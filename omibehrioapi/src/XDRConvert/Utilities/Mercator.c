/******************************************************************************
PURPOSE: Mercator.c - Define Mercator Conformal Conic Projectors ADT.
NOTES:   See Projection.h and Mercator.h. Formulations from the USGS PROJ Lib.
HISTORY: 2004/10, Todd Plessel Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <math.h>  /* For M_PI, sqrt(), sin(), cos(), atan2(), log(), pow().*/
#include <float.h> /* For DBL_MAX. */

#include <Assertions.h>    /* For PRE*(), POST*(), CHECK*(), IN_RANGE().     */
#include <BasicNumerics.h> /* For Integer, Real, isNan(), withinTolerance(). */
#include <Memory.h>        /* For NEW_ZERO(), ZERO_OBJECT(), FREE_ZERO().    */
#include <Mercator.h>      /* For public interface.                          */

/*================================== TYPES ==================================*/

struct MercatorPrivate {
  Real majorSemiaxis;    /* Mean equitorial radius in meters 6370000.0. */
  Real minorSemiaxis;    /* Mean polar      radius in meters 6370000.0. */
  Real centralLongitude; /* Projects to zero, e.g., -100.0 degrees.     */
  Real falseEasting;     /* Skew offset in meters, e.g., 0.0.           */
  Real falseNorthing;    /* Skew offset in meters, e.g., 0.0.           */
  Real eccentricity;     /* Of ellipsoid approximation of planet.       */
  Real lambda0;          /* Central longitude in radians.               */
};

/*=========================== FORWARD DECLARATIONS ==========================*/

/* Commands: */

static void free__( Mercator* self );

static void setEllipsoid( Mercator* self,
                          Real majorSemiaxis, Real minorSemiaxis );

static void setFalseEasting( Mercator* self, Real falseEasting );

static void setFalseNorthing( Mercator* self, Real falseNorthing );

static void project( const Mercator* self, Real longitude, Real latitude,
                     Real* x, Real* y );

static void unproject( const Mercator* self, Real x, Real y,
                       Real* longitude, Real* latitude );

/* Queries: */

static Integer invariant( const Mercator* self );
static Integer equal( const Mercator* self, const Mercator* other );
static Mercator* clone( const Mercator* self );

static void ellipsoid( const Mercator* self,
                       Real* majorSemiaxis, Real* minorSemiaxis );

static Real falseEasting( const Mercator* self );
static Real falseNorthing( const Mercator* self );
static Real centralLongitude( const Mercator* self );
static Real centralLatitude( const Mercator* self );
static const char* name( const Mercator* self );

/* Helpers: */

static void assignMembers( Mercator* self );
static Integer hasMembers( const Mercator* self );
static void computeDerivedTerms( Mercator* self );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: newMercator - Construct a Mercator projector.
INPUTS:  Real majorSemiaxis    Mean equitorial radius in meters 6370000.0.
         Real minorSemiaxis    Mean polar      radius in meters 6370000.0.
         Real centralLongitude Projects to zero, e.g., -100.0 degrees.
         Real falseEasting     Skew offset in meters, e.g., 0.0.
         Real falseNorthing    Skew offset in meters, e.g., 0.0.
RETURNS: Mercator* else 0 if failed and failureMessage() called.
******************************************************************************/

Mercator* newMercator( Real newMajorSemiaxis,    Real newMinorSemiaxis,
                       Real newCentralLongitude,
                       Real newFalseEasting,     Real newFalseNorthing ) {

  PRE04( isValidEllipsoid( newMajorSemiaxis, newMinorSemiaxis ),
         isValidLongitude( newCentralLongitude ),
         ! isNan( newFalseEasting ),
         ! isNan( newFalseNorthing ) );

  Mercator* result = 0;
  MercatorPrivate* data = NEW_ZERO( MercatorPrivate, 1 );

  if ( data ) {
    result = NEW_ZERO( Mercator, 1 );

    if ( ! result ) {
      FREE( data );
    } else {
      data->majorSemiaxis    = newMajorSemiaxis;
      data->minorSemiaxis    = newMinorSemiaxis;
      data->centralLongitude = newCentralLongitude;
      data->falseEasting     = newFalseEasting;
      data->falseNorthing    = newFalseNorthing;
      result->data = data;
      assignMembers( result );
      computeDerivedTerms( result );
    }
  }

  POST0( IMPLIES( result,
                  AND3( result->invariant( result ),
                        result->falseEasting( result ) == newFalseEasting,
                        result->falseNorthing( result) == newFalseNorthing)));
  return result;
}



/******************************************************************************
PURPOSE: free - Destruct a Mercator.
INPUTS:  Mercator* self  Object to destruct.
NOTE:    Use FREE_OBJECT( object ) instead since it zeros argument.
******************************************************************************/

static void free__( Mercator* self ) {
  PRE( self );
  FREE_ZERO( self->data );
  FREE_ZERO( self );
  POST0( self == 0 );
}



/******************************************************************************
PURPOSE: setEllipsoid - Set the ellipsoid approximation of planet.
INPUTS:  Real     newMajorSemiaxis  Mean equitorial radius in meters.
         Real     newMinorSemiaxis  Mean polar      radius in meters.
OUTPUTS: Mercator* self  Object to update.
******************************************************************************/

static void setEllipsoid( Mercator* self,
                          Real newMajorSemiaxis, Real newMinorSemiaxis ) {

  PRE2( self, isValidEllipsoid( newMajorSemiaxis, newMinorSemiaxis ) );
  self->data->majorSemiaxis = newMajorSemiaxis;
  self->data->minorSemiaxis = newMinorSemiaxis;
  computeDerivedTerms( self );
}



/******************************************************************************
PURPOSE: setFalseEasting - Set the projected x offset in meters.
INPUTS:  Real     newFalseEasting  Projected x offset in meters.
OUTPUTS: Mercator* self  Object to update.
******************************************************************************/

static void setFalseEasting( Mercator* self, Real newFalseEasting ) {

  PRE2( self, ! isNan( newFalseEasting ) );
  self->data->falseEasting = newFalseEasting;
  POST( aboutEqual( self->falseEasting( self ), newFalseEasting ) );
}



/******************************************************************************
PURPOSE: setFalseNorthing - Set the projected y offset in meters.
INPUTS:  Real     newFalseNorthing  Projected y offset in meters.
OUTPUTS: Mercator* self  Object to update.
******************************************************************************/

static void setFalseNorthing( Mercator* self, Real newFalseNorthing ) {

  PRE2( self, ! isNan( newFalseNorthing ) );
  self->data->falseNorthing = newFalseNorthing;
  POST( aboutEqual( self->falseNorthing( self ), newFalseNorthing ) );
}



/******************************************************************************
PURPOSE: project - Project a point.
INPUTS:  const Mercator* self  Projector.
         Real longitude  E.g., -78.7268.
         Real latitude   E.g., 35.9611.
OUTPUTS: Real* x         Projected longitude.
         Real* y         Projected latitude.
******************************************************************************/

static void project( const Mercator* self, Real longitude, Real latitude,
                     Real* x, Real* y ) {

  PRE4( isValidLongitude( longitude ), isValidLatitude( latitude ), x, y );

  const MercatorPrivate* const data = self->data;
  Real lambda = radians( longitude );
  Real phi    = radians( latitude  );
  Real lambdaDelta = 0.0;  /* Radians from central latitude. */

  /*
   * If phi is too near a pole tweak it so that projecting
   * succeeds and unprojecting yields original latitude.
   */

  if ( ! IN_RANGE( phi, -PI_OVER_2 + TOLERANCE, PI_OVER_2 - TOLERANCE ) ) {
    phi = phi + TOLERANCE * -SIGN( phi );
  }


  /*
   * If lambda is too near +/-180 longitude tweak it so that projecting
   * succeeds and unprojecting yields original longitude.
   */

  if ( ! IN_RANGE( lambda, -M_PI + TOLERANCE, M_PI - TOLERANCE ) ) {
    lambda = lambda + TOLERANCE * -SIGN( lambda );
  }

  for ( lambdaDelta = lambda - data->lambda0; fabs( lambdaDelta ) > M_PI; ) {

    if ( lambdaDelta < 0.0 ) {
      lambdaDelta = lambdaDelta + M_PI + M_PI;
    } else {
      lambdaDelta = lambdaDelta - M_PI - M_PI;
    }
  }

  *x = lambdaDelta * data->majorSemiaxis + data->falseEasting;

  if ( data->eccentricity == 0.0 ) { /* Sphere: */
    *y = log( tan( PI_OVER_4 + phi * 0.5 ) );
  } else { /* Ellipsoid: */
    *y = -log( tsfn( phi, sin( phi ), data->eccentricity ) );
  }

  *y = *y * data->majorSemiaxis + data->falseNorthing;

  POST2( ! isNan( *x ), ! isNan( *y ) );
}



/******************************************************************************
PURPOSE: unproject - Unproject a point.
INPUTS:  const Mercator* self   Projector.
         Real x           X-coordinate of point to unproject.
         Real y           X-coordinate of point to unproject.
OUTPUTS: Real* longitude  Unprojected x.
         Real* latitude   Unprojected y.
******************************************************************************/

static void unproject( const Mercator* self, Real x, Real y,
                       Real* longitude, Real* latitude ) {

  PRE4( ! isNan( x ), ! isNan( y ), longitude, latitude );

  const MercatorPrivate* const data = self->data;
  Real oneOverMajorSemiaxis = 1.0 / data->majorSemiaxis;
  Real xp = x - data->falseEasting;
  Real yp = y - data->falseNorthing;
  Real expYP = exp( -yp * oneOverMajorSemiaxis );
  Real lambda = xp * oneOverMajorSemiaxis + data->lambda0;
  Real phi =
    data->eccentricity == 0.0 ? PI_OVER_2 - 2.0 * atan( expYP )
    :                     phi2Iterate( expYP, data->eccentricity );

  *longitude = degrees( lambda );
  *latitude = degrees( phi );
  CHECK( fabs( *longitude ) < REAL_MAX );

  while ( *longitude < -180.0 ) {
    *longitude += 360.0;
  }

  while ( *longitude > 180.0 ) {
    *longitude -= 360.0;
  }

  POST2( isValidLongitude( *longitude ), isValidLatitude( *latitude ) );
}



/******************************************************************************
PURPOSE: invariant - Class invariant.
INPUTS:  const Mercator* self  Projector to check.
RETURNS: Integer 1 if valid, else 0.
NOTE:    If this query ever returns 0 then there is a defect in the code.
******************************************************************************/

static Integer invariant( const Mercator* self ) {

  const MercatorPrivate* const data = self ? self->data : 0;
  const Integer result =
    AND8( self,
          hasMembers( self ),
          self->data,
          data == self->data,
          isValidEllipsoid( data->majorSemiaxis, data->minorSemiaxis ),
          isValidLongitude( data->centralLongitude ),
          ! isNan( data->falseEasting ),
          ! isNan( data->falseNorthing ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: equal - Is self functionally equivalent to other?
INPUTS:  const Mercator* self   Projector to compare.
         const Mercator* other  Projector to compare.
RETURNS: Integer 1 if equal, else 0.
******************************************************************************/

static Integer equal( const Mercator* self, const Mercator* other ) {

  PRE2( other, other->invariant( other ) );

  const MercatorPrivate* const data = self->data;
  const MercatorPrivate* const otherData = other->data;
  const Integer result =
    AND5( aboutEqual( data->majorSemiaxis,    otherData->majorSemiaxis ),
          aboutEqual( data->minorSemiaxis,    otherData->minorSemiaxis ),
          aboutEqual( data->centralLongitude, otherData->centralLongitude ),
          aboutEqual( data->falseEasting,     otherData->falseEasting ),
          aboutEqual( data->falseNorthing,    otherData->falseNorthing ) );

  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: clone - Yield a clone - a NEW functionally equivalent Mercator.
INPUTS:  const Mercator* self   Projector to clone.
RETURNS: Mercator* result clone of self.
******************************************************************************/

static Mercator* clone( const Mercator* self ) {

  PRE( self );

  const MercatorPrivate* const data = self->data;
  Mercator* const result =
    newMercator( data->majorSemiaxis, data->minorSemiaxis,
                 data->centralLongitude,
                 data->falseEasting, data->falseNorthing );

  POST( IMPLIES( result,
                 AND2( result->invariant( result ),
                       result->equal( result, self ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: ellipsoid - The ellipsoid axes of the planet approximation.
INPUTS:  const Mercator* self   Projector to query.
OUTPUTS: Real* majorSemiaxis   Mean equitorial radius.
OUTPUTS: Real* minorSemiaxis   Mean polar      radius.
******************************************************************************/

static void ellipsoid( const Mercator* self,
                       Real* majorSemiaxis, Real* minorSemiaxis ) {

  PRE3( self, majorSemiaxis, minorSemiaxis );

  const MercatorPrivate* const data = self->data;
  *majorSemiaxis = data->majorSemiaxis;
  *minorSemiaxis = data->majorSemiaxis;

  POST( isValidEllipsoid( *majorSemiaxis, *minorSemiaxis ) );
}



/******************************************************************************
PURPOSE: falseEasting - Projected x offset in meters.
INPUTS:  const Mercator* self   Projector to query.
RETURNS: Real result  Projected x offset in meters.
******************************************************************************/

static Real falseEasting( const Mercator* self ) {

  PRE( self );

  const Real result = self->data->falseEasting;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: falseNorthing - Projected y offset in meters.
INPUTS:  const Mercator* self   Projector to query.
RETURNS: Real result  Projected y offset in meters.
******************************************************************************/

static Real falseNorthing( const Mercator* self ) {

  PRE( self );

  const Real result = self->data->falseNorthing;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLongitude - Longitude of center of projection.
INPUTS:  const Mercator* self   Projector to query.
RETURNS: Real result longitude of center of projection.
******************************************************************************/

static Real centralLongitude( const Mercator* self ) {

  PRE( self );

  const Real result = self->data->centralLongitude;

  POST( isValidLongitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLatitude - Latitude of center of projection.
INPUTS:  const Mercator* self   Projector to query.
RETURNS: Real result latitude of center of projection.
******************************************************************************/

static Real centralLatitude( const Mercator* self ) {

  PRE( self );

  const Real result = 0.0;

  POST( result == 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: name - Name of projection.
INPUTS:  const Mercator* self   Projector to query.
RETURNS: const char*  Name of projection: "Mercator".
******************************************************************************/

static const char* name( const Mercator* self ) {

  PRE( self );

  const char* const result = "Mercator";

  POST( ! strcmp( result, "Mercator" ) );
  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: assignMembers - Assign pointers to member functions.
OUTPUTS: Mercator* self  Mercator to initialize.
******************************************************************************/

static void assignMembers( Mercator* self ) {

  PRE0( self );

  self->free             = free__;
  self->setEllipsoid     = setEllipsoid;
  self->setFalseEasting  = setFalseEasting;
  self->setFalseNorthing = setFalseNorthing;
  self->project          = project;
  self->unproject        = unproject;
  self->invariant        = invariant;
  self->equal            = equal;
  self->clone            = clone;
  self->ellipsoid        = ellipsoid;
  self->falseEasting     = falseEasting;
  self->falseNorthing    = falseNorthing;
  self->centralLongitude = centralLongitude;
  self->centralLatitude  = centralLatitude;
  self->name             = name;

  POST0( hasMembers( self ) );
}



/******************************************************************************
PURPOSE: hasMembers - Does self have all pointers to member functions?
INPUTS:  const Mercator* self  Mercator to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer hasMembers( const Mercator* self ) {

  PRE0( self );

  const Integer result =
    AND15( self->free             == free__,
           self->setEllipsoid     == setEllipsoid,
           self->setFalseEasting  == setFalseEasting,
           self->setFalseNorthing == setFalseNorthing,
           self->project          == project,
           self->unproject        == unproject,
           self->invariant        == invariant,
           self->equal            == equal,
           self->clone            == clone,
           self->ellipsoid        == ellipsoid,
           self->falseEasting     == falseEasting,
           self->falseNorthing    == falseNorthing,
           self->centralLongitude == centralLongitude,
           self->centralLatitude  == centralLatitude,
           self->name             == name );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: computeDerivedTerms - Compute trigonometry terms independent of
         longitude/latitude of projection point.
OUTPUTS: Mercator* self  Mercator to (re)initialize.
******************************************************************************/

static void computeDerivedTerms( Mercator* self ) {

  PRE0( self );

  MercatorPrivate* const data = self->data;
  const Real eccentricity0 = data->majorSemiaxis == data->minorSemiaxis ? 0.0 :
    safeQuotient( sqrt( safeDifference( SQUARE( data->majorSemiaxis ),
                                        SQUARE( data->minorSemiaxis ) ) ),
                  data->majorSemiaxis );

  const Real eccentricity1 = eccentricity0 > 1.0 ? 1.0 : eccentricity0;
  data->eccentricity = eccentricity1;
  data->lambda0 = radians( data->centralLongitude );

  POST3( ! isNan( data->eccentricity ),
         IN_RANGE( data->eccentricity, 0.0, 1.0 ),
         ! isNan( data->lambda0 ) );
}





/******************************************************************************
PURPOSE: Lambert.c - Define Lambert Conformal Conic Projectors ADT.
NOTES:   See Projection.h and Lambert.h. Formulations from the USGS PROJ Lib.
HISTORY: 2004/10, Todd Plessel Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <math.h>  /* For M_PI, sqrt(), sin(), cos(), atan2(), log(), pow().*/
#include <float.h> /* For DBL_MAX. */

#include <Assertions.h>    /* For PRE*(), POST*(), CHECK*(), IN_RANGE().     */
#include <BasicNumerics.h> /* For Integer, Real, isNan(), withinTolerance(). */
#include <Memory.h>        /* For NEW_ZERO(), ZERO_OBJECT(), FREE_ZERO().    */
#include <Lambert.h>       /* For public interface.                          */

/*================================== TYPES ==================================*/

struct LambertPrivate {
  Real majorSemiaxis;    /* Mean equitorial radius in meters 6370000.0. */
  Real minorSemiaxis;    /* Mean polar      radius in meters 6370000.0. */
  Real lowerLatitude;    /* Lower tangent in degrees, e.g., 30.0.       */
  Real upperLatitude;    /* Upper tangent in degrees, e.g., 60.0.       */
  Real centralLongitude; /* Projects to zero, e.g., -100.0 degrees.     */
  Real centralLatitude;  /* Projects to zero, e.g., 40.0 degrees.       */
  Real falseEasting;     /* Skew offset in meters, e.g., 0.0.           */
  Real falseNorthing;    /* Skew offset in meters, e.g., 0.0.           */
  Real eccentricity;     /* Of ellipsoid approximation of planet.       */
  Real lambda0;          /* Central longitude in radians.               */
  Real rho0;             /* See USGS PROJ Library.                      */
  Real n;                /* See USGS PROJ Library.                      */
  Real c;                /* See USGS PROJ Library.                      */
};

/*=========================== FORWARD DECLARATIONS ==========================*/

/* Commands: */

static void free__( Lambert* self );

static void setEllipsoid( Lambert* self,
                          Real majorSemiaxis, Real minorSemiaxis );

static void setFalseEasting( Lambert* self, Real falseEasting );

static void setFalseNorthing( Lambert* self, Real falseNorthing );

static void project( const Lambert* self, Real longitude, Real latitude,
                     Real* x, Real* y );

static void unproject( const Lambert* self, Real x, Real y,
                       Real* longitude, Real* latitude );

/* Queries: */

static Integer invariant( const Lambert* self );
static Integer equal( const Lambert* self, const Lambert* other );
static Lambert* clone( const Lambert* self );

static void ellipsoid( const Lambert* self,
                       Real* majorSemiaxis, Real* minorSemiaxis );

static Real falseEasting( const Lambert* self );
static Real falseNorthing( const Lambert* self );
static Real lowerLatitude( const Lambert* self );
static Real upperLatitude( const Lambert* self );
static Real centralLongitude( const Lambert* self );
static Real centralLatitude( const Lambert* self );
static const char* name( const Lambert* self );

/* Helpers: */

static void assignMembers( Lambert* self );
static Integer hasMembers( const Lambert* self );
static void computeDerivedTerms( Lambert* self );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: newLambert - Construct a Lambert projector.
INPUTS:  Real majorSemiaxis    Mean equitorial radius in meters 6370000.0.
         Real minorSemiaxis    Mean polar      radius in meters 6370000.0.
         Real lowerLatitude    Lower tangent in degrees, e.g., 30.0.
         Real upperLatitude    Upper tangent in degrees, e.g., 60.0.
         Real centralLongitude Projects to zero, e.g., -100.0 degrees.
         Real centralLatitude  Projects to zero, e.g., 40.0 degrees.
         Real falseEasting     Skew offset in meters, e.g., 0.0.
         Real falseNorthing    Skew offset in meters, e.g., 0.0.
RETURNS: Lambert* else 0 if failed and failureMessage() called.
******************************************************************************/

Lambert* newLambert( Real newMajorSemiaxis,    Real newMinorSemiaxis,
                     Real newLowerLatitude,    Real newUpperLatitude,
                     Real newCentralLongitude, Real newCentralLatitude,
                     Real newFalseEasting,     Real newFalseNorthing ) {

  PRE012( isValidEllipsoid( newMajorSemiaxis, newMinorSemiaxis ),
          isValidLatitude( newLowerLatitude ),
          isValidLatitude( newUpperLatitude ),
          isValidLongitude( newCentralLongitude ),
          isValidLatitude( newCentralLatitude ),
          newLowerLatitude <= newUpperLatitude,
          SIGN( newLowerLatitude ) == SIGN( newUpperLatitude ),
          IMPLIES_ELSE( newLowerLatitude >= 0.0,
                        IN_RANGE( newLowerLatitude, 1.0, 89.0 ),
                        IN_RANGE( newLowerLatitude, -89.0, -1.0 ) ),
          IMPLIES_ELSE( newUpperLatitude >= 0.0,
                        IN_RANGE( newUpperLatitude, 1.0, 89.0 ),
                        IN_RANGE( newUpperLatitude, -89.0, -1.0 ) ),
          IN_RANGE( newCentralLatitude, -89.0, 89.0 ),
          ! isNan( newFalseEasting ),
          ! isNan( newFalseNorthing ) );

  Lambert* result = 0;
  LambertPrivate* data = NEW_ZERO( LambertPrivate, 1 );

  if ( data ) {
    result = NEW_ZERO( Lambert, 1 );

    if ( ! result ) {
      FREE( data );
    } else {
      data->majorSemiaxis    = newMajorSemiaxis;
      data->minorSemiaxis    = newMinorSemiaxis;
      data->lowerLatitude    = newLowerLatitude;
      data->upperLatitude    = newUpperLatitude;
      data->centralLongitude = newCentralLongitude;
      data->centralLatitude  = newCentralLatitude;
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
PURPOSE: free - Destruct a Lambert.
INPUTS:  Lambert* self  Object to destruct.
NOTE:    Use FREE_OBJECT( object ) instead since it zeros argument.
******************************************************************************/

static void free__( Lambert* self ) {
  PRE( self );
  FREE_ZERO( self->data );
  FREE_ZERO( self );
  POST0( self == 0 );
}



/******************************************************************************
PURPOSE: setEllipsoid - Set the ellipsoid approximation of planet.
INPUTS:  Real     newMajorSemiaxis  Mean equitorial radius in meters.
         Real     newMinorSemiaxis  Mean polar      radius in meters.
OUTPUTS: Lambert* self  Object to update.
******************************************************************************/

static void setEllipsoid( Lambert* self,
                          Real newMajorSemiaxis, Real newMinorSemiaxis ) {

  PRE2( self, isValidEllipsoid( newMajorSemiaxis, newMinorSemiaxis ) );
  self->data->majorSemiaxis = newMajorSemiaxis;
  self->data->minorSemiaxis = newMinorSemiaxis;
  computeDerivedTerms( self );
}



/******************************************************************************
PURPOSE: setFalseEasting - Set the projected x offset in meters.
INPUTS:  Real     newFalseEasting  Projected x offset in meters.
OUTPUTS: Lambert* self  Object to update.
******************************************************************************/

static void setFalseEasting( Lambert* self, Real newFalseEasting ) {

  PRE2( self, ! isNan( newFalseEasting ) );
  self->data->falseEasting = newFalseEasting;
  POST( aboutEqual( self->falseEasting( self ), newFalseEasting ) );
}



/******************************************************************************
PURPOSE: setFalseNorthing - Set the projected y offset in meters.
INPUTS:  Real     newFalseNorthing  Projected y offset in meters.
OUTPUTS: Lambert* self  Object to update.
******************************************************************************/

static void setFalseNorthing( Lambert* self, Real newFalseNorthing ) {

  PRE2( self, ! isNan( newFalseNorthing ) );
  self->data->falseNorthing = newFalseNorthing;
  POST( aboutEqual( self->falseNorthing( self ), newFalseNorthing ) );
}



/******************************************************************************
PURPOSE: project - Project a point.
INPUTS:  const Lambert* self   Projector.
         Real longitude  E.g., -78.7268.
         Real latitude   E.g., 35.9611.
OUTPUTS: Real* x         Projected longitude.
         Real* y         Projected latitude.
******************************************************************************/

static void project( const Lambert* self, Real longitude, Real latitude,
                     Real* x, Real* y ) {

  PRE4( isValidLongitude( longitude ), isValidLatitude( latitude ), x, y );

  const LambertPrivate* const data = self->data;
  Real lambda = radians( longitude );
  Real phi    = radians( latitude  );
  Real rho = 0.0;    /* See USGS PROJ Library. */
  Real lambdaDelta = 0.0;  /* Radians from central latitude. */
  Real nLambdaDelta = 0.0; /* n times lambdaDelta. */

  /*
   * If phi is too near a pole tweak it so that projecting
   * succeeds and unprojecting yields original longitude
   * (instead of central longitude).
   */

  if ( ! IN_RANGE( phi, -PI_OVER_2 + PROJECTION_TOLERANCE,
                         PI_OVER_2 - PROJECTION_TOLERANCE ) ) {
    phi = phi + sqrt( PROJECTION_TOLERANCE ) * -SIGN( phi );
  }

  rho = data->c * pow( tsfn( phi, sin( phi ), data->eccentricity ), data->n );

  /*
   * If lambda is too near +/-180 longitude tweak it so that projecting
   * succeeds and unprojecting yields original longitude
   * (instead of central longitude).
   */

  if ( ! IN_RANGE( lambda, -M_PI + PROJECTION_TOLERANCE,
                            M_PI - PROJECTION_TOLERANCE ) ) {
    lambda = lambda + sqrt( PROJECTION_TOLERANCE ) * -SIGN( lambda );
  }

  for ( lambdaDelta = lambda - data->lambda0; fabs( lambdaDelta ) > M_PI; ) {

    if ( lambdaDelta < 0.0 ) {
      lambdaDelta = lambdaDelta + M_PI + M_PI;
    } else {
      lambdaDelta = lambdaDelta - M_PI - M_PI;
    }
  }

  nLambdaDelta = data->n * lambdaDelta;
  *x = rho          * sin( nLambdaDelta ) * data->majorSemiaxis +
       data->falseEasting;
  *y = ( data->rho0 - rho * cos( nLambdaDelta)) * data->majorSemiaxis +
       data->falseNorthing;

  POST2( ! isNan( *x ), ! isNan( *y ) );
}



/******************************************************************************
PURPOSE: unproject - Unproject a point.
INPUTS:  const Lambert* self    Projector.
         Real x           X-coordinate of point to unproject.
         Real y           X-coordinate of point to unproject.
OUTPUTS: Real* longitude  Unprojected x.
         Real* latitude   Unprojected y.
******************************************************************************/

static void unproject( const Lambert* self, Real x, Real y,
                       Real* longitude, Real* latitude ) {

  PRE4( ! isNan( x ), ! isNan( y ), longitude, latitude );

  const LambertPrivate* const data = self->data;
  Real oneOverMajorSemiaxis = 1.0 / data->majorSemiaxis;
  Real xp = ( x - data->falseEasting  ) * oneOverMajorSemiaxis;
  Real yp = ( y - data->falseNorthing ) * oneOverMajorSemiaxis;
  Real ypDelta = data->rho0 - yp;/*Dist from yp to central lat (in radians)*/
  Real rho     = hypot( xp, ypDelta );
  Real lambda  = 0.0;       /* Radians of longitude. */
  Real phi     = PI_OVER_2; /* Radians of latitude.  */

  if ( rho != 0.0 ) {

    if ( data->n < 0.0 ) {
      rho = -rho;
      xp = -xp;
      ypDelta = -ypDelta;
    }

    CHECK3( data->c != 0.0, data->n != 0.0, rho != 0.0 );

    if ( data->eccentricity == 0.0 ) { /* Sphere: */
      phi = 2.0 * atan( pow( data->c / rho, 1.0 / data->n ) ) - PI_OVER_2;
    } else { /* Ellipsoid: */
      phi = phi2Iterate( pow( rho / data->c, 1.0 / data->n ),
                         data->eccentricity );
    }

    lambda = atan2( xp, ypDelta ) / data->n;
  } else if ( data->n < 0.0 ) {
    phi = -PI_OVER_2;
  }

  lambda += data->lambda0;
  *longitude = degrees( lambda );
  *latitude = degrees( phi );
  CHECK( fabs( *longitude ) < DBL_MAX );

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
INPUTS:  const Lambert* self  Projector to check.
RETURNS: Integer 1 if valid, else 0.
NOTE:    If this query ever returns 0 then there is a defect in the code.
******************************************************************************/

static Integer invariant( const Lambert* self ) {

  const LambertPrivate* const data = self ? self->data : 0;
  const Integer result =
    AND15( self, hasMembers( self ), self->data, data == self->data,
           isValidEllipsoid( data->majorSemiaxis, data->minorSemiaxis ),
           isValidLatitude( data->lowerLatitude ),
           isValidLatitude( data->upperLatitude ),
           isValidLatitude( data->centralLatitude ),
           isValidLongitude( data->centralLongitude ),
           data->lowerLatitude <= data->upperLatitude,
           SIGN( data->lowerLatitude ) == SIGN( data->upperLatitude ),
           IMPLIES_ELSE( data->lowerLatitude >= 0.0,
                         IN_RANGE( data->lowerLatitude, 1.0, 89.0 ),
                         IN_RANGE( data->lowerLatitude, -89.0, -1.0 ) ),
           IMPLIES_ELSE( data->upperLatitude >= 0.0,
                         IN_RANGE( data->upperLatitude, 1.0, 89.0 ),
                         IN_RANGE( data->upperLatitude, -89.0, -1.0 ) ),
           ! isNan( data->falseEasting ),
           ! isNan( data->falseNorthing ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: equal - Is self functionally equivalent to other?
INPUTS:  const Lambert* self   Projector to compare.
         const Lambert* other  Projector to compare.
RETURNS: Integer 1 if equal, else 0.
******************************************************************************/

static Integer equal( const Lambert* self, const Lambert* other ) {

  PRE2( other, other->invariant( other ) );

  const LambertPrivate* const data = self->data;
  const LambertPrivate* const otherData = other->data;
  const Integer result =
    AND8( aboutEqual( data->majorSemiaxis,    otherData->majorSemiaxis ),
          aboutEqual( data->minorSemiaxis,    otherData->minorSemiaxis ),
          aboutEqual( data->lowerLatitude,    otherData->lowerLatitude ),
          aboutEqual( data->upperLatitude,    otherData->upperLatitude ),
          aboutEqual( data->centralLatitude,  otherData->centralLatitude ),
          aboutEqual( data->centralLongitude, otherData->centralLongitude ),
          aboutEqual( data->falseEasting,     otherData->falseEasting ),
          aboutEqual( data->falseNorthing,    otherData->falseNorthing ) );

  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: clone - Yield a clone - a NEW functionally equivalent Lambert.
INPUTS:  const Lambert* self   Projector to clone.
RETURNS: Lambert* result clone of self.
******************************************************************************/

static Lambert* clone( const Lambert* self ) {

  PRE( self );

  const LambertPrivate* const data = self->data;
  Lambert* const result =
    newLambert( data->majorSemiaxis,    data->minorSemiaxis,
                data->lowerLatitude,    data->upperLatitude,
                data->centralLongitude, data->centralLatitude,
                data->falseEasting,     data->falseNorthing );

  POST( IMPLIES( result,
                 AND2( result->invariant( result ),
                       result->equal( result, self ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: ellipsoid - The ellipsoid axes of the planet approximation.
INPUTS:  const Lambert* self   Projector to query.
OUTPUTS: Real* majorSemiaxis   Mean equitorial radius.
OUTPUTS: Real* minorSemiaxis   Mean polar      radius.
******************************************************************************/

static void ellipsoid( const Lambert* self,
                       Real* majorSemiaxis, Real* minorSemiaxis ) {

  PRE3( self, majorSemiaxis, minorSemiaxis );

  const LambertPrivate* const data = self->data;
  *majorSemiaxis = data->majorSemiaxis;
  *minorSemiaxis = data->majorSemiaxis;

  POST( isValidEllipsoid( *majorSemiaxis, *minorSemiaxis ) );
}



/******************************************************************************
PURPOSE: falseEasting - Projected x offset in meters.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: Real result  Projected x offset in meters.
******************************************************************************/

static Real falseEasting( const Lambert* self ) {

  PRE( self );

  const Real result = self->data->falseEasting;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: falseNorthing - Projected y offset in meters.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: Real result  Projected y offset in meters.
******************************************************************************/

static Real falseNorthing( const Lambert* self ) {

  PRE( self );

  const Real result = self->data->falseNorthing;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLongitude - Longitude of center of projection.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: Real result longitude of center of projection.
******************************************************************************/

static Real centralLongitude( const Lambert* self ) {

  PRE( self );

  const Real result = self->data->centralLongitude;

  POST( isValidLongitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLatitude - Latitude of center of projection.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: Real result latitude of center of projection.
******************************************************************************/

static Real centralLatitude( const Lambert* self ) {

  PRE( self );

  const Real result = self->data->centralLatitude;

  POST( isValidLatitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: lowerLatitude - Lower latitude of secant plane.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: Real result lower latitude of secant plane.
******************************************************************************/

static Real lowerLatitude( const Lambert* self ) {

  PRE( self );

  const Real result = self->data->lowerLatitude;

  POST( isValidLatitude( result ) );
  return result;
}




/******************************************************************************
PURPOSE: name - Name of projection.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: const char*  Name of projection: "Lambert".
******************************************************************************/

static const char* name( const Lambert* self ) {

  PRE( self );

  const char* const result = "Lambert";

  POST( ! strcmp( result, "Lambert" ) );
  return result;
}




/******************************************************************************
PURPOSE: upperLatitude - Upper latitude of secant plane.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: Real result upper latitude of secant plane.
******************************************************************************/

static Real upperLatitude( const Lambert* self ) {

  PRE( self );

  const Real result = self->data->upperLatitude;

  POST2( isValidLatitude( result ), result >= self->lowerLatitude( self ) );
  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: assignMembers - Assign pointers to member functions.
OUTPUTS: Lambert* self  Lambert to initialize.
******************************************************************************/

static void assignMembers( Lambert* self ) {

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

  self->lowerLatitude    = lowerLatitude;
  self->upperLatitude    = upperLatitude;

  POST0( hasMembers( self ) );
}



/******************************************************************************
PURPOSE: hasMembers - Does self have all pointers to member functions?
INPUTS:  const Lambert* self  Lambert to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer hasMembers( const Lambert* self ) {

  PRE0( self );

  const Integer result =
    AND17( self->free             == free__,
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
           self->name             == name,
           self->lowerLatitude    == lowerLatitude,
           self->upperLatitude    == upperLatitude );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: computeDerivedTerms - Compute trigonometry terms independent of
         longitude/latitude of projection point.
OUTPUTS: Lambert* self  Lambert to (re)initialize.
******************************************************************************/

static void computeDerivedTerms( Lambert* self ) {

  PRE0( self );

  LambertPrivate* const data = self->data;
  const Real eccentricity0 = data->majorSemiaxis == data->minorSemiaxis ? 0.0 :
    safeQuotient( sqrt( safeDifference( SQUARE( data->majorSemiaxis ),
                                        SQUARE( data->minorSemiaxis ) ) ),
                  data->majorSemiaxis );

  const Real eccentricity1 = eccentricity0 > 1.0 ? 1.0 : eccentricity0;
  const Real eccentricitySquared = SQUARE( eccentricity1 );
  const Real phi0 = radians( data->centralLatitude );
  const Real phi1 = radians( data->lowerLatitude );
  const Real phi2 = radians( data->upperLatitude );
  const Real sinePhi1   = sin( phi1 );
  const Real cosinePhi1 = cos( phi1 );
  const Real sinePhi2   = sin( phi2 );
  const Real cosinePhi2 = cos( phi2 );
  /* Are lower/upperLatitude about equal? */
  const Integer isTangent = phi1 + PROJECTION_TOLERANCE >= phi2;

  data->eccentricity = eccentricity1;
  data->lambda0 = radians( data->centralLongitude );
  data->n = sinePhi1;

  if ( eccentricitySquared != 0.0 ) { /* Ellipsoid planet: */
    const Real m1  = msfn( sinePhi1, cosinePhi1, eccentricitySquared );
    const Real ml1 = tsfn( phi1, sinePhi1, data->eccentricity );

    if ( ! isTangent ) { /* Secant form: */
      const Real numerator =
        log( m1 / msfn( sinePhi2, cosinePhi2, eccentricitySquared ) );
      const Real denominator =
        log( ml1 / tsfn( phi2, sinePhi2, data->eccentricity ) );
      CHECK( denominator != 0.0 );
      data->n = numerator / denominator;
    }

    CHECK( data->n != 0.0 );
    data->c = m1 * pow( ml1, -data->n ) / data->n;

    if ( fabs( fabs( phi0 ) - PI_OVER_2 ) < PROJECTION_TOLERANCE ) {
      data->rho0 = 0.0;
    } else {
      data->rho0 = data->c * pow( tsfn( phi0, sin( phi0 ), data->eccentricity),
                             data->n );
    }
  } else { /* Sphere planet: */
    const Real denominator = tan( PI_OVER_4 + 0.5 * phi1 );

    if ( ! isTangent ) { /* Secant form: */
      CHECK( ! aboutEqual( fabs( phi1 ), PI_OVER_2 ) ); /* Not near pole. */
      CHECK( ! aboutEqual( fabs( phi2 ), PI_OVER_2 ) ); /* Not near pole. */
      CHECK( cosinePhi1 != 0.0 );
      CHECK( cosinePhi2 != 0.0 );
      CHECK( tan( PI_OVER_4 + 0.5 * phi2 ) != 0.0 );
      CHECK( denominator != 0.0 );
      data->n = log( cosinePhi1 / cosinePhi2 ) /
                log( tan( PI_OVER_4 + 0.5 * phi2 ) / denominator );
    }

    data->c = cosinePhi1 * pow( denominator, data->n ) / data->n;

    if ( fabs( fabs( phi0 ) - PI_OVER_2 ) < PROJECTION_TOLERANCE ) {
      data->rho0 = 0.0;
    } else {
      data->rho0 = data->c * pow( tan( PI_OVER_4 + 0.5 * phi0 ), -data->n );
    }
  }

  POST6( ! isNan( data->eccentricity ),
         IN_RANGE( data->eccentricity, 0.0, 1.0 ),
         ! isNan( data->lambda0 ),
         ! isNan( data->rho0 ),
         ! isNan( data->n ),
         ! isNan( data->c ) );
}




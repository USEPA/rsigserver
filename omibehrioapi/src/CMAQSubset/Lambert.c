/******************************************************************************
PURPOSE: Lambert.c - Define Lambert Conformal Conic Projectors ADT.
NOTES:   See Projection.h and Lambert.h. Formulations from the USGS PROJ Lib.
HISTORY: 2004-10-01 plessel.todd@epa.gov Created based on C++ version.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <math.h>  /* For M_PI, sqrt(), sin(), cos(), atan2(), log(), pow().*/
#include <float.h> /* For DBL_MAX. */

#include <Assertions.h> /* For PRE*(), POST*(), CHECK*(), IN_RANGE().     */
#include <Utilities.h>  /* For NEW().                                     */
#include <Lambert.h>    /* For public interface.                          */

/*================================== TYPES ==================================*/

struct LambertPrivate {
  double majorSemiaxis;    /* Mean equitorial radius in meters 6370000.0. */
  double minorSemiaxis;    /* Mean polar      radius in meters 6370000.0. */
  double lowerLatitude;    /* Lower tangent in degrees, e.g., 30.0.       */
  double upperLatitude;    /* Upper tangent in degrees, e.g., 60.0.       */
  double centralLongitude; /* Projects to zero, e.g., -100.0 degrees.     */
  double centralLatitude;  /* Projects to zero, e.g., 40.0 degrees.       */
  double falseEasting;     /* Skew offset in meters, e.g., 0.0.           */
  double falseNorthing;    /* Skew offset in meters, e.g., 0.0.           */
  double eccentricity;     /* Of ellipsoid approximation of planet.       */
  double lambda0;          /* Central longitude in radians.               */
  double rho0;             /* See USGS PROJ Library.                      */
  double n;                /* See USGS PROJ Library.                      */
  double c;                /* See USGS PROJ Library.                      */
};

/*=========================== FORWARD DECLARATIONS ==========================*/

/* Commands: */

static void free__( Lambert* self );

static void setEllipsoid( Lambert* self,
                          double majorSemiaxis, double minorSemiaxis );

static void setFalseEasting( Lambert* self, double falseEasting );

static void setFalseNorthing( Lambert* self, double falseNorthing );

static void project( const Lambert* self, double longitude, double latitude,
                     double* x, double* y );

static void unproject( const Lambert* self, double x, double y,
                       double* longitude, double* latitude );

/* Queries: */

static int invariant( const Lambert* self );
static int equal( const Lambert* self, const Lambert* other );
static Lambert* clone( const Lambert* self );

static void ellipsoid( const Lambert* self,
                       double* majorSemiaxis, double* minorSemiaxis );

static double falseEasting( const Lambert* self );
static double falseNorthing( const Lambert* self );
static double lowerLatitude( const Lambert* self );
static double upperLatitude( const Lambert* self );
static double centralLongitude( const Lambert* self );
static double centralLatitude( const Lambert* self );
static const char* name( const Lambert* self );

/* Helpers: */

static void assignMembers( Lambert* self );
static int hasMembers( const Lambert* self );
static void computeDerivedTerms( Lambert* self );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: newLambert - Construct a Lambert projector.
INPUTS:  double majorSemiaxis    Mean equitorial radius in meters 6370000.0.
         double minorSemiaxis    Mean polar      radius in meters 6370000.0.
         double lowerLatitude    Lower tangent in degrees, e.g., 30.0.
         double upperLatitude    Upper tangent in degrees, e.g., 60.0.
         double centralLongitude Projects to zero, e.g., -100.0 degrees.
         double centralLatitude  Projects to zero, e.g., 40.0 degrees.
         double falseEasting     Skew offset in meters, e.g., 0.0.
         double falseNorthing    Skew offset in meters, e.g., 0.0.
RETURNS: Lambert* else 0 if failed and failureMessage() called.
******************************************************************************/

Lambert* newLambert( double newMajorSemiaxis,    double newMinorSemiaxis,
                     double newLowerLatitude,    double newUpperLatitude,
                     double newCentralLongitude, double newCentralLatitude,
                     double newFalseEasting,     double newFalseNorthing ) {

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
  LambertPrivate* data = NEW( LambertPrivate, 1 );

  if ( data ) {
    result = NEW( Lambert, 1 );

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
INPUTS:  double newMajorSemiaxis  Mean equitorial radius in meters.
         double newMinorSemiaxis  Mean polar      radius in meters.
OUTPUTS: Lambert* self  Object to update.
******************************************************************************/

static void setEllipsoid( Lambert* self,
                          double newMajorSemiaxis, double newMinorSemiaxis ) {

  PRE2( self, isValidEllipsoid( newMajorSemiaxis, newMinorSemiaxis ) );
  self->data->majorSemiaxis = newMajorSemiaxis;
  self->data->minorSemiaxis = newMinorSemiaxis;
  computeDerivedTerms( self );
}



/******************************************************************************
PURPOSE: setFalseEasting - Set the projected x offset in meters.
INPUTS:  double newFalseEasting  Projected x offset in meters.
OUTPUTS: Lambert* self  Object to update.
******************************************************************************/

static void setFalseEasting( Lambert* self, double newFalseEasting ) {

  PRE2( self, ! isNan( newFalseEasting ) );
  self->data->falseEasting = newFalseEasting;
  POST( aboutEqual( self->falseEasting( self ), newFalseEasting ) );
}



/******************************************************************************
PURPOSE: setFalseNorthing - Set the projected y offset in meters.
INPUTS:  double newFalseNorthing  Projected y offset in meters.
OUTPUTS: Lambert* self  Object to update.
******************************************************************************/

static void setFalseNorthing( Lambert* self, double newFalseNorthing ) {

  PRE2( self, ! isNan( newFalseNorthing ) );
  self->data->falseNorthing = newFalseNorthing;
  POST( aboutEqual( self->falseNorthing( self ), newFalseNorthing ) );
}



/******************************************************************************
PURPOSE: project - Project a point.
INPUTS:  const Lambert* self   Projector.
         double longitude  E.g., -78.7268.
         double latitude   E.g., 35.9611.
OUTPUTS: double* x         Projected longitude.
         double* y         Projected latitude.
******************************************************************************/

static void project( const Lambert* self, double longitude, double latitude,
                     double* x, double* y ) {

  PRE4( isValidLongitude( longitude ), isValidLatitude( latitude ), x, y );

  const LambertPrivate* const data = self->data;
  double lambda = radians( longitude );
  double phi    = radians( latitude  );
  double rho = 0.0;    /* See USGS PROJ Library. */
  double lambdaDelta = 0.0;  /* Radians from central latitude. */
  double nLambdaDelta = 0.0; /* n times lambdaDelta. */

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
         double x           X-coordinate of point to unproject.
         double y           X-coordinate of point to unproject.
OUTPUTS: double* longitude  Unprojected x.
         double* latitude   Unprojected y.
******************************************************************************/

static void unproject( const Lambert* self, double x, double y,
                       double* longitude, double* latitude ) {

  PRE4( ! isNan( x ), ! isNan( y ), longitude, latitude );

  const LambertPrivate* const data = self->data;
  double oneOverMajorSemiaxis = 1.0 / data->majorSemiaxis;
  double xp = ( x - data->falseEasting  ) * oneOverMajorSemiaxis;
  double yp = ( y - data->falseNorthing ) * oneOverMajorSemiaxis;
  double ypDelta = data->rho0 - yp;/*Dist from yp to central lat (in radians)*/
  double rho     = hypot( xp, ypDelta );
  double lambda  = 0.0;       /* Radians of longitude. */
  double phi     = PI_OVER_2; /* Radians of latitude.  */

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
RETURNS: int 1 if valid, else 0.
NOTE:    If this query ever returns 0 then there is a defect in the code.
******************************************************************************/

static int invariant( const Lambert* self ) {

  const LambertPrivate* const data = self ? self->data : 0;
  const int result =
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
RETURNS: int 1 if equal, else 0.
******************************************************************************/

static int equal( const Lambert* self, const Lambert* other ) {

  PRE2( other, other->invariant( other ) );

  const LambertPrivate* const data = self->data;
  const LambertPrivate* const otherData = other->data;
  const int result =
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
OUTPUTS: double* majorSemiaxis   Mean equitorial radius.
OUTPUTS: double* minorSemiaxis   Mean polar      radius.
******************************************************************************/

static void ellipsoid( const Lambert* self,
                       double* majorSemiaxis, double* minorSemiaxis ) {

  PRE3( self, majorSemiaxis, minorSemiaxis );

  const LambertPrivate* const data = self->data;
  *majorSemiaxis = data->majorSemiaxis;
  *minorSemiaxis = data->majorSemiaxis;

  POST( isValidEllipsoid( *majorSemiaxis, *minorSemiaxis ) );
}



/******************************************************************************
PURPOSE: falseEasting - Projected x offset in meters.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: double result  Projected x offset in meters.
******************************************************************************/

static double falseEasting( const Lambert* self ) {

  PRE( self );

  const double result = self->data->falseEasting;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: falseNorthing - Projected y offset in meters.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: double result  Projected y offset in meters.
******************************************************************************/

static double falseNorthing( const Lambert* self ) {

  PRE( self );

  const double result = self->data->falseNorthing;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLongitude - Longitude of center of projection.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: double result longitude of center of projection.
******************************************************************************/

static double centralLongitude( const Lambert* self ) {

  PRE( self );

  const double result = self->data->centralLongitude;

  POST( isValidLongitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLatitude - Latitude of center of projection.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: double result latitude of center of projection.
******************************************************************************/

static double centralLatitude( const Lambert* self ) {

  PRE( self );

  const double result = self->data->centralLatitude;

  POST( isValidLatitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: lowerLatitude - Lower latitude of secant plane.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: double result lower latitude of secant plane.
******************************************************************************/

static double lowerLatitude( const Lambert* self ) {

  PRE( self );

  const double result = self->data->lowerLatitude;

  POST( isValidLatitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: upperLatitude - Upper latitude of secant plane.
INPUTS:  const Lambert* self   Projector to query.
RETURNS: double result upper latitude of secant plane.
******************************************************************************/

static double upperLatitude( const Lambert* self ) {

  PRE( self );

  const double result = self->data->upperLatitude;

  POST2( isValidLatitude( result ), result >= self->lowerLatitude( self ) );
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
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int hasMembers( const Lambert* self ) {

  PRE0( self );

  const int result =
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
  const double eccentricity0 = data->majorSemiaxis == data->minorSemiaxis ? 0.0 :
    safeQuotient( sqrt( safeDifference( SQUARE( data->majorSemiaxis ),
                                        SQUARE( data->minorSemiaxis ) ) ),
                  data->majorSemiaxis );

  const double eccentricity1 = eccentricity0 > 1.0 ? 1.0 : eccentricity0;
  const double eccentricitySquared = SQUARE( eccentricity1 );
  const double phi0 = radians( data->centralLatitude );
  const double phi1 = radians( data->lowerLatitude );
  const double phi2 = radians( data->upperLatitude );
  const double sinePhi1   = sin( phi1 );
  const double cosinePhi1 = cos( phi1 );
  const double sinePhi2   = sin( phi2 );
  const double cosinePhi2 = cos( phi2 );
  /* Are lower/upperLatitude about equal? */
  const int isTangent = phi1 + PROJECTION_TOLERANCE >= phi2;

  data->eccentricity = eccentricity1;
  data->lambda0 = radians( data->centralLongitude );
  data->n = sinePhi1;

  if ( eccentricitySquared != 0.0 ) { /* Ellipsoid planet: */
    const double m1  = msfn( sinePhi1, cosinePhi1, eccentricitySquared );
    const double ml1 = tsfn( phi1, sinePhi1, data->eccentricity );

    if ( ! isTangent ) { /* Secant form: */
      const double numerator =
        log( m1 / msfn( sinePhi2, cosinePhi2, eccentricitySquared ) );
      const double denominator =
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
    const double denominator = tan( PI_OVER_4 + 0.5 * phi1 );

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




/******************************************************************************
PURPOSE: Stereographic.c - Define Stereographic Projectors ADT.
NOTES:   See Projection.h and Stereographic.h. Formulations from USGS PROJ Lib.
HISTORY: 2004-10-01 plessel.todd@epa.gov Created based on C++ version.
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h> /* For stderr, fprintf(). */
#endif

#include <math.h>  /* For M_PI, sqrt(), sin(), cos(), atan2(), log(), pow().*/
#include <float.h> /* For DBL_MAX. */

#include <Assertions.h>    /* For PRE*(), POST*(), CHECK*(), IN_RANGE().     */
#include <Utilities.h>     /* For NEW().                                     */
#include <Stereographic.h> /* For public interface.                          */

/*================================== TYPES ==================================*/

enum { NORTH_POLE, SOUTH_POLE, EQUITORIAL, OBLIQUE }; /* Projection subtype. */
#define IS_VALID_SUBTYPE( subtype ) IN_RANGE( subtype, NORTH_POLE, OBLIQUE )

struct StereographicPrivate {
  double majorSemiaxis;    /* Mean equitorial radius in meters 6370000.0. */
  double minorSemiaxis;    /* Mean polar      radius in meters 6370000.0. */
  double falseEasting;     /* Skew offset in meters, e.g., 0.0.           */
  double falseNorthing;    /* Skew offset in meters, e.g., 0.0.           */
  double centralLongitude; /* Projects to zero, e.g., -98.0 degrees.     */
  double centralLatitude;  /* Projects to zero, e.g., 90.0 degrees.       */
  double secantLatitude;   /* Secant in degrees, e.g., 45.0.             */
  double eccentricity;     /* Of ellipsoid approximation of planet.       */
  double lambda0;          /* Central longitude in radians.               */
  double phi0;             /* Secant latitude in radians.                */
  double sineX1;           /* See USGS PROJ Library.                      */
  double cosineX1;         /* See USGS PROJ Library.                      */
  double akm1;             /* See USGS PROJ Library.                      */
  double projectedCenterX; /* Projected central longitude.                */
  double projectedCenterY; /* Projected central latitude.                 */
  int subtype;       /* EQUITORIAL, etc.                            */
  int initialized;   /* Initialized this object?                    */
};

/*=========================== FORWARD DECLARATIONS ==========================*/

/* Commands: */

static void free__( Stereographic* self );

static void setEllipsoid( Stereographic* self,
                          double majorSemiaxis, double minorSemiaxis );

static void setFalseEasting( Stereographic* self, double falseEasting );

static void setFalseNorthing( Stereographic* self, double falseNorthing );

static void project( const Stereographic* self, double longitude, double latitude,
                     double* x, double* y );

static void unproject( const Stereographic* self, double x, double y,
                       double* longitude, double* latitude );

/* Queries: */

static int invariant( const Stereographic* self );
static int equal( const Stereographic* self, const Stereographic* other );
static Stereographic* clone( const Stereographic* self );

static void ellipsoid( const Stereographic* self,
                       double* majorSemiaxis, double* minorSemiaxis );

static double falseEasting( const Stereographic* self );
static double falseNorthing( const Stereographic* self );
static double centralLongitude( const Stereographic* self );
static double centralLatitude( const Stereographic* self );
static double secantLatitude( const Stereographic* self );
static const char* name( const Stereographic* self );

/* Helpers: */

static void assignMembers( Stereographic* self );
static int hasMembers( const Stereographic* self );
static void computeDerivedTerms( Stereographic* self );

static void projectEllipsoid( const StereographicPrivate* data,
                              double phi,
                              double sineLambda,
                              double cosineLambda,
                              double sinePhi,
                              double* x, double* y );

static void projectSphere( const StereographicPrivate* data,
                           double phi,
                           double sineLambda,
                           double cosineLambda,
                           double sinePhi,
                           double* x, double* y );

static void unprojectEllipsoid( const StereographicPrivate* data,
                                double xp, double yp, double rho,
                                double* lambda, double* phi );

static void unprojectSphere( const StereographicPrivate* data,
                             double xp, double yp, double rho,
                             double* lambda, double* phi );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: newStereographic - Construct a Stereographic projector.
INPUTS:  double majorSemiaxis    Mean equitorial radius in meters 6370000.0.
         double minorSemiaxis    Mean polar      radius in meters 6370000.0.
         double centralLongitude Projects to zero, e.g., -98.0 degrees.
         double centralLatitude  Projects to zero, e.g., 90.0 degrees.
         double secantLatitude   Secant in degrees, e.g., 45.0.
         double falseEasting     Skew offset in meters, e.g., 0.0.
         double falseNorthing    Skew offset in meters, e.g., 0.0.
RETURNS: Stereographic* else 0 if failed and failureMessage() called.
******************************************************************************/

Stereographic* newStereographic( double newMajorSemiaxis,
                                 double newMinorSemiaxis,
                                 double newCentralLongitude,
                                 double newCentralLatitude,
                                 double newSecantLatitude,
                                 double newFalseEasting,
                                 double newFalseNorthing ) {

  PRE06( isValidEllipsoid( newMajorSemiaxis, newMinorSemiaxis ),
         isValidLongitude( newCentralLongitude ),
         isValidLatitude( newCentralLatitude ),
         isValidLatitude( newSecantLatitude ),
         ! isNan( newFalseEasting ),
         ! isNan( newFalseNorthing ) );

  Stereographic* result = 0;
  StereographicPrivate* data = NEW( StereographicPrivate, 1 );

  if ( data ) {
    result = NEW( Stereographic, 1 );

    if ( ! result ) {
      FREE( data );
    } else {
      data->majorSemiaxis    = newMajorSemiaxis;
      data->minorSemiaxis    = newMinorSemiaxis;
      data->falseEasting     = newFalseEasting;
      data->falseNorthing    = newFalseNorthing;
      data->centralLongitude = newCentralLongitude;
      data->centralLatitude  = newCentralLatitude;
      data->secantLatitude   = newSecantLatitude;
      result->data = data;
      assignMembers( result );
      computeDerivedTerms( result );
      result->data->initialized = 1;
    }
  }

  POST0( IMPLIES( result,
                  AND3( result->invariant( result ),
                        result->falseEasting( result ) == newFalseEasting,
                        result->falseNorthing( result) == newFalseNorthing)));
  return result;
}



/******************************************************************************
PURPOSE: free - Destruct a Stereographic.
INPUTS:  Stereographic* self  Object to destruct.
NOTE:    Use FREE_OBJECT( object ) instead since it zeros argument.
******************************************************************************/

static void free__( Stereographic* self ) {
  PRE( self );
  FREE_ZERO( self->data );
  FREE_ZERO( self );
  POST0( self == 0 );
}



/******************************************************************************
PURPOSE: setEllipsoid - Set the ellipsoid approximation of planet.
INPUTS:  double     newMajorSemiaxis  Mean equitorial radius in meters.
         double     newMinorSemiaxis  Mean polar      radius in meters.
OUTPUTS: Stereographic* self  Object to update.
******************************************************************************/

static void setEllipsoid( Stereographic* self,
                          double newMajorSemiaxis, double newMinorSemiaxis ) {

  PRE2( self, isValidEllipsoid( newMajorSemiaxis, newMinorSemiaxis ) );
  self->data->majorSemiaxis = newMajorSemiaxis;
  self->data->minorSemiaxis = newMinorSemiaxis;
  computeDerivedTerms( self );
}



/******************************************************************************
PURPOSE: setFalseEasting - Set the projected x offset in meters.
INPUTS:  double newFalseEasting  Projected x offset in meters.
OUTPUTS: Stereographic* self  Object to update.
******************************************************************************/

static void setFalseEasting( Stereographic* self, double newFalseEasting ) {

  PRE2( self, ! isNan( newFalseEasting ) );
  self->data->falseEasting = newFalseEasting;
  POST( aboutEqual( self->falseEasting( self ), newFalseEasting ) );
}



/******************************************************************************
PURPOSE: setFalseNorthing - Set the projected y offset in meters.
INPUTS:  double newFalseNorthing  Projected y offset in meters.
OUTPUTS: Stereographic* self  Object to update.
******************************************************************************/

static void setFalseNorthing( Stereographic* self, double newFalseNorthing ) {

  PRE2( self, ! isNan( newFalseNorthing ) );
  self->data->falseNorthing = newFalseNorthing;
  POST( aboutEqual( self->falseNorthing( self ), newFalseNorthing ) );
}



/******************************************************************************
PURPOSE: project - Project a point.
INPUTS:  project(Stereographic* self   Projector.
         double longitude  E.g., -78.7268.
         double latitude   E.g., 35.9611.
OUTPUTS: double* x         Projected longitude.
         double* y         Projected latitude.
******************************************************************************/

static void project( const Stereographic* self,
                     double longitude, double latitude,
                     double* x, double* y ) {

  PRE4( isValidLongitude( longitude ), isValidLatitude( latitude ), x, y );

  const StereographicPrivate* const data = self->data;
  const double lambda0 = radians( longitude );
  const double phi0    = radians( latitude );
  const double lambda  =
    OR2( ! data->initialized,
         IN_RANGE( lambda0, -M_PI + PROJECTION_TOLERANCE,
                             M_PI - PROJECTION_TOLERANCE ) ) ?
      lambda0 - data->lambda0
    : SIGN( lambda0 ) * ( M_PI - PROJECTION_TOLERANCE ) - data->lambda0;

  const double phi =
    OR2( ! data->initialized,
         IN_RANGE( phi0, -PI_OVER_2 + PROJECTION_TOLERANCE,
                          PI_OVER_2 - PROJECTION_TOLERANCE ) ) ?
      phi0
    : SIGN( phi0 ) * ( PI_OVER_2 - PROJECTION_TOLERANCE );

  const double sineLambda   = sin( lambda );
  const double cosineLambda = cos( lambda );
  const double sinePhi      = sin( phi );
  CHECK( IS_VALID_SUBTYPE( data->subtype ) );
  *x = *y = 0.0;
  DEBUG( fprintf( stderr, "init = %d\n", data->initialized ); )
  DEBUG( if ( ! data->initialized ) {
           fprintf( stderr, "lambda0 = %lf, phi0 = %lf, akm1 = %lf, "
                    "pCX = %lf, pCY = %lf\n",
                    data->lambda0, data->phi0, data->akm1,
                    data->projectedCenterX, data->projectedCenterY );
          } )
  DEBUG( fprintf( stderr, "p(lon, lat) = (%28.20lf, %28.20lf)\n",
                  longitude, latitude ); )
  DEBUG( fprintf( stderr, "p lambda = %28.20lf, phi = %28.20lf\n",
                  lambda, phi ); )

  if ( data->eccentricity ) { /* Ellipsoid: */
    projectEllipsoid( data, phi, sineLambda, cosineLambda, sinePhi, x, y );
  } else { /* Sphere: */
    projectSphere( data, phi, sineLambda, cosineLambda, sinePhi, x, y );
  }

  *x = *x * data->majorSemiaxis + data->falseEasting;
  *y = *y * data->majorSemiaxis + data->falseNorthing;
  DEBUG( fprintf( stderr, "(x, y) = (%28.20lf, %28.20lf)\n", *x, *y ); )
  POST2( ! isNan( *x ), ! isNan( *y ) );
}



/******************************************************************************
PURPOSE: unproject - Unproject a point.
INPUTS:  const Stereographic* self    Projector.
         double x           X-coordinate of point to unproject.
         double y           X-coordinate of point to unproject.
OUTPUTS: double* longitude  Unprojected x.
         double* latitude   Unprojected y.
******************************************************************************/

static void unproject( const Stereographic* self, double x, double y,
                       double* longitude, double* latitude ) {

  PRE4( ! isNan( x ), ! isNan( y ), longitude, latitude );

  const StereographicPrivate* const data = self->data;
  double oneOverMajorSemiaxis = 1.0 / data->majorSemiaxis;
  double xp = ( x + data->projectedCenterX - data->falseEasting ) *
            oneOverMajorSemiaxis;
  double yp = ( y + data->projectedCenterY - data->falseNorthing ) *
            oneOverMajorSemiaxis;
  const double rho = hypot( xp, yp );
  double lambda = 0.0; /* Radians of longitude. */
  double phi    = 0.0; /* Radians of latitude.  */
  CHECK( IS_VALID_SUBTYPE( data->subtype ) );
  DEBUG( fprintf( stderr, "u(x, y) = (%28.20lf, %28.20lf)\n", x, y ); )

  if ( data->eccentricity ) { /* Ellipsoid: */
    unprojectEllipsoid( data, xp, yp, rho, &lambda, &phi );
  } else { /* Sphere: */
    unprojectSphere( data, xp, yp, rho, &lambda, &phi );
  }

  lambda += data->lambda0;
  *longitude = degrees( lambda );
  *latitude  = degrees( phi );
  CHECK( fabs( *longitude ) < DBL_MAX );

  while ( *longitude < -180.0 ) {
    *longitude += 360.0;
  }

  while ( *longitude > 180.0 ) {
    *longitude -= 360.0;
  }

  DEBUG( fprintf( stderr, "(lon, lat) = (%28.20lf, %28.20lf)\n",
                  *longitude, *latitude ); )
  POST2( isValidLongitude( *longitude ), isValidLatitude( *latitude ) );
}



/******************************************************************************
PURPOSE: invariant - Class invariant.
INPUTS:  const Stereographic* self  Projector to check.
RETURNS: int 1 if valid, else 0.
NOTE:    If this query ever returns 0 then there is a defect in the code.
******************************************************************************/

static int invariant( const Stereographic* self ) {

  const StereographicPrivate* const data = self ? self->data : 0;
  const int result =
    AND10( self,
           hasMembers( self ),
           self->data,
           data == self->data,
           isValidEllipsoid( data->majorSemiaxis, data->minorSemiaxis ),
           ! isNan( data->falseEasting ),
           ! isNan( data->falseNorthing ),
           isValidLongitude( data->centralLongitude ),
           isValidLatitude( data->centralLatitude ),
           isValidLatitude( data->secantLatitude ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: equal - Is self functionally equivalent to other?
INPUTS:  const Stereographic* self   Projector to compare.
         const Stereographic* other  Projector to compare.
RETURNS: int 1 if equal, else 0.
******************************************************************************/

static int equal( const Stereographic* self, const Stereographic* other ) {

  PRE2( other, other->invariant( other ) );

  const StereographicPrivate* const data = self->data;
  const StereographicPrivate* const otherData = other->data;
  const int result =
    AND7( aboutEqual( data->majorSemiaxis,    otherData->majorSemiaxis ),
          aboutEqual( data->minorSemiaxis,    otherData->minorSemiaxis ),
          aboutEqual( data->falseEasting,     otherData->falseEasting ),
          aboutEqual( data->falseNorthing,    otherData->falseNorthing ),
          aboutEqual( data->centralLongitude, otherData->centralLongitude ),
          aboutEqual( data->centralLatitude,  otherData->centralLatitude ),
          aboutEqual( data->secantLatitude,   otherData->secantLatitude ) );

  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: clone - Yield a clone - a NEW functionally equivalent Stereographic.
INPUTS:  const Stereographic* self   Projector to clone.
RETURNS: Stereographic* result clone of self.
******************************************************************************/

static Stereographic* clone( const Stereographic* self ) {

  PRE( self );

  const StereographicPrivate* const data = self->data;
  Stereographic* const result =
    newStereographic( data->majorSemiaxis,    data->minorSemiaxis,
                      data->centralLongitude, data->centralLatitude,
                      data->secantLatitude,
                      data->falseEasting,     data->falseNorthing );

  POST( IMPLIES( result,
                 AND2( result->invariant( result ),
                       result->equal( result, self ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: ellipsoid - The ellipsoid axes of the planet approximation.
INPUTS:  const Stereographic* self   Projector to query.
OUTPUTS: double* majorSemiaxis   Mean equitorial radius.
OUTPUTS: double* minorSemiaxis   Mean polar      radius.
******************************************************************************/

static void ellipsoid( const Stereographic* self,
                       double* majorSemiaxis, double* minorSemiaxis ) {

  PRE3( self, majorSemiaxis, minorSemiaxis );

  const StereographicPrivate* const data = self->data;
  *majorSemiaxis = data->majorSemiaxis;
  *minorSemiaxis = data->majorSemiaxis;

  POST( isValidEllipsoid( *majorSemiaxis, *minorSemiaxis ) );
}



/******************************************************************************
PURPOSE: falseEasting - Projected x offset in meters.
INPUTS:  const Stereographic* self   Projector to query.
RETURNS: double result  Projected x offset in meters.
******************************************************************************/

static double falseEasting( const Stereographic* self ) {

  PRE( self );

  const double result = self->data->falseEasting;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: falseNorthing - Projected y offset in meters.
INPUTS:  const Stereographic* self   Projector to query.
RETURNS: double result  Projected y offset in meters.
******************************************************************************/

static double falseNorthing( const Stereographic* self ) {

  PRE( self );

  const double result = self->data->falseNorthing;

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLongitude - Longitude of center of projection.
INPUTS:  const Stereographic* self   Projector to query.
RETURNS: double result longitude of center of projection.
******************************************************************************/

static double centralLongitude( const Stereographic* self ) {

  PRE( self );

  const double result = self->data->centralLongitude;

  POST( isValidLongitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: centralLatitude - Latitude of center of projection.
INPUTS:  const Stereographic* self   Projector to query.
RETURNS: double result latitude of center of projection.
******************************************************************************/

static double centralLatitude( const Stereographic* self ) {

  PRE( self );

  const double result = self->data->centralLatitude;

  POST( isValidLatitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: secantLatitude - Latitude of secant plane.
INPUTS:  const Stereographic* self   Projector to query.
RETURNS: double result lower latitude of secant plane.
******************************************************************************/

static double secantLatitude( const Stereographic* self ) {

  PRE( self );

  const double result = self->data->secantLatitude;

  POST( isValidLatitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: name - Name of projection.
INPUTS:  const Stereographic* self   Projector to query.
RETURNS: const char*  Name of projection: "Stereographic".
******************************************************************************/

static const char* name( const Stereographic* self ) {

  PRE( self );

  const char* const result = "Stereographic";

  POST( ! strcmp( result, "Stereographic" ) );
  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: assignMembers - Assign pointers to member functions.
OUTPUTS: Stereographic* self  Stereographic to initialize.
******************************************************************************/

static void assignMembers( Stereographic* self ) {

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
  self->secantLatitude   = secantLatitude;

  POST0( hasMembers( self ) );
}



/******************************************************************************
PURPOSE: hasMembers - Does self have all pointers to member functions?
INPUTS:  const Stereographic* self  Stereographic to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int hasMembers( const Stereographic* self ) {

  PRE0( self );

  const int result =
    AND16( self->free             == free__,
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
           self->secantLatitude   == secantLatitude );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: computeDerivedTerms - Compute trigonometry terms independent of
         longitude/latitude of projection point.
OUTPUTS: Stereographic* self  Stereographic to (re)initialize.
******************************************************************************/

static void computeDerivedTerms( Stereographic* self ) {

  PRE0( self );

  StereographicPrivate* const data = self->data;
  const double eccentricity0 = data->majorSemiaxis == data->minorSemiaxis ? 0.0 :
    safeQuotient( sqrt( safeDifference( SQUARE( data->majorSemiaxis ),
                                        SQUARE( data->minorSemiaxis ) ) ),
                  data->majorSemiaxis );

  const double eccentricity1 = eccentricity0 > 1.0 ? 1.0 : eccentricity0;
  const double phits = fabs( radians( data->secantLatitude ) );
  const double k0 = ( 1.0 + sin( phits ) ) * 0.5;
  double t = 0.0;

  data->eccentricity = eccentricity1;
  data->lambda0 = radians( data->centralLongitude );
  data->phi0    = radians( data->centralLatitude  );
  t = fabs( data->phi0 );

  if ( fabs( t - PI_OVER_2 ) < PROJECTION_TOLERANCE ) {
    data->subtype = data->phi0 < 0.0 ? SOUTH_POLE : NORTH_POLE;
  } else {
    data->subtype = t > PROJECTION_TOLERANCE ? OBLIQUE : EQUITORIAL;
  }

  if ( data->eccentricity != 0.0 ) { /* Ellipsoid planet: */

    switch ( data->subtype ) {
    case EQUITORIAL:
      data->akm1 = k0 + k0;
      break;
    case OBLIQUE:
      {
        double X = 0.0;
        t = sin( data->phi0 );
        X = atan( ssfn( data->phi0, t, data->eccentricity ) );
        X += X;
        X -= PI_OVER_2;
        t *= data->eccentricity;
        CHECK( SQUARE( t ) < 1.0 );
        data->akm1 = ( k0 + k0 ) * cos( data->phi0 ) / sqrt( 1.0 - SQUARE( t ) );
        data->sineX1 = sin( X );
        data->cosineX1 = cos( X );
      }
      break;
    default:
      CHECK( IN3( data->subtype, NORTH_POLE, SOUTH_POLE ) );

      if ( fabs( phits - PI_OVER_2 ) < PROJECTION_TOLERANCE ) {
        const double X  = 1.0 + data->eccentricity;
        const double X2 = 1.0 - data->eccentricity;
        CHECK( sqrt( pow( X, X ) * pow( X2, X2 ) ) > 0.0 );
        data->akm1 = ( k0 + k0 ) / sqrt( pow( X, X ) * pow( X2, X2 ) );
      } else {
        t = sin( phits );
        data->akm1 = cos( phits ) / tsfn( phits, t, data->eccentricity );
        t *= data->eccentricity;
        CHECK2( SQUARE( t ) < 1.0, sqrt( 1.0 - SQUARE( t ) ) > 0.0 );
        data->akm1 /= sqrt( 1.0 - SQUARE( t ) );
      }

      break;
    }
  } else { /* Sphere planet: */

    switch ( data->subtype ) {
    case EQUITORIAL:
      data->akm1 = k0 + k0;
      break;
    case OBLIQUE:
      data->sineX1   = sin( data->phi0 );
      data->cosineX1 = cos( data->phi0 );
      data->akm1 = k0 + k0;
      break;
    default:
      CHECK( IN3( data->subtype, NORTH_POLE, SOUTH_POLE ) );
      data->akm1 = fabs( phits - PI_OVER_2 ) >= PROJECTION_TOLERANCE ?
                     cos( phits ) / tan( PI_OVER_4 - 0.5 * phits )
                   : k0 + k0 ;
      DEBUG( fprintf( stderr, "phits = %lf, k0 = %lf, "
                      "data->eccentricity = %lf, "
                      "data->lambda0 = %lf, "
                      "data->phi0 = %lf, "
                      "data->sineX1 = %lf, "
                      "data->cosineX1 = %lf, "
                      "data->akm1 = %lf, "
                      "data->subtype = %d\n",
                      phits, k0,
                      data->eccentricity,
                      data->lambda0,
                      data->phi0,
                      data->sineX1,
                      data->cosineX1,
                      data->akm1,
                      data->subtype ); )

      break;
    }
  }

  project( self, data->centralLongitude, data->centralLatitude,
           &data->projectedCenterX, &data->projectedCenterY );

  POST010( ! isNan( data->eccentricity ),
           IN_RANGE( data->eccentricity, 0.0, 1.0 ),
           ! isNan( data->lambda0 ),
           ! isNan( data->phi0 ),
           ! isNan( data->sineX1 ),
           ! isNan( data->cosineX1 ),
           IN_RANGE( data->sineX1, -1.0, 1.0 ),
           IN_RANGE( data->cosineX1, -1.0, 1.0 ),
           ! isNan( data->akm1 ),
           IS_VALID_SUBTYPE( data->subtype ) );
}



/******************************************************************************
PURPOSE: projectEllipsoid - Project a point on an ellipsoid planet.
INPUTS:  const StereographicPrivate* data  Object data.
         double phi           Angle in radians of adjusted latitude to project.
         double sineLambda    Sine of lambda.
         double cosineLambda  Cosine of lambda.
         double sinePhi       Sine of phi.
OUTPUTS: double* x            Projected longitude.
         double* y            Projected latitude.
******************************************************************************/

static void projectEllipsoid( const StereographicPrivate* data,
                              double phi,
                              double sineLambda, double cosineLambda, double sinePhi,
                              double* x, double* y ) {

  PRE011( data, data->eccentricity,
          ! isNan( phi ),
          ! isNan( sineLambda ), ! isNan( cosineLambda ), ! isNan( sinePhi ),
          aboutEqual( sinePhi, sin( phi ) ),
          x, y, *x == 0.0, *y == 0.0 );

  switch ( data->subtype ) {
  case OBLIQUE: {
    const double X =
      2.0 * atan( ssfn( phi, sinePhi, data->eccentricity ) ) - PI_OVER_2;
    const double sineX   = sin( X );
    const double cosineX = cos( X );
    const double A = data->akm1 /
      ( data->cosineX1 * ( 1.0 + data->sineX1 * sineX +
                           data->cosineX1 * cosineX * cosineLambda ) );
    *x = A * cosineX;
    *y = A * ( data->cosineX1 * sineX -
               data->sineX1 * cosineX * cosineLambda );
    }

    break;
  case EQUITORIAL: {
    const double X =
      2.0 * atan( ssfn( phi, sinePhi, data->eccentricity ) ) - PI_OVER_2;
    const double sineX   = sin( X );
    const double cosineX = cos( X );
    const double A = 2.0 * data->akm1 / ( 1.0 + cosineX * cosineLambda );
    *x = A * cosineX;
    *y = A * sineX;
    }

    break;
  case SOUTH_POLE:
    *x = data->akm1 * tsfn( -phi, -sinePhi, data->eccentricity );
    *y = *x * cosineLambda;
    break;
  default:
    CHECK( data->subtype == NORTH_POLE );
    *x = data->akm1 * tsfn( phi, sinePhi, data->eccentricity );
    *y = -(*x) * cosineLambda;
   break;
  }

  *x *= sineLambda;

  POST02( ! isNan( *x ), ! isNan( *y ) );
}



/******************************************************************************
PURPOSE: projectSphere - Project a point on a sphere planet.
INPUTS:  const StereographicPrivate* data  Object data.
         double phi           Angle in radians of adjusted latitude to project.
         double sineLambda    Sine of lambda.
         double cosineLambda  Cosine of lambda.
         double sinePhi       Sine of phi.
OUTPUTS: double* x            Projected longitude.
         double* y            Projected latitude.
******************************************************************************/

static void projectSphere( const StereographicPrivate* data,
                           double phi,
                           double sineLambda, double cosineLambda, double sinePhi,
                           double* x, double* y ) {

  PRE011( data, data->eccentricity == 0.0,
          ! isNan( phi ),
          ! isNan( sineLambda ), ! isNan( cosineLambda ), ! isNan( sinePhi ),
          aboutEqual( sinePhi, sin( phi ) ),
          x, y, *x == 0.0, *y == 0.0 );

  switch ( data->subtype ) {
  case EQUITORIAL:

    {
      const double cosinePhi = cos( phi );
      *y = 1.0 + cosinePhi * cosineLambda;

      if ( *y != 0.0) {
        *y = data->akm1 / *y;
        *x = *y * cosinePhi * sineLambda;
        *y *= sinePhi;
      }
    }

    break;
  case OBLIQUE:

    {
      const double cosinePhi = cos( phi );
      *y = 1.0 +
           data->sineX1 * sinePhi + data->cosineX1 * cosinePhi * cosineLambda;

      if ( *y != 0.0) {
        *y = data->akm1 / *y;
        *x = *y * cosinePhi * sineLambda;
        *y *=
           data->cosineX1 * sinePhi - data->sineX1 * cosinePhi * cosineLambda;
      }
    }

    break;
  case NORTH_POLE:

    if ( fabs( phi - PI_OVER_2 ) >= PROJECTION_TOLERANCE ) {
      *y = data->akm1 * tan( PI_OVER_4 + 0.5 * -phi );
      *x = sineLambda * *y;
      *y *= -cosineLambda;
    }

    break;
  default:
    CHECK( data->subtype == SOUTH_POLE );

    if ( fabs( phi - PI_OVER_2 ) >= PROJECTION_TOLERANCE ) {
      *y = data->akm1 * tan( PI_OVER_4 + 0.5 * phi );
      *x = sineLambda * *y;
      *y *= cosineLambda;
    }

    DEBUG( fprintf( stderr, "ps phi = %lf, sineLambda = %lf, "
                    "cosineLambda = %lf, sinePhi = %lf, "
                    "data->akm1 = %lf, "
                    "*x = %lf, *y = %lf\n",
                    phi, sineLambda, cosineLambda, sinePhi,
                    data->akm1, *x, *y ); )

    break;
  }

  POST02( ! isNan( *x ), ! isNan( *y ) );
}



/******************************************************************************
PURPOSE: unprojectEllipsoid - Unproject a point on an ellipsoid planet.
INPUTS:  const StereographicData* data    Object data.
         double xp       Adjusted x-coordinate of point to unproject.
         double yp       Adjusted y-coordinate of point to unproject.
         double rho      Hypoteneuse of xp, yp.
OUTPUTS: double* lambda  Unprojected xp.
         double* phi     Unprojected yp.
******************************************************************************/

static void unprojectEllipsoid( const StereographicPrivate* data,
                                double xp, double yp, double rho,
                                double* lambda, double* phi ) {

  PRE09( data, data->eccentricity,
         ! isNan( xp ), ! isNan( yp ), ! isNan( rho ),
         lambda, phi, *lambda == 0.0, *phi == 0.0 );

  double sinePhi   = 0.0;
  double cosinePhi = 0.0;
  double phiL      = 0.0;
  double tp        = 0.0;
  double halfPi    = PI_OVER_2;
  double halfEccentricity = 0.5 * data->eccentricity;
  int iteration = 0;

  switch ( data->subtype ) {
  case EQUITORIAL:
    /* Fall through: */
  case OBLIQUE:
    tp = atan2( rho * data->cosineX1, data->akm1 );
    tp += tp;
    cosinePhi = cos( tp );
    sinePhi   = sin( tp );
    CHECK( rho != 0.0 );
    phiL = asin( cosinePhi * data->sineX1 +
                 ( yp * sinePhi * data->cosineX1 / rho ) );
    tp = tan( 0.5 * ( PI_OVER_2 + phiL ) );
    xp *= sinePhi;
    yp = rho * data->cosineX1 * cosinePhi - yp * data->sineX1 * sinePhi;
    break;
  case NORTH_POLE:
    yp = -yp;
    /* Fall through: */
  default: /* SOUTH_POLE: */
    CHECK( data->akm1 != 0.0 );
    tp = -rho / data->akm1;
    phiL = PI_OVER_2 - 2.0 * atan( tp );
    halfPi = -halfPi;
    halfEccentricity = -halfEccentricity;
    break;
  }

  do {
    sinePhi = data->eccentricity * sin( phiL );
    CHECK( sinePhi != 1.0 );
    *phi = 2.0 * atan( tp * pow( ( 1.0 + sinePhi ) / ( 1.0 - sinePhi ),
                                 halfEccentricity ) )
           - halfPi;

    if ( fabs( phiL - *phi ) < CONVERGENCE_TOLERANCE ) {

      if (data->subtype == SOUTH_POLE ) {
        *phi = -*phi;
      }

      iteration = MAXIMUM_ITERATIONS;
    }

    ++iteration;
    phiL = *phi;
  } while ( iteration < MAXIMUM_ITERATIONS );

  if ( iteration == MAXIMUM_ITERATIONS ) {
    *lambda = *phi = 0.0; /* Failed to converge! */
  } else {
    *lambda = IS_ZERO2( xp, yp ) ? 0.0 : atan2( xp, yp );
  }

  POST02( ! isNan( *lambda ), ! isNan( *phi ) );
}



/******************************************************************************
PURPOSE: unprojectSphere - Unproject a point on a sphere planet.
INPUTS:  const StereographicData* data    Object data.
         double xp           Adjusted x-coordinate of point to unproject.
         double yp           Adjusted y-coordinate of point to unproject.
         double rho          Hypoteneuse of xp, yp.
OUTPUTS: double* lambda      Unprojected xp.
         double* phi         Unprojected yp.
******************************************************************************/

static void unprojectSphere( const StereographicPrivate* data,
                             double xp, double yp, double rho,
                             double* lambda, double* phi ) {

  PRE09( data, data->eccentricity == 0.0,
         ! isNan( xp ), ! isNan( yp ), ! isNan( rho ),
         lambda, phi, *lambda == 0.0, *phi == 0.0 );

  double c = 2.0 * atan( rho / data->akm1 );
  const double cosineC = cos( c );

  DEBUG( fprintf( stderr, "akm1 = %lf, c = %lf\n", data->akm1, c ); )

  switch ( data->subtype ) {
  case EQUITORIAL:

    {
      const double sineC = sin( c );

      if ( fabs( rho ) > PROJECTION_TOLERANCE ) {
        *phi = asin( yp * sineC / rho );
      }

      if ( OR2( cosineC, xp ) ) {
        *lambda = atan2( xp * sineC, cosineC * rho );
      }
    }

    break;
  case OBLIQUE:

    {
      const double sineC = sin( c );

      if ( fabs( rho ) <= PROJECTION_TOLERANCE ) {
        *phi = data->phi0;
      } else {
        *phi = asin(cosineC * data->sineX1 + yp * sineC * data->cosineX1 /rho);
      }

      c = cosineC - data->sineX1 * sin( *phi );

      if ( OR2( c, xp ) ) {
        *lambda = atan2( xp * sineC * data->cosineX1, c * rho );
      }
    }

    break;
  case NORTH_POLE:
    yp = -yp;
    /* Fall through. */
  default:
    CHECK( IN3( data->subtype, NORTH_POLE, SOUTH_POLE ) );

    if ( fabs( rho ) <= PROJECTION_TOLERANCE ) {
      *phi = data->phi0;
    } else {
      *phi = asin( data->subtype == SOUTH_POLE ? -cosineC : cosineC );
    }

    *lambda = IS_ZERO2( xp, yp ) ? 0.0 : atan2( xp, yp );
    DEBUG( fprintf( stderr, "us rho = %lf, phi = %lf, *lambda = %lf,"
                    "xp = %lf, yp = %lf, cosineC = %lf\n",
                    rho, *phi, *lambda, xp, yp, cosineC ); )
    break;
  }

  POST02( ! isNan( *lambda ), ! isNan( *phi ) );
}





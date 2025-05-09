
/* lambert.c - Lambert Conformal Conic projections derived from USGS PROJ. */

#include <math.h>  /* For M_PI, sqrt(), sin(), cos(), atan2(), log(), pow().*/
#include <float.h> /* For DBL_MAX. */

#include <Assertions.h>  /* For PRE*(), POST*(), CHECK*(), IN_RANGE(). */
#include <projections.h> /* For is_valid_latitude(), etc. */
#include <lambert.h>     /* For public interface. */

/*============================ PRIVATE VARIABLES ============================*/

/* Inputs: */

static double major_semiaxis;    /* Mean equitorial radius in meters 6370997*/
static double minor_semiaxis;    /* Mean polar      radius in meters 6370997*/
static double lower_latitude;    /* Lower tangent in degrees, e.g., 30.0.   */
static double upper_latitude;    /* Upper tangent in degrees, e.g., 60.0.   */
static double central_latitude;  /* Projects to zero, e.g., 49.0 degrees.   */
static double central_longitude; /* Projects to zero, e.g., -121.0 degrees. */
static double false_easting;     /* Skew offset in meters, e.g., 0.0.       */
static double false_northing;    /* Skew offset in meters, e.g., 0.0.       */

/* Derived terms: */

static double eccentricity;      /* Of ellipsoid approximation of planet.   */
static double lambda0;           /* Central longitude in radians.           */
static double rho0;              /* See USGS PROJ Library.                  */
static double n;                 /* See USGS PROJ Library.                  */
static double c;                 /* See USGS PROJ Library.                  */
static int initialized;          /* Has client called initialize_lambert()? */

/*=========================== FORWARD DECLARATIONS ==========================*/

static void recompute_derived_terms( void );
static double phi2_iterate( double ts, double the_eccentricity );

/*============================= PUBLIC FUNCTIONS ============================*/


void initialize_lambert(
  double new_major_semiaxis,   double new_minor_semiaxis,
  double new_lower_latitude,   double new_upper_latitude,
  double new_central_latitude, double new_central_longitude,
  double new_false_easting,    double new_false_northing ) {

  PRE012( is_valid_ellipsoid( new_major_semiaxis, new_minor_semiaxis ),
          is_valid_latitude( new_lower_latitude ),
          is_valid_latitude( new_upper_latitude ),
          is_valid_latitude( new_central_latitude ),
          is_valid_longitude( new_central_longitude ),
          new_lower_latitude <= new_upper_latitude,
          SIGN( new_lower_latitude ) == SIGN( new_upper_latitude ),
          IN_RANGE( new_lower_latitude,
                    SIGN( new_lower_latitude ) * 1.0,
                    SIGN( new_lower_latitude ) * 89.0 ),
          IN_RANGE( new_upper_latitude,
                    SIGN( new_upper_latitude ) * 1.0,
                    SIGN( new_upper_latitude ) * 89.0 ),
          IN_RANGE( new_central_latitude, -89.0, 89.0 ),
          ! is_nan( new_false_easting ),
          ! is_nan( new_false_northing ) );

  major_semiaxis    = new_major_semiaxis;
  minor_semiaxis    = new_minor_semiaxis;
  lower_latitude    = new_lower_latitude;
  upper_latitude    = new_upper_latitude;
  central_latitude  = new_central_latitude;
  central_longitude = new_central_longitude;
  false_easting     = new_false_easting;
  false_northing    = new_false_northing;

  recompute_derived_terms();
}



void project_lambert( double longitude, double latitude,
                      double* x, double* y ) {

  PRE05( is_valid_longitude( longitude ), is_valid_latitude( latitude ),
         x, y, initialized );

  double lambda = to_radians( longitude );
  double phi    = to_radians( latitude  );
  double rho;            /* See USGS PROJ Library. */
  double lambda_delta;   /* Radians from central latitude. */
  double n_lambda_delta; /* n times lambda_delta. */

  /*
   * If phi is too near a pole tweak it so that projecting
   * succeeds and unprojecting yields original longitude
   * (instead of central longitude).
   */

  if ( ! IN_RANGE( phi, -PI_OVER_2 + TOLERANCE, PI_OVER_2 - TOLERANCE ) ) {
    phi = phi + TOLERANCE * -SIGN( phi );
  }

  rho = c * pow( tsfn( phi, sin( phi ), eccentricity ), n );

  /*
   * If lambda is too near +/-180 longitude tweak it so that projecting
   * succeeds and unprojecting yields original longitude
   * (instead of central longitude).
   */

  if ( ! IN_RANGE( phi, -M_PI + TOLERANCE, M_PI - TOLERANCE ) ) {
    lambda = lambda + SQUARE( TOLERANCE ) * -SIGN( lambda );
  }

  for ( lambda_delta = lambda - lambda0; fabs( lambda_delta ) > M_PI; ) {

    if ( lambda_delta < 0.0 ) {
      lambda_delta = lambda_delta + M_PI + M_PI;
    } else {
      lambda_delta = lambda_delta - M_PI - M_PI;
    }
  }

  n_lambda_delta = n * lambda_delta;
  *x = rho          * sin( n_lambda_delta ) * major_semiaxis + false_easting;
  *y = ( rho0 - rho * cos( n_lambda_delta)) * major_semiaxis + false_northing;

  POST02( ! is_nan( *x ), ! is_nan( *y ) );
}



void unproject_lambert( double x, double y,
                        double* longitude, double* latitude ) {

  PRE05( ! is_nan( x ), ! is_nan( y ), longitude, latitude, initialized );
  double one_over_major_semiaxis = 1.0 / major_semiaxis;
  double xp = ( x - false_easting  ) * one_over_major_semiaxis;
  double yp = ( y - false_northing ) * one_over_major_semiaxis;
  double yp_delta = rho0 - yp;/*Dist from yp to central latitude (in radians)*/
  double rho      = hypot( xp, yp_delta );
  double lambda   = 0.0;       /* Radians of longitude. */
  double phi      = PI_OVER_2; /* Radians of latitude.  */

  if ( rho != 0.0 ) {

    if ( n < 0.0 ) {
      rho = -rho;
      xp = -xp;
      yp_delta = -yp_delta;
    }

    CHECK3( c != 0.0, n != 0.0, rho != 0.0 );

    if ( eccentricity == 0.0 ) { /* Sphere: */
      phi = 2.0 * atan( pow( c / rho, 1.0 / n ) ) - PI_OVER_2;
    } else { /* Ellipsoid: */
      phi = phi2_iterate( pow( rho / c, 1.0 / n ), eccentricity );
    }

    lambda = atan2( xp, yp_delta ) / n;
  } else if ( n < 0.0 ) {
    phi = -PI_OVER_2;
  }

  lambda += lambda0;
  *longitude = to_degrees( lambda );
  *latitude = to_degrees( phi );
  CHECK( fabs( *longitude ) < DBL_MAX );

  while ( *longitude < -180.0 ) {
    *longitude += 360.0;
  }

  while ( *longitude > 180.0 ) {
    *longitude -= 360.0;
  }

  POST02( is_valid_longitude( *longitude ), is_valid_latitude( *latitude ) );
}



void lambert_center( double* the_central_longitude,
                     double* the_central_latitude ) {
  
  PRE03( the_central_longitude, the_central_latitude, initialized );
  *the_central_longitude = central_longitude;
  *the_central_latitude  = central_latitude;
  POST02( is_valid_longitude_latitude( *the_central_longitude,
                                       *the_central_latitude ),
          IN_RANGE( *the_central_latitude, -89.0, 89.0 ) );
}



void lambert_tangents( double* the_lower_latitude,
                       double* the_upper_latitude ) {
  
  PRE03( the_lower_latitude, the_upper_latitude, initialized );
  *the_lower_latitude = lower_latitude;
  *the_upper_latitude = upper_latitude;
  POST06( is_valid_latitude( *the_lower_latitude ),
          is_valid_latitude( *the_upper_latitude ),
          *the_lower_latitude <= *the_upper_latitude,
          SIGN( *the_lower_latitude ) == SIGN( *the_upper_latitude ),
          IN_RANGE( *the_lower_latitude,
                    SIGN( *the_lower_latitude ) * 1.0,
                    SIGN( *the_lower_latitude ) * 89.0 ),
          IN_RANGE( *the_upper_latitude,
                    SIGN( *the_upper_latitude ) * 1.0,
                    SIGN( *the_upper_latitude ) * 89.0 ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



static void recompute_derived_terms( void ) {

  const double eccentricity0 = major_semiaxis == minor_semiaxis ? 0.0 :
    safe_quotient( sqrt( safe_difference( SQUARE( major_semiaxis ),
                                          SQUARE( minor_semiaxis ) ) ),
                   major_semiaxis );

  const double eccentricity1 = eccentricity0 > 1.0 ? 1.0 : eccentricity0;
  const double eccentricity_squared = SQUARE( eccentricity1 );
  const double phi0 = to_radians( central_latitude );
  const double phi1 = to_radians( lower_latitude );
  const double phi2 = to_radians( upper_latitude );
  const double sine_phi1   = sin( phi1 );
  const double cosine_phi1 = cos( phi1 );
  const double sine_phi2   = sin( phi2 );
  const double cosine_phi2 = cos( phi2 );
  /* Are lower/upper_latitude about equal? */
  const int is_tangent = phi1 + TOLERANCE >= phi2;

  eccentricity = eccentricity1;
  lambda0 = to_radians( central_longitude );
  n = sine_phi1;

  if ( eccentricity_squared != 0.0 ) { /* Ellipsoid planet: */
    const double m1  = msfn( sine_phi1, cosine_phi1, eccentricity_squared );
    const double ml1 = tsfn( phi1, sine_phi1, eccentricity );

    if ( ! is_tangent ) { /* Secant form: */
      const double numerator =
        log( m1 / msfn( sine_phi2, cosine_phi2, eccentricity_squared ) );
      const double denominator =
        log( ml1 / tsfn( phi2, sine_phi2, eccentricity ) );
      CHECK( denominator != 0.0 );
      n = numerator / denominator;
    }

    CHECK( n != 0.0 );
    c = m1 * pow( ml1, -n ) / n;

    if ( fabs( fabs( phi0 ) - PI_OVER_2 ) < TOLERANCE ) {
      rho0 = 0.0;
    } else {
      rho0 = c * pow( tsfn( phi0, sin( phi0 ), eccentricity ), n );
    }
  } else { /* Sphere planet: */
    const double denominator = tan( PI_OVER_4 + 0.5 * phi1 );

    if ( ! is_tangent ) { /* Secant form: */
      CHECK( ! about_equal( fabs( phi1 ), PI_OVER_2 ) ); /* Not near pole. */
      CHECK( ! about_equal( fabs( phi2 ), PI_OVER_2 ) ); /* Not near pole. */
      CHECK( cosine_phi1 != 0.0 );
      CHECK( cosine_phi2 != 0.0 );
      CHECK( tan( PI_OVER_4 + 0.5 * phi2 ) != 0.0 );
      CHECK( denominator != 0.0 );
      n = log( cosine_phi1 / cosine_phi2 ) /
          log( tan( PI_OVER_4 + 0.5 * phi2 ) / denominator );
    }

    c = cosine_phi1 * pow( denominator, n ) / n;

    if ( fabs( fabs( phi0 ) - PI_OVER_2 ) < TOLERANCE ) {
      rho0 = 0.0;
    } else {
      rho0 = c * pow( tan( PI_OVER_4 + 0.5 * phi0 ), -n );
    }
  }

  initialized = 1;

  POST07( initialized,
          ! is_nan( eccentricity ), IN_RANGE( eccentricity, 0.0, 1.0 ),
          ! is_nan( lambda0 ), ! is_nan( rho0 ), ! is_nan( n ), ! is_nan( c));
}



static double phi2_iterate( double ts, double the_eccentricity ) {

  PRE0( IN_RANGE( the_eccentricity, 0.0, 1.0 ) );
  const int maximum_iterations       = MAXIMUM_ITERATIONS;
  const double convergence_tolerance = CONVERGENCE_TOLERANCE;
  const double half_eccentricity     = the_eccentricity * 0.5;
  double delta_phi = 0.0;
  int iteration = 0;
  double result = PI_OVER_2 - 2.0 * atan( ts );

  do {
    const double con = the_eccentricity * sin( result );
    CHECK( con != -1.0 );
    delta_phi =
      PI_OVER_2 -
      2.0 * atan( ts * pow( ( 1.0 - con ) / ( 1.0 + con ), half_eccentricity))
      - result;
    result += delta_phi;
    ++iteration;
  } while ( fabs( delta_phi ) >= convergence_tolerance &&
          iteration < maximum_iterations );

  POST0( ! is_nan( result ) );
  return result;
}




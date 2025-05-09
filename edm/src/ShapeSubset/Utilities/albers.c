
/*
albers.c - Albers Conformal Conic projections derived from USGS PROJ.
http://mathworld.wolfram.com/AlbersEqual-AreaConicProjection.html
*/

#include <math.h>  /* For M_PI, sqrt(), sin(), hypot(), atan2(), log(), pow()*/
#include <float.h> /* For DBL_MAX. */

#include <Assertions.h>  /* For PRE*(), POST*(), CHECK*(), IN_RANGE(). */
#include <projections.h> /* For is_valid_latitude(), etc. */
#include <albers.h>      /* For public interface. */

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
static double one_minus_eccentricity_squared;
static double lambda0;           /* Central longitude in radians.           */
static double rho0;              /* See USGS PROJ Library.                  */
static double n;                 /* See USGS PROJ Library.                  */
static double n2;                /* See USGS PROJ Library.                  */
static double c;                 /* See USGS PROJ Library.                  */
static double ec;                /* See USGS PROJ Library.                  */
static double dd;                /* See USGS PROJ Library.                  */
static int initialized;          /* Has client called initialize_albers()?  */

/*=========================== FORWARD DECLARATIONS ==========================*/

static void recompute_derived_terms( void );

static double phi1_iterate( double phi, double the_eccentricity,
                            double the_one_minus_eccentricity_squared );

/*============================= PUBLIC FUNCTIONS ============================*/


void initialize_albers(
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



void project_albers( double longitude, double latitude,
                      double* x, double* y ) {

  PRE05( is_valid_longitude( longitude ), is_valid_latitude( latitude ),
         x, y, initialized );

  double lambda = to_radians( longitude );
  double phi    = to_radians( latitude  );
  double rho = 0.0;            /* See USGS PROJ Library. */
  double lambda_delta = 0.0;   /* Radians from central latitude. */
  double n_lambda_delta = 0.0; /* n times lambda_delta. */

  /*
   * If phi is too near a pole tweak it so that projecting
   * succeeds and unprojecting yields original longitude
   * (instead of central longitude).
   */

  if ( ! IN_RANGE( phi, -PI_OVER_2 + TOLERANCE, PI_OVER_2 - TOLERANCE ) ) {
    phi = phi + TOLERANCE * -SIGN( phi );
  }

  rho = c - n * qsfn( sin(phi), eccentricity, one_minus_eccentricity_squared );
  CHECK( rho >= 0.0 );
  rho = sqrt( rho ) * dd;

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



void unproject_albers( double x, double y,
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

    CHECK4( c != 0.0, n != 0.0, rho != 0.0, dd != 0.0 );
    phi = rho / dd;

    if ( eccentricity != 0.0 ) { /* Ellipsoid: */
      phi = ( c - phi * phi ) / n;

      if ( fabs( ec - fabs( phi ) ) > TOLERANCE ) {
        phi = phi1_iterate( phi, eccentricity, one_minus_eccentricity_squared);
      } else {
        phi = phi < 0.0 ? -PI_OVER_2 : PI_OVER_2;
      }
    } else { /* Sphere: */
      phi = ( c - SQUARE( phi ) ) / n2;

      if ( fabs( phi ) < 1.0 ) {
        phi = asin( phi );
      } else {
        phi = phi < 0.0 ? -PI_OVER_2 : PI_OVER_2;
      }
    }

    lambda = atan2( xp, yp_delta ) / n;
  } else {
    phi = n > 0.0 ? PI_OVER_2 : -PI_OVER_2;
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



void albers_center( double* the_central_longitude,
                    double* the_central_latitude ) {

  PRE03( the_central_longitude, the_central_latitude, initialized );
  *the_central_longitude = central_longitude;
  *the_central_latitude  = central_latitude;
  POST02( is_valid_longitude_latitude( *the_central_longitude,
                                       *the_central_latitude ),
          IN_RANGE( *the_central_latitude, -89.0, 89.0 ) );
}



void albers_tangents( double* the_lower_latitude,
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
  const double sine_phi0   = sin( phi0 );
  const double sine_phi1   = sin( phi1 );
  const double cosine_phi1 = cos( phi1 );
  const double sine_phi2   = sin( phi2 );
  const double cosine_phi2 = cos( phi2 );
  /* Are lower/upper_latitude about equal? */
  const int is_tangent = phi1 + TOLERANCE >= phi2;

  eccentricity = eccentricity1;
  one_minus_eccentricity_squared = 1.0 - eccentricity_squared;
  lambda0 = to_radians( central_longitude );
  n = sine_phi1;

  if ( eccentricity_squared != 0.0 ) { /* Ellipsoid planet: */
    const double m1  = msfn( sine_phi1, cosine_phi1, eccentricity_squared );
    const double ml1 =
      qsfn( sine_phi1, eccentricity, one_minus_eccentricity_squared );

    if ( ! is_tangent ) { /* Secant form: */
      const double m2  = msfn( sine_phi2, cosine_phi2, eccentricity_squared );
      const double ml2 =
        qsfn( sine_phi2, eccentricity, one_minus_eccentricity_squared );
      CHECK( ml1 != ml2 );
      n = ( SQUARE( m1 ) - SQUARE( m2 ) ) / ( ml2 - ml1 );
    }

    CHECK2( n != 0.0, eccentricity != 0.0 );
    ec =
      1.0 - 0.5 * one_minus_eccentricity_squared *
      log( ( 1.0 - eccentricity ) / ( 1.0 + eccentricity ) ) / eccentricity;
    c = SQUARE( m1 ) + n * ml1;
    dd = 1.0 / n;
    rho0 =
      dd * sqrt( c - n *
                 qsfn( sine_phi0, eccentricity, one_minus_eccentricity_squared ) );
  } else { /* Sphere planet: */

    if ( ! is_tangent ) { /* Secant form: */
      n = 0.5 * ( n + sine_phi2 );
    }

    CHECK( ! about_equal( fabs( phi1 ), PI_OVER_2 ) ); /* Not near pole. */
    CHECK( ! about_equal( fabs( phi2 ), PI_OVER_2 ) ); /* Not near pole. */
    CHECK( cosine_phi1 != 0.0 );
    CHECK( cosine_phi2 != 0.0 );
    n2 = n + n;
    c = SQUARE( cosine_phi1 ) + n2 * sine_phi1;
    CHECK2( n != 0.0, c > n2 * sine_phi0 );
    dd = 1.0 / n;
    rho0 = dd * sqrt( c - n2 * sine_phi0 );
  }

  initialized = 1;

  POST012( initialized,
           ! is_nan( eccentricity ),
           IN_RANGE( eccentricity, 0.0, 1.0 ),
           IN_RANGE( one_minus_eccentricity_squared, 0.0, 1.0 ),
           about_equal( one_minus_eccentricity_squared,
                        1.0 - SQUARE( eccentricity ) ),
           ! is_nan( lambda0 ),
           ! is_nan( rho0 ),
           ! is_nan( n ),
           ! is_nan( n2 ),
           ! is_nan( c ),
           ! is_nan( dd ),
           ! is_nan( ec ) );
}



static double phi1_iterate( double phi, double the_eccentricity,
                            double the_one_minus_eccentricity_squared ) {

  PRE02( IN_RANGE( the_eccentricity, 0.0, 1.0 ),
         IN_RANGE( the_one_minus_eccentricity_squared, 0.0, 1.0 ) );
  double result = asin( 0.5 * phi );

  if ( the_eccentricity > PROJECTION_TOLERANCE ) {
    const int maximum_iterations       = MAXIMUM_ITERATIONS;
    const double convergence_tolerance = CONVERGENCE_TOLERANCE;
    double delta_phi = 0.0;
    int iteration = 0;

    do {
      const double sine_phi = sin( result );
      const double cosine_phi = cos( result );
      const double con = the_eccentricity * sine_phi;
      const double com = 1.0 - SQUARE( con );
      CHECK4( cosine_phi != 0.0, con != -1.0, com != 0.0,
              the_one_minus_eccentricity_squared != 0.0 );
      delta_phi =
        0.5 * SQUARE( com ) / cosine_phi *
        ( phi / the_one_minus_eccentricity_squared -
          sine_phi / com +
          0.5 / the_eccentricity * log( ( 1.0 - con ) / ( 1.0 + con ) ) );
      result += delta_phi;
      ++iteration;
    } while ( fabs( delta_phi ) >= convergence_tolerance &&
              iteration < maximum_iterations );
  }

  POST0( ! is_nan( result ) );
  return result;
}




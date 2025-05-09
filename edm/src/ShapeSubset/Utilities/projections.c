
/* projections.c - Common code for projecting/unprojecting. */

#include <float.h> /* For DBL_MIN, DBLMAX. */
#include <math.h>  /* For M_PI, fabs().    */

#include <projections.h> /* For public interface. */

#include <Assertions.h> /* For PRE0*(), POST0*(), IN_RANGE(), IS_BOOL(). */


static int within_tolerance( double x, double y, double tolerance );


int is_valid_ellipsoid( double major_semiaxis, double minor_semiaxis ) {
  const int result =
    AND7( ! is_nan( major_semiaxis ),
          ! is_nan( minor_semiaxis ),
          major_semiaxis > 0.0,
          minor_semiaxis > 0.0,
          major_semiaxis >= minor_semiaxis,
          SQUARE( major_semiaxis ) > 0.0,
          SQUARE( minor_semiaxis ) > 0.0 );
  POST0( IS_BOOL( result ) );
  return result;
}


int is_valid_longitude( double longitude ) {
  const int result =
    AND2( ! is_nan( longitude ), IN_RANGE( longitude, -180.0, 180.0 ) );
  POST0( IS_BOOL( result ) );
  return result;
}


int is_valid_latitude( double latitude ) {
  const int result =
    AND2( ! is_nan( latitude ), IN_RANGE( latitude, -90.0, 90.0 ) );
  POST0( IS_BOOL( result ) );
  return result;
}


int is_valid_longitude_latitude( double longitude, double latitude ) {
  return is_valid_longitude( longitude ) && is_valid_latitude( latitude );
}


int is_valid_longitudes_and_latitudes( long long count,
                                       const double longitudes[],
                                       const double latitudes[] ) {
  PRE03( count > 0, longitudes, latitudes );
  int result = 0;
  long long index = 0;

  do {

    if ( ! AND2( is_valid_longitude( longitudes[ index ] ),
                 is_valid_latitude(  latitudes[ index ] ) ) ) {
      index = count;
    }

    ++index;
  } while ( index < count );

  result = index == count;
  POST0( IS_BOOL( result ) );
  return result;
}


int is_nan( double value ) {
  const double copy = value;
  const int result = (copy != value);
  POST0( IS_BOOL( result ) );
  return result;
}


double to_radians( double the_degrees ) {
  PRE0( ! is_nan( the_degrees ) );
  const double result = the_degrees * ( M_PI / 180.0 );
  POST03( ! is_nan( result ),
          OR2( SIGN( result ) == SIGN( the_degrees ), result == 0.0 ),
          fabs( result ) <= fabs( the_degrees ) );
  return result;
}


double to_degrees( double the_radians ) {
  PRE0( ! is_nan( the_radians ) );
  const double result = the_radians * ( 180.0 / M_PI );
  POST03( ! is_nan( result ),
          OR2( SIGN( result ) == SIGN( the_radians ), result == 0.0 ),
          fabs( result ) >= fabs( the_radians ) );
  return result;

}


double safe_difference( double left, double right ) {
  PRE02( ! is_nan( left ), ! is_nan( right ) );
  const double result = left == right ? 0.0 : left - right;
  POST02( ! is_nan( result ), IMPLIES( left == right, result == 0.0 ) );
  return result;
}


double safe_quotient( double numerator, double denominator ) {
  PRE03( ! is_nan( numerator ), ! is_nan( denominator ), denominator != 0.0 );
  const double result = numerator   ==  0.0 ? 0.0
                      : denominator ==  1.0 ? numerator
                      : denominator == -1.0 ? -numerator
                      : numerator == denominator ? 1.0
                      : numerator == -denominator ? -1.0
                      : numerator / denominator;
  POST04( ! is_nan( result ),
          IMPLIES( numerator ==  0.0,         result ==  0.0 ),
          IMPLIES( numerator ==  denominator, result ==  1.0 ),
          IMPLIES( numerator == -denominator, result == -1.0 ) );
  return result;
}


int about_equal( double x, double y ) {
  return within_tolerance( x, y, TOLERANCE );
}


/*
Derived from Squassabia, Alberto, "Comparing Floats", C++ Report, Vol 12,
No 2, February 2000, pp 30-32, 39. SIGS Publications.
*/

static int within_tolerance( double x, double y, double tolerance ) {
  PRE02( ! is_nan( tolerance ), tolerance <= 0.1 );

  /* First try bitwise comparison (handles nans): */

  const double* const xp = &x;
  const double* const yp = &y;
  const long long* const ixp = (const long long*) xp;
  const long long* const iyp = (const long long*) yp;
  const long long ix = *ixp;
  const long long iy = *iyp;
  int result = ix == iy;

  if ( result == 0 ) {
    if ( x == 0.0 ) {
      result = IN_RANGE( y, -tolerance, tolerance ); /* Close enough to 0? */
    } else if ( y == 0.0 ) {
      result = IN_RANGE( x, -tolerance, tolerance ); /* Close enough to 0? */
    } else if ( IN_RANGE( x, y - tolerance, y + tolerance)) {/*Or each other?*/
      result = 1;
    } else if ( IN_RANGE( y, x - tolerance, x + tolerance)) {/*Or each other?*/
      result = 1;
    } else { /* Ratio handles cases of large values differing in last digits.*/
      const double ax = fabs( x );
      const double ay = fabs( y );

      if ( AND2( ay < 1.0, ax > ay * DBL_MAX ) ) { /* Avoid overflow. */
        result = 0;
      } else if ( AND2( ay > 1.0, ax < ay * DBL_MIN ) ) { /* Avoid underflow.*/
        result = 0;
      } else {
        const double ratio = x / y;
        result = IN_RANGE( ratio, 1.0 - tolerance, 1.0 + tolerance );
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: ssfn - See USGS PROJ Library.
INPUTS:  double phi      Angle in radians.
         double sine_phi  Sine of phi.
         double ellipsoid_eccentricity  Of planet approximation.
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double ssfn( double phi, double sine_phi, double ellipsoid_eccentricity ) {

  PRE07( ! is_nan( phi ),
         ! is_nan( sine_phi ),
         ! is_nan( ellipsoid_eccentricity ),
         within_tolerance( sine_phi, sin( phi ), PROJECTION_TOLERANCE ),
         sine_phi > -1.0,
         sine_phi < 1.0,
         IN_RANGE( ellipsoid_eccentricity, 0.0, 1.0 ) );

  const double eccentricity_sine_phi = ellipsoid_eccentricity * sine_phi;
  const double exponent = ellipsoid_eccentricity * 0.5;
  const double factor1 = tan( ( PI_OVER_2 + phi ) * 0.5 );
  const double factor2 = pow( ( ( 1.0 - eccentricity_sine_phi ) /
                                ( 1.0 + eccentricity_sine_phi ) ),
                              exponent );
  const double result = factor1 * factor2;

  POST0( ! is_nan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: msfn - See USGS PROJ Library.
INPUTS:  double sine_phi    Sine of phi.
         double cosine_phi  Cosine of phi.
         double eccentricity_squared  Of planet approximation.
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double msfn(double sine_phi, double cosine_phi, double eccentricity_squared) {

  PRE012( ! is_nan( sine_phi ),
          ! is_nan( cosine_phi ),
          ! is_nan( eccentricity_squared ),
          within_tolerance( sine_phi, sqrt( 1.0 - SQUARE( cosine_phi ) ),
                           PROJECTION_TOLERANCE ),
          sine_phi   > -1.0,
          sine_phi   < 1.0,
          cosine_phi > -1.0,
          cosine_phi < 1.0,
          cosine_phi != 0.0,
          IN_RANGE( eccentricity_squared, 0.0, 1.0 ),
          eccentricity_squared * sine_phi * sine_phi < 1.0,
          sqrt( 1.0 - eccentricity_squared * SQUARE( sine_phi ) ) != 0.0 );

  const double result =
    cosine_phi / sqrt( 1.0 - eccentricity_squared * SQUARE( sine_phi ) );

  POST02( ! is_nan( result ), result != 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: tsfn - See USGS PROJ Library.
INPUTS:  double phi       Angle in radians.
         double sine_phi  Sine of phi.
         double ellipsoid_eccentricity  Of planet approximation.
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double tsfn( double phi, double sine_phi, double ellipsoid_eccentricity ) {

  PRE010( ! is_nan( phi ),
          ! is_nan( sine_phi ),
          ! is_nan( ellipsoid_eccentricity ),
          within_tolerance( sine_phi, sin( phi ), PROJECTION_TOLERANCE ),
          sine_phi > -1.0,
          sine_phi < 1.0,
          IN_RANGE( ellipsoid_eccentricity, 0.0, 1.0 ),
          tan( ( PI_OVER_2 - phi ) * 0.5 ) != 0.0,
          fabs( ellipsoid_eccentricity * sine_phi ) != 1.0,
          ( 1.0 + ellipsoid_eccentricity * sine_phi ) != 0.0 );

  const double eccentricity_sine_phi = ellipsoid_eccentricity * sine_phi;
  const double exponent = ellipsoid_eccentricity * 0.5;
  const double numerator = tan( ( PI_OVER_2 - phi ) * 0.5 );
  const double denominator = pow( ( ( 1.0 - eccentricity_sine_phi ) /
                                    ( 1.0 + eccentricity_sine_phi ) ),
                                  exponent );
  const double result = numerator / denominator;

  POST02( ! is_nan( result ), result != 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: qsfn - See USGS PROJ Library.
INPUTS:  double sine_phi                         Sine of phi.
         double ellipsoid_eccentricity           Of planet approximation.
         double one_minus_ellipsoid_eccentricity 1 - ellipsoid_eccentricity^2.
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double qsfn( double sine_phi, double ellipsoid_eccentricity,
             double one_minus_ellipsoid_eccentricity_squared ) {

  PRE08( ! is_nan( sine_phi ),
         ! is_nan( ellipsoid_eccentricity ),
         ! is_nan( one_minus_ellipsoid_eccentricity_squared ),
         sine_phi > -1.0,
         sine_phi < 1.0,
         IN_RANGE( ellipsoid_eccentricity, 0.0, 1.0 ),
         IN_RANGE( one_minus_ellipsoid_eccentricity_squared, 0.0, 1.0 ),
         within_tolerance( one_minus_ellipsoid_eccentricity_squared,
                           1.0 - SQUARE( ellipsoid_eccentricity ),
                           PROJECTION_TOLERANCE ) );

  double result = 0.0;

  if ( ellipsoid_eccentricity < PROJECTION_TOLERANCE ) {
    result = sine_phi + sine_phi;
  } else {
    const double con = ellipsoid_eccentricity * sine_phi;
    CHECK3( con != 1.0, con != -1.0, ellipsoid_eccentricity != 0.0 );
    result =
      one_minus_ellipsoid_eccentricity_squared *
        ( sine_phi / ( 1.0 - SQUARE( con ) ) -
          0.5 / ellipsoid_eccentricity *
          log( ( 1.0 - con ) / ( 1.0 + con ) ) );
  }

  POST0( ! is_nan( result ) );
  return result;
}




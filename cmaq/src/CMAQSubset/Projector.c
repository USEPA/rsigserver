/******************************************************************************
PURPOSE: Projector.c - Define Cartographic projectors ADT ABC.
NOTES:   See Lambert.h for example usage.
HISTORY: 2004-10-01 plessel.todd@epa.gov Created based on C++ version.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */
#include <float.h>  /* For DBL_MIN, DBLMAX. */
#include <math.h>   /* For M_PI, fabs(), sqrt(), sin(), tan(), atan(), pow().*/

#include <Assertions.h>    /* For PRE0*(), POST0*(), IN_RANGE(), IS_BOOL().  */
#include <Projector.h>     /* For public interface.                          */

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: isNan - Is x a NaN (Not a Number)?
INPUTS:  double x - The 64-bit value to test.
RETURNS: int 1 if x is a NaN, else 0.
******************************************************************************/

int isNan( double x ) {
  const double copy = x;
  const int result = (copy != x);
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeDifference - NaN-free difference.
INPUTS:  double left  - The left value of difference.
         double right - The right value of difference.
RETURNS: double difference of arguments.
******************************************************************************/

double safeDifference( double x, double y ) {
  PRE02( ! isNan( x ), ! isNan( y ) );
  const double result = x == y ? 0.0 : x - y;
  POST02( ! isNan( result ), IMPLIES( x == y, result == 0.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: safeQuotient - NaN-free quotient.
INPUTS:  double numerator    The numerator   of ratio.
         double denominator  The denominator of ratio.
RETURNS: double quotient of arguments.
******************************************************************************/

double safeQuotient( double numerator, double denominator ) {
  PRE03( ! isNan( numerator ), ! isNan( denominator ), denominator != 0.0 );
  const double result =   numerator   ==  0.0 ? 0.0
                        : denominator ==  1.0 ? numerator
                        : denominator == -1.0 ? -numerator
                        : numerator == denominator ? 1.0
                        : numerator == -denominator ? -1.0
                        : numerator / denominator;
  POST04( ! isNan( result ),
          IMPLIES( numerator ==  0.0,         result ==  0.0 ),
          IMPLIES( numerator ==  denominator, result ==  1.0 ),
          IMPLIES( numerator == -denominator, result == -1.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: withinTolerance - Do x and y differ by less than (non-negative,
         finite) tolerance (i.e., fabs( x - y ) <= tolerance) or, for large
         values, differ only in digits beyond number of significant digits in
         tolerance? (E.g., tolerance = 1e-6 and 1.0000000001e30 and
         1.0000000002e30 are considered equal.)
INPUTS:  double x - first value to compare.
         double y - second value to compare.
         double tolerance - tolerance threshold (e.g., 1-e6).
RETURNS: int 1 if x and y differ (in significant digits) by less than
         tolerance, else 0.
NOTES:   This function is commutative:
           withinTolerance( x, y, tolerance ) ==
           withinTolerance( y, x, tolerance )
         but not transitive:
           (withinTolerance( x, y, tolerance ) &&
            withinTolerance( x, z, tolerance ))
         does not imply
            withinTolerance( y, z, tolerance )
         (due to pivoting, e.g., about zero:
         if x == 0.0 and y = -tolerance and z == +tolerance,
         x ~ y and x ~ z but y !~ z)
         See: Squassabia, Alberto, "Comparing Floats", C++ Report, Vol 12,
         No 2, February 2000, pp 30-32, 39. SIGS Publications.
******************************************************************************/

int withinTolerance( double x, double y, double tolerance ) {

  PRE02( ! isNan( tolerance ), tolerance <= 0.1 );
  assert_static( sizeof (double) == sizeof (long long) );

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
PURPOSE: aboutEqual - Is withinTolerance( x, y, TOLERANCE )?
INPUTS:  double x - first value to compare.
         double y - second value to compare.
RETURNS: int 1 if x and y differ (in significant digits) by less than
         TOLERANCE, else 0.
******************************************************************************/

int aboutEqual( double x, double y ) {
  const int result = withinTolerance( x, y, TOLERANCE );
  POST0( result == withinTolerance( x, y, TOLERANCE ) );
  return result;
}



/******************************************************************************
PURPOSE: radians - Radians of degrees.
INPUTS:  double theDegrees - The degrees to convert.
RETURNS: double radians of degrees.
******************************************************************************/

double radians( double theDegrees ) {
  PRE0( ! isNan( theDegrees ) );
  const double result = theDegrees * ( M_PI / 180.0 );
  POST03( ! isNan( result ),
          OR2( SIGN( result ) == SIGN( theDegrees ), result == 0.0 ),
          fabs( result ) <= fabs( theDegrees ) );
  return result;
}



/******************************************************************************
PURPOSE: degrees - Degrees of radians.
INPUTS:  double theDegrees - The degrees to convert.
RETURNS: double radians of degrees.
******************************************************************************/

double degrees( double theRadians ) {
  PRE0( ! isNan( theRadians ) );
  const double result = theRadians * ( 180.0 / M_PI );
  POST03( ! isNan( result ),
          OR2( SIGN( result ) == SIGN( theRadians ), result == 0.0 ),
          fabs( result ) >= fabs( theRadians ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidEllipsoid - Do the arguments define a valid ellipsoid?
INPUTS:  double majorSemiaxis Mean equitorial radius of planet approximation.
         double minorSemiaxis Mean polar      radius of planet approximation.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidEllipsoid( double majorSemiaxis, double minorSemiaxis ) {
  const int result =
    AND7( ! isNan( majorSemiaxis ),
          ! isNan( minorSemiaxis ),
          majorSemiaxis > 0.0,
          minorSemiaxis > 0.0,
          majorSemiaxis >= minorSemiaxis,
          SQUARE( majorSemiaxis ) > 0.0,
          SQUARE( minorSemiaxis ) > 0.0 );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidLongitude - Is the argument a valid longitude?
INPUTS:  double longitude In degrees.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidLongitude( double longitude ) {
  const int result = IN_RANGE( longitude, -180.0, 180.0 );
  POST0( result == IN_RANGE( longitude, -180.0, 180.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidLatitude - Is the argument a valid latitude?
INPUTS:  double latitude In degrees.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidLatitude( double latitude ) {
  const int result = IN_RANGE( latitude, -90.0, 90.0 );
  POST0( result == IN_RANGE( latitude, -90.0, 90.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidLongitudeLatitude - Are the arguments a valid longitude
         latitude point?
INPUTS:  double longitude In degrees.
         double latitude  In degrees.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidLongitudeLatitude( double longitude, double latitude ) {
  const int result =
    AND2( isValidLongitude( longitude ), isValidLatitude( latitude ) );
  POST0( result ==
         AND2( isValidLongitude( longitude ), isValidLatitude( latitude ) ) );
  return result;
}



/******************************************************************************
PURPOSE: validLongitudesAndLatitudes - Are longitudes and latitudes valid?
INPUTS:  size_t count            Number of points.
         const double longitudes[]  Longitudes to check.
         const double latitudes[]   Latitudes to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int validLongitudesAndLatitudes( size_t count,
                                 const double longitudes[],
                                 const double latitudes[] ) {
  PRE03( count > 0, longitudes, latitudes );
  int result = 0;
  int index = 0;

  do {

    if ( ! isValidLongitudeLatitude( longitudes[ index ], latitudes[ index])) {
      index = count;
    }

    ++index;
  } while ( index < count );

  result = index == count;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: ssfn - See USGS PROJ Library.
INPUTS:  double phi      Angle in radians.
         double sinePhi  Sine of phi.
         double ellipsoidEccentricity  Of planet approximation.
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double ssfn( double phi, double sinePhi, double ellipsoidEccentricity ) {

  PRE04( withinTolerance( sinePhi, sin( phi ), PROJECTION_TOLERANCE ),
         sinePhi > -1.0,
         sinePhi < 1.0,
         IN_RANGE( ellipsoidEccentricity, 0.0, 1.0 ) );

  const double eccentricitySinePhi = ellipsoidEccentricity * sinePhi;
  const double exponent = ellipsoidEccentricity * 0.5;
  const double factor1 = tan( ( PI_OVER_2 + phi ) * 0.5 );
  const double factor2 = pow( ( ( 1.0 - eccentricitySinePhi ) /
                                ( 1.0 + eccentricitySinePhi ) ),
                              exponent );
  const double result = factor1 * factor2;

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: msfn - See USGS PROJ Library.
INPUTS:  double sinePhi    Sine of phi.
         double cosinePhi  Cosine of phi.
         double eccentricitySquared  Of planet approximation.
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double msfn( double sinePhi, double cosinePhi, double eccentricitySquared ) {

  PRE09( withinTolerance( sinePhi, sqrt( 1.0 - SQUARE( cosinePhi ) ),
                          PROJECTION_TOLERANCE ),
         sinePhi   > -1.0,
         sinePhi   < 1.0,
         cosinePhi > -1.0,
         cosinePhi < 1.0,
         cosinePhi != 0.0,
         IN_RANGE( eccentricitySquared, 0.0, 1.0 ),
         eccentricitySquared * sinePhi * sinePhi < 1.0,
         sqrt( 1.0 - eccentricitySquared * SQUARE( sinePhi ) ) != 0.0 );

  const double result =
    cosinePhi / sqrt( 1.0 - eccentricitySquared * SQUARE( sinePhi ) );

  POST02( ! isNan( result ), result != 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: qsfn - See USGS PROJ Library.
INPUTS:  double sinePhi                         Sine of phi.
         double ellipsoidEccentricity           Of planet approximation.
         double oneMinusEllipsoidEccentricitySquared 1-ellipsoidEccentricity^2
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double qsfn( double sinePhi, double ellipsoidEccentricity,
             double oneMinusEllipsoidEccentricitySquared ) {

  PRE08( ! isNan( sinePhi ),
         ! isNan( ellipsoidEccentricity ),
         ! isNan( oneMinusEllipsoidEccentricitySquared ),
         sinePhi > -1.0,
         sinePhi < 1.0,
         IN_RANGE( ellipsoidEccentricity, 0.0, 1.0 ),
         IN_RANGE( oneMinusEllipsoidEccentricitySquared, 0.0, 1.0 ),
         withinTolerance( oneMinusEllipsoidEccentricitySquared,
                          1.0 - SQUARE( ellipsoidEccentricity ),
                          PROJECTION_TOLERANCE ) );

  double result = 0.0;

  if ( ellipsoidEccentricity < PROJECTION_TOLERANCE ) {
    result = sinePhi + sinePhi;
  } else {
    const double con = ellipsoidEccentricity * sinePhi;
    CHECK3( con != 1.0, con != -1.0, ellipsoidEccentricity != 0.0 );
    result =
      oneMinusEllipsoidEccentricitySquared *
        ( sinePhi / ( 1.0 - SQUARE( con ) ) -
          0.5 / ellipsoidEccentricity *
          log( ( 1.0 - con ) / ( 1.0 + con ) ) );
  }

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: tsfn - See USGS PROJ Library.
INPUTS:  double phi      Angle in radians.
         double sinePhi  Sine of phi.
         double ellipsoidEccentricity  Of planet approximation.
RETURNS: double See USGS PROJ Library.
******************************************************************************/

double tsfn( double phi, double sinePhi, double ellipsoidEccentricity ) {

  PRE07( withinTolerance( sinePhi, sin( phi ), PROJECTION_TOLERANCE ),
         sinePhi > -1.0,
         sinePhi < 1.0,
         IN_RANGE( ellipsoidEccentricity, 0.0, 1.0 ),
         tan( ( PI_OVER_2 - phi ) * 0.5 ) != 0.0,
         fabs( ellipsoidEccentricity * sinePhi ) != 1.0,
         ( 1.0 + ellipsoidEccentricity * sinePhi ) != 0.0 );

  const double eccentricitySinePhi = ellipsoidEccentricity * sinePhi;
  const double exponent = ellipsoidEccentricity * 0.5;
  const double numerator = tan( ( PI_OVER_2 - phi ) * 0.5 );
  const double denominator = pow( ( ( 1.0 - eccentricitySinePhi ) /
                                    ( 1.0 + eccentricitySinePhi ) ),
                                  exponent );
  const double result = numerator / denominator;

  POST02( ! isNan( result ), result != 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: phi1Iterate - Iterate on unprojected y coordinate.
RETURNS: double converged phi.
******************************************************************************/

double phi1Iterate( double phi,
                    double eccentricity,
                    double oneMinusEccentricitySquared ) {

  PRE02( IN_RANGE( eccentricity, 0.0, 1.0 ),
         IN_RANGE( oneMinusEccentricitySquared, 0.0, 1.0 ) );

  double result = asin( 0.5 * phi );

  if ( eccentricity > PROJECTION_TOLERANCE ) {
    const int maximumIterations       = MAXIMUM_ITERATIONS;
    const double convergenceTolerance = CONVERGENCE_TOLERANCE;
    double deltaPhi = 0.0;
    int iteration = 0;

    do {
      const double sinePhi = sin( result );
      const double cosinePhi = cos( result );
      const double con = eccentricity * sinePhi;
      const double com = 1.0 - SQUARE( con );
      CHECK4( cosinePhi != 0.0, con != -1.0, com != 0.0,
              oneMinusEccentricitySquared != 0.0 );
      deltaPhi =
        0.5 * SQUARE( com ) / cosinePhi *
        ( phi / oneMinusEccentricitySquared -
          sinePhi / com +
          0.5 / eccentricity * log( ( 1.0 - con ) / ( 1.0 + con ) ) );
      result += deltaPhi;
      ++iteration;
    } while ( fabs( deltaPhi ) >= convergenceTolerance &&
              iteration < maximumIterations );
  }

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: phi2Iterate - Iterate on unprojected y coordinate.
RETURNS: double converged phi.
******************************************************************************/

double phi2Iterate( double ts, double theEccentricity ) {

  PRE0( IN_RANGE( theEccentricity, 0.0, 1.0 ) );

  const int maximumIterations = MAXIMUM_ITERATIONS;
  const double convergenceTolerance = CONVERGENCE_TOLERANCE;
  const double halfEccentricity     = theEccentricity * 0.5;
  double deltaPhi = 0.0;
  int iteration = 0;
  double result = PI_OVER_2 - 2.0 * atan( ts );

  do {
    const double con = theEccentricity * sin( result );
    CHECK( con != -1.0 );
    deltaPhi =
      PI_OVER_2 -
      2.0 * atan( ts * pow( ( 1.0 - con ) / ( 1.0 + con ), halfEccentricity))
      - result;
    result += deltaPhi;
    ++iteration;
  } while ( AND2( fabs( deltaPhi ) >= convergenceTolerance,
                  iteration < maximumIterations ) );

  POST0( ! isNan( result ) );
  return result;
}





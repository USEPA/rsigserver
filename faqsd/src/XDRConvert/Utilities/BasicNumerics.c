/******************************************************************************
PURPOSE: BasicNumerics.c - Defines routines for manipulating 64-bit types.

NOTES:   See BasicNumerics.h.

HISTORY: 2001/04 Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <limits.h> /* For DBL_MAX, DBL_MIN.             */
#include <float.h>  /* For DBL_MAX, DBL_MIN on DEC OSF1. */
#include <math.h>   /* For fabs().                       */
#include <stdlib.h> /* For drand48() or random() to "implement" drand48. */
#include <ctype.h>  /* For isspace().                    */
#include <errno.h>  /* For errno.                        */

#include <Assertions.h>    /* For PRE0*(), POST0*().  */
#include <BasicNumerics.h> /* For public definitions. */

/*================================= MACROS ==================================*/

#ifdef MIN
#undef MIN
#endif
#define MIN(a,b) ((a)<(b)?(a):(b))

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: isSignedChar - Is value representable without loss as a signed char?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isSignedChar( Integer value ) {
  const Integer result = IN_RANGE( value, SCHAR_MIN, SCHAR_MAX );
  POST0( result == IN_RANGE( value, SCHAR_MIN, SCHAR_MAX ) );
  return result;
}



/******************************************************************************
PURPOSE: isUnsignedChar - Is value representable without loss as an unsigned
         char?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isUnsignedChar( Integer value ) {
  const Integer result = IN_RANGE( value, 0, UCHAR_MAX );
  POST0( result == IN_RANGE( value, 0, UCHAR_MAX ) );
  return result;
}



/******************************************************************************
PURPOSE: isSignedShort - Is value representable without loss as a signed
         short?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isSignedShort( Integer value ) {
  const Integer result = IN_RANGE( value, SHRT_MIN, SHRT_MAX );
  POST0( result == IN_RANGE( value, SHRT_MIN, SHRT_MAX ) );
  return result;
}



/******************************************************************************
PURPOSE: isUnsignedShort - Is value representable without loss as an unsigned
         short?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isUnsignedShort( Integer value ) {
  const Integer result = AND2( value >= 0,
                               OR2(sizeof (unsigned short) >= sizeof (Integer),
                                   value <= USHRT_MAX ) );
  POST0( result == AND2( value >= 0,
                         OR2( sizeof (unsigned short) >= sizeof (Integer),
                              value <= USHRT_MAX ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isSignedInt - Is value representable without loss as a signed int?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isSignedInt( Integer value ) {
  const Integer result = IN_RANGE( value, INT_MIN, INT_MAX );
  POST0( result == IN_RANGE( value, INT_MIN, INT_MAX ) );
  return result;
}



/******************************************************************************
PURPOSE: isUnsignedInt - Is value representable without loss as an unsigned
         int?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isUnsignedInt( Integer value ) {
  const Integer result = AND2( value >= 0,
                               OR2( sizeof (unsigned int) >= sizeof (Integer),
                                    value <= UINT_MAX ) );
  POST0( result == AND2( value >= 0,
                         OR2( sizeof (unsigned int) >= sizeof (Integer),
                              value <= UINT_MAX ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isSignedLong - Is value representable without loss as a signed long?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isSignedLong( Integer value ) {
  const Integer result = IN_RANGE( value, LONG_MIN, LONG_MAX );
  POST0( result == IN_RANGE( value, LONG_MIN, LONG_MAX ) );
  return result;
}



/******************************************************************************
PURPOSE: isUnsignedLong - Is value representable without loss as an unsigned
         long?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isUnsignedLong( Integer value ) {
  const Integer result = AND2( value >= 0,
                               OR2(sizeof (unsigned long) >= sizeof (Integer),
                                   value <= ULONG_MAX ) );
  POST0( result == AND2( value >= 0,
                         OR2( sizeof (unsigned long) >= sizeof (Integer),
                              value <= ULONG_MAX ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isSizet - Is value representable without loss as a size_t?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isSizet( Integer value ) {
  const Integer result = AND2( value >= 0,
                               OR2( sizeof (size_t) >= sizeof (Integer),
                                    value <= SIZET_MAX ) );
  POST0( result == AND2( value >= 0,
                         OR2( sizeof (size_t) >= sizeof (Integer),
                              value <= SIZET_MAX ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isSignedLongLong - Is value representable without loss as a
         signed long long?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isSignedLongLong( Integer value ) {
  const Integer result = IN_RANGE( value, LONGLONG_MIN, LONGLONG_MAX );
  POST0( result == IN_RANGE( value, LONGLONG_MIN, LONGLONG_MAX ) );
  return result;
}



/******************************************************************************
PURPOSE: isUnsignedLongLong - Is value representable without loss as an
         unsigned long long?
INPUTS:  Integer value - The 64-bit value to test.
RETURNS: Integer 1 if representable, else 0.
******************************************************************************/

Integer isUnsignedLongLong( Integer value ) {
  const Integer result =
    AND2( value >= 0,
          OR2( sizeof (unsigned long long) >= sizeof (Integer),
               value <= ULONGLONG_MAX ) );
  POST0( result == AND2( value >= 0,
                         OR2( sizeof (unsigned long long) >= sizeof (Integer),
                              value <= ULONGLONG_MAX ) ) );
  return result;
}



/******************************************************************************
PURPOSE: withinTolerance - Do x and y differ by less than (non-negative,
         finite) tolerance (i.e., fabs( x - y ) <= tolerance) or, for large
         values, differ only in digits beyond number of significant digits in
         tolerance? (E.g., tolerance = 1e-6 and 1.0000000001e30 and
         1.0000000002e30 are considered equal.)
INPUTS:  Real x - first value to compare.
         Real y - second value to compare.
         Real tolerance - tolerance threshold (e.g., 1-e6).
RETURNS: Integer 1 if x and y differ (in significant digits) by less than
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

Integer withinTolerance( Real x, Real y, Real tolerance ) {
  PRE02( ! isNan( tolerance ), tolerance <= 0.1 );

  /* First try bitwise comparison (handles nans): */

  const Real* const xp = &x;
  const Real* const yp = &y;
  const Integer* const ixp = (const Integer*) xp;
  const Integer* const iyp = (const Integer*) yp;
  const Integer ix = *ixp;
  const Integer iy = *iyp;
  Integer result = ix == iy;

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
      const Real ax = fabs( x );
      const Real ay = fabs( y );

      if ( AND2( ay < 1.0, ax > ay * DBL_MAX ) ) { /* Avoid overflow. */
        result = 0;
      } else if ( AND2( ay > 1.0, ax < ay * DBL_MIN ) ) { /* Avoid underflow.*/
        result = 0;
      } else {
        const Real ratio = x / y;
        result = IN_RANGE( ratio, 1.0 - tolerance, 1.0 + tolerance );
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: aboutEqual - Is withinTolerance( x, y, TOLERANCE )?
INPUTS:  Real x - first value to compare.
         Real y - second value to compare.
RETURNS: Integer 1 if x and y differ (in significant digits) by less than
         TOLERANCE, else 0.
******************************************************************************/

Integer aboutEqual( Real x, Real y ) {
  const Integer result = withinTolerance( x, y, TOLERANCE );
  POST0( result == withinTolerance( x, y, TOLERANCE ) );
  return result;
}



/******************************************************************************
PURPOSE: isNan - Is x a NaN (Not a Number)?
INPUTS:  Real x - The 64-bit value to test.
RETURNS: Integer 1 if x is a NaN, else 0.
******************************************************************************/

Integer isNan( Real x ) {
  const Real copy = x;
  const Integer result = (copy != x);
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: filterNan - If x is a NaN return 0 else return x.
INPUTS:  Real x - The 64-bit value to test.
RETURNS: Real 0 if isNan( x ), else x.
******************************************************************************/

Real filterNan( Real x ) {
  const Real result = isNan( x ) ? 0.0 : x;
  POST0( IMPLIES_ELSE( isNan( x ), result == 0.0, result == x ) );
  return result;
}



/******************************************************************************
PURPOSE: absI - Absolute-value of x (or INTEGER_MAX if x == INTEGER_MIN).
INPUTS:  Integer x - The 64-bit value to test.
RETURNS: Integer |x|.
******************************************************************************/

Integer absI( Integer x ) {
  const Integer result = x >= 0 ? x : x == INTEGER_MIN ? INTEGER_MAX : -x;
  POST02( result >= 0,
          IMPLIES_ELSE( x >= 0,
                        result == x,
                        IMPLIES_ELSE( x == INTEGER_MIN,
                                      result == INTEGER_MAX,
                                      result == -x ) ) );
  return result;
}



/******************************************************************************
PURPOSE: toInteger - Integer value of string if within range [lower, upper].
INPUTS:  const char* string  - The string to convert.
         Integer lower  - The lower limit of valid range.
         Integer upper  - The upper limit of valid range.
OUTPUTS: Integer* ok    - Does string represent an integer in [lower,upper]?
RETURNS: Integer value of string within range [lower, upper], else 0.
NOTES:   Unlike atoI(), strtoI() and scanf(), this routine rejects strings
         that would overflow or that contain non-digit characters (except
         an optional leading sign) or lack digit characters, or contain
         multiple whitespace-separated tokens.
******************************************************************************/

Integer toInteger( const char* string, Integer lower, Integer upper,
                   Integer* ok ) {
  PRE03( string, lower <= upper, ok );
  Integer result = 0;
  *ok = 0;
  errno = 0; /*strtoll sets errno upon failure but won't clear it on success!*/

  if ( *string ) {
    char* terminator = 0;
    const Integer convertionResult = strtoll( string, &terminator, 10 );

    if ( AND2( terminator, terminator != string ) ) {

      while ( isspace( *terminator ) ) { /* Skip trailing whitespace. */
        ++terminator;
      }

      *ok = AND3( errno == 0, /* No overflow. */
                  *terminator == '\0', /* No remaining characters. */
                  IN_RANGE( convertionResult, lower, upper ) );

      if ( *ok ) {
        result = convertionResult;
      }
    }
  }

  POST0( IMPLIES_ELSE( *ok, IN_RANGE( result, lower, upper ), result == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: toReal - Real value of string if within range [lower, upper].
INPUTS:  const char* string  - The string to convert.
         Real lower  - The lower limit of valid range.
         Real upper  - The upper limit of valid range.
OUTPUTS: Integer* ok    - Does string represent an integer in [lower,upper]?
RETURNS: Real value of string within range [lower, upper], else 0.0.
******************************************************************************/

Real toReal( const char* string, Real lower, Real upper, Integer* ok ) {
  PRE05( string, ! isNan( lower ), ! isNan( upper ), lower <= upper, ok );
  Real result = 0;
  *ok = 0;
  errno = 0; /*strtod sets errno upon failure but won't clear it on success!*/

  if ( *string ) {
    char* terminator = 0;
    const Real convertionResult = strtod( string, &terminator );

    if ( AND2( terminator, terminator != string ) ) {

      while ( isspace( *terminator ) ) { /* Skip trailing whitespace. */
        ++terminator;
      }

      *ok = AND2( *terminator == '\0', /* No remaining characters. */
                  IN_RANGE( convertionResult, lower, upper ) );

      if ( *ok ) {
        result = convertionResult;
      }
    }
  }

  POST0( IMPLIES_ELSE( *ok, IN_RANGE( result, lower, upper ), result == 0.0 ));
  return result;
}



/******************************************************************************
PURPOSE: isInfinity - Is the value so large that its recipricol is zero?
INPUTS:  Real x - The value to test.
RETURNS: Integer 1 if recipricol is zero, else 0.
******************************************************************************/

Integer isInfinity( Real x ) {
  const Integer result =
    AND3( ! isNan( x ), x > 0.0, OR2( x > REAL_MAX, 1.0 / x == 0.0 ) );
  POST0( result == AND3( ! isNan( x ),
                         x > 0.0,
                         OR2( x > REAL_MAX, 1.0 / x == 0.0 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isMinusInfinity - Is the value so small that its recipricol is zero?
INPUTS:  Real x - The value to test.
RETURNS: Integer 1 if recipricol is zero, else 0.
******************************************************************************/

Integer isMinusInfinity( Real x ) {
  const Integer result =
    AND3( ! isNan( x ), x < 0.0, OR2( x < -REAL_MAX, 1.0 / x == 0.0 ) );
  POST0( result == AND3( ! isNan( x ),
                         x < 0.0,
                         OR2( x < -REAL_MAX, 1.0 / x == 0.0 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isFinite - Is the value and its recipricol
         non-NaN/+Infinity/-Infinity?
INPUTS:  Real x - The value to test.
RETURNS: Integer 1 if finite, else 0.
******************************************************************************/

Integer isFinite( Real x ) {
  const Integer result =
    AND3( ! isNan( x ), ! isInfinity( x ), ! isMinusInfinity( x ) );
  POST0( result == AND3( ! isNan(x), ! isInfinity(x), ! isMinusInfinity(x) ));
  return result;
}



/******************************************************************************
PURPOSE: radians - Radians of degrees.
INPUTS:  Real theDegrees - The degrees to convert.
RETURNS: Real radians of degrees.
******************************************************************************/

Real radians( Real theDegrees ) {
  PRE0( ! isNan( theDegrees ) );
  const Real result = theDegrees * ( M_PI / 180.0 );
  POST03( ! isNan( result ),
          OR2( SIGN( result ) == SIGN( theDegrees ), result == 0.0 ),
          fabs( result ) <= fabs( theDegrees ) );
  return result;
}



/******************************************************************************
PURPOSE: degrees - Degrees of radians.
INPUTS:  Real theDegrees - The degrees to convert.
RETURNS: Real radians of degrees.
******************************************************************************/

Real degrees( Real theRadians ) {
  PRE0( ! isNan( theRadians ) );
  const Real result = theRadians * ( 180.0 / M_PI );
  POST03( ! isNan( result ),
          OR2( SIGN( result ) == SIGN( theRadians ), result == 0.0 ),
          fabs( result ) >= fabs( theRadians ) );
  return result;
}



/******************************************************************************
PURPOSE: safeSum - NaN-free sum.
INPUTS:  Real x - The first value to sum.
         Real y - The second value to sum.
RETURNS: Real sum of arguments.
******************************************************************************/

Real safeSum( Real x, Real y ) {
  PRE02( ! isNan( x ), ! isNan( y ) );
  const Real result = x == -y ? 0.0 : x + y;
  POST02( ! isNan( result ), IMPLIES( x == -y, result == 0.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: safeSum3 - NaN-free sum.
INPUTS:  Real a1 - The first value to sum.
         Real a2 - The second value to sum.
         Real a3 - The third value to sum.
RETURNS: Real sum of arguments.
******************************************************************************/

Real safeSum3( Real a1, Real a2, Real a3 ) {
  PRE03( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ) );
  const Real result = safeSum( safeSum( a1, a2 ), a3 );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeSum4 - NaN-free sum.
INPUTS:  Real a1 - The first value to sum.
         Real a2 - The second value to sum.
         Real a3 - The third value to sum.
         Real a4 - The fourth value to sum.
RETURNS: Real sum of arguments.
******************************************************************************/

Real safeSum4( Real a1, Real a2, Real a3, Real a4 ) {
  PRE04( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ) );
  const Real result = safeSum( safeSum( a1, a2 ), safeSum( a3, a4 ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeSum5 - NaN-free sum.
INPUTS:  Real a1 - The first value to sum.
         Real a2 - The second value to sum.
         Real a3 - The third value to sum.
         Real a4 - The fourth value to sum.
         Real a5 - The fifth value to sum.
RETURNS: Real sum of arguments.
******************************************************************************/

Real safeSum5( Real a1, Real a2, Real a3, Real a4, Real a5 ) {
  PRE05( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ) );
  const Real result = safeSum( safeSum( a1, a2 ),
                               safeSum( safeSum( a3, a4 ), a5 ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeSum6 - NaN-free sum.
INPUTS:  Real a1 - The first value to sum.
         Real a2 - The second value to sum.
         Real a3 - The third value to sum.
         Real a4 - The fourth value to sum.
         Real a5 - The fifth value to sum.
         Real a6 - The sixth value to sum.
RETURNS: Real sum of arguments.
******************************************************************************/

Real safeSum6( Real a1, Real a2, Real a3, Real a4, Real a5, Real a6 ) {
  PRE06( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ), ! isNan( a6 ) );
  const Real result = safeSum( safeSum( a1, a2 ),
                               safeSum( safeSum( a3, a4 ),
                                        safeSum( a5, a6 ) ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeSum7 - NaN-free sum.
INPUTS:  Real a1 - The first value to sum.
         Real a2 - The second value to sum.
         Real a3 - The third value to sum.
         Real a4 - The fourth value to sum.
         Real a5 - The fifth value to sum.
         Real a6 - The sixth value to sum.
         Real a7 - The seventh value to sum.
RETURNS: Real sum of arguments.
******************************************************************************/

Real safeSum7( Real a1, Real a2, Real a3, Real a4, Real a5, Real a6,
                      Real a7 ) {
  PRE07( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ), ! isNan( a6 ), ! isNan( a7 ) );
  const Real result = safeSum( safeSum( a1, a2 ),
                               safeSum( safeSum( a3, a4 ),
                                        safeSum( safeSum( a5, a6 ), a7 ) ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeSum8 - NaN-free sum.
INPUTS:  Real a1 - The first value to sum.
         Real a2 - The second value to sum.
         Real a3 - The third value to sum.
         Real a4 - The fourth value to sum.
         Real a5 - The fifth value to sum.
         Real a6 - The sixth value to sum.
         Real a7 - The seventh value to sum.
         Real a8 - The eighth value to sum.
RETURNS: Real sum of arguments.
******************************************************************************/

Real safeSum8( Real a1, Real a2, Real a3, Real a4, Real a5, Real a6,
                      Real a7, Real a8 ) {
  PRE08( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ), ! isNan( a6 ), ! isNan( a7 ), ! isNan( a8 ) );
  const Real result = safeSum( safeSum( a1, a2 ),
                               safeSum( safeSum( safeSum( a3, a4 ), a5 ),
                                        safeSum( safeSum( a6, a7 ), a8 ) ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeDifference - NaN-free difference.
INPUTS:  Real left  - The left value of difference.
         Real right - The right value of difference.
RETURNS: Real difference of arguments.
******************************************************************************/

Real safeDifference( Real x, Real y ) {
  PRE02( ! isNan( x ), ! isNan( y ) );
  const Real result = x == y ? 0.0 : x - y;
  POST02( ! isNan( result ), IMPLIES( x == y, result == 0.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: safeProduct - NaN-free product.
INPUTS:  Real x - The first value of product.
         Real y - The second value of product.
RETURNS: Real product of arguments.
******************************************************************************/

Real safeProduct( Real x, Real y ) {
  PRE02( ! isNan( x ), ! isNan( y ) );
  const Real result = ( x == 0.0 || y == 0.0 ) ? 0.0 : x * y;
  POST02( ! isNan( result ),
          IMPLIES( OR2( x == 0.0, y == 0.0 ), result == 0.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: safeProduct3 - NaN-free product.
INPUTS:  Real a1 - The first value to product.
         Real a2 - The second value to product.
         Real a3 - The third value to product.
RETURNS: Real product of arguments.
******************************************************************************/

Real safeProduct3( Real a1, Real a2, Real a3 ) {
  PRE03( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ) );
  const Real result = safeProduct( safeProduct( a1, a2 ), a3 );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeProduct4 - NaN-free product.
INPUTS:  Real a1 - The first value to product.
         Real a2 - The second value to product.
         Real a3 - The third value to product.
         Real a4 - The fourth value to product.
RETURNS: Real product of arguments.
******************************************************************************/

Real safeProduct4( Real a1, Real a2, Real a3, Real a4 ) {
  PRE04( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ) );
  const Real result = safeProduct( safeProduct( a1, a2 ),
                                   safeProduct( a3, a4 ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeProduct5 - NaN-free product.
INPUTS:  Real a1 - The first value to product.
         Real a2 - The second value to product.
         Real a3 - The third value to product.
         Real a4 - The fourth value to product.
         Real a5 - The fifth value to product.
RETURNS: Real product of arguments.
******************************************************************************/

Real safeProduct5( Real a1, Real a2, Real a3, Real a4, Real a5 ) {
  PRE05( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ) );
  const Real result = safeProduct( safeProduct( a1, a2 ),
                                   safeProduct( safeProduct( a3, a4 ), a5 ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeProduct6 - NaN-free product.
INPUTS:  Real a1 - The first value to product.
         Real a2 - The second value to product.
         Real a3 - The third value to product.
         Real a4 - The fourth value to product.
         Real a5 - The fifth value to product.
         Real a6 - The sixth value to product.
RETURNS: Real product of arguments.
******************************************************************************/

Real safeProduct6( Real a1, Real a2, Real a3, Real a4, Real a5,
                          Real a6 ) {
  PRE06( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ), ! isNan( a6 ) );
  const Real result = safeProduct( safeProduct( a1, a2 ),
                                   safeProduct( safeProduct( a3, a4 ),
                                                safeProduct( a5, a6 ) ) );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeProduct7 - NaN-free product.
INPUTS:  Real a1 - The first value to product.
         Real a2 - The second value to product.
         Real a3 - The third value to product.
         Real a4 - The fourth value to product.
         Real a5 - The fifth value to product.
         Real a6 - The sixth value to product.
         Real a7 - The seventh value to product.
RETURNS: Real product of arguments.
******************************************************************************/

Real safeProduct7( Real a1, Real a2, Real a3, Real a4, Real a5,
                          Real a6, Real a7 ) {
  PRE07( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ), ! isNan( a6 ), ! isNan( a7 ) );
  const Real result = safeProduct( safeProduct( a1, a2 ),
                                   safeProduct( safeProduct( a3, a4 ),
                                                safeProduct3( a5, a6, a7 ) ));
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeProduct8 - NaN-free product.
INPUTS:  Real a1 - The first value to product.
         Real a2 - The second value to product.
         Real a3 - The third value to product.
         Real a4 - The fourth value to product.
         Real a5 - The fifth value to product.
         Real a6 - The sixth value to product.
         Real a7 - The seventh value to product.
         Real a8 - The eighth value to product.
RETURNS: Real product of arguments.
******************************************************************************/

Real safeProduct8( Real a1, Real a2, Real a3, Real a4, Real a5,
                          Real a6, Real a7, Real a8 ) {
  PRE08( ! isNan( a1 ), ! isNan( a2 ), ! isNan( a3 ), ! isNan( a4 ),
         ! isNan( a5 ), ! isNan( a6 ), ! isNan( a7 ), ! isNan( a8 ) );
  const Real result = safeProduct( safeProduct( a1, a2 ),
                                   safeProduct( safeProduct3( a3, a4, a5 ),
                                                safeProduct3( a6, a7, a8 ) ));
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: safeQuotient - NaN-free quotient.
INPUTS:  Real numerator - The numerator value of ratio.
         Real denominator - The denominator value of ratio.
RETURNS: Real quotient of arguments.
******************************************************************************/

Real safeQuotient( Real numerator, Real denominator ) {
  PRE03( ! isNan( numerator ), ! isNan( denominator ), denominator != 0.0 );
  const Real result =   numerator   ==  0.0 ? 0.0
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
PURPOSE: kahanSum - Kahan summation of array elements.
INPUTS:  const Real* items - The array of items to sum.
         Integer count - The number of items in the array.
RETURNS: Real The minimal-round-off-error sum of the array items.
******************************************************************************/

Real kahanSum( const Real* items, Integer count ) {
  PRE03( items,
         IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ),
         isNanFree( items, count ) );
  Real result     = 0.0; /* The sum of the values in the array. */
  Real correction = 0.0; /* Kahan corrector subtracts each round-off error. */
  Integer index;         /* The index on each value in the array. */

  for ( index = 0; index < count; ++index ) {
    const Real nextTerm           = items[ index ];
    const Real correctedNextTerm  = nextTerm - correction;
    const Real newSum             = result + correctedNextTerm;
    correction = ( newSum - result ) - correctedNextTerm;
    result = newSum;
  }

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sumI - Summation of array elements.
INPUTS:  const Integer* items - The array of items to sum.
         Integer count - The number of items in the array.
RETURNS: Integer The sum of the array items.
******************************************************************************/

Integer sumI( const Integer* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer result = 0; /* The sum of the values in the array. */
  Integer index;      /* The index on each value in the array. */

  for ( index = 0; index < count; ++index ) {
    result += items[ index ];
  }

  return result;
}



/******************************************************************************
PURPOSE: filterNans - Replace NaN array items with zero.
INPUTS:  Real*   items - The array of items to check.
         Integer count - The number of items in the array.
OUTPUTS: Real*   items - The NaN-free array items.
******************************************************************************/

void filterNans( Real* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer index;

  for ( index = 0; index < count; ++index ) {
    items[ index ] = filterNan( items[ index ] );
  }

  POST0( isNanFree( items, count ) );
}



/******************************************************************************
PURPOSE: isNanFree - Verify that an array contains no NaNs (Not A Number).
INPUTS:  const Real* items - The array of items to check.
         Integer     count - The number of items in the array.
RETURNS: Integer 1 if no NaNs found, else 0.
******************************************************************************/

Integer isNanFree( const Real* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer result = 1;
  Integer index;

  for ( index = 0; AND2( result, index < count ); ++index ) {
    result = ! isNan( items[ index ] );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: allFinite - Are all array items finite?
INPUTS:  const Real* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Integer 1 if all items are finite, else 0.
******************************************************************************/

Integer allFinite( const Real* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer result = 1;
  Integer index;

  for ( index = 0; AND2( result, index < count ); ++index ) {
    result = isFinite( items[ index ] );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: increasing - Are all array items increasing?
INPUTS:  const Real* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Integer 1 if all items are increasing, else 0.
******************************************************************************/

Integer increasing( const Real* items, Integer count ) {
  PRE03( items,
         IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ),
         isNanFree( items, count ) );
  Integer result = 1;
  Integer index;

  for ( index = 1; AND2( result, index < count ); ++index ) {
    result = items[ index - 1 ] < items[ index ];
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: decreasing - Are all array items decreasing?
INPUTS:  const Real* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Integer 1 if all items are decreasing, else 0.
******************************************************************************/

Integer decreasing( const Real* items, Integer count ) {
  PRE03( items,
         IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ),
         isNanFree( items, count ) );
  Integer result = 1;
  Integer index;

  for ( index = 1; AND2( result, index < count ); ++index ) {
    result = items[ index - 1 ] > items[ index ];
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: allZero - Are all array items zero?
INPUTS:  const Real* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Integer 1 if all items are zero, else 0.
******************************************************************************/

Integer allZero( const Real* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer result = 1;
  Integer index;

  for ( index = 0; AND2( result, index < count ); ++index ) {
    result = items[ index ] == 0.0;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: allZeroI - Are all array items zero?
INPUTS:  const Integer* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Integer 1 if all items are zero, else 0.
******************************************************************************/

Integer allZeroI( const Integer* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer result = 1;
  Integer index;

  for ( index = 0; AND2( result, index < count ); ++index ) {
    result = items[ index ] == 0;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: maximumItem - Largest value in the array.
INPUTS:  const Real* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Real largest value in the array.
******************************************************************************/

Real maximumItem( const Real* items, Integer count ) {
  PRE03( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ),
         isNanFree( items, count ) );
  Real result = items[ 0 ];
  Integer index;

  for ( index = 1; AND2( result <= REAL_MAX, index < count ); ++index ) {
    const Real value = items[ index ];

    if ( value > result ) {
      result = value;
    }
  }

  POST03( ! isNan( result ),
          result >= items[ 0 ],
          result >= items[ count - 1 ] );
  return result;
}



/******************************************************************************
PURPOSE: minimumItem - Smallest value in the array.
INPUTS:  const Real* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Real smallest value in the array.
******************************************************************************/

Real minimumItem( const Real* items, Integer count ) {
  PRE03( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ),
         isNanFree( items, count ) );
  Real result = items[ 0 ];
  Integer index;

  for ( index = 1; AND2( result >= -REAL_MAX, index < count ); ++index ) {
    const Real value = items[ index ];

    if ( value < result ) {
      result = value;
    }
  }

  POST03( ! isNan( result ),
          result <= items[ 0 ],
          result <= items[ count - 1 ] );
  return result;
}



/******************************************************************************
PURPOSE: maximumItemI - Largest value in the array.
INPUTS:  const Integer* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Integer largest value in the array.
******************************************************************************/

Integer maximumItemI( const Integer* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer result = items[ 0 ];
  Integer index;

  for ( index = 1; AND2( result < INTEGER_MAX, index < count ); ++index ) {
    const Integer value = items[ index ];

    if ( value > result ) {
      result = value;
    }
  }

  POST02( result >= items[ 0 ], result >= items[ count - 1 ] );
  return result;
}



/******************************************************************************
PURPOSE: minimumItemI - Smallest value in the array.
INPUTS:  const Integer* items - The array of items to check.
         Integer count - The number of items in the array.
RETURNS: Real smallest value in the array.
******************************************************************************/

Integer minimumItemI( const Integer* items, Integer count ) {
  PRE02( items, IN_RANGE( count, 1, INTEGER_MAX / sizeof items[ 0 ] ) );
  Integer result = items[ 0 ];
  Integer index;

  for ( index = 1; AND2( result > INTEGER_MIN, index < count ); ++index ) {
    const Integer value = items[ index ];

    if ( value < result ) {
      result = value;
    }
  }

  POST02( result <= items[ 0 ], result <= items[ count - 1 ] );
  return result;
}



/******************************************************************************
PURPOSE: reverseItems - Reverse the order of the items in an array.
INPUTS:  Real    array[ count ] - Array of values.
         Integer count          - Number of elements in the array.
OUTPUTS: Real    array[ count ] - Array of values in reverse order.
******************************************************************************/

void reverseItems( Real array[], Integer count ) {
  PRE02( array, IN_RANGE( count, 1, INTEGER_MAX / sizeof array[ 0 ] ) );
  const Integer halfCount = count / 2;
  const Integer count1 = count - 1;
  Integer index = 0;

  for ( index = 0; index < halfCount; ++index ) {
    const Integer other = count1 - index;
    const Real temp = array[ index ];
    array[ index ] = array[ other ];
    array[ other ] = temp;
  }

}



/******************************************************************************
PURPOSE: fillRandom - Initialize an array with a pseudo-random sequence of
         Reals uniformly distributed in the range [0, 1].
INPUTS:  Integer count          - Number of elements in the array.
OUTPUTS: Real    array[ count ] - Array of pseudo-random values.
******************************************************************************/

void fillRandom( Real array[], Integer count ) {
  PRE02( array, IN_RANGE( count, 1, INTEGER_MAX / sizeof array[ 0 ] ) );
  Integer index;

  for ( index = 0; index < count; ++index ) {
    array[ index ] = drand48();
  }

  POST03( isNanFree( array, count ),
          minimumItem( array, count ) >= 0,
          maximumItem( array, count ) <= 1.0 );
}



/******************************************************************************
PURPOSE: fillRandomI - Initialize an array with a pseudo-random sequence of
         Integers uniformly distributed in the range [low, high].
INPUTS:  Integer count          - Number of elements in the array.
         Integer low            - Minimum value to generate.
         Integer high           - Maximum value to generate.
OUTPUTS: Integer array[ count ] - Array of pseudo-random values.
******************************************************************************/

void fillRandomI(Integer array[], Integer count, Integer low, Integer high) {
  PRE04( array,
         IN_RANGE( count, 1, INTEGER_MAX / sizeof array[ 0 ] ),
         high - low >= 0,
         high - low < INTEGER_MAX );
  Integer index;

  for ( index = 0; index < count; ++index ) {
    array[ index ] = randomInteger( low, high );
  }

  POST02( minimumItemI( array, count ) >= low,
          maximumItemI( array, count ) <= high );
}



/******************************************************************************
PURPOSE: shuffle - Pseudo-randomly permute (reorder) an array's items.
INPUTS:  Integer count          - Number of elements in the array.
         Real    array[ count ] - Array of values.
OUTPUTS: Real    array[ count ] - Array of shuffled values.
******************************************************************************/

void shuffle( Real array[], Integer count ) {
  PRE03( array, IN_RANGE( count, 1, INTEGER_MAX / sizeof array[ 0 ] ),
         isNanFree( array, count ) );
  CHECKING( const Real OLD( minimumItem ) = minimumItem( array, count ); )
  CHECKING( const Real OLD( maximumItem ) = maximumItem( array, count ); )
  const Integer last = count - 1;
  Integer index;

  for ( index = 0; index < last; ++index ) {
    const Integer other = randomInteger( index + 1, last );
    const Real swapped = array[ index ];
    CHECK2( other != index, IN_RANGE( other, 0, count - 1 ) );
    array[ index ] = array[ other ];
    array[ other ] = swapped;
  }

  POST02( minimumItem( array, count ) == OLD( minimumItem ),
          maximumItem( array, count ) == OLD( maximumItem ) );
}



/******************************************************************************
PURPOSE: shuffleI - Pseudo-randomly permute (reorder) an array's items.
INPUTS:  Integer count          - Number of elements in the array.
         Integer array[ count ] - Array of values.
OUTPUTS: Integer array[ count ] - Array of shuffled values.
******************************************************************************/

void shuffleI( Integer array[], Integer count ) {
  PRE02( array, IN_RANGE( count, 1, INTEGER_MAX / sizeof array[ 0 ] ) );
  CHECKING( const Integer OLD( minimumItem ) = minimumItemI( array, count ); )
  CHECKING( const Integer OLD( maximumItem ) = maximumItemI( array, count ); )
  const Integer last = count - 1;
  Integer index;

  for ( index = 0; index < last; ++index ) {
    const Integer other = randomInteger( index + 1, last );
    const Integer swapped = array[ index ];
    CHECK2( other != index, IN_RANGE( other, 0, count - 1 ) );
    array[ index ] = array[ other ];
    array[ other ] = swapped;
  }

  POST02( minimumItemI( array, count ) == OLD( minimumItem ),
          maximumItemI( array, count ) == OLD( maximumItem ) );
}



/******************************************************************************
PURPOSE: isSorted - Determine if an array of Reals is sorted in ascending
         order.
INPUTS:  const Real array[ count ] - Array of values to check.
         Integer    count          - Number of elements in the array.
RETURNS: Integer 1 if sorted, else 0.
******************************************************************************/

Integer isSorted( const Real array[], Integer count ) {
  PRE03( array,
         IN_RANGE( count, 1, INTEGER_MAX / sizeof array[ 0 ] ),
         isNanFree( array, count ) );
  Integer result = 1;
  Integer index;

  for ( index = 1; AND2( result, index < count ); ++index ) {
    result = array[ index - 1 ] <= array[ index ];
  }

  POST03( IS_BOOL( result ),
          IMPLIES( result, array[ 0 ] <= array[ count - 1 ] ),
          IMPLIES( array[ 0 ] > array[ count - 1 ], ! result ) );
  return result;
}



/******************************************************************************
PURPOSE: isSortedI - Determine if an array of Integers is sorted in ascending
         order.
INPUTS:  const Integer array[ count ] - Array of values to check.
         Integer       count          - Number of elements in the array.
RETURNS: Integer 1 if sorted, else 0.
******************************************************************************/

Integer isSortedI( const Integer array[], Integer count ) {
  PRE02( array, IN_RANGE( count, 1, INTEGER_MAX / sizeof array[ 0 ] ) );
  Integer result = 1;
  Integer index;

  for ( index = 1; AND2( result, index < count ); ++index ) {
    result = array[ index - 1 ] <= array[ index ];
  }

  POST03( IS_BOOL( result ),
          IMPLIES( result, array[ 0 ] <= array[ count - 1 ] ),
          IMPLIES( array[ 0 ] > array[ count - 1 ], ! result ) );
  return result;
}



/******************************************************************************
PURPOSE: shellsort - Sort an array of Reals using the Shellsort algorithm.
INPUTS:  Real    array[ count ] - Array of values to sort.
         Integer count          - Number of elements in the array.
OUTPUTS: Real    array[ count ] - Array of values sorted in ascending order.
NOTES:   Shellsort sorts 'count' items using the diminishing increment sequence
         [6078832729528464400, ..., 121, 40, 13, 4, 1] resulting in about
         'count' ^ 1.25 comparisons < count * log2(count) for count < 65,536
         - i.e., theoretically more efficient than Quicksort for up to 2^16
         items. In practice, processor time is less than a qsort()-based
         implementation for up to about 10^7 items. Also, Shellsort is stable
         (unlike Quicksort) - meaning items that compare equal remain in their
         original relative order.
         See Sedgewick, Robert, "Algorithms", Second Edition, pp 107-111,
         Addision Wesley, Reading MA, 1988.
******************************************************************************/

void shellsort( Real array[], Integer count ) {
  PRE03( array,
         IN_RANGE( count, 1,
                          MIN( INTEGER_MAX / sizeof array[ 0 ],
                               INTEGER_CONSTANT(6078832729528464400) - 1 ) ),
         isNanFree( array, count ) );

  Integer h;

  /*
   * Initialize diminishing increment 'h' to value in sequence
   * [..., 121, 40, 13, 4, 1] that bounds 'count' and results in
   * about 'count' ^ 1.25 number of comparisons.
   * This takes just log3(count) iterations.
   */

  for ( h = 4; h <= count; h = 3 * h + 1 ) {
    CHECK( h > 0 ); /* next */
  }

  do {
    Integer i;
    h /= 3;
    CHECK( h > 0 );

    for ( i = h; i < count; ++i ) {
      const Real value = array[ i ];
      Integer j   = i;
      Integer jh  = j - h;
      Real    ajh = array[ jh ];
      CHECK2( IN_RANGE( i, 1, count - 1 ), IN_RANGE( jh, 0, count - 1 ) );

      while ( value < ajh ) {
        CHECK( IN_RANGE( j, 0, count - 1 ) );
        array[ j ] = ajh;
        j = jh;

        if ( j < h ) {
          ajh = value; /* Done. */
        } else {
          jh = j - h;
          ajh = array[ jh ];
        }
      }

      CHECK( IN_RANGE( j, 0, count - 1 ) );
      array[ j ] = value;
    }

  } while ( h != 1 );

  POST0( isSorted( array, count ) );
}



/******************************************************************************
PURPOSE: shellsortI - Sort an array of Integers using the Shellsort algorithm.
INPUTS:  Integer array[ count ] - Array of values to sort.
         Integer count          - Number of elements in the array.
OUTPUTS: Integer array[ count ] - Array of values sorted in ascending order.
NOTES:   Shellsort sorts 'count' items using the diminishing increment sequence
         [6078832729528464400, ..., 121, 40, 13, 4, 1] resulting in about
         'count' ^ 1.25 comparisons < count * log2(count) for count < 65,536
         - i.e., theoretically more efficient than Quicksort for up to 2^16
         items. In practice, processor time is less than a qsort()-based
         implementation for up to about 10^7 items. Also, Shellsort is stable
         (unlike Quicksort) - meaning items that compare equal remain in their
         original relative order.
         See Sedgewick, Robert, "Algorithms", Second Edition, pp 107-111,
         Addision Wesley, Reading MA, 1988.
******************************************************************************/

void shellsortI( Integer array[], Integer count ) {
  PRE02( array,
         IN_RANGE( count, 1,
                          MIN( INTEGER_MAX / sizeof array[ 0 ],
                               INTEGER_CONSTANT(6078832729528464400) - 1 ) ));
  Integer h;

  /*
   * Initialize diminishing increment 'h' to value in sequence
   * [..., 121, 40, 13, 4, 1] that bounds 'count' and results in
   * about 'count' ^ 1.25 number of comparisons.
   * This takes just log3(count) iterations.
   */

  for ( h = 4; h <= count; h = 3 * h + 1 ) {
    CHECK( h > 0 ); /* next */
  }

  do {
    Integer i;
    h /= 3;
    CHECK( h > 0 );

    for ( i = h; i < count; ++i ) {
      const Integer value = array[ i ];
      Integer j   = i;
      Integer jh  = j - h;
      Integer ajh = array[ jh ];
      CHECK2( IN_RANGE( i, 1, count - 1 ), IN_RANGE( jh, 0, count - 1 ) );

      while ( value < ajh ) {
        CHECK( IN_RANGE( j, 0, count - 1 ) );
        array[ j ] = ajh;
        j = jh;

        if ( j < h ) {
          ajh = value; /* Done. */
        } else {
          jh = j - h;
          ajh = array[ jh ];
        }
      }

      CHECK( IN_RANGE( j, 0, count - 1 ) );
      array[ j ] = value;
    }

  } while ( h != 1 );

  POST0( isSortedI( array, count ) );
}



/******************************************************************************
PURPOSE: randomInteger() - Yield a pesudo-randomly selected integer within the
         range [low, high].
INPUTS:  Integer low  Minimum value of range.
         Integer high Maximum value of range.
RETURNS: Integer within the range [low, high].
******************************************************************************/

Integer randomInteger( Integer low, Integer high ) {
  PRE02( high - low >= 0, high - low < INTEGER_MAX );
  const Integer range = high - low;
  const Integer result = drand48() * range + low;
  POST0( IN_RANGE( result, low, high ) );
  return result;
}



/******************************************************************************
PURPOSE: rotate4ByteWordIfLittleEndian() - Rotate 4-bytes of value if on a
         little-endian platform.
INPUTS:  void* value  4-byte value to rotate.
OUTPUTS: void* value  Rotated value.
******************************************************************************/

void rotate4ByteWordIfLittleEndian( void* value ) {

#if IS_LITTLE_ENDIAN

  PRE0( value );
  assert_static( sizeof (int) == 4 );
  int* const ivalue = value;
  const int value4 = *ivalue;
  const int result =
    ( value4 & 0xff000000 ) >> 24 |
    ( value4 & 0x00ff0000 ) >>  8 |
    ( value4 & 0x0000ff00 ) <<  8 |
    ( value4 & 0x000000ff ) << 24;
  *ivalue = result;

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: rotate4ByteArrayIfLittleEndian() - Rotate 4-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array    Array of 4-byte values to rotate.
         Integer count  Number of items in array.
OUTPUTS: void* array    Array of rotated values.
******************************************************************************/

void rotate4ByteArrayIfLittleEndian( void* array, Integer count ) {

#if IS_LITTLE_ENDIAN

  PRE02( array, count > 0 );
  assert_static( sizeof (int) == 4 );
  int* const array4 = array;
  Integer index = 0;

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const int value = array4[ index ];
    const int newValue =
      ( value & 0xff000000 ) >> 24 |
      ( value & 0x00ff0000 ) >>  8 |
      ( value & 0x0000ff00 ) <<  8 |
      ( value & 0x000000ff ) << 24;
    array4[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: rotate8ByteWordIfLittleEndian() - Rotate 8-bytes of value if on a
         little-endian platform.
INPUTS:  void* value  8-byte value to rotate.
OUTPUTS: void* value  Rotated value.
******************************************************************************/

void rotate8ByteWordIfLittleEndian( void* value ) {

#if IS_LITTLE_ENDIAN

  PRE0( value );
  assert_static( sizeof (Integer) == 8 );
  Integer* const ivalue = value;
  const Integer value8 = *ivalue;
  const Integer result =
    ( value8 & INTEGER_CONSTANT( 0xff00000000000000 ) ) >> 56 |
    ( value8 & INTEGER_CONSTANT( 0x00ff000000000000 ) ) >> 40 |
    ( value8 & INTEGER_CONSTANT( 0x0000ff0000000000 ) ) >> 24 |
    ( value8 & INTEGER_CONSTANT( 0x000000ff00000000 ) ) >>  8 |
    ( value8 & INTEGER_CONSTANT( 0x00000000ff000000 ) ) <<  8 |
    ( value8 & INTEGER_CONSTANT( 0x0000000000ff0000 ) ) << 24 |
    ( value8 & INTEGER_CONSTANT( 0x000000000000ff00 ) ) << 40 |
    ( value8 & INTEGER_CONSTANT( 0x00000000000000ff ) ) << 56;
  *ivalue = result;

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: rotate8ByteArrayIfLittleEndian() - Rotate 8-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array    Array of 8-byte values to rotate.
         Integer count  Number of items in array.
OUTPUTS: void* array    Array of rotated values.
******************************************************************************/

void rotate8ByteArrayIfLittleEndian( void* array, Integer count ) {

#if IS_LITTLE_ENDIAN

  PRE02( array, count > 0 );
  assert_static( sizeof (Integer) == 8 );
  Integer* const array8 = array;
  Integer index = 0;

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const Integer value = array8[ index ];
    const Integer newValue =
    ( value & INTEGER_CONSTANT( 0xff00000000000000 ) ) >> 56 |
    ( value & INTEGER_CONSTANT( 0x00ff000000000000 ) ) >> 40 |
    ( value & INTEGER_CONSTANT( 0x0000ff0000000000 ) ) >> 24 |
    ( value & INTEGER_CONSTANT( 0x000000ff00000000 ) ) >>  8 |
    ( value & INTEGER_CONSTANT( 0x00000000ff000000 ) ) <<  8 |
    ( value & INTEGER_CONSTANT( 0x0000000000ff0000 ) ) << 24 |
    ( value & INTEGER_CONSTANT( 0x000000000000ff00 ) ) << 40 |
    ( value & INTEGER_CONSTANT( 0x00000000000000ff ) ) << 56;
    array8[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: expand32BitValues - Copy/expand 32-bit floating-point values to
         64-bit values.
INPUTS:  Real* array    Array of 32-bit values to expand in-place.
         Integer count  Number of values in array.
OUTPUTS: Real* array    Expanded array of 64-bit values.
******************************************************************************/

void expand32BitValues( Real* array, Integer count ) {
  PRE02( array, count > 0 );
  const float* const farray = (float*) array;
  const float* source = farray + count; /* 1 past last element. */
  Real* destination   = array  + count; /* 1 past last element. */

  do {
    *--destination = *--source; /* Expand 32-bits to 64-bits. */
  } while ( destination != array );
}



/******************************************************************************
PURPOSE: compress64BitValues - Copy/compress 64-bit floating-point values to
         32-bit values.
INPUTS:  Real* array    Array of 64-bit values to expand in-place.
         Integer count  Number of values in array.
OUTPUTS: Real* array    Compressed array of 32-bit values.
NOTES:   Values in the array are overwritten to the first half and will have
         lost precision and will be clamped to 32-bit range. UGLY.
******************************************************************************/

void compress64BitValues( Real* array, Integer count ) {
  PRE02( array, count > 0 );
  float* destination = (float*) array;

  while ( count-- ) {
    const float value = CLAMPED_TO_RANGE( *array, -FLT_MAX, FLT_MAX );
    *destination++ = value;
    ++array;
  }
}



/******************************************************************************
PURPOSE: windDirectionAndSpeed - Compute direction and speed from
         u/v-components.
INPUTS:  const Real windU      Eastward component (m/s).
         const Real windV      Northward component (m/s).
OUTPUTS: Real* windDirection   Bering angle the wind comes from.
         Real* windSpeed       Magnitude of wind vector (m/s).
******************************************************************************/

void windDirectionAndSpeed( const Real windU, const Real windV,
                            Real* windDirection, Real* windSpeed ) {
  PRE04( isFinite( windU ), isFinite( windV ), windDirection, windSpeed );
  const Real speed = hypot( windU, windV );
  const Real angleRadians = atan2( windV, windU );
  const Real angleDegrees0 = degrees( angleRadians );
  const Real angleDegrees =
    angleDegrees0 < 0.0 ? angleDegrees0 + 360.0 : angleDegrees0; 
  const Real direction0 = 270.0 - angleDegrees;
  const Real direction = direction0 < 0.0 ? direction0 + 360.0 : direction0;
  *windDirection = direction;
  *windSpeed = speed;
  POST02( IN_RANGE( *windDirection, 0.0, 360.0 ), *windSpeed >= 0.0 );
}



/******************************************************************************
PURPOSE: windUV - Compute wind u/v-components from direction and speed.
INPUTS:  const Real windDirection   Bering angle the wind comes from.
         const Real windSpeed       Magnitude of wind vector (m/s).
OUTPUTS: Real* windU                Eastward component (m/s).
         Real* windV                Northward component (m/s).
******************************************************************************/

void windUV( const Real windDirection, const Real windSpeed,
             Real* windU, Real* windV ) {
  POST04( IN_RANGE( windDirection, 0.0, 360.0 ), windSpeed >= 0.0,
          windU, windV );
  const Real direction = 270.0 - windDirection;
  const Real angleDegrees =
    direction < 0.0   ? direction + 360.0
    : direction > 360.0 ? direction - 360.0
    : direction;
  const Real angleRadians = radians( angleDegrees );
  const Real u = cos( angleRadians );
  const Real v = sin( angleRadians );
  *windU = windSpeed * u;
  *windV = windSpeed * v;
  POST02( isFinite( *windU ), isFinite( *windV ) );
}




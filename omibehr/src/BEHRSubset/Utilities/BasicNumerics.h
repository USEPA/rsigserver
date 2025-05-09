
#ifndef BASICNUMERICS_H
#define BASICNUMERICS_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: BasicNumerics.h - Defines typedefs for basic integer and floating-
         point types as exactly 64-bit values that facilitate passing
         between routines written in different languages (e.g., C/C++,
         FORTRAN, Java and Eiffel) and I/O transmission/storage without loss
         of precision.

NOTES:   This can only work on platforms that support 64-bit integers
         and 64-bit reals. Also, unsigned integers are not supported since
         the 64-bit integer is defined as a typedef and most other languages
         don't support unsigned integers. No assumption is made that the
         64-bit floating-point type is IEEE 754 compatible.
         For example, it may not support nan/inf/-inf and may also use an
         incompatible number of bits for the mantissa or characteristic as
         found on Crays, for example.

         Provisions are also made to alias the various string-to-integer
         and integer-to-string routines such as atoll() and strtoll(),
         printf()/scanf() format specification strings and variations of
         LONGLONG_MAX, etc., for example:

#include <stdio.h>         // For printf().

#include <Assertions.h>    // For CHECK().
#include <BasicNumerics.h> // For INTEGER_MAX, Integer, Real, atoI(), strtoR()

int main( void ) {

  const Integer maximumInteger      = INTEGER_MAX;
  const Integer minimumInteger      = INTEGER_MIN;
  const Real    maximumReal         = REAL_MAX;
  const Real    minimumReal         = -REAL_MAX;
  const Real    minimumPositiveReal = REAL_MIN;
  const Real    epsilonReal         = REAL_EPSILON;

  Integer i = INTEGER_MIN;
  Real    r = -REAL_MAX;

  printf( "Integer range: [ %"INTEGER_FORMAT", %"INTEGER_FORMAT" ]\n",
          minimumInteger, maximumInteger );

  printf( "Real range: [ %"REAL_E_FORMAT", %"REAL_E_FORMAT" ]\n",
          minimumReal, maximumReal );

  printf( "Real epsilon: %"REAL_E_FORMAT"\n", epsilonReal );

  printf( "i = %"INTEGER_FORMAT", and r = %"REAL_E_FORMAT"\n", i, r );
  i = atoI( "9223372036854775807" );
  r = atoR( "3.14159265358979323846" );
  printf( "i = %"INTEGER_FORMAT", and r = %"REAL_F_FORMAT"\n", i, r );
  i = strtoI( "9223372036854775807", 0, 10 );
  r = strtoR( "3", 0 );
  printf( "i = %"INTEGER_FORMAT", and r = %10.2"REAL_G_FORMAT"\n", i, r );
  sscanf( "-9223372036854775808 -3.14159265358979323846E+2",
          "%"INTEGER_FORMAT" %"REAL_F_FORMAT, &i, &r );

  // Check that the implementation is correct:

  CHECK( sizeof (Integer)    == 8           );
  CHECK( sizeof (Real)       == 8           );
  CHECK( maximumInteger      == INTEGER_MAX );
  CHECK( minimumInteger      == INTEGER_MIN );
  CHECK( maximumReal         == REAL_MAX    );
  CHECK( minimumReal         == -REAL_MAX   );
  CHECK( minimumPositiveReal == REAL_MIN    );
  CHECK( i == INTEGER_CONSTANT( -9223372036854775808 ) );
  CHECK( epsilonReal         == REAL_EPSILON );

  // Use macro to append 'LL' or 'L' as a suffix on all Integer constants:

  return i == INTEGER_CONSTANT( -9223372036854775808 );
}


HISTORY: 2001/04, Todd Plessel, EPA/LM, Created (based on NumericTypes.h).
STATUS:  reviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <limits.h> /* For DBL_MAX, DBL_MIN.             */
#include <float.h>  /* For DBL_MAX, DBL_MIN on DEC OSF1. */
#include <math.h>   /* For fabs().                       */
#include <stdlib.h> /* For atoll(), strtoll().           */

/*================================= MACROS ==================================*/

/* Macro to append 'L' or 'LL' to form appropriate constant */

#if _CRAY || __alpha
#define INTEGER_CONSTANT(x) x##L
#define ULONGLONG_CONSTANT(x) x##UL
#else
#define INTEGER_CONSTANT(x) x##LL
#define ULONGLONG_CONSTANT(x) x##ULL
#endif

/* Analogous to LONGLONG_MAX (2^31) and LONGLONG_MIN -(2^31)-1: */

#define INTEGER_MAX INTEGER_CONSTANT(9223372036854775807)
#define INTEGER_MIN \
  (INTEGER_CONSTANT(-9223372036854775807)-INTEGER_CONSTANT(1))

#ifndef LONGLONG_MAX
#define LONGLONG_MAX INTEGER_MAX
#endif

#ifndef LONGLONG_MIN
#define LONGLONG_MIN INTEGER_MIN
#endif

#ifndef ULONGLONG_MAX
#define ULONGLONG_MAX ULONGLONG_CONSTANT(18446744073709551615)
#endif


#define SIZET_MAX ULONG_MAX

#define MIN_OF_SIZET_OR_INTEGER_MAX \
  ((sizeof (size_t) < 8 && SIZET_MAX < INTEGER_MAX) ? SIZET_MAX : INTEGER_MAX)

/* Analogous to DBL_MAX and DBL_MIN. (REAL_MIN is smallest positive value): */

#define REAL_MAX DBL_MAX
#define REAL_MIN DBL_MIN
#define REAL_EPSILON DBL_EPSILON

/* Recommended tolerance for most floating-point operations. */

#define TOLERANCE 1.0e-6

/* Value clamped to range [low, high]. */

#define CLAMPED_TO_RANGE( value, low, high ) \
((value) < (low) ? (low) : (value) > (high) ? (high) : (value))

#define SIGN(x) ((x) < 0 ? -1 : 1)

/* For testing Integer increment: */

#define NEXT_STRICTLY_POSITIVE_INTEGER(i) ( (i) == INTEGER_MAX ? 1 : (i) + 1 )

/*
 * To facilitate correct and optimized I/O, determine if the platform is native
 * IEEE/MSB/XDR.
 */

#ifndef _IEEE
#if defined(__alpha) || defined(_CRAY) || \
    defined(__i386__) || defined(__i486__) || \
    defined(__i586__) || defined(__i686__) || \
    defined(__ia64__) || defined(__x86_64__)
#define _IEEE 0
#else
#define _IEEE 1
#endif
#endif

/*
 * Is the platform big-endian (MSB: most significant byte first) or
 * little-endian (LSB: least significant byte first)?
 */

#ifndef IS_LITTLE_ENDIAN
#if ( \
    ( defined(__BYTE_ORDER__) && \
      defined(__ORDER_LITTLE_ENDIAN__) && \
      __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ ) || \
    defined(__x86_64__) || \
    defined(__i386__)   || \
    defined(__i486__)   || \
    defined(__i586__)   || \
    defined(__i686__)   || \
    defined(__alpha)    || \
    defined(__ia64__)   || \
    defined(__ARMEL__)  || \
    defined(__MIPSEL)   || \
    defined(_WIN32)     || \
    defined(_WIN64)     || \
    defined(_M_IX86)    || \
    defined(_M_AMD64)      \
  )
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif
#endif


/* Is the platform native XDR? */

#if _IEEE == 0 || IS_LITTLE_ENDIAN == 1
#define IS_NATIVE_XDR 0
#else
#define IS_NATIVE_XDR 1
#endif

/* These are defined by the XDR standard: */

#define SIZEOF_XDR_int    4
#define SIZEOF_XDR_long   4
#define SIZEOF_XDR_float  4
#define SIZEOF_XDR_double 8

/* See (at compile time) if the native size matches the size of the XDR type: */
#define USES_NATIVE_XDR(type_) \
  (IS_NATIVE_XDR && sizeof (type_) == SIZEOF_XDR_##type_)


/*
 * Facilities for string-to-integer conversions:
 *   atoI(), strtoI(), atoR(), strtoR()
 * and printf()/scanf() format specifications:
 *   INTEGER_FORMAT, REAL_FORMAT
 */

/*
 * Select the appropriate strtoll() implementation and 64-bit integer format
 * specification for printf/scanf:
 */

#if ( __osf__ || _CRAY )
/*
 * strtoll() is missing and "%lld" is invalid,
 * but strtol() is 64-bits as is "%ld" so just use them instead.
 */
#define strtoll strtol
#define INTEGER_FORMAT "ld"
#elif _WIN32
#define INTEGER_FORMAT "I64d"
#elif __OPENNT
/*
 * strtoll() is missing,
 * but strtoq() exists in libc.a and works for 64-bit integers (typedef quad_t)
 * however, it is not declared in a header so declare it & alias it to strtoll.
 * scanf( "%lld", &a_long_long ) is broken (as is "%qd") for 64-bit integers.
 * Therefore, instead of scanf/sscanf/fscanf() of 64-bit integers,
 * use (aliased) strtoll(), reading a string first if necessary.
 * (Note: printf( "%lld", a_long_long ) works.)
 */
extern long long strtoq( const char* s, char** unused1, int unused2 );
#define strtoll strtoq
#define SCANF_IS_BROKEN 1
#elif __APPLE__
/*
 * "%lld" does not work for scanf of 64-bit integers, but "%qd" does so use it.
 * strtoll() is missing but strtoq() works so alias it.
 */
#define INTEGER_FORMAT "qd"
#define strtoll strtoq
#else
#define INTEGER_FORMAT "lld"
#endif


/* Define atoll() in terms of strtoll() for all platforms: */

#define atoll(s) strtoll((s),0,10)

/*
 * Now define the custom versions - atoI, atoR(), strtoI(), strtoR() -
 * used for 64-bit Integer and 64-bit Real types:
 */

#define atoI atoll
#define strtoI strtoll
#define atoR atof
#define strtoR strtod
#define REAL_F_FORMAT "lf"
#define REAL_E_FORMAT "le"
#define REAL_G_FORMAT "lg"

/*================================== TYPES ==================================*/

/*
 * Finally, here are the definitions of the 64-bit integer and
 * 64-bit floating-point types to be used in place of the various
 * ambiguously-sized C basic types. (Using the typedefs allows for easy
 * redefinition on new platforms that do support 64-bit integers and
 * floating-point values, but use different C types to do so.)
 *
 * Replace all usage (in public interfaces at least) of:
 *   signed/unsigned short, int, long, long long, and size_t with Integer and
 *   float, double and long double with Real.
 */

typedef long long Integer; /* Exactly 64-bits on all supported platforms. */
typedef double Real;       /* Exactly 64-bits on all supported platforms. */

/*================================ FUNCTIONS ================================*/

extern Integer isSignedChar(       Integer value );
extern Integer isUnsignedChar(     Integer value );
extern Integer isSignedShort(      Integer value );
extern Integer isUnsignedShort(    Integer value );
extern Integer isSignedInt(        Integer value );
extern Integer isUnsignedInt(      Integer value );
extern Integer isSignedLong(       Integer value );
extern Integer isUnsignedLong(     Integer value );
extern Integer isSizet(            Integer value );
extern Integer isSignedLongLong(   Integer value );
extern Integer isUnsignedLongLong( Integer value );
extern Integer absI( Integer x );

extern Integer toInteger( const char* string, Integer lower, Integer upper,
                          Integer* ok );

extern Real toReal( const char* string, Real lower, Real upper, Integer* ok );

extern Integer isInfinity( Real x );
extern Integer isMinusInfinity( Real x );
extern Integer isFinite( Real x );
extern Integer isNan(  Real x );
extern Real filterNan( Real x );
extern Integer withinTolerance( Real x, Real y, Real tolerance );
extern Integer aboutEqual( Real x, Real y );

extern Real radians( Real theDegrees );
extern Real degrees( Real theRadians );

extern Real safeSum( Real a1, Real a2 );
extern Real safeSum3( Real a1, Real a2, Real a3 );
extern Real safeSum4( Real a1, Real a2, Real a3, Real a4 );
extern Real safeSum5( Real a1, Real a2, Real a3, Real a4, Real a5 );
extern Real safeSum6( Real a1, Real a2, Real a3, Real a4, Real a5, Real a6 );
extern Real safeSum7( Real a1, Real a2, Real a3, Real a4, Real a5, Real a6,
                      Real a7 );
extern Real safeSum8( Real a1, Real a2, Real a3, Real a4, Real a5, Real a6,
                      Real a7, Real a8 );

extern Real safeDifference( Real left, Real right );

extern Real safeProduct( Real a1, Real a2 );
extern Real safeProduct3( Real a1, Real a2, Real a3 );
extern Real safeProduct4( Real a1, Real a2, Real a3, Real a4 );
extern Real safeProduct5( Real a1, Real a2, Real a3, Real a4, Real a5 );

extern Real safeProduct6( Real a1, Real a2, Real a3, Real a4, Real a5,
                          Real a6 );

extern Real safeProduct7( Real a1, Real a2, Real a3, Real a4, Real a5,
                          Real a6, Real a7 );

extern Real safeProduct8( Real a1, Real a2, Real a3, Real a4, Real a5,
                          Real a6, Real a7, Real a8 );

extern Real safeQuotient( Real numerator, Real denominator );

extern Real kahanSum( const Real* items, Integer count );
extern void filterNans( Real* items, Integer count );
extern Integer isNanFree( const Real* items, Integer count );
extern Integer allFinite( const Real* items, Integer count );
extern Integer allZero( const Real* items, Integer count );
extern Integer allZeroI( const Integer* items, Integer count );
extern Integer allNull( const void* items[], Integer count );
extern Integer nonNull( const void* items[], Integer count );
extern Real maximumItem( const Real* items, Integer count );
extern Real minimumItem( const Real* items, Integer count );
extern Integer maximumItemI( const Integer* items, Integer count );
extern Integer minimumItemI( const Integer* items, Integer count );
extern Integer decreasing( const Real* items, Integer count );
extern Integer increasing( const Real* items, Integer count );
extern Integer sumI( const Integer* items, Integer count );

extern void fillRandom( Real array[], Integer count );
extern void fillRandomI( Integer array[], Integer count,
                         Integer low, Integer high );
extern void shuffle( Real array[], Integer count );
extern void shuffleI( Integer array[], Integer count );
extern Integer isSorted( const Real array[], Integer count );
extern Integer isSortedI( const Integer array[], Integer count );
extern void shellsort( Real array[], Integer count );
extern void shellsortI( Integer array[], Integer count );
extern Integer randomInteger( Integer low, Integer high );
extern void rotate8ByteArrayIfLittleEndian( void* array, Integer count );
extern void rotate4ByteWordIfLittleEndian( void* value );
extern void rotate4ByteArrayIfLittleEndian( void* array, Integer count );
extern void rotate2ByteWordIfLittleEndian( void* value );
extern void rotate2ByteArrayIfLittleEndian( void* array, Integer count );
extern void compress64BitValues( Real* array, Integer count );

extern Integer isValidArgs( Integer argc, const char* argv[] );

#ifdef __cplusplus
}
#endif

#endif /* BASICNUMERICS_H */


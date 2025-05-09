
#ifndef UTILITIES_H
#define UTILITIES_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Utilities.h - Declare some general-purpose reusable routines.

HISTORY: 2017-01-02 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================= MACROS ==================================*/

#define MISSING_VALUE (-9999.0)

#define OR2(a,b) ((a)||(b))
#define OR3(a,b,c) ((a)||(b)||(c))
#define OR4(a,b,c,d) ((a)||(b)||(c)||(d))
#define OR5(a,b,c,d,e) ((a)||(b)||(c)||(d)||(e))
  
#define AND2(a,b) ((a)&&(b))
#define AND3(a,b,c) ((a)&&(b)&&(c))
#define AND4(a,b,c,d) ((a)&&(b)&&(c)&&(d))
#define AND5(a,b,c,d,e) ((a)&&(b)&&(c)&&(d)&&(e))
#define AND6(a,b,c,d,e,f) ((a)&&(b)&&(c)&&(d)&&(e)&&(f))
  
#define NON_ZERO7(a,b,c,d,e,f,g) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g))
#define IS_ZERO2(a,b) ((a)==0&&(b)==0)
#define IS_ZERO7(a,b,c,d,e,f,g) \
((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0)

#define IMPLIES(p,c) (!(p)||(c))
#define IMPLIES_ELSE(p,c1,c2) (((p)&&(c1))||((!(p))&&(c2)))

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

#define CLAMPED_TO_RANGE( value, low, high ) \
((value) < (low) ? (low) : (value) > (high) ? (high) : (value))
  
#define IS_BOOL(x) ((x)==0||(x)==1)
#define SIGN(x) ((x)<0?-1:1)

#define IN3(x,a,b) ((x)==(a)||(x)==(b))
#define IN4(x,a,b,c) ((x)==(a)||(x)==(b)||(x)==(c))
#define IN5(x,a,b,c,d) ((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d))

#define IN7(x,a,b,c,d,e,f) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f))

#define IN9(x,a,b,c,d,e,f,g,h) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||\
(x)==(h))

#define IN10(x,a,b,c,d,e,f,g,h,i) \
((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f)||(x)==(g)||\
(x)==(h)||(x)==(i))

#ifdef DEBUGGING
#define DEBUG( s ) s
#else
#define DEBUG(unused)
#endif

/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };

typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE,LATITUDE ][ MINIMUM,MAXIMUM ]*/

/*================================ FUNCTIONS ================================*/

extern int isValidLongitude( const double longitude );

extern int isValidLatitude( const double latitude );

extern int isValidElevation( const double elevation );

extern int clampInvalidCoordinates( const size_t points,
                                    double longitudes[],
                                    double latitudes[] );

extern void computeBounds( const size_t count,
                           const double longitudes[],
                           const double latitudes[],
                           Bounds bounds );

extern int compactPointsInSubset( const Bounds domain,
                                  const double minimumElevation,
                                  const double maximumElevation,
                                  const size_t points,
                                  const size_t levels,
                                  double timestamps[],
                                  double longitudes[],
                                  double latitudes[],
                                  double elevations[],
                                  double values[],
                                  double thicknesses[],
                                  size_t* subsetPoints,
                                  size_t* subsetLevels );

extern int isLeapYear( const int yyyy );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int isValidYYYYMMDDHHMM( const long long yyyymmddhhmm );

extern int isValidYYYYDDDHHMM( const long long yyyymmddhhmm );

extern long long convertTimestamp( const long long yyyymmddhhmm );

extern long long offsetTimestamp( const long long yyyydddhhmm,
                                  const long long hours );

extern int isValidBounds( const Bounds bounds );

extern int boundsOverlap( const Bounds a, const Bounds b );

extern void rotate8ByteArrayIfLittleEndian( const size_t count, void* array );

extern void reverseLevels( const size_t points, const size_t levels,
                           double array[] );

extern void expandInt8( const size_t count, double array[] );

extern void expandUint8( const size_t count, double array[] );

extern void expandInt16( const size_t count, double array[] );

extern void expandUint16( const size_t count, double array[] );

extern void expandInt32( const size_t count, double array[] );

extern void expandUint32( const size_t count, double array[] );

extern void expandReals( const size_t count, double array[] );

extern void scaleValues( const double factor, const size_t count,
                         double array[] );

extern void offsetValues( const double offset, const size_t count,
                          double array[] );

extern void copyVectorComponent( const size_t count, const size_t components,
                                 const size_t component,
                                 const double input[], double output[] );

extern void copyMaximumComponent( const size_t count, const size_t components,
                                  const double input[], double output[] );

extern void copyMeanComponents( const size_t count, const size_t components,
                                const double input[], double output[] );

extern size_t fileSize( const char* name );

extern char* readFile( const char* name, size_t* length );

extern void controlMToSpace( char* string );

extern void spacesToUnderscores( char* string );

extern size_t linesInString( const char* string );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */




#ifndef UTILITIES_H
#define UTILITIES_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Utilities.h - Declare some general-purpose reusable routines.

HISTORY: 2017-10-14 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================= MACROS ==================================*/

#define MISSING_VALUE (-9.999e36)

#ifndef NDEBUG
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };
enum { COLUMN, ROW };
enum { FIRST, LAST };

typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE,LATITUDE ][ MINIMUM,MAXIMUM ]*/

/*================================ FUNCTIONS ================================*/

extern int clampInvalidCoordinates( const size_t points,
                                    double longitudes[], double latitudes[] );

extern size_t pointsInDomain( const Bounds domain,
                              const size_t points,
                              const double longitudes[],
                              const double latitudes[],
                              const double values[],
                              unsigned char mask[] );
  
extern void computeCorners( const size_t rows,
                            const size_t columns,
                            const double longitudes[],
                            const double latitudes[],
                            double longitudesSW[],
                            double longitudesSE[],
                            double longitudesNW[],
                            double longitudesNE[],
                            double latitudesSW[],
                            double latitudesSE[],
                            double latitudesNW[],
                            double latitudesNE[] );

extern void clampLongitudes( const double longitude,
                             double* longitude1,
                             double* longitude2,
                             double* longitude3,
                             double* longitude4 );

extern void replicateRows( const size_t columns,
                           const size_t rows,
                           double data[] );

extern void replicateColumns( const size_t rows,
                              const size_t columns,
                              double data[] );

extern void transpose( const size_t rows, const size_t columns,
                       double data[], double temp[] );

extern int isLeapYear( const int yyyy );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int isValidYYYYMMDDHHMM( const long long yyyymmddhhmm );

extern int isValidYYYYDDDHHMM( const long long yyyymmddhhmm );

extern long long convertTimestamp( const long long yyyymmddhhmm );

extern int incrementHours( const int yyyymmddhh, const int hours );

extern int hoursUntil( const int yyyymmddhh1, const int yyyymmddhh2 );

extern int isValidBounds( const Bounds bounds );

extern int boundsOverlap( const Bounds a, const Bounds b );

extern void fillArray( const double value, const size_t count, double array[] );

extern int writeArray( const double array[], const size_t count,
                       const int times );

extern int fillAndWriteArray( const double value, const size_t count,
                              double array[], const int times );

extern void rotate8ByteArrayIfLittleEndian( void* array, const size_t count);

extern size_t fileSize( const char* name );

extern char* readFile( const char* name, size_t* length );

extern void controlMToSpace( char* string );

extern size_t linesInString( const char* string );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */



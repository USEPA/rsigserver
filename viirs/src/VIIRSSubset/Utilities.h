
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

#define MISSING_VALUE (-9999.0)

/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };

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

extern int isLeapYear( const int yyyy );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int isValidYYYYMMDDHHMM( const long long yyyymmddhhmm );

extern int isValidYYYYDDDHHMM( const long long yyyymmddhhmm );

extern long long convertTimestamp( const long long yyyymmddhhmm );

extern int isValidBounds( const Bounds bounds );

extern int boundsOverlap( const Bounds a, const Bounds b );

extern void expand32BitReals( const size_t count, double* array );

extern void rotate8ByteArrayIfLittleEndian( void* array, const size_t count);

extern size_t fileSize( const char* name );

extern char* readFile( const char* name );

extern size_t linesInString( const char* string );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */




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
#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { COLUMN, ROW };
enum { MINIMUM, MAXIMUM };

typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE,LATITUDE ][ MINIMUM,MAXIMUM ]*/

/*================================ FUNCTIONS ================================*/

extern size_t pointsInDomain( const Bounds domain,
                              const size_t points,
                              const double longitudes[],
                              const double latitudes[],
                              unsigned char mask[] );
  
extern int isLeapYear( const int yyyy );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int incrementHours( const int yyyymmddhh, const int hours );

extern int isValidBounds( const Bounds bounds );

extern int boundsOverlap( const Bounds a, const Bounds b );

extern int subsetIndicesByBounds( const Bounds bounds,
                                  const size_t rows,
                                  const size_t columns,
                                  const double longitudes[],
                                  const double latitudes[],
                                  size_t* firstRow,
                                  size_t* lastRow,
                                  size_t* firstColumn,
                                  size_t* lastColumn );

extern void doublesToFloats( double* array, const size_t count);

extern void rotate4ByteArrayIfLittleEndian( void* array, const size_t count);

extern void rotate8ByteArrayIfLittleEndian( void* array, const size_t count);

extern void fillArray( const double value, const size_t count, double array[] );

extern size_t fileSize( const char* name );

extern char* readFile( const char* name, size_t* length );

extern void controlMToSpace( char* string );

extern size_t linesInString( const char* string );


#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */



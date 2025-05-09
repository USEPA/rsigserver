
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

#include <stdlib.h> /* For size_t, malloc(), free(). */
#include <stdio.h>  /* For FILE. */

/*================================== MACROS =================================*/

extern void* newMemory( size_t, size_t );
#define NEW( type, count ) ((type*) newMemory( (count), sizeof (type) ))
#define FREE( p ) ( ( (p) ? free(p) : (void) 0 ), (p) = 0 )

#define MIN(a,b) ((a)<(b)?(a):(b))

#define COUNT_IN_RANGE( first, last ) ( (last) - (first) + 1 )

/* Value clamped to range [low, high]. */

#define CLAMPED_TO_RANGE( value, low, high ) \
  ((value) < (low) ? (low) : (value) > (high) ? (high) : (value))

/* From M3IO Library: */

#define BADVAL3 (-9.999E36)
#define AMISS3 (-9.000E36)
#define IS_VALID_VALUE( value ) ( (value) > AMISS3 && (value) < -AMISS3 )
  
/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };
enum { COLUMN, ROW, LAYER, TIME, DIMENSIONS };
enum { TIMESTAMP_LENGTH = 24 }; /* YYYY-MM-DDTHH:MM:SS-ZZZZ */

typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE,LATITUDE ][ MINIMUM,MAXIMUM ]*/

/* Command-line options: */

enum {
  NO_TYPE, FILE_TYPE, DIRECTORY_TYPE, STRING_TYPE, ENUM_TYPE,
  INT_TYPE, INTEGER64_TYPE, REAL64_TYPE, YYYYMMDDHH_TYPE, BOUNDS_TYPE
};

typedef struct {
  const char* name;     /* E.g., -tmpdir. */
  const int   required; /* Is option mandatory? */
  const int   type;     /* Type of values. E.g., DIRECTORY_TYPE. */
  const int   count;    /* Number of values for option. E.g., 0, -2, 4. */
                        /* -n means any count in the range [1, n]. */
  const void* range;    /* Minimum, maximum valid values to accept. */
  const char* valids;   /* String of space-delimied valid values to accept. */
  int         parsed;   /* Number of values parsed including name. */
  void*       values;   /* Array of count values to store. */
} Option;

/*================================ FUNCTIONS ================================*/

extern void* newMemory( size_t count, size_t sizeEach ); /* Called by NEW() */

extern int inRange( const size_t count, const double values[],
                    const double minimum, const double maximum );

extern int inRangeF( const size_t count, const float values[],
                     const double minimum, const double maximum );

extern int valuesIncrease( const size_t count, const double values[] );

extern int valuesDecrease( const size_t count, const float values[] );

/* Parsing: */

extern void underscoreToSpace( char* string );

extern int indexOfWord( const char* word, const char* words );

extern char* paddedString( const char* const input, const size_t length,
                           char output[] );

extern int parseOptions( const int argc, char* argv[], const int count,
                         Option options[] );

extern int parseOption( const int argc, char* argv[],
                        int* const arg, Option* const option );

extern int parseOptionValue( const char* argument, const int valueIndex,
                             Option* const option );

extern int isDouble( const char* const string );

extern void eraseBackslashR( char* string );

/* Date-time routines: */

extern int isLeapYear( const int yyyy );

extern int daysInMonth( const int yyyy, const int mm );

extern int isValidYYYYDDD( const int yyyyddd );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int isValidHHMMSS( const int hhmmss );

extern int toYYYYMMDD( const int yyyyddd );

extern int toYYYYDDD( const int yyyymmdd );

extern void nowUTC( int* const yyyy,
                    int* const ddd,
                    int* const hh,
                    int* const mm,
                    int* const ss );

extern int incrementHours( const int yyyymmddhh, const int hours );

extern int hoursInRange( const int yyyymmddhh1, const int yyyymmddhh2 );

extern int timestepsUntil( const int yyyymmddhh1,
                           const int yyyymmddhh2,
                           const int hours );

/* I/O: */

extern void rotate4ByteArrayIfLittleEndian( const size_t count, void* array );

extern int writeFloats( const size_t count, float array[], FILE* output );

extern int readFloats(  const size_t count, float array[], FILE* output );

extern int isDirectory( const char* name );

extern int fileDateUTC( const char* name );

extern size_t fileSize( const char* name );

extern int streamFile( const char* const name );

extern int printWorkingDirectory( void );
  
extern int isNetCDFFile( const char* const name );

extern int printDirectoryListing( const char* const name );

/* Geometry: */

extern int isValidBounds( const Bounds bounds );

extern int boundsOverlap( const Bounds a, const Bounds b );

extern int boundsSubsumes( const Bounds a, const Bounds b );

extern int clipPolygon( const int discardDegenerates,
                        const double clipXMin,
                        const double clipYMin,
                        const double clipXMax,
                        const double clipYMax,
                        const int count,
                        const double x[],
                        const double y[],
                        double cx[],
                        double cy[] );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */



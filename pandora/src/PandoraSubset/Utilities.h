
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

/*================================== MACROS =================================*/

#define IN3(x,a,b ) ((x)==(a)||(x)==(b))
#define IN4(x,a,b,c ) ((x)==(a)||(x)==(b)||(x)==(c))
#define IN5(x,a,b,c,d ) ((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d))

/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };

typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE,LATITUDE ][ MINIMUM,MAXIMUM ]*/

extern void* newMemory( size_t count, size_t  sizeEach );
#define NEW( type, count ) ((type*) newMemory( (count), sizeof (type) ))
#define FREE( p ) ( ( (p) ? free(p) : (void) 0 ), (p) = 0 )

/* Command-line options: */

enum {
  NO_TYPE, FILE_TYPE, DIRECTORY_TYPE, STRING_TYPE, ENUM_TYPE,
  INT_TYPE, INTEGER64_TYPE, REAL64_TYPE, YYYYMMDDHHMMSS_TYPE, BOUNDS_TYPE
};

typedef struct {
  const char* name;     /* E.g., -tmpdir. */
  const int   required; /* Is option mandatory? */
  const int   type;     /* Type of values. E.g., DIRECTORY_TYPE. */
  const int   count;    /* Number of values for option. E.g., 0, 1, 4. */
  const void* range;    /* Minimum, maximum valid values to accept. */
  const char* valids;   /* String of space-delimied valid values to accept. */
  int         parsed;   /* Was value parsed? */
  void*       values;   /* Array of count values to store. */
} Option;

/*================================ FUNCTIONS ================================*/

extern int parseOptions( const int argc, char* argv[], const int count,
                         Option options[] );

extern int parseOption( const int argc, char* argv[],
                        int* const arg, Option* const option );

extern int parseOptionValue( const char* argument, const int valueIndex,
                             Option* const option );

extern int isLeapYear( const int yyyy );

extern int isValidYYYYMMDDHHMMSS( const long long yyyymmddhhmmss );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int isValidBounds( const Bounds bounds );

extern void rotate8ByteArrayIfLittleEndian( void* array, const size_t count);

extern int isDirectory( const char* name );

extern size_t fileSize( const char* name );

extern int readFile( const char* name, size_t* length, char** content );

extern int streamBytes( FILE* const file, size_t bytes );

extern int nextLine( char** line, char** end );

extern int indexOfWord( const char* word, const char* words );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */



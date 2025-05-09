
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

/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };

typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE,LATITUDE ][ MINIMUM,MAXIMUM ]*/

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

extern int parseOptions( int argc, char* argv[], int count, Option options[] );

extern int isLeapYear( const int yyyy );

extern int isValidYYYYMMDDHHMMSS( const long long yyyymmddhhmmss );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int isValidBounds( const Bounds bounds );

extern void rotate8ByteArrayIfLittleEndian( void* array, const size_t count);

extern int isDirectory( const char* name );

extern size_t fileSize( const char* name );

extern int readFile( const char* name, size_t* length, size_t* lines,
                     char** content );

extern size_t controlMToSpace( char* string );

extern int indexOfWord( const char* word, const char* words );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */



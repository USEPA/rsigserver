
#ifndef UTILITIES_H
#define UTILITIES_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Utilities.h - Group related header files for client convenience.
NOTES:   See included header files and their documentation for details.
HISTORY: 2004/07 Todd Plessel EPA/LM Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Assertions.h>
#include <BasicNumerics.h>
#include <DateTime.h>
#include <Failure.h>
#include <Memory.h>
#include <Stream.h>
#include <Projector.h>
#include <Lambert.h>
#include <Mercator.h>
#include <Stereographic.h>
#include <Grid.h>
#include <VoidList.h>
#include <RegridQuadrilaterals.h>

static const Real invalid = -9999.0; /* Invalid data value. */

enum { X, Y, Z };
enum { COLUMN, ROW, LAYER };
enum { LONGITUDE, LATITUDE, ELEVATION };
enum { MINIMUM, MAXIMUM };
enum { FIRST, LAST };

typedef Real Bounds[ 2 ][ 2 ]; /* [ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ]*/

enum { FILE_NAME_LENGTH = 256 };
typedef char FileName[ FILE_NAME_LENGTH ]; /* pathed file name. */


extern Integer isValidBounds( const Bounds bounds );

extern Integer pointInDomain( Real longitude, Real latitude,
                              const Bounds domain );

extern Integer overlap( const Bounds a, const Bounds b );

extern Integer isValidArgs( Integer argc, const char* argv[] );

extern void checkForTest( int* argc, char* argv[] );

extern const char* parseArgument2( Integer argc, char* argv[],
                                  const char* option, Integer* arg );

extern Integer parseBounds( Integer argc, char* argv[], Integer* arg,
                            Bounds bounds );

extern Integer parseTimestampAndHours( Integer argc, char* argv[],
                                       Integer* arg,
                                       Integer* timestamp, Integer* hours );

extern Integer indexOfString( const char* string,
                              const char* const strings[], Integer count );

extern void lowercase( char string[] );

extern void uppercase( char string[] );

extern Integer fileExists( const char* name );

extern Integer fileSize( const char* name );

extern char* readFile( const char* name, Integer* length );

extern void controlMToSpace( char* string );

extern void trimTrailingWhitespace( size_t size, char string[] );

extern Integer linesInString( const char* string );

extern const char* skipLines( const char* string, Integer lines );

extern char* startOfLineWithWord( const char* string, const char* word );

extern char* findLine( char* string, const char* tag, char** end );

extern Integer binIndex( Integer index, Integer bins, const Integer counts[] );

extern void appendNote( RegriddedNote regriddedNote, const Note note );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */


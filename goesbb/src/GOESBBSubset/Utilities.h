
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

/*================================ FUNCTIONS ================================*/

extern int isLeapYear( const int yyyy );

extern int isValidYYYYMMDDHH( const int yyyymmddhh );

extern int incrementHours( const int yyyymmddhh, const int hours );

extern int isValidBounds( const Bounds bounds );

extern void expand32BitReals( const size_t count, double* array );

extern void rotate8ByteArrayIfLittleEndian( void* array, const size_t count);

extern size_t fileSize( const char* name );

extern char* readFile( const char* name, size_t* length );

extern void controlMToSpace( char* string );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */



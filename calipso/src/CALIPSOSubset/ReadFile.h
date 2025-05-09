
#ifndef READFILE_H
#define READFILE_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadFile.h - Declare simple to use wrapper routines to read data from
         CALIPSO HDF files.

NOTES:   Uses HDF libraries and libs they depend on (z, jpeg).

HISTORY: 2017-01-02 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include "Utilities.h" /* For Bounds. */

/*================================ FUNCTIONS ================================*/

extern int openFile( const char* const fileName );

extern void closeFile( const int file );

extern int readFileBounds( const int file, Bounds bounds );

extern int fileVariableExists( const int file, const char* const variable );

extern int readVariableDimensions( const int file, const char* const variable,
                                   int* const rank, int dimensions[ 32 ] );

extern int readFileData( const int file, const char* const variable,
                         const int rank, const int dimensions[],
                         char units[ 80 ], double data[] );

extern int readFileVData( const int file, const char* const variable,
                          const int count, double data[] );

#ifdef __cplusplus
}
#endif

#endif /* READFILE_H */




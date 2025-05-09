
#ifndef READFILE_H
#define READFILE_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare simple to use wrapper routines to read data from
         ceilometer HDF5 files.

NOTES:   Uses NetCDF4 library and libs it depends on (HDF5, curl, z, dl).

HISTORY: 2021-07-12 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================ FUNCTIONS ================================*/


extern int openFile( const char* const fileName );

extern void closeFile( const int file );

extern int readVariableId( const int file, const char* const variable );

extern int readFileAttribute( const int file, const char* const attribute,
                              double* value );

extern int readVariableDimensions( const int file, const int id,
                                   size_t* dimension0, size_t* dimension1 );

extern int readFileData( const int file, const int id,
                         const size_t dimension0, const size_t dimension1,
                         double data[] );

#ifdef __cplusplus
}
#endif

#endif /* READFILE_H */




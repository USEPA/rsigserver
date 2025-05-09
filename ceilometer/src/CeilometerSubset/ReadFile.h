
#ifndef READFILE_H
#define READFILE_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare simple to use wrapper routines to read data from
         ceilometer HDF5 files.

NOTES:   Uses HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-11-29 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================ FUNCTIONS ================================*/


extern int openFile( const char* const fileName );

extern void closeFile( const int file );

extern int openDataset( const int file, const char* const dataset );

extern void closeDataset( const int dataset );

extern int readFileAttribute( const int file, const char* const attribute,
                              double* value );

extern int readDatasetDimensions( const int dataset,
                                  size_t* dimension0, size_t* dimension1 );

extern int readFileData( const int dataset,
                         const size_t dimension0,
                         const size_t dimension1,
                         double data[] );

extern int readFileDataIntegers( const int dataset,
                                 const size_t dimension0,
                                 const size_t dimension1,
                                 long long data[] );

#ifdef __cplusplus
}
#endif

#endif /* READFILE_H */




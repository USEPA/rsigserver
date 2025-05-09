
#ifndef READFILE_H
#define READFILE_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare simple to use wrapper routines to read data from
         OMI-BEHR HDF5 files.

NOTES:   Uses HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-12-06 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================ CONSTANTS ================================*/

#define MISSING_VALUE (-9999.0)

/*================================ FUNCTIONS ================================*/


extern int openFile( const char* const fileName );

extern void closeFile( const int file );

extern int swathsInFile( const int file );

extern int readDimensions( const int file, const int swath,
                           size_t* rows, size_t* columns );

extern int readDataset( const int file,
                        const int swath,
                        const size_t rows,
                        const size_t columns,
                        const char* const variable,
                        char units[ 80],
                        double data[] );

#ifdef __cplusplus
}
#endif

#endif /* READFILE_H */




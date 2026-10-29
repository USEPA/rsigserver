
#ifndef READDATA_H
#define READDATA_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare simple to use wrapper routines to read data from
         VIIRS NetCDF4 files.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-10-14 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

#include "Utilities.h" /* For Bounds. */

/*================================ FUNCTIONS ================================*/

extern int openFile( const char* const fileName );

extern void closeFile( const int file );

extern int readFileBounds( const int file, Bounds bounds );

extern int readFileDimensions( const int file,
                               size_t* const rows, size_t* const columns );

extern int readFileData( const int file, const char* const variable,
                         const size_t rows, const size_t columns,
                         const int qcMinimum,
                         const int groundPixelMinimum,
                         const int groundPixelMaximum,
                         const double maximumCloudFraction,
                         const int allowNegativeCounts,
                         char units[ 80 ], double data[] );

#ifdef __cplusplus
}
#endif

#endif /* READDATA_H */




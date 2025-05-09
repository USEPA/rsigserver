
#ifndef READDATA_H
#define READDATA_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare simple to use wrapper routines to read data from
         NetCDF4 files.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2019-06-25 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================ FUNCTIONS ================================*/

extern int openFile( const char* const fileName );

extern void closeFile( const int file );

extern int readFileDimensions( const int file,
                               const char* const product,
                               const char* const variable,
                               size_t* const rows,
                               size_t* const columns );

extern size_t readFileData( const int file,
                            const char* const product,
                            const char* const variable,
                            const size_t rows,
                            const size_t columns,
                            const size_t gridSubsetIndices[2][2],
                            const int qcMinimum,
                            const double maximumCloudFraction,
                            const double maximumSolarZenithAngle,
                            const int allowNegativeCounts,
                            char units[ 80 ],
                            double data[],
                            double scratch[] );

#ifdef __cplusplus
}
#endif

#endif /* READDATA_H */




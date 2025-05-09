
#ifndef READFILE_H
#define READFILE_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare simple to use wrapper routines to read data from
         OMI-AURA HDF5-EOS files.

NOTES:   Uses HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-12-06 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================ CONSTANTS ================================*/

#define MISSING_VALUE (-9.999e36)

/*================================ FUNCTIONS ================================*/


extern int openFile( const char* const fileName );

extern void closeFile( const int file );

extern int readDimensions( const int file, const char* const product,
                           size_t* rows, size_t* columns );

extern size_t readDataset( const int file,
                           const size_t rows,
                           const size_t columns,
                           const char* const product,
                           const char* const variable,
                           const double maximumCloudFraction,
                           const double maximumSolarZenithAngle,
                           const int allowNegativeCounts,
                           char units[ 80],
                           double data[],
                           double temp[] );

#ifdef __cplusplus
}
#endif

#endif /* READFILE_H */




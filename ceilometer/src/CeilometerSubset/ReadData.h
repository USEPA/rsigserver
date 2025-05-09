
#ifndef READDATA_H
#define READDATA_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare routine for reading ceilometer data and
         subsetting to domain.

NOTES:   Uses HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-11-29 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include "Utilities.h" /* For Bounds. */

/*================================ FUNCTIONS ================================*/

extern size_t readSubsetCeilometerData( const char* const fileName,
                                        const char* const variable,
                                        const Bounds domain,
                                        const Integer firstTimestamp,
                                        const Integer lastTimestamp,
                                        char units[ 80 ],
                                        double* longitude,
                                        double* latitude,
                                        double* elevation,
                                        double** data,
                                        double** elevations,
                                        double** timestamps );

#ifdef __cplusplus
}
#endif

#endif /* READDATA_H */




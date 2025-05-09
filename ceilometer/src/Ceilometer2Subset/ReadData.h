
#ifndef READDATA_H
#define READDATA_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare routine for reading ceilometer data and
         subsetting to domain.

NOTES:   Uses NetCDF4 library and libs it depends on (HDF5, curl, z, dl).

HISTORY: 2022-07-12 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include "Utilities.h" /* For Bounds. */

/*================================ CONSTANTS ================================*/

#define MISSING_VALUE (-9999.0)

#define MINIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE 0.0
#define MAXIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE 540000.0

/* Meters above mean sea level: */

#define MINIMUM_SURFACE_ELEVATION (-500.0)
#define MAXIMUM_SURFACE_ELEVATION 1e4
#define MAXIMUM_ELEVATION 1e5

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




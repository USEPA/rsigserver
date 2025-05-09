
#ifndef READDATA_H
#define READDATA_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare routine to read and filter/process CALIPSO data.

NOTES:   Uses routines from ReadFile.c and Utilities.c.

HISTORY: 2017-01-02 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

#include "Utilities.h" /* For Bounds. */

/*================================== TYPES ==================================*/

/* Type of CALIPSO file: */

enum {
  CALIPSO_L1, CALIPSO_L2_05KMAPRO, CALIPSO_L2_05KMCPRO,
  CALIPSO_L2_05KMALAY, CALIPSO_L2_05KMCLAY,
  CALIPSO_L2_01KMCLAY, CALIPSO_L2_333MCLAY, CALIPSO_L2_VFM,
  CALIPSO_FILE_TYPES
};

#define IS_CALIPSO( type ) \
  IN9( (type), \
       CALIPSO_L1, CALIPSO_L2_05KMAPRO, CALIPSO_L2_05KMCPRO, \
       CALIPSO_L2_05KMALAY, CALIPSO_L2_05KMCLAY, \
       CALIPSO_L2_01KMCLAY, CALIPSO_L2_333MCLAY, CALIPSO_L2_VFM )

#define IS_LAYERED( type ) \
  IN5( (type), \
       CALIPSO_L2_05KMALAY, CALIPSO_L2_05KMCLAY, CALIPSO_L2_01KMCLAY, \
       CALIPSO_L2_333MCLAY )

/*================================ FUNCTIONS ================================*/

extern int typeOfCALIPSOFile( const char* const fileName );

extern int readCALIPSOVariableDimensions( const int file,
                                          const char* const variable,
                                          size_t* const points,
                                          size_t* const levels );

extern int readCALIPSOData( const int file, const int fileType,
                            const char* const variable,
                            const size_t points, const size_t levels,
                            const double minimumCAD,
                            const double maximumUncertainty,
                            char units[ 80 ],
                            double timestamps[],
                            double longitudes[], double latitudes[],
                            double elevations[], double thicknesses[],
                            double data[] );

extern void aggregateCALIPSOData( const size_t points,
                                  const size_t levels,
                                  const size_t window,
                                  const size_t targetLevels,
                                  double timestamps[],
                                  double longitudes[],
                                  double latitudes[],
                                  double elevations[],
                                  double data[],
                                  size_t* const subsetPoints,
                                  size_t* const subsetLevels );

#ifdef __cplusplus
}
#endif

#endif /* READDATA_H */




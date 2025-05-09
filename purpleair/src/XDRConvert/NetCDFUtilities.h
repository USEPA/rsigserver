#ifndef NETCDFUTILITIES_H
#define NETCDFUTILITIES_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: NetCDFUtilities.h - Declare convenience routines for NetCDF files.

NOTES:

HISTORY: 2007/12, plessel.todd@epa.gov, Created.

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <BasicNumerics.h> /* For Integer, Real. */

/*================================ FUNCTIONS ================================*/


extern Integer createLongitudeAndLatitude( Integer file,
                                           Integer dimensionality,
                                           const Integer dimensions[] );

extern Integer writeExtraAttributes( Integer file, const Real domain[2][2],
                                     Integer dimensionId );

extern Integer writeTimeData( Integer file, Integer count, Integer dims,
                              Integer useBothDims,
                              const Integer* timestamps,
                              const Integer* dimensions,
                              Real* buffer );

extern Integer writeTimeData1( const Integer file, const Integer count,
                               const int yyyyddd[], const int hhmmss[],
                               const float fhour[] );

extern Integer createNetCDFFile( const char* fileName,
                                 Integer create64BitFile );

extern Integer createDimensions( Integer file, Integer count,
                                 const char* const names[],
                                 const Integer values[],
                                 Integer ids[] );

extern Integer writeStandardContents( Integer file, const char* history,
                                      const UTCTimestamp timestamp,
                                      Integer timeDimension,
                                      Integer timesteps, Integer writeTime );

extern Integer createVariable( Integer file, const char* name,
                               const char* units,
                               Integer type, Integer hasMissingValues,
                               Integer dimensionality,
                               const Integer* dimensionIds );

extern Integer writeIntegerAttribute( Integer file, const char* name,
                                      Integer value );

extern Integer writeRealAttribute( Integer file, Integer id, Integer type,
                                   const char* name, Real value );

extern Integer writeTextAttribute( Integer file, Integer id, const char* name,
                                   const char* value );

extern Integer writeRealArrayAttribute( Integer file, Integer type,
                                        const char* name,
                                        const Real* values, Integer count );

extern Integer writeAllData( Integer file, const char* variableName,
                             Integer dimension1, Integer dimension2,
                             Integer dimension3, Integer dimension4,
                             Real* data );

extern Integer writeAllIntData( Integer file, const char* variableName,
                                Integer dimension1, Integer dimension2,
                                Integer dimension3, Integer dimension4,
                                Integer* data );

extern Integer writeAllCharData( Integer file, const char* variableName,
                                 Integer count, Integer length,
                                 const char data[] );

extern Integer writeSomeData( Integer file, const char* variableName,
                              Integer timestep,
                              Integer dimension1, Integer dimension2,
                              Integer dimension3, Integer dimension4,
                              Real* data );

extern Integer writeSomeIntData( Integer file, const char* variableName,
                                 Integer timestep,
                                 Integer dimension1, Integer dimension2,
                                 Integer dimension3, Integer dimension4,
                                 const int* data );

extern Integer writeSomeIntegerData( Integer file, const char* variableName,
                                     Integer timestep,
                                     Integer dimension1, Integer dimension2,
                                     Integer dimension3, Integer dimension4,
                                     Integer* data );

#ifdef __cplusplus
}
#endif

#endif /* NETCDFUTILITIES_H */




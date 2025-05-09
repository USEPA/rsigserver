/******************************************************************************
PURPOSE: ReadData.c - Implement routines for reading ceilometer data and
         subsetting to domain.

NOTES:   Uses HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-11-16 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(), snprintf(). */
#include <string.h> /* For memset(), memcpy(), strcmp(), strstr(). */

#include <Utilities.h> /* For LONGITUDE, Bounds, expand32BitReals(). */
#include "ReadFile.h"  /* For openFile(), readFileData(), openDataset(). */
#include "ReadData.h"  /* For public interface. */

/*================================ CONSTANTS ================================*/

#define MISSING_VALUE (-9999.0)

#define MAXIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE 540000.0

/* Meters above mean sea level: */

#define MINIMUM_SURFACE_ELEVATION (-500.0)
#define MAXIMUM_SURFACE_ELEVATION 1e4
#define MAXIMUM_ELEVATION 1e5

/*=========================== FORWARD DECLARATIONS ==========================*/

static int readCeilometerLocation( const int file,
                                   double* longitude,
                                   double* latitude,
                                   double* elevation );

static int readCeilometerData( const int dataset,
                               const size_t timesteps,
                               const size_t levels,
                               const char* const datasetName,
                               const int isMixingLayerHeight,
                               double data[] );

static int readCeilometerElevations( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
                                     const char* const datasetName,
                                     const double surfaceElevation,
                                     const double minimumElevation,
                                     const double maximumElevation,
                                     const double data[],
                                     double elevations[] );

static int readCeilometerTimestamps( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
                                     const char* const datasetName,
                                     const Integer firstTimestamp,
                                     const Integer lastTimestamp,
                                     double timestamps[] );

static size_t subsetCeilometerData( const size_t count,
                                    double data[],
                                    double elevations[],
                                    double timestamps[] );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: readSubsetCeilometerData - Read a subset of the Ceilometer data that
         are valid and within the subset domain and time-range.
INPUTS:  const char* const* fileName  File to open/read.
         const char* const variable   Name of data variable to read.
         const Bounds domain          Lon-lat-elevation of data subset.
         const Integer firstTimestamp First YYYYMMDDHHMMSS of subset.
         const Integer lastTimestamp  Last  YYYYMMDDHHMMSS of subset.
OUTPUTS: char units[ 80 ]             Units of variable.
         double* longitude            Longitude of site.
         double* latitude             Latitude  of site.
         double* elevation            Elevation of site.
         double** data                Allocated data[ result ].
         double** elevations          Allocated elevations[ result ].
         double** timestamps          Allocated timestamps[ result ].
RETURNS: size_t number of valid data points within domain or 0 if failed.
******************************************************************************/

size_t readSubsetCeilometerData( const char* const fileName,
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
                                 double** timestamps ) {
  size_t result = 0;
  int file = -1;

  assert( fileName ); assert( *fileName );
  assert( variable ); assert( *variable );
  assert( isValidBounds( (const double (*)[2]) domain ) );
  assert( isValidYYYYMMDDHHMMSS( firstTimestamp ) );
  assert( isValidYYYYMMDDHHMMSS( lastTimestamp ) );
  assert( firstTimestamp <= lastTimestamp );
  assert( units );
  assert( longitude );
  assert( latitude );
  assert( elevation );
  assert( data );
  assert( elevations );
  assert( timestamps );

  *longitude = *latitude = *elevation = 0.0;
  *data = *elevations = *timestamps = 0;
  units[ 0 ] = '\0';

  file = openFile( fileName );

  if ( file != -1 ) {
    int ok = readCeilometerLocation( file, longitude, latitude, elevation );

    if ( ok ) {
      const int inDomain =
        AND2( IN_RANGE( *longitude,
                        domain[ LONGITUDE ][ MINIMUM ],
                        domain[ LONGITUDE ][ MAXIMUM ] ),
              IN_RANGE( *latitude,
                        domain[ LATITUDE ][ MINIMUM ],
                        domain[ LATITUDE ][ MAXIMUM ] ) );

      if ( inDomain ) {
        int dataset = -1;
        const char* const datasetPath =
          strstr( variable, "full_" ) == variable ? "/FULL_PROFILE"
          : "/BLVIEW_PROFILE";
        char datasetName[ 256 ] = "";

        /*
         * mixing_layer_height (timesteps x 1) is a pseudo-variable name for
         * the first level of aerosol_layer_heights (timesteps x 3 levels).
         */

        const int isMixingLayerHeight =
          ! strcmp( variable, "mixing_layer_height" );
        const char* const variable0 =
          isMixingLayerHeight ? "aerosol_layer_heights" : variable;

        memset( datasetName, 0, sizeof datasetName );
        snprintf( datasetName, sizeof datasetName / sizeof *datasetName,
                  "%s/%s", datasetPath, variable0 );

        DEBUG( fprintf( stderr, "readSubsetCeilometerData(): "
                        "datasetName = '%s', isMixingLayerHeight = %d\n",
                        datasetName, isMixingLayerHeight ); )

        dataset = openDataset( file, datasetName );

        if ( dataset != -1 ) {
          size_t timesteps = 0;
          size_t levels = 0;
          ok = readDatasetDimensions( dataset, &timesteps, &levels );

          if ( ok ) {
            const size_t count = timesteps * levels;
            *data = NEW( double, count );
            *elevations = *data ? NEW( double, count ) : 0;
            *timestamps = *elevations ? NEW( double, count ) : 0;
            ok = *timestamps != 0;

            if ( ok ) {
              ok = readCeilometerData( dataset, timesteps, levels,
                                       datasetName, isMixingLayerHeight, *data);

              if ( ok ) {
                ok = readCeilometerElevations( file, timesteps, levels,
                                               datasetName,
                                               *elevation,
                                               domain[ ELEVATION ][ MINIMUM ],
                                               domain[ ELEVATION ][ MAXIMUM ],
                                               *data, *elevations );

                if ( ok ) {
                  ok = readCeilometerTimestamps( file, timesteps, levels,
                                                 datasetName,
                                                 firstTimestamp,
                                                 lastTimestamp,
                                                 *timestamps );

                  if ( ok ) {
                   result =
                      subsetCeilometerData( timesteps * levels,
                                            *data, *elevations, *timestamps );

                    if ( result ) {
                      const int isHeight = strstr( variable, "_height" ) != 0;

                      if ( isHeight ) {
                        units[ 0 ] = 'm';
                        units[ 1 ] = '\0';
                      } else {
                        strcpy( units, "10^-9/m/sr" );
                      }
                    }
                  }
                }
              }
            }
          }

          closeDataset( dataset ), dataset = -1;
        }
      }
    }

    closeFile( file ), file = -1;
  }

  if ( ! result ) {
    FREE( *data );
    FREE( *elevations );
    FREE( *timestamps );
    *longitude = *latitude = *elevation = 0.0;
    units[ 0 ] = '\0';
  }

  assert( IMPLIES_ELSE( result,
                        AND13( *units,
                               isValidLongitude( *longitude ),
                               isValidLatitude( *latitude ),
                               IN_RANGE( *elevation, -500.0, 1e5 ),
                               *data,
                               minimumItem( *data, result ) > MISSING_VALUE,
                               maximumItem( *data, result ) <=
                                 MAXIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE,
                               *elevations,
                               minimumItem( *elevations, result ) >= -500.0,
                               maximumItem( *elevations, result ) <= 1e5,
                               *timestamps,
                               isValidYYYYMMDDHHMMSS( (Integer)
                                 minimumItem( *timestamps, result ) ),
                               isValidYYYYMMDDHHMMSS( (Integer)
                                 maximumItem( *timestamps, result ) ) ),
                        IS_ZERO7( *units, *longitude, *latitude, *elevation,
                                  *data, *elevations, *timestamps ) ) );
  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: readCeilometerLocation - Read ceilometer location coordinates.
INPUTS:  const int file      File to read.
OUTPUTS: double* longitude   Longitude of ceilometer.
         double* latitude    Latitude  of ceilometer.
         double* elevation   Elevation of ceilometer.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int readCeilometerLocation( const int file,
                                   double* longitude,
                                   double* latitude,
                                   double* elevation ) {
  int result = 0;

  assert( file > -1 );
  assert( longitude ); assert( latitude ); assert( elevation );

  *longitude = *latitude = *elevation = MISSING_VALUE;

  result = readFileAttribute( file, "/LOCATION/longitude", longitude );

  if ( result ) {
    result = isValidLongitude( *longitude );

    if ( ! result ) {
      fprintf( stderr, "\a\nFailed to read valid file longitude (%lf).\n",
               *longitude );
    } else {
      result = readFileAttribute( file, "/LOCATION/latitude", latitude );

      if ( result ) {
        result = isValidLatitude( *latitude );

        if ( ! result ) {
          fprintf( stderr, "\a\nFailed to read valid file latitude (%lf).\n",
                   *latitude );
        } else {
          result = readFileAttribute( file, "/LOCATION/elevation",
                                      elevation );

          if ( result ) {
            result = IN_RANGE( *elevation,
                               MINIMUM_SURFACE_ELEVATION,
                               MAXIMUM_SURFACE_ELEVATION );
            if ( ! result ) {
              fprintf( stderr,
                       "\a\nFailed to read valid file elevation (%lf).\n",
                       *elevation );
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    *longitude = *latitude = *elevation = MISSING_VALUE;
  }

  assert( IS_BOOL( result ) );
  assert( IMPLIES_ELSE( result,
                        AND3( isValidLongitude( *longitude ),
                              isValidLatitude( *latitude ),
                              IN_RANGE( *elevation,
                                        MINIMUM_SURFACE_ELEVATION,
                                        MAXIMUM_SURFACE_ELEVATION ) ),
                       AND3( *longitude == MISSING_VALUE,
                             *latitude  == MISSING_VALUE,
                             *elevation == MISSING_VALUE ) ) );

  return result;
}



/******************************************************************************
PURPOSE: readCeilometerData - Read ceilometer variable data.
INPUTS:  const int dataset                  Id of dataset to read.
         const size_t timesteps             Number of timesteps to match/read.
         const size_t levels                Number of levels to match/read.
         const char* const datasetName      Pathed name of dataset to read.
         const int isMixingLayerHeight      Is variable mixing_layer_height?
OUTPUTS: double data[ timesteps * levels ]  Data read.
RETURNS: int 1 if successful and at least one valid value was read, else 0.
******************************************************************************/

static int readCeilometerData( const int dataset,
                               const size_t timesteps,
                               const size_t levels,
                               const char* const datasetName,
                               const int isMixingLayerHeight,
                               double data[] ) {
  int result = 0;

  assert( dataset > -1 );
  assert( timesteps );
  assert( levels );
  assert( datasetName ); assert( *datasetName );
  assert( IS_BOOL( isMixingLayerHeight ) );
  assert( data );

  result = readFileData( dataset, timesteps, levels, data );

  if ( result ) { /* Filter-out invalid data values: */
    const int isHeight = strstr( datasetName, "_height" ) != 0;
    const double minimum = isHeight ? MINIMUM_SURFACE_ELEVATION : 0.0;
    const double maximum =
      isHeight ? MAXIMUM_ELEVATION
      : MAXIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE;
    double* output = data;
    const size_t count = timesteps * levels;
    const size_t stride = isMixingLayerHeight ? 3 : 1;
    size_t index = 0;

    DEBUG( double minval = data[ 0 ]; double maxval = minval;
    fprintf( stderr, "%s [%lu x %lu] = %lu, filter range = [%lf, %lf]\n",
              datasetName, timesteps, levels, count, minimum, maximum ); )

    result = 0;

    for ( index = 0; index < count; index += stride, ++output ) {
      const double value = data[ index ];

      DEBUG( if ( index < 15 || index > count - 9 )
               fprintf( stderr, "  data[ %4lu ] = %lf\n", index, value ); )
#ifdef DEBUGGING
      if ( isNan( minval ) ) minval = maxval = value;
      else if ( value < minval ) minval = value;
      else if ( value > maxval ) maxval = value;
#endif

      if ( ! IN_RANGE( value, minimum, maximum ) ) {
        *output = MISSING_VALUE;
      } else {
        *output = value;
        result = 1;
      }
    }

    if ( isMixingLayerHeight ) { /* Filter remaining data levels: */
      const double* const end = data + count;

      while ( output < end ) {
        *output++ = MISSING_VALUE;
      }
    }

    DEBUG( fprintf( stderr,
                    "filtered data[%lu x %lu] range = [%lf, %lf], "
                    "result = %d\n",
                    timesteps, isMixingLayerHeight ? 1 : levels,
                    minval, maxval, result ); )
  }

  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readCeilometerElevations - Read or compute ceilometer elevations
         appropriate for the specified data variable.
INPUTS:  const int file                 Id of file to read.
         const size_t timesteps         Number of data timesteps.
         const size_t levels            Number of data levels.
         const char* const datasetName  Pathed name of dataset.
         const double surfaceElevation  Elevation of site in meters above sea. 
         const double minimumElevation  Minimum elevation of data subset. 
         const double maximumElevation  Maximum elevation of data subset. 
         const double data[ timesteps * levels ]  Variable data read.
OUTPUTS: double elevations[ timesteps * levels ]  Elevations for variable data.
RETURNS: int 1 if successful and at least one valid value was read, else 0.
******************************************************************************/

static int readCeilometerElevations( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
                                     const char* const datasetName,
                                     const double surfaceElevation,
                                     const double minimumElevation,
                                     const double maximumElevation,
                                     const double data[],
                                     double elevations[] ) {
  int result = 0;

  assert( file > -1 ); assert( timesteps ); assert( levels );
  assert( datasetName ); assert( *datasetName );
  assert( IN_RANGE( surfaceElevation,
                    MINIMUM_SURFACE_ELEVATION, MAXIMUM_SURFACE_ELEVATION ) );
  assert( IN_RANGE( minimumElevation,
                    MINIMUM_SURFACE_ELEVATION, MAXIMUM_ELEVATION ) );
  assert( IN_RANGE( maximumElevation,
                    minimumElevation, MAXIMUM_ELEVATION ) );
  assert( data ); assert( elevations );

  if ( strstr( datasetName, "_heights" ) ) {
    const size_t count = timesteps * levels;
    size_t index = 0;

    for ( index = 0; index < count; ++index ) {
      elevations[ index ] = data[ index ] + surfaceElevation;
    }

    result = 1;
  } else { /* Read "altitude" and expand it to all timesteps: */
    char elevationDatasetName[ 256 ] = "";
    char* slash = 0;
    memset( elevationDatasetName, 0, sizeof elevationDatasetName );
    strncpy( elevationDatasetName, datasetName,
             sizeof elevationDatasetName / sizeof *elevationDatasetName - 1 );
    slash = strrchr( elevationDatasetName, '/' );

    if ( slash ) {
      int dataset = -1;
      const size_t length =
        sizeof elevationDatasetName / sizeof *elevationDatasetName -
        (size_t) ( slash - elevationDatasetName );
      strncpy( slash, "/altitude", length - 1 );

      assert( elevationDatasetName[ sizeof elevationDatasetName /
                                    sizeof *elevationDatasetName - 1] == '\0');

      dataset = openDataset( file, elevationDatasetName );

      if ( dataset > -1 ) {
        result = readFileData( dataset, levels, 0, elevations );
        closeDataset( dataset ), dataset = -1;

        /*
         * Convert from meters above ground to meters above mean sea level
         * (filtering invalid values)
         * then copy to all timesteps:
         */

        if ( result ) {
          size_t level = 0;
          result = 0;

          DEBUG( fprintf( stderr, "readCeilometerElevations(): "
                          "surfaceElevation = %lf, "
                          "minimumElevation = %lf, maximumElevation = %lf\n",
                          surfaceElevation,
                          minimumElevation, maximumElevation ); )

          for ( level = 0; level < levels; ++level ) {
            const double elevationAboveMeanSeaLevel =
              ( elevations[ level ] += surfaceElevation );

            if ( ! IN_RANGE( elevationAboveMeanSeaLevel,
                             minimumElevation, maximumElevation ) ) {
              elevations[ level ] = MISSING_VALUE;
            } else {
              result = 1;
            }

            DEBUG( if ( level < 5 || level > levels - 5 )
                     fprintf( stderr, "  elevation[ %lu ] = %lf\n",
                              level, elevations[ level ] ); )
          }

          /* If valid elevation was read then copy to remaining timesteps: */

          if ( result ) {
            double* copy = elevations + levels;
            const size_t bytes = levels * sizeof (double);
            size_t timestep = 1;

            for ( ; timestep < timesteps; ++timestep, copy += levels ) {
              memcpy( copy, elevations, bytes );
            }
          }
        }
      }
    }
  }

  DEBUG( fprintf( stderr, "readCeilometerElevations() result = %d\n", result);)
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readCeilometerTimestamps - Read ceilometer timestamps.
INPUTS:  const int file                 File to read.
         const size_t timesteps         Number of timesteps to match/read.
         const size_t levels            Number of levels to match/read.
         const char* const datasetName  Pathed name of dataset to read.
         const Integer firstTimestamp   First YYYYMMDDHHMMSS of subset.
         const Integer lastTimestamp    Last  YYYYMMDDHHMMSS of subset.
OUTPUTS: double timestamps[ timesteps * levels ]  YYYYMMDDHHMMSS values.
RETURNS: int 1 if successful and at least one valid in-subset value was read,
         else 0.
******************************************************************************/

static int readCeilometerTimestamps( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
                                     const char* const datasetName,
                                     const Integer firstTimestamp,
                                     const Integer lastTimestamp,
                                     double timestamps[] ) {

  int result = 0;
  char yyyymmddDatasetName[ 256 ] = "";
  char hhmmssDatasetName[ 256 ] = "";
  char* yyyymmddSlash = 0;
  char* hhmmssSlash = 0;
  long long* yyyymmdd_hhmmss = 0;

  assert( file > -1 ); assert( timesteps ); assert( levels );
  assert( datasetName ); assert( *datasetName );
  assert( isValidYYYYMMDDHHMMSS( firstTimestamp ) );
  assert( isValidYYYYMMDDHHMMSS( lastTimestamp ) );
  assert( firstTimestamp <= lastTimestamp );
  assert( timestamps );

  yyyymmdd_hhmmss = NEW( long long, timesteps * 2 );

  if ( yyyymmdd_hhmmss ) {
    long long* const yyyymmdd0 = yyyymmdd_hhmmss;
    long long* const hhmmss0 = yyyymmdd0 + timesteps;

    /* Read yyyymmdd and hhmmss arrays and concatenate their values: */

    memset( yyyymmddDatasetName, 0, sizeof yyyymmddDatasetName );
    memset( hhmmssDatasetName, 0, sizeof hhmmssDatasetName );
    strncpy( yyyymmddDatasetName, datasetName,
             sizeof yyyymmddDatasetName / sizeof *yyyymmddDatasetName - 1 );
    strncpy( hhmmssDatasetName, datasetName,
             sizeof hhmmssDatasetName / sizeof *hhmmssDatasetName - 1 );
    yyyymmddSlash = strrchr( yyyymmddDatasetName, '/' );
    hhmmssSlash = strrchr( hhmmssDatasetName, '/' );

    if ( AND2( yyyymmddSlash, hhmmssSlash ) ) {
      int yyyymmddDatasetId = -1;
      int hhmmssDatasetId = -1;
      const size_t yyyymmddLength =
        sizeof yyyymmddDatasetName / sizeof *yyyymmddDatasetName -
          (size_t) ( yyyymmddSlash - yyyymmddDatasetName );
      const size_t hhmmssLength =
        sizeof hhmmssDatasetName / sizeof *hhmmssDatasetName -
          (size_t) ( hhmmssSlash - hhmmssDatasetName );
      strncpy( yyyymmddSlash, "/yyyymmdd", yyyymmddLength - 1 );
      strncpy( hhmmssSlash, "/hhmmss", hhmmssLength - 1 );

      assert( yyyymmddDatasetName[ sizeof yyyymmddDatasetName /
                                   sizeof *yyyymmddDatasetName - 1 ] == '\0' );

      assert( hhmmssDatasetName[ sizeof hhmmssDatasetName /
                                 sizeof *hhmmssDatasetName - 1 ] == '\0' );

      yyyymmddDatasetId = openDataset( file, yyyymmddDatasetName );

      if ( yyyymmddDatasetId > -1 ) {
        result =
          readFileDataIntegers( yyyymmddDatasetId, timesteps, 0, yyyymmdd0 );
        closeDataset( yyyymmddDatasetId ), yyyymmddDatasetId = -1;
      }
  
      if ( result ) {
        result = 0;
        hhmmssDatasetId = openDataset( file, hhmmssDatasetName );

        if ( hhmmssDatasetId > -1 ) {
          result =
            readFileDataIntegers( hhmmssDatasetId, timesteps, 0, hhmmss0 );
          closeDataset( hhmmssDatasetId ), hhmmssDatasetId = -1;
        }
      }

      if ( result ) { /* Check and concatenate values: */
        size_t validCount = 0;
        size_t timestep = 0;

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const long long yyyymmdd = yyyymmdd0[ timestep ];
          const long long hhmmss = hhmmss0[ timestep ];
          const long long yyyymmddhhmmss = yyyymmdd * 1000000LL + hhmmss;
          int ok = isValidYYYYMMDDHHMMSS( yyyymmddhhmmss );

          if ( timestep > 0 ) {
            ok = yyyymmddhhmmss >= (long long) timestamps[ timestep - 1 ];
          }

          if ( ok ) {
            ok = IN_RANGE( yyyymmddhhmmss, firstTimestamp, lastTimestamp );
          }

          if ( ok ) {
            timestamps[ timestep ] = yyyymmddhhmmss;
            ++validCount;
          } else {
            timestamps[ timestep ] = MISSING_VALUE;
          }

          DEBUG( if ( timestep < 5 || timestep > timesteps - 5 )
                   fprintf( stderr,
                            "  yyyymmddhhmmss[ %lu ] = %lld, ok = %d\n",
                              timestep,
                              (long long) timestamps[ timestep ],
                              ok ); )
        }

        /* Copy timestamps to all vertical levels: */

        if ( validCount ) {
          double* output = timestamps + timesteps;

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double value = timestamps[ timestep ];
            size_t count = levels - 1;

            while ( count-- ) {
              *output++ = value;
            }
          }
        }

        DEBUG( fprintf( stderr, "readCeilometerTimestamps() validCount = %lu\n",
                        validCount ); )
        result = validCount > 0;
      }
    }
 
    FREE( yyyymmdd_hhmmss );
  }

  DEBUG( fprintf( stderr, "readCeilometerTimestamps() result = %d\n", result);)
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: subsetCeilometerData - Subset/filter ceilometer data by
         MISSING_VALUE.
INPUTS:  const size_t count          Number of data values.
         double data[       count ]  Data to be subset.
         double elevations[ count ]  Elevations to be subset.
         double timestamps[ count ]  Timestamps to be subset.
OUTPUTS: double data[       count ]  Subsetted data
         double elevations[ count ]  Subsetted elevations
         double timestamps[ count ]  Subsetted timestamps.
RETURNS: size_t Number of valid/filtered subset points if successful, else 0.
******************************************************************************/

static size_t subsetCeilometerData( const size_t count,
                                    double data[],
                                    double elevations[],
                                    double timestamps[] ) {

  size_t result = 0;
  size_t validCount = 0;
  size_t index = 0;

  assert( count );
  assert( data );
  assert( elevations );
  assert( timestamps );

  /*
   * If any indexed value of data, elevations or timestamps is MISSING_VALUE
   * then set all indexed values to MISSING_VALUE
   * else copy valid values to front part of arrays:
   */

  for ( index = 0; index < count; ++index ) {
    const double value     = data[ index ];
    const double elevation = elevations[ index ];
    const double timestamp = timestamps[ index ];
    const int ok =
      AND3( value     > MISSING_VALUE,
            elevation > MISSING_VALUE,
            timestamp > MISSING_VALUE );

    if ( ! ok ) {
      data[ index ] = elevations[ index ] = timestamps[ index] = MISSING_VALUE;
    } else {

      if ( index != validCount ) {
        data[ validCount ] = value;
        elevations[ validCount ] = elevation;
        timestamps[ validCount ] = timestamp;
      }

      DEBUG( if ( validCount < 10 )
               fprintf( stderr, "data[ %lu ] = %lf\n",
                        validCount, data[ validCount ] ); )
      ++validCount;
    }
  }

  result = validCount;

  DEBUG( fprintf( stderr,
                  "subsetCeilometerData() result = %lu, data = [%lf ... %lf]\n",
                  result, data[ 0 ], data[ result - 1 ] ); )
  assert( minimumItem( data, result ) > MISSING_VALUE );
  return result;
}




/******************************************************************************
PURPOSE: ReadData.c - Implement routines for reading ceilometer data and
         subsetting to domain.

NOTES:   Uses NetCDF4 library and libs it depends on (HDF5, curl, z, dl).

HISTORY: 2022-07-12 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(), snprintf(). */
#include <string.h> /* For memset(), memcpy(), strcmp(), strstr(). */

#include <Utilities.h> /* For LONGITUDE, Bounds, expand32BitReals(). */
#include "ReadFile.h"  /* For openFile(), readFileData(). */
#include "ReadData.h"  /* For public interface. */

/*=========================== FORWARD DECLARATIONS ==========================*/

static int readCeilometerLocation( const int file,
                                   double* longitude,
                                   double* latitude,
                                   double* elevation );

static int readCeilometerData( const int file,
                               const int id,
                               const size_t timesteps,
                               const size_t levels,
                               const double validMinimum,
                               const double validMaximum,
                               double data[] );

static int readCeilometerElevations( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
                                     const double surfaceElevation,
                                     const double minimumElevation,
                                     const double maximumElevation,
                                     const int isHeight,
                                     const double data[],
                                     double elevations[] );

static int readCeilometerTimestamps( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
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
        const char* const fileVariable =
          ! strcmp( variable, "mixing_layer_height" ) ? "cML"
          : ! strcmp( variable, "cloud_base_heights" ) ? "Cloud_base_heights"
          : variable;
        const int id = readVariableId( file, fileVariable );

        DEBUG( fprintf( stderr, "readSubsetCeilometerData(): "
                        "fileVariable = '%s', id = %d\n", fileVariable, id ); )

        if ( id != -1 ) {
          size_t timesteps = 0;
          size_t levels = 0;
          ok = readVariableDimensions( file, id, &timesteps, &levels );

          if ( ok ) {
            const size_t count = timesteps * levels;
            *data = NEW( double, count );
            *elevations = *data ? NEW( double, count ) : 0;
            *timestamps = *elevations ? NEW( double, count ) : 0;
            ok = *timestamps != 0;

            if ( ok ) {
              const int isHeight = strstr( variable, "scatter" ) == 0;
              double validMinimum = 0.0;
              double validMaximum = 0.0;

              /* UGLY: The units and valid range should be a file attribute! */

              if ( isHeight ) {
                units[ 0 ] = 'm';
                units[ 1 ] = '\0';
                validMinimum = MINIMUM_SURFACE_ELEVATION;
                validMaximum = MAXIMUM_ELEVATION;
              } else {
                strcpy( units, "/km/sr" );
                validMinimum = MINIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE;
                validMaximum = MAXIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE;
              }

              ok = readCeilometerData( file, id, timesteps, levels,
                                       validMinimum, validMaximum, *data );

              if ( ok ) {
                ok = readCeilometerElevations( file, timesteps, levels,
                                               *elevation,
                                               domain[ ELEVATION ][ MINIMUM ],
                                               domain[ ELEVATION ][ MAXIMUM ],
                                               isHeight,
                                               *data, *elevations );

                if ( ok ) {
                  ok = readCeilometerTimestamps( file, timesteps, levels,
                                                 firstTimestamp,
                                                 lastTimestamp,
                                                 *timestamps );

                  if ( ok ) {
                   result =
                      subsetCeilometerData( timesteps * levels,
                                            *data, *elevations, *timestamps );
                  }
                }
              }
            }
          }
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
                               IN_RANGE( *elevation, MINIMUM_SURFACE_ELEVATION,
                                         MAXIMUM_ELEVATION ),
                               *data,
                               minimumItem( *data, result ) > MISSING_VALUE,
                               maximumItem( *data, result ) <=
                                 MAXIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE,
                               *elevations,
                               minimumItem( *elevations, result ) >=
                                 MINIMUM_SURFACE_ELEVATION,
                               maximumItem( *elevations, result ) <=
                                 MAXIMUM_ELEVATION,
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

  result = readFileAttribute( file, "Longitude", longitude );

  if ( result ) {
    result = isValidLongitude( *longitude );

    if ( ! result ) {
      fprintf( stderr, "\a\nFailed to read valid file longitude (%lf).\n",
               *longitude );
    } else {
      result = readFileAttribute( file, "Latitude", latitude );

      if ( result ) {
        result = isValidLatitude( *latitude );

        if ( ! result ) {
          fprintf( stderr, "\a\nFailed to read valid file latitude (%lf).\n",
                   *latitude );
        } else {
          result = readFileAttribute( file, "Elevation ASL", elevation);/*UGLY*/

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
INPUTS:  const int file                     Id of file to read.
         const int id                       Id of variable to read.
         const size_t timesteps             Number of timesteps to match/read.
         const size_t levels                Number of levels to match/read.
         const double validMinimum          Minimum data value that is valid.
         const double validMaximum          Maximum data value that is valid.
OUTPUTS: double data[ timesteps * levels ]  Data read.
RETURNS: int 1 if successful and at least one valid value was read, else 0.
******************************************************************************/

static int readCeilometerData( const int file,
                               const int id,
                               const size_t timesteps,
                               const size_t levels,
                               const double validMinimum,
                               const double validMaximum,
                               double data[] ) {
  int result = 0;
  DEBUG( size_t filteredCount = 0; )

  assert( file > -1 );
  assert( id > -1 );
  assert( timesteps );
  assert( levels );
  assert( data );
  assert( validMinimum < validMaximum );

  result = readFileData( file, id, timesteps, levels, data );

  if ( result ) { /* Filter-out invalid data values: */
    const size_t count = timesteps * levels;
    size_t index = 0;

    DEBUG( fprintf( stderr, "data[%lu x %lu] = %lu, filter = [%lf, %lf]\n",
                    timesteps, levels, count,
                    validMinimum, validMaximum ); )

    result = 0;

    for ( index = 0; index < count; ++index ) {
      const double value = data[ index ];

      if ( ! IN_RANGE( value, validMinimum, validMaximum ) ) {
        data[ index ] = MISSING_VALUE;
        DEBUG( ++filteredCount; )
      } else {
        result = 1;
      }
    }
  }

  DEBUG( fprintf( stderr, "readCeilometerData() filteredCount = %lu, "
                  "data = [%f ... %f], result = %d\n",
                  filteredCount,
                  data[ 0 ], data[ timesteps * levels - 1 ],
                  result ); )
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readCeilometerElevations - Read or compute ceilometer elevations
         appropriate for the specified data variable.
INPUTS:  const int file                 Id of file to read.
         const size_t timesteps         Number of data timesteps.
         const size_t levels            Number of data levels.
         const double surfaceElevation  Elevation of site in meters above sea.
         const double minimumElevation  Minimum elevation of data subset. 
         const double maximumElevation  Maximum elevation of data subset.
         const int isHeight             Is the data heights?
         const double data[ timesteps * levels ]  Variable data read.
OUTPUTS: double elevations[ timesteps * levels ]  Elevations for variable data.
RETURNS: int 1 if successful and at least one valid value was read, else 0.
******************************************************************************/

static int readCeilometerElevations( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
                                     const double surfaceElevation,
                                     const double minimumElevation,
                                     const double maximumElevation,
                                     const int isHeight,
                                     const double data[],
                                     double elevations[] ) {
  int result = 0;
  DEBUG( size_t filteredCount = 0; )

  assert( file > -1 ); assert( timesteps ); assert( levels );
  assert( IN_RANGE( surfaceElevation,
                    MINIMUM_SURFACE_ELEVATION, MAXIMUM_SURFACE_ELEVATION ) );
  assert( IN_RANGE( minimumElevation,
                    MINIMUM_SURFACE_ELEVATION, MAXIMUM_ELEVATION ) );
  assert( IN_RANGE( maximumElevation,
                    minimumElevation, MAXIMUM_ELEVATION ) );
  assert( data ); assert( elevations );

  if ( isHeight ) {
    const size_t count = timesteps * levels;
    size_t index = 0;

    for ( index = 0; index < count; ++index ) {
      elevations[ index ] = data[ index ] + surfaceElevation;
    }

    result = 1;
  } else { /* Read "range" and expand it to all timesteps: */
    const int rangeId = readVariableId( file, "range" );

    if ( rangeId > -1 ) {
      result = readFileData( file, rangeId, levels, 1, elevations );

      /*
       * Convert from meters above ground to meters above mean sea level
       * (filtering invalid values) then copy to all timesteps:
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
            DEBUG( ++filteredCount; )
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

  DEBUG( fprintf( stderr, "readCeilometerElevations() filteredCount = %lu, "
                  "elevations = [%f ... %f], result = %d\n",
                  filteredCount,
                  elevations[ 0 ], elevations[ timesteps * levels - 1 ],
                  result ); )
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readCeilometerTimestamps - Read ceilometer timestamps.
INPUTS:  const int file                 File to read.
         const size_t timesteps         Number of timesteps to match/read.
         const size_t levels            Number of levels to match/read.
         const Integer firstTimestamp   First YYYYMMDDHHMMSS of subset.
         const Integer lastTimestamp    Last  YYYYMMDDHHMMSS of subset.
OUTPUTS: double timestamps[ timesteps * levels ]  YYYYMMDDHHMMSS values.
RETURNS: int 1 if successful and at least one valid in-subset value was read,
         else 0.
******************************************************************************/

static int readCeilometerTimestamps( const int file,
                                     const size_t timesteps,
                                     const size_t levels,
                                     const Integer firstTimestamp,
                                     const Integer lastTimestamp,
                                     double timestamps[] ) {

  int result = 0;
  const int id = readVariableId( file, "time" );
  DEBUG( size_t filteredCount = 0; )

  if ( id > -1 ) {
    result = readFileData( file, id, timesteps, 1, timestamps );

    if ( result ) { /* Convert seconds since 1970 to yyyymmddhhmmss: */
      Integer validCount = 0;
      size_t timestep = 0;

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        const double seconds = timestamps[ timestep ];
        const Integer iseconds = (Integer) seconds;
        const Integer maximumSeconds =
          ( 2100 - 1970 ) * 365LL * 24 * 60 * 60;
        int ok = 0;

        if ( IN_RANGE( iseconds, 0, maximumSeconds ) ) {
          const Integer yyyymmddhhmmss = fromSeconds( iseconds );
          ok = isValidYYYYMMDDHHMMSS( yyyymmddhhmmss );

          if ( ok ) {

            if ( timestep > 0 ) {
              ok = yyyymmddhhmmss >= (Integer) timestamps[ timestep - 1 ];
            }

            if ( ok ) {
              ok = IN_RANGE( yyyymmddhhmmss, firstTimestamp, lastTimestamp );
            }

            if ( ok ) {
              timestamps[ timestep ] = yyyymmddhhmmss;
              ++validCount;
            }
          }
        }

        if ( ! ok )  {
          timestamps[ timestep ] = MISSING_VALUE;
          DEBUG( ++filteredCount; )
        }

        DEBUG( if ( timestep < 5 || timestep > timesteps - 5 )
                 fprintf( stderr,
                          "  yyyymmddhhmmss[ %3lu ] = %lld, ok = %d\n",
                          timestep,
                          (Integer) timestamps[ timestep ],
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

      DEBUG( fprintf( stderr, "readCeilometerTimestamps() validCount = %lld\n",
                      validCount ); )
      result = validCount > 0;
    }
  }

  DEBUG( fprintf( stderr, "readCeilometerTimestamps() filteredCount = %lu, "
                  "timestamps = [%f ... %f], result = %d\n",
                  filteredCount,
                  timestamps[ 0 ], timestamps[ timesteps * levels - 1 ],
                  result ); )
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
  DEBUG( size_t filteredCount = 0; )

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
      DEBUG( ++filteredCount; )
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
                  "subsetCeilometerData() filteredCount = %lu, "
                  "data = [%lf ... %lf], result = %lu\n",
                  filteredCount,
                  data[ 0 ], data[ result - 1 ],
                  result ); )
  assert( result == 0 || minimumItem( data, result ) > MISSING_VALUE );
  return result;
}




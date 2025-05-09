/******************************************************************************
PURPOSE: ReadData.c - routine to read and filter/process CALIPSO data.

NOTES:   Uses routines in ReadFile.c and Utilities.c.

HISTORY: 2017-01-02 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For strcmp(), strstr(), strncpy(). */

#include "Utilities.h" /* For MISSING_VALUE, LONGITUDE, IN_RANGE(), Bounds. */
#include "ReadFile.h"  /* For readFileData(), readFileVData(). */
#include "ReadData.h"  /* For public interface. */

/*=========================== FORWARD DECLARATIONS ==========================*/

static int readCALIPSOVariable( const int file, const int fileType,
                                const char* const variable,
                                const size_t points, const size_t levels,
                                char units[ 80 ], double data[] );

static int readCALIPSOTimestamps( const int file, const size_t points,
                                  double timestamps[] );

static int readCALIPSOCoordinates( const int file, const int fileType,
                                   const size_t points,
                                   double longitudes[], double latitudes[] );

static int readCALIPSOElevations( const int file, const int fileType,
                                  const size_t points, const size_t levels,
                                  double elevations[], double thicknesses[] );

static int filterCALIPSOData( const int file, const int fileType,
                              const char* const variable,
                              const size_t points, const size_t levels,
                              const double minimumCADScore,
                              const double maximumUncertainty,
                              const double elevations[], double data[] );

static void copyWorstCADScore( const size_t count, const double input[],
                               double output[] );

static const char* uncertaintyVariable( const char* const variable );

static void filterDataNearSurface( const size_t points, const size_t levels,
                                   const double elevations[], double data[] );

static int filterDataByQC( const int file, const int fileType,
                           const double dataMinimum, const double dataMaximum,
                           const double missingValue,
                           const char* const qcVariable, const size_t qcLevels,
                           const double qcMinimum, const double qcMaximum,
                           const unsigned int mask, const int propagateDown,
                           const size_t points, const size_t levels,
                           double data[] );

static void propagateBadUncertaintyDownward( const size_t points,
                                             const size_t levels,
                                             double uncertainty[] );

static void propagateWorstCADScoreDownward( const size_t points,
                                            const size_t levels,
                                            double score[] );

static void propagateWorstQCFlagDownward( const unsigned int mask,
                                          const size_t points,
                                          const size_t levels,
                                          double qc[] );

static size_t filterData( const double dataMinimum, const double dataMaximum,
                          const double missingValue,
                          const size_t points, const size_t levels,
                          const size_t qcLevels,
                          const double qcFlags[],
                          const double qcMinimum, const double qcMaximum,
                          const unsigned int mask,
                          double data[] );

static size_t computeAggregateLevels( const size_t levels,
                                      const size_t targetLevels,
                                      size_t* const stride );

static size_t computeStride( const size_t count, const size_t target );

static void aggregateDataAndElevations( const size_t points,
                                        const size_t levels,
                                        const size_t point,
                                        const size_t level,
                                        const size_t width,
                                        const size_t height,
                                        const double data[],
                                        const double elevations[],
                                        double* const meanDatum,
                                        double* const meanElevation );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: typeOfCALIPSOFile - Type of CALIPSO file, given name.
INPUTS:  const char* const fileName      Name of file.
RETURNS: int type of CALIPSO file or -1 if invalid.
******************************************************************************/

int typeOfCALIPSOFile( const char* const fileName ) {
  int result = -1;
  assert( fileName ); assert( *fileName );

  if ( strstr( fileName, "CAL_LID_L1" ) ) {
    result = CALIPSO_L1;
  } else if ( strstr( fileName, "CAL_LID_L2_05kmAPro" ) ) {
      result = CALIPSO_L2_05KMAPRO;
  } else if ( strstr( fileName, "CAL_LID_L2_05kmCPro" ) ) {
      result = CALIPSO_L2_05KMCPRO;
  } else if ( strstr( fileName, "CAL_LID_L2_05kmALay" ) ) {
      result = CALIPSO_L2_05KMALAY;
  } else if ( strstr( fileName, "CAL_LID_L2_05kmCLay" ) ) {
      result = CALIPSO_L2_05KMCLAY;
  } else if ( strstr( fileName, "CAL_LID_L2_01kmCLay" ) ) {
      result = CALIPSO_L2_01KMCLAY;
  } else if ( strstr( fileName, "CAL_LID_L2_333mCLay" ) ) {
      result = CALIPSO_L2_333MCLAY;
  } else if ( strstr( fileName, "CAL_LID_L2_VFM" ) ) {
      result = CALIPSO_L2_VFM;
  }

  return result;
}



/******************************************************************************
PURPOSE: readCALIPSOVariableDimensions - Read CALIPSO variable dimensions.
INPUTS:  const int file             ID of file to read.
         const char* const variable Name of variable to read.
OUTPUTS: size_t* const points       Number of ground points of data.
         size_t* const levels       Number of vertical levels of data.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readCALIPSOVariableDimensions( const int file,
                                   const char* const variable,
                                   size_t* const points,
                                   size_t* const levels ) {

  int result = 0;

  assert( file >= 0 ); assert( variable ); assert( *variable );
  assert( points ); assert( levels );

  *points = *levels = 0;

  {
    const size_t length = strlen( variable );
    const int isVector =
      AND3( length > 2,
            IN4( variable[ length - 1 ], 'X', 'Y', 'Z' ),
            variable[ length - 2 ] == '_' );
    enum { NAME_LENGTH = 255 };
    char variable2[ NAME_LENGTH + 1 ] = "";
    int rank = 0;
    int dimensions[ 32 ] = { 0, 0 };

    memset( dimensions, 0, sizeof dimensions );
    memset( variable2, 0, sizeof variable2 );
    strncpy( variable2, variable,
             sizeof variable2 / sizeof *variable2 - 1 );
    
    if ( isVector ) {
      variable2[ length - 1 ] = '\0'; /* Erase component _X. */
      variable2[ length - 2 ] = '\0';
    }

    result = readVariableDimensions( file, variable2, &rank, dimensions );

    if ( result ) {

      if ( ! IN3( rank, 2, 3 ) ) {
        fprintf( stderr, "\a\n\nInvalid rank for variable '%s' "
                 "(rank = %d, expected 2 or 3).\n",
                 variable2, rank );
        result = 0;
      } else {

        DEBUG( fprintf( stderr, "%d x %d isVector = %d\n",
                        dimensions[ 0 ], dimensions[ 1 ], isVector ); )
        *points = dimensions[ 0 ];
        *levels = isVector ? 1 : dimensions[ 1 ];
      }
    }
  }

  assert( IMPLIES_ELSE( result,
                        AND2( *points, *levels ),
                        IS_ZERO2( *points, *levels ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readCALIPSOData - Read, filter and process CALIPSO data for variable.
INPUTS:  const int file             ID of file to read.
         const int fileType         CALIPSO_L1 ... CALIPSO_L2_VFM.
         const char* const variable Name of variable to read.
         const size_t points        Number of ground points of data.
         const size_t levels        Number of vertical levels of data.
         const double minimumCAD    Minimum cloud/aerosol discrimination score.
                                    E.g., 20 accepts score in range [20, 100].
         const double maximumUncertainty  Maximum uncertainty acceptable.
                                    E.g., 99 accepts absolute uncertainty <= 99.
                                    Units are the same as the data variable.
OUTPUTS: char units[ 80 ]              Units of variable.
         double timestamps[  points ]  Timestamps of data points.
         double longitudes[  points ]  Longitudes of data points.
         double latitudes[   points ]  Latitudes  of data points.
         double elevations[  points * levels ] Elevations (mAMSL) of data points
         double thicknesses[ points * levels ] 0 or layer thicknesses of data.
         double data[        points * levels ] Values of data points.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readCALIPSOData( const int file, const int fileType,
                     const char* const variable,
                     const size_t points, const size_t levels,
                     const double minimumCAD,
                     const double maximumUncertainty,
                     char units[ 80 ],
                     double timestamps[],
                     double longitudes[], double latitudes[],
                     double elevations[], double thicknesses[],
                     double data[] ) {

  int result = 0;

  assert( file >= 0 ); assert( IS_CALIPSO( fileType ) );
  assert( IMPLIES( AND2( IS_LAYERED( fileType ), levels > 1 ), thicknesses ) );
  assert( variable ); assert( *variable );
  assert( points ); assert( levels );
  assert( units );
  assert( timestamps ); assert( longitudes ); assert( latitudes );
  assert( elevations ); assert( data );

  DEBUG( fprintf( stderr, "readCALIPSOData( %s: %lu x %lu )\n",
                  variable, points, levels ) );

  result = readCALIPSOTimestamps( file, points, timestamps );

  if ( result ) {
    result =
      readCALIPSOCoordinates( file, fileType, points, longitudes, latitudes );

    if ( result ) {
      result =
        readCALIPSOElevations( file, fileType, points, levels,
                               elevations, thicknesses );

      if ( result ) {
        result =
          readCALIPSOVariable( file, fileType, variable, points, levels,
                               units, data );
        if ( result ) {
          result = filterCALIPSOData( file, fileType, variable,
                                      points, levels,
                                      minimumCAD, maximumUncertainty,
                                      elevations, data );
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: aggregateCALIPSOData - Aggregate CALIPSO_L1 data from 333m to 5km
         ground point resolution using window = 15.
INPUTS:  const size_t points                  Number of ground points.
         const size_t levels                  Number of elevation levels.
         const size_t window                  # of ground points to aggregate.
         const size_t targetLevels            Desired number of vertical levels
                                              after aggregation.
         double timestamps[ points ]          Timetamps to aggregate.
         double longitudes[ points ]          Longitudes to aggregate.
         double latitudes[  points ]          Latitudes to aggregate.
         double elevations[ points * levels ] Elevations to aggregate.
         double data[ points * levels ]       Data to aggregate.
OUTPUTS: double timestamps[ subsetPoints ] Aggregated timetamps.
         double longitudes[ subsetPoints ] Aggregated longitudes.
         double latitudes[  subsetPoints ] Aggregated latitudes.
         double elevations[ subsetPoints * subsetLevels ] Aggregated elevations.
         double data[ subsetPoints * subsetLevels ] Aggregated data.
         size_t* const subsetPoints      Number of aggregated ground points.
         size_t* const levels            Number of aggregated elevation levels.
******************************************************************************/

void aggregateCALIPSOData( const size_t points,
                           const size_t levels,
                           const size_t window,
                           const size_t targetLevels,
                           double timestamps[],
                           double longitudes[],
                           double latitudes[],
                           double elevations[],
                           double data[],
                           size_t* const subsetPoints,
                           size_t* const subsetLevels ) {

  assert( points ); assert( levels ); assert( window ); assert( targetLevels );
  assert( timestamps ); assert( longitudes ); assert( latitudes );
  assert( elevations ); assert( data );
  assert( subsetPoints ); assert( subsetLevels );

  *subsetPoints = points;
  *subsetLevels = levels;

  if ( window > 1 ) {
    size_t levelStride = 0;
    const size_t aggregateLevels =
      computeAggregateLevels( levels, targetLevels, &levelStride );
    const size_t remainder = ( points % window ) != 0;
    const size_t aggregatePoints = points / window + remainder;

    DEBUG( fprintf( stderr, "levels = %lu, targetLevels = %lu, "
                    "levelStride = %lu, aggregateLevels = %lu\n",
                    levels, targetLevels, levelStride, aggregateLevels ); )
    DEBUG( fprintf( stderr, "points = %lu, window = %lu, "
                    "aggregatePoints = %lu\n",
                    points, window, aggregatePoints ); )

    size_t point = 0;
    size_t level = 0;
    size_t index = 0;
    size_t index2 = 0;

    for ( point = 0, index = 0; point < points; point += window, ++index ) {
      const size_t width = point + window < points ? window : points - point;
      const size_t middlePoint = point + width / 2;
      timestamps[ index ] = timestamps[ middlePoint ];
      longitudes[ index ] = longitudes[ middlePoint ];
      latitudes[  index ] = latitudes[  middlePoint ];

      for ( level = 0; level < levels; level += levelStride, ++index2 ) {
        const size_t height =
          level + levelStride < levels ? levelStride : levels - level;
        double meanDatum = 0.0;
        double meanElevation = 0.0;
        aggregateDataAndElevations( points, levels, point, level,
                                    width, height, data, elevations,
                                    &meanDatum, &meanElevation );
        assert( level < levels ); assert( index2 <= points * levels );
        data[       index2 ] = meanDatum;
        elevations[ index2 ] = meanElevation;
      }
    }

    DEBUG( fprintf( stderr, "aggregatePoints = %lu, aggregateLevels = %lu, "
                    "levelStride = %lu, index = %lu, index2 = %lu\n",
                    aggregatePoints, aggregateLevels,
                    levelStride, index, index2 ); )

    assert( index == aggregatePoints );
    assert( index2 == aggregatePoints * aggregateLevels );
    *subsetPoints = aggregatePoints;
    *subsetLevels = aggregateLevels;
  }
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: readCALIPSOVariable - Read CALIPSO variable.
INPUTS:  const int file             ID of file to read.
         const int fileType         CALIPSO_L1 ... CALIPSO_L2_VFM.
         const char* const variable Name of variable to read.
                                    E.g., "Total_Attenuated_Backscatter_532" or
                                    "Surface_Wind_X" X-component of wind.
         const size_t points        Number of ground points of data.
         const size_t levels        Number of vertical levels of data.
OUTPUTS: char units[ 80 ]           Units of variable.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int readCALIPSOVariable( const int file, const int fileType,
                                const char* const variable,
                                const size_t points, const size_t levels,
                                char units[ 80 ], double data[]  ) {

  int result = 0;

  assert( file >= 0 ); assert( IS_CALIPSO( fileType ) );
  assert( variable ); assert( *variable );
  assert( points ); assert( levels );
  assert( units );  assert( data );

  DEBUG( fprintf( stderr, "readCALIPSOVariable( %s, %lu x %lu )\n",
                  variable, points, levels ) );

  {
    const char* const xp = strstr( variable, "_X" );
    const char* const yp = strstr( variable, "_Y" );
    const char* const zp = strstr( variable, "_Z" );
    const int isVectorX = AND2( xp, xp[ 2 ] == '\0' );
    const int isVectorY = AND2( yp, yp[ 2 ] == '\0' );
    const int isVectorZ = AND2( zp, zp[ 2 ] == '\0' );
    const int isVector = OR3( isVectorX, isVectorY, isVectorZ );
    int rank = 0;
    int dimensions[ 32 ];
    enum { NAME_LENGTH = 255 };
    char variable2[ NAME_LENGTH + 1 ] = "";
    const size_t nameLength0 = strlen( variable ) - 2 * isVector; /* Trim _X */
    const size_t nameLength =
      nameLength0 < NAME_LENGTH ? nameLength0 : NAME_LENGTH;
    memset( dimensions, 0, sizeof dimensions );
    memset( variable2, 0, sizeof variable2 );
    strncpy( variable2, variable, nameLength );
    result = readVariableDimensions( file, variable2, &rank, dimensions );

    if ( result ) {
      const int matchedDimensions =
        AND3( rank == 2, dimensions[ 0 ] == points, dimensions[ 1 ] == levels );
      const size_t count =
        dimensions[ 0 ] * dimensions[ 1 ] * ( rank == 3 ? dimensions[ 2 ] : 1 );
      const size_t bytes = count * sizeof (double);
      double* buffer = matchedDimensions ? data : malloc( bytes );

      DEBUG( fprintf( stderr,
                      "data = %p, buffer = %p, bytes = %lu, rank = %d\n",
                      data, buffer, bytes, rank ); )

      if ( ! buffer ) {
        fprintf( stderr, "\a\n\nFailed to allocate %lu bytes "
                 "to complete requested operation.\n", bytes );
      } else {
        result =
          readFileData( file, variable2, rank, dimensions, units, buffer );

        if ( result ) {

          if ( rank == 3 ) {

            if ( dimensions[ 2 ] > 1 ) {

             if ( isVector ) {
                const int length = strstr( variable2, "Surface_Wind" ) ? 2 : 3;
                const int vectorComponent = isVectorY + isVectorZ + isVectorZ;
                copyVectorComponent( dimensions[ 0 ], length, vectorComponent,
                                     buffer, data );
              } else if ( ! strcmp( variable2, "CAD_Score" ) ) {
                copyWorstCADScore( count, buffer, data );
              } else if ( ! strcmp( variable2,
                                    "Atmospheric_Volume_Description" ) ) {
                /* Pick 2nd value of pair (classification of lower bin): */
                copyVectorComponent( dimensions[ 0 ] * dimensions[ 1 ],
                                     dimensions[ 2 ], 1, buffer, data );
              } else {
                copyMaximumComponent( dimensions[ 0 ] * dimensions[ 1 ],
                                      dimensions[ 2 ], buffer, data );
              }
            }
          } else {
            assert( rank == 2 );

            if ( dimensions[ 1 ] > 1 ) {

              if ( ! strcmp( variable2, "Profile_ID" ) ) {
                copyVectorComponent( dimensions[ 0 ], dimensions[ 1 ], 0,
                                     buffer, data );
              } else if ( ! strcmp( variable2, "Lidar_Surface_Elevation" ) ) {
                copyMeanComponents( dimensions[ 0 ], dimensions[ 1 ],
                                   buffer, data ); /* (top + bottom) / 2. */
              } else if ( OR4( ! strcmp( variable2, "Profile_UTC_Time" ),
                               ! strcmp( variable2, "Profile_Time" ),
                               ! strcmp( variable2, "Latitude" ),
                               ! strcmp( variable2, "Longitude" ) ) ) {
                copyVectorComponent( dimensions[ 0 ], dimensions[ 1 ], 1,
                                     buffer, data );
              } else if ( OR2( ! strcmp( variable2,
                                         "Surface_Elevation_Statistics"),
                               ! strcmp( variable2,
                                         "DEM_Surface_Elevation" ) ) ) {
                const int meanComponent = dimensions[ 1 ] / 2;
                copyVectorComponent( dimensions[ 0 ], dimensions[ 1 ],
                                     meanComponent, buffer, data );
              } else if ( isVector ) {
                const int length = strstr( variable2, "Surface_Wind" ) ? 2 : 3;
                const int vectorComponent = isVectorY + isVectorZ + isVectorZ;
                copyVectorComponent( dimensions[ 0 ], length, vectorComponent,
                                     buffer, data );
              }
            }
          }

          if ( levels > 1 ) {
            reverseLevels( points, levels, data );
          }
        }

        if ( buffer != data ) {
          free( buffer ), buffer = 0;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readCALIPSOTimestamps - Read CALIPSO profile timestamps.
INPUTS:  const int file               ID of file to read.
         const size_t points          Number of ground points of data.
OUTPUTS: double timestamps[ points ]  Timestamps of data points.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int readCALIPSOTimestamps( const int file, const size_t points,
                                  double timestamps[] ) {

  const char* const variable = "Profile_UTC_Time";
  int result = 0;
  int rank = 0;
  int dimensions[ 32 ];
  memset( dimensions, 0, sizeof dimensions );

  assert( file >= 0 );
  assert( points );
  assert( timestamps );

  result = readVariableDimensions( file, variable, &rank, dimensions );

  if ( result ) {
    result = 0;

    if ( ! AND2( rank == 2, dimensions[ 0 ] == points ) ) {
      fprintf( stderr, "\a\n\nInvalid dimensions of variable %s.\n", variable );
    } else {
      const int components = dimensions[ 1 ];
      char unused[ 80 ] = "";

      if ( components == 1 ) {
        result = readFileData(file, variable, 2, dimensions, unused, timestamps);
      } else if ( components == 3 ) {
        const size_t bytes = points * components * sizeof (double);
        double* buffer = malloc( bytes );

        if ( ! buffer ) {
          fprintf( stderr, "\a\n\nFailed to allocate %lu bytes "
                  "to complete requested operation.\n", bytes );
        } else {
          dimensions[ 1 ] = components;
          result = readFileData( file, variable, 2, dimensions, unused, buffer);

          if ( result ) {
            copyVectorComponent( points, components, 1, buffer, timestamps );
          }
        }
      } else {
        fprintf( stderr, "\a\n\nInvalid dimensions for variable %s: "
                 "components = %d (expected 1 or 3).\n", variable, components );
      }

      if ( result ) {
        const double y2k = 20000000.0; /* Convert yymmdd.f to yyyymmdd.f */
        offsetValues(  y2k, points, timestamps );
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readCALIPSOCoordinates - Read CALIPSO longitudes and latitudes.
INPUTS:  const int file             ID of file to read.
         const int fileType         CALIPSO_L1 ... CALIPSO_L2_VFM.
         const size_t points        Number of ground points of data.
OUTPUTS: double longitudes[  points ]  Longitudes of data points.
         double latitudes[   points ]  Latitudes  of data points.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/


static int readCALIPSOCoordinates( const int file, const int fileType,
                                   const size_t points,
                                   double longitudes[], double latitudes[] ) {
  int result = 0;
  int rank = 0;
  int dimensions[ 32 ] = { 0, 0, 0 };

  assert( file >= 0 ); assert( IS_CALIPSO( fileType ) );
  assert( points ); assert( longitudes ); assert( latitudes );

  result = readVariableDimensions( file, "Longitude", &rank, dimensions );

  if ( AND2( rank == 2, dimensions[ 0 ] == points ) ) {

    if ( dimensions[ 1 ] == 1 ) { /* One value per ground point: */
      char unused[ 80 ] = "";
      result =
        readFileData( file, "Longitude", 2, dimensions, unused, longitudes );

      if ( result ) {
        result =
          readFileData( file, "Latitude", 2, dimensions, unused, latitudes );
      }

    } else if ( dimensions[ 1 ] == 3 ) { /* Copy middle of 3 coordinates: */
      const size_t bytes = points * 3 * sizeof (double);
      double* buffer = malloc( bytes );

      if ( ! buffer ) {
        fprintf( stderr, "\a\n\nFailed to allocate %lu bytes "
                 "to complete requested operation.\n", bytes );
      } else {
        char unused[ 80 ] = "";
        result =
          readFileData( file, "Longitude", 2, dimensions, unused, buffer );

        if ( result ) {
          copyVectorComponent( points, 3, 1, buffer, longitudes );
          result =
            readFileData( file, "Latitude", 2, dimensions, unused, buffer );

          if ( result ) {
            copyVectorComponent( points, 3, 1, buffer, latitudes );
          }
        }
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\a\n\nFailed to read valid coordinates.\n" );
  }

  return result;
}



/******************************************************************************
PURPOSE: readCALIPSOElevations - Read CALIPSO elevations and thicknesses.
INPUTS:  const int file             ID of file to read.
         const int fileType         CALIPSO_L1 ... CALIPSO_L2_VFM.
         const size_t points        Number of ground points of data.
         const size_t levels        Number of vertical points of data.
OUTPUTS: double elevations[  points * levels ]  Elevations (mAMSL) of data.
         double thicknesses[ points * levels ]  0 or layer thicknesses of data.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int readCALIPSOElevations( const int file, const int fileType,
                                  const size_t points, const size_t levels,
                                  double elevations[], double thicknesses[] ) {

  const double kilometersToMeters = 1000.0;
  int foundElevations = 0;
  int result = 0;

  assert( file >= 0 ); assert( IS_CALIPSO( fileType ) );
  assert( IMPLIES( AND2( IS_LAYERED( fileType ), levels > 1 ), thicknesses ) );
  assert( points ); assert( levels ); assert( elevations );

  DEBUG(fprintf(stderr, "readCALIPSOElevtions( %lu x %lu ):\n", points,levels));

  if ( AND2( levels > 1, IS_LAYERED( fileType ) ) ) {
    int dimensions[ 2 ] = { 0, 0 };
    char unused[ 80 ] = "";
    dimensions[ 0 ] = points;
    dimensions[ 1 ] = levels;
    foundElevations = 1;
    result =
      readFileData( file, "Layer_Base_Altitude", 2, dimensions, unused,
                    elevations );

    if ( result ) {
      result =
        readFileData( file, "Layer_Top_Altitude", 2, dimensions, unused,
                      thicknesses );

      if ( result ) {
        size_t point = 0;

#pragma omp parallel for

        for ( point = 0; point < points; ++point ) {
          const size_t offset = point * levels;
          size_t level = 0;
          size_t level2 = levels;

          /* Reorder level data as surface-to-sky: */

          for ( level = 0; level < level2--; ++level ) {
            const size_t index1 = offset + level;
            const size_t index2 = offset + level2;
            const double swapElevationTemp = elevations[  index1 ];
            const double swapThicknessTemp = thicknesses[ index1 ];
            elevations[  index1 ] = elevations[ index2 ];
            elevations[  index2 ] = swapElevationTemp;
            thicknesses[ index1 ] = thicknesses[ index2 ];
            thicknesses[ index2 ] = swapThicknessTemp;
          }

          /* Compute elevations and thicknesses: */

          for ( level = 0; level < levels; ++level ) {
            const size_t index = offset + level;
            const double kmBottom = elevations[  index ];
            const double kmTop    = thicknesses[ index ];

            if ( OR2( kmBottom < 0.0, kmTop <= kmBottom ) ) {
              elevations[  index ] = 0.0;
              thicknesses[ index ] = 0.0;
            } else {
              const double theBottom = kmBottom * kilometersToMeters;
              const double theTop    = kmTop    * kilometersToMeters;
              const double thickness = theTop - theBottom;
              const double middle    = ( theTop + theBottom ) * 0.5;
              elevations[  index ] = middle;
              thicknesses[ index ] = thickness;
            }
          }
        }
      }
    }
  } else {

    /*
     * Search for possible name of variable containing surface elevation in km.
     * The first dimension must match points.
     * The second dimension may be 1 or greater (vector).
     * If vector of values, compute the mean of 2-component vectors
     * else select component #2 (0, 1, 2, ...) as the mean.
     */

    const char* const elevationVariables[] = { /* Ordered by preference: */
      "Surface_Elevation",
      "Surface_Elevation_Statistics",
      "DEM_Surface_Elevation",
      "Lidar_Surface_Elevation"
    };
    const size_t elevationVariableCount =
      sizeof elevationVariables / sizeof *elevationVariables;
    int found = 0;
    char unused[ 80 ] = "";
    size_t variableIndex = 0;

    for ( variableIndex = 0;
          AND2( ! found, variableIndex < elevationVariableCount );
          ++variableIndex ) {
      const char* const elevationVariable = elevationVariables[ variableIndex ];
      found = fileVariableExists( file, elevationVariable );

      if ( found ) {
        int rank = 0;
        int dimensions[ 32 ];
        memset( dimensions, 0, sizeof dimensions );
        foundElevations = 1;
        result =
          readVariableDimensions(file, elevationVariable, &rank, dimensions);

        if ( AND4( result, rank == 2, dimensions[ 0 ] == points,
                   dimensions[ 1 ] >= 1 ) ) {
          const int components = dimensions[ 1 ];

          if ( components == 1 ) {
            result =
              readFileData( file, elevationVariable, 2, dimensions, unused,
                            elevations );
          } else {
            const size_t bytes = points * components * sizeof (double);
            double* buffer = malloc( bytes );

            if ( ! buffer ) {
              fprintf( stderr, "\a\n\nFailed to allocate %lu bytes "
                       "to complete requested operation.\n", bytes );
            } else {
              result =
                readFileData( file, elevationVariable, 2, dimensions, unused,
                              buffer );

              if ( result ) {

                if ( components == 2 ) {
                  copyMeanComponents( points, components, buffer, elevations );
                } else {
                  const int meanComponent = 2;
                  assert( components > meanComponent );
                  copyVectorComponent( points, components, meanComponent,
                                       buffer, elevations );
                }
              }

              free( buffer ), buffer = 0;
            }
          }
        }
      }
    }

    /* VFM files lack elevations so just set elevations to sea-level (HACK): */

    if ( AND2( ! foundElevations, fileType == CALIPSO_L2_VFM ) ) {
      memset( elevations, 0, points * levels * sizeof *elevations );
      result = 1;
    }

    if ( result ) {
      scaleValues( kilometersToMeters, points, elevations );

      if ( levels > 1 ) { /* Must read altitudes and convert/replicate: */
        const size_t bytes = levels * sizeof (double);
        double* altitudes = malloc( bytes );
        result = altitudes != 0;

        if ( ! altitudes ) {
          fprintf( stderr, "\a\n\nFailed to allocate %lu bytes "
                  "to complete requested operation.\n", bytes );
        } else {
          const char* const altitudesVariable =
            levels == 33 ? "Met_Data_Altitudes" : "Lidar_Data_Altitudes";
          result = readFileVData( file, altitudesVariable, levels, altitudes );

          /*
           * Given surface elevations (m): elevations[ points ]
           * and altitudes[ levels ] (km, sky-to-surface order),
           * compute elevations[ points ][ levels ] (m, surface-to-sky order):
           */

          if ( result ) {
            size_t point = points;

            while ( point-- ) {
              const double surfaceElevation = elevations[ point ];
              size_t level = levels;

              while ( level-- ) {
                const double altitude = altitudes[ level ] * kilometersToMeters;
                const double elevation = altitude;
                const size_t index = point * levels + ( levels - level - 1 );

                if ( elevation < surfaceElevation ) {
                  elevations[ index ] = surfaceElevation;
                } else {
                  elevations[ index ] = elevation;
                }
              }
            }
          }

          free( altitudes ), altitudes = 0;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: filterCALIPSOData - Filter CALIPSO data by applicable QC variables and
         hard-coded rules.
INPUTS:  const int file             ID of file to read.
         const int fileType         CALIPSO_L1 ... CALIPSO_L2_VFM.
         const char* const variable Name of data variable to filter.
         const size_t points        Number of ground points of data.
         const size_t levels        Number of vertical points of data.
         const double minimumCAD    Minimum cloud/aerosol discrimination score.
         E.g., 20 accepts score in range [20, 100].
         const double maximumUncertainty  Maximum uncertainty acceptable.
         E.g., 99 accepts absolute uncertainty <= 99.
         Units are the same as the data variable.
         const double elevations[ points * levels ] Elevations (mAMSL) of data.
         double data[ points * levels ]             Variable data to filter.
OUTPUTS: double data[ points * levels ]             Filtered variable data.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int filterCALIPSOData( const int file, const int fileType,
                              const char* const variable,
                              const size_t points, const size_t levels,
                              double minimumCAD, double maximumUncertainty,
                              const double elevations[], double data[] ) {

  typedef struct {
    int type;                 /* CALIPSO_L1, CALIPSO_L2_05KMAPRO, etc. */
    const char* variableName; /* Name of data variable to filter. */
    double dataMinimum;       /* Minimum valid value for this variable. */
    double dataMaximum;       /* Maximum valid value for this variable. */
    double missingValue;      /* Missing value for this variable. */
    const char* qcVariable;   /* Name of data variable to filter by. */
    size_t qcLevels;          /* if 0 then use variable levels. */
    double qcMinimum;         /* Minimum QC value to accept. */
    double qcMaximum;         /* Maximum QC value to accept. */
    unsigned int mask;        /* If non-0 then & with unsigned int QC and */
                              /* accept if 0, else apply qcMinimum,qcMaximum.*/
    int propagateDown;        /* Propagate the bad/worst qc value from the */
                              /* above layers downward before filtering? */
  } FilterParameters;

  /*
   * Variables named in this table get filtered by their named QC variables
   * according to recommendations of the NASA Langley CALIPSO Team 2012-02-09.
   */

  static const FilterParameters filterTable[] = {

    /* CALIPSO_L1: */

    {
      CALIPSO_L1,
      "Total_Attenuated_Backscatter_532", -0.075, 2.5, MISSING_VALUE,
      "QC_Flag", 1, 0.0, 0.0, 0xfffffc3f, 0 /* Clear bits 6-9. */
    },
#if 0
    /* Applying QC_Flag_2 yields almost no data! So disable this filter. */
    {
      CALIPSO_L1,
      "Total_Attenuated_Backscatter_532", -0.075, 2.5, MISSING_VALUE,
      "QC_Flag_2", 1, 0.0, 0.0, 0, 0
    },
#endif
    {
      CALIPSO_L1,
      "Perpendicular_Attenuated_Backscatter_532", -0.075, 1.5, MISSING_VALUE,
      "QC_Flag", 1, 0.0, 0.0, 0xfffffc3f, 0 /* Clear bits 6-9. */
    },
#if 0
    /* Applying QC_Flag_2 yields almost no data! So disable this filter. */
    {
      CALIPSO_L1,
      "Perpendicular_Attenuated_Backscatter_532", -0.075, 1.5, MISSING_VALUE,
      "QC_Flag_2", 1, 0.0, 0.0, 0, 0
    },
#endif
    {
      CALIPSO_L1,
      "Attenuated_Backscatter_1064", -0.075, 2.5, MISSING_VALUE,
      "QC_Flag", 1, 0.0, 0.0, 0xfffffc3f, 0 /* Clear bits 6-9. */
    },
#if 0
    /* Applying QC_Flag_2 yields almost no data! So disable this filter. */
    {
      CALIPSO_L1,
      "Attenuated_Backscatter_1064", -0.075, 2.5, MISSING_VALUE,
      "QC_Flag_2", 1, 0.0, 0.0, 0, 0
    },
#endif
    {
      CALIPSO_L1,
      "Depolarization_Gain_Ratio_532", 0.0, 2.0, MISSING_VALUE,
      0, 0, 0.0, 0.0, 0, 0
    },

    /* CALIPSO_L2_05KMAPRO: */

    {
      CALIPSO_L2_05KMAPRO,
      "Extinction_Coefficient_532", -0.2, 2.5, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Extinction_Coefficient_532", -0.2, 2.5, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Extinction_Coefficient_1064", -0.2, 2.5, MISSING_VALUE,
      "Extinction_QC_Flag_1064", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4.*/
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Extinction_Coefficient_1064", -0.2, 2.5, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Total_Backscatter_Coefficient_532", -0.01, 0.125, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Total_Backscatter_Coefficient_532", -0.01, 0.125, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Perpendicular_Backscatter_Coefficient_532", -0.01, 0.025, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Perpendicular_Backscatter_Coefficient_532", -0.01, 0.025, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Backscatter_Coefficient_1064", -0.01, 0.075, MISSING_VALUE,
      "Extinction_QC_Flag_1064", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4.*/
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Backscatter_Coefficient_1064", -0.01, 0.075,
      MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Particulate_Depolarization_Ratio_Profile_532", -0.05, 0.8,MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4.*/
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Particulate_Depolarization_Ratio_Profile_532", -0.05, 0.8,MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Aerosol_Layer_Fraction", 0.0, 30.0, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO, /* Filtered-out use CALIPSO_L2_05KMCPRO instead. */
      "Cloud_Layer_Fraction", 0.0, 30.0, MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 0
    },
#if 0
    /*
     * If *_QC_Flag_* and CAD_Score are applied (propagating worst value down)
     * then all Column_* data is filtered-out! So disable this filter.
     */
    {
      CALIPSO_L2_05KMAPRO,
      "Column_Optical_Depth_Aerosols_532", 0.0, 5.0, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 1 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Column_Optical_Depth_Aerosols_532", 0.0, 5.0, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 1
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Column_Optical_Depth_Aerosols_1064", 0.0, 5.0, MISSING_VALUE,
      "Extinction_QC_Flag_1064", 0, 0.0, 0.0, 0xffffffec, 1 /* Clear 0,1,4.*/
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Column_Optical_Depth_Aerosols_1064", 0.0, 5.0, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 1
    },
    {
      CALIPSO_L2_05KMAPRO, /* Filtered-out use CALIPSO_L2_05KMCPRO instead. */
      "Column_Optical_Depth_Cloud_532", 0.0, 5.0, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 1 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMAPRO, /* Filtered-out use CALIPSO_L2_05KMCPRO instead. */
      "Column_Optical_Depth_Cloud_532", 0.0, 5.0, MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 1
    },
#else
    /*
     * Instead, just filter by data range (then Uncertainty)
     * so we can get Column_* data...
     */
    {
      CALIPSO_L2_05KMAPRO,
      "Column_Optical_Depth_Aerosols_532", 0.0, 5.0, MISSING_VALUE,
      0, 0, 0.0, 0.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO,
      "Column_Optical_Depth_Aerosols_1064", 0.0, 5.0, MISSING_VALUE,
      0, 0, 0.0, 0.0, 0, 0
    },
    {
      CALIPSO_L2_05KMAPRO, /* Filtered-out use CALIPSO_L2_05KMCPRO instead. */
      "Column_Optical_Depth_Cloud_532", 0.0, 5.0, MISSING_VALUE,
      0, 0, 0.0, 0.0, 0, 0
    },
#endif

    /* CALIPSO_L2_05KMCPRO: */

    {
      CALIPSO_L2_05KMCPRO,
      "Extinction_Coefficient_532", -0.2, 2.5, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Extinction_Coefficient_532", -0.2, 2.5, MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 0
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Total_Backscatter_Coefficient_532", -0.01, 0.125, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Total_Backscatter_Coefficient_532", -0.01, 0.125, MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 0
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Perpendicular_Backscatter_Coefficient_532", -0.01, 0.025, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Perpendicular_Backscatter_Coefficient_532", -0.01, 0.025, MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 0
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Particulate_Depolarization_Ratio_Profile_532", 0.0, 1.0, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 0 /* Clear 0,1,4.*/
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Particulate_Depolarization_Ratio_Profile_532", -0.05, 0.8,MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 0
    },
    {
      CALIPSO_L2_05KMCPRO, /* Filtered-out use CALIPSO_L2_05KMAPRO instead. */
      "Aerosol_Layer_Fraction", 0.0, 30.0, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 0
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Cloud_Layer_Fraction", 0.0, 30.0, MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 0
    },
#if 0
    /*
     * If *_QC_Flag_* and CAD_Score are applied (propagating worst value down)
     * then all Column_* data is filtered-out! So disable this filter.
     */
    {
      CALIPSO_L2_05KMCPRO,
      "Column_Optical_Depth_Cloud_532", 0.0, 5.0, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 1 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMCPRO,
      "Column_Optical_Depth_Cloud_532", 0.0, 5.0, MISSING_VALUE,
      "CAD_Score", 0, 20.0, 100.0, 0, 1
    },
    {
      CALIPSO_L2_05KMCPRO, /* Filtered-out use CALIPSO_L2_05KMAPRO instead. */
      "Column_Optical_Depth_Aerosols_532", 0.0, 5.0, MISSING_VALUE,
      "Extinction_QC_Flag_532", 0, 0.0, 0.0, 0xffffffec, 1 /* Clear 0,1,4. */
    },
    {
      CALIPSO_L2_05KMCPRO, /* Filtered-out use CALIPSO_L2_05KMAPRO instead. */
      "Column_Optical_Depth_Aerosols_532", 0.0, 5.0, MISSING_VALUE,
      "CAD_Score", 0, -100.0, -20.0, 0, 1
    },
#else
    /*
     * Instead, just filter by data range (then Uncertainty)
     * so we can get Column_* data...
     */
    {
      CALIPSO_L2_05KMCPRO,
      "Column_Optical_Depth_Cloud_532", 0.0, 5.0, MISSING_VALUE,
      0, 0, 0.0, 0.0, 0, 0
    },
    {
      CALIPSO_L2_05KMCPRO, /* Filtered-out use CALIPSO_L2_05KMAPRO instead. */
      "Column_Optical_Depth_Aerosols_532", 0.0, 5.0, MISSING_VALUE,
      0, 0, 0.0, 0.0, 0, 0
    },
#endif

    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } /* End of table. */
  };

  int result = 1;
  int index = 0;

  assert( file >= 0 ); assert( IS_CALIPSO( fileType ) );
  assert( variable ); assert( *variable );
  assert( points ); assert( levels ); assert( elevations ); assert( data );

  /* Search table for variable name and apply each filter variable found: */

  for ( index = 0; AND2( result, filterTable[ index ].variableName); ++index ) {
    const FilterParameters* const entry = filterTable + index;
    const int found =
      AND2( fileType == entry->type, ! strcmp( variable, entry->variableName ));

    if ( found ) {
      const int qcLevels = entry->qcLevels ? entry->qcLevels : levels;
      double qcMinimum = entry->qcMinimum;
      double qcMaximum = entry->qcMaximum;
      const int isCAD =
        AND2( entry->qcVariable, ! strcmp( entry->qcVariable, "CAD_Score" ) );

      if ( isCAD ) {

        /* Adjust QC range by minimumCADScore w/ appropriate sign: */

        if ( qcMinimum < 0.0 ) {
          qcMaximum = -minimumCAD;
        } else {
          qcMinimum = minimumCAD;
        }
      }

      result =
        filterDataByQC( file, fileType,
                        entry->dataMinimum, entry->dataMaximum,
                        entry->missingValue,
                        entry->qcVariable,
                        qcLevels, qcMinimum, qcMaximum, entry->mask,
                        entry->propagateDown, points, levels, data );
    }
  }

  if ( result ) {

    /*
     * Filter points with associated absolute uncertainty > maximumUncertainty.
     */

    const char* const uncertainty = uncertaintyVariable( variable );

    if ( *uncertainty ) {
      const int found = fileVariableExists( file, uncertainty );

      if ( found ) {
          result = filterDataByQC( file, fileType,
                                   -1000.0, 1000.0, MISSING_VALUE,
                                   uncertainty, levels, 0.0,
                                   maximumUncertainty, 0, 0,
                                   points, levels, data );
      }
    }

    if ( result ) {

      /*
       * 2012-02-09 Per NASA Langley CALIPSO Team:
       * L2 profile data within 180m of surface (first 3 layers)
       * is possibly invalid so filter it out!
       */

      if ( AND2( levels > 1,
                 ! OR2( fileType == CALIPSO_L1, IS_LAYERED( fileType ) ) ) ) {
        filterDataNearSurface( points, levels, elevations, data );
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: copyWorstCADScore - Copy an array of the worst of each pair of CAD
         scores.
INPUTS:  const size_t  count          Number of scalar values in array.
         const double input[ count ]  Array with sequence of pairs of values.
OUTPUTS: double output[ count ]       Array of worst scalar component values.
NOTES:   CAD scores are in the range [-100, 100] with closest to +/-100
         (i.e., furthest from 0) considered best. Some values may be -9999.0.
******************************************************************************/

static void copyWorstCADScore( const size_t count, const double input[],
                               double output[] ) {

  const size_t points = count / 2;
  size_t point = 0;
  assert( count >= 2 ); assert( count % 2 == 0 );
  assert( input ); assert( output );

  /* Copy the worst component forward: */

  for ( point = 0; point < points; ++point ) {
    const size_t index = point + point;
    const double score1 = input[ index ];
    const double score2 = input[ index + 1 ];
    const double absScore1 = score1 < 0.0 ? -score1 : score1;
    const double absScore2 = score2 < 0.0 ? -score2 : score2;
    double worst = score1;

    if ( absScore1 > 100.0 ) {

      if ( absScore2 > 100.0 ) {

        if ( absScore2 > absScore1 ) {
          worst = score2;
        }
      }
    } else if ( absScore2 > 100.0 ) {
      worst = score2;
    } else if ( absScore2 < absScore1 ) {
      worst = score2;
    }

    output[ point ] = worst;
  }
}



/******************************************************************************
PURPOSE: uncertaintyVariable - Name of uncertainty variable or 0 if none.
INPUTS:  const char* const variable  Variable name.
RETURNS: const char* uncertainty variable name, else "".
******************************************************************************/

static const char* uncertaintyVariable( const char* const variable ) {
  const char* const variables =
    " Depolarization_Gain_Ratio_532 "
    /*    " Calibration_Constant_532 " */
    /*  " Calibration_Constant_1064 " */
    " Column_Optical_Depth_Cloud_532 "
    " Column_Optical_Depth_Aerosols_532 "
    " Column_Optical_Depth_Aerosols_1064 "
    " Column_Optical_Depth_Stratospheric_532 "
    " Column_Optical_Depth_Stratospheric_1064 "
    " Total_Backscatter_Coefficient_532 "
    " Perpendicular_Backscatter_Coefficient_532 "
    " Paticulate_Depolarization_Ratio_Profile_532 "
    " Extinction_Coefficient_532 "
    " Backscatter_Coefficient_1064 "
    " Extinction_Coefficient_1064 "
    " Parallel_Column_Reflectance_532 "
    " Perpendiculal_Column_Reflectance_532 "
    " Integrated_Attenuated_Backscatter_532 "
    " Integrated_Attenuated_Backscatter_1064 "
    " Integrated_Volume_Depolarization_Ratio "
    " Integrated_Attenuated_Total_Color_Ratio "
    " Measured_Two_Way_Transmittance_532 "
    " Normalization_Constant_532 "
    " Feature_Optical_Depth_532 "
    " Feature_Optical_Depth_1064 "
    " Integrated_Particulate_Color_Ratio "
    " Integrated_Particulate_Depolarization_Ratio "
    " Cirrus_Shape_Parameter "
    " Ice_Water_Path "
  ;
  static char result[ 80 ] = "";
  assert( variable ); assert( *variable );
  memset( result, 0, sizeof result );

  if ( strstr( variables, variable ) ) {
    strncpy( result, variable, sizeof result / sizeof *result );
    {
      char* const tag1 = strstr( result, "_532" );
      char* const tag2 = tag1 ? tag1 : strstr( result, "_1064" );

      if ( tag2 ) {
        *tag2 = '\0';
      }

      strcat( result, "_Uncertainty" );

      if ( tag2 ) {
        const char* const tag3 = strstr( variable, "_532" );
        const char* const tag4 = tag3 ? tag3 : strstr( variable, "_1064" );
        assert( tag4 ); assert( *tag4 );
        strcat( result, tag4 );
      }
    }
  }

  assert( strlen( result ) < sizeof result / sizeof *result );
  return result;
}



/******************************************************************************
PURPOSE: filterDataNearSurface - Filter-out data in/below 180m of surface!
INPUTS:  const size_t points                         Number of ground points.
         const size_t levels                         Number of vertical levels.
         const double elevations[ points * levels ]  Elevations of data (mAMSL)
OUTPUTS: double data[ points * levels ]              Filtered data.
******************************************************************************/

static void filterDataNearSurface( const size_t points, const size_t levels,
                                   const double elevations[], double data[] ) {
  const double nearSurface = 180.0; /* meters. */
  size_t point = 0;

  assert( points > 0 ); assert( levels > 0 ); assert( elevations );
  assert( data );

  DEBUG( fprintf( stderr,
                 "filterDataNearSurface"
                 "(points = %lu, levels = %lu, elevations = %p, data = %p)\n",
                 points, levels, elevations, data ); )

  for ( point = 0; point < points; ++point ) {
    const size_t pointOffset = point * levels;
    const double surfaceElevation = elevations[ pointOffset ];
    size_t level = 0;

    for ( level = 0; level < levels; ++level ) {
      const size_t index = pointOffset + level;
      const double elevation = elevations[ index ];
      assert( IN_RANGE( index, 0, points * levels - 1 ) );

      if ( elevation <= surfaceElevation + nearSurface ) {
        data[ index ] = MISSING_VALUE;
      } else {
        level = levels; /* Above filter limit so stop looping. */
      }
    }
  }
}



/******************************************************************************
PURPOSE: filterDataByQC - Filter data based on data threshold and QC data.
INPUTS:  const int file                ID of HDF file.
         const int fileType            CALIPSO_L1 ... CALIPSO_L2_VFM.
         const double dataMinimum      Minimum valid data value.
         const double dataMaximum      Maximum valid data value.
         const double missingValue     Value for filtered data.
         const char* const qcVariable  Optional: Name of QC variable to read.
         const size_t qcLevels         Optional: Number of qc flag levels.
         const double qcMinimum        Optional: Minimum acceptable qc flag value.
         const double qcMaximum        Optional: Maximum acceptable qc flag value.
         const unsigned int mask       If non-0 then & with unsigned int QC and
                                       accept if 0,
                                       else apply qcMinimum,qmMaximum.
         const int propagateDown       Propagate the bad/worst qc value from the
                                       above layers downward before filtering?
         const size_t points           Number of ground points.
         const size_t levels           Number of vertical levels.
OUTPUTS: double data[ points * levels ]  Filtered data.
RETURNS: int 1 if there are unfiltered points remaining, else 0.
******************************************************************************/

static int filterDataByQC( const int file, const int fileType,
                           const double dataMinimum, const double dataMaximum,
                           const double missingValue,
                           const char* const qcVariable, const size_t qcLevels,
                           const double qcMinimum, const double qcMaximum,
                           const unsigned int mask, const int propagateDown,
                           const size_t points, const size_t levels,
                           double data[] ) {

  int result = 1;

  assert( file > -1 );
  assert( dataMinimum <= dataMaximum );
  assert( IMPLIES( qcVariable, OR2( qcLevels == 1, qcLevels == levels ) ) );
  assert( IS_BOOL( propagateDown ) );
  assert( data ); assert( points > 0 ); assert( levels > 0 );

  DEBUG( fprintf( stderr,
                  "qcVariable = '%s', mask = %x, propagateDown = %d\n",
                  qcVariable ? qcVariable : "null",
                  mask, propagateDown); )

  {
    char unused[ 80 ] = "";
    int rank = 0;
    int dimensions[ 32 ] = { 0, 0, 0 };
    double* qcFlags = 0;

    if ( qcVariable ) {
      result = readVariableDimensions( file, qcVariable, &rank, dimensions );
      result = AND3( result, IN3( rank, 2, 3), dimensions[ 0 ] == points );

      if ( result ) {
        const size_t qcLevelsRead = dimensions[ 1 ];
        const size_t bytes = points * qcLevelsRead * sizeof (double);
        qcFlags = malloc( bytes );
        result = qcFlags != 0;

        if ( ! qcFlags ) {
          fprintf( stderr, "\a\n\nFailed to allocate %lu bytes "
                   "to complete requested operation.\n", bytes );
        } else {
          result =
            readCALIPSOVariable( file, fileType, qcVariable,
                                 points, qcLevelsRead, unused, qcFlags );

          if ( result ) {

            if ( propagateDown ) {

              if ( strstr( qcVariable, "_Uncertainty" ) ) {
                propagateBadUncertaintyDownward( points, qcLevelsRead, qcFlags);
              } else if ( ! strcmp( qcVariable, "CAD_Score" ) ) {
                propagateWorstCADScoreDownward( points, qcLevelsRead, qcFlags );
              } else if ( strstr( qcVariable, "QC_Flag" ) ) {
                propagateWorstQCFlagDownward( mask, points, qcLevelsRead,
                                              qcFlags );
              }
            }
          }
        }
      }
    }

    if ( result ) {
      const size_t count =
        filterData( dataMinimum, dataMaximum, missingValue,
                    points, levels, qcLevels,
                    qcFlags, qcMinimum, qcMaximum, mask, data );
      result = count > 0;
    }

    if ( qcFlags ) {
      free( qcFlags ), qcFlags = 0;
    }
  }

  assert( IS_BOOL( result ) );
  DEBUG( fprintf( stderr, "filterDataByQC() returning %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: propagateBadUncertaintyDownward - Propagate the bad uncertainty value
         if found, downward to all lower layers within the column.
INPUTS:  const size_t points                    Number of data ground points.
         const size_t levels                    Number of data elevation levels
         double uncertainty[ points * levels ]  Uncertainty values to check.
OUTPUTS: double uncertainty[ points * levels ]  Propagated uncertainty.
******************************************************************************/

static void propagateBadUncertaintyDownward( const size_t points,
                                             const size_t levels,
                                             double uncertainty[] ) {
  size_t point = 0;
  const double badUncertainty = 99.99;

  assert( points > 0 ); assert( levels > 0 ); assert( uncertainty );

  DEBUG( fprintf( stderr,
                  "propagateBadUncertaintyDownward"
                  "( points = %lu, levels = %lu, uncertainty = %p )\n",
                  points, levels, uncertainty ); )

  for ( point = 0; point < points; ++point ) {
    const size_t pointOffset = point * levels;
    size_t level = 0;
    int foundBadUncertainty = 0;

    for ( level = levels; level-- > 0; ) {
      const size_t index = pointOffset + level;
      const double value = uncertainty[ index ];
      assert( IN_RANGE( index, 0, points * levels - 1 ) );

      if ( value == badUncertainty ) {
        foundBadUncertainty = 1;
      }

      if ( foundBadUncertainty ) {
        uncertainty[ index ] = badUncertainty;
      }
    }
  }
}



/******************************************************************************
PURPOSE: propagateWorstCADScoreDownward - Propagate the worst CAD_Score value
         downward to all lower layers within the column.
INPUTS:  const size_t points              Number of data ground points.
         const size_t levels              Number of data elevation levels.
         double score[ points * levels ]  CAD score values to check.
OUTPUTS: double score[ points * levels ]  Propagated cad score values.
******************************************************************************/

static void propagateWorstCADScoreDownward( const size_t points,
                                            const size_t levels,
                                            double score[] ) {
  size_t point = 0;

  assert( points > 0 ); assert( levels > 0 ); assert( score );

  DEBUG( fprintf( stderr,
                  "propagateWorstCADScoreDownward"
                  "( points = %lu, levels = %lu, score = %p )\n",
                  points, levels, score ); )

  for ( point = 0; point < points; ++point ) {
    const size_t pointOffset = point * levels;
    int initialized = 0;
    double worst = 0.0;
    size_t level = 0;

    for ( level = levels; level-- > 0; ) {
      const size_t index = pointOffset + level;
      const double value = score[ index ];

      if ( ! initialized ) {
        worst = value;
        initialized = 1;
      } else {
        const double absScore = value < 0.0 ? -value : value;
        const double absWorst = worst < 0.0 ? -worst : worst;

        if ( absWorst > 100.0 ) {

          if ( absScore > 100.0 ) {

            if ( absScore > absWorst ) {
              worst = value;
            }
          }
        } else if ( absScore > 100.0 ) {
          worst = value;
        } else if ( absScore < absWorst ) {
          worst = value;
        }

        score[ index ] = worst;
      }
    }
  }
}



/******************************************************************************
PURPOSE: propagateWorstQCFlagDownward - Propagate the worst QC_Flag value
         downward to all lower layers within the column.
INPUTS:  const unsigned int mask       If non-0 then & with unsigned int QC and
                                       accept if 0,
                                       else apply qcMinimum,qcMaximum.
         const size_t points           Number of data ground points.
         const size_t levels           Number of data elevation levels.
         double qc[ points * levels ]  QC values to check.
OUTPUTS: double qc[ points * levels ]  QC values propgated.
******************************************************************************/

static void propagateWorstQCFlagDownward( const unsigned int mask,
                                          const size_t points,
                                          const size_t levels,
                                          double qc[] ) {
  size_t point = 0;

  assert( points > 0 ); assert( levels > 0 ); assert( qc );

  DEBUG( fprintf( stderr,
                  "propagateWorstQCFlagDownward"
                  "( mask = %x, points = %lu, levels = %lu, qc = %p )\n",
                  mask, points, levels, qc ); )

  for ( point = 0; point < points; ++point ) {
    const size_t pointOffset = point * levels;
    int initialized = 0;
    double worst = 0.0;
    size_t level = 0;

    for ( level = levels; level-- > 0; ) {
      const size_t index = pointOffset + level;
      const double value = qc[ index ];

      if ( ! initialized ) {
        worst = value; /* Lowest 4 bytes only. */
        initialized = 1;
      } else if ( mask ) {
        const unsigned int intQCValue = value;         /* Lowest 4 bytes only */
        const unsigned int ivalue = intQCValue & mask; /* Clear some bits. */
        const unsigned int intWorstValue = worst;      /* Lowest 4 bytes only */
        const unsigned int worstValue = intWorstValue & mask; /* Clear bits.*/

        if ( ivalue > worstValue ) {
          worst = value;
        } else {
          qc[ index ] = worst;
        }
      } else if ( value > worst ) {
        worst = value;
      } else {
        qc[ index ] = worst;
      }
    }
  }
}



/******************************************************************************
PURPOSE: filterData - Filter data based on data threshold and QC flag.
INPUTS: const double dataMinimum        Minimum valid data value.
        const double dataMaximum        Maximum valid data value.
        const double missingValue       Value for filtered data.
        const size_t points             Number of data ground points.
        const size_t levels             Number of data elevation levels.
        const size_t qcLevels           Optional: Number of qc flag levels.
        const double qcFlags[ points * qcLevels ]     Optional:  qcFlags.
        const double qcMinimum          Optional: Minimum acceptable qc value.
        const double qcMaximum          Optional: Maximum acceptable qc value.
        const unsigned int mask         If non-0 then & with unsigned int QC and
                                        accept if 0,
                                        else apply qcMinimum,qmMaximum.
         double data[ points * levels ] Data to filter.
OUTPUTS: double data[ points * levels ] Filtered data.
RETURNS: size_t number of unfiltered data points.
******************************************************************************/

static size_t filterData( const double dataMinimum, const double dataMaximum,
                          const double missingValue,
                          const size_t points, const size_t levels,
                          const size_t qcLevels,
                          const double qcFlags[],
                          const double qcMinimum, const double qcMaximum,
                          const unsigned int mask, double data[] ) {
  size_t result = 0;
  size_t point = 0;
  assert( data ); assert( points > 0 ); assert( levels > 0 );
  assert( dataMinimum <= dataMaximum );
  assert( IMPLIES( qcFlags, OR2( qcLevels == 1, qcLevels == levels ) ) );

  DEBUG( fprintf( stderr, "filterData(  dataMinimum = %lg, dataMaximum = %lg, "
                          "missingValue = %lg, "
                          "points = %lu, levels = %lu, qcLevels = %lu, "
                          "qcFlags = %p, qcMinimum = %lg, qcMaximum = %lg, "
                          "mask = %x, data = %p, )\n",
                  dataMinimum, dataMaximum, missingValue,
                  points, levels, qcLevels, qcFlags, qcMinimum, qcMaximum,
                  mask, data ); )

  for ( point = 0; point < points; ++point ) {
    const size_t pointOffset = point * levels;
    size_t level = 0;

    for ( level = 0; level < levels; ++level ) {
      const size_t dataIndex = pointOffset + level;
      const double dataValue = data[ dataIndex ];
      assert( IN_RANGE( dataIndex, 0, points * levels - 1 ) );

      if ( dataValue != missingValue ) {

        if ( ! IN_RANGE( dataValue, dataMinimum, dataMaximum ) ) {
          data[ dataIndex ] = missingValue;
        } else if ( qcFlags ) {
          const size_t qcLevel = qcLevels == 1 ? 0 : level;
          const size_t qcIndex = point * qcLevels + qcLevel;
          const double qcValue = qcFlags[ qcIndex ];
          assert( IN_RANGE( qcIndex, 0, points * qcLevels - 1 ) );

          if ( mask ) {
            const unsigned int intQCValue = qcValue; /* Lowest 4 bytes only. */
            const unsigned int value = intQCValue & mask; /* Clear some bits.*/

            if ( value != 0 ) {
              data[ dataIndex ] = missingValue;
            } else {
              ++result;
            }
          } else if ( ! IN_RANGE( qcValue, qcMinimum, qcMaximum ) ) {
            data[ dataIndex ] = missingValue;
          } else {
            ++result;
          }

        } else {
          ++result;
        }
      }
    }
  }

  DEBUG( fprintf( stderr, "unfiltered points = %lu\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: computeAggregateLevels - Compute number and stride of aggregated
         levels.
INPUTS:  const size_t levels        Number of levels.
         const size_t targetLevels  Desired number of aggregate levels.
OUTPUTS: size_t* const stride       Stride to yield result.
RETURNS: size_t number of aggregate levels.
******************************************************************************/

static size_t computeAggregateLevels( const size_t levels,
                                      const size_t targetLevels,
                                      size_t* const stride ) {
  size_t result = levels;
  assert( levels > 0 ); assert( targetLevels > 0 ); assert( stride );
  *stride = 1;

  if ( targetLevels < levels ) {
    size_t remainder = 0;
    *stride = computeStride( levels, targetLevels );
    remainder = ( levels % *stride ) != 0;
    result = levels / *stride + remainder;
  }

  assert( result > 0 ); assert( result <= levels );
  return result;
}



/******************************************************************************
PURPOSE: computeStride - Compute indexing stride that yields approximately
         target count.
INPUTS:  const size_t count   Number of points.
         const size_t target  Desired count. 0 = all points.
RETURNS: size_t Stride that yields approximately target points.
******************************************************************************/

static size_t computeStride( const size_t count, const size_t target ) {
  size_t result = 1;
  assert( count );

  if ( AND2( target, count > target ) ) {
    const double fcount = count;
    result = fcount / target + 0.5; /* Rounded. */
  }

  assert( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: aggregateDataAndElevations - Compute aggregated (mean) subsetData
         and subsetElevations at given point and level.
INPUTS:  const size_t points       Number of ground points.
         const size_t levels       Number of vertical levels.
         const size_t point        Current point (lower-left corner).
         const size_t level        Current level (bottom).
         const size_t width        # of ground points in aggregate rectangle.
         const size_t height       # of level  points in aggregate rectangle.
         const double data[ points * levels ]        Data to aggregate.
         const double elevations[ points * levels ]  Elevations to aggregate.
OUTPUTS: double* const meanDatum   Mean of all data in rectangle with
                                   lower-left corner at (point, level) and
                                   width stride and height levelStride.
         double* const meanElevation  Mean of all elevations in same rectangle.
******************************************************************************/

static void aggregateDataAndElevations( const size_t points,
                                        const size_t levels,
                                        const size_t point,
                                        const size_t level,
                                        const size_t width,
                                        const size_t height,
                                        const double data[],
                                        const double elevations[],
                                        double* const meanDatum,
                                        double* const meanElevation ) {

  const size_t endPoint = point + width;
  const size_t endLevel = level + height;
  size_t pointIndex = 0;
  size_t dataCount = 0;
  size_t elevationCount = 0;
  double elevationSum = 0.0;
  double dataSum = 0.0;

  assert( data ); assert( elevations );
  assert( points > 0 ); assert( levels > 0 );
  assert( point < points ); assert( level < levels );
  assert( width > 0 ); assert( width <= points );
  assert( point + width <= points );
  assert( height > 0 ); assert( height <= levels );
  assert( level + height <= levels );
  assert( meanDatum ); assert( meanElevation );

  for ( pointIndex = point; pointIndex < endPoint; ++pointIndex ) {
    const size_t pointOffset = pointIndex * levels;
    size_t levelIndex = 0;

    for ( levelIndex = level; levelIndex < endLevel; ++levelIndex ) {
      const size_t index = pointOffset + levelIndex;
      const double datum = data[ index ];
      const double elevation = elevations[ index ];
      elevationSum += elevation;
      ++elevationCount;

      DEBUG( if ( point == 0 && level == 0 )
               fprintf( stderr, "index = %lu, datum = %lf, elevation = %lf\n",
                        index, datum, elevation ); )

      if ( datum != MISSING_VALUE ) {
        dataSum += datum;
        ++dataCount;
      }
    }
  }

  assert( elevationCount );
  *meanElevation = elevationSum / elevationCount;

  if ( dataCount ) {
    *meanDatum = dataSum / dataCount;
  } else {
    *meanDatum = MISSING_VALUE;
  }

  DEBUG( if ( point == 0 && level == 0 )
           fprintf( stderr, "aggregate: [%lu .. %lu) [%lu .. %lu) "
                    "meanDatum (of %lu valid values) = %lf, "
                    "meanElevation (of %lu values) = %lf\n",
                    point, endPoint, level, endLevel,
                    dataCount, *meanDatum, elevationCount, *meanElevation ); )
}





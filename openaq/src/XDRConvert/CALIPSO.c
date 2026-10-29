
/******************************************************************************
PURPOSE: CALIPSO.c - Define routines for processing CALIPSO data.

NOTES:

HISTORY: 2007/12 plessel.todd@epa.gov, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h>  /* For stderr, fprintf(). */
#endif
#include <string.h> /* For memset().  */

#include <netcdf.h> /* For nc_close().  */

#include <Utilities.h> /* For PRE*(), NEW_ZERO(), Integer, Real, Stream. */

#include <Helpers.h>         /* For Name. */
#include <NetCDFUtilities.h> /* For createNetCDFFile(). */
#include <M3IO.h>            /* For writeM3IOHeader(). */
#include <Parameters.h>      /* For Parameters. */

/*================================== TYPES ==================================*/

enum { POINT, LEVEL, DIMENSIONS };

typedef struct {
  Integer      variables;     /* Profile_Time,Longitude,Latitude,Elevation,V.*/
  Integer      timesteps;     /* Number of timesteps. */
  Integer      profiles;      /* Number of profile fly-overs to process. */
  Integer      points;        /* Sum of all profile ground points. */
  Integer      maximumPoints; /* Number of points in largest profile. */
  Integer      levels;        /* Profile vertical levels. 1, 33, 583. */
  Real         domain[2][2];  /* domain[LONGITUDE LATITUDE][MINIMUM MAXIMUM] */
  UTCTimestamp timestamp;
  Line         note;          /* File note/description. */
  Name*        variable;      /* variable[ variables ]. E.g., "Temperature". */
  Name*        units;         /* units[ variables ]. E.g., "C". */
  Integer*     timestamps;    /* timestamps[ profiles ]. */
  Integer*     dimensions;    /* dimensions[ profiles ][ POINT LEVEL ]. */
  Real*        data;          /* variable[ maximumPoints ][ levels ]. */
  /* Regrid data: */
  Integer      totalRegriddedPoints;
  Integer*     outputPoints;  /* outputPoints[ timesteps ]. */
  Real*        longitudes;    /*longitudes[2*maxPoints]*/
  Real*        latitudes;     /*latitudes[ 2*maxPoints]*/
  Real*        elevations;    /*elevations[2*maxPoints*levels]*/
  Real*        thickness;     /*thickness[2*maxPoints*levels]*/
  Real*        gridLongitudes;/* gridLongitudes[ totalRegriddedPoints ]. */
  Real*        gridLatitudes; /* gridLatitudes[ totalRegriddedPoints ]. */
  Real*        gridElevations; /* Z in meters AMSL for gridData[] */
  Integer*     columns;       /* columns[ totalRegriddedPoints ]. */
  Integer*     rows;          /* rows[ totalRegriddedPoints ]. */
  Integer*     layers;        /* layers[ totalRegriddedPoints * gridLayers].*/
  Real*        gridData;      /* gridData[totalRegriddedPoints * gridLayers].*/
  Integer      gridLayers;    /* Number of grid layers. */
} CALIPSO;

typedef Integer (*Writer)( CALIPSO* calipso, const Parameters* parameters );

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateCALIPSO( CALIPSO* calipso );

static Integer isValidCALIPSO( const CALIPSO* calipso );

static Integer readXDRHeader( Stream* input, CALIPSO* calipso );

static Integer readXDRData( Stream* input, CALIPSO* calipso );

static Integer readRegriddedXDR( Stream* input, CALIPSO* calipso );

static Integer readRegriddedXDRData( Stream* input, const int version,
                                     CALIPSO* calipso );

static void computeLayers( CALIPSO* calipso );

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    CALIPSO* calipso );

static void countCALIPSOPoints( CALIPSO* calipso );

static Writer dispatcher( Integer format, Integer regrid );

static Integer writeASCII( CALIPSO* calipso, const Parameters* parameters );

static Integer writeASCIIData( CALIPSO* calipso, Stream* input,
                               Stream* output );

static Integer writeCOARDS( CALIPSO* calipso, const Parameters* parameters );

static Integer writeCOARDSHeader( Integer file, const CALIPSO* calipso );

static Integer writeCOARDSData( Stream* input, Integer file,
                                const CALIPSO* calipso );

static Integer writeRegriddedXDR( CALIPSO* calipso,
                                  const Parameters* parameters );

static Integer writeRegriddedASCII( CALIPSO* calipso,
                                    const Parameters* parameters );

static Integer writeRegriddedCOARDS( CALIPSO* calipso,
                                     const Parameters* parameters );

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const CALIPSO* calipso );

static Integer writeRegriddedCOARDSData( Integer file,
                                         CALIPSO* calipso,
                                         const Parameters* parameters );

static Integer writeRegriddedIOAPI( CALIPSO* calipso,
                                    const Parameters* parameters );

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const CALIPSO* calipso,
                                          const Grid* grid );

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const CALIPSO* calipso,
                                        const Grid* grid );

static void regridCALIPSO( Stream* input, Integer method, Grid* grid,
                           CALIPSO* calipso );

static Integer readProfileDataForTimestamp( Integer yyyydddhh00,
                                            Stream* input,
                                            CALIPSO* calipso,
                                            Integer* points );

static Integer skipProfileDataBeforeTimestamp( Integer yyyydddhh00,
                                               const CALIPSO* calipso,
                                               Stream* stream );

static Integer readProfileData( Stream* input,
                                Integer variables, Integer points,
                                Integer levels,
                                Real longitudes[], Real latitudes[],
                                Real elevations[], Real thickness[],
                                Real data[] );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateCALIPSO - Read input and write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  Parameters updated by translation.
******************************************************************************/

void translateCALIPSO( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok,
         parameters->input->ok( parameters->input ) );

  CALIPSO calipso;
  ZERO_OBJECT( &calipso );
  parameters->ok = 0;

  if ( readXDRHeader( parameters->input, &calipso ) ) {
    Writer writer = dispatcher( parameters->format, parameters->regrid );

    if ( ! writer ) {
      failureMessage( "Invalid/unsupported format/regrid specification." );
    } else if ( parameters->regrid ) {
      regridCALIPSO( parameters->input, parameters->regrid, parameters->grid,
                     &calipso );

      if ( calipso.totalRegriddedPoints == 0 ) {
        failureMessage( "No points projected onto the grid." );
      } else {

        if ( parameters->aggregationTimesteps ) {
          const Integer dataVariable = calipso.variables - 1;
          Integer totalOutputPoints = 0;
          const Integer aggregatedTimesteps =
            aggregateData( parameters->aggregationTimesteps,
                           0,
                           calipso.timesteps,
                           calipso.outputPoints,
                           calipso.gridLongitudes,
                           calipso.gridLatitudes,
                           calipso.elevations,
                           calipso.columns,
                           calipso.rows,
                           calipso.layers,
                           calipso.gridData,
                           0,
                           &totalOutputPoints );
          calipso.timesteps = aggregatedTimesteps;
          calipso.totalRegriddedPoints = totalOutputPoints;

          if ( AND2( parameters->aggregationTimesteps == 24,
                     ! OR2( strstr( calipso.variable[ dataVariable ], "daily" ),
                            strstr( calipso.variable[ dataVariable], "DAILY")))) {
            Name dailyName = "";
            memset( dailyName, 0, sizeof dailyName );
            snprintf( dailyName, sizeof dailyName / sizeof *dailyName,
                      "daily_%s",
                      calipso.variable[ dataVariable ] );
            strncpy( calipso.variable[ dataVariable ], dailyName,
                     sizeof dailyName / sizeof *dailyName );
          }
        }

        parameters->ok = writer( &calipso, parameters );
      }
    } else {
      parameters->ok = writer( &calipso, parameters );
    }
  }

  deallocateCALIPSO( &calipso );
  POST0( isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: compareRegriddedCALIPSO - Read REGRIDDED-CALIPSO input,
         compare it to CMAQ XDR data & write it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedCALIPSO( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok, parameters->input->ok( parameters->input ) );

  if ( ! AND2( parameters->compareFunction, parameters->data ) ) {
    failureMessage( "Invalid input for comparing." );
    parameters->ok = 0;
  } else {
    CALIPSO calipso;
    ZERO_OBJECT( &calipso );
    parameters->ok = 0;

    if ( readRegriddedXDR( parameters->input, &calipso ) ) {
      compareFunctionNameUnits( parameters->compareFunction,
                                parameters->convertFunction,
                                calipso.variable[ 0 ],
                                calipso.units[ 0 ],
                                parameters->variable,
                                parameters->units );

      if ( compareRegriddedXDR( parameters, &calipso ) ) {
        Writer writer = dispatcher( parameters->format, 1 );
        CHECK( writer );

        if ( calipso.totalRegriddedPoints == 0 ) {
          failureMessage( "No points projected onto the grid." );
        } else {
          parameters->ok = writer( &calipso, parameters );
        }
      }
    }

    deallocateCALIPSO( &calipso );
  }

  POST0( isValidParameters( parameters ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateCALIPSO - Deallocate contents of calipso structure.
INPUTS: CALIPSO* calipso Structure to deallocate contents of.
******************************************************************************/

static void deallocateCALIPSO( CALIPSO* calipso ) {
  PRE0( calipso );
  FREE( calipso->variable );
  FREE( calipso->timestamps );
  FREE( calipso->dimensions );
  FREE( calipso->data );
  ZERO_OBJECT( calipso );
}



/******************************************************************************
PURPOSE: isValidCALIPSO - Check calipso structure.
INPUTS: const CALIPSO* calipso Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidCALIPSO( const CALIPSO* calipso ) {
  const Integer result =
    AND10( calipso,
           calipso->note[ 0 ],
           isValidUTCTimestamp( calipso->timestamp ),
           calipso->variable,
           calipso->units,
           calipso->variable[ 0 ],
           calipso->variable[ calipso->variables - 1 ],
           calipso->units[ 0 ],
           calipso->units[ calipso->variables - 1 ],
           IMPLIES_ELSE( calipso->totalRegriddedPoints == 0,
             AND16( calipso->variables >= 5,
                    GT_ZERO3( calipso->profiles,
                              calipso->levels,
                              calipso->points ),
                    isValidLongitudeLatitude(
                                     calipso->domain[ LONGITUDE ][ MINIMUM ],
                                     calipso->domain[ LATITUDE  ][ MINIMUM ] ),
                    isValidLongitudeLatitude(
                                     calipso->domain[ LONGITUDE ][ MAXIMUM ],
                                     calipso->domain[ LATITUDE  ][ MAXIMUM ] ),
                    calipso->domain[ LONGITUDE ][ MINIMUM ] <=
                      calipso->domain[ LONGITUDE ][ MAXIMUM ],
                    calipso->domain[ LATITUDE  ][ MINIMUM ] <=
                      calipso->domain[ LATITUDE  ][ MAXIMUM ],
                    calipso->timestamps,
                    isValidTimestamp( calipso->timestamps[ 0 ] ),
                    isValidTimestamp(calipso->timestamps[calipso->profiles-1]),
                    calipso->timestamps[ calipso->profiles - 1 ] >=
                      calipso->timestamps[ 0 ],
                    calipso->dimensions,
                    calipso->dimensions[ LEVEL ] > 0,
                    calipso->dimensions[ POINT ] > 0,
                    calipso->dimensions[(calipso->profiles-1) * 2 + LEVEL] > 0,
                    calipso->dimensions[(calipso->profiles-1) * 2 + POINT] > 0,
                    calipso->data ),
             AND16( calipso->totalRegriddedPoints > 0,
                    calipso->timesteps > 0,
                    calipso->outputPoints,
                    minimumItemI( calipso->outputPoints,
                                  calipso->timesteps ) >= 0,
                    calipso->columns,
                    calipso->rows,
                    calipso->gridLayers > 0,
                    calipso->gridLongitudes,
                    calipso->gridLatitudes,
                    calipso->gridElevations,
                    calipso->gridData,
                    minimumItemI( calipso->columns,
                                  calipso->totalRegriddedPoints ) > 0,
                    minimumItemI( calipso->rows,
                                  calipso->totalRegriddedPoints ) > 0,
                    isNanFree( calipso->gridElevations,
                               calipso->totalRegriddedPoints *
                                 calipso->gridLayers ),
                    isNanFree( calipso->gridData,
                               calipso->totalRegriddedPoints *
                                 calipso->gridLayers ),
                    validLongitudesAndLatitudes(
                      calipso->totalRegriddedPoints,
                      calipso->gridLongitudes,
                      calipso->gridLatitudes ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRHeader - Read input and initialize calipso structure.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: CALIPSO* calipso Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDRHeader( Stream* input, CALIPSO* calipso ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         calipso, calipso->variable == 0, calipso->data == 0 );

  Integer result = 0;

  input->readString( input, calipso->note, COUNT( calipso->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( calipso->note );

    if ( readTimestamp( input, calipso->timestamp ) ) {
      Integer dimensions[ 3 ] = { 0, 0, 0 };

      if ( readDimensions( input, COUNT( dimensions ), dimensions ) ) {

        if ( dimensions[ 0 ] < 5 ) {
          failureMessage( "Invalid input data: variables must be 5." );
        } else {
          calipso->variables = dimensions[ 0 ];
          calipso->timesteps = dimensions[ 1 ];
          calipso->profiles  = dimensions[ 2 ];
          calipso->variable  = NEW_ZERO( Name, calipso->variables * 2 );

          if ( calipso->variable ) {
            calipso->units = calipso->variable + calipso->variables;

            if ( readVariablesAndUnits( input, calipso->variables,
                                        calipso->variable, calipso->units ) ) {

              if ( readDomain( input, calipso->domain ) ) {

                if ( skipInputLines( input, 4 ) ) {
                  result = readXDRData( input, calipso );
                }
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid CALIPSO data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidCALIPSO( calipso ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read initial binary data from input.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: CALIPSO* calipso Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDRData( Stream* input, CALIPSO* calipso ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         calipso, calipso->variables >= 5, calipso->profiles > 0,
         calipso->timestamps == 0, calipso->data == 0 );

  Integer result = 0;
  calipso->timestamps = NEW_ZERO( Integer, calipso->profiles );

  if ( calipso->timestamps ) {
    input->read64BitIntegers( input, calipso->timestamps, calipso->profiles );

    if ( AND3( input->ok( input ),
               isValidTimestamp( calipso->timestamps[ 0 ] ),
               isValidTimestamp( calipso->timestamps[calipso->profiles-1]))) {
      const Integer boundsCount = calipso->profiles * 2 * 2;
      Real* bounds = NEW_ZERO( Real, boundsCount );

      if ( bounds ) {
        input->read64BitReals(input, bounds, boundsCount );
        FREE( bounds ); /* Read and ignored bounds. */

        if ( input->ok( input ) ) {
          const Integer dimensionsCount = calipso->profiles * 2;
          calipso->dimensions = NEW_ZERO( Integer, dimensionsCount );

          if ( calipso->dimensions ) {
            input->read64BitIntegers( input, calipso->dimensions,
                                      dimensionsCount );

            /* Sum ground points per profile: */

            if ( input->ok( input ) ) {
              calipso->levels = calipso->dimensions[ LEVEL ];
              countCALIPSOPoints( calipso );

              if ( calipso->points > 0 ) {
                const Integer dataCount =
                  calipso->maximumPoints * calipso->levels;
                calipso->data = NEW_ZERO( Real, dataCount ); /* Largest only */
                result = isValidCALIPSO( calipso );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid CALIPSO data." );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, isValidCALIPSO( calipso ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDR - Read REGRIDDED-CALIPSO and initialize calipso.
INPUTS:  Stream* input     Input stream read from and parse.
OUTPUTS: CALIPSO* calipso  Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
NOTES:   Input data format is:

REGRIDDED-CALIPSO 2.0
http://eosweb.larc.nasa.gov/PRODOCS/calipso/table_calipso.html,CALIPSOSubset,\
 XDRConvert
2006-07-03T00:00:00-0000
# timesteps layers
24 24
# Variable name:
Total_Attenuated_Backscatter_532
# Variable units:
per_kilometer_per_steradian
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
 33 45 40 -97 6.37e+06 6.37e+06
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[25]:
 279 240 -1.008e+06 -1.62e+06 12000 12000 2 10000 1 0.995 0.99 0.98 0.97 \
 0.96 0.94 0.92 0.9 0.88 0.86 0.84 0.82 0.8 0.77 0.74 0.7 0.65 0.6 0.5 \
 0.4 0.3 0.2 0.1 0
# MSB 64-bit integers points[timesteps] and
# IEEE-754 64-bit reals longitudes[timesteps][points] and
# IEEE-754 64-bit reals latitudes[timesteps][points] and
# IEEE-754 64-bit reals elevations[timesteps][points][layers] and
# MSB 64-bit integers columns[timesteps][points] and
# MSB 64-bit integers rows[timesteps][points] and
# MSB 64-bit integers layers[timesteps][points][layers] and
# IEEE-754 64-bit reals data[timesteps][points][layers]:

******************************************************************************/

static Integer readRegriddedXDR( Stream* input, CALIPSO* calipso ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         calipso, calipso->variable == 0, calipso->data == 0,
         calipso->gridData == 0 );

  Integer result = 0;
  input->readString( input, calipso->note, COUNT( calipso->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( calipso->note );

    if ( readTimestamp( input, calipso->timestamp ) ) {
      Integer dimensions[ 2 ] = { 0, 0 };

      if ( readDimensions( input, 2, dimensions ) ) {
        calipso->timesteps  = dimensions[ 0 ];
        calipso->gridLayers = dimensions[ 1 ];
        calipso->timestamps = NEW_ZERO( Integer, calipso->timesteps );

        if ( calipso->timestamps ) {
          Integer timestamp = fromUTCTimestamp( calipso->timestamp );
          Integer timestep = 0;

          for ( timestep = 0; timestep < calipso->timesteps; ++timestep ) {
            calipso->timestamps[ timestep ] = timestamp;
            incrementTimestamp( &timestamp );
          }

          calipso->variables = 1;
          calipso->variable = NEW_ZERO( Name, 2 );

          if ( calipso->variable ) {
            calipso->units = calipso->variable + calipso->variables;

            if ( readVariablesAndUnits( input, calipso->variables,
                                        calipso->variable, calipso->units ) ) {
              char line[ 256 ] = "";
              int count = 7;
              int version = 1;
              memset( line, 0, sizeof line );
              input->readString( input, line, sizeof line / sizeof *line - 1 );

              if ( strcmp( line,
                           "# MSB 64-bit integers points[timesteps] and\n" )) {
                count += 4 + 1; /* Skip 4 line projection/grid + layers. */
                version = 2;
              }

              if ( skipInputLines( input, count - 1 ) ) {
                result = readRegriddedXDRData( input, version, calipso );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid CALIPSO data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidCALIPSO( calipso ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRData - Read regridded binary array data from input.
INPUTS:  Stream* input       Input stream read from and parse.
         const int version   1 or 2.
OUTPUTS: CALIPSO* calipso  Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Input data format is:

# MSB 64-bit integers points[timesteps] and
# IEEE-754 64-bit reals longitudes[timesteps][points] and
# IEEE-754 64-bit reals latitudes[timesteps][points] and
# IEEE-754 64-bit reals elevations[timesteps][points][layers] and
# MSB 64-bit integers columns[timesteps][points] and
# MSB 64-bit integers rows[timesteps][points] and
# MSB 64-bit integers layers[timesteps][points][layers] and
# IEEE-754 64-bit reals data[timesteps][points][layers]:

******************************************************************************/

static Integer readRegriddedXDRData( Stream* input, const int version,
                                     CALIPSO* calipso ) {

  PRE010( input, input->ok( input ), input->isReadable( input ),
          IN3( version, 1, 2 ),
          calipso, calipso->timesteps > 0, calipso->gridLayers > 0,
          calipso->variables == 1,
          calipso->profiles == 0, calipso->data == 0 );

  Integer result = 0;
  calipso->outputPoints = NEW_ZERO( Integer, calipso->timesteps );

  if ( calipso->outputPoints ) {
    input->read64BitIntegers( input, calipso->outputPoints,
                              calipso->timesteps );

    if ( input->ok( input ) ) {
      const Integer count = calipso->totalRegriddedPoints =
        sumI( calipso->outputPoints, calipso->timesteps );

      if ( count > 0 ) {
        const Integer countVertical = count * calipso->gridLayers;
        calipso->gridLongitudes =
          NEW_ZERO( Real, count * 4 + countVertical * 3 );

        if ( calipso->outputPoints ) {
          calipso->gridLatitudes  = calipso->gridLongitudes + count;
          calipso->gridElevations = calipso->gridLatitudes + count;
          calipso->columns        = (Integer*)
                                    calipso->gridElevations + countVertical;
          calipso->rows           = calipso->columns + count;
          calipso->layers         = calipso->rows + count;
          calipso->gridData       = (Real*) calipso->layers + countVertical;

          input->read64BitReals( input, calipso->gridLongitudes,
                                 count * 2 +
                                 ( version > 1 ) * countVertical );

          if ( input->ok( input ) ) {
            input->read64BitIntegers( input, calipso->columns,
                                      count * 2 +
                                      ( version > 1 ) * countVertical );

            if ( input->ok( input ) ) {
              input->read64BitReals( input, calipso->gridData, countVertical );

              if ( version == 1 ) {
                input->read64BitReals( input, calipso->gridElevations,
                                       countVertical );
              }

              if ( input->ok( input ) ) {

                if ( version == 1 ) {
                  computeLayers( calipso );
                }

                result = isValidCALIPSO( calipso );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid CALIPSO data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidCALIPSO( calipso ) ) );
  return result;
}



/******************************************************************************
PURPOSE: computeLayers - Compute layers for regridded data.
INPUTS:  CALIPSO* calipso  calipso->layers.
OUTPUTS: CALIPSO* calipso  Initialized calipso->layers.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static void computeLayers( CALIPSO* calipso ) {
  PRE02( calipso->totalRegriddedPoints, calipso->layers );
  const Integer points = calipso->totalRegriddedPoints;
  const Integer layers = calipso->gridLayers;
  Integer index = 0;
  Integer point = 0;

  for ( point = 0; point < points; ++point ) {
    Integer layer = 0;

    for ( layer = 0; layer < layers; ++layer ) {
      calipso->layers[ index++ ] = layer + 1;
    }
  }
}



/******************************************************************************
PURPOSE: compareRegriddedXDR - Compare Regridded data with CMAQ data.
INPUTS:  const Parameters* parameters  CMAQ data to compare to.
OUTPUTS: CALIPSO* calipso            Updated calipso->data.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    CALIPSO* calipso ) {

  PRE05( parameters, isValidParameters( parameters ),
         parameters->compareFunction,
         calipso, isValidCALIPSO( calipso ) );

  Integer result = 0;

  if ( ! AND2( ! strcmp( parameters->timestamp, calipso->timestamp ),
               parameters->timesteps == calipso->timesteps ) ) {
    failureMessage( "Mismatched time steps (%s %lld)"
                    " for comparison to CMAQ data (%s %lld).",
                    calipso->timestamp, calipso->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else {
    Real* const calipsoData             = calipso->gridData;
    const Integer* const calipsoRows    = calipso->rows;
    const Integer* const calipsoColumns = calipso->columns;
    const Integer* const calipsoPoints  = calipso->outputPoints;
    const Real* const cmaqData = parameters->data;
    CompareFunction comparer   = parameters->compareFunction;
    const Integer timesteps    = parameters->timesteps;
    const Integer firstLayer   = parameters->firstLayer;
    const Integer lastLayer    = parameters->lastLayer;
    const Integer firstRow     = parameters->firstRow;
    const Integer lastRow      = parameters->lastRow;
    const Integer firstColumn  = parameters->firstColumn;
    const Integer lastColumn   = parameters->lastColumn;
    const Integer layers       = lastLayer  - firstLayer  + 1;
    const Integer rows         = lastRow    - firstRow    + 1;
    const Integer columns      = lastColumn - firstColumn + 1;
    const Integer rowsTimesColumns = rows * columns;
    const Integer layersTimesRowsTimesColumns = layers * rowsTimesColumns;
    Integer timestep = 0;
    Integer calipsoDataIndex = 0;
    Integer calipsoPointIndex = 0;

    DEBUG( fprintf( stderr,
                    "grid L,R,C ranges: %lld..%lld %lld..%lld %lld..%lld\n",
                    firstLayer, lastLayer, firstColumn, lastColumn,
                    firstRow, lastRow ); )

    DEBUG( fprintf( stderr, "timesteps = %lld\n", timesteps ); )

    for ( timestep = 0; timestep < timesteps; ++timestep ) {
      const Integer points = calipsoPoints[ timestep ];
      const Integer timestepOffset = timestep * layersTimesRowsTimesColumns;
      Integer point = 0;

      DEBUG( fprintf( stderr, "timestep = %lld, points = %lld\n",
                      timestep, points ); )

      for ( point = 0; point < points; ++point, ++calipsoPointIndex ) {
        const Integer calipsoRow    = calipsoRows[    calipsoPointIndex ];
        const Integer calipsoColumn = calipsoColumns[ calipsoPointIndex ];
        Integer layer = 0;

        for ( layer = 0; layer < layers; ++layer, ++calipsoDataIndex ) {
          const Integer calipsoLayer = layer + 1;
          DEBUG( fprintf( stderr, "  @[ %lld, %lld, %lld ]: ",
                          calipsoLayer, calipsoRow, calipsoColumn ); )

          if ( AND3( IN_RANGE( calipsoLayer, firstLayer, lastLayer ),
                     IN_RANGE( calipsoRow, firstRow, lastRow ),
                     IN_RANGE( calipsoColumn, firstColumn, lastColumn ) ) ) {
            const Integer calipsoLayer0  = layer;
            const Integer calipsoRow0    = calipsoRow    - firstRow;
            const Integer calipsoColumn0 = calipsoColumn - firstColumn;
            const Integer dataIndex =
              timestepOffset + calipsoLayer0 * rowsTimesColumns +
              calipsoRow0 * columns + calipsoColumn0;
            CHECK4( IN_RANGE( calipsoLayer0, 0, layers - 1 ),
                    IN_RANGE( calipsoRow0, 0, rows - 1 ),
                    IN_RANGE( calipsoColumn0, 0, columns - 1 ),
                    IN_RANGE( dataIndex, 0,
                              timesteps * layers * rows * columns - 1 ) );
            const Real calipsoDatum = calipsoData[ calipsoDataIndex ];
            const Real cmaqDatum = cmaqData[ dataIndex ];
            const Real comparedDatum = comparer( calipsoDatum, cmaqDatum );
            calipsoData[ calipsoDataIndex ] = comparedDatum;
            result = 1;
            DEBUG( fprintf( stderr, "f(%lf, %lf) -> %lf\n",
                            calipsoDatum, cmaqDatum, comparedDatum ); )
          } else {
            calipsoData[ calipsoDataIndex ] = -9999.0;
            DEBUG( fprintf( stderr, "-9999\n" ); )
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "No points in output." );
  }

  POST02( IS_BOOL( result ), isValidCALIPSO( calipso ) );
  return result;
}



/******************************************************************************
PURPOSE: countCALIPSOPoints - Compute sum and max of profile ground points.
INPUTS:  CALIPSO* calipso  Partially initialized structure.
OUTPUTS: CALIPSO* calipso  calipso->maximumPoints is maximum profile points.
******************************************************************************/

static void countCALIPSOPoints( CALIPSO* calipso ) {

  PRE05( calipso, calipso->profiles > 0, calipso->dimensions,
         calipso->points == 0, calipso->maximumPoints == 0 );

  Integer profile = 0;

  do {
    const Integer profilePoints = calipso->dimensions[ profile * 2 + POINT ];

    if ( profilePoints > 0 ) {
      calipso->points += profilePoints;

      if ( profilePoints > calipso->maximumPoints ) {
        calipso->maximumPoints = profilePoints;
      }
    } else {
      calipso->points = calipso->maximumPoints = 0;
      profile = calipso->profiles;
    }

    ++profile;
  } while ( profile < calipso->profiles );


  POST02( calipso->points >= 0,
          IMPLIES_ELSE( calipso->points > 0,
                        calipso->maximumPoints > 0,
                        calipso->maximumPoints == 0 ) );
}



/******************************************************************************
PURPOSE: dispatcher - Look-up and return a writer for the given format/regrid.
INPUTS:  Integer format  E.g., FORMAT_XDR.
         Integer regrid  E.g., 0 or AGGREGATE_MEAN.
RETURNS: Writer a writer routine or 0 if none found.
******************************************************************************/

static Writer dispatcher( Integer format, Integer regrid ) {

  PRE02( IS_VALID_FORMAT( format ),
         OR2( regrid == 0, IS_VALID_AGGREGATE_METHOD( regrid ) ) );

  typedef struct {
    Integer format;         /* FORMAT_XDR, etc. */
    Writer writer;          /* Routine that writes data in this format. */
    Writer regriddedWriter; /* Routine that writes regridded data in format. */
  } Entry;

  static const Entry writers[] = {
    { FORMAT_XDR,    0,           writeRegriddedXDR    },
    { FORMAT_ASCII,  writeASCII,  writeRegriddedASCII  },
    { FORMAT_COARDS, writeCOARDS, writeRegriddedCOARDS },
    { FORMAT_IOAPI,  0,           writeRegriddedIOAPI  },
    { -1, 0, 0 }
  };

  Writer result = 0;
  const Integer count = COUNT( writers );
  Integer index = 0;

  do {
    const Entry* const entry = writers + index;

    if ( entry->format == -1 ) {
      index = count;
    } else if ( entry->format == format ) {
      result = regrid == 0 ? entry->writer : entry->regriddedWriter;
      index = count;
    }

    ++index;
  } while ( index < count );

  return result;
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format output.
INPUTS:  CALIPSO* calipso   Structure to write.
         const Parameters*  Parameters input stream to read from.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeASCII( CALIPSO* calipso, const Parameters* parameters ) {

  PRE03( isValidCALIPSO( calipso ),
         isValidParameters( parameters ),
         parameters->input->ok( parameters->input ) );

  Integer result = 0;
  const Integer dataSize =
    3 * calipso->maximumPoints +
    ( calipso->variables - 3 ) * calipso->maximumPoints * calipso->levels;

  /*
   * First reallocate calipso->data to be large enough to hold
   * all data for the largest profile, data[ variables ][ points ][ levels ].
   * Note that the first 3 variables (Profile_Time, Longitude, Latitude) are
   * only surface variables (i.e., have only 1 level).
   */

  FREE( calipso->data );
  calipso->data = NEW_ZERO( Real, dataSize );

  if ( calipso->data ) {
    Stream* output = newFileStream( "-stdout", "wb" );

    if ( output ) {
      int variable = 0;

      output->writeString( output, 
                           "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)"
                           "\tELEVATION(m)\t%s(%s)\t",
                           calipso->variable[ 0 ], calipso->units[ 0 ] );

      output->writeString( output, "%s(%s)",
                           calipso->variable[ 4 ], calipso->units[ 4 ] );

      for ( variable = 5; variable < calipso->variables; ++variable ) {
        output->writeString( output, "\t%s(%s)",
                            calipso->variable[ variable ],
                            calipso->units[ variable ] );
      }

      output->writeString( output, "\n" );

      if ( output->ok( output ) ) {
        result = writeASCIIData( calipso, parameters->input, output );
      }

      FREE_OBJECT( output );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCIIData - Write ASCII-format data lines.
INPUTS:  const CALIPSO* calipso  Structure to write.
         Stream* input           Stream to read from.
         Stream* output          Stream to write to.
RETURNS: Integer 1 if successful, else 0 and failureMesage is called.
NOTES:   Requires calipso->data[ variables ][ maximumPoints ][ levels ].
******************************************************************************/

static Integer writeASCIIData( CALIPSO* calipso, Stream* input,
                               Stream* output ) {

  PRE05( isValidCALIPSO( calipso ), input, input->isReadable( input ),
         output, output->isWritable( output ) );

  Integer result = 0;
  const Integer variables = calipso->variables;
  const char* const dataFormat =
    "%s\t%28.6"REAL_F_FORMAT"\t%28.6"REAL_F_FORMAT"\t%28.6"REAL_F_FORMAT
    "\t%28.6"REAL_F_FORMAT"\t%28.6"REAL_F_FORMAT;
  const Integer profiles = calipso->profiles;
  Integer profile = 0;

  /* Write data points for each timestep, profile, point, level, variable: */

  do {
    const Integer profile2 = profile + profile;
    const Integer groundPoints = calipso->dimensions[ profile2 + POINT ];
    const Integer levels       = calipso->dimensions[ profile2 + LEVEL ];
    const Integer profilePoints = groundPoints * levels;
    const Integer profileSize =
      3 * groundPoints + ( variables - 3 ) * profilePoints;
    const Integer profileTimestamp = calipso->timestamps[ profile ];
    UTCTimestamp profileUTCTimestamp = "";
    toUTCTimestamp( profileTimestamp, profileUTCTimestamp );

    CHECK( IN_RANGE( profileSize, 0,
                     calipso->variables * calipso->maximumPoints *
                     calipso->levels ) );

    /*
     * The first 3 variables (Profile_Time, Longitude, Latitude) are only
     * ground points (levels = 1) whereas the rest of the variables
     * (Elevation, Total_Attenuated_Backscatter_532) have levels >= 1.
     */

    input->read64BitReals( input, calipso->data, profileSize );

    if ( ! input->ok( input ) ) {
      profile = profiles;
    } else {
      const Real* const values = calipso->data;
      Integer value = 0;

      do {
        const Integer profileTimeOffset = value / levels;
        const Integer longitudeOffset   = profileTimeOffset + groundPoints;
        const Integer latitudeOffset    = longitudeOffset   + groundPoints;
        const Integer elevationOffset   = 3 * groundPoints  + value;
        const Integer dataOffset        = elevationOffset   + profilePoints;
        const Real profileTime          = values[ profileTimeOffset ];
        const Real longitude            = values[ longitudeOffset ];
        const Real latitude             = values[ latitudeOffset ];
        const Real elevation            = values[ elevationOffset ];
        const Real datum = MAX( values[ dataOffset ], -9999.0 );
        Integer v = 0;

        CHECK( dataOffset < 3 * groundPoints + 2 * groundPoints * levels );
        /* late. */

        output->writeString( output, dataFormat,
                             profileUTCTimestamp,
                             longitude,
                             latitude,
                             elevation,
                             profileTime,
                             datum );


        for ( v = 5; AND2( v < calipso->variables, output->ok(output)); ++v ) {
          const Integer offset = dataOffset + profilePoints * ( v - 4 );
          const Real datum2 = MAX( values[ offset ], -9999.0 );
          output->writeString( output, "\t%28.6"REAL_F_FORMAT, datum2 );
        }

        if ( output->ok( output ) ) {
          output->writeString( output, "\n" );          
        }

        if ( ! output->ok( output ) ) {
          profile = profiles;
          value = profilePoints;
        }

        ++value;
      } while ( value < profilePoints );
    }

    ++profile;
  } while ( profile < profiles );

  result = AND2( input->ok( input ), output->ok( output ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDS - Read input and write it in COARDS format to output.
INPUTS:  CALIPSO* calipso              Object to translate.
         const Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  Parameters updated by translation.
******************************************************************************/

static Integer writeCOARDS( CALIPSO* calipso, const Parameters* parameters ) {

  PRE04( calipso, isValidCALIPSO( calipso ),
         isValidParameters( parameters ),
         parameters->input->ok( parameters->input ) );

  Integer result = 0;
  const Integer surfaceVariables = 5;
  /* time, longitude, latitude, yyyyddd, hhmmss. */
  const Integer levelVariables = calipso->variables - 3;
  const Integer fileSizeEstimate =
    surfaceVariables * calipso->points * BYTES_PER_NETCDF_FLOAT +
    levelVariables * calipso->points * calipso->levels * BYTES_PER_NETCDF_FLOAT
    + 10000; /* header/extra. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {

    if ( writeCOARDSHeader( file, calipso ) ) {
      result = writeCOARDSData( parameters->input, file, calipso );
    }

    nc_close( file );
    file = -1;
  }

  POST03( IS_BOOL( result ),
          isValidCALIPSO( calipso ), isValidParameters( parameters ) );

  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSHeader - Write calipso header info to outputFileName.
INPUTS:  Integer file NetCDF file to write to.
         const CALIPSO* calipso Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSHeader( Integer file, const CALIPSO* calipso ) {

  PRE02( file != -1, isValidCALIPSO( calipso ) );

  Integer result = 0;
  const char* const names[ DIMENSIONS ] = { "points", "levels" };
  Integer dimensionIds[ DIMENSIONS ] = { -1, -1 };
  Integer dimensions[   DIMENSIONS ] = { 0, 0 };
  dimensions[ POINT ] = calipso->points;
  dimensions[ LEVEL ] = calipso->levels;

  if ( createDimensions( file, DIMENSIONS, names, dimensions, dimensionIds)) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createLongitudeAndLatitude( file, 1, dimensionIds ) ) {
        char timeUnits[ 80 ] = "";
        strcpy( timeUnits, calipso->units[ 0 ] );
        underscoreToSpace( timeUnits );

        if ( createVariable( file, calipso->variable[ 0 ], timeUnits,
                             NC_DOUBLE, 0, 1, dimensionIds ) != -1 ) {
          const Integer variables = calipso->variables;
          Integer index = 3;

          do {
            const char* const units =
              strcmp( calipso->units[ index ], "-"   ) == 0 ? "none" :
              strcmp( calipso->units[ index ], "deg" ) == 0 ? "degrees" :
              calipso->units[ index ];
              const char* const name =
                strcmp( calipso->variable[ index ], "Elevation" ) ?
                calipso->variable[ index ] : "elevation";
            const Integer ok =
              createVariable( file, name, units, NC_FLOAT,
                              strcmp( name, "elevation" ) != 0, 2,
                              dimensionIds ) != -1;

             if ( ! ok ) {
               index = variables;
             }

             ++index;
          } while ( index < variables );

          if ( index == variables ) {

            if ( writeExtraAttributes( file, (const Real (*)[2]) calipso->domain,
                                       dimensionIds[ POINT ] ) ) {
              UTCTimestamp timestamp = "";
              const char* const history =
                "http://eosweb.larc.nasa.gov/HORDERBIN/HTML_Start.cgi/"
                ",CALIPSOSubset,XDRConvert";
              toUTCTimestamp( calipso->timestamps[ 0 ], timestamp );
              result = writeStandardContents( file, history, timestamp,
                                              dimensionIds[ POINT ],
                                              calipso->points, 0 );
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSData - Write calipso data info to outputFileName.
INPUTS:  Stream* input    Stream to read data arrays from.
         Integer file     NetCDF file to write to.
         const CALIPSO* calipso Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSData( Stream* input, Integer file,
                                 const CALIPSO* calipso ) {

  PRE05( input, input->ok( input ), input->isReadable( input ),
         file != -1, isValidCALIPSO( calipso ) );

  Integer result = 0;
  const Integer profiles  = calipso->profiles;
  const Integer levels    = calipso->levels;
  const Integer variables = calipso->variables;
  Integer profile = 0;
  Integer variable = 0;
  Integer offset = 0;

  DEBUG( fprintf( stderr, "profiles: %lld levels: %lld variables: %lld\n",
                  profiles, levels, variables ); )

  do {
    const Integer points = calipso->dimensions[ profile * 2 + POINT ];
    variable = 0;

    DEBUG( fprintf( stderr, "points: %lld\n", points ); )

    do {
      const char* const variableName =
        variable == 0 ? calipso->variable[0] :
        variable == 1 ? "longitude" :
        variable == 2 ? "latitude" :
        variable == 3 ? "elevation" :
        calipso->variable[ variable ];
      const Integer count = variable < 3 ? points : points * levels;
      const Integer dimension1 = points;
      const Integer dimension2 = variable < 3 ? 1 : levels;

      DEBUG( fprintf( stderr, "%s count: %lld offset: %lld dims: %lld %lld\n",
                      variableName, count, offset, dimension1, dimension2 ); )

      input->read64BitReals( input, calipso->data, count );

      if ( ! AND2( input->ok( input ),
                   writeSomeData( file, variableName, offset,
                                  dimension1, dimension2, 1, 1,
                                  calipso->data ) ) ) {
        profile = profiles;
        variable = variables;
      }

      ++variable;
    } while ( variable < variables );

    offset += points;
    ++profile;
  } while ( profile < profiles );

  result = AND3( profile == profiles,
                 variable == variables,
                 writeTimeData( file, calipso->profiles, 2, 0,
                                calipso->timestamps, calipso->dimensions,
                                calipso->data ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  CALIPSO* calipso              Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedXDR( CALIPSO* calipso,
                                  const Parameters* parameters) {

  PRE02( isValidCALIPSO( calipso ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = calipso->timesteps;
    const Integer points    = calipso->totalRegriddedPoints;
    const Integer layers    = calipso->gridLayers;
    const Integer variableIndex = calipso->variables >= 5 ? 4 : 0;
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
    Name variable = "";
    aggregateName( calipso->variable[ variableIndex ], hoursPerTimestep,
                   variable );

    output->writeString( output,
                         "REGRIDDED-CALIPSO 2.0\n"
                         "%s,XDRConvert\n"
                         "%s\n"
                         "# timesteps layers\n"
                         "%"INTEGER_FORMAT" %"INTEGER_FORMAT"\n"
                         "# Variable name:\n%s\n"
                         "# Variable units:\n%s\n",
                         calipso->note,
                         calipso->timestamp,
                         timesteps, layers,
                         variable,
                         calipso->units[ variableIndex ] );

    writeProjectionAndGrid( parameters->grid, output );

    output->writeString( output,
          "# MSB 64-bit integers points[timesteps] and\n"
          "# IEEE-754 64-bit reals longitudes[timesteps][points] and\n"
          "# IEEE-754 64-bit reals latitudes[timesteps][points] and\n"
          "# IEEE-754 64-bit reals elevations[timesteps][points][layers] and\n"
          "# MSB 64-bit integers columns[timesteps][points] and\n"
          "# MSB 64-bit integers rows[timesteps][points] and\n"
          "# MSB 64-bit integers layers[timesteps][points][layers] and\n"
          "# IEEE-754 64-bit reals data[timesteps][points][layers]:\n" );

    if ( output->ok( output ) ) {
      output->write64BitIntegers( output, calipso->outputPoints, timesteps );

      if ( output->ok( output ) ) {
        output->write64BitReals( output, calipso->gridLongitudes, points );

        if ( output->ok( output ) ) {
          output->write64BitReals( output, calipso->gridLatitudes, points );

          if ( output->ok( output ) ) {
            output->write64BitReals( output, calipso->gridElevations,
                                     points * layers );

            if ( output->ok( output ) ) {
              output->write64BitIntegers( output, calipso->columns, points);

              if ( output->ok( output ) ) {
                output->write64BitIntegers( output, calipso->rows, points );

                if ( output->ok( output ) ) {
                  output->write64BitIntegers( output, calipso->layers,
                                              points * layers );

                  if ( output->ok( output ) ) {
                    output->write64BitReals( output, calipso->gridData,
                                             points * layers );

                    result = output->ok( output );
                  }
                }
              }
            }
          }
        }
      }
    }

    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedASCII - Write regridded ASCII-format data.
INPUTS:  CALIPSO* calipso              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedASCII( CALIPSO* calipso,
                                    const Parameters* parameters ) {

  PRE03( isValidCALIPSO( calipso ), IN3( calipso->variables, 1, 5 ),
         isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const char* const headerStart =
      "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tELEVATION(m)"
      "\tCOLUMN(-)\tROW(-)\tLAYER(-)";
    const char* const headerFormat = "\t%s(%s)\n";
    const char* const dataFormat =
      "%s\t%10.4lf\t%10.4lf\t%10.4lf"
      "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT
      "\t%10.4lf\n";

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {
      const Integer variableIndex = calipso->variables >= 5 ? 4 : 0;
      const Integer hoursPerTimestep =
        parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
      Name variable = "";
      aggregateName( calipso->variable[ variableIndex ], hoursPerTimestep,
                     variable );

      output->writeString( output, headerFormat,
                           variable, calipso->units[ variableIndex ] );

      if ( output->ok( output ) ) {
        const Integer layers    = calipso->gridLayers;
        const Integer timesteps = calipso->timesteps;
        const Real* longitudes  = calipso->gridLongitudes;
        const Real* latitudes   = calipso->gridLatitudes;
        const Integer* columns  = calipso->columns;
        const Integer* rows     = calipso->rows;
        const Real* data        = calipso->gridData;
        const Real* elevations  = calipso->gridElevations;
        Integer timestep = 0;
        Integer yyyydddhh00 =
          ( fromUTCTimestamp( calipso->timestamp ) / 100 ) * 100;
        UTCTimestamp timestamp;

        /* Write data rows: */

        do {
          const Integer points = calipso->outputPoints[ timestep ];
          Integer point = 0;
          toUTCTimestamp( yyyydddhh00, timestamp );

          for ( point = 0; point < points; ++point ) {
            const Real longitude = *longitudes++;
            const Real latitude  = *latitudes++;
            const Integer column = *columns++;
            const Integer row    = *rows++;
            Integer layer = 0;

            do {
              const Real value0    = *data++;
              const Real value     = MAX( value0, -9999.0 );
              const Real elevation = *elevations++;

              output->writeString( output, dataFormat,
                                   timestamp, longitude, latitude, elevation,
                                   column, row, layer + 1, value );

              if ( ! output->ok( output ) ) {
                layer = layers;
                point = points;
                timestep = timesteps;
              }

              ++layer;
            } while ( layer < layers );
          }

          yyyydddhh00 = offsetTimestamp( yyyydddhh00, hoursPerTimestep );
          ++timestep;
        } while ( timestep < timesteps );
      }
    }

    result = output->ok( output );
    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDS - Write regridded COARDS-format data.
INPUTS:  CALIPSO* calipso              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDS( CALIPSO* calipso,
                                     const Parameters* parameters ) {

  PRE02( isValidCALIPSO( calipso ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer points          = calipso->totalRegriddedPoints;
  const Integer levels          = calipso->gridLayers;
  const Integer headerBytes     = 10000;
  const Integer bytesPerDatum   = 4;
  const Integer groundVariables = 5; /* lon[ points ], lat, col, row, time. */
  const Integer pointsTimesLevelVariables = 2; /* data, elevations. */
  const Integer groundSize      = points * groundVariables;
  const Integer pointsTimesLevelSize =
    points * levels * pointsTimesLevelVariables;
  const Integer fileSizeEstimate =
    ( groundSize + pointsTimesLevelSize ) * bytesPerDatum + headerBytes;

  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedCOARDSHeader( file, hoursPerTimestep, calipso ) ) {
      result = writeRegriddedCOARDSData( file, calipso, parameters );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDSHeader - Write header to file.
INPUTS:  Integer file            NetCDF file to write to.
         Integer hoursPerTimestep  E.g., 1, 24, 744, etc.
         const CALIPSO* calipso  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const CALIPSO* calipso ) {

  PRE03( file != -1, hoursPerTimestep > 0, isValidCALIPSO( calipso ) );

  Integer result = 0;
  const Integer variableIndex = calipso->variables >= 5 ? 4 : 0;
  const char* const dimensionNames[ 2 ] = { "points", "levels" };
  Integer dimensionIds[ 2 ] = { -1, -1 };
  Integer dimensions[ 2 ] = { 0, 0 };
  dimensions[ 0 ] = calipso->totalRegriddedPoints;
  dimensions[ 1 ] = calipso->gridLayers;

  DEBUG( fprintf( stderr, "dims: %lld %lld\n", dimensions[0], dimensions[1]);)

  if ( createDimensions( file, 2, dimensionNames, dimensions, dimensionIds)) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "column", "-",
                           NC_INT, 0, 1, &dimensionIds[ 0 ] ) != -1 ) {

        if ( createVariable( file, "row", "-",
                             NC_INT, 0, 1, &dimensionIds[ 0 ] ) != -1 ) {

          if ( createLongitudeAndLatitude( file, 1, &dimensionIds[ 0 ] ) ) {

            if ( createVariable( file, "elevation", "m",
                                 NC_FLOAT, 0, 2, dimensionIds ) != -1 ) {
              Name variable = "";
              aggregateName( calipso->variable[ variableIndex ],
                             hoursPerTimestep, variable );

              if ( createVariable( file,
                                   variable,
                                   calipso->units[ variableIndex ],
                                   NC_FLOAT, 1, 2, dimensionIds ) != -1 ) {

                Line history = "";
                appendToLine( history, calipso->note );
                appendToLine( history, ",XDRConvert" );

                result = writeStandardContents( file, history,
                                                calipso->timestamp,
                                                dimensionIds[ 0 ], 0, 0 );
              }
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDSData - Write COARDS-format data to file.
INPUTS:  Integer file            NetCDF file to write to.
         CALIPSO* calipso  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSData( Integer file,
                                         CALIPSO* calipso,
                                         const Parameters* parameters ) {

  PRE03( file != -1, isValidCALIPSO( calipso ),
         isValidParameters( parameters ) );

  Integer result = 0;
  const Integer points = calipso->totalRegriddedPoints;
  const Integer levels = calipso->gridLayers;

  if ( writeAllIntData( file, "column", points, 1, 1, 1,
                        calipso->columns ) ) {

    if ( writeAllIntData( file, "row", points, 1, 1, 1,
                          calipso->rows ) ) {

      if ( writeAllData( file, "longitude", points, 1, 1, 1,
                         calipso->gridLongitudes ) ) {

        if ( writeAllData( file, "latitude", points, 1, 1, 1,
                           calipso->gridLatitudes ) ) {

          if ( writeAllData( file, "elevation",
                             points, levels, 1, 1, calipso->gridElevations )) {
            const Integer variableIndex = calipso->variables >= 5 ? 4 : 0;
            const Integer hoursPerTimestep =
              parameters->aggregationTimesteps ? parameters->aggregationTimesteps
              : 1;
            Name variable = "";
            aggregateName( calipso->variable[ variableIndex ],
                           hoursPerTimestep, variable );

            replaceMissingValues( points * levels, calipso->gridData );

            if ( writeAllData( file, variable, points, levels, 1, 1,
                               calipso->gridData ) ) {
              timeData( calipso->timesteps, hoursPerTimestep, points,
                        calipso->outputPoints, calipso->gridData );

              result = writeAllData( file, "time", points, 1, 1, 1,
                                     calipso->gridData );
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPI - Write regridded IOAPI-format data.
INPUTS:  CALIPSO* calipso              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPI( CALIPSO* calipso,
                                    const Parameters* parameters ) {

  PRE02( isValidCALIPSO( calipso ), isValidParameters( parameters ) );

  Integer result = 0;
  const Grid* const grid = parameters->grid;
  const Integer timesteps = calipso->timesteps;
  const Integer layers   = grid->layers( grid );
  const Integer rows     = grid->rows( grid );
  const Integer columns  = grid->columns( grid );
  const Integer headerBytes   = 10000;
  const Integer bytesPerDatum = 4;
  const Integer variables     = 4; /* lon, lat, elv, data. */
  const Integer dataSize =
    variables * timesteps * layers * rows * columns * bytesPerDatum;
  const Integer timeSize = timesteps * variables * 2 * bytesPerDatum;
  const Integer fileSizeEstimate =
    dataSize + timeSize + headerBytes;
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedIOAPIHeader( file, hoursPerTimestep,
                                    calipso, parameters->grid ) ) {
      result =
        writeRegriddedIOAPIData( file, hoursPerTimestep, calipso,
                                 parameters->grid );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIHeader - Write header to file.
INPUTS:  Integer file            NetCDF file to write to.
         Integer hoursPerTimestep  Hours per timestep: 1, 24, etc.
         const CALIPSO* calipso  Structure to write.
         const Grid* grid        Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const CALIPSO* calipso,
                                          const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0,
         isValidCALIPSO( calipso ), grid, grid->invariant( grid ));

  Integer result = 0;
  const Integer layers = MIN( calipso->gridLayers, grid->layers( grid ) );
  enum { VARIABLES = 4 }; /* LONGITUDE, LATITUDE, ELEVATION, calipso. */
  Name variableNames[ VARIABLES ] = {
    "LONGITUDE", "LATITUDE", "ELEVATION", "calipso"
  };
  Name variableUnits[ VARIABLES ] = { "deg", "deg", "m", "-" };
  const Integer timestamp = fromUTCTimestamp( calipso->timestamp );
  const Integer variableIndex = calipso->variables >= 5 ? 4 : 0;
  Line history = "";
  appendToLine( history, calipso->note );
  appendToLine( history, ",XDRConvert" );

  aggregateName( calipso->variable[ variableIndex ], hoursPerTimestep,
                 variableNames[ VARIABLES - 1 ] );
  variableNames[ VARIABLES - 1 ][ 15 ] = '\0';
  strncpy( variableUnits[ VARIABLES - 1 ],
           calipso->units[ variableIndex ], 16 );
  variableUnits[ VARIABLES - 1 ][ 16 ] = '\0';
  uppercase( variableNames[ VARIABLES - 1 ] );
  lowercase( variableUnits[ VARIABLES - 1 ] );

  result = writeM3IOHeader( file, calipso->timesteps, hoursPerTimestep,
                            timestamp,
                            VARIABLES, layers,
                            (const Name*) variableNames,
                            (const Name*) variableUnits,
                            history, grid );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIData - Write IOAPI-format data to file.
INPUTS:  Integer file              NetCDF file to write to.
         Integer hoursPerTimestep  E.g., 1, 24, 744, etc.
         const CALIPSO* calipso    Structure to write.
         const Grid* grid          Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const CALIPSO* calipso,
                                        const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0, isValidCALIPSO( calipso ),
         grid, grid->invariant( grid ) );

  Integer result = 0;
  const Integer layers    = MIN( calipso->gridLayers, grid->layers( grid ) );
  const Integer rows      = grid->rows( grid );
  const Integer columns   = grid->columns( grid );
  const Integer cells     = layers * rows * columns;
  const Integer variables = 2; /* elv, data. */
  Real* gridElevations    = NEW_ZERO( Real, cells * variables );
  Real* expandedGridData  = gridElevations ? gridElevations + cells : 0;

  if ( gridElevations ) {
    const Integer variableIndex = calipso->variables >= 5 ? 4 : 0;
    const Integer timesteps = calipso->timesteps;
    Integer timestep = 0;
    Integer offset   = 0;
    Name variable    = "";
    aggregateName( calipso->variable[ variableIndex ], hoursPerTimestep,
                   variable );
    variable[ 15 ] = '\0';
    uppercase( variable );

    if ( writeM3IOGrid( grid, timesteps, layers, file ) ) {

      do {
        const Integer points = calipso->outputPoints[ timestep ];
        const Integer* const rowData    = calipso->rows    + offset;
        const Integer* const columnData = calipso->columns + offset;
        const Integer offset2 = offset * layers;
        const Real* const elevationData = calipso->gridElevations + offset2;
        const Real* const regriddedData = calipso->gridData       + offset2;
        const Real scale = 1.0;
        Integer ok = 0;

        DEBUG(fprintf(stderr,"writing timestep %"INTEGER_FORMAT"\n",timestep);)

        copyDataToGrid( points, rowData, columnData, elevationData,
                        scale, layers, rows, columns, gridElevations );

        if ( writeM3IOData( file, "ELEVATION",
                            timestep, layers, rows, columns,
                            gridElevations ) ) {

          copyDataToGrid( points, rowData, columnData, regriddedData,
                          scale, layers, rows, columns, expandedGridData );

          ok = writeM3IOData( file, variable,
                              timestep, layers, rows, columns,
                              expandedGridData );
        }

        if ( ! ok ) {
          timestep = timesteps;
        }

        offset += points;
        ++timestep;
      } while ( timestep < timesteps );

      result = timestep == timesteps;
    }

    FREE( gridElevations );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: regridCALIPSO - Regrid data.
INPUTS:  Stream* input     Stream to read CALIPSO data from.
         Integer method    E.g., AGGREGATE_MEAN.
         Grid* grid        Grid to project and aggregate points into.
         CALIPSO* calipso  Data to regrid.
OUTPUTS: CALIPSO* calipso  Regridded data.
******************************************************************************/

static void regridCALIPSO( Stream* input, Integer method, Grid* grid,
                           CALIPSO* calipso ) {

  PRE08( input, input->isReadable( input ),
         IS_VALID_AGGREGATE_METHOD( method ),
         grid, grid->invariant( grid ),
         isValidCALIPSO( calipso ),
         calipso->totalRegriddedPoints == 0,
         calipso->longitudes == 0 );

  const Integer variables = calipso->variables;

  if ( variables >= 5 ) {
    const Integer gridColumns   = grid->columns( grid );
    const Integer gridRows      = grid->rows( grid );
    const Integer gridLayers    = calipso->levels > 1 ? grid->layers(grid) : 1;
    const Integer timesteps     = calipso->timesteps;
    const Integer profiles      = calipso->profiles;
    const Integer maximumPoints = calipso->maximumPoints;/*Longest profile*/
    const Integer levels        = calipso->levels;
    /* profile: elevation, data, optional thickness. */
    const Integer levelVariables        = variables - 3;
    const Integer inputGroundVariables  = 2; /* profile: lon, lat. */
    const Integer outputGroundVariables = 4; /* grid: lon, lat, col, row. */
    const Integer layerVariables        = 3; /* grid: layer, data, elevation.*/
    const Integer maximumProfilesPerTimestep  = 2; /* 45 minute half orbit. */
    const Integer inputGroundSize = maximumProfilesPerTimestep * maximumPoints;
    const Integer maximumRegriddedPoints =
      maximumProfilesPerTimestep * ( gridColumns + gridRows ); /*Worst-case?*/
    const Integer outputGroundSize =
      MIN( profiles, timesteps ) * maximumRegriddedPoints;
    const Integer levelSize         = inputGroundSize * levels;
    const Integer layerSize         = outputGroundSize * gridLayers;
    const Integer processGroundSize = inputGroundSize + outputGroundSize;
    const Integer processLayerSize  = levelSize + layerSize;
    const Integer inputDataSize =
      levelVariables * levelSize + inputGroundVariables * inputGroundSize;
    const Integer outputDataSize =
      layerVariables * processLayerSize +
      outputGroundVariables * processGroundSize +
      timesteps;
    const Integer dataSize = inputDataSize + outputDataSize;
    FREE( calipso->data );
    calipso->data = NEW_ZERO( Real, dataSize );

    DEBUG( fprintf( stderr, "levelSize = %lld, inputGroundSize = %lld\n",
                    levelSize, inputGroundSize ); )

    DEBUG( fprintf( stderr, "layerSize = %lld, outputGroundSize = %lld\n",
                    layerSize, outputGroundSize ); )

    DEBUG( fprintf( stderr, "inputDataSize = %lld, outputDataSize = %lld\n",
                    inputDataSize, outputDataSize ); )

    DEBUG( fprintf( stderr, "processGroundSize = %lld\n",
                    processGroundSize ); )

    if ( calipso->data ) {
      Integer totalRegriddedPoints = 0;
      Integer timestep = 0;
      Integer yyyydddhh00 =
        ( fromUTCTimestamp( calipso->timestamp ) / 100 ) * 100;

      /* Input data: calipso->data,elevations,longitudes,latitudes: */

      calipso->elevations = calipso->data + levelSize;
      calipso->thickness = variables <= 5 ? 0 : calipso->elevations +levelSize;
      calipso->longitudes =
        calipso->thickness ? calipso->thickness + levelSize :
        calipso->elevations + levelSize;
      calipso->latitudes = calipso->longitudes + inputGroundSize;

      /*
       * Output data: calipso->gridLongitudes,GridLatitudes,columns,rows,
       * gridData,gridElevations,outputPoints:
       */

      calipso->gridLongitudes = calipso->latitudes      + inputGroundSize;
      calipso->gridLatitudes  = calipso->gridLongitudes + processGroundSize;
      calipso->columns        = (Integer*)
                                calipso->gridLatitudes  + processGroundSize;
      calipso->rows           = calipso->columns        + processGroundSize;
      calipso->layers         = calipso->rows           + processGroundSize;
      calipso->gridData       = (Real*) calipso->layers + processLayerSize;
      calipso->gridElevations = calipso->gridData       + processLayerSize;
      calipso->outputPoints   = (Integer*)
                                calipso->gridElevations + processLayerSize;

      /* For these discrete variables use AGGREGATE_NEAREST method: */

      if ( strstr( "Atmospheric_Volume_Description "
                   "Feature_Classification_Flags "
                   "Profile_ID "
                   "Land_Water_Mask "
                   "Cirrus_Shape_Parameter "
                   "Horizontal_Averaging "
                   "Number_Layers_Found "
                   "Layer_Base_Extended "
                   "IGBP_Surface_Type "
                   "NSIDC_Surface_Type "
                   "Frame_Number "
                   "Lidar_Mode "
                   "Lidar_Submode "
                   "Surface_Elevation_Detection_Frequency "
                   "FeatureFinderQC "
                   "QC_FLag "
                   "QC_Flag_2 "
                   "Day_Night_Flag "
                   "Frame_Number "
                   "Lidar_Mode "
                   "Lidar_Submode "
                   "Opacity_Flag "
                   "Cirrus_Shape_Parameter "
                   "Cirrus_Shape_Parameter_Uncertainty "
                   "Cirrus_Shape_Parameter_Invalid_Points "
                   "Number_Layers_Found "
                   "Samples_Averaged "
                   "FeatureFinderQC "
                   "CAD_Score "
                   "ExtinctionQC_532 "
                   "ExtinctionQC_1064 "
                   "Lidar_Ratio_532_Selection_Method "
                   "Lidar_Ratio_1064_Selection_Method ",
                   calipso->variable[ 4 ] ) ) {
        method = AGGREGATE_NEAREST;
      }

      DEBUG( fprintf( stderr,
                      "method: %lld, variables = %lld, levelVariables = %lld, "
                      "calipso->thickness = %p, "
                      "inputDataSize = %lld, outputDataSize = %lld\n",
                      method, variables, levelVariables, calipso->thickness,
                      inputDataSize, outputDataSize ); )

      if ( skipProfileDataBeforeTimestamp( yyyydddhh00, calipso, input ) ) {

        do {
          Integer inputPoints = 0;

          if ( ! readProfileDataForTimestamp( yyyydddhh00, input, calipso,
                                              &inputPoints ) ) {
            timestep = timesteps;
          } else if ( inputPoints > 0 ) {
            const Integer layerOffset = totalRegriddedPoints * gridLayers;
            Integer outputPoints = 0;
            const Real minimumValidValue =
              strstr( calipso->variable[ 4 ], "emperature" ) ? -120.0 : 0.0;

            DEBUG( fprintf( stderr,
                            "inputPoints = %lld, minimumValidValue = %lf\n",
                             inputPoints, minimumValidValue ); )
            CHECK2( inputPoints <= inputGroundSize,
                    totalRegriddedPoints < outputGroundSize );

            grid->regrid( grid, method, minimumValidValue, inputPoints,
                          calipso->levels,
                          calipso->longitudes, calipso->latitudes,
                          calipso->levels > 1 ? calipso->elevations : 0,
#if 0
/* UNIMPLEMENTED: thickness-based regriding in Utilities/Grid.c,regrid() */
                          calipso->thickness,
#endif
                          calipso->data,
                          0, /* No input vector data. */
                          0, /* No notes. */
                          &outputPoints,
                          calipso->columns        + totalRegriddedPoints,
                          calipso->rows           + totalRegriddedPoints,
                          0,
                          calipso->gridLongitudes + totalRegriddedPoints,
                          calipso->gridLatitudes  + totalRegriddedPoints,
                          calipso->levels > 1 ?
                            calipso->gridElevations + layerOffset : 0,
                          calipso->gridData       + layerOffset,
                          0, /* No output vector data. */
                          0 /* no regriddedNotes */ );

            calipso->outputPoints[ timestep ] = outputPoints;
            totalRegriddedPoints += outputPoints;
          }

          incrementTimestamp( &yyyydddhh00 );
          ++timestep;
        } while ( timestep < timesteps );

        if ( totalRegriddedPoints ) {
          calipso->totalRegriddedPoints = totalRegriddedPoints;
          calipso->gridLayers = gridLayers;
          computeLayers( calipso );
        }
      }
    }
  }

  POST02( calipso->totalRegriddedPoints >= 0,
          IMPLIES( calipso->totalRegriddedPoints > 0,
                   AND10( IN_RANGE( minimumItemI( calipso->outputPoints,
                                                  calipso->timesteps ),
                                    0, calipso->totalRegriddedPoints ),
                          IN_RANGE( maximumItemI( calipso->outputPoints,
                                                  calipso->timesteps ),
                                    1, calipso->totalRegriddedPoints ),
                          IN_RANGE( minimumItemI(calipso->columns,
                                                calipso->totalRegriddedPoints),
                                                 1, grid->columns( grid ) ),
                          IN_RANGE( maximumItemI(calipso->columns,
                                                calipso->totalRegriddedPoints),
                                                 1, grid->columns( grid ) ),
                          IN_RANGE( minimumItemI(calipso->rows,
                                                calipso->totalRegriddedPoints),
                                                 1, grid->rows( grid ) ),
                          IN_RANGE( maximumItemI(calipso->rows,
                                                calipso->totalRegriddedPoints),
                                                 1, grid->rows( grid ) ),
                          validLongitudesAndLatitudes(
                                                calipso->totalRegriddedPoints,
                                                       calipso->gridLongitudes,
                                                       calipso->gridLatitudes),
                          isNanFree( calipso->elevations,
                                     calipso->totalRegriddedPoints ),
                          isNanFree( calipso->gridData,
                                     calipso->totalRegriddedPoints *
                                     grid->layers( grid ) ),
                          isNanFree( calipso->gridElevations,
                                     calipso->totalRegriddedPoints *
                                     grid->layers( grid ) ) ) ) );
}



/******************************************************************************
PURPOSE: skipProfileDataBeforeTimestamp - Skip data before a given timestamp.
INPUTS:  Integer  yyyydddhh00    Timestamp to compare.
         const CALIPSO* calipso  CALIPSO structure.
         Stream*  stream         Stream to read/skip CALIPSO data from.
OUTPUTS: Stream*  stream         Stream to read/skip CALIPSO data from.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer skipProfileDataBeforeTimestamp( Integer yyyydddhh00,
                                               const CALIPSO* calipso,
                                               Stream* stream ) {

  PRE05( isValidTimestamp( yyyydddhh00 ), calipso, isValidCALIPSO( calipso ),
          stream, stream->isReadable( stream ) );

  Integer result = 1;
  const Integer profiles = calipso->profiles;
  Integer profile = 0;
  Integer done = 0;
  enum { BUFFER_SIZE = 1024 * 1024 };
  char buffer[ BUFFER_SIZE ] = "";

  do {
    const Integer yyyydddhhmm = calipso->timestamps[ profile ];
    const Integer yyyydddhh002 = ( yyyydddhhmm / 100 ) * 100;

    DEBUG( fprintf( stderr, "%lld <? %lld\n", yyyydddhh002, yyyydddhh00 ); )

    if ( yyyydddhh002 < yyyydddhh00 ) {
      const Integer points = calipso->dimensions[ profile * 2 + POINT ];
      const Integer levels = calipso->dimensions[ profile * 2 + LEVEL ];
      const Integer hasThickness = calipso->thickness != 0;
      const Integer dataBytes =
        ( points * 3 + /* Profile_UTC_Time, Longitude, Latitude */
          points * levels * ( 2 + hasThickness ) ) * sizeof (Real);
      Integer bytesRemaining = dataBytes;
      Integer done = 0;
      DEBUG( fprintf( stderr, "skiping dataBytes = %lld\n", dataBytes ); )

      do {
        const Integer bytesToReadNow = MIN( bytesRemaining, BUFFER_SIZE );
        stream->readBytes( stream, buffer, bytesToReadNow );

        if ( stream->ok( stream ) ) {
          bytesRemaining -= bytesToReadNow;
          done = bytesRemaining == 0;
        } else {
          result = 0;
          done = 1;
          failureMessage( "Failed to read/skip %lld bytes of profile data.",
                          dataBytes );
        }
      } while ( ! done );

    } else {
      done = 1;
    }

    ++profile;
  } while ( AND2( profile < profiles, ! done ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readProfileDataForTimestamp - Read all data for a given timestamp
         for regridding.
INPUTS:  Integer  yyyydddhh00  Timestamp to compare.
         Stream*  input        Stream to read CALIPSO data from.
         CALIPSO* calipso      CALIPSO structure.
OUTPUTS: Integer* points       Number of points read for this timestamp.
         CALIPSO* calipso      data, longitudes, latitudes, elevations.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readProfileDataForTimestamp( Integer yyyydddhh00,
                                            Stream* input,
                                            CALIPSO* calipso,
                                            Integer* points ) {

  PRE010( isValidTimestamp( yyyydddhh00 ),
          input, input->isReadable( input ),
          calipso, isValidCALIPSO( calipso ), points,
          calipso->longitudes, calipso->latitudes, calipso->elevations,
          calipso->data );

  const Integer profiles = calipso->profiles;
  Integer result  = 0;
  Integer profile = 0;
  Integer groundOffset = 0;
  Integer levelOffset  = 0;
  *points = 0;

  do {
    const Integer yyyydddhhmm = calipso->timestamps[ profile ];
    const Integer yyyydddhh002 = ( yyyydddhhmm / 100 ) * 100;

    DEBUG( fprintf( stderr, "%lld =?= %lld\n", yyyydddhh002, yyyydddhh00 ); )

    if ( yyyydddhh002 == yyyydddhh00 ) {
      const Integer groundPoints = calipso->dimensions[ profile * 2 + POINT ];
      const Integer levels       = calipso->dimensions[ profile * 2 + LEVEL ];
      const Integer dataSize     = groundPoints * levels;
      Real* const thickness =
        calipso->thickness ? calipso->thickness + levelOffset : 0;

      DEBUG( fprintf( stderr, "groundPoints = %lld, levels = %lld\n",
                      groundPoints, levels ); )

      if ( ! readProfileData( input, calipso->variables, groundPoints, levels,
                              calipso->longitudes + groundOffset,
                              calipso->latitudes  + groundOffset,
                              calipso->elevations + levelOffset,
                              thickness,
                              calipso->data       + levelOffset ) ) {
        profile = profiles;
      } else {
        *points      += groundPoints;
        groundOffset += groundPoints;
        levelOffset  += dataSize;
      }
    }

    ++profile;
  } while ( profile < profiles );

  result = profile == profiles;

  POST03( IS_BOOL( result ),
          *points >= 0,
          IMPLIES( AND2( result, *points > 0 ),
                   AND3( validLongitudesAndLatitudes( *points,
                                                      calipso->longitudes,
                                                      calipso->latitudes ),
                         isNanFree( calipso->elevations,
                                    *points * calipso->levels ),
                         isNanFree( calipso->data,
                                    *points * calipso->levels ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readProfileData - Read all variables of a profile for regridding.
INPUTS:  Stream* input               Stream to read CALIPSO data from.
         Integer variables           Number of profile variables.
         Integer points              Number of profile points to read.
         Integer levels              Number of profile levels to read.
OUTPUTS: Real* longitudes[ points ]  Longitude data.
         Real* latitudes[ points ]   Latitude data.
         Real* elevations[ points * levels ]  Elevation data.
         Real* thickness[ points * levels ]   Optional thickness data.
         Real* data[ points * levels ]        5th data variable. Skip rest.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readProfileData( Stream* input,
                                Integer variables, Integer points,
                                Integer levels,
                                Real longitudes[], Real latitudes[],
                                Real elevations[], Real thickness[],
                                Real data[] ) {

  PRE09( input, input->isReadable( input ), variables >= 5, points > 0,
         levels > 0,
         longitudes,latitudes, elevations, data );

  Integer result = 0;
  const Integer dataSize = points * levels;
  Real* skipData =
    OR2( variables > 6, AND2( variables == 6, elevations == 0 ) ) ?
    NEW_ZERO( Real, dataSize ) : 0;

  DEBUG( fprintf( stderr, "variables = %lld, points = %lld, levels = %lld\n",
                  variables, points, levels ); )

  if ( IMPLIES( OR2( variables > 6, AND2( variables == 6, elevations == 0 ) ),
                skipData ) ) {
    Integer variable = 0;

    do {
      Real* const output =
        variable == 0 ? longitudes : /* Read/ignore profile_time. */
        variable == 1 ? longitudes :
        variable == 2 ? latitudes  :
        variable == 3 ? elevations :
        variable == 4 ? data :
        variable == 5 ? ( thickness != 0 ? thickness : skipData ) :
        skipData;
      const Integer count =
        IN3( output, longitudes, latitudes ) ? points : dataSize;
      CHECK( output );
      input->read64BitReals( input, output, count );
      DEBUG( fprintf( stderr,
                      "%lld: count = %lld, output = %p [%lf, ..., %lf]\n",
                      variable, count, output, output[ 0 ],
                      output[ count - 1 ] ); );
      ++variable;
    } while ( AND2( input->ok( input ), variable < variables ) );

    if ( AND2( input->ok( input ), variable == variables ) ) {

      if ( ! validLongitudesAndLatitudes( points, longitudes, latitudes ) ) {
        failureMessage( "Input data has invalid longitudes and/or latitudes.");
      } else if ( ! isNanFree( elevations, points * levels ) ) {
        failureMessage( "Input data has invalid elevations." );
      } else if ( ! isNanFree( data, points * levels ) ) {
        failureMessage( "Input data is invalid." );
      } else {
        result = 1;
      }
    }
  }

  FREE( skipData );

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND4( validLongitudesAndLatitudes( points,
                                                      longitudes, latitudes ),
                         isNanFree( elevations, points * levels ),
                         IMPLIES( thickness,
                                  isNanFree( thickness, points * levels ) ),
                         isNanFree( data, points * levels ) ) ) );
  return result;
}





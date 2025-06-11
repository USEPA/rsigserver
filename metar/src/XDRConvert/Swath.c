
/******************************************************************************
PURPOSE: Data.c - Define routines for processing satellite-measured swath data.

NOTES:

HISTORY: 2008-02 plessel.todd@epa.gov
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h>  /* For stderr, fprintf(). */
#include <time.h>  /* For time_t, time(), ctime(). */
#endif
#include <string.h> /* For memset(), memcpy(). */

#include <netcdf.h> /* For nc_close(). */

#include <Utilities.h> /* For PRE*(), NEW_ZERO(), Integer, Real, Stream. */

#include <Helpers.h>         /* For Name, timeData(). */
#include <M3IO.h>            /* For writeM3IOHeader(), writeM3IOData(). */
#include <NetCDFUtilities.h> /* For createNetCDFFile(). */
#include <Parameters.h>      /* For Parameters. */

/*================================= MACROS ==================================*/

/* Write regridded cell counts to final regridded output? */

#define OUTPUT_CELL_COUNTS 0
enum { OUTPUT_REGRID_VARIABLES = 6 }; /* lon,lat,col,row,counts,data. */

/*================================== TYPES ==================================*/

typedef struct {
  Integer      variables;       /* 3 = lon,lat,aod or 4 = lon,lat,time,aod. */
  Integer      timesteps;       /* E.g., 24. */
  Integer      scans;           /* E.g., 35 half-hour daylight scans. */
  Integer      totalPoints;     /* Sum of points[ scan ]. */
  Integer      maximumPoints;   /* Number of scan points in a hour. */
  Real         domain[ 2 ][ 2]; /*domain[LONGITUDE LATITUDE][MINIMUM MAXIMUM]*/
  UTCTimestamp timestamp;
  Line         note;            /* File note/description. */
  Name*        variable;        /* variable[ variables ]. E.g., "aod". */
  Name*        units;           /* units[ variables ]. E.g., "-". */
  Integer*     timestamps;      /* timestamps[ scans ] yyyydddhhmm. */
  Integer*     points;          /* points[ scans ] */
  Real*        data;            /* data[ variable ][ points_s ] */
  /* Regrid data: */
  Integer      totalRegriddedPoints;
  Integer*     outputPoints;    /* outputPoints[ timesteps ]. */
  Real*        longitudes;      /*longitudes[MIN(scans,timesteps)*maxPoints]*/
  Real*        latitudes;       /*latitudes[ MIN(scans,timesteps)*maxPoints]*/
  Real*        gridLongitudes;  /* gridLongitudes[ totalRegriddedPoints ]. */
  Real*        gridLatitudes;   /* gridLatitudes[ totalRegriddedPoints ]. */
  Integer*     columns;         /* columns[ totalRegriddedPoints ]. */
  Integer*     rows;            /* rows[ totalRegriddedPoints ]. */
  Integer*     counts;          /* counts[ totalRegriddedPoints ]. */
  Real*        weights;         /* weights[ totalRegriddedPoints ]. */
  Real*        gridData;        /* gridData[ totalRegriddedPoints ]. */
} Data;

typedef Integer (*Writer)( Data* data, const Parameters* parameters );

typedef struct {
  Integer format;         /* FORMAT_XDR, etc. */
  Writer writer;          /* Routine that writes data in this format. */
  Writer regriddedWriter; /* Routine that writes regridded data in format. */
} Entry;

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateData( Data* data );

static Integer isValidData( const Data* data );

static Writer dispatcher( Integer format, Integer regrid );

static Integer readXDR( Stream* input, Data* data );

static Integer readXDRData( Stream* input, Data* data );

static Integer readRegriddedXDR( Stream* input, Data* data );

static Integer readRegriddedXDRData( const Integer hasCounts,
                                     Stream* input, Data* data );

static Integer compareRegriddedXDR( const Parameters* parameters, Data* data );

static void countDataPoints( Data* data );

static Integer dataVariableIndex( const Data* data );

static Integer writeASCII( Data* data, const Parameters* parameters );

static void writeASCIIHeader( const Data* data, Stream* output );

static Integer writeASCIIData( Data* data, Stream* input, Stream* output );

static Integer writeCOARDS( Data* data, const Parameters* parameters );

static Integer writeCOARDSHeader( Integer file, const Data* data );

static Integer writeCOARDSData( Integer file, Stream* input,
                                const Data* data );

static Integer writeRegriddedXDR( Data* data, const Parameters* unused );

static Integer copyRegriddedXDRDataFromTempFile( Data* data,
                                                 const Parameters* parameters,
                                                 Stream* output );

static Integer writeRegriddedASCII( Data* data, const Parameters* unused );

static Integer writeRegriddedCOARDS( Data* data,
                                     const Parameters* parameters );

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Data* data );

static Integer writeRegriddedCOARDSData( Integer file, const Data* data,
                                         const Parameters* parameters );

static Integer writeBufferedRegriddedCOARDSData( Integer file,
                                                 Data* data,
                                                 const Parameters* parameters );

static Integer writeRegriddedIOAPI(Data* data, const Parameters* parameters);

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Data* data,
                                          const Grid* grid );

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        Data* data,
                                        const Parameters* parameters );

static void regridDataWithCorners( Parameters* const parameters, Data* data );

static void appendRegriddedData( Stream* const output,
                                 const size_t count,
                                 const size_t cellCounts[],
                                 const double cellData[],
                                 const double cellLongitudes[],
                                 const double cellLatitudes[],
                                 const size_t cellColumns[],
                                 const size_t cellRows[] );

static Integer copyRegriddedXDRDataFromTempFile( Data* data,
                                                 const Parameters* parameters,
                                                 Stream* output );

static Integer readRegriddedXDRDataFromTempFile( Stream* tempFile,
                                                 const Integer timestep,
                                                 Data* data );

static Integer seekToTimestep( const Integer timestep, const Data* data,
                               Stream* tempFile );

static void regridData( Stream* input, Integer method, Grid* grid, Data* data);

static Integer readScanDataForTimestamp( Integer yyyydddhh00,
                                         Stream* input,
                                         Data* data,
                                         Real longitudesSW[],
                                         Real longitudesSE[],
                                         Real longitudesNW[],
                                         Real longitudesNE[],
                                         Real latitudesSW[],
                                         Real latitudesSE[],
                                         Real latitudesNW[],
                                         Real latitudesNE[],
                                         Integer* points );

static Integer readScanData( Stream* input, Integer variables, Integer points,
                             Real longitudes[], Real latitudes[],
                             Real data[],
                             Real longitudesSW[],
                             Real longitudesSE[],
                             Real longitudesNW[],
                             Real longitudesNE[],
                             Real latitudesSW[],
                             Real latitudesSE[],
                             Real latitudesNW[],
                             Real latitudesNE[] );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateSwath - Read input and write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void translateSwath( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok,
         parameters->input->ok( parameters->input ) );

  Data data;
  ZERO_OBJECT( &data );
  parameters->ok = 0;

  if ( readXDR( parameters->input, &data ) ) {
    Writer writer = dispatcher( parameters->format, parameters->regrid );

    if ( ! writer ) {
      failureMessage( "Invalid/unsupported format/regrid specification." );
    } else if ( parameters->regrid ) {
      const Integer hasCorners = IN3( data.variables, 11, 12 );
      DEBUG( { time_t s = time(0);
               fprintf( stderr, "%sregridData()...\n", ctime(&s) ); } )
      if ( hasCorners ) { /* Has corners so regrid and aggregate. */
        regridDataWithCorners( parameters, &data ); /* Also aggregates. */
      } else {
        regridData( parameters->input, parameters->regrid, parameters->grid,
                    &data );
      }

      if ( data.totalRegriddedPoints == 0 ) {
        failureMessage( "No points projected onto the grid." );
      } else {

        if ( parameters->aggregationTimesteps ) {
          const Integer dataVariable = dataVariableIndex( &data );

          if ( ! hasCorners ) { /* Not aggregated already. */
            Integer totalOutputPoints = 0;

            DEBUG( time_t s = time(0); int unused_ =
              fprintf( stderr, "%s aggregate( aggregationTimesteps = %lld )\n",
                       ctime(&s), parameters->aggregationTimesteps ); )

            const Integer aggregatedTimesteps =
              aggregateData( parameters->aggregationTimesteps,
                             0,
                             data.timesteps,
                             data.outputPoints,
                             data.gridLongitudes,
                             data.gridLatitudes,
                             0,
                             data.columns,
                             data.rows,
                             0,
                             data.gridData,
                             0,
                             &totalOutputPoints );

            data.timesteps = aggregatedTimesteps;
            data.totalRegriddedPoints = totalOutputPoints;

            DEBUG( fprintf( stderr,
                            "timesteps = %lld, totalRegriddedPoints = %lld\n",
                            data.timesteps, data.totalRegriddedPoints ); )
          }

          if ( AND2( parameters->aggregationTimesteps == 24,
                     ! OR2( strstr( data.variable[ dataVariable ], "daily" ),
                            strstr( data.variable[ dataVariable], "DAILY")))) {
            Name dailyName = "";
            memset( dailyName, 0, sizeof dailyName );
            snprintf( dailyName, sizeof dailyName / sizeof *dailyName,
                      "daily_%s",
                      data.variable[ dataVariable ] );
            strncpy( data.variable[ dataVariable ], dailyName,
                     sizeof dailyName / sizeof *dailyName );
          }
        }

        DEBUG( {time_t s=time(0);fprintf(stderr,"%swriter()...\n",ctime(&s));})
        parameters->ok = writer( &data, parameters );
      }
    } else {
      parameters->ok = writer( &data, parameters );
    }
  }

  deallocateData( &data );
  DEBUG( { time_t s = time(0); fprintf( stderr, "%sdone\n", ctime(&s) ); } )
  POST0( isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: compareRegriddedSwath - Read REGRIDDED-Swath input, compare it to
         CMAQ XDR data and hasCornerse it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedSwath( Parameters* parameters ) {
    
  PRE03( isValidParameters( parameters ),
         parameters->ok, parameters->input->ok( parameters->input ) );

  if ( ! AND3( ! parameters->regrid,
               OR2( parameters->compareFunction, parameters->convertFunction ),
               parameters->data ) ) {
    failureMessage( "Invalid input for comparing." );
    parameters->ok = 0;
  } else {
    Data data;
    ZERO_OBJECT( &data );
    parameters->ok = 0;

    if ( readRegriddedXDR( parameters->input, &data ) ) {
      compareFunctionNameUnits( parameters->compareFunction,
                                parameters->convertFunction,
                                data.variable[ 0 ], data.units[ 0 ], 
                                parameters->variable, 
                                parameters->units );

      if ( compareRegriddedXDR( parameters, &data ) ) {
        Writer writer = dispatcher( parameters->format, 1 );
        CHECK( writer );

        if ( data.totalRegriddedPoints == 0 ) {
          failureMessage( "No points projected onto the grid." );
        } else {
          parameters->ok = writer( &data, parameters );
        }
      }
    }

    deallocateData( &data );
  }

  POST0( isValidParameters( parameters ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateData - Deallocate contents of data structure.
INPUTS:  Data* data Structure to deallocate contents of.
******************************************************************************/

static void deallocateData( Data* data ) {
  PRE0( data );

  if ( data->data == 0 ) { /* Read Regridded XDR compare case. */
    FREE( data->gridLongitudes );
  }

  FREE( data->outputPoints );
  FREE( data->variable );
  FREE( data->timestamps );
  FREE( data->data );
  ZERO_OBJECT( data );
}



/******************************************************************************
PURPOSE: isValidData - Check data structure.
INPUTS:  const Data* data Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidData( const Data* data ) {
  const Integer result =
    AND15( data,
           data->note[ 0 ],
           isValidUTCTimestamp( data->timestamp ),
           data->timesteps > 0,
           data->variables > 0,
           isValidLongitudeLatitude( data->domain[ LONGITUDE ][ MINIMUM ],
                                     data->domain[ LATITUDE  ][ MINIMUM ] ),
           isValidLongitudeLatitude( data->domain[ LONGITUDE ][ MAXIMUM ],
                                     data->domain[ LATITUDE  ][ MAXIMUM ] ),
           data->domain[ LONGITUDE ][ MINIMUM ] <=
             data->domain[ LONGITUDE ][ MAXIMUM ],
           data->domain[ LATITUDE  ][ MINIMUM ] <=
             data->domain[ LATITUDE  ][ MAXIMUM ],
           data->variable,
           data->units,
           data->variable[ 0 ],
           data->units[ 0 ],
           IMPLIES( GT_ZERO2( data->scans, data->totalPoints ),
             AND10( data->maximumPoints >=
                      maximumItemI( data->points, data->scans ),
                    data->timestamps,
                    isValidTimestamp( data->timestamps[ 0 ] ),
                    isValidTimestamp( data->timestamps[ data->scans - 1 ] ),
                    data->timestamps[ data->scans - 1 ] >= data->timestamps[0],
                    data->points,
                    data->points[ 0 ] > 0,
                    data->points[ data->scans - 1 ] > 0,
                    data->data,
                    data->totalRegriddedPoints >= 0 ) ),
          IMPLIES( data->totalRegriddedPoints > 0,
                   AND7( data->outputPoints,
                         minimumItemI( data->outputPoints, data->timesteps )
                           >= 0,
                         data->columns,
                         data->rows,
                         data->gridLongitudes,
                         data->gridLatitudes,
                         data->gridData ) ) );

  return result;
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
PURPOSE: readXDR - Read input and initialize data structure.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Data* data Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDR( Stream* input, Data* data ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         data, data->variable == 0, data->data == 0 );

  Integer result = 0;

  input->readString( input, data->note, COUNT( data->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( data->note );

    if ( readTimestamp( input, data->timestamp ) ) {
      Integer dimensions[ 3 ] = { 0, 0, 0 };

      if ( readDimensions( input, COUNT( dimensions ), dimensions ) ) {
        data->variables = dimensions[ 0 ];
        data->timesteps = dimensions[ 1 ];
        data->scans     = dimensions[ 2 ];
        data->variable  = NEW_ZERO( Name, data->variables * 2 );

        if ( data->variable ) {
          data->units = data->variable + data->variables;

          if ( readVariablesAndUnits( input, data->variables,
                                      data->variable, data->units ) ) {

            if ( readDomain( input, data->domain ) ) {

              if ( skipInputLines( input, 3 ) ) {
                result = readXDRData( input, data );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Data data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read initial binary data from input.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Data* data Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDRData( Stream* input, Data* data ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         data, IN5( data->variables, 3, 4, 11, 12 ), data->scans > 0,
         data->timestamps == 0, data->data == 0 );

  Integer result = 0;
  data->timestamps = NEW_ZERO( Integer, data->scans * 2 );

  if ( data->timestamps ) {
    data->points = data->timestamps + data->scans;
    input->read64BitIntegers( input, data->timestamps, data->scans );

    if ( AND3( input->ok( input ),
               isValidTimestamp( data->timestamps[ 0 ] ),
               isValidTimestamp( data->timestamps[ data->scans - 1 ] ) ) ) {

      input->read64BitIntegers( input, data->points, data->scans );

      /* Sum rows x columns per scan: */

      if ( input->ok( input ) ) {
        countDataPoints( data );

        if ( data->totalPoints > 0 ) {
          data->data = NEW_ZERO( Real, data->maximumPoints );
          result = isValidData( data );
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Data data." );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDR - Read REGRIDDED-Swath and initialize data.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Data* data Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readRegriddedXDR( Stream* input, Data* data ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         data, data->variable == 0, data->data == 0, data->gridData == 0 );

  Integer result = 0;

  input->readString( input, data->note, COUNT( data->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( data->note );

    if ( readTimestamp( input, data->timestamp ) ) {

      if ( readDimensions( input, 1, &data->timesteps ) ) {
        data->timestamps = NEW_ZERO( Integer, data->timesteps );

        if ( data->timestamps ) {
          Integer timestamp = fromUTCTimestamp( data->timestamp );
          Integer timestep = 0;

          for ( timestep = 0; timestep < data->timesteps; ++timestep ) {
            data->timestamps[ timestep ] = timestamp;
            incrementTimestamp( &timestamp );
          }

          data->variables = 1;
          data->variable = NEW_ZERO( Name, 2 );

          if ( data->variable ) {
            data->units = data->variable + 1;

            if ( readVariablesAndUnits( input, data->variables,
                                        data->variable, data->units ) ) {
              Integer hasCounts = 0;
              char line[ 256 ] = "";
              memset( line, 0, sizeof line );

              do {
                input->readString(input, line, sizeof line / sizeof *line - 1);
              } while ( AND2( input->ok( input ),
                              strcmp( line,
                      "# MSB 64-bit integers rows[timesteps][points] and\n")));

              if ( input->ok( input) ) {
                input->readString(input, line, sizeof line / sizeof *line - 1);

                if ( input->ok( input) ) {

                  if ( ! strcmp( line,
                    "# MSB 64-bit integers counts[timesteps][points] and\n")) {
                    hasCounts = 1;
                    input->readString( input, line,
                                       sizeof line / sizeof *line - 1 );
                  }

                  if ( ! strcmp( line,
                      "# IEEE-754 64-bit reals data[timesteps][points]:\n" )) {
                    result = readRegriddedXDRData( hasCounts, input, data );
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Data data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRData - Read regridded data from input.
INPUTS:  const Integer hasCounts   Read counts?
         Stream* input             Input stream read from and parse.
OUTPUTS: Data* data Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readRegriddedXDRData(  const Integer hasCounts,
                                      Stream* input, Data* data ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         data, data->timesteps > 0, data->variables > 0,
         data->scans == 0, data->data == 0 );

  Integer result = 0;

  data->outputPoints = NEW_ZERO( Integer, data->timesteps );

  if ( data->outputPoints ) {
    input->read64BitIntegers( input, data->outputPoints, data->timesteps );

    if ( input->ok( input ) ) {
      const Integer count = data->totalRegriddedPoints =
        sumI( data->outputPoints, data->timesteps );

      if ( count > 0 ) {
        data->gridLongitudes = NEW_ZERO( Real, count * OUTPUT_REGRID_VARIABLES);

        if ( data->outputPoints ) {
          data->gridLatitudes = data->gridLongitudes + count;
          data->columns       = (Integer*) data->gridLatitudes + count;
          data->rows          = data->columns + count;
          data->counts        = data->rows + count;
          data->gridData      = (Real*) data->counts + count;

          input->read64BitReals( input, data->gridLongitudes, count * 2 );

          if ( input->ok( input ) ) {
            input->read64BitIntegers( input, data->columns,
                                      count * ( 2 + hasCounts ) );

            if ( input->ok( input ) ) {
              input->read64BitReals( input, data->gridData, count );

              if ( input->ok( input ) ) {
                result = isValidData( data );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Data data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: compareRegriddedXDR - Compare Regridded data with CMAQ data.
INPUTS:  const Parameters* parameters  CMAQ data to compare to.
OUTPUTS: Data* data                    Updated data->data.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    Data* data ) {

  PRE05( parameters, isValidParameters( parameters ),
         XOR2( parameters->compareFunction, parameters->convertFunction ),
         data, isValidData( data ) );

  const Integer isDaily = parameters->timesteps * 24 == data->timesteps;
  Integer result = 0;

  if ( ! OR2( isDaily,
              AND2( ! strcmp( parameters->timestamp, data->timestamp ),
                    parameters->timesteps == data->timesteps ) ) ) {
    failureMessage( "Mismatched time steps (%s %lld)"
                    " for comparison to CMAQ data (%s %lld).",
                    data->timestamp, data->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else {
    const Integer timesteps = data->timesteps;
    Real* const regriddedData        = data->gridData;
    const Integer* const dataRows    = data->rows;
    const Integer* const dataColumns = data->columns;
    const Integer* const dataPoints  = data->outputPoints;
    const Real* const cmaqData = parameters->data;
    const Real* const cmaqData2 = parameters->data2;
    CompareFunction comparer   = parameters->compareFunction;
    ConvertFunction converter  = parameters->convertFunction;
    const Integer firstRow     = parameters->firstRow;
    const Integer lastRow      = parameters->lastRow;
    const Integer firstColumn  = parameters->firstColumn;
    const Integer lastColumn   = parameters->lastColumn;
    const Integer rows         = lastRow    - firstRow    + 1;
    const Integer columns      = lastColumn - firstColumn + 1;
    const Integer rowsTimesColumns = rows * columns;
    Integer timestep = 0;
    Integer cmaqTimestep = 0;
    Integer dataIndex = 0;

    DEBUG( fprintf( stderr, "timesteps = %lld, rows = %lld, columns = %lld\n",
                    timesteps, rows, columns ); )

    for ( timestep = 0; timestep < timesteps;
          ++timestep, cmaqTimestep = ( isDaily ? timestep / 24 : timestep ) ) {
      const Integer points = dataPoints[ timestep ];
      const Integer cmaqTimestepOffset = cmaqTimestep * rowsTimesColumns;
      Integer point = 0;

      DEBUG( fprintf( stderr, "timestep = %lld, "
                      "cmaqTimestep = %lld, cmaqTimestepOffset = %lld, "
                      "points = %lld\n",
                      timestep, cmaqTimestep, cmaqTimestepOffset, points ); )


      for ( point = 0; point < points; ++point, ++dataIndex ) {
        const Integer dataRow    = dataRows[    dataIndex ];
        const Integer dataColumn = dataColumns[ dataIndex ];
        DEBUG( fprintf( stderr, "  @[ %lld, %lld ]: ", dataRow, dataColumn);)
        CHECK( IN_RANGE( dataIndex, 0, data->totalRegriddedPoints - 1 ) );

        if ( AND2( IN_RANGE( dataRow, firstRow, lastRow ),
                   IN_RANGE( dataColumn, firstColumn, lastColumn ) ) ) {
          const Integer cmaqRow0    = dataRow    - firstRow;
          const Integer cmaqColumn0 = dataColumn - firstColumn;
          const Integer cmaqIndex =
            cmaqTimestepOffset + cmaqRow0 * columns + cmaqColumn0;
          const Real dataDatum = regriddedData[ dataIndex ];
          const Real cmaqDatum = cmaqData[ cmaqIndex ];
          const Real newDatum =
            comparer ? comparer( dataDatum, cmaqDatum )
            : converter( dataDatum, cmaqDatum, cmaqData2[ cmaqIndex ] );
          CHECK3( IN_RANGE( cmaqRow0, 0, rows - 1 ),
                  IN_RANGE( cmaqColumn0, 0, columns - 1 ),
                  IN_RANGE( cmaqIndex, 0, timesteps * rows * columns - 1 ) );
          regriddedData[ dataIndex ] = newDatum;
          result = 1;
          DEBUG( fprintf( stderr, "@[dataIndex = %lld, cmaqIndex = %lld, "
                          "cmaqRow0 = %lld, cmaqColumn0 = %lld] ",
                          dataIndex, cmaqIndex, cmaqRow0, cmaqColumn0 ); )
          DEBUG( if ( comparer ) fprintf( stderr, "f(%lf, %lf) -> %lf\n",
                                          dataDatum, cmaqDatum, newDatum );
                else fprintf( stderr, "f(%lf, %lf, %lf) -> %lf\n",
                              dataDatum, cmaqDatum, cmaqData2[ cmaqIndex ],
                              newDatum ); )
        } else {
          regriddedData[ dataIndex ] = -9999.0;
          DEBUG( fprintf( stderr, "-9999\n" ); )
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "No points in output." );
  }

  POST02( IS_BOOL( result ), isValidData( data ) );
  return result;
}



/******************************************************************************
PURPOSE: countDataPoints - Compute sum of all scan points and also the maximum
         number of scan points in an hour.
INPUTS:  Data* data  Partially initialized structure.
OUTPUTS: Data* data  totalPoints = Sum of all scan points.
                     maximumPoints = Maximum number of scan points in an hour.
******************************************************************************/

static void countDataPoints( Data* data ) {

  PRE05( data, data->scans > 0, data->points, data->totalPoints == 0,
         data->maximumPoints == 0 );

  const Integer scans = data->scans;
  Integer scan = 0;

  data->totalPoints = 0;
  data->maximumPoints = 0;

  do {
    const Integer scanPoints = data->points[ scan ];
    const Integer yyyydddhhmm = data->timestamps[ scan ];
    const Integer yyyydddhh = yyyydddhhmm / 100;
    Integer otherScan = 0;
    Integer scanPointsInHour = scanPoints;

    for ( otherScan = scan + 1; otherScan < scans; ++otherScan ) {
      const Integer yyyydddhhmm2 = data->timestamps[ otherScan ];
      const Integer yyyydddhh2 = yyyydddhhmm2 / 100;

      if ( yyyydddhh2 == yyyydddhh ) {
        const Integer scanPoints2 = data->points[ otherScan ];
        scanPointsInHour += scanPoints2;
      }
    }

    data->totalPoints += scanPoints;
    data->maximumPoints = MAX( data->maximumPoints, scanPointsInHour );

    ++scan;
  } while ( scan < scans );

  POST02( data->totalPoints >= 0,
          IMPLIES_ELSE( data->totalPoints > 0,
                        AND2( data->maximumPoints > 0,
                              data->maximumPoints >=
                                maximumItemI( data->points, data->scans ) ),
                        data->maximumPoints == 0 ) );
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format output.
INPUTS:  Data* data       Structure to write.
         const Parameters*  Parameters input stream to read from.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeASCII( Data* data, const Parameters* parameters ) {

  PRE03( isValidData( data ),
         isValidParameters( parameters ),
         parameters->input->ok( parameters->input ) );

  Integer result = 0;

  /*
   * First reallocate data->data to be large enough to hold
   * all variables of the largest scan:
   */

  FREE( data->data );
  data->data = NEW_ZERO( Real, data->maximumPoints * data->variables );

  if ( data->data ) {
    Stream* output = newFileStream( "-stdout", "wb" );

    if ( output ) {
      writeASCIIHeader( data, output );

      if ( output->ok( output ) ) {
        result = writeASCIIData( data, parameters->input, output );
      }

      FREE_OBJECT( output );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCIIHeader - Write ASCII-format header line.
INPUTS:  const Data* data  Structure to write.
         Stream* output      Stream to write to.
******************************************************************************/

static void writeASCIIHeader( const Data* data, Stream* output ) {

  PRE03( isValidData( data ), output, output->isWritable( output ) );

  const char* const headerStart =
    "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)";
  const char* const headerFormat = "\t%s(%s)";

  /* Write header row: */

  output->writeString( output, headerStart );

  if ( output->ok( output ) ) {
    const Integer variables = data->variables;
    Integer variable = 2;

    do {
      output->writeString( output, headerFormat,
                           data->variable[ variable ],
                           data->units[ variable ] );

      if ( ! output->ok( output ) ) {
        variable = variables;
      }

      ++variable;
    } while ( variable < variables );

    if ( output->ok( output ) ) {
      output->writeString( output, "\n" );
    }
  }
}



/******************************************************************************
PURPOSE: writeASCIIData - Write ASCII-format data lines.
INPUTS:  const Data* data  Structure to write.
         Stream* input       Stream to read from.
         Stream* output      Stream to write to.
RETURNS: Integer 1 if successful, else 0 and failureMesage is called.
******************************************************************************/

static Integer writeASCIIData( Data* data, Stream* input, Stream* output ) {

  PRE05( isValidData( data ), input, input->isReadable( input ),
         output, output->isWritable( output ) );

  Integer result = 0;
  const Integer variables = data->variables;
  const char* const dataFormat = "\t%28.12e";
  const Integer scans = data->scans;
  Integer scan = 0;

  /* Write data rows: */

  do {
    const Integer scanPoints = data->points[ scan ];
    const Integer scanSize = variables * scanPoints;

    input->read64BitReals( input, data->data, scanSize );

    if ( ! input->ok( input ) ) {
      scan = scans;
    } else {
      Real* scanData = data->data;
      Integer point = 0;
      const Integer timestamp = data->timestamps[ scan ];
      UTCTimestamp timestampString;
      toUTCTimestamp( timestamp, timestampString );

      do {
        output->writeString( output, timestampString ); /* Begin row. */

        if ( output->ok( output ) ) {
          Integer variable = 0;

          do {
            const Real datum = *( scanData + variable * scanPoints );
            output->writeString( output, dataFormat, datum );

            if ( ! output->ok( output ) ) {
              variable = variables;
              point = scanPoints;
              scan = scans;
            }

            ++variable;
          } while ( variable < variables );

          if ( output->ok( output ) ) {
            output->writeString( output, "\n" ); /* End of row. */
          }
        }

        ++scanData;
        ++point;
      } while ( point < scanPoints );
    }

    ++scan;
  } while ( scan < scans );

  result = AND2( input->ok( input ), output->ok( output ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDS - Write COARDS-format data.
INPUTS:  Data* data                  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDS( Data* data, const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
      data->variables * data->totalPoints * 4 + /* variables(points). */
      data->totalPoints * 3 * 4 + 2000;
      /* yyyyddd(points),hhmmss(points),time(points) + header/extra. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file = createNetCDFFile(parameters->netcdfFileName, create64BitFile);

  if ( file != -1 ) {

    if ( writeCOARDSHeader( file, data ) ) {
      result = writeCOARDSData( file, parameters->input, data );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSHeader - Write header to file.
INPUTS:  Integer file      NetCDF file to write to.
         const Data* data  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSHeader( Integer file, const Data* data ) {

  PRE02( file != -1, isValidData( data ) );

  Integer result = 0;
  const char* const name = "points";
  Integer dimensionId = -1;

  if ( createDimensions( file, 1, &name, &data->totalPoints, &dimensionId ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {
        const Integer variables = data->variables;
        Integer index = 2;

        do {
          const char* const units =
            strcmp( data->units[ index ], "-"   ) == 0 ? "none" :
            strcmp( data->units[ index ], "deg" ) == 0 ? "degrees" :
            data->units[ index ];
          const Integer ok =
            createVariable( file, data->variable[ index ], units,
                            NC_FLOAT, 1, 1, &dimensionId ) != -1;

           if ( ! ok ) {
             index = variables;
           }

           ++index;
        } while ( index < variables );

        if ( index == variables ) {

          if ( writeExtraAttributes( file, (const Real (*)[2]) data->domain,
                                     dimensionId ) ) {
            UTCTimestamp timestamp = "";
            Line history = "";
            appendToLine( history, data->note );
            appendToLine( history, ",XDRConvert" );
            toUTCTimestamp( data->timestamps[ 0 ], timestamp );
            result = writeStandardContents( file, history,
                                            timestamp,
                                            dimensionId, data->totalPoints, 0);
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSData - Write COARDS-format data to file.
INPUTS:  Integer file        NetCDF file to write to.
         Stream* input       Stream to read from.
         const Data* data  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSData(Integer file, Stream* input, const Data* data) {

  PRE05( input, input->ok( input ), input->isReadable( input ),
         file != -1, isValidData( data ) );

  Integer result = 0;
  const Integer scans     = data->scans;
  const Integer variables = data->variables;
  Integer scan = 0;
  Integer variable = 0;
  Integer start = 0;

  do {
    const Integer count = data->points[ scan ];
    variable = 0;

    do {
      const char* const variableName =
        variable == 0 ? "longitude" :
        variable == 1 ? "latitude" :
        data->variable[ variable ];

      input->read64BitReals( input, data->data, count );

      if ( ! AND2( input->ok( input ),
                   writeSomeData( file, variableName, start,
                                  count, 1, 1, 1, data->data ) ) ) {
        scan = scans;
        variable = variables;
      }

      ++variable;
    } while ( variable < variables );

    start += count;
    ++scan;
  } while ( scan < scans );

  result = AND3( scan == scans, variable == variables,
                 writeTimeData( file, data->scans, 1, 0,
                                data->timestamps, data->points, data->data ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: dataVariableIndex - 0-based index of data variable (not coordinates or
         time).
INPUTS:  const Data* data  Structure to check.
RETURNS: Integer 1 index.
******************************************************************************/

static Integer dataVariableIndex( const Data* data ) {
  PRE0( isValidData( data ) );

  const Integer hasCorners =
    ! strcmp( data->variable[ data->variables - 1 ], "Latitude_NE" );

  const Integer result =
    hasCorners ? data->variables - 9 : data->variables - 1;

  POST05( IN_RANGE( result, 0, data->variables - 1 ),
          ! strstr( data->variable[ result ], "Longitude" ),
          ! strstr( data->variable[ result ], "Latitude" ),
          strcmp( data->variable[ result ], "Scan_Start_Time" ),
          strcmp( data->variable[ result ], "Profile_UTC_Time" ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  Data* data          Structure to write.
         const Parameters*isValidData(
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedXDR( Data* data, const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = data->timesteps;
    const Integer points = data->totalRegriddedPoints;
    const Integer index = dataVariableIndex( data );
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps
      : 1;
    Name variable = "";
    aggregateName( data->variable[ index ], hoursPerTimestep, variable );

    output->writeString( output,
                         "REGRIDDED-Swath 2.0\n"
                         "%s,XDRConvert\n"
                         "%s\n"
                         "# timesteps\n%"INTEGER_FORMAT"\n"
                         "# Variable name:\n%s\n"
                         "# Variable units:\n%s\n",
                         data->note,
                         data->timestamp,
                         timesteps,
                         variable,
                         data->units[ index ] );
  
    writeProjectionAndGrid( parameters->grid, output );

#if OUTPUT_CELL_COUNTS
    output->writeString( output,
                 "# MSB 64-bit integers points[timesteps] and\n"
                 "# IEEE-754 64-bit reals longitudes[timesteps][points] and\n"
                 "# IEEE-754 64-bit reals latitudes[timesteps][points] and\n"
                 "# MSB 64-bit integers columns[timesteps][points] and\n"
                 "# MSB 64-bit integers rows[timesteps][points] and\n"
                 "# MSB 64-bit integers counts[timesteps][points] and\n"
                 "# IEEE-754 64-bit reals data[timesteps][points]:\n" );
#else
    output->writeString( output,
                 "# MSB 64-bit integers points[timesteps] and\n"
                 "# IEEE-754 64-bit reals longitudes[timesteps][points] and\n"
                 "# IEEE-754 64-bit reals latitudes[timesteps][points] and\n"
                 "# MSB 64-bit integers columns[timesteps][points] and\n"
                 "# MSB 64-bit integers rows[timesteps][points] and\n"
                 "# IEEE-754 64-bit reals data[timesteps][points]:\n" );
#endif
    if ( output->ok( output ) ) {
      const Integer hasCorners = IN3( data->variables, 11, 12 );
      output->write64BitIntegers( output, data->outputPoints, timesteps );

      if ( hasCorners ) {
        result = copyRegriddedXDRDataFromTempFile( data, parameters, output );
      } else {

        if ( output->ok( output ) ) {
          output->write64BitReals( output, data->gridLongitudes, points );

          if ( output->ok( output ) ) {
            output->write64BitReals( output, data->gridLatitudes, points );

            if ( output->ok( output ) ) {
              output->write64BitIntegers( output, data->columns, points);

              if ( output->ok( output ) ) {
                output->write64BitIntegers( output, data->rows, points );
#if OUTPUT_CELL_COUNTS
                if ( output->ok( output ) ) {
                  output->write64BitIntegers( output, data->counts, points );
                }
#endif
                if ( output->ok( output ) ) {
                  output->write64BitReals( output, data->gridData, points );
                  result = output->ok( output );
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
PURPOSE: copyRegriddedXDRDataFromTempFile - Copy regridded XDR-format data from
         tempFile to output.
INPUTS:  Data* data                    Structure to write.
         const Parameters* parameters  parameters->regridFileName.
         Stream* output                Stream to write to.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer copyRegriddedXDRDataFromTempFile( Data* data,
                                                 const Parameters* parameters,
                                                 Stream* output ) {

  PRE05( isValidData( data ),
         isValidParameters( parameters ),
         output,
         output->ok( output ),
         IN3( data->variables, 11, 12 ) );

  Integer result = 0;

  /*
   * Copy each block of data from tempFile:
   * (with timesteps = 1, points = outputPoints[ timestep ] ):
   * # IEEE-754 64-bit reals longitudes[timesteps][points]
   * # IEEE-754 64-bit reals latitudes[timesteps][points]
   * # MSB 64-bit integers columns[timesteps][points]
   * # MSB 64-bit integers rows[timesteps][points]
   * # MSB 64-bit integers counts[timesteps][points]
   * # IEEE-754 64-bit reals data[timesteps][points]
   */

  Stream* tempFile = newFileStream( parameters->regridFileName, "rb" );

  if ( tempFile ) {
    void* buffer = data->gridData;
    const Integer tempFileVariables = OUTPUT_REGRID_VARIABLES;
    /* lon,lat,col,row,counts,data. */
    Integer variable = 0;

    for ( variable = 0;
          AND3( tempFile->ok( tempFile ), output->ok( output ),
                variable < tempFileVariables ); ++variable ) {
      const Integer timesteps = data->timesteps;
      Integer timestep = 0;
      tempFile->seekFromStart( tempFile, 0 ); /* Rewind for each variable. */

      for ( timestep = 0;
            AND3( tempFile->ok( tempFile ), output->ok( output ),
                  timestep < timesteps ); ++timestep ) {
        const Integer points = data->outputPoints[ timestep ];

        if ( points > 0 ) {
          const Integer bytes = points * 8;
          CHECK( bytes > 0 );
          const Integer seekBytes = variable * bytes;
          CHECK( seekBytes >= 0 );
          tempFile->seekFromCurrent( tempFile, seekBytes ); /* Seek variable */

          if ( tempFile->ok( tempFile ) ) {
            tempFile->readBytes( tempFile, buffer, bytes );

            if ( tempFile->ok( tempFile ) ) {
#if OUTPUT_CELL_COUNTS
              output->writeBytes( output, buffer, bytes );
#else
              if ( variable != 4 ) {
                output->writeBytes( output, buffer, bytes );
              }
#endif

              /* Skip over remaining variables for this timestep: */

              if ( output->ok( output ) ) {
                const Integer skipBytes =
                  ( tempFileVariables - variable - 1 ) * bytes;
                tempFile->seekFromCurrent( tempFile, skipBytes );
              }
            }
          }
        }
      }
    }

    result = AND2( tempFile->ok( tempFile ), output->ok( output ) );
    FREE_OBJECT( tempFile );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRDataFromTempFile - Read regridded XDR-format data from
         tempFile into data structure.
INPUTS:  Stream* tempFile         Temp regrid file to read from.
         const Integer timestep   Timestep to read.
OUTPUTS: Data* data               data->gridLongitudes,
                                  data->gridLatitudes,
                                  data->columns,
                                  data->rows,
                                  data->counts,
                                  data->gridData
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readRegriddedXDRDataFromTempFile( Stream* tempFile,
                                                 const Integer timestep,
                                                 Data* data ) {

  PRE012( tempFile,
          tempFile->invariant( tempFile ),
          tempFile->isReadable( tempFile ),
          tempFile->isSeekable( tempFile ),
          tempFile->ok( tempFile ),
          isValidData( data ),
          data->totalRegriddedPoints > 0,
          data->data,
          IN3( data->variables, 11, 12 ),
          IN_RANGE( timestep, 0, data->timesteps - 1 ),
          data->outputPoints,
          data->outputPoints[ timestep ] > 0 );

  Integer result = 0;

  /*
   * Read a block of data from temp regrid file for the given timestep:
   * (with timesteps = 1, points = outputPoints[ timestep ] ):
   * # IEEE-754 64-bit reals longitudes[timesteps][points]
   * # IEEE-754 64-bit reals latitudes[timesteps][points]
   * # MSB 64-bit integers columns[timesteps][points]
   * # MSB 64-bit integers rows[timesteps][points]
   * # MSB 64-bit integers counts[timesteps][points]
   * # IEEE-754 64-bit reals data[timesteps][points]
   */

  if ( seekToTimestep( timestep, data, tempFile ) ) {
    void* const pointers[ OUTPUT_REGRID_VARIABLES ] = {
      data->gridLongitudes,
      data->gridLatitudes,
      data->columns,
      data->rows,
      data->counts,
      data->gridData
    };
    const Integer points = data->outputPoints[ timestep ];
    const Integer bytes = points * 8; /* All data values are 64-bit. */
    Integer variable = 0;
    CHECK( sizeof pointers / sizeof *pointers == OUTPUT_REGRID_VARIABLES );

    for ( variable = 0;
          AND2( tempFile->ok( tempFile ), variable < OUTPUT_REGRID_VARIABLES );
          ++variable ) {
      void* const buffer = pointers[ variable ];
      CHECK( buffer );
      tempFile->readBytes( tempFile, buffer, bytes );

      if ( tempFile->ok( tempFile ) ) {
        rotate8ByteArrayIfLittleEndian( buffer, points );
      }
    }
  }

  result = tempFile->ok( tempFile );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: seekToTimestep - Seek to the given timestep in the temp regrid file.
INPUTS:  const Integer timestep  Timestep to read.
         const Data* data        data->outputPoints[ data->timesteps ].
OUTPUTS: Stream* tempFile        Temp regrid file.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer seekToTimestep( const Integer timestep, const Data* data,
                               Stream* tempFile ) {

  PRE011( data,
          isValidData( data ),
          data->totalRegriddedPoints > 0,
          data->data,
          IN3( data->variables, 11, 12 ),
          IN_RANGE( timestep, 0, data->timesteps - 1 ),
          tempFile,
          tempFile->invariant( tempFile ),
          tempFile->isReadable( tempFile ),
          tempFile->isSeekable( tempFile ),
          tempFile->ok( tempFile ) );

  /*
   * Seek to a block of data from temp regrid file for the given timestep:
   * (with timesteps = 1, points = outputPoints[ timestep ] ):
   * # IEEE-754 64-bit reals longitudes[timesteps][points]
   * # IEEE-754 64-bit reals latitudes[timesteps][points]
   * # MSB 64-bit integers columns[timesteps][points]
   * # MSB 64-bit integers rows[timesteps][points]
   * # MSB 64-bit integers counts[timesteps][points]
   * # IEEE-754 64-bit reals data[timesteps][points]
   */

  Integer result = 0;
  Integer pointSum = 0;
  Integer thisTimestep = 0;

  for ( thisTimestep = 0; thisTimestep < timestep; ++thisTimestep ) {
    const Integer points = data->outputPoints[ thisTimestep ];
    pointSum += points;
  }

  {
    const Integer bytes = pointSum * OUTPUT_REGRID_VARIABLES * 8; /* 64-bit */
    tempFile->seekFromStart( tempFile, bytes );
  }

  result = tempFile->ok( tempFile );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedASCII - Write regridded ASCII-format data.
INPUTS:  Data* data         Structure to write.
         const Parameters*  Structure of parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedASCII( Data* data, const Parameters* parameters ) {

  PRE03( isValidData( data ), data->variables > 0,
         isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const char* const headerStart =
      "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)"
      "\tCOLUMN(-)\tROW(-)\tCOUNT(-)";
    const char* const headerFormat = "\t%s(%s)\n";
    const char* const dataFormat =
      "%s\t%10.5f\t%10.5f"
      "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%28.12e\n";
    const Integer hasCorners = IN3( data->variables, 11, 12 );
    Stream* tempFile =
      hasCorners ? newFileStream( parameters->regridFileName, "rb" ) : 0;
    Integer ok = IMPLIES( hasCorners, tempFile );

    /* Write header row: */

    if ( ok ) {
      output->writeString( output, headerStart );
      ok = output->ok( output );
    }

    if ( ok ) {
      const Integer index = dataVariableIndex( data );
      const Integer hoursPerTimestep =
        parameters->aggregationTimesteps ? parameters->aggregationTimesteps
        : 1;
      Name variable = "";
      aggregateName( data->variable[ index ], hoursPerTimestep, variable );
      output->writeString( output, headerFormat,
                           variable, data->units[ index ] );

      if ( output->ok( output ) ) {
        const Integer timesteps = data->timesteps;
        const Real* longitudes  = data->gridLongitudes;
        const Real* latitudes   = data->gridLatitudes;
        const Integer* columns  = data->columns;
        const Integer* rows     = data->rows;
        const Integer* counts   = data->counts;
        const Real* gridData    = data->gridData;
        const Integer hoursPerTimestep =
          parameters->aggregationTimesteps ? parameters->aggregationTimesteps
          : 1;
        Integer timestep = 0;
        Integer yyyydddhh00 =
          ( fromUTCTimestamp( data->timestamp ) / 100 ) * 100;
        UTCTimestamp timestamp;

        /* Write data rows: */

        do {
          const Integer points = data->outputPoints[ timestep ];

          if ( points > 0 ) {
            Integer point = 0;
            toUTCTimestamp( yyyydddhh00, timestamp );

            if ( hasCorners ) {
              ok = readRegriddedXDRDataFromTempFile( tempFile, timestep, data );
            }

            if ( ok ) {
              for ( point = 0; point < points; ++point ) {
                const Real longitude = longitudes[ point ];
                const Real latitude  = latitudes[ point ];
                const Integer column = columns[ point ];
                const Integer row    = rows[ point ];
                const Integer count  = counts[ point ];
                const Real value     = gridData[ point ];

                if ( count > 0 ) {
                  output->writeString( output, dataFormat,
                                       timestamp, longitude, latitude,
                                       column, row, count, value );

                  if ( ! output->ok( output ) ) {
                    point = points;
                    timestep = timesteps;
                  }
                }
              }
            }
          }

          yyyydddhh00 = offsetTimestamp( yyyydddhh00, hoursPerTimestep );
          ++timestep;
        } while ( AND2( ok, timestep < timesteps ) );
      }
    }

    result = AND2( ok, output->ok( output ) );
    FREE_OBJECT( tempFile );
    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDS - Write regridded COARDS-format data.
INPUTS:  Data* data                    Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDS( Data* data,
                                     const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    data->totalRegriddedPoints * ( OUTPUT_REGRID_VARIABLES + 2 ) * 4 + 10000;
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps
      : 1;

    if ( writeRegriddedCOARDSHeader( file, hoursPerTimestep, data ) ) {
      const Integer hasCorners = IN3( data->variables, 11, 12 );

      if ( hasCorners ) {
        result = writeBufferedRegriddedCOARDSData( file, data, parameters );
      } else {
        result = writeRegriddedCOARDSData( file, data, parameters );
      }
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDSHeader - Write header to file.
INPUTS:  Integer file              NetCDF file to write to.
         Integer hoursPerTimestep  E.g., 1 or 24.
         const Data* data          Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Data* data ) {

  PRE03( file != -1, hoursPerTimestep > 0, isValidData( data ) );

  const char* const dimensionName = "points";
  Integer dimensionId = -1;
  const Integer dimension = data->totalRegriddedPoints;
  Integer result = 0;

  if ( createDimensions( file, 1, &dimensionName, &dimension, &dimensionId ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "column", "-",  NC_INT, 0, 1, &dimensionId )
           != -1 ) {

        if ( createVariable( file, "row", "-",  NC_INT, 0, 1, &dimensionId )
             != -1 ) {

          if ( createVariable( file, "count", "-",  NC_INT, 0, 1, &dimensionId )
               != -1 ) {

            if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {
              const Integer index = dataVariableIndex( data );
              Name variable = "";
              aggregateName(data->variable[index], hoursPerTimestep, variable);

              if ( createVariable( file, variable, data->units[ index ],
                                   NC_FLOAT, 1, 1, &dimensionId ) != -1 ) {
                UTCTimestamp timestamp;
                Line history = "";
                appendToLine( history, data->note );
                appendToLine( history, ",XDRConvert" );
                toUTCTimestamp( data->timestamps[ 0 ], timestamp );

                result = writeStandardContents( file, history, timestamp,
                                                dimensionId, 0, 0 );
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
INPUTS:  Integer file                  NetCDF file to write to.
         const Data* data              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSData( Integer file, const Data* data,
                                         const Parameters* parameters ) {

  PRE05( file != -1, isValidData( data ), data->variables < 11,
         data->totalRegriddedPoints > 0, isValidParameters( parameters ) );

  const Integer count = data->totalRegriddedPoints;
  const Integer index = dataVariableIndex( data );
  Integer result =
    writeAllIntData( file, "column", count, 1, 1, 1, data->columns );

  if ( result ) {
    result = writeAllIntData( file, "row", count, 1, 1, 1, data->rows );

    if ( result ) {
      result = writeAllIntData( file, "count", count, 1, 1, 1, data->counts );

      if ( result ) {
        result =
          writeAllData( file, "longitude", count, 1,1,1, data->gridLongitudes);

        if ( result ) {
          result =
             writeAllData( file, "latitude", count, 1,1,1, data->gridLatitudes);

          if ( result ) {
            const Integer hoursPerTimestep =
              parameters->aggregationTimesteps ?
              parameters->aggregationTimesteps
              : 1;
            Name variable = "";
            aggregateName( data->variable[ index ], hoursPerTimestep, variable);
            result =
              writeAllData( file, variable, count, 1,1,1, data->gridData );

            if ( result ) {
              const Integer hoursPerTimestep =
                parameters->aggregationTimesteps ?
                  parameters->aggregationTimesteps
                : 1;

              timeData( data->timesteps, hoursPerTimestep, count,
                        data->outputPoints, data->gridData );

              result =
                writeAllData( file, "time", count, 1, 1, 1, data->gridData );
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
PURPOSE: writeBufferedRegriddedCOARDSData - Write COARDS-format data to file
         by buffering each timestep of data from the temp regrid file.
INPUTS:  Integer file                  NetCDF file to write to.
         const Data* data              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeBufferedRegriddedCOARDSData( Integer file,
                                                 Data* data,
                                                const Parameters* parameters) {

  PRE05( file != -1, isValidData( data ), IN3( data->variables, 11, 12 ),
         data->totalRegriddedPoints > 0, isValidParameters( parameters ) );

  Stream* tempFile = newFileStream( parameters->regridFileName, "rb" );
  Integer result = tempFile != 0;

  if ( result ) {
    const Integer timesteps = data->timesteps;
    Integer timestep = 0;
    Integer offset = 0;
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
    Name variable = "";
    const Integer index = dataVariableIndex( data );
    aggregateName( data->variable[ index ], hoursPerTimestep, variable );

    for ( timestep = 0; AND2( result, timestep < timesteps ); ++timestep ) {
      const Integer count = data->outputPoints[ timestep ];

      if ( count ) {
        result = readRegriddedXDRDataFromTempFile( tempFile, timestep, data );

        if ( result ) {
          result =
            writeSomeIntegerData( file, "column", offset, count, 1, 1, 1,
                                  data->columns );
          if ( result ) {
            result =
              writeSomeIntegerData( file, "row", offset, count, 1, 1, 1,
                                    data->rows );

            if ( result ) {
              result =
                writeSomeIntegerData( file, "count", offset, count, 1, 1, 1,
                                      data->counts );

              if ( result ) {
                result =
                  writeSomeData( file, "longitude", offset, count, 1, 1, 1,
                                 data->gridLongitudes );



                if ( result ) {
                  result =
                    writeSomeData( file, "latitude", offset, count, 1, 1, 1,
                                   data->gridLatitudes );

                  if ( result ) {
                    result =
                      writeSomeData( file, variable, offset, count, 1, 1, 1,
                                     data->gridData );
                  }
                }
              }
            }
          }
        }

        offset += count;
      }
    }

    if ( result ) {
      const Integer count = data->totalRegriddedPoints;
      Real* allTimeData = NEW_ZERO( Real, count );
      result = allTimeData != 0;

      if ( result ) {
        timeData( data->timesteps, hoursPerTimestep, count,
                  data->outputPoints, allTimeData );
        result = writeAllData( file, "time", count, 1, 1, 1, allTimeData );
      }

      FREE( allTimeData );
    }

    FREE_OBJECT( tempFile );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPI - Write regridded IOAPI-format data.
INPUTS:  Data* data                    Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPI( Data* data,
                                    const Parameters* parameters ) {

  PRE03( isValidData( data ), isValidParameters( parameters ),
         parameters->grid );

  Integer result = 0;
  const Integer layers = 1;
  const Integer rows = parameters->grid->rows( parameters->grid );
  const Integer columns = parameters->grid->columns( parameters->grid );
  const Integer fileSizeEstimate =
    data->timesteps * layers * rows * columns * 4 * 4 + 10000;
    /* lon,lat,count,var */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedIOAPIHeader( file, hoursPerTimestep,
                                    data, parameters->grid ) ) {
      result =
        writeRegriddedIOAPIData( file, hoursPerTimestep, data,
                                 parameters );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIHeader - Write header to file.
INPUTS:  Integer file              NetCDF file to write to.
         Integer hoursPerTimestep  Hours per timestep: 1, 24, etc.
         const Data* data          Structure to write.
         const Grid* grid          Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Data* data,
                                          const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0,
         isValidData( data ), grid, grid->invariant( grid ) );

  Integer result = 0;
  enum { VARIABLES = 4 }; /* LONGITUDE, LATITUDE, COUNTS, data. */
  const Integer index = dataVariableIndex( data );
  Name variableNames[ VARIABLES ] = { "LONGITUDE", "LATITUDE", "COUNT", "data"};
  Name variableUnits[ VARIABLES ] = { "deg", "deg", "-", "-" };
  const Integer firstTimestamp = fromUTCTimestamp( data->timestamp );
  Line history = "";
  appendToLine( history, data->note );
  appendToLine( history, ",XDRConvert" );
  aggregateName( data->variable[ index ], hoursPerTimestep,
                 variableNames[ VARIABLES - 1 ] );
  variableNames[ VARIABLES - 1 ][ 15 ] = '\0';
  strncpy( variableUnits[ VARIABLES - 1 ], data->units[ index ], 16 );
  variableUnits[ VARIABLES - 1 ][ 16 ] = '\0';
  uppercase( variableNames[ VARIABLES - 1 ] );
  lowercase( variableUnits[ VARIABLES - 1 ] );

  result = writeM3IOHeader( file, data->timesteps, hoursPerTimestep,
                            firstTimestamp, VARIABLES, 1,
                            (const Name*) variableNames,
                            (const Name*) variableUnits,
                            history, grid );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIData - Write IOAPI-format data to file.
INPUTS:  Integer file                  NetCDF file to write to.
         Integer hoursPerTimestep      E.g., 1 or 24.
         const Data* data              Structure to write.
         const Parameters* parameters  Parameters info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        Data* data,
                                        const Parameters* parameters ) {

  PRE04( file != -1, hoursPerTimestep > 0, isValidData( data ),
         isValidParameters( parameters ) );

  Integer result = 0;
  const Grid* const grid = parameters->grid;
  const Integer rows     = grid->rows( grid );
  const Integer columns  = grid->columns( grid );
  const Integer cells    = rows * columns;
  Real* expandedGridData = NEW_ZERO( Real, cells );

  if ( expandedGridData ) {
    const Integer hasCorners = IN3( data->variables, 11, 12 );
    Stream* tempFile =
      hasCorners ? newFileStream( parameters->regridFileName, "rb" ) : 0;
    result = IMPLIES( hasCorners, tempFile );

    if ( result ) {
      const Integer timesteps = data->timesteps;
      const Integer layers = 1;
      const Integer index = dataVariableIndex( data );
      Integer timestep = 0;
      Integer offset = 0;
      Name variable = "";
      aggregateName( data->variable[ index ], hoursPerTimestep, variable );
      variable[ 15 ] = '\0';
      uppercase( variable );
      result = writeM3IOGrid( grid, timesteps, layers, file );

      for ( timestep = 0; AND2( result, timestep < timesteps ); ++timestep ) {
        const Integer points = data->outputPoints[ timestep ];
        const Integer offset2 = hasCorners ? 0 : offset;

        if ( hasCorners ) {

          if ( points > 0 ) {
            result = readRegriddedXDRDataFromTempFile(tempFile, timestep, data);
          } else {
            Integer cell = 0;

            for ( cell = 0; cell < cells; ++cell ) {
              data->counts[ cell ] = IMISS3;
              data->gridData[ cell ] = BADVAL3;
            }
          }
        }

        if ( result ) {
          Integer* const iExpandedGridData = (Integer*) expandedGridData;
          copyIntDataToGrid( points,
                             data->rows    + offset2,
                             data->columns + offset2,
                             data->counts  + offset2,
                             1, rows, columns, iExpandedGridData );

          result =
            writeM3IOData( file, "COUNT", timestep, 1, rows, columns,
                           iExpandedGridData );

          if ( result ) {
            copyDataToGrid( points,
                            data->rows     + offset2,
                            data->columns  + offset2,
                            data->gridData + offset2,
                            1.0, 1, rows, columns, expandedGridData );

            result =
              writeM3IOData( file, variable, timestep, 1, rows, columns,
                             expandedGridData );
          }
        }

        offset += points;
      }

      FREE_OBJECT( tempFile );
      FREE( expandedGridData );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: regridDataWithCorners - Regrid and mean-aggregate quadrilateral data.
 INPUTS:  Parameters* parameters   Parameters:
                                      aggregationTimesteps, method, grid, etc.
         Data* data                Data to regrid.
OUTPUTS: Data* data                Regridded data attributes:
           data->totalRegriddedPoints, data->outputPoints
           data->longitudes, data->latitudes,
           data->gridLongitudes, data->gridLatitudes,
           data->columns, data->rows, data->gridData.
******************************************************************************/

static void regridDataWithCorners( Parameters* const parameters, Data* data ) {

  PRE05( isValidParameters( parameters ),
         parameters->input->ok( parameters->input ),
         parameters->regridFileName,
         IN3( data->variables, 11, 12 ),
         data->totalRegriddedPoints == 0 );
  assert_static( sizeof (size_t) == sizeof (Integer ) );

  const Integer dataVariable = dataVariableIndex( data );
  const Real defaultMinimumValidValue =
    OR4( ! strcmp( data->units[ dataVariable ], "m" ),
         ! strcmp( data->units[ dataVariable ], "m/s" ),
         ! strcmp( data->units[ dataVariable ], "degC" ),
         ! strcmp( data->units[ dataVariable ], "degrees" ) ) ?
    -900.0 : 0.0;
  const double minimumValidValue =
    AND2( parameters->minimumValidValue > AMISS3,
          ! strcmp( data->units[ dataVariable ], "molecules/cm2" ) ) ?
      parameters->minimumValidValue
    : defaultMinimumValidValue;

  Integer ok = 0;

  /*
   * create temp output file
   * allocate and clear a buffer to hold one hour of input and regrid data
   * For each hourly timestep:
   *   read swath data for the hour
   *   project and/or reorder quad vertices
   *   bin the swath data into grid cells
   *   if timestep is an aggregation output timestep then
   *     compute mean of grid cell values
   *     write regridded data to temp output file
   *     clear cell counts/data
   *   end
   * end
   * close temp output file
   */

  /* Create temp file for regridded data: */

  const char* const tempFileName = parameters->regridFileName;
  Stream* tempFile = newFileStream( parameters->regridFileName, "wb" );

  if ( ! tempFile ) {
    failureMessage( "Failed to create temporary file for regrid data.\n",
                    tempFileName );
  } else {

    /*
     * Allocate a buffer that can hold
     * input data for one hour + output data for one hour.
     * Within the buffer data pointers are ordered:
     * data[ maximumPoints ]
     * longitudes[ maximumPoints ]
     * longitudesSW[ maximumPoints ]
     * longitudesSE[ maximumPoints ]
     * longitudesNW[ maximumPoints ]
     * longitudesNE[ maximumPoints ]
     * latitudes[ maximumPoints ]
     * latitudesSW[ maximumPoints ]
     * latitudesSE[ maximumPoints ]
     * latitudesNW[ maximumPoints ]
     * latitudesNE[ maximumPoints ]
     * vx[ maximumPoints * 4 ]
     * vx[ maximumPoints * 4 ]
     * gridLongitudes[ cellCount ]
     * gridLatitudes[ cellCount ]
     * gridColumns[ cellCount ]
     * gridRows[ cellCount ]
     * gridCounts[ cellCount ]
     * gridWeights[ cellCount ]
     * gridData[ cellCount ]
     * Also, separately allocate:
     * outputPoints[ timestep ]
     */

    const Integer variables = data->variables;
    Stream* input = parameters->input;
    const Integer aggregationTimesteps = /* If not given, default to hourly. */
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
    const int weighted = parameters->regrid == AGGREGATE_WEIGHTED;
    Grid* const grid = parameters->grid;
    const Integer timesteps  = data->timesteps;
    const Integer inputSize  = data->maximumPoints;
    const Integer rows = grid->rows( grid );
    const Integer columns = grid->columns( grid );
    const Integer cellCount = rows * columns;
    const Integer outputSize = cellCount;
    const Integer inputVariables  = variables; /* data,lon,lat,2*corners */
    const Integer vertexCoordinates  = 8; /* vxSW,vxSE,vxNE,vxNW, vySW..vyNW*/
    const Integer outputVariables = OUTPUT_REGRID_VARIABLES + 1;
    /* glon,glat,col,row,counts,weights,gdata */
    const Integer dataSize =
      inputVariables    * inputSize +
      vertexCoordinates * inputSize +
      outputVariables   * outputSize;

    DEBUG( fprintf( stderr, "NEW_ZERO( Real, dataSize = %lld ), "
                    "inputSize = %lld, outputSize = %lld, "
                    "inputVariables = %lld, outputVariables = %lld, "
                    "timesteps = %lld, aggregationTimesteps = %lld\n"
                    "totalPoints = %lld, cellCount = %lld\n",
                    dataSize, inputSize, outputSize,
                    inputVariables, outputVariables,
                    timesteps, aggregationTimesteps,
                    data->totalPoints, cellCount ); )

    FREE( data->data );
    FREE( data->outputPoints );
    data->outputPoints = NEW_ZERO( Integer, timesteps );
    data->data = data->outputPoints ? NEW_ZERO( Real, dataSize ) : 0;

    if ( data->data ) {
      const double gridXMinimum = grid->westEdge( grid );
      const double gridYMinimum = grid->southEdge( grid );
      const double cellWidth    = grid->cellWidth( grid );
      const double cellHeight   = grid->cellHeight( grid );
      Integer binnedSomePoints = 0;
      Integer totalRegriddedPoints = 0;
      Integer timestep = 0;
      Integer outputTimestep = 0;
      Integer yyyydddhh00 = ( fromUTCTimestamp( data->timestamp ) / 100) * 100;
      double* longitudesSW = 0;
      double* longitudesSE = 0;
      double* longitudesNW = 0;
      double* longitudesNE = 0;
      double* latitudesSW  = 0;
      double* latitudesSE  = 0;
      double* latitudesNW  = 0;
      double* latitudesNE  = 0;
      double* vx           = 0; /* Projected counter-clockwise x-vertices. */
      double* vy           = 0; /* Projected counter-clockwise y-vertices. */
      Projector* const projector = grid->projector( grid );

      data->longitudes     = data->data                     + inputSize;
      longitudesSW         = data->longitudes               + inputSize;
      longitudesSE         = longitudesSW                   + inputSize;
      longitudesNW         = longitudesSE                   + inputSize;
      longitudesNE         = longitudesNW                   + inputSize;
      data->latitudes      = longitudesNE                   + inputSize;
      latitudesSW          = data->latitudes                + inputSize;
      latitudesSE          = latitudesSW                    + inputSize;
      latitudesNW          = latitudesSE                    + inputSize;
      latitudesNE          = latitudesNW                    + inputSize;
      vx                   = latitudesNE                    + inputSize;
      vy                   = vx                             + inputSize * 4;
      data->gridLongitudes = vy                             + inputSize * 4;
      data->gridLatitudes  = data->gridLongitudes           + outputSize;
      data->columns        = (Integer*) data->gridLatitudes + outputSize;
      data->rows           = data->columns                  + outputSize;
      data->counts         = data->rows                     + outputSize;
      data->weights        = (Real*) data->counts           + outputSize;
      data->gridData       = data->weights                  + outputSize;

      /* For each hourly timestep: */

      do {

        /* Read swath data for the hour: */

        Integer inputPoints = 0;
        Integer outputPoints = 0;
        DEBUG( fprintf( stderr, "readScanDataForTimestamp yyyydddhh00 = %lld\n",
                        yyyydddhh00 ); )

        ok = readScanDataForTimestamp( yyyydddhh00, input, data,
                                       longitudesSW, longitudesSE,
                                       longitudesNW, longitudesNE,
                                       latitudesSW,  latitudesSE,
                                       latitudesNW,  latitudesNE,
                                       &inputPoints );

        if ( ok ) {

          if ( inputPoints > 0 ) {
            Integer binnedPoints = 0;

            CHECK( totalRegriddedPoints + inputPoints <= data->totalPoints );

            DEBUG( fprintf( stderr, "read swath data:"
                            "inputPoints = %lld, \n"
                            "longitudesSW = [%f ... %f]\n"
                            "longitudesSE = [%f ... %f]\n"
                            "longitudesNW = [%f ... %f]\n"
                            "longitudesNE = [%f ... %f]\n"
                            "latitudesSW  = [%f ... %f]\n"
                            "latitudesSE  = [%f ... %f]\n"
                            "latitudesNW  = [%f ... %f]\n"
                            "latitudesNE  = [%f ... %f]\n"
                            "data         = [%f ... %f]\n",
                            inputPoints,
                            longitudesSW[ 0 ], longitudesSW[ inputPoints - 1 ],
                            longitudesSE[ 0 ], longitudesSE[ inputPoints - 1 ],
                            longitudesNW[ 0 ], longitudesNW[ inputPoints - 1 ],
                            longitudesNE[ 0 ], longitudesNE[ inputPoints - 1 ],
                            latitudesSW[ 0 ], latitudesSW[ inputPoints - 1 ],
                            latitudesSE[ 0 ], latitudesSE[ inputPoints - 1 ],
                            latitudesNW[ 0 ], latitudesNW[ inputPoints - 1 ],
                            latitudesNE[ 0 ], latitudesNE[ inputPoints - 1 ],
                            data->data[ 0 ], data->data[ inputPoints - 1 ] ); )

            /* Project and/or reorder quad lonlat vertices into vx, vy: */

            projectAndOrReorderQuadrilateralVertices( (const size_t) inputPoints,
                                                      longitudesSW,
                                                      longitudesSE,
                                                      longitudesNW,
                                                      longitudesNE,
                                                      latitudesSW,
                                                      latitudesSE,
                                                      latitudesNW,
                                                      latitudesNE,
                                                      projector,
                                                      (ProjectFunction)
                                                        (projector ?
                                                         projector->project : 0),
                                                      vx, vy );

            DEBUG( fprintf( stderr, "projected/reordered quad vertices:"
                            "vx = [%f ... %f], "
                            "vy = [%f ... %f]\n",
                            vx[ 0 ], vx[ inputPoints - 1 ],
                            vy[ 0 ], vy[ inputPoints - 1 ] ); )

            /* Bin the swath data into grid cells: */

            binnedPoints =
              binQuadrilateralData( minimumValidValue,
                                    inputPoints,
                                    data->data,
                                    vx, vy,
                                    rows, columns,
                                    gridXMinimum, gridYMinimum,
                                    cellWidth, cellHeight,
                                    (size_t*) data->counts,
                                    weighted ? data->weights : 0,
                                    data->gridData );

            DEBUG( fprintf( stderr, "binned quad data: "
                            "weighted = %d, binnedPoints = %lld, "
                            "data->counts = [%lld ... %lld], "
                            "data->weights = [%f ... %f], "
                            "data->gridData = [%f ... %f]\n",
                            weighted, binnedPoints,
                            data->counts[0], data->counts[ cellCount - 1 ],
                            data->weights[0], data->weights[ cellCount - 1 ],
                            data->gridData[0], data->gridData[cellCount - 1]);)

            if ( binnedPoints ) {
              binnedSomePoints = 1;
            }
          } /* End if inputPoints > 0 */

          /*
           * if timestep is an aggregation timestep then
           *   compute mean of grid cell values
           *   write regridded data to temp output file
           *   clear cell counts/data
           * end
           */

          if ( ( timestep + 1 ) % aggregationTimesteps == 0 ) {
            outputPoints = ! binnedSomePoints ? 0 :
              computeCellMeans( minimumValidValue, cellCount,
                                (size_t*) data->counts,
                                weighted ? data->weights : 0,
                                data->gridData );

            DEBUG(fprintf(stderr, "mean outputPoints = %lld\n", outputPoints);)

            if ( outputPoints ) {

              /*
               * Compute compact arrays of
               * data->gridLongitudes, data->gridLatitudes,
               * data->columns, data->rows and
               * data->counts, data->gridData
               * so afterwards, all these arrays have length outputPoints.
               */

              compactCells( projector,
                           (UnprojectFunction)
                             ( projector ? projector->unproject : 0 ),
                           columns,
                           rows,
                           gridXMinimum,
                           gridYMinimum,
                           cellWidth,
                           cellHeight,
                           outputPoints,
                           (size_t*) data->counts,
                           data->gridData,
                           data->gridLongitudes,
                           data->gridLatitudes,
                           (size_t*) data->columns,
                           (size_t*) data->rows );

            DEBUG( if ( outputPoints > 0 )
                fprintf( stderr,
                        "compactCells:\n"
                        "  columns        = [%lld ... %lld]\n"
                        "  rows           = [%lld ... %lld]\n"
                        "  gridLongitudes = [%f ... %f]\n"
                        "  gridLatitudes  = [%f ... %f]\n"
                        "  gridData       = [%f ... %f]\n"
                        "  counts         = [%lld ... %lld]\n",
                        data->columns[ 0 ], data->columns[ outputPoints - 1 ],
                        data->rows[    0 ], data->rows[    outputPoints - 1 ],
                        data->gridLongitudes[ 0 ],
                        data->gridLongitudes[ outputPoints - 1 ],
                        data->gridLatitudes[ 0 ],
                        data->gridLatitudes[ outputPoints - 1 ],
                        data->gridData[ 0 ],
                        data->gridData[ outputPoints - 1 ],
                        data->counts[ 0 ],
                        data->counts[ outputPoints - 1 ] ); )

              /* Write regridded output data to temp file: */

              appendRegriddedData( tempFile,
                                   outputPoints,
                                   (size_t*) data->counts,
                                   data->gridData,
                                   data->gridLongitudes,
                                   data->gridLatitudes,
                                   (size_t*) data->columns,
                                   (size_t*) data->rows );
              ok = tempFile->ok( tempFile );

              DEBUG( fprintf( stderr, "appended data to temp regrid file, "
                              "ok = %lld\n", ok ); )

              /* Clear data: */

              binnedSomePoints = 0;
              memset( data->data, 0, dataSize * sizeof data->data[0] );

            } /* End if outputPoints > 0 */

            /* Record regridded point count for this output timestep: */

            data->outputPoints[ outputTimestep ] = outputPoints;
            ++outputTimestep;
            totalRegriddedPoints += outputPoints;

          } /* If aggregation output timestep. */
        } /* If read data did not fail. */

        incrementTimestamp( &yyyydddhh00 );
        ++timestep;
      } while ( AND2( timestep < timesteps, ok ) );

      data->totalRegriddedPoints = totalRegriddedPoints;
      data->timesteps = outputTimestep;
    } /* If allocated data. */

    tempFile->flush( tempFile );

    if ( ! ok ) {
      data->totalRegriddedPoints = 0;
    }

    FREE_OBJECT( tempFile );
  }

  DEBUG( fprintf( stderr, "data->totalRegriddedPoints = %lld\n",
                  data->totalRegriddedPoints ); )

  POST02( data->totalRegriddedPoints >= 0,
          IMPLIES( data->totalRegriddedPoints > 0,
                   AND2( minimumItemI( data->outputPoints, data->timesteps )
                           >= 0,
                         maximumItemI( data->outputPoints, data->timesteps )
                           <= data->totalRegriddedPoints ) ) );
}



/******************************************************************************
PURPOSE: appendRegriddedData - Append regridded data to a temp file.
INPUTS:  Stream* output                        Stream to append data to.
         const size_t count                    Number of values to write.
         const size_t cellCounts[ count ]      Cell counts to write.
         const double cellData[ count ]        Cell data to write.
         const double cellLongitudes[ count ]  Cell center longitudes to write.
         const double cellLatitudes[  count ]  Cell center latitudes  to write.
         const size_t cellColumns[    count ]  Cell columns (1-based) to write.
         const size_t cellRows[       count ]  Cell rows    (1-based) to write.
OUTPUTS: Stream* output                        Binary data appended to file.
******************************************************************************/

static void appendRegriddedData( Stream* const output,
                                 const size_t count,
                                 const size_t cellCounts[],
                                 const double cellData[],
                                 const double cellLongitudes[],
                                 const double cellLatitudes[],
                                 const size_t cellColumns[],
                                 const size_t cellRows[] ) {

  PRE014( output,
          output->ok( output ),
          count,
          cellCounts,
          cellData,
          cellLongitudes,
          cellLatitudes,
          cellColumns,
          cellRows,
          validLongitudesAndLatitudes( (Integer) count,
                                       cellLongitudes, cellLatitudes ),
          minimumItemI( (const Integer*) cellCounts,
                        (Integer) count ) > 0,
          minimumItemI( (const Integer*) cellColumns,
                        (Integer) count ) > 0,
          minimumItemI( (const Integer*) cellRows,
                        (Integer) count ) > 0,
          isNanFree( cellData, (Integer) count ) );

  /*
   * Write a block of data (with timesteps = 1, points = count):
   * # IEEE-754 64-bit reals longitudes[timesteps][points]
   * # IEEE-754 64-bit reals latitudes[timesteps][points]
   * # MSB 64-bit integers columns[timesteps][points]
   * # MSB 64-bit integers rows[timesteps][points]
   * # MSB 64-bit integers counts[timesteps][points]
   * # IEEE-754 64-bit reals data[timesteps][points]
   */

  output->write64BitReals( output, cellLongitudes, count );

  if ( output->ok( output ) ) {
    output->write64BitReals( output, cellLatitudes, count );

    if ( output->ok( output ) ) {
      output->write64BitIntegers( output, (const Integer*) cellColumns, count);

      if ( output->ok( output ) ) {
        output->write64BitIntegers( output, (const Integer*) cellRows, count );

        if ( output->ok( output ) ) {
          output->write64BitIntegers( output, (const Integer*) cellCounts,
                                      count );

          if ( output->ok( output ) ) {
            output->write64BitReals( output, cellData, count );
          }
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: regridData - Regrid data.
INPUTS:  Stream* input   Stream to read Data data from.
         Integer method  Regrid method, AGGREGATE_MEAN, AGGREGATE_WEIGHTED.
         Grid* grid      Grid to project and aggregate points into.
         Data* data      Data to regrid.
OUTPUTS: Data* data      Regridded data attributes:
           data->totalRegriddedPoints, data->outputPoints
           data->longitudes, data->latitudes,
           data->gridLongitudes, data->gridLatitudes,
           data->columns, data->rows, data->gridData.
******************************************************************************/

static void regridData( Stream* input, Integer method, Grid* grid, Data* data) {

  PRE08( input, input->isReadable( input ),
         IS_VALID_AGGREGATE_METHOD( method ),
         grid,
         grid->invariant( grid ),
         isValidData( data ),
         data->totalRegriddedPoints == 0,
         data->longitudes == 0 );

  const Integer variables = data->variables;

  if ( IN5( variables, 3, 4, 11, 12 ) ) {
          const Integer dataVariable = dataVariableIndex( data );
    const Real minimumValidValue =
      OR4( ! strcmp( data->units[ dataVariable ], "m" ),
           ! strcmp( data->units[ dataVariable ], "m/s" ),
           ! strcmp( data->units[ dataVariable ], "degC" ),
           ! strcmp( data->units[ dataVariable ], "degrees" ) ) ?
      -900.0 : 0.0;

    /*
     * Allocate the smallest buffer that can hold:
     *   input  data for one hour +
     *   output data for all hours.
     */

    const Integer hasCorners = IN3( variables, 11, 12 );
    const Integer timesteps  = data->timesteps;
    const Integer inputSize  = data->maximumPoints;
    const Integer totalGridCells =
      timesteps * grid->rows( grid ) * grid->columns( grid );
    const Integer outputSize =
      ! hasCorners ? MIN( data->totalPoints, totalGridCells )
      : totalGridCells;
    const Integer inputVariables  = variables;
    const Integer outputVariables = OUTPUT_REGRID_VARIABLES;
    /* glon,glat,col,row,count,gdata */
    const Integer dataSize =
      inputSize * inputVariables + outputSize * outputVariables;

    DEBUG( fprintf( stderr, "NEW_ZERO( Real, dataSize = %lld ), "
                    "inputSize = %lld, outputSize = %lld, "
                    "inputVariables = %lld, outputVariables = %lld, "
                    "timesteps = %lld, hasCorners = %lld, "
                    "totalPoints = %lld, totalGridCells = %lld\n",
                    dataSize, inputSize, outputSize,
                    inputVariables, outputVariables,
                    timesteps, hasCorners, data->totalPoints, totalGridCells ); )

    FREE( data->data );
    FREE( data->outputPoints );
    data->outputPoints = NEW_ZERO( Integer, timesteps );
    data->data = data->outputPoints ? NEW_ZERO( Real, dataSize ) : 0;

    if ( data->data ) {
      Integer totalRegriddedPoints = 0;
      Integer timestep = 0;
      Integer yyyydddhh00   =
        ( fromUTCTimestamp( data->timestamp ) / 100 ) * 100;
      const Integer inputSize2 =
        variables >= 11 ? inputSize * OUTPUT_REGRID_VARIABLES : inputSize;
      double* longitudesSW = 0;
      double* longitudesSE = 0;
      double* longitudesNW = 0;
      double* longitudesNE = 0;
      double* latitudesSW  = 0;
      double* latitudesSE  = 0;
      double* latitudesNW  = 0;
      double* latitudesNE  = 0;

      data->longitudes     = data->data                     + inputSize;
      data->latitudes      = data->longitudes               + inputSize2;
      data->gridLongitudes = data->latitudes                + inputSize2;
      data->gridLatitudes  = data->gridLongitudes           + outputSize;
      data->columns        = (Integer*) data->gridLatitudes + outputSize;
      data->rows           = data->columns                  + outputSize;
      data->counts         = data->rows                     + outputSize;
      data->gridData       = (Real*) data->counts           + outputSize;

      if ( IN3( variables, 11, 12 ) ) {
        longitudesSW = data->longitudes + inputSize;
        longitudesSE = longitudesSW     + inputSize;
        longitudesNW = longitudesSE     + inputSize;
        longitudesNE = longitudesNW     + inputSize;
        latitudesSW  = data->latitudes  + inputSize;
        latitudesSE  = latitudesSW      + inputSize;
        latitudesNW  = latitudesSE      + inputSize;
        latitudesNE  = latitudesNW      + inputSize;
      }

      do {
        Integer inputPoints = 0;
        Integer outputPoints = 0;
        DEBUG( fprintf( stderr, "readScanDataForTimestamp yyyydddhh00 = %lld\n",
                        yyyydddhh00 ); )

        if ( AND2( readScanDataForTimestamp( yyyydddhh00, input, data,
                                             longitudesSW, longitudesSE,
                                             longitudesNW, longitudesNE,
                                             latitudesSW,  latitudesSE,
                                             latitudesNW,  latitudesNE,
                                             &inputPoints ),
                   inputPoints > 0 ) ) {

          CHECK( totalRegriddedPoints + inputPoints <= data->totalPoints );

          if ( IN3( variables, 3, 4 ) ) {
            DEBUG( fprintf( stderr, "grid->regrid( inputPoints = %lld, ... ) "
                           "totalRegriddedPoints = %lld\n",
                            inputPoints, totalRegriddedPoints ); )
            grid->regrid( grid, method, minimumValidValue, inputPoints, 1,
                          data->longitudes, data->latitudes, 0, data->data,
                          0, /* No input vector data. */
                          0, /* No notes. */
                          &outputPoints,
                          data->columns        + totalRegriddedPoints,
                          data->rows           + totalRegriddedPoints,
                          0,
                          data->gridLongitudes + totalRegriddedPoints,
                          data->gridLatitudes  + totalRegriddedPoints,
                          0,
                          data->gridData       + totalRegriddedPoints,
                          0, /* No output vector data. */
                          0  /* No regriddedNotes. */ );
          } else { /* Regrid using swath quadrilateral corners: */
            CHECK( IN3( variables, 11, 12 ) );

            DEBUG( fprintf( stderr, "grid->regridSwath( grid, method = %lld, "
                           "inputPoints = %lld, \n"
                           "longitudesSW = [%f ... %f]\n"
                           "longitudesSE = [%f ... %f]\n"
                           "longitudesNW = [%f ... %f]\n"
                           "longitudesNE = [%f ... %f]\n"
                           "latitudesSW  = [%f ... %f]\n"
                           "latitudesSE  = [%f ... %f]\n"
                           "latitudesNW  = [%f ... %f]\n"
                           "latitudesNE  = [%f ... %f]\n"
                           "data         = [%f ... %f]\n"
                           "... )\n",
                            method,
                            inputPoints,
                           longitudesSW[ 0 ], longitudesSW[ inputPoints - 1 ],
                           longitudesSE[ 0 ], longitudesSE[ inputPoints - 1 ],
                           longitudesNW[ 0 ], longitudesNW[ inputPoints - 1 ],
                           longitudesNE[ 0 ], longitudesNE[ inputPoints - 1 ],
                           latitudesSW[ 0 ], latitudesSW[ inputPoints - 1 ],
                           latitudesSE[ 0 ], latitudesSE[ inputPoints - 1 ],
                           latitudesNW[ 0 ], latitudesNW[ inputPoints - 1 ],
                           latitudesNE[ 0 ], latitudesNE[ inputPoints - 1 ],
                           data->data[ 0 ], data->data[ inputPoints - 1 ] ); )

            grid->regridSwath( grid, method, minimumValidValue, inputPoints,
                               longitudesSW,
                               longitudesSE,
                               longitudesNW,
                               longitudesNE,
                               latitudesSW,
                               latitudesSE,
                               latitudesNW,
                               latitudesNE,
                               data->data,
                               &outputPoints,
                               data->columns        + totalRegriddedPoints,
                               data->rows           + totalRegriddedPoints,
                               data->gridLongitudes + totalRegriddedPoints,
                               data->gridLatitudes  + totalRegriddedPoints,
                               data->gridData       + totalRegriddedPoints );
          }

          DEBUG( fprintf( stderr, "outputPoints = %lld\n", outputPoints ); )
          DEBUG( if ( outputPoints > 0 )
                fprintf( stderr,
                        "columns = [%lld ... %lld]\n"
                        "rows    = [%lld ... %lld]\n"
                        "gridLongitudes = [%f ... %f]\n"
                        "gridLatitudes  = [%f ... %f]\n"
                        "gridData       = [%f ... %f]\n",
                        data->columns[ 0 ], data->columns[ outputPoints - 1 ],
                        data->rows[    0 ], data->rows[    outputPoints - 1 ],
                        data->gridLongitudes[ 0 ],
                        data->gridLongitudes[ outputPoints - 1 ],
                        data->gridLatitudes[ 0 ],
                        data->gridLatitudes[ outputPoints - 1 ],
                        data->gridData[ 0 ],
                        data->gridData[ outputPoints - 1 ] ); )
        }

        data->outputPoints[ timestep ] = outputPoints;
        totalRegriddedPoints += outputPoints;
        incrementTimestamp( &yyyydddhh00 );
        ++timestep;
      } while ( timestep < timesteps );

      data->totalRegriddedPoints = totalRegriddedPoints;
    }
  }

  DEBUG( fprintf( stderr, "data->totalRegriddedPoints = %lld\n",
                  data->totalRegriddedPoints ); )

  POST02( data->totalRegriddedPoints >= 0,
          IMPLIES( data->totalRegriddedPoints > 0,
                   AND8( IN_RANGE( minimumItemI( data->outputPoints,
                                                 data->timesteps ),
                                   0, data->totalRegriddedPoints ),
                         IN_RANGE( maximumItemI( data->outputPoints,
                                                 data->timesteps ),
                                   1, data->totalRegriddedPoints ),
                         IN_RANGE( minimumItemI(data->columns,
                                                data->totalRegriddedPoints),
                                                1, grid->columns( grid ) ),
                         IN_RANGE( maximumItemI(data->columns,
                                                data->totalRegriddedPoints),
                                                1, grid->columns( grid ) ),
                         IN_RANGE( minimumItemI(data->rows,
                                                data->totalRegriddedPoints),
                                                1, grid->rows( grid ) ),
                         IN_RANGE( maximumItemI(data->rows,
                                                data->totalRegriddedPoints),
                                                1, grid->rows( grid ) ),
                         validLongitudesAndLatitudes(
                                                  data->totalRegriddedPoints,
                                                      data->gridLongitudes,
                                                      data->gridLatitudes ),
                         isNanFree( data->gridData,
                                    data->totalRegriddedPoints ) ) ) );
}



/******************************************************************************
PURPOSE: readScanDataForTimestamp - Read all data for a given timestamp for
         regridding.
INPUTS:  Integer  yyyydddhh00  Timestamp to compare.
         Stream*  input        Stream to read Data data from.
         Data*    data         Data structure.
OUTPUTS: Integer* points       Number of points read for this timestamp.
         Data*    data         data->data, longitudes, latitudes.
                               If non-0:
                               longitudesSW, longitudesSE,
                               longitudesNW, longitudesNE,
                               latitudesSW, latitudesSE,
                               latitudesNW, latitudesNE.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readScanDataForTimestamp( Integer yyyydddhh00, Stream* input,
                                         Data* data,
                                         Real longitudesSW[],
                                         Real longitudesSE[],
                                         Real longitudesNW[],
                                         Real longitudesNE[],
                                         Real latitudesSW[],
                                         Real latitudesSE[],
                                         Real latitudesNW[],
                                         Real latitudesNE[],
                                         Integer* points ) {

  PRE011( isValidTimestamp( yyyydddhh00 ),
          input, input->isReadable( input ),
          data, isValidData( data ), points,
          data->longitudes, data->latitudes, data->data,
          IN5( data->variables, 3, 4, 11, 12 ),
          IMPLIES_ELSE( IN3( data->variables, 11, 12 ),
                        NON_ZERO8( longitudesSW, longitudesSE,
                                   longitudesNW, longitudesNE,
                                   latitudesSW,  latitudesSE,
                                   latitudesNW,  latitudesNE ),
                        IS_ZERO8(  longitudesSW, longitudesSE,
                                   longitudesNW, longitudesNE,
                                   latitudesSW,  latitudesSE,
                                   latitudesNW,  latitudesNE ) ) );

  const Integer scans = data->scans;
  const Integer hasCorners = longitudesSW != 0;
  Integer result = 0;
  Integer scan = 0;
  Integer offset = 0;
  *points = 0;

  do {
    const Integer yyyydddhhmm = data->timestamps[ scan ];
    const Integer yyyydddhh002 = ( yyyydddhhmm / 100 ) * 100;

    if ( yyyydddhh002 == yyyydddhh00 ) {
      const Integer count = data->points[ scan ];
      Integer ok = 0;

      if ( hasCorners ) {
        DEBUG( fprintf( stderr, "  %lld reading %lld values of data into offset = %lld\n",
                       yyyydddhh00, count, offset ); )
        ok = readScanData( input, data->variables, count,
                           data->longitudes + offset,
                           data->latitudes  + offset,
                           data->data       + offset,
                           longitudesSW     + offset,
                           longitudesSE     + offset,
                           longitudesNW     + offset,
                           longitudesNE     + offset,
                           latitudesSW      + offset,
                           latitudesSE      + offset,
                           latitudesNW      + offset,
                           latitudesNE      + offset );

      } else {
        ok = readScanData( input, data->variables, count,
                           data->longitudes + offset,
                           data->latitudes  + offset,
                           data->data       + offset,
                           0, 0, 0, 0, 0, 0, 0, 0 );
      }

      if ( ! ok ) {
        scan = scans; /* Stop looping. */
      } else {
        *points += count;
        offset  += count;
      }
    }

    ++scan;
  } while ( scan < scans );

  result = scan == scans;

  POST03( IS_BOOL( result ),
          *points >= 0,
          IMPLIES( AND2( result, *points > 0 ),
                   AND2( validLongitudesAndLatitudes( *points,
                                                      data->longitudes,
                                                      data->latitudes ),
                         isNanFree( data->data, *points ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readScanData - Read all variables for a given scan for regridding.
INPUTS:  Stream* input               Stream to read data from.
         Integer variables           Number of scan variables.
         Integer points              Number of scan points to read.
OUTPUTS: Real longitudes[ points ]  Longitude data.
         Real latitudes[ points ]   Latitude data.
         Real data[ points ]        3rd data variable. Others are skipped.
         Real longitudesSW[ points ]  0 or longitudes of SW vertex.
         Real longitudesSE[ points ]  0 or longitudes of SE vertex.
         Real longitudesNW[ points ]  0 or longitudes of NW vertex.
         Real longitudesNE[ points ]  0 or longitudes of NE vertex.
         Real latitudesSW[  points ]  0 or latitudes  of SW vertex.
         Real latitudesSE[  points ]  0 or latitudes  of SE vertex.
         Real latitudesNW[  points ]  0 or latitudes  of NW vertex.
         Real latitudesNE[  points ]  0 or latitudes  of NE vertex.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readScanData( Stream* input, Integer variables, Integer points,
                             Real longitudes[], Real latitudes[],
                             Real data[],
                             Real longitudesSW[], Real longitudesSE[],
                             Real longitudesNW[], Real longitudesNE[],
                             Real latitudesSW[],  Real latitudesSE[],
                             Real latitudesNW[],  Real latitudesNE[] ) {

  PRE08( input, input->isReadable( input ), IN5( variables, 3, 4, 11, 12 ),
         points > 0,
         longitudes, latitudes, data,
         IMPLIES_ELSE( IN3( variables, 11, 12 ),
                       NON_ZERO8( longitudesSW, longitudesSE,
                                  longitudesNW, longitudesNE,
                                  latitudesSW,  latitudesSE,
                                  latitudesNW,  latitudesNE ),
                       IS_ZERO8(  longitudesSW, longitudesSE,
                                  longitudesNW, longitudesNE,
                                  latitudesSW,  latitudesSE,
                                  latitudesNW,  latitudesNE ) ) );

  Integer result = 0;
  Integer variable = 0;
  enum { CORNER_VARIABLES = 8 }; /* 4 corners (SW,SE,NW,NE) * 2 (lon,lat). */
  const Integer nonCornerVariables =
    variables >= 11 ? variables - CORNER_VARIABLES : variables;

  /* Read 3 or 4 non-corner variables: */

  do {
    Real* const output =
      variable == 0 ? longitudes :
      variable == 1 ? latitudes  :
      data; /* Read Scan_Start_Time into data then overwrite it with AOD. */

    input->read64BitReals( input, output, points );
    DEBUG( fprintf( stderr, "    v = %lld, read [%f ... %f]\n",
                    variable, output[ 0 ], output[ points - 1 ] ); )
    ++variable;
  } while ( AND2( input->ok( input ), variable < nonCornerVariables ) );

  result = input->ok( input );

  if ( AND2( result, variables >= 11 ) ) { /* Read 8 corner variables: */
    Real* cornerCoordinates[ CORNER_VARIABLES ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    cornerCoordinates[ 0 ] = longitudesSW;
    cornerCoordinates[ 1 ] = longitudesSE;
    cornerCoordinates[ 2 ] = longitudesNW;
    cornerCoordinates[ 3 ] = longitudesNE;
    cornerCoordinates[ 4 ] = latitudesSW;
    cornerCoordinates[ 5 ] = latitudesSE;
    cornerCoordinates[ 6 ] = latitudesNW;
    cornerCoordinates[ 7 ] = latitudesNE;
    variable = 0;

    do {
      Real* const output = cornerCoordinates[ variable ];
      input->read64BitReals( input, output, points );
      DEBUG( fprintf( stderr, "    v = %lld, readlonlat [%f ... %f]\n",
                      variable, output[ 0 ], output[ points - 1 ] ); )
      ++variable;
    } while ( AND2( input->ok( input ), variable < CORNER_VARIABLES ) );

    result = input->ok( input );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND2( validLongitudesAndLatitudes( points,
                                                      longitudes, latitudes ),
                         isNanFree( data, points ) ) ) );
  return result;
}





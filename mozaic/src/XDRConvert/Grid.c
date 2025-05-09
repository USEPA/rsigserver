
/******************************************************************************
PURPOSE: Grid.c - Define routines for processing Grid point data.

NOTES:   XDR-format Grid point data looks like this:

Grid 1.0
http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,HRRRSubset
2020-02-17T00:00:00-0000
# Dimensions: timesteps variables rows columns:
24 2 1059 1799
# Variable names:
wind_u wind_v
# Variable units:
m/s m/s
# IEEE-754 64-bit reals longitudes[rows][columns] and
# IEEE-754 64-bit reals latitudes[rows][columns] and
# IEEE-754 64-bit reals data[timesteps][variables][rows][columns]:

Regridded data looks like this:

REGRIDDED-Grid 1.0
http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,XDRConvert
2008-07-03T00:00:00-0000
# timesteps
24
# Variable name:
wind_u wind_v
# Variable units:
m/s m/s
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000.000000 6370000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
459 299 -2556000.000000 -1728000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 64-bit integers points[timesteps] and
# IEEE-754 64-bit reals longitudes[timesteps][points] and
# IEEE-754 64-bit reals latitudes[timesteps][points] and
# MSB 64-bit integers columns[timesteps][points] and
# MSB 64-bit integers rows[timesteps][points] and
# IEEE-754 64-bit reals data[timesteps][points]:

HISTORY: 2020-03-24 plessel.todd@epa.gov
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h>  /* For stderr, fprintf(). */
#endif
#include <string.h> /* For memset(), memcpy(). */

#include <netcdf.h> /* For nc_close(). */

#include <Utilities.h> /* For PRE*(), NEW_ZERO(), Integer, Real, Stream. */

#include <Helpers.h>         /* For Name, timeData(). */
#include <M3IO.h>            /* For writeM3IOHeader(), writeM3IOData(). */
#include <NetCDFUtilities.h> /* For createNetCDFFile(). */
#include <Parameters.h>      /* For Parameters. */

/*================================== TYPES ==================================*/

typedef struct {
  Stream* input;                  /* Input stream to read per timestep. */
  Line     note;                  /* URL to data source. */
  UTCTimestamp startingTimestamp; /* Starting timestamp. */
  Integer  timesteps;             /* Number of hours in time range. */
  Integer  variables;             /* 4+: timestamp, longitude, latitude, pm25*/
          /* Optional: elevation and/or 2d vector variable (wind_u, wind_v). */
  Integer  rows;                  /* Number of rows of data. */
  Integer  columns;               /* Number of columns of data. */
  Name*    variable;              /* variable[ variables ]. */
  Name*    units;                 /* units[ variables ]. */
  Real*    longitudes;            /* longitudes[ rows * columns ]. */
  Real*    latitudes;             /* latitudes[  rows * columns ]. */
  Real*    data;                  /* data[ variables * rows * columns ]. */
  /* Regrid data: */
  Integer  totalRegriddedPoints;  /* Total number of output points. */
  Integer* outputPoints;          /* outputPoints[ timestep ]. */
  Real*    gridLongitudes;        /* gridLongitudes[ totalRegriddedPoints ]. */
  Real*    gridLatitudes;         /* gridLatitudes[ totalRegriddedPoints ]. */
  Integer* gridColumns;           /* gridColumns[ totalRegriddedPoints ]. */
  Integer* gridRows;              /* gridRows[ totalRegriddedPoints ]. */
  Real*    gridData;              /* gridData[ totalRegriddedPoints ]. */
} Data;

typedef Integer (*Writer)( const Data* data, const Parameters* parameters );

typedef struct {
  Integer format;         /* FORMAT_XDR, etc. */
  Writer writer;          /* Routine that writes data in this format. */
  Writer regriddedWriter; /* Routine that writes regridded data in format. */
} Entry;


/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateData( Data* data );

static Integer isValidData( const Data* data );

static Integer isVectorVariable( const Data* data );

static Writer dispatcher( Integer format, Integer regrid );

static Integer readXDR( Stream* input, Data* data );

static Integer readXDRData( Stream* input, Data* data );

static Integer readRegriddedXDR( Stream* input, Data* data );

static Integer readRegriddedXDRData( Stream* input, Data* data );

static Integer readVariablesAndUnits2( Stream* input, Data* data );

static Integer compareRegriddedXDR( const Parameters* parameters, Data* data );

static Integer writeASCII( const Data* data, const Parameters* unused );

static Integer writeCOARDS( const Data* data, const Parameters* parameters );

static Integer writeCOARDSHeader( Integer file, const Data* data );

static Integer writeCOARDSData( Integer file, const Data* data );

static Integer writeRegriddedXDR( const Data* data, const Parameters* unused );

static Integer writeRegriddedASCII(const Data* data, const Parameters* unused);

static Integer writeRegriddedCOARDS( const Data* data,
                                     const Parameters* parameters );

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Data* data);

static Integer writeRegriddedCOARDSData( Integer file, const Data* data,
                                         const Parameters* parameters );

static Integer writeRegriddedIOAPI( const Data* data,
                                    const Parameters* parameters );

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Data* data,
                                          const Grid* grid );

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Data* data,
                                        const Grid* grid );

static void regridData( Integer method, Grid* grid, Data* data );

static void computeBounds( const Integer count, const Real longitudes[],
                           const Real latitudes[], Bounds bounds );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateGrid - Read input and write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void translateGrid( Parameters* parameters ) {

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
      regridData( parameters->regrid, parameters->grid, &data );

      if ( data.totalRegriddedPoints == 0 ) {
        failureMessage( "No points projected onto the grid." );
      } else {

        if ( parameters->aggregationTimesteps ) {
          const Integer isVector2 = isVectorVariable( &data );
          Integer dataVariable = data.variables - 1;
          Integer totalOutputPoints = 0;
          const Integer aggregatedTimesteps =
            aggregateData( parameters->aggregationTimesteps,
                           isVector2,
                           data.timesteps,
                           data.outputPoints,
                           data.gridLongitudes,
                           data.gridLatitudes,
                           0,
                           data.gridColumns,
                           data.gridRows,
                           0,
                           data.gridData,
                           0,
                           &totalOutputPoints );
          data.timesteps = aggregatedTimesteps;
          data.totalRegriddedPoints = totalOutputPoints;

          if ( AND2( parameters->aggregationTimesteps == 24,
                     ! OR2( strstr( data.variable[ dataVariable ], "daily" ),
                            strstr( data.variable[ dataVariable], "DAILY")))) {
            Integer count = 1 + isVector2;

            while ( count-- ) {
              Name dailyName = "";
              memset( dailyName, 0, sizeof dailyName );
              snprintf( dailyName, sizeof dailyName / sizeof *dailyName,
                        "daily_%s",
                        data.variable[ dataVariable ] );
              strncpy( data.variable[ dataVariable ], dailyName,
                       sizeof dailyName / sizeof *dailyName );
              --dataVariable;
            }
          }
        }

        parameters->ok = writer( &data, parameters );
      }
    } else {
      parameters->ok = writer( &data, parameters );
    }
  }

  deallocateData( &data );
  POST0( isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: compareRegriddedGrid - Read regridded-grid input, compare it to
         CMAQ XDR data and write it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedGrid( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok, parameters->input->ok( parameters->input ) );

  if ( ! AND3( ! parameters->regrid, parameters->compareFunction,
               parameters->data ) ) {
    failureMessage( "Invalid input for comparing." );
    parameters->ok = 0;
  } else {
    Data data;
    ZERO_OBJECT( &data );
    parameters->ok = 0;

    DEBUG( fprintf( stderr, "compareRegriddedData()\n" ); )

    if ( readRegriddedXDR( parameters->input, &data ) ) {
      compareFunctionNameUnits( parameters->compareFunction,
                                parameters->convertFunction,
                                data.variable[ 3 ], data.units[ 3 ],
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
PURPOSE: deallocateData - Deallocate contents of point structure.
INPUTS:  Data* data Structure to deallocate contents of.
******************************************************************************/

static void deallocateData( Data* data ) {
  PRE0( data );
  FREE( data->variable );
  FREE( data->longitudes );
  FREE( data->data );
  FREE( data->outputPoints );
  FREE( data->gridLongitudes );
  ZERO_OBJECT( data );
}



/******************************************************************************
PURPOSE: isValidData - Check data structure.
INPUTS:  const Data* data Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidData( const Data* data ) {
  const Integer result =
    AND14( data,
           data->note[ 0 ],
           isValidUTCTimestamp( data->startingTimestamp ),
           data->timesteps > 0,
           IN_RANGE( data->variables, 1, 2 ),
           data->variable,
           data->variable[ 0 ],
           data->variable[ data->variables - 1 ],
           data->units,
           data->units[ 0 ],
           data->units[ data->variables - 1 ],
           IMPLIES( isVectorVariable( data ),
                    AND3( data->variables == 2,
                          ! strcmp( data->units[ data->variables - 2], "m/s" ),
                          ! strcmp( data->units[ data->variables - 1], "m/s"))),
           IMPLIES_ELSE( data->rows > 0,
                         AND2( data->columns > 0, data->data ),
                         IS_ZERO2( data->columns, data->data ) ),
           IMPLIES( data->totalRegriddedPoints > 0,
                    AND6( data->outputPoints,
                          data->gridLongitudes,
                          data->gridLatitudes,
                          data->gridColumns,
                          data->gridRows,
                          data->gridData ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isVectorVariable - Is the data variable a 2d wind vector?
INPUTS:  const Data* data Structure to check.
RETURNS: Integer 1 if vector, else 0.
******************************************************************************/

static Integer isVectorVariable( const Data* data ) {

  PRE05( data, data->variables > 0, data->variable,
         data->variable[ 0 ], data->variable[ data->variables ] );

  const Integer result =
    AND2( data->variables >= 2,
          OR2( AND2( ! strcmp( data->variable[ data->variables - 2 ], "windU" ),
                     ! strcmp( data->variable[ data->variables - 1 ], "windV" ) ),
               AND2( ! strcmp( data->variable[ data->variables - 2 ], "wind_u" ),
                     ! strcmp( data->variable[ data->variables - 1 ], "wind_v" ) ) ) );

  POST0( IS_BOOL( result ) );
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
PURPOSE: readXDR - Read XDR-format input and initialize data.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Data* data     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   File looks like this:
Grid 1.0
http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,HRRRSubset
2020-02-17T00:00:00-0000
# Dimensions: timesteps variables rows columns:
24 2 rows columns
# Variable names:
wind_u wind_v
# Variable units:
m/s m/s
# IEEE-754 64-bit reals longitudes[rows][columns]
# IEEE-754 64-bit reals latitudes[rows][columns]
# IEEE-754 64-bit reals data[timesteps][variables][rows][columns]:
******************************************************************************/

static Integer readXDR( Stream* input, Data* data ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         data, data->rows == 0, data->columns == 0, data->data == 0 );

  data->input = input; /* Must save for buffered reading per timestep later. */
  input->readString( input, data->note, COUNT( data->note ) );
  Integer result = input->ok( input );

  if ( result ) {
    removeTrailingNewline( data->note );
    result = readTimestamp( input, data->startingTimestamp );

    if ( result ) {
      Integer dimensions[ 4 ] = { 0, 0, 0, 0 };
      result = readDimensions( input, COUNT( dimensions ), dimensions );

      if ( result ) {
        data->timesteps = dimensions[ 0 ];
        data->variables = dimensions[ 1 ];
        data->rows      = dimensions[ 2 ];
        data->columns   = dimensions[ 3 ];
        data->variable = NEW_ZERO( Name, data->variables * 2 );
        result = data->variable != 0;

        if ( result ) {
          data->units = data->variable + data->variables;
          result = readVariablesAndUnits( input, data->variables,
                                          data->variable, data->units );

          if ( result ) {
            const char* const lines[] = {
              "# IEEE-754 64-bit reals longitudes[rows][columns] and\n",
              "# IEEE-754 64-bit reals latitudes[rows][columns] and\n",
              "# IEEE-754 64-bit reals data[timesteps][variables][rows][columns]:\n"
            };
            result =
              AND3( readMatchedLine( input, lines[ 0 ] ),
                    readMatchedLine( input, lines[ 1 ] ),
                    readMatchedLine( input, lines[ 2 ] ) );

            if ( result ) {
              result = readXDRData( input, data );
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Grid data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read XDR-format data and initialize data.
INPUTS:  Stream* input  Input stream read from and parse.
         Data* data     data->variables,timesteps,stations.
OUTPUTS: Data* data     data->longitudes, data->latitudes,
                        allocated data->data for 1 timestep of data.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:  Data looks like this:
# IEEE-754 64-bit reals longitudes[rows][columns] and
# IEEE-754 64-bit reals latitudes[rows][columns] and
# IEEE-754 64-bit reals data[timesteps][variables][rows][columns]:
******************************************************************************/

static Integer readXDRData( Stream* input, Data* data ) {

  PRE010( input, input->ok( input ), input->isReadable( input ),
          data, IN3( data->variables, 1, 2 ), data->rows > 0, data->columns > 0,
          data->variable, data->units, data->data == 0 );

  Integer result = 0;
  const Integer points = data->rows * data->columns;
  data->longitudes = NEW_ZERO( Real, 2 * points );

  result = data->longitudes != 0;

  if ( result ) {
    data->latitudes = data->longitudes + points;
    const Integer dataPoints = data->variables * points; /* Only 1 timestep. */
    data->data = NEW_ZERO( Real, dataPoints );

    result = data->data != 0;

    if ( result ) {
      input->read64BitReals( input, data->longitudes, 2 * points );
      result = input->ok( input );

      if ( result ) {
        result =
          validLongitudesAndLatitudes( points,
                                       data->longitudes,
                                       data->latitudes );

        if ( ! result ) {
          failureMessage( "Read invalid longitude-latitude coordinates." );
        } else {
          result = isValidData( data );
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Grid data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDR - Read Regridded XDR-format input & initialize data.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Data* data     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Input data format is:

REGRIDDED-Grid 1.0
http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,XDRConvert
2008-07-03T00:00:00-0000
# timesteps
24
# Variable name:
wind_u wind_v
# Variable units:
m/s m/s
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000.000000 6370000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
459 299 -2556000.000000 -1728000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 64-bit integers points[timesteps] and
# IEEE-754 64-bit reals longitudes[timesteps][points] and
# IEEE-754 64-bit reals latitudes[timesteps][points] and
# MSB 64-bit integers columns[timesteps][points] and
# MSB 64-bit integers rows[timesteps][points] and
# IEEE-754 64-bit reals data[timesteps][points]:

******************************************************************************/

static Integer readRegriddedXDR( Stream* input, Data* data ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         data, data->rows == 0, data->columns == 0, data->data == 0 );

  Integer result = 0;
  data->input = input; /* Must save for buffered reading per timestep later. */
  input->readString( input, data->note, COUNT( data->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( data->note );

    if ( readTimestamp( input, data->startingTimestamp ) ) {

      if ( readDimensions( input, 1, &data->timesteps ) ) {

        if ( readVariablesAndUnits2( input, data ) ) {

          if ( skipInputLines( input, 10 ) ) {
            data->outputPoints = NEW_ZERO( Integer, data->timesteps );

            if ( data->outputPoints ) {
              input->read64BitIntegers( input, data->outputPoints,
                                        data->timesteps );

              if ( input->ok( input ) ) {
                data->totalRegriddedPoints =
                  sum( data->timesteps, data->outputPoints );
                result = readRegriddedXDRData( input, data );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid REGRIDDED-Grid data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRData - Read Regridded XDR-format input & init data.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Data* data     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:  Data looks like this:
# MSB 64-bit integers points[timesteps] and
# IEEE-754 64-bit reals longitudes[timesteps][points] and
# IEEE-754 64-bit reals latitudes[timesteps][points] and
# MSB 64-bit integers columns[timesteps][points] and
# MSB 64-bit integers rows[timesteps][points] and
# IEEE-754 64-bit reals data[timesteps][points]:
******************************************************************************/

static Integer readRegriddedXDRData( Stream* input, Data* data ) {

  PRE09( input, input->ok( input ), input->isReadable( input ),
         data, data->rows == 0, data->columns == 0, data->data == 0,
         data->gridData == 0, data->timesteps > 0 );

  Integer result = 0;
  const Integer count = data->totalRegriddedPoints;

  if ( count > 0 ) {
    const Integer isVector = isVectorVariable( data );
    data->gridLongitudes = NEW_ZERO( Real, count * ( 5 + isVector ) );

    if ( data->outputPoints ) {
      data->gridLatitudes = data->gridLongitudes + count;
      data->gridColumns   = (Integer*) ( data->gridLatitudes + count );
      data->gridRows      = data->gridColumns + count;
      data->gridData      = (Real*) ( data->gridRows + count );

      input->read64BitReals( input, data->gridLongitudes, count * 2 );

      if ( input->ok( input ) ) {
        input->read64BitIntegers( input, data->gridColumns, count * 2 );

        if ( input->ok( input ) ) {
          const size_t count2 = isVector ? count + count : count;
          input->read64BitReals( input, data->gridData, count2 );

          if ( input->ok( input ) ) {
            result = isValidData( data );
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid REGRIDDED-Point data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readVariablesAndUnits2 - Read 1 (e.g., temperature) or 2
         (wind_u, wind_v) sets of variables and units.
INPUTS:  Stream* input Input stream read from and parse.
OUTPUTS: Data* data    Allocated data->variables,variable,units.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
NOTES:   Input data format is:
# Variable name:
wind_u wind_v
# Variable units:
m/s m/s
******************************************************************************/

static Integer readVariablesAndUnits2( Stream* input, Data* data ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         data, data->variables == 0, data->variable == 0 );

  Integer result = 0;
  char line[ 256 ] = "";
  memset( line, 0, sizeof line );
  input->readString( input, line, sizeof line / sizeof *line - 1 );

  if ( OR2( ! strcmp( line, "# Variable name:\n" ),
            ! strcmp( line, "# Variable names:\n" ) ) ) {
    input->readString( input, line, sizeof line / sizeof *line - 1 );
    data->variables = wordsInString( line );

    if ( IN3( data->variables, 1, 2 ) ) {
      data->variable = NEW_ZERO( Name, data->variables * 2 );

      if ( data->variable ) {
        data->units = data->variable + data->variables;

        if ( data->variables == 1 ) {
          result = sscanf( line, "%s\n", data->variable[ 0 ] ) == 1;
        } else {
          result = sscanf( line, "%s %s\n", data->variable[ 0 ],
                           data->variable[ 1 ] ) == 2;
        }

        if ( result ) {
          result = 0;
          input->readString( input, line, sizeof line / sizeof *line - 1 );

          if ( ! strcmp( line, "# Variable units:\n" ) ) {
            input->readString( input, line, sizeof line / sizeof *line - 1 );

            if ( wordsInString( line ) == data->variables ) {

              if ( data->variables == 1 ) {
                result = sscanf( line, "%s\n", data->units[ 0 ] ) == 1;
              } else {
                result = sscanf( line, "%s %s\n", data->units[ 0 ],
                                 data->units[ 1 ] ) == 2;
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    failureMessage( "Invalid Grid header (variables/units)." );
    data->variables = 0;
    FREE( data->variable );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND7( IN3( data->variables, 1, 2 ),
                         data->variable,
                         data->variable[ 0 ][ 0 ],
                         data->variable[ data->variables - 1 ][ 0 ],
                         data->units,
                         data->units[ 0 ][ 0 ],
                         data->units[ data->variables - 1 ][ 0 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: compareRegriddedXDR - Compare Regridded data with CMAQ data.
INPUTS:  const Parameters* parameters  CMAQ data to compare to.
OUTPUTS: Data* data                    Updated data->data.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static Integer compareRegriddedXDR( const Parameters* parameters, Data* data) {

  PRE05( parameters, isValidParameters( parameters ),
         parameters->compareFunction,
         data, isValidData( data ) );

  Integer result = 0;

  if ( ! AND2( ! strcmp( parameters->timestamp, data->startingTimestamp ),
               parameters->timesteps == data->timesteps ) ) {
    failureMessage( "Mismatched time steps (%s %lld)"
                    " for comparison to CMAQ data (%s %lld).",
                    data->startingTimestamp, data->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else {
    Real* const pointData             = data->gridData;
    const Integer* const pointRows    = data->gridRows;
    const Integer* const pointColumns = data->gridColumns;
    const Integer* const pointsPerTimestep = data->outputPoints;
    const Real* const cmaqData = parameters->data;
    CompareFunction comparer   = parameters->compareFunction;
    const Integer timesteps    = parameters->timesteps;
    const Integer firstRow     = parameters->firstRow;
    const Integer lastRow      = parameters->lastRow;
    const Integer firstColumn  = parameters->firstColumn;
    const Integer lastColumn   = parameters->lastColumn;
    const Integer firstLayer   = parameters->firstLayer;
    const Integer lastLayer    = parameters->lastLayer;
    const Integer layers       = lastLayer  - firstLayer  + 1;
    const Integer rows         = lastRow    - firstRow    + 1;
    const Integer columns      = lastColumn - firstColumn + 1;
    const Integer rowsTimesColumns = rows * columns;
    const Integer layersTimesRowsTimesColumns = layers * rowsTimesColumns;
    Integer timestep = 0;
    Integer pointIndex = 0;

    DEBUG( fprintf( stderr, "timesteps = %lld\n", timesteps ); )

    for ( timestep = 0; timestep < timesteps; ++timestep ) {
      const Integer points = pointsPerTimestep[ timestep ];
      const Integer timestepOffset = timestep * layersTimesRowsTimesColumns;
      Integer point = 0;

      DEBUG( fprintf( stderr, "timestep = %lld, points = %lld\n",
                      timestep, points ); )

      for ( point = 0; point < points; ++point, ++pointIndex ) {
        const Integer pointLayer  = 1;
        const Integer pointRow    = pointRows[    pointIndex ];
        const Integer pointColumn = pointColumns[ pointIndex ];
        DEBUG( fprintf( stderr, "  @[ %lld, %lld, %lld ]: ",
                        pointLayer, pointRow, pointColumn ); )

        if ( AND3( IN_RANGE( pointLayer, firstLayer, lastLayer ),
                   IN_RANGE( pointRow, firstRow, lastRow ),
                   IN_RANGE( pointColumn, firstColumn, lastColumn ) ) ) {
          const Integer pointLayer0  = pointLayer  - firstLayer;
          const Integer pointRow0    = pointRow    - firstRow;
          const Integer pointColumn0 = pointColumn - firstColumn;
          const Integer dataIndex =
            timestepOffset + pointLayer0 * rowsTimesColumns +
            pointRow0 * columns + pointColumn0;
          CHECK4( IN_RANGE( pointLayer0,  0, layers  - 1 ),
                  IN_RANGE( pointRow0,    0, rows    - 1 ),
                  IN_RANGE( pointColumn0, 0, columns - 1 ),
                  IN_RANGE( dataIndex, 0,
                            timesteps * layers * rows * columns - 1 ) );
          const Real pointDatum = pointData[ pointIndex ];
          const Real cmaqDatum = cmaqData[ dataIndex ];
          const Real comparedDatum = comparer( pointDatum, cmaqDatum );
          pointData[ pointIndex ] = comparedDatum;
          result = 1;
          DEBUG( fprintf( stderr, "f(%lf, %lf) -> %lf\n",
                          pointDatum, cmaqDatum, comparedDatum ); )
        } else {
          pointData[ pointIndex ] = -9999.0;
          DEBUG( fprintf( stderr, "-9999\n" ); )
        }
      }
    }
  }

  POST02( IS_BOOL( result ), isValidData( data ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format data.
INPUTS:  const Data* data              Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeASCII( const Data* data, const Parameters* unused ) {

  PRE0( isValidData( data ) );

  Integer result = 0;
  Stream* const input = data->input;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer isVector = isVectorVariable( data );
    const char* const headerStart =
      "Timestamp(UTC)\tLongitude(deg)\tLatitude(deg)";

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {

      if ( isVector ) {
        const char* const headerFormat = "\t%s(%s)\t%s(%s)\n";
        output->writeString( output, headerFormat,
                             data->variable[ 0 ], data->units[ 0 ],
                             data->variable[ 1 ], data->units[ 1 ] );
      } else {
        const char* const headerFormat = "\t%s(%s)\n";
        output->writeString( output, headerFormat,
                             data->variable[ 0 ], data->units[ 0 ] );
      }

      if ( output->ok( output ) ) {
        const Integer pointCount = data->rows * data->columns;
        const Real* const longitudes = data->longitudes;
        const Real* const latitudes  = data->latitudes;
        const Real* const values = data->data;
        const Real* const values2 = isVector ? values + pointCount : 0;
        const Integer timesteps = data->timesteps;
        Integer timestep = 0;
        Integer yyyydddhhmm = fromUTCTimestamp( data->startingTimestamp );

        for ( timestep = 0;
              AND3( input->ok( input ), output->ok( output ),
                    timestep < timesteps );
              ++timestep, incrementTimestamp( &yyyydddhhmm ) ) {

          input->read64BitReals( input, data->data,
                                 pointCount * ( 1 + isVector ) );

          if ( input->ok( input ) ) {
            Integer pointIndex = 0;

            for ( pointIndex = 0;
                  AND2( output->ok( output ), pointIndex < pointCount );
                  ++pointIndex ) {
              const Real longitude = longitudes[ pointIndex ];
              const Real latitude  = latitudes[  pointIndex ];
              const Real value     = values[     pointIndex ];
              UTCTimestamp timestamp = "";
              toUTCTimestamp( yyyydddhhmm, timestamp );

              if ( isVector ) {
                const char* const dataFormat =
                  "%s\t%10.4lf\t%10.4lf\t%10.4lf\t%10.4lf\n";
                const Real value2 = values2[ pointIndex ];
                output->writeString( output, dataFormat,
                                     timestamp, longitude, latitude,
                                     value, value2 );
              } else { /* Not vector: */
                const char* const dataFormat =
                  "%s\t%10.4lf\t%10.4lf\t%10.4lf\n";
                output->writeString( output, dataFormat,
                                     timestamp, longitude, latitude, value );
              }
            }
          }
        }
      }
    }

    result = AND2( input->ok( input), output->ok( output ) );
    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDS - Write COARDS-format data.
INPUTS:  const Data* data              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDS( const Data* data,
                            const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    data->timesteps * data->rows * data->columns * ( 2 + data->variables ) * 4
    + 10000;
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file = createNetCDFFile( parameters->netcdfFileName, create64BitFile);

  if ( file != -1 ) {

    if ( writeCOARDSHeader( file, data ) ) {
      result = writeCOARDSData( file, data );
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
  const char* const names[ 3 ] = { "timesteps", "latitude", "longitude" };
  Integer dimensionIds[ 3 ] = { -1, -1, -1 };
  Integer sizes[ 3 ] = { 0, 0, 0 };
  sizes[ 0 ] = data->timesteps;
  sizes[ 1 ] = data->rows;
  sizes[ 2 ] = data->columns;

  if ( createDimensions( file, 3, names, sizes, dimensionIds ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createLongitudeAndLatitude( file, 2, dimensionIds + 1 ) ) {

        if ( createVariable( file, data->variable[ 0 ], data->units[ 0 ],
                             NC_FLOAT, 1, 3, dimensionIds ) != -1 ) {
          const Integer isVector = isVectorVariable( data );

          if ( OR2( ! isVector,
                    createVariable( file, data->variable[ 1 ],
                                    data->units[ 1 ],
                                    NC_FLOAT, 1, 3, dimensionIds ) != -1 ) ) {
            const Integer points = data->rows * data->columns;
            Bounds bounds;
            computeBounds( points, data->longitudes, data->latitudes, bounds );

            if ( writeExtraAttributes( file, (const Real (*)[2]) bounds,
                                       dimensionIds[ 0 ] ) ) {
              Line history = "";
              appendToLine( history, data->note );
              appendToLine( history, ",XDRConvert" );
              result =
                writeStandardContents( file, history,
                                       data->startingTimestamp,
                                       dimensionIds[ 0 ], data->timesteps, 1 );
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
PURPOSE: writeCOARDSData - Write COARDS-format data to file.
INPUTS:  Integer file      NetCDF file to write to.
         const Data* data  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSData( Integer file, const Data* data ) {

  PRE02( file != -1, isValidData( data ) );

  Integer result = 0;
  const Integer isVector = isVectorVariable( data );
  const Integer count = data->rows * data->columns;

  result = writeAllData( file, "longitude", data->rows, data->columns, 1, 1,
                         data->longitudes );

  if ( result ) {
    result = writeAllData( file, "latitude", data->rows, data->columns, 1, 1,
                           data->latitudes );

    if ( result ) {
      Stream* const input = data->input;
      const Integer timesteps = data->timesteps;
      Integer timestep = 0;

      for ( timestep = 0;
            AND3( result, input->ok( input ), timestep < timesteps );
            ++timestep ) {

        input->read64BitReals( input, data->data,
                               data->rows * data->columns * ( 1 + isVector ) );

        if ( input->ok( input ) ) {
          result =
            writeSomeData( file, data->variable[ 0 ], timestep,
                           1, data->rows, data->columns, 1, data->data);

          if ( result ) {

            if ( isVector ) {
              result =
                writeSomeData( file, data->variable[ 1 ], timestep,
                               1, data->rows, data->columns, 1,
                               data->data + count );
            }

            if ( result ) { /* Create and write yyyyddd and hhmmss data: */
              int* buffer = NEW_ZERO( int, timesteps * 3 );
              result = buffer != 0;

              if ( buffer ) {
                int* const yyyyddd = buffer;
                int* const hhmmss = buffer + timesteps;
                float* fhour = (float*) ( hhmmss + timesteps );
                Integer yyyydddhhmm =
                  fromUTCTimestamp( data->startingTimestamp );

                for ( timestep = 0; timestep < timesteps; ++timestep ) {
                  yyyyddd[ timestep ] = yyyydddhhmm / 10000LL;
                  hhmmss[  timestep ] = yyyydddhhmm % 10000LL * 100LL;
                  fhour[   timestep ] = timestep;
                  incrementTimestamp( &yyyydddhhmm );
                }

                result =
                  writeTimeData1( file, timesteps, yyyyddd, hhmmss, fhour );

                FREE( buffer );
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
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  const Data* data              Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Output data format is:

REGRIDDED-Grid 1.0
http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,XDRConvert
2008-07-03T00:00:00-0000
# timesteps
24
# Variable name:
wind_u wind_v
# Variable units:
m/s m/s
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000.000000 6370000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
459 299 -2556000.000000 -1728000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 64-bit integers points[timesteps] and
# IEEE-754 64-bit reals longitudes[timesteps][points] and
# IEEE-754 64-bit reals latitudes[timesteps][points] and
# MSB 64-bit integers columns[timesteps][points] and
# MSB 64-bit integers rows[timesteps][points] and
# IEEE-754 64-bit reals data[timesteps][points]:

******************************************************************************/

static Integer writeRegriddedXDR( const Data* data,
                                  const Parameters* parameters) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = data->timesteps;
    const Integer isVector  = isVectorVariable( data );
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps
      : 1;
    Name variable = "";
    aggregateName( data->variable[ 0 ], hoursPerTimestep, variable );

    if ( isVector ) {
      Name variable2 = "";
      aggregateName( data->variable[ 1 ], hoursPerTimestep, variable2 );
      output->writeString( output,
                           "REGRIDDED-Grid 1.0\n"
                           "%s,XDRConvert\n"
                           "%s\n"
                           "# timesteps\n%"INTEGER_FORMAT"\n"
                           "# Variable name:\n%s %s\n"
                           "# Variable units:\n%s %s\n",
                           data->note,
                           data->startingTimestamp,
                           timesteps,
                           variable,
                           variable2,
                           data->units[ 0 ],
                           data->units[ 1 ] );
    } else {
      output->writeString( output,
                           "REGRIDDED-Grid 1.0\n"
                           "%s,XDRConvert\n"
                           "%s\n"
                           "# timesteps\n%"INTEGER_FORMAT"\n"
                           "# Variable name:\n%s\n"
                           "# Variable units:\n%s\n",
                           data->note,
                           data->startingTimestamp,
                           timesteps,
                           variable,
                           data->units[ 0 ] );
    }

    if ( output->ok( output ) ) {
      writeProjectionAndGrid( parameters->grid, output );

      if ( output->ok( output ) ) {
        output->writeString( output,
          "# MSB 64-bit integers points[timesteps] and\n"
          "# IEEE-754 64-bit reals longitudes[timesteps][points] and\n"
          "# IEEE-754 64-bit reals latitudes[timesteps][points] and\n"
          "# MSB 64-bit integers columns[timesteps][points] and\n"
          "# MSB 64-bit integers rows[timesteps][points] and\n"
          "# IEEE-754 64-bit reals data[timesteps][points]:\n" );

        if ( output->ok( output ) ) {
          output->write64BitIntegers( output, data->outputPoints, timesteps );

          if ( output->ok( output ) ) {
            const Integer points = data->totalRegriddedPoints;
            output->write64BitReals( output, data->gridLongitudes, points );

            if ( output->ok( output ) ) {
              output->write64BitReals( output, data->gridLatitudes, points );

              if ( output->ok( output ) ) {
                output->write64BitIntegers( output, data->gridColumns, points );

                if ( output->ok( output ) ) {
                  output->write64BitIntegers( output, data->gridRows, points );

                  if ( output->ok( output ) ) {
                    const Integer points2 =
                      isVector ? points + points : points;
                    output->write64BitReals( output, data->gridData, points2 );
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
INPUTS:  const Data* data              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedASCII( const Data* data,
                                    const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer isVector = isVectorVariable( data );
    const char* const headerStart =
      "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tCOLUMN(-)\tROW(-)";

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {
      const Integer hoursPerTimestep =
        parameters->aggregationTimesteps ? parameters->aggregationTimesteps
        : 1;
      Name variable = "";
      aggregateName( data->variable[ 0 ], hoursPerTimestep, variable );

      if ( isVector ) {
        const char* const headerFormat = "\t%s(%s)\t%s(%s)\n";
        Name variable2 = "";
        aggregateName( data->variable[ 1 ], hoursPerTimestep, variable2 );
        output->writeString( output, headerFormat,
                             variable,
                             data->units[    0 ],
                             variable2,
                             data->units[    1 ] );
      } else {
        const char* const headerFormat = "\t%s(%s)\n";
        output->writeString( output, headerFormat,
                             variable,
                             data->units[    0 ] );
      }

      if ( output->ok( output ) ) {
        const Integer timesteps = data->timesteps;
        const Real* longitudes  = data->gridLongitudes;
        const Real* latitudes   = data->gridLatitudes;
        const Integer* columns  = data->gridColumns;
        const Integer* rows     = data->gridRows;
        const Real* values      = data->gridData;
        const Real* values2 = isVector ? values + data->totalRegriddedPoints :0;
        Integer timestep = 0;
        Integer yyyydddhhmm = fromUTCTimestamp( data->startingTimestamp );

        /* Write data rows: */

        do {
          const Integer points = data->outputPoints[ timestep ];
          Integer point = 0;
          UTCTimestamp timestamp = "";
          toUTCTimestamp( yyyydddhhmm, timestamp );

          for ( point = 0; point < points; ++point ) {
            const Real longitude = *longitudes++;
            const Real latitude  = *latitudes++;
            const Integer column = *columns++;
            const Integer row    = *rows++;
            const Real value     = *values++;

            if ( isVector ) {
              const Real value2 = *values2++;
              const char* const dataFormat =
                "%s\t%10.4lf\t%10.4lf"
                "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT
                "\t%10.4lf\t%10.4lf\n";
              output->writeString( output, dataFormat,
                                   timestamp, longitude, latitude,
                                   column, row, value, value2 );
            } else {
              const char* const dataFormat =
                "%s\t%10.4lf\t%10.4lf"
                "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%10.4lf\n";
              output->writeString( output, dataFormat,
                                   timestamp, longitude, latitude,
                                   column, row, value );
            }

            if ( ! output->ok( output ) ) {
              point = points;
              timestep = timesteps;
            }
          }

          yyyydddhhmm = offsetTimestamp( yyyydddhhmm, hoursPerTimestep );
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
INPUTS:  const Data* data              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDS( const Data* data,
                                     const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    data->totalRegriddedPoints * 9 * 4 + 10000;
    /* lon, lat, col, row, time, var, var2, + hdr. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedCOARDSHeader( file, hoursPerTimestep, data ) ) {
      result =
        writeRegriddedCOARDSData( file, data,  parameters );
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
         Integer hoursPerTimestep  Hours per timestep: 1, 24, etc.
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

  if ( createDimensions( file, 1, &dimensionName, &dimension, &dimensionId )) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "column", "-",  NC_INT, 0, 1, &dimensionId )
           != -1 ) {

        if ( createVariable( file, "row", "-", NC_INT, 0, 1, &dimensionId )
             != -1 ) {

          if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {
            const Integer isVector = isVectorVariable( data );
            Name variable = "";
            aggregateName( data->variable[ 0 ], hoursPerTimestep, variable );
            result = createVariable( file, variable, data->units[ 0 ], NC_FLOAT,
                                     1, 1, &dimensionId ) != -1;
            if ( result ) {

              if ( isVector ) {
                aggregateName( data->variable[ 1 ], hoursPerTimestep, variable );
                result =
                  createVariable( file, variable, data->units[ 1 ], NC_FLOAT,
                                  1, 1, &dimensionId ) != -1;
              }

              if ( result ) {
                Line history = "";
                appendToLine( history, data->note );
                appendToLine( history, ",XDRConvert" );
                result =
                  writeStandardContents( file, history, data->startingTimestamp,
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
INPUTS:  Integer file      NetCDF file to write to.
         const Data* data  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSData( Integer file, const Data* data,
                                         const Parameters* parameters ) {

  PRE03( file != -1, isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer count = data->totalRegriddedPoints;

  if ( writeAllIntData( file, "column", count, 1, 1, 1, data->gridColumns ) ) {

    if ( writeAllIntData( file, "row", count, 1, 1, 1, data->gridRows ) ) {

      if ( writeAllData( file, "longitude", count, 1, 1, 1,
                           data->gridLongitudes ) ) {

        if ( writeAllData( file, "latitude", count, 1, 1, 1,
                           data->gridLatitudes ) ) {
          const Integer hoursPerTimestep =
            parameters->aggregationTimesteps ? parameters->aggregationTimesteps
            : 1;
          Name variable = "";
          aggregateName( data->variable[ 0 ], hoursPerTimestep, variable );

          if ( writeAllData( file, variable, count, 1, 1, 1, data->gridData)) {
            const Integer isVector = isVectorVariable( data );
            result = 1;

            if ( isVector ) {
              aggregateName( data->variable[ 1 ], hoursPerTimestep, variable );
              result =
                writeAllData( file, variable, count, 1, 1, 1,
                              data->gridData + count );
            }

            if ( result ) {
              timeData( data->timesteps, hoursPerTimestep, count,
                        data->outputPoints, data->gridData );

              result = writeAllData( file, "time", count, 1, 1, 1,
                                     data->gridData );
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
INPUTS:  const Data* data              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPI( const Data* data,
                                    const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    data->totalRegriddedPoints * 5 * 4 + 10000; /* lon, lat, var, var2 + hdr */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedIOAPIHeader( file, hoursPerTimestep,
                                    data, parameters->grid ) ) {
      result = writeRegriddedIOAPIData( file, hoursPerTimestep,
                                        data, parameters->grid );
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

  PRE05( file != -1, hoursPerTimestep,
         isValidData( data ), grid, grid->invariant( grid ) );

  Integer result = 0;
  enum { VARIABLES = 3 }; /* LONGITUDE, LATITUDE, var. */
  Name variableNames[ VARIABLES + 1 ] =
    { "LONGITUDE", "LATITUDE", "var", "WIND_V" };
  Name variableUnits[ VARIABLES + 1 ] = { "deg", "deg", "m/s", "m/s" };
  const Integer firstTimestamp = fromUTCTimestamp( data->startingTimestamp );
  const Integer isVector = isVectorVariable( data );
  Line history = "";
  appendToLine( history, data->note );
  appendToLine( history, ",XDRConvert" );
  aggregateName( data->variable[ 0 ], hoursPerTimestep, variableNames[ 2 ] );
  variableNames[ 2 ][ 15 ] = '\0';
  strncpy( variableUnits[ VARIABLES - 1 ], data->units[ 0 ], 16 );
  variableUnits[ 2 ][ 16 ] = '\0';
  uppercase( variableNames[ 2 ] );
  lowercase( variableUnits[ 2 ] );

  if ( isVector ) {
    aggregateName( data->variable[ 1 ], hoursPerTimestep, variableNames[ 3 ] );
    variableNames[ 3 ][ 15 ] = '\0';
    strncpy( variableUnits[ 3 ], data->units[ 1 ], 16 );
    variableUnits[ 3 ][ 16 ] = '\0';
    uppercase( variableNames[ 3 ] );
    lowercase( variableUnits[ 3 ] );
  }

  result = writeM3IOHeader( file, data->timesteps, hoursPerTimestep,
                            firstTimestamp,
                            VARIABLES + isVector,
                            1,
                            (const Name*) variableNames,
                            (const Name*) variableUnits,
                            history, grid );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIData - Write IOAPI-format data to file.
INPUTS:  Integer file              NetCDF file to write to.
         Integer hoursPerTimestep  Hours per timestep: 1, 24, etc.
         const Data* data          Structure to write.
         const Grid* grid          Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Data* data,
                                        const Grid* grid ) {

  PRE04( file != -1, isValidData( data ), grid, grid->invariant( grid ) );

  Integer result = 0;
  const Integer isVector = isVectorVariable( data );
  const Integer layers  = 1;
  const Integer rows    = grid->rows( grid );
  const Integer columns = grid->columns( grid );
  const Integer cells   = layers * rows * columns;
  Real* gridData = NEW_ZERO( Real, cells  );

  if ( gridData ) {
    const Integer timesteps = data->timesteps;
    Integer offset = 0;
    Integer timestep = 0;

    if ( writeM3IOGrid( grid, timesteps, layers, file ) ) {

      do {
        const Integer points = data->outputPoints[ timestep ];
        Name variable = "";
        aggregateName( data->variable[ 0 ], hoursPerTimestep, variable );
        variable[ 15 ] = '\0';
        uppercase( variable );
        copyDataToGrid( points,
                          data->gridRows     + offset,
                          data->gridColumns  + offset,
                          data->gridData     + offset,
                          1.0, layers, rows, columns, gridData );

        if ( ! writeM3IOData( file, variable,
                              timestep, layers, rows, columns, gridData ) ) {
          timestep = timesteps;
        } else {

          if ( isVector ) {
            const Real* gridData2 =
              data->gridData + data->totalRegriddedPoints;
            memset( variable, 0, sizeof variable );
            strncpy( variable, data->variable[ 1 ], 16 );
            variable[ 16 ] = '\0';
            uppercase( variable );
            copyDataToGrid( points,
                            data->gridRows    + offset,
                            data->gridColumns + offset,
                            gridData2         + offset,
                            1.0, layers, rows, columns, gridData );
          }

          if ( ! writeM3IOData( file, variable, timestep, layers, rows,
                                columns, gridData ) ) {
            timestep = timesteps;
          }
        }

        offset += points;
        ++timestep;
      } while ( timestep < timesteps );

      result = timestep == timesteps;
    }

    FREE( gridData );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: regridData - Regrid data.
INPUTS:  Integer method  E.g., AGGREGATE_MEAN.
         Grid* grid      Grid to project and aggregate points into.
         Data* data      Data to regrid.
OUTPUTS: Data* data      Regridded data.
******************************************************************************/

static void regridData( Integer method, Grid* grid, Data* data ) {

  PRE05( IS_VALID_AGGREGATE_METHOD( method ),
         grid,
         grid->invariant( grid ),
         isValidData( data ),
         data->totalRegriddedPoints == 0 );

  Integer totalRegriddedPoints = 0;
  const Integer isVector = isVectorVariable( data );
  const Integer timesteps = data->timesteps;
  const Integer points = data->rows * data->columns;
  data->outputPoints = NEW_ZERO( Integer, timesteps );

  if ( data->outputPoints ) {
    const Integer outputVariables = 5 + isVector;
    /* lon,lat,col,row,var + var2 */
    const Integer gridPoints = grid->rows( grid ) * grid->columns( grid );
    const Integer outputPoints = MIN( points, gridPoints );
    const Integer variableSize = timesteps * outputPoints;
    const Integer outputSize = variableSize * outputVariables;
    data->gridLongitudes = NEW_ZERO( Real, outputSize );

    if ( data->gridLongitudes ) {
      Stream* const input = data->input;
      Integer timestep = 0;
      Real* gridData2 = 0;
      data->gridLatitudes = data->gridLongitudes + variableSize;
      data->gridColumns = (Integer*) ( data->gridLatitudes + variableSize );
      data->gridRows = data->gridColumns + variableSize;
      data->gridData = (Real*) ( data->gridRows + variableSize );

      if ( isVector ) {
        gridData2 = data->gridData + variableSize;
      }

      for ( timestep = 0;
            AND2( input->ok( input ), timestep < timesteps );
            ++timestep ) {

        input->read64BitReals( input, data->data, points * ( 1 + isVector ) );

        if ( input->ok( input ) ) {
          Integer outputPoints = 0;
          const Real minimumValidValue = -500.0;

          grid->regrid( grid, method, minimumValidValue, points, 1,
                        data->longitudes, data->latitudes,
                        0, /* No elevations. */
                        data->data,
                        isVector ? data->data + points : 0,
                        0, /* No notes. */
                        &outputPoints,
                        data->gridColumns    + totalRegriddedPoints,
                        data->gridRows       + totalRegriddedPoints,
                        0, /* No layers. */
                        data->gridLongitudes + totalRegriddedPoints,
                        data->gridLatitudes  + totalRegriddedPoints,
                        0, /* No gridElevations. */
                        data->gridData       + totalRegriddedPoints,
                        isVector ? gridData2 + totalRegriddedPoints : 0,
                        0 /* No regriddedNotes. */ );

          data->outputPoints[ timestep ] = outputPoints;
          totalRegriddedPoints += outputPoints;
        }
      }
    }
  }

  data->totalRegriddedPoints = totalRegriddedPoints;

  POST02( data->totalRegriddedPoints >= 0,
          IMPLIES( data->totalRegriddedPoints > 0,
                   AND6( IN_RANGE( minimumItemI( data->gridColumns,
                                                 data->totalRegriddedPoints ),
                                                 1, grid->columns( grid ) ),
                         IN_RANGE( maximumItemI( data->gridColumns,
                                                 data->totalRegriddedPoints ),
                                                 1, grid->columns( grid ) ),
                         IN_RANGE( minimumItemI( data->gridRows,
                                                 data->totalRegriddedPoints ),
                                                 1, grid->rows( grid ) ),
                         IN_RANGE( maximumItemI( data->gridRows,
                                                 data->totalRegriddedPoints ),
                                                 1, grid->rows( grid ) ),
                         validLongitudesAndLatitudes(data->totalRegriddedPoints,
                                                     data->gridLongitudes,
                                                     data->gridLatitudes ),
                         isNanFree( data->gridData,
                                    data->totalRegriddedPoints
                                    * ( 1 + isVector ) ) ) ) );
}



/******************************************************************************
PURPOSE: computeBounds - Compute 2d bounds of coordinates.
INPUTS:  const Integer count             Number of points.
         const Real longitudes[ count ]  Longitudes of points.
         const Real latitudes[  count ]  Latitudes of points.
OUTPUTS: Bounds bounds                   Bounds of coordinates.
******************************************************************************/

static void computeBounds( const Integer count, const Real longitudes[],
                           const Real latitudes[], Bounds bounds ) {

  PRE04( count > 0, longitudes, latitudes, bounds );
  Real longitudeMinimum = *longitudes;
  Real longitudeMaximum = longitudeMinimum;
  Real latitudeMinimum = *latitudes;
  Real latitudeMaximum = latitudeMinimum;
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {
    const Real longitude = longitudes[ index ];
    const Real latitude  = latitudes[  index ];

    if ( longitude < longitudeMinimum ) {
      longitudeMinimum = longitude;
    } else if ( longitude > longitudeMaximum ) {
      longitudeMaximum = longitude;
    }

    if ( latitude < latitudeMinimum ) {
      latitudeMinimum = latitude;
    } else if ( latitude > latitudeMaximum ) {
      latitudeMaximum = latitude;
    }
  }

  bounds[ LONGITUDE ][ MINIMUM ] = longitudeMinimum;
  bounds[ LONGITUDE ][ MAXIMUM ] = longitudeMaximum;
  bounds[ LATITUDE  ][ MINIMUM ] = latitudeMinimum;
  bounds[ LATITUDE  ][ MAXIMUM ] = latitudeMaximum;

  POST0( isValidBounds( (const Real(*)[2]) bounds ) );
}




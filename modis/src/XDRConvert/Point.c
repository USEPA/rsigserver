
/******************************************************************************
PURPOSE: Point.c - Define routines for processing Point data.

NOTES:   XDR-format Point data looks like this:

Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/
2008-07-03T00:00:00-0000 2008-07-03T23:59:59-0000
# Dimensions: variables points
3 72
# Variable names:
timestamp longitude latitude pm25
# Variable units:
yyyymmddhhmmss deg deg metric_tons
# IEEE-754 64-bit reals data[variables][points]:

or, with elevations:

Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/
2008-07-03T00:00:00-0000 2008-07-03T23:59:59-0000
# Dimensions: variables points
4 72
# Variable names:
timestamp longitude latitude elevation pm25
# Variable units:
yyyymmddhhmmss deg deg m metric_tons
# IEEE-754 64-bit reals data[variables][points# IEEE-754 64-bit reals data[variables][points]:

or, with elevations and vector2 variable:

Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/
2008-07-03T00:00:00-0000 2008-07-03T23:59:59-0000
# Dimensions: variables points
6 72
# Variable names:
timestamp longitude latitude elevation wind_u wind_v
# Variable units:
yyyymmddhhmmss deg deg m m/s m/s
# IEEE-754 64-bit reals data[variables][points]:

or, with elevations and vector2 variable and notes-per-point:

Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/
2008-07-03T00:00:00-0000 2008-07-03T23:59:59-0000
# Dimensions: variables points
6 72
# Variable names:
timestamp longitude latitude elevation wind_u wind_v
# Variable units:
yyyymmddhhmmss deg deg m m/s m/s
# char notes[points][80] and
# IEEE-754 64-bit reals data[variables][points# IEEE-754 64-bit reals data[variables][points]:

Regridded data looks like this:

REGRIDDED-Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/,XDRConvert
2008-07-03T00:00:00-0000
# timesteps
24
# Variable name:
pm25
# Variable units:
metric_tons
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000.000000 6370000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
459 299 -2556000.000000 -1728000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 32-bit integers points[timesteps] and
# IEEE-754 32-bit reals longitudes[timesteps][points] and
# IEEE-754 32-bit reals latitudes[timesteps][points] and
# MSB 32-bit integers columns[timesteps][points] and
# MSB 32-bit integers rows[timesteps][points] and
# IEEE-754 32-bit reals data[timesteps][points]:

or, with elevation:

REGRIDDED-Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/,XDRConvert
2008-07-03T00:00:00-0000
# timesteps
24
# Variable name:
pm25
# Variable units:
metric_tons
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000.000000 6370000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
459 299 -2556000.000000 -1728000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 32-bit integers points[timesteps] and
# IEEE-754 32-bit reals longitudes[timesteps][points] and
# IEEE-754 32-bit reals latitudes[timesteps][points] and
# IEEE-754 32-bit reals elevations[timesteps][points] and
# MSB 32-bit integers columns[timesteps][points] and
# MSB 32-bit integers rows[timesteps][points] and
# MSB 32-bit integers layers[timesteps][points] and
# IEEE-754 32-bit reals data[timesteps][points]:

HISTORY: 2018-10-27 plessel.todd@epa.gov
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
  Line     note;                  /* URL to data source. */
  UTCTimestamp startingTimestamp; /* Starting timestamp. */
  UTCTimestamp endingTimestamp;   /* Ending timestamp. */
  Integer  timesteps;             /* Number of hours in time range. */
  Integer  variables;             /* 4+: timestamp, longitude, latitude, pm25*/
          /* Optional: elevation and/or 2d vector variable (wind_u, wind_v). */
  Integer  points;                /* Number of data points. */
  Name*    variable;              /* variable[ variables ]. */
  Name*    units;                 /* units[ variables ]. */
  Note*    notes;                 /* 0 or notes[ points ]. */
  Real*    data;                  /* data[ variables * points ]. */
  /* Regrid data: */
  Integer  totalRegriddedPoints;  /* Total number of output points. */
  Integer* outputPoints;          /* outputPoints[ timestep ]. */
  Real*    gridLongitudes;        /* gridLongitudes[ totalRegriddedPoints ]. */
  Real*    gridLatitudes;         /* gridLatitudes[ totalRegriddedPoints ]. */
  Real*    gridElevations;        /* gridElevations[ totalRegriddedPoints ]. */
  Integer* columns;               /* columns[ totalRegriddedPoints ]. */
  Integer* rows;                  /* rows[ totalRegriddedPoints ]. */
  Integer* layers;                /* layers[ totalRegriddedPoints ]. */
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

static Integer writeConvertedTimeData( const Integer file,
                                       const UTCTimestamp startingTimestamp,
                                       const Integer count,
                                       const Real* const timestamps );

static Integer writeRegriddedXDR( const Data* data,
                                  const Parameters* parameters );

static Integer writeRegriddedASCII( const Data* data,
                                    const Parameters* parameters );

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

static Integer pointsMatchingHour( const Integer yyyydddhh,
                                   const Data* const data,
                                   const Integer start );

static void computeBounds( const Integer count, const Real longitudes[],
                           const Real latitudes[], Bounds bounds );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translatePoint - Read input and write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void translatePoint( Parameters* parameters ) {

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
                           data.gridElevations,
                           data.columns,
                           data.rows,
                           data.layers,
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
PURPOSE: compareRegriddedPoint - Read regridded-point input, compare it to
         CMAQ XDR data and write it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedPoint( Parameters* parameters ) {

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
  FREE( data->data );
  FREE( data->notes );
  FREE( data->outputPoints );
  FREE( data->gridLongitudes );
  ZERO_OBJECT( data );
}



/******************************************************************************
PURPOSE: isValidData - Check point structure.
INPUTS:  const Data* data Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidData( const Data* data ) {
  const Integer result =
    AND22( data,
           data->note[ 0 ],
           isValidUTCTimestamp( data->startingTimestamp ),
           data->timesteps > 0,
           IN_RANGE( data->variables, 4, 7 ),
           data->variable,
           data->variable[ 0 ],
           data->variable[ data->variables - 1 ],
           data->units,
           data->units[ 0 ],
           data->units[ data->variables - 1 ],
           OR2( ! strcmp( data->variable[ 0 ], "Timestamp" ),
                ! strcmp( data->variable[ 0 ], "timestamp" ) ),
           OR2( ! strcmp( data->variable[ 1 ], "Longitude" ),
                ! strcmp( data->variable[ 1 ], "longitude" ) ),
           OR2( ! strcmp( data->variable[ 2 ], "Latitude" ),
                ! strcmp( data->variable[ 2 ], "latitude" ) ),
           ! strcmp( data->units[ 0 ], "yyyymmddhhmmss" ),
           ! strcmp( data->units[ 1 ], "deg" ),
           ! strcmp( data->units[ 2 ], "deg" ),
           IMPLIES( OR2( ! strcmp( data->variable[ 3 ], "Elevation" ),
                         ! strcmp( data->variable[ 3 ], "elevation" ) ),
                    AND2( data->variables > 4,
                          ! strcmp( data->units[ 3 ], "m" ) ) ),
           IMPLIES( isVectorVariable( data ),
                    AND3( data->variables > 4,
                          ! strcmp( data->units[ data->variables - 2], "m/s" ),
                          ! strcmp( data->units[ data->variables - 1], "m/s"))),
           IMPLIES_ELSE( data->points > 0,
                         data->data,
                         data->data == 0 ),
           IMPLIES( data->notes,
                    AND2( data->notes[ 0 ], data->notes[ data->points - 1 ] ) ),
           IMPLIES( data->totalRegriddedPoints > 0,
                    AND7( data->outputPoints,
                          data->gridLongitudes,
                          data->gridLatitudes,
                          IMPLIES(OR2(! strcmp(data->variable[3], "Elevation"),
                                      ! strcmp(data->variable[3], "elevation")),
                                   data->gridElevations ),
                          data->columns,
                          data->rows,
                          data->gridData ) ) );

if ( result == 0 ) abort(); /* TEMP HACK */

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
Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/
2008-07-03T00:00:00-0000 2008-07-03T23:59:59-0000
# Dimensions: variables points
4 72
# Variable names:
timestamp longitude latitude pm25
# Variable units:
yyyymmddhhmmss deg deg metric_tons
# IEEE-754 64-bit reals data[variables][points]:
******************************************************************************/

static Integer readXDR( Stream* input, Data* data ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         data, data->points == 0, data->data == 0 );

  input->readString( input, data->note, COUNT( data->note ) );
  Integer result = input->ok( input );

  if ( result ) {
    removeTrailingNewline( data->note );
    result = readTimestamps( input, data->startingTimestamp,
                             data->endingTimestamp );

    if ( result ) {
      Integer dimensions[ 2 ] = { 0, 0 };
      data->timesteps =
        hoursInRange( data->startingTimestamp, data->endingTimestamp );
      result = readDimensions( input, COUNT( dimensions ), dimensions );

      if ( result ) {
        data->variables = dimensions[ 0 ];
        data->points    = dimensions[ 1 ];
        data->variable = NEW_ZERO( Name, data->variables * 2 );
        result = data->variables != 0;

        if ( result ) {
          data->units = data->variable + data->variables;
          result = readVariablesAndUnits( input, data->variables,
                                          data->variable, data->units );

          if ( AND2( result, data->variables >= 4 ) ) {
            result =
              OR2( AND3( ! strcmp( data->variable[ 0 ], "Timestamp" ),
                         ! strcmp( data->variable[ 1 ], "Longitude" ),
                         ! strcmp( data->variable[ 2 ], "Latitude" ) ),
                   AND3( ! strcmp( data->variable[ 0 ], "timestamp" ),
                         ! strcmp( data->variable[ 1 ], "longitude" ),
                         ! strcmp( data->variable[ 2 ], "latitude" ) ) );

            if ( result ) {
              char line[ 80 ] = "";
              memset( line, 0, sizeof line );
              input->readString( input, line, sizeof line / sizeof *line );
              result = input->ok( input );

              if ( result ) {

                if ( ! strcmp( line, "# char notes[points][80] and\n" ) ) {
                  data->notes = NEW_ZERO( Note, data->points );
                  result = data->notes != 0;

                  if ( result ) {
                    input->readString( input, line, sizeof line / sizeof *line);
                    result = input->ok( input );
                  }
                }

                if ( result ) {
                  result = ! strcmp( line,
                          "# IEEE-754 64-bit reals data[variables][points]:\n" );

                  if ( result ) {
                    result = readXDRData( input, data );
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
    failureMessage( "Invalid Point header." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidData( data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read XDR-format data and initialize data.
INPUTS:  Stream* input  Input stream read from and parse.
         Data* data     data->variables,timesteps,stations.
OUTPUTS: Data* data     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:  Data looks like this:
# IEEE-754 64-bit reals data[variables][points]:
or:
# char notes[points][80] and
# IEEE-754 64-bit reals data[variables][points]:
******************************************************************************/

static Integer readXDRData( Stream* input, Data* data ) {

  PRE09( input, input->ok( input ), input->isReadable( input ),
         data, data->variables > 3, data->points > 0,
         data->variable, data->units, data->data == 0 );

  Integer result = 0;
  const Integer count = data->variables * data->points;
  data->data = NEW_ZERO( Real, count );
  result = data->data != 0;

  if ( result ) {

    if ( data->notes ) {
      readNotes( input, data->points, data->notes );
      result = input->ok( input );
    }

    if ( result ) {
      input->read64BitReals( input, data->data, count );
      result = input->ok( input );

      if ( result ) {
        result = isValidData( data );
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Point data." );
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

REGRIDDED-Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/,XDRConvert
2008-07-03T00:00:00-0000
# timesteps
24
# Variable name:
pm25
# Variable units:
metric_tons
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000.000000 6370000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
459 299 -2556000.000000 -1728000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 32-bit integers points[timesteps] and
# IEEE-754 32-bit reals longitudes[timesteps][points] and
# IEEE-754 32-bit reals latitudes[timesteps][points] and
# MSB 32-bit integers columns[timesteps][points] and
# MSB 32-bit integers rows[timesteps][points] and
# IEEE-754 32-bit reals data[timesteps][points]:

******************************************************************************/

static Integer readRegriddedXDR( Stream* input, Data* data ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         data, data->points == 0, data->data == 0 );

  Integer result = 0;
  input->readString( input, data->note, COUNT( data->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( data->note );

    if ( readTimestamp( input, data->startingTimestamp ) ) {

      if ( readDimensions( input, 1, &data->timesteps ) ) {

        if ( readVariablesAndUnits2( input, data ) ) {

          if ( skipInputLines( input, 10 ) ) {
            data->outputPoints = NEW_ZERO( Integer, data->timesteps );

            if ( data->outputPoints ) {
              input->read32BitIntegers( input, data->outputPoints,
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
    failureMessage( "Invalid REGRIDDED-Point data." );
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
# MSB 32-bit integers points[timesteps] and
# IEEE-754 32-bit reals longitudes[timesteps][points] and
# IEEE-754 32-bit reals latitudes[timesteps][points] and
# MSB 32-bit integers columns[timesteps][points] and
# MSB 32-bit integers rows[timesteps][points] and
# IEEE-754 32-bit reals data[timesteps][points]:
******************************************************************************/

static Integer readRegriddedXDRData( Stream* input, Data* data ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         data, data->points == 0, data->data == 0, data->gridData == 0,
         data->timesteps > 0 );

  Integer result = 0;
  const Integer count = data->totalRegriddedPoints;

  if ( count > 0 ) {
    const Integer isVector = isVectorVariable( data );
    data->gridLongitudes = NEW_ZERO( Real, count * ( 5 + isVector ) );

    if ( data->outputPoints ) {
      data->gridLatitudes = data->gridLongitudes + count;
      data->columns   = (Integer*) data->gridLatitudes + count;
      data->rows      = data->columns + count;
      data->gridData  = (Real*) data->rows + count;

      input->read32BitReals( input, data->gridLongitudes, count * 2 );

      if ( input->ok( input ) ) {
        input->read32BitIntegers( input, data->columns, count * 2 );

        if ( input->ok( input ) ) {
          const size_t count2 = isVector ? count + count : count;
          input->read32BitReals( input, data->gridData, count2 );

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
PURPOSE: readVariablesAndUnits2 - Read 1 (e.g., ozone) or 2 (windU windV)
         sets of variables and units.
INPUTS:  Stream* input Input stream read from and parse.
OUTPUTS: Data* data    Allocated data->variables,variable,units.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
NOTES:   Input data format is:
# Variable name:
windU windV
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
    failureMessage( "Invalid Point header (variables/units)." );
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
    const Integer* const pointLayers  = data->layers;
    const Integer* const pointRows    = data->rows;
    const Integer* const pointColumns = data->columns;
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
        const Integer pointLayer  = pointLayers ? pointLayers[ pointIndex ] : 1;
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
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer hasElevation =
      OR2( ! strcmp( data->variable[ 3 ], "Elevation" ),
           ! strcmp( data->variable[ 3 ], "elevation" ) );
    const char* const headerStart =
      hasElevation ?
        "Timestamp(UTC)\tLongitude(deg)\tLatitude(deg)\tElevation(m)"
      : "Timestamp(UTC)\tLongitude(deg)\tLatitude(deg)";
    Integer variable = 0;

    /* Write header row: */

    output->writeString( output, headerStart );

    for ( variable = 3 + hasElevation;
          AND2( output->ok( output ), variable < data->variables );
          ++variable ) {
        const char* const headerFormat = "\t%s(%s)";
        output->writeString( output, headerFormat,
                             data->variable[ variable ],
                             data->units[ variable ] );
    }

    if ( output->ok( output ) ) {

      if ( data->notes ) {
        output->writeString( output, "\tNotes(-)\n" );
      } else {
        output->writeString( output, "\n" );
      }

      /* Write data rows: */

      if ( output->ok( output ) ) {
        const Integer pointCount = data->points;
        const Real* const timestamps = data->data;
        const Real* const longitudes = timestamps + pointCount;
        const Real* const latitudes  = longitudes + pointCount;
        const Real* const elevations = hasElevation ? latitudes + pointCount :0;
        Integer pointIndex = 0;

        for ( pointIndex = 0;
              AND2( output->ok( output ), pointIndex < pointCount );
              ++pointIndex ) {
          const Integer yyyymmddhhmmss = (Integer) timestamps[ pointIndex ];
          const Real elevation = elevations ? elevations[ pointIndex ] : 0.0;
          const Real longitude = longitudes[ pointIndex ];
          const Real latitude  = latitudes[  pointIndex ];
          UTCTimestamp timestamp;
          toUTCTimestamp2( yyyymmddhhmmss, timestamp );

          output->writeString( output, "%s\t%10.5lf\t%10.5lf",
                              timestamp, longitude, latitude );

          if ( output->ok( output ) ) {

            if ( hasElevation ) {
              output->writeString( output, "\t%10.5lf", elevation );
            }

            for ( variable = 3 + hasElevation;
                 AND2( output->ok( output ), variable < data->variables );
                 ++variable ) {
              const Integer index = variable * pointCount + pointIndex;
              const Real value = data->data[ index ];
              output->writeString( output, "\t%10.5lf", value );
            }

            if ( output->ok( output ) ) {

              if ( data->notes ) {
                output->writeString( output, "\t%-80s\n",
                                     data->notes[ pointIndex ] );
              } else {
                output->writeString( output, "\n" );
              }
            }
          }
        }
      }
    }

    result = output->ok( output );
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
  const Integer hasNotes = data->notes != 0;
  const Integer fileSizeEstimate =
    data->points * data->variables * 4 +
    data->points * hasNotes * sizeof (Note) + 10000; /* header. */
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
  const Integer hasNotes = data->notes != 0;
  const char* const names[ 2 ] = { "points", "length" };
  Integer dimensionIds[ 2 ] = { -1, -1 };
  Integer sizes[ 2 ] = { 0, 0 };
  sizes[ 0 ] = data->points;

  if ( hasNotes ) {
    sizes[ 1 ] = sizeof (Note) / sizeof (char);
  }

  if ( createDimensions( file, 1 + hasNotes, names, sizes, dimensionIds ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createLongitudeAndLatitude( file, 1, dimensionIds ) ) {
        int variable = 0;
        int ok = 1;

        for ( variable = 3; AND2(ok, variable < data->variables); ++variable) {
          const int isInteger =
            OR2( AND2( ! strcmp( data->variable[ variable ], "id" ),
                       ! strcmp( data->units[ variable ], "-" ) ),
                 AND2( ! strcmp( data->variable[ variable ], "count" ),
                       ! strcmp( data->units[ variable ], "-" ) ) );
          ok = createVariable( file, data->variable[ variable ],
                               data->units[ variable ],
                               isInteger ? NC_INT : NC_FLOAT,
                               1, 1, dimensionIds ) != -1;
        }

        if ( ok ) {

          if ( OR2( ! hasNotes,
                    createVariable( file, "notes", "-",
                                    NC_CHAR, 0, 2, dimensionIds ) != -1 ) ) {
            const Integer points = data->points;
            const Real* longitudes = data->data + points;
            const Real* latitudes = longitudes  + points;
            Bounds bounds;
            computeBounds( points, longitudes, latitudes, bounds );

            if ( writeExtraAttributes( file, (const Real (*)[2]) bounds,
                                       dimensionIds[ 0 ] ) ) {
              Line history = "";
              appendToLine( history, data->note );
              appendToLine( history, ",XDRConvert" );
              result = writeStandardContents( file, history,
                                              data->startingTimestamp,
                                              dimensionIds[ 0 ],
                                              data->points, 0 );
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
  const Integer count = data->points;
  Real* const timestamps = data->data;
  Real* const longitudes = timestamps + count;
  Real* const latitudes  = longitudes + count;

  result = writeAllData( file, "longitude", count, 1, 1, 1, longitudes );

  if ( result ) {
    result = writeAllData( file, "latitude", count, 1, 1, 1, latitudes );

    if ( result ) {
      Real* values = latitudes + count;
      int variable = 0;
      result = 1;

      for ( variable = 3; AND2( result, variable < data->variables );
           ++variable, values += count ) {
        const int isInteger =
          OR2( AND2( ! strcmp( data->variable[ variable ], "id" ),
                     ! strcmp( data->units[ variable ], "-" ) ),
               AND2( ! strcmp( data->variable[ variable ], "count" ),
                     ! strcmp( data->units[ variable ], "-" ) ) );

        if ( isInteger ) { /* Convert 64-reals to 64-bit integers in-place: */
          Integer* ivalues = (Integer*) values;
          Integer index = 0;

          for ( index = 0; index < count; ++index ) {
            const Real value = values[ index ];
            const Integer ivalue = (Integer) value;
            ivalues[ index ] = ivalue;
          }

          result = writeAllIntData( file, data->variable[ variable ], count,
                                    1, 1, 1, ivalues );
        } else {
          result = writeAllData( file, data->variable[ variable ], count,
                               1, 1, 1, values );
        }
      }

      if ( result ) {
        result = writeConvertedTimeData( file, data->startingTimestamp,
                                         count, timestamps );

        if ( result ) {

          if ( data->notes ) {
            result = writeAllCharData( file, "notes",
                                       count, sizeof (Note) / sizeof (char),
                                       data->notes[0] );
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeConvertedTimeData - Write COARDS-format time data to file.
INPUTS:  const Integer file            NetCDF file to write to.
         const UTCTimestamp startingTimestamp  Starting timestamp of subset.
         const Integer count           Number of timestamps.
         const Real* const timestamps  Timestamps to convert and write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeConvertedTimeData( const Integer file,
                                       const UTCTimestamp startingTimestamp,
                                       const Integer count,
                                       const Real* const timestamps ) {
  PRE06( file != -1, isValidUTCTimestamp( startingTimestamp ),
         count > 0, timestamps,
         isValidYYYYMMDDHHMMSS( timestamps[ 0 ] ),
         isValidYYYYMMDDHHMMSS( timestamps[ count - 1 ] ) );

  int* yyyyddd = NEW_ZERO( int, count + count );
  float* fhour = yyyyddd ? NEW_ZERO( float, count ) : 0;
  Integer result = 0;

  if ( fhour != 0 ) {
    int* hhmmss = yyyyddd + count;
    const Integer yyyydddhhmmStart = fromUTCTimestamp( startingTimestamp );
    Integer index = 0;

    for ( index = 0; index < count; ++index ) {
      const Real timestamp = timestamps[ index ];
      const Integer yyyymmddhhmmss = (Integer) timestamp;

      if ( ! isValidYYYYMMDDHHMMSS( yyyymmddhhmmss ) ) {
        failureMessage( "Invalid timestamp read: %lld.", yyyymmddhhmmss );
        index = count;
      } else {
        const Integer yyyymmdd = yyyymmddhhmmss / 1000000;
        const Integer hhmmss0  = yyyymmddhhmmss % 1000000;
        CHECK2( isValidYearMonthDay( yyyymmdd ), isValidTime( hhmmss0 ) );

        {
          const Integer yyyyddd0 = convertYearMonthDay( yyyymmdd );
          CHECK( isValidDate( yyyyddd0 ) );
          yyyyddd[ index ] = (int) yyyyddd0;
          hhmmss[  index ] = (int) hhmmss0;

          {
            const Integer hhmm = hhmmss0 / 100;
            const Integer yyyydddhhmm =  yyyyddd0 * 10000 + hhmm;
            CHECK( isValidTimestamp( yyyydddhhmm ) );

            {
              const Real fraction =
                fractionalHours( yyyydddhhmmStart, yyyydddhhmm );
              fhour[ index ] = (float) fraction;
            }
          }
        }
      }
    }

    if ( index == count ) {
      result = writeTimeData1( file, count, yyyyddd, hhmmss, fhour );
    }
  }

  FREE( yyyyddd );
  FREE( fhour );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  const Data* data              Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Output data format is:

REGRIDDED-Point 1.0
ftp://ftp.orbit.nesdis.noaa.gov/pub/smcd/xzhang/PM25/,XDRConvert
2008-07-03T00:00:00-0000
# timesteps
24
# Variable name:
pm25
# Variable units:
metric_tons
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000.000000 6370000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
459 299 -2556000.000000 -1728000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 32-bit integers points[timesteps] and
# IEEE-754 32-bit reals longitudes[timesteps][points] and
# IEEE-754 32-bit reals latitudes[timesteps][points] and
# MSB 32-bit integers columns[timesteps][points] and
# MSB 32-bit integers rows[timesteps][points] and
# IEEE-754 32-bit reals data[timesteps][points]:

******************************************************************************/

static Integer writeRegriddedXDR( const Data* data,
                                  const Parameters* parameters ) {

  PRE02( isValidData( data ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = data->timesteps;
    const Integer hasElevation =
      OR2( ! strcmp( data->variable[ 3 ], "Elevation" ),
           ! strcmp( data->variable[ 3 ], "elevation" ) );
    const Integer isVector  = isVectorVariable( data );
    const Integer variableIndex0 = data->variables - 1 - isVector;
    const Integer variableIndex =
      variableIndex0 >= 0 ? variableIndex0 : 0;
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
    Name variable = "";
    aggregateName( data->variable[ variableIndex ], hoursPerTimestep, variable);

    if ( isVector ) {
     Name variable2 = "";
      aggregateName( data->variable[ variableIndex + 1  ], hoursPerTimestep,
                     variable2 );
      output->writeString( output,
                           "REGRIDDED-Point 1.0\n"
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
                           data->units[ variableIndex ],
                           data->units[ variableIndex + 1 ] );
    } else {
      output->writeString( output,
                           "REGRIDDED-Point 1.0\n"
                           "%s,XDRConvert\n"
                           "%s\n"
                           "# timesteps\n%"INTEGER_FORMAT"\n"
                           "# Variable name:\n%s\n"
                           "# Variable units:\n%s\n",
                           data->note,
                           data->startingTimestamp,
                           timesteps,
                           variable,
                           data->units[ variableIndex ] );
    }

    if ( output->ok( output ) ) {
      writeProjectionAndGrid( parameters->grid, output );

      if ( output->ok( output ) ) {

        if ( hasElevation ) {
          output->writeString( output,
                  "# MSB 32-bit integers points[timesteps] and\n"
                  "# IEEE-754 32-bit reals longitudes[timesteps][points] and\n"
                  "# IEEE-754 32-bit reals latitudes[timesteps][points] and\n"
                  "# IEEE-754 32-bit reals elevations[timesteps][points] and\n"
                  "# MSB 32-bit integers columns[timesteps][points] and\n"
                  "# MSB 32-bit integers rows[timesteps][points] and\n"
                  "# MSB 32-bit integers layers[timesteps][points] and\n"
                  "# IEEE-754 32-bit reals data[timesteps][points]:\n" );
        } else {
          output->writeString( output,
                  "# MSB 32-bit integers points[timesteps] and\n"
                  "# IEEE-754 32-bit reals longitudes[timesteps][points] and\n"
                  "# IEEE-754 32-bit reals latitudes[timesteps][points] and\n"
                  "# MSB 32-bit integers columns[timesteps][points] and\n"
                  "# MSB 32-bit integers rows[timesteps][points] and\n"
                  "# IEEE-754 32-bit reals data[timesteps][points]:\n" );
        }

        if ( output->ok( output ) ) {
          output->write32BitIntegers( output, data->outputPoints, timesteps );

          if ( output->ok( output ) ) {
            const Integer points = data->totalRegriddedPoints;
            output->write32BitReals( output, data->gridLongitudes, points );

            if ( output->ok( output ) ) {
              output->write32BitReals( output, data->gridLatitudes, points );

              if ( output->ok( output ) ) {

                if ( hasElevation ) {
                  output->write32BitReals(output, data->gridElevations, points);
                }

                if ( output->ok( output ) ) {
                  output->write32BitIntegers( output, data->columns, points );

                  if ( output->ok( output ) ) {
                    output->write32BitIntegers( output, data->rows, points );

                    if ( output->ok( output ) ) {

                      if ( hasElevation ) {
                        output->write32BitIntegers(output, data->layers,points);
                      }

                      if ( output->ok( output ) ) {
                        const Integer points2 =
                          isVector ? points + points : points;
                        output->write32BitReals(output, data->gridData,points2);
                        result = output->ok( output );
                      }
                    }
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
    const Integer hasElevation =
      OR2( ! strcmp( data->variable[ 3 ], "Elevation" ),
           ! strcmp( data->variable[ 3 ], "elevation" ) );
    const Integer isVector = isVectorVariable( data );
    const Integer variableIndex0 = data->variables - 1 - isVector;
    const Integer variableIndex =
      variableIndex0 >= 0 ? variableIndex0 : 0;
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
    const char* const headerStart =
      hasElevation ?
        "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tELEVATION(m)"
        "\tCOLUMN(-)\tROW(-)\tLAYER(-)"
      : "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)"
        "\tCOLUMN(-)\tROW(-)";
    Name variable = "";
    aggregateName( data->variable[ variableIndex ], hoursPerTimestep, variable);

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {

      if ( isVector ) {
        const char* const headerFormat = "\t%s(%s)\t%s(%s)\n";
        Name variable2 = "";
        aggregateName( data->variable[ variableIndex + 1 ], hoursPerTimestep,
                       variable2 );

        output->writeString( output, headerFormat,
                             variable,  data->units[ variableIndex ],
                             variable2, data->units[ variableIndex + 1 ] );
      } else {
        const char* const headerFormat = "\t%s(%s)\n";
        output->writeString( output, headerFormat, variable, data->units );
      }

      if ( output->ok( output ) ) {
        const Integer timesteps = data->timesteps;
        const Real* longitudes  = data->gridLongitudes;
        const Real* latitudes   = data->gridLatitudes;
        const Real* elevations  = data->gridElevations;
        const Integer* columns  = data->columns;
        const Integer* rows     = data->rows;
        const Integer* layers   = data->layers;
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

              if ( hasElevation ) {
                const Real elevation = *elevations++;
                const Integer layer  = *layers++;
                const char* const dataFormat =
                  "%s\t%10.5lf\t%10.5lf\t%10.5lf"
                  "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT
                  "\t%9"INTEGER_FORMAT"\t%10.5lf\t%10.5lf\n";
                output->writeString( output, dataFormat,
                                     timestamp, longitude, latitude, elevation,
                                     column, row, layer, value, value2 );
              } else {
                const char* const dataFormat =
                  "%s\t%10.5lf\t%10.5lf"
                  "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT
                  "\t%10.5lf\t%10.5lf\n";
                output->writeString( output, dataFormat,
                                     timestamp, longitude, latitude,
                                     column, row, value, value2 );
              }
            } else {

              if ( hasElevation ) {
                const Real elevation = *elevations++;
                const Integer layer  = *layers++;
                const char* const dataFormat =
                  "%s\t%10.5lf\t%10.5lf\t%10.5lf"
                  "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT
                  "\t%9"INTEGER_FORMAT"\t%10.5lf\n";
                output->writeString( output, dataFormat,
                                     timestamp, longitude, latitude, elevation,
                                     column, row, layer, value );
              } else {
                const char* const dataFormat =
                  "%s\t%10.5lf\t%10.5lf"
                  "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%10.5lf\n";
                output->writeString( output, dataFormat,
                                     timestamp, longitude, latitude,
                                     column, row, value );
              }
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
      result = writeRegriddedCOARDSData( file, data, parameters );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDSHeader - Write header to file.
INPUTS:  Integer file               NetCDF file to write to.
         Integer hoursPerTimestep   E.g., 1, 24, 744, etc.
         const Data* data           Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Data* data ) {

  PRE03( file != -1, hoursPerTimestep > 0, isValidData( data ) );

  Integer result = 0;
  const char* const dimensionName = "points";
  Integer dimensionId = -1;
  const Integer dimension = data->totalRegriddedPoints;

  if ( createDimensions( file, 1, &dimensionName, &dimension, &dimensionId)) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "column", "-",
                           NC_INT, 0, 1, &dimensionId ) != -1 ) {

        if ( createVariable( file, "row", "-",
                             NC_INT, 0, 1, &dimensionId ) != -1 ) {
          const Integer hasElevation = data->gridElevations != 0;

          if ( OR2( ! hasElevation,
                    createVariable( file, "layer", "-",
                                    NC_INT, 0, 1, &dimensionId ) != -1 ) ) {

            if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {

              if ( OR2( ! hasElevation,
                        createVariable( file, "elevation", "m",
                                        NC_INT, 0, 1, &dimensionId ) != -1 ) ) {
                const Integer isVector = isVectorVariable( data );
                const Integer variableIndex0 = data->variables - 1 - isVector;
                const Integer variableIndex =
                  variableIndex0 >= 0 ? variableIndex0 : 0;
                Name variable = "";
                aggregateName( data->variable[ variableIndex ], hoursPerTimestep,
                               variable );

                if ( createVariable( file, variable,
                                     data->units[ variableIndex ],
                                     NC_FLOAT, 1, 1, &dimensionId ) != -1 ) {

                  result = 1;

                  if ( isVector ) {
                    Name variable2 = "";
                    aggregateName( data->variable[ variableIndex + 1 ],
                                   hoursPerTimestep, variable2 );
                    result =
                      createVariable( file, variable2,
                                      data->units[ variableIndex + 1 ],
                                      NC_FLOAT, 1, 1, &dimensionId ) != -1;
                  }

                  if ( result ) {
                    Line history = "";
                    appendToLine( history, data->note );
                    appendToLine( history, ",XDRConvert" );
                    result = writeStandardContents( file, history,
                                                    data->startingTimestamp,
                                                    dimensionId, 0, 0 );
                  }
                }
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

  if ( writeAllIntData( file, "column", count, 1, 1, 1, data->columns ) ) {

    if ( writeAllIntData( file, "row", count, 1, 1, 1, data->rows ) ) {
      const Integer hasElevation = data->gridElevations != 0;

      if ( OR2( ! hasElevation,
                writeAllIntData( file, "layer", count, 1,1,1, data->layers ))) {

        if ( writeAllData( file, "longitude", count, 1, 1, 1,
                           data->gridLongitudes ) ) {

          if ( writeAllData( file, "latitude", count, 1, 1, 1,
                             data->gridLatitudes ) ) {

            if ( OR2( ! hasElevation,
                      writeAllData( file, "elevation", count, 1, 1, 1,
                                    data->gridElevations ) ) ) {
              const Integer isVector = isVectorVariable( data );
              const Integer variableIndex0 = data->variables - 1 - isVector;
              const Integer variableIndex =
                variableIndex0 >= 0 ? variableIndex0 : 0;
              const Integer hoursPerTimestep =
                parameters->aggregationTimesteps ?
                  parameters->aggregationTimesteps
                : 1;
              Name variable = "";
              aggregateName( data->variable[ variableIndex ], hoursPerTimestep,
                             variable );

              if ( writeAllData( file, variable, count, 1, 1, 1,
                                 data->gridData ) ) {
                const Integer isVector = isVectorVariable( data );
                result = 1;

                if ( isVector ) {
                  Name variable2 = "";
                  aggregateName( data->variable[ variableIndex + 1 ],
                                 hoursPerTimestep, variable2 );
                  result =
                    writeAllData( file, variable2, count, 1,1,1,
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
      result =
        writeRegriddedIOAPIData(file, hoursPerTimestep, data, parameters->grid);
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
  const Integer hasElevation = grid->layers( grid ) > 1;
  const Integer layers = hasElevation ? grid->layers( grid ) : 1;
  enum { VARIABLES = 3 }; /* LONGITUDE, LATITUDE, var. */
  Name variableNames[ VARIABLES + 2 ] =
    { "LONGITUDE", "LATITUDE", "ELEVATION", "var", "WindV" };
  Name variableUnits[ VARIABLES + 2 ] = { "deg", "deg", "m", "m/s", "m/s" };
  const Integer firstTimestamp = fromUTCTimestamp( data->startingTimestamp );
  Integer variable = VARIABLES - 1 + hasElevation;
  const Integer isVector = isVectorVariable( data );
  const Integer variableIndex0 = data->variables - 1 - isVector;
  const Integer variableIndex = variableIndex0 >= 0 ? variableIndex0 : 0;
  Line history = "";
  appendToLine( history, data->note );
  appendToLine( history, ",XDRConvert" );
  aggregateName( data->variable[ variableIndex ], hoursPerTimestep,
                 variableNames[ variable ] );
  variableNames[ variable ][ 15 ] = '\0';
  strncpy( variableUnits[ variable ], data->units[ variableIndex ], 16 );
  variableUnits[ variable ][ 16 ] = '\0';
  uppercase( variableNames[ variable ] );
  lowercase( variableUnits[ variable ] );

  if ( isVector ) {
    ++variable;
    aggregateName( data->variable[ variableIndex + 1 ], hoursPerTimestep,
                   variableNames[ variable ] );
    variableNames[ variable ][ 15 ] = '\0';
    strncpy( variableUnits[ variable ], data->units[ variableIndex + 1 ], 16 );
    variableUnits[ variable ][ 16 ] = '\0';
    uppercase( variableNames[ variable ] );
    lowercase( variableUnits[ variable ] );
  }

  result = writeM3IOHeader( file, data->timesteps, hoursPerTimestep,
                            firstTimestamp,
                            VARIABLES + hasElevation + isVector,
                            layers,
                            (const Name*) variableNames,
                            (const Name*) variableUnits,
                            history, grid );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIData - Write IOAPI-format data to file.
INPUTS:  Integer file      NetCDF file to write to.
         Integer hoursPerTimestep   E.g., 1, 24, 744, etc.
         const Data* data  Structure to write.
         const Grid* grid  Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Data* data,
                                        const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0, isValidData( data ),
         grid, grid->invariant( grid ) );

  Integer result = 0;
  const Integer hasElevation = data->layers != 0;
  const Integer isVector = isVectorVariable( data );
  const Integer layers  = hasElevation ? grid->layers( grid ) : 1;
  const Integer rows    = grid->rows( grid );
  const Integer columns = grid->columns( grid );
  const Integer cells   = layers * rows * columns;
  Real* gridData = NEW_ZERO( Real, cells  );

  if ( gridData ) {
    const Integer timesteps = data->timesteps;
    Integer offset = 0;
    Integer timestep = 0;

    if ( writeM3IOGrid( grid, timesteps, layers, file ) ) {
      const Integer variableIndex0 = data->variables - 1 - isVector;
      const Integer variableIndex = variableIndex0 >= 0 ? variableIndex0 : 0;
      Name variable = "";
      Name variable2 = "";
      aggregateName( data->variable[ variableIndex ], hoursPerTimestep,
                     variable );
      variable[ 15 ] = '\0';
      uppercase( variable );

      if ( isVector ) {
        aggregateName( data->variable[ variableIndex + 1 ], hoursPerTimestep,
                       variable2 );
        variable2[ 15 ] = '\0';
        uppercase( variable2 );
      }

      do {
        const Integer points = data->outputPoints[ timestep ];

        if ( hasElevation ) {
          copyDataToGrid3( points,
                           data->layers   + offset,
                           data->rows     + offset,
                           data->columns  + offset,
                           data->gridData + offset,
                           1.0, layers, rows, columns, gridData );
        } else {
          copyDataToGrid( points,
                          data->rows     + offset,
                          data->columns  + offset,
                          data->gridData + offset,
                          1.0, layers, rows, columns, gridData );
        }

        if ( ! writeM3IOData( file, variable,
                              timestep, layers, rows, columns, gridData ) ) {
          timestep = timesteps;
        } else {

          if ( isVector ) {
            const Real* gridData2 =
              data->gridData + data->totalRegriddedPoints;

            if ( hasElevation ) {
              copyDataToGrid3( points,
                               data->layers    + offset,
                               data->rows      + offset,
                               data->columns   + offset,
                               gridData2       + offset,
                               1.0, layers, rows, columns, gridData );
            } else {
              copyDataToGrid( points,
                              data->rows    + offset,
                              data->columns + offset,
                              gridData2     + offset,
                              1.0, layers, rows, columns, gridData );
            }

            if ( ! writeM3IOData( file, variable2, timestep, layers, rows,
                                  columns, gridData ) ) {
              timestep = timesteps;
            }
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
  const Integer hasElevation = ! strcmp( data->variable[ 3 ], "elevation" );
  const Integer isVector = isVectorVariable( data );
  const Integer timesteps = data->timesteps;
  const Integer points = data->points;
  data->outputPoints = NEW_ZERO( Integer, timesteps );

  if ( data->outputPoints ) {
    const Integer outputVariables = 5 + hasElevation * 2 + isVector;
    /* lon,lat,col,row,var + elv, lay, + var2 */
    const Integer outputSize = outputVariables * points;
    data->gridLongitudes = NEW_ZERO( Real, outputSize );

    if ( data->gridLongitudes ) {
      const Integer yyyydddhhmm = fromUTCTimestamp( data->startingTimestamp );
      Integer yyyydddhh00 = yyyydddhhmm / 100 * 100;
      Integer timestep = 0;
      Integer index = 0;
      const Real* const longitudes = data->data + points;
      const Real* const latitudes  = longitudes + points;
      const Real* const elevations = hasElevation ? latitudes + points : 0;
      const Real* const values =
        data->data + ( data->variables - 1 - isVector ) * points;
      const Real* const values2 = isVector ? values + points : 0;
      Real* gridData2 = 0;
      data->gridLatitudes = data->gridLongitudes + points;
      data->gridElevations = hasElevation ? data->gridLatitudes + points : 0;
      data->columns = (Integer*)
        ( hasElevation ? data->gridElevations + points
        : data->gridLatitudes + points );
      data->rows = data->columns + points;
      data->layers = hasElevation ? data->rows + points : 0;
      data->gridData = (Real*)
        ( hasElevation ? data->layers + points : data->rows + points );

      if ( isVector ) {
        gridData2 = data->gridData + points;
      }

      for ( timestep = 0;
            AND2( timestep < timesteps, index < points );
            ++timestep ) {
        const Integer yyyydddhh = yyyydddhh00 / 100;
        const Integer count = pointsMatchingHour( yyyydddhh, data, index );
        Integer outputPoints = 0;

        if ( count ) {
          const Real minimumValidValue = -500.0;

          grid->regrid( grid, method, minimumValidValue, count, 1,
                        longitudes + index, latitudes + index,
                        hasElevation ? elevations + index : 0,
                        values + index,
                        isVector ? values2 + index : 0,
                        0, /* No notes. */
                        &outputPoints,
                        data->columns        + totalRegriddedPoints,
                        data->rows           + totalRegriddedPoints,
                        hasElevation ? data->layers + totalRegriddedPoints : 0,
                        data->gridLongitudes + totalRegriddedPoints,
                        data->gridLatitudes  + totalRegriddedPoints,
                        hasElevation ?
                          data->gridElevations + totalRegriddedPoints : 0,
                        data->gridData       + totalRegriddedPoints,
                        isVector ? gridData2 + totalRegriddedPoints : 0,
                        0 /* No regriddedNotes. */ );

          index += count;
        }

        data->outputPoints[ timestep ] = outputPoints;
        totalRegriddedPoints += outputPoints;
        incrementTimestamp( &yyyydddhh00 );
      }
    }
  }

  data->totalRegriddedPoints = totalRegriddedPoints;

  POST02( data->totalRegriddedPoints >= 0,
          IMPLIES( data->totalRegriddedPoints > 0,
                   AND8( IN_RANGE( minimumItemI(data->columns,
                                                data->totalRegriddedPoints) ,
                                                1, grid->columns( grid ) ),
                         IN_RANGE( maximumItemI(data->columns,
                                                data->totalRegriddedPoints ),
                                                1, grid->columns( grid ) ),
                         IN_RANGE( minimumItemI(data->rows,
                                                data->totalRegriddedPoints ),
                                                1, grid->rows( grid ) ),
                         IN_RANGE( maximumItemI(data->rows,
                                                data->totalRegriddedPoints ),
                                                1, grid->rows( grid ) ),
                         IMPLIES( data->layers,
                                  IN_RANGE( minimumItemI(data->layers,
                                                data->totalRegriddedPoints ),
                                                1, grid->layers( grid ) ) ),
                         IMPLIES( data->layers,
                                  IN_RANGE( maximumItemI(data->layers,
                                                data->totalRegriddedPoints ),
                                                1, grid->layers( grid ) ) ),
                         validLongitudesAndLatitudes(data->totalRegriddedPoints,
                                                      data->gridLongitudes,
                                                      data->gridLatitudes ),
                         isNanFree( data->gridData,
                                    data->totalRegriddedPoints
                                    * ( 1 + isVector ) ) ) ) );
}



/******************************************************************************
PURPOSE: pointsMatchingHour - Count number of data points match the given hour.
INPUTS:  const Integer yyyydddhh     Hour to match.
         const Data* const data      Data to check.
         const Integer start         Index to start search.
RETURNS: Integer number of poitms matching hour.
NOTES:   Requires data->data timestamps in form yyyymmddhhmmss and sorted.
******************************************************************************/

static Integer pointsMatchingHour( const Integer yyyydddhh,
                                   const Data* const data,
                                   const Integer start ) {

  PRE05( isValidTimestamp( yyyydddhh * 100 ),
         data, data->data, data->variables > 0,
         IN_RANGE( start, 0, data->points - 1 ) );

  Integer result = 0;
  const Integer points = data->points;
  const Real* const timestamps = data->data;
  Integer point = start;

  for ( point = 0; point < points; ++point ) {
    const Real timestamp = timestamps[ point ];
    const Integer yyyymmddhhmmss = (Integer) timestamp;

    if ( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss ) ) {
      const Integer hh = yyyymmddhhmmss / 10000 % 100;
      const Integer yyyymmdd = yyyymmddhhmmss / 1000000;
      const Integer yyyyddd = convertYearMonthDay( yyyymmdd );
      const Integer yyyydddhh0 = yyyyddd * 100 + hh;

      if ( yyyydddhh0 == yyyydddhh ) {
        ++result;
      } else if ( yyyydddhh0 > yyyydddhh ) {
        point = points; /* Stop looping. */
      }
    }
  }

  POST0( IN_RANGE( result, 0, data->points - start ) );
  return result;
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





/******************************************************************************
PURPOSE: Aircraft.c - Define routines for processing Aircraft data.

NOTES:   Input aircraft data in XDR format is as follows:

Aircraft 2.0
http://mozaic.aero.obs-mip.fr/web/,MOZAICSubset
2006-07-03T00:00:00-0000 2006-07-04T23:59:59-0000
# Subset bounds: <min_lon> <min_lat> <max_lon> <max_lat>:
-180 -90 180 90
# Dimensions: variables points tracks:
5 371 2
# Variable names:
 timestamp longitude latitude elevation ozone
# Variable units:
 yyyymmddhhmmss deg deg m ppmV
# char notes[tracks][80] and
# IEEE-754 64-bit reals bounds[tracks][2=lon,lat][2=min,max] and
# MSB 64-bit integers points[tracks] and
# IEEE-754 64-bit reals data_1[points_1][variables]
  ... data_T[points_T][variables]:


HISTORY: 2010/02 plessel.todd@epa.gov, Created.
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

/* MOZAIC variables must be > 4 and exactly 5 if regridding: */

enum {
  AIRCRAFT_TIMESTAMP, AIRCRAFT_LONGITUDE, AIRCRAFT_LATITUDE,AIRCRAFT_ELEVATION,
  IMPLICIT_VARIABLES = 4
};

typedef struct {
  /* Input data: */
  Line         description;    /* File description. */
  UTCTimestamp firstTimestamp; /* Earliest timestamp of data. */
  UTCTimestamp lastTimestamp;  /* Latest   timestamp of data. */
  Bounds       bounds;         /*bounds[LONGITUDE LATITUDE][MINIMUM MAXIMUM]*/
  Integer      variables;      /* E.g., 4+1 = timestamp,lon,lat,elv,ozone. */
  Integer      totalPoints;    /* Sum of points[ track ]. */
  Integer      tracks;         /* E.g., 2 aircraft tracks. */
  Note*        notes;          /* notes[ tracks ]. flight:from->to. */
  Name*        variable;       /* variable[ variables ]. E.g., "ozone"*/
  Name*        units;          /* units[ variables ]. E.g., "ppb". */
  Integer*     points;         /* points[ tracks ] */
  Real*        data;           /* data_1[ variables ][ points_1 ] ... */
                               /* data_T[ variables ][ points_T ] */
  /* Regridded data: */
  Integer      totalRegriddedPoints; /* Total number of regridded points. */
  Integer      timesteps;      /* Hours in regridded output. */
  Integer*     timestamps;     /* timestamps[ timesteps ]. */
  Integer*     outputPoints;   /* outputPoints[ timesteps ]. */
  Real*        longitudes;     /*longitudes[MIN(tracks,timesteps)*maxPoints]*/
  Real*        latitudes;      /*latitudes[ MIN(tracks,timesteps)*maxPoints]*/
  Real*        elevations;     /*elevations[MIN(tracks,timesteps)*maxPoints]*/
  Real*        gridLongitudes; /* gridLongitudes[ totalRegriddedPoints ]. */
  Real*        gridLatitudes;  /* gridLatitudes[ totalRegriddedPoints ]. */
  Real*        gridElevations; /* gridElevations[ totalRegriddedPoints ]. */
  Integer*     columns;        /* columns[ totalRegriddedPoints ]. */
  Integer*     rows;           /* rows[ totalRegriddedPoints ]. */
  Integer*     layers;         /* layers[ totalRegriddedPoints ]. */
  Real*        copyData;       /* copyData[ totalPoints ]. */
  Real*        gridData;       /* gridData[ totalRegriddedPoints ]. */
  Note*        copyNotes;      /* copyNotes[ totalPoints ]. */
  RegriddedNote* regriddedNotes; /* regriddedNotes[ totalRegriddedPoints ]. */
} Aircraft;

typedef Integer (*Writer)( Aircraft* aircraft, const Parameters* parameters );

typedef struct {
  Integer format;         /* FORMAT_XDR, etc. */
  Writer writer;          /* Routine that writes data in this format. */
  Writer regriddedWriter; /* Routine that writes regridded data in format. */
} Entry;

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateAircraft( Aircraft* aircraft );

static Integer isValidAircraft( const Aircraft* aircraft );

static Integer isVectorVariable( const Aircraft* aircraft );

static Writer dispatcher( Integer format, Integer regrid );

static Integer readXDR( Stream* input, Aircraft* aircraft );

static Integer readXDRData( Stream* input, Aircraft* aircraft );

static Integer readRegriddedXDR( Stream* input, Aircraft* aircraft );

static Integer readVariablesAndUnits2( Stream* input, Aircraft* aircraft );

static Integer readRegriddedXDRData( Stream* input, Aircraft* aircraft );

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    Aircraft* aircraft );

static Integer writeASCII( Aircraft* aircraft, const Parameters* parameters );

static void writeASCIIHeader( const Aircraft* aircraft, Stream* output );

static Integer writeASCIIData( Aircraft* aircraft, Stream* output );

static Integer writeCOARDS( Aircraft* aircraft, const Parameters* parameters );

static Integer writeCOARDSHeader( Integer file, const Aircraft* aircraft );

static Integer writeCOARDSData( Integer file, Aircraft* aircraft );


static Integer writeRegriddedXDR( Aircraft* aircraft,
                                  const Parameters* unused );

static Integer writeRegriddedASCII( Aircraft* aircraft,
                                    const Parameters* unused );

static Integer writeRegriddedCOARDS( Aircraft* aircraft,
                                     const Parameters* parameters );

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Aircraft* aircraft );

static Integer writeRegriddedCOARDSData( Integer file,
                                         const Aircraft* aircraft,
                                         const Parameters* parameters );

static Integer writeRegriddedIOAPI( Aircraft* aircraft,
                                    const Parameters* parameters) ;

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Aircraft* aircraft,
                                          const Grid* grid );

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Aircraft* aircraft,
                                        const Grid* grid );

static void regridAircraft( Integer method, Grid* grid, Aircraft* aircraft );


static Integer copyDataForTimestamp( Integer yyyydddhh00, Aircraft* aircraft );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateAircraft - Read input & write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void translateAircraft( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok,
         parameters->input->ok( parameters->input ) );

  Aircraft aircraft;
  ZERO_OBJECT( &aircraft );
  parameters->ok = 0;

  if ( readXDR( parameters->input, &aircraft ) ) {
    Writer writer = dispatcher( parameters->format, parameters->regrid );

    if ( ! writer ) {
      failureMessage( "Invalid/unsupported format/regrid specification." );
    } else if ( parameters->regrid ) {
      regridAircraft( parameters->regrid, parameters->grid, &aircraft );

      if ( aircraft.totalRegriddedPoints == 0 ) {
        failureMessage( "No points projected onto the grid." );
      } else {

        if ( parameters->aggregationTimesteps ) {
          const Integer dataVariable = aircraft.variables - 1;
          Integer totalOutputPoints = 0;
          const Integer aggregatedTimesteps =
            aggregateData( parameters->aggregationTimesteps,
                           0,
                           aircraft.timesteps,
                           aircraft.outputPoints,
                           aircraft.gridLongitudes,
                           aircraft.gridLatitudes,
                           aircraft.elevations,
                           aircraft.columns,
                           aircraft.rows,
                           aircraft.layers,
                           aircraft.gridData,
                           aircraft.regriddedNotes,
                           &totalOutputPoints );
          aircraft.timesteps = aggregatedTimesteps;
          aircraft.totalRegriddedPoints = totalOutputPoints;

          if ( AND2( parameters->aggregationTimesteps == 24,
                     ! OR2( strstr( aircraft.variable[ dataVariable ], "daily" ),
                            strstr( aircraft.variable[ dataVariable], "DAILY")))) {
            Name dailyName = "";
            memset( dailyName, 0, sizeof dailyName );
            snprintf( dailyName, sizeof dailyName / sizeof *dailyName,
                      "daily_%s",
                      aircraft.variable[ dataVariable ] );
            strncpy( aircraft.variable[ dataVariable ], dailyName,
                     sizeof dailyName / sizeof *dailyName );
          }
        }

        parameters->ok = writer( &aircraft, parameters );
      }
    } else {
      parameters->ok = writer( &aircraft, parameters );
    }
  }

  deallocateAircraft( &aircraft );
  POST0( isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: compareRegriddedAircraft - Read REGRIDDED-Aircraft input,
         compare it to CMAQ XDR data & write it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedAircraft( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok, parameters->input->ok( parameters->input ) );

  if ( ! AND2( parameters->compareFunction,
               parameters->data ) ) {
    failureMessage( "Invalid input for comparing." );
    parameters->ok = 0;
  } else {
    Aircraft aircraft;
    ZERO_OBJECT( &aircraft );
    parameters->ok = 0;

    if ( readRegriddedXDR( parameters->input, &aircraft ) ) {
      compareFunctionNameUnits( parameters->compareFunction,
                                parameters->convertFunction,
                                aircraft.variable[ 0 ],
                                aircraft.units[ 0 ],
                                parameters->variable,
                                parameters->units );

      if ( compareRegriddedXDR( parameters, &aircraft ) ) {
        Writer writer = dispatcher( parameters->format, 1 );
        CHECK( writer );

        if ( aircraft.totalRegriddedPoints == 0 ) {
          failureMessage( "No points projected onto the grid." );
        } else {
          parameters->ok = writer( &aircraft, parameters );
        }
      }
    }

    deallocateAircraft( &aircraft );
  }

  POST0( isValidParameters( parameters ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateAircraft - Deallocate contents of aircraft structure.
INPUTS:  Aircraft* aircraft Structure to deallocate contents of.
******************************************************************************/

static void deallocateAircraft( Aircraft* aircraft ) {
  PRE0( aircraft );
  FREE( aircraft->variable );
  FREE( aircraft->points );
  FREE( aircraft->notes );
  FREE( aircraft->data );
  FREE( aircraft->longitudes );
  FREE( aircraft->copyNotes );
  FREE( aircraft->regriddedNotes );
  ZERO_OBJECT( aircraft );
}



/******************************************************************************
PURPOSE: isValidAircraft - Check aircraft structure.
INPUTS:  const Aircraft* aircraft Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidAircraft( const Aircraft* aircraft ) {
  const Integer result =
    AND6( aircraft,
          aircraft->description[ 0 ],
          isValidUTCTimestamp( aircraft->firstTimestamp ),
          IMPLIES_ELSE( aircraft->variables > IMPLICIT_VARIABLES,
                        AND2( aircraft->variable[ IMPLICIT_VARIABLES ],
                              aircraft->units[ IMPLICIT_VARIABLES ] ),
                        AND2( aircraft->variable[ 0 ],
                              aircraft->units[ 0 ] ) ),
          IMPLIES( GT_ZERO2( aircraft->tracks, aircraft->totalPoints ),
            AND13( isValidUTCTimestamp( aircraft->lastTimestamp ),
                   aircraft->variables > IMPLICIT_VARIABLES,
                   isValidBounds( (const Real (*)[2]) aircraft->bounds ),
                   aircraft->variable,
                   aircraft->units,
                   aircraft->points,
                   minimumItemI( aircraft->points, aircraft->tracks ) > 0,
                   aircraft->notes,
                   aircraft->notes[ 0 ],
                   aircraft->notes[ aircraft->tracks - 1 ],
                   aircraft->data,
                   isNanFree( aircraft->data,
                              aircraft->variables * aircraft->totalPoints ),
                   aircraft->totalRegriddedPoints >= 0 ) ),
          IMPLIES( aircraft->totalRegriddedPoints > 0,
                   AND19( aircraft->timesteps > 0,
                          aircraft->outputPoints,
                          minimumItemI( aircraft->outputPoints,
                                        aircraft->timesteps ) >= 0,
                          aircraft->regriddedNotes,
                          aircraft->regriddedNotes[ 0 ],
                          aircraft->regriddedNotes[
                            aircraft->totalRegriddedPoints - 1 ],
                          aircraft->columns,
                          aircraft->rows,
                          aircraft->layers,
                          aircraft->gridLongitudes,
                          aircraft->gridLatitudes,
                          aircraft->gridElevations,
                          aircraft->gridData,
                          minimumItemI( aircraft->columns,
                                        aircraft->totalRegriddedPoints ) > 0,
                          minimumItemI( aircraft->rows,
                                        aircraft->totalRegriddedPoints ) > 0,
                          minimumItemI( aircraft->layers,
                                        aircraft->totalRegriddedPoints ) > 0,
                          isNanFree( aircraft->gridElevations,
                                     aircraft->totalRegriddedPoints ),
                          isNanFree( aircraft->gridData,
                                     aircraft->totalRegriddedPoints ),
                          validLongitudesAndLatitudes(
                            aircraft->totalRegriddedPoints,
                            aircraft->gridLongitudes,
                            aircraft->gridLatitudes ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isVectorVariable - Is the data variable a 2d wind vector?
INPUTS:  const Aircraft* aircraft Structure to check.
RETURNS: Integer 1 if vector, else 0.
******************************************************************************/

static Integer isVectorVariable( const Aircraft* aircraft ) {

  PRE05( aircraft, aircraft->variables > 0, aircraft->variable,
         aircraft->variable[ 0 ], aircraft->variable[ aircraft->variables ] );

    const Integer result =
      OR2( AND3( aircraft->variables == 2,
                 ! strcmp( aircraft->variable[ 0 ], "wind_u" ),
                 ! strcmp( aircraft->variable[ 1 ], "wind_v" ) ),
           AND3( aircraft->variables == IMPLICIT_VARIABLES + 2,
                 ! strcmp( aircraft->variable[ IMPLICIT_VARIABLES ],
                           "wind_u" ),
                 ! strcmp( aircraft->variable[ IMPLICIT_VARIABLES + 1 ],
                           "wind_v" ) ) );

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
PURPOSE: readXDR - Read input and initialize aircraft structure.
INPUTS:  Stream* input      Input stream read from and parse.
OUTPUTS: Aircraft* aircraft Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Input aircraft data in XDR format is:
Aircraft 2.0
http://mozaic.aero.obs-mip.fr/web/,MOZAICSubset
2006-07-03T00:00:00-0000 2006-07-04T23:59:59-0000
# Subset bounds: <min_lon> <min_lat> <max_lon> <max_lat>:
-180 -90 180 90
# Dimensions: variables points tracks:
5 371 2
# Variable names:
 timestamp longitude latitude elevation ozone
# Variable units:
 yyyymmddhhmmss deg deg m ppmV
# char notes[tracks][80] and
# IEEE-754 64-bit reals bounds[tracks][2=lon,lat][2=min,max] and
# MSB 64-bit integers points[tracks] and
# IEEE-754 64-bit reals data_1[points_1][variables]
  ... data_T[points_T][variables]:

******************************************************************************/

static Integer readXDR( Stream* input, Aircraft* aircraft ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         aircraft, aircraft->variable == 0, aircraft->data == 0 );

  Integer result = 0;

  input->readString( input, aircraft->description,
                     COUNT( aircraft->description ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( aircraft->description );

    if ( readTimestamps( input,
                         aircraft->firstTimestamp,
                         aircraft->lastTimestamp ) ) {

      if ( readDomain( input, aircraft->bounds ) ) {
        Integer dimensions[ 3 ] = { 0, 0, 0 };

        if ( readDimensions( input, COUNT( dimensions ), dimensions ) ) {
          aircraft->variables   = dimensions[ 0 ];
          aircraft->totalPoints = dimensions[ 1 ];
          aircraft->tracks      = dimensions[ 2 ];
          aircraft->variable    = NEW_ZERO( Name, aircraft->variables * 2 );

          if ( aircraft->variable ) {
            aircraft->units = aircraft->variable + aircraft->variables;

            if ( readVariablesAndUnits( input, aircraft->variables,
                                        aircraft->variable,
                                        aircraft->units ) ) {

              if ( skipInputLines( input, 4 ) ) {
                result = readXDRData( input, aircraft );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Aircraft data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidAircraft( aircraft ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read binary data from input.
INPUTS:  Stream* input      Input stream read from and parse.
OUTPUTS: Aircraft* aircraft Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDRData( Stream* input, Aircraft* aircraft ) {

  PRE010( input,
          input->ok( input ),
          input->isReadable( input ),
          aircraft,
          aircraft->variables > IMPLICIT_VARIABLES,
          aircraft->tracks > 0,
          aircraft->totalPoints > 0,
          aircraft->variables * aircraft->totalPoints > 0,
          aircraft->points == 0,
          aircraft->data == 0 );

  Integer result = 0;
  const Integer dataSize = aircraft->variables * aircraft->totalPoints;
  aircraft->data = NEW_ZERO( Real, dataSize );

  if ( aircraft->data ) {
    aircraft->points = NEW_ZERO( Integer, aircraft->tracks );

    if ( aircraft->points ) {
      aircraft->notes = NEW_ZERO( Note, aircraft->tracks );

      if ( aircraft->notes ) {
        readNotes( input, aircraft->tracks, aircraft->notes );

        if ( input->ok( input ) ) {
          const Integer boundsSize = aircraft->tracks * 2 * 2;
          CHECK( boundsSize < dataSize );

          input->read64BitReals( input, aircraft->data, boundsSize ); /*Skip*/

          if ( input->ok( input ) ) {
            input->read64BitIntegers( input, aircraft->points,
                                      aircraft->tracks );

            if ( input->ok( input ) ) {
              input->read64BitReals( input, aircraft->data, dataSize );

              if ( input->ok( input )  ) {
                result = isValidAircraft( aircraft );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Aircraft data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidAircraft( aircraft ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDR - Read REGRIDDED-Aircraft and initialize aircraft.
INPUTS:  Stream* input      Input stream read from and parse.
OUTPUTS: Aircraft* aircraft Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
NOTES:   Input data format is:

REGRIDDED-Aircraft 3.0
http://mozaic.aero.obs-mip.fr/web/,MOZAICSubset,XDRConvert
2006-07-03T00:00:00-0000
# timesteps
24
# Variable name:
ozone
# Variable units:
ppmV
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
 33 45 40 -97 6.37e+06 6.37e+06
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[25]:
 279 240 -1.008e+06 -1.62e+06 12000 12000 2 10000 1 0.995 0.99 0.98 0.97 \
 0.96 0.94 0.92 0.9 0.88 0.86 0.84 0.82 0.8 0.77 0.74 0.7 0.65 0.6 0.5 \
 0.4 0.3 0.2 0.1 0
# MSB 64-bit integers points[timesteps] and
# char notes[points][256] and
# IEEE-754 64-bit reals longitudes[points] and
# IEEE-754 64-bit reals latitudes[points] and
# IEEE-754 64-bit reals elevations[points] and
# MSB 64-bit integers columns[points] and
# MSB 64-bit integers rows[points] and
# MSB 64-bit integers layers[points] and
# IEEE-754 64-bit reals data[points]:

******************************************************************************/

static Integer readRegriddedXDR( Stream* input, Aircraft* aircraft ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         aircraft, aircraft->variable == 0, aircraft->data == 0,
         aircraft->gridData == 0 );

  Integer result = 0;
  input->readString( input, aircraft->description,
                     COUNT( aircraft->description ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( aircraft->description );

    if ( readTimestamp( input, aircraft->firstTimestamp ) ) {

      if ( readDimensions( input, 1, &aircraft->timesteps ) ) {
        aircraft->timestamps = NEW_ZERO( Integer, aircraft->timesteps );

        if ( aircraft->timestamps ) {
          Integer timestamp = fromUTCTimestamp( aircraft->firstTimestamp );
          Integer timestep = 0;

          for ( timestep = 0; timestep < aircraft->timesteps; ++timestep ) {
            aircraft->timestamps[ timestep ] = timestamp;
            incrementTimestamp( &timestamp );
          }

          if ( readVariablesAndUnits2( input, aircraft ) ) {
            char line[ 256 ] = "";
            int count = 9;
            memset( line, 0, sizeof line );
            input->readString( input, line, sizeof line / sizeof *line - 1 );

            if ( strcmp( line,
                         "# MSB 64-bit integers points[timesteps] and\n" )) {
              count += 4; /* Skip 4 line projection/grid. */
            }

            if ( skipInputLines( input, count - 1 ) ) {
              result = readRegriddedXDRData( input, aircraft );
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Aircraft data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidAircraft( aircraft ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readVariablesAndUnits2 - Read 1 (e.g., ozone) or 2 (wind_u wind_v)
         sets of variables and units.
INPUTS:  Stream* input       Input stream read from and parse.
OUTPUTS: Aircraft* aircraft  Allocated aircraft->variables,variable,units.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
NOTES:   Input data format is:
# Variable name:
wind_u wind_v
# Variable units:
m/s m/s
******************************************************************************/

static Integer readVariablesAndUnits2( Stream* input, Aircraft* aircraft ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         aircraft, aircraft->variables == 0, aircraft->variable == 0 );

  Integer result = 0;
  char line[ 256 ] = "";
  memset( line, 0, sizeof line );
  input->readString( input, line, sizeof line / sizeof *line - 1 );

  if ( ! strcmp( line, "# Variable name:\n" ) ) {
    input->readString( input, line, sizeof line / sizeof *line - 1 );
    aircraft->variables = wordsInString( line );

    if ( IN3( aircraft->variables, 1, 2 ) ) {
      aircraft->variable = NEW_ZERO( Name, aircraft->variables * 2 );

      if ( aircraft->variable ) {
        aircraft->units = aircraft->variable + aircraft->variables;

        if ( aircraft->variables == 1 ) {
          result = sscanf( line, "%s\n", aircraft->variable[ 0 ] ) == 1;
        } else {
          result = sscanf( line, "%s %s\n", aircraft->variable[ 0 ],
                           aircraft->variable[ 1 ] ) == 2;
        }

        if ( result ) {
          result = 0;
          input->readString( input, line, sizeof line / sizeof *line - 1 );

          if ( ! strcmp( line, "# Variable units:\n" ) ) {
            input->readString( input, line, sizeof line / sizeof *line - 1 );

            if ( wordsInString( line ) == aircraft->variables ) {

              if ( aircraft->variables == 1 ) {
                result = sscanf( line, "%s\n", aircraft->units[ 0 ] ) == 1;
              } else {
                result = sscanf( line, "%s %s\n", aircraft->units[ 0 ],
                                 aircraft->units[ 1 ] ) == 2;
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    failureMessage( "Invalid Aircraft header (variables/units)." );
    aircraft->variables = 0;
    FREE( aircraft->variable );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND7( IN3( aircraft->variables, 1, 2 ),
                         aircraft->variable,
                         aircraft->variable[ 0 ][ 0 ],
                         aircraft->variable[ aircraft->variables - 1 ][ 0 ],
                         aircraft->units,
                         aircraft->units[ 0 ][ 0 ],
                         aircraft->units[ aircraft->variables - 1 ][ 0 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRData - Read regridded binary array data from input.
INPUTS:  Stream* input       Input stream read from and parse.
OUTPUTS: Aircraft* aircraft  Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Input data format is:

# MSB 64-bit integers points[timesteps] and
# char notes[points][256] and
# IEEE-754 64-bit reals longitudes[points] and
# IEEE-754 64-bit reals latitudes[points] and
# IEEE-754 64-bit reals elevations[points] and
# MSB 64-bit integers columns[points] and
# MSB 64-bit integers rows[points] and
# MSB 64-bit integers layers[points] and
# IEEE-754 64-bit reals data[points]:

******************************************************************************/

static Integer readRegriddedXDRData( Stream* input, Aircraft* aircraft ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         aircraft, aircraft->timesteps > 0,
         IN3( aircraft->variables, 1, 2 ),
         aircraft->tracks == 0, aircraft->data == 0 );

  Integer result = 0;
  aircraft->outputPoints = NEW_ZERO( Integer, aircraft->timesteps );

  if ( aircraft->outputPoints ) {
    input->read64BitIntegers( input, aircraft->outputPoints,
                              aircraft->timesteps );

    if ( input->ok( input ) ) {
      const Integer count = aircraft->totalRegriddedPoints =
        sumI( aircraft->outputPoints, aircraft->timesteps );

      if ( count > 0 ) {
        aircraft->regriddedNotes = NEW_ZERO( RegriddedNote, count );

        if ( aircraft->regriddedNotes ) {
          const Integer isVector2 = isVectorVariable( aircraft );
          aircraft->gridLongitudes = NEW_ZERO( Real, count * ( 7 + isVector2 ));

          if ( aircraft->outputPoints ) {
            aircraft->gridLatitudes = aircraft->gridLongitudes + count;
            aircraft->gridElevations = aircraft->gridLatitudes + count;
            aircraft->columns = (Integer*) aircraft->gridElevations + count;
            aircraft->rows = aircraft->columns + count;
            aircraft->layers = aircraft->rows + count;
            aircraft->gridData = (Real*) aircraft->layers + count;

            readRegriddedNotes( input, count, aircraft->regriddedNotes );

            if ( input->ok( input ) ) {
              input->read64BitReals( input, aircraft->gridLongitudes, count*3);

              if ( input->ok( input ) ) {
                input->read64BitIntegers( input, aircraft->columns, count * 3);

                if ( input->ok( input ) ) {
                  const Integer count2 = isVector2 ? count + count : count;
                  input->read64BitReals( input, aircraft->gridData, count2 );

                  if ( input->ok( input ) ) {
                    result = isValidAircraft( aircraft );
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
    failureMessage( "Invalid Aircraft data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidAircraft( aircraft ) ) );
  return result;
}



/******************************************************************************
PURPOSE: compareRegriddedXDR - Compare Regridded data with CMAQ data.
INPUTS:  const Parameters* parameters  CMAQ data to compare to.
OUTPUTS: Aircraft* aircraft            Updated aircraft->data.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    Aircraft* aircraft ) {

  PRE05( parameters, isValidParameters( parameters ),
         parameters->compareFunction,
         aircraft, isValidAircraft( aircraft ) );

  Integer result = 0;

  if ( ! AND2( ! strcmp( parameters->timestamp, aircraft->firstTimestamp ),
               parameters->timesteps == aircraft->timesteps ) ) {
    failureMessage( "Mismatched time steps (%s %lld)"
                    " for comparison to CMAQ data (%s %lld).",
                    aircraft->firstTimestamp, aircraft->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else {
    Real* const aircraftData             = aircraft->gridData;
    const Integer* const aircraftLayers  = aircraft->layers;
    const Integer* const aircraftRows    = aircraft->rows;
    const Integer* const aircraftColumns = aircraft->columns;
    const Integer* const aircraftPoints  = aircraft->outputPoints;
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
    Integer aircraftIndex = 0;

    DEBUG( fprintf( stderr,
                    "grid L,R,C ranges: %lld..%lld %lld..%lld %lld..%lld\n",
                    firstLayer, lastLayer, firstColumn, lastColumn,
                    firstRow, lastRow ); )

    DEBUG( fprintf( stderr, "timesteps = %lld\n", timesteps ); )

    for ( timestep = 0; timestep < timesteps; ++timestep ) {
      const Integer points = aircraftPoints[ timestep ];
      const Integer timestepOffset = timestep * layersTimesRowsTimesColumns;
      Integer point = 0;

      DEBUG( fprintf( stderr, "timestep = %lld, points = %lld\n",
                      timestep, points ); )

      for ( point = 0; point < points; ++point, ++aircraftIndex ) {
        const Integer aircraftLayer  = aircraftLayers[  aircraftIndex ];
        const Integer aircraftRow    = aircraftRows[    aircraftIndex ];
        const Integer aircraftColumn = aircraftColumns[ aircraftIndex ];
        DEBUG( fprintf( stderr, "  @[ %lld, %lld, %lld ]: ",
                        aircraftLayer, aircraftRow, aircraftColumn ); )

        if ( AND3( IN_RANGE( aircraftLayer, firstLayer, lastLayer ),
                   IN_RANGE( aircraftRow, firstRow, lastRow ),
                   IN_RANGE( aircraftColumn, firstColumn, lastColumn ) ) ) {
          const Integer aircraftLayer0  = aircraftLayer  - firstLayer;
          const Integer aircraftRow0    = aircraftRow    - firstRow;
          const Integer aircraftColumn0 = aircraftColumn - firstColumn;
          const Integer dataIndex =
            timestepOffset + aircraftLayer0 * rowsTimesColumns +
            aircraftRow0 * columns + aircraftColumn0;
          CHECK4( IN_RANGE( aircraftLayer0, 0, layers - 1 ),
                  IN_RANGE( aircraftRow0, 0, rows - 1 ),
                  IN_RANGE( aircraftColumn0, 0, columns - 1 ),
                  IN_RANGE( dataIndex, 0,
                            timesteps * layers * rows * columns - 1 ) );
          const Real aircraftDatum = aircraftData[ aircraftIndex ];
          const Real cmaqDatum = cmaqData[ dataIndex ];
          const Real comparedDatum = comparer( aircraftDatum, cmaqDatum );
          aircraftData[ aircraftIndex ] = comparedDatum;
          result = 1;
          DEBUG( fprintf( stderr, "f(%lf, %lf) -> %lf\n",
                          aircraftDatum, cmaqDatum, comparedDatum ); )
        } else {
          aircraftData[ aircraftIndex ] = -9999.0;
          DEBUG( fprintf( stderr, "-9999\n" ); )
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "No points in output." );
  }

  POST02( IS_BOOL( result ), isValidAircraft( aircraft ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format output.
INPUTS:  Aircraft* aircraft  Structure to write.
         const Parameters*   Parameters input stream to read from.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeASCII( Aircraft* aircraft, const Parameters* parameters ) {

  PRE03( isValidAircraft( aircraft ),
         isValidParameters( parameters ),
         parameters->input->ok( parameters->input ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    writeASCIIHeader( aircraft, output );

    if ( output->ok( output ) ) {
      result = writeASCIIData( aircraft, output );
    }

    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCIIHeader - Write ASCII-format header line.
INPUTS:  const Aircraft* aircraft  Structure to write.
         Stream* output            Stream to write to.
******************************************************************************/

static void writeASCIIHeader( const Aircraft* aircraft, Stream* output ) {

  PRE03( isValidAircraft( aircraft ), output, output->isWritable( output ) );

  const char* const headerStart =
    "timestamp(UTC)\tlongitude(deg)\tlatitude(deg)\televation(m)";
  const char* const headerFormat = "\t%s(%s)";

  /* Write header row: */

  output->writeString( output, headerStart );

  if ( output->ok( output ) ) {
    const Integer variables = aircraft->variables;
    Integer variable = IMPLICIT_VARIABLES;

    do {
      output->writeString( output, headerFormat,
                           aircraft->variable[ variable ],
                           aircraft->units[ variable ] );

      if ( ! output->ok( output ) ) {
        variable = variables;
      }

      ++variable;
    } while ( variable < variables );

    if ( output->ok( output ) ) {
      output->writeString( output, "\tnote(-)\n" );
    }
  }
}



/******************************************************************************
PURPOSE: writeASCIIData - Write ASCII-format data lines.
INPUTS:  const Aircraft* aircraft  Structure to write.
OUTPUTS: Stream* output            Stream to write to.
RETURNS: Integer 1 if successful, else 0 and failureMesage is called.
******************************************************************************/

static Integer writeASCIIData( Aircraft* aircraft, Stream* output ) {

  PRE03( isValidAircraft( aircraft ), output, output->isWritable( output ) );

  Integer result = 0;
  const char* const dataFormat = "\t%28.6"REAL_F_FORMAT;
  const Integer variables = aircraft->variables;
  const Integer tracks    = aircraft->tracks;
  const Real* trackData   = aircraft->data;
  Integer track = 0;

  /* Write data rows: */

  for ( track = 0; track < tracks; ++track ) {
    const Integer trackPoints = aircraft->points[ track ];
    Integer point = 0;

    for ( point = 0; point < trackPoints; ++point, trackData += variables ) {
      Integer variable = 0;
      const Integer timestamp = trackData[ 0 ];
      UTCTimestamp timestampString;
      CHECK( isValidYYYYMMDDHHMMSS( timestamp ) );
      toUTCTimestamp2( timestamp, timestampString );
      output->writeString( output, timestampString ); /* Begin row. */

      for ( variable = 1;
            AND2( output->ok( output ), variable < variables );
            ++variable ) {
        const Real datum = trackData[ variable ];
        output->writeString( output, dataFormat, datum );
      }

      if ( output->ok( output ) ) { /* End row: */
        char format[ 80 ] = "";
        memset( format, 0, sizeof format );
        sprintf( format, "\t%%-%lus\n", sizeof aircraft->notes[ 0 ] - 1 );
        output->writeString( output, format, aircraft->notes[ track ] );
      }

      if ( ! output->ok( output ) ) {
        point = trackPoints;
        track = tracks;
      }
    }
  }

  result = output->ok( output );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDS - Write COARDS-format data.
INPUTS:  Aircraft* aircraft            Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDS( Aircraft* aircraft, const Parameters* parameters) {

  PRE02( isValidAircraft( aircraft ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
      aircraft->variables * aircraft->totalPoints * 4 + /*variables(points).*/
      aircraft->totalPoints * 3 * 4 + 2000;
      /* yyyyddd(points),hhmmss(points),time(points) + header/extra. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file = createNetCDFFile(parameters->netcdfFileName, create64BitFile);

  if ( file != -1 ) {

    if ( writeCOARDSHeader( file, aircraft ) ) {
      result = writeCOARDSData( file, aircraft );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSHeader - Write header to file.
INPUTS:  Integer file              NetCDF file to write to.
         const Aircraft* aircraft  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSHeader( Integer file, const Aircraft* aircraft ) {

  PRE02( file != -1, isValidAircraft( aircraft ) );

  Integer result = 0;
  enum { DIMENSIONS = 3 };
  const char* const names[ DIMENSIONS ] = { "points", "tracks", "length" };
  Integer sizes[ DIMENSIONS ] = { 0, 0, 0 };
  Integer dimensionIds[ 3 ] = { -1, -1, -1 };
  sizes[ 0 ] = aircraft->totalPoints;
  sizes[ 1 ] = aircraft->tracks;
  sizes[ 2 ] = sizeof aircraft->notes[ 0 ] / sizeof (char);

  if ( createDimensions( file, 3, names, sizes, dimensionIds ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "notes", "-", NC_CHAR, 0, 2,
                           dimensionIds + 1 ) != -1 ) {

        if ( createLongitudeAndLatitude( file, 1, dimensionIds ) ) {
          const Integer variables = aircraft->variables;
          Integer index = AIRCRAFT_ELEVATION; /* Only write elevation and data.*/

          do {
            const char* const units =
              strcmp( aircraft->units[ index ], "-"   ) == 0 ? "none" :
              strcmp( aircraft->units[ index ], "deg" ) == 0 ? "degrees" :
              aircraft->units[ index ];
            const Integer ok =
              createVariable( file, aircraft->variable[ index ], units,
                              NC_FLOAT, 1, 1, dimensionIds ) != -1;

            if ( ! ok ) {
              index = variables;
            }

            ++index;
          } while ( index < variables );

          if ( index == variables ) {

            if ( writeExtraAttributes( file,
                                       (const Real (*)[2]) aircraft->bounds,
                                       dimensionIds[ 0 ] ) ) {
                UTCTimestamp timestamp = "";
                Line history = "";
                appendToLine( history, aircraft->description );
                appendToLine( history, ",XDRConvert" );
                toUTCTimestamp2( aircraft->data[ 0 ], timestamp );
                result =
                  writeStandardContents( file, history, timestamp,
                                         dimensionIds[ 0 ],
                                         aircraft->totalPoints, 0 );
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
INPUTS:  Integer file        NetCDF file to write to.
         Aircraft* aircraft  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSData( Integer file, Aircraft* aircraft ) {

  PRE02( file != -1, isValidAircraft( aircraft ) );

  Integer result = 0;
  const Integer totalPoints = aircraft->totalPoints;
  const Integer variables   = aircraft->variables;
  const Integer timesteps   = aircraft->totalPoints;

  /* Generate date-time data in aircraft->timestamps[]: */

  CHECK2( aircraft->copyData == 0, aircraft->timestamps == 0 );
  aircraft->copyData = NEW_ZERO( Real, totalPoints * variables );
  aircraft->timestamps =
    aircraft->copyData ? NEW_ZERO( Integer, 3 * timesteps) :0;

  if ( aircraft->timestamps ) {
    const Integer yyyydddhhmmStart =
      fromUTCTimestamp( aircraft->firstTimestamp );
    Integer* yyyyddd = aircraft->timestamps;
    Integer* hhmmss  = yyyyddd + timesteps;
    Real* time       = (Real*) hhmmss + timesteps;
    Integer variable = 0;
    Real* input = aircraft->data;
    Real* output = aircraft->copyData;
    Integer point = 0;

    /* Copy data to timestamps and non-interleaved copyData: */

    for ( point = 0; point < totalPoints; ++point, ++output ) {
      const Integer yyyymmddhhmmss = input[ AIRCRAFT_TIMESTAMP ];
      const Integer yyyymmdd0 = yyyymmddhhmmss / 1000000;
      const Integer hhmmss0   = yyyymmddhhmmss % 1000000;
      const Integer yyyyddd0  = convertYearMonthDay( yyyymmdd0 );
      const Integer yyyydddhhmmss00 = yyyyddd0 * 1000000 + hhmmss0;
      const Integer yyyydddhhmmNow = yyyydddhhmmss00 / 100;
      const Real fractionalTime =
        fractionalHours( yyyydddhhmmStart, yyyydddhhmmNow );
      yyyyddd[ point ] = yyyyddd0;
      hhmmss[  point ] = hhmmss0;
      time[    point ] = fractionalTime;

      for ( variable = 0; variable < variables; ++variable ) {
        output[ variable * totalPoints ] = *input++;
      }
    }

    /* Write and ruin copyData: */

    for ( variable = AIRCRAFT_LONGITUDE; variable < variables; ++variable ) {
      const char* const variableName = aircraft->variable[ variable ];

      if ( ! writeAllData( file, variableName, totalPoints, 1, 1, 1,
                           aircraft->copyData + variable * totalPoints ) ) {
        variable = variables;
      }
    }

    /* Write and ruin yyyyddd: */

    result = variable == variables;

    result = AND2( result,
                   writeAllIntData( file, "yyyyddd", totalPoints, 1, 1, 1,
                                yyyyddd ) );

    result = AND2( result,
                   writeAllIntData( file, "hhmmss", totalPoints, 1, 1, 1,
                                    hhmmss ) );

    result = AND2( result,
                   writeAllData( file, "time", totalPoints, 1, 1, 1, time ) );

    if ( result ) {
      const Integer bufferLength =
        aircraft->tracks * sizeof aircraft->notes[ 0 ];
      char* buffer = NEW_ZERO( char, bufferLength + 1 );
      result = buffer != 0;

      if ( result ) {
        expandNotes( aircraft->tracks, aircraft->notes, buffer );
        result = writeAllCharData( file, "notes", aircraft->tracks,
                                   sizeof aircraft->notes[ 0 ], buffer );
        FREE( buffer );
      }
    }
  }

  FREE( aircraft->copyData );
  FREE( aircraft->timestamps );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  Aircraft* aircraft  Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.

NOTES:   Output data format is:

REGRIDDED-Aircraft 3.0
http://mozaic.aero.obs-mip.fr/web/,MOZAICSubset,XDRConvert
2006-07-03T00:00:00-0000
# timesteps
24
# Variable name:
ozone
# Variable units:
ppmV
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
 33 45 40 -97 6.37e+06 6.37e+06
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[25]:
 279 240 -1.008e+06 -1.62e+06 12000 12000 2 10000 1 0.995 0.99 0.98 0.97 \
 0.96 0.94 0.92 0.9 0.88 0.86 0.84 0.82 0.8 0.77 0.74 0.7 0.65 0.6 0.5 \
 0.4 0.3 0.2 0.1 0
# MSB 64-bit integers points[timesteps] and
# char notes[points][256] and
# IEEE-754 64-bit reals longitudes[points] and
# IEEE-754 64-bit reals latitudes[points] and
# IEEE-754 64-bit reals elevations[points] and
# MSB 64-bit integers columns[points] and
# MSB 64-bit integers rows[points] and
# MSB 64-bit integers layers[points] and
# IEEE-754 64-bit reals data[points]:

******************************************************************************/

static Integer writeRegriddedXDR( Aircraft* aircraft,
                                  const Parameters* parameters ) {

  PRE02( isValidAircraft( aircraft ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = aircraft->timesteps;
    const Integer points    = aircraft->totalRegriddedPoints;
    const Integer isVector2 = isVectorVariable( aircraft );
    const Integer variableIndex =
      aircraft->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
    Name variable = "";
    aggregateName( aircraft->variable[ variableIndex ], hoursPerTimestep,
                   variable );

    output->writeString( output,
                         "REGRIDDED-Aircraft 3.0\n"
                         "%s,XDRConvert\n"
                         "%s\n"
                         "# timesteps\n%"INTEGER_FORMAT"\n",
                         aircraft->description,
                         aircraft->firstTimestamp,
                         timesteps );

    if ( ! isVector2 ) {
      output->writeString( output,
                           "# Variable name:\n%s\n"
                           "# Variable units:\n%s\n",
                           variable,
                           aircraft->units[ variableIndex ] );
    } else {
      Name variable2 = "";
      aggregateName( aircraft->variable[ variableIndex + 1 ], hoursPerTimestep,
                     variable2 );
      output->writeString( output,
                           "# Variable name:\n%s %s\n"
                           "# Variable units:\n%s %s\n",
                           variable, variable2,
                           aircraft->units[ variableIndex ],
                           aircraft->units[ variableIndex + 1 ] );
    }

    if ( output->ok( output ) ) {
      writeProjectionAndGrid( parameters->grid, output );

      if ( output->ok( output ) ) {
        output->writeString( output,
                             "# MSB 64-bit integers points[timesteps] and\n"
                             "# char notes[points][256] and\n"
                             "# IEEE-754 64-bit reals longitudes[points] and\n"
                             "# IEEE-754 64-bit reals latitudes[points] and\n"
                             "# IEEE-754 64-bit reals elevations[points] and\n"
                             "# MSB 64-bit integers columns[points] and\n"
                             "# MSB 64-bit integers rows[points] and\n"
                             "# MSB 64-bit integers layers[points] and\n"
                             "# IEEE-754 64-bit reals data[points]:\n" );

        if ( output->ok( output ) ) {
          output->write64BitIntegers( output, aircraft->outputPoints, timesteps );

          if ( output->ok( output ) ) {
            writeRegriddedNotes( output, points,
                            (const RegriddedNote*) aircraft->regriddedNotes );

            if ( output->ok( output ) ) {
              output->write64BitReals(output, aircraft->gridLongitudes, points);

              if ( output->ok( output ) ) {
                output->write64BitReals(output, aircraft->gridLatitudes, points);

                if ( output->ok( output ) ) {
                  output->write64BitReals( output, aircraft->gridElevations,
                                           points );

                  if ( output->ok( output ) ) {
                    output->write64BitIntegers( output, aircraft->columns,
                                                points );

                    if ( output->ok( output ) ) {
                      output->write64BitIntegers( output, aircraft->rows,
                                                  points );

                      if ( output->ok( output ) ) {
                        output->write64BitIntegers( output, aircraft->layers,
                                                    points );

                        if ( output->ok( output ) ) {
                          const Integer points2 =
                            isVector2 ? points + points : points;
                          output->write64BitReals( output, aircraft->gridData,
                                                   points2 );
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
    }

    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedASCII - Write regridded ASCII-format data.
INPUTS:  Aircraft* aircraft            Structure to write.
         const Parameters* parameters  Structure of options.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedASCII( Aircraft* aircraft,
                                    const Parameters* parameters ) {

  PRE03( isValidAircraft( aircraft ), aircraft->variables > 0,
        isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer isVector2 = isVectorVariable( aircraft );
    const char* const headerStart =
      "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tELEVATION(m)"
      "\tCOLUMN(-)\tROW(-)\tLAYER(-)";

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {
      const Integer variableIndex =
        aircraft->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
      const Integer hoursPerTimestep =
        parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
      Name variable = "";
      aggregateName( aircraft->variable[ variableIndex ], hoursPerTimestep,
                     variable );

      if ( isVector2 ) {
        const char* const headerFormat = "\t%s(%s)\t%s(%s)\tnote(-)\n";
        Name variable2 = "";
        aggregateName( aircraft->variable[ variableIndex + 1], hoursPerTimestep,
                       variable2 );

        output->writeString( output, headerFormat,
                             variable,
                             aircraft->units[ variableIndex ],
                             variable2,
                             aircraft->units[ variableIndex + 1 ] );

      } else {
        const char* const headerFormat = "\t%s(%s)\tnote(-)\n";
        output->writeString( output, headerFormat,
                             variable,
                             aircraft->units[ variableIndex ] );
      }

      if ( output->ok( output ) ) {
        const Integer timesteps = aircraft->timesteps;
        const Real* longitudes  = aircraft->gridLongitudes;
        const Real* latitudes   = aircraft->gridLatitudes;
        const Real* elevations  = aircraft->gridElevations;
        const Integer* columns  = aircraft->columns;
        const Integer* rows     = aircraft->rows;
        const Integer* layers   = aircraft->layers;
        const Real* data        = aircraft->gridData;
        const Real* data2 =
          isVector2 ? aircraft->gridData + aircraft->totalRegriddedPoints : 0;
        const RegriddedNote* regriddedNotes =
          (const RegriddedNote*) aircraft->regriddedNotes;
        const Integer hoursPerTimestep =
          parameters->aggregationTimesteps ? parameters->aggregationTimesteps
          : 1;
        Integer timestep = 0;
        Integer yyyydddhh00 =
          ( fromUTCTimestamp( aircraft->firstTimestamp ) / 100 ) * 100;
        UTCTimestamp timestamp = "";

        /* Write data rows: */

        do {
          const Integer points = aircraft->outputPoints[ timestep ];
          Integer point = 0;
          toUTCTimestamp( yyyydddhh00, timestamp );

          for ( point = 0; point < points; ++point ) {
            const Real longitude = *longitudes++;
            const Real latitude  = *latitudes++;
            const Real elevation = *elevations++;
            const Integer column = *columns++;
            const Integer row    = *rows++;
            const Integer layer  = *layers++;
            const Real value     = *data++;

            if ( isVector2 ) {
              const char* const dataFormat =
                "%s\t%10.4lf\t%10.4lf\t%10.4lf"
                "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT
                "\t%10.4lf\t%10.4lf";
              const Real value2 = *data2++;
              output->writeString( output, dataFormat,
                                   timestamp, longitude, latitude, elevation,
                                   column, row, layer, value, value2 );
            } else {
              const char* const dataFormat =
                "%s\t%10.4lf\t%10.4lf\t%10.4lf"
                "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT
                "\t%10.4lf";
              output->writeString( output, dataFormat,
                                   timestamp, longitude, latitude, elevation,
                                   column, row, layer, value );
            }

            if ( output->ok( output ) ) { /* End row: */
              char format[ 80 ] = "";
              memset( format, 0, sizeof format );
              sprintf( format, "\t%%-%lus\n", sizeof *regriddedNotes - 1 );
              output->writeString( output, format, *regriddedNotes );
              ++regriddedNotes;
            }

            if ( ! output->ok( output ) ) {
              point = points;
              timestep = timesteps;
            }
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
INPUTS:  Aircraft* aircraft            Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDS( Aircraft* aircraft,
                                     const Parameters* parameters ) {

  PRE02( isValidAircraft( aircraft ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    aircraft->totalRegriddedPoints * ( 7 + 1 ) * 4 + 10000;
    /* lon, lat, elv, col, row, lay, time, hdr. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedCOARDSHeader( file, hoursPerTimestep, aircraft ) ) {
      result = writeRegriddedCOARDSData( file, aircraft, parameters );
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
         const Aircraft* aircraft  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Aircraft* aircraft ) {

  PRE03( file != -1, hoursPerTimestep > 0, isValidAircraft( aircraft ) );

  Integer result = 0;
  const char* const names[ 2 ] = { "points", "length" };
  Integer sizes[ 2 ] = { 0, 0 };
  Integer dimensionIds[ 2 ] = { -1, -1 };
  sizes[ 0 ] = aircraft->totalRegriddedPoints;
  sizes[ 1 ] = sizeof aircraft->regriddedNotes[ 0 ] / sizeof (char);

  if ( createDimensions( file, 2, names, sizes, dimensionIds ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "notes", "-", NC_CHAR, 0, 2, dimensionIds )
           != -1 ) {

        if ( createVariable( file, "column", "-",
                             NC_INT, 0, 1, dimensionIds ) != -1 ) {

          if ( createVariable( file, "row", "-",
                               NC_INT, 0, 1, dimensionIds ) != -1 ) {

            if ( createVariable( file, "layer", "-",
                                 NC_INT, 0, 1, dimensionIds ) != -1 ) {

              if ( createLongitudeAndLatitude( file, 1, dimensionIds ) ) {

                if ( createVariable( file, "elevation", "-",
                                     NC_FLOAT, 0, 1, dimensionIds ) != -1) {
                  const Integer variableIndex =
                    aircraft->variables > IMPLICIT_VARIABLES ?
                      IMPLICIT_VARIABLES
                    : 0;
                  const Integer isVector2 = isVectorVariable( aircraft );
                  Name variable = "";
                  aggregateName( aircraft->variable[ variableIndex ],
                                 hoursPerTimestep, variable );

                  if ( createVariable( file,
                                       variable,
                                       aircraft->units[ variableIndex ],
                                       NC_FLOAT, 1, 1, dimensionIds ) != -1 ) {
                    result = 1;

                    if ( isVector2 ) {
                      Name variable2 = "";
                      aggregateName( aircraft->variable[ variableIndex + 1 ],
                                     hoursPerTimestep,
                                     variable2 );
                      result =
                        createVariable( file, variable2,
                                        aircraft->units[ variableIndex + 1 ],
                                        NC_FLOAT, 1, 1, dimensionIds ) != -1;
                    }

                    if ( result ) {
                      UTCTimestamp timestamp = "";
                      Line history = "";
                      appendToLine( history, aircraft->description );
                      appendToLine( history, ",XDRConvert" );
                      toUTCTimestamp( aircraft->timestamps[ 0 ], timestamp );

                      result =
                        writeStandardContents( file, history, timestamp,
                                               dimensionIds[ 0 ], 0, 0 );
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

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDSData - Write COARDS-format data to file.
INPUTS:  Integer file              NetCDF file to write to.
         const Aircraft* aircraft  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSData( Integer file,
                                         const Aircraft* aircraft,
                                         const Parameters* parameters ) {

  PRE03( file != -1, isValidAircraft( aircraft ),
         isValidParameters( parameters ) );

  Integer result = 0;
  const Integer count = aircraft->totalRegriddedPoints;
  const Integer isVector2 = isVectorVariable( aircraft );

  if ( writeAllIntData( file, "column", count, 1, 1, 1,
                        aircraft->columns ) ) {

    if ( writeAllIntData( file, "row", count, 1, 1, 1,
                          aircraft->rows ) ) {

      if ( writeAllIntData( file, "layer", count, 1, 1, 1,
                            aircraft->layers ) ) {

        if ( writeAllData( file, "longitude", count, 1, 1, 1,
                            aircraft->gridLongitudes ) ) {

          if ( writeAllData( file, "latitude", count, 1, 1, 1,
                              aircraft->gridLatitudes ) ) {

            if ( writeAllData( file, "elevation", count, 1, 1, 1,
                                aircraft->gridElevations ) ) {
              const Integer variableIndex =
                aircraft->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES
                : 0;
              const Integer hoursPerTimestep =
                parameters->aggregationTimesteps ?
                  parameters->aggregationTimesteps
                : 1;
              Name variable = "";
              aggregateName( aircraft->variable[ variableIndex ],
                             hoursPerTimestep, variable );

              if ( writeAllData( file, variable, count, 1, 1, 1,
                                 aircraft->gridData ) ) {
                result = 1;

                if ( isVector2 ) {
                  Name variable2 = "";
                  aggregateName( aircraft->variable[ variableIndex + 1 ],
                                 hoursPerTimestep, variable2 );
                  result =
                    writeAllData( file, variable2, count, 1, 1, 1,
                                  aircraft->gridData + count );
                }

                if ( result ) {
                  const Integer hoursPerTimestep =
                    parameters->aggregationTimesteps ?
                      parameters->aggregationTimesteps
                    : 1;

                  timeData( aircraft->timesteps, hoursPerTimestep, count,
                            aircraft->outputPoints, aircraft->gridData );

                  if ( writeAllData( file, "time", count, 1, 1, 1,
                                     aircraft->gridData ) ) {
                    const Integer bufferLength =
                      count * sizeof aircraft->regriddedNotes[ 0 ];
                    char* buffer = NEW_ZERO( char, bufferLength + 1 );

                    if ( buffer ) {
                      expandRegriddedNotes( count, aircraft->regriddedNotes,
                                            buffer );
                      result = writeAllCharData( file, "notes", count,
                                          sizeof aircraft->regriddedNotes[ 0 ],
                                                buffer );
                      FREE( buffer );
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

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPI - Write regridded IOAPI-format data.
INPUTS:  Aircraft* aircraft                  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPI( Aircraft* aircraft,
                                    const Parameters* parameters ) {

  PRE02( isValidAircraft( aircraft ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    aircraft->totalRegriddedPoints * 5 * 4 + 10000; /* lon,lat,elv,var, hdr */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedIOAPIHeader( file, hoursPerTimestep,
                                    aircraft, parameters->grid ) ) {
      result = writeRegriddedIOAPIData( file, hoursPerTimestep, aircraft,
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
INPUTS:  Integer file              NetCDF file to write to.
         Integer hoursPerTimestep  Hours per timestep: 1, 24, etc.
         const Aircraft* aircraft  Structure to write.
         const Grid* grid  Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Aircraft* aircraft,
                                          const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0,
         isValidAircraft( aircraft ),
         grid, grid->invariant( grid ) );

  Integer result = 0;
  enum { VARIABLES = 4 }; /* LONGITUDE, LATITUDE, ELEVATION, aircraft. */
  Name variableNames[ VARIABLES + 1 ] = {
    "LONGITUDE", "LATITUDE", "ELEVATION", "aircraft", "aircraft2"
  };
  Name variableUnits[ VARIABLES + 1 ] = { "deg", "deg", "m", "-", "m/s" };
  const Integer layers = grid->layers( grid );
  const Integer firstTimestamp = fromUTCTimestamp( aircraft->firstTimestamp );
  const Integer variableIndex =
    aircraft->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
  const Integer isVector2 = isVectorVariable( aircraft );
  Line history = "";
  appendToLine( history, aircraft->description );
  appendToLine( history, ",XDRConvert" );
  aggregateName( aircraft->variable[ variableIndex ], hoursPerTimestep,
                 variableNames[ VARIABLES - 1 ] );
  variableNames[ VARIABLES - 1 ][ 15 ] = '\0';
  strncpy( variableUnits[ VARIABLES - 1 ],
           aircraft->units[ variableIndex ], 16 );
  variableUnits[ VARIABLES - 1 ][ 16 ] = '\0';
  uppercase( variableNames[ VARIABLES - 1 ] );
  lowercase( variableUnits[ VARIABLES - 1 ] );

  if ( isVector2 ) {
    aggregateName( aircraft->variable[ variableIndex + 1 ], hoursPerTimestep,
                   variableNames[ VARIABLES ] );
    variableNames[ VARIABLES ][ 15 ] = '\0';
    strncpy( variableUnits[ VARIABLES ],
             aircraft->units[ variableIndex + 1 ], 16 );
    variableUnits[ VARIABLES ][ 16 ] = '\0';
    uppercase( variableNames[ VARIABLES ] );
    lowercase( variableUnits[ VARIABLES ] );
  }

  result = writeM3IOHeader( file, aircraft->timesteps, hoursPerTimestep,
                            firstTimestamp,
                            VARIABLES + isVector2, layers,
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
         const Aircraft* aircraft  Structure to write.
         const Grid* grid  Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Aircraft* aircraft,
                                        const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0, isValidAircraft( aircraft ),
         grid, grid->invariant( grid ) );

  Integer result = 0;
  const Integer layers   = grid->layers( grid );
  const Integer rows     = grid->rows( grid );
  const Integer columns  = grid->columns( grid );
  const Integer cells    = layers * rows * columns;
  Real* expandedGridData = NEW_ZERO( Real, cells );

  if ( expandedGridData ) {
    const Integer timesteps = aircraft->timesteps;
    const Real scale = 1.0;
    Integer timestep = 0;
    Integer offset = 0;

    if ( writeM3IOGrid( grid, timesteps, layers, file ) ) {

      do {
        const Integer points = aircraft->outputPoints[ timestep ];
        Name variable = "";
        memset( variable, 0, sizeof variable );
        strncpy( variable, "ELEVATION", 16 );
        uppercase( variable );

        copyDataToGrid3( points,
                         aircraft->layers         + offset,
                         aircraft->rows           + offset,
                         aircraft->columns        + offset,
                         aircraft->gridElevations + offset,
                         scale, layers, rows, columns,
                         expandedGridData );

        if ( ! writeM3IOData( file, variable,
                              timestep, layers, rows, columns,
                              expandedGridData ) ) {
          timestep = timesteps;
        } else {
          const Integer variableIndex =
            aircraft->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
          aggregateName( aircraft->variable[ variableIndex ], hoursPerTimestep,
                         variable );
          variable[ 15 ] = '\0';
          uppercase( variable );

          copyDataToGrid3( points,
                           aircraft->layers   + offset,
                           aircraft->rows     + offset,
                           aircraft->columns  + offset,
                           aircraft->gridData + offset,
                           scale, layers, rows, columns,
                           expandedGridData );

          if ( ! writeM3IOData( file, variable,
                                timestep, layers, rows, columns,
                                expandedGridData ) ) {
            timestep = timesteps;
          } else {
            const Integer isVector2 = isVectorVariable( aircraft );

            if ( isVector2 ) {
              const Real* gridData2 =
                aircraft->gridData + aircraft->totalRegriddedPoints;
              Name variable2 = "";
              aggregateName( aircraft->variable[ variableIndex + 1 ],
                             hoursPerTimestep, variable2 );
              variable2[ 15 ] = '\0';
              uppercase( variable2 );

              copyDataToGrid3( points,
                               aircraft->layers   + offset,
                               aircraft->rows     + offset,
                               aircraft->columns  + offset,
                               gridData2 + offset,
                               scale, layers, rows, columns,
                               expandedGridData );

              if ( ! writeM3IOData( file, variable2,
                                    timestep, layers, rows, columns,
                                    expandedGridData ) ) {
                timestep = timesteps;
              }
            }
          }
        }

        offset += points;
        ++timestep;
      } while ( timestep < timesteps );

      result = timestep == timesteps;
    }

    FREE( expandedGridData );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: regridAircraft - Regrid data.
INPUTS:  Integer method      E.g., AGGREGATE_MEAN.
         Grid* grid          Grid to project and aggregate points into.
         Aircraft* aircraft  Data to regrid.
OUTPUTS: Aircraft* aircraft  Regridded data.
******************************************************************************/

static void regridAircraft( Integer method, Grid* grid, Aircraft* aircraft ) {

  PRE08( IS_VALID_AGGREGATE_METHOD( method ),
         grid,
         grid->invariant( grid ),
         isValidAircraft( aircraft ),
         aircraft->totalRegriddedPoints == 0,
         aircraft->longitudes == 0,
         aircraft->copyNotes == 0,
         aircraft->regriddedNotes == 0 );

  const Integer variables = aircraft->variables;

  if ( IN3( variables, IMPLICIT_VARIABLES + 1, IMPLICIT_VARIABLES + 2 ) ) {
    const Integer timesteps =
      hoursInRange( aircraft->firstTimestamp, aircraft->lastTimestamp );
    const Integer isVector2 = isVectorVariable( aircraft );
    const Integer inputVariables = 4 + isVector2; /* lon, lat, elv, dat, dat2.*/
    const Integer inputSize      = aircraft->totalPoints; /*All points in 1hr*/
    const Integer inputDataSize  = inputVariables * inputSize;
    const Integer outputVariables =
      7 + isVector2; /* glon,glat,gelv,col,row,lay,gdata,gdata2. */
    const Integer outputSize      = inputSize; /* 1 point in each hr? */
    const Integer outputDataSize  = outputVariables * outputSize;
    const Integer dataSize  = inputDataSize + outputDataSize + timesteps * 2;

    aircraft->copyNotes = NEW_ZERO( Note, aircraft->totalPoints );
    aircraft->regriddedNotes =
      aircraft->copyNotes ? NEW_ZERO( RegriddedNote, outputSize ) : 0;
    aircraft->longitudes =
      aircraft->regriddedNotes ? NEW_ZERO( Real, dataSize ) : 0;

    if ( aircraft->longitudes ) {
      Integer totalRegriddedPoints = 0;
      Integer timestep = 0;
      Integer yyyydddhh00 =
        ( fromUTCTimestamp( aircraft->firstTimestamp ) / 100 ) * 100;
      aircraft->latitudes      = aircraft->longitudes     + inputSize;
      aircraft->elevations     = aircraft->latitudes      + inputSize;
      aircraft->copyData       = aircraft->elevations     + inputSize;
      aircraft->gridLongitudes =
        aircraft->copyData + inputSize + inputSize * isVector2;
      aircraft->gridLatitudes  = aircraft->gridLongitudes + outputSize;
      aircraft->gridElevations = aircraft->gridLatitudes  + outputSize;
      aircraft->gridData       = aircraft->gridElevations + outputSize;
      aircraft->columns        = (Integer*)
        aircraft->gridData + outputSize + outputSize * isVector2;
      aircraft->rows           = aircraft->columns        + outputSize;
      aircraft->layers         = aircraft->rows           + outputSize;
      aircraft->outputPoints   = aircraft->layers         + outputSize;
      aircraft->timestamps     = aircraft->outputPoints   + timesteps;
      aircraft->timesteps      = timesteps;

      do {
        const Real* const copyData2 =
          isVector2 ? aircraft->copyData + inputSize : 0;
        Real* const gridData2 =
          isVector2 ? aircraft->gridData + outputSize : 0;
        const Integer inputPoints =
          copyDataForTimestamp( yyyydddhh00, aircraft );
        aircraft->timestamps[ timestep ] = yyyydddhh00;

        if ( inputPoints ) {
          const Real minimumValidValue = -900.0;
          Integer outputPoints = 0;

          DEBUG( fprintf( stderr, "regridding %lld points for %lld...\n",
                          inputPoints, yyyydddhh00 ); )
          CHECK( totalRegriddedPoints < outputSize );


          grid->regrid( grid, method, minimumValidValue, inputPoints, 1,
                        aircraft->longitudes, aircraft->latitudes,
                        aircraft->elevations, aircraft->copyData,
                        isVector2 ? copyData2 : 0,
                        (const Note*) aircraft->copyNotes,
                        &outputPoints,
                        aircraft->columns        + totalRegriddedPoints,
                        aircraft->rows           + totalRegriddedPoints,
                        aircraft->layers         + totalRegriddedPoints,
                        aircraft->gridLongitudes + totalRegriddedPoints,
                        aircraft->gridLatitudes  + totalRegriddedPoints,
                        aircraft->gridElevations + totalRegriddedPoints,
                        aircraft->gridData       + totalRegriddedPoints,
                        isVector2 ? gridData2 + totalRegriddedPoints : 0,
                        aircraft->regriddedNotes + totalRegriddedPoints );

          aircraft->outputPoints[ timestep ] = outputPoints;
          totalRegriddedPoints += outputPoints;

          DEBUG( fprintf( stderr, "outputPoints = %lld,"
                          " totalRegriddedPoints = %lld\n",
                          outputPoints, totalRegriddedPoints ); )
        }

        incrementTimestamp( &yyyydddhh00 );
        ++timestep;
      } while ( timestep < timesteps );

      aircraft->totalRegriddedPoints = totalRegriddedPoints;

      if ( isVector2 ) { /* Append gridData2 onto end of gridData: */
        const Real* const gridData2 = aircraft->gridData + outputSize;
        memcpy( aircraft->gridData + totalRegriddedPoints,
               gridData2,
               totalRegriddedPoints * sizeof (Real) );
      }
    }
  }

  POST02( aircraft->totalRegriddedPoints >= 0,
          IMPLIES( aircraft->totalRegriddedPoints > 0,
            AND10( IN_RANGE( minimumItemI( aircraft->outputPoints,
                                           aircraft->timesteps ),
                             0, aircraft->totalRegriddedPoints ),
                   IN_RANGE( maximumItemI( aircraft->outputPoints,
                                           aircraft->timesteps ),
                                           1, aircraft->totalRegriddedPoints ),
                   IN_RANGE( minimumItemI( aircraft->columns,
                                           aircraft->totalRegriddedPoints),
                                           1, grid->columns( grid ) ),
                   IN_RANGE( maximumItemI( aircraft->columns,
                                           aircraft->totalRegriddedPoints),
                                           1, grid->columns( grid ) ),
                   IN_RANGE( minimumItemI( aircraft->rows,
                                           aircraft->totalRegriddedPoints),
                                           1, grid->rows( grid ) ),
                   IN_RANGE( maximumItemI( aircraft->rows,
                                           aircraft->totalRegriddedPoints),
                                           1, grid->rows( grid ) ),
                   IN_RANGE( minimumItemI( aircraft->layers,
                                           aircraft->totalRegriddedPoints),
                                           1, grid->layers( grid ) ),
                   IN_RANGE( maximumItemI( aircraft->layers,
                                           aircraft->totalRegriddedPoints),
                                           1, grid->layers( grid ) ),
                   validLongitudesAndLatitudes( aircraft->totalRegriddedPoints,
                                                aircraft->gridLongitudes,
                                                aircraft->gridLatitudes ),
                   isNanFree( aircraft->gridData,
                              aircraft->totalRegriddedPoints ) ) ) );
}



/******************************************************************************
PURPOSE: copyDataForTimestamp - Copy data for given regrid timestamp.
INPUTS:  Integer yyyydddhh00  Timestamp to copy data for.
         Aircraft* aircraft   aircraft->data data to copy.
OUTPUTS: Aircraft* aircraft   aircraft->longitudes, latitudes, elevations,
                              copyData, copyNotes.
RETURNS: Integer number of points copied for the timestamp.
******************************************************************************/

static Integer copyDataForTimestamp(Integer yyyydddhh00, Aircraft* aircraft) {

  PRE03( isValidTimestamp( yyyydddhh00 ),
         isValidAircraft( aircraft ),
         aircraft->copyNotes );

  const Integer isVector2 = isVectorVariable( aircraft );
  const Integer variables = aircraft->variables;
  const Integer points = aircraft->totalPoints;
  const Real* data = aircraft->data;
  Real* longitudes = aircraft->longitudes;
  Real* latitudes  = aircraft->latitudes;
  Real* elevations = aircraft->elevations;
  Real* copyData   = aircraft->copyData;
  Real* copyData2  = isVector2 ? aircraft->copyData + points : 0;
  Note* notes      = aircraft->notes;
  Note* copyNotes  = aircraft->copyNotes;
  Integer point = 0;
  Integer result = 0;

  for ( point = 0; point < points; ++point, data += variables ) {
    const Integer pointTimestamp = data[ AIRCRAFT_TIMESTAMP ];
    Integer timestamp = 0;
    UTCTimestamp timestampString = "";
    toUTCTimestamp2( pointTimestamp, timestampString );
    timestamp = fromUTCTimestamp( timestampString );
    timestamp /= 100;
    timestamp *= 100;

    if ( timestamp == yyyydddhh00 ) {
      const Real longitude = data[ AIRCRAFT_LONGITUDE ];
      const Real latitude  = data[ AIRCRAFT_LATITUDE ];
      const Real elevation = data[ AIRCRAFT_ELEVATION ];
      const Real datum     = data[ IMPLICIT_VARIABLES ];
      const Integer track =
        binIndex( point, aircraft->tracks, aircraft->points );
      CHECK5( IN_RANGE( track, 0, aircraft->tracks - 1 ),
              isValidLongitude( longitude ), isValidLatitude( latitude ),
              ! isNan( elevation ), ! isNan( datum ) );
      *longitudes++ = longitude;
      *latitudes++  = latitude;
      *elevations++ = elevation;
      *copyData++   = datum;

      if ( isVector2 ) {
        const Real datum2 = data[ IMPLICIT_VARIABLES + 1 ];
        *copyData2++ = datum2;
      }

      strcpy( *copyNotes, notes[ track ] );
      ++copyNotes;
      ++result;
    }
  }

  POST02( result >= 0,
         IMPLIES( result > 0,
                  AND4( validLongitudesAndLatitudes( result,
                                                     aircraft->longitudes,
                                                     aircraft->latitudes ),
                        isNanFree( aircraft->elevations, result ),
                        isNanFree( aircraft->copyData, result ),
                        aircraft->copyNotes[ 0 ][ 0 ] ) ) );
  return result;
}



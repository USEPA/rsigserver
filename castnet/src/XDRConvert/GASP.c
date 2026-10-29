
/******************************************************************************
PURPOSE: GASP.c - Define routines for processing GASP data.

NOTES:

HISTORY: 2008/02 plessel.todd@epa.gov, Created.
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
  Integer      variables;       /* E.g., 3 = longitude, latitude, aod. */
  Integer      timesteps;       /* E.g., 24. */
  Integer      scans;           /* E.g., 35 half-hour daylight scans. */
  Integer      totalPoints;     /* Sum of points[ scan ]. */
  Integer      maximumPoints;   /* Number of points in largest scan. */
  Real         domain[ 2 ][ 2]; /*domain[LONGITUDE LATITUDE][MINIMUM MAXIMUM]*/
  UTCTimestamp timestamp;
  Line         note;            /* File note/description. */
  Name*        variable;        /* variable[ variables ]. E.g., "aod"*/
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
  Real*        gridData;        /* gridData[ totalRegriddedPoints ]. */
} GASP;

typedef Integer (*Writer)( GASP* gasp, const Parameters* parameters );

typedef struct {
  Integer format;         /* FORMAT_XDR, etc. */
  Writer writer;          /* Routine that writes data in this format. */
  Writer regriddedWriter; /* Routine that writes regridded data in format. */
} Entry;

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateGASP( GASP* gasp );

static Integer isValidGASP( const GASP* gasp );

static Writer dispatcher( Integer format, Integer regrid );

static Integer readXDR( Stream* input, GASP* gasp );

static Integer readXDRData( Stream* input, GASP* gasp );

static Integer readRegriddedXDR( Stream* input, GASP* gasp );

static Integer readRegriddedXDRData( Stream* input, GASP* gasp );

static Integer compareRegriddedXDR( const Parameters* parameters, GASP* gasp );

static void countGASPPoints( GASP* gasp );

static Integer writeASCII( GASP* gasp, const Parameters* parameters );

static void writeASCIIHeader( const GASP* gasp, Stream* output );

static Integer writeASCIIData( GASP* gasp, Stream* input, Stream* output );

static Integer writeCOARDS( GASP* gasp, const Parameters* parameters );

static Integer writeCOARDSHeader( Integer file, const GASP* gasp );

static Integer writeCOARDSData( Integer file, Stream* input,
                                const GASP* gasp );


static Integer writeRegriddedXDR( GASP* gasp, const Parameters* unused );

static Integer writeRegriddedASCII( GASP* gasp, const Parameters* unused );

static Integer writeRegriddedCOARDS( GASP* gasp,
                                     const Parameters* parameters );

static Integer writeRegriddedCOARDSHeader( Integer file, const GASP* gasp );

static Integer writeRegriddedCOARDSData( Integer file, const GASP* gasp );

static Integer writeRegriddedIOAPI(GASP* gasp, const Parameters* parameters);

static Integer writeRegriddedIOAPIHeader( Integer file, const GASP* gasp,
                                          const Grid* grid );

static Integer writeRegriddedIOAPIData( Integer file, const GASP* gasp,
                                        const Grid* grid );

static Integer writeRegriddedMCMC( GASP* gasp, const Parameters* unused );


static void regridGASP( Stream* input, Integer method, Grid* grid, GASP* gasp);

static Integer readScanDataForTimestamp( Integer yyyydddhh00,
                                         Stream* input,
                                         GASP* gasp,
                                         Integer* points );

static Integer readScanData( Stream* input, Integer variables, Integer points,
                             Real longitudes[], Real latitudes[],
                             Real data[] );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateGASP - Read input and write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void translateGASP( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok,
         parameters->input->ok( parameters->input ) );

  GASP gasp;
  ZERO_OBJECT( &gasp );
  parameters->ok = 0;

  if ( readXDR( parameters->input, &gasp ) ) {
    Writer writer = dispatcher( parameters->format, parameters->regrid );

    if ( ! writer ) {
      failureMessage( "Invalid/unsupported format/regrid specification." );
    } else if ( parameters->regrid ) {
      regridGASP( parameters->input, parameters->regrid, parameters->grid,
                  &gasp );

      if ( gasp.totalRegriddedPoints == 0 ) {
        failureMessage( "No points projected onto the grid." );
      } else {
        parameters->ok = writer( &gasp, parameters );
      }
    } else {
      parameters->ok = writer( &gasp, parameters );
    }
  }

  deallocateGASP( &gasp );
  POST0( isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: compareRegriddedGASP - Read REGRIDDED-GASP input, compare it to
         CMAQ XDR data and write it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedGASP( Parameters* parameters ) {
    
  PRE03( isValidParameters( parameters ),
         parameters->ok, parameters->input->ok( parameters->input ) );

  if ( ! AND3( ! parameters->regrid,
               OR2( parameters->compareFunction, parameters->convertFunction ),
               parameters->data ) ) {
    failureMessage( "Invalid input for comparing." );
    parameters->ok = 0;
  } else {
    GASP gasp;
    ZERO_OBJECT( &gasp );
    parameters->ok = 0;

    if ( readRegriddedXDR( parameters->input, &gasp ) ) {
      compareFunctionNameUnits( parameters->compareFunction,
                                parameters->convertFunction,
                                gasp.variable[ 0 ], gasp.units[ 0 ], 
                                parameters->variable, 
                                parameters->units );

      if ( compareRegriddedXDR( parameters, &gasp ) ) {
        Writer writer = dispatcher( parameters->format, 1 );
        CHECK( writer );

        if ( gasp.totalRegriddedPoints == 0 ) {
          failureMessage( "No points projected onto the grid." );
        } else {
          parameters->ok = writer( &gasp, parameters );
        }
      }
    }

    deallocateGASP( &gasp );
  }

  POST0( isValidParameters( parameters ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateGASP - Deallocate contents of gasp structure.
INPUTS:  GASP* gasp Structure to deallocate contents of.
******************************************************************************/

static void deallocateGASP( GASP* gasp ) {
  PRE0( gasp );

  if ( gasp->data == 0 ) { /* Read Regridded XDR compare case. */
    FREE( gasp->outputPoints );
    FREE( gasp->gridLongitudes );    
  }

  FREE( gasp->variable );
  FREE( gasp->timestamps );
  FREE( gasp->data );
  ZERO_OBJECT( gasp );
}



/******************************************************************************
PURPOSE: isValidGASP - Check gasp structure.
INPUTS:  const GASP* gasp Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidGASP( const GASP* gasp ) {
  const Integer result =
    AND14( gasp,
           gasp->note[ 0 ],
           isValidUTCTimestamp( gasp->timestamp ),
           GT_ZERO2( gasp->variables, gasp->timesteps ),
           isValidLongitudeLatitude( gasp->domain[ LONGITUDE ][ MINIMUM ],
                                     gasp->domain[ LATITUDE  ][ MINIMUM ] ),
           isValidLongitudeLatitude( gasp->domain[ LONGITUDE ][ MAXIMUM ],
                                     gasp->domain[ LATITUDE  ][ MAXIMUM ] ),
           gasp->domain[ LONGITUDE ][ MINIMUM ] <=
             gasp->domain[ LONGITUDE ][ MAXIMUM ],
           gasp->domain[ LATITUDE  ][ MINIMUM ] <=
             gasp->domain[ LATITUDE  ][ MAXIMUM ],
           gasp->variable,
           gasp->units,
           gasp->variable[ 0 ],
           gasp->units[ 0 ],
           IMPLIES( GT_ZERO2( gasp->scans, gasp->totalPoints ),
             AND12( gasp->maximumPoints ==
                      maximumItemI( gasp->points, gasp->scans ),
                    gasp->timestamps,
                    isValidTimestamp( gasp->timestamps[ 0 ] ),
                    isValidTimestamp( gasp->timestamps[ gasp->scans - 1 ] ),
                    gasp->timestamps[ gasp->scans - 1 ] >= gasp->timestamps[0],
                    gasp->points,
                    gasp->points[ 0 ] > 0,
                    gasp->points[ gasp->scans - 1 ] > 0,
                    gasp->data,
                    isNanFree( gasp->data,
                               gasp->maximumPoints ),
                    IMPLIES( gasp->longitudes,
                             AND2( gasp->latitudes,
                                   validLongitudesAndLatitudes(
                                     gasp->totalPoints,
                                     gasp->longitudes, gasp->latitudes ) ) ),
                    gasp->totalRegriddedPoints >= 0 ) ),
          IMPLIES( gasp->totalRegriddedPoints > 0,
                   AND11( gasp->outputPoints,
                          minimumItemI( gasp->outputPoints, gasp->timesteps )
                            >= 0,
                          gasp->columns,
                          gasp->rows,
                          gasp->gridLongitudes,
                          gasp->gridLatitudes,
                          gasp->gridData,
                          minimumItemI( gasp->columns,
                                        gasp->totalRegriddedPoints ) > 0,
                          minimumItemI( gasp->rows,
                                        gasp->totalRegriddedPoints ) > 0,
                          isNanFree( gasp->gridData,
                                     gasp->totalRegriddedPoints ),
                          validLongitudesAndLatitudes(
                            gasp->totalRegriddedPoints,
                            gasp->gridLongitudes, gasp->gridLatitudes ))));
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
    { FORMAT_MCMC,   0,           writeRegriddedMCMC   },
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
PURPOSE: readXDR - Read input and initialize gasp structure.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: GASP* gasp Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDR( Stream* input, GASP* gasp ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         gasp, gasp->variable == 0, gasp->data == 0 );

  Integer result = 0;

  input->readString( input, gasp->note, COUNT( gasp->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( gasp->note );

    if ( readTimestamp( input, gasp->timestamp ) ) {
      Integer dimensions[ 3 ] = { 0, 0, 0 };

      if ( readDimensions( input, COUNT( dimensions ), dimensions ) ) {
        gasp->variables = dimensions[ 0 ];
        gasp->timesteps = dimensions[ 1 ];
        gasp->scans     = dimensions[ 2 ];
        gasp->variable  = NEW_ZERO( Name, gasp->variables * 2 );

        if ( gasp->variable ) {
          gasp->units = gasp->variable + gasp->variables;

          if ( readVariablesAndUnits( input, gasp->variables,
                                      gasp->variable, gasp->units ) ) {

            if ( readDomain( input, gasp->domain ) ) {

              if ( skipInputLines( input, 3 ) ) {
                result = readXDRData( input, gasp );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid GASP data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidGASP( gasp ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read initial binary data from input.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: GASP* gasp Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDRData( Stream* input, GASP* gasp ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         gasp, gasp->variables > 0, gasp->scans > 0,
         gasp->timestamps == 0, gasp->data == 0 );

  Integer result = 0;
  gasp->timestamps = NEW_ZERO( Integer, gasp->scans * 2 );

  if ( gasp->timestamps ) {
    gasp->points = gasp->timestamps + gasp->scans;
    input->read64BitIntegers( input, gasp->timestamps, gasp->scans );

    if ( AND3( input->ok( input ),
               isValidTimestamp( gasp->timestamps[ 0 ] ),
               isValidTimestamp( gasp->timestamps[ gasp->scans - 1 ] ) ) ) {

      input->read64BitIntegers( input, gasp->points, gasp->scans );

      /* Sum rows x columns per scan: */

      if ( input->ok( input ) ) {
        countGASPPoints( gasp );

        if ( gasp->totalPoints > 0 ) {
          gasp->data = NEW_ZERO( Real, gasp->maximumPoints );
          result = isValidGASP( gasp );
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid GASP data." );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, isValidGASP( gasp ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDR - Read REGRIDDED-GASP and initialize gasp.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: GASP* gasp Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readRegriddedXDR( Stream* input, GASP* gasp ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         gasp, gasp->variable == 0, gasp->data == 0, gasp->gridData == 0 );

  Integer result = 0;

  input->readString( input, gasp->note, COUNT( gasp->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( gasp->note );

    if ( readTimestamp( input, gasp->timestamp ) ) {

      if ( readDimensions( input, 1, &gasp->timesteps ) ) {
        gasp->timestamps = NEW_ZERO( Integer, gasp->timesteps );

        if ( gasp->timestamps ) {
          Integer timestamp = fromUTCTimestamp( gasp->timestamp );
          Integer timestep = 0;

          for ( timestep = 0; timestep < gasp->timesteps; ++timestep ) {
            gasp->timestamps[ timestep ] = timestamp;
            incrementTimestamp( &timestamp );
          }

          gasp->variables = 1;
          gasp->variable = NEW_ZERO( Name, 2 );

          if ( gasp->variable ) {
            gasp->units = gasp->variable + 1;

            if ( readVariablesAndUnits( input, gasp->variables,
                                        gasp->variable, gasp->units ) ) {
              char line[ 256 ] = "";
              int count = 6;
              memset( line, 0, sizeof line );
              input->readString( input, line, sizeof line / sizeof *line - 1 );

              if ( strcmp( line,
                           "# MSB 64-bit integers points[timesteps] and\n" )) {
                count += 4; /* Skip 4 line projection/grid. */
              }

              if ( skipInputLines( input, count - 1 ) ) {
                result = readRegriddedXDRData( input, gasp );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid GASP data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidGASP( gasp ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRData - Read regridded data from input.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: GASP* gasp Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readRegriddedXDRData( Stream* input, GASP* gasp ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         gasp, gasp->timesteps > 0, gasp->variables > 0,
         gasp->scans == 0, gasp->data == 0 );

  Integer result = 0;

  gasp->outputPoints = NEW_ZERO( Integer, gasp->timesteps );

  if ( gasp->outputPoints ) {
    input->read64BitIntegers( input, gasp->outputPoints, gasp->timesteps );

    if ( input->ok( input ) ) {
      const Integer count = gasp->totalRegriddedPoints =
        sumI( gasp->outputPoints, gasp->timesteps );

      if ( count > 0 ) {
        gasp->gridLongitudes = NEW_ZERO( Real, count * 5 );

        if ( gasp->outputPoints ) {
          gasp->gridLatitudes = gasp->gridLongitudes + count;
          gasp->columns       = (Integer*) gasp->gridLatitudes + count;
          gasp->rows          = gasp->columns + count;
          gasp->gridData      = (Real*) gasp->rows + count;

          input->read64BitReals( input, gasp->gridLongitudes, count * 2 );

          if ( input->ok( input ) ) {
            input->read64BitIntegers( input, gasp->columns, count * 2 );

            if ( input->ok( input ) ) {
              input->read64BitReals( input, gasp->gridData, count );

              if ( input->ok( input ) ) {
                result = isValidGASP( gasp );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid GASP data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidGASP( gasp ) ) );
  return result;
}



/******************************************************************************
PURPOSE: compareRegriddedXDR - Compare Regridded data with CMAQ data.
INPUTS:  const Parameters* parameters  CMAQ data to compare to.
OUTPUTS: GASP* gasp                  Updated gasp->data.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    GASP* gasp ) {

  PRE05( parameters, isValidParameters( parameters ),
         OR2( parameters->compareFunction, parameters->convertFunction ),
         gasp, isValidGASP( gasp ) );

  Integer result = 0;

  if ( ! IMPLIES( parameters->compareFunction,
                  AND2( ! strcmp( parameters->timestamp, gasp->timestamp ),
                        parameters->timesteps == gasp->timesteps ) ) ) {
    failureMessage( "Mismatched comparison timesteps: "
                    "satellite (%s %lld) vs model (%s %lld).",
                    gasp->timestamp, gasp->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else if ( ! IMPLIES( parameters->convertFunction,
                    gasp->timesteps == parameters->timesteps * 24 ) ) {
    failureMessage( "Mismatched conversion timesteps: "
                    "satellite (%s %lld) vs model (%s %lld).",
                    gasp->timestamp, gasp->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else {
    const Integer daily = gasp->timesteps == parameters->timesteps * 24;
    const Integer timesteps          = gasp->timesteps;
    Real* const gaspData             = gasp->gridData;
    const Integer* const gaspRows    = gasp->rows;
    const Integer* const gaspColumns = gasp->columns;
    const Integer* const gaspPoints  = gasp->outputPoints;
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
    Integer cmaqTimestep = 0;
    Integer timestep = 0;
    Integer gaspIndex = 0;

    DEBUG( fprintf( stderr, "timesteps = %lld\n", timesteps ); )

    for ( timestep = 0; timestep < timesteps; ++timestep ) {
      const Integer points = gaspPoints[ timestep ];
      const Integer cmaqTimestep = daily ? timestep / 24 : timestep;
      const Integer timestepOffset = cmaqTimestep * rowsTimesColumns;
      Integer point = 0;

      DEBUG( fprintf( stderr, "cmaqTimestep = %lld, "
                      "timestep = %lld, points = %lld\n",
                      cmaqTimestep, timestep, points ); )

      for ( point = 0; point < points; ++point, ++gaspIndex ) {
        const Integer gaspRow    = gaspRows[    gaspIndex ];
        const Integer gaspColumn = gaspColumns[ gaspIndex ];
        DEBUG( fprintf( stderr, "  @[ %lld, %lld ]: ", gaspRow, gaspColumn);)

        if ( AND2( IN_RANGE( gaspRow, firstRow, lastRow ),
                   IN_RANGE( gaspColumn, firstColumn, lastColumn ) ) ) {
          const Integer gaspRow0    = gaspRow    - firstRow;
          const Integer gaspColumn0 = gaspColumn - firstColumn;
          const Integer cmaqIndex =
            timestepOffset + gaspRow0 * columns + gaspColumn0;
          CHECK3( IN_RANGE( gaspRow0, 0, rows - 1 ),
                  IN_RANGE( gaspColumn0, 0, columns - 1 ),
                  IN_RANGE( cmaqIndex, 0,
                            parameters->timesteps * rows * columns - 1 ) );
          const Real gaspDatum = gaspData[ gaspIndex ];
          const Real cmaqDatum = cmaqData[ cmaqIndex ];
          const Real newDatum =
            comparer ? comparer( gaspDatum, cmaqDatum )
            : converter( gaspDatum, cmaqDatum, cmaqData2[ cmaqIndex ] );
          gaspData[ gaspIndex ] = newDatum;
          result = 1;
          DEBUG( if ( comparer ) fprintf( stderr, "f(%lf, %lf) -> %lf\n",
                                          gaspDatum, cmaqDatum, newDatum );
                else fprintf( stderr, "f(%lf, %lf, %lf) -> %lf\n",
                              gaspDatum, cmaqDatum, cmaqData2[ cmaqIndex ],
                              newDatum ); )
        } else {
          gaspData[ gaspIndex ] = -9999.0;
          DEBUG( fprintf( stderr, "-9999\n" ); )
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "No points in output." );
  }

  POST02( IS_BOOL( result ), isValidGASP( gasp ) );
  return result;
}



/******************************************************************************
PURPOSE: countGASPPoints - Compute sum and maximum of product of scan
         rows x columns.
INPUTS:  GASP* gasp      Partially initialized structure.
OUTPUTS: GASP* gasp      points = Sum of each scan rows x columns.
                         maximumPoints = Maximum scan rows x columns.
******************************************************************************/

static void countGASPPoints( GASP* gasp ) {

  PRE05( gasp, gasp->scans > 0, gasp->points, gasp->totalPoints == 0,
         gasp->maximumPoints == 0 );

  Integer scan = 0;

  do {
    const Integer scanPoints = gasp->points[ scan ];
    gasp->totalPoints += scanPoints;
    gasp->maximumPoints = MAX( gasp->maximumPoints, scanPoints );
    ++scan;
  } while ( scan < gasp->scans );

  POST02( gasp->totalPoints >= 0,
          IMPLIES_ELSE( gasp->totalPoints > 0,
                        gasp->maximumPoints > 0,
                        gasp->maximumPoints == 0 ) );
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format output.
INPUTS:  GASP* gasp       Structure to write.
         const Parameters*  Parameters input stream to read from.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeASCII( GASP* gasp, const Parameters* parameters ) {

  PRE03( isValidGASP( gasp ),
         isValidParameters( parameters ),
         parameters->input->ok( parameters->input ) );

  Integer result = 0;

  /*
   * First reallocate gasp->data to be large enough to hold
   * all variables of the largest scan:
   */

  FREE( gasp->data );
  gasp->data = NEW_ZERO( Real, gasp->maximumPoints * gasp->variables );

  if ( gasp->data ) {
    Stream* output = newFileStream( "-stdout", "wb" );

    if ( output ) {
      writeASCIIHeader( gasp, output );

      if ( output->ok( output ) ) {
        result = writeASCIIData( gasp, parameters->input, output );
      }

      FREE_OBJECT( output );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCIIHeader - Write ASCII-format header line.
INPUTS:  const GASP* gasp  Structure to write.
         Stream* output      Stream to write to.
******************************************************************************/

static void writeASCIIHeader( const GASP* gasp, Stream* output ) {

  PRE03( isValidGASP( gasp ), output, output->isWritable( output ) );

  const char* const headerStart =
    "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)";
  const char* const headerFormat = "\t%s(%s)";

  /* Write header row: */

  output->writeString( output, headerStart );

  if ( output->ok( output ) ) {
    const Integer variables = gasp->variables;
    Integer variable = 2;

    do {
      output->writeString( output, headerFormat,
                           gasp->variable[ variable ],
                           gasp->units[ variable ] );

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
INPUTS:  const GASP* gasp  Structure to write.
         Stream* input       Stream to read from.
         Stream* output      Stream to write to.
RETURNS: Integer 1 if successful, else 0 and failureMesage is called.
******************************************************************************/

static Integer writeASCIIData( GASP* gasp, Stream* input, Stream* output ) {

  PRE05( isValidGASP( gasp ), input, input->isReadable( input ),
         output, output->isWritable( output ) );

  Integer result = 0;
  const Integer variables = gasp->variables;
  const char* const dataFormat = "\t%28.6"REAL_F_FORMAT;
  const Integer scans = gasp->scans;
  Integer scan = 0;

  /* Write data rows: */

  do {
    const Integer scanPoints = gasp->points[ scan ];
    const Integer scanSize = variables * scanPoints;

    input->read64BitReals( input, gasp->data, scanSize );

    if ( ! input->ok( input ) ) {
      scan = scans;
    } else {
      Real* scanData = gasp->data;
      Integer point = 0;
      const Integer timestamp = gasp->timestamps[ scan ];
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
INPUTS:  GASP* gasp                  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDS( GASP* gasp, const Parameters* parameters ) {

  PRE02( isValidGASP( gasp ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
      gasp->variables * gasp->totalPoints * 4 + /* variables(points). */
      gasp->totalPoints * 3 * 4 + 2000;
      /* yyyyddd(points),hhmmss(points),time(points) + header/extra. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file = createNetCDFFile(parameters->netcdfFileName, create64BitFile);

  if ( file != -1 ) {

    if ( writeCOARDSHeader( file, gasp ) ) {
      result = writeCOARDSData( file, parameters->input, gasp );
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
         const GASP* gasp  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSHeader( Integer file, const GASP* gasp ) {

  PRE02( file != -1, isValidGASP( gasp ) );

  Integer result = 0;
  const char* const name = "points";
  Integer dimensionId = -1;

  if ( createDimensions( file, 1, &name, &gasp->totalPoints, &dimensionId ) ) {

    if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {
      const Integer variables = gasp->variables;
      Integer index = 2;

      do {
        const char* const units =
          strcmp( gasp->units[ index ], "-"   ) == 0 ? "none" :
          strcmp( gasp->units[ index ], "deg" ) == 0 ? "degrees" :
          gasp->units[ index ];
        const Integer ok =
          createVariable( file, gasp->variable[ index ], units,
                          NC_FLOAT, 1, 1, &dimensionId ) != -1;

         if ( ! ok ) {
           index = variables;
         }

         ++index;
      } while ( index < variables );

      if ( index == variables ) {

        if ( writeExtraAttributes( file, (const Real (*)[2]) gasp->domain,
                                   dimensionId ) ) {
          UTCTimestamp timestamp = "";
          Line history = "";
          appendToLine( history, gasp->note );
          appendToLine( history, ",XDRConvert" );
          toUTCTimestamp( gasp->timestamps[ 0 ], timestamp );
          result = writeStandardContents( file, history,
                                          timestamp,
                                          dimensionId, gasp->totalPoints, 0);
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
         const GASP* gasp  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSData(Integer file, Stream* input, const GASP* gasp) {

  PRE05( input, input->ok( input ), input->isReadable( input ),
         file != -1, isValidGASP( gasp ) );

  Integer result = 0;
  const Integer scans     = gasp->scans;
  const Integer variables = gasp->variables;
  Integer scan = 0;
  Integer variable = 0;
  Integer start = 0;

  do {
    const Integer count = gasp->points[ scan ];
    variable = 0;

    do {
      const char* const variableName =
        variable == 0 ? "longitude" :
        variable == 1 ? "latitude" :
        gasp->variable[ variable ];

      input->read64BitReals( input, gasp->data, count );

      if ( ! AND2( input->ok( input ),
                   writeSomeData( file, variableName, start,
                                  count, 1, 1, 1, gasp->data ) ) ) {
        scan = scans;
        variable = variables;
      }

      ++variable;
    } while ( variable < variables );

    start += count;
    ++scan;
  } while ( scan < scans );

  result = AND3( scan == scans, variable == variables,
                 writeTimeData( file, gasp->scans, 1, 0,
                                gasp->timestamps, gasp->points, gasp->data ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  GASP* gasp       Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedXDR( GASP* gasp, const Parameters* parameters ) {

  PRE02( isValidGASP( gasp ),isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = gasp->timesteps;
    const Integer points    = gasp->totalRegriddedPoints;
    const Integer index     = gasp->variables >= 3 ? 2 : 0;

    output->writeString( output,
                         "REGRIDDED-GASP 2.0\n"
                         "%s,XDRConvert\n"
                         "%s\n"
                         "# timesteps\n%"INTEGER_FORMAT"\n"
                         "# Variable name:\n%s\n"
                         "# Variable units:\n%s\n",
                         gasp->note,
                         gasp->timestamp,
                         timesteps,
                         gasp->variable[ index ],
                         gasp->units[ index ] );
  
    writeProjectionAndGrid( parameters->grid, output );

    output->writeString( output,
                 "# MSB 64-bit integers points[timesteps] and\n"
                 "# IEEE-754 64-bit reals longitudes[timesteps][points] and\n"
                 "# IEEE-754 64-bit reals latitudes[timesteps][points] and\n"
                 "# MSB 64-bit integers columns[timesteps][points] and\n"
                 "# MSB 64-bit integers rows[timesteps][points] and\n"
                        "# IEEE-754 64-bit reals data[timesteps][points]:\n" );
                         

    if ( output->ok( output ) ) {
      output->write64BitIntegers( output, gasp->outputPoints, timesteps );

      if ( output->ok( output ) ) {
        output->write64BitReals( output, gasp->gridLongitudes, points );

        if ( output->ok( output ) ) {
          output->write64BitReals( output, gasp->gridLatitudes, points );

          if ( output->ok( output ) ) {
            output->write64BitIntegers( output, gasp->columns, points);

            if ( output->ok( output ) ) {
              output->write64BitIntegers( output, gasp->rows, points );

              if ( output->ok( output ) ) {
                output->write64BitReals( output, gasp->gridData, points );
                result = output->ok( output );
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
INPUTS:  GASP* gasp       Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedASCII( GASP* gasp, const Parameters* unused ) {

  PRE02( isValidGASP( gasp ), gasp->variables > 0 );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const char* const headerStart =
      "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)"
      "\tCOLUMN(-)\tROW(-)";
    const char* const headerFormat = "\t%s(%s)\n";
    const char* const dataFormat =
      "%s\t%10.4lf\t%10.4lf"
      "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%10.4lf\n";

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {
      const Integer index = gasp->variables >= 3 ? 2 : 0;
      output->writeString( output, headerFormat,
                           gasp->variable[ index ], gasp->units[ index ] );

      if ( output->ok( output ) ) {
        const Integer timesteps = gasp->timesteps;
        const Real* longitudes  = gasp->gridLongitudes;
        const Real* latitudes   = gasp->gridLatitudes;
        const Integer* columns  = gasp->columns;
        const Integer* rows     = gasp->rows;
        const Real* data        = gasp->gridData;

        Integer timestep = 0;
        Integer yyyydddhh00 =
          ( fromUTCTimestamp( gasp->timestamp ) / 100 ) * 100;
        UTCTimestamp timestamp;

        /* Write data rows: */

        do {
          const Integer points = gasp->outputPoints[ timestep ];
          Integer point = 0;
          toUTCTimestamp( yyyydddhh00, timestamp );

          for ( point = 0; point < points; ++point ) {
            const Real longitude = *longitudes++;
            const Real latitude  = *latitudes++;
            const Integer column = *columns++;
            const Integer row    = *rows++;
            const Real value     = *data++;

            output->writeString( output, dataFormat,
                                 timestamp, longitude, latitude,
                                 column, row, value );

            if ( ! output->ok( output ) ) {
              point = points;
              timestep = timesteps;
            }
          }

          incrementTimestamp( &yyyydddhh00 );
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
INPUTS:  GASP* gasp                  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDS( GASP* gasp,
                                     const Parameters* parameters ) {

  PRE02( isValidGASP( gasp ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    gasp->totalRegriddedPoints * 5 * 4 + 10000; /*lon,lat,col,row,time,hdr.*/
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {

    if ( writeRegriddedCOARDSHeader( file, gasp ) ) {
      result = writeRegriddedCOARDSData( file, gasp );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedCOARDSHeader - Write header to file.
INPUTS:  Integer file      NetCDF file to write to.
         const GASP* gasp  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           const GASP* gasp ) {

  PRE02( file != -1, isValidGASP( gasp ) );

  Integer result = 0;
  const char* const dimensionName = "points";
  Integer dimensionId = -1;
  const Integer dimension = gasp->totalRegriddedPoints;

  if ( createDimensions( file, 1, &dimensionName, &dimension, &dimensionId)) {

    if ( createVariable( file, "column", "-",
                         NC_INT, 0, 1, &dimensionId ) != -1 ) {

      if ( createVariable( file, "row", "-",
                           NC_INT, 0, 1, &dimensionId ) != -1 ) {

        if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {
          const Integer index = gasp->variables >= 3 ? 2 : 0;

          if ( createVariable( file,
                               gasp->variable[ index ], gasp->units[ index ],
                               NC_FLOAT, 1, 1, &dimensionId ) != -1 ) {

            UTCTimestamp timestamp;
            Line history = "";
            appendToLine( history, gasp->note );
            appendToLine( history, ",XDRConvert" );
            toUTCTimestamp( gasp->timestamps[ 0 ], timestamp );

            result = writeStandardContents( file, history,
                                            timestamp,
                                            dimensionId, 0, 0 );
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
         const GASP* gasp  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSData( Integer file, const GASP* gasp ) {

  PRE02( file != -1, isValidGASP( gasp ) );

  Integer result = 0;
  const Integer count = gasp->totalRegriddedPoints;

  if ( writeAllIntData( file, "column", count, 1, 1, 1,
                        gasp->columns ) ) {

    if ( writeAllIntData( file, "row", count, 1, 1, 1,
                          gasp->rows ) ) {

      if ( writeAllData( file, "longitude", count, 1, 1, 1,
                         gasp->gridLongitudes ) ) {

        if ( writeAllData( file, "latitude", count, 1, 1, 1,
                           gasp->gridLatitudes ) ) {

          const Integer index = gasp->variables >= 3 ? 2 : 0;

          if ( writeAllData( file, gasp->variable[ index ], count, 1, 1, 1,
                             gasp->gridData ) ) {

            timeData( gasp->timesteps, count, gasp->outputPoints,
                      gasp->gridData );

            result = writeAllData( file, "time", count, 1, 1, 1,
                                   gasp->gridData );
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
INPUTS:  GASP* gasp                  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPI( GASP* gasp,
                                    const Parameters* parameters ) {

  PRE02( isValidGASP( gasp ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    gasp->totalRegriddedPoints * 3 * 4 + 10000; /* lon, lat, var, hdr .*/
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {

    if ( writeRegriddedIOAPIHeader( file, gasp, parameters->grid ) ) {
      result = writeRegriddedIOAPIData( file, gasp, parameters->grid );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIHeader - Write header to file.
INPUTS:  Integer file      NetCDF file to write to.
         const GASP* gasp  Structure to write.
         const Grid* grid  Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          const GASP* gasp,
                                          const Grid* grid ) {

  PRE04( file != -1, isValidGASP( gasp ), grid, grid->invariant( grid ) );

  Integer result = 0;
  enum { VARIABLES = 3 }; /* LONGITUDE, LATITUDE, gasp. */
  const Integer index = gasp->variables >= 3 ? 2 : 0;
  Name variableNames[ VARIABLES ] = { "LONGITUDE", "LATITUDE", "gasp" };
  Name variableUnits[ VARIABLES ] = { "deg", "deg", "-" };
  const Integer firstTimestamp = fromUTCTimestamp( gasp->timestamp );
  Line history = "";
  appendToLine( history, gasp->note );
  appendToLine( history, ",XDRConvert" );
  strncpy( variableNames[ VARIABLES - 1 ], gasp->variable[ index ], 16 );
  variableNames[ VARIABLES - 1 ][ 16 ] = '\0';
  strncpy( variableUnits[ VARIABLES - 1 ], gasp->units[ index ], 16 );
  variableUnits[ VARIABLES - 1 ][ 16 ] = '\0';
  uppercase( variableNames[ VARIABLES - 1 ] );
  lowercase( variableUnits[ VARIABLES - 1 ] );

  result = writeM3IOHeader( file, gasp->timesteps, firstTimestamp, VARIABLES,1,
                            (const Name*) variableNames,
                            (const Name*) variableUnits,
                            history, grid );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIData - Write IOAPI-format data to file.
INPUTS:  Integer file      NetCDF file to write to.
         const GASP* gasp  Structure to write.
         const Grid* grid  Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file, const GASP* gasp,
                                        const Grid* grid ) {

  PRE04( file != -1, isValidGASP( gasp ), grid, grid->invariant( grid ) );

  Integer result = 0;
  const Integer rows     = grid->rows( grid );
  const Integer columns  = grid->columns( grid );
  const Integer cells    = rows * columns;
  Real* expandedGridData = NEW_ZERO( Real, cells );

  if ( expandedGridData ) {
    const Integer timesteps = gasp->timesteps;
    const Integer layers = 1;
    const Real scale = 1.0;
    const Integer index = gasp->variables >= 3 ? 2 : 0;
    Integer timestep = 0;
    Integer offset = 0;
    Name variable = "";
    memset( variable, 0, sizeof variable );
    strncpy( variable, gasp->variable[ index ], 16 );
    uppercase( variable );

    if ( writeM3IOGrid( grid, timesteps, layers, file ) ) {

      do {
        const Integer points = gasp->outputPoints[ timestep ];

        copyDataToGrid( points,
                        gasp->rows + offset,
                        gasp->columns + offset,
                        gasp->gridData + offset,
                        scale, 1, rows, columns,
                        expandedGridData );

        if ( ! writeM3IOData( file, variable,
                              timestep, 1, rows, columns,
                              expandedGridData ) ) {
          timestep = timesteps;
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
PURPOSE: writeRegriddedMCMC - Write regridded MCMC-format data.
INPUTS:  GASP* gasp       Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedMCMC( GASP* gasp,
                                   const Parameters* unused ) {

  PRE0( isValidGASP( gasp ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {

    /* Write header row: */

    output->writeString( output, "time,xcell,ycell,%s\n", gasp->variable[3] );

    if ( output->ok( output ) ) {
      const char* const dataFormat =
        "%5"INTEGER_FORMAT",%5"INTEGER_FORMAT",%5"INTEGER_FORMAT","
        "%28.18"REAL_F_FORMAT"\n";
      const Integer timesteps = gasp->timesteps;
      const Integer* columns  = gasp->columns;
      const Integer* rows     = gasp->rows;
      const Real* data        = gasp->gridData;
      Integer timestep = 0;

      /* Write data rows: */

      do {
        const Integer points = gasp->outputPoints[ timestep ];
        Integer point = 0;

        for ( point = 0; point < points; ++point ) {
          const Integer column = *columns++;
          const Integer row    = *rows++;
          const Real value     = *data++;

          output->writeString( output, dataFormat,
                               timestep + 1, column, row, value );

          if ( ! output->ok( output ) ) {
            point = points;
            timestep = timesteps;
          }
        }

        ++timestep;
      } while ( timestep < timesteps );
    }

    result = output->ok( output );
    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: regridGASP - Regrid data.
INPUTS:  Stream* input    Stream to read GASP data from.
         Integer method  E.g., AGGREGATE_MEAN.
         Grid* grid      Grid to project and aggregate points into.
         GASP* gasp    Data to regrid.
OUTPUTS: GASP* gasp    Regridded data.
******************************************************************************/

static void regridGASP(Stream* input, Integer method, Grid* grid, GASP* gasp) {

  PRE08( input, input->isReadable( input ),
         IS_VALID_AGGREGATE_METHOD( method ),
         grid,
         grid->invariant( grid ),
         isValidGASP( gasp ),
         gasp->totalRegriddedPoints == 0,
         gasp->longitudes == 0 );

  const Integer variables = gasp->variables;

  if ( variables >= 3 ) {
    const Integer timesteps     = gasp->timesteps;
    const Integer scans         = gasp->scans;
    const Integer maximumPoints = gasp->maximumPoints;
    const Integer outputSize    = MIN( scans, timesteps ) * maximumPoints;
    const Integer regridVariables = 8; /* lon,lat,dat,glon,glat,col,row,gdata*/
    const Integer dataSize      = outputSize * regridVariables + timesteps;
    FREE( gasp->data );
    gasp->data = NEW_ZERO( Real, dataSize );

    if ( gasp->data ) {
      Integer totalRegriddedPoints = 0;
      Integer timestep = 0;
      Integer yyyydddhh00   =
        ( fromUTCTimestamp( gasp->timestamp ) / 100 ) * 100;
      gasp->longitudes     = gasp->data                     + outputSize;
      gasp->latitudes      = gasp->longitudes               + outputSize;
      gasp->gridLongitudes = gasp->latitudes                + outputSize;
      gasp->gridLatitudes  = gasp->gridLongitudes           + outputSize;
      gasp->columns        = (Integer*) gasp->gridLatitudes + outputSize;
      gasp->rows           = gasp->columns                  + outputSize;
      gasp->gridData       = (Real*) gasp->rows             + outputSize;
      gasp->outputPoints   = (Integer*) gasp->gridData      + outputSize;

      do {
        Integer inputPoints = 0;

        if ( ! readScanDataForTimestamp( yyyydddhh00, input, gasp,
                                         &inputPoints ) ) {
          timestep = timesteps;
        } else if ( inputPoints > 0 ) {
          Integer outputPoints = 0;
          const Real minimumValidValue =
            OR3( ! strcmp( gasp->variable[3], "Optical_Depth_Land_And_Ocean"),
                 ! strcmp( gasp->variable[3], "Cloud_Optical_Thickness" ),
                 ! strcmp( gasp->variable[3], "Total_Ozone" ) ) ?
            0.0 : -900.0;

          grid->regrid( grid, method, minimumValidValue, inputPoints, 1,
                        gasp->longitudes, gasp->latitudes, 0, gasp->data,
                        0, /* no notes */
                        &outputPoints,
                        gasp->columns        + totalRegriddedPoints,
                        gasp->rows           + totalRegriddedPoints,
                        0,
                        gasp->gridLongitudes + totalRegriddedPoints,
                        gasp->gridLatitudes  + totalRegriddedPoints,
                        0,
                        gasp->gridData       + totalRegriddedPoints,
                        0 /* No regriddedNotes. */ );

          gasp->outputPoints[ timestep ] = outputPoints;
          totalRegriddedPoints += outputPoints;
        }

        incrementTimestamp( &yyyydddhh00 );
        ++timestep;
      } while ( timestep < timesteps );

      if ( timestep == timesteps ) {
        gasp->totalRegriddedPoints = totalRegriddedPoints;
      }
    }
  }

  POST02( gasp->totalRegriddedPoints >= 0,
          IMPLIES( gasp->totalRegriddedPoints > 0,
                   AND8( IN_RANGE( minimumItemI( gasp->outputPoints,
                                                 gasp->timesteps ),
                                   0, gasp->totalRegriddedPoints ),
                         IN_RANGE( maximumItemI( gasp->outputPoints,
                                                 gasp->timesteps ),
                                   1, gasp->totalRegriddedPoints ),
                         IN_RANGE( minimumItemI(gasp->columns,
                                                gasp->totalRegriddedPoints),
                                                1, grid->columns( grid ) ),
                         IN_RANGE( maximumItemI(gasp->columns,
                                                gasp->totalRegriddedPoints),
                                                1, grid->columns( grid ) ),
                         IN_RANGE( minimumItemI(gasp->rows,
                                                gasp->totalRegriddedPoints),
                                                1, grid->rows( grid ) ),
                         IN_RANGE( maximumItemI(gasp->rows,
                                                gasp->totalRegriddedPoints),
                                                1, grid->rows( grid ) ),
                         validLongitudesAndLatitudes(
                                                  gasp->totalRegriddedPoints,
                                                      gasp->gridLongitudes,
                                                      gasp->gridLatitudes ),
                         isNanFree( gasp->gridData,
                                    gasp->totalRegriddedPoints ) ) ) );
}



/******************************************************************************
PURPOSE: readScanDataForTimestamp - Read all data for a given timestamp for regridding.
INPUTS:  Integer  yyyydddhh00  Timestamp to compare.
         Stream*  input        Stream to read GASP data from.
         GASP*   gasp          GASP structure.
OUTPUTS: Integer* points       Number of points read for this timestamp.
         GASP*   gasp          gasp->data, longitudes, latitudes.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readScanDataForTimestamp( Integer yyyydddhh00, Stream* input,
                                         GASP* gasp, Integer* points ) {

  PRE09( isValidTimestamp( yyyydddhh00 ),
         input, input->isReadable( input ),
         gasp, isValidGASP( gasp ), points,
         gasp->longitudes, gasp->latitudes, gasp->data );

  const Integer scans = gasp->scans;
  Integer result = 0;
  Integer scan = 0;
  Integer offset = 0;
  *points = 0;

  do {
    const Integer yyyydddhhmm = gasp->timestamps[ scan ];
    const Integer yyyydddhh002 = ( yyyydddhhmm / 100 ) * 100;

    if ( yyyydddhh002 == yyyydddhh00 ) {
      const Integer count = gasp->points[ scan ];

      if ( ! readScanData( input, gasp->variables, count,
                            gasp->longitudes + offset,
                            gasp->latitudes  + offset,
                            gasp->data       + offset ) ) {
        scan = scans;
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
                                                      gasp->longitudes,
                                                      gasp->latitudes ),
                         isNanFree( gasp->data, *points ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readScanData - Read all variables for a given scan for regridding.
INPUTS:  Stream* input                Stream to read GASP data from.
         Integer variables           Number of scan variables.
         Integer points              Number of scan points to read.
OUTPUTS: Real* longitudes[ points ]  Longitude data.
         Real* latitudes[ points ]   Latitude data.
         Real* data[ points ]        3rd data variable. Others are skipped.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readScanData( Stream* input, Integer variables, Integer points,
                             Real longitudes[], Real latitudes[],
                             Real data[] ) {

  PRE07( input, input->isReadable( input ), variables >= 3, points > 0,
         longitudes, latitudes, data );

  Integer result = 0;
  Real* skipData = variables > 3 ? NEW_ZERO( Real, points ) : 0;

  if ( OR2( variables == 3, skipData ) ) {
    Integer variable = 0;

    do {
      Real* const output =
        variable == 0 ? longitudes :
        variable == 1 ? latitudes  :
        variable >  2 ? skipData   :
        data;

      input->read64BitReals( input, output, points );

      ++variable;
    } while ( AND2( input->ok( input ), variable < variables ) );

    result = input->ok( input );
  }

  FREE( skipData );

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND2( validLongitudesAndLatitudes( points,
                                                      longitudes, latitudes ),
                         isNanFree( data, points ) ) ) );
  return result;
}






/******************************************************************************
PURPOSE: Site.c - Define routines for processing Site data.

NOTES:

HISTORY: 2007/12 plessel.todd@epa.gov, Created.
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
  Line     note;
  UTCTimestamp timestamp;
  Integer  variables;            /* 1 (e.g., ozone) or 2 (windU, windV). */
  Integer  timesteps;
  Integer  stations;
  Integer* ids;                  /* ids[ stations ]. */
  Name*    variable;             /* variable[ variables ]. ["windU", "windV"]*/
  Name*    units;                /* units[ variables ]. ["m/s", "m/s"].*/
  Note*    notes;                /* notes[ stations ]. */
  Real*    sites;           /* sites[ stations * 2 = LONGITUDE, LATITUDE ]. */
  Real*    data;            /* data[ variables * timesteps * stations ]. */
  /* Regrid data: */
  Integer  totalRegriddedPoints; /* Total number of aggregated points. */
  Real*    stationLongitudes;    /* stationLongitudes[ station ]. */
  Real*    stationLatitudes;     /* stationLatitudes[ station ]. */
  Integer* stationColumns;       /* stationColumns[ station ]. */
  Integer* stationRows;          /* stationRows[ station ]. */
  Real*    stationXOffsets;      /* stationXOffsets[ station ]. */
  Real*    stationYOffsets;      /* stationYOffsets[ station ]. */
  Integer* outputColumns;        /* outputColumns[ outputPoint ]. */
  Integer* outputRows;           /* outputRows[ outputPoint ]. */
  Real*    outputLongitudes;     /* outputLongitudes[ outputPoint ]. */
  Real*    outputLatitudes;      /* outputLatitudes[ outputPoint ]. */
  Real*    outputData;           /* outputData[ outputPoint ]. */
  Integer* outputPoints;         /* outputPoints[ timestep ]. */
  Real     scale;                /* scale outputData by this factor. */
} Site;

typedef Integer (*Writer)( const Site* site, const Parameters* parameters );

typedef struct {
  Integer format;         /* FORMAT_XDR, etc. */
  Writer writer;          /* Routine that writes data in this format. */
  Writer regriddedWriter; /* Routine that writes regridded data in format. */
} Entry;


/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateSite( Site* site );

static Integer isValidSite( const Site* site );

static Integer isVectorVariable( const Site* site );

static Writer dispatcher( Integer format, Integer regrid );

static Integer readXDR( Stream* input, Site* site );

static Integer readXDRData( Stream* input, Site* site );

static Integer readRegriddedXDR( Stream* input, Site* site );

static Integer readRegriddedXDRData( Stream* input, Site* site );

static Integer readVariablesAndUnits2( Stream* input, Site* site );

static Integer compareRegriddedXDR( const Parameters* parameters, Site* site );

static Integer writeASCII( const Site* site, const Parameters* unused );

static Integer writeCOARDS( const Site* site, const Parameters* parameters );

static Integer writeCOARDSHeader( Integer file, const Site* site );

static Integer writeCOARDSData( Integer file, const Site* site );


static Integer writeRegriddedXDR( const Site* site,
                                  const Parameters* parameters );

static Integer writeRegriddedASCII( const Site* site,
                                    const Parameters* parameters );

static Integer writeRegriddedCOARDS( const Site* site,
                                     const Parameters* parameters );

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Site* site) ;

static Integer writeRegriddedCOARDSData( Integer file, const Site* site,
                                         const Parameters* parameters );

static Integer writeRegriddedIOAPI( const Site* site,
                                    const Parameters* parameters );

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Site* site,
                                          const Grid* grid );

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Site* site,
                                        const Grid* grid );

static void regridSite( Integer method, Grid* grid, Site* site );

static void compactIntegerData( Integer timesteps, Integer stations,
                                Integer totalPoints, const Integer points[],
                                const Integer input[], Integer output[] );

static void compactRealData( Integer timesteps, Integer stations,
                             Integer totalPoints, const Integer points[],
                             const Real input[], Real output[] );

static void copy2( Integer count, const Real input[],
                   Real output1[], Real output2[] );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateSite - Read input and write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void translateSite( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok,
         parameters->input->ok( parameters->input ) );

  Site site;
  ZERO_OBJECT( &site );
  site.scale = 1.0;
  parameters->ok = 0;

  if ( readXDR( parameters->input, &site ) ) {
    Writer writer = dispatcher( parameters->format, parameters->regrid );

    if ( ! writer ) {
      failureMessage( "Invalid/unsupported format/regrid specification." );
    } else if ( parameters->regrid ) {
      regridSite( parameters->regrid, parameters->grid, &site );

      if ( site.totalRegriddedPoints == 0 ) {
        failureMessage( "No points projected onto the grid." );
      } else {

        if ( parameters->aggregationTimesteps ) {
          const Integer isVector2 = isVectorVariable( &site );
          Integer dataVariable = site.variables - 1;
          Integer totalOutputPoints = 0;
          const Integer aggregatedTimesteps =
            aggregateData( parameters->aggregationTimesteps,
                           isVector2,
                           site.timesteps,
                           site.outputPoints,
                           site.outputLongitudes,
                           site.outputLatitudes,
                           0,
                           site.outputColumns,
                           site.outputRows,
                           0,
                           site.outputData,
                           0,
                           &totalOutputPoints );
          site.timesteps = aggregatedTimesteps;
          site.totalRegriddedPoints = totalOutputPoints;

          if ( AND2( parameters->aggregationTimesteps == 24,
                     ! OR2( strstr( site.variable[ dataVariable ], "daily" ),
                            strstr( site.variable[ dataVariable], "DAILY")))) {
            Integer count = 1 + isVector2;

            while ( count-- ) {
              Name dailyName = "";
              memset( dailyName, 0, sizeof dailyName );
              snprintf( dailyName, sizeof dailyName / sizeof *dailyName,
                        "daily_%s",
                        site.variable[ dataVariable ] );
              strncpy( site.variable[ dataVariable ], dailyName,
                       sizeof dailyName / sizeof *dailyName );
              --dataVariable;
            }
          }
        }

        parameters->ok = writer( &site, parameters );
      }
    } else {
      parameters->ok = writer( &site, parameters );
    }
  }

  deallocateSite( &site );
  POST0( isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: compareRegriddedSite - Read regridded-site input, compare it to
         CMAQ XDR data and write it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedSite( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok, parameters->input->ok( parameters->input ) );

  if ( ! AND3( ! parameters->regrid, parameters->compareFunction,
               parameters->data ) ) {
    failureMessage( "Invalid input for comparing." );
    parameters->ok = 0;
  } else {
    Site site;
    ZERO_OBJECT( &site );
    site.scale = 1.0;
    parameters->ok = 0;

    DEBUG( fprintf( stderr, "compareRegriddedSite()\n" ); )

    if ( readRegriddedXDR( parameters->input, &site ) ) {
      compareFunctionNameUnits( parameters->compareFunction,
                                parameters->convertFunction,
                                site.variable[ 0 ], site.units[ 0 ], 
                                parameters->variable,
                                parameters->units );

      if ( compareRegriddedXDR( parameters, &site ) ) {
        Writer writer = dispatcher( parameters->format, 1 );
        CHECK( writer );

        if ( site.totalRegriddedPoints == 0 ) {
          failureMessage( "No points projected onto the grid." );
        } else {
          parameters->ok = writer( &site, parameters );
        }
      }
    }

    deallocateSite( &site );
  }

  POST0( isValidParameters( parameters ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateSite - Deallocate contents of site structure.
INPUTS:  Site* site Structure to deallocate contents of.
******************************************************************************/

static void deallocateSite( Site* site ) {
  PRE0( site );

  if ( ! site->stationLongitudes ) { /* Read Regridded XDR compare case. */
    FREE( site->outputPoints );
    FREE( site->outputLongitudes );
  }

  FREE( site->variable );
  FREE( site->notes );
  FREE( site->ids );
  FREE( site->sites );
  FREE( site->stationLongitudes );
  ZERO_OBJECT( site );
}



/******************************************************************************
PURPOSE: isValidSite - Check site structure.
INPUTS:  const Site* site Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidSite( const Site* site ) {
  const Integer result =
    AND12( site,
           site->note[ 0 ],
           isValidUTCTimestamp( site->timestamp ),
           site->variables > 0,
           site->variable,
           site->variable[ 0 ],
           site->variable[ site->variables - 1 ],
           site->units,
           site->units[ 0 ],
           site->units[ site->variables - 1 ],
           site->timesteps > 0,
           IMPLIES_ELSE( site->stations > 0,
             AND8( site->notes,
                   site->notes[ 0 ],
                   site->notes[ site->stations - 1 ],
                   site->ids,
                   site->sites,
                   site->data,
                   ! isNan( site->scale ),
                   IMPLIES_ELSE( site->stationLongitudes,
                     AND12( site->stationLatitudes,
                            site->stationColumns,
                            site->stationRows,
                            site->stationXOffsets,
                            site->stationYOffsets,
                            site->outputColumns,
                            site->outputRows,
                            site->outputLongitudes,
                            site->outputLatitudes,
                            site->outputData,
                            site->outputPoints,
                            site->totalRegriddedPoints > 0 ),
                     IS_ZERO12( site->stationLatitudes,
                                site->stationColumns,
                                site->stationRows,
                                site->stationXOffsets,
                                site->stationYOffsets,
                                site->outputColumns,
                                site->outputRows,
                                site->outputLongitudes,
                                site->outputLatitudes,
                                site->outputData,
                                site->outputPoints,
                                site->totalRegriddedPoints ) ) ),
              AND7( site->totalRegriddedPoints > 0,
                    site->outputPoints,
                    site->outputLongitudes,
                    site->outputLatitudes,
                    site->outputColumns,
                    site->outputRows,
                    site->outputData ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isVectorVariable - Is the data variable a 2d wind vector?
INPUTS:  const Site* site Structure to check.
RETURNS: Integer 1 if vector, else 0.
******************************************************************************/

static Integer isVectorVariable( const Site* site ) {

  PRE05( site, site->variables > 0, site->variable,
         site->variable[ 0 ], site->variable[ site->variables ] );

  const Integer result =
    AND2( site->variables >= 2,
          OR2( AND2( ! strcmp( site->variable[ site->variables - 2 ], "windU" ),
                     ! strcmp( site->variable[ site->variables - 1 ], "windV" ) ),
               AND2( ! strcmp( site->variable[ site->variables - 2 ], "wind_u" ),
                     ! strcmp( site->variable[ site->variables - 1 ], "wind_v" ) ) ) );

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
PURPOSE: readXDR - Read XDR-format input and initialize site.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Site* site     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   File looks like this:
SITE 2.0
http://airnow.gov/,reformat_airnow_obs,SiteSubset
2005-08-26T00:00:00-0000
# data dimensions: timesteps stations
3 4
# Variable names:
OZONE
# Variable units:
ppb
# char notes[stations][80] and
# MSB 64-bit integers ids[stations] and
# IEEE-754 64-bit reals sites[stations][2=<longitude,latitude>] and
# IEEE-754 64-bit reals data[timesteps][stations]:
******************************************************************************/

static Integer readXDR( Stream* input, Site* site ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         site, site->sites == 0, site->data == 0 );

  Integer result = 0;
  input->readString( input, site->note, COUNT( site->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( site->note );

    if ( readTimestamp( input, site->timestamp ) ) {
      Integer dimensions[ 2 ] = { 0, 0 };

      if ( readDimensions( input, COUNT( dimensions ), dimensions ) ) {
        site->timesteps = dimensions[ 0 ];
        site->stations  = dimensions[ 1 ];

        if ( readVariablesAndUnits2( input, site ) ) {

          if ( skipInputLines( input, 4 ) ) {
            result = readXDRData( input, site );
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Site data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidSite( site ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read XDR-format data and initialize site.
INPUTS:  Stream* input  Input stream read from and parse.
         Site* site     site->variables,timesteps,stations.
OUTPUTS: Site* site     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDRData( Stream* input, Site* site ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         site, site->sites == 0, site->units[ 0 ],
         site->ids == 0, site->data == 0 );

  Integer result = 0;
  site->ids = NEW_ZERO( Integer, site->stations );

  if ( site->ids ) {
    const Integer siteCount = site->stations * 2;
    const Integer dataCount = site->variables * site->timesteps * site->stations;
    const Integer count     = siteCount + dataCount;
    site->sites = NEW_ZERO( Real, count );

    if ( site->sites ) {
      site->notes = NEW_ZERO( Note, site->stations );

      if ( site->notes ) {
        readNotes( input, site->stations, site->notes );

        if ( input->ok( input ) ) {
          input->read64BitIntegers( input, site->ids, site->stations );

          if ( input->ok( input ) ) {
            site->data = site->sites + siteCount;
            input->read64BitReals( input, site->sites, count );

            if ( input->ok( input ) ) {
              result = isValidSite( site );
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Site data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidSite( site ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDR - Read Regridded XDR-format input & initialize site.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Site* site     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Input data format is:

REGRIDDED-SITE 2.0
http://airnow.gov/,reformat_airnow_obs,SiteSubset,XDRConvert
2005-08-26T00:00:00-0000
# timesteps
3
# Variable name:
OZONE
# Variable units:
ppb
# lonlat projection: major_semiaxis minor_semiaxis
6.37e+06 6.37e+06
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[35]:
282 140 -130.000000 23.500000 0.250000 0.250000 5 5000.000000 0 40.385 80.9298 121.636 162.504 244.735 327.633 411.21 495.479 580.452 666.142 752.561 839.725 1016.34 1196.1 1379.14 1565.58 1755.57 2047.53 2348.35 2764.29 3310.41 3889.16 4505.11 5163.9 5872.63 6640.37 7479.09 8404.94 9440.67 10619.8 11995.2 13657.8 15789.1 18855.2
# MSB 32-bit integers points[timesteps] and
# IEEE-754 32-bit reals longitudes[timesteps][points] and
# IEEE-754 32-bit reals latitudes[timesteps][points] and
# MSB 32-bit integers columns[timesteps][points] and
# MSB 32-bit integers rows[timesteps][points] and
# IEEE-754 32-bit reals data[timesteps][points]:

******************************************************************************/

static Integer readRegriddedXDR( Stream* input, Site* site ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         site, site->sites == 0, site->data == 0 );

  Integer result = 0;
  input->readString( input, site->note, COUNT( site->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( site->note );

    if ( readTimestamp( input, site->timestamp ) ) {

      if ( readDimensions( input, 1, &site->timesteps ) ) {

        if ( readVariablesAndUnits2( input, site ) ) {

          if ( strcmp( site->units[ 0 ], "ppb" ) == 0 ) {
            site->scale = 0.001;
            strcpy( site->units[ 0 ], "ppm" ); /* Match CMAQ. */
          }

          {
            char line[ 256 ] = "";
            int count = 6;
            memset( line, 0, sizeof line );
            input->readString( input, line, sizeof line / sizeof *line - 1 );

            if ( strcmp( line,
                         "# MSB 32-bit integers points[timesteps] and\n" ) ) {
              count += 4; /* Skip 4 line projection/grid. */
            }

            if ( skipInputLines( input, count - 1 ) ) {
              site->outputPoints = NEW_ZERO( Integer, site->timesteps );

              if ( site->outputPoints ) {
                input->read32BitIntegers( input, site->outputPoints,
                                          site->timesteps );

                if ( input->ok( input ) ) {
                  site->totalRegriddedPoints =
                    sum( site->timesteps, site->outputPoints );
                  result = readRegriddedXDRData( input, site );
                }
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Site data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidSite( site ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRData - Read Regridded XDR-format input & init site.
INPUTS:  Stream* input  Input stream read from and parse.
OUTPUTS: Site* site     Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readRegriddedXDRData( Stream* input, Site* site ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         site, site->sites == 0, site->data == 0, site->outputData == 0,
         site->timesteps > 0 );

  Integer result = 0;
  const Integer count = site->totalRegriddedPoints;

  if ( count > 0 ) {
    const Integer isVector = isVectorVariable( site );
    site->outputLongitudes = NEW_ZERO( Real, count * ( 5 + isVector ) );

    if ( site->outputPoints ) {
      site->outputLatitudes = site->outputLongitudes + count;
      site->outputColumns   = (Integer*) site->outputLatitudes + count;
      site->outputRows      = site->outputColumns + count;
      site->outputData      = (Real*) site->outputRows + count;

      input->read32BitReals( input, site->outputLongitudes, count * 2 );

      if ( input->ok( input ) ) {
        input->read32BitIntegers( input, site->outputColumns, count * 2 );

        if ( input->ok( input ) ) {
          const size_t count2 = isVector ? count + count : count;
          input->read32BitReals( input, site->outputData, count2 );

          if ( input->ok( input ) ) {

            if ( site->scale ) {
              scale( site->scale, count2, site->outputData );
            }

            result = isValidSite( site );
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Site data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidSite( site ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readVariablesAndUnits2 - Read 1 (e.g., ozone) or 2 (windU windV)
         sets of variables and units.
INPUTS:  Stream* input Input stream read from and parse.
OUTPUTS: Site* site    Allocated site->variables,variable,units.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
NOTES:   Input data format is:
# Variable name:
windU windV
# Variable units:
m/s m/s
******************************************************************************/

static Integer readVariablesAndUnits2( Stream* input, Site* site ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         site, site->variables == 0, site->variable == 0 );

  Integer result = 0;
  char line[ 256 ] = "";
  memset( line, 0, sizeof line );
  input->readString( input, line, sizeof line / sizeof *line - 1 );

  if ( OR2( ! strcmp( line, "# Variable name:\n" ),
            ! strcmp( line, "# Variable names:\n" ) ) ) {
    input->readString( input, line, sizeof line / sizeof *line - 1 );
    site->variables = wordsInString( line );

    if ( IN3( site->variables, 1, 2 ) ) {
      site->variable = NEW_ZERO( Name, site->variables * 2 );

      if ( site->variable ) {
        site->units = site->variable + site->variables;

        if ( site->variables == 1 ) {
          result = sscanf( line, "%s\n", site->variable[ 0 ] ) == 1;
        } else {
          result = sscanf( line, "%s %s\n", site->variable[ 0 ],
                           site->variable[ 1 ] ) == 2;
        }

        if ( result ) {
          result = 0;
          input->readString( input, line, sizeof line / sizeof *line - 1 );

          if ( ! strcmp( line, "# Variable units:\n" ) ) {
            input->readString( input, line, sizeof line / sizeof *line - 1 );

            if ( wordsInString( line ) == site->variables ) {

              if ( site->variables == 1 ) {
                result = sscanf( line, "%s\n", site->units[ 0 ] ) == 1;
              } else {
                result = sscanf( line, "%s %s\n", site->units[ 0 ],
                                 site->units[ 1 ] ) == 2;
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    failureMessage( "Invalid SITE header (variables/units)." );
    site->variables = 0;
    FREE( site->variable );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND7( IN3( site->variables, 1, 2 ),
                         site->variable,
                         site->variable[ 0 ][ 0 ],
                         site->variable[ site->variables - 1 ][ 0 ],
                         site->units,
                         site->units[ 0 ][ 0 ],
                         site->units[ site->variables - 1 ][ 0 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: compareRegriddedXDR - Compare Regridded data with CMAQ data.
INPUTS:  const Parameters* parameters  CMAQ data to compare to.
OUTPUTS: Site* site                    Updated site->data.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static Integer compareRegriddedXDR( const Parameters* parameters, Site* site) {

  PRE05( parameters, isValidParameters( parameters ),
         parameters->compareFunction,
         site, isValidSite( site ) );

  Integer result = 0;

  if ( ! AND2( ! strcmp( parameters->timestamp, site->timestamp ),
               parameters->timesteps == site->timesteps ) ) {
    failureMessage( "Mismatched time steps (%s %lld)"
                    " for comparison to CMAQ data (%s %lld).",
                    site->timestamp, site->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else {
    Real* const siteData             = site->outputData;
    const Integer* const siteRows    = site->outputRows;
    const Integer* const siteColumns = site->outputColumns;
    const Integer* const sitePoints  = site->outputPoints;
    const Real* const cmaqData = parameters->data;
    CompareFunction comparer   = parameters->compareFunction;
    const Integer timesteps    = parameters->timesteps;
    const Integer firstRow     = parameters->firstRow;
    const Integer lastRow      = parameters->lastRow;
    const Integer firstColumn  = parameters->firstColumn;
    const Integer lastColumn   = parameters->lastColumn;
    const Integer rows         = lastRow    - firstRow    + 1;
    const Integer columns      = lastColumn - firstColumn + 1;
    const Integer rowsTimesColumns = rows * columns;
    Integer timestep = 0;
    Integer siteIndex = 0;

    DEBUG( fprintf( stderr, "timesteps = %lld\n", timesteps ); )

    for ( timestep = 0; timestep < timesteps; ++timestep ) {
      const Integer points = sitePoints[ timestep ];
      const Integer timestepOffset = timestep * rowsTimesColumns;
      Integer point = 0;

      DEBUG( fprintf( stderr, "timestep = %lld, points = %lld\n",
                      timestep, points ); )

      for ( point = 0; point < points; ++point, ++siteIndex ) {
        const Integer siteRow    = siteRows[    siteIndex ];
        const Integer siteColumn = siteColumns[ siteIndex ];
        DEBUG( fprintf( stderr, "  @[ %lld, %lld ]: ", siteRow, siteColumn ); )

        if ( AND2( IN_RANGE( siteRow, firstRow, lastRow ),
                   IN_RANGE( siteColumn, firstColumn, lastColumn ) ) ) {
          const Integer siteRow0    = siteRow    - firstRow;
          const Integer siteColumn0 = siteColumn - firstColumn;
          const Integer dataIndex =
            timestepOffset + siteRow0 * columns + siteColumn0;
          CHECK3( IN_RANGE( siteRow0, 0, rows - 1 ),
                  IN_RANGE( siteColumn0, 0, columns - 1 ),
                  IN_RANGE( dataIndex, 0, timesteps * rows * columns - 1 ) );
          const Real siteDatum = siteData[ siteIndex ];
          const Real cmaqDatum = cmaqData[ dataIndex ];
          const Real comparedDatum = comparer( siteDatum, cmaqDatum );
          siteData[ siteIndex ] = comparedDatum;
          result = 1;
          DEBUG( fprintf( stderr, "f(%lf, %lf) -> %lf\n",
                          siteDatum, cmaqDatum, comparedDatum ); )
        } else {
          siteData[ siteIndex ] = -9999.0;
          DEBUG( fprintf( stderr, "-9999\n" ); )
        }
      }
    }
  }

  POST02( IS_BOOL( result ), isValidSite( site ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format data.
INPUTS:  const Site* site              Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeASCII( const Site* site, const Parameters* unused ) {

  PRE0( isValidSite( site ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer isVector = isVectorVariable( site );
    const char* const headerStart =
      "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tSTATION(-)";

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {

      if ( isVector ) {
        const char* const headerFormat = "\t%s(%s)\t%s(%s)\n";
        output->writeString( output, headerFormat,
                             site->variable[ 0 ], site->units[ 0 ],
                             site->variable[ 1 ], site->units[ 1 ] );
      } else {
        const char* const headerFormat = "\t%s(%s)\n";
        output->writeString( output, headerFormat,
                             site->variable[ 0 ], site->units[ 0 ] );
      }

      if ( output->ok( output ) ) {
        const Integer timesteps = site->timesteps;
        const Integer stations = site->stations;
        const Integer totalPoints = timesteps * stations;
        Integer timestep = 0;
        Integer yyyydddhhmm = fromUTCTimestamp( site->timestamp );
        UTCTimestamp timestamp;

        /* Write data rows: */

        do {
          Integer station = 0;
          toUTCTimestamp( yyyydddhhmm, timestamp );

          do {
            const Integer id     = site->ids[ station ];
            const Integer station2 = station + station;
            const Real longitude = site->sites[ station2 + LONGITUDE ];
            const Real latitude  = site->sites[ station2 + LATITUDE  ];
            const Integer index  = timestep * stations + station;
            const Real data      = site->data[ index ];

            if ( isVector ) {
              const char* const dataFormat =
                "%s\t%10.5f\t%10.5f\t%20"INTEGER_FORMAT"\t%10.5f\t%10.5f\n";
              const Integer index2 = index + totalPoints;
              const Real data2 = site->data[ index2 ];
              output->writeString( output, dataFormat,
                                   timestamp, longitude, latitude, id,
                                   data, data2 );
           } else {
              const char* const dataFormat =
                "%s\t%10.5f\t%10.5f\t%20"INTEGER_FORMAT"\t%10.5f\n";
              output->writeString( output, dataFormat,
                                   timestamp, longitude, latitude, id, data );
            }

            if ( ! output->ok( output ) ) {
              station  = stations;
              timestep = timesteps;
            }

            ++station;
          } while ( station < stations );

          incrementTimestamp( &yyyydddhhmm );
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
PURPOSE: writeCOARDS - Write COARDS-format data.
INPUTS:  const Site* site              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDS( const Site* site,
                            const Parameters* parameters ) {

  PRE02( isValidSite( site ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    site->stations * 2 * 4 + /* longitude( station ), latitude( station ).*/
    site->variables * site->stations * site->timesteps * 4 +
    /* variable( station, time ).*/
    site->timesteps * 4 + /* time( time ). */
    1000; /* header. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file = createNetCDFFile( parameters->netcdfFileName, create64BitFile);

  if ( file != -1 ) {

    if ( writeCOARDSHeader( file, site ) ) {
      result = writeCOARDSData( file, site );
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
         const Site* site  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSHeader( Integer file, const Site* site ) {

  PRE02( file != -1, isValidSite( site ) );

  Integer result = 0;
  enum { TIME, STATION, DIMENSIONS };
  const char* const names[ DIMENSIONS ] = { "time", "station" };
  Integer dimensionIds[ DIMENSIONS ] = { -1, -1 };
  Integer dimensions[   DIMENSIONS ] = { 0, 0 };
  dimensions[ TIME    ] = site->timesteps;
  dimensions[ STATION ] = site->stations;

  if ( createDimensions( file, DIMENSIONS, names, dimensions,
                         dimensionIds ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "station_id", "-",
                           NC_INT, 0, 1, dimensionIds + STATION ) != -1 ) {

        if ( createLongitudeAndLatitude( file, 1, dimensionIds + STATION ) ) {

          if ( createVariable( file, site->variable[ 0 ], site->units[ 0 ],
                               NC_FLOAT, 1, 2, dimensionIds ) != -1 ) {
            const Integer isVector = isVectorVariable( site );

            if ( OR2( ! isVector,
                      createVariable( file, site->variable[1], site->units[1],
                                      NC_FLOAT, 1, 2, dimensionIds ) != -1 )) {
              Line history = "";
              appendToLine( history, site->note );
              appendToLine( history, ",XDRConvert" );
              result = writeStandardContents( file, history, site->timestamp,
                                              dimensionIds[ TIME ],
                                              site->timesteps, 1 );
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
         const Site* site  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSData( Integer file, const Site* site ) {

  PRE02( file != -1, isValidSite( site ) );

  Integer result = 0;
  const Integer stations = site->stations;

  if ( writeAllIntData( file, "station_id", stations, 1, 1, 1, site->ids ) ) {
    const Integer timesteps = site->timesteps;
    const Integer count = stations * MAX( timesteps, 2 );
    Real* data = NEW_ZERO( Real, count );

    if ( data ) {
      Integer index = 0;

      do {
        const Integer index2 = index + index;
        data[ index ] = site->sites[ index2 ];
        ++index;
      } while ( index < stations );

      if ( writeAllData( file, "longitude", stations, 1, 1, 1, data ) ) {
        index = 0;

        do {
          const Integer index2 = index + index;
          data[ index ] = site->sites[ index2 + 1 ];
          ++index;
        } while ( index < stations );

        if ( writeAllData( file, "latitude", stations, 1, 1, 1, data ) ) {
          const Integer totalPoints = timesteps * stations;
          const Integer variableBytes = totalPoints * sizeof *data;
          Integer variableOffset = 0;
          Integer variable = 0;

          /*
           * Copy variable data to data[timestep][station] because call to
           * writeAllData() will ruin it.
           */

          do {
            memcpy( data, site->data + variableOffset, variableBytes );
            result = writeAllData( file, site->variable[ variable ],
                                   timesteps, stations, 1, 1, data );
            ++variable;
            variableOffset += totalPoints;
          } while ( AND2( result, variable < site->variables ) );
        }
      }

      FREE( data );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  const Site* site              Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Output data format is:

REGRIDDED-SITE 2.0
http://www.epa.gov/ttn/airs/aqsdatamart,SiteSubset,XDRConvert
2001-08-26T00:00:00-0000
# timesteps
13
# Variable name:
ozone
# Variable units:
ppb
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6371000.000000 6371000.000000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
279 240 -1008000.000000 -1620000.000000 12000.000000 12000.000000 2 10000.000000 1 0.995
# MSB 32-bit integers points[timesteps] and
# IEEE-754 32-bit reals longitudes[timesteps][points] and
# IEEE-754 32-bit reals latitudes[timesteps][points] and
# MSB 32-bit integers columns[timesteps][points] and
# MSB 32-bit integers rows[timesteps][points] and
# IEEE-754 32-bit reals data[timesteps][points]:

******************************************************************************/

static Integer writeRegriddedXDR( const Site* site,
                                  const Parameters* parameters) {

  PRE02( isValidSite( site ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = site->timesteps;
    const Integer points    = site->totalRegriddedPoints;
    const Integer isVector = isVectorVariable( site );
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps
      : 1;
    Name variable = "";
    aggregateName( site->variable[ 0 ], hoursPerTimestep, variable );

    if ( isVector ) {
      Name variable2 = "";
      aggregateName( site->variable[ 1 ], hoursPerTimestep, variable2 );

      output->writeString( output,
                           "REGRIDDED-SITE 2.0\n"
                           "%s,XDRConvert\n"
                           "%s\n"
                           "# timesteps\n%"INTEGER_FORMAT"\n"
                           "# Variable name:\n%s %s\n"
                           "# Variable units:\n%s %s\n",
                           site->note,
                           site->timestamp,
                           timesteps,
                           site->variable, variable2,
                           site->units[ 0 ], site->units[ 1 ] );

    } else {
      output->writeString( output,
                           "REGRIDDED-SITE 2.0\n"
                           "%s,XDRConvert\n"
                           "%s\n"
                           "# timesteps\n%"INTEGER_FORMAT"\n"
                           "# Variable name:\n%s\n"
                           "# Variable units:\n%s\n",
                           site->note,
                           site->timestamp,
                           timesteps,
                           variable,
                           site->units[ 0 ] );
    }

    if ( output->ok( output ) ) {
      writeProjectionAndGrid( parameters->grid, output );

      if ( output->ok( output ) ) {
        output->writeString( output,
                    "# MSB 32-bit integers points[timesteps] and\n"
                    "# IEEE-754 32-bit reals longitudes[timesteps][points] and\n"
                    "# IEEE-754 32-bit reals latitudes[timesteps][points] and\n"
                    "# MSB 32-bit integers columns[timesteps][points] and\n"
                    "# MSB 32-bit integers rows[timesteps][points] and\n"
                    "# IEEE-754 32-bit reals data[timesteps][points]:\n" );

        if ( output->ok( output ) ) {
          output->write32BitIntegers( output, site->outputPoints, timesteps );

          if ( output->ok( output ) ) {
            output->write32BitReals( output, site->outputLongitudes, points );

            if ( output->ok( output ) ) {
              output->write32BitReals( output, site->outputLatitudes, points );

              if ( output->ok( output ) ) {
                output->write32BitIntegers( output, site->outputColumns, points);

                if ( output->ok( output ) ) {
                  output->write32BitIntegers( output, site->outputRows, points);

                  if ( output->ok( output ) ) {
                    const Integer points2 =
                      isVector ? points + points : points;
                    output->write32BitReals( output, site->outputData, points2);
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
INPUTS:  const Site* site              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedASCII( const Site* site,
                                    const Parameters* parameters ) {

  PRE02( isValidSite( site ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer isVector = isVectorVariable( site );
    const char* const headerStart =
      "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)"
      "\tCOLUMN(-)\tROW(-)";

    /* Write header row: */

    output->writeString( output, headerStart );

    if ( output->ok( output ) ) {
      const Integer hoursPerTimestep =
        parameters->aggregationTimesteps ? parameters->aggregationTimesteps
        : 1;
      Name variable = "";
      aggregateName( site->variable[ 0 ], hoursPerTimestep, variable );

      if ( isVector ) {
        const char* const headerFormat = "\t%s(%s)\t%s(%s)\n";
        Name variable2 = "";
        aggregateName( site->variable[ 1 ], hoursPerTimestep, variable2 );
        output->writeString( output, headerFormat,
                             variable,  site->units[ 0 ],
                             variable2, site->units[ 1 ] );
      } else {
        const char* const headerFormat = "\t%s(%s)\n";
        output->writeString( output, headerFormat,
                             variable, site->units[ 0 ] );
      }

      if ( output->ok( output ) ) {
        const Integer timesteps = site->timesteps;
        const Real* longitudes  = site->outputLongitudes;
        const Real* latitudes   = site->outputLatitudes;
        const Integer* columns  = site->outputColumns;
        const Integer* rows     = site->outputRows;
        const Real* data        = site->outputData;
        const Real* data2 =
          isVector ? site->outputData + site->totalRegriddedPoints : 0;
        const Integer hoursPerTimestep =
          parameters->aggregationTimesteps ? parameters->aggregationTimesteps
          : 1;
        Integer timestep = 0;
        Integer yyyydddhhmm = fromUTCTimestamp( site->timestamp );
        UTCTimestamp timestamp;

        /* Write data rows: */

        do {
          const Integer points = site->outputPoints[ timestep ];
          Integer point = 0;
          toUTCTimestamp( yyyydddhhmm, timestamp );

          for ( point = 0; point < points; ++point ) {
            const Real longitude = *longitudes++;
            const Real latitude  = *latitudes++;
            const Integer column = *columns++;
            const Integer row    = *rows++;
            const Real value     = *data++;

            if ( isVector ) {
              const char* const dataFormat =
                "%s\t%10.4lf\t%10.4lf"
                "\t%9"INTEGER_FORMAT"\t%9"INTEGER_FORMAT"\t%10.4lf\t%10.4lf\n";
              const Real value2 = *data2++;
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
INPUTS:  const Site* site              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDS( const Site* site,
                                     const Parameters* parameters ) {

  PRE02( isValidSite( site ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    site->totalRegriddedPoints * 7 * 4 + 10000;
    /* lon, lat, col, row, time, var, var2, + hdr. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedCOARDSHeader( file, hoursPerTimestep, site ) ) {
      result = writeRegriddedCOARDSData( file, site, parameters );
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
         Integer hoursPerTimestep  E.g., 1, 24, 744, etc.
         const Site* site          Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Site* site ) {

  PRE03( file != -1, hoursPerTimestep > 0, isValidSite( site ) );

  Integer result = 0;
  const char* const dimensionName = "points";
  Integer dimensionId = -1;
  const Integer dimension = site->totalRegriddedPoints;

  if ( createDimensions( file, 1, &dimensionName, &dimension, &dimensionId)) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "column", "-",
                           NC_INT, 0, 1, &dimensionId ) != -1 ) {

        if ( createVariable( file, "row", "-",
                             NC_INT, 0, 1, &dimensionId ) != -1 ) {

          if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {
            Name variable = "";
            aggregateName( site->variable[ 0 ], hoursPerTimestep, variable );

            if ( createVariable( file, variable, site->units[ 0 ],
                                 NC_FLOAT, 1, 1, &dimensionId ) != -1 ) {
              const Integer isVector = isVectorVariable( site );
              result = 1;

              if ( isVector ) {
                Name variable2 = "";
                aggregateName( site->variable[ 1 ], hoursPerTimestep, variable2);
                result =
                  createVariable( file, variable2, site->units[ 1 ],
                                  NC_FLOAT, 1, 1, &dimensionId ) != -1;
              }

              if ( result ) {
                Line history = "";
                appendToLine( history, site->note );
                appendToLine( history, ",XDRConvert" );
                result = writeStandardContents( file, history, site->timestamp,
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
         const Site* site  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSData( Integer file, const Site* site,
                                        const Parameters* parameters ) {

  PRE03( file != -1, isValidSite( site ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer count = site->totalRegriddedPoints;

  if ( writeAllIntData( file, "column", count, 1, 1, 1,
                        site->outputColumns ) ) {

    if ( writeAllIntData( file, "row", count, 1, 1, 1,
                          site->outputRows ) ) {

      if ( writeAllData( file, "longitude", count, 1, 1, 1,
                         site->outputLongitudes ) ) {

        if ( writeAllData( file, "latitude", count, 1, 1, 1,
                           site->outputLatitudes ) ) {
          const Integer hoursPerTimestep =
          parameters->aggregationTimesteps ? parameters->aggregationTimesteps
            : 1;
          Name variable = "";
          aggregateName( site->variable[ 0 ], hoursPerTimestep, variable );

          if ( writeAllData( file, variable, count, 1, 1, 1,
                             site->outputData ) ) {
            const Integer isVector = isVectorVariable( site );
            result = 1;

            if ( isVector ) {
              Name variable2 = "";
              aggregateName( site->variable[ 1 ], hoursPerTimestep, variable2);
              result =
                writeAllData( file, variable2, count, 1, 1, 1,
                              site->outputData + count );
            }

            if ( result ) {
              timeData( site->timesteps, hoursPerTimestep, count,
                        site->outputPoints, site->outputData );
              result = writeAllData( file, "time", count, 1, 1, 1,
                                     site->outputData );
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
INPUTS:  const Site* site              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPI( const Site* site,
                                    const Parameters* parameters ) {

  PRE02( isValidSite( site ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    site->totalRegriddedPoints * 4 * 4 + 10000; /* lon, lat, var, var2 + hdr */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedIOAPIHeader( file, hoursPerTimestep,
                                    site, parameters->grid ) ) {
      result =
        writeRegriddedIOAPIData(file, hoursPerTimestep, site, parameters->grid);
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
         const Site* site          Structure to write.
         const Grid* grid          Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Site* site,
                                          const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0,
         isValidSite( site ), grid, grid->invariant( grid ) );

  Integer result = 0;
  enum { VARIABLES = 3 }; /* LONGITUDE, LATITUDE, site. */
  Name variableNames[ VARIABLES + 1 ] = { "LONGITUDE", "LATITUDE", "site" "windV" };
  Name variableUnits[ VARIABLES + 1 ] = { "deg", "deg", "-" "m/s"};
  const Integer firstTimestamp = fromUTCTimestamp( site->timestamp );
  const Integer isVector = isVectorVariable( site );
  Line history = "";
  appendToLine( history, site->note );
  appendToLine( history, ",XDRConvert" );
  aggregateName( site->variable[ 0 ], hoursPerTimestep,
                 variableNames[ VARIABLES - 1 ] );
  variableNames[ VARIABLES - 1 ][ 15 ] = '\0';
  strncpy( variableUnits[ VARIABLES - 1 ], site->units[ 0 ], 16 );
  variableUnits[ VARIABLES - 1 ][ 16 ] = '\0';
  uppercase( variableNames[ VARIABLES - 1 ] );
  lowercase( variableUnits[ VARIABLES - 1 ] );

  if ( strcmp( variableUnits[ VARIABLES - 1 ], "ppb" ) == 0 ) {
    Site* hack = (Site*) site; /* UGLY HACK: const-cast-away! */
    hack->scale = 0.001;
    strcpy( variableUnits[ VARIABLES - 1 ], "ppm" ); /* Match CMAQ units. */
  }

  if ( isVector ) {
    aggregateName( site->variable[ 1 ], hoursPerTimestep,
                   variableNames[ VARIABLES ] );
    variableNames[ VARIABLES ][ 15 ] = '\0';
    strncpy( variableUnits[ VARIABLES ], site->units[ 1 ], 16 );
    variableUnits[ VARIABLES ][ 16 ] = '\0';
    uppercase( variableNames[ VARIABLES ] );
    lowercase( variableUnits[ VARIABLES ] );
  }

  result = writeM3IOHeader( file, site->timesteps, hoursPerTimestep,
                            firstTimestamp,
                            VARIABLES + isVector, 1,
                            (const Name*) variableNames,
                            (const Name*) variableUnits,
                            history, grid );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedIOAPIData - Write IOAPI-format data to file.
INPUTS:  Integer file      NetCDF file to write to.
         Integer hoursPerTimestep  Hours per timestep: 1, 24, etc.
         const Site* site  Structure to write.
         const Grid* grid  Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Site* site,
                                        const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep, isValidSite( site ),
         grid, grid->invariant( grid ) );

  Integer result = 0;
  const Integer layers  = 1;
  const Integer rows    = grid->rows( grid );
  const Integer columns = grid->columns( grid );
  const Integer cells   = layers * rows * columns;
  Real* gridData  = NEW_ZERO( Real, cells  );

  if ( gridData ) {
    const Integer timesteps = site->timesteps;
    Integer offset = 0;
    Integer timestep = 0;

    if ( writeM3IOGrid( grid, timesteps, layers, file ) ) {

      do {
        const Integer points = site->outputPoints[ timestep ];
        Name variable = "";
        aggregateName( site->variable[ 0 ], hoursPerTimestep, variable );
        variable[ 15 ] = '\0';
        uppercase( variable );

        copyDataToGrid( points,
                        site->outputRows    + offset,
                        site->outputColumns + offset,
                        site->outputData    + offset,
                        site->scale, layers, rows, columns, gridData );

        if ( ! writeM3IOData( file, variable,
                              timestep, layers, rows, columns, gridData ) ) {
          timestep = timesteps;
        } else {
          const Integer isVector = isVectorVariable( site );

          if ( isVector ) {
            const Real* outputData2 =
              site->outputData + site->totalRegriddedPoints;
            Name variable2 = "";
            aggregateName( site->variable[ 1 ], hoursPerTimestep, variable2 );
            variable2[ 15 ] = '\0';
            uppercase( variable2 );

            copyDataToGrid( points,
                            site->outputRows    + offset,
                            site->outputColumns + offset,
                            outputData2         + offset,
                            site->scale, layers, rows, columns, gridData );

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
PURPOSE: regridSite - Regrid data.
INPUTS:  Integer method  E.g., AGGREGATE_MEAN.
         Grid* grid      Grid to project and aggregate points into.
         Site* site      Data to regrid.
OUTPUTS: Site* site      Regridded data.
******************************************************************************/

static void regridSite( Integer method, Grid* grid, Site* site ) {

  PRE05( IS_VALID_AGGREGATE_METHOD( method ),
         grid,
         grid->invariant( grid ),
         isValidSite( site ),
         site->totalRegriddedPoints == 0 );

  Integer totalRegriddedPoints = 0;
  const Integer isVector = isVectorVariable( site );
  const Integer timesteps = site->timesteps;
  const Integer stations = site->stations;
  const Integer inputCount = stations * 8;/*lon,lat,col,row,xOff,yOff,glon/at*/
  const Integer maximumOutputPoints = timesteps * stations;
  const Integer outputCount = maximumOutputPoints * ( 5 + isVector );
  /* lon, lat, col, row, data, data2 */
  const Integer count = inputCount + outputCount + timesteps;
  site->stationLongitudes = NEW_ZERO( Real, count );

  if ( site->stationLongitudes ) {
    Real* gridLongitudes = 0;
    Real* gridLatitudes  = 0;
    site->stationLatitudes = site->stationLongitudes + stations;
    site->stationColumns   = (Integer*) site->stationLatitudes + stations;
    site->stationRows      = site->stationColumns + stations;
    site->stationXOffsets  = (Real*) site->stationRows + stations;
    site->stationYOffsets  = site->stationXOffsets + stations;
    gridLongitudes         = site->stationYOffsets + stations;
    gridLatitudes          = gridLongitudes + stations;
    site->outputColumns    = (Integer*) gridLatitudes + stations;
    site->outputRows       = site->outputColumns + maximumOutputPoints;
    site->outputLongitudes = (Real*) site->outputRows + maximumOutputPoints;
    site->outputLatitudes  = site->outputLongitudes + maximumOutputPoints;
    site->outputData       = site->outputLatitudes + maximumOutputPoints;
    site->outputPoints     = (Integer*)
      site->outputData + maximumOutputPoints * ( 1 + isVector );

    copy2( stations, site->sites,
           site->stationLongitudes, site->stationLatitudes );

    grid->projectXY( grid, stations,
                     site->stationLongitudes, site->stationLatitudes,
                     site->outputPoints,
                     site->stationColumns, site->stationRows,
                     site->stationXOffsets, site->stationYOffsets,
                     gridLongitudes, gridLatitudes );

    if ( site->outputPoints[ 0 ] ) {
      const Integer projectedSiteCount = site->outputPoints[ 0 ];
      const Real minimumValidValue = 0.0;
      const Real* const xOffsets = site->stationXOffsets;
      const Real* const yOffsets = site->stationYOffsets;
      Integer timestep = 0;

#pragma omp parallel for reduction( + : totalRegriddedPoints )

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        Integer* const outputPoints = site->outputPoints + timestep;
        const Integer offset        = timestep * stations;
        const Real* const inputData = site->data             + offset;
        const Real* const inputData2 =
          isVector ? inputData + maximumOutputPoints : 0;
        Integer* const columns      = site->outputColumns    + offset;
        Integer* const rows         = site->outputRows       + offset;
        Real* const longitudes      = site->outputLongitudes + offset;
        Real* const latitudes       = site->outputLatitudes  + offset;
        Real* const outputData      = site->outputData       + offset;
        Real* const outputData2 =
          isVector ? outputData + maximumOutputPoints : 0;

        *outputPoints = projectedSiteCount;
        memcpy( columns, site->stationColumns, stations * sizeof *columns   );
        memcpy( rows,    site->stationRows,    stations * sizeof *rows      );
        memcpy( longitudes, gridLongitudes,    stations * sizeof *longitudes);
        memcpy( latitudes,  gridLatitudes,     stations * sizeof *latitudes );

        grid->aggregate( grid, method, minimumValidValue, stations,
                         columns, rows, xOffsets, yOffsets,
                         longitudes, latitudes, 1, 0, inputData, inputData2,
                         outputPoints, outputData, outputData2, 0 );

        totalRegriddedPoints += *outputPoints;
      }

      /* Compact the output arrays so all values are contiguous: */

      compactIntegerData( timesteps, stations, totalRegriddedPoints,
                          site->outputPoints,
                          site->outputColumns, site->outputColumns );

      compactIntegerData( timesteps, stations, totalRegriddedPoints,
                          site->outputPoints,
                          site->outputRows, site->outputRows );

      compactRealData( timesteps, stations, totalRegriddedPoints,
                       site->outputPoints,
                       site->outputLongitudes, site->outputLongitudes );

      compactRealData( timesteps, stations, totalRegriddedPoints,
                       site->outputPoints,
                       site->outputLatitudes, site->outputLatitudes );

      compactRealData( timesteps, stations, totalRegriddedPoints,
                       site->outputPoints,
                       site->outputData, site->outputData );

      if ( isVector ) {
        compactRealData( timesteps, stations, totalRegriddedPoints,
                         site->outputPoints,
                         site->outputData + maximumOutputPoints,
                         site->outputData + totalRegriddedPoints );
      }
    }
  }

  site->totalRegriddedPoints = totalRegriddedPoints;

  POST02( site->totalRegriddedPoints >= 0,
          IMPLIES( site->totalRegriddedPoints > 0,
                   AND6( IN_RANGE( minimumItemI(site->outputColumns,
                                                site->totalRegriddedPoints),
                                                1, grid->columns( grid ) ),
                         IN_RANGE( maximumItemI(site->outputColumns,
                                                site->totalRegriddedPoints),
                                                1, grid->columns( grid ) ),
                         IN_RANGE( minimumItemI(site->outputRows,
                                                site->totalRegriddedPoints),
                                                1, grid->rows( grid ) ),
                         IN_RANGE( maximumItemI(site->outputRows,
                                                site->totalRegriddedPoints),
                                                1, grid->rows( grid ) ),
                         validLongitudesAndLatitudes(site->totalRegriddedPoints,
                                                      site->outputLongitudes,
                                                      site->outputLatitudes ),
                         isNanFree( site->outputData,
                                    site->totalRegriddedPoints
                                    * ( 1 + isVector ) ) ) ) );
}



/******************************************************************************
PURPOSE: compactIntegerData - Compact integer data as contiguous storage.
INPUTS:  Integer timesteps                       Number of timesteps.
         Integer stations                        Number of stations.
         Integer totalPoints                     Total number of points.
         const Integer points[ timesteps ]       Number of points per timstep.
         const Integer input[ timesteps * stations ] Data to read.
OUTPUTS: Integer output[ totalPoints ]           Contiguous data.
******************************************************************************/

static void compactIntegerData( Integer timesteps, Integer stations,
                                Integer totalPoints, const Integer points[],
                                const Integer input[], Integer output[] ) {

  PRE08( GT_ZERO3( timesteps, stations, totalPoints ),
         totalPoints <= timesteps * stations,
         points,
         minimumItemI( points, timesteps ) >= 0,
         maximumItemI( points, timesteps ) <= totalPoints,
         input,
         minimumItemI( input, timesteps * stations ) >= 0,
         output );

  Integer outputIndex = 0;
  Integer timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    Integer inputIndex = timestep * stations;
    Integer count = points[ timestep ];

    while ( count-- ) {
      output[ outputIndex++ ] = input[ inputIndex++ ];
    }
  }

  CHECK( outputIndex == totalPoints );

  POST0( minimumItemI( output, totalPoints ) >= 1 );
}



/******************************************************************************
PURPOSE: compactRealData - Compact real data as contiguous storage.
INPUTS:  Integer timesteps                        Number of timesteps.
         Integer stations                         Number of stations.
         Integer totalPoints                      Total number of points.
         const Integer points[ timesteps ]        Number of points per timstep.
         const Real input[ timesteps * stations ] Data to read.
OUTPUTS: Real output[ totalPoints ]               Contiguous data.
******************************************************************************/

static void compactRealData( Integer timesteps, Integer stations,
                             Integer totalPoints, const Integer points[],
                             const Real input[], Real output[] ) {

  PRE07( GT_ZERO3( timesteps, stations, totalPoints ),
         points,
         minimumItemI( points, timesteps ) >= 0,
         maximumItemI( points, timesteps ) <= totalPoints,
         input,
         isNanFree( input, timesteps * stations ),
         output );

  Integer outputIndex = 0;
  Integer timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    Integer inputIndex = timestep * stations;
    Integer count = points[ timestep ];

    while ( count-- ) {
      output[ outputIndex++ ] = input[ inputIndex++ ];
    }
  }

  CHECK( outputIndex == totalPoints );

  POST0( isNanFree( output, totalPoints ) );
}



/******************************************************************************
PURPOSE: copy2 - Copy 2d data to 1d arrays.
INPUTS:  Integer count                  Number of pairs.
         const Real input[ count * 2 ]  Pairs to copy.
OUTPUTS: Real output1[ count ]          1st item per pair.
         Real output2[ count ]          2nd item per pair.
******************************************************************************/

static void copy2( Integer count, const Real input[],
                   Real output1[], Real output2[] ) {

  PRE05( count > 0, input, isNanFree( input, count * 2 ), output1, output2);

  const Integer count2 = count + count;
  Integer index = 0;

  for ( index = 0; index < count2; index += 2 ) {
    const Integer outputIndex = index / 2;
    output1[ outputIndex ] = input[ index ];
    output2[ outputIndex ] = input[ index + 1 ];
  }

  POST04( output1[ 0 ] == input[ 0 ], output2[ 0 ] == input[ 1 ],
          output1[ count - 1 ] == input[ 2 * count - 2 ],
          output2[ count - 1 ] == input[ 2 * count - 1 ] );
}








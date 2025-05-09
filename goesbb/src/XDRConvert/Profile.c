
/******************************************************************************
PURPOSE: Profile.c - Define routines for processing Profile data.

NOTES:   Input profile data in XDR format is as follows:

Profile 2.0
http://www.esrl.noaa.gov/gmd/grad/neubrew/,NEUBrewSubset
2011-07-12T00:00:00-0000 2011-07-13T23:59:59-0000
# Subset domain: <min_lon> <min_lat> <max_lon> <max_lat>:
-90 25 -60 50
# Dimensions: variables profiles:
6 6
# Variable names:
timestamp id longitude latitude elevation ozone
# Variable units:
yyyymmddhhmmss - deg deg m molecules/cm3
# char notes[profiles][80] and
# MSB 64-bit integers points[profiles] and
# IEEE-754 64-bit reals data_1[variables][points_1] ...
 data_P[variables][points_T]:
 
HISTORY: 2011/08 plessel.todd@epa.gov, Created.
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

enum {
  DATA_TIMESTAMP, DATA_ID, DATA_LONGITUDE, DATA_LATITUDE, DATA_ELEVATION,
  DATA_OZONE, IMPLICIT_VARIABLES = 5
};

typedef struct {
  /* Input data: */
  Line         note;           /* File note/description. */
  UTCTimestamp firstTimestamp; /* Earliest timestamp of data. */
  UTCTimestamp lastTimestamp;  /* Latest   timestamp of data. */
  Bounds       bounds;         /*bounds[LONGITUDE LATITUDE][MINIMUM MAXIMUM]*/
  Integer      variables;      /* 6 = timestamp,id,lon,lat,elv,ozone. */
  Integer      totalPoints;    /* Sum of points[ profile ]. */
  Integer      profiles;       /* E.g., 2 profile profiles. */
  Name*        variable;       /* variable[ variables ]. E.g., "ozone"*/
  Name*        units;          /* units[ variables ]. E.g., "molecules/cm3". */
  Integer*     points;         /* points[ profiles ] */
  Real*        data;           /* data_1[ variables ][ points_1 ] ... */
                               /* data_P[ variables ][ points_P ] */
  /* Regridded data: */
  Integer      totalRegriddedPoints; /* Total number of regridded points. */
  Integer      timesteps;      /* Hours in regridded output. */
  Integer*     timestamps;     /* timestamps[ timesteps ]. */
  Integer*     outputPoints;   /* outputPoints[ timesteps ]. */
  Real*        longitudes;    /*longitudes[MIN(profiles,timesteps)*maxPoints]*/
  Real*        latitudes;     /*latitudes[ MIN(profiles,timesteps)*maxPoints]*/
  Real*        elevations;    /*elevations[MIN(profiles,timesteps)*maxPoints]*/
  Real*        gridLongitudes; /* gridLongitudes[ totalRegriddedPoints ]. */
  Real*        gridLatitudes;  /* gridLatitudes[ totalRegriddedPoints ]. */
  Real*        gridElevations; /* gridElevations[ totalRegriddedPoints ]. */
  Integer*     columns;        /* columns[ totalRegriddedPoints ]. */
  Integer*     rows;           /* rows[ totalRegriddedPoints ]. */
  Integer*     layers;         /* layers[ totalRegriddedPoints ]. */
  Real*        copyData;       /* copyData[ totalPoints ]. */
  Real*        gridData;       /* gridData[ totalRegriddedPoints ]. */
} Profile;

typedef Integer (*Writer)( Profile* profile, const Parameters* parameters );

typedef struct {
  Integer format;         /* FORMAT_XDR, etc. */
  Writer writer;          /* Routine that writes data in this format. */
  Writer regriddedWriter; /* Routine that writes regridded data in format. */
} Entry;

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateProfile( Profile* profile );

static Integer isValidProfile( const Profile* profile );

static Writer dispatcher( Integer format, Integer regrid );

static Integer readXDR( Stream* input, Profile* profile );

static Integer readXDRData( Stream* input, Profile* profile );

static Integer readRegriddedXDR( Stream* input, Profile* profile );

static Integer readRegriddedXDRData( Stream* input, Profile* profile );

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    Profile* profile );

static Integer writeASCII( Profile* profile, const Parameters* parameters );

static void writeASCIIHeader( const Profile* profile, Stream* output );

static Integer writeASCIIData( Profile* profile, Stream* output );

static Integer writeCOARDS( Profile* profile, const Parameters* parameters );

static Integer writeCOARDSHeader( Integer file, const Profile* profile );

static Integer writeCOARDSData( Integer file, Profile* profile );


static Integer writeRegriddedXDR( Profile* profile,
                                  const Parameters* parameters );

static Integer writeRegriddedASCII( Profile* profile,
                                    const Parameters* parameters );

static Integer writeRegriddedCOARDS( Profile* profile,
                                     const Parameters* parameters );

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Profile* profile );

static Integer writeRegriddedCOARDSData( Integer file,
                                         const Profile* profile,
                                         const Parameters* parameters );

static Integer writeRegriddedIOAPI( Profile* profile,
                                    const Parameters* parameters );

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Profile* profile,
                                          const Grid* grid );

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Profile* profile,
                                        const Grid* grid );

static void regridProfile( Integer method, Grid* grid, Profile* profile );


static Integer copyDataForTimestamp( Integer yyyydddhh00, Profile* profile );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateProfile - Read input & write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void translateProfile( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok,
         parameters->input->ok( parameters->input ) );

  Profile profile;
  ZERO_OBJECT( &profile );
  parameters->ok = 0;

  if ( readXDR( parameters->input, &profile ) ) {
    Writer writer = dispatcher( parameters->format, parameters->regrid );

    if ( ! writer ) {
      failureMessage( "Invalid/unsupported format/regrid specification." );
    } else if ( parameters->regrid ) {
      regridProfile( parameters->regrid, parameters->grid, &profile );

      if ( profile.totalRegriddedPoints == 0 ) {
        failureMessage( "No points projected onto the grid." );
      } else {

        if ( parameters->aggregationTimesteps ) {
          const Integer dataVariable = profile.variables - 1;
          Integer totalOutputPoints = 0;
          const Integer aggregatedTimesteps =
            aggregateData( parameters->aggregationTimesteps,
                           0,
                           profile.timesteps,
                           profile.outputPoints,
                           profile.gridLongitudes,
                           profile.gridLatitudes,
                           profile.gridElevations,
                           profile.columns,
                           profile.rows,
                           profile.layers,
                           profile.gridData,
                           0,
                           &totalOutputPoints );
          profile.timesteps = aggregatedTimesteps;
          profile.totalRegriddedPoints = totalOutputPoints;

          if ( AND2( parameters->aggregationTimesteps == 24,
                     ! OR2( strstr( profile.variable[ dataVariable ], "daily" ),
                            strstr( profile.variable[ dataVariable], "DAILY")))) {
            Name dailyName = "";
            memset( dailyName, 0, sizeof dailyName );
            snprintf( dailyName, sizeof dailyName / sizeof *dailyName,
                      "daily_%s",
                      profile.variable[ dataVariable ] );
            strncpy( profile.variable[ dataVariable ], dailyName,
                     sizeof dailyName / sizeof *dailyName );
          }
        }

        parameters->ok = writer( &profile, parameters );
      }
    } else {
      parameters->ok = writer( &profile, parameters );
    }
  }

  deallocateProfile( &profile );
  POST0( isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: compareRegriddedProfile - Read REGRIDDED-Profile input,
         compare it to CMAQ XDR data & write it in the given format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  parameters->ok updated by translation.
******************************************************************************/

void compareRegriddedProfile( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok, parameters->input->ok( parameters->input ) );

  if ( ! AND2( parameters->compareFunction, parameters->data ) ) {
    failureMessage( "Invalid input for comparing." );
    parameters->ok = 0;
  } else {
    Profile profile;
    ZERO_OBJECT( &profile );
    parameters->ok = 0;

    if ( readRegriddedXDR( parameters->input, &profile ) ) {
      compareFunctionNameUnits( parameters->compareFunction,
                                parameters->convertFunction,
                                profile.variable[ 0 ], profile.units[ 0 ],
                                parameters->variable, parameters->units );

      if ( compareRegriddedXDR( parameters, &profile ) ) {
        Writer writer = dispatcher( parameters->format, 1 );
        CHECK( writer );

        if ( profile.totalRegriddedPoints == 0 ) {
          failureMessage( "No points projected onto the grid." );
        } else {
          parameters->ok = writer( &profile, parameters );
        }
      }
    }

    deallocateProfile( &profile );
  }

  POST0( isValidParameters( parameters ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateProfile - Deallocate contents of profile structure.
INPUTS:  Profile* profile Structure to deallocate contents of.
******************************************************************************/

static void deallocateProfile( Profile* profile ) {
  PRE0( profile );
  FREE( profile->variable );
  FREE( profile->points );
  FREE( profile->data );
  FREE( profile->longitudes );
  ZERO_OBJECT( profile );
}



/******************************************************************************
PURPOSE: isValidProfile - Check profile structure.
INPUTS:  const Profile* profile Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidProfile( const Profile* profile ) {
  Integer result =
    AND6( profile,
          profile->note[ 0 ],
          isValidUTCTimestamp( profile->firstTimestamp ),
          IMPLIES_ELSE( profile->variables > IMPLICIT_VARIABLES,
                        AND2( profile->variable[ IMPLICIT_VARIABLES ],
                              profile->units[ IMPLICIT_VARIABLES ] ),
                        AND2( profile->variable[ 0 ],
                              profile->units[ 0 ] ) ),
          IMPLIES( GT_ZERO2( profile->profiles, profile->totalPoints ),
            AND10( isValidUTCTimestamp( profile->lastTimestamp ),
                   profile->variables > IMPLICIT_VARIABLES,
                   isValidBounds( (const Real (*)[2]) profile->bounds ),
                   profile->variable,
                   profile->units,
                   profile->points,
                   minimumItemI( profile->points, profile->profiles ) > 0,
                   profile->data,
                   isNanFree( profile->data,
                              profile->variables * profile->totalPoints ),
                   profile->totalRegriddedPoints >= 0 ) ),
          IMPLIES( profile->totalRegriddedPoints > 0,
                   AND16( profile->timesteps > 0,
                          profile->outputPoints,
                          minimumItemI( profile->outputPoints,
                                        profile->timesteps ) >= 0,
                          profile->columns,
                          profile->rows,
                          profile->layers,
                          profile->gridLongitudes,
                          profile->gridLatitudes,
                          profile->gridElevations,
                          profile->gridData,
                          minimumItemI( profile->columns,
                                        profile->totalRegriddedPoints ) > 0,
                          minimumItemI( profile->rows,
                                        profile->totalRegriddedPoints ) > 0,
                          minimumItemI( profile->layers,
                                        profile->totalRegriddedPoints ) > 0,
                          isNanFree( profile->gridElevations,
                                     profile->totalRegriddedPoints ),
                          isNanFree( profile->gridData,
                                     profile->totalRegriddedPoints ),
                          validLongitudesAndLatitudes(
                            profile->totalRegriddedPoints,
                            profile->gridLongitudes,
                            profile->gridLatitudes ) ) ) );

  if ( AND2( result, GT_ZERO2( profile->profiles, profile->totalPoints ) ) ) {
    const Integer variables = profile->variables;
    const Integer profiles = profile->profiles;
    const Real* data = profile->data;
    Integer profileIndex = 0;

    for ( profileIndex = 0; AND2( result, profileIndex < profiles );
          ++profileIndex ) {
      const Integer points = profile->points[ profileIndex ];
      const Real* const timestamps = data + DATA_TIMESTAMP * points;
      const Real* const ids        = data + DATA_ID        * points;
      const Real* const longitudes = data + DATA_LONGITUDE * points;
      const Real* const latitudes  = data + DATA_LATITUDE  * points;
      const Real* const elevations = data + DATA_ELEVATION * points;
      const Real* const ozones     = data + DATA_OZONE     * points;
      Integer point = 0;
      data += variables * points;

      for ( point = 0; AND2( result, point < points ); ++point ) {
        const Integer timestamp = (Integer) timestamps[ point ];
        const Integer id        = (Integer) ids[ point ];
        const Real longitude    = longitudes[ point ];
        const Real latitude     = latitudes[ point ];
        const Real elevation    = elevations[ point ];
        const Real ozone        = ozones[ point ];
        result = isValidYYYYMMDDHHMMSS( timestamp );
        result = AND2( result, id > 0 );
        result = AND2( result, isValidLongitude( longitude ) );
        result = AND2( result, isValidLatitude( latitude ) );
        result = AND3( result, ! isNan( elevation ),
                       IN_RANGE( elevation, -500.0, 110000.0 ) );
        result = AND3( result, ! isNan( ozone ), ozone >= 0.0 );
      }
    }
  }

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
PURPOSE: readXDR - Read input and initialize profile structure.
INPUTS:  Stream* input      Input stream read from and parse.
OUTPUTS: Profile* profile Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Input profile data in XDR format is:

Profile 2.0
http://www.esrl.noaa.gov/gmd/grad/neubrew/,NEUBrewSubset
2011-07-12T00:00:00-0000 2011-07-13T23:59:59-0000
# Subset domain: <min_lon> <min_lat> <max_lon> <max_lat>:
-90 25 -60 50
# Dimensions: variables profiles:
6 6
# Variable names:
timestamp id longitude latitude elevation ozone
# Variable units:
yyyymmddhhmmss - deg deg m molecules/cm3
# char notes[profiles][80] and
# MSB 64-bit integers points[profiles] and
# IEEE-754 64-bit reals data_1[variables][points_1] ...
 data_P[variables][points_T]:

******************************************************************************/

static Integer readXDR( Stream* input, Profile* profile ) {

  PRE06( input, input->ok( input ), input->isReadable( input ),
         profile, profile->variable == 0, profile->data == 0 );

  Integer result = 0;

  input->readString( input, profile->note, COUNT( profile->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( profile->note );

    if ( readTimestamps( input,
                         profile->firstTimestamp,
                         profile->lastTimestamp ) ) {

      if ( readDomain( input, profile->bounds ) ) {
        Integer dimensions[ 2 ] = { 0, 0 };

        if ( readDimensions( input, COUNT( dimensions ), dimensions ) ) {
          profile->variables = dimensions[ 0 ];
          profile->profiles  = dimensions[ 1 ];
          profile->variable  = NEW_ZERO( Name, profile->variables * 2 );

          if ( profile->variable ) {
            profile->units = profile->variable + profile->variables;

            if ( readVariablesAndUnits( input, profile->variables,
                                        profile->variable,
                                        profile->units ) ) {

              if ( skipInputLines( input, 3 + profile->profiles ) ) {
                result = readXDRData( input, profile );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Profile data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidProfile( profile ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRData - Read binary data from input.
INPUTS:  Stream* input      Input stream read from and parse.
OUTPUTS: Profile* profile Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDRData( Stream* input, Profile* profile ) {

  PRE09( input,
         input->ok( input ),
         input->isReadable( input ),
         profile,
         profile->variables > IMPLICIT_VARIABLES,
         profile->profiles > 0,
         profile->totalPoints == 0,
         profile->points == 0,
         profile->data == 0 );

  Integer result = 0;
  profile->points = NEW_ZERO( Integer, profile->profiles );

  if ( profile->points ) {
    input->read64BitIntegers( input, profile->points, profile->profiles );

    if ( input->ok( input ) ) {
      const Integer totalPoints = sum( profile->profiles, profile->points );
      const Integer dataSize = profile->variables * totalPoints;
      profile->data = NEW_ZERO( Real, dataSize );

      if ( profile->data ) {
        profile->totalPoints = totalPoints;
        input->read64BitReals( input, profile->data, dataSize );

        if ( input->ok( input )  ) {
          result = isValidProfile( profile );
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Profile data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidProfile( profile ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDR - Read REGRIDDED-Profile and initialize profile.
INPUTS:  Stream* input      Input stream read from and parse.
OUTPUTS: Profile* profile Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
NOTES:   Input data format is:

REGRIDDED-Profile 2.0
http://www.esrl.noaa.gov/gmd/grad/neubrew/,NEUBrewSubset,XDRConvert
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
# IEEE-754 64-bit reals longitudes[points] and
# IEEE-754 64-bit reals latitudes[points] and
# IEEE-754 64-bit reals elevations[points] and
# MSB 64-bit integers columns[points] and
# MSB 64-bit integers rows[points] and
# MSB 64-bit integers layers[points] and
# IEEE-754 64-bit reals data[points]:

******************************************************************************/

static Integer readRegriddedXDR( Stream* input, Profile* profile ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         profile, profile->variable == 0, profile->data == 0,
         profile->gridData == 0 );

  Integer result = 0;
  input->readString( input, profile->note, COUNT( profile->note ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( profile->note );

    if ( readTimestamp( input, profile->firstTimestamp ) ) {

      if ( readDimensions( input, 1, &profile->timesteps ) ) {
        profile->timestamps = NEW_ZERO( Integer, profile->timesteps );

        if ( profile->timestamps ) {
          Integer timestamp = fromUTCTimestamp( profile->firstTimestamp );
          Integer timestep = 0;

          for ( timestep = 0; timestep < profile->timesteps; ++timestep ) {
            profile->timestamps[ timestep ] = timestamp;
            incrementTimestamp( &timestamp );
          }

          profile->variables = 1;
          profile->variable = NEW_ZERO( Name, 2 );

          if ( profile->variable ) {
            profile->units = profile->variable + profile->variables;

            if ( readVariablesAndUnits( input, profile->variables,
                                        profile->variable, profile->units)) {
              char line[ 256 ] = "";
              int count = 8;
              memset( line, 0, sizeof line );
              input->readString( input, line, sizeof line / sizeof *line - 1 );

              if ( strcmp( line,
                           "# MSB 64-bit integers points[timesteps] and\n" )) {
                count += 4; /* Skip 4 line projection/grid. */
              }

              if ( skipInputLines( input, count - 1 ) ) {
                result = readRegriddedXDRData( input, profile );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Profile data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidProfile( profile ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readRegriddedXDRData - Read regridded binary array data from input.
INPUTS:  Stream* input       Input stream read from and parse.
OUTPUTS: Profile* profile  Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Input data format is:

# MSB 64-bit integers points[timesteps] and
# IEEE-754 64-bit reals longitudes[points] and
# IEEE-754 64-bit reals latitudes[points] and
# IEEE-754 64-bit reals elevations[points] and
# MSB 64-bit integers columns[points] and
# MSB 64-bit integers rows[points] and
# MSB 64-bit integers layers[points] and
# IEEE-754 64-bit reals data[points]:

******************************************************************************/

static Integer readRegriddedXDRData( Stream* input, Profile* profile ) {

  PRE08( input, input->ok( input ), input->isReadable( input ),
         profile, profile->timesteps > 0,
         profile->variables == 1,
         profile->profiles == 0, profile->data == 0 );

  Integer result = 0;
  profile->outputPoints = NEW_ZERO( Integer, profile->timesteps );

  if ( profile->outputPoints ) {
    input->read64BitIntegers( input, profile->outputPoints,
                              profile->timesteps );

    if ( input->ok( input ) ) {
      const Integer count = profile->totalRegriddedPoints =
        sumI( profile->outputPoints, profile->timesteps );

      if ( count > 0 ) {
        profile->gridLongitudes = NEW_ZERO( Real, count * 7 );

        if ( profile->outputPoints ) {
          profile->gridLatitudes = profile->gridLongitudes + count;
          profile->gridElevations = profile->gridLatitudes + count;
          profile->columns = (Integer*) profile->gridElevations + count;
          profile->rows = profile->columns + count;
          profile->layers = profile->rows + count;
          profile->gridData = (Real*) profile->layers + count;

          input->read64BitReals( input, profile->gridLongitudes, count * 3 );

          if ( input->ok( input ) ) {
            input->read64BitIntegers( input, profile->columns, count * 3 );

            if ( input->ok( input ) ) {
              input->read64BitReals( input, profile->gridData, count );

              if ( input->ok( input ) ) {
                result = isValidProfile( profile );
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid Profile data." );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidProfile( profile ) ) );
  return result;
}



/******************************************************************************
PURPOSE: compareRegriddedXDR - Compare Regridded data with CMAQ data.
INPUTS:  const Parameters* parameters  CMAQ data to compare to.
OUTPUTS: Profile* profile            Updated profile->data.
RETURNS: Integer 1 if comparable, else 0 and failureMessage() called.
******************************************************************************/

static Integer compareRegriddedXDR( const Parameters* parameters,
                                    Profile* profile ) {

  PRE05( parameters, isValidParameters( parameters ),
         parameters->compareFunction,
         profile, isValidProfile( profile ) );

  Integer result = 0;
  
  if ( ! AND2( ! strcmp( parameters->timestamp, profile->firstTimestamp ),
               parameters->timesteps == profile->timesteps ) ) {
    failureMessage( "Mismatched time steps (%s %lld)"
                    " for comparison to CMAQ data (%s %lld).",
                    profile->firstTimestamp, profile->timesteps,
                    parameters->timestamp, parameters->timesteps );
  } else {
    Real* const profileData             = profile->gridData;
    const Integer* const profileLayers  = profile->layers;
    const Integer* const profileRows    = profile->rows;
    const Integer* const profileColumns = profile->columns;
    const Integer* const profilePoints  = profile->outputPoints;
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
    Integer profileIndex = 0;

    DEBUG( fprintf( stderr,
                    "grid L,R,C ranges: %lld..%lld %lld..%lld %lld..%lld\n",
                    firstLayer, lastLayer, firstColumn, lastColumn,
                   firstRow, lastRow ); )

    DEBUG( fprintf( stderr, "timesteps = %lld\n", timesteps ); )

    for ( timestep = 0; timestep < timesteps; ++timestep ) {
      const Integer points = profilePoints[ timestep ];
      const Integer timestepOffset = timestep * layersTimesRowsTimesColumns;
      Integer point = 0;

      DEBUG( fprintf( stderr, "timestep = %lld, points = %lld\n",
                      timestep, points ); )

      for ( point = 0; point < points; ++point, ++profileIndex ) {
        const Integer profileLayer  = profileLayers[  profileIndex ];
        const Integer profileRow    = profileRows[    profileIndex ];
        const Integer profileColumn = profileColumns[ profileIndex ];
        DEBUG( fprintf( stderr, "  @[ %lld, %lld, %lld ]: ",
                        profileLayer, profileRow, profileColumn ); )

        if ( AND3( IN_RANGE( profileLayer, firstLayer, lastLayer ),
                   IN_RANGE( profileRow, firstRow, lastRow ),
                   IN_RANGE( profileColumn, firstColumn, lastColumn ) ) ) {
          const Integer profileLayer0  = profileLayer  - firstLayer;
          const Integer profileRow0    = profileRow    - firstRow;
          const Integer profileColumn0 = profileColumn - firstColumn;
          const Integer dataIndex =
            timestepOffset + profileLayer0 * rowsTimesColumns +
            profileRow0 * columns + profileColumn0;
          const Real profileDatum = profileData[ profileIndex ];
          const Real cmaqDatum = cmaqData[ dataIndex ];
          const Real comparedDatum = comparer( profileDatum, cmaqDatum );
          CHECK4( IN_RANGE( profileLayer0, 0, layers - 1 ),
                 IN_RANGE( profileRow0, 0, rows - 1 ),
                 IN_RANGE( profileColumn0, 0, columns - 1 ),
                 IN_RANGE( dataIndex, 0,
                          timesteps * layers * rows * columns - 1 ) );
          profileData[ profileIndex ] = comparedDatum;
          result = 1;
          DEBUG( fprintf( stderr, "f(%lf, %lf) -> %lf\n",
                          profileDatum, cmaqDatum, comparedDatum ); )
        } else {
          profileData[ profileIndex ] = -9999.0;
          DEBUG( fprintf( stderr, "-9999\n" ); )
        }
      }
    }
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "No points in output." );
  }

  POST02( IS_BOOL( result ), isValidProfile( profile ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format output.
INPUTS:  Profile* profile  Structure to write.
         const Parameters*   Parameters input stream to read from.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeASCII( Profile* profile, const Parameters* parameters ) {

  PRE03( isValidProfile( profile ),
         isValidParameters( parameters ),
         parameters->input->ok( parameters->input ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    writeASCIIHeader( profile, output );

    if ( output->ok( output ) ) {
      result = writeASCIIData( profile, output );
    }

    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeASCIIHeader - Write ASCII-format header line.
INPUTS:  const Profile* profile  Structure to write.
         Stream* output          Stream to write to.
******************************************************************************/

static void writeASCIIHeader( const Profile* profile, Stream* output ) {

  PRE03( isValidProfile( profile ), output, output->isWritable( output ) );

  const char* const headerStart =
    "timestamp(UTC)\tid(-)\tlongitude(deg)\tlatitude(deg)\televation(m)";
  const char* const headerFormat = "\t%s(%s)";

  /* Write header row: */

  output->writeString( output, headerStart );

  if ( output->ok( output ) ) {
    const Integer variables = profile->variables;
    Integer variable = IMPLICIT_VARIABLES;

    do {
      output->writeString( output, headerFormat,
                           profile->variable[ variable ],
                           profile->units[ variable ] );

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
INPUTS:  const Profile* profile  Structure to write.
OUTPUTS: Stream* output            Stream to write to.
RETURNS: Integer 1 if successful, else 0 and failureMesage is called.
******************************************************************************/

static Integer writeASCIIData( Profile* profile, Stream* output ) {

  PRE03( isValidProfile( profile ), output, output->isWritable( output ) );

  Integer result = 0;
  const char* const dataFormat = "\t%28.6"REAL_F_FORMAT;
  const Integer variables = profile->variables;
  const Integer profiles  = profile->profiles;
  const Real* profileData = profile->data;
  Integer p = 0;
  Integer theProfile = 0;

  /* Write data rows: */

  for ( theProfile = 0; theProfile < profiles; ++theProfile ) {
    const Integer profilePoints = profile->points[ theProfile ];
    Integer point = 0;

    for ( point = 0; point < profilePoints; ++point ) {
      Integer variable = 0;
      const Integer offset = p + point;
      const Integer timestamp = (Integer) profileData[ offset ];
      UTCTimestamp timestampString = "";
      CHECK( isValidYYYYMMDDHHMMSS( timestamp ) );
      toUTCTimestamp2( timestamp, timestampString );
      output->writeString( output, timestampString ); /* Begin row. */

      if ( output->ok( output ) ) {
        output->writeString( output, "\t%10lld", (Integer)
                             profileData[ offset + DATA_ID * profilePoints ] );

        for ( variable = 2;
              AND2( output->ok( output ), variable < variables );
              ++variable ) {
          const Real datum = profileData[ offset + variable * profilePoints ];
          output->writeString( output, dataFormat, datum );
        }

        if ( output->ok( output ) ) {
          output->writeString( output, "\n" ); /* End row. */
        }

        if ( ! output->ok( output ) ) {
          point = profilePoints;
          theProfile = profiles;
        }
      }
    }

    p += variables * profilePoints;
  }

  CHECK( p == profile->variables * profile->totalPoints );
  result = output->ok( output );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDS - Write COARDS-format data.
INPUTS:  Profile* profile              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDS( Profile* profile, const Parameters* parameters) {

  PRE02( isValidProfile( profile ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
      profile->variables * profile->totalPoints * 4 + /*variables(points).*/
      profile->totalPoints * 3 * 4 + 2000;
      /* yyyyddd(points),hhmmss(points),time(points) + header/extra. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file = createNetCDFFile(parameters->netcdfFileName, create64BitFile);

  if ( file != -1 ) {

    if ( writeCOARDSHeader( file, profile ) ) {
      result = writeCOARDSData( file, profile );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSHeader - Write header to file.
INPUTS:  Integer file            NetCDF file to write to.
         const Profile* profile  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSHeader( Integer file, const Profile* profile ) {

  PRE02( file != -1, isValidProfile( profile ) );

  Integer result = 0;
  const char* const name = "points";
  Integer dimensionId = -1;

  if ( createDimensions( file, 1, &name, &profile->totalPoints,
                         &dimensionId ) ) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {

        if ( createVariable( file, profile->variable[ DATA_ID ], "none",
                             NC_INT, 0, 1, &dimensionId ) != -1 ) {
          const Integer variables = profile->variables;
          Integer index = DATA_ELEVATION; /* Only write elevation and data. */

          do {
            const char* const units =
              strcmp( profile->units[ index ], "-"   ) == 0 ? "none" :
              strcmp( profile->units[ index ], "deg" ) == 0 ? "degrees" :
              profile->units[ index ];
            const Integer ok =
              createVariable( file, profile->variable[ index ], units,
                              NC_FLOAT, 1, 1, &dimensionId ) != -1;

            if ( ! ok ) {
              index = variables;
            }

            ++index;
          } while ( index < variables );

          if ( index == variables ) {

            if ( writeExtraAttributes( file, (const Real (*)[2]) profile->bounds,
                                       dimensionId ) ) {
              UTCTimestamp timestamp = "";
              Line history = "";
              appendToLine( history, profile->note );
              appendToLine( history, ",XDRConvert" );
              toUTCTimestamp2( profile->data[ 0 ], timestamp );
              result =
                writeStandardContents( file, history, timestamp,
                                       dimensionId, profile->totalPoints, 0 );
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
         Profile* profile  Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSData( Integer file, Profile* profile ) {

  PRE02( file != -1, isValidProfile( profile ) );

  Integer result = 0;
  const Integer variables   = profile->variables;
  const Integer totalPoints = profile->totalPoints;
  const Integer timesteps   = totalPoints;
  const Integer dataCount   = variables * totalPoints;

  /* Generate date-time data in profile->timestamps[]: */

  CHECK2( profile->copyData == 0, profile->timestamps == 0 );
  profile->copyData = NEW_ZERO( Real, dataCount );
  profile->timestamps = 
    profile->copyData ? NEW_ZERO( Integer, 3 * timesteps ) :0;

  if ( profile->timestamps ) {
    const Integer yyyydddhhmmStart =
      fromUTCTimestamp( profile->firstTimestamp );
    Integer* yyyyddd = profile->timestamps;
    Integer* hhmmss  = yyyyddd + timesteps;
    Real* time       = (Real*) hhmmss + timesteps;
    Integer variable = 0;
    Real* input = profile->data;
    Real* output = profile->copyData;
    Integer timestep = 0;
    Integer p = 0;
    Integer profileIndex = 0;

    /* Copy data to timestamp arrays and non-interleaved copyData: */

    for ( profileIndex = 0; profileIndex < profile->profiles; ++profileIndex) {
      const Integer profilePoints = profile->points[ profileIndex ];

      for ( variable = 0; variable < variables; ++variable ) {
        const Integer outputOffset = variable * totalPoints + p;
        Integer profilePoint = 0;

        for ( profilePoint = 0; profilePoint < profilePoints;
              ++profilePoint, ++input ) {
          const Real inputValue = *input;
          output[ outputOffset + profilePoint ] = inputValue;

          if ( variable == DATA_ID ) {
            Integer* const integerOutput =
              (Integer*) output + outputOffset + profilePoint;
            const Integer integerInputValue = (Integer) inputValue;
            *integerOutput = integerInputValue;
          }

          if ( variable == DATA_TIMESTAMP ) {
            const Integer yyyymmddhhmmss = (Integer) inputValue;
            const Integer yyyymmdd0 = yyyymmddhhmmss / 1000000;
            const Integer hhmmss0   = yyyymmddhhmmss % 1000000;
            const Integer yyyyddd0  = convertYearMonthDay( yyyymmdd0 );
            const Integer yyyydddhhmmss00 = yyyyddd0 * 1000000 + hhmmss0;
            const Integer yyyydddhhmmNow = yyyydddhhmmss00 / 100;
            const Real fractionalTime =
              fractionalHours( yyyydddhhmmStart, yyyydddhhmmNow );
            CHECK( IN_RANGE( timestep, 0, totalPoints - 1 ) );
            yyyyddd[ timestep ] = yyyyddd0;
            hhmmss[  timestep ] = hhmmss0;
            time[    timestep ] = fractionalTime;
            ++timestep;
          }
        }
      }

      p += profilePoints;
    }

    CHECK( timestep == totalPoints );

    /* Write and ruin copyData: */

    for ( variable = DATA_LONGITUDE; variable < variables; ++variable ) {
      const char* const variableName = profile->variable[ variable ];

      if ( ! writeAllData( file, variableName, totalPoints, 1, 1, 1,
                           profile->copyData + variable * totalPoints ) ) {
        variable = variables;
      }
    }

    /* Write and ruin yyyyddd: */

    result = variable == variables;
    
    result =
      AND2( result,
            writeAllIntData( file, "yyyyddd", totalPoints, 1, 1, 1, yyyyddd ));
      
    result =
      AND2( result,
            writeAllIntData( file, "hhmmss", totalPoints, 1, 1, 1, hhmmss ) );

    result =
      AND2( result, writeAllData( file, "time", totalPoints, 1, 1, 1, time ) );

    result =
      AND2( result,
            writeAllIntData( file, "id", totalPoints, 1, 1, 1,
                      (Integer*) (profile->copyData + DATA_ID * totalPoints)));
    
  }

  FREE( profile->copyData );
  FREE( profile->timestamps );

  POST02( IS_BOOL( result ), isValidProfile( profile ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRegriddedXDR - Write regridded XDR-format data.
INPUTS:  Profile* profile  Structure to write.
         const Parameters*
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
 
NOTES:   Output data format is:

REGRIDDED-Profile 1.0
http://www.esrl.noaa.gov/gmd/grad/neubrew/,NEUBrewSubset,XDRConvert
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
# IEEE-754 64-bit reals longitudes[points] and
# IEEE-754 64-bit reals latitudes[points] and
# IEEE-754 64-bit reals elevations[points] and
# MSB 64-bit integers columns[points] and
# MSB 64-bit integers rows[points] and
# MSB 64-bit integers layers[points] and
# IEEE-754 64-bit reals data[points]:

******************************************************************************/

static Integer writeRegriddedXDR( Profile* profile,
                                  const Parameters* parameters ) {

  PRE02( isValidProfile( profile ), isValidParameters( parameters ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer timesteps = profile->timesteps;
    const Integer points    = profile->totalRegriddedPoints;
    const Integer variableIndex =
      profile->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
    Name variable = "";
    aggregateName( profile->variable[ variableIndex ], hoursPerTimestep,
                   variable );

    output->writeString( output,
                         "REGRIDDED-Profile 2.0\n"
                         "%s,XDRConvert\n"
                         "%s\n"
                         "# timesteps\n%"INTEGER_FORMAT"\n"
                         "# Variable name:\n%s\n"
                         "# Variable units:\n%s\n",
                         profile->note,
                         profile->firstTimestamp,
                         timesteps,
                         variable,
                         profile->units[ variableIndex ] );

    writeProjectionAndGrid( parameters->grid, output );

    output->writeString( output,
                        "# MSB 64-bit integers points[timesteps] and\n"
                        "# IEEE-754 64-bit reals longitudes[points] and\n"
                        "# IEEE-754 64-bit reals latitudes[points] and\n"
                        "# IEEE-754 64-bit reals elevations[points] and\n"
                        "# MSB 64-bit integers columns[points] and\n"
                        "# MSB 64-bit integers rows[points] and\n"
                        "# MSB 64-bit integers layers[points] and\n"
                        "# IEEE-754 64-bit reals data[points]:\n" );

    if ( output->ok( output ) ) {
      output->write64BitIntegers( output, profile->outputPoints, timesteps );

      if ( output->ok( output ) ) {
        output->write64BitReals( output, profile->gridLongitudes, points );

        if ( output->ok( output ) ) {
          output->write64BitReals( output, profile->gridLatitudes, points );

          if ( output->ok( output ) ) {
            output->write64BitReals( output, profile->gridElevations, points);

            if ( output->ok( output ) ) {
              output->write64BitIntegers( output, profile->columns, points );

              if ( output->ok( output ) ) {
                output->write64BitIntegers( output, profile->rows, points );

                if ( output->ok( output ) ) {
                  output->write64BitIntegers(output, profile->layers, points);

                  if ( output->ok( output ) ) {
                    output->write64BitReals( output, profile->gridData,
                                             points );
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
INPUTS:  Profile* profile  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedASCII( Profile* profile,
                                    const Parameters* parameters ) {

  PRE03( isValidProfile( profile ), profile->variables > 0,
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
      const Integer variableIndex =
        profile->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
      const Integer hoursPerTimestep =
        parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;
      Name variable = "";
      aggregateName( profile->variable[ variableIndex ], hoursPerTimestep,
                     variable );

      output->writeString( output, headerFormat,
                           variable, profile->units[ variableIndex ] );

      if ( output->ok( output ) ) {
        const Integer timesteps = profile->timesteps;
        const Real* longitudes  = profile->gridLongitudes;
        const Real* latitudes   = profile->gridLatitudes;
        const Real* elevations  = profile->gridElevations;
        const Integer* columns  = profile->columns;
        const Integer* rows     = profile->rows;
        const Integer* layers   = profile->layers;
        const Real* data        = profile->gridData;
        const Integer hoursPerTimestep =
          parameters->aggregationTimesteps ? parameters->aggregationTimesteps
          : 1;
        Integer timestep = 0;
        Integer yyyydddhh00 =
          ( fromUTCTimestamp( profile->firstTimestamp ) / 100 ) * 100;
        UTCTimestamp timestamp = "";

        /* Write data rows: */

        do {
          const Integer points = profile->outputPoints[ timestep ];
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

            output->writeString( output, dataFormat,
                                 timestamp, longitude, latitude, elevation,
                                 column, row, layer, value );

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
INPUTS:  Profile* profile              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDS( Profile* profile,
                                     const Parameters* parameters ) {

  PRE02( isValidProfile( profile ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    profile->totalRegriddedPoints * 7 * 4 + 10000;
    /* lon, lat, elv, col, row, lay, time, hdr. */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedCOARDSHeader( file, hoursPerTimestep, profile ) ) {
      result = writeRegriddedCOARDSData( file, profile, parameters );
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
         const Profile* profile    Structure to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSHeader( Integer file,
                                           Integer hoursPerTimestep,
                                           const Profile* profile ) {

  PRE03( file != -1, hoursPerTimestep > 0, isValidProfile( profile ) );

  Integer result = 0;
  const char* const dimensionName = "points";
  Integer dimensionId = -1;
  const Integer dimension = profile->totalRegriddedPoints;

  if ( createDimensions( file, 1, &dimensionName, &dimension, &dimensionId)) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createVariable( file, "column", "-",
                           NC_INT, 0, 1, &dimensionId ) != -1 ) {

        if ( createVariable( file, "row", "-",
                             NC_INT, 0, 1, &dimensionId ) != -1 ) {

          if ( createVariable( file, "layer", "-",
                               NC_INT, 0, 1, &dimensionId ) != -1 ) {

            if ( createLongitudeAndLatitude( file, 1, &dimensionId ) ) {

              if ( createVariable( file, "elevation", "-",
                                   NC_FLOAT, 0, 1, &dimensionId ) != -1 ) {
                const Integer variableIndex =
                  profile->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES
                  : 0;
                Name variable = "";
                aggregateName( profile->variable[ variableIndex ],
                               hoursPerTimestep, variable );

                if ( createVariable( file,
                                     variable,
                                     profile->units[ variableIndex ],
                                     NC_FLOAT, 1, 1, &dimensionId ) != -1 ) {

                  UTCTimestamp timestamp = "";
                  Line history = "";
                  appendToLine( history, profile->note );
                  appendToLine( history, ",XDRConvert" );
                  toUTCTimestamp( profile->timestamps[ 0 ], timestamp );

                  result = writeStandardContents( file, history,
                                                  timestamp,
                                                  dimensionId, 0, 0 );
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
INPUTS:  Integer file            NetCDF file to write to.
         const Profile* profile  Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedCOARDSData( Integer file,
                                         const Profile* profile,
                                         const Parameters* parameters ) {

  PRE03( file != -1, isValidProfile( profile ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer count = profile->totalRegriddedPoints;

  if ( writeAllIntData( file, "column", count, 1, 1, 1,
                        profile->columns ) ) {

    if ( writeAllIntData( file, "row", count, 1, 1, 1,
                          profile->rows ) ) {

      if ( writeAllIntData( file, "layer", count, 1, 1, 1,
                            profile->layers ) ) {

        if ( writeAllData( file, "longitude", count, 1, 1, 1,
                            profile->gridLongitudes ) ) {

          if ( writeAllData( file, "latitude", count, 1, 1, 1,
                              profile->gridLatitudes ) ) {

            if ( writeAllData( file, "elevation", count, 1, 1, 1,
                                profile->gridElevations ) ) {
              const Integer variableIndex =
                profile->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES
                : 0;
              const Integer hoursPerTimestep =
                parameters->aggregationTimesteps ?
                  parameters->aggregationTimesteps
                : 1;
              Name variable = "";
              aggregateName( profile->variable[ variableIndex ],
                             hoursPerTimestep, variable );

              if ( writeAllData( file, variable, count, 1, 1, 1,
                                 profile->gridData ) ) {
                timeData( profile->timesteps, hoursPerTimestep, count,
                          profile->outputPoints, profile->gridData );
                
                result = writeAllData( file, "time", count, 1, 1, 1,
                                       profile->gridData );
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
INPUTS:  Profile* profile              Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPI( Profile* profile,
                                    const Parameters* parameters ) {

  PRE02( isValidProfile( profile ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    profile->totalRegriddedPoints * 4 * 4 + 10000; /* lon,lat,elv,var, hdr */
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {
    const Integer hoursPerTimestep =
      parameters->aggregationTimesteps ? parameters->aggregationTimesteps : 1;

    if ( writeRegriddedIOAPIHeader( file, hoursPerTimestep,
                                    profile, parameters->grid ) ) {
      result = writeRegriddedIOAPIData( file, hoursPerTimestep, profile,
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
         const Profile* profile    Structure to write.
         const Grid* grid          Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIHeader( Integer file,
                                          Integer hoursPerTimestep,
                                          const Profile* profile,
                                          const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0,
         isValidProfile( profile ),
         grid, grid->invariant( grid ) );

  Integer result = 0;
  enum { VARIABLES = 4 }; /* LONGITUDE, LATITUDE, ELEVATION, profile. */
  Name variableNames[ VARIABLES ] = {
    "LONGITUDE", "LATITUDE", "ELEVATION", "profile"
  };
  Name variableUnits[ VARIABLES ] = { "deg", "deg", "m", "-" };
  const Integer layers = grid->layers( grid );
  const Integer firstTimestamp = fromUTCTimestamp( profile->firstTimestamp );
  const Integer variableIndex =
    profile->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
  Line history = "";
  appendToLine( history, profile->note );
  appendToLine( history, ",XDRConvert" );
  aggregateName( profile->variable[ variableIndex ], hoursPerTimestep,
                 variableNames[ VARIABLES - 1 ] );
  variableNames[ VARIABLES - 1 ][ 15 ] = '\0';
  strncpy( variableUnits[ VARIABLES - 1 ],
           profile->units[ variableIndex ], 16 );
  variableUnits[ VARIABLES - 1 ][ 16 ] = '\0';
  uppercase( variableNames[ VARIABLES - 1 ] );
  lowercase( variableUnits[ VARIABLES - 1 ] );

  result = writeM3IOHeader( file, profile->timesteps, hoursPerTimestep,
                            firstTimestamp,
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
         const Profile* profile    Structure to write.
         const Grid* grid          Grid info.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeRegriddedIOAPIData( Integer file,
                                        Integer hoursPerTimestep,
                                        const Profile* profile,
                                        const Grid* grid ) {

  PRE05( file != -1, hoursPerTimestep > 0, isValidProfile( profile ),
         grid, grid->invariant( grid ) );

  Integer result = 0;
  const Integer layers   = grid->layers( grid );
  const Integer rows     = grid->rows( grid );
  const Integer columns  = grid->columns( grid );
  const Integer cells    = layers * rows * columns;
  Real* expandedGridData = NEW_ZERO( Real, cells );

  if ( expandedGridData ) {
    const Integer timesteps = profile->timesteps;
    const Real scale = 1.0;
    Integer timestep = 0;
    Integer offset = 0;

    if ( writeM3IOGrid( grid, timesteps, layers, file ) ) {

      do {
        const Integer points = profile->outputPoints[ timestep ];
        Name variable = "";
        memset( variable, 0, sizeof variable );
        strncpy( variable, "ELEVATION", 16 );
        uppercase( variable );

        copyDataToGrid3( points,
                         profile->layers         + offset,
                         profile->rows           + offset,
                         profile->columns        + offset,
                         profile->gridElevations + offset,
                         scale, layers, rows, columns,
                         expandedGridData );

        if ( ! writeM3IOData( file, variable,
                              timestep, layers, rows, columns,
                              expandedGridData ) ) {
          timestep = timesteps;
        } else {
          const Integer variableIndex =
            profile->variables > IMPLICIT_VARIABLES ? IMPLICIT_VARIABLES : 0;
          Name variable = "";
          aggregateName( profile->variable[ variableIndex ], hoursPerTimestep,
                         variable );
          variable[ 15 ] = '\0';
          uppercase( variable );

          copyDataToGrid3( points,
                           profile->layers   + offset,
                           profile->rows     + offset,
                           profile->columns  + offset,
                           profile->gridData + offset,
                           scale, layers, rows, columns,
                           expandedGridData );

          if ( ! writeM3IOData( file, variable,
                                timestep, layers, rows, columns,
                                expandedGridData ) ) {
            timestep = timesteps;
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
PURPOSE: regridProfile - Regrid data.
INPUTS:  Integer method    E.g., AGGREGATE_MEAN.
         Grid* grid        Grid to project and aggregate points into.
         Profile* profile  Data to regrid.
OUTPUTS: Profile* profile  Regridded data.
******************************************************************************/

static void regridProfile( Integer method, Grid* grid, Profile* profile ) {

  PRE06( IS_VALID_AGGREGATE_METHOD( method ),
         grid,
         grid->invariant( grid ),
         isValidProfile( profile ),
         profile->totalRegriddedPoints == 0,
         profile->longitudes == 0 );

  const Integer variables = profile->variables;

  if ( variables == IMPLICIT_VARIABLES + 1 ) {
    const Integer timesteps = hoursInRange( profile->firstTimestamp,
                                            profile->lastTimestamp );
    const Integer inputVariables = 4; /* lon, lat, elv, dat. */
    const Integer inputSize      = profile->totalPoints; /*All points in 1hr*/
    const Integer inputDataSize  = inputVariables * inputSize;
    const Integer outputVariables = 7; /* glon,glat,gelv,col,row,lay,gdata. */
    const Integer outputSize      = inputSize; /* 1 point in each hr? */
    const Integer outputDataSize  = outputVariables * outputSize;
    const Integer dataSize  = inputDataSize + outputDataSize + timesteps * 2;
    profile->longitudes = NEW_ZERO( Real, dataSize );

    if ( profile->longitudes ) {
      Integer totalRegriddedPoints = 0;
      Integer timestep = 0;
      Integer yyyydddhh00 =
        ( fromUTCTimestamp( profile->firstTimestamp ) / 100 ) * 100;
      profile->latitudes      = profile->longitudes     + inputSize;
      profile->elevations     = profile->latitudes      + inputSize;
      profile->copyData       = profile->elevations     + inputSize;
      profile->gridLongitudes = profile->copyData       + inputSize;
      profile->gridLatitudes  = profile->gridLongitudes + outputSize;
      profile->gridElevations = profile->gridLatitudes  + outputSize;
      profile->gridData       = profile->gridElevations + outputSize;
      profile->columns        = (Integer*)
                                 profile->gridData      + outputSize;
      profile->rows           = profile->columns        + outputSize;
      profile->layers         = profile->rows           + outputSize;
      profile->outputPoints   = profile->layers         + outputSize;
      profile->timestamps     = profile->outputPoints   + timesteps;
      profile->timesteps      = timesteps;

      do {
        const Integer inputPoints =
          copyDataForTimestamp( yyyydddhh00, profile );
          profile->timestamps[ timestep ] = yyyydddhh00;

        if ( inputPoints ) {
          Integer outputPoints = 0;
          const Real minimumValidValue = -900.0;

          DEBUG( fprintf( stderr, "regridding %lld points for %lld...\n",
                          inputPoints, yyyydddhh00 ); )
          
          grid->regrid( grid, method, minimumValidValue, inputPoints, 1,
                        profile->longitudes, profile->latitudes,
                        profile->elevations, profile->copyData,
                        0, /* No input vector data. */
                        0, /* No notes. */
                        &outputPoints,
                        profile->columns        + totalRegriddedPoints,
                        profile->rows           + totalRegriddedPoints,
                        profile->layers         + totalRegriddedPoints,
                        profile->gridLongitudes + totalRegriddedPoints,
                        profile->gridLatitudes  + totalRegriddedPoints,
                        profile->gridElevations + totalRegriddedPoints,
                        profile->gridData       + totalRegriddedPoints,
                        0, /* No output vector data. */
                        0 /* No regriddedNotes. */ );

          profile->outputPoints[ timestep ] = outputPoints;
          totalRegriddedPoints += outputPoints;

          DEBUG( fprintf( stderr, "outputPoints = %lld,"
                          " totalRegriddedPoints = %lld\n",
                          outputPoints, totalRegriddedPoints ); )
        }

        incrementTimestamp( &yyyydddhh00 );
        ++timestep;
      } while ( timestep < timesteps );

      profile->totalRegriddedPoints = totalRegriddedPoints;
    }
  }

  POST02( profile->totalRegriddedPoints >= 0,
          IMPLIES( profile->totalRegriddedPoints > 0,
            AND10( IN_RANGE( minimumItemI( profile->outputPoints,
                                           profile->timesteps ),
                             0, profile->totalRegriddedPoints ),
                   IN_RANGE( maximumItemI( profile->outputPoints,
                                           profile->timesteps ),
                                           1, profile->totalRegriddedPoints ),
                   IN_RANGE( minimumItemI( profile->columns,
                                           profile->totalRegriddedPoints),
                                           1, grid->columns( grid ) ),
                   IN_RANGE( maximumItemI( profile->columns,
                                           profile->totalRegriddedPoints),
                                           1, grid->columns( grid ) ),
                   IN_RANGE( minimumItemI( profile->rows,
                                           profile->totalRegriddedPoints),
                                           1, grid->rows( grid ) ),
                   IN_RANGE( maximumItemI( profile->rows,
                                           profile->totalRegriddedPoints),
                                           1, grid->rows( grid ) ),
                   IN_RANGE( minimumItemI( profile->layers,
                                           profile->totalRegriddedPoints),
                                           1, grid->layers( grid ) ),
                   IN_RANGE( maximumItemI( profile->layers,
                                           profile->totalRegriddedPoints),
                                           1, grid->layers( grid ) ),
                   validLongitudesAndLatitudes( profile->totalRegriddedPoints,
                                                profile->gridLongitudes,
                                                profile->gridLatitudes ),
                   isNanFree( profile->gridData,
                              profile->totalRegriddedPoints ) ) ) );
}



/******************************************************************************
PURPOSE: copyDataForTimestamp - Copy data for given regrid timestamp.
INPUTS:  Integer yyyydddhh00  Timestamp to copy data for.
         Profile* profile     profile->data data to copy.
OUTPUTS: Profile* profile     profile->longitudes, latitudes, elevations,
                              copyData.
RETURNS: Integer number of points copied for the timestamp.
******************************************************************************/

static Integer copyDataForTimestamp(Integer yyyydddhh00, Profile* profile) {

  PRE02( isValidTimestamp( yyyydddhh00 ), isValidProfile( profile ) );
  const Integer profiles = profile->profiles;
  const Integer variables = profile->variables;
  const Real* data = profile->data;
  Real* longitudes = profile->longitudes;
  Real* latitudes  = profile->latitudes;
  Real* elevations = profile->elevations;
  Real* copyData   = profile->copyData;
  Integer p = 0;
  Integer theProfile = 0;
  Integer result = 0;

  for ( theProfile = 0; theProfile < profiles; ++theProfile ) {
    const Integer profilePoints = profile->points[ theProfile ];
    Integer point = 0;

    for ( point = 0; point < profilePoints; ++point ) {
      const Integer offset = p + point;
      const Integer pointTimestamp = (Integer) data[ offset ];
      Integer timestamp = 0;
      UTCTimestamp timestampString = "";
      CHECK( isValidYYYYMMDDHHMMSS( pointTimestamp ) );
      toUTCTimestamp2( pointTimestamp, timestampString );
      timestamp = fromUTCTimestamp( timestampString );
      timestamp /= 100;
      timestamp *= 100;

      if ( timestamp == yyyydddhh00 ) {
        const Real longitude = data[ offset + DATA_LONGITUDE * profilePoints ];
        const Real latitude  = data[ offset + DATA_LATITUDE  * profilePoints ];
        const Real elevation = data[ offset + DATA_ELEVATION * profilePoints ];
        const Real ozone     = data[ offset + DATA_OZONE     * profilePoints ];
        CHECK4( isValidLongitude( longitude ), isValidLatitude( latitude ),
                IN_RANGE( elevation, -500.0, 1e6 ), ozone >= 0.0 );
        *longitudes++ = longitude;
        *latitudes++  = latitude;
        *elevations++ = elevation;
        *copyData++   = ozone;
        ++result;
      }
    }

    p += variables * profilePoints;
  }

  CHECK( p == profile->variables * profile->totalPoints );
  POST02( result >= 0,
         IMPLIES( result > 0,
                  AND3( validLongitudesAndLatitudes( result,
                                                     profile->longitudes,
                                                     profile->latitudes ),
                        isNanFree( profile->elevations, result ),
                        isNanFree( profile->copyData, result ) ) ) );
  return result;
}



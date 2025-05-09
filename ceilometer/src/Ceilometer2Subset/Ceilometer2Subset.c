
/******************************************************************************
PURPOSE: Ceilometer2Subset.c - Read a set of a Ceilometer2 files, subset the
         data to a bounds (longitude-latitude rectangle) and write it to
         stdout as XDR (IEEE-754) format binary.

NOTES:   http://www.vaisala.com/en/products/ceilometers/Pages/CL51.aspx

HISTORY: 2022-07-12 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>     /* For printf(). */
#include <string.h>    /* For strlen(). */
#include <ctype.h>     /* For isdigit(). */

#include <Utilities.h> /* For PRE0*(), NEW_ZERO(), Stream, VoidList, etc. */
#include "ReadData.h"  /* For MAXIMUM_ELEVATION, readSubsetCeilometerData(). */

/*================================== TYPES ==================================*/

enum { VARIABLES = 6 }; /* timestamp, id, longitude, latitude, elevation, var*/

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;         /* File listing Ceilometer files to read. */
  const char* description;      /* User-supplied description. */
  const char* variable;         /* Data variable name. */
  Bounds      bounds;           /*bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]*/
  Real        elevationRange[ 2 ]; /* elevationRange[ MINIMUM, MAXIMUM ]. */
  Integer     firstTimestamp;   /* YYYYMMDDHHMMSS of subset. */
  Integer     lastTimestamp;    /* YYYYMMDDHHMMSS of subset. */
} Arguments;


static Integer isValidElevationRange( const Real elevationRange[ 2 ] ) {
  const Integer result =
  AND3( isNanFree( elevationRange, 2 ),
        IN_RANGE( elevationRange[ MINIMUM ],
                  MINIMUM_SURFACE_ELEVATION, MAXIMUM_ELEVATION ),
        IN_RANGE( elevationRange[ MAXIMUM ],
                  elevationRange[ MINIMUM ], MAXIMUM_ELEVATION ) );
  POST0( IS_BOOL( result ) );
  return result;
}


#ifndef NO_ASSERTIONS

/* Arguments invariant: */

static Integer isValidArguments( const Arguments* arguments ) {
  const Integer result =
    AND13( arguments,
           arguments->listFile,
           arguments->listFile[ 0 ],
           arguments->description,
           arguments->description[ 0 ],
           arguments->variable,
           arguments->variable[ 0 ],
           isValidBounds( arguments->bounds ),
           isNanFree( arguments->elevationRange, 2 ),
           isValidElevationRange( arguments->elevationRange ),
           isValidYYYYMMDDHHMMSS( arguments->firstTimestamp ),
           isValidYYYYMMDDHHMMSS( arguments->lastTimestamp ),
           arguments->firstTimestamp <= arguments->lastTimestamp );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Profile: Subsetted result of reading a Ceilometer LIDAR data file. */

typedef struct {
  Integer timestamp; /* YYYYDDDHHMMSS. */
  Integer points;    /* # of filtered data points (timesteps x levels).*/
  Real* data;        /* data[ points ]. */
  Real* elevations;  /* elevations[points]. Meters above mean sea level*/
  Real* timestamps;  /* timestamps[ points ]. yyyymmddhhmmss. */
  Real longitude;    /* Of station. */
  Real latitude;     /* Of station. */
  Real elevation;    /* Of station. Meters above mean sea level. */
  Real id;           /* Constructed from hashing part of data file name.*/
  char units[ 80 ];  /* Of data variable. */
  char notes[ 80 ];  /* Site name, location, id parsed from file name. */
} Profile;


/* Profile destructor: */

static void deallocateProfile( void* vprofile ) {
  PRE0( vprofile );
  Profile* profile = (Profile*) vprofile;
  FREE( profile->data );
  FREE( profile->elevations );
  FREE( profile->timestamps );
  ZERO_OBJECT( profile );
  POST0( vprofile );
}

#ifndef NO_ASSERTIONS

/* Profile invariant: */

static Integer isValidProfile( const Profile* profile ) {
  Integer result =
    AND13( profile,
           isValidYYYYMMDDHHMMSS( profile->timestamp ),
           IN_RANGE( profile->points, 1, INTEGER_MAX / VARIABLES ),
           profile->data,
           profile->elevations,
           profile->timestamps,
           isNanFree( profile->data, profile->points ),
           isValidLongitude( profile->longitude ),
           isValidLatitude( profile->latitude ),
           IN_RANGE( profile->elevation,
                     MINIMUM_SURFACE_ELEVATION, MAXIMUM_ELEVATION ),
           profile->id > 0,
           profile->units[ 0 ],
           profile->notes[ 0 ] );

  if ( result ) {
    const Integer points         = profile->points;
    const Real* const data       = profile->data;
    const Real* const timestamps = profile->timestamps;
    const Real* const elevations = profile->elevations;
    const Real minimum =
      profile->units[ 0 ] == 'm' ? MINIMUM_SURFACE_ELEVATION : 0.0;
    const Real maximum =
      profile->units[ 0 ] == 'm' ? MAXIMUM_ELEVATION
      : MAXIMUM_VALID_CEILOMETER_BACKSCATTER_VALUE;
    Integer point = 0;

    for ( point = 0; AND2( result, point < points ); ++point ) {
      const Integer elevation = elevations[ point ];
      const Integer timestamp = (Integer) timestamps[ point ];
      const Real value        = data[ point ];
      result = isValidYYYYMMDDHHMMSS( timestamp );
      result = AND2( result, IN_RANGE( value, minimum, maximum ) );
      result = AND2( result, IN_RANGE( elevation, MINIMUM_SURFACE_ELEVATION,
                                       MAXIMUM_ELEVATION ) );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */

#ifdef DEBUGGING

static void printProfile( const Profile* profile ) {
  PRE0( profile );

  fprintf( stderr,
            "\nprofile: id = %lld (%s) @ (%lg, %lg, %lf), %lld, points = %lld",
            (Integer) profile->id, profile->notes,
            profile->longitude, profile->latitude,
            profile->elevation, profile->timestamp, profile->points );

  if ( AND3( profile->data, profile->elevations, profile->timestamps ) ) {
    const Real* const data = profile->data;
    const Real* const elevations = profile->elevations;
    const Real* const timestamps = profile->timestamps;
    const Integer points = profile->points;

    fprintf( stderr, "\n  data:" );

    if ( points >= 3 ) {
      fprintf( stderr, " %lg %lg %lg", data[ 0 ] , data[ 1 ] , data[ 2 ] );
      fprintf( stderr, " ... %lg %lg %lg",
               data[ points - 3 ], data[ points - 2 ] , data[ points - 1 ] );
    }

    fprintf( stderr, "\n  elevations:" );

    if ( points >= 3 ) {
      fprintf( stderr, " %lg %lg %lg",
               elevations[ 0 ] , elevations[ 1 ] , elevations[ 2 ] );
      fprintf( stderr, " ... %lg %lg %lg",
               elevations[ points - 3 ], elevations[ points - 2 ] ,
               elevations[ points - 1 ] );
    }

    fprintf( stderr, "\n  timestamps:" );

    if ( points >= 3 ) {
      fprintf( stderr, " %lld %lld %lld",
               (Integer) timestamps[ 0 ],
               (Integer) timestamps[ 1 ],
               (Integer) timestamps[ 2 ] );
      fprintf( stderr, " ... %lld %lld %lld",
               (Integer) timestamps[ points - 3 ] ,
               (Integer) timestamps[ points - 2 ] ,
               (Integer) timestamps[ points - 1 ] );
    }
  }

  fprintf( stderr, "\n\n" );
}

#endif /* DEBUGGING */


/* Data type: */

typedef struct {
  Arguments arguments; /* User-supplied (command-line) arguments. */
  VoidList* profiles;  /* List of subsetted profiles. */
  Integer   ok;        /* Did last command succeed? */
} Data;

/* Data destructor: */

static void deallocateData( Data* data ) {
  PRE0( data );
  FREE_OBJECT( data->profiles ); /* Calls deallocateProfile(). */
  ZERO_OBJECT( data );
  POST0( data );
}

#ifndef NO_ASSERTIONS

/* Data invariant: */

static Integer isValidData( const Data* data ) {
  Integer result =
    AND5( data,
          isValidArguments( &data->arguments ),
          data->profiles,
          data->profiles->count( data->profiles ),
          IS_BOOL( data->ok ) );

  if ( result ) {
    const VoidList* const profiles = data->profiles;
    const Integer profileCount = profiles->count( profiles );
    Integer index = 0;

    do {
      const Profile* const profile = profiles->item( profiles, index );
      result = isValidProfile( profile );
      ++index;
    } while ( AND2( result, index < profileCount ) );
  }

  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */


/*========================== FORWARD DECLARATIONS ===========================*/

static void printUsage( const char* programName );

static Integer parseArguments( Integer argc, char* argv[],
                               Arguments* arguments );

static void initializeArguments( Arguments* arguments );

static Integer parseOptionalArguments( Integer argc, char* argv[],
                                       Integer* arg,
                                       Arguments* arguments );

static Integer parsedElevationRange( Integer argc, char* argv[], Integer* arg,
                                     Real range[ 2 ] );

static void readData( Data* data );

static Profile* readCeilometerFile( const char* fileName,
                                    const char* variableName,
                                    const Integer firstTimestamp,
                                    const Integer lastTimestamp,
                                    const Bounds bounds,
                                    const Real elevationRange[ 2 ] );

static Integer instrumentId( const char* fileName );

static void appendOrFreeProfile( Data* data, Profile* profile );

static Integer parseFileTimestamp( const char* fileName );

static void parseFileNotes( const char* fileName, const Integer id,
                            char notes[ 80 ] );

static void writeData( Data* data );

static void writeHeader( Data* data, Stream* output );

static void writeXDR( Data* data, Stream* output );

static void writeProfileNotes( const VoidList* profiles, Stream* output );

static void writeProfilePoints( const VoidList* profiles, Stream* output );

static Integer writeProfileData( const VoidList* profiles, Stream* output );

static Integer maximumProfilePoints( const VoidList* profiles );

static void copyValue( const Integer count, Real array[], const Real value );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Read a subset of a list Ceilometer files and write it to stdout
         in XDR format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Integer ok = 0;

  if ( ! isValidArgs( argc, (const char**) argv ) ) {
    failureMessage( "Invalid command-line arguments." );
    printUsage( argv[ 0 ] );
  } else {
    Data data;
    ZERO_OBJECT( &data );
    checkForTest( &argc, argv ); /* Check for and remove any -test arguments.*/
    data.ok = parseArguments( argc, argv, &data.arguments );

    if ( data.ok ) {
      readData( &data ); /* From list file named in arguments. */

      if ( data.ok ) {
        writeData( &data ); /* To stdout. */
      }
    }

    ok = AND2( data.ok, failureCount() == 0 );
    deallocateData( &data );
  }

  POST0( IS_BOOL( ok ) );
  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* programName  Name of program.
******************************************************************************/

static void printUsage( const char* programName ) {
  PRE0( programName );
  fprintf( stderr,
           "\n\n%s - Read a set of Ceilometer files and extract profile\n",
           programName );
  fprintf( stderr, "data for selected variable" );
  fprintf( stderr, " subsetted by time range, lon-lat rectangle, and "
           "elevation range.\n" );
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "  -files <listFile> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -time <yyyymmddhhmmss> <yyyymmddhhmmss> \\\n" );
  fprintf( stderr, "  -variable name \\\n" );
  fprintf( stderr, "  [ -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, "  [ -elevation <maximum_elevation> <maximum_elevation> ] "
                   "\\\n\n\n" );
  fprintf( stderr, "Note: times are in UTC (GMT), " );
  fprintf( stderr, "and elevations are in meters above mean sea level.\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files /data/files.txt \\\n" );
  fprintf( stderr, "-desc "
           "http://www.vaisala.com/en/products/ceilometers/Pages/CL51.aspx,"
           "CeilometerSubset \\\n" );
  fprintf( stderr, "-time 20211007000000 20211007235959 \\\n" );
  fprintf( stderr, "-variable attenuated_backscatter \\\n" );
  fprintf( stderr, "-domain -130 25 -60 50 \\\n" );
  fprintf( stderr, "-elevation 0 20000 > subset.xdr\n\n" );
  fprintf( stderr,"Subset of data up to 20,000m for July 17, 2014 over USA\n");
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays\n" );
  fprintf( stderr, "For example:\n" );
  fprintf( stderr, "Profile 2.0\n" );
  fprintf( stderr,
           "http://www.vaisala.com/en/products/ceilometers/Pages/CL51.aspx,"
           "CeilometerSubset\n");
  fprintf( stderr, "2014-07-17T00:00:00-0000 2014-07-17T23:59:59-0000\n" );
  fprintf( stderr, "# Subset domain: <min_lon> <min_lat> <max_lon> <max_lat>");
  fprintf( stderr, ":\n");
  fprintf( stderr, "-130 25 -60 50\n" );
  fprintf( stderr, "# Dimensions: variables profiles:\n" );
  fprintf( stderr, "6 1000000\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "timestamp id longitude latitude elevation "
                   "attenuated_backscatter\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss - deg deg m molecules/m3\n" );
  fprintf( stderr, "# char notes[profiles][80] and\n" );
  fprintf( stderr, "# MSB 64-bit integers points[profiles] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals" );
  fprintf( stderr, " data_1[variables][points_1] ..." );
  fprintf( stderr, " data_P[variables][points_1]:\n" );
  fprintf( stderr, "<binary data arrays here>\n\n\n" );
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: Integer 1 if successful, else 0 and failureMessage() then printUsage()
         are called.
******************************************************************************/

static Integer parseArguments( Integer argc, char* argv[],
                               Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  Integer result = 0;
  ZERO_OBJECT( arguments );
  initializeArguments( arguments );

  if ( ! IN_RANGE( argc, 10, 18 ) ) {
    failureMessage( "Invalid/insufficient command line arguments.");
  } else {
    Integer arg = 1; /* Argument to parse next, updated by parseArgument2(). */
    arguments->listFile = parseArgument2( argc, argv, "-files", &arg );

    if ( arguments->listFile ) {
      arguments->description = parseArgument2( argc, argv, "-desc", &arg );

      if ( arguments->description ) {

        if ( AND2( ! strcmp( argv[ arg ], "-time" ),
                   parseTimeRange( argv[ arg + 1 ], argv[ arg + 2 ],
                             &arguments->firstTimestamp,
                             &arguments->lastTimestamp ) ) ) {
          arg += 3;
          {
            const char* const variableName =
              parseArgument2( argc, argv, "-variable", &arg );
            result = variableName != 0;

            if ( ! result ) {
              failureMessage( "Invalid variable name spcified '%s'\n",
                              variableName );
            } else {
              arguments->variable = variableName;
              result =
                AND2( result,
                      parseOptionalArguments( argc, argv, &arg, arguments ) );
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    ZERO_OBJECT( arguments );
    printUsage( argv[ 0 ] );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments )));
  return result;
}



/******************************************************************************
PURPOSE: initializeArguments - Initialize arguments.
INPUTS:  Arguments* arguments  Arguments to initialize.
OUTPUTS: Arguments* arguments  Arguments initialized.
******************************************************************************/

static void initializeArguments( Arguments* arguments ) {
  PRE0( arguments );
  ZERO_OBJECT( arguments );
  arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
  arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
  arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
  arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
  arguments->elevationRange[ MINIMUM ] = MINIMUM_SURFACE_ELEVATION;
  arguments->elevationRange[ MAXIMUM ] = MAXIMUM_ELEVATION;
}



/******************************************************************************
PURPOSE: parseOptionalArguments - Parse optional command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
         Integer* arg          Index of current argument to parse.
OUTPUTS: Integer* arg          Index of next argument to parse.
         Arguments* arguments  Updated arguments.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer parseOptionalArguments( Integer argc, char* argv[],
                                       Integer* arg,
                                       Arguments* arguments ) {

  PRE04( isValidArgs( argc, (const char**) argv ), arg, *arg > 0, arguments );

  Integer result = 1;
  Integer parsedBounds = 0;
  Integer parsedElevation = 0;

  while ( *arg < argc ) {

    if ( AND2( ! strcmp( argv[ *arg ], "-domain" ), ! parsedBounds )) {
      parsedBounds = 1;
      result = parseBounds( argc, argv, arg, arguments->bounds );
    } else if ( AND2( ! strcmp( argv[ *arg ], "-elevation" ),
                      ! parsedElevation ) ) {
      parsedElevation = 1;
      result = parsedElevationRange( argc, argv, arg,
                                     arguments->elevationRange );
    } else {
      failureMessage( "Invalid/redundant command-line argument: %s.",
                      argv[ *arg ] );
      result = 0;
      *arg = argc;
    }
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments )));
  return result;
}



/******************************************************************************
PURPOSE: parsedElevationRange - Parse command-line arguments for -elevation.
INPUTS:  Integer argc     Number of command-line arguments.
         char* argv[]     Command-line argument strings.
         Integer* arg     Index of argument to parse.
OUTPUTS: Integer* arg     Index of next argument to parse.
         Real range[ 2 ]  range[ MINIMUM, MAXIMUM ].
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer parsedElevationRange( Integer argc, char* argv[], Integer* arg,
                                     Real range[ 2 ] ) {

  PRE06( isValidArgs( argc, (const char**) argv ),
         argc > 0, arg, IN_RANGE( *arg, 1, argc - 1 ),
         ! strcmp( argv[ *arg ], "-elevation" ), range );
  CHECKING( Integer OLD( arg ) = *arg; )

  Integer result = 0;
  range[ MINIMUM ] = range[ MAXIMUM ] = 0.0;

  if ( *arg + 2 >= argc ) {
    failureMessage( "Missing parameter to command-line argument -elevation." );
  } else {
    char* end = 0;
    const Real minimum = strtod( argv[ *arg + 1 ], &end );

    if ( ! AND2( end != argv[ *arg + 1 ],
                 IN_RANGE( minimum,
                           MINIMUM_SURFACE_ELEVATION, MAXIMUM_ELEVATION ) ) ) {
      failureMessage( "Invalid command-line argument '%s'.", argv[ *arg + 1 ]);
    } else {
      const Real maximum = strtod( argv[ *arg + 2 ], &end );

      if ( ! AND2( end != argv[ *arg + 2 ],
                   IN_RANGE( maximum, minimum, MAXIMUM_ELEVATION ) ) ) {
        failureMessage( "Invalid command-line argument '%s'.", argv[*arg + 2]);
      } else {
        range[ MINIMUM ] = minimum;
        range[ MAXIMUM ] = maximum;
        *arg += 3;
        result = 1;
      }
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        AND3( IN_RANGE( range[ MINIMUM ],
                                        MINIMUM_SURFACE_ELEVATION,
                                        MAXIMUM_ELEVATION ),
                              IN_RANGE( range[ MAXIMUM ],
                                        range[ MINIMUM ], MAXIMUM_ELEVATION ),
                              *arg == OLD( arg ) + 3 ),
                        AND2( IS_ZERO2( range[ MINIMUM ], range[ MAXIMUM ] ),
                              *arg == OLD( arg ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readData - Read data from Ceilometer files and subset it by time,
         lon-lat box and elevation range.
INPUTS:  Data* data  data->arguments->listfile is the file containing the
         names of Ceilometer files to read and subset.
OUTPTUS: Data* data  data->profiles is the list of subset data to write.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void readData( Data* data ) {

  PRE04( data, data->ok, isValidArguments( &data->arguments ),
         data->profiles == 0 );

  Stream* listFile = newFileStream( data->arguments.listFile, "r" );
  data->ok = listFile != 0;

  if ( listFile ) {
    const Integer firstTimestamp = data->arguments.firstTimestamp;
    const Integer lastTimestamp  = data->arguments.lastTimestamp;
    DEBUG( fprintf( stderr, "readData for time range: %lld ... %lld\n",
                    firstTimestamp, lastTimestamp ); )

    /* For each file, read a subset of it into profile and append to list: */

    do {
      FileName fileName = "";
      listFile->readWord( listFile, fileName,
                          sizeof fileName / sizeof *fileName );
      data->ok = listFile->ok( listFile );

      if ( data->ok ) {
        DEBUG( fprintf( stderr, "reading Ceilometer file %s\n", fileName ); )
        Profile* const profile =
          readCeilometerFile( fileName, data->arguments.variable,
                              firstTimestamp, lastTimestamp,
                              (const Real (*)[2]) data->arguments.bounds,
                              data->arguments.elevationRange );
        data->ok = failureCount() == 0;
        
        if ( AND2( data->ok, profile ) ) {
          DEBUG( printProfile( profile ); )
          appendOrFreeProfile( data, profile );
        }
      }

      listFile->readString( listFile, fileName, 2 ); /* Read '\n'. */
    } while ( AND2( data->ok, ! listFile->isAtEnd( listFile ) ) );

    FREE_OBJECT( listFile );
  }

  if ( AND2( data->ok,
             OR2( data->profiles == 0,
                  data->profiles->count( data->profiles ) == 0 ) ) ) {
    failureMessage( "No profiles were in the subset." );
    data->ok = 0;
  }

  POST02( isValidArguments( &data->arguments ),
          IMPLIES( data->ok, isValidData( data ) ) );
}



/******************************************************************************
PURPOSE: appendOrFreeProfile - Append profile onto list or free it.
INPUTS:  Data* data  data->arguments->profiles List to create and/or append to.
         Profile* profiles[ 2 ]  Profiles to append or free.
OUTPTUS: Data* data              data->profiles  Created/appended list.
                                 data->ok
         Profile* profile        Appended or freed profile.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         profiles are deallocated and failureMessage() is called.
******************************************************************************/

static void appendOrFreeProfile( Data* data, Profile* profile ) {

  PRE03( data, IS_BOOL( data->ok ), profile );

  if ( data->ok ) {

    /* Create list if needed: */

    if ( data->profiles == 0 ) {
      data->profiles = newVoidList( deallocateProfile, 0 );
      data->ok = data->profiles != 0;
    }
  }

  if ( data->ok ) {
    data->profiles->insert( data->profiles, profile, LAST_ITEM );
    data->ok = data->profiles->ok( data->profiles );
  }

  if ( ! data->ok ) {
    deallocateProfile( profile );
    FREE_ZERO( profile );
  }

  POST0( IMPLIES( data->ok,
                  AND2( data->profiles,
                        data->profiles->has( data->profiles, profile ) ) ) );
}



/******************************************************************************
PURPOSE: readCeilometerFile - Read a subset of data from Ceilometer file.
INPUTS:  const char* fileName            Name of Ceilometer file to read.
         const char* variableName        Name of Ceilometer variable to read.
         const Integer firstTimestamp    Beginning timestamp of subset.
         const Integer lastTimestamp     Ending timestamp of subset.
         const Bounds bounds             Lon-lat bounds of subset.
         const Real elevationRange[ 2 ]  elevationRange[ MINIMUM, MAXIMUM ].
RETURNS: Profile* profile                Profile of subset data or
                                         0 if no data in subset.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.
******************************************************************************/

static Profile* readCeilometerFile( const char* fileName,
                                    const char* variableName,
                                    const Integer firstTimestamp,
                                    const Integer lastTimestamp,
                                    const Bounds bounds,
                                    const Real elevationRange[ 2 ] ) {

  PRE08( fileName, *fileName,
         variableName, *variableName,
         isValidYYYYMMDDHHMMSS( firstTimestamp ),
         isValidYYYYMMDDHHMMSS( lastTimestamp ),
         isValidBounds( bounds ),
         isValidElevationRange( elevationRange ) );

  Profile* result = 0;
  const Integer id = instrumentId( fileName );

  if ( id ) {
    const Integer dayBeforeFirstTimestamp = previousDay( firstTimestamp );
    const Integer profileTimestamp = parseFileTimestamp( fileName );

    if ( IN_RANGE( profileTimestamp, dayBeforeFirstTimestamp, lastTimestamp)) {
      double* data = 0;
      double* elevations = 0;
      double* timestamps = 0;
      double longitude = 0.0;
      double latitude  = 0.0;
      double elevation = 0.0;
      Real domain[ 3 ][ 2 ];
      char units[ 80 ] = "";
      size_t points = 0;
      domain[ LONGITUDE ][ MINIMUM ] = bounds[ LONGITUDE ][ MINIMUM ];
      domain[ LONGITUDE ][ MAXIMUM ] = bounds[ LONGITUDE ][ MAXIMUM ];
      domain[ LATITUDE  ][ MINIMUM ] = bounds[ LATITUDE  ][ MINIMUM ];
      domain[ LATITUDE  ][ MAXIMUM ] = bounds[ LATITUDE  ][ MAXIMUM ];
      domain[ ELEVATION ][ MINIMUM ] = bounds[ ELEVATION ][ MINIMUM ];
      domain[ ELEVATION ][ MAXIMUM ] = bounds[ ELEVATION ][ MAXIMUM ];

      /*
       * Read site location and a subset of valid data within domain and
       * time-range:
       */

      points =
        readSubsetCeilometerData( fileName, variableName,
                                  (const double (*)[2]) domain,
                                  firstTimestamp, lastTimestamp,
                                  units,
                                  &longitude, &latitude, &elevation,
                                  &data, &elevations, &timestamps );
      DEBUG( fprintf( stderr,
                      "Filter %lu points by date-time range = [%lld, %lld]\n",
                      points, firstTimestamp, lastTimestamp ); )

      if ( points ) {
        result = NEW( Profile, 1 );

        if ( result ) { /* Copy subset points to minimize cummulatve memory: */
          result->data = NEW( Real, points );
          result->elevations = result->data ? NEW( Real, points ) : 0;
          result->timestamps = result->elevations ? NEW( Real, points ) : 0;

          if ( ! result->timestamps ) {
            deallocateProfile( result );
            result = 0;
          } else {
            size_t point = 0;
            result->timestamp = profileTimestamp;
            result->points    = points;
            result->longitude = longitude;
            result->latitude  = latitude;
            result->elevation = elevation;
            result->id = id;
            memcpy( result->units, units, sizeof result->units );
            parseFileNotes( fileName, id, result->notes );

            /* Copy filtered data, elevations, timestamps: */

            for ( point = 0; point < points; ++point ) {
              const double timestamp = timestamps[ point ];
              const double value = data[ point ];
              const double elevation = elevations[ point ];
              result->data[ point ] = value;
              result->elevations[ point ] = elevation;
              result->timestamps[ point ] = timestamp;
            }
          }
        }

        FREE( data );
        FREE( elevations );
        FREE( timestamps );
      }
    }
  }

  DEBUG( fprintf( stderr, "readCeilometerFile result = %p\n", result ); )
  POST0( IMPLIES( result, isValidProfile( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: instrumentId - Parse instrument id from the name of a Ceilometer file.
INPUTS:  const char* fileName  Name of a Ceilometer file.
                               testdata/CL51_STMA_20211007.nc
RETURNS: Integer id if successful, else 0 and failureMessage() called.
NOTES:   Unlike NEUBrew files (which contain numeric ids) just hash part of
         file name before timestamp.
******************************************************************************/

static Integer instrumentId( const char* fileName ) {
  PRE03( fileName, *fileName, strlen( fileName ) > 12 );
  Integer result = 0;
  const Integer length = strlen( fileName );
  const Integer end = length - 12; /* Before timestamp. */
  Integer index = 0;

  for ( index = 0; index < end; ++index ) {
    result += fileName[ index ];
  }

  if ( result <= 0 ) {
    failureMessage( "Invalid Ceilometer file name '%s'.", fileName );
    result = 0;
  }

  DEBUG( fprintf( stderr, "instrumentId = %lld\n", result ); )
  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: parseFileTimestamp - Parse timestamp of a Ceilometer file data.
INPUTS:  const char* fileName   Name of Ceilometer file.
                                testdata/CL51_STMA_20211007.nc
RETURNS: Integer valid YYYYMMDDHHMMSS, else 0 and failureMessage() called.
******************************************************************************/

static Integer parseFileTimestamp( const char* fileName ) {

  PRE03( fileName, *fileName, strlen( fileName ) > 12 );

  Integer result = 0;
  const size_t offset = strlen( fileName ) - 11;
  const char* const yyyymmdd = fileName + offset;
  result = atoI( yyyymmdd );
  result *= 1000000; /* 000000 hhmmss. */
  DEBUG( fprintf( stderr, "file timestamp '%s' = %lld\n", yyyymmdd, result ); )

  if ( ! isValidYYYYMMDDHHMMSS( result ) ) {
      failureMessage( "Invalid/missing timestamp in Ceilometer file." );
      result = 0;
  }

  POST0( OR2( result == 0, isValidYYYYMMDDHHMMSS( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseFileNotes - Parse site/instrument name of a Ceilometer file.
INPUTS:  const char* fileName   Name of Ceilometer file.
                                testdata/CL51_STMA_20211007.nc
         const Integer id       File hashed id to use if parsing fails.
OUTPUTS: char notes[ 80 ]       Site name / instrument.
******************************************************************************/

static void parseFileNotes( const char* fileName, const Integer id,
                            char notes[ 80 ] ) {

  PRE04( fileName, *fileName, id > 0, notes );

  char* slash = strrchr( fileName, '/' );

  if ( slash ) {
    fileName = slash + 1; /* Skip past any path. */
  }

  memset( notes, 0, 80 * sizeof (char) );

  {
    size_t length = strlen( fileName );

    if ( length > 12 ) {
      int index = 0;
      length -= 12;

      if ( fileName[ length - 1 ] == '_' ) {
        --length;
      }

      if ( length > 79 ) {
        length = 79;
      }

      for ( index = 0; index < length; ++index ) {
        notes[ index ] = fileName[ index ];
      }
    }
  }

  if ( *notes == '\0' ) {
    snprintf( notes, 79, "%lld", id );
  }

  DEBUG( fprintf( stderr, "ParsedFileNotes(): notes = '%s'\n", notes ); )
  POST0( strlen( notes ) <= 79 );
}



/******************************************************************************
PURPOSE: writeData - Write subsetted profile data to stdout.
INPUTS:  Data* data  Data to write.
OUTPTUS: Data* data  data->ok.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void writeData( Data* data ) {
  PRE0( isValidData( data ) );
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    writeHeader( data, output );

    if ( data->ok ) {
      writeXDR( data, output );
    }

    FREE_OBJECT( output );
  }

  POST0( isValidData( data ) );
}



/******************************************************************************
PURPOSE: writeHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data      Data to write to stream.
         Stream* output  Stream to write data to.
OUTPUTS: Data* data      data->ok set to indicate success or failure.
         Stream* output  Stream to written to.
******************************************************************************/

static void writeHeader( Data* data, Stream* output ) {

  PRE08( isValidData( data ),
         data->ok,
         output,
         output->invariant( output ),
         output->isWritable( output ),
         data->ok == output->ok( output ),
         data->profiles,
         data->profiles->count( data->profiles ) > 0 );

  const Arguments* const arguments = &data->arguments;
  const VoidList*  const profiles = data->profiles;
  const Profile* const profile = profiles->item( profiles, 0 );
  UTCTimestamp firstTimestamp;
  UTCTimestamp lastTimestamp;
  toUTCTimestamp2( arguments->firstTimestamp, firstTimestamp );
  toUTCTimestamp2( arguments->lastTimestamp,  lastTimestamp );

  output->writeString( output,
                       "Profile 2.0\n"
                       "%s\n"
                       "%s %s\n"
                       "# Subset domain:"
                       " <min_lon> <min_lat> <max_lon> <max_lat>:\n"
                       "%lg %lg %lg %lg\n"
                       "# Dimensions: variables profiles:\n"
                       "%lld %lld\n"
                       "# Variable names:\n"
                       "timestamp id longitude latitude elevation %s\n"
                       "# Variable units:\n"
                       "yyyymmddhhmmss - deg deg m %s\n"
                       "# char notes[profiles][80] and\n"
                       "# MSB 64-bit integers points[profiles] and\n"
                       "# IEEE-754 64-bit reals"
                       " data_1[variables][points_1] ..."
                       " data_P[variables][points_P]:\n",
                       arguments->description,
                       firstTimestamp, lastTimestamp,
                       arguments->bounds[ LONGITUDE ][ MINIMUM ],
                       arguments->bounds[ LATITUDE  ][ MINIMUM ],
                       arguments->bounds[ LONGITUDE ][ MAXIMUM ],
                       arguments->bounds[ LATITUDE  ][ MAXIMUM ],
                       VARIABLES, profiles->count( profiles ),
                       arguments->variable, profile->units );

  data->ok = output->ok( output );

  POST04( isValidData( data ), output->invariant( output ),
          output->isWritable( output ), data->ok == output->ok( output ) );
}



/******************************************************************************
PURPOSE: writeXDR - Write XDR format data arrays of subset to stdout.
INPUTS:  Data* data      Data to write to stream.
         Stream* output  Stream to write data to.
OUTPUTS: Data* data      data->ok set to indicate success or failure.
         Stream* output  Stream to written to.
******************************************************************************/

static void writeXDR( Data* data, Stream* output ) {

  PRE06( isValidData( data ), data->ok, output, output->invariant( output ),
         output->isWritable( output ), output->ok( output ) );

  const VoidList* profiles = data->profiles;
  writeProfileNotes( profiles, output );
  writeProfilePoints( profiles, output );
  data->ok = output->ok( output );

  if ( data->ok ) {
    data->ok = writeProfileData( profiles, output );
  }

  POST04( isValidData( data ), output->invariant( output ),
          output->isWritable( output ),
          IMPLIES( data->ok, output->ok( output ) ) );
}



/******************************************************************************
PURPOSE: writeProfileNotes - Write profile notes.
INPUTS:  const VoidList* profiles   List of Profile*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeProfileNotes( const VoidList* profiles, Stream* output ) {

  PRE07( profiles, profiles->invariant( profiles ),
         profiles->count( profiles ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer profileCount = profiles->count( profiles );
  Integer profileIndex = 0;

  do {
    const Profile* const profile = profiles->item( profiles, profileIndex );
    CHECK( isValidProfile( profile ) );
    output->writeString( output, "%-79s\n", profile->notes );
    ++profileIndex;
  } while ( AND2( output->ok( output ), profileIndex < profileCount ) );

  POST04( profiles->invariant( profiles ), profiles->count( profiles ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: writeProfilePoints - Write MSB 64-bit integer profile subset point
         counts.
INPUTS:  const VoidList* profiles   List of Profile*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeProfilePoints( const VoidList* profiles, Stream* output ) {

  PRE07( profiles, profiles->invariant( profiles ),
         profiles->count( profiles ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer profileCount = profiles->count( profiles );
  Integer profileIndex = 0;

  do {
    const Profile* const profile = profiles->item( profiles, profileIndex );
    CHECK( isValidProfile( profile ) );
    output->write64BitInteger( output, profile->points );
    ++profileIndex;
  } while ( AND2( output->ok( output ), profileIndex < profileCount ) );

  POST04( profiles->invariant( profiles ), profiles->count( profiles ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: writeProfileData - Write 64-bit IEEE-754 profile subset variable data.
INPUTS:  const VoidList* profiles   List of Profile*.
         Stream*         output  Stream to write to.
RETURNS: Integer 1 if successful, else 0 and print failure message to stderr.
******************************************************************************/

static Integer writeProfileData( const VoidList* profiles, Stream* output ) {

  PRE07( profiles, profiles->invariant( profiles ),
         profiles->count( profiles ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer profileCount = profiles->count( profiles );
  const Integer tempSize = maximumProfilePoints( profiles );
  Real* temp = NEW( Real, tempSize ); /* To hold replicated ids, etc. */
  Integer result = 0;

  if ( temp ) {
    Integer profileIndex = 0;

    do {
      const Profile* const profile = profiles->item( profiles, profileIndex );
      CHECK( isValidProfile( profile ) );
      output->write64BitReals( output, profile->timestamps, profile->points );

      if ( output->ok( output ) ) {
        copyValue( profile->points, temp, profile->id );
        output->write64BitReals( output, temp, profile->points );

        if ( output->ok( output ) ) {
          copyValue( profile->points, temp, profile->longitude );
          output->write64BitReals( output, temp, profile->points );

          if ( output->ok( output ) ) {
            copyValue( profile->points, temp, profile->latitude );
            output->write64BitReals( output, temp, profile->points );

            if ( output->ok( output ) ) {
              output->write64BitReals( output, profile->elevations,
                                       profile->points );

              if ( output->ok( output ) ) {
                output->write64BitReals( output, profile->data,
                                         profile->points );
              }
            }
          }
        }
      }

      ++profileIndex;
    } while ( AND2( output->ok( output ), profileIndex < profileCount ) );

    FREE( temp );
    result = output->ok( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: copyValue - Copy value to array.
INPUTS:  const Integer count  Size of array.
         const Real value     Value to copy.
OUTPUTS: Real array[ count ]  Array filled with value.
******************************************************************************/

static void copyValue( const Integer count, Real array[], const Real value ) {
  PRE03( count > 0, array, ! isNan( value ) );
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {
    array[ index ] = value;
  }

  POST02( array[ 0 ] == value, array[ count - 1 ] == value );
}



/******************************************************************************
PURPOSE: maximumProfilePoints - Maximum number of profile points.
INPUTS:  const VoidList* profiles  List of profiles to check.
RETURNS: Integer number of points in the largest profile.
******************************************************************************/

static Integer maximumProfilePoints( const VoidList* profiles ) {
  PRE03( profiles,
         profiles->invariant( profiles ),
         profiles->count( profiles ) > 0 );
  const Integer count = profiles->count( profiles );
  Integer result = 0;
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {
    const Profile* const profile = (const Profile*)
      profiles->item( profiles, index );
    const Integer points = profile->points;

    if ( points > result ) {
      result = points;
    }
  }

  POST0( result > 0 );
  return result;
}




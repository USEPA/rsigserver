
/******************************************************************************
PURPOSE: MOZAICSubset.c - Read a set of a MOZAIC files, subset the data to a
         bounds (longitude-latitude rectangle) and write it to
         stdout as XDR (IEEE-754) format binary.

NOTES:   Uses libUtilities.a (../../libs/Utilities).

HISTORY: 2010-01-26 plessel.todd@epa.gov, Created.
STATUS: unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>     /* For printf(). */
#include <math.h>      /* For hypot(), atan2(). */
#include <string.h>    /* For strlen(). */
#include <ctype.h>     /* For isdigit(), isalnum(), isspace(), isprint(). */

#include <Utilities.h> /* For PRE0*(), NEW_ZERO(), Stream, VoidList, Note. */

/*================================== TYPES ==================================*/

/* MOZAIC variables: */

enum {
  AIRCRAFT_TIMESTAMP, AIRCRAFT_LONGITUDE, AIRCRAFT_LATITUDE,AIRCRAFT_ELEVATION,
  RADIO_ALTITUDE, PRESSURE, TEMPERATURE, AIR_SPEED, GROUND_SPEED,
  WIND_U, WIND_V, OZONE, H2O_STATIC_TEMPERATURE,
  RELATIVE_HUMIDITY, H2O, CO, NOY, NO, NOX, VARIABLES,
  IMPLICIT_VARIABLES = 4
};

static const char* const variableNames[ VARIABLES ] = {
  "timestamp", "longitude", "latitude", "elevation",
  "radio_altitude", "pressure", "temperature", "air_speed", "ground_speed",
  "wind_u", "wind_v", "ozone", "h2o_static_temperature",
  "relative_humidity", "h2o", "co", "noy", "no", "nox"
};

static const char* const variableUnits[ VARIABLES ] = {
  "yyyymmddhhmmss", "deg", "deg", "m",
  "m", "Pa", "C", "m/s", "m/s",
  "m/s", "m/s", "ppmV", "C",
  "%", "g/kg", "ppmV", "ppmV", "ppmV", "ppmV"
};

/* Scale factors (data multipliers) to convert units to those above: */

#define ppb2ppm 0.001

static const Real unitScales[ VARIABLES ] = {
  1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, 1.0, 1.0, 1.0,
  1.0, 1.0, ppb2ppm, 1.0,
  1.0, 1.0, ppb2ppm, ppb2ppm, ppb2ppm, ppb2ppm
};

#undef ppb2ppm

static void scaleData( Real data[ VARIABLES ] ) {
  Integer index = 0;

  for ( index = 0; index < VARIABLES; ++index ) {
    data[ index ] *= unitScales[ index ];
  }
}

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;              /* File listing MOZAIC files to read. */
  const char* description;           /* User-supplied description. */
  Bounds      bounds;          /*bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]*/
  Integer     firstTimestamp;        /* YYYYMMDDHHMMSS of subset. */
  Integer     lastTimestamp;         /* YYYYMMDDHHMMSS of subset. */
  Integer     selected[ VARIABLES ]; /* User-specified output variables. */
} Arguments;

#ifndef NO_ASSERTIONS

/* Arguments invariant: */

static Integer isValidArguments( const Arguments* arguments ) {
  const Integer result =
    AND12( arguments,
           arguments->listFile,
           arguments->listFile[ 0 ],
           arguments->description,
           arguments->description[ 0 ],
           isValidBounds( arguments->bounds ),
           isValidYYYYMMDDHHMMSS( arguments->firstTimestamp ),
           isValidYYYYMMDDHHMMSS( arguments->lastTimestamp ),
           arguments->firstTimestamp <= arguments->lastTimestamp,
           minimumItemI( arguments->selected, VARIABLES ) >= 0,
           maximumItemI( arguments->selected, VARIABLES ) == 1,
           IN_RANGE( sumI( arguments->selected, VARIABLES ),
                     IMPLICIT_VARIABLES + 1, VARIABLES ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Track: Result of reading a subset of a MOZAIC aircraft data file. */


typedef struct {
  Note note;              /* E.g., "MD20060703014:FRANKFURT->ATLANTA". */
  Bounds  bounds;         /* bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  Integer firstTimestamp; /* YYYYDDDHHMMSS. */
  Integer lastTimestamp;  /* YYYYDDDHHMMSS. */
  Integer variables;      /* Selected variables <= VARIABLES. */
  Integer points;         /* Number of data points. */
  Real*   data;           /* data[ variables ][ points ]. */
} Track;


/* Track destructor: */

static void deallocateTrack( void* vtrack ) {
  PRE0( vtrack );
  Track* track = (Track*) vtrack;
  FREE( track->data );
  ZERO_OBJECT( track );
  POST0( vtrack );
}

#ifndef NO_ASSERTIONS

/* Track invariant: */

static Integer isValidTrack( const Track* track ) {
  const Integer result =
    AND12( track,
           isalnum( track->note[ 0 ] ),
           track->note[ NOTE_LENGTH ] == '\0',
           strchr( track->note, ':' ),
           strstr( track->note, "->" ),
           strlen( track->note ) == NOTE_LENGTH,
           isValidYYYYMMDDHHMMSS( track->firstTimestamp ),
           isValidYYYYMMDDHHMMSS( track->lastTimestamp ),
           track->firstTimestamp <= track->lastTimestamp,
           IN_RANGE( track->variables, IMPLICIT_VARIABLES + 1, VARIABLES ),
           IN_RANGE( track->points, 1, INTEGER_MAX / VARIABLES ),
           isNanFree( track->data + track->points, /* Skip timestamp. */
                      ( track->variables - 1 ) * track->points ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Data type: */

typedef struct {
  Arguments arguments; /* User-supplied (command-line) arguments. */
  VoidList* tracks;    /* List of subsetted tracks. */
  Integer   ok;        /* Did last command succeed? */
} Data;

/* Data destructor: */

static void deallocateData( Data* data ) {
  PRE0( data );
  FREE_OBJECT( data->tracks ); /* Calls deallocateTrack(). */
  ZERO_OBJECT( data );
  POST0( data );
}

#ifndef NO_ASSERTIONS

/* Data invariant: */

static Integer isValidData( const Data* data ) {
  Integer result =
    AND5( data,
          isValidArguments( &data->arguments ),
          data->tracks,
          data->tracks->invariant( data->tracks ),
          IS_BOOL( data->ok ) );

  if ( result ) {
    const VoidList* const tracks = data->tracks;
    const Integer trackCount = tracks->count( tracks );
    Integer index = 0;

    do {
      const Track* const track = tracks->item( tracks, index );
      result = isValidTrack( track );
      ++index;
    } while ( AND2( result, index < trackCount ) );
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

static Integer parseVariables( Integer argc, char* argv[], Integer* arg,
                               Integer selected[ VARIABLES ] );

static void readData( Data* data );

static Track* readMOZAICFile( const char* fileName,
                              const Integer firstTimestamp,
                              const Integer lastTimestamp,
                              const Integer selected[ VARIABLES ],
                              const Bounds bounds );

static Integer isProfileFile( const char* fileName, const char* fileData,
                              Integer* firstTimestamp, Integer* lastTimestamp,
                              Integer* seconds1, Integer* seconds2 );

static Integer isValidDataLine( const Integer selected[ VARIABLES ],
                                const Real variables[ VARIABLES ],
                                Integer relativeHumidityValidity,
                                Integer noValidity );

static Integer parseDataLine( const char* headerLine,
                              const char* dataLine,
                              Integer points,
                              Integer firstTimestamp,
                              Integer lastTimestamp,
                              Integer profile,
                              Integer profileFirstTimestamp,
                              Integer profileLastTimestamp,
                              Integer seconds1,
                              Integer seconds2,
                              const Bounds bounds,
                              const Integer selected[ VARIABLES ],
                              Real variables[ VARIABLES ] );

static Integer profileTimestamp( Integer point, Integer points,
                                 Integer profile,
                                 Integer firstTimestamp, Integer lastTimestamp,
                                 Integer seconds1, Integer seconds2 );

static Integer timestampOfFileName( const char* fileName );

static void parseNote( const char* fileData, Note note );

static Track* copySubsetData( const Real* data,
                              const Integer subsetVariables,
                              const Integer subsetPoints,
                              const Integer reverse,
                              const Note note );

static Integer totalSubsetPoints( const VoidList* tracks );

static void writeData( Data* data );

static void writeHeader( Data* data, Stream* output );

static void writeXDR( Data* data, Stream* output );

static void writeTrackNotes( const VoidList* tracks, Stream* output );

static void writeTrackBounds( const VoidList* tracks, Stream* output );

static void writeTrackPoints( const VoidList* tracks, Stream* output );

static void writeTrackData( const VoidList* tracks, Stream* output );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Read a subset of a list MOZAIC files and write it to stdout in
         XDR format.
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
    checkForTest( &argc, argv ); /* Check for and remove any -test arguments.*/
    ZERO_OBJECT( &data );
    data.ok = parseArguments( argc, argv, &data.arguments );

    if ( data.ok ) {
      readData( &data ); /* From list file named in arguments. */

      if ( data.ok ) {
        writeData( &data ); /* To stdout. */
      }
    }

    ok = data.ok;
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
  Integer variable = 0;
  fprintf( stderr,
           "\n\n%s - Read a set of MOZAIC files and extract track\n",
           programName );
  fprintf( stderr, "data for selected variables" );
  fprintf( stderr, " subsetted by a lon-lat rectangle.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "  -files <listFile> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -time <yyyymmddhhmmss> <yyyymmddhhmmss> \\\n" );
  fprintf( stderr, "  [ -variable" );

  for ( variable = IMPLICIT_VARIABLES; variable < VARIABLES; ++variable ) {
    fprintf( stderr, " %s", variableNames[ variable ] );
  }

  fprintf( stderr, " ] \\\n" );
  fprintf( stderr, "  [ -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ] \\\n" );
  fprintf( stderr, "\n\n" );
  fprintf( stderr, "Note: times are in UTC (GMT)\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files /mozaic/data/files.txt \\\n" );
  fprintf(stderr,"-desc http://mozaic.aero.obs-mip.fr/web/,MOZAICSubset \\\n");
  fprintf( stderr, "-time 20060703000000 20060703235959 \\\n" );
  fprintf( stderr, "-variable ozone \\\n" );
  fprintf( stderr, "-domain -84 33 -82 34 > subset.xdr\n\n" );
  fprintf( stderr, "Subset of data for July 3, 2006 near Atlanta, GA, USA\n");
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays\n" );
  fprintf( stderr, "For example:\n" );
  fprintf( stderr, "AIRCRAFT 2.0\n" );
  fprintf( stderr, "http://mozaic.aero.obs-mip.fr/web/,MOZAICSubset\n" );
  fprintf( stderr, "2006-07-03T00:00:00-0000 2006-07-03T23:59:59-0000\n" );
  fprintf( stderr, "# Subset domain: <min_lon> <min_lat> <max_lon> <max_lat>");
  fprintf( stderr, ":\n");
  fprintf( stderr, "-85 33 -82 34\n" );
  fprintf( stderr, "# Dimensions: variables points tracks:\n" );
  fprintf( stderr, "5 48 2\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "timestamp longitude latitude elevation ozone\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss deg deg m ppb\n" );
  fprintf( stderr, "# char notes[tracks][80] and" );
  fprintf( stderr, "# IEEE-754 64-bit reals bounds[tracks][2=lon,lat][2=min,max]" );
  fprintf( stderr, " and\n" );
  fprintf( stderr, "# MSB 64-bit integers points[tracks] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals" );
  fprintf( stderr, " data_1[points_1][variables] ..." );
  fprintf( stderr, " data_T[points_T][variables]:\n" );
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

  if ( ! IN_RANGE( argc, 8, 14 + VARIABLES - IMPLICIT_VARIABLES ) ) {
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
          result = parseOptionalArguments( argc, argv, &arg, arguments );
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
  Integer variable = 0;
  ZERO_OBJECT( arguments );
  arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
  arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
  arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
  arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;

  for ( variable = 0; variable < VARIABLES; ++variable ) {
    arguments->selected[ variable ] = 1;
  }
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
  Integer parsedVariable = 0;
  Integer parsedBounds   = 0;

  while ( *arg < argc ) {

    if ( ( AND2( ! strcmp( argv[ *arg ], "-variable" ), ! parsedVariable ))) {
      parsedVariable = 1;
      result = parseVariables( argc, argv, arg, arguments->selected );
    } else if ( AND2( ! strcmp( argv[ *arg ], "-domain" ), ! parsedBounds )) {
      parsedBounds = 1;
      result = parseBounds( argc, argv, arg, arguments->bounds );
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
PURPOSE: parseVariables - Parse command-line arguments for -variable.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
         Integer* arg          Index of argument to parse.
OUTPUTS: Integer* arg          Index of next argument to parse.
         Integer selected[ VARIABLES ]  Set to 1 if named, else 0.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer parseVariables( Integer argc, char* argv[], Integer* arg,
                               Integer selected[ VARIABLES ] ) {

  PRE06( isValidArgs( argc, (const char**) argv ),
         argc > 0, arg, IN_RANGE( *arg, 1, argc - 1 ),
         ! strcmp( argv[ *arg ], "-variable" ), selected );
  CHECKING( Integer OLD( arg ) = *arg; )

  Integer result = 0;
  memset( selected + IMPLICIT_VARIABLES, 0,
          ( VARIABLES - IMPLICIT_VARIABLES ) * sizeof (Integer) );

  if ( *arg + 1 >= argc ) {
    failureMessage( "Missing parameter to command-line argument -variables." );
  } else {
    Integer ok = 1;
    ++*arg;

    while ( AND4( ok, *arg < argc, argv[ *arg ][ 0 ],
                  argv[ *arg ][ 0 ] != '-' ) ) {
      const char* const variableName = argv[ *arg ];
      const Integer variable =
        indexOfString( variableName, variableNames, VARIABLES );

      if ( OR3( variable == -1, variable < IMPLICIT_VARIABLES,
                selected[ variable ] ) ) {
        failureMessage( "Invalid/redundant variable name %s.", variableName );
        ok = 0;
      } else {
        selected[ variable ] = 1;
        ++*arg;
        result = 1;
      }
    }
  }

  if ( ! result ) {
    memset( selected, 0, VARIABLES * sizeof (Integer) );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, AND4( minimumItemI( selected, VARIABLES ) >= 0,
                                 maximumItemI( selected, VARIABLES ) == 1,
                                 sumI( selected, VARIABLES ) >
                                   IMPLICIT_VARIABLES,
                                 *arg == OLD( arg ) +
                                         1 + sumI( selected, VARIABLES ) -
                                         IMPLICIT_VARIABLES ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readData - Read data from MOZAIC files and subset it by time,
         lon-lat box and selected variables.
INPUTS:  Data* data  data->arguments->listfile is the file containing the
         names of MOZAIC files to read and subset.
OUTPTUS: Data* data  data->tracks is the list of subset data to write.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void readData( Data* data ) {

  PRE04( data, data->ok, isValidArguments( &data->arguments ),
         data->tracks == 0 );

  Stream* listFile = newFileStream( data->arguments.listFile, "r" );
  data->ok = listFile != 0;

  if ( listFile ) {
    const Integer firstTimestamp = data->arguments.firstTimestamp;
    const Integer lastTimestamp  = data->arguments.lastTimestamp;
    const Integer fileTimestamp  = previousDay( firstTimestamp );

    /* For each file, read a subset of it into track and append to list: */

    do {
      FileName fileName = "";
      listFile->readWord( listFile, fileName,
                          sizeof fileName / sizeof *fileName );
      data->ok = listFile->ok( listFile );

      if ( data->ok ) {
        const Integer currentTimestamp = timestampOfFileName( fileName );
        DEBUG( fprintf( stderr, "listing MOZAIC file %s\n", fileName ); )
        data->ok = currentTimestamp > 0;

        if ( ! data->ok ) {
          failureMessage( "Invalid MOZAIC file %s.", fileName );
        } else if (IN_RANGE(currentTimestamp, fileTimestamp, lastTimestamp)) {
          Track* track =
            readMOZAICFile( fileName, firstTimestamp, lastTimestamp,
                            data->arguments.selected,
                            (const Real (*)[2]) data->arguments.bounds );

          if ( track ) {

            if ( data->tracks == 0 ) { /* Create list if needed: */
              data->tracks = newVoidList( deallocateTrack, 0 );
              data->ok = data->tracks != 0;
            }

            if ( data->tracks ) { /* Append subsetted track to list: */
              data->tracks->insert( data->tracks, track, LAST_ITEM );
              data->ok = data->tracks->ok( data->tracks );
            }

            if ( ! data->ok ) {
              deallocateTrack( track );
              FREE( track );
            }
          }
        }
      }

      listFile->readString( listFile, fileName, 2 ); /* Read '\n'. */
    } while ( AND2( data->ok, ! listFile->isAtEnd( listFile ) ) );

    FREE_OBJECT( listFile );
  }

  if ( AND2( data->ok, data->tracks == 0 ) ) {
    failureMessage( "No tracks were in the subset." );
    data->ok = 0;
  }

  POST02( isValidArguments( &data->arguments ),
          IMPLIES( data->ok, isValidData( data ) ) );
}



/******************************************************************************
PURPOSE: readMOZAICFile - Read a subset of track data from MOZAIC file.
INPUTS:  const char* fileName          Name of compressed MOZAIC file to read.
         const Integer firstTimestamp  Beginning timestamp of subset.
         const Integer lastTimestamp   Ending timestamp of subset.
         const Integer selected[ VARIABLES ]  Selected variables of subset.
         const Bounds bounds                  Lon-lat bounds of subset.
RETURNS: Track* track  Track of subset data or 0 if no data in subset.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.
******************************************************************************/

static Track* readMOZAICFile( const char* fileName,
                              const Integer firstTimestamp,
                              const Integer lastTimestamp,
                              const Integer selected[ VARIABLES ],
                              const Bounds bounds ) {

  PRE08( fileName,
         isValidYYYYMMDDHHMMSS( firstTimestamp ),
         isValidYYYYMMDDHHMMSS( lastTimestamp ),
         selected,
         minimumItemI( selected, VARIABLES ) >= 0,
         maximumItemI( selected, VARIABLES ) == 1,
         sumI( selected, VARIABLES ) > IMPLICIT_VARIABLES,
         isValidBounds( bounds ) );

  Track* result = 0;
  Integer length = 0;
  char* fileData = readFile( fileName, &length );

  if ( fileData ) {
    Integer profileFirstTimestamp = 0;
    Integer profileLastTimestamp = 0;
    Integer seconds1 = 0;
    Integer seconds2 = 0;
    const Integer profile =
      isProfileFile( fileName, fileData,
                     &profileFirstTimestamp, &profileLastTimestamp,
                     &seconds1, &seconds2 );

    if ( IMPLIES( profile,
                  AND2( IN_RANGE( profileFirstTimestamp,
                                  firstTimestamp, lastTimestamp ),
                        IN_RANGE( profileLastTimestamp,
                                  firstTimestamp, lastTimestamp ) ) ) ) {
      const Integer headerLines = profile ? 5 : 3;
      const Integer points = linesInString( fileData ) - headerLines;
      DEBUG( fprintf( stderr, "Reading MOZAIC file %s, points = %lld\n",
                      fileName, points ); )

      if ( points > 0 ) {
        const Integer subsetVariables = sumI( selected, VARIABLES );
        Integer subsetPoints = 0;
        Real* data = NEW_ZERO( Real, subsetVariables * points );

        if ( data ) {
          Real* output = data;
          const char* const headerLine = skipLines( fileData, headerLines - 1);
          const char* dataLine = skipLines( fileData, headerLines );

          while ( dataLine ) {
            Real variables[ VARIABLES ];

            if ( parseDataLine( headerLine, dataLine, points,
                                firstTimestamp, lastTimestamp,
                                profile, profileFirstTimestamp,
                                profileLastTimestamp, seconds1, seconds2,
                                bounds, selected, variables ) ) {
              Integer variable = 0;

              for ( variable = 0; variable < VARIABLES; ++variable ) {

                if ( selected[ variable ] ) {
                  *output++ = variables[ variable ];
                }
              }

              ++subsetPoints;
            }

            dataLine = skipLines( dataLine, 1 );
          }

          if ( subsetPoints ) {
            Note note = "";
            memset( note, 0, sizeof note );
            parseNote( fileData, note );
            result =
              copySubsetData( data, subsetVariables, subsetPoints,
                              profile == -1, note );
          }

          FREE( data );
        }
      }
    }

    FREE( fileData );
  }

  POST0( IMPLIES( result, isValidTrack( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isProfileFile - Is it a profile (ascending or descending) MOZIAC file?
INPUTS:  const char* fileName     Name of MOZAIC file.
         const char* fileData     String containing file contents.
OUTPUTS: Integer* firstTimestamp  Beginning timestamp of profile file.
         Integer* lastTimestamp   Ending timestamp of profile file.
         Integer* seconds1        Total seconds from Jan 1 of year of
                                  firstTimestamp to firstTimestamp.
         Integer* seconds2        Total seconds from Jan 1 of year of
                                  firstTimestamp to lastTimestamp.
RETURNS: Integer 0 = non-profile, 1 = ascent, -1 = descent.
******************************************************************************/

static Integer isProfileFile( const char* fileName, const char* fileData,
                              Integer* firstTimestamp, Integer* lastTimestamp,
                              Integer* seconds1, Integer* seconds2 ) {

  Integer result = 0;

  /*
   Check if file name contains a path e.g.,
   /data/MOZAIC/data/2006/MP2006/MD20060703014.txt
   then parse timestamp from header that looks like:
   LA2006-04-17,KFA2007-07-15,CNRM2005-06-09,Version5.4
   flight: MA20060417015 NL=070 SP=010000 SHANGHAI->VIENNA
   START at: Lat 31.13,Lon  121.82,date 20060417,time 015141
   END   at: Lat 32.63,Lon  119.45,date 20060417,time 022319
  */

  const char* name = strrchr( fileName, '/' );
  const char* const start = strstr( fileData, "START " );
  const char* const end   = strstr( fileData, "END " );

  /* Skip directory, if present: */

  if ( name ) {
    ++name;
  } else {
    name = fileName;
  }

  *firstTimestamp = *lastTimestamp = *seconds1 = *seconds2 = 0;

  if ( AND6( name, *name == 'M', strlen( name ) > 9,
             start, end, start < end ) ) {
    result = name[ 1 ] == 'A' ? 1 : name[ 1 ] == 'D' ? -1 : 0;

    if ( result ) {
      const char* const startDate = strstr( start, "ate " ); /* Date or date */
      const char* const startTime = strstr( start, "ime " ); /* Time or time */
      const char* const endDate = strstr( end, "ate " );
      const char* const endTime = strstr( end, "ime " );

      if ( AND6( startDate, startTime, endDate, endTime,
                 startDate < startTime, endDate < endTime ) ) {
        const Integer yyyymmdd1 = atoI( startDate + 4 );
        const Integer hhmmss1   = atoI( startTime + 4 );
        const Integer yyyymmdd2 = atoI( endDate + 4 );
        const Integer hhmmss2   = atoI( endTime + 4 );
        const Integer yyyymmddhhmmss1 = yyyymmdd1 * 1000000 + hhmmss1;
        const Integer yyyymmddhhmmss2 = yyyymmdd2 * 1000000 + hhmmss2;

        if ( AND2( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss1 ),
                   isValidYYYYMMDDHHMMSS( yyyymmddhhmmss2 ) ) ) {

          if ( AND2( result == 1, yyyymmddhhmmss1 <= yyyymmddhhmmss2 ) ) {
            *firstTimestamp = yyyymmddhhmmss1;
            *lastTimestamp  = yyyymmddhhmmss2;
          } else if ( AND2( result == -1,
                            yyyymmddhhmmss2 <= yyyymmddhhmmss1 ) ) {
            *firstTimestamp = yyyymmddhhmmss2;
            *lastTimestamp  = yyyymmddhhmmss1;
          }

          if ( *firstTimestamp ) {
            totalSeconds( *firstTimestamp, *lastTimestamp, seconds1, seconds2);
            DEBUG( fprintf( stderr, "%lld ... %lld seconds: %lld ... %lld\n",
                   *firstTimestamp, *lastTimestamp, *seconds1, *seconds2 ); )
          }
        }
      }
    }
  }

  if ( *seconds1 == 0 ) {
    result = *firstTimestamp = *lastTimestamp = *seconds2 = 0;
  }

  POST02( IN_RANGE( result, -1, 1 ),
          IMPLIES_ELSE( result,
                        AND5( isValidYYYYMMDDHHMMSS( *firstTimestamp ),
                              isValidYYYYMMDDHHMMSS( *lastTimestamp ),
                              *firstTimestamp <= *lastTimestamp,
                              *seconds1 >= 0,
                              *seconds2 >= *seconds1 ),
                        IS_ZERO4( *firstTimestamp, *lastTimestamp,
                                  *seconds1, *seconds2 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseDataLine - Parse, validate and subset a line of MOZAIC file data.
INPUTS:  const char* headerLine   Line of header from MOZAIC file.
         const char* dataLine     Line of data from MOZAIC file.
         Integer points           Number of data points in file.
         Integer firstTimestamp   Beginning timestamp of subset.
         Integer lastTimestamp    Ending timestamp of subset.
         Integer profile          0 = non-profile, 1 = ascent, -1 = descent.
         Integer profileFirstTimestamp  Beginning timestamp of profile.
         Integer profileLastTimestamp   Ending timestamp of profile.
         Integer seconds1               Total seconds from Jan 1 yyyy of
                                        firstTimestamp to firstTimestamp.
         Integer seconds2               Total seconds from Jan 1 yyyy of
                                        firstTimestamp to lastTimestamp.
         const Bounds bounds      Lon-lat bounds of subset.
         const Integer selected[ VARIABLES ]  Selected variables of subset.
OUTPUTS: Real variables[ VARIABLES ]  Data variable values if in subset.
RETURNS: Integer 1 if valid and in subset, else 0.
******************************************************************************/

static Integer parseDataLine( const char* headerLine,
                              const char* dataLine,
                              Integer points,
                              Integer firstTimestamp,
                              Integer lastTimestamp,
                              Integer profile,
                              Integer profileFirstTimestamp,
                              Integer profileLastTimestamp,
                              Integer seconds1,
                              Integer seconds2,
                              const Bounds bounds,
                              const Integer selected[ VARIABLES ],
                              Real variables[ VARIABLES ] ) {

  PRE014( headerLine,
          dataLine,
          points > 0,
          isValidYYYYMMDDHHMMSS( firstTimestamp ),
          isValidYYYYMMDDHHMMSS( lastTimestamp ),
          firstTimestamp <= lastTimestamp,
          IN4( profile, 0, 1, -1 ),
          IMPLIES_ELSE( profile,
                        AND5( isValidYYYYMMDDHHMMSS( profileFirstTimestamp ),
                              isValidYYYYMMDDHHMMSS( profileLastTimestamp ),
                              profileFirstTimestamp <= profileLastTimestamp,
                              seconds1 >= 0,
                              seconds2 >= seconds1 ),
                        IS_ZERO4( profileFirstTimestamp,
                                  profileLastTimestamp,
                                  seconds1,
                                  seconds2 ) ),
          isValidBounds( bounds ),
          selected,
          minimumItemI( selected, VARIABLES ) >= 0,
          maximumItemI( selected, VARIABLES ) == 1,
          sumI( selected, VARIABLES ) > IMPLICIT_VARIABLES,
          variables );

  Integer result = 0;
  Integer yyyymmdd = 0;
  Integer hhmmss   = 0;
  Integer relativeHumidityValidity = 0;
  Integer noValidity = 0;

  /*
   new columns:
    0 Level or Date
    1 Level_Altitude or Time
    2 Latitude
    3 Longitude
    4 Baro_Altitude
    5 Radio_Altitude
    6 GPS_Altitude
    7 Pressure
    8 Aircraft_Air_Speed
    9 Aircraft_Ground_Speed
    10 Aircraft_Static_Temperature
    11 Air_Temperature
    12 Zonal_Wind
    13 Meridian_Wind
    14 O3
    15 H2O_Static_Temperature
    16 Relative_Humidity
    17 Relative_Humidity_Validity
    18 Relative_Humidity_Accuracy
    19 H2O
    20 CO
    21 NOy
    22 NO
    23 NOx
    24 NOy_Uncertainty
    25 NOy_Validity
   */

  const char* const startsWith =
    profile ? "Level Level_Altitude "
    : "Date Time Latitude Longitude Baro_Altitude";
  const Integer isNewFormat =
    ! strncmp( headerLine, startsWith, strlen( startsWith ) );
  const void* const unused_ =
    memset( variables, 0, VARIABLES * sizeof *variables );
  Real windDirection = 0.0;
  Real windSpeed = 0.0;
  const Integer ok =
    isNewFormat ?
      sscanf( dataLine,
              "%lld %lld %lf %lf %lf %lf %*s %lf %lf %lf %*s %lf %lf %lf "
              "%lf %lf %lf %lld %*s %lf %lf %lf %lf %lf %*s %lld\n",
              &yyyymmdd, &hhmmss,
              &variables[ AIRCRAFT_LATITUDE ],
              &variables[ AIRCRAFT_LONGITUDE ],
              &variables[ AIRCRAFT_ELEVATION ],
              &variables[ RADIO_ALTITUDE ],
              /* skip GPS_Altitude */
              &variables[ PRESSURE ],
              &variables[ AIR_SPEED ],
              &variables[ GROUND_SPEED ],
              /* Skip Aircraft_Static_Temperature */
              &variables[ TEMPERATURE ],
              &variables[ WIND_U ],
              &variables[ WIND_V ],
              &variables[ OZONE ],
              &variables[ H2O_STATIC_TEMPERATURE ],
              &variables[ RELATIVE_HUMIDITY ],
              &relativeHumidityValidity,
              /* Skip Relative_Humidity_Accuracy */
              &variables[ H2O ],
              &variables[ CO ],
              &variables[ NOY ],
              &variables[ NO ],
              &variables[ NOX ],
              /* Skip NOy_Uncertainty */
              &noValidity ) == 22
    : sscanf( dataLine,
              "%lld %lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
              "%lld %*s %lf %lf %lf %lf %lf %*s %lld\n",
              &yyyymmdd, &hhmmss,
              &variables[ AIRCRAFT_LATITUDE ],
              &variables[ AIRCRAFT_LONGITUDE ],
              &variables[ AIRCRAFT_ELEVATION ],
              &variables[ RADIO_ALTITUDE ],
              &variables[ PRESSURE ],
              &variables[ TEMPERATURE ],
              &variables[ AIR_SPEED ],
              &variables[ GROUND_SPEED ],
              &windDirection,
              &windSpeed,
              &variables[ OZONE ],
              &variables[ H2O_STATIC_TEMPERATURE ],
              &variables[ RELATIVE_HUMIDITY ],
              &relativeHumidityValidity,
              &variables[ H2O ],
              &variables[ CO ],
              &variables[ NOY ],
              &variables[ NO ],
              &variables[ NOX ],
              &noValidity ) == 22;

  if ( ok ) {
    const Integer point = profile ? yyyymmdd : 0; /* Value in first column. */
    const Integer timestamp =
      profile == 0 ? yyyymmdd * 1000000 + hhmmss
      : profileTimestamp( point, points, profile,
                          profileFirstTimestamp, profileLastTimestamp,
                          seconds1, seconds2 );

    variables[ AIRCRAFT_TIMESTAMP ] = timestamp;

#if 0
    /* Don't use profile mid-level altitude. */

    if ( profile ) {
      variables[ AIRCRAFT_ELEVATION ] = hhmmss; /* Value in 2nd column. */
      DEBUG( fprintf( stderr, "%lld ... %lld = "
                      "%lld ... %lld at %lld of %lld (profile = %lld)\n",
                      profileFirstTimestamp, profileLastTimestamp,
                      seconds1, seconds2, point, points, profile ); )
    }
#endif

    DEBUG2( fprintf( stderr, "%lld @ (%lg, %lg, %lg) in %lld ... %lld?\n",
                    timestamp,
                    variables[ AIRCRAFT_LONGITUDE ],
                    variables[ AIRCRAFT_LATITUDE ],
                    variables[ AIRCRAFT_ELEVATION ],
                    firstTimestamp, lastTimestamp ); )

    if ( AND4( isValidYYYYMMDDHHMMSS( timestamp ),
               IN_RANGE( timestamp, firstTimestamp, lastTimestamp ),
               IN_RANGE( variables[ AIRCRAFT_LATITUDE ],
                         bounds[ LATITUDE ][ MINIMUM ],
                         bounds[ LATITUDE ][ MAXIMUM ] ),
               IN_RANGE( variables[ AIRCRAFT_LONGITUDE ],
                         bounds[ LONGITUDE ][ MINIMUM ],
                         bounds[ LONGITUDE ][ MAXIMUM ] ) ) ) {

      /*
       * variables[ AIRCRAFT_ELEVATION ] contains baro_altitude
       * (a pressure-based estimate of height in meters above mean sea level).
       * If this estimated elevation is <= 2,000m then
       * sum radio_altitude (height above ground) + terrain_height
       * (from 2km global file)
       * to get (a more accurate) elevation in meters above mean sea level
       * and use this sum instead of baro_alititude
       * if it differs from baro_altitude by less than 1000m.
       * Note: One cannot reliably use radio_altitude above about 2,000m
       * since the value in the MOZAIC flight file is bogus
       * (approximately repeated from last valid value) during cruising
       * and apparently the radio altimeter instrument has too weak a signal
       * above about 2,500m to be of any use anyway.
       */

      if ( variables[ RADIO_ALTITUDE ] <= 2000.0 ) {
        const Real longitude = variables[ AIRCRAFT_LONGITUDE ];
        const Real latitude  = variables[ AIRCRAFT_LATITUDE ];
        extern float elevationAt( float longitude, float latitude );
        const Real surfaceElevation0 = elevationAt( longitude, latitude );
        const Real surfaceElevation =
          surfaceElevation0 > 0.0 ? surfaceElevation0 : 0.0;
        const Real heightAboveGround = variables[ RADIO_ALTITUDE ];
        const Real elevation = surfaceElevation + heightAboveGround;

        if ( fabs( elevation - variables[ AIRCRAFT_ELEVATION ] ) < 1000.0 ) {
          DEBUG2( fprintf( stderr, "Adjusting elevation: "
                           "At (%lg, %lg) z = %lg -> %lg + %lg = %lg\n",
                           longitude, latitude,
                           variables[ AIRCRAFT_ELEVATION ],
                           surfaceElevation, heightAboveGround, elevation ); )
          variables[ AIRCRAFT_ELEVATION ] = elevation;
        }
      }

      scaleData( variables );
      DEBUG( fprintf( stderr, "variables[ NOY ] = %lf\n", variables[ NOY ] ); )
      result = isValidDataLine( selected, variables,
                                relativeHumidityValidity, noValidity );

      if ( AND2( result, ! isNewFormat ) ) {
        windUV( windDirection, windSpeed,
                variables + WIND_U, variables + WIND_V );
      }

      DEBUG( fprintf( stderr, "in subset, relativeHumidityValidity = %lld, "
                              "noValidity = %lld, valid = %lld\n",
                              relativeHumidityValidity, noValidity, result ); )
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: profileTimestamp - Compute linearly interpolated timestamp of a point
         on a profile.
INPUTS:  Integer point           1-based point number in profile.
         Integer points          Number of points in profile.
         Integer profile         -1 = descending, 1 = ascending.
         Integer firstTimestamp  Beginning (earliest) timestamp of profile.
         Integer lastTimestamp   Ending (latest) timestamp of profile.
         Integer seconds1        Total seconds from Jan 1 yyyy of
                                 firstTimestamp to firstTimestamp.
         Integer seconds2        Total seconds from Jan 1 yyyy of
                                 firstTimestamp to lastTimestamp.
RETURNS: Integer yyyymmddhhmmss of point or 0 if invalid point.
******************************************************************************/

static Integer profileTimestamp( Integer point, Integer points,
                                 Integer profile,
                                 Integer firstTimestamp, Integer lastTimestamp,
                                 Integer seconds1, Integer seconds2 ) {

  PRE08( isValidYYYYMMDDHHMMSS( firstTimestamp ),
         isValidYYYYMMDDHHMMSS( lastTimestamp ),
         firstTimestamp <= lastTimestamp,
         seconds1 >= 0,
         seconds2 >= seconds1,
         IN3( profile, -1, 1 ),
         seconds1,
         seconds2 );

  Integer result = 0;

  if ( IN_RANGE( point, 1, points ) ) {

    if ( point == 1 ) {

      if ( profile == -1 ) {
        result = lastTimestamp;
      } else {
        result = firstTimestamp;
      }
    } else if ( point == points ) {

      if ( profile == -1 ) {
        result = firstTimestamp;
      } else {
        result = lastTimestamp;
      }
    } else {
      const Integer secondsDifference = seconds2 - seconds1;
      const Real interpolation = ( point - 1.0 ) / ( points - 1.0 );
      const Real interpolatedSeconds = secondsDifference * interpolation;
      const Integer targetSeconds =
        profile == 1 ?
          (Integer) ( seconds1 + interpolatedSeconds )
        : (Integer) ( seconds2 - interpolatedSeconds );
      result =
        timestampOfTargetSeconds( firstTimestamp, seconds1, targetSeconds );

      CHECK5( secondsDifference >= 0,
              IN_RANGE( interpolation, 0.0, 1.0 ),
              IN_RANGE( interpolatedSeconds, 0.0, secondsDifference ),
              IN_RANGE( targetSeconds, seconds1, seconds2 ),
              IN_RANGE( result, firstTimestamp, lastTimestamp ) );
    }
  }

  POST0( IMPLIES_ELSE( IN_RANGE( point, 1, points ),
                       AND2( isValidYYYYMMDDHHMMSS( result ),
                             IN_RANGE( result, firstTimestamp, lastTimestamp)),
                       result == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidDataLine - Validate a line of MOZAIC file data.
INPUTS:  const Integer selected[ VARIABLES ]  Selected variables of subset.
         const Real variables[ VARIABLES ]    Data variable values.
         Integer relativeHumidityValidity     RH validity flag.
         Integer noValidity                   NO validity flag.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidDataLine( const Integer selected[ VARIABLES ],
                                const Real variables[ VARIABLES ],
                                Integer relativeHumidityValidity,
                                Integer noValidity ) {

  PRE05( selected,
         minimumItemI( selected, VARIABLES ) >= 0,
         maximumItemI( selected, VARIABLES ) == 1,
         sumI( selected, VARIABLES ) > IMPLICIT_VARIABLES,
         variables );

  const Integer timestamp = variables[ AIRCRAFT_TIMESTAMP ];
  const Integer result =
    AND19( isValidYYYYMMDDHHMMSS( timestamp ),
           isValidLongitude( variables[ AIRCRAFT_LONGITUDE ] ),
           isValidLatitude( variables[ AIRCRAFT_LATITUDE ] ),
           IN_RANGE( variables[ AIRCRAFT_ELEVATION ], -500.0, 1e6 ),
           IMPLIES( selected[ RADIO_ALTITUDE ],
                    IN_RANGE( variables[ RADIO_ALTITUDE ], -500.0, 1e6 ) ),
           IMPLIES( selected[ PRESSURE ],
                    IN_RANGE( variables[ PRESSURE ], 1e3, 1e6 ) ),
           IMPLIES( selected[ TEMPERATURE ],
                    IN_RANGE( variables[ TEMPERATURE ], -90.0, 100.0 ) ),
           IMPLIES( selected[ AIR_SPEED ],
                    IN_RANGE( variables[ AIR_SPEED ], 0.0, 500.0 ) ),
           IMPLIES( selected[ GROUND_SPEED ],
                    IN_RANGE( variables[ GROUND_SPEED ], 0.0, 500.0 ) ),
           IMPLIES( selected[ WIND_U ],
                    IN_RANGE( variables[ WIND_U ], -200.0, 200.0 ) ),
           IMPLIES( selected[ WIND_V ],
                    IN_RANGE( variables[ WIND_V ], -200.0, 200.0 ) ),
           IMPLIES( selected[ OZONE ],
                    IN_RANGE( variables[ OZONE ], 0.0, 10.0 ) ),
           IMPLIES( selected[ H2O_STATIC_TEMPERATURE ],
                    IN_RANGE( variables[ H2O_STATIC_TEMPERATURE ],
                              -90.0, 100.0 ) ),
           IMPLIES( selected[ RELATIVE_HUMIDITY ],
                    AND2( IN_RANGE( variables[ RELATIVE_HUMIDITY ],
                                    0.0, 100.0 ),
                          relativeHumidityValidity != 0 ) ),
           IMPLIES( selected[ H2O ],
                    IN_RANGE( variables[ H2O ], 0.0, 1000.0 ) ),
           IMPLIES( selected[ CO ],
                    IN_RANGE( variables[ CO ], 0.0, 1000.0 ) ),
           IMPLIES( selected[ NOY ],
                    AND2( IN_RANGE( variables[ NOY ], 0.0, 1000.0 ),
                          noValidity != 0 ) ),
           IMPLIES( selected[ NO ],
                    AND2( IN_RANGE( variables[ NO ], 0.0, 1000.0 ),
                          noValidity != 0 ) ),
           IMPLIES( selected[ NOX ],
                    AND2( IN_RANGE( variables[ NOX ], 0.0, 1000.0 ),
                          noValidity != 0 ) ) );

  POST0( IS_BOOL( result ) );

  return result;
}



/******************************************************************************
PURPOSE: timestampOfFileName - Timestamp of MOZAIC file name.
INPUTS:  const char* fileName  Name of MOZAIC file to parse. E.g.,
                               "testdata/M20060703014.txt"
                               "testdata/MA20060703014.txt"
                               "testdata/MD20060703014.txt"
RETURNS: Integer YYYYDDDHHMMSS. E.g., 20060703000000, or 0 if failed and
         failureMessage() is called.
******************************************************************************/

static Integer timestampOfFileName( const char* fileName ) {
  PRE0( fileName );
  Integer result = 0;
  const char* name = strrchr( fileName, '/' );

  /* Skip directory, if present: */

  if ( name ) {
    ++name;
  } else {
    name = fileName;
  }

  /* Parse timestamp: */

  if ( AND3( name, *name == 'M', strlen( name ) > 9 ) ) {
    char timestamp[ 8 + 1 ] = "";
    memset( timestamp, 0, sizeof timestamp );

    if ( isdigit( name[ 1 ] ) ) { /* E.g., M20060703014.txt. */
      strncpy( timestamp, name + 1, 8 ); /* E.g., "20060703". */
    } else if ( IN3( name[ 1 ], 'A', 'D' ) ) { /* E.g., MD20060703014.txt. */
      strncpy( timestamp, name + 2, 8 ); /* E.g., "20060703". */
    }

    CHECK( timestamp[ 8 ] == '\0' );
    result = atoI( timestamp );
    result *= 1000000; /* HHMMSS = 000000. */

    if ( ! isValidYYYYMMDDHHMMSS( result ) ) {
      failureMessage( "Invalid timestamp %s.", timestamp );
      result = 0;
    }
  }

  POST0( IMPLIES( result, isValidYYYYMMDDHHMMSS( result ) ) );
  return result;
}


/******************************************************************************
PURPOSE: parseNote - Parse flight info (line 2) of MOZAIC file.
INPUTS:  const char* fileData  Contents of MOZAIC file to parse.
OUTPUTS: Note note  "MD20060703014:FRANKFURT->ATLANTA".
NOTES:   fileData first two lines to parse should look like:
         LA2006-07-03,KFA2008-01-31,CNRM2005-06-09,Version5.4
         flight: MD20060703014 NL=083 SP=111000 FRANKFURT->ATLANTA
******************************************************************************/

static void parseNote( const char* fileData, Note note ) {
  PRE04( fileData, *fileData, strchr( fileData, '\n' ), note );
  const char* c = strchr( fileData, '\n' );
  assert_static( NOTE_LENGTH > 20 );
  memset( note, 0, sizeof (Note) );

  if ( c ) {

    while ( AND2( *c, *c != ' ' ) ) {
      ++c;
    }

    if ( *c == ' '  ) {
      int index = 0;
      ++c;

      /* Copy flight number, e.g., "MD20060703014": */

      while ( AND2( index < NOTE_LENGTH / 2, isalnum( *c ) ) ) {
        note[ index ] = *c;
        ++index;
        ++c;
      }

      if ( AND2( index < NOTE_LENGTH / 2, *c == ' ' ) ) {
        note[ index ] = ':';
        ++index;
        c = strchr( c + 1, ' '  );

        if ( c ) {
          c = strchr( c + 1, ' '  );

          if ( c ) {
            ++c;

            while ( AND3( index < NOTE_LENGTH, *c, *c != '\n' ) ) {

              if ( isprint( *c ) ) {

                if ( isspace( *c ) ) {
                  note[ index ] = '_';
                } else {
                  note[ index ] = *c;
                }

                ++index;
              }

              ++c;
            }

            while ( index < NOTE_LENGTH ) {
              note[ index ] = ' ';
              ++index;
            }

            note[ NOTE_LENGTH ] = '\0';
          }
        }
      }
    }
  }

  if ( ! AND5( isalnum( *note ),
               note[ NOTE_LENGTH ] == '\0',
               strchr( note, ':' ),
               strstr( note, "->" ),
               strlen( note ) == NOTE_LENGTH ) ) {
    char format[ 80 ] = "";
    sprintf( format, "%%-%ds", NOTE_LENGTH );
    sprintf( note, format, "flight?:from?->to?" );
  }

  POST05( isalnum( *note ), note[ NOTE_LENGTH ] == '\0', strchr( note, ':' ),
          strstr( note, "->" ), strlen( note ) == NOTE_LENGTH );
}



/******************************************************************************
PURPOSE: copySubsetData - Copy subset data into a Track.
INPUTS:  const Real* data               data[ subsetVariables * subsetPoints ].
         const Integer subsetVariables  Number of variables in subset.
         const Integer subsetPoints     Number of points in subset.
         const Integer reverse          Reverse the data line order?
         const Note note                Note on track.
RETURNS: Track* if successful, else 0 if no points are in the subset or
         there was a memory allocation failure and failureMessage() is called.
NOTES:   Use deallocateTrack() then FREE() on returned result when finished.
******************************************************************************/

static Track* copySubsetData( const Real* data,
                              const Integer subsetVariables,
                              const Integer subsetPoints,
                              const Integer reverse,
                              const Note note ) {

  PRE08( data,
         subsetVariables > IMPLICIT_VARIABLES,
         subsetPoints > 0,
         subsetVariables * subsetPoints > 0,
         isNanFree( data, subsetVariables * subsetPoints ),
         note, isalnum( *note ), note[ NOTE_LENGTH ] == '\0' );

  Track* result = NEW_ZERO( Track, 1 );

  if ( result ) {
    result->data = NEW_ZERO( Real, subsetVariables * subsetPoints );

    if ( ! result->data ) {
      FREE( result );
    } else {
      const Integer factor = reverse ? -1 : 1;
      const Real* input = data;
      Real* output =
        result->data + reverse * ( subsetVariables * ( subsetPoints - 1 ) );
      const Integer subsetVariableBytes = subsetVariables * sizeof (Real);
      Integer point = 0;
      const Integer timestamp = input[ AIRCRAFT_TIMESTAMP ];
      result->firstTimestamp = timestamp;
      result->lastTimestamp  = timestamp;
      result->variables = subsetVariables;
      result->points    = subsetPoints;
      result->bounds[ LONGITUDE ][ MINIMUM ] = input[ AIRCRAFT_LONGITUDE ];
      result->bounds[ LONGITUDE ][ MAXIMUM ] = input[ AIRCRAFT_LONGITUDE ];
      result->bounds[ LATITUDE  ][ MINIMUM ] = input[ AIRCRAFT_LATITUDE ];
      result->bounds[ LATITUDE  ][ MAXIMUM ] = input[ AIRCRAFT_LATITUDE ];
      strncpy( result->note, note, NOTE_LENGTH );
      result->note[ NOTE_LENGTH ] = '\0';
      CHECK( isalnum( result->note[ 0 ] ) );

      for ( point = 0; point < subsetPoints;
            ++point, input += subsetVariables,
            output += subsetVariables * factor ) {
        const Integer timestamp = input[ AIRCRAFT_TIMESTAMP ];
        const Real longitude    = input[ AIRCRAFT_LONGITUDE ];
        const Real latitude     = input[ AIRCRAFT_LATITUDE  ];

        if ( timestamp < result->firstTimestamp ) {
          result->firstTimestamp = timestamp;
        } else if ( timestamp > result->lastTimestamp ) {
          result->lastTimestamp = timestamp;
        }

        if ( longitude < result->bounds[ LONGITUDE ][ MINIMUM ] ) {
          result->bounds[ LONGITUDE ][ MINIMUM ] = longitude;
        } else if ( longitude > result->bounds[ LONGITUDE ][ MAXIMUM ] ) {
          result->bounds[ LONGITUDE ][ MAXIMUM ] = longitude;
        }

        if ( latitude < result->bounds[ LATITUDE ][ MINIMUM ] ) {
          result->bounds[ LATITUDE ][ MINIMUM ] = latitude;
        } else if ( latitude > result->bounds[ LATITUDE ][ MAXIMUM ] ) {
          result->bounds[ LATITUDE ][ MAXIMUM ] = latitude;
        }

        CHECK4( IN_RANGE( timestamp,
                          result->firstTimestamp, result->lastTimestamp ),
                IN_RANGE( longitude,
                          result->bounds[ LONGITUDE ][ MINIMUM ],
                          result->bounds[ LONGITUDE ][ MAXIMUM ] ),
                IN_RANGE( latitude,
                          result->bounds[ LATITUDE ][ MINIMUM ],
                          result->bounds[ LATITUDE ][ MAXIMUM ] ),
                IN_RANGE( input[ AIRCRAFT_ELEVATION ], -500.0, 1e6 ) );

        memcpy( output, input, subsetVariableBytes );
      }
    }
  }

  POST0( IMPLIES( result, isValidTrack( result ) ) );

  return result;
}



/******************************************************************************
PURPOSE: totalSubsetPoints - Count the total number of subset points over all
         tracks.
INPUTS:  const VoidList* tracks   List of subsetted tracks.
RETURNS: Integer number of points.
******************************************************************************/

static Integer totalSubsetPoints( const VoidList* tracks ) {

  PRE02( tracks, tracks->invariant( tracks ) );

  const Integer count = tracks->count( tracks );
  Integer result = 0;
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {
    const Track* const track = (const Track*) tracks->item( tracks, index );
    result += track->points;
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: writeData - Write subsetted track data to stdout.
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

  PRE06( isValidData( data ), data->ok, output, output->invariant( output ),
         output->isWritable( output ), data->ok == output->ok( output ) );

  const Arguments* const arguments = &data->arguments;
  const VoidList*  const tracks = data->tracks;
  const Track* const firstTrack = tracks->item( tracks, 0 );
  const Integer variables = firstTrack->variables;
  const Integer totalPoints = totalSubsetPoints( tracks );
  Integer variable = 0;
  UTCTimestamp firstTimestamp;
  UTCTimestamp lastTimestamp;
  toUTCTimestamp2( arguments->firstTimestamp, firstTimestamp );
  toUTCTimestamp2( arguments->lastTimestamp,  lastTimestamp );
  CHECK( variables > 0 );

  output->writeString( output,
                       "Aircraft 2.0\n%s\n%s %s\n"
                       "# Subset domain:"
                       " <min_lon> <min_lat> <max_lon> <max_lat>:\n"
                       "%g %g %g %g\n"
                       "# Dimensions: variables points tracks:\n"
                       "%lld %lld %lld\n"
                       "# Variable names:\n",
                       arguments->description, firstTimestamp, lastTimestamp,
                       arguments->bounds[ LONGITUDE ][ MINIMUM ],
                       arguments->bounds[ LATITUDE  ][ MINIMUM ],
                       arguments->bounds[ LONGITUDE ][ MAXIMUM ],
                       arguments->bounds[ LATITUDE  ][ MAXIMUM ],
                       variables, totalPoints, tracks->count( tracks ) );

  for ( variable = 0;
        AND2( output->ok( output ), variable < VARIABLES );
        ++variable ) {

    if ( arguments->selected[ variable ] ) {
      const char* const variableName = variableNames[ variable ];
      CHECK2( variableName, variableName[ 0 ] );
      output->writeString( output, " %s", variableName );
    }
  }

  if ( output->ok( output ) ) {
    output->writeString( output, "\n# Variable units:\n" );
  }

  for ( variable = 0;
        AND2( output->ok( output ), variable < VARIABLES );
        ++variable ) {

    if ( arguments->selected[ variable ] ) {
      const char* const units = variableUnits[ variable ];
      CHECK2( units, units[ 0 ] );
      output->writeString( output, " %s", units );
    }
  }

  if ( output->ok( output ) ) {
    output->writeString( output,
                         "\n# char notes[tracks][%d] and\n"
                         "# IEEE-754 64-bit reals "
                         "bounds[tracks][2=lon,lat][2=min,max] and\n"
                         "# MSB 64-bit integers points[tracks] and\n"
                         "# IEEE-754 64-bit reals"
                         " data_1[points_1][variables] ..."
                         " data_T[points_T][variables]:\n",
                         NOTE_LENGTH + 1,
                         arguments->bounds[ LONGITUDE ][ MINIMUM ],
                         arguments->bounds[ LATITUDE  ][ MINIMUM ],
                         arguments->bounds[ LONGITUDE ][ MAXIMUM ],
                         arguments->bounds[ LATITUDE  ][ MAXIMUM ] );
  }

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

  const VoidList* tracks = data->tracks;
  writeTrackNotes( tracks, output );
  data->ok = output->ok( output );

  if ( data->ok ) {
    writeTrackBounds( tracks, output );
    data->ok = output->ok( output );

    if ( data->ok ) {
      writeTrackPoints( tracks, output );
      data->ok = output->ok( output );

      if ( data->ok ) {
        writeTrackData( tracks, output );
        data->ok = output->ok( output );
      }
    }
  }

  POST04( isValidData( data ), output->invariant( output ),
          output->isWritable( output ), data->ok == output->ok( output ) );
}



/******************************************************************************
PURPOSE: writeTrackNotes - Write track notes.
INPUTS:  const VoidList* tracks   List of Track*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeTrackNotes( const VoidList* tracks, Stream* output ) {

  PRE07( tracks, tracks->invariant( tracks ), tracks->count( tracks ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer trackCount = tracks->count( tracks );
  Integer trackIndex = 0;

  do {
    const Track* const track = tracks->item( tracks, trackIndex );
    char format[ 80 ] = "";
    memset( format, 0, sizeof format );
    sprintf( format, "%%-%ds\n", NOTE_LENGTH );
    CHECK( format[ sizeof format / sizeof *format - 1 ] == '\0' );
    CHECK( isValidTrack( track ) );
    output->writeString( output, format, track->note );
    ++trackIndex;
  } while ( AND2( output->ok( output ), trackIndex < trackCount ) );

  POST04( tracks->invariant( tracks ), tracks->count( tracks ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: writeTrackBounds - Write MSB 64-bit integer track bounds.
INPUTS:  const VoidList* tracks   List of Track*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeTrackBounds( const VoidList* tracks, Stream* output ) {

  PRE07( tracks, tracks->invariant( tracks ), tracks->count( tracks ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer trackCount = tracks->count( tracks );
  Integer trackIndex = 0;

  do {
    const Track* const track = tracks->item( tracks, trackIndex );
    CHECK( isValidTrack( track ) );
    output->write64BitReals( output, &track->bounds[0][0], 4 );
    ++trackIndex;
  } while ( AND2( output->ok( output ), trackIndex < trackCount ) );

  POST04( tracks->invariant( tracks ), tracks->count( tracks ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: writeTrackPoints - Write MSB 64-bit integer track subset point counts.
INPUTS:  const VoidList* tracks   List of Track*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeTrackPoints( const VoidList* tracks, Stream* output ) {

  PRE07( tracks, tracks->invariant( tracks ), tracks->count( tracks ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer trackCount = tracks->count( tracks );
  Integer trackIndex = 0;

  do {
    const Track* const track = tracks->item( tracks, trackIndex );
    CHECK( isValidTrack( track ) );
    output->write64BitInteger( output, track->points );
    ++trackIndex;
  } while ( AND2( output->ok( output ), trackIndex < trackCount ) );

  POST04( tracks->invariant( tracks ), tracks->count( tracks ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: writeTrackData - Write 64-bit IEEE-754 track subset variable data.
INPUTS:  const VoidList* tracks   List of Track*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeTrackData( const VoidList* tracks, Stream* output ) {

  PRE07( tracks, tracks->invariant( tracks ), tracks->count( tracks ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer trackCount = tracks->count( tracks );
  Integer trackIndex = 0;

  do {
    const Track* const track = tracks->item( tracks, trackIndex );
    CHECK( isValidTrack( track ) );
    output->write64BitReals( output, track->data,
                             track->variables * track->points );
    ++trackIndex;
  } while ( AND2( output->ok( output ), trackIndex < trackCount ) );

  POST04( tracks->invariant( tracks ), tracks->count( tracks ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}




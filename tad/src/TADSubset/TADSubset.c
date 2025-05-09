
/******************************************************************************
PURPOSE: TADSubset.c - Read a set of a TAD aircraft measurement files,
         subset the data to a bounds (longitude-latitude rectangle) and
         write it to stdout as XDR (IEEE-754) format binary.

NOTES:   Uses libUtilities.a (../../libs/Utilities).
         https://tad.larc.nasa.gov/


INPUTS: data files look like this
        time-sorted tab-delimited sample 405 data line file with 3-line header:

INTEX-B Flight 15-20 samping over pacific						
405						
Timestamp(UTC_Start)	LATITUDE(degN)	LONGITUDE(degE)	PressureAtltitude(km)	CO(ppbv)		
2006-04-21T17:07:30Z	48.0132	-122.332	0.890735	168.37		
2006-04-21T17:08:30Z	47.9626	-122.376	1.05256	235.787		
...
2006-04-22T00:46:30Z	47.8575	-122.31	0.461381	149.125		
2006-04-22T00:47:30Z	47.8835	-122.282	0.210875	153.744		

OUTPUTS: XDR (portable binary) format file with 15-line ASCII header:

Aircraft 2.0
https://tad.larc.nasa.gov/,TADSubset
2006-04-21T00:00:00-0000 2006-04-22T23:59:59-0000
# Subset domain: <min_lon> <min_lat> <max_lon> <max_lat>:
-125 45 -120 50
# Dimensions: variables points tracks:
4 394 1
# Variable names:
timestamp longitude latitude elevation co
# Variable units:
yyyymmddhhmmss deg deg m ppbv
# char notes[tracks][80] and
# IEEE-754 64-bit reals bounds[tracks][2=lon,lat][2=min,max] and
# MSB 64-bit integers points[tracks] and
# IEEE-754 64-bit reals data_1[points_1][variables] ... data_T[points_T][variables]:
INTEX-B Flight 15-20 samping over pacific
-1.220000000000000e+01
-1.230000000000000e+00
4.7100000000000000e+01
4.7900000000000000e+01
394
2.0060421170730000e+13
-1.223320000000000e+01
4.8013200000000000e+01
8.9073500000000000e+01
1.6837000000000000e+03
2.0060421170830000e+13
-1.223760000000000e+01
4.7962000000000000e+01
1.0525600000000000e+01
2.3578700000000000e+03
...
2.0060422004630000e+13
-1.223100000000000e+01
4.7857500000000000e+01
4.6138100000000000e+01
1.4912500000000000e+03
2.0060422004730000e+13
-1.222820000000000e+01
4.7883500000000000e+01
2.1087500000000000e+01
1.5374400000000000e+03

HISTORY: 2015-02-19 plessel.todd@epa.gov Created
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>     /* For printf(), snprintf(). */
#include <string.h>    /* For strlen(). */
#include <ctype.h>     /* For isdigit(), isalnum(), isspace(), isprint(). */

#include <Utilities.h> /* For PRE0*(), NEW_ZERO(), Stream, VoidList, Note. */

/*================================== TYPES ==================================*/

enum { VARIABLES = 5 };/* timestamp, longitude, latitude, elevation, ozone */

static const double MISSING = -9999.0;

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;       /* File listing TAD files to read. */
  const char* description;    /* User-supplied description. */
  const char* variable;       /* User-supplied variable name. */
  Note units;                 /* Units of variable (from TAD file header). */
  Bounds      bounds;         /* bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM].*/
  long long   firstTimestamp; /* YYYYMMDDHHMMSS of subset. */
  long long   lastTimestamp;  /* YYYYMMDDHHMMSS of subset. */
} Arguments;

#ifndef NO_ASSERTIONS

/* Arguments invariant: */

static int isValidArguments( const Arguments* arguments ) {
  const int result =
    AND13( arguments,
           arguments->listFile,
           arguments->listFile[ 0 ],
           arguments->description,
           arguments->description[ 0 ],
           arguments->variable,
           isalpha( arguments->variable[ 0 ] ),
           strlen( arguments->variable ) < 32,
           arguments->units,
           isValidBounds( arguments->bounds ),
           isValidYYYYMMDDHHMMSS( arguments->firstTimestamp ),
           isValidYYYYMMDDHHMMSS( arguments->lastTimestamp ),
           arguments->firstTimestamp <= arguments->lastTimestamp );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Track: Result of reading a subset of a TAD aircraft data file. */

typedef struct {
  Note note;     /* "INTEX-B Flight 15-20 samping over pacific". */
  Bounds bounds; /* bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  int points;    /* Number of data points. */
  double* data;  /* data[ 5 ][ points ]. */
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

static int isValidTrack( const Track* track ) {
  const int result =
    AND7( track,
          track->note[ 0 ],
          track->note[ NOTE_LENGTH ] == '\0',
          strlen( track->note ) == NOTE_LENGTH,
          isValidBounds( track->bounds ),
          track->points >= 1,
          isNanFree( track->data, 5 * track->points ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Data type: */

typedef struct {
  Arguments arguments; /* User-supplied (command-line) arguments. */
  VoidList* tracks;    /* List of subsetted tracks. */
  int ok;              /* Did last command succeed? */
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

static int isValidData( const Data* data ) {
  int result =
    AND5( data,
          isValidArguments( &data->arguments ),
          data->tracks,
          data->tracks->invariant( data->tracks ),
          IS_BOOL( data->ok ) );

  if ( result ) {
    const VoidList* const tracks = data->tracks;
    const int trackCount = tracks->count( tracks );
    int index = 0;

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

static int parseArguments( int argc, char* argv[], Arguments* arguments );

static void initializeArguments( Arguments* arguments );

static void readData( Data* data );

static Track* readTADFile( const char* fileName, Arguments* arguments );

static void updateBounds( const double values[ 5 ], int* initialized,
                          Bounds bounds );

static void updateRange( const double value,
                         double* minimum, double* maximum );

static const char* parseHeaderLines( const char* fileData,
                                     const char* variable,
                                     int* points,
                                     int* longitudeColumn, int* variableColumn,
                                     Note note, Note units );

static void parseHeaderNote( const char* fileData, Note note );

static char readDelimiter( const char* dataLine );

static int headerColumn( const char* headerLine, const char* variable,
                         Note units );

static int parseDataLine( const char* dataLine,
                          const char delimiter,
                          const long long firstTimestamp,
                          const long long lastTimestamp,
                          const Bounds bounds,
                          const int longitudeColumn,
                          const int variableColumn,
                          const Note units,
                          double dataValues[ 5 ] );

static double parseColumnValue( const char* dataLine, const char delimiter,
                                const int column );

static int totalSubsetPoints( const VoidList* tracks );

static void writeData( Data* data );

static void writeHeader( Data* data, Stream* output );

static void writeXDR( Data* data, Stream* output );

static void writeTrackNotes( const VoidList* tracks, Stream* output );

static void writeTrackBounds( const VoidList* tracks, Stream* output );

static void writeTrackPoints( const VoidList* tracks, Stream* output );

static void writeTrackData( const VoidList* tracks, Stream* output );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Read a subset of a list TAD files and write it to stdout in
         XDR format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int ok = 0;

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
  fprintf( stderr,
           "\n\n%s - Read a set of TAD files and extract track\n",
           programName );
  fprintf( stderr, "data for selected variables" );
  fprintf( stderr, " subsetted by a lon-lat rectangle.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "  -files <listFile> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -time <yyyymmddhhmmss> <yyyymmddhhmmss> \\\n" );
  fprintf( stderr, "  -variable <name>\\n" );
  fprintf( stderr, "  [ -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ] \\\n" );
  fprintf( stderr, "\n\n" );
  fprintf( stderr, "Note: times are in UTC (GMT)\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files /data/files.txt \\\n" );
  fprintf( stderr, "-desc https://tad.larc.nasa.gov/,TADSubset \\\n" );
  fprintf( stderr, "-time 20060421000000 20060422235959 \\\n" );
  fprintf( stderr, "-variable co \\\n" );
  fprintf( stderr, "-domain -125 45 -120 50 > subset.xdr\n\n" );
  fprintf( stderr, "Subset of data for April 21, 2006 ovwer the Pacific\n");
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays\n" );
  fprintf( stderr, "For example:\n" );
  fprintf( stderr, "AIRCRAFT 2.0\n" );
  fprintf( stderr, "https://tad.larc.nasa.gov/,TADSubset\n" );
  fprintf( stderr, "2006-04-21T00:00:00-0000 2006-04-22T23:59:59-0000\n" );
  fprintf( stderr, "# Subset domain: <min_lon> <min_lat> <max_lon> <max_lat>");
  fprintf( stderr, ":\n");
  fprintf( stderr, "-125 45 -120 50\n" );
  fprintf( stderr, "# Dimensions: variables points tracks:\n" );
  fprintf( stderr, "4 48 2\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "timestamp longitude latitude elevation co\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss deg deg m ppbv\n" );
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
INPUTS:  int argc              Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: int 1 if successful, else 0 and failureMessage() then printUsage()
         are called.
******************************************************************************/

static int parseArguments( int argc, char* argv[],
                           Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  int result = 0;
  ZERO_OBJECT( arguments );
  initializeArguments( arguments );

  if ( ! IN3( argc, 10, 15 ) ) {
    failureMessage( "Invalid/insufficient command line arguments." );
  } else {
    long long arg = 1; /* Argument to parse next, updated by parseArgument2()*/
    arguments->listFile = parseArgument2( argc, argv, "-files", &arg );

    if ( arguments->listFile ) {
      arguments->description = parseArgument2( argc, argv, "-desc", &arg );

      if ( arguments->description ) {

        if ( AND2( ! strcmp( argv[ arg ], "-time" ),
                   parseTimeRange( argv[ arg + 1 ], argv[ arg + 2 ],
                                   &arguments->firstTimestamp,
                                   &arguments->lastTimestamp ) ) ) {
          arg += 3;
          arguments->variable = parseArgument2( argc, argv, "-variable", &arg);

          if ( arguments->variable ) {
            result = 1;

            if ( AND2( arg < argc, ! strcmp( argv[ arg ], "-domain" ) ) ) {
              result = parseBounds( argc, argv, &arg, arguments->bounds );
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
  arguments->units[ 0 ] = '-'; /* Placeholder until input file read. */
}



/******************************************************************************
PURPOSE: readData - Read data from TAD files and subset it by time,
         lon-lat box and selected variables.
INPUTS:  Data* data  data->arguments->listfile is the file containing the
         names of TAD files to read and subset.
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

    /* For each file, read a subset of it into track and append to list: */

    do {
      FileName fileName = "";
      listFile->readWord( listFile, fileName,
                          sizeof fileName / sizeof *fileName );
      data->ok = listFile->ok( listFile );

      if ( data->ok ) {
        Track* track = readTADFile( fileName, &data->arguments );

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
PURPOSE: readTADFile - Read a subset of track data from TAD file.
INPUTS:  const char* fileName  Name of compressed TAD file to read.
         Arguments* arguments  Input arguments.
OUTPUTS: Arguments* arguments  Updated arguments->units.
RETURNS: Track* track  Track of subset data or 0 if no data in subset.
NOTES:   If unsuccessful then failureMessage() is called and 0 is returned.
******************************************************************************/

static Track* readTADFile( const char* fileName, Arguments* arguments ) {

  PRE02( fileName, isValidArguments( arguments ) );

  Track* result = 0;
  long long length = 0;
  char* fileData = readFile( fileName, &length );

  if ( fileData ) {
    int points = 0;
    int variableColumn = 0;
    int longitudeColumn = 0;
    Note units = "";
    Note note = "";
    const char* dataLine = 0;
    memset( note, 0, sizeof note );
    memset( units, 0, sizeof units );

    dataLine =
      parseHeaderLines( fileData, arguments->variable,
                        &points, &longitudeColumn, &variableColumn,
                        note, units );
    DEBUG( fprintf( stderr, "Reading TAD file %s, points = %d, "
                    "longitudeColumn = %d, variableColumn = %d, "
                    "note = '%s', units = '%s'\n",
                    fileName, points,
                    longitudeColumn, variableColumn, note, units ); )

    if ( *units ) {
      memcpy( arguments->units, units, sizeof units );
    }

    if ( points > 0 ) {
      const char delimiter = readDelimiter( dataLine );

      if ( delimiter ) {
        double* data = NEW_ZERO( double, VARIABLES * points );

        if ( data ) {
          int initialized = 0;
          Bounds bounds;
          double* output = data;
          int subsetPoints = 0;

          while ( dataLine ) {
            double values[ VARIABLES ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

            if ( parseDataLine( dataLine, delimiter,
                                arguments->firstTimestamp,
                                arguments->lastTimestamp,
                                (const double (*)[2]) arguments->bounds,
                                longitudeColumn, variableColumn, units,
                                values ) ) {
              int variable = 0;
              updateBounds( values, &initialized, bounds );

              for ( variable = 0; variable < VARIABLES; ++variable ) {
                *output++ = values[ variable ];
              }

              ++subsetPoints;
            }

            dataLine = skipLines( dataLine, 1 );
          }

          if ( subsetPoints ) {
            result = NEW_ZERO( Track, 1 );

            if ( result ) {
              strncpy( result->note, note, NOTE_LENGTH );
              result->note[ NOTE_LENGTH ] = '\0';
              memcpy( result->bounds, bounds, sizeof bounds );
              result->points = subsetPoints;
    
              if ( subsetPoints < points / 2 ) { /* Copy data if small subset*/
                result->data = NEW_ZERO( double, VARIABLES * subsetPoints );
              }

              if ( result->data ) {
                memcpy( result->data, data,
                        VARIABLES * subsetPoints * sizeof *data );
              } else { /* Transfer ownership of data: */
                result->data = data;
                data = 0;
              }
            }
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
PURPOSE: updateBounds - Update longitude latitude bounds.
INPUTS:  const double values[ 5 ] Current values parsed.
         int* initialized         Has bounds been initialized?
         Bounds bounds            Current data bounds.
OUTPUTS: int* initialized         Updated.
         Bounds bounds            Updated data bounds.
******************************************************************************/

static void updateBounds( const double values[ 5 ], int* initialized,
                          Bounds bounds ) {

  PRE03( values, initialized, bounds );
  const double longitude = values[ 1 ];
  const double latitude  = values[ 2 ];

  if ( *initialized == 0 ) {
    bounds[ LONGITUDE ][ MINIMUM ] = longitude;
    bounds[ LONGITUDE ][ MAXIMUM ] = longitude;
    bounds[ LATITUDE  ][ MINIMUM ] = latitude;
    bounds[ LATITUDE  ][ MAXIMUM ] = latitude;
    *initialized = 1;
  } else {
    updateRange( longitude,
                 &bounds[ LONGITUDE ][ MINIMUM ],
                 &bounds[ LONGITUDE ][ MAXIMUM ] );
    updateRange( latitude,
                 &bounds[ LATITUDE ][ MINIMUM ],
                 &bounds[ LATITUDE ][ MAXIMUM ] );
  }

  POST03( *initialized == 1,
          IN_RANGE( values[ 1 ],
                   bounds[ LONGITUDE ][ MINIMUM ],
                   bounds[ LONGITUDE ][ MAXIMUM ] ),
          IN_RANGE( values[ 2 ],
                   bounds[ LATITUDE ][ MINIMUM ],
                   bounds[ LATITUDE ][ MAXIMUM ] ) );
}



/******************************************************************************
PURPOSE: updateRange - Update data range.
INPUTS:  const double value     Current value.
         const double* minimum  Current minimum value.
         const double* maximum  Current maximum value.
OUTPUTS: const double* minimum  Updated minimum value.
         const double* maximum  Updated maximum value.
******************************************************************************/

static void updateRange( const double value,
                         double* minimum, double* maximum ) {

  PRE04( value >= MISSING, minimum, maximum, *minimum <= *maximum );

  if ( value < *minimum ) {
    *minimum = value;
  } else if ( value > *maximum ) {
    *maximum = value;
  }

  POST0( IN_RANGE( value, *minimum, *maximum ) );
}



/******************************************************************************
PURPOSE: parseHeaderLines - Parse header lines of TAD file data.
INPUTS:  const char* fileData   All lines from TAD file.
         const char* variable   Name of variable.
OUTPUTS: int* points            Number of data lines.
         int* longitudeColumn   0-based column containing longitude.
         int* variableColumn    0-based column containing variable.
         Note note              Header note.
         Note units             Units of variable.
RETURNS: char* data lines from TAD file or 0 if invalid fileData.
******************************************************************************/

static const char* parseHeaderLines( const char* fileData,
                                     const char* variable,
                                     int* points,
                                     int* longitudeColumn, int* variableColumn,
                                     Note note, Note units ) {

  PRE08( fileData, variable, isalpha( *variable ),
         points, longitudeColumn, variableColumn, note, units );

  const int headerLines = 3;
  const char* result = skipLines( fileData, headerLines );

  *points = *longitudeColumn = *variableColumn = 0;
  note[ 0 ] = units[ 0 ] = '\0';
  parseHeaderNote( fileData, note );

  if ( AND3( result,
             sscanf( fileData, "%*[^\n]\n%d\n", points ) == 1,
             *points > 0 ) ) {
    const char* headerLine = skipLines( fileData, 2 );

    if ( headerLine ) {
      *longitudeColumn = headerColumn( headerLine, "longitude", 0 );

      if ( *longitudeColumn > 0 ) {
        *variableColumn = headerColumn( headerLine, variable, units );
      }
    }
  }

  if ( *variableColumn == -1 ) {
    *variableColumn = *longitudeColumn = 0;
    *points = 0;
    result = 0;
    memset( note, 0, sizeof (Note ) );
    memset( units, 0, sizeof (Note ) );
  }

  POST0( IMPLIES_ELSE( result,
                       AND5( *points > 0,
                             *longitudeColumn > 0,
                             *variableColumn > *longitudeColumn,
                             strlen( note ) == NOTE_LENGTH,
                             *units ),
                       AND5( *points == 0,
                             *longitudeColumn == 0,
                             *variableColumn == 0,
                             note[ 0 ] == '\0',
                             units[ 0 ] == '\0' ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseHeaderNote - Parse note on first line of TAD file.
INPUTS:  const char* fileData   All lines from TAD file.
OUTPUTS: Note note              Header note.
******************************************************************************/

static void parseHeaderNote( const char* fileData, Note note ) {

  PRE02( fileData, note );

  const char* input = fileData;
  int index = 0;
  int done = 0;

  note[ 0 ] = '\0';

  while ( AND2( ! done, index < NOTE_LENGTH ) ) {
    char ch = *input++;
    done = OR2( IN4( ch, '\0', '\n', '\r' ), ! isprint( ch ) );

    if ( ! done ) {

      if ( ch == '\t' ) { /* Replace tabs with spaces: */
        ch = ' ';
      }

      /* Skip leading delimiters: */

      if ( IMPLIES( IN3( ch, ' ', ',' ), index > 0 ) ) {
        note[ index ] = ch;
        ++index;
      }
    }
  }

  for ( ; index < NOTE_LENGTH; ++index ) {
    note[ index ] = ' ';
  }

  /* Trim trailing commas: */

  for ( index = NOTE_LENGTH - 1; index > 0; --index ) {

    if ( AND2( note[ index ] == ',', note[ index - 1 ] == ',' ) ) {
      note[ index ] = ' ';
    }
  }

  if ( isspace( *note ) ) { /* Replace empty note with "_": */
    *note = '_';
  }

  POST02( strlen( note ) == NOTE_LENGTH, ! isspace( *note ) );
}



/******************************************************************************
PURPOSE: readDelimiter - Read delimiter from first data line of TAD file.
INPUTS:  const char* dataLine   Data line from TAD file.
RETURNS: char delimiter separating columns or 0 if invalid and failure printed.
******************************************************************************/

static char readDelimiter( const char* dataLine ) {

  PRE0( dataLine );

  char result = 0;
  const char* p = dataLine;

  do {
    const char c = *p++;

    if ( IN6( c, '\t', ',', '\n', '\r', '\0' ) ) {
      result = c;
    }

  } while ( result == 0 );

  if ( ! IN3( result, '\t', ',' ) ) {
    failureMessage( "Invalid/missing TAD input file data column delimiters." );
    result = 0;
  }

  POST0( IMPLIES( result, IN3( result, '\t', ',' ) ) );
  return result;
}



/******************************************************************************
PURPOSE: headerColumn - Return the 0-based column number for variable and
         optionally its units.
INPUTS:  const char* headerLine  Header line from TAD file.
         const char* variable    Name of variable.
OUTPUTS: Note units              0 or units of variable.
RETURNS: int 0-based column number of variable or -1 if not found.
NOTES:   Will match to either all uppercase or all lowercase variable.
******************************************************************************/

static int headerColumn( const char* headerLine, const char* variable,
                         Note units ) {

  PRE03( headerLine, variable, isalpha( *variable ) );

  int result = -1;
  char delimiter = '\t';
  int ok = strchr( headerLine, delimiter ) != 0;

  if ( ! ok ) {
    delimiter = ',';
    ok = strchr( headerLine, delimiter ) != 0;
  }

  if ( ok ) {
    const char* c = 0;
    char word[ 80 ] = "";
    memset( word, 0, sizeof word );
    snprintf( word, sizeof word / sizeof *word - 1,
             "%c%s(", delimiter, variable );
    uppercase( word );
    c = strstr( headerLine, word );

    if ( ! c ) {
      lowercase( word );
      c = strstr( headerLine, word );
      ok = c != 0;
    }

    if ( ok ) {
      int column = 0;
      const char* p = 0;

      for ( p = headerLine; p <= c; ++p ) {

        if ( *p == delimiter ) {
          ++column;
        }
      }

      result = column;

      if ( units ) { /* Copy units of variable. */
        int index = 0;
        ok = 1;
        c += strlen( word ); /* Skip to inside parentheses. */
        memset( units, 0, sizeof (Note) );

        for ( index = 0; index < NOTE_LENGTH; ++index ) {
          char ch = ok ? c[ index ] : '\0';
          ok = AND3( ok, isprint( ch ), ! IN5( ch, '\0', ')', '\r', '\n' ) );

          if ( AND2( ok, IN3( ch, ' ', '\t' ) ) ) {
            ch = '_'; /* Replace spaces in units name with underscores. */
          }

          if ( ok ) {
            units[ index ] = ch;
          }
        }

        if ( units[ 0 ] == '\0' ) {
          result = 0;
        }
      }
    }
  }

  if ( result == -1 ) {
    failureMessage( "Variable %s not available in TAD input file.", variable );

    if ( units ) {
      memset( units, 0, sizeof (Note) );
    }
  }

  POST0( IMPLIES_ELSE( result >= 0,
                       IMPLIES( units, *units ),
                       IMPLIES( units, units[ 0 ] == '\0' ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseDataLine - Parse, validate and subset a line of TAD file data.
INPUTS:  const char* dataLine            Line of data from TAD file.
         const char delimiter            Delimiter character, e.g., '\t', ','.
         const long long firstTimestamp  Beginning timestamp of subset.
         const long long lastTimestamp   Ending timestamp of subset.
         const Bounds bounds             Lon-lat bounds of subset.
         const int longitudeColumn       0-based column number of longitude.
         const int variableColumn        0-based column number of variable.
         const Note units                Units of variable.
OUTPUTS: double dataValues[ 5 ]          Data variable values if in subset.
RETURNS: int 1 if valid and in subset, else 0.
******************************************************************************/

static int parseDataLine( const char* dataLine,
                          const char delimiter,
                          const long long firstTimestamp,
                          const long long lastTimestamp,
                          const Bounds bounds,
                          const int longitudeColumn,
                          const int variableColumn,
                          const Note units,
                          double dataValues[ 5 ] ) {

  PRE09( dataLine, delimiter,
         isValidYYYYMMDDHHMMSS( firstTimestamp ),
         isValidYYYYMMDDHHMMSS( lastTimestamp ),
         firstTimestamp <= lastTimestamp,
         isValidBounds( bounds ),
         longitudeColumn > 0,
         variableColumn > longitudeColumn,
         dataValues );

  int result = 0;
  UTCTimestamp utcTimestamp = "";
  memset( utcTimestamp, 0, sizeof utcTimestamp );
  result = sscanf( dataLine, "%20s", utcTimestamp ) == 1;

  if ( result ) { /* Convert 'Z' suffix to "-0000" */
    utcTimestamp[ 19 ] = '-';
    utcTimestamp[ 20 ] = '0';
    utcTimestamp[ 21 ] = '0';
    utcTimestamp[ 22 ] = '0';
    utcTimestamp[ 23 ] = '0';
    result = isValidUTCTimestamp( utcTimestamp );

    if ( result ) {
      const long long yyyydddhhmm = fromUTCTimestamp( utcTimestamp );
      const long long timestamp = toYYYYMMDDHHMMSS( yyyydddhhmm );
      result = IN_RANGE( timestamp, firstTimestamp, lastTimestamp );

      if ( result ) {
        const int latitudeColumn =
          longitudeColumn == 1 ? 2 : longitudeColumn - 1;
        const int elevationColumn = 3;
        double value = 0.0;
        dataValues[ 0 ] = timestamp;
        value = parseColumnValue( dataLine, delimiter, longitudeColumn );
        result = IN_RANGE( value,
                           bounds[ LONGITUDE ][ MINIMUM ],
                           bounds[ LONGITUDE ][ MAXIMUM ] );

        if ( result ) {
          dataValues[ 1 ] = value;
          value = parseColumnValue( dataLine, delimiter, latitudeColumn );
          result = IN_RANGE( value,
                             bounds[ LATITUDE ][ MINIMUM ],
                             bounds[ LATITUDE ][ MAXIMUM ] );

          if ( result ) {
            const double minimumValidElevation = -500.0; /* Meters. */
            const double maximumValidElevation = 1e6;    /* Meters. */
            const double km_to_m = 1000.0;
            dataValues[ 2 ] = value;
            value = parseColumnValue( dataLine, delimiter, elevationColumn );
            value *= km_to_m;
            result = IN_RANGE( value,
                               minimumValidElevation, maximumValidElevation );

            if ( result ) {
              double minimumValidValue = 0.0;
              double maximumValidValue = 1e6;
              dataValues[ 3 ] = value;

              if ( *units == '%' ) {
                maximumValidValue = 100.0;
              } else if ( ! strcmp( units, "m/s" ) ) {
                minimumValidValue = -500.0;
                maximumValidValue =  500.0;
              } else if ( OR3( ! strcmp( units, "C" ),
                               ! strcmp( units, "degC" ),
                               ! strcmp( units, "degreesC" ) ) ) {
                minimumValidValue = -100.0;
                maximumValidValue =  100.0;
              }

              value = parseColumnValue( dataLine, delimiter, variableColumn );
              result = IN_RANGE( value, minimumValidValue, maximumValidValue );

              if ( result ) {
                dataValues[ 4 ] = value;
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    memset( dataValues, 0, VARIABLES * sizeof *dataValues );
  }

  POST02( IS_BOOL( result ),
         IMPLIES_ELSE( result,
                       AND5( dataValues[ 0 ],
                             IN_RANGE( dataValues[ 1 ],
                                       bounds[ LONGITUDE ][ MINIMUM ],
                                       bounds[ LONGITUDE ][ MAXIMUM ] ),
                             IN_RANGE( dataValues[ 2 ],
                                       bounds[ LATITUDE ][ MINIMUM ],
                                       bounds[ LATITUDE ][ MAXIMUM ] ),
                             IN_RANGE( dataValues[ 3 ], -500.0, 1e6 ),
                             dataValues[ 4 ] > MISSING ),
                      IS_ZERO5( dataValues[ 0 ], dataValues[ 1 ],
                                dataValues[ 2 ], dataValues[ 3 ],
                                dataValues[ 4 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseColumnValue - Parse data value at given column.
INPUTS:  const char* dataLine  Data line to parse.
         const char delimiter  Delimiter between column values. '\t' or ','.
         const int column      0-based column number to parse.
RETURNS: double value of data at column or MISSING if invalid.
******************************************************************************/

static double parseColumnValue( const char* dataLine, const char delimiter,
                                const int column ) {

  PRE03( dataLine, delimiter, column > 0 );
  double result = MISSING;
  const char* word = 0;
  int skippedColumns = 1;

  for ( word = strchr( dataLine, delimiter );
        AND3( word, *word, skippedColumns < column );
        ++skippedColumns ) {
    word = strchr( word + 1, delimiter );
  }

  if ( AND3( skippedColumns == column, word, *word == delimiter ) ) {
    const int ok =
      AND3( sscanf( word + 1, "%lf", &result ) == 1,
            ! isNan( result ),
           result >= MISSING );

    if ( ! ok ) {
      result = MISSING;
    }
  }

  POST0( result >= MISSING );
  return result;
}



/******************************************************************************
PURPOSE: totalSubsetPoints - Count the total number of subset points over all
         tracks.
INPUTS:  const VoidList* tracks   List of subsetted tracks.
RETURNS: int number of points.
******************************************************************************/

static int totalSubsetPoints( const VoidList* tracks ) {

  PRE02( tracks, tracks->invariant( tracks ) );

  const int count = tracks->count( tracks );
  int result = 0;
  int index = 0;

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
  const int totalPoints = totalSubsetPoints( tracks );
  UTCTimestamp firstTimestamp;
  UTCTimestamp lastTimestamp;
  toUTCTimestamp2( arguments->firstTimestamp, firstTimestamp );
  toUTCTimestamp2( arguments->lastTimestamp,  lastTimestamp );

  output->writeString( output,
                       "Aircraft 2.0\n%s\n%s %s\n"
                       "# Subset domain:"
                       " <min_lon> <min_lat> <max_lon> <max_lat>:\n"
                       "%lf %lf %lf %lf\n"
                       "# Dimensions: variables points tracks:\n"
                       "%d %d %lld\n"
                       "# Variable names:\n"
                       "timestamp longitude latitude elevation %s\n"
                       "# Variable units:\n"
                       "yyyymmddhhmmss deg deg m %s\n"
                       "# char notes[tracks][%d] and\n"
                       "# IEEE-754 64-bit reals "
                       "bounds[tracks][2=lon,lat][2=min,max] and\n"
                       "# MSB 64-bit integers points[tracks] and\n"
                       "# IEEE-754 64-bit reals data_1[points_1][variables] "
                       "... data_T[points_T][variables]:\n",
                       arguments->description, firstTimestamp, lastTimestamp,
                       arguments->bounds[ LONGITUDE ][ MINIMUM ],
                       arguments->bounds[ LATITUDE  ][ MINIMUM ],
                       arguments->bounds[ LONGITUDE ][ MAXIMUM ],
                       arguments->bounds[ LATITUDE  ][ MAXIMUM ],
                       VARIABLES, totalPoints, tracks->count( tracks ),
                       data->arguments.variable, data->arguments.units,
                       NOTE_LENGTH + 1 );

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

  const int trackCount = tracks->count( tracks );
  int trackIndex = 0;

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

  const int trackCount = tracks->count( tracks );
  int trackIndex = 0;

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

  const int trackCount = tracks->count( tracks );
  int trackIndex = 0;

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

  const int trackCount = tracks->count( tracks );
  int trackIndex = 0;

  do {
    const Track* const track = tracks->item( tracks, trackIndex );
    CHECK( isValidTrack( track ) );
    output->write64BitReals( output, track->data,
                             VARIABLES * track->points );
    ++trackIndex;
  } while ( AND2( output->ok( output ), trackIndex < trackCount ) );

  POST04( tracks->invariant( tracks ), tracks->count( tracks ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}




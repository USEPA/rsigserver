
/******************************************************************************
PURPOSE: Utilities.c - Utility routines.
HISTORY: 2009/10/28 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <string.h>    /* For strcmp(). */
#include <ctype.h>     /* For isspace(). */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

#include <Assertions.h>    /* For PRE*(), POST*(), CHECK*(), AND*(). */
#include <BasicNumerics.h> /* For Integer.*/
#include <Failure.h>       /* For failureSetProgramName(). */
#include <Memory.h>        /* For setCountDownToFailMemory(). */
#include <Utilities.h>     /* For public interface. */

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: isValidBounds - CHeck validity of bounds object.
INPUTS:  const Bounds bounds  Object to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidBounds( const Bounds bounds) {

  const Integer result =
    AND5( bounds,
          IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ),
          IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
                    bounds[ LONGITUDE ][ MINIMUM ], 180.0 ),
          IN_RANGE( bounds[ LATITUDE ][ MINIMUM ], -90.0, 90.0 ),
          IN_RANGE( bounds[ LATITUDE ][ MAXIMUM ],
                   bounds[ LATITUDE ][ MINIMUM ], 90.0 ) );

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND5( bounds,
                         IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ],
                                   -180.0, 180.0 ),
                         IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
                                   bounds[ LONGITUDE ][ MINIMUM ], 180.0 ),
                         IN_RANGE( bounds[ LATITUDE ][ MINIMUM ], -90.0, 90.0),
                         IN_RANGE( bounds[ LATITUDE ][ MAXIMUM ],
                                   bounds[ LATITUDE ][ MINIMUM ], 90.0 ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: pointInDomain - Is the (longitude, latitude) point inside domain?
INPUTS:  Real longitude      Longitude of point.
         Real latitude       Latitude of point.
         const Bounds domain domain[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
RETURNS: Integer 1 if inside, else 0.
******************************************************************************/

Integer pointInDomain( Real longitude, Real latitude, const Bounds domain ) {

  PRE03( isValidLongitude( longitude ),
         isValidLatitude( latitude ),
         isValidBounds( domain ) );

  const Integer result =
    AND2( IN_RANGE( longitude,
                    domain[ LONGITUDE ][ MINIMUM ],
                    domain[ LONGITUDE ][ MAXIMUM ] ),
          IN_RANGE( latitude,
                    domain[ LATITUDE ][ MINIMUM ],
                    domain[ LATITUDE ][ MAXIMUM ] ) );

  POST02( IS_BOOL( result ),
          result == AND2( IN_RANGE( longitude,
                                    domain[ LONGITUDE ][ MINIMUM ],
                                    domain[ LONGITUDE ][ MAXIMUM ] ),
                          IN_RANGE( latitude,
                                    domain[ LATITUDE ][ MINIMUM ],
                                    domain[ LATITUDE ][ MAXIMUM ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: overlap - Do rectangles overlap?
INPUTS:  const Bounds a  a[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
         const Bounds b  b[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
RETURNS: Integer 1 if overlap, else 0.
******************************************************************************/

Integer overlap( const Bounds a, const Bounds b ) {
  PRE04( a, b, isValidBounds( a ), isValidBounds( b ) );
  const Integer outside =
    OR4( a[ LATITUDE  ][ MINIMUM ] > b[ LATITUDE  ][ MAXIMUM ],
         a[ LATITUDE  ][ MAXIMUM ] < b[ LATITUDE  ][ MINIMUM ],
         a[ LONGITUDE ][ MINIMUM ] > b[ LONGITUDE ][ MAXIMUM ],
         a[ LONGITUDE ][ MAXIMUM ] < b[ LONGITUDE ][ MINIMUM ] );
  const Integer result = ! outside;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidArgs() - Are each of the strings non-zero and non-zero-length?
INPUTS:  Integer argc              Number of strings to check.
         const char* argv[ argc ]  Strings to check.
RETURNS: Integer 1 if non-zero, else 0.
NOTES:   Also set program name for failureMessage().
******************************************************************************/

Integer isValidArgs( Integer argc, const char* argv[] ) {
  Integer result = 0;

  if ( AND2( IN_RANGE( argc, 1, INT_MAX - 1 ), argv != 0 ) ) {
    Integer index = 0;

    do {

      if ( OR2( argv[ index ] == 0, argv[ index ][ 0 ] == '\0' ) ) {
        index = argc;
      }

      ++index;
    } while ( index < argc );

    result = index == argc;
    failureSetProgramName( argv[ 0 ] );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: checkForTest - Check for and set-up for testing.
INPUTS:  int* argc  Number of command-line arguments.
         char* argv[]   Command-line argument strings.
OUTPUTS: int* argc  Possibly reduced count.
NOTES:   Used to set-up simulated memory allocation failures to attain
         coverage of failure-handling code.
******************************************************************************/

void checkForTest( int* argc, char* argv[] ) {

  PRE0( isValidArgs( *argc, (const char**) argv ) );

  if ( AND3( *argc >= 3, strcmp( argv[ *argc - 2 ], "-test" ) == 0,
             atoI( argv[ *argc - 1 ] ) > 0 ) ) {
    const Integer count = atoI( argv[ *argc - 1 ] );
    *argc -= 2; /* Remove last two arguments: -test 3 */
    setCountDownToFailMemory( count ); /* Simulate failure on count NEW(). */
  }

  POST0( *argc > 0 );
}



/******************************************************************************
PURPOSE: parseArgument2 - Parse next 2 command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
         const char* option    Command string to match, e.g., "-desc".
         Integer* arg          Index of argument to parse.
OUTPUTS: Integer* arg          Index of next argument to parse.
RETURNS: const char* argument after option if successful, else 0 and
        failureMessage() is called.
******************************************************************************/

const char* parseArgument2( Integer argc, char* argv[], const char* option,
                            Integer* arg ) {

  PRE05( isValidArgs( argc, (const char**) argv ),
         option, *option, arg, *arg > 0 );
  CHECKING( const Integer OLD( arg ) = *arg; )
  const char* result = 0;

  if ( *arg + 1 >= argc ) {
    failureMessage( "Invalid/missing command-line arguments: %s.", option );
  } else if ( strcmp( argv[ *arg ], option ) ) {
    failureMessage( "Invalid command-line argument: %s (expected %s).",
                    argv[ *arg ], option );
  } else if ( ! AND3( argv[ *arg + 1 ], argv[ *arg + 1 ][ 0 ],
                      argv[ *arg + 1 ][ 0 ] != '-' ) ) {
    failureMessage( "Invalid/missing parameter to command-line argument: %s.",
                    option );
  } else {
    result = argv[ *arg + 1 ];
    *arg += 2;
  }

  POST0( IMPLIES( result, AND3( *result, *result != '-',
                                *arg == OLD( arg ) + 2 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseBounds - Parse command-line arguments for -bounds.
INPUTS:  Integer argc   Number of command-line arguments.
         char* argv[]   Command-line argument strings.
         Integer* arg   Index of argument to parse.
OUTPUTS: Integer* arg   Index of next argument to parse.
         Bounds bounds  bounds[ LONGITUDE LATITUDE ][ MINIMUM MAXIMUM ]
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

Integer parseBounds( Integer argc, char* argv[], Integer* arg, Bounds bounds) {

  PRE05( isValidArgs( argc, (const char**) argv ),
         argc > 0, arg, *arg > 0, bounds );
  CHECKING( const Integer OLD( arg ) = *arg; )

  Integer result = 0;

  if ( *arg + 4 >= argc ) {
    failureMessage( "Invalid/missing command-line arguments: -bounds." );
  } else {
    bounds[ LONGITUDE ][ MINIMUM ] = atoR( argv[ *arg + 1 ] );
    bounds[ LATITUDE  ][ MINIMUM ] = atoR( argv[ *arg + 2 ] );
    bounds[ LONGITUDE ][ MAXIMUM ] = atoR( argv[ *arg + 3 ] );
    bounds[ LATITUDE  ][ MAXIMUM ] = atoR( argv[ *arg + 4 ] );
    result = isValidBounds( (const Real(*)[2]) bounds );

    if ( result ) {
      *arg += 5;
    } else {
      failureMessage( "Invalid bounds specified [%lg %lg %lg %lg].\n",
               bounds[ LONGITUDE ][ MINIMUM ],
               bounds[ LATITUDE  ][ MINIMUM ],
               bounds[ LONGITUDE ][ MAXIMUM ],
               bounds[ LATITUDE  ][ MAXIMUM ] );
      memset( bounds, 0, sizeof (Bounds) );
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, AND2( isValidBounds( (const Real(*)[2]) bounds ),
                                 *arg == OLD( arg ) + 5 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseTimestampAndHours - Parse command-line arguments for
         -timestamp YYYYMMDDHH -hours n.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
         Integer* arg          Index of argument to parse.
OUTPUTS: Integer* arg          Index of next argument to parse.
         Integer* yyyydddhh    Timestamp.
         Integer* hours        Number of hours.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

Integer parseTimestampAndHours( Integer argc, char* argv[], Integer* arg,
                                Integer* timestamp, Integer* hours ) {

  PRE06( isValidArgs( argc, (const char**) argv ),
         argc > 0, arg, *arg > 0, timestamp, hours );
  CHECKING( Integer OLD( arg ) = *arg; )
  Integer result = 0;

  if ( *arg + 3 >= argc ) {
    failureMessage( "Invalid/missing command-line arguments: "
                    "-timestamp and -hours." );
  } else if ( strcmp( argv[ *arg ], "-timestamp" ) ) {
    failureMessage("Invalid command-line argument: %s (expecting -timestamp).",
                    argv[ *arg ] );
  } else if ( strcmp( argv[ *arg + 2 ], "-hours" ) ) {
    failureMessage("Invalid command-line argument: %s (expecting -hours).",
                    argv[ *arg ] );
  } else {
    const Integer yyyymmddhh = atoI( argv[ *arg + 1 ] );
    const Integer yyyymmdd   = yyyymmddhh / 100;
    const Integer hh         = yyyymmddhh % 100;

    if ( ! AND2( isValidYearMonthDay( yyyymmdd ), IN_RANGE( hh, 0, 23 ) ) ) {
      failureMessage( "Invalid command-line parameter for -timestamp: %s.",
                      argv[ *arg + 1 ] );
    } else {
      const Integer hh00 = hh * 100;
      const Integer yyyyddd = convertYearMonthDay( yyyymmdd );
      *timestamp = yyyyddd;
      *timestamp *= 10000; /* Make space to append HHMM. */
      *timestamp += hh00;
      *hours = atoI( argv[ *arg + 3 ] );

      if ( *hours <= 0 ) {
        failureMessage( "Invalid command-line parameter for -hours: %s.",
                        argv[ *arg + 3 ] );
      } else {
        *arg += 4;
        result = 1;
      }
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, AND3( isValidTimestamp( *timestamp ),
                                 *hours > 0,
                                 *arg = OLD( arg ) + 4 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: indexOfString - Index of string in strings[] or -1 if not present.
INPUTS:  const char* string           String to search for.
         const char* const strings[]  Strings to search.
         Integer count                    Size of strings[].
RETURNS: Integer index of string in strings[], else -1 if not present.
******************************************************************************/

Integer indexOfString( const char* string, const char* const strings[],
                       Integer count ) {

  PRE08( string, *string, strings, strings[ 0 ], *strings[ 0 ],
         count > 0, strings[ count - 1 ], *strings[ count - 1 ] );

  Integer result = -1;
  Integer index = 0;

  do {

    if ( strcmp( string, strings[ index ] ) == 0 ) {
      result = index;
      index = count;
    }

    ++index;
  } while ( index < count );

  POST0( OR2( result == -1,
             AND2( result >= 0,
                   strcmp( string, strings[ result ] ) == 0 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: lowercase - Convert a string to lowercase.
INPUTS:  char string[]  String to convert.
OUTPUTS: char string[]  Lowercase string.
******************************************************************************/

void lowercase( char string[] ) {
  PRE0( string );
  char* s = string;

  for ( s = string; *s; ++s ) {
    *s = tolower( *s );
  }
}



/******************************************************************************
PURPOSE: uppercase - Convert a string to uppercase.
INPUTS:  char string[]  String to convert.
OUTPUTS: char string[]  Uppercase string.
******************************************************************************/

void uppercase( char string[] ) {
  PRE0( string );
  char* s = string;

  for ( s = string; *s; ++s ) {
    *s = toupper( *s );
  }
}



/******************************************************************************
PURPOSE: fileExists - Determine named file exists.
INPUTS:  const char* name  Name of file to check for.
RETURNS: Integer 1 if it exists, else 0.
******************************************************************************/

Integer fileExists( const char* name ) {
  PRE0( name );
  struct stat buf;
  const Integer result = stat( name, &buf ) != -1;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: fileSize - Determine size of named file.
INPUTS:  const char* name  Name of file to examine.
RETURNS: Integer Size, in bytes, of named file, else 0 if failed.
******************************************************************************/

Integer fileSize( const char* name ) {
  PRE0( name );
  Integer result = 0;
  struct stat buf;

  if ( stat( name, &buf ) == -1 ) {
    failureMessage( "Failed to determine size of file '%s'.\n", name );
  } else {
    result = buf.st_size;

    if ( result < 0 ) {
      failureMessage( "Negative size of file '%s'.\n", name );
      result = 0;
    }
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: readFile - Read named file into memory and return it as an allocated
         string.
INPUTS:  const char* name  Name of file to examine.
OUTPUTS: Integer* length   Length of string.
RETURNS: char* string contents of file
         (with any '\r' characters converted to ' '),
         else 0 if failed and failureMessage() is called.
******************************************************************************/

char* readFile( const char* name, Integer* length ) {
  PRE02( name, length );
  char* result = 0;
  *length = fileSize( name ) / sizeof (char);

  if ( *length > 0 ) {
    result = NEW_ZERO( char, *length + 1 );

    if ( ! result ) {
      *length = 0;
    } else {
      FILE* file = fopen( name, "rb" );

      if ( file ) {
        const Integer itemsRead =
          fread( result, *length * sizeof (char), 1, file );

        if ( itemsRead != 1 ) {
          failureMessage( "Failed to read entire file '%s'.\n", name );
          FREE( result );
          *length = 0;
        } else {
          result[ *length ] = '\0'; /* Terminate string. */
          controlMToSpace( result );
        }

        fclose( file );
        file = 0;
      }
    }
  }

  POST0( IMPLIES_ELSE( result, *length > 0, *length == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: controlMToSpace - Convert any '\r' characters to ' '.
INPUTS:  char* string  String to filter.
OUTPUTS: char* string  Filtered string.
******************************************************************************/

void controlMToSpace( char* string ) {
  PRE0( string );
  char* s = 0;

  for ( s = string; *s; ++s ) {

    if ( *s == '\r' ) {
      *s = ' ';
    }
  }

  POST0( ! strchr( string, '\r' ) );
}



/******************************************************************************
PURPOSE: trimTrailingWhitespace - CHange trailing whitespace characters to '\0'
INPUTS:  size_t size          Size of string.
         char string[ size ]  String to filter.
OUTPUTS: char string[ size ]  Filtered string.
******************************************************************************/

void trimTrailingWhitespace( size_t size, char string[] ) {
  PRE02( size, string );
  char* s = string + size - 1;

  while ( AND2( s != string, isspace( *s ) ) ) {
    *s = '\0';
    --s;
  }
}



/******************************************************************************
PURPOSE: linesInString - Count number of lines in a string.
INPUTS:  const char* string  String to scan.
RETURNS: Integer number of lines in string.
******************************************************************************/

Integer linesInString( const char* string ) {
  PRE0( string );
  Integer result = 0;
  const char* c = 0;

  for ( c = string; *c; ++c ) {
    result += *c == '\n';
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: skipLines - Skip lines in a string.
INPUTS:  const char* string  String to scan.
         Integer lines       Number of lines in string to skip.
RETURNS: const char* Pointer into string after skipped lines, or 0 if at end.
******************************************************************************/

const char* skipLines( const char* string, Integer lines ) {
  PRE02( string, lines > 0 );
  const char* result = string;
  Integer count = 0;

  for ( result = string; AND2( *result, count < lines ); ++result ) {

    if ( *result == '\n' ) {
      ++count;
    }
  }

  if ( count != lines ) {
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: startOfLineWithWord - Search for word in string and if found, return
         pointer to first character in the line containing the word, else 0.
INPUTS:  const char* string   Multi-line string to search.
         const char* word*    Word to search for in string.
RETURNS: char* pointer to first character of line in string containing word,
               else 0 if not found.
******************************************************************************/

char* startOfLineWithWord( const char* string, const char* word ) {

  PRE02( string, word );

  char* result = strstr( string, word );

  if ( result ) { /* Found so back-up to beginning of line: */

    while ( AND2( result > string, *result != '\n' ) ) {
      --result;
    }

    if ( *result == '\n' ) {
      ++result;
    }
  }

  POST0( OR2( result == 0, result >= string ) );
  return result;
}



/******************************************************************************
PURPOSE: findLine - Find tag in string and skip tag, terminate end of line.
INPUTS:  char* string      String to search.
         const char* tag   String to find.
OUTPUTS: char** end        Pointer to end of line (changed to '\0').
RETURNS: char* pointer past tag if found, else 0.
******************************************************************************/

char* findLine( char* string, const char* tag, char** end ) {

  PRE03( string, tag, end );

  char* result = strstr( string, tag );
  *end = 0;

  if ( result ) {
    *end = strchr( result, '\n' );
    result += strlen( tag );

    if ( *end ) {
      **end = '\0'; /* Terminate to avoid slow strlen() in strtod(). */
    } else {
      result = 0;
    }
  }

  POST0( OR2( AND2( result != 0, **end == '\0' ), IS_ZERO2( result, *end ) ) );
  return result;
}



/******************************************************************************
PURPOSE: binIndex - Linear index into sequence of 'bins' of various counts.
INPUTS:  Integer index                0-based index.
         Integer bins                 Number of bins.
         const Integer counts[ bins ] Counts per bin.
RETURNS: Integer 0-based index into counts[] of index, else -1 if out-of-range.
******************************************************************************/

Integer binIndex( Integer index, Integer bins, const Integer counts[] ) {
  PRE04( index >= 0, bins > 0, counts, minimumItemI( counts, bins ) >= 0 );
  Integer result = -1;
  Integer bin = 0;
  Integer counter = 0;

  for ( bin = 0; bin < bins; ++bin ) {
    counter += counts[ bin ];

    if ( index < counter ) {
      result = bin;
      bin = bins;
    }
  }

  POST0( IN_RANGE( result, -1, bins - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: appendNote - Append comma-separated note to regriddedNote
         if not already present.
INPUTS:  const Note note             Note to append (space-permitting).
OUTPUTS: RegriddedNote regriddedNote Regridded note with note appended.
RETURNS: Integer 0-based index into counts[] of index, else -1 if out-of-range.
******************************************************************************/

void appendNote( RegriddedNote regriddedNote, const Note note ) {
  PRE03( regriddedNote, note, note[ 0 ] );

  if ( ! strstr( regriddedNote, note ) ) { /* If not already present: */
    Integer index = 0;

    /* Get index of end of regriddedNote: */

    while ( AND2( index < REGRIDDED_NOTE_LENGTH, regriddedNote[ index ] ) ) {
      ++index;
    }

    if ( index < REGRIDDED_NOTE_LENGTH ) { /* If there is space: */
      Integer readIndex = 0;

      if ( index > 0 ) { /* Append comma delimiter: */
        regriddedNote[ index ] = ',';
        ++index;
      }

      /* Append note in available space: */

      while ( AND2( index < REGRIDDED_NOTE_LENGTH, note[ readIndex ] ) ) {
        CHECK2( IN_RANGE( index, 0, REGRIDDED_NOTE_LENGTH ),
                IN_RANGE( readIndex, 0, NOTE_LENGTH ) );
        regriddedNote[ index ] = note[ readIndex ];
        ++index;
        ++readIndex;
      }
    }
  }

  regriddedNote[ REGRIDDED_NOTE_LENGTH ] = '\0';
  POST02( regriddedNote[ 0 ], regriddedNote[ REGRIDDED_NOTE_LENGTH ] == '\0' );
}





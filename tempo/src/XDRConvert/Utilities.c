/******************************************************************************
PURPOSE: Utilities.c - Some general-purpose reusable routines.

HISTORY: 2017-10-14 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/



/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, stderr, fprintf(). */
#include <stdlib.h>    /* For malloc(), free(), atoll(), setenv(), unsetenv()*/
#include <string.h>    /* For memset(). */
#include <ctype.h>     /* For isalnum(), isprint(). */
#include <time.h>      /* For time_t, struct_tm, mktime(), tzset(). */
#include <dirent.h>    /* For opendir(), closedir(). */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

#include "Utilities.h" /* For public interface. */

/*================================= MACROS ==================================*/

/*
 * Is the platform big-endian (MSB: most significant byte first) or
 * little-endian (LSB: least significant byte first)?
 */

#if defined(__alpha) || \
    defined(__i386__) || defined(__i486__) || \
    defined(__i586__) || defined(__i686__) || \
    defined(__ia64__) || defined(__x86_64__)
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

/*============================ GLOBAL CONSTANTS =============================*/

/*
 * 30 days hath September, April, June and November, all the rest have 31,
 * except February which has either 28 or 29 (on a leap year).
 */

static const int daysPerMonth[ 2 ][ 12 ] =
{
  { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
  { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
};

/*================================ FUNCTIONS ================================*/


static time_t my_timegm( struct tm *s ) {
  time_t seconds = 0;
  char* tz = getenv( "TZ" );
  setenv( "TZ", "UTC", 1 );
  tzset();
  seconds = mktime( s );

  if ( tz ) {
    setenv( "TZ", tz, 1 );
  } else {
    unsetenv( "TZ" );
  }

  tzset();
  return seconds;
}



/******************************************************************************
PURPOSE: isLeapYear - Is the yyyy a leap year (i.e., has 366 days)?
INPUTS:  const int yyyy year to check.
RETURNS: int 1 if leap year else 0.
******************************************************************************/

int isLeapYear( const int yyyy ) {
  return yyyy % 4 == 0 && ( yyyy % 100 != 0 || yyyy % 400 == 0 );
}



/******************************************************************************
PURPOSE: secondsSince1970 - UNIX seconds since 1970 from yyyymmddhhmmss.
INPUTS:  const pong long yyyymmddhhmmss  Timestamp.
RETURNS: int seconds since 1970.
******************************************************************************/

int secondsSince1970( const long long yyyymmddhhmmss ) {
  const int yyyy = (int) ( yyyymmddhhmmss / 10000000000LL );
  const int mm   = (int) ( yyyymmddhhmmss / 100000000LL % 100LL );
  const int dd   = (int) ( yyyymmddhhmmss / 1000000LL   % 100LL );
  const int hh   = (int) ( yyyymmddhhmmss / 10000LL     % 100LL );
  const int mm2  = (int) ( yyyymmddhhmmss / 100LL       % 100LL );
  const int ss   = (int) ( yyyymmddhhmmss               % 100LL );
  long seconds = 0;
  struct tm s;
  memset( &s, 0, sizeof s );
  s.tm_year = yyyy - 1900;
  s.tm_mon = mm - 1;
  s.tm_mday = dd;
  s.tm_hour = hh;
  s.tm_min = mm2;
  s.tm_sec = ss;
  s.tm_zone = "UTC";
  seconds = my_timegm( &s );
  return seconds;
}



/******************************************************************************
PURPOSE: fromSeconds - Compute YYYYMMDDHHMMSS from UNIX seconds since 1970.
INPUTS:  const int seconds  Seconds since 1970.
RETURNS: long long yyyymmddhhmmss.
******************************************************************************/

long long secondsToYYYYMMDDHHMMSS( const int seconds ) {
  long long yyyymmddhhmmss = 0;
  const time_t seconds2 = seconds;
  const struct tm* const s = gmtime( &seconds2 );

  if ( s ) {
    const int yyyy = s->tm_year + 1900;
    const int mm   = s->tm_mon + 1;
    const int dd   = s->tm_mday;
    const int hh   = s->tm_hour;
    const int mm2  = s->tm_min;
    const int ss   = s->tm_sec;
    yyyymmddhhmmss = yyyy;
    yyyymmddhhmmss *= 100;
    yyyymmddhhmmss += mm;
    yyyymmddhhmmss *= 100;
    yyyymmddhhmmss += dd;
    yyyymmddhhmmss *= 100;
    yyyymmddhhmmss += hh;
    yyyymmddhhmmss *= 100;
    yyyymmddhhmmss += mm2;
    yyyymmddhhmmss *= 100;
    yyyymmddhhmmss += ss;
  }

  return yyyymmddhhmmss;
}



/******************************************************************************
PURPOSE: isValidYYYYMMDDHHMMSS - Is the timestamp valid?
INPUTS:  const int yyyymmddhhmmss.
RETURNS: int 1 if valid else 0.
******************************************************************************/

int isValidYYYYMMDDHHMMSS( const long long yyyymmddhhmmss ) {
  const int yyyymmddhh = yyyymmddhhmmss / 10000;
  int result = isValidYYYYMMDDHH( yyyymmddhh );

  if ( result ) {
    const int mm = yyyymmddhhmmss / 100 % 100;
    const int ss = yyyymmddhhmmss % 100;
    result = mm >= 0 && mm <= 59 && ss >= 0 && ss <= 59;
  }

  return result;
}



/******************************************************************************
PURPOSE: isValidYYYYMMDDHH - Is the timestamp valid?
INPUTS:  const int yyyymmddhh.
RETURNS: int 1 if valid else 0.
******************************************************************************/

int isValidYYYYMMDDHH( const int yyyymmddhh ) {
  const int yyyy = yyyymmddhh / 1000000;
  const int mm   = yyyymmddhh / 10000 % 100;
  const int dd   = yyyymmddhh / 100 % 100;
  const int hh   = yyyymmddhh % 100;
  const int leap = isLeapYear( yyyy );
  const int result =
    IN_RANGE( yyyy, 1900, 9999 ) &&
    IN_RANGE( mm, 1, 12 ) &&
    IN_RANGE( dd, 1, daysPerMonth[ leap ][ mm - 1 ] ) &&
    IN_RANGE( hh, 0, 23 );
  return result;
}



/******************************************************************************
PURPOSE: isValidBounds - CHeck validity of bounds object.
INPUTS:  const Bounds bounds  Object to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

int isValidBounds( const Bounds bounds) {
  const int result =
    bounds != 0 &&
    IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ) &&
    IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
              bounds[ LONGITUDE ][ MINIMUM ], 180.0 ) &&
    IN_RANGE( bounds[ LATITUDE  ][ MINIMUM ], -90.0, 90.0 ) &&
    IN_RANGE( bounds[ LATITUDE  ][ MAXIMUM ],
              bounds[ LATITUDE  ][ MINIMUM ], 90.0 );
  return result;
}



/******************************************************************************
PURPOSE: rotate8ByteArrayIfLittleEndian() - Rotate 8-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array         Array of 8-byte values to rotate.
         const size_t count  Number of items in array.
OUTPUTS: void* array         Array of rotated values.
******************************************************************************/

void rotate8ByteArrayIfLittleEndian( void* array, const size_t count) {

#if IS_LITTLE_ENDIAN

  long long* const array8 = array;
  long long index = 0;
  assert( array ); assert( count > 0 );
  assert( sizeof (long long) == 8 );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const long long value = array8[ index ];
    const long long newValue =
    ( value & 0xff00000000000000LL ) >> 56 |
    ( value & 0x00ff000000000000LL ) >> 40 |
    ( value & 0x0000ff0000000000LL ) >> 24 |
    ( value & 0x000000ff00000000LL ) >>  8 |
    ( value & 0x00000000ff000000LL ) <<  8 |
    ( value & 0x0000000000ff0000LL ) << 24 |
    ( value & 0x000000000000ff00LL ) << 40 |
    ( value & 0x00000000000000ffLL ) << 56;
    array8[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: isDirectory - Determine if name is a directory.
INPUTS:  const char* name  Name of directory to examine.
RETURNS: int 1 if directory, else 0.
******************************************************************************/

int isDirectory( const char* name ) {
  int result = 0;
  DIR* directory = 0;
  assert( name );
  directory = opendir( name );
  result = directory != 0;

  if ( directory ) {
    closedir( directory );
    directory = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: fileSize - Determine size of named file.
INPUTS:  const char* name  Name of file to examine.
RETURNS: size_t size, in bytes, of named file, else 0 if failed.
******************************************************************************/

size_t fileSize( const char* name ) {
  size_t result = 0;
  struct stat buf;
  assert( name );

  if ( stat( name, &buf ) == -1 ) {
    fprintf( stderr, "\nFailed to determine size of file '%s'.\n", name );
  } else if ( buf.st_size > 0 ) {
    result = buf.st_size;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFile - Read named file into memory.
INPUTS:  const char* name  Name of file to examine.
         size_t* length   Current length of content.
         char** content   Currently allocated buffer of given length or 0.
OUTPUTS: size_t* length   Length of content if increased.
         char** content   Contents of file.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   All '\r' characters are converted to ' '.
******************************************************************************/

int readFile( const char* name, size_t* length, char** content ) {
  int result = 0;
  size_t fileLength = 0;
  assert( name ); assert( *name );
  assert( length ); assert( content );
  assert( *length == 0 || *content );
  fileLength = fileSize( name ) / sizeof (char);

  if ( fileLength > 0 ) {

    if ( fileLength > *length ) {
      const size_t bytes = ( fileLength + 1 ) * sizeof (char);

      if ( *content ) {
        free( *content ), *content = 0;
        *length = 0;
      }

      *content = malloc( bytes );

      if ( ! *content ) {
        fprintf( stderr,
              "\nCan't allocate %lu bytes to complete the requested action.\n",
                 bytes );
      } else {
        *length = fileLength;
      }
    }

    if ( *content ) {
      FILE* file = fopen( name, "rb" );

      if ( file ) {
        const size_t itemsRead =
          fread( *content, sizeof (char), fileLength, file );

        if ( itemsRead != fileLength ) {
          fprintf( stderr, "\nFailed to read entire file '%s'.\n", name );
          free( *content ), *content = 0;
          *length = 0;
        } else {
          (*content)[ fileLength ] = '\0'; /* Terminate string. */
          result = 1;
        }

        fclose( file );
        file = 0;
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: eraseQuotedCommasAndQuotesAndBrackets - Convert any ',' characters to
         ' ' when they appear inside double-quotes (") and also replace
         double-quotes and '[]' characters with ' '.
INPUTS:  char* string  String to filter.
OUTPUTS: char* string  Filtered string.
NOTES:   This routine is a profiled hotspot so it needs to be highly optimized.
******************************************************************************/

void eraseQuotedCommasAndQuotesAndBrackets( char* string ) {
  unsigned char quoted = 0; /* Inside double-quoted string? */
  assert( string );
loop:
  switch ( *string ) {
  case '"':
    quoted ^= 1; /* Toggle flag. */
    /* lint -fallthrough */
  case '[':
  case ']':
    *string = ' ';
    break;
  case ',':
    if ( quoted ) {
      *string = ' ';
    }
    break;
  case '\0':
    return;
  default:
    break;
  }
  ++string;
  goto loop;
}



/******************************************************************************
PURPOSE: nextLine - Extract a line of text from a string.
INPUTS:  char** line  String to parse.
OUTPUTS: char** line  0-terminated line with no leading whitespace.
         char** end   End of line within string.
RETURNS: int 1 if extracted line is non-empty.
******************************************************************************/

int nextLine( char** line, char** end ) {
  int result = 0;
  char* s = 0;
  char* e = 0;
  assert( line ); assert( end );
  s = *line;

  /* Trim leading whitespace from line: */

  while ( isspace( *s ) ) {
    ++s;
  }

  *line = s;
  *end = s;

  if ( *s ) {
    result = 1; /* Non-empty line remains. */

    /* Terminate line: */

    e = strchr( s, '\n' );

    if ( e ) {
      *e = '\0';
      *end = e;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: indexOfWord - Find 0-based index of word in string of single
         space-delimited words.
INPUTS:  const char* word   String to find.
         const char* words  String to search.
RETURNS: int 0-based index of word in words, else -1.
******************************************************************************/

int indexOfWord( const char* word, const char* words ) {
  int result = -1;
  const char* s = 0;
  size_t wordLength = 0;

  assert( word );  assert( isalnum( *word  ) ); assert( ! strchr( word, ' ' ));
  assert( words ); assert( isalnum( *words ) ); assert(   strchr( words, ' '));
  assert( ! strstr( words, "  " ) ); /* No consecutive spaces. */

  wordLength = strlen( word );
  s = strstr( words, word );

  while ( s ) {
    const char c = s[ wordLength ];
    const int found =
      ( c == ' ' || c == '\0' ) && ( s == words || *( s - 1 ) == ' ' );

    if ( found ) { /* Count words (spaces) before. */
      result = 0;

      while ( s > words ) {

        if ( *s == ' ' ) {
          ++result;
        }

        --s;
      }

      s = 0; /* Stop looping. */
    } else {
      s = strstr( s + wordLength, word );
    }
  }

  assert( result >= -1 );
  assert( result == -1 || strstr( words, word ) );
  return result;
}



/******************************************************************************
PURPOSE: parseOptions - Parse command-line options.
INPUTS:  const int argc           Number of command-line arguments.
         char argv[]              Command-line argument strings.
         const int count          Number of options.
         Option options[ count ]  Option: name, required, type, count, range.
OUTPUTS: Option options[ count ]  Option: parsed, values.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int parseOptions( const int argc, char* argv[], const int count,
                  Option options[] ) {
  int result = 1;
  int arg = 1;
  int option = 0;
  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] ); assert( count > 0 ); assert( options );

  for ( option = 0; option < count; ++option ) {
    options[ option ].parsed = 0;
  }

  while ( result && arg < argc ) {
    const char* const argument = argv[ arg ];

    for ( option = 0; result && option < count; ++option ) {
      Option* const thisOption = options + option;
      const char* name = thisOption->name;

      if ( ! strcmp( argument, name ) ) {
        result = parseOption( argc, argv, &arg, thisOption );
        option = count; /* Found so stop looping. */
      }
    }

    if ( option == count ) { /* Not found. */
      fprintf( stderr, "%s: Invalid command-line option %s.\n",
               argv[ 0 ], argument );
      result = 0;
    }
  } /* End loop on argv[]. */

  /* Check that required arguments were provided: */

  for ( option = 0; result && option < count; ++option ) {
    Option* const  thisOption = options + option;
    result = ! thisOption->required || thisOption->parsed;

    if ( ! result ) {
      fprintf( stderr, "%s: Missing required command-line option %s.\n",
               argv[ 0 ], thisOption->name );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: parseOption - Parse a command-line option.
INPUTS:  const int argc        Number of command-line arguments.
         char argv[]           Command-line argument strings.
         int* const arg        Number of argument being parsed.
         Option* const option  Option: name, required, type, count, range.
OUTPUTS: int* const arg        Increased number of argument being parsed.
         Option* const option  Option: values.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int parseOption( const int argc, char* argv[],
                 int* const arg, Option* const option ) {
  int result = 0;
  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] );
  assert( arg ), assert( *arg > 0 ); assert( *arg < argc );
  assert( option ); assert( option->count >= 0 );

  {
    const char* argument = argv[ *arg ];
    const int valueCount = option->count;

    if ( option->parsed ) {
      fprintf( stderr, "%s: Redundant command-line option %s.\n",
               argv[ 0 ], argument );
    } else if ( *arg + valueCount >= argc ) {
      fprintf( stderr,
               "%s: Require %d values for command-line option %s.\n",
               argv[ 0 ], valueCount, argument );
    } else {
      void* const values = option->values;
      option->parsed = 1;

      ++*arg; /* Skip argument. */

      if ( values && valueCount > 0 ) {
        const int type = option->type;
        int valueIndex = 0;

        for ( valueIndex = 0, result = 1;
              result && *arg < argc && valueIndex < valueCount;
              ++valueIndex, ++*arg ) {
          argument = argv[ *arg ];
          result = parseOptionValue( argument, valueIndex, option );
        }

        if ( result ) { /* Perform additional type-specific checks: */

          if ( type == YYYYMMDDHHMMSS_TYPE && valueCount == 2 ) {
            long long* const ivalues = values;
            const long long yyyymmddhhmmss1 = ivalues[ 0 ];
            const long long yyyymmddhhmmss2 = ivalues[ 1 ];
            result = yyyymmddhhmmss1 <= yyyymmddhhmmss2;
          } else if ( type == BOUNDS_TYPE ) {
            result = 0;

            if ( valueCount >= 4 ) {
              double* rvalues = values;
              const double latitudeMinimum  = rvalues[ 1 ];
              const double longitudeMaximum = rvalues[ 2 ];
              rvalues[ 1 ] = longitudeMaximum; /* Reorder to 2d array. */
              rvalues[ 2 ] = latitudeMinimum;
              result = isValidBounds( (const double (*)[2]) rvalues );
            }
          }
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: parseOptionValue - Parse a command-line option value.
INPUTS:  const char* const argument  Command-line argument string.
         const int valueIndex        0-based index of option value being parsed.
         Option* const option  Option: type, values, range.
OUTPUTS: Option* const option  Option: values[ valueIndex ].
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int parseOptionValue( const char* const argument, const int valueIndex,
                      Option* const option ) {
  int result = 0;
  assert( argument ); assert( argument[ 0 ] );
  assert( valueIndex >= 0 );
  assert( option ); assert( option->values ); assert( option->type >= 0 );

  {
    void* const values = option->values;

    switch ( option->type ) {
    case FILE_TYPE:
      {
        const char** svalues = values;
        const char** svalue = svalues + valueIndex;
        *svalue = argument;
        result = fileSize( *svalue ) > 0;
      }
      break;
    case DIRECTORY_TYPE:
      {
        const char** svalues = values;
        const char** svalue = svalues + valueIndex;
        *svalue = argument;
        result = isDirectory( *svalue );
      }
      break;
    case STRING_TYPE:
      {
        const char** svalues = values;
        const char** svalue = svalues + valueIndex;
        *svalue = argument;
        result = isprint( *svalue[ 0 ] );
      }
      break;
    case ENUM_TYPE:
      {
        int* const ivalues = values;
        int* const ivalue = ivalues + valueIndex;
        result = isalnum( *argument ) && ! strchr( argument, ' ' );

        if ( result ) {
          *ivalue = indexOfWord( argument, option->valids );
          result = *ivalue >= 0;
        }
      }
      break;
    case INT_TYPE:
      {
        int* const ivalues = values;
        int* const ivalue = ivalues + valueIndex;
        const int* const range = option->range;
        char* end = 0;
        *ivalue = strtol( argument, &end, 10 );
        result = end != argument;

        if ( range ) {
          const int minimum = range[ 0 ];
          const int maximum = range[ 1 ];
          result = *ivalue >= minimum && *ivalue <= maximum;
        }
      }
      break;
    case INTEGER64_TYPE:
      {
        long long* const ivalues = values;
        long long* const ivalue = ivalues + valueIndex;
        const long long* const range = option->range;
        char* end = 0;
        *ivalue = strtoll( argument, &end, 10 );
        result = end != argument;

        if ( range ) {
          const long long minimum = range[ 0 ];
          const long long maximum = range[ 1 ];
          result = *ivalue >= minimum && *ivalue <= maximum;
        }
      }
      break;
    case REAL64_TYPE:
      {
        double* const rvalues = values;
        double* const rvalue = rvalues + valueIndex;
        const double* const range = option->range;
        char* end = 0;
        *rvalue = strtod( argument, &end );
        result = end != argument;

        if ( range ) {
          const double minimum = range[ 0 ];
          const double maximum = range[ 1 ];
          result = *rvalue >= minimum && *rvalue <= maximum;
        }
      }
      break;
    case YYYYMMDDHHMMSS_TYPE:
      {
        long long* const ivalues = values;
        long long* const ivalue = ivalues + valueIndex;
        char* end = 0;
        *ivalue = strtoll( argument, &end, 10 );
        result = end != argument;

        if ( result ) {
          result = isValidYYYYMMDDHHMMSS( *ivalue );
        }
      }
      break;
    case BOUNDS_TYPE:
      {
        double* const rvalues = values;
        double* const rvalue = rvalues + valueIndex;
        char* end = 0;
        *rvalue = strtod( argument, &end );
        result = end != argument;
      }
      break;
    default:
      break;
    }
  }

  return result;
}




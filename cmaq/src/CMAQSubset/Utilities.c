/******************************************************************************
PURPOSE: Utilities.c - Some general-purpose reusable routines.

HISTORY: 2017-10-14 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/



/*================================ INCLUDES =================================*/

#include <stdio.h>     /* For FILE, stderr, fprintf(). */
#include <stdlib.h>    /* For malloc(), free(), atoll(), qsort(), setenv()*/
#include <string.h>    /* For memset(). */
#include <ctype.h>     /* For isalnum(), isprint(). */
#include <float.h>     /* For DBL_MAX. */
#include <time.h>      /* For time_t, struct_tm, gmtime_r(). */
#include <unistd.h>    /* For getcwd(). */
#include <dirent.h>    /* For opendir(), closedir(). */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */
#ifndef _WIN32
#include <sys/utsname.h> /* For uname(). */
#endif

#include "Assertions.h" /* For PRE*(), POST*(), DEBUG(). */
#include "Utilities.h"  /* For public interface. */

/*================================= MACROS ==================================*/

/*
 * Is the platform big-endian (MSB: most significant byte first) or
 * little-endian (LSB: least significant byte first)?
 */

#ifndef IS_LITTLE_ENDIAN

#if \
  ( \
    ( defined(__BYTE_ORDER__) && \
      defined(__ORDER_LITTLE_ENDIAN__) && \
      __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ ) || \
    defined(__x86_64__) || \
    defined(__i386__)   || \
    defined(__i486__)   || \
    defined(__i586__)   || \
    defined(__i686__)   || \
    defined(__alpha)    || \
    defined(__ia64__)   || \
    defined(__ARMEL__)  || \
    defined(__MIPSEL)   || \
    defined(_WIN32)     || \
    defined(_WIN64)     || \
    defined(_M_IX86)    || \
    defined(_M_AMD64)      \
  )
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif

#endif

#ifdef _WIN32
#ifdef _WIN64
/*
 * HACK non-standard gmtime_s() exists but with a different argument order
 * and return logic, but we can use it with adjustments:
 */
#define gmtime_r( clock, result ) ( gmtime_s( result, clock ) ? 0 : result )
#else
/*
 * gmtime_s() is missing in mingw32 so just use gmtime().
 * According to this reference there is no actual race-condition:
 * https://sourceforge.net/p/mingw-w64/mailman/message/21523522/
 * "Please note that plain old gmtime() and localtime() *are* thread-safe
 * in Microsoft's C library. They use a thread-local buffer."
 */
static struct tm* gmtime_r( const time_t* clock, struct tm* result ) {
  const struct tm* const result0 = gmtime( clock ); /* HACK: not thread-safe!*/
  const int ok = result0 != 0;
  if ( result0 ) memcpy( result, result0, sizeof (struct tm) );
  return ok ? result : 0;
}
#endif
#endif


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



/******************************************************************************
PURPOSE: newMemory - Implements NEW() macro by calling the standard malloc()
         and, if successful, zeros the memory else upon failure prints a
         message to stderr.
INPUTS:  size_t count     The number of items to allocate.
         size_t sizeEach  The number of bytes per item.
RETURNS: void*  The resulting address, zeroed-out if successful,
                              else returns 0 if unsuccessful.
NOTES:   If malloc() returns 0 then a failure message is printer to stderr.
******************************************************************************/

void* newMemory( size_t count, size_t sizeEach ) {
  PRE02( count > 0, sizeEach > 0 );
  void* address = 0;
  const size_t bytes = count * sizeEach;

  /* Don't attempt to allocate if too large to represent: */

  if ( AND3( bytes > 0, bytes >= count, bytes >= sizeEach ) ) {
    address = malloc( bytes );
  }

  if ( address == 0 ) {
    fprintf ( stderr,
             "\nCan't allocate %lu bytes to complete the requested action.\n",
             bytes );
  } else {
    memset( address, 0, bytes );
  }

  return address;
}



/******************************************************************************
PURPOSE: inRange - Are values within range [minimum, maximum]?
INPUTS:  const size_ count             Number of values.
         const double values[ count ]  Values to check.
         const double minimum          Minimum value to check for.
         const double maximum          Maximum value to check for.
RETURNS: int 1 if all values are within range, else 0.
******************************************************************************/

int inRange( const size_t count, const double values[],
             const double minimum, const double maximum ) {
  PRE02( count, values );
  int result = 1;
  size_t counter = count;

  while ( AND2( result, counter-- ) ) {
    const double value = *values++;
    result = IN_RANGE( value, minimum, maximum );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: inRangeF - Are values within range [minimum, maximum]?
INPUTS:  const size_ count             Number of values.
         const float values[ count ]  Values to check.
         const double minimum          Minimum value to check for.
         const double maximum          Maximum value to check for.
RETURNS: int 1 if all values are within range, else 0.
******************************************************************************/

int inRangeF( const size_t count, const float values[],
             const double minimum, const double maximum ) {
  PRE02( count, values );
  int result = 1;
  size_t counter = count;

  while ( AND2( result, counter-- ) ) {
    const double value = *values++;
    result = IN_RANGE( value, minimum, maximum );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: valuesIncrease - Are values[ i ] < values[ i + 1 ]?
INPUTS:  const size_ count             Number of values.
         const double values[ count ]  Values to check.
RETURNS: int 1 if increasing, else 0.
******************************************************************************/

int valuesIncrease( const size_t count, const double values[] ) {
  PRE02( count, values );
  int result = 1;
  size_t index = 1;

  for ( ; AND2( result, index < count ); ++index ) {
    result = values[ index ] > values[ index - 1 ];
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: valuesDecrease - Are values[ i ] > values[ i + 1 ]?
INPUTS:  const size_ count             Number of values.
         const float values[ count ]  Values to check.
RETURNS: int 1 if increasing, else 0.
******************************************************************************/

int valuesDecrease( const size_t count, const float values[] ) {
  PRE02( count, values );
  int result = 1;
  size_t index = 1;

  for ( ; AND2( result, index < count ); ++index ) {
    result = values[ index ] < values[ index - 1 ];
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: underscoreToSpace - Convert any '_' characters to ' '.
INPUTS:  char* string  String to convert.
OUTPUTS: char* string  Converted string.
******************************************************************************/

void underscoreToSpace( char* string ) {
  PRE0( string );

  while ( *string ) {

    if ( *string == '_' ) {
      *string = ' ';
    }

    ++string;
  }
}



/******************************************************************************
PURPOSE: indexOfWord - Find 0-based index of word in string of single
         space-delimited words.
INPUTS:  const char* word   String to find.
         const char* words  String to search.
RETURNS: int 0-based index of word in words, else -1.
******************************************************************************/

int indexOfWord( const char* word, const char* words ) {

  PRE07( word, isalnum( *word ), ! strchr( word, ' ' ),
         words, isalnum( *words ), strchr( words, ' '),
         ! strstr( words, "  ") );

  int result = -1;
  const char* s = 0;
  size_t wordLength = 0;
  wordLength = strlen( word );
  s = strstr( words, word );

  while ( s ) {
    const char c = s[ wordLength ];
    const int found =
      AND2( IN3( c, ' ', '\0' ), OR2( s == words, *( s - 1 ) == ' ' ) );

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

  POST02( result >= -1, OR2( result == -1, strstr( words, word ) ) );
  return result;
}



/******************************************************************************
PURPOSE: paddedString - Append spaces to a string up to a specified length.
INPUTS:  const char* const input    String to copy.
         const size_t length        Length of padded output.
OUTPUTS: char output[ length + 1 ]  Padded, null-terminated string.
RETURNS: char* output.
******************************************************************************/

char* paddedString( const char* const input, const size_t length,
                    char output[] ) {
  PRE04( input, length, output, input != output );
  CHECKING( const size_t inputLength = strlen( input ); )
  const char* source = input;
  char* result = output;
  size_t index = 0;

  while (  AND2( *source, index < length ) ) {
    *result++ = *source++;
    ++index;
  }

  while ( index < length ) {
    *result++ = ' ';
    ++index;
  }

  *result = '\0';
  result = output;
  CHECK( strlen( input ) == inputLength );
  POST04( result == output,
          output[ length ] == '\0',
          strlen( output ) == length,
          IMPLIES_ELSE( length > strlen( input ),
                        output[ length - 1 ] == ' ',
                        output[ length - 1 ] == input[ index - 1 ] ) );
  return result;
}



/******************************************************************************
PURPOSE: parseOptions - Parse command-line options.
INPUTS:  const int argc           Number of command-line arguments.
         char argv[]              Command-line argument strings.
         const int count          Number of options.
         Option options[ count ]  Option: name, required, type, count, range.
OUTPUTS: Option options[ count ]  Option: parsed, values.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int parseOptions( const int argc, char* argv[], const int count,
                  Option options[] ) {

  PRE06( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], count > 0, options );

  int result = 1;
  int arg = 1;
  int option = 0;

  for ( option = 0; option < count; ++option ) {
    options[ option ].parsed = 0;
  }

  while ( AND2( result, arg < argc ) ) {
    const char* const argument = argv[ arg ];

    for ( option = 0; AND2( result, option < count ); ++option ) {
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

  for ( option = 0; AND2( result, option < count ); ++option ) {
    Option* const thisOption = options + option;
    result = IMPLIES( thisOption->required, thisOption->parsed );

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
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int parseOption( const int argc, char* argv[],
                 int* const arg, Option* const option ) {

  PRE08( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], arg , *arg > 0,
         *arg < argc, option );

  int result = 0;

  {
    const char* argument = argv[ *arg ];
    const int upToCount = option->count < 0;
    const int valueCount = option->count < 0 ? -option->count : option->count;

    if ( option->parsed ) {
      fprintf( stderr, "%s: Redundant command-line option %s.\n",
               argv[ 0 ], argument );
    } else if ( AND3( upToCount, option->required, *arg + 1 >= argc ) ) {
      fprintf( stderr,
               "%s: Require at least 1 value for command-line option %s.\n",
               argv[ 0 ], argument );
    } else if ( AND2( ! upToCount, *arg + valueCount >= argc ) ) {
      fprintf( stderr,
               "%s: Require %d values for command-line option %s.\n",
               argv[ 0 ], valueCount, argument );
    } else {
      void* const values = option->values;
      option->parsed = 1;
      ++*arg; /* Skip argument. */
      result = valueCount == 0;

      if ( AND3( values, valueCount > 0, *arg < argc ) ) {
        int valueIndex = 0;
        int isViableValue = 0;

        do {
          const char* const optionValue = argv[ *arg ];
          isViableValue =
            OR2( optionValue[ 0 ] != '-', isDouble( optionValue ) );
          result = OR2( upToCount, isViableValue );

          if ( ! result ) {
            fprintf( stderr,
               "%s: Require value for command-line option %s.\n",
               argv[ 0 ], option->name );
          } else if ( isViableValue ) {
            result = parseOptionValue( optionValue, valueIndex, option );

            if ( ! result ) {
              fprintf( stderr,
                 "%s: Invalid value '%s' for command-line option %s.\n",
                 argv[ 0 ], optionValue , option->name );
            } else {
              option->parsed += 1;
              ++*arg;
              ++valueIndex;
            }
          }

        } while ( AND4( result, *arg < argc, isViableValue,
                        valueIndex < valueCount ) );

        if ( result ) {
          result =
             OR2( AND2( upToCount == 0, valueIndex == valueCount ),
                  AND2( upToCount == 1, valueIndex <= valueCount ) );

          if ( ! result ) {
            fprintf( stderr,
                     "%s: Invalid value count for command-line option %s.\n",
                     argv[ 0 ], option->name );
          } else { /* Perform additional type-specific checks: */
            const int type = option->type;
            const int parsedValueCount = valueIndex;

            if ( AND2( type == YYYYMMDDHH_TYPE, parsedValueCount == 2 ) ) {
              const int* const ivalues = values;
              const int yyyymmddhh1 = ivalues[ 0 ];
              const int yyyymmddhh2 = ivalues[ 1 ];
              result = yyyymmddhh1 <= yyyymmddhh2;

              if ( ! result ) {
                fprintf( stderr,
                         "%s: Require %d <= %d for command-line option %s.\n",
                         argv[ 0 ], yyyymmddhh1, yyyymmddhh2, option->name );
              }
            } else if ( type == BOUNDS_TYPE ) {
              result = 0;

              if ( parsedValueCount >= 4 ) {
                double* rvalues = values;
                const double latitudeMinimum  = rvalues[ 1 ];
                const double longitudeMaximum = rvalues[ 2 ];
                rvalues[ 1 ] = longitudeMaximum; /* Reorder to 2d array. */
                rvalues[ 2 ] = latitudeMinimum;
                result = isValidBounds( (const double (*)[2]) rvalues );

                if ( ! result ) {
                  fprintf( stderr,
                           "%s: Invalid bounds for command-line option %s.\n",
                           argv[ 0 ], option->name );
                }
              }
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
INPUTS:  const char* const optionValue  Command-line option value string.
         const int valueIndex        0-based index of option value being parsed.
         Option* const option  Option: type, values, range.
OUTPUTS: Option* const option  Option: values[ valueIndex ].
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int parseOptionValue( const char* const optionValue, const int valueIndex,
                      Option* const option ) {

  PRE06( optionValue, optionValue[ 0 ], valueIndex >= 0, option, option->values,
        option->type >= 0 );

  int result = 0;

  {
    void* const values = option->values;

    switch ( option->type ) {
    case FILE_TYPE:
      {
        const char** svalues = values;
        const char** svalue = svalues + valueIndex;
        *svalue = optionValue;
        result = fileSize( *svalue ) > 0;
      }
      break;
    case DIRECTORY_TYPE:
      {
        const char** svalues = values;
        const char** svalue = svalues + valueIndex;
        *svalue = optionValue;
        result = isDirectory( *svalue );
      }
      break;
    case STRING_TYPE:
      {
        const char** svalues = values;
        const char** svalue = svalues + valueIndex;
        *svalue = optionValue;
        result = isprint( *svalue[ 0 ] );
      }
      break;
    case ENUM_TYPE:
      {
        int* const ivalues = values;
        int* const ivalue = ivalues + valueIndex;
        result = AND2( isalnum( *optionValue ), ! strchr( optionValue, ' ' ) );

        if ( result ) {
          *ivalue = indexOfWord( optionValue, option->valids );
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
        *ivalue = strtol( optionValue, &end, 10 );
        result = end != optionValue;

        if ( range ) {
          const int minimum = range[ 0 ];
          const int maximum = range[ 1 ];
          result = IN_RANGE( *ivalue, minimum, maximum );
        }
      }
      break;
    case INTEGER64_TYPE:
      {
        long long* const ivalues = values;
        long long* const ivalue = ivalues + valueIndex;
        const long long* const range = option->range;
        char* end = 0;
        *ivalue = strtoll( optionValue, &end, 10 );
        result = end != optionValue;

        if ( range ) {
          const long long minimum = range[ 0 ];
          const long long maximum = range[ 1 ];
          result = IN_RANGE( *ivalue, minimum, maximum );
        }
      }
      break;
    case REAL64_TYPE:
      {
        double* const rvalues = values;
        double* const rvalue = rvalues + valueIndex;
        const double* const range = option->range;
        char* end = 0;
        *rvalue = strtod( optionValue, &end );
        result = end != optionValue;

        if ( range ) {
          const double minimum = range[ 0 ];
          const double maximum = range[ 1 ];
          result = IN_RANGE( *rvalue, minimum, maximum );
        }
      }
      break;
    case YYYYMMDDHH_TYPE:
      {
        int* const ivalues = values;
        int* const ivalue = ivalues + valueIndex;
        char* end = 0;
        *ivalue = strtol( optionValue, &end, 10 );
        result = end != optionValue;

        if ( result ) {
          result = isValidYYYYMMDDHH( *ivalue );
        }
      }
      break;
    case BOUNDS_TYPE:
      {
        double* const rvalues = values;
        double* const rvalue = rvalues + valueIndex;
        char* end = 0;
        *rvalue = strtod( optionValue, &end );
        result = end != optionValue;
      }
      break;
    default:
      break;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: isDouble - Is the string parsable as a finite double?
INPUTS:  const char* const string  String to parse.
RETURNS: int 1 if string is parsable as a finite double else 0.
******************************************************************************/

int isDouble( const char* const string ) {
  char* end = 0;
  const double value = strtod( string, &end );
  const int result =
    AND2( end != string, IN_RANGE( value, -DBL_MAX, DBL_MAX ) );
  POST0( IS_BOOL( result ) );
  return result;
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
PURPOSE: isValidYYYYDDD - Is the timestamp valid?
INPUTS:  const int yyyyddd.
RETURNS: int 1 if valid else 0.
******************************************************************************/

int isValidYYYYDDD( const int yyyyddd ) {
  const int yyyy = yyyyddd / 1000;
  const int ddd  = yyyyddd % 1000;
  const int leap = isLeapYear( yyyy );
  const int result =
    AND2( IN_RANGE( yyyy, 1900, 9999 ), IN_RANGE( ddd, 1, 365 + leap ) );
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
    AND4( IN_RANGE( yyyy, 1900, 9999 ),
          IN_RANGE( mm, 1, 12 ),
          IN_RANGE( dd, 1, daysPerMonth[ leap ][ mm - 1 ] ),
          IN_RANGE( hh, 0, 23 ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidHHMMSS - Is the timestamp valid?
INPUTS:  const int hhmmss.
RETURNS: int 1 if valid else 0.
******************************************************************************/

int isValidHHMMSS( const int hhmmss ) {
  const int hh = hhmmss / 10000;
  const int mm = hhmmss / 100 % 100;
  const int ss = hhmmss % 100;
  const int result =
    AND3( IN_RANGE( hh, 0, 23 ), IN_RANGE( mm, 0, 59 ), IN_RANGE( ss, 0, 59 ));
  return result;
}



/******************************************************************************
PURPOSE: toYYYYMMDD - Convert yyyyddd to yyyymmdd.
INPUTS:  const int yyyyddd.
RETURNS: int yyyymmdd.
******************************************************************************/

int toYYYYMMDD( const int yyyyddd ) {
  PRE0( isValidYYYYDDD( yyyyddd ) );
  const int yyyy = yyyyddd / 1000;
  int ddd        = yyyyddd % 1000;
  const int leap = isLeapYear( yyyy );
  int mm = 0;
  int done = 0;
  int result = 0;
  CHECK2( IS_BOOL( leap ), IN_RANGE( ddd, 1, 366 ) );

  do {
    const int monthDays = daysPerMonth[ leap ][ mm ];
    CHECK( IN_RANGE( monthDays, 28, 31 ) );

    if ( ddd > monthDays ) {
      ddd -= monthDays;
    } else {
      done = 1;
    }

    ++mm;
    CHECK( OR2( done, mm < 12 ) );
  } while ( done == 0 );

  CHECK2( IN_RANGE( mm, 1, 12 ),
          IN_RANGE( ddd, 1, daysPerMonth[ leap ][ mm - 1 ] ) );
  result = ( yyyy * 100 + mm ) * 100 + ddd;
  POST0( isValidYYYYMMDDHH( result * 100 ) );
  return result;
}



/******************************************************************************
PURPOSE: toYYYYDDD - Convert yyyymmdd to yyyyddd.
INPUTS:  const int yyyymmdd.
RETURNS: int yyyyddd.
******************************************************************************/

int toYYYYDDD( const int yyyymmdd ) {
  PRE0( isValidYYYYMMDDHH( yyyymmdd * 100 ) );
  const int yyyy = yyyymmdd / 10000;
  const int mm0  = yyyymmdd / 100 % 100 - 1;
  const int dd   = yyyymmdd % 100;
  const int leap = isLeapYear( yyyy );
  int result = 0;
  int month = 0;

  for ( month = 0; month < mm0; ++month ) {
    CHECK( IN_RANGE( month, 0, 11 ) );
    result += daysPerMonth[ leap ][ month ];
  }

  result += dd;
  result += yyyy * 1000;
  POST0( isValidYYYYDDD( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nowUTC - Current timestamp in UTC.
OUTPUTS: int* const yyyy  Year.
         int* const ddd   Day of year.
         int* const hh    Hour of day.
         int* const mm    Minute of hour.
         int* const ss    Second of minute.
******************************************************************************/

void nowUTC( int* const yyyy,
             int* const ddd,
             int* const hh,
             int* const mm,
             int* const ss ) {
  PRE05( yyyy, ddd, hh, mm, ss );
  time_t time0 = time( 0 );
  struct tm timeInfo;

  if ( gmtime_r( &time0, &timeInfo ) ) {
    *yyyy = timeInfo.tm_year + 1900;
    *ddd  = timeInfo.tm_yday + 1;
    *hh   = timeInfo.tm_hour;
    *mm   = timeInfo.tm_min;
    *ss   = timeInfo.tm_sec;
  } else {
    *yyyy = 1970;
    *ddd  = 1;
    *hh = 0;
    *mm = 0;
    *ss = 0;
  }

  POST02( isValidYYYYDDD( *yyyy * 1000 + *ddd ),
          isValidHHMMSS( ( *hh * 100 + *mm ) * 100 + *ss ) );
}



/******************************************************************************
PURPOSE: daysInMonth - Number of days in given month.
INPUTS:  const int yyyy  Year.
         const int mm    Month.
RETURNS: int number of days in yyyymm.
******************************************************************************/

int daysInMonth( const int yyyy, const int mm ) {
  PRE02( IN_RANGE( yyyy, 1900, 9999 ), IN_RANGE( mm, 1, 12 ) );
  const int leap = mm == 2 ? isLeapYear( yyyy ) : 0;
  const int result = daysPerMonth[ leap ][ mm - 1 ];
  POST0( IN_RANGE( result, 1, 31 ) );
  return result;
}



/******************************************************************************
PURPOSE: incrementHours - Increment yyyymmddhh by hours.
INPUTS:  const int yyyymmddhh   Timestamp to increment.
         const int hours        Number of hours to increment.
RETURNS: int incremented timestamp.
******************************************************************************/

int incrementHours( const int yyyymmddhh, const int hours ) {
  PRE02( isValidYYYYMMDDHH( yyyymmddhh ), hours >= 0 );
  int result = yyyymmddhh;

  if ( hours > 0 ) {
    int hour = 0;
    int yyyy = yyyymmddhh / 1000000;
    int mm   = yyyymmddhh / 10000 % 100;
    int dd   = yyyymmddhh / 100 % 100;
    int hh   = yyyymmddhh % 100;
    int leap = isLeapYear( yyyy );

    for ( hour = 0; hour < hours; ++hour ) {
      ++hh;

      if ( hh > 23 ) {
        hh = 0;
        ++dd;

        if ( dd > daysPerMonth[ leap ][ mm - 1 ] ) {
          dd = 1;
          ++mm;

          if ( mm > 12 ) {
            mm = 1;
            ++yyyy;
            leap = isLeapYear( yyyy );
          }
        }
      }
    }

    result = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
  }

  POST0( isValidYYYYMMDDHH( result ) );
  return result;
}



/******************************************************************************
PURPOSE: hoursInRange - Number of hours in range[yyyymmddhh1, yyyymmddhh2].
INPUTS:  const int yyyymmddhh1  First timestamp.
         const int yyyymmddhh1  Last  timestamp.
RETURNS: int number of hours in time range.
******************************************************************************/

int hoursInRange( const int yyyymmddhh1, const int yyyymmddhh2 ) {

  PRE03( isValidYYYYMMDDHH( yyyymmddhh1 ), isValidYYYYMMDDHH( yyyymmddhh2 ),
         yyyymmddhh1 <= yyyymmddhh2 );

  int result = 1;
  int yyyymmddhh = yyyymmddhh1;

  while ( yyyymmddhh < yyyymmddhh2 ) {
    yyyymmddhh = incrementHours( yyyymmddhh, 1 );
    ++result;
  }

  POST0( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: timestepsUntil - Timesteps from yyyymmddhh1 until yyyymmddhh2.
INPUTS:  const int yyyymmddhh1  First timestamp.
         const int yyyymmddhh1  Last  timestamp.
         const int hours        Hours per timestep.
RETURNS: int number of timesteps from first to last.
******************************************************************************/

int timestepsUntil( const int yyyymmddhh1,
                    const int yyyymmddhh2,
                    const int hours ) {

  PRE04( isValidYYYYMMDDHH( yyyymmddhh1 ), isValidYYYYMMDDHH( yyyymmddhh2 ),
         yyyymmddhh1 <= yyyymmddhh2, hours > 0 );

  int result = 0;
  int yyyymmddhh = yyyymmddhh1;

  while ( yyyymmddhh < yyyymmddhh2 ) {
    yyyymmddhh = incrementHours( yyyymmddhh, hours );
    ++result;
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: rotate4ByteArrayIfLittleEndian() - Rotate 4-bytes of each array item
         if on a little-endian platform.
INPUTS:  const size_t count  Number of items in array.
         void* array         Array of 4-byte values to rotate.
OUTPUTS: void* array         Array of rotated values.
******************************************************************************/

void rotate4ByteArrayIfLittleEndian( const size_t count, void* array ) {

#if IS_LITTLE_ENDIAN

  PRE02( count, array );
  assert_static( sizeof (int) == 4 );
  int* const array4 = array;
  long long index = 0; /* OpenMP requires a signed type for the loop index. */

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const int value = array4[ index ];
    const int newValue =
      ( value & 0xff000000 ) >> 24 |
      ( value & 0x00ff0000 ) >>  8 |
      ( value & 0x0000ff00 ) <<  8 |
      ( value & 0x000000ff ) << 24;
    array4[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: writeFloats() - Write big-endian float array to a file/stream.
INPUTS:  const size_t count    Number of items in array.
         float array[ count ]  Array of floats to write.
OUTPUTS: FILE* const output    File/stream to write to.
         float array[ count ]  Array of floats in big-endian format.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int writeFloats( const size_t count, float array[], FILE* const output ) {
  PRE03( count, array, output );
  int result = 0;
  rotate4ByteArrayIfLittleEndian( count, array );
  result = fwrite( array, 4, count, output ) == count;

  DEBUG( fprintf( stderr, "writeFloats( count = %lu ), result = %d\n",
                  count, result ); )

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write %lu bytes\n", 4 * count );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readFloats() - Read big-endian float array to a file/stream.
INPUTS:  const size_t count    Number of items in array.
OUTPUTS: FILE* const output    File/stream to read from.
         float array[ count ]  Array of floats in native format.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int readFloats( const size_t count, float array[], FILE* const output ) {
  PRE03( count, array, output );
  int result = 0;
  result = fread( array, 4, count, output ) == count;
  rotate4ByteArrayIfLittleEndian( count, array );

  DEBUG( fprintf( stderr, "readFloats( count = %lu ), result = %d\n",
                  count, result ); )

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read %lu bytes\n", 4 * count );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isDirectory - Determine if name is a directory.
INPUTS:  const char* name  Name of directory to examine.
RETURNS: int 1 if directory, else 0.
******************************************************************************/

int isDirectory( const char* name ) {
  PRE0( name );
  int result = 0;
  DIR* directory = 0;
  directory = opendir( name );
  result = directory != 0;

  if ( directory ) {
    closedir( directory );
    directory = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: fileDateUTC - UTC date of a file.
INPUTS:  const char* fileName  File to check.
RETURNS: int yyyymmdd of file.
******************************************************************************/

int fileDateUTC( const char* fileName ) {
  PRE02( fileName, *fileName );
  int result = 19000101;
  struct stat fileInfo;

  if ( stat( fileName, &fileInfo ) == 0 ) {
    const time_t seconds = fileInfo.st_mtime;
    struct tm timestamp;

    if ( gmtime_r( &seconds, &timestamp ) ) {
      const int yyyy = timestamp.tm_year + 1900;
      const int mm   = timestamp.tm_mon + 1;
      const int dd   = timestamp.tm_mday;
      const int yyyymmdd = yyyy * 10000 + mm * 100 + dd;

      if ( isValidYYYYMMDDHH( yyyymmdd * 100 ) ) {
        result = yyyymmdd;
      }
    }
  }

  POST0( isValidYYYYMMDDHH( result * 100 ) );
  return result;
}



/******************************************************************************
PURPOSE: fileSize - Determine size of named file.
INPUTS:  const char* name  Name of file to examine.
RETURNS: size_t size, in bytes, of named file, else 0 if failed.
******************************************************************************/

size_t fileSize( const char* name ) {
  PRE0( name );
  size_t result = 0;
  struct stat buf;

  if ( stat( name, &buf ) == -1 ) {
    fprintf( stderr, "\nFailed to determine size of file '%s'.\n", name );
  } else if ( buf.st_size > 0 ) {
    result = buf.st_size;
  }

  return result;
}



/******************************************************************************
PURPOSE: streamFile - Write bytes of named file to stdout.
INPUTS:  const char* name  Name of file to stream.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int streamFile( const char* const name ) {
  const size_t bufferSize = 1024 * 1024; /* 1MB. */
  void* buffer = NEW( char, bufferSize );
  int result = 0;

  if ( buffer ) {
    FILE* file = fopen( name, "rb" );

    if ( file ) {

      do {
        const size_t bytesRead = fread( buffer, 1, bufferSize, file );

        if ( bytesRead > 0 ) {
          result = fwrite( buffer, 1, bytesRead, stdout ) == bytesRead;
        }

       } while ( AND2( result, ! feof( file ) ) );

      fclose( file ), file = 0;
    }

    FREE( buffer );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: printWorkingDirectory - Print working directory to stdout.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int printWorkingDirectory( void ) {
  int result = 0;
  char path[ 256 ] = "";
  memset( path, 0, sizeof path );

  if ( getcwd( path, sizeof path / sizeof *path - 1 ) ) {
    result = puts( path ) != EOF;
  }

  return result;
}



/******************************************************************************
PURPOSE: printDirectoryListing - Print sub-directories & NetCDF files to stdout
INPUTS:  const char* name  Name of directory to list.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int printDirectoryListing( const char* const name ) {
  enum { NAME_SIZE = 256 };
  typedef char FileName[ NAME_SIZE ];
  enum { MAXIMUM_FILES = 4096 };
  static FileName listing[ MAXIMUM_FILES ]; /* 1MB. */
  DIR* directory = opendir( name );
  int result = 0;

  if ( directory ) {
#ifdef _WIN32
    const char slash = '\\';
#else
    const char slash = '/';
#endif
    struct dirent* entry = 0;
    size_t count = 0;
    memset( listing, 0, sizeof listing );
    result = 1;

    while ( count < MAXIMUM_FILES &&
            ( entry = readdir( directory ) ) != 0 ) {

      /* Only consider regular files, directories or links to these: */

#ifdef _WIN32
      /* On Windows, struct dirent lacks attribute d_type! */
      const int considerEntry = entry->d_name[ 0 ] != '.';
#else
      const int considerEntry = entry->d_name[ 0 ] != '.' &&
        ( entry->d_type == DT_DIR || entry->d_type == DT_REG ||
          entry->d_type == DT_LNK );
#endif

      if ( considerEntry ) {
        char path[ 256 ] = "";
        memset( path, 0, sizeof path );
        snprintf( path, sizeof path / sizeof *path,
                 "%s%c%s", name, slash, entry->d_name );

        {
#ifdef _WIN32
          const int isSubdirectory = isDirectory( path );
#else
          const int isSubdirectory =
            entry->d_type == DT_DIR || isDirectory( path );
#endif

          if ( isSubdirectory ) {
            snprintf( listing[ count ], sizeof *listing / sizeof (char),
                 "%s%c", entry->d_name, slash );
            ++count;
          } else if ( isNetCDFFile( path ) ) {
            strncpy( listing[ count ], entry->d_name,
                     sizeof *listing / sizeof (char) - 1 );
            ++count;
          }
        }
      }
    }

    closedir( directory ), directory = 0;

    /* Sort listing, ignoring case: */

    qsort( listing, count, sizeof *listing,
           (int (*)(const void*, const void*)) strcasecmp );

    /* Print listing: */

    printf( "..%c\n", slash ); /* Print parent directory first. */

    {
      size_t index = 0;

      for ( ; index < count; ++index ) {
        puts( listing[ index ] );
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: isNetCDFFile - Is file a NetCDF file?
INPUTS:  const char* name  Name of file to check.
RETURNS: int 1 if NetCDF file, else 0.
NOTES:   Does not check files that are on tape since this would be too slow.
         This routine only checks that the first 4 bytes of a file match one of
         CDF1, CDF2, 211HDF.
******************************************************************************/

int isNetCDFFile( const char* const name ) {
  int result = 0;
  struct stat buffer;
  memset( &buffer, 0, sizeof buffer );

  if ( stat( name, &buffer ) == 0 ) {

    if ( buffer.st_size > 10000 ) { /* Minimal CMAQ file header size. */
#ifdef _WIN32
      const int onTape = 0;
#else
      const int onTape = buffer.st_blocks == 0;
#endif

      if ( onTape ) { /* 0 blocks means on tape */
        result = 1; /* which is too slow to demigrate/read/check so allow. */
      } else { /* Read and check the first 4 bytes of file: */
        FILE* file = fopen( name, "rb" );

        if ( file ) {
          unsigned char bytes[ 5 ] = "";
          const int ok = fread( &bytes, 4, 1, file ) == 1;

          if ( ok ) {
            result =
              OR2( AND4( bytes[ 0 ] == 'C',
                         bytes[ 1 ] == 'D',
                         bytes[ 2 ] == 'F',
                         IN3( bytes[ 3 ], 1, 2 ) ),
                   AND4( bytes[ 0 ] == 0211,
                         bytes[ 1 ] == 'H',
                         bytes[ 2 ] == 'D',
                         bytes[ 3 ] == 'F'));
          }

          fclose( file ), file = 0;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: isValidBounds - Check validity of bounds object.
INPUTS:  const Bounds bounds  Object to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

int isValidBounds( const Bounds bounds) {
  const int result =
    AND5( bounds != 0,
          IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ),
          IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
                    bounds[ LONGITUDE ][ MINIMUM ], 180.0 ),
          IN_RANGE( bounds[ LATITUDE  ][ MINIMUM ], -90.0, 90.0 ),
          IN_RANGE( bounds[ LATITUDE  ][ MAXIMUM ],
                    bounds[ LATITUDE  ][ MINIMUM ], 90.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: boundsOverlap - Do rectangles overlap?
INPUTS:  const Bounds a  a[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
         const Bounds b  b[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
RETURNS: int 1 if overlap, else 0.
******************************************************************************/

int boundsOverlap( const Bounds a, const Bounds b ) {
  PRE02( isValidBounds( a ), isValidBounds( b ) );
  const int outside =
    OR4( a[ LATITUDE  ][ MINIMUM ] > b[ LATITUDE  ][ MAXIMUM ],
         a[ LATITUDE  ][ MAXIMUM ] < b[ LATITUDE  ][ MINIMUM ],
         a[ LONGITUDE ][ MINIMUM ] > b[ LONGITUDE ][ MAXIMUM ],
         a[ LONGITUDE ][ MAXIMUM ] < b[ LONGITUDE ][ MINIMUM ] );
  const int result = ! outside;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: boundsSubsumes - Does a subsume b. I.e., is b completely inside a?
INPUTS:  const Bounds a  a[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
         const Bounds b  b[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
RETURNS: int 1 if a subsumes b, else 0.
******************************************************************************/

int boundsSubsumes( const Bounds a, const Bounds b ) {
  PRE02( isValidBounds( a ), isValidBounds( b ) );
  double minimum    = a[ LONGITUDE ][ MINIMUM ];
  double maximum    = a[ LONGITUDE ][ MAXIMUM ];
  double coordinate = b[ LONGITUDE ][ MINIMUM ];
  int result = IN_RANGE( coordinate, minimum, maximum );

  if ( result ) {
    coordinate = b[ LONGITUDE ][ MAXIMUM ];
    result = IN_RANGE( coordinate, minimum, maximum );

    if ( result ) {
      minimum    = a[ LATITUDE ][ MINIMUM ];
      maximum    = a[ LATITUDE ][ MAXIMUM ];
      coordinate = b[ LATITUDE ][ MINIMUM ];
      result = IN_RANGE( coordinate, minimum, maximum );

      if ( result ) {
        coordinate = b[ LATITUDE ][ MAXIMUM ];
        result = IN_RANGE( coordinate, minimum, maximum );
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: areaOfTriangle - Area of triangle with vertices
         (x1, y1), (x2, y2), (x3, y3).
INPUTS:  double x1   X-Coordinate of 1st vertex of triangle.
         double y1   Y-Coordinate of 1st vertex of triangle.
         double x2   X-Coordinate of 2nd vertex of triangle.
         double y2   Y-Coordinate of 2nd vertex of triangle.
         double x3   X-Coordinate of 3rd vertex of triangle.
         double y3   Y-Coordinate of 3rd vertex of triangle.
RETURNS: double area.
******************************************************************************/

double areaOfTriangle( double x1, double y1,
                       double x2, double y2,
                       double x3, double y3 ) {
  const double a = x1 - x3;
  const double b = y1 - y3;
  const double c = x2 - x3;
  const double d = y2 - y3;
  const double triangleArea = 0.5 * ( a * d - b * c );
  const double result = triangleArea < 0.0 ? -triangleArea : triangleArea;
  return result;
}



/******************************************************************************
PURPOSE: clipPolygon - Clip polygon to an axis-aligned rectangle
         and return the number of vertices in clipped polygon.
INPUTS:  const int discardDegenerates  Discard degenerate triangles?
                                       Slower, but useful if rendering.
         const double clipXMin  X-coordinate of lower-left  corner of clip rect
         const double clipYMin  Y-coordinate of lower-left  corner of clip rect
         const double clipXMax  X-coordinate of upper-right corner of clip rect
         const double clipYMax  Y-coordinate of upper-right corner of clip rect
         const int count          Number of input vertices.
         const double x[ count ]  X-coordinates of input polygon to clip.
         const double y[ count ]  Y-coordinates of input polygon to clip.
         double cx[ 2 * count + 2 ] Storage for 2 * count + 2 X-coordinates.
         double cy[ 2 * count + 2 ] Storage for 2 * count + 2 Y-coordinates.
OUTPUTS: double cx[ result ]      X-coordinates of vertices of clipped poly.
         double cy[ result ]      Y-coordinates of vertices of clipped poly.
RETURNS: int number of vertices in clipped polygon.
NOTES:   Uses the Liang-Barsky polygon clipping algorithm. (Fastest known.)
         "An Analysis and Algorithm for Polygon Clipping",
         You-Dong Liang and Brian Barsky, UC Berkeley,
         CACM Vol 26 No. 11, November 1983.
         https://www.longsteve.com/fixmybugs/?page_id=210
******************************************************************************/

int clipPolygon( const int discardDegenerates,
                 const double clipXMin,
                 const double clipYMin,
                 const double clipXMax,
                 const double clipYMax,
                 const int count,
                 const double x[],
                 const double y[],
                 double cx[],
                 double cy[] ) {
  int result = 0;
  const double inf = DBL_MAX;
  double xIn   = 0.0; /* X-coordinate of entry point. */
  double yIn   = 0.0; /* Y-coordinate of entry point. */
  double xOut  = 0.0; /* X-coordinate of exit point. */
  double yOut  = 0.0; /* Y-coordinate of exit point. */
  double tInX  = 0.0; /* Parameterized X-coordinate of entry intersection. */
  double tInY  = 0.0; /* Parameterized Y-coordinate of entry intersection. */
  double tOutX = 0.0; /* Parameterized X-coordinate of exit intersection. */
  double tOutY = 0.0; /* Parameterized Y-coordinate of exit intersection. */
  int vertex = 0;

  for ( vertex = 0; vertex < count; ++vertex ) {
    const int vertexp1 = vertex + 1;
    const int vertex1 = vertexp1 < count ? vertexp1 : 0;
    const double vx = x[ vertex ];
    const double vy = y[ vertex ];
    const double deltaX = x[ vertex1 ] - vx; /* Edge direction. */
    const double deltaY = y[ vertex1 ] - vy;
    const double oneOverDeltaX = deltaX ? 1.0 / deltaX : 0.0;
    const double oneOverDeltaY = deltaY ? 1.0 / deltaY : 0.0;
    double tOut1 = 0.0;
    double tOut2 = 0.0;
    CHECK( result < count + count + 2 );

    /*
     * Determine which bounding lines for the clip window the containing line
     * hits first:
     */

    if ( deltaX > 0.0 || ( deltaX == 0.0 && vx > clipXMax ) ) {
      xIn  = clipXMin;
      xOut = clipXMax;
    } else {
      xIn  = clipXMax;
      xOut = clipXMin;
    }

    if ( deltaY > 0.0 || ( deltaY == 0.0 && vy > clipYMax ) ) {
      yIn  = clipYMin;
      yOut = clipYMax;
    } else {
      yIn  = clipYMax;
      yOut = clipYMin;
    }

    /* Find the t values for the x and y exit points: */

    if ( deltaX != 0.0 ) {
      tOutX = ( xOut - vx ) * oneOverDeltaX;
    } else if ( vx <= clipXMax && clipXMin <= vx ) {
      tOutX = inf;
    } else {
      tOutX = -inf;
    }

    if ( deltaY != 0.0 ) {
      tOutY = ( yOut - vy ) * oneOverDeltaY;
    } else if ( vy <= clipYMax && clipYMin <= vy ) {
      tOutY = inf;
    } else {
      tOutY = -inf;
    }

    /* Set tOut1 = min( tOutX, tOutY ) and tOut2 = max( tOutX, tOutY ): */

    if ( tOutX < tOutY ) {
      tOut1 = tOutX;
      tOut2 = tOutY;
    } else {
      tOut1 = tOutY;
      tOut2 = tOutX;
    }

    if ( tOut2 > 0.0 ) {
      double tIn2 = 0.0;

      if ( deltaX != 0.0 ) {
        tInX = ( xIn - vx ) * oneOverDeltaX;
      } else {
        tInX = -inf;
      }

      if ( deltaY != 0.0 ) {
        tInY = ( yIn - vy ) * oneOverDeltaY;
      } else {
        tInY = -inf;
      }

      /* Set tIn2 = max( tInX, tInY ): */

      if ( tInX < tInY ) {
        tIn2 = tInY;
      } else {
        tIn2 = tInX;
      }

      if ( tOut1 < tIn2 ) { /* No visible segment. */

        if ( 0.0 < tOut1 && tOut1 <= 1.0 ) {
          CHECK( result < count + count );

          /* Line crosses over intermediate corner region. */

          if ( tInX < tInY ) {
            cx[ result ] = xOut;
            cy[ result ] = yIn;
          } else {
            cx[ result ] = xIn;
            cy[ result ] = yOut;
          }

          ++result;
        }
      } else { /* Line crosses through window: */

        if ( 0.0 < tOut1 && tIn2 <= 1.0 ) {

          if ( 0.0 <= tIn2 ) { /* Visible segment: */
            CHECK( result < count + count );

            if ( tInX > tInY ) {
              cx[ result ] = xIn;
              cy[ result ] = vy + ( tInX * deltaY );
            } else {
              cx[ result ] = vx + ( tInY * deltaX );
              cy[ result ] = yIn;
            }

            ++result;
          }

          CHECK( result < count + count );

          if ( 1.0 >= tOut1 ) {

            if ( tOutX < tOutY ) {
              cx[ result ] = xOut;
              cy[ result ] = vy + ( tOutX * deltaY );
            } else {
              cx[ result ] = vx + ( tOutY * deltaX );
              cy[ result ] = yOut;
            }

            ++result;
          } else {
            cx[ result ] = x[ vertex1 ];
            cy[ result ] = y[ vertex1 ];
            ++result;
          }
        }
      }

      if ( 0.0 < tOut2 && tOut2 <= 1.0 ) {
        CHECK( result < count + count );
        cx[ result ] = xOut;
        cy[ result ] = yOut;
        ++result;
      }
    }
  }

  /*
   * The above algorithm can generate 5-vertex 'line' or 'hat' polygons: _/\_
   * where the last 3 vertices are colinear
   * which yields a degenerate 'triangle' (i.e., with 0 area).
   * Here we discard the last 2 verticies in such cases.
   */

  if ( discardDegenerates && result == 5 ) {
    int twice = 2; /* Check twice in case of 5-vertex 'line'. */

    do {

      if ( result >= 3 ) {
        const size_t count_3 = result - 3;
        const size_t count_2 = result - 2;
        const size_t count_1 = result - 1;

        const double lastTriangleArea =
          areaOfTriangle( cx[ count_3 ], cy[ count_3 ],
                          cx[ count_2 ], cy[ count_2 ],
                          cx[ count_1 ], cy[ count_1 ] );

        DEBUG( fprintf( stderr, "lastTriangleArea = %f\n", lastTriangleArea );)

        if ( lastTriangleArea == 0.0 ) {
          result -= 2;
        }
      }
    } while (--twice );
  }

  /* Always discard any result less than a triangle. */

  if ( result < 3 ) {
    result = 0;
  }

  POST0( OR2( result == 0, IN_RANGE( result, 3, count + count + 2 ) ) );
  return result;
}




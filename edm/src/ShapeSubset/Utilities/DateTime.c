/******************************************************************************
PURPOSE: DateTime.c - Define routines for date/time computation.
NOTES:   See DataTime.h.
HISTORY: 2004/10, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For stderr, snprintf(). */
#include <stdlib.h> /* For abs(). */
#include <string.h> /* For strlen(), memset(), memcpy(). */
#include <time.h>   /* For gmtime_r(). */
#include <sys/stat.h> /* For stat(). */

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

#include <Assertions.h> /* For PRE0*(), POST0*(), AND?(). */
#include <DateTime.h> /* For public interface. */

/*============================= PUBLIC FUNCTIONS ============================*/


/*
 * 30 days hath September, April, June and November, all the rest have 31,
 * except February which has either 28 or 29 (on a leap year).
 */

static const Integer daysPerMonth[ 2 ][ 12 ] =
{
  { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
  { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
};


/*============================= PUBLIC FUNCTIONS ============================*/




/******************************************************************************
PURPOSE: isValidDate - Is the given date valid YYYYDDD format?
INPUTS:  Integer yyyyddd  The date to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidDate( Integer yyyyddd ) {
  const Integer yyyy = yyyyddd / 1000;
  const Integer ddd  = yyyyddd % 1000;
  const Integer result = AND3( yyyy >= 1800, IN_RANGE( ddd, 1, 366 ),
                               IMPLIES( ddd == 366, isLeapYear( yyyy ) ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidTime - Is the given time valid HHMMSS format?
INPUTS:  Integer hhmmss  The time to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidTime( Integer hhmmss ) {
  const Integer hh     = hhmmss / 10000;
  const Integer mm     = ( hhmmss / 100 ) % 100;
  const Integer ss     = hhmmss % 100;
  const Integer result =
    AND3( IN_RANGE( hh, 0, 23 ), IN_RANGE( mm, 0, 59 ), IN_RANGE( ss, 0, 59));
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidTimestepSize - Is the given time valid *HHMMSS format?
INPUTS:  Integer hhmmss  The timestep size to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidTimestepSize( Integer hhmmss ) {
  const Integer hh     = hhmmss / 10000;
  const Integer mm     = ( hhmmss / 100 ) % 100;
  const Integer ss     = hhmmss % 100;
  const Integer result =
    AND4( hhmmss > 0, hh >= 0, IN_RANGE( mm, 0, 59 ), IN_RANGE( ss, 0, 59 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isLeapYear - Does the year have 366 days?
INPUTS:  Integer yyyy  YYYY.
RETURNS: Integer 1 if yyyy is a leap year, else 0.
******************************************************************************/

Integer isLeapYear( Integer yyyy ) {
  PRE0( isValidDate( yyyy * 1000 + 1 ) );
  const Integer result = yyyy % 4 == 0 && (yyyy % 100 != 0 || yyyy % 400 == 0);
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: timestepsInRange - Number of one-hour timesteps in range, inclusive.
INPUTS:  Integer firstDate  YYYYDDD
         Integer firstTime  HHMMSS
         Integer lastDate   YYYYDDD
         Integer lastTime   HHMMSS
RETURNS: Integer Number of one-hour timesteps in range, inclusive.
******************************************************************************/

Integer timestepsInRange( Integer firstDate, Integer firstTime,
                          Integer lastDate,  Integer lastTime ) {

  PRE06( isValidDate( firstDate ), isValidTime( firstTime ),
         isValidDate( lastDate  ), isValidTime( lastTime  ),
         firstDate <= lastDate,
         IMPLIES( firstDate == lastDate, firstTime <= lastTime ) );

  Integer yyyyddd = firstDate;
  Integer hhmmss  = firstTime;
  Integer result  = 1;

  while ( ! AND2( yyyyddd == lastDate, hhmmss == lastTime ) ) {
    incrementOneHour( &yyyyddd, &hhmmss );
    ++result;
  }

  POST0( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: monthAndDay - Month [1..12] and day of month [1..31] of yyyyddd.
INPUTS:  Integer yyyyddd      Year and day of year.
OUTPUTS: Integer* month       Month of year [1..12].
         Integer* dayOfMonth  Day of month [1..31].
******************************************************************************/

void monthAndDay( Integer yyyyddd, Integer* month, Integer* dayOfMonth ) {

  PRE03( isValidDate( yyyyddd ), month, dayOfMonth );

  enum { MONTHS = 12 };
  Integer daysPerMonth[ MONTHS ] = {
    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
  };

  const Integer yyyy = yyyyddd / 1000;
  const Integer ddd  = yyyyddd % 1000;
  Integer dddRemaining = ddd;
  Integer days = 0, thisMonth = 0;
  *month = *dayOfMonth = 0;
  daysPerMonth[ 1 ] += isLeapYear( yyyy );

  for ( thisMonth = 0; thisMonth < MONTHS; ++thisMonth ) {
    const Integer daysThisMonth = daysPerMonth[ thisMonth ];
    days += daysThisMonth;

    if ( days >= ddd ) {
      *month = thisMonth + 1;
      *dayOfMonth = dddRemaining;
      thisMonth = MONTHS;
    }

    dddRemaining -= daysThisMonth;
  }

  POST02( IN_RANGE( *month, 1, 12 ), IN_RANGE( *dayOfMonth, 1, 31 ) );
}



/******************************************************************************
PURPOSE: incrementOneHour - Increment date/time by one hour.
INPUTS:  Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
OUTPUTS: Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
******************************************************************************/

void incrementOneHour( Integer* yyyyddd, Integer* hhmmss ) {

  PRE04( yyyyddd, hhmmss, isValidDate( *yyyyddd ), isValidTime( *hhmmss ) );

  const Integer oneHour     =  10000; /* HHMMSS = 01:00:00. */
  const Integer maximumTime = 235959; /* HHMMSS = 23:59:59. */
  *hhmmss += oneHour;

  if ( *hhmmss > maximumTime ) {
    const Integer ss = *hhmmss % 100;
    const Integer mm = ( *hhmmss / 100 ) % 100;
    Integer ddd;
    *hhmmss = mm * 100 + ss;
    *yyyyddd += 1;
    ddd = *yyyyddd % 1000;

    if ( ddd > 365 ) {
      const Integer yyyy = *yyyyddd / 1000;
      const Integer daysInYear = 365 + isLeapYear( yyyy );

      if ( ddd > daysInYear ) {
        *yyyyddd = ( yyyy + 1 ) * 1000 + 1; /* Next year, first day. */
      }
    }
  }

  POST02( isValidDate( *yyyyddd ), isValidTime( *hhmmss ) );
}



/******************************************************************************
PURPOSE: decrementOneHour - Decrement date/time by one hour.
INPUTS:  Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
OUTPUTS: Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
******************************************************************************/

void decrementOneHour( Integer* yyyyddd, Integer* hhmmss ) {

  PRE04( yyyyddd, hhmmss, isValidDate( *yyyyddd ), isValidTime( *hhmmss ) );

  const Integer oneHour = 10000; /* HHMMSS = 01:00:00. */
  *hhmmss -= oneHour;

  if ( *hhmmss < 0 ) {
    const Integer lastHour = 230000; /* HHMMSS = 23:00:00. */
    Integer ss, mm, ddd;
    *hhmmss += oneHour;
    ss = *hhmmss % 100;
    mm = ( *hhmmss / 100 ) % 100;
    *hhmmss = lastHour + mm * 100 + ss;
    *yyyyddd -= 1;
    ddd = *yyyyddd % 1000;

    if ( ddd < 1 ) {
      const Integer yyyy = *yyyyddd / 1000;
      const Integer daysInYear = 365 + isLeapYear( yyyy );
      *yyyyddd = ( yyyy - 1 ) * 1000 + daysInYear; /* Previous year last day*/
    }
  }

  POST02( isValidDate( *yyyyddd ), isValidTime( *hhmmss ) );
}



/******************************************************************************
PURPOSE: incrementTime - Increment date/time by step.
INPUTS:  Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
         Integer  step     ...HHMMSS.
OUTPUTS: Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
******************************************************************************/

void incrementTime( Integer* yyyyddd, Integer* hhmmss, Integer step ) {

  PRE05( yyyyddd, hhmmss, isValidDate( *yyyyddd ), isValidTime( *hhmmss ),
         isValidTimestepSize( step ) );

  Integer hours = step / 10000; /* ...hh0000 */

  while ( hours-- ) {
    incrementOneHour( yyyyddd, hhmmss );
  }

  step %= 10000; /* 00mmss. */

  if ( step ) {
    const Integer stepSS = step % 100;
    Integer ss = *hhmmss % 100;
    Integer mm = ( *hhmmss / 100 ) % 100;
    Integer hh = 0;
    ss += stepSS;
    mm += ss / 60;
    ss %= 60;
    hh = mm / 60;
    mm %= 60;

    if ( hh > 0 ) {
      incrementOneHour( yyyyddd, hhmmss );
      hh = *hhmmss / 10000;
    }

    *hhmmss = hh * 10000 + mm * 100 + ss;
  }

  POST02( isValidDate( *yyyyddd ), isValidTime( *hhmmss ) );
}



/******************************************************************************
PURPOSE: decrementTime - Decrement date/time by step.
INPUTS:  Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
         Integer  step     HHMMSS.
OUTPUTS: Integer* yyyyddd  YYYYDDD.
         Integer* hhmmss   HHMMSS.
******************************************************************************/

void decrementTime( Integer* yyyyddd, Integer* hhmmss, Integer step ) {

  PRE05( yyyyddd, hhmmss, isValidDate( *yyyyddd ), isValidTime( *hhmmss ),
         isValidTimestepSize( step ) );

  Integer hours = 1 + step / 10000; /* ...hh0000 */

  do {
    decrementOneHour( yyyyddd, hhmmss );
    --hours;
  } while ( hours );

  step %= 10000; /* 00mmss. */

  if ( step == 0 ) {
    incrementOneHour( yyyyddd, hhmmss );
  } else {
    incrementTime( yyyyddd, hhmmss, step );
  }

  POST02( isValidDate( *yyyyddd ), isValidTime( *hhmmss ) );
}



/******************************************************************************
PURPOSE: isValidUTCTimestamp - Is string a valid UTCTimestamp?
INPUTS:  const char* string  The string to examine.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidUTCTimestamp( const char* string ) {
  Integer result = 0;

  if ( AND2( string, strlen( string ) == UTC_TIMESTAMP_LENGTH ) ) {
    int yyyy = 0, mo = 0, dd = 0, hh = 0, mm = 0, ss = 0, zone = 0;

    if ( sscanf( string, "%04d-%02d-%02dT%02d:%02d:%02d-%04d",
                 &yyyy, &mo, &dd, &hh, &mm, &ss, &zone ) == 7 ) {
      const Integer yyyymmdd = yyyy * 10000 + mo * 100 + dd;
      result = isValidYearMonthDay( yyyymmdd );
      result = AND6( result,
                     IN_RANGE( hh, 0, 23 ),
                     IN_RANGE( mm, 0, 59 ),
                     IN_RANGE( ss, 0, 59 ),
                     IN_RANGE( zone / 100, -23, 23 ),
                     IN_RANGE( abs( zone % 100 ), 0, 59 ) );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: toUTCTimestamp - Convert timestamp to ISO UTC string format.
INPUTS:  Integer yyyydddhhmm  Timestamp to convert.
OUTPUTS: UTCTimestamp string  ISO UTC string format for the timestamp.
******************************************************************************/

void toUTCTimestamp( Integer yyyydddhhmm, UTCTimestamp string ) {
  PRE02( isValidTimestamp( yyyydddhhmm ), string );
  const Integer mm      = yyyydddhhmm % 100;
  const Integer hh      = yyyydddhhmm / 100 % 100;
  const Integer yyyyddd = yyyydddhhmm / 10000;
  const Integer yyyy    = yyyyddd / 1000;
  Integer mo = 0;
  Integer dd = 0;
  monthAndDay( yyyyddd, &mo, &dd );
  snprintf( string, sizeof (UTCTimestamp) / sizeof (char),
            "%04"INTEGER_FORMAT"-%02"INTEGER_FORMAT
            "-%02"INTEGER_FORMAT
            "T%02"INTEGER_FORMAT
            ":%02"INTEGER_FORMAT":00-0000",
            yyyy, mo, dd, hh, mm );
  POST0( strlen( string ) == UTC_TIMESTAMP_LENGTH );
}



/******************************************************************************
PURPOSE: fromUTCTimestamp - Convert ISO UTC string to integer.
INPUTS:  const UTCTimestamp string  ISO UTC string to convert.
RETURNS: Integer yyyydddhhmm  Integer timestamp.
******************************************************************************/

Integer fromUTCTimestamp( const UTCTimestamp string ) {

  PRE0( isValidUTCTimestamp( (const char*) string ) );

  Integer result = 0;
  int yyyy = 1900, mo = 1, dd = 1, hh = 0, mm = 0, ss = 0, zone = 0;

  if ( sscanf( string, "%04d-%02d-%02dT%02d:%02d:%02d-%04d",
               &yyyy, &mo, &dd, &hh, &mm, &ss, &zone ) == 7 ) {
    result = yyyy;
    result *= 100;
    result += mo;
    result *= 100;
    result += dd;
    result = convertYearMonthDay( result );
    result *= 100;
    result += hh;
    result *= 100;
    result += mm;
  }

  POST0( isValidTimestamp( result ) );
  return result;
}



/******************************************************************************
PURPOSE: parseTimestamp - Parse string timestamp into its integer value.
INPUTS:  const char* string      String representation of timestamp.
OUTPUTS: Integer*    yyyydddhh00 Value of string as yyyydddhh00.
RETURNS: Integer 1 if successful, else 0.
******************************************************************************/

Integer parseTimestamp( const char* string, Integer* yyyydddhh00 ) {

  PRE02( string, yyyydddhh00 );

  const Integer yyyymmddhh = atoI( string );
  const Integer yyyymmdd       = yyyymmddhh / 100;
  const Integer hh             = yyyymmddhh % 100;
  const Integer result =
    AND2( IN_RANGE( hh, 0, 23 ), isValidYearMonthDay( yyyymmdd ) );

  *yyyydddhh00 = 0;

  if ( ! result ) {
    fprintf( stderr, "\a\n\nInvalid timestamp specified '%s'.\n", string );
  } else {
    const Integer yyyyddd   = convertYearMonthDay( yyyymmdd );
    const Integer yyyydddhh = yyyyddd * 100 + hh;
    *yyyydddhh00 = yyyydddhh * 100;
  }

  POST02( IS_BOOL( result ),
         IMPLIES_ELSE( result,
                       isValidTimestamp( *yyyydddhh00 ), *yyyydddhh00 == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidTimestamp - Is the timestamp valid?
INPUTS:  Integer yyyydddhhmm The timestamp to examine.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidTimestamp( Integer yyyydddhhmm ) {
  const Integer yyyy = yyyydddhhmm / 10000000;
  const Integer ddd  = yyyydddhhmm / 10000 % 1000;
  const Integer hh   = yyyydddhhmm / 100 % 100;
  const Integer mm   = yyyydddhhmm % 100;
  const Integer result =
    AND4( IN_RANGE( yyyy, 1800, 9999 ),
          IN_RANGE( ddd, 1, 365 + isLeapYear( yyyy ) ),
          IN_RANGE( hh, 0, 23 ),
          IN_RANGE( mm, 0, 59 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidYearMonthDay - Is the date valid?
INPUTS:  Integer yyyymmdd The date to examine.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidYearMonthDay( Integer yyyymmdd ) {
  const Integer yyyy = yyyymmdd / 10000;
  const Integer mm   = yyyymmdd / 100 % 100;
  const Integer dd   = yyyymmdd % 100;
  const Integer result =
    AND3( IN_RANGE( yyyy, 1800, 9999 ),
          IN_RANGE( mm, 1, 12 ),
          IN_RANGE( dd, 1, daysPerMonth[ isLeapYear( yyyy ) ][ mm - 1 ] ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: convertYearMonthDay - Convert date from YYYYMMDD to YYYYDDD.
INPUTS:  Integer yyyymmdd  Year, month, day to convert.
RETURNS: Integer yyyyddd equivalent day.
******************************************************************************/

Integer convertYearMonthDay( Integer yyyymmdd ) {

  PRE03( yyyymmdd / 10000 > 1000,
        IN_RANGE(  yyyymmdd / 100 % 100, 1, 12 ),
        IN_RANGE(  yyyymmdd % 100, 1,
                   daysPerMonth[ isLeapYear( yyyymmdd / 10000 ) ]
                   [ yyyymmdd / 100 % 100 - 1] ) );

  const Integer yyyy = yyyymmdd / 10000;
  const Integer mm0  = yyyymmdd / 100 % 100 - 1;
  const Integer dd   = yyyymmdd % 100;
  const Integer leap = isLeapYear( yyyy );
  Integer result = 0;
  Integer month = 0;

  for ( month = 0; month < mm0; ++month ) {
    result += daysPerMonth[ leap ][ month ];
  }

  result += dd;
  result += yyyy * 1000;

  POST03( result / 1000 == yyyymmdd / 10000,
          result % 1000 > 0,
          result % 1000 <= 365 + isLeapYear( yyyymmdd / 10000 ) );
  return result;
}



/******************************************************************************
PURPOSE: incrementTimestamp - Increment timestamp by one hour.
INPUTS:  Integer* yyyydddhhmm  Timestamp to increment.
OUTPUTS: Integer* yyyydddhhmm  Inremented timestamp.
******************************************************************************/

void incrementTimestamp( Integer* yyyydddhhmm ) {
  PRE02( yyyydddhhmm, isValidTimestamp( *yyyydddhhmm ) );
  const Integer mm = *yyyydddhhmm % 100;
  Integer hh = *yyyydddhhmm / 100 % 100;
  ++hh;

  if ( hh < 24 ) { /* Just update the hour: */
    *yyyydddhhmm = *yyyydddhhmm / 10000 * 10000 + hh * 100 + mm;
  } else { /* Reset the hour and increment the day: */
    Integer yyyy = *yyyydddhhmm / 10000000;
    Integer ddd  = *yyyydddhhmm / 10000 % 1000;
    hh = 0;
    ++ddd;

    if ( AND2( ddd > 365, ddd > 365 + isLeapYear( yyyy ) ) ) {
      ddd = 1; /* First day,    */
      ++yyyy;  /* of next year. */
    }

    *yyyydddhhmm = yyyy * 10000000 + ddd * 10000 + hh * 00 + mm;
  }

  POST0( isValidTimestamp( *yyyydddhhmm ) );
}



/******************************************************************************
PURPOSE: decrementTimestamp - Decrement timestamp by one hour.
INPUTS:  Integer* yyyydddhhmm  Timestamp to increment.
OUTPUTS: Integer* yyyydddhhmm  Deremented timestamp.
******************************************************************************/

void decrementTimestamp( Integer* yyyydddhhmm ) {
  PRE02( yyyydddhhmm, isValidTimestamp( *yyyydddhhmm ) );
  const Integer mm = *yyyydddhhmm % 100;
  Integer hh = *yyyydddhhmm / 100 % 100;
  --hh;

  if ( hh >= 0 ) { /* Just update the hour: */
    *yyyydddhhmm = *yyyydddhhmm / 10000 * 10000 + hh * 100 + mm;
  } else { /* Reset the hour and decrement the day: */
    Integer yyyy = *yyyydddhhmm / 10000000;
    Integer ddd  = *yyyydddhhmm / 10000 % 1000;
    hh = 23;
    --ddd;

    if ( ddd == 0 ) {
      --yyyy; /* Previous year. */
      ddd = 365 + isLeapYear( yyyy );
    }

    *yyyydddhhmm = yyyy * 10000000 + ddd * 10000 + hh * 100 + mm;
  }

  POST0( isValidTimestamp( *yyyydddhhmm ) );
}



/******************************************************************************
PURPOSE: offsetTimestamp - Compute timestamp + offset.
INPUTS:  Integer yyyydddhhmm  Initial timestamp.
         Integer hours        + or - hours to offset.
RETURNS: Integer yyyydddhhmm + hours.
******************************************************************************/

Integer offsetTimestamp( Integer yyyydddhhmm, Integer hours ) {
  PRE0( isValidTimestamp( yyyydddhhmm ) );
  Integer result = yyyydddhhmm;
  Integer hour = hours;

  while ( hour < 0 ) {
    decrementTimestamp( &result );
    ++hour;
  }

  while ( hour > 0 ) {
    incrementTimestamp( &result );
    --hour;
  }

  POST02( isValidTimestamp( result ),
          OR3( AND2( hours == 0, result == yyyydddhhmm ),
               AND2( hours < 0,  result <  yyyydddhhmm ),
               AND2( hours > 0,  result >  yyyydddhhmm ) ) );
  return result;
}



/******************************************************************************
PURPOSE: nowUTC - Current timestamp in UTC.
RETURNS: Integer yyyydddhhmm.
******************************************************************************/

Integer nowUTC( void ) {
  const time_t clock = time( 0 );
  struct tm timeInfo;
  Integer result = INTEGER_CONSTANT( 19000010000 );

  if ( gmtime_r( &clock, &timeInfo ) ) {
    result = timeInfo.tm_year + 1900;
    result = result * 1000 + timeInfo.tm_yday + 1;
    result = result *  100 + timeInfo.tm_hour;
    result = result *  100 + timeInfo.tm_min;
  }

  POST0( isValidTimestamp( result ) );
  return result;
}




/******************************************************************************
PURPOSE: fileDateUTC - UTC date when named file was last modified.
INPUTS:  const char* fileName  Name of file to check.
RETURNS: int yyyymmdd.
******************************************************************************/

int fileDateUTC( const char* fileName ) {
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

      if ( isValidYearMonthDay( yyyymmdd ) ) {
        result = yyyymmdd;
      }
    }
  }

  POST0( isValidYearMonthDay( result ) );
  return result;
}



/******************************************************************************
PURPOSE: timeZoneOffset - Hour offset (subtracted) from UTC.
INPUTS:  const char* timeZoneName  E.g., "EST".
RETURNS: int offset from UTC (e.g., 4).
******************************************************************************/

int timeZoneOffset( const char* timeZoneName ) {
  PRE02( timeZoneName, *timeZoneName );
  typedef struct { const char* name; int offset; } TimeZoneEntry;
  static const TimeZoneEntry timeZones[] = {
    { "AST", 4 },
    { "EST", 5 },
    { "EDT", 4 },
    { "CST", 6 },
    { "CDT", 5 },
    { "MST", 7 },
    { "MDT", 6 },
    { "PST", 8 },
    { "PDT", 7 },
    { "AKST", 9 },
    { "AKDT", 8 },
    { "HAST", 10 },
    { "HASDT", 9 },
    { 0, 0 }
  };
  int result = 0;
  int index = 0;

  for ( index = 0; AND2( ! result, timeZones[ index ].name ); ++index ) {

    if ( ! strcmp( timeZoneName, timeZones[ index ].name ) ) {
      result = timeZones[ index ].offset;
    }
  }

  POST0( IN_RANGE( result, 0, 23 ) );
  return result;
}


/******************************************************************************
PURPOSE: DateTime.c - Define routines for date/time computation.
NOTES:   See DataTime.h.
HISTORY: 2004/10, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For stderr, sprintf(). */
#include <stdlib.h> /* For abs(). */
#include <string.h> /* For strlen(). */
#include <time.h>   /* For gmtime_r(). */

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
  const Integer result = AND3( yyyy >= 1950, IN_RANGE( ddd, 1, 366 ),
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
PURPOSE: daysInMonth - Number of days in month of year.
INPUTS:  Integer year  YYYY.
         Integer month 1..12.
RETURNS: Integer number of days in month of year.
******************************************************************************/

Integer daysInMonth( Integer year, Integer month ) {
  PRE02( isValidDate( year * 1000 + 1 ), IN_RANGE( month, 1, 12 ) );
  const Integer leap = isLeapYear( year );
  const Integer result = daysPerMonth[ leap ][ month - 1 ];
  POST0( IN_RANGE( result, 1, 31 ) );
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
  sprintf( string, "%04"INTEGER_FORMAT"-%02"INTEGER_FORMAT
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
    AND4( IN_RANGE( yyyy, 1900, 9999 ),
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
    AND3( IN_RANGE( yyyy, 1900, 9999 ),
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
PURPOSE: fromSeconds - Current timestamp in UTC from seconds since 1970.
INPUTS:  const Integer seconds   Seconds since 1970.
RETURNS: Integer yyyymmddhhmmss.
******************************************************************************/

Integer fromSeconds( const Integer seconds ) {
  const time_t clock = (time_t) seconds;
  struct tm timeInfo;
  Integer result = INTEGER_CONSTANT( 19000101000000 );

  if ( gmtime_r( &clock, &timeInfo ) ) {
    result = timeInfo.tm_year + 1900;
    result = result * 100 + timeInfo.tm_mon + 1;
    result = result * 100 + timeInfo.tm_mday;
    result = result * 100 + timeInfo.tm_hour;
    result = result * 100 + timeInfo.tm_min;
    result = result * 100 + timeInfo.tm_sec;
  }

  POST0( isValidYYYYMMDDHHMMSS( result ) );
  return result;
}



/******************************************************************************
PURPOSE: toYYYYMMDDHHMMSS - Convert yyyydddhhmm to yyyymmddhhmmss.
INPUTS:  Integer yyyydddhhmm  Day of year format.
RETURNS: Integer yyyymmddhhmmss Year month day format.
******************************************************************************/

Integer toYYYYMMDDHHMMSS( Integer yyyydddhhmm ) {
  PRE0( isValidTimestamp( yyyydddhhmm ) );
  const Integer yyyyddd = yyyydddhhmm / 10000;
  const Integer hhmm    = yyyydddhhmm % 10000;
  Integer yyyy = yyyyddd / 1000;
  Integer mo = 0;
  Integer dd = 0;
  Integer result = 0;
  monthAndDay( yyyyddd, &mo, &dd );
  result =
    yyyy * 10000000000LL + mo * 100000000LL + dd * 1000000LL + hhmm * 100LL;
  POST0( isValidYYYYMMDDHHMMSS( result ) );
  return result;
}



/******************************************************************************
PURPOSE: hoursInRange - Number of hours from first to last timestamps.
INPUTS:  const UTCTimestamp first  First timestamp.
         const UTCTimestamp last   Last  timestamp.
RETURNS: Integer number of hours in [first, last].
******************************************************************************/

Integer hoursInRange( const UTCTimestamp first, const UTCTimestamp last ) {
  PRE03( isValidUTCTimestamp( first ),
         isValidUTCTimestamp( last ),
         fromUTCTimestamp( first ) <= fromUTCTimestamp( last ) );
  const Integer lastTimestamp = fromUTCTimestamp( last ) / 100 * 100;
  Integer timestamp = fromUTCTimestamp( first ) / 100 * 100;
  Integer result = 1;

  while ( timestamp < lastTimestamp ) {
    incrementTimestamp( &timestamp );
    ++result;
  }

  POST0( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: daysInRange - Number of days from first to last timestamps.
INPUTS:  const UTCTimestamp first  First timestamp.
         const UTCTimestamp last   Last  timestamp.
RETURNS: Integer number of days in [first, last].
******************************************************************************/

Integer daysInRange( const UTCTimestamp first, const UTCTimestamp last ) {
  PRE03( isValidUTCTimestamp( first ),
         isValidUTCTimestamp( last ),
         fromUTCTimestamp( first ) <= fromUTCTimestamp( last ) );
  const Integer lastTimestamp = fromUTCTimestamp( last ) / 10000 * 10000;
  Integer timestamp = fromUTCTimestamp( first ) / 10000 * 10000;
  Integer result = 1;

  while ( timestamp < lastTimestamp ) {
    int hours = 24;

    while ( hours-- ) {
      incrementTimestamp( &timestamp );
    }

    ++result;
  }

  POST0( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: isValidYYYYMMDDHHMMSS - Is the timestamp valid?
INPUTS:  Integer yyyymmddhhmmss.
******************************************************************************/

Integer isValidYYYYMMDDHHMMSS( Integer yyyymmddhhmmss ) {
  const Integer yyyy = yyyymmddhhmmss / 10000000000LL;
  const Integer mo = yyyymmddhhmmss / 100000000 % 100;
  const Integer dd = yyyymmddhhmmss / 1000000 % 100;
  const Integer hh = yyyymmddhhmmss / 10000 % 100;
  const Integer mm = yyyymmddhhmmss / 100 % 100;
  const Integer ss = yyyymmddhhmmss % 100;
  const Integer result =
    AND6( IN_RANGE( yyyy, 1900, 9999 ),
          IN_RANGE( mo, 1, 12 ),
          IN_RANGE( dd, 1, daysPerMonth[ isLeapYear( yyyy ) ][ mo - 1 ] ),
          IN_RANGE( hh, 0, 23 ),
          IN_RANGE( mm, 0, 59 ),
          IN_RANGE( ss, 0, 59 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: parseTimeRange - Parse string into its two integer timestamps.
INPUTS:  const char* string1         String representation of 1st timestamp.
         const char* string2         String representation of 2nd timestamp.
OUTPUTS: Integer*    yyyymmddhhmmss1 1st value in string.
         Integer*    yyyymmddhhmmss2 2nd value in string.
RETURNS: Integer 1 if successful, else 0.
******************************************************************************/

Integer parseTimeRange( const char* string1, const char* string2,
                        Integer* yyyymmddhhmmss1, Integer* yyyymmddhhmmss2 ) {

  PRE04( string1, string2, yyyymmddhhmmss1, yyyymmddhhmmss2 );

  const Integer result =
    AND5( sscanf( string1, "%lld", yyyymmddhhmmss1 ) == 1,
          sscanf( string2, "%lld", yyyymmddhhmmss2 ) == 1,
          isValidYYYYMMDDHHMMSS( *yyyymmddhhmmss1 ),
          isValidYYYYMMDDHHMMSS( *yyyymmddhhmmss2 ),
          *yyyymmddhhmmss1 <= *yyyymmddhhmmss2 );

  if ( ! result ) {
    *yyyymmddhhmmss1 = *yyyymmddhhmmss2 = 0;
  }

  POST02( IS_BOOL( result ),
         IMPLIES_ELSE( result,
                       AND3( isValidYYYYMMDDHHMMSS( *yyyymmddhhmmss1 ),
                             isValidYYYYMMDDHHMMSS( *yyyymmddhhmmss2 ),
                             *yyyymmddhhmmss1 <= *yyyymmddhhmmss2 ),
                       IS_ZERO2( *yyyymmddhhmmss1, *yyyymmddhhmmss2 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: totalSeconds - Compute total seconds from January 1 of year of
         first timestamp to first and second timestamps.
INPUTS:  Integer yyyymmddhhmmss1  1st timestamp.
         Integer yyyymmddhhmmss2  2nd timestamp.
OUTPUTS: Integer* seconds1  Seconds from year of 1st timestamp to 1st timestamp.
         Integer* seconds2  Seconds from year of 1st timestamp to 2nd timestamp.
RETURNS: Integer 1 if successful, else 0.
******************************************************************************/

void totalSeconds( Integer yyyymmddhhmmss1, Integer yyyymmddhhmmss2,
                   Integer* seconds1, Integer* seconds2 ) {

  PRE05( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss1 ),
         isValidYYYYMMDDHHMMSS( yyyymmddhhmmss2 ),
         yyyymmddhhmmss1 <= yyyymmddhhmmss2,
         seconds1,
         seconds2 );

  const Integer yyyy1 = yyyymmddhhmmss1 / 10000000000LL;
  Integer mo1         = yyyymmddhhmmss1 / 100000000 % 100;
  Integer dd1         = yyyymmddhhmmss1 / 1000000 % 100;
  const Integer hh1   = yyyymmddhhmmss1 / 10000 % 100;
  const Integer mm1   = yyyymmddhhmmss1 / 100 % 100;
  const Integer ss1   = yyyymmddhhmmss1 % 100;

  Integer yyyy2       = yyyymmddhhmmss2 / 10000000000LL;
  Integer mo2         = yyyymmddhhmmss2 / 100000000 % 100;
  Integer dd2         = yyyymmddhhmmss2 / 1000000 % 100;
  const Integer hh2   = yyyymmddhhmmss2 / 10000 % 100;
  const Integer mm2   = yyyymmddhhmmss2 / 100 % 100;
  const Integer ss2   = yyyymmddhhmmss2 % 100;

  /* Compute seconds1 = total seconds since yyyy of yyyymmddhhmmss1: */

  *seconds1 =
    ss1 + mm1 * SECONDS_PER_MINUTE +
    hh1 * MINUTES_PER_HOUR * SECONDS_PER_MINUTE;

  while ( --dd1 ) {
    *seconds1 += HOURS_PER_DAY * MINUTES_PER_HOUR * SECONDS_PER_MINUTE;
  }

  while ( --mo1 ) {
    *seconds1 +=
      daysInMonth( yyyy1, mo1 ) *
      HOURS_PER_DAY * MINUTES_PER_HOUR * SECONDS_PER_MINUTE;
  }

  /* Compute seconds2 = total seconds since yyyy of yyyymmddhhmmss2: */

  *seconds2 =
    ss2 + mm2 * SECONDS_PER_MINUTE +
    hh2 * MINUTES_PER_HOUR * SECONDS_PER_MINUTE;

  while ( --dd2 ) {
    *seconds2 += HOURS_PER_DAY * MINUTES_PER_HOUR * SECONDS_PER_MINUTE;
  }

  while ( --mo2 ) {
    *seconds2 +=
      daysInMonth( yyyy2, mo2 ) *
      HOURS_PER_DAY * MINUTES_PER_HOUR * SECONDS_PER_MINUTE;
  }

  while ( yyyy2 > yyyy1 ) {
    const Integer leap2 = isLeapYear( --yyyy2 );
    *seconds2 +=
      ( 365 + leap2 ) * HOURS_PER_DAY * MINUTES_PER_HOUR * SECONDS_PER_MINUTE;
  }

  POST02( *seconds1 >= 0, *seconds2 >= *seconds1 );
}



/******************************************************************************
PURPOSE: timestampOfTargetSeconds - Compute timestamp (yyyymmddhhmmss) that is
         yyyymddhhmmss (seconds + targetSeconds).
INPUTS:  Integer yyyymmddhhmmss  Starting timestamp.
         Integer seconds         Total seconds from Jan 1 of year of
                                 yyyymmddhhmmss.
         Integer targetSeconds   Total seconds from Jan 1 of year of
                                 yyyymmddhhmmss of target timestamp.
RETURNS: Integer yyyymmddhhmmss of targetSeconds.
******************************************************************************/

Integer timestampOfTargetSeconds( Integer yyyymmddhhmmss,
                                  Integer seconds, Integer targetSeconds ) {

  PRE03( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss ), seconds >= 0,
         targetSeconds >= seconds );

  Integer result = yyyymmddhhmmss;
  Integer secondsDifference = targetSeconds - seconds;

  if ( secondsDifference > 0 ) {
    Integer yyyy = yyyymmddhhmmss / 10000000000LL;
    Integer mo   = yyyymmddhhmmss / 100000000 % 100;
    Integer dd   = yyyymmddhhmmss / 1000000 % 100;
    Integer hh   = yyyymmddhhmmss / 10000 % 100;
    Integer mm   = yyyymmddhhmmss / 100 % 100;
    Integer ss   = yyyymmddhhmmss % 100;

    while ( secondsDifference-- ) {
      ++ss;

      if ( ss == SECONDS_PER_MINUTE ) {
        ss = 0;
        ++mm;

        if ( mm == MINUTES_PER_HOUR ) {
          mm = 0;
          ++hh;

          if ( hh == HOURS_PER_DAY ) {
            hh = 0;
            ++dd;
            CHECK( IN_RANGE( mo, 1, MONTHS_PER_YEAR ) );

            if ( dd > daysInMonth( yyyy, mo ) ) {
              dd = 1;
              ++mo;

              if ( mo > MONTHS_PER_YEAR ) {
                mo = 1;
                ++yyyy;
              }
            }
          }
        }
      }
    }

    CHECK6( yyyy >= 1900,
            IN_RANGE( mo, 1, 12 ),
            IN_RANGE( dd, 1, daysInMonth( yyyy, mo ) ),
            IN_RANGE( hh, 0, 23 ),
            IN_RANGE( mm, 0, 59 ),
            IN_RANGE( ss, 0, 59 ) );

    result =
      yyyy * 10000000000LL + mo * 100000000 + dd * 1000000 +
      hh * 10000 + mm * 100 + ss;
  }

  POST02( isValidYYYYMMDDHHMMSS( result ), result >= yyyymmddhhmmss );
  return result;
}



/******************************************************************************
PURPOSE: toUTCTimestamp2 - Convert timestamp to ISO UTC string format.
INPUTS:  Integer yyyymmddhhmmss  Timestamp to convert.
OUTPUTS: UTCTimestamp string  ISO UTC string format for the timestamp.
******************************************************************************/

void toUTCTimestamp2( Integer yyyymmddhhmmss, UTCTimestamp string ) {
  PRE02( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss ), string );
  const Integer yyyy = yyyymmddhhmmss / 10000000000LL;
  const Integer mo   = yyyymmddhhmmss / 100000000 % 100;
  const Integer dd   = yyyymmddhhmmss / 1000000 % 100;
  const Integer hh   = yyyymmddhhmmss / 10000 % 100;
  const Integer mm   = yyyymmddhhmmss / 100 % 100;
  const Integer ss   = yyyymmddhhmmss % 100;

  sprintf( string, "%04"INTEGER_FORMAT"-%02"INTEGER_FORMAT
                   "-%02"INTEGER_FORMAT
                   "T%02"INTEGER_FORMAT
                   ":%02"INTEGER_FORMAT":%02"INTEGER_FORMAT"-0000",
           yyyy, mo, dd, hh, mm, ss );
  POST0( strlen( string ) == UTC_TIMESTAMP_LENGTH );
}



/******************************************************************************
PURPOSE: previousDay - Timestamp of previous day.
INPUTS:  Integer yyyymmddhhmmss  Current day.
RETURNS: Integer timestamp of previous day.
******************************************************************************/

Integer previousDay( Integer yyyymmddhhmmss ) {
  PRE0( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss ) );
  const Integer hhmmss = yyyymmddhhmmss % 1000000;
  Integer yyyy = yyyymmddhhmmss / 10000000000LL;
  Integer mo = yyyymmddhhmmss / 100000000 % 100;
  Integer dd = yyyymmddhhmmss / 1000000 % 100;
  Integer result = 0;

  --dd;
  
  if ( dd < 1 ) {
    --mo;
    
    if ( mo < 1 ) {
      --yyyy;
      mo = 12;
    }
    
    dd = daysInMonth( yyyy, mo );
  }

  result = yyyy * 10000000000LL + mo * 100000000 + dd * 1000000 + hhmmss;
  POST02( isValidYYYYMMDDHHMMSS( result ), result < yyyymmddhhmmss );
  return result;
}



/******************************************************************************
PURPOSE: incrementHour - Increment timestamp by 1 hour.
INPUTS:  int* yyyymmddhh Current timestamp.
OUTPUTS: int* yyyymmddhh Incremented timestamp.
******************************************************************************/

void incrementHour( int* yyyymmddhh ) {
  PRE03( yyyymmddhh, isValidYearMonthDay( *yyyymmddhh / 100 ),
         IN_RANGE( *yyyymmddhh % 100, 0, 23 ) );
  int yyyy = *yyyymmddhh / 1000000;
  int mm   = *yyyymmddhh / 10000 % 100;
  int dd   = *yyyymmddhh / 100 % 100;
  int hh   = *yyyymmddhh % 100;

  ++hh;

  if ( hh > 23 ) {
    hh = 0;
    ++dd;

    if ( dd > daysPerMonth[ isLeapYear( yyyy ) ][ mm - 1 ] ) {
      dd = 1;
      ++mm;

      if ( mm > 12 ) {
        mm = 1;
        ++yyyy;
      }
    }
  }

  *yyyymmddhh = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;

  POST02( isValidYearMonthDay( *yyyymmddhh / 100 ),
          IN_RANGE( *yyyymmddhh % 100, 0, 23 ) );
}


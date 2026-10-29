
#ifndef DATETIME_H
#define DATETIME_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: DateTime.h - Declare routines for date/time computation.

NOTES:   Commands:
           void incrementOneHour( Integer* yyyyddd, Integer* hhmmss );
           void decrementOneHour( Integer* yyyyddd, Integer* hhmmss );
         Queries:
           Integer isValidDate( Integer yyyyddd );
           Integer isValidTime( Integer hhmmss );
           Integer isLeapYear( Integer yyyy );
           Integer timestepsInRange( Integer firstDate, Integer firstTime,
                                     Integer lastDate,  Integer lastTime );
           void monthAndDay( Integer yyyyddd,
                             Integer* month, Integer* dayOfMonth );

HISTORY: 2004/10 Todd Plessel EPA/LM Created.

STATUS:  unreviewed, untested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <BasicNumerics.h> /* For Integer. */

/*================================ FUNCTIONS ================================*/

enum {
  MONTHS_PER_YEAR = 12, SECONDS_PER_MINUTE = 60, MINUTES_PER_HOUR = 60,
  HOURS_PER_DAY = 24
};

extern Integer isValidDate( Integer yyyyddd );

extern Integer isValidTime( Integer hhmmss );

extern Integer isValidTimestepSize( Integer hhmmss );

extern Integer isLeapYear( Integer yyyy );

extern Integer daysInMonth( Integer year, Integer month );

extern Integer timestepsInRange( Integer firstDate, Integer firstTime,
                                 Integer lastDate,  Integer lastTime );

extern void monthAndDay( Integer yyyyddd, Integer* month, Integer* dayOfMonth);

extern void incrementOneHour( Integer* yyyyddd, Integer* hhmmss );

extern void decrementOneHour( Integer* yyyyddd, Integer* hhmmss );

extern void incrementTime( Integer* yyyyddd, Integer* hhmmss, Integer step );

extern void decrementTime( Integer* yyyyddd, Integer* hhmmss, Integer step );



/* New date/time routines: */

enum { UTC_TIMESTAMP_LENGTH = 24 }; /* YYYY-MM-DDTHH:MM:SS-ZZZZ */
typedef char UTCTimestamp[ UTC_TIMESTAMP_LENGTH + 1 ];

extern Integer isValidUTCTimestamp( const char* string );
extern void toUTCTimestamp( Integer yyyydddhhmm, UTCTimestamp string );
extern Integer fromUTCTimestamp( const UTCTimestamp string );
extern void incrementTimestamp( Integer* yyyydddhhmm );
extern void decrementTimestamp( Integer* yyyydddhhmm );
extern Integer offsetTimestamp( Integer yyyydddhhmm, Integer hours );
extern Integer parseTimestamp( const char* string, Integer* yyyydddhh00 );
extern Integer isValidTimestamp( Integer yyyydddhhmm );
extern Integer isValidYearMonthDay( Integer yyyymmdd );
extern Integer convertYearMonthDay( Integer yyyymmdd );
extern Integer nowUTC( void ); /* Returns YYYYDDDHHMM. */
extern Integer fromSeconds( const Integer seconds ); /*Returns YYYYMMDDHHMMSS*/

extern Integer hoursInRange(const UTCTimestamp first, const UTCTimestamp last);
  
extern Integer daysInRange(const UTCTimestamp first, const UTCTimestamp last);

extern Integer toYYYYMMDDHHMMSS( Integer yyyydddhhmm );

extern Integer isValidYYYYMMDDHHMMSS( Integer yyyymmddhhmmss );

extern Integer parseTimeRange( const char* string1, const char* string2,
                               Integer* yyyymmddhhmmss1,
                               Integer* yyyymmddhhmmss2 );

extern void totalSeconds( Integer yyyymmddhhmmss1, Integer yyyymmddhhmmss2,
                          Integer* seconds1, Integer* seconds2 );

extern Integer timestampOfTargetSeconds( Integer yyyymmddhhmmss,
                                         Integer seconds,
                                         Integer targetSeconds );

extern void toUTCTimestamp2( Integer yyyymmddhhmmss, UTCTimestamp string );

extern Integer previousDay( Integer yyyymmddhhmmss );

extern void incrementHour( int* yyyymmddhh );

#ifdef __cplusplus
}
#endif

#endif /* DATETIME_H */



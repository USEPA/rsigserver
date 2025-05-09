
/******************************************************************************
PURPOSE: Failure.c - Defines routines for setting user-defined failure
         handlers.
NOTES:   The presence of the 'userFailureHandler' static 'global' variable
         means that multi-thread applications cannot have thread-specific
         handlers.
HISTORY: 1993/03 Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>        /* For fprintf() and vfprintf().                   */
#include <string.h>       /* For strcat().                                   */
#include <errno.h>        /* For errno and strerror().                       */
#include <limits.h>       /* For UINT_MAX.                                   */
#include <stdarg.h>       /* For va_args.                                    */

#include <Assertions.h>    /* For PRE0*(), POST0*(), DEBUG().                */
#include <BasicNumerics.h> /* For INTEGER_MAX, INTEGER_FORMAT, Integer.      */
#include <Failure.h>       /* For the routine definitions.                   */

#if __sgi||__sun||__osf__||_CRAY||_AIX||__OPENNT||__linux__||__APPLE__||\
_WIN32||__CYGWIN__
#define USE_MY_VNSPRINTF
#define vnsprintf my_vnsprintf
static int my_vnsprintf( char* buffer, size_t buffer_size, const char* format,
                         va_list args );
#endif

#if defined( __linux__ ) || defined( __APPLE__ ) || defined( _WIN32 ) || \
    defined( __CYGWIN__ )
#define STDERR_NON_STATIC 1
#else
#define STDERR_NON_STATIC 0
#endif

#if STDERR_NON_STATIC
#define ASSIGN_FILE(unused) 0
#define ENSURE_ASSIGNED_OR_STDERR(f) ((f) == 0 ? (f) = stderr : (f))
#define ENSURE_ASSIGNED_OR_STDOUT(f) ((f) == 0 ? (f) = stdout : (f))
#else
#define ASSIGN_FILE(f) f
#define ENSURE_ASSIGNED_OR_STDERR(f) f
#define ENSURE_ASSIGNED_OR_STDOUT(f) f
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

/* No-op macro used for readability in PRE0/POST0 of failureLogFile(). */
#define IS_OPEN_FOR_WRITING(f) 1

/*============================ PRIVATE VARIABLES ============================*/

static FailureHandler userFailureHandler; /* Points to user's handler.     */
static Integer totalFailureCount;         /* Count of all failures so far. */
static Integer logging = 1;               /* Failure printing turned on?   */
static Integer calling = 1;               /* Failure callback turned on?   */
static Integer ringing = 1;               /* Terminal bell ringing on?     */
static Integer verbose = 1;               /* Verbose message adornment on? */
static FILE* theFailureLogFile = ASSIGN_FILE(stderr); /* Default messages. */
static FILE* theInfoLogFile    = ASSIGN_FILE(stdout); /* Default messages. */
static const char* programName = 0;       /* Optional name of program.     */

/*============================ PUBLIC FUNCTIONS =============================*/


/******************************************************************************
PURPOSE: failureCount - Get the total number of failures that have occurred
         since the process started. Initally zero.
RETURNS: Integer The current failure count.
******************************************************************************/

Integer failureCount( void ) {
  const Integer result = totalFailureCount;
  POST0( result >= 0 );
  return result;
}


/******************************************************************************
PURPOSE: failureProgramName - The name of the program used to prefix messages.
         If not set by failureSetProgramName() then
         it defaults to 0 and is not used to prefix messages.
RETURNS: const char* The name of the program (e.g., argv[0]) or 0 if not set.
******************************************************************************/

const char* failureProgramName( void ) {
  const char* const result = programName;
  return result;
}


/******************************************************************************
PURPOSE: failureLoggingEnabled - Determine if failure logging (printing to
         log files) is enabled or not. On by default.
RETURNS: Integer 1 if enabled else 0.
******************************************************************************/

Integer failureLoggingEnabled( void ) {
  return logging;
}


/******************************************************************************
PURPOSE: failureCallingEnabled - Determine if failure calling (of optional
         client handler) is enabled or not. On by default.
RETURNS: Integer 1 if enabled else 0.
******************************************************************************/

Integer failureCallingEnabled( void ) {
  return calling;
}


/******************************************************************************
PURPOSE: failureRingingEnabled - Determine if ringing of the terminal bell
         upon failureMessage() calls is enabled or not. On by default.
RETURNS: Integer 1 if enabled else 0.
******************************************************************************/

Integer failureRingingEnabled( void ) {
  return ringing;
}


/******************************************************************************
PURPOSE: failureVerboseEnabled - Determine if verbose message adornments
         are enabled or not. On by default.
RETURNS: Integer 1 if enabled else 0.
******************************************************************************/

Integer failureVerboseEnabled( void ) {
  return verbose;
}


/******************************************************************************
PURPOSE: failureHandler - Returns the current failure handler (if any).
RETURNS: FailureHandler  The current failure handler or 0 if none exists.
NOTES:   Handler is an optional client-supplied routine that will be called
         (if it exists) when failures are encountered and reported via
         failureMessage() AND failureCallingEnabled() is true (the default).
******************************************************************************/

FailureHandler failureHandler( void ) {
  return userFailureHandler;
}


/******************************************************************************
PURPOSE: failureLogFile - The opened file used to print failure messages to.
         (Default is stderr.)
RETURNS: FILE* The file used to print messages to.
******************************************************************************/

FILE* failureLogFile( void ) {
  FILE* const result = ENSURE_ASSIGNED_OR_STDERR( theFailureLogFile );
  POST02( result, IS_OPEN_FOR_WRITING( result ) );
  return result;
}


/******************************************************************************
PURPOSE: infoLogFile - The opened file used to print info messages to.
         (Default is stdout.)
RETURNS: FILE* The file used to print messages to.
******************************************************************************/

FILE* infoLogFile( void ) {
  FILE* const result = ENSURE_ASSIGNED_OR_STDOUT( theInfoLogFile );
  POST02( result, IS_OPEN_FOR_WRITING( result ) );
  return result;
}


/******************************************************************************
PURPOSE: failureSetHandler - Sets a new user-specified failure handler.
INPUTS:  FailureHandler newFailureHandler  Failure handler routine or 0.
NOTES:   Handler is an optional client-supplied routine that will be called
         (if it exists) when failures are encountered and reported via
         failureMessage() AND failureCallingEnabled() is true (the default).
******************************************************************************/

void failureSetHandler( FailureHandler newFailureHandler ) {
  userFailureHandler = newFailureHandler;
  POST0( failureHandler() == newFailureHandler );
}


/******************************************************************************
PURPOSE: failureSetProgramName - Set the name of the program used to prefix
         messages. Default is 0. When 0 no prefixing is done.
RETURNS: const char* The name of the program (e.g., argv[0]) or 0 if not set.
******************************************************************************/

void failureSetProgramName( const char* newProgramName ) {
  PRE02( newProgramName, strlen( newProgramName ) < 256 );
  programName = newProgramName;
  POST0( programName == newProgramName );
}


/******************************************************************************
PURPOSE: failureEnableLogging - Enable failure logging (printing of messages
         to log file).
NOTE:    Subsequent failures will be printed to log file.
******************************************************************************/

void failureEnableLogging( void ) {
  logging = 1;
  POST0( failureLoggingEnabled() );
}


/******************************************************************************
PURPOSE: failureDisableLogging - Disable failure reporting.
NOTE:    Failures will not be reported, that is, printed to the log file,
         but will still be counted.
******************************************************************************/

void failureDisableLogging( void ) {
  logging = 0;
  POST0( ! failureLoggingEnabled() );
}


/******************************************************************************
PURPOSE: failureEnableCalling - Enable callbacks to client handler.
NOTE:    Client handler will be called when subsequent failures arise.
******************************************************************************/

void failureEnableCalling( void ) {
  calling = 1;
  POST0( failureCallingEnabled() );
}


/******************************************************************************
PURPOSE: failureDisableCalling - Disable callbacks to client handler.
NOTE:    Client handler will not be called when subsequent failures arise,
         but will still be counted.
******************************************************************************/

void failureDisableCalling( void ) {
  calling = 0;
  POST0( ! failureCallingEnabled() );
}


/******************************************************************************
PURPOSE: failureEnableRinging - Enable ringing of the terminal bell when
         failureMessage() is called. On by default.
NOTE:    No effect unless failureLoggingEnabled() != 0.
******************************************************************************/

void failureEnableRinging( void ) {
  ringing = 1;
  POST0( failureRingingEnabled() );
}


/******************************************************************************
PURPOSE: failureDisableRinging - Disable ringing of the terminal bell when
         failureMessage() is called. On by default.
******************************************************************************/

void failureDisableRinging( void ) {
  ringing = 0;
  POST0( ! failureRingingEnabled() );
}


/******************************************************************************
PURPOSE: failureEnableVerbose - Enable construction of verbose/adorned messages
         when failureMessage() is called. On by default.
******************************************************************************/

void failureEnableVerbose( void ) {
  verbose = 1;
  POST0( failureVerboseEnabled() );
}


/******************************************************************************
PURPOSE: failureDisableVerbose - Disable construction of verbose/adorned
         messages when failureMessage() is called. On by default.
******************************************************************************/

void failureDisableVerbose( void ) {
  verbose = 0;
  POST0( ! failureVerboseEnabled() );
}


/******************************************************************************
PURPOSE: failureSetLogFile - Set the opened file used to print failure messages
         to. (Default is stderr.)
INPUTS:  FILE* The file opened for writing.
******************************************************************************/

void failureSetLogFile( FILE* newLogFile ) {
  PRE02( newLogFile, IS_OPEN_FOR_WRITING( newLogFile ) );
  theFailureLogFile = newLogFile;
  POST0( failureLogFile() == newLogFile );
}


/******************************************************************************
PURPOSE: infoSetLogFile - Set the opened file used to print info messages to.
         (Default is stdout.)
INPUTS:  FILE* The file opened for writing.
******************************************************************************/

void infoSetLogFile( FILE* newLogFile ) {
  PRE02( newLogFile, IS_OPEN_FOR_WRITING( newLogFile ) );
  theInfoLogFile = newLogFile;
  POST0( infoLogFile() == newLogFile );
}


/******************************************************************************
PURPOSE: failureMessage - Prints the given failure message (plus strerror())
         to the failureLogFile and then invokes the user failure handler
         (if any).
INPUTS:  const char* message  A string like those used in printf().
         ...                  Other arguments implied by '%' in message.
NOTES:   Implemented in terms of vsprintf (man vsprintf).
         For Integers use INTEGER_FORMAT macro from BasicNumerics.h:
         const Integer i = INTEGER_CONSTANT(12345678901234567);
         failureMessage( "example: i = %30"INTEGER_FORMAT" bytes.\n", i );
         Also, see BasicNumerics.h for REAL_G_FORMAT, etc.
******************************************************************************/

void failureMessage( const char* message, ... ) {
  /* REMEMBER_F( Integer, failureCount ); */

  /*
   * errnoCopy holds current value to pass to user after clearing the real
   * errno to protect against user longjumps.
   */

  const Integer errnoCopy = errno;

  enum { MAX_MESSAGE_SIZE = 1024 };
  char expandedMessage[ MAX_MESSAGE_SIZE + 1 ]; /* Buffer for full message. */

  const char* const prefix = "I'm sorry: ";
  const char* const reason = "\nPossible reason: ";
  const char* const default_explanation = "Invalid user or data input.";
  const char* const system_explanation =
    "Temporary system resource acquisition/access/usage failure:\n";
  const char* const errno_explanation = strerror( errno );
  const char* const advice =
    "\nSee console window for possible details "
    "then perhaps try operation again.";

  const size_t length_of_system_explanation  = strlen( system_explanation );
  const size_t length_of_errno_explanation   = strlen( errno_explanation );
  const size_t length_of_default_explanation = strlen( default_explanation );

  const size_t length_of_explanation =
    MAX( length_of_system_explanation + length_of_errno_explanation,
         length_of_default_explanation );

  const size_t length_of_additional_messages =
    strlen( reason ) + length_of_explanation + strlen( advice );


  /* Increment global failure count, skipping 0 when wrapping around: */

  totalFailureCount = NEXT_STRICTLY_POSITIVE_INTEGER( totalFailureCount );

  /* Begin message with a standard prefix (with or without programName): */

  expandedMessage[ 0 ] = '\0';

  if ( NON_ZERO2( programName, *programName ) ) {
    const size_t programNameLength = strlen( programName );
    const size_t maximumProgramNameLength =
      MAX( programNameLength, MAX_MESSAGE_SIZE / 4 );

    CHECK( strlen( prefix ) + 2 < MAX_MESSAGE_SIZE / 4 );
    strncpy( expandedMessage, programName, maximumProgramNameLength );
    strcat( expandedMessage, ": " );
  }

  if ( verbose ) {
    strcat( expandedMessage, prefix ); /* Add standard prefix. */
  }

  if ( AND2( message, *message ) ) {
    const size_t maximum_length = MAX_MESSAGE_SIZE -
      strlen( expandedMessage ) - length_of_additional_messages;
    va_list args;              /* For stdarg magic. */

    va_start( args, message ); /* Begin stdarg magic. */
    CHECK( MAX_MESSAGE_SIZE >
           strlen( expandedMessage ) + length_of_additional_messages );
    vnsprintf( expandedMessage + strlen( expandedMessage ),
               maximum_length, message, args );
    va_end( args );            /* End of stdarg magic. */
  }

  if ( verbose ) {

    /* Add the reason and advice: */

    strcat( expandedMessage, reason );

    if ( errno ) {
      strcat( expandedMessage, system_explanation );
      strcat( expandedMessage, errno_explanation );
    } else {
      strcat( expandedMessage, default_explanation );
    }

    strcat( expandedMessage, advice );
  }

  CHECK( strlen( expandedMessage ) < MAX_MESSAGE_SIZE ); /* HACK too late... */

  if ( failureLoggingEnabled() ) {

    /* Print newlines, ring bell and print the expanded message. */

    fprintf( ENSURE_ASSIGNED_OR_STDERR( theFailureLogFile ), "\n\n" );
    fprintf( theFailureLogFile, "%s", expandedMessage );

    /* Print the current failure count. */

    fprintf( theFailureLogFile, "\n(program failure # %"INTEGER_FORMAT")\n\n",
             totalFailureCount );

    if ( ringing ) {
      fprintf( stderr, "\a\n" );
    }
  }

  errno = 0; /* Clear global now in case userFailureHandler never returns. */

  if ( failureCallingEnabled() ) {

    /* Finally, call the user's failure handler routine (if it exists). */

    if ( userFailureHandler ) {
      userFailureHandler( errnoCopy, expandedMessage );
    }
  }

  /*
  POST0( failureCount() == NEXT_STRICTLY_POSITIVE_INTEGER( OLD(failureCount)));
  */
}


/******************************************************************************
PURPOSE: infoMessage - Prints the given info message to the failureLogFile and
         then invokes the user failure handler (if any).
INPUTS:  const char* message  A string like those used in printf().
         ...                  Other arguments implied by '%' in message.
NOTES:   Implemented in terms of vsprintf (man vsprintf).
         For Integers use INTEGER_FORMAT macro from BasicNumerics.h:
         const Integer i = INTEGER_CONSTANT(12345678901234567);
         infoMessage( "example: i = %30"INTEGER_FORMAT" bytes.\n", i );
         Also, see BasicNumerics.h for REAL_G_FORMAT, etc.
******************************************************************************/

void infoMessage( const char* message, ... ) {

  enum { MAX_MESSAGE_SIZE = 1024 };
  char expandedMessage[ MAX_MESSAGE_SIZE + 1 ]; /* Buffer for full message. */
  const char* const prefix = "Info: ";
  REMEMBER_F( Integer, failureCount );

  /* Begin message with a standard prefix (with or without programName): */

  expandedMessage[ 0 ] = '\0';

  if ( NON_ZERO2( programName, *programName ) ) {
    const size_t programNameLength = strlen( programName );
    const size_t maximumProgramNameLength =
      MAX( programNameLength, MAX_MESSAGE_SIZE / 4 );

    CHECK( strlen( prefix ) + 2 < MAX_MESSAGE_SIZE / 4 );
    strncpy( expandedMessage, programName, maximumProgramNameLength );
    strcat( expandedMessage, ": " );
  }

  if ( verbose ) {
    strcat( expandedMessage, prefix );
  }

  if ( AND2( message, *message ) ) {

    const size_t maximum_length = MAX_MESSAGE_SIZE - strlen( expandedMessage );
    va_list args;              /* For stdarg magic. */

    va_start( args, message ); /* Begin stdarg magic. */
    CHECK( MAX_MESSAGE_SIZE > strlen( expandedMessage ) );
    vnsprintf( expandedMessage + strlen( expandedMessage ),
               maximum_length, message, args );
    va_end( args );            /* End of stdarg magic. */
  }

  CHECK( strlen( expandedMessage ) < MAX_MESSAGE_SIZE ); /* HACK too late... */

  if ( failureLoggingEnabled() ) {

     /* Print the expanded message. */

    fprintf( ENSURE_ASSIGNED_OR_STDOUT(theInfoLogFile), "%s", expandedMessage );
  }

  if ( failureCallingEnabled() ) {

    /* Finally, call the user's failure handler routine (if it exists). */

    if ( userFailureHandler ) {
      userFailureHandler( 0, expandedMessage );
    }
  }

  POST0( failureCount() == OLD( failureCount ) );
}


/*============================ PRIVATE FUNCTIONS ============================*/


#ifdef USE_MY_VNSPRINTF

/******************************************************************************
PURPOSE: my_vnsprintf - Implements vnsprintf() for platforms that lack it.
INPUTS:  size_t buffer_size   Size of output buffer.
         const char* format   A string like those used in printf().
         va_list args         Other arguments implied by '%' in format.
OUTPUTS: const char* buffer   String initialized by format and args, or "".
RETURNS: int Number of characters written to buffer - i.e., its length.
NOTES:   Implemented in terms of vfprintf() to /dev/null HACK.
******************************************************************************/

static int my_vnsprintf( char* buffer, size_t buffer_size, const char* format,
                         va_list args ) {
  PRE04( buffer, buffer_size, format, *format );
  int result = 0;
  FILE* test_sink = fopen( "/dev/null", "w" );
  *buffer = '\0';

  if ( test_sink ) {
#ifdef __x86_64__
    /*
     * __x86_64__ BUG: must copy args otherwise after passing it to vfprintf,
     * apparently its contents are altered, causing its second passing to
     * vsprintf below to crash!
     */
    int required_buffer_size = 0;
    va_list args_copy;
    va_copy( args_copy, args ); /* UGLY: this is a non-portable macro! */
    required_buffer_size = vfprintf( test_sink, format, args_copy );
#else
    const int required_buffer_size = vfprintf( test_sink, format, args );
#endif
    fclose( test_sink );
    test_sink = 0;

    if ( IN_RANGE( required_buffer_size, 1, buffer_size - 1 ) ) {
      result = vsprintf( buffer, format, args );
    }
  }

  POST0( result >= 0 );
  return result;
}

#endif /* ifdef USE_MY_VNSPRINTF */


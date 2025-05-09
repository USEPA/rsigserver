
#ifndef FAILURE_H
#define FAILURE_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Failure.h - Declares routines for global failure logging and handling.
NOTES:   Example usage:

           #include <stdio.h>

           #include <Memory.h>
           #include <Failure.h>

           static void yourFailureHandler( Integer failureCode,
                                           const char* message ) {
             printf( "Got failure code %"INTEGER_FORMAT".\n", failureCode );
             printf( "Got message '%s'\n", message );
           }

           int main( int argc, char* argv[] ) {
             FailureHandler oldHandler = failureHandler();
             FILE* nullFile = 0;
             Integer* p = 0;

             failureSetProgramName( argv[ 0 ] );
             failureSetHandler( yourFailureHandler );
             p = NEW( Integer, 1000000000 );
             FREE( p );
             nullFile = fopen( "/dev/null", "w" );

             if ( nullFile ) {
               failureSetLogFile( nullFile );
               failureMessage( "Something is amuck with p = %p.", p );
               fclose( nullFile );
               nullFile = 0;
             }

             failureSetHandler( oldHandler );
             return 0;
           }

         $ cc -n32 -mips3 -I../../../../include -L../../../../lib/IRIX \
              -o example example.c -lFailure.debug -lMemeory.debug

         $ example


         Correctness vs Robustness or Defects vs Failures:

         Errors are human mistakes. When the humans are programmers, these
         mistakes are often manifest as defects in the software they write -
         incorrect or missing code fragments. When these programs are executed
         their defects may cause faults - undesired behavior (e.g., hangs,
         crashes, or worse, unnoticed data corruption). The realm of
         errors/defects/faults is excellently, if partially, addressed by
         Design By Contract - specification assertions optionally checkable by
         macros. Use Assertions.h for this purpose. It will allow the bulk
         of (library) code not dealing with the 'external world' (e.g.,
         files, users) to express correctness conditions without the
         inefficient complex redundancy of 'defensive programming'.

         For the minority of code that does interface with the 'external world'
         of data/inputs/outputs/events, the issue is robustness - the
         ability to detect and appropriately handle anomalies that are not
         due to defects and are beyond the control of the given program(mer).
         This library is used to support one aspect of robustness, namely,
         reporting failures due to detected 'external factors'.

         Failures are of two kinds:

           I.  Due to resource acquisition/access/usage failures.
           II. Due to invalid (user/data) input.

         Examples of type I are: failure to successfully allocate a segment of
         memory, open a file or seek beyond the start of a file (defect?).

         Type I failure messages have the following structure:

         <program-name><generic-apology><description-of-failure>
         <generic-reason-for-failure>
         <specific-reason-for-failure>
         <generic-suggestion-of-how-to-proceed>
         <current-total-failure-count>

         When the generic parts are filled-in it appears as:

         [<programName> :]I'm sorry: <expanded-message>
         Reason: Temporary system resource acquisition/access/usage failure:
         <strerror>
         See console window for possible details then perhaps try operation
         again.
         (program failure # <N>)

         <program-name> is established by calling failureSetProgramName(argv[0])
         otherwise it defaults to 0 and is not printed. Users prefer polite
         software thus the apology. errno is used to elaborate on the
         underlying (low-level) reason for the failure.

         Example:

         FailureTest: I'm sorry: Couldn't open file '/a/non-existent/file'.
         Reason: Temporary system resource acquisition/access/usage failure:
         No such file or directory
         See console window for possible details then perhaps try operation
         again.
         (program failure # 20)

         In type II failures, errno is irrelevant since the cause of the
         failure is higher-level - typically mis-typed input data values or
         inappropriately-formatted input data files.

         An example of a type II failure message is:

         Importer: I'm sorry: Can't read 88000 integers from file 'bogus'.
         Reason: Invalid user or data input.
         See console window for possible details then perhaps try operation
         again.
         (program failure # 13)


         Basic usage:

         To log failures call failureMessage() as if it were printf().
         It accepts the same format-string (with %-specifiers) and
         optional arguments (stdarg). failureMessage( "like this %d.", 1 );

         By default it will print to stderr.
         Use failureSetLogFile( yourOpenedFile ) to redirect subsequent
         messages elsewhere or call failureDisableLogging() to supress
         such logging altogether. Subsequent calls to failureMessage()
         will still be counted though. The current count is available
         from failureCount().

         Orthogonal to logging is callbacks. You can install a callback
         handler routine that will be called with each occurrence of the
         expanded message resulting from failureMessage() and also the
         'errno' code (man perror strerror) or 0 if type II failure.

           void yourFailureHandler( Integer failureCode,
                                    const char* message ) {
             ... do whatever ...
           }

           failureSetHandler( yourFailureHandler );

         Since 'do whatever' may include calls to abort() or exit() or
         long-jump to main() (not recommended for 'robust' applications)
         it is important that 'library' routines 'clean-up state' (e.g.,
         restore invariant) before calling failureMessage() since there
         may not be a chance to do so afterwards. Conversely, such
         library routines should not be written with the assumption that
         failureMessage() never returns. That is, robust routines that
         call failureMessage() must be designed to be correct regardless
         if failureMessage() returns or not.
         (Strousrtup's Acquisition-Is-Initialization idiom is appropriate
         to help ensure this.)

         Sometimes library routines call other library routines which may
         also call failureMessage(). In such cases, it is often better to
         temporarily suppress logging and callbacks until such caller
         routines have had a chance to handle/summarize failures itself
         so that client code is not called with multiple failures.
         This can be accomplished with:

           const Integer currentFailureCount = failureCount();
           const Integer wereLogging = failureLoggingEnabled();
           const Integer wereCalling = failureCallingEnabled();
           failureDisableLogging();
           failureDisableCalling();

           callRoutineThatMightFail();
           callAnotherRoutineThatMightFail();

           if ( wereLogging ) {
             failureEnableLogging();
           }

           if ( wereCalling ) {
             failureEnableCalling();
           }

           if ( failureCount() != currentFailureCount ) {
             failureMessage( "No can-do." );
           }

         Finally, informational messages can be displayed using
         infoMessage() (also like printf()) and similarly re-directed
         from stdout (the default) to another file using infoSetLogFile().
         Informational message have the form:

         [<program-name> :]Info :<expanded-message>

         An example of an info message is:

         Export: Info: file = 'wonderful_result', timestep = 128

HISTORY: 1993/03 Todd Plessel EPA/LM Created.
STATUS:  reviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h> /* For FILE. */

#include <BasicNumerics.h> /* For Integer. */

/*================================== TYPES ==================================*/

typedef void (*FailureHandler)( Integer failureCode, const char* message );
/* PRE3( failureCode >= 0, message, *message ); */

/*================================ FUNCTIONS ================================*/

/* Queries: */

extern Integer failureCount( void );
extern Integer failureLoggingEnabled( void );
extern Integer failureCallingEnabled( void );
extern Integer failureRingingEnabled( void );
extern Integer failureVerboseEnabled( void );
extern FailureHandler failureHandler( void );
extern FILE* failureLogFile( void );
extern const char* failureProgramName( void );
extern FILE* infoLogFile( void );

/* Commands: */

extern void failureSetHandler( FailureHandler newFailureHandler );
extern void failureEnableLogging( void );
extern void failureDisableLogging( void );
extern void failureEnableCalling( void );
extern void failureDisableCalling( void );
extern void failureEnableRinging( void );
extern void failureDisableRinging( void );
extern void failureEnableVerbose( void );
extern void failureDisableVerbose( void );
extern void failureSetLogFile( FILE* newLogFile );
extern void failureSetProgramName( const char* newProgramName );
extern void infoSetLogFile( FILE* newLogFile );
extern void failureMessage( const char* message, ... );
extern void infoMessage( const char* message, ... );

#ifdef __cplusplus
}
#endif

#endif /* FAILURE_H */


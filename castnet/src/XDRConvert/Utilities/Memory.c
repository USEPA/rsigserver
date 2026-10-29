
/******************************************************************************
PURPOSE: Memory.c - Defines memory allocation routines.
NOTES:   Calls the Standard C Library malloc(), realloc() and free() routines.
HISTORY: 1993/03, Todd Plessel, EPA/LM, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h>        /* For malloc(), free(), realloc(). */
#include <string.h>        /* For memset().                    */

#include <Assertions.h>    /* For PRE(), POST(), DEBUG().  */
#include <BasicNumerics.h> /* For Integer, isSizet().      */
#include <Failure.h>       /* For failureMessage().        */
#include <Memory.h>        /* For the routine definitions. */

/*============================ PRIVATE VARIABLES ============================*/

static Integer failureCountDown = 0; /* UGLY: GLOBAL variable in library! */

/*
malloc BUG: Must work-around non-robust malloc() which either hangs or
crashes rather than returns 0 when unable to allocate the requested bytes.
*/

static unsigned long long available_bytes( void );

static int available( unsigned long long bytes ) {
  const unsigned long long safe_minimum_free = 5 * 1024 * 1024; /* 5MB. */
  const unsigned long long free_bytes = available_bytes();
  const int result =
    free_bytes > bytes && free_bytes - bytes > safe_minimum_free;
  return result;
}

/*============================ PUBLIC FUNCTIONS =============================*/


/******************************************************************************
PURPOSE: freeMemory - Implements the FREE() macro by calling free().
         This provides a checking and debugging wrapper around free().
INPUTS:  void* address  The address to free.
******************************************************************************/

void freeMemory( void* address ) {
  PRE0( address );
  char* p = (char*) address;
  DEBUG( fprintf( stderr, "freeMemory( address = %p )\n", address ); )
  *p = '\0'; /* Zero first byte to reveal dangling pointer in client. */
  free( address );
}


/******************************************************************************
PURPOSE: newMemory - Implements NEW() macro by calling the standard malloc()
         and, upon failure failureMessage() is called.
         This provides a checking and debugging wrapper around malloc().
INPUTS:  Integer count     The number of items to allocate.
         Integer sizeEach  The number of bytes per item.
         Integer zeroIt    Zero-out the memory?
RETURNS: void*  The resulting address, or 0 if unsuccessful.
NOTES:   If malloc() returns 0 then failureMessage() is invoked.
******************************************************************************/

void* newMemory( Integer count, Integer sizeEach, Integer zeroIt ) {
  PRE03( count > 0, sizeEach > 0, IS_BOOL( zeroIt ) );
  void* address = 0;
  const Integer bytes = count * sizeEach;
  Integer forceFailure = 0;

  DEBUG( fprintf( stderr, "newMemory( count = %"INTEGER_FORMAT", "
                  "sizeEach = %"INTEGER_FORMAT", "
                  "zeroIt = %"INTEGER_FORMAT" ) "
                  "(failureCountDown = %"INTEGER_FORMAT" ) ",
                  count, sizeEach, zeroIt, failureCountDown ); )

  /* Decide if simulating failure or not: */

  if ( failureCountDown > 0 ) { /* If we're counting down then... */
    --failureCountDown;
    forceFailure = failureCountDown == 0; /* Fail iff 0 is reached. */
  }

  /* Don't attempt to allocate if too large to represent: */

  if ( AND8( ! forceFailure, isSizet( count ), isSizet( sizeEach ),
             bytes > 0, bytes >= count, bytes >= sizeEach, isSizet( bytes ),
             available( bytes ) ) ) {
    address = malloc( bytes );
    DEBUG( fprintf( stderr, " yields address = %p\n", address ); )
  }

  if ( address == 0 ) {
    failureMessage( "Can't allocate %"INTEGER_FORMAT" bytes to complete"
                    " the requested action.", bytes );
  } else if ( zeroIt ) {
    memset( address, 0, bytes );
  }

  return address;
}


/******************************************************************************
PURPOSE: resizeMemory - Implements the RESIZE() macro by calling the standard
         realloc() and, upon failure failureMessage() is called.
         This provides an easier-to-use substitute for the notoriously
         problematic realloc().
INPUTS:  Integer typeSize        Size of type in bytes.
         const void** existingAddress  Pointer to pointer to a block
                                 to be reallocated (or pointer to 0 to
                                 allocate a new block).
         Integer* existingCount  Pointer to # of items in existing block.
         Integer deltaCount      The extra (or fewer) number of items to
                                 grow (or shink) the existing block by.
                                 If deltaCount == -existingCount then the
                                 existing block is freed (and
                                 *existingAddress == 0).
         Integer zeroDelta       Zero-out extra memory (if deltaCount > 0)?
OUTPUTS: const void** existingAddress  Possibly changed address of a block of
                                       size existingCount + deltaCount items.
         Integer* existingCount  Pointer to new # of items in existing block.
RETURNS: Integer 1 if successful, else 0.
NOTES:   If realloc() fails then failureMessage() is invoked and
         existingAddress and existingCount are unchanged.

         Example:

           const Integer initialEmployeeCount = 100;
           Employee* employees = NEW_ZERO( Employee, initialEmployeeCount );

           if ( employees ) {

             Integer employeeCount = initialEmployeeCount;

             ... do stuff with employees ...

             if ( RESIZE_NEW( Employee, employees, employeeCount, 10 ) ) {

               (employeeCount is now 110 and
               employees may be a different address)

               ... do stuff with additional employees ...

               if ( RESIZE_NEW( Employee, employees, employeeCount,
                                -employeeCount / 2 ) ) {

                 (after 50% down-sizing, employeeCount is now 55)
                 ... do stuff with remaining employees ...
               }

             }
           }

           FREE( employees );

******************************************************************************/

Integer resizeMemory( Integer typeSize,
                      const void** existingAddress, Integer* existingCount,
                      Integer deltaCount, Integer zeroDelta ) {

  PRE08( typeSize > 0,
         existingAddress,
         existingCount,
         *existingCount >= 0,
         IMPLIES( *existingAddress,      *existingCount >  0 ),
         IMPLIES( *existingAddress == 0, *existingCount == 0 ),
         *existingCount + deltaCount >= 0,
         IS_BOOL( zeroDelta ) );

  CHECKING( const void* const OLD( existingAddress ) = *existingAddress; )
  CHECKING( const Integer     OLD( existingCount   ) = *existingCount;   )

  Integer ok = 1;
  const Integer oldBytes = ( *existingCount              ) * typeSize;
  const Integer newBytes = ( *existingCount + deltaCount ) * typeSize;

  DEBUG( fprintf( stderr, "resizeMemory( existingAddress = %p (%p), "
                          "existingCount = %p "
                          "(%"INTEGER_FORMAT"), "
                          "deltaCount = %"INTEGER_FORMAT", "
                          "zeroDelta = %"INTEGER_FORMAT" )",
                          existingAddress, *existingAddress,
                          existingCount, *existingCount,
                          deltaCount, zeroDelta ); )

  DEBUG( fprintf( stderr, " newBytes = %"INTEGER_FORMAT" (%lu)\n",
                  newBytes, (size_t) newBytes ); )

  if ( newBytes == 0 ) {
    /* Must use a temporary to free the const-cast-awayed pointer. */
    void* theExistingAddress = (void*) *existingAddress;
    FREE( theExistingAddress );
    *existingAddress = 0;
    *existingCount   = 0;
  } else if ( newBytes != oldBytes ) {
    void* newAddress = 0;
    Integer forceFailure = 0;

    /* Decide if simulating failure or not: */

    if ( failureCountDown > 0 ) { /* If we're counting down then... */
      --failureCountDown;
      forceFailure = failureCountDown == 0; /* Fail iff 0 is reached. */
    }

    if ( AND3( ! forceFailure, isSizet( newBytes ), available( newBytes ) ) ) {
      /* Don't attempt if too large to represent. */
      /* HACK: avoid calling realloc with pointer-to-0 - it is not portable. */

      newAddress = *existingAddress == 0 ? malloc( newBytes )
                   : realloc( (void*) *existingAddress, newBytes );

      DEBUG( fprintf( stderr, " yields newAddress = %p\n", newAddress ); )
    }

    if ( newAddress == 0 ) {
      ok = 0;
      failureMessage( "Can't re-allocate %"INTEGER_FORMAT" bytes to complete"
                      " the requested action.", newBytes );
    } else {
      *existingAddress = newAddress;
      newAddress = 0;

      if ( AND2( deltaCount > 0, zeroDelta ) ) {
        const Integer deltaBytes    = deltaCount     * typeSize;
        const Integer existingBytes = *existingCount * typeSize;
        char* const baseAddress     = (char*) *existingAddress;
        char* const extraPortion    = baseAddress +
                                      existingBytes / sizeof (char);

        DEBUG( fprintf( stderr, "Zeroing %"INTEGER_FORMAT" (%lu) bytes at"
                        " extraPortion = %p\n",
                        deltaBytes, (size_t) deltaBytes, extraPortion ); )

        memset( extraPortion, 0, deltaBytes );
      }

      *existingCount += deltaCount;
    }
  }

  POST05( *existingCount >= 0,
          IMPLIES( *existingAddress,      *existingCount >  0 ),
          IMPLIES( *existingAddress == 0, *existingCount == 0 ),
          IMPLIES( OR2( deltaCount == 0, ! ok ),
                   AND2( *existingAddress == OLD( existingAddress ),
                         *existingCount   == OLD( existingCount   ) ) ),
          IMPLIES( ok, *existingCount == OLD( existingCount ) + deltaCount ) );

  return ok;
}



/******************************************************************************
PURPOSE: setCountDownToFailMemory - Sets a count-down that will cause the
         'countDown'-th call to memory allocation to simulate failure
         to facilitate writing test code that will automatically achieve
         coverage of failure-handling paths.
INPUTS:  Integer countDown  The number of the subsequent call to memory
                            allocation that will 'fail'.
                            E.g., if countDown is 1 then the next call to
                            newMemory() or resizeMemory() will fail
                            (but subsequent calls will not simulate failure).
                            If 2 then the next call will not simulate failure
                            but the one after that will. Etc.
NOTES:   UGLY: This routine is not thread-safe.
******************************************************************************/

void setCountDownToFailMemory( Integer countDown ) {
  PRE0( countDown >= 0 );
  failureCountDown = countDown;
}



/*============================ PRIVATE FUNCTIONS ============================*/



#if defined( __APPLE__ )

#include <mach/mach.h>

static unsigned long long available_bytes( void ) {
  unsigned long long result = 0;
  const mach_port_t host_priv_port = mach_host_self();
  vm_size_t page_size = 0;

  if ( host_page_size( host_priv_port, &page_size ) == KERN_SUCCESS &&
       page_size > 0 ) {
    vm_statistics_data_t vm_stat;
    unsigned int host_count = HOST_VM_INFO_COUNT;

    if ( host_statistics( host_priv_port, HOST_VM_INFO,
                          (host_info_t) &vm_stat, &host_count )
         == KERN_SUCCESS && vm_stat.free_count > 0 ) {
      result = vm_stat.free_count * page_size;
    }
  }

  DEBUG( fprintf( stderr, "available_bytes = %llu\n", result ); )
  return result;
}

#elif defined( __sun )

#include <unistd.h> /* For sysconf(). */

static unsigned long long available_bytes( void ) {
  const unsigned long long page_size = sysconf( _SC_PAGESIZE );
  const unsigned long long pages = sysconf( _SC_AVPHYS_PAGES );
  const unsigned long long result = page_size * pages;
  DEBUG( fprintf( stderr, "available_bytes = %llu\n", result ); )
  return result;
}

#elif defined( __linux__ )

#include <unistd.h> /* For sysconf(). */

static unsigned long long available_bytes( void ) {
  const unsigned long long page_size = sysconf( _SC_PAGESIZE );
/*
  Due to problematic memory buffering in Linux,
  _SC_AVPHYS_PAGES greatly underestimates actual free memory
  so instead overestimate by using _SC_PHYS_PAGES (installed memory)
  and rely on malloc returning 0 when the overestimate exceeds the actual
  (unknowable?) free bytes.
  const unsigned long long pages = sysconf( _SC_AVPHYS_PAGES );
*/
  const unsigned long long pages = sysconf( _SC_PHYS_PAGES );
  const unsigned long long result = page_size * pages;
  DEBUG( fprintf( stderr, "available_bytes = %llu\n", result ); )
  return result;
}

#else

/* Unknown so return the maximum and let malloc return 0 if it is exceeded. */

#include <limits.h> /* For UULONG_MAX. */

static unsigned long long available_bytes( void ) { return ULLONG_MAX; }

#endif



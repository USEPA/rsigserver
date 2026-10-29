
#ifndef MEMORY_H
#define MEMORY_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Memory.h - Defines macros for memory allocation and declares
         related routines.
NOTES:   Implementation is based on calls the standard malloc() and free().
         Usage:
           instead of:
             #include <errno.h>
             #include <stdlib.h>
             Integer* p = (Integer*) malloc( n * sizeof (Integer) );
             if ( p ) { free( p ); p = 0; }
             else myFailureHandler( errno );
           use:
             #include <Failure.h>
             #include <Memory.h>
             failureHandler( myFailureHandler );
             Integer* p = NEW( Integer, n );
             FREE( p );
         Uses the failureMessage() routine from Failure.c
HISTORY: 1993/03, Todd Plessel, EPA/LM, Created.
STATUS:  reviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <string.h>        /* For memset(). */
#include <BasicNumerics.h> /* For Integer.  */

/*================================= MACROS ==================================*/

/* Allocate memory large enough to hold 'count' items of given 'type': */

#define NEW( type, count ) ((type*) newMemory( (count), sizeof (type), 0 ))

/* Allocate memory large enough to hold 'count' items of size 'size': */

#define NEW_SIZE( size, count ) newMemory( (count), (size), 0 )

/* Versions that allocate then zero-out the memory: */

#define NEW_ZERO( type, count ) \
((type*) newMemory( (count), sizeof (type), 1 ))

#define NEW_SIZE_ZERO( size, count ) newMemory( (count), (size), 1 )

/*
 * Re-size an existing allocated block by deltaCount:
 * Example:
 *
 *   Foo* fooArray = NEW( Foo, 100 );
 *
 *   if ( fooArray ) {
 *
 *     Integer fooCount = 100;
 *
 *     ...grow by 10 more Foo...
 *
 *     if ( RESIZE( &fooArray, &fooCount, 10 ) ) {
 *       ...fooCount is now 110...
 *     }
 *
 *     FREE( fooArray );
 *   }
 */

#define RESIZE( addressOfExistingAddress, addressOfExistingCount, deltaCount ) \
resizeMemory( sizeof **(addressOfExistingAddress), \
(const void**) (addressOfExistingAddress), \
(addressOfExistingCount), (deltaCount), 0 )

/* Re-size an existing allocated block and zero-out any extra portion: */

#define RESIZE_ZERO( addressOfExistingAddress, addressOfExistingCount, \
deltaCount ) \
resizeMemory( sizeof **(addressOfExistingAddress), \
(const void**) (addressOfExistingAddress), \
(addressOfExistingCount), (deltaCount), 1 )

/*
 * Re-size an existing allocated block by deltaCount:
 * Example:
 *
 *   const Integer elementSize = type == INTEGER ? sizeof (int) :sizeof (float);
 *   void* array = NEW_SIZE( elementSize, 100 );
 *
 *   if ( array ) {
 *
 *     Integer arrayCount = 100;
 *
 *     ...grow by 10 more ints or floats...
 *
 *     if ( RESIZE_SIZE( &array, elementSize, &arrayCount, 10 ) ) {
 *       ...arrayCount is now 110...
 *     }
 *
 *     FREE( array );
 *   }
 */

#define RESIZE_SIZE( addressOfExistingAddress, elementSize, \
addressOfExistingCount, deltaCount ) \
resizeMemory( elementSize, (const void**) (addressOfExistingAddress), \
(addressOfExistingCount), (deltaCount), 0 )

#define RESIZE_SIZE_ZERO( addressOfExistingAddress, elementSize, \
addressOfExistingCount, deltaCount ) \
resizeMemory( elementSize, (const void**) (addressOfExistingAddress), \
(addressOfExistingCount), (deltaCount), 1 )

/* Zero the memory of the pointed-to structure: */

#define ZERO_OBJECT(p) (((p) ? memset( (p), 0, sizeof *(p) ) : (void*) 0))

/* Deallocate the memory and zero the pointer: */

#define FREE( p ) ( ( (p) ? freeMemory(p) : (void) 0 ), (p) = 0 )

/* Zero the memory, deallocate the memory and zero the pointer: */

#define FREE_ZERO( p ) \
( ( (p) ? memset( (p), 0, sizeof *(p) ), freeMemory(p) : (void) 0 ), (p) = 0 )

/* Macro version of destructor zeros pointer to avoid dangling references: */

#define FREE_OBJECT( p ) (((p) ? (p)->free(p) : (void) 0), (p) = 0)

/*=============================== FUNCTIONS =================================*/

extern void* newMemory( Integer count, Integer sizeEach, Integer zeroIt );
extern void  freeMemory( void* address );
extern Integer resizeMemory( Integer typeSize,
                             const void** existingAddress,
                             Integer* existingCount,
                             Integer deltaCount,
                             Integer zeroExtra );
extern void setCountDownToFailMemory( Integer countdown );

#ifdef __cplusplus
}
#endif

#endif /* MEMORY_H */


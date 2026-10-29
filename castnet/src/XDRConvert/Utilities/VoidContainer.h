#ifndef VOIDCONTAINER_H
#define VOIDCONTAINER_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: VoidContainer.h - Declares types and routines for using containers of
         dynamically allocated void pointers (which the containers 'owns'
         after insertion).

NOTES:   This is an "interface" / pure virtual ABC (Abstract Base Class) for
         "descendent classes".

         Commands:

         FREE_OBJECT( VoidContainer* self )
         void free( VoidContainer* self )
         void removeAll( VoidContainer* self )
         void apply( const VoidContainer* self, VoidVisitor visitor )

         Queries:

         Integer invariant( const VoidContainer* self )
         Integer ok( const VoidContainer* self )
         Integer has( const VoidContainer* self, const void* itemPointer )
         Integer equal( const VoidContainer* self, const VoidContainer* other )
         Integer count( const VoidContainer* self )
         VoidComparer comparer( const VoidContainer* self )
         VoidVisitor deleter( const VoidContainer* self )

HISTORY: 2003/03: Todd Plessel, EPA/LM, Created.
STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <BasicNumerics.h> /* For Integer.                         */
#include <Memory.h>        /* To export FREE_OBJECT() for clients. */

/*================================== TYPES ==================================*/

typedef void (*VoidVisitor)( void* item );
/* PRE0( item ) */

typedef Integer (*VoidComparer)( const void* item1, const void* item2 );
/* PRE02( item2, item2 ); POST0( IN4( result, -1, 0, 1 ) ) */

typedef struct VoidContainer VoidContainer;

/*= PRIVATE =*/

#define DECLARE_VOIDCONTAINER_MEMBERS( Type ) \
  void (*free)( Type* self ); \
  void (*removeAll)( Type* self ); \
  void (*apply)( const Type* self, VoidVisitor visitor ); \
  Integer (*invariant)( const Type* self ); \
  Integer (*ok)( const Type* self ); \
  Integer (*has)( const Type* self, const void* itemPointer ); \
  Integer (*equal)( const Type* self, const Type* other ); \
  Integer (*count)( const Type* self ); \
  VoidComparer (*comparer)( const Type* self ); \
  VoidVisitor (*deleter)( const Type* self )

struct VoidContainer {
  DECLARE_VOIDCONTAINER_MEMBERS( VoidContainer );
};



#ifdef __cplusplus
}
#endif

#endif /* VOIDCONTAINER_H */




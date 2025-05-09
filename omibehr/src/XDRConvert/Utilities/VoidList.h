#ifndef VOIDLIST_H
#define VOIDLIST_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: VoidList.h - Declares types and routines for using doubly-linked
         lists of dynamically allocated (with NEW()) void pointers
         (which the list 'owns' after insertion).
         It supports O(n) traversal in either direction and O(1) insertion,
         deletion and access to items at either end and adjacent to the most
         recently accessed item (when assertions checking is disabled).

NOTES:   Commands:

         VoidList* newVoidList( VoidVisitor deleterCallback,
                                VoidComparer comparerCallback )
         FREE_OBJECT( VoidList* self )
         void free( VoidList* self )
         void removeAll( VoidList* self )
         void apply( const VoidList* self, VoidVisitor visitor )
         void insert( VoidList* self, void* anItem, Integer index )
         void remove( VoidList* self, Integer index )
         void replace( VoidList* self, void* anItem, Integer index )
         void sort( VoidList* self )

         Queries:

         Integer invariant( const VoidList* self )
         Integer ok( const VoidList* self )
         Integer has( const VoidList* self, const void* itemPointer )
         Integer equal( const VoidList* self, const VoidList* other )
         Integer count( const VoidList* self )
         VoidComparer comparer( const VoidList* self )
         VoidVisitor deleter( const VoidList* self )
         const void* item( const VoidList* self, Integer index )
         Integer sorted( const VoidList* self )
         Integer index( const VoidList* self, const void* anItem )



         Example usage:

         Creates, prints, edits and frees a list of numbers.

         #include <stdio.h>
         #include <VoidList.h>


         typedef struct {
           char*   lexeme;
           Integer value;
         } Number;


         static Number* newNumber( Integer value ) {
           Number* number = NEW_ZERO( Number, 1 );

           if ( number ) {
             number->value  = value;
             number->lexeme = NEW_ZERO( char, 80 );

             if ( number->lexeme ) {
               sprintf( number->lexeme, "%"INTEGER_FORMAT, value );
             } else {
               FREE( number );
             }
           }

           return number;
         }


         static void freeNumber( Number* number ) {
           FREE( number->lexeme );
           number->value = 0;
         }


         static Integer compareNumbers( const Number* a, const Number* b ) {
           const Integer result = a == b ?  0 :
                                  a == 0 ? -1 :
                                  b == 0 ?  1 :
                                  a->value < b->value ? -1 :
                                  a->value > b->value ?  1 :
                                  0;
           return result;
         }


         static void printNumber( const Number* number ) {
           printf( "number->lexeme = '%s', ", number->lexeme );
           printf( "number->value  = %"INTEGER_FORMAT"\n", number->value );
         }



         int main( void ) {
           VoidList* list = newVoidList( (VoidVisitor) freeNumber,
                                         (VoidComparer) compareNumbers );
           Integer ok = list != 0;

           if ( list ) {
             Integer index = 0;

             for ( index = 0; ok && index < 5; ++index ) {
               Number* number = newNumber( ( index + 1 ) * 10 + 1 );
               ok = number != 0;

               if ( ok ) {
                 list->insert( list, number, index );
                 ok = list->ok( list );
               }

               if ( ! ok && number != 0 ) {
                 freeNumber( number );
                 FREE( number );
               }
             }

             list->apply( list, (VoidVisitor) printNumber );

             for ( index = 0; index < list->count( list ); index += 2 ) {
               list->remove( list, index );
             }

             if ( list->count( list ) ) {
               const Number* const last = list->item( list, LAST_ITEM );
               Number* other = newNumber( last->value - 20 );

               if ( other ) {
                 const Integer index = list->index( list, other );

                 if ( index != NOT_FOUND ) {
                   list->replace( list, other, LAST_ITEM );
                 } else {
                   freeNumber( other );
                   FREE( other );
                 }
               }
             }

             if ( list->count( list ) ) {
               const VoidContainer* container = (const VoidContainer*) list;
               container->apply( container, (VoidVisitor) printNumber );
             }
           }

           FREE_OBJECT( list );
           return ! ok;
         }


         cc -64 -mips4 -xansi -fullwarn -g \
            -I/usr/include -I/usr/include/sys -I../../../../../include -I.. \
            -o test_list test_list.c \
            -L/usr/lib64/mips4 -L/usr/lib64 -L../../../../../lib/IRIX64 \
            -lVoidList.debug \
            -lMemory.debug -lFailure.debug -lBasicNumerics.debug \
            -lm -lc

         test_list

         number->lexeme = '11', number->value  = 11
         number->lexeme = '21', number->value  = 21
         number->lexeme = '31', number->value  = 31
         number->lexeme = '41', number->value  = 41
         number->lexeme = '51', number->value  = 51
         number->lexeme = '21', number->value  = 21
         number->lexeme = '31', number->value  = 31
         number->lexeme = '31', number->value  = 31

HISTORY: 1996/06: Todd Plessel, EPA/LM, Created.
STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <VoidContainer.h> /* For DECLARE_VOIDCONTAINER_MEMBERS. */

/*================================== TYPES ==================================*/

enum { LAST_ITEM = -1 }; /* To access the last item without using count(). */
enum { NOT_FOUND = -1 }; /* Returned by index() if item is not found.      */

typedef struct VoidList VoidList;

/* Constructor: */

extern VoidList* newVoidList( VoidVisitor deleterCallback,
                              VoidComparer comparerCallback );

/* Destructor: use macro FREE_OBJECT( yourVoidList ); */

/*= PRIVATE =*/

typedef struct VoidListPrivate VoidListPrivate;

#define DECLARE_VOIDLIST_MEMBERS( Type ) \
  DECLARE_VOIDCONTAINER_MEMBERS( Type ); \
  void (*insert)( VoidList* self, void* anItem, Integer index ); \
  void (*remove)( VoidList* self, Integer index ); \
  void (*replace)( VoidList* self, void* anItem, Integer index ); \
  const void* (*item)( const VoidList* self, Integer index ); \
  void (*sort)( VoidList* self ); \
  Integer (*sorted)( const VoidList* self ); \
  Integer (*index)( const VoidList* self, const void* anItem ); \
  VoidListPrivate* p

struct VoidList {
  DECLARE_VOIDLIST_MEMBERS( VoidList );
};



#ifdef __cplusplus
}
#endif

#endif /* VOIDLIST_H */




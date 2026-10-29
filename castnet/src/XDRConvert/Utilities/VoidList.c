
/******************************************************************************
PURPOSE: VoidList.c - Implements doubly-linked lists of owned void pointers.
NOTES:   See VoidList.h.
HISTORY: 1996/06: Todd Plessel, EPA/LM, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>         /* For stderr, fprintf().                   */
#include <string.h>        /* For memset().                            */

#include <Assertions.h>    /* For PRE(), POST(), CHECK(), NON_ZERO*(). */
#include <BasicNumerics.h> /* For Integer.                             */
#include <Memory.h>        /* For macros NEW_ZERO(), FREE_ZERO().      */
#include <VoidList.h>      /* For public interface.                    */

/*================================== TYPES ==================================*/

typedef struct VoidListNode VoidListNode;

struct VoidListNode {
  void*         item;     /* Pointer to client's NEW'd data item.   */
  VoidListNode* next;     /* Pointer to the next node.              */
  VoidListNode* previous; /* Pointer to the previous node.          */
};

struct VoidListPrivate {
  VoidListNode*    head;       /* Points to the first node.              */
  VoidListNode*    tail;       /* Points to the last  node.              */
  VoidListNode*    cache;      /* Points to the last-accessed node.      */
  Integer          cacheIndex; /* Index of  the last-accessed node.      */
  Integer          count;      /* The number of items in the list.       */
  Integer          ok;         /* Did last command succeed?              */
  VoidVisitor  deleter;        /* Optional client callback invoked on    */
                               /* each non-0 stored item before freeing. */
  VoidComparer comparer;   /* Optional client callback invoked by index. */
};

/*========================== FORWARD DECLARATIONS ===========================*/

/* Member functions: */

/* Queries: */

static Integer invariant( const VoidList* self );
static Integer ok(        const VoidList* self );
static Integer sorted(    const VoidList* self );
static Integer has(       const VoidList* self, const void* item );
static Integer equal(     const VoidList* self, const VoidList* other );
static Integer count(     const VoidList* self );
static Integer index__(   const VoidList* self, const void* anItem );
static const void* item(  const VoidList* self, Integer index );
static VoidComparer comparer( const VoidList* self );
static VoidVisitor  deleter(  const VoidList* self );

/* Commands: */

static void free__(    VoidList* self );
static void insert(    VoidList* self, void* anItem, Integer index);
static void remove__(  VoidList* self, Integer index );
static void removeAll( VoidList* self );
static void replace(   VoidList* self, void* anItem, Integer index);
static void sort(      VoidList* self );
static void apply( const VoidList* self, VoidVisitor visitor );

/* Helper functions: */

static Integer hasPrivateDataAndMemberFunctions( const VoidList* self );
static void swapPrivateContents( VoidList* a, VoidList* b );
static VoidListNode* newNode( void* item );
static void freeNode( VoidListNode* nodeToFree, VoidVisitor deleter );
static VoidListNode* indexedNode( const VoidList* self, Integer index );
static VoidListNode* cachedNode( const VoidList* self, Integer index );
static void setCache( const VoidList* self, Integer index );
static void insertNode( VoidList* self, VoidListNode* node, Integer index );
static void removeNode( VoidList* self, VoidListNode* node );
static void moveNode( VoidList* source, VoidListNode* node, VoidList* result );
static void mergeSublists( VoidList* source, VoidList* result, Integer length);

#ifdef DEBUGGING
static void print( const VoidList* self );
#endif /* DEBUGGING */



/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: newVoidList - Allocate a new empty VoidList.
INPUTS:  VoidVisitor deleterCallback  Optional client callback called on each
                                      stored item prior to freeing it.
                                      The routine should not free its argument.
                                      Its purpose is to release resources
                                      (e.g., sub-object pointers) held by the
                                      object.
                                      VoidStack will call FREE_ZERO() on the
                                      argument.
         VoidComparer comparerCallback  Optional client callback called on each
                                        stored item during equal().
                                        Semantics like strcmp().
RETURNS: VoidList* New VoidList or 0 if unsuccessful (due to memory allocation
                    failure in which case failureMessage() is called).
******************************************************************************/

VoidList* newVoidList( VoidVisitor deleterCallback,
                       VoidComparer comparerCallback ) {
  VoidList* result = NEW_ZERO( VoidList, 1 );

  if ( result ) {
    VoidListPrivate* const p = NEW_ZERO( VoidListPrivate, 1 );

    if ( p == 0 ) {
      FREE( result );
    } else {
      p->cacheIndex = NOT_FOUND;
      p->ok         = 1;
      p->deleter    = deleterCallback;
      p->comparer   = comparerCallback;

      result->free      = free__;
      result->insert    = insert;
      result->remove    = remove__;
      result->removeAll = removeAll;
      result->replace   = replace;
      result->sort      = sort;
      result->apply     = apply;
      result->invariant = invariant;
      result->ok        = ok;
      result->sorted    = sorted;
      result->has       = has;
      result->equal     = equal;
      result->count     = count;
      result->index     = index__;
      result->item      = item;
      result->comparer  = comparer;
      result->deleter   = deleter;

      result->p = p;
    }
  }

  POST0( IMPLIES( result,
                  AND5( result->invariant( result ),
                        result->ok( result ),
                        result->count( result ) == 0,
                        result->comparer( result ) == comparerCallback,
                        result->deleter( result ) == deleterCallback ) ) );
  return result;
}



/******************************************************************************
PURPOSE: free - Deallocate a VoidList.
INPUTS:  VoidList* self  The list to free.
NOTES:   If the deleter callback exists, it is called on each item
         prior to freeing it with FREE_ZERO(). Then the list itself is freed.
         Called by macro FREE_OBJECT( self ) (which also sets self = 0 to
         avoid dangling pointers in client code).
******************************************************************************/

static void free__( VoidList* self ) {
  PRE( self );
  VoidListPrivate* p = self->p;

  while ( p->head ) {
    VoidListNode* const nodeToFree = p->head;
    p->head = p->head->next;
    freeNode( nodeToFree, p->deleter );
  }

  FREE_ZERO( p );
  FREE_ZERO( self );
  POST0( self == 0 );
}



/******************************************************************************
PURPOSE: invariant - Checks if the object is in a valid state.
INPUTS:  const VoidList* self  The object to check.
RETURNS: Integer 1 if non-0 and valid, else 0.
NOTES:   If this routine ever returns a value other than 1, then this indicates
         that the code contains a defect. Used in PRE() and POST().
******************************************************************************/

static Integer invariant( const VoidList* self ) {
  Integer result = hasPrivateDataAndMemberFunctions( self );

  DEBUG( if ( self ) print( self ); )

  if ( result ) {
    const VoidListPrivate* const p = self->p;

    if ( p->count == 0 ) {
      result = AND2( IS_ZERO3( p->head, p->tail, p->cache ),
                     p->cacheIndex == NOT_FOUND );
    } else { /* Count the nodes. */
      Integer theCount   = 0;
      Integer foundTail  = 0; /* Encountered tail        while traversing? */
      Integer foundCache = 0; /* Encountered cached node while traversing? */
      Integer foundCacheAt = NOT_FOUND; /* Index of cached node found. */
      const VoidListNode* node         = p->head;
      const VoidListNode* previousNode = 0;

      while ( node ) { /* Check for incorrectly-linked nodes: */

        if ( AND6( node->item     != 0,
                   node->next     != p->head,
                   node->previous != p->tail,
                   node->next     != node,
                   node->previous != node,
                   node->previous == previousNode ) ) {

          if ( node == p->tail ) {
            ++foundTail;
          }

          if ( node == p->cache ) {
            ++foundCache;
            foundCacheAt = theCount;
          }

          previousNode = node;
          node         = node->next;
          ++theCount;

          if ( theCount > p->count ) { /* Defect: cycle in list! */
            node = 0;
          }
        } else { /* Defect: incorrectly linked node! */
          node = 0;
        }
      }

      result = AND8( IS_BOOL( p->ok ),
                     theCount          == p->count,
                     foundTail         == 1,
                     p->head->previous == 0,
                     p->tail->next     == 0,
                     p->head->item     != 0,
                     p->tail->item     != 0,
                     IMPLIES_ELSE( p->cache,
                                   AND2( foundCache   == 1,
                                         foundCacheAt == p->cacheIndex ),
                                   AND2( foundCache   == 0,
                                         foundCacheAt == NOT_FOUND ) ) );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: ok - Did the last command succeed?
INPUTS:  const VoidList* self  The list to check.
RETURNS: Integer 1 or 0.
******************************************************************************/

static Integer ok( const VoidList* self ) {
  PRE( self );
  const Integer result = self->p->ok;
  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: count - The number of items stored in the list.
INPUTS:  const VoidList* self  The list to check.
RETURNS: Integer The number of items in the list.
NOTES:   This operation is O(1) (when assertions checking is disabled).
******************************************************************************/

static Integer count( const VoidList* self ) {
  PRE( self );
  const Integer result = self->p->count;
  POST( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: item - Item at the given index in the list.
INPUTS:  const VoidList* self  The list to scan.
         Integer         index The index of the item.
RETURNS: const void* The item at the given index.
NOTES:   If index is a middle item then the index is cached so that subsequent
         accesses to this or adjacent items are O(1) (when assertions
         checking is disabled).
******************************************************************************/

static const void* item( const VoidList* self, Integer index ) {
  PRE2( self,
        OR2( index == LAST_ITEM,
             IN_RANGE( index, 0, self->count( self ) - 1 ) ) );
  const Integer theIndex = index == LAST_ITEM ? self->p->count - 1 : index;
  const void* const result = indexedNode( self, theIndex )->item;
  POST( result );
  return result;
}



/******************************************************************************
PURPOSE: has - Does the list already contain this pointed to item?
INPUTS:  const VoidList* self         The list to scan.
         const void*     itemPointer  The item pointer to search for.
RETURNS: Integer 1 if found else 0.
******************************************************************************/

static Integer has( const VoidList* self, const void* itemPointer ) {
  PRE2( self, itemPointer );
  Integer result = 0;
  const VoidListNode* node = self->p->head;

  for ( ; AND2( node, result == 0 ); node = node->next ) {
    result = node->item == itemPointer;
  }

  POST( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: index - Index of first matching item in the list.
INPUTS:  const VoidList* self  The list to scan.
         const void*     item  The item to match (via comparer()).
RETURNS: Integer  The index of the first occurrence of item or NOT_FOUND.
NOTES:   The client comparer callback must exist.
         This routine is an O(n) linear search.
         If the resulting index is a middle item then the index is cached so
         subsequent accesses via item() to the resulting indexed item
         (or adjacent items) are O(1) (when assertions checking is disabled).
******************************************************************************/

static Integer index__( const VoidList* self, const void* item ) {
  PRE3( self, item, self->comparer( self ) );
  Integer result = 0;
  VoidListPrivate* const p = self->p;
  VoidListNode* node = p->head;
  VoidListNode* foundNode = 0;
  VoidComparer comparer = p->comparer;

  /* Search the list: */

  while ( node ) {

    if ( comparer( node->item, item ) == 0 ) {
      foundNode = node;
      node = 0;
    } else {
      node = node->next,
      ++result;
    }
  }

  if ( foundNode ) {

    /* If the node found is not the head or tail then cache it: */

    if ( ! IN3( result, 0, p->count - 1 ) ) {
      p->cache      = foundNode;
      p->cacheIndex = result;
    }
  } else {
    result = NOT_FOUND;
  }

/* Define this macro for an encapsulated readable checkable postcondition. */

#define ITEM_IS_CACHED( index ) \
AND2( self->p->cacheIndex == index, self->p->cache == foundNode )

  POST( IMPLIES( result != NOT_FOUND,
                 AND3( IN_RANGE( result, 0, self->count( self ) - 1 ),
                       self->comparer(self)(self->item(self,result),item) == 0,
                       IMPLIES( ! IN3( result, 0, self->count( self ) - 1 ),
                                ITEM_IS_CACHED( result ) ) ) ) );

  return result;
}



/******************************************************************************
PURPOSE: comparer - Client-supplied item comparison function used by index()
         and equal().
INPUTS:  const VoidList* self The list to query.
RETURNS: VoidComparer  The comparison function or 0 if it was not set.
NOTES:   The client compare callback is set during object construction and has
         semantics like strcmp().
******************************************************************************/

static VoidComparer comparer( const VoidList* self ) {
  PRE( self );
  VoidComparer result = self->p->comparer;
  return result;
}



/******************************************************************************
PURPOSE: deleter - Client-supplied deleter function called on removed items
         before they are freed (with FREE_ZERO()).
INPUTS:  const VoidList* self The list to query.
RETURNS: VoidVisitor  The deleter function or 0 if it was not set.
NOTES:   The client deleter callback is set during object construction.
******************************************************************************/

static VoidVisitor deleter( const VoidList* self ) {
  PRE( self );
  VoidVisitor result = self->p->deleter;
  return result;
}



/******************************************************************************
PURPOSE: equal - Do two lists have equivalent items according to comparer()?
INPUTS:  const VoidList* self  The 1st list to compare.
         const VoidList* other The 2nd list to compare.
RETURNS: Integer 1 or 0.
NOTES:   The client comparer callback must be set during object construction
         and must be the same for both lists.
******************************************************************************/

static Integer equal( const VoidList* self, const VoidList* other ) {
  PRE5( self,
        other,
        other->invariant( other ),
        self->comparer( self ),
        self->comparer( self ) == other->comparer( other ) );
  Integer result = self == other;

  if ( ! result ) {
    const VoidListPrivate* const p = self->p;
    const VoidListPrivate* const other_p = other->p;
    result = p->count == other_p->count;

    if ( AND2( result, p->count ) ) { /* Same non-0 count so compare items: */
      const VoidListNode* selfNode  = p->head;
      const VoidListNode* otherNode = other_p->head;
      VoidComparer comparer         = p->comparer;

      do {
        result    = comparer( selfNode->item, otherNode->item ) == 0;
        selfNode  = selfNode->next;
        otherNode = otherNode->next;
      } while ( AND2( result, selfNode ) );
    }
  }

  POST3( IS_BOOL( result ),
         IMPLIES( self == other, result ),
         IMPLIES( result, self->count( self ) == other->count( other ) ) );
  return result;
}



/******************************************************************************
PURPOSE: sorted - Are the items ordered by comparer()?.
INPUTS:  const VoidList* self  The list to check.
RETURNS: Integer 1 or 0.
******************************************************************************/

static Integer sorted( const VoidList* self ) {
  PRE2( self, self->comparer( self ) );
  const VoidListPrivate* const p = self->p;
  Integer result = p->count < 2;

  if ( ! result ) {
    const VoidListNode* previousNode = p->head;
    const VoidListNode* node         = previousNode->next;
    VoidComparer comparer            = p->comparer;

    do {
      result       = comparer( previousNode->item, node->item ) <= 0;
      previousNode = node;
      node         = node->next;
    } while ( AND2( result, node ) );
  }

  POST3( IS_BOOL( result ),
         IMPLIES( self->count( self ) < 2, result ),
         IMPLIES( AND2( result, self->count( self ) >= 2 ),
                  self->comparer( self )( self->item( self, 0 ),
                                          self->item( self, LAST_ITEM)) <= 0));
  return result;
}



/******************************************************************************
PURPOSE: sort - Sort the items in ascending order by comparer().
INPUTS:  VoidList* self  The list to sort.
NOTES:   Implements an O(n*log2(n)) bottom-up (non-recursive) mergesort
         in-place (in constant space) scanning the list in strictly sequential
         order and moving only links between self and a temporary list.
******************************************************************************/

static void sort( VoidList* self ) {
  PRE2( self, self->comparer( self ) );
  CHECKING( const Integer OLD(count) = self->p->count; )
  const Integer theCount = self->p->count;
  Integer length = 1;          /* Of each sorted sub-list to merge.   */
  VoidListPrivate tempPrivate; /* On stack so no allocation failures. */
  VoidList temp = *self;       /* Copy member functions.              */
  temp.p = &tempPrivate;       /* But don't share private data.       */
  ZERO_OBJECT( &tempPrivate ); /* Initialize temp's private data:     */
  tempPrivate.cacheIndex = NOT_FOUND;
  tempPrivate.ok         = 1;
  tempPrivate.deleter    = self->p->deleter;
  tempPrivate.comparer   = self->p->comparer;

  for ( length = 1; length < theCount; length += length ) {
    mergeSublists( self, &temp, length ); /* Split self & merge it into temp.*/
    swapPrivateContents( self, &temp );   /* Move nodes from temp to self.   */
  }

  self->p->ok = 1;
  POST3( self->ok( self ),
         self->count( self ) == OLD(count),
         self->sorted( self ) );
}



/******************************************************************************
PURPOSE: apply - Apply visitor to each of the items in order.
INPUTS:  const VoidList* self    The list to visit.
         VoidListvisitor visitor The routine to apply (call) passing the item.
NOTES:   Example:

         static void printNumber( const Number* number ) {
           fprintf( stderr, "number->lexeme = '%s', ", number->lexeme );
           fprintf( stderr, "number->value  = %"INTEGER_FORMAT"\n", number->value );
         }
         ...
         list->apply( list, (VoidVisitor) printNumber ); // Print numbers.
******************************************************************************/

static void apply( const VoidList* self, VoidVisitor visitor ) {
  PRE2( self, visitor );
  CHECKING( const Integer OLD(count) = self->p->count; )
  VoidListPrivate* const p = self->p;
  const VoidListNode* node = p->head;

  for ( ; node; node = node->next ) {
    visitor( node->item );
  }

  p->ok = 1;
  POST2( self->ok( self ), self->count( self ) == OLD(count) );
}



/******************************************************************************
PURPOSE: insert - Insert an item at the given index.
INPUTS:  VoidList* self  The list to insert into.
         void*     item  The item to insert.
                         This must be a pointer obtained from NEW() (or
                         NEW_ZERO()) since it will eventually be freed using
                         FREE_ZERO() when the node is deleted (either replaced,
                         removed or upon list destruction).
         Integer   index The index to insert the item at.
NOTES:   Use self->ok( self ) to determine if insertion was successful.
         If unsuccessful (due to memory allocation failure) then
         failureMessage() is called.
         Insertion at either end of the list or at or adjacent to a
         cached node is O(1) (when assertions checking is disabled).
******************************************************************************/

static void insert( VoidList* self, void* item, Integer index ) {
  PRE4( self,
        item,
        ! self->has( self, item ),
        OR2( index == LAST_ITEM,
             IN_RANGE( index, 0, self->count( self ) ) ) );
  CHECKING( const Integer OLD(count) = self->p->count; )
  VoidListPrivate* const p = self->p;
  const Integer theIndex = index == LAST_ITEM ? p->count : index;
  VoidListNode* const theNewNode = newNode( item ); /* Has zeroed links. */
  p->ok = theNewNode != 0;

  if ( p->ok ) {
    insertNode( self, theNewNode, theIndex );
  }

  POST( IMPLIES_ELSE( self->ok( self ),
                      AND2( self->count( self ) == OLD(count) + 1,
                            self->item( self, index ) == item ),
                      self->count( self ) == OLD(count) ) );
}



/******************************************************************************
PURPOSE: remove - Remove the item at the given index.
INPUTS:  VoidList* self  The list to remove from.
         Integer   index The index of the item to remove.
NOTES:   If the deleter callback exists, it is called on item prior to freeing
         it with FREE_ZERO(). Removal at either end of the list is O(1).
******************************************************************************/

static void remove__( VoidList* self, Integer index ) {
  PRE3( self,
        self->count( self ),
        OR2( index == LAST_ITEM,
             IN_RANGE( index, 0, self->count( self ) - 1 ) ) );
  CHECKING( const Integer OLD(count) = self->p->count; )
  VoidListPrivate* const p = self->p;
  const Integer theIndex = index == LAST_ITEM ? p->count - 1 : index;
  VoidListNode* const node = indexedNode( self, theIndex );
  removeNode( self, node );
  p->ok = 1;
  freeNode( node, p->deleter );
  POST2( self->ok( self ), self->count( self ) == OLD(count) - 1 );
}



/******************************************************************************
PURPOSE: removeAll - Remove all items (if any).
INPUTS:  VoidList* self  The list to remove from.
NOTES:   If the deleter callback exists, it is called on each item prior to
         freeing it with FREE_ZERO().
******************************************************************************/

static void removeAll( VoidList* self ) {
  PRE( self );
  VoidListPrivate* const p = self->p;

  while ( p->count ) {
    VoidListNode* const node = p->head;
    removeNode( self, node );
    freeNode( node, p->deleter );
  }

  p->ok = 1;
  POST2( self->ok( self ), self->count( self ) == 0 );
}



/******************************************************************************
PURPOSE: replace - Replace the item at the given index.
INPUTS:  VoidList* self  The list to replace an item in.
         void*     item  The replacement item.
                         This must be a pointer obtained from NEW() (or
                         NEW_ZERO()) since it will eventually be freed using
                         FREE_ZERO() when the node is destroyed (either
                         replaced, removed or upon list destruction).
         Integer   index The index of the item to replace.
NOTES:   If the deleter callback exists, it is called on the item to be
         replaced prior to freeing it with FREE().
******************************************************************************/

static void replace( VoidList* self, void* item, Integer index ) {
  PRE5( self,
        item,
        ! self->has( self, item ),
        self->count( self ),
        OR2( index == LAST_ITEM,
             IN_RANGE( index, 0, self->count( self ) - 1 ) ) );
  CHECKING( const Integer OLD(count) = self->p->count; )
  VoidListPrivate* const p = self->p;
  const Integer theIndex = index == LAST_ITEM ? p->count - 1 : index;
  VoidVisitor deleter = p->deleter;
  VoidListNode* const nodeToReplace = indexedNode( self, theIndex );
  void* itemToReplace = nodeToReplace->item;

  if ( deleter ) {
    CHECK( itemToReplace );
    deleter( itemToReplace );
  }

  FREE( itemToReplace );
  nodeToReplace->item = item;
  p->ok = 1;
  POST3( self->ok( self ),
         self->count( self ) == OLD(count),
         self->item( self, index ) == item );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: hasPrivateDataAndMemberFunctions - Is object structurally valid?
INPUTS:  const VoidList* self  The list to verify.
RETURNS: Integer 1 if valid, else 0.
NOTES:   Used by invariant().
******************************************************************************/

static Integer hasPrivateDataAndMemberFunctions( const VoidList* self ) {
  const Integer result =
    AND21( self,
           self->p,
           self->p->count >= 0,
           IS_BOOL( self->p->ok ),
           self->free      == free__,
           self->insert    == insert,
           self->remove    == remove__,
           self->removeAll == removeAll,
           self->replace   == replace,
           self->sort      == sort,
           self->apply     == apply,
           self->invariant == invariant,
           self->ok        == ok,
           self->sorted    == sorted,
           self->has       == has,
           self->equal     == equal,
           self->count     == count,
           self->index     == index__,
           self->item      == item,
           self->comparer  == comparer,
           self->deleter   == deleter );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: swapPrivateContents - Swap the private contents of two lists.
INPUTS:  VoidList* a  1st list operand.
         VoidList* b  2nd list operand.
OUTPUTS: VoidList* a  1st list operand.
         VoidList* b  2nd list operand.
******************************************************************************/

static void swapPrivateContents( VoidList* a, VoidList* b ) {
  PRE02( a, b );
  const VoidListPrivate temp = *a->p;
  *a->p = *b->p;
  *b->p = temp;
}



/******************************************************************************
PURPOSE: newNode - Return an allocated initialized VoidListNode.
INPUTS:  void* item     The node data item to initialize to.
RETURNS: VoidListNode*  The newly allocated node (or 0 if failed).
NOTES:   If unsuccessful (due to memory allocation failure) then
         failureMessage() is called.
******************************************************************************/

static VoidListNode* newNode( void* item ) {
  PRE0( item );
  VoidListNode* const result = NEW_ZERO( VoidListNode, 1 );

  if ( result ) {
    result->item = item;
  }

  POST0( IMPLIES( result,
                  AND2( result->item == item,
                        IS_ZERO2( result->next, result->previous ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: freeNode - Deallocate a VoidListNode.
INPUTS:  VoidListNode* nodeToFree  The node to free.
         VoidVisitor deleter   Optional client's destroy callback.
NOTES:   If the deleter callback exists, it is called on this item prior to
         freeing it with FREE().
******************************************************************************/

static void freeNode( VoidListNode* nodeToFree, VoidVisitor deleter ) {
  PRE0( nodeToFree );
  void* itemOfNodeToFree = nodeToFree->item;

  if ( deleter ) {
    CHECK( itemOfNodeToFree );
    deleter( itemOfNodeToFree );
  }

  FREE( itemOfNodeToFree );
  FREE_ZERO( nodeToFree );
  POST0( nodeToFree == 0 );
}



/******************************************************************************
PURPOSE: indexedNode - Return the indexed node of a list.
INPUTS:  const VoidList* self  The list to scan.
         Integer         index The index of the node to return.
RETURNS: VoidListNode*   The node at the given index.
NOTES:   Caches the indexed node if necessary.
******************************************************************************/

static VoidListNode* indexedNode( const VoidList* self, Integer index ) {
  PRE2( self, IN_RANGE( index, 0, self->count( self ) - 1 ) );
  VoidListPrivate* const p = self->p;
  VoidListNode* result = cachedNode( self, index );

  if ( result == 0 ) {
    setCache( self, index );
    result = p->cache;
  }

  POST7( result,
         result->item,
         IMPLIES( index == 0,                  result == p->head ),
         IMPLIES( index == self->p->count - 1, result == self->p->tail ),
         IMPLIES( index == p->cacheIndex,      result == p->cache ),
         IMPLIES( AND2( p->cacheIndex >= 0, index == p->cacheIndex + 1 ),
                  result == p->cache->next ),
         IMPLIES( AND2( p->cacheIndex >= 0, index == p->cacheIndex - 1 ),
                  result == p->cache->previous ) );

  return result;
}



/******************************************************************************
PURPOSE: cachedNode - Returns a node at the given index if it is cached or
         otherwise obtainable without iteration otherwise returns 0.
INPUTS:  const VoidList* self  The list to check.
         Integer         index The index of the node to return.
RETURNS: VoidListNode*   The cached node or 0 if not available.
NOTES:   Helper for indexedNode().
******************************************************************************/

static VoidListNode* cachedNode( const VoidList* self, Integer index ) {
  PRE2( self, IN_RANGE( index, 0, self->count( self ) - 1 ) );
  VoidListPrivate* const p = self->p;
  VoidListNode* result = 0;

  if ( index == 0 ) {
    result = p->head;
  } else if ( index == p->count - 1 ) {
    result = p->tail;
  } else if ( p->cache ) {
    VoidListNode* const cache = p->cache;
    const Integer cacheIndex = p->cacheIndex;

    if ( index == cacheIndex ) {
      result = cache;
    } else if ( index == cacheIndex + 1 ) {
      result = cache->next;
    } else if ( index == cacheIndex - 1 ) {
      result = cache->previous;
    }
  }

  POST3( IMPLIES( result, result->item ),
         IMPLIES( index == 0,                  result == p->head ),
         IMPLIES( index == self->p->count - 1, result == self->p->tail ) );

  return result;
}



/******************************************************************************
PURPOSE: setCache - Sets a list's cache to the given index.
INPUTS:  const VoidList* self  The list to set the cache of.
         Integer         index The index of the node to set the cache to.
OUTPUTS: const VoidList* self  The list with cache set at the given index.
NOTES:   Helper for indexedNode().
******************************************************************************/

static void setCache( const VoidList* self, Integer index ) {
  PRE2( self, IN_RANGE( index, 0, self->count( self ) - 1 ) );
  VoidListPrivate* const p = self->p;
  const Integer distanceToFirst = index;
  const Integer distanceToLast  = p->count - 1 - index;
  const Integer distanceToCache = index > p->cacheIndex ?
                                  index - p->cacheIndex :
                                  p->cacheIndex - index;

  if ( distanceToFirst < distanceToLast ) {

    if ( distanceToFirst < distanceToCache ) {
      p->cache      = p->head;
      p->cacheIndex = 0;
    }
  } else if ( distanceToLast < distanceToCache ) {
    p->cache      = p->tail;
    p->cacheIndex = p->count - 1;
  }

  while ( p->cacheIndex < index ) {
    p->cache = p->cache->next;
    ++p->cacheIndex;
  }

  while ( p->cacheIndex > index ) {
    p->cache = p->cache->previous;
    --p->cacheIndex;
  }

  POST( p->cacheIndex == index );
}



/******************************************************************************
PURPOSE: insertNode - Insert the node at the given index.
INPUTS:  VoidList*     self  The list to insert into.
         VoidListNode* node  The node to insert.
         Integer       index Index to insert the node at.
NOTES:   Cache is preserved when inserting at either end of the list otherwise
         when inserting in the middle, the inserted node is cached.
******************************************************************************/

static void insertNode( VoidList* self, VoidListNode* node, Integer index ) {
  PRE06( self, node, node->next == 0, node->previous == 0, node->item,
         IN_RANGE( index, 0, self->p->count ) );
  CHECKING( const Integer OLD( count ) = self->p->count; )
  VoidListPrivate* const p = self->p;

  if ( index == 0 ) { /* Insert the new node at the head of the list: */

    if ( p->head ) {
      node->next = p->head;
      p->head->previous = node;
    }

    p->head = node;

    if ( p->tail == 0 ) {
      p->tail = node;
    }

    if ( p->cache ) {
      p->cache = p->cache->previous;
      CHECK( p->cache );
    }
  } else if ( index == p->count ) { /* Insert at tail of (non-empty) list: */
    p->tail->next  = node;
    node->previous = p->tail;
    p->tail        = node;
  } else { /* Insert in the middle of the list (with two or more nodes): */
    VoidListNode* const previousNode = indexedNode( self, index - 1 );
    node->next           = previousNode->next;
    previousNode->next   = node;
    node->previous       = previousNode;
    node->next->previous = node;
    p->cache             = node;
    p->cacheIndex        = index;
  }

  ++p->count; /* Count the new node. */

  POST0( self->p->count == OLD( count ) + 1 );
}



/******************************************************************************
PURPOSE: removeNode - Unlink the node from the list.
INPUTS:  VoidList* self      The list to remove from.
         VoidListNode* node  The node to remove.
NOTES:   Clears the cache (since to maintain it would require iteration).
******************************************************************************/

static void removeNode( VoidList* self, VoidListNode* node ) {
  PRE02( self, node );
  CHECKING( const Integer OLD( count ) = self->p->count; )
  CHECKING( const void* OLD( nodeItem ) = node->item; )
  VoidListPrivate* const p = self->p;
  VoidListNode* const previousNode = node->previous;

  /* Adjust head, tail and count and clear the cache: */

  if ( node == p->head ) {
    p->head = p->head->next;
  }

  if ( node == p->tail ) {
    p->tail = previousNode;
  }

  p->cache = 0;
  p->cacheIndex = NOT_FOUND;
  --p->count;

  /* Unlink the node to remove: */

  if ( previousNode ) {
    previousNode->next = node->next;
  }

  if ( node->next ) {
    node->next->previous = node->previous;
  }

  node->next     = 0;
  node->previous = 0;

  POST06( node->item     == OLD( nodeItem ),
          node->next     == 0,
          node->previous == 0,
          self->p->count == OLD( count ) - 1,
          self->p->cache == 0,
          self->p->cacheIndex == NOT_FOUND );
}



/******************************************************************************
PURPOSE: moveNode - Unlink the node from source and append it to result.
INPUTS:  VoidList* source    The list to remove from.
         VoidListNode* node  The node to remove from source.
OUTPUTS: VoidList* source    The list with node removed.
         VoidList* result    The list with node appended.
NOTES:   Helper for mergeSublists().
******************************************************************************/

static void moveNode(VoidList* source, VoidListNode* node, VoidList* result) {
  PRE06( source,
         source->invariant( source ),
         node,
         node->item,
         result,
         result->invariant( result ) );
  CHECKING( const Integer OLD( sourceCount ) = source->p->count; )
  CHECKING( const Integer OLD( resultCount ) = result->p->count; )
  CHECKING( const void* OLD( nodeItem ) = node->item; )

  removeNode( source, node );
  insertNode( result, node, result->p->count );

  POST06( source->invariant( source ),
          result->invariant( result ),
          node == result->p->tail,
          node->item == OLD( nodeItem ),
          source->p->count == OLD( sourceCount ) - 1,
          result->p->count == OLD( resultCount ) + 1 );
}



/******************************************************************************
PURPOSE: mergeSublists - Treat 'source' as a sequence of pairs of sorted
         sublists of (at most) 'length' items and merge their nodes (by moving
         links only) onto the list 'result'.
INPUTS:  VoidList* source  The list to split and merge.
         Integer   length  The length of the sublists of source.
         VoidList* result  The empty list.
OUTPUTS: VoidList* source  The empty list.
         VoidList* result  The list of merged nodes extracted from source.
******************************************************************************/

static void mergeSublists(VoidList* source, VoidList* result, Integer length) {
  PRE08( source, source->invariant( source ),
         source->count( source ) >= 2,
         source->comparer( source ),
         result, result->invariant( result ),
         result->count( result ) == 0,
         IN_RANGE( length, 1, source->count( source ) ) );
  CHECKING( const Integer OLD(sourceCount) = source->p->count; )
  VoidListPrivate* const source_p = source->p;
  VoidComparer comparer = source->p->comparer;

  do {
    VoidListNode* node1 = source_p->head;
    VoidListNode* node2 = indexedNode( source, length );
    Integer count1 = 0, count2 = 0;

    do {
      CHECK7( node1, node2, node1 != node2, node1 == source_p->head,
              node1->item, node2->item, node1->item != node2->item );

      if ( comparer( node1->item, node2->item ) <= 0 ) {
        VoidListNode* const nodeToMove = node1;
        node1 = node1->next; /* Advance list1 and append previous node. */
        moveNode( source, nodeToMove, result );
        ++count1; /* Increment count of nodes removed from the first sublist.*/
      } else {
        VoidListNode* const nodeToMove = node2;
        node2 = node2->next; /* Advance list2 and append previous node. */
        moveNode( source, nodeToMove, result );

        if ( node2 == 0 ) {
          count2 = length; /* Reached end of (possibly short) 2nd sublist. */
        } else {
          ++count2; /* Increment count of nodes moved from the 2nd sublist. */
        }
      }
    } while ( AND2( count1 < length, count2 < length ) );

    /* Append any remaining nodes from source sublist1: */

    for ( ; count1 < length; ++count1 ) {
      moveNode( source, source_p->head, result );
    }

    /* Append any remaining nodes from source sublist2: */

    for ( ; AND2( count2 < length, length < source_p->count ); ++count2 ) {
      moveNode( source, source_p->head, result );
    }
  } while ( length < source_p->count );

  /* Append any remainder of the list (no 2nd sublist to compare): */

  while ( source_p->head ) {
    moveNode( source, source_p->head, result );
  }

  POST04( source->invariant( source ),
          source->count( source ) == 0,
          result->invariant( result ),
          result->count( result ) == OLD(sourceCount) );
}



#ifdef DEBUGGING

/******************************************************************************
PURPOSE: print - Prints a VoidList to stderr for debugging purposes.
INPUTS:  VoidList* self  The list to print.
******************************************************************************/

static void print( const VoidList* self ) {
  PRE0( self );
  const VoidListPrivate* const p = self->p;

  fprintf( stderr, "\nself      = %p\n", (void*) self );
  fprintf( stderr, "  free      = %p == %p\n",
           (void*) self->free, (void*) free__ );
  fprintf( stderr, "  insert    = %p == %p\n",
           (void*) self->insert, (void*) insert  );
  fprintf( stderr, "  remove    = %p == %p\n",
           (void*) self->remove, (void*)remove__ );
  fprintf( stderr, "  removeAll = %p == %p\n",
           (void*)self->removeAll, (void*) removeAll );
  fprintf( stderr, "  replace   = %p == %p\n",
           (void*) self->replace, (void*) replace );
  fprintf( stderr, "  sort      = %p == %p\n",
           (void*) self->sort, (void*) sort );
  fprintf( stderr, "  apply     = %p == %p\n",
           (void*) self->apply, (void*) apply );
  fprintf( stderr, "  invariant = %p == %p\n",
           (void*)self->invariant, (void*) invariant );
  fprintf( stderr, "  ok        = %p == %p\n",
           (void*) self->ok, (void*) ok );
  fprintf( stderr, "  sorted    = %p == %p\n",
           (void*) self->sorted, (void*) sorted );
  fprintf( stderr, "  has       = %p == %p\n",
           (void*) self->has, (void*) has );
  fprintf( stderr, "  equal     = %p == %p\n",
           (void*) self->equal, (void*) equal );
  fprintf( stderr, "  count     = %p == %p\n",
           (void*) self->count, (void*) count );
  fprintf( stderr, "  index     = %p == %p\n",
           (void*) self->index, (void*) index__ );
  fprintf( stderr, "  item      = %p == %p\n",
           (void*) self->item, (void*) item );
  fprintf( stderr, "  comparer  = %p == %p\n",
           (void*) self->comparer, (void*) comparer );
  fprintf( stderr, "  deleter   = %p == %p\n",
           (void*) self->deleter, (void*) deleter );

  if ( p ) {
    const VoidListNode* node         = p->head;
    const VoidListNode* previousNode = 0;
    Integer count = 0;
    fprintf( stderr, "  p            = %p\n", (void*) p           );
    fprintf( stderr, "    deleter    = %p\n", (void*) p->deleter  );
    fprintf( stderr, "    comparer   = %p\n", (void*) p->comparer );
    fprintf( stderr, "    head       = %p\n", (void*) p->head     );
    fprintf( stderr, "    tail       = %p\n", (void*) p->tail     );
    fprintf( stderr, "    cache      = %p\n", (void*) p->cache    );
    fprintf( stderr, "    cacheIndex = %"INTEGER_FORMAT"\n", p->cacheIndex );
    fprintf( stderr, "    count      = %"INTEGER_FORMAT"\n", p->count      );
    fprintf( stderr, "    ok         = %"INTEGER_FORMAT"\n", p->ok         );

    while ( node ) {

      if ( node->previous ) {
        fprintf( stderr, "<-" );
      }

      fprintf( stderr, "%p", (void*) node );

      if ( node == p->head ) {
        fprintf( stderr, "(h)" );
      }

      if ( node == p->tail ) {
        fprintf( stderr, "(t)" );
      }

      if ( node == p->cache ) {
        fprintf( stderr, "(c)" );
      }

      if ( node->next ) {
        fprintf( stderr, "->" );
      }

      if ( AND4( node->item     != 0,
                 node->next     != node,
                 node->previous != node,
                 node->previous == previousNode ) ) {
        previousNode = node;
        node         = node->next;
        ++count;

        if ( count > p->count ) {
          fprintf( stderr, "\nDEFECT: cycle in list!\n" );
          node = 0;
        }
      } else {
        fprintf( stderr, "\nDEFECT: incorrectly linked node!\n" );
        node = 0;
      }
    }
  }

  fprintf( stderr, "\n\n" );
}

#endif /* DEBUGGING */





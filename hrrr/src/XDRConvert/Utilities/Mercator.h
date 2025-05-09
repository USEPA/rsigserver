
#ifndef MERCATOR_H
#define MERCATOR_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Mercator.h - Declare Mercator Conformal Conic Projectors ADT.

NOTES:   Commands:

           Mercator* newMercator( Real newMajorSemiaxis, Real newMinorSemiaxis,
                                  Real newCentralLongitude,
                                  Real newFalseEasting, Real newFalseNorthing);

           void free( Mercator* self );

           void setEllipsoid( Mercator* self,
                              Real majorSemiaxis, Real minorSemiaxis );

           void setFalseEasting( Mercator* self, Real falseEasting );

           void setFalseNorthing( Mercator* self, Real falseNorthing );

           void project( const Mercator* self, Real longitude, Real latitude,
                         Real* x, Real* y );

           void unproject( const Mercator* self, Real x, Real y,
                           Real* longitude, Real* latitude );

         Queries:

           Integer invariant( const Mercator* self );
           Integer equal( const Mercator* self, const Mercator* other );
           Mercator* clone( const Mercator* self );

           void ellipsoid( const Mercator* self,
                           Real* majorSemiaxis, Real* minorSemiaxis );

           Real falseEasting( const Mercator* self );
           Real falseNorthing( const Mercator* self );
           Real centralLongitude( const Mercator* self );
           const char* name( const Lambert* self );


         Example:

           #include <stdio.h>
           #include <Mercator.h>

           int main( void ) {
             Projector* projector = (Projector*)
               newMercator( 6370000.0, 6370000.0, -100.0, 0.0, 0.0 );

             if ( projector ) {
               Real x = 0.0, y = 0.0, longitude = 0.0, latitude = 0.0;
               projector->project( projector, -78.7268, 35.9611, &x, &y );
               printf( "%"REAL_F_FORMAT", %"REAL_F_FORMAT"\n", x, y );
               projector->unproject( projector, x, y, &longitude, &latitude );
               printf( "%"REAL_F_FORMAT", %"REAL_F_FORMAT"\n",
                       longitude, latitude );
               FREE_OBJECT( projector );
             }

             return 0;
           }

           cc -64 -mips4 -ansi -fullwarn -g -o example example.c \
              -I/usr/include -I../../../../../include \
              -L../../../../../lib/IRIX64 \
              -lMercator.debug -lProjector.debug \
              -lMemory.debug -lFailure.debug -lBasicNumerics.debug -lm

           example
           1852180.851293, -189978.517654
           -78.726800, 35.961100

HISTORY: 2004/10 Todd Plessel EPA/LM Created.

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Projector.h> /* For DECLARE_PROJECTOR_MEMBERS. */

/*================================== TYPES ==================================*/

typedef struct Mercator Mercator;

/* Constructor: */

extern Mercator* newMercator( Real newMajorSemiaxis,    Real newMinorSemiaxis,
                              Real newCentralLongitude,
                              Real newFalseEasting,     Real newFalseNorthing);

/* Destructor: use macro FREE_OBJECT( yourProjector ); */

/*= PRIVATE =*/

typedef struct MercatorPrivate MercatorPrivate;

#define DECLARE_MERCATOR_MEMBERS( Type ) \
  DECLARE_PROJECTOR_MEMBERS( Type ); \
  MercatorPrivate* data

struct Mercator {
  DECLARE_MERCATOR_MEMBERS( Mercator );
};


#ifdef __cplusplus
}
#endif

#endif /* MERCATOR_H */





#ifndef LAMBERT_H
#define LAMBERT_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Lambert.h - Declare Lambert Conformal Conic Projectors ADT.

NOTES:   Commands:

           Lambert* newLambert( Real newMajorSemiaxis, Real newMinorSemiaxis,
                                Real newLowerLatitude, Real newUpperLatitude,
                                Real newCentralLongitude,
                                Real newCentralLatitude,
                                Real newFalseEasting, Real newFalseNorthing );

           void free( Lambert* self );

           void setEllipsoid( Lambert* self,
                              Real majorSemiaxis, Real minorSemiaxis );

           void setFalseEasting( Lambert* self, Real falseEasting );

           void setFalseNorthing( Lambert* self, Real falseNorthing );

           void project( const Lambert* self, Real longitude, Real latitude,
                         Real* x, Real* y );

           void unproject( const Lambert* self, Real x, Real y,
                           Real* longitude, Real* latitude );

         Queries:

           Integer invariant( const Lambert* self );
           Integer equal( const Lambert* self, const Lambert* other );
           Lambert* clone( const Lambert* self );

           void ellipsoid( const Lambert* self,
                           Real* majorSemiaxis, Real* minorSemiaxis );

           Real falseEasting( const Lambert* self );
           Real falseNorthing( const Lambert* self );
           Real lowerLatitude( const Lambert* self );
           Real upperLatitude( const Lambert* self );
           Real centralLongitude( const Lambert* self );
           Real centralLatitude( const Lambert* self );
           const char* name( const Lambert* self );


         Example:

           #include <stdio.h>
           #include <Lambert.h>

           int main( void ) {
             Projector* projector = (Projector*)
               newLambert( 6370000.0, 6370000.0,
                           30.0, 60.0, -100.0, 40.0, 0.0, 0.0 );

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
              -lLambert.debug -lProjector.debug \
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

typedef struct Lambert Lambert;

/* Constructor: */

extern Lambert* newLambert( Real newMajorSemiaxis,    Real newMinorSemiaxis,
                            Real newLowerLatitude,    Real newUpperLatitude,
                            Real newCentralLongitude, Real newCentralLatitude,
                            Real newFalseEasting,     Real newFalseNorthing );

/* Destructor: use macro FREE_OBJECT( yourProjector ); */

/*= PRIVATE =*/

typedef struct LambertPrivate LambertPrivate;

#define DECLARE_LAMBERT_MEMBERS( Type ) \
  DECLARE_PROJECTOR_MEMBERS( Type ); \
  Real (*lowerLatitude)( const Type* self ); \
  Real (*upperLatitude)( const Type* self ); \
  LambertPrivate* data

struct Lambert {
  DECLARE_LAMBERT_MEMBERS( Lambert );
};


#ifdef __cplusplus
}
#endif

#endif /* LAMBERT_H */




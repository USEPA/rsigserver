
#ifndef MERCATOR_H
#define MERCATOR_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Mercator.h - Declare Mercator Projectors ADT.

NOTES:   Commands:

           Mercator* newMercator( double newMajorSemiaxis,
                                  double newMinorSemiaxis,
                                  double newCentralLongitude,
                                  double newFalseEasting,
                                  double newFalseNorthing );

           void free( Mercator* self );

           void setEllipsoid( Mercator* self,
                              double majorSemiaxis, double minorSemiaxis );

           void setFalseEasting( Mercator* self, double falseEasting );

           void setFalseNorthing( Mercator* self, double falseNorthing );

           void project( const Mercator* self, double longitude, double latitude,
                         double* x, double* y );

           void unproject( const Mercator* self, double x, double y,
                           double* longitude, double* latitude );

         Queries:

           int invariant( const Mercator* self );
           int equal( const Mercator* self, const Mercator* other );
           Mercator* clone( const Mercator* self );

           void ellipsoid( const Mercator* self,
                           double* majorSemiaxis, double* minorSemiaxis );

           double falseEasting( const Mercator* self );
           double falseNorthing( const Mercator* self );
           double centralLongitude( const Mercator* self );
           const char* name( const Lambert* self );


         Example:

           #include <stdio.h>
           #include <Mercator.h>

           int main( void ) {
             Projector* projector = (Projector*)
               newMercator( 6370000.0, 6370000.0, -100.0, 0.0, 0.0 );

             if ( projector ) {
               double x = 0.0, y = 0.0, longitude = 0.0, latitude = 0.0;
               projector->project( projector, -78.7268, 35.9611, &x, &y );
               printf( "%f, %f\n", x, y );
               projector->unproject( projector, x, y, &longitude, &latitude );
               printf( "%f, %f\n",
                       longitude, latitude );
               FREE_OBJECT( projector );
             }

             return 0;
           }

           cc -o example example.c Mercator.c Projector.c -lm

           example
           1852180.851293, -189978.517654
           -78.726800, 35.961100

HISTORY: 2004-10-01 plessel.todd@epa.gov

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Projector.h> /* For macro DECLARE_PROJECTOR_MEMBERS(). */

/*================================== TYPES ==================================*/

typedef struct Mercator Mercator;

/* Constructor: */

extern Mercator* newMercator( double newMajorSemiaxis,    double newMinorSemiaxis,
                              double newCentralLongitude,
                              double newFalseEasting,     double newFalseNorthing);

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




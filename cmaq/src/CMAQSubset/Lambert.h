
#ifndef LAMBERT_H
#define LAMBERT_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Lambert.h - Declare Lambert Conformal Conic Projectors ADT.

NOTES:   Commands:

           Lambert* newLambert( double newMajorSemiaxis,
                                double newMinorSemiaxis,
                                double newLowerLatitude,
                                double newUpperLatitude,
                                double newCentralLongitude,
                                double newCentralLatitude,
                                double newFalseEasting,
                                double newFalseNorthing );

           void free( Lambert* self );

           void setEllipsoid( Lambert* self,
                              double majorSemiaxis, double minorSemiaxis );

           void setFalseEasting( Lambert* self, double falseEasting );

           void setFalseNorthing( Lambert* self, double falseNorthing );

           void project( const Lambert* self, double longitude, double latitude,
                         double* x, double* y );

           void unproject( const Lambert* self, double x, double y,
                           double* longitude, double* latitude );

         Queries:

           int invariant( const Lambert* self );
           int equal( const Lambert* self, const Lambert* other );
           Lambert* clone( const Lambert* self );

           void ellipsoid( const Lambert* self,
                           double* majorSemiaxis, double* minorSemiaxis );

           double falseEasting( const Lambert* self );
           double falseNorthing( const Lambert* self );
           double lowerLatitude( const Lambert* self );
           double upperLatitude( const Lambert* self );
           double centralLongitude( const Lambert* self );
           double centralLatitude( const Lambert* self );
           const char* name( const Lambert* self );


         Example:

           #include <stdio.h>
           #include <Lambert.h>

           int main( void ) {
             Projector* projector = (Projector*)
               newLambert( 6370000.0, 6370000.0,
                           30.0, 60.0, -100.0, 40.0, 0.0, 0.0 );

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

           cc -o example example.c Lambert.c Projector.c -lm

           example
           1852180.851293, -189978.517654
           -78.726800, 35.961100

HISTORY: 2004-10-01 plessel.todd@epa.gov

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Projector.h> /* For macro DECLARE_PROJECTOR_MEMBERS(). */

/*================================== TYPES ==================================*/

typedef struct Lambert Lambert;

/* Constructor: */

extern Lambert* newLambert( double newMajorSemiaxis,    double newMinorSemiaxis,
                            double newLowerLatitude,    double newUpperLatitude,
                            double newCentralLongitude, double newCentralLatitude,
                            double newFalseEasting,     double newFalseNorthing );

/* Destructor: use macro FREE_OBJECT( yourProjector ); */

/*= PRIVATE =*/

typedef struct LambertPrivate LambertPrivate;

#define DECLARE_LAMBERT_MEMBERS( Type ) \
  DECLARE_PROJECTOR_MEMBERS( Type ); \
  double (*lowerLatitude)( const Type* self ); \
  double (*upperLatitude)( const Type* self ); \
  LambertPrivate* data

struct Lambert {
  DECLARE_LAMBERT_MEMBERS( Lambert );
};


#ifdef __cplusplus
}
#endif

#endif /* LAMBERT_H */




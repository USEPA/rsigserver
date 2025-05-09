
#ifndef ALBERS_H
#define ALBERS_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Albers.h - Declare Albers Conformal Conic Projectors ADT.

NOTES:   Commands:

           Albers* newAlbers( double newMajorSemiaxis, double newMinorSemiaxis,
                              double newLowerLatitude, double newUpperLatitude,
                              double newCentralLongitude,
                              double newCentralLatitude,
                              double newFalseEasting, double newFalseNorthing);

           void free( Albers* self );

           void setEllipsoid( Albers* self,
                              double majorSemiaxis, double minorSemiaxis );

           void setFalseEasting( Albers* self, double falseEasting );

           void setFalseNorthing( Albers* self, double falseNorthing );

           void project( const Albers* self, double longitude, double latitude,
                         double* x, double* y );

           void unproject( const Albers* self, double x, double y,
                           double* longitude, double* latitude );

         Queries:

           int invariant( const Albers* self );
           int equal( const Albers* self, const Albers* other );
           Albers* clone( const Albers* self );

           void ellipsoid( const Albers* self,
                           double* majorSemiaxis, double* minorSemiaxis );

           double falseEasting( const Albers* self );
           double falseNorthing( const Albers* self );
           double lowerLatitude( const Albers* self );
           double upperLatitude( const Albers* self );
           double centralLongitude( const Albers* self );
           double centralLatitude( const Albers* self );
           const char* name( const Albers* self );


         Example:

           #include <stdio.h>
           #include <Albers.h>

           int main( void ) {
             Projector* projector = (Projector*)
               newAlbers( 6370000.0, 6370000.0,
                          30.0, 60.0, -100.0, 40.0, 0.0, 0.0 );

             if ( projector ) {
               double x = 0.0, y = 0.0, longitude = 0.0, latitude = 0.0;
               projector->project( projector, -78.7268, 35.9611, &x, &y );
               printf( "%f, %f\n", x, y );
               projector->unproject( projector, x, y, &longitude, &latitude );
               printf( "%f, %f\n", longitude, latitude );
               FREE_OBJECT( projector );
             }

             return 0;
           }

           cc -o example example.c Albers.c Projector.c -lm

           example
           1852180.851293, -189978.517654
           -78.726800, 35.961100

HISTORY: 2004-10-01 plessel.todd@epa.gov

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Projector.h> /* For macro DECLARE_PROJECTOR_MEMBERS(). */

/*================================== TYPES ==================================*/

typedef struct Albers Albers;

/* Constructor: */

extern Albers* newAlbers( double newMajorSemiaxis,   double newMinorSemiaxis,
                          double newLowerLatitude,   double newUpperLatitude,
                          double newCentralLongitude,double newCentralLatitude,
                          double newFalseEasting,    double newFalseNorthing);

/* Destructor: use macro FREE_OBJECT( yourProjector ); */

/*= PRIVATE =*/

typedef struct AlbersPrivate AlbersPrivate;

#define DECLARE_ALBERS_MEMBERS( Type ) \
  DECLARE_PROJECTOR_MEMBERS( Type ); \
  double (*lowerLatitude)( const Type* self ); \
  double (*upperLatitude)( const Type* self ); \
  AlbersPrivate* data

struct Albers {
  DECLARE_ALBERS_MEMBERS( Albers );
};


#ifdef __cplusplus
}
#endif

#endif /* ALBERS_H */




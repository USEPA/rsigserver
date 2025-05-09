
#ifndef STEREOGRAPHIC_H
#define STEREOGRAPHIC_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Stereographic.h - Declare Stereographic Projectors ADT.

NOTES:   Commands:

           Stereographic* newStereographic( double newMajorSemiaxis,
                                            double newMinorSemiaxis,
                                            double newCentralLongitude,
                                            double newCentralLatitude,
                                            double newSecantLatitude,
                                            double newFalseEasting,
                                            double newFalseNorthing );

           void free( Stereographic* self );

           void setEllipsoid( Stereographic* self,
                              double majorSemiaxis, double minorSemiaxis );

           void setFalseEasting( Stereographic* self, double falseEasting );

           void setFalseNorthing( Stereographic* self, double falseNorthing );

           void project( const Stereographic* self,
                         double longitude, double latitude,
                         double* x, double* y );

           void unproject( const Stereographic* self, double x, double y,
                           double* longitude, double* latitude );

         Queries:

           int invariant( const Stereographic* self );
           int equal( const Stereographic* self,
                          const Stereographic* other );
           Stereographic* clone( const Stereographic* self );

           void ellipsoid( const Stereographic* self,
                           double* majorSemiaxis, double* minorSemiaxis );

           double falseEasting( const Stereographic* self );
           double falseNorthing( const Stereographic* self );
           double centralLongitude( const Stereographic* self );
           double centralLatitude( const Stereographic* self );
           double secantLatitude( const Stereographic* self );
           double centralLongitude( const Stereographic* self );
           double centralLatitude( const Stereographic* self );
           const char* name( const Stereographic* self );

         Example:

           #include <stdio.h>
           #include <Stereographic.h>

           int main( void ) {
             Projector* projector =
             (Projector*) newStereographic( 6370997.0, 6370997.0,
                                            -98.0, 90.0, 45.0, 0.0, 0.0 );

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

           cc -o example example.c Stereographic.c Projector.c -lm

HISTORY: 2004-10-01 plessel.todd@epa.gov

STATUS:  unreviewed.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Projector.h> /* For macro DECLARE_PROJECTOR_MEMBERS(). */

/*================================== TYPES ==================================*/

typedef struct Stereographic Stereographic;

/* Constructor: */

extern Stereographic* newStereographic( double newMajorSemiaxis,
                                        double newMinorSemiaxis,
                                        double newCentralLongitude,
                                        double newCentralLatitude,
                                        double newSecantLatitude,
                                        double newFalseEasting,
                                        double newFalseNorthing );

/* Destructor: use macro FREE_OBJECT( yourProjector ); */

/*= PRIVATE =*/

typedef struct StereographicPrivate StereographicPrivate;

#define DECLARE_STEREOGRAPHIC_MEMBERS( Type ) \
  DECLARE_PROJECTOR_MEMBERS( Type ); \
  double (*secantLatitude)( const Type* self ); \
  StereographicPrivate* data

struct Stereographic {
  DECLARE_STEREOGRAPHIC_MEMBERS( Stereographic );
};


#ifdef __cplusplus
}
#endif

#endif /* STEREOGRAPHIC_H */





#ifndef STEREOGRAPHIC_H
#define STEREOGRAPHIC_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Stereographic.h - Declare Stereographic Projectors ADT.

NOTES:   Commands:

           Stereographic* newStereographic( Real newMajorSemiaxis,
                                            Real newMinorSemiaxis,
                                            Real newCentralLongitude,
                                            Real newCentralLatitude,
                                            Real newSecantLatitude,
                                            Real newFalseEasting,
                                            Real newFalseNorthing );

           void free( Stereographic* self );

           void setEllipsoid( Stereographic* self,
                              Real majorSemiaxis, Real minorSemiaxis );

           void setFalseEasting( Stereographic* self, Real falseEasting );

           void setFalseNorthing( Stereographic* self, Real falseNorthing );

           void project( const Stereographic* self,
                         Real longitude, Real latitude,
                         Real* x, Real* y );

           void unproject( const Stereographic* self, Real x, Real y,
                           Real* longitude, Real* latitude );

         Queries:

           Integer invariant( const Stereographic* self );
           Integer equal( const Stereographic* self,
                          const Stereographic* other );
           Stereographic* clone( const Stereographic* self );

           void ellipsoid( const Stereographic* self,
                           Real* majorSemiaxis, Real* minorSemiaxis );

           Real falseEasting( const Stereographic* self );
           Real falseNorthing( const Stereographic* self );
           Real centralLongitude( const Stereographic* self );
           Real centralLatitude( const Stereographic* self );
           Real secantLatitude( const Stereographic* self );
           Real centralLongitude( const Stereographic* self );
           Real centralLatitude( const Stereographic* self );
           const char* name( const Stereographic* self );

         Example:

           #include <stdio.h>
           #include <Stereographic.h>

           int main( void ) {
             Projector* projector =
             (Projector*) newStereographic( 6370997.0, 6370997.0,
                                            -98.0, 90.0, 45.0, 0.0, 0.0 );

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

HISTORY: 2004/10 Todd Plessel  EPA/LM Created.

STATUS:  unreviewed.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Projector.h> /* For DECLARE_PROJECTOR_MEMBERS. */

/*================================== TYPES ==================================*/

typedef struct Stereographic Stereographic;

/* Constructor: */

extern Stereographic* newStereographic( Real newMajorSemiaxis,
                                        Real newMinorSemiaxis,
                                        Real newCentralLongitude,
                                        Real newCentralLatitude,
                                        Real newSecantLatitude,
                                        Real newFalseEasting,
                                        Real newFalseNorthing );

/* Destructor: use macro FREE_OBJECT( yourProjector ); */

/*= PRIVATE =*/

typedef struct StereographicPrivate StereographicPrivate;

#define DECLARE_STEREOGRAPHIC_MEMBERS( Type ) \
  DECLARE_PROJECTOR_MEMBERS( Type ); \
  Real (*secantLatitude)( const Type* self ); \
  StereographicPrivate* data

struct Stereographic {
  DECLARE_STEREOGRAPHIC_MEMBERS( Stereographic );
};


#ifdef __cplusplus
}
#endif

#endif /* STEREOGRAPHIC_H */




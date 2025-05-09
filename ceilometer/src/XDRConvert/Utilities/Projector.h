
#ifndef PROJECTOR_H
#define PROJECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Projector.h - Declare Cartographic projectors ADT ABC.

NOTES:   Commands:

           void free( Projector* self );

           void setEllipsoid( Projector* self,
                              Real majorSemiaxis, Real minorSemiaxis );

           void setFalseEasting( Projector* self, Real falseEasting );

           void setFalseNorthing( Projector* self, Real falseNorthing );

           void project( const Projector* self, Real longitude, Real latitude,
                         Real* x, Real* y );

           void unproject( const Projector* self, Real x, Real y,
                           Real* longitude, Real* latitude );

         Queries:

           Integer invariant( const Projector* self );

           Integer equal( const Projector* self, const Projector* other );

           Projector* clone( const Projector* self );

           void ellipsoid( const Projector* self,
                           Real* majorSemiaxis, Real* minorSemiaxis );

           Real falseEasting( const Projector* self );

           Real falseNorthing( const Projector* self );

           Real centralLongitude( const Projector* self );

           Real centralLatitude( const Projector* self );

           const char* name( const Projector* self );
 
         See Lambert.h for example usage.

HISTORY: 2004/10 Todd Plessel EPA/LM Created.

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <BasicNumerics.h> /* For Real, Integer. */
#include <Memory.h> /* For macro FREE_OBJECT() used by clients. */

/*================================= MACROS ==================================*/

#define PI_OVER_2 1.57079632679489661923
#define PI_OVER_4 0.78539816339744830962
#define PROJECTION_TOLERANCE 1e-10
#define CONVERGENCE_TOLERANCE 1e-12
enum { MAXIMUM_ITERATIONS = 15 };
#define AIRY_1830_MAJOR_SEMIAXIS 6377563.4
#define AIRY_1830_MINOR_SEMIAXIS 6356256.9
#define MODIFIED_AIRY_MAJOR_SEMIAXIS 6377340.2
#define MODIFIED_AIRY_MINOR_SEMIAXIS 6356034.4
#define ANDRAE_1876_MAJOR_SEMIAXIS 6377104.4
#define ANDRAE_1876_MINOR_SEMIAXIS 6355847.4
#define APPLIED_PHYSICS_1965_MAJOR_SEMIAXIS 6378137.0
#define APPLIED_PHYSICS_1965_MINOR_SEMIAXIS 6356751.8
#define AUSTRALIAN_NATL_SA_1969_MAJOR_SEMIAXIS 6378160.0
#define AUSTRALIAN_NATL_SA_1969_MINOR_SEMIAXIS 6356774.7
#define BESSEL_1841_MAJOR_SEMIAXIS 6377397.2
#define BESSEL_1841_MINOR_SEMIAXIS 6356079.0
#define BESSEL_NAMIBIA_1841_MAJOR_SEMIAXIS 6377483.9
#define BESSEL_NAMIBIA_1841_MINOR_SEMIAXIS 6356165.4
#define CLARKE_1866_MAJOR_SEMIAXIS 6378206.4
#define CLARKE_1866_MINOR_SEMIAXIS 6356583.8
#define CLARKE_1880_MAJOR_SEMIAXIS 6378249.1
#define CLARKE_1880_MINOR_SEMIAXIS 6356515.0
#define COMM_DES_POIDS_ET_MESURES_1799_MAJOR_SEMIAXIS 6375738.7
#define COMM_DES_POIDS_ET_MESURES_1799_MINOR_SEMIAXIS 6356666.2
#define DELAMBRE_1810_BELGIUM_MAJOR_SEMIAXIS 6376428.0
#define DELAMBRE_1810_BELGIUM_MINOR_SEMIAXIS 6355957.9
#define ENGELIS_1985_MAJOR_SEMIAXIS 6378136.1
#define ENGELIS_1985_MINOR_SEMIAXIS 6356751.3
#define EVEREST_1830_MAJOR_SEMIAXIS 6377276.3
#define EVEREST_1830_MINOR_SEMIAXIS 6356075.4
#define EVEREST_1948_MAJOR_SEMIAXIS 6377304.1
#define EVEREST_1948_MINOR_SEMIAXIS 6356103.0
#define EVEREST_1956_MAJOR_SEMIAXIS 6377301.2
#define EVEREST_1956_MINOR_SEMIAXIS 6356100.2
#define EVEREST_1969_MAJOR_SEMIAXIS 6377295.7
#define EVEREST_1969_MINOR_SEMIAXIS 6356094.7
#define EVEREST_SABAH_SARAWAK_MAJOR_SEMIAXIS 6377298.6
#define EVEREST_SABAH_SARAWAK_MINOR_SEMIAXIS 6356097.6
#define FISCHER_MERCURY_DATUM_1960_MAJOR_SEMIAXIS 6378166.0
#define FISCHER_MERCURY_DATUM_1960_MINOR_SEMIAXIS 6356784.3
#define MODIFIED_FISCHER_1960_MAJOR_SEMIAXIS 6378155.0
#define MODIFIED_FISCHER_1960_MINOR_SEMIAXIS 6356773.3
#define FISCHER_1968_MAJOR_SEMIAXIS 6378150.0
#define FISCHER_1968_MINOR_SEMIAXIS 6356768.3
#define GRS_IUGG_1967_MAJOR_SEMIAXIS 6378160.0
#define GRS_IUGG_1967_MINOR_SEMIAXIS 6352363.3
#define GRS_IUGG_1980_MAJOR_SEMIAXIS 6378137.0
#define GRS_IUGG_1980_MINOR_SEMIAXIS 6356752.3
#define HELMERT_1906_MAJOR_SEMIAXIS 6378200.0
#define HELMERT_1906_MINOR_SEMIAXIS 6356818.2
#define HOUGH_MAJOR_SEMIAXIS 6378270.0
#define HOUGH_MINOR_SEMIAXIS 6356794.3
#define IAU_1976_MAJOR_SEMIAXIS 6378140.0
#define IAU_1976_MINOR_SEMIAXIS 6356755.3
#define INTL_HAYFORD_1909_MAJOR_SEMIAXIS 6378388.0
#define INTL_HAYFORD_1909_MINOR_SEMIAXIS 6356911.9
#define KRASSOVSKY_1942_MAJOR_SEMIAXIS 6378245.0
#define KRASSOVSKY_1942_MINOR_SEMIAXIS 6356863.0
#define KAULA_1961_MAJOR_SEMIAXIS 6378163.0
#define KAULA_1961_MINOR_SEMIAXIS 6356777.0
#define LERCH_1979_MAJOR_SEMIAXIS 6378139.0
#define LERCH_1979_MINOR_SEMIAXIS 6356754.3
#define MAUPERTIUS_1738_MAJOR_SEMIAXIS 6397300.0
#define MAUPERTIUS_1738_MINOR_SEMIAXIS 6363806.3
#define MERIT_1983_MAJOR_SEMIAXIS 6378137.0
#define MERIT_1983_MINOR_SEMIAXIS 6356752.3
#define NAVAL_WEAPONS_LAB_1965_MAJOR_SEMIAXIS 6378145.0
#define NAVAL_WEAPONS_LAB_1965_MINOR_SEMIAXIS 6356759.8
#define NEW_INTERNATIONAL_1967_MAJOR_SEMIAXIS 6378157.5
#define NEW_INTERNATIONAL_1967_MINOR_SEMIAXIS 6356772.2
#define PLESSIS_1817_MAJOR_SEMIAXIS 6376523.0
#define PLESSIS_1817_MINOR_SEMIAXIS 6355863.0
#define SGS_1985_MAJOR_SEMIAXIS 6378136.0
#define SGS_1985_MINOR_SEMIAXIS 6356751.3
#define SOUTHEAST_ASIA_MAJOR_SEMIAXIS 6378155.0
#define SOUTHEAST_ASIA_MINOR_SEMIAXIS 6356773.0
#define WALBECK_MAJOR_SEMIAXIS 6376896.0
#define WALBECK_MINOR_SEMIAXIS 6355835.0
#define WGS_1960_MAJOR_SEMIAXIS 6378165.0
#define WGS_1960_MINOR_SEMIAXIS 6356783.3
#define WGS_1966_MAJOR_SEMIAXIS 6378145.0
#define WGS_1966_MINOR_SEMIAXIS 6356759.8
#define WGS_1972_MAJOR_SEMIAXIS 6378135.0
#define WGS_1972_MINOR_SEMIAXIS 6356750.5
#define WGS_1984_MAJOR_SEMIAXIS 6378137.0
#define WGS_1984_MINOR_SEMIAXIS 6356752.3
#define MM5_RADIUS 6370997.0
#define MCIDAS_RADIUS 6371230.0
#define MOON_RADIUS 1738000.0
#define MARS_MAJOR_SEMIAXIS 3394500.0
#define MARS_MINOR_SEMIAXIS 3376400.0
#define VENUS_RADIUS 6051000.0

#define SQUARE(x) ((x)*(x))

/*================================== TYPES ==================================*/

typedef struct Projector Projector;

/*= PRIVATE =*/

#define DECLARE_PROJECTOR_MEMBERS( Type ) \
  void (*free)( Type* self ); \
  void (*setEllipsoid)( Type* self, Real majorSemiaxis, Real minorSemiaxis); \
  void (*setFalseEasting)( Type* self, Real falseEasting ); \
  void (*setFalseNorthing)( Type* self, Real falseNorthing ); \
  void (*project)( const Type* self, Real longitude, Real latitude, \
                   Real* x,Real* y ); \
  void (*unproject)( const Type* self, Real x, Real y, \
                     Real* longitude, Real* latitude ); \
  Integer (*invariant)( const Type* self ); \
  Integer (*equal)( const Type* self, const Type* other ); \
  Type* (*clone)( const Type* self ); \
  void (*ellipsoid)( const Type* self, \
                     Real* majorSemiaxis, Real* minorSemiaxis ); \
  Real (*falseEasting)( const Type* self ); \
  Real (*falseNorthing)( const Type* self ); \
  Real (*centralLongitude)( const Type* self ); \
  Real (*centralLatitude)( const Type* self ); \
  const char* (*name)( const Type* self );

struct Projector {
  DECLARE_PROJECTOR_MEMBERS( Projector );
};



/*================================ FUNCTIONS ================================*/

extern Integer isValidEllipsoid( Real majorSemiaxis, Real minorSemiaxis );

extern Integer isValidLongitude( Real longitude );
extern Integer isValidLatitude( Real latitude );
extern Integer isValidLongitudeLatitude( Real longitude, Real latitude );

extern Integer validLongitudesAndLatitudes( Integer count,
                                            const Real longitudes[],
                                            const Real latitudes[] );

extern Real latitudeWGS84( Real latitudeSphere );
extern Real latitudeSphere( Real latitudeWGS84 );

extern Real msfn( Real sinePhi, Real cosinePhi,
                  Real eccentricitySquared );

extern Real tsfn( Real phi, Real sinePhi, Real ellipsoidEccentricity );
extern Real ssfn( Real phi, Real sinePhi, Real ellipsoidEccentricity );
extern Real phi2Iterate( Real ts, Real theEccentricity );


#ifdef __cplusplus
}
#endif

#endif /* PROJECTOR_H */



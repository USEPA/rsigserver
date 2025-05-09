
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
                              double majorSemiaxis, double minorSemiaxis );

           void setFalseEasting( Projector* self, double falseEasting );

           void setFalseNorthing( Projector* self, double falseNorthing );

           void project( const Projector* self, double longitude, double latitude,
                         double* x, double* y );

           void unproject( const Projector* self, double x, double y,
                           double* longitude, double* latitude );

         Queries:

           int invariant( const Projector* self );

           int equal( const Projector* self, const Projector* other );

           Projector* clone( const Projector* self );

           void ellipsoid( const Projector* self,
                           double* majorSemiaxis, double* minorSemiaxis );

           double falseEasting( const Projector* self );

           double falseNorthing( const Projector* self );

           double centralLongitude( const Projector* self );

           double centralLatitude( const Projector* self );

           const char* name( const Projector* self );
 
         See Lambert.h for example usage.

HISTORY: 2004-10-01 plessel.todd@epa.gov

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h>  /* For size_t. */
#include <stdlib.h>  /* For free(). */
#include <string.h>  /* For memset(). */

/*================================= MACROS ==================================*/

#define PI_OVER_2 1.57079632679489661923
#define PI_OVER_4 0.78539816339744830962
#define TOLERANCE 1e-6
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
#ifndef SIGN
#define SIGN(x) ((x) < 0 ? -1 : 1)
#endif

#define FREE_ZERO( p ) \
( ( (p) ? memset( (p), 0, sizeof *(p) ), free(p) : (void) 0 ), (p) = 0 )

#define FREE_OBJECT( p ) (((p) ? (p)->free(p) : (void) 0), (p) = 0)

/*================================== TYPES ==================================*/

typedef struct Projector Projector;

/*= PRIVATE =*/

#define DECLARE_PROJECTOR_MEMBERS( Type ) \
  void (*free)( Type* self ); \
  void (*setEllipsoid)( Type* self, double majorSemiaxis, double minorSemiaxis); \
  void (*setFalseEasting)( Type* self, double falseEasting ); \
  void (*setFalseNorthing)( Type* self, double falseNorthing ); \
  void (*project)( const Type* self, double longitude, double latitude, \
                   double* x,double* y ); \
  void (*unproject)( const Type* self, double x, double y, \
                     double* longitude, double* latitude ); \
  int (*invariant)( const Type* self ); \
  int (*equal)( const Type* self, const Type* other ); \
  Type* (*clone)( const Type* self ); \
  void (*ellipsoid)( const Type* self, \
                     double* majorSemiaxis, double* minorSemiaxis ); \
  double (*falseEasting)( const Type* self ); \
  double (*falseNorthing)( const Type* self ); \
  double (*centralLongitude)( const Type* self ); \
  double (*centralLatitude)( const Type* self ); \
  const char* (*name)( const Type* self );

struct Projector {
  DECLARE_PROJECTOR_MEMBERS( Projector );
};



/*================================ FUNCTIONS ================================*/

extern int isNan( double x );
extern double safeDifference( double x, double y );
extern double safeQuotient( double numerator, double denominator );
extern int withinTolerance( double x, double y, double tolerance );
extern int aboutEqual( double x, double y );

extern double radians( double theDegrees );
extern double degrees( double theRadians );

extern int isValidEllipsoid( double majorSemiaxis, double minorSemiaxis );

extern int isValidLongitude( double longitude );
extern int isValidLatitude( double latitude );
extern int isValidLongitudeLatitude( double longitude, double latitude );

extern int validLongitudesAndLatitudes( size_t count,
                                        const double longitudes[],
                                        const double latitudes[] );

extern double msfn( double sinePhi, double cosinePhi,
                    double eccentricitySquared );

extern double qsfn( double sinePhi, double ellipsoidEccentricity,
                    double oneMinusEllipsoidEccentricitySquared );

extern double tsfn( double phi, double sinePhi, double ellipsoidEccentricity );

extern double ssfn( double phi, double sinePhi, double ellipsoidEccentricity );

extern double phi1Iterate( double phi,
                           double the_eccentricity,
                           double the_one_minus_eccentricity_squared );

extern double phi2Iterate( double ts, double theEccentricity );


#ifdef __cplusplus
}
#endif

#endif /* PROJECTOR_H */



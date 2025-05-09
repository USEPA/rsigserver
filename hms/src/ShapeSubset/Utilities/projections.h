
#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef TOLERANCE
#define TOLERANCE 1e-6
#endif
#define PROJECTION_TOLERANCE 1e-10
#define CONVERGENCE_TOLERANCE 1e-12
enum { MAXIMUM_ITERATIONS = 15 };

#define PI_OVER_2 1.57079632679489661923
#define PI_OVER_4 0.78539816339744830962
#define SIGN(x) ((x) < 0 ? -1 : 1)
#define SQUARE(x) ((x)*(x))

int is_valid_ellipsoid( double major_semiaxis, double minor_semiaxis );
int is_valid_longitude( double longitude );
int is_valid_latitude( double latitude );
int is_valid_longitude_latitude( double longitude, double latitude );
int is_valid_longitudes_and_latitudes( long long count,
                                       const double longitudes[],
                                       const double latitudes[] );

int is_nan( double value );
double to_radians( double the_degrees );
double to_degrees( double the_radians );
double safe_difference( double left, double right );
double safe_quotient( double numerator, double denominator );
int about_equal( double x, double y );

double ssfn( double phi, double sine_phi, double ellipsoid_eccentricity );
double msfn( double sine_phi, double cosine_phi, double eccentricity_squared );
double tsfn( double phi, double sine_phi, double ellipsoid_eccentricity );
double qsfn( double sine_phi, double ellipsoid_eccentricity,
             double one_minus_ellipsoid_eccentricity_squared );


#ifdef __cplusplus
}
#endif

#endif /* PROJECTIONS_H */


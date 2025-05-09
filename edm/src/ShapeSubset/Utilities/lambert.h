
#ifndef LAMBERT_H
#define LAMBERT_H

#ifdef __cplusplus
extern "C" {
#endif

/* lambert.h - declare C-callable routines for Lambert unprojecting. */

void initialize_lambert( double new_major_semiaxis,
                         double new_minor_semiaxis,
                         double new_lower_latitude,
                         double new_upper_latitude,
                         double new_central_latitude,
                         double new_central_longitude,
                         double new_false_easting,
                         double new_false_northing );

void project_lambert( double longitude, double latitude, double* x, double* y );

void unproject_lambert( double x, double y,
                        double* longitude, double* latitude );

void lambert_center( double* central_longitude, double* central_latitude );

void lambert_tangents( double* lower_latitude, double* upper_latitude );

#ifdef __cplusplus
}
#endif

#endif /* LAMBERT_H */



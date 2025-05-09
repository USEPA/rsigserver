
#ifndef ALBERS_H
#define ALBERS_H

#ifdef __cplusplus
extern "C" {
#endif

/* albers.h - declare C-callable routines for Albers unprojecting. */

void initialize_albers( double new_major_semiaxis,
                        double new_minor_semiaxis,
                        double new_lower_latitude,
                        double new_upper_latitude,
                        double new_central_latitude,
                        double new_central_longitude,
                        double new_false_easting,
                        double new_false_northing );

void project_albers( double longitude, double latitude, double* x, double* y );

void unproject_albers( double x, double y,
                       double* longitude, double* latitude );

void albers_center( double* central_longitude, double* central_latitude );

void albers_tangents( double* lower_latitude, double* upper_latitude );

#ifdef __cplusplus
}
#endif

#endif /* ALBERS_H */



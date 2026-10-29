/******************************************************************************
PURPOSE: Projector.c - Define Cartographic projectors ADT ABC.
NOTES:   See Lambert.h for example usage.
HISTORY: 2004/10, Created based on C++ version.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <float.h> /* For DBL_MIN, DBLMAX. */
#include <math.h>  /* For M_PI, fabs(), tan(), atan(). */

#include <Assertions.h>    /* For PRE0*(), POST0*(), IN_RANGE(), IS_BOOL().  */
#include <BasicNumerics.h> /* For Integer, Real, isNan(), withinTolerance(). */
#include <Projector.h>     /* For public interface.                          */

/* Adjust latitude to/from WGS84/Sphere? */

#define ADJUST_LATITUDE 0

/*============================ PUBLIC FUNCTIONS =============================*/


/******************************************************************************
PURPOSE: isValidEllipsoid - Do the arguments define a valid ellipsoid?
INPUTS:  Real majorSemiaxis Mean equitorial radius of planet approximation.
         Real minorSemiaxis Mean polar      radius of planet approximation.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidEllipsoid( Real majorSemiaxis, Real minorSemiaxis ) {
  const Integer result =
    AND7( ! isNan( majorSemiaxis ),
          ! isNan( minorSemiaxis ),
          majorSemiaxis > 0.0,
          minorSemiaxis > 0.0,
          majorSemiaxis >= minorSemiaxis,
          SQUARE( majorSemiaxis ) > 0.0,
          SQUARE( minorSemiaxis ) > 0.0 );
  POST0( IS_BOOL( result ) );
  return result;
}


/******************************************************************************
PURPOSE: isValidLongitude - Is the argument a valid longitude?
INPUTS:  Real longitude In degrees.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidLongitude( Real longitude ) {
  const Integer result =
    AND2( ! isNan( longitude ), IN_RANGE( longitude, -180.0, 180.0 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidLatitude - Is the argument a valid latitude?
INPUTS:  Real latitude In degrees.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidLatitude( Real latitude ) {
  const Integer result =
    AND2( ! isNan( latitude ), IN_RANGE( latitude, -90.0, 90.0 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidLongitudeLatitude - Are the arguments a valid longitude
         latitude point?
INPUTS:  Real longitude In degrees.
         Real latitude  In degrees.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidLongitudeLatitude( Real longitude, Real latitude ) {
  const Integer result =
    AND2( isValidLongitude( longitude ), isValidLatitude( latitude ) );
  POST0( result ==
         AND2( isValidLongitude( longitude ), isValidLatitude( latitude ) ) );
  return result;
}



/******************************************************************************
PURPOSE: validLongitudesAndLatitudes - Are longitudes and latitudes valid?
INPUTS:  Integer count            Number of points.
         const Real longitudes[]  Longitudes to check.
         const Real latitudes[]   Latitudes to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer validLongitudesAndLatitudes( Integer count,
                                     const Real longitudes[],
                                     const Real latitudes[] ) {
  PRE03( count > 0, longitudes, latitudes );
  Integer result = 0;
  Integer index = 0;

  do {

    if ( ! isValidLongitudeLatitude( longitudes[ index ],
                                     latitudes[ index ] ) ) {
      index = count;
    }

    ++index;
  } while ( index < count );

  result = index == count;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: latitudeWGS84 - Convert latitude on a sphere to latitude on a
         WGS84/GRS80/NAD83 spheroid.
INPUTS:  Real latitudeSphere   Latitude on a sphere.
RETURNS: Real equivalent latitude on a WGS84 spheroid.
NOTES:   http://en.wikipedia.org/wiki/Latitude#Geocentric_latitude.
******************************************************************************/

Real latitudeWGS84( Real latitudeSphere ) {
  PRE0( isValidLatitude( latitudeSphere ) );
#if ADJUST_LATITUDE
  const Real inverseWGS84AxisRatioSquared = 1.006739496756587;
  const Real latitudeSphereRadians = radians( latitudeSphere );
  const Real latitudeWGS84Radians =
    atan( tan( latitudeSphereRadians ) * inverseWGS84AxisRatioSquared );
  const Real result = degrees( latitudeWGS84Radians );
#else
  const Real result = latitudeSphere;
#endif
  POST02( isValidLatitude( result ),
          IMPLIES_ELSE( OR3( aboutEqual( latitudeSphere, 0.0 ),
                             aboutEqual( latitudeSphere, -90.0 ),
                             aboutEqual( latitudeSphere, 90.0 ) ),
                        aboutEqual( result, latitudeSphere ),
                        fabs( result - latitudeSphere ) < 0.1925 ) );                   
  return result;
}



/******************************************************************************
PURPOSE: latitudeSphere - Convert latitude on a  WGS84/GRS80/NAD83 spheroid to
         a latitude on a sphere.
INPUTS:  Real latitudeSphere   Latitude on a sphere.
RETURNS: Real equivalent latitude on a WGS84 spheroid.
NOTES:   http://en.wikipedia.org/wiki/Latitude#Geocentric_latitude.
******************************************************************************/

Real latitudeSphere( Real latitudeWGS84 ) {
  PRE0( isValidLatitude( latitudeWGS84 ) );
#if ADJUST_LATITUDE
  const Real WGS84AxisRatioSquared = 0.9933056199957391;
  const Real latitudeWGS84Radians = radians( latitudeWGS84 );
  const Real latitudeSphereRadians =
    atan( tan( latitudeWGS84Radians ) * WGS84AxisRatioSquared );
  const Real result = degrees( latitudeSphereRadians );
#else
  const Real result = latitudeWGS84;
#endif
  POST02( isValidLatitude( result ),
          IMPLIES_ELSE( OR3( aboutEqual( latitudeWGS84, 0.0 ),
                             aboutEqual( latitudeWGS84, -90.0 ),
                             aboutEqual( latitudeWGS84, 90.0 ) ),
                        aboutEqual( result, latitudeWGS84 ),
                        fabs( result - latitudeWGS84 ) < 0.1925 ) );                   
  return result;
}



/******************************************************************************
PURPOSE: ssfn - See USGS PROJ Library.
INPUTS:  Real phi      Angle in radians.
         Real sinePhi  Sine of phi.
         Real ellipsoidEccentricity  Of planet approximation.
RETURNS: Real See USGS PROJ Library.
******************************************************************************/

Real ssfn( Real phi, Real sinePhi, Real ellipsoidEccentricity ) {

  PRE07( ! isNan( phi ),
         ! isNan( sinePhi ),
         ! isNan( ellipsoidEccentricity ),
         withinTolerance( sinePhi, sin( phi ), PROJECTION_TOLERANCE ),
         sinePhi > -1.0,
         sinePhi < 1.0,
         IN_RANGE( ellipsoidEccentricity, 0.0, 1.0 ) );

  const Real eccentricitySinePhi = ellipsoidEccentricity * sinePhi;
  const Real exponent = ellipsoidEccentricity * 0.5;
  const Real factor1 = tan( ( PI_OVER_2 + phi ) * 0.5 );
  const Real factor2 = pow( ( ( 1.0 - eccentricitySinePhi ) /
                              ( 1.0 + eccentricitySinePhi ) ),
                            exponent );
  const Real result = factor1 * factor2;

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: msfn - See USGS PROJ Library.
INPUTS:  Real sinePhi    Sine of phi.
         Real cosinePhi  Cosine of phi.
         Real eccentricitySquared  Of planet approximation.
RETURNS: Real See USGS PROJ Library.
******************************************************************************/

Real msfn(Real sinePhi, Real cosinePhi, Real eccentricitySquared) {

  PRE012( ! isNan( sinePhi ),
          ! isNan( cosinePhi ),
          ! isNan( eccentricitySquared ),
          withinTolerance( sinePhi, sqrt( 1.0 - SQUARE( cosinePhi ) ),
                           PROJECTION_TOLERANCE ),
          sinePhi   > -1.0,
          sinePhi   < 1.0,
          cosinePhi > -1.0,
          cosinePhi < 1.0,
          cosinePhi != 0.0,
          IN_RANGE( eccentricitySquared, 0.0, 1.0 ),
          eccentricitySquared * sinePhi * sinePhi < 1.0,
          sqrt( 1.0 - eccentricitySquared * SQUARE( sinePhi ) ) != 0.0 );

  const Real result = cosinePhi /
               sqrt( 1.0 - eccentricitySquared * SQUARE( sinePhi ) );

  POST02( ! isNan( result ), result != 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: tsfn - See USGS PROJ Library.
INPUTS:  Real phi      Angle in radians.
         Real sinePhi  Sine of phi.
         Real ellipsoidEccentricity  Of planet approximation.
RETURNS: Real See USGS PROJ Library.
******************************************************************************/

Real tsfn( Real phi, Real sinePhi, Real ellipsoidEccentricity ) {

  PRE010( ! isNan( phi ),
          ! isNan( sinePhi ),
          ! isNan( ellipsoidEccentricity ),
          withinTolerance( sinePhi, sin( phi ), PROJECTION_TOLERANCE ),
          sinePhi > -1.0,
          sinePhi < 1.0,
          IN_RANGE( ellipsoidEccentricity, 0.0, 1.0 ),
          tan( ( PI_OVER_2 - phi ) * 0.5 ) != 0.0,
          fabs( ellipsoidEccentricity * sinePhi ) != 1.0,
          ( 1.0 + ellipsoidEccentricity * sinePhi ) != 0.0 );

  const Real eccentricitySinePhi = ellipsoidEccentricity * sinePhi;
  const Real exponent = ellipsoidEccentricity * 0.5;
  const Real numerator = tan( ( PI_OVER_2 - phi ) * 0.5 );
  const Real denominator = pow( ( ( 1.0 - eccentricitySinePhi ) /
                                  ( 1.0 + eccentricitySinePhi ) ),
                                exponent );
  const Real result = numerator / denominator;

  POST02( ! isNan( result ), result != 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: phi2Iterate - Iterate on unprojected y coordinate.
RETURNS: Real converged phi.
******************************************************************************/

Real phi2Iterate( Real ts, Real theEccentricity ) {

  PRE0( IN_RANGE( theEccentricity, 0.0, 1.0 ) );

  const Integer maximumIterations = MAXIMUM_ITERATIONS;
  const Real convergenceTolerance = CONVERGENCE_TOLERANCE;
  const Real halfEccentricity     = theEccentricity * 0.5;
  Real deltaPhi = 0.0;
  Integer iteration = 0;
  Real result = PI_OVER_2 - 2.0 * atan( ts );

  do {
    const Real con = theEccentricity * sin( result );
    CHECK( con != -1.0 );
    deltaPhi =
      PI_OVER_2 -
      2.0 * atan( ts * pow( ( 1.0 - con ) / ( 1.0 + con ), halfEccentricity))
      - result;
    result += deltaPhi;
    ++iteration;
  } while ( AND2( fabs( deltaPhi ) >= convergenceTolerance,
                  iteration < maximumIterations ) );

  POST0( ! isNan( result ) );
  return result;
}





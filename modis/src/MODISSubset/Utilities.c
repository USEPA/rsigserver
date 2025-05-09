/******************************************************************************
PURPOSE: Utilities.c - Some general-purpose reusable routines.

HISTORY: 2017-12-19 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/



/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, stderr, fprintf(). */
#include <stdlib.h>    /* For malloc(), free(), abs(). */
#include <string.h>    /* For strlen(). */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

#include "Utilities.h" /* For public interface. */

/*================================= MACROS ==================================*/

/*
 * Is the platform big-endian (MSB: most significant byte first) or
 * little-endian (LSB: least significant byte first)?
 */

#if defined(__alpha) || \
    defined(__i386__) || defined(__i486__) || \
    defined(__i586__) || defined(__i686__) || \
    defined(__ia64__) || defined(__x86_64__)
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif

/*============================ GLOBAL CONSTANTS =============================*/

static const double EDGE = 179.99; /* Clamp longitude to EDGE if cross +/-180*/

/*
 * 30 days hath September, April, June and November, all the rest have 31,
 * except February which has either 28 or 29 (on a leap year).
 */

static const int daysPerMonth[ 2 ][ 12 ] = {
  { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
  { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
};

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: isValidLongitude - Is the argument a valid longitude?
INPUTS:  const double longitude In degrees.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidLongitude( const double longitude ) {
  const int result = IN_RANGE( longitude, -180.0, 180.0 );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidLatitude - Is the argument a valid latitude?
INPUTS:  const double latitude In degrees.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidLatitude( const double latitude ) {
  const int result = IN_RANGE( latitude, -90.0, 90.0 );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: clampInvalidCoordinates - Clamp invalid longitude-latitude points.
INPUTS:  const size_t points          Points of data.
         double longitudes[ points ]  Longitudes to check.
         double latitudes[  points ]  Latitudes to check.
INPUTS:  double longitudes[ points ]  Clamped Longitudes.
         double latitudes[  points ]  Clamped Latitudes.
RETURNS: int 1 if at least one valid point was found.
******************************************************************************/

int clampInvalidCoordinates(  const size_t points,
                              double longitudes[],
                              double latitudes[] ) {

  int result = 0;

  assert( points );
  assert( longitudes ); assert( latitudes );

  {
    size_t point = 0;
    size_t validPoint = 0; /* Index of first valid point: */

    for ( point = 0; point < points; ++point ) {
      const double longitude = longitudes[ point ];

      if ( IN_RANGE( longitude, -180.0, 180.0 ) ) {
        const double latitude = latitudes[ point ];

        if ( IN_RANGE( latitude, -90.0, 90.0 ) ) {
          result = 1;
          validPoint = point;
          point = points; /* Stop looping. */
        }
      }
    }

    if ( result ) {
      const double longitude = longitudes[ validPoint ];
      const double latitude  = latitudes[  validPoint ];

      /* Clamp all previous points to this first valid point: */

      for ( point = 0; point < validPoint; ++point ) {
        longitudes[ point ] = longitude;
        latitudes[  point ] = latitude;
      }

      /* Clamp all remaining points to the previous valid point: */

      for ( ; point < points; ++point ) {

        if ( IN_RANGE( longitudes[ point ], -180.0, 180.0 ) &&
             IN_RANGE( latitudes[  point ],  -90.0,  90.0 ) ) {
          validPoint = point;
        } else {
          longitudes[ point ] = longitudes[ validPoint ];
          latitudes[  point ] = latitudes[  validPoint ];
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: pointsInSubset - Compute number of points in subset based on domain
         and assign MISSING_VALUE to data values that are outside the subset.
INPUTS:  const Bounds domain                  Domain to subset data to.
         const size_t points                  Points of data.
         const double longitudes[   points ]  Longitudes to check.
         const double latitudes[    points ]  Latitudes to check.
         const double longitudesSW[ points ]  0 or corners to check.
         const double longitudesSE[ points ]  0 or corners to check.
         const double longitudesNW[ points ]  0 or corners to check.
         const double longitudesNE[ points ]  0 or corners to check.
         const double latitudesSW[  points ]  0 or corners to check.
         const double latitudesSE[  points ]  0 or corners to check.
         const double latitudesNW[  points ]  0 or corners to check.
         const double latitudesNE[  points ]  0 or corners to check.
OUTPUTS: double values[ points ]              Values filtered by subset.
RETURNS: size_t number of points in subset.
******************************************************************************/

size_t pointsInSubset( const Bounds domain,
                       const size_t points,
                       const double longitudes[],
                       const double latitudes[],
                       double values[],
                       const double longitudesSW[],
                       const double longitudesSE[],
                       const double longitudesNW[],
                       const double longitudesNE[],
                       const double latitudesSW[],
                       const double latitudesSE[],
                       const double latitudesNW[],
                       const double latitudesNE[] ) {

  size_t result = 0;

  assert( isValidBounds( domain ) );
  assert( points ); assert( longitudes ); assert( latitudes ); assert( values );
  assert( IMPLIES_ELSE( longitudesSW,
                        NON_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                   latitudesSW, latitudesSE,
                                   latitudesNW, latitudesNE ),
                        IS_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                  latitudesSW, latitudesSE,
                                  latitudesNW, latitudesNE ) ) );


  {
    const double longitudeMinimum = domain[ LONGITUDE ][ MINIMUM ];
    const double longitudeMaximum = domain[ LONGITUDE ][ MAXIMUM ];
    const double latitudeMinimum  = domain[ LATITUDE  ][ MINIMUM ];
    const double latitudeMaximum  = domain[ LATITUDE  ][ MAXIMUM ];
    size_t point = 0;

    for ( point = 0; point < points; ++point ) {
      const double value     = values[     point ];
      const double longitude = longitudes[ point ];
      const double latitude  = latitudes[  point ];
      int valid =
        value > MISSING_VALUE &&
        IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ) &&
        IN_RANGE( latitude,  latitudeMinimum,  latitudeMaximum );

      if ( valid && longitudesSW ) { /* Check for degenerate cells: */
        const double longitudeSW = longitudesSW[ point ];
        const double longitudeSE = longitudesSE[ point ];
        const double longitudeNW = longitudesNW[ point ];
        const double longitudeNE = longitudesNE[ point ];
        const double latitudeSW  = latitudesSW[  point ];
        const double latitudeSE  = latitudesSE[  point ];
        const double latitudeNW  = latitudesNW[  point ];
        const double latitudeNE  = latitudesNE[  point ];

        valid =
          longitudeSW != longitude   &&
          longitudeSW != longitudeSE &&
          longitudeSW != longitudeNW &&
          longitudeSW != longitudeNE &&
          longitudeSE != longitude   &&
          longitudeSE != longitudeNW &&
          longitudeSE != longitudeNE &&
          longitudeNW != longitude   &&
          longitudeNW != longitudeNE &&
          longitudeNE != longitude   &&
          latitudeSW  != latitude    &&
          latitudeSW  != latitudeSE  &&
          latitudeSW  != latitudeNW  &&
          latitudeSW  != latitudeNE  &&
          latitudeSE  != latitude    &&
          latitudeSE  != latitudeNW  &&
          latitudeSE  != latitudeNE  &&
          latitudeNW  != latitude    &&
          latitudeNW  != latitudeNE  &&
          latitudeNE  != latitude;
      }

      if ( valid ) {
        ++result;
      } else {
        values[ point ] = MISSING_VALUE;
      }
    }
  }

  assert( result <= points );
  return result;
}



/******************************************************************************
PURPOSE: compactSubsetData - Copy valid data points to the first subsetPoints
         elements of the arrays.
INPUTS:  const size_t subsetPoints      Number of valid points.
         const size_t points            Number of total points.
         double longitudes[   points ]  Longitudes to compact.
         double latitudes[    points ]  Latitudes to compact.
         double values[       points ]  Values to compact.
         double longitudesSW[ points ]  0 or corners to compact.
         double longitudesSE[ points ]  0 or corners to compact.
         double longitudesNW[ points ]  0 or corners to compact.
         double longitudesNE[ points ]  0 or corners to compact.
         double latitudesSW[  points ]  0 or corners to compact.
         double latitudesSE[  points ]  0 or corners to compact.
         double latitudesNW[  points ]  0 or corners to compact.
         double latitudesNE[  points ]  0 or corners to compact.
OUTPUTS: double longitudes[   subsetPoints ]  Longitudes after compacting.
         double latitudes[    subsetPoints ]  Latitudes after compacting.
         double values[       subsetPoints ]  Values after compacting.
         double longitudesSW[ subsetPoints ]  0 or corners after compacting.
         double longitudesSE[ subsetPoints ]  0 or corners after compacting.
         double longitudesNW[ subsetPoints ]  0 or corners after compacting.
         double longitudesNE[ subsetPoints ]  0 or corners after compacting.
         double latitudesSW[  subsetPoints ]  0 or corners after compacting.
         double latitudesSE[  subsetPoints ]  0 or corners after compacting.
         double latitudesNW[  subsetPoints ]  0 or corners after compacting.
         double latitudesNE[  subsetPoints ]  0 or corners after compacting.
******************************************************************************/

void compactSubsetData( const size_t subsetPoints,
                        const size_t points,
                        double longitudes[],
                        double latitudes[],
                        double values[],
                        double longitudesSW[],
                        double longitudesSE[],
                        double longitudesNW[],
                        double longitudesNE[],
                        double latitudesSW[],
                        double latitudesSE[],
                        double latitudesNW[],
                        double latitudesNE[] ) {

  size_t input = 0;
  size_t output = 0;

  assert( subsetPoints ); assert( points > subsetPoints );
  assert( longitudes ); assert( latitudes ); assert( values );
  assert( IMPLIES_ELSE( longitudesSW,
                        NON_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                   latitudesSW, latitudesSE,
                                   latitudesNW, latitudesNE ),
                        IS_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                  latitudesSW, latitudesSE,
                                  latitudesNW, latitudesNE ) ) );


  for ( input = 0; input < points; ++input ) {
    const double value = values[ input ];
    const int valid = value > MISSING_VALUE;
    assert( output <= input );

    if ( valid ) {

      if ( output < input ) { /* Copy valid values to front of arrays: */
        longitudes[ output ] = longitudes[ input ];
        latitudes[  output ] = latitudes[  input ];
        values[     output ] = value;

        if ( longitudesSW ) {
          longitudesSW[ output ] = longitudesSW[ input ];
          longitudesSE[ output ] = longitudesSE[ input ];
          longitudesNW[ output ] = longitudesNW[ input ];
          longitudesNE[ output ] = longitudesNE[ input ];
          latitudesSW[  output ] = latitudesSW[  input ];
          latitudesSE[  output ] = latitudesSE[  input ];
          latitudesNW[  output ] = latitudesNW[  input ];
          latitudesNE[  output ] = latitudesNE[  input ];
        }
      }

      ++output;
    }
  }

  assert( output == subsetPoints );
  assert( isValidLongitude( longitudes[ 0 ] ) );
  assert( isValidLongitude( longitudes[ subsetPoints - 1 ] ) );
  assert( isValidLatitude( latitudes[ 0 ] ) );
  assert( isValidLatitude( latitudes[ subsetPoints - 1 ] ) );
  assert( values[ 0 ] > MISSING_VALUE );
  assert( values[ subsetPoints - 1 ] > MISSING_VALUE );
  assert( IMPLIES( longitudesSW,
                   AND2( isValidLongitude( longitudesSW[ 0 ] ),
                         isValidLongitude( longitudesSW[ subsetPoints - 1] ))));
  assert( IMPLIES( longitudesSE,
                   AND2( isValidLongitude( longitudesSE[ 0 ] ),
                         isValidLongitude( longitudesSE[ subsetPoints - 1] ))));
  assert( IMPLIES( longitudesNW,
                   AND2( isValidLongitude( longitudesNW[ 0 ] ),
                         isValidLongitude( longitudesNW[ subsetPoints - 1] ))));
  assert( IMPLIES( longitudesNE,
                   AND2( isValidLongitude( longitudesNE[ 0 ] ),
                         isValidLongitude( longitudesNE[ subsetPoints - 1] ))));
  assert( IMPLIES( latitudesSW,
                   AND2( isValidLatitude( latitudesSW[ 0 ] ),
                         isValidLatitude( latitudesSW[ subsetPoints - 1 ] ))));
  assert( IMPLIES( latitudesSE,
                   AND2( isValidLatitude( latitudesSE[ 0 ] ),
                         isValidLatitude( latitudesSE[ subsetPoints - 1 ] ))));
  assert( IMPLIES( latitudesNW,
                   AND2( isValidLatitude( latitudesNW[ 0 ] ),
                         isValidLatitude( latitudesNW[ subsetPoints - 1 ] ))));
  assert( IMPLIES( latitudesNE,
                   AND2( isValidLatitude( latitudesNE[ 0 ] ),
                         isValidLatitude( latitudesNE[ subsetPoints - 1 ] ))));
}



/******************************************************************************
PURPOSE: computeCorners - Compute corner vertices given quadrillateral centers.
         const size_t rows                            Rows of data.
         const size_t columns                         Columns of data.
         const double longitudes[   rows * columns ]  Longitudes of centers.
         const double latitudes[    rows * columns ]  Latitudes of centers.
OUTPUTS: double longitudesSW[ rows * columns ]  Cell corner vertices.
         double longitudesSE[ rows * columns ]  Cell corner vertices.
         double longitudesNW[ rows * columns ]  Cell corner vertices.
         double longitudesNE[ rows * columns ]  Cell corner vertices.
         double latitudesSW[  rows * columns ]  Cell corner vertices.
         double latitudesSE[  rows * columns ]  Cell corner vertices.
         double latitudesNW[  rows * columns ]  Cell corner vertices.
         double latitudesNE[  rows * columns ]  Cell corner vertices.
******************************************************************************/

void computeCorners( const size_t rows,
                     const size_t columns,
                     const double longitudes[],
                     const double latitudes[],
                     double longitudesSW[],
                     double longitudesSE[],
                     double longitudesNW[],
                     double longitudesNE[],
                     double latitudesSW[],
                     double latitudesSE[],
                     double latitudesNW[],
                     double latitudesNE[] ) {

  const size_t rows_1    = rows - 1;
  const size_t columns_1 = columns - 1;
  const size_t cells     = rows * columns;
  size_t cell = 0;
  size_t index = 0;

  assert( rows != 0 ); assert( columns != 0 );
  assert( longitudes ); assert( latitudes );
  assert( longitudesSW ); assert( longitudesSE );
  assert( longitudesNW ); assert( longitudesNE );
  assert( latitudesSW ); assert( latitudesSE );
  assert( latitudesNW ); assert( latitudesNE );

#ifndef NDEBUG

  /* Init corners to MISSING_VALUE to ensure they are written exactly once: */

  for ( cell = 0; cell < cells; ++cell ) {
    longitudesSW[ cell ] = MISSING_VALUE;
    longitudesSE[ cell ] = MISSING_VALUE;
    longitudesNW[ cell ] = MISSING_VALUE;
    longitudesNE[ cell ] = MISSING_VALUE;
    latitudesSW[  cell ] = MISSING_VALUE;
    latitudesSE[  cell ] = MISSING_VALUE;
    latitudesNW[  cell ] = MISSING_VALUE;
    latitudesNE[  cell ] = MISSING_VALUE;
  }

#endif

  if ( OR2( rows < 2, columns < 2 ) ) {

    /* Copy all center values to the corners in such degenerate cases: */

#pragma omp parallel for

    for ( cell = 0; cell < cells; ++cell ) {
      longitudesSW[ cell ] =
      longitudesSE[ cell ] =
      longitudesNW[ cell ] =
      longitudesNE[ cell ] = longitudes[ cell ];
      latitudesSW[ cell ] =
      latitudesSE[ cell ] =
      latitudesNW[ cell ] =
      latitudesNE[ cell ] = latitudes[ cell ];
    }

  } else { /* Linearly interpolate and extrapolate the corner points: */
    size_t row    = 0;
    size_t column = 0;

    /*
     * First compute linearly interpolated corners of all interior cells:
     * Note: rows increase north to south and columns increase west to east.
     */

#pragma omp parallel for private( column )

    for ( row = 0; row < rows_1; ++row ) {
      const size_t rowOffset = row * columns;
      const size_t nextRowOffset = rowOffset + columns;

      /* Interior row, interior columns: */

      for ( column = 0; column < columns_1; ++column ) {
        const size_t thisIndex         = rowOffset + column;
        const size_t nextColumn        = thisIndex + 1;
        const size_t nextRow           = nextRowOffset + column;
        const size_t nextRowNextColumn = nextRow + 1;

        const double longitude            = longitudes[ thisIndex ];
        double nextColumnLongitude        = longitudes[ nextColumn ];
        double nextRowLongitude           = longitudes[ nextRow ];
        double nextRowNextColumnLongitude = longitudes[ nextRowNextColumn ];

        const double latitude                  = latitudes[ thisIndex ];
        const double nextColumnLatitude        = latitudes[ nextColumn ];
        const double nextRowLatitude           = latitudes[ nextRow ];
        const double nextRowNextColumnLatitude = latitudes[ nextRowNextColumn];

        if ( OR2( longitude < -179.0, longitude > 179.0 ) ) {
          clampLongitudes( longitude,
                           &nextColumnLongitude,
                           &nextRowLongitude,
                           &nextRowNextColumnLongitude,
                           &nextRowNextColumnLongitude );
        }

        {
          const double interpolatedLongitude = 0.25 *
            ( longitude + nextColumnLongitude + nextRowLongitude +
              nextRowNextColumnLongitude );

          const double interpolatedLatitude = 0.25 *
            ( latitude + nextColumnLatitude + nextRowLatitude +
              nextRowNextColumnLatitude );

          assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                           SIGN( interpolatedLongitude ) == SIGN( longitude)));

          assert( longitudesNE[ thisIndex         ] == MISSING_VALUE );
          assert( longitudesNW[ nextColumn        ] == MISSING_VALUE );
          assert( longitudesSE[ nextRow           ] == MISSING_VALUE );
          assert( longitudesSW[ nextRowNextColumn ] == MISSING_VALUE );
          assert( latitudesNE[  thisIndex         ] == MISSING_VALUE );
          assert( latitudesNW[  nextColumn        ] == MISSING_VALUE );
          assert( latitudesSE[  nextRow           ] == MISSING_VALUE );
          assert( latitudesSW[  nextRowNextColumn ] == MISSING_VALUE );

          longitudesNE[ thisIndex         ] = interpolatedLongitude;
          longitudesNW[ nextColumn        ] = interpolatedLongitude;
          longitudesSE[ nextRow           ] = interpolatedLongitude;
          longitudesSW[ nextRowNextColumn ] = interpolatedLongitude;

          latitudesNE[ thisIndex         ] = interpolatedLatitude;
          latitudesNW[ nextColumn        ] = interpolatedLatitude;
          latitudesSE[ nextRow           ] = interpolatedLatitude;
          latitudesSW[ nextRowNextColumn ] = interpolatedLatitude;
        }
      } /* End loop on interior columns. */

    } /* End parallel loop on interior rows. */

    /* Serial region (not worth parallelizing): */

    /* Last row, interior columns (extrapolated top edge): */

    for ( column = 1, index = rows_1 * columns + 1; column < columns;
          ++column, ++index ) {
      const size_t previousColumn = index - 1;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousColumnLongitude = longitudes[ previousColumn ];
      const int signPreviousColumnLongitude = SIGN( previousColumnLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousColumnLongitude == signLongitude );

      assert( longitudesNW[ index          ] == MISSING_VALUE );
      assert( longitudesNE[ previousColumn ] == MISSING_VALUE );
      assert( latitudesNW[  index          ] == MISSING_VALUE );
      assert( latitudesNE[  previousColumn ] == MISSING_VALUE );

      if ( ! ok ) { /* longitude to previousColumnLongitude crosses +/-180: */
        longitudesNW[ index          ] = signLongitude * EDGE;
        longitudesNE[ previousColumn ] = signPreviousColumnLongitude * EDGE;
        latitudesNW[  index          ] = latitude;
        latitudesNE[  previousColumn ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesSW[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousColumnLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

        const double previousColumnLatitude = latitudes[ previousColumn ];
        const double midpointLatitude =
          0.5 * ( latitude + previousColumnLatitude );
        const double interpolatedLatitude = latitudesSW[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesNW[ index          ] = extrapolatedLongitude;
        longitudesNE[ previousColumn ] = extrapolatedLongitude;
        latitudesNW[  index          ] = extrapolatedLatitude;
        latitudesNE[  previousColumn ] = extrapolatedLatitude;
      }
    }

    /* First row, interior columns (extrapolated bottom edge): */

    for ( column = 1, index = 1; column < columns; ++column, ++index ) {
      const size_t previousColumn = index - 1;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousColumnLongitude = longitudes[ previousColumn ];
      const int signPreviousColumnLongitude = SIGN( previousColumnLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousColumnLongitude == signLongitude );

      assert( longitudesSW[ index          ] == MISSING_VALUE );
      assert( longitudesSE[ previousColumn ] == MISSING_VALUE );
      assert( latitudesSW[  index          ] == MISSING_VALUE );
      assert( latitudesSE[  previousColumn ] == MISSING_VALUE );

      if ( ! ok ) { /* longitude to previousColumnLongitude crosses +/-180: */
        longitudesSW[ index          ] = signLongitude * EDGE;
        longitudesSE[ previousColumn ] = signPreviousColumnLongitude * EDGE;
        latitudesSW[  index          ] = latitude;
        latitudesSE[  previousColumn ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesNW[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousColumnLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          IN_RANGE( extrapolatedLongitude0, -180.0, 180.0 ) ?
            extrapolatedLongitude0
          : signLongitude * EDGE;

        const double previousColumnLatitude = latitudes[ previousColumn ];
        const double midpointLatitude =
          0.5 * ( latitude + previousColumnLatitude );
        const double interpolatedLatitude = latitudesNW[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesSW[ index          ] = extrapolatedLongitude;
        longitudesSE[ previousColumn ] = extrapolatedLongitude;
        latitudesSW[  index          ] = extrapolatedLatitude;
        latitudesSE[  previousColumn ] = extrapolatedLatitude;
      }
    }

    /* First column, interior rows (extrapolated left edge, except corners): */

    for ( row = 1, index = columns; row < rows; ++row, index += columns ) {
      const size_t previousRow = index - columns;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousRowLongitude = longitudes[ previousRow ];
      const int signPreviousRowLongitude = SIGN( previousRowLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousRowLongitude == signLongitude );

      assert( longitudesSW[ index       ] == MISSING_VALUE );
      assert( longitudesNW[ previousRow ] == MISSING_VALUE );
      assert( latitudesSW[  index       ] == MISSING_VALUE );
      assert( latitudesNW[  previousRow ] == MISSING_VALUE );

      if ( ! ok ) { /* longitude to previousRowLongitude crosses +/-180: */
        longitudesSW[ index       ] = signLongitude * EDGE;
        longitudesNW[ previousRow ] = signPreviousRowLongitude * EDGE;
        latitudesSW[  index       ] = latitude;
        latitudesNW[  previousRow ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesSE[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousRowLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

        const double previousRowLatitude = latitudes[ previousRow ];
        const double midpointLatitude =
          0.5 * ( latitude + previousRowLatitude );
        const double interpolatedLatitude = latitudesSE[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesSW[ index       ] = extrapolatedLongitude;
        longitudesNW[ previousRow ] = extrapolatedLongitude;
        latitudesSW[  index       ] = extrapolatedLatitude;
        latitudesNW[  previousRow ] = extrapolatedLatitude;
      }
    }

    /* Last column, interior rows (extrapolated right edge, except corners): */

    for ( row = 1, index = columns + columns - 1;
          row < rows; ++row, index += columns ) {
      const size_t previousRow = index - columns;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousRowLongitude = longitudes[ previousRow ];
      const int signPreviousRowLongitude = SIGN( previousRowLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousRowLongitude == signLongitude );

      assert( longitudesSE[ index       ] == MISSING_VALUE );
      assert( longitudesNE[ previousRow ] == MISSING_VALUE );
      assert( latitudesSE[  index       ] == MISSING_VALUE );
      assert( latitudesNE[  previousRow ] == MISSING_VALUE );

#ifdef DEBUGGING
if ( IN_RANGE( latitude, 30.7, 32.2 ) )
fprintf( stderr, "row %5lu index %5lu (%lf %lf): prl = %lf ok = %d\n",
row, index, longitude, latitude, previousRowLongitude, ok );
#endif

      if ( ! ok ) { /* longitude to previousRowLongitude crosses +/-180: */
        longitudesSE[ index       ] = signLongitude * EDGE;
        longitudesNE[ previousRow ] = signPreviousRowLongitude * EDGE;
        latitudesSE[  index       ] = latitude;
        latitudesNE[  previousRow ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesSW[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousRowLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

        const double previousRowLatitude = latitudes[ previousRow ];
        const double midpointLatitude = 0.5 * (latitude + previousRowLatitude);
        const double interpolatedLatitude = latitudesSW[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesSE[ index       ] = extrapolatedLongitude;
        longitudesNE[ previousRow ] = extrapolatedLongitude;
        latitudesSE[  index       ] = extrapolatedLatitude;
        latitudesNE[  previousRow ] = extrapolatedLatitude;

#ifdef DEBUGGING
if ( IN_RANGE( latitude, 30.7, 32.2 ) )
fprintf( stderr, "  previousRow   %5lu (%lf %lf) [(%lf %lf)-(%lf %lf) (%lf %lf)-(%lf %lf)] ok2 = %d\n",
previousRow,
longitude, latitude,
longitudesNW[ index ], latitudesNW[ index ],
longitudesNE[ index ], latitudesNE[ index ],
longitudesSW[ index ], latitudesSW[ index ],
longitudesSE[ index ], latitudesSE[ index ],
ok2 );
#endif

      }
    }

DEBUG( fprintf( stderr, "\n" ); )

    /* First row, first column cell (extrapolated bottom-left corner): */

    {
      const double latitude               = latitudes[ 0 ];
      const double diagonalLatitude       = latitudesNE[ 0 ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ 0 ];
      const double diagonalLongitude      = longitudesNE[ 0 ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesSW[ 0 ] == MISSING_VALUE );
      assert( latitudesSW[  0 ] == MISSING_VALUE );

      longitudesSW[ 0 ] = extrapolatedLongitude;
      latitudesSW[  0 ] = extrapolatedLatitude;
    }

    /* First row, last column cell (extrapolated bottom-right corner): */

    {
      const double latitude               = latitudes[ columns_1 ];
      const double diagonalLatitude       = latitudesNW[ columns_1 ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude   =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ columns_1 ];
      const double diagonalLongitude      = longitudesNW[ columns_1 ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesSE[ columns_1 ] == MISSING_VALUE );
      assert( latitudesSE[  columns_1 ] == MISSING_VALUE );

      longitudesSE[ columns_1 ] = extrapolatedLongitude;
      latitudesSE[  columns_1 ] = extrapolatedLatitude;
    }

    /* Last row, first column cell (extrapolated top-left corner): */

    index = cells - columns;

    {
      const double latitude               = latitudes[ index ];
      const double diagonalLatitude       = latitudesSE[ index ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude   =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ index ];
      const double diagonalLongitude      = longitudesSE[ index ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesNW[ index ] == MISSING_VALUE );
      assert( latitudesNW[  index ] == MISSING_VALUE );

      longitudesNW[ index ] = extrapolatedLongitude;
      latitudesNW[  index ] = extrapolatedLatitude;
    }

    /* Last row, last column cell (extrapolated top-right corner): */

    index = cells - 1;

    {
      const double latitude               = latitudes[ index ];
      const double diagonalLatitude       = latitudesSW[ index ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude   =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ index ];
      const double diagonalLongitude      = longitudesSW[ index ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesNE[ index ] == MISSING_VALUE );
      assert( latitudesNE[  index ] == MISSING_VALUE );

      longitudesNE[ index ] = extrapolatedLongitude;
      latitudesNE[  index ] = extrapolatedLatitude;
    }

#ifndef NDEBUG

    /* Check corners to ensure they were all computed: */

    for ( cell = 0; cell < cells; ++cell ) {
      assert( longitudesSW[ cell ] != MISSING_VALUE );
      assert( longitudesSE[ cell ] != MISSING_VALUE );
      assert( longitudesNW[ cell ] != MISSING_VALUE );
      assert( longitudesNE[ cell ] != MISSING_VALUE );
      assert( latitudesSW[  cell ] != MISSING_VALUE );
      assert( latitudesSE[  cell ] != MISSING_VALUE );
      assert( latitudesNW[  cell ] != MISSING_VALUE );
      assert( latitudesNE[  cell ] != MISSING_VALUE );
    }

#endif

    /* Clamp any out-of-range values: */

#pragma omp parallel for

    for ( cell = 0; cell < cells; ++cell ) {
      const double longitude = longitudes[ cell ];

      longitudesNW[cell] = CLAMPED_TO_RANGE(longitudesNW[cell], -180.0,180.0);
      longitudesSW[cell] = CLAMPED_TO_RANGE(longitudesSW[cell], -180.0,180.0);
      longitudesSE[cell] = CLAMPED_TO_RANGE(longitudesSE[cell], -180.0,180.0);
      longitudesNE[cell] = CLAMPED_TO_RANGE(longitudesNE[cell], -180.0,180.0);

      latitudesNW[ cell ] = CLAMPED_TO_RANGE( latitudesNW[ cell ], -90.0,90.0);
      latitudesSW[ cell ] = CLAMPED_TO_RANGE( latitudesSW[ cell ], -90.0,90.0);
      latitudesSE[ cell ] = CLAMPED_TO_RANGE( latitudesSE[ cell ], -90.0,90.0);
      latitudesNE[ cell ] = CLAMPED_TO_RANGE( latitudesNE[ cell ], -90.0,90.0);

      if ( OR2( longitude < -179.0, longitude > 179.0 ) ) {
        clampLongitudes( longitude,
                         &longitudesNW[cell],
                         &longitudesSW[cell],
                         &longitudesSE[cell],
                         &longitudesNE[cell] );
      }

      /* Check for bogus (stretched) cells and collapse them to the center: */

      {
        const double maximumDistance = 3.0; /* Degrees. */
        const double longitudeNW = longitudesNW[ cell ];
        const double longitudeSW = longitudesSW[ cell ];
        const double longitudeSE = longitudesSE[ cell ];
        const double longitudeNE = longitudesNE[ cell ];
        double distanceNW =
          longitudeNW < longitude ? longitude - longitudeNW
          : longitudeNW - longitude;
        double distanceSW =
          longitudeSW < longitude ? longitude - longitudeSW
          : longitudeSW - longitude;
        double distanceNE =
          longitudeNE < longitude ? longitude - longitudeNE
          : longitudeNE - longitude;
        double distanceSE =
          longitudeSE < longitude ? longitude - longitudeSE
          : longitudeSE - longitude;
        const int bogusCell =
          distanceNW > maximumDistance ||
          distanceNE > maximumDistance ||
          distanceSW > maximumDistance ||
          distanceSE > maximumDistance;

        if ( bogusCell ) { /* Collapse bogus cells to the cell center point: */
          longitudesSW[ cell ] =
          longitudesSE[ cell ] =
          longitudesNW[ cell ] =
          longitudesNE[ cell ] = longitude;
          latitudesSW[ cell ] =
          latitudesSE[ cell ] =
          latitudesNW[ cell ] =
          latitudesNE[ cell ] = latitudes[ cell ];
        }
      }
    }

  } /* End else non-degenerate cases: */
}



/******************************************************************************
PURPOSE: clampLongitudes - Clamp cell longitudes and match sign of first one.
INPUTS:  const double longitude First longitude center point.
         double* longitude1     Longitude to check/clamp.
         double* longitude2     Longitude to check/clamp.
         double* longitude3     Longitude to check/clamp.
         double* longitude4     Longitude to check/clamp.
OUTPUTS: double* longitude1     Clamped longitude with same sign as longitude.
         double* longitude2     Clamped longitude with same sign as longitude.
         double* longitude3     Clamped longitude with same sign as longitude.
         double* longitude4     Clamped longitude with same sign as longitude.
******************************************************************************/

void clampLongitudes( const double longitude,
                      double* longitude1,
                      double* longitude2,
                      double* longitude3,
                      double* longitude4 ) {

  assert( longitude1 );
  assert( longitude2 );
  assert( longitude3 );
  assert( longitude4 );

  if ( longitude < -179.0 ) {

    if ( *longitude1 >= 0.0 ) {
      *longitude1 = -EDGE;
    }

    if ( *longitude2 >= 0.0 ) {
      *longitude2 = -EDGE;
    }

    if ( *longitude3 >= 0.0 ) {
      *longitude3 = -EDGE;
    }

    if ( *longitude4 >= 0.0 ) {
      *longitude4 = -EDGE;
    }

  } else if ( longitude > 179.0 ) {

    if ( *longitude1 <= 0.0 ) {
      *longitude1 = EDGE;
    }

    if ( *longitude2 <= 0.0 ) {
      *longitude2 = EDGE;
    }

    if ( *longitude3 <= 0.0 ) {
      *longitude3 = EDGE;
    }

    if ( *longitude4 <= 0.0 ) {
      *longitude4 = EDGE;
    }
  }

  assert( SIGN( *longitude1 ) == SIGN( longitude ) );
  assert( SIGN( *longitude2 ) == SIGN( longitude ) );
  assert( SIGN( *longitude3 ) == SIGN( longitude ) );
  assert( SIGN( *longitude4 ) == SIGN( longitude ) );
}



/******************************************************************************
PURPOSE: isLeapYear - Is the yyyy a leap year (i.e., has 366 days)?
INPUTS:  const int yyyy year to check.
RETURNS: int 1 if leap year else 0.
******************************************************************************/

int isLeapYear( const int yyyy ) {
  const int result =
    ( yyyy % 4 == 0 && ( yyyy % 100 != 0 || yyyy % 400 == 0 ) );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidYYYYMMDDHH - Is the timestamp valid?
INPUTS:  const int yyyymmddhh.
RETURNS: int 1 if valid else 0.
******************************************************************************/

int isValidYYYYMMDDHH( const int yyyymmddhh ) {
  const int yyyy = yyyymmddhh / 1000000;
  const int mm   = yyyymmddhh / 10000 % 100;
  const int dd   = yyyymmddhh / 100 % 100;
  const int hh   = yyyymmddhh % 100;
  const int leap = isLeapYear( yyyy );
  const int result =
    IN_RANGE( yyyy, 1900, 9999 ) &&
    IN_RANGE( mm, 1, 12 ) &&
    IN_RANGE( dd, 1, daysPerMonth[ leap ][ mm - 1 ] ) &&
    IN_RANGE( hh, 0, 23 );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidYYYYMMDDHHMM - Is the timestamp valid?
INPUTS:  const long long yyyymmddhhmm.
RETURNS: int 1 if valid else 0.
******************************************************************************/

int isValidYYYYMMDDHHMM( const long long yyyymmddhhmm ) {
  const int yyyymmddhh = (int) ( yyyymmddhhmm / 100 );
  const int mm = (int) ( yyyymmddhhmm % 100 );
  const int result =
    isValidYYYYMMDDHH( yyyymmddhh ) && IN_RANGE( mm, 0, 59 );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidYYYYDDDHHMM - Is the timestamp valid?
INPUTS:  const long long yyyymmddhhmm.
RETURNS: int 1 if valid else 0.
******************************************************************************/

int isValidYYYYDDDHHMM( const long long yyyydddhhmm ) {
  const int yyyy = yyyydddhhmm / 10000000;
  const int ddd  = yyyydddhhmm / 10000 % 1000;
  const int hh   = yyyydddhhmm / 100 % 100;
  const int mm   = yyyydddhhmm % 100;
  const int leap = isLeapYear( yyyy );
  const int result =
    IN_RANGE( yyyy, 1900, 9999 ) &&
    IN_RANGE( ddd, 1, 365 + leap ) &&
    IN_RANGE( hh, 0, 23 ) &&
    IN_RANGE( mm, 0, 59 );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: convertTimestamp - Convert yyyymmddhhmm to yyyydddhhmm
INPUTS:  const long long yyyymmddhhmm  Timestamp to convert.
RETURNS: long long       yyyydddhhmm   Converted timestamp.
******************************************************************************/

long long convertTimestamp( const long long yyyymmddhhmm ) {
  long long result = 0;
  assert( isValidYYYYMMDDHHMM( yyyymmddhhmm ) );

  {
    const int yyyy = yyyymmddhhmm / 100000000;
    const int mo   = yyyymmddhhmm / 1000000 % 100;
    const int dd   = yyyymmddhhmm / 10000 % 100;
    const int hh   = yyyymmddhhmm / 100 % 100;
    const int mm   = yyyymmddhhmm % 100;
    const int leap = isLeapYear( yyyy );
    int ddd = dd;
    int month = 0;

    for ( month = 1; month < mo; ++month ) {
      ddd += daysPerMonth[ leap ][ month - 1 ];
    }

    result = yyyy;
    result *= 1000;
    result += ddd;
    result *= 100;
    result += hh;
    result *= 100;
    result += mm;
  }

  assert( isValidYYYYDDDHHMM( result ) );
  return result;
}



/******************************************************************************
PURPOSE: offsetTimestamp - Compute timestamp + hours.
INPUTS:  const long long yyyydddhhmm  Initial timestamp.
         const long long hours        hours to offset.
RETURNS: long long yyyydddhhmm + hours.
******************************************************************************/

long long offsetTimestamp(const long long yyyydddhhmm, const long long hours) {

  long long result = yyyydddhhmm;

  assert( isValidYYYYDDDHHMM( yyyydddhhmm ) ); assert( hours >= 0 );

  if ( hours > 0 ) {
    const int mm = yyyydddhhmm % 100;
    int yyyy = yyyydddhhmm / 10000000;
    int ddd  = yyyydddhhmm / 10000 % 1000;
    int hh   = yyyydddhhmm / 100 % 100;
    long long hour = hours;

    while ( hour-- ) {
      ++hh;

      if ( hh > 23 ) {
        hh = 0;
        ++ddd;

        if ( ddd > 365 ) {
          const int leap = isLeapYear( yyyy );

          if ( ddd > 365 + leap ) {
            ddd = 1;
            ++yyyy;
          }
        }
      }
    }

    result = yyyy;
    result *= 1000;
    result += ddd;
    result *= 100;
    result += hh;
    result *= 100;
    result += mm;
  }

  assert( isValidYYYYDDDHHMM( result ) );
  assert( IMPLIES_ELSE( hours == 0,
                        result == yyyydddhhmm,
                        result >  yyyydddhhmm ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidBounds - CHeck validity of bounds object.
INPUTS:  const Bounds bounds  Object to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

int isValidBounds( const Bounds bounds) {
  const int result =
    bounds != 0 &&
    IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ) &&
    IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
              bounds[ LONGITUDE ][ MINIMUM ], 180.0 ) &&
    IN_RANGE( bounds[ LATITUDE  ][ MINIMUM ], -90.0, 90.0 ) &&
    IN_RANGE( bounds[ LATITUDE  ][ MAXIMUM ],
              bounds[ LATITUDE  ][ MINIMUM ], 90.0 );
  return result;
}



/******************************************************************************
PURPOSE: boundsOverlap - Do the given bounds overlap?
INPUTS:  const Bounds a  1st bounds to compare.
         const Bounds b  2nd bounds to compare.
RETURNS: int 1 if overlap, else 0.
******************************************************************************/

int boundsOverlap( const Bounds a, const Bounds b ) {
  int result = 0;
  assert( isValidBounds( a ) ); assert( isValidBounds( b ) );
  {
    const int outside =
      a[ LATITUDE  ][ MINIMUM ] > b[ LATITUDE  ][ MAXIMUM ] ||
      a[ LATITUDE  ][ MAXIMUM ] < b[ LATITUDE  ][ MINIMUM ] ||
      a[ LONGITUDE ][ MINIMUM ] > b[ LONGITUDE ][ MAXIMUM ] ||
      a[ LONGITUDE ][ MAXIMUM ] < b[ LONGITUDE ][ MINIMUM ];
    result = ! outside;
  }

  return result;
}



/******************************************************************************
PURPOSE: rotate8ByteArrayIfLittleEndian() - Rotate 8-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array         Array of 8-byte values to rotate.
         const size_t count  Number of items in array.
OUTPUTS: void* array         Array of rotated values.
******************************************************************************/

void rotate8ByteArrayIfLittleEndian( void* array, const size_t count) {

#if IS_LITTLE_ENDIAN

  long long* const array8 = array;
  long long index = 0;
  assert( array ); assert( count > 0 );
  assert( sizeof (long long) == 8 );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const long long value = array8[ index ];
    const long long newValue =
    ( value & 0xff00000000000000LL ) >> 56 |
    ( value & 0x00ff000000000000LL ) >> 40 |
    ( value & 0x0000ff0000000000LL ) >> 24 |
    ( value & 0x000000ff00000000LL ) >>  8 |
    ( value & 0x00000000ff000000LL ) <<  8 |
    ( value & 0x0000000000ff0000LL ) << 24 |
    ( value & 0x000000000000ff00LL ) << 40 |
    ( value & 0x00000000000000ffLL ) << 56;
    array8[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: fileSize - Determine size of named file.
INPUTS:  const char* name  Name of file to examine.
RETURNS: size_t size, in bytes, of named file, else 0 if failed.
******************************************************************************/

size_t fileSize( const char* name ) {
  size_t result = 0;
  struct stat buf;
  assert( name );

  if ( stat( name, &buf ) == -1 ) {
    fprintf( stderr, "\nFailed to determine size of file '%s'.\n", name );
  } else {
    result = buf.st_size;

    if ( result < 1 ) {
      fprintf( stderr, "\nNegative size of file '%s'.\n", name );
      result = 0;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readFile - Read named file into memory and return it as an allocated
         string.
INPUTS:  const char* name  Name of file to examine.
OUTPUTS: size_t* length   Length of string.
RETURNS: char* string contents of file
         (with any '\r' characters converted to ' '),
         else 0 and a failure message is printed to stderr.
******************************************************************************/

char* readFile( const char* name, size_t* length ) {
  char* result = 0;
  assert( name ); assert( length );
  *length = fileSize( name ) / sizeof (char);

  if ( *length > 0 ) {
    const size_t bytes = ( *length + 1 ) * sizeof (char);
    result = malloc( bytes );

    if ( ! result ) {
      *length = 0;
      fprintf( stderr,
              "\nCan't allocate %lu bytes to complete the requested action.\n",
               bytes );
    } else {
      FILE* file = fopen( name, "rb" );

      if ( file ) {
        const size_t itemsRead = fread( result, sizeof (char), *length, file );

        if ( itemsRead != *length ) {
          fprintf( stderr, "\nFailed to read entire file '%s'.\n", name );
          free( result );
          result = 0;
          *length = 0;
        } else {
          result[ *length ] = '\0'; /* Terminate string. */
          controlMToSpace( result );
        }

        fclose( file );
        file = 0;
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: controlMToSpace - Convert any '\r' characters to ' '.
INPUTS:  char* string  String to filter.
OUTPUTS: char* string  Filtered string.
******************************************************************************/

void controlMToSpace( char* string ) {
  char* s = 0;
  assert( string );

  for ( s = string; *s; ++s ) {

    if ( *s == '\r' ) {
      *s = ' ';
    }
  }
}



/******************************************************************************
PURPOSE: spacesToUnderscores - Convert any ' ' characters to '_'.
INPUTS:  char* string  String to filter.
OUTPUTS: char* string  Filtered string.
******************************************************************************/

void spacesToUnderscores( char* string ) {
  char* s = 0;
  assert( string );

  for ( s = string; *s; ++s ) {

    if ( *s == ' ' ) {
      *s = '_';
    }
  }
}



/******************************************************************************
PURPOSE: linesInString - Count number of lines in a string.
INPUTS:  const char* string  String to scan.
RETURNS: size_t number of lines in string.
******************************************************************************/

size_t linesInString( const char* string ) {
  size_t result = 0;
  const char* c = 0;
  assert( string );

  for ( c = string; *c; ++c ) {
    result += *c == '\n';
  }

  return result;
}



/******************************************************************************
PURPOSE: Utilities.c - Some general-purpose reusable routines.

HISTORY: 2017-01-02 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/



/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, stderr, fprintf(). */
#include <stdlib.h>    /* For malloc(), free(), abs(). */
#include <float.h>     /* For DBL_MAX. */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

#include "Utilities.h" /* For public interface. */

/*================================= MACROS ==================================*/

/*
 * Is the platform big-endian (MSB: most significant byte first) or
 * little-endian (LSB: least significant byte first)?
 */

#ifndef IS_LITTLE_ENDIAN
#define IS_LITTLE_ENDIAN \
  ( \
    ( defined(__BYTE_ORDER__) && \
      defined(__ORDER_LITTLE_ENDIAN__) && \
      __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ ) || \
    defined(__x86_64__) || \
    defined(__i386__)   || \
    defined(__i486__)   || \
    defined(__i586__)   || \
    defined(__i686__)   || \
    defined(__alpha)    || \
    defined(__ia64__)   || \
    defined(__ARMEL__)  || \
    defined(__MIPSEL)   || \
    defined(_WIN32)     || \
    defined(_WIN64)     || \
    defined(_M_IX86)    || \
    defined(_M_AMD64)      \
  )
#endif

/*============================ GLOBAL CONSTANTS =============================*/

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
PURPOSE: isValidElevation - Is the argument a valid elevation?
INPUTS:  const double elevation  In meters above mean sea level.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidElevation( const double elevation ) {
  const int result = IN_RANGE( elevation, -500.0, 1e5 );
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
PURPOSE: computeBounds - Copmute range of longitudes and latitudes.
INPUTS:  const size_t count                  Number of lon-lat values.
         const double longitudes[ count ]    Longitudes to check.
         const double latitudes[  count ]    Latitudes  to check.
OUTPUTS: Bounds bounds                       Bounds of lon-lat values.
******************************************************************************/

void computeBounds( const size_t count,
                    const double longitudes[], const double latitudes[],
                    Bounds bounds ) {

  size_t index = 0;
  double longitudeMinimum = 0.0;
  double longitudeMaximum = 0.0;
  double latitudeMinimum  = 0.0;
  double latitudeMaximum  = 0.0;
  assert( count ); assert( longitudes ); assert( latitudes ); assert( bounds );

  longitudeMinimum = longitudeMaximum = *longitudes;
  latitudeMinimum  = latitudeMaximum  = *latitudes;

  for ( index = 1; index < count; ++index ) {
    const double longitude = longitudes[ index ];
    const double latitude  = latitudes[  index ];

    if ( longitude < longitudeMinimum ) {
      longitudeMinimum = longitude;
    } else if ( longitude > longitudeMaximum ) {
      longitudeMaximum = longitude;
    }

    if ( latitude < latitudeMinimum ) {
      latitudeMinimum = latitude;
    } else if ( latitude > latitudeMaximum ) {
      latitudeMaximum = latitude;
    }
  }

  bounds[ LONGITUDE ][ MINIMUM ] = longitudeMinimum;
  bounds[ LONGITUDE ][ MAXIMUM ] = longitudeMaximum;
  bounds[ LATITUDE  ][ MINIMUM ] = latitudeMinimum;
  bounds[ LATITUDE  ][ MAXIMUM ] = latitudeMaximum;

  DEBUG( fprintf( stderr, "scan bounds: [%lf %lf][%lf %lf]\n",
                  bounds[ LONGITUDE ][ MINIMUM ],
                  bounds[ LONGITUDE ][ MAXIMUM ],
                  bounds[ LATITUDE  ][ MINIMUM ],
                  bounds[ LATITUDE  ][ MAXIMUM ] ) );
  assert( isValidBounds( (const double (*)[2]) bounds ) );
}



/******************************************************************************
PURPOSE: compactPointsInSubset - Compact data inside domain and output
         reduced dimensions.
INPUTS:  const Bounds domain                  Domain to subset data to.
         const double minimumElevation        Minimum elevation in subset.
         const double maximumElevation        Maximum elevation in subset.
         const size_t points                  Ground points of data.
         const size_t levels                  Vertical levels of data.
         double timestamps[  points ]           Timestamps before compacting.
         double longitudes[  points ]           Longitudes before compacting.
         double latitudes[   points ]           Latitudes  before compacting.
         double elevations[  points * levels ]  Elevations before compacting.
         double values[      points * levels]   Values before compacting.
         double thicknesses[ points * levels]   0 or thicknesses before compact
OUTPUTS: double timestamps[  subsetPoints ]     Timestamps compacted.
         double longitudes[  subsetPoints ]     Longitudes compacted.
         double latitudes[   subsetPoints ]     Latitudes  compacted.
         double elevations[  subsetPoints * subsetLevels ] Elevations compacted
         double values[      subsetPoints * subsetLevels]  Values compacted.
         double thicknesses[ subsetPoints * subsetLevels]  Thicknesses compacted
         size_t* subsetPoints  Number of ground points in subset domain.
         size_t* subsetLevels  Number of level points in subset elevation range
RETURNS: int 1 if there is at least one valid point in the subset, else 0.
******************************************************************************/

int compactPointsInSubset( const Bounds domain,
                           const double minimumElevation,
                           const double maximumElevation,
                           const size_t points,
                           const size_t levels,
                           double timestamps[],
                           double longitudes[],
                           double latitudes[],
                           double elevations[],
                           double values[],
                           double thicknesses[],
                           size_t* subsetPoints,
                           size_t* subsetLevels ) {

  int result = 0;
  size_t level = 0;

  assert( isValidBounds( domain ) );
  assert( isValidElevation( minimumElevation ) );
  assert( isValidElevation( maximumElevation ) );
  assert( minimumElevation <= maximumElevation );
  assert( points ); assert( levels );
  assert( timestamps ); assert( longitudes ); assert( latitudes );
  assert( elevations ); assert( values );
  assert( subsetPoints ); assert( subsetLevels );

  *subsetPoints = *subsetLevels = 0;

  /*
   * Every ground point has the same set of elevation values.
   * Compute elevations levels in subset elevation range:
   */

  while ( AND2( level < levels, elevations[ level ] < minimumElevation ) ) {
    ++level;
  }

  if ( level < levels ) {
    const double longitudeMinimum = domain[ LONGITUDE ][ MINIMUM ];
    const double longitudeMaximum = domain[ LONGITUDE ][ MAXIMUM ];
    const double latitudeMinimum  = domain[ LATITUDE  ][ MINIMUM ];
    const double latitudeMaximum  = domain[ LATITUDE  ][ MAXIMUM ];
    size_t point = 0;
    size_t index = 0;
    size_t output = 0;
    size_t output2 = 0;
    size_t firstElevationIndex = level;
    size_t lastElevationIndex  = level;

    do {
      lastElevationIndex = level;
      *subsetLevels += 1;
      ++level;
    } while ( AND2( level < levels, elevations[ level ] <= maximumElevation ) );

    assert( *subsetLevels );
    assert( lastElevationIndex < levels );
    assert( firstElevationIndex <= lastElevationIndex );

    /* For each gound point in domain, copy/compact data and elevations: */

    for ( point = 0, index = 0; point < points; ++point ) {
      const double longitude = longitudes[ point ];
      const double latitude  = latitudes[  point ];
      const int inDomain =
        AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
              IN_RANGE( latitude,  latitudeMinimum,  latitudeMaximum ) );

      if ( ! inDomain ) {
        index += levels;
      } else {
        *subsetPoints += 1;

        if ( output < point ) {
          timestamps[ output ] = timestamps[ point ];
          longitudes[ output ] = longitudes[ point ];
          latitudes[  output ] = latitudes[  point ];
        }

        ++output;

        index += firstElevationIndex;

        for ( level = firstElevationIndex; level <= lastElevationIndex;
              ++level, ++index ) {
      
          if ( output2 < index ) {
            values[ output2 ] = values[ index ];
            elevations[ output2 ] = elevations[ index ];

            if ( thicknesses ) {
              thicknesses[ output2 ] = thicknesses[ index ];
            }
          }

          ++output2;
        }

        index += levels - index - 1; /* Skip out-of-bounds levels. */
        assert( ( index + 1 ) % levels == 0 );
      }
    }
  }

  result = *subsetPoints != 0;
  assert( IS_BOOL( result ) );
  assert( *subsetPoints <= points );
  assert( *subsetLevels <= levels );
  return result;
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
INPUTS:  const size_t count  Number of items in array.
         void* array         Array of 8-byte values to rotate.
OUTPUTS: void* array         Array of rotated values.
******************************************************************************/

void rotate8ByteArrayIfLittleEndian( const size_t count, void* array ) {

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
PURPOSE: reverseLevels() - Reverse the order of the data levels.
INPUTS:  const size_t points  Number of ground points of data.
         const size_t levels  Number of vertical points per ground point.
OUTPUTS: double array[ points * levels ]  Array of reversed values.
******************************************************************************/

void reverseLevels( const size_t points, const size_t levels, double array[] ) {

  size_t point = 0;
  size_t levels2 = levels / 2;

  assert( points ); assert( levels ); assert( array );

#pragma omp parallel for

  for ( point = 0; point < points; ++point ) {
    const size_t offset = point * levels;
    size_t level = 0;

    for ( level = 0; level < levels2; ++level ) {
      const size_t index1 = offset + level;
      const size_t index2 = offset + levels - ( level + 1 );
      const double swapTemp = array[ index1 ];
      array[ index1 ] = array[ index2 ];
      array[ index2 ] = swapTemp;
    }
  }
}




/******************************************************************************
PURPOSE: expandInt8 - Expand signed 8-bit integer values in first half of
         array to doubles, in place.
INPUTS:  size_t count           Number of items in array.
         double array[ count ]  Array with int8 values stored in lower half.
OUTPUTS: double array[ count ]  Array with double values throughout.
******************************************************************************/

void expandInt8( size_t count, double array[] ) {
  signed char* iarray = (signed char*) array;
  assert( array ); assert( count >= 1 );
  array  += count;
  iarray += count;

  while ( count-- ) {
    *--array = *--iarray;
  }
}



/******************************************************************************
PURPOSE: expandUint8 - Expand unsigned 8-bit integer values in first half of
         array to doubles, in place.
INPUTS:  size_t count           Number of items in array.
         double array[ count ]  Array with int8 values stored in lower half.
OUTPUTS: double array[ count ]  Array with double values throughout.
******************************************************************************/

void expandUint8( size_t count, double array[] ) {
  unsigned char* iarray = (unsigned char*) array;
  assert( array ); assert( count >= 1 );
  array  += count;
  iarray += count;

  while ( count-- ) {
    *--array = *--iarray;
  }
}



/******************************************************************************
PURPOSE: expandInt16 - Expand signed 16-bit integer values in first half of
         array to doubles, in place.
INPUTS:  size_t count           Number of items in array.
         double array[ count ]  Array with int8 values stored in lower half.
OUTPUTS: double array[ count ]  Array with double values throughout.
******************************************************************************/

void expandInt16( size_t count, double array[] ) {
  signed short* iarray = (signed short*) array;
  assert( array ); assert( count >= 1 );
  array  += count;
  iarray += count;

  while ( count-- ) {
    *--array = *--iarray;
  }
}



/******************************************************************************
PURPOSE: expandUint16 - Expand unsigned 16-bit integer values in first half of
         array to doubles, in place.
INPUTS:  size_t count           Number of items in array.
         double array[ count ]  Array with int8 values stored in lower half.
OUTPUTS: double array[ count ]  Array with double values throughout.
******************************************************************************/

void expandUint16( size_t count, double array[] ) {
  unsigned short* iarray = (unsigned short*) array;
  assert( array ); assert( count >= 1 );
  array  += count;
  iarray += count;

  while ( count-- ) {
    *--array = *--iarray;
  }
}



/******************************************************************************
PURPOSE: expandInt32 - Expand signed 32-bit integer values in first half of
         array to doubles, in place.
INPUTS:  size_t count           Number of items in array.
         double array[ count ]  Array with int8 values stored in lower half.
OUTPUTS: double array[ count ]  Array with double values throughout.
******************************************************************************/

void expandInt32( size_t count, double array[] ) {
  signed int* iarray = (signed int*) array;
  assert( array ); assert( count >= 1 );
  array  += count;
  iarray += count;

  while ( count-- ) {
    *--array = *--iarray;
  }
}



/******************************************************************************
PURPOSE: expandUint32 - Expand unsigned 32-bit integer values in first half of
         array to doubles, in place.
INPUTS:  size_t count           Number of items in array.
         double array[ count ]  Array with int8 values stored in lower half.
OUTPUTS: double array[ count ]  Array with double values throughout.
******************************************************************************/

void expandUint32( size_t count, double array[] ) {
  unsigned int* iarray = (unsigned int*) array;
  assert( array ); assert( count >= 1 );
  array  += count;
  iarray += count;

  while ( count-- ) {
    *--array = *--iarray;
  }
}


/******************************************************************************
PURPOSE: expandReals - Expand float values in first half of array to doubles,
         in place.
INPUTS:  size_t count           Number of items in array.
         double array[ count ]  Array with int8 values stored in lower half.
OUTPUTS: double array[ count ]  Array with double values throughout.
******************************************************************************/

void expandReals( size_t count, double array[] ) {
  float* farray = (float*) array;
  assert( array ); assert( count >= 1 );
  array  += count;
  farray += count;

  while ( count-- ) {
    *--array = *--farray;
  }
}



/******************************************************************************
PURPOSE: scaleValues - Multiply each array element by factor.
INPUTS:  const double factor    Factor to multiply array elements by.
         const size_t count     Number of elements in array.
         double array[ count ]  Array of values.
OUTPUTS: double array[ count ]  Array of values scaled by factor.
******************************************************************************/

void scaleValues( const double factor, const size_t count, double array[] ) {
  size_t index = 0;
  assert( factor ); /* Non-zero multiplier. */
  assert( IN_RANGE( factor, -DBL_MAX, DBL_MAX ) ); /* Finite, non-NAN. */
  assert( count ); assert( array );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    array[ index ] *= factor;
  }
}



/******************************************************************************
 PURPOSE: offsetValues - Add an offset to each array element.
 INPUTS:  const size_t  count          Number of elements in array.
          double array[ count ]  Array of values.
 OUTPUTS: double array[ count ]  Array of values scaled by factor.
 ******************************************************************************/

void offsetValues( const double offset, const size_t count, double array[] ) {
  size_t index = 0;
  assert( offset ); /* Non-zero multiplier. */
  assert( IN_RANGE( offset, -DBL_MAX, DBL_MAX ) ); /* Finite, non-NAN. */
  assert( count ); assert( array );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    array[ index ] += offset;
  }
}



/******************************************************************************
PURPOSE: copyVectorComponent - Copy a component of an array of vectors.
INPUTS:  const size_t count          Number of vectors in input array.
         const size_t components     Number of components per vector.
         const size_t component      Components to copy.
         const double input[ count * components ]  Array of vectors to read.
OUTPUTS: double output[ count ]                    Array of scalars copied.
******************************************************************************/

void copyVectorComponent( const size_t count, const size_t components,
                          const size_t component,
                          const double input[], double output[] ) {

  size_t index = 0;
  size_t index2 = component;
  assert( count ); assert( components );
  assert( IN_RANGE( component, 0, components - 1 ) );
  assert( input ); assert( output );

/* #pragma omp parallel for private( index2 ) */

  for ( index = 0, index2 = component; index < count;
        ++index, index2 += components ) {
    output[ index ] = input[ index2 ];
  }
}



/******************************************************************************
PURPOSE: copyMaximumComponent - Copy the maximum of component of an array of
         vectors.
INPUTS:  const size_t count                   Number of vectors in input array.
         const size_t components              Number of components per vector.
         const double input[ count * components ]  Array of vectors to read.
OUTPUTS: double output[ count ]           Array of scalars copied.
******************************************************************************/

void copyMaximumComponent( const size_t count, const size_t components,
                           const double input[], double output[] ) {

  size_t index = 0;
  assert( count ); assert( components );
  assert( input ); assert( output );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const size_t offset = index * components;
    size_t component = 0;
    double maximum = MISSING_VALUE;

    for ( component = 0; component < components; ++component ) {
      const double value = input[ offset + component ];

      if ( value > maximum ) {
        maximum = value;
      }
    }

    output[ index ] = maximum;
  }
}



/******************************************************************************
PURPOSE: copyMeanComponents - Copy mean of components of an array of vectors.
INPUTS:  const size_t count               Number of vectors in input array.
         const size_t components          Number of components per vector.
const double input[ count * components ]  Array of vectors to read.
OUTPUTS: double output[ count ]           Array of scalars copied.
******************************************************************************/

void copyMeanComponents( const size_t count, const size_t components,
                         const double input[], double output[] ) {

  size_t index = 0;
  assert( count ); assert( components );
  assert( input ); assert( output );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const size_t offset = index * components;
    size_t validCount = 0;
    size_t component = 0;
    double sum = 0.0;

    for ( component = 0; component < components; ++component ) {
      const double value = input[ offset + component ];

      if ( value > MISSING_VALUE ) {
        sum += value;
        ++validCount;
      }
    }

    if ( validCount ) {
      output[ index ] = sum / validCount;
    } else {
      output[ index ] = MISSING_VALUE;
    }
  }
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



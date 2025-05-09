/******************************************************************************
PURPOSE: Utilities.c - Some general-purpose reusable routines.

HISTORY: 2017-10-14 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/



/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, stderr, fprintf(). */
#include <stdlib.h>    /* For malloc(), free(). */
#include <string.h>    /* For strlen(). */
#include <float.h>     /* For FLT_MAX. */
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

#define AND2(a, b) ( (a) && (b) )
#define AND4(a, b, c, d) ( (a) && (b) && (c) && (d) )

/*============================ GLOBAL CONSTANTS =============================*/

/*
 * 30 days hath September, April, June and November, all the rest have 31,
 * except February which has either 28 or 29 (on a leap year).
 */

static const int daysPerMonth[ 2 ][ 12 ] =
{
  { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
  { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
};

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: pointsInDomain - Compute number of points in subset domain and output
         mask.
INPUTS:  const Bounds domain                Domain to subset data to.
         const size_t points                Points of data.
         const double longitudes[ points ]  Longitudes to check.
         const double latitudes[  points ]  Latitudes to check.
OUTPUTS: unsigned char mask[points] mask[point] = 1 if point in domain, else 0.
RETURNS: size_t number of points in subset.
******************************************************************************/

size_t pointsInDomain( const Bounds domain,
                       const size_t points,
                       const double longitudes[],
                       const double latitudes[],
                       unsigned char mask[] ) {
  size_t result = 0;

  assert( isValidBounds( domain ) ); assert( points );
  assert( longitudes ); assert( latitudes );
  assert( mask );

  {
    const double longitudeMinimum = domain[ LONGITUDE ][ MINIMUM ];
    const double longitudeMaximum = domain[ LONGITUDE ][ MAXIMUM ];
    const double latitudeMinimum  = domain[ LATITUDE  ][ MINIMUM ];
    const double latitudeMaximum  = domain[ LATITUDE  ][ MAXIMUM ];
    size_t point = 0;

    for ( point = 0; point < points; ++point ) {
      mask[ point ] = 0;
      const double longitude = longitudes[ point ];

      if ( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ) ) {
        const double latitude = latitudes[ point ];

        if ( IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) ) {
          mask[ point ] = 1;
          ++result;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: isLeapYear - Is the yyyy a leap year (i.e., has 366 days)?
INPUTS:  const int yyyy year to check.
RETURNS: int 1 if leap year else 0.
******************************************************************************/

int isLeapYear( const int yyyy ) {
  return yyyy % 4 == 0 && ( yyyy % 100 != 0 || yyyy % 400 == 0 );
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
  return result;
}



/******************************************************************************
PURPOSE: incrementHours - Increment yyyymmddhh by hours.
INPUTS:  const int yyyymmddhh   Timestamp to increment.
         const int hours        Number of hours to increment.
RETURNS: int incremented timestamp.
******************************************************************************/

int incrementHours( const int yyyymmddhh, const int hours ) {
  int result = yyyymmddhh;
  assert( isValidYYYYMMDDHH( yyyymmddhh ) ); assert( hours >= 0 );

  if ( hours > 0 ) {
    int hour = 0;
    int yyyy = yyyymmddhh / 1000000;
    int mm   = yyyymmddhh / 10000 % 100;
    int dd   = yyyymmddhh / 100 % 100;
    int hh   = yyyymmddhh % 100;
    int leap = isLeapYear( yyyy );

    for ( hour = 0; hour < hours; ++hour ) {
      ++hh;

      if ( hh > 23 ) {
        hh = 0;
        ++dd;

        if ( dd > daysPerMonth[ leap ][ mm - 1 ] ) {
          dd = 1;
          ++mm;

          if ( mm > 12 ) {
            mm = 1;
            ++yyyy;
            leap = isLeapYear( yyyy );
          }
        }
      }
    }

    result = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
  }

  assert( isValidYYYYMMDDHH( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidBounds - CHeck validity of bounds object.
INPUTS:  const Bounds bounds  Object to check.
RETURNS: size_t 1 if valid, else 0.
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
PURPOSE: subsetIndicesByBounds - Subset row and column indices by bounds.
INPUTS:  const Bounds bounds    Longitude-latitude bounds of subset.
         const size_t rows     Number of row points.
         const size_t columns  Number of column points.
         const double longitudes[ rows ][ columns ]  Longitudes.
         const double latitudes[  rows ][ columns ]  Latitudes.
OUTPUTS: size_t* firstRow      Index of first row of subset.
         size_t* lastRow       Index of last  row of subset.
         size_t* firstColumn   Index of first column of subset.
         size_t* lastColumn    Index of last  column of subset.
RETURNS: int 1 if a non-empty subset exists, else 0 and outputs are 0.
******************************************************************************/

int subsetIndicesByBounds( const Bounds bounds,
                           const size_t rows,
                           const size_t columns,
                           const double longitudes[],
                           const double latitudes[],
                           size_t* firstRow,
                           size_t* lastRow,
                           size_t* firstColumn,
                           size_t* lastColumn ) {

  int result = 0;

  assert( isValidBounds( bounds ) );
  assert( rows ); assert( columns );
  assert( IN_RANGE( longitudes[ 0 ], -180.0, 180.0 ) ),
  assert( IN_RANGE( longitudes[ rows * columns - 1 ], -180.0, 180.0 ) ),
  assert( IN_RANGE( latitudes[ 0 ], -90.0, 90.0 ) ),
  assert( IN_RANGE( latitudes[ rows * columns - 1 ], -90.0, 90.0 ) ),
  assert( firstRow ); assert( lastRow );
  assert( firstColumn ); assert( lastColumn );

  {
    const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
    const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
    const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
    const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
    size_t theFirstRow    = 0;
    size_t theLastRow     = 0;
    size_t theFirstColumn = 0;
    size_t theLastColumn  = 0;
    size_t row            = 0;
    size_t column         = 0;

    *firstRow = *lastRow = *firstColumn = *lastColumn = 0;

    /* Loop forward through rows to find first subset row: */

    for ( row = 0; row < rows; ++row ) {
      const size_t rowOffset = row * columns;
      assert( rowOffset + columns - 1 < rows * columns );

      for ( column = 0; column < columns; ++column ) {
        const size_t index = rowOffset + column;
        const double longitude = longitudes[ index ];
        const double latitude  = latitudes[ index ];
        const size_t inside =
          AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

        if ( inside ) {
          theFirstRow = theLastRow = row;
          result = 1;
          row = rows;
          column = columns;
        }
      }
    }

    if ( result ) {
      assert( IN_RANGE( theFirstRow, 0, rows - 1 ) );

      /* Loop backward through rows to find last subset row: */

      for ( row = rows; row > theFirstRow; ) {
        const size_t rowOffset = --row * columns;
        assert( rowOffset + columns - 1 < rows * columns );

        for ( column = 0; column < columns; ++column ) {
          const size_t index = rowOffset + column;
          const double longitude = longitudes[ index ];
          const double latitude  = latitudes[ index ];
          const size_t inside =
            AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                  IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

          if ( inside ) {
            theLastRow = row;
            row = 0;
            column = columns;
          }
        }
      }

      assert( IN_RANGE( theLastRow, theFirstRow, rows - 1 ) );

      /* Loop forward through columns to find first subset column: */

      for ( column = 0; column < columns; ++column ) {

        for ( row = theFirstRow; row <= theLastRow; ++row ) {
          const size_t index = row * columns + column;
          const double longitude = longitudes[ index ];
          const double latitude  = latitudes[ index ];
          const size_t inside =
            AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                  IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

          if ( inside ) {
            theLastColumn = theFirstColumn = column;
            row = rows;
            column = columns;
          }
        }
      }

      assert( IN_RANGE( theFirstColumn, 0, columns - 1 ) );

      /* Loop backward through columns to find last subset column: */

      for ( column = columns; column > theFirstColumn; ) {
        --column;

        for ( row = theFirstRow; row <= theLastRow; ++row ) {
          const size_t index = row * columns + column;
          const double longitude = longitudes[ index ];
          const double latitude  = latitudes[ index ];
          const size_t inside =
            AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                  IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

          if ( inside ) {
            theLastColumn = column;
            row = rows;
            column = 0;
          }
        }
      }

      *firstRow    = theFirstRow;
      *lastRow     = theLastRow;
      *firstColumn = theFirstColumn;
      *lastColumn  = theLastColumn;
    }
  }

  assert( result == 0 ||
          AND4( IN_RANGE( *firstRow, 0, rows - 1 ),
                IN_RANGE( *lastRow, *firstRow, rows - 1 ),
                IN_RANGE( *firstColumn, 0, columns - 1 ),
                IN_RANGE( *lastColumn, *firstColumn, columns - 1 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: doublesToFloats - Copy 64-bit floating-point values to 32-bit values.
INPUTS:  double* array       Array of 64-bit values to expand in-place.
         const size_t count  Number of values in array.
OUTPUTS: float* array        Compressed array of 32-bit values.
NOTES:   Values in the array are overwritten to the first half and will have
         lost precision and will be clamped to 32-bit range. UGLY.
******************************************************************************/

void doublesToFloats( double* array, const size_t count ) {
  float* destination = (float*) array;
  size_t counter = count;

  while ( counter-- ) {
    float value = *array;

    if ( ! ( value > -FLT_MAX ) ) { /* Clamps NaNs and Infs. */
      value = -FLT_MAX;
    } else if ( value > FLT_MAX ) {
      value = FLT_MAX;
    }

    *destination++ = value;
    ++array;
  }
}



/******************************************************************************
PURPOSE: rotate4ByteArrayIfLittleEndian() - Rotate 4-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array         Array of 8-byte values to rotate.
         const size_t count  Number of items in array.
OUTPUTS: void* array         Array of rotated values.
******************************************************************************/

void rotate4ByteArrayIfLittleEndian( void* array, const size_t count ) {

#if IS_LITTLE_ENDIAN

  unsigned int* const array4 = array;
  long long index = 0;
  assert( array ); assert( count > 0 );
  assert( sizeof (int) == 4 );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const unsigned int value = array4[ index ];
    const unsigned int newValue =
    ( value & 0xff000000 ) >> 24 |
    ( value & 0x00ff0000 ) >>  8 |
    ( value & 0x0000ff00 ) <<  8 |
    ( value & 0x000000ff ) << 24;
    array4[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/******************************************************************************
PURPOSE: rotate8ByteArrayIfLittleEndian() - Rotate 8-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array         Array of 8-byte values to rotate.
         const size_t count  Number of items in array.
OUTPUTS: void* array         Array of rotated values.
******************************************************************************/

void rotate8ByteArrayIfLittleEndian( void* array, const size_t count ) {

#if IS_LITTLE_ENDIAN

  unsigned long long* const array8 = array;
  long long index = 0;
  assert( array ); assert( count > 0 );
  assert( sizeof (long long) == 8 );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const unsigned long long value = array8[ index ];
    const unsigned long long newValue =
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
PURPOSE: fillArray() - Write value to each element of an array.
INPUTS:  const double value     Value to write.
         const size_t count     Number of items in array.
OUTPUTS: double array[ count ]  Array initialized to value.
******************************************************************************/

void fillArray( const double value, const size_t count, double array[] ) {
  long long index = 0;
  assert( array ); assert( count > 0 );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    array[ index ] = value;
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




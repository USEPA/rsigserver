
/******************************************************************************
PURPOSE: elevation.c - Read a sequence of longitude-latitude pairs from stdin
         and write their corresponding elevation (in meters above mean sea
         level) to stdout.
NOTES:   Latitude is assumed to be on a WGS84 spheroid.
         The following files are read from /data/land_use/
         grid_surface_-180_-90_-90_0.bin
         grid_surface_-180_-90_0_90.bin
         grid_surface_-90_0_-90_0.bin
         grid_surface_-90_0_0_90.bin
         grid_surface_0_90_-90_0.bin
         grid_surface_0_90_0_90.bin
         grid_surface_90_180_-90_0.bin
         grid_surface_90_180_0_90.bin
HISTORY: 2011-10-19 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For FILE, stdin, stdout, fscanf(), fprintf(), fread().*/
#include <stdlib.h> /* For malloc(), free(), atexit(). */
#include <string.h> /* For memset(). */
#include <limits.h> /* For ULONG_MAX. */
#include <sys/stat.h> /* For struct stat, stat(). */

/*================================== MACROS =================================*/

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

/* Compile-time assertion: */

#define assert_static(a) extern int unused_assert_static_[ (a) ? 1 : -1 ]

/* Logic macros: */

#define IS_BOOL(c) ((c)==0||(c)==1)
#define IMPLIES_ELSE(p,c1,c2) (((p)&&(c1))||((!(p))&&(c2)))
#define OR2(a,b) ((a)||(b))
#define AND2(a,b) ((a)&&(b))
#define AND6(a,b,c,d,e,f) ((a)&&(b)&&(c)&&(d)&&(e)&&(f))
#define AND7(a,b,c,d,e,f,g) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g))
#define AND8(a,b,c,d,e,f,g,h) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h))
#define IS_ZERO6(a,b,c,d,e,f) \
  ((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0)
#define IS_ZERO8(a,b,c,d,e,f,g,h) \
  ((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0)
#define IN_RANGE( value, minimum , maximum ) \
  ( (value) >= (minimum) && (value) <= (maximum) )
#define CLAMPED_TO_RANGE( value, low, high ) \
  ((value) < (low) ? (low) : (value) > (high) ? (high) : (value))

/* Optionally include debug statements if -DDEBUGGING: */

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(s)
#endif

/* Memory allocation/deallocation: */

#define FREE( p ) ( ( (p) ? free(p) : (void) 0 ), (p) = 0 )

/*================================== TYPES ==================================*/

typedef struct {
  const char* name;
  FILE* file;
  int rows;
  int columns;
  float longitudeMinimum;
  float longitudeMaximum;
  float latitudeMinimum;
  float latitudeMaximum;
  float* data;
} ElevationFile;

/*============================= GLOBAL VARIABLES ============================*/

static int registeredCloseAtExit = 0;

#ifdef USE_CWD
#define DIRECTORY "."
#else
#define DIRECTORY "/data/land_use"
#endif

static ElevationFile elevationFiles[] = {
  {
    DIRECTORY"/grid_surface_-180_-90_-90_0.bin",
    0, 0, 0, -180.0, -90.0, -90.0, 0.0, 0
  },
  {
    DIRECTORY"/grid_surface_-180_-90_0_90.bin",
    0, 0, 0, -180.0, -90.0, 0.0, 90.0, 0
  },
  {
    DIRECTORY"/grid_surface_-90_0_-90_0.bin",
    0, 0, 0, -90.0, 0.0, -90.0, 0.0, 0
  },
  {
    DIRECTORY"/grid_surface_-90_0_0_90.bin",
    0, 0, 0, -90.0, 0.0, 0.0, 90.0, 0
  },
  {
    DIRECTORY"/grid_surface_0_90_-90_0.bin",
    0, 0, 0, 0.0, 90.0, -90.0, 0.0, 0
  },
  {
    DIRECTORY"/grid_surface_0_90_0_90.bin",
    0, 0, 0, 0.0, 90.0, 0.0, 90.0, 0
  },
  {
    DIRECTORY"/grid_surface_90_180_-90_0.bin",
    0, 0, 0, 90.0, 180.0, -90.0, 0.0, 0
  },
  {
    DIRECTORY"/grid_surface_90_180_0_90.bin",
    0, 0, 0, 90.0, 180.0, 0.0, 90.0, 0
  }
};

/*=========================== FORWARD DECLARATIONS ==========================*/

static void closeElevationFiles( void );

static int selectElevationFile( float longitude, float latitude );

static int readElevationFile( int index );

static int readElevationFileHeader( FILE* file, int* rows, int* columns,
                                    float* longitudeMinimum,
                                    float* longitudeMaximum,
                                    float* latitudeMinimum,
                                    float* latitudeMaximum );

static void* allocate( size_t bytesEach, size_t count );

static void rotate4ByteArrayIfLittleEndian( void* array, size_t count );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: elevation - Elevation at a longitude-latitude point.
INPUTS:  float longitude  Longitude of point.
         float latitude   Latitude (WGS84) of point.
RETURNS: float elevation in meters above/below mean sea level.
******************************************************************************/

float elevationAt( float longitude, float latitude ) {
  float result = 0.0;
  struct stat unused;
  int index = 0;
  assert( IN_RANGE( longitude, -180.0, 180.0 ) );
  assert( IN_RANGE( latitude, -90.0, 90.0 ) );
  index = selectElevationFile( longitude, latitude );

  if ( stat( elevationFiles[ index ].name, &unused ) == 0 ) { /* File exists.*/

    if ( ! elevationFiles[ index ].file ) {
      readElevationFile( index ); /* Avoid unneeded alloc/read: */

      if ( ! registeredCloseAtExit ) {
        registeredCloseAtExit = atexit( closeElevationFiles ) == 0;
      }
    }

    if ( elevationFiles[ index ].data ) {
      const ElevationFile* const elevationFile = elevationFiles + index;
      const float longitudeFraction =
        ( longitude - elevationFile->longitudeMinimum ) /
        ( elevationFile->longitudeMaximum - elevationFile->longitudeMinimum );
      const float latitudeFraction =
        ( latitude - elevationFile->latitudeMinimum ) /
        ( elevationFile->latitudeMaximum - elevationFile->latitudeMinimum );
      const int columns = elevationFile->columns;
      const int rows = elevationFile->rows;
      const int column0 = longitudeFraction * columns;
      const int row0 = latitudeFraction * rows;
      const int column = CLAMPED_TO_RANGE( column0, 0, columns - 1 );
      const int row = CLAMPED_TO_RANGE( row0, 0, rows - 1 );
      const int dataIndex = ( rows - 1 - row ) * columns + column;
      result = elevationFile->data[ dataIndex ];
      DEBUG( fprintf( stderr, "%s [%d, %d](%f, %f) -> %f\n",
                      elevationFile->name,
                      column, row, longitude, latitude, result ); )
    }
  }

  assert( IN_RANGE( result, -11000.0, 9000.0 ) );
  return result;
}



/*============================= PRIVATE FUNCTIONS ===========================*/



/******************************************************************************
PURPOSE: closeElevationFiles - Close all elevation files.
******************************************************************************/

static void closeElevationFiles( void ) {
  const int count = sizeof elevationFiles / sizeof *elevationFiles;
  int index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( elevationFiles[ index ].file ) {
      fclose( elevationFiles[ index ].file );
      elevationFiles[ index ].file = 0;
    }

    FREE( elevationFiles[ index ].data );
  }

  memset( elevationFiles, 0, sizeof elevationFiles );
}



/******************************************************************************
PURPOSE: selectElevationFile - Return index of elevation file for a given
         longitude-latitude point.
INPUTS:  float longitude  Longitude of point.
         float latitude   Latitude of point.
RETURNS: int index of elevation file for point.
******************************************************************************/

static int selectElevationFile( float longitude, float latitude ) {
  const int count = sizeof elevationFiles / sizeof *elevationFiles;
  int index = 0;
  int result = -1;
  assert( IN_RANGE( longitude, -180.0, 180.0 ) );
  assert( IN_RANGE( latitude, -90.0, 90.0 ) );

  for ( index = 0; AND2( result == -1, index < count ); ++index ) {

    if ( AND2( IN_RANGE( longitude,
                         elevationFiles[ index ].longitudeMinimum,
                         elevationFiles[ index ].longitudeMaximum ),
               IN_RANGE( latitude,
                         elevationFiles[ index ].latitudeMinimum,
                         elevationFiles[ index ].latitudeMaximum ) ) ) {
      result = index;
    }
  }

  assert( IN_RANGE( result, 0,
                    sizeof elevationFiles / sizeof *elevationFiles - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: readElevationFile - Read indexed elevation file.
INPUTS:  int index  Index of file to read.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readElevationFile( int index ) {
  int result = 0;
  assert( IN_RANGE( index, 0,
                    sizeof elevationFiles / sizeof *elevationFiles - 1 ) );
  assert( elevationFiles[ index ].file == 0 );
  assert( elevationFiles[ index ].data == 0 );

  elevationFiles[ index ].file = fopen( elevationFiles[ index ].name, "rb" );

  if ( elevationFiles[ index ].file ) {

    if ( readElevationFileHeader( elevationFiles[ index ].file,
                                  &elevationFiles[ index ].rows,
                                  &elevationFiles[ index ].columns,
                                  &elevationFiles[ index ].longitudeMinimum,
                                  &elevationFiles[ index ].longitudeMaximum,
                                  &elevationFiles[ index ].latitudeMinimum,
                                  &elevationFiles[ index ].latitudeMaximum )) {
      const int count =
        elevationFiles[ index ].rows * elevationFiles[ index ].columns;
      elevationFiles[ index ].data = allocate( sizeof (float), count );

      if ( elevationFiles[ index ].data ) {

        if ( fread( elevationFiles[ index ].data, sizeof (float), count,
                    elevationFiles[ index ].file ) != count ) {
          fprintf( stderr, "Failed to read %d 4-byte float elevations.\n",
                   count );
        } else {
          rotate4ByteArrayIfLittleEndian( elevationFiles[ index ].data, count);
          result = 1;
          DEBUG( fprintf( stderr, "Read elevation file index = %d\n", index );)
        }
      }
    }
  }

  if ( ! result ) {

    if ( elevationFiles[ index ].file ) {
      fclose( elevationFiles[ index ].file );
      elevationFiles[ index ].file = 0;
    }
    
    FREE( elevationFiles[ index ].data );
    memset( elevationFiles + index, 0, sizeof (ElevationFile) );
  }

  assert( IS_BOOL( result ) );
  assert( IMPLIES_ELSE( result,
                        AND8( elevationFiles[ index ].file,
                              elevationFiles[ index ].rows > 0,
                              elevationFiles[ index ].columns > 0,
                              IN_RANGE( elevationFiles[index].longitudeMinimum,
                                        -180.0, 180.0 ),
                              IN_RANGE( elevationFiles[index].longitudeMaximum,
                                        elevationFiles[index].longitudeMinimum,
                                        180.0 ),
                              IN_RANGE( elevationFiles[index].latitudeMinimum,
                                        -90.0, 90.0 ),
                              IN_RANGE( elevationFiles[index].latitudeMaximum,
                                        elevationFiles[index].latitudeMinimum,
                                        90.0 ),
                              elevationFiles[ index ].data ),
                        IS_ZERO8( elevationFiles[ index ].file,
                                  elevationFiles[ index ].rows,
                                  elevationFiles[ index ].columns,
                                  elevationFiles[ index ].longitudeMinimum,
                                  elevationFiles[ index ].longitudeMaximum,
                                  elevationFiles[ index ].latitudeMinimum,
                                  elevationFiles[ index ].latitudeMaximum,
                                  elevationFiles[ index ].data ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readElevationFileHeader - Read 4 line ASCII header of elevation file.
OUTPUTS: FILE* file     File to read.
         int* rows      Number of grid rows.
         int* columns   Number of grid columns.
         float* longitudeMinimum  Longitude minimum of grid domain.
         float* longitudeMaximum  Longitude maximum of grid domain.
         float* latitudeMinimum  Latitude minimum of grid domain.
         float* latitudeMaximum  Latitude maximum of grid domain.
RETURNS: int 1 if successful, else 0 and a failure message is printed.
NOTES: Header line looks like:
Content-type: application/octet-stream; charset=iso-8859-1
# Dimensions: rows columns lonmin lonmax latmin latmax
780        598  -84.4962262273  -66.810974227   24.4000371   47.467757100
# IEEE-754 32-bit float data[rows][columns]:
******************************************************************************/

static int readElevationFileHeader( FILE* file, int* rows, int* columns,
                                    float* longitudeMinimum,
                                    float* longitudeMaximum,
                                    float* latitudeMinimum,
                                    float* latitudeMaximum ) {

  int result = 0;
  char typeString[ 16 ] = "";
  memset( typeString, 0, sizeof typeString );
  assert( rows ); assert( columns );
  assert( longitudeMinimum ); assert( longitudeMaximum );
  assert( latitudeMinimum ); assert( latitudeMaximum );

  *rows = *columns = 0;
  *longitudeMinimum = *longitudeMaximum = 0.0;
  *latitudeMaximum = *latitudeMaximum = 0.0;

  result = fscanf( file,
                   "%*[^\n]\n"
                   "%*[^\n]\n"
                   "%d %d %f %f %f %f\n"
                   "%*s %15s %*[^\n]%*c",
                   rows, columns,
                   longitudeMinimum, longitudeMaximum,
                   latitudeMinimum, latitudeMaximum,
                   typeString ) == 7;

  if ( ! result ) {
    perror( "Failed to read 4-line ASCII header because" );
  } else {
    result = AND7( *rows > 0,
                   *columns > 0,
                   IN_RANGE( *longitudeMinimum, -180.0, 180.0 ),
                   IN_RANGE( *latitudeMinimum, -90.0, 90.0 ),
                   IN_RANGE( *longitudeMaximum, *longitudeMinimum, 180.0 ),
                   IN_RANGE( *latitudeMaximum, *longitudeMinimum, 90.0 ),
                   OR2( ! strcmp( typeString, "signed" ),
                        ! strcmp( typeString, "IEEE-754" ) ) );

    if ( ! result ) {
      fprintf( stderr,
               "Read invalid 4-line ASCII header:\n"
               "rows %d columns "
               "%d lonmin %f latmin %f lonmax %f latmax %f type %s\n",
               *rows, *columns,
               *longitudeMinimum, *latitudeMinimum,
               *longitudeMaximum, *latitudeMaximum,
               typeString );
      *rows = *columns = 0;
      *longitudeMinimum = *longitudeMaximum = 0.0;
      *latitudeMaximum = *latitudeMaximum = 0.0;
    }
  }

  assert( IS_BOOL( result ) );
  assert( IMPLIES_ELSE( result,
                        AND6( *rows > 0,
                              *columns > 0,
                              IN_RANGE( *longitudeMinimum,
                                        -180.0, 180.0 ),
                              IN_RANGE( *latitudeMinimum,
                                        -90.0, 90.0 ),
                              IN_RANGE( *longitudeMaximum,
                                        *longitudeMinimum, 180.0 ),
                              IN_RANGE( *latitudeMaximum,
                                        *latitudeMinimum, 90.0 ) ),
                        IS_ZERO6( *rows,
                                  *columns,
                                  *longitudeMinimum,
                                  *longitudeMaximum,
                                  *latitudeMinimum,
                                  *latitudeMaximum ) ) );
  return result;
}



/******************************************************************************
PURPOSE: allocate - Allocates memory by calling malloc() to allocate
         (then zero) an array of count items, each of size bytesEach.
         If allocation fails then a message is printed to stderr.
INPUTS:  size_t sizeEach  Number of bytes per item.
         size_t count     Number of items to allocate.
RETURNS: void*  The resulting address (zeroed), or 0 if unsuccessful.
******************************************************************************/

static void* allocate( size_t bytesEach, size_t count ) {
  void* result = 0;
  assert( bytesEach > 0 ); assert( count > 0 );

  if ( bytesEach > ULONG_MAX / count || count > ULONG_MAX / bytesEach ) {
    fprintf( stderr, "\a\nCannot allocate %lu x %lu bytes ", count, bytesEach);
  } else {
    const size_t bytes = bytesEach * count;
    result = malloc( bytes );

    if ( ! result ) {
      fprintf(stderr, "\a\nCannot allocate %lu x %lu bytes ", count,bytesEach);
      perror( "because" );
    } else {
      memset( result, 0, bytes );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: rotate4ByteArrayIfLittleEndian() - Rotate 4-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array    Array of 4-byte values to rotate.
         size_t count   Number of items in array.
OUTPUTS: void* array    Array of rotated values.
******************************************************************************/

static void rotate4ByteArrayIfLittleEndian( void* array, size_t count ) {

#if IS_LITTLE_ENDIAN

  int* const array4 = array;
  size_t index = 0;
  assert_static( sizeof (int) == 4 );
  assert( array ); assert( count > 0 );

  for ( index = 0; index < count; ++index ) {
    const int value = array4[ index ];
    const int newValue =
      ( value & 0xff000000 ) >> 24 |
      ( value & 0x00ff0000 ) >>  8 |
      ( value & 0x0000ff00 ) <<  8 |
      ( value & 0x000000ff ) << 24;
    array4[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}

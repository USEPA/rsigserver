
/******************************************************************************
PURPOSE: GridSubset.c - Subset and aggregate binary grid surface (*.bin) files.
NOTES:   Row order is north-to-south (like ESRI ASCII Grid files).
         Compile:
           cc -DNDEBUG -o GridSubset GridSubset.c
         Run:
           GridSubset input.bin \
                      [-time first_timestep/stamp last_timestep/stamp ] \
                      [-subset lonmin latmin lonmax latmax ] \
                      [-aggregate size mean | mode ] \
             > output.bin
         Example:
           GridSubset grid_surface_nlcd2001_gulf.bin \
                      -subset -90 28 -85 32 -aggregate 2048 mode \
                      > subset.bin
           head -4 subset.bin
HISTORY: 2011-02-14 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For FILE, stderr, fscanf(), fopen(), fseek(). */
#include <stdlib.h> /* For malloc(), free(), atof(). */
#include <string.h> /* For memset(). */
#include <limits.h> /* For ULONG_MAX. */

/*================================== MACROS =================================*/

/* Compile-time assertion: */

#define assert_static(a) extern int unused_assert_static_[ (a) ? 1 : -1 ]

/* Logic macros: */

#define IS_BOOL(c) ((c)==0||(c)==1)
#define IMPLIES(p,c) (!(p)||(c))
#define IMPLIES_ELSE(p,c1,c2) (((p)&&(c1))||((!(p))&&(c2)))
#define IN3(x,a,b) ((x)==(a)||(x)==(b))
#define IN7(x,a,b,c,d,e,f) \
  ((x)==(a)||(x)==(b)||(x)==(c)||(x)==(d)||(x)==(e)||(x)==(f))
#define OR2(a,b) ((a)||(b))
#define OR3(a,b,c) ((a)||(b)||(c))
#define OR4(a,b,c,d) ((a)||(b)||(c)||(d))
#define AND2(a,b) ((a)&&(b))
#define AND3(a,b,c) ((a)&&(b)&&(c))
#define AND4(a,b,c,d) ((a)&&(b)&&(c)&&(d))
#define AND5(a,b,c,d,e) ((a)&&(b)&&(c)&&(d)&&(e))
#define AND6(a,b,c,d,e,f) ((a)&&(b)&&(c)&&(d)&&(e)&&(f))
#define AND7(a,b,c,d,e,f,g) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g))
#define AND8(a,b,c,d,e,f,g,h) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h))
#define AND9(a,b,c,d,e,f,g,h,i) ((a)&&(b)&&(c)&&(d)&&(e)&&(f)&&(g)&&(h)&&(i))

#define IS_ZERO6(a,b,c,d,e,f) \
  ((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0)

#define IS_ZERO7(a,b,c,d,e,f,g) \
  ((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0)

#define IS_ZERO8(a,b,c,d,e,f,g,h) \
  ((a)==0&&(b)==0&&(c)==0&&(d)==0&&(e)==0&&(f)==0&&(g)==0&&(h)==0)

#define IN_RANGE( value, minimum , maximum ) \
  ( (value) >= (minimum) && (value) <= (maximum) )

#define CLAMPED_TO_RANGE( value, low, high ) \
((value) < (low) ? (low) : (value) > (high) ? (high) : (value))

#ifndef MIN
#define MIN(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

/* Optionally include debug statements if -DDEBUGGING: */

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(s)
#endif

/* Optionally include debug statements if -DDEBUGGING2: */

#ifdef DEBUGGING2
#define DEBUG2(s) s
#else
#define DEBUG2(s)
#endif

/* Memory allocation/deallocation: */

#define FREE( p ) ( ( (p) ? free(p) : (void) 0 ), (p) = 0 )

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

/*=================================== TYPES =================================*/

#define F_MISSING (-9999.0)
enum { MISSING = -99 };          /* Missing data value. */
enum { MEAN, MODE };             /* Grid cell data aggregation methods. */
enum { MINIMUM, MAXIMUM };
enum { LONGITUDE, LATITUDE };
enum { COLUMN, ROW };
enum { FIRST, LAST };            /* 0-based indices. */
enum { MAXIMUM_TIMESTAMPS = 24 * 366 }; /* Hours in a leap year. */
typedef char LongName[ 63 + 1 ];
typedef char Units[ 15 + 1 ];
typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE LATITUDE ][MINIMUM MAXIMUM].*/
typedef int Range[ 2 ][ 2 ];     /* [ COLUMN ROW ][ FIRST LAST ]. */

enum { BYTE_TYPE, UINT16_TYPE, FLOAT_TYPE, DATA_TYPES }; /* Cell value type. */

#define IS_VALID_TYPE( type ) IN_RANGE( type, 0, DATA_TYPES - 1 )

#define WORD_SIZE(type) \
  ((type) == FLOAT_TYPE ? sizeof (float) \
   : (type) == UINT16_TYPE ? sizeof (short) : sizeof (char))

/*=========================== FORWARD DECLARATIONS ==========================*/


static void usage( const char* program );

static int parseArguments( int argc, const char* argv[],
                           const char** inputFileName,
                           int timestepSubset[ 2 ],
                           Bounds subset, int* dimension, int* method );

static int isValidArgs( int argc, const char* argv[] );

static int readHeader( FILE* file, int* timesteps, int* rows, int* columns,
                       Bounds bounds, int* type,
                       LongName name, Units units, int timestamps[] );

static void convertTimestepSubset( int timestepSubset[ 2 ],
                                   const int timesteps,
                                   const int timestamps[] );

static int computeSubset( int rows, int columns, Bounds domain,
                          Bounds subset, int dimension, int* size,
                          int* subsetRows, int* subsetColumns, Range range );

static void subsetIndices( const double clip[ 2 ],
                           double range[ 2 ], int count, int indices[ 2 ] );

static void computeStride( int rows, int columns, int dimension,
                           int* subsetRows, int* subsetColumns, int* size,
                           Range range, Bounds subset );

static void aggregate( int size, int columns, int firstColumn,
                       int outputColumns, int method, int type,
                       const void* input, void* output );

static int indexOfMaximum( const int array[], int count );

static void* allocate( size_t bytesEach, size_t count );

static int isValidBounds( const Bounds bounds );

static void rotate4ByteArrayIfLittleEndian( void* array, long long count );

static void rotate2ByteArrayIfLittleEndian( void* array, long long count );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: main - Main routine.
INPUTS:  int argc     Number of command-line arguments.
         char* argv[]  Command-line arguments.
RETURNS: int 0 if successful, else 1 and failure message is printed.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int ok = 0;
  const char* inputFileName = 0;
  int timestepSubset[ 2 ] = { MISSING, MISSING };
  int dimension = 0;
  int method = 0;
  Bounds subset = { { 0.0, 0.0 }, { 0.0, 0.0 } };

  if ( parseArguments( argc, (const char**) argv, &inputFileName,
                       timestepSubset, subset, &dimension, &method ) ) {
    FILE* inputFile = fopen( inputFileName, "rb" );

    if ( ! inputFile ) {
      fprintf( stderr, "Failed to open input file %s ", inputFileName );
      perror( "because" );
    } else {
      int type = 0;
      int timesteps = 0;
      int rows = 0;
      int columns = 0;
      Bounds domain = { { 0.0, 0.0 }, { 0.0, 0.0 } };
      LongName name = "";
      Units units = "";
      int timestamps[ MAXIMUM_TIMESTAMPS ];
      memset( timestamps, 0, sizeof timestamps );
      memset( name, 0, sizeof name );
      memset( units, 0, sizeof units );

      if ( readHeader( inputFile, &timesteps, &rows, &columns, domain, &type,
                       name, units, timestamps ) ) {
        int subsetRows = 0;
        int subsetColumns = 0;
        int size = 0;
        Range range = { { 0, 0 }, { 0, 0 } };

        if ( type == FLOAT_TYPE ) {
          method = MEAN;
        }

        if ( timestepSubset[ FIRST ] == MISSING ) { /* If unspecified, */
          timestepSubset[ FIRST ] = 0;
          timestepSubset[ LAST  ] = timesteps - 1; /* then use full range. */
        } else if ( timestepSubset[ FIRST ] > 1900010100 ) { /* YYYYMMDDHH: */
          convertTimestepSubset( timestepSubset, timesteps, timestamps );
        }

        DEBUG( fprintf( stderr, "timestepSubset[] = [%d %d], dimension = %d\n",
                        timestepSubset[ FIRST ], timestepSubset[ LAST ],
                        dimension ); )

        if ( AND3( IN_RANGE( timestepSubset[ FIRST ], 0, timesteps - 1 ),
                   IN_RANGE( timestepSubset[ LAST ], timestepSubset[ FIRST ],
                             timesteps - 1 ),
                   computeSubset( rows, columns, domain, subset, dimension,
                             &size, &subsetRows, &subsetColumns, range ) ) ) {
          const int wordSize = WORD_SIZE( type );
          const long timestepSkipBytes =
            (const long) timestepSubset[ FIRST ] * rows * columns * wordSize;
          const long skipRowBytes =
            (const long) range[ ROW ][ FIRST ] * columns * wordSize;
          timesteps = 1 + timestepSubset[ LAST ] - timestepSubset[ FIRST ];

          DEBUG( fprintf( stderr, "size = %d\n", size ); )
          DEBUG( fprintf( stderr, "wordSize = %d\n", wordSize ); )
          DEBUG( fprintf( stderr, "columns = %d\n", columns ); )
          DEBUG( fprintf( stderr, "range[ ROW ][ FIRST ] = %d\n",
                          range[ ROW ][ FIRST ] ); )
          DEBUG( fprintf( stderr, "timestepSkipBytes = %ld\n",
                          timestepSkipBytes ); )
          DEBUG( fprintf( stderr, "skipRowBytes = %ld\n", skipRowBytes ); )
          assert( skipRowBytes >= 0 );

          if ( fseek( inputFile, timestepSkipBytes + skipRowBytes, SEEK_CUR)) {
            fprintf( stderr, "Failed to seek %ld bytes into file %s ",
                     timestepSkipBytes + skipRowBytes, inputFileName );
            perror( "because" );
          } else {
            const size_t inputDataSize = (size_t) size * columns;
            const size_t outputDataSize = (size_t) subsetColumns;
            const size_t dataSize = inputDataSize + outputDataSize;
            signed char* data = allocate( wordSize, dataSize );

            DEBUG( fprintf( stderr,
                            "inputDataSize = %lu, outputDataSize = %lu\n",
                            inputDataSize, outputDataSize ); )

            if ( data ) {
              const long remainingRows = rows - ( range[ ROW ][ LAST ] + size);
              const long remainingRowBytes =
                (const long) remainingRows * columns * wordSize;
              const long seekRowBytes =
                (const long) remainingRowBytes + skipRowBytes;
              signed char* const outputData =
                data + inputDataSize * wordSize / sizeof (char);
              const int firstColumn = range[ COLUMN ][ FIRST ];
              int row = 0;
              int timestep = 0;
              const char* const timestepsDimension =
                timesteps > 1 ? "[timesteps]" : "";
              char timestepsValue[ 64 ] = "";
              memset( timestepsValue, 0, sizeof timestepsValue );
              DEBUG( fprintf( stderr, "seekRowBytes = %ld\n", seekRowBytes );)

              if ( timesteps > 1 ) {
                snprintf( timestepsValue,
                          sizeof timestepsValue / sizeof timestepsValue[ 0 ],
                          "%10d", timesteps );
              }

              /* Write ASCII header: */

              printf( "Content-type: application/octet-stream; "
                      "charset=iso-8859-1\n" ); /* For streaming over web.*/

              if ( *name ) { /* Write variable name units of 7-line header: */
                printf( "# variable units:\n%s %s\n", name, units );
              }

              printf( "# Dimensions: "
                      "%srows columns lonmin lonmax latmin latmax\n",
                      timesteps > 1 ? "timesteps " : "" );
              printf( "%s%10d %10d %22.17lf %22.17lf %22.17lf %22.17lf\n",
                      timestepsValue,
                      subsetRows, subsetColumns,
                      subset[ LONGITUDE ][ MINIMUM ],
                      subset[ LONGITUDE ][ MAXIMUM ],
                      subset[ LATITUDE  ][ MINIMUM ],
                      subset[ LATITUDE  ][ MAXIMUM ] );

              if ( timestamps[ 0 ] ) { /* For 7-line header: */
                printf( "# char yyyymmddhh[timesteps][1]] and\n" );
              }

              if ( type == FLOAT_TYPE ) {
                printf( "# IEEE-754 32-bit float data%s[rows][columns]:\n",
                        timestepsDimension);
              } else if ( type == UINT16_TYPE ) {
                printf( "# MSB 16-bit uint16 data%s[rows][columns]:\n",
                        timestepsDimension );
              } else {
                printf( "# signed char data%s[rows][columns]:\n",
                        timestepsDimension );
              }

              if ( timestamps[ 0 ] ) { /* Write timestamps: */

                for ( timestep = timestepSubset[ FIRST ];
                      timestep <= timestepSubset[ LAST ];
                      ++timestep ) {
                  printf( "%d\n", timestamps[ timestep ] );
                }
              }

              ok = 1;

              for ( timestep = 0; AND2(ok, timestep < timesteps); ++timestep) {

                for ( row = range[ ROW ][ FIRST ];
                      AND2( ok, row <= range[ ROW ][ LAST ] );
                      row += size ) {
                  DEBUG( fprintf( stderr,
                                 "Processing timestep %d, rows [%6d %6d]...\n",
                                  timestep, row, row + size - 1 ); )
                  ok = fread( data, inputDataSize * wordSize, 1, inputFile )
                       == 1;

                  if ( ! ok ) {
                    fprintf( stderr, "Failed to read %lu bytes of row data "
                             "from file %s ", inputDataSize, inputFileName );
                    perror( "because" );
                  } else {

                    if ( type == FLOAT_TYPE ) {
                      rotate4ByteArrayIfLittleEndian( data, inputDataSize );
                    } else if ( type == UINT16_TYPE ) {
                      rotate2ByteArrayIfLittleEndian( data, inputDataSize );
                    }

#ifdef DEBUGGING
                    if ( type == FLOAT_TYPE ) {
                      const float* const fdata = (const float*) data;
                      fprintf( stderr, "input data: %f %f %f ... %f\n",
                               fdata[ 0 ],
                               inputDataSize >= 2 ? fdata[ 1 ] : 0.0,
                               inputDataSize >= 3 ? fdata[ 2 ] : 0.0,
                               fdata[ inputDataSize - 1 ] );
                    } else if ( type == UINT16_TYPE ) {
                      const unsigned short* sdata =
                        (const unsigned short*) data;
                      fprintf( stderr, "input data: %d %d %d ... %d\n",
                               sdata[ 0 ],
                               inputDataSize >= 2 ? sdata[ 1 ] : 0,
                               inputDataSize >= 3 ? sdata[ 2 ] : 0,
                               sdata[ inputDataSize - 1 ] );
                    } else {
                      fprintf( stderr, "input data: %d %d %d ... %d\n",
                               data[ 0 ],
                               inputDataSize >= 2 ? data[ 1 ] : 0,
                               inputDataSize >= 3 ? data[ 2 ] : 0,
                               data[ inputDataSize - 1 ] );
                    }
#endif

                    aggregate( size, columns, firstColumn, subsetColumns,
                               method, type, data, outputData );

                    if ( type == FLOAT_TYPE ) {
                      rotate4ByteArrayIfLittleEndian( outputData,
                                                      outputDataSize );
                    } else if ( type == UINT16_TYPE ) {
                      rotate2ByteArrayIfLittleEndian( outputData,
                                                      outputDataSize );
                    }

                    ok = fwrite( outputData, outputDataSize * wordSize, 1,
                                 stdout ) == 1;
                  }

                  if ( ! ok ) {
                    perror( "Failed to write row data because" );
                  }
                } /* Next row. */

                if ( AND2( ok, timestep + 1 < timesteps ) ) {
                  ok = fseek( inputFile, seekRowBytes, SEEK_CUR ) == 0;
        
                  if ( ! ok ) {
                    fprintf( stderr,
                             "Failed to seek %ld bytes into file %s ",
                             seekRowBytes, inputFileName );
                    perror( "because" );                        
                  }
                }
              } /* Next timestep. */

              FREE( data );
            }
          }
        }
      }

      fclose( inputFile ), inputFile = 0;
    }
  }

  return ! ok;
}



/*============================= PRIVATE FUNCTIONS ===========================*/



/******************************************************************************
PURPOSE: usage - Print program usage.
OUTPUTS: const char* program  Name of program.
******************************************************************************/

static void usage( const char* program ) {
  assert( program );
  fprintf( stderr, "\n%s - Subset and aggregate binary grid surface "
           "(*.bin) files.\n", program );
  fprintf( stderr, "usage: %s input.bin ", program );
  fprintf( stderr, "[-time first_timestep/stamp last_timestep/stamp] \\\n" );
  fprintf( stderr, "[-subset lonmin latmin lonmax latmax ] \\\n" );
  fprintf( stderr, "[-aggregate size mean | mode ] > output.bin\n" );
  fprintf( stderr, "example: %s grid_surface_nlcd2001_gulf.bin \\\n", program);
  fprintf( stderr, "-subset -90 28 -85 32 -aggregate 2048 mode \\\n" );
  fprintf( stderr, "> subset.bin\n" );
  fprintf( stderr, "head -4 subset.bin\n" );
  fprintf( stderr, "Notes:\n" );
  fprintf( stderr, "Row order is north-to-south like ASCII Grids.\n\n" );
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
OUTPUTS: int argc                    Number of command-line arguments.
         const char argv[]           Command-line arguments.
OUTPUTS: const char** inputFileName  Name of grid file to read.
         int timestepSubset[ 2 ]     0-based first and last timesteps or
                                     yyyymmddhh format first and last
                                     timestamps to read.
         Bounds subset               Subset to apply to grid file.
         int* dimension              Maximum result dimension.
         int* method                 Aggregation method: MEAN or MODE.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int parseArguments( int argc, const char* argv[],
                           const char** inputFileName,
                           int timestepSubset[ 2 ],
                           Bounds subset, int* dimension, int* method ) {
  int result = 0;
  assert( timestepSubset );
  assert( subset ); assert( dimension ); assert( method );

  /* Initialize outputs: */

  timestepSubset[ FIRST ] = timestepSubset[ LAST ] = MISSING;
  subset[ LONGITUDE ][ MINIMUM ] = -180.0;
  subset[ LONGITUDE ][ MAXIMUM ] =  180.0;
  subset[ LATITUDE  ][ MINIMUM ] =  -90.0;
  subset[ LATITUDE  ][ MAXIMUM ] =   90.0;
  *dimension = INT_MAX;
  *method = MEAN;

  if ( AND2( IN7( argc, 2, 5, 7, 8, 10, 13 ),
             isValidArgs( argc, argv ) ) ) {
    int arg = 2;
    *inputFileName = argv[ 1 ];
    result = 1;

    while ( AND2( result, arg < argc ) ) {

      if ( AND2( ! strcmp( argv[ arg ], "-time" ), arg + 2 < argc ) ) {
        timestepSubset[ FIRST ] = atoi( argv[ arg + 1 ] );
        timestepSubset[ LAST  ] = atoi( argv[ arg + 2 ] );
        result = AND2( timestepSubset[ FIRST ] >= 0,
                       timestepSubset[ LAST ] >= timestepSubset[ FIRST ] );

        if ( ! result ) {
          fprintf( stderr, "Invalid -time arguments.\n" );
        } else {
          arg += 3;
        }
      } else if ( AND2( ! strcmp( argv[ arg ], "-subset" ), arg + 4 < argc )) {
        subset[ LONGITUDE ][ MINIMUM ] = atof( argv[ arg + 1 ] );
        subset[ LATITUDE  ][ MINIMUM ] = atof( argv[ arg + 2 ] );
        subset[ LONGITUDE ][ MAXIMUM ] = atof( argv[ arg + 3 ] );
        subset[ LATITUDE  ][ MAXIMUM ] = atof( argv[ arg + 4 ] );
        result = isValidBounds( (const double(*)[2]) subset );

        if ( ! result ) {
          fprintf( stderr, "Invalid -subset arguments.\n" );
        } else {
          arg += 5;
        }
      } else if ( AND2( ! strcmp( argv[ arg ], "-aggregate" ),
                        arg + 2 < argc ) ) {
        *dimension = atoi( argv[ arg + 1 ] );
        *method = -1;

        if ( ! strcmp( argv[ arg + 2 ], "mean" ) ) {
          *method = MEAN;
        } else if ( ! strcmp( argv[ arg + 2 ], "mode" ) ) {
          *method = MODE;
        }

        result = AND2( *dimension > 0, IN3( *method, MEAN, MODE ) );

        if ( ! result ) {
          fprintf( stderr, "Invalid -aggregate arguments.\n" );
        } else {
          arg += 3;
        }
      } else {
        result = 0;
        fprintf( stderr, "Invalid/incomplete arguments.\n" );
      }
    }
  }

  if ( ! result ) {
    timestepSubset[ FIRST ] = timestepSubset[ LAST ] = 0;
    subset[ LONGITUDE ][ MINIMUM ] = 0.0;
    subset[ LONGITUDE ][ MAXIMUM ] = 0.0;
    subset[ LATITUDE  ][ MINIMUM ] = 0.0;
    subset[ LATITUDE  ][ MAXIMUM ] = 0.0;
    *dimension = 0;
    *method = 0;
    usage( argv[ 0 ] );
  }

  assert( IS_BOOL( result ) );
  assert( IMPLIES_ELSE( result,
                        AND4( isValidBounds( (const double (*)[2]) subset ),
                              IMPLIES_ELSE( timestepSubset[ FIRST ] == MISSING,
                                            timestepSubset[ LAST  ] == MISSING,
                                            AND2( timestepSubset[ FIRST ] >= 0,
                                                  timestepSubset[ LAST  ] >=
                                                  timestepSubset[ FIRST ] ) ),
                              *dimension >= 0,
                              IN3( *method, MEAN, MODE ) ),
                        IS_ZERO8( timestepSubset[ FIRST ],
                                  timestepSubset[ LAST ],
                                  *dimension, *method,
                                  subset[ 0 ][ 0 ], subset [ 0 ][ 1 ],
                                  subset[ 1 ][ 0 ], subset [ 1 ][ 1 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidArgs() - Are each of the strings non-zero and non-zero-length?
INPUTS:  int argc                  Number of strings to check.
         const char* argv[ argc ]  Strings to check.
RETURNS: int 1 if non-zero, else 0.
******************************************************************************/

static int isValidArgs( int argc, const char* argv[] ) {
  int result = 0;

  if ( AND2( IN_RANGE( argc, 1, INT_MAX - 1 ), argv != 0 ) ) {
    int index = 0;

    do {

      if ( OR2( argv[ index ] == 0, argv[ index ][ 0 ] == '\0' ) ) {
        index = argc;
      }

      ++index;
    } while ( index < argc );

    result = index == argc;
  }

  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readHeader - Read 4 line ASCII header of binary Grid file from stdin.
OUTPUTS: FILE* file     File to read.
         int* timesteps Number of timesteps.
         int* rows      Number of grid rows.
         int* columns   Number of grid columns.
         Bounds domain  Domain of grid.
         int* type      Type of cell value: BYTE_TYPE, FLOAT_TYPE, etc.
         LongName name  Variable name or "".
         Units    units Variable units or "".
         int timestamps[ MAXIMUM_TIMESTAMPS ] yyyymmddhh or 0.
RETURNS: int 1 if successful, else 0 and a failure message is printed.
NOTES: Header line looks like:
Content-type: application/octet-stream; charset=iso-8859-1
# Dimensions: rows columns lonmin lonmax latmin latmax
780        598  -84.4962262273  -66.810974227   24.4000371   47.467757100
# IEEE-754 32-bit float data[rows][columns]:
******************************************************************************/

static int readHeader( FILE* file, int* timesteps, int* rows, int* columns,
                       Bounds domain, int* type, LongName name, Units units,
                       int timestamps[ MAXIMUM_TIMESTAMPS ] ) {

  int result = 0;
  int hasTimestep = 0;
  char typeString[ 16 ] = "";
  char line[ 256 ] = "";
  memset( typeString, 0, sizeof typeString );
  memset( line, 0, sizeof line );
  assert( rows ); assert( columns );
  assert( domain ); assert( type );
  assert( name ); assert( units ); assert( timestamps );

  *timesteps = *rows = *columns = *type = 0;
  memset( domain, 0, sizeof (Bounds) );
  *name = *units = '\0';
  timestamps[ 0 ] = 0;

  if ( fgets( line, sizeof line / sizeof *line, file ) ) {

    if ( fgets( line, sizeof line / sizeof *line, file ) ) {
      const int hasVariable = strstr( line, "# variable units:" ) != 0;
      result = 1;

      if ( hasVariable ) {
        result = fscanf( file, "%63s %15s\n", name, units ) == 2;
        name[ sizeof (LongName) / sizeof *name  - 1 ] = '\0';
        units[ sizeof (Units)   / sizeof *units - 1 ] = '\0';

        if ( AND3( result, *name, *units ) ) { /* Read dimensions line: */
          result = fgets( line, sizeof line / sizeof *line, file ) != 0;
        }
      }

      if ( result ) {
        hasTimestep = strstr( line, ": timesteps rows " ) != 0;

        if ( hasTimestep ) {
          result = fscanf( file,
                      "%d %d %d %lf %lf %lf %lf\n"
                      "%*s %15s %*[^\n]%*c",
                      timesteps, rows, columns,
                      &domain[ LONGITUDE ][ MINIMUM ],
                      &domain[ LONGITUDE ][ MAXIMUM ],
                      &domain[ LATITUDE  ][ MINIMUM ],
                      &domain[ LATITUDE  ][ MAXIMUM ],
                      typeString ) == 8;
        } else {
          *timesteps = 1;
          result = fscanf( file,
                      "%d %d %lf %lf %lf %lf\n"
                      "%*s %15s %*[^\n]%*c",
                      rows, columns,
                      &domain[ LONGITUDE ][ MINIMUM ],
                      &domain[ LONGITUDE ][ MAXIMUM ],
                      &domain[ LATITUDE  ][ MINIMUM ],
                      &domain[ LATITUDE  ][ MAXIMUM ],
                      typeString ) == 7;
        }
      }

      if ( AND3( result, *name, ! strcmp( typeString, "char" ) ) ) {
        int timestep = 0;
        result = fscanf( file, "%*s %15s %*[^\n]%*c", typeString ) == 1;

        /* Read yyyymmddhh array: */

        for ( timestep = 0; AND2( result, timestep < *timesteps); ++timestep) {
          int yyyymmddhh = 0;
          result = fscanf( file, "%d\n", &yyyymmddhh ) == 1;
          result = AND2( result,
                         IN_RANGE( yyyymmddhh, 1970010100, 2100123123 ) );
          timestamps[ timestep ] = yyyymmddhh;
        }
      }
    }
  }

  DEBUG( fprintf( stderr, "timesteps = %d\n", *timesteps ); )
  DEBUG( fprintf( stderr, "rows = %d\n", *rows ); )
  DEBUG( fprintf( stderr, "columns = %d\n", *columns ); )
  DEBUG( fprintf( stderr, "longitudeMinimum = %28.17lf\n",
                  domain[ LONGITUDE ][ MINIMUM ] ); )
  DEBUG( fprintf( stderr, "longitudeMaximum = %28.17lf\n",
                  domain[ LONGITUDE ][ MAXIMUM ] ); )
  DEBUG( fprintf( stderr, "latitudeMinimum  = %28.17lf\n",
                  domain[ LATITUDE ][ MINIMUM ] ); )
  DEBUG( fprintf( stderr, "latitudeMaximum  = %28.17lf\n",
                  domain[ LATITUDE ][ MAXIMUM ] ); )
  DEBUG( fprintf( stderr, "type = %d\n", *type ); )
  DEBUG( fprintf( stderr, "name = '%s'\n", name ); )  
  DEBUG( fprintf( stderr, "units = '%s'\n", units ); )  
  DEBUG( fprintf( stderr, "timestamps[ 0 ] = %d\n", timestamps[ 0 ] ); )  

  if ( ! result ) {
    perror( "Failed to read ASCII header because" );
  } else {
    result = AND8( *timesteps > 0,
                   *rows > 0,
                   *columns > 0,
                   IN_RANGE( domain[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ),
                   IN_RANGE( domain[ LATITUDE  ][ MINIMUM ], -90.0, 90.0 ),
                   IN_RANGE( domain[ LONGITUDE ][ MAXIMUM ],
                             domain[ LONGITUDE ][ MINIMUM ], 180.0 ),
                   IN_RANGE( domain[ LATITUDE  ][ MAXIMUM ],
                             domain[ LATITUDE  ][ MINIMUM ], 90.0 ),
                   OR3( ! strcmp( typeString, "signed" ),
                        ! strcmp( typeString, "MSB" ),
                        ! strcmp( typeString, "IEEE-754" ) ) );

    if ( ! result ) {
      fprintf( stderr,
               "Read invalid ASCII header:\n"
               "timesteps %d rows %d columns %d "
               "lonmin %lf latmin %lf lonmax %lf latmax %lf type %s"
               "timestamps[0] = %d\n",
               *timesteps, *rows, *columns,
               domain[ LONGITUDE ][ MINIMUM ], domain[ LATITUDE ][ MINIMUM ],
               domain[ LONGITUDE ][ MAXIMUM ], domain[ LATITUDE ][ MAXIMUM ],
               typeString, timestamps[ 0 ] );
      *timesteps = *rows = *columns = *type = 0;
      memset( domain, 0, sizeof (Bounds) );
      *name = *units = '\0';
      timestamps[ 0 ] = 0;
    }
  }

  *type =
    ! strcmp( typeString, "IEEE-754" ) ? FLOAT_TYPE
    : ! strcmp( typeString, "MSB" ) ? UINT16_TYPE
    : BYTE_TYPE;
  DEBUG( fprintf( stderr, "result = %d\n", result ); )

  assert( IS_BOOL( result ) );
  assert( IS_VALID_TYPE( *type ) );
  assert( IMPLIES_ELSE( result,
                        AND8( *timesteps > 0,
                              *rows > 0,
                              *columns > 0,
                              IN_RANGE( domain[ LONGITUDE ][ MINIMUM ],
                                        -180.0, 180.0 ),
                              IN_RANGE( domain[ LATITUDE ][ MINIMUM ],
                                        -90.0, 90.0 ),
                              IN_RANGE( domain[ LONGITUDE ][ MAXIMUM ],
                                        domain[ LONGITUDE ][ MINIMUM ],
                                        180.0 ),
                              IN_RANGE( domain[ LATITUDE ][ MAXIMUM ],
                                        domain[ LATITUDE ][ MINIMUM ],
                                        90.0 ),
                              IMPLIES_ELSE( *name, *units, ! *units ) ),
                        IS_ZERO7( *timesteps, *rows, *columns,
                                  domain[ LONGITUDE ][ MINIMUM ],
                                  domain[ LONGITUDE ][ MAXIMUM ],
                                  domain[ LATITUDE ][ MINIMUM ],
                                  domain[ LATITUDE ][ MAXIMUM ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: convertTimestepSubset - Convert yyyymmddhh format first and last
         subset timestamps to 0-based timestep indices.
INPUTS:  int timestepSubset[ 2 ]            First, last yyyymmddhh timestamps.
         const int timesteps                Number of timesteps.
         const int timestamps[ timesteps ]  Increasing yyyymmddhh timestamps.
OUTPUTS: int timestepSubset[ 2 ]            First, last 0-based timesteps
                                            else MISSING if not in time range.
******************************************************************************/

static void convertTimestepSubset( int timestepSubset[ 2 ],
                                   const int timesteps,
                                   const int timestamps[] ) {

  assert( timestepSubset ); assert( timesteps > 0 ); assert( timestamps );
  assert( IN_RANGE( timestepSubset[ FIRST ], 1970010100, 2100123123 ) );
  assert( IN_RANGE( timestepSubset[ LAST  ], timestepSubset[ FIRST ],
                    2100123123 ) );
  assert( IN_RANGE( timestamps[ 0 ], 1970010100, 2100123123 ) );
  assert( IN_RANGE( timestamps[ timesteps - 1 ], timestamps[ 0 ],
                    2100123123 ) );

  if ( OR2( timestepSubset[ LAST  ] < timestamps[ 0 ],
            timestepSubset[ FIRST ] > timestamps[ timesteps - 1 ] ) ) {
    timestepSubset[ FIRST ] = timestepSubset[ LAST ] = MISSING; /* Not found*/
  } else {
    const int yyyymmddhh1 = timestepSubset[ FIRST ];
    const int yyyymmddhh2 = timestepSubset[ LAST  ];
    int yyyymmddhh = 0;
    int index = 0;

    /* Find first timestep: */

    do {
      yyyymmddhh = timestamps[ index ];

      if ( yyyymmddhh == yyyymmddhh1 ) {
        timestepSubset[ FIRST ] = index;
        index = timesteps; /* Stop looping. */
      } else if ( yyyymmddhh > yyyymmddhh1 ) {
        timestepSubset[ FIRST ] = index > 0 ? index - 1 : index;
        index = timesteps; /* Stop looping. */
      }

      ++index;
    }  while  ( index < timesteps );

    assert( IN_RANGE( timestepSubset[ FIRST ], 0, timesteps - 1 ) );

    /* Find last timestep: */

    index = timestepSubset[ LAST ] = timestepSubset[ FIRST ];

    do {
      yyyymmddhh = timestamps[ index ];

      if ( yyyymmddhh <= yyyymmddhh2 ) {
        timestepSubset[ LAST ] = index;
      } else {
        index = timesteps; /* Stop looping. */
      }

      ++index;
    }  while ( index < timesteps );

    assert( IN_RANGE( timestepSubset[ LAST ], timestepSubset[ FIRST ],
                      timesteps - 1 ) );
  }
}


/******************************************************************************
PURPOSE: computeSubset - Compute aggregated subset of grid data.
INPUTS:  int rows             Number of grid rows.
         int columns          Number of grid columns.
         Bounds domain        Grid domain.
         Bounds subset        Target subset domain.
         int dimension        Maximum target output dimension.
OUTPUTS: Bounds domain        Adjusted grid domain.
         Bounds subset        Possibly adjusted subset domain.
         int* size            Aggregation size (number of rows/columns to agg).
         int* subsetRows      Number of output rows.
         int* subsetColumns   Number of output columns.
         Range range          0-based indices of subset.
RETURNS: int 1 if a non-empty subset, else 0.
******************************************************************************/

static int computeSubset( int rows, int columns, Bounds domain,
                          Bounds subset, int dimension, int* size,
                          int* subsetRows, int* subsetColumns, Range range ) {

  int result = 0;

  assert( rows > 0 ); assert( columns > 0 );
  assert( isValidBounds( (const double (*)[2]) domain ) );
  assert( isValidBounds( (const double (*)[2]) subset ) );
  assert( dimension > 0 );
  assert( size ); assert( subsetRows ); assert( subsetColumns );
  assert( range );

  /* Initialize outputs: */

  *size = *subsetRows = *subsetColumns = 0;
  memset( range, 0, sizeof (Range) );

  /* Subset rows: */

  subsetIndices( subset[ LATITUDE ], domain[ LATITUDE ], rows, range[ ROW ] );

  if ( range[ ROW ][ 0 ] != -1 ) { /* Intersected: */

    /* Convert to north-to-south row-order: */

    const int swapTemp = range[ ROW ][ FIRST ];
    range[ ROW ][ FIRST ] = range[ ROW ][ LAST ];
    range[ ROW ][ LAST  ] = swapTemp;
    range[ ROW ][ FIRST ] = rows - 1 - range[ ROW ][ FIRST ];
    range[ ROW ][ LAST  ] = rows - 1 - range[ ROW ][ LAST ];

    /* Subset columns: */

    subsetIndices( subset[ LONGITUDE ], domain[ LONGITUDE ],
                   columns, range[ COLUMN ] );

    DEBUG( fprintf( stderr, "initial range: columns [%d %d] rows [%d %d]\n",
                    range[ COLUMN ][ FIRST ], range[ COLUMN ][ LAST ],
                    range[ ROW ][ FIRST ], range[ ROW ][ LAST ] ); )

    if ( range[ COLUMN ][ 0 ] != -1 ) { /* Intersected: */
      result = 1;
      *subsetColumns = range[ COLUMN ][ LAST ] - range[ COLUMN ][ FIRST ] + 1;
      *subsetRows    = range[ ROW    ][ LAST ] - range[ ROW    ][ FIRST ] + 1;
      *size = 1;
      memcpy( subset, domain, sizeof (Bounds) );

      /* Ensure subset size does not exceed maximum target dimension: */

      computeStride( rows, columns, dimension,
                     subsetRows, subsetColumns, size, range, subset );
    }
  }

  if ( ! result ) { /* Clear outputs: */
    *size = *subsetRows = *subsetColumns = 0;
    memset( range, 0, sizeof (Range) );
  }

  DEBUG( fprintf( stderr, "computeSubset: result = %d, "
                  "*size = %d, *subsetRows = %d, *subsetColumns = %d, "
                  "range = columns [%d %d] rows [%d %d]\n",
                  result, *size, *subsetRows, *subsetColumns,
                  range[ COLUMN ][ FIRST ], range[ COLUMN ][ LAST ],
                  range[ ROW ][ FIRST ], range[ ROW ][ LAST ] ); )

  assert( IS_BOOL( result ) );
  assert( isValidBounds( (const double (*)[2]) domain ) );
  assert( isValidBounds( (const double (*)[2]) subset ) );
  assert( IMPLIES_ELSE( result,
                        AND9( *size > 0,
                              IN_RANGE( range[ COLUMN ][FIRST], 0, columns-1),
                              IN_RANGE( range[ COLUMN ][LAST],
                                        range[ COLUMN ][ FIRST ], columns - 1),
                              IN_RANGE( range[ ROW ][ FIRST ], 0, rows - 1 ),
                              IN_RANGE( range[ ROW ][ LAST ],
                                        range[ ROW ][ FIRST ], rows - 1 ),
                              IN_RANGE( *subsetRows, 1, rows ),
                              *subsetRows ==
                                ( range[ ROW ][ LAST ] + *size - 1 -
                                  range[ ROW ][ FIRST ] + 1 ) / *size,
                              IN_RANGE( *subsetColumns, 1, columns ),
                              *subsetColumns ==
                                ( range[ COLUMN ][ LAST ] + *size - 1 -
                                  range[ COLUMN ][ FIRST ] + 1 ) / *size ),
                        IS_ZERO7( *size, *subsetRows, *subsetColumns,
                                  range[ 0 ][ 0 ], range[ 0 ][ 1 ],
                                  range[ 1 ][ 0 ], range[ 1 ][ 1 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: subsetIndices - Subset image range based on clip and compute indices.
INPUTS:  const double clip[ 2 ]   clip [ MINIMUM MAXIMUM ] lon-lat bounds.
         double range[ 2 ]        range[ MINIMUM MAXIMUM ] of unclipped image.
         int count                Pixels in image along subset dimension.
OUTPUTS: double range[ 2 ]        range[ MINIMUM MAXIMUM ] of clipped image.
         int indices[ 2 ]         indices[ MINIMUM MAXIMUM ] of clipped image.
******************************************************************************/

static void subsetIndices( const double clip[ 2 ],
                           double range[ 2 ], int count, int indices[ 2 ] ) {

  const double clipMinimum = clip ? clip[ MINIMUM ] : 0.0;
  const double clipMaximum = clip ? clip[ MAXIMUM ] : 0.0;
  const double rangeMinimum = range ? range[ MINIMUM ] : 0.0;
  const double rangeMaximum = range ? range[ MAXIMUM ] : 0.0;
  const double clipRange = clipMaximum - clipMinimum;
  const double rangeRange = rangeMaximum - rangeMinimum;
  const double TOO_SMALL = 1e-6;
  const double margin = 1.0;

  assert( clip );
  assert( range );
  assert( count > 0 );
  assert( indices );

  indices[ MINIMUM ] = indices[ MAXIMUM ] = -1;

  DEBUG( fprintf( stderr,
                  "subsetIndices input:  "
                  "range[%lf %lf] clip(%lf %lf) count = %d\n",
                  range[ MINIMUM ], range[ MAXIMUM ],
                  clip[ MINIMUM ], clip[ MAXIMUM ], count ); )

  if ( ! OR4( clipRange < TOO_SMALL, rangeRange < TOO_SMALL,
              rangeMaximum < clipMinimum, rangeMinimum > clipMaximum ) ) {
    const double scale = 1.0 / rangeRange;
    const double rangeIncrement = rangeRange / count;

    DEBUG( fprintf( stderr, "scale = %lf, rangeIncrement = %lf\n",
                    scale, rangeIncrement ); )

    if ( rangeMinimum > clipMinimum ) {
      indices[ MINIMUM ] = 0;
      DEBUG( fprintf( stderr, "set indices[ MINIMUM ] = 0\n" ); )
    } else {
      const double interpolation = ( clipMinimum - rangeMinimum ) * scale;
      DEBUG( fprintf( stderr, "interpolation = %lf\n", interpolation ); )
      indices[ MINIMUM ] = (int) ( interpolation * count - 0.5 );
      DEBUG( fprintf( stderr, "indices[ MINIMUM ] = %d\n", indices[ MINIMUM]);)
      indices[ MINIMUM ] = CLAMPED_TO_RANGE( indices[ MINIMUM ], 0, count - 1);
      DEBUG( fprintf( stderr, "indices[ MINIMUM ] = %d\n", indices[ MINIMUM]);)
      range[ MINIMUM ] += indices[ MINIMUM ] * rangeIncrement;
      DEBUG( fprintf( stderr,
                      "indices[ MINIMUM ] = %d, range[ MINIMUM ] = %lf\n",
                      indices[ MINIMUM ], range[ MINIMUM ] ); )
    }

    if ( rangeMaximum < clipMaximum ) {
      indices[ MAXIMUM ] = count - 1;
      DEBUG( fprintf( stderr, "set indices[ MAXIMUM ] = count - 1\n" ); )
    } else {
      const double interpolation = ( clipMaximum - rangeMinimum ) * scale;
      DEBUG( fprintf( stderr, "interpolation = %lf\n", interpolation ); )
      indices[ MAXIMUM ] = (int) ( interpolation * count + 0.5 );
      DEBUG( fprintf( stderr, "indices[ MAXIMUM ] = %d\n", indices[ MAXIMUM]);)
      indices[ MAXIMUM ] =
        CLAMPED_TO_RANGE( indices[ MAXIMUM ], indices[ MINIMUM ], count - 1 );
      DEBUG( fprintf( stderr,
                      "indices[ MAXIMUM ] = %d\n", indices[ MAXIMUM ] ); )
    }

    range[ MAXIMUM ] = range[ MINIMUM ] +
      ( indices[ MAXIMUM ] - indices[ MINIMUM ] ) * rangeIncrement;

    DEBUG( fprintf( stderr, "range[ MAXIMUM ] = %lf\n", range[ MAXIMUM ] ); )

    if ( ! AND4( IN_RANGE( range[ MINIMUM ],
                           clip[ MINIMUM ] - margin, clip[ MAXIMUM ] + margin ),
                 IN_RANGE( range[ MAXIMUM ],
                           range[ MINIMUM ], clip[ MAXIMUM ] + margin ),
                 IN_RANGE( indices[ MINIMUM ], 0, count - 1 ),
                 IN_RANGE( indices[ MAXIMUM ],
                           indices[ MINIMUM ], count - 1 ) ) ) {

      DEBUG( fprintf( stderr, "set indices[] = -1\n" ); )
      DEBUG( fprintf( stderr, "because:\n"
                      "IN_RANGE( range[ MINIMUM ], "
                      "clip[ MINIMUM ] - 0.1, clip[ MAXIMUM ] + margin ) = %d,\n"
                      "IN_RANGE( range[ MAXIMUM ], "
                      "range[ MINIMUM ], clip[ MAXIMUM ] + margin ) = %d,\n"
                      "IN_RANGE( indices[ MINIMUM ], 0, count - 1 ) = %d,\n"
                      "IN_RANGE( indices[ MAXIMUM ],"
                      "indices[ MINIMUM ], count - 1 ) = %d\n",
                      IN_RANGE( range[ MINIMUM ],
                                clip[ MINIMUM ] - margin, clip[ MAXIMUM ] + margin ),
                      IN_RANGE( range[ MAXIMUM ],
                                range[ MINIMUM ], clip[ MAXIMUM ] + margin ),
                      IN_RANGE( indices[ MINIMUM ], 0, count - 1 ),
                      IN_RANGE( indices[ MAXIMUM ],
                                indices[ MINIMUM ], count - 1 ) ) ; )

      indices[ MINIMUM ] = indices[ MAXIMUM ] = -1;
    }
  }

  DEBUG( fprintf( stderr,
                  "subsetIndices output: "
                  "range[%lf %lf] indices[%d %d]\n",
                  range[ MINIMUM ], range[ MAXIMUM ],
                  indices[ MINIMUM ], indices[ MAXIMUM ] ); )

  assert( IMPLIES( ! AND2( indices[ MINIMUM ] == -1, indices[ MAXIMUM ] == -1),
                     AND4( IN_RANGE( range[ MINIMUM ],
                                     clip[MINIMUM] - margin, clip[MAXIMUM] + margin),
                           IN_RANGE( range[ MAXIMUM ],
                                     range[ MINIMUM ], clip[ MAXIMUM ] + margin ),
                           IN_RANGE( indices[ MINIMUM ], 0, count - 1 ),
                           IN_RANGE( indices[ MAXIMUM ],
                                     indices[ MINIMUM ], count - 1 ) ) ) );
}



/******************************************************************************
PURPOSE: computeStride - Compute stride so subsetDimension does not exceed
         maximum target dimension.
INPUTS:  int rows             Number of grid rows.
         int columns          Number of grid columns.
         Bounds subset        Target subset domain.
         int dimension        Maximum target output dimension.
OUTPUTS: Bounds domain        Adjusted grid domain.
         Bounds subset        Possibly adjusted subset domain.
         int* size            Aggregation size (# of rows/columns to agg).
         int* subsetRows      Number of output rows.
         int* subsetColumns   Number of output columns.
         Range range          0-based indices of subset.
RETURNS: int 1 if a non-empty subset, else 0.
*****************************************************************************/

static void computeStride( int rows, int columns, int dimension,
                           int* subsetRows, int* subsetColumns, int* size,
                           Range range, Bounds subset ) {

  assert( rows > 0 ); assert( columns > 0 ); assert( dimension > 0 );
  assert( subsetRows );
  assert( IN_RANGE( *subsetRows, 1, rows ) );
  assert( subsetColumns );
  assert( IN_RANGE( *subsetColumns, 1, columns ) );
  assert( size );
  assert( range );
  assert( IN_RANGE( range[ COLUMN ][FIRST], 0, columns - 1 ) );
  assert( IN_RANGE( range[ COLUMN ][LAST],
                    range[ COLUMN ][ FIRST ], columns - 1 ) );
  assert( IN_RANGE( range[ ROW ][ FIRST ], 0, rows - 1 ) );
  assert( IN_RANGE( range[ ROW ][ LAST ], range[ ROW ][ FIRST ], rows - 1 ) );
  assert( *subsetRows == ( range[ ROW ][ LAST ] - range[ ROW ][ FIRST ] + 1 ));
  assert( *subsetColumns ==
          ( range[ COLUMN ][ LAST ] - range[ COLUMN ][ FIRST ] + 1 ) );
  assert( isValidBounds( (const double (*)[2]) subset ) );

  DEBUG( fprintf( stderr, "computeStride: before dimension = %d "
                 "*size = %d, *subsetRows = %d, *subsetColumns = %d, "
                 "range = columns [%d %d] rows [%d %d] "
                 "subset = [%g %g] [ %g %g]\n",
                 dimension, *size, *subsetRows, *subsetColumns,
                 range[ COLUMN ][ FIRST ], range[ COLUMN ][ LAST ],
                 range[ ROW ][ FIRST ], range[ ROW ][ LAST ],
                 subset[ LONGITUDE ][ MINIMUM ], subset[ LONGITUDE ][ MAXIMUM],
                 subset[ LATITUDE ][ MINIMUM ], subset[ LATITUDE ][ MAXIMUM]);)

  {
    const double cellWidth =
      ( subset[ LONGITUDE ][ MAXIMUM ] - subset[ LONGITUDE ][ MINIMUM ] ) /
      *subsetColumns;

    const double cellHeight =
      ( subset[ LATITUDE ][ MAXIMUM ] - subset[ LATITUDE ][ MINIMUM ] ) /
      *subsetRows;

    const int maximumSubsetDimension = MAX( *subsetColumns, *subsetRows );
    int adjustedSize = 0;
    *size = 1;

    if ( maximumSubsetDimension > dimension ) {
      *size = maximumSubsetDimension / dimension +
             (maximumSubsetDimension % dimension != 0);

      if ( *size > 1 ) {
        const int minimumSubsetDimension = MIN( *subsetColumns, *subsetRows );
        *size = CLAMPED_TO_RANGE( *size, 1, minimumSubsetDimension );
        *subsetColumns /= *size;
        *subsetRows    /= *size;
        adjustedSize = 1;
      }
    }

    DEBUG( fprintf( stderr, "  *size = %d\n", *size ); )

    if ( adjustedSize ) {
      int last = 0;
      int count = 0;

      /* Adjust column range: */

      for ( last =  range[ COLUMN ][ FIRST ], count = 1;
            last < range[ COLUMN ][ LAST  ] + *size;
            last += *size, ++count ) {
      }

      DEBUG(fprintf(stderr, "  column: count = %d, last = %d\n", count,last);)

      if ( last == range[ COLUMN ][ LAST ] + 1 ) { // Aligned to subset:
        last -= *size;
        range[ COLUMN ][ LAST ] = last;
        DEBUG( fprintf( stderr, "  aligned last column = %d\n", last ); )
      } else { // Non-aligned so adjust subset:

        while ( last + *size > columns ) {
          last -= *size;
        }

        if ( last < range[ COLUMN ][ FIRST ] ) {
          last = range[ COLUMN ][ FIRST ];
        }

        DEBUG( fprintf( stderr, "  non-aligned last column = %d\n", last ); )

        assert( last >= range[ COLUMN ][ FIRST ] ); assert( last < columns );
        range[ COLUMN ][ LAST ] = last;
        *subsetColumns =
          ( range[ COLUMN ][ LAST ] + *size - 1 - range[ COLUMN ][ FIRST ] + 1)
          / *size;
        subset[ LONGITUDE ][ MAXIMUM ] =
          subset[ LONGITUDE ][ MINIMUM ] + cellWidth * *subsetColumns * *size;
      }

      /* Adjust row range: */

      for ( last =  range[ ROW ][ FIRST ], count = 1;
            last < range[ ROW ][ LAST  ] + *size;
            last += *size, ++count ) {
      }

      DEBUG( fprintf( stderr, "  row: count = %d, last = %d\n", count, last );)

      if ( last == range[ ROW ][ LAST ] + 1 ) { // Aligned to subset:
        last -= *size;
        range[ ROW ][ LAST ] = last;
        DEBUG( fprintf( stderr, "  aligned last row = %d\n", last ); )
      } else { // Non-aligned so adjust subset:

        while ( last + *size > rows ) {
          last -= *size;
        }

        if ( last < range[ ROW ][ FIRST ] ) {
          last = range[ ROW ][ FIRST ];
        }

        DEBUG( fprintf( stderr, "  non-aligned last row = %d\n", last ); )

        assert( last >= range[ ROW ][ FIRST ] ); assert( last < rows );
        range[ ROW ][ LAST ] = last;
        *subsetRows =
          ( range[ ROW ][ LAST ] + *size - 1 - range[ ROW ][ FIRST ] + 1 )
          / *size;
        subset[ LATITUDE ][ MINIMUM ] =
          subset[ LATITUDE ][ MAXIMUM ] - cellHeight * *subsetRows * *size;
      }
    }
  }

  DEBUG( fprintf( stderr, "computeStride: after "
                 "*size = %d, *subsetRows = %d, *subsetColumns = %d, "
                 "range = columns [%d %d] rows [%d %d] "
                 "subset = [%g %g] [ %g %g]\n",
                 *size, *subsetRows, *subsetColumns,
                 range[ COLUMN ][ FIRST ], range[ COLUMN ][ LAST ],
                 range[ ROW ][ FIRST ], range[ ROW ][ LAST ],
                 subset[ LONGITUDE ][ MINIMUM ], subset[ LONGITUDE ][ MAXIMUM],
                 subset[ LATITUDE ][ MINIMUM ], subset[ LATITUDE ][ MAXIMUM]);)

  assert( IN_RANGE( *subsetRows, 1, rows ) );
  assert( IN_RANGE( *subsetColumns, 1, columns ) );
  assert( *size >= 1 ); assert( *size <= rows ); assert( *size <= columns );
  assert( IN_RANGE( range[ COLUMN ][FIRST], 0, columns - 1 ) );
  assert( IN_RANGE( range[ COLUMN ][LAST],
                   range[ COLUMN ][ FIRST ], columns - 1 ) );
  assert( IN_RANGE( range[ ROW ][ FIRST ], 0, rows - 1 ) );
  assert( IN_RANGE( range[ ROW ][ LAST ], range[ ROW ][ FIRST ], rows - 1 ) );
  assert( *subsetRows ==
          ( range[ ROW ][ LAST ] + *size - 1 - range[ ROW ][ FIRST ] + 1 )
          / *size );
  assert( *subsetColumns ==
          ( range[ COLUMN ][ LAST ] + *size - 1 - range[ COLUMN ][ FIRST ] + 1)
          / *size );
  assert( isValidBounds( (const double (*)[2]) subset ) );
}



/******************************************************************************
PURPOSE: aggregate - Aggregates a set of row data.
INPUTS:  int size           Number of input rows / columns to aggregate.
         int columns        Number of input columns.
         int firstColumn    0-based index of first column to aggregate.
         int outputColumns  Number of aggregated output columns.
         int method         MEAN or MODE.
         int type           BYTE_TYPE, FLOAT_TYPE, etc.
         const void* input  Row data input[ size * columns ] to aggregate.
OUTPUTS: void* output       Aggregated output[ size ] data row.
NOTES:   Not thread-safe due to use of static variables.
******************************************************************************/

static void aggregate( int size, int columns, int firstColumn,
                       int outputColumns, int method, int type,
                       const void* input, void* output ) {

  const signed char* const cinput = input;
  signed char* const coutput = output;
  const unsigned short* const sinput = input;
  unsigned short* const soutput = output;
  const float* finput = input;
  float* const foutput = output;
  int startColumn = firstColumn;
  int outputColumn = 0;
  int counts[ 256 ];
  static int scounts[ 65536 ]; /* static because too large for the stack. */
  memset( counts, 0, sizeof counts );

  assert( size > 0 ); assert( columns > 0 );
  assert( IN_RANGE( firstColumn, 0, columns - 1 ) );
  assert( IN_RANGE( outputColumns, 1, columns ) );
  assert( IN3( method, MEAN, MODE ) );
  assert( IS_VALID_TYPE( type ) );
  assert( IMPLIES( method == MODE, type != FLOAT_TYPE ) );
  assert( input ); assert( output );

  DEBUG( fprintf( stderr, "aggregate: size = %d, columns = %d, "
                  "firstColumn = %d, outputColumns = %d, method = %d, "
                  "type = %d\n",
                  size, columns, firstColumn, outputColumns, method, type);)

  for ( outputColumn = 0; outputColumn < outputColumns;
        ++outputColumn, startColumn += size ) {
    int row = 0;
    int missingCount = 0;
    int count = 0;
    double mean = 0.0;

    if ( method == MODE ) {

      if ( type == UINT16_TYPE ) {
        memset( scounts, 0, sizeof scounts );
      } else {
        memset( counts, 0, sizeof counts );
      }
    }

    for ( row = 0; row < size; ++row ) {
      const int endColumn = startColumn + size;
      const int rowOffset = row * columns;
      int column = 0;

      for ( column = startColumn; column < endColumn; ++column ) {
        const int index = rowOffset + column;
        const signed char cvalue = type == BYTE_TYPE ? cinput[ index ] : 0;
        const unsigned short svalue = type == UINT16_TYPE ? sinput[ index] : 0;
        const float fvalue = type == FLOAT_TYPE ? finput[ index ] : 0.0;

        DEBUG2( fprintf( stderr, "input[%d] = (%d, %d, %f)\n",
                         index, cvalue, svalue, fvalue ); )

        if ( method == MODE ) {

          if ( type == UINT16_TYPE ) {
            const int valueIndex = svalue;
            assert( IN_RANGE( valueIndex, 0, 65535 ) );
            scounts[ valueIndex ] += 1;
          } else {
            const unsigned char valueIndex = (unsigned char) cvalue;
            assert( type == BYTE_TYPE );
            counts[ valueIndex ] += 1;
          }
        } else {
          assert( method == MEAN );

          if ( cvalue == MISSING || fvalue == F_MISSING ||
               ( type == UINT16_TYPE && svalue == 0 ) ) { /* 0 = missing. */
            ++missingCount;
          } else {
            const float value =
              type == BYTE_TYPE ? cvalue
              : type == UINT16_TYPE ? svalue
              : fvalue;
            const int count1 = count + 1;
            mean = ( count * mean + value ) / count1;
            count = count1;
          }
        }
      } /* for column. */
    } /* for row. */

    if ( method == MODE ) {

      if ( type == UINT16_TYPE ) {
        const int mode =
          indexOfMaximum( scounts, sizeof scounts / sizeof *scounts );
        assert( IN_RANGE( mode, 0, 65535 ) );
        soutput[ outputColumn ] = mode;
      } else {
        const int mode =
          indexOfMaximum( counts, sizeof counts / sizeof *counts );
        assert( type == BYTE_TYPE );
        coutput[ outputColumn ] = mode;
      }
    } else {
      DEBUG2( fprintf( stderr, "mean = %lf, count = %d, missingCount = %d\n",
                       mean, count, missingCount ); )
      assert( method == MEAN );

      if ( type == FLOAT_TYPE ) {
        foutput[ outputColumn ] = missingCount > count ? F_MISSING : mean;
      } else if ( type == UINT16_TYPE ) {
        const int m = missingCount > count ? 0 : (int) ( mean + 0.5 );
        assert( IN_RANGE( m, 0, 65535 ) );
        soutput[ outputColumn ] = m;
      } else {
        const int m = missingCount > count ? MISSING : (int) ( mean + 0.5 );
        assert( type == BYTE_TYPE );
        coutput[ outputColumn ] = m;
      }
    }

    DEBUG2( fprintf( stderr, "output[%d] = (%d, %d, %f)\n",
                     outputColumn,
                     coutput[ outputColumn ],
                     soutput[ outputColumn ],
                     foutput[ outputColumn ] ); )
  }
}



/******************************************************************************
PURPOSE: indexOfMaximum - 0-based index of maximum value in array.
INPUTS:  const int array[ count ]  Array to scan.
         int count                 Number of in array.
RETURNS: int index of largest item in array.
******************************************************************************/

static int indexOfMaximum( const int array[], int count ) {
  int result = 0;
  int index = 0;
  int maximum = 0;
  assert( array ); assert( count > 0 );
  maximum = array[ 0 ];

  for ( index = 1; index < count; ++index ) {
    const int value = array[ index ];

    if ( value > maximum ) {
      maximum = value;
      result = index;
    }
  }

  assert( IN_RANGE( result, 0, count - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: allocate - Allocate memory by calling malloc() to allocate
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
PURPOSE: isValidBounds - Is bounds a valid lon-lat rectangle?
INPUTS:  const Bounds bounds  Bounds to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidBounds( const Bounds bounds ) {
  const int result =
    AND5( bounds,
          IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ),
          IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
                    bounds[ LONGITUDE ][ MINIMUM ], 180.0 ),
          IN_RANGE( bounds[ LATITUDE ][ MINIMUM ], -90.0, 90.0 ),
          IN_RANGE( bounds[ LATITUDE ][ MAXIMUM ],
                    bounds[ LATITUDE ][ MINIMUM ], 90.0 ) );
  assert( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: rotate4ByteArrayIfLittleEndian() - Rotate 4-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array      Array of 4-byte values to rotate.
         long long count  Number of items in array.
OUTPUTS: void* array      Array of rotated values.
******************************************************************************/

static void rotate4ByteArrayIfLittleEndian( void* array, long long count ) {

#if IS_LITTLE_ENDIAN

  int* const array4 = array;
  long long index = 0;
  assert_static( sizeof (int) == 4 );
  assert( array ); assert( count > 0 );

#pragma omp parallel for

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



/******************************************************************************
PURPOSE: rotate2ByteArrayIfLittleEndian() - Rotate 2-bytes of each array item
         if on a little-endian platform.
INPUTS:  void* array      Array of 2-byte values to rotate.
         long long count  Number of items in array.
OUTPUTS: void* array      Array of rotated values.
******************************************************************************/

static void rotate2ByteArrayIfLittleEndian( void* array, long long count ) {

#if IS_LITTLE_ENDIAN

  unsigned short* const array2 = array;
  long long index = 0;
  assert_static( sizeof (unsigned short) == 2 );
  assert( array ); assert( count > 0 );

  for ( index = 0; index < count; ++index ) {
    const unsigned short value = array2[ index ];
    const unsigned short newValue =
      ( value & 0x00ff ) <<  8 |
      ( value & 0xff00 ) >>  8;
    array2[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}



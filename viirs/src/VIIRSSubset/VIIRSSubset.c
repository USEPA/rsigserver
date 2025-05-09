/******************************************************************************
PURPOSE: VIIRSSubset.c - Extract a lon-lat subset of data from a list of
         VIIRS NetCDF4 files and write it to stdout as XDR binary format.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).
         Compile:
         gcc -Wall -g -o VIIRS2Subset VIIRS2Subset.c Utilities.c ReadData.c \
                   -I../../../include/NetCDF4 \
                   -L../../../lib/$platform \
                   -lnetcdf4 -lhdf5_hl -lhdf5 -lcurl -lz -ldl -lm -lc

         Usage:
         VIIRS2Subset -files <listfile> \
                      -tmpdir <temp_directory> \
                      -desc "description text" \
                      -timestamp <yyyymmddhh> -hours <count> \
                      -variable <name> \
                       -domain <minimum_longitude> <minimum_latitude> \
                              <maximum_longitude> <maximum_latitude> \
                      [-corners]
                      [-minimumQuality [medium|high]] (default is high)

          Example:
          ../../../bin/$platform/VIIRS2Subset \
          -files testdata/file_list \
          -tmpdir testdata \
          -desc "https://ncc.nesdis.noaa.gov/VIIRS/,VIIRS2Subset" \
          -timestamp 2017080100 -hours 24 \
          -variable AOD550 \
          -domain -75 35 -70 36 \
          -corners \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 14-line ASCII header followed by binary 64-bit big-endian arrays:

Swath 2.0
https://ncc.nesdis.noaa.gov/VIIRS/,VIIRS2Subset
2013-06-15T20:00:00-0000
# Dimensions: variables timesteps scans:
11 24 2
# Variable names:
Longitude Latitude AOD550 \
 LongitudeSW LongitudeSE LongitudeNW LongitudeNE \
 LatitudeSW LatitudeSE LatitudeNW LatitudeNE
# Variable units:
deg deg - deg deg deg deg deg deg deg deg
# Domain: <min_lon> <min_lat> <max_lon> <max_lat>
-90 35 -50 40
# MSB 64-bit integers (yyyydddhhmm) timestamps[scans] and
# MSB 64-bit integers points[scans] and
# IEEE-754 64-bit reals data_1[variables][points_1] ... \
 data_S[variables][points_S]:
<big-endian binary format arrays>
20172131804
20172131806
5
122
-7.1847106933593750e+01
-7.1855308532714844e+01
 ...
3.5999182701110840e+01
3.5997957229614258e+01


HISTORY: 2017-10-14 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, printf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strchr(). */
#include <stdlib.h>    /* For malloc(), free(). */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For LONGITUDE, Bounds, expand32BitReals(). */
#include "ReadData.h"  /* For openFile(), swathInDomain(), readFileData(). */

/*================================= MACROS =================================*/

/* Name of temporary file created in -tmpdir will have PID appended: */

#define TEMP_FILE_NAME "junk_VIIRS2Subset"

/*================================== TYPES ==================================*/

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;    /* File containing list of VIIRS files to read. */
  const char* tmpdir;      /* Name of directory to write temp files. */
  const char* description; /* User-supplied description. */
  const char* variable;    /* Name of variable to read. */
  const char* minimumQuality; /* "medium" or "high". Default is "high". */
  Bounds      domain; /* Subset domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  int         yyyymmddhh;  /* First timestamp of subset. */
  int         hours;       /* Number of hours in subset. */
  int         corners;     /* Compute interpolated lon-lat corner points?*/
} Arguments;

/* Data type: */

typedef struct {
  Arguments  arguments;    /* User-supplied (command-line) arguments. */
  FileName   tempFileName; /* Name of temp file of output subset data. */
  FILE*      tempFile;     /* Temp file of output subset data. */
  long long* yyyydddhhmm;  /* Timestamp per output subset scan. */
  long long* points;       /* Number of points per output subset scan. */
  int        scans;        /* Number of output subset scans. */
  int        ok;           /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocate( Data* data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static void readData( Data* const data );

static char* readListFileAndAllocateTimestampsAndPoints( Data* const data );

static void allocateTimestampsAndPoints( const char* const listFileContent,
                                         Data* const data );

static int fileInDomain( Data* const data,
                         const char* const fileName,
                         int* const file,
                         long long* const yyyymmddhhmm,
                         size_t* const rows,
                         size_t* const columns,
                         size_t* const size,
                         int* const changedDimensions );

static void readCoordinatesAndValues( Data* const data,
                                      const int file,
                                      const size_t rows,
                                      const size_t columns,
                                      double* const longitudes,
                                      double* const latitudes,
                                      double* const values );

static long long fileTimestamp( const char* const fileName );

static void writeDataSubset( Data* const data,
                             const long long yyyymmddhhmm,
                             const size_t rows,
                             const size_t columns,
                             const double longitudes[],
                             const double latitudes[],
                             const double values[],
                             unsigned char mask[],
                             double longitudesSW[],
                             double longitudesSE[],
                             double longitudesNW[],
                             double longitudesNE[],
                             double latitudesSW[],
                             double latitudesSE[],
                             double latitudesNW[],
                             double latitudesNE[] );

static void writeSubset( Data* const data,
                         const long long  yyyymmddhhmm,
                         const size_t subsetPoints,
                         const size_t points,
                         const unsigned char mask[],
                         const double longitudes[],
                         const double latitudes[],
                         const double values[],
                         const double longitudesSW[],
                         const double longitudesSE[],
                         const double longitudesNW[],
                         const double longitudesNE[],
                         const double latitudesSW[],
                         const double latitudesSE[],
                         const double latitudesNW[],
                         const double latitudesNE[] );

static void streamData( Data* const data );

static void streamHeader( const Data* const data );


/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of VIIRS files and write
         it to stdout in XDR format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int ok = 0;
  Data data;
  memset( &data, 0, sizeof data );
  data.ok = parseArguments( argc, argv, &data.arguments );

  if ( ! data.ok ) {
    printUsage( argv[ 0 ] );
  } else {
    readData( &data ); /* Read subset of VIIRS files and write temp files. */

    if ( data.ok && data.scans ) {
      streamData( &data ); /* Write header and temp file to stdout & rm temp.*/
      ok = data.ok;
    }
  }

  deallocate( &data );
  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocate - Deallocate data.
INPUTS:  Data* data  Data to deallocate.
******************************************************************************/

static void deallocate( Data* data ) {
  assert( data );

  if ( data->points ) {
    free( data->points );
    data->points = 0;
  }

  if ( data->yyyydddhhmm ) {
    free( data->yyyydddhhmm );
    data->yyyydddhhmm = 0;
  }

  if ( data->tempFile ) {
    fclose( data->tempFile );
    data->tempFile = 0;
  }

  if ( data->tempFileName[ 0 ] ) {
    unlink( data->tempFileName );
    memset( data->tempFileName, 0, sizeof data->tempFileName );
  }
}



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* name  Name of program.
******************************************************************************/

static void printUsage( const char* name ) {
  assert( name ); assert( *name );
  fprintf( stderr,
           "\a\n\n%s - Extract a lon-lat subset of data from a list of\n"
          "VIIRS NetCDF4 files and write it to stdout as XDR binary format.\n",
           name );
  fprintf( stderr, "Data is subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "  -files <listfile> \\\n" );
  fprintf( stderr, "  -tmpdir <temp_directory> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, "  -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, "  [-corners]\n" );
  fprintf( stderr, "  [-minimumQuality [medium | high] (defualt is high)\n\n" );
  fprintf( stderr, "Note:\ntimestamp is in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory were a transient file is " );
  fprintf ( stderr, "written.\nIt should have enough disk space (1TB).\n" );
  fprintf( stderr, "-corners option will output 8 additional variables:\n" );
  fprintf( stderr, "  Longitude_SW Longitude_SE Longitude_NW Longitude_NE\n");
  fprintf( stderr, "  Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "that are the linearly interpolated " );
  fprintf( stderr, "(and edge extrapolated)\n" );
  fprintf( stderr, "corner points for each center-pixel point.\n\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-files JPSSRR_AOD_files \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc "
                   "\"http://www.star.nesdis.noaa.gov/smcd/emb/viirs_aerosol/"
                   ",VIIRS2Subset\" \\\n");
  fprintf( stderr, "-timestamp 2017081500 -hours 24 \\\n" );
  fprintf( stderr, "-variable AOD550 \\\n" );
  fprintf( stderr, "-domain -126 25 -65 50 -corners > subset.xdr\n\n" );
  fprintf( stderr, "AOD over US on August 15, 2013.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Swath 2.0\n" );
  fprintf( stderr, "http://www.star.nesdis.noaa.gov/smcd/emb/viirs_aerosol\n");
  fprintf( stderr, "2017-08-15T00:00:00-0000\n" );
  fprintf( stderr, "# Dimensions: variables timesteps scans:\n" );
  fprintf( stderr, "11 24 2\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Longitude Latitude AOD550 "
                   "Longitude_SW Longitude_SE Longitude_NW Longitude_NE"
                   "Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "deg deg - deg deg deg deg deg deg deg deg\n" );
  fprintf( stderr, "# Domain: <min_lon> <min_lat> <max_lon> <max_lat>\n" );
  fprintf( stderr, "-126 25 -65 50\n" );
  fprintf( stderr, "# MSB 64-bit integers (yyyydddhhmm)" );
  fprintf( stderr, " timestamps[scans] and\n" );
  fprintf( stderr, "# MSB 64-bit integers points[scans] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals" );
  fprintf( stderr, " data_1[variables][points_1] ..." );
  fprintf( stderr, " data_S[variables][points_S]:\n" );
  fprintf( stderr, "<big-endian binary format arrays>\n" );
  fprintf( stderr, "20172131804\n" );
  fprintf( stderr, "20172131806\n" );
  fprintf( stderr, "5\n" );
  fprintf( stderr, "122\n" );
  fprintf( stderr, "-7.1847106933593750e+01\n" );
  fprintf( stderr, "-7.1855308532714844e+01\n" );
  fprintf( stderr, " ...\n" );
  fprintf( stderr, "3.5999182701110840e+01\n" );
  fprintf( stderr, "3.5997957229614258e+01\n" );
  fprintf( stderr, "\n\n\n");
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int parseArguments(int argc, char* argv[], Arguments* const arguments) {
  int result = 0;
  int arg = 0;
  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] ); assert( arguments );
  memset( arguments, 0, sizeof (Arguments) );
  arguments->domain[ LONGITUDE ][ MINIMUM ] = -180.0;
  arguments->domain[ LONGITUDE ][ MAXIMUM ] =  180.0;
  arguments->domain[ LATITUDE  ][ MINIMUM ] =  -90.0;
  arguments->domain[ LATITUDE  ][ MAXIMUM ] =   90.0;
  arguments->minimumQuality = "high";

  result = argc == 18 || argc == 19 || argc == 20 || argc == 21;

  for ( arg = 1; result && arg < argc; ++arg ) {

    if ( ! strcmp( argv[ arg ], "-files" ) ) {
      ++arg;
      arguments->listFile = argv[ arg ];
    } else if ( ! strcmp( argv[ arg ], "-tmpdir" ) ) {
      ++arg;
      arguments->tmpdir = argv[ arg ];
    } else if ( ! strcmp( argv[ arg ], "-desc" ) ) {
      ++arg;
      arguments->description = argv[ arg ];
    } else if ( ! strcmp( argv[ arg ], "-timestamp" ) ) {
      ++arg;
      arguments->yyyymmddhh = atoi( argv[ arg ] );
      result = isValidYYYYMMDDHH( arguments->yyyymmddhh );
    } else if ( ! strcmp( argv[ arg ], "-hours" ) ) {
      ++arg;
      arguments->hours = atoi( argv[ arg ] );
      result = arguments->hours > 0;
    } else if ( ! strcmp( argv[ arg ], "-variable" ) ) {
      ++arg;
      arguments->variable = argv[ arg ];
      result = arguments->variable[ 0 ] != '\0';
    } else if ( ! strcmp( argv[ arg ], "-domain" ) ) {
      double value = 0.0;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->domain[ LONGITUDE ][ MINIMUM ] = value;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->domain[ LATITUDE ][ MINIMUM ] = value;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->domain[ LONGITUDE ][ MAXIMUM ] = value;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->domain[ LATITUDE ][ MAXIMUM ] = value;
      result = isValidBounds( (const double (*)[2]) arguments->domain );
    } else if ( ! strcmp( argv[ arg ], "-corners" ) ) {
      arguments->corners = 1;
    } else if ( ! strcmp( argv[ arg ], "-minimumQuality" ) ) {
      ++arg;

      if ( arg < argc &&
           ( ! strcmp( argv[ arg ], "medium" ) ||
             ! strcmp( argv[ arg ], "high" ) ) ) {
        arguments->minimumQuality = argv[ arg ];
      } else {
        result = 0;
      }
    } else {
      result = 0;
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nInvalid/insufficient command-line arguments.\n" );
  }

  return result;
}



/******************************************************************************
PURPOSE: readData - Read swath data from each listed VIIRS file and
         write the lon-lat subset of data to the temporary file.
INPUTS:  Data* data  Data description to read.
OUTPUTS: Data* data  data->yyyydddhhmm, points, scans, ok, tempFile = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  char* listFileContent = readListFileAndAllocateTimestampsAndPoints( data );
  int wroteSomeData = 0;
  data->ok = listFileContent != 0;

  if ( data->ok ) {
    size_t rows          = 0;
    size_t columns       = 0;
    size_t size          = 0;
    double* buffer       = 0;
    double* longitudes   = 0;
    double* latitudes    = 0;
    double* values       = 0;
    double* longitudesSW = 0;
    double* longitudesSE = 0;
    double* longitudesNW = 0;
    double* longitudesNE = 0;
    double* latitudesSW  = 0;
    double* latitudesSE  = 0;
    double* latitudesNW  = 0;
    double* latitudesNE  = 0;
    unsigned char* mask  = 0;
    char* fileName = listFileContent;

    /* Get each line of list file. It is the VIIRS data file to read: */

    do {
      char* const newline = strchr( fileName, '\n' );

      if ( newline ) {
        *newline = '\0';

        {
          int file = 0;
          long long yyyymmddhhmm = 0;
          int changedDimensions = 0;
          const int inSubset =
            fileInDomain( data, fileName,
                          &file, &yyyymmddhhmm, &rows, &columns, &size,
                          &changedDimensions );

          if ( inSubset ) {

            if ( changedDimensions ) {
              const int corners = arguments->corners;
              const size_t variables = 3 + 8 * corners;
              const size_t maskSize = size * sizeof (unsigned char);
              const size_t dataSize = variables * size * sizeof (double);
              const size_t bytes = dataSize + maskSize;

              if ( buffer ) {
                free( buffer );
                buffer = 0;
              }

              buffer = malloc( bytes );
              data->ok = buffer != 0;

              if ( buffer ) {
                longitudes   = buffer;
                latitudes    = longitudes + size;
                values       = latitudes  + size;
                longitudesSW = corners ? values       + size : 0;
                longitudesSE = corners ? longitudesSW + size : 0;
                longitudesNW = corners ? longitudesSE + size : 0;
                longitudesNE = corners ? longitudesNW + size : 0;
                latitudesSW  = corners ? longitudesNE + size : 0;
                latitudesSE  = corners ? latitudesSW  + size : 0;
                latitudesNW  = corners ? latitudesSE  + size : 0;
                latitudesNE  = corners ? latitudesNW  + size : 0;
                mask = (unsigned char*)
                  ( corners ? latitudesNE + size : values + size );
              } else {
                fprintf( stderr,
                       "\nCan't allocate %lu bytes "
                       "to complete the requested action.\n", bytes );
              }
            }

            if ( data->ok ) {
              readCoordinatesAndValues( data, file, rows, columns,
                                        longitudes, latitudes, values );
            }

            closeFile( file );
            file = -1;

            if ( data->ok ) {
              writeDataSubset( data, yyyymmddhhmm, rows, columns,
                               longitudes, latitudes, values, mask,
                               longitudesSW, longitudesSE,
                               longitudesNW, longitudesNE,
                               latitudesSW, latitudesSE,
                               latitudesNW, latitudesNE );

              if ( data->ok ) {
                wroteSomeData = 1;
              }
            }
          }

          fileName = newline + 1;
        }
      } else {
        fileName = 0;
      }
    } while ( fileName ); /* End loop on listFile. */

    if ( listFileContent ) {
      free( listFileContent );
      listFileContent = 0;
    }

    if ( data->tempFile ) { /* Done writing to temp file so close it: */
      fclose( data->tempFile );
      data->tempFile = 0;
    }
  }

  data->ok = wroteSomeData;
}



/******************************************************************************
PURPOSE: readListFileAndAllocateTimestampsAndPoints - Read list file and
         return its contents as a string and allocate timestamps and points
         arrays with length equal to lines in the list file.
INPUTS:  Data* const data  data->arguments.listFile.
OUTPUTS: Data* const data  data->ok, yyyydddhhmm, points allocated.
RETURNS: char* Allocated string containing lines of list file or 0 if failed.
******************************************************************************/

static char* readListFileAndAllocateTimestampsAndPoints( Data* const data ) {
  char* result = 0;
  assert( data ); assert( data->arguments.listFile );
  assert( data->arguments.listFile[ 0 ] );
  assert( data->yyyydddhhmm == 0 );
  assert( data->points == 0 );

  {
    result = readFile( data->arguments.listFile );
    data->ok = result != 0;

    if ( result ) {
      allocateTimestampsAndPoints( result, data );

      if ( ! data->ok ) {
        free( result );
        result = 0;
      }
    }
  }

  assert( ( result &&
            data->ok == 1 && data->yyyydddhhmm != 0 && data->points != 0 ) ||
          ( ! result &&
            data->ok == 0 && data->yyyydddhhmm == 0 && data->points == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: allocateTimestampsAndPoints - Allocate arrays for timestamps
         and points per scan.
INPUTS:  const char* const listFileContent  String content of list file.
         Data* const data                   Data to allocate.
OUTPUTS: Data* const data  data->ok, data->yyyydddhhmm, points allocated.
******************************************************************************/

static void allocateTimestampsAndPoints( const char* const listFileContent,
                                         Data* const data ) {
  assert( listFileContent );
  assert( data );
  assert( data->arguments.listFile );
  assert( data->yyyydddhhmm == 0 );
  assert( data->points == 0 );

  {
    const size_t lines = linesInString( listFileContent );
    data->ok = lines != 0;

    if ( ! lines ) {
      fprintf(stderr, "\nInvalid list file '%s'.\n", data->arguments.listFile);
    } else {
      const size_t bytes = lines * sizeof (long long);
      data->yyyydddhhmm = malloc( bytes );
      data->points = data->yyyydddhhmm ? malloc( bytes ) : 0;
      data->ok = data->points != 0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nCan't allocate %lu bytes to complete the "
                 "requested action.\n", bytes );

        if ( data->yyyydddhhmm ) {
          free( data->yyyydddhhmm );
          data->yyyydddhhmm = 0;
        }
      }
    }
  }

  assert( ( data->ok == 1 && data->yyyydddhhmm != 0 && data->points != 0 ) ||
          ( data->ok == 0 && data->yyyydddhhmm == 0 && data->points == 0 ) );
}



/******************************************************************************
PURPOSE: fileInDomain - Is the file in domain? If so get its timestamp and
         dimensions.
INPUTS:  Data* const data              data->arguments.domain.
         const char* const fileName    Name of file to check bounds/dims.
         size_t* const rows            Rows of data in previous file or 0.
         size_t* const columns         Columns of data in previous file or 0.
OUTPUTS: Data* const data              data->ok.
         int* const file               NetCDF file id of file.
         long long* const yyyymmddhhmm Timestamp of file.
         size_t* const rows            Rows of data in file.
         size_t* const columns         Columns of data in file.
         size_t* const size            rows * columns.
         int* const changedDimensions  1 if rows or columns changed.
RETURNS: int 1 if in domain else 0.
******************************************************************************/

static int fileInDomain( Data* const data,
                         const char* const fileName,
                         int* const file,
                         long long* const yyyymmddhhmm,
                         size_t* const rows,
                         size_t* const columns,
                         size_t* const size,
                         int* const changedDimensions ) {

  int result = 0;
  assert( data );
  assert( fileName ); assert( *fileName ); assert( file );
  assert( yyyymmddhhmm );
  assert( rows ); assert( columns ); assert( size );
  assert( changedDimensions );

  *yyyymmddhhmm = fileTimestamp( fileName );
  data->ok = *yyyymmddhhmm != 0;
  *changedDimensions = 0;

  if ( data->ok ) {
    *file = openFile( fileName );
    data->ok = *file != -1;

    if ( *file != -1 ) {
      result = swathInDomain( *file,
                              (const double (*)[2]) data->arguments.domain );

      if ( result ) {
        size_t rows2 = 0;
        size_t columns2 = 0;
        data->ok = readFileDimensions( *file, &rows2, &columns2 );

        if ( data->ok ) {

          if ( rows2 != *rows || columns2 != *columns ) {
            *rows = rows2;
            *columns = columns2;
            *size = rows2 * columns2;
            *changedDimensions = 1;
          }
        }
      }
    }
  }

  *size = *rows * *columns;
  result = result && *size;

  if ( ! ( data->ok && result ) ) {
    result = 0;

    if ( *file != -1 ) {
      closeFile( *file );
      *file = -1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readCoordinatesAndValues - Read lon-lats and variable data.
INPUTS:  Data* const data              data->arguments.variable
         const int file                NetCDF file id of file to read.
         const size_t rows             Rows of data to read.
         const size_t columns          Columns of data to read.
OUTPUTS: Data* const data              data->ok.
         double* const longitudes[ rows * columns ]  Longitudes read.
         double* const latitudes[  rows * columns ]  Latitudes read.
         double* const values[     rows * columns ]  Values read.
******************************************************************************/

static void readCoordinatesAndValues( Data* const data,
                                      const int file,
                                      const size_t rows,
                                      const size_t columns,
                                      double* const longitudes,
                                      double* const latitudes,
                                      double* const values ) {

  assert( data ); assert( data->ok ); assert( data->arguments.variable );
  assert( file > -1 );
  assert( rows != 0 ); assert( columns != 0 );
  assert( longitudes ); assert( latitudes ); assert( values );

  data->ok = readFileData( file, "Longitude", 0, rows, columns, longitudes );

  if ( data->ok ) {
    data->ok = readFileData( file, "Latitude", 0, rows, columns, latitudes );

    if ( data->ok ) {
      const int allowMediumQuality =
        ! strcmp( data->arguments.minimumQuality, "medium" );
      data->ok = readFileData( file, data->arguments.variable,
                               allowMediumQuality, rows, columns, values );
    }
  }
}



/******************************************************************************
PURPOSE: fileTimestamp - Timestamp of file.
INPUTS:  const char* const fileName   Name of VIIRS file.
RETURNS: long long yyyymmddhhmm of file or 0 if invalid and message on stderr.
******************************************************************************/

static long long fileTimestamp( const char* const fileName ) {
  long long result = 0;
  const char* const tag = "_s";
  const char* s = strstr( fileName, tag );

  if ( s ) {
    int c = 0;
    s += strlen( tag );

    for ( c = 0; *s && c < 12; ++c, ++s ) {
      result = result * 10 + *s - '0';
    }
  }

  if ( ! isValidYYYYMMDDHHMM( result ) ) {
    fprintf( stderr, "\nInvalid file name timestamp '%s'.\n", fileName );
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: writeDataSubset - Write subset of data (and corners) to temp file.
INPUTS:  Data* const data                     Data description.
         const long long yyyymmddhhmm         File timestamp.
         const size_t rows                    Rows of input data points.
         const size_t columns                 Columns of input data points.
         const double longitudes[   points ]  Longitudes to write.
         const double latitudes[    points ]  Latitudes to write.
         const double values[       points ]  Values to write.
OUTPUTS: const unsigned char mask[  points ]  Subset mask.
         double longitudesSW[ points ]  0 or computed corners to written.
         double longitudesSE[ points ]  0 or computed corners to written.
         double longitudesNW[ points ]  0 or computed corners to written.
         double longitudesNE[ points ]  0 or computed corners to written.
         double latitudesSW[  points ]  0 or computed corners to written.
         double latitudesSE[  points ]  0 or computed corners to written.
         double latitudesNW[  points ]  0 or computed corners to written.
         double latitudesNE[  points ]  0 or computed corners to written.
OUTPUTS: Data* data  data->ok.
******************************************************************************/

static void writeDataSubset( Data* const data,
                             const long long yyyymmddhhmm,
                             const size_t rows,
                             const size_t columns,
                             const double longitudes[],
                             const double latitudes[],
                             const double values[],
                             unsigned char mask[],
                             double longitudesSW[],
                             double longitudesSE[],
                             double longitudesNW[],
                             double longitudesNE[],
                             double latitudesSW[],
                             double latitudesSE[],
                             double latitudesNW[],
                             double latitudesNE[] ) {

  assert( data ); assert( data->ok );
  assert( isValidYYYYMMDDHHMM( yyyymmddhhmm ) );
  assert( rows != 0 ); assert( columns != 0 ); assert( mask );
  assert( longitudes ); assert( latitudes ); assert( values );

  {
    const size_t points = rows * columns;
    const size_t subsetPoints =
      pointsInDomain( (const double (*)[2]) data->arguments.domain,
                      points, longitudes, latitudes, values, mask );

    if ( subsetPoints ) {

      if ( longitudesSW ) {
        computeCorners( rows, columns, longitudes, latitudes,
                        longitudesSW, longitudesSE,
                        longitudesNW, longitudesNE,
                        latitudesSW, latitudesSE,
                        latitudesNW, latitudesNE );
      }

      writeSubset( data, yyyymmddhhmm, subsetPoints, points,
                   mask, longitudes, latitudes, values,
                   longitudesSW, longitudesSE,
                   longitudesNW, longitudesNE,
                   latitudesSW, latitudesSE,
                   latitudesNW, latitudesNE );
    }
  }
}



/******************************************************************************
PURPOSE: writeSubset - Store timestamps and subset point counts and write
         subset of data to temp file.
INPUTS:  Data* const data                     Data description.
         const long long yyyymmddhhmm         File timestamp.
         const size_t subsetPoints            Points in subset domain.
         const size_t points                  Input data points.
         const unsigned char mask[  points ]  Subset mask.
         const double longitudes[   points ]  Longitudes to write.
         const double latitudes[    points ]  Latitudes to write.
         const double values[       points ]  Values to write.
         const double longitudesSW[ points ]  0 or corners to write.
         const double longitudesSE[ points ]  0 or corners to write.
         const double longitudesNW[ points ]  0 or corners to write.
         const double longitudesNE[ points ]  0 or corners to write.
         const double latitudesSW[  points ]  0 or corners to write.
         const double latitudesSE[  points ]  0 or corners to write.
         const double latitudesNW[  points ]  0 or corners to write.
         const double latitudesNE[  points ]  0 or corners to write.
OUTPUTS: Data* data  data->ok.
******************************************************************************/

static void writeSubset( Data* const data,
                         const long long yyyymmddhhmm,
                         const size_t subsetPoints,
                         const size_t points,
                         const unsigned char mask[],
                         const double longitudes[],
                         const double latitudes[],
                         const double values[],
                         const double longitudesSW[],
                         const double longitudesSE[],
                         const double longitudesNW[],
                         const double longitudesNE[],
                         const double latitudesSW[],
                         const double latitudesSE[],
                         const double latitudesNW[],
                         const double latitudesNE[] ) {

  assert( data );
  assert( isValidYYYYMMDDHHMM( yyyymmddhhmm ) );
  assert( subsetPoints != 0 ); assert( points != 0 ); assert( mask );
  assert( longitudes ); assert( latitudes ); assert( values );

  /* Open temp file for writing if it does not yet exist: */

  if ( data->tempFile == 0 ) {
    const int pid = getpid();
    memset( data->tempFileName, 0, sizeof (FileName) );
    snprintf( data->tempFileName,
              sizeof (FileName) / sizeof (char) - 1,
              "%s/%s.%04d", data->arguments.tmpdir, TEMP_FILE_NAME, pid );
    data->tempFile = fopen( data->tempFileName, "wb" );

    if ( ! data->tempFile ) {
      fprintf( stderr, "\nCan't create temporary output file '%s'.\n",
               data->tempFileName );
      data->ok = 0;
    }

  }

  if ( data->ok ) { /* Append to arrays of timestamps and points: */
    const long long yyyydddhhmm = convertTimestamp( yyyymmddhhmm );
    data->yyyydddhhmm[ data->scans ] = yyyydddhhmm;
    data->points[ data->scans ] = (long long) subsetPoints;
    data->scans += 1;
  }

  if ( data->ok ) { /* Write subset data to temp file: */
    const size_t variables = 3 + 8 * data->arguments.corners;
    const size_t count = variables * subsetPoints;
    const size_t bytes = count * sizeof (double);
    double* buffer = malloc( bytes );
    data->ok = buffer != 0;

    if ( ! buffer ) {
      fprintf( stderr,
              "\nCan't allocate %lu bytes to complete the requested action.\n",
               bytes );
    } else {
      const int corners = longitudesSW != 0;
      double* const subsetLongitudes = buffer;
      double* const subsetLatitudes  = subsetLongitudes + subsetPoints;
      double* const subsetValues     = subsetLatitudes  + subsetPoints;
      double* const subsetLongitudesSW =
        corners ? subsetValues + subsetPoints : 0;
      double* const subsetLongitudesSE =
        corners ? subsetLongitudesSW + subsetPoints : 0;
      double* const subsetLongitudesNW =
        corners ? subsetLongitudesSE + subsetPoints : 0;
      double* const subsetLongitudesNE =
        corners ? subsetLongitudesNW + subsetPoints : 0;
      double* const subsetLatitudesSW =
        corners ? subsetLongitudesNE + subsetPoints : 0;
      double* const subsetLatitudesSE =
        corners ? subsetLatitudesSW + subsetPoints : 0;
      double* const subsetLatitudesNW =
        corners ? subsetLatitudesSE + subsetPoints : 0;
      double* const subsetLatitudesNE =
        corners ? subsetLatitudesNW + subsetPoints : 0;
      size_t outputPoints = 0;
      size_t point = 0;

      for ( point = 0; point < points; ++point ) {
        const int m = mask[ point ];

        if ( m ) {
          subsetLongitudes[ outputPoints ] = longitudes[ point ];
          subsetLatitudes[  outputPoints ] = latitudes[  point ];
          subsetValues[     outputPoints ] = values[     point ];

          if ( corners ) {
            subsetLongitudesSW[ outputPoints ] = longitudesSW[ point ];
            subsetLongitudesSE[ outputPoints ] = longitudesSE[ point ];
            subsetLongitudesNW[ outputPoints ] = longitudesNW[ point ];
            subsetLongitudesNE[ outputPoints ] = longitudesNE[ point ];
            subsetLatitudesSW[  outputPoints ] = latitudesSW[  point ];
            subsetLatitudesSE[  outputPoints ] = latitudesSE[  point ];
            subsetLatitudesNW[  outputPoints ] = latitudesNW[  point ];
            subsetLatitudesNE[  outputPoints ] = latitudesNE[  point ];
          }

          ++outputPoints;
        }
      }

      assert( outputPoints == subsetPoints );
      rotate8ByteArrayIfLittleEndian( buffer, count );

      data->ok =
        fwrite( buffer, sizeof *buffer, count, data->tempFile ) == count;

      if ( ! data->ok ) {
        fprintf( stderr, "\nFailed to write subset data to temp file '%s'.\n",
                 data->tempFileName );
      }

      free( buffer );
      buffer = 0;
    }
  }
}



/******************************************************************************
PURPOSE: streamData - Write ASCII header and XDR binary data (content of temp
         file) to stdout.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* const data  data->ok, tempFile = 0 (closed and removed).
******************************************************************************/

static void streamData( Data* const data ) {
  const size_t bytes = 256 * 1024 * 1024;
  void* buffer = malloc( bytes );
  assert( data );
  assert( data->tempFileName[ 0 ] );
  assert( data->tempFile == 0 ); /* Temp file is closed after writing it. */
  data->ok = buffer != 0;

  if ( ! buffer ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes to complete the requested action.\n",
             bytes );
  } else {
    data->tempFile = fopen( data->tempFileName, "rb" );
    data->ok = data->tempFile != 0;

    if ( ! data->ok ) {
      fprintf( stderr, "\nCan't open temp data file '%s' for reading.\n",
               data->tempFileName );
    } else {
      streamHeader( data );
      rotate8ByteArrayIfLittleEndian( data->yyyydddhhmm, data->scans );
      rotate8ByteArrayIfLittleEndian( data->points, data->scans );
      data->ok =
        fwrite( data->yyyydddhhmm, sizeof data->yyyydddhhmm[ 0 ], data->scans,
                stdout ) == data->scans;
      data->ok = data->ok &&
        fwrite( data->points, sizeof data->points[ 0 ], data->scans, stdout )
        == data->scans;
      rotate8ByteArrayIfLittleEndian( data->yyyydddhhmm, data->scans );
      rotate8ByteArrayIfLittleEndian( data->points, data->scans );

      while ( data->ok && ! feof( data->tempFile ) ) {
        const size_t bytesRead = fread( buffer, 1, bytes, data->tempFile );

        if ( bytesRead ) {
          const size_t bytesWritten = fwrite( buffer, 1, bytesRead, stdout );
          data->ok = bytesWritten == bytesRead;
        }
      }
    }

    free( buffer );
    buffer = 0;
  }

  if ( ! data->ok ) {
    fprintf( stderr, "\nFailed to stream subset data from temp file '%s'.\n",
             data->tempFileName );
  }

  if ( data->tempFile ) {
    fclose( data->tempFile );
    data->tempFile = 0;
  }

  unlink( data->tempFileName );
}



/******************************************************************************
PURPOSE: streamHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
******************************************************************************/

static void streamHeader( const Data* const data ) {
  const Arguments* const arguments = &data->arguments;
  const int variables = 3 + arguments->corners * 8;
  const int yyyymmddhh = arguments->yyyymmddhh;
  const int yyyy       = yyyymmddhh / 1000000;
  const int mm         = yyyymmddhh / 10000 % 100;
  const int dd         = yyyymmddhh / 100 % 100;
  const int hh         = yyyymmddhh % 100;

  printf( "Swath 2.0\n%s\n%04d-%02d-%02dT%02d:00:00-0000\n",
          arguments->description, yyyy, mm, dd, hh );
  printf( "# Dimensions: variables timesteps scans:\n%d %d %d\n",
          variables, arguments->hours, data->scans );
  printf( "# Variable names:\n" );
  printf( "Longitude Latitude %s", arguments->variable );

  if ( variables == 11 ) {
    printf( " Longitude_SW Longitude_SE Longitude_NW Longitude_NE"
            " Latitude_SW Latitude_SE Latitude_NW Latitude_NE" );
  }

  printf( "\n# Variable units:\ndeg deg -" );

  if ( variables == 11 ) {
    printf( " deg deg deg deg deg deg deg deg" );
  }

  printf( "\n# Domain: <min_lon> <min_lat> <max_lon> <max_lat>\n%g %g %g %g\n",
          arguments->domain[ LONGITUDE ][ MINIMUM ],
          arguments->domain[ LATITUDE  ][ MINIMUM ],
          arguments->domain[ LONGITUDE ][ MAXIMUM ],
          arguments->domain[ LATITUDE  ][ MAXIMUM ] );
  printf( "# MSB 64-bit integers (yyyydddhhmm)" );
  printf( " timestamps[scans] and\n" );
  printf( "# MSB 64-bit integers points[scans] and\n" );
  printf( "# IEEE-754 64-bit reals" );
  printf( " data_1[variables][points_1] ..." );
  printf( " data_S[variables][points_S]:\n" );
}




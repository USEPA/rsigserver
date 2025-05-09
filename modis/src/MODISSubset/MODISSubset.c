/******************************************************************************
PURPOSE: MODISSubset.c - Extract a lon-lat subset of data from a list of
         MODIS HDF-EOS files and write it to stdout as XDR binary format.

NOTES:   Uses HDF-EOS libraries and libs they depend on (z, jpeg).
         Compile:
         gcc -Wall -g -o MODISSubset MODISSubset.c Utilities.c ReadData.c \
                   -L../../../lib/$platform \
                   -lhdfeos -lmfhdf -ldf -lsz -lz -ljpeg -lm

         Usage:
         MODISSubset -files <listfile> \
                     -tmpdir <temp_directory> \
                     -desc "description text" \
                     -timestamp <yyyymmddhh> -hours <count> \
                     -variable <name> \
                     -domain <minimum_longitude> <minimum_latitude> \
                             <maximum_longitude> <maximum_latitude> \
                     [-corners]

          Example:
          ../../../bin/$platform/MODISSubset \
          -files testdata/files \
          -tmpdir testdata \
          -desc \
 "https://modwebsrv.modaps.eosdis.nasa.gov/cgi-bin/RSIGservice,MODISSubset" \
          -timestamp 2013061500 -hours 24 \
          -variable Optical_Depth_Land_And_Ocean \
          -domain -75 35 -70 36 \
          -corners \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 14-line ASCII header followed by binary 64-bit big-endian arrays:

Swath 2.0
https://modwebsrv.modaps.eosdis.nasa.gov/cgi-bin/RSIGservice,MODISSubset
2013-06-15T00:00:00-0000
# Dimensions: variables timesteps scans:
11 24 2
# Variable names:
Longitude Latitude Optical_Depth_Land_And_Ocean \
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
20130661848
20130661942
5
122
-7.1847106933593750e+01
-7.1855308532714844e+01
 ...
3.5999182701110840e+01
3.5997957229614258e+01


HISTORY: 2005-11-01 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, printf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <stdlib.h>    /* For malloc(), free(). */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For MISSING_VALUE, IN_RANGE(), LONGITUDE, Bounds. */
#include "ReadData.h"  /* For openFile(), readFileBounds(), readFileData()*/

/*================================= MACROS =================================*/

/* Name of temporary file created in -tmpdir will have PID appended: */

#define TEMP_FILE_NAME "junk_MODISSubset"

/*================================== TYPES ==================================*/

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;    /* File containing list of MODIS files to read. */
  const char* tmpdir;      /* Name of directory to write temp files. */
  const char* description; /* User-supplied description. */
  const char* variable;    /* Name of variable to read. */
  Bounds      domain; /* Subset domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  int         yyyymmddhh;  /* First timestamp of subset. */
  int         hours;       /* Number of hours in subset. */
  int         corners;     /* Compute interpolated lon-lat corner points?*/
} Arguments;

/* Data type: */

typedef struct {
  Arguments  arguments;    /* User-supplied (command-line) arguments. */
  char units[ 80 ];        /* Units of variable. */
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

static long long swathFileTimestamp( const char* const fileName );

static void readFileInfo( Data* const data,
                          const char* const fileName,
                          int* const file,
                          long long* const yyyydddhhmm,
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

static void writeSubsetData( Data* const data,
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
                             double latitudesNE[] );

static void streamData( Data* const data );

static void streamHeader( const Data* const data );

static void streamSwathTimestamps( Data* const data );

static void streamSwathPoints( Data* const data );

static void streamTempFile( Data* const data );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of MODIS files and write
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
    readData( &data ); /* Read subset of MODIS files and write temp files. */

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

static void printUsage( const char* programName ) {
  assert( programName );
  fprintf( stderr, "\a\n\n%s - Read a set of MODIS files and extract swath\n",
           programName );
  fprintf( stderr, "data subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable(s).\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "  -files <file> \\\n" );
  fprintf( stderr, "  -tmpdir <temp_directory> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, "  -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, "  -corners \n\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies where to write temporary files.\n" );
  fprintf( stderr, "-corners option will output 8 additional variables:\n" );
  fprintf( stderr, "  Longitude_SW Longitude_SE Longitude_NW Longitude_NE\n");
  fprintf( stderr, "  Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "that are the linearly interpolated " );
  fprintf( stderr, "(and edge extrapolated)\n" );
  fprintf( stderr, "corner points for each center-pixel point.\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example #1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files testdata/files \\\n" );
  fprintf( stderr, "-tmpdir testdata \\\n" );
  fprintf( stderr, "-desc "
           "\"https://modwebsrv.modaps.eosdis.nasa.gov/cgi-bin/RSIGservice,"
           "MODISSubset\" \\\n");
  fprintf( stderr, "-timestamp 2013061500 -hours 24 \\\n" );
  fprintf( stderr, "-variable Optical_Depth_Land_And_Ocean \\\n" );
  fprintf( stderr, "-domain -126 25 -65 50 > subset.xdr\n\n" );
  fprintf( stderr, "AOD over US on June 15, 2013.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n" );
  fprintf( stderr, "Swath 2.0\n" );
  fprintf( stderr, "http://www.star.nesdis.noaa.gov/smcd/emb/MODIS_aerosol\n");
  fprintf( stderr, "2013-06-15T00:00:00-0000\n" );
  fprintf( stderr, "# Dimensions: variables timesteps scans:\n" );
  fprintf( stderr, "4 24 2\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Longitude Latitude AerosolOpticalDepth_at_555nm\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "deg deg YYYYDDDHHMM -\n" );
  fprintf( stderr, "# Domain: <min_lon> <min_lat> <max_lon> <max_lat>\n" );
  fprintf( stderr, "-126 25 -65 50\n" );
  fprintf( stderr, "# MSB 64-bit integers (yyyydddhhmm)" );
  fprintf( stderr, " timestamps[scans] and\n" );
  fprintf( stderr, "# MSB 64-bit integers points[scans] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals" );
  fprintf( stderr, " data_1[variables][points_1] ..." );
  fprintf( stderr, " data_S[variables][points_S]:\n" );
  fprintf( stderr, "<binary data arrays here>\n\n\n");
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

  result = argc == 18 || argc == 19;

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
PURPOSE: readData - Read swath data from each listed MODIS file and
         write the lon-lat subset of data to the temporary file.
INPUTS:  Data* data  Data description to read.
OUTPUTS: Data* data  data->yyyydddhhmm, points, scans, ok, tempFile = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  const int corners = arguments->corners;
  char* listFileContent = readListFileAndAllocateTimestampsAndPoints( data );
  int wroteSomeData = 0;
  char* fileName       = 0;
  char* end            = 0;
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
  data->ok = 0;

  /* Get each line of list file. It is the MODIS data file to read: */

  for ( fileName = strtok_r( listFileContent, "\n", &end );
        fileName;
        fileName = strtok_r( 0, "\n", &end ) ) {
    int file = 0;
    long long yyyydddhhmm = 0;
    int changedDimensions = 0;

    readFileInfo( data, fileName, &file, &yyyydddhhmm, &rows, &columns,
                  &size, &changedDimensions );

    if ( data->ok ) {
      const size_t points = size;

      if ( changedDimensions ) {
        const size_t variables = 3 + 8 * corners;
        const size_t dataSize = variables * size;
        const size_t bytes = dataSize * sizeof (double);

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

      closeFile( file ), file = -1;

      if ( data->ok ) {
        size_t subsetPoints = 0;

        if ( corners ) {
          computeCorners( rows, columns,
                          longitudes,   latitudes,
                          longitudesSW, longitudesSE,
                          longitudesNW, longitudesNE,
                          latitudesSW,  latitudesSE,
                          latitudesNW,  latitudesNE );
        }

        subsetPoints =
          pointsInSubset( (const double (*)[2]) arguments->domain,
                          points, longitudes, latitudes, values,
                          longitudesSW, longitudesSE,
                          longitudesNW, longitudesNE,
                          latitudesSW,  latitudesSE,
                          latitudesNW,  latitudesNE );

        DEBUG( fprintf( stderr, "subsetPoints = %lu\n", subsetPoints ); )

        if ( subsetPoints ) {

          if ( subsetPoints < points ) {
            compactSubsetData( subsetPoints, points,
                               longitudes, latitudes, values,
                               longitudesSW, longitudesSE,
                               longitudesNW, longitudesNE,
                               latitudesSW,  latitudesSE,
                               latitudesNW,  latitudesNE );
          }

          writeSubsetData( data, subsetPoints,
                           longitudes, latitudes, values,
                           longitudesSW, longitudesSE,
                           longitudesNW, longitudesNE,
                           latitudesSW,  latitudesSE,
                           latitudesNW,  latitudesNE );

          if ( data->ok ) {
            data->yyyydddhhmm[ data->scans ] = yyyydddhhmm;
            data->points[ data->scans ] = subsetPoints;
            DEBUG( fprintf( stderr, "scan %d: %lld %lld\n",
                            wroteSomeData, data->yyyydddhhmm[ data->scans ],
                            data->points[ data->scans ] ); )
            data->scans += 1;
            wroteSomeData = 1;
          }
        }
      }
    }
  } /* End loop on listFile. */

  free( listFileContent );
  listFileContent = 0;

  if ( data->tempFile ) { /* Done writing to temp file so close it: */
    fclose( data->tempFile );
    data->tempFile = 0;
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
    size_t length = 0;
    result = readFile( data->arguments.listFile, &length );
    data->ok = result != 0;

    if ( result ) {
      allocateTimestampsAndPoints( result, data );

      if ( ! data->ok ) {
        free( result );
        result = 0;
      }
    }
  }

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
      }
    }
  }
}



/******************************************************************************
PURPOSE: swathFileTimestamp - Timestamp of swath file.
INPUTS:  const char* const fileName   Name of MODIS file.
RETURNS: long long yyyydddhhmm of file or 0 if failed and message on stderr.
NOTES:   File names look like:
         MOD04_L2.A2006100.0425.005.2006243075120.hdf
         MOD04_3K.A2008185.1815.006.2015029004104.hdf
******************************************************************************/

static long long swathFileTimestamp( const char* const fileName ) {
  long long result = 0;
  const char* const slash = strrchr( fileName, '/' );
  const char* name = slash ? slash + 1 : fileName;
  const char* tag = "_L2.A";
  const char* s = strstr( name, tag );

  if ( s == 0 ) {
    tag = "_3K.A";
    s = strstr( fileName, tag );
  }

  if ( s ) {
    int c = 0;
    s += strlen( tag );

    /* Parse YYYYDDD: */

    for ( c = 0; *s && *s != '.' && c < 7; ++c, ++s ) {
      result = result * 10 + *s - '0';
    }

    /* Parse HHMM: */

    if ( *s == '.' ) {
      ++s;

      for ( c = 0; *s && *s != '.' && c < 4; ++c, ++s ) {
        result = result * 10 + *s - '0';
      }
    }
  }

  if ( ! isValidYYYYDDDHHMM( result ) ) {
    fprintf( stderr, "\nInvalid file name timestamp '%s'.\n", fileName );
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileInfo - Parse file timestamp, open it and read bounds and,
         if in subset, read dimensions.
INPUTS:  Data* const data              data->arguments.domain.
         const char* const fileName    Name of data file to open.
         size_t* const rows            Rows of data in previous file or 0.
         size_t* const columns         Columns of data in previous file or 0.
OUTPUTS: Data* const data              data->ok.
         int* const file               HDF file id of data file.
         long long* const yyyydddhhmm  Timestamp of file (if in subset range).
         size_t* const rows            Rows of data in file.
         size_t* const columns         Columns of data in file.
         size_t* const size            rows * columns.
         int* const changedDimensions  1 if rows or columns changed.
******************************************************************************/

static void readFileInfo( Data* const data,
                          const char* const fileName,
                          int* const file,
                          long long* const yyyydddhhmm,
                          size_t* const rows,
                          size_t* const columns,
                          size_t* const size,
                          int* const changedDimensions ) {

  assert( data );
  assert( isValidBounds( (const double (*)[2]) data->arguments.domain ) );
  assert( fileName ); assert( *fileName );  assert( file );
  assert( yyyydddhhmm ); assert( rows ); assert( columns ); assert( size );
  assert( changedDimensions );

  *file = -1;
  *yyyydddhhmm = swathFileTimestamp( fileName );
  data->ok = *yyyydddhhmm != 0;
  *changedDimensions = 0;

  if ( data->ok ) {
    const Arguments* const arguments = &( data->arguments );
    const long long firstTimestamp =
      convertTimestamp( arguments->yyyymmddhh * 100LL );
    const long long lastTimestamp =
      offsetTimestamp( firstTimestamp, arguments->hours );
    data->ok = IN_RANGE( *yyyydddhhmm, firstTimestamp, lastTimestamp );

    if ( data->ok ) {
      *file = openFile( fileName );
      data->ok = *file != -1;

      if ( data->ok ) {
        Bounds bounds = { { -180.0, 180.0 }, { -90.0, 90.0 } };
        data->ok = readFileBounds( *file, bounds );

        if ( data->ok ) {
          data->ok =
            boundsOverlap( (const double (*)[2]) bounds,
                           (const double (*)[2]) arguments->domain );

          if ( data->ok ) {
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
    }
  }

  if ( ! data->ok ) {

    if ( *file != -1 ) {
      closeFile( *file ), *file = -1;
    }
  }
}



/******************************************************************************
PURPOSE: readCoordinatesAndValues - Read lon-lats and variable data.
INPUTS:  Data* const data              data->arguments.variable
         const int file                Data file id of file to read.
         const size_t rows             Rows of data to read.
         const size_t columns          Columns of data to read.
OUTPUTS: Data* const data              data->ok, data->units[ 80 ]              Units of variable.
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

  char unused[ 80 ] = "";

  assert( data ); assert( data->ok ); assert( data->arguments.variable );
  assert( file > -1 );
  assert( rows != 0 ); assert( columns != 0 );
  assert( longitudes ); assert( latitudes ); assert( values );

  data->ok =
    readFileData( file, "Longitude", rows, columns, unused, longitudes );

  if ( data->ok ) {
    data->ok =
      readFileData( file, "Latitude", rows, columns, unused, latitudes ) ;

    if ( data->ok ) {
      data->ok =
        clampInvalidCoordinates( rows * columns, longitudes, latitudes );

      if ( data->ok ) {
        data->ok = readFileData( file, data->arguments.variable, rows, columns,
                                 data->units, values );
      }
    }
  }
}



/******************************************************************************
PURPOSE: writeSubsetData - Write subset of data to temp file.
INPUTS:  Data* const data                     Data description.
         const size_t points                  Points in subset domain.
         double longitudes[   points ]  Longitudes to write.
         double latitudes[    points ]  Latitudes to write.
         double values[       points ]  Values to write.
         double longitudesSW[ points ]  0 or corners to write.
         double longitudesSE[ points ]  0 or corners to write.
         double longitudesNW[ points ]  0 or corners to write.
         double longitudesNE[ points ]  0 or corners to write.
         double latitudesSW[  points ]  0 or corners to write.
         double latitudesSE[  points ]  0 or corners to write.
         double latitudesNW[  points ]  0 or corners to write.
         double latitudesNE[  points ]  0 or corners to write.
OUTPUTS: Data* data  data->ok, ruined arrays.
******************************************************************************/

static void writeSubsetData( Data* const data,
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

  assert( data ); assert( data->ok );
  assert( points );
  assert( longitudes ); assert( latitudes ); assert( values );
  assert( IMPLIES_ELSE( longitudesSW,
                        NON_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                   latitudesSW, latitudesSE,
                                   latitudesNW, latitudesNE ),
                        IS_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                  latitudesSW, latitudesSE,
                                  latitudesNW, latitudesNE ) ) );

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

  if ( data->ok ) {
    rotate8ByteArrayIfLittleEndian( longitudes, points );
    rotate8ByteArrayIfLittleEndian( latitudes, points );
    rotate8ByteArrayIfLittleEndian( values, points );

    if ( longitudesSW ) {
      rotate8ByteArrayIfLittleEndian( longitudesSW, points );
      rotate8ByteArrayIfLittleEndian( longitudesSE, points );
      rotate8ByteArrayIfLittleEndian( longitudesNW, points );
      rotate8ByteArrayIfLittleEndian( longitudesNE, points );
      rotate8ByteArrayIfLittleEndian( latitudesSW, points );
      rotate8ByteArrayIfLittleEndian( latitudesSE, points );
      rotate8ByteArrayIfLittleEndian( latitudesNW, points );
      rotate8ByteArrayIfLittleEndian( latitudesNE, points );
    }

    data->ok =
      fwrite(longitudes, sizeof *longitudes, points, data->tempFile) == points;

    if ( data->ok ) {
      data->ok =
        fwrite(latitudes, sizeof *latitudes, points, data->tempFile) == points;

      if ( data->ok ) {
        data->ok =
          fwrite( values, sizeof *values, points, data->tempFile ) == points;

        if ( data->ok ) {

          if ( longitudesSW ) {
            data->ok =
              fwrite( longitudesSW, sizeof *longitudesSW, points,
                      data->tempFile ) == points;

            if ( data->ok ) {
              data->ok =
                fwrite( longitudesSE, sizeof *longitudesSE, points,
                        data->tempFile ) == points;

              if ( data->ok ) {
                data->ok =
                  fwrite( longitudesNW, sizeof *longitudesNW, points,
                          data->tempFile ) == points;

                if ( data->ok ) {
                  data->ok =
                    fwrite( longitudesNE, sizeof *longitudesNE, points,
                            data->tempFile ) == points;

                  if ( data->ok ) {
                    data->ok =
                      fwrite( latitudesSW, sizeof *latitudesSW, points,
                              data->tempFile ) == points;

                    if ( data->ok ) {
                      data->ok =
                        fwrite( latitudesSE, sizeof *latitudesSE, points,
                                data->tempFile ) == points;

                      if ( data->ok ) {
                        data->ok =
                          fwrite( latitudesNW, sizeof *latitudesNW, points,
                                  data->tempFile ) == points;

                        if ( data->ok ) {
                          data->ok =
                            fwrite( latitudesNE, sizeof *latitudesNE, points,
                                    data->tempFile ) == points;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: streamData - Write XDR format subset data to stdout.
INPUTS:  Data* data     Data to write to stream.
******************************************************************************/

static void streamData( Data* const data ) {

  assert( data ); assert( data->ok ); assert( ! data->tempFile );

  streamHeader( data );
  streamSwathTimestamps( data );

  if ( data->ok ) {
    streamSwathPoints( data );

    if ( data->ok ) {
      streamTempFile( data );
    }
  }
}



/******************************************************************************
PURPOSE: streamHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
******************************************************************************/

static void streamHeader( const Data* const data ) {

  assert( data ); assert( data->ok );

  {
    const Arguments* const arguments = &data->arguments;
    const int variables = 3 + arguments->corners * 8;
    const int scans = data->scans;
    const int yyyymmddhh = arguments->yyyymmddhh;
    const int yyyy       = yyyymmddhh / 1000000;
    const int mm         = yyyymmddhh / 10000 % 100;
    const int dd         = yyyymmddhh / 100 % 100;
    const int hh         = yyyymmddhh % 100;

    printf( "Swath 2.0\n%s\n%04d-%02d-%02dT%02d:00:00-0000\n",
            arguments->description, yyyy, mm, dd, hh );
    printf( "# Dimensions: variables timesteps scans:\n%d %d %d\n",
            variables, arguments->hours, scans );
    printf( "# Variable names:\n" );
    printf( "Longitude Latitude %s", arguments->variable );

    if ( variables > 3 ) {
      printf( " Longitude_SW Longitude_SE Longitude_NW Longitude_NE"
              " Latitude_SW Latitude_SE Latitude_NW Latitude_NE" );
    }

    printf( "\n# Variable units:\ndeg deg %s", data->units );

    if ( variables > 3 ) {
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
}



/******************************************************************************
PURPOSE: streamSwathTimestamps - Write MSB 64-bit integer swath timestamps to
         stdout.
INPUTS:  Data* const data  Data to stream.
OUTPUTS: Data* const data  data->ok.
******************************************************************************/

static void streamSwathTimestamps( Data* const data ) {

  assert( data ); assert( data->scans > 0 );
  assert( data->yyyydddhhmm );
  assert( isValidYYYYDDDHHMM( data->yyyydddhhmm[ 0 ] ) );
  assert( isValidYYYYDDDHHMM( data->yyyydddhhmm[ data->scans - 1 ] ) );

  rotate8ByteArrayIfLittleEndian( data->yyyydddhhmm, data->scans );
  data->ok =
    fwrite( data->yyyydddhhmm, sizeof data->yyyydddhhmm[ 0 ], data->scans,
            stdout ) == data->scans;
  rotate8ByteArrayIfLittleEndian( data->yyyydddhhmm, data->scans );

  if ( ! data->ok ) {
    fprintf( stderr, "\a\nFailed to stream subset swath timestamps.\n" );
  }
}



/******************************************************************************
PURPOSE: streamSwathPoints - Write MSB 64-bit integer swath subset point counts
INPUTS:  Data* const data  Data to stream.
OUTPUTS: Data* const data  data->ok.
******************************************************************************/

static void streamSwathPoints( Data* const data ) {

  assert( data ); assert( data->scans > 0 );
  assert( data->points );
  assert( data->points[ 0 ] > 0 );
  assert( data->points[ data->scans - 1 ] > 0 );

  rotate8ByteArrayIfLittleEndian( data->points, data->scans );
  data->ok =
    fwrite( data->points, sizeof data->points[ 0 ], data->scans, stdout )
    == data->scans;
  rotate8ByteArrayIfLittleEndian( data->points, data->scans );

  if ( ! data->ok ) {
    fprintf( stderr, "\a\nFailed to stream subset swath point counts.\n" );
  }
}



/******************************************************************************
PURPOSE: streamTempFile - Stream XDR binary subset data (content of tempfile)
         to stdout.
INPUTS:  Data* const data  data->tempFileName.
OUTPUTS: Data* const data  data->ok, tempFile = 0 (closed).
******************************************************************************/

static void streamTempFile( Data* const data ) {

  const size_t bytes = 256 * 1024 * 1024;
  void* buffer = malloc( bytes );
  data->ok = buffer != 0;

  assert( data ); assert( data->ok ); assert( data->tempFileName[ 0 ] );
  assert( data->tempFile == 0 ); /* tempFile closed after writing to it. */


  if ( ! buffer ) {
    fprintf( stderr,
            "\a\nCan't allocate %lu bytes to complete the requested action.\n",
             bytes );
  } else {
    data->tempFile = fopen( data->tempFileName, "rb" );
    data->ok = data->tempFile != 0;

    if ( ! data->ok ) {
      fprintf( stderr, "\a\nCan't open temp data file '%s' for reading.\n",
                       data->tempFileName );
    } else {

      do {
        const size_t bytesRead = fread( buffer, 1, bytes, data->tempFile );

        if ( bytesRead ) {
          const size_t bytesWritten = fwrite( buffer, 1, bytesRead, stdout );
          data->ok = bytesWritten == bytesRead;
        }
      } while ( data->ok && ! feof( data->tempFile ) );
    }

    free( buffer ), buffer = 0;
  }

  if ( ! data->ok ) {
    fprintf( stderr, "\a\nFailed to stream subset data from temp file '%s'.\n",
             data->tempFileName );
  }

  if ( data->tempFile ) {
    fclose( data->tempFile );
    data->tempFile = 0;
  }
}




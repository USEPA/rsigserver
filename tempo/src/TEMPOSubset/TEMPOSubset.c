/******************************************************************************
PURPOSE: TEMPOSubset.c - Extract a lon-lat subset of data from a list of
         TEMPO NetCDF4 files and write it to stdout as XDR binary format.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).
         Compile:
         gcc -Wall -g -o TEMPOSubset TEMPOSubset.c Utilities.c ReadData.c \
                   -I../../../include/NetCDF4 \
                   -L../../../lib/$platform \
                   -lnetcdf4 -lhdf5_hl -lhdf5 -lcurl -lz -ldl -lm -lc

         Usage:
         TEMPOSubset -files <listfile> \
                      -tmpdir <temp_directory> \
                      -desc "description text" \
                      -timestamp <yyyymmddhh> -hours <count> \
                      -variable <name> \
                      -domain <minimum_longitude> <minimum_latitude> \
                              <maximum_longitude> <maximum_latitude> \
                      [-minimumQuality value] \
                      [-maximumCloudFraction value] \
                      [-maximumSolarZenithAngle value] \
                      [-allowNegativeCounts] \
                      [-corners]

          Example:
          ../../../bin/$platform/TEMPOSubset \
          -files testdata/file_list \
          -tmpdir testdata \
          -desc "http://tempo.si.edu/,TEMPOSubset" \
          -timestamp 2023101700 -hours 24 \
          -variable no2 \
          -domain -130 20 -60 50 \
          -corners \
          > subset.xdr

          Outputs to stdout a data stream in the following formats.
          For L2 data,
          a 14-line ASCII header followed by binary 64-bit big-endian arrays:

Swath 2.0
http://tempo.si.edu/,TEMPOSubset
2023-10-17T00:00:00-0000
# Dimensions: variables timesteps scans:
11 24 2
# Variable names:
Longitude Latitude no2 \
 LongitudeSW LongitudeSE LongitudeNW LongitudeNE \
 LatitudeSW LatitudeSE LatitudeNW LatitudeNE
# Variable units:
deg deg molecules/cm2 deg deg deg deg deg deg deg deg
# Domain: <min_lon> <min_lat> <max_lon> <max_lat>
-130 20 -60 50
# MSB 64-bit integers (yyyydddhhmm) timestamps[scans] and
# MSB 64-bit integers points[scans] and
# IEEE-754 64-bit reals data_1[variables][points_1] ... \
 data_S[variables][points_S]:
<big-endian binary format arrays>
20232831659
20232831705
5
122
-7.1847106933593750e+01
-7.1855308532714844e+01
 ...
3.5999182701110840e+01
3.5997957229614258e+01

L3 data (lonlat grid) uses 64-bit CMAQ XDR format (no lons,lats,elevations):

SUBSET 9.0 CMAQ
TEMPO_L3
http://tempo.si.edu/,TEMPOSubset
2024-05-13T00:00:00-0000
# data dimensions: timesteps variables layers rows columns:
24 1 1 13 12
# subset indices (0-based time, 1-based layer/row/column): \
 first-timestep last-timestep first-layer last-layer first-row last-row \
 first-column last-column:
0 23 1 1 43 55 139 150
# Variable names:
nitrogen_dioxide_vertical_column_total
# Variable units:
molecules/cm2
# lonlat projection: major_semiaxis minor_semiaxis
6370000 6370000
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:
7750 2950 -167.99 -14.01 0.02 0.02 6 40000 0 1
# IEEE-754 64-bit reals data[variables][timesteps][layers][rows][columns]:

HISTORY: 2019-06-25 plessel.todd@epa.gov
STATUS:  unreviewd tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, printf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <stdlib.h>    /* For malloc(), free(), strtod(), atof(), atoi(). */
#include <math.h>      /* For fabs(). */
#include <ctype.h>     /* For isdigit(). */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For LONGITUDE, Bounds. */
#include "ReadData.h"  /* For openFile(), swathInDomain(), readFileData(). */

/*================================= MACROS =================================*/

/* Name of temporary file created in -tmpdir will have PID appended: */

#define TEMP_FILE_NAME "junk_TEMPOSubset"

/*
 * For L3 grid data, there can be more than one file for a given hour,
 * e.g., T12, so the output will be the mean (of non-missing values) per cell:
 * TEMPO_NO2_L3_V02_20240513T000143Z_S017.nc
 * TEMPO_NO2_L3_V02_20240513T004148Z_S018.nc
 * TEMPO_NO2_L3_V03_20240513T104103Z_S001.nc
 * TEMPO_NO2_L3_V03_20240513T112108Z_S002.nc
 * TEMPO_NO2_L3_V03_20240513T120113Z_S003.nc
 * TEMPO_NO2_L3_V03_20240513T124118Z_S004.nc
 * TEMPO_NO2_L3_V03_20240513T132123Z_S005.nc
 * TEMPO_NO2_L3_V03_20240513T140128Z_S006.nc
 * TEMPO_NO2_L3_V03_20240513T150128Z_S007.nc
 * TEMPO_NO2_L3_V03_20240513T160128Z_S008.nc
 * TEMPO_NO2_L3_V03_20240513T170128Z_S009.nc
 * TEMPO_NO2_L3_V03_20240513T180128Z_S010.nc
 * TEMPO_NO2_L3_V03_20240513T190128Z_S011.nc
 * TEMPO_NO2_L3_V03_20240513T200128Z_S012.nc
 * TEMPO_NO2_L3_V03_20240513T210128Z_S013.nc
 * TEMPO_NO2_L3_V03_20240513T220128Z_S014.nc
 * TEMPO_NO2_L3_V03_20240513T224133Z_S015.nc
 * TEMPO_NO2_L3_V03_20240513T232139Z_S016.nc
 */

/*================================== TYPES ==================================*/

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;       /* File containing list of TEMPO files to read*/
  const char* tmpdir;         /* Name of directory to write temp files. */
  const char* description;    /* User-supplied description. */
  const char* variable;       /* Name of variable to read. */
  Bounds      domain; /* Subset domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  int         yyyymmddhh;     /* First timestamp of subset. */
  int         hours;          /* Number of hours in subset. */
  int         minimumQuality; /* Minimum quality filter [0, 2]. */
  double      maximumCloudFraction; /* Maximum acceptable cloud fraction. */
  double      maximumSolarZenithAngle; /* Max acceptable solar zenith angle. */
  int         allowNegativeCounts;     /* Allow negative molecules/cm2? */
  int         corners;        /* Compute interpolated lon-lat corner points?*/
} Arguments;

/*
 * GridInfo is for L3 gridded data. It describes the full grid of the file.
 * It is initialized from the first file and subsequent files that don't match
 * are skipped.
 */

typedef struct {
  size_t columns;
  size_t rows;
  double west;
  double south;
  double cellWidth;
  double cellHeight;
} GridInfo;

/* Data type: */

typedef struct {
  Arguments   arguments;    /* User-supplied (command-line) arguments. */
  const char* product;      /* E.g., "NO2_L2" or "L3". */
  const char* variable;     /* E.g., "no2". */
  char        units[ 80 ];  /* Units of variable. */
  FileName    tempFileName; /* Name of temp file of output subset data. */
  FILE*       tempFile;     /* Temp file of output subset data. */
  long long*  yyyydddhhmm;  /* Timestamp per output subset scan. */
  long long*  points;       /* Number of points per output subset scan. */
  int         isL3;         /* Is this a set of L3 grid files? */
  GridInfo    gridInfo;     /* Of L3 gridded data. */
  size_t      gridSubsetIndices[2][2]; /* [COLUMN, ROW][FIRST, LAST]. 1-based*/
  double*     gridSubsetValues; /* For L3: subset variable read buffer. */
  double*     gridSubsetScratch; /* For L3: subset scratch read buffer. */
  double*     gridSubsetMeans; /* For L3: subset hourly mean output buffer. */
  unsigned char* gridSubsetCounts; /* For L3 hourly counts per cell for mean.*/
  int         gridSubsetMeans_yyyymmddhh; /* timestamp of gridSubsetMeans. */
  int         scans;        /* Number of output subset scans. */
  int         ok;           /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static int isValidArguments( const Arguments* const arguments ) {
  const int result =
    arguments != 0 &&
    arguments->listFile &&
    arguments->listFile[ 0 ] &&
    arguments->description &&
    arguments->description[ 0 ] &&
    arguments->variable &&
    arguments->variable[ 0 ] &&
    isValidBounds( (const double (*)[2]) arguments->domain ) &&
    isValidYYYYMMDDHH( arguments->yyyymmddhh ) &&
    arguments->hours > 0 &&
    arguments->minimumQuality >= 0 &&
    arguments->minimumQuality <= 2 &&
    arguments->maximumCloudFraction >= 0.0 &&
    arguments->maximumCloudFraction <= 1.0 &&
    arguments->maximumSolarZenithAngle >= 0.0 &&
    arguments->maximumSolarZenithAngle <= 90.0 &&
    ( arguments->allowNegativeCounts == 0 ||
      arguments->allowNegativeCounts == 1 ) &&
    ( arguments->corners == 0 || arguments->corners == 1 );
  return result;
}

static void deallocate( Data* data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static void readData( Data* const data );

static size_t initializeGridInfo( const int file,
                                  const size_t rows,
                                  const size_t columns,
                                  Data* const data );

static int matchesGridInfo( const int file,
                            const size_t rows,
                            const size_t columns,
                            Data* const data );

static void subsetCellIndex( const double coordinate, const double minimum,
                             const double cellSize, const size_t count,
                             size_t* index );

static char* readListFileAndAllocateScanArrays( Data* const data );

static void allocateScanArrays( const char* const listFileContent,
                                Data* const data );

static long long swathFileTimestamp( const char* const fileName );

static const char* swathFileProduct( const char* const fileName );

static const char* swathFileVariable( const char* const fileName );

static int readFileInfo( const char* const fileName,
                         const char** product,
                         const char** variable,
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
                                      double* const values,
                                      double* const scratch );

static void writeGridHeader( const Data* const data );

static void writeSubsetGridData( Data* const data,
                                 const int yyyymmddhh );

static void writeDataBeforeYYYYMMDDHH( Data* const data, const int yyyymmddhh);

static void updateMeans( Data* const data, const int reinitialize );

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
PURPOSE: main - Extract a subset of data from a list of TEMPO files and write
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
    readData( &data ); /* Read subset of TEMPO files and write temp files. */

    if ( data.isL3 ) {
      ok = data.ok;
    } else if ( data.ok && data.scans ) {
      streamData( &data ); /* Write header and temp file to stdout & rm temp.*/
      ok = data.ok;
    }
  }

  deallocate( &data );
  DEBUG( fprintf( stderr, "%s exiting main with value %d\n\n", argv[0], ! ok);)
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
    free( data->points ), data->points = 0;
  }

  if ( data->yyyydddhhmm ) {
    free( data->yyyydddhhmm ), data->yyyydddhhmm = 0;
  }

  if ( data->gridSubsetValues ) {
    free( data->gridSubsetValues ), data->gridSubsetValues = 0;
  }

  if ( data->gridSubsetScratch ) {
    free( data->gridSubsetScratch ), data->gridSubsetScratch = 0;
  }

  if ( data->gridSubsetMeans ) {
    free( data->gridSubsetMeans ), data->gridSubsetMeans = 0;
  }

  if ( data->gridSubsetCounts ) {
    free( data->gridSubsetCounts ), data->gridSubsetCounts = 0;
  }

  if ( data->tempFile ) {
    fclose( data->tempFile ), data->tempFile = 0;
  }

  if ( data->tempFileName[ 0 ] ) {
    unlink( data->tempFileName );
    memset( data->tempFileName, 0, sizeof data->tempFileName );
  }

  memset( data, 0, sizeof (Data) );
}



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* name  Name of program.
******************************************************************************/

static void printUsage( const char* name ) {
  assert( name ); assert( *name );
  fprintf( stderr,
           "\a\n\n%s - Extract a lon-lat subset of data from a list of\n"
          "TEMPO NetCDF4 files and write it to stdout as XDR binary format.\n",
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
  fprintf( stderr, "  [-minimumQuality value]\\\n" );
  fprintf( stderr, "  [-maximumCloudFraction value]\\\n" );
  fprintf( stderr, "  [-maximumSolarZenithAngle value]\\\n" );
  fprintf( stderr, "  [-allowNegativeCounts]\\\n" );
  fprintf( stderr, "  [-corners]\n\n" );
  fprintf( stderr, "Note:\ntimestamp is in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory were a transient file is " );
  fprintf( stderr, "written.\nIt should have enough disk space (1TB).\n" );
  fprintf( stderr, "-minimumQuality option filters-out values less than the " );
  fprintf( stderr, "specified allowed minimum to accept:\n" );
  fprintf( stderr, "[0 = only 'normal' values allowed,\n" );
  fprintf( stderr, "1 = either 'normal' or 'suspect' values are allowed,\n" );
  fprintf( stderr, "2 = no quality filtering is applied].\n" );
  fprintf( stderr, "Default is 0.\n" );
  fprintf( stderr, "-maximumCloudFraction option filter-out values greater " );
  fprintf( stderr, "than the specified value [0.0, 1.0]. Default is 1.0.\n");
  fprintf( stderr, "-maximumSolarZenithAngle option filter-out values greater " );
  fprintf( stderr, "than the specified value [0.0, 90.0]. Default is 90.0.\n");
  fprintf( stderr, "-allowNegativeCounts will allow negative counts of "
                   "molecules/cm2 (non-physical).\n" );
  fprintf( stderr, "-corners option will output 8 additional variables:\n" );
  fprintf( stderr, "  Longitude_SW Longitude_SE Longitude_NW Longitude_NE\n");
  fprintf( stderr, "  Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "that are the linearly interpolated " );
  fprintf( stderr, "(and edge extrapolated)\n" );
  fprintf( stderr, "corner points for each center-pixel point.\n\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-files vnpaerdt_files \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc "
                   "\"http://tempo.si.edu/,TEMPOSubset\" \\\n");
  fprintf( stderr, "-timestamp 2023101700 -hours 24 \\\n" );
  fprintf( stderr, "-variable vertical_column_total \\\n" );
  fprintf( stderr, "-domain -126 25 -65 50 -corners > subset.xdr\n\n" );
  fprintf( stderr, "AOD over US on Novembr 28, 2017.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Swath 2.0\n" );
  fprintf( stderr, "http://tempo.si.edu/,TEMPOSubset\n");
  fprintf( stderr, "2023-10-17T00:00:00-0000\n" );
  fprintf( stderr, "# Dimensions: variables timesteps scans:\n" );
  fprintf( stderr, "11 24 2\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Longitude Latitude no2_vertical_column_total "
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
  fprintf( stderr, "20232831659\n" );
  fprintf( stderr, "20232831705\n" );
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
  arguments->minimumQuality = 0; /* Range: [0, 2], default 0 = normal. */
  arguments->maximumCloudFraction = 1.0; /* Allow any fraction of clouds. */
  arguments->maximumSolarZenithAngle = 90.0; /* Allow any solar zenith angle */

  result = argc >= 18 && argc <= 26;

  for ( arg = 1; result && arg < argc; ++arg ) {

    if ( ! strcmp( argv[ arg ], "-files" ) && arg + 1 < argc ) {
      ++arg;
      arguments->listFile = argv[ arg ];
    } else if ( ! strcmp( argv[ arg ], "-tmpdir" ) && arg + 1 < argc ) {
      ++arg;
      arguments->tmpdir = argv[ arg ];
    } else if ( ! strcmp( argv[ arg ], "-desc" ) && arg + 1 < argc ) {
      ++arg;
      arguments->description = argv[ arg ];
    } else if ( ! strcmp( argv[ arg ], "-timestamp" ) && arg + 1 < argc ) {
      ++arg;
      arguments->yyyymmddhh = atoi( argv[ arg ] );
      result = isValidYYYYMMDDHH( arguments->yyyymmddhh );
    } else if ( ! strcmp( argv[ arg ], "-hours" ) && arg + 1 < argc ) {
      ++arg;
      arguments->hours = atoi( argv[ arg ] );
      result = arguments->hours > 0;
    } else if ( ! strcmp( argv[ arg ], "-variable" ) && arg + 1 < argc ) {
      ++arg;
      arguments->variable = argv[ arg ];
      result = arguments->variable[ 0 ] != '\0' &&
               arguments->variable[ 0 ] != '-';
    } else if ( ! strcmp( argv[ arg ], "-domain" ) && arg + 4 < argc ) {
      double value = 0.0;
      char* end = 0;
      ++arg;
      value = strtod( argv[ arg ], &end );
      result += end != argv[ arg ];
      arguments->domain[ LONGITUDE ][ MINIMUM ] = value;
      ++arg;
      value = strtod( argv[ arg ], &end );
      result += end != argv[ arg ];
      arguments->domain[ LATITUDE ][ MINIMUM ] = value;
      ++arg;
      value = strtod( argv[ arg ], &end );
      result += end != argv[ arg ];
      arguments->domain[ LONGITUDE ][ MAXIMUM ] = value;
      ++arg;
      value = strtod( argv[ arg ], &end );
      result += end != argv[ arg ];
      arguments->domain[ LATITUDE ][ MAXIMUM ] = value;
      result = result == 5 &&
               isValidBounds( (const double (*)[2]) arguments->domain );
    } else if ( ! strcmp( argv[ arg ], "-minimumQuality") &&
                arg + 1 < argc && isdigit( argv[ arg + 1 ][ 0 ] ) ) {
      ++arg;
      arguments->minimumQuality = atoi( argv[ arg ] );
      result =
        arguments->minimumQuality >= 0 && arguments->minimumQuality <= 2;
    } else if ( ! strcmp( argv[ arg ], "-maximumCloudFraction" ) &&
                arg + 1 < argc && isdigit( argv[ arg + 1 ][ 0 ] ) ) {
      ++arg;
      arguments->maximumCloudFraction = atof( argv[ arg ] );
      result =
        arguments->maximumCloudFraction >= 0.0 &&
        arguments->maximumCloudFraction <= 1.0;
    } else if ( ! strcmp( argv[ arg ], "-maximumSolarZenithAngle" ) &&
                arg + 1 < argc && isdigit( argv[ arg + 1 ][ 0 ] ) ) {
      ++arg;
      arguments->maximumSolarZenithAngle = atof( argv[ arg ] );
      result =
        arguments->maximumSolarZenithAngle >= 0.0 &&
        arguments->maximumSolarZenithAngle <= 90.0;
    } else if ( ! strcmp( argv[ arg ], "-allowNegativeCounts" ) ) {
      arguments->allowNegativeCounts = 1;
    } else if ( ! strcmp( argv[ arg ], "-corners" ) ) {
      arguments->corners = 1;
    } else {
      result = 0;
    }
  }

  result = result && isValidArguments( arguments );

  if ( ! result ) {
    fprintf( stderr, "\nInvalid/insufficient command-line arguments.\n" );
  }

  return result;
}



/******************************************************************************
PURPOSE: readData - Read swath data from each listed TEMPO file and
         write the lon-lat subset of data to the temporary file.
INPUTS:  Data* data  Data description to read.
OUTPUTS: Data* data  data->yyyydddhhmm, points, cellWidthHeight, scans, ok,
                     tempFile = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  char* listFileContent = readListFileAndAllocateScanArrays( data );
  int wroteSomeData = 0;
  data->ok = listFileContent != 0;

  if ( data->ok ) {
    char* fileName       = 0;
    char* end            = 0;
    size_t rows          = 0;
    size_t columns       = 0;
    size_t size          = 0;
    double* buffer       = 0;
    double* longitudes   = 0;
    double* latitudes    = 0;
    double* values       = 0;
    double* scratch      = 0;
    double* longitudesSW = 0;
    double* longitudesSE = 0;
    double* longitudesNW = 0;
    double* longitudesNE = 0;
    double* latitudesSW  = 0;
    double* latitudesSE  = 0;
    double* latitudesNW  = 0;
    double* latitudesNE  = 0;
    unsigned char* mask  = 0;
    int wroteGridHeader = 0;

    /* If corners option is used then process gridded L3 files as swaths. */

    data->isL3 =
      ! arguments->corners &&
      strstr( listFileContent, "_L3_V" ) != 0 &&
      strstr( listFileContent, "_L2_V" ) == 0 &&
      strstr( listFileContent, "_PM25_L3_" ) == 0;

    /* Get each line of list file. It is the TEMPO data file to read: */

    for ( fileName = strtok_r( listFileContent, "\n", &end );
          fileName;
          fileName = strtok_r( 0, "\n", &end ) ) {
      int file = 0;
      long long yyyymmddhhmm = 0;
      int changedDimensions = 0;

      if ( strstr( fileName, "_PM25_L3_V" ) ||
           strstr( fileName, "_ADP_L2_V" ) ) {
        data->variable = arguments->variable;
      }

      data->ok =
        readFileInfo( fileName,
                      &(data->product),
                      &(data->variable),
                      &file, &yyyymmddhhmm, &rows, &columns, &size,
                      &changedDimensions );

      if ( data->ok && ! strcmp( data->product, "PM25_L3" ) ) {
        data->variable = arguments->variable;
      }

      DEBUG( fprintf( stderr,
                     "\n%s %lld %s %s isL3 = %d, %lu x %lu = %lu "
                     "changed = %d, ok = %d\n",
                     fileName, yyyymmddhhmm, data->product, data->variable,
                     data->isL3, rows, columns, size, changedDimensions,
                     data->ok ); )

      if ( data->ok ) {

        if ( data->isL3 ) {

          if ( data->gridInfo.rows == 0 ) {
            data->ok = initializeGridInfo( file, rows, columns, data );
          } else {
            data->ok =
              rows == data->gridInfo.rows && columns == data->gridInfo.columns;

            if ( data->ok ) {
              data->ok = matchesGridInfo( file, rows, columns, data );
            }

            if ( ! data->ok ) {
              fprintf( stderr, "Skipping file with unmatched grid '%s'.\n",
                       fileName );
            }
          }
        } else if ( changedDimensions ) {
          const int corners = arguments->corners;
          const size_t variables = 4 + 8 * corners; /* lon,lat,values,scratch*/
          const size_t dataSize = variables * size * sizeof (double);
          const size_t maskSize = size * sizeof (unsigned char);
          const size_t bytes = dataSize + maskSize;

          if ( buffer ) {
            free( buffer ), buffer = 0;
          }

          buffer = malloc( bytes );
          data->ok = buffer != 0;

          if ( buffer ) {
            memset( buffer, 0, bytes );
            longitudes = buffer;
            latitudes  = longitudes + size;
            values     = latitudes  + size;
            scratch    = values + size;

            if ( corners ) {
              longitudesSW = scratch      + size;
              longitudesSE = longitudesSW + size;
              longitudesNW = longitudesSE + size;
              longitudesNE = longitudesNW + size;
              latitudesSW  = longitudesNE + size;
              latitudesSE  = latitudesSW  + size;
              latitudesNW  = latitudesSE  + size;
              latitudesNE  = latitudesNW  + size;
              mask = (unsigned char*) ( latitudesNE + size );
            } else {
              mask = (unsigned char*) ( scratch + size );
            }
          } else {
            fprintf( stderr,
                     "\nCan't allocate %lu bytes "
                     "to complete the requested action.\n", bytes );
          }
        }

        if ( data->ok ) {
          readCoordinatesAndValues( data, file, rows, columns,
                                    longitudes, latitudes, values, scratch );
        }

        closeFile( file );
        file = -1;

        if ( data->ok ) {

          if ( data->isL3 ) {

            if ( ! wroteGridHeader ) {
              writeGridHeader( data );
              wroteGridHeader = 1;
            }

            writeSubsetGridData( data, yyyymmddhhmm / 100 );
          } else {
            const size_t subsetPoints =
              pointsInDomain( (const double (*)[2]) data->arguments.domain,
                              size, longitudes, latitudes, values, mask );
            DEBUG( fprintf( stderr, "subsetPoints = %lu\n", subsetPoints ); )

            if ( subsetPoints ) {

              if ( longitudesSW ) {
                computeCorners( rows, columns, longitudes, latitudes,
                                longitudesSW, longitudesSE,
                                longitudesNW, longitudesNE,
                                latitudesSW, latitudesSE,
                                latitudesNW, latitudesNE );
              }

              writeSubset( data, yyyymmddhhmm, subsetPoints, size,
                           mask, longitudes, latitudes, values,
                           longitudesSW, longitudesSE,
                           longitudesNW, longitudesNE,
                           latitudesSW, latitudesSE,
                           latitudesNW, latitudesNE );
            }
          }

          if ( data->ok ) {
            wroteSomeData = 1;
          }
        }
      } /* If ok */
    } /* End loop on listFile. */

    free( listFileContent ), listFileContent = 0;

    if ( buffer ) {
      free( buffer ), buffer = 0;
    }

    if ( data->tempFile ) { /* Done writing to temp file so close it: */
      fclose( data->tempFile ), data->tempFile = 0;
    }
  }

  data->ok = wroteSomeData;

  DEBUG(fprintf(stderr, "\nEnd of file processing, data->ok = %d\n",data->ok);)

  if ( data->ok && data->isL3 ) {
    const int end_yyyymmddhh =
      incrementHours( arguments->yyyymmddhh, arguments->hours );
    writeDataBeforeYYYYMMDDHH( data, end_yyyymmddhh );
  }
}



/******************************************************************************
PURPOSE: initializeGridInfo - Initialize grid info and buffers.
INPUTS:  const int file           File to read coordinates from.
         const size_t rows        File grid rows.
         const size_t columns     File grid columns.
         Data* const data         data->arguments->domain.
OUTPUTS: Data* const data         data->gridInfo
                                  data->gridSubsetIndices
                                  data->gridSubsetValues
                                  data->gridSubsetScratch
                                  data->gridSubsetMeans
                                  data->gridSubsetCounts.
RETURNS: size_t number of subset points in data->gridSubsetIndices domain.
******************************************************************************/

static size_t initializeGridInfo( const int file,
                                  const size_t rows,
                                  const size_t columns,
                                  Data* const data ) {

  size_t result = 0;
  size_t count = rows > columns ? rows : columns;
  size_t bytes = count * sizeof (double);
  double* buffer = 0;

  assert( rows ); assert( columns ); assert( data ); assert( data->ok );
  assert( data->isL3 ); assert( data->gridInfo.rows == 0 );
  assert( data->gridSubsetValues == 0 );
  assert( data->gridSubsetScratch == 0 );
  assert( data->gridSubsetMeans == 0 );
  assert( data->gridSubsetCounts == 0 );

  buffer = malloc( bytes );

  if ( ! buffer ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes to complete the requested action.\n",
             bytes );
  } else {
    char unused[ 80 ] = "";
    size_t indices[ 2 ][ 2 ] = { { 1, 1 }, { 1, 1 } };
    indices[ COLUMN ][ LAST ] = columns;
    memset( buffer, 0, bytes );
    data->ok =
      readFileData( file, data->product, "longitude", 1, columns,
                    (const size_t (*)[2]) indices,
                    0, 1.0, 90.0, 0, unused, buffer, 0 ) > 0;

    if ( data->ok ) {
      const double west = buffer[ 0 ];
      const double east = buffer[ columns - 1 ];
      const double cellWidth = ( east - west ) / ( columns - 1 );
      data->ok =
        cellWidth > 0.0 &&
        west >= -180.0 && west + ( columns - 1 ) * cellWidth <= 180.0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nFailed to read valid longitudes.\n" );
      } else {
        const Arguments* const arguments = &data->arguments;
        const double longitudeMinimum =
          arguments->domain[ LONGITUDE ][ MINIMUM ];
        const double longitudeMaximum =
          arguments->domain[ LONGITUDE ][ MAXIMUM ];
        int outside = longitudeMaximum < west || longitudeMinimum > east;

        if ( ! outside ) {
          indices[ COLUMN ][ LAST ] = 1;
          indices[ ROW ][ LAST ] = rows;
          data->ok =
            readFileData( file, data->product, "latitude", rows, 1,
                          (const size_t (*)[2]) indices,
                          0, 1.0, 90.0, 0, unused, buffer, 0 ) > 0;

          if ( data->ok ) {
            const double south = buffer[ 0 ];
            const double north = buffer[ rows - 1 ];
            const double cellHeight = ( north - south ) / ( rows - 1 );
            data->ok =
              cellHeight > 0.0 &&
              south >= -90.0 && south + ( rows - 1 ) * cellHeight <= 90.0;

            if ( ! data->ok ) {
              fprintf( stderr, "\nFailed to read valid latitudes.\n" );
            } else {
              const double latitudeMinimum =
                arguments->domain[ LATITUDE ][ MINIMUM ];
              const double latitudeMaximum =
                arguments->domain[ LATITUDE ][ MAXIMUM ];
              outside = latitudeMaximum < south || latitudeMinimum > north;

              if ( ! outside ) {
                size_t firstColumn = 1;
                size_t lastColumn = columns;
                size_t firstRow = 1;
                size_t lastRow = rows;

                subsetCellIndex( longitudeMinimum, west, cellWidth, columns,
                                 &firstColumn );
                subsetCellIndex( longitudeMaximum, west, cellWidth, columns,
                                 &lastColumn );
                subsetCellIndex( latitudeMinimum, south, cellHeight, rows,
                                 &firstRow );
                subsetCellIndex( latitudeMaximum, south, cellHeight, rows,
                                 &lastRow );

                assert( firstColumn >= 1 );
                assert( firstColumn <= columns );
                assert( lastColumn >= firstColumn );
                assert( lastColumn <= columns );

                assert( firstRow >= 1 );
                assert( firstRow <= rows );
                assert( lastRow >= firstRow );
                assert( lastRow <= rows );

                data->gridSubsetIndices[ COLUMN ][ FIRST ] = firstColumn;
                data->gridSubsetIndices[ COLUMN ][ LAST  ] = lastColumn;
                data->gridSubsetIndices[ ROW    ][ FIRST ] = firstRow;
                data->gridSubsetIndices[ ROW    ][ LAST  ] = lastRow;

                data->gridInfo.rows       = rows;
                data->gridInfo.columns    = columns;
                data->gridInfo.cellWidth  = cellWidth;
                data->gridInfo.cellHeight = cellHeight;
                data->gridInfo.west       = west;
                data->gridInfo.south      = south;

                result =
                  ( 1 + lastRow - firstRow ) * ( 1 + lastColumn - firstColumn);

                DEBUG( fprintf( stderr,
                                "Grid subset: columns = [%lu, %lu] of %lu, "
                                "rows [%lu, %lu] of %lu, "
                                "origin (%f, %f), cell size %f x %f\n",
                                data->gridSubsetIndices[ COLUMN ][ FIRST ],
                                data->gridSubsetIndices[ COLUMN ][ LAST  ],
                                data->gridInfo.columns,
                                data->gridSubsetIndices[ ROW    ][ FIRST ],
                                data->gridSubsetIndices[ ROW    ][ LAST  ],
                                data->gridInfo.rows,
                                data->gridInfo.west,
                                data->gridInfo.south,
                                data->gridInfo.cellWidth,
                                data->gridInfo.cellHeight ); )
              }
            }
          }
        }
      }
    }
  }

  assert( buffer );
  free( buffer ), buffer = 0;

  /* Allocate buffers needed to read and process L3 subset grid data: */

  if ( result ) {
    data->ok = 0;
    bytes = result * sizeof (double);
    data->gridSubsetValues = malloc( bytes );

    if ( ! data->gridSubsetValues ) {
      fprintf( stderr,
               "\nCan't allocate %lu bytes to complete the requested action.\n",
               bytes );
    } else {
      data->gridSubsetScratch = malloc( bytes );

      if ( ! data->gridSubsetScratch ) {
        fprintf( stderr,
                 "\nCan't allocate %lu bytes "
                 "to complete the requested action.\n",
                 bytes );
      } else {
        data->gridSubsetMeans = malloc( bytes );

        if ( ! data->gridSubsetMeans ) {
          fprintf( stderr,
                   "\nCan't allocate %lu bytes "
                   "to complete the requested action.\n",
                   bytes );
        } else {
          bytes = result * sizeof (unsigned char);
          data->gridSubsetCounts = malloc( bytes );

          if ( ! data->gridSubsetCounts ) {
            fprintf( stderr,
                     "\nCan't allocate %lu bytes "
                     "to complete the requested action.\n",
                     bytes );
          } else {
            data->ok = 1;
          }
        }
      }
    }
  }

  if ( ! data->ok ) {
    data->gridInfo.columns = data->gridInfo.rows = 0;

    if ( data->gridSubsetValues ) {
      free( data->gridSubsetValues ), data->gridSubsetValues = 0;
    }

    if ( data->gridSubsetScratch ) {
      free( data->gridSubsetScratch ), data->gridSubsetScratch = 0;
    }

    if ( data->gridSubsetMeans ) {
      free( data->gridSubsetMeans ), data->gridSubsetMeans = 0;
    }

    if ( data->gridSubsetCounts ) {
      free( data->gridSubsetCounts ), data->gridSubsetCounts = 0;
    }

    result = 0;
  }

  data->gridSubsetMeans_yyyymmddhh = 0;

  assert( data->ok == 0 || (
          data->gridInfo.columns == columns &&
          data->gridInfo.rows == rows &&
          data->gridInfo.west >= -180.0 &&
          data->gridInfo.west <=  180.0 &&
          data->gridInfo.south >= -90.0 &&
          data->gridInfo.south <=  90.0 &&
          data->gridInfo.cellWidth > 0.0 &&
          data->gridInfo.cellHeight > 0.0 &&
          data->gridSubsetIndices[ COLUMN ][ FIRST ] >= 1 &&
          data->gridSubsetIndices[ COLUMN ][ LAST  ] >=
          data->gridSubsetIndices[ COLUMN ][ FIRST ] &&
          data->gridSubsetIndices[ COLUMN ][ LAST  ] <= columns &&
          data->gridSubsetIndices[ ROW    ][ FIRST ] >= 1 &&
          data->gridSubsetIndices[ ROW    ][ LAST  ] >=
          data->gridSubsetIndices[ ROW    ][ FIRST ] &&
          data->gridSubsetIndices[ ROW    ][ LAST  ] <= rows &&
          data->gridInfo.west  + ( columns - 1 ) * data->gridInfo.cellWidth
            <= 180.0 &&
          data->gridInfo.south + ( rows    - 1 ) * data->gridInfo.cellHeight
            <=  90.0 &&
          data->gridSubsetValues  != 0 &&
          data->gridSubsetScratch != 0 &&
          data->gridSubsetMeans   != 0 &&
          data->gridSubsetCounts  != 0 &&
          data->gridSubsetMeans_yyyymmddhh == 0 &&
          result == ( 1 + data->gridSubsetIndices[ ROW    ][ LAST  ] -
                          data->gridSubsetIndices[ ROW    ][ FIRST  ] ) *
                    ( 1 + data->gridSubsetIndices[ COLUMN ][ LAST  ] -
                          data->gridSubsetIndices[ COLUMN ][ FIRST  ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: matchesGridInfo - Does file grid match data->gridInfo?
INPUTS:  const int file        File to read coordinates from.
         const size_t rows     File grid rows.
         const size_t columns  File grid columns.
         Data* const data      data->gridInfo
OUTPUTS: Data* const data      data->ok
RETURNS: int 1 if matched, else 0.
******************************************************************************/

static int matchesGridInfo( const int file,
                            const size_t rows,
                            const size_t columns,
                            Data* const data ) {

  const double tolerance = 1e6;
  size_t result = 0;
  size_t count = rows > columns ? rows : columns;
  size_t bytes = count * sizeof (double);
  double* buffer = 0;

  assert( rows ); assert( columns ); assert( data ); assert( data->ok );
  assert( data->isL3 ); assert( data->gridInfo.rows );

  data->ok = 0;
  buffer = malloc( bytes );

  if ( ! buffer ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes to complete the requested action.\n",
             bytes );
  } else {
    char unused[ 80 ] = "";
    size_t indices[ 2 ][ 2 ] = { { 1, 1 }, { 1, 1 } };
    indices[ COLUMN ][ LAST ] = columns;
    memset( buffer, 0, bytes );
    data->ok =
      readFileData( file, data->product, "longitude", 1, columns,
                    (const size_t (*)[2]) indices,
                    0, 1.0, 90.0, 0, unused, buffer, 0 ) > 0;

    if ( data->ok ) {
      const double west = buffer[ 0 ];
      const double east = buffer[ columns - 1 ];
      const double cellWidth = ( east - west ) / ( columns - 1 );
      data->ok =
        cellWidth > 0.0 &&
        west >= -180.0 && west + ( columns - 1 ) * cellWidth <= 180.0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nFailed to read valid longitudes.\n" );
      } else {
        data->ok =
          fabs( west - data->gridInfo.west ) < tolerance &&
          fabs( cellWidth - data->gridInfo.cellWidth ) < tolerance;

        if ( ! data->ok ) {
          fprintf( stderr, "\nRead unmatched longitudes.\n" );
        } else {
          indices[ COLUMN ][ LAST ] = 1;
          indices[ ROW ][ LAST ] = rows;
          data->ok =
            readFileData( file, data->product, "latitude", rows, 1,
                          (const size_t (*)[2]) indices,
                          0, 1.0, 90.0, 0, unused, buffer, 0 ) > 0;

          if ( data->ok ) {
            const double south = buffer[ 0 ];
            const double north = buffer[ rows - 1 ];
            const double cellHeight = ( north - south ) / ( rows - 1 );
            data->ok =
              cellHeight > 0.0 &&
              south >= -90.0 && south + ( rows - 1 ) * cellHeight <= 90.0;

            if ( ! data->ok ) {
              fprintf( stderr, "\nFailed to read valid latitudes.\n" );
            } else {
              data->ok =
                fabs( south - data->gridInfo.south ) < tolerance &&
                fabs( cellHeight - data->gridInfo.cellHeight ) < tolerance;

              if ( ! data->ok ) {
                fprintf( stderr, "\nRead unmatched latitudes.\n" );
              } else {
                result = 1;
              }
            }
          }
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: subsetCellIndex - Subset 1-based cell index to grid dimension.
INPUTS:  const double coordinate  Coordinate to compute index of.
         const double minimum     Minimum cell dimension coordinate.
         const double cellSize    Grid cell length in dimension.
         const size_t count       Number of grid cells along dimension.
         size_t* index            Cell index initialized to [1, count].
OUTPUTS: size_t* index            Updated 1-based cell index of coordinate.
******************************************************************************/


static void subsetCellIndex( const double coordinate, const double minimum,
                             const double cellSize, const size_t count,
                             size_t* index ) {

  assert( cellSize > 0.0 ); assert( count );
  assert( index ); assert( *index >= 1 ); assert( *index <= count );

  {
    const double delta = coordinate - minimum;

    if ( delta >= 0.0 ) {
      *index = (int) ( delta / cellSize ) + 1;

      if ( *index > count ) {
        *index = count;
      }
    }
  }

  assert( *index >= 1 ); assert( *index <= count );
}



/******************************************************************************
PURPOSE: readListFileAndAllocateScanArrays - Read list file and
         return its contents as a string and allocate timestamps, points and
         cellWidthHeight arrays with length equal to lines in the list file.
INPUTS:  Data* const data  data->arguments.listFile.
OUTPUTS: Data* const data  data->ok, yyyydddhhmm, points, cellWidthHeight
                           allocated.
RETURNS: char* Allocated string containing lines of list file or 0 if failed.
******************************************************************************/

static char* readListFileAndAllocateScanArrays( Data* const data ) {
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
      allocateScanArrays( result, data );

      if ( ! data->ok ) {
        free( result );
        result = 0;
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: allocateScanArrays - Allocate per-scan arrays for timestamps, points.
INPUTS:  const char* const listFileContent  String content of list file.
         Data* const data                   Data to allocate.
OUTPUTS: Data* const data  data->ok, data->yyyydddhhmm, points, cellWidthHeight
                           allocated.
******************************************************************************/

static void allocateScanArrays( const char* const listFileContent,
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

      if ( data->ok ) {
        memset( data->yyyydddhhmm, 0, bytes );
        memset( data->points, 0, bytes );
      } else {
        fprintf( stderr, "\nCan't allocate %lu bytes to complete the "
                 "requested action.\n", bytes );
      }
    }
  }
}



/******************************************************************************
PURPOSE: swathFileTimestamp - Timestamp of swath file.
INPUTS:  const char* const fileName   Name of TEMPO file.
RETURNS: long long yyyymmddhhmm of file or 0 if failed and message on stderr.
NOTES: File names look like:
S5P_OFFL_L2__NO2____20171128T163259_20171128T181628_00657_03_001108_20171220T145115.nc
******************************************************************************/

static long long swathFileTimestamp( const char* const fileName ) {
  long long result = 0;
  const char* slash = strchr( fileName, '/' );
  const char* const name = slash ? slash + 1 : fileName;
  const char* const tag = "_20";
  const char* s = strstr( name, tag );

  if ( s ) {
   int c = 0;
   ++s;

    /* Parse YYYYMMDD: */

    for ( c = 0; *s && *s != '_' && c < 8; ++c, ++s ) {
      result = result * 10 + *s - '0';
    }

    /* Parse HHMM: */

    if ( *s == 'T' ) {
      ++s;

      for ( c = 0; *s && *s != '_' && c < 4; ++c, ++s ) {
        result = result * 10 + *s - '0';
      }
    }
  }

  if ( ! isValidYYYYMMDDHHMM( result ) ) {
    fprintf( stderr, "\nInvalid file name timestamp '%s'.\n", fileName );
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: swathFileProduct - Product name of swath file.
INPUTS:  const char* const fileName   Name of TEMPO file.
RETURNS: const char* product of file or 0 if failed and message on stderr.
NOTES: File names look like:
       TEMPO_HCHO_L2_V01_20130715T165956Z_S002G01.nc
       so result is "HCHO_L2" or 0 if failed.
******************************************************************************/

static const char* swathFileProduct( const char* const fileName ) {
  const char* result = 0;
  const char* slash = strchr( fileName, '/' );
  const char* const name = slash ? slash + 1 : fileName;

  if ( strstr( name, "NO2_L2_" ) ) {
    result = "NO2_L2";
  } else if ( strstr( name, "HCHO_L2_" ) ) {
    result = "HCHO_L2";
  } else if ( strstr( name, "O3TOT_L2_" ) ) {
    result = "O3TOT_L2";
  } else if ( strstr( name, "CLDO4_L2_" ) ) {
    result = "CLDO4_L2";
  } else if ( strstr( name, "NO2_L3_" ) ) {
    result = "NO2_L3";
  } else if ( strstr( name, "HCHO_L3_" ) ) {
    result = "HCHO_L3";
  } else if ( strstr( name, "O3TOT_L3_" ) ) {
    result = "O3TOT_L3";
  } else if ( strstr( name, "CLDO4_L3_" ) ) {
    result = "CLDO4_L3";
  } else if ( strstr( name, "AODALH_L2_" ) ) {
    result = "AODALH_L2";
  } else if ( strstr( name, "PM25_L3_" ) ) {
    result = "PM25_L3";
  } else if ( strstr( name, "_ADP_L2_V" ) ) {
    result = "ADP_L2";
  } else {
    fprintf( stderr, "\nInvalid/unsupported file name '%s'.\n", fileName );
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: swathFileVariable - Variable name of swath file.
INPUTS:  const char* const fileName   Name of TEMPO file.
RETURNS: const char* variable of file or 0 if failed and message on stderr.
NOTES:   File names look like:
         TEMPO_CLDO4_L2_V01_20231017T111336Z_S001G01.nc
         TEMPO_HCHO_L2_V01_20231017T111336Z_S001G01.nc
         TEMPO_NO2_L2_V01_20231017T111336Z_S001G01.nc
         TEMPO_O3TOT_L2_V01_20231017T111336Z_S001G01.nc
         so result is "hcho", "no2", etc. or 0 if failed.
******************************************************************************/

static const char* swathFileVariable( const char* const fileName ) {
  const char* result = 0;
  const char* slash = strchr( fileName, '/' );
  const char* const name = slash ? slash + 1 : fileName;

  if ( strstr( name, "NO2" ) ) {
    result = "no2";
  } else if ( strstr( name, "HCHO" ) ) {
    result = "hcho";
  } else if ( strstr( name, "O3" ) ) {
    result = "o3";
  } else if ( strstr( name, "CLDO4" ) ) {
    result = "cloud";
  } else if ( strstr( name, "AODALH" ) ) {
    result = "aod";
  } else if ( strstr( name, "PM25" ) ) {
    result = "pm25";
  } else if ( strstr( name, "ADP" ) ) {
    result = "smoke";
  } else {
    fprintf( stderr, "\nInvalid/unsupported file name '%s'.\n", fileName );
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileInfo - Parse file name timestamp, product and variable and
         open it and read dimensions.
INPUTS:  const char* const fileName    Name of file to check bounds/dims.
         size_t* const rows            Rows of data in previous file or 0.
         size_t* const columns         Columns of data in previous file or 0.
         const char** variable         E.g., "pm25_ge".
OUTPUTS: const char** product          E.g., "NO2_L2" or "L3".
         const char** variable         E.g., "no2".
         int* const file               NetCDF file id of file.
         long long* const yyyymmddhhmm Timestamp of file.
         size_t* const rows            Rows of data in file.
         size_t* const columns         Columns of data in file.
         size_t* const size            rows * columns.
         int* const changedDimensions  1 if rows or columns changed.
RETURNS: int 1 if in domain else 0.
******************************************************************************/

static int readFileInfo( const char* const fileName,
                         const char** product,
                         const char** variable,
                         int* const file,
                         long long* const yyyymmddhhmm,
                         size_t* const rows,
                         size_t* const columns,
                         size_t* const size,
                         int* const changedDimensions ) {

  int result = 0;
  assert( fileName ); assert( *fileName ); assert( product ), assert( file );
  assert( yyyymmddhhmm );
  assert( rows ); assert( columns ); assert( size );
  assert( changedDimensions );

  *yyyymmddhhmm = swathFileTimestamp( fileName );
  *product = 0;
  *changedDimensions = 0;

  if ( *yyyymmddhhmm ) {
    *product = swathFileProduct( fileName );

    if ( *product ) {

      if ( strcmp( *product, "PM25_L3" ) && strcmp( *product, "ADP_L2" ) ) {
        *variable = swathFileVariable( fileName );
      }

      if ( *variable ) {
        *file = openFile( fileName );

        if ( *file != -1 ) {
          size_t rows2 = 0;
          size_t columns2 = 0;
          result =
            readFileDimensions( *file, *product, *variable, &rows2, &columns2 );

          if ( result ) {

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

  *size = *rows * *columns;
  result = result && *size;

  if ( ! result ) {

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
OUTPUTS: Data* const data              data->ok, data->units[ 80 ]
         double* const longitudes[ rows * columns ]  Longitudes read.
         double* const latitudes[  rows * columns ]  Latitudes read.
         double* const values[     rows * columns ]  Values read.
         double* const scratch[    rows * columns ]  Temp buffer for reorder.
******************************************************************************/

static void readCoordinatesAndValues( Data* const data,
                                      const int file,
                                      const size_t rows,
                                      const size_t columns,
                                      double* const longitudes,
                                      double* const latitudes,
                                      double* const values,
                                      double* const scratch ) {

  char unused[ 80 ] = "";

  assert( data ); assert( data->ok ); assert( data->arguments.variable );
  assert( file > -1 );
  assert( rows ); assert( columns );

  if ( ! data->isL3 ) {
    const char* longitude = "longitude";
    const char* latitude = "latitude";
    assert( longitudes ); assert( latitudes );
    assert( longitudes != latitudes );
    assert( values ); assert( scratch );
    assert( values != scratch );

    if ( ! strcmp( data->product, "PM25_L3" ) ) {

      if ( strstr( data->variable, "_ge" ) ) {
        longitude = "lon_ge";
        latitude = "lat_ge";
      } else {
        longitude = "lon_gw";
        latitude = "lat_gw";
      }
    }

    data->ok =
      readFileData( file, data->product, longitude, rows, columns,
                    (const size_t (*)[2]) data->gridSubsetIndices,
                    0, 1.0, 90.0, 0, unused,
                    longitudes, scratch ) > 0;

    if ( data->ok ) {
      data->ok =
        readFileData( file, data->product, latitude, rows, columns,
                      (const size_t (*)[2]) data->gridSubsetIndices,
                      0, 1.0, 90.0, 0, unused,
                      latitudes, scratch ) > 0;

      if ( data->ok ) {
        data->ok =
          clampInvalidCoordinates( rows * columns, longitudes, latitudes );
      }
    }
  }

  if ( data->ok ) {
    const size_t readRows =
      data->isL3 ?
        1 + data->gridSubsetIndices[ ROW ][ LAST  ] -
            data->gridSubsetIndices[ ROW ][ FIRST ]
      : rows;
    const size_t readColumns =
      data->isL3 ?
        1 + data->gridSubsetIndices[ COLUMN ][ LAST  ] -
            data->gridSubsetIndices[ COLUMN ][ FIRST ]
      : columns;

    data->ok =
      readFileData( file,
                    data->product,
                    data->arguments.variable,
                    readRows,
                    readColumns,
                    (const size_t (*)[2]) data->gridSubsetIndices,
                    data->arguments.minimumQuality,
                    data->arguments.maximumCloudFraction,
                    data->arguments.maximumSolarZenithAngle,
                    data->arguments.allowNegativeCounts,
                    data->units,
                    data->isL3 ? data->gridSubsetValues : values,
                    data->isL3 ? data->gridSubsetScratch : scratch ) > 0;
  }
}



/******************************************************************************
PURPOSE: writeGridHeader - Write grid header to stdout.
INPUTS:  const Data* const data  Data description.
******************************************************************************/

static void writeGridHeader( const Data* const data ) {
  const Arguments* const arguments = &data->arguments;
  const int hours = arguments->hours;
  const int yyyymmddhh = arguments->yyyymmddhh;
  const int yyyy       = yyyymmddhh / 1000000;
  const int mm         = yyyymmddhh / 10000 % 100;
  const int dd         = yyyymmddhh / 100 % 100;
  const int hh         = yyyymmddhh % 100;
  const int firstRow    = data->gridSubsetIndices[ ROW ][ FIRST ];
  const int lastRow     = data->gridSubsetIndices[ ROW ][ LAST ];
  const int firstColumn = data->gridSubsetIndices[ COLUMN ][ FIRST ];
  const int lastColumn  = data->gridSubsetIndices[ COLUMN ][ LAST ];
  const int subsetRows = 1 + lastRow - firstRow;
  const int subsetColumns = 1 + lastColumn - firstColumn;

  printf( "SUBSET 9.0 CMAQ\n" );
  printf( "TEMPO_L3\n" );
  printf( "http://tempo.si.edu/,TEMPOSubset\n" );
  printf( "%4d-%02d-%02dT%02d:00:00-0000\n", yyyy, mm, dd, hh );
  printf( "# data dimensions: timesteps variables layers rows columns:\n" );
  printf( "%d 1 1 %d %d\n", hours, subsetRows, subsetColumns );
  printf( "# subset indices (0-based time, 1-based layer/row/column):" );
  printf( " first-timestep last-timestep first-layer last-layer" );
  printf( " first-row last-row" );
  printf( " first-column last-column:\n" );
  printf( "0 %d 1 1 %d %d %d %d\n",
          hours - 1, firstRow, lastRow, firstColumn, lastColumn );
  printf( "# Variable names:\n" );
  printf( "%s\n", arguments->variable );
  printf( "# Variable units:\n" );
  printf( "%s\n", data->units );
  printf( "# lonlat projection: major_semiaxis minor_semiaxis\n" );
  printf( "6370000 6370000\n" );
  printf( "# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[2]:\n" );
  printf( "%lu %lu %f %f %f %f 6 40000 0 1\n",
          data->gridInfo.columns, data->gridInfo.rows,
          data->gridInfo.west, data->gridInfo.south,
          data->gridInfo.cellWidth, data->gridInfo.cellHeight );
  printf( "# IEEE-754 64-bit reals data[variables][timesteps][layers][rows][columns]:\n" );
}



/******************************************************************************
PURPOSE: writeSubsetGridData - Write subset grid data to buffer or stdout.
INPUTS:  Data* const data            Data description.
         const int yyyymmddhh        Timestamp of values[].
OUTPUTS: Data* const data->buffer, data->gridSubsetMeans_yyyymmddhh, data->mask.
NOTES:   First all data before yyyymmddhh is written - including missing values
         for previous timesteps. Then the buffer is initialized or updated
         with these values.
******************************************************************************/

static void writeSubsetGridData( Data* const data,
                                 const int yyyymmddhh ) {

  assert( data ); assert( data->ok ); assert( data->isL3 );
  assert( data->gridInfo.rows >= 1 );
  assert( data->gridSubsetIndices[ 0 ][ 0 ] >= 1 );
  assert( data->gridSubsetValues );
  assert( data->gridSubsetMeans );
  assert( data->gridSubsetCounts );
  assert( isValidYYYYMMDDHHMM( yyyymmddhh * 100LL ) );

  writeDataBeforeYYYYMMDDHH( data, yyyymmddhh );

  if ( data->ok ) {

    if ( data->gridSubsetMeans_yyyymmddhh == 0 ) {
      DEBUG( fprintf( stderr, "updateMeans overwrite  for %d\n", yyyymmddhh );)
      updateMeans( data, 1 );
      data->gridSubsetMeans_yyyymmddhh = yyyymmddhh;
    } else {
      assert( data->gridSubsetMeans_yyyymmddhh == yyyymmddhh );
      DEBUG( fprintf( stderr, "updateMeans accumulate for %d\n", yyyymmddhh );)
      updateMeans( data, 0 );
    }
  }
}



/******************************************************************************
PURPOSE: writeDataBeforeYYYYMMDDHH - Write subset grid data before timestamp
         to data->gridSubsetMeans_yyyymmddhh or stdout.
INPUTS:  Data* const data      data->arguments.yyyymmddhh,
                               data->gridSubsetMeans_yyyymmddhh,
                               data->gridSubsetIndices.
         const int yyyymmddhh  Timestamp to write data before.
OUTPUTS: Data* const Data      data->ok, data->gridSubsetMeans,
                               data->gridSubsetMeans_yyyymmddhh
******************************************************************************/

static void writeDataBeforeYYYYMMDDHH( Data* const data,
                                       const int yyyymmddhh ) {

  const size_t subsetPoints =
    ( 1 + data->gridSubsetIndices[ ROW    ][ LAST  ] -
          data->gridSubsetIndices[ ROW    ][ FIRST ] ) *
    ( 1 + data->gridSubsetIndices[ COLUMN ][ LAST  ] -
          data->gridSubsetIndices[ COLUMN ][ FIRST ] );
  int hours = 0;
  assert( data->ok );
  assert( data->isL3 );
  assert( data->gridInfo.rows );
  assert( data->gridSubsetIndices[0][0] );
  assert( data->gridSubsetMeans );

  DEBUG( fprintf( stderr, "writeDataBeforeYYYYMMDDHH( %d ) with "
                  "data->gridSubsetMeans_yyyymmddhh = %d, "
                  "subsetPoints = %lu, data->ok = %d\n",
                  yyyymmddhh, data->gridSubsetMeans_yyyymmddhh,
                  subsetPoints, data->ok ); )

  if ( data->gridSubsetMeans_yyyymmddhh == 0 ) {
    /* assert( yyyymmddhh >= data->arguments.yyyymmddhh ); */
    hours = hoursUntil( data->arguments.yyyymmddhh, yyyymmddhh );

    DEBUG( fprintf( stderr, "  hours = %d\n", hours ); )

    if ( hours > 0 ) {
      data->ok =
        fillAndWriteArray( MISSING_VALUE, subsetPoints, data->gridSubsetScratch,
                           hours );
      DEBUG( fprintf( stderr, "  wrote %d hours of MISSING_VALUE\n", hours ); )
    }
  } else {
    /* assert( yyyymmddhh >= data->gridSubsetMeans_yyyymmddhh ); */
    hours = hoursUntil( data->gridSubsetMeans_yyyymmddhh, yyyymmddhh );

    DEBUG( fprintf( stderr, "  hours = %d\n", hours ); )

    if ( hours > 0 ) {
      rotate8ByteArrayIfLittleEndian( data->gridSubsetMeans, subsetPoints );
      data->ok = writeArray( data->gridSubsetMeans, subsetPoints, 1 );
      DEBUG( fprintf( stderr, "  wrote 1 hour of data for %d\n",
                      data->gridSubsetMeans_yyyymmddhh ); )
      --hours;

      if ( data->ok && hours > 0 ) {
        data->ok = fillAndWriteArray( MISSING_VALUE, subsetPoints,
                                      data->gridSubsetScratch, hours );
        DEBUG( fprintf( stderr, "  wrote %d hours of MISSING_VALUE\n", hours);)
      }

      data->gridSubsetMeans_yyyymmddhh = 0;
    }
  }

  DEBUG( fprintf( stderr, "writeDataBeforeYYYYMMDDHH( %d ) returning with "
                  "data->ok = %d\n", yyyymmddhh, data->ok ); )
}



/******************************************************************************
PURPOSE: updateMeans - Update gridSubsetMeans with gridSubsetValues.
INPUTS:  Data* const data        data->gridSubsetIndices,
                                 data->gridSubsetValues,
                                 data->gridSubsetMeans.
                                 data->gridSubsetCounts.
         const int reinitialize  If 1 and value is missing then store
                                 MISSING_VALUE.
OUTPUTS: Data* const data        data->gridSubsetMeans updated means.
                                 data->gridSubsetCounts updated mean counts.
******************************************************************************/

static void updateMeans( Data* const data, const int reinitialize ) {

  assert( data->ok );
  assert( data->isL3 );
  assert( data->gridSubsetIndices[0][0] );
  assert( data->gridSubsetValues );
  assert( data->gridSubsetMeans );
  assert( data->gridSubsetCounts );

  {
    const size_t subsetPoints =
      ( 1 + data->gridSubsetIndices[ ROW    ][ LAST  ] -
            data->gridSubsetIndices[ ROW    ][ FIRST ] ) *
      ( 1 + data->gridSubsetIndices[ COLUMN ][ LAST  ] -
            data->gridSubsetIndices[ COLUMN ][ FIRST ] );
    double* values = data->gridSubsetValues;
    double* means = data->gridSubsetMeans;
    unsigned char* counts = data->gridSubsetCounts;
    size_t point = 0;

    for ( point = 0; point < subsetPoints; ++point ) {
      const double value  = values[ point ];
      double mean         = means[ point ];
      unsigned char count = counts[ point ];

      if ( reinitialize || count == 0 ) {

        if ( value > MISSING_VALUE ) {
          means[ point ] = value;
          counts[ point] = 1;
        } else {
          means[ point ] = MISSING_VALUE;
          counts[ point] = 0;
        }
      } else if ( value > MISSING_VALUE ) {

        if ( count < 255 ) {
          assert( mean > MISSING_VALUE );
          mean *= count;
          mean += value;
          ++count;
          mean /= count;
          means[ point ] = mean;
          counts[ point ] = count;
        }
      }
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
  assert( ! data->isL3 );
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

  if ( data->ok ) { /* Append to arrays of timestamps, points: */
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
      memset( buffer, 0, bytes );

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
  assert( ! data->isL3 );
  data->ok = buffer != 0;

  if ( ! buffer ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes to complete the requested action.\n",
             bytes );
  } else {
    memset( buffer, 0, bytes );
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

  /* HACK: Prepend variable name so "column_amount" is "no2_column_amount": */

  const char* const has_o3 = strstr( arguments->variable, "column_amount_o3" );
  const char* const has_column_ =
    has_o3 ? "" : strstr( arguments->variable, "column_" );
  const char* const prefix = has_column_ ? data->variable : "";
  const char* const underscore = has_column_ ? "_" : "";

  assert( ! data->isL3 );

  printf( "Swath 2.0\n%s\n%04d-%02d-%02dT%02d:00:00-0000\n",
          arguments->description, yyyy, mm, dd, hh );
  printf( "# Dimensions: variables timesteps scans:\n%d %d %d\n",
          variables, arguments->hours, data->scans );
  printf( "# Variable names:\n" );
  printf( "Longitude Latitude %s%s%s",
          prefix, underscore, arguments->variable );

  if ( variables == 11 ) {
    printf( " Longitude_SW Longitude_SE Longitude_NW Longitude_NE"
            " Latitude_SW Latitude_SE Latitude_NW Latitude_NE" );
  }

  printf( "\n# Variable units:\ndeg deg %s", data->units );

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




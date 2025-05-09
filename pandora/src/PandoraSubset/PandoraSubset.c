/******************************************************************************
PURPOSE: PandoraSubset.c - Extract a lon-lat subset of data from a list of
         Pandora files and write it to stdout as XDR binary format.

NOTES:   This program calls system( "/usr/bin/sort ..." ).

         Compile:
         gcc -Wall -g -o PandoraSubset PandoraSubset.c Utilities.c \
                   -L../../../lib/$platform \
                   -lm -lc

         Usage:
         PandoraSubset \
           -files <listfile> \
           -tmpdir <temp_directory> \
           -desc "description text" \
           -timerange <yyyymmddhhmmss> <yyyymmddhhmmss> \
           -variable <name> \
           -bounds <minimum_longitude> <minimum_latitude> \
                   <maximum_longitude> <maximum_latitude> \
           -format ascii | xdr \
           [-minimum_quality high | medium | low] (default is high) \
           [-aggregate hourly | daily | monthly | all] (default is none)
           [-minimum_aggregation_count_percentage 0-100] (default 75)

          Example:
          ../../../bin/$platform/PandoraSubset \
          -files testdata/file_list \
          -tmpdir testdata \
          -desc "https://data.pandonia-global-network.org,PandoraSubset" \
          -timerange 20190910000000 20190910235959 \
          -variable nitrogen_dioxide_total_vertical_column_amount \
          -bounds -130 24 -65 50 \
          -format xdr \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 14-line ASCII header followed by binary 64-bit big-endian arrays:

Profile 2.0
https://data.pandonia-global-network.org,PandoraSubset
2019-09-10T00:00:00-0000 2019-09-10T23:59:59-0000
# Subset domain:
-130.0 24.0 -65.0 50.0
# Dimensions: variables profiles:
6 2
# Variable names:
timestamp id longitude latitude elevation nitrogen_dioxide_total_vertical_column_amount
# Variable units:
yyyymmddhhmmss - deg deg m mol/m2
# char notes[profiles][80] and
# MSB 64-bit integers points[profiles] and
# IEEE-754 64-bit reals data_1[variables][points_1] ... data_P[variables][points_P]:
<big-endian binary format array>

HISTORY: 2023-02-01 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/


/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, printf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <ctype.h>     /* For isalpha(). */
#include <stdlib.h>    /* For malloc(), free(), atoll(), atof(), system(). */
#include <limits.h>    /* For INT_MAX. */
#include <float.h>     /* For DBL_MAX. */
#include <math.h>      /* For exp(). */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For LONGITUDE, Bounds, readFile(). */

/*================================= MACROS =================================*/

/* Name of temp per-variable files created in -tmpdir with PID appended: */

#define TEMP_FILE_NAME "junk_PandoraSubset"

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

/*================================ CONSTANTS ================================*/

#define MISSING_VALUE (-999.0)

/*
 * Percentage of maximum number of available values per sensor,
 * per aggregation period. E.g., if aggregating hourly and a sensor reports a
 * value every 2 minutes then the maximum number of values from
 * that sensor would be 60 / 2 = 30. So 75% of 30 = 22.5. Rounding up = 23.
 * So omit that sensor from that hour if it reports less than 23 valid values.
 */

static const double defaultMinimumAggregationCountPercentage = 75.0;

/*---------------------------------------------------------------------------*/

/* Output vars: timestamp, id, longitude, latitude, elevation, var: */

enum { VARIABLES = 6 };

enum { TEMP_FILES = VARIABLES };

static const char* const tempFileNames[ TEMP_FILES ] = {
  "timestamp", "id", "longitude", "latitude", "elevation", "data"
};

enum { FORMAT_XDR, FORMAT_ASCII };
#define FORMAT_STRING "xdr ascii"

enum {
  AGGREGATE_NONE, AGGREGATE_ALL, AGGREGATE_HOURLY, AGGREGATE_DAILY,
  AGGREGATE_MONTHLY
};
#define AGGREGATE_STRING "none all hourly daily monthly"

enum {
  TIMESTAMP_INDEX,
  VARIABLE_INDEX,
  VARIABLE_STRIDE,
  ELEVATION_INDEX,
  ELEVATION_STRIDE,
  L2_QUALITY_FLAG_INDEX,
  COLUMN_INDICES
};

enum { HIGH_QUALITY, MEDIUM_QUALITY, LOW_QUALITY };
#define QUALITY_STRING "high medium low"
#define IS_HIGH_QUALITY(flag) ((flag) == 0 || (flag) == 10)
#define IS_MEDIUM_QUALITY(flag) ((flag) == 1 || (flag) == 11)

/*
 * Maximum number of vertical levels an instrument could measure.
 * Elevated measures beyond this limit will be filtered-out.
 */

enum { MAXIMUM_LEVELS = 256 };

/*================================== TYPES ==================================*/

#include "ColumnInfoTable.h" /* For ColumnInfo type, columnInfoTable[]. */

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

enum { NOTE_LENGTH = 79 };
typedef char Note[ NOTE_LENGTH + 1 ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;       /* File containing list of input files to read.*/
  const char* tmpdir;         /* Name of directory to write temp files. */
  const char* description;    /* User-supplied description. */
  const char* variable;       /* Name of variable to read. */
  Bounds      bounds; /* Subset bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  double      elevationRange[ 2 ]; /* elevationRange[ MINIMUM, MAXIMUM ]. */
  long long   yyyymmddhhmmss[ 2 ]; /* Beginning/ending timestamp of subset. */
  double      minimumAggregationCountPercentage; /* Default is 75%. */
  int         aggregate;       /* AGGREGATE_NONE/ALL/HOURLY/DAILY/MONTHLY. */
  int         minimumQuality; /* HIGH_QUALITY, MEDIUM_QUALITY, LOW_QUALITY. */
  int         format;         /* FORMAT_XDR, FORMAT_ASCII. */
} Arguments;



/*
 * Profile: subsetted result of reading a Pandora data file.
 *
 * # char notes[profiles][80] and
 * # MSB 64-bit integers points[profiles] and
 */

typedef struct {
  long long points; /* # of filtered data points (timesteps x levels). */
  long long id;     /* Constructed from part of data file name. */
  double longitude; /* Of site. */
  double latitude;  /* Of site. */
  double elevation; /* Meters above mean sea level of site. */
  Note note;        /* Site name, location, id parsed from file name. */
} ProfileInfo;


/* Simple linked-list: */

typedef struct Node { void* data; struct Node* next; } Node;

static Node* appendListData( Node* node, void* data ) {
  Node* newNode = NEW( Node, 1 );
  assert( data );

  if ( newNode ) {
    newNode->data = data;
    newNode->next = 0;

    if ( node ) {

      while ( node->next ) {
        node = node->next;
      }

      node->next = newNode;
    }
  }

  return newNode;
}

static void deallocateList( Node* node ) {

  while ( node ) {
    Node* const next = node->next;
    node->next = 0;
    FREE( node->data );
    FREE( node );
    node = next;
  }
}


/* Data type: */

typedef struct {
  Arguments   arguments;       /* User-supplied (command-line) arguments. */
  const char* units;           /* Units of output variable. E.g., "ug/m3". */
  char        fileType[ 32 ];  /* E.g., "_L2Tot_rnvs1p1-7.txt". */
  FileName    tempFileNames[ VARIABLES ]; /* Pathed names of temp files. */
  FILE*       tempFiles[ VARIABLES ];     /* Temp files of output subset data*/
  char*       inputBuffer;                /* Current file content. */
  size_t      bufferSize;                 /* Allocated size of inputBuffer. */
  int         columnInfoIndex; /* Index of fileType/variable in columnInfo[].*/
  Note        note;            /* Site description. */
  size_t      points;          /* Number of valid data points in subset. */
  Node*       profileInfoList; /* List of subsetted profile info. */
  int         ok;              /* Did last command succeed? */
} Data;

/* Data destructor: */

static void deallocateData( Data* data ) {
  assert( data );
  FREE( data->inputBuffer );
  deallocateList( data->profileInfoList );
  memset( data, 0, sizeof *data );
}


/*========================== FORWARD DECLARATIONS ===========================*/


static void createTempFiles( Data* const data );

static void removeTempFiles( Data* const data );

static void openTempFiles( Data* const data );

static void closeTempFiles( Data* const data );


static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static int fileInSubset( const char* const fileName,
                         const Bounds bounds,
                         const double bottom,
                         const double top,
                         const long long yyyymmddhhmmss0,
                         const long long yyyymmddhhmmss1,
                         double* const longitude,
                         double* const latitude,
                         double* const elevation,
                         long long* const id,
                         Note note );

static double parseDoubleFromHeader( const char* header, const char* match );

static double parseTimestampFromHeader( const char* header, const char* match );

static int parseNoteFromHeader( const char* header, Note note);

static const char* parseFileType( const char* const fileName );

static int getVariableIndex( const char* const fileType,
                             const char* const variable );

static void readData( Data* const data );

static void extractSubset( Data* data, ProfileInfo* const profileInfo );

static int parseColumnValues( char* line,
                              const long long yyyymmddhhmmss0,
                              const long long yyyymmddhhmmss1,
                              const double siteElevation,
                              const double minimumElevation,
                              const double maximumElevation,
                              const int minimumQuality,
                              const int columnInfoIndex,
                              long long* yyyymmddhhmmss,
                              double measures[],
                              double elevations[] );

static void initializeAggregated( const int count,
                                  const double measures[],
                                  const double elevations[],
                                  double aggregatedMeasures[],
                                  double aggregatedElevations[],
                                  int aggregatedCounts[] );

static void aggregateData( const int count,
                           const double measures[],
                           const double elevations[],
                           double aggregatedMeasures[],
                           double aggregatedElevations[],
                           int aggregatedCounts[] );

static int aggregationLevel( const double elevation,
                             double* const aggregatedElevation );

static void writeTempData( Data* const data,
                           const long long yyyymmddhhmmss,
                           const long long id,
                           const double longitude,
                           const double latitude,
                           const double elevation,
                           const double measure,
                           const Note note );

static void streamXDRHeader( const Data* const data );

static void streamXDRData( Data* const data );



/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of Pandora files and write
         it to stdout in XDR format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int result = 0;
  Data data;
  memset( &data, 0, sizeof data );
  data.ok = parseArguments( argc, argv, &data.arguments );

  if ( ! data.ok ) {
    printUsage( argv[ 0 ] );
  } else {
    createTempFiles( &data );

    if ( data.ok ) {
      readData( &data ); /* Read input data and write temp files or ASCII. */
      closeTempFiles( &data );

      if ( data.ok && data.points > 0 &&
           data.arguments.format == FORMAT_XDR ) {
        streamXDRHeader( &data );
        streamXDRData( &data );
      }
    }

    removeTempFiles( &data );
    result = data.ok && data.points > 0;

    if ( ! result ) {
      fprintf( stderr, "\n%s: No points were in the subset.\n", argv[ 0 ] );
    }

    deallocateData( &data );
  }

  return ! result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: createTempFiles - Create temp output files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->ok, data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void createTempFiles( Data* const data ) {
  assert( data );

  {
    const int pid = getpid();
    size_t index = 0;
    data->ok = 1;

    for ( index = 0; data->ok && index < TEMP_FILES; ++index ) {
      memset( data->tempFileNames[ index ], 0, sizeof (FileName) );
      assert( data->arguments.tmpdir ); assert( data->arguments.tmpdir[ 0 ] );
      assert( tempFileNames[ index ] ); assert( tempFileNames[ index ][ 0 ] );
      snprintf( data->tempFileNames[ index ],
                sizeof (FileName) / sizeof (char) - 1,
                "%s/%s_%s.%d",
                data->arguments.tmpdir, TEMP_FILE_NAME,
                tempFileNames[ index ], pid );
      data->tempFiles[ index ] = fopen( data->tempFileNames[ index ], "wb" );

      if ( ! data->tempFiles[ index ] ) {
        fprintf( stderr, "\nCan't create temporary output file '%s'.\n",
                 data->tempFileNames[ index ] );
        data->ok = 0;
      }
    }

    if ( ! data->ok ) {
      closeTempFiles( data );
      removeTempFiles( data );
    }
  }
}



/******************************************************************************
PURPOSE: removeTempFiles - Close and remove temp files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void removeTempFiles( Data* const data ) {
  const size_t count =
    data ? sizeof data->tempFiles / sizeof data->tempFiles[0] : 0;
  size_t index = 0;
  assert( data );

  for ( index = 0; index < count; ++index ) {

    if ( data->tempFiles[ index ] ) {
      fclose( data->tempFiles[ index ] );
      data->tempFiles[ index ] = 0;
    }

    if ( data->tempFileNames[ index ] && data->tempFileNames[ index ][ 0 ] ) {
      unlink( data->tempFileNames[ index ] );
      memset( data->tempFileNames[ index ],
              0, sizeof data->tempFileNames[ index ] );
    }
  }
}



/******************************************************************************
PURPOSE: openTempFiles - Open temp files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void openTempFiles( Data* const data ) {
  const size_t count =
    data ? sizeof data->tempFiles / sizeof data->tempFiles[0] : 0;
  size_t index = 0;
  assert( data );


  for ( index = 0; index < count; ++index ) {
    assert( ! data->tempFiles[ index ] );
    data->tempFiles[ index ] = fopen( data->tempFileNames[ index ], "rb" );

    if ( ! data->tempFiles[ index ] ) {
      fprintf( stderr, "\nCan't open temporary output file '%s'.\n",
               data->tempFileNames[ index ] );
      data->ok = 0;
    }
  }
}



/******************************************************************************
PURPOSE: closeTempFiles - Close temp files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void closeTempFiles( Data* const data ) {
  const size_t count =
    data ? sizeof data->tempFiles / sizeof data->tempFiles[0] : 0;
  size_t index = 0;
  assert( data );


  for ( index = 0; index < count; ++index ) {

    if ( data->tempFiles[ index ] ) {
      fclose( data->tempFiles[ index ] );
      data->tempFiles[ index ] = 0;
    }
  }
}



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* name  Name of program.
******************************************************************************/

static void printUsage( const char* name ) {
  assert( name ); assert( *name );
  fprintf( stderr,
           "\n%s - Extract a subset of data from a time-sorted list of\n"
          "Pandora files and write it to stdout in XDR binary format.\n",
           name );
  fprintf( stderr, "Data is subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "  -files <listfile> \\\n" );
  fprintf( stderr, " [-tmpdir <temp_directory>] (default = .)\\\n" );
  fprintf( stderr, " [-desc \"description text\"] (default = pandonia)\\\n" );
  fprintf( stderr, "  -timerange <yyyymmddhhmmss> <yyyymmddhhmmss> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, " [-bounds <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ]\\\n" );
  fprintf( stderr, " [-elevation <minimum_> <maximum>] "
                   "(default = -500 100000)\\\n" );
  fprintf( stderr, " [-format ascii | xdr] (default = xdr)\\\n" );
  fprintf( stderr, " [-minimum_quality high | medium ] (default = high) \\\n" );
  fprintf( stderr, " [-aggregate none | hourly | daily | monthly | all] "
                   "(default = none) \\\n" );
  fprintf( stderr, " [-minimum_aggregation_count_percentage 0-100 "
                   "(default = 75)\n\n" );
  fprintf( stderr, "Note:\ntimes are in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory were temp files are " );
  fprintf ( stderr, "written.\nIt should have enough disk space (100GB).\n" );
  fprintf( stderr, "\nExample:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-files file_list \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc "
                   "\"https://data.pandonia-global-network.org,"
                   "PandoraSubset\" \\\n");
  fprintf( stderr, "-timerange 20190910000000 20190910235959 \\\n" );
  fprintf( stderr,
          "-variable nitrogen_dioxide_total_vertical_column_amount \\\n" );
  fprintf( stderr, "-bounds -130 24 -65 50 \\\n" );
  fprintf( stderr, "-format xdr \\\n" );
  fprintf( stderr, "> subset.xdr\n\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Profile 2.0\n" );
  fprintf( stderr, "https://data.pandonia-global-network.org,PandoraSubset\n");
  fprintf( stderr, "2019-09-10T00:00:00-0000 2019-09-10T23:59:59-0000\n" );
  fprintf( stderr, "# Subset domain <min_lon> <min_lat> <max_lon> <max_lat>:\n");
  fprintf( stderr, "-130 24 -65 50\n" );
  fprintf( stderr, "# Dimensions: variables profiles\n" );
  fprintf( stderr, "6 20\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "timestamp id longitude latitude elevation "
           "nitrogen_dioxide_total_vertical_column_amount\n");
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss - deg deg m mol/cm2\n" );
  fprintf( stderr, "# char notes[profiles][80] and\n" );
  fprintf( stderr, "# MSB 64-bit integers points[profiles] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals "
                   "data_1[variables][points_1] ..."
                   "data_P[variables][points_P]:\n" );
  fprintf( stderr, "<big-endian binary format array>\n" );
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
  int result = argc >= 8;
  static const double elevationRange[ 2 ] = {
    MINIMUM_VALID_SURFACE_ELEVATION_METERS, MAXIMUM_VALID_ELEVATION_METERS
  };
  static const double zero_100[ 2 ] = { 0.0, 100.0 };
  static Option options[] = {
    { "-files",              1, FILE_TYPE,           1, 0, 0, 0, 0 },
    { "-tmpdir",             0, DIRECTORY_TYPE,      1, 0, 0, 0, 0 },
    { "-desc",               0, STRING_TYPE,         1, 0, 0, 0, 0 },
    { "-variable",           1, ENUM_TYPE,           1, 0, 0, 0, 0 },
    { "-timerange",          1, YYYYMMDDHHMMSS_TYPE, 2, 0, 0, 0, 0 },
    { "-format",             0, ENUM_TYPE,           1, 0, FORMAT_STRING, 0, 0},
    { "-bounds",             0, BOUNDS_TYPE,         4, 0, 0, 0, 0 },
    { "-elevation",          0, REAL64_TYPE,         2, elevationRange, 0,0,0 },
    { "-minimum_quality",    0, ENUM_TYPE,           1, 0, QUALITY_STRING, 0,0},
    { "-aggregate",          0, ENUM_TYPE,           1, 0,AGGREGATE_STRING,0,0},
    { "-minimum_aggregation_count_percentage",
                             0, REAL64_TYPE,         1, zero_100, 0, 0, 0 }
  };
  enum { MAXIMUM_VARIABLE_NAME_LENGTH = 80 };
  char variableNames[ sizeof columnInfoTable / sizeof columnInfoTable[0] *
                      MAXIMUM_VARIABLE_NAME_LENGTH ] = "";
  int variableIndex = -1;

  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] ); assert( arguments );

  if ( result ) {

    /* Initialize arguments to defaults: */

    memset( arguments, 0, sizeof (Arguments) );
    arguments->tmpdir = ".";
    arguments->description =
      "https://data.pandonia-global-network.org,PandoraSubset";
    arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
    arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
    arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
    arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
    arguments->elevationRange[ MINIMUM ] =
      MINIMUM_VALID_SURFACE_ELEVATION_METERS;
    arguments->elevationRange[ MAXIMUM ] = MAXIMUM_VALID_ELEVATION_METERS;
    arguments->minimumQuality = HIGH_QUALITY;
    arguments->minimumAggregationCountPercentage =
      defaultMinimumAggregationCountPercentage;

    /* Init space-delimited string of valid variable names for parse check: */

    {
      const size_t variables =
        sizeof columnInfoTable / sizeof columnInfoTable[ 0 ];
      size_t variable = 0;
      memset( variableNames, 0, sizeof variableNames );

      for ( variable = 0; variable < variables; ++variable ) {
        assert( strlen( variableNames ) + MAXIMUM_VARIABLE_NAME_LENGTH + 2
                < sizeof variableNames / sizeof *variableNames );
        assert( ! strchr( columnInfoTable[ variable ].name, ' ' ) );
        strncat( variableNames, columnInfoTable[ variable ].name,
                 MAXIMUM_VARIABLE_NAME_LENGTH );
        strcat( variableNames, " " );
      }
    }

    /* Finish initializing non-compile-time-constant parts of options: */

    options[ 0 ].values = &arguments->listFile;
    options[ 1 ].values = &arguments->tmpdir;
    options[ 2 ].values = &arguments->description;
    options[ 3 ].valids = variableNames; /* Space-delimited var names. */
    options[ 3 ].values = &variableIndex; /* Temp int to get index of name. */
    options[ 4 ].values = &arguments->yyyymmddhhmmss[ 0 ];
    options[ 5 ].values = &arguments->format;
    options[ 6 ].values = &arguments->bounds[ 0 ][ 0 ];
    options[ 7 ].values = &arguments->elevationRange[ 0 ];
    options[ 8 ].values = &arguments->minimumQuality;
    options[ 9 ].values = &arguments->aggregate;
    options[ 10].values = &arguments->minimumAggregationCountPercentage;

    result =
      parseOptions( argc, argv, sizeof options / sizeof *options, options );

    if ( result ) { /* Initialize argument->variable to static string name: */
      arguments->variable = columnInfoTable[ variableIndex ].name;
    }

  } else {
    fprintf( stderr, "\n%s: Invalid/insufficient command-line arguments.\n",
             argv[ 0 ] );
  }

  assert( result == 0 ||
         ( arguments->listFile && arguments->listFile[ 0 ] &&
           arguments->tmpdir && arguments->tmpdir[ 0 ] &&
           arguments->description && arguments->description[ 0 ] &&
           arguments->variable && arguments->variable[ 0 ] &&
           ( arguments->format == FORMAT_XDR ||
             arguments->format == FORMAT_ASCII ) &&
           isValidYYYYMMDDHHMMSS( arguments->yyyymmddhhmmss[ 0 ] ) &&
           isValidYYYYMMDDHHMMSS( arguments->yyyymmddhhmmss[ 1 ] ) &&
           arguments->yyyymmddhhmmss[ 0 ] <= arguments->yyyymmddhhmmss[ 1 ] &&
           isValidBounds( (const double (*)[2]) arguments->bounds ) &&
           IN4( arguments->minimumQuality,
                HIGH_QUALITY, MEDIUM_QUALITY, LOW_QUALITY ) &&
           IN_RANGE( arguments->minimumAggregationCountPercentage,0.0,100.0)));

  return result;
}



/******************************************************************************
PURPOSE: readData - Read data from each listed data file and
         write the subset of data to the temporary files.
INPUTS:  Data* data  Data to read.
OUTPUTS: Data* data  data->points, ok.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  size_t length = 0;
  char* listFileContent = 0;
  data->ok = readFile( arguments->listFile, &length, &listFileContent );

  if ( data->ok ) {
    int wroteSomeData = 0;
    char* inputFileName = 0;
    char* end = 0;

    /* Get each line of list file. It is the data file to read: */

    for ( inputFileName = strtok_r( listFileContent, "\n", &end );
          inputFileName;
          inputFileName = strtok_r( 0, "\n", &end ) ) {

      DEBUG( fprintf( stderr, "Processing Pandora file %s\n", inputFileName );)

      /* Check that file type matches: */

      const char* fileType = parseFileType( inputFileName );

      if ( fileType ) {
        DEBUG( fprintf( stderr, "File type = %s\n", fileType ); )

        if ( data->fileType[ 0 ] == '\0' ) { /* Initialize. */
          memset( data->fileType, 0, sizeof data->fileType );
          strncpy( data->fileType, fileType, sizeof data->fileType - 1 );

          /* Get index of variable in columnInfoTable[]: */

          data->columnInfoIndex =
            getVariableIndex( fileType, arguments->variable );
          data->ok = data->columnInfoIndex > -1;

          if ( ! data->ok ) {
            fprintf( stderr, "\nInvalid variable '%s' for file '%s'\n",
                     arguments->variable, fileType );
            break;
          } else {
            data->units = columnInfoTable[ data->columnInfoIndex ].units;
          }
        }

        if ( ! strcmp( fileType, data->fileType ) ) {

          /* Read file header and check if instrument is within subset: */

          double longitude = 0.0;
          double latitude  = 0.0;
          double elevation = 0.0;
          long long id = 0;
          Note note = "";
          const int inSubset =
            fileInSubset( inputFileName,
                          arguments->bounds,
                          arguments->elevationRange[ MINIMUM ],
                          arguments->elevationRange[ MAXIMUM ],
                          arguments->yyyymmddhhmmss[0],
                          arguments->yyyymmddhhmmss[1],
                          &longitude,
                          &latitude,
                          &elevation,
                          &id,
                          note);

          DEBUG( fprintf( stderr, "inSubset = %d, "
                          "longitude = %f, latitude = %f, elevation = %f, "
                          "id = %lld, note = '%s'\n",
                          inSubset, longitude, latitude, elevation, id, note);)

          if ( inSubset ) {
            data->ok =
              readFile( inputFileName, &data->bufferSize, &data->inputBuffer );

            if ( data->ok ) {
              ProfileInfo* profileInfo = NEW( ProfileInfo, 1 );
              data->ok = profileInfo != 0;

              if ( profileInfo ) {
                profileInfo->longitude = longitude;
                profileInfo->latitude = latitude;
                profileInfo->elevation = elevation;
                profileInfo->id = id;
                strncpy( profileInfo->note, note, sizeof (Note) );

                /* Write data to temp file: */

                extractSubset( data, profileInfo );

                /* If data was extracted then append it else free it: */

                if ( profileInfo->points > 0 ) {
                  Node* newNode =
                    appendListData( data->profileInfoList, profileInfo );
                  data->ok = newNode != 0;

                  if ( newNode == 0 ) {
                    FREE( profileInfo );
                    break;
                  } else if ( data->profileInfoList == 0 ) {
                    data->profileInfoList = newNode;
                  }

                } else {
                  FREE( profileInfo );
                }

                if ( data->ok ) {
                  wroteSomeData = 1;
                }
              } else {
                break;
              }
            }

            if ( ! data->ok ) {
              fprintf( stderr, "\nOmitting invalid file %s\n", inputFileName);
            }
          }
        } else {
          fprintf( stderr, "\nOmitting unmatched file %s\n", inputFileName );
        }
      } else {
        fprintf( stderr, "\nOmitting unknown type of file %s\n", inputFileName);
      }
    } /* End loop on listFile. */

    free( listFileContent ), listFileContent = 0;
    data->ok = wroteSomeData;
  }

  DEBUG( fprintf( stderr,
                  "End of readData() data->points = %lu, data->ok = %d\n",
                  data->points, data->ok ); )
}



/******************************************************************************
PURPOSE: fileInSubset - Check if file's instrument is located within bounds
         and has data within the given time range.
INPUTS:  const char* const fileName  Pathed name of Pandora file.
         const Bounds bounds         Lon-lat bounds to check against.
         const double bottom         Lower elevation (m).
         const double top            Upper elevation (m).
         const int yyyymmddhhmmss0   Beginning timestamp.
         const int yyyymmddhhmmss1   Ending timestamp.
OUTPUTS: double* const longitude     Longitude of site.
         double* const latitude      Latitude of site.
         double* const elevation     Elevation of site (meters above sea level).
         long long id                Integer id of site.
         Note note                   Note includes file and site description.
RETURNS: int 1 if within subset else 0.
NOTES:   Data file headers look like:
File name: Pandora147s1_BronxNY_L2Tot_rnvs0p1-7.txt
File generation date: 20200621T093155Z
Data description: Level 2 total columns file
Data file version: rnvs0p1-7
Data product status: Nitrogen dioxide data are disused
Local principal investigator: Jim Szykman
Network principal investigator: Alexander Cede
Instrument type: Pandora
Instrument number: 147
Spectrometer number: 1
Processing software version used: BlickP v1.7.16
Full location name: New York Botanical Garden NYSDEC (USEPA AQS ID 36-005-0133)
Short location name: BronxNY
Country of location: USA
Location latitude [deg]: 40.8679
Location longitude [deg]: -73.8781
Location altitude [m]: 31
Data start time: 20190910T144614Z
Data end time: 20191228T195136Z
Data caveats: None
-------------------------------------------------------------------------------
******************************************************************************/

static int fileInSubset( const char* const fileName,
                         const Bounds bounds,
                         const double bottom,
                         const double top,
                         const long long yyyymmddhhmmss0,
                         const long long yyyymmddhhmmss1,
                         double* const longitude,
                         double* const latitude,
                         double* const elevation,
                         long long* const id,
                         Note note ) {
  int result = 0;

  assert( fileName ); assert( *fileName );
  assert( isValidBounds( bounds ) );
  assert( IN_RANGE( bottom,
                    MINIMUM_VALID_SURFACE_ELEVATION_METERS,
                    MAXIMUM_VALID_ELEVATION_METERS ) );
  assert( IN_RANGE( top, bottom, MAXIMUM_VALID_ELEVATION_METERS ) );
  assert( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss0 ) );
  assert( isValidYYYYMMDDHHMMSS( yyyymmddhhmmss1 ) );
  assert( yyyymmddhhmmss0 <= yyyymmddhhmmss1 );
  assert( longitude ); assert( latitude ); assert( elevation );
  assert( id ); assert( note );

  {
    FILE* file = fopen( fileName, "rb" );

    if ( file ) {
      enum { SIZE = 2048 }; /* Enough to read file header. */
      char buffer[ SIZE ] = "";
      int ok = fread( buffer, sizeof (char), SIZE, file ) == SIZE;

      if ( ok ) {
        *latitude = parseDoubleFromHeader( buffer, "\nLocation latitude" );
        result =
          IN_RANGE( *latitude,
                    bounds[ LATITUDE ][ MINIMUM ],
                    bounds[ LATITUDE ][ MAXIMUM ] );

        if ( result ) {
          *longitude = parseDoubleFromHeader( buffer, "\nLocation longitude");
          result =
            IN_RANGE( *longitude,
                      bounds[ LONGITUDE ][ MINIMUM ],
                      bounds[ LONGITUDE ][ MAXIMUM ] );

          if ( result ) {
            *elevation = parseDoubleFromHeader(buffer, "\nLocation altitude");
            result = IN_RANGE( *elevation, bottom, top );

            if ( result ) {
              const long long start =
                parseTimestampFromHeader( buffer, "\nData start time" );
              result = start <= yyyymmddhhmmss1;

              if ( result ) {
                const long long end =
                  parseTimestampFromHeader( buffer, "\nData end time" );
                result = end == LLONG_MIN || end >= yyyymmddhhmmss0;


                if ( result ) {
                  result = parseNoteFromHeader( buffer, note );

                  if ( result ) {
                    const char* const slash = strrchr( fileName, '/' );
                    const char* const name = slash ? slash + 1 : fileName;
                    const char* const tag = "Pandora";
                    const char* const found = strstr( name, tag );

                    if ( found ) {
                      *id = atoll( found + strlen( tag ) );
                      result = *id > 0;
                    }
                  }
                }
              }
            }
          }
        }
      }

      fclose( file ), file = 0;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: parseDoubleFromHeader - Parse a double value from data file header.
INPUTS:  const char* header  File header.
         const char* match   Line to match. E.g., "\nLocation latitude".
RETURNS: double parsed value or -DBL_MAX if unsuccessful.
NOTES:   Data file headers look like:
File name: Pandora147s1_BronxNY_L2Tot_rnvs0p1-7.txt
File generation date: 20200621T093155Z
Data description: Level 2 total columns file
Data file version: rnvs0p1-7
Data product status: Nitrogen dioxide data are disused
Local principal investigator: Jim Szykman
Network principal investigator: Alexander Cede
Instrument type: Pandora
Instrument number: 147
Spectrometer number: 1
Processing software version used: BlickP v1.7.16
Full location name: New York Botanical Garden NYSDEC (USEPA AQS ID 36-005-0133)
Short location name: BronxNY
Country of location: USA
Location latitude [deg]: 40.8679
Location longitude [deg]: -73.8781
Location altitude [m]: 31
Data start time: 20190910T144614Z
Data end time: 20191228T195136Z
Data caveats: None
-------------------------------------------------------------------------------
******************************************************************************/

static double parseDoubleFromHeader( const char* header, const char* match ) {

  double result = -DBL_MAX;
  assert( header ); assert( match ); assert( *match );
  header = strstr( header, match );

  if ( header ) {
    header = strchr( header, ':' );

    if ( header ) {
      result = atof( header + 1 );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: parseTimestampFromHeader - Parse a timestamp from data file header.
INPUTS:  const char* header  File header.
         const char* match   Line to match. E.g., "\nData start time".
RETURNS: long long parsed value or LLONG_MIN if unsuccessful.
NOTES:   Data file headers look like:
File name: Pandora147s1_BronxNY_L2Tot_rnvs0p1-7.txt
File generation date: 20200621T093155Z
Data description: Level 2 total columns file
Data file version: rnvs0p1-7
Data product status: Nitrogen dioxide data are disused
Local principal investigator: Jim Szykman
Network principal investigator: Alexander Cede
Instrument type: Pandora
Instrument number: 147
Spectrometer number: 1
Processing software version used: BlickP v1.7.16
Full location name: New York Botanical Garden NYSDEC (USEPA AQS ID 36-005-0133)
Short location name: BronxNY
Country of location: USA
Location latitude [deg]: 40.8679
Location longitude [deg]: -73.8781
Location altitude [m]: 31
Data start time: 20190910T144614Z
Data end time: 20191228T195136Z
Data caveats: None
-------------------------------------------------------------------------------
******************************************************************************/

static double parseTimestampFromHeader( const char* header, const char* match) {

  long long result = LLONG_MIN;
  assert( header ); assert( match ); assert( *match );
  header = strstr( header, match );

  if ( header ) {
    header = strchr( header, ':' );

    if ( header ) {
      long long yyyymmdd = atoi( header + 1 );
      yyyymmdd *= 1000000;

      if ( isValidYYYYMMDDHHMMSS( yyyymmdd ) ) {
        header = strchr( header, 'T' );

        if ( header ) {
          const long long hhmmss = atoi( header + 1 );

          if ( IN_RANGE( hhmmss, 0, 235959 ) ) {
            result = yyyymmdd + hhmmss;
          }
        }
      }
    }
  }

  assert( result == LLONG_MIN || isValidYYYYMMDDHHMMSS( result ) );
  return result;
}



/******************************************************************************
PURPOSE: parseNoteFromHeader - Parse a file name and site description from
         data file header.
INPUTS:  const char* header  File header.
OUTPUTS: Note note           File name and site description. E.g.,:
"Pandora147s1_BronxNY_L2Tot_rnvs0p1-7;Garden NYSDEC (USEPA AQS ID 36-005-0133)"
RETURNS: int 1 if successful, else 0.
NOTES:   Data file headers look like:
File name: Pandora147s1_BronxNY_L2Tot_rnvs0p1-7.txt
File generation date: 20200621T093155Z
Data description: Level 2 total columns file
Data file version: rnvs0p1-7
Data product status: Nitrogen dioxide data are disused
Local principal investigator: Jim Szykman
Network principal investigator: Alexander Cede
Instrument type: Pandora
Instrument number: 147
Spectrometer number: 1
Processing software version used: BlickP v1.7.16
Full location name: New York Botanical Garden NYSDEC (USEPA AQS ID 36-005-0133)
Short location name: BronxNY
Country of location: USA
Location latitude [deg]: 40.8679
Location longitude [deg]: -73.8781
Location altitude [m]: 31
Data start time: 20190910T144614Z
Data end time: 20191228T195136Z
Data caveats: None
-------------------------------------------------------------------------------
******************************************************************************/

static int parseNoteFromHeader( const char* header, Note note) {
  int result = 0;
  assert( header ); assert( note );
  memset( note, 0, sizeof (Note) );

  {
    const int last = sizeof (Note) - 1;
    const char* match = "File name: ";
    const char* s = strstr( header, match );

    if ( s ) {
      int skip = strlen( match );
      int length = 0;
      char c = 0;
      s += skip;

      /* Copy file name (except .txt extension) to first part of note: */

      for ( c = *s;
            length < last && ( isalnum( c ) || IN3( c, '_', '-' ) );
            ++length, ++s, c = *s ) {
        note[ length ] = c;
      }

      if ( length < last ) {

        /* Append semi-colon delimiter: */

        note[ length ] = ';';
        ++length;

        if ( length < last ) {

          /* Copy last part of site description to remaining part of note: */

          match = "Full location name: ";
          s = strstr( header, match );

          if ( s ) {
            skip = strlen( match );
            s += skip;

            {
              const char* s2 = s;

              /* Find end of site description: */

              for ( c = *s2; ! IN4( c, '\r', '\n', '\0' ); c = *s2 ) {
                ++s2;
              }

              {
                int remainder = last - length;

                /* Find end portion of description that will fit: */

                while ( s2 != s && remainder-- ) {
                  --s2;
                }

                /* If in the middle of a word then skip to next word: */

                for ( c = *s2; ! IN5( c, ' ', '\r', '\n', '\0' ); c = *s2 ) {
                  ++s2;
                }

                /* Append end of site description: */

                for ( c = *s2;
                      length < last && ! IN4( c, '\r', '\n', '\0' );
                      ++length, ++s2, c = *s2 ) {
                  note[ length ] = c;
                }
              }
            }
          }
        }
      }
    }
  }

  result = *note != '\0';
  assert( note[ sizeof (Note) - 1 ] == '\0' );
  assert( result == 0 || *note );
  return result;
}



/******************************************************************************
PURPOSE: parseFileType - Parse file type from file name.
INPUTS:  const char* const fileName  Name of input file.
RETURNS: const char* file type. E.g., "L2Tot_rnvs1p1-7.txt" else 0 if invalid.
******************************************************************************/

static const char* parseFileType( const char* const fileName ) {
  const char* result = 0;
  const char* s = 0;

  assert( fileName ); assert( *fileName );

  s = strrchr( fileName, '/' );

  if ( s ) {
    ++s;
  } else {
    s = fileName;
  }

  s = strstr( s, "_L2" );

  if ( s ) {
    result = s;
  }

  return result;
}



/******************************************************************************
PURPOSE: getVariableIndex - Get index into columnInfoTable[] of matching file
         type and variable name.
INPUTS:  const char* const fileType  Type of file. E.g., "_L2_rnvh3p1-8.txt".
         const char* const variable  Variable name. E.g., "temperature".
RETURNS: int 0-based index or -1 if not found.
******************************************************************************/

static int getVariableIndex( const char* const fileType,
                             const char* const variable ) {

  const size_t count = sizeof columnInfoTable / sizeof columnInfoTable[ 0 ];
  int result = 0;
  assert( variable ); assert( variable[0] );
  assert( fileType ); assert( fileType[0] );

  for ( result = 0; result < count; ++result ) {
    const ColumnInfo* const columnInfo = columnInfoTable + result;

    if ( ! strcmp( columnInfo->type, fileType ) ) {

      if ( ! strcmp( columnInfo->name, variable ) ) {
        break;
      }
    }
  }

  if ( result == count ) {
    result = -1;
  }

  return result;
}



/******************************************************************************
PURPOSE: extractSubset - Parse variable column data from inputBuffer and write
         to temp file if XDR or stdout if ASCII format.
INPUTS:  Data* data  data->columnInfoIndex  Index into columnInfoTable[] of
                                            variable.
                     data->inputBuffer Content of input data file.
         ProfileInfo* const profileInfo  Info on profile being read.
OUTPUTS: Data* data  data->tempFiles[ TEMP_FILES ]  Write variables .
                     data->points  Increased if AGGREGSTE_NONE, else unchanged.
                     data->ok = 1 if successful, else 0.
         ProfileInfo* const profileInfo   Initialized profileInfo->points.
******************************************************************************/

static void extractSubset( Data* data, ProfileInfo* const profileInfo ) {
  const char* const tag = "\n-------";
  char* beginSection = 0; /* Beginning of section to parse. */
  assert( data );
  assert( data->arguments.variable ); assert( data->arguments.variable[ 0 ] );
  assert( data->inputBuffer );
  assert( data->tempFiles[ 0 ] );
  data->ok = 0;

  beginSection = strstr( data->inputBuffer, tag );

  if ( beginSection ) {
    beginSection = strstr( beginSection + 1, tag );

    if ( beginSection ) {
      beginSection = strchr( beginSection + 1, '\n' );

      if ( beginSection ) {
        ++beginSection;

        {
          const Arguments* const arguments = &data->arguments;
          const long long yyyymmddhhmmss0 = arguments->yyyymmddhhmmss[ 0 ];
          const long long yyyymmddhhmmss1 = arguments->yyyymmddhhmmss[ 1 ];
          long long yyyymmddhhmmss = 0; /* Timestamp of line read. */
          long long yyyymmddhh0000 = 0; /* Current aggregation timestamp*/
          const double minimumElevation = arguments->elevationRange[ 0 ];
          const double maximumElevation = arguments->elevationRange[ 1 ];
          const int aggregate = arguments->aggregate;
          const int minimumQuality = arguments->minimumQuality;
          const long long id = profileInfo->id;
          const double longitude = profileInfo->longitude;
          const double latitude = profileInfo->latitude;
          const double siteElevation = profileInfo->elevation;
          char* line = beginSection;
          char* endLine = 0;
          double measures[ MAXIMUM_LEVELS ];
          double elevations[ MAXIMUM_LEVELS ];
          double aggregatedMeasures[ MAXIMUM_LEVELS ];
          double aggregatedElevations[ MAXIMUM_LEVELS ];
          int aggregatedCounts[ MAXIMUM_LEVELS ];
          int initialized = 0;
          int point = 0;
          int level = 0;

          memset( measures, 0, sizeof measures);
          memset( elevations, 0, sizeof elevations);
          memset( aggregatedMeasures, 0, sizeof aggregatedMeasures);
          memset( aggregatedElevations, 0, sizeof aggregatedElevations);
          memset( aggregatedCounts, 0, sizeof aggregatedCounts);

          DEBUG( fprintf( stderr,
                         "columnInfoIndex = %d: "
                         "columnInfo: index = %d, indexStride = %d, "
                         "elevationIndex = %d, elevationStride = %d, "
                         "filterIndex = %d\n",
                         data->columnInfoIndex,
                         columnInfoTable[data->columnInfoIndex].index,
                         columnInfoTable[data->columnInfoIndex].indexStride,
                         columnInfoTable[data->columnInfoIndex].elevationIndex,
                         columnInfoTable[data->columnInfoIndex].elevationStride,
                         columnInfoTable[data->columnInfoIndex].filterIndex);)

          /* Read each line until EOF or timestamp is beyond yyyymmddhhmmss1:*/

          for ( data->ok = 1;
                data->ok && nextLine( &line, &endLine ) &&
                yyyymmddhhmmss <= yyyymmddhhmmss1;
                *endLine = '\n', line = endLine + 1 ) {
            const int pointsInSubset =
              parseColumnValues( line, yyyymmddhhmmss0, yyyymmddhhmmss1,
                                 siteElevation,
                                 minimumElevation, maximumElevation,
                                 minimumQuality, data->columnInfoIndex,
                                 &yyyymmddhhmmss, measures, elevations );

            if ( pointsInSubset > 0 ) {
              DEBUG( fprintf( stderr,
                              "pointsInSubset = %d, yyyymmddhhmmss = %lld\n",
                              pointsInSubset, yyyymmddhhmmss );)

              if ( aggregate != AGGREGATE_NONE ) {

                if ( ! initialized ) {
                  initialized = 1;
                  yyyymmddhh0000 = yyyymmddhhmmss;
                  initializeAggregated( pointsInSubset, measures, elevations,
                                        aggregatedMeasures,
                                        aggregatedElevations,
                                        aggregatedCounts );
                } else {
                  const int doAggregate =
                    aggregate == AGGREGATE_ALL ||
                    ( aggregate == AGGREGATE_HOURLY &&
                      yyyymmddhhmmss / 10000 == yyyymmddhh0000 / 10000 ) ||
                    ( aggregate == AGGREGATE_DAILY &&
                      yyyymmddhhmmss / 1000000 == yyyymmddhh0000 / 1000000 ) ||
                    ( aggregate == AGGREGATE_MONTHLY &&
                      yyyymmddhhmmss / 100000000 == yyyymmddhh0000 / 100000000);

                  if ( doAggregate ) {
                    aggregateData( pointsInSubset, measures, elevations,
                                   aggregatedMeasures,
                                   aggregatedElevations,
                                   aggregatedCounts );
                  } else {

                    /* Write aggregated values: */

                    for ( level = 0; level < MAXIMUM_LEVELS; ++level ) {
                      const int count = aggregatedCounts[ level ];

                      if ( count ) {
                        const double aggregatedMeasure =
                          aggregatedMeasures[ level ];
                        const double aggregatedElevation =
                          aggregatedElevations[ level ];
                        assert( aggregatedMeasure != MISSING_VALUE );
                        assert( aggregatedElevation >= minimumElevation );
                        assert( aggregatedElevation <= maximumElevation );

                        writeTempData( data,
                                       yyyymmddhh0000, id,
                                       longitude, latitude,
                                       aggregatedElevation,
                                       aggregatedMeasure,
                                       profileInfo->note );
                        profileInfo->points += data->ok;

                        DEBUG( fprintf( stderr,
                                        "wrote aggregatedMeasures[%d] = %e "
                                        "@ %f\n",
                                        level, aggregatedMeasure,
                                        aggregatedElevation ); )

                      }
                    }

                    /* Re-initialize aggregated values to last value read: */

                    yyyymmddhh0000 = yyyymmddhhmmss;
                    initializeAggregated( pointsInSubset, measures, elevations,
                                          aggregatedMeasures,
                                          aggregatedElevations,
                                          aggregatedCounts );

                  }
                }
              } else { /* Write non-aggregated values: */

                for ( point = 0; point < pointsInSubset; ++point ) {
                  const double measure = measures[ point ];
                  const double elevation = elevations[ point ];
                  assert( measure != MISSING_VALUE );
                  assert( elevation >= minimumElevation );
                  assert( elevation <= maximumElevation );

                  writeTempData( data,
                                 yyyymmddhhmmss, id,
                                 longitude, latitude,
                                 elevation,
                                 measure,
                                 profileInfo->note );
                  profileInfo->points += data->ok;

                  DEBUG( fprintf( stderr, "wrote measures[%d] = %e @ %f\n",
                                  point, measure, elevation ); )
                }
              }
            }
          }

          if ( aggregate != AGGREGATE_NONE ) { /* Write final values: */

            for ( level = 0; level < MAXIMUM_LEVELS; ++level ) {
              const int count = aggregatedCounts[ level ];

              if ( count ) {
                const double aggregatedMeasure =
                  aggregatedMeasures[ level ];
                const double aggregatedElevation =
                  aggregatedElevations[ level ];
                writeTempData( data,
                               yyyymmddhh0000, id,
                               longitude, latitude,
                               aggregatedElevation,
                               aggregatedMeasure,
                               profileInfo->note );
                profileInfo->points += data->ok;

                DEBUG( fprintf( stderr,
                                "wrote aggregatedMeasures[%d] = %e "
                                "@ %f\n",
                                level, aggregatedMeasure,
                                aggregatedElevation ); )

              }
            }
          }
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: parseColumnValues - Parse column values from a line.
INPUTS:  const char* line                  Data line to parse.
         const long long yyyymmddhhmmss0   Beginning timestamp of subset.
         const long long yyyymmddhhmmss1   Ending    timestamp of subset.
         const double siteElevation        Site meters above mean sea level.
         const double minimumElevation     Minimum elevation of subset.
         const double maximumElevation     Maximum elevation of subset.
         const int minimumQuality          HIGH_QUALITY, etc.
         const int columnInfoIndex       Index of variable in columnInfoTable[]
OUTPUTS: long long* yyyymmddhhmmss         Read timestamp.
         double measures[MAXIMUM_LEVELS]   Read filtered values.
         double elevations[MAXIMUM_LEVELS] Read elevations (m) above sea level.
RETURNS: int number of values (after filtering) in subset.
******************************************************************************/

static int parseColumnValues( char* line,
                              const long long yyyymmddhhmmss0,
                              const long long yyyymmddhhmmss1,
                              const double siteElevation,
                              const double minimumElevation,
                              const double maximumElevation,
                              const int minimumQuality,
                              const int columnInfoIndex,
                              long long* yyyymmddhhmmss,
                              double measures[],
                              double elevations[] ) {
  int result = 0;
  int elevationCount = 0;
  int ok = 0;

  assert( line );
  assert( yyyymmddhhmmss0 > 0 );
  assert( yyyymmddhhmmss1 > yyyymmddhhmmss0 );
  assert( IN_RANGE( siteElevation,
                    MINIMUM_VALID_SURFACE_ELEVATION_METERS,
                    MAXIMUM_VALID_SURFACE_ELEVATION_METERS ) );
  assert( IN_RANGE( minimumElevation,
                    MINIMUM_VALID_SURFACE_ELEVATION_METERS,
                    MAXIMUM_VALID_ELEVATION_METERS ) );
  assert( IN_RANGE( maximumElevation,
                    minimumElevation,
                    MAXIMUM_VALID_ELEVATION_METERS ) );
  assert( IN4( minimumQuality, HIGH_QUALITY, MEDIUM_QUALITY, LOW_QUALITY ) );
  assert( columnInfoIndex >= 0 );
  assert( columnInfoIndex < sizeof columnInfoTable / sizeof *columnInfoTable );
  assert( yyyymmddhhmmss );
  assert( measures );
  assert( elevations );

  /* Initialize outputs: */

  *yyyymmddhhmmss = 0;

  for ( result = 0; result < MAXIMUM_LEVELS; ++result ) {
    measures[ result ] = MISSING_VALUE;
    elevations[ result ] = MISSING_VALUE;
  }

  result = 0;

  {
    const ColumnInfo* const columnInfo = columnInfoTable + columnInfoIndex;
    const int index = columnInfo->index;
    const int stride = columnInfo->indexStride;
    const int filter = columnInfo->filterIndex;
    const int elevationIndex = columnInfo->elevationIndex;
    const int elevationStride = columnInfo->elevationStride;
    int column = 0;
    char* word = line;
    char* end = 0;

    for ( ok = 1, column = 1, end = strchr( word + 1, ' ' );
          ok && word;
          word = end, end = word ? strchr( word + 1, ' ' ) : 0, ++column ) {
      char* end2 = 0;

      if ( column == 1 ) { /* Read timestamp like 20210304T154639.1Z */
        const long long yyyymmdd = strtoll( word, &end2, 10 );
        ok = *end2 == 'T';

        if ( ok ) {
          const char* const rest = end2 + 1;
          const long long hhmmss = strtoll( rest, &end2, 10 );
          ok = end2 != rest;

          if ( ok ) {
            *yyyymmddhhmmss = yyyymmdd * 1000000 + hhmmss;
            ok =
              isValidYYYYMMDDHHMMSS( *yyyymmddhhmmss ) &&
              IN_RANGE( *yyyymmddhhmmss, yyyymmddhhmmss0, yyyymmddhhmmss1 );

            DEBUG( if(ok) fprintf( stderr, "Read timetamp = %lld\n",
                                   *yyyymmddhhmmss ); )
          }
        }
      } else {

        if ( column == index ||
             ( stride > 0 && column > index &&
               ( column - index ) % stride == 0 ) ) { /* Read measure: */
          double measure = strtod( word, &end2 );
          ok = end2 != word;
          DEBUG( fprintf( stderr, "Read measure = %e\n", measure ); )

          if ( ok ) {
            measure += columnInfo->offset;
            measure *= columnInfo->scale;

            if ( columnInfo->converter ) {
              measure = columnInfo->converter( measure );
            }

            ok = IN_RANGE( measure, columnInfo->minimum, columnInfo->maximum );


            DEBUG( fprintf( stderr, "Final measure = %e, ok = %d\n",
                            measure, ok  ); )

            if ( ok ) {

              if ( result < MAXIMUM_LEVELS ) {
                measures[ result ] = measure;
                DEBUG( fprintf( stderr, "Stored measures[ %d ] = %e\n", result,
                                measures[ result ]  ); )
                ++result;
              }
            }
          }
        } else if ( column == elevationIndex ||
                    ( elevationStride > 0 && column > elevationIndex &&
                      ( column - elevationIndex ) % elevationStride == 0 ) ) {

          /* Read elevation and convert km to m and add siteElevation (m): */

          double elevation = strtod( word, &end2 );
          ok = end2 != word;
          DEBUG( fprintf( stderr, "Read elevation = %f\n", elevation ); )

          if ( ok ) {
            elevation *= KM_TO_M;
            elevation += siteElevation;
            ok = IN_RANGE( elevation, minimumElevation, maximumElevation );
            DEBUG( fprintf( stderr, "Final elevation = %f, ok = %d\n",
                            elevation, ok  ); )

            if ( ok ) {

              if ( elevationCount < MAXIMUM_LEVELS ) {
                elevations[ elevationCount ] = elevation;
                DEBUG( fprintf( stderr, "Stored elevations[ %d ] = %f\n",
                                elevationCount,
                                elevations[ elevationCount ]  ); )
                ++elevationCount;
              }
            }
          }
        } else if ( column == filter ) {
          const long value = strtol( word, &end2, 10 );
          ok = end2 != word;

          if ( ok && minimumQuality != LOW_QUALITY ) {
            ok = IS_HIGH_QUALITY( value );

            if ( ! ok && minimumQuality == MEDIUM_QUALITY ) {
              ok = IS_MEDIUM_QUALITY( value );
            }
          }

          DEBUG( fprintf( stderr, "filter value = %ld, ok = %d\n",
                          value, ok ); )
        }
      }
    }
  }

  /* If no (above surface) elevations were read then store siteElevation: */

  if ( ok ) {
    int point = 0;

    for ( point = 0; point < result; ++point ) {

      if ( measures[ point ] != MISSING_VALUE &&
           elevations[ point ] == MISSING_VALUE ) {
        elevations[ point ] = siteElevation;
      }
    }
  } else {
    result = 0;
  }


  return result;
}


/******************************************************************************
PURPOSE: aggregationLevel - Compute the aggregation level index (0-based) of
         a given elevation and return the aggregation elevation.
INPUTS:  const double elevation Elevation to bin.
OUTPUTS: double* const aggregatedElevation  Aggregated elevation.
RETURNS: int index or -1 if outside aggregation elevation range.
******************************************************************************/

static int aggregationLevel( const double elevation,
                             double* const aggregatedElevation ) {
  int result = -1;

  assert( elevation >= MINIMUM_VALID_SURFACE_ELEVATION_METERS );
  assert( elevation <= MAXIMUM_VALID_ELEVATION_METERS );
  assert( aggregatedElevation );

  if ( elevation <= MAXIMUM_AGGREGATION_ELEVATION ) {
    static const double recipricol = 1.0 /
    ( MAXIMUM_AGGREGATION_ELEVATION - MINIMUM_VALID_SURFACE_ELEVATION_METERS );
    double normalized =
      ( elevation - MINIMUM_VALID_SURFACE_ELEVATION_METERS ) * recipricol;
    result = normalized * AGGREGATION_LEVELS;

    if ( result >= AGGREGATION_LEVELS ) {
      --result;
    }

    *aggregatedElevation = MINIMUM_VALID_SURFACE_ELEVATION_METERS +
      ( result + 1 ) * AGGREGATION_LEVEL_THICKNESS_METERS;
  }

  assert( result == -1 || IN_RANGE( result, 0, MAXIMUM_LEVELS - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: initializeAggregated - Initialize aggregated data to given measures.
INPUTS:  const int count                   Number of measures/elevations.
         const double measures[ count ]    Measures to initialize to.
         const double elevations[ count ]  Elevations to bin.
OUTPUTS: double aggregatedMeasures[ MAXIMUM_LEVELS ]    Initialized to measures
         double aggregatedElevations[ MAXIMUM_LEVELS ]  Initialized.
         int aggregatedCounts[ MAXIMUM_LEVELS ]         Initialized.
******************************************************************************/

static void initializeAggregated( const int count,
                                  const double measures[],
                                  const double elevations[],
                                  double aggregatedMeasures[],
                                  double aggregatedElevations[],
                                  int aggregatedCounts[] ) {

  int index = 0;
  assert( count > 0 ); assert( count <= MAXIMUM_LEVELS );
  assert( measures );
  assert( measures[ 0 ] != MISSING_VALUE );
  assert( measures[ count - 1 ] != MISSING_VALUE );
  assert( elevations );
  assert( elevations[ 0 ] != MISSING_VALUE );
  assert( elevations[ count - 1 ] != MISSING_VALUE );
  assert( aggregatedMeasures );
  assert( aggregatedElevations );
  assert( aggregatedCounts );

  for ( index = 0; index < MAXIMUM_LEVELS; ++index ) {
    aggregatedMeasures[ index ] = MISSING_VALUE;
    aggregatedElevations[ index ] = MISSING_VALUE;
    aggregatedCounts[ index ] = 0;
  }

  for ( index = 0; index < count; ++index ) {
    const double elevation = elevations[ index ];
    double aggregatedElevation = 0.0;
    const int level = aggregationLevel( elevation, &aggregatedElevation );

    if ( level >= 0 ) {
      const double measure = measures[ index ];
      assert( measure != MISSING_VALUE );
      assert( level < MAXIMUM_LEVELS );
      aggregatedMeasures[ level ] = measure;
      aggregatedElevations[ level ] = aggregatedElevation;
      aggregatedCounts[ level ] = 1;

      DEBUG( fprintf( stderr, "initialized aggregatedMeasures[%d] = %e @ %f\n",
                      level, aggregatedMeasures[ level ],
                      aggregatedElevations[ level ] ); )

    }
  }
}



/******************************************************************************
PURPOSE: aggregateData - Aggregate measures.
INPUTS:  const int count                   Number of measures/elevations.
         const double measures[ count ]    Measures to initialize to.
         const double elevations[ count ]  Elevations to bin.
OUTPUTS: double aggregatedMeasures[ MAXIMUM_LEVELS ]    Aggregated to measures
         double aggregatedElevations[ MAXIMUM_LEVELS ]  Aggregated.
         int aggregatedCounts[ MAXIMUM_LEVELS ]         Aggregated.
******************************************************************************/

static void aggregateData( const int count,
                           const double measures[],
                           const double elevations[],
                           double aggregatedMeasures[],
                           double aggregatedElevations[],
                           int aggregatedCounts[] ) {

  int index = 0;
  assert( count > 0 ); assert( count <= MAXIMUM_LEVELS );
  assert( measures );
  assert( measures[ 0 ] != MISSING_VALUE );
  assert( measures[ count - 1 ] != MISSING_VALUE );
  assert( elevations );
  assert( elevations[ 0 ] != MISSING_VALUE );
  assert( elevations[ count - 1 ] != MISSING_VALUE );
  assert( aggregatedMeasures );
  assert( aggregatedElevations );
  assert( aggregatedCounts );

  for ( index = 0; index < count; ++index ) {
    const double elevation = elevations[ index ];
    double aggregatedElevation = 0.0;
    const int level = aggregationLevel( elevation, &aggregatedElevation );
    assert( elevation != MISSING_VALUE );

    if ( level >= 0 ) {
      const double measure = measures[ index ];
      const int aggregatedCount1 = aggregatedCounts[ level ] + 1;
      assert( measure != MISSING_VALUE );
      assert( IN_RANGE( aggregatedElevation,
                        MINIMUM_VALID_SURFACE_ELEVATION_METERS,
                        MAXIMUM_AGGREGATION_ELEVATION ) );
      assert( level < MAXIMUM_LEVELS );

      if ( aggregatedCount1 == 1 ) {
        aggregatedMeasures[ level ] = measure;
        aggregatedElevations[ level ] = aggregatedElevation;
      } else {
        assert( aggregatedMeasures[ level ] != MISSING_VALUE );
        aggregatedMeasures[ level ] *= aggregatedCounts[ level ];
        aggregatedMeasures[ level ] += measure;
        aggregatedMeasures[ level ] /= aggregatedCount1;
        assert( aggregatedElevations[ level ] == aggregatedElevation );
      }

      aggregatedCounts[ level ] = aggregatedCount1;

      DEBUG( fprintf( stderr, "aggregated aggregatedMeasures[%d] = %e @ %f\n",
                      level, aggregatedMeasures[ level ],
                      aggregatedElevations[ level ] ); )
    }
  }
}



/******************************************************************************
PURPOSE: writeTempData - Write data to per-variable temp files if XDR format or
         write ASCII format data.
INPUTS:  Data* data  data->tempFiles[ TEMP_FILES ].
                     data->points        Previous number of output points.
         const long long yyyymmddhhmmss  Timestamp of data point.
         const long long id              Id of data point.
         const double longitude          Longitude of data point.
         const double latitude           Latitude of data point.
         const double elevation          Elevation of data point.
         const double measure            Measure of data point.
         const Note note                 Note of profile.
OUTPUTS: Data* data  data->tempFiles[ TEMP_AGGREGATED_FILE ] Data appended.
                     data->points  Increased number of output points.
                     data->ok = 1 if successful, else 0.
******************************************************************************/

static void writeTempData( Data* const data,
                           const long long yyyymmddhhmmss,
                           const long long id,
                           const double longitude,
                           const double latitude,
                           const double elevation,
                           const double measure,
                           const Note note ) {

  assert( data ); assert( data->ok );

  if ( data->arguments.format == FORMAT_ASCII ) {

    if ( data->points == 0 ) {
      const Arguments* const arguments = &data->arguments;
      const char* const aggregation =
        arguments->aggregate == AGGREGATE_HOURLY ? "_hourly"
        : arguments->aggregate == AGGREGATE_DAILY ? "_daily"
        : arguments->aggregate == AGGREGATE_MONTHLY ? "_monthly"
        : arguments->aggregate == AGGREGATE_ALL ? "_mean"
        : "";

      printf( "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tELEVATION(m)\t"
              "STATION(-)\t%s%s(%s)\tNOTE\n",
              arguments->variable, aggregation, data->units );
    }

    {
      const int yyyy = yyyymmddhhmmss / 10000000000;
      const int mm = yyyymmddhhmmss / 100000000 % 100;
      const int dd = yyyymmddhhmmss / 1000000 % 100;
      const int hh = yyyymmddhhmmss / 10000 % 100;
      const int min = yyyymmddhhmmss / 100 % 100;
      const int ss = yyyymmddhhmmss % 100;
      printf( "%04d-%02d-%02dT%02d:%02d:%02d-0000\t%f\t%f\t%e\t%lld\t%0.8e\t%s\n",
              yyyy, mm, dd, hh, min, ss,
              longitude, latitude, elevation, id, measure, note );
    }

    data->points += 1;
  } else {
    double dvalues[ 6 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    dvalues[ 0 ] = yyyymmddhhmmss;
    dvalues[ 1 ] = id;
    dvalues[ 2 ] = longitude;
    dvalues[ 3 ] = latitude;
    dvalues[ 4 ] = elevation;
    dvalues[ 5 ] = measure;
    rotate8ByteArrayIfLittleEndian( &dvalues, sizeof dvalues / sizeof *dvalues);
    assert( data->tempFiles[ 0 ] );
    assert( data->tempFiles[ 1 ] );
    assert( data->tempFiles[ 2 ] );
    assert( data->tempFiles[ 3 ] );
    assert( data->tempFiles[ 4 ] );
    assert( data->tempFiles[ 5 ] );
    data->ok =
      fwrite( &dvalues[ 0 ], 8, 1, data->tempFiles[ 0 ] ) == 1 &&
      fwrite( &dvalues[ 1 ], 8, 1, data->tempFiles[ 1 ] ) == 1 &&
      fwrite( &dvalues[ 2 ], 8, 1, data->tempFiles[ 2 ] ) == 1 &&
      fwrite( &dvalues[ 3 ], 8, 1, data->tempFiles[ 3 ] ) == 1 &&
      fwrite( &dvalues[ 4 ], 8, 1, data->tempFiles[ 4 ] ) == 1 &&
      fwrite( &dvalues[ 5 ], 8, 1, data->tempFiles[ 5 ] ) == 1;
    data->points += data->ok;
  }
}




/******************************************************************************
PURPOSE: streamXDRHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
NOTES: Header llooks like this:
Profile 2.0
https://data.pandonia-global-network.org,PandoraSubset
2019-09-10T00:00:00-0000 2019-09-10T23:59:59-0000
# Subset domain <min_lon> <min_lat> <max_lon> <max_lat>:
-130.0 24.0 -65.0 50.0
# Dimensions: variables profiles:
6 2
# Variable names:
timestamp id longitude latitude elevation nitrogen_dioxide_total_vertical_column_amount
# Variable units:
yyyymmddhhmmss - deg deg m mol/m2
# char notes[profiles][80] and
# MSB 64-bit integers points[profiles] and
# IEEE-754 64-bit reals data_1[variables][points_1] ... data_P[variables][points_P]:
<big-endian binary format arrays follow>
******************************************************************************/

static void streamXDRHeader( const Data* const data ) {
  const Arguments* const arguments = &data->arguments;

  /* Append to name to indicate aggregation: */

  const char* const aggregation =
    arguments->aggregate == AGGREGATE_HOURLY ? "_hourly"
    : arguments->aggregate == AGGREGATE_DAILY ? "_daily"
    : arguments->aggregate == AGGREGATE_MONTHLY ? "_monthly"
    : arguments->aggregate == AGGREGATE_ALL ? "_mean"
    : "";

  {
    const long long yyyymmddhhmmss1 = arguments->yyyymmddhhmmss[ 0 ];
    const long long yyyymmddhhmmss2 = arguments->yyyymmddhhmmss[ 1 ];
    const int yyyy1 = yyyymmddhhmmss1 / 10000000000;
    const int mm1   = yyyymmddhhmmss1 / 100000000 % 100;
    const int dd1   = yyyymmddhhmmss1 / 1000000 % 100;
    const int hh1   = yyyymmddhhmmss1 / 10000 % 100;
    const int mm21  = yyyymmddhhmmss1 / 100 % 100;
    const int ss1   = yyyymmddhhmmss1 % 100;
    const int yyyy2 = yyyymmddhhmmss2 / 10000000000;
    const int mm2   = yyyymmddhhmmss2 / 100000000 % 100;
    const int dd2   = yyyymmddhhmmss2 / 1000000 % 100;
    const int hh2   = yyyymmddhhmmss2 / 10000 % 100;
    const int mm22  = yyyymmddhhmmss2 / 100 % 100;
    const int ss2   = yyyymmddhhmmss2 % 100;
    int profileCount = 0;
    const Node* node = data->profileInfoList;

    while ( node ) {
      ++profileCount;
      node = node->next;
    }

    printf( "Profile 2.0\n"
            "%s\n"
            "%04d-%02d-%02dT%02d:%02d:%02d-0000 "
            "%04d-%02d-%02dT%02d:%02d:%02d-0000\n"
            "# Subset domain <min_lon> <min_lat> <max_lon> <max_lat>:\n"
            "%g %g %g %g\n",
            arguments->description,
            yyyy1, mm1, dd1, hh1, mm21, ss1,
            yyyy2, mm2, dd2, hh2, mm22, ss2,
            arguments->bounds[ LONGITUDE ][ MINIMUM ],
            arguments->bounds[ LATITUDE  ][ MINIMUM ],
            arguments->bounds[ LONGITUDE ][ MAXIMUM ],
            arguments->bounds[ LATITUDE  ][ MAXIMUM ] );
    printf( "# Dimensions: variables profiles:\n%d %d\n",
            VARIABLES, profileCount );
    printf( "# Variable names:\n" );
    printf( "timestamp id longitude latitude elevation %s%s\n",
            arguments->variable, aggregation );
    printf( "# Variable units:\nyyyymmddhhmmss - deg deg m %s\n", data->units);
    printf( "# char notes[profiles][80] and\n" );
    printf( "# MSB 64-bit integers points[profiles] and\n" );
    printf( "# IEEE-754 64-bit reals data_1[variables][points_1] ... "
            "data_P[variables][points_P]:\n" );
  }
}



/******************************************************************************
PURPOSE: streamXDRData - Write final XDR binary data to stdout.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* const data  data->ok.
NOTES: XDR format data looks like this:
          char notes[profiles][80] and
          MSB 64-bit integers points[profiles] and
          IEEE-754 64-bit reals data_1[variables][points_1] ...
          data_P[variables][points_P]:
******************************************************************************/

static void streamXDRData( Data* const data ) {

  const Node* node = 0;
  int fileIndex = 0;

  assert( data ); assert( data->ok ); assert( data->points > 0 );

  /* Temp file names are initialized but temp files are closed after writing:*/

  assert( data->tempFileNames[ 0 ][ 0 ] );
  assert( data->tempFiles[     0 ] == 0 );
  assert( data->tempFileNames[ 1 ][ 0 ] );
  assert( data->tempFiles[     1 ] == 0 );
  assert( data->tempFileNames[ 2 ][ 0 ] );
  assert( data->tempFiles[     2 ] == 0 );
  assert( data->tempFileNames[ 3 ][ 0 ] );
  assert( data->tempFiles[     3 ] == 0 );
  assert( data->tempFileNames[ 4 ][ 0 ] );
  assert( data->tempFiles[     4 ] == 0 );
  assert( data->tempFileNames[ 5 ][ 0 ] );
  assert( data->tempFiles[     5 ] == 0 );

  /* Write notes: */

  for ( node = data->profileInfoList; data->ok && node; node = node->next) {
    const ProfileInfo* const profileInfo = (const ProfileInfo*) node->data;
    assert( profileInfo ); assert( profileInfo->note[ 0 ] );
    assert( strlen( profileInfo->note ) < 80 );
    data->ok = printf( "%-79s\n", profileInfo->note ) == 80;
  }

  /* Write profile point counts: */

  for ( node = data->profileInfoList; data->ok && node; node = node->next) {
    const ProfileInfo* const profileInfo = (const ProfileInfo*) node->data;
    long long count = profileInfo->points;
    rotate8ByteArrayIfLittleEndian( &count, 1 );
    data->ok = fwrite( &count, 8, 1, stdout ) == 1;
  }

  /* Write each profile's variable temp file content to stdout: */

  openTempFiles( data );

  for ( node = data->profileInfoList; data->ok && node; node = node->next) {
    const ProfileInfo* const profileInfo = (const ProfileInfo*) node->data;
    const size_t bytes = profileInfo->points * 8;

    for ( fileIndex = 0; data->ok && fileIndex < VARIABLES; ++fileIndex ) {
      data->ok = streamBytes( data->tempFiles[ fileIndex ], bytes );
    }
  }

  closeTempFiles( data );

  if ( ! data->ok ) {
    fprintf( stderr, "\nFailed to stream all subset data to stdout.\n" );
  }
}




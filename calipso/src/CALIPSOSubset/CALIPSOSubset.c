/******************************************************************************
PURPOSE: CALIPSOSubset.c - Extract a lon-lat subset of data from a list of
         CALIPSO HDF files and write it to stdout as XDR binary format.

NOTES:   Uses HDF libraries and libs they depend on.

         Compile with debugging enabled:
         gcc -Wall -g -o CALIPSOSubset \
             CALIPSOSubset.c Utilities.c ReadData.c \
             -L../../../lib/$platform \
             -lhdfeos -lmfhdf -ldf -lsz -lz -ljpeg -lm

         Compile optimized:
         gcc -Wall -DNDEBUG -O -o CALIPSOSubset \
             CALIPSOSubset.c Utilities.c ReadData.c \
             -L../../../lib/$platform \
             -lhdfeos -lmfhdf -ldf -lsz -lz -ljpeg -lm

         Usage:
         CALIPSOSubset -files <listfile> \
                     -tmpdir <temp_directory> \
                     -desc "description text" \
                     -timestamp <yyyymmddhh> -hours <count> \
                     -variable <name> \
                     -domain <minimum_longitude> <minimum_latitude> \
                             <maximum_longitude> <maximum_latitude> \
                     [-elevation <minimum_elevation> <maximum_elevation>] \
                     [-minimumCAD <value>] \
                     [-maximumUncertainty <value>]

          Example:
          ../../../bin/$platform/CALIPSOSubset \
          -files testdata/files \
          -tmpdir testdata \
          -desc \
 https://eosweb.larc.nasa.gov/project/calipso/calipso_table,CALIPSOSubset \
          -timestamp 2006080600 -hours 24 \
          -variable Extinction_Coefficient_532 \
          -domain -110 35 -75 36 \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 15-line ASCII header followed by binary 64-bit big-endian arrays:

CALIPSO 1.0
https://eosweb.larc.nasa.gov/project/calipso/calipso_table,CALIPSOSubset
2006-08-06T00:00:00-0000
# Dimensions: variables timesteps profiles:
6 24 4
# Variable names:
Profile_UTC_Time Longitude Latitude Elevation Extinction_Coefficient_532 Thickness
# Variable units:
yyyymmdd.f deg deg m - m
# Domain: <min_lon> <min_lat> <max_lon> <max_lat>
-110 35 -75 36
# MSB 64-bit integers (yyyydddhhmm) profile_timestamps[profiles] and
# IEEE-754 64-bit reals profile_bounds[profiles][2=<lon,lat>][2=<min,max>] \
   and
# MSB 64-bit integers profile_dimensions[profiles][2=<points,levels>] and
# IEEE-754 64-bit reals profile_data_1[variables][points_1][levels] ... \
   profile_data_S[variables][points_S][levels]:
<big-endian binary format arrays:>
20061860735
20061860914
20061861815
20061861954
-8.3809753417968750e+01
-8.3527038574218750e+01
3.5009464263916016e+01
3.5990005493164062e+01
-1.0852429962158203e+02
-1.0825477600097656e+02
3.5032341003417969e+01
3.5969318389892578e+01
-8.1201416015625000e+01
-8.0931808471679688e+01
3.5020053863525391e+01
3.5956905364990234e+01
-1.0592688751220703e+02
-1.0564466094970703e+02
3.5001144409179688e+01
3.5982353210449219e+01
23
8
22
8
22
8
23
8
2.0060705320214204e+07
2.0060705320222814e+07
...
-9.9990000000000000e+03
-9.9990000000000000e+03

Profile_UTC_Time Longitude Latitude are only for ground points -
i.e., implicitly dimensioned with levels = 1. UGLY.

HISTORY: 2007-02-01 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h>   /* For assert(). */
#include <stdio.h>    /* For FILE, printf(), snprintf(). */
#include <string.h>   /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <stdlib.h>   /* For malloc(), free(). */
#include <unistd.h>   /* For unlink(), getpid() */

#include "Utilities.h"/* For MISSING_VALUE, IN_RANGE(), LONGITUDE, Bounds. */
#include "ReadFile.h" /* For openFile(), readFileBounds(). */
#include "ReadData.h" /* For readCALIPSOData(). */

/*================================= MACROS =================================*/

/* Name of temporary file created in -tmpdir will have PID appended: */

#define TEMP_FILE_NAME "junk_CALIPSOSubset"

/*================================== TYPES ==================================*/

enum { CALIPSO_L1_AGGREGATION_WINDOW = 15 }; /* 333m to 5km ground points. */
enum { CALIPSO_L1_AGGREGATION_TARGET_LEVELS = 100 }; /* Curtain height pixels.*/

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;    /* File containing list of CALIPSO files to read. */
  const char* tmpdir;      /* Name of directory to write temp files. */
  const char* description; /* User-supplied description. */
  const char* variable;    /* Name of variable to read. */
  Bounds      domain; /* Subset domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  double      elevationRange[ 2 ]; /* Subset elevation[MINIMUM,MAXIMUM]. */
  double      minimumCAD;  /* Minimum cloud/aerosol discrimination score. */
                           /* E.g., 20 accepts score in range [20, 100]. */
  double      maximumUncertainty; /* Maximum uncertainty acceptable. */
                           /* E.g., 99 accepts absolute uncertainty <= 99. */
                           /* Units are the same as the data variable. */
  int         yyyymmddhh;  /* First timestamp of subset. */
  int         hours;       /* Number of hours in subset. */
} Arguments;

/* Data type: */

typedef struct {
  Arguments  arguments;    /* User-supplied (command-line) arguments. */
  char units[ 80 ];        /* Units of variable. */
  FileName   tempFileName; /* Name of temp file of output subset data. */
  FILE*      tempFile;     /* Temp file of output subset data. */
  long long* yyyydddhhmm;  /* Timestamp per output subset scan. */
  long long* pointsAndLevels; /* # of ground points and levels per subset scan*/
                              /* pointsAndLevels[scan][2=ground,level]. */
  Bounds*    bounds; /* bounds[scan][LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  int fileType;            /* CALIPSO_L1 ... CALIPSO_L2_VFM. */
  int hasThickness;        /* Does the variable have thickness data? */
  int        scans;        /* Number of output subset scans. */
  int        ok;           /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocate( Data* data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static void readData( Data* const data );

static char* readListFileAndAllocateScanMetadata( Data* const data );

static void allocateScanMetadata( const char* const listFileContent,
                                  Data* const data );

static long long dataFileTimestamp( const char* const fileName );

static void readFileInfo( Data* const data,
                          const char* const fileName,
                          int* const file,
                          long long* const yyyydddhhmm,
                          size_t* const points,
                          size_t* const levels,
                          size_t* const size,
                          int* const changedDimensions );

static void writeSubsetData( Data* const data,
                             const size_t points,
                             const size_t levels,
                             double timestamps[],
                             double longitudes[],
                             double latitudes[],
                             double elevations[],
                             double thicknesses[],
                             double values[] );

static void streamData( Data* const data );

static void streamHeader( const Data* const data );

static void streamSwathTimestamps( Data* const data );

static void streamSwathBounds( Data* const data );

static void streamSwathPointsAndLevels( Data* const data );

static void streamTempFile( Data* const data );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of CALIPSO files and write
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
    readData( &data ); /* Read subset of CALIPSO files and write temp files. */

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

  if ( data->tempFile ) {
    fclose( data->tempFile );
    data->tempFile = 0;
  }

  if ( data->tempFileName[ 0 ] ) {
    unlink( data->tempFileName );
    memset( data->tempFileName, 0, sizeof data->tempFileName );
  }

  if ( data->yyyydddhhmm ) {
    free( data->yyyydddhhmm );
    data->yyyydddhhmm = 0;
  }

  if ( data->pointsAndLevels ) {
    free( data->pointsAndLevels );
    data->pointsAndLevels = 0;
  }

  if ( data->bounds ) {
    free( data->bounds );
    data->bounds = 0;
  }

  memset( data, 0, sizeof *data );
}



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* name  Name of program.
******************************************************************************/

static void printUsage( const char* programName ) {
  assert( programName );
  fprintf( stderr, "\a\n\n%s - Read a set of CALIPSO files and extract swath\n",
           programName );
  fprintf( stderr, "data subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable(s).\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "  -files <file> \\\n" );
  fprintf( stderr, "  -tmpdir <temp_directory> \\\n" );
  fprintf( stderr, "  -desc description \\\n" );
  fprintf( stderr, "  -timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, "  -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, "  [-elevation <minimum_elevation> <maximum_elevation>] \\\n");
  fprintf( stderr, "  [-minimumCAD <value>] \\\n" );
  fprintf( stderr, "  [-maximumUncertainty <value>]\n\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies where to write temporary files.\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example #1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files testdata/files2.txt \\\n" );
  fprintf( stderr, "-tmpdir testdata \\\n" );
  fprintf( stderr, "-desc "
           "https://eosweb.larc.nasa.gov/project/calipso/calipso_table,"
           "CALIPSOSubset \\\n");
  fprintf( stderr, "-timestamp 2006070500 -hours 24 \\\n" );
  fprintf( stderr, "-variable Extinction_Coefficient_532 \\\n" );
  fprintf( stderr, "-domain -110 35 -75 36 -elevation 0 16000 "
           "> subset.xdr\n\n" );
  fprintf( stderr, "AOD over part of US on July 5, 2014 not more than 16km "
          "above the surface.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n" );
  fprintf( stderr, "CALIPSO 1.0\n" );
  fprintf( stderr, "https://eosweb.larc.nasa.gov/project/calipso/calipso_table,"
           "CALIPSOSubset\n");
  fprintf( stderr, "2014-07-05T00:00:00-0000\n" );
  fprintf( stderr, "# Dimensions: variables timesteps profiles:\n" );
  fprintf( stderr, "6 24 4\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Profile_UTC_Time Longitude Latitude Elevation "
                     "Extinction_Coefficient_532 Thickness\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmdd.f deg deg m - m\n" );
  fprintf( stderr, "# Domain: <min_lon> <min_lat> <max_lon> <max_lat>\n" );
  fprintf( stderr, "-110 35 -75 36\n" );
  fprintf( stderr, "# MSB 64-bit integers (yyyydddhhmm)" );
  fprintf( stderr, " profile_timestamps[profiles] and\n" );
  fprintf( stderr, "IEEE-754 64-bit reals profile_bounds[profiles]"
                   "[2=<lon,lat>][2=<min,max>] and\n" );
  fprintf( stderr, "# MSB 64-bit integers"
                   " profile_dimensions[profiles][2=<points,levels>] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals" );
  fprintf( stderr, " profile_data_1[variables][points_1][levels] ..." );
  fprintf( stderr, " profile_data_S[variables][points_S][levels]:\n" );
  fprintf( stderr, "<binary data arrays here>\n\n");
  fprintf( stderr, "Note: Profile_UTC_Time Longitude Latitude are" );
  fprintf( stderr, " only for ground points -\n" );
  fprintf( stderr, "i.e., implicitly dimensioned with levels = 1. UGLY.\n\n\n");
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
  arguments->elevationRange[ MINIMUM ] = -500.0;
  arguments->elevationRange[ MAXIMUM ] = 1e5;
  arguments->minimumCAD = 20.0;
  arguments->maximumUncertainty = 99.0;

  result = IN7( argc, 18, 20, 21, 22, 23, 25 );

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
    } else if ( ! strcmp( argv[ arg ], "-elevation" ) ) {
      double value = 0.0;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->elevationRange[ MINIMUM ] = value;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->elevationRange[ MAXIMUM ] = value;
      result =
        AND3( isValidElevation( arguments->elevationRange[ MINIMUM ] ),
              isValidElevation( arguments->elevationRange[ MAXIMUM ] ),
              arguments->elevationRange[ MINIMUM ] <=
                arguments->elevationRange[ MAXIMUM ] );
    } else if ( ! strcmp( argv[ arg ], "-minimumCAD" ) ) {
      double value = 0.0;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->minimumCAD = value;
      result = IN_RANGE( arguments->minimumCAD, 0.0, 100.0 );
    } else if ( ! strcmp( argv[ arg ], "-maximumUncertainty" ) ) {
      double value = 0.0;
      ++arg;
      value = atof( argv[ arg ] );
      arguments->maximumUncertainty = value;
      result = IN_RANGE( arguments->maximumUncertainty, 0.0, 99.0 );
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
PURPOSE: readData - Read swath data from each listed CALIPSO file and
         write the lon-lat subset of data to the temporary file.
INPUTS:  Data* data  Data description to read.
OUTPUTS: Data* data  data->yyyydddhhmm, points, scans, ok, tempFile = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  char* listFileContent = readListFileAndAllocateScanMetadata( data );
  int wroteSomeData    = 0;
  char* fileName       = 0;
  char* end            = 0;
  size_t points        = 0;
  size_t levels        = 0;
  size_t size          = 0;
  double* buffer       = 0;
  double* timestamps   = 0;
  double* longitudes   = 0;
  double* latitudes    = 0;
  double* elevations   = 0;
  double* values       = 0;
  double* thicknesses  = 0;
  data->ok = 0;

  /* Get each line of list file. It is the CALIPSO data file to read: */

  for ( fileName = strtok_r( listFileContent, "\n", &end );
        fileName;
        fileName = strtok_r( 0, "\n", &end ) ) {
    int file = 0;
    long long yyyydddhhmm = 0;
    int changedDimensions = 0;

    readFileInfo( data, fileName, &file, &yyyydddhhmm, &points, &levels,
                  &size, &changedDimensions );

    if ( data->ok ) {

      if ( changedDimensions ) {
        const size_t groundVariables = 3; /* time,lon,lat. */
        const size_t levelVariables = 2 + data->hasThickness; /*elv,aod+thick*/
        const size_t dataSize =
          groundVariables * points + levelVariables * points * levels;
        const size_t bytes = dataSize * sizeof (double);

        if ( buffer ) {
          free( buffer );
          buffer = 0;
        }

        buffer = malloc( bytes );
        data->ok = buffer != 0;

        if ( buffer ) {
          timestamps  = buffer;
          longitudes  = timestamps + points;
          latitudes   = longitudes + points;
          elevations  = latitudes  + points;
          values      = elevations + size;
          thicknesses = data->hasThickness ? values + size : 0;
        } else {
          fprintf( stderr,
                   "\nCan't allocate %lu bytes "
                   "to complete the requested action.\n", bytes );
        }
      }

      if ( data->ok ) {
        data->ok =
          readCALIPSOData( file, data->fileType, arguments->variable,
                           points, levels,
                           arguments->minimumCAD,
                           arguments->maximumUncertainty,
                           data->units,
                           timestamps, longitudes, latitudes,
                           elevations, thicknesses, values );
      }

      closeFile( file ), file = -1;

      if ( data->ok ) {
        size_t subsetPoints = 0;
        size_t subsetLevels = 0;
        data->ok =
          compactPointsInSubset( (const double (*)[2]) arguments->domain,
                                 arguments->elevationRange[ MINIMUM ],
                                 arguments->elevationRange[ MAXIMUM ],
                                 points, levels,
                                 timestamps, longitudes, latitudes,
                                 elevations, values, thicknesses,
                                 &subsetPoints, &subsetLevels );

        DEBUG( fprintf( stderr,
                        "ok = %d, subsetPoints = %lu, subsetLevels = %lu\n",
                        data->ok, subsetPoints, subsetLevels ); )

        if ( data->ok ) {

          if ( data->fileType == CALIPSO_L1 ) {
            assert( thicknesses == 0 ); /* L1 data does not use thicknesses. */
            aggregateCALIPSOData( subsetPoints, subsetLevels,
                                  CALIPSO_L1_AGGREGATION_WINDOW,
                                  CALIPSO_L1_AGGREGATION_TARGET_LEVELS,
                                  timestamps, longitudes, latitudes,
                                  elevations, values,
                                  &subsetPoints, &subsetLevels );
          }

          {
            const size_t scans2 = data->scans * 2;
            data->yyyydddhhmm[ data->scans ] = yyyydddhhmm;
            data->pointsAndLevels[ scans2     ] = subsetPoints;
            data->pointsAndLevels[ scans2 + 1 ] = subsetLevels;
            computeBounds( subsetPoints, longitudes, latitudes,
                           data->bounds[ data->scans ] );
            DEBUG( fprintf( stderr, "output scan %d: %lld %lld x %lld\n",
                            data->scans, data->yyyydddhhmm[ data->scans ],
                            data->pointsAndLevels[ scans2 ],
                            data->pointsAndLevels[ scans2 + 1 ] ); )
            data->scans += 1;
          }

          writeSubsetData( data, subsetPoints, subsetLevels,
                           timestamps, longitudes, latitudes,
                           elevations, values, thicknesses );

          if ( data->ok ) {
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
PURPOSE: readListFileAndAllocateScanMetadata - Read list file and
         return its contents as a string and allocate timestamps, bounds and
         pointsAndLevels arrays with length equal to lines in the list file.
INPUTS:  Data* const data  data->arguments.listFile.
OUTPUTS: Data* const data  data->ok, yyyydddhhmm, bounds, pointsAndLevels
                                     allocated.
RETURNS: char* Allocated string containing lines of list file or 0 if failed.
******************************************************************************/

static char* readListFileAndAllocateScanMetadata( Data* const data ) {
  char* result = 0;
  assert( data ); assert( data->arguments.listFile );
  assert( data->arguments.listFile[ 0 ] );
  assert( data->yyyydddhhmm == 0 );
  assert( data->pointsAndLevels == 0 );

  {
    size_t length = 0;
    result = readFile( data->arguments.listFile, &length );
    data->ok = result != 0;

    if ( result ) {
      allocateScanMetadata( result, data );

      if ( ! data->ok ) {
        free( result );
        result = 0;
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: allocateScanMetadata - Allocate arrays for timestamps, bounds
         and pointsAndLevels per scan.
INPUTS:  const char* const listFileContent  String content of list file.
         Data* const data                   Data to allocate.
OUTPUTS: Data* const data  data->ok, data->yyyydddhhmm, bounds, pointsAndLevels
                           allocated.
******************************************************************************/

static void allocateScanMetadata( const char* const listFileContent,
                                         Data* const data ) {
  assert( listFileContent );
  assert( data );
  assert( data->arguments.listFile );
  assert( data->yyyydddhhmm == 0 );
  assert( data->pointsAndLevels == 0 );

  {
    const size_t lines = linesInString( listFileContent );
    data->ok = lines != 0;

    if ( ! lines ) {
      fprintf(stderr, "\nInvalid list file '%s'.\n", data->arguments.listFile);
    } else {
      size_t bytes = lines * sizeof (long long);
      data->yyyydddhhmm = malloc( bytes );
      data->pointsAndLevels = data->yyyydddhhmm ? malloc( bytes * 2 ) : 0;
      data->ok = data->pointsAndLevels != 0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nCan't allocate %lu bytes to complete the "
                 "requested action.\n", bytes * 3 );
      } else {
        bytes = lines * sizeof (Bounds);
        data->bounds = malloc( bytes );
        data->ok = data->bounds != 0;

        if ( ! data->ok ) {
          fprintf( stderr, "\nCan't allocate %lu bytes to complete the "
                   "requested action.\n", bytes );
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: dataFileTimestamp - Timestamp of CALIPSO file.
INPUTS:  const char* const fileName   Name of CALIPSO file.
RETURNS: long long yyyydddhhmm of file or 0 if failed and message on stderr.
NOTES:   File names look like: CAL_LID_L1-Prov-V1-10.2006-07-04T23-21-01ZN.hdf
******************************************************************************/

static long long dataFileTimestamp( const char* const fileName ) {
  long long result = 0;
  const size_t length = fileName ? strlen( fileName ) : 0;
  assert( fileName ), assert( *fileName );

  if ( length >= 14 ) {
    const size_t start = length - 25;
    const char* const fileNameTimestamp = fileName + start;
    const long yyyy = atoi( fileNameTimestamp );
    const long mo   = atoi( fileNameTimestamp + 5 );
    const long dd   = atoi( fileNameTimestamp + 8 );
    const long hh   = atoi( fileNameTimestamp + 11 );
    const long mm   = atoi( fileNameTimestamp + 14 );
    const long ss   = atoi( fileNameTimestamp + 17 );

    if ( AND6( IN_RANGE( yyyy, 1900, 3000 ),
               IN_RANGE( mo, 1, 12 ),
               IN_RANGE( dd, 1, 31 ),
               IN_RANGE( hh, 0, 23 ),
               IN_RANGE( mm, 0, 59 ),
               IN_RANGE( ss, 0, 59 ) ) ) {
      long long yyyymmddhhmm = yyyy;
      yyyymmddhhmm *= 100;
      yyyymmddhhmm += mo;
      yyyymmddhhmm *= 100;
      yyyymmddhhmm += dd;
      yyyymmddhhmm *= 100;
      yyyymmddhhmm += hh;
      yyyymmddhhmm *= 100;
      yyyymmddhhmm += mm;

      if ( isValidYYYYMMDDHHMM( yyyymmddhhmm ) ) {
        result = convertTimestamp( yyyymmddhhmm );
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
         size_t* const points            points of data in previous file or 0.
         size_t* const levels          levels of data in previous file or 0.
OUTPUTS: Data* const data              data->ok, data->hasThickness.
         int* const file               HDF file id of data file.
         long long* const yyyydddhhmm  Timestamp of file (if in subset range).
         size_t* const points          Ground points of data in file.
         size_t* const levels          Vetical levels of data in file.
         size_t* const size            points * levels.
         int* const changedDimensions  1 if points or levels changed.
******************************************************************************/

static void readFileInfo( Data* const data,
                          const char* const fileName,
                          int* const file,
                          long long* const yyyydddhhmm,
                          size_t* const points,
                          size_t* const levels,
                          size_t* const size,
                          int* const changedDimensions ) {

  assert( data );
  assert( isValidBounds( (const double (*)[2]) data->arguments.domain ) );
  assert( fileName ); assert( *fileName );  assert( file );
  assert( yyyydddhhmm ); assert( points ); assert( levels ); assert( size );
  assert( changedDimensions );

  *file = -1;
  *yyyydddhhmm = dataFileTimestamp( fileName );
  data->ok = *yyyydddhhmm != 0;
  data->hasThickness = 0;
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
            size_t points2 = 0;
            size_t levels2 = 0;
            data->ok =
              readCALIPSOVariableDimensions( *file, arguments->variable,
                                             &points2, &levels2 );

            if ( data->ok ) {
              data->fileType = typeOfCALIPSOFile( fileName );
              data->hasThickness =
              AND2( IS_LAYERED( data->fileType ), levels2 > 1 );

              if ( OR2( points2 != *points, levels2 != *levels ) ) {
                *points = points2;
                *levels = levels2;
                *size = *points * *levels;
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
PURPOSE: writeSubsetData - Write subset of data to temp file.
INPUTS:  Data* const data                     Data description.
         const size_t points                  Points in subset domain.
         double* const timestamps[  points ]      Timestamps (yyyymmdd.f) read.
         double* const longitudes[  points ]           Longitudes read.
         double* const latitudes[   points ]           Latitudes read.
         double* const elevations[  points * levels ]  Elevations read.
         double* const thicknesses[ points * levels ] If non-0 Thicknesses read
         double* const values[      points * levels ]  Values read.
OUTPUTS: Data* data  data->ok, ruined arrays.
******************************************************************************/

static void writeSubsetData( Data* const data,
                             const size_t points,
                             const size_t levels,
                             double timestamps[],
                             double longitudes[],
                             double latitudes[],
                             double elevations[],
                             double values[],
                             double thicknesses[] ) {

  assert( data ); assert( data->ok );
  assert( points ); assert( levels );
  assert( timestamps ); assert( longitudes ); assert( latitudes );
  assert( elevations ); assert( values );

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
    const size_t pointsTimesLevels = points * levels;
    rotate8ByteArrayIfLittleEndian( points, timestamps );
    rotate8ByteArrayIfLittleEndian( points, longitudes );
    rotate8ByteArrayIfLittleEndian( points, latitudes );
    rotate8ByteArrayIfLittleEndian( pointsTimesLevels, elevations );
    rotate8ByteArrayIfLittleEndian( pointsTimesLevels, values );

    if ( thicknesses ) {
      rotate8ByteArrayIfLittleEndian( pointsTimesLevels, thicknesses );
    }

    data->ok =
      fwrite( timestamps, sizeof *timestamps, points, data->tempFile) == points;

    if ( data->ok ) {
      data->ok =
        fwrite( longitudes, sizeof *longitudes, points, data->tempFile )
        == points;

      if ( data->ok ) {
        data->ok =
          fwrite( latitudes, sizeof *latitudes, points, data->tempFile )
          == points;

        if ( data->ok ) {
          data->ok =
          fwrite( elevations, sizeof *elevations, pointsTimesLevels,
                  data->tempFile ) == pointsTimesLevels;

          if ( data->ok ) {
            data->ok =
              fwrite( values, sizeof *values, pointsTimesLevels, data->tempFile)
              == pointsTimesLevels;

            if ( data->ok ) {

              if ( thicknesses ) {
                data->ok =
                  fwrite( thicknesses, sizeof *thicknesses, pointsTimesLevels,
                          data->tempFile ) == pointsTimesLevels;
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
    streamSwathBounds( data );

    if ( data->ok ) {
      streamSwathPointsAndLevels( data );

      if ( data->ok ) {
        streamTempFile( data );
      }
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
    const int hasThickness = data->hasThickness;
    const int variables = 5 + hasThickness;
    const int scans = data->scans;
    const int yyyymmddhh = arguments->yyyymmddhh;
    const int yyyy       = yyyymmddhh / 1000000;
    const int mm         = yyyymmddhh / 10000 % 100;
    const int dd         = yyyymmddhh / 100 % 100;
    const int hh         = yyyymmddhh % 100;
    printf( "CALIPSO 1.0\n%s\n%04d-%02d-%02dT%02d:00:00-0000\n",
            arguments->description, yyyy, mm, dd, hh );
    printf( "# Dimensions: variables timesteps profiles:\n%d %d %d\n",
            variables, arguments->hours, scans );
    printf( "# Variable names:\n" );
    printf( "Profile_UTC_Time Longitude Latitude Elevation %s%s\n",
            arguments->variable, hasThickness ? " Thickness" : "" );
    printf( "# Variable units:\nyyyymmdd.f deg deg m %s%s\n",
            data->units, hasThickness ? " m" : "" );
    printf( "# Domain: <min_lon> <min_lat> <max_lon> <max_lat>\n%g %g %g %g\n",
            arguments->domain[ LONGITUDE ][ MINIMUM ],
            arguments->domain[ LATITUDE  ][ MINIMUM ],
            arguments->domain[ LONGITUDE ][ MAXIMUM ],
            arguments->domain[ LATITUDE  ][ MAXIMUM ] );
    printf( "# MSB 64-bit integers (yyyydddhhmm)" );
    printf( " profile_timestamps[profiles] and\n" );
    printf( "# IEEE-754 64-bit reals profile_bounds[profiles]"
            "[2=<lon,lat>][2=<min,max>] and\n" );
    printf( "# MSB 64-bit integers profile_dimensions[profiles]"
            "[2=<points,levels>] and\n" );
    printf( "# IEEE-754 64-bit reals" );
    printf( " profile_data_1[variables][points_1][levels] ..." );
    printf( " profile_data_S[variables][points_S][levels]:\n" );
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

  rotate8ByteArrayIfLittleEndian( data->scans, data->yyyydddhhmm );
  data->ok =
    fwrite( data->yyyydddhhmm, sizeof data->yyyydddhhmm[ 0 ], data->scans,
            stdout ) == data->scans;
  rotate8ByteArrayIfLittleEndian( data->scans, data->yyyydddhhmm );

  if ( ! data->ok ) {
    fprintf( stderr, "\a\nFailed to stream subset swath timestamps.\n" );
  }
}



/******************************************************************************
 PURPOSE: streamSwathBounds - Write IEEE-754 64-bit real swath subset bounds.
 INPUTS:  Data* const data  Data to stream.
 OUTPUTS: Data* const data  data->ok.
 ******************************************************************************/

static void streamSwathBounds( Data* const data ) {

  assert( data ); assert( data->scans > 0 );
  assert( isValidBounds( (const double (*)[2]) data->bounds[ 0 ] ) );
  assert( isValidBounds( (const double (*)[2]) data->bounds[ data->scans - 1]));

  rotate8ByteArrayIfLittleEndian( data->scans * 4, data->bounds );
  data->ok =
    fwrite( data->bounds, sizeof data->bounds[ 0 ],
            data->scans, stdout ) == data->scans;
  rotate8ByteArrayIfLittleEndian( data->scans * 4, data->bounds );

  if ( ! data->ok ) {
    fprintf( stderr, "\a\nFailed to stream subset swath bounds.\n" );
  }
}



/******************************************************************************
 PURPOSE: streamSwathPointsAndLevels - Write MSB 64-bit integer swath subset point counts
 INPUTS:  Data* const data  Data to stream.
 OUTPUTS: Data* const data  data->ok.
 ******************************************************************************/

static void streamSwathPointsAndLevels( Data* const data ) {

  assert( data ); assert( data->scans > 0 );
  assert( data->pointsAndLevels );
  assert( data->pointsAndLevels[ 0 ] > 0 );
  assert( data->pointsAndLevels[ data->scans * 2 - 1 ] > 0 );

  rotate8ByteArrayIfLittleEndian( data->scans * 2, data->pointsAndLevels );
  data->ok =
    fwrite( data->pointsAndLevels, sizeof data->pointsAndLevels[ 0 ],
            data->scans * 2, stdout ) == data->scans * 2;
  rotate8ByteArrayIfLittleEndian( data->scans * 2, data->pointsAndLevels );

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




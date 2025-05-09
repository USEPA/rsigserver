/******************************************************************************
PURPOSE: PointSubset.c - Extract a time and lon-lat subset of data from a list
         of Point xdr files and write it to stdout as XDR binary format.

NOTES:   The input files are created by the ../PurpleAirSubset program.

         Compile:
         gcc -Wall -g -o PointSubset PointSubset.c Utilities.c \
                   -L../../../lib/$platform \
                   -lm -lc

         Usage:
         PointSubset \
           -files <listfile> \
           -bounds <minimum_longitude> <minimum_latitude> \
                   <maximum_longitude> <maximum_latitude> \
           [-timerange <yyyymmddhhmmss> <yyyymmddhhmmss> ] \
           [-sensor sensor_id] \

          Example:
          ../../../bin/$platform/PointSubset \
          -files testdata/file_list \
          -tmpdir testdata \
          -bounds -124 49 -123 50 \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 10-line ASCII header followed by binary 64-bit big-endian arrays:

Point 1.0
https://api.purpleair.com,PurpleAirSubset,PointSubset
2020-12-02T00:00:00-0000 2020-12-02T23:59:59-0000
# Dimensions: variables points
6 459
# Variable names:
timestamp longitude latitude elevation id pm25
# Variable units:
yyyymmddhhmmss deg deg m - ug/m3
# char notes[points][80] and
# IEEE-754 64-bit reals data[variables][points]:
<big-endian binary format array>

HISTORY: 2020-12-29 plessel.todd@epa.gov
STATUS:  unreviewed untested
******************************************************************************/


/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, printf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <stdlib.h>    /* For malloc(), free(), system(). */
#include <limits.h>    /* For INT_MAX. */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For LONGITUDE, Bounds. */

/*================================= MACROS =================================*/

/* Name of temp per-variable files created in -tmpdir with PID appended: */

#define TEMP_FILE_NAME "junk_PointSubset"

/*================================== TYPES ==================================*/

enum { VARIABLES = 7 };/*timestamp,longitude,latitude,elevation,id,count,pm25*/
enum { TEMP_FILES = 1 + VARIABLES }; /* notes + variables. */
enum { HEADER_LINES = 11 }; /* Input file header lines. */
static const char* const tempFileNames[ TEMP_FILES ] = {
  "notes", "timestamp", "longitude", "latitude", "elevation", "id",
  "count", "data"
};

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

enum { NOTE_LENGTH = 79 };
typedef char Note[ NOTE_LENGTH + 1 ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile; /* File containing list of input files to read.*/
  const char* tmpdir;   /* Name of directory to write temp files. */
  Bounds      bounds;   /* Subset bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]*/
  long long   yyyymmddhhmmss[ 2 ]; /* Beginning/ending timestamp of subset. */
  int         sensor;   /* 0 = all. > 0 for specific sensor id. */
} Arguments;

/* Data type: */

typedef struct {
  Arguments arguments;                   /* Command-line arguments. */
  FileName  tempFileNames[ TEMP_FILES ]; /* Pathed names of temp files. */
  FILE*     tempFiles[ TEMP_FILES ];     /* Temp files of output subset data.*/
  char*     buffer;                      /* Buffer for reading data. */
  size_t    bufferSize;                  /* Size, in bytes, of buffer. */
  Note      header[ HEADER_LINES ];      /* Input file header lines. */
  long long yyyymmddhhmmss[ 2 ];         /* First and last file timestamps. */
  int variables;                         /* 6 or 7 if counts. */
  size_t    points;                      /* Number of data points in subset. */
  int       ok;                          /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static void removeTempFiles( Data* const data );

static void createTempFiles( Data* const data );

static void closeTempFiles( Data* const data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static void checkAndReallocateBuffer( const size_t bytes, Data* const data );

static void readData( Data* const data );

static int readHeader( FILE* inputFile,
                       int* variables, size_t* points,
                       long long yyyymmddhhmmss[], Note header[] );

static int parseTimestamps( const char* line, long long yyyymmddhhmmss[] );

static void extractSubset( FILE* inputFile,
                           const size_t points, Data* data );

static void streamData( Data* const data );

static void streamHeader( Data* const data );


/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of Point files and write
         it to stdout in XDR format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Data data;
  memset( &data, 0, sizeof data );
  data.ok = parseArguments( argc, argv, &data.arguments );

  if ( ! data.ok ) {
    printUsage( argv[ 0 ] );
  } else {
    createTempFiles( &data );

    if ( data.ok ) {
      const size_t bytes = 2 * 1024 * 1024; /* 2MB should avoid reallocations*/
      checkAndReallocateBuffer( bytes, &data );

      if ( data.ok ) {
        readData( &data ); /* Read input data and write subset to temp files. */

        if ( data.ok ) {
          streamData( &data ); /* Write header & temp files to stdout. */
        }
      }
    }

    if ( ! data.ok ) {
      fprintf( stderr, "\n%s: No points were in the subset.\n", argv[ 0 ] );
    }
  }

  if ( data.buffer ) {
    free( data.buffer ), data.buffer = 0;
    data.bufferSize = 0;
  }

  removeTempFiles( &data );
  data.ok = data.ok && data.points > 0;
  return ! data.ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* name  Name of program.
******************************************************************************/

static void printUsage( const char* name ) {
  assert( name ); assert( *name );
  fprintf( stderr,
           "\n%s - Extract a subset of data from a time-sorted list of\n"
          "Point xdr files and write it to stdout in XDR binary format.\n",
           name );
  fprintf( stderr, "Data is subsetted by lon-lat rectangle.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "  -files <listfile> \\\n" );
  fprintf( stderr, "  -tmpdir <temp_directory> \\\n" );
  fprintf( stderr, "  -bounds <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, " [-timerange <yyyymmddhhmmss> <yyyymmddhhmmss> ]] \\\n" );
  fprintf( stderr, " [-sensor sensor_id] (subset to specific sensor id)\n\n" );
  fprintf( stderr, "-tmpdir specifies a directory were temp files are " );
  fprintf ( stderr, "written.\nIt should have enough disk space (1TB).\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-files file_list \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-bounds -124 49 -123 50 \\\n" );
  fprintf( stderr, "> subset.xdr\n\n" );
  fprintf( stderr, "Daily corrected PM2.5 over BC on December 2, 2020.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Point 1.0\n" );
  fprintf( stderr, "https://api.purpleair.com,PointSubset\n" );
  fprintf( stderr, "2020-12-02T00:00:00-0000 2020-12-02T23:59:59-0000\n" );
  fprintf( stderr, "# Dimensions: variables points\n" );
  fprintf( stderr, "6 20\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "timestamp longitude latitude elevation id pm25_corrected\n");
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss deg deg m - ug/m3\n" );
  fprintf( stderr, "# char notes[points][80] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals data[variables][points]:\n" );
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
  int result = argc >= 5;
  static const int one_int_max[ 2 ] = { 1, INT_MAX };
  static Option options[] = {
    { "-files",  1, FILE_TYPE,      1, 0, 0, 0, 0 },
    { "-tmpdir", 1, DIRECTORY_TYPE, 1, 0, 0, 0, 0 },
    { "-bounds", 0, BOUNDS_TYPE,    4, 0, 0, 0, 0 },
    { "-sensor", 0, INT_TYPE,       1, one_int_max, 0, 0, 0 },
    { "-timerange", 0, YYYYMMDDHHMMSS_TYPE, 2, 0, 0, 0, 0 }
  };

  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] ); assert( arguments );

  if ( result ) {

    /* Initialize arguments to defaults: */

    memset( arguments, 0, sizeof (Arguments) );
    arguments->tmpdir = ".";
    arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
    arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
    arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
    arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;

    /* Finish initializing non-compile-time-constant parts of options: */

    options[ 0 ].values = &arguments->listFile;
    options[ 1 ].values = &arguments->tmpdir;
    options[ 2 ].values = &arguments->bounds[ 0 ][ 0 ];
    options[ 3 ].values = &arguments->sensor;
    options[ 4 ].values = &arguments->yyyymmddhhmmss[ 0 ];

    result =
      parseOptions( argc, argv, sizeof options / sizeof *options, options );
  } else {
    fprintf( stderr, "\n%s: Invalid/insufficient command-line arguments.\n",
             argv[ 0 ] );
  }

  assert( result == 0 ||
         ( arguments->listFile && arguments->listFile[ 0 ] &&
           arguments->tmpdir && arguments->tmpdir[ 0 ] &&
           isValidBounds( (const double (*)[2]) arguments->bounds ) &&
           arguments->sensor >= 0 ) );

  return result;
}



/******************************************************************************
PURPOSE: removeTempFiles - Close and remove temp files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void removeTempFiles( Data* const data ) {
  const size_t count =
    data ? sizeof data->tempFiles / sizeof data->tempFiles[ 0 ] : 0;
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
PURPOSE: createTempFiles - Create temp output files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->ok, data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void createTempFiles( Data* const data ) {
  const size_t count =
    data ? sizeof data->tempFiles / sizeof data->tempFiles[ 0 ] : 0;
  const int pid = getpid();
  size_t index = 0;
  assert( data );
  data->ok = 1;

  for ( index = 0; data->ok && index < count; ++index ) {
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
PURPOSE: checkAndReallocateBuffer - Reallocate buffer if too small for points.
INPUTS:  const size_t bytes  Number of bytes needed for buffer.
         Data* const data    data->bufferSize  Bytes allocated for buffer.
                             data->buffer      Allocated buffer.
OUTPUTS: Data* const data  data->bufferSize    Size of buffer.
                           data->buffer        Allocated  buffer.
                           data->ok 1 if successful, else 0.
******************************************************************************/

static void checkAndReallocateBuffer( const size_t bytes, Data* const data ) {
  assert( bytes );
  assert( data );

  if ( bytes > data->bufferSize ) {
    free( data->buffer ), data->buffer = 0;
    data->bufferSize = 0;
  }

  data->bufferSize = bytes;
  data->buffer = malloc( data->bufferSize );
  data->ok = data->buffer != 0;

  if ( ! data->ok ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes "
             "to complete the requested action.\n", data->bufferSize );
    data->bufferSize = 0;
  } else {
    memset( data->buffer, 0, data->bufferSize );
  }

  assert( ( data->ok == 0 && data->buffer == 0 && data->bufferSize == 0 ) ||
          ( data->ok == 1 && data->buffer != 0 && data->bufferSize == bytes));
}



/******************************************************************************
PURPOSE: readData - Read data from each listed data file and
         write the subset of data to the temporary files.
INPUTS:  Data* data  Data to read.
OUTPUTS: Data* data  data->points, ok, tempFiles[] = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  size_t length = 0;
  size_t lines = 0;
  char* listFileContent = 0;
  data->ok = readFile( arguments->listFile, &length, &lines, &listFileContent);

  if ( data->ok ) {
    char* inputFileName = 0;
    char* end = 0;
    int fileCount = 0;
    long long yyyymmddhhmmss1 = 0; /* First file header first timestamp. */

    /* Get each line of list file. It is the data file to read: */

    for ( inputFileName = strtok_r( listFileContent, "\n", &end );
          inputFileName && data->ok;
          inputFileName = strtok_r( 0, "\n", &end ) ) {

      if ( fileSize( inputFileName ) ) { /* Ignore empty files. */
        int variables = 0;
        size_t points = 0;
        FILE* inputFile = fopen( inputFileName, "rb" );
        data->ok = inputFile != 0;

        if ( ! inputFile ) {
          fprintf( stderr, "\nFailed to open file '%s' for reading.\n",
                   inputFileName );
        } else {
          data->ok =
            readHeader( inputFile, &variables, &points,
                        data->yyyymmddhhmmss, data->header );

          if ( data->ok ) {

            if ( data->variables == 0 ) {
              data->variables = variables;
            } else if ( variables != data->variables ) {
              fprintf( stderr, "\nMismatched variable count "
                       "(actual %d, expected %d) in file '%s'.\n",
                       variables, data->variables, inputFileName );
              data->ok = 0;
            }

            if ( data->ok ) {

              if ( fileCount == 0 ) {
                yyyymmddhhmmss1 = data->yyyymmddhhmmss[ 0 ];
              }

              /* Appends to temp files: */

              extractSubset( inputFile, points, data );
            }
          }

          fclose( inputFile ), inputFile = 0;
          ++fileCount;
        }
      }
    } /* End loop on listFile. */

    data->yyyymmddhhmmss[ 0 ] = yyyymmddhhmmss1;
    free( listFileContent ), listFileContent = 0;
    data->ok = data->ok && data->points > 0;
  }

  closeTempFiles( data );
}



/******************************************************************************
PURPOSE: readHeader - Read Point xdr file header.
INPUTS:  FILE* inputFile                File to read header from.
OUTPUTS: int* variables                 Number of data variables in file.
         size_t* points                 Number of data points in file.
         long long yyyymmddhhmmss[ 2 ]  Header first and last timestamps.
         Note header[ HEADER_LINES ]    Header lines read.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int readHeader( FILE* inputFile, int* variables, size_t* points,
                       long long yyyymmddhhmmss[], Note header[] ) {
  int ok = 1;
  int index = 0;
  Note line = "";
  const size_t length = sizeof line / sizeof *line;
  assert( inputFile ); assert( yyyymmddhhmmss ); assert( header );
  memset( line, 0, sizeof line );
  *points = 0;

  for ( index = 0; ok && index < HEADER_LINES; ++index ) {
    ok = fgets( line, length, inputFile ) != 0;

    if ( ok ) {

      if ( index == 2 ) {
        ok = parseTimestamps( line, yyyymmddhhmmss );
      } else if ( index == 4 ) {
        ok = sscanf( line, "%d %lu\n", variables, points ) == 2
             && ( *variables == 6 || *variables == 7 ) && *points != 0;

        if ( ! ok ) {
          fprintf( stderr,
                   "\nFailed to read valid point count in header line '%s'\n",
                   line );
        }
      } else if ( header[ index ][ 0 ] ) {
        ok = ! strcmp( line, header[ index ] );

        if ( ! ok ) {
          fprintf( stderr, "\nFailed to read valid header line '%s'\n", line );
        }

      } else {
        memset( header[ index ], 0, sizeof (Note) );
        strncpy( header[ index ], line, length - 1 );
      }
    }
  }

  return ok;
}



/******************************************************************************
PURPOSE: parseTimestamps - Parse header timestamps.
INPUTS:  const char* line      Header line to parse.
OUTPUTS: long long yyyymmddhhmmss[ 2 ]  Header first and last timestamps.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int parseTimestamps( const char* line, long long yyyymmddhhmmss[] ) {
  int ok = 0;
  int yyyy1 = 0;
  int mm1   = 0;
  int dd1   = 0;
  int hh1   = 0;
  int min1  = 0;
  int ss1   = 0;
  int yyyy2 = 0;
  int mm2   = 0;
  int dd2   = 0;
  int hh2   = 0;
  int min2  = 0;
  int ss2   = 0;
  assert( line ); assert( yyyymmddhhmmss );
  yyyymmddhhmmss[ 0 ] = yyyymmddhhmmss[ 1 ] = 0;
  ok = sscanf( line, "%d-%d-%dT%d:%d:%d-0000 %d-%d-%dT%d:%d:%d-0000\n",
               &yyyy1, &mm1, &dd1, &hh1, &min1, &ss1,
               &yyyy2, &mm2, &dd2, &hh2, &min2, &ss2 ) == 12;

  if ( ok ) {
    yyyymmddhhmmss[ 0 ] = yyyy1;
    yyyymmddhhmmss[ 0 ] *= 100;
    yyyymmddhhmmss[ 0 ] += mm1;
    yyyymmddhhmmss[ 0 ] *= 100;
    yyyymmddhhmmss[ 0 ] += dd1;
    yyyymmddhhmmss[ 0 ] *= 100;
    yyyymmddhhmmss[ 0 ] += hh1;
    yyyymmddhhmmss[ 0 ] *= 100;
    yyyymmddhhmmss[ 0 ] += min1;
    yyyymmddhhmmss[ 0 ] *= 100;
    yyyymmddhhmmss[ 0 ] += ss1;
    ok = isValidYYYYMMDDHHMMSS( yyyymmddhhmmss[ 0 ] );

    if ( ok ) {
      yyyymmddhhmmss[ 1 ] = yyyy2;
      yyyymmddhhmmss[ 1 ] *= 100;
      yyyymmddhhmmss[ 1 ] += mm2;
      yyyymmddhhmmss[ 1 ] *= 100;
      yyyymmddhhmmss[ 1 ] += dd2;
      yyyymmddhhmmss[ 1 ] *= 100;
      yyyymmddhhmmss[ 1 ] += hh2;
      yyyymmddhhmmss[ 1 ] *= 100;
      yyyymmddhhmmss[ 1 ] += min2;
      yyyymmddhhmmss[ 1 ] *= 100;
      yyyymmddhhmmss[ 1 ] += ss2;
      ok = isValidYYYYMMDDHHMMSS( yyyymmddhhmmss[ 1 ] );
    }
  }

  if ( ! ok ) {
    fprintf( stderr,
             "\nFailed to read valid timestamps in header line '%s'\n", line );
  }

  return ok;
}



/******************************************************************************
PURPOSE: extractSubset - Extract data within bounds and write to temp files.
INPUTS:  FILE* inputFile                     File to read from.
         const size_t points                 Number of data points to read.
         Data* data  data->arguments.bounds  Lon-lat bounds to subset by.
OUTPUTS: Data* data  data->tempFiles[]  Write variable data to these files.
                     data->ok = 1 if successful, else 0.
******************************************************************************/

static void extractSubset( FILE* inputFile, const size_t points, Data* data ) {
  const int variables = data ? data->variables : 0;
  size_t bytes = points * sizeof (Note) + data->variables * points * sizeof (double);
  assert( inputFile ); assert( points );
  assert( data );
  assert( data->variables == 6 || data->variables == 7 );
  assert( data->tempFiles[ 0 ] );
  checkAndReallocateBuffer( bytes, data );

  if ( data->ok ) {
    char* const buffer = data->buffer;
    data->ok = fread( buffer, 1, bytes, inputFile ) == bytes;

    if ( ! data->ok ) {
      fprintf( stderr, "\nFailed to read %lu bytes of input file.\n", bytes );
    } else {
      const Arguments* const arguments = &( data->arguments );
      const int sensor = arguments->sensor;
      const double longitudeMinimum = arguments->bounds[ LONGITUDE ][ MINIMUM ];
      const double longitudeMaximum = arguments->bounds[ LONGITUDE ][ MAXIMUM ];
      const double latitudeMinimum  = arguments->bounds[ LATITUDE  ][ MINIMUM ];
      const double latitudeMaximum  = arguments->bounds[ LATITUDE  ][ MAXIMUM ];
      const long long yyyymmddhhmmss0 = arguments->yyyymmddhhmmss[ 0 ];
      const long long yyyymmddhhmmss1 = arguments->yyyymmddhhmmss[ 1 ];
      Note* notes = (Note*) buffer;
      double* const dp = (double*) ( buffer + points * sizeof (Note) );
      double* const timestamps = dp;
      double* const longitudes = timestamps + points;
      double* const latitudes  = longitudes + points;
      double* const elevations = latitudes  + points;
      double* const ids        = elevations + points;
      double* const counts     = variables == 7 ? ids + points : 0;
      double* const measures   = counts ? counts + points : ids + points;
      const int filtering = /* Check to avoid comparisons. */
        sensor != 0 ||
        longitudeMinimum > -180.0 ||
        longitudeMaximum <  180.0 ||
        latitudeMinimum  >  -90.0 ||
        latitudeMaximum  <   90.0 ||
        yyyymmddhhmmss0 != 0;
      size_t outputPoints = filtering ? 0 : points;

      if ( filtering ) {
        size_t point = 0;
        rotate8ByteArrayIfLittleEndian( timestamps, variables * points );

        for ( point = 0; point < points; ++point ) {

          if ( sensor == 0 || ids[ point ] == sensor ) {
            const double longitude = longitudes[ point ];

            if ( longitude >= longitudeMinimum &&
                 longitude <= longitudeMaximum) {
              const double latitude = latitudes[ point ];

              if ( latitude >= latitudeMinimum &&
                   latitude <= latitudeMaximum ) {
                const double timestamp = timestamps[ point ];
                const long long yyyymmddhhmmss = (long long) timestamp;

                if ( yyyymmddhhmmss0 == 0 ||
                     ( yyyymmddhhmmss >= yyyymmddhhmmss0 &&
                       yyyymmddhhmmss <= yyyymmddhhmmss1 ) ) {

                  if ( outputPoints < point ) {
                    const double elevation = elevations[ point ];
                    const double id        = ids[        point ];
                    const double count     = counts ? counts[ point ] : 0.0;
                    const double measure   = measures[   point ];
                    timestamps[ outputPoints ] = timestamp;
                    longitudes[ outputPoints ] = longitude;
                    latitudes[  outputPoints ] = latitude;
                    elevations[ outputPoints ] = elevation;
                    ids[        outputPoints ] = id;

                    if ( counts ) {
                      counts[   outputPoints ] = count;
                    }

                    measures[   outputPoints ] = measure;
                    memcpy( notes[ outputPoints ], notes[ point ],
                            sizeof (Note) );
                  }

                  ++outputPoints;
                }
              }
            }
          }
        }
      }

      if ( outputPoints ) {

        if ( filtering ) {
          rotate8ByteArrayIfLittleEndian( timestamps, outputPoints );
          rotate8ByteArrayIfLittleEndian( longitudes, outputPoints );
          rotate8ByteArrayIfLittleEndian( latitudes,  outputPoints );
          rotate8ByteArrayIfLittleEndian( elevations, outputPoints );
          rotate8ByteArrayIfLittleEndian( ids,        outputPoints );

          if ( counts ) {
            rotate8ByteArrayIfLittleEndian( counts,   outputPoints );
          }

          rotate8ByteArrayIfLittleEndian( measures,   outputPoints );
        }

        bytes = outputPoints * sizeof (Note);
        data->ok = fwrite( notes, 1, bytes, data->tempFiles[ 0 ] ) == bytes;

        if ( data->ok ) {
          bytes = outputPoints * sizeof (double);
          data->ok =
            fwrite( timestamps, 1, bytes, data->tempFiles[ 1 ] ) == bytes &&
            fwrite( longitudes, 1, bytes, data->tempFiles[ 2 ] ) == bytes &&
            fwrite( latitudes,  1, bytes, data->tempFiles[ 3 ] ) == bytes &&
            fwrite( elevations, 1, bytes, data->tempFiles[ 4 ] ) == bytes &&
            fwrite( ids,        1, bytes, data->tempFiles[ 5 ] ) == bytes &&
            ( counts == 0 ||
              fwrite( counts,   1, bytes, data->tempFiles[ 6 ] ) == bytes ) &&
            fwrite( measures,   1, bytes, data->tempFiles[ 7 ] ) == bytes;

          if ( data->ok ) {
            data->points += outputPoints;
          }
        }
      }

      if ( ! data->ok ) {
        bytes =
          outputPoints * sizeof (Note) +
          variables * outputPoints * sizeof (double);
        fprintf( stderr,
                 "\nFailed to write %lu bytes of subset data to temp files.\n",
                bytes );
      } else {
        data->variables = variables;
      }
    }
  }
}



/******************************************************************************
PURPOSE: streamData - Write final content of temp files to stdout.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* const data  data->ok, tempFiles[] = 0 (closed and removed).
******************************************************************************/

static void streamData( Data* const data ) {

  assert( data ); assert( data->ok ); assert( data->points > 0 );
  assert( data->buffer ); assert( data->bufferSize );

  /* Temp file names are initialized but temp files are closed after writing:*/

  assert( data->tempFileNames[ 0 ][ 0 ] );
  assert( data->tempFileNames[ TEMP_FILES - 1 ][ 0 ] );
  assert( data->tempFiles[ 0 ] == 0 );
  assert( data->tempFiles[ TEMP_FILES - 1 ] == 0 );

  {
    int fileIndex = 0;
    streamHeader( data );

    for ( fileIndex = 0; data->ok && fileIndex < TEMP_FILES; ++fileIndex ) {
      const char* const tempFileName = data->tempFileNames[ fileIndex ];
      FILE* tempFile = data->tempFiles[fileIndex] = fopen( tempFileName, "rb" );
      data->ok = tempFile != 0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nCan't open temp data file '%s' for reading.\n",
                 tempFileName );
      } else {
        const size_t bytes = data->bufferSize;
        char* const buffer = data->buffer;

        while ( data->ok && ! feof( tempFile ) ) {
          const size_t bytesRead = fread( buffer, 1, bytes, tempFile );

          if ( bytesRead ) {
            const size_t bytesWritten = fwrite( buffer, 1, bytesRead, stdout );
            data->ok = bytesWritten == bytesRead;
          }
        }

        if ( ! data->ok ) {
          fprintf( stderr,
                   "\nFailed to stream subset data from temp file '%s'.\n",
                   tempFileName );
        }
      }
    }
  }

  removeTempFiles( data );
}



/******************************************************************************
PURPOSE: streamHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
******************************************************************************/

static void streamHeader( Data* const data ) {
  int line = 0;
  assert( data );
  assert( data->header[ 0 ][ 0 ] );
  assert( data->header[ HEADER_LINES - 1 ][ 0 ] );

  for ( line = 0; line < HEADER_LINES; ++line ) {
    char* const headerLine = data->header[ line ];

    if ( line == 1 ) {
      char* s = strrchr( headerLine, '\n' );

      if ( s ) {
        *s = '\0';
      }

      fputs( headerLine, stdout );
      fputs( ",PointSubset\n", stdout );
    } else if ( line == 2 ) {
      const Arguments* const arguments = &( data->arguments );
      const long long yyyymmddhhmmssFirst = arguments->yyyymmddhhmmss[ 0 ];
      const long long yyyymmddhhmmssLast = arguments->yyyymmddhhmmss[ 1 ];

      const long long yyyymmddhhmmss1 =
        yyyymmddhhmmssFirst ? yyyymmddhhmmssFirst : data->yyyymmddhhmmss[ 0 ];
      const long long yyyymmddhhmmss2 =
        yyyymmddhhmmssLast ? yyyymmddhhmmssLast : data->yyyymmddhhmmss[ 1 ];
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
      printf( "%04d-%02d-%02dT%02d:%02d:%02d-0000 "
              "%04d-%02d-%02dT%02d:%02d:%02d-0000\n",
              yyyy1, mm1, dd1, hh1, mm21, ss1,
              yyyy2, mm2, dd2, hh2, mm22, ss2 );
    } else if ( line == 4 ) {
      printf( "%d %lu\n", data->variables, data->points );
    } else {
      fputs( headerLine, stdout );
    }
  }
}




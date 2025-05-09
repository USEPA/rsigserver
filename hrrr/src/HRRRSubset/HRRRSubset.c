/******************************************************************************
PURPOSE: HRRRSubset.c - Extract a lon-lat subset of data from a list of
         HRRR GRIB2 files and write it to stdout as XDR binary format.

 NOTES:  HRRR model variables are described here:
         https://home.chpc.utah.edu/~u0553130/Brian_Blaylock/HRRR_archive/\
         hrrr_sfc_table_f00-f01.html

         Uses GRIB2 Library.

         Compile:
         gcc -Wall -g -o HRRRSubset HRRRSubset.c Utilities.c ReadData.c \
                   -I../../../include \
                   -L../../../lib/$platform \
                   -lgrib2 -lm

         Usage:
         HRRRSubset -lonlats <lonlatfile> -files <listfile> \
                      -tmpdir <temp_directory> \
                      -desc "description text" \
                      -timestamp <yyyymmddhh> -hours <count> \
                      -variable <name> \
                      -units <name> \
                      -domain <minimum_longitude> <minimum_latitude> \
                              <maximum_longitude> <maximum_latitude> \
                      [-is_vector2] \
                      [-cmaq] \

          Example:
          ../../../bin/$platform/HRRRSubset \
          -lonlats testdata/HRRR_lonlat.bin \
          -files testdata/file_list \
          -tmpdir testdata \
          -desc "http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,HRRRSubset" \
          -timestamp 2020021700 -hours 24 \
          -variable wind_10m \
          -units m/s \
          -domain -75 35 -70 36 \
          -is_vector2 \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 12-line ASCII header followed by binary 64-bit big-endian arrays:

Grid 1.0
http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,HRRRSubset
2020-02-17T00:00:00-0000
# Dimensions: timesteps variables rows columns:
24 2 1059 1799
# Variable names:
wind_10m_u wind_10m_v
# Variable units:
m/s m/s
# IEEE-754 64-bit reals longitudes[rows][columns] and
# IEEE-754 64-bit reals latitudes[rows][columns] and
# IEEE-754 64-bit reals data[timesteps][variables][rows][columns]:
<big-endian binary format arrays>
-1.2272000000000000e+02
-1.2269333000000000e+02
-1.2266665999999999e+02
...
-6.0990960000000001e+01
-6.0954389999999997e+01
-6.0917839999999998e+01
2.1138000000000002e+01
2.1144990000000000e+01
2.1151969999999999e+01
...
4.7862949999999998e+01
4.7852670000000003e+01
4.7842379999999999e+01
-1.234567890123456e+00
...

HISTORY: 2020-02-21 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For FILE, printf(), snprintf(). */
#include <string.h> /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <ctype.h>  /* For isspace(), isprint(), isalpha(), isdigit(). */
#include <stdlib.h> /* For malloc(), free(), atof(). */
#include <float.h>  /* For MAX_FLT. */
#include <unistd.h> /* For unlink(), getpid() */

#include "Utilities.h" /* For LONGITUDE, Bounds, isValidYYYYMMDDHH(). */
#include "ReadData.h"  /* For readData(). */

/*================================= MACROS =================================*/

/* Name of temporary file created in -tmpdir will have PID appended: */

#define TEMP_FILE_NAME "junk_HRRRSubset"

#define AND2(a, b) ( (a) && (b) )
#define AND3(a, b, c) ( (a) && (b) && (c) )

/*================================== TYPES ==================================*/

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* lonlatFile;  /* File containing lonlat coordinates to read. */
  const char* listFile;    /* File containing list of HRRR files to read. */
  const char* tmpdir;      /* Name of directory to write temp files. */
  const char* description; /* User-supplied description. */
  const char* variable;    /* Name of variable to read.  E.g.,"wind_10m". */
  const char* units;       /* Variable units. E.g., "m/s". */
  Bounds      domain; /* Subset domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  int         yyyymmddhh;  /* First timestamp of subset. */
  int         hours;       /* Number of hours in subset. */
  int         isVector2;   /* Is variable a 2D vector? */
  int         cmaq;        /* Output CMAQ implicit grid XDR format? */
} Arguments;

/* Data type: */

typedef struct {
  Arguments  arguments;     /* User-supplied (command-line) arguments. */
  FileName   tempFileNames[2]; /* Name of temp file of output subset data. */
  FILE*      tempFiles[2];     /* Temp file of output subset data. */
  size_t     rows;          /* Number of rows    in grid. */
  size_t     columns;       /* Number of columns in grid. */
  size_t     subsetIndices[ 2 ][ 2 ]; /* [ COLUMN,ROW ][ MINIMUM,MAXIMUM]. */
  double*    longitudesLatitudes; /* longitudesLatitudes[2 * rows * columns] */
  double*    buffer;        /* buffer[ rows * columns * 4 ] */
                            /* for one timestep of wind_10m_u, wind_10m_v */
                            /* plus write temp space. */
  int        ok;            /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocate( Data* data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static void readDataFiles( Data* const data );

static int dataFileTimestamp( const char* const fileName );

static double* readCoordinatesFile( const char* const fileName,
                                    size_t* const rows, size_t* const columns);

static void writeDataSubset( Data* const data );

static void streamData( Data* const data );

static void streamHeader( const Data* const data );


/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of HRRR files and write
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
    readDataFiles( &data ); /* Read each subset of HRRR files, write temp file*/

    if ( data.ok ) {
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

  if ( data->longitudesLatitudes ) {
    free( data->longitudesLatitudes );
    data->longitudesLatitudes = 0;
  }

  if ( data->buffer ) {
    free( data->buffer );
    data->buffer = 0;
  }

  if ( data->tempFiles[ 0 ] ) {
    fclose( data->tempFiles[ 0 ] );
    data->tempFiles[ 0 ] = 0;
  }

  if ( data->tempFiles[ 1 ] ) {
    fclose( data->tempFiles[ 1 ] );
    data->tempFiles[ 1 ] = 0;
  }

  if ( data->tempFileNames[ 0 ][ 0 ] ) {
    unlink( data->tempFileNames[ 0 ] );
    memset( data->tempFileNames[ 0 ], 0, sizeof data->tempFileNames[ 0 ] );
  }

  if ( data->tempFileNames[ 1 ][ 0 ] ) {
    unlink( data->tempFileNames[ 1 ] );
    memset( data->tempFileNames[ 1 ], 0, sizeof data->tempFileNames[ 1 ] );
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
          "HRRR GRIB2 files and write it to stdout as XDR binary format.\n",
           name );
  fprintf( stderr, "Data is subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "  -lonlats <lonlatfile> \\\n" );
  fprintf( stderr, "  -files <listfile> \\\n" );
  fprintf( stderr, "  -tmpdir <temp_directory> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, "  -units <name> \\\n" );
  fprintf( stderr, "  -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, "  [-is_vector2] \\\n" );
  fprintf( stderr, "  [-cmaq] \n\n" );
  fprintf( stderr, "Note:\ntimestamp is in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory were a transient file is " );
  fprintf ( stderr, "written.\nIt should have enough disk space (1TB).\n\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-lonlats testdata/HRRR_lonlat.bin \\\n");
  fprintf( stderr, "-files testdata/filelist \\\n");
  fprintf( stderr, "-tmpdir testdata \\\n");
  fprintf( stderr, "-desc "
                   "\"http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/"
                   ",HRRRSubset\" \\\n");
  fprintf( stderr, "-timestamp 2020021700 -hours 24 \\\n" );
  fprintf( stderr, "-variable wind_10m \\\n" );
  fprintf( stderr, "-units m/s \\\n" );
  fprintf( stderr, "-domain -126 25 -65 50 \\\n" );
  fprintf( stderr, "-is_vector2 \\\n" );
  fprintf( stderr, "> subset.xdr\n\n" );
  fprintf( stderr, "HRRR modeled 2D wind at 10m above ground over US on "
           "February 17, 2020.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr,
           "Grid 1.0\n"
           "http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/,HRRRSubset\n"
           "2020-02-17T00:00:00-0000\n"
           "# Dimensions: timesteps variables rows columns:\n"
           "24 2 rows columns\n"
           "# Variable names:\n"
           "wind_10m_u wind_10m_v\n"
           "# Variable units:\n"
           "m/s m/s\n"
           "# IEEE-754 64-bit reals longitudes[rows][columns] and\n"
           "# IEEE-754 64-bit reals latitudes[rows][columns] and\n"
           "# IEEE-754 64-bit reals data[timesteps][variables][rows][columns]:\n"
           "<big-endian binary format arrays>\n"
           "-1.2272000000000000e+02\n"
           "-1.2269333000000000e+02\n"
           "-1.2266665999999999e+02\n"
           "...\n"
           "-6.0990960000000001e+01\n"
           "-6.0954389999999997e+01\n"
           "-6.0917839999999998e+01\n"
           "2.1138000000000002e+01\n"
           "2.1144990000000000e+01\n"
           "2.1151969999999999e+01\n"
           "...\n"
           "4.7862949999999998e+01\n"
           "4.7852670000000003e+01\n"
           "4.7842379999999999e+01\n"
           "-1.234567890123456e+00\n"
           "...\n" );
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

  result = argc >= 22 && argc <= 24;

  for ( arg = 1; result && arg < argc; ++arg ) {

    if ( ! strcmp( argv[ arg ], "-lonlats" ) && arg + 1 < argc ) {
      ++arg;
      arguments->lonlatFile = argv[ arg ];
      result = arguments->lonlatFile[ 0 ] != '\0';
    } else if ( ! strcmp( argv[ arg ], "-files" ) && arg + 1 < argc ) {
      ++arg;
      arguments->listFile = argv[ arg ];
      result = arguments->listFile[ 0 ] != '\0';
    } else if ( ! strcmp( argv[ arg ], "-tmpdir" ) && arg + 1 < argc ) {
      ++arg;
      arguments->tmpdir = argv[ arg ];
      result = arguments->tmpdir[ 0 ] != '\0';
    } else if ( ! strcmp( argv[ arg ], "-desc" ) && arg + 1 < argc ) {
      ++arg;
      arguments->description = argv[ arg ];
      result = arguments->description[ 0 ] != '\0';
    } else if ( ! strcmp( argv[ arg ], "-timestamp" ) && arg + 1 < argc ) {
      ++arg;
      arguments->yyyymmddhh = atoi( argv[ arg ] );
      result = isValidYYYYMMDDHH( arguments->yyyymmddhh );
    } else if ( ! strcmp( argv[ arg ], "-hours" ) && arg + 1 < argc ) {
      ++arg;
      arguments->hours = atoi( argv[ arg ] );
      result = arguments->hours > 0;
    } else if ( ! strcmp( argv[ arg ], "-variable" ) && arg + 1 < argc ) {
      int index = 0;
      ++arg;
      arguments->variable = argv[ arg ];
      result = isalpha( arguments->variable[ 0 ] );

      for ( index = 1; arguments->variable[ index ] && result; ++index  ) {
        const char ch = arguments->variable[ index ];
        result = isalpha( ch ) || isdigit( ch ) || ch == '_';
      }
    } else if ( ! strcmp( argv[ arg ], "-units" ) && arg + 1 < argc ) {
      int index = 0;
      ++arg;
      arguments->units = argv[ arg ];
      result =
        isprint( arguments->units[ 0 ] ) && ! isspace( arguments->units[ 0 ] );

      for ( index = 1; arguments->units[ index ] && result; ++index  ) {
        const char ch = arguments->units[ index ];
        result = isprint( ch ) && ! isspace( ch );
      }
    } else if ( ! strcmp( argv[ arg ], "-is_vector2" ) && arg < argc ) {
      arguments->isVector2 = 1;
    } else if ( ! strcmp( argv[ arg ], "-cmaq" ) && arg < argc ) {
      arguments->cmaq = 1;
    } else if ( ! strcmp( argv[ arg ], "-domain" ) && arg + 4 < argc ) {
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
PURPOSE: readDataFiles - Read data from each listed HRRR file and
         write the lon-lat subset of data to the temporary file.
INPUTS:  Data* data  Data description to read.
OUTPUTS: Data* data  data->yyyydddhhmm, points, scans, ok, tempFiles = 0.
******************************************************************************/

static void readDataFiles( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  data->longitudesLatitudes =
    readCoordinatesFile( arguments->lonlatFile, &data->rows, &data->columns );
  data->ok = data->longitudesLatitudes != 0;

  if ( data->ok ) {
    const size_t rows    = data->rows;
    const size_t columns = data->columns;
    const size_t count = rows * columns;
    const double* const longitudes = data->longitudesLatitudes;
    const double* const latitudes  = longitudes + count;
    int wroteSomeData = 0;

    data->ok =
      subsetIndicesByBounds( (const double (*)[2]) data->arguments.domain,
                             rows, columns, longitudes, latitudes,
                             &data->subsetIndices[ ROW ][ MINIMUM ],
                             &data->subsetIndices[ ROW ][ MAXIMUM ],
                             &data->subsetIndices[ COLUMN ][ MINIMUM ],
                             &data->subsetIndices[ COLUMN ][ MAXIMUM ] );

    if ( data->ok ) {
      const size_t count2 = count + count;
      const size_t count4 = count2 + count2;
      data->buffer = (double*) malloc( count4 * sizeof (double) );
      data->ok = data->buffer != 0;

      if ( ! data->ok ) {
          fprintf( stderr,
               "\nFailed to allocate %lu bytes "
               "to complete the requested action.\n",
               count4 * sizeof (double) );
      } else {
        size_t unused_ = 0;
        char* listFileContent = readFile( arguments->listFile, &unused_ );
        data->ok = listFileContent != 0;

        if ( data->ok ) {
          const int hoursPerTimestep = 1;
          char* fileName = 0;
          char* end      = 0;
          int yyyymmddhh = arguments->yyyymmddhh;
          const int hours = arguments->hours;
          int hoursWritten = 0;

          /* Get each line of list file. It is the HRRR data file to read: */

          for ( fileName = strtok_r( listFileContent, "\n", &end );
                AND2( data->ok, fileName );
                fileName = strtok_r( 0, "\n", &end ) ) {
            const int fileYYYYMMDDHH = dataFileTimestamp( fileName );

            /* Write missing values for hours before first data file: */

            if ( yyyymmddhh < fileYYYYMMDDHH ) {
              fillArray( MISSING_VALUE, count2, data->buffer );

              while ( AND3( yyyymmddhh < fileYYYYMMDDHH,
                            hoursWritten < hours,
                            data->ok ) ) {
                writeDataSubset( data ); /* Sets data->ok. */
                ++hoursWritten;
                yyyymmddhh = incrementHours( yyyymmddhh, hoursPerTimestep );
              }
            }

            /*
             * If data file timestamp matches hour then
             * try to read data, if failed substitute missing values,
             * then write data:
             */

            if ( AND3( data->ok,
                       fileYYYYMMDDHH == yyyymmddhh,
                       hoursWritten < hours ) ) {
              const int readSomeData =
                readData( fileName, arguments->isVector2, count, data->buffer );

              if ( ! readSomeData ) {
                fillArray( MISSING_VALUE, count2, data->buffer );
              }

              writeDataSubset( data ); /* Sets data->ok. */

              if ( AND2( readSomeData, data->ok ) ) {
                wroteSomeData = 1;
              }

              ++hoursWritten;
              yyyymmddhh = incrementHours( yyyymmddhh, hoursPerTimestep );
            }
          } /* End loop on listFile. */

          /* Write missing values for hours after last data file: */

          if ( AND3( hoursWritten < hours, wroteSomeData, data->ok ) ) {
            fillArray( MISSING_VALUE, count2, data->buffer );

            while ( AND2( hoursWritten < hours, data->ok ) ) {
              writeDataSubset( data ); /* Sets data->ok. */
              ++hoursWritten;
            }
          }

          free( listFileContent );
          listFileContent = 0;
        }
      }

      /* Done writing tempFiles so close them: */

      if ( data->tempFiles[ 0 ] ) {
        fclose( data->tempFiles[ 0 ] );
        data->tempFiles[ 0 ] = 0;
      }

      if ( data->tempFiles[ 1 ] ) {
        fclose( data->tempFiles[ 1 ] );
        data->tempFiles[ 1 ] = 0;
      }
    }

    data->ok = wroteSomeData;
  }
}



/******************************************************************************
PURPOSE: dataFileTimestamp - Timestamp of file.
INPUTS:  const char* const fileName   Name of HRRR file.
RETURNS: int yyyymmddhh of file or 0 if failed and message on stderr.
******************************************************************************/

static int dataFileTimestamp( const char* const fileName ) {
  int result = 0;
  const char* s = strrchr( fileName, '/' );

  if ( s ) {
    int c = 0;
    ++s;

    /* Parse YYYYMMDD: */

    for ( c = 0; *s && *s != '_' && c < 8; ++c, ++s ) {
      result = result * 10 + *s - '0';
    }

    /* Parse YYYYMMDD: */

    s = strstr( s, ".t" );

    if ( s ) {
      s += 2;

      for ( c = 0; *s && *s != 'z' && c < 2; ++c, ++s ) {
        result = result * 10 + *s - '0';
      }
    }
  }

  if ( ! isValidYYYYMMDDHH( result ) ) {
    fprintf( stderr, "\nInvalid file name timestamp '%s'.\n", fileName );
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readCoordinatesFile - Read longitude-latitude coordinates file.
INPUTS:  const char* const fileName  Name of file to read.
OUTPUTS: size_t* const rows          Number of rows of points.
         size_t* const columns       Number of columns of points.
RETURNS: double* array of coordinates: array[ ROWS * COLUMNS ].
******************************************************************************/

static double* readCoordinatesFile( const char* const fileName,
                                   size_t* const rows, size_t* const columns) {

  double* result = 0;
  int ok = 0;
  assert( fileName ); assert( *fileName ); assert( rows ); assert( columns );

  {
    FILE* inputFile = fopen( fileName, "rb" );
    ok = inputFile != 0;

    if ( ok ) {
      const char* const format =
        "Content-type: application/octet-stream; charset=iso-8859-1\n"
        "# dimensions: variables rows columns\n"
        "2 %lu %lu\n"
        "# variable names:\n"
        "longitude latitude\n"
        "# variable units:\n"
        "deg deg\n"
        "# IEEE-754 64-bit real data[variables][rows][columns]:\n";
      ok = fscanf( inputFile, format, rows, columns ) == 2 &&
           rows && columns;

      if ( ok ) {
        const size_t count = *rows * *columns;
        const size_t count2 = count + count;
        result = (double*) malloc( count2 * sizeof (double) );
        ok = result != 0;

        if ( ! ok ) {
          fprintf( stderr,
               "\nFailed to allocate %lu bytes "
               "to complete the requested action.\n",
               count2 * sizeof (double) );
        } else {
          ok = fread( result, sizeof *result, count2, inputFile ) == count2;
          rotate8ByteArrayIfLittleEndian( result, count2 );

          if ( ok ) {
            double* const longitudes = result;
            double* const latitudes  = longitudes + count;
            size_t point = 0;

            for ( point = 0; ok && point < count; ++point ) {
              const double longitude = longitudes[ point ];
              const double latitude  = latitudes[  point ];
              ok =
                IN_RANGE( longitude, -180.0, 180.0 ) &&
                IN_RANGE( latitude,   -90.0,  90.0 );
            }
          }
        }
      }

      fclose( inputFile );
      inputFile = 0;
    }
  }

  if ( ! ok ) {
    fprintf( stderr, "\nInvalid coordinates file '%s'.\n", fileName );
    free( result );
    result = 0;
    *rows = *columns = 0;
  }

  assert( ( result == 0 && *rows == 0 && *columns == 0 ) ||
          ( result != 0 && *rows && *columns &&
            IN_RANGE( result[ 0 ], -180.0, 180.0 ) &&
            IN_RANGE( result[ *rows * *columns - 1 ], -180.0, 180.0 ) &&
            IN_RANGE( result[ *rows * *columns ], -90.0, 90.0 ) &&
            IN_RANGE( result[ *rows * *columns * 2 - 1 ], -90.0, 90.0 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: writeDataSubset - Write one hour subset of data to temp file.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* data  data->ok.
******************************************************************************/

static void writeDataSubset( Data* const data ) {
  int writeCoordinates = 0;
  assert( data ); assert( data->ok );
  assert( data->longitudesLatitudes ); assert( data->buffer );

  /* Open temp file for writing if it does not yet exist: */

  if ( data->tempFiles[ 0 ] == 0 ) {
    const int pid = getpid();
    const int files = 1 + data->arguments.isVector2;
    int file = 0;
    assert( data->tempFiles[ 1 ] == 0 );
    writeCoordinates = data->arguments.cmaq == 0;

    for ( file = 0; data->ok && file < files; ++file ) {
      memset( data->tempFileNames[ file ], 0, sizeof (FileName) );
      snprintf( data->tempFileNames[ file ],
                sizeof (FileName) / sizeof (char) - 1,
               "%s/%s%d.%04d",
                data->arguments.tmpdir, TEMP_FILE_NAME, file, pid );
      data->tempFiles[ file ] = fopen( data->tempFileNames[ file ], "wb" );
      data->ok = data->tempFiles[ file ] != 0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nCan't create temporary output file '%s'.\n",
                data->tempFileNames[ file ] );
      }
    }
  }

  if ( data->ok ) { /* Write subset data to temp file: */
    const size_t firstColumn = data->subsetIndices[ COLUMN ][ MINIMUM ];
    const size_t lastColumn  = data->subsetIndices[ COLUMN ][ MAXIMUM ];
    const size_t firstRow    = data->subsetIndices[ ROW ][ MINIMUM ];
    const size_t lastRow     = data->subsetIndices[ ROW ][ MAXIMUM ];
    const size_t subsetPoints =
      ( 1 + lastRow - firstRow ) * ( 1 + lastColumn - firstColumn );
    const size_t subsetPoints2 = subsetPoints + subsetPoints;
    const size_t columns = data->columns;
    const size_t count = data->rows * columns;
    const size_t count2 = count + count;
    const double* const longitudes = data->longitudesLatitudes;
    double* buffer = data->buffer + count2; /* Temp memory for rotate. */

    if ( writeCoordinates ) {
      const double* const latitudes = longitudes + count;
      double* subsetLongitudes = buffer;
      double* subsetLatitudes  = subsetLongitudes + subsetPoints;
      size_t row = 0;

      for ( row = firstRow; row <= lastRow; ++row ) {
        const size_t rowOffset = row * columns;
        size_t column = 0;

        for ( column = firstColumn; column <= lastColumn; ++column ) {
          const size_t index = rowOffset + column;
          const double longitude = longitudes[ index ];
          const double latitude  = latitudes[  index ];
          assert( IN_RANGE( longitude, -180.0, 180.0 ) );
          assert( IN_RANGE( latitude,   -90.0,  90.0 ) );
          *subsetLongitudes++ = longitude;
          *subsetLatitudes++  = latitude;
        }
      }

      rotate8ByteArrayIfLittleEndian( buffer, subsetPoints2 );

      data->ok =
        fwrite( buffer, sizeof *buffer, subsetPoints2, data->tempFiles[0 ] )
        == subsetPoints2;

      if ( ! data->ok ) {
        fprintf( stderr,
                 "\nFailed to write subset coordinates to temp file '%s'.\n",
                 data->tempFileNames[ 0 ] );
      }
    }

    if ( data->ok ) { /* Copy subset data to temp part of buffer: */
      const double* const variable1 = data->buffer;
      const double* const variable2 =
        data->arguments.isVector2 ? variable1 + count : 0;
      double* subsetVariable1 = buffer;
      double* subsetVariable2 = variable2 ? subsetVariable1 + subsetPoints : 0;
      const size_t variableSubsetPoints =
        variable2 ? subsetPoints2 : subsetPoints;
      const size_t wordSize = data->arguments.cmaq == 0 ? 8 : 4;
      size_t row = 0;

      for ( row = firstRow; row <= lastRow; ++row ) {
        const size_t rowOffset = row * columns;
        size_t column = 0;

        for ( column = firstColumn; column <= lastColumn; ++column ) {
          const size_t index = rowOffset + column;
          *subsetVariable1++ = variable1[ index ];

          if ( variable2 ) {
            *subsetVariable2++ = variable2[ index ];
          }
        }
      }

      /* Write subset data at temp part of buffer: */

      if ( wordSize == 8 ) { /* GRID XDR format 64-bit data: */
        rotate8ByteArrayIfLittleEndian( buffer, variableSubsetPoints );
        data->ok =
          fwrite( buffer, wordSize, variableSubsetPoints, data->tempFiles[ 0 ])
          == variableSubsetPoints;

        if ( ! data->ok ) {
          fprintf( stderr, "\nFailed to write subset data to temp file '%s'.\n",
                   data->tempFileNames[ 0 ] );
        }
      } else { /* CMAQ XDR format 32-bit data: */
        const int files = 1 + data->arguments.isVector2;
        int file = 0;
        float* fbuffer = (float*) buffer;
        doublesToFloats( buffer, variableSubsetPoints );
        rotate4ByteArrayIfLittleEndian( fbuffer, variableSubsetPoints );

        for ( file = 0; data->ok && file < files; ++file,
              fbuffer += subsetPoints ) {
          data->ok =
            fwrite( fbuffer, wordSize, subsetPoints, data->tempFiles[ file ] )
            == subsetPoints;

          if ( ! data->ok ) {
            fprintf( stderr,
                     "\nFailed to write subset data to temp file '%s'.\n",
                     data->tempFileNames[ file ] );
          }
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: streamData - Stream content of binary temp file to stdout.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* const data  data->ok, data->buffer contents overwritten.
******************************************************************************/

static void streamData( Data* const data ) {
  const int files = data ? 1 + data->arguments.isVector2 : 0;
  int file = 0;
  assert( data ); assert( data->ok ); assert( data->buffer );
  assert( data->tempFileNames[ 0 ][ 0 ] );
  assert( data->tempFiles[ 0 ] == 0 ); /* Temp file is closed after writing.*/
  assert( data->tempFiles[ 1 ] == 0 ); /* Temp file is closed after writing.*/

  /* Open temp files: */

  for ( file = 0; data->ok && file < files; ++file ) {
    assert( data->tempFileNames[ file ][ 0 ] );
    data->tempFiles[ file ] = fopen( data->tempFileNames[ file ], "rb" );
    data->ok = data->tempFiles[ file ] != 0;

    if ( ! data->ok ) {
      fprintf( stderr, "\nCan't open temp data file '%s' for reading.\n",
                data->tempFileNames[ file ] );
    }
  }

  if ( data->ok ) { /* Stream temp files to stdout: */
    double* buffer = data->buffer;
    const size_t bufferBytes = data->rows * data->columns * 4 * sizeof *buffer;
    streamHeader( data );

    for ( file = 0; data->ok && file < files; ++file ) {

      while ( data->ok && ! feof( data->tempFiles[ file ] ) ) {
        const size_t bytesRead =
          fread( buffer, 1, bufferBytes, data->tempFiles[ file ] );

        if ( bytesRead ) {
          const size_t bytesWritten = fwrite( buffer, 1, bytesRead, stdout );
          data->ok = bytesWritten == bytesRead;
        }
      }

      if ( ! data->ok ) {
        fprintf( stderr,
                 "\nFailed to stream subset data from temp file '%s'.\n",
                 data->tempFileNames[ file ] );
      }
    }
  }
}



/******************************************************************************
PURPOSE: streamHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
******************************************************************************/

static void streamHeader( const Data* const data ) {
  const Arguments* const arguments = &data->arguments;
  const int yyyymmddhh = arguments->yyyymmddhh;
  const int yyyy       = yyyymmddhh / 1000000;
  const int mm         = yyyymmddhh / 10000 % 100;
  const int dd         = yyyymmddhh / 100 % 100;
  const int hh         = yyyymmddhh % 100;
  const size_t subsetRows =
    1 +
    data->subsetIndices[ ROW ][ MAXIMUM ] -
    data->subsetIndices[ ROW ][ MINIMUM ];
  const size_t subsetColumns =
    1 +
    data->subsetIndices[ COLUMN ][ MAXIMUM ] -
    data->subsetIndices[ COLUMN ][ MINIMUM ];

  if ( arguments->cmaq ) {
    printf( "SUBSET 9.0 CMAQ\nHRRR\n%s\n"
            "%04d-%02d-%02dT%02d:00:00-0000\n",
            arguments->description, yyyy, mm, dd, hh );

    printf( "# data dimensions: timesteps variables layers rows columns:\n"
            "%d %d 1 %lu %lu\n",
            arguments->hours, 1 + arguments->isVector2,
            subsetRows, subsetColumns );

    printf( "# subset indices (0-based time, 1-based layer/row/column): "
            "first-timestep last-timestep "
            "first-layer last-layer "
            "first-row last-row "
            "first-column last-column:\n"
            "0 %d 1 1 %lu %lu %lu %lu\n",
            arguments->hours - 1,
            1 + data->subsetIndices[ ROW ][ MINIMUM ],
            1 + data->subsetIndices[ ROW ][ MAXIMUM ],
            1 + data->subsetIndices[ COLUMN ][ MINIMUM ],
            1 + data->subsetIndices[ COLUMN ][ MAXIMUM ] );

    if ( arguments->isVector2 ) {
      printf( "# Variable names:\n%s_u %s_v\n",
              arguments->variable, arguments->variable );
      printf( "# Variable units:\n%s %s\n", arguments->units, arguments->units);
    } else {
      printf( "# Variable names:\n%s\n", arguments->variable );
      printf( "# Variable units:\n%s\n", arguments->units );
    }

    printf( "# lcc projection: "
            "lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis\n"
            "38.5 38.5 38.5 -97.5 6370000 6370000\n" );

    printf( "# Grid: ncols nrows xorig yorig xcell ycell "
            "vgtyp vgtop vglvls[2]:\n"
            "%lu %lu "
            "-2698552.865923 -1588499.88061 3000 3000 2 10000 1 0.995\n",
            data->columns, data->rows );

    printf( "# IEEE-754 32-bit reals "
            "data[variables][timesteps][layers][rows][columns]:\n" );
  } else {
    printf( "Grid 1.0\n%s\n%04d-%02d-%02dT%02d:00:00-0000\n",
            arguments->description, yyyy, mm, dd, hh );
    printf( "# Dimensions: timesteps variables rows columns:\n%d %d %lu %lu\n",
            arguments->hours, 1 + arguments->isVector2,
            subsetRows, subsetColumns );

    if ( arguments->isVector2 ) {
      printf( "# Variable names:\n%s_u %s_v\n",
              arguments->variable, arguments->variable );
      printf( "# Variable units:\n%s %s\n", arguments->units, arguments->units);
    } else {
      printf( "# Variable names:\n%s\n", arguments->variable );
      printf( "# Variable units:\n%s\n", arguments->units );
    }

    printf( "# IEEE-754 64-bit reals longitudes[rows][columns] and\n" );
    printf( "# IEEE-754 64-bit reals latitudes[rows][columns] and\n" );
    printf( "# IEEE-754 64-bit reals "
            "data[timesteps][variables][rows][columns]:\n" );
  }
}




/******************************************************************************
PURPOSE: OpenAQSubset.c - Extract a lon-lat subset of data from a list of
         OpenAQ files and write it to stdout as XDR binary format.

NOTES:   This program calls system( "/usr/bin/sort ..." ).

         Compile:
         gcc -Wall -g -o OpenAQSubset OpenAQSubset.c Utilities.c \
                   -L../../../lib/$platform \
                   -lz -lm -lc

         Usage:
         OpenAQSubset \
           [-files | -info] <listfile> \
           [-tmpdir <temp_directory>] (default .)\
           [-desc "description text"] (default URL) \
           -timerange <yyyymmddhh> <yyyymmddhh> \
           -variable <name> \
           [-bounds <minimum_longitude> <minimum_latitude> \
                    <maximum_longitude> <maximum_latitude>] \
           [-format ascii | xdr] (default xdr) \
           [-sensor_type all | \
                         ABCD | Aernode | AIO | AirGradient | BAM | BEACO2N | \
                         Clarity | Government | HabitatMap | Miri | RAMP \
                         Senstate] (default is all)
           [-aggregate none | hourly | daily | all] (default is hourly)

          Example:
          ../../../bin/$platform/OpenAQSubset \
          -files testdata/file_list \
          -tmpdir testdata \
          -desc "https://docs.openaq.org/aws/about,OpenAQSubset" \
          -timerange 2026060100 2026060023 \
          -variable pm25 \
          -bounds -75 41 -74 42 \
          -format xdr \
          -aggregate hourly \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 10-line ASCII header followed by binary 64-bit big-endian arrays:

Point 1.0
https://docs.openaq.org/aws/about,OpenAQSubset
2026-06-01T00:00:00-0000 2026-06-02T23:59:59-0000
# Dimensions: variables points
5 459
# Variable names:
timestamp longitude latitude id pm25
# Variable units:
yyyymmddhhmmss deg deg - ug/m3
# char notes[points][80] and
# IEEE-754 64-bit reals data[variables][points]:
<big-endian binary format array>

HISTORY: 2026-06-06 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/


/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, sscanf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <stdlib.h>    /* For malloc(), free(), atoll(), atof(), system(). */
#include <limits.h>    /* For INT_MAX, LLONG_MAX. */
#include <float.h>     /* For DBL_MAX. */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For Bounds, readFile(), readCompressedFile(). */

/*================================= MACROS =================================*/


/* Name of temp per-variable files created in -tmpdir with PID appended: */

#define TEMP_FILE_NAME "junk_OpenAQSubset"

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

/*================================ CONSTANTS ================================*/


/*---------------------------------------------------------------------------*/

/* Output vars: timestamp, longitude, latitude, id, pm25: */

enum {
  TIMESTAMP_INDEX,
  LONGITUDE_INDEX,
  LATITUDE_INDEX,
  ID_INDEX,
  VARIABLE_INDEX,
  VARIABLES
};

enum { TEMP_FILE_1 = VARIABLES, TEMP_FILE_2, TEMP_AGGREGATED_FILE, TEMP_FILES};

static const char* const tempFileNames[ TEMP_FILES ] = {
  "timestamp", "longitude", "latitude", "id", "data",
  "temp1", "temp2", "aggregated"
};

enum { FORMAT_XDR, FORMAT_ASCII };
#define FORMAT_STRING "xdr ascii"

enum { SENSOR_TYPE_ABCD, SENSOR_TYPE_AERNODE, SENSOR_TYPE_AIO,
  SENSOR_TYPE_AIRGRADIENT, SENSOR_TYPE_BAM, SENSOR_TYPE_BEACO2N,
  SENSOR_TYPE_CLARITY, SENSOR_TYPE_GOVERNMENT, SENSOR_TYPE_HABITATMAP,
  SENSOR_TYPE_MIRI, SENSOR_TYPE_RAMP, SENSOR_TYPE_SENSTATE
};
#define SENSOR_TYPE_STRING "ABCD Aernode AIO AirGradient BAM BEACO2N Clarity Government HabitatMap Miri RAMP, Senstate"

enum {
  AGGREGATE_NONE, AGGREGATE_ALL, AGGREGATE_DAILY, AGGREGATE_HOURLY
};
#define AGGREGATE_STRING "none all daily hourly"



/*================================== TYPES ==================================*/


/* List of known variables: */

typedef struct {
  const char* name;          /* Name of output variable. "pm25". */
  const char* units;         /* Units of output variable. "ug/m3". */
  const double validMinimum; /* Minimum valid value. */
  const double validMaximum; /* Maximum valid value. */
} ColumnInfo;

static ColumnInfo columnInfo[] = {
  { "id",                "-",               1, LLONG_MAX },
  { "description",       "-",               0, 0 },
  { "timestamp",         "UTC",             0, 0 },
  { "longitude",         "deg",             -180.0, 180.0 },
  { "latitude",          "deg",             -90.0, 90.0 },

  { "bc",                "ug/m3",           0, 1e6 },
  { "bc_370",            "ng/m3",           0, 1e9 },
  { "bc_375",            "ng/m3",           0, 1e9 },
  { "bc_470",            "ng/m3",           0, 1e9 },
  { "bc_528",            "ng/m3",           0, 1e9 },
  { "bc_625",            "ng/m3",           0, 1e9 },
  { "bc_880",            "ng/m3",           0, 1e9 },
  { "co",                "ppm",             0, 1e6 },
  { "co2",               "ppm",             0, 1e6 },
  { "humidity",          "%",               0, 100 },
  { "no",                "ppb",             0, 1e6 },
  { "no2",               "ppb",             0, 1e6 },
  { "nox",               "ppb",             0, 1e6 },
  { "o3",                "ppb",             0, 1e6 },
  { "pm1",               "ug/m3",           0, 1e6 },
  { "pm4",               "ug/m3",           0, 1e6 },
  { "pm10",              "ug/m3",           0, 1e6 },
  { "pm25",              "ug/m3",           0, 1e6,},
  { "pressure",          "hpa",             800, 1500 },
  { "so2",               "ppb",             0, 1e6 },
  { "temperature",       "C",               -50, 100 },
  { "um003",             "particles/cm3",   0, 1e6 },
  { "wind_direction",    "deg",             0, 360 },
  { "wind_speed",        "m/s",             0, 50 }
};

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
  long long   yyyymmddhh[ 2 ]; /* Beginning/ending timestamp of subset. */
  int         format;          /* FORMAT_XDR, FORMAT_ASCII. */
  const char* sensorType;      /* 0 or "Clarity", "Government", etc. */
  int         aggregate;       /* AGGREGATE_NONE/ALL/DAILY/HOURLY. */
} Arguments;

/* Data type: */

typedef struct {
  Arguments   arguments;     /* User-supplied (command-line) arguments. */
  const char* units;         /* Units of output variable. E.g., "ug/m3". */
  double validMinimum;       /* Minimum valid value for variable. */
  double validMaximum;       /* Maximum valid value for variable. */
  FileName    tempFileNames[ TEMP_FILES ]; /* Pathed names of temp files. */
  FILE*       tempFiles[ TEMP_FILES ]; /* Temp files of output subset data. */
  size_t      bufferSize;    /* Byte size of input/output buffer. */
  size_t      inputLength;   /* Length of input buffer. */
  char*       inputBuffer;   /* Holds current input file content. */
  char*       outputBuffer;  /* Holds an output subset data content. */
  double      columnValues[ VARIABLES ]; /* Value of required columns. */
  Note        note;          /* Column description. */
  size_t      points;        /* Number of valid data points in subset. */
  int         ok;            /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static int printInfo( Data* const data, const char* const listFile );

static int updateInfo( char* content,
                       int* id,
                       double* longitude,
                       double* latitude,
                       int* yyyymmdd1,
                       int* yyyymmdd2,
                       const int variablesSize,
                       char* variables,
                       Note description );

static int parseColumnInfo( char* line,
                            int* id,
                            double* longitude,
                            double* latitude,
                            int* yyyymmdd1,
                            int* yyyymmdd2,
                            const int variablesSize,
                            char* variables,
                            Note description );

static void removeTempFiles( Data* const data );

static void createTempFiles( Data* const data );

static void closeTempFiles( Data* const data );

static void openVariableTempFiles( Data* const data );

static void allocateBuffers( Data* const data );

static void deallocateBuffers( Data* const data );

static void reallocateOutputBuffer( Data* const data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static void setVariableUnitsAndRange( Data* const data );

static int unitsMatch( const char* const string, const char* const units,
                       double* scale);

static void readData( Data* const data );

static void extractSubset( Data* data );

static int parseColumnValues( char* line,
                              const long long yyyymmddhh1,
                              const long long yyyymmddhh2,
                              const double west,
                              const double east,
                              const double south,
                              const double north,
                              const char* const variable,
                              const char* const units,
                              const double validMinimum,
                              const double validMaximum,
                              double columnValues[],
                              Note note );

static void aggregateData( Data* data );

static void reformatData( Data* data );

static void writeTempVariableFiles( Data* data );

static void sortTempData( Data* data );

static void streamData( Data* const data );

static void streamHeader( const Data* const data );


/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of OpenAQ files and write
         it to stdout in XDR format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Data data;
  memset( &data, 0, sizeof data );

  if ( argc == 3 && ! strcmp( argv[ 1 ], "-info" ) ) {
    data.ok = fileSize( argv[ 2 ] ) > 0;

    if ( data.ok ) {
      allocateBuffers( &data );
      data.points = printInfo( &data, argv[ 2 ] );
      deallocateBuffers( &data );
    }
  } else {
    data.ok = parseArguments( argc, argv, &data.arguments );

    if ( ! data.ok ) {
      printUsage( argv[ 0 ] );
    } else {
      setVariableUnitsAndRange( &data );
      createTempFiles( &data );

      if ( data.ok ) {
        allocateBuffers( &data );

        if ( data.ok ) {
          readData( &data ); /* Read input data and write temp files. */

          if ( data.ok && data.points > 0 ) {
            reformatData( &data );

            if ( data.ok ) {
              streamData( &data ); /* Write header & temp files to stdout. */
            }
          }
        }

        deallocateBuffers( &data );
      }

      removeTempFiles( &data );

      if ( ! data.ok ) {
        fprintf( stderr, "\n%s: No points were in the subset.\n", argv[ 0 ] );
      }
    }
  }

  data.ok = data.ok && data.points > 0;
  return ! data.ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: printInfo - Read all listed data files and print a site info line.
INPUTS:  Data* data                  Data with allocated inputBuffer.
         const char* const listFile  List of data files to read.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int printInfo( Data* const data, const char* const listFile ) {
  int parsed = 0;
  size_t length = 0;
  char* listFileContent = 0;
  int ok = readFile( listFile, &length, &listFileContent );

  if ( ok ) {
    int id = 0;
    double longitude = 0.0;
    double latitude = 0.0;
    int yyyymmdd1 = 0;
    int yyyymmdd2 = 0;
    enum { VARIALBES_SIZE = 4096 };
    char variables[ VARIALBES_SIZE ] = "";
    Note description = "";
    char* inputFileName = 0;
    char* end = 0;
    memset( variables, 0, sizeof variables );
    memset( description, 0, sizeof description );

    /* Get each line of list file. It is the data file to read: */

    for ( inputFileName = strtok_r( listFileContent, "\n", &end );
          inputFileName;
          inputFileName = strtok_r( 0, "\n", &end ) ) {

      const size_t contentLength =
        readCompressedFile( inputFileName,
                            &data->inputLength, &data->inputBuffer );

      ok = contentLength > 0;

      if ( ! ok ) {
        fprintf( stderr, "\nOmitting invalid file %s\n", inputFileName );
      } else {
        eraseQuotedCommasAndQuotesAndBrackets( data->inputBuffer ); /*hotspot*/
        ok = updateInfo( data->inputBuffer,
                         &id, &longitude, &latitude, &yyyymmdd1, &yyyymmdd2,
                         VARIALBES_SIZE, variables, description );
        if ( ok ) {
          parsed = 1;
        }
      }
    } /* End loop on listFile. */

    if ( parsed ) {
      printf( "%d\t%f\t%f\t%d\t%d\t%s\t%s\n",
              id, longitude, latitude, yyyymmdd1, yyyymmdd2,
              variables, description );
    }

    free( listFileContent ), listFileContent = 0;
  }

  return parsed;
}




/******************************************************************************
PURPOSE: updateInfo - Parse required column data from inputBuffer and update
         info.
INPUTS:  char* content  Content of data file.
OUTPUTS: int* id             Site id (initialize if 0 else match).
         double* longitude   Longitude of site.
         double* latitude    Latitude of site.
         int* yyyymmdd1      Earliest date of data.
         int* yyyymmdd2      Latest date of data.
         int variablesSize   Size of variables array.
         char* variables     Space-delimited list of each variable(unit).
         Note description    Site description (filtered of problematic chars).
RETURNS: 1 if successful, else 0.
******************************************************************************/

static int updateInfo( char* content,
                       int* id,
                       double* longitude,
                       double* latitude,
                       int* yyyymmdd1,
                       int* yyyymmdd2,
                       const int variablesSize,
                       char* variables,
                       Note description ) {
  int result = 0;
  assert( content );
  assert( id );
  assert( longitude );
  assert( latitude );
  assert( yyyymmdd1 );
  assert( yyyymmdd2 );
  assert( variablesSize );
  assert( variables );
  assert( description );

  {
    char* line = 0;
    char* end = 0;
    int skip = 1;

    /* Get each line but skip the first header line */

    for ( line = strtok_r( content, "\n", &end );
          line;
          line = strtok_r( 0, "\n", &end ) ) {

      if ( skip ) {
        skip = 0;
      } else {
        const int ok =
          parseColumnInfo( line,
                           id, longitude, latitude, yyyymmdd1, yyyymmdd2,
                           variablesSize, variables, description );

        if ( ok ) {
          result = 1;
        }
      }
    }
  }

  return result;
}


/******************************************************************************
PURPOSE: parseColumnInfo - Parse column info from a line and update values.
INPUTS:  char* line          Data line to parse.
OUTPUTS: int* id             Site id to initialize or match.
         double* longitude   Longitude of site.
         double* latitude    Latitude of site.
         int* yyyymmdd1      Earliest date of data.
         int* yyyymmdd2      Latest date of data.
         int variablesSize   Size of variables array.
         char* variables     Space-delimited list of each variable(unit).
         Note description     Site description (filtered of problematic chars).
RETURNS: int 1 successful, else 0.
******************************************************************************/

static int parseColumnInfo( char* line,
                            int* id,
                            double* longitude,
                            double* latitude,
                            int* yyyymmdd1,
                            int* yyyymmdd2,
                            const int variablesSize,
                            char* variables,
                            Note description ) {
  int ok = 0;

  assert( line );
  assert( id );
  assert( longitude );
  assert( latitude );
  assert( yyyymmdd1 );
  assert( yyyymmdd2 );
  assert( variablesSize > sizeof (Note) );
  assert( variables );
  assert( description );

  {
    char* value = 0;
    char* endValue = 0;
    int column = 0;
    Note parsedVariable = "";
    memset( parsedVariable, 0, sizeof parsedVariable );
    memset( description, 0, sizeof (Note) );
    description[ 0 ] = '-';

    for ( value = strtok_r( line, ",", &endValue ), column = 0;
          value;
          value = strtok_r( 0, ",", &endValue ), ++column ) {

      if ( *value ) {

        if ( column == 0 ) {
          long long parsedId = 0;
          ok = parseLongLong( value, 1, LLONG_MAX, &parsedId );

          if ( ok ) {

            if ( *id == 0 ) {
              *id = parsedId;
            } else {
              ok = parsedId == *id;

              if ( ! ok ) {
                break;
              }
            }
          }
        } else if ( column == 2 ) {
          strncpy( description, value, sizeof (Note) / sizeof (char) - 1 );
          filterNotes( description, sizeof (Note) ); /* Erase problem chars! */
        } else if ( column == 3 ) {
          long long timestamp = 0;
          ok = parseTimestamp( value, 2000010100, 2147123123, &timestamp );

          if ( ok ) {
            timestamp /= 100; /* Remove hh. */

            if ( *yyyymmdd1 == 0 ) {
              *yyyymmdd1 = *yyyymmdd2 = timestamp;
            } else if ( timestamp < *yyyymmdd1 ) {
              *yyyymmdd1 = timestamp;
            } else if ( timestamp > *yyyymmdd2 ) {
              *yyyymmdd2 = timestamp;
            }
          }
        } else if ( column == 4 ) {
          double parsedLatitude = 0.0;
          ok = parseDouble( value, -90.0, 90.0, &parsedLatitude );

          if ( ok ) {
            *latitude = parsedLatitude;
          }
        } else if ( column == 5 ) {
          double parsedLongitude = 0.0;
          ok = parseDouble( value, -180.0, 180.0, &parsedLongitude );

          if ( ok ) {
            *longitude = parsedLongitude;
          }
        } else if ( column == 6 ) {
          char* space = strrchr( value, ' ' );
          memset( parsedVariable, 0, sizeof parsedVariable );

          /* Trim any leading or trailing space: */

          if ( space ) {
            *space = '\0';

            if ( *value == ' ' ) {
              ++value;
            }
          }

          strncpy( parsedVariable, value, sizeof parsedVariable );
        } else if ( column == 7 ) {
          char* space = strrchr( value, ' ' );

          /* Trim any leading or trailing space: */

          if ( space ) {
            *space = '\0';

            if ( *value == ' ' ) {
              ++value;
            }
          }

          strncat( parsedVariable, "(",
                   sizeof parsedVariable - strlen( parsedVariable ) - 1 );
          strncat( parsedVariable, value,
                   sizeof parsedVariable - strlen( parsedVariable ) - 1 );
          strncat( parsedVariable, ")",
                   sizeof parsedVariable - strlen( parsedVariable ) - 1 );

          if ( ! strstr( variables, parsedVariable ) ) {
            strncat( variables, parsedVariable,
                     variablesSize - strlen( variables ) - 2 );
            strncat( variables, " ",
                     variablesSize - strlen( variables ) - 2 );
          }
        }
      }
    }
  }

  ok =
    *id > 0 &&
    IN_RANGE( *longitude, -180.0, 180.0 ) &&
    IN_RANGE( *latitude, -90.0, 90.0 ) &&
    isValidYYYYMMDD( *yyyymmdd1 ) &&
    isValidYYYYMMDD( *yyyymmdd2 ) &&
    *yyyymmdd2 >= *yyyymmdd1 &&
    strstr( variables, ")" ) &&
    *description;

  return ok;
}



/******************************************************************************
PURPOSE: setVariableUnitsAndRange - Set variable units and valid range.
INPUTS:  Data* const data  data->arguments->variable
OUTPUTS: Data* const data  data->units, data->validMinimum, data->validMaximum.
******************************************************************************/

static void setVariableUnitsAndRange( Data* const data ) {
  const size_t count = sizeof columnInfo / sizeof *columnInfo;
  const char* const variableName = data->arguments.variable;
  size_t variableIndex = 0;

  for ( ; variableIndex < count; ++variableIndex ) {
    const ColumnInfo* const entry = columnInfo + variableIndex;

    if ( ! strcmp( entry->name, variableName ) ) {
      data->units = entry->units;
      data->validMinimum = entry->validMinimum;
      data->validMaximum = entry->validMaximum;
      variableIndex = count;
    }
  }
}



/******************************************************************************
PURPOSE: unitsMatch - Check if string matches units.
INPUTS:  const char* const string      String to parse.
         const char* const units       Units to match.
         double* scale                 Scale factor for units conversion.
RETURNS: int 1 if matched, else 0.
NOTES: To see the ASCII code for the problematic extended ASCII characters in
the data file, copy and echo the units string | od -b
\302 \265 = mu character
\302 \263 = superscript 3
******************************************************************************/

static int unitsMatch( const char* const string, const char* const units,
                       double* scale ) {
  int result = 0;
  assert( string ); assert( units ); assert( *units ); assert( scale );
  *scale = 1.0;

  if ( ! strcmp( units, "ug/m3" ) ) {
    result = string[ 0 ] == ' ' &&
             string[ 1 ] == '\302' &&
             string[ 2 ] == '\265' &&
             string[ 3 ] == 'g' &&
             string[ 4 ] == '/' &&
             string[ 5 ] == 'm' &&
             string[ 6 ] == '\302' &&
             string[ 7 ] == '\263' &&
             string[ 8 ] == ' ';
  } else if ( ! strcmp( units, "ng/m3" ) ) {
    result = string[ 0 ] == ' ' &&
             string[ 1 ] == 'n' &&
             string[ 2 ] == 'g' &&
             string[ 3 ] == '/' &&
             string[ 4 ] == 'm' &&
             string[ 5 ] == '\302' &&
             string[ 6 ] == '\263' &&
             string[ 7 ] == ' ';
  } else if ( ! strcmp( units, "ppm" ) ) {
    result = ! strcmp( string, " ppm " ) || ! strcmp( string, "ppm" );

    if ( ! result ) {
      result = ! strcmp( string, " ppb " ) || ! strcmp( string, "ppb" );

      if ( result ) {
        *scale = 0.001; /* ppb to ppm. */
      }
    }
  } else if ( ! strcmp( units, "ppb" ) ) {
    result = ! strcmp( string, " ppb " ) || ! strcmp( string, "ppb" );

    if ( ! result ) {
      result = ! strcmp( string, " ppm " ) || ! strcmp( string, "ppm" );

      if ( result ) {
        *scale = 1000.0; /* ppm to ppb. */
      }
    }
  } else if ( ! strcmp( units, "%" ) ) {
    result = ! strcmp( string, " % " );
  } else if ( ! strcmp( units, "hpa" ) ) {
    result = ! strcmp( string, " hpa " );
  } else if ( ! strcmp( units, "C" ) ) {
    result = ! strcmp( string, " c " );
  } else if ( ! strcmp( units, "particles/cm3" ) ) {
    result = strlen( string ) == 16 &&
             strstr( string, " particles/cm" ) == string &&
             string[ 13 ] == '\302' &&
             string[ 14 ] == '\263' &&
             string[ 15 ] == ' ';
  } else if ( ! strcmp( units, "deg" ) ) {
    result = ! strcmp( string, " deg " );
  } else if ( ! strcmp( units, "m/s" ) ) {
    result = ! strcmp( string, " m/s " );
  }

  return result;
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
PURPOSE: openVariableTempFiles - Open temp variable output files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->ok, data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void openVariableTempFiles( Data* const data ) {
  size_t index = 0;
  assert( data );

  for ( index = 0; data->ok && index < VARIABLES; ++index ) {
    assert( data->tempFileNames[ index ] );
    assert( data->tempFileNames[ index ][ 0 ] );
    data->tempFiles[ index ] = fopen( data->tempFileNames[ index ], "wb" );

    if ( ! data->tempFiles[ index ] ) {
      fprintf( stderr, "\nCan't create temporary output file '%s'.\n",
               data->tempFileNames[ index ] );
      data->ok = 0;
    }
  }
}



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
PURPOSE: allocateBuffers - Allocate input/output buffers.
OUTPUTS: Data* const data  data->inputLength  Allocated length of input buffer.
                           data->inputBuffer  Allocated input buffer.
                           data->outputBuffer Allocated output buffer.
                           data->ok 1 if successful, else 0.
******************************************************************************/

static void allocateBuffers( Data* const data ) {
  assert( data );
  assert( data->inputBuffer  == 0 );
  assert( data->outputBuffer == 0 );
  assert( data->inputLength == 0 );
  data->bufferSize = 100 * 1024 * sizeof (char);
  data->inputBuffer = malloc( data->bufferSize );
  data->outputBuffer = data->inputBuffer ? malloc( data->bufferSize ) : 0;
  data->ok = data->outputBuffer != 0;

  if ( ! data->ok ) {
      fprintf( stderr,
               "\nCan't allocate %lu bytes "
               "to complete the requested action.\n", 2 * data->bufferSize );

    if ( data->inputBuffer ) {
      free( data->inputBuffer ), data->inputBuffer = 0;
    }

    data->bufferSize = 0;
  } else {
    memset( data->inputBuffer,  0, data->bufferSize );
    memset( data->outputBuffer, 0, data->bufferSize );
    data->inputLength = data->bufferSize;
    DEBUG( fprintf( stderr, "allocateBuffers(): bufferSize = %lu\n",
                    data->bufferSize ); )
  }
}



/******************************************************************************
PURPOSE: deallocateBuffers - Deallocate input/output buffers.
OUTPUTS: Data* const data  data->inputLength  0.
                           data->inputBuffer  0.
                           data->outputBuffer 0.
******************************************************************************/

static void deallocateBuffers( Data* const data ) {
  assert( data );

  if ( data->inputBuffer ) {
    free( data->inputBuffer ), data->inputBuffer = 0;
  }

  if ( data->outputBuffer ) {
    free( data->outputBuffer ), data->outputBuffer = 0;
  }

  data->bufferSize = 0;
  data->inputLength = 0;

  assert( data->inputBuffer == 0 );
  assert( data->outputBuffer == 0 );
  assert( data->bufferSize  == 0 );
  assert( data->inputLength == 0 );
}



/******************************************************************************
PURPOSE: reallocateOutputBuffer - Reallocate output buffer.
INPUTS:  Data* const data  data->inputLength  Allocated length of input buffer.
                           data->outputBuffer Allocated output buffer.
OUTPUTS: Data* const data  data->bufferSize   Size of output buffer.
                           data->outputBuffer Allocated output buffer.
                           data->ok 1 if successful, else 0.
******************************************************************************/

static void reallocateOutputBuffer( Data* const data ) {
  assert( data );
  assert( data->inputBuffer );
  assert( data->inputLength * sizeof (char) > data->bufferSize );

  if ( data->outputBuffer ) {
    free( data->outputBuffer ), data->outputBuffer = 0;
  }

  data->bufferSize = data->inputLength * sizeof (char);
  data->outputBuffer = malloc( data->bufferSize );
  data->ok = data->outputBuffer != 0;

  if ( ! data->ok ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes "
             "to complete the requested action.\n", data->bufferSize );
    data->bufferSize = 0;
  } else {
    memset( data->outputBuffer, 0, data->bufferSize );

    DEBUG( fprintf( stderr, "reallocateOutputBuffer(): bufferSize = %lu\n",
                    data->bufferSize ); )

  }

  assert( (data->ok == 0 && data->outputBuffer == 0 && data->bufferSize == 0) ||
          (data->ok == 1 && data->outputBuffer != 0 &&
            data->bufferSize == data->inputLength * sizeof (char) ) );
}



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* name  Name of program.
******************************************************************************/

static void printUsage( const char* name ) {
  assert( name ); assert( *name );
  fprintf( stderr,
           "\n%s - Extract a subset of data from a time-sorted list of\n"
          "OpenAQ files and write it to stdout in XDR binary format.\n",
           name );
  fprintf( stderr, "Data is subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "  [-files | -info] <listfile> \\\n" );
  fprintf( stderr, "  [-tmpdir <temp_directory>] (default .) \\\n" );
  fprintf( stderr, "  [-desc \"description text\"] (default URL) \\\n" );
  fprintf( stderr, "  -timerange <yyyymmddhh> <yyyymmddhh> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, " [-bounds <minimum_longitude> <minimum_latitude>]" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, " [-format ascii | xdr] (default xdr) \\\n" );
  fprintf( stderr, " [-sensor_type all" );
  fprintf( stderr, " | ABCD | Aernode | AIO | AirGradient" );
  fprintf( stderr, " | BAM | BEACO2N | Clarity" );
  fprintf( stderr, " | Government | HabitatMap | Miri" );
  fprintf( stderr, " | RAMP | Senstate] (default is all)\n" );
  fprintf( stderr, " [-aggregate none | hourly | daily | all] " );
  fprintf( stderr, "(default = hourly) \\\n" );
  fprintf( stderr, "Note:\ntimes are in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory were temp files are " );
  fprintf( stderr, "written.\nIt should have enough disk space (1TB).\n" );
  fprintf( stderr, "-info is used to print a summary site info line.\n" );
  fprintf( stderr, "\nExample:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-files file_list \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc "
                   "\"https://docs.openaq.org/aws/about,OpenAQSubset\" \\\n");
  fprintf( stderr, "-timerange 2026060100 2026060223 \\\n" );
  fprintf( stderr, "-variable pm25 \\\n" );
  fprintf( stderr, "-bounds -75 40 -74 41 \\\n" );
  fprintf( stderr, "-format xdr \\\n" );
  fprintf( stderr, "> subset.xdr\n\n" );
  fprintf( stderr, "Hourly corrected PM2.5 over NYC on June 1-2, 2026.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Point 1.0\n" );
  fprintf( stderr, "https://docs.openaq.org/aws/about,OpenAQSubset\n" );
  fprintf( stderr, "2026-06-01T00:00:00-0000 2026-06-02T23:59:59-0000\n" );
  fprintf( stderr, "# Dimensions: variables points\n" );
  fprintf( stderr, "6 20\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "timestamp longitude latitude id pm25\n");
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss deg deg - ug/m3\n" );
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
  int result = argc >= 8;
  static Option options[] = {
    { "-files",        1, FILE_TYPE,        1, 0, 0, 0, 0 },
    { "-tmpdir",       0, DIRECTORY_TYPE,   1, 0, 0, 0, 0 },
    { "-desc",         0, STRING_TYPE,      1, 0, 0, 0, 0 },
    { "-variable",     1, ENUM_TYPE,        1, 0, 0, 0, 0 },
    { "-timerange",    1, YYYYMMDDHH_TYPE,  2, 0, 0, 0, 0 },
    { "-format",       0, ENUM_TYPE,        1, 0, FORMAT_STRING, 0, 0 },
    { "-bounds",       0, BOUNDS_TYPE,      4, 0, 0, 0, 0 },
    { "-sensor_type",  0, STRING_TYPE,      1, 0, SENSOR_TYPE_STRING, 0, 0 },
    { "-aggregate",    0, ENUM_TYPE,        1, 0, AGGREGATE_STRING, 0, 0 }
  };
  char variableNames[ sizeof columnInfo / sizeof columnInfo[ 0 ] * 32 ] = "";
  int variableIndex = -1;

  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] ); assert( arguments );

  if ( result ) {

    /* Initialize arguments to defaults: */

    memset( arguments, 0, sizeof (Arguments) );
    arguments->tmpdir = ".";
    arguments->description = "https://docs.openaq.org/aws/about,OpenAQSubset";
    arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
    arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
    arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
    arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
    arguments->aggregate = AGGREGATE_HOURLY;

    /* Init space-delimited string of valid variable names for parse check: */

    {
      const size_t variables = sizeof columnInfo / sizeof columnInfo[ 0 ];
      size_t variable = 0;
      memset( variableNames, 0, sizeof variableNames );

      for ( variable = 0; variable < variables; ++variable ) {
        assert( strlen( variableNames ) + 32
                < sizeof variableNames / sizeof *variableNames );
        strncat( variableNames, columnInfo[ variable ].name, 30 );
        strcat( variableNames, " " );
      }
    }

    /* Finish initializing non-compile-time-constant parts of options: */

    options[ 0 ].values = &arguments->listFile;
    options[ 1 ].values = &arguments->tmpdir;
    options[ 2 ].values = &arguments->description;
    options[ 3 ].valids = variableNames; /* Space-delimited var names. */
    options[ 3 ].values = &variableIndex; /* Temp int to get index of name. */
    options[ 4 ].values = &arguments->yyyymmddhh[ 0 ];
    options[ 5 ].values = &arguments->format;
    options[ 6 ].values = &arguments->bounds[ 0 ][ 0 ];
    options[ 7 ].values = &arguments->sensorType;
    options[ 8 ].values = &arguments->aggregate;

    result =
      parseOptions( argc, argv, sizeof options / sizeof *options, options );

    if ( result ) { /* Initialize argument->variable to static string name: */
      arguments->variable = columnInfo[ variableIndex ].name;
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
           isValidYYYYMMDDHH( arguments->yyyymmddhh[ 0 ] ) &&
           isValidYYYYMMDDHH( arguments->yyyymmddhh[ 1 ] ) &&
           arguments->yyyymmddhh[ 0 ] <= arguments->yyyymmddhh[ 1 ] &&
           isValidBounds( (const double (*)[2]) arguments->bounds ) ) );

  return result;
}



/******************************************************************************
PURPOSE: readData - Read data from each listed data file and
         write the subset of data to the temporary files.
INPUTS:  Data* data  Data to read.
OUTPUTS: Data* data  data->points, ok, tempFiles[] = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  const int aggregate = arguments->aggregate;
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

      const size_t contentLength =
        readCompressedFile( inputFileName,
                            &data->inputLength, &data->inputBuffer );

      data->ok = contentLength > 0;

      DEBUG( fprintf( stderr,
                      "%s, contentLength = %lu, data->inputLength = %lu\n",
                      inputFileName, contentLength, data->inputLength ); )

      if ( data->ok ) {

        if ( data->inputLength * sizeof (char) > data->bufferSize ) {
          reallocateOutputBuffer( data );
        }

        if ( data->ok ) {
          extractSubset( data ); /* Appends to temp file 1. */

          if ( data->ok ) {
            wroteSomeData = 1;
          } else {
            fprintf( stderr, "\nOmitting invalid file %s\n", inputFileName );
          }
        }
      }

      data->ok = 1; /* Ignore invalid files. */
    } /* End loop on listFile. */

    free( listFileContent ), listFileContent = 0;
    data->ok = wroteSomeData;

    if ( wroteSomeData ) {

      /* Sort temp file 1 into temp file 2 and empty temp file 1: */

      sortTempData( data );

      if ( data->ok && aggregate != AGGREGATE_NONE ) {

        /* Aggregate sorted temp file 2 appending onto temp aggregate file: */

        aggregateData( data );
        DEBUG( fprintf( stderr, "after final aggregateData() data->ok = %d\n",
                        data->ok ); )
      }

      closeTempFiles( data );
    }
  }

  DEBUG( fprintf( stderr,
                  "End of readData() data->points = %lu, data->ok = %d\n",
                  data->points, data->ok ); )
}



/******************************************************************************
PURPOSE: extractSubset - Parse required column data from inputBuffer and write
         to outputBuffer then to temp file 1.
INPUTS:  Data* data  data->columnIndices[]  0-based indices of required columns
                     data->inputBuffer Content of input data file.
OUTPUTS: Data* data  data->outputBuffer  Content of extract.
                     data->tempFiles[ TEMP_FILE_1 ]  Write outputBuffer to it.
                     data->points  Increased if AGGREGSTE_NONE, else unchanged.
                     data->ok = 1 if successful, else 0.
NOTES:   For performance, both the input data file and the subset output are
         written to internal string buffers and the output buffer is written
         to temp file 1 at the end of extracting subset for the current file.
******************************************************************************/

static void extractSubset( Data* data ) {
  assert( data );
  assert( data->arguments.variable ); assert( data->arguments.variable[ 0 ] );
  assert( data->inputBuffer );
  assert( data->tempFiles[ TEMP_FILE_1 ] );

  {
    const Arguments* const arguments = &data->arguments;
    const char* const variable = arguments->variable;
    const long long yyyymmddhh1 = arguments->yyyymmddhh[ 0 ];
    const long long yyyymmddhh2 = arguments->yyyymmddhh[ 1 ];
    const double west  = arguments->bounds[ LONGITUDE ][ MINIMUM ];
    const double east  = arguments->bounds[ LONGITUDE ][ MAXIMUM ];
    const double south = arguments->bounds[ LATITUDE  ][ MINIMUM ];
    const double north = arguments->bounds[ LATITUDE  ][ MAXIMUM ];
    double* const columnValues  = data->columnValues;
    char* output = data->outputBuffer;
    const char* const outputEnd = output + data->bufferSize / sizeof (char);
    char* line = 0;
    char* end = 0;
    int previous_id = 0;
    long long previous_yyyymmddhh = 0;
    const char* const sensorType =
      arguments->sensorType ? arguments->sensorType : "";
    const char* const sensorTypeBlank = arguments->sensorType ? " " : "";
    const char* const sensorTypeDelimiter = arguments->sensorType ? ";" : "";

    eraseQuotedCommasAndQuotesAndBrackets( data->inputBuffer ); /*perf hotspot*/
    memset( data->outputBuffer, 0, data->bufferSize );

    /* Get each line */

    for ( line = strtok_r( data->inputBuffer, "\n", &end );
          line;
          line = strtok_r( 0, "\n", &end ) ) {

      int inSubset =
        parseColumnValues( line,
                           yyyymmddhh1, yyyymmddhh2, west, east, south, north,
                           variable, data->units,
                           data->validMinimum, data->validMaximum,
                           columnValues, data->note );

      /* Append subset data to string buffer for speed: */

      if ( inSubset ) {
        const int id               = columnValues[ ID_INDEX ];
        const long long yyyymmddhh = columnValues[ TIMESTAMP_INDEX ];

        if ( id == previous_id && yyyymmddhh == previous_yyyymmddhh ) {
          inSubset = 0; /* Skip extra measures within the hour! */
        } else {
          previous_id = id;
          previous_yyyymmddhh = yyyymmddhh;
        }

        if ( inSubset ) {
          const double longitude     = columnValues[ LONGITUDE_INDEX ];
          const double latitude      = columnValues[ LATITUDE_INDEX ];
          const double measure       = columnValues[ VARIABLE_INDEX ];
          const char* const description = data->note;
          enum { OUTPUT_LINE_LENGTH = 255 };
          char outputLine[ OUTPUT_LINE_LENGTH + 1 ] = "";
          const size_t length =
            snprintf( outputLine, sizeof outputLine / sizeof *outputLine,
                      "%d,%lld,%f,%f,%f,%s%s%s%s\n",
                      id, yyyymmddhh, longitude, latitude,
                      measure,
                      sensorTypeBlank,
                      sensorType,
                      sensorTypeDelimiter,
                      description );

          /* If line is not truncated and line fits in outputBuffer, append: */

          if ( length < sizeof outputLine / sizeof *outputLine &&
               output + length < outputEnd ) {
            strcpy( output, outputLine );
            output += length;

            if ( data->arguments.aggregate == AGGREGATE_NONE ) {
              data->points += 1; /* Update number of output points. */
            }
          }
        }
      }
    } /* End of parsing data lines. */

    if ( data->outputBuffer[ 0 ] ) { /* Write string buffer to temp file: */
      data->ok =
        fputs( data->outputBuffer, data->tempFiles[ TEMP_FILE_1 ] ) != EOF;
      data->outputBuffer[ 0 ] = '\0';
    }
  }
}



/******************************************************************************
PURPOSE: parseColumnValues - Parse column values from a line.
INPUTS:  char* line                      Data line to parse.
         const long long yyyymmddhh1     Beginning time of subset.
         const long long yyyymmddhh2     Ending    time of subset.
         const double west               West bound of subset.
         const double east               East bound of subset.
         const double south              South bound of subset.
         const double north              North bound of subset.
         const char* const variable      Varaible to read.
         const char* const units         Units of measure.
         const double validMinimum       Minimum valid value of measure.
         const double validMaximum       Maximum valid value of measure.
OUTPUTS: double columnValues[ VARIABLES ]  Parsed column values.
         Note note                       Description of sensor location
RETURNS: int 1 if line values are in subset, else 0.
******************************************************************************/

static int parseColumnValues( char* line,
                              const long long yyyymmddhh1,
                              const long long yyyymmddhh2,
                              const double west,
                              const double east,
                              const double south,
                              const double north,
                              const char* const variable,
                              const char* const units,
                              const double validMinimum,
                              const double validMaximum,
                              double columnValues[],
                              Note note ) {
  int inSubset = 0;
  int column= 0;

  assert( line );
  assert( isValidYYYYMMDDHH( yyyymmddhh1 ) );
  assert( isValidYYYYMMDDHH( yyyymmddhh2 ) );
  assert( yyyymmddhh2 >= yyyymmddhh1 );
  assert( IN_RANGE( west, -180.0, 180.0 ) );
  assert( IN_RANGE( east, west, 180.0 ) );
  assert( IN_RANGE( south, -90.0, 90.0 ) );
  assert( IN_RANGE( north, south, 90.0 ) );
  assert( variable ); assert( *variable );
  assert( units ); assert( *units );
  assert( columnValues );
  assert( note );

  {
    char* value = 0;
    char* endValue = 0;
    double scale = 1.0;

    for ( column = 0; column < VARIABLES; ++column ) {
      columnValues[ column ] = MISSING_VALUE;
    }

    memset( note, 0, sizeof (Note) );
    note[ 0 ] = '-';

    for ( value = strtok_r( line, ",", &endValue ), inSubset = 1, column = 0;
          value && inSubset;
          value = strtok_r( 0, ",", &endValue ), ++column ) {

      if ( *value ) {

        if ( column == 0 ) {
          long long id = 0;
          inSubset = parseLongLong( value, 1, LLONG_MAX, &id );
          columnValues[ ID_INDEX ] = id;
        } else if ( column == 2 ) {
          strncpy( note, value, sizeof (Note) / sizeof (char) - 1 );
          filterNotes( note, sizeof (Note) ); /* Erase problematic chars! */
        } else if ( column == 3 ) {
          long long timestamp = 0;
          inSubset =
            parseTimestamp( value, yyyymmddhh1, yyyymmddhh2, &timestamp );
          columnValues[ TIMESTAMP_INDEX ] = timestamp;
        } else if ( column == 4 ) {
          inSubset =
            parseDouble( value, south, north, columnValues + LATITUDE_INDEX );
        } else if ( column == 5 ) {
          inSubset =
            parseDouble( value, west, east, columnValues + LONGITUDE_INDEX );
        } else if ( column == 6 ) {

          if ( ! strcmp( variable, "humidity" ) ) {
            inSubset = strstr( value, variable ) != 0; /* relativehumidity. */
          } else {

            if ( *value == ' ' ) {
              const char* const value1 = value + 1;
              inSubset = strstr( value1, variable ) == value1;

              if ( inSubset ) {
                inSubset = value1[ strlen( variable ) ] == ' ';
              }
            } else {
              inSubset = ! strcmp( value, variable );
            }
          }
        } else if ( column == 7 ) {
          inSubset = unitsMatch( value, units, &scale );
        } else if ( column == 8 ) {

          if ( scale == 1.0 ) {
            inSubset = parseDouble( value, validMinimum, validMaximum,
                                    columnValues + VARIABLE_INDEX );
          } else {
            double measure = 0.0;
            inSubset = parseDouble( value, -DBL_MAX, DBL_MAX, &measure );

            if ( inSubset ) {
              measure *= scale;
              inSubset = IN_RANGE( measure, validMinimum, validMaximum );

              if ( inSubset ) {
                columnValues[ VARIABLE_INDEX ] = measure;
              }
            }
          }
        }

        DEBUG( fprintf( stderr, "column %d, value = '%s', inSubset = %d\n",
                        column, value, inSubset ); )
      }
    }
  }

  for ( column = 0; column < VARIABLES && inSubset; ++column ) {
    inSubset = columnValues[ column ] != MISSING_VALUE;
  }

  return inSubset;
}



/******************************************************************************
PURPOSE: aggregateData - Compute per-id means of temp subset data extracted
         thus far.
INPUTS:  Data* data  data->tempFiles[ TEMP_FILE_2 ]  id-time sorted temp data
                                                     to aggregate.
                     data->points  Previous number of output points.
OUTPUTS: Data* data  data->tempFiles[ TEMP_AGGREGATED_FILE ]
                     Aggregated data appended.
                     data->points  Increased number of output points.
                     data->ok = 1 if successful, else 0.
******************************************************************************/

static void aggregateData( Data* data ) {
  assert( data );
  assert( data->tempFiles[ TEMP_FILE_2 ] == 0 );
  assert( data->tempFiles[ TEMP_AGGREGATED_FILE ] );
  assert( data->arguments.aggregate != AGGREGATE_NONE );
  data->tempFiles[ TEMP_FILE_2 ] =
    fopen( data->tempFileNames[ TEMP_FILE_2 ], "rb" );
  data->ok = data->tempFiles[ TEMP_FILE_2 ] != 0;

  if ( data->ok ) {
    const int aggregate = data->arguments.aggregate;
    FILE* const input = data->tempFiles[ TEMP_FILE_2 ];
    FILE* const output = data->tempFiles[ TEMP_AGGREGATED_FILE ];
    enum { LINE_LENGTH = 255 };
    char line[ LINE_LENGTH + 1 ] = "";
    int id = 0;
    long long timestamp = 0;
    double longitude = 0.0;
    double latitude  = 0.0;
    double measure   = 0.0;
    double sum = 0.0;
    long long count = 0;
    long long previousTimestamp = 0;
    long long firstTimestamp = 0;

    DEBUG( fprintf( stderr, "start aggregateData() data->points = %lu\n",
                    data->points ); )

    memset( line, 0, sizeof line );
    memset( data->note, 0, sizeof data->note );

    /*
     * Read each line and parse column values comparing id.
     * Optimized for most common case of consecutive id and
     * use atof() which is faster that sscanf():
     */

    while ( data->ok &&
            fgets( line, sizeof line, input ) != 0 ) {
      int reinitialize = 0;
      const int thisId = atoi( line ); /* atoi() is faster than sscanf(). */

      if ( thisId == id ) { /* Handle most frequent case first: aggregate sum*/

        /* Parse measure which appears after penultimate comma: */

        char* comma = strrchr( line, ',' );
        int commaCount = comma != 0;

        while ( comma > line && commaCount < 2 ) {
          --comma;

          if ( *comma == ',' ) {
            ++commaCount;
          }
        }

        data->ok = commaCount == 2 && comma > line;

        if ( data->ok ) {
          measure = atof( comma + 1 );

          if ( aggregate == AGGREGATE_ALL ) {
            sum += measure; /* Aggregate. */
            ++count;
          } else { /* Must read timestamp which is after the first comma: */
            comma = strchr( line, ',' );

            if ( comma ) {
              timestamp = atoll( comma + 1 );

              if ( aggregate == AGGREGATE_HOURLY ||
                   aggregate == AGGREGATE_DAILY ) {
                const int matched =
                  aggregate == AGGREGATE_HOURLY ?
                    timestamp == previousTimestamp
                  : timestamp / 100 == previousTimestamp / 100;

                if ( matched ) {
                  sum += measure; /* Aggregate. */
                  ++count;
                } else {
                  reinitialize = 1;
                }
              }
            }
          }
        }
      }

      if ( thisId != id || reinitialize ) {

        if ( id != 0 && count > 0 ) {

          /* Output mean: */

          const double mean = sum / count;
          const long long outputTimestamp =
            aggregate == AGGREGATE_DAILY ? previousTimestamp / 100 * 100
            : aggregate == AGGREGATE_ALL ? firstTimestamp / 100 * 100
            : previousTimestamp;

          data->ok =
            fprintf( output, "%d,%lld,%lf,%lf,%lf,%s\n",
                     id, outputTimestamp, longitude, latitude, mean,
                     data->note ) > 6;
          data->points += data->ok;
        }

        /* Re-initialize: */

        data->ok =
          sscanf( line, "%d,%lld,%lf,%lf,%lf",
                  &id, &timestamp, &longitude, &latitude, &measure ) == 5;
        sum = measure;
        count = 1;
        previousTimestamp = timestamp;

        if ( firstTimestamp == 0 ) {
          firstTimestamp = timestamp;
          DEBUG( fprintf( stderr, "firstTimestamp = %lld\n", firstTimestamp ); )
        }

        /* Parse description (appearing after last comma) into data->note: */

        {
          char* comma = strrchr( line, ',' );
          data->ok = comma != 0;

          if ( data->ok ) {
            char* newline = 0;
            strncpy( data->note, comma + 1,
                     sizeof data->note / sizeof data->note[ 0 ] - 1 );
            newline = strchr( data->note, '\n' );

            if ( newline ) {
              *newline = '\0';
            }
          }
        }
      }
    }

    DEBUG( fprintf( stderr, "  count = %lld\n", count ); )

    /* Write final mean: */

    if ( data->ok && count > 0 ) {
      const double mean = sum / count;
      const long long outputTimestamp =
        aggregate == AGGREGATE_DAILY ? previousTimestamp / 100 * 100
        : aggregate == AGGREGATE_ALL ? firstTimestamp / 100 * 100
        : previousTimestamp;

      data->ok =
        fprintf( output, "%d,%lld,%lf,%lf,%lf,%s\n",
                 id, outputTimestamp, longitude, latitude, mean, data->note )
         > 6;
      data->points += data->ok; /* Update number of output points. */
    }
  }

  DEBUG( fprintf( stderr,
                  "end aggregateData() data->points = %lu, data->ok = %d\n",
                  data->points, data->ok ); )
}



/******************************************************************************
PURPOSE: reformatData - Reformat final subset temp file to output format.
INPUTS:  Data* data  data->tempFiles[ TEMP_FILE_2 ]  id-time sorted temp data
                                                     to reformat or
                     data->tempFiles[ TEMP_AGGREGATED_FILE ]
                     Aggregated data to reformat.
                     data->points  Number of output points.
OUTPUTS: Data* data  data->tempFiles[ TEMP_FILE_1 ] if ascii else
                     data->tempFiles[ VARIABLES ] if xdr.
                     Reformatted data.
                     data->ok = 1 if successful, else 0.
NOTES: If data was aggregated then read from the temp aggregated file
       else read from temp file 2.
       Write intermediate unsorted reformatted file to temp file 1.
       Then sort temp file 1 into temp file 2.
       If format is ascii, the final result is temp file 2,
       else for xdr format the final results will be in per-variable temp files
******************************************************************************/

static void reformatData( Data* data ) {
  FILE* input = 0;
  int inputIndex = 0;
  assert( data );
  assert( data->ok );
  assert( data->points > 0 );
  assert( data->tempFiles[ TEMP_FILE_1 ] == 0 );
  assert( data->tempFiles[ TEMP_FILE_2 ] == 0 );
  assert( data->tempFiles[ TEMP_AGGREGATED_FILE ] == 0 );

  inputIndex =
    data->arguments.aggregate == AGGREGATE_NONE ? TEMP_FILE_2
    : TEMP_AGGREGATED_FILE;

  input =
    data->tempFiles[ inputIndex ] =
    fopen( data->tempFileNames[ inputIndex ], "rb" );

  data->ok = input != 0;

  if ( input ) {
    FILE* output =
      data->tempFiles[ TEMP_FILE_1 ] =
      fopen( data->tempFileNames[ TEMP_FILE_1 ], "wb" );
    data->ok = output != 0;

    if ( output ) {
      const int isASCII = data->arguments.format == FORMAT_ASCII;
      enum { LINE_LENGTH = 255 };
      char line[ LINE_LENGTH + 1 ] = "";
      size_t points = 0;
      int id = -1;
      long long yyyymmddhh = 0;
      double longitude = 0.0;
      double latitude  = 0.0;
      double measure   = 0.0;

      memset( line, 0, sizeof line );
      memset( data->note, 0, sizeof data->note );

      /* Read each line: */

      while ( data->ok &&
              fgets( line, sizeof line, input ) != 0 ) {

        data->ok =
          sscanf( line, "%d,%lld,%lf,%lf,%lf",
                  &id, &yyyymmddhh, &longitude, &latitude, &measure )
          == VARIABLES;

        if ( data->ok ) {

          /* Parse description (after last comma) into data->note: */

          char* comma = strrchr( line, ',' );
          data->ok = comma != 0;

          if ( data->ok ) {
            const long long yyyymmddhhmmss = yyyymmddhh * 10000;
            char* newline = 0;
            strncpy( data->note, comma + 2, /* Skip space replacing quote. */
                     sizeof data->note / sizeof data->note[ 0 ] - 1 );
            newline = strchr( data->note, '\n' );

            if ( newline ) {
              *newline = '\0';
            }

            /* Write reformatted line to temp file 1: */

            if ( isASCII ) {
              const int yyyy = yyyymmddhhmmss / 10000000000LL;
              const int mm   = yyyymmddhhmmss / 100000000 % 100;
              const int dd   = yyyymmddhhmmss / 1000000 % 100;
              const int hh   = yyyymmddhhmmss / 10000 % 100;
              const int mm2  = yyyymmddhhmmss / 100 % 100;
              const int ss   = yyyymmddhhmmss % 100;
              data->ok =
                fprintf( output,
                         "%d-%02d-%02dT%02d:%02d:%02d-0000"
                         "\t%f\t%f\t%d\t%f\t%s\n",
                         yyyy, mm, dd, hh, mm2, ss,
                         longitude, latitude, id, measure,
                         data->note ) > 11;
            } else {
              data->ok =
                fprintf( output,
                         "%lld,%f,%f,%d,%f,%s\n",
                         yyyymmddhhmmss,
                         longitude, latitude, id, measure,
                         data->note ) > 6;
            }

            points += data->ok;
          }
        }
      }

      /* data->ok = points == data->points; */

      data->points = points;
      data->ok = data->points > 0;

      if ( data->ok ) {

        /* Sort temp file 1 into temp file 2 and empty temp file 1: */

        sortTempData( data );

        if ( data->ok && ! isASCII ) { /* Read temp 2. Write per-var temps: */
          writeTempVariableFiles( data );
        }
      }
    }
  }

  closeTempFiles( data );
}



/******************************************************************************
PURPOSE: writeTempVariableFiles - Read temp file 2 and write per-variable
         temp files in big-endian binary format.
INPUTS:  Data* data  data->tempFiles[ TEMP_FILE_2 ]  time sorted temp data
                                                     to read.
                     data->points  Number of output points.
OUTPUTS: Data* data  data->tempFiles[ VARIABLES ]  Temp files to write.
                     data->ok = 1 if successful, else 0.
******************************************************************************/

static void writeTempVariableFiles( Data* data ) {
  enum { LINE_LENGTH = 255 };
  char line[ LINE_LENGTH + 1 ] = "";
  double longitude = 0.0;
  double latitude  = 0.0;
  double measure   = 0.0;
  long long yyyymmddhh = 0;
  int id = 0;
  size_t points = 0;

  assert( data );
  assert( data->ok );
  assert( data->points > 0 );
  assert( data->tempFiles[ TEMP_FILE_1 ] );

  memset( line, 0, sizeof line );
  memset( data->note, 0, sizeof data->note );

  /* Read each line: */

  openVariableTempFiles( data );

  if ( data->ok ) {
    FILE* notesOutput = data->tempFiles[ TEMP_FILE_1 ];
    FILE* input =
      data->tempFiles[ TEMP_FILE_2 ] =
      fopen( data->tempFileNames[ TEMP_FILE_2 ], "rb" );
      data->ok = input != 0;

    if ( input ) {
      FILE** tempFiles = data->tempFiles;

      /* Read each line: */

      while ( data->ok &&
              fgets( line, sizeof line, input ) != 0 ) {
        data->ok =
          sscanf( line, "%lld,%lf,%lf,%d,%lf",
                  &yyyymmddhh,
                  &longitude, &latitude,
                  &id, &measure ) == 5;

        /* Parse description (after last comma) into data->note: */

        if ( data->ok ) {
          char* comma = strrchr( line, ',' );
          data->ok = comma != 0;

          if ( data->ok ) {
            double words[ VARIABLES ] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
            char* newline = 0;
            strncpy( data->note, comma + 1,
                     sizeof data->note / sizeof data->note[ 0 ] - 1 );
            newline = strchr( data->note, '\n' );

            if ( newline ) {
              *newline = '\0';
            }

            /* Write note to temp file 1: */

            data->ok = fprintf( notesOutput, "%-79s\n", data->note ) == 80;

            if ( data->ok ) { /* Write data to per-variable temp files: */
              words[ 0 ] = yyyymmddhh;
              words[ 1 ] = longitude;
              words[ 2 ] = latitude;
              words[ 3 ] = id;
              words[ 4 ] = measure;
              rotate8ByteArrayIfLittleEndian( words, VARIABLES );
              data->ok =
                fwrite( &words[ 0 ], 8, 1, tempFiles[ 0 ] ) == 1 &&
                fwrite( &words[ 1 ], 8, 1, tempFiles[ 1 ] ) == 1 &&
                fwrite( &words[ 2 ], 8, 1, tempFiles[ 2 ] ) == 1 &&
                fwrite( &words[ 3 ], 8, 1, tempFiles[ 3 ] ) == 1 &&
                fwrite( &words[ 4 ], 8, 1, tempFiles[ 4 ] ) == 1;
              points += data->ok;
            }
          }
        }
      }
    }
  }

  data->ok = points == data->points;
}




/******************************************************************************
PURPOSE: sortTempData - Sort temp file 1 into temp file 2 then empty temp 1.
INPUTS:  Data* data  data->tempFileNames[ TEMP_FILE_1 ]  Name of file to sort.
                     data->tempFileNames[ TEMP_FILE_2 ]  Name of sorted file.
OUTPUTS: Data* data  data->tempFileNames[ TEMP_FILE_1 ]  Empty file.
                     data->tempFileNames[ TEMP_FILE_2 ]  Sorted closed file.
                     data->tempFiles[ TEMP_FILE_1 ]  Rewound.
                     data->tempFiles[ TEMP_FILE_2 ]  Closed.
                     data->ok = 1 if successful, else 0.
NOTES:   HACK: invokes UNIX /usr/bin/sort!
******************************************************************************/

static void sortTempData( Data* data ) {
  char command[ 1024 ] = "";
  assert( data );
  assert( data->tempFileNames ); assert( data->tempFileNames[ 0 ] );

  /* Close temp files to flush content: */

  if ( data->tempFiles[ TEMP_FILE_1 ] ) {
    fclose( data->tempFiles[ TEMP_FILE_1 ] );
    data->tempFiles[ TEMP_FILE_1 ] = 0;
  }

  if ( data->tempFiles[ TEMP_FILE_2 ] ) {
    fclose( data->tempFiles[ TEMP_FILE_2 ] );
    data->tempFiles[ TEMP_FILE_2 ] = 0;
  }

  /* Sort temp file 1 into temp file 2 (overwriting it): */

  memset( command, 0, sizeof command );
  snprintf( command, sizeof command / sizeof *command - 1,
           "/usr/bin/sort -T %s -o %s %s",
            data->arguments.tmpdir,
            data->tempFileNames[ TEMP_FILE_2 ],
            data->tempFileNames[ TEMP_FILE_1 ] );
  data->ok = system( command ) == 0;

  if ( data->ok ) { /* Remove temp file 1 then reopen it for more appending: */
    data->ok = unlink( data->tempFileNames[ TEMP_FILE_1 ] ) == 0;
    data->tempFiles[ TEMP_FILE_1 ] =
      fopen( data->tempFileNames[ TEMP_FILE_1 ], "wb" );
    data->ok = data->tempFiles[ TEMP_FILE_1 ] != 0;
  } else {
    fprintf( stderr, "\nFailed to sort file %s into %s\n",
             data->tempFileNames[ TEMP_FILE_1 ],
             data->tempFileNames[ TEMP_FILE_2 ] );
  }
}



/******************************************************************************
PURPOSE: streamData - Write final content of temp files to stdout.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* const data  data->ok, tempFiles[] = 0 (closed and removed).
******************************************************************************/

static void streamData( Data* const data ) {

  assert( data ); assert( data->ok ); assert( data->points > 0 );
  assert( data->outputBuffer ); assert( data->bufferSize );

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

  {
    const size_t bytes = data->bufferSize;
    void* buffer = data->outputBuffer;
    const int isASCII = data->arguments.format == FORMAT_ASCII;
    const int first = isASCII ? TEMP_FILE_2 : -1; /* -1 = notes first. */
    const int last  = isASCII ? first : VARIABLES - 1;
    int fileIndex = 0;
    int index = 0;
    streamHeader( data );

    for ( index = first; data->ok && index <= last; ++index ) {
      fileIndex = index < 0 ? TEMP_FILE_1 : index;

      data->tempFiles[ fileIndex ] =
        fopen( data->tempFileNames[ fileIndex ], "rb" );
      data->ok = data->tempFiles[ fileIndex ] != 0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nCan't open temp data file '%s' for reading.\n",
                 data->tempFileNames[ fileIndex ] );
      } else {

        while ( data->ok && ! feof( data->tempFiles[ fileIndex ] ) ) {
          const size_t bytesRead =
            fread( buffer, 1, bytes, data->tempFiles[ fileIndex ] );

          if ( bytesRead ) {
            const size_t bytesWritten = fwrite(buffer, 1, bytesRead, stdout);
            data->ok = bytesWritten == bytesRead;
          }
        }
      }
    }

    if ( ! data->ok ) {
      fprintf( stderr, "\nFailed to stream subset data from temp file '%s'.\n",
               data->tempFileNames[ fileIndex ] );
    }
  }
}



/******************************************************************************
PURPOSE: streamHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
******************************************************************************/

static void streamHeader( const Data* const data ) {
  const Arguments* const arguments = &data->arguments;

  /* Append to name to indicate aggregation: */

  const char* const aggregation =
    arguments->aggregate == AGGREGATE_DAILY ? "_daily"
    : arguments->aggregate == AGGREGATE_ALL ? "_mean"
    : "";

  if ( arguments->format == FORMAT_ASCII ) {
    printf( "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\t"
            "STATION(-)\t%s%s(%s)\tNOTE\n",
            arguments->variable, aggregation, data->units );
  } else {
    const long long yyyymmddhhmmss1 = arguments->yyyymmddhh[ 0 ] * 10000;
    const long long yyyymmddhhmmss2 = arguments->yyyymmddhh[ 1 ] * 10000;
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
    printf( "Point 1.0\n"
            "%s\n"
            "%04d-%02d-%02dT%02d:%02d:%02d-0000 "
            "%04d-%02d-%02dT%02d:%02d:%02d-0000\n",
            arguments->description,
            yyyy1, mm1, dd1, hh1, mm21, ss1,
            yyyy2, mm2, dd2, hh2, mm22, ss2 );
    printf( "# Dimensions: variables points:\n%d %lu\n",
            VARIABLES, data->points );
    printf( "# Variable names:\n" );
    printf( "timestamp longitude latitude id %s%s\n",
            arguments->variable, aggregation );
    printf( "# Variable units:\nyyyymmddhhmmss deg deg - %s\n",
            data->units );
    printf( "# char notes[points][80] and\n" );
    printf( "# IEEE-754 64-bit reals data[variables][points]:\n" );
  }
}




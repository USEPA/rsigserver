/******************************************************************************
PURPOSE: GOESBBSubset.c - Extract a lon-lat subset of data from a list of
         GOESBB NetCDF files and write it to stdout as XDR binary format.

NOTES:   Uses NetCDF 3 library.
         Compile:
         gcc -Wall -g -o GOESBBSubset GOESBBSubset.c Utilities.c ReadData.c \
                   -I../../../include/NetCDF4 \
                   -L../../../lib/$platform \
                   -lNetCDF -lm -lc

         Usage:
         GOESBBSubset -files <listfile> \
                      -tmpdir <temp_directory> \
                      -desc "description text" \
                      -timestamp <yyyymmddhh> -hours <count> \
                      -variable <name> \
                      -domain <minimum_longitude> <minimum_latitude> \
                              <maximum_longitude> <maximum_latitude>

          Example:
          ../../../bin/$platform/GOESBBSubset \
          -files testdata/file_list \
          -tmpdir testdata \
  -desc "https://satepsanone.nesdis.noaa.gov/pub/FIRE/BBEP-geo,GOESBBSubset" \
          -timestamp 2018110900 -hours 24 \
          -variable PM2.5_emission \
          -domain -130 20 -60 50 \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 10-line ASCII header followed by binary 64-bit big-endian arrays:

Point 1.0
https://satepsanone.nesdis.noaa.gov/pub/FIRE/BBEP-geo,GOESBBSubset
2018-11-09T00:00:00-0000 2018-11-09T23:59:59-0000
# Dimensions: variables points
4 24
# Variable names:
Timestamp Longitude Latitude PM25_emission
# Variable units:
yyyymmddhhmmss deg deg kg
# IEEE-754 64-bit reals data[variables][points]:
<big-endian binary format array>

HISTORY: 2018-11-11 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/


/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, printf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <stdlib.h>    /* For malloc(), free(). */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For LONGITUDE, Bounds, expand32BitReals(). */
#include "ReadData.h"  /* For openFile(), readFileData(). */

/*================================= MACROS =================================*/

/* Name of temporary file created in -tmpdir will have PID appended: */

#define TEMP_FILE_NAME "junk_GOESBBSubset"

/*================================== TYPES ==================================*/

enum { VARIABLES = 4 }; /* Timestamp, Longitude, Latitude, PM25_emission. */

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;       /* File containing list of GOESBB files to read*/
  const char* tmpdir;         /* Name of directory to write temp files. */
  const char* description;    /* User-supplied description. */
  const char* variable;       /* Name of variable to read. */
  Bounds      domain; /* Subset domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  int         yyyymmddhh;     /* First timestamp of subset. */
  int         hours;          /* Number of hours in subset. */
} Arguments;

/* Data type: */

typedef struct {
  Arguments  arguments;     /* User-supplied (command-line) arguments. */
  char units[ 80 ];         /* Units of variable. */
  FileName   tempFileNames[ VARIABLES ];
  /* Pathed names of temp files of output subset data. */
  FILE*      tempFiles[ 4 ]; /* Temp files of output subset data. */
  int        points;        /* Number of valid data points in subset. */
  int        ok;            /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static void removeTempFiles( Data* const data );

static void createTempFiles( Data* const data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static void readData( Data* const data );

static int dataFileTimestamp( const char* const fileName );

static int readFileInfo( const char* const fileName,
                         int* const file, int* const yyyymmddhh,
                         size_t* const timesteps, size_t* const points,
                         int* const changedDimensions );

static void readCoordinatesAndValues( Data* const data,
                                      const int file,
                                      const int yyyymmddhh,
                                      const size_t timesteps,
                                      const size_t points,
                                      double* const timestamps,
                                      double* const longitudes,
                                      double* const latitudes,
                                      double* const values );

static void replicateData( const size_t timesteps, const size_t points,
                           double data[] );

static void reorderData( const size_t timesteps, const size_t points,
                         const double input[], double output[] );

static void computeTimestamps( const int yyyymmddhh,
                               const size_t timesteps, const size_t points,
                               double timestamps[] );

static void writeSubset( Data* const data,
                         const size_t timesteps,
                         const size_t points,
                         double timestamps[],
                         double longitudes[],
                         double latitudes[],
                         double values[] );

static size_t subsetData( Data* const data,
                          const size_t timesteps,
                          const size_t points,
                          double timestamps[],
                          double longitudes[],
                          double latitudes[],
                          double values[] );

static void streamData( Data* const data );

static void streamHeader( const Data* const data );


/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of GOESBB files and write
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
    readData( &data ); /* Read subset of GOESBB files and write temp files. */

    if ( data.ok && data.points > 0 ) {
      streamData( &data ); /* Write header and temp files to stdout & rm temps*/
    }
  }

  removeTempFiles( &data );
  data.ok = data.ok && data.points > 0;
  return ! data.ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: removeTempFiles - Close and remove temp files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void removeTempFiles( Data* const data ) {
  int variable = 0;
  assert( data );


  for ( variable = 0; variable < VARIABLES; ++variable ) {

    if ( data->tempFiles[ variable ] ) {
      fclose( data->tempFiles[ variable ] );
      data->tempFiles[ variable ] = 0;
    }

    unlink( data->tempFileNames[ variable ] );
    memset( data->tempFileNames[ variable ],
            0, sizeof data->tempFileNames[ variable ] );
  }
}



/******************************************************************************
PURPOSE: createTempFiles - Create temp output files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->ok, data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void createTempFiles( Data* const data ) {
  assert( data );
   /* Temp file names are null. */
  assert( data->tempFiles[ 0 ] == 0 );
  assert( data->tempFiles[ 1 ] == 0 );
  assert( data->tempFiles[ 2 ] == 0 );
  assert( data->tempFiles[ 3 ] == 0 );

  {
    const int pid = getpid();
    const char* const variableNames[ VARIABLES ] = {
      "Timestamp", "Longitude", "Latitude", "Data"
    };
    int variable = 0;

    for ( variable = 0; data->ok && variable < VARIABLES; ++variable ) {
      memset( data->tempFileNames[ variable ], 0, sizeof (FileName) );
      snprintf( data->tempFileNames[ variable ],
                sizeof (FileName) / sizeof (char) - 1,
                "%s/%s_%s.%04d",
                data->arguments.tmpdir, TEMP_FILE_NAME,
                variableNames[ variable ], pid );
      data->tempFiles[ variable ] =
        fopen( data->tempFileNames[ variable ], "wb" );

      if ( ! data->tempFiles[ variable ] ) {
        fprintf( stderr, "\nCan't create temporary output file '%s'.\n",
                 data->tempFileNames[ variable ] );
        data->ok = 0;
      }
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
           "\a\n\n%s - Extract a lon-lat subset of data from a list of\n"
          "GOESBB NetCDF files and write it to stdout as XDR binary format.\n",
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
  fprintf( stderr, " <maximum_longitude> <maximum_latitude>\n\n" );
  fprintf( stderr, "Note:\ntimestamp is in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory were a transient file is " );
  fprintf ( stderr, "written.\nIt should have enough disk space (1TB).\n" );
  fprintf( stderr, "Example:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-files file_list \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc "
                   "\"https://satepsanone.nesdis.noaa.gov/pub/FIRE/BBEP-geo,"
                   ",GOESBBSubset\" \\\n");
  fprintf( stderr, "-timestamp 2018110900 -hours 24 \\\n" );
  fprintf( stderr, "-variable PM25_emission \\\n" );
  fprintf( stderr, "-domain -126 25 -65 50 > subset.xdr\n\n" );
  fprintf( stderr, "PM25 over US on Novembr 9, 2018.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Point 1.0\n" );
  fprintf( stderr, "https://satepsanone.nesdis.noaa.gov/pub/FIRE/BBEP-geo,GOESBBSubset\n" );
  fprintf( stderr, "2018-11-09T00:00:00-0000 2018-11-09T23:59:59-0000\n" );
  fprintf( stderr, "# Dimensions: variables points\n" );
  fprintf( stderr, "4 24\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Timestamp Longitude Latitude PM25_emission\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss deg deg kg\n" );
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
  int result = 0;
  int arg = 0;
  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] ); assert( arguments );
  memset( arguments, 0, sizeof (Arguments) );
  arguments->domain[ LONGITUDE ][ MINIMUM ] = -180.0;
  arguments->domain[ LONGITUDE ][ MAXIMUM ] =  180.0;
  arguments->domain[ LATITUDE  ][ MINIMUM ] =  -90.0;
  arguments->domain[ LATITUDE  ][ MAXIMUM ] =   90.0;

  result = argc == 13 || argc == 18;

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
PURPOSE: readData - Read data from each listed data file and
         write the subset of data to the temporary file.
INPUTS:  Data* data  Data description to read.
OUTPUTS: Data* data  data->points, ok, tempFiles[] = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  size_t length = 0;
  char* listFileContent = readFile( arguments->listFile, &length );
  int wroteSomeData = 0;
  data->ok = listFileContent != 0;

  if ( data->ok ) {
    char* fileName       = 0;
    char* end            = 0;
    size_t timesteps     = 0;
    size_t points        = 0;
    double* buffer       = 0;
    double* timestamps   = 0;
    double* longitudes   = 0;
    double* latitudes    = 0;
    double* values       = 0;

    /* Get each line of list file. It is the GOESBB data file to read: */

    for ( fileName = strtok_r( listFileContent, "\n", &end );
          fileName;
          fileName = strtok_r( 0, "\n", &end ) ) {
      int file = 0;
      int yyyymmddhh = 0;
      int changedDimensions = 0;
      int ok = readFileInfo( fileName,
                             &file, &yyyymmddhh, &timesteps, &points,
                             &changedDimensions );

      if ( ok ) {

        if ( changedDimensions ) {
          const size_t variables = 4; /* timestamp, longitude, latitude, value*/
          const size_t size =  points * timesteps;
          const size_t bytes = variables * size * sizeof (double);

          if ( buffer ) {
            free( buffer );
            buffer = 0;
          }

          buffer = (double*) malloc( bytes );
          data->ok = buffer != 0;

          if ( buffer ) {
            timestamps = buffer;
            longitudes = timestamps + size;
            latitudes  = longitudes + size;
            values     = latitudes  + size;
          } else {
            fprintf( stderr,
                   "\nCan't allocate %lu bytes "
                   "to complete the requested action.\n", bytes );
          }
        }

        if ( data->ok ) {
          readCoordinatesAndValues( data, file, yyyymmddhh,
                                    timesteps, points,
                                    timestamps, longitudes, latitudes, values );
        }

        closeFile( file );
        file = -1;

        if ( data->ok ) {
          writeSubset( data, timesteps, points,
                       timestamps, longitudes, latitudes, values );

          if ( data->ok ) {
            wroteSomeData = 1;
          }
        }
      }
    } /* End loop on listFile. */

    free( listFileContent );
    listFileContent = 0;

    { /* Done writing to temp files so close them: */
      int variable = 0;

      for ( variable = 0; variable < VARIABLES; ++variable ) {

        if ( data->tempFiles[ variable ] ) {
          fclose( data->tempFiles[ variable ] );
          data->tempFiles[ variable ] = 0;
        }
      }
    }
  }

  data->ok = wroteSomeData;
}





/******************************************************************************
PURPOSE: dataFileTimestamp - Timestamp of data file.
INPUTS:  const char* const fileName   Name of GOESBB file.
RETURNS: int yyyymmddhh of file or 0 if failed and message on stderr.
NOTES: Data file names are of the form:
         biomass_burning_YYYYMMDD_HH_HH.nc
       For example,
         biomass_burning_20181109_00_05.nc
       which yields 2018110900.
******************************************************************************/

static int dataFileTimestamp( const char* const fileName ) {
  int result = 0;
  const char* const tag = "biomass_burning_";
  const char* s = strstr( fileName, tag );

  if ( s ) {
    int c = 0;
    s += strlen( tag );

    /* Parse YYYYMMDD: */

    for ( c = 0; *s && *s != '_' && c < 8; ++c, ++s ) {
      result = result * 10 + *s - '0';
    }

    /* Parse HH: */

    if ( *s == '_' ) {
      ++s;

      for ( c = 0; *s && *s != '_' && c < 2; ++c, ++s ) {
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
PURPOSE: readFileInfo - Parse file timestamp, open it and read bounds and
         dimensions.
INPUTS:  const char* const fileName    Name of file to check bounds/dims.
         size_t* const timesteps       Timesteps of data in previous file or 0.
         size_t* const points          Points of data in previous file or 0.
OUTPUTS: int* const file               NetCDF file id of file.
         long long* const yyyymmddhh   Timestamp of file.
         size_t* const timesteps       Timesteps of data in file.
         size_t* const points          Points of data in file.
         int* const changedDimensions  1 if timesteps or points changed.
RETURNS: int 1 if in domain else 0.
******************************************************************************/

static int readFileInfo( const char* const fileName,
                         int* const file, int* const yyyymmddhh,
                         size_t* const timesteps, size_t* const points,
                         int* const changedDimensions ) {

  int result = 0;
  assert( fileName ); assert( *fileName ); assert( file );
  assert( yyyymmddhh );
  assert( timesteps ); assert( points ); assert( changedDimensions );

  *yyyymmddhh = dataFileTimestamp( fileName );
  result = *yyyymmddhh != 0;
  *changedDimensions = 0;

  if ( result ) {
    *file = openFile( fileName );
    result = *file != -1;

    if ( result ) {
      size_t timesteps2 = 0;
      size_t points2 = 0;
      result = readFileDimensions( *file, &timesteps2, &points2 );

      if ( result ) {

        if ( timesteps2 != *timesteps || points2 != *points ) {
          *timesteps = timesteps2;
          *points = points2;
          *changedDimensions = 1;
        }
      }
    }
  }

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
         const int yyyymmddhh          Timestamp of data file.
         const size_t timesteps        Timesteps of data to read.
         const size_t points           Points of data to read.
OUTPUTS: Data* const data              data->ok, data->units[ 80 ]
         double* const timestamps[ timesteps * points ]  Timestamps of data
                                                         points.
         double* const longitudes[ timesteps * points ]  Longitudes read.
         double* const latitudes[  timesteps * points ]  Latitudes read.
         double* const values[     timesteps * points ]  Values read.
******************************************************************************/

static void readCoordinatesAndValues( Data* const data,
                                      const int file,
                                      const int yyyymmddhh,
                                      const size_t timesteps,
                                      const size_t points,
                                      double* const timestamps,
                                      double* const longitudes,
                                      double* const latitudes,
                                      double* const values ) {

  char unused[ 80 ] = "";

  assert( data ); assert( data->ok ); assert( data->arguments.variable );
  assert( file > -1 ); assert( isValidYYYYMMDDHH( yyyymmddhh ) );
  assert( timesteps ); assert( points ); assert( timestamps );
  assert( longitudes ); assert( latitudes ); assert( values );

  data->ok =
    readFileData( file, "Longitude", 0, points, unused, longitudes );

  if ( data->ok ) {
    replicateData( timesteps, points, longitudes );
    data->ok =
      readFileData( file, "Latitude", 0, points, unused, latitudes ) ;

    if ( data->ok ) {
      replicateData( timesteps, points, latitudes );

      if ( data->ok ) {
        const char* const variable = data->arguments.variable;
        size_t timesteps2 = 0;
        size_t points2 = 0;
        data->ok = readVariableDimensions( file, variable,
                                           &timesteps2, &points2 );

        if ( data->ok ) {

          data->ok =
            ( timesteps2 == 0 || timesteps2 == timesteps ) && points2 == points;

          if ( ! data->ok ) {
            fprintf( stderr,
                     "\nUnmatched variable dimensions: "
                     "%s [%lu %lu] expected [%lu %lu].\n",
                     variable, timesteps2, points2, timesteps, points );
          } else if ( timesteps2 == 0 ) { /* 1D variable: */
            data->ok =
              readFileData( file, variable, 0, points, data->units, values );

            if ( data->ok ) {
              replicateData( timesteps, points, values );
            }
          } else { /* 2D variable: */
            data->ok = readFileData( file, variable,
                                     timesteps, points,
                                     data->units, timestamps );

            if ( data->ok ) {
              reorderData( timesteps, points, timestamps, values );
            }
          }

          if ( data->ok ) {
            computeTimestamps( yyyymmddhh, timesteps, points, timestamps );
          }
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: replicateData - Replicate data values for each timestep.
INPUTS:  const size_t timesteps             Timesteps of data.
         const size_t points                Points of data.
         const double data[ points ]        Data to replicate.
OUTPUTS: double data[ timesteps * points ]  Replicated data.
******************************************************************************/

static void replicateData( const size_t timesteps, const size_t points,
                           double data[] ) {

  assert( timesteps ); assert( points ); assert( data );

  {
    double* input  = data + points;
    double* output = data + timesteps * points;
    size_t point = points;

    while ( point-- ) {
      const double value = *--input;
      size_t timestep = timesteps;

      while ( timestep-- ) {
        *--output = value;
      }
    }
  }
}



/******************************************************************************
PURPOSE: reorderData - Reorder 2D array.
INPUTS:  const size_t timesteps              Timesteps of data.
         const size_t points                 Points of data.
         const double input[ points * timesteps ]  Data to reorder.
OUTPUTS: double output[ timesteps * points ]       Reordered data.
******************************************************************************/

static void reorderData( const size_t timesteps, const size_t points,
                         const double input[], double output[] ) {

  size_t timestep = 0;
  double* outputPointer = output;

  assert( timesteps ); assert( points );
  assert( input ); assert( output ); assert( input != output );

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    size_t point = 0;

    for ( point = 0; point < points; ++point ) {
      const size_t index = point * timesteps + timestep;
      const double value = input[ index ];
      *outputPointer++ = value;
    }
  }
}



/******************************************************************************
PURPOSE: computeTimestamps - Compute timestamp yyyymmddhhmmss from starting
         timestamp yyyymmddhh and number of timesteps and points to replicate.
INPUTS:  const int yyyymmddhh                    First timestamp.
         const size_t timesteps                  Timesteps (hours) of data.
         const size_t points                     Points of data.
OUTPUTS: double timestamps[ timesteps * points ] yyyymmddhhmmss per point.
******************************************************************************/

static void computeTimestamps( const int yyyymmddhh,
                               const size_t timesteps, const size_t points,
                               double timestamps[] ) {

  size_t timestep = 0;
  double* outputPointer = timestamps;
  int yyyymmddhh2 = yyyymmddhh;

  assert( isValidYYYYMMDDHH( yyyymmddhh ) );
  assert( timesteps ); assert( points ); assert( timestamps );

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    const double yyyymmddhhmmss = incrementHours( yyyymmddhh2, 1 ) * 10000.0;
    size_t point = 0;

    for ( point = 0; point < points; ++point ) {
      *outputPointer++ = yyyymmddhhmmss;
    }
  }
}



/******************************************************************************
PURPOSE: writeSubset - Write subset of data to temp file.
INPUTS:  Data* const data                    Data description.
         const size_t timesteps              Timesteps of input data.
         const size_t points                 Points of input data.
         double timestamps[ timesteps * points ]  Timestamps to write.
         double longitudes[ timesteps * points ]  Longitudes to write.
         double latitudes[  timesteps * points ]  Latitudes to write.
         double values[     timesteps * points ]  Values to write.
OUTPUTS: Data* data  data->ok, data->tempFileNames, data->tempFiles.
         double timestamps[ timesteps * points ]  Ruined data.
         double longitudes[ timesteps * points ]  Ruined data.
         double latitudes[  timesteps * points ]  Ruined data.
         double values[     timesteps * points ]  Ruined data.
******************************************************************************/

static void writeSubset( Data* const data,
                         const size_t timesteps,
                         const size_t points,
                         double timestamps[],
                         double longitudes[],
                         double latitudes[],
                         double values[] ) {

  size_t subsetCount = 0;

  assert( data );
  assert( timesteps != 0 ); assert( points != 0 ); assert( timestamps );
  assert( longitudes ); assert( latitudes ); assert( values );

  /* Store compact sequence of valid data that is within time/domain subset: */

  subsetCount = subsetData( data, timesteps, points, timestamps,
                            longitudes, latitudes, values );

  if ( subsetCount ) { /* Write subset of data: */
    data->points += subsetCount;

    /* Open temp files for writing if they do not yet exist: */

    if ( data->tempFiles[ 0 ] == 0 ) {
      createTempFiles( data );
    }

    if ( data->ok ) { /* Write subset data to temp file: */
      int variable = 0;
      rotate8ByteArrayIfLittleEndian( timestamps, subsetCount );
      data->ok =
        fwrite( timestamps, sizeof *timestamps, subsetCount,
                data->tempFiles[ variable ] ) == subsetCount;

      if ( data->ok ) {
        ++variable;
        rotate8ByteArrayIfLittleEndian( longitudes, subsetCount );
        data->ok =
          fwrite( longitudes, sizeof *longitudes, subsetCount,
                  data->tempFiles[variable ] ) == subsetCount;

        if ( data->ok ) {
          ++variable;
          rotate8ByteArrayIfLittleEndian( latitudes, subsetCount );
          data->ok =
            fwrite( latitudes, sizeof *latitudes, subsetCount,
                    data->tempFiles[ variable ] ) == subsetCount;

          if ( data->ok ) {
            ++variable;
            rotate8ByteArrayIfLittleEndian( values, subsetCount );
            data->ok =
              fwrite( values, sizeof *values, subsetCount,
                      data->tempFiles[ variable ] ) == subsetCount;
          }
        }
      }

      if ( ! data->ok ) {
        fprintf( stderr, "\nFailed to write subset data to temp files '%s'.\n",
                 data->tempFileNames[ variable ] );
      }
    }
  }
}



/******************************************************************************
PURPOSE: subsetData - Subset/filter and compact data matching subset.
INPUTS:  Data* const data                    Data description.
         const size_t timesteps              Timesteps of input data.
         const size_t points                 Points of input data.
         double timestamps[ timesteps * points ]  Timestamps to filter.
         double longitudes[ timesteps * points ]  Longitudes to filter.
         double latitudes[  timesteps * points ]  Latitudes to filter.
         double values[     timesteps * points ]  Values to filter.
OUTPUTS: Data* data  data->ok, data->tempFileNames, data->tempFiles.
         double timestamps[ result ]  Compact subset data.
         double longitudes[ result ]  Compact subset data.
         double latitudes[  result ]  Compact subset data.
         double values[     result ]  Compact subset data.
RETURNS: size_t number of data values in subset.
******************************************************************************/

static size_t subsetData( Data* const data,
                          const size_t timesteps,
                          const size_t points,
                          double timestamps[],
                          double longitudes[],
                          double latitudes[],
                          double values[] ) {

  size_t result = 0;

  assert( data );
  assert( timesteps != 0 ); assert( points != 0 ); assert( timestamps );
  assert( longitudes ); assert( latitudes ); assert( values );

  /* Store compact sequence of valid data that is within time/domain subset: */

  {
    const double minimumValidValue = 0.0;
    const Arguments* const arguments = &data->arguments;
    const double longitudeMinimum = arguments->domain[ LONGITUDE ][ MINIMUM ];
    const double longitudeMaximum = arguments->domain[ LONGITUDE ][ MAXIMUM ];
    const double latitudeMinimum  = arguments->domain[ LATITUDE  ][ MINIMUM ];
    const double latitudeMaximum  = arguments->domain[ LATITUDE  ][ MAXIMUM ];
    const int yyyymmddhhFirst = arguments->yyyymmddhh;
    const int yyyymmddhhLast =
      incrementHours( yyyymmddhhFirst, arguments->hours - 1 );
    const size_t count = timesteps * points;
    size_t index = 0;

    for ( index = 0; index < count; ++index ) {
      const double yyyymmddhhmmss = timestamps[ index ];
      const int yyyymmddhh = (int) ( yyyymmddhhmmss / 10000.0 );

      if ( yyyymmddhh >= yyyymmddhhFirst && yyyymmddhh <= yyyymmddhhLast ) {
        const double longitude = longitudes[ index ];

        if ( longitude >= longitudeMinimum && longitude <= longitudeMaximum ) {
          const double latitude = latitudes[ index ];

          if ( latitude >= latitudeMinimum && latitude <= latitudeMaximum ) {
            const double value = values[ index ];

            if ( value >= minimumValidValue ) {

              if ( index > result ) {
                timestamps[ result ] = yyyymmddhhmmss;
                longitudes[ result ] = longitude;
                latitudes[  result ] = latitude;
                values[     result ] = value;
                ++result;
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
PURPOSE: streamData - Write content of temp files to stdout.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* const data  data->ok, tempFiles[] = 0 (closed and removed).
******************************************************************************/

static void streamData( Data* const data ) {
  const size_t bytes = 1024 * 1024;
  void* buffer = malloc( bytes );
  int variable = 0;

  assert( data ); assert( data->ok );
   /* Temp file names are initialized but temp files are closed after writing*/
  assert( data->tempFileNames[ 0 ][ 0 ] );
  assert( data->tempFiles[     0 ] == 0 );
  assert( data->tempFileNames[ 1 ][ 0 ] );
  assert( data->tempFiles[     1 ] == 0 );
  assert( data->tempFileNames[ 2 ][ 0 ] );
  assert( data->tempFiles[     2 ] == 0 );
  assert( data->tempFileNames[ 3 ][ 0 ] );
  assert( data->tempFiles[     3 ] == 0 );

  data->ok = buffer != 0;

  if ( ! buffer ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes to complete the requested action.\n",
             bytes );
  } else {
    streamHeader( data );

    for ( variable = 0; data->ok && variable < VARIABLES; ++variable ) {
      data->tempFiles[ variable ] =
        fopen( data->tempFileNames[ variable ], "rb" );
      data->ok = data->tempFiles[ variable ] != 0;

      if ( ! data->ok ) {
        fprintf( stderr, "\nCan't open temp data file '%s' for reading.\n",
                 data->tempFileNames[ variable ] );
      } else {

        while ( data->ok && ! feof( data->tempFiles[ variable ] ) ) {
          const size_t bytesRead =
            fread( buffer, 1, bytes, data->tempFiles[ variable ] );

          if ( bytesRead ) {
            const size_t bytesWritten = fwrite( buffer, 1, bytesRead, stdout );
            data->ok = bytesWritten == bytesRead;
          }
        }
      }
    }

    if ( ! data->ok ) {
      fprintf( stderr, "\nFailed to stream subset data from temp file '%s'.\n",
               data->tempFileNames[ variable ] );
    }

    free( buffer );
    buffer = 0;
  }

  removeTempFiles( data );
}



/******************************************************************************
PURPOSE: streamHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
******************************************************************************/

static void streamHeader( const Data* const data ) {
  const Arguments* const arguments = &data->arguments;
  const int variables = 4;
  const int yyyymmddhh1 = arguments->yyyymmddhh;
  const int yyyymmddhh2 = incrementHours( yyyymmddhh1, arguments->hours - 1 );
  const int yyyy1       = yyyymmddhh1 / 1000000;
  const int mm1         = yyyymmddhh1 / 10000 % 100;
  const int dd1         = yyyymmddhh1 / 100 % 100;
  const int hh1         = yyyymmddhh1 % 100;
  const int yyyy2       = yyyymmddhh2 / 1000000;
  const int mm2         = yyyymmddhh2 / 10000 % 100;
  const int dd2         = yyyymmddhh2 / 100 % 100;
  const int hh2         = yyyymmddhh2 % 100;
  const char* input = arguments->variable;
  enum { LENGTH = 79 };
  char variableName[ LENGTH + 1 ] = "";
  memset( variableName, 0, sizeof variableName );

  if ( strchr( arguments->variable, '.' ) ) {
    int index = 0;

    do {

      if ( *input != '.' ) {
        variableName[ index ] = *input;
        ++index;
      }

      ++input;
    } while ( index < LENGTH && *input != '\0' );
  } else {
    strncpy( variableName, arguments->variable, LENGTH );
  }

  printf( "Point 1.0\n"
          "%s\n"
          "%04d-%02d-%02dT%02d:00:00-0000 %04d-%02d-%02dT%02d:00:00-0000\n",
          arguments->description, yyyy1, mm1, dd1, hh1, yyyy2, mm2, dd2, hh2 );
  printf( "# Dimensions: variables points:\n%d %d\n",
          variables, data->points );
  printf( "# Variable names:\n" );
  printf( "Timestamp Longitude Latitude %s\n", variableName );
  printf( "# Variable units:\nyyyymmddhhmmss deg deg %s\n", data->units );
  printf( "# IEEE-754 64-bit reals data[variables][points]:\n" );
}




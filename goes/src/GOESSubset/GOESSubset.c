
/******************************************************************************
PURPOSE: GOESSubset.c - Read a set of a GOES files, subset the scans to a
         bounds (longitude-latitude rectangle),
         optionally aggregate to a daily mean and
         write it (with optionally computed coordinate corners)
         to stdout as XDR (IEEE-754) format binary.

NOTES:   Uses libz.a (see ../../libs/MODIS/zlib) and
         libUtilities.a (../../libs/Utilities).

         gcc -m64 -fopenmp \
             -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -g -DNO_ASSERTIONS \
             -I../../../include -I. -o ../../../bin/$platform/GOESSubset \
             GOESSubset.c -L/usr/local/lib/x86_64 \
             -L../../../lib/$platform -lUtilities -lz -lm

HISTORY: 2015-05-08 plessel.todd@epa.gov
STATUS: unreviewed untested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For stderr, FILE, fprintf(). */
#include <string.h> /* For strlen(), memset(). */
#include <ctype.h>  /* For isdigit(), isprint(). */

#include <Utilities.h> /* For PRE0*(), NEW_ZERO(), Stream, VoidList, etc. */

/* Z Library routines used (simply prototyped) to read compressed files: */

extern void* gzopen( const char* file_name, const char* mode );
extern int gzread( void* file, void* data, unsigned int bytes );
extern int gzclose( void* file );
extern const char* gzerror( void* file, int* unused );

static int readCompressedFile( const char* fileName,
                               const int bufferSize, void* buffer ) {
  void* inputFile = gzopen( fileName, "rb" );
  int result = 0;
  memset( buffer, 0, bufferSize );

  DEBUG( fprintf( stderr, "Reading file %s\n", fileName ); )

  if ( ! inputFile ) {
    int unused = 0;
    failureMessage( "Failed to open file %s for reading because %s.",
                    fileName, gzerror( inputFile, &unused ) );
  } else {
    const int bytesRead = gzread( inputFile, buffer, bufferSize - 1);
    CHECK( bytesRead <= bufferSize - 1 );
    result = bytesRead > 0;

    if ( ! result ) {
      int unused = 0;
      failureMessage( "Failed to read up to %d bytes from file %s"
                      "because %s.",
                      bufferSize - 1, fileName, gzerror( inputFile, &unused ));
    }

    gzclose( inputFile ), inputFile = 0;
  }

  DEBUG( fprintf( stderr, "result = %d\n", result ); )
  return result;
}

/*================================== TYPES ==================================*/

#define MISSING (-9999.0)

/* For -corners option: */

enum{ VARIABLE_SIZE = 32 };
typedef char Variable[ VARIABLE_SIZE ]; /* E.g., "INSL" ... "PARM". */

/* Input command-line arguments: */

typedef struct {
  const char* listFile;       /* File listing time-sorted GOES files to read */
  const char* description;    /* User-supplied description. */
  const char* units;          /* User-supplied variable units. Not in file! */
  Integer     firstTimestamp; /* YYYYDDDHHMM of subset. */
  Integer     hours;          /* Number of hours in subset. */
  Integer     daily;          /* Compute daily mean of variable? */
  Integer     corners;        /* Compute grid cell corners? */
  Bounds      bounds;         /* bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM] */
} Arguments;

#ifndef NO_ASSERTIONS

/* Arguments invariant: */

static Integer isValidArguments( const Arguments* arguments ) {
  const Integer result =
    AND12( arguments,
           arguments->listFile,
           arguments->listFile[ 0 ],
           arguments->description,
           arguments->description[ 0 ],
           arguments->units,
           arguments->units[ 0 ],
           isValidTimestamp( arguments->firstTimestamp ),
           arguments->hours > 0,
           IS_BOOL( arguments->daily ),
           IS_BOOL( arguments->corners ),
           isValidBounds( arguments->bounds ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Scan: Result of reading a GOES data file. */

typedef struct {
  Integer yyyydddhhmm; /* UTC timestamp of file containing this scan. */
  FileName coordinatesFileName; /* File containing longitudes and latitudes. */
  Integer rows;        /* Number of data rows. */
  Integer columns;     /* Number of data columns. */
  Real*   data;        /* data[ rows ][ columns ]. */
} Scan;

/* Scan destructor: */

static void deallocateScan( Scan* const scan ) {
  PRE0( scan );
  FREE( scan->data );
  ZERO_OBJECT( scan );
  POST0( scan );
}

#ifndef NO_ASSERTIONS

/* Scan invariant: */

static Integer isValidScan( const Scan* scan ) {
  const Integer result =
    AND7( scan,
          isValidTimestamp( scan->yyyydddhhmm ),
          isprint( scan->coordinatesFileName[ 0 ] ),
          scan->rows > 0,
          scan->columns > 0,
          scan->rows * scan->columns > 0,
          isNanFree( scan->data, scan->rows * scan->columns ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* SubsettedScan: After bounds subsetting and data filtering. */

typedef struct {
  Integer timestamp; /* YYYYDDDHHMM of file containing scan. */
  Integer variables; /* 3 or 11: Longitude, Latitude, INSL, Longitude_SE, ...*/
  Integer points;    /* Number of subsetted points. */
  Real*   data;      /* Subsetted data[ variables ][ points ]. */
} SubsettedScan;

/* SubsettedScan destructor (void* argument due to callback from VoidList): */

static void deallocateSubsettedScan( void* argument ) {
  PRE0( argument );
  SubsettedScan* const subsettedScan = argument;
  FREE( subsettedScan->data );
  ZERO_OBJECT( subsettedScan );
  POST0( argument );
}

#ifndef NO_ASSERTIONS

/* SubsettedScan invariant: */

static Integer isValidSubsettedScan( const SubsettedScan* subsettedScan ) {
  const Integer result =
    AND7( subsettedScan,
          isValidTimestamp( subsettedScan->timestamp ),
          IN3( subsettedScan->variables, 3, 11 ),
          subsettedScan->points > 0,
          validLongitudesAndLatitudes( subsettedScan->points,
                                       subsettedScan->data,
                                       subsettedScan->data +
                                         subsettedScan->points ),
          minimumItem( subsettedScan->data + 2 * subsettedScan->points,
                       subsettedScan->points ) >= 0.0,
          IMPLIES( subsettedScan->variables == 11,
                   validLongitudesAndLatitudes( 4 * subsettedScan->points,
                                       subsettedScan->data +
                                         3 * subsettedScan->points,
                                       subsettedScan->data +
                                         7 * subsettedScan->points ) ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Data type: */

typedef struct {
  Integer   ok;                /* Did last command succeed? */
  Integer   rows;              /* Number of coordinate rows. */
  Integer   columns;           /* Number of coordinate columns. */
  Integer   subset[ 2 ][ 2 ];  /* subset[ ROW COLUMN ][ MINIMUM MAXIMUM ]. */
  Variable  variable;          /* Name of data variable. "INSL" ... "PARM". */
  Real*     longitudes;        /* longitudes[ rows ][ columns ]. */
  Real*     latitudes;         /* latitudes[  rows ][ columns ]. */
  Real*     corners;/*corners[LONGITUDE LATITUDE][SW SE NW NE][rows][columns]*/
  Integer   bufferSize;        /* sizeof buffer. */
  char*     buffer;            /* buffer[bufferSize] uncompressed file bytes.*/
  Arguments arguments;         /* User-supplied (command-line) arguments. */
  Scan      scan;              /* Scan of current file being read/processed. */
  VoidList* subsettedScans;    /* List of subsetted scans already read. */
} Data;

/* Data destructor: */

static void deallocateData( Data* data ) {
  PRE0( data );
  FREE( data->longitudes );
  FREE( data->latitudes );
  FREE( data->corners );
  FREE( data->buffer );
  deallocateScan( &data->scan );
  FREE_OBJECT( data->subsettedScans ); /* Calls deallocateSubsettedScan(). */
  ZERO_OBJECT( data );
  POST0( data );
}

#ifndef NO_ASSERTIONS

/* Data invariant: */

static Integer isValidData( const Data* data ) {
  Integer result =
    AND20( data,
           isValidArguments( &data->arguments ),
           isprint( data->variable[ 0 ] ),
           data->rows > 0,
           data->columns > 0,
           data->rows * data->columns > 0,
           IN_RANGE( data->subset[ ROW ][ MINIMUM ], 0, data->rows - 1 ),
           IN_RANGE( data->subset[ ROW ][ MAXIMUM ],
                     data->subset[ ROW ][ MINIMUM ], data->rows - 1 ),
           IN_RANGE( data->subset[ COLUMN][MINIMUM], 0, data->columns - 1),
           IN_RANGE( data->subset[ COLUMN ][ MAXIMUM ],
                     data->subset[ COLUMN ][ MINIMUM ], data->columns - 1),
           data->longitudes,
           data->latitudes,
           validLongitudesAndLatitudes( data->rows * data->columns,
                                        data->longitudes, data->latitudes ),
           IMPLIES_ELSE( data->arguments.corners,
                         AND2( data->corners,
                               validLongitudesAndLatitudes(
                                  4 * data->rows * data->columns,
                                  data->corners,
                                  data->corners +
                                    4 * data->rows * data->columns ) ),
                         ! data->corners ),
           data->bufferSize > 0,
           data->buffer,
           isValidScan( &data->scan ),
           data->subsettedScans,
           data->subsettedScans->invariant( data->subsettedScans ),
           IS_BOOL( data->ok ) );

  if ( result ) {
    const VoidList* const scans = data->subsettedScans;
    const Integer scanCount = scans->count( scans );
    Integer index = 0;

    do {
      const SubsettedScan* const subsettedScan = scans->item( scans, index );
      result = isValidSubsettedScan( subsettedScan );
      ++index;
    } while ( AND2( result, index < scanCount ) );
  }

  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */


/*========================== FORWARD DECLARATIONS ===========================*/

static void copyWord( char* output, const char* input, const size_t length ) {

  if ( input && output && length ) {
    size_t index = 0;

    while ( isspace( *input ) ) {
      ++input;
    }

    while ( *input && ! isspace( *input ) && index < length ) {
      output[ index ] = *input++;
      ++index;
    }

    while ( index < length ) {
      output[ index ] = '\0';
      ++index;
    }
  }
}

static const char* readAndSkipReal( const char* string, Real* value ) {
  const char* result = string;

  if ( result && value ) {
    *value = MISSING;

    while ( isspace( *result ) ) {
      ++result;
    }

    *value = atof( result );

    while ( ! isspace( *result ) ) {
      ++result;
    }

    if ( ! ( *value >= MISSING ) ) {
      *value = 0.0;
      result = 0;
    }
  }

  return result;
}

static void printUsage( void );

static Integer parseArguments( Integer argc, char* argv[],
                               Arguments* arguments );

static void initializeArguments( Arguments* arguments );

static Integer parseOptionalArguments( Integer argc, char* argv[],
                                       Integer* arg,
                                       Arguments* arguments );

static void readData( Data* const data );

static void appendDailyMeans( const Integer yyyyddd0000,
                              const Integer counts[], const Real means[],
                              Data* data );

static void readGOESDataFile( const char* fileName, Data* data );

static Integer findCoordinatesFile( const char* const dataFileName,
                                    FileName coordinatesFileName );

static void readGOESCoordinatesFile( Data* data );

static const char* parseHeader( const char* buffer,
                                Integer* rows, Integer* columns,
                                Variable variable,
                                FileName coordinatesFileName );

static Integer parseCount( const char* header, const char* name );

static const char* parseWord( const char* header, const char* name );

static Integer parseData( const char* buffer, const Integer count,
                          Real data[], Real data2[] );

static Integer timestampOfFileName( const char* fileName );

static SubsettedScan* subsetScan( const Bounds bounds,
                                  const Integer subset[ 2 ][ 2 ],
                                  const Real longitudes[],
                                  const Real latitudes[],
                                  const Real corners[],
                                  Scan* scan );

static Integer subsetScanCount( const Bounds bounds,
                                const Real longitudes[],
                                const Real latitudes[],
                                Integer subset[ 2 ][ 2 ],
                                Scan* scan );

static Integer subsetIndicesByBounds( const Bounds bounds,
                                      const Integer rows,
                                      const Integer columns,
                                      const Real longitudes[],
                                      const Real latitudes[],
                                      Integer* firstRow,
                                      Integer* lastRow,
                                      Integer* firstColumn,
                                      Integer* lastColumn );

static void reduceSubset( const Integer rows, const Integer columns,
                          const Real data[],
                          Integer* firstRow, Integer* lastRow,
                          Integer* firstColumn, Integer* lastColumn );

static void copySubsetLongitudesAndLatitudes( const Integer rows,
                                              const Integer columns,
                                              const Real longitudes[],
                                              const Real latitudes[],
                                              const Real corners[],
                                              const Real data[],
                                              Integer points,
                                              Integer firstRow,
                                              Integer lastRow,
                                              Integer firstColumn,
                                              Integer lastColumn,
                                              Real subsetLongitudes[],
                                              Real subsetLatitudes[],
                                              Real subsetLongitudesSW[],
                                              Real subsetLongitudesSE[],
                                              Real subsetLongitudesNW[],
                                              Real subsetLongitudesNE[],
                                              Real subsetLatitudesSW[],
                                              Real subsetLatitudesSE[],
                                              Real subsetLatitudesNW[],
                                              Real subsetLatitudesNE[] );

static void copySubsetData( const Integer rows, const Integer columns,
                            const Real data[],
                            Integer firstRow, Integer lastRow,
                            Integer firstColumn, Integer lastColumn,
                            Real output[] );

static void computeMean( const Data* const data,
                         Integer counts[], Real means[] );

static Integer meanPoints( const Data* data, const Integer counts[] );

static void copyMeanData( const Integer rows,
                          const Integer columns,
                          const Integer subset[ 2 ][ 2 ],
                          const Integer points,
                          const Integer counts[],
                          const Real means[],
                          const Real longitudes[],
                          const Real latitudes[],
                          const Real corners[],
                          Real subsetLongitudes[],
                          Real subsetLatitudes[],
                          Real subsetData[],
                          Real subsetLongitudesSW[],
                          Real subsetLongitudesSE[],
                          Real subsetLongitudesNW[],
                          Real subsetLongitudesNE[],
                          Real subsetLatitudesSW[],
                          Real subsetLatitudesSE[],
                          Real subsetLatitudesNW[],
                          Real subsetLatitudesNE[] );

static void writeData( Data* data );

static void writeHeader( Data* data, Stream* output );

static void writeXDR( Data* data, Stream* output );

static void writeScanTimestamps( const VoidList* scans, Stream* output );

static void writeScanPoints( const VoidList* scans, Stream* output );

static void writeScanData( const VoidList* scans, Stream* output );

static void computeCorners( const Integer rows, const Integer columns,
                            const Real longitudes[], const Real latitudes[],
                            Real corners[] );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Read a subset of a GOES file and write it to stdout in
         XDR or ASCII format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Integer ok = 0;

  if ( ! isValidArgs( argc, (const char**) argv ) ) {
    failureMessage( "Invalid command-line arguments." );
    printUsage();
  } else {
    Data data;
    memset( &data, 0, sizeof data );
    checkForTest( &argc, argv ); /* Check for and remove any -test arguments.*/
    data.ok = parseArguments( argc, argv, &data.arguments );

    if ( data.ok ) {
      const Integer bufferSize = 100 * 1024 * 1024; /*Decompressed file bytes*/
      data.buffer = NEW_ZERO( char, bufferSize );
      data.ok = data.buffer != 0;

      if ( data.ok ) {
        data.bufferSize = bufferSize;
        readData( &data ); /* From list file named in arguments. */

        if ( data.ok ) {
          writeData( &data ); /* To stdout. */
        }
      }
    }

    ok = data.ok;
    deallocateData( &data );
  }

  ok = AND2( ok, failureCount() == 0 );
  POST0( IS_BOOL( ok ) );
  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
******************************************************************************/

static void printUsage( void ) {
  fprintf( stderr,
           "\n\nGOESSubset - Read a set of GOES files and extract scan\n" );
  fprintf( stderr, "data subsetted by a lon-lat rectangle and" );
  fprintf( stderr, " filtered by variable ranges.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "GOESSubset \\\n" );
  fprintf( stderr, "  -files <listFile> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -units \"units\" \\\n" );
  fprintf( stderr, "  -timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "  [ -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ]\\\n" );
  fprintf( stderr, "  [ -daily ] \\\n" );
  fprintf( stderr, "  [ -corners ]\\\n" );
  fprintf( stderr, "\n\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "  -daily computes daily mean of filtered data\n" );
  fprintf( stderr, "  -corners option will output 8 additional variables:\n" );
  fprintf( stderr, "  Longitude_SW Longitude_SE Longitude_NW Longitude_NE\n");
  fprintf( stderr, "  Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "that are the linearly interpolated " );
  fprintf( stderr, "(and edge extrapolated)\n" );
  fprintf( stderr, "corner points for each center-pixel point.\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example #1:\n\n" );
  fprintf( stderr, "GOESSubset \\\n" );
  fprintf( stderr, "-files testdata/insl_files.txt \\\n" );
  fprintf( stderr,"-desc http://www.nsstc.uah.edu/nsstc/,GOESSubset \\\n");
  fprintf( stderr, "-units W/m2 \\\n" );
  fprintf( stderr, "-timestamp 2013090100 -hours 24 \\\n" );
  fprintf( stderr, "-domain -76 34 -74 36 > subset.xdr\n\n" );
  fprintf( stderr, "Subset of data for September 1, 2013 near Raleigh, NC, USA\n");
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays\n" );
  fprintf( stderr, "For example:\n" );
  fprintf( stderr, "Swath 2.0\n" );
  fprintf( stderr, "http://www.nsstc.uah.edu/nsstc/,GOESSubset\n" );
  fprintf( stderr, "2013-09-01T00:00:00-0000\n" );
  fprintf( stderr, "# Dimensions: variables timesteps scans:\n" );
  fprintf( stderr, "3 24 25\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Longitude Latitude INSL\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "deg deg W/m2\n" );
  fprintf( stderr, "# Domain: <min_lon> <min_lat> <max_lon> <max_lat>\n" );
  fprintf( stderr, "-76 34 -74 36\n" );
  fprintf( stderr, "# MSB 64-bit integers (yyyydddhhmmss)" );
  fprintf( stderr, " timestamps[scans] and\n" );
  fprintf( stderr, "# MSB 64-bit integers" );
  fprintf( stderr, " points[scans] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals" );
  fprintf( stderr, " data_1[variables][points_1] ..." );
  fprintf( stderr, " data_S[variables][points_S]:\n" );
  fprintf( stderr, "<binary data arrays here>\n\n\n" );
  fprintf( stderr, "Example #2:\n\n" );
  fprintf( stderr, "GOESSubset \\\n" );
  fprintf( stderr, "-files /data/tmp/files.txt \\\n" );
  fprintf( stderr, "-desc http://www.nsstc.uah.edu/nsstc/,GOESSubset \\\n");
  fprintf( stderr, "-units W/m2 \\\n" );
  fprintf( stderr, "-timestamp 2013090100 -hours 48 \\\n" );
  fprintf( stderr, "-daily -corners \\\n" );
  fprintf( stderr, "-domain -76 34 -74 36 > subset.xdr\n\n" );
  fprintf( stderr, "Computes daily mean of filtered data with corners.\n" );
  fprintf( stderr, "\n\n\n");
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: Integer 1 if successful, else 0 and failureMessage() then printUsage()
         are called.
******************************************************************************/

static Integer parseArguments( Integer argc, char* argv[],
                               Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  Integer result = 0;
  ZERO_OBJECT( arguments );
  initializeArguments( arguments );

  if ( ! IN_RANGE( argc, 11, 49 ) ) {
    failureMessage( "Invalid/insufficient/redundant command line arguments.");
  } else {
    Integer arg = 1; /* Argument to parse next, updated by parseArgument2(). */
    arguments->listFile = parseArgument2( argc, argv, "-files", &arg );

    if ( arguments->listFile ) {
      arguments->description = parseArgument2( argc, argv, "-desc", &arg );

      if ( arguments->description ) {
        arguments->units = parseArgument2( argc, argv, "-units", &arg );

        if ( arguments->units ) {

          if ( parseTimestampAndHours( argc, argv, &arg,
                                       &arguments->firstTimestamp,
                                       &arguments->hours ) ) {
            result = parseOptionalArguments( argc, argv, &arg, arguments );
          }
        }
      }
    }
  }

  if ( ! result ) {
    ZERO_OBJECT( arguments );
    printUsage();
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments )));
  return result;
}



/******************************************************************************
PURPOSE: initializeArguments - Initialize arguments.
INPUTS:  Arguments* arguments  Arguments to initialize.
OUTPUTS: Arguments* arguments  Arguments initialized.
******************************************************************************/

static void initializeArguments( Arguments* arguments ) {
  PRE0( arguments );
  ZERO_OBJECT( arguments );
  arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
  arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
  arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
  arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
}



/******************************************************************************
PURPOSE: parseOptionalArguments - Parse optional command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
         Integer* arg          Index of current argument to parse.
OUTPUTS: Integer* arg          Index of next argument to parse.
         Arguments* arguments  Updated arguments.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer parseOptionalArguments( Integer argc, char* argv[],
                                       Integer* arg,
                                       Arguments* arguments ) {

  PRE04( isValidArgs( argc, (const char**) argv ), arg, *arg > 0, arguments );

  Integer result = 1;
  Integer parsedBounds = 0;

  while ( *arg < argc ) {

    if ( AND2( ! strcmp( argv[ *arg ], "-domain" ), ! parsedBounds ) ) {
      parsedBounds = 1;
      result = parseBounds( argc, argv, arg, arguments->bounds );
    } else if ( AND2( ! strcmp( argv[ *arg ], "-corners" ),
                      ! arguments->corners ) ) {
      *arg += 1; /* Skip "-corners" option. */
      arguments->corners = 1;
    } else if ( AND2( ! strcmp( argv[ *arg ], "-daily" ),
                      ! arguments->daily ) ) {
      *arg += 1; /* Skip "-daily" option. */
      arguments->daily = 1;
    } else {
      failureMessage( "Invalid/redundant command-line argument: %s.",
                      argv[ *arg ] );
      result = 0;
      *arg = argc;
    }
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments )));
  return result;
}



/******************************************************************************
PURPOSE: readData - Read scan data from GOES files and subset it by time,
         lon-lat and data filtering.
INPUTS:  Data* data  data->arguments->listfile is the file containing the
         names of GOES files to read and subset by data->arguments->bounds.
OUTPTUS: Data* data  data->subsettedScans is the list of subset data to write.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void readData( Data* data ) {

  PRE05( data, data->ok, isValidArguments( &data->arguments ),
         data->scan.yyyydddhhmm == 0, data->subsettedScans == 0 );

  Stream* listFile = newFileStream( data->arguments.listFile, "r" );
  const Integer computeDailyMean = data->arguments.daily;
  Integer rows    = 0;
  Integer columns = 0;
  Integer* counts = 0;
  Real* means     = 0;
  Integer yyyyddd = 0;
  data->ok = 0;

  if ( listFile ) {
    const Integer firstTimestamp = data->arguments.firstTimestamp;
    const Integer lastTimestamp =
      offsetTimestamp( firstTimestamp, data->arguments.hours );
    Integer previousTimestamp = 0;

    /*
     * For each GOES data file,
     *   read its corresponding lonlat file and dimensions,
     *   optionally compute corners,
     *   read data into scan, subset scan by bounds, filter and append to list:
     */

    do {
      FileName fileName = "";
      memset( fileName, 0, sizeof fileName );
      listFile->readWord( listFile, fileName,
                          sizeof fileName / sizeof *fileName );

      if ( listFile->ok( listFile ) ) {
        const Integer currentTimestamp = timestampOfFileName( fileName );
        DEBUG( fprintf( stderr, "listing GOES file %s\n", fileName ); )

        if ( ! AND2( currentTimestamp > 0,
                     IMPLIES( previousTimestamp,
                              currentTimestamp > previousTimestamp ) ) ) {
          failureMessage( "Invalid/unordered GOES file %s.", fileName );
        } else if (IN_RANGE(currentTimestamp, firstTimestamp, lastTimestamp)) {
          readGOESDataFile( fileName, data );

          if ( data->ok ) {

            /* Read coordinates if scan dimensions do not match lonlats: */

            if ( OR2( data->scan.rows    != data->rows,
                      data->scan.columns != data->columns ) ) {
              readGOESCoordinatesFile( data );

              if ( data->ok ) {
                data->ok =
                  subsetIndicesByBounds( (const Real (*)[2])
                                           data->arguments.bounds,
                                         data->rows,
                                         data->columns,
                                         data->longitudes,
                                         data->latitudes,
                                         &data->subset[ ROW ][ MINIMUM ],
                                         &data->subset[ ROW ][ MAXIMUM ],
                                         &data->subset[ COLUMN ][ MINIMUM ],
                                         &data->subset[ COLUMN ][ MAXIMUM ] );
              }
            }

            if ( data->ok ) {

              if ( computeDailyMean ) {
                const Integer scanDay = data->scan.yyyydddhhmm / 10000;
                const Integer points = data->scan.rows * data->scan.columns;

                /*
                 * Reallocate if first time
                 * or if dimensions changed (discarding previous means):
                 */

                if ( OR2( rows != data->rows, columns != data->columns ) ) {
                  FREE( counts );
                  FREE( means );
                  rows = columns = 0;
                  counts = NEW_ZERO( Integer, points );
                  means  = counts ? NEW_ZERO( Real, points ) : 0;
                  data->ok = means != 0;

                  if ( ! data->ok ) {
                    FREE( counts );
                  } else {
                    rows    = data->scan.rows;
                    columns = data->scan.columns;
                    yyyyddd = scanDay;
                  }
                } else if ( yyyyddd != scanDay ) { /* Write previous day: */
                  CHECK4( rows == data->rows, columns == data->columns,
                          counts, means );
                  appendDailyMeans( yyyyddd * 10000, counts, means, data );
                  memset( counts, 0, points * sizeof counts[ 0 ] );
                  memset( means,  0, points * sizeof means[ 0 ] );
                  yyyyddd = scanDay;
                }

                if ( data->ok ) {
                  CHECK4( rows    == data->scan.rows,
                          columns == data->scan.columns,
                          counts, means );
                  computeMean( data, counts, means );
                }

              } else {
                SubsettedScan* subsettedScan =
                  subsetScan( (const Real (*)[2]) data->arguments.bounds,
                              (const Integer (*)[2]) data->subset,
                              data->longitudes, data->latitudes, data->corners,
                              &data->scan );

                previousTimestamp = currentTimestamp;

                if ( subsettedScan ) { /* Non-empty subset remains: */

                  if ( data->subsettedScans == 0 ) { /* Create list if needed*/
                    data->subsettedScans =
                      newVoidList( deallocateSubsettedScan, 0 );
                  }

                  if ( data->subsettedScans) { /* Append subset to list: */
                    data->subsettedScans->insert( data->subsettedScans,
                                                  subsettedScan, LAST_ITEM );

                    if ( data->subsettedScans->ok( data->subsettedScans ) ) {
                      data->ok = 1;
                    }
                  }
                }
              }
            }
          }
        }
      }

      listFile->readString( listFile, fileName, 2 ); /* Read '\n'. */
    } while ( ! listFile->isAtEnd( listFile ) );

    FREE_OBJECT( listFile );
  }

  if ( data->ok ) {

    if ( computeDailyMean ) {
      CHECK5( yyyyddd,
              rows == data->rows, columns == data->columns, counts, means );

      if ( data->subsettedScans == 0 ) {
        appendDailyMeans( yyyyddd * 10000, counts, means, data );
      } else {
        const SubsettedScan* const lastSubsettedScan =
          data->subsettedScans->item( data->subsettedScans, LAST_ITEM );
        const Integer lastScanDay = lastSubsettedScan->timestamp / 10000;

        if ( yyyyddd != lastScanDay ) {
          appendDailyMeans( yyyyddd * 10000, counts, means, data );
        }
      }

    } else if ( data->subsettedScans == 0 ) {
      failureMessage( "No scans were in the subset." );
      data->ok = 0;
    }
  }

  FREE( counts );
  FREE( means );
  POST02( isValidArguments( &data->arguments ),
          IMPLIES( data->ok, isValidData( data ) ) );
}



/******************************************************************************
PURPOSE: appendDailyMeans - Append a SubsettedScan of the daily means.
INPUTS:  const Integer yyyyddd0000  Timestamp of day.
         const Integer counts[ rows ][ columns ]  Counts of means.
         const Real    means[ rows ][ columns ]   Running means.
         Data* data                 Data.
OUTPTUS: Data* data  data->subsettedScans the list of subset data to write.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void appendDailyMeans( const Integer yyyyddd0000,
                              const Integer counts[], const Real means[],
                              Data* data ) {

  PRE04( isValidTimestamp( yyyyddd0000 ),
         data, data->ok, isValidArguments( &data->arguments ) );

  const Integer points = meanPoints( data, counts );
  CHECKING( const Integer OLD( count ) =
    data->subsettedScans ? data->subsettedScans->count( data->subsettedScans )
    : 0 );
  data->ok = 0;

  if ( points == 0 )  {
    failureMessage( "No scans were in the subset." );
  } else {
    SubsettedScan* subsettedScan = NEW_ZERO( SubsettedScan, 1 );

    if ( subsettedScan ) {
      const Real* const corners = data->corners;
      subsettedScan->timestamp = yyyyddd0000;
      subsettedScan->variables = 3 + ( corners != 0 ) * 2 * 4;
      subsettedScan->points = points;
      subsettedScan->data = NEW_ZERO( Real, subsettedScan->variables * points);

      if ( subsettedScan->data ) {
        Real* const subsetLongitudes   = subsettedScan->data;
        Real* const subsetLatitudes    = subsetLongitudes + points;
        Real* const subsetData         = subsetLatitudes  + points;
        Real* const subsetLongitudesSW = corners ? subsetData + points : 0;
        Real* const subsetLongitudesSE = corners ? subsetData + points * 2 : 0;
        Real* const subsetLongitudesNW = corners ? subsetData + points * 3 : 0;
        Real* const subsetLongitudesNE = corners ? subsetData + points * 4 : 0;
        Real* const subsetLatitudesSW  = corners ? subsetData + points * 5 : 0;
        Real* const subsetLatitudesSE  = corners ? subsetData + points * 6 : 0;
        Real* const subsetLatitudesNW  = corners ? subsetData + points * 7 : 0;
        Real* const subsetLatitudesNE  = corners ? subsetData + points * 8 : 0;
        copyMeanData( data->rows, data->columns,
                      (const Integer (*)[2]) data->subset,
                      points,
                      counts, means,
                      data->longitudes, data->latitudes, data->corners,
                      subsetLongitudes, subsetLatitudes, subsetData,
                      subsetLongitudesSW, subsetLongitudesSE,
                      subsetLongitudesNW, subsetLongitudesNE,
                      subsetLatitudesSW, subsetLatitudesSE,
                      subsetLatitudesNW, subsetLatitudesNE );
        CHECK( isValidSubsettedScan( subsettedScan ) );

        if ( data->subsettedScans == 0 ) {
          data->subsettedScans = newVoidList( deallocateSubsettedScan, 0 );
          data->ok = data->subsettedScans != 0;
        }

        if ( data->subsettedScans ) { /* Append subsetted scan to list: */
          data->subsettedScans->insert( data->subsettedScans,
                                        subsettedScan, LAST_ITEM );
          data->ok = data->subsettedScans->ok( data->subsettedScans );
        }
      }

      if ( ! data->ok ) {
        deallocateSubsettedScan( subsettedScan );
        FREE( subsettedScan );
      }
    }
  }

  POST02( isValidArguments( &data->arguments ),
          IMPLIES( data->ok,
                   AND2( isValidData( data ),
                         data->subsettedScans->count( data->subsettedScans )
                           == OLD( count ) + 1 ) ) );
}



/******************************************************************************
PURPOSE: readGOESDataFile - Read scan data from GOES data file.
INPUTS:  const char* fileName  Name of compressed GOES file to read.
         Data* data            data->buffer, data->bufferSize, data->variable
OUTPTUS: Data* data            data->buffer, data->bufferSize, data->variable
                               data->scan.yyyydddhhmm  YYYYDDDHHMM of file
                               data->rows, data->columns,
                               data->scan.data[ rows ][ columns ] scan data.
NOTES:  If successful data->ok is 1
        else data->ok is 0 and failureMessage() is called.
******************************************************************************/

static void readGOESDataFile( const char* fileName, Data* data ) {

  PRE04( fileName, data,
         IN_RANGE(  data->bufferSize, 10LL * 1024LL * 1024LL, INT_MAX / 2LL ),
         data->buffer );

  data->ok = readCompressedFile( fileName, (const int) data->bufferSize,
                                 (void*) data->buffer );

  if ( data->ok ) { /* Parse header and data: */
    Scan* const scan = &data->scan;
    const char* const dataPointer =
      parseHeader( data->buffer,
                   &scan->rows, &scan->columns, data->variable,
                   scan->coordinatesFileName );
    data->ok = dataPointer != 0;

    if ( data->ok ) {
      const Integer points = scan->rows * scan->columns;
      data->scan.data = NEW_ZERO( Real, points );
      data->ok = data->scan.data != 0;

      if ( data->ok ) {
        data->ok = parseData( dataPointer, points, scan->data, 0 );

        if ( data->ok ) {
          const Integer yyyydddhhmm = timestampOfFileName( fileName );
          CHECK( isValidTimestamp( yyyydddhhmm ) );
          scan->yyyydddhhmm = yyyydddhhmm;
          data->ok = findCoordinatesFile( fileName, scan->coordinatesFileName);
        }
      }
    }
  }

  POST0( IMPLIES( data->ok, isValidScan( &data->scan ) ) );
}



/******************************************************************************
PURPOSE: findCoordinatesFile - Find the correctly-pathed name of coordinates
         file that matches the given data file.
INPUTS:  const char* dataFileName  Name of GOES data file.
         FileName coordinatesFileName  Name of coordinates file.
OUTPTUS: FileName coordinatesFileName  Correctly-pathed name of file.
NOTES:  Integer 1 if file is found else 0 and failureMessage() is called.
******************************************************************************/

static Integer findCoordinatesFile( const char* const dataFileName,
                                    FileName coordinatesFileName ) {

  PRE04( dataFileName, *dataFileName,
         coordinatesFileName, *coordinatesFileName );

  Integer result = fileExists( coordinatesFileName );

  if ( ! result ) { /* Try compressed version of file name: */
    FileName name = "";
    memset( name, 0, sizeof name );
    strncpy( name, coordinatesFileName, sizeof name - 4 );
    CHECK( name[ sizeof name / sizeof *name - 4 ] == '\0' );
    strncat( name, ".gz", 3 );
    CHECK( name[ sizeof name / sizeof *name - 1 ] == '\0' );
    result = fileExists( name );

    if ( result ) {
      strncpy( coordinatesFileName, name,
               sizeof (FileName) / sizeof *coordinatesFileName - 1 );
    } else { /* Try adding the same path as the data file: */
      const char* const slash = strrchr( dataFileName, '/' );

      if ( slash ) {
        const char* c = dataFileName;
        Integer index = 0;
        FileName path = "";
        const Integer last = sizeof path / sizeof *path - 1;
        memset( path, 0, sizeof path );

        while ( AND2( c <= slash, index < last ) ) {
          path[ index ] = *c;
          ++index;
          ++c;
        }

        CHECK( path[ sizeof path / sizeof *path - 1 ] == '\0' );

        {
          FileName pathedName = "";
          const Integer pathLength = strlen( path );
          memset( pathedName, 0, sizeof pathedName / sizeof *pathedName );
          CHECK( pathLength < sizeof path / sizeof *path );

          /* Try path with compressed file name: */

          strncpy( pathedName, path, sizeof pathedName / sizeof *pathedName );
          CHECK(pathedName[sizeof pathedName / sizeof *pathedName -1] == '\0');
          strncat( pathedName, name,
                   sizeof pathedName /
                   sizeof *pathedName - strlen( pathedName ) );
          pathedName[ sizeof pathedName / sizeof *pathedName - 1 ] = '\0';
          result = fileExists( pathedName );

          if ( result ) {
            strncpy( coordinatesFileName, pathedName,
              sizeof (FileName) / sizeof *coordinatesFileName - 1 );
          } else { /* Try path with uncompressed file name: */
            strncpy( pathedName, path, sizeof pathedName / sizeof *pathedName );
            CHECK(pathedName[sizeof pathedName/sizeof *pathedName -1] == '\0');
            strncat( pathedName, coordinatesFileName,
                     sizeof pathedName /
                     sizeof *pathedName - strlen( pathedName ) );
            pathedName[ sizeof pathedName / sizeof *pathedName - 1 ] = '\0';
            result = fileExists( pathedName );

            if ( result ) {
              strncpy( coordinatesFileName, pathedName,
                       sizeof (FileName) / sizeof *coordinatesFileName - 1 );
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    failureMessage( "Failed to find coordinates file '%s'.",
                    coordinatesFileName );
  }

  CHECK( coordinatesFileName[
          sizeof (FileName) / sizeof *coordinatesFileName - 1 ]
         == '\0' );

  DEBUG( fprintf( stderr, "coordinatesFileName = '%s'.\n",
                  coordinatesFileName ); )

  POST02( IS_BOOL( result ),
          IMPLIES( result, fileExists( coordinatesFileName ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readGOESCoordinatesFile - Read scan lonlats from GOES lonlat file.
INPUTS:  const char* fileName  Name of compressed GOES file to read.
         Data* data            data->buffer, data->bufferSize,
OUTPTUS: Data* data            data->buffer, data->bufferSize,
                               data->scan.yyyydddhhmm  YYYYDDDHHMM of file
                               data->rows, data->columns,
                               data->scan.data[ rows ][ columns ] scan data.
NOTES:  If successful data->ok is 1
        else data->ok is 0 and failureMessage() is called.
******************************************************************************/

static void readGOESCoordinatesFile( Data* data ) {

  PRE04( data, data->scan.coordinatesFileName,
         IN_RANGE(  data->bufferSize, 64LL * 1024LL * 1024LL, INT_MAX / 2LL ),
         data->buffer );

  data->ok = readCompressedFile( data->scan.coordinatesFileName,
                                 (const int) data->bufferSize,
                                 (void*) data->buffer );

  if ( data->ok ) { /* Parse header and data: */
    const char* const dataPointer =
      parseHeader( data->buffer, &data->rows, &data->columns, 0, 0 );
    data->ok = dataPointer != 0;

    if ( data->ok ) {
      const Integer points = data->rows * data->columns;
      FREE( data->longitudes );
      FREE( data->latitudes );
      FREE( data->corners );
      data->longitudes = NEW_ZERO( Real, points );
      data->latitudes = data->longitudes ? NEW_ZERO( Real, points ) : 0;
      data->ok = data->latitudes != 0;

      if ( data->ok ) {
        data->ok =
          parseData( dataPointer, data->rows * data->columns,
                     data->longitudes, data->latitudes );

        if ( data->ok ) {

          if ( data->arguments.corners ) {
            data->corners = NEW_ZERO( Real, 2 * 4 * points );
            data->ok = data->corners != 0;

            if ( data->ok ) {
              computeCorners( data->rows, data->columns,
                              data->longitudes, data->latitudes,
                              data->corners );
            }
          }
        }
      }
    }
  }

  POST0( IMPLIES( data->ok, isValidScan( &data->scan ) ) );
}



/******************************************************************************
PURPOSE: parseHeader - Parse dimensions from a GOES header.
INPUTS:  const char* buffer  String contents of a GOES file.
OUTPTUS: Integer* rows       Number of data rows.
         Integer* columns    Number of data columns.
         Variable variable   Or 0. If empty then init else verify match.
         FileName coordinatesFileName  If not 0 then
                                       Name of corresponding coordinates file.
RETURNS: const char* pointer into buffer of start of data or 0 if failed to
         parse valid header content and failureMessage() is called.
******************************************************************************/

static const char* parseHeader( const char* buffer,
                                Integer* rows, Integer* columns,
                                Variable variable,
                                FileName coordinatesFileName ) {
  PRE03( buffer, rows, columns );
  const char* result = 0;
  const Integer headerLines = parseCount( buffer, "hdr_lines" );

  if ( headerLines > 0 ) {
    *columns = parseCount( buffer, "NX" );

    if ( *columns > 0 ) {
      *rows = parseCount( buffer, "NY" );

      if ( *rows * *columns > 0 ) {
        result = skipLines( buffer, headerLines );

        if ( result ) {

          if ( variable ) {
            const char* const name = parseWord( buffer, "product" );

            if ( name ) {

              if ( *variable == '\0' ) { /* Initialize variable: */
                memset( variable, 0, sizeof (Variable) );
                copyWord( variable, name,
                          sizeof (Variable) / sizeof *variable - 1 );
              } else if ( strncmp( variable, name, strlen( variable ) ) ) {
                failureMessage( "Data file product name does not match "
                                "expected variable name '%s'.",
                                variable ); 
                result = 0;
              }
            }
          }

          if ( AND2( result, coordinatesFileName ) ) {
            const char* const name = parseWord( buffer, "georef" );

            if ( name ) {
              memset( coordinatesFileName, 0, sizeof (FileName) );
              copyWord( coordinatesFileName, name,
                        sizeof (FileName) / sizeof *coordinatesFileName - 1 );
            } else {
              result = 0;
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    *rows = *columns = 0;

    if ( variable ) {
      memset( variable, 0, sizeof (Variable) );
    }

    if ( coordinatesFileName ) {
      memset( coordinatesFileName, 0, sizeof (FileName) );
    }
  }

  POST0( IMPLIES_ELSE( result,
                       AND3( GT_ZERO3( *rows, *columns, *rows * *columns ),
                             IMPLIES( coordinatesFileName,
                                      *coordinatesFileName ),
                             IMPLIES( variable, *variable ) ),
                       AND3( IS_ZERO2( *rows, *columns ),
                             IMPLIES( coordinatesFileName,
                                      *coordinatesFileName == '\0' ),
                             IMPLIES( variable, *variable == '\0' ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseCount - Read an integer value > 0 for given name in a string.
INPUTS:  const char* header   GOES header string to parse.
         const char* name     Name of value to parse.
RETURNS: Integer value > 0 or 0 if failed and failureMessage() is called.
******************************************************************************/

static Integer parseCount( const char* header, const char* name ) {
  PRE03( header, name, *name );
  Integer result = 0;
  const char* const line = strstr( header, name );

  if ( line ) {
    const Integer ok =
      AND2( sscanf( line, "%*s : %lld", &result ) == 1, result > 0 );

    if ( ! ok ) {
      result = 0;
    }
  }

  if ( ! result ) {
    failureMessage( "Failed to read valid count for %s in header.", name );
  }

  return result;
}



/******************************************************************************
PURPOSE: parseWord - Find a word value for given name in a string.
INPUTS:  const char* header   GOES header string to parse.
         const char* name     Name of value to parse.
RETURNS: const char* start of word or 0 if failed and failureMessage() called.
******************************************************************************/

static const char* parseWord( const char* header, const char* name ) {
  PRE03( header, name, *name );
  const char* result = 0;
  const char* const line = strstr( header, name );

  if ( line ) {
    const char* const colon = strchr( line, ':' );

    if ( colon ) {
      result = colon + 1;

      while ( isspace( *result ) ) {
        ++result;
      }
    }
  }

  if ( ! result ) {
    failureMessage( "Failed to read valid word for %s in header.", name );
  }

  return result;
}



/******************************************************************************
PURPOSE: parseData - Read real data values from GOES file content string.
INPUTS:  const char* content  GOES data string to parse.
         const Integer count  Number of values to parse.
OUTPUTS: Real data[ count ]   Data values read.
         Real data2[ count ]  0 or 2nd data values read (if lonlats).
RETURNS: Integer 1 if valid values or 0 if failed and failureMessage() called.
******************************************************************************/

static Integer parseData( const char* buffer, const Integer count,
                          Real data[], Real data2[] ) {
  PRE03( buffer, count > 0, data );
  Integer result = 0;
  Integer read = 0;
  const char* line = buffer;

  if ( data2 ) {

    do {
#if 0
      /* HACK: sscanf() calls strlen() which is prohibitively slow! */

      result =
        AND3( sscanf( line, "%*s %*s %lf %lf\n",
                      data + read, data2 + read ) == 2,
              isValidLongitude( data[ read ] ),
              isValidLatitude( data2[ read ] ) );
      read += result;
      line = strchr( line + 1, '\n' );
#else
      line = readAndSkipReal( line, data + read );
      line = readAndSkipReal( line, data + read );
      line = readAndSkipReal( line, data2 + read );
      line = readAndSkipReal( line, data + read );
      result = AND2( isValidLongitude( data[ read ] ),
                     isValidLatitude( data2[ read ] ) );
      read += result;
#endif
    } while ( AND3( result, line, read < count ) ); 

  } else {

    do {
      data[ read ] = atof( line );
      result = data[ read ] >= MISSING;
#if 1
      /* HACK: some files contain Nans! Replace them with MISSING: */

      if ( ! result ) {
        data[ read ] = MISSING;
        result = 1;
      }
#endif
      read += result;
      line = strchr( line + 1, '\n' );
    } while ( AND3( result, line, read < count ) ); 
  }

  if ( ! result ) {
    memset( data, 0, count * sizeof *data );

    if ( data2 ) {
      memset( data2, 0, count * sizeof *data2 );      
    }
  }

  POST0( IMPLIES( result,
                  IMPLIES_ELSE( data2,
                                validLongitudesAndLatitudes(count, data,data2),
                                isNanFree( data, count ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: timestampOfFileName - Timestamp of GOES file name.
INPUTS:  const char* fileName  Name of GOES file to parse. E.g.,
                               "testdata/goes_snd_TSKN_v00_201309012345Z.txt"
RETURNS: Integer YYYYDDDHHMM. E.g., 20132442345, or 0 if failed and
         failureMessage() is called.
******************************************************************************/

static Integer timestampOfFileName( const char* fileName ) {
  PRE0( fileName );
  Integer result = 0;
  const char* timestamp = strrchr( fileName, '_' );

  if ( timestamp ) {
    int yyyy = 0;
    int mo = 0;
    int dd = 0;
    int hh = 0;
    int mm = 0;

    if ( sscanf( timestamp, "_%4d%2d%2d%2d%2dZ",
                 &yyyy, &mo, &dd, &hh, &mm ) == 5 ) {
      const Integer yyyymmdd = yyyy * 10000 + mo * 100 + dd;

      if ( AND3( isValidYearMonthDay( yyyymmdd ),
                 IN_RANGE( hh, 0, 23 ), IN_RANGE( mm, 0, 59 ) ) ) {
        const Integer yyyyddd = convertYearMonthDay( yyyymmdd );
        result = yyyyddd * 10000 + hh * 100 + mm;
      }
    }
  }

  if ( ! result ) {
    failureMessage( "Invalid file name %s.", fileName );
  }

  POST0( IMPLIES( result, isValidTimestamp( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: subsetScan - Subset scan by bounds and data filtering.
INPUTS:  const Bounds bounds         Subset longitude-latitude bounds.
         const Integer subset[2][2]  Subset indices[ ROW COLUMN][MIN/MAXIMUM].
         const Real longitudes[rows][columns]     Longitudes.
         const Real latitudes[rows][columns]      Latitudes.
         const Real corners[2][4][rows][columns]  lonlat corners or 0.
         Scan* scan                               Uncompressed scan.
OUTPUTS: Scan* scan  scan->data[ rows ][ columns ] Zero'd by filter.
RETURNS: SubsettedScan* if successful, else 0 if no points are in the subset or
         there was a memory allocation failure and failureMessage() is called.
NOTES:   Must use FREE() on returned result when finished.
******************************************************************************/

static SubsettedScan* subsetScan( const Bounds bounds,
                                  const Integer subset[ 2 ][ 2 ],
                                  const Real longitudes[],
                                  const Real latitudes[],
                                  const Real corners[],
                                  Scan* scan ) {

  PRE010( isValidBounds( bounds ),
          isValidScan( scan ),
          IN_RANGE( subset[ ROW ][ MINIMUM ], 0, scan->rows - 1 ),
          IN_RANGE( subset[ ROW ][ MAXIMUM ],
                    subset[ ROW ][ MINIMUM ], scan->rows - 1 ),
          IN_RANGE( subset[ COLUMN ][ MINIMUM ], 0, scan->columns - 1 ),
          IN_RANGE( subset[ COLUMN ][ MAXIMUM ],
                    subset[ COLUMN ][ MINIMUM ], scan->columns - 1 ),
          longitudes,
          latitudes,
          validLongitudesAndLatitudes( scan->rows * scan->columns,
                                       longitudes, latitudes ),
          IMPLIES( corners,
                   validLongitudesAndLatitudes( scan->rows * scan->columns * 4,
                                                corners,
                                                corners +
                                                  scan->rows *
                                                  scan->columns * 4 ) ) );

  SubsettedScan* result = 0;
  Integer filteredSubset[ 2 ][ 2 ] = { { 0, 0 }, { 0, 0 } };
  const void* const unused_ =
    memcpy( filteredSubset, subset, sizeof filteredSubset );
  const Integer points =
    subsetScanCount( bounds, longitudes, latitudes, filteredSubset, scan );

  if ( points ) { /* If not all points are filtered-out, copy subset of them:*/
    result = NEW_ZERO( SubsettedScan, 1 );

    if ( result ) {
      const Integer variables = 3 + ( corners != 0 ) * 4 * 2; /*lon,lat,var,*/
      CHECK3( variables >= 3, points > 0, variables * points >= 3 );
      result->timestamp = scan->yyyydddhhmm;
      result->variables = variables;
      result->points    = points;
      result->data = NEW_ZERO( Real, variables * points );

      if ( result->data ) {
        Real* const subsetLongitudes   = result->data;
        Real* const subsetLatitudes    = subsetLongitudes + points;
        Real* const subsetData         = subsetLatitudes  + points;
        Real* const subsetLongitudesSW = corners ? subsetData + points : 0;
        Real* const subsetLongitudesSE =
          corners ? subsetLongitudesSW + points : 0;
        Real* const subsetLongitudesNW =
          corners ? subsetLongitudesSE + points : 0;
        Real* const subsetLongitudesNE =
          corners ? subsetLongitudesNW + points : 0;
        Real* const subsetLatitudesSW  =
          corners ? subsetLongitudesNE + points : 0;
        Real* const subsetLatitudesSE  =
          corners ? subsetLatitudesSW + points : 0;
        Real* const subsetLatitudesNW  =
          corners ? subsetLatitudesSE + points : 0;
        Real* const subsetLatitudesNE  =
          corners ? subsetLatitudesNW + points : 0;
        const Integer firstRow    = filteredSubset[ ROW    ][ MINIMUM ];
        const Integer lastRow     = filteredSubset[ ROW    ][ MAXIMUM ];
        const Integer firstColumn = filteredSubset[ COLUMN ][ MINIMUM ];
        const Integer lastColumn  = filteredSubset[ COLUMN ][ MAXIMUM ];

        CHECK( IN_RANGE( ( 1 + lastRow - firstRow ) *
                         ( 1 + lastColumn - firstColumn ),
                          points, scan->rows * scan->columns ) );

        copySubsetLongitudesAndLatitudes( scan->rows,
                                          scan->columns,
                                          longitudes,
                                          latitudes,
                                          corners,
                                          scan->data,
                                          points,
                                          firstRow, lastRow,
                                          firstColumn, lastColumn,
                                          subsetLongitudes,
                                          subsetLatitudes,
                                          subsetLongitudesSW,
                                          subsetLongitudesSE,
                                          subsetLongitudesNW,
                                          subsetLongitudesNE,
                                          subsetLatitudesSW,
                                          subsetLatitudesSE,
                                          subsetLatitudesNW,
                                          subsetLatitudesNE );

        copySubsetData( scan->rows, scan->columns, scan->data,
                        firstRow, lastRow, firstColumn, lastColumn,
                        subsetData );
      }
    }
  }

  POST0( IMPLIES( result,
                  AND6( isValidScan( scan ),
                        isValidSubsettedScan( result ),
                        minimumItem( result->data, result->points ) >=
                          bounds[ LONGITUDE ][ MINIMUM ],
                        maximumItem( result->data, result->points ) <=
                          bounds[ LONGITUDE ][ MAXIMUM ],
                        minimumItem( result->data + result->points,
                                     result->points ) >=
                          bounds[ LATITUDE ][ MINIMUM ],
                        maximumItem( result->data + result->points,
                                     result->points ) <=
                          bounds[ LATITUDE ][ MAXIMUM ] ) ) );

  return result;
}



/******************************************************************************
PURPOSE: subsetScanCount - Count scan points subsetted by indices and ranges.
INPUTS:  const Bounds bounds                  Subset lon-lat bounds.
         const Real ranges[ VARIABLES ][ 2 ]  Data ranges to filter by.
         Integer subset[ 2 ][ 2 ]             subset[ROW COLUMN][MIN/MAXIMUM].
         const Real longitudes[ scan->rows ][ scan->rows ]  Longitudes.
         const Real latitudes[ scan->rows ][ scan->rows ]   Latitudes.
         Scan* scan                           Uncompressed & decoded scan.
OUTPUTS: Integer subset[ 2 ][ 2 ]  Subset indices reduced by bounds/filtering.
RETURNS: Integer number of points remaining after subsetting/filtering.
******************************************************************************/

static Integer subsetScanCount( const Bounds bounds,
                                const Real longitudes[],
                                const Real latitudes[],
                                Integer subset[ 2 ][ 2 ],
                                Scan* scan ) {

  PRE08( isValidBounds( bounds ),
         isValidScan( scan ),
         subset,
         IN_RANGE( subset[ ROW ][ MINIMUM ], 0, scan->rows - 1 ),
         IN_RANGE( subset[ ROW ][ MAXIMUM ],
                   subset[ ROW ][ MINIMUM ], scan->rows - 1 ),
         IN_RANGE( subset[ COLUMN ][ MINIMUM ], 0, scan->columns - 1 ),
         IN_RANGE( subset[ COLUMN ][ MAXIMUM ],
                   subset[ COLUMN ][ MINIMUM ], scan->columns - 1 ),
         validLongitudesAndLatitudes( scan->rows * scan->columns,
                                      longitudes, latitudes ) );

  Real* const data = scan->data;
  const Integer rows          = scan->rows;
  const Integer columns       = scan->columns;
  const Integer firstRow      = subset[ ROW    ][ MINIMUM ];
  const Integer lastRow       = subset[ ROW    ][ MAXIMUM ];
  const Integer firstColumn   = subset[ COLUMN ][ MINIMUM ];
  const Integer lastColumn    = subset[ COLUMN ][ MAXIMUM ];
  const Real longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const Real longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const Real latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const Real latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  Integer result = 0;
  Integer row = 0;

  /* Compute reduced mask in parallel: */

#pragma omp parallel for reduction( + : result )

  for ( row = firstRow; row <= lastRow; ++row ) {
    const Integer rowOffset = row * columns;
    Integer column = 0;
    CHECK( rowOffset + lastColumn < scan->rows * scan->columns );

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Real longitude = longitudes[ index ];
      const Real latitude  = latitudes[ index ];
      const Real value = data[ index ];
      const Integer output =
        AND3( value >= 0.0,
              IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
              IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );
      result += output;
      data[ index ] = output ? value : MISSING;
    }
  }

  /* Compute reduced subset sequentially: */

  if ( result ) {
    reduceSubset( rows, columns, data,
                  &subset[ ROW ][ MINIMUM ],
                  &subset[ ROW ][ MAXIMUM ],
                  &subset[ COLUMN ][ MINIMUM ],
                  &subset[ COLUMN ][ MAXIMUM ] );

    CHECK4( IN_RANGE( subset[ ROW ][ MINIMUM ], firstRow, lastRow ),
            IN_RANGE( subset[ ROW ][ MAXIMUM ],
                      subset[ ROW ][ MINIMUM ], lastRow ),
            IN_RANGE( subset[ COLUMN ][ MINIMUM ], firstColumn, lastColumn ),
            IN_RANGE( subset[ COLUMN ][ MAXIMUM ],
                      subset[ COLUMN ][ MINIMUM ], lastColumn ) );
  }

  POST06( result >= 0,
          isValidScan( scan ),
          IN_RANGE( subset[ ROW ][ MINIMUM ], 0, scan->rows - 1 ),
          IN_RANGE( subset[ ROW ][ MAXIMUM ],
                    subset[ ROW ][ MINIMUM ], scan->rows - 1 ),
          IN_RANGE( subset[ COLUMN ][ MINIMUM ], 0, scan->columns - 1 ),
          IN_RANGE( subset[ COLUMN ][ MAXIMUM ],
                    subset[ COLUMN ][ MINIMUM ], scan->columns - 1 ) );

  return result;
}



/******************************************************************************
PURPOSE: subsetIndicesByBounds - Subset row and column indices by bounds.
INPUTS:  const Bounds bounds    Longitude-latitude bounds of subset.
         const Integer rows     Number of row points.
         const Integer columns  Number of column points.
         const Real longitudes[ rows ][ columns ]  Longitudes.
         const Real latitudes[  rows ][ columns ]  Latitudes.
OUTPUTS: Integer* firstRow      Index of first row of subset.
         Integer* lastRow       Index of last  row of subset.
         Integer* firstColumn   Index of first column of subset.
         Integer* lastColumn    Index of last  column of subset.
RETURNS: Integer 1 if a non-empty subset exists, else 0 and outputs are -1.
******************************************************************************/

static Integer subsetIndicesByBounds( const Bounds bounds,
                                      const Integer rows,
                                      const Integer columns,
                                      const Real longitudes[],
                                      const Real latitudes[],
                                      Integer* firstRow,
                                      Integer* lastRow,
                                      Integer* firstColumn,
                                      Integer* lastColumn ) {

  PRE010( bounds, isValidBounds( bounds ),
          rows > 0, columns > 0, rows * columns > 0,
          validLongitudesAndLatitudes( rows * columns, longitudes, latitudes ),
          firstRow, lastRow, firstColumn, lastColumn );

  const Real longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const Real longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const Real latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const Real latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  Integer theFirstRow    = 0;
  Integer theLastRow     = 0;
  Integer theFirstColumn = 0;
  Integer result         = 0;
  Integer row            = 0;
  Integer column         = 0;

  *firstRow = *lastRow = *firstColumn = *lastColumn = -1;

  /* Loop forward through rows to find first subset row: */

  for ( row = 0; row < rows; ++row ) {
    const Integer rowOffset = row * columns;
    CHECK( rowOffset + columns - 1 < rows * columns );

    for ( column = 0; column < columns; ++column ) {
      const Integer index = rowOffset + column;
      const Real longitude = longitudes[ index ];
      const Real latitude  = latitudes[ index ];
      const Integer inside =
        AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
              IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

      if ( inside ) {
        *firstRow = theFirstRow = row;
        result = 1;
        row = rows;
        column = columns;
      }
    }
  }

  if ( result ) {

    /* Loop backward through rows to find last subset row: */

    for ( row = rows - 1; row >= theFirstRow; --row ) {
      const Integer rowOffset = row * columns;
      CHECK( rowOffset + columns - 1 < rows * columns );

      for ( column = 0; column < columns; ++column ) {
        const Integer index = rowOffset + column;
        const Real longitude = longitudes[ index ];
        const Real latitude  = latitudes[ index ];
        const Integer inside =
          AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

        if ( inside ) {
          *lastRow = theLastRow = row;
          row = 0;
          column = columns;
        }
      }
    }

    CHECK2( IN_RANGE( *firstRow, 0, rows - 1 ),
            IN_RANGE( *lastRow, *firstRow, rows - 1 ) );


    /* Loop forward through columns to find first subset column: */

    for ( column = 0; column < columns; ++column ) {

      for ( row = theFirstRow; row <= theLastRow; ++row ) {
        const Integer index = row * columns + column;
        const Real longitude = longitudes[ index ];
        const Real latitude  = latitudes[ index ];
        const Integer inside =
          AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

        if ( inside ) {
          *firstColumn = theFirstColumn = column;
          row = rows;
          column = columns;
        }
      }
    }

    /* Loop backward through columns to find last subset column: */

    for ( column = columns - 1; column >= theFirstColumn; --column ) {

      for ( row = theFirstRow; row <= theLastRow; ++row ) {
        const Integer index = row * columns + column;
        const Real longitude = longitudes[ index ];
        const Real latitude  = latitudes[ index ];
        const Integer inside =
          AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

        if ( inside ) {
          *lastColumn = column;
          row = rows;
          column = 0;
        }
      }
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                       AND4( IN_RANGE( *firstRow, 0, rows - 1 ),
                             IN_RANGE( *lastRow, *firstRow, rows - 1 ),
                             IN_RANGE( *firstColumn, 0, columns - 1 ),
                             IN_RANGE( *lastColumn, *firstColumn, columns -1)),
                       AND4( *firstRow == -1, *lastRow == -1,
                             *firstColumn == -1, *lastColumn == -1 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: reduceSubset - Reduce subset row and column indices by MISSING data.
INPUTS:  const Real data[ rows ][ columns ]  Subset-modified data.
         Integer*   firstRow     Index of first row of bounds subset.
         Integer*   lastRow      Index of alst row of bounds subset.
         Integer*   firstColumn  Index of first column of bounds subset.
         Integer*   lastColumn   Index of last column of bounds subset.
OUTPUTS: Integer*   firstRow     Index of first row of mask-filtered subset.
         Integer*   lastRow      Index of last row of mask-filtered subset.
         Integer*   firstColumn  Index of first column of mask-filtered subset.
         Integer*   lastColumn   Index of last column of mask-filtered subset.
******************************************************************************/

static void reduceSubset( const Integer rows, const Integer columns,
                          const Real data[],
                          Integer* firstRow, Integer* lastRow,
                          Integer* firstColumn, Integer* lastColumn) {

  PRE012( rows > 0,
          columns > 0,
          rows * columns > 0,
          isNanFree( data, rows * columns ),
          firstRow,
          lastRow,
          firstColumn,
          lastColumn,
          IN_RANGE( *firstRow, 0, rows - 1 ),
          IN_RANGE( *lastRow, *firstRow, rows - 1 ),
          IN_RANGE( *firstColumn, 0, columns - 1 ),
          IN_RANGE( *lastColumn, *firstColumn, columns - 1 ) );

  CHECKING( const Integer OLD( firstRow ) = *firstRow; )
  CHECKING( const Integer OLD( lastRow ) = *lastRow; )
  CHECKING( const Integer OLD( firstColumn ) = *firstColumn; )
  CHECKING( const Integer OLD( lastColumn ) = *lastColumn; )
  Integer theFirstRow    = *firstRow;
  Integer theLastRow     = *lastRow;
  Integer theFirstColumn = *firstColumn;
  Integer theLastColumn  = *lastColumn;
  Integer row = 0;
  Integer column = 0;

  /* Loop forward through rows to find first subset row: */

  for ( row = theFirstRow; row <= theLastRow; ++row ) {
    const Integer rowOffset = row * columns;
    CHECK( rowOffset + theLastColumn < rows * columns );

    for ( column = theFirstColumn; column <= theLastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Real value = data[ index ];

      if ( value >= 0.0 ) {
        theFirstRow = row;
        row = rows;
        column = columns;
      }
    }
  }

  /* Loop backward through rows to find last subset row: */

  for ( row = theLastRow; row >= theFirstRow; --row ) {
    const Integer rowOffset = row * columns;
    CHECK( rowOffset + theLastColumn < rows * columns );

    for ( column = theFirstColumn; column <= theLastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Real value = data[ index ];

      if ( value >= 0.0 ) {
        theLastRow = row;
        row = 0;
        column = columns;
      }
    }
  }

  CHECK2( IN_RANGE( theFirstRow, *firstRow, *lastRow ),
          IN_RANGE( theLastRow, theFirstRow, *lastRow ) );

  /* Loop forward through columns to find first subset column: */

  for ( column = theFirstColumn; column <= theLastColumn; ++column ) {

    for ( row = theFirstRow; row <= theLastRow; ++row ) {
      const Integer index = row * columns + column;
      const Real value = data[ index ];

      if ( value >= 0.0 ) {
        theFirstColumn = column;
        row = rows;
        column = columns;
      }
    }
  }

  /* Loop backward through columns to find last subset column: */

  for ( column = theLastColumn; column >= theFirstColumn; --column ) {

    for ( row = theFirstRow; row <= theLastRow; ++row ) {
      const Integer index = row * columns + column;
      const Real value = data[ index ];

      if ( value >= 0.0 ) {
        theLastColumn = column;
        row = rows;
        column = 0;
      }
    }
  }

  *firstRow    = theFirstRow;
  *lastRow     = theLastRow;
  *firstColumn = theFirstColumn;
  *lastColumn  = theLastColumn;

  POST04( IN_RANGE( *firstRow, OLD( firstRow ), OLD( lastRow ) ),
          IN_RANGE( *lastRow, *firstRow, OLD( lastRow ) ),
          IN_RANGE( *firstColumn, OLD( firstColumn ), OLD( lastColumn ) ),
          IN_RANGE( *lastColumn, *firstColumn, OLD( lastColumn ) ) );
}



/******************************************************************************
PURPOSE: copySubsetLongitudesAndLatitudes - Copy subsetted/filtered longitudes
         and latitudes.
INPUTS:  const Integer rows                  Number of data rows.
         const Integer columns               Number of data columns.
         const Real longitudes[ rows ][ columns ]  Data longitudes.
         const Real latitudes[  rows ][ columns ]  Data latitudes.
         const Real corners[ 2 ][ 4 ][  rows ][ columns ]  Lonlat corners or 0.
         const Real data[ rows ][ columns ]  Subset-modified data.
         Integer    points       Number of subset points.
         Integer    firstRow     Index of first row of subset.
         Integer    lastRow      Index of alst row os subset.
         Integer    firstColumn  Index of first column of subset.
         Integer    lastColumn   Index of last column of subset.
OUTPUTS: Real subsetLongitudes[ points ]  Longitudes to append to.
         Real subsetLatitudes[ points ]   Latitudes to append to.
         Real subsetLongitudesSW[ points ]  Optional: for corners or 0.
         Real subsetLongitudesSE[ points ]  Optional: for corners or 0.
         Real subsetLongitudesNW[ points ]  Optional: for corners or 0.
         Real subsetLongitudesNE[ points ]  Optional: for corners or 0.
         Real subsetLatitudesSW[ points ]   Optional: for corners or 0.
         Real subsetLatitudesSE[ points ]   Optional: for corners or 0.
         Real subsetLatitudesNW[ points ]   Optional: for corners or 0.
         Real subsetLatitudesNE[ points ]   Optional: for corners or 0.
NOTES:   Uses static global longitudesLatitudes.
******************************************************************************/

static void copySubsetLongitudesAndLatitudes( const Integer rows,
                                              const Integer columns,
                                              const Real longitudes[],
                                              const Real latitudes[],
                                              const Real corners[],
                                              const Real data[],
                                              Integer points,
                                              Integer firstRow,
                                              Integer lastRow,
                                              Integer firstColumn,
                                              Integer lastColumn,
                                              Real subsetLongitudes[],
                                              Real subsetLatitudes[],
                                              Real subsetLongitudesSW[],
                                              Real subsetLongitudesSE[],
                                              Real subsetLongitudesNW[],
                                              Real subsetLongitudesNE[],
                                              Real subsetLatitudesSW[],
                                              Real subsetLatitudesSE[],
                                              Real subsetLatitudesNW[],
                                              Real subsetLatitudesNE[] ) {

  PRE016( rows > 0,
          columns > 0,
          rows * columns > 0,
          longitudes,
          latitudes,
          validLongitudesAndLatitudes( rows * columns, longitudes, latitudes ),
          IMPLIES( corners,
                   validLongitudesAndLatitudes( rows * columns * 4,
                                                corners,
                                                corners + rows * columns * 4)),
          isNanFree( data, rows * columns ),
          IN_RANGE( points, 1,
                    (1 + lastRow - firstRow) * (1 + lastColumn - firstColumn)),
          IN_RANGE( firstRow, 0, rows - 1 ),
          IN_RANGE( lastRow, firstRow, rows - 1 ),
          IN_RANGE( firstColumn, 0, columns - 1 ),
          IN_RANGE( lastColumn, firstColumn, columns - 1 ),
          subsetLongitudes,
          subsetLatitudes,
          IMPLIES_ELSE( subsetLongitudesSW,
                        NON_ZERO7( subsetLongitudesSE,
                                   subsetLongitudesNW,
                                   subsetLongitudesNE,
                                   subsetLatitudesSW,
                                   subsetLatitudesSE,
                                   subsetLatitudesNW,
                                   subsetLatitudesNE ),
                        IS_ZERO7( subsetLongitudesSE,
                                  subsetLongitudesNW,
                                  subsetLongitudesNE,
                                  subsetLatitudesSW,
                                  subsetLatitudesSE,
                                  subsetLatitudesNW,
                                  subsetLatitudesNE ) ) );

  CHECKING( Integer count = 0; )
  Integer row = 0;
  Integer output = 0;
  const Integer rowsTimesColumns = rows * columns;
  const Real* const longitudesSW = corners;
  const Real* const longitudesSE = corners ? corners + rowsTimesColumns : 0;
  const Real* const longitudesNW = corners ? corners + rowsTimesColumns * 2 :0;
  const Real* const longitudesNE = corners ? corners + rowsTimesColumns * 3 :0;
  const Real* const latitudesSW  = corners ? corners + rowsTimesColumns * 4 :0;
  const Real* const latitudesSE  = corners ? corners + rowsTimesColumns * 5 :0;
  const Real* const latitudesNW  = corners ? corners + rowsTimesColumns * 6 :0;
  const Real* const latitudesNE  = corners ? corners + rowsTimesColumns * 7 :0;

  for ( row = firstRow; row <= lastRow; ++row ) {
    const Integer rowOffset = row * columns;
    Integer column = 0;
    CHECK( rowOffset + lastColumn < rows * columns );

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Real value = data[ index ];

      if ( value >= 0.0 ) {
        const Real longitude = longitudes[ index ];
        const Real latitude  = latitudes[ index ];
        CHECKING( ++count; )
        CHECK3( IN_RANGE( longitude, -180.0, 180.0 ),
                IN_RANGE( latitude, -90.0, 90.0 ),
                IN_RANGE( count, 1, points ) );
        subsetLongitudes[ output ] = longitude;
        subsetLatitudes[  output ] = latitude;

        if ( subsetLongitudesSW ) {
          subsetLongitudesSW[ output ] = longitudesSW[ index ];
          subsetLongitudesSE[ output ] = longitudesSE[ index ];
          subsetLongitudesNW[ output ] = longitudesNW[ index ];
          subsetLongitudesNE[ output ] = longitudesNE[ index ];
          subsetLatitudesSW[  output ] = latitudesSW[ index ];
          subsetLatitudesSE[  output ] = latitudesSE[ index ];
          subsetLatitudesNW[  output ] = latitudesNW[ index ];
          subsetLatitudesNE[  output ] = latitudesNE[ index ];
        }

        ++output;
      }
    }
  }

  POST02( validLongitudesAndLatitudes( points,
                                       subsetLongitudes, subsetLatitudes ),
         IMPLIES( subsetLongitudesSW,
                  AND4( validLongitudesAndLatitudes( points,
                                                     subsetLongitudesSW,
                                                     subsetLatitudesSW ),
                        validLongitudesAndLatitudes( points,
                                                     subsetLongitudesSE,
                                                     subsetLatitudesSE ),
                        validLongitudesAndLatitudes( points,
                                                     subsetLongitudesNW,
                                                     subsetLatitudesNW ),
                        validLongitudesAndLatitudes( points,
                                                     subsetLongitudesNE,
                                                     subsetLatitudesNE ) ) ) );
}



/******************************************************************************
PURPOSE: copySubsetData - Copy subsetted/filtered variable data.
INPUTS:  const Integer rows                  Number of data rows.
         const Integer columns               Number of data columns.
         const Real data[ rows ][ columns ]  Data to copy.
         Integer     firstRow     Index of first row of subset.
         Integer     lastRow      Index of alst row os subset.
         Integer     firstColumn  Index of first column of subset.
         Integer     lastColumn   Index of last column of subset.
         Real        output[ points ] Array of subset data to append data to.
******************************************************************************/

static void copySubsetData( const Integer rows, const Integer columns,
                            const Real data[],
                            Integer firstRow, Integer lastRow,
                            Integer firstColumn, Integer lastColumn,
                            Real output[] ) {

  PRE010( rows > 0,
          columns > 0,
          rows * columns > 0,
          data,
          isNanFree( data, rows * columns ),
          IN_RANGE( firstRow, 0, rows - 1 ),
          IN_RANGE( lastRow, firstRow, rows - 1 ),
          IN_RANGE( firstColumn, 0, columns - 1 ),
          IN_RANGE( lastColumn, firstColumn, columns - 1 ),
          output );

  CHECKING( Integer points = 0; )
  Real* values = output;
  Integer row = 0;

  for ( row = firstRow; row <= lastRow; ++row ) {
    const Integer rowOffset = row * columns;
    Integer column = 0;
    CHECK( rowOffset + lastColumn < rows * columns );

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Real value = data[ index ];

      if ( value >= 0.0 ) {
        CHECKING(  ++points; )
        CHECK( IN_RANGE( points, 1,
                        ( 1 + lastRow - firstRow ) *
                        ( 1 + lastColumn - firstColumn ) ) );
        *values++ = value;
      }
    }
  }

  POST0( isNanFree( output, points ) );
}



/******************************************************************************
PURPOSE: computeMean - Compute running mean of subset of scan points.
INPUTS:  const Data* const data  Data structure with scan, subset, etc.
         Integer counts[ rows ][ columns ]    Current mean counts.
         Real    values[ rows ][ columns ]    Current means.
OUTPUTS: Integer counts[ rows ][ columns ]    Updated mean counts.
         Real    means[  rows ][ columns ]    Updated means.
******************************************************************************/

static void computeMean( const Data* const data,
                         Integer counts[], Real means[] ) {

  PRE014( data,
          data->rows > 0,
          data->columns > 0,
          data->rows * data->columns > 0,
          data->scan.data,
          isValidScan( &data->scan ),
          IN_RANGE( data->subset[ ROW ][ MINIMUM ], 0, data->rows - 1 ),
          IN_RANGE( data->subset[ ROW ][ MAXIMUM ],
                    data->subset[ ROW ][ MINIMUM ], data->rows - 1 ),
          IN_RANGE( data->subset[ COLUMN][MINIMUM], 0, data->columns - 1),
          IN_RANGE( data->subset[ COLUMN ][ MAXIMUM ],
                    data->subset[ COLUMN ][ MINIMUM ], data->columns - 1),
          counts,
          means,
          minimumItemI( counts, data->rows * data->columns ) >= 0,
          minimumItem(  means,  data->rows * data->columns ) >= 0.0 );

  const Real* const scanData = data->scan.data;
  const Integer columns      = data->columns;
  const Integer firstRow     = data->subset[ ROW    ][ MINIMUM ];
  const Integer lastRow      = data->subset[ ROW    ][ MAXIMUM ];
  const Integer firstColumn  = data->subset[ COLUMN ][ MINIMUM ];
  const Integer lastColumn   = data->subset[ COLUMN ][ MAXIMUM ];
  Integer row = 0;

#pragma omp parallel for

  for ( row = firstRow; row <= lastRow; ++row ) {
    const Integer rowOffset = row * columns;
    Integer column = 0;
    CHECK( rowOffset + lastColumn < data->rows * data->columns );

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Real value = scanData[ index ];

      if ( value >= 0.0 ) { /* All valid GOES data is non-negative. */
        const Real mean = means[ index ];
        const Integer count = counts[ index ];
        means[ index ] = ( count * mean + value ) / ( count + 1 );
        counts[ index ] += 1;
      }
    }
  }

  POST02( minimumItemI( counts, data->rows * data->columns ) >= 0,
          minimumItem(  means,  data->rows * data->columns ) >= 0.0 );
}



/******************************************************************************
PURPOSE: meanPoints - Count number of non-zero counts.
INPUTS:  const Data* data                    Data structure with subset.
         Integer counts[ rows ][ columns ]   Final mean counts.
RETURNS: Integer number of non-zero counts.
******************************************************************************/

static Integer meanPoints( const Data* data, const Integer counts[] ) {

  PRE010( data,
          data->ok,
          data->rows > 0,
          data->columns > 0,
          data->rows * data->columns > 0,
          IN_RANGE( data->subset[ ROW ][ MINIMUM ], 0, data->rows - 1 ),
          IN_RANGE( data->subset[ ROW ][ MAXIMUM ],
                    data->subset[ ROW ][ MINIMUM ], data->rows - 1 ),
          IN_RANGE( data->subset[ COLUMN ][ MINIMUM ], 0, data->columns - 1 ),
          IN_RANGE( data->subset[ COLUMN ][ MAXIMUM ],
                    data->subset[ COLUMN ][ MINIMUM ], data->columns - 1 ),
          counts );

  const Integer columns = data->columns;
  const Integer firstRow    = data->subset[ ROW    ][ MINIMUM ];
  const Integer lastRow     = data->subset[ ROW    ][ MAXIMUM ];
  const Integer firstColumn = data->subset[ COLUMN ][ MINIMUM ];
  const Integer lastColumn  = data->subset[ COLUMN ][ MAXIMUM ];
  Integer row = 0;
  Integer result = 0;

  for ( row = firstRow; row <= lastRow; ++row ) {
    const Integer rowOffset = row * columns;
    Integer column = 0;
    CHECK( rowOffset + lastColumn < data->rows * data->columns );

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Integer count = counts[ index ];
      result += count != 0;
    }
  }

  POST0( IN_RANGE( result, 0,
                   ( 1 + data->subset[ ROW    ][ MAXIMUM ] -
                         data->subset[ ROW    ][ MINIMUM ] ) *
                   ( 1 + data->subset[ COLUMN ][ MAXIMUM ] -
                         data->subset[ COLUMN ][ MINIMUM ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: copyMeanData - Copy subset of mean data.
INPUTS:  const Integer rows                 Number of data rows.
         const Integer columns              Number of data columns.
         const Integer subset[ 2 ][ 2 ]     subset[ROW COLUMN][MIN/MAXIMUM].
         const Integer points               Subset points.
         const Integer counts[ rows ][ columns ]  Final counts.
         const Real    means[  rows ][ columns ]  Final means.
         const Real longitudes[  rows ][ columns ]  Longitudes.
         const Real latitudes[   rows ][ columns ]  Latitudes.
         const Real corners[ 2 ][ 4 ][ rows ][ columns ]    Lon-lat corners.
OUTPUTS: Real subsetLongitudes[ points ]    Subset longitudes.
         Real subsetLatitudes[  points ]    Subset longitudes.
         Real subsetData[       points ]    Subset data, e.g., AOD.
         Real subsetLongitudesSW[ points ]  Optional: for corners or 0.
         Real subsetLongitudesSE[ points ]  Optional: for corners or 0.
         Real subsetLongitudesNW[ points ]  Optional: for corners or 0.
         Real subsetLongitudesNE[ points ]  Optional: for corners or 0.
         Real subsetLatitudesSW[ points ]   Optional: for corners or 0.
         Real subsetLatitudesSE[ points ]   Optional: for corners or 0.
         Real subsetLatitudesNW[ points ]   Optional: for corners or 0.
         Real subsetLatitudesNE[ points ]   Optional: for corners or 0.
******************************************************************************/

static void copyMeanData( const Integer rows,
                          const Integer columns,
                          const Integer subset[ 2 ][ 2 ],
                          const Integer points,
                          const Integer counts[],
                          const Real means[],
                          const Real longitudes[],
                          const Real latitudes[],
                          const Real corners[],
                          Real subsetLongitudes[],
                          Real subsetLatitudes[],
                          Real subsetData[],
                          Real subsetLongitudesSW[],
                          Real subsetLongitudesSE[],
                          Real subsetLongitudesNW[],
                          Real subsetLongitudesNE[],
                          Real subsetLatitudesSW[],
                          Real subsetLatitudesSE[],
                          Real subsetLatitudesNW[],
                          Real subsetLatitudesNE[] ) {

  PRE016( rows > 0,
          columns > 0,
          rows * columns > 0,
          subset,
          counts,
          means,
          IN_RANGE( subset[ ROW ][ MINIMUM ], 0, rows - 1 ),
          IN_RANGE( subset[ ROW ][ MAXIMUM ],
                    subset[ ROW ][ MINIMUM ], rows - 1 ),
          IN_RANGE( subset[ COLUMN ][ MINIMUM ], 0, columns - 1 ),
          IN_RANGE( subset[ COLUMN ][ MAXIMUM ],
                    subset[ COLUMN ][ MINIMUM ], columns - 1 ),
          IN_RANGE( points, 1, rows * columns ),
          IN_RANGE( points, 1,
                    ( 1 + subset[ ROW ][ MAXIMUM ] - subset[ROW][MINIMUM] ) *
                    ( 1 + subset[COLUMN][MAXIMUM ] - subset[COLUMN][MINIMUM])),
          subsetLongitudes,
          subsetLatitudes,
          subsetData,
          IMPLIES_ELSE( subsetLongitudesSW,
                        NON_ZERO7( subsetLongitudesSE,
                                   subsetLongitudesNW,
                                   subsetLongitudesNE,
                                   subsetLatitudesSW,
                                   subsetLatitudesSE,
                                   subsetLatitudesNW,
                                   subsetLatitudesNE ),
                        IS_ZERO7( subsetLongitudesSE,
                                  subsetLongitudesNW,
                                  subsetLongitudesNE,
                                  subsetLatitudesSW,
                                  subsetLatitudesSE,
                                  subsetLatitudesNW,
                                  subsetLatitudesNE ) ) );

  const Integer firstRow    = subset[ ROW    ][ MINIMUM ];
  const Integer lastRow     = subset[ ROW    ][ MAXIMUM ];
  const Integer firstColumn = subset[ COLUMN ][ MINIMUM ];
  const Integer lastColumn  = subset[ COLUMN ][ MAXIMUM ];
  const Integer rowsTimesColumns = rows * columns;
  const Real* const longitudesSW = corners;
  const Real* const longitudesSE = corners + rowsTimesColumns;
  const Real* const longitudesNW = corners + rowsTimesColumns * 2;
  const Real* const longitudesNE = corners + rowsTimesColumns * 3;
  const Real* const latitudesSW  = corners + rowsTimesColumns * 4;
  const Real* const latitudesSE  = corners + rowsTimesColumns * 5;
  const Real* const latitudesNW  = corners + rowsTimesColumns * 6;
  const Real* const latitudesNE  = corners + rowsTimesColumns * 7;
  Integer row = 0;
  Integer output = 0;

  for ( row = firstRow, output = 0; row <= lastRow; ++row ) {
    const Integer rowOffset = row * columns;
    Integer column = 0;
    CHECK( rowOffset + lastColumn < rows * columns );

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Integer index = rowOffset + column;
      const Integer count = counts[ index ];
      CHECK( IN_RANGE( output, 0, points - 1 ) );

      if ( count ) {
        const Real longitude = longitudes[ index ];
        const Real latitude  = latitudes[ index ];
        const Real mean      = means[ index ];
        CHECK3( isValidLongitude( longitude ),
                isValidLatitude( latitude ),
                mean >= 0.0 );
        subsetLongitudes[ output ] = longitude;
        subsetLatitudes[  output ] = latitude;
        subsetData[       output ] = mean;

        if ( subsetLongitudesSW ) {
          const Real longitudeSW = longitudesSW[ index ];
          const Real longitudeSE = longitudesSE[ index ];
          const Real longitudeNW = longitudesNW[ index ];
          const Real longitudeNE = longitudesNE[ index ];
          const Real latitudeSW  = latitudesSW[  index ];
          const Real latitudeSE  = latitudesSE[  index ];
          const Real latitudeNW  = latitudesNW[  index ];
          const Real latitudeNE  = latitudesNE[  index ];
          CHECK8( isValidLongitude( longitudeSW ),
                  isValidLongitude( longitudeSE ),
                  isValidLongitude( longitudeNE ),
                  isValidLongitude( longitudeNW ),
                  isValidLatitude( latitudeSW ),
                  isValidLatitude( latitudeSE ),
                  isValidLatitude( latitudeNE ),
                  isValidLatitude( latitudeNW ) );
          subsetLongitudesSW[ output ] = longitudeSW;
          subsetLongitudesSE[ output ] = longitudeSE;
          subsetLongitudesNW[ output ] = longitudeNW;
          subsetLongitudesNE[ output ] = longitudeNE;
          subsetLatitudesSW[  output ] = latitudeSW;
          subsetLatitudesSE[  output ] = latitudeSE;
          subsetLatitudesNW[  output ] = latitudeNW;
          subsetLatitudesNE[  output ] = latitudeNE;
        }

        ++output;
      }
    }
  }

  CHECK( output == points );
  POST03( validLongitudesAndLatitudes( points,
                                       subsetLongitudes, subsetLatitudes ),
          isNanFree( subsetData, points ),
          IMPLIES( subsetLongitudesSW,
                   AND4( validLongitudesAndLatitudes( points,
                                                      subsetLongitudesSW,
                                                      subsetLatitudesSW ),
                         validLongitudesAndLatitudes( points,
                                                      subsetLongitudesSE,
                                                      subsetLatitudesSE ),
                         validLongitudesAndLatitudes( points,
                                                      subsetLongitudesNW,
                                                      subsetLatitudesNW ),
                         validLongitudesAndLatitudes( points,
                                                      subsetLongitudesNE,
                                                      subsetLatitudesNE ) ) ));
}



/******************************************************************************
PURPOSE: writeData - Write subsetted scan data to stdout.
INPUTS:  Data* data  Data to write.
OUTPTUS: Data* data  data->ok.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void writeData( Data* data ) {
  PRE0( isValidData( data ) );
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    writeHeader( data, output );

    if ( data->ok ) {
      writeXDR( data, output );
    }

    FREE_OBJECT( output );
  }

  POST0( isValidData( data ) );
}



/******************************************************************************
PURPOSE: writeHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data      Data to write to stream.
         Stream* output  Stream to write data to.
OUTPUTS: Data* data      data->ok set to indicate success or failure.
         Stream* output  Stream to written to.
******************************************************************************/

static void writeHeader( Data* data, Stream* output ) {

  PRE06( isValidData( data ), data->ok,
         output, output->invariant( output ), output->ok( output ),
         output->isWritable( output ) );

  const Arguments* const arguments = &data->arguments;
  const VoidList*  const scans = data->subsettedScans;
  const SubsettedScan* const firstScan = scans->item( scans, 0 );
  const Integer variables = firstScan->variables;
  const Integer writeCorners = arguments->corners;
  const char* const daily = arguments->daily ? "daily_" : "";
  UTCTimestamp timestamp;
  toUTCTimestamp( arguments->firstTimestamp, timestamp );
  CHECK( variables > 0 );

  output->writeString( output,
                       "Swath 2.0\n%s\n%s\n"
                       "# Dimensions: variables timesteps scans:\n"
                       "%"INTEGER_FORMAT" %"INTEGER_FORMAT" %"
                       INTEGER_FORMAT"\n"
                       "# Variable names:\nLongitude Latitude %s%s",
                       arguments->description, timestamp,
                       variables,
                       arguments->hours, scans->count( scans ),
                       daily, data->variable );

  if ( writeCorners ) {
    output->writeString( output,
                         " Longitude_SW Longitude_SE"
                         " Longitude_NW Longitude_NE"
                         " Latitude_SW Latitude_SE"
                         " Latitude_NW Latitude_NE" );
  }

  if ( output->ok( output ) ) {
    output->writeString( output, "\n# Variable units:\ndeg deg %s",
                         arguments->units );
  }

  if ( writeCorners ) {
    output->writeString( output, " deg deg deg deg deg deg deg deg" );
  }

  if ( output->ok( output ) ) {
    output->writeString( output,
                         "\n# Domain: <min_lon> <min_lat> <max_lon> <max_lat>"
                         "\n%g %g %g %g\n"
                         "# MSB 64-bit integers (yyyydddhhmm)"
                         " timestamps[scans] and\n"
                         "# MSB 64-bit integers points[scans] and\n"
                         "# IEEE-754 64-bit reals"
                         " data_1[variables][points_1] ..."
                         " data_S[variables][points_S]:\n",
                         arguments->bounds[ LONGITUDE ][ MINIMUM ],
                         arguments->bounds[ LATITUDE  ][ MINIMUM ],
                         arguments->bounds[ LONGITUDE ][ MAXIMUM ],
                         arguments->bounds[ LATITUDE  ][ MAXIMUM ] );
  }

  data->ok = output->ok( output );

  POST04( isValidData( data ), output->invariant( output ),
          output->isWritable( output ), data->ok == output->ok( output ) );
}



/******************************************************************************
PURPOSE: writeXDR - Write XDR format data arrays of subset to stdout.
INPUTS:  Data* data      Data to write to stream.
         Stream* output  Stream to write data to.
OUTPUTS: Data* data      data->ok set to indicate success or failure.
         Stream* output  Stream to written to.
******************************************************************************/

static void writeXDR( Data* data, Stream* output ) {

  PRE06( isValidData( data ), data->ok, output, output->invariant( output ),
         output->isWritable( output ), output->ok( output ) );

  const VoidList* scans = data->subsettedScans;
  writeScanTimestamps( scans, output );
  data->ok = output->ok( output );

  if ( data->ok ) {
    writeScanPoints( scans, output );
    data->ok = output->ok( output );

    if ( data->ok ) {
      writeScanData( scans, output );
      data->ok = output->ok( output );
    }
  }

  POST04( isValidData( data ), output->invariant( output ),
          output->isWritable( output ), data->ok == output->ok( output ) );
}



/******************************************************************************
PURPOSE: writeScanTimestamps - Write MSB 64-bit integer scan timestamps.
INPUTS:  const VoidList* scans   List of SubsettedScan*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeScanTimestamps( const VoidList* scans, Stream* output ) {

  PRE07( scans, scans->invariant( scans ), scans->count( scans ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer scanCount = scans->count( scans );
  Integer scanIndex = 0;

  do {
    const SubsettedScan* const scan = scans->item( scans, scanIndex );
    CHECK( isValidSubsettedScan( scan ) );
    output->write64BitInteger( output, scan->timestamp );
    ++scanIndex;
  } while ( AND2( output->ok( output ), scanIndex < scanCount ) );

  POST04( scans->invariant( scans ), scans->count( scans ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: writeScanPoints - Write MSB 64-bit integer scan subset point counts.
INPUTS:  const VoidList* scans   List of SubsettedScan*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeScanPoints( const VoidList* scans, Stream* output ) {

  PRE07( scans, scans->invariant( scans ), scans->count( scans ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer scanCount = scans->count( scans );
  Integer scanIndex = 0;

  do {
    const SubsettedScan* const scan = scans->item( scans, scanIndex );
    CHECK( isValidSubsettedScan( scan ) );
    output->write64BitInteger( output, scan->points );
    ++scanIndex;
  } while ( AND2( output->ok( output ), scanIndex < scanCount ) );

  POST04( scans->invariant( scans ), scans->count( scans ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: writeScanData - Write 64-bit IEEE-754 scan subset variable data.
INPUTS:  const VoidList* scans   List of SubsettedScan*.
         Stream*         output  Stream to write to.
OUTPUTS: Stream*         output->ok( output ) indicates success or failure.
******************************************************************************/

static void writeScanData( const VoidList* scans, Stream* output ) {

  PRE07( scans, scans->invariant( scans ), scans->count( scans ) > 0,
         output, output->invariant( output ), output->isWritable( output ),
         output->ok( output ) );

  const Integer scanCount = scans->count( scans );
  Integer scanIndex = 0;

  do {
    const SubsettedScan* const scan = scans->item( scans, scanIndex );
    CHECK( isValidSubsettedScan( scan ) );
    output->write64BitReals( output, scan->data,
                             scan->variables * scan->points );
    ++scanIndex;
  } while ( AND2( output->ok( output ), scanIndex < scanCount ) );

  POST04( scans->invariant( scans ), scans->count( scans ) > 0,
          output->invariant( output ), output->isWritable( output ) );
}



/******************************************************************************
PURPOSE: computeCorners - Compute and store the corner variables
         (longitude_sw, ... , latitude_ne) for each center/pixel.
INPUTS:  const Integer rows                       Number of rows of cells.
         const Integer columns                    Number of columns of cells.
         const Real longitudes[ rows * columns ]  Longitudes of cell centers.
         const Real latitudes[  rows * columns ]  Latitudes  of cell centers.
OUTPUTS: Real corners[ 8 * rows * columns ]       longitudes_sw
                                                  longitudes_se
                                                  longitudes_nw
                                                  longitudes_ne
                                                  latitudes_sw
                                                  latitudes_se
                                                  latitudes_nw
                                                  latitudes_ne
NOTES:   Uses linear interpolation and extrapolation to the edges.
******************************************************************************/

static void computeCorners( const Integer rows, const Integer columns,
                            const Real longitudes[], const Real latitudes[],
                            Real corners[] ) {

  PRE07( rows > 0, columns > 0, rows * columns > 0, longitudes, latitudes,
         corners,
         validLongitudesAndLatitudes( rows * columns, longitudes, latitudes ));

  const Integer rows_1      = rows - 1;
  const Integer columns_1   = columns - 1;
  const Integer cells       = rows * columns;
  Real* const longitudes_sw = corners;
  Real* const longitudes_se = longitudes_sw + cells;
  Real* const longitudes_nw = longitudes_se + cells;
  Real* const longitudes_ne = longitudes_nw + cells;
  Real* const latitudes_sw  = longitudes_ne + cells;
  Real* const latitudes_se  = latitudes_sw + cells;
  Real* const latitudes_nw  = latitudes_se + cells;
  Real* const latitudes_ne  = latitudes_nw + cells;
  Integer cell = 0;
  Integer index = 0;

  if ( OR2( rows < 2, columns < 2 ) ) {

    /* Copy all center values to the corners in such degenerate cases: */

#pragma omp parallel for

    for ( cell = 0; cell < cells; ++cell ) {
      longitudes_sw[ cell ] =
      longitudes_se[ cell ] =
      longitudes_nw[ cell ] =
      longitudes_ne[ cell ] = longitudes[ cell ];
      latitudes_sw[ cell ] =
      latitudes_se[ cell ] =
      latitudes_nw[ cell ] =
      latitudes_ne[ cell ] = latitudes[ cell ];
    }

  } else { /* Linearly interpolate and extrapolate the corner points: */
    Integer row    = 0;
    Integer column = 0;

    /*
     * First compute linearly interpolated corners of all interior cells:
     * Note: rows increase north to south and columns increase west to east.
     */

#pragma omp parallel for private( column )

    for ( row = 0; row < rows_1; ++row ) {
      const Integer rowOffset     = row * columns;
      const Integer nextRowOffset = rowOffset + columns;

      /* Interior row, interior columns: */

      for ( column = 0; column < columns_1; ++column ) {
        const Integer thisIndex         = rowOffset + column;
        const Integer nextColumn        = thisIndex + 1;
        const Integer nextRow           = nextRowOffset + column;
        const Integer nextRowNextColumn = nextRow + 1;

        const Real longitude                  = longitudes[ thisIndex ];
        const Real nextColumnLongitude        = longitudes[ nextColumn ];
        const Real nextRowLongitude           = longitudes[ nextRow ];
        const Real nextRowNextColumnLongitude = longitudes[ nextRowNextColumn];

        const Real latitude                  = latitudes[ thisIndex ];
        const Real nextColumnLatitude        = latitudes[ nextColumn ];
        const Real nextRowLatitude           = latitudes[ nextRow ];
        const Real nextRowNextColumnLatitude = latitudes[ nextRowNextColumn ];

        const Real interpolatedLongitude = 0.25 *
          ( longitude + nextColumnLongitude + nextRowLongitude +
            nextRowNextColumnLongitude );

        const Real interpolatedLatitude = 0.25 *
          ( latitude + nextColumnLatitude + nextRowLatitude +
            nextRowNextColumnLatitude );

        longitudes_ne[ thisIndex         ] = interpolatedLongitude;
        longitudes_nw[ nextColumn        ] = interpolatedLongitude;
        longitudes_se[ nextRow           ] = interpolatedLongitude;
        longitudes_sw[ nextRowNextColumn ] = interpolatedLongitude;

        latitudes_ne[ thisIndex         ] = interpolatedLatitude;
        latitudes_nw[ nextColumn        ] = interpolatedLatitude;
        latitudes_se[ nextRow           ] = interpolatedLatitude;
        latitudes_sw[ nextRowNextColumn ] = interpolatedLatitude;
      } /* End loop on interior columns. */

    } /* End parallel loop on interior rows. */

    /* Serial region (not worth parallelizing): */

    /* Last row, interior columns (extrapolated top edge): */

    for ( column = 1, index = rows_1 * columns + 1; column < columns;
          ++column, ++index ) {
      const Integer previousColumn = index - 1;

      const Real longitude = longitudes[ index ];
      const Real previousColumnLongitude = longitudes[ previousColumn ];
      const Real midpointLongitude =
        0.5 * ( longitude + previousColumnLongitude );
      const Real interpolatedLongitude = longitudes_sw[ index ];
      const Real longitudeDifference =
        midpointLongitude - interpolatedLongitude;
      const Real extrapolatedLongitude =
        midpointLongitude + longitudeDifference;

      const Real latitude = latitudes[ index ];
      const Real previousColumnLatitude = latitudes[ previousColumn ];
      const Real midpointLatitude = 0.5 * ( latitude + previousColumnLatitude);
      const Real interpolatedLatitude = latitudes_sw[ index ];
      const Real latitudeDifference = midpointLatitude - interpolatedLatitude;
      const Real extrapolatedLatitude = midpointLatitude + latitudeDifference;

      longitudes_nw[ index          ] = extrapolatedLongitude;
      longitudes_ne[ previousColumn ] = extrapolatedLongitude;

      latitudes_nw[ index          ] = extrapolatedLatitude;
      latitudes_ne[ previousColumn ] = extrapolatedLatitude;
    }

    /* First row, interior columns (extrapolated bottom edge): */

    for ( column = 1, index = 1; column < columns; ++column, ++index ) {
      const Integer previousColumn = index - 1;

      const Real longitude = longitudes[ index ];
      const Real previousColumnLongitude = longitudes[ previousColumn ];
      const Real midpointLongitude =
        0.5 * ( longitude + previousColumnLongitude );
      const Real interpolatedLongitude = longitudes_nw[ index ];
      const Real longitudeDifference =
        midpointLongitude - interpolatedLongitude;
      const Real extrapolatedLongitude =
        midpointLongitude + longitudeDifference;

      const Real latitude = latitudes[ index ];
      const Real previousColumnLatitude = latitudes[ previousColumn ];
      const Real midpointLatitude = 0.5 * ( latitude + previousColumnLatitude);
      const Real interpolatedLatitude = latitudes_nw[ index ];
      const Real latitudeDifference = midpointLatitude - interpolatedLatitude;
      const Real extrapolatedLatitude = midpointLatitude + latitudeDifference;

      longitudes_sw[ index          ] = extrapolatedLongitude;
      longitudes_se[ previousColumn ] = extrapolatedLongitude;

      latitudes_sw[ index          ] = extrapolatedLatitude;
      latitudes_se[ previousColumn ] = extrapolatedLatitude;
    }

    /* First column, interior rows (extrapolated left edge, except corners): */

    for ( row = 1, index = columns; row < rows; ++row, index += columns ) {
      const Integer previousRow = index - columns;

      const Real longitude = longitudes[ index ];
      const Real previousRowLongitude = longitudes[ previousRow ];
      const Real midpointLongitude =
        0.5 * ( longitude + previousRowLongitude );
      const Real interpolatedLongitude = longitudes_se[ index ];
      const Real longitudeDifference =
        midpointLongitude - interpolatedLongitude;
      const Real extrapolatedLongitude =
        midpointLongitude + longitudeDifference;

      const Real latitude = latitudes[ index ];
      const Real previousRowLatitude = latitudes[ previousRow ];
      const Real midpointLatitude = 0.5 * ( latitude + previousRowLatitude );
      const Real interpolatedLatitude = latitudes_se[ index ];
      const Real latitudeDifference = midpointLatitude - interpolatedLatitude;
      const Real extrapolatedLatitude = midpointLatitude + latitudeDifference;

      longitudes_sw[ index       ] = extrapolatedLongitude;
      longitudes_nw[ previousRow ] = extrapolatedLongitude;

      latitudes_sw[ index       ] = extrapolatedLatitude;
      latitudes_nw[ previousRow ] = extrapolatedLatitude;
    }

    /* Last column, interior rows (extrapolated right edge, except corners): */

    for ( row = 1, index = columns + columns - 1;
          row < rows; ++row, index += columns ) {
      const Integer previousRow = index - columns;

      const Real longitude = longitudes[ index ];
      const Real previousRowLongitude = longitudes[ previousRow ];
      const Real midpointLongitude =
        0.5 * ( longitude + previousRowLongitude );
      const Real interpolatedLongitude = longitudes_sw[ index ];
      const Real longitudeDifference =
        midpointLongitude - interpolatedLongitude;
      const Real extrapolatedLongitude =
        midpointLongitude + longitudeDifference;

      const Real latitude = latitudes[ index ];
      const Real previousRowLatitude = latitudes[ previousRow ];
      const Real midpointLatitude = 0.5 * ( latitude + previousRowLatitude );
      const Real interpolatedLatitude = latitudes_sw[ index ];
      const Real latitudeDifference = midpointLatitude - interpolatedLatitude;
      const Real extrapolatedLatitude = midpointLatitude + latitudeDifference;

      longitudes_se[ index       ] = extrapolatedLongitude;
      longitudes_ne[ previousRow ] = extrapolatedLongitude;

      latitudes_se[ index       ] = extrapolatedLatitude;
      latitudes_ne[ previousRow ] = extrapolatedLatitude;
    }

    /* First row, first column cell (extrapolated bottom-left corner): */

    {
      const Real longitude             = longitudes[ 0 ];
      const Real latitude              = latitudes[ 0 ];
      const Real diagonalLongitude     = longitudes_ne[ 0 ];
      const Real diagonalLatitude      = latitudes_ne[ 0 ];
      const Real longitudeDifference   = longitude - diagonalLongitude;
      const Real latitudeDifference    = latitude  - diagonalLatitude;
      const Real extrapolatedLongitude = longitude + longitudeDifference;
      const Real extrapolatedLatitude  = latitude  + latitudeDifference;
      longitudes_sw[ 0 ]                = extrapolatedLongitude;
      latitudes_sw[  0 ]                = extrapolatedLatitude;
    }

    /* First row, last column cell (extrapolated bottom-right corner): */

    {
      const Real longitude             = longitudes[ columns_1 ];
      const Real latitude              = latitudes[ columns_1 ];
      const Real diagonalLongitude     = longitudes_nw[ columns_1 ];
      const Real diagonalLatitude      = latitudes_nw[ columns_1 ];
      const Real longitudeDifference   = longitude - diagonalLongitude;
      const Real latitudeDifference    = latitude  - diagonalLatitude;
      const Real extrapolatedLongitude = longitude + longitudeDifference;
      const Real extrapolatedLatitude  = latitude  + latitudeDifference;
      longitudes_se[ columns_1 ]        = extrapolatedLongitude;
      latitudes_se[  columns_1 ]        = extrapolatedLatitude;
    }

    /* Last row, first column cell (extrapolated top-left corner): */

    index = cells - columns;

    {
      const Real longitude             = longitudes[ index ];
      const Real latitude              = latitudes[ index ];
      const Real diagonalLongitude     = longitudes_se[ index ];
      const Real diagonalLatitude      = latitudes_se[ index ];
      const Real longitudeDifference   = longitude - diagonalLongitude;
      const Real latitudeDifference    = latitude  - diagonalLatitude;
      const Real extrapolatedLongitude = longitude + longitudeDifference;
      const Real extrapolatedLatitude  = latitude  + latitudeDifference;
      longitudes_nw[ index ]            = extrapolatedLongitude;
      latitudes_nw[  index ]            = extrapolatedLatitude;
    }

    /* Last row, last column cell (extrapolated top-right corner): */

    index = cells - 1;

    {
      const Real longitude             = longitudes[ index ];
      const Real latitude              = latitudes[ index ];
      const Real diagonalLongitude     = longitudes_sw[ index ];
      const Real diagonalLatitude      = latitudes_sw[ index ];
      const Real longitudeDifference   = longitude - diagonalLongitude;
      const Real latitudeDifference    = latitude  - diagonalLatitude;
      const Real extrapolatedLongitude = longitude + longitudeDifference;
      const Real extrapolatedLatitude  = latitude  + latitudeDifference;
      longitudes_ne[ index ]            = extrapolatedLongitude;
      latitudes_ne[  index ]            = extrapolatedLatitude;
    }

    /* Clamp any out-of-range values: */

#pragma omp parallel for

    for ( cell = 0; cell < cells; ++cell ) {
      longitudes_nw[cell] = CLAMPED_TO_RANGE(longitudes_nw[cell],-180.0,180.0);
      longitudes_sw[cell] = CLAMPED_TO_RANGE(longitudes_sw[cell],-180.0,180.0);
      longitudes_se[cell] = CLAMPED_TO_RANGE(longitudes_se[cell],-180.0,180.0);
      longitudes_ne[cell] = CLAMPED_TO_RANGE(longitudes_ne[cell],-180.0,180.0);
      latitudes_nw[cell] = CLAMPED_TO_RANGE( latitudes_nw[cell], -90.0, 90.0 );
      latitudes_sw[cell] = CLAMPED_TO_RANGE( latitudes_sw[cell], -90.0, 90.0 );
      latitudes_se[cell] = CLAMPED_TO_RANGE( latitudes_se[cell], -90.0, 90.0 );
      latitudes_ne[cell] = CLAMPED_TO_RANGE( latitudes_ne[cell], -90.0, 90.0 );
    }

  } /* End else non-degenerate cases: */

  POST0( validLongitudesAndLatitudes( 4 * rows * columns,
                                      corners, corners + 4 * rows * columns ));
}


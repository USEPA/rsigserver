
/******************************************************************************
PURPOSE: OMISubset.c - Read a set of a OMI-AURA files, subset the swath to a
         domain (longitude-latitude rectangle) and write it to
         stdout as XDR (IEEE-754) format binary.

NOTES:   Uses HDF5 Library and libUtilities.a (../../libs/Utilities).

HISTORY: 2016-09-12 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h>    /* For macro assert(). */
#include <stdio.h>     /* For FILE, stderr, fprintf(), fopen(), fwrite(). */
#include <string.h>    /* For strchr(). */
#include <unistd.h>    /* For unlink(), getpid() */

#include <Utilities.h> /* For PRE0*(), NEW_ZERO(), Stream, VoidList, etc. */
#include "ReadFile.h"  /* For openFile(), swathsInFile(), readDataset(). */

/*================================== TYPES ==================================*/

enum { NAME_LENGTH = 80 };
typedef char Name[ NAME_LENGTH + 1 ]; /* Variable name. */

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;       /* File containing list of OMI files to read*/
  const char* tmpdir;         /* Name of directory to write temp files. */
  const char* description;    /* User-supplied description. */
  const char* variableName;   /* Selected variable name. */
  Integer     firstTimestamp; /* YYYYDDDHH00 of subset. */
  Integer     timesteps;      /* Number of hours in subset. */
  Integer     corners;        /* Compute interpolated lon-lat corner points? */
  Bounds      domain;         /* domain[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM] */
  double      maximumCloudFraction; /* Maximum acceptable cloud fraction. */
  double      maximumSolarZenithAngle; /* Max acceptable solar zenith angle. */
  int         allowNegativeCounts;     /* Allow negative molecules/cm2? */
} Arguments;

/* Swath info: */

typedef struct {
  Integer timestamp; /* YYYYDDDHHMM of file containing swath. */
  Integer points;    /* Number of points in domain-subsetted filtered swath.*/
} SwathInfo;

/* Data type: */

typedef struct {
  Arguments arguments;    /* User-supplied (command-line) arguments. */
  FileName  tempFileName; /* Name of temp file of output subset data. */
  FILE*     tempFile;     /* Temp file of output subset data. */
  VoidList* swaths;       /* swaths[ swath ]->SwathInfo. */
  Name      units;        /* Units of arguments->variableName. */
  size_t    rows;         /* Rows of unsubsetted data. */
  size_t    columns;      /* Columns of unsubsetted data. */
  double*   buffer;       /* Allocated buffer[ variables * rows * columns ]; */
  Integer   ok;           /* Did last command succeed? */
} Data;

/*================================= CONSTANTS ================================*/

/* Name of temporary file created in -tmpdir will have PID appended: */

#define TEMP_FILE_NAME "junk_OMISubset"

static const double EDGE = 179.99; /* Clamp longitude to EDGE if cross +/- 180*/

/*========================== FORWARD DECLARATIONS ===========================*/

/* Commands: */

CHECKING( static Integer isValidArguments( const Arguments* arguments ); )

CHECKING( static Integer isValidSwathInfo( const SwathInfo* swathInfo ); )

CHECKING( static Integer isValidData( const Data* data ); )

static void deallocateData( Data* data );

static void deallocateSwathInfo( void* argument );

static void printUsage( const char* programName );

static Integer parseArguments( Integer argc, char* argv[],
                               Arguments* arguments );

static Integer parseDomain( Integer argc, char* argv[], Arguments* arguments );

static Integer parseCorners( Integer argc, char* argv[], Arguments* arguments);

static Integer parseMaximumCloudFraction( Integer argc, char* argv[],
                                          Arguments* arguments );

static Integer parseMaximumSolarZenithAngle( Integer argc, char* argv[],
                                             Arguments* arguments );

static Integer parseAllowNegativeCounts( Integer argc, char* argv[],
                                         Arguments* arguments );

static Integer parseFileTimestamp( const char* fileName );

static void readData( Data* data );

static SwathInfo* readSwathData( Data* data,
                                 const int file,
                                 const char* const product,
                                 const Integer timestamp );

static int clampInvalidCoordinates( const size_t points,
                                    double longitudes[],
                                    double latitudes[] );

static size_t pointsInSubset( const Bounds domain,
                              const size_t points,
                              const double longitudes[],
                              const double latitudes[],
                              double values[],
                              const double longitudesSW[],
                              const double longitudesSE[],
                              const double longitudesNW[],
                              const double longitudesNE[],
                              const double latitudesSW[],
                              const double latitudesSE[],
                              const double latitudesNW[],
                              const double latitudesNE[] );

static void compactSubsetData( const size_t subsetPoints,
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

static void computeCorners( const size_t rows,
                            const size_t columns,
                            const double longitudes[],
                            const double latitudes[],
                            double longitudesSW[],
                            double longitudesSE[],
                            double longitudesNW[],
                            double longitudesNE[],
                            double latitudesSW[],
                            double latitudesSE[],
                            double latitudesNW[],
                            double latitudesNE[] );

static void clampLongitudes( const double longitude,
                             double* longitude1,
                             double* longitude2,
                             double* longitude3,
                             double* longitude4 );

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

static void streamData( Data* data );

static void streamHeader( const Data* const data );

static int streamSwathTimestamps( const VoidList* swaths );

static int streamSwathPoints( const VoidList* swaths );

static void streamTempFile( Data* data );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Read a subset of a OMI file and write it to stdout in
         XDR or ASCII format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Integer ok = 0;

  if ( ! isValidArgs( argc, (const char**) argv ) ) {
    failureMessage( "Invalid command-line arguments." );
  } else {
    Data data;
    ZERO_OBJECT( &data );
    checkForTest( &argc, argv );
    data.ok = parseArguments( argc, argv, &data.arguments );

    if ( data.ok ) {
      readData( &data );

      if ( data.ok ) {
        streamData( &data );
      }
    }

    ok = data.ok;
    deallocateData( &data );
  }

  POST0( IS_BOOL( ok ) );
  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



#ifndef NO_ASSERTIONS



/******************************************************************************
PURPOSE: isValidArguments - Is arguments in a valid state?
INPUTS:  const Arguments* arguments  Structure to examine.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidArguments( const Arguments* arguments ) {
  const Integer result =
    AND16( arguments,
           arguments->listFile,
           arguments->listFile[ 0 ],
           arguments->tmpdir,
           arguments->tmpdir[ 0 ],
           arguments->description,
           arguments->description[ 0 ],
           isValidTimestamp( arguments->firstTimestamp ),
           arguments->timesteps > 0,
           IS_BOOL( arguments->corners ),
           isValidBounds( arguments->domain ),
           arguments->maximumCloudFraction >= 0.0,
           arguments->maximumCloudFraction <= 1.0,
           arguments->maximumSolarZenithAngle >= 0.0,
           arguments->maximumSolarZenithAngle <= 90.0,
           IS_BOOL( arguments->allowNegativeCounts ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidSwathInfo - Is swathInfo in a valid state?
INPUTS:  const SwathInfo* swathInfo  Structure to examine.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidSwathInfo( const SwathInfo* swathInfo ) {
  const Integer result =
    AND3( swathInfo,
          swathInfo->points > 0,
          isValidTimestamp( swathInfo->timestamp ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidData - Is data in a valid state?
INPUTS:  const Data* data  Structure to examine.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidData( const Data* data ) {
  Integer result =
    AND9( data,
          isValidArguments( &data->arguments ),
          IMPLIES( data->tempFile,
                   AND2( data->tempFileName, data->tempFileName[ 0 ] ) ),
          data->swaths,
          data->swaths->invariant( data->swaths ),
          data->units[ 0 ],
          ! strchr( data->units, ' ' ),
          IMPLIES_ELSE( data->buffer,
                        AND2( data->rows, data->columns ),
                        IS_ZERO2( data->rows, data->columns ) ),
          IS_BOOL( data->ok ) );

  if ( result ) {
    const VoidList* const swaths = data->swaths;
    const Integer swathCount = swaths->count( swaths );
    Integer index = 0;

    do {
      const SwathInfo* const swathInfo = swaths->item( swaths, index );

      if ( ! isValidSwathInfo( swathInfo ) ) {
        result = 0;
        index = swathCount;
      }

      ++index;
    } while ( index < swathCount );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



#endif /* ! defined( NO_ASSERTIONS ) */



/******************************************************************************
PURPOSE: deallocateData - Deallocate and zero contents of data.
INPUTS:  Data* data Data structure to examine and deallocate/zero members.
OUTPUTS: Data* data Data structure with deallocated/zeroed members.
******************************************************************************/

static void deallocateData( Data* data ) {
  PRE0( data );
  FREE_OBJECT( data->swaths );

  if ( data->tempFile ) {
    fclose( data->tempFile ), data->tempFile = 0;
  }

  if ( data->tempFileName[ 0 ] ) {
    unlink( data->tempFileName );
    data->tempFileName[ 0 ] = '\0';
  }

  FREE( data->buffer );
  ZERO_OBJECT( data );
  POST0( data );
}



/******************************************************************************
PURPOSE: deallocateSwathInfo - Deallocate and zero contents of swathInfo.
INPUTS:  void* argument  Data structure to examine and deallocate/zero members.
OUTPUTS: void* argument  Data structure with deallocated/zeroed members.
NOTES:   Called by VoidList, thus the void* argument.
******************************************************************************/

static void deallocateSwathInfo( void* argument ) {
  PRE0( argument );
  SwathInfo* const swathInfo = argument;
  ZERO_OBJECT( swathInfo );
  POST0( argument );
}



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* programName  Name of program.
******************************************************************************/

static void printUsage( const char* programName ) {
  PRE0( programName );
  fprintf( stderr, "\a\n\n%s - Read a set of OMI files and extract swath\n",
           programName );
  fprintf( stderr, "data subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable(s).\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "  -files <listfile> \\\n" );
  fprintf( stderr, "  -tmpdir <temporary_directory> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, "  -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, "  [-maximumCloudFraction value]\\\n" );
  fprintf( stderr, "  [-maximumSolarZenithAngle value]\\\n" );
  fprintf( stderr, "  [-allowNegativeCounts]\\\n" );
  fprintf( stderr, "  -corners \\n\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory to write transient file:\n");
  fprintf( stderr, "-maximumCloudFraction option filter-out values greater " );
  fprintf( stderr, "than the specified value [0.0, 1.0]. Default is 1.0.\n");
  fprintf(stderr, "-maximumSolarZenithAngle option filter-out values greater ");
  fprintf( stderr, "than the specified value [0.0, 90.0]. Default is 90.0.\n");
  fprintf( stderr, "-allowNegativeCounts will allow negative counts of "
                   "molecules/cm2 (non-physical).\n" );
  fprintf( stderr, "-corners option will output 8 additional variables:\n" );
  fprintf( stderr, "  Longitude_SW Longitude_SE Longitude_NW Longitude_NE\n");
  fprintf( stderr, "  Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "that are the linearly interpolated " );
  fprintf( stderr, "(and edge extrapolated)\n" );
  fprintf( stderr, "corner points for each center-pixel point.\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example #1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files file_list.txt \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr,
    "-desc \"https://disc.gdfc.nasa.gov/datasets/OMNO2_004/summary\" \\\n");
  fprintf( stderr, "-timestamp 2005080100 -hours 24 \\\n" );
  fprintf( stderr, "-variable ColumnAmountNO2\\\n" );
  fprintf( stderr, "-domain -51 35 -50 36 -corners > subset.xdr\n\n" );
  fprintf( stderr, "NO2 on August 1 2005.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Swath 2.0\n" );
  fprintf( stderr, "https://disc.gdfc.nasa.gov/datasets/OMNO2_004/summary\n" );
  fprintf( stderr, "2005-08-01T00:00:00-0000\n" );
  fprintf( stderr, "# Dimensions: variables timesteps scans:\n" );
  fprintf( stderr, "11 24 2\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Longitude Latitude ColumnAmountNO2 Longitude_SW Longitude_SE Longitude_NW Longitude_NE Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "deg deg molecules/cm2 deg deg deg deg deg deg deg deg\n" );
  fprintf( stderr, "# Domain: <min_lon> <min_lat> <max_lon> <max_lat>\n" );
  fprintf( stderr, "-51 35 -50 36\n" );
  fprintf( stderr, "# MSB 64-bit integers (yyyydddhhmm) timestamps[scans] and\n" );
  fprintf( stderr, "# MSB 64-bit integers points[scans] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals data_1[variables][points_1] ... data_S[variables][points_S]:\n" );
  fprintf( stderr, "<binary data arrays here>\n\n\n");
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: Integer 1 if successful, else 0.
******************************************************************************/

static Integer parseArguments( Integer argc, char* argv[],
                               Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  Integer result = 0;
  Integer arg = 1;
  ZERO_OBJECT( arguments );
  arguments->maximumCloudFraction = 1.0; /* Allow any fraction of clouds. */
  arguments->maximumSolarZenithAngle = 90.0; /* Allow any solar zenith angle */
  arguments->listFile = parseArgument2( argc, argv, "-files", &arg );

  if ( arguments->listFile ) {
    arguments->tmpdir = parseArgument2( argc, argv, "-tmpdir", &arg );

    if ( arguments->tmpdir ) {
      arguments->description = parseArgument2( argc, argv, "-desc", &arg );

      if ( arguments->description ) {

        if ( parseTimestampAndHours( argc, argv, &arg,
                                     &arguments->firstTimestamp,
                                     &arguments->timesteps ) ) {
          arguments->variableName =
            parseArgument2( argc, argv, "-variable", &arg );

          if ( arguments->variableName ) {

            if ( parseDomain( argc, argv, arguments ) ) {

              if ( parseCorners( argc, argv, arguments ) ) {

                if ( parseMaximumCloudFraction( argc, argv, arguments ) ) {

                  if ( parseMaximumSolarZenithAngle( argc, argv, arguments ) ) {

                    if ( parseAllowNegativeCounts( argc, argv, arguments ) ) {
                      result = 1;
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

  if ( ! result ) {
    ZERO_OBJECT( arguments );
    printUsage( argv[ 0 ] );
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments )));
  return result;
}



/******************************************************************************
PURPOSE: parseDomain - Parse command-line arguments for domain, if specified,
         else initialize domain to entire Earth.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: Integer 1 if successful, else 0.
******************************************************************************/

static Integer parseDomain(Integer argc, char* argv[], Arguments* arguments) {

  PRE03( isValidArgs( argc, (const char**) argv ), argc >= 10, arguments );

  Integer result = 0;
  Integer arg = 10;

  while ( AND2( arg < argc, strcmp( argv[ arg ], "-domain" ) ) ) {
    ++arg;
  }

  if ( arg == argc ) { /* Not specified so just use default: */
    arguments->domain[ LONGITUDE ][ MINIMUM ] = -180.0;
    arguments->domain[ LONGITUDE ][ MAXIMUM ] =  180.0;
    arguments->domain[ LATITUDE  ][ MINIMUM ] =  -90.0;
    arguments->domain[ LATITUDE  ][ MAXIMUM ] =   90.0;
    result = 1;
  } else if ( arg + 4 < argc ) {
    arguments->domain[ LONGITUDE ][ MINIMUM ] = atoR( argv[ arg + 1 ] );
    arguments->domain[ LATITUDE  ][ MINIMUM ] = atoR( argv[ arg + 2 ] );
    arguments->domain[ LONGITUDE ][ MAXIMUM ] = atoR( argv[ arg + 3 ] );
    arguments->domain[ LATITUDE  ][ MAXIMUM ] = atoR( argv[ arg + 4 ] );

    if ( ! isValidBounds( (const Real(*)[2]) arguments->domain ) ) {
      failureMessage( "\a\n\nInvalid domain specified [%lg %lg %lg %lg].\n",
               arguments->domain[ LONGITUDE ][ MINIMUM ],
               arguments->domain[ LATITUDE  ][ MINIMUM ],
               arguments->domain[ LONGITUDE ][ MAXIMUM ],
               arguments->domain[ LATITUDE  ][ MAXIMUM ] );
    } else {
      result = 1;
    }
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments )));
  return result;
}



/******************************************************************************
PURPOSE: parseCorners - Parse command-line argument -corners, if specified.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: Integer 1.
******************************************************************************/

static Integer parseCorners(Integer argc, char* argv[], Arguments* arguments) {

  PRE03( isValidArgs( argc, (const char**) argv ), argc >= 10, arguments );

  Integer result = 1;
  Integer arg = 10;

  while ( AND2( arg < argc, strcmp( argv[ arg ], "-corners" ) ) ) {
    ++arg;
  }

  arguments->corners = arg < argc;

  POST02( result == 1, isValidArguments( arguments ) );
  return result;
}



/******************************************************************************
PURPOSE: parseMaximumCloudFraction - Parse command-line argument
         -maximumCloudFraction, if specified.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  arguments->maximumCloudFraction.
RETURNS: Integer 1.
******************************************************************************/

static Integer parseMaximumCloudFraction( Integer argc, char* argv[],
                                          Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  Integer result = 1;
  Integer arg = 1;

  while ( AND2( arg < argc, strcmp( argv[ arg ], "-maximumCloudFraction" ))) {
    ++arg;
  }

  if ( arg < argc - 1 ) {
    const double value = atof( argv[ arg + 1 ] );

    if ( IN_RANGE( value, 0.0, 1.0 ) ) {
      arguments->maximumCloudFraction = value;
    } else {
      failureMessage( "\a\n\nInvalid maximumCloudFraction specified: %g.\n",
                      value );
      result = 0;
    }
  }

  POST0( OR2( result == 0,
              IN_RANGE( arguments->maximumCloudFraction, 0.0, 1.0 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseMaximumSolarZenithAngle - Parse command-line argument
         -maximumSolarZenithAngle, if specified.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  arguments->maximumSolarZenithAngle
RETURNS: Integer 1.
******************************************************************************/

static Integer parseMaximumSolarZenithAngle( Integer argc, char* argv[],
                                             Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  Integer result = 1;
  Integer arg = 1;

  while ( AND2( arg < argc, strcmp( argv[ arg ], "-maximumSolarZenithAngle" ))) {
    ++arg;
  }

  if ( arg < argc - 1 ) {
    const double value = atof( argv[ arg + 1 ] );

    if ( IN_RANGE( value, 0.0, 90.0 ) ) {
      arguments->maximumSolarZenithAngle = value;
    } else {
      failureMessage( "\a\n\nInvalid maximumSolarZenithAngle specified: %g.\n",
                      value );
      result = 0;
    }
  }

  POST0( OR2( result == 0,
              IN_RANGE( arguments->maximumSolarZenithAngle, 0.0, 90.0 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseMaximumAllowNegativeCounts - Parse command-line argument
         -allowNegativeCounts, if specified.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  arguments->allowNegativeCounts
RETURNS: Integer 1.
******************************************************************************/

static Integer parseAllowNegativeCounts( Integer argc, char* argv[],
                                         Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  Integer result = 1;
  Integer arg = 1;

  while ( AND2( arg < argc, strcmp( argv[ arg ], "-allowNegativeCounts" ))) {
    ++arg;
  }

  if ( arg < argc ) {
    arguments->allowNegativeCounts = 1;
  }

  POST0( IS_BOOL( arguments->allowNegativeCounts ) );
  return result;
}



/******************************************************************************
PURPOSE: parseFileTimestamp - Parse timestamp from file name.
INPUTS:  const char* const fileName  Name of data file to parse.
RETURNS: Integer yyyydddhh00  Timestamp of file.
******************************************************************************/

static Integer parseFileTimestamp( const char* fileName ) {
  PRE02( fileName, *fileName );
  Integer result = 0;
  const char* const slash = strrchr( fileName, '/' );
  const char* const tag =
    strstr( fileName, "OMI-Aura_L2-OMNO2_" ) ? "OMI-Aura_L2-OMNO2_"
    : strstr( fileName, "OMI-Aura_L2-OMTO3_" ) ? "OMI-Aura_L2-OMTO3_"
    : strstr( fileName, "OMI-Aura_L2-OMCLDRR_" ) ? "OMI-Aura_L2-OMCLDRR_"
    : "";
  const char* part = strstr( slash ? slash + 1 : fileName, tag );

  if ( part ) {
    int yyyy = 0;
    int mm = 0;
    int dd = 0;
    int hh = 0;
    int parsed = 0;
    part += strlen( tag );
    parsed = sscanf( part, "%04dm%02d%02dt%02d\n", &yyyy, &mm, &dd, &hh );

    if ( parsed == 4 ) {
      const Integer yyyymmdd = yyyy * 10000 + mm * 100 + dd;
      const Integer yyyyddd = convertYearMonthDay( yyyymmdd );
      result = yyyyddd * 100 + hh;
      result *= 100; /* Append 00 minutes. */
    }
  }

  POST0( OR2( result == 0, isValidTimestamp( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readData - Read each OMI file for selected time range and subset
         their swaths to the specified domain.
INPUTS:  Data* data  Data to update.
OUTPUTS: Data* data  Update data->swaths.
******************************************************************************/

static void readData( Data* data ) {

  PRE04( data, data->ok, data->swaths == 0,
         isValidArguments( &data->arguments ) );

  const Arguments* const arguments = &data->arguments;
  const Integer firstTimestamp = arguments->firstTimestamp;
  const Integer lastTimestamp =
    offsetTimestamp( firstTimestamp, arguments->timesteps );
  FILE* listFile = fopen( arguments->listFile, "r" );

  if ( listFile ) {
    data->swaths = newVoidList( deallocateSwathInfo, 0 );

    if ( data->swaths ) {
      char fileName[ 256 ] = "";
      memset( fileName, 0, sizeof fileName );

      while ( fgets( fileName, sizeof fileName / sizeof *fileName, listFile )) {
        char* const newline = strchr( fileName, '\n' );
        int file = 0;
        Integer fileTimestamp = 0;

        if ( newline ) {
          *newline = '\0';
        }

        fileTimestamp = parseFileTimestamp( fileName );
        DEBUG( fprintf( stderr, "%s %lld\n", fileName, fileTimestamp ); )

        if ( IN_RANGE( fileTimestamp, firstTimestamp, lastTimestamp ) ) {
          file = openFile( fileName );

          if ( file != -1 ) {

            /*
             * Read, filter and subset swath data and
             * write subset data to temporary file in tmpdir
             * then store timestamp and number of subset points in swaths list.
             */

            const char* const product =
              strstr( fileName, "OMI-Aura_L2-OMNO2_" ) ? "OMNO2"
              : strstr( fileName, "OMI-Aura_L2-OMTO3_" ) ? "OMTO3"
              : strstr( fileName, "OMI-Aura_L2-OMCLDRR_" ) ? "OMCLDRR"
              : "";
            SwathInfo* swathInfo =
              readSwathData( data, file, product, fileTimestamp );
            data->ok = 1;

            if ( swathInfo ) {
              data->swaths->insert( data->swaths, swathInfo, LAST_ITEM );
            }

            closeFile( file ), file = -1;
          }
        }
      }
    }

    fclose( listFile ), listFile = 0;
  }

  if ( data->tempFile ) { /* Done writing to temp file so close it: */
    fclose( data->tempFile ), data->tempFile = 0;
  }

  data->ok = AND2( data->swaths, data->swaths->count( data->swaths ) > 0 );
  POST0( IMPLIES( data->ok, isValidData( data ) ) );
}



/******************************************************************************
PURPOSE: readSwathData - Read a data file, write subset data to a temporary
         file and return timestamp and number of subset points.
INPUTS:  Data*  data                data->arguments.
         const int file             File to read.
         const char* const product  Name of file type: OMNO2, OMCLDRR, etc.
         const Integer  timestamp   File timestamp (yyyydddhhmm.
OUTPUTS: Data* data   data->ok, data->units.
RETURNS: SwathInfo*  Initialized swath info from file, else 0 if unsuccessful.
******************************************************************************/

static SwathInfo* readSwathData( Data* data,
                                 const int file,
                                 const char* const product,
                                 const Integer timestamp ) {

  PRE08( data, data->ok, data->swaths,
         isValidArguments( &data->arguments ),
         file > -1, product, *product, isValidTimestamp( timestamp ) );

  const Arguments* const arguments = &data->arguments;
  SwathInfo* result = 0;
  size_t rows = 0;
  size_t columns = 0;
  data->ok = readDimensions( file, product, &rows, &columns );

  DEBUG( fprintf( stderr, "readSwathData ( file = %d, product = '%s' ), "
                  "rows = %lu, columns = %lu\n",
                  file, product, rows, columns ); )

  if ( data->ok ) {
    const size_t variables =
      3 + arguments->corners * 8; /* lon,lat,var + corners. */
    const size_t points = rows * columns;
    const size_t size = ( variables + 1 ) * points;
    const int changedDimensions =
      OR2( rows != data->rows, columns != data->columns );

    if ( changedDimensions ) {
      FREE( data->buffer );
      data->buffer = NEW( double, size );

      if ( data->buffer ) {
        data->rows = rows;
        data->columns = columns;
      } else {
        data->rows = 0;
        data->columns = 0;
        data->ok = 0;
      }
    }

    if ( data->buffer ) {
      const int corners = arguments->corners;
      double* const longitudes = data->buffer;
      double* const latitudes  = longitudes + points;
      double* const values     = latitudes  + points;
      double* const temp       = values     + points;
      double* const longitudesSW = corners ? temp         + points : 0;
      double* const longitudesSE = corners ? longitudesSW + points : 0;
      double* const longitudesNW = corners ? longitudesSE + points : 0;
      double* const longitudesNE = corners ? longitudesNW + points : 0;
      double* const latitudesSW  = corners ? longitudesNE + points : 0;
      double* const latitudesSE  = corners ? latitudesSW  + points : 0;
      double* const latitudesNW  = corners ? latitudesSE  + points : 0;
      double* const latitudesNE  = corners ? latitudesNW  + points : 0;
      Name unused = "";
      size_t subsetPoints =
        readDataset( file, rows, columns, product, "Longitude", 1.0, 90.0, 0,
                     unused, longitudes, temp );

      DEBUG( fprintf( stderr, "Longitude subsetPoints = %lu\n", subsetPoints);)

      if ( subsetPoints ) {
        subsetPoints =
          readDataset( file, rows, columns, product, "Latitude", 1.0, 90.0, 0,
                       unused, latitudes, temp );

        DEBUG(fprintf( stderr, "Latitude subsetPoints = %lu\n", subsetPoints);)

        if ( subsetPoints ) {
          subsetPoints = clampInvalidCoordinates(points, longitudes, latitudes);

          if ( subsetPoints ) {
            subsetPoints =
              readDataset( file, rows, columns, product,
                           arguments->variableName,
                           arguments->maximumCloudFraction,
                           arguments->maximumSolarZenithAngle,
                           arguments->allowNegativeCounts,
                           data->units,
                           values, temp );

            DEBUG( fprintf( stderr, "%s (%g, %g, %d) subsetPoints = %lu\n",
                            arguments->variableName,
                            arguments->maximumCloudFraction,
                            arguments->maximumSolarZenithAngle,
                            arguments->allowNegativeCounts,
                            subsetPoints ); )

            if ( subsetPoints ) {

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

              DEBUG(fprintf( stderr, "  subsetPoints = %lu\n", subsetPoints );)

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
                  result = NEW_ZERO( SwathInfo, 1 );
                  data->ok = result != 0;
                  DEBUG(fprintf(stderr, "  swathInfo      = %p\n", result);)

                  if ( result ) {
                    result->timestamp = timestamp;
                    result->points    = subsetPoints;

                    DEBUG( fprintf( stderr,
                                    "  timestamp = %lld, points = %lld, "
                                    "units = '%s'\n",
                                    result->timestamp, result->points,
                                    data->units ); )
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if ( AND2( ! data->ok, result ) ) {
    deallocateSwathInfo( result );
    FREE( result );
  }

  DEBUG( fprintf( stderr, "data->ok = %lld, result = %p\n", data->ok, result );)
  POST0( IMPLIES( result, AND2( isValidSwathInfo( result ), data->ok ) ) );
  return result;
}



/******************************************************************************
PURPOSE: clampInvalidCoordinates - Clamp invalid longitude-latitude points.
INPUTS:  const size_t points          Points of data.
         double longitudes[ points ]  Longitudes to check.
         double latitudes[  points ]  Latitudes to check.
INPUTS:  double longitudes[ points ]  Clamped Longitudes.
         double latitudes[  points ]  Clamped Latitudes.
RETURNS: int 1 if at least one valid point was found.
******************************************************************************/

static int clampInvalidCoordinates(  const size_t points,
                                     double longitudes[],
                                     double latitudes[] ) {

  int result = 0;

  assert( points );
  assert( longitudes ); assert( latitudes );

  {
    size_t point = 0;
    size_t validPoint = 0; /* Index of first valid point: */

    for ( point = 0; point < points; ++point ) {
      const double longitude = longitudes[ point ];

      if ( IN_RANGE( longitude, -180.0, 180.0 ) ) {
        const double latitude = latitudes[ point ];

        if ( IN_RANGE( latitude, -90.0, 90.0 ) ) {
          result = 1;
          validPoint = point;
          point = points; /* Stop looping. */
        }
      }
    }

    if ( result ) {
      const double longitude = longitudes[ validPoint ];
      const double latitude  = latitudes[  validPoint ];

      /* Clamp all previous points to this first valid point: */

      for ( point = 0; point < validPoint; ++point ) {
        longitudes[ point ] = longitude;
        latitudes[  point ] = latitude;
      }

      /* Clamp all remaining points to the previous valid point: */

      for ( ; point < points; ++point ) {

        if ( IN_RANGE( longitudes[ point ], -180.0, 180.0 ) &&
             IN_RANGE( latitudes[  point ],  -90.0,  90.0 ) ) {
          validPoint = point;
        } else {
          longitudes[ point ] = longitudes[ validPoint ];
          latitudes[  point ] = latitudes[  validPoint ];
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: pointsInSubset - Compute number of points in subset based on domain
         and assign MISSING_VALUE to data values that are outside the subset.
INPUTS:  const Bounds domain                  Domain to subset data to.
         const size_t points                  Points of data.
         const double longitudes[   points ]  Longitudes to check.
         const double latitudes[    points ]  Latitudes to check.
         const double longitudesSW[ points ]  0 or corners to check.
         const double longitudesSE[ points ]  0 or corners to check.
         const double longitudesNW[ points ]  0 or corners to check.
         const double longitudesNE[ points ]  0 or corners to check.
         const double latitudesSW[  points ]  0 or corners to check.
         const double latitudesSE[  points ]  0 or corners to check.
         const double latitudesNW[  points ]  0 or corners to check.
         const double latitudesNE[  points ]  0 or corners to check.
OUTPUTS: double values[ points ]              Values filtered by subset.
RETURNS: size_t number of points in subset.
******************************************************************************/

static size_t pointsInSubset( const Bounds domain,
                              const size_t points,
                              const double longitudes[],
                              const double latitudes[],
                              double values[],
                              const double longitudesSW[],
                              const double longitudesSE[],
                              const double longitudesNW[],
                              const double longitudesNE[],
                              const double latitudesSW[],
                              const double latitudesSE[],
                              const double latitudesNW[],
                              const double latitudesNE[] ) {

  PRE06( isValidBounds( domain ), points, longitudes, latitudes, values,
         IMPLIES_ELSE( longitudesSW,
                       NON_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                  latitudesSW, latitudesSE,
                                  latitudesNW, latitudesNE ),
                       IS_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                 latitudesSW, latitudesSE,
                                 latitudesNW, latitudesNE ) ) );

  size_t result = 0;
  const double longitudeMinimum = domain[ LONGITUDE ][ MINIMUM ];
  const double longitudeMaximum = domain[ LONGITUDE ][ MAXIMUM ];
  const double latitudeMinimum  = domain[ LATITUDE  ][ MINIMUM ];
  const double latitudeMaximum  = domain[ LATITUDE  ][ MAXIMUM ];
  size_t point = 0;

  for ( point = 0; point < points; ++point ) {
    const double value     = values[     point ];
    const double longitude = longitudes[ point ];
    const double latitude  = latitudes[  point ];
    int valid =
      AND3( value > MISSING_VALUE,
            IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
            IN_RANGE( latitude,  latitudeMinimum,  latitudeMaximum ) );

    if ( AND2( valid, longitudesSW ) ) { /* Check for degenerate cells: */
      const double longitudeSW = longitudesSW[ point ];
      const double longitudeSE = longitudesSE[ point ];
      const double longitudeNW = longitudesNW[ point ];
      const double longitudeNE = longitudesNE[ point ];
      const double latitudeSW  = latitudesSW[  point ];
      const double latitudeSE  = latitudesSE[  point ];
      const double latitudeNW  = latitudesNW[  point ];
      const double latitudeNE  = latitudesNE[  point ];

      valid =
        AND20( longitudeSW != longitude,
               longitudeSW != longitudeSE,
               longitudeSW != longitudeNW,
               longitudeSW != longitudeNE,
               longitudeSE != longitude,
               longitudeSE != longitudeNW,
               longitudeSE != longitudeNE,
               longitudeNW != longitude,
               longitudeNW != longitudeNE,
               longitudeNE != longitude,
               latitudeSW != latitude,
               latitudeSW != latitudeSE,
               latitudeSW != latitudeNW,
               latitudeSW != latitudeNE,
               latitudeSE != latitude,
               latitudeSE != latitudeNW,
               latitudeSE != latitudeNE,
               latitudeNW != latitude,
               latitudeNW != latitudeNE,
               latitudeNE != latitude );
    }

    if ( valid ) {
      ++result;
    } else {
      values[ point ] = MISSING_VALUE;
    }
  }

  POST0( result <= points );
  return result;
}



/******************************************************************************
PURPOSE: compactSubsetData - Copy valid data points to the first subsetPoints
         elements of the arrays.
INPUTS:  const size_t subsetPoints      Number of valid points.
         const size_t points            Number of total points.
         double longitudes[   points ]  Longitudes to compact.
         double latitudes[    points ]  Latitudes to compact.
         double values[       points ]  Values to compact.
         double longitudesSW[ points ]  0 or corners to compact.
         double longitudesSE[ points ]  0 or corners to compact.
         double longitudesNW[ points ]  0 or corners to compact.
         double longitudesNE[ points ]  0 or corners to compact.
         double latitudesSW[  points ]  0 or corners to compact.
         double latitudesSE[  points ]  0 or corners to compact.
         double latitudesNW[  points ]  0 or corners to compact.
         double latitudesNE[  points ]  0 or corners to compact.
OUTPUTS: double longitudes[   subsetPoints ]  Longitudes after compacting.
         double latitudes[    subsetPoints ]  Latitudes after compacting.
         double values[       subsetPoints ]  Values after compacting.
         double longitudesSW[ subsetPoints ]  0 or corners after compacting.
         double longitudesSE[ subsetPoints ]  0 or corners after compacting.
         double longitudesNW[ subsetPoints ]  0 or corners after compacting.
         double longitudesNE[ subsetPoints ]  0 or corners after compacting.
         double latitudesSW[  subsetPoints ]  0 or corners after compacting.
         double latitudesSE[  subsetPoints ]  0 or corners after compacting.
         double latitudesNW[  subsetPoints ]  0 or corners after compacting.
         double latitudesNE[  subsetPoints ]  0 or corners after compacting.
******************************************************************************/

static void compactSubsetData( const size_t subsetPoints,
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

  PRE06( subsetPoints, points > subsetPoints,
         longitudes, latitudes, values,
         IMPLIES_ELSE( longitudesSW,
                       NON_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                  latitudesSW, latitudesSE,
                                  latitudesNW, latitudesNE ),
                       IS_ZERO7( longitudesSE, longitudesNW, longitudesNE,
                                 latitudesSW, latitudesSE,
                                 latitudesNW, latitudesNE ) ) );

  size_t input = 0;
  size_t output = 0;

  for ( input = 0; input < points; ++input ) {
    const double value = values[ input ];
    const int valid = value > MISSING_VALUE;
    CHECK( output <= input );

    if ( valid ) {

      if ( output < input ) { /* Copy valid values to front of arrays: */
        longitudes[ output ] = longitudes[ input ];
        latitudes[  output ] = latitudes[  input ];
        values[     output ] = value;

        if ( longitudesSW ) {
          longitudesSW[ output ] = longitudesSW[ input ];
          longitudesSE[ output ] = longitudesSE[ input ];
          longitudesNW[ output ] = longitudesNW[ input ];
          longitudesNE[ output ] = longitudesNE[ input ];
          latitudesSW[  output ] = latitudesSW[  input ];
          latitudesSE[  output ] = latitudesSE[  input ];
          latitudesNW[  output ] = latitudesNW[  input ];
          latitudesNE[  output ] = latitudesNE[  input ];
        }
      }

      ++output;
    }
  }

  CHECK( output == subsetPoints );
  POST014( isValidLongitude( longitudes[ 0 ] ),
           isValidLongitude( longitudes[ subsetPoints - 1 ] ),
           isValidLatitude( latitudes[ 0 ] ),
           isValidLatitude( latitudes[ subsetPoints - 1 ] ),
           values[ 0 ] > MISSING_VALUE,
           values[ subsetPoints - 1 ] > MISSING_VALUE,
           IMPLIES( longitudesSW,
                    AND2( isValidLongitude( longitudesSW[ 0 ] ),
                          isValidLongitude( longitudesSW[ subsetPoints - 1]))),
           IMPLIES( longitudesSE,
                    AND2( isValidLongitude( longitudesSE[ 0 ] ),
                          isValidLongitude( longitudesSE[ subsetPoints - 1]))),
           IMPLIES( longitudesNW,
                    AND2( isValidLongitude( longitudesNW[ 0 ] ),
                          isValidLongitude( longitudesNW[ subsetPoints - 1]))),
           IMPLIES( longitudesNE,
                    AND2( isValidLongitude( longitudesNE[ 0 ] ),
                          isValidLongitude( longitudesNE[ subsetPoints - 1]))),
           IMPLIES( latitudesSW,
                    AND2( isValidLatitude( latitudesSW[ 0 ] ),
                          isValidLatitude( latitudesSW[ subsetPoints - 1 ] ))),
           IMPLIES( latitudesSE,
                    AND2( isValidLatitude( latitudesSE[ 0 ] ),
                          isValidLatitude( latitudesSE[ subsetPoints - 1 ] ))),
           IMPLIES( latitudesNW,
                    AND2( isValidLatitude( latitudesNW[ 0 ] ),
                          isValidLatitude( latitudesNW[ subsetPoints - 1 ]))),
           IMPLIES( latitudesNE,
                    AND2( isValidLatitude( latitudesNE[ 0 ] ),
                          isValidLatitude( latitudesNE[ subsetPoints - 1 ]))));
}



/******************************************************************************
PURPOSE: computeCorners - Compute corner vertices given quadrillateral centers.
         const size_t rows                            Rows of data.
         const size_t columns                         Columns of data.
         const double longitudes[   rows * columns ]  Longitudes of centers.
         const double latitudes[    rows * columns ]  Latitudes of centers.
OUTPUTS: double longitudesSW[ rows * columns ]  Cell corner vertices.
         double longitudesSE[ rows * columns ]  Cell corner vertices.
         double longitudesNW[ rows * columns ]  Cell corner vertices.
         double longitudesNE[ rows * columns ]  Cell corner vertices.
         double latitudesSW[  rows * columns ]  Cell corner vertices.
         double latitudesSE[  rows * columns ]  Cell corner vertices.
         double latitudesNW[  rows * columns ]  Cell corner vertices.
         double latitudesNE[  rows * columns ]  Cell corner vertices.
******************************************************************************/

static void computeCorners( const size_t rows,
                            const size_t columns,
                            const double longitudes[],
                            const double latitudes[],
                            double longitudesSW[],
                            double longitudesSE[],
                            double longitudesNW[],
                            double longitudesNE[],
                            double latitudesSW[],
                            double latitudesSE[],
                            double latitudesNW[],
                            double latitudesNE[] ) {

  const size_t rows_1    = rows - 1;
  const size_t columns_1 = columns - 1;
  const size_t cells     = rows * columns;
  size_t cell = 0;
  size_t index = 0;

  assert( rows != 0 ); assert( columns != 0 );
  assert( longitudes ); assert( latitudes );
  assert( longitudesSW ); assert( longitudesSE );
  assert( longitudesNW ); assert( longitudesNE );
  assert( latitudesSW ); assert( latitudesSE );
  assert( latitudesNW ); assert( latitudesNE );

#ifndef NDEBUG

  /* Init corners to MISSING_VALUE to ensure they are written exactly once: */

  for ( cell = 0; cell < cells; ++cell ) {
    longitudesSW[ cell ] = MISSING_VALUE;
    longitudesSE[ cell ] = MISSING_VALUE;
    longitudesNW[ cell ] = MISSING_VALUE;
    longitudesNE[ cell ] = MISSING_VALUE;
    latitudesSW[  cell ] = MISSING_VALUE;
    latitudesSE[  cell ] = MISSING_VALUE;
    latitudesNW[  cell ] = MISSING_VALUE;
    latitudesNE[  cell ] = MISSING_VALUE;
  }

#endif

  if ( OR2( rows < 2, columns < 2 ) ) {

    /* Copy all center values to the corners in such degenerate cases: */

#pragma omp parallel for

    for ( cell = 0; cell < cells; ++cell ) {
      longitudesSW[ cell ] =
      longitudesSE[ cell ] =
      longitudesNW[ cell ] =
      longitudesNE[ cell ] = longitudes[ cell ];
      latitudesSW[ cell ] =
      latitudesSE[ cell ] =
      latitudesNW[ cell ] =
      latitudesNE[ cell ] = latitudes[ cell ];
    }

  } else { /* Linearly interpolate and extrapolate the corner points: */
    size_t row    = 0;
    size_t column = 0;

    /*
     * First compute linearly interpolated corners of all interior cells:
     * Note: rows increase north to south and columns increase west to east.
     */

#pragma omp parallel for private( column )

    for ( row = 0; row < rows_1; ++row ) {
      const size_t rowOffset = row * columns;
      const size_t nextRowOffset = rowOffset + columns;

      /* Interior row, interior columns: */

      for ( column = 0; column < columns_1; ++column ) {
        const size_t thisIndex         = rowOffset + column;
        const size_t nextColumn        = thisIndex + 1;
        const size_t nextRow           = nextRowOffset + column;
        const size_t nextRowNextColumn = nextRow + 1;

        const double longitude            = longitudes[ thisIndex ];
        double nextColumnLongitude        = longitudes[ nextColumn ];
        double nextRowLongitude           = longitudes[ nextRow ];
        double nextRowNextColumnLongitude = longitudes[ nextRowNextColumn ];

        const double latitude                  = latitudes[ thisIndex ];
        const double nextColumnLatitude        = latitudes[ nextColumn ];
        const double nextRowLatitude           = latitudes[ nextRow ];
        const double nextRowNextColumnLatitude = latitudes[ nextRowNextColumn];

        if ( OR2( longitude < -179.0, longitude > 179.0 ) ) {
          clampLongitudes( longitude,
                           &nextColumnLongitude,
                           &nextRowLongitude,
                           &nextRowNextColumnLongitude,
                           &nextRowNextColumnLongitude );
        }

        {
          const double interpolatedLongitude = 0.25 *
            ( longitude + nextColumnLongitude + nextRowLongitude +
              nextRowNextColumnLongitude );

          const double interpolatedLatitude = 0.25 *
            ( latitude + nextColumnLatitude + nextRowLatitude +
              nextRowNextColumnLatitude );

          assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                           SIGN( interpolatedLongitude ) == SIGN( longitude)));

          assert( longitudesNE[ thisIndex         ] == MISSING_VALUE );
          assert( longitudesNW[ nextColumn        ] == MISSING_VALUE );
          assert( longitudesSE[ nextRow           ] == MISSING_VALUE );
          assert( longitudesSW[ nextRowNextColumn ] == MISSING_VALUE );
          assert( latitudesNE[  thisIndex         ] == MISSING_VALUE );
          assert( latitudesNW[  nextColumn        ] == MISSING_VALUE );
          assert( latitudesSE[  nextRow           ] == MISSING_VALUE );
          assert( latitudesSW[  nextRowNextColumn ] == MISSING_VALUE );

          longitudesNE[ thisIndex         ] = interpolatedLongitude;
          longitudesNW[ nextColumn        ] = interpolatedLongitude;
          longitudesSE[ nextRow           ] = interpolatedLongitude;
          longitudesSW[ nextRowNextColumn ] = interpolatedLongitude;

          latitudesNE[ thisIndex         ] = interpolatedLatitude;
          latitudesNW[ nextColumn        ] = interpolatedLatitude;
          latitudesSE[ nextRow           ] = interpolatedLatitude;
          latitudesSW[ nextRowNextColumn ] = interpolatedLatitude;
        }
      } /* End loop on interior columns. */

    } /* End parallel loop on interior rows. */

    /* Serial region (not worth parallelizing): */

    /* Last row, interior columns (extrapolated top edge): */

    for ( column = 1, index = rows_1 * columns + 1; column < columns;
          ++column, ++index ) {
      const size_t previousColumn = index - 1;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousColumnLongitude = longitudes[ previousColumn ];
      const int signPreviousColumnLongitude = SIGN( previousColumnLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousColumnLongitude == signLongitude );

      assert( longitudesNW[ index          ] == MISSING_VALUE );
      assert( longitudesNE[ previousColumn ] == MISSING_VALUE );
      assert( latitudesNW[  index          ] == MISSING_VALUE );
      assert( latitudesNE[  previousColumn ] == MISSING_VALUE );

      if ( ! ok ) { /* longitude to previousColumnLongitude crosses +/-180: */
        longitudesNW[ index          ] = signLongitude * EDGE;
        longitudesNE[ previousColumn ] = signPreviousColumnLongitude * EDGE;
        latitudesNW[  index          ] = latitude;
        latitudesNE[  previousColumn ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesSW[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousColumnLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

        const double previousColumnLatitude = latitudes[ previousColumn ];
        const double midpointLatitude =
          0.5 * ( latitude + previousColumnLatitude );
        const double interpolatedLatitude = latitudesSW[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesNW[ index          ] = extrapolatedLongitude;
        longitudesNE[ previousColumn ] = extrapolatedLongitude;
        latitudesNW[  index          ] = extrapolatedLatitude;
        latitudesNE[  previousColumn ] = extrapolatedLatitude;
      }
    }

    /* First row, interior columns (extrapolated bottom edge): */

    for ( column = 1, index = 1; column < columns; ++column, ++index ) {
      const size_t previousColumn = index - 1;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousColumnLongitude = longitudes[ previousColumn ];
      const int signPreviousColumnLongitude = SIGN( previousColumnLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousColumnLongitude == signLongitude );

      assert( longitudesSW[ index          ] == MISSING_VALUE );
      assert( longitudesSE[ previousColumn ] == MISSING_VALUE );
      assert( latitudesSW[  index          ] == MISSING_VALUE );
      assert( latitudesSE[  previousColumn ] == MISSING_VALUE );

      if ( ! ok ) { /* longitude to previousColumnLongitude crosses +/-180: */
        longitudesSW[ index          ] = signLongitude * EDGE;
        longitudesSE[ previousColumn ] = signPreviousColumnLongitude * EDGE;
        latitudesSW[  index          ] = latitude;
        latitudesSE[  previousColumn ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesNW[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousColumnLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          IN_RANGE( extrapolatedLongitude0, -180.0, 180.0 ) ?
            extrapolatedLongitude0
          : signLongitude * EDGE;

        const double previousColumnLatitude = latitudes[ previousColumn ];
        const double midpointLatitude =
          0.5 * ( latitude + previousColumnLatitude );
        const double interpolatedLatitude = latitudesNW[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesSW[ index          ] = extrapolatedLongitude;
        longitudesSE[ previousColumn ] = extrapolatedLongitude;
        latitudesSW[  index          ] = extrapolatedLatitude;
        latitudesSE[  previousColumn ] = extrapolatedLatitude;
      }
    }

    /* First column, interior rows (extrapolated left edge, except corners): */

    for ( row = 1, index = columns; row < rows; ++row, index += columns ) {
      const size_t previousRow = index - columns;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousRowLongitude = longitudes[ previousRow ];
      const int signPreviousRowLongitude = SIGN( previousRowLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousRowLongitude == signLongitude );

      assert( longitudesSW[ index       ] == MISSING_VALUE );
      assert( longitudesNW[ previousRow ] == MISSING_VALUE );
      assert( latitudesSW[  index       ] == MISSING_VALUE );
      assert( latitudesNW[  previousRow ] == MISSING_VALUE );

      if ( ! ok ) { /* longitude to previousRowLongitude crosses +/-180: */
        longitudesSW[ index       ] = signLongitude * EDGE;
        longitudesNW[ previousRow ] = signPreviousRowLongitude * EDGE;
        latitudesSW[  index       ] = latitude;
        latitudesNW[  previousRow ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesSE[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousRowLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

        const double previousRowLatitude = latitudes[ previousRow ];
        const double midpointLatitude =
          0.5 * ( latitude + previousRowLatitude );
        const double interpolatedLatitude = latitudesSE[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesSW[ index       ] = extrapolatedLongitude;
        longitudesNW[ previousRow ] = extrapolatedLongitude;
        latitudesSW[  index       ] = extrapolatedLatitude;
        latitudesNW[  previousRow ] = extrapolatedLatitude;
      }
    }

    /* Last column, interior rows (extrapolated right edge, except corners): */

    for ( row = 1, index = columns + columns - 1;
          row < rows; ++row, index += columns ) {
      const size_t previousRow = index - columns;

      const double latitude = latitudes[ index ];

      const double longitude = longitudes[ index ];
      const int closeToEdge = OR2( longitude < -179.0, longitude > 179.0 );
      const int signLongitude = SIGN( longitude );
      const double previousRowLongitude = longitudes[ previousRow ];
      const int signPreviousRowLongitude = SIGN( previousRowLongitude );
      const int ok =
        IMPLIES( closeToEdge, signPreviousRowLongitude == signLongitude );

      assert( longitudesSE[ index       ] == MISSING_VALUE );
      assert( longitudesNE[ previousRow ] == MISSING_VALUE );
      assert( latitudesSE[  index       ] == MISSING_VALUE );
      assert( latitudesNE[  previousRow ] == MISSING_VALUE );

      if ( ! ok ) { /* longitude to previousRowLongitude crosses +/-180: */
        longitudesSE[ index       ] = signLongitude * EDGE;
        longitudesNE[ previousRow ] = signPreviousRowLongitude * EDGE;
        latitudesSE[  index       ] = latitude;
        latitudesNE[  previousRow ] = latitude;
      } else {
        const double interpolatedLongitude0 = longitudesSW[ index ];
        const int signInterpolatedLongitude0 = SIGN( interpolatedLongitude0 );
        const int ok2 =
          IMPLIES( closeToEdge, signInterpolatedLongitude0 == signLongitude );
        const double interpolatedLongitude =
          ok2 ? interpolatedLongitude0 : signLongitude * EDGE;
        const double midpointLongitude =
          0.5 * ( longitude + previousRowLongitude );
        const double longitudeDifference =
          midpointLongitude - interpolatedLongitude;
        const double extrapolatedLongitude0 =
          midpointLongitude + longitudeDifference;
        const double extrapolatedLongitude =
          CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

        const double previousRowLatitude = latitudes[ previousRow ];
        const double midpointLatitude = 0.5 * (latitude + previousRowLatitude);
        const double interpolatedLatitude = latitudesSW[ index ];
        const double latitudeDifference =
          midpointLatitude - interpolatedLatitude;
        const double extrapolatedLatitude0 =
          midpointLatitude + latitudeDifference;
        const double extrapolatedLatitude =
          CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

        assert( IMPLIES( closeToEdge,
                         SIGN( interpolatedLongitude ) == signLongitude ) );
        assert( IMPLIES( closeToEdge,
                         SIGN( extrapolatedLongitude ) == signLongitude ) );

        longitudesSE[ index       ] = extrapolatedLongitude;
        longitudesNE[ previousRow ] = extrapolatedLongitude;
        latitudesSE[  index       ] = extrapolatedLatitude;
        latitudesNE[  previousRow ] = extrapolatedLatitude;
      }
    }

    /* First row, first column cell (extrapolated bottom-left corner): */

    {
      const double latitude               = latitudes[ 0 ];
      const double diagonalLatitude       = latitudesNE[ 0 ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ 0 ];
      const double diagonalLongitude      = longitudesNE[ 0 ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesSW[ 0 ] == MISSING_VALUE );
      assert( latitudesSW[  0 ] == MISSING_VALUE );

      longitudesSW[ 0 ] = extrapolatedLongitude;
      latitudesSW[  0 ] = extrapolatedLatitude;
    }

    /* First row, last column cell (extrapolated bottom-right corner): */

    {
      const double latitude               = latitudes[ columns_1 ];
      const double diagonalLatitude       = latitudesNW[ columns_1 ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude   =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ columns_1 ];
      const double diagonalLongitude      = longitudesNW[ columns_1 ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesSE[ columns_1 ] == MISSING_VALUE );
      assert( latitudesSE[  columns_1 ] == MISSING_VALUE );

      longitudesSE[ columns_1 ] = extrapolatedLongitude;
      latitudesSE[  columns_1 ] = extrapolatedLatitude;
    }

    /* Last row, first column cell (extrapolated top-left corner): */

    index = cells - columns;

    {
      const double latitude               = latitudes[ index ];
      const double diagonalLatitude       = latitudesSE[ index ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude   =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ index ];
      const double diagonalLongitude      = longitudesSE[ index ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesNW[ index ] == MISSING_VALUE );
      assert( latitudesNW[  index ] == MISSING_VALUE );

      longitudesNW[ index ] = extrapolatedLongitude;
      latitudesNW[  index ] = extrapolatedLatitude;
    }

    /* Last row, last column cell (extrapolated top-right corner): */

    index = cells - 1;

    {
      const double latitude               = latitudes[ index ];
      const double diagonalLatitude       = latitudesSW[ index ];
      const double latitudeDifference     = latitude  - diagonalLatitude;
      const double extrapolatedLatitude0  = latitude  + latitudeDifference;
      const double extrapolatedLatitude   =
        CLAMPED_TO_RANGE( extrapolatedLatitude0, -90.0, 90.0 );

      const double longitude              = longitudes[ index ];
      const double diagonalLongitude      = longitudesSW[ index ];
      const double longitudeDifference    = longitude - diagonalLongitude;
      const double extrapolatedLongitude0 = longitude + longitudeDifference;
      const double extrapolatedLongitude  =
        CLAMPED_TO_RANGE( extrapolatedLongitude0, -180.0, 180.0 );

      assert( IMPLIES( OR2( longitude < -179.0, longitude > 179.0 ),
                       SIGN( extrapolatedLongitude ) == SIGN( longitude ) ) );

      assert( longitudesNE[ index ] == MISSING_VALUE );
      assert( latitudesNE[  index ] == MISSING_VALUE );

      longitudesNE[ index ] = extrapolatedLongitude;
      latitudesNE[  index ] = extrapolatedLatitude;
    }

    /* Clamp any out-of-range values: */

#pragma omp parallel for

    for ( cell = 0; cell < cells; ++cell ) {
      const double longitude = longitudes[ cell ];

      longitudesNW[cell] = CLAMPED_TO_RANGE(longitudesNW[cell], -180.0,180.0);
      longitudesSW[cell] = CLAMPED_TO_RANGE(longitudesSW[cell], -180.0,180.0);
      longitudesSE[cell] = CLAMPED_TO_RANGE(longitudesSE[cell], -180.0,180.0);
      longitudesNE[cell] = CLAMPED_TO_RANGE(longitudesNE[cell], -180.0,180.0);

      latitudesNW[ cell ] = CLAMPED_TO_RANGE( latitudesNW[ cell ], -90.0,90.0);
      latitudesSW[ cell ] = CLAMPED_TO_RANGE( latitudesSW[ cell ], -90.0,90.0);
      latitudesSE[ cell ] = CLAMPED_TO_RANGE( latitudesSE[ cell ], -90.0,90.0);
      latitudesNE[ cell ] = CLAMPED_TO_RANGE( latitudesNE[ cell ], -90.0,90.0);

      if ( OR2( longitude < -179.0, longitude > 179.0 ) ) {
        clampLongitudes( longitude,
                         &longitudesNW[cell],
                         &longitudesSW[cell],
                         &longitudesSE[cell],
                         &longitudesNE[cell] );
      }

      /* Check for bogus (stretched) cells and collapse them to the center: */

      {
        const double maximumDistance = 3.0; /* Degrees. */
        const double longitudeNW = longitudesNW[ cell ];
        const double longitudeSW = longitudesSW[ cell ];
        const double longitudeSE = longitudesSE[ cell ];
        const double longitudeNE = longitudesNE[ cell ];
        double distanceNW =
          longitudeNW < longitude ? longitude - longitudeNW
          : longitudeNW - longitude;
        double distanceSW =
          longitudeSW < longitude ? longitude - longitudeSW
          : longitudeSW - longitude;
        double distanceNE =
          longitudeNE < longitude ? longitude - longitudeNE
          : longitudeNE - longitude;
        double distanceSE =
          longitudeSE < longitude ? longitude - longitudeSE
          : longitudeSE - longitude;
        const int bogusCell =
          distanceNW > maximumDistance ||
          distanceNE > maximumDistance ||
          distanceSW > maximumDistance ||
          distanceSE > maximumDistance;

        if ( bogusCell ) { /* Collapse bogus cells to the cell center point: */
          longitudesSW[ cell ] =
          longitudesSE[ cell ] =
          longitudesNW[ cell ] =
          longitudesNE[ cell ] = longitude;
          latitudesSW[ cell ] =
          latitudesSE[ cell ] =
          latitudesNW[ cell ] =
          latitudesNE[ cell ] = latitudes[ cell ];
        }
      }
    }

  } /* End else non-degenerate cases: */
}



/******************************************************************************
PURPOSE: clampLongitudes - Clamp cell longitudes and match sign of first one.
INPUTS:  const double longitude First longitude center point.
         double* longitude1     Longitude to check/clamp.
         double* longitude2     Longitude to check/clamp.
         double* longitude3     Longitude to check/clamp.
         double* longitude4     Longitude to check/clamp.
OUTPUTS: double* longitude1     Clamped longitude with same sign as longitude.
         double* longitude2     Clamped longitude with same sign as longitude.
         double* longitude3     Clamped longitude with same sign as longitude.
         double* longitude4     Clamped longitude with same sign as longitude.
******************************************************************************/

static void clampLongitudes( const double longitude,
                             double* longitude1,
                             double* longitude2,
                             double* longitude3,
                             double* longitude4 ) {

  assert( longitude1 );
  assert( longitude2 );
  assert( longitude3 );
  assert( longitude4 );

  if ( longitude < -179.0 ) {

    if ( *longitude1 >= 0.0 ) {
      *longitude1 = -EDGE;
    }

    if ( *longitude2 >= 0.0 ) {
      *longitude2 = -EDGE;
    }

    if ( *longitude3 >= 0.0 ) {
      *longitude3 = -EDGE;
    }

    if ( *longitude4 >= 0.0 ) {
      *longitude4 = -EDGE;
    }

  } else if ( longitude > 179.0 ) {

    if ( *longitude1 <= 0.0 ) {
      *longitude1 = EDGE;
    }

    if ( *longitude2 <= 0.0 ) {
      *longitude2 = EDGE;
    }

    if ( *longitude3 <= 0.0 ) {
      *longitude3 = EDGE;
    }

    if ( *longitude4 <= 0.0 ) {
      *longitude4 = EDGE;
    }
  }

  assert( SIGN( *longitude1 ) == SIGN( longitude ) );
  assert( SIGN( *longitude2 ) == SIGN( longitude ) );
  assert( SIGN( *longitude3 ) == SIGN( longitude ) );
  assert( SIGN( *longitude4 ) == SIGN( longitude ) );
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

  PRE07( data, data->ok, points, longitudes, latitudes, values,
        IMPLIES_ELSE( longitudesSW,
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

  PRE03( isValidData( data ), data->ok, ! data->tempFile );

  streamHeader( data );
  data->ok = streamSwathTimestamps( data->swaths );

  if ( data->ok ) {
    data->ok = streamSwathPoints( data->swaths );

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

  PRE02( isValidData( data ), data->ok );

  const Arguments* const arguments = &data->arguments;
  const int variables = 3 + arguments->corners * 8;
  const Integer scans = data->swaths->count( data->swaths );
  UTCTimestamp timestamp;
  toUTCTimestamp( arguments->firstTimestamp, timestamp );

  printf( "Swath 2.0\n%s\n%s\n", arguments->description, timestamp );
  printf( "# Dimensions: variables timesteps scans:\n%d %lld %lld\n",
          variables, arguments->timesteps, scans );
  printf( "# Variable names:\n" );
  printf( "Longitude Latitude %s", arguments->variableName );

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



/******************************************************************************
PURPOSE: streamSwathTimestamps - Write MSB 64-bit integer swath timestamps to
         stdout.
INPUTS:  const VoidList* swaths  List of SwathInfo*.
RETURNS: int 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static int streamSwathTimestamps( const VoidList* swaths ) {

  PRE03( swaths, swaths->invariant( swaths ), swaths->count( swaths ) > 0 );

  int result = 0;
  const Integer swathCount = swaths->count( swaths );
  Integer swathIndex = 0;

  do {
    const SwathInfo* const swathInfo = swaths->item( swaths, swathIndex );
    CHECK( swathInfo );
    {
      Integer word = swathInfo->timestamp;
      rotate8ByteArrayIfLittleEndian( &word, 1 );
      result = fwrite( &word, sizeof word, 1, stdout ) == 1;
    }

    ++swathIndex;
  } while ( AND2( result, swathIndex < swathCount ) );

  if ( ! result ) {
    failureMessage( "Failed to stream subset swath timestamps." );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: streamSwathPoints - Write MSB 64-bit integer swath subset point counts
         to stdout.
INPUTS:  const VoidList* swaths  List of SwathInfo*.
RETURNS: int 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static int streamSwathPoints( const VoidList* swaths ) {

  PRE03( swaths, swaths->invariant( swaths ), swaths->count( swaths ) > 0 );

  int result = 0;
  const Integer swathCount = swaths->count( swaths );
  Integer swathIndex = 0;

  do {
    const SwathInfo* const swathInfo = swaths->item( swaths, swathIndex );
    CHECK( swathInfo );
    {
      Integer word = swathInfo->points;
      rotate8ByteArrayIfLittleEndian( &word, 1 );
      result = fwrite( &word, sizeof word, 1, stdout ) == 1;
    }

    ++swathIndex;
  } while ( AND2( result, swathIndex < swathCount ) );

  if ( ! result ) {
    failureMessage( "Failed to stream subset swath point counts." );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: streamTempFile - Stream XDR binary subset data (content of tempfile)
         to stdout.
INPUTS:  Data* const data  data->tempFileName.
OUTPUTS: Data* const data  data->ok, tempFile = 0 (closed).
******************************************************************************/

static void streamTempFile( Data* const data ) {

  PRE04( isValidData( data ), data->ok, data->tempFileName[ 0 ],
         data->tempFile == 0 ); /* tempFile closed after writing to it. */

  const size_t bytes = 256 * 1024 * 1024;
  void* buffer = NEW( char*, bytes );
  data->ok = buffer != 0;

  if ( buffer ) {
    data->tempFile = fopen( data->tempFileName, "rb" );
    data->ok = data->tempFile != 0;

    if ( ! data->ok ) {
      failureMessage( "Can't open temp data file '%s' for reading.",
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

    FREE( buffer );
  }

  if ( ! data->ok ) {
    failureMessage( "Failed to stream subset data from temp file '%s'.",
                    data->tempFileName );
  }

  if ( data->tempFile ) {
    fclose( data->tempFile );
    data->tempFile = 0;
  }
}




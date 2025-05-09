
/******************************************************************************
PURPOSE: METARSubset.c - Read a subset of a NWS METAR NetCDF station files
         and write it to stdout as XDR (IEEE-754) format binary or ASCII
         tab-delimited spreadsheet.

NOTES:   To compile: ./makeit.
         To run:
         ../../../bin/$platform/METARSubset \
           -files test/filelist \
           -xdr \
           -desc https://madis.noaa.gov/,METARSubset \
           -time 20160224000000 20160224015959 \
           -variable temperature \
           -units C -scale 1.0 -offset -273.15 -min -50 -max 50 \
           -domain -76 35 -75 40 > subset.xdr
           head -12 subset.xdr
SITE 2.0
https://madis.noaa.gov/,METARSubset
2016-02-24T00:00:00-0000
# data dimensions: timesteps stations
2 16
# Variable names:
temperature
# Variable units:
C
# char notes[stations][80] and
# MSB 64-bit integers ids[stations] and
# IEEE-754 64-bit reals sites[stations][2=<longitude,latitude>] and
# IEEE-754 64-bit reals data[timesteps][stations]:

         ../../../bin/$platform/METARSubset \
           -files test/filelist \
           -ascii \
           -desc https://madis.noaa.gov/,METARSubset \
           -time 20160224000000 20160224015959 \
           -variable temperature \
           -units C -scale 1.0 -offset -273.15 -min -50 -max 50 \
           -domain -76 35 -75 40

Timestamp(UTC)	LONGITUDE(deg)	LATITUDE(deg)	STATION(-)	temperature(C)
2016-02-24T00:00:00-0000	  -75.4800	   37.9300	000072402	   10.0000
2016-02-24T00:00:00-0000	  -75.2500	   39.8700	000072408	    5.0000
2016-02-24T00:00:00-0000	  -75.9000	   35.6800	000074695	   13.0000
2016-02-24T00:00:00-0000	  -75.0700	   39.3700	750703937	    7.0000
2016-02-24T00:00:00-0000	  -75.1200	   38.3000	751203830	    9.0000
2016-02-24T00:00:00-0000	  -75.3700	   38.6800	753703868	    8.0000
2016-02-24T00:00:00-0000	  -75.4700	   39.1300	754703913	    7.0000
2016-02-24T00:00:00-0000	  -75.5200	   38.3300	755203833	    8.0000
2016-02-24T00:00:00-0000	  -75.6000	   39.6800	756003968	    4.0000
2016-02-24T00:00:00-0000	  -75.6200	   35.2200	756203522	   12.0000
2016-02-24T00:00:00-0000	  -75.6700	   36.0200	756703602	    9.0000
2016-02-24T00:00:00-0000	  -75.7000	   35.9200	757003592	   10.0000
2016-02-24T00:00:00-0000	  -75.7700	   37.6500	757703765	    6.0000
2016-02-24T00:00:00-0000	  -75.8700	   39.9800	758703998	    3.0000
2016-02-24T00:00:00-0000	  -75.9300	   38.8000	759303880	    8.0000
2016-02-24T00:00:00-0000	  -76.0000	   37.8300	760003783	    6.0000
2016-02-24T01:00:00-0000	  -75.4800	   37.9300	000072402	   10.0000
2016-02-24T01:00:00-0000	  -75.2500	   39.8700	000072408	    4.0000
2016-02-24T01:00:00-0000	  -75.9000	   35.6800	000074695	   11.0000
2016-02-24T01:00:00-0000	  -75.0700	   39.3700	750703937	    6.0000
2016-02-24T01:00:00-0000	  -75.1200	   38.3000	751203830	    9.0000
2016-02-24T01:00:00-0000	  -75.3700	   38.6800	753703868	    7.0000
2016-02-24T01:00:00-0000	  -75.4700	   39.1300	754703913	    6.0000
2016-02-24T01:00:00-0000	  -75.5200	   38.3300	755203833	    9.0000
2016-02-24T01:00:00-0000	  -75.6000	   39.6800	756003968	    4.0000
2016-02-24T01:00:00-0000	  -75.6200	   35.2200	756203522	   13.0000
2016-02-24T01:00:00-0000	  -75.6700	   36.0200	756703602	    9.0000
2016-02-24T01:00:00-0000	  -75.7000	   35.9200	757003592	   10.0000
2016-02-24T01:00:00-0000	  -75.7700	   37.6500	757703765	    6.0000
2016-02-24T01:00:00-0000	  -75.8700	   39.9800	758703998	    3.0000
2016-02-24T01:00:00-0000	  -75.9300	   38.8000	759303880	    8.0000
2016-02-24T01:00:00-0000	  -76.0000	   37.8300	760003783	    6.0000

HISTORY: 2016-02-25 plessel.todd@epa.gov
STATUS: unreviewed untested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For malloc(), free(), qsort(). */
#include <stdio.h>  /* For FILE*, fwrite(). */
#include <string.h> /* For strrchr(). */
#include <math.h>   /* For sin(), cos(). */
#include <float.h>  /* For DBL_MAX. */

#include <netcdf.h>  /* For NC_* nc_*(). */

#include <Utilities.h> /* For PRE*(), Integer, VoidList, etc. */

/*================================== TYPES ==================================*/

/* Aggregation modes (note AGGREGATE_NONE outputs hourly data (like input)): */

enum { AGGREGATE_NONE, AGGREGATE_ALL, AGGREGATE_DAILY, AGGREGATE_MODES };

#define IS_VALID_AGGREGATE_MODE( mode ) IN_RANGE( mode, 0, AGGREGATE_MODES - 1)

static const char* const aggregateModes[ AGGREGATE_MODES ] = {
  "none", "all", "daily"
};

/* Output formats: */

enum { OUTPUT_XDR, OUTPUT_ASCII, OUTPUT_FORMATS };

#define IS_VALID_OUTPUT_FORMAT(format) IN_RANGE(format, 0, OUTPUT_FORMATS - 1)

enum { NAME_LENGTH = 5, LOCATION_NAME_LENGTH = 24 };

static const char* const outputFormats[ OUTPUT_FORMATS ] = { "-xdr", "-ascii"};

/* User-supplied command-line arguments: */

typedef struct {
  const char*  listFile;       /* File containing list of files to read. */
  const char*  description;    /* User-supplied description. */
  const char*  variable;       /* Name of variable to read. */
  const char*  units;          /* Units of variable to read. */
  double       scale;          /* Scale to multiply data by. */
  double       offset;         /* Value to add to data value. */
  double       minimum;        /* Minimum valid value of variable data. */
  double       maximum;        /* Maximum valid value of variable data. */
  Bounds       bounds;         /* bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]*/
  Integer      firstTimestamp; /* YYYYMMDDHHMMSS of subset. */
  Integer      lastTimestamp;  /* YYYYMMDDHHMMSS of subset. */
  int          outputFormat;   /* OUTPUT_XDR, OUTPUT_ASCII. */
  int          aggregate;      /* AGGREGATE_NONE/ALL/DAILY. */
} Arguments;

#ifndef NO_ASSERTIONS

static int isValidArguments( const Arguments* arguments ) {
  const int result =
    AND20( arguments,
           arguments->listFile,
           arguments->listFile[ 0 ],
           arguments->description,
           arguments->description[ 0 ],
           arguments->variable,
           arguments->variable[ 0 ],
           arguments->units,
           arguments->units[ 0 ],
           IN_RANGE( arguments->scale,  -DBL_MAX, DBL_MAX ),
           arguments->scale != 0.0,
           IN_RANGE( arguments->offset, -DBL_MAX, DBL_MAX ),
           IN_RANGE( arguments->minimum, -DBL_MAX, DBL_MAX ),
           IN_RANGE( arguments->maximum, arguments->minimum, DBL_MAX ),
           isValidBounds( arguments->bounds ),
           isValidYYYYMMDDHHMMSS( arguments->firstTimestamp ),
           isValidYYYYMMDDHHMMSS( arguments->lastTimestamp ),
           arguments->firstTimestamp <= arguments->lastTimestamp,
           IS_VALID_OUTPUT_FORMAT( arguments->outputFormat ),
           IS_VALID_AGGREGATE_MODE( arguments->aggregate ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* NO_ASSERTIONS */

/* Filtered/subsetted file data: */

typedef struct {
  Integer timestamp; /* YYYYMMDDHHMMSS of file. */
  int* ids;          /* Station ids[ count ] (WMO or generated). */
  float* longitudes; /* Station longitudes[ count ]. */
  float* latitudes;  /* Station latitudes[ count ]. */
  float* values;     /* Station values[ count ] of specified variable. */
  float* values2;    /* Station values[ count ] of wind 2nd component or 0. */
  Note* notes;       /* Station notes[ count ] name/location/description. */
  int count;         /* Length of above arrays. */
} SubsetData;

static void deallocateSubsetData( SubsetData* subsetData ) {
  PRE0( subsetData );
  FREE( subsetData->ids );
  FREE( subsetData->longitudes );
  FREE( subsetData->latitudes );
  FREE( subsetData->values );
  FREE( subsetData->values2 );
  FREE( subsetData->notes );
  ZERO_OBJECT( subsetData );
  POST0( IS_ZERO7( subsetData->timestamp,
                   subsetData->ids,
                   subsetData->longitudes,
                   subsetData->latitudes,
                   subsetData->values,
                   subsetData->values2,
                   subsetData->count ) );
}

/* Data output: */

typedef struct {
  Arguments arguments;      /* Command-line arguments. */
  VoidList* subsetDataList; /* List of filtered/subsetted file data. */
  int timesteps;            /* Number of timesteps in subset time range*/
  int stations;             /* Number of unique stations. */
  int* ids;                 /* ids[ stations ] station ids. */
  float* lonlats;           /* lonlats[ stations ][ LONGITUDE LATITUDE ]. */
  float* values;            /* values[ timesteps ][ stations ]. */
  float* values2;           /* values2[timesteps ][ stations ] of windV or 0.*/
  Note* notes;              /* notes[ stations ] station name/location/desc. */
  int ok;                   /* Did last command succeed? */
} Data;

static void deallocateData( Data* data ) {
  PRE0( data );
  FREE_OBJECT( data->subsetDataList );
  FREE( data->ids );
  FREE( data->lonlats );
  FREE( data->values );
  FREE( data->values2 );
  ZERO_OBJECT( data );
  POST0( data );
}

#ifndef NO_ASSERTIONS

static int isValidData( const Data* data ) {
  const int result =
  AND8( data,
        isValidArguments( &data->arguments ),
        IMPLIES( data->subsetDataList,
                 data->subsetDataList->invariant( data->subsetDataList ) ),
        IMPLIES( data->stations > 0, data->ids ),
        IMPLIES( data->stations > 0,
                 AND3( data->notes, data->notes[0],
                       data->notes[ data->stations - 1 ] ) ),
        IMPLIES( data->timesteps > 0, AND2( data->lonlats, data->values ) ),
        IMPLIES( AND2( ! strcmp( data->arguments.variable, "wind" ),
                       data->values ),
                 data->values2 ),
        IS_BOOL( data->ok ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* NO_ASSERTIONS */

static const float missingValue = -9999.0; /* Missing or invalid data value.*/

assert_static( sizeof (long long) == 8 );
assert_static( sizeof (float)     == 4 );
assert_static( sizeof (int)       == 4 );

/*========================== FORWARD DECLARATIONS ===========================*/

static void printUsage( const char* programName );

static int parseArguments( int argc, char* argv[], Arguments* arguments );

static void deallocateData( Data* data );

static void readData( Data* data );

static void readDataFile( const char* const fileName, Data* data );

static int convertWind( const int count,
                        const float windDirection[], const float windSpeed[],
                        float windU[], float windV[] );

static int computeRelativeHumidity( const int count,
                                    const float temperature[],
                                    const float dewPointTemperature[],
                                    float relativeHumidity[] );

static int filterData( const Arguments* const arguments,
                       const int count,
                       int ids[], float longitudes[], float latitudes[],
                       float values[], float values2[],
                       const Note stationNames[], const Note locationNames[],
                       Note notes[]);

static void appendData( const Integer timestamp,
                        const int count, const int subsetCount,
                        const int ids[],
                        const float longitudes[], const float latitudes[],
                        const float values[], const float values2[],
                        const Note notes[],
                        Data* data );

static void consolidateData( Data* data );

static void allocateConsolidatedData( Data* data );

static void copyConsolidatedData( Data* data );

static void copyUniqueStations( Data* data );

static void findStationLongitudeLatitudeNote( const Data* data, const int id,
                                              float* longitude, float* latitude,
                                              Note note );

static float findStationData( const Data* data, const Integer timestamp,
                              const int id, float* value2 );

static int intCompare( const void* pa, const void* pb ) {
  const int* const pai = (const int*) pa;
  const int* const pbi = (const int*) pb;
  const int a = *pai;
  const int b = *pbi;
  const int result = a - b;
  return result;
}

static void writeHeader( const Data* data );

static void writeXDR( Data* data );

static void writeASCII( const Data* data );

static int generateId( const float longitude, const float latitude );

static Integer fileTimestamp( const char* const fileName );

static int openNetCDFFile( const char* const fileName );

static int stationsInNetCDFFile( const int file );

static int readArray( const int file, const int count, const int type,
                      const char* const name, void* values );

static int readStringArray( const int file, const int count, const int length,
                            const char* const name, Note notes[] );

/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Read a subset of data and write it to stdout in
         XDR or ASCII format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int ok = 0;

  if ( ! isValidArgs( argc, (const char**) argv ) ) {
    failureMessage( "Invalid command-line arguments." );
  } else {
    Data data;
    ZERO_OBJECT( &data );

    if ( parseArguments( argc, argv, &data.arguments ) ) {
      readData( &data );

      if ( data.ok ) {

        if ( data.arguments.outputFormat == OUTPUT_XDR ) {
          writeXDR( &data );
        } else {
          CHECK( data.arguments.outputFormat == OUTPUT_ASCII );
          writeASCII( &data );
        }

        ok = data.ok;
      }
    }

    deallocateData( &data );
  }

  POST0( IS_BOOL( ok ) );
  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* programName  Name of program.
******************************************************************************/

static void printUsage( const char* programName ) {
  PRE0( programName );
  fprintf( stderr, "\n\n\n%s - Read a subset of NWS METAR station files\n",
           programName );
  fprintf( stderr, "and write it to stdout in XDR or ASCII format.\n" );
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files <file_name> \\\n" );
  fprintf( stderr, "-xdr | -ascii \\\n" );
  fprintf( stderr, "-desc \"description text\" \\\n" );
  fprintf( stderr, "-time <yyyymmddhhmmss> <yyyymmddhhmmss> \\\n" );
  fprintf( stderr, "-variable <name> -units <name> \\\n" );
  fprintf( stderr, "[-scale <number> -offset <number>]\\\n" );
  fprintf( stderr, "-min <number> -max <number>\\\n" );
  fprintf( stderr, "[ -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ]\n" );
  fprintf( stderr, "[ -aggregate none|all|daily ] \\\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example 1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files /data/tmp/metar_file_list.txt \\\n" );
  fprintf( stderr, "-xdr -desc https://madis.noaa.gov/,METARSubset \\\n" );
  fprintf( stderr, "-time 20160224000000 20160224015959 \\\n" );
  fprintf( stderr, "-variable temperature -units C \\\n" );
  fprintf( stderr, "-scale 1 -offset -273.15 -min -50 -max 50 \\\n" );
  fprintf( stderr, "-domain -76 35 -75 40 > subset.xdr\n" );
  fprintf( stderr, "\n\nOutputs an ASCII header followed by binary arrays:\n");
  fprintf( stderr, "SITE 2.0\n" );
  fprintf( stderr, "https://madis.noaa.gov/,METARSubset\n" );
  fprintf( stderr, "2016-02-24T00:00:00-0000\n" );
  fprintf( stderr, "# data dimensions: timesteps stations:\n" );
  fprintf( stderr, "2 16\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "temperature\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "C\n" );
  fprintf( stderr, "# char notes[stations][80] and\n" );
  fprintf( stderr, "# MSB 64-bit integers ids[stations] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals " );
  fprintf( stderr, "sites[stations][2=<longitude,latitude>] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals data[timesteps][stations]:\n");
  fprintf( stderr, "<binary data arrays here>\n\n\n");
  fprintf( stderr, "Example 2:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files /data/tmp/metar_file_list.txt \\\n" );
  fprintf( stderr, "-ascii -desc https://madis.noaa.gov/,METARSubset \\\n" );
  fprintf( stderr, "-time 20160224000000 20160224015959 \\\n" );
  fprintf( stderr, "-variable temperature -units C \\\n" );
  fprintf( stderr, "-scale 1 -offset -273.15 -min -50 -max 50 \\\n" );
  fprintf( stderr, "-domain -76 35 -75 40 > subset.xdr\n" );
  fprintf( stderr, "\n\nOutputs an ASCII spreadsheet (tab-delimited):\n" );
  fprintf( stderr,
  "Timestamp(UTC)	LONGITUDE(deg)	LATITUDE(deg)	STATION(-)	temperature(C)\n");
  fprintf( stderr,
  "2016-02-24T00:00:00-0000	  -75.4800	   37.9300	000072402	   10.0000\n" );
  fprintf( stderr,
  "2016-02-24T00:00:00-0000	  -75.2500	   39.8700	000072408	    5.0000\n" );
  fprintf( stderr, "...\n" );
  fprintf( stderr, "\n\n\n");
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  int argc              Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: int 1 if successful, else 0.
NOTES:   Command-line options look like this:
         METARSubset
         -files /data/tmp/metarserver.32345
         -xdr
         -desc https://madis.noaa.gov/,METARSubset
         -time 20160224000000 20160224235959
         -variable temperature
         -units C
         -scale 1
         -offset -273.15
         -min -50
         -max 50
         -bounds -130 20 -60 50
         -aggregate daily
******************************************************************************/

static int parseArguments( int argc, char* argv[], Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  int result = 0;

  checkForTest( &argc, argv );
  ZERO_OBJECT( arguments );
  arguments->scale = 1.0;
  arguments->offset = 0.0;
  arguments->minimum = -DBL_MAX;
  arguments->minimum = DBL_MAX;
  arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
  arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
  arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
  arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
  arguments->outputFormat = -1;

  if ( IN_RANGE( argc, 13, 28 ) ) {
    int arg = 1;

    if ( ! strcmp( argv[ arg ], "-files" ) ) {
      ++arg;
      arguments->listFile = (const char*) argv[ arg ];
      ++arg;

      arguments->outputFormat =
        indexOfString( argv[ arg ], outputFormats,
                       sizeof outputFormats / sizeof *outputFormats );

      if ( ! IS_VALID_OUTPUT_FORMAT( arguments->outputFormat ) ) {
        failureMessage( "Invalid output format '%s'.", argv[ arg ] );
      } else {
        ++arg;

        if ( ! strcmp( argv[ arg ], "-desc" ) ) {
          ++arg;
          arguments->description = argv[ arg ];
          ++arg;

          if ( ! strcmp( argv[ arg ], "-time" ) ) {
            ++arg;

            if ( parseTimeRange( argv[ arg ], argv[ arg + 1 ],
                                 &arguments->firstTimestamp,
                                 &arguments->lastTimestamp ) ) {
              arg += 2;

              if ( ! strcmp( argv[ arg ], "-variable" ) ) {
                ++arg;
                arguments->variable = argv[ arg ];
                ++arg;

                if ( ! strcmp( argv[ arg ], "-units" ) ) {
                  ++arg;
                  arguments->units = argv[ arg ];
                  ++arg;
                  result = 1;

                  /* Parse optional arguments: */

                  while ( AND2( result, arg < argc ) ) {

                    if ( AND2( ! strcmp( argv[ arg ], "-scale" ),
                               arg + 1 < argc ) ) {
                      ++arg;
                      arguments->scale = atof( argv[ arg ] );

                      if ( OR2( arguments->scale == 0.0,
                                ! IN_RANGE( arguments->scale,
                                            -DBL_MAX, DBL_MAX ) ) ) {
                        failureMessage( "Invalid scale: '%s'.", argv[ arg ] );
                        result = 0;
                      }

                    } else if ( AND2( ! strcmp( argv[ arg ], "-offset" ),
                               arg + 1 < argc ) ) {
                      ++arg;
                      arguments->offset = atof( argv[ arg ] );

                      if ( ! IN_RANGE( arguments->offset, -DBL_MAX, DBL_MAX)) {
                        failureMessage( "Invalid offset: '%s'.", argv[ arg ] );
                        result = 0;
                      }

                    } else if ( AND2( ! strcmp( argv[ arg ], "-min" ),
                               arg + 1 < argc ) ) {
                      ++arg;
                      arguments->minimum = atof( argv[ arg ] );
                    } else if ( AND2( ! strcmp( argv[ arg ], "-max" ),
                               arg + 1 < argc ) ) {
                      ++arg;
                      arguments->maximum = atof( argv[ arg ] );
                    } else if ( ! strcmp( argv[ arg ], "-domain" ) ) {
                      Integer skip = arg;
                      result =
                        parseBounds( argc, argv, &skip, arguments->bounds );
                      arg = skip - 1;
                    } else if ( AND2( ! strcmp( argv[ arg ], "-aggregate" ),
                                      arg + 1 < argc ) ) {
                      ++arg;
                      arguments->aggregate =
                        indexOfString( argv[ arg ], aggregateModes,
                                       sizeof aggregateModes /
                                       sizeof *aggregateModes );

                      if ( ! IS_VALID_AGGREGATE_MODE( arguments->aggregate ) ) {
                        failureMessage( "Invalid aggregate mode '%s'.",
                                        argv[ arg ] );
                        result = 0;
                      }

                    } else {
                      failureMessage( "Invalid option: '%s'.", argv[ arg ] );
                      result = 0;
                    }

                    ++arg;
                  }

                  if ( result ) {
                    result = AND2( IN_RANGE( arguments->minimum,
                                             -DBL_MAX, DBL_MAX ),
                                   IN_RANGE( arguments->maximum,
                                             arguments->minimum,
                                             DBL_MAX ) );

                    if ( ! result ) {
                      failureMessage( "Invalid -min -max options: "
                                      "%lf, %lf.",
                                      arguments->minimum,
                                      arguments->maximum );
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
PURPOSE: readData - Read data files.
INPUTS:  Data* data  data->listFile is name of file containing time-sorted
                     list of data files to read.
OUTPUTS: Data* data  data->subsetDataList is list of subset data read.
                     data->ok = 1 if at least one valid subset value was read.
******************************************************************************/

static void readData( Data* data ) {
  PRE02( data, isValidArguments( &data->arguments ) );
  FILE* listFile = fopen( data->arguments.listFile, "r" );

  if ( listFile ) {
    char fileName[ 256 ] = "";
    memset( fileName, 0, sizeof fileName );

    while ( fgets( fileName, sizeof fileName / sizeof *fileName, listFile ) ) {
      char* const newline = strchr( fileName, '\n' );

      if ( newline ) {
        *newline = '\0';
      }

      readDataFile( fileName, data );
    }

    fclose( listFile ), listFile = 0;

    if ( data->ok ) {
      consolidateData( data );
    }
  }

  POST0( isValidData( data ) );
}



/******************************************************************************
PURPOSE: readDataFile - Read named data file.
INPUTS:  const char* const fileName  Name of data file to read.
         Data* data  data->listFile is name of file containing time-sorted
                     list of data files to read.
OUTPUTS: Data* data  data->subsetDataList is list of filtered/subset data read.
                     data->ok = 1 if at least one valid subset value was read.
******************************************************************************/

static void readDataFile( const char* const fileName, Data* data ) {
  PRE04( fileName, fileName[ 0 ], data, isValidArguments( &data->arguments ));

  const Arguments* const arguments = &data->arguments;
  const Integer firstTimestamp = arguments->firstTimestamp;
  const Integer lastTimestamp  = arguments->lastTimestamp;
  const Integer timestamp = fileTimestamp( fileName );

  DEBUG( fprintf( stderr, "readDataFile( fileName = '%s' ) timestamp = %lld\n",
                  fileName, timestamp ); )

  if ( AND2( timestamp, IN_RANGE( timestamp, firstTimestamp, lastTimestamp))) {
    int file = openNetCDFFile( fileName );

    if ( file > -1 ) {
      const int stations = stationsInNetCDFFile( file );

      if ( stations > 0 ) {
        const int isWind = ! strcmp( arguments->variable, "wind" );
        const int isRelativeHumidity =
          ! strcmp( arguments->variable, "relativeHumidity" );
        float* buffer =
          NEW_ZERO( float, stations * ( 4 + isWind + isRelativeHumidity ) );
        Note* notes = buffer ? NEW_ZERO( Note, stations * 3 ) : 0;

        if ( notes ) {
          int* const ids          = (int*) buffer;
          float* const longitudes = buffer     + stations;
          float* const latitudes  = longitudes + stations;
          float* const values     = latitudes  + stations;
          float* values2 =
            OR2( isWind, isRelativeHumidity ) ? values + stations : 0;
          Note* const stationNames = notes + stations;
          Note* const locationNames = stationNames + stations;

          if ( readArray( file, stations, NC_INT, "wmoId", ids ) ) {

            if ( readArray( file, stations, NC_FLOAT, "longitude",
                            longitudes ) ) {

              if ( readArray( file, stations, NC_FLOAT, "latitude",
                              latitudes ) ) {
                int ok = 0;

                if ( isWind ) {
                  ok = readArray( file, stations, NC_FLOAT, "windDir", values );

                  if ( ok ) {
                    ok = readArray( file, stations, NC_FLOAT, "windSpeed",
                                    values2 );

                    if ( ok ) {
                      ok = convertWind( stations, values, values2,
                                        values, values2 );
                    }
                  }
                } else if ( isRelativeHumidity ) { /* Compute pseudo-variable:*/
                  ok = readArray( file, stations, NC_FLOAT, "temperature", values );

                  if ( ok ) {
                    ok = readArray( file, stations, NC_FLOAT, "dewpoint",
                                    values2 );

                    if ( ok ) {
                      ok = computeRelativeHumidity( stations, values, values2,
                                                    values );
                      values2 = 0; /* values now contains relativeHumidity. */
                    }
                  }
                } else {
                  ok = readArray( file, stations, NC_FLOAT, arguments->variable,
                                  values );
                }

                if ( ok ) {
                  ok = readStringArray( file, stations,
                                        NAME_LENGTH, "stationName",
                                        stationNames );

                  if ( ok ) {
                    ok = readStringArray( file, stations,
                                          LOCATION_NAME_LENGTH, "locationName",
                                          locationNames );

                    if ( ok ) {
                      const int validCount =
                        filterData( arguments, stations,
                                    ids, longitudes, latitudes, values, values2,
                                    (const Note*) stationNames,
                                    (const Note*) locationNames, notes );

                      DEBUG( fprintf( stderr, "filtered count = %d\n",
                                        validCount ); )

                      if ( validCount ) {
                        appendData( timestamp, stations, validCount,
                                    ids, longitudes, latitudes, values,
                                    values2, (const Note*) notes, data );
                      }
                    }
                  }
                }
              }
            }
          }
        }

        FREE( buffer );
        FREE( notes );
      }

      nc_close( file ), file = -1;
    }
  }

  POST0( isValidData( data ) );
}



/******************************************************************************
PURPOSE: convertWind - Convert wind direction/speed to windU, windV.
INPUTS:  const int count                     Number of elements in arrays.
         const float windDirection[ count ]  Wind direction (bearing).
         const float windSpeed[ count ]      Wind speed (m/s).
OUTPUTS: float windU[ count ]                Wind eastward component (m/s).
         float windV[ count ]                Wind northward component (m/s).
RETURNS: int number of valid/converted pairs of values or 0.
******************************************************************************/

static int convertWind( const int count,
                        const float windDirection[], const float windSpeed[],
                        float windU[], float windV[] ) {

  PRE05( count > 0, windDirection, windSpeed, windU, windV);
  int result = 0;
  int index = 0;

#pragma omp parallel for reduction( + : result )

  for ( index = 0; index < count; ++index ) {
    const double direction = windDirection[ index ];
    const double speed     = windSpeed[     index ];
    const int ok =
      AND2( IN_RANGE( direction, 0.0, 360.0 ), IN_RANGE( speed, 0.0, 150.0 ) );

    if ( ok ) {
      const double direction0 = 270.0 - direction;
      const double angleDegrees =
          direction0 <   0.0 ? direction0 + 360.0
        : direction0 > 360.0 ? direction0 - 360.0
        : direction0;
      const double angleRadians = radians( angleDegrees );
      const double u = cos( angleRadians );
      const double v = sin( angleRadians );
      const double scaledU = speed * u;
      const double scaledV = speed * v;
      CHECK2( IN_RANGE( scaledU, -150.0, 150.0 ),
              IN_RANGE( scaledV, -150.0, 150.0 ) );
      windU[ index ] = scaledU;
      windV[ index ] = scaledV;
      ++result;
    } else {
      windU[ index ] = missingValue;
      windV[ index ] = missingValue;
    }
  }

  POST0( IN_RANGE( result, 0, count ) );
  return result;
}



/******************************************************************************
PURPOSE: computeRelativeHumidity - Compute relative humidity (%) from
         temperature (C) and dew-point temperature (C).
INPUTS:  const int count                   Number of elements in arrays.
         const float temperature[ count ]  Air temperature (K).
         const float dewPointTemperature[ count ]  Dewpoint temperature (K).
OUTPUTS: float relativeHumidity[ count ]   Relative humidity (%).
RETURNS: int number of valid/converted pairs of values or 0.
NOTES:   Uses Magnus approximation which is considered valid for
         air temperature in range [0, 60] (degrees C)
         dew point air temperature in range [0, 50] (degrees C)
         https://en.wikipedia.org/wiki/Dew_point
******************************************************************************/

static int computeRelativeHumidity( const int count,
                                    const float temperature[],
                                    const float dewPointTemperature[],
                                    float relativeHumidity[] ) {

  PRE04( count > 0, temperature, dewPointTemperature, relativeHumidity);
  const double kelvinToCelsius = -273.15;
  const double b = 17.625;
  const double c = 243.04;
  int result = 0;
  int index = 0;

#pragma omp parallel for reduction( + : result )

  for ( index = 0; index < count; ++index ) {
    const double theTemperature = temperature[ index ] + kelvinToCelsius;
    const double theDewPointTemperature =
      dewPointTemperature[ index ] + kelvinToCelsius;
    const int ok =
      AND2( IN_RANGE( theTemperature, 0.0, 60.0 ),
            IN_RANGE( theDewPointTemperature, 0.0, 50.0 ) );

    if ( ok ) {
      const double numerator =
        c * b * ( theDewPointTemperature - theTemperature );
      const double denominator =
        ( c + theTemperature ) * ( c + theDewPointTemperature );
      const double theRelativeHumidity = 100.0 * exp( numerator / denominator );
      relativeHumidity[ index ] = theRelativeHumidity;
      ++result;
    } else {
      relativeHumidity[ index ] = missingValue;
    }
  }

  POST0( IN_RANGE( result, 0, count ) );
  return result;
}



/******************************************************************************
PURPOSE: filterData - Filter data arrays by subset by overwriting values
         outside the subset with missingValue.
INPUTS:  const Arguments* arguments   Subset arguments.
         const int count              Length of arrays.
         int ids[ count ]             Station Ids.
         float longitudes[ count ]    Station longitudes.
         float latitudes[ count ]     Station latitudes.
         float values[ count ]        Station data values.
         float values2[ count ]       Station data 2nd component values or 0.
         const Note stationNames[ count ]   Station names.
         const Note locationNames[ count ]  Station location names.
OUTPUTS: int ids[ count ]             Reduced station Ids.
         float longitudes[ count ]    Reduced station longitudes.
         float latitudes[ count ]     Reduced station latitudes.
         float values[ count ]        Reduced station data values.
         float values2[ count ]       Reduced station data values2 if not 0.
         Note notes[ count ]          Reduced station notes.
RETURNS: int Number of values inside subset.
NOTES:   Uses WMO Id of station unless it is not available in which case
         an id is generated from the station longitude and latitude.
         Notes are constructed from id+stationName+locationName.
******************************************************************************/

static int filterData( const Arguments* const arguments,
                       const int count,
                       int ids[], float longitudes[], float latitudes[],
                       float values[], float values2[],
                       const Note stationNames[], const Note locationNames[],
                       Note notes[] ) {

  PRE09( isValidArguments( arguments ),
         count > 0, ids, longitudes, latitudes, values,
         stationNames, locationNames, notes );

  const double scale = arguments->scale;
  const double offset = arguments->offset;
  const double minimum = arguments->minimum;
  const double maximum = arguments->maximum;
  const double minimumLongitude = arguments->bounds[ LONGITUDE ][ MINIMUM ];
  const double maximumLongitude = arguments->bounds[ LONGITUDE ][ MAXIMUM ];
  const double minimumLatitude  = arguments->bounds[ LATITUDE  ][ MINIMUM ];
  const double maximumLatitude  = arguments->bounds[ LATITUDE  ][ MAXIMUM ];
  int result = 0;
  int index = 0;

#pragma omp parallel for reduction( + : result )

  for ( index = 0; index < count; ++index ) {
    const float longitude = longitudes[ index ];
    const float latitude  = latitudes[ index ];
    const float value     = values[ index ];
    const float value2    = values2 ? values2[ index ] : value;
    const float convertedValue  = value  * scale + offset;
    const float convertedValue2 = value2 * scale + offset;
    const int id0         = ids[ index ];
    const int id = id0 > 0 ? id0 : generateId( longitude, latitude );
    const Note* const stationName = stationNames + index;
    const Note* const locationName = locationNames + index;
    Note* const note = notes + index;
    const int valid =
      AND5( id >= 0,
            IN_RANGE( longitude, minimumLongitude, maximumLongitude ),
            IN_RANGE( latitude,  minimumLatitude,  maximumLatitude ),
            IN_RANGE( convertedValue,  minimum, maximum ),
            IN_RANGE( convertedValue2, minimum, maximum ) );

    if ( valid ) {
      const char* const name = stationName[0][0] ? *stationName : "no_name";
      const char* const location =
        locationName[0][0] ? *locationName : "no_location";
      ids[ index ] = id;
      values[ index ] = convertedValue;

      if ( values2 ) {
        values2[ index ] = convertedValue2;
      }

      snprintf( *note, sizeof (Note) / sizeof (char), "%d-%s-%s",
                id, name, location );
      ++result;
    } else {
      ids[ index ]        = (int) missingValue;
      longitudes[ index ] = missingValue;
      latitudes[ index ]  = missingValue;
      values[ index ]     = missingValue;

      if ( values2 ) {
        values2[ index ] = missingValue;
      }

      strncpy( *note, "None", sizeof (Note) / sizeof (char) - 1 );
    }
  }

  POST0( IN_RANGE( result, 0, count ) );
  return result;
}



/******************************************************************************
PURPOSE: appendData - Append filtered/subset data to list.
INPUTS:  const Integer timestamp          YYYYMMDDHHMMSS of file data.
         const int count                  Number of elements in arrays.
         const int subsetCount            Number of valid values in arrays.
         const int ids[ count ]           Station ids.
         const float longitudes[ count ]  Station longitudes.
         const float latitudes[ count ]   Station latitudes.
         const float values[ count ]      Station variable values.
         const float values2[ count ]     Station variable 2nd component or 0.
         const Note notes[ count ]              Station notes.
OUTPUTS: Data* data  data->subsetDataList appended with filtered/subset data.
                     data->ok = 1 if successful, else data->ok = 0 and
                     failureMessage() was called.
******************************************************************************/

static void appendData( const Integer timestamp,
                        const int count, const int subsetCount,
                        const int ids[],
                        const float longitudes[], const float latitudes[],
                        const float values[], const float values2[],
                        const Note notes[],
                        Data* data ) {

  PRE010( isValidYYYYMMDDHHMMSS( timestamp ),
          count > 0,
          IN_RANGE( subsetCount, 1, count ),
          ids, longitudes, latitudes, values, notes, data,
          isValidData( data ) );

  SubsetData* subsetData = NEW_ZERO( SubsetData, 1 );
  data->ok = 0;

  if ( subsetData ) {
    subsetData->notes = NEW_ZERO( Note, subsetCount );
    subsetData->ids = subsetData->notes ? NEW_ZERO( int, subsetCount ) : 0;
    subsetData->longitudes =
      subsetData->ids ? NEW_ZERO( float, subsetCount ) : 0;
    subsetData->latitudes =
      subsetData->longitudes ? NEW_ZERO( float, subsetCount ) : 0;
    subsetData->values =
      subsetData->latitudes ? NEW_ZERO( float, subsetCount ) : 0;
    subsetData->values2 =
      AND2( values2, subsetData->values ) ? NEW_ZERO( float, subsetCount ) : 0;

    if ( AND2( subsetData->values, IMPLIES( values2, subsetData->values2 ) ) ) {
      int index = 0;
      int writeIndex = 0;

      for ( index = 0; index < count; ++index ) {

        if ( ids[ index ] >= 0 ) {
          subsetData->ids[        writeIndex ] = ids[        index ];
          subsetData->longitudes[ writeIndex ] = longitudes[ index ];
          subsetData->latitudes[  writeIndex ] = latitudes[  index ];
          subsetData->values[     writeIndex ] = values[     index ];

          if ( values2 ) {
            subsetData->values2[  writeIndex ] = values2[    index ];
          }

          strncpy( subsetData->notes[ writeIndex ], notes[ index ],
                   sizeof (Note) / sizeof (char) );
          ++writeIndex;
        }
      }

      CHECK( writeIndex == subsetCount );
      subsetData->count = writeIndex;
      subsetData->timestamp = timestamp;

      if ( ! data->subsetDataList ) {
        data->subsetDataList =
          newVoidList( (VoidVisitor) deallocateSubsetData, 0 );
      }

      if ( data->subsetDataList ) {
        data->subsetDataList->insert( data->subsetDataList, subsetData,
                                      LAST_ITEM );
        data->ok = data->subsetDataList->ok( data->subsetDataList );
      }
    }

    if ( ! data->ok ) {
      deallocateSubsetData( subsetData );
      FREE( subsetData );
    }
  }

  POST0( isValidData( data ) );
}



/******************************************************************************
PURPOSE: consolidateData - Copy data from each subset to contiguous arrays
         inserting missingValue for stations not reporting for all timesteps.
INPUTS:  Data* data  data->subsetDataList  Subset data.
OUTPUTS: Data* data  data->timesteps  Number of timesteps in subset.
                     data->stations   Number of unique stations.
                     data->ids[ stations ]  Station ids.
                     data->lonlats[ stations ][ 2 ]  Station lon-lats.
                     data->values[ timesteps ][ stations ]  Station data.
                     data->values2[timesteps ][ stations ]  Station data2 or 0
******************************************************************************/

static void consolidateData( Data* data ) {

  PRE04( isValidData( data ), data->subsetDataList,
         IS_ZERO5( data->stations, data->timesteps,
                   data->ids, data->lonlats, data->values ),
         data->ok );

  allocateConsolidatedData( data );

  if ( data->ok ) {
    copyConsolidatedData( data );
  }

  POST02( isValidData( data ),
          IMPLIES_ELSE( data->ok,
                        AND6( data->subsetDataList == 0,
                              data->stations > 0, data->timesteps > 0,
                              data->ids, data->lonlats, data->values ),
                        AND2( data->subsetDataList,
                              IS_ZERO5( data->stations, data->timesteps,
                                        data->ids, data->lonlats,
                                        data->values ) ) ) );
}



/******************************************************************************
PURPOSE: allocateConsolidatedData - Allocate consolidated data arrays.
INPUTS:  Data* data  data->arguments, data->subsetDataList
OUTPUTS: Data* data  data->timesteps  Number of timesteps in subset.
                     data->stations   Number of unique stations.
                     data->ids[ stations ]  Sorted unique station ids.
                     data->lonlats[ stations ][ 2 ]  Allocated station lon-lats
                     data->values[ timesteps ][ stations ]  Initialized missing
******************************************************************************/

static void allocateConsolidatedData( Data* data ) {

  PRE04( isValidData( data ), data->subsetDataList,
         IS_ZERO5( data->stations, data->timesteps,
                   data->ids, data->lonlats, data->values ),
         data->ok );

  /* Allocate consolidated data: */

  UTCTimestamp firstTimestamp = "";
  UTCTimestamp lastTimestamp = "";
  toUTCTimestamp2( data->arguments.firstTimestamp, firstTimestamp );
  toUTCTimestamp2( data->arguments.lastTimestamp, lastTimestamp );

  {
    const int aggregate = data->arguments.aggregate;
    const int timesteps =
      aggregate == AGGREGATE_ALL ? 1
      : aggregate == AGGREGATE_DAILY ?
        daysInRange( firstTimestamp, lastTimestamp )
      : hoursInRange( firstTimestamp, lastTimestamp );
    copyUniqueStations( data );

    if ( data->ok ) {
      const int stations = data->stations;
      const int isWind = ! strcmp( data->arguments.variable, "wind" );
      CHECK4( timesteps > 0, stations > 0, data->ids, data->lonlats );
      data->values = NEW_ZERO( float, timesteps * stations );
      data->values2 =
        AND2( data->values, isWind ) ? NEW_ZERO( float, timesteps * stations )
        : 0;
      data->ok = AND2( data->values, IMPLIES( isWind, data->values2 ) );

      if ( ! data->ok ) {
        data->stations = 0;
        FREE( data->ids );
        FREE( data->lonlats );
        FREE( data->values );
        FREE( data->values2 );
      } else { /* Initialize values to missingValue: */
        const int count = timesteps * stations;
        int index = 0;

        for ( index = 0; index < count; ++index ) {
          data->values[ index ] = missingValue;

          if ( data->values2 ) {
            data->values2[ index ] = missingValue;
          }
        }

        data->timesteps = timesteps;
      }
    }
  }

  POST02( isValidData( data ),
          IMPLIES_ELSE( data->ok,
                        AND8( data->subsetDataList,
                              data->stations > 0, data->timesteps > 0,
                              data->ids, data->lonlats, data->values,
                              data->values[ 0 ] == missingValue,
                              data->values[data->timesteps * data->stations -1]
                                 == missingValue ),
                        AND2( data->subsetDataList,
                              IS_ZERO5( data->stations, data->timesteps,
                                        data->ids, data->lonlats,
                                        data->values ) ) ) );
}



/******************************************************************************
PURPOSE: copyConsolidatedData - Copy subset data to consolidated data arrays.
INPUTS:  Data* data  data->arguments, data->subsetDataList
OUTPUTS: Data* data  data->ids[ stations ]  Station ids.
                     data->lonlats[ stations ][ 2 ]  Station lon-lats.
                     data->values[ timesteps ][ stations ]  Station data.
                     data->values2[timesteps ][ stations ]  Station data2.
******************************************************************************/

static void copyConsolidatedData( Data* data ) {

  PRE010( isValidData( data ), data->subsetDataList,
          data->stations > 0, data->timesteps > 0,
          data->ids, data->lonlats,
          data->values, data->values[ 0 ] == missingValue,
          data->values[ data->timesteps * data->stations - 1 ] == missingValue,
          data->ok );

  const Arguments* const arguments = &data->arguments;
  const int aggregate = arguments->aggregate;
  const int stations  = data->stations;
  int* ids = data->ids;
  float* values = data->values;
  float* values2 = data->values2;
  int timestep = 0;
  int index = 0;
  int count = 0;
  int count1 = 1;
  int timesteps = 0;
  Integer yyyydddhhmm = 0;
  Integer yyyydddhhmm2 = 0;
  UTCTimestamp firstTimestamp = "";
  UTCTimestamp lastTimestamp = "";

  toUTCTimestamp2( arguments->firstTimestamp, firstTimestamp );
  toUTCTimestamp2( arguments->lastTimestamp, lastTimestamp );
  timesteps = hoursInRange( firstTimestamp, lastTimestamp );
  yyyydddhhmm = fromUTCTimestamp( firstTimestamp );
  yyyydddhhmm2 = yyyydddhhmm;

  for ( timestep = 0; timestep < timesteps; ++timestep,
        incrementTimestamp( &yyyydddhhmm ) ) {
    int station = 0;

    for ( station = 0; station < stations; ++station, ++index ) {
      const int id = ids[ station ];
      const Integer yyyymmddhhmmss = toYYYYMMDDHHMMSS( yyyydddhhmm );
      float value2 = 0.0;
      const float value = findStationData( data, yyyymmddhhmmss, id, &value2 );
      CHECK( IN_RANGE( index, 0, data->timesteps * data->stations - 1 ) );

      if ( aggregate == AGGREGATE_NONE ) {
        values[ index ] = value;

        if ( values2 ) {
          values2[ index ] = value2;
        }

      } else if ( value != missingValue ) {

        /* Re-initialize or compute running mean: */

        if ( values[ index ] == missingValue ) {
          values[ index ] = value;
        } else {
          values[ index ] = ( count * values[ index ] + value ) / count1;
        }

        if ( values2 ) {

          if ( values2[ index ] == missingValue ) {
            values2[ index ] = value2;
          } else {
            values2[ index ] = ( count * values2[ index ] + value2 ) / count1;
          }
        }
      }
    }

    if ( aggregate == AGGREGATE_ALL ) {
      index = 0;
      ++count;
      ++count1;
    } else if ( aggregate == AGGREGATE_DAILY ) {
      const Integer ddd  = yyyydddhhmm  / 10000 % 1000;
      const Integer ddd2 = yyyydddhhmm2 / 10000 % 1000;
      const int sameDay = ddd == ddd2;

      if ( sameDay ) {
        index -= stations;
        ++count;
        ++count1;
      } else {
        count = 0;
        count1 = 1;
        yyyydddhhmm2 = yyyydddhhmm;
      }
    }
  }

  FREE_OBJECT( data->subsetDataList );

  POST08( isValidData( data ), data->ok,
          data->subsetDataList == 0,
          data->stations > 0, data->timesteps > 0,
          data->ids, data->lonlats, data->values );
}



/******************************************************************************
PURPOSE: copyUniqueStations - Copy unique station ids and lon-lats subset data.
INPUTS:  Data* data  data->subsetDataList.
OUTPTUS: Data* data  data->stations, data->ids, data->lonlats, data->ok.
******************************************************************************/

static void copyUniqueStations( Data* data ) {

  PRE06( isValidData( data ), data->subsetDataList, data->ok,
         data->stations == 0, data->ids == 0, data->lonlats == 0 );

  /* Sum all subset counts: */

  const VoidList* const subsetDataList = data->subsetDataList;
  const int subsetCount = subsetDataList->count( subsetDataList );
  int counts = 0;
  int subset = 0;
  data->ok = 0;

  for ( subset = 0; subset < subsetCount; ++subset ) {
    const SubsetData* const subsetData =
      subsetDataList->item( subsetDataList, subset );
    counts += subsetData->count;
  }

  CHECK2( counts >= 0, counts >= subsetCount );

  /* Allocate and copy a temporary array to hold all subset ids: */

  {
    int* ids = NEW_ZERO( int, counts );

    if ( ids ) {
      int offset = 0;

      for ( subset = 0; subset < subsetCount; ++subset ) {
        const SubsetData* const subsetData =
          subsetDataList->item( subsetDataList, subset );
        const int count = subsetData->count;
        memcpy( ids + offset, subsetData->ids, count * sizeof *ids );
        offset += count;
      }

      /* Sort temporary array of ids: */

      qsort( ids, counts, sizeof *ids, intCompare );

      /* Scan temporary array of ids and count unique ones: */

      {
        int id = 0;
        int index = 0;
        int uniqueCount = 0;

        for ( index = 0; index < counts; ++index ) {
          const int thisId = ids[ index ];

          if ( thisId != id ) {
            id = thisId;
            ++uniqueCount;
          }
        }

        /* Copy sorted unique ids to data->ids[]: */

        data->ids = NEW_ZERO( int, uniqueCount );

        if ( data->ids ) {
          int writeIndex = 0;
          id = 0;

          for ( index = 0; index < counts; ++index ) {
            const int thisId = ids[ index ];

            if ( thisId != id ) {
              id = thisId;
              CHECK( IN_RANGE( writeIndex, 0, uniqueCount - 1 ) );
              data->ids[ writeIndex ] = id;
              ++writeIndex;
            }
          }

          data->stations = uniqueCount;

          /* Copy sorted unique station lonlats to data->lonlats[]: */

          FREE( ids );
          data->lonlats = NEW_ZERO( float, uniqueCount * 2 );
          data->notes = data->lonlats ? NEW_ZERO( Note, uniqueCount ) : 0;

          if ( data->notes ) {

            for ( index = 0; index < uniqueCount; ++index ) {
              const int lonlatIndex = index + index;
              float longitude = 0.0;
              float latitude = 0.0;
              id = data->ids[ index ];
              findStationLongitudeLatitudeNote( data, id,
                                                &longitude, &latitude,
                                                data->notes[ index ] );
              data->lonlats[ lonlatIndex     ] = longitude;
              data->lonlats[ lonlatIndex + 1 ] = latitude;
            }

            data->ok = 1;
          }
        }
      }

      FREE( ids );
    }
  }

  if ( ! data->ok ) {
    FREE( data->ids );
    FREE( data->lonlats );
    data->stations = 0;
  }

  POST03( isValidData( data ),
          data->subsetDataList,
          IMPLIES_ELSE( data->ok,
                        AND7( data->stations > 0,
                              data->ids,
                              data->ids[ 0 ] > 0,
                              IMPLIES( data->stations > 1,
                                data->ids[ data->stations - 1] > data->ids[0]),
                              data->lonlats,
                              isValidLongitudeLatitude( data->lonlats[ 0 ],
                                                        data->lonlats[ 1 ] ),
                              isValidLongitudeLatitude(
                                data->lonlats[ 2 * data->stations - 2 ],
                                data->lonlats[ 2 * data->stations - 1 ] ) ),
                        IS_ZERO3( data->stations, data->ids, data->lonlats )));
}



/******************************************************************************
PURPOSE: findStationLongitudeLatitudeNote - Find station lon-lat and note
         in subset data.
INPUTS:  const Data* data  data->subsetDataList.
         const int id      Id of station to find.
OUTPTUS: float* longitude  Longitude of station with given id.
         float* latitude   Latitude  of station with given id.
         Note note         Note of station with given id.
******************************************************************************/

static void findStationLongitudeLatitudeNote( const Data* data, const int id,
                                              float* longitude, float* latitude,
                                              Note note ) {

  PRE06( isValidData( data ), data->subsetDataList, id > 0,
         longitude, latitude, note );

  const VoidList* const subsetDataList = data->subsetDataList;
  const int subsetCount = subsetDataList->count( subsetDataList );
  int subset = 0;
  *longitude = *latitude = missingValue;

  for ( subset = 0; subset < subsetCount; ++subset ) {
    const SubsetData* const subsetData =
      subsetDataList->item( subsetDataList, subset );
    const int count = subsetData->count;
    int index = 0;

    for ( index = 0; index < count; ++index ) {
      const int thisId = subsetData->ids[ index ];

      if ( thisId == id ) {
        *longitude = subsetData->longitudes[ index ];
        *latitude  = subsetData->latitudes[  index ];
        strncpy( note, subsetData->notes[ index ],
                 sizeof (Note) / sizeof (char) );
        index = count;
        subset = subsetCount;
      }
    }
  }

  POST0( isValidLongitudeLatitude( *longitude, *latitude ) );
}



/******************************************************************************
PURPOSE: findStationData - Find station timestamped data in subset data.
INPUTS:  const Data* data  data->subsetDataList.
         const Integer timestamp  Timestamp of data to find.
         const int id             Id of station to find.
OUTPUTS: float* value2            2nd component or missingValue if not found.
RETURNS: float Value or missingValue if not found.
******************************************************************************/

static float findStationData( const Data* data, const Integer timestamp,
                              const int id, float* value2 ) {

  PRE05( isValidData( data ), data->subsetDataList, id > 0,
         isValidYYYYMMDDHHMMSS( timestamp ), value2 );

  float result = missingValue;
  const VoidList* const subsetDataList = data->subsetDataList;
  const int subsetCount = subsetDataList->count( subsetDataList );
  int subset = 0;

  for ( subset = 0; subset < subsetCount; ++subset ) {
    const SubsetData* const subsetData =
      subsetDataList->item( subsetDataList, subset );

    if ( subsetData->timestamp == timestamp ) {
      const int count = subsetData->count;
      int index = 0;

      for ( index = 0; index < count; ++index ) {
        const int thisId = subsetData->ids[ index ];

        if ( thisId == id ) {
          result = subsetData->values[ index ];

          if ( subsetData->values2 ) {
            *value2 = subsetData->values2[ index ];
          } else {
            *value2 = missingValue;
          }

          index = count;
          subset = subsetCount;
        }
      }
    }
  }

  POST0( OR2( IN_RANGE( result,
                        data->arguments.minimum, data->arguments.maximum ),
              result == missingValue ) );
  return result;
}



/******************************************************************************
PURPOSE: writeHeader - Write ASCII header of subset to stdout.
INPUTS:  const Data*  data  Data to write
******************************************************************************/

static void writeHeader( const Data* data ) {

  PRE03( isValidData( data ), data->stations, data->ok );
  const int isWind = ! strcmp( data->arguments.variable, "wind" );
  UTCTimestamp firstTimestamp = "";
  toUTCTimestamp2( data->arguments.firstTimestamp, firstTimestamp );
  printf( "SITE 2.0\n" );
  printf( "%s\n", data->arguments.description );
  printf( "%s\n", firstTimestamp );
  printf( "# data dimensions: timesteps stations\n" );
  printf( "%d %d\n", data->timesteps, data->stations );
  printf( "# Variable names:\n" );

  if ( isWind ) {
    printf( "windU windV\n" );
  } else {
    printf( "%s\n", data->arguments.variable );
  }

  printf( "# Variable units:\n" );

  if ( isWind ) {
    printf( "%s %s\n", data->arguments.units, data->arguments.units );
  } else {
    printf( "%s\n", data->arguments.units );
  }

  printf( "# char notes[stations][80] and\n" );
  printf( "# MSB 64-bit integers ids[stations] and\n" );
  printf( "# IEEE-754 64-bit reals sites[stations][2=<longitude,latitude>]" );
  printf( " and\n");
  printf( "# IEEE-754 64-bit reals data[timesteps][stations]:\n");

  POST0( data->ok );
}



/******************************************************************************
PURPOSE: writeXDR - Write XDR format output of subset to stdout.
INPUTS:  INPUTS: Data* data  Data to write
******************************************************************************/

static void writeXDR( Data* data ) {

  PRE02( isValidData( data ), data->ok );

  const int timesteps = data->timesteps;
  const int stations = data->stations;
  const int timestepsTimesStations = timesteps * stations;
  const int stations2 = stations + stations;
  const int count =
    timestepsTimesStations > stations2 ? timestepsTimesStations : stations2;
  double* buffer = NEW_ZERO( double, count );

  if ( buffer ) {
    long long* ibuffer = (long long*) buffer;
    int station = 0;

    writeHeader( data );

    for ( station = 0; station < stations; ++station ) {
      printf( "%-79s\n", data->notes[ station ] );
      ibuffer[ station ] = data->ids[ station ];
    }

    rotate8ByteArrayIfLittleEndian( ibuffer, stations );
    data->ok = fwrite( ibuffer, sizeof *ibuffer, stations, stdout) == stations;

    if ( data->ok ) {

      for ( station = 0; station < stations2; station += 2 ) {
        buffer[ station     ] = data->lonlats[ station ];
        buffer[ station + 1 ] = data->lonlats[ station + 1 ];
      }

      rotate8ByteArrayIfLittleEndian( buffer, stations2 );
      data->ok =
        fwrite( buffer, sizeof *buffer, stations2, stdout ) == stations2;

      if ( data->ok ) {
        int index = 0;

        for ( index = 0; index < timestepsTimesStations; ++index ) {
          buffer[ index ] = data->values[ index ];
        }

        rotate8ByteArrayIfLittleEndian( buffer, timestepsTimesStations );
        data->ok =
          fwrite( buffer, sizeof *buffer, timestepsTimesStations, stdout )
          == timestepsTimesStations;

        if ( data->ok ) {

          if ( data->values2 ) {

            for ( index = 0; index < timestepsTimesStations; ++index ) {
              buffer[ index ] = data->values2[ index ];
            }

            rotate8ByteArrayIfLittleEndian( buffer, timestepsTimesStations );
            data->ok =
              fwrite( buffer, sizeof *buffer, timestepsTimesStations, stdout )
              == timestepsTimesStations;
          }
        }
      }
    }

    FREE( buffer );
  }

  POST0( IS_BOOL( data->ok ) );
}



/******************************************************************************
PURPOSE: writeASCII - Write subset data to stdout in ASCII format
         (spreadsheet, tab-separated text values).
INPUTS:  const Data* data  Data to write.
******************************************************************************/

static void writeASCII( const Data* data ) {

  PRE03( isValidData( data ), data->values, data->ok );
  const int isWind = ! strcmp( data->arguments.variable, "wind" );
  const char* const headerStart =
    "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tSTATION(-)";
  const char* const headerFormat1 = "\t%s(%s)\tNOTE\n";
  const char* const headerFormat2 = "\t%s(%s)\t%s(%s)\tNOTE\n";
  const char* const dataFormat1 =
    "%s\t%10.5f\t%10.5f\t%09d\t%12.5f\t%46s\n";
  const char* const dataFormat2 =
    "%s\t%10.5f\t%10.5f\t%09d\t%12.5f\t%12.5f\t%46s\n";
  const int* ids = data->ids;
  const float* lonlats = data->lonlats;
  const float* values = data->values;
  const float* values2 = data->values2;
  int timesteps    = data->timesteps;
  int stations     = data->stations;
  int timestep = 0;
  Integer timestamp = 0;
  UTCTimestamp timestampString;
  Note note = "";

  toUTCTimestamp2( data->arguments.firstTimestamp, timestampString );
  timestamp = fromUTCTimestamp( timestampString );

  /* Write header row: */

  printf( "%s", headerStart );

  if ( isWind ) {
    printf( headerFormat2,
            "windU", data->arguments.units,
            "windV", data->arguments.units );
  } else {
    printf( headerFormat1, data->arguments.variable, data->arguments.units );
  }

  /* Write data rows: */

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    int station = 0;

    for ( station = 0; station < stations; ++station ) {
      const int station2 = station + station;
      const int id = ids[ station ];
      const float longitude = lonlats[ station2 ];
      const float latitude  = lonlats[ station2 + 1 ];
      const float value  = *values++;
      const float value2 = values2 ? *values2++ : missingValue;
      strncpy( note, data->notes[ station ], 45 );
      note[ 45 ] = '\0';
      CHECK5( id > 0,
              IN_RANGE( longitude,
                        data->arguments.bounds[ LONGITUDE ][ MINIMUM ],
                        data->arguments.bounds[ LONGITUDE ][ MAXIMUM ] ),
              IN_RANGE( latitude,
                        data->arguments.bounds[ LATITUDE ][ MINIMUM ],
                        data->arguments.bounds[ LATITUDE ][ MAXIMUM ] ),
              OR2( value == missingValue,
                   IN_RANGE( value,
                             data->arguments.minimum,
                             data->arguments.maximum ) ),
              OR2( value2 == missingValue,
                   IN_RANGE( value2,
                             data->arguments.minimum,
                             data->arguments.maximum ) ) );

      if ( isWind ) {
        printf( dataFormat2, timestampString, longitude, latitude, id, value,
                value2, note );
      } else {
        printf( dataFormat1, timestampString, longitude, latitude, id, value,
                note );
      }
    }

    incrementTimestamp( &timestamp );
    toUTCTimestamp( timestamp, timestampString );
  }

  POST0( data->ok );
}



/******************************************************************************
PURPOSE: fileTimestamp - YYYYMMDDHHMMSS of file.
INPUTS:  const char* fileName  Name of file to open.
RETURNS: Integer YYYYMMDDHHMMSS of file if successful, else 0 if failed and
         failureMessage() called.
******************************************************************************/

static Integer fileTimestamp( const char* const fileName ) {
  PRE02( fileName, *fileName );
  Integer result = 0;
  const char* const slash = strrchr( fileName, '/' );
  const char* const start = slash ? slash + 1 : fileName;
  const char* const underscore = strchr( start, '_' );

  if ( underscore ) {
    const int yyyymmdd = atoi( start );
    const int hhmm     = atoi( underscore + 1 );
    result = yyyymmdd;
    result *= 10000;
    result += hhmm;
    result *= 100;

    if ( ! isValidYYYYMMDDHHMMSS( result ) ) {
      result = 0;
    }
  }

  POST0( OR2( result == 0, isValidYYYYMMDDHHMMSS( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: generateId - Generate station id from its location.
INPUTS:  const float longitude  Station longitude.
         const float latitude   Station latitude.
RETURNS: int station id.
******************************************************************************/

static int generateId( const float longitude, const float latitude ) {
  PRE0( isValidLongitudeLatitude( longitude, latitude ) );
  const float pLongitude = longitude < 0.0 ? -longitude : longitude;
  const float pLatitude = latitude < 0.0 ? -latitude : latitude;
  const int iLongitude =
    pLongitude >= 100.0 ? (int) ( pLongitude * 100.0 + 0.5 )
    : (int) ( pLongitude * 1000.0 + 0.5 ); /* First 5 sigdig. */
  const int iLatitude  = (int) ( pLatitude * 100.0 + 0.5 ); /* First 4 sigdig*/
  const int result = iLongitude * 10000 + iLatitude; /* 9 sigdig. */
  POST0( IN_RANGE( result, 0, 999999999 ) );
  return result;
}



/* NetCDF convenience wrappers: */



/******************************************************************************
PURPOSE: openNetCDFFile - Open a NetCDFFile for reading.
INPUTS:  const char* fileName Name of file to open.
RETURNS: int NetCDF file ID if successful, else -1 if failed and
         failureMessage() called.
******************************************************************************/

static int openNetCDFFile( const char* const fileName ) {
  PRE02( fileName, *fileName );
  int result = -1;
  const int status = nc_open( fileName, NC_NOWRITE, &result );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't open file '%s' because %s.", fileName, message );
    result = -1;
  } else if ( result < 0 ) {
    nc_close( result );
    failureMessage( "Invalid id for file '%s'.", fileName );
    result = -1;
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: stationsInNetCDFFile - How many stations are in the file?
INPUTS:  const int file  NetCDF file ID.
RETURNS: int > 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

static int stationsInNetCDFFile( const int file ) {
  PRE0( file >= 0 );
  int result = 0;
  int dimensionId = -1;
  int status = nc_inq_dimid( file, "recNum", &dimensionId );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine number of stations because %s.",
                    message );
  } else {
    size_t value = 0;
    status = nc_inq_dimlen( file, dimensionId, &value );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't determine number of stations because %s.",
                      message );
    } else if ( value < 1 ) {
      failureMessage( "Invalid number of stations (%lu).", value );
    } else {
      result = value;
    }
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: readArray - Read 1D array.
INPUTS:  const int file          NetCDF file ID.
         const int count         Number of values in array.
         const int type          NC_INT or NC_FLOAT.
         const char* const name  Name of array to read.
OUTPUTS: void* values            Array of values read.
RETURNS: int 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

static int readArray( const int file, const int count, const int type,
                      const char* const name, void* values ) {
  PRE06( file >= 0, count > 0, IN3( type, NC_INT, NC_FLOAT ), name, *name,
         values );
  int result = 0;
  int variableId = -1;
  int status = nc_inq_varid( file, name, &variableId );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    name, message );
  } else {
    int rank = 0;
    int dimIds[ 32 ];
    nc_type nctype = (nc_type) 0;
    status = nc_inq_var( file, variableId, 0, &nctype, &rank, dimIds, 0 );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't determine info on variable '%s' because %s.",
                       name, message );
    } else if ( OR2( rank != 1, (int) nctype != type ) ) {
      failureMessage( "Mismatched rank/type of variable '%s'.", name );
    } else {
      char dimName[ 32 ] = "";
      size_t size = 0;
      status = nc_inq_dim( file, dimIds[ 0 ], dimName, &size );

      if ( status != NC_NOERR ) {
        const char* const message = nc_strerror( status );
        failureMessage( "Can't determine size of dimension '%s' because %s.",
                       dimName, message );
      } else if ( (int) size != count ) {
        failureMessage( "Mismatched size of variable '%s'.", name );
      } else {
        const size_t starts[ 1 ] = { 0 };
        size_t counts[ 1 ] = { 0 };
        counts[ 0 ] = size;

        if ( nctype == NC_INT ) {
          status = nc_get_vara_int( file, variableId, starts, counts, values );
        } else {
          status = nc_get_vara_float( file, variableId, starts, counts,values);
        }

        if ( status != NC_NOERR ) {
          const char* const message = nc_strerror( status );
          failureMessage( "Can't read variable '%s' because %s.",
                         name, message );
        } else {
          result = 1;
        }
      }
    }
  }

  POST0( result >= 0 );
  return result;
}




/******************************************************************************
PURPOSE: readStringArray - Read 1D array of strings.
INPUTS:  const int file          NetCDF file ID.
         const int count         Number of values in array.
         const int length        Maximum string size.
         const char* const name  Name of array to read.
OUTPUTS: Note notes[]            Array of notes read.
RETURNS: int 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

static int readStringArray( const int file, const int count, const int length,
                            const char* const name, Note notes[] ) {
  PRE07( file >= 0, count > 0, length > 0,
         length < sizeof (Note) / sizeof (char) - 1,
         name, *name, notes );
  int result = 0;
  int variableId = -1;
  int status = nc_inq_varid( file, name, &variableId );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    name, message );
  } else {
    int rank = 0;
    int dimIds[ 32 ];
    nc_type nctype = (nc_type) 0;
    status = nc_inq_var( file, variableId, 0, &nctype, &rank, dimIds, 0 );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't determine info on variable '%s' because %s.",
                       name, message );
    } else if ( OR2( rank != 2, nctype != NC_CHAR ) ) {
      failureMessage( "Mismatched rank/type of variable '%s'.", name );
    } else {
      char dimName[ 32 ] = "";
      size_t size = 0;
      status = nc_inq_dim( file, dimIds[ 0 ], dimName, &size );

      if ( status != NC_NOERR ) {
        const char* const message = nc_strerror( status );
        failureMessage( "Can't determine size of dimension '%s' because %s.",
                       dimName, message );
      } else if ( (int) size != count ) {
        failureMessage( "Mismatched size of variable '%s'.", name );
      } else {
        size_t maximumLength = 0;
        status = nc_inq_dim( file, dimIds[ 1 ], dimName, &maximumLength );

        if ( status != NC_NOERR ) {
          const char* const message = nc_strerror( status );
          failureMessage("Can't determine length of dimension '%s' because %s.",
                         dimName, message );
        } else if ( (int) maximumLength != length ) {
          failureMessage( "Mismatched length of variable '%s'.", name );
        } else {
          const size_t starts[ 2 ] = { 0, 0 };
          size_t counts[ 2 ] = { 0, 0 };
          counts[ 0 ] = size;
          counts[ 1 ] = maximumLength;
          status =
            nc_get_vara_text( file, variableId, starts, counts, notes[0] );

          if ( status != NC_NOERR ) {
            const char* const message = nc_strerror( status );
            failureMessage( "Can't read variable '%s' because %s.",
                           name, message );
          } else {
            result = 1;
          }
        }
      }
    }
  }

  POST0( result >= 0 );
  return result;
}




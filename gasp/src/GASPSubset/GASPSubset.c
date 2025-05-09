
/******************************************************************************
PURPOSE: GASPSubset.c - Read a set of a GASP files, subset the scans to a
         bounds (longitude-latitude rectangle) and write it to
         stdout as XDR (IEEE-754) format binary.

NOTES:   Uses libz.a (see ../../libs/MODIS/zlib) and
         libUtilities.a (../../libs/Utilities).
         Compile 3 different executables - one for each GOES grid:

         gcc -m64 -fopenmp \
             -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -g -DNO_ASSERTIONS \
             -I../../../include -I. -o ../../../bin/$platform/GASPSubset \
             GASPSubset.c -L/usr/local/lib/x86_64 \
             -L../../../lib/$platform -lUtilities.debug -lz -lm

         gcc -m64 -fopenmp -DGOES_13 \
             -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -g -DNO_ASSERTIONS \
             -I../../../include -I. -o ../../../bin/$platform/GASPSubset13 \
             GASPSubset.c -L/usr/local/lib/x86_64 \
             -L../../../lib/$platform -lUtilities.debug -lz -lm

         gcc -m64 -fopenmp -DGOES_13NEW \
             -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -g -DNO_ASSERTIONS \
             -I../../../include -I. -o ../../../bin/$platform/GASPSubset13new \
             GASPSubset.c -L/usr/local/lib/x86_64 \
             -L../../../lib/$platform -lUtilities.debug -lz -lm

         gcc -m64 -fopenmp -DGOES_15 \
             -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -g -DNO_ASSERTIONS \
             -I../../../include -I. -o ../../../bin/$platform/GASPSubset15 \
             GASPSubset.c -L/usr/local/lib/x86_64 \
             -L../../../lib/$platform -lUtilities.debug -lz -lm

HISTORY: 2009-10-26 plessel.todd@epa.gov, Created.
STATUS: unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>     /* For printf(). */
#include <string.h>    /* For strlen(). */
#include <ctype.h>     /* For isdigit(). */

#include <Utilities.h> /* For PRE0*(), NEW_ZERO(), Stream, VoidList, etc. */

/* Z Library routines used (simply prototyped): */

extern void* gzopen( const char* file_name, const char* mode );
extern int gzread( void* file, void* data, unsigned int bytes );
extern int gzclose( void* file );
extern const char* gzerror( void* file, int* unused );

/*================================== TYPES ==================================*/

#if defined( GASP_15 )
enum { ROWS = 962, COLUMNS = 2800 }; /* GOES-15 West [2017-09-06, present]. */
#elif defined( GASP_13NEW )
enum { ROWS = 880, COLUMNS = 2126 }; /* GOES-13 [2012-10-18, 2017-09-05]. */
#elif defined( GASP_13 )
enum { ROWS = 880, COLUMNS = 2128 }; /* GOES-13 [2010-05-03, 2012-10-17]. */
#else
enum { ROWS = 850, COLUMNS = 2000 }; /* GOES-8 [2006-06-01, 2010-05-02]. */
#endif

/* No SCA before YYYYDDDHHMMSS = 20101230000 = 2010-05-03 00:00:00: */

static const Integer timestampWithSCA = 20101230000;

static Real longitudesLatitudes[ 2 ][ ROWS ][ COLUMNS ]; /* Read from file. */

/* For -daily option: */

static Integer counts[ ROWS ][ COLUMNS ];
static Real    means[ ROWS ][ COLUMNS ];

/* For -corners option: */

enum { SW, SE, NW, NE };

static Real longitudeLatitudeCorners[ 2 ][ 4 ][ ROWS ][ COLUMNS ];

/* GASP variables: SCA is only in files after 2010-05-03. */

enum { AOD, MSK, CLS, STD, SFC, CH1, MOS, CLD, SIG, SCA, VARIABLES };

static const char* const variableNames[ VARIABLES ] = {
  "aod", "msk", "cls", "std", "sfc", "ch1", "mos", "cld", "sig", "sca"
};

static const char* const variableRanges[ VARIABLES ] = {
  "-aod_range", "-msk_range", "-cls_range", "-std_range", "-sfc_range",
  "-ch1_range", "-mos_range", "-cld_range", "-sig_range", "-sca_range"
};


/* User-supplied command-line arguments: */

typedef struct {
  const char*  listFile;        /* File listing GASP files to read. */
  const char*  description;     /* User-supplied description. */
  Integer      firstTimestamp;  /* YYYYDDDHHMM of subset. */
  Integer      hours;           /* Number of hours in subset. */
  Integer      daily;           /* Compute daily mean AOD? */
  Integer      corners;         /* Compute grid cell corners? */
  Integer      selected[ VARIABLES ];    /* User-specified output variables. */
  Real         ranges[ VARIABLES ][ 2 ]; /* User-specified filter ranges. */
  Bounds       bounds;          /*bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]*/
} Arguments;

#ifndef NO_ASSERTIONS

/* Arguments invariant: */

static Integer isValidArguments( const Arguments* arguments ) {
  const Integer result =
    AND15( arguments,
           arguments->listFile,
           arguments->listFile[ 0 ],
           arguments->description,
           arguments->description[ 0 ],
           isValidTimestamp( arguments->firstTimestamp ),
           arguments->hours > 0,
           IS_BOOL( arguments->daily ),
           IS_BOOL( arguments->corners ),
           minimumItemI( arguments->selected, VARIABLES ) >= 0,
           maximumItemI( arguments->selected, VARIABLES ) == 1,
           isNanFree( &arguments->ranges[0][0], VARIABLES * 2 ),
           arguments->ranges[0][ MINIMUM ] <= arguments->ranges[0][ MAXIMUM ],
           arguments->ranges[ VARIABLES - 1 ][ MINIMUM ] <=
             arguments->ranges[ VARIABLES - 1 ][ MAXIMUM ],
           isValidBounds( arguments->bounds ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Scan: Result of reading a GASP file. */

typedef struct {
  Integer       timestamp;   /* YYYYDDDHHMM of file containing scan. */
  unsigned char byteData[ VARIABLES ][ ROWS ][ COLUMNS ]; /* Uncompressed. */
  Real          data[     VARIABLES ][ ROWS ][ COLUMNS ]; /* Decoded. */
} Scan;

#ifndef NO_ASSERTIONS

/* Scan invariant: */

static Integer isValidScan( const Scan* scan ) {
  const Integer result =
    AND3( scan,
          isValidTimestamp( scan->timestamp ),
          isNanFree( &scan->data[0][0][0], VARIABLES * ROWS * COLUMNS ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* SubsettedScan: After time/bounds subsetting and variable range filtering. */

typedef struct {
  Integer hasCorners;        /* Are 8 lon/lat corner arrays at end of data? */
  Integer timestamp;         /* YYYYDDDHHMM of file containing scan. */
  Integer variables;         /* longitude, latitude, aod, ... */
  Integer points;            /* Number of subsetted points. */
  Integer indices[ 2 ][ 2 ]; /* indices[ ROW COLUMN ][ MINIMUM MAXIMUM ]. */
  Real*   data;              /* Subsetted data[ variables ][ points ]. */
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
    AND11( subsettedScan,
          IS_BOOL( subsettedScan->hasCorners ),
           isValidTimestamp( subsettedScan->timestamp ),
           IN_RANGE( subsettedScan->variables, 3, 2 + VARIABLES - 1 ),
           IN_RANGE( subsettedScan->points, 1, ROWS * COLUMNS ),
           IN_RANGE( subsettedScan->indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
           IN_RANGE( subsettedScan->indices[ ROW ][ MAXIMUM ],
                     subsettedScan->indices[ ROW ][ MINIMUM ], ROWS - 1 ),
           IN_RANGE( subsettedScan->indices[ COLUMN][MINIMUM], 0, COLUMNS - 1),
           IN_RANGE( subsettedScan->indices[ COLUMN ][ MAXIMUM ],
                     subsettedScan->indices[ COLUMN ][ MINIMUM ], COLUMNS - 1),
           validLongitudesAndLatitudes( subsettedScan->points,
                                        subsettedScan->data,
                                        subsettedScan->data +
                                          subsettedScan->points ),
           isNanFree( subsettedScan->data,
                      subsettedScan->variables * subsettedScan->points ) );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */



/* Data type: */

typedef struct {
  Arguments arguments;         /* User-supplied (command-line) arguments. */
  Scan      scan;              /* Scan of current file being read/processed. */
  VoidList* subsettedScans;    /* List of processes/subsetted scans. */
  Integer   indices[ 2 ][ 2 ]; /* indices[ ROW COLUMN ][ MINIMUM MAXIMUM ]. */
  Integer   ok;                /* Did last command succeed? */
} Data;

/* Data destructor: */

static void deallocateData( Data* data ) {
  PRE0( data );
  FREE_OBJECT( data->subsettedScans ); /* Calls deallocateSubsettedScan(). */
  ZERO_OBJECT( data );
  POST0( data );
}

#ifndef NO_ASSERTIONS

/* Data invariant: */

static Integer isValidData( const Data* data ) {
  Integer result =
    AND10( data,
           isValidArguments( &data->arguments ),
           isValidScan( &data->scan ),
           data->subsettedScans,
           data->subsettedScans->invariant( data->subsettedScans ),
           IN_RANGE( data->indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
           IN_RANGE( data->indices[ ROW ][ MAXIMUM ],
                     data->indices[ ROW ][ MINIMUM ], ROWS - 1 ),
           IN_RANGE( data->indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
           IN_RANGE( data->indices[ COLUMN ][ MAXIMUM ],
                     data->indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ),
           IS_BOOL( data->ok ) );

  if ( result ) {
    const VoidList* const scans = data->subsettedScans;
    const Integer scanCount = scans->count( scans );
    Integer index = 0;

    do {
      const SubsettedScan* const subsettedScan = scans->item( scans, index );
      result = AND6( result,
                     isValidSubsettedScan( subsettedScan ),
                     IN_RANGE( subsettedScan->indices[ ROW ][ MINIMUM ],
                               data->indices[ ROW ][ MINIMUM ],
                               data->indices[ ROW ][ MAXIMUM ] ),
                     IN_RANGE( subsettedScan->indices[ ROW ][ MAXIMUM ],
                               subsettedScan->indices[ ROW ][ MINIMUM ],
                               data->indices[ ROW ][ MAXIMUM ] ),
                     IN_RANGE( subsettedScan->indices[ COLUMN ][ MINIMUM ],
                               data->indices[ COLUMN ][ MINIMUM ],
                               data->indices[ COLUMN ][ MAXIMUM ] ),
                     IN_RANGE( subsettedScan->indices[ COLUMN ][ MAXIMUM ],
                               subsettedScan->indices[ COLUMN ][ MINIMUM ],
                              data->indices[ COLUMN ][ MAXIMUM ] ) );

      ++index;
    } while ( index < scanCount );
  }

  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* ! defined( NO_ASSERTIONS ) */

static Data data; /* Too large to declare as a local variable of main(). */

/*========================== FORWARD DECLARATIONS ===========================*/

static void printUsage( const char* programName );

static Integer parseArguments( Integer argc, char* argv[],
                               Arguments* arguments );

static void initializeArguments( Arguments* arguments );

static Integer parseOptionalArguments( Integer argc, char* argv[],
                                       Integer* arg,
                                       Arguments* arguments );

static Integer parseVariables( Integer argc, char* argv[], Integer* arg,
                               Integer selected[ VARIABLES ] );

static Integer parseRange( Integer argc, char* argv[], Integer* arg,
                           Real range[ 2 ] );

static Integer readLongitudeLatitudeFile( const char* fileName );

static void readData( Data* data );

static void appendDailyMeans( const Integer yyyyddd0000, Data* data );

static Integer readGASPFile( const char* fileName, const Integer indices[2][2],
                             Scan* scan );

static Integer timestampOfFileName( const char* fileName );

static SubsettedScan* subsetScan( const Bounds bounds,
                                  const Integer indices[ 2 ][ 2 ],
                                  const Real ranges[ VARIABLES ][ 2 ],
                                  const Integer selected[ VARIABLES ],
                                  const Integer hasCorners,
                                  Scan* scan );

static Integer subsetScanCount( const Bounds bounds,
                                const Real ranges[ VARIABLES ][ 2 ],
                                Integer indices[ 2 ][ 2 ],
                                Real data[ VARIABLES ][ ROWS ][ COLUMNS ] );

static Integer subsetIndicesByBounds( const Bounds bounds,
                                      Integer* firstRow,
                                      Integer* lastRow,
                                      Integer* firstColumn,
                                      Integer* lastColumn );

static void subsetIndicesByMask( const Real msk[ ROWS ][ COLUMNS ],
                                 Integer* firstRow, Integer* lastRow,
                                 Integer* firstColumn, Integer* lastColumn );

static void copySubsetLongitudesAndLatitudes( const Real msk[ ROWS ][ COLUMNS],
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

static void copySubsetData( const Real data[ VARIABLES ][ ROWS ][ COLUMNS ],
                            const Integer selected[ VARIABLES ],
                            Integer firstRow, Integer lastRow,
                            Integer firstColumn, Integer lastColumn,
                            Real output[] );

static void computeMean( const Bounds bounds,
                         const Real ranges[ VARIABLES ][ 2 ],
                         const Integer indices[ 2 ][ 2 ],
                         Real data[ VARIABLES ][ ROWS ][ COLUMNS ] );

static Integer meanPoints( const Integer indices[ 2 ][ 2 ] );

static void copyMeanData( const Integer indices[ 2 ][ 2 ],
                          const Integer points,
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

static void clampLongitudes( const double longitude,
                             double* nextColumnLongitude,
                             double* nextRowLongitude,
                             double* nextRowNextColumnLongitude );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Read a subset of a GASP file and write it to stdout in
         XDR or ASCII format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Integer ok = 0;

  if ( ! isValidArgs( argc, (const char**) argv ) ) {
    failureMessage( "Invalid command-line arguments." );
    printUsage( argv[ 0 ] );
  } else {
    checkForTest( &argc, argv ); /* Check for and remove any -test arguments.*/
    ZERO_OBJECT( &data );
    data.ok = parseArguments( argc, argv, &data.arguments );

    if ( data.ok ) {
      data.ok =
        subsetIndicesByBounds( (const Real (*)[2]) data.arguments.bounds,
                               &data.indices[ ROW ][ MINIMUM ],
                               &data.indices[ ROW ][ MAXIMUM ],
                               &data.indices[ COLUMN ][ MINIMUM ],
                               &data.indices[ COLUMN ][ MAXIMUM ] );

      if ( data.ok ) {

        if ( data.arguments.corners ) {
          DEBUG( fprintf( stderr, "calling computeCorners()...\n" ); )
          computeCorners( ROWS, COLUMNS,
                          &longitudesLatitudes[ LONGITUDE ][ 0 ][ 0 ],
                          &longitudesLatitudes[ LATITUDE ][ 0 ][ 0 ],
                          &longitudeLatitudeCorners[ 0 ][ 0 ][ 0 ][ 0 ] );
          DEBUG( fprintf( stderr, "finished computeCorners().\n" ); )
        }

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
INPUTS:  const char* programName  Name of program.
******************************************************************************/

static void printUsage( const char* programName ) {
  PRE0( programName );
  fprintf( stderr,
           "\n\n%s - Read a set of GASP files and extract scan\n",
           programName );
  fprintf( stderr, "data subsetted by a lon-lat rectangle and optionally" );
  fprintf( stderr, " filtered by variable ranges.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "  -lonlats <lonlat_file> \\\n" );
  fprintf( stderr, "  -files <listFile> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -timestamp <yyyymmddhh> -hours <count> \\\n" );
  fprintf( stderr, "  [ -daily | " );
  fprintf( stderr, "-variable aod std cld cls sfc ch1 sig mos ] \\\n" );
  fprintf( stderr, "  [ -domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ] " );
  fprintf( stderr, "[ -corners ]\\\n" );
  fprintf( stderr, "  [ -aod_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -std_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -cld_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -cls_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -sfc_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -ch1_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -sig_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -mos_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "  [ -sca_range <minimum> <maximum> ] \\\n");
  fprintf( stderr, "\n\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "Available variables (unitless) are:\n" );
  fprintf( stderr, "  aod: Aerosol Optical Depth [-0.5 2.05]\n" );
  fprintf( stderr, "  std: Standard deviation of aod [0 2.55]\n" );
  fprintf( stderr, "  cld: Cloud flag: 1 is cloudless, 0 is clouded [0 1]\n");
  fprintf( stderr, "  cls: Sum of cld for 25 surrounding pixels [0 25]\n");
  fprintf( stderr, "  sfc: Surface reflectivities [-0.1 0.41]\n" );
  fprintf( stderr, "  ch1: Channel 1 (visible) reflectance [0 0.425]\n" );
  fprintf( stderr, "  sig: Aerosol signal [-0.5 0.52]\n" );
  fprintf( stderr, "  mos: 28-day composite visible image [0 0.425]\n" );
  fprintf( stderr, "  sca: scattering angle of image [0 180]\n" );
  fprintf( stderr, "  -daily computes daily mean of filtered aod\n" );
  fprintf( stderr, "-corners option will output 8 additional variables:\n" );
  fprintf( stderr, "  Longitude_SW Longitude_SE Longitude_NW Longitude_NE\n");
  fprintf( stderr, "  Latitude_SW Latitude_SE Latitude_NW Latitude_NE\n" );
  fprintf( stderr, "that are the linearly interpolated " );
  fprintf( stderr, "(and edge extrapolated)\n" );
  fprintf( stderr, "corner points for each center-pixel point.\n" );
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );
  fprintf( stderr, "Example #1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-lonlats /gasp/data/goes12/lonlats.bin \\\n" );
  fprintf( stderr, "-files /gasp/data/goes12/files.txt \\\n" );
  fprintf( stderr,"-desc http://www.ssd.noaa.gov/PS/FIRE/GASP/gasp.html \\\n");
  fprintf( stderr, "-timestamp 2008062200 -hours 24 \\\n" );
  fprintf( stderr, "-variable aod \\\n" );
  fprintf( stderr, "-domain -76 34 -74 36 > subset.xdr\n\n" );
  fprintf( stderr, "Subset of data for June 22, 2008 near Raleigh, NC, USA\n");
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays\n" );
  fprintf( stderr, "For example:\n" );
  fprintf( stderr, "Swath 2.0\n" );
  fprintf( stderr, "http://www.ssd.noaa.gov/PS/FIRE/GASP/gasp.html" );
  fprintf( stderr, ",GASPSubset\n" );
  fprintf( stderr, "2008-06-22T00:00:00-0000\n" );
  fprintf( stderr, "# Dimensions: variables timesteps scans:\n" );
  fprintf( stderr, "3 24 25\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "Longitude Latitude aod\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "deg deg -\n" );
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
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-lonlats /gasp/data/goes12/lonlats.bin \\\n" );
  fprintf( stderr, "-files /gasp/data/goes12/files.txt \\\n" );
  fprintf( stderr,"-desc http://www.ssd.noaa.gov/PS/FIRE/GASP/gasp.html \\\n");
  fprintf( stderr, "-timestamp 2008062100 -hours 48 \\\n" );
  fprintf( stderr, "-variable aod \\\n" );
  fprintf( stderr, "-aod_range 0 1 \\\n" );
  fprintf( stderr, "-std_range 0 0.3 \\\n" );
  fprintf( stderr, "-cld_range 1 1 \\\n" );
  fprintf( stderr, "-cls_range 15 25 \\\n" );
  fprintf( stderr, "-sfc_range 0.005 0.41 \\\n" );
  fprintf( stderr, "-ch1_range 0.001 0.425 \\\n" );
  fprintf( stderr, "-sig_range 0.01 0.52 \\\n" );
  fprintf( stderr, "-domain -76 34 -74 36 > subset.xdr\n\n" );
  fprintf( stderr, "Like above but includes data filtering ranges.\n" );
  fprintf( stderr, "Example #3:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-lonlats /gasp/data/goes12/lonlats.bin \\\n" );
  fprintf( stderr, "-files /gasp/data/goes12/files.txt \\\n" );
  fprintf( stderr,"-desc http://www.ssd.noaa.gov/PS/FIRE/GASP/gasp.html \\\n");
  fprintf( stderr, "-timestamp 2008062100 -hours 48 \\\n" );
  fprintf( stderr, "-daily \\\n" );
  fprintf( stderr, "-domain -76 34 -74 36 > subset.xdr\n\n" );
  fprintf( stderr, "Computes daily mean of filtered AOD.\n" );
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
    const char* const longitudeLatitudeFile =
      parseArgument2( argc, argv, "-lonlats", &arg );

    if ( longitudeLatitudeFile ) {

      if ( readLongitudeLatitudeFile( longitudeLatitudeFile ) ) {
        arguments->listFile = parseArgument2( argc, argv, "-files", &arg );

        if ( arguments->listFile ) {
          arguments->description = parseArgument2( argc, argv, "-desc", &arg );

          if ( arguments->description ) {

            if ( parseTimestampAndHours( argc, argv, &arg,
                                         &arguments->firstTimestamp,
                                         &arguments->hours ) ) {
              result = parseOptionalArguments( argc, argv, &arg, arguments );
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
PURPOSE: initializeArguments - Initialize arguments.
INPUTS:  Arguments* arguments  Arguments to initialize.
OUTPUTS: Arguments* arguments  Arguments initialized.
******************************************************************************/

static void initializeArguments( Arguments* arguments ) {
  PRE0( arguments );
  ZERO_OBJECT( arguments );
  arguments->ranges[ AOD ][ MINIMUM ] = -0.5;
  arguments->ranges[ AOD ][ MAXIMUM ] = 2.05;
  arguments->ranges[ MSK ][ MINIMUM ] = 1.0;
  arguments->ranges[ MSK ][ MAXIMUM ] = 1.0;
  arguments->ranges[ CLS ][ MINIMUM ] = 0.0;
  arguments->ranges[ CLS ][ MAXIMUM ] = 25.0;
  arguments->ranges[ STD ][ MINIMUM ] = 0.0;
  arguments->ranges[ STD ][ MAXIMUM ] = 2.55;
  arguments->ranges[ SFC ][ MINIMUM ] = -0.1;
  arguments->ranges[ SFC ][ MAXIMUM ] = 0.41;
  arguments->ranges[ CH1 ][ MINIMUM ] = 0.0;
  arguments->ranges[ CH1 ][ MAXIMUM ] = 0.425;
  arguments->ranges[ MOS ][ MINIMUM ] = 0.0;
  arguments->ranges[ MOS ][ MAXIMUM ] = 0.425;
  arguments->ranges[ CLD ][ MINIMUM ] = 0.0;
  arguments->ranges[ CLD ][ MAXIMUM ] = 1.0;
  arguments->ranges[ SIG ][ MINIMUM ] = -0.5;
  arguments->ranges[ SIG ][ MAXIMUM ] = 0.52;
  arguments->ranges[ SCA ][ MINIMUM ] = 0.0;
  arguments->ranges[ SCA ][ MAXIMUM ] = 255.0;
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
  Integer parsedVariable = 0;
  Integer parsedBounds   = 0;
  Integer parsedVariableRange[ VARIABLES ];
  memset( parsedVariableRange, 0, sizeof parsedVariableRange );

  while ( *arg < argc ) {
    const Integer variable =
      indexOfString( argv[ *arg ], variableRanges, VARIABLES );

    if ( AND3( ! IN3( variable, -1, MSK ),
               ! parsedVariableRange[ variable ],
               ! arguments->daily ) ) {
      result = parseRange( argc, argv, arg,
                           &arguments->ranges[ variable ][ 0 ] );
      parsedVariableRange[ variable ] = 1;
    } else if ( AND3( ! strcmp( argv[ *arg ], "-variable" ),
                      ! parsedVariable,
                      ! arguments->daily ) ) {
      parsedVariable = 1;
      result = parseVariables( argc, argv, arg,
                               arguments->selected );
    } else if ( AND3( ! strcmp( argv[ *arg ], "-daily" ),
                      ! parsedVariable,
                      ! arguments->daily ) ) {
      const Integer timestamp = arguments->firstTimestamp;
      *arg += 1; /* Skip "-daily" option. */
      arguments->daily = 1; /* Use hard-coded filter ranges from IDL code*/
      arguments->selected[ AOD ] = 1;

      if ( timestamp < timestampWithSCA ) {
        arguments->ranges[ AOD ][ MINIMUM ] = 0.0;
        arguments->ranges[ AOD ][ MAXIMUM ] = 2.05;
        arguments->ranges[ MSK ][ MINIMUM ] = 1.0;
        arguments->ranges[ MSK ][ MAXIMUM ] = 1.0;
        arguments->ranges[ CLS ][ MINIMUM ] = 15.0;
        arguments->ranges[ CLS ][ MAXIMUM ] = 25.0;
        arguments->ranges[ STD ][ MINIMUM ] = 0.0;
        arguments->ranges[ STD ][ MAXIMUM ] = 0.3 - 1e-6; /* About < 0.3 */
        arguments->ranges[ SFC ][ MINIMUM ] = 0.005 + 1e-6;
        arguments->ranges[ SFC ][ MAXIMUM ] = 0.15 - 1e-6;
        arguments->ranges[ CH1 ][ MINIMUM ] = 0.0 + 1e-6;
        arguments->ranges[ CH1 ][ MAXIMUM ] = 0.425;
        arguments->ranges[ MOS ][ MINIMUM ] = 0.0;
        arguments->ranges[ MOS ][ MAXIMUM ] = 0.425;
        arguments->ranges[ CLD ][ MINIMUM ] = 1.0;
        arguments->ranges[ CLD ][ MAXIMUM ] = 1.0;
        arguments->ranges[ SIG ][ MINIMUM ] = 0.01 + 1e-6;
        arguments->ranges[ SIG ][ MAXIMUM ] = 0.52;
        arguments->ranges[ SCA ][ MINIMUM ] = 0.0;
        arguments->ranges[ SCA ][ MAXIMUM ] = 255.0;
      } else {
        arguments->ranges[ AOD ][ MINIMUM ] = 0.0;
        arguments->ranges[ AOD ][ MAXIMUM ] = 2.05;
        arguments->ranges[ MSK ][ MINIMUM ] = 1.0;
        arguments->ranges[ MSK ][ MAXIMUM ] = 1.0;
        arguments->ranges[ CLS ][ MINIMUM ] = 25.0;
        arguments->ranges[ CLS ][ MAXIMUM ] = 25.0;
        arguments->ranges[ STD ][ MINIMUM ] = 0.0;
        arguments->ranges[ STD ][ MAXIMUM ] = 0.2 - 1e-6; /* About < 0.2 */
        arguments->ranges[ SFC ][ MINIMUM ] = 0.005 + 1e-6;
        arguments->ranges[ SFC ][ MAXIMUM ] = 0.15 - 1e-6;
        arguments->ranges[ CH1 ][ MINIMUM ] = 0.0 + 1e-6;
        arguments->ranges[ CH1 ][ MAXIMUM ] = 0.425;
        arguments->ranges[ MOS ][ MINIMUM ] = 0.0;
        arguments->ranges[ MOS ][ MAXIMUM ] = 0.425;
        arguments->ranges[ CLD ][ MINIMUM ] = 1.0;
        arguments->ranges[ CLD ][ MAXIMUM ] = 1.0;
        arguments->ranges[ SIG ][ MINIMUM ] = 0.01 + 1e-6;
        arguments->ranges[ SIG ][ MAXIMUM ] = 0.52;
        arguments->ranges[ SCA ][ MINIMUM ] = 70.0;
        arguments->ranges[ SCA ][ MAXIMUM ] = 170.0;
      }
    } else if ( AND2( ! strcmp( argv[ *arg ], "-domain" ), ! parsedBounds ) ) {
      parsedBounds = 1;
      result = parseBounds( argc, argv, arg, arguments->bounds );
    } else if ( AND2( ! strcmp( argv[ *arg ], "-corners" ),
                      ! arguments->corners ) ) {
      *arg += 1; /* Skip "-corners" option. */
      arguments->corners = 1;
    } else {
      failureMessage( "Invalid/redundant command-line argument: %s.",
                      argv[ *arg ] );
      result = 0;
      *arg = argc;
    }
  }

  if ( AND3( result, ! parsedVariable, ! arguments->daily ) ) {
    arguments->selected[ AOD ] = 1;
    arguments->selected[ CLS ] = 1;
    arguments->selected[ STD ] = 1;
    arguments->selected[ SFC ] = 1;
    arguments->selected[ CH1 ] = 1;
    arguments->selected[ MOS ] = 1;
    arguments->selected[ CLD ] = 1;
    arguments->selected[ SIG ] = 1;
    arguments->selected[ SCA ] = arguments->firstTimestamp >= timestampWithSCA;
  }

  if ( arguments->firstTimestamp < timestampWithSCA ) {
    arguments->ranges[ SCA ][ MINIMUM ] = 0.0;
    arguments->ranges[ SCA ][ MAXIMUM ] = 255.0;
  }

  POST02( IS_BOOL( result ), IMPLIES( result, isValidArguments( arguments )));
  return result;
}



/******************************************************************************
PURPOSE: parseRange - Parse command-line arguments for -aod_range min max, etc.
INPUTS:  Integer argc     Number of command-line arguments.
         char* argv[]     Command-line argument strings.
         Integer* arg     Index of argument to parse.
OUTPUTS: Integer* arg     Index of next argument to parse.
         Real range[ 2 ]  Initialized MINIMUM MAXIMUM.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer parseRange( Integer argc, char* argv[], Integer* arg,
                           Real range[ 2 ] ) {

  PRE06( isValidArgs( argc, (const char**) argv ),
         argc > 0, arg, IN_RANGE( *arg, 1, argc - 1 ),
         indexOfString( argv[ *arg ], variableRanges, VARIABLES ) >= 0,
         range );
  CHECKING( Integer OLD( arg ) = *arg; )
  Integer result = 0;

  if ( *arg + 2 >= argc ) {
    failureMessage( "Missing parameters to command-line argument %s.",
                    argv[ *arg ] );
  } else {
    ++*arg; /* Skip -aod_range. */
    range[ MINIMUM ] = atoR( argv[ *arg ] );
    ++*arg;
    range[ MAXIMUM ] = atoR( argv[ *arg ] );

    if ( range[ MAXIMUM ] < range[ MINIMUM ] ) {
      failureMessage( "Invalid 2nd (maximum) parameter to command-line "
                      "argument %s.", argv[ *arg ] );
    } else {
      ++*arg;
      result = 1;
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, AND4( ! isNan( range[ MINIMUM ] ),
                                 ! isNan( range[ MAXIMUM ] ),
                                 range[ MINIMUM ] <= range[ MAXIMUM ],
                                 *arg == OLD( arg ) + 3 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseVariables - Parse command-line arguments for -variable.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
         Integer* arg          Index of argument to parse.
OUTPUTS: Integer* arg          Index of next argument to parse.
         Integer selected[ VARIABLES ]  Set to 1 if named, else 0.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer parseVariables( Integer argc, char* argv[], Integer* arg,
                               Integer selected[ VARIABLES ] ) {

  PRE06( isValidArgs( argc, (const char**) argv ),
         argc > 0, arg, IN_RANGE( *arg, 1, argc - 1 ),
         ! strcmp( argv[ *arg ], "-variable" ), selected );
  CHECKING( Integer OLD( arg ) = *arg; )

  Integer result = 0;
  memset( selected, 0, VARIABLES * sizeof (Integer) );

  if ( *arg + 1 >= argc ) {
    failureMessage( "Missing parameter to command-line argument -variables." );
  } else {
    Integer ok = 1;
    ++*arg;

    while ( AND4( ok, *arg < argc, argv[ *arg ][ 0 ],
                  argv[ *arg ][ 0 ] != '-' ) ) {
      const char* const variableName = argv[ *arg ];
      const Integer variable =
        indexOfString( variableName, variableNames, VARIABLES );

      if ( OR3( variable == -1, variable == MSK, selected[ variable ] ) ) {
        failureMessage( "Invalid/redundant variable name %s.", variableName );
        ok = 0;
      } else {
        selected[ variable ] = 1;
        ++*arg;
        result = 1;
      }
    }
  }

  if ( ! result ) {
    memset( selected, 0, VARIABLES * sizeof (Integer) );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result, AND3( minimumItemI( selected, VARIABLES ) >= 0,
                                 maximumItemI( selected, VARIABLES ) == 1,
                                 *arg == OLD( arg ) +
                                         1 + sumI( selected, VARIABLES ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readLongitudeLatitudeFile - Read longitude-latitude file.
INPUTS:  const char* fileName  Name of Longitude-latitude .bin file.
OUTPUTS: Real longitudesLatitudes[ 2 ][ ROWS ][ COLUMNS ].
RETURNS: int 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer readLongitudeLatitudeFile( const char* fileName ) {
  PRE02( fileName, *fileName );
  Integer result  = 0;
  FILE* file = fopen( fileName, "rb" );

  if ( ! file ) {
    failureMessage( "Failed to open lon-lat file %s.", fileName );
  } else {
    const Integer fileSize = 2 * ROWS * COLUMNS * sizeof (float);
    Integer rows    = 0;
    Integer columns = 0;

    if ( fscanf( file, "%*[^\n]\n%*[^\n]\n%lld %lld\n%*[^\n]%*c",
                 &rows, &columns ) != 2 ) {
      failureMessage( "Invalid header in lon-lat file %s.", fileName );
    } else if ( ! AND2( rows == ROWS, columns == COLUMNS ) ) {
      failureMessage( "Unmatched row/column dimensions in lon-lat file %s.",
                      fileName );
    } else if ( fread( longitudesLatitudes, fileSize, 1, file ) != 1 ) {
      failureMessage( "Failed to read coordinates from lon-lat file %s.",
                      fileName );
    } else {
      const Integer count = 2 * rows * columns;
      rotate4ByteArrayIfLittleEndian( longitudesLatitudes, count );
      expand32BitValues( &longitudesLatitudes[0][0][0], count );
      CHECK(validLongitudesAndLatitudes(ROWS * COLUMNS,
                                        &longitudesLatitudes[LONGITUDE][0][0],
                                        &longitudesLatitudes[LATITUDE][0][0]));
      result = 1;
    }

    fclose( file ), file = 0;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readData - Read scan data from GASP files and subset it by time,
         lon-lat and variable ranges.
INPUTS:  Data* data  data->arguments->listfile is the file containing the
         names of GASP files to read and subset by data->arguments->bounds and
         data->arguments->ranges.
OUTPTUS: Data* data  data->subsettedScans is the list of subset data to write.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void readData( Data* data ) {

  PRE05( data, data->ok, isValidArguments( &data->arguments ),
         data->scan.timestamp == 0, data->subsettedScans == 0 );

  Stream* listFile = newFileStream( data->arguments.listFile, "r" );
  Integer yyyyddd = 0;
  data->ok = 0;

  if ( listFile ) {
    const Integer firstTimestamp = data->arguments.firstTimestamp;
    const Integer lastTimestamp =
      offsetTimestamp( firstTimestamp, data->arguments.hours );
    Integer previousTimestamp = 0;

    /* For each GASP file, read it into scan, subset scan and append to list:*/

    do {
      FileName fileName = "";
      listFile->readWord( listFile, fileName,
                          sizeof fileName / sizeof *fileName );

      if ( listFile->ok( listFile ) ) {
        const Integer currentTimestamp = timestampOfFileName( fileName );
        DEBUG( fprintf( stderr, "listing GASP file %s\n", fileName ); )

        if ( ! AND2( currentTimestamp > 0,
                     IMPLIES( previousTimestamp,
                              currentTimestamp > previousTimestamp ) ) ) {
          failureMessage( "Invalid/unordered GASP file %s.", fileName );
        } else if (IN_RANGE(currentTimestamp, firstTimestamp, lastTimestamp)) {

          if ( readGASPFile( fileName, (const Integer (*)[2]) data->indices,
                             &data->scan ) ) {

            if ( data->arguments.daily ) {
              const int scanDay = data->scan.timestamp / 10000;

              if ( yyyyddd == 0 ) {
                memset( counts, 0, sizeof counts );
                memset( means,  0, sizeof means );
                yyyyddd = scanDay;
              } else if ( yyyyddd != scanDay ) {
                appendDailyMeans( yyyyddd * 10000, data );
                memset( counts, 0, sizeof counts );
                memset( means,  0, sizeof means );
                yyyyddd = scanDay;
              }

              /* Update running daily mean with this sub-hourly data: */

              computeMean( (const Real (*)[2]) data->arguments.bounds,
                           (const Real (*)[2]) data->arguments.ranges,
                           (const Integer (*)[2]) data->indices,
                           data->scan.data );
              data->ok = 1;
            } else {
              SubsettedScan* subsettedScan =
                subsetScan( (const Real (*)[2]) data->arguments.bounds,
                            (const Integer (*)[2]) data->indices,
                            (const Real (*)[2]) data->arguments.ranges,
                            data->arguments.selected,
                            data->arguments.corners,
                            &data->scan );

              previousTimestamp = currentTimestamp;

              if ( subsettedScan ) { /* Non-empty subset remains: */

                if ( data->subsettedScans == 0 ) { /* Create list if needed: */
                  data->subsettedScans =
                    newVoidList( deallocateSubsettedScan, 0 );
                }

                if ( data->subsettedScans) { /*Append subsetted scan to list:*/
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

      listFile->readString( listFile, fileName, 2 ); /* Read '\n'. */
    } while ( ! listFile->isAtEnd( listFile ) );

    FREE_OBJECT( listFile );
  }

  if ( data->ok ) {

    if ( data->arguments.daily ) {
      CHECK( yyyyddd );

      if ( data->subsettedScans == 0 ) {
        appendDailyMeans( yyyyddd * 10000, data );
      } else {
        const SubsettedScan* const lastSubsettedScan =
          data->subsettedScans->item( data->subsettedScans, LAST_ITEM );
        const Integer lastScanDay = lastSubsettedScan->timestamp / 10000;

        if ( yyyyddd != lastScanDay ) {
          appendDailyMeans( yyyyddd * 10000, data );
          data->ok =
            AND2( data->subsettedScans,
                  data->subsettedScans->count( data->subsettedScans ) > 0 );
        }
      }

    } else if ( data->subsettedScans == 0 ) {
      failureMessage( "No scans were in the subset." );
      data->ok = 0;
    }
  }

  POST02( isValidArguments( &data->arguments ),
          IMPLIES( data->ok, isValidData( data ) ) );
}



/******************************************************************************
PURPOSE: appendDailyMeans - Append a SubsettedScan of the daily means.
INPUTS:  const integer yyyyddd0000  Timestamp of first scan.
         Data* data                 Data.
OUTPTUS: Data* data  data->subsettedScans the list of subset data to write.
NOTES:   If successful then data->ok == 1, else data->ok == 0 and
         failureMessage() is called.
******************************************************************************/

static void appendDailyMeans( const Integer yyyyddd0000, Data* data ) {

  PRE04( isValidTimestamp( yyyyddd0000 ),
         data, data->ok, isValidArguments( &data->arguments ) );

  const Integer points = meanPoints((const Integer (*)[2]) data->indices);
  CHECKING( const Integer OLD( count ) =
    data->subsettedScans ? data->subsettedScans->count( data->subsettedScans )
    : 0 );
  data->ok = 0;

  if ( points > 0 )  {
    SubsettedScan* subsettedScan = NEW_ZERO( SubsettedScan, 1 );

    if ( subsettedScan ) {
      const Integer corners = data->arguments.corners;
      subsettedScan->hasCorners = corners;
      subsettedScan->timestamp = yyyyddd0000;
      subsettedScan->variables = 3; /* Longitude, latitude, AOD. */
      memcpy( &subsettedScan->indices[0][0], &data->indices[0][0],
              sizeof data->indices );
      subsettedScan->points = points;
      subsettedScan->data =
        NEW_ZERO( Real,
                  ( subsettedScan->variables + subsettedScan->hasCorners * 8 )
                  * points );

      if ( subsettedScan->data ) {
        Real* const subsetLongitudes   = subsettedScan->data;
        Real* const subsetLatitudes    = subsetLongitudes + points;
        Real* const subsetData         = subsetLatitudes  + points;
        Real* const subsetLongitudesSW =
          corners ? subsettedScan->data + subsettedScan->variables * points :0;
        Real* const subsetLongitudesSE = corners ? subsetData + points * 2 : 0;
        Real* const subsetLongitudesNW = corners ? subsetData + points * 3 : 0;
        Real* const subsetLongitudesNE = corners ? subsetData + points * 4 : 0;
        Real* const subsetLatitudesSW  = corners ? subsetData + points * 5 : 0;
        Real* const subsetLatitudesSE  = corners ? subsetData + points * 6 : 0;
        Real* const subsetLatitudesNW  = corners ? subsetData + points * 7 : 0;
        Real* const subsetLatitudesNE  = corners ? subsetData + points * 8 : 0;
        copyMeanData( (const Integer (*)[2]) data->indices, points,
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
PURPOSE: readGASPFile - Read and decode a subset of scan data from GASP file.
INPUTS:  const char* fileName  Name of compressed GASP file to read.
         const indices[ 2 ][ 2 ]  Subset indices[ROW COLUMN][MINIMUM MAXIMUM].
OUTPTUS: Scan* scan  scan->timestamp  YYYYDDDHHMM of file
                     scan->byteData[ VARIABLES ][ ROWS ][ COLUMNS ]
                       uncompressed byte data
                     scan->data[ VARIABLES ][ ROWS ][ COLUMNS ] decoded data.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer readGASPFile( const char* fileName, const Integer indices[2][2],
                             Scan* scan ) {

  PRE07( fileName,
         indices,
         IN_RANGE( indices[ ROW    ][ MINIMUM ], 0, ROWS - 1 ),
         IN_RANGE( indices[ ROW    ][ MAXIMUM ],
                   indices[ ROW    ][ MINIMUM ], ROWS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                   indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ),
         scan );

  Integer result = 0;
  void* inputFile = gzopen( fileName, "rb" );
  int unused = 0;

  DEBUG( fprintf( stderr, "Reading GASP file %s [%lld %lld] [%lld %lld]\n",
                  fileName,
                  indices[ ROW ][ MINIMUM ], indices[ ROW ][ MAXIMUM ],
                  indices[ COLUMN ][ MINIMUM ], indices[ COLUMN ][ MAXIMUM ]);)

  if ( ! inputFile ) {
    failureMessage( "Failed to open GASP file %s for reading because %s.",
                    fileName, gzerror( inputFile, &unused ) );
  } else {
    const Integer yyyydddhhmmss = timestampOfFileName( fileName );
    Integer ok = 0;

    if ( isValidTimestamp( yyyydddhhmmss ) ) {
      const Integer sizeOfVariable = sizeof scan->byteData / VARIABLES;
      const Integer bytesToRead =
        sizeof scan->byteData -
        sizeOfVariable * ( yyyydddhhmmss < timestampWithSCA );
      const Integer bytesRead = /* The most time-consuming routine: gzread().*/
        gzread( inputFile, &scan->byteData[0][0][0], bytesToRead );
      ok = bytesRead == bytesToRead;

      if ( ! ok ) {
        failureMessage( "Failed to read %lld bytes from GASP file %s"
                        "because %s.",
                        bytesToRead, fileName,
                        gzerror( inputFile, &unused ) );
      }
    }

    if ( ok ) { /* Decode data: */
      const Real one600th = 1.0 / 600.0;
      const Integer firstRow    = indices[ ROW    ][ MINIMUM ];
      const Integer lastRow     = indices[ ROW    ][ MAXIMUM ];
      const Integer firstColumn = indices[ COLUMN ][ MINIMUM ];
      const Integer lastColumn  = indices[ COLUMN ][ MAXIMUM ];
      Integer row = 0;

      /* Decode subset of data in parallel: */

#pragma omp parallel for

      for ( row = firstRow; row <= lastRow; ++row ) {
        Integer column = 0;

        for ( column = firstColumn; column <= lastColumn; ++column ) {
          scan->data[ AOD ][ row ][ column ] =
            scan->byteData[ AOD ][ row ][ column ] * 0.01 - 0.5;
          scan->data[ MSK ][ row ][ column ] =
            scan->byteData[ MSK ][ row ][ column ];
          scan->data[ CLS ][ row ][ column ] =
            scan->byteData[ CLS ][ row ][ column ];
          scan->data[ STD ][ row ][ column ] =
            scan->byteData[ STD ][ row ][ column ] * 0.01;
          scan->data[ SFC ][ row ][ column ] =
            scan->byteData[ SFC ][ row ][ column ] * 0.002 - 0.1;
          scan->data[ CH1 ][ row ][ column ] =
            scan->byteData[ CH1 ][ row ][ column ] * one600th;
          scan->data[ MOS ][ row ][ column ] =
            scan->byteData[ MOS ][ row ][ column ] * one600th;
          scan->data[ CLD ][ row ][ column ] =
            scan->byteData[ CLD ][ row ][ column ];
          scan->data[ SIG ][ row ][ column ] =
            scan->byteData[ SIG ][ row ][ column ] * 0.004 - 0.5;
          scan->data[ SCA ][ row ][ column ] =
            scan->byteData[ SCA ][ row ][ column ];

          DEBUG( if ( row - firstRow < 3 && column - firstColumn < 3 )
                   fprintf( stderr,
                            "Raw subset data:\n"
                            " AOD: %u MSK: %u CLS: %u STD: %u"
                            " SFC: %u CH1: %u MOS: %u CLD: %u"
                            " SIG: %u SCA: %u\n"
                            "Decoded subset data:\n"
                            " AOD: %lf MSK: %lf CLS: %lf STD: %lf"
                            " SFC: %lf CH1: %lf MOS: %lf CLD: %lf"
                            " SIG: %lf SCA: %lf\n",
                            scan->byteData[ AOD ][ row ][ column ],
                            scan->byteData[ MSK ][ row ][ column ],
                            scan->byteData[ CLS ][ row ][ column ],
                            scan->byteData[ STD ][ row ][ column ],
                            scan->byteData[ SFC ][ row ][ column ],
                            scan->byteData[ CH1 ][ row ][ column ],
                            scan->byteData[ MOS ][ row ][ column ],
                            scan->byteData[ CLD ][ row ][ column ],
                            scan->byteData[ SIG ][ row ][ column ],
                            scan->byteData[ SCA ][ row ][ column ],
                            scan->data[ AOD ][ row ][ column ],
                            scan->data[ MSK ][ row ][ column ],
                            scan->data[ CLS ][ row ][ column ],
                            scan->data[ STD ][ row ][ column ],
                            scan->data[ SFC ][ row ][ column ],
                            scan->data[ CH1 ][ row ][ column ],
                            scan->data[ MOS ][ row ][ column ],
                            scan->data[ CLD ][ row ][ column ],
                            scan->data[ SIG ][ row ][ column ],
                            scan->data[ SCA ][ row ][ column ] ); )
        }
      }

      scan->timestamp = timestampOfFileName( fileName );
      result = scan->timestamp > 0;
    }

    gzclose( inputFile ), inputFile = 0;
  }

  DEBUG( fprintf( stderr, "%lld\n", result ); )
  POST02( IS_BOOL( result ), IMPLIES( result, isValidScan( scan ) ) );
  return result;
}



/******************************************************************************
PURPOSE: timestampOfFileName - Timestamp of GASP file name.
INPUTS:  const char* fileName  Name of GASP file to parse. E.g.,
                               "testdata/2008174231515_i16_US.all.aod.gz"
RETURNS: Integer YYYYDDDHHMM. E.g., 20081742315, or 0 if failed and
         failureMessage() is called.
******************************************************************************/

static Integer timestampOfFileName( const char* fileName ) {
  PRE0( fileName );
  Integer result = 0;
  const char* name = strrchr( fileName, '/' );

  /* Skip directory, if present: */

  if ( name ) {
    ++name;
  } else {
    name = fileName;
  }

  /* Parse timestamp: */

  if ( AND3( name, isdigit( *name ), strlen( name ) > 11 ) ) {
    char timestamp[ 11 + 1 ] = "";
    memset( timestamp, 0, sizeof timestamp );
    strncpy( timestamp, name, 11 ); /* E.g., "20081742315". */
    result = atoI( timestamp );

    if ( ! isValidTimestamp( result ) ) {
      failureMessage( "Invalid timestamp %s.", timestamp );
      result = 0;
    }
  }

  POST0( IMPLIES( result, isValidTimestamp( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: subsetScan - Filter and subset scan by bounds and variable ranges.
INPUTS:  const Bounds bounds          Subset longitude-latitude bounds.
         const Integer indices[2][2]  Subset indices[ ROW COLUMN][MIN/MAXIMUM].
         const Real ranges[VARIABLES][2]     Data ranges to filter by.
         const Integer selected[VARAIBLES]   Flags: 1=output, 0=skip.
         const Integer corners               Include corners in data?
         Scan* scan                          Uncompressed/decoded scan.
OUTPUTS: Scan* scan   scan->data[ MSK ][ ROWS ][COLUMNS] 0'd by filter.
RETURNS: SubsettedScan* if successful, else 0 if no points are in the subset or
         there was a memory allocation failure and failureMessage() is called.
NOTES:   Must use FREE() on returned result when finished.
******************************************************************************/

static SubsettedScan* subsetScan( const Bounds bounds,
                                  const Integer indices[ 2 ][ 2 ],
                                  const Real ranges[ VARIABLES ][ 2 ],
                                  const Integer selected[ VARIABLES ],
                                  const Integer corners,
                                  Scan* scan ) {

  PRE012( isValidBounds( bounds ),
          isValidScan( scan ),
          IN_RANGE( indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
          IN_RANGE( indices[ ROW ][ MAXIMUM ],
                    indices[ ROW ][ MINIMUM ], ROWS - 1 ),
          IN_RANGE( indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
          IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                    indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ),
          isNanFree( &ranges[0][0], VARIABLES * 2 ),
          selected,
          minimumItemI( selected, VARIABLES ) >= 0,
          maximumItemI( selected, VARIABLES ) == 1,
          sumI( selected, VARIABLES ) > 0,
          IS_BOOL( corners ) );

  SubsettedScan* result = 0;
  Integer filteredIndices[ 2 ][ 2 ] = { { 0, 0 }, { 0, 0 } };
  const void* const unused_ =
    memcpy( filteredIndices, indices, sizeof filteredIndices );
  const Integer points =
    subsetScanCount( bounds, (const Real (*)[2]) ranges,
                     filteredIndices, scan->data );

  if ( points ) {
    result = NEW_ZERO( SubsettedScan, 1 );

    if ( result ) {
      const Integer variables = 2 + sumI( selected, VARIABLES );/*lon,lat,aod*/
      const Integer subsetCount = ( variables + corners * 8 ) * points;
      CHECK3( variables >= 3, points > 0, subsetCount >= 3 );
      result->hasCorners = corners;
      result->timestamp = scan->timestamp;
      result->variables = variables;
      result->points    = points;
      memcpy( result->indices, filteredIndices, sizeof result->indices );
      result->data = NEW_ZERO( Real, subsetCount );

      if ( result->data ) {
        Real* const subsetLongitudes   = result->data;
        Real* const subsetLatitudes    = subsetLongitudes + points;
        Real* const subsetData         = subsetLatitudes  + points;
        Real* const subsetLongitudesSW =
          corners ? result->data + variables * points : 0;
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
        const Integer firstRow         = filteredIndices[ ROW    ][ MINIMUM ];
        const Integer lastRow          = filteredIndices[ ROW    ][ MAXIMUM ];
        const Integer firstColumn      = filteredIndices[ COLUMN ][ MINIMUM ];
        const Integer lastColumn       = filteredIndices[ COLUMN ][ MAXIMUM ];

        CHECK( IN_RANGE( ( 1 + lastRow - firstRow ) *
                         ( 1 + lastColumn - firstColumn ),
                          points, ROWS * COLUMNS ) );

        copySubsetLongitudesAndLatitudes( (const Real (*)[ COLUMNS ])
                                          scan->data[ MSK ],
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

        copySubsetData( (const Real (*)[ROWS][COLUMNS]) scan->data,
                        selected, firstRow, lastRow, firstColumn, lastColumn,
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
         Integer indices[ 2 ][ 2 ]            indices[ROW COLUMN][MIN/MAXIMUM].
         Real data[ VARIABLES ][ ROWS ][ COLUMNS ] Uncompressed & decoded scan.
OUTPUTS: Integer indices[ 2 ][ 2 ]  Subset indices reduced by computed mask.
         Real data[ MSK ][ ROWS ][COLUMNS] Zereoed by filter.
RETURNS: Integer number of points remaining after subsetting/filtering.
******************************************************************************/

static Integer subsetScanCount( const Bounds bounds,
                                const Real ranges[ VARIABLES ][ 2 ],
                                Integer indices[ 2 ][ 2 ],
                                Real data[ VARIABLES ][ ROWS ][ COLUMNS ] ) {

  PRE08( isValidBounds( bounds ),
         isNanFree( &data[0][0][0], VARIABLES * ROWS * COLUMNS ),
         indices,
         IN_RANGE( indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
         IN_RANGE( indices[ ROW ][ MAXIMUM ],
                   indices[ ROW ][ MINIMUM ], ROWS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                   indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ),
         isNanFree( &ranges[0][0], VARIABLES * 2 ) );

  const Integer firstRow    = indices[ ROW    ][ MINIMUM ];
  const Integer lastRow     = indices[ ROW    ][ MAXIMUM ];
  const Integer firstColumn = indices[ COLUMN ][ MINIMUM ];
  const Integer lastColumn  = indices[ COLUMN ][ MAXIMUM ];
  const Real longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const Real longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const Real latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const Real latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const Real aodMinimum = ranges[ AOD ][ MINIMUM ];
  const Real aodMaximum = ranges[ AOD ][ MAXIMUM ];
  const Real clsMinimum = ranges[ CLS ][ MINIMUM ];
  const Real clsMaximum = ranges[ CLS ][ MAXIMUM ];
  const Real stdMinimum = ranges[ STD ][ MINIMUM ];
  const Real stdMaximum = ranges[ STD ][ MAXIMUM ];
  const Real sfcMinimum = ranges[ SFC ][ MINIMUM ];
  const Real sfcMaximum = ranges[ SFC ][ MAXIMUM ];
  const Real ch1Minimum = ranges[ CH1 ][ MINIMUM ];
  const Real ch1Maximum = ranges[ CH1 ][ MAXIMUM ];
  const Real mosMinimum = ranges[ MOS ][ MINIMUM ];
  const Real mosMaximum = ranges[ MOS ][ MAXIMUM ];
  const Real cldMinimum = ranges[ CLD ][ MINIMUM ];
  const Real cldMaximum = ranges[ CLD ][ MAXIMUM ];
  const Real sigMinimum = ranges[ SIG ][ MINIMUM ];
  const Real sigMaximum = ranges[ SIG ][ MAXIMUM ];
  const Real scaMinimum = ranges[ SCA ][ MINIMUM ];
  const Real scaMaximum = ranges[ SCA ][ MAXIMUM ];
  Integer result = 0;
  Integer row = 0;

  /* Compute reduced mask in parallel: */

#pragma omp parallel for reduction( + : result )

  for ( row = firstRow; row <= lastRow; ++row ) {
    Integer column = 0;

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Real longitude = longitudesLatitudes[ LONGITUDE ][ row ][ column ];
      const Real latitude  = longitudesLatitudes[ LATITUDE  ][ row ][ column ];
      const Real aod = data[ AOD ][ row ][ column ];
      const Real msk = data[ MSK ][ row ][ column ];
      const Real cls = data[ CLS ][ row ][ column ];
      const Real std = data[ STD ][ row ][ column ];
      const Real sfc = data[ SFC ][ row ][ column ];
      const Real ch1 = data[ CH1 ][ row ][ column ];
      const Real mos = data[ MOS ][ row ][ column ];
      const Real cld = data[ CLD ][ row ][ column ];
      const Real sig = data[ SIG ][ row ][ column ];
      const Real sca = data[ SCA ][ row ][ column ];
      const Integer output =
        AND12( msk == 1.0,
               IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
               IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ),
               IN_RANGE( aod, aodMinimum, aodMaximum ),
               IN_RANGE( cls, clsMinimum, clsMaximum ),
               IN_RANGE( std, stdMinimum, stdMaximum ),
               IN_RANGE( sfc, sfcMinimum, sfcMaximum ),
               IN_RANGE( ch1, ch1Minimum, ch1Maximum ),
               IN_RANGE( mos, mosMinimum, mosMaximum ),
               IN_RANGE( cld, cldMinimum, cldMaximum ),
               IN_RANGE( sig, sigMinimum, sigMaximum ),
               IN_RANGE( sca, scaMinimum, scaMaximum ) );
      data[ MSK ][ row ][ column ] = output; /* Zero msk for next pass.*/
      result += output;
    }
  }

  /* Compute reduced indices sequentially: */

  if ( result ) {
    subsetIndicesByMask( (const Real (*)[ COLUMNS ]) data[ MSK ],
                         &indices[ ROW ][ MINIMUM ],
                         &indices[ ROW ][ MAXIMUM ],
                         &indices[ COLUMN ][ MINIMUM ],
                         &indices[ COLUMN ][ MAXIMUM ] );

    CHECK4( IN_RANGE( indices[ ROW ][ MINIMUM ], firstRow, lastRow ),
            IN_RANGE( indices[ ROW ][ MAXIMUM ],
                      indices[ ROW ][ MINIMUM ], lastRow ),
            IN_RANGE( indices[ COLUMN ][ MINIMUM ], firstColumn, lastColumn ),
            IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                      indices[ COLUMN ][ MINIMUM ], lastColumn ) );
  }

  POST08( result >= 0,
          isNanFree( &data[0][0][0], VARIABLES * ROWS * COLUMNS ),
          IN_RANGE( indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
          IN_RANGE( indices[ ROW ][ MAXIMUM ],
                    indices[ ROW ][ MINIMUM ], ROWS - 1 ),
          IN_RANGE( indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
          IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                    indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ),
          minimumItem( &data[ MSK ][0][0], ROWS * COLUMNS ) >= 0.0,
          maximumItem( &data[ MSK ][0][0], ROWS * COLUMNS ) <= 1.0 );

  return result;
}



/******************************************************************************
PURPOSE: subsetIndicesByBounds - Subset row and column indices by bounds.
INPUTS:  const Bounds bounds     longitude-latitude bounds of subset.
OUTPUTS: Integer*   firstRow     Index of first row of subset.
         Integer*   lastRow      Index of last row of subset.
         Integer*   firstColumn  Index of first column of subset.
         Integer*   lastColumn   Index of last column of subset.
RETURNS: Integer 1 if a non-empty subset exists, else 0 and outputs are -1.
NOTES:   Uses global array longitudesLatitudes.
******************************************************************************/

static Integer subsetIndicesByBounds( const Bounds bounds,
                                      Integer* firstRow,
                                      Integer* lastRow,
                                      Integer* firstColumn,
                                      Integer* lastColumn ) {

  PRE07( bounds, isValidBounds( bounds ),
         firstRow, lastRow, firstColumn, lastColumn,
         validLongitudesAndLatitudes( ROWS * COLUMNS,
           &longitudesLatitudes[ LONGITUDE ][ 0 ][ 0 ],
           &longitudesLatitudes[ LATITUDE  ][ 0 ][ 0 ] ) );

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

  for ( row = 0; row < ROWS; ++row ) {

    for ( column = 0; column < COLUMNS; ++column ) {
      const Real longitude = longitudesLatitudes[ LONGITUDE ][ row ][ column ];
      const Real latitude  = longitudesLatitudes[ LATITUDE  ][ row ][ column ];
      const Integer inside =
        AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
              IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

      if ( inside ) {
        *firstRow = theFirstRow = row;
        result = 1;
        row = ROWS;
        column = COLUMNS;
      }
    }
  }

  if ( result ) {

    /* Loop backward through rows to find last subset row: */

    for ( row = ROWS - 1; row >= theFirstRow; --row ) {

      for ( column = 0; column < COLUMNS; ++column ) {
        const Real longitude = longitudesLatitudes[ LONGITUDE ][ row ][column];
        const Real latitude  = longitudesLatitudes[ LATITUDE  ][ row ][column];
        const Integer inside =
          AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

        if ( inside ) {
          *lastRow = theLastRow = row;
          row = 0;
          column = COLUMNS;
        }
      }
    }

    CHECK2( IN_RANGE( *firstRow, 0, ROWS - 1 ),
            IN_RANGE( *lastRow, *firstRow, ROWS - 1 ) );


    /* Loop forward through columns to find first subset column: */

    for ( column = 0; column < COLUMNS; ++column ) {

      for ( row = theFirstRow; row <= theLastRow; ++row ) {
        const Real longitude = longitudesLatitudes[LONGITUDE][ row ][ column ];
        const Real latitude  = longitudesLatitudes[LATITUDE ][ row ][ column ];
        const Integer inside =
          AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

        if ( inside ) {
          *firstColumn = theFirstColumn = column;
          row = ROWS;
          column = COLUMNS;
        }
      }
    }

    /* Loop backward through columns to find last subset column: */

    for ( column = COLUMNS - 1; column >= theFirstColumn; --column ) {

      for ( row = theFirstRow; row <= theLastRow; ++row ) {
        const Real longitude = longitudesLatitudes[LONGITUDE][ row ][ column ];
        const Real latitude  = longitudesLatitudes[LATITUDE ][ row ][ column ];
        const Integer inside =
          AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
                IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ) );

        if ( inside ) {
          *lastColumn = column;
          row = ROWS;
          column = 0;
        }
      }
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                       AND4( IN_RANGE( *firstRow, 0, ROWS - 1 ),
                             IN_RANGE( *lastRow, *firstRow, ROWS - 1 ),
                             IN_RANGE( *firstColumn, 0, COLUMNS - 1 ),
                             IN_RANGE( *lastColumn, *firstColumn, COLUMNS -1)),
                       AND4( *firstRow == -1, *lastRow == -1,
                             *firstColumn == -1, *lastColumn == -1 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: subsetIndicesByMask - Subset row and column indices by mask flag.
INPUTS:  const Real msk[ ROWS ][ COLUMNS ]  Subset-modified mask flag.
         Integer*   firstRow     Index of first row of bounds subset.
         Integer*   lastRow      Index of alst row of bounds subset.
         Integer*   firstColumn  Index of first column of bounds subset.
         Integer*   lastColumn   Index of last column of bounds subset.
OUTPUTS: Integer*   firstRow     Index of first row of mask-filtered subset.
         Integer*   lastRow      Index of last row of mask-filtered subset.
         Integer*   firstColumn  Index of first column of mask-filtered subset.
         Integer*   lastColumn   Index of last column of mask-filtered subset.
******************************************************************************/

static void subsetIndicesByMask( const Real msk[ ROWS ][ COLUMNS ],
                                 Integer* firstRow, Integer* lastRow,
                                 Integer* firstColumn, Integer* lastColumn) {

  PRE011( msk,
          minimumItem( &msk[0][0], ROWS * COLUMNS ) >= 0.0,
          maximumItem( &msk[0][0], ROWS * COLUMNS ) <= 1.0,
          firstRow,
          lastRow,
          firstColumn,
          lastColumn,
          IN_RANGE( *firstRow, 0, ROWS - 1 ),
          IN_RANGE( *lastRow, *firstRow, ROWS - 1 ),
          IN_RANGE( *firstColumn, 0, COLUMNS - 1 ),
          IN_RANGE( *lastColumn, *firstColumn, COLUMNS - 1 ) );

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

    for ( column = theFirstColumn; column <= theLastColumn; ++column ) {
      const Real mask = msk[ row ][ column ];

      if ( mask ) {
        theFirstRow = row;
        row = ROWS;
        column = COLUMNS;
      }
    }
  }

  /* Loop backward through rows to find last subset row: */

  for ( row = theLastRow; row >= theFirstRow; --row ) {

    for ( column = theFirstColumn; column <= theLastColumn; ++column ) {
      const Real mask = msk[ row ][ column ];

      if ( mask ) {
        theLastRow = row;
        row = 0;
        column = COLUMNS;
      }
    }
  }

  CHECK2( IN_RANGE( theFirstRow, *firstRow, *lastRow ),
          IN_RANGE( theLastRow, theFirstRow, *lastRow ) );

  /* Loop forward through columns to find first subset column: */

  for ( column = theFirstColumn; column <= theLastColumn; ++column ) {

    for ( row = theFirstRow; row <= theLastRow; ++row ) {
      const Real mask = msk[ row ][ column ];

      if ( mask ) {
        theFirstColumn = column;
        row = ROWS;
        column = COLUMNS;
      }
    }
  }

  /* Loop backward through columns to find last subset column: */

  for ( column = theLastColumn; column >= theFirstColumn; --column ) {

    for ( row = theFirstRow; row <= theLastRow; ++row ) {
      const Real mask = msk[ row ][ column ];

      if ( mask ) {
        theLastColumn = column;
        row = ROWS;
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
INPUTS:  const Real msk[ ROWS ][ COLUMNS ]  Subset-modified mask flag.
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

static void copySubsetLongitudesAndLatitudes( const Real msk[ ROWS ][ COLUMNS],
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

  PRE010( minimumItem( &msk[0][0], ROWS * COLUMNS ) >= 0.0,
          maximumItem( &msk[0][0], ROWS * COLUMNS ) <= 1.0,
          IN_RANGE( points, 1,
                    (1 + lastRow - firstRow) * (1 + lastColumn - firstColumn)),
          IN_RANGE( firstRow, 0, ROWS - 1 ),
          IN_RANGE( lastRow, firstRow, ROWS - 1 ),
          IN_RANGE( firstColumn, 0, COLUMNS - 1 ),
          IN_RANGE( lastColumn, firstColumn, COLUMNS - 1 ),
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
  Integer index = 0;

  for ( row = firstRow; row <= lastRow; ++row ) {
    Integer column = 0;

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Real mask = msk[ row ][ column ];

      if ( mask ) {
        const Real longitude = longitudesLatitudes[ LONGITUDE][ row ][ column];
        const Real latitude  = longitudesLatitudes[ LATITUDE ][ row ][ column];
        CHECKING( ++count; )
        CHECK3( IN_RANGE( longitude, -180.0, 180.0 ),
                IN_RANGE( latitude, -90.0, 90.0 ),
                IN_RANGE( count, 1, points ) );
        subsetLongitudes[ index ] = longitude;
        subsetLatitudes[  index ] = latitude;

        if ( subsetLongitudesSW ) {
          subsetLongitudesSW[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ SW ][ row ][ column ];
          subsetLongitudesSE[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ SE ][ row ][ column ];
          subsetLongitudesNW[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ NW ][ row ][ column ];
          subsetLongitudesNE[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ NE ][ row ][ column ];
          subsetLatitudesSW[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ SW ][ row ][ column ];
          subsetLatitudesSE[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ SE ][ row ][ column ];
          subsetLatitudesNW[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ NW ][ row ][ column ];
          subsetLatitudesNE[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ NE ][ row ][ column ];
        }

        ++index;
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
INPUTS:  const Real data[ VARIABLES ][ ROWS ][ COLUMNS ]  Decoded data to copy.
         const Integer selected[ VARIABLES ]  Flags for each variable selected.
         Integer     firstRow     Index of first row of subset.
         Integer     lastRow      Index of alst row os subset.
         Integer     firstColumn  Index of first column of subset.
         Integer     lastColumn   Index of last column of subset.
         Real        output[ variables * points ] Array to append data to.
******************************************************************************/

static void copySubsetData( const Real data[ VARIABLES ][ ROWS ][ COLUMNS ],
                            const Integer selected[ VARIABLES ],
                            Integer firstRow, Integer lastRow,
                            Integer firstColumn, Integer lastColumn,
                            Real output[] ) {

  PRE012( data,
          isNanFree( &data[0][0][0], VARIABLES * ROWS * COLUMNS ),
          minimumItem( &data[ MSK ][0][0], ROWS * COLUMNS ) >= 0.0,
          maximumItem( &data[ MSK ][0][0], ROWS * COLUMNS ) <= 1.0,
          selected,
          minimumItemI( selected, VARIABLES ) >= 0,
          maximumItemI( selected, VARIABLES ) == 1,
          IN_RANGE( firstRow, 0, ROWS - 1 ),
          IN_RANGE( lastRow, firstRow, ROWS - 1 ),
          IN_RANGE( firstColumn, 0, COLUMNS - 1 ),
          IN_RANGE( lastColumn, firstColumn, COLUMNS - 1 ),
          output );

  CHECKING( const Integer variables = sumI( selected, VARIABLES ); )
  CHECKING( Integer points = 0; )
  Integer variable = 0;
  Real* values = output;

  for ( variable = 0; variable < VARIABLES; ++variable ) {

    if ( selected[ variable ] ) {
      Integer row = 0;

      for ( row = firstRow; row <= lastRow; ++row ) {
        Integer column = 0;

        for ( column = firstColumn; column <= lastColumn; ++column ) {
          const Real mask = data[ MSK ][ row ][ column ];

          if ( mask ) {
            const Real value = data[ variable ][ row ][ column ];
            CHECKING(  ++points; )
            CHECK2( ! isNan( value ),
                    IN_RANGE( points, 1,
                              variables *
                              ( 1 + lastRow - firstRow ) *
                              ( 1 + lastColumn - firstColumn ) ) );
            *values++ = value;
          }
        }
      }
    }
  }

  POST0( isNanFree( output, points ) );
}



/******************************************************************************
PURPOSE: computeMean - Compute mean of scan points subsetted by indices/ranges.
INPUTS:  const Bounds bounds                  Subset lon-lat bounds.
         const Real ranges[ VARIABLES ][ 2 ]  Data ranges to filter by.
         const Integer indices[ 2 ][ 2 ]      indices[ROW COLUMN][MIN/MAXIMUM].
         const Real data[ VARIABLES ][ ROWS ][ COLUMNS ]  Decoded scan data.
         Real longitudes[ ROWS ][ COLUMNS ]   Read from a file.
         Real latitudes[ ROWS ][ COLUMNS ]    Read from a file.
         Real data[ VARIABLES ][ ROWS ][ COLUMNS ]  Decoded data.
         Integer counts[ ROWS ][ COLUMNS ]    Current mean counts.
         Real values[ ROWS ][COLUMNS]         Current means.
OUTPUTS: Integer counts[ ROWS ][ COLUMNS ]    Updated mean counts.
         Real values[ ROWS ][COLUMNS]         Updated means.
******************************************************************************/

static void computeMean( const Bounds bounds,
                         const Real ranges[ VARIABLES ][ 2 ],
                         const Integer indices[ 2 ][ 2 ],
                         Real data[ VARIABLES ][ ROWS ][ COLUMNS ] ) {

  PRE07( isValidBounds( bounds ),
         indices,
         IN_RANGE( indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
         IN_RANGE( indices[ ROW ][ MAXIMUM ],
                   indices[ ROW ][ MINIMUM ], ROWS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                   indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ),
         isNanFree( &ranges[0][0], VARIABLES * 2 ) );

  const Integer firstRow    = indices[ ROW    ][ MINIMUM ];
  const Integer lastRow     = indices[ ROW    ][ MAXIMUM ];
  const Integer firstColumn = indices[ COLUMN ][ MINIMUM ];
  const Integer lastColumn  = indices[ COLUMN ][ MAXIMUM ];
  const Real longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const Real longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const Real latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const Real latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const Real aodMinimum = ranges[ AOD ][ MINIMUM ];
  const Real aodMaximum = ranges[ AOD ][ MAXIMUM ];
  const Real clsMinimum = ranges[ CLS ][ MINIMUM ];
  const Real clsMaximum = ranges[ CLS ][ MAXIMUM ];
  const Real stdMinimum = ranges[ STD ][ MINIMUM ];
  const Real stdMaximum = ranges[ STD ][ MAXIMUM ];
  const Real sfcMinimum = ranges[ SFC ][ MINIMUM ];
  const Real sfcMaximum = ranges[ SFC ][ MAXIMUM ];
  const Real ch1Minimum = ranges[ CH1 ][ MINIMUM ];
  const Real ch1Maximum = ranges[ CH1 ][ MAXIMUM ];
  const Real mosMinimum = ranges[ MOS ][ MINIMUM ];
  const Real mosMaximum = ranges[ MOS ][ MAXIMUM ];
  const Real cldMinimum = ranges[ CLD ][ MINIMUM ];
  const Real cldMaximum = ranges[ CLD ][ MAXIMUM ];
  const Real sigMinimum = ranges[ SIG ][ MINIMUM ];
  const Real sigMaximum = ranges[ SIG ][ MAXIMUM ];
  const Real scaMinimum = ranges[ SCA ][ MINIMUM ];
  const Real scaMaximum = ranges[ SCA ][ MAXIMUM ];
  Integer row = 0;

  /* Compute reduced mask in parallel: */

#pragma omp parallel for

  for ( row = firstRow; row <= lastRow; ++row ) {
    Integer column = 0;

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      const Real longitude = longitudesLatitudes[ LONGITUDE ][ row ][ column ];
      const Real latitude  = longitudesLatitudes[ LATITUDE  ][ row ][ column ];
      const Real aod = data[ AOD ][ row ][ column ];
      const Real msk = data[ MSK ][ row ][ column ];
      const Real cls = data[ CLS ][ row ][ column ];
      const Real std = data[ STD ][ row ][ column ];
      const Real sfc = data[ SFC ][ row ][ column ];
      const Real ch1 = data[ CH1 ][ row ][ column ];
      const Real mos = data[ MOS ][ row ][ column ];
      const Real cld = data[ CLD ][ row ][ column ];
      const Real sig = data[ SIG ][ row ][ column ];
      const Real sca = data[ SCA ][ row ][ column ];
      const Integer output =
        AND12( msk == 1.0,
               IN_RANGE( longitude, longitudeMinimum, longitudeMaximum ),
               IN_RANGE( latitude, latitudeMinimum, latitudeMaximum ),
               IN_RANGE( aod, aodMinimum, aodMaximum ),
               IN_RANGE( cls, clsMinimum, clsMaximum ),
               IN_RANGE( std, stdMinimum, stdMaximum ),
               IN_RANGE( sfc, sfcMinimum, sfcMaximum ),
               IN_RANGE( ch1, ch1Minimum, ch1Maximum ),
               IN_RANGE( mos, mosMinimum, mosMaximum ),
               IN_RANGE( cld, cldMinimum, cldMaximum ),
               IN_RANGE( sig, sigMinimum, sigMaximum ),
               IN_RANGE( sca, scaMinimum, scaMaximum ) );

      if ( output ) {
        const Real mean = means[ row ][ column ];
        const Integer count = counts[ row ][ column ];
        means[ row ][ column ] = ( count * mean + aod ) / ( count + 1 );
        counts[ row ][ column ] += 1;
      }
    }
  }

  POST0( isNanFree( &means[0][0], ROWS * COLUMNS ) );
}



/******************************************************************************
PURPOSE: meanPoints - Count and return number of non-zero counts[][].
INPUTS:  const Integer indices[ 2 ][ 2 ]      indices[ROW COLUMN][MIN/MAXIMUM].
         Integer counts[ ROWS ][ COLUMNS ]    Final mean counts.
RETURNS: Integer number of non-zero counts[][].
******************************************************************************/

static Integer meanPoints( const Integer indices[ 2 ][ 2 ] ) {

  PRE05( indices,
         IN_RANGE( indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
         IN_RANGE( indices[ ROW ][ MAXIMUM ],
                   indices[ ROW ][ MINIMUM ], ROWS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
         IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                   indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ) );

  const Integer firstRow    = indices[ ROW    ][ MINIMUM ];
  const Integer lastRow     = indices[ ROW    ][ MAXIMUM ];
  const Integer firstColumn = indices[ COLUMN ][ MINIMUM ];
  const Integer lastColumn  = indices[ COLUMN ][ MAXIMUM ];
  Integer row = 0;
  Integer result = 0;

  for ( row = firstRow; row <= lastRow; ++row ) {
    Integer column = 0;

    for ( column = firstColumn; column <= lastColumn; ++column ) {
      result += counts[ row ][ column ] != 0;
    }
  }

  POST0( IN_RANGE( result, 0, ROWS * COLUMNS ) );
  return result;
}



/******************************************************************************
PURPOSE: copyMeanData - Copy subset of mean data.
INPUTS:  const Integer indices[ 2 ][ 2 ]    indices[ROW COLUMN][MIN/MAXIMUM].
         const Integer points               Subset points.
         Integer counts[ ROWS ][ COLUMNS ]  Final mean counts.
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

static void copyMeanData( const Integer indices[ 2 ][ 2 ],
                          const Integer points,
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

  PRE011( indices,
          IN_RANGE( indices[ ROW ][ MINIMUM ], 0, ROWS - 1 ),
          IN_RANGE( indices[ ROW ][ MAXIMUM ],
                    indices[ ROW ][ MINIMUM ], ROWS - 1 ),
          IN_RANGE( indices[ COLUMN ][ MINIMUM ], 0, COLUMNS - 1 ),
          IN_RANGE( indices[ COLUMN ][ MAXIMUM ],
                    indices[ COLUMN ][ MINIMUM ], COLUMNS - 1 ),
          IN_RANGE( points, 1, ROWS * COLUMNS ),
          IN_RANGE( points, 1,
                    ( 1 + indices[ ROW ][ MAXIMUM ] - indices[ROW][MINIMUM] ) *
                    (1 + indices[COLUMN][MAXIMUM] - indices[COLUMN][MINIMUM])),
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

  const Integer firstRow    = indices[ ROW    ][ MINIMUM ];
  const Integer lastRow     = indices[ ROW    ][ MAXIMUM ];
  const Integer firstColumn = indices[ COLUMN ][ MINIMUM ];
  const Integer lastColumn  = indices[ COLUMN ][ MAXIMUM ];
  Integer row = 0;
  Integer index = 0;

  for ( row = firstRow; row <= lastRow; ++row ) {
    Integer column = 0;

    for ( column = firstColumn; column <= lastColumn; ++column ) {

      if ( counts[ row ][ column ] ) {
        const Real longitude =
          longitudesLatitudes[ LONGITUDE ][ row ][ column ];
        const Real latitude = longitudesLatitudes[ LATITUDE ][ row ][ column ];
        const Real mean = means[ row ][ column ];
        CHECK( IN_RANGE( index, 0, points - 1 ) );
        subsetLongitudes[ index ] = longitude;
        subsetLatitudes[  index ] = latitude;
        subsetData[       index ] = mean;

        if ( subsetLongitudesSW ) {
          subsetLongitudesSW[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ SW ][ row ][ column ];
          subsetLongitudesSE[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ SE ][ row ][ column ];
          subsetLongitudesNW[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ NW ][ row ][ column ];
          subsetLongitudesNE[ index ] =
            longitudeLatitudeCorners[ LONGITUDE ][ NE ][ row ][ column ];
          subsetLatitudesSW[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ SW ][ row ][ column ];
          subsetLatitudesSE[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ SE ][ row ][ column ];
          subsetLatitudesNW[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ NW ][ row ][ column ];
          subsetLatitudesNE[ index ] =
            longitudeLatitudeCorners[ LATITUDE ][ NE ][ row ][ column ];
        }

        ++index;
      }
    }
  }

  CHECK( index == points );
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

  PRE06( isValidData( data ), data->ok, output, output->invariant( output ),
         output->isWritable( output ), data->ok == output->ok( output ) );

  const Arguments* const arguments = &data->arguments;
  const VoidList*  const scans = data->subsettedScans;
  const SubsettedScan* const firstScan = scans->item( scans, 0 );
  const Integer variables = firstScan->variables;
  const Integer writeCorners = arguments->corners;
  const char* const daily = arguments->daily ? "daily_" : "";
  Integer variable = 0;
  UTCTimestamp timestamp;
  toUTCTimestamp( arguments->firstTimestamp, timestamp );
  CHECK( variables > 0 );

  output->writeString( output,
                       "Swath 2.0\n%s\n%s\n"
                       "# Dimensions: variables timesteps scans:\n"
                       "%"INTEGER_FORMAT" %"INTEGER_FORMAT" %"
                       INTEGER_FORMAT"\n"
                       "# Variable names:\nLongitude Latitude",
                       arguments->description, timestamp,
                       variables + writeCorners * 8,
                       arguments->hours, scans->count( scans ) );

  for ( variable = 0;
        AND2( output->ok( output ), variable < VARIABLES );
        ++variable ) {

    if ( arguments->selected[ variable ] ) {
      const char* const variableName = variableNames[ variable ];
      CHECK2( variableName, variableName[ 0 ] );
      output->writeString( output, " %s%s", daily, variableName );
    }
  }

  if ( writeCorners ) {
    output->writeString( output,
                         " Longitude_SW Longitude_SE"
                         " Longitude_NW Longitude_NE"
                         " Latitude_SW Latitude_SE"
                         " Latitude_NW Latitude_NE" );
  }

  if ( output->ok( output ) ) {
    output->writeString( output, "\n# Variable units:\ndeg deg" );
  }

  for ( variable = 0;
        AND2( output->ok( output ), variable < VARIABLES );
        ++variable ) {

    if ( arguments->selected[ variable ] ) {

      if ( variable == SCA ) {
        output->writeString( output, " deg" );
      } else {
        output->writeString( output, " -" );
      }
    }
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
                             ( scan->variables + scan->hasCorners * 8 )
                             * scan->points );
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



/******************************************************************************
PURPOSE: clampLongitudes - Clamp cell longitudes to match sign of first one.
INPUTS:  const double longitude              First longitude center point.
         double* nextColumnLongitude         Longitude to check/clamp.
         double* nextRowLongitude            Longitude to check/clamp.
         double* nextRowNextColumnLongitude  Longitude to check/clamp.
OUTPUTS: double* nextColumnLongitude         Signed to match longitude.
         double* nextRowLongitude            Signed to match longitude.
         double* nextRowNextColumnLongitude  Signed to match longitude.
******************************************************************************/

static void clampLongitudes( const double longitude,
                             double* nextColumnLongitude,
                             double* nextRowLongitude,
                             double* nextRowNextColumnLongitude ) {

  PRE03( nextColumnLongitude, nextRowLongitude, nextRowNextColumnLongitude );

  if ( longitude < -179.0 ) {

    if ( *nextColumnLongitude > 0.0 ) {
      *nextColumnLongitude = -*nextColumnLongitude;
    }

    if ( *nextRowLongitude > 0.0 ) {
      *nextRowLongitude = -*nextRowLongitude;
    }

    if ( *nextRowNextColumnLongitude > 0.0 ) {
      *nextRowNextColumnLongitude = -*nextRowNextColumnLongitude;
    }
  } else if ( longitude > 179.0 ) {

    if ( *nextColumnLongitude < 0.0 ) {
      *nextColumnLongitude = -*nextColumnLongitude;
    }

    if ( *nextRowLongitude < 0.0 ) {
      *nextRowLongitude = -*nextRowLongitude;
    }

    if ( *nextRowNextColumnLongitude < 0.0 ) {
      *nextRowNextColumnLongitude = -*nextRowNextColumnLongitude;
    }
  }

  POST03( SIGN( longitude ) == SIGN( *nextColumnLongitude ),
          SIGN( longitude ) == SIGN( *nextRowLongitude ),
          SIGN( longitude ) == SIGN( *nextRowNextColumnLongitude ) );
}




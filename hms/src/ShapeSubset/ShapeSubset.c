
/******************************************************************************
PURPOSE: ShapeSubset.c - Subset by lon-lat polygon Shapefiles (shp,shx,dbf).
NOTES:   Uses Shapefile and GPC free, open-source libraries.
         http://shapelib.maptools.org/shp_api.html
         http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
         Compile (with debugging):
           gcc -Wall -m64 -g -DDEBUGGING -I../../../../code/include \
               -o ShapeSubset ShapeSubset.c \
               -L../../../../code/lib/Darwin \
               -lUtilities.DEBUG -lShapefile -lGPC
         Compile (with optimization):
           gcc -Wall -m64 -O -DNO_ASSERTIONS -I../../../../code/include \
               -o ShapeSubset ShapeSubset.c \
               -L../../../../code/lib/Darwin -lUtilities -lShapefile -lGPC
         Run:
           ShapeSubset -input input [input2 ...] \
                       -subset minlon minlat maxlon maxlat \
                       [-mindist mindist] \
                       [-huc 109000107] \
                       [-comid 4570503] \
                       [-estcode INDI] \
                       [-time yyyymmdd1 yyyymmdd2] \
                       [-dbf_only] \
                       -output output > output.bin
           Note: mindist is the 'minimum adjacent vertex distance' so that
           vertices closer than this (in either x or y) are merged to sparse
           the resulting clipped polygons. Use -mindist 0 to disable sparse.
         Examples:
           ShapeSubset -input temperature_gulf -subset -90 28 -85 35 \
                       -mindist 0.001 -output temperature > output.bin
           unbin output.bin
           ls -l temperature.???
           open temperature.shp
           ShapeSubset -input sediment_nca_gulf -subset -90 28 -85 35 \
                       -dbf_only -output sediment > sediment.bin
           unbin sediment.bin
           dbfdump sediment.dbf
HISTORY: 2011-08-25 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For FILE, stderr, fprintf(). */
#include <stdlib.h> /* For malloc(), free(). */
#include <string.h> /* For memset(), strrchr. */
#include <stdlib.h> /* For malloc(), free(), qsort(). */
#include <dirent.h> /* For DIR, opendir(), closedir(). */

#include <Assertions.h>    /* For PRE*(), POST*(). */
#include <BasicNumerics.h> /* For isValidArgs(). */
#include <Utilities.h>     /* For Bounds, isValidBounds(). */
#include <DateTime.h>      /* For isValidYearMonthDay(). */
#include <Shapefile.h>     /* For PolygonShape, readAndClipShapes(). */

/*=========================== FORWARD DECLARATIONS ==========================*/

static int intCompare( const void* a, const void* b ) {
  const int* const pa = (const int*) a;
  const int* const pb = (const int*) b;
  const int ia = *pa;
  const int ib = *pb;
  return ia - ib;
}

static int isDirectory( const char* const name ) {
  DIR* directory = name ? opendir( name ) : 0;
  const int result = directory ? closedir( directory ) == 0 : 0;
  return result;
}

static void usage( const char* program );

static int parseArguments( int argc, char* argv[],
                           const char** baseOutputFileName,
                           const char** csvFileName,
                           Bounds subset, long long* huc, int* comid,
                           const char** estcode,
                           double* mindist, int* yyyymmdd1, int* yyyymmdd2,
                           int* dbfOnly, int* noStream );

static int writeMaskedSubsetDBF( const char* const baseInputFileName,
                                 const char* const csvFileName,
                                 const int dbfOnly,
                                 const long long huc,
                                 const int comid,
                                 const char* const estcode,
                                 const Bounds subset,
                                 const int yyyymmdd1,
                                 const int yyyymmdd2,
                                 char* const mask,
                                 DBFHandle outputDBF,
                                 int* wroteDBF );

static int writeMaskedSubsetSHP( const char* const baseInputFileName,
                                 const char* const baseOutputFileName,
                                 const double minimumAdjacentVertexDistance,
                                 const int clipShapes,
                                 const Bounds subset,
                                 char* const mask,
                                 int* const offset,
                                 DBFHandle outputDBF,
                                 SHPHandle* const outputSHP );

static int writeMaskedSubsetCSV( const char* const baseInputFileName,
                                 const char* const baseOutputFileName,
                                 const char* const csvFileName,
                                 const char* const estcode,
                                 const int yyyymmdd1,
                                 const int yyyymmdd2,
                                 char* const mask );

static int parseUniqueSiteIds( const char* csvData,
                               const int lines,
                               int siteIds[] );

static int getSiteIdColumn( const char* csv_data );

static const char* getColumnValue( const char* line, const int column );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: main - Main routine.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line arguments.
RETURNS: int 0 if successful, else 1 and failure message is printed.
******************************************************************************/

int main( int argc, char* argv[] ) {
  int ok = 0;
  const char* baseOutputFileName = 0;
  const char* csvFileName = 0;
  Bounds subset = { { 0.0, 0.0 }, { 0.0, 0.0 } };
  double minimumAdjacentVertexDistance = 0.0;
  long long huc = 0; /* 14-digit integer HUC id to filter/subset shp/dbf. */
  int comid = 0; /* 10-digit integer COMID to filter/subset shp/dbf. */
  int yyyymmdd1 = 0; /* Optional 1st timestamp to filter by. */
  int yyyymmdd2 = 0; /* Optional 2nd timestamp to filter by. */
  int dbfOnly = 0;
  int noStream = 0;
  int inputFileIndex = 0;
  const char* estcode = 0;
  const int inputFileCount =
    parseArguments( argc, argv, &baseOutputFileName, &csvFileName,
                    subset, &huc, &comid, &estcode,
                    &minimumAdjacentVertexDistance,
                    &yyyymmdd1, &yyyymmdd2, &dbfOnly, &noStream );

  /*
   * Read 1 or more input dbf and possibly matching shp file too and
   * write subset output to a dbf file and sometimes also a shp file:
   */

  if ( inputFileCount > 0 ) {
    DBFHandle outputDBF = DBFCreate( baseOutputFileName );

    if ( ! outputDBF ) {
      fprintf( stderr, "Failed to create temporary dbf file '%s'.\n",
               baseOutputFileName );
    } else {
      int offset = 0;
      SHPHandle outputSHP = 0;
      DEBUG( fprintf( stderr, "inputFileCount = %d, "
                      "baseOutputFileName = %s, csvFileName = '%s', "
                      "subset = [%lg %lg][%lg %lg], "
                      "yyyymmdd1 = %d, yyyymmdd2 = %d, "
                      "huc = %lld, comid = %d, estcode = '%s', "
                      "minimumAdjacentVertexDistance = %lf, "
                      "dbfOnly = %d, noStream = %d\n",
                      inputFileCount, baseOutputFileName,
                      csvFileName ? csvFileName : "",
                      subset[ LONGITUDE ][ MINIMUM ],
                      subset[ LONGITUDE ][ MAXIMUM ],
                      subset[ LATITUDE  ][ MINIMUM ],
                      subset[ LATITUDE  ][ MAXIMUM ],
                      yyyymmdd1, yyyymmdd2,
                      huc, comid, estcode ? estcode : "",
                      minimumAdjacentVertexDistance,
                      dbfOnly, noStream ); )

      for ( inputFileIndex = 0, ok = 1;
            AND2( ok, inputFileIndex < inputFileCount );
            ++inputFileIndex ) {
        const char* const baseInputFileName = argv[ 2 + inputFileIndex ];
        const int rows = getRowsDBF( baseInputFileName );
        char* mask = rows > 0 ? malloc( rows * sizeof (char) ) : 0;
        ok = mask != 0;

        if ( mask ) {
          int wroteDBF = 0;
          const int clipShapes =
            AND5( ! dbfOnly, huc == 0, comid == 0, estcode == 0,
                  ! isPointType( baseInputFileName ));
          DEBUG( fprintf( stderr, "clipShapes = %d\n", clipShapes ); )

          memset( mask, 1, rows * sizeof (char) );

          if ( ! clipShapes ) {

            /* Write subset dbf file and/or set mask based on various filters:*/

            ok = writeMaskedSubsetDBF( baseInputFileName, csvFileName,
                                       dbfOnly, huc, comid, estcode,
                                       (const double (*)[2]) subset,
                                       yyyymmdd1, yyyymmdd2,
                                       mask, outputDBF, &wroteDBF );

            CHECK( IMPLIES( AND2( ok, dbfOnly ), wroteDBF ) );
            CHECK( IMPLIES( wroteDBF, offset == 0 ) );
          }

          if ( ok ) {

            if ( ! dbfOnly ) {

              /*
               * Write subset shp file and maybe set mask to subset clip.
               * Some usages require deferring output dbf until after mask
               * has been set based on clipping shp to subset. In such cases,
               * the DBF file will be written by the below routine.
               */

              const int polygonCount =
                writeMaskedSubsetSHP( baseInputFileName,
                                      baseOutputFileName,
                                      minimumAdjacentVertexDistance,
                                      clipShapes,
                                      (const double (*)[2]) subset,
                                      mask, &offset, outputDBF, &outputSHP );

              ok = OR2( polygonCount > 0,
                        strstr( baseInputFileName, "hms_smoke" ) );
            }

            if ( ok ) {

              /*
               * If csvFileName specified then write subset csv based on
               * various id or time filtering.
               */

              if ( csvFileName ) {
                const int alreadyProcessedCSV =
                  OR7( strstr( baseInputFileName, "hspf_" ),
                       strstr( baseInputFileName, "swat_" ),
                       strstr( baseInputFileName, "swmm_" ),
                       strstr( baseInputFileName, "esat1_" ),
                       strstr( baseInputFileName, "esat2_" ),
                       strstr( baseInputFileName, "seagrass" ),
                       AND3( comid == 0, estcode,
                             OR2( strstr( baseInputFileName, "seagrass" ),
                                  strstr( baseInputFileName,
                                         "flowlines_puget_sound_watershed"))));

                /*
                 * WSM csv is a temporary file that was already subsetted by
                 * landuseserver routine subset_csv_file()
                 * so no need to write another, just include it in the streamed
                 * output.
                 */

                if ( ! alreadyProcessedCSV ) {
                  const int isSubsetFlowlines =
                    OR2( strstr( baseInputFileName, "flowlines" ),
                         strstr( baseInputFileName, "RBEROST" ) );
                  const int valuesWritten =
                    writeMaskedSubsetCSV( baseInputFileName,
                                          baseOutputFileName, csvFileName,
                                          estcode, yyyymmdd1, yyyymmdd2,
                                          mask );
                  ok = OR2( isSubsetFlowlines, valuesWritten > 0 );
                }
              }
            }
          }

          free( mask ), mask = 0;
        } /* end if ( mask ) */
      } /* Loop on input files. */

      if ( outputSHP ) {
        SHPClose( outputSHP ), outputSHP = 0;
      }

      DBFClose( outputDBF ), outputDBF = 0;
    } /* outputDBF created. */

    if ( ok ) {

      if ( ! noStream ) {

        /*
         * Stream set of subset file bytes in 'bin' format (with one line
         * header) to stdout:
         */

        const char* const slash = strrchr( baseOutputFileName, '/' );
        const char* const unpathedName = slash ? slash + 1: baseOutputFileName;
        ok = streamShapefiles( baseOutputFileName, unpathedName, dbfOnly,
                               csvFileName != 0 );
        removeShapefiles( baseOutputFileName );
      }
    } else {
      fprintf( stderr, "%s: No shapes were in the subset.\n", argv[ 0 ] );
    }
  }

  DEBUG( fprintf( stderr, "main returning %d\n", ! ok ); )
  return ! ok;
}



/*============================= PRIVATE FUNCTIONS ===========================*/



/******************************************************************************
PURPOSE: usage - Print program usage.
OUTPUTS: const char* program  Name of program.
******************************************************************************/

static void usage( const char* program ) {
  PRE02( program, *program );
  fprintf( stderr, "\n%s - ", program );
  fprintf( stderr, "Subset by lon-lat or HUC Shapefiles (shp,shx,dbf).\n" );
  fprintf( stderr, "usage: %s -input input [input2 ...] ", program );
  fprintf( stderr, " [-subset minlon minlat maxlon maxlat] | [-huc id] " );
  fprintf( stderr, " [-estcode code] " );
  fprintf( stderr, " [-mindist mindist] [-time yyyymmdd1 yyyymmdd2] " );
  fprintf( stderr, " [-dbf_only] " );
  fprintf( stderr, " [-csv file] -output output " );
  fprintf( stderr, " [-no_stream] | > output.bin\n" );
  fprintf( stderr, "\n" );
  fprintf( stderr, "Examples:\n" );
  fprintf( stderr, "%s temperature_annual_gulf -subset -90 28 -85 35 ",
           program );
  fprintf( stderr, "-mindist 0.007 -output temperature > output.bin\n" );
  fprintf( stderr, "head -1 output.bin\n" );
  fprintf( stderr, "unbin output.bin\n" );
  fprintf( stderr, "ls -l temperature.???\n" );
  fprintf( stderr, "After line 1 is a concatenation of shx,shp,dbf.\n\n" );
  fprintf( stderr,
           "%s -input sediment_nca_gulf -subset -90 28 -85 35 ", program );
  fprintf( stderr, "-dbf_only -output sediment > sediment.bin\n" );
  fprintf( stderr, "unbin sediment.bin\n" );
  fprintf( stderr, "dbfdump sediment.dbf\n\n" );
  fprintf( stderr,
           "%s -input /data/land_use/wmost_charles3_huc10_point ", program );
  fprintf( stderr, "-huc 109000107 -dbf_only " );
  fprintf( stderr, "-csv /data/tmp/wmost_charles3_loadings_tn.csv " );
  fprintf( stderr, "-output /data/tmp/wmost_charles3_huc10_point_loadings_tn");
  fprintf( stderr, " -no_stream\n");
  fprintf( stderr,
           "ls -l /data/tmp/wmost_charles3_huc10_point_loadings_tn.???\n" );
  fprintf( stderr, "\n" );
  fprintf( stderr,
           "%s -input /data/land_use/stream_discharge_nhd_line_gulf ", program);
  fprintf( stderr, "-subset -85 25 -80 30 -time 20080701 20080930 " );
  fprintf( stderr, "-estcode ALLI " );
  fprintf( stderr, "-csv /data/land_use/discharge/monthly/gulf " );
  fprintf( stderr, "-output /data/tmp/monthly_stream_discharge_nhd_line");
  fprintf( stderr, "\n" );
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
OUTPUTS: int argc                         Number of command-line arguments.
         char argv[]                      Command-line arguments.
OUTPUTS: const char** baseOutputFileName  Base name of shp/dbf file to write.
         const char** csvFileName         Base name of csv file to stream.
         Bounds subset                    Subset to apply to grid file.
         long long* huc                   HUC_ID to match instead of subset.
         int* comid                       COMID to match instead of subset.
         const char** estcode             ESTCODE to match instead of subset.
         double* mindist                  Minimum adjacent vertex distance.
         int* yyyymmdd1                   1st timestamp of subset or 0.
         int* yyyymmdd2                   2nd timestamp of subset or 0.
         int* dbfOnly                     Only process DBF file?
         int* noStream                    Don't stream files to stdout?
RETURNS: int input file count >= 1 if successful, else 0 and a failure message
         is printed to stderr.
******************************************************************************/

static int parseArguments( int argc, char* argv[],
                           const char** baseOutputFileName,
                           const char** csvFileName,
                           Bounds subset, long long* huc, int* comid,
                           const char** estcode,
                           double* mindist, int* yyyymmdd1, int* yyyymmdd2,
                           int* dbfOnly, int* noStream ) {

  PRE010( baseOutputFileName, csvFileName,
          subset, subset, estcode, huc, mindist, yyyymmdd1, yyyymmdd2, dbfOnly);

  int result = isValidArgs( argc, (const char**) argv );
  int inputFileCount = 0;
  int parsedSubset = 0;
  int parsedTime = 0;
  int parsedMindist = 0;
  int arg = 1;
  *baseOutputFileName = 0;
  *csvFileName = 0;
  subset[ LONGITUDE ][ MINIMUM ] = -180.0;
  subset[ LONGITUDE ][ MAXIMUM ] =  180.0;
  subset[ LATITUDE  ][ MINIMUM ] =  -90.0;
  subset[ LATITUDE  ][ MAXIMUM ] =   90.0;
  *huc = 0;
  *comid = 0;
  *mindist = 0.0;
  *yyyymmdd1 = 0;
  *yyyymmdd2 = 0;
  *dbfOnly = 0;
  *noStream = 0;
  *estcode = 0;

  while ( AND2( result, arg < argc ) ) {

    if ( AND4( ! strcmp( argv[ arg ], "-input" ),
               arg + 1 < argc,
               argv[ arg + 1 ][ 0 ] != '-',
               inputFileCount == 0 ) ) {
      ++arg;

      while ( AND2( arg + 1 < argc, argv[ arg ][ 0 ] != '-' ) ) {
        ++inputFileCount;
        ++arg;
      }

    } else if ( AND2( ! strcmp( argv[ arg ], "-output" ), arg + 1 < argc ) ) {

      if ( *baseOutputFileName ) {
        fprintf( stderr, "Redundant -output arguments.\n" );
        result = 0;
      } else {
        *baseOutputFileName = argv[ arg + 1 ];
        arg += 2;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-csv" ), arg + 1 < argc ) ) {

      if ( *csvFileName ) {
        fprintf( stderr, "Redundant -csv arguments.\n" );
        result = 0;
      } else {
        *csvFileName = argv[ arg + 1 ];
        arg += 2;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-dbf_only" ), arg < argc ) ) {

      if ( *dbfOnly ) {
        fprintf( stderr, "Redundant -dbf_only arguments.\n" );
        result = 0;
      } else {
        *dbfOnly = 1;
        ++arg;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-mindist" ), arg + 1 < argc ) ) {

      if ( parsedMindist ) {
        fprintf( stderr, "Redundant -mindist arguments.\n" );
        result = 0;
      } else {
        *mindist = atof( argv[ arg + 1 ] );

        if ( *mindist >= 0.0 ) {
          arg += 2;
          parsedMindist = 1;
        } else {
          result = 0;
        }
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-no_stream" ), arg < argc ) ) {

      if ( *noStream ) {
        fprintf( stderr, "Redundant -no_stream arguments.\n" );
        result = 0;
      } else {
        *noStream = 1;
        ++arg;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-estcode" ), arg + 1 < argc ) ) {

      if ( *huc ) {
        fprintf( stderr, "Redundant -estcode arguments.\n" );
        result = 0;
      } else {
        *estcode = argv[ arg + 1 ];
        arg += 2;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-huc" ), arg + 1 < argc ) ) {

      if ( *huc ) {
        fprintf( stderr, "Redundant -huc arguments.\n" );
        result = 0;
      } else if ( parsedSubset ) {
        fprintf( stderr, "Conflicting usage: -huc / -subset arguments.\n" );
        result = 0;
      } else {
        *huc = atoll( argv[ arg + 1 ] );
        arg += 2;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-comid" ), arg + 1 < argc ) ) {

      if ( *comid ) {
        fprintf( stderr, "Redundant -comid arguments.\n" );
        result = 0;
      } else if ( parsedSubset ) {
        fprintf( stderr, "Conflicting usage: -comid / -subset arguments.\n" );
        result = 0;
      } else {
        *comid = atoll( argv[ arg + 1 ] );
        arg += 2;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-subset" ), arg + 4 < argc ) ) {
      subset[ LONGITUDE ][ MINIMUM ] = atof( argv[ arg + 1 ] );
      subset[ LATITUDE  ][ MINIMUM ] = atof( argv[ arg + 2 ] );
      subset[ LONGITUDE ][ MAXIMUM ] = atof( argv[ arg + 3 ] );
      subset[ LATITUDE  ][ MAXIMUM ] = atof( argv[ arg + 4 ] );

      if ( parsedSubset ) {
        fprintf( stderr, "Redundant -subset arguments.\n" );
        result = 0;
      } else if ( *huc ) {
        fprintf( stderr, "Conflicting usage: -huc / -subset arguments.\n" );
        result = 0;
      } else if ( ! isValidBounds( (const double(*)[2]) subset ) ) {
        fprintf( stderr, "Invalid -subset arguments.\n" );
        result = 0;
      } else {
        parsedSubset = 1;
        arg += 5;
      }
    } else if ( AND2( ! strcmp( argv[ arg ], "-time" ), arg + 2 < argc ) ) {
      *yyyymmdd1 = atoi( argv[ arg + 1 ] );
      *yyyymmdd2 = atoi( argv[ arg + 2 ] );

      if ( parsedTime ) {
        fprintf( stderr, "Redundant -time arguments.\n" );
        result = 0;
      } else if ( *huc ) {
        fprintf( stderr, "Conflicting usage: -huc / -time arguments.\n" );
        result = 0;
      } else if ( ! AND3( isValidYearMonthDay( *yyyymmdd1 ),
                          isValidYearMonthDay( *yyyymmdd2 ),
                          *yyyymmdd1 <= *yyyymmdd2 ) ) {
        fprintf( stderr, "Invalid -time arguments.\n" );
        result = 0;
      } else {
        parsedTime = 1;
        arg += 3;
      }
    } else {
      result = 0;
    }
  }

  if ( result ) {
    result = AND5( inputFileCount > 0,
                   *baseOutputFileName,
                   *baseOutputFileName[ 0 ],
                   strcmp( *baseOutputFileName, argv[ 0 ] ),
                   IMPLIES( *csvFileName,
                            AND3( *csvFileName[ 0 ],
                                  strcmp( *csvFileName, argv[ 0 ] ),
                                  strcmp( *csvFileName,*baseOutputFileName))));
  }

  if ( ! result ) {
    *baseOutputFileName = 0;
    *csvFileName = 0;
    memset( subset, 0, sizeof (Bounds) );
    *huc = 0;
    *comid = 0;
    *mindist = 0.0;
    *yyyymmdd1 = 0;
    *yyyymmdd2 = 0;
    *dbfOnly = 0;
    *estcode = 0;
    usage( argv[ 0 ] );
  } else {
    result = inputFileCount;
  }

  POST02( result >= 0,
          IMPLIES_ELSE( result,
                        AND4( isValidBounds( (const double (*)[2]) subset ),
                              IN_RANGE( *mindist, 0.0, 1.0 ),
                              IS_BOOL( *dbfOnly ),
                              IMPLIES( *yyyymmdd1,
                                       AND3( isValidYearMonthDay( *yyyymmdd1 ),
                                             isValidYearMonthDay( *yyyymmdd2 ),
                                             *yyyymmdd1 <= *yyyymmdd2 ) ) ),
                        IS_ZERO9( subset[ 0 ][ 0 ], subset [ 0 ][ 1 ],
                                  subset[ 1 ][ 0 ], subset [ 1 ][ 1 ],
                                  *estcode,
                                  *mindist, *yyyymmdd1, *yyyymmdd2, *dbfOnly)));
  return result;
}



/******************************************************************************
PURPOSE: writeMaskedSubsetDBF - Write masked subset dbf file.
INPUTS:  const char* const baseInputFileName  Name of input dbf/shp file.
         const char* const csvFileName        0 or csv file/dir name.
         const int dbfOnly                    Write only DBF file?
         const long long huc                  0 or 14-digit HUC to filter by.
         const int comid                      0 or COMID to filter by.
         const char* const estcode            ESTCODE to filter by.
         const Bounds subset                  Subset lon-lat bounds.
         const int yyyymmdd1                  0 or start date to subset by.
         const int yyyymmdd2                  0 or end date to subset by.
         char* const mask                     mask per row (all 1).
OUTPUTS: char* const mask                     mask per row (0 or 1).
         DBFHandle outputDBF                  DBF file to write.
         int* wroteDBF                        Was DBF file written to?
RETURNS: int 1 if ok, else 0 and failure message is printed to stderr.
******************************************************************************/

static int writeMaskedSubsetDBF( const char* const baseInputFileName,
                                 const char* const csvFileName,
                                 const int dbfOnly,
                                 const long long huc,
                                 const int comid,
                                 const char* const estcode,
                                 const Bounds subset,
                                 const int yyyymmdd1,
                                 const int yyyymmdd2,
                                 char* const mask,
                                 DBFHandle outputDBF,
                                 int* wroteDBF ) {
  int ok = 1;
  const int rows = getRowsDBF( baseInputFileName );
  const int seagrass = strstr( baseInputFileName, "seagrass_" ) != 0;
  const int esat =
    OR2( strstr( baseInputFileName, "esat1_" ),
         strstr( baseInputFileName, "esat2_" ) );
  const int subsetFlowlines =
    OR2( strstr( baseInputFileName, "flowlines" ),
         strstr( baseInputFileName, "RBEROST" ) );
  const int isCSVDirectory = isDirectory( csvFileName );
  *wroteDBF = 0;

  DEBUG( fprintf( stderr, "huc = %lld, seagrass = %d, esat = %d, "
                  "subsetFlowlines = %d, "
                  "csvFileName = '%s', isCSVDirectory = %d, yyyymmdd1 = %d\n",
                  huc, seagrass, esat, subsetFlowlines, csvFileName,
                  isCSVDirectory,
                  yyyymmdd1 ); )

  if ( subsetFlowlines ) {
    ShapeData* shapeData = readDBF( baseInputFileName );
    ok = shapeData != 0;

    if ( shapeData ) {

      if ( comid > 0 ) { /* Subset by mask of rows upstream of comid: */
        int subsetRows = subsetByCOMID( shapeData, comid, mask );
        ok = subsetRows > 0;
        DEBUG( fprintf( stderr, "subsetByCOMID() = %d\n", subsetRows ); )
        deallocateShapeData( shapeData ), shapeData = 0;

        if ( ok ) {
          const Bounds global = { { -180.0, 180.0 }, { -90.0, 90.0 } };
          subsetRows =
            writeSubsetDBF( baseInputFileName, global, 0, 0, 0, 0, rows, mask,
                            outputDBF );
          ok = subsetRows > 0;
          *wroteDBF = ok;
          DEBUG( fprintf( stderr, "comid = %d, writeSubsetDBF() = %d\n",
                          comid, subsetRows ); )
        }
      } else if ( estcode ) {
        const int subsetRows =
          writeSubsetDBF( baseInputFileName, (const double (*)[2]) subset,
                          0, estcode, 0, 0, rows, mask, outputDBF );
        ok = subsetRows > 0;
        *wroteDBF = ok;
        DEBUG( fprintf( stderr, "estcode = '%s', writeSubsetDBF() = %d\n",
                        estcode, subsetRows );)
      }
    }
  } else if ( seagrass ) {
     const int subsetRows =
      writeSubsetDBF( baseInputFileName, (const double (*)[2]) subset,
                      0, estcode, 0, 0, rows, mask, outputDBF );
    ok = subsetRows > 0;
    *wroteDBF = ok;
    DEBUG( fprintf( stderr, "writeSubsetDBF() = %d\n", subsetRows );)
  } else if ( AND2( esat, csvFileName ) ) {

    /* Read the csv file and get a list of unique siteIds to filter DBF by: */

    size_t length = 0;
    size_t lines = 0;
    char* csvData = readFile( csvFileName, &length, &lines );
    ok = 0;

    if ( csvData ) {

      if ( lines < INT_MAX ) {
        int* siteIds = malloc( lines * sizeof (int) );

        if ( ! siteIds ) {
          fprintf( stderr, "\nFailed to allocate enough memory "
                   "to complete the requested operation.\n" );
        } else {
          const int uniqueCount = parseUniqueSiteIds( csvData, lines, siteIds);
          DEBUG( fprintf( stderr, "uniqueCount = %d\n", uniqueCount );)

          if ( uniqueCount > 0 ) {
            const int subsetRows =
              writeSubsetDBF( baseInputFileName, (const double (*)[2]) subset,
                              0,0, uniqueCount, siteIds, rows, mask, outputDBF);
            ok = subsetRows > 0;
            *wroteDBF = ok;
            DEBUG( fprintf( stderr, "writeSubsetDBF() = %d\n", subsetRows );)
          }
        }

        free( siteIds ), siteIds = 0;
      }

      free( csvData ), csvData = 0;
    }
  } else if ( AND3( dbfOnly, yyyymmdd1, ! csvFileName ) ) {
    const int subsetRows =
      subsetDBFByTimeAndBoundsOrEstcode( baseInputFileName,
                                         yyyymmdd1, yyyymmdd2,
                                         (const double (*)[2]) subset,
                                         estcode, rows, mask );
    ok = subsetRows > 0;
    DEBUG( fprintf( stderr, "%d = subsetDBFByTimeAndBoundsOrEstcode()\n",
                    subsetRows ); )

    if ( ok ) {
      const char* const estcode2 = estcode ? estcode : "all";
      const int subsetRowsWritten =
        writeSubsetDBF( baseInputFileName, (const double (*)[2]) subset,
                        0, estcode2, 0, 0, rows, mask, outputDBF );
      ok = subsetRowsWritten == subsetRows;
      *wroteDBF = ok;
      DEBUG( fprintf( stderr, "writeSubsetDBF() = %d\n", subsetRowsWritten ); )
    }
  } else if ( OR2( dbfOnly, csvFileName ) ) {

    /* Subset by huc/estcode/bounds: */

    if ( AND3( huc == 0, isCSVDirectory,
               OR2( ! estcode, ! strcmp( estcode, "all" ) ) ) ) {
      const int subsetShapes =
        computeShapeSubsetBoundsMask( baseInputFileName,
                                      (const double (*)[2]) subset,
                                      rows, mask );
      ok = subsetShapes > 0;
      DEBUG( fprintf( stderr, "subsetShapes = %d\n", subsetShapes ); )
    }

    if ( ok ) {
      const char* const estcode2 =
        ! isCSVDirectory ? 0
        :  estcode ? estcode
        : "all";
      const int subsetRows =
        writeSubsetDBF( baseInputFileName, (const double (*)[2]) subset,
                        huc, estcode2, 0, 0, rows, mask, outputDBF );
      ok = subsetRows > 0;
      *wroteDBF = ok;
      DEBUG( fprintf( stderr,
                      "estcode2 = '%s', writeSubsetDBF() = %d\n",
                      estcode2 ? estcode2 : "", subsetRows ); )
    }
  } else if ( yyyymmdd1 ) { /* Subset by time: */
    const int subsetRows =
      subsetDBFByTime( baseInputFileName, yyyymmdd1, yyyymmdd2, rows,
                       mask );
    ok = subsetRows > 0;
    DEBUG( fprintf( stderr, "%d = subsetDBFByTime()\n", subsetRows ); )
  } else if ( huc ) { /* Filter by HUC. */
    const Bounds global = { { -180.0, 180.0 }, { -90.0, 90.0 } };
    const int subsetRows =
      writeSubsetDBF( baseInputFileName, global, huc, 0, 0, 0, rows, mask,
                      outputDBF );
    ok = subsetRows > 0;
    *wroteDBF = ok;
    DEBUG( fprintf( stderr,
                    "huc = %lld, writeSubsetDBF() = %d subsetRows, ok = %d\n",
                    huc, subsetRows, ok ); )
  }

  DEBUG( fprintf( stderr,
                  "writeMaskedSubsetDBF(), returning %d, wroteDBF = %d\n",
                  ok, *wroteDBF ); )
  return ok;
}



/******************************************************************************
PURPOSE: writeMaskedSubsetSHP - Write masked subset shp (and maybe dbf) file.
INPUTS:  const char* const baseInputFileName  Name of input dbf/shp file.
         const char* const baseOutputFileName  Name of output dbf/shp file.
         const double minimumAdjacentVertexDistance  For polygons.
         const int clipShapes                 1 = clip shapes by mask & bounds,
                                              0 = just copy shapes per mask.
         const Bounds subset                  Subset lon-lat bounds.
         char* const mask                     mask per row (all 1).
         int* const offset                    Offset row to start writing.
         DBFHandle outputDBF                  DBF file to possibly write.
OUTPUTS: char* const mask                     mask per row (0 or 1).
         int* const offset                    Advanced by result.
         SHPHandle* const outputSHP           SHP file to write.
RETURNS: int polygonCount written if ok, else 0 and failure message is printed.
******************************************************************************/

static int writeMaskedSubsetSHP( const char* const baseInputFileName,
                                 const char* const baseOutputFileName,
                                 const double minimumAdjacentVertexDistance,
                                 const int clipShapes,
                                 const Bounds subset,
                                 char* const mask,
                                 int* const offset,
                                 DBFHandle outputDBF,
                                 SHPHandle* const outputSHP ) {
  int result = 0;
  int ok = 0;
  const int rows = getRowsDBF( baseInputFileName );

  DEBUG( fprintf( stderr, "writeMaskedSubsetSHP(): clipShapes = %d\n",
                  clipShapes ); )

  if ( ! clipShapes ) {

    /* Write unclipped shapes per mask. */

    result = copyMaskedShapes( baseInputFileName, baseOutputFileName,
                               outputSHP, rows, mask );
    DEBUG( fprintf( stderr, "  copyMaskedShapes() result = %d\n", result ); )
  } else { /* Write clipped unmasked shapes to subset: */
    int polygonCount = 0;
    int isPolyline = 0;
    PolygonShape* polygons =
      readAndClipShapes( baseInputFileName,
                         (const double (*)[2]) subset,
                         minimumAdjacentVertexDistance, mask,
                         &polygonCount, &isPolyline );

    DEBUG( fprintf( stderr, "  readAndClipShapes() polygonCount = %d\n",
                    polygonCount ); )

    if ( polygons ) {

      if ( *outputSHP == 0 ) {
        const int type = isPolyline ? SHPT_ARC : SHPT_POLYGON;
        *outputSHP = SHPCreate( baseOutputFileName, type );

        if ( *outputSHP == 0 ) {
          fprintf(stderr, "Failed to create temporary shp file '%s'.\n",
                          baseOutputFileName );
        }
      }

      if ( *outputSHP ) {
        ok = writePolygonsToShapefile( *outputSHP, isPolyline,
                                       polygonCount, polygons );

        if ( ok ) { /* Write filtered DBF based on clipped polygon id/area: */
          ok = writePolygonDBF( baseInputFileName, outputDBF, *offset,
                                polygonCount, mask, polygons );

          DEBUG( fprintf( stderr, "  writePolygonDBF() ok = %d\n", ok ); )

          if ( ok ) {
            result = polygonCount;
          }
        }
      }

      deallocatePolygons( polygonCount, polygons ), polygons = 0;
    }
  }

  *offset += result;
  DEBUG( fprintf( stderr, "writeMaskedSubsetSHP(), returning %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: writeMaskedSubsetCSV - Write masked subset csv file.
INPUTS:  const char* const baseInputFileName  Name of input dbf/shp file.
         const char* const baseOutputFileName  Name of output dbf/shp file.
         const char* const csvFileName        0 or csv file/dir name.
         const char* const estcode            ESTCODE to filter by.
         const int yyyymmdd1                  0 or start date to subset by.
         const int yyyymmdd2                  0 or end date to subset by.
         char* const mask                     mask per row.
OUTPUTS: char* const mask                     mask per row.
RETURNS: int values written if ok, else 0.
******************************************************************************/

static int writeMaskedSubsetCSV( const char* const baseInputFileName,
                                 const char* const baseOutputFileName,
                                 const char* const csvFileName,
                                 const char* const estcode,
                                 const int yyyymmdd1,
                                 const int yyyymmdd2,
                                 char* const mask ) {
  int result = 0;
  const int isSubsetFlowlines =
    OR2( strstr( baseInputFileName, "flowlines" ),
         strstr( baseInputFileName, "RBEROST" ) );
  const int isCSVDirectory = isDirectory( csvFileName );

  DEBUG( fprintf( stderr, "isSubsetFlowlines = %d, isCSVDirectory = %d\n",
                  isSubsetFlowlines, isCSVDirectory ); )

  if ( AND2( isSubsetFlowlines, csvFileName ) ) {

    /* Write subset of csv lines matching named id column in dbf: */

    const char* const columnName = "COMID";
    const int allowEmptyOutputCSV = 1;
    char outputCSVFileName[ 256 ] = "";
    memset( outputCSVFileName, 0, sizeof outputCSVFileName );
    strncpy( outputCSVFileName, baseOutputFileName,
             sizeof outputCSVFileName / sizeof *outputCSVFileName - 5 );
    strcat( outputCSVFileName, ".csv" );
    result =
      writeSubsetCSVFileById( baseInputFileName, csvFileName,
                              outputCSVFileName, columnName,
                              allowEmptyOutputCSV,
                              mask );
    DEBUG( fprintf( stderr, "writeSubsetCSVFileById() = %d\n", result ); )
  } else if ( AND2( isCSVDirectory, yyyymmdd1 ) ) {

    /* Write masked time-subset csv: */

    const char* const inputDBFFileName = baseInputFileName;
    const char* const inputCSVDirectory = csvFileName;
    char outputCSVFileName[ 256 ] = "";
    memset( outputCSVFileName, 0, sizeof outputCSVFileName );
    strncpy( outputCSVFileName, baseOutputFileName,
             sizeof outputCSVFileName / sizeof *outputCSVFileName - 5 );
    strcat( outputCSVFileName, ".csv" );
    result =
      writeSubsetCSVFile( inputDBFFileName, inputCSVDirectory,
                          outputCSVFileName,
                          yyyymmdd1, yyyymmdd2,
                          AND2( estcode, strcmp( estcode, "any" ) ) ?
                            estcode : 0,
                          mask );
    DEBUG( fprintf( stderr, "writeSubsetCSVFile() = %d\n", result ); )
  }

  DEBUG( fprintf( stderr, "writeMaskedSubsetCSV(), returning %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: parseUniqueSiteIds - Parse csv data for unique site ids.
INPUTS:  const char* const csvData  Content of a csv file to parse.
         const int lines            Lines of csv data (including header line).
OUTPUTS: int sideIds[ lines ]       Sorted array of unique site ids.
RETURNS: int number of unique site ids stored in sideIds[].
******************************************************************************/

static int parseUniqueSiteIds( const char* csvData,
                               const int lines,
                               int siteIds[] ) {

  PRE03( csvData, lines > 0, siteIds );
  int result = 0;
  const char newline = '\n';
  const int column = getSiteIdColumn( csvData );
  memset( siteIds, 0, lines * sizeof siteIds[0] );

  DEBUG( fprintf( stderr, "parseUniqueSiteIds()... column = %d, lines = %d\n",
                  column, lines ); )

  if ( column > 0 ) {

    /* Parse and store all non-0 ids: */

    const char* line = strchr( csvData, newline ); /* Skip header. */
    int previousId = 0;

    while ( AND2( result < lines, line ) ) {
      const char* const value = getColumnValue( line, column );

      if ( value ) {
        const int id = atoi( value );

        if ( AND2( id > 0, id != previousId ) ) {
          siteIds[ result ] = id;
          ++result;
          previousId = id;
        }
      }

      line = strchr( line + 1, newline );
    }

    /* Sort and uniq: */

    DEBUG( fprintf( stderr, "  temp result = %d\n", result ); )

    if ( result > 1 ) {
      int output = 0;
      int input = 1;
      DEBUG( fprintf( stderr, "  qsort()...\n" ); )
      qsort( siteIds, result, sizeof siteIds[0], intCompare );
      DEBUG( fprintf( stderr, "  uniq...\n" ); )

      while ( input < result ) {

        if ( siteIds[ input ] > siteIds[ output ] ) {
          ++output;
          siteIds[ output ] = siteIds[ input ];
        }

        ++input;
      }

      ++output;
      result = output;

      DEBUG( fprintf( stderr, "  uniq result = %d\n", result ); )

      /* Zero any trailing values: */

      for ( ; output < lines; ++output ) {
        siteIds[ output ] = 0;
      }
    }
  }

  DEBUG( fprintf( stderr, "  result = %d\n", result ); )
  POST02( result < lines,
          IMPLIES( result > 1, siteIds[ result - 1 ] > siteIds[ result - 2 ]));
  return result;
}



/******************************************************************************
PURPOSE: getSiteIdColumn - Parse csv data for 0-based index of siteId column.
INPUTS:  const char* csvData   Content of a csv file to parse.
RETURNS: int 0-based index of site id column.
******************************************************************************/

static int getSiteIdColumn( const char* csv_data ) {
  PRE0( csv_data );
  int result = 0;
  char header[ 4096 ] = "";
  const char* s = 0;
  memset( header, 0, sizeof header );
  strncpy( header, csv_data, sizeof header / sizeof *header - 1 );
  s = strstr( header, ",Site_Id(" );

  if ( ! s ) {
    s = strstr( header, ",SITE_ID(" );
  }

  if ( ! s ) {
    s = strstr( header, ",site_id(" );
  }

  if ( ! s ) {
    s = strstr( header, ",SiteId(" );
  }

  if ( s ) {
    const char* c = header;
    result = 1;

    while ( c != s ) {

      if ( *c == ',' ) {
        ++result;
      }

      ++c;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: getColumnValue - Parse csv data for value of specified column.
INPUTS:  const char* line  Line to parse.
RETURNS: const char* value of specified column.
******************************************************************************/

static const char* getColumnValue( const char* line, const int column ) {
  const char* result = 0;
  const char* next = line;
  int count = 0;

  do {
    next = strchr( next, ',' );

    if ( next ) {
      ++next;
      ++count;
    }

  } while ( next && count < column );

  if ( count == column ) {
    result = next;
  }

  return result;
}


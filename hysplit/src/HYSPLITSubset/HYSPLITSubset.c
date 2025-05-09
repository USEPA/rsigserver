
/******************************************************************************
PURPOSE: HYSPLITSubset.c - Read a subset of a HYSPLIT PM25 files
         and write it to stdout as one ASCII tab-delimited line.

NOTES:   To compile: ./makeit.
         To run:
         ../../../bin/$platform/HYSPLITSubset \
           -files test/filelist \
           -timestamp 2008022800 -hours 72 \
           -domain -76 35 -75 40 -scale 1e9

HISTORY: 2017-01-25 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For malloc(), free(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For strcmp(). */
#include <float.h>  /* For DBL_MAX. */

#include <Utilities.h> /* For PRE*(), FREE(), readFile(), etc. */

/*================================== TYPES ==================================*/

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;   /* File containing list of files to read. */
  double      scale;      /* Scale to multiply data by. */
  Bounds      bounds;     /* bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]*/
  int         yyyymmddhh; /* Starting timestamp of subset. */
  int         hours;      /* hours in subset. */
} Arguments;

#ifndef NO_ASSERTIONS

static int isValidArguments( const Arguments* arguments ) {
  const int result =
    AND9( arguments,
          arguments->listFile,
          arguments->listFile[ 0 ],
          IN_RANGE( arguments->scale, -DBL_MAX, DBL_MAX ),
          arguments->scale != 0.0,
          isValidBounds( arguments->bounds ),
          isValidYearMonthDay( arguments->yyyymmddhh / 100 ),
          IN_RANGE( arguments->yyyymmddhh % 100, 0, 23 ),
          arguments->hours > 0 );
  POST0( IS_BOOL( result ) );
  return result;
}

#endif /* NO_ASSERTIONS */

/* Subsetted/scaled file data: */

typedef struct {
  double* values;    /* Hourly PM25 values. */
  double* distances; /* Manhatten distance from center of arguments->bounds. */
  size_t count;      /* Length of above arrays == arguments->hours. */
} SubsetData;

static int initializeSubsetData( const int hours,
                                 SubsetData* const subsetData ) {
  PRE02( hours > 0, subsetData );
  int result = 0;
  subsetData->distances = 0;
  subsetData->count = 0;
  subsetData->values = NEW_ZERO( double, hours * 2 );

  if ( subsetData->values ) {
    subsetData->distances = subsetData->values + hours;
    subsetData->count = hours;
    result = 1;

    {
      double* distances = subsetData->distances;
      size_t hour = 0;

      for ( hour = 0; hour < (size_t) hours; ++hour ) {
        distances[ hour ] = DBL_MAX;
      }
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        AND7( subsetData->count > 0,
                              subsetData->values,
                              subsetData->distances ==
                                subsetData->values + subsetData->count,
                              subsetData->values[ 0 ] == 0.0,
                              subsetData->values[subsetData->count - 1] == 0.0,
                              subsetData->distances[ 0 ] == DBL_MAX,
                              subsetData->distances[ subsetData->count - 1 ] ==
                                DBL_MAX ),
                        IS_ZERO3( subsetData->values,
                                  subsetData->distances,
                                  subsetData->count ) ) );
  return result;
}

static void deallocateSubsetData( SubsetData* subsetData ) {
  PRE0( subsetData );
  FREE( subsetData->values );
  /* FREE( subsetData->distances ); distance is a pointer into values. */
  ZERO_OBJECT( subsetData );
  POST0( IS_ZERO3( subsetData->values,
                   subsetData->distances,
                   subsetData->count ) );
}

#ifndef NO_ASSERTIONS

static int isValidSubsetData( const SubsetData* subsetData ) {
  const int result =
  AND5( subsetData,
        subsetData->count,
        subsetData->values,
        subsetData->distances,
        subsetData->distances == subsetData->values + subsetData->count );
  POST0( IS_BOOL( result ) );
  return result;

}

#endif /* NO_ASSERTIONS */


/*========================== FORWARD DECLARATIONS ===========================*/

static void printUsage( const char* programName );

static int parseArguments( int argc, char* argv[], Arguments* arguments );

static int parseTimestampHours( const int arg, char* argv[],
                                int* yyyymmddhh, int* hours );

static int readData( const Arguments* const arguments,
                     SubsetData* const subsetData );

static char* readDataFile( const char* const fileName, size_t* lines );

static void processData( char* const fileData, const size_t lines,
                         const Arguments* const arguments,
                         SubsetData* const subsetData );

static size_t findMatchedLine( char* const fileData, const size_t lines,
                               const size_t line, const int yyyymmddhh,
                               const Bounds bounds,
                               double* const longitude,
                               double* const latitude,
                               double* const value );

static size_t lineLength( const char* const string );

static void writeData( const SubsetData* const subsetData );

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
    Arguments arguments;
    ZERO_OBJECT( &arguments );
    ok = parseArguments( argc, argv, &arguments );

    if ( ok ) {
      SubsetData subsetData;
      ok = initializeSubsetData( arguments.hours, &subsetData );
      
      if ( ok ) {
        ok = readData( &arguments, &subsetData );

        if ( ok ) {
          writeData( &subsetData );
        }
      }

      deallocateSubsetData( &subsetData );
    }
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
  fprintf( stderr, "\n\n\n%s - Read a subset of HYSPLIT PM25 files\n",
           programName );
  fprintf( stderr, "and write it to stdout in ASCII format.\n" );
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files <file_name> \\\n" );
  fprintf( stderr, "-timestamp <yyyymmddhh> -hours <hours> \\\n" );
  fprintf( stderr, "-domain <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, "[-scale <number>]\\\n" );
  fprintf( stderr, "Note: timestamp is in UTC (GMT)\n" );
  fprintf( stderr, "\nExample 1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-files test/file_list \\\n" );
  fprintf( stderr, "-timestamp 2008022800 -hours 12 \\\n" );
  fprintf( stderr, "-domain -111 31 -110 32 -scale 1e9\n" );
  fprintf( stderr, "\nPrints tab-delimited hourly PM25 values (ug/m3) \n");
  fprintf( stderr, "of nearest HYSPLIT point to center of domain:\n" );
  fprintf( stderr, "0\t0\t0\t0\t0\t0\t0\t0.05425\t0.02862\t0.08325\t0.04003\t0.004173\n" );
  fprintf( stderr, "\n\n\n");
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  int argc              Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: int 1 if successful, else 0.
NOTES:   Command-line options look like this:
         HYSPLITSubset
         -files /data/tmp/HYSPLITserver.32345
         -time 2008022800 -hours 72
         -domain -76 35 -75 36
         -scale 1e9
******************************************************************************/

static int parseArguments( int argc, char* argv[], Arguments* arguments ) {

  PRE02( isValidArgs( argc, (const char**) argv ), arguments );

  int result = 0;

  checkForTest( &argc, argv );
  ZERO_OBJECT( arguments );
  arguments->scale = 1.0;

  if ( IN_RANGE( argc, 12, 14 ) ) {
    int arg = 1;

    if ( ! strcmp( argv[ arg ], "-files" ) ) {
      ++arg;
      arguments->listFile = (const char*) argv[ arg ];
      ++arg;

      if ( parseTimestampHours( arg, argv,
                                &arguments->yyyymmddhh,
                                &arguments->hours ) ) {
        arg += 4;

        if ( ! strcmp( argv[ arg ], "-domain" ) ) {
          Integer skip = arg;
          result = parseBounds( argc, argv, &skip, arguments->bounds );
          arg = skip;
          
          /* Parse optional arguments: */

          if ( AND2( ! strcmp( argv[ arg ], "-scale" ), arg + 1 < argc ) ) {
            ++arg;
            arguments->scale = atof( argv[ arg ] );

            if ( OR2( arguments->scale == 0.0,
                      ! IN_RANGE( arguments->scale, -DBL_MAX, DBL_MAX ) ) ) {
              failureMessage( "Invalid scale: '%s'.", argv[ arg ] );
              result = 0;
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
PURPOSE: parseTimestampHours - Parse timestamp and hours options.
INPUTS:  const int arg   0-based index of -timestamp option in argv[].
         char* argv[]    Command line options.
OUTPUTS: int* yyyydddhh  Timestamp.
         int* hours      Number of hours.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int parseTimestampHours( const int arg, char* argv[],
                                int* yyyymmddhh, int* hours ) {

  PRE04( arg > 0, argv, yyyymmddhh, hours );
  int result = 0;
  *yyyymmddhh = *hours = 0;

  if ( AND5( ! strcmp( argv[ arg ], "-timestamp" ),
             argv[ arg + 1 ],
             strlen( argv[ arg + 1 ] ) == 10,
             ! strcmp( argv[ arg + 2 ], "-hours" ),
             argv[ arg + 3 ] ) ) {
    *yyyymmddhh = atoi( argv[ arg + 1 ] );
    *hours = atoi( argv[ arg + 3 ] );
    result =
      AND3( isValidYearMonthDay( *yyyymmddhh / 100 ),
            IN_RANGE( *yyyymmddhh % 100, 0, 23 ),
            *hours > 0 );

    if ( ! result ) {
      fprintf( stderr, "\a\n\nInvalid timestamp '%s' hours '%s' specified.\n",
               argv[ arg + 1 ], argv[ arg + 3 ] );
      *yyyymmddhh = *hours = 0;
    }
  }

  POST02( IS_BOOL( result ),
         IMPLIES_ELSE( result,
                       AND3( isValidYearMonthDay( *yyyymmddhh / 100 ),
                             IN_RANGE( *yyyymmddhh % 100, 0, 23 ),
                             *hours > 0 ),
                       IS_ZERO2( *yyyymmddhh, *hours ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readData - Read data files and extract subset data for output.
INPUTS:  const Arguments* const arguments  Command-line options.
OUTPUTS: SubsetData* const subsetData      Extracted subset data values.
RETURNS: int 1 if successful else 0.
******************************************************************************/

static int readData( const Arguments* const arguments,
                     SubsetData* const subsetData ) {

  PRE02( isValidArguments( arguments ), isValidSubsetData( subsetData ) );
  int result = 0;
  FILE* listFile = fopen( arguments->listFile, "r" );

  if ( listFile ) {
    char fileName[ 256 ] = "";
    memset( fileName, 0, sizeof fileName );

    while ( fgets( fileName, sizeof fileName / sizeof *fileName, listFile ) ) {
      char* const newline = strchr( fileName, '\n' );

      if ( newline ) {
        *newline = '\0';
      }

      DEBUG( fprintf( stderr, "Reading data file: %s\n", fileName ); )

      {
        size_t lines = 0;
        char* fileData = readDataFile( fileName, &lines );

        if ( fileData ) {
          result = 1;
          processData( fileData, lines, arguments, subsetData );
          FREE( fileData );
        }
      }
    }

    fclose( listFile ), listFile = 0;
  }

  POST02( IS_BOOL( result ), isValidSubsetData( subsetData ) );
  return result;
}



/******************************************************************************
PURPOSE: readDataFile - Read named data file as an array of strings.
INPUTS:  const char* const fileName  Name of data file to read.
OUTPUTS: size_t* lines  Number of strings/lines in file.
RETURNS: char* Allocated array of string lines.
******************************************************************************/

static char* readDataFile( const char* const fileName, size_t* lines ) {
  PRE03( fileName, *fileName, lines );
  Integer length = 0;
  char* result = readFile( fileName, &length );
  *lines = 0;

  if ( result ) {
    *lines = (size_t) linesInString( result );
  }

  POST0( IMPLIES_ELSE( result, *lines, *lines == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: processData - Extract subset of file data into subsetData.
INPUTS:  char* const fileData              Array of strings/lines of data.
         const size_t lines                Number of lines in fileData[].
         const Arguments* const arguments  Command-line options.
OUTPUTS: char* const fileData              Array of 0-terminated strings/lines.
         SubsetData* const subsetData      Subset/scaled values.
******************************************************************************/

static void processData( char* const fileData, const size_t lines,
                         const Arguments* const arguments,
                         SubsetData* const subsetData ) {
  PRE04( fileData, lines,
         isValidArguments( arguments ), isValidSubsetData( subsetData ) );

  const double longitude0 = 0.5 *
    ( arguments->bounds[ LONGITUDE ][ MINIMUM ] +
      arguments->bounds[ LONGITUDE ][ MAXIMUM ] );
  const double latitude0 = 0.5 *
    ( arguments->bounds[ LATITUDE ][ MINIMUM ] +
      arguments->bounds[ LATITUDE ][ MAXIMUM ] );
  double* const values = subsetData->values;
  double* const distances = subsetData->distances;
  const double scale = arguments->scale;
  const int hours = arguments->hours;
  int hour = 0;
  int yyyymmddhh = arguments->yyyymmddhh;
  size_t line0 = 1; /* Skip header line. */

  DEBUG( fprintf( stderr, "processData: lines = %lu, yyyymmddhh = %d, "
                  "hours = %d, longitude0 = %le, latitude0 = %le, "
                  "scale = %le\n",
                  lines, yyyymmddhh, hours, longitude0, latitude0, scale ); )

  for ( hour = 0; hour < hours; ++hour ) {
    size_t line = line0; /* Skip earlier lines. */
    double distance = distances[ hour ];

    DEBUG( fprintf( stderr, "  %d distance = %lg\n", yyyymmddhh, distance ); )

    if ( distance > 0.0 ) { /* It is possible that there is a closer point. */

      do {
        double longitude = 0.0;
        double latitude = 0.0;
        double value = 0.0;
        line = findMatchedLine( fileData, lines, line, yyyymmddhh,
                                arguments->bounds,
                                &longitude, &latitude, &value );

        if ( line < lines ) {
          const double longitudeDistance =
            longitude < longitude0 ? longitude0 - longitude
            : longitude - longitude0;
          CHECK( longitudeDistance >= 0.0 );

          if ( line0 == 1 ) {
            line0 = line; /* Index of subset hour 0. Subsequent hours are > */
          }

          if ( longitudeDistance < distance ) {
            const double latitudeDistance =
              latitude < latitude0 ? latitude0 - latitude
              : latitude - latitude0;
            const double thisDistance = longitudeDistance + latitudeDistance;
            CHECK3( latitudeDistance >= 0.0, thisDistance >= 0.0,
                    distance >= 0.0 );

            if ( thisDistance < distance ) {
              const double scaledValue = scale * value;
              values[ hour ] = scaledValue;
              distances[ hour ] = distance = thisDistance;
              DEBUG( fprintf( stderr, "  updated for hour = %d: "
                              "value = %le distance = %le\n",
                              hour, values[ hour ], distances[ hour ] ); )
            }
          }
        }

        ++line;
      } while ( line < lines );
    }

    DEBUG( fprintf( stderr, "  Final updated for hour = %d: "
                    "value = %le distance = %le\n",
                    hour, values[ hour ], distances[ hour ] ); )

    incrementHour( &yyyymmddhh );
  }

  POST02( isValidArguments( arguments ), isValidSubsetData( subsetData ) );
}



/******************************************************************************
PURPOSE: findMatchedLine - Search for the first line in fileData matching
         yyyymmddhh and within bounds and output its coordinates and value.
INPUTS:  char* const fileData              Array of strings/lines of data.
         const size_t lines                Number of lines in fileData[].
         const size_t line                 First line index to search from.
         const int yyyymmddhh              Timestamp to match.
         const Bounds bounds               Longitude-latitude bounds to match.
OUTPUTS: char* const fileData              Array of 0-terminated strings/lines.
         double* const longitude           Matched longitude or 0.
         double* const latitude            Matched latitude or 0.
         double* const value               Matched value or 0.
RETURNS: int index >= line of first matched line or else lines if no match.
NOTES:   The fileData has constant-length lines that look like:
YEAR, MO, DA, HR,     LAT,      LON,  PM2500100
2008,  1,  1,  7, 24.3000,-120.4000, 0.5017E-11
2008,  1,  1,  7, 24.4500,-120.4000, 0.5017E-11
...
******************************************************************************/

static size_t findMatchedLine( char* const fileData, const size_t lines,
                               const size_t line, const int yyyymmddhh,
                               const Bounds bounds,
                               double* const longitude,
                               double* const latitude,
                               double* const value ) {

  PRE09( fileData, lines, line < lines,
         isValidYearMonthDay( yyyymmddhh / 100 ),
         IN_RANGE( yyyymmddhh % 100, 0, 23 ),
         isValidBounds( bounds ),
         longitude, latitude, value );

  size_t result = lines; /* Default indicating 'not found'. */
  const size_t eachLineLength = lineLength( fileData ); /* Assume constant! */
  char* dataLine = fileData + line * eachLineLength;
  const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const int yyyy0 = yyyymmddhh / 1000000;
  const int mm0   = yyyymmddhh / 10000 % 100;
  const int dd0   = yyyymmddhh / 100 % 100;
  const int hh0   = yyyymmddhh % 100;
  size_t index = 0;

  DEBUG( fprintf( stderr, "    findMatchedLine(): "
                  "line   = %lu, eachLineLength = %lu\n",
                  line, eachLineLength ); )

  for ( index = line; index < lines; ++index, dataLine += eachLineLength ) {

    CHECK( dataLine <= fileData + ( lines - 1 ) * eachLineLength );

    /*
     * This is the most time-consuming section of code in this program
     * (because it is the inner-most loop, 16 seconds, 80% of total runtime),
     * so optimize it as follows:
     * 1. Terminate the line string for faster parsing by atoi/f()
     *    (apparently some implementations call strlen()).
     * 2. Assume fixed locations of commas as shown below
     *    012345678901234567890123456789012345678901234567
     *    2008,  1,  1,  7, 24.3000,-120.4000, 0.5017E-11
     * 3. Minimize the number of parsing calls to slow atoi/f().
     *    Note data lines are date-time-ordered so match each time component
     *    then match each coordinate.
     */

    dataLine[ eachLineLength - 1 ] = '\0';

    {
      const int yyyy = atoi( dataLine );

      if ( yyyy > yyyy0 ) {
        index = lines; /* Past target time so stop looping. */
      } else if ( yyyy == yyyy0 ) {
        const int mm = atoi( dataLine + 5 );

        if ( mm > mm0 ) {
          index = lines; /* Past target time so stop looping. */
        } else if ( mm == mm0 ) {
          const int dd = atoi( dataLine + 9 );

          if ( dd > dd0 ) {
            index = lines; /* Past target time so stop looping. */
          } else if ( dd == dd0 ) {
            const int hh = atoi( dataLine + 13 );

            if ( hh > hh0 ) {
              index = lines; /* Past target time so stop looping. */
            } else if ( hh == hh0 ) {
              *latitude = atof( dataLine + 17 );

              if ( IN_RANGE( *latitude, latitudeMinimum,  latitudeMaximum ) ) {
                *longitude = atof( dataLine + 26 );

                if ( IN_RANGE( *longitude,
                               longitudeMinimum, longitudeMaximum ) ) {
                  *value = atof( dataLine + 36 );

                  if ( IN_RANGE( *value, -DBL_MIN, DBL_MAX ) ) {
                    result = index; /* Found 1st time-matched point in bounds*/
                    index = lines; /* Stop looping. */
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if ( result == lines ) {
    *longitude = *latitude = *value = 0.0;
  }

  DEBUG( fprintf( stderr, "    findMatchedLine(): result = %lu\n", result ); )

  POST0( IMPLIES_ELSE( IN_RANGE( result, line, lines - 1 ),
                       AND4( isValidLongitudeLatitude( *longitude, *latitude ),
                             IN_RANGE( *longitude,
                                       bounds[ LONGITUDE ][ MINIMUM ],
                                       bounds[ LONGITUDE ][ MAXIMUM ] ),
                             IN_RANGE( *latitude,
                                       bounds[ LATITUDE ][ MINIMUM ],
                                       bounds[ LATITUDE ][ MAXIMUM ] ),
                             IN_RANGE( *value, -DBL_MIN, DBL_MAX ) ),
                       AND2( result == lines,
                             IS_ZERO3( *longitude, *latitude, *value ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: lineLength - Length of line string (including \n).
INPUTS:  const char* const string  String to check.
RETURNS: size_t length of line.
******************************************************************************/

static size_t lineLength( const char* const string ) {
  PRE0( string );
  const char* s = string;
  size_t result = 0;

  while ( AND2( *s, *s != '\n' ) ) {
    ++result;
    ++s;
  }

  result += *s == '\n';
  return result;
}



/******************************************************************************
PURPOSE: writeData - Write subset of file data to stdout.
INPUTS:  const SubsetData* const subsetData  Array of values to write.
******************************************************************************/

static void writeData( const SubsetData* const subsetData ) {
  PRE0( isValidSubsetData( subsetData ) );
  const size_t count = subsetData->count;
  const double* const values = subsetData->values;
  size_t index = 1;

  printf( "%g", values[ 0 ] );

  for ( index = 1; index < count; ++index ) {
    const double value = values[ index ];
    printf( "\t%g", value );
  }

  printf( "\n" );
  POST0( isValidSubsetData( subsetData ) );
}



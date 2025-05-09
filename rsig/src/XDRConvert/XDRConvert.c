
/******************************************************************************
PURPOSE: XDRConvert.c - Read a stream of data in XDR-format and write it to
         stdout in various formats including: NetCDF COARDS,
         regridded NetCDF COARDS, regridded NetCDF IOAPI, regridded XDR,
         regridded ASCII.

NOTES:   For a description of NetCDF COARDS Conventions see:
         http://ferret.wrc.noaa.gov/noaa_coop/coop_cdf_profile.html

         To implement translation of a new input format it is not necessary to
         edit any existing routines, only insert new entries into the
         translators table (data) and declare and define a new routine
         in the appropriate .h and .c files and include the .h file below.
         Just follow the steps documented in the code. (Search for STEP)
         Then add a test input file and update runit and testit scripts.

HISTORY: 2007-12 plessel.todd@epa.gov

STATUS:  unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h>     /* For stderr, fprintf(). */
#endif
#include <string.h>    /* For strcmp(). */
#include <unistd.h>    /* For unlink(), getpid(). */
#include <sys/stat.h>  /* For struct stat, stat(). */

#include <Utilities.h>  /* For PRE*(), NEW_ZERO(), Integer, Real, Stream. */

#include <Helpers.h>    /* For COUNT(), indexOfString(). */
#include <Parameters.h> /* For Parameters, isValidParameters(). */

/*================================== TYPES ==================================*/

typedef void (*Translator)( Parameters* );

typedef struct {
  const char* const tag; /* First line of input stream identifies its type. */
  Translator translator; /* Function to translate format. */
} Entry;

/*============================= GLOBAL VARIABLES ============================*/

extern void translateCMAQ(         Parameters* );
extern void translatePoint(        Parameters* );
extern void translateSite(         Parameters* );
extern void translateSwath(        Parameters* );
extern void translateCALIPSO(      Parameters* );
extern void translateAircraft(     Parameters* );
extern void translateProfile(      Parameters* );
extern void translateGrid(         Parameters* );
extern void compareRegriddedPoint( Parameters* );
extern void compareRegriddedSite(  Parameters* );
extern void compareRegriddedSwath( Parameters* );
extern void compareRegriddedCALIPSO(  Parameters* );
extern void compareRegriddedAircraft( Parameters* );
extern void compareRegriddedProfile(  Parameters* );
extern void compareRegriddedGrid(     Parameters* );

/* STEP 1. Declare new translator routine above this line. */

static const Entry translators[] = {
  { "SUBSET 9.0 CMAQ\n",        translateCMAQ            },
  { "Point 1.0\n",              translatePoint           },
  { "SITE 2.0\n",               translateSite            },
  { "Swath 2.0\n",              translateSwath           },
  { "CALIPSO 1.0\n",            translateCALIPSO         },
  { "Aircraft 2.0\n",           translateAircraft        },
  { "Profile 2.0\n",            translateProfile         },
  { "Grid 1.0\n",               translateGrid            },
  { "REGRIDDED-Point 1.0\n",    compareRegriddedPoint    },
  { "REGRIDDED-SITE 2.0\n",     compareRegriddedSite     },
  { "REGRIDDED-Swath 2.0\n",    compareRegriddedSwath    },
  { "REGRIDDED-CALIPSO 2.0\n",  compareRegriddedCALIPSO  },
  { "REGRIDDED-Aircraft 3.0\n", compareRegriddedAircraft },
  { "REGRIDDED-Profile 2.0\n",  compareRegriddedProfile  },
  { "REGRIDDED-Grid 1.0\n",     compareRegriddedGrid  },

/* STEP 2. Add new translator entry above this line. */

  { 0, 0 }
};

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateParameters( Parameters* parameters );

static void printUsage( const char* programName );

static void temporaryFileName( const char* const tmpdir,
                               const char* const tag,
                               char result[ 256 ] );

static void parseParameters( int argc, char* argv[],
                             Parameters* parameters );

static const char* guessTmpDir( void );

static void parseTmpdir( int argc, char* argv[], Integer* arg,
                         Parameters* parameters );

static void parseRegrid( int argc, char* argv[], Integer* arg,
                         Parameters* parameters );

static void parseFormat( int argc, char* argv[], Integer* arg,
                         Parameters* parameters );

static void parseAggregate( int argc, char* argv[], Integer* arg,
                            Parameters* parameters );

static void parseCompare( int argc, char* argv[], Integer* arg,
                          Parameters* parameters );

static Translator findTranslator( const char* line );

static void readCMAQXDR( const char* fileName, Parameters* parameters );

static void readCMAQXDRHeader( Stream* input, Parameters* parameters );

static void readCMAQXDRData( Stream* input, Parameters* parameters );

static Projector* parseProjectorFromXDRHeader( Stream* input, Integer* ok );

static Grid* parseGridFromXDRHeader( Stream* input, Projector* projector );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Read a stream of data in XDR-format and write it to
         stdout in various formats including: NetCDF COARDS,
         regridded NetCDF COARDS, regridded NetCDF IOAPI, regridded XDR,
         regridded ASCII.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Integer ok = 0;

  if ( ! isValidArgs( argc, (const char**) argv ) ) {
    failureMessage( "Invalid command-line arguments." );
  } else {
    Parameters parameters;
    ZERO_OBJECT( &parameters );
    checkForTest( &argc, argv ); /* Check for and remove any -test argument. */
    parseParameters( argc, argv, &parameters );

    if ( parameters.ok ) {
      parameters.input = newFileStream( "-stdin", "rb" );

      if ( parameters.input ) {

        if ( AND4( parameters.format == FORMAT_XDR,
                   parameters.regrid == 0,
                   parameters.compareFunction == 0,
                   parameters.convertFunction == 0 ) ) {
          ok = parameters.ok = copyToStdout( parameters.input );
        } else {
          char line[ 80 ] = ""; /* Read first line to determine input type. */
          parameters.input->readString( parameters.input, line, COUNT( line ));

          if ( parameters.input->ok( parameters.input ) ) {
            Translator translator = findTranslator( line ); /*For input type*/

            if ( translator ) {
              parameters.ok = 1;

              if ( IN3( parameters.format, FORMAT_COARDS, FORMAT_IOAPI ) ) {
                temporaryFileName( parameters.temporaryDirectory, "netcdf",
                                   parameters.netcdfFileName );
              }

              if ( parameters.regrid ) {
                  temporaryFileName( parameters.temporaryDirectory, "regrid",
                                     parameters.regridFileName );
              }

              if ( parameters.ok ) {
                DEBUG( fprintf( stderr, "Translating:\n" ); )
                translator( &parameters ); /* Call translator for input. */

                if ( AND2( parameters.ok, parameters.netcdfFileName[ 0 ] ) ) {

                  /* Read and copy temporary NetCDF file to stdout: */

                  DEBUG( fprintf( stderr, "Streaming:\n" ); )
                  parameters.ok = streamFile( parameters.netcdfFileName );
                }
              }

              ok = parameters.ok;
            }
          }
        }
      }
    }

    deallocateParameters( &parameters );
  }

  if ( AND2( ! ok, ! failureCount() ) ) {
    failureMessage( "No points in output." );
  }

  DEBUG( fprintf( stderr, "XDRConvert is returning %d\n", ! ok ); )

  return ! ok;
}



/******************************************************************************
PURPOSE: isValidParameters - Is parameters structure valid?
INPUTS:  const Parameters* parameters  Parameters to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

Integer isValidParameters( const Parameters* parameters ) {

  const Integer result =
    AND8( parameters,
          parameters->input,
          parameters->input->isReadable( parameters->input ),
          IS_VALID_FORMAT( parameters->format ),
          IMPLIES( parameters->regrid,
                   AND4( IS_VALID_AGGREGATE_METHOD( parameters->regrid ),
                         parameters->grid->invariant( parameters->grid ),
                         parameters->temporaryDirectory,
                         parameters->temporaryDirectory[ 0 ] ) ),
          IMPLIES_ELSE( OR2( parameters->format == FORMAT_COARDS,
                             parameters->format == FORMAT_IOAPI ),
                        parameters->netcdfFileName[ 0 ] != '\0',
                        parameters->netcdfFileName[ 0 ] == '\0' ),
          IS_BOOL( parameters->ok ),
          IMPLIES( OR2( parameters->compareFunction,
                        parameters->convertFunction ),
                   AND11( parameters->timestamp,
                          parameters->timesteps > 0,
                          parameters->aggregationTimesteps >= 0,
                          parameters->firstLayer > 0,
                          parameters->lastLayer >= parameters->firstLayer,
                          parameters->firstRow > 0,
                          parameters->lastRow >= parameters->firstRow,
                          parameters->firstColumn > 0,
                          parameters->lastColumn >= parameters->firstColumn,
                          parameters->data,
                          IMPLIES_ELSE( parameters->compareFunction,
                                        IS_ZERO2( parameters->convertFunction,
                                                  parameters->data2 ),
                                        AND2( parameters->convertFunction,
                                              parameters->data2 ) ) ) ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateParameters - Deallocate contents of parameters structure.
INPUTS:  Parameters* parameters  Parameters to deallocate contents of.
******************************************************************************/

static void deallocateParameters( Parameters* parameters ) {

  PRE0( parameters );

  FREE_OBJECT( parameters->input );
  FREE_OBJECT( parameters->grid );

#ifndef DEBUGGING

  if ( parameters->netcdfFileName[ 0 ] ) {
    unlink( parameters->netcdfFileName ); /* Remove temporary file. */
  }

  if ( parameters->regridFileName[ 0 ] ) {
    unlink( parameters->regridFileName ); /* Remove temporary file. */
  }

#endif

  FREE( parameters->data );
  ZERO_OBJECT( parameters );

  POST0( ! isValidParameters( parameters ) );
}



/******************************************************************************
PURPOSE: temporaryFileName - Generate and return a unique file name to be used
         for the temporary NetCDF file.
INPUTS:  const char* const tmpdir  Name of a writable directory for temp files.
         const char* const tag     Tag name for file.
OUTPUTS: char result[ 256 ]        Temporary file name.
******************************************************************************/

static void temporaryFileName( const char* const tmpdir,
                               const char* const tag,
                               char result[ 256 ] ) {
  PRE03( tag, *tag, result );
  const char* const ext = ! strcmp( tag, "netcdf" ) ? ".nc" : ".bin";
  const Integer pid = getpid();
  memset( result, 0, 256 * sizeof (char) );
  snprintf( result, 256 * sizeof (char),
            "%s/rsig_temp_%s_%lld%s", tmpdir, tag, pid, ext );
  POST02( result, *result );
}



/******************************************************************************
PURPOSE: printUsage - Print program printUsage instructions to stderr.
INPUTS:  const char* programName  Name of this program.
******************************************************************************/

static void printUsage( const char* programName ) {

  PRE02( programName, *programName );

  fprintf( stderr, "\n\n%s - Read a stream of XDR-format data from stdin\n",
           programName );
  fprintf( stderr, "and write it to stdout in various formats including:\n" );
  fprintf( stderr, "  XDR\n  ASCII\n  NetCDF-COARDS\n" );
  fprintf( stderr, "  Regridded XDR\n  Regridded ASCII\n" );
  fprintf( stderr, "  Regridded NetCDF-COARDS\n" );
  fprintf( stderr, "  Regridded NetCDF-IOAPI\n\n" );
  fprintf( stderr, "Usage: %s", programName );
  fprintf( stderr, " [-tmpdir temporary_directory ]\n" );
  fprintf( stderr, " [-regrid nearest | mean | weighted\n" );
  fprintf( stderr, "  [ -lambert lower_latitude upper_latitude" );
  fprintf( stderr, " central_longitude central_latitude |\n" );
  fprintf( stderr, "    -mercator central_longitude |\n" );
  fprintf( stderr, "    -stereographic central_longitude central_latitude" );
  fprintf( stderr, " secant_latitude \n" );
  fprintf( stderr, "    -lonlat ]\n" );
  fprintf( stderr, "  -ellipsoid major_semiaxis minor_semiaxis\n" );
  fprintf( stderr, "  -grid columns rows west_edge south_edge" );
  fprintf( stderr, " cell_width cell_height\n" );
  fprintf( stderr, "  [-layers layers type top_pressure" );
  fprintf( stderr, " level_1 level_2 ... level_layers+1 g R A T0s P00]]\n" );
  fprintf( stderr, " [-compare difference | absolute_difference | " );
  fprintf( stderr, "percent_difference | ratio | convert compare_file]\n" );
  fprintf( stderr, " [-aggregate timesteps ]\n" );
  fprintf( stderr, " -xdr | -ascii | -ioapi | -coards | -mcmc\n\n" );
  fprintf( stderr, "Note the following constants are from MM5:\n");
  fprintf( stderr, "g   = 9.81     Gravitational force m/s^2.\n" );
  fprintf( stderr, "R   = 287.04   Gas constant J/kg/K = m^3/s/K.\n" );
  fprintf( stderr, "A   = 50.0     Atmospheric lapse rate in K/kg.\n" );
  fprintf( stderr, "T0s = 290.0    Reference surface temperature in K.\n" );
  fprintf( stderr, "P00 = 100000.0 Reference surface pressure in Pa.\n" );
  fprintf( stderr, "\nexamples:\n\n" );
  fprintf( stderr, "  cat airnow.xdr | %s -coards", programName );
  fprintf( stderr, " > airnow.nc ; ncdump airnow.nc | more\n\n" );
  fprintf( stderr, "  cat modis.xdr | %s", programName );
  fprintf( stderr, " -regrid mean -lambert 33 45 -97 40" );
  fprintf( stderr, " -ellipsoid 6370000 6370000" );
  fprintf( stderr, " -grid 268 259 -420000 -1716000 12000 12000 -ioapi" );
  fprintf( stderr, " > modis.nc ; ncdump modis.nc | more\n\n" );
  fprintf( stderr, "  cat calipso.xdr | %s", programName );
  fprintf( stderr, " -regrid mean -lambert 33 45 -97 40" );
  fprintf( stderr, " -ellipsoid 6370000 6370000" );
  fprintf( stderr, " -grid 268 259 -420000 -1716000 12000 12000" );
  fprintf( stderr, " -layers 22 2 10000" );
  fprintf( stderr, " 1.0 0.995 0.988 0.979 0.97 0.96 0.938 0.914" );
  fprintf( stderr, " 0.889 0.862 0.834 0.804 0.774 0.743 0.694 0.644" );
  fprintf( stderr, " 0.592 0.502 0.408 0.311 0.21 0.106 0.0" );
  fprintf( stderr, " 9.81 287.04 50.0 290.0 100000.0" );
  fprintf( stderr, " -ioapi" );
  fprintf( stderr, " > calipso.nc ; ncdump calipso.nc | more\n\n" );
  fprintf( stderr, "  cat airnow.xdr | %s -coards", programName );
  fprintf( stderr, " > airnow.nc ; ncdump airnow.nc | more\n\n" );
  fprintf( stderr, "  cat modis.xdr | %s", programName );
  fprintf( stderr, " -regrid mean -stereographic -98 90 45" );
  fprintf( stderr, " -ellipsoid 6370000 6370000" );
  fprintf( stderr, " -grid 137 137 -7398000 -7398000 108000 108000 -ioapi" );
  fprintf( stderr, " -aggregate 24 > modis.nc ; ncdump modis.nc | more\n\n" );
}



/******************************************************************************
PURPOSE: parseParameters - Parse command-line arguments into parameters.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
OUTPUTS: Parameters* parameters  Updated parameters.
******************************************************************************/

static void parseParameters( int argc, char* argv[], Parameters* parameters ) {

  PRE02( isValidArgs( argc, (const char**) argv ), parameters );

  if ( argc < 2 ) {
    failureMessage( "Missing command-line arguments." );
  } else {
    Integer argument = 1; /* Index of argument to parse. */
    parameters->ok = 1;
    parseTmpdir( argc, argv, &argument, parameters );

    if ( parameters->ok ) {
      parseRegrid( argc, argv, &argument, parameters );

      if ( parameters->ok ) {
        parseAggregate( argc, argv, &argument, parameters );

        if ( parameters->ok ) {
          parseCompare( argc, argv, &argument, parameters );

          if ( parameters->ok ) {
            parseFormat( argc, argv, &argument, parameters );
          }
        }
      }
    }
  }

  if ( ! parameters->ok ) {
    printUsage( argv[ 0 ] );
  }

  POST02( IS_BOOL( parameters->ok ),
          IMPLIES( parameters->ok,
                   AND2( IS_VALID_FORMAT( parameters->format ),
                         IMPLIES( parameters->regrid,
                           AND2( IS_VALID_AGGREGATE_METHOD(parameters->regrid),
                                 parameters->grid->invariant(
                                   parameters->grid ) ) ) ) ) );
}



/******************************************************************************
PURPOSE: guessTmpdir - Guess a usable temporary default directory.
RETURNS: static const char* default name of temporary directory.
******************************************************************************/

static const char* guessTmpDir( void ) {
  static const char* const table[] = {
    "/data/tmp",
    "/data/rsig/tmp",
    "/local_proc",
    "/home/tplessel/RSIG/tmp",
    "testdata"
  };
  const size_t count = sizeof table / sizeof *table;
  size_t index = 0;
  const char* result = 0;
  struct stat unused_;

  for ( index = 0; index < count; ++index ) {
    const char* const name = table[ index ];

    if ( ! stat( name, &unused_ ) ) {
      result = name;
      index = count;
    }
  }

  if ( ! result ) {
    result = "/tmp";
  }

  POST0( stat( result, &unused_ ) == 0 );
  return result;
}



/******************************************************************************
PURPOSE: parseTmpdir - Parse tmpdir command-line option into parameters.
INPUTS:  int argc           Number of command-line arguments.
         char* argv[]       Command-line argument strings.
         Integer* argument  Index of argument to parse.
OUTPUTS: Integer* argument  Updated index of argument to parse.
         Parameters* parameters  Updated parameters->temporaryDirectory.
******************************************************************************/

static void parseTmpdir( int argc, char* argv[], Integer* argument,
                         Parameters* parameters ) {

  PRE05( isValidArgs( argc, (const char**) argv ),
         argument, *argument > 0,
         parameters, parameters->ok );

  struct stat unused_;

  if ( AND2( *argument + 2 < argc,
             ! strcmp( argv[ *argument ], "-tmpdir" ) ) ) {
    parameters->temporaryDirectory = argv[ *argument + 1 ];
    parameters->ok =
      AND3( parameters->temporaryDirectory,
            parameters->temporaryDirectory[ 0 ],
            stat( parameters->temporaryDirectory, &unused_ ) == 0 );

    if ( ! parameters->ok ) {
      failureMessage( "Invalid -tmpdir argument '%s'\n",
                      argv[ *argument + 1 ] );
    } else {
      *argument += 2;
    }
  } else {
    parameters->temporaryDirectory = guessTmpDir();
    parameters->ok = stat( parameters->temporaryDirectory, &unused_ ) == 0;

    if ( ! parameters->ok ) {
      failureMessage( "No valid temporary directory found.\n" );
    }
  }

  POST03( *argument <= argc,
          IS_BOOL( parameters->ok ),
          IMPLIES( parameters->ok,
                   AND2( parameters->temporaryDirectory,
                         stat(parameters->temporaryDirectory, &unused_) ==0)));
}



/******************************************************************************
PURPOSE: parseRegrid - Parse optional regrid command-line arguments into
         parameters.
INPUTS:  int argc           Number of command-line arguments.
         char* argv[]       Command-line argument strings.
         Integer* argument  Index of argument to parse.
OUTPUTS: Integer* argument  Updated index of argument to parse.
         Parameters* parameters  Updated parameters.
******************************************************************************/

static void parseRegrid( int argc, char* argv[], Integer* argument,
                         Parameters* parameters ) {

  PRE06( isValidArgs( argc, (const char**) argv ),
         argument, *argument > 0, *argument <= argc,
         parameters, parameters->ok );

  if ( AND2( *argument + 1 < argc,
             ! strcmp( argv[ *argument ], "-regrid" ) ) ) {
    const char* const method = argv[ *argument + 1 ];
    const char* const methods[] = { "nearest", "mean", "weighted" };
    const Integer count = COUNT( methods );
    const Integer index = indexOfString( method, methods, count );
    parameters->ok = index != -1;

    if ( ! parameters->ok ) {
      failureMessage( "Invalid regrid method '%s'\n", method );
    } else {
      int isLonLat = 0;
      Projector* projector = 0;
      *argument += 2;
      parameters->regrid = index + 1;

      if ( ! strcmp( argv[ *argument ], "-lonlat" ) ) {
        Real majorSemiaxis = 0.0;
        Real minorSemiaxis = 0.0;
        Integer ok = 0;
        isLonLat = 1;
        *argument += 1;
        parseEllipsoid( argc, argv, argument,
                        &majorSemiaxis, &minorSemiaxis,
                        &ok );
      } else {
        projector = parseProjection( argc, argv, argument );
        parameters->ok = projector != 0;
      }

      if ( OR2( isLonLat, projector ) ) {
        parameters->grid = parseGrid( argc, argv, argument, projector );
        parameters->ok = parameters->grid != 0;

        if ( parameters->grid ) {
          projector = 0;
        }
      }
    }
  }

  POST03( *argument <= argc,
          IS_BOOL( parameters->ok ),
          IMPLIES( parameters->ok,
                   IMPLIES( parameters->regrid,
                           AND2( IS_VALID_AGGREGATE_METHOD(parameters->regrid),
                                 parameters->grid->invariant(
                                   parameters->grid ) ) ) ) );
}



/******************************************************************************
PURPOSE: parseAggregate - Parse aggregate command-line arguments into
         parameters.
INPUTS:  int argc           Number of command-line arguments.
         char* argv[]       Command-line argument strings.
         Integer* argument  Index of argument to parse.
OUTPUTS: Integer* argument  Updated index of argument to parse.
         Parameters* parameters  Updated parameters.
******************************************************************************/

static void parseAggregate( int argc, char* argv[], Integer* argument,
                            Parameters* parameters ) {

  PRE05( isValidArgs( argc, (const char**) argv ),
         argument, *argument > 0,
         parameters, parameters->ok );

  if ( AND2( *argument + 2 < argc,
             ! strcmp( argv[ *argument ], "-aggregate" ) ) ) {
    parameters->aggregationTimesteps = atoi( argv[ *argument + 1 ] );
    parameters->ok = parameters->aggregationTimesteps >= 0;

    if ( ! parameters->ok ) {
      failureMessage( "Invalid -aggregate argument '%s'\n",
                      argv[ *argument + 1 ] );
    } else {
      *argument += 2;
    }
  }

  POST03( *argument <= argc,
          IS_BOOL( parameters->ok ),
          IMPLIES( parameters->ok, parameters->aggregationTimesteps >= 0 ) );
}



/******************************************************************************
PURPOSE: parseFormat - Parse format command-line arguments into
         parameters.
INPUTS:  int argc           Number of command-line arguments.
         char* argv[]       Command-line argument strings.
         Integer* argument  Index of argument to parse.
OUTPUTS: Integer* argument  Updated index of argument to parse.
         Parameters* parameters  Updated parameters.
******************************************************************************/

static void parseFormat( int argc, char* argv[], Integer* argument,
                         Parameters* parameters ) {

  PRE05( isValidArgs( argc, (const char**) argv ),
         argument, *argument > 0,
         parameters, parameters->ok );

  if ( *argument < argc ) {
    const char* const format = argv[ *argument ];
    const char* const formats[] = {
      "-xdr", "-ascii", "-coards", "-ioapi", "-mcmc"
    };
    const Integer count = COUNT( formats );
    const Integer index = indexOfString( format, formats, count );
    parameters->ok = index != -1;

    if ( ! parameters->ok ) {
      failureMessage( "Invalid format '%s'\n", format );
    } else {
      ++*argument;
      parameters->format = index;
    }
  } else {
    failureMessage( "Missing format argument\n" );
    parameters->ok = 0;
  }

  POST03( *argument <= argc,
          IS_BOOL( parameters->ok ),
          IMPLIES( parameters->ok,
                   IS_VALID_FORMAT( parameters->format ) ) );
}



/******************************************************************************
PURPOSE: parseCompare - Parse compare command-line arguments into
         parameters.
INPUTS:  int argc           Number of command-line arguments.
         char* argv[]       Command-line argument strings.
         Integer* argument  Index of argument to parse.
OUTPUTS: Integer* argument  Updated index of argument to parse.
         Parameters* parameters  Updated parameters.
******************************************************************************/

static void parseCompare( int argc, char* argv[], Integer* argument,
                          Parameters* parameters ) {

  PRE05( isValidArgs( argc, (const char**) argv ),
         argument, *argument > 0,
         parameters, parameters->ok );

  if ( AND2( *argument + 2 < argc,
             ! strcmp( argv[ *argument ], "-compare" ) ) ) {
    const char* const operator = argv[ *argument + 1 ];
    parameters->compareFunction = compareFunction( operator );
    parameters->convertFunction = convertFunction( operator );
    parameters->ok =
      OR2( parameters->compareFunction, parameters->convertFunction );

    if ( ! parameters->ok ) {
      failureMessage( "Invalid compare operator '%s'\n", operator );
    } else {
      readCMAQXDR( argv[ *argument + 2 ], parameters );
      *argument += 3;
    }
  }

  POST03( *argument <= argc,
          IS_BOOL( parameters->ok ),
          IMPLIES( AND2( parameters->ok,
                         OR2( parameters->compareFunction,
                              parameters->convertFunction ) ),
                   AND8( parameters->timestamp,
                         parameters->timesteps > 0,
                         parameters->firstRow > 0,
                         parameters->lastRow >= parameters->firstRow,
                         parameters->firstColumn > 0,
                         parameters->lastColumn >= parameters->firstColumn,
                         parameters->data,
                         IMPLIES( parameters->convertFunction,
                           AND2( parameters->compareFunction == 0,
                                 parameters->data2 ) ) ) ) );
}



/******************************************************************************
PURPOSE: findTranslator - Read a stream of data in XDR-format and write it to
         stdout in another format.
INPUTS:  const char* line First line of input used to determine translator.
RETURNS: Translator if successful, else 0.
******************************************************************************/

static Translator findTranslator( const char* line ) {

  PRE0( line );

  Translator result = 0;
  const Integer count = COUNT( translators ) - 1;
  Integer index = 0;

  do {

    if ( strcmp( line, translators[ index ].tag ) == 0 ) {
      result = translators[ index ].translator;
      index = count;
    }

    ++index;
  } while ( index < count );

  if ( ! result ) {
    failureMessage( "Invalid input data format '%s'.", line );
  }

  return result;
}



/******************************************************************************
PURPOSE: readCMAQXDR - Read CMAQ XDR-format file for comparing.
INPUTS:  const Char* const fileName CMAQ XDR-format file to read.
OUTPUTS: Parameters* parameters  Updated parameters->data.
******************************************************************************/

static void readCMAQXDR( const char* fileName, Parameters* parameters ) {

  PRE04( fileName, parameters, parameters->ok,
         OR2( parameters->compareFunction, parameters->convertFunction ) );

  Stream* input = newFileStream( fileName, "rb" );
  parameters->ok = input != 0;

  if ( input ) {
    readCMAQXDRHeader( input, parameters );

    if ( parameters->ok ) {
      readCMAQXDRData( input, parameters );
    }

    FREE_OBJECT( input );
  }

  POST02( IS_BOOL( parameters->ok ),
          IMPLIES( parameters->ok,
                   AND9( OR2( parameters->compareFunction,
                              parameters->convertFunction ),
                         parameters->timestamp,
                         parameters->timesteps > 0,
                         parameters->firstRow > 0,
                         parameters->lastRow >= parameters->firstRow,
                         parameters->firstColumn > 0,
                         parameters->lastColumn >= parameters->firstColumn,
                         parameters->data,
                         IMPLIES( parameters->convertFunction,
                           AND2( parameters->compareFunction == 0,
                                  parameters->data2 ) ) ) ) );

}



/******************************************************************************
PURPOSE: readCMAQXDRHeader - Read CMAQ XDR-format file header for comparing.
INPUTS:  Stream* input           Stream to read from.
OUTPUTS: Parameters* parameters  Updated parameters->timestamp...lastColumn.
******************************************************************************/

static void readCMAQXDRHeader( Stream* input, Parameters* parameters ) {

  PRE07( input, input->invariant( input ), input->isReadable( input ),
         input->ok( input ),
         parameters, parameters->ok,
         OR2( parameters->compareFunction, parameters->convertFunction ) );


  char line[ 1024 ] = "";
  const Integer size = COUNT( line );
  parameters->ok = 0;
  input->readString( input, line, size );

  if ( input->ok( input ) ) {

    if ( strcmp( line, "SUBSET 9.0 CMAQ\n" ) ) {
      failureMessage( "Invalid CMAQ XDR file." );
    } else {

      if ( skipInputLines( input, 2 ) ) {

        if ( readTimestamp( input, parameters->timestamp ) ) {
          Integer dimensions[ 5 ] = { 0, 0, 0, 0, 0 };

          if ( readDimensions( input, COUNT(dimensions), dimensions ) ) {
            const Integer timesteps = dimensions[ 0 ];
            const Integer variables = dimensions[ 1 ];
            const Integer layers    = dimensions[ 2 ];
            const Integer rows      = dimensions[ 3 ];
            const Integer columns   = dimensions[ 4 ];
            parameters->timesteps   = timesteps;

            if ( ! AND5( parameters->timesteps > 0,
                         IN_RANGE( variables, 4, 5 ),
                         layers > 0, rows > 0, columns > 0 ) ) {
              failureMessage( "Invalid dimensions in CMAQ XDR file." );
            } else {
              Integer subset[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };

              if ( readSubsetIndices( input, subset ) ) {
                Name cmaqVariables[ 5 ];
                Name cmaqUnits[ 5 ];
                memset( cmaqVariables, 0, sizeof cmaqVariables );
                memset( cmaqUnits, 0, sizeof cmaqUnits );
                parameters->firstLayer  = subset[ 2 ];
                parameters->lastLayer   = subset[ 3 ];
                parameters->firstRow    = subset[ 4 ];
                parameters->lastRow     = subset[ 5 ];
                parameters->firstColumn = subset[ 6 ];
                parameters->lastColumn  = subset[ 7 ];

                if ( ! AND4( parameters->timesteps ==
                               subset[ 1 ] - subset[ 0 ] + 1,
                             layers ==
                                parameters->lastLayer -
                                parameters->firstLayer + 1,
                             rows ==
                                parameters->lastRow - parameters->firstRow + 1,
                             columns ==
                                parameters->lastColumn -
                                parameters->firstColumn + 1 ) ) {
                  failureMessage( "Invalid subset indices in CMAQ XDR file." );
                } else if ( readVariablesAndUnits( input, variables,
                                                   cmaqVariables,
                                                   cmaqUnits ) ) {
                  Integer ok = 0;
                  Projector* projector =
                    parseProjectorFromXDRHeader( input, &ok );
                  memcpy( parameters->variable, cmaqVariables[ 3 ],
                          sizeof parameters->variable );
                  memcpy( parameters->units, cmaqUnits[ 3 ],
                          sizeof parameters->units );

                  if ( ok ) {
                    parameters->grid =
                      parseGridFromXDRHeader( input, projector );

                    if ( parameters->grid ) {
                      projector = 0; /* Transfered ownership to grid. */
                    } else {
                      FREE_OBJECT( projector );
                    }

                    parameters->ok =
                      AND2( parameters->grid, skipInputLines( input, 1 ) );
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  POST02( IS_BOOL( parameters->ok ),
      IMPLIES( parameters->ok,
               AND9( OR2( parameters->compareFunction,
                          parameters->convertFunction ),
                     parameters->timestamp,
                     parameters->timesteps > 0,
                     parameters->firstLayer > 0,
                     parameters->lastLayer >= parameters->firstLayer,
                     parameters->firstRow > 0,
                     parameters->lastRow >= parameters->firstRow,
                     parameters->firstColumn > 0,
                     parameters->lastColumn >= parameters->firstColumn ) ) );
}



/******************************************************************************
PURPOSE: readCMAQXDRData - Read CMAQ XDR-format file data for comparing.
INPUTS:  Stream* input           Stream to read from.
OUTPUTS: Parameters* parameters  Updated parameters->data.
******************************************************************************/

static void readCMAQXDRData( Stream* input, Parameters* parameters ) {

  PRE07( input, input->invariant( input ), input->isReadable( input ),
         input->ok( input ),
         parameters, parameters->ok,
         OR2( parameters->compareFunction, parameters->convertFunction ) );

  const Integer layers  = parameters->lastLayer - parameters->firstLayer + 1;
  const Integer rows    = parameters->lastRow - parameters->firstRow + 1;
  const Integer columns = parameters->lastColumn - parameters->firstColumn + 1;
  const Integer count   = parameters->timesteps * layers * rows * columns;
  const Integer count2  = parameters->convertFunction ? count + count : count;
  parameters->data = NEW_ZERO( Real, count2 );
  parameters->ok = 0;

  if ( parameters->data ) {
    Integer variables = 4;
    parameters->data2 = count2 > count ? parameters->data + count : 0;

    /* Read/skip over Longitudes, Latitudes, Elevations then read Data: */

    do {
      input->read32BitReals( input, parameters->data, count );
      --variables;
    } while ( AND2( input->ok( input ), variables ) );

    parameters->ok =
      AND2( input->ok( input ), isNanFree( parameters->data, count ) );

    if ( AND2( parameters->ok, parameters->data2 ) ) { /* Read 2nd data var. */
      input->read32BitReals( input, parameters->data2, count );
      parameters->ok =
        AND2( input->ok( input ), isNanFree( parameters->data2, count ) );
    }
  }

  if ( ! parameters->ok ) {
    FREE( parameters->data );
    parameters->data2 = 0;
  }

  POST02( IS_BOOL( parameters->ok ),
      IMPLIES( parameters->ok,
               AND12( OR2( parameters->compareFunction,
                           parameters->convertFunction ),
                      parameters->timestamp,
                      parameters->timesteps > 0,
                      parameters->firstLayer > 0,
                      parameters->lastLayer >= parameters->firstLayer,
                      parameters->firstRow > 0,
                      parameters->lastRow >= parameters->firstRow,
                      parameters->firstColumn > 0,
                      parameters->lastColumn >= parameters->firstColumn,
                      parameters->data,
                      isNanFree( parameters->data, count ),
                      IMPLIES( parameters->convertFunction,
                               AND3( parameters->compareFunction == 0,
                                     parameters->data2 != 0,
                                     isNanFree( parameters->data2, count))))));
}




/******************************************************************************
PURPOSE: parseProjectorFromXDRHeader - Read CMAQ XDR-format header projection.
INPUTS:  Stream* input           Stream to read from.
OUTPUTS: Integer* ok             Parsed ok?
RETURNS: Projector* Created from projection parameters in XDR header
         or 0 and failureMessage() called.
******************************************************************************/

static Projector* parseProjectorFromXDRHeader( Stream* input, Integer* ok ) {

  PRE05( input, input->invariant( input ), input->isReadable( input ),
         input->ok( input ), ok );

  Projector* result = 0;
  char line[ 1024 ] = "";
  const Integer size = COUNT( line );
  input->readString( input, line, size );
  *ok = 0;

  if ( input->ok( input ) ) {
    const char* const lambertString =
      "# lcc projection:"
      " lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis\n";
    const char* const mercatorString =
      "# mercator projection:"
      " lon_0 major_semiaxis minor_semiaxis\n";
    const char* const stereographicString =
      "# stereographic projection:"
      " lon_0 lat_0 lat_sec major_semiaxis minor_semiaxis\n";
    const char* const lonLatString =
      "# lonlat projection: major_semiaxis minor_semiaxis\n";
    Integer isLonLat = 0;
    Real majorSemiaxis = 0.0;
    Real minorSemiaxis = 0.0;
    Real centralLongitude = 0.0;
    Real centralLatitude = 0.0;

    if ( ! strcmp( line, lambertString ) ) {
      Real lowerLatitude = 0.0;
      Real upperLatitude = 0.0;
      input->readString( input, line, size );

      if ( AND12( input->ok( input ),
                  sscanf( line, "%lf %lf %lf %lf %lf %lf",
                          &lowerLatitude, &upperLatitude,
                          &centralLatitude, &centralLongitude,
                          &majorSemiaxis, &minorSemiaxis ) == 6,
                  isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                  isValidLatitude( lowerLatitude ),
                  isValidLatitude( upperLatitude ),
                  isValidLongitude( centralLongitude ),
                  isValidLatitude( centralLatitude ),
                  lowerLatitude <= upperLatitude,
                  SIGN( lowerLatitude ) == SIGN( upperLatitude ),
                  IMPLIES_ELSE( lowerLatitude >= 0.0,
                                IN_RANGE( lowerLatitude, 1.0, 89.0 ),
                                IN_RANGE( lowerLatitude, -89.0, -1.0 ) ),
                  IMPLIES_ELSE( upperLatitude >= 0.0,
                                IN_RANGE( upperLatitude, 1.0, 89.0 ),
                                IN_RANGE( upperLatitude, -89.0, -1.0 ) ),
                  IN_RANGE( centralLatitude, -89.0, 89.0 ) ) ) {

        result = (Projector*)
          newLambert( majorSemiaxis, minorSemiaxis,
                      lowerLatitude, upperLatitude,
                      centralLongitude, centralLatitude, 0.0, 0.0 );
      }
    } else if ( ! strcmp( line, mercatorString ) ) {
      input->readString( input, line, size );

      if ( AND4( input->ok( input ),
                 sscanf( line, "%lf %lf %lf",
                         &centralLongitude,
                         &majorSemiaxis, &minorSemiaxis ) == 3,
                 isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                 isValidLongitude( centralLongitude ) ) ) {

        result = (Projector*)
          newMercator( majorSemiaxis, minorSemiaxis,
                       centralLongitude, 0.0, 0.0 );
      }
    } else if ( ! strcmp( line, stereographicString ) ) {
      Real secantLatitude = 0.0;
      input->readString( input, line, size );

      if ( AND6( input->ok( input ),
                 sscanf( line, "%lf %lf %lf %lf %lf",
                         &centralLatitude, &centralLongitude,
                         &secantLatitude,
                         &majorSemiaxis, &minorSemiaxis ) == 5,
                 isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                 isValidLongitude( centralLongitude ),
                 isValidLatitude( centralLatitude ),
                 isValidLatitude( secantLatitude ) ) ) {

        result = (Projector*)
          newStereographic( majorSemiaxis, minorSemiaxis,
                            centralLongitude, centralLatitude,
                            secantLatitude, 0.0, 0.0 );
      }
    } else if ( ! strcmp( line, lonLatString ) ) {
      input->readString( input, line, size );

      if ( AND3( input->ok( input ),
                 sscanf( line, "%lf %lf", &majorSemiaxis, &minorSemiaxis) == 2,
                 isValidEllipsoid( majorSemiaxis, minorSemiaxis ) ) ) {
        isLonLat = 1;
      }
    }

    *ok = XOR2( isLonLat, result );

    if ( AND2( failureCount() == 0, ! *ok ) ) {
      failureMessage( "Invalid/unsupported projection: '%s'.", line );
    }
  }

  POST03( IMPLIES( result, input->ok( input ) ),
          IMPLIES( ! input->ok( input ), ! result ),
          IS_BOOL( *ok ) );

  return result;
}



/******************************************************************************
PURPOSE: parseGridFromXDRHeader - Read CMAQ XDR-format header proj/grid info.
INPUTS:  Stream* input           Stream to read from.
         Projector* projector    Projector for grid.
OUTPUTS: Projector* projector    0 iff result is non-NULL.
RETURNS: Grid* Grid created from projection and grid parameters in XDR header
         or 0 and failureMessage() called.
******************************************************************************/

static Grid* parseGridFromXDRHeader( Stream* input, Projector* projector ) {

  PRE06( input, input->invariant( input ), input->isReadable( input ),
         input->ok( input ), projector, projector->invariant( projector ) );

  Grid* result = 0;
  char line[ 1024 ] = "";
  const Integer size = COUNT( line );
  input->readString( input, line, size );

  if ( input->ok( input ) ) {
    const char* const gridString =
      "# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[";
    const Integer length = strlen( gridString );

    if ( ! strncmp( line, gridString, length ) ) {
      enum { MAXIMUM_LEVELS = 1000 };
      const Integer levels = atoI( line + length );
      Integer columns  = 0;
      Integer rows     = 0;
      Real westEdge    = 0.0;
      Real southEdge   = 0.0;
      Real cellWidth   = 0.0;
      Real cellHeight  = 0.0;
      Integer type     = 0;
      Real topPressure = 0.0;
      input->readString( input, line, size );

      if ( AND15( IN_RANGE( levels, 2, MAXIMUM_LEVELS ),
                  input->ok( input ),
                  sscanf( line, "%lld %lld %lf %lf %lf %lf %lld %lf",
                          &columns, &rows,
                          &westEdge, &southEdge,
                          &cellWidth, &cellHeight,
                          &type, &topPressure ) == 8,
                  columns > 0,
                  rows > 0,
                  columns * rows > 0,
                  ! isNan( westEdge ),
                  ! isNan( southEdge ),
                  ! isNan( cellWidth ),
                  ! isNan( cellHeight ),
                  cellWidth > 0.0,
                  cellHeight > 0.0,
                  IS_VALID_VERTICAL_GRID_TYPE( type ),
                  ! isNan( topPressure ),
                  topPressure > 0.0 ) ) {
        Real sigmaLevels[ MAXIMUM_LEVELS ];
        Integer level = 0;
        Integer skip = 0;
        const char* word = line;
        const Real g   = 9.81; /* HACK: should be in the CMAQ.xdr header! */
        const Real R   = 287.04;
        const Real A   = 50.0;
        const Real T0s = 290.0;
        const Real P00 = 100000.0;
        memset( sigmaLevels, 0, sizeof sigmaLevels );

        for ( skip = 0; AND2( word, skip < 8 ); ++skip ) {
          word = strchr( word + 1, ' ' );
        }

        for ( level = 0; AND2( word, level < levels ); ++level ) {
          sigmaLevels[ level ] = atoR( word );
          word = strchr( word + 1, ' ' );
        }

        result = newGrid( projector, columns, rows, westEdge, southEdge,
                          cellWidth, cellHeight,
                          levels - 1, type, topPressure, sigmaLevels,
                          g, R, A, T0s, P00 );
      }
    }

    if ( result ) {
      projector = 0; /* Grid result now owns the projector. */
    } else if ( failureCount() == 0 ) {
      failureMessage( "Invalid/unsupported grid: '%s'.", line );
    }
  }

  POST02( IMPLIES( result, AND2( input->ok( input ), projector == 0 ) ),
          IMPLIES( ! input->ok( input ), AND2( ! result, projector != 0 ) ) );

  return result;
}





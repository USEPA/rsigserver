
/******************************************************************************
PURPOSE: CMAQ.c - Define routines for processing CMAQ data.

NOTES:

HISTORY: 2007/12 plessel.todd@epa.gov, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h>  /* For stderr, fprintf(), tempnam(). */
#endif
#include <string.h> /* For memset(), memcpy().  */

#include <netcdf.h> /* For nc_close().  */

#include <Utilities.h> /* For PRE*(), NEW_ZERO(), Integer, Real, Stream. */

#include <Helpers.h>         /* For Name. */
#include <NetCDFUtilities.h> /* For createNetCDFFile(). */
#include <M3IO.h>            /* For writeM3IOHeader(). */
#include <Parameters.h>      /* For Parameters. */

/*================================== TYPES ==================================*/

typedef struct {
  UTCTimestamp timestamp;
  Name  grid;           /* E.g., "M_02_99BRACE". */
  Line  description;    /* E.g., "CMAQ CONC during Katrina". */
  Name* variable;       /* variable[ variables ]. E.g., "O3". */
  Name* units;          /* units[ variables ]. E.g., "ppmV". */
  Integer timesteps;
  Integer variables;
  Integer layers;
  Integer rows;
  Integer columns;
  Integer is64Bit;      /* Is input XDR file 64-bit data? */
  Real* data;           /* data[ layers ][ rows ][ columns ]. */
} CMAQ;

typedef Integer (*Writer)( CMAQ*, const Parameters* );

/*========================== FORWARD DECLARATIONS ===========================*/

static void deallocateCMAQ( CMAQ* cmaq );

static Integer isValidCMAQ( const CMAQ* cmaq );

static Integer readXDR( Parameters* parameters, CMAQ* cmaq );

static Projector* readProjector( Stream* input, int* isLonlat );

static Grid* readGrid( Stream* input, Projector* projector,
                       Integer columns, Integer rows, Integer layers );

static Integer parseLevels( const int vgtyp, const char* string, Integer count,
                            Real values[] );

static Integer writeASCII( CMAQ* cmaq, const Parameters* parameters );

static Integer writeCOARDS( CMAQ* cmaq, const Parameters* parameters );

static Integer writeCOARDSHeader( const CMAQ* cmaq, Integer file );

static Integer writeNetCDFData( const CMAQ* cmaq,
                                const Parameters* parameters,
                                Integer file );

static Integer writeComputedCoordinates( const Integer writeLonLats,
                                         const Integer writeElevations,
                                         const CMAQ* cmaq,
                                         const Parameters* parameters,
                                         Integer file );

static Integer writeIOAPI( CMAQ* cmaq, const Parameters* parameters );

static Integer writeIOAPIHeader( const CMAQ* cmaq,
                                 const Parameters* parameters,
                                 Integer file );

static Writer dispatcher( Integer format, Integer regrid );

static void computeLonLats( const Parameters* parameters,
                            const CMAQ* cmaq,
                            Real longitudes[], Real latitudes[] );

static void computeElevations( const Parameters* parameters,
                               const CMAQ* cmaq,
                               Real elevations[] );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: translateCMAQ - Read input and write it in another format to output.
INPUTS:  Parameters* parameters  Parameters used for translation.
OUTPUTS: Parameters* parameters  Parameters updated by translation.
******************************************************************************/

void translateCMAQ( Parameters* parameters ) {

  PRE03( isValidParameters( parameters ),
         parameters->ok,
         parameters->input->ok( parameters->input ) );

  CMAQ cmaq;
  ZERO_OBJECT( &cmaq );
  parameters->ok = 0;

  if ( readXDR( parameters, &cmaq ) ) {
    Writer writer = dispatcher( parameters->format, parameters->regrid );

    if ( ! writer ) {
      failureMessage( "Invalid/unsupported format/regrid specification." );
    } else {
      parameters->ok = writer( &cmaq, parameters );
    }
  }

  deallocateCMAQ( &cmaq );
  POST02( isValidParameters( parameters ),
          IMPLIES( AND2( parameters->ok, parameters->format == FORMAT_IOAPI ),
                   parameters->grid ) );
}



/*============================= PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: deallocateCMAQ - Deallocate contents of cmaq structure.
INPUTS: CMAQ* cmaq Structure to deallocate contents of.
******************************************************************************/

static void deallocateCMAQ( CMAQ* cmaq ) {
  PRE0( cmaq );
  FREE( cmaq->variable );
  FREE( cmaq->data );
  ZERO_OBJECT( cmaq );
}



/******************************************************************************
PURPOSE: isValidCMAQ - Check cmaq structure.
INPUTS: const CMAQ* cmaq Structure to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer isValidCMAQ( const CMAQ* cmaq ) {
  const Integer result =
    AND10( cmaq,
           cmaq->grid[ 0 ],
           cmaq->description[ 0 ],
           strlen( cmaq->timestamp ) == 24,
           cmaq->variable, cmaq->units,
           cmaq->variable[ 0 ], cmaq->units[ 0 ],
           GT_ZERO5( cmaq->timesteps, cmaq->variables, cmaq->layers,
                     cmaq->rows, cmaq->columns ),
           cmaq->data );
  return result;
}



/******************************************************************************
PURPOSE: readXDR - Read input and initialize cmaq structure.
INPUTS:  Parameters* parameters  Input parameters.
OUTPUTS: CMAQ* cmaq Structure to initialize.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer readXDR( Parameters* parameters, CMAQ* cmaq ) {

  PRE06( isValidParameters( parameters ),
         parameters->input->ok( parameters->input ),
         parameters->input->isReadable( parameters->input ),
         cmaq,
         cmaq->variable == 0,
         cmaq->data == 0 );

  Integer result = 0;
  Integer ok = 0;
  Stream* input = parameters->input;

  input->readString( input, cmaq->grid, COUNT( cmaq->grid ) );

  if ( input->ok( input ) ) {
    removeTrailingNewline( cmaq->grid );
    input->readString( input, cmaq->description, COUNT( cmaq->description ) );

    if ( input->ok( input ) ) {
      removeTrailingNewline( cmaq->description );

      if ( readTimestamp( input, cmaq->timestamp ) ) {
        Integer dimensions[ 5 ] = { 0, 0, 0, 0, 0 };

        if ( readDimensions( input, COUNT( dimensions ), dimensions ) ) {
          Integer subsetIndices[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };
          memset( subsetIndices, 0, sizeof subsetIndices );
          cmaq->timesteps = dimensions[ 0 ];
          cmaq->variables = dimensions[ 1 ];
          cmaq->layers    = dimensions[ 2 ];
          cmaq->rows      = dimensions[ 3 ];
          cmaq->columns   = dimensions[ 4 ];

          if ( readSubsetIndices( input, subsetIndices ) ) {

            if ( ! AND4( cmaq->timesteps ==
                           subsetIndices[ 1 ] - subsetIndices[ 0 ] + 1,
                         cmaq->layers ==
                           subsetIndices[ 3 ] -
                           subsetIndices[ 2 ] + 1,
                         cmaq->rows ==
                           subsetIndices[ 5 ] -
                           subsetIndices[ 4 ] + 1,
                         cmaq->columns ==
                           subsetIndices[ 7 ] -
                           subsetIndices[ 6 ] + 1 ) ) {
              failureMessage( "Invalid subset indices in CMAQ XDR file." );
            } else {
              parameters->firstLayer  = subsetIndices[ 2 ];
              parameters->lastLayer   = subsetIndices[ 3 ];
              parameters->firstRow    = subsetIndices[ 4 ];
              parameters->lastRow     = subsetIndices[ 5 ];
              parameters->firstColumn = subsetIndices[ 6 ];
              parameters->lastColumn  = subsetIndices[ 7 ];
              cmaq->variable = NEW_ZERO( Name, cmaq->variables * 2 );

              if ( cmaq->variable ) {
                cmaq->units = cmaq->variable + cmaq->variables;

                if ( readVariablesAndUnits( input, cmaq->variables,
                                            cmaq->variable, cmaq->units ) ) {


                  if ( AND2( parameters->regrid == 0,
                             parameters->grid   == 0 ) ) {
                    int isLonlat = 0;
                    Projector* projector = readProjector( input, &isLonlat );

                    if ( OR2( projector, isLonlat ) ) {
                      parameters->grid =
                        readGrid( input, projector,
                                  cmaq->columns, cmaq->rows, cmaq->layers );

                      if ( ! parameters->grid ) {
                        FREE( projector );
                      } else {
                        const char* const lines[] = {
                          "# IEEE-754 32-bit reals data[variables][timesteps][layers][rows][columns]:\n",
                          "# IEEE-754 64-bit reals data[variables][timesteps][layers][rows][columns]:\n"
                        };

                        const Integer match =
                          readMatchedLine2( input, lines[ 0], lines[ 1 ] );

                        if ( match ) {
                          cmaq->is64Bit = match == 2;
                          ok = 1;
                        }
                      }
                    }
                  }
                } else {
                  ok = skipInputLines( input, 5 );
                }
              }
            }
          }
        }
      }
    }
  }

  if ( ok ) {
    const Integer count = cmaq->layers * cmaq->rows * cmaq->columns;
    cmaq->data = NEW_ZERO( Real, count );
    result = isValidCMAQ( cmaq );
  }

  if ( AND2( ! result, failureCount() == 0 ) ) {
    failureMessage( "Invalid CMAQ data." );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   IMPLIES( parameters->format == FORMAT_IOAPI,
                          AND2( parameters->grid,
                          parameters->grid->invariant( parameters->grid )))));
  return result;
}



/******************************************************************************
PURPOSE: readProjector - Read projector from input XDR stream.
INPUTS:  Stream* input  Input XDR stream to read.
OUTPUTS: int* isLonlat  1 if lonlat projection.
RETURNS: Projection* if successful, else 0.
NOTES:   Expects to read two lines like:

# lonlat projection: major_semiaxis minor_semiaxis
6370000 6370000

# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
33 45 40 -97 6370000 6370000

# stereographic projection: lon_0 lat_0 lat_sec major_semiaxis minor_semiaxis
-98 90 45 6370000 6370000

# mercator projection: lon_0 major_semiaxis minor_semiaxis
-97 6370000 6370000

******************************************************************************/

static Projector* readProjector( Stream* input, int* isLonlat ) {

  PRE04( input, input->ok( input ), input->isReadable( input ), isLonlat );

  Projector* result = 0;
  Line line = "";
  *isLonlat = 0;

  input->readString( input, line, COUNT( line ) );

  if ( input->ok( input ) ) {

    if ( strncmp( line, "# lcc ", 6 ) == 0 ) {
      input->readString( input, line, COUNT( line ) );

      if ( input->ok( input ) ) {
        Real lowerLatitude    = 0.0;
        Real upperLatitude    = 0.0;
        Real centralLatitude  = 0.0;
        Real centralLongitude = 0.0;
        Real majorSemiaxis    = 0.0;
        Real minorSemiaxis    = 0.0;

        if ( ! AND11( sscanf( line, "%lf %lf %lf %lf %lf %lf",
                              &lowerLatitude, &upperLatitude, &centralLatitude,
                              &centralLongitude,
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
          failureMessage( "Invalid Lambert parameters '%s'.", line );
        } else {
          result = (Projector*)
            newLambert( majorSemiaxis, minorSemiaxis,
                        lowerLatitude, upperLatitude,
                        centralLongitude, centralLatitude, 0.0, 0.0 );
        }
      }
    } else if ( strncmp( line, "# stereographic ", 16 ) == 0 ) {
      input->readString( input, line, COUNT( line ) );

      if ( input->ok( input ) ) {
        Real centralLatitude  = 0.0;
        Real centralLongitude = 0.0;
        Real secantLatitude   = 0.0;
        Real majorSemiaxis    = 0.0;
        Real minorSemiaxis    = 0.0;

        if ( ! AND5( sscanf( line, "%lf %lf %lf %lf %lf",
                             &centralLongitude, &centralLatitude,
                             &secantLatitude,
                             &majorSemiaxis, &minorSemiaxis ) == 5,
                     isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                     isValidLongitude( centralLongitude ),
                     isValidLatitude( centralLatitude ),
                     isValidLatitude( secantLatitude ) ) ) {
          failureMessage( "Invalid Stereographic parameters '%s'.", line );
        } else {
          result = (Projector*)
            newStereographic( majorSemiaxis, minorSemiaxis,
                              centralLongitude, centralLatitude,
                              secantLatitude, 0.0, 0.0 );
        }
      }
    } else if ( strncmp( line, "# mercator ", 11 ) == 0 ) {
      input->readString( input, line, COUNT( line ) );

      if ( input->ok( input ) ) {
        Real centralLongitude = 0.0;
        Real majorSemiaxis    = 0.0;
        Real minorSemiaxis    = 0.0;

        if ( ! AND3( sscanf( line, "%lf %lf %lf",
                             &centralLongitude,
                             &majorSemiaxis, &minorSemiaxis ) == 3,
                     isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                     isValidLongitude( centralLongitude ) ) ) {
          failureMessage( "Invalid Mercator parameters '%s'.", line );
        } else {
          result = (Projector*)
            newMercator( majorSemiaxis, minorSemiaxis,
                         centralLongitude, 0.0, 0.0 );
        }
      }

    } else if ( strncmp( line, "# lonlat ", 9 ) == 0 ) {
      input->readString( input, line, COUNT( line ) );

      if ( ! input->ok( input ) ) {
        failureMessage( "Invalid lonlat parameters '%s'.", line );
      } else {
        *isLonlat = 1;
      }

    } else {
      failureMessage( "Invalid/unsupported projection '%s'.", line );
    }
  }

  POST0( IMPLIES( result, result->invariant( result ) ) );

  return result;
}



/******************************************************************************
PURPOSE: readGrid - Read grid description from input XDR stream.
INPUTS:  Stream* input         Input XDR stream to read.
         Projector* projector  Projector for grid or 0 for lonlat grid.
         Integer columns       Number of columns of grid cells.
         Integer rows          Number of rows    of grid cells.
         Integer layers        Number of layers  of grid cells.
RETURNS: Grid* grid if successful, else 0 and failureMessage is called.
NOTES:   Expects to read two lines like:
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[23]:
268 259 -420000 -1.716e+06 12000 12000 \
2 10000 1 0.995 0.988 0.979 0.97 0.96 0.938 0.914 0.889 0.862 0.834 \
0.804 0.774 0.743 0.694 0.644 0.592 0.502 0.408 0.311 0.21 0.106 0
******************************************************************************/

static Grid* readGrid( Stream* input, Projector* projector,
                       Integer columns, Integer rows, Integer layers ) {

  PRE07( input, input->ok( input ), input->isReadable( input ),
         IMPLIES( projector, projector->invariant( projector ) ),
         columns > 0, rows > 0, layers > 0 );

  Grid* result = 0;
  Line line = "";
  const char* const expectedLine =
    "# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[";

  input->readString( input, line, COUNT( line ) );

  if ( input->ok( input ) ) {

    if ( strncmp( line, expectedLine, strlen( expectedLine ) ) != 0 ) {
      failureMessage( "Invalid/unsupported grid '%s'.", line );
    } else {
      input->readString( input, line, COUNT( line ) );

      if ( input->ok( input ) ) {
        Integer ncols = 0;
        Integer nrows = 0;
        Real xorig = 0.0;
        Real yorig = 0.0;
        Real xcell = 0.0;
        Real ycell = 0.0;
        Integer vgtyp = 0;
        Real vgtop = 0.0;
        Real* vglvls = NEW_ZERO( Real, layers + 1 );

        if ( vglvls ) {

          if ( ! AND14( sscanf( line, "%"INTEGER_FORMAT" %"INTEGER_FORMAT
                                " %lf %lf %lf %lf %"INTEGER_FORMAT" %lf",
                                &ncols, &nrows,
                                &xorig, &yorig, &xcell, &ycell,
                                &vgtyp, &vgtop ) == 8,
                        ncols > 0,
                        nrows > 0,
                        ncols * nrows > 0,
                        ! isNan( xorig ),
                        ! isNan( yorig ),
                        ! isNan( xcell ),
                        ! isNan( ycell ),
                        xcell > 0.0,
                        ycell > 0.0,
                        IMPLIES( projector == 0,
                                 AND6( IN_RANGE( xorig, -180.0, 180.0 ),
                                       IN_RANGE( yorig, -90.0, 90.0 ),
                                       xcell <= 360.0,
                                       ycell <= 180.0,
                                       xorig + ncols * xcell <= 180.0,
                                       yorig + nrows * ycell <= 90.0 ) ),
                        IS_VALID_VERTICAL_GRID_TYPE( vgtyp ),
                        ! isNan( vgtop ),
                        vgtop > 0.0 ) ) {
            failureMessage( "Invalid grid '%s'.", line );
          } else {
            const char* const levelValues = skipWords( line, 8 );

            if ( levelValues == 0 ) {
              failureMessage( "Invalid grid '%s'.", line );
            } else if ( parseLevels( vgtyp, levelValues, layers + 1, vglvls)) {
              const double g = 9.81;       /* Gravitational force m/s^2. */
              const double R = 287.04;     /* J/kg/K = m^3/s/K. */
              const double A = 50.0;       /* Atmospheric lapse rate in K/kg*/
              const double T0s = 290.0;    /* Surface temperature in Kelvin.*/
              const double P00 = 100000.0; /* Surface pressure in Pascals. */

              result = newGrid( projector, ncols, nrows, /* columns, rows, */
                                xorig, yorig, xcell, ycell,
                                layers, vgtyp, vgtop, vglvls,
                                g, R, A, T0s, P00 );

              if ( result ) {
                projector = 0; /* Transfered ownership to Grid result. */
              }
            }
          }

          FREE( vglvls );
        }
      }
    }
  }

  POST0( IMPLIES( result, result->invariant( result ) ) );

  return result;
}



/******************************************************************************
PURPOSE: parseLevels - Parse level values from a string.
INPUTS:  const int vgtyp       Vertical grid type.
         const char* string    String to parse.
         Integer count         Number of values to parse.
OUTPUTS: Real values[ count ]  Parsed values.
RETURNS: 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer parseLevels( const int vgtyp, const char* string, Integer count,
                            Real values[] ) {

  PRE03( string, count > 0, values );

  Integer result = 0;
  Integer index  = 0;
  Real previous  = IN5( vgtyp, VGSGPH3, VGSGPN3, VGPRES3, VGWRFEM) ? 1.1 : -1.0;
  const char* word = string;

  for ( index = 0;
        AND3( word, *word, index < count );
        ++index, word = skipWords( word, 1 ) ) {
    const Real value = atoR( word );

    if ( OR2( value < 0.0,
              IMPLIES_ELSE( IN5( vgtyp, VGSGPH3, VGSGPN3, VGPRES3, VGWRFEM ),
                                 value >= previous,
                                 value <= previous ) ) ) {
      failureMessage( "Invalid level value %lf.", value );
      index = count; /* Stop looping. */
    } else {
      values[ index ] = value;
      previous = value;
    }
  }

  result = index == count;

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND3( isNanFree( values, count ),
                         minimumItem( values, count ) >= 0.0,
                         maximumItem( values, count ) <= 1.0 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: writeASCII - Write ASCII-format data.
INPUTS:  CMAQ*             cmaq        Structure to write.
         const Parameters* parameters  parameters->input.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   Reallocates cmaq->data to full size of data! Could fail if large.
******************************************************************************/

static Integer writeASCII( CMAQ* cmaq, const Parameters* parameters ) {

  PRE02( isValidCMAQ( cmaq ), isValidParameters( parameters ) );

  Integer result = 0;
  const char* const dataFormat = "\t%28.18"REAL_E_FORMAT;
  const Integer dataFormatLength = 30;
  const Integer variables    = cmaq->variables;
  const Integer timesteps    = cmaq->timesteps;
  const Integer layers       = cmaq->layers;
  const Integer rows         = cmaq->rows;
  const Integer columns      = cmaq->columns;
  const Integer layerSize    = rows * columns;
  const Integer timestepSize = layers * layerSize;
  const Integer variableSize = timesteps * timestepSize;
  const Integer dataSize     = variables * variableSize;
  const Integer outputVariables =
    variables <= 2 ? 1 + 3 + variables : 1 + variables;
  const Integer bufferSize = layerSize * outputVariables * dataFormatLength;
  FREE( cmaq->data );
  cmaq->data = NEW_ZERO( Real, dataSize ); /* Good luck... */

  if ( cmaq->data ) {
    Real* const data = cmaq->data;
    char* buffer = NEW_ZERO( char, bufferSize );

    if ( buffer ) {

      if ( cmaq->is64Bit ) {
        parameters->input->read64BitReals( parameters->input, data, dataSize );
      } else {
        parameters->input->read32BitReals( parameters->input, data, dataSize );
      }

      if ( parameters->input->ok( parameters->input ) ) {
        Stream* output = newFileStream( "-stdout", "wb" );

        if ( output ) {
          Integer variable = 0;

          /* Write header row: */

          output->writeString( output, "Timestamp(UTC)" );

          if ( variables <= 2 ) {
            output->writeString( output,
                          "\tLONGITUDE(deg)\tLATITUDE(deg)\tELEVATION(deg)" );
          }

          for ( variable = 0;
                AND2( variable < variables, output->ok( output ) );
                ++variable ) {
            output->writeString( output, "\t%s(%s)",
                                 cmaq->variable[ variable ],
                                 cmaq->units[ variable ] );
          }

          if ( output->ok( output ) ) {
            Integer timestep = 0;
            Integer yyyydddhhmm = fromUTCTimestamp( cmaq->timestamp );
            UTCTimestamp timestamp;
            const Integer timestampLength = COUNT( timestamp ) - 1;
            const Grid* const grid = parameters->grid;
            output->writeString( output, "\n" ); /* End of header line. */

            /* Write data[V][T][L][R][C] as data[T][L][R][C][V]: */

            for ( timestep = 0;
                  AND2( output->ok( output ), timestep < timesteps );
                  incrementTimestamp( &yyyydddhhmm ), ++timestep ) {

              Integer layer = 0;
              toUTCTimestamp( yyyydddhhmm, timestamp );

              /* Write one layer worth of buffered data: */

              for ( layer = 0;
                    AND2( output->ok( output ), layer < layers ); ++layer ) {
                Integer row = 0;
                char* outputBuffer = buffer;
                memset( buffer, 0, bufferSize );

                for ( row = 0; row < rows; ++row ) {
                  Integer column = 0;

                  for ( column = 0; column < columns; ++column ) {
                    CHECK( strlen( timestamp ) == timestampLength );
                    strncpy( outputBuffer, timestamp, timestampLength );
                    outputBuffer += timestampLength;

                    if ( variables <= 2 ) { /* Compute and write coordinates:*/
                      const double elevation = layer;
                      const double latitude =
                        grid->southEdge( grid ) +
                        ( parameters->firstRow + row - 1 ) *
                        grid->cellHeight( grid );
                      const double longitude =
                        grid->westEdge( grid ) +
                        ( parameters->firstColumn + column - 1 ) *
                        grid->cellWidth( grid );
                      char word[ 40 ] = "";
                      Integer count = 0;
                      count = sprintf( word, dataFormat, longitude );
                      CHECK( count < sizeof word / sizeof *word );
                      CHECK( strlen( word ) == count );
                      strncpy( outputBuffer, word, count );
                      outputBuffer += count;
                      count = sprintf( word, dataFormat, latitude );
                      CHECK( count < sizeof word / sizeof *word );
                      CHECK( strlen( word ) == count );
                      strncpy( outputBuffer, word, count );
                      outputBuffer += count;
                      count = sprintf( word, dataFormat, elevation );
                      CHECK( count < sizeof word / sizeof *word );
                      CHECK( strlen( word ) == count );
                      strncpy( outputBuffer, word, count );
                      outputBuffer += count;
                    }

                    for ( variable = 0; variable < variables; ++variable ) {
                      Integer dataIndex =
                        variable * variableSize + timestep * timestepSize +
                        layer * layerSize + row * columns + column;
                      char word[ 40 ] = "";
                      Integer count = 0;
                      CHECK( IN_RANGE( dataIndex, 0, dataSize - 1 ) );
                      count = sprintf( word, dataFormat, data[ dataIndex ] );
                      CHECK( count < sizeof word / sizeof *word );
                      CHECK( strlen( word ) == count );
                      strncpy( outputBuffer, word, count );
                      outputBuffer += count;
                    }

                    strcpy( outputBuffer, "\n" ); /* End of spreadsheet row. */
                    ++outputBuffer;
                  }
                }

                /* Write buffered output to stream: */

                CHECK( strlen( buffer ) < bufferSize );
                output->writeString( output, buffer );
              }
            }

            result = output->ok( output );
          }

          FREE_OBJECT( output );
        }
      }

      FREE( buffer );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDS - Write COARDS-format data.
INPUTS:  CMAQ* cmaq                    Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDS( CMAQ* cmaq, const Parameters* parameters ) {

  PRE02( isValidCMAQ( cmaq ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    cmaq->timesteps * cmaq->variables *
    cmaq->layers * cmaq->rows * cmaq->columns * 4 +
    cmaq->timesteps * 4 +
    5000;
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {

    if ( writeCOARDSHeader( cmaq, file ) ) {
      result = writeNetCDFData( cmaq, parameters, file );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSHeader - Write cmaq header info to outputFileName.
INPUTS:  const CMAQ* cmaq  Structure to write.
OUTPUTS: Integer file      NetCDF fle to write to.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeCOARDSHeader( const CMAQ* cmaq, Integer file ) {

  PRE02( isValidCMAQ( cmaq ), file != -1 );

  Integer result = 0;
  enum { TIME, Z, Y, X, DIMENSIONS };
  const char* const names[ DIMENSIONS ] = {
    "time", "elevation", "latitude", "longitude"
  };
  Integer dimensionIds[ DIMENSIONS ] = { -1, -1, -1, -1 };
  Integer dimensions[   DIMENSIONS ] = { 0, 0, 0, 0 };
  dimensions[ TIME ] = cmaq->timesteps;
  dimensions[ Z    ] = cmaq->layers;
  dimensions[ Y    ] = cmaq->rows;
  dimensions[ X    ] = cmaq->columns;

  if ( createDimensions( file, DIMENSIONS, names, dimensions, dimensionIds )) {

    if ( createCRSVariable( file ) != -1 ) {

      if ( createLongitudeAndLatitude( file, DIMENSIONS, dimensionIds ) ) {
        const Integer height =
          createVariable( file, "height", "meters", NC_FLOAT, 1, DIMENSIONS,
                          dimensionIds );

        if ( height != -1 ) {

          if ( writeTextAttribute( file, height, "positive", "up" ) ) {
            const Integer variables = cmaq->variables;
            const Integer noLonlats = cmaq->variables < 3;
            const Integer noElevations =
              OR2( AND3( noLonlats == 1,
                         strcmp( cmaq->variable[0], "ELEVATION" ),
                         strcmp( cmaq->variable[0], "elevation" ) ),
                  AND3( noLonlats == 0,
                        strcmp(cmaq->variable[2], "ELEVATION" ),
                        strcmp(cmaq->variable[2], "elevation" ) ) );
            Integer index = 2 * ( noLonlats == 0 ) + ( noElevations == 0 );

            do {
              const Integer ok =
                createVariable( file, cmaq->variable[ index ],
                                cmaq->units[ index ], NC_FLOAT, 1,
                                DIMENSIONS, dimensionIds ) != -1;

              if ( ! ok ) {
                index = variables;
              }

              ++index;
            } while ( index < variables );

            if ( index == variables ) {

              if ( writeTextAttribute( file, NC_GLOBAL, "grid", cmaq->grid )) {

                if ( writeTextAttribute( file, NC_GLOBAL, "description",
                                         cmaq->description ) ) {
                  const char* const history =
                    "http://www.ncep.noaa.gov,EPA-RTP,"
                    ",CMAQSubset,XDRConvert";

                  result = writeStandardContents( file, history,
                                                  cmaq->timestamp,
                                                  dimensionIds[ TIME ],
                                                  cmaq->timesteps, 1 );
                }
              }
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeNetCDFData - Write COARDS or IOAPI format cmaq data.
INPUTS:  Stream* input                Stream to read data arrays from.
         const CMAQ* cmaq             Structure to write.
         const Parameters* parameters Structure to query.
OUTPUTS: Integer file                 NetCDF file to write to.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeNetCDFData( const CMAQ* cmaq,
                                const Parameters* parameters,
                                Integer file ) {

  PRE07( isValidCMAQ( cmaq ), isValidParameters( parameters ),
        IN3( parameters->format, FORMAT_COARDS, FORMAT_IOAPI ),
        parameters->input, parameters->input->ok( parameters->input ),
        parameters->input->isReadable( parameters->input ),
        file != -1 );

  const Integer isCOARDS = parameters->format == FORMAT_COARDS;
  const Integer noLonlats = cmaq->variables < 3;
  const Integer noElevations =
    OR2( AND3( noLonlats == 1,
               strcmp(cmaq->variable[0], "ELEVATION" ),
               strcmp(cmaq->variable[0], "elevation" ) ),
        AND3( noLonlats == 0,
              strcmp(cmaq->variable[2], "ELEVATION" ),
              strcmp(cmaq->variable[2], "elevation" ) ) );
  const Integer variables = cmaq->variables;
  const Integer timesteps = cmaq->timesteps;
  const Integer layers    = cmaq->layers;
  const Integer rows      = cmaq->rows;
  const Integer columns   = cmaq->columns;
  const Integer count = layers * rows * columns;
  Integer result = 1;

  if ( OR2( noLonlats, noElevations ) ) {
    result =
      writeComputedCoordinates(noLonlats, noElevations, cmaq, parameters, file);
  }

  if ( result ) {
    Integer variable = 0;

    do {
      const char* const variableName = cmaq->variable[ variable ];
      const char* const outputVariableName =
        isCOARDS ?
          ( ! strcmp( variableName, "LATITUDE" ) ? "latitude"
          : ! strcmp( variableName, "LONGITUDE" ) ? "longitude"
          : ! strcmp( variableName, "ELEVATION" ) ? "height"
          : ! strcmp( variableName, "elevation" ) ? "height"
          : variableName )
        :
          ( ! strcmp( variableName, "latitude" ) ? "LATITUDE"
          : ! strcmp( variableName, "longitude" ) ? "LONGITUDE"
          : ! strcmp( variableName, "elevation" ) ? "ELEVATION"
          : variableName );

      Stream* const input = parameters->input;
      Integer timestep = 0;

      do {

        if ( cmaq->is64Bit ) {
          input->read64BitReals( input, cmaq->data, count );
        } else {
          input->read32BitReals( input, cmaq->data, count );
        }

        result = input->ok( input );

        if ( result ) {
          result = writeSomeData( file, outputVariableName, timestep,
                                  1, layers, rows, columns, cmaq->data );
        }

        if ( ! result ) {
          timestep = timesteps;
          variable = variables;
        }

        ++timestep;
      } while ( timestep < timesteps );

      ++variable;
    } while ( variable < variables );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeComputedCoordinates - Compute and write longitude, latitude,
         elevation to a NetCDF file.
INPUTS:  const Integer writeLonLats     Write longitudes, latitudes?
         const Integer writeElevations  Write elevations?
         const CMAQ* cmaq               Structure to query.
         const Parameters* parameters   Structure to query.
OUTPUTS: Integer file     NetCDF file to write to.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeComputedCoordinates( const Integer writeLonLats,
                                         const Integer writeElevations,
                                         const CMAQ* cmaq,
                                         const Parameters* parameters,
                                         Integer file ) {

  PRE05( IS_BOOL( writeLonLats ), IS_BOOL( writeElevations ),
         isValidCMAQ( cmaq ), isValidParameters( parameters ), file != -1 );

  const Integer timesteps = cmaq->timesteps;
  const Integer layers    = cmaq->layers;
  const Integer rows      = cmaq->rows;
  const Integer columns   = cmaq->columns;
  const Integer count = layers * rows * columns;
  Integer timestep = 0;
  Integer result = 1;

  if ( writeLonLats ) {
    Real* lonlats = NEW_ZERO( Real, 2 * count );
    result = lonlats != 0;

    if ( result ) {
      Real* const longitudes = lonlats;
      Real* const latitudes  = longitudes + count;
      computeLonLats( parameters, cmaq, longitudes, latitudes );
      timestep = 0;

      do {

        /* Must use cmaq->data as a temporary since writeSomeData() ruins it.*/

        memcpy( cmaq->data, longitudes, count * sizeof *longitudes );
        result = writeSomeData( file,
                                parameters->format == FORMAT_COARDS ?
                                  "longitude" : "LONGITUDE",
                                timestep,
                                1, layers, rows, columns, cmaq->data );

        if ( result ) {
          memcpy( cmaq->data, latitudes, count * sizeof *latitudes );
          result = writeSomeData( file,
                                  parameters->format == FORMAT_COARDS ?
                                    "latitude" : "LATITUDE",
                                  timestep,
                                  1, layers, rows, columns, cmaq->data );
        }

        ++timestep;
      } while ( result && timestep < timesteps );

      FREE( lonlats );
    }
  }

  if ( result ) {

    if ( writeElevations ) {
      Real* elevations = NEW_ZERO( Real, count );
      result = elevations != 0;

      if ( result ) {
        computeElevations( parameters, cmaq, elevations );
        timestep = 0;

        do {
          memcpy( cmaq->data, elevations, count * sizeof *elevations );
          result = writeSomeData( file,
                                  parameters->format == FORMAT_COARDS ?
                                    "height" : "ELEVATION",
                                  timestep,
                                  1, layers, rows, columns, cmaq->data );
          ++timestep;
        } while ( result && timestep < timesteps );

        FREE( elevations );
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeIOAPI - Write IOAPI-format data.
INPUTS:  CMAQ* cmaq                    Structure to write.
         const Parameters* parameters  Parameters.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeIOAPI( CMAQ* cmaq, const Parameters* parameters ) {

  PRE02( isValidCMAQ( cmaq ), isValidParameters( parameters ) );

  Integer result = 0;
  const Integer fileSizeEstimate =
    cmaq->timesteps * cmaq->variables *
    cmaq->layers * cmaq->rows * cmaq->columns * 4 + 10000;
  const Integer create64BitFile = fileSizeEstimate > TWO_GB;
  Integer file =
    createNetCDFFile( parameters->netcdfFileName, create64BitFile );

  if ( file != -1 ) {

    if ( writeIOAPIHeader( cmaq, parameters, file ) ) {
      result = writeNetCDFData( cmaq, parameters, file );
    }

    nc_close( file );
    file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeIOAPIHeader - Write header to file.
INPUTS:  const CMAQ* cmaq  Structure to write.
         const Parameters* parameters  Parameters structure.
OUTPUTS: Integer file      NetCDF file to write to.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

static Integer writeIOAPIHeader( const CMAQ* cmaq,
                                 const Parameters* parameters,
                                 Integer file ) {

  PRE03( isValidCMAQ( cmaq ), isValidParameters( parameters ), file != -1 );

  Integer result = 0;
  const Integer noLonlats = cmaq->variables < 3;
  const Integer noElevations =
    OR2( AND3( noLonlats == 1,
               strcmp(cmaq->variable[0], "ELEVATION" ),
               strcmp(cmaq->variable[0], "elevation" ) ),
        AND3( noLonlats == 0,
              strcmp(cmaq->variable[2], "ELEVATION" ),
              strcmp(cmaq->variable[2], "elevation" ) ) );
  Name variableNames[ 120 ];
  Name variableUnits[ 120 ];
  Integer variableIndex = 0;
  const Integer variables = cmaq->variables;
  Integer variable = 0;
  const Integer firstTimestamp = fromUTCTimestamp( cmaq->timestamp );
  Line history = "";

  appendToLine( history, cmaq->description );
  appendToLine( history, ",XDRConvert" );
  memset( variableNames, 0, sizeof variableNames );
  memset( variableUnits, 0, sizeof variableNames );

  if ( noLonlats ) {
    memcpy( variableNames[ 0 ], "LONGITUDE", sizeof (Name) );
    memcpy( variableNames[ 1 ], "LATITUDE", sizeof (Name) );
    variableIndex = 2;
  }

  if ( noElevations ) {
    memcpy( variableNames[ 2 ], "ELEVATION", sizeof (Name) );
    variableIndex = 3;
  }

  for ( variable = 0; variable < variables; ++variable, ++variableIndex ) {
    memcpy( variableNames[ variableIndex ], cmaq->variable[ variable ],
            sizeof (Name) );
    memcpy( variableUnits[ variableIndex ], cmaq->units[ variable ],
            sizeof (Name) );
  }

  {
    Grid* subsetGrid =
      newSubsetGrid( parameters->grid,
                     parameters->firstLayer - 1, parameters->lastLayer - 1,
                     parameters->firstRow - 1, parameters->lastRow - 1,
                     parameters->firstColumn - 1, parameters->lastColumn - 1 );

    result = subsetGrid != 0;

    if ( result ) {
      result = writeM3IOHeader( file, cmaq->timesteps, 1, firstTimestamp,
                                variableIndex, cmaq->layers,
                                variableNames, variableUnits,
                                history, subsetGrid );
      FREE_OBJECT( subsetGrid );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: dispatcher - Look-up and return a writer for the given format/regrid.
INPUTS:  Integer format  E.g., FORMAT_XDR.
         Integer regrid  E.g., 0 or AGGREGATE_MEAN.
RETURNS: Writer a writer routine or 0 if none found.
******************************************************************************/

static Writer dispatcher( Integer format, Integer regrid ) {

  PRE02( IS_VALID_FORMAT( format ),
         OR2( regrid == 0, IS_VALID_AGGREGATE_METHOD( regrid ) ) );

  Writer result = 0;

  if ( ! regrid ) {
    typedef struct {
      Integer format; /* FORMAT_XDR, etc. */
      Writer writer;  /* Routine that writes data in this format. */
    } Entry;

    static const Entry writers[] = {
      { FORMAT_XDR,    0,          },
      { FORMAT_ASCII,  writeASCII  },
      { FORMAT_COARDS, writeCOARDS },
      { FORMAT_IOAPI,  writeIOAPI  },
      { -1, 0 }
    };

    const Integer count = COUNT( writers );
    Integer index = 0;

    do {
      const Entry* const entry = writers + index;

      if ( entry->format == -1 ) {
        index = count;
      } else if ( entry->format == format ) {
        result = entry->writer;
        index = count;
      }

      ++index;
    } while ( index < count );
  }

  return result;
}



/******************************************************************************
PURPOSE: computeLonLats - Compute longitude, latitude coordinates.
INPUTS:  const Parameters* parameters  Structure to query.
         const CMAQ* cmaq              Structure to query.
OUTPUTS: Real longitudes[ cmaq->layers * cmaq->rows * cmaq->columns ]
         Real latitudes[  cmaq->layers * cmaq->rows * cmaq->columns ]
******************************************************************************/

static void computeLonLats( const Parameters* parameters,
                            const CMAQ* cmaq,
                            Real longitudes[], Real latitudes[] ) {

  PRE04( isValidParameters( parameters ), isValidCMAQ( cmaq ),
         longitudes, latitudes );

  const Integer layers  = cmaq->layers;
  const Integer rows    = cmaq->rows;
  const Integer columns = cmaq->columns;
  const Integer layerSize = rows * columns;
  const Grid* const grid = parameters->grid;
  const Integer lastRow = parameters->lastRow;
  const Integer firstColumn = parameters->firstColumn;
  const Integer lastColumn = parameters->lastColumn;
  Integer row = parameters->firstRow - 1;
  Integer layer = 0;
  Integer index = 0;

  for ( ; row < lastRow; ++row ) {
    Integer column = firstColumn - 1;

    for ( ; column < lastColumn; ++column, ++index ) {
      const Real longitude = grid->longitude( grid, row, column );
      const Real latitude = grid->latitude( grid, row, column );
      longitudes[ index ] = longitude;
      latitudes[  index ] = latitude;
    }
  }

  /* Copy to all layers: */

  for ( layer = 0, index = 0; layer < layers; ++layer, index += layerSize ) {
    memcpy( longitudes + index, longitudes, layerSize * sizeof *longitudes );
    memcpy( latitudes  + index, latitudes,  layerSize * sizeof *latitudes );
  }
}



/******************************************************************************
PURPOSE: computeElevations - Compute elevation coordinates.
INPUTS:  const Parameters* parameters  Structure to query.
         const CMAQ* cmaq              Structure to query.
OUTPUTS: Real elevations[ cmaq->layers * cmaq->rows * cmaq->columns ]
******************************************************************************/

static void computeElevations( const Parameters* parameters,
                               const CMAQ* cmaq,
                               Real elevations[] ) {

  PRE03( isValidParameters( parameters ), isValidCMAQ( cmaq ), elevations );

  const Integer layers  = cmaq->layers;
  const Integer rows    = cmaq->rows;
  const Integer columns = cmaq->columns;
  const Integer layerSize = rows * columns;
  Integer layer = 0;
  const Grid* const grid = parameters->grid;
  Real* output = elevations;

  for ( layer = 0; layer < layers; ++layer ) {
    Integer counter = layerSize;
    const Real elevation = grid->elevation( grid, layer );

    while ( counter-- ) {
      *output++ = elevation;
    }
  }

  POST02( IN_RANGE( elevations[ 0 ], -500.0, 100000.0 ),
          IN_RANGE( elevations[ cmaq->layers * cmaq->rows * cmaq->columns - 1 ],
                    elevations[ 0 ], 100000.0 ) );
}




/******************************************************************************
PURPOSE: Helpers.c - Define general utility helper routines.

NOTES:   See Helpers.h.

HISTORY: 2007/12 plessel.todd@epa.gov, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <ctype.h>  /* For toupper(), tolower(). */
#include <string.h> /* For memset().  */
#include <float.h>  /* For FLT_MAX.  */
#include <stdlib.h> /* For qsort().  */

#include <Helpers.h> /* For public interface. */

/*================================= GLOBALS =================================*/

/* Compare operators: */

static Real difference( Real a, Real b ) {
  PRE02( ! isNan( a ), ! isNan( b ) );
  const Real safeResult = safeDifference( a, b );
  const Real result = CLAMPED_TO_RANGE( safeResult, -FLT_MAX, FLT_MAX );
  POST0( ! isNan( result ) );
  return result;
}

static Real absoluteDifference( Real a, Real b ) {
  PRE02( ! isNan( a ), ! isNan( b ) );
  const Real safeResult = safeDifference( a, b );
  const Real absResult = safeResult < 0.0 ? -safeResult : safeResult;
  const Real result = absResult > FLT_MAX ? FLT_MAX : absResult;
  POST02( ! isNan( result ), result >= 0.0 );
  return result;
}

static Real percentDifference( Real a, Real b ) {
  PRE02( ! isNan( a ), ! isNan( b ) );
  const Real diffResult = safeDifference( a, b );
  const Real numerator = diffResult < 0.0 ? -diffResult : diffResult;
  const Real denominator = safeSum( a, b );
  const Real safeResult =
    numerator == 0.0 ? 0.0
    : denominator == 0.0 ? FLT_MAX
    : 200.0 * safeQuotient( numerator, denominator );
  const Real result = CLAMPED_TO_RANGE( safeResult, -FLT_MAX, FLT_MAX );
  POST0( ! isNan( result ) );
  return result;
}

static Real ratio( Real a, Real b ) {
  PRE02( ! isNan( a ), ! isNan( b ) );
  const Real safeResult =
    a == 0.0 ? 0.0
    : b == 0.0 ? ( a < 0.0 ? -FLT_MAX : FLT_MAX )
    : safeQuotient( a, b );
  const Real result = CLAMPED_TO_RANGE( safeResult, -FLT_MAX, FLT_MAX );
  POST0( ! isNan( result ) );
  return result;
}

static Real replace( Real a, Real b ) {
  PRE02( ! isNan( a ), ! isNan( b ) );
  const Real result = b;
  POST0( result == b );
  return result;
}

static Real scaledOffset( Real input, Real scale, Real offset ) {
  PRE03( ! isNan( input ), ! isNan( scale ), ! isNan( offset ) );
  const Real missing = -9999.0;
  const Real safeResult =
    AND3( input > missing, scale > missing, offset > missing ) ?
  safeSum( safeProduct( input, scale ), offset ) : missing;
  const Real result = CLAMPED_TO_RANGE( safeResult, -FLT_MAX, FLT_MAX );
  POST0( ! isNan( result ) );
  return result;
}

typedef struct {
  const char* name;
  const char* variable;      /* Optional: altered variable name. */
  const char* units;         /* Optional: altered variable units. */
  CompareFunction comparer;
  ConvertFunction converter;
} CompareEntry;

static const CompareEntry compareEntries[] = {
  { "difference", 0, 0, difference, 0 },
  { "absolute_difference", 0, 0, absoluteDifference, 0 },
  { "percent_difference", 0, "%", percentDifference, 0 },
  { "ratio", 0, "-", ratio, 0 },
  { "replace", 0, 0, replace, 0 },
  { "convert", "PM25", "ug/m3", 0, scaledOffset },
};



/* Grid cell info for aggregating: */

typedef struct {
  Integer index;
  Integer count;
  Real value;
  Real value2;
} Cell;



/* (Re)initialize cells: */

static void initializeCells( const Integer count, Cell cells[] ) {
  Integer index = 0;

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    Cell* const cell = cells + index;
    cell->index = -1;
    cell->count = 0;
    cell->value = 0.0;
    cell->value2 = 0.0;
  }
}



/* Called by qsort(): */

static int compareCells( const void* va, const void* vb ) {
  const Cell* const a = va;
  const Cell* const b = vb;
  const int result = a->index - b->index;
  return result;
}



/* Copy mean values from cells to output arrays: */

static Integer copyAggregatedData( const Integer cellCount,
                                   Cell cells[],
                                   const Integer offset,
                                   Real longitudes[],
                                   Real latitudes[],
                                   Real elevations[],
                                   Integer columns[],
                                   Integer rows[],
                                   Integer layers[],
                                   Real data[],
                                   Real data2[],
                                   RegriddedNote notes[] ) {

  PRE08( cellCount > 0, cells, longitudes, latitudes, columns, rows, data,
         IMPLIES_ELSE( elevations, layers, ! layers ) );

  Integer result = 0;
  Integer output = offset;
  Integer cellIndex = 0; /* Index of first cell with aggregated data: */

  /* Sort cells by index to output in original input order: */

  qsort( cells, cellCount, sizeof *cells, compareCells );

  while ( cellIndex < cellCount && cells[ cellIndex ].index == -1 ) {
    ++cellIndex;
  }

  while ( cellIndex < cellCount ) {
    Cell* const cell = cells + cellIndex;
    const Integer input = cell->index;
    CHECK( input >= output );
    cell->value /= cell->count;
    data[ output ] = cell->value;
    longitudes[ output ] = longitudes[ input ];
    latitudes[  output ] = latitudes[  input ];
    columns[    output ] = columns[    input ];
    rows[       output ] = rows[       input ];

    if ( layers ) {
      layers[ output ] = layers[ input ];
      elevations[ output ] = elevations[ input ];
    }

    if ( data2 ) {
      cell->value2 /= cell->count;
      data2[ output ] = cell->value2;
    }

    if ( notes ) {
      memcpy( notes + output, notes + input, sizeof *notes );
    }

    ++result;
    ++output;
    ++cellIndex;
  }

  POST0( IN_RANGE( result, 1, cellCount ) );
  return result;
}



/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: compareFunction - Return comparer function given its name.
INPUTS:  const char* name  Name of function, e.g., "ratio".
RETURNS: CompareFunction or 0 if not found.
******************************************************************************/

CompareFunction compareFunction( const char* name ) {
  PRE0( name );
  CompareFunction result = 0;
  const Integer count = COUNT( compareEntries );
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( ! strcmp( name, compareEntries[ index ].name ) ) {
      result = compareEntries[ index ].comparer;
      index = count;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: convertFunction - Return converter function given its name.
INPUTS:  const char* name  Name of function, e.g., "convert".
RETURNS: ConvertFunction or 0 if not found.
******************************************************************************/

ConvertFunction convertFunction( const char* name ) {
  PRE0( name );
  ConvertFunction result = 0;
  const Integer count = COUNT( compareEntries );
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( ! strcmp( name, compareEntries[ index ].name ) ) {
      result = compareEntries[ index ].converter;
      index = count;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: compareFunctionNameUnits - Update units if needed for compare.
INPUTS:  CompareFunction comparer  E.g., ratio.
         ConvertFunction converter E.g., scaledOffset.
         Name name                 E.g., "ozone".
         Name units                E.g., "ppm".
         Name otherName            E.g., "O3".
         Name otherUnits           E.g., "ppmV".
OUTPUTS: Name name                 E.g., "ozone_ratio".
         Name units                E.g., "-".
******************************************************************************/

void compareFunctionNameUnits( CompareFunction comparer,
                               ConvertFunction converter,
                               Name name, Name units,
                               Name otherName, Name otherUnits ) {
  PRE09( OR2( comparer, converter ),
         name, *name, units, *units, otherName, *otherName,
         otherUnits, *otherUnits );
  const Integer count = COUNT( compareEntries );
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {
    const CompareEntry* const entry = compareEntries + index;

    if ( OR2( AND2( comparer,  comparer  == entry->comparer ),
              AND2( converter, converter == entry->converter ) ) ) {
      const Integer length = sizeof (Name) / sizeof *name - strlen( name ) - 2;

      if ( entry->variable ) {
        CHECK2( entry->variable[ 0 ],
                strlen( entry->variable ) < sizeof (Name) / sizeof *name );
        strncpy( name, entry->variable, sizeof (Name) / sizeof *name );
      } else if ( length > 0 ) {
        CHECK3( entry->name, entry->name[ 0 ],
                strlen( name ) + length < sizeof (Name) / sizeof *name );
        strcat( name, "_" );
        strncat( name, entry->name, length );
      }

      if ( entry->units ) {
        CHECK2( entry->units[ 0 ],
                strlen( entry->units ) < sizeof (Name) / sizeof *units );
        strncpy( units, entry->units, sizeof (Name) / sizeof *name );
      } else if ( comparer == replace ) {
        strcpy( name, otherName );
        strcpy( units, otherUnits );
      }

      index = count;
    }
  }
}



/******************************************************************************
PURPOSE: sum - Return sum of data items.
INPUTS:  Integer count  Number of data items.
         const Integer data[ count ]  Data items to sum.
RETURNS: Integer sum of data items.
******************************************************************************/

Integer sum( Integer count, const Integer data[] ) {
  PRE02( count > 0, data );
  Integer result = 0;
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {
    result += data[ index ];
  }

  return result;
}



/******************************************************************************
PURPOSE: scale - Scale data items.
INPUTS:  Real factor    Factor to scale data items by.
         Integer count  Number of data items.
OUTPUTS: Real data[ count ]  Data items multiplied by factor.
******************************************************************************/

void scale( Real factor, Integer count, Real data[] ) {
  PRE04( ! isNan( factor ), count > 0, data, isNanFree( data, count ) );
  Integer index = 0;

  for ( index = 0; index < count; ++index ) {
    data[ index ] *= factor;
  }

  POST0( isNanFree( data, count ) );
}



/******************************************************************************
PURPOSE: streamFile - Read/stream a file to stdout.
INPUTS:  const char* fileName  File to read/stream.
RETURNS: 1 if successful, else 0 and failureMessage() called.
******************************************************************************/

Integer streamFile( const char* fileName ) {

  PRE02( fileName, *fileName );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer bufferSize = 10 * 1024 * 1024; /* 10MB. */
    void* buffer = NEW_ZERO( char, bufferSize );

    if ( buffer ) {
      Stream* input = newFileStream( fileName, "rb" );

      if ( input ) {

        if ( output->ok( output ) ) {

          do {
            Integer bytesRead = 0;
            DEBUG( fprintf( stderr, "reading 10MB...\n" ); )
            input->readUpToNBytes( input, buffer, bufferSize, &bytesRead );
            result = input->ok( input );

            if ( AND2( result, bytesRead > 0 ) ) {
              DEBUG( fprintf( stderr, "writing 10MB...\n" ); )
              output->writeBytes( output, buffer, bytesRead );
              result = output->ok( output );
            }

          } while ( AND2( result, ! input->isAtEnd( input ) ) );
        }

        FREE_OBJECT( input );
      }

      FREE( buffer );
    }

    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: copyToStdout - Read from input and write to stdout.
INPUTS:  Stream* input  The stream to read from.
RETURNS: Integer 1 if successful, else 0.
******************************************************************************/

Integer copyToStdout( Stream* input ) {

  PRE02( input, input->isReadable( input ) );

  Integer result = 0;
  Stream* output = newFileStream( "-stdout", "wb" );

  if ( output ) {
    const Integer bufferSize = 1024 * 1024;
    void* buffer = NEW_ZERO( char, bufferSize );

    if ( buffer ) {

      do {
        Integer bytesRead = 0;
        input->readUpToNBytes( input, buffer, bufferSize, &bytesRead );
        result = input->ok( input );

        if ( AND2( result, bytesRead > 0 ) ) {
          output->writeBytes( output, buffer, bytesRead );
          result = output->ok( output );
        }

      } while ( AND2( result, ! input->isAtEnd( input ) ) );

      FREE( buffer );
    }

    FREE_OBJECT( output );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: skipInputLines - Skip lines of input.
INPUTS:  Stream* input  Stream to read from.
         Integer count  Number of lines to read/skip.
OUTPUTS: Stream* input  Stream read from.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer skipInputLines( Stream* input, Integer count ) {
  PRE03( input, input->isReadable( input ), count > 0 );
  char line[ 1024 ] = "";
  Integer index = 0;
  Integer result = 0;

  do {
    input->readString( input, line, COUNT( line ) );

    if ( ! input->ok( input ) ) {
      index = count;
    }

    ++index;
  } while ( index < count );

  result = index == count;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: aggregateName - Generate aggregate name of variable.
INPUTS:  const char* const inputVariableName  Name of variable.
         const Integer hoursPerTimestep       Number of hours per aggregation.
OUTPUTS: Name outputVariableName              Output name. E.g., daily_no2.
******************************************************************************/

void aggregateName( const char* const inputVariableName,
                    const Integer hoursPerTimestep,
                    Name outputVariableName ) {
  PRE04( inputVariableName, *inputVariableName, hoursPerTimestep > 0,
         outputVariableName );
  const char* const prefix0 =
    hoursPerTimestep == 1 ? "" /* Default is hourly. No prefix. */
    : hoursPerTimestep == 24 ? "daily_"
    : hoursPerTimestep <= 31 * 24 ? "monthly_"
    : hoursPerTimestep <= 3 * 31 * 24 ? "seasonal_"
    : IN_RANGE( hoursPerTimestep, 365 * 24, 366 * 24 ) ? "yearly_"
    : "";
  const char* const prefix =
    ! strstr( inputVariableName, prefix0 ) ? prefix0 : ""; /* Avoid duplicate*/
  const size_t length = strlen( prefix );
  memset( outputVariableName, 0, sizeof (Name) );

  if ( *prefix ) {
    strncpy( outputVariableName, prefix, sizeof (Name) / sizeof (char) - 1 );
  }

  strncpy( outputVariableName + length, inputVariableName,
            sizeof (Name) / sizeof (char) - length - 1 );
  CHECK( outputVariableName[ sizeof (Name) / sizeof (char) - 1 ] == '\0' );

  /* HACKs to shorten name: */

  if ( strlen( outputVariableName ) > 15 ) {
    char* replace = 0;

    if ( ( replace = strstr( outputVariableName, "nitrogendioxide" ) ) ) {
      ++replace; /* Skip 'n'. */
      *replace++ = 'o';
      *replace++ = '2';
      *replace = '\0'; /* no2. */
    } else if ( ( replace = strstr( outputVariableName, "carbonmonoxide" ) ) ) {
      ++replace; /* Skip 'c'. */
      *replace++ = 'o';
      *replace = '\0'; /* co. */
    } else if ( ( replace = strstr( outputVariableName, "formaldehyde" ) ) ) {
      *replace++ = 'h';
      *replace++ = 'c';
      *replace++ = 'h';
      *replace++ = 'o';
      *replace = '\0'; /* hcho. */
    }

/*    outputVariableName[ 15 ] = '\0'; */
  }

  POST03( outputVariableName[ 0 ],
          outputVariableName[ sizeof (Name) / sizeof (char) - 1 ] == '\0',
          strlen( outputVariableName ) <= 80 /* 15 */ );
}



/******************************************************************************
PURPOSE: readMatchedLine - Read a line of input and check if it matches pattern.
INPUTS:  Stream* input              Stream to read from.
         const char* const pattern  Line to match.
OUTPUTS: Stream* input              Stream read from.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readMatchedLine( Stream* input, const char* const pattern ) {
  PRE04( input, input->isReadable( input ), pattern, *pattern );
  char line[ 1024 ] = "";
  Integer result = 0;
  input->readString( input, line, COUNT( line ) );
  result = input->ok( input );

  if ( result ) {
    result = ! strcmp( line, pattern );

    if ( ! result ) {
      failureMessage( "Invalid line in file.\n'%s'\nexpected '%s'\n",
                      line, pattern);
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readMatchedLine2 - Read a line of input and check if it matches
         pattern1 or pattern2.
INPUTS:  Stream* input               Stream to read from.
         const char* const pattern1  1st line to match.
         const char* const pattern2  Else 2nd line to match.
OUTPUTS: Stream* input              Stream read from.
RETURNS: Integer 1 if matches pattern1, 2 if matches pattern2, else 0 and
         failureMessage() was called.
******************************************************************************/

Integer readMatchedLine2( Stream* input,
                          const char* const pattern1,
                          const char* const pattern2 ) {
  PRE06( input, input->isReadable( input ),
         pattern1, *pattern1 , pattern2, *pattern2 );
  char line[ 1024 ] = "";
  Integer result = 0;
  input->readString( input, line, COUNT( line ) );
  result = input->ok( input );

  if ( result ) {
    result = ! strcmp( line, pattern1 );

    if ( ! result ) {
      result = ! strcmp( line, pattern2 );
      result += result;
    }

    if ( ! result ) {
      failureMessage( "Invalid line in file.\n'%s'\nexpected '%s' or '%s'\n",
                      line, pattern1, pattern2 );
    }
  }

  POST0( IN4( result, 0, 1, 2 ) );
  return result;
}



/******************************************************************************
PURPOSE: skipWords - Skip over a given number of words.
INPUTS:  const char* string  String to scan.
         Integer count       Number of words to skip.
RETURNS: const char*  if skipped, else 0.
******************************************************************************/

const char* skipWords( const char* string, Integer count ) {

  PRE02( string, count > 0 );

  const char* result = string;
  Integer counter = 0;

  do {

    /* Skip any blanks before a possible word: */

    while ( isspace( *result ) ) {
      ++result;
    }

    if ( *result ) {
      ++counter; /* Count the word. */

      /* Skip the counted word: */

      do {
        ++result;
      } while ( AND2( *result, ! isspace( *result ) ) );
    }

  } while ( AND2( counter < count, *result ) );

  if ( *result ) {

    /* Skip any blanks after the last word: */

    while ( isspace( *result ) ) {
      ++result;
    }
  }

  if ( OR2( *result == '\0', counter != count ) ) {
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readTimestamp - Read/check/return a timestamp from input.
INPUTS:  Stream* input           Stream to read from.
OUTPUTS: Stream* input           Stream read from.
         UTCTimestamp timestamp  Timestamp read: "2005-08-26T23:00:00-0000".
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readTimestamp( Stream* input, UTCTimestamp timestamp ) {
  PRE03( input, input->isReadable( input ), timestamp );
  Integer result = 0;
  input->readString( input, timestamp, sizeof (UTCTimestamp) / sizeof (char));

  if ( AND2( input->ok( input ), skipInputLines( input, 1 ) ) ) {
    result = isValidUTCTimestamp( timestamp );

    if ( ! result ) {
      failureMessage( "Invalid timestamp read '%s'.", timestamp );
    }
  }

  if ( ! result ) {
    memset( timestamp, 0, sizeof (UTCTimestamp) );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        isValidUTCTimestamp( timestamp ),
                        ! timestamp[ 0 ] ) );
  return result;
}



/******************************************************************************
PURPOSE: readTimestamps - Read/check/return first/last timestamp from input.
INPUTS:  Stream* input           Stream to read from.
OUTPUTS: Stream* input           Stream read from.
         UTCTimestamp firstTimestamp  E.g., "2006-07-03T00:00:00-0000".
         UTCTimestamp lastTimestamp   E.g., "2006-07-04T23:59:59-0000".
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readTimestamps( Stream* input,
                        UTCTimestamp firstTimestamp,
                        UTCTimestamp lastTimestamp ) {

  PRE04( input, input->isReadable( input ), firstTimestamp, lastTimestamp );
  Integer result = 0;
  input->readWord( input, firstTimestamp,
                   sizeof (UTCTimestamp) / sizeof (char) );

  if ( input->ok( input ) ) {

    if ( ! isValidUTCTimestamp( firstTimestamp ) ) {
      failureMessage( "Invalid timestamp read '%s'.", firstTimestamp );
    } else {
      input->readWord( input, lastTimestamp,
                       sizeof (UTCTimestamp) / sizeof (char) );

      if ( input->ok( input ) ) {

        if ( ! isValidUTCTimestamp( lastTimestamp ) ) {
          failureMessage( "Invalid timestamp read '%s'.", lastTimestamp );
        } else if ( strcmp( firstTimestamp, lastTimestamp ) > 0 ) {
          failureMessage( "Unordered timestamps read '%s' to '%s'.",
                          firstTimestamp, lastTimestamp );
        } else {
          result = skipInputLines( input, 1 );
        }
      }
    }
  }

  if ( ! result ) {
    memset( firstTimestamp, 0, sizeof (UTCTimestamp) );
    memset( lastTimestamp,  0, sizeof (UTCTimestamp) );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        AND3( isValidUTCTimestamp( firstTimestamp ),
                              isValidUTCTimestamp( lastTimestamp ),
                              strcmp( firstTimestamp, lastTimestamp ) <= 0 ),
                        AND2( firstTimestamp[ 0 ] == '\0',
                              lastTimestamp[  0 ] == '\0' ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readDimensions - Read/check dimensions from input.
INPUTS:  Stream* input  Stream to read from.
         Integer count  Number of dimensions to read.
OUTPUTS: Stream* input  Stream read from.
         Integer dimensions[]  Dimensions read or 0.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readDimensions( Stream* input, Integer count, Integer dimensions[] ) {

  PRE04( input, input->isReadable( input ), count > 0, dimensions );

  Integer result = 0;

  if ( skipInputLines( input, 1 ) ) { /* Skip '# Dimensions: comment line. */
    char word[ 40 ] = "";
    Integer index = 0;

    do {
      input->readWord( input, word, COUNT( word ) );

      if ( ! input->ok( input ) ) {
        index = count;
      } else {
        dimensions[ index ] = atoI( word );

        if ( dimensions[ index ] < 1 ) {
          failureMessage( "Invalid dimension read '%s'.", word );
          index = count;
        }
      }

      ++index;
    } while ( index < count );

    result = index == count;

    if ( result ) {
      input->readString( input, word, COUNT( word ) ); /* Newline. */
      result = input->ok( input );
    }
  }

  if ( ! result ) {
    memset( dimensions, 0, count * sizeof *dimensions );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        GT_ZERO2( dimensions[ 0 ], dimensions[ count - 1 ] ),
                        IS_ZERO2( dimensions[ 0 ], dimensions[ count - 1 ])));
  return result;
}



/******************************************************************************
PURPOSE: readSubsetIndices - Read/check subset indices from input.
INPUTS:  Stream* input  Stream to read from.
OUTPUTS: Stream* input  Stream read from.
         Integer indices[ 8 ]  Indices read or 0.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readSubsetIndices( Stream* input, Integer indices[ 8 ] ) {

  PRE03( input, input->isReadable( input ), indices );

  const Integer count = 8;
  Integer minimum = 0;
  Integer result = 0;

  if ( skipInputLines( input, 1 ) ) { /* Skip '# subset indices: line. */
    char word[ 40 ] = "";
    Integer index = 0;

    do {
      input->readWord( input, word, COUNT( word ) );

      if ( ! input->ok( input ) ) {
        index = count;
      } else {
        indices[ index ] = atoI( word );

        if ( indices[ index ] < minimum ) {
          failureMessage( "Invalid subset index read '%s'.", word );
          index = count;
        }
      }

      ++index;

      if ( index > 1 ) {
        minimum = 1;
      }

    } while ( index < count );

    result = index == count;

    if ( result ) {
      input->readString( input, word, COUNT( word ) ); /* Newline. */
      result = input->ok( input );
    }
  }

  if ( ! result ) {
    memset( indices, 0, count * sizeof *indices );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        AND2( GE_ZERO2( indices[ 0 ], indices[ 1 ] ),
                              GT_ZERO6( indices[ 2 ], indices[ 3 ],
                                        indices[ 4 ], indices[ 5 ],
                                        indices[ 6 ], indices[ 7 ] ) ),
                        IS_ZERO8( indices[ 0 ], indices[ 1 ],
                                  indices[ 2 ], indices[ 3 ],
                                  indices[ 4 ], indices[ 5 ],
                                  indices[ 6 ], indices[ 7 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readVariablesAndUnits - Read/check variables and units from input.
INPUTS:  Stream* input  Stream to read from.
         Integer count  Number of variables to read.
OUTPUTS: Stream* input  Stream read from.
         Name    variables[]  Variables read.
         Name    units[]      Units read.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readVariablesAndUnits( Stream* input, Integer count,
                               Name variables[], Name units[] ) {

  PRE05( input, input->isReadable( input ), count > 0, variables, units );

  Integer result = 0;

  if ( skipInputLines( input, 1 ) ) { /* Skip '# Variable names:' line. */
    const Integer length = COUNT( *variables );
    Integer index = 0;

    do {
      input->readWord( input, variables[ index ], length );

      if ( ! input->ok( input ) ) {
        index = count;
      }

      ++index;
    } while ( index < count );

    if ( index == count ) {

      if ( skipInputLines( input, 2 ) ) { /* Skip '\n' '# Variable units:'*/
        index = 0;

        do {
          input->readWord( input, units[ index ], length );

          if ( ! input->ok( input ) ) {
            index = count;
          }

          ++index;
        } while ( index < count );

        result = index == count;
      }
    }
  }

  if ( result ) {
    char word[ 8 ] = "";
    input->readString( input, word, COUNT( word ) ); /* Newline. */
    result = input->ok( input );
  }

  if ( ! result ) {
    memset( variables, 0, count * sizeof *variables );
    memset( units,     0, count * sizeof *units );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        AND4( variables[ 0 ][ 0 ],
                              variables[ count - 1 ][ 0 ],
                              units[ 0 ][ 0 ],
                              units[ count - 1 ][ 0 ] ),
                        IS_ZERO4( variables[ 0 ][ 0 ],
                                  variables[ count - 1 ][ 0 ],
                                  units[ 0 ][ 0 ],
                                  units[ count - 1 ][ 0 ] ) ) );

  return result;
}



/******************************************************************************
PURPOSE: readDomain - Read/check domain from input.
INPUTS:  Stream* input             Stream to read from.
OUTPUTS: Stream* input             Stream read from.
         Real    domain[ 2 ][ 2 ]  Domain read.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readDomain( Stream* input, Real domain[ 2 ][ 2 ] ) {

  PRE03( input, input->isReadable( input ), domain );

  const char* const message = "Invalid domain value";
  Integer result = skipInputLines( input, 1 ); /* Skip '# Domain:' line. */

  result = AND2( result,
                 readReal( input, -180.0, 180.0, message,
                           &domain[ LONGITUDE ][ MINIMUM ] ) );

  result = AND2( result,
                 readReal( input, -90.0, 90.0, message,
                           &domain[ LATITUDE ][ MINIMUM ] ) );

  result = AND2( result,
                 readReal( input, domain[ LONGITUDE ][ MINIMUM ], 180.0,
                           message, &domain[ LONGITUDE ][ MAXIMUM ] ) );

  result = AND2( result,
                 readReal( input, domain[ LATITUDE ][ MINIMUM ], 90.0,
                           message, &domain[ LATITUDE ][ MAXIMUM ] ) );

  result = AND2( result, skipInputLines( input, 1 ) ); /* Read '\n'. */

  if ( ! result ) {
    memset( domain, 0, 2 * 2 * sizeof domain[ 0 ][ 0 ] );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        AND4( IN_RANGE( domain[ LONGITUDE ][ MINIMUM ],
                                        -180.0, 180.0 ),
                              IN_RANGE( domain[ LATITUDE ][ MINIMUM ],
                                        -90.0, 90.0 ),
                              IN_RANGE( domain[ LONGITUDE ][ MAXIMUM ],
                                        domain[ LONGITUDE ][ MINIMUM ],
                                        180.0 ),
                              IN_RANGE( domain[ LATITUDE ][ MAXIMUM ],
                                        domain[ LATITUDE ][ MINIMUM ],
                                        90.0 ) ),
                        IS_ZERO4( domain[ LONGITUDE ][ MINIMUM ],
                                  domain[ LONGITUDE ][ MAXIMUM ],
                                  domain[ LATITUDE  ][ MINIMUM ],
                                  domain[ LATITUDE  ][ MAXIMUM ] ) ) );

  return result;
}



/******************************************************************************
PURPOSE: readReal - Read/check real value from input.
INPUTS:  Stream* input        Stream to read from.
         Real minimum         Minimum valid value.
         Real maximum         Maximum valid value.
         const char* message  Failure message if invalid value is read.
OUTPUTS: Stream* input  Stream read from.
         Real*   value  Value read or 0 if failed.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer readReal( Stream* input, Real minimum, Real maximum,
                  const char* message, Real* value ) {

  PRE05( input, input->isReadable( input ), minimum <= maximum, message,
         value );

  Integer result = 0;
  char word[ 40 ] = "";
  *value = 0.0;

  input->readWord( input, word, COUNT( word ) );

  if ( input->ok( input ) ) {
    const Real valueRead = atoR( word );

    if ( IN_RANGE( valueRead, minimum, maximum ) ) {
      *value = valueRead;
      result = 1;
    } else {
      failureMessage( "%s '%s'.", message, word );
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        IN_RANGE( *value, minimum, maximum ),
                        *value == 0.0 ) );
  return result;
}



/******************************************************************************
PURPOSE: compress64BitIntegerValues - Copy/compress 64-bit floating-point values to
         32-bit values.
INPUTS:  Integer* array  Array of 64-bit values to expand in-place.
         Integer count   Number of values in array.
OUTPUTS: Integer* array  Compressed array of 32-bit values.
NOTES:   Values in the array are overwritten to the first half and will have
         lost precision and will be clamped to 32-bit range. UGLY.
******************************************************************************/

void compress64BitIntegerValues( Integer* array, Integer count ) {
  PRE02( array, count > 0 );
  int* destination = (int*) array;

  while ( count-- ) {
    const int value = CLAMPED_TO_RANGE( *array, INT_MIN, INT_MAX );
    *destination++ = value;
    ++array;
  }
}



/******************************************************************************
PURPOSE: replicateRealValue - Copy value to each item of output.
         32-bit values.
INPUTS:  Integer count  Number of values in array.
         Real value     value to assign to each item of array.
OUTPUTS: Real* array    Array initialized with value.
******************************************************************************/

void replicateRealValue( Real* array, Integer count, Real value ) {
  PRE03( array, count > 0, ! isNan( value ) );
  Integer index = 0;

  do {
    array[ index ] = value;
    ++index;
  } while ( index < count );

  POST02( minimumItem( array, count ) == value,
          minimumItem( array, count ) == value );

}



/******************************************************************************
PURPOSE: replicateIntValue - Copy value to each item of output.
         32-bit values.
INPUTS:  Integer count    Number of values in array.
         int value    value to assign to each item of array.
OUTPUTS: int* array   Array initialized with value.
******************************************************************************/

void replicateIntValue( int* array, Integer count, int value ) {
  PRE02( array, count > 0 );
  Integer index = 0;

  do {
    array[ index ] = value;
    ++index;
  } while ( index < count );

  POST02( array[ 0 ] == value, array[ count - 1 ] == value );
}



/******************************************************************************
PURPOSE: replaceMissingValues - Replace values < -9999.0 with -9999.0.
INPUTS:  Integer count    Number of values in array.
         Real    values[] Values to scan.
OUTPUTS: Real    values[] Updated values.
******************************************************************************/

void replaceMissingValues( Integer count, Real values[] ) {

  PRE03( count > 0, values, isNanFree( values, count ) );

  Integer index = 0;

  do {

    if ( values[ index ] < -9999.0 ) {
      values[ index ] = -9999.0;
    }

    ++index;
  } while ( index < count );

  POST02( isNanFree( values, count ), minimumItem( values, count ) >= -9999.0);
}



/******************************************************************************
PURPOSE: fractionalHours - Fractional hours from yyyydddhhmmStart to
         yyyydddhhmmNow.
INPUTS:  Integer yyyydddhhmmStart  First timestamp.
         Integer yyyydddhhmmNow    Current timestamp.
RETURNS: Real fractional hours from start to now.
******************************************************************************/

Real fractionalHours(Integer yyyydddhhmmStart, Integer yyyydddhhmmNow) {

  PRE03( isValidTimestamp( yyyydddhhmmStart ),
         isValidTimestamp( yyyydddhhmmNow ),
          yyyydddhhmmStart <= yyyydddhhmmNow );

  Real result = 0.0;
  Integer yyyydddhhmm = yyyydddhhmmStart;
  Integer hours = 0;
  Integer minutes = 0;

  while ( yyyydddhhmm < yyyydddhhmmNow ) {
    incrementTimestamp( &yyyydddhhmm );
    ++hours;
  }

  if ( yyyydddhhmm > yyyydddhhmmNow ) {
    minutes = yyyydddhhmmNow % 100;
    minutes -= yyyydddhhmm % 100;

    if ( minutes < 0 ) {
      minutes = -minutes;
    }

    --hours;
  }

  result = hours + minutes / 60.0;
  POST0( result >= 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: underscoreToSpace - Convert underscores to spaces.
INPUTS:  char* string  String to convert.
OUTPUTS: char* string  Converted string.
******************************************************************************/

void underscoreToSpace( char* string ) {
  PRE0( string );
  char* s = string;

  while ( *s ) {

    if ( *s == '_' ) {
      *s = ' ';
    }

    ++s;
  }

  POST0( strchr( string, '_' ) == 0 );
}



/******************************************************************************
PURPOSE: removeTrailingNewline - Remove possible '\n' from end of string.
INPUTS:  char* string  String to convert.
OUTPUTS: char* string  Converted string.
******************************************************************************/

void removeTrailingNewline( char* string ) {
  PRE0( string );
  Integer last = strlen( string ) - 1;

  while ( AND2( last >= 0, string[ last ] == '\n' ) ) {
    string[ last ] = '\0';
    --last;
  }
}



/******************************************************************************
PURPOSE: appendToLine - Append to line, if possible.
INPUTS:  Line  line           Line to append to.
         const const* string  String to append.
OUTPUTS: Line  line           Line with string appended, if there is space.
******************************************************************************/

void appendToLine( Line line, const char* string ) {
  PRE02( line, string );
  const size_t lengthLimit = sizeof (Line) / sizeof *line;
  const size_t lineLength = strlen( string );
  const size_t stringLength = strlen( string );

  if ( lineLength + stringLength < lengthLimit ) {
    strncat( line, string, stringLength );
  }

  line[ lengthLimit - 1 ] = '\0';
  POST0( strlen( line ) < sizeof (Line) / sizeof *line );
}



/******************************************************************************
PURPOSE: expandString - Copy a string to padded/truncated form.
INPUTS:  const char* source  String to copy.
         Integer length      Length to truncate/pad copy.
OUTPUTS: char* copy          Padded/truncated copy.
******************************************************************************/

void expandString( char* copy, const char* source, Integer length ) {
  PRE04( copy, source, *source, length > 0 );
  Integer index = 0;

  do {
    CHECK2( index < length, source[ index ] );
    copy[ index ] = source[ index ];
    ++index;
  } while ( AND2( index < length, source[ index ] ) );

  while ( index < length ) {
    copy[ index ] = ' ';
    ++index;
  }

  CHECK( index == length );
  copy[ index ] = '\0';

  POST0( strlen( copy ) == length );
}



/******************************************************************************
PURPOSE: wordsInString - Count whitespace-delimited words in string.
INPUTS:  const char* const string  String to scan.
RETURNS: Integer Number of words in string.
******************************************************************************/

Integer wordsInString( const char* const string ) {
  PRE0( string );
  Integer result = 0;
  Integer index = 0;
  Integer inWord = 0;

  /* Skip leading space: */

  while ( isspace( string[ index ] ) ) {
    ++index;
  }

  /* Count remaining words (non-space sequences): */

  while ( string[ index ] ) {
    const char c = string[ index ];

    if ( isspace( c ) ) {
      inWord = 0;
    } else if ( ! inWord ) {
      inWord = 1;
      ++result;
    }

    ++index;
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: timeData - Expand time data into contiguous storage.
INPUTS:  Integer timesteps                  Number of timesteps.
         Integer hoursPerTimestep           Hours per timestep: 1, 24, etc.
         Integer totalPoints                Total number of points.
         const Integer points[ timesteps ]  Number of points per timestep.
OUTPUTS: Real output[ totalPoints ]         0-based sequential time (hour).
******************************************************************************/

void timeData( Integer timesteps, Integer hoursPerTimestep, Integer totalPoints,
               const Integer points[], Real output[] ) {

  PRE05( GT_ZERO3( timesteps, totalPoints, hoursPerTimestep ),
         points,
         minimumItemI( points, timesteps ) >= 0,
         maximumItemI( points, timesteps ) <= totalPoints,
         output );

  Integer outputIndex = 0;
  Integer timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    const Integer hours = timestep * hoursPerTimestep;
    Integer count = points[ timestep ];

    while ( count-- ) {
      output[ outputIndex++ ] = hours;
    }
  }

  CHECK( outputIndex == totalPoints );

  POST03( isNanFree( output, totalPoints ),
          minimumItem( output, totalPoints ) >= 0.0,
          maximumItem( output, totalPoints ) <= timesteps * hoursPerTimestep );
}



/******************************************************************************
PURPOSE: readNotes - Read track notes from a stream.
INPUTS:  Stream* input        Stream to read from.
         Integer count        Number of notes to read (e.g., tracks).
OUTPUTS: Note notes[ count ]  Track notes. "MD20060703014:FRANKFURT->ATLANTA".
******************************************************************************/

void readNotes( Stream* input, Integer count, Note notes[] ) {
  PRE05( input, input->ok( input ), input->isReadable( input ), count > 0,
         notes );
  const size_t size = sizeof *notes;
  int index = 0;

  for ( index = 0; AND2( index < count, input->ok( input ) ); ++index ) {
    input->readBytes( input, notes[ index ], size ); /* Read '\n'. */
    notes[ index ][ size - 1 ] = '\0'; /* Change '\n' to '\0'. */
    trimTrailingWhitespace( size, notes[ index ] );

    if ( notes[ index ][ 0 ] == '\0' ) {
      notes[ index ][ 0 ] = '?';
    }

    CHECK( notes[ index ][ size - 1 ] == '\0' );
  }

  POST0( IMPLIES( input->ok( input ),
                  AND4( strlen( notes[ 0 ] ) > 0,
                        strlen( notes[ 0 ] ) <
                          sizeof notes[ 0 ] / sizeof *notes[ 0 ],
                        strlen( notes[ count - 1 ] ) > 0,
                        strlen( notes[ count - 1 ] ) <
                          sizeof notes[ 0 ] / sizeof *notes[ 0 ] ) ) );
}



/******************************************************************************
PURPOSE: readRegriddedNotes - Read regridded track notes from a stream.
INPUTS:  Stream* input         Stream to read from.
         Integer count         Number of notes to read
                              (e.g., totalRegriddedPoints).
OUTPUTS: RegriddedNote regriddedNotes[ count ]  Regridded track notes. E.g.,
                               "MD20060703014:FRANKFURT->ATLANTA,...".
******************************************************************************/

void readRegriddedNotes( Stream* input, Integer count,
                         RegriddedNote regriddedNotes[] ) {
  PRE05( input, input->ok( input ), input->isReadable( input ), count > 0,
         regriddedNotes );
  const size_t size = sizeof *regriddedNotes;
  int index = 0;

  for ( index = 0; AND2( index < count, input->ok( input ) ); ++index ) {
    input->readBytes( input, regriddedNotes[ index ], size ); /* Read '\n'. */
    regriddedNotes[ index ][ size - 1 ] = '\0'; /* Change '\n' to '\0'. */
    trimTrailingWhitespace( size, regriddedNotes[ index ] );

    if ( regriddedNotes[ index ][ 0 ] == '\0' ) {
      regriddedNotes[ index ][ 0 ] = '?';
    }

    CHECK( regriddedNotes[ index ][ size - 1 ] == '\0' );
  }

  POST0( IMPLIES( input->ok( input ),
                  AND4( strlen( regriddedNotes[ 0 ] ) > 0,
                        strlen( regriddedNotes[ 0 ] ) <
                          sizeof regriddedNotes[0] / sizeof *regriddedNotes[0],
                        strlen( regriddedNotes[ count - 1 ] ) > 0,
                        strlen( regriddedNotes[ count - 1 ] ) <
                        sizeof regriddedNotes[0]/sizeof *regriddedNotes[0])));
}



/******************************************************************************
PURPOSE: writeRegriddedNotes - Write regridded track notes to a stream.
INPUTS:  Stream* output        Stream to write to.
         Integer count         Number of regridded notes to write
                               (e.g., totalRegriddedPoints).
         const RegriddedNote regriddedNotes[ count ]  Regridded track notes.
                               E.g., "MD20060703014:FRANKFURT->ATLANTA,...".
******************************************************************************/

void writeRegriddedNotes( Stream* output, Integer count,
                          const RegriddedNote regriddedNotes[] ) {
  PRE05( output, output->ok( output ), output->isWritable( output ), count > 0,
         regriddedNotes );
  int index = 0;

  for ( index = 0; AND2( index < count, output->ok( output ) ); ++index ) {
    char format[ 80 ] = "";
    memset( format, 0, sizeof format );
    sprintf( format, "%%-%lus\n", sizeof regriddedNotes[ 0 ] - 1 );
    CHECK( format[ sizeof format / sizeof *format - 1 ] == '\0' );
    output->writeString( output, format, regriddedNotes[ index ] );
  }
}



/******************************************************************************
PURPOSE: expandNotes - Expand track notes to a buffer.
INPUTS:  Integer count         Number of notes to expand
         const Note notes[ count ]  Track notes.
                               E.g., "MD20060703014:FRANKFURT->ATLANTA,...".
OUTPUTS: char buffer[ count * sizeof (Note) + 1 ]  Expanded '\0' to ' '.
******************************************************************************/

void expandNotes( Integer count, Note notes[], char buffer[] ) {
  PRE03( count > 0, notes, buffer );
  const Integer size = sizeof notes[ 0 ];
  char* b = buffer;
  int index = 0;

  for ( index = 0; index < count; ++index, b += size ) {
    char format[ 80 ] = "";
    memset( format, 0, sizeof format );
    sprintf( format, "%%-%llds", size);
    CHECK( format[ sizeof format / sizeof *format - 1 ] == '\0' );
    sprintf( b, format, notes[ index ] );
    CHECK( buffer[ count * sizeof *notes ] == '\0' );
  }

  POST02( buffer[ count * sizeof *notes ] == '\0',
          strlen( buffer ) == count * sizeof *notes );
}



/******************************************************************************
PURPOSE: expandRegriddedNotes - Expand regridded track notes to a buffer.
INPUTS:  Integer count         Number of regridded notes to expand
         const RegriddedNote regriddedNotes[ count ]  Regridded track notes.
                               E.g., "MD20060703014:FRANKFURT->ATLANTA,...".
OUTPUTS: char buffer[count * sizeof (RegriddedNote) + 1]  Expanded '\0' to ' '.
******************************************************************************/

void expandRegriddedNotes( Integer count,
                           RegriddedNote regriddedNotes[],
                           char buffer[] ) {
  PRE03( count > 0, regriddedNotes, buffer );
  const Integer size = sizeof regriddedNotes[ 0 ];
  char* b = buffer;
  int index = 0;

  for ( index = 0; index < count; ++index, b += size ) {
    char format[ 80 ] = "";
    memset( format, 0, sizeof format );
    sprintf( format, "%%-%llds", size );
    CHECK( format[ sizeof format / sizeof *format - 1 ] == '\0' );
    sprintf( b, format, regriddedNotes[ index ] );
    CHECK( buffer[ count * sizeof *regriddedNotes ] == '\0' );
  }

  POST02( buffer[ count * sizeof *regriddedNotes ] == '\0',
          strlen( buffer ) == count * sizeof *regriddedNotes );
}



/******************************************************************************
PURPOSE: aggregateData - Mean-aggregate regridded data points.
INPUTS:  Integer timestepsPerAggregation  Number of timesteps for each
                                          aggregation. E.g., 24 for daily or
                                          totalpoints = sum(points[]).
         Integer isVector2                Are there 2 data variables?
         Integer timesteps                Number of timesteps of data.
         Integer points[ timesteps ]      Number of data points per timestep.
         Real longitudes[  totalPoints ]  Longitudes of points.
         Real latitudes[   totalPoints ]  Latitudes  of points.
         Real elevations[  totalPoints ]  Elevations of points or 0.
         Integer columns[  totalPoints ]  1-based grid columns of points.
         Integer rows[     totalPoints ]  1-based grid rows    of points.
         Integer layers[   totalPoints ]  1-based grid layers of points or 0.
         Real data[        totalPoints ]  Data of points.
OUTPUTS: Integer* totalOutputPoints       Total number of aggregated/output
                                            points = sum(points[]).
         Integer points[   result ] Reduced points per output timestep.

         Real longitudes[  totalOutputPoints ]  Longitudes of aggregated points
         Real latitudes[   totalOutputPoints ]  Latitudes  of aggregated points
         Real elevations[  totalOutputPoints ]  Elevations of aggregated points
                                                or 0.
         Integer columns[  totalOutputPoints ]  1-based columns of aggregated
                                                points.
         Integer rows[     totalOutputPoints ]  1-based rows of aggregated
                                                points.
         Integer layers[   totalOutputPoints ]  1-based layers of aggregated
                                                points or 0.
         Real data[        totalOutputPoints ]  Aggregated data.
         RegriddedNote notes[ totalOutputPoints ] Notes of aggregated data or 0
RETURNS: Integer total number of output timesteps or 0 if failed and a message
         is printed to stderr.
******************************************************************************/

Integer aggregateData( Integer timestepsPerAggregation,
                       Integer isVector2,
                       Integer timesteps,
                       Integer points[],
                       Real longitudes[],
                       Real latitudes[],
                       Real elevations[],
                       Integer columns[],
                       Integer rows[],
                       Integer layers[],
                       Real data[],
                       RegriddedNote notes[],
                       Integer* totalOutputPoints ) {

  PRE011( timestepsPerAggregation > 0, IS_BOOL( isVector2 ), timesteps > 0,
          points, longitudes, latitudes, columns, rows, data,
          IMPLIES_ELSE( elevations, layers, ! layers ),
          totalOutputPoints );

  Integer result = 0;
  Integer output = 0;
  Integer input = 0;
  Integer timestep = 0;
  const Integer totalRegriddedPoints = sum( timesteps, points );
  const Integer maximumLayers =
    layers ? maximumItemI( layers, totalRegriddedPoints ) : 1;
  const Integer maximumRows = maximumItemI( rows, totalRegriddedPoints );
  const Integer maximumColumns = maximumItemI( columns, totalRegriddedPoints );
  const Integer maximumRowsTimesMaximumColumns = maximumRows * maximumColumns;
  const Integer cellCount = maximumLayers * maximumRowsTimesMaximumColumns;
  Cell* cells = NEW_ZERO( Cell, cellCount );
  DEBUG( fprintf( stderr,
                "maximumRows = %lld, maximumColumns = %lld, cellCount = %lld\n",
                 maximumRows, maximumColumns, cellCount ); )

  if ( cells ) {
    Real* data2 = isVector2 ? data + totalRegriddedPoints : 0;

    while ( timestep < timesteps ) {
      const Integer timestepsRemaining = timesteps - timestep;
      const Integer timestepsToAggregate =
        MIN( timestepsPerAggregation, timestepsRemaining );
      const Integer inputPoints = sum( timestepsToAggregate, points + timestep);
      const Integer inputEnd = input + inputPoints;
      const Integer output0 = output;
      CHECK( inputEnd <= totalRegriddedPoints );
      DEBUG( fprintf( stderr, "  inputPoints = %lld, initializeCells()...\n",
                     inputPoints ); )
      initializeCells( cellCount, cells );
      DEBUG( fprintf( stderr, "  mapping to cells...\n" ); )

      /* Map each point within the timestep to a grid cell and aggregate: */

      for ( ; input < inputEnd; ++input ) {
        const Integer column = columns[ input ] - 1;
        const Integer row    = rows[    input ] - 1;
        const Integer layer  = layers ? layers[ input ] - 1: 0;
        const Integer cellIndex =
          layer * maximumRowsTimesMaximumColumns + row * maximumColumns +column;
        Cell* const cell = cells + cellIndex;
        CHECK( cellIndex < cellCount );
        cell->count += 1;
        cell->value += data[ input ];

        if ( data2 ) {
          cell->value2 +=  data2[ input ];
        }

        if ( cell->index == -1 ) {
          cell->index = input;
          ++output;
        }
      }

      DEBUG( fprintf( stderr, "  copyAggregatedData( output0 = %lld )...\n",
                      output0 ); )

      CHECKING( const Integer pointsInTimestep = )
        copyAggregatedData( cellCount, cells,
                            output0, longitudes, latitudes, elevations,
                            columns, rows, layers, data, data2, notes );
      CHECK2( pointsInTimestep == output, output < totalRegriddedPoints );
      points[ result ] = output - output0;
      ++result;
      timestep += timestepsPerAggregation;
    }
  }

  FREE( cells );
  *totalOutputPoints = output;

  DEBUG( fprintf( stderr, "result = %lld, totalOutputPoints = %lld\n",
                  result, *totalOutputPoints ); )
  POST02( result == timesteps / timestepsPerAggregation +
                    ( ( timesteps % timestepsPerAggregation ) != 0 ),
          *totalOutputPoints > 0 );
  return result;
}



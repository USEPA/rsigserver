#ifndef HELPERS_H
#define HELPERS_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Helpers.h - Declare general utility helper routines.

NOTES:

HISTORY: 2007/12, plessel.todd@epa.gov, Created.

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Utilities.h> /* For Integer, Real, UTCTimestamp, Stream. */

/*================================== MACROS =================================*/

#ifdef MAX
#undef MAX
#endif
#define MAX( a, b ) ( (a) > (b) ? (a) : (b) )

#ifdef MIN
#undef MIN
#endif
#define MIN( a, b ) ( (a) < (b) ? (a) : (b) )

#define COUNT( array ) ( sizeof (array) / sizeof *(array) )

enum { TWO_GB = 2147483647, BYTES_PER_NETCDF_FLOAT = 4 };

/*================================== TYPES ==================================*/

typedef char Name[ 80 ];  /* Ozone, ppb. */
typedef char Line[ 256 ]; /* Longer description. */

typedef Real (*CompareFunction)( Real, Real );
typedef Real (*ConvertFunction)( Real, Real, Real );

/*================================ FUNCTIONS ================================*/

extern CompareFunction compareFunction( const char* name );

extern ConvertFunction convertFunction( const char* name );

extern void compareFunctionNameUnits( CompareFunction comparer,
                                      ConvertFunction converter,
                                      Name name, Name units,
                                      Name otherName, Name otherUnits );

extern Integer sum( Integer count, const Integer data[] );
  
extern void scale( Real factor, Integer count, Real data[] );
  
extern Integer copyToStdout( Stream* input );

extern Integer streamFile( const char* fileName );

extern Integer readMatchedLine( Stream* input, const char* const pattern );

extern Integer readMatchedLine2( Stream* input,
                                 const char* const pattern1,
                                 const char* const pattern2 );

extern Integer skipInputLines( Stream* input, Integer count );

extern const char* skipWords( const char* string, Integer count );

extern void aggregateName( const char* const inputVariableName,
                           const Integer hoursPerTimestep,
                           Name outputVariableName );

extern Integer readTimestamp( Stream* input, UTCTimestamp timestamp );

extern Integer readTimestamps( Stream* input,
                               UTCTimestamp firstTimestamp,
                               UTCTimestamp lastTimestamp );

extern Integer readDimensions( Stream* input, Integer count,
                               Integer dimensions[] );

extern Integer readSubsetIndices( Stream* input, Integer indices[ 8 ] );

extern Integer readVariablesAndUnits( Stream* input, Integer count,
                                      Name variables[], Name units[] );

extern Integer readDomain( Stream* input, Real domain[ 2 ][ 2 ] );

extern Integer readReal( Stream* input, Real minimum, Real maximum,
                         const char* message, Real* value );

extern void compress64BitIntegerValues( Integer* array, Integer count );

extern void replicateRealValue( Real* array, Integer count, Real value );

extern void replicateIntValue( int* array, Integer count, int value );

extern void replaceMissingValues( Integer count, Real data[] );

extern Real fractionalHours(Integer yyyydddhhmmStart, Integer yyyydddhhmmNow);

extern void uppercase( char* string );

extern void lowercase( char* string );

extern void underscoreToSpace( char* string );

extern void removeTrailingNewline( char* string );

extern void appendToLine( Line line, const char* string );

extern void expandString( char* copy, const char* source, Integer length );

extern Integer wordsInString( const char* const string );

extern void timeData( Integer timesteps,
                      Integer hoursPerTimestep,
                      Integer totalPoints,
                      const Integer points[],
                      Real output[] );

extern void readNotes( Stream* input, Integer count, Note notes[] );

extern void readRegriddedNotes( Stream* input, Integer count,
                                RegriddedNote regriddedNotes[] );

extern void writeRegriddedNotes( Stream* output, Integer count,
                                 const RegriddedNote regriddedNotes[] );

extern void expandNotes( Integer count, Note notes[], char buffer[] );

extern void expandRegriddedNotes( Integer count,
                                  RegriddedNote regriddedNotes[],
                                  char buffer[] );

extern Integer aggregateData( Integer timestepsPerAggregation,
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
                              Integer* totalOutputPoints );

#ifdef __cplusplus
}
#endif

#endif /* HELPERS_H */




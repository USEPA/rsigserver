/******************************************************************************
PURPOSE: ReadData.c - Simple to use wrapper routine to read data data from
                      HRRR .grib2 files.

NOTES:   Uses grib2 library.

HISTORY: 2020-02-21 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(),fopen(),fseek(),fread(),fclose()*/
#include <stdlib.h> /* For size_t, malloc(), free(). */

#include <grib2.h> /* For gribfield, seekgb(), g2_info/getfld/free(). */

#include "ReadData.h" /* For public interface. */

/*================================== MACROS =================================*/

#define MISSING_VALUE (-9999.0)
#define IN_RANGE( x, low, high ) ( (low) <= (x) && (x) <= (high) )

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: readData - Read data-U/V data.
INPUTS:  const char* const fileName    Name of file to read.
         const int isVector2           Is this a 2D vector variable?
         const size_t count            Number of data points to read.
OUTPUTS: double data[ count * ( 1 + isVector2 ) ]  Subset data points.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readData( const char* const fileName,
              const int isVector2,
              const size_t count,
              double data[] ) {

  const double validMinimum = -1e30;
  const double validMaximum =  1e30;

  int result = 0;
  FILE* inputFile = 0;

  assert( fileName ); assert( *fileName ); assert( count ); assert( data );

  inputFile = fopen( fileName, "rb" );

  if ( ! inputFile ) {
    fprintf( stderr, "\nFailed to open file '%s' for reading.\n", fileName );
  } else {
    double* output = data;
    long iseek = 0;
    int messageCount = 1 + isVector2;

    /* Read 1 or 2 'messages' - 1st has datau field, 2nd has datav field: */

    for ( result = 1; result && messageCount-- ; ) {
      long messageSkip = 0;
      long messageLength = 0;

      seekgb( inputFile, iseek, 32000, &messageSkip, &messageLength );
      result = messageLength > 0;

      if ( result ) {
        unsigned char* message = (unsigned char *)
          malloc( messageLength * sizeof *message );
        result = message != 0;

        if ( ! message ) {
          fprintf( stderr, "\nFailed to allocate %lu bytes.\n",
                   messageLength * sizeof *message );
        } else {
          result = fseek( inputFile, messageSkip, SEEK_SET ) == 0;

          if ( ! result ) {
            fprintf( stderr, "\nFailed to seek %ld bytes into file '%s'.\n",
                     messageSkip, fileName );
          } else {
            const size_t bytesRead =
              fread( message, sizeof *message, messageLength, inputFile );
            result = bytesRead == messageLength * sizeof *message;

            if ( ! result ) {
              fprintf( stderr, "\nFailed to read %lu bytes from file '%s'.\n",
                       messageSkip, fileName );
            } else {
              iseek = messageSkip + messageLength;
              long unused_1[  3 ] = { 0, 0, 0 };
              long unused_2[ 13 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
              long unused_3 = 0;
              long fieldCount = 0;

              result =
                g2_info( message, unused_1, unused_2, &fieldCount, &unused_3 )
                == 0 && fieldCount == 1;

              if ( ! result ) {
               fprintf( stderr, "\nInvalid info in file '%s'.\n", fileName );
              } else {
                gribfield* gribField = 0;
                result = g2_getfld( message, 1, 1, 1, &gribField ) == 0 &&
                         gribField->ngrdpts == count && gribField->fld != 0;

                if ( ! result ) {
                  fprintf( stderr, "\nInvalid field in file '%s'.\n", fileName);
                } else if ( output <= data + count ) {
                  const float* const input = gribField->fld;
                  size_t index = 0;

                  for ( index = 0; index < count; ++index ) {
                    double value = input[ index ];

                    /* Replace NaN/invalid values with MISSING_VALUE. */

                    if ( ( gribField->bmap != 0 &&
                           gribField->bmap[ index ] == 0 ) ||
                         ! IN_RANGE( value, validMinimum, validMaximum ) ) {
                      value = MISSING_VALUE;
                    }

                    output[ index ] = value;
                  }

                  output += count;
                }

                if ( gribField ) {
                  g2_free( gribField ), gribField = 0;
                }
              }
            }
          }

          free( message ), message = 0;
        }
      }
    }

    fclose( inputFile ), inputFile = 0;
  }

  return result;
}


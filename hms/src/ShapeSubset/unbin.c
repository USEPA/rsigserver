/*
unbin.c - Unpack binary output from ShapeSubset.
gcc -m64 -Wall -O -o unbin unbin.c
usage: unbin file.bin ; ls -last *.shx *.shp *.dbf
example: unbin subset.bin ; ls -last *.shx *.shp *.dbf
*/

#include <stdio.h> /* For FILE, fopen(), fread(), fclose(), fscanf(). */
#include <string.h> /* For strncpy(), strncat(). */
#include <ctype.h> /* For isspace(). */

#ifdef MIN
#undef MIN
#endif
#define MIN(a, b) ((a) < (b) ? (a) : (b))

enum { BUFFER_SIZE = 1024 * 1024 };
static char buffer[ BUFFER_SIZE ] = "";

static int copyFileBytes( FILE* inputFile, const char* outputFileName,
                          size_t bytes ) {
  int result = 0;
  FILE* outputFile = fopen( outputFileName, "wb" );

  if ( outputFile ) {
    size_t bytesRemaining = bytes;

    do {
      const size_t bytesToRead = MIN( bytesRemaining, BUFFER_SIZE );
      const size_t bytesRead = fread( buffer, 1, bytesToRead, inputFile );

      if ( bytesRead ) {
        const size_t bytesWritten = fwrite( buffer, 1, bytesRead, outputFile );
        result = bytesWritten == bytesRead;
        bytesRemaining -= bytesWritten;
      } else {
        result = 0;
      }

    } while ( result && bytesRemaining );

    fclose( outputFile ), outputFile = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "Failed to copy %lu bytes to output file %s ",
            bytes, outputFileName );
    perror( "because" );
  }

  return result;
}



static int wordsInString( const char* s ) {
  int result = 0;

  while ( *s ) {

    while ( isspace( *s ) ) {
      ++s;
    }

    if ( *s ) {
      ++result;

      while ( *s && ! isspace( *s ) ) {
        ++s;
      }
    }
  }

  return result;
}



static int unpackShapeBinFile( const char* const inputFileName ) {

  FILE* inputFile = fopen( inputFileName, "rb" );
  int ok = 0;

  if ( inputFile ) {
    enum { FILE_NAME_LENGTH = 256, HEADER_LINE_LENGTH = 512 };
    char baseOutputFileName[ FILE_NAME_LENGTH ] = "";
    char headerLine[ HEADER_LINE_LENGTH ] = "";
    int shxBytes = 0;
    int shpBytes = 0;
    int dbfBytes = 0;
    int csvBytes = 0;
    ok = 0;
    memset( headerLine, 0, sizeof headerLine );
    memset( baseOutputFileName, 0, sizeof baseOutputFileName );

    if ( fgets( headerLine, sizeof headerLine / sizeof *headerLine,
                inputFile ) ) {
      const int words = wordsInString( headerLine );

      if ( ( ( words == 4 &&
               sscanf( headerLine, "%s %d %d %d\n",
                       baseOutputFileName, &shxBytes, &shpBytes, &dbfBytes )
                  == words )
             || ( words == 5 &&
                  sscanf( headerLine, "%s %d %d %d %d\n",
                          baseOutputFileName, &shxBytes, &shpBytes, &dbfBytes,
                          &csvBytes ) == words ) ) &&
           baseOutputFileName[ 0 ] &&
           shxBytes >= 0 && shpBytes >= 0 && dbfBytes > 0 &&
          ( words == 4 || csvBytes > 0 ) ) {
        char outputFileName[ FILE_NAME_LENGTH ] = "";
        memset( outputFileName, 0, sizeof outputFileName );
        strncpy( outputFileName, baseOutputFileName, FILE_NAME_LENGTH - 5 );
        strncat( outputFileName, ".shx", 4 );

        if ( shxBytes == 0 ||
             copyFileBytes( inputFile, outputFileName, shxBytes ) ) {
          strncpy( outputFileName, baseOutputFileName, FILE_NAME_LENGTH - 5 );
          strncat( outputFileName, ".shp", 4 );

          if ( shpBytes == 0 ||
               copyFileBytes( inputFile, outputFileName, shpBytes ) ) {
            strncpy( outputFileName, baseOutputFileName,FILE_NAME_LENGTH -5 );
            strncat( outputFileName, ".dbf", 4 );
            ok = copyFileBytes( inputFile, outputFileName, dbfBytes );
          }

          if ( ok && csvBytes ) {
            strncpy( outputFileName, baseOutputFileName,FILE_NAME_LENGTH -5 );
            strncat( outputFileName, ".csv", 4 );
            ok = copyFileBytes( inputFile, outputFileName, csvBytes );
          }
        }
      } else {
        fprintf( stderr, "Invalid header line in file %s\n", inputFileName );
      }
    }

    fclose( inputFile ), inputFile = 0;
  } else {
    fprintf( stderr, "Failed to open input file %s\n", inputFileName );
  }

  return ok;
}


int main( int argc, char* argv[] ) {
  int ok = 0;

  if ( argc == 2 && argv[ 1 ] != argv[ 0 ] ) {
    ok = unpackShapeBinFile( argv[ 1 ] );
  } else {
    fprintf( stderr, "\n%s - Unpack binary output from ShapeSubset.\n",
             argv[ 0 ] );
    fprintf( stderr, "usage: %s file.bin ", argv[ 0 ] );
    fprintf( stderr, "; ls -last *.shx *.shp *.dbf | head -3\n" );
    fprintf( stderr, "example: %s subset.bin ", argv[ 0 ] );
    fprintf( stderr, "; ls -last *.shx *.shp *.dbf | head -3\n" );
  }

  return ! ok;
}

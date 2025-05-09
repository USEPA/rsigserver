/*
convert_lonlats_file.c - Convert HRRR_latlon.h5 file to HRRR_lonlat.bin
gcc -Wall -g -o convert_lonlats_file convert_lonlats_file.c \
    -I../../../../include/NetCDF4 \
    -L../../../../lib/$platform \
    -lnetcdf4 -lhdf5_hl -lhdf5 -lcurl -lz -ldl -lm -lc
 convert_lonlats_file HRRR_latlon.h5 HRRR_lonlat.bin
 echo $?
 head -8 HRRR_lonlat.bin
 tail +9 HRRR_lonlat.bin | fdd conv=real8-ascii conv=swab8 | more

 2020-02-21 plessel.todd@epa.gov
*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(), fopen(), fwrite(), fclose(). */
#include <string.h> /* For strcmp(), strlen(). */

#include <netcdf.h> /* For NC*, nc_*(). */


static int openFile( const char* const fileName ) {
  int result = 0;
  int status = 0;
  assert( fileName );
  status = nc_open( fileName, NC_NOWRITE, &result );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to open NetCDF file for reading because: %s\n",
             message);
    result = -1;
  }

  return result;
}


static void closeFile( int file ) {
  const int status = nc_close( file );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to close NetCDF file because: %s\n",
             message );
  }
}


static int readFileData( const int file, const char* const variable,
                         const size_t rows, const size_t columns,
                         double data[] ) {

  int result = 0;
  int id = -1;
  int status = nc_inq_varid( file, variable, &id );

  if ( status == NC_NOERR ) {
    size_t start[ 2 ] = { 0, 0 };
    size_t count[ 2 ] = { 0, 0 };
    count[ 0 ] = rows;
    count[ 1 ] = columns;
    status = nc_get_vara_double( file, id, start, count, data );

    if ( status == NC_NOERR ) {
      result = 1;
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "%s\n", message );
  }

  return result;
}


/*
 * Is the platform big-endian (MSB: most significant byte first) or
 * little-endian (LSB: least significant byte first)?
 */

#if defined(__alpha) || \
    defined(__i386__) || defined(__i486__) || \
    defined(__i586__) || defined(__i686__) || \
    defined(__ia64__) || defined(__x86_64__)
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif


static void rotate8ByteArrayIfLittleEndian( void* array, const size_t count) {

#if IS_LITTLE_ENDIAN

  long long* const array8 = array;
  long long index = 0;
  assert( array ); assert( count > 0 );
  assert( sizeof (long long) == 8 );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const long long value = array8[ index ];
    const long long newValue =
    ( value & 0xff00000000000000LL ) >> 56 |
    ( value & 0x00ff000000000000LL ) >> 40 |
    ( value & 0x0000ff0000000000LL ) >> 24 |
    ( value & 0x000000ff00000000LL ) >>  8 |
    ( value & 0x00000000ff000000LL ) <<  8 |
    ( value & 0x0000000000ff0000LL ) << 24 |
    ( value & 0x000000000000ff00LL ) << 40 |
    ( value & 0x00000000000000ffLL ) << 56;
    array8[ index ] = newValue;
  }

#endif /* IS_LITTLE_ENDIAN */

}


enum { ROWS = 1059, COLUMNS = 1799, COUNT = ROWS * COLUMNS };
static double lonlats[ COUNT * 2 ];

int main( int argc, char* argv[] ) {
  int ok = 0;

  if ( argc == 3 &&
      strcmp( argv[ 1 ], argv[ 0 ] ) &&
      strcmp( argv[ 2 ], argv[ 0 ] ) &&
      strcmp( argv[ 1 ], argv[ 2 ] ) ) {
    const char* const inputFileName = argv[ 1 ];
    int inputFile = openFile( inputFileName );

    if ( inputFile >= 0 ) {
      ok = readFileData( inputFile, "longitude", ROWS, COLUMNS, lonlats );

      if ( ok ) {
        ok = readFileData( inputFile, "latitude", ROWS, COLUMNS, lonlats + COUNT );

        if ( ok ) {
          const char* const outputFileName = argv[ 2 ];
          FILE* outputFile = fopen( outputFileName, "wb" );

          if ( outputFile ) {
            const char* const header =
              "Content-type: application/octet-stream; charset=iso-8859-1\n"
              "# dimensions: variables rows columns\n"
              "2 %d %d\n"
              "# variable names:\n"
              "longitude latitude\n"
              "# variable units:\n"
              "deg deg\n"
              "# IEEE-754 64-bit real data[variables][rows][columns]:\n";
            const size_t headerLength = strlen( header );
            ok = fprintf( outputFile, header, ROWS, COLUMNS ) >= headerLength - 10;

            if ( ok ) {
              rotate8ByteArrayIfLittleEndian( lonlats, COUNT * 2 );
              ok = fwrite( lonlats, sizeof lonlats, 1, outputFile ) == 1;
            }

            fclose( outputFile ), outputFile = 0;
          }
        }
      }

      closeFile( inputFile ), inputFile = -1;
    }
  }

  return ! ok;
}




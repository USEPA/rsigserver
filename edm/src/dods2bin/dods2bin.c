/*
dods2bin.c - C program to convert DODS grid data to bin format.
2012-09-19 plessel.todd@epa.gov
To compile:
$ cc -o dods2bin dods2bin.c
Usage:
dods2bin name units minimum maximum < data.dods > data.bin
data outside the range [minimum, maximum] is mapped to -9999.
Usage examples:

$ wget -q -O - \
'http://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMHsstd1day.dods?sst\
[(2006-12-20T12:00:00Z):1:(2006-12-23T12:00:00Z)]\
[(0.0):1:(0.0)][(35):1:(36.0)][(285):1:(286)]' > sst.dods
$ head -12 sst.dods
Dataset {
  GRID {
    ARRAY:
      Float32 sst[time = 4][altitude = 1][latitude = 25][longitude = 25];
    MAPS:
      Float64 time[time = 4];
      Float64 altitude[altitude = 1];
      Float64 latitude[latitude = 25];
      Float64 longitude[longitude = 25];
  } sst;
} erdMHsstd1day;

$ dods2bin sst C 0 50 < sst.dods > sst.bin
$ head -7 sst.bin
Content-type: application/octet-stream; charset=iso-8859-1
# variable units:
sst C
# dimensions: timesteps z rows columns lonmin lonmax latmin latmax
4       0        25         25  -84.02304222733599204  -81.98243622733599523   32.97649710028599657   34.04116110028599707
# char yyyymmddhh[timesteps][11] and
# IEEE-754 32-bit float data[timesteps][rows][columns]:

Input formats can have 32 or 64 bit data and 32 or 64 bit coordinates and
include or exclude z dimension:

Float32 chlorophyll[time = 1][altitude = 1][latitude = 241][longitude = 241];
Float64 time[time = 1];
Float64 altitude[altitude = 1];
Float64 latitude[latitude = 241];
Float64 longitude[longitude = 241];

Float32 temperature_an[time = 2][depth = 1][latitude = 2][longitude = 2];
Float64 time[time = 2];
Float32 depth[depth = 1];
Float32 latitude[latitude = 2];
Float32 longitude[longitude = 2];

Float64 analysed_sst[time = 3][latitude = 92][longitude = 92];
Float64 time[time = 3];
Float32 latitude[latitude = 92];
Float32 longitude[longitude = 92];

Note: each of the DODS input arrays are MSB/big-endian IEEE-754 32 or 64-bit
format and are preceeded by two 4-byte words (MSB integer) containing the
count of array items that follows.
*/


#include <stdio.h>  /* For FILE, stdoerr, scanf(), printf(), fwrite(). */
#include <stdlib.h> /* For atof(), malloc(), free(). */
#include <string.h> /* For strcmp(), memset(). */


#if defined(__alpha) || \
    defined(__i386__) || defined(__i486__) || \
    defined(__i586__) || defined(__i686__) || \
    defined(__ia64__) || defined(__x86_64__)
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif

#define IN_RANGE( x, lower, upper ) ( (x) >= (lower) && (x) <= (upper) )

typedef char YYYYMMDDHH[ 11 ];

/* Forward declarations: */

static void usage( const char* program );

static int parse_header( int* timesteps, int* rows, int* columns,
                         int* data_bits, int* coordinate_bits,
                         int* has_layers, int* negative_z );

static int process_data(const int timesteps, const int rows, const int columns,
                        const int data_bits, const int coordinate_bits,
                        const int has_layers, const int negative_z,
                        const int swap_rows,
                        const char* const name, const char* const units,
                        const double minimum, const double maximum );

static float* read_data( const int timesteps,
                         const int rows, const int columns,
                         const int bits, const int swap_rows,
                         const double minimum, const double maximum );

static YYYYMMDDHH* read_time( const int timesteps );

static double read_z( const int bits );

static int read_coordinates( const size_t count, const int bits,
                             double* minimum, double* maximum );

static int write_output( const int timesteps,
                         const int rows, const int columns,
                         const char* const name, const char* const units,
                         YYYYMMDDHH yyyymmddhh[],
                         const double z,
                         const double longitude_minimum,
                         const double longitude_maximum,
                         const double latitude_minimum,
                         const double latitude_maximum,
                         float data[] );

static int skip_8_bytes( void );

static void rotate_4_byte_words_if_little_endian( const size_t count,
                                                  void* array );

static void rotate_8_byte_words_if_little_endian( const size_t count,
                                                  void* array );

static void swap_data_rows( float array[], const int timesteps,
                            const int rows, const int columns );

static int hours_of_seconds( const double seconds );

static int yyyymmddhh_of_hours( const int hours_tai );

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours );

static void decrement_yyyymmddhh( int* yyyymmddhh, const int hours );

static int days_in_month( int year, int month );



/* Read input DODS data from stdin, write output bin data to stdout: */

int main( int argc, char* argv[] ) {
  int result = 0;

  if ( argc != 6 ) {
    usage( argv[ 0 ] );
  } else {
    const char* const name  = argv[ 1 ];
    const char* const units = argv[ 2 ];
    const double minimum    = atof( argv[ 3 ] );
    const double maximum    = atof( argv[ 4 ] );
    const int swap_rows     = argv[ 5 ][ 0 ] == '1';
    int timesteps = 0;
    int rows = 0;
    int columns = 0;
    int data_bits = 0;
    int coordinate_bits = 0;
    int has_layers = 0;
    int negative_z = 0;

    if ( parse_header( &timesteps, &rows, &columns,
                       &data_bits, &coordinate_bits,
                       &has_layers, &negative_z ) ) {
      result =
        process_data( timesteps, rows, columns, data_bits, coordinate_bits,
                      has_layers, negative_z, swap_rows,
                      name, units, minimum, maximum );
    }
  }

  return ! result;
}



/* Print program usage instructions to stderr: */

static void usage( const char* program ) {
  fprintf( stderr, "\n%s - Convert DODS grid data to bin format.\n", program );
  fprintf( stderr, "usage: %s name units minimum maximum swap_rows"
                   "< input.dods > output.bin\n", program );
  fprintf( stderr, "example: %s sst C 0 50 0 < sst.dods > sst.bin\n", program );
  fprintf( stderr, "head -7 sst.bin\n\n" );
}



/* Read/check DODS header: */

static int parse_header( int* timesteps, int* rows, int* columns,
                         int* data_bits, int* coordinate_bits,
                         int* has_layers, int* negative_z ) {
  char layer_name[ 32 ] = "";
  char temp_name[ 32 ] = "";
  int layers = 0;
  int temp1 = 0;
  int temp2 = 0;
  int temp3 = 0;
  int result = 0;
  memset( layer_name, 0, sizeof layer_name );
  memset( temp_name, 0, sizeof temp_name );

  *timesteps = *rows = *columns = *data_bits = *coordinate_bits = 0;
  *has_layers = *negative_z = 0;

  result =
    scanf( "Dataset {\n"
              "  GRID {\n"
              "    ARRAY:\n"
              "      Float%d %*s = %d%31s = %d]",
          data_bits, timesteps, layer_name, &layers ) == 4;

  result = result &&
           ( *data_bits == 32 || *data_bits == 64 ) &&
           *timesteps > 0 && *layer_name && layers >= 0;

  if ( result ) {

    if ( ! strcmp( layer_name, "][latitude" ) ) { /* Read rows as layers. */
      *rows = layers;
      layers = 1;
      result = scanf( "%*s = %d];\n", columns ) == 1;
    } else { /* Read some named layer dimension: */
      *has_layers = 1;
      *negative_z = ! strcmp( layer_name, "][depth" );
      result =
        layers == 1 && scanf( "%*s = %d%*s = %d];\n", rows, columns ) == 2;
    }

    result = result && *rows > 0 && *columns > 0;

    result = result &&
      scanf( "    MAPS:\n"
             "      Float64 %*s = %d];\n", &temp1 ) == 1;
    result = result && temp1 == *timesteps;

    if ( result ) {

      if ( *has_layers ) {
        result =
          scanf("      Float%d %31s = %d];\n", &temp1, temp_name, &temp2) == 3;
        result = result &&
                 ( temp1 == 32 || temp1 == 64 ) &&
                 ! strncmp(temp_name, layer_name + 2, strlen(layer_name +2)) &&
                 temp2 == layers;
      }

      if ( result ) {
        result =
          scanf( "      Float%d %*s = %d];\n"
                 "      Float%d %*s = %d];\n"
                 "%*[^\n]\n"
                 "%*[^\n]\n"
                 "%31s\n",
                 coordinate_bits, &temp1, &temp2, &temp3, temp_name ) == 5;
        result = result &&
                 ( *coordinate_bits == 32 || *coordinate_bits == 64 ) &&
                 temp1 == *rows && temp2 == *coordinate_bits &&
                 temp3 == *columns && ! strcmp( temp_name, "Data:" );
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nRead invalid input DODS header.\n" );
  }

  return result;
}



/* Read/write data: */

static int process_data(const int timesteps, const int rows, const int columns,
                        const int data_bits, const int coordinate_bits,
                        const int has_layers, const int negative_z,
                        const int swap_rows,
                        const char* const name, const char* const units,
                        const double minimum, const double maximum ) {
  int result = 0;
  float* data =
    read_data(timesteps, rows, columns, data_bits, swap_rows, minimum, maximum);

  if ( data ) {
    YYYYMMDDHH* yyyymmddhh = read_time( timesteps );

    if ( yyyymmddhh ) {
      double z = 0.0;
      double longitude_minimum = 0.0;
      double longitude_maximum = 0.0;
      double latitude_minimum = 0.0;
      double latitude_maximum = 0.0;

      if ( has_layers ) {
        z = read_z( coordinate_bits );

        if ( negative_z && z > 0.0) {
          z = -z;
        }
      }

      if ( read_coordinates( rows, coordinate_bits,
                             &latitude_minimum, &latitude_maximum ) ) {

        if ( read_coordinates( columns, coordinate_bits,
                               &longitude_minimum, &longitude_maximum ) ) {
          result =
            write_output( timesteps, rows, columns, name, units, yyyymmddhh, z,
                          longitude_minimum, longitude_maximum,
                          latitude_minimum, latitude_maximum, data );
        }
      }

      free( yyyymmddhh ), yyyymmddhh = 0;
    }

    free( data ), data = 0;
  }


  return result;

}



/* Read/convert/filter/return data: */

static float* read_data( const int timesteps,
                         const int rows, const int columns,
                         const int bits, const int swap_rows,
                         const double minimum, const double maximum ) {
  float* result = 0;
  const size_t count = timesteps * rows * columns;
  const size_t size = bits == 32 ? 4 : 8;
  void* data = malloc( count * size );

  if ( data ) {
    memset( data, 0, count * size );

    if ( skip_8_bytes() && fread( data, size, count, stdin ) == count ) {

      if ( bits == 32 ) {
        float* fdata = data;
        size_t index = 0;
        rotate_4_byte_words_if_little_endian( count, data );

        /* Map invalid values to -9999: */

        for ( index = 0; index < count; ++index ) {
          const double value = fdata[ index ];

          if ( ! IN_RANGE( value, minimum, maximum ) ) {
            fdata[ index ] = -9999.0;
          }
        }

      } else { /* Convert 64-bit input data to 32-bit: */
        float*  fdata = data;
        double* ddata = data;
        size_t index = 0;
        rotate_8_byte_words_if_little_endian( count, data );

        /* Map invalid values to -9999: */

        for ( index = 0; index < count; ++index ) {
          double value = ddata[ index ];

          if ( ! IN_RANGE( value, minimum, maximum ) ) {
            value = -9999.0;
          }

          fdata[ index ] = value; /* Copy to lower 4-bytes of 8-byte word. */
        }

        memset( fdata + count, 0, count * sizeof (float) ); /* Zero unused. */
      }
    }

    if ( swap_rows ) {
      swap_data_rows( (float*) data, timesteps, rows, columns );
    }
  }

  result = data;
  data = 0;

  if ( ! result ) {
    fprintf( stderr, "\nFailed to allocate/read/convert data.\n" );
    perror( "because" );
  }

  return result;
}



/* Read/convert array of TAI seconds (non-const-delta) to UTC date + hour: */

static YYYYMMDDHH* read_time( const int timesteps ) {
  YYYYMMDDHH* result = 0;
  double* seconds_tai = malloc( timesteps * sizeof (double) );

  if ( seconds_tai ) {
    memset( seconds_tai, 0, timesteps * sizeof (double) );

    if ( skip_8_bytes() &&
         fread( seconds_tai, sizeof (double), timesteps, stdin) == timesteps) {
      rotate_8_byte_words_if_little_endian( timesteps, seconds_tai );
      result = malloc( timesteps * sizeof (YYYYMMDDHH) );

      if ( result ) {
        int timestep = 0;
        memset( result, 0, timesteps * sizeof (YYYYMMDDHH) );

        for ( timestep = 0; timestep < timesteps;  ++timestep ) {
          const double seconds = seconds_tai[ timestep ];
          const int hours = hours_of_seconds( seconds );
          const int yyyymmddhh = yyyymmddhh_of_hours( hours );
          sprintf( result[ timestep ], "%010d", yyyymmddhh );
        }
      }
    }

    free( seconds_tai ), seconds_tai = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to allocate/read/convert time.\n" );
    perror( "because" );
  }

  return result;
}



/* Read z data: */

static double read_z( const int bits ) {
  double result = 0.0;
  int ok = 0;

  if ( bits == 32 ) {
    float z = 0.0;

    if ( skip_8_bytes() && fread( &z, sizeof (float), 1, stdin ) == 1 ) {
      rotate_4_byte_words_if_little_endian( 1, &z );
      result = z;
      ok = 1;
    }

  } else if ( skip_8_bytes() &&
              fread( &result, sizeof (double), 1, stdin ) == 1 ) {
    rotate_8_byte_words_if_little_endian( 1, &result );
      ok = 1;
  }

  if ( ! ok ) {
    fprintf( stderr, "\nFailed to read z coordinate " );
    perror( "because" );
  }

  return result;
}



/* Read coordinate data: */

static int read_coordinates( const size_t count, const int bits,
                             double* minimum, double* maximum ) {
  int result = 0;
  const size_t size = bits == 32 ? 4 : 8;
  void* data = malloc( count * size );
  *minimum = *maximum = 0.0;

  if ( data ) {
    memset( data, 0, count * size );

    if ( skip_8_bytes() && fread( data, size, count, stdin ) == count ) {
      double first = 0.0;
      double last  = 0.0;

      if ( bits == 32 ) {
        float* fdata = data;
        rotate_4_byte_words_if_little_endian( count, data );
        first = fdata[ 0 ];
        last  = fdata[ count - 1 ];
      } else {
        double* fdata = data;
        rotate_8_byte_words_if_little_endian( count, data );
        first = fdata[ 0 ];
        last  = fdata[ count - 1 ];
      }

      /* Convert longitudes in range [0, 360] to [-180, 180]: */

      if ( first > 180.0 ) {
        first -= 360.0;
      }

      if ( last > 180.0 ) {
        last -= 360.0;
      }

      if ( first > last ) {
        const double swap_temp = first;
        first = last;
        last  = swap_temp;
      }

      if ( count > 1 ) {
        const double half_delta = 0.5 * ( last - first ) / ( count - 1 );
        *minimum = first - half_delta;
        *maximum = last  + half_delta;
        result = 1;
      } else {
        *minimum = first;
        *maximum = last;
        result = 1;
      }
    }

    free( data ), data = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to allocate/read/convert coordinates.\n" );
    perror( "because" );
  }

  return result;
}



/* Write bin-format ASCII header and binary data to stdout: */

static int write_output( const int timesteps,
                         const int rows, const int columns,
                         const char* const name, const char* const units,
                         YYYYMMDDHH yyyymmddhh[],
                         const double z,
                         const double longitude_minimum,
                         double longitude_maximum,
                         const double latitude_minimum,
                         double latitude_maximum,
                         float data[] ) {
  const size_t count = timesteps * rows * columns;
  int result = 0;
  int timestep = 0;

  /* Write ASCII header: */

  printf( "Content-type: application/octet-stream; charset=iso-8859-1\n" );
  printf( "# variable units:\n" );
  printf( "%s %s\n", name, units );
  printf( "# dimensions: "
          "timesteps z rows columns lonmin lonmax latmin latmax\n" );
  printf( "%-5d %5g %10d %10d %24.18f %24.18f %24.18f %24.18f\n",
          timesteps, z, rows, columns, longitude_minimum, longitude_maximum,
          latitude_minimum, latitude_maximum );
  printf( "# char yyyymmddhh[timesteps][11] and\n" );
  printf( "# IEEE-754 32-bit float data[timesteps][rows][columns]:\n" );

  /* Write timestamps: */

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    printf( "%s\n", yyyymmddhh[ timestep ] );
  }

  /* Write binary MSB 32-bit IEEE-754 data: */

  rotate_4_byte_words_if_little_endian( count, data ); /* Big-endian output. */
  result = fwrite( data, sizeof (float), count, stdout ) == count;

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write all %lu bytes of data.\n",
             count * sizeof (float) );
    perror( "because" );
  }

  return result;
}



/* Read and ignore 8 bytes from stdin: */

static int skip_8_bytes( void ) {
  double unused = 0.0;
  const int result = fread( &unused, 8, 1, stdin ) == 1;
  return result;
}



/* On little-endian platforms, change endianess of an array of 4-byte words: */

static void rotate_4_byte_words_if_little_endian( const size_t count,
                                                  void* array ) {

#if IS_LITTLE_ENDIAN

  int* const array4 = array;
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const unsigned int value = array4[ index ];
    const unsigned int new_value =
      ( value & 0xff000000 ) >> 24 |
      ( value & 0x00ff0000 ) >>  8 |
      ( value & 0x0000ff00 ) <<  8 |
      ( value & 0x000000ff ) << 24;
    array4[ index ] = new_value;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/* On little-endian platforms, change endianess of an array of 8-byte words: */

static void rotate_8_byte_words_if_little_endian( const size_t count,
                                                  void* array ) {
#if IS_LITTLE_ENDIAN

  unsigned long long* const array8 = (unsigned long long*) array;
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const unsigned long long value = array8[ index ];
    const unsigned long long new_value =
    ( value & 0xff00000000000000ULL ) >> 56 |
    ( value & 0x00ff000000000000ULL ) >> 40 |
    ( value & 0x0000ff0000000000ULL ) >> 24 |
    ( value & 0x000000ff00000000ULL ) >>  8 |
    ( value & 0x00000000ff000000ULL ) <<  8 |
    ( value & 0x0000000000ff0000ULL ) << 24 |
    ( value & 0x000000000000ff00ULL ) << 40 |
    ( value & 0x00000000000000ffULL ) << 56;
    array8[ index ] = new_value;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/*
 * swap_data_rows - Swap data rows in-place.
 * float array[ timesteps ][ rows ][ columns ]
 */

static void swap_data_rows( float array[], const int timesteps,
                            const int rows, const int columns ) {
  float* a = array;
  const int timestep_offset = rows * columns;
  int timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep, a += timestep_offset ) {
    int lower_row_index = 0;
    int upper_row_index = rows - 1;

    for ( ; lower_row_index < upper_row_index;
         ++lower_row_index, --upper_row_index ) {
      float* const lower_row_data = a + lower_row_index * columns;
      float* const upper_row_data = a + upper_row_index * columns;
      int column = 0;

      for ( column = 0; column < columns; ++column ) {
        float temp = lower_row_data[ column ];
        lower_row_data[ column ] = upper_row_data[ column ];
        upper_row_data[ column ] = temp;
      }
    }
  }
}



/* convert seconds TAI to whole hours: */

static int hours_of_seconds( const double seconds ) {
  const int result = seconds / 60.0 / 60.0;
  return result;
}



/* Convert hours TAI to UTC data+hh: */

static int yyyymmddhh_of_hours( const int hours_tai ) {
  int result = 1970010100; /* 1970-01-01T00. */

  if ( hours_tai < 0 ) {
    decrement_yyyymmddhh( &result, hours_tai );
  } else {
    increment_yyyymmddhh( &result, hours_tai );
  }

  if ( result / 1000000 <= 0 ) { /* If year < 1 make it 1970. */
    result += 1970000000;
  }

  return result;
}



/* Increment yyyymmddhh by hours: */

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours ) {
  int yyyy = *yyyymmddhh / 1000000;
  int mm   = *yyyymmddhh / 10000 % 100;
  int dd   = *yyyymmddhh / 100 % 100;
  int hh   = *yyyymmddhh % 100;
  int hour = 0;

  for ( hour = 0; hour < hours; ++hour ) {
    ++hh;

    if ( hh > 23 ) {
      hh = 0;
      ++dd;

      if ( dd > 28 && dd > days_in_month( yyyy, mm ) ) {
        dd = 1;
        ++mm;

        if ( mm > 12 ) {
          mm = 1;
          ++yyyy;
        }
      }
    }
  }

  *yyyymmddhh = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
}



/* Decrement yyyymmddhh by hours: */

static void decrement_yyyymmddhh( int* yyyymmddhh, const int hours ) {
  int yyyy = *yyyymmddhh / 1000000;
  int mm   = *yyyymmddhh / 10000 % 100;
  int dd   = *yyyymmddhh / 100 % 100;
  int hh   = *yyyymmddhh % 100;
  int hour = hours;

  while ( hour++ ) {
    --hh;

    if ( hh < 0 ) {
      hh = 23;
      --dd;

      if ( dd < 1 ) {
        --mm;

        if ( mm < 1 ) {
          mm = 12;
          --yyyy;
        }

        dd = days_in_month( yyyy, mm );
      }
    }
  }

  *yyyymmddhh = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
}



/* Number of days in year/month: */

static int days_in_month( int year, int month ) {
  static const int days_per_month[ 2 ][ 12 ] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Leap year. */
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
  };
  const int leap =
    month != 2 ? 0 : year % 4 == 0 && ! ( year % 100 == 0 && year % 400 != 0 );
  const int result = days_per_month[ leap ][ month - 1 ];
  return result;
}




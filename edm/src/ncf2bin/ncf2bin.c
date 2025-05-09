/*
ncf2bin.c - C program to convert ncf grid data to bin format.
2024-02-22 plessel.todd@epa.gov
To compile:
 $ gcc -Wall -O -I. -o ncf2bin ncf2bin.c libnetcdf.a ; strip ncf2bin
Usage:
ncf2bin list_file name units minimum maximum west south east north " \
 > output.bin
data outside the range [minimum, maximum] is mapped to -9999.
Usage examples:

wget -q --header "Authorization: Bearer $TOKEN" -O - \
  'https://opendap.earthdata.nasa.gov/collections/C2491756421-POCLOUD/\
   granules/Q2012033.L3m_DAY_SCI_V5.0_SSS_1deg.dap.nc?\
   dap4.ce=/l3m_data%5B120:1:121%5D%5B89:1:90%5D' > salinity_20120202.nc

ncdump salinity_20120202.nc | more

echo salinity_20120202.nc > list_file

ncf2bin list_file salinity PSU 0 50 -90.5 30.5 -89.5 31.5 \
  > salinity_20120202.bin

head -7 salinity_20120202.bin
Content-type: application/octet-stream; charset=iso-8859-1
# variable units:
salinity PSU
# dimensions: timesteps z rows columns lonmin lonmax latmin latmax
1 0 2 2 -90.5 -89.5 30.5 31.5
# char yyyymmddhh[timesteps][11] and
# IEEE-754 32-bit float data[timesteps][rows][columns]:

*/

#include <stdio.h>     /* For stderr, fprintf(), fwrite(). */
#include <stdlib.h>    /* For strtod(), malloc(), free(). */
#include <string.h>    /* For strcmp(), memset(). */
#include <ctype.h>     /* For isalpha(). */
#include <limits.h>    /* For INT_MAX. */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

#include <netcdf.h> /* For nc_*(). */

#ifndef IS_LITTLE_ENDIAN
#define IS_LITTLE_ENDIAN \
  ( \
    ( defined(__BYTE_ORDER__) && \
      defined(__ORDER_LITTLE_ENDIAN__) && \
      __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ ) || \
    defined(__x86_64__) || \
    defined(__i386__)   || \
    defined(__i486__)   || \
    defined(__i586__)   || \
    defined(__i686__)   || \
    defined(__alpha)    || \
    defined(__ia64__)   || \
    defined(__ARMEL__)  || \
    defined(__MIPSEL)   || \
    defined(_WIN32)     || \
    defined(_WIN64)     || \
    defined(_M_IX86)    || \
    defined(_M_AMD64)      \
  )
#endif

#define IN_RANGE( x, lower, upper ) ( (x) >= (lower) && (x) <= (upper) )

typedef struct {
  char* input_file_names;
  const char* name;
  const char* units;
  double minimum;
  double maximum;
  double west;
  double south;
  double east;
  double north;
  int timesteps;
  int rows;
  int columns;
  int* yyyymmddhh; /* yyyymmddhh[ data->timesteps ]. */
  float* data;     /* data[ data->timesteps * data->rows * data->columns ]. */
} Data;

static void deallocate_data( Data* const data ) {

  if ( data ) {

    if ( data->input_file_names ) {
      free( data->input_file_names ), data->input_file_names = 0;
    }

    if ( data->yyyymmddhh ) {
      free( data->yyyymmddhh ), data->yyyymmddhh = 0;
    }

    if ( data->data ) {
      free( data->data ), data->data = 0;
    }
  }
}

/* Forward declarations: */

static void usage( const char* program );

static int parse_arguments( int argc, char* argv[], Data* const data );

static char* read_file( const char* const file_name, int* const line_count );

static size_t file_size( const char* file_name );

static size_t control_m_to_space( char* string );

static int parse_timestamp( const char* const file_name );

static int is_yyyymmddhh( const int yyyymmddhh );

static int days_in_month( const int year, const int month );

static int read_all_data( Data* data );

static int read_file_dimensions( const int file_id,
                                 int* const rows,
                                 int* const columns );

static int read_file_data( const int file_id,
                           const int rows,
                           const int columns,
                           const double minimum,
                           const double maximum,
                           float* const data );

static int write_output( Data* const data );

static void rotate_4_byte_words_if_little_endian( const size_t count,
                                                  void* array );

static void swap_data_rows( float array[], const int rows, const int columns );



/* Read data from list of NetCDF files then write bin-format data to stdout: */

int main( int argc, char* argv[] ) {
  int result = 0;
  Data data;
  memset( &data, 0, sizeof data );

  if ( ! parse_arguments( argc, argv, &data) ) {
    usage( argv[ 0 ] );
  } else {
    result = read_all_data( &data );

    if ( result ) {
      result = write_output( &data );
    }
  }

  deallocate_data( &data );
  return ! result;
}



/* Print program usage instructions to stderr: */

static void usage( const char* program ) {
  fprintf( stderr, "\n%s - Convert NetCDF 2D grid data to bin format.\n",
           program );
  fprintf( stderr, "usage: %s list_file name units minimum maximum "
          "west south east north > output_file\n", program );
  fprintf( stderr, "example: %s list_file salinity PSU 0 50 "
           "-90.5 30.5 -89.5 31.5 > salinity_daily_20120202.bin\n", program );
  fprintf( stderr, "head -7 salinity_daily_20120202.bin\n\n" );
}



/* Read and check command-line arguments: */

static int parse_arguments( int argc, char* argv[], Data* const data ) {
  int result = 0;
  memset( data, 0, sizeof *data );

  if ( argc == 10 ) {
    data->input_file_names = read_file( argv[ 1 ], &data->timesteps );

    if ( data->timesteps > 0 ) {
      const size_t bytes = data->timesteps * sizeof (int);
      data->yyyymmddhh = malloc( bytes );

      if ( ! data->yyyymmddhh ) {
        fprintf( stderr, "\nFailed to allocate %lu bytes for timestamps.\n",
                 bytes );
      } else {
        memset( data->yyyymmddhh, 0, bytes );
        data->name = argv[ 2 ];

        if ( isalpha( data->name[ 0 ] ) ) {
          data->units = argv[ 3 ];

          if ( isalpha( data->units[ 0 ] ) ) {
            char* end = 0;
            data->minimum = strtod( argv[ 4 ], &end );

            if ( end != argv[ 4 ] ) {
              data->maximum = strtod( argv[ 5 ], &end );

              if ( end != argv[ 5 ] && data->maximum >= data->minimum ) {
                data->west = strtod( argv[ 6 ], &end );

                if ( end != argv[ 6 ] && IN_RANGE( data->west, -180.0, 180.0)) {
                  data->south = strtod( argv[ 7 ], &end );

                  if ( end != argv[ 7 ] && IN_RANGE( data->south, -90.0, 90.0)) {
                    data->east = strtod( argv[ 8 ], &end );

                    if ( end != argv[ 8 ] &&
                         IN_RANGE( data->east, data->west, 180.0 ) ) {
                      data->north = strtod( argv[ 9 ], &end );

                      if ( end != argv[ 9 ] &&
                           IN_RANGE( data->north, data->south, 90.0 ) ) {
                        result = 1;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return result;
}



/* Read file into a string and get line count: */

static char* read_file( const char* const file_name, int* const line_count ) {
  char* result = 0;
  const size_t size = file_size( file_name ) / sizeof (char);
  *line_count = 0;

  if ( size > 0 ) {
    const size_t bytes = ( size + 1 ) * sizeof (char);
    result = malloc( bytes );

    if ( ! result ) {
      fprintf( stderr,
              "\nCan't allocate %lu bytes to complete the requested action.\n",
               bytes );
    } else {
      FILE* file = fopen( file_name, "rb" );
      memset( result, 0, bytes );

      if ( file ) {
        const size_t itemsRead = fread( result, sizeof (char), size, file );

        if ( itemsRead != size ) {
          fprintf( stderr, "\nFailed to read entire file '%s'.\n", file_name );
          free( result ), result = 0;
        } else {
          result[ size ] = '\0'; /* Terminate string. */
          *line_count = control_m_to_space( result );
        }

        fclose( file );
        file = 0;
      }
    }
  }

  return result;
}



/* Get size in bytes of named file: */

static size_t file_size( const char* file_name ) {
  size_t result = 0;
  struct stat buf;

  if ( stat( file_name, &buf ) == -1 ) {
    fprintf( stderr, "\nFailed to determine size of file '%s'.\n", file_name );
  } else {
    result = buf.st_size;

    if ( result < 1 ) {
      fprintf( stderr, "\nNegative size of file '%s'.\n", file_name );
      result = 0;
    }
  }

  return result;
}



/* Change any control-m (\r) characters to space and return line count: */

static size_t control_m_to_space( char* string ) {
  size_t result = 0;
  char* s = string;

  while ( *s ) {
    const char c = *s;

    if ( c == '\r' ) {
      *s = ' ';
    } else if ( c == '\n' ) {
      ++result;
    }

    ++s;
  }

  return result;
}



/* Get timestamp yyyymmddhh from file name: */

static int parse_timestamp( const char* const file_name ) {
  int result = 0;
  const char* underscore = strrchr( file_name, '_' );

  if ( underscore ) {
    result = atoi( underscore + 1 );

    if ( result >= 1900 ) {

      if ( result <= 2147 ) {
        result = result * 1000000 + 10100;
      } else if ( result < 300012 ) {
        result = result * 10000 + 100;
      } else if ( result < 2147123123 ) {
        result *= 100;
      }
    }
  }

  if ( ! is_yyyymmddhh( result ) ) {
    fprintf( stderr, "\nInvalid timestamp (%d) of file '%s'.\n",
             result, file_name );
    result = 0;
  }

  return result;
}



/* Is yyyymmddhh valid? */

static int is_yyyymmddhh( const int yyyymmddhh ) {
  int result = 0;

  if ( yyyymmddhh <= 2147123123 ) {
    const int yyyy = yyyymmddhh / 1000000;

    if ( IN_RANGE( yyyy, 1900, 2147 ) ) {
      const int mm = yyyymmddhh / 10000 % 100;

      if ( IN_RANGE( mm , 1, 12 ) ) {
        const int dd = yyyymmddhh / 100 % 100;
        const int month_days = days_in_month( yyyy, mm );

        if ( IN_RANGE( dd, 1, month_days ) ) {
          const int hh = yyyymmddhh % 100;
          result = IN_RANGE( hh , 0, 23 );
        }
      }
    }
  }

  return result;
}



/* Get number of days in yyyy-mm: */

static int days_in_month( const int year, const int month ) {
  static const int days_per_month[ 2 ][ 12 ] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
  };
  const int leap =
    month != 2 ? 0 : year % 4 == 0 && ! ( year % 100 == 0 && year % 400 != 0 );
  const int result =
    IN_RANGE( month, 1, 12 ) ? days_per_month[ leap ][ month - 1 ] : 0;
  return result;
}



/* Read all input data into memory: */

static int read_all_data( Data* const data ) {
  int result = 0;
  int ok = 1;
  int timestep = 0;
  int timestep_offset = 0;
  const char* input_file_name = 0;
  char* temp = 0;

  for ( input_file_name = strtok_r( data->input_file_names, "\n", &temp );
        input_file_name && ok && timestep < data->timesteps;
        input_file_name = strtok_r( 0, "\n", &temp ),
         ++timestep, timestep_offset += data->rows * data->columns ) {
    data->yyyymmddhh[ timestep ] = parse_timestamp( input_file_name );

    if ( data->yyyymmddhh[ timestep ] ) {
      int input_file = -1;
      int status = nc_open( input_file_name, NC_NOWRITE, &input_file );

      if ( status != NC_NOERR ) {
        fprintf( stderr, "\nFailed to open file '%s' for reading.\n",
                input_file_name );
        fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
        ok = 0;
      } else {
        int rows = 0;
        int columns = 0;
        ok = read_file_dimensions( input_file, &rows, &columns );

        if ( ok ) {

          if ( data->rows == 0 ) {
            const size_t bytes = data->timesteps * rows * columns;
            data->rows = rows;
            data->columns = columns;
            data->data = malloc( bytes );

            if ( ! data->data ) {
              fprintf( stderr, "\nFailed to allocate %lu bytes for data.\n",
                       bytes );
              ok = 0;
            } else if ( rows != data->rows || columns != data->columns ) {
              fprintf( stderr, "\nUnmatched dimensions (%d x %d) in file %s.\n",
                       rows, columns, input_file_name );
              ok = 0;
            } else {
              memset( data->data, 0, bytes );
              ok = read_file_data( input_file, rows, columns,
                                   data->minimum, data->maximum,
                                   data->data + timestep_offset );
              swap_data_rows( data->data + timestep_offset, rows, columns );
              ++result;
            }
          }
        }

        nc_close( input_file ), input_file = -1;
      }
    }
  }

  return result;
}



/* Get/check file dimensions: row and column: */

static int read_file_dimensions( const int file_id,
                                 int* const rows,
                                 int* const columns ) {
  int result = 0;
  int dimensions = 0;
  int variables = 0;
  int attributes = 0;
  int unlimited = 0;
  int status =
    nc_inq( file_id, &dimensions, &variables, &attributes, &unlimited );
  *rows = *columns = 0;

  if ( status != NC_NOERR ) {
    fprintf( stderr, "\nFailed to read dimensions of file.\n" );
    fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
  } else if ( dimensions != 2  ) {
    fprintf( stderr, "\nInvalid dimensions (%d) of file.\n",
             dimensions );
  } else if ( variables != 1 ) {
    fprintf( stderr, "\nInvalid variables (%d) of file.\n",
             variables );
  } else {
    size_t dimension0 = 0;
    status = nc_inq_dimlen(file_id, 0, &dimension0 );

    if ( status != NC_NOERR ) {
      fprintf( stderr, "\nFailed to read first dimension of file.\n" );
      fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
    } else if ( dimension0 > INT_MAX ) {
      fprintf( stderr, "\nFile dimension %lu is too large.\n", dimension0 );
    } else {
      size_t dimension1 = 0;
      status = nc_inq_dimlen(file_id, 1, &dimension1 );

      if ( status != NC_NOERR ) {
        fprintf( stderr, "\nFailed to read second dimension of file.\n" );
        fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
      } else if ( dimension1 > INT_MAX / dimension0 ) {
        fprintf( stderr, "\nFile dimension %lu is too large.\n", dimension1 );
      } else {
        *rows = (int) dimension0;
        *columns = (int) dimension1;
        result = 1;
      }
    }
  }

  return result;
}



/* Read file data and return number of values in range [minimum, maximum]: */

static int read_file_data( const int file_id,
                           const int rows,
                           const int columns,
                           const double minimum,
                           const double maximum,
                           float* const data ) {
  int result = 0;
  int status = 0;
  const size_t starts[ 2 ] = { 0, 0 };
  size_t counts[ 2 ] = { 0, 0 };
  counts[ 0 ] = rows;
  counts[ 1 ] = columns;
  status = nc_get_vara_float( file_id, 0, starts, counts, data );

  if ( status != NC_NOERR ) {
    fprintf( stderr, "\nFailed to read file data.\n" );
    fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
  } else {
    const int count = rows * columns;
    int index = 0;

    for ( index = 0; index < count; ++index ) {
      const float value = data[ index ];

      if ( ! IN_RANGE( value, minimum, maximum ) ) {
        data[ index ] = -9999.0;
      } else {
        ++result;
      }
    }
  }

  return result;
}



/* Write bin-format ASCII header and binary data to stdout: */

static int write_output( Data* const data ) {
  int result = 0;
  const int timesteps = data->timesteps;
  const size_t count = timesteps * data->rows * data->columns;
  int timestep = 0;

  /* Write ASCII header: */

  printf( "Content-type: application/octet-stream; charset=iso-8859-1\n" );
  printf( "# variable units:\n" );
  printf( "%s %s\n", data->name, data->units );
  printf( "# dimensions: "
          "timesteps rows columns lonmin lonmax latmin latmax\n" );
  printf( "%-5d %10d %10d %24.18f %24.18f %24.18f %24.18f\n",
          timesteps, data->rows, data->columns,
          data->west, data->east, data->south, data->north );
  printf( "# char yyyymmddhh[timesteps][11] and\n" );
  printf( "# IEEE-754 32-bit float data[timesteps][rows][columns]:\n" );

  /* Write timestamps: */

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    printf( "%d\n", data->yyyymmddhh[ timestep ] );
  }

  /* Write binary MSB 32-bit IEEE-754 data: */

  rotate_4_byte_words_if_little_endian( count, data->data );
  result = fwrite( data->data, sizeof data->data[0], count, stdout ) == count;

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write all %lu bytes of data.\n",
             count * sizeof (float) );
    perror( "because" );
  }

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



/*
 * swap_data_rows - Swap data rows in-place.
 * float array[ timesteps ][ rows ][ columns ]
 */

static void swap_data_rows( float array[], const int rows, const int columns) {
  int lower_row_index = 0;
  int upper_row_index = rows - 1;

  for ( ; lower_row_index < upper_row_index;
       ++lower_row_index, --upper_row_index ) {
    float* const lower_row_data = array + lower_row_index * columns;
    float* const upper_row_data = array + upper_row_index * columns;
    int column = 0;

    for ( column = 0; column < columns; ++column ) {
      float temp = lower_row_data[ column ];
      lower_row_data[ column ] = upper_row_data[ column ];
      upper_row_data[ column ] = temp;
    }
  }
}




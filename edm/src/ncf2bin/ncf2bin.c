/*
ncf2bin.c - C program to convert NetCDF grid data to bin format.
2024-02-22 plessel.todd@epa.gov
To compile:
 $ gcc -Wall -O -I. -o ncf2bin ncf2bin.c libnetcdf.a ; strip ncf2bin
Usage:
ncf2bin list_file name units minimum maximum west south east north [depth] [yyyymm]" \
 > output.bin
data values outside the range [minimum, maximum] are mapped to -9999.
Usage example: Get April mean water temperature near Florida:

ln -s /data/land_use/woa_decav_t04_04.nc /data/tmp/temp.$$_202304.nc
echo /data/tmp/temp.$$_202304.nc >  /data/tmp/temp_file_list.$$

ncf2bin /data/tmp/temp_file_list.$$ temperature C 0 50 -85 28 -84 29 100 > temperature.bin

/bin/rm /data/tmp/temp.$$_202304.nc /data/tmp/temp_file_list.$$

head -8 junk ; tail +9 junk | fdd conv=real4-ascii conv=swab4
Content-type: application/octet-stream; charset=iso-8859-1
# variable units:
temperature C
# dimensions: timesteps rows columns lonmin lonmax latmin latmax
1              5          5   -85.250   -84.000    27.750    29.000
# char yyyymmddhh[timesteps][11] and
# IEEE-754 32-bit float data[timesteps][rows][columns]:
2023040100
1.9742109298706055e+01
1.9575889587402344e+01
1.9179010391235352e+01
-9.9990000000000000e+03
-9.9990000000000000e+03
1.9492290496826172e+01
1.9620990753173828e+01
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
1.9073808670043945e+01
1.9300899505615234e+01
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
1.8354900360107422e+01
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03
-9.9990000000000000e+03

*/

#include <stdio.h>     /* For stderr, fprintf(), fwrite(). */
#include <stdlib.h>    /* For strtod(), malloc(), free(). */
#include <string.h>    /* For strcmp(), memset(), strstr(), strchr(). */
#include <ctype.h>     /* For isalpha(). */
#include <limits.h>    /* For INT_MAX. */
#include <assert.h>    /* For macro assert(). */
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

/* Map specified command-line variable name to woa*.nc file variable: */

static const char* const woa_file_variable_names[] = {
  "nitrate=n_an",
  "oxygen=o_an",
  "oxygen_saturation=O_an",
  "phosphate=p_an",
  "salinity=s_an",
  "silicate=i_an",
  "temperature=t_an"
};

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
  double west_edge;
  double south_edge;
  double east_edge;
  double north_edge;
  double depth;
  int yyyymm;
  int subset_rows;
  int subset_columns;
  int timesteps;
  int rows;
  int columns;
  int layers;
  int* yyyymmddhh; /* yyyymmddhh[ data->timesteps ]. */
  float* values;   /* values[ data->timesteps * data->rows * data->columns ] */
  float* depths;   /* depths[ data->layers ] */
} Data;

static void deallocate_data( Data* const data ) {

  if ( data ) {

    if ( data->input_file_names ) {
      free( data->input_file_names ), data->input_file_names = 0;
    }

    if ( data->yyyymmddhh ) {
      free( data->yyyymmddhh ), data->yyyymmddhh = 0;
    }

    if ( data->values ) {
      free( data->values ), data->values = 0;
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
                                 int* const columns,
                                 int* layers );

static int read_file_data( const int file_id,
                           const int timestep_offset,
                           Data* const data );

static int lookup_variable_id( const int file_id, const char* const name );

static int compute_subset( const int file_id, Data* const data,
                           size_t starts[], size_t counts[] );

static size_t depth_index( const int file_id, const double depth,
                           const int layers, float depths[] );

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
          "west south east north [depth] > output_file\n", program );
  fprintf( stderr, "example: %s list_file salinity PSU 0 50 "
           "-90.5 30.5 -89.5 31.5 > salinity_daily_20120202.bin\n", program );
  fprintf( stderr, "head -7 salinity_daily_20120202.bin\n\n" );
}



/* Read and check command-line arguments: */

static int parse_arguments( int argc, char* argv[], Data* const data ) {
  int result = 0;
  memset( data, 0, sizeof *data );

  if ( argc >= 10 && argc <= 12) {
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

          if ( data->units[ 0 ] == '%' || isalpha( data->units[ 0 ] ) ) {
            char* end = 0;
            data->minimum = strtod( argv[ 4 ], &end );

            if ( end != argv[ 4 ] ) {
              data->maximum = strtod( argv[ 5 ], &end );

              if ( end != argv[ 5 ] && data->maximum >= data->minimum ) {
                data->west = strtod( argv[ 6 ], &end );

                if ( end != argv[ 6 ] && IN_RANGE( data->west, -180.0, 180.0)) {
                  data->south = strtod( argv[ 7 ], &end );

                  if ( end != argv[ 7 ] && IN_RANGE(data->south, -90.0, 90.0)) {
                    data->east = strtod( argv[ 8 ], &end );

                    if ( end != argv[ 8 ] &&
                         IN_RANGE( data->east, data->west, 180.0 ) ) {
                      data->north = strtod( argv[ 9 ], &end );

                      if ( end != argv[ 9 ] &&
                           IN_RANGE( data->north, data->south, 90.0 ) ) {

                        if ( argc >= 11 ) {
                          data->depth = strtod( argv[ 10 ], &end );

                          if ( end != argv[ 10 ] &&
                               IN_RANGE( data->depth, 0.0, 11000.0 ) ) {

                            if ( argc == 12 ) {
                              data->yyyymm = atoi( argv[ 11 ] );

                              if ( data->yyyymm >= 1900 ) {

                                if ( data->yyyymm <= 2147 ) {
                                  const int yyyy = data->yyyymm;
                                  result =
                                    is_yyyymmddhh( yyyy * 1000000 +
                                                   10100 );
                                } else if ( data->yyyymm <= 214712 ) {
                                  result =
                                    is_yyyymmddhh( data->yyyymm * 10000 + 100 );
                                }
                              }
                            } else {
                              result = 1;
                            }
                          }
                        } else {
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
      } else if ( result < 214712 ) {
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
  const int is_monthly = IN_RANGE( data->yyyymm, 190001, 214712 );
  int yyyymmddhh =
    data->yyyymm == 0 ? 0
    : is_monthly ? data->yyyymm * 10000 + 100
    : data->yyyymm * 1000000 + 10100;

  for ( input_file_name = strtok_r( data->input_file_names, "\n", &temp );
        input_file_name && ok && timestep < data->timesteps;
        input_file_name = strtok_r( 0, "\n", &temp ),
         ++timestep,
         timestep_offset += data->subset_rows * data->subset_columns ) {

    if ( yyyymmddhh ) {
      data->yyyymmddhh[ timestep ] = yyyymmddhh;
    } else {
      data->yyyymmddhh[ timestep ] = parse_timestamp( input_file_name );
    }

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
        int layers = 0;
        ok = read_file_dimensions( input_file, &rows, &columns, &layers );

        if ( ok ) {

          if ( data->rows == 0 ) {
            const size_t bytes =
              ( data->timesteps * rows * columns + layers * 2 ) * sizeof (float);
            data->rows = rows;
            data->columns = columns;
            data->layers = layers;
            data->values = malloc( bytes );

            if ( ! data->values ) {
              fprintf( stderr, "\nFailed to allocate %lu bytes for data.\n",
                       bytes );
              ok = 0;
            } else {
              memset( data->values, 0, bytes );
              data->depths = layers ?
                data->values + data->timesteps * rows * columns : 0;
            }
          }

          if ( ok ) {

            if ( rows != data->rows || columns != data->columns ||
                 layers != data->layers ) {
              fprintf( stderr,
                       "Unmatched dimensions (%d x %d x %d) in file %s.\n",
                       rows, columns, layers, input_file_name );
              ok = 0;
            } else {
              ok = read_file_data( input_file, timestep_offset, data );


              if ( ! layers ) {
                swap_data_rows( data->values + timestep_offset,
                                rows, columns );
              }

              ++result;
            }
          }
        }

        nc_close( input_file ), input_file = -1;
      }
    }

    if ( yyyymmddhh ) {

      if ( is_monthly ) { /* Increment yyyymm. */
        int yyyy = yyyymmddhh / 1000000;
        int mm   = yyyymmddhh / 10000 % 100;

        ++mm;

        if ( mm > 12 ) {
          mm = 1;
          ++yyyy;
        }

        yyyymmddhh = yyyy * 1000000 + mm * 10000 + 100;
      } else {
        yyyymmddhh += 1000000; /* Increment yyyy. */
      }
    }
  } /* Loop on files. */

  return result;
}



/* Get/check file dimensions: row and column and possible layers: */

static int read_file_dimensions( const int file_id,
                                 int* const rows,
                                 int* const columns,
                                 int* const layers ) {
  int result = 0;
  int dimensions = 0;
  int variables = 0;
  int attributes = 0;
  int unlimited = 0;
  int status =
    nc_inq( file_id, &dimensions, &variables, &attributes, &unlimited );
  *rows = *columns = *layers = 0;

  if ( status != NC_NOERR ) {
    fprintf( stderr, "\nFailed to read dimensions of file.\n" );
    fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
  } else {
    const int is_woa_file = dimensions == 5 && variables >= 2;

    if ( ! is_woa_file && dimensions != 2  ) {
      fprintf( stderr, "\nInvalid dimensions (%d) of file.\n",
               dimensions );
    } else if ( ! is_woa_file && variables != 1 ) {
      fprintf( stderr, "\nInvalid variables (%d) of file.\n",
               variables );
    } else {
      size_t dimension0 = 0;
      status = nc_inq_dimlen(file_id, is_woa_file, &dimension0 );

      if ( status != NC_NOERR ) {
        fprintf( stderr, "\nFailed to read first dimension of file.\n" );
        fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
      } else if ( dimension0 > INT_MAX ) {
        fprintf( stderr, "\nFile dimension %lu is too large.\n", dimension0 );
      } else {
        size_t dimension1 = 0;
        status = nc_inq_dimlen(file_id, is_woa_file + 1, &dimension1 );

        if ( status != NC_NOERR ) {
          fprintf( stderr, "\nFailed to read second dimension of file.\n" );
          fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
        } else if ( dimension1 > INT_MAX / dimension0 ) {
          fprintf( stderr, "\nFile dimension %lu is too large.\n", dimension1 );
        } else {
          size_t dimension2 = 0;
          *rows = (int) dimension0;
          *columns = (int) dimension1;

          if ( is_woa_file ) {
            status = nc_inq_dimlen(file_id, is_woa_file + 2, &dimension2 );

            if ( status != NC_NOERR ) {
              fprintf( stderr, "\nFailed to read second dimension of file.\n" );
              fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
            } else {
              *layers = (int) dimension2;
              result = 1;
            }

          } else {
            result = 1;
          }
        }
      }
    }
  }

  return result;
}



/* Read file data and return number of values in range [minimum, maximum]: */

static int read_file_data( const int file_id,
                           const int timestep_offset,
                           Data* const data ) {

  int result = 0;
  float* const values = data->values + timestep_offset;
  const int variable_id =
    data->layers ? lookup_variable_id( file_id, data->name ) : 0;

  if ( variable_id >= 0 ) {
    size_t starts[ 4 ] = { 0, 0, 0, 0 };
    size_t counts[ 4 ] = { 0, 0, 0, 0 };
    int count = 0;

    if ( ! data->layers ) {
      counts[ 0 ] = data->rows;
      counts[ 1 ] = data->columns;
      data->subset_rows = data->rows;
      data->subset_columns = data->columns;
    } else {
      compute_subset( file_id, data, starts, counts );
    }

    count = data->subset_rows * data->subset_columns;

    if ( count ) {
      const int status =
        nc_get_vara_float( file_id, variable_id, starts, counts, values );

      if ( status != NC_NOERR ) {
        fprintf( stderr, "\nFailed to read file data.\n" );
        fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
      } else { /* Filter invalid values: */
        const double minimum = data->minimum;
        const double maximum = data->maximum;
        int index = 0;

        for ( index = 0; index < count; ++index ) {
          const float value = values[ index ];

          if ( ! IN_RANGE( value, minimum, maximum ) ) {
            values[ index ] = -9999.0;
          } else {
            ++result;
          }
        }
      }
    }
  }

  return result;
}



/* Get NetCDF variable id of variable corresponding to WCS name: */

static int lookup_variable_id( const int file_id, const char* const name ) {
  static const size_t count =
    sizeof woa_file_variable_names / sizeof *woa_file_variable_names;
  size_t index = 0;
  int result = -1;

  while ( index < count &&
          ! strstr( woa_file_variable_names[ index ], name ) ) {
    ++index;
  }

  if ( index < count ) {
    int status = 0;
    const char* file_variable_name =
      strrchr( woa_file_variable_names[ index ], '=' );
    assert( file_variable_name );
    ++file_variable_name;
    assert( *file_variable_name );
    status = nc_inq_varid( file_id, file_variable_name, &result );

    if ( status != NC_NOERR ) {
      fprintf( stderr, "\nFailed to get file variable %s id.\n",
               file_variable_name );
      fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
      result = -1;
    }
  } else {
    fprintf( stderr, "\nInvalid variable %s.\n", name );
  }

  return result;
}



/*
 * Compute subset indices of global grid within bounds and
 * output:
 *   data->subset_rows
 *   data->subset_columns
 *   data->west_edge
 *   data->east_edge
 *   data->south_edge
 *   data->north_edge
 * and return number of 2D grid subset points.
 */

static int compute_subset( const int file_id, Data* const data,
                           size_t starts[], size_t counts[] ) {
  const double west = data->west;
  const double east = data->east;
  const double south = data->south;
  const double north = data->north;
  const double depth = data->depth;
  const double cell_width  = 360.0 / data->columns;
  const double cell_height = 180.0 / data->rows;
  double edge = -90.0 + cell_height;
  int index = 0;
  int result = 0;

  /* Data array is float t_an(time, depth, lat, lon) */

  /* Time: */

  starts[ 0 ] = 0;
  counts[ 0 ] = 1;

  /* Depth index: */

  starts[ 1 ] = 0;
  counts[ 1 ] = 1;

  if ( depth > 0.0 && data->layers > 0 &&
       data->layers * 2 <= data->rows * data->columns ) {
    starts[ 1 ] =
      depth_index( file_id, depth, data->layers, data->depths );
  }

  /* Compute subset rows: */

  for ( index = 0; edge < south; ++index ) {
    edge += cell_height;
  }

  starts[ 2 ] = index;
  data->south_edge = edge - cell_height;

  for ( ; edge < north; ++index ) {
    edge += cell_height;
  }

  counts[ 2 ] = 1 + index - starts[ 2 ];
  data->subset_rows = counts[ 2 ];
  data->north_edge = edge;

  /* Compute subset columns: */

  edge = -180.0 + cell_width;

  for ( index = 0; edge < west; ++index ) {
    edge += cell_width;
  }

  starts[ 3 ] = index;
  data->west_edge = edge - cell_width;

  for ( ; edge < east; ++index ) {
    edge += cell_width;
  }

  counts[ 3 ] = 1 + index - starts[ 3 ];
  data->subset_columns = counts[ 3 ];
  data->east_edge = edge;

  result = counts[ 2 ] * counts[ 3 ];
  return result;
}



/* Compute layer index of depth: */

static size_t depth_index( const int file_id, const double depth,
                           const int layers, float depths[] ) {
  size_t result = 0;
  int id = -1;
  int status = nc_inq_varid( file_id, "depth_bnds", &id );

  if ( status == NC_NOERR ) {
    const size_t starts[ 2 ] = { 0, 0 };
    size_t counts[ 2 ] = { 0, 0 };
    counts[ 0 ] = layers;
    counts[ 1 ] = 2;
    status = nc_get_vara_float( file_id, id, starts, counts, depths );

    if ( status == NC_NOERR ) {

      while ( result < layers && depth > depths[ result + result + 1 ] ) {
        ++result;
      }

      if ( result == layers ) {
        result = layers - 1;
      }
    }
  }

  if ( status != NC_NOERR ) {
    fprintf( stderr, "\nFailed to read file depth data.\n" );
    fprintf( stderr, "Because: %s\n", nc_strerror( status ) );
  }

  return result;
}



/* Write bin-format ASCII header and binary data to stdout: */

static int write_output( Data* const data ) {
  int result = 0;
  const int timesteps = data->timesteps;
  const size_t count = timesteps * data->subset_rows * data->subset_columns;
  int timestep = 0;

  /* Write ASCII header: */

  printf( "Content-type: application/octet-stream; charset=iso-8859-1\n" );
  printf( "# variable units:\n" );
  printf( "%s %s\n", data->name, data->units );
  printf( "# dimensions: "
          "timesteps rows columns lonmin lonmax latmin latmax\n" );
  printf( "%-5d %10d %10d %24.18f %24.18f %24.18f %24.18f\n",
          timesteps, data->subset_rows, data->subset_columns,
          data->west_edge, data->east_edge,
          data->south_edge, data->north_edge );
  printf( "# char yyyymmddhh[timesteps][11] and\n" );
  printf( "# IEEE-754 32-bit float data[timesteps][rows][columns]:\n" );

  /* Write timestamps: */

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    printf( "%d\n", data->yyyymmddhh[ timestep ] );
  }

  /* Write binary MSB 32-bit IEEE-754 data: */

  rotate_4_byte_words_if_little_endian( count, data->values );
  result =
    fwrite( data->values, sizeof data->values[0], count, stdout ) == count;

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




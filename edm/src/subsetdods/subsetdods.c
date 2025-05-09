/*
subsetdods.c - C program to read a sequence of DODS gridded 1-degree global
data, subset it to a given lonlat box and output it in bin format.
2013-04-29 plessel.todd@epa.gov
To compile:
$ cc -o subsetdods subsetdods.c
Usage:
subsetdods variable units minimum maximum lonmin latmin lonmax latmax \
  startdate hours_per_timestep timesteps input_files > output.bin
data outside the range [minimum, maximum] is mapped to -9999.
Usage examples:

$ curl -k -g --silent --retry 0 -L --tcp-nodelay --max-time 0 \
 'https://opendap.jpl.nasa.gov/opendap/hyrax/allData/\
 aquarius/L3/mapped/V5/daily/SCI/2013/088/\
 Q2013088.L3m_DAY_SCI_V5.0_SSS_1deg.bz2.dods?/SSS[0:1:179][0:1:359]' \
 > daily_salinity_20130329.dods

$ curl -k -g --silent --retry 0 -L --tcp-nodelay --max-time 0 \
 'https://opendap.jpl.nasa.gov/opendap/hyrax/allData/\
 aquarius/L3/mapped/V5/daily/SCI/2013/089/\
 Q2013089.L3m_DAY_SCI_V5.0_SSS_1deg.bz2.dods?/SSS[0:1:179][0:1:359]' \
 > daily_salinity_20130330.dods

$ curl -k -g --silent --retry 0 -L --tcp-nodelay --max-time 0 \
 'https://opendap.jpl.nasa.gov/opendap/hyrax/allData/\
 aquarius/L3/mapped/V5/daily/SCI/2013/090/\
 Q2013090.L3m_DAY_SCI_V5.0_SSS_1deg.bz2.dods?/SSS[0:1:179][0:1:359]' \
 > daily_salinity_20130331.dods

$ head -4 daily_salinity_20130329.dods
Dataset {
    Float32 /SSS[180][360];
} Q2013088.L3m_DAY_SCI_V5.0_SSS_1deg;
Data:

$ ls -1 test/daily_salinity_*.dods > test/input_files

$ subsetdods salinity PSU 0 50 -78 35 -70 40 20130329 24 3 test/input_files > \
  test/daily_salinity.bin

$ head -7 daily_salinity.bin
Content-type: application/octet-stream; charset=iso-8859-1
# variable units:
sst C
# dimensions: timesteps z rows columns lonmin lonmax latmin latmax
3       0        25         25  \
 -84.02304222733599204  -81.98243622733599523   \
 32.97649710028599657   34.04116110028599707
# char yyyymmddhh[timesteps][11] and
# IEEE-754 32-bit float data[timesteps][rows][columns]:

Note: DODS input arrays are MSB/big-endian IEEE-754 32-bit floating-point
format and are preceeded by two 4-byte words (MSB integer) containing the
count of array items that follows.
*/


#include <stdio.h>  /* For FILE, stdoerr, scanf(), printf(), fwrite(). */
#include <stdlib.h> /* For atof(), malloc(), free(). */
#include <string.h> /* For strcmp(), memset(). */
#include <math.h>   /* For floor(), ceil(). */


#if defined(__alpha) || \
    defined(__i386__) || defined(__i486__) || \
    defined(__i586__) || defined(__i686__) || \
    defined(__ia64__) || defined(__x86_64__)
#define IS_LITTLE_ENDIAN 1
#else
#define IS_LITTLE_ENDIAN 0
#endif

#define IN_RANGE( x, lower, upper ) ( (x) >= (lower) && (x) <= (upper) )

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

#define MISSING (-9999.0)

/* Forward declarations: */

static void usage( const char* program );

static int parse_options( const int argc, const char* argv[],
                          const char** variable, const char** units,
                          double* minimum, double* maximum,
                          double* longitude_minimum, double* latitude_minimum,
                          double* longitude_maximum, double* latitude_maximum,
                          int* yyyymmdd, int* hours, int* timesteps,
                          const char** files );

static int read_header( FILE* file );

static int read_data( FILE* file );

static int copy_subset_data( const int first_row, const int last_row,
                             const int first_column, const int last_column,
                             const double minimum, const double maximum,
                             float subset_data[] );

static int write_output( const int timesteps,
                         const int yyyymmdd, const int hours,
                         const int rows, const int columns,
                         const char* const name, const char* const units,
                         const double longitude_minimum,
                         const double longitude_maximum,
                         const double latitude_minimum,
                         const double latitude_maximum,
                         float subset_data[] );

static int skip_8_bytes( FILE* file );

static void reverse_row_order( const size_t timesteps,
                               const size_t rows,
                               const size_t columns,
                               float array[] );

static void reverse_4_byte_words_if_little_endian( const size_t count,
                                                  void* array );

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours );

static int days_in_month( const int year, const int month );

static int clamped_to_range( const int value,
                             const int minimum, const int maximum );

static float data[ 180 ][ 360 ]; /* data[ rows ][ columns ]. */



/* Read input DODS data from stdin, write output bin data to stdout: */

int main( int argc, char* argv[] ) {
  int ok = 0;
  const char* variable = 0;
  const char* units = 0;
  double minimum = 0.0;
  double maximum = 0.0;
  double longitude_minimum = 0.0;
  double latitude_minimum = 0.0;
  double longitude_maximum = 0.0;
  double latitude_maximum = 0.0;
  int yyyymmdd = 0;
  int hours = 0;
  int timesteps = 0;
  const char* files = 0;

  if ( argc != 13 ) {
    usage( argv[ 0 ] );
  } else if ( parse_options( argc, (const char**) argv, &variable, &units,
                             &minimum, &maximum,
                             &longitude_minimum, &latitude_minimum,
                             &longitude_maximum, &latitude_maximum,
                             &yyyymmdd, &hours, &timesteps,
                             &files ) ) {

    /* Expand lon-lat bounds to be whole degrees since the data grid is: */

    longitude_minimum =
      clamped_to_range( (int) floor( longitude_minimum ), -180, 179 );
    longitude_maximum =
      clamped_to_range( (int) ceil( longitude_maximum ), -179, 180 );

    latitude_minimum =
      clamped_to_range( (int) floor( latitude_minimum ), -90, 89 );
    latitude_maximum =
      clamped_to_range( (int) ceil( latitude_maximum ), -89, 90 );

    /* Ensure there is at least one grid cell width, height: */

    if ( longitude_minimum == longitude_maximum ) {

      if ( longitude_minimum < -179.0 ) {
        ++longitude_maximum;
      } else {
        --longitude_minimum;
      }
    }

    if ( latitude_minimum == latitude_maximum ) {

      if ( latitude_minimum < -89.0 ) {
        ++latitude_maximum;
      } else {
        --latitude_minimum;
      }
    }

    /*
     * Compute row and column indices.
     * Note rows are north to south and columns are west to east:
     */

    {
      const int first_row = (int) ( 90.0 - latitude_maximum );
      const int latitude_difference =
        (int) ( latitude_maximum - latitude_minimum );
      const int rows = clamped_to_range( latitude_difference, 1, 180 );
      const int last_row = first_row + rows - 1;

      const int first_column = (int) ( 180.0 + longitude_minimum );
      const int longitude_difference =
        (int) ( longitude_maximum - longitude_minimum );
      const int columns =
        clamped_to_range( longitude_difference, 1, 360 );
      const int last_column = first_column + columns - 1;

      const size_t rows_times_columns = rows * columns;
      const size_t subset_bytes =
        timesteps * rows_times_columns * sizeof (float);
      float* subset_data = malloc( subset_bytes );

      DEBUG( fprintf( stderr, "%lu bytes: %d x %d = "
                      "row[%d %d] col[%d %d] lat(%lg %lg) lon(%lg %lg)\n",
                      subset_bytes, rows, columns,
                      first_row, last_row, first_column, last_column,
                      latitude_minimum, latitude_maximum,
                      longitude_minimum, longitude_maximum ); )

      if ( ! subset_data ) {
        fprintf( stderr, "\a\n\nI'm sorry: "
                 "Failed to allocate %lu bytes "
                 "to complete the requested operation.\n\n", subset_bytes );
      } else {
        FILE* inputs = fopen( files, "r" );

        if ( inputs ) {
          int timestep = 0;
          char file_name[ 256 ] = "";
          memset( file_name, 0, sizeof file_name );
          memset( subset_data, 0, subset_bytes );

          for ( timestep = 0; timestep < timesteps; ++timestep ) {

            if ( fgets( file_name,
                        sizeof file_name / sizeof *file_name, inputs ) ) {
              char* const newline = strrchr( file_name, '\n' );

              if ( newline ) {
                *newline = '\0';
              }

              {
                FILE* input = fopen( file_name, "rb" );
                DEBUG( fprintf( stderr, "file_name = '%s', input = %p\n",
                                file_name, input ); )

                if ( input ) {
                  const size_t offset = timestep * rows_times_columns;
                  int read_valid_subset_data = 0;

                  if ( read_header( input ) ) {

                    if ( read_data( input ) ) {
                      read_valid_subset_data =
                        copy_subset_data( first_row, last_row,
                                          first_column, last_column,
                                          minimum, maximum,
                                          subset_data + offset );
                    }
                  }

                  DEBUG( fprintf( stderr, "read_valid_subset_data = %d\n",
                                  read_valid_subset_data ); )

                  if ( ! read_valid_subset_data ) {

                    /* Fill 'dropped-out' subset with MISSING: */

                    size_t count = rows_times_columns;
                    float* output = subset_data + offset;

                    while ( count-- ) {
                      *output++ = MISSING;
                    }
                  } else {
                    ok = 1;
                  }

                  fclose( input ), input = 0;
                }
              }
            }
          }

          fclose( inputs ), inputs = 0;
        }

        if ( ok ) {
          write_output( timesteps, yyyymmdd, hours, rows, columns,
                        variable, units,
                        longitude_minimum, longitude_maximum,
                        latitude_minimum, latitude_maximum,
                        subset_data );
        }

        free( subset_data ), subset_data = 0;
      }
    }
  }

  return ! ok;
}



/* Print program usage instructions to stderr: */

static void usage( const char* program ) {
  fprintf( stderr, "\n%s - Read a sequence of DODS gridded 1-degree global "
                   "data,\nsubset it to a given lonlat box and output it "
                   "in bin format.\n", program );
  fprintf( stderr, "usage: %s variable units minimum maximum "
                   "lonmin latmin lonmax latmax "
                   "startdate hours_per_timestep timesteps input_files "
                   "> output.bin\n", program );
  fprintf( stderr, "example: %s salinity PSU 0 50 -78 35 -70 40 20130329 24 3 "
                   "input_files > salinity.bin\n", program );
  fprintf( stderr, "head -7 salinity.bin\n\n" );
}



/* Read/check command-line options: */

static int parse_options( const int argc, const char* argv[],
                          const char** variable, const char** units,
                          double* minimum, double* maximum,
                          double* longitude_minimum, double* latitude_minimum,
                          double* longitude_maximum, double* latitude_maximum,
                          int* yyyymmdd, int* hours, int* timesteps,
                          const char** files ) {

  int result = 0;
  *variable = *units = 0;
  *minimum = *maximum = 0.0;
  *longitude_minimum = *longitude_maximum = 0.0;
  *latitude_minimum = *latitude_maximum = 0.0;
  *yyyymmdd = *hours = *timesteps = 0;
  *files = 0;

  if ( argc == 13 ) {
    *variable = argv[ 1 ];

    if ( *variable ) {
      *units = argv[ 2 ];

      if ( *units ) {
        *minimum = atof( argv[ 3 ] );
        *maximum = atof( argv[ 4 ] );

        if ( *maximum > *minimum ) {
          *longitude_minimum = atof( argv[ 5 ] );

          if ( IN_RANGE( *longitude_minimum, -180.0, 180.0 ) ) {
            *latitude_minimum = atof( argv[ 6 ] );

            if ( IN_RANGE( *latitude_minimum, -90.0, 90.0 ) ) {
              *longitude_maximum = atof( argv[ 7 ] );

              if ( IN_RANGE( *longitude_maximum, *longitude_minimum, 180.0 )) {
                *latitude_maximum = atof( argv[ 8 ] );

                if ( IN_RANGE( *latitude_maximum, *latitude_minimum, 90.0 ) ) {
                  *yyyymmdd = atoi( argv[ 9 ] );

                  {
                    const int yyyy = *yyyymmdd / 10000;
                    const int mm = *yyyymmdd / 100 % 100;
                    const int dd = *yyyymmdd % 100;

                    if ( IN_RANGE( yyyy, 1900, 3000 )&&
                         IN_RANGE( mm, 1, 12 ) &&
                         IN_RANGE( dd, 1, days_in_month( yyyy, mm ) ) ) {
                      *hours = atoi( argv[ 10 ] );

                      if ( *hours >= 24 ) {
                        *timesteps = atoi( argv[ 11 ] );

                        if ( *timesteps >= 1 ) {
                          *files = argv[ 12 ];
                          result = *files[ 0 ] != '\0';
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

  if ( ! result ) {
    fprintf( stderr, "\nInvalid command-line options.\n" );
  }

  return result;
}



/* Read/check DODS header: */

static int read_header( FILE* file ) {
  int result = 0;
  enum { LENGTH = 255 };
  char line[ LENGTH + 1 ] = "";
  memset( line, 0, sizeof line );

  if ( fgets( line, LENGTH, file ) ) {

    if ( ! strcmp( line, "Dataset {\n" ) ) {

      if ( fgets( line, LENGTH, file ) ) {

        if ( ! strcmp( line, "    Float32 l3m_data[180][360];\n" ) ||
             ! strcmp( line, "    Float32 SSS[180][360];\n" ) ) {

          if ( fgets( line, LENGTH, file ) ) {

            if ( strstr( line, ".L3m_DAY_SCI_V5.0_SSS_1deg;\n" ) ||
                 strstr( line, ".L3m_7D_SCI_V5.0_SSS_1deg;\n" ) ||
                 strstr( line, ".L3m_MO_SCI_V5.0_SSS_1deg;\n" ) ||
                 strstr( line, ".L3m_YR_SCI_V5.0_SSS_1deg;\n" ) ) {

              if ( fgets( line, LENGTH, file ) ) {
                result = ! strcmp( line, "Data:\n" );
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nInvalid input DODS header.\n" );
  }

  return result;
}



/* Read input file data into data[]: */

static int read_data( FILE* file ) {
  int result = skip_8_bytes( file );

  if ( result ) {
    const size_t count = sizeof data / sizeof data[ 0 ][ 0 ];
    result = fread( data, sizeof (float), count, file ) == count;
    DEBUG( fprintf( stderr, "Read %lu floats = %d.\n", count, result ); )

    if ( result ) {
      reverse_4_byte_words_if_little_endian( count, data );
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read all DODS file data.\n" );
  }

  return result;
}



/* Copy subset of data[] into subset_data[]: */

static int copy_subset_data( const int first_row, const int last_row,
                             const int first_column, const int last_column,
                             const double minimum, const double maximum,
                             float subset_data[] ) {
  int result = 0;
  int row = 0;
  float* output = subset_data;

  for ( row = first_row; row <= last_row; ++row ) {
    int column = 0;

    for ( column = first_column; column <= last_column; ++column ) {
      const float value = data[ row ][ column ];
      DEBUG( fprintf( stderr, " %g", value ); )

      if ( IN_RANGE( value, minimum, maximum ) ) {
        *output++ = value;
        result = 1;
      } else {
        *output++ = MISSING;
      }
    }
  }

  return result;
}



/* Write bin-format ASCII header and binary data to stdout: */

static int write_output( const int timesteps,
                         const int yyyymmdd, const int hours,
                         const int rows, const int columns,
                         const char* const name, const char* const units,
                         const double longitude_minimum,
                         const double longitude_maximum,
                         const double latitude_minimum,
                         double latitude_maximum,
                         float subset_data[] ) {
  const size_t count = timesteps * rows * columns;
  int result = 0;
  int timestep = 0;
  int yyyy = yyyymmdd / 10000;
  int mm = yyyymmdd / 100 % 100;
  int yyyymmddhh = yyyymmdd * 100;
  const double z = 0.0;

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

    if ( hours == 8760 ) { /* Yearly: */
      yyyymmddhh = yyyy * 1000000 + 10100; /* January 1. */
      ++yyyy;
    } else if ( hours == 744 ) { /* Monthly: */
      yyyymmddhh = yyyy * 1000000 + mm * 10000 + 100; /* First day of month. */
      ++mm;

      if ( mm > 12 ) {
        mm = 1;
        ++yyyy;
      }
    }

    printf( "%010d\n", yyyymmddhh );

    if ( hours != 8760 && hours != 744 ) { /* Weekly or daily: */
      increment_yyyymmddhh( &yyyymmddhh, hours );
    }
  }

  /* Write binary MSB 32-bit IEEE-754 data: */

  reverse_row_order( timesteps, rows, columns, subset_data ); /* south->north*/
  reverse_4_byte_words_if_little_endian( count, subset_data ); /* Big-endian.*/
  result = fwrite( subset_data, sizeof (float), count, stdout ) == count;

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write all %lu bytes of data.\n",
             count * sizeof (float) );
    perror( "because" );
  }

  return result;
}



/* Read and ignore 8 bytes from file: */

static int skip_8_bytes( FILE* file ) {
  double unused = 0.0;
  const int result = fread( &unused, 8, 1, file ) == 1;
  return result;
}



/* Reverse row order from north->south to south->north: */

static void reverse_row_order( const size_t timesteps,
                               const size_t rows,
                               const size_t columns,
                               float array[] ) {

  const size_t rows_times_columns = rows * columns;
  const size_t rows_minus_1_times_columns = ( rows - 1 ) * columns;
  size_t timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    float* north_row = array + timestep * rows_times_columns;
    float* south_row = north_row + rows_minus_1_times_columns;

    for ( ; north_row < south_row; north_row += columns, south_row -= columns ) {
      size_t column = 0;

      for ( column = 0; column < columns; ++column ) {
        const float swap_temp = north_row[ column ];
        north_row[ column ] = south_row[ column ];
        south_row[ column ] = swap_temp;
      }
    }
  }
}




/* On little-endian platforms, change endianess of an array of 4-byte words: */

static void reverse_4_byte_words_if_little_endian( const size_t count,
                                                   void* array ) {

#if IS_LITTLE_ENDIAN

  unsigned int* const array4 = array;
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



/* Number of days in year/month: */

static int days_in_month( const int year, const int month ) {
  static const int days_per_month[ 2 ][ 12 ] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
  };
  const int leap =
    month != 2 ? 0 : year % 4 == 0 && ! ( year % 100 == 0 && year % 400 != 0 );
  const int result = days_per_month[ leap ][ month - 1 ];
  return result;
}



/* Return value clamped to range [minimum, maximum]: */

static int clamped_to_range( const int value,
                             const int minimum, const int maximum ) {
  const int result =
    value < minimum ? minimum : value > maximum ? maximum : value;
  return result;
}




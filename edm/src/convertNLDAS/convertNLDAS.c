/*
convertNLDAS.c - C program to read a sequence of NetCDF gridded NLDAS files
and convert it to a bin format grid file and write it to stdout.
2016-08-03 plessel.todd@epa.gov

To compile:
  cc -I../../../include -L../../../lib/$platform \
     -o convertNLDAS convertNLDAS.c -lNetCDF -lm

Usage:
convertNLDAS variable units minimum maximum \
  yyyymmddhh hours_per_timestep timesteps input_files > output.bin

Data outside the range [minimum, maximum] is mapped to -9999.

Example:

echo 'machine urs.earthdata.nasa.gov login rsig password Rsig4444' > \
test/temp_netrc

touch test/junk ; /bin/rm -f test/junk*

/rsig/current/code/bin/Linux.x86_64/curl -k --silent --retry 0 -L \
  --max-redirs 10 --tcp-nodelay --max-time 3600 \
  --netrc-file nldas_netrc -c nldas_cookie -b nldas_cookie \
 'https://hydro1.gesdisc.eosdis.nasa.gov/thredds/ncss/grid/NLDAS_aggregation/\
 NLDAS_FORA0125_H.2.0/NLDAS_FORA0125_H.2.0_Aggregation_2016.ncml?\
 time_start=2016-07-29T00:00:00Z&time_end=2016-07-29T00:59:59Z\
 &horizStride=1&addLatLon=true&accept=netcdf3\
 &north=40.125&west=-76&east=-74.875&south=39&var=Wind_E,Wind_N' \
 > test/wind00.nc

/rsig/current/code/bin/Linux.x86_64/curl -k --silent --retry 0 -L \
  --max-redirs 10 --tcp-nodelay --max-time 3600 \
  --netrc-file nldas_netrc -c nldas_cookie -b nldas_cookie \
 'https://hydro1.gesdisc.eosdis.nasa.gov/thredds/ncss/grid/NLDAS_aggregation/\
 NLDAS_FORA0125_H.2.0/NLDAS_FORA0125_H.2.0_Aggregation_2016.ncml?\
 time_start=2016-07-29T01:00:00Z&time_end=2016-07-29T01:59:59Z\
 &horizStride=1&addLatLon=true&accept=netcdf3\
 &north=40.125&west=-76&east=-74.875&south=39&var=Wind_E,Wind_N' \
 > test/wind01.nc

/rsig/current/code/bin/Linux.x86_64/curl -k --silent --retry 0 -L \
  --max-redirs 10 --tcp-nodelay --max-time 3600 \
  --netrc-file nldas_netrc -c nldas_cookie -b nldas_cookie \
 'https://hydro1.gesdisc.eosdis.nasa.gov/thredds/ncss/grid/NLDAS_aggregation/\
 NLDAS_FORA0125_H.2.0/NLDAS_FORA0125_H.2.0_Aggregation_2016.ncml?\
 time_start=2016-07-29T02:00:00Z&time_end=2016-07-29T02:59:59Z\
 &horizStride=1&addLatLon=true&accept=netcdf3\
 &north=40.125&west=-76&east=-74.875&south=39&var=Wind_E,Wind_N' \
 > test/wind02.nc

ls -1 test/wind*.nc > test/wind_files

convertNLDAS wind m/s -500 500 2016072900 1 3 test/wind_files > test/wind.bin

$ head -7 wind.bin
Content-type: application/octet-stream; charset=iso-8859-1
# variable units:
wind m/s
# dimensions: timesteps components rows columns lonmin lonmax latmin latmax
2       3       9          9  -76.0  -74.875   39.0   40.125
# char yyyymmddhh[timesteps][11] and
# IEEE-754 32-bit float data[timesteps][components][rows][columns]:
<binary data follows>

*/


#include <stdio.h>  /* For FILE, stderr, scanf(), printf(), fwrite(). */
#include <stdlib.h> /* For atof(), malloc(), free(). */
#include <string.h> /* For strcmp(), memset(). */

#include <netcdf.h> /* For NC* nc*(). */


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

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

#define MISSING (-9999.0)

#define FREE( p ) ( ( (p) ? free(p) : (void) 0 ), (p) = 0 )
#define NEW( type, count ) allocate( sizeof (type) * (count) )

static void* allocate( size_t bytes ) {
  void* result = 0;
  result = malloc( bytes );

  if ( result ) {
    memset( result, 0, bytes );
  } else {
    fprintf( stderr, "\nFailed to allocate %lu bytes "
             "to complete the requested action.\n", bytes );
  }

  return result;
}



/*
 * Translate variable name into NLDAS NetCDF file variable name.
 * See nldasserver for list of parameter / file variable names.
 */

static const char* lookup_file_variable_name( const char* const variable ) {
  typedef struct {
    const char* const variable;
    const char* const file_variable;
  } Entry;
  const Entry table[] = {
    { "precipitation",
      /* "Total_Precipitation_surface_1_Hour_Accumulation" */
      /* "Total_Precipitation_surface_Accumulation" */
      "Rainf"
    },
    { "convective_available_potential_energy",
      /* "Convective_Available_Potential_Energy"
      "_layer_between_two_pressure_difference_from_ground_layer" */
      "CAPE"
    },
    { "convective_precipitation",
      /* "Convective_Precipitation_surface_1_Hour_Accumulation" */
      /* "Convective_Precipitation_surface_Accumulation" */
      "CRainf_frac"
    },
    { "long_wave_radiation_flux",
      /* "Downward_Longwave_Radiation_Flux_surface" */
      "LWdown"
    },
    { "short_wave_radiation_flux",
      /* "Downward_Shortwave_Radiation_Flux_surface" */
      "SWdown"
    },
    { "potential_evaporation",
      /* "Potential_Evaporation_surface_1_Hour_Accumulation" */
      /* "Potential_Evaporation_surface_Accumulation" */
      "PotEvap"
    },
    { "pressure",
      /* "Pressure_surface" */
      "PSurf"
    },
    { "humidity",
      /* "Specific_Humidity_height_above_ground" */
      "Qair"
    },
    { "temperature",
      /* "Temperature_height_above_ground" */
      "Tair"
    },
    { "u_wind",
      /* "u-component_of_wind_height_above_ground" */
      "Wind_E"
    },
    { "v_wind",
      /* "v-component_of_wind_height_above_ground" */
      "Wind_N"
    },
    { 0, 0 }
  };
  const size_t count = sizeof table / sizeof *table;
  const char* result = 0;
  size_t index = 0;

  do {
    const Entry* const entry = table + index;

    if ( entry->variable && ! strcmp( entry->variable, variable ) ) {
      result = entry->file_variable;
    }

    ++index;
  } while ( index < count && result == 0 );

  if ( result == 0 ) {
    fprintf( stderr, "\n\nFailure: unknown variable '%s'.\n\n", variable );
  }

  return result;
}



/* Forward declarations: */

static void usage( const char* program );

static int parse_options( const int argc, const char* argv[],
                          const char** variable, const char** units,
                          double* minimum, double* maximum,
                          int* yyyymmddhh, int* hours_per_timestep,
                          int* timesteps,
                          const char** files );

static int process_files( const char* const variable, const char* const units,
                          const double minimum, const double maximum,
                          const int yyyymmddhh, const int hours_per_timestep,
                          const int timesteps, const char* files );

static void write_header( const char* const variable, const char* const units,
                          const int yyyymmddhh, const int hours_per_timestep,
                          const int timesteps,
                          const int components,
                          const int rows, const int columns,
                          const double west, const double east,
                          const double south, const double north );

static int write_data( const size_t count, float data[] );

static int read_file_dimensions( const int file,
                                 int* rows, int* columns,
                                 double* west, double* east,
                                 double* south, double* north );

static int read_data( const int file, const char* const variable,
                      const double minimum, const double maximum,
                      const int rows, const int columns, float data[] );

static void fill_data( const size_t count, const float value, float data[] );

static void filter_data( const size_t count,
                         const double minimum, const double maximum,
                         float data[] );

static void reverse_4_byte_words_if_little_endian( const size_t count,
                                                   void* array );

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours );

static int days_in_month( const int year, const int month );





/* Read arguments, process input data files and write data to stdout: */

int main( int argc, char* argv[] ) {
  const char* variable = 0;
  const char* units = 0;
  double minimum = 0.0;
  double maximum = 0.0;
  int yyyymmddhh = 0;
  int hours_per_timestep = 0;
  int timesteps = 0;
  const char* files = 0;
  int result =
    parse_options( argc, (const char**) argv,
                   &variable, &units, &minimum, &maximum,
                   &yyyymmddhh, &hours_per_timestep, &timesteps, &files );

  if ( result ) {
    result =
      process_files( variable, units, minimum, maximum,
                     yyyymmddhh, hours_per_timestep, timesteps, files );
  }

  return ! result;
}



/* Print program usage instructions to stderr: */

static void usage( const char* program ) {
  fprintf( stderr,
           "\n%s - Read a sequence of NetCDF gridded NLDAS files\n"
           "and convert it to bin format and write it to stdout.\n",
           program );
  fprintf( stderr, "usage: %s variable units minimum maximum "
                   "yyyymmddhh hours_per_timestep timesteps input_files "
                   "> output.bin\n", program );
  fprintf( stderr,
           "example: %s wind m/s -500 500 2016072900 1 3 "
           "test/wind_files > test/wind.bin\n", program );
  fprintf( stderr, "head -7 test/wind.bin\n\n" );
}



/* Read/check command-line options: */

static int parse_options( const int argc, const char* argv[],
                          const char** variable, const char** units,
                          double* minimum, double* maximum,
                          int* yyyymmddhh, int* hours_per_timestep,
                          int* timesteps,
                          const char** files ) {

  int result = 0;
  *variable = *units = 0;
  *minimum = *maximum = 0.0;
  *yyyymmddhh = *hours_per_timestep = *timesteps = 0;
  *files = 0;

  if ( argc == 9 ) {
    *variable = argv[ 1 ];

    if ( *variable ) {
      *units = argv[ 2 ];

      if ( *units ) {
        *minimum = atof( argv[ 3 ] );
        *maximum = atof( argv[ 4 ] );

        if ( *maximum > *minimum ) {
          *yyyymmddhh = atoi( argv[ 5 ] );

          {
            const int yyyy = *yyyymmddhh / 1000000;
            const int mm = *yyyymmddhh / 10000 % 100;
            const int dd = *yyyymmddhh / 100 % 100;
            const int hh = *yyyymmddhh % 100;

            if ( IN_RANGE( yyyy, 1900, 3000 )&&
                 IN_RANGE( mm, 1, 12 ) &&
                 IN_RANGE( dd, 1, days_in_month( yyyy, mm ) ) &&
                 IN_RANGE( hh, 0, 23 ) ) {
              *hours_per_timestep = atoi( argv[ 6 ] );

              if ( *hours_per_timestep >= 1 ) {
                *timesteps = atoi( argv[ 7 ] );

                if ( *timesteps >= 1 ) {
                  *files = argv[ 8 ];
                  result = *files[ 0 ] != '\0';
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
    usage( argv[ 0 ] );
  }

  return result;
}



/* Process each input file and write the converted data to stdout: */

static int process_files( const char* const variable, const char* const units,
                          const double minimum, const double maximum,
                          const int yyyymmddhh, const int hours_per_timestep,
                          const int timesteps, const char* files ) {
  int result = 0;
  const int components = ! strcmp( variable, "wind" ) ? 2 : 1;
  FILE* inputs = fopen( files, "r" );

  if ( inputs ) {
    double west = 0.0;
    double east = 0.0;
    double south = 0.0;
    double north = 0.0;
    int rows = 0;
    int columns = 0;
    float* data = 0;
    int timestep = 0;
    size_t timestepSize = 0;
    size_t componentSize = 0;
    size_t dataSize = 0;

    for ( timestep = 0; timestep < timesteps; ++timestep ) {
      int ok = 0;
      char file_name[ 256 ] = "";
      memset( file_name, 0, sizeof file_name );

      if ( fgets( file_name, sizeof file_name / sizeof *file_name, inputs ) ) {
        char* const newline = strrchr( file_name, '\n' );
        int input = -1;
        int status = 0;

        if ( newline ) {
          *newline = '\0';
        }

        status = nc_open( file_name, NC_NOWRITE, &input );

        DEBUG( fprintf( stderr, "file_name = '%s', input = %d\n",
                        file_name, input ); )

        if ( status != NC_NOERR ) {
          const char* const message = nc_strerror( status );
          fprintf( stderr, "Can't open file '%s' because %s.\n",
                   file_name, message );
        } else {

          if ( data == 0 ) {

            if ( read_file_dimensions( input, &rows, &columns,
                                       &west, &east, &south, &north ) ) {
              timestepSize = rows * columns;
              componentSize = timesteps * timestepSize;
              dataSize = components * componentSize;
              data = NEW( float, dataSize );

              if ( data ) {
                fill_data( dataSize, MISSING, data );
              }
            }
          }

          if ( data ) {
            const char* variable0 = variable;
            float* const udata = data + timestep * timestepSize;

            if ( ! strcmp( variable, "wind" ) ) {
              variable0 = "u_wind";
            }

            ok =
              read_data( input, variable0, minimum, maximum, rows, columns,
                         udata );

            if ( ok ) {

              if ( ! strcmp( variable, "wind" ) ) {
                float* const vdata = udata + componentSize;
                ok =
                  read_data( input, "v_wind", minimum, maximum, rows, columns,
                             vdata );
              }

              if ( ok ) {
                result = 1; /* At least one file was successfully processed.*/
              }
            }
          }

          nc_close( input );
          input = -1;
        }
      }
    }

    if ( result ) {
      write_header( variable, units,
                    yyyymmddhh, hours_per_timestep, timesteps,
                    components, rows, columns,
                    west, east, south, north );
      result = write_data( dataSize, data );
    }

    FREE( data );
  }

  return result;
}



/* Write grid header and timestamps to stdout: */

static void write_header( const char* const variable, const char* const units,
                          const int yyyymmddhh, const int hours_per_timestep,
                          const int timesteps,
                          const int components,
                          const int rows, const int columns,
                          const double west, const double east,
                          const double south, const double north ) {

  int timestep = 0;
  int timestamp = yyyymmddhh;
  printf( "Content-type: application/octet-stream; charset=iso-8859-1\n" );
  printf( "# variable units:\n%s %s\n", variable, units );
  printf( "# dimensions: components timesteps rows columns "
          "lonmin lonmax latmin latmax\n"
          "%-5d %5d %10d %10d %24.18f %24.18f %24.18f %24.18f\n",
          components, timesteps, rows, columns, west, east, south, north );
  printf( "# char yyyymmddhh[timesteps][11] and\n" );
  printf( "# IEEE-754 32-bit float "
          "data[components][timesteps][rows][columns]:\n" );

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    printf( "%010d\n", timestamp );
    increment_yyyymmddhh( &timestamp, hours_per_timestep );
  }
}



/* Write one timestep-worth of grid data to stdout: */

static int write_data( const size_t count, float data[] ) {

  /* Write binary MSB 32-bit IEEE-754 data: */

  int result = 0;
  reverse_4_byte_words_if_little_endian( count, data ); /* Big-endian.*/
  result = fwrite( data, sizeof (float), count, stdout ) == count;

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write all %lu bytes of data.\n",
             count * sizeof (float) );
    perror( "because" );
  }

  return result;
}



/* Read NetCDF file dimensions and lon-lat extent: */

static int read_file_dimensions( const int file,
                                 int* rows, int* columns,
                                 double* west, double* east,
                                 double* south, double* north ) {
  const double default_cell_size = 0.125; /* If degenerate use 1/8 degree. */
  int result = 1;
  int index = 0;

  for ( index = 0; result && index < 2; ++index ) {
    int id = -1;
    int status = 0;

    if ( index == 0 ) {
      status = nc_inq_varid( file, "lon", &id );

      if ( status != NC_NOERR ) {
        status = nc_inq_varid( file, "longitude", &id );
      }
    } else {
      status = nc_inq_varid( file, "lat", &id );

      if ( status != NC_NOERR ) {
        status = nc_inq_varid( file, "latitude", &id );
      }
    }

    result = 0;

    if ( status == NC_NOERR ) {
      int rank = 0;
      int dim_ids[ 32 ];
      size_t size = 0;
      float data[ 2 ] = { 0.0, 0.0 };
      nc_type nctype = (nc_type) 0;
      status = nc_inq_var( file, id, 0, &nctype, &rank, dim_ids, 0 );

      if ( status == NC_NOERR && nctype == NC_FLOAT && rank == 1 ) {
        status = nc_inq_dimlen( file, dim_ids[ 0 ], &size );

        if ( status == NC_NOERR && size > 0 ) {
          const size_t start[ 1 ] = { 0 };
          size_t count[ 1 ] = { 2 };

          if ( size == 1 ) {
            count[ 0 ] = 1;
          }

          status = nc_get_vara_float( file, id, start, count, data );

          if ( size == 1 ) {
            data[ 1 ] = data[ 0 ] + default_cell_size;
          }

          if ( status == NC_NOERR && data[ 0 ] < data[ 1 ] ) {
            const double delta = data[ 1 ] - data[ 0 ];
            const double half_delta = delta * 0.5;
            const double minimum = data[ 0 ] - half_delta;
            const double maximum = minimum + size * delta;

            if ( index == 0 &&
                 IN_RANGE( minimum, -180.0, 180.0 ) &&
                 IN_RANGE( maximum, minimum, 180.0 ) ) {
              *columns = size;
              *west = minimum;
              *east = maximum;
              result = 1;
            } else if ( IN_RANGE( minimum, -90.0, 90.0 ) &&
                        IN_RANGE( maximum, minimum, 90.0 ) ) {
              *rows = size;
              *south = minimum;
              *north = maximum;
              result = 1;
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    *rows = *columns = 0;
    *west = *east = *south = *north = MISSING;
  }

  DEBUG( fprintf( stderr, "read dimensions: [%d %d], [%lf %f] [%lf %lf]\n",
                  *columns, *rows, *west, *east, *south, *north ); )
  return result;
}



/* Read one timestep worth of input file data into data[]: */

static int read_data( const int file, const char* const variable,
                      const double minimum, const double maximum,
                      const int rows, const int columns, float data[] ) {
  int result = 0;
  const char* const file_variable = lookup_file_variable_name( variable );
  const size_t data_size = rows * columns;

  fill_data( data_size, MISSING, data );

  if ( file_variable ) {
    int id = -1;
    int status = nc_inq_varid( file, file_variable, &id );

    if ( status == NC_NOERR ) {
      int rank = 0;
      int dim_ids[ 32 ];
      nc_type nctype = (nc_type) 0;
      status = nc_inq_var( file, id, 0, &nctype, &rank, dim_ids, 0 );

      if ( status == NC_NOERR && nctype == NC_FLOAT &&
           ( rank == 3 || rank == 4 ) ) {
        const size_t start[ 4 ] = { 0, 0, 0, 0 };
        size_t count[ 4 ] = { 1, 1, 0, 0 };
        count[ 2 ] = rows;
        count[ 3 ] = columns;
        status =
          nc_get_vara_float( file, id, start,
                             rank == 4 ? count : count + 1, data );
      }
    }

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "Can't read NetCDF file variable '%s' because %s.\n",
               variable, message );
    } else {
      result = 1;
    }
  }

  filter_data( data_size, minimum, maximum, data );
  return result;
}



/* Fill data[ count ] with MISSING: */

static void fill_data( const size_t count, const float value, float data[] ) {
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    data[ index ] = value;
  }
}



/* Filter data[ count ] by [minimum, maximum], assigning MISSING otherwise: */

static void filter_data( const size_t count,
                         const double minimum, const double maximum,
                         float data[] ) {
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const float value = data[ index ];

    if ( ! IN_RANGE( value, minimum, maximum ) ) {
      data[ index ] = MISSING;      
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





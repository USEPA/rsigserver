/******************************************************************************
PURPOSE: bounds_filter.c - Filter a list of VIIRS NetCDF4 files by bounds.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).
         Compile:
         gcc -Wall -DNDEBUG -O -o bounds_filter bounds_filter.c \
                   -I../../../include/NetCDF4 \
                   -lnetcdf4 -lhdf5_hl -lhdf5 -lcurl -lz -ldl -lm -lc
         strip bounds_filter

         Usage:
         bounds_filter -files <listfile> \
                      -domain <minimum_longitude> <minimum_latitude> \
                              <maximum_longitude> <maximum_latitude> \

          Example:
          bounds_filter \
          -files testdata/file_list -domain -75 35 -70 36

          Prints a filtered list of files whose swath bounds intersect domain.

HISTORY: 2021-12-18 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h>    /* For macro assert(). */
#include <stdio.h>     /* For FILE, stderr, fprintf(), puts(). */
#include <string.h>    /* For memset(), strcmp(), strchr(). */
#include <stdlib.h>    /* For malloc(), free(), strtod(). */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

/* Declare just the part of NetCDF Library needed: */

enum { NC_NOERR = 0, NC_NOWRITE = 0, NC_GLOBAL = -1 };
extern int nc_open( const char* path, int mode, int* id );
extern int nc_close( int id );
extern int nc_get_att_float(int ncid, int varid, const char* name, float* val);
extern const char* nc_strerror( int ncerr );

/*================================== MACROS =================================*/

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

/*================================== TYPES ==================================*/

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };
typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE,LATITUDE ][ MINIMUM,MAXIMUM ]*/

/*========================== FORWARD DECLARATIONS ===========================*/

static void print_usage( const char* const name );
static int parse_arguments( int argc, char* argv[],
                            const char** list_file, Bounds domain );
static int process_files( const char* const list_file, const Bounds domain );
static int read_file_bounds( const char* const file_name, Bounds bounds );
static int read_viirs_file_bounds( const int file, Bounds bounds );
static int read_tropomi_file_bounds( const int file, Bounds bounds );
static int is_valid_bounds( const Bounds bounds );
static int bounds_overlap( const Bounds a, const Bounds b );
static char* read_file( const char* name );
static size_t file_size( const char* name );

/*================================ FUNCTIONS ================================*/


int main( int argc, char* argv[] ) {
  const char* list_file = 0;
  Bounds domain = { { 0.0, 0.0 }, { 0.0, 0.0 } };
  int ok = parse_arguments( argc, argv, &list_file, domain );

  if ( ! ok ) {
    print_usage( argv[ 0 ] );
  } else {
    ok = process_files( list_file, (const double (*)[2]) domain ) > 0;
  }

  return ! ok;
}



static void print_usage( const char* name ) {
  assert( name ); assert( *name );
  fprintf( stderr,
           "\a\n\n%s - "
           "Filter a list of VIIRS/TROPOMI L2 NetCDF4 files by domain.\n",
           name );
  fprintf( stderr, "\nUsage:\n%s -files <listfile> "
           "  -domain <minimum_longitude> <minimum_latitude> "
           "<maximum_longitude> <maximum_latitude>\n\n", name );
  fprintf( stderr, "Example:\n\n%s -files testdate/file_list "
           "-domain -59.5 40 -59 41\n", name );
  fprintf( stderr,
           "testdata/JRR-AOD_v1r1_npp_"
           "s201708011806009_e201708011807251_c201708011857490.nc\n" );
}



static int parse_arguments( int argc, char* argv[],
                            const char** list_file, Bounds domain ) {
  int result = 0;
  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] );
  assert( list_file );
  assert( domain );
  *list_file = 0;
  memset( domain, 0, sizeof (Bounds) );

  result =
    argc == 8 &&
    ! strcmp( argv[ 1 ], "-files" ) &&
    argv[ 2 ] &&
    argv[ 2 ][ 0 ] &&
    ! strcmp( argv[ 3 ], "-domain" );

  if ( result ) {
    char* end = 0;
    *list_file = argv[ 2 ];
    domain[ LONGITUDE ][ MINIMUM ] = strtod( argv[ 4 ], &end );
    result = end != argv[ 4 ];

    if ( result ) {
      domain[ LATITUDE ][ MINIMUM ] = strtod( argv[ 5 ], &end );
      result = end != argv[ 5 ];

      if ( result ) {
        domain[ LONGITUDE ][ MAXIMUM ] = strtod( argv[ 6 ], &end );
        result = end != argv[ 6 ];

        if ( result ) {
          domain[ LATITUDE ][ MAXIMUM ] = strtod( argv[ 7 ], &end );
          result = end != argv[ 7 ];

          if ( result ) {
            result = is_valid_bounds( (const double (*)[2]) domain );
          }
        }
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nInvalid/insufficient command-line arguments.\n" );
  }

  return result;
}



static int process_files( const char* const list_file, const Bounds domain ) {
  int result = 0;
  assert( list_file ); assert( is_valid_bounds( domain ) );

  DEBUG( fprintf( stderr, "domain = [%f %f][%f %f]\n",
                  domain[ LONGITUDE ][ MINIMUM ],
                  domain[ LONGITUDE ][ MAXIMUM ],
                  domain[ LATITUDE  ][ MINIMUM ],
                  domain[ LATITUDE  ][ MAXIMUM ] ); )

  {
    char* file_list = read_file( list_file );

    if ( file_list ) {
      char* file_name = file_list;

      do {
        char* const newline = strchr( file_name, '\n' );

        if ( newline ) {
          *newline = '\0';

          {
            Bounds file_bounds = { { 0.0, 0.0 }, { 0.0, 0.0 } };
            int ok = read_file_bounds( file_name, file_bounds );

            if ( ok ) {
              const int file_in_domain =
                bounds_overlap( domain, (const double (*)[2]) file_bounds );

              if ( file_in_domain ) {
                puts( file_name );
                ++result;
              }
            }
          }

          file_name = newline + 1;
        } else {
          file_name = 0;
        }

      } while ( file_name );

      free( file_list ), file_list = 0;
    }
  }

  return result;
}



static int read_file_bounds( const char* const file_name, Bounds bounds ) {
  int result = 0;
  assert( file_name ); assert( bounds );
  memset( bounds, 0, sizeof (Bounds) );

  {
    int file = -1;
    int status = nc_open( file_name, NC_NOWRITE, &file );

    if ( status == NC_NOERR ) {

      if ( strstr( file_name, "JRR-AOD_" ) ) {
        result = read_viirs_file_bounds( file, bounds );
      } else {
        result = read_tropomi_file_bounds( file, bounds );
      }

      nc_close( file );
    } else {
      const char* const message = nc_strerror( status );
      fprintf( stderr,
               "Failed to open NetCDF file %s for reading because: %s\n",
               file_name, message );
      result = 0;
    }
  }

  return result;
}



static int read_viirs_file_bounds( const int file, Bounds bounds ) {
  int result = 0;
  assert( file >= 0 ); assert( bounds );
  memset( bounds, 0, sizeof (Bounds) );

  {
    const char* const attributes[ 8 ] = {
      "geospatial_first_scanline_first_fov_lon",
      "geospatial_first_scanline_last_fov_lon",
      "geospatial_last_scanline_first_fov_lon",
      "geospatial_last_scanline_last_fov_lon",
      "geospatial_first_scanline_first_fov_lat",
      "geospatial_first_scanline_last_fov_lat",
      "geospatial_last_scanline_first_fov_lat",
      "geospatial_last_scanline_last_fov_lat"
    };
    float values[ 8 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    const int count = sizeof attributes / sizeof *attributes;
    int index = 0;
    int status = NC_NOERR;
    result = 0;

    for ( index = 0; index < count && status == NC_NOERR; ++index ) {
      status =
        nc_get_att_float( file, NC_GLOBAL, attributes[ index ],
                          values + index );
    }

    if ( status == NC_NOERR ) {
      bounds[ LONGITUDE ][ MINIMUM ] =
      bounds[ LONGITUDE ][ MAXIMUM ] = values[ 0 ];
      bounds[ LATITUDE  ][ MINIMUM ] =
      bounds[ LATITUDE  ][ MAXIMUM ] = values[ 4 ];

      for ( index = 1; index < 4; ++index ) {
        const float longitude = values[ index ];
        const float latitude  = values[ index + 4 ];

        if ( longitude < bounds[ LONGITUDE ][ MINIMUM ] ) {
          bounds[ LONGITUDE ][ MINIMUM ] = longitude;
        } else if ( longitude > bounds[ LONGITUDE ][ MAXIMUM ] ) {
          bounds[ LONGITUDE ][ MAXIMUM ] = longitude;
        }

        if ( latitude < bounds[ LATITUDE ][ MINIMUM ] ) {
          bounds[ LATITUDE ][ MINIMUM ] = latitude;
        } else if ( latitude > bounds[ LATITUDE ][ MAXIMUM ] ) {
          bounds[ LATITUDE ][ MAXIMUM ] = latitude;
        }
      }

      /*
       * Some VIIRS files contain some bogus (-999.3) longitude/latitudes
       * on the edges of the swath!
       * This is indicated by the geospatial bounds having invalid range.
       * Allow such files and expect the Subset program to clamp/filter.
       * Here just clamp invalid coordinates to valid range.
       */

      DEBUG( fprintf( stderr, "read file bounds = [%f %f][%f %f]\n",
                      bounds[ LONGITUDE ][ MINIMUM ],
                      bounds[ LONGITUDE ][ MAXIMUM ],
                      bounds[ LATITUDE  ][ MINIMUM ],
                      bounds[ LATITUDE  ][ MAXIMUM ] ); )

      if ( ! IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ) ) {
        bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
      }

      if ( ! IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
                       bounds[ LONGITUDE ][ MINIMUM ], 180.0 ) ) {
        bounds[ LONGITUDE ][ MAXIMUM ] = 180.0;
      }

      if ( ! IN_RANGE( bounds[ LATITUDE ][ MINIMUM ], -90.0, 90.0 ) ) {
        bounds[ LATITUDE ][ MINIMUM ] = -90.0;
      }

      if ( ! IN_RANGE( bounds[ LATITUDE ][ MAXIMUM ],
                       bounds[ LATITUDE ][ MINIMUM ], 90.0 ) ) {
        bounds[ LATITUDE ][ MAXIMUM ] = 90.0;
      }

      result = 1;
    }

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr,
               "Failed to read NetCDF file %d geospatial bounds because: %s\n",
               file, message );
      result = 0;
    }
  }

  assert( ! result || is_valid_bounds( (const double (*)[2]) bounds ) );
  return result;
}



/* Also works for TEMPO files. */

static int read_tropomi_file_bounds( const int file, Bounds bounds ) {
  int result = 0;
  assert( file >= 0 ); assert( bounds );
  memset( bounds, 0, sizeof (Bounds) );

  {
    float longitude_minimum = 0.0;
    float longitude_maximum = 0.0;
    float latitude_minimum  = 0.0;
    float latitude_maximum  = 0.0;
    int status =
      nc_get_att_float( file, NC_GLOBAL, "geospatial_lon_min",
                        &longitude_minimum );

    if ( status == NC_NOERR ) {
      status =
        nc_get_att_float( file, NC_GLOBAL, "geospatial_lon_max",
                          &longitude_maximum );

      if ( status == NC_NOERR ) {
        status =
          nc_get_att_float( file, NC_GLOBAL, "geospatial_lat_min",
                            &latitude_minimum );

        if ( status == NC_NOERR ) {
          status =
            nc_get_att_float( file, NC_GLOBAL, "geospatial_lat_max",
                              &latitude_maximum );
        }
      }
    }

    if ( status == NC_NOERR ) {
      DEBUG( fprintf( stderr, "read file bounds = [%f %f][%f %f]\n",
                      longitude_minimum, longitude_maximum,
                      latitude_minimum, latitude_maximum ); )

      /* Sometimes the attributes are not ordered min <= max so fix here: */

      if ( longitude_minimum > longitude_maximum ) {
        const float swap = longitude_minimum;
        longitude_minimum = longitude_maximum;
        longitude_maximum = swap;
      }

      if ( latitude_minimum > latitude_maximum ) {
        const float swap = latitude_minimum;
        latitude_minimum = latitude_maximum;
        latitude_maximum = swap;
      }

      /* Clamp/expand to valid range: */

      if ( ! IN_RANGE( longitude_minimum, -180.0, 180.0 ) ) {
        longitude_minimum = -180.0;
      }

      if ( ! IN_RANGE( longitude_maximum, longitude_minimum, 180.0 ) ) {
        longitude_maximum = 180.0;
      }

      if ( ! IN_RANGE( latitude_minimum, -90.0, 90.0 ) ) {
        latitude_minimum = -90.0;
      }

      if ( ! IN_RANGE( latitude_maximum, latitude_minimum, 90.0 ) ) {
        latitude_maximum = 90.0;
      }

      bounds[ LONGITUDE ][ MINIMUM ] = longitude_minimum;
      bounds[ LONGITUDE ][ MAXIMUM ] = longitude_maximum;
      bounds[ LATITUDE  ][ MINIMUM ] = latitude_minimum;
      bounds[ LATITUDE  ][ MAXIMUM ] = latitude_maximum;
      result = 1;
    } else {
      const char* const message = nc_strerror( status );
      fprintf( stderr,
               "Failed to read NetCDF file %d geospatial bounds because: %s\n",
               file, message );
      result = 0;
    }
  }

  assert( ! result || is_valid_bounds( (const double (*)[2]) bounds ) );
  return result;
}



static int is_valid_bounds( const Bounds bounds ) {
  const int result =
    bounds != 0 &&
    IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ) &&
    IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
              bounds[ LONGITUDE ][ MINIMUM ], 180.0 ) &&
    IN_RANGE( bounds[ LATITUDE ][ MINIMUM ], -90.0, 90.0 ) &&
    IN_RANGE( bounds[ LATITUDE ][ MAXIMUM ],
              bounds[ LATITUDE ][ MINIMUM ], 90.0 );
  return result;
}



static int bounds_overlap( const Bounds a, const Bounds b ) {
  int result = 0;
  assert( is_valid_bounds( a ) ); assert( is_valid_bounds( b ) );

  {
    const int outside =
      a[ LATITUDE  ][ MINIMUM ] > b[ LATITUDE  ][ MAXIMUM ] ||
      a[ LATITUDE  ][ MAXIMUM ] < b[ LATITUDE  ][ MINIMUM ] ||
      a[ LONGITUDE ][ MINIMUM ] > b[ LONGITUDE ][ MAXIMUM ] ||
      a[ LONGITUDE ][ MAXIMUM ] < b[ LONGITUDE ][ MINIMUM ];

    result = ! outside;
  }

  return result;
}



static char* read_file( const char* name ) {
  char* result = 0;
  assert( name );

  {
    const size_t length = file_size( name ) / sizeof (char);

    if ( length > 0 ) {
      const size_t bytes = ( length + 1 ) * sizeof (char);
      result = malloc( bytes );

      if ( ! result ) {
        fprintf( stderr,
                 "\nCan't allocate %lu bytes "
                 "to complete the requested action.\n",
                 bytes );
      } else {
        FILE* file = fopen( name, "rb" );

        if ( file ) {
          const size_t items_read = fread( result, sizeof (char), length, file);

          if ( items_read != length ) {
            fprintf( stderr, "\nFailed to read entire file '%s'.\n", name );
            free( result );
            result = 0;
          } else {
            result[ length ] = '\0'; /* Terminate string. */
          }

          fclose( file );
          file = 0;
        }
      }
    }
  }

  return result;
}



static size_t file_size( const char* name ) {
  size_t result = 0;
  struct stat buf;
  assert( name );

  if ( stat( name, &buf ) == -1 ) {
    fprintf( stderr, "\nFailed to determine size of file '%s'.\n", name );
  } else {
    result = buf.st_size;

    if ( result < 1 ) {
      fprintf( stderr, "\nNegative size of file '%s'.\n", name );
      result = 0;
    }
  }

  return result;
}





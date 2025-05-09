/*
read_lonlats.c - Read GOES longitudes/latitudes from a file.
gcc -Wall -g -o read_lonlats read_lonlats.c
time read_lonlats data/sport_goesImager_latLon_20130919.txt
*/

#include <stdio.h> /* For fgets(), fscanf(), sscanf(). */
#include <stdlib.h> /* For malloc(), free(). */
#include <string.h> /* For memset(), strlen(), strncmp(), strchr(). */


static int parse_arguments( const int argc, char* argv[],
                            const char** lonlat_file );

static double* read_lonlats( const char* const file_name,
                             int* rows, int* columns );

static int skip_lines( FILE* file, const int lines );

static int read_count( FILE* file, const char* const tag, int* line );

static int read_coordinate_pair( FILE* file,
                                 double* longitude, double* latitude );


int main( int argc, char* argv[] ) {
  int ok = 0;
  const char* lonlat_file = 0;

  if ( parse_arguments( argc, argv, &lonlat_file ) ) {
    int rows = 0;
    int columns = 0;
    double* lonlats = read_lonlats( lonlat_file, &rows, &columns );

    if ( lonlats ) {
      printf( "coordinates: %d rows x %d columns [%lf ... %lf][%lf ... %lf]\n",
              rows, columns,
              lonlats[ 0 ], lonlats[ rows * columns - 1 ],
              lonlats[ rows * columns ], lonlats[ 2 * rows * columns - 1 ] );
      free( lonlats ), lonlats = 0;
      ok = 1;
    }
  }

  return ! ok;
}


static int parse_arguments( const int argc, char* argv[],
                            const char** lonlat_file ) {
  int result = 0;
  *lonlat_file = 0;

  if ( argc == 2 && argv && argv[ 0 ] && argv[ 1 ] && argv[ 1 ][ 0 ] &&
       strcmp( argv[ 0 ], argv[ 1 ] ) ) {
    *lonlat_file = argv[ 1 ];
  } else {
    fprintf( stderr, "\nUsage: read_lonlats lonlat_file\n" );
    fprintf( stderr, "Example: read_lonlats "
             "data/sport_goesImager_latLon_20130919.txt\n\n" );
  }

  result = *lonlat_file != 0;
  return result;
}


static double* read_lonlats( const char* const file_name,
                             int* rows, int* columns ) {
  double* result = 0;
  FILE* file = fopen( file_name, "r" );
  *columns = *rows = 0;

  if ( file ) {
    int line = 0;
    const int header_lines = read_count( file, "hdr_lines", &line );

    if ( header_lines ) {
      *columns = read_count( file, "NX", &line );

      if ( *columns ) {
        *rows = read_count( file, "NY", &line );

        if ( *rows && *rows * *columns > 0 ) {
          const int points = *rows * *columns;
          result = malloc( points * 2 * sizeof *result );

          if ( result ) {
            int ok = skip_lines( file, header_lines - line );
            int point = 0;
            double* const longitudes = result;
            double* const latitudes  = result + points;

            for ( point = 0; ok && point < points; ++point ) {
              ok = read_coordinate_pair( file,
                                         longitudes + point, latitudes + point);
            }

            if ( ! ok ) {
              free( result ), result = 0;
            }
          }
        }
      }
    } 

    fclose( file ), file = 0;
  }

  return result;
}


static int skip_lines( FILE* file, const int lines ) {
  int result = 0;
  int count = 0;

  do {
    char line[ 256 ] = "";
    memset( line, 0, sizeof line );

    if ( ! fgets( line, sizeof line / sizeof *line, file ) ) {
      count = lines;
    }

    ++count;
  } while ( count < lines );

  result = count == lines;
  return result;
}


static int read_count( FILE* file, const char* const tag, int* line ) {
  int result = 0;
  int done = 0;
  const int tag_length = strlen( tag );

  do {
    char buffer[ 256 ] = "";
    memset( buffer, 0, sizeof buffer );

    if ( ! fgets( buffer, sizeof buffer / sizeof *buffer, file ) ) {
      done = 1;
    } else {
      const int found_line = ! strncmp( buffer, tag, tag_length );
      *line += 1;

      if ( found_line ) {
        const char* const colon = strchr( buffer, ':' );

        if ( colon ) {
          const int ok = sscanf( colon + 1, "%d", &result ) == 1 && result > 0;

          if ( ! ok ) {
            result = 0;
          }
        }

        done = 1;
      }
    }

  } while ( ! done );

  return result;
}


static int read_coordinate_pair( FILE* file,
                                 double* longitude, double* latitude ) {
  const int result =
    fscanf( file, "%*s %*s %lf %lf\n", latitude, longitude ) == 2 &&
    *longitude >= -180.0 && *longitude <= 180.0 &&
    *latitude >= -90.0 && *latitude <= 90.0;

  if ( ! result ) {
    *longitude = *latitude = 0.0;
  }

  return result;
}



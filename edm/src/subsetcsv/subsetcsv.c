/*
subsetcsv.c - C program to extract time subset of rows from a WMOST csv file.
2017-02-16 plessel.todd@epa.gov
To compile:
$ cc -o subsetcsv subsetcsv.c
Usage:
subsetcsv input_csv_file yyyymmddhh1 yyyymmddhh2 output_csv_file \
 [-subset lonmin lonmax latmin latmax] [-layer Surface|Middle|Bottom]
Example:
subsetcsv /data/land_use/hspf_charles3_loadings_tn_1961_2015.csv \
          1987123100 1988010223 \
          /data/tmp/hspf_charles3_loadings_tn.csv
subsetcsv /data/land_use/esat1_water_quality_2018-2020_gulf.csv \
          2020030100 2020033023 \
          /data/tmp/esat1_water_quality.csv -subset -83 25 -81 30 -layer Surface
*/


#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For FILE, stdoerr, scanf(), printf(), fwrite(). */
#include <stdlib.h> /* For atoi(), malloc(), free(), strtol(). */
#include <string.h> /* For strcmp(), memset(). */
#include <ctype.h> /* For isdigit(). */
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */

#define IN_RANGE( x, lower, upper ) ( (x) >= (lower) && (x) <= (upper) )

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

enum { LONGITUDE, LATITUDE };
enum { MINIMUM, MAXIMUM };
typedef double Bounds[ 2 ][ 2 ]; /* [LONGITUDE,LATITUDE][MINIMUM,MAXIMUM] */

/* Forward declarations: */

static void usage( const char* program );

static int parse_arguments( const int argc, const char* argv[],
                            const char** input_file_name,
                            int* yyyymmddhh1, int* yyyymmddhh2,
                            int* no_hru_columns,
                            const char** output_file_name );

static int parse_options( int* argc, const char* argv[],
                          Bounds bounds,
                          const char** layer );

static char* read_file( const char* file_name );

static size_t file_size( const char* file_name );

static int write_file( const char* file_name,
                       const int yyyymmddh1,
                       const int yyyymmddhh2,
                       const int no_hru_columns,
                       const Bounds bounds,
                       const char* const layer,
                       const char* csv_data );

static const char* write_header_line( FILE* file, const char* line,
                                      const int no_hru_columns,
                                      int* first_columns, int* last_columns );

static const char* process_line( FILE* file, const char* line,
                                 const int yyyymmddhh1,
                                 const int yyyymmddhh2,
                                 const int first_columns,
                                 const int last_columns,
                                 const int longitude_column,
                                 const int latitude_column,
                                 const int layer_column,
                                 const Bounds bounds,
                                 const char* const layer );

static const char* parse_timestamp( const char* line, int* yyyymmddhh );

static void get_longutude_latitude_columns( const char* csv_data,
                                            int* longitude_column,
                                            int* latitude_column );

static int get_layer_column( const char* csv_data );

static int is_valid_yyyymmddhh( const int yyyymmddhh );

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours );

static void decrement_yyyymmddhh( int* yyyymmddhh, const int hours );

static int days_in_month( const int year, const int month );

static int in_bounds( const char* line,
                      const int longitude_column,
                      const int latitude_column,
                      const Bounds bounds );

static int matches_layer( const char* line, const int layer_column,
                          const char* layer );

static const char* get_column_value( const char* line, const int column );


/* Main routine: */

int main( int argc, char* argv[] ) {
  int ok = 0;
  const char* input_file_name = 0;
  const char* output_file_name = 0;
  int yyyymmddhh1 = 0;
  int yyyymmddhh2 = 0;
  int no_hru_columns = 0;
  Bounds bounds = { { 0.0, 0.0 }, { 0.0, 0.0 } };
  const char* layer = 0;

  if ( ! parse_options( &argc, (const char**) argv, bounds, &layer ) ||
       ! parse_arguments( argc, (const char**) argv,
                          &input_file_name, &yyyymmddhh1, &yyyymmddhh2,
                          &no_hru_columns, &output_file_name ) ) {
    usage( argv[ 0 ] );
  } else {
    char* csv_data = read_file( input_file_name );

    if ( csv_data ) {
      ok = write_file( output_file_name, yyyymmddhh1, yyyymmddhh2,
                       no_hru_columns,
                       (const double(*)[2]) bounds, layer, csv_data );
      free( csv_data ), csv_data = 0;
    }
  }

  return ! ok;
}



/* Print program usage instructions to stderr: */

static void usage( const char* program ) {
  assert( program );
  fprintf( stderr, "\n%s - "
           "Extract time subset of rows from a WMOST csv file.\n", program );
  fprintf( stderr, "usage: %s input_csv_file yyyymmddhh1 yyyymmddhh2 "
           "[-no_hru_columns] output_csv_file\n", program );
  fprintf( stderr, "example: %s "
           "/data/land_use/hsfp_charles3_loadings_tn_1961_2015.csv "
           "1987123100 1988010223 "
           "/data/tmp/hspf_charles3_loadings_tn.csv\n\n", program );
  fprintf( stderr, "example: %s "
           "/data/land_use/esat1_water_quality_2018-2020_gulf.csv "
           "2020030100 2020033023 "
           "/data/tmp/esat1_water_quality.csv "
           "-subset -83 25 -81 30 -layer Surface\n\n", program );
}



/* Read/check optional command-line arguments: */

static int parse_options( int* argc, const char* argv[],
                          Bounds bounds, const char** layer ) {

  int ok = 1;

  assert( argc ); assert( argv ); assert( bounds ); assert( layer );

  *layer = 0;
  memset( bounds, 0, sizeof (Bounds) );

  if ( *argc > 2 &&
       argv[0] && argv[0][0] &&
       argv[1] && argv[1][0] &&
       argv[2] && argv[2][0] ) {

    if ( ! strcmp( argv[ *argc - 2 ], "-layer" ) ) {
      *layer = argv[ *argc - 1 ];
       *argc -= 2;
    }

    if ( *argc > 5 ) {

      if ( ! strcmp( argv[ *argc - 5 ], "-subset" ) ) {
        char* end = 0;
        const double longitude_minimum = strtod( argv[ *argc - 4 ], &end );
        ok =
          end != argv[ *argc - 4 ] &&
          longitude_minimum >= -180.0 &&
          longitude_minimum <=  180.0;

        if ( ok ) {
        const double latitude_minimum = strtod( argv[ *argc - 3 ], &end );
            ok =
            end != argv[ *argc - 3 ] &&
            latitude_minimum >= -90.0 &&
            latitude_minimum <=  90.0;

          if ( ok ) {
            const double longitude_maximum = strtod( argv[ *argc - 2 ], &end );
            ok =
              end != argv[ *argc - 2 ] &&
              longitude_maximum >  longitude_minimum &&
              longitude_maximum <= 180.0;

            if ( ok ) {
              const double latitude_maximum = strtod( argv[ *argc - 1 ], &end );
              ok =
                end != argv[ *argc - 1 ] &&
                latitude_maximum >  latitude_minimum &&
                latitude_maximum <=  90.0;

              if ( ok ) {
                bounds[ LONGITUDE ][ MINIMUM ] = longitude_minimum;
                bounds[ LONGITUDE ][ MAXIMUM ] = longitude_maximum;
                bounds[ LATITUDE  ][ MINIMUM ] = latitude_minimum;
                bounds[ LATITUDE  ][ MAXIMUM ] = latitude_maximum;
                *argc -= 5;
              }
            }
          }
        }
      }
    }
  }

  return ok;
}



/* Read/check command-line arguments: */

static int parse_arguments( const int argc, const char* argv[],
                            const char** input_file_name,
                            int* yyyymmddhh1, int* yyyymmddhh2,
                            int* no_hru_columns,
                            const char** output_file_name ) {

  int result = 0;
  assert( input_file_name ); assert( yyyymmddhh1 ); assert( yyyymmddhh2 );
  assert( no_hru_columns ); assert( output_file_name );

  *input_file_name = *output_file_name = 0;
  *no_hru_columns = *yyyymmddhh1 = *yyyymmddhh2 = 0;

  if ( ( argc == 5 || argc == 6 ) &&
       argv && argv[0] && argv[1] && argv[2] && argv[3] && argv[4] &&
       argv[0][0] && argv[1][0] && argv[2][0] && argv[3][0] && argv[4][0] &&
       ( argc == 5 || ( argv[5] && argv[5][0] ) ) ) {
    *input_file_name = argv[ 1 ];
    *yyyymmddhh1 = atoi( argv[ 2 ] );
    *yyyymmddhh2 = atoi( argv[ 3 ] );
    *no_hru_columns =
      argc == 6 && strcmp( argv[ 4 ], "-no_hru_columns" ) == 0;
    *output_file_name = argv[ argc - 1 ];
    result = strcmp( *input_file_name, *output_file_name ) != 0;
    result = result && is_valid_yyyymmddhh( *yyyymmddhh1 );
    result = result && is_valid_yyyymmddhh( *yyyymmddhh2 );
    result = result && *yyyymmddhh1 <= *yyyymmddhh2;
  }

  if ( ! result ) {
    fprintf( stderr, "\nInvalid command-line arguments.\n" );
    *input_file_name = *output_file_name = 0;
    *no_hru_columns = *yyyymmddhh1 = *yyyymmddhh2 = 0;
  }

  assert( ( result == 0 && *input_file_name == 0 && *output_file_name == 0 &&
                           *yyyymmddhh1 == 0 && *yyyymmddhh2 == 0 &&
                           *no_hru_columns == 0 ) ||
          ( result == 1 && *input_file_name && *output_file_name &&
                        strcmp( *input_file_name, *output_file_name ) != 0 &&
                        is_valid_yyyymmddhh( *yyyymmddhh1 ) &&
                        is_valid_yyyymmddhh( *yyyymmddhh2 ) &&
                        *yyyymmddhh1 <= *yyyymmddhh2 &&
                        ( *no_hru_columns == 0 ||
                          *no_hru_columns == 1 ) ) );
  return result;
}



/* Read file and return as a string: */

static char* read_file( const char* name ) {
  char* result = 0;
  assert( name ); assert( *name );

  {
    const size_t length = file_size( name ) / sizeof (char);

    if ( length ) {
      result = (char*) malloc( length + 1 );

      if ( result ) {
        FILE* file = fopen( name, "r" );

        if ( file ) {
          size_t count = 0;
          count = fread( result, length * sizeof (char), 1, file );

          if ( count != 1 ) {
            fprintf( stderr, "Failed to read entire file '%s'.\n", name );
            free( result ), result = 0;
          } else {
            result[ length ] = '\0'; /* Terminate string. */
          }

          fclose( file ), file = 0;
        }
      }
    }
  }

  return result;
}



/* Get size of file, in bytes: */

static size_t file_size( const char* name ) {
  size_t result = 0;
  struct stat buf;
  assert( name ); assert( *name );

  if ( stat( name, &buf ) == -1 ) {
    fprintf( stderr, "Failed to determine size of file '%s'.\n", name );
  } else {
    result = buf.st_size;
  }

  return result;
}


/* Write csv_data rows within the time range to the output file: */

static int write_file( const char* file_name, const int yyyymmddhh1,
                       const int yyyymmddhh2, const int no_hru_columns,
                       const Bounds bounds,
                       const char* const layer,
                       const char* csv_data ) {
  int result = 0;
  FILE* file = 0;
  assert( file_name ); assert( *file_name );
  assert( is_valid_yyyymmddhh( yyyymmddhh1 ) );
  assert( is_valid_yyyymmddhh( yyyymmddhh2 ) );
  assert( yyyymmddhh1 <= yyyymmddhh2 );
  assert( no_hru_columns == 0 || no_hru_columns == 1 );
  assert( csv_data ); assert( *csv_data );

  file = fopen( file_name, "w" );

  if ( file ) {
    int first_columns = 0;
    int last_columns = 0;
    const char* line =
      write_header_line( file, csv_data, no_hru_columns,
                         &first_columns, &last_columns );

    if ( line ) {
      int longitude_column = 0;
      int latitude_column  = 0;
      int layer_column     = 0;
      int ok = 1;

      if ( bounds[ LONGITUDE ][ MINIMUM ] < bounds[ LONGITUDE ][ MAXIMUM ] ) {
        get_longutude_latitude_columns( csv_data,
                                        &longitude_column,
                                        &latitude_column );
        ok = longitude_column > 0;
      }

      if ( layer ) {
        layer_column = get_layer_column( csv_data );
        ok = layer_column > 0;
      }

      if ( ok ) {

        while ( line ) {
          line =
            process_line( file, line, yyyymmddhh1, yyyymmddhh2,
                          first_columns, last_columns,
                          longitude_column, latitude_column, layer_column,
                          bounds, layer );

          if ( line ) {
            result = 1;
          }
        }
      }
    }

    fclose( file ), file = 0;
  }

  return result;
}



/*
 * Write header line to a file and return next line or 0 if done/failed.
 * If no_hru_columns == 1 and numbered columns are parsed (and skipped),
 * *first_columns is the number of non-numbered columns in the first part and
 * *last_columns  is the number of non-numbered columns in the last part.
 */

static const char* write_header_line( FILE* file, const char* line,
                                      const int no_hru_columns,
                                      int* first_columns, int* last_columns ) {
  const char* result = 0;
  size_t length = 0; /* Length of header line to skip over for result. */
  size_t length1 = 0; /* Length of first part of header line to write. */
  size_t length2 = 0; /* Length of last  part of header line to write. */
  size_t comma_index = 0; /* Index in line of comma. */
  int commas = 0;     /* Number of commas on line. */
  assert( file ); assert( line );
  assert( no_hru_columns == 0 || no_hru_columns == 1 );
  assert( first_columns ); assert( last_columns );
  *first_columns = *last_columns = 0;

  while ( line[ length ] != '\0' && line[ length ] != '\n' ) {

    if ( no_hru_columns ) {
      const int is_comma = line[ length ] == ',';

      if ( is_comma ) {
        const int is_numbered_column =
          length > 0 && isdigit( line[ length - 1 ] );

        if ( is_numbered_column ) {

          if ( length1 == 0 ) {
            length1 = comma_index;   
            *first_columns = commas;
          }
        } else if ( length1 > 0 && length2 == 0 ) {
          length2 = comma_index;
          *last_columns = commas;
        }

        comma_index = length;
        ++commas;
      }
    }

    ++length;
  }

  if ( length ) {
    const int has_newline = line[ length ] == '\n';
    length += has_newline;

    if ( length1 == 0 ) {

      if ( fwrite( line, sizeof (char), length, file ) == length ) {
        result = line + length;
      }
    } else {
      const size_t offset = length2; /* Last comma_index. */
      assert( length1 > 0 );
      assert( length2 > length1 ); assert( length2 < length );
      length2 = length - offset;
      *last_columns = 1 + commas - *last_columns;

      if ( fwrite( line, sizeof (char), length1, file ) == length1  ) {

        if ( fwrite( line + offset, sizeof (char), length2, file) == length2) {
          result = line + length;
        }
      }
    }
  }

  assert( *first_columns >= 0 ); assert( *last_columns >= 0 );
  return result;
}



/* Find 0-based column-index of longitude and latitude: */

static void get_longutude_latitude_columns( const char* csv_data,
                                            int* longitude_column,
                                            int* latitude_column ) {


  assert( csv_data ); assert( longitude_column ); assert( latitude_column );

  {
    char header[ 4096 ] = "";
    const char* s = 0;
    *longitude_column = *latitude_column = 0;
    memset( header, 0, sizeof header );
    strncpy( header, csv_data, sizeof header / sizeof *header - 1 );
    s = strstr( header, ",Longitude(" );

    if ( ! s ) {
      s = strstr( header, ",longitude(" );
    }

    if ( ! s ) {
      s = strstr( header, ",LONGITUDE(" );
    }

    if ( s ) {
      const char* c = header;
      int column = 1;

      while ( c != s ) {

        if ( *c == ',' ) {
          ++column;
        }

        ++c;
      }

      *longitude_column = column;

      s = strstr( header, ",Latitude(" );

      if ( ! s ) {
        s = strstr( header, ",latitude(" );
      }

      if ( ! s ) {
        s = strstr( header, ",LATITUDE(" );
      }

      if ( s ) {
        c = header;
        column = 1;

        while ( c != s ) {

          if ( *c == ',' ) {
            ++column;
          }

          ++c;
        }

        *latitude_column = column;
      } else {
        *longitude_column = 0;
      }
    }
  }
}



/* Find 0-based column-index of layer: */

static int get_layer_column( const char* csv_data ) {
  int result = 0;
  assert( csv_data );

  {
    char header[ 4096 ] = "";
    const char* s = 0;
    memset( header, 0, sizeof header );
    strncpy( header, csv_data, sizeof header / sizeof *header - 1 );
    s = strstr( header, ",Layer(" );

    if ( ! s ) {
      s = strstr( header, ",layer(" );
    }

    if ( ! s ) {
      s = strstr( header, ",LAYER(" );
    }

    if ( s ) {
      const char* c = header;
      result = 1;

      while ( c != s ) {

        if ( *c == ',' ) {
          ++result;
        }

        ++c;
      }
    }
  }

  return result;
}



/*
 * Parse the timestamp and write line to file if within time range and
 * return pointer to next line or 0 if past yyyymmddhh2 (or write failure):
 */

static const char* process_line( FILE* file, const char* line,
                                 const int yyyymmddhh1,
                                 const int yyyymmddhh2,
                                 const int first_columns,
                                 const int last_columns,
                                 const int longitude_column,
                                 const int latitude_column,
                                 const int layer_column,
                                 const Bounds bounds,
                                 const char* const layer ) {

  const char* result = 0; /* Default is to stop processing lines. */

  assert( file );
  assert( line );
  assert( is_valid_yyyymmddhh( yyyymmddhh1 ) );
  assert( is_valid_yyyymmddhh( yyyymmddhh2 ) );
  assert( yyyymmddhh1 <= yyyymmddhh2 );
  assert( first_columns >= 0 ); assert( last_columns >= 0 );

  if ( *( line - 1 ) == '\n' ) {
    int yyyymmddhh = 0;
    const char* word = parse_timestamp( line, &yyyymmddhh );
    int skip_line = yyyymmddhh < yyyymmddhh1 || yyyymmddhh > yyyymmddhh2;

    if ( ! skip_line ) {

      if ( longitude_column ) { /* Filter by bounds: */
        skip_line =
          ! in_bounds( line, longitude_column, latitude_column, bounds );
      }

      if ( ! skip_line ) {

        if ( layer ) { /* Filter by layer: */
          skip_line = ! matches_layer( line, layer_column, layer );
        }
      }
    }

    if ( skip_line ) { /* Skip line. */
      result = word;

      while ( *result != '\0' && *result != '\n' ) {
        ++result;
      }

      if ( *result == '\n' ) {
        ++result;
      } else {
        result = 0; /* No more lines. */
      }

    } else { /* Write */
      size_t length = 0;
      size_t length1 = 0; /* Length of first part. */
      size_t comma_index = 0;
      int commas = 0;
      result = line;

      do {

        if ( first_columns ) {
          const int is_comma = *result == ',';

          if ( is_comma ) {
            comma_index = length;
            ++commas;

            if ( commas == first_columns ) {
              length1 = comma_index;
            }
          }
        }

        ++result;
        ++length;
      } while ( *result != '\0' && *result != '\n' );

      if ( *result == '\n' ) {
        ++result;
        ++length;
      } else {
        result = 0; /* No more lines. */
      }

      if ( length1 == 0 ) { /* Write whole line: */

        if ( fwrite( line, sizeof (char), length, file ) != length ) {
          result = 0; /* Write failure. */
        }
      } else { /* Write first and last columns only: */
        size_t index = comma_index;
        size_t length2 = 0;
        int count = 1;

        while ( index > 0 && count < last_columns ) {
          --index;

          if ( line[ index ] == ',' ) {
            ++count;
          }
        }

        length2 = length - index;

        if ( fwrite( line, sizeof (char),  length1, file ) != length1 ) {
          result = 0; /* Write failure. */
        } else if ( fwrite( line + index, sizeof (char), length2, file )
                    != length2 ) {
          result = 0; /* Write failure. */
        }
      }
    }
  }

  return result;
}



/* Parse the timestamp and return pointer to next character to parse: */

static const char* parse_timestamp( const char* line, int* yyyymmddhh ) {
  char* next = 0;
  assert( yyyymmddhh ); assert( line );
  *yyyymmddhh = 0;

  if ( *( line - 1 ) == '\n' ) {
    const char* word = line;
    const int first = strtol( word, &next, 10 );

    if ( IN_RANGE( first, 1, 12 ) && *next == '/' ) { /* m/d/yyyy */
      const int mm = first;
      word = next + 1;

      {
        const int dd = strtol( word, &next, 10 );

        if ( IN_RANGE( dd, 1, 31 ) && *next == '/' ) {
          word = next + 1;

          {
            const int yyyy = strtol( word, &next, 10 );

            if ( IN_RANGE( yyyy, 1900, 3000 ) && *next == ' ' ) {
              word = next + 1;

              {
                const int hh = strtol( word, &next, 10 );

                if ( IN_RANGE( hh, 0, 23 ) && *next == ':' ) {
                  word = next + 1;

                  {
                    const int min = strtol( word, &next, 10 );

                    if ( IN_RANGE( min, 0, 59 ) && *next == ',' ) {
                      word = next + 1;

                      {
                        const int utc_offset = strtol( word, &next, 10 );

                        if ( IN_RANGE( utc_offset, -23, 23 ) && *next == ',') {
                          *yyyymmddhh =
                            yyyy * 1000000 + mm * 10000 + dd * 100 + hh;

                          if ( is_valid_yyyymmddhh( *yyyymmddhh ) ) {

                            if ( utc_offset < 0 ) {
                              increment_yyyymmddhh( yyyymmddhh, -utc_offset );
                            } else if ( utc_offset > 0 ) {
                              decrement_yyyymmddhh( yyyymmddhh, utc_offset );
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
    } else if ( IN_RANGE( first, 1900, 3000 ) && *next == '-' ) { /*YYYY-MM-DD*/
      const int yyyy = first;
      word = next + 1;

      {
        const int mm = strtol( word, &next, 10 );

        if ( IN_RANGE( mm, 1, 12 ) && *next == '-' ) {
          word = next + 1;

          {
            const int dd = strtol( word, &next, 10 );

            if ( IN_RANGE( dd, 1, 31 ) && ( *next == 'T' || *next == ' ') ) {
              word = next + 1;

              {
                const int hh = strtol( word, &next, 10 );

                if ( IN_RANGE( hh, 0, 23 ) ) {

                  while ( *next != '\0' && *next != '\n' && *next != ',' ) {
                    ++next;
                  }

                  *yyyymmddhh =  yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
                }
              }
            }
          }
        }
      }
    }
  }

  if ( ! is_valid_yyyymmddhh( *yyyymmddhh ) ) {
    *yyyymmddhh = 0;
  }

  return next;
}



/* Is yyyymmddhh valid? */


static int is_valid_yyyymmddhh( const int yyyymmddhh ) {
  const int yyyy = yyyymmddhh / 1000000;
  const int mm   = yyyymmddhh / 10000 % 100;
  const int dd   = yyyymmddhh / 100 % 100;
  const int hh   = yyyymmddhh % 100;
  const int result =
    IN_RANGE( yyyy, 1900, 3000 ) &&
    IN_RANGE( mm, 1, 12 ) &&
    IN_RANGE( dd, 1, days_in_month( yyyy, mm ) ) &&
    IN_RANGE( hh, 0, 23 );
  return result;
}



/* Increment yyyymmddhh by hours: */

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours ) {
  int yyyy = *yyyymmddhh / 1000000;
  int mm   = *yyyymmddhh / 10000 % 100;
  int dd   = *yyyymmddhh / 100 % 100;
  int hh   = *yyyymmddhh % 100;
  int hour = 0;
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) );

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
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) );
}



/* Decrement yyyymmddhh by hours: */

static void decrement_yyyymmddhh( int* yyyymmddhh, const int hours ) {
  int yyyy = *yyyymmddhh / 1000000;
  int mm   = *yyyymmddhh / 10000 % 100;
  int dd   = *yyyymmddhh / 100 % 100;
  int hh   = *yyyymmddhh % 100;
  int hour = 0;
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) );

  for ( hour = 0; hour < hours; ++hour ) {
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
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) );
}



/* Number of days in year/month: */

static int days_in_month( const int year, const int month ) {
  int result = 0;
  assert( IN_RANGE( year, 1900, 3000 ) ); assert( IN_RANGE( month, 1, 12 ) );

  {
    static const int days_per_month[ 2 ][ 12 ] = {
      { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
      { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
    };
    const int leap =
      month != 2 ? 0 : year % 4 == 0 && ! ( year % 100 == 0 && year % 400 != 0 );
    result = days_per_month[ leap ][ month - 1 ];
  }

  return result;
}



/* Check if line has longitude latitude within bounds: */

static int in_bounds( const char* line,
                      const int longitude_column,
                      const int latitude_column,
                      const Bounds bounds ) {
  int result = 0;
  const char* column_value = get_column_value( line, longitude_column );

  if ( column_value ) {
    char* end = 0;
    const double longitude = strtod( column_value, &end );

    if ( end != column_value ) {
      result =
        longitude >= bounds[ LONGITUDE ][ MINIMUM ] &&
        longitude <= bounds[ LONGITUDE ][ MAXIMUM ];
    }
  }

  if ( result ) {
    result = 0;
    column_value = get_column_value( line, latitude_column );

    if ( column_value ) {
      char* end = 0;
      const double latitude = strtod( column_value, &end );

      if ( end != column_value ) {
        result =
          latitude >= bounds[ LATITUDE ][ MINIMUM ] &&
          latitude <= bounds[ LATITUDE ][ MAXIMUM ];
      }
    }
  }

  return result;
}



/* Check if line has matched layer: */

static int matches_layer( const char* line, const int layer_column,
                          const char* layer ) {
  int result = 0;
  const char* const column_value = get_column_value( line, layer_column );

  if ( column_value ) {
    const size_t length = strlen( layer );

    if ( length ) {
      result = ! strncmp( column_value, layer, length );
    }
  }

  return result;
}




/* Get column value: */

static const char* get_column_value( const char* line, const int column ) {
  const char* result = 0;
  const char* next = line;
  int count = 0;

  do {
    next = strchr( next, ',' );

    if ( next ) {
      ++next;
      ++count;
    }

  } while ( next && count < column );

  if ( count == column ) {
    result = next;
  }

  return result;
}




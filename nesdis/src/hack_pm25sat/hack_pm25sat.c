/*
hack_pm25sat.c - Make pm25satYYYYMMDD.nc files M3IO-compliant.
gcc -Wall -g -I. -o hack_pm25sat hack_pm25sat.c libnetcdf.a
hack_pm25sat pm25sat*.nc
ls -l pm25sat*.ncf
*/

#include <stdio.h>  /* For stderr, fprintf(), snprintf(). */
#include <string.h> /* For memset(). */
#include <stdlib.h> /* For malloc(), free(). */
#include <assert.h> /* For assert(). */

#include <netcdf.h> /* For NC_*, nc_*(). */



static int open( const char* const name, const char rw );

static int end_definition( const int file );

static int read_dimensions( const int file,
                            size_t* const rows, size_t* const columns );

static int create_dimension( int file, const char* name, size_t size,
                             int unlimited );

static int file_dimensions( const int file );

static int file_attributes( const int file );

static int copy_file_metadata( const int input, const int output );

static int copy_dimensions( const int input, const int output );

static int copy_attributes( const int input, const int output );

static int create_variable_definitions( const int output );

static int create_tflag_data( int input, int output );

static int copy_variable_data( const int input_file,
                               const size_t rows,
                               const size_t columns,
                               float* buffer,
                               const int output_file );



int main( int argc, char* argv[] ) {
  int ok = argc >= 2 && argv[ 1 ][ 0 ];

  if ( ! ok ) {
    fprintf( stderr, "usage: %s files\n", argv[ 0 ] );
    fprintf( stderr, "example: %s pm25sat/*/pm25sat*.nc\n", argv[ 0 ] );
    fprintf( stderr, "ls -l pm25sat/*/pm25sat*.ncf\n" );
  } else {
    int arg = 0;

    for ( arg = 1; ok && arg < argc; ++arg ) {
      const char* const input_file_name = argv[ arg ];
      int input_file = open( input_file_name, 'r' );

      if ( input_file != -1 ) {
        int output_file = -1;
        char output_file_name [ 256 ] = "";
        memset( output_file_name, 0, sizeof output_file_name );
        snprintf( output_file_name,
                  sizeof output_file_name / sizeof *output_file_name - 1,
                  "%sf", input_file_name );
        output_file = open( output_file_name, 'w' );

        if ( output_file != -1 ) {
          size_t rows = 0;
          size_t columns = 0;
          ok = read_dimensions( input_file, &rows, &columns );

          if ( ok ) {
            const size_t bytes = rows * columns * sizeof (float);
            float* buffer = malloc( bytes );
            ok = buffer != 0;

            if ( ok ) {
              memset( buffer, 0, bytes );
              ok = copy_file_metadata( input_file, output_file );

              if ( ok ) {
                ok = create_tflag_data( input_file, output_file );

                if ( ok ) {
                  ok = copy_variable_data( input_file, rows, columns,
                                           buffer, output_file );
                }
              }

              free( buffer ), buffer = 0;
            }
          }
        }

        nc_close( output_file ), output_file = -1;
      }

      nc_close( input_file ), input_file = -1;
    }
  }

  return ! ok;
}



static int open( const char* const name, const char rw ) {
  int result = -1;
  assert( name ); assert( rw == 'r' || rw == 'w' );

  {
    const int status =
      rw == 'r' ? nc_open( name, NC_NOWRITE, &result )
      : nc_create( name, NC_CLOBBER, &result );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't open file '%s' because %s.\n", name, message );
      result = -1;
    } else if ( result < 0 ) {
      nc_close( result );
      fprintf( stderr, "\nInvalid id (%d) for file '%s'.\n", result, name );
      result = -1;
    }
  }

  assert( result >= -1 );
  return result;
}



static int read_dimensions( const int file,
                            size_t* const rows, size_t* const columns ) {
  int result = 0;
  assert( file >= 0 ); assert( rows ); assert( columns );

  *rows = *columns = 0;

  {
    int status = NC_NOERR;
    const int dimensions = file_dimensions( file );
    int dimension = 0;

    for ( dimension = 0; dimension < dimensions; ++dimension ) {
      char name[ NC_MAX_NAME + 1 ] = "";
      size_t dimension_size = 0;
      status = nc_inq_dim( file, dimension, name, &dimension_size );

      if ( status == NC_NOERR ) {

        if ( ! strcmp( name, "ROW" ) ) {
          *rows = dimension_size;
        } else if ( ! strcmp( name, "COL" ) ) {
          *columns = dimension_size;
        }
      }
    }

    if ( status != NC_NOERR || *rows == 0 || *columns == 0 ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nFailed to read valid file dimension because %s.\n",
               message );
      *rows = *columns = 0;
    } else {
      result = 1;
    }
  }

  assert( result == 0 || result == 1 );
  return result;
}



static int file_dimensions( const int file ) {
  int result = 0;
  assert( file >= 0 );

  {
    const int status = nc_inq_ndims( file, &result );

    if ( status != NC_NOERR || result < 1 ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't query file dimensions because %s.\n", message);
      result = 0;
    }
  }

  assert( result >= 0 );
  return result;
}



static int file_attributes( const int file ) {
  int result = 0;
  assert( file >= 0 );

  {
    const int status = nc_inq_natts( file, &result );

    if ( status != NC_NOERR || result < 1 ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't read attributes in file because %s.\n",
               message );
      result = 0;
    }
  }

  assert( result >= 0 );
  return result;
}



static int copy_file_metadata( const int input, const int output ) {
  int result = 0;
  assert( input ); assert( output );

  if ( copy_dimensions( input, output ) ) {

    if ( copy_attributes( input, output ) ) {

      if ( create_variable_definitions( output ) ) {
        result = end_definition( output );
      }
    }
  }

  return result;
}



static int copy_dimensions( const int input, const int output ) {
  int result = 0;
  assert( input >= 0 ); assert( output >= 0 );

  {
    size_t rows = 0;
    size_t columns = 0;
    result = read_dimensions( input, &rows, &columns );

    if ( result ) {
      result = create_dimension( output, "TSTEP", 1, 1 ) >= 0;

      if ( result ) {
        result = create_dimension( output, "DATE-TIME", 2, 0 ) >= 0;

        if ( result ) {
          result = create_dimension( output, "LAY", 1, 0 ) >= 0;

          if ( result ) {
            result = create_dimension( output, "VAR", 2, 0 ) >= 0;

            if ( result ) {
              result = create_dimension( output, "ROW", rows, 0 ) >= 0;

              if ( result ) {
                result = create_dimension( output, "COL", columns, 0 ) >= 0;
              }
            }
          }
        }
      }
    }
  }

  assert( result == 0 || result == 1 );
  return result;
}



static int copy_attributes( const int input, const int output ) {
  int result = 0;
  assert( input >= 0 ); assert( output >= 0 );

  {
    const int attributes = file_attributes( input );
    int attribute = 0;
    int status = NC_NOERR;

    for ( attribute = 0; status == NC_NOERR && attribute < attributes;
          ++attribute ) {
      char name[ NC_MAX_NAME ] = "";
      status = nc_inq_attname( input, NC_GLOBAL, attribute, name );

      if ( status == NC_NOERR ) {

        if ( ! strcmp( name, "VGLVLS" ) ) {
          const float values[ 2 ] = { 1.0, 0.995 };
          status =
            nc_put_att_float( output, NC_GLOBAL, name, NC_FLOAT, 2, values );
        } else {
          status = nc_copy_att( input, NC_GLOBAL, name, output, NC_GLOBAL );
        }
      }
    }

    if ( status == NC_NOERR ) {
      status =
        nc_put_att_text( output, NC_GLOBAL, "GDNAM", 16, "4US1            " );

      if ( status == NC_NOERR ) {
        status =
          nc_put_att_text( output, NC_GLOBAL, "UPNAM", 16, "hack_pm25sat    ");

        if ( status == NC_NOERR ) {
          status =
            nc_put_att_text( output, NC_GLOBAL, "VAR-LIST", 32,
                             "DAILY_PM25      DAILY_PM25UNC   " );

          if ( status == NC_NOERR ) {
            enum { MXDLEN3 = 80, MXDESC3 = 60 };
            enum { FILE_DESCRIPTION_LENGTH = MXDLEN3 * MXDESC3 };
            char filedesc[ FILE_DESCRIPTION_LENGTH + 1 ] = "";
            const char* const content = "NOAA pm25sat product";
            memset( filedesc, 0, sizeof filedesc );
            snprintf( filedesc, FILE_DESCRIPTION_LENGTH, "%-4800s", content );
            status =
              nc_put_att_text( output, NC_GLOBAL, "FILEDESC",
                               FILE_DESCRIPTION_LENGTH, filedesc );

            if ( status == NC_NOERR ) {
              const char* const history = "NOAA, US EPA"
                "                                  "
                "                                  ";
              status =
                nc_put_att_text( output, NC_GLOBAL, "HISTORY", 80, history );
            }
          }
        }
      }
    }

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nFailed to copy file attribute (%d)"
               " in file because %s.\n", attribute, message );

    } else {
      result = 1;
    }
  }

  assert( result == 0 || result == 1 );
  return result;
}



static int create_variable_definitions( const int output ) {
  int result = 0;
  int dimension_ids[ 4 ] = { 0, 3, 1, -1 }; /* TSTEP, VAR, DATE-TIME. */
  int id = -1;
  int status = 0;
  assert( output >= 0 );

  /* Create correct TFLAG variable definition: */

  status = nc_def_var( output, "TFLAG", NC_INT, 3, dimension_ids, &id );

  if ( status == NC_NOERR ) {
    status = nc_put_att_text( output, id, "units", 16, "<YYYYDDD,HHMMSS>");

    if ( status == NC_NOERR ) {
      status =
        nc_put_att_text( output, id, "long_name", 16, "TFLAG           " );

      if ( status == NC_NOERR ) {
        const char* const var_desc =
          "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS"
          "                                ";
        status = nc_put_att_text( output, id, "var_desc", 80, var_desc );
      }
    }
  }

 /* Create correct PM25 variable definition: */

  if ( status == NC_NOERR ) {
    dimension_ids[ 0 ] = 0; /* TSTEP */
    dimension_ids[ 1 ] = 2; /* LAY */
    dimension_ids[ 2 ] = 4; /* ROW */
    dimension_ids[ 3 ] = 5; /* COL */
    status = nc_def_var( output, "DAILY_PM25", NC_FLOAT, 4, dimension_ids, &id);

    if ( status == NC_NOERR ) {
      status =
        nc_put_att_text( output, id, "long_name", 16, "PM25 conc.     " );

      if ( status == NC_NOERR ) {
        status =
          nc_put_att_text( output, id, "units", 16, "ug/m3           ");

        if ( status == NC_NOERR ) {
          const char* const var_desc =
            "Fine particulates not more than 2.5 microns in diameter."
            "                        ";
          status = nc_put_att_text( output, id, "var_desc", 80, var_desc );
        }
      }
    }
  }

 /* Create correct PM25UNC variable definition: */

  if ( status == NC_NOERR ) {
    status =
      nc_def_var( output, "DAILY_PM25UNC", NC_FLOAT, 4, dimension_ids, &id );

    if ( status == NC_NOERR ) {
      status =
        nc_put_att_text( output, id, "long_name", 16, "PM25 uncertainty" );

      if ( status == NC_NOERR ) {
        status =
          nc_put_att_text( output, id, "units", 16, "ug/m3           ");

        if ( status == NC_NOERR ) {
          const char* const var_desc =
            "Uncertainty of fine particulates not more than 2.5 microns in diameter."
            "         ";
          status = nc_put_att_text( output, id, "var_desc", 80, var_desc );
        }
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't define variable because %s.\n", message );
  } else {
    result = 1;
  }


  assert( result == 0 || result == 1 );
  return result;
}



static int create_dimension( int file, const char* name, size_t size,
                             int unlimited ) {
  int result = -1;
  assert( file >= 0 ); assert( name ); assert( size > 0 );
  assert( unlimited == 0 || unlimited == 1 );

  {
    const size_t length = unlimited ? NC_UNLIMITED : size;
    const int status = nc_def_dim( file, name, length, &result );

    if ( status != NC_NOERR || result < 0 ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't define file dimension because %s.\n", message);
      result = -1;
    }
  }

  assert( result >= -1 );
  return result;
}



static int end_definition( const int file ) {
  int result = 0;
  assert( file >= 0 );

  {
    const int status = nc_enddef( file );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't end definition because %s.\n", message );
    } else {
      result = 1;
    }
  }

  assert( result == 0 || result == 1 );
  return result;
}



static int create_tflag_data( int input, int output ) {
  int result = 0;
  int status = 0;
  int yyyyddd = 0;
  assert( input >= 0 ); assert( output >= 0 );

  /* Read SDATE and use it to create TFLAG values: */

  status = nc_get_att_int( input, NC_GLOBAL, "SDATE", &yyyyddd );

  if ( status == NC_NOERR ) {
    const size_t starts[ 3 ] = { 0, 0, 0 };
    const size_t counts[ 3 ] = { 1, 2, 2 }; /* TSTEP, VAR, 2. */
    int data[ 4 ] = { 0, 0, 0, 0 }; /* YYYYDDD, HHMMSS. */
    data[ 0 ] = yyyyddd;
    data[ 1 ] = 0;
    data[ 2 ] = yyyyddd;
    data[ 3 ] = 0;
    status = nc_put_vara_int( output, 0, starts, counts, data );
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't copy variable (0) TFLAG because %s.\n", message);
  } else {
    result = 1;
  }

  return result;
}


static int copy_variable_data( const int input_file,
                               const size_t rows,
                               const size_t columns,
                               float* buffer,
                               const int output_file ) {

  int result = 0;
  int status = NC_NOERR;
  int variable = 0;
  const size_t count = rows * columns;
  const size_t starts[ 4 ] = { 0, 0, 0, 0 };
  size_t counts[ 4 ] = { 1, 1, 0, 0 }; /* TSTEP, LAY, ROW, COL. */

  assert( input_file >= 0 );
  assert( rows ); assert( columns );
  assert( buffer );
  assert( output_file >= 0 );

  counts[ 2 ] = rows;
  counts[ 3 ] = columns;

  status = nc_get_vara_float( input_file, 3, starts, counts, buffer );

  if ( status == NC_NOERR ) {
    const double BADVAL3 = -9.999e+36;
    size_t index = 0;

    for ( index = 0; index < count; ++index ) {
      const float value = buffer[ index ];

      /* Replace invalid values (ug/m3) with BADVAL3: */

      if ( value < 0.0 ) {
        buffer[ index ] = BADVAL3;
      }
    }

    status = nc_put_vara_float( output_file, 1, starts, counts, buffer );
  }

  if ( status == NC_NOERR ) {
    status = nc_get_vara_float( input_file, 5, starts, counts, buffer );

    if ( status == NC_NOERR ) {
      const double BADVAL3 = -9.999e+36;
      size_t index = 0;

      for ( index = 0; index < count; ++index ) {
        const float value = buffer[ index ];

        /* Replace invalid values (ug/m3) with BADVAL3: */

        if ( value < 0.0 ) {
          buffer[ index ] = BADVAL3;
        }
      }

      status = nc_put_vara_float( output_file, 2, starts, counts, buffer );
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't copy variable (%d) because %s.\n",
             variable, message );
  } else {
    result = 1;
  }

  return result;
}






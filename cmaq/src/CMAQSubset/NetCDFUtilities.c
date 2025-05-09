
/******************************************************************************
PURPOSE: NetCDFUtilities.c - Define convenience routines for NetCDF files.
NOTES:   For a description of NetCDF COARDS Conventions see:
         http://ferret.wrc.noaa.gov/noaa_coop/coop_cdf_profile.html
HISTORY: 2007-12-01 plessel.todd@epa.gov
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For stderr, fprintf(), puts(). */
#include <string.h> /* For memset(), strcmp(). */
#include <float.h>  /* For FLT_MAX, DBL_MAX. */
#include <limits.h>  /* For INT_MAX. */

#include <netcdf.h> /* For NC_*, nc_*. */

#include <Assertions.h>      /* For PRE*(), POST*(), DEBUG(), IN_RANGE(). */
#include <Utilities.h>       /* For TIME, LAYER, ROW, COLUMN. */
#include <NetCDFUtilities.h> /* For public interface. */

static int isNaN( double x ) {
  const double copy = x;
  const int result = (copy != x);
  return result;
}

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: printM3IOVariables - Print names of M3IO float variables to stdout
         one variable per line.
INPUTS:  const char* const name  Name of NetCDF file to examine.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int printM3IOVariables( const char* const name ) {
  PRE02( name, *name );
  enum { NAME_SIZE = 32 };
  typedef char VariableName[ NAME_SIZE ];
  enum { MAXIMUM_VARIABLES = 512 };
  static VariableName listing[ MAXIMUM_VARIABLES ]; /* 16KB. */
  size_t count = 0;
  int result = 0;
  const int file = openNetCDFFile( name, 'r' );

  if ( file ) {
    int variables = 0;
    int status = nc_inq_nvars( file, &variables );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't read file '%s' variable count because %s.\n",
               name, message );
    } else {
      int variable = 0;
      memset( listing, 0, sizeof listing );

      do {
        int type = 0;
        int rank = 0;
        int dimensions[ 32 ];
        char variableName[ 256 ] = "";
        char units[ 256 ] = "";
        memset( dimensions, 0, sizeof dimensions );
        memset( variableName, 0, sizeof variableName );
        memset( units, 0, sizeof units );
        result =
          getNetCDFVariableInfo( file, variable, variableName, &type, &rank,
                                 dimensions, units, 0 );

        if ( AND2( isNetCDFFloat( type ), rank == 4 ) ) {
          strncpy( listing[ count ], variableName,
                   sizeof *listing / sizeof (char) - 1 );
          ++count;
        }

        ++variable;
      } while ( AND2( result, variable < variables ) );
    }

    closeNetCDFFile( file );

    /* Sort listing: */

    qsort( listing, count, sizeof *listing,
           (int (*)(const void*, const void*)) strcasecmp );

    /* Print sorted listing: */

    {
      size_t index = 0;

      for ( ; index < count; ++index ) {
        puts( listing[ index ] );
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFFile - Create a NetCDF file.
INPUTS:  const char* const name  Name of NetCDF file to create.
RETURNS: int >= 0 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int createNetCDFFile( const char* const name ) {
  PRE02( name, *name );
  int result = -1;
  const int status = nc_create( name, NC_CLOBBER | NC_SHARE, &result );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't create file '%s' because %s.\n", name, message );
    result = -1;
  } else if ( result < 0 ) {
    nc_close( result );
    fprintf( stderr, "\nInvalid id (%d) for file '%s'.\n", result, name );
    result = -1;
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: openNetCDFFile - Open an existing NetCDF file for reading or writing.
INPUTS:  const char* const name  Name of NetCDF file to open.
         const char rw           'r' to read, 'w' to write.
RETURNS: int >= 0 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int openNetCDFFile( const char* const name, const char rw ) {
  PRE02( name, *name );
  int result = -1;
  const int status =
    nc_open( name, rw == 'r' ? NC_NOWRITE | NC_SHARE : NC_WRITE | NC_SHARE, &result );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't open file '%s' because %s.\n", name, message );
    result = -1;
  } else if ( result < 0 ) {
    nc_close( result );
    fprintf( stderr, "\nInvalid id (%d) for file '%s'.\n", result, name );
    result = -1;
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: closeNetCDFFile - Close a NetCDF file.
INPUTS:  const int id  Id of opened NetCDF file to close.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int closeNetCDFFile( const int id ) {
  PRE0( id >= 0 );
  int status = nc_sync( id );
  int result = status == NC_NOERR;

  if ( result ) {
    status = nc_close( id );
    result = status == NC_NOERR;
  }

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't flush/close file because %s.\n", message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: flushNetCDFFile - Flush write buffers to a NetCDF file.
INPUTS:  const int id  Id of opened NetCDF file to flush.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int flushNetCDFFile( const int id ) {
  PRE0( id >= 0 );
  const int status = nc_sync( id );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't flush file because %s.\n", message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: endNetCDFHeader - End NetCDF file header definitions.
INPUTS:  const int id  Id of writable NetCDF file.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int endNetCDFHeader( const int id ) {
  PRE0( id >= 0 );
  const int status = nc_enddef( id );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't end definition because %s.\n", message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFDimension - Get size of named NetCDF dimension.
INPUTS:  const int file  Id of NetCDF file.
         const char* const name  Name of dimension.
RETURNS: int size if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int getNetCDFDimension( const int file, const char* const name ) {
  PRE0( file >= 0 );
  int id = -1;
  int status = nc_inq_dimid( file, name, &id );
  int result = 0;

  if ( status == NC_NOERR ) {
    size_t size = 0;
    status = nc_inq_dimlen( file, id, &size );

    if ( status == NC_NOERR ) {

      if ( size > INT_MAX ) {
        fprintf( stderr, "\nInvalid size (%lu) for dimension '%s'.\n",
                 size, name );
      } else {
        result = size;
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't get dimension '%s' size because %s.\n",
             name, message );
    result = 0;
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: checkNetCDFVariableId - Check if named NetCDF variable if it exists.
INPUTS:  const int file  Id of NetCDF file.
         const char* const name  Name of variable.
RETURNS: int id if ok, else -1.
******************************************************************************/

int checkNetCDFVariableId( const int file, const char* const name ) {
  PRE0( file >= 0 );
  int result = -1;
  int status = nc_inq_varid( file, name, &result );

  DEBUG( fprintf( stderr, "%d = nc_inq_varid( file %d, %s, result = %d )\n",
         status, file, name, result ); )

  if ( status == NC_NOERR ) {

    if ( result < 0 ) {
      fprintf( stderr, "\nInvalid file (id = %d) id (%d) for variable '%s'.\n",
               file, result, name );
      result = -1;
    }
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFVariableId - Get id of named NetCDF variable.
INPUTS:  const int file  Id of NetCDF file.
         const char* const name  Name of variable.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int getNetCDFVariableId( const int file, const char* const name ) {
  PRE0( file >= 0 );
  int result = -1;
  int status = nc_inq_varid( file, name, &result );

  DEBUG( fprintf( stderr, "%d = nc_inq_varid( file %d, %s, result = %d )\n",
         status, file, name, result ); )

  if ( status == NC_NOERR ) {

    if ( result < 0 ) {
      fprintf( stderr, "\nInvalid file (id = %d) id (%d) for variable '%s'.\n",
               file, result, name );
      result = -1;
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr,
             "\nCan't get file (id = %d) variable '%s' id because %s.\n",
             file, name, message );
    result = -1;
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: isNetCDFFloat - Is NetCDF type float?
INPUTS:  const int type  NetCDF type.
RETURNS: int 1 if float else 0.
******************************************************************************/

int isNetCDFFloat( const int type ) {
  const int result = type == NC_FLOAT;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: getM3IOVariableDimensions - Get dimension of M3IO variable.
INPUTS:  const int file  NetCDF file id.
         const char* const variable  Name of variable.
OUTPUTS: int dims[ 4 ]  dims[ COLUMN, ROW, LAYER, TIME ].
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getM3IOVariableDimensions( const int file,
                               const char* const variable,
                               int dims[ 4 ] ) {

  PRE04( file >= 0, variable, *variable, dims );

  int result = 0;
  const int id = getNetCDFVariableId( file, variable );
  memset( dims, 0, 4 * sizeof (int) );

  if ( id >= 0 ) {
    nc_type type = 0;
    int rank = 0;
    int dimids[ NC_MAX_DIMS ];
    int status = nc_inq_var( file, id, 0, &type, &rank, dimids, 0 );

    if ( status == NC_NOERR ) {

      if ( type != NC_FLOAT ) {
        fprintf( stderr, "\nInvalid type (%d) of variable '%s'.\n",
                 (int) type, variable );
      } else if ( rank != 4 ) {
        fprintf( stderr, "\nInvalid rank (%d) of variable '%s'.\n",
                 rank, variable );
      } else {
        size_t size = 0;
        status = nc_inq_dimlen( file, dimids[ 0 ], &size );

        if ( status == NC_NOERR ) {
          dims[ TIME ] = size;
          status = nc_inq_dimlen( file, dimids[ 1 ], &size );

          if ( status == NC_NOERR ) {
            dims[ LAYER ] = size;
            status = nc_inq_dimlen( file, dimids[ 2 ], &size );

            if ( status == NC_NOERR ) {
              dims[ ROW ] = size;
              status = nc_inq_dimlen( file, dimids[ 3 ], &size );

              if ( status == NC_NOERR ) {
                dims[ COLUMN ] = size;

                if ( ! GT_ZERO4( dims[ 0 ], dims[ 1 ], dims[ 2 ], dims[ 3 ])) {
                  fprintf( stderr, "\nInvalid dimensions (%d, %d, %d, %d) "
                           "of variable '%s'.\n",
                           dims[ 0 ], dims[ 1 ], dims[ 2 ], dims[ 3 ],
                           variable );
                } else {
                  result = 1;
                }
              }
            }
          }
        }
      }
    }

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't get variable '%s' id because %s.\n",
               variable, message );
    }
  }

  if ( ! result ) {
    memset( dims, 0, 4 * sizeof (int) );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        GT_ZERO4( dims[ 0 ], dims[ 1 ], dims[ 2 ], dims[ 3 ] ),
                        IS_ZERO4( dims[ 0 ], dims[ 1 ], dims[ 2 ], dims[ 3])));
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFVariableInfo - Get metadata info of NetCDF variable.
INPUTS:  const int file        NetCDF file id.
         const int id          Id of variable.
OUTPUTS: char name[ 128 ]      0 or Name of variable.
         int* type             Type of variable.
         int* rank             Number of dimensions of variable.
         int dimensions[ 32 ]  Dimensions of variable.
         char units[ 128 ]     0 or Units of variable.
         char description[ 128 ] 0 or description of variable.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getNetCDFVariableInfo( const int file, const int id,
                           char name[], int* type, int* rank,
                           int dimensions[], char units[],
                           char description[] ) {

  PRE05( file >= 0, id >= 0, type, rank, dimensions );

  nc_type ntype = 0;
  char name0[ 128 ] = "";
  int dimids[ NC_MAX_DIMS ];
  int status = nc_inq_var( file, id, name0, &ntype, rank, dimids, 0 );
  int result = AND2( status == NC_NOERR, *name0 );
  int dimension = 0;
  *type = (int) ntype;

  for ( ; AND2( result, dimension < *rank ); ++dimension ) {
    size_t size = 0;
    status = nc_inq_dimlen( file, dimids[ dimension ], &size );
    result = AND2( status == NC_NOERR, size > 0 );
    dimensions[ dimension ] = size;

    if ( size == 0 ) {
      fprintf( stderr,
               "\nInvalid dimension #%d of variable %s.\n", dimension, name0 );
    }
  }

  if ( AND2( result, units ) ) {
    result = getNetCDFStringAttribute( file, id, "units", 128, units );
  }

  if ( AND2( result, description ) ) {
    result = getNetCDFStringAttribute( file, id, "var_desc", 128, description );
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't get variable %s info because %s.\n",
             name0, message );
  }

  if ( name ) {
    strncpy( name, name0, sizeof name0 );
  }

  /* Ensure non-empty units output: */

  if ( AND2( units, IN3( *units, '\0', ' ' ) ) ) {
    units[ 0 ] = '-';
    units[ 1 ] = '\0';

  }

  /* Ensure non-empty description output: */

  if ( AND2( description, IN3( *description, '\0', ' ' ) ) ) {
    snprintf( description, sizeof name0 / sizeof *name0,
              "Variable %s", name0 );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND3( GT_ZERO3( *rank, dimensions[ 0 ],
                                   dimensions[ *rank - 1 ] ),
                         OR2( description == 0,
                              AND2( *description, strlen( description) < 128)),
                         OR2( name == 0,
                              AND2( *name, strlen( name ) < 128 ) ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFStringAttribute - Get named attribute string value.
INPUTS:  const int file          NetCDF file id.
         const int id            Id  of variable or -1 if global.
         const char* const name  Name of attribute.
         const size_t size       Size of value[].
OUTPUTS: char value[ size ]      Value of attribute.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getNetCDFStringAttribute( const int file, const int id,
                              const char* const name,
                              const size_t size, char value[] ) {

  PRE06( file >= 0, id >= -1, name, *name, size, value );

  int result = 0;
  size_t length = 0;
  int status = nc_inq_attlen( file, id, name, &length );
  memset( value, 0, size );

  if ( AND2( status == NC_NOERR, IN_RANGE( length, 1, size - 1 ) ) ) {
    status = nc_get_att_text( file, id, name, value );
    result = status == NC_NOERR;
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't get attribute '%s' value because %s.\n",
             name, message );
  } else { /* Trim trailing spaces: */
    size_t index = strlen( value );
    int done = 0;

    while ( AND2( index--, done == 0 ) ) {

      if ( value[ index ] == ' ' ) {
        value[ index ] = '\0';
      } else {
        done = 1;
      }
    }
  }

  POST02( IS_BOOL( result ), value[ size - 1 ] == '\0' );
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFDoubleAttribute - Get named attribute double value.
INPUTS:  const int file          NetCDF file id.
         const char* const name  Name of attribute.
OUTPUTS: double* const value     Value of attribute.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getNetCDFDoubleAttribute( const int file, const char* const name,
                              double* const value ) {

  PRE04( file >= 0, name, *name, value );

  int result = 0;
  int status = 0;
  *value = 0.0;
  status = nc_get_att_double( file, NC_GLOBAL, name, value );
  result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't get attribute '%s' value because %s.\n",
             name, message );
  } else if ( isNaN( *value ) ) {
    *value = BADVAL3;
  } else if ( *value < -FLT_MAX ) {
    *value = -FLT_MAX;
  } else if ( *value > FLT_MAX ) {
    *value = FLT_MAX;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFFloatAttribute - Get named attribute float value.
INPUTS:  const int file          NetCDF file id.
         const char* const name  Name of attribute.
OUTPUTS: float* const value     Value of attribute.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getNetCDFFloatAttribute( const int file, const char* const name,
                             float* const value ) {

  PRE04( file >= 0, name, *name, value );

  int result = 0;
  int status = 0;
  *value = 0.0;
  status = nc_get_att_float( file, NC_GLOBAL, name, value );
  result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't get attribute '%s' value because %s.\n",
             name, message );
  } else if ( isNaN( *value ) ) {
    *value = BADVAL3;
  } else if ( *value < -FLT_MAX ) {
    *value = -FLT_MAX;
  } else if ( *value > FLT_MAX ) {
    *value = FLT_MAX;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFIntAttribute - Get named attribute int value.
INPUTS:  const int file          NetCDF file id.
         const char* const name  Name of attribute.
OUTPUTS: int* const value        Value of attribute.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getNetCDFIntAttribute( const int file, const char* const name,
                           int* const value ) {

  PRE04( file >= 0, name, *name, value );

  int result = 0;
  int status = 0;
  *value = 0;
  status = nc_get_att_int( file, NC_GLOBAL, name, value );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't get attribute '%s' value because %s.\n",
             name, message );
  } else {
    result = 1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: getNetCDFDoubleAttribute - Get named attribute double value.
INPUTS:  const int file          NetCDF file id.
         const char* const name  Name of attribute.
         const int count         Number of values.
OUTPUTS: float values[ count ]   Values of attribute.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getNetCDFFloatArrayAttribute( const int file, const char* const name,
                                  const int count, float values[] ) {

  PRE05( file >= 0, name, *name, count > 0, values );

  int result = 0;
  int status = 0;
  nc_type type = 0;
  size_t length = 0;
  memset( values, 0, count * sizeof *values );
  status = nc_inq_att( file, NC_GLOBAL, name, &type, &length );
  result = AND3( status == NC_NOERR, type == NC_FLOAT, length == count );

  if ( result ) {
    status = nc_get_att_float( file, NC_GLOBAL, name, values );
    result = status == NC_NOERR;
  }

  if ( result ) {
    int index = 0;

    for ( ; AND2( result, index < count ); ++index ) {
      const float value = values[ index ];
      result = IN_RANGE( value, -FLT_MAX, FLT_MAX );
    }
  }

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't get attribute '%s' values because %s.\n",
             name, message );
    memset( values, 0, count * sizeof *values );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: getM3IOFileTimeRange - Get first/last timestamp of file data.
INPUTS:  const int file                 NetCDF file id.
OUTPUTS: int* const yyyymmddhh1  First timestamp.
         int* const yyyymmddhh2  Last  timestamp.
         int* timesteps          Number of timesteps.
         int* hoursPerTimestep   Number of hours per timestep.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int getM3IOFileTimeRange( const int file,
                          int* const yyyymmddhh1,
                          int* const yyyymmddhh2,
                          int* timesteps,
                          int* hoursPerTimestep ) {

  PRE05( file >= 0, yyyymmddhh1, yyyymmddhh2, timesteps, hoursPerTimestep );

  int result = 0;
  *timesteps = getNetCDFDimension( file, "TSTEP" );

  if ( *timesteps == 0 ) {
    *timesteps = 1;
  }

  result = *timesteps > 0;

  if ( result ) {
    int yyyyddd = 0;
    result = getNetCDFIntAttribute( file, "SDATE", &yyyyddd );

    if ( AND2( result, isValidYYYYDDD( yyyyddd ) ) ) {
      int yyyymmdd = toYYYYMMDD( yyyyddd );
      int hhmmss = 0;
      result = getNetCDFIntAttribute( file, "STIME", &hhmmss );

      if ( AND2( result, isValidHHMMSS( hhmmss ) ) ) {
        int hh = hhmmss / 10000;
        *yyyymmddhh2 = *yyyymmddhh1 = yyyymmdd * 100 + hh;
        hhmmss = 0;
        result = getNetCDFIntAttribute( file, "TSTEP", &hhmmss );

        if ( result ) {
          hh = hhmmss / 10000;

          if ( hh <= 0 ) {
            hh = 1;
          }

          *hoursPerTimestep = hh;
          hh *= *timesteps;

          if ( hh > 0 ) {
            *yyyymmddhh2 = incrementHours( *yyyymmddhh1, hh - 1 );
          }
        }
      }
    }
  }

  if ( ! result ) {
    fprintf(stderr, "\nFailed to read valid file time dimension/attributes.\n");
    *yyyymmddhh1 = *yyyymmddhh2 = *timesteps = *hoursPerTimestep = 0;
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                         AND5( isValidYYYYMMDDHH( *yyyymmddhh1 ),
                               isValidYYYYMMDDHH( *yyyymmddhh2 ),
                               *yyyymmddhh1 <= *yyyymmddhh2,
                               *timesteps > 0,
                               *hoursPerTimestep > 0 ),
                         IS_ZERO4( *yyyymmddhh1, *yyyymmddhh2,
                                   *timesteps, *hoursPerTimestep ) ) );
  return result;
}



/******************************************************************************
PURPOSE: readM3IOVariable - Read a subset of variable data.
INPUTS:  const int file                 NetCDF file id.
OUTPUTS: int* const yyyymmddhh1  First timestamp.
         int* const yyyymmddhh2  Last  timestamp.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

int readM3IOVariable( const int file,
                      const int id,
                      const int time0,   const int time1,
                      const int layer0,  const int layer1,
                      const int row0,    const int row1,
                      const int column0, const int column1,
                      float array[] ) {

  PRE06( GE_ZERO10( file, id, time0, time1, layer0, layer1, row0, row1,
                    column0, column1 ),
         time0 <= time1, layer0 <= layer1, row0 <= row1, column0 <= column1,
         array );

  int result = 0;
  int status = 0;
  size_t starts[ 4 ] = { 0, 0, 0, 0 };
  size_t counts[ 4 ] = { 0, 0, 0, 0 };
  starts[ 0 ] = time0;
  starts[ 1 ] = layer0;
  starts[ 2 ] = row0;
  starts[ 3 ] = column0;
  counts[ 0 ] = COUNT_IN_RANGE( time0, time1 );
  counts[ 1 ] = COUNT_IN_RANGE( layer0, layer1 );
  counts[ 2 ] = COUNT_IN_RANGE( row0, row1 );
  counts[ 3 ] = COUNT_IN_RANGE( column0, column1 );
  status = nc_get_vara_float( file, id, starts, counts, array );

  DEBUG( fprintf( stderr,
                  "starts = [%lu %lu %lu %lu], counts = [%lu %lu %lu %lu]\n",
                  starts[0], starts[1], starts[2], starts[3],
                  counts[0], counts[1], counts[2], counts[3] ); )

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't read subset of variable data because %s.\n",
             message );
  } else { /* Filter-out invalid values: */
    const size_t count = counts[ 0 ] * counts[ 1 ] * counts[ 2 ] * counts[ 3 ];
    size_t index = 0;

    for ( ; index < count; ++index ) {
      const float value = array[ index ];

      if ( ! IS_VALID_VALUE( value ) ) {
        array[ index ] = BADVAL3;
      }
    }

    result = 1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFDimension - Create a dimension.
INPUTS:  const int file          NetCDF file id.
         const char* const name  Name of dimension.
         const int size          Size of dimension or 0 if UNLIMITED.
OUTPUTS: int* const id           Id of dimension.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int createNetCDFDimension( const int file, const char* const name,
                           const int size, int* const id ) {

  PRE05( file >= 0, name, *name, size >= 0, id );

  const int status = nc_def_dim( file, name, size, id );
  const int result = AND2( status == NC_NOERR, *id >= 0 );

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't create dimension '%s' of size %d because %s.\n",
             name, size, message );
    *id = -1;
  }

  POST02( IS_BOOL( result ), IMPLIES_ELSE( result, *id >= 0, *id == -1 ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFVariable - Create a variable.
INPUTS:  const int file            NetCDF file id.
         const char* const name    Name of variable.
         const char* const units   Units of variable.
         const char* const attributeName   0 or optional attribute name.
         const char* const attributeVale   0 or optional attribute value.
         const int rank            Number of dimensions of variable.
         const int dimids[ rank ]  Dimension ids of variable (left to right).
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int createNetCDFVariable( const int file,
                          const char* const name,
                          const char* const units,
                          const char* const attributeName,
                          const char* const attributeValue,
                          const int rank,
                          const int dimids[] ) {

  PRE08( file >= 0, name, *name, units, *units,
         IMPLIES_ELSE( attributeName,
                       AND3( *attributeName, attributeValue, *attributeValue),
                       attributeValue == 0 ),
         IN_RANGE( rank, 1, 32 ),
         dimids );

  int result = -1;
  const int isTFLAG = ! strcmp( name, "TFLAG" );
  const int isInt =
    OR3( isTFLAG,
         ! strcmp( name, "yyyymmdd" ),
         ! strcmp( name, "hhmmss" ) );
  int status =
    nc_def_var( file, name, isInt ? NC_INT : NC_FLOAT, rank, dimids, &result );

  if ( AND2( status == NC_NOERR, result >= 0 ) ) {
    char long_name[ 17 ] = ""; /* Detect and add M3IO long_name attribute. */
    memset( long_name, 0, sizeof long_name );

    if ( AND2( attributeName, ! strcmp( attributeName, "var_desc" ) ) ) {
      snprintf(long_name, sizeof long_name / sizeof *long_name, "%-16s", name);
    }

    if ( AND2( long_name[ 0 ], ! isTFLAG ) ) {
      status = nc_put_att_text( file, result, "long_name", 16, long_name );
    }

    if ( status == NC_NOERR ) {
      status = nc_put_att_text( file, result, "units", strlen( units ), units );

      if ( status == NC_NOERR ) {

        if ( AND2( long_name[ 0 ], isTFLAG ) ) {
          status = nc_put_att_text( file, result, "long_name", 16, long_name );
        }

        if ( attributeName ) {
          status = nc_put_att_text( file, result, attributeName,
                                    strlen( attributeValue ), attributeValue );
        }
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't create variable '%s' because %s.\n",
             name, message );
    result = -1;
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: copyNetCDFAttribute - Copy a global string attribute.
INPUTS:  const int file           NetCDF file id.
         const char* const name   Name of attribute to copy.
OUTPUTS: const int output         Output NetCDF file id.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int copyNetCDFAttribute( const int file,
                         const char* const name,
                         const int output ) {

  PRE04( file >= 0, name, *name, output >= 0 );

  const int status = nc_copy_att( file, -1, name, output, -1 );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't create variable '%s' because %s.\n",
             name, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFStringAttribute - Create a string attribute.
INPUTS:  const int file           NetCDF file id.
         const int id             NetCDF id or -1 if global.
         const char* const name   Name of attribute
         const char* const value  Value of attribute.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int createNetCDFStringAttribute( const int file,
                                 const int id,
                                 const char* const name,
                                 const char* const value ) {

  PRE06( file >= 0, id >= -1, name, *name, value, *value );

  const int status = nc_put_att_text( file, id, name, strlen( value ), value );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr,
             "\nCan't create attribute '%s' with value '%s' because %s.\n",
             name, value, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFIntAttribute - Create an int global attribute.
INPUTS:  const int file           NetCDF file id.
         const char* const name   Name of attribute
         const int value          Value of attribute.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int createNetCDFIntAttribute( const int file,
                              const char* const name,
                              const int value ) {

  PRE03( file >= 0, name, *name );

  const int status = nc_put_att_int( file, NC_GLOBAL, name, NC_INT, 1, &value );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr,
             "\nCan't create attribute '%s' with value %d because %s.\n",
             name, value, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFDoubleAttribute - Create a double global attribute.
INPUTS:  const int file           NetCDF file id.
         const char* const name   Name of attribute
         const double value       Value of attribute.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int createNetCDFDoubleAttribute( const int file,
                                 const char* const name,
                                 const double value ) {

  PRE03( file >= 0, name, *name );

  const int status =
    nc_put_att_double( file, NC_GLOBAL, name, NC_DOUBLE, 1, &value );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr,
             "\nCan't create attribute '%s' with value %g because %s.\n",
             name, value, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFFloatAttribute - Create a float global attribute.
INPUTS:  const int file           NetCDF file id.
         const char* const name   Name of attribute
         const float value          Value of attribute.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int createNetCDFFloatAttribute( const int file,
                                const char* const name,
                                const float value ) {

  PRE03( file >= 0, name, *name );

  const int status =
    nc_put_att_float( file, NC_GLOBAL, name, NC_FLOAT, 1, &value );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr,
             "\nCan't create attribute '%s' with value %g because %s.\n",
             name, value, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createNetCDFFloatArrayAttribute - Create a float array global
         attribute.
INPUTS:  const int file           NetCDF file id.
         const char* const name   Name of attribute
         const int value          Value of attribute.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/


int createNetCDFFloatArrayAttribute( const int file,
                                     const char* const name,
                                     const int count,
                                     const float values[] ) {

  PRE05( file >= 0, name, *name, count > 0, values );

  const int status =
    nc_put_att_float( file, NC_GLOBAL, name, NC_FLOAT, count, values );
  const int result = status == NC_NOERR;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "\nCan't create attribute '%s' because %s.\n",
             name, message );
  }

  POST0( IS_BOOL( result ) );
  return result;

}



/******************************************************************************
PURPOSE: writeCOARDS2DVariable - Write 2D variable data.
INPUTS:  const int file            NetCDF file id.
         const char* const name    Name of variable.
         const int rows            Number of rows.
         const int columns         Number of columns.
         const float array[ rows * columns ]  Subset data to write.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int writeCOARDS2DVariable( const int file,
                           const char* const name,
                           const int rows,
                           const int columns,
                           const float array[] ) {

  PRE05( file >= 0, name, *name, GT_ZERO2( rows, columns ), array );

  int result = 0;
  int id = getNetCDFVariableId( file, name );

  if ( id >= 0 ) {
    int status = 0;
    size_t starts[ 3 ] = { 0, 0 };
    size_t counts[ 3 ] = { 0, 0 };
    counts[ 0 ] = rows;
    counts[ 1 ] = columns;
    status = nc_put_vara_float( file, id, starts, counts, array );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr,
               "\nCan't write variable '%s' because %s.\n",
               name, message );
    } else {
      result = 1;
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSTimeVariables - Write COARDS standard time data for
         timestep.
INPUTS:  const int file            NetCDF file id.
         const int timestep        0-based timestep.
         const int yyyymmddhh      Timestamp of timestep.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int writeCOARDSTimeVariables( const int file,
                              const int timestep,
                              const int yyyymmddhh ) {

  PRE03( file >= 0, timestep >= 0, isValidYYYYMMDDHH( yyyymmddhh ) );

  int result = 0;
  int id = getNetCDFVariableId( file, "time" );

  if ( id >= 0 ) {
    const size_t start = timestep;
    const size_t count = 1;
    const float value = timestep;
    int status = nc_put_vara_float( file, id, &start, &count, &value );

    if ( status == NC_NOERR) {
      id = getNetCDFVariableId( file, "yyyymmdd" );

      if ( id >= 0 ) {
        const int yyyymmdd = yyyymmddhh / 100;
        status = nc_put_vara_int( file, id, &start, &count, &yyyymmdd );

        if ( status == NC_NOERR ) {
          id = getNetCDFVariableId( file, "hhmmss" );

          if ( id >= 0 ) {
            const int hhmmss = yyyymmddhh % 100 * 10000;
            status = nc_put_vara_int( file, id, &start, &count, &hhmmss );
          }
        }
      }
    }

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "\nCan't write time variables because %s.\n", message );
    } else {
      result = 1;
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeM3IOVariable - Write a single timestep of a variable's data.
INPUTS:  const int file            NetCDF file id.
         const char* const name    Name of variable.
         const int timestep        0-based timestep index.
         const int layers          Number of layers.
         const int rows            Number of rows.
         const int columns         Number of columns.
         const float array[ layers ][ rows ][ columns ]  Subset data to write.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int writeM3IOVariable( const int file,
                       const char* const name,
                       const int timestep,
                       const int layers,
                       const int rows,
                       const int columns,
                       const float array[] ) {

  PRE06( file >= 0, name, *name, timestep >= 0,
         GT_ZERO3( layers, rows, columns ), array );

  int result = 0;
  int id = getNetCDFVariableId( file, name );

  if ( id >= 0 ) {
    int status = 0;
    size_t starts[ 4 ] = { 0, 0, 0, 0 };
    size_t counts[ 4 ] = { 1, 0, 0, 0 };
    starts[ 0 ] = timestep;
    counts[ 1 ] = layers;
    counts[ 2 ] = rows;
    counts[ 3 ] = columns;
    status = nc_put_vara_float( file, id, starts, counts, array );

    DEBUG( fprintf( stderr,
                   "writeM3IOVariable( file = %d, name = %s, timestep = %d, "
                   "layers = %d, rows = %d, columns = %d, "
                   "array = [%f ... %f]\n",
                   file, name, timestep, layers, rows, columns,
                   array[0], array[ layers * rows * columns - 1 ] ); )

    DEBUG( fprintf( stderr,
                    "starts = [%lu, %lu, %lu, %lu], "
                    "counts = [%lu, %lu, %lu, %lu]\n",
                    starts[0], starts[1], starts[2], starts[3],
                    counts[0], counts[1], counts[2], counts[3] ); )

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr,
               "\nCan't write timestep %d of variable '%s' because %s.\n",
               timestep, name, message );
    } else {
      result = 1;
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeTFLAGVariable - Write a single timestep of a TFLAG data.
INPUTS:  const int file            NetCDF file id.
         const int timestep        0-based timestep index.
         const int variable        Number of variables.
         const int datetime        DATE-TIME dimension.
         const int array[ layers ][ rows ][ columns ]  Subset data to write.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int writeTFLAGVariable( const int file,
                        const int timestep,
                        const int variables,
                        const int datetime,
                        const int array[] ) {

  PRE05( file >= 0, timestep >= 0, variables > 0, datetime > 0, array );

  int result = 0;
  const char* const name = "TFLAG";
  int id = getNetCDFVariableId( file, name );

  if ( id >= 0 ) {
    int status = 0;
    size_t starts[ 3 ] = { 0, 0, 0 };
    size_t counts[ 3 ] = { 1, 0, 0 };
    starts[ 0 ] = timestep;
    counts[ 1 ] = variables;
    counts[ 2 ] = datetime;
    status = nc_put_vara_int( file, id, starts, counts, array );

    DEBUG( fprintf( stderr,
                   "writeTFLAGVariable( file = %d, timestep = %d, "
                   "variables = %d, datetime = %d, "
                   "array = [%d ... %d]\n",
                   file, timestep, variables, datetime,
                   array[0], array[ variables * datetime - 1 ] ); )

    DEBUG( fprintf( stderr,
                    "starts = [%lu, %lu, %lu], "
                    "counts = [%lu, %lu, %lu]\n",
                    starts[0], starts[1], starts[2],
                    counts[0], counts[1], counts[2] ); )

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr,
               "\nCan't write timestep %d of variable '%s' because %s.\n",
               timestep, name, message );
    } else {
      result = 1;
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



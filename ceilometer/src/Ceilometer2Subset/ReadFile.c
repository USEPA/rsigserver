/******************************************************************************
PURPOSE: ReadFile.c - Simple to use wrapper routines to read data from
                      ceilometer HDF5 files.

NOTES:   Uses NetCDF4 library and libs it depends on (HDF5, curl, z, dl).

HISTORY: 2022-07-12 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h>  /* For assert(). */
#include <stdio.h>   /* For stderr, fprintf(). */
#include <stdlib.h>  /* For atof(). */
#include <string.h>  /* For memset(). */

#include <netcdf.h> /* For nc_*(). */

#include "ReadFile.h"  /* For public interface. */

/*=========================== FORWARD DECLARATIONS ==========================*/

static int matchedVariableDimensions( const int file,
                                      const int id,
                                      const size_t dimension0,
                                      const size_t dimension1 );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: openFile - Open NetCDF file for reading.
INPUTS:  const char* const fileName   Name of NetCDF file to open for reading.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int openFile( const char* const fileName ) {
  int result = -1;
  int status = -1;
  assert( fileName ); assert( *fileName );
  status = nc_open( fileName, NC_NOWRITE, &result );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to open NetCDF file for reading because: %s\n",
             message );
    result = -1;
  }

  assert( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: closeFile - Close NetCDF file.
INPUTS:  const int file  NetCDF file to close.
******************************************************************************/

void closeFile( const int file ) {
  assert( file > -1 );
  nc_close( file );
}



/******************************************************************************
PURPOSE: readVariableId - Get id of named variable.
INPUTS:  const int file               File to read.
         const char* const variable   Name of variable.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int readVariableId( const int file, const char* const variable ) {
  int result = -1;
  int status = -1;
  assert( file > -1 ); assert( variable ); assert( *variable );
  status = nc_inq_varid( file, variable, &result );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to get id of variable '%s' because: %s\n",
             variable, message );
    result = -1;
  }

  assert( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: readFileAttribute - Read global attribute.
INPUTS:  const int file               Id of file to read.
         const char* const attribute  Name of global attribute.
OUTPUTS: double* value                Value of attribute read.
RETURNS: int 1 if successful, else 0 and failure message on stderr.
******************************************************************************/

int readFileAttribute( const int file, const char* const attribute,
                       double* value ) {
  int result = 0;
  int status = -1;
  int id = -1;

  assert( file > -1 ); assert( attribute ); assert( *attribute );
  assert( value );

  status = nc_inq_attid( file, NC_GLOBAL, attribute, &id );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to get id of attribute '%s' because: %s\n",
             attribute, message );
  } else {
    nc_type type = 0;
    status = nc_inq_atttype( file, NC_GLOBAL, attribute, &type );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "Failed to type of attribute '%s' because: %s\n",
               attribute, message );
    } else {

      if ( type == NC_DOUBLE ) {
        status = nc_get_att_double( file, NC_GLOBAL, attribute, value );

        if ( status != NC_NOERR ) {
          const char* const message = nc_strerror( status );
          fprintf( stderr, "Failed to read double attribute '%s' because: %s\n",
                  attribute, message );
        } else {
          result = 1;
        }
      } else if ( type == NC_CHAR ) { /* UGLY */
        char svalue[ NC_MAX_NAME + 1 ] = "";
        memset( svalue, 0, sizeof svalue );
        status = nc_get_att_text( file, NC_GLOBAL, attribute, svalue );

        if ( status != NC_NOERR ) {
          const char* const message = nc_strerror( status );
          fprintf( stderr, "Failed to read string attribute '%s' because: %s\n",
                  attribute, message );
        } else {
          *value = atof( svalue );
          result = 1;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readVariableDimensions - Read variable dimensions.
INPUTS:  const int file      File to query.
         const int id        Id of variable to query.
OUTPUTS: size_t* dimension0  1st dimension of variable.
         size_t* dimension1  2nd dimension of variable or 1 if rank 1.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readVariableDimensions( const int file,
                            const int id,
                            size_t* const dimension0,
                            size_t* const dimension1 ) {
  int result = 0;

  assert( file > -1 ); assert( id > -1 );
  assert( dimension0 ); assert( dimension1 );

  *dimension0 = *dimension1 = 0;

  {
    int status = 0;
    int rank = 0;
    int dimensionIds[ NC_MAX_DIMS ];
    memset( dimensionIds, 0, sizeof dimensionIds );
    status = nc_inq_var( file, id, 0, 0, &rank, dimensionIds, 0 );

    if ( status == NC_NOERR ) {

      if ( rank == 1 ) {
        status = nc_inq_dimlen( file, dimensionIds[ 0 ], dimension0 );
        *dimension1 = 1;
      } else if ( rank == 2 ) {
        status = nc_inq_dimlen( file, dimensionIds[ 0 ], dimension0 );

        if ( status == NC_NOERR ) {
          status = nc_inq_dimlen( file, dimensionIds[ 1 ], dimension1 );
        }
      } else {
        fprintf( stderr, "Failed to read valid rank: %d\n", rank );
      }
    }

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "Failed to read valid dimensions because: %s\n", message);
    } else if ( *dimension0 < 1 ) {
      fprintf(stderr, "Failed to read valid 1st dimension: %lu\n", *dimension0);
    } else if ( rank == 2 && *dimension1 < 1 ) {
      fprintf(stderr, "Failed to read valid 2nd dimension: %lu\n", *dimension1);
    } else {
      result = 1;
    }
  }

  if ( ! result ) {
    *dimension0 = *dimension1 = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file data.
INPUTS:  const int file            ID of file to read.
         const int id              ID of varaible to read.
         const size_t dimension0   1st dimension of data to read.
         const size_t dimension1   2nd dimension of data to read or 0.
OUTPUTS: double data[ dimension0 * dimension1 ]  Data read.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileData( const int file,
                  const int id,
                  const size_t dimension0,
                  const size_t dimension1,
                  double data[] ) {

  int result = 0;

  assert( file > -1 ); assert( id > -1 ); assert( dimension0 ); assert( data );

  result = matchedVariableDimensions( file, id, dimension0, dimension1 );

  if ( result ) {
    nc_type type = 0;
    int status = nc_inq_vartype( file, id, &type );

    if ( status == NC_NOERR ) {
      const size_t starts[ 2 ] = { 0, 0 };
      size_t counts[ 2 ] = { 0, 0 };
      counts[ 0 ] = dimension0;
      counts[ 1 ] = dimension1;

      if ( type == NC_FLOAT ) {
        status = nc_get_vara_float( file, id, starts, counts, (float*) data );

        if ( status == NC_NOERR ) { /* Expand floats to doubles: */
          size_t count = dimension0 * ( dimension1 > 0 ? dimension1 : 1 );
          const float* input = (float*) data;
          double* output = data;
          input += count;
          output += count;

          while ( count-- ) {
            *--output = *--input;
          }

          result = 1;
        }
      } else if ( type == NC_DOUBLE ) {
        status = nc_get_vara_double( file, id, starts, counts, data );
        result = status == NC_NOERR;
      }
    }

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "Failed to read valid data because: %s\n", message );
    }
  }

  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: matchedVariableDimensions - Does the variable have the expected
         dimensions?
INPUTS:  const int file          ID of file to check.
         const int id            ID of variable to check.
         const int dimension0    1st dimension.
         const int dimension1    2nd dimension or 0 if rank 1.
RETURNS: int 1 if matched, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int matchedVariableDimensions( const int file,
                                      const int id,
                                      const size_t dimension0,
                                      const size_t dimension1 ) {

  int result = 0;
  size_t readDimension0 = 0;
  size_t readDimension1 = 0;

  assert( file > -1 ); assert( id > -1 );
  assert( dimension0 > 0 || dimension1 > 0 );

  result = readVariableDimensions( file, id, &readDimension0, &readDimension1 );

  if ( result ) {

    if ( readDimension0 != dimension0 ) {
      fprintf( stderr, "Invalid/mismatched dims[0] in dataset: "
               "%lu (expected %lu)\n", readDimension0, dimension0 );
    } else if ( readDimension1 != dimension1 ) {
      fprintf( stderr, "Invalid/mismatched dims[1] in dataset: "
               "%lu (expected %lu)\n", readDimension1, dimension1 );
    } else {
      result = 1;
    }
  }

  return result;
}



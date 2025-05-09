/******************************************************************************
PURPOSE: ReadData.c - Simple to use wrapper routines to read data from GOESBB
                      NetCDF files.

NOTES:   Uses NetCDF library.

HISTORY: 2018-11-12 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For memset(). */

#include <netcdf.h> /* For NC*, nc_*(). */

#include "Utilities.h" /* For LONGITUDE, Bounds, expand32BitReals(). */
#include "ReadData.h"  /* For public interface. */

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: openFile - Open NetCDF file for reading.
INPUTS:  const char* const fileName   Name of NetCDF file to open for reading.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int openFile( const char* const fileName ) {
  int result = 0;
  int status = 0;
  assert( fileName );
  status = nc_open( fileName, NC_NOWRITE, &result );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to open NetCDF file for reading because: %s\n",
             message);
    result = -1;
  }

  return result;
}



/******************************************************************************
PURPOSE: closeFile - Close NetCDF file.
INPUTS:  int file NetCDF file to close.
******************************************************************************/

void closeFile( int file ) {
  const int status = nc_close( file );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to close NetCDF file because: %s\n",
             message );
  }
}



/******************************************************************************
PURPOSE: readFileDimensions - Read file dimensions.
INPUTS:  const int file     NetCDF file ID to read.
OUTPUTS: size_t* timesteps  Timesteps.
         size_t* points     Points.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileDimensions( const int file,
                        size_t* const timesteps, size_t* const points ) {
  int result = 0;
  int id = 0;
  int status = nc_inq_dimid( file, "x", &id );
  assert( timesteps ); assert( points );
  *timesteps = *points = 0;

  if ( status == NC_NOERR ) {
    status = nc_inq_dimlen( file, id, timesteps );

    if ( status == NC_NOERR ) {
      status = nc_inq_dimid( file, "y", &id );

      if ( status == NC_NOERR ) {
        status = nc_inq_dimlen( file, id, points );
      }
    }
  }

  result = status == NC_NOERR && *timesteps != 0 && *points != 0;

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read valid dimensions because: %s\n", message);
  }

  if ( ! result ) {
    *timesteps = *points = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readVariableDimensions - Read variable dimensions.
INPUTS:  const int file              NetCDF file ID to read.
         const char* const variable  Name of variable to check.
OUTPUTS: size_t* timesteps           Timesteps or 0.
         size_t* points              Points.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readVariableDimensions( const int file,
                            const char* const variable,
                            size_t* const timesteps,
                            size_t* const points ) {
  int result = 0;
  int id = 0;
  int status = 0;

  assert( file >= 0 ); assert( variable ); assert( *variable );
  assert( timesteps ); assert( points );

  *timesteps = *points = 0;
  status = nc_inq_varid( file, variable, &id );

  if ( status == NC_NOERR ) {
    int rank = 0;
    int dimensionIds[ NC_MAX_DIMS ];
    memset( dimensionIds, 0, sizeof dimensionIds );
    status = nc_inq_var( file, id, 0, 0, &rank, dimensionIds, 0 );

    if ( status == NC_NOERR ) {

      if ( rank == 1 ) {
        status = nc_inq_dimlen( file, dimensionIds[ 0 ], points );
      } else if ( rank == 2 ) {
        status = nc_inq_dimlen( file, dimensionIds[ 0 ], points );

        if ( status == NC_NOERR ) {
          status = nc_inq_dimlen( file, dimensionIds[ 1 ], timesteps );
        }
      }
    }
  }

  result = status == NC_NOERR && *points != 0;

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read valid dimensions because: %s\n", message);
  }

  if ( ! result ) {
    *timesteps = *points = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file data.
INPUTS:  const int file                 NetCDF file ID to read.
         const char* const variable     Name of variable to read.
         const size_t timesteps         Timesteps of variable to read or 0.
         const size_t points            Points of variable to read.
OUTPUTS: char units[ 80 ]               Units of variable.
         double data[ points * timesteps ]  Data for variable.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileData( const int file, const char* const variable,
                  const size_t timesteps, const size_t points,
                  char units[ 80 ], double data[] ) {

  int result = 0;
  int id = -1;
  int status = 0;

  assert( file >= 0 ); assert( variable ); assert( *variable );
  assert( points != 0 );
  assert( units ); assert( data );

  status = nc_inq_varid( file, variable, &id );

  if ( status == NC_NOERR ) {
    float* fdata = (float*) data;
    size_t start[ 2 ] = { 0, 0 };
    size_t count[ 2 ] = { 0, 0 };
    count[ 0 ] = points;
    count[ 1 ] = timesteps;
    status = nc_get_vara_float( file, id, start, count, fdata );

    if ( status == NC_NOERR ) {
      const size_t size = timesteps ? timesteps * points : points;
      expand32BitReals( size, data );
      status = nc_get_att_text( file, id, "units", units );
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read variable %s data/units because: %s\n",
             variable, message );
  } else {
    result = 1;
  }

  return result;
}




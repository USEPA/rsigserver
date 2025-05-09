/******************************************************************************
PURPOSE: ReadData.c - Simple to use wrapper routines to read data from VIIRS
                      NetCDF4 files.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-10-14 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For strcmp(), strstr(). */
#include <stdlib.h> /* For atoi(). */
#include <ctype.h>  /* For isdigit(). */
#include <float.h>  /* For DBL_MAX. */

#include <netcdf.h> /* For NC*, nc_*(). */

#include "Utilities.h" /* For LONGITUDE, Bounds, expand32BitReals(). */
#include "ReadData.h"  /* For public interface. */

/*=========================== FORWARD DECLARATIONS ===========================*/


static int clampInvalidValues( const size_t count,
                               const double validMinimum,
                               const double validMaximum,
                               float data[] );

static int parseHighQualityValue( const int file, const int id );

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
    fprintf( stderr, "Failed to open NetCDF file %s for reading because: %s\n",
             fileName, message);
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
PURPOSE: swathInDomain - Do the file longitude-latitude coordinates overlap
         domain?
INPUTS:  const int file       NetCDF id of input file.
         const Bounds domain  Lon-lat bounds of domain of subset of interest.
RETURNS: int 1 if file lon-lat coordinates overlap domain, else 0.
******************************************************************************/

int swathInDomain( const int file, const Bounds domain ) {
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
  int status = NC_NOERR;
  int index = 0;
  int result = 0;
  assert( file >= 0 ); assert( domain );

  for ( index = 0; index < count && status == NC_NOERR; ++index ) {
    status =
      nc_get_att_float( file, NC_GLOBAL, attributes[ index ], values + index );
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "%s\n", message );
  } else {
    Bounds bounds = { { 0.0, 0.0 }, { 0.0, 0.0 } };
    bounds[ LONGITUDE ][ MINIMUM ] = bounds[ LONGITUDE ][ MAXIMUM] = values[0];
    bounds[ LATITUDE  ][ MINIMUM ] = bounds[ LATITUDE  ][ MAXIMUM] = values[4];

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
     * Some VIIRS files contain bogus (-999.3) longitude/latitude values!
     * This is indicated by the geospatial bounds having invalid range.
     * Rather than reject the file in such cases,
     * proceed and later clamp these invalid values to the nearest valid value
     * in case such invalid values only appear at the edge of the swath rather
     * than in the middle somewhere (which would make it impossible to correctly
     * compute the corner vertices via geometric dual).
     */

    if ( isValidBounds( (const double (*)[2]) bounds ) ) {
      result = boundsOverlap( domain, (const double (*)[2]) bounds );
    } else {
      result = 1; /* Proceed anyway in case clamped points intersect domain. */
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileDimensions - Read file swath row/column dimensions.
INPUTS:  const int file   NetCDF file ID to read.
OUTPUTS: size_t* rows     Rows in swath.
         size_t* columns  Columns in swath.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileDimensions( const int file,
                        size_t* const rows, size_t* const columns ) {
  int result = 0;
  int id = 0;
  int status = nc_inq_dimid( file, "Columns", &id );
  assert( rows ); assert( columns );
  *rows = *columns = 0;

  if ( status == NC_NOERR ) {
    status = nc_inq_dimlen( file, id, columns );

    if ( status == NC_NOERR ) {
      status = nc_inq_dimid( file, "Rows", &id );

      if ( status == NC_NOERR ) {
        status = nc_inq_dimlen( file, id, rows );
      }
    }
  }

  result = status == NC_NOERR && rows != 0 && columns != 0;

  if ( ! result ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read valid dimensions because: %s\n", message );
    *columns = *rows = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file swath data.
INPUTS:  const int file                 NetCDF file ID to read.
         const char* const variable     Name of variable to read.
         const int allowMediumQuality   Allow medium-quality data?
         const size_t rows              Rows of variable to read.
         const size_t columns           Columns of variable to read.
OUTPUTS: double data[ rows * columns ]  Swath data for variable.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
NOTES:   Data is filtered by QC flag - i.e., set to MISSING_VALUE if
         QC != high = (3 before 2018-02-13 or 0 afterwards).
******************************************************************************/

int readFileData( const int file, const char* const variable,
                  const int allowMediumQuality,
                  const size_t rows, const size_t columns,
                  double data[] ) {
  int result = 0;
  int id = 0;
  int status = nc_inq_varid( file, variable, &id );
  assert( file >= 0 ); assert( variable ); assert( *variable );
  assert( allowMediumQuality == 0 || allowMediumQuality == 1 );
  assert( rows != 0 ); assert( columns != 0 ); assert( data );

  if ( status == NC_NOERR ) {
    float* const fdata = (float*) data;
    size_t start[ 2 ] = { 0, 0 };
    size_t count[ 2 ] = { 0, 0 };
    count[ 0 ] = rows;
    count[ 1 ] = columns;
    status = nc_get_vara_float( file, id, start, count, fdata );

    if ( status == NC_NOERR ) {
      const size_t points = rows * columns;

      if ( ! strcmp( variable, "Longitude" ) ) {
        result = clampInvalidValues( points, -180.0, 180.0, fdata );
      } else if ( ! strcmp( variable, "Latitude" ) ) {
        result = clampInvalidValues( points, -90.0, 90.0, fdata );
      } else {
        const int isAOD = ! strcmp( variable, "AOD550" );
        const int isAngstromExponent =
          strstr( variable, "AngsExp" ) == variable;
        size_t validPoints = 0;
        unsigned char* const bdata =
          (isAOD || isAngstromExponent) ? (unsigned char*) (fdata + points) :0;
        int highQuality = 0;

        if ( bdata ) { /* Read quality flag: */
          unsigned char* const bdata = (unsigned char*) ( fdata + points );
          status = nc_inq_varid( file, "QCAll", &id );

          if ( status == NC_NOERR ) {
            status = nc_get_vara_uchar( file, id, start, count, bdata );

            if ( status == NC_NOERR ) {
              highQuality = parseHighQualityValue( file, id );
            }
          }
        }

        /* Filter out-of-range or non-high-quality data: */

        if ( status == NC_NOERR ) {
          const double validMinimum =
            isAOD ? -0.05
            : isAngstromExponent ? -1.0
            : -DBL_MAX;
          const double validMaximum =
            isAOD ? 5.0
            : isAngstromExponent ? 3.0
            : DBL_MAX;
          size_t point = 0;

          for ( point = 0; point < points; ++point ) {
            const float value = fdata[ point ];
            const int isValid = value >= validMinimum && value <= validMaximum;

            if ( ! isValid ) {
              fdata[ point ] = MISSING_VALUE;
            } else if ( bdata != 0 ) { /* Check quality flag: */
              const int qualityValue = bdata[ point ];

              if ( qualityValue == highQuality ) {
                ++validPoints;
              } else if ( ! allowMediumQuality ) {
                fdata[ point ] = MISSING_VALUE; /* Filter non-high-quality. */
              } else if ( qualityValue == 1 ) { /* 1 is Medium quality value.*/
                ++validPoints;
              } else {
                fdata[ point ] = MISSING_VALUE; /* Filter non-high-quality. */
              }
            } else {
              ++validPoints;
            }
          }
        }

        /*
        fprintf( stderr,
                 "%-10s points = %lu, validPoints = %lu, highQuality = %d, "
                 "allowMediumQuality = %d\n",
                 variable, points, validPoints, highQuality,
                 allowMediumQuality );
        */

        result = validPoints != 0; /* At least one data value must be valid. */
      }

      if ( result ) {
        expand32BitReals( points, data );
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "%s\n", message );
  }

  return result;
}



/*============================= PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: clampInvalidValues - Clamp invalid data values to the most recent
         valid value.
INPUTS:  const size_t count         Number of data values.
         const double validMinimum  Minimum valid value.
         const double validMaximum  Maximum valid value.
         float data[ count ]        Data values to clamp.
OUTPUTS: float data[ count ]        Data values clamped.
RETURNS: int 1 if at least one valid value remains else 0.
******************************************************************************/

static int clampInvalidValues( const size_t count,
                               const double validMinimum,
                               const double validMaximum,
                               float data[] ) {
  int result = 0;
  size_t index = 0;
  float validValue = MISSING_VALUE; /* Initialize to invalid. */

  assert( count ); assert( validMinimum <= validMaximum );
  assert( MISSING_VALUE < validMinimum );
  assert( data );

  /* Find first valid value: */

  while ( index < count && validValue != MISSING_VALUE ) {
    const float value = data[ index ];
    const int isValid = value >= validMinimum && value <= validMaximum;

    if ( isValid ) {
      validValue = value; /* Store first valid value. */
    }

    ++index;
  }

  /* Clamp remaining invalid values to nearest valid value: */

  while ( index < count ) {
    const float value = data[ index ];
    const int isValid = value >= validMinimum && value <= validMaximum;

    if ( isValid ) {
      validValue = value; /* Update to nearest valid value. */
    } else {
      data[ index ] = validValue; /* Overwrite invalid value to nearest valid.*/
    }

    ++index;
  }

  result = validValue != MISSING_VALUE; /* Was there at least one valid value?*/
  return result;
}



/******************************************************************************
PURPOSE: parseHighQuality - Determine the value used to represent high-quality
         data.
INPUTS:  const int file  NetCDF file ID to read.
         const int id    NetCDF variable ID (for QCAll) to read.
RETURNS: int value used to represent high-quality data retrieval.
NOTES:   high = 3 before 2018-02-13 or 0 afterwards.
******************************************************************************/

static int parseHighQualityValue( const int file, const int id ) {
  const char* const qualityAttribute = "long_name";
  char qualityText[ 256 ] = "";
  size_t length = 0;
  int status = 0;
  int result = 0;

  assert( file >= 0 ); assert( id >= 0 );

  /*
   * Parse QCAll::long_name attribute text to determine the value used
   * to represent 'high' quality:
   */

  status = nc_inq_att( file, id, qualityAttribute, 0, &length );

  if ( status == NC_NOERR && length > 0 &&
       length < sizeof qualityText / sizeof *qualityText ) {
    status = nc_get_att_text( file, id, qualityAttribute, qualityText );

    if ( status == NC_NOERR ) {
      const char* c = strstr( qualityText, "high" );

      if ( c ) {
        --c;

        while ( c > qualityText && ( *c == ' ' || *c == ':' ) ) {
          --c;
        }

        while ( c > qualityText && isdigit( *c ) ) {
          --c;
        }

        ++c;
        result = atoi( c );
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "%s\n", message );
  }

  return result;
}




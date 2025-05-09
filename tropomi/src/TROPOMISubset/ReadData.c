/******************************************************************************
PURPOSE: ReadData.c - Simple to use wrapper routines to read data from TROPOMI
                      NetCDF4 files.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2018-04-04 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For strcmp(), strstr(). */
#include <float.h>  /* For DBL_MAX. */

#include <netcdf.h> /* For NC*, nc_*(). */

#include "Utilities.h" /* For LONGITUDE, Bounds, expand32BitReals(). */
#include "ReadData.h"  /* For public interface. */

/*================================== MACROS =================================*/

#define MIN( a, b ) ( ( a ) < ( b ) ? ( a ) : ( b ) )
#define MAX( a, b ) ( ( a ) > ( b ) ? ( a ) : ( b ) )
#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

#ifdef DEBUGGING
#define DEBUG( s ) s
#else
#define DEBUG(unused)
#endif

/*=========================== FORWARD DECLARATIONS ==========================*/

static size_t filterByQC( const int group,
                          const char* const qcVariable,
                          const size_t starts[],
                          const size_t counts[],
                          const size_t points,
                          const int qcMinimum,
                          float data[],
                          unsigned char qcData[] );

static size_t filterByGroundPixelRange( const int group,
                                        const char* const groundPixelVariable,
                                        const size_t starts[],
                                        const size_t counts[],
                                        const size_t points,
                                        const int groundPixelMinimum,
                                        const int groundPixelMaximum,
                                        float data[],
                                        int groundPixelData[] );

static size_t filterByCloudFraction( const int group,
                                     const char* const cloudVariable,
                                     const size_t starts[],
                                     const size_t counts[],
                                     const size_t points,
                                     const double maximumCloudFraction,
                                     float data[],
                                     float cloudData[] );

static size_t expandAndFilterData( const size_t points,
                                   const char* const units,
                                   const double scale,
                                   const double validMinimum,
                                   const double validMaximum,
                                   double data[] );

static int getNestedGroupId( const int group,
                             const char* const subgroups[],
                             const size_t count );

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

  DEBUG( fprintf( stderr, "%s\n", fileName ); )
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
PURPOSE: readFileBounds - Read file lon-lat bounds (extent of swath).
INPUTS:  const int file   NetCDF file ID to read.
OUTPUTS: Bounds bounds    Lon-lat bounds of swath.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileBounds( const int file, Bounds bounds ) {
  int result = 0;
  int status = 0;
  float longitudeMinimum = 0.0;
  float longitudeMaximum = 0.0;
  float latitudeMinimum  = 0.0;
  float latitudeMaximum  = 0.0;

  bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
  bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
  bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
  bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;

  status =
    nc_get_att_float( file, NC_GLOBAL, "geospatial_lon_min",
                      &longitudeMinimum );

  if ( status == NC_NOERR ) {
    status =
      nc_get_att_float( file, NC_GLOBAL, "geospatial_lon_max",
                        &longitudeMaximum );

    if ( status == NC_NOERR ) {
      status =
        nc_get_att_float( file, NC_GLOBAL, "geospatial_lat_min",
                          &latitudeMinimum );

      if ( status == NC_NOERR ) {
        status =
          nc_get_att_float( file, NC_GLOBAL, "geospatial_lat_max",
                            &latitudeMaximum );

        /* Files contain out-of-order longitudes so order them here: */

        bounds[ LONGITUDE ][ MINIMUM ] =
          MIN( longitudeMinimum, longitudeMaximum );
        bounds[ LONGITUDE ][ MAXIMUM ] =
          MAX( longitudeMinimum, longitudeMaximum );
        bounds[ LATITUDE  ][ MINIMUM ] =
          MIN( latitudeMinimum, latitudeMaximum );
        bounds[ LATITUDE  ][ MAXIMUM ] =
          MAX( latitudeMinimum, latitudeMaximum );
      }
    }
  }


  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read valid bounds because: %s\n", message);
  } else if ( ! isValidBounds( (const double (*)[2]) bounds ) ) {
    fprintf( stderr, "Failed to read valid bounds\n" );
  } else {
    result = 1;
  }

  if ( ! result ) {
    bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
    bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
    bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
    bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
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
  int gid = 0;
  int id = 0;
  int status = nc_inq_ncid( file, "PRODUCT", &gid );

  if ( status == NC_NOERR ) {
    status = nc_inq_dimid( gid, "scanline", &id );
    assert( rows ); assert( columns );
    *rows = *columns = 0;

    if ( status == NC_NOERR ) {
      status = nc_inq_dimlen( gid, id, rows );

      if ( status == NC_NOERR ) {
        status = nc_inq_dimid( gid, "ground_pixel", &id );

        if ( status == NC_NOERR ) {
          status = nc_inq_dimlen( gid, id, columns );
        }
      }
    }
  }

  result = status == NC_NOERR && *rows != 0 && *columns != 0;

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read valid dimensions because: %s\n", message);
  }

  if ( ! result ) {
    *columns = *rows = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file swath data.
INPUTS:  const int file                 NetCDF file ID to read.
         const char* const variable     Name of variable to read.
         const size_t rows              Rows of variable to read.
         const size_t columns           Columns of variable to read.
         const int qcMinimum            Minimum acceptable data quality [0,100]
         const int groundPixelMinimum   Minimum acceptable ground pixel or -1.
         const int groundPixelMaximum   Maximum acceptable ground pixel or -1.
         const double maximumCloudFraction  Maximum acceptable cloud fraction.
         const int allowNegativeCounts  Allow negative molecules/cm2? UGLY.
OUTPUTS: char units[ 80 ]               Units of variable.
         double data[ rows * columns ]  Swath data for variable.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
NOTES:   Data is filtered by qa_value - i.e., set to MISSING_VALUE
         if qa_value != 1 (full quality data).
******************************************************************************/

int readFileData( const int file,
                  const char* const variable,
                  const size_t rows,
                  const size_t columns,
                  const int qcMinimum,
                  const int groundPixelMinimum,
                  const int groundPixelMaximum,
                  const double maximumCloudFraction,
                  const int allowNegativeCounts,
                  char units[ 80 ],
                  double data[] ) {

  typedef struct {
    const char* const group;
    const char* const name;
    const char* const units;
    const double scale;
    const double validMinimum;
    const double validMaximum;
    const char* const group2;
    const char* const qcVariable;
    const char* const filterVariable;
  } Entry;
  const Entry table[] = {

    /* L2 NO2: */

    { "PRODUCT", "longitude", "deg", 1.0, -180.0, 180.0, 0, 0, 0 },

    { "PRODUCT", "latitude",  "deg", 1.0,  -90.0,  90.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/GEOLOCATIONS", "solar_zenith_angle",
      "deg", 1.0,    0.0, 180.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/GEOLOCATIONS", "solar_azimuth_angle",
      "deg", 1.0, -180.0, 180.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/GEOLOCATIONS", "viewing_zenith_angle",
      "deg", 1.0,    0.0, 180.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/GEOLOCATIONS", "viewing_azimuth_angle",
      "deg", 1.0, -180.0, 180.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "surface_altitude",
      "m", 1.0, -500.0, 10000.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "surface_albedo",
      "-", 1.0, 0.0, 1.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "surface_pressure",
      "hPa", 1e-2, 0.0, 1500.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "cloud_fraction_crb",
      "-", 1.0, 0.0, 1.0, 0, 0, 0 },

    { "PRODUCT", "nitrogendioxide_tropospheric_column", "molecules/cm2",
       6.022141e+19, 0.0, 1e30,
      "PRODUCT", "qa_value", "cloud_fraction_crb_nitrogendioxide_window" },

    { "PRODUCT", "nitrogendioxide_tropospheric_column_precision",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30,
      "PRODUCT", "qa_value", "cloud_fraction_crb_nitrogendioxide_window" },

    { "PRODUCT", "nitrogendioxide_tropospheric_column_precision_kernel",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30,
      "PRODUCT", "qa_value", "cloud_fraction_crb_nitrogendioxide_window" },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS", "nitrogendioxide_total_column",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30,
      "PRODUCT", "qa_value", "cloud_fraction_crb_nitrogendioxide_window" },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
      "nitrogendioxide_total_column_precision",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30,
      "PRODUCT", "qa_value", "cloud_fraction_crb_nitrogendioxide_window" },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
      "nitrogendioxide_total_column_precision_kernel",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30,
      "PRODUCT", "qa_value", "cloud_fraction_crb_nitrogendioxide_window" },

    { "PRODUCT", "air_mass_factor_troposphere", "-", 1.0, 0.0, 1e30,
      "PRODUCT", "qa_value", 0 },

    { "PRODUCT", "air_mass_factor_total",       "-", 1.0, 0.0, 1e30,
      "PRODUCT", "qa_value", 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
      "nitrogendioxide_stratospheric_column",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30,
      "PRODUCT", "qa_value", "cloud_fraction_crb_nitrogendioxide_window" },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
      "cloud_fraction_crb_nitrogendioxide_window",
      "-", 1.0, 0.0, 1.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS",
      "cloud_radiance_fraction_nitrogendioxide_window",
      "-", 1.0, 0.0, 1.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA",
      "surface_albedo_nitrogendioxide_window",
      "-", 1.0, 0.0, 1.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "surface_pressure",
      "hPa", 1e-2, 0.0, 1500.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "surface_classification",
      "-", 1.0, 0.0, 249.0, 0, 0, 0 },

    { "PRODUCT", "qa_value", "-", 0.01, 0.0, 1.0, 0, 0, 0 },

    /* RPRO L2 NO2 files have these additional variables: */

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "eastward_wind",
      "m/s", 1.0, -500.0, 500.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/INPUT_DATA", "northward_wind",
      "m/s", 1.0, -500.0, 500.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD",
      "o22cld_cloud_height_crb",
      "m", 1.0, -500.0, 1e5, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO",
      "fresco_cloud_pressure_crb",
      "hPa", 1e-2, 0.0, 1500.0, 0, 0, 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS", "air_mass_factor_cloudy",
      "-", 1.0, 0.0, 1e30,
      "PRODUCT", "qa_value", 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS", "air_mass_factor_clear",
      "-", 1.0, 0.0, 1e30,
      "PRODUCT", "qa_value", 0 },


    /* L2 HCHO: */

    { "PRODUCT", "formaldehyde_tropospheric_vertical_column", "molecules/cm2",
       6.022141e+19, 0.0, 1e30, "PRODUCT", "qa_value", 0 },

    { "PRODUCT", "formaldehyde_tropospheric_vertical_column_precision",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30, "PRODUCT", "qa_value", 0 },

    /* L2 CO: */

    { "PRODUCT", "carbonmonoxide_total_column", "molecules/cm2",
      6.022141e+19, 0.0, 1e30, "PRODUCT", "qa_value", 0 },

    { "PRODUCT", "carbonmonoxide_total_column_precision", "molecules/cm2",
      6.022141e+19, 0.0, 1e30, "PRODUCT", "qa_value", 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS", "water_total_column",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30, "PRODUCT", "qa_value", 0 },

    { "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS", "water_total_column_precision",
      "molecules/cm2", 6.022141e+19, 0.0, 1e30, "PRODUCT", "qa_value", 0 },

    /* L2 CH4: */

    { "PRODUCT", "methane_mixing_ratio",                "-", 1e-9, 0.0, 1.0,
      "PRODUCT", "qa_value", 0 },

    { "PRODUCT", "methane_mixing_ratio_precision",      "-", 1e-9, 0.0, 1.0,
      "PRODUCT", "qa_value", 0 },

    { "PRODUCT", "methane_mixing_ratio_bias_corrected", "-", 1e-9, 0.0, 1.0,
      "PRODUCT", "qa_value", 0 },

 };
  const size_t variables = sizeof table / sizeof *table;
  int status = 0;
  int result = 0;
  size_t validValues = 0; /* Number of non-filtered values read. */
  size_t index = 0;

  assert( file >= 0 ); assert( variable ); assert( *variable );
  assert( rows != 0 ); assert( columns != 0 );
  assert( qcMinimum >= 0 ); assert( qcMinimum <= 100 );
  assert( maximumCloudFraction >= 0.0 ); assert( maximumCloudFraction <= 1.0 );
  assert( allowNegativeCounts == 0 || allowNegativeCounts == 1 );
  assert( units ); assert( data );

  while ( index < variables && strcmp( variable, table[ index ].name ) ) {
    ++index;
  }

  DEBUG( fprintf( stderr, "%s index = %lu\n", variable, index ); )

  if ( index < variables ) {
    const Entry* const entry = table + index;
    int id = -1;
    int gid = -1; /* Groups are more trouble than they're worth! */
    int gid2 = -1;
    int gid0 = file;
    char* rest = 0;
    char* name = 0;
    char groupName[80] = "";
    memset( groupName, 0, sizeof groupName );
    strncpy( groupName, entry->group, sizeof groupName / sizeof *groupName );

    for ( name = strtok_r( groupName, "/", &rest );
          name;
          name = strtok_r( 0, "/", &rest ) ) {
      status = nc_inq_ncid( gid0, name, &gid );
      gid0 = gid;

      if ( gid2 == -1 ) {
        gid2 = gid;
      }
    }

    strcpy( units, entry->units );

    if ( status == NC_NOERR ) {
      status = nc_inq_varid( gid, variable, &id );

      if ( status == NC_NOERR ) {
        const size_t points = rows * columns;
        float* fdata = (float*) data;
        size_t starts[ 3 ] = { 0, 0, 0 };
        size_t counts[ 3 ] = { 1, 0, 0 };
        counts[ 1 ] = rows;
        counts[ 2 ] = columns;
        status = nc_get_vara_float( gid, id, starts, counts, fdata );

        if ( status == NC_NOERR ) {
          const char* const qcVariable = qcMinimum == 0 ? 0 : entry->qcVariable;
          validValues = points;
          DEBUG( fprintf( stderr, "points = %lu\n", points ); )

          if ( qcVariable ) {
            unsigned char* const qcData = (unsigned char*) ( fdata + points );
            validValues =
              filterByQC( gid2, qcVariable, starts, counts, points,
                          qcMinimum, fdata, qcData );
            DEBUG( fprintf( stderr, "After filterByQC(%s), validValues = %lu\n",
                            qcVariable, validValues ); )
          }

          if ( validValues ) {
            const char* const groundPixelVariable =
              groundPixelMinimum == -1 ? 0 : "ground_pixel";

            if ( groundPixelVariable ) {
              int groundPixelGroupId = -1;
              status = nc_inq_ncid( file, "PRODUCT", &groundPixelGroupId );

              if ( status == NC_NOERR ) {
                int* const groundPixelData = (int*) ( fdata + points );
                validValues =
                  filterByGroundPixelRange( groundPixelGroupId,
                                            groundPixelVariable,
                                            starts, counts, points,
                                            groundPixelMinimum,
                                            groundPixelMaximum,
                                            fdata, groundPixelData );
                DEBUG( fprintf( stderr,
                                "After filterByGroundPixelRange(%s), "
                                "validValues = %lu\n",
                                groundPixelVariable, validValues ); )
              }
            }

            if ( status == NC_NOERR && validValues ) {
              const char* const filterVariable =
                maximumCloudFraction == 1.0 ? 0 : entry->filterVariable;

              if ( filterVariable ) {
                float* const filterData = fdata + points;
                validValues =
                  filterByCloudFraction( gid2, filterVariable, starts, counts,
                                         points, maximumCloudFraction, fdata,
                                         filterData );
                DEBUG( fprintf( stderr, "After filterByCloudFraction(%s), "
                                "validValues = %lu\n",
                                filterVariable, validValues ); )
              }

              if ( validValues ) {
                const double scale = entry->scale;
                const double validMinimum =
                  allowNegativeCounts &&
                  ! strcmp( entry->units, "molecules/cm2" ) ? -1e29
                  : entry->validMinimum;
                const double validMaximum = entry->validMaximum;
                validValues =
                  expandAndFilterData( points, units, scale,
                                       validMinimum, validMaximum, data );
                DEBUG( fprintf( stderr,
                                "After expandAndFilterData(%s, %f, [%e, %e]), "
                                "validValues = %lu\n",
                                units, scale, validMinimum, validMaximum,
                                validValues ); )
              }
            }
          }
        }
      }
    }
  }

  if ( status ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read variable %s because: %s\n",
             variable, message );
  }

  result = validValues > 0;
  DEBUG( fprintf( stderr, "result = %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: filterByQC - Read qc array and use it to filter data.
INPUTS:  const int group               NetCDF file group ID of qcVariable.
         const char* const qcVariable  Name of qc variable to read.
         const size_t starts[]         Start dimension indices of variable.
         const size_t counts[]         Dimensions of variable.
         const size_t points           Total number of variable data points.
         const int qcMinimum           Minimum acceptable data quality [0,100]
         float data[ points ]          Swath data to filter.
OUTPUTS: float data[ points ]          Fitered swath data.
         unsigned char qcData[ points ]  Storage to hold qcData.
RETURNS: size_t number of unfiltered points,
         else 0 and a failure message is printed to stderr.
******************************************************************************/

static size_t filterByQC( const int group,
                          const char* const qcVariable,
                          const size_t starts[],
                          const size_t counts[],
                          const size_t points,
                          const int qcMinimum,
                          float data[],
                          unsigned char qcData[] ) {

  const int qcMaximum = 100;
  size_t result = 0;
  int qcId = -1;
  int status = 0;

  assert( group >= 0 ); assert( qcVariable ); assert( *qcVariable );
  assert( starts ); assert( counts );
  assert( points );
  assert( IN_RANGE( qcMinimum, 0, 100 ) );
  assert( data ); assert( qcData );

  status = nc_inq_varid( group, qcVariable, &qcId );

  if ( status == NC_NOERR ) {
    status = nc_get_vara_uchar( group, qcId, starts, counts, qcData );

    if ( status == NC_NOERR ) {
      size_t point = 0;

      for ( point = 0; point < points; ++point ) {
        const unsigned char qc = qcData[ point ];

        if ( IN_RANGE( qc, qcMinimum, qcMaximum ) ) {
          ++result;
        } else {
          data[ point ] = MISSING_VALUE;
        }
      }
    }
  }

  if ( status ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read qc variable %s because: %s\n",
             qcVariable, message );
  }

  return result;
}



/******************************************************************************
PURPOSE: filterByGroundPixelRange - Read ground pixel array and use it to
         filter data.
INPUTS:  const int group           NetCDF file group ID of groundPixelVariable.
         const char* const groundPixelVariable  Name of ground pixel variable.
         const size_t starts[]        Start dimension indices of variable.
         const size_t counts[]        Dimensions of variable.
         const size_t points          Total number of variable data points.
         const int groundPixelMinimum   Minimum acceptable ground pixel value.
         const int groundPixelMaximum   Maximum acceptable ground pixel value.
OUTPUTS: float data[ points ]           Fitered swath data.
         int groundPixelData[ points ]  Storage to hold ground pixel data.
RETURNS: size_t number of unfiltered points,
         else 0 and a failure message is printed to stderr.
******************************************************************************/

static size_t filterByGroundPixelRange( const int group,
                                        const char* const groundPixelVariable,
                                        const size_t starts[],
                                        const size_t counts[],
                                        const size_t points,
                                        const int groundPixelMinimum,
                                        const int groundPixelMaximum,
                                        float data[],
                                        int groundPixelData[] ) {

  size_t result = 0;
  int id = -1;
  int status = 0;

  assert( group >= 0 );
  assert( groundPixelVariable ); assert( *groundPixelVariable );
  assert( starts ); assert( counts );
  assert( points );
  assert( groundPixelMinimum >= 0 );
  assert( groundPixelMaximum >= groundPixelMinimum );
  assert( data ); assert( groundPixelData );

  status = nc_inq_varid( group, groundPixelVariable, &id );

  if ( status == NC_NOERR ) {
    status =
      nc_get_vara_int( group, id, starts + 2, counts + 2, groundPixelData );

    if ( status == NC_NOERR ) {
      const size_t groundPoints = counts[ 2 ];
      size_t point = 0;

      for ( point = 0; point < points; ++point ) {
        const size_t groundPixelIndex = point % groundPoints;
        const int groundPixelValue = groundPixelData[ groundPixelIndex ];

        if ( IN_RANGE( groundPixelValue,
                       groundPixelMinimum, groundPixelMaximum ) ) {
          ++result;
        } else {
          data[ point ] = MISSING_VALUE;
        }
      }
    }
  }

  if ( status ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read ground pixel variable %s because: %s\n",
             groundPixelVariable, message );
  }

  return result;
}



/******************************************************************************
PURPOSE: filterByCloudFraction - Read qc array and use it to filter data.
INPUTS:  const int group                  NetCDF file group ID of variable.
         const char* const cloudVariable  Name of cloud variable to read.
         const size_t starts[]            Start dimension indices of variable.
         const size_t counts[]            Dimensions of variable.
         const size_t points              Total number of variable data points.
         const double maximumCloudFraction  Maximum acceptable cloud fraction.
         float data[ points ]             Swath data to filter.
OUTPUTS: float cloudData[ points ]        Storage to hold cloudData.
RETURNS: size_t number of unfiltered points,
         else 0 and a failure message is printed to stderr.
******************************************************************************/

static size_t filterByCloudFraction( const int group,
                                     const char* const cloudVariable,
                                     const size_t starts[],
                                     const size_t counts[],
                                     const size_t points,
                                     const double maximumCloudFraction,
                                     float data[],
                                     float cloudData[] ) {

  const char* const subgroups[] = {
    "SUPPORT_DATA",
    "DETAILED_RESULTS"
  };
  const double minimumCloudFraction = 0.0;
  size_t result = 0;
  int groupId = 0;
  int cloudId = -1;
  int status = 0;

  assert( group >= 0 ); assert( cloudVariable ); assert( *cloudVariable );
  assert( starts ); assert( counts );
  assert( points );
  assert( IN_RANGE( maximumCloudFraction, 0.0, 1.0 ) );
  assert( data ); assert( cloudData );

  groupId =
    getNestedGroupId(group, subgroups, sizeof subgroups / sizeof subgroups[0]);

  if ( groupId ) {
    status = nc_inq_varid( groupId, cloudVariable, &cloudId );

    if ( status == NC_NOERR ) {
      status = nc_get_vara_float( groupId, cloudId, starts, counts, cloudData);

      if ( status == NC_NOERR ) {
        size_t point = 0;

        for ( point = 0; point < points; ++point ) {

          if ( data[ point ] != MISSING_VALUE ) {
            const double cloudFraction = cloudData[ point ];

            if ( IN_RANGE( cloudFraction,
                           minimumCloudFraction, maximumCloudFraction ) ) {
              ++result;
            } else {
              data[ point ] = MISSING_VALUE;
            }
          }
        }
      }
    }
  }

  if ( status ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read cloud variable %s because: %s\n",
             cloudVariable, message );
  }

  return result;
}



/******************************************************************************
PURPOSE: expandAndFilterData - Expand float data to double and filter by
         valid range and convert values per units.
INPUTS:  const size_t points          Total number of variable data points.
         const char* const units      Units of values.
         const double scale           Scale data value.
         const double validMinimum    Minimum acceptable data value.
         const double validMaximum    Maximum acceptable data value.
         double data[ points ]        Float swath data to expand/filter.
OUTPUTS: double data[ points ]        Expanded/filtered data.
RETURNS: size_t number of unfiltered points,
         else 0 and a failure message is printed to stderr.
******************************************************************************/

static size_t expandAndFilterData( const size_t points,
                                   const char* const units,
                                   const double scale,
                                   const double validMinimum,
                                   const double validMaximum,
                                   double data[] ) {

  size_t result = 0;

  assert( points ); assert( units ); assert( *units );
  assert( scale > 0.0 );
  assert( validMinimum <= validMaximum );
  assert( data );

  {
    float* const fdata = (float*) data;
    size_t index = points;

    /* Loop backward to avoid overwrite when expanding 4 bytes to 8 bytes: */

    while ( index-- ) {
      const double value = fdata[ index ];

      if ( value != MISSING_VALUE ) {
        const double convertedValue = value * scale;

        if ( IN_RANGE( convertedValue, validMinimum, validMaximum ) ) {
          data[ index ] = convertedValue;
          ++result;
        } else {
          data[ index ] = MISSING_VALUE;
        }
      } else {
        data[ index ] = MISSING_VALUE; /* Store 64-bit missing value. */
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: getNestedGroupId - Get NetCDF group id of nested group.
INPUTS:  const int group              NetCDF group ID to start from.
         const char* const subgroups[] Names of nested subgroups.
         const size_t count            Number of subgroups.
RETURNS: int NetCDF id of nested subgroup,
         else 0 and a failure message is printed to stderr.
NOTES:   Nested groups are _much_ more trouble than they're worth!
******************************************************************************/

static int getNestedGroupId( const int group,
                             const char* const subgroups[],
                             const size_t count ) {

  int result = group;
  int status = 0;
  size_t index = 0;

  assert( group >= 0 ); assert( count );
  assert( subgroups ); assert( subgroups[ 0 ] );
  assert( subgroups[ count - 1 ] );

  for ( index = 0; status == 0 && index < count; ++index ) {
    int nestedGroupId = -1;
    status = nc_inq_ncid( result, subgroups[ index ] , &nestedGroupId );
    result = nestedGroupId;
  }

  if ( status ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read nested groups because: %s\n",
             message );
  }

  return result;
}



/******************************************************************************
PURPOSE: ReadData.c - Simple to use wrapper routines to read data from TEMPO
                      NetCDF4 files.

NOTES:   Uses NetCDF4 and HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2019-06-25 plessel.todd@epa.gov
STATUS:  unreviewed untested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For strcmp(), strstr(). */
#include <float.h>  /* For DBL_MAX. */

#include <netcdf.h> /* For NC*, nc_*(). */

#include "Utilities.h" /* For MISSING_VALUE. */
#include "ReadData.h"  /* For public interface. */

/*========================== FORWARD DECLARATIONS ===========================*/

static int findVariable( const int parentId, const char* const name,
                         int* gid, int* rank, nc_type* type );

static size_t readAndExpandData( const int gid,
                                 const int id,
                                 const nc_type type,
                                 const size_t starts[],
                                 const size_t counts[],
                                 const double validMinimum,
                                 const double validMaximum,
                                 double data[] );

static size_t sumDataAndScratch( const int gid,
                                 const int id,
                                 const size_t starts[],
                                 const size_t counts[],
                                 const double validMinimum,
                                 const double validMaximum,
                                 double scratch[],
                                 double data[] );

static size_t filterDataByQC( const int file,
                              const char* const qcVariable,
                              const size_t starts[],
                              const size_t counts[],
                              const int qcMinimum,
                              const int qcMaximum,
                              const int shift,
                              double scratch[],
                              double data[] );

static size_t filterDataByAuxiliaryVariable( const int file,
                                            const char* const auxiliaryVariable,
                                             const size_t starts[],
                                             const size_t counts[],
                                             const double maximum,
                                             double scratch[],
                                             double data[] );

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
PURPOSE: readFileDimensions - Read file swath row/column dimensions.
INPUTS:  const int file            NetCDF file ID to read.
         const char* const product   "L2" or "L3".
         const char* const variable  Variable name. "pm25_ge"
OUTPUTS: size_t* rows              Rows in swath.
         size_t* columns           Columns in swath.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileDimensions( const int file,
                        const char* const product,
                        const char* const variable,
                        size_t* const rows,
                        size_t* const columns ) {
  int result = 0;
  int id = 0;
  int status = NC_NOERR - 1;
  assert( file > -1 );
  assert( product ); assert( *product );
  assert( strstr( product, "L2" ) || strstr( product, "L3" ) );
  assert( rows ); assert( columns );
  *rows = *columns = 0;

  if ( strstr( product, "L2" ) ) {
    status = nc_inq_dimid( file, "mirror_step", &id );
  } else if ( strstr( product, "PM25_L3" ) ) {

    if ( strstr( variable, "_ge" ) ) {
      status = nc_inq_dimid( file, "xdim_ge", &id );
    } else {
      status = nc_inq_dimid( file, "xdim_gw", &id );
    }
  } else {
    status = nc_inq_dimid( file, "longitude", &id );
  }

  if ( status == NC_NOERR ) {
    status = nc_inq_dimlen( file, id, columns );

    if ( status == NC_NOERR ) {

      if ( strstr( product, "L2" )  ) {
        status = nc_inq_dimid( file, "xtrack", &id );
      } else if ( strstr( product, "PM25_L3" ) ) {

        if ( strstr( variable, "_ge" ) ) {
          status = nc_inq_dimid( file, "ydim_ge", &id );
        } else {
          status = nc_inq_dimid( file, "ydim_gw", &id );
        }
      } else {
        status = nc_inq_dimid( file, "latitude", &id );
      }

      if ( status == NC_NOERR ) {
        status = nc_inq_dimlen( file, id, rows );
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
         const char* const product      File product e.g., "NO2_L2" or "L3".
         const char* const variable     E.g., "vertical_column_total".
         const size_t rows              Rows of variable to read.
         const size_t columns           Columns of variable to read.
         const size_t gridSubsetIndices[2][2]  L3 grid subset 1-based indices.
         const int minimumQuality       Minimum acceptable data quality
                                        0 = normal, 1 = suspect, 2 = bad.
         const double maximumCloudFraction  Maximum allowed cloud fraction.
         const double maximumCSolarZenithAngle  Max allowed solar zenith angle.
         const int allowNegativeCounts  Allow negative molecules/cm2? UGLY.
         double scratch[ rows * columns ]   Temp buffer for reordering.
OUTPUTS: char units[ 80 ]               Units of variable.
         double data[ rows * columns ]  Swath data for variable.
RETURNS: size_t number of unfiltered points,
         else 0 and a failure message is printed to stderr.
NOTES:   Data is filtered by qa_value - i.e., set to MISSING_VALUE
         if qa_value != 1 (full quality data).
******************************************************************************/

size_t readFileData( const int file,
                     const char* const product,
                     const char* const variable,
                     const size_t rows,
                     const size_t columns,
                     const size_t gridSubsetIndices[2][2],
                     const int minimumQuality,
                     const double maximumCloudFraction,
                     const double maximumSolarZenithAngle,
                     const int allowNegativeCounts,
                     char units[ 80 ],
                     double data[],
                     double scratch[] ) {

  typedef struct {
    const char* const product;    /* E.g., "NO2_L2". */
    const char* const group;      /* E.g., "product". */
    const char* const name;       /* E.g., "vertical_column_total". */
    const char* const units;      /* E.g., "molecules/cm2". */
    const double validMinimum;    /* E.g., 0.0. */
    const double validMaximum;    /* E.g., 1e30. */
    const char* const qcVariable; /* E.g., "main_data_quality_flag". */
    const char* const cloudFractionVariable; /* "eff_cloud_fraction" or "fc" */
    const char* const solarZenithAngleVariable; /* 0 or "solar_zenith_angle" */
  } Entry;
  const Entry table[] = {

    /* NO2_L2: ------------------------------------------------------------- */

    { "NO2_L2", "geolocation", "longitude",              "deg", -180.0, 180.0, 0, 0, 0},
    { "NO2_L2", "geolocation", "latitude",               "deg",  -90.0,  90.0, 0, 0, 0},
    { "NO2_L2", "geolocation", "solar_zenith_angle",     "deg",    0.0,  90.0, 0, 0, 0},
    { "NO2_L2", "geolocation", "solar_azimuth_angle",    "deg", -180.0, 180.0, 0, 0, 0},
    { "NO2_L2", "geolocation", "viewing_zenith_angle",   "deg",    0.0,  90.0, 0, 0, 0},
    { "NO2_L2", "geolocation", "viewing_azimuth_angle",  "deg", -180.0, 180.0, 0, 0, 0},
    { "NO2_L2", "geolocation", "relative_azimuth_angle", "deg", -180.0, 180.0, 0, 0, 0},

    /*
     * vertical_column_sum is a pseudo-variable to be computed as the sum of
     * vertical_column_troposphere and vertical_column_stratosphere:
     */

    { "NO2_L2", "product", "vertical_column_sum", "molecules/cm2", 0.0, 1e30,
      "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },

    { "NO2_L2", "product", "vertical_column_total", "molecules/cm2", 0.0, 1e30,
      "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L2", "product", "vertical_column_total_uncertainty", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L2", "product", "vertical_column_troposphere", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L2", "product", "vertical_column_troposphere_uncertainty",
      "molecules/cm2", 0.0, 1e30,
      "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L2", "product", "vertical_column_stratosphere", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },

    { "NO2_L2", "support_data", "fitted_slant_column", "molecules/cm2",
      0.0, 1e30, 0, 0, 0 },
    { "NO2_L2", "support_data", "fitted_slant_column_uncertainty",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L2", "support_data", "snow_ice_fraction",    "-",       0.0,  1.0, 0, 0, 0 },
    { "NO2_L2", "support_data", "terrain_height",        "m",  -1000.0, 10000.0, 0, 0, 0},
    { "NO2_L2", "support_data", "ground_pixel_quality_flag", "-",  0.0,  1e30,   0, 0, 0},
    { "NO2_L2", "support_data", "tropopause_pressure", "hPa",      0.0,  1030.0, 0, 0, 0},
    { "NO2_L2", "support_data", "surface_pressure",    "hPa",      0.0,  1030.0, 0, 0, 0},
    { "NO2_L2", "support_data", "albedo", "-", 0.0, 1.0, 0, 0, 0 },
    { "NO2_L2", "support_data", "amf_total",           "-",  0.0, 1e30, 0, 0, 0},
    { "NO2_L2", "support_data", "amf_diagnostic_flag", "-", -2.0, 500.0, 0, 0, 0 },
    { "NO2_L2", "support_data", "eff_cloud_fraction", "-", 0.0, 1.0, 0, 0, 0 },
    { "NO2_L2", "support_data", "amf_cloud_fraction", "-", 0.0, 1.0, 0, 0, 0 },
    { "NO2_L2", "support_data", "amf_cloud_pressure", "hPa", 0.0, 1030.0, 0, 0, 0 },
    { "NO2_L2", "support_data", "amf_troposphere", "-",       0.0, 1e30, 0, 0, 0},
    { "NO2_L2", "support_data", "amf_stratosphere", "-",      0.0, 1e30, 0, 0, 0},

    { "NO2_L2", "qa_statistics", "fit_rms_residual", "-", 0.0, 0.01, 0, 0, 0 },
    { "NO2_L2", "qa_statistics", "fit_convergence_flag", "-", -10.0, 12344.0, 0, 0, 0 },



    /* HCHO_L2: ------------------------------------------------------------ */

    { "HCHO_L2", "geolocation", "longitude",              "deg", -180.0,  180.0, 0,  0, 0},
    { "HCHO_L2", "geolocation", "latitude",               "deg",  -90.0,   90.0, 0,  0, 0},
    { "HCHO_L2", "geolocation", "solar_zenith_angle",     "deg",    0.0,   90.0, 0,  0, 0},
    { "HCHO_L2", "geolocation", "solar_azimuth_angle",    "deg", -180.0,  180.0, 0,  0, 0},
    { "HCHO_L2", "geolocation", "viewing_zenith_angle",   "deg",    0.0,   90.0, 0,  0, 0},
    { "HCHO_L2", "geolocation", "viewing_azimuth_angle",  "deg", -180.0,  180.0, 0,  0, 0},
    { "HCHO_L2", "geolocation", "relative_azimuth_angle", "deg", -180.0,  180.0, 0,  0, 0},

    { "HCHO_L2", "product", "vertical_column", "molecules/cm2", 0.0, 1e30,
                 "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "HCHO_L2", "product", "vertical_column_uncertainty", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },

    { "HCHO_L2", "support_data", "fitted_slant_column", "molecules/cm2",
      0.0, 1e30, 0, 0, 0 },
    { "HCHO_L2", "support_data", "fitted_slant_column_uncertainty",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "HCHO_L2", "support_data", "snow_ice_fraction",    "-",       0.0,  1.0, 0, 0, 0 },
    { "HCHO_L2", "support_data", "terrain_height",        "m",  -1000.0, 10000.0, 0, 0, 0},
    { "HCHO_L2", "support_data", "ground_pixel_quality_flag", "-", 0.0, 1e30, 0, 0, 0 },
    { "HCHO_L2", "support_data", "surface_pressure", "hPa", 0.0, 1030.0, 0, 0, 0 },
    { "HCHO_L2", "support_data", "albedo", "-", 0.0, 1.0, 0, 0, 0 },
    { "HCHO_L2", "support_data", "amf", "-", 0.0, 1e30, 0, 0, 0 },
    { "HCHO_L2", "support_data", "amf_diagnostic_flag", "-", -2.0, 500.0, 0, 0, 0 },
    { "HCHO_L2", "support_data", "eff_cloud_fraction", "-", 0.0, 1.0, 0, 0, 0 },
    { "HCHO_L2", "support_data", "amf_cloud_fraction", "-", 0.0, 1.0, 0, 0, 0 },
    { "HCHO_L2", "support_data", "amf_cloud_pressure", "hPa", 0.0, 1030.0, 0, 0, 0 },

    { "HCHO_L2", "qa_statistics", "fit_rms_residual", "-", 0.0, 0.01, 0, 0, 0 },
    { "HCHO_L2", "qa_statistics", "fit_convergence_flag", "-", -10.0, 12344.0, 0, 0, 0 },



    /* O3TOT_L2: ----------------------------------------------------------- */

    { "O3TOT_L2", "geolocation", "longitude",              "deg", -180.0, 180.0, 0, 0, 0},
    { "O3TOT_L2", "geolocation", "latitude",               "deg",  -90.0,  90.0, 0, 0, 0},
    { "O3TOT_L2", "geolocation", "solar_zenith_angle",     "deg",    0.0,  90.0, 0, 0, 0},
    { "O3TOT_L2", "geolocation", "solar_azimuth_angle",    "deg", -180.0, 180.0, 0, 0, 0},
    { "O3TOT_L2", "geolocation", "viewing_zenith_angle",   "deg",    0.0,  90.0, 0, 0, 0},
    { "O3TOT_L2", "geolocation", "viewing_azimuth_angle",  "deg", -180.0, 180.0, 0, 0, 0},
    { "O3TOT_L2", "geolocation", "relative_azimuth_angle", "deg", -180.0, 180.0, 0, 0, 0},
    { "O3TOT_L2", "geolocation", "terrain_height", "m", -500.0,  10000.0, 0, 0, 0 },

    { "O3TOT_L2", "product", "column_amount_o3", "DU", 0.0, 700.0,
      "quality_flag", "fc", "solar_zenith_angle" },
    { "O3TOT_L2", "product", "radiative_cloud_frac", "-", 0.0, 1.0, 0, 0, 0 },
    { "O3TOT_L2", "product", "fc", "-", 0.0, 1.0, 0, 0, 0 },
    { "O3TOT_L2", "product", "o3_below_cloud", "DU", 0.0, 700.0,
      "quality_flag", "fc", "solar_zenith_angle" },
    { "O3TOT_L2", "product", "quality_flag", "-", 0.0, 32768.0, 0, 0, 0 },
    { "O3TOT_L2", "product", "so2_index", "-", -300.0, 300.0,
      0, "fc", "solar_zenith_angle" },
    { "O3TOT_L2", "product", "uv_aerosol_index", "-", -30.0, 30.0,
      0, "fc", "solar_zenith_angle" },

    { "O3TOT_L2", "support_data", "ground_pixel_quality_flag", "-", 0.0, 32768.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "lut_wavelength", "nm", 300.0, 400.0, 0, 0, 0},
    { "O3TOT_L2", "support_data", "cloud_pressure", "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "terrain_pressure", "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "algorithm_flags", "-", 0.0, 13.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "radiance_bpix_flag_accepted", "-",
      0.0, 32768.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "surface_reflectivity_at_331nm",
      "-", -15.0, 115.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "surface_reflectivity_at_360nm",
      "-", -15.0, 115.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "step1_o3", "DU", 0.0, 700.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "step2_o3", "DU", 0.0, 700.0, 0, 0, 0 },
    { "O3TOT_L2", "support_data", "cal_adjustment", "-", -10.0, 10.0, 0, 0, 0 },



    /* CLDO4_L2: ----------------------------------------------------------- */

    { "CLDO4_L2", "geolocation", "longitude",              "deg", -180.0, 180.0, 0, 0, 0},
    { "CLDO4_L2", "geolocation", "latitude",               "deg",  -90.0,  90.0, 0, 0, 0},
    { "CLDO4_L2", "geolocation", "solar_zenith_angle",     "deg",    0.0,  90.0, 0, 0, 0},
    { "CLDO4_L2", "geolocation", "solar_azimuth_angle",    "deg", -180.0, 180.0, 0, 0, 0},
    { "CLDO4_L2", "geolocation", "viewing_zenith_angle",   "deg",    0.0,  90.0, 0, 0, 0},
    { "CLDO4_L2", "geolocation", "viewing_azimuth_angle",  "deg", -180.0, 180.0, 0, 0, 0},
    { "CLDO4_L2", "geolocation", "relative_azimuth_angle", "deg", -180.0, 180.0, 0, 0, 0},

    { "CLDO4_L2", "product", "CloudRadianceFraction440", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L2", "product", "CloudRadianceFraction466", "-", 0.0, 1.0, 0, 0, 0},
    { "CLDO4_L2", "product", "cloud_fraction", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L2", "product", "cloud_pressure", "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "CLDO4_L2", "product", "processing_quality_flag", "-", 0.0, 32767.0, 0, 0, 0 },

    { "CLDO4_L2", "qa_statistics", "fit_convergence_flag", "-", -10.0, 12344.0, 0, 0, 0},
    { "CLDO4_L2", "qa_statistics", "fit_rms_residual", "-", 0.0, 0.01, 0, 0, 0 },

    { "CLDO4_L2", "support_data", "GLER440", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "GLER466", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "SCD_MainDataQualityFlags", "-", 0.0, 2.0, 0, 0, 0},
    { "CLDO4_L2", "support_data", "SceneLER440", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "SceneLER466", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "ScenePressure", "hPa", 0.0, 1500.0, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "surface_pressure", "hPa", 0.0, 1500.0, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "fitted_slant_column", "molecules2/cm5",
      0.0, 1e30, "SCD_MainDataQualityFlags", 0, 0 },
    { "CLDO4_L2", "support_data", "fitted_slant_column_uncertainty",
      "molecules2/cm5", 0.0, 1e30, "SCD_MainDataQualityFlags", 0, 0 },
    { "CLDO4_L2", "support_data", "ground_pixel_quality_flag", "-", 0.0, 1e30, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "snow_ice_fraction", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L2", "support_data", "terrain_height", "m", -500.0, 10000.0, 0, 0, 0 },



    /* NO2_L3: ------------------------------------------------------------- */

    { "NO2_L3", 0, "longitude", "deg", -180.0, 180.0, 0, 0, 0 },
    { "NO2_L3", 0, "latitude",  "deg",  -90.0,  90.0, 0, 0, 0 },
    { "NO2_L3", 0, "weight",   "km2",   0.0,  1e30, 0, 0, 0 },

    { "NO2_L3", "geolocation", "solar_zenith_angle",     "deg",    0.0,  90.0, 0, 0, 0},
    { "NO2_L3", "geolocation", "viewing_zenith_angle",   "deg",    0.0,  90.0, 0, 0, 0},
    { "NO2_L3", "geolocation", "relative_azimuth_angle", "deg", -180.0, 180.0, 0, 0, 0},

    /*
     * vertical_column_sum is a pseudo-variable to be computed as the sum of
     * vertical_column_troposphere and vertical_column_stratosphere:
     */

    { "NO2_L3", "product", "vertical_column_sum", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },

    { "NO2_L3", "product", "vertical_column_troposphere", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L3", "product", "vertical_column_troposphere_uncertainty",
      "molecules/cm2", 0.0, 1e30,
      "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L3", "product", "vertical_column_stratosphere", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L3", "product", "vertical_column_total", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L3", "product", "vertical_column_total_uncertainty", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "NO2_L3", "product", "main_data_quality_flag", "-", 0.0, 2.0, 0, 0, 0 },

    { "NO2_L3", "qa_statistics", "num_vertical_column_troposphere_samples",
      "-", 0.0, 1e10, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "min_vertical_column_troposphere_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "max_vertical_column_troposphere_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "num_vertical_column_troposphere_uncertainty_samples",
      "-", 0.0, 1e10, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "min_vertical_column_troposphere_uncertainty_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "max_vertical_column_troposphere_uncertainty_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "num_vertical_column_stratosphere_samples",
      "-", 0.0, 1e10, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "min_vertical_column_stratosphere_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "max_vertical_column_stratosphere_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "num_vertical_column_total_samples",
      "-", 0.0, 1e10, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "min_vertical_column_total_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "qa_statistics", "max_vertical_column_total_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },

    { "NO2_L3", "support_data", "surface_pressure",    "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "NO2_L3", "support_data", "terrain_height",      "m", -500.0, 10000.0, 0, 0, 0},
    { "NO2_L3", "support_data", "snow_ice_fraction",   "-",  0.0, 1.0, 0, 0, 0 },
    { "NO2_L3", "support_data", "fitted_slant_column", "molecules/cm2",
      0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "support_data", "fitted_slant_column_uncertainty",
                "molecules/cm2",  0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "support_data", "albedo",              "-",  0.0, 1.0, 0, 0, 0 },
    { "NO2_L3", "support_data", "tropopause_pressure", "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "NO2_L3", "support_data", "amf_total",           "-", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "support_data", "eff_cloud_fraction",  "-", 0.0, 1.0, 0, 0, 0 },
    { "NO2_L3", "support_data", "amf_cloud_fraction",  "-", 0.0, 1.0, 0, 0, 0 },
    { "NO2_L3", "support_data", "amf_cloud_pressure",  "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "NO2_L3", "support_data", "amf_troposphere",     "-", 0.0, 1e30, 0, 0, 0 },
    { "NO2_L3", "support_data", "amf_stratosphere",    "-", 0.0, 1e30, 0, 0, 0 },



    /* HCHO_L3: ------------------------------------------------------------ */

    { "HCHO_L3", 0, "longitude", "deg", -180.0, 180.0, 0, 0, 0 },
    { "HCHO_L3", 0, "latitude",  "deg",  -90.0,  90.0, 0, 0, 0 },
    { "HCHO_L3", 0, "weight",   "km2",   0.0,  1e30, 0, 0, 0 },

    { "HCHO_L3", "geolocation", "solar_zenith_angle",     "deg",    0.0,  90.0, 0, 0, 0},
    { "HCHO_L3", "geolocation", "viewing_zenith_angle",   "deg",    0.0,  90.0, 0, 0, 0},
    { "HCHO_L3", "geolocation", "relative_azimuth_angle", "deg", -180.0, 180.0, 0, 0, 0},

    { "HCHO_L3", "product", "vertical_column", "molecules/cm2",
      0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "HCHO_L3", "product", "vertical_column_uncertainty",
      "molecules/cm2", 0.0, 1e30, "main_data_quality_flag", "eff_cloud_fraction", "solar_zenith_angle" },
    { "HCHO_L3", "product", "main_data_quality_flag", "-", 0.0, 2.0, 0, 0, 0 },

    { "HCHO_L3", "qa_statistics", "num_vertical_column_total_samples",
      "-", 0.0, 1e10, 0, 0, 0 },
    { "HCHO_L3", "qa_statistics", "min_vertical_column_total_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },
    { "HCHO_L3", "qa_statistics", "max_vertical_column_total_sample",
      "molecules/cm2", 0.0, 1e30, 0, 0, 0 },

    { "HCHO_L3", "support_data", "surface_pressure",    "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "HCHO_L3", "support_data", "terrain_height",      "m", -500.0, 10000.0, 0, 0, 0},
    { "HCHO_L3", "support_data", "snow_ice_fraction",   "-",  0.0, 1.0, 0, 0, 0 },
    { "HCHO_L3", "support_data", "fitted_slant_column", "molecules/cm2",
      0.0, 1e30, 0, 0, 0 },
    { "HCHO_L3", "support_data", "fitted_slant_column_uncertainty",
                "molecules/cm2",  0.0, 1e30, 0, 0, 0 },
    { "HCHO_L3", "support_data", "albedo",              "-",  0.0, 1.0, 0, 0, 0 },
    { "HCHO_L3", "support_data", "amf",                 "-", 0.0, 1e30, 0, 0, 0 },
    { "HCHO_L3", "support_data", "eff_cloud_fraction",  "-", 0.0, 1.0, 0, 0, 0 },
    { "HCHO_L3", "support_data", "amf_cloud_fraction",  "-", 0.0, 1.0, 0, 0, 0 },
    { "HCHO_L3", "support_data", "amf_cloud_pressure",  "hPa", 0.0, 1200.0, 0, 0, 0 },



    /* O3TOT_L3: ----------------------------------------------------------- */

    { "O3TOT_L3", 0, "longitude", "deg", -180.0, 180.0, 0, 0, 0 },
    { "O3TOT_L3", 0, "latitude",  "deg",  -90.0,  90.0, 0, 0, 0 },
    { "O3TOT_L3", 0, "weight",   "km2",   0.0,  1e30, 0, 0, 0 },

    { "O3TOT_L3", "geolocation", "solar_zenith_angle",     "deg",    0.0,  90.0, 0, 0, 0},
    { "O3TOT_L3", "geolocation", "viewing_zenith_angle",   "deg",    0.0,  90.0, 0, 0, 0},
    { "O3TOT_L3", "geolocation", "relative_azimuth_angle", "deg", -180.0, 180.0, 0, 0, 0},
    { "O3TOT_L3", "geolocation", "terrain_height", "m", -500.0, 10000.0, 0, 0, 0 },

    { "O3TOT_L3", "product", "column_amount_o3", "DU", 0.0, 700.0,
      0, "fc", "solar_zenith_angle" },
    { "O3TOT_L3", "product", "radiative_cloud_frac", "-", 0.0, 1.0, 0, 0, 0 },
    { "O3TOT_L3", "product", "fc", "-", 0.0, 1.0, 0, 0, 0 },
    { "O3TOT_L3", "product", "o3_below_cloud", "DU", 0.0, 100.0,
      0, "fc", "solar_zenith_angle" },
    { "O3TOT_L3", "product", "so2_index", "-", -300.0, 300.0,
      0, "fc", "solar_zenith_angle" },
    { "O3TOT_L3", "product", "uv_aerosol_index", "-", -30.0, 30.0,
      0, "fc", "solar_zenith_angle" },

    { "O3TOT_L3", "qa_statistics", "num_column_samples", "-", 0.0, 1e30, 0, 0, 0 },
    { "O3TOT_L3", "qa_statistics", "min_column_sample", "DU", 0.0, 700.0, 0, 0, 0 },
    { "O3TOT_L3", "qa_statistics", "max_column_sample", "DU", 0.0, 700.0, 0, 0, 0 },

    { "O3TOT_L3", "support_data", "cloud_pressure",    "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "O3TOT_L3", "support_data", "terrain_pressure",  "hPa", 0.0, 1200.0, 0, 0, 0 },


    /* CLDO4_L3: ----------------------------------------------------------- */

    { "CLDO4_L3", 0, "longitude", "deg", -180.0, 180.0, 0, 0, 0 },
    { "CLDO4_L3", 0, "latitude",  "deg",  -90.0,  90.0, 0, 0, 0 },
    { "CLDO4_L3", 0, "weight",   "km2",   0.0,  1e30, 0, 0, 0 },

    { "CLDO4_L3", "geolocation", "solar_zenith_angle",     "deg",    0.0,  90.0, 0, 0, 0},
    { "CLDO4_L3", "geolocation", "viewing_zenith_angle",   "deg",    0.0,  90.0, 0, 0, 0},
    { "CLDO4_L3", "geolocation", "relative_azimuth_angle", "deg", -180.0, 180.0, 0, 0, 0},

    { "CLDO4_L3", "product", "cloud_pressure", "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "CLDO4_L3", "product", "cloud_fraction", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L3", "product", "CloudRadianceFraction440", "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L3", "product", "CloudRadianceFraction466", "-", 0.0, 1.0, 0, 0, 0 },

    { "CLDO4_L3", "support_data", "surface_pressure",  "hPa", 0.0, 1200.0, 0, 0, 0 },
    { "CLDO4_L3", "support_data", "GLER440",           "-", 0.0, 1.0, 0, 0, 0 },
    { "CLDO4_L3", "support_data", "GLER446",           "-", 0.0, 1.0, 0, 0, 0 },


    /* AODALH_L2: -----------------------------------------------------------*/

    { "AODALH_L2", "geolocation", "longitude", "deg", -180.0, 180.0, 0, 0, 0},
    { "AODALH_L2", "geolocation", "latitude",  "deg",  -90.0,  90.0, 0, 0, 0},

    { "AODALH_L2", "product", "aod550", "-", -0.05, 10.0, "dqf", 0, 0 },
    { "AODALH_L2", "product", "alh", "km", -0.5, 15.0, "dqf", 0, 0 },
    { "AODALH_L2", "product", "aermodel", "-", 0.0, 3.0, 0, 0, 0 },

    { "AODALH_L2", "quality_diagnostic_flags", "lwmask", "-", 0.0, 3.0, 0,0,0 },
    { "AODALH_L2", "quality_diagnostic_flags", "qctest", "-", 0.0, 63.0, 0,0,0},
    { "AODALH_L2", "quality_diagnostic_flags", "dqf", "-", 0.0, 3.0, 0,0,0},


    /* PM25_L3: -----------------------------------------------------------*/

    { "PM25_L3", "geolocation", "lon_ge", "deg", -180.0, 180.0, 0, 0, 0},
    { "PM25_L3", "geolocation", "lat_ge",  "deg",  -90.0,  90.0, 0, 0, 0},
    { "PM25_L3", "geolocation", "lon_gw", "deg", -180.0, 180.0, 0, 0, 0},
    { "PM25_L3", "geolocation", "lat_gw",  "deg",  -90.0,  90.0, 0, 0, 0},

    { "PM25_L3", "product", "pm25sat_ge", "-", 0.0, 1000.0, 0, 0, 0 },
    { "PM25_L3", "product", "pm25sat_gw", "-", 0.0, 1000.0, 0, 0, 0 },

    { "PM25_L3", "support_data", "aod_ge", "-", -0.05, 5.0, 0,0,0 },
    { "PM25_L3", "support_data", "aod_gw", "-", -0.05, 5.0, 0,0,0 },
    { "PM25_L3", "support_data", "alh_ge", "km", -0.5, 15.0, 0,0,0 },
    { "PM25_L3", "support_data", "alh_gw", "km", -0.5, 15.0, 0,0,0 },
    { "PM25_L3", "support_data", "slope_aod_ge", "ug/m3", -1e3, 1e3, 0,0,0 },
    { "PM25_L3", "support_data", "slope_aod_gw", "ug/m3", -1e3, 1e3, 0,0,0 },
    { "PM25_L3", "support_data", "slope_alh_ge", "ug/m3/km", -1e3, 1e3, 0,0,0 },
    { "PM25_L3", "support_data", "slope_alh_gw", "ug/m2/km", -1e3, 1e3, 0,0,0 },
    { "PM25_L3", "support_data", "intercept_ge", "ug/m3", -1e3, 1e3, 0,0,0 },
    { "PM25_L3", "support_data", "intercept_gw", "ug/m3", -1e3, 1e3, 0,0,0 },
    { "PM25_L3", "support_data", "count_aod_ge", "-", 0.0, 60.0, 0,0,0 },
    { "PM25_L3", "support_data", "count_aod_gw", "-", 0.0, 60.0, 0,0,0 },


    /* ADP_L2: -----------------------------------------------------------*/

    { "ADP_L2", "geolocation", "longitude", "deg", -180.0, 180.0, 0, 0, 0},
    { "ADP_L2", "geolocation", "latitude", "deg", -180.0, 180.0, 0, 0, 0},

    { "ADP_L2", "product", "smoke", "-", 0.0, 1.0, "qc_flag", 0, 0 },
    { "ADP_L2", "product", "dust", "-", 0.0, 1.0, "qc_flag", 0, 0 },
    { "ADP_L2", "product", "cloud", "-", 0.0, 1.0, 0, 0, 0 },
    { "ADP_L2", "product", "nuc", "-", 0.0, 1.0, "qc_flag", 0, 0 },
    { "ADP_L2", "product", "snowice", "-", 0.0, 1.0, 0, 0, 0 },
    { "ADP_L2", "product", "saai", "-", 0.0, 30.0, 0, 0, 0 },
    { "ADP_L2", "product", "dsdi", "-", -50.0, 50.0, 0, 0, 0 },
    { "ADP_L2", "product", "deepblue_aai", "-", -30.0, 30.0, 0, 0, 0 },
    { "ADP_L2", "product", "uv_aai", "-", -50.0, 50.0, 0, 0, 0 },

    { "ADP_L2", "quality_diagnostic_flags", "qc_flag", "-", -128.0, 127.0,
      0, 0, 0 },
    { "ADP_L2", "quality_diagnostic_flags", "pqi1", "-", -128.0, 127.0,
      0, 0, 0 },
    { "ADP_L2", "quality_diagnostic_flags", "pqi2", "-", -128.0, 127.0,
      0, 0, 0 },
    { "ADP_L2", "quality_diagnostic_flags", "pqi3", "-", -128.0, 127.0,
      0, 0, 0 },
    { "ADP_L2", "quality_diagnostic_flags", "pqi4", "-", -128.0, 127.0,
      0, 0, 0 },
    { "ADP_L2", "quality_diagnostic_flags", "std_dev_410nm", "-", 0.0, 10.0,
      0, 0, 0 },
    { "ADP_L2", "quality_diagnostic_flags", "std_dev_865nm", "-", 0.0, 10.0,
      0, 0, 0 },
    { "ADP_L2", "quality_diagnostic_flags", "std_dev_2210nm", "-", 0.0, 10.0,
      0, 0, 0 },

    { 0, 0, 0, 0, 0.0, 0.0, 0, 0, 0 } /* End of table. */
  };
  const size_t variables = sizeof table / sizeof *table - 1;
  const int qcMinimum = 0;
  const int qcMaximum = minimumQuality; /* [0 = normal, 1 = suspect, 2 = bad]*/
  size_t result = 0;
  size_t index = 0;

  assert( file >= 0 );
  assert( product ), assert( *product );
  assert( strstr( product, "L2" ) || strstr( product, "L3" ) );
  assert( variable ); assert( *variable );
  assert( rows != 0 ); assert( columns != 0 );
  assert( qcMinimum >= 0 ); assert( qcMinimum <= qcMaximum );
  assert( maximumCloudFraction >= 0.0 ); assert( maximumCloudFraction <= 1.0 );
  assert( maximumSolarZenithAngle >= 0.0 );
  assert( maximumSolarZenithAngle <= 90.0 );
  assert( allowNegativeCounts == 0 || allowNegativeCounts == 1 );
  assert( units );
  assert( data );

  while ( index < variables &&
          ! ( strcmp( product, table[ index ].product ) == 0 &&
              strcmp( variable, table[ index ].name ) == 0 ) ) {
    ++index;
  }

  if ( index < variables ) {
    const Entry* const entry = table + index;
    int gid = -1;
    int rank = 0;
    nc_type type = -1;
    int id = -1;
    int gid2 = -1;
    int id2 = -1;
    strcpy( units, entry->units );

    if ( ! strcmp( variable, "vertical_column_sum" ) ) {
      id = findVariable( file, "vertical_column_troposphere",
                         &gid, &rank, &type );

      if ( id >= 0 && type == NC_DOUBLE ) {
        id2 = findVariable( file, "vertical_column_stratosphere",
                            &gid2, &rank, &type );

        if ( id2 < 0 || type != NC_DOUBLE ) {
          id = -1;
        }
      }

    } else {
      id = findVariable( file, variable, &gid, &rank, &type );
    }

    if ( id >= 0 ) {
      const double validMinimum =
        allowNegativeCounts &&
        ! strcmp( entry->units, "molecules/cm2" ) ? -1e29
        : entry->validMinimum;
      const double validMaximum = entry->validMaximum;
      const int isL3 = strstr( product, "L3" ) != 0 &&
                       strstr( product, "PM25_L3" ) == 0;
      const int isL3Longitude = isL3 && ! strcmp( variable, "longitude" );
      const int isL3Latitude  = isL3 && ! strcmp( variable, "latitude" );
      const int isGridSubset = isL3 && gridSubsetIndices[ 0 ][ 0 ] != 0;
      size_t starts[ 3 ] = { 0, 0, 0 };
      size_t counts[ 3 ] = { 1, 1, 1 };

      assert( ! isGridSubset ||
              columns == 1 + gridSubsetIndices[ COLUMN ][ LAST ] -
                             gridSubsetIndices[ COLUMN ][ FIRST ] );

      assert( ! isGridSubset ||
              rows == 1 + gridSubsetIndices[ ROW ][ LAST ] -
                          gridSubsetIndices[ ROW ][ FIRST ] );

      if ( isL3 ) {

        if ( isL3Longitude ) {
          counts[ 0 ] = columns;
          starts[ 0 ] =
            isGridSubset ? gridSubsetIndices[ COLUMN ][ FIRST ] - 1
            : 0;
        } else if ( isL3Latitude ) {
          counts[ 0 ] = rows;
          starts[ 0 ] =
            isGridSubset ? gridSubsetIndices[ ROW ][ FIRST ] - 1
            : 0;
        } else {
          counts[ 0 ] = rows;
          counts[ 1 ] = columns;
          starts[ 0 ] =
            isGridSubset ? gridSubsetIndices[ ROW ][ FIRST ] - 1
            : 0;
          starts[ 1 ] =
            isGridSubset ? gridSubsetIndices[ COLUMN ][ FIRST ] - 1
            : 0;
        }
      } else if ( ! strcmp( product, "PM25_L3" ) ) {
        counts[ 0 ] = rows;
        counts[ 1 ] = columns;
      } else {
        counts[ 0 ] = columns;
        counts[ 1 ] = rows;
      }

      /*
       * HACK: L3_*_V03 files have some variables with useless first dimension
       * time = 1. E.g., terrain_height(time, latitude, longitude)
       * so shift the subset:
       */

      if ( rank == 3 ) {
        counts[ 2 ] = counts[ 1 ];
        counts[ 1 ] = counts[ 0 ];
        counts[ 0 ] = 1;
        starts[ 2 ] = starts[ 1 ];
        starts[ 1 ] = starts[ 0 ];
        starts[ 0 ] = 0;
      }

      DEBUG( fprintf(stderr, "\nReadFileData(): %s, [%g, %g] starts = [%lu %lu %lu], counts = [%lu %lu %lu]\n",
                       entry->name, validMinimum, validMaximum,
                       starts[ 0 ], starts[ 1 ], starts[ 2 ],
                       counts[ 0 ], counts[ 1 ], counts[ 2 ] ); )

      /* Read data: */

      memset( data, 0, rows * columns * sizeof (double) );
      result = readAndExpandData( gid, id, type, starts, counts,
                                  validMinimum, validMaximum, data );

      if ( result ) {

        /* Replicate L3 1D coordinate arrays to 2D for consistency w/ L2: */

        if ( ! isGridSubset ) {

          if ( isL3Longitude ) {
            replicateRows( columns, rows, data );
          } else if ( isL3Latitude ) {
            replicateColumns( rows, columns, data );
          }
        }

        /* If vertical_column_sum then sum troposphere + stratosphere: */

        if ( ! strcmp( variable, "vertical_column_sum" ) ) {
          assert( type == NC_DOUBLE ); assert( scratch );
          result = sumDataAndScratch( gid2, id2, starts, counts,
                                      validMinimum, validMaximum,
                                      scratch, data );
        }

        if ( result ) {
          const char* const qcVariable =
            minimumQuality == 2 ? 0 : entry->qcVariable;
          const char* const cloudFractionVariable =
            maximumCloudFraction == 1.0 ? 0 : entry->cloudFractionVariable;
          const char* const solarZenithAngleVariable =
            maximumSolarZenithAngle == 90.0 ? 0
            : entry->solarZenithAngleVariable;

          if ( qcVariable ) { /* Read qc variable and use it to filter data:*/

            /*
             * TEMPO-ABI_ADP_Users_Guide_V1_20250107.pdf, Table 4 byte qc_flag
             * bits 2-3 = QC_SMOKE_CONFIDENCE
             * bits 4-5 = QC_DUST_CONFIDENCE
             * bits 6-7 = QC_NUC_CONFIDENCE
             * 00 = high, 01 = low, 10 = medium, 11 = bad/missing.
             * To extract the 2-bits, shift right then & 0x3.
             */

            const unsigned char shift =
              ! strcmp( variable, "smoke") ? 2
              : ! strcmp( variable, "dust") ? 4
              : ! strcmp( variable, "nuc") ? 6
              : 0;
            assert( scratch );
            result = filterDataByQC( file, qcVariable, starts, counts,
                                     qcMinimum, qcMaximum,
                                     shift, scratch, data );
          }

          /* Read cloud fraction variable and use it to filter data: */

          if ( result && cloudFractionVariable ) {
            assert( scratch );
            result =
              filterDataByAuxiliaryVariable( file, cloudFractionVariable,
                                             starts, counts,
                                             maximumCloudFraction,
                                             scratch, data );
          }

          /* Read solar zenith angle variable and use it to filter data: */

          if ( result && solarZenithAngleVariable ) {
            assert( scratch );
            result =
              filterDataByAuxiliaryVariable( file, solarZenithAngleVariable,
                                             starts, counts,
                                             maximumSolarZenithAngle,
                                             scratch, data );
          }

          /* Reorder L2 data from data[columns * rows] to data[rows * columns] */

          if ( result && ! isL3 && strcmp( product, "PM25_L3" ) ) {
            transpose( rows, columns, data, scratch );
          }
        }
      }
    }
  }

  DEBUG( fprintf( stderr, "readFileData() returning result = %lu\n\n", result );)

  return result;
}



/******************************************************************************
PURPOSE: findVariable - Find id and type of named variable.
INPUTS:  const int parentId      Id of file or group to search
                                 (non-tail recursively depth-first).
         const char* const name  Name of variable to find.
OUTPUTS: int* gid                Variable's group id.
         int* rank               Rank of variable.
         nc_type* type           Type of variable.
RETURNS: int id >= 0 if ok, else -1.
NOTES:   This is a non-tail recursive (stack-dangerous) depth-first search
         since NetCDF4/HDF files can have nested groups.
         Groups are not worth their complexity!
******************************************************************************/

static int findVariable( const int parentId, const char* const name,
                         int* gid, int* rank, nc_type* type ) {
  int result = -1;
  int status = nc_inq_varid( *gid = parentId, name, &result );

  if ( status == NC_NOERR && result > -1 ) {
    status = nc_inq_vartype( *gid, result, type );

    if ( status == NC_NOERR ) {
      status = nc_inq_varndims( *gid, result, rank );

#ifndef NDEBUG
      if ( status == NC_NOERR && *rank > 0 ) {
        char groupName[ NC_MAX_NAME ];
        memset( groupName, 0, sizeof groupName );
        status = nc_inq_grpname( *gid, groupName );

        if ( status == NC_NOERR ) {
          fprintf( stderr, "found group %s variable %s rank %d type %d\n",
                   groupName, name, *rank, (int) *type );
        }
      }
#endif
    }

    if ( status != NC_NOERR ) {
      result = -1;
      *gid = -1;
      *rank = 0;
    }

  } else {
    int count = 0;
    int ids[ NC_MAX_VARS ];
    result = -1;
    memset( ids, 0, sizeof ids );
    status = nc_inq_grps( parentId, &count, ids );

    if ( status == NC_NOERR ) {
      int index = 0;

      do {
        result = findVariable( ids[ index ], name, gid, rank, type );
        ++index;
      } while ( result == -1 && index < count );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readAndExpandData - Read named variable and expand to 64-bit doubles.
INPUTS:  const int gid               Group of variable.
         const int id                Id of variable.
         const nc_type type          Type of variable.
         const size_t starts[]       Subset start 0-based indices to read.
         const size_t counts[]       Subset counts to read.
         const double validMinimum   Valid minimum of data and data.
         const double validMaximum   Valid maximum of data and data.
OUTPUTS: double data[]               Data values read and expanded.
RETURNS: size_t number of unfiltered points.
******************************************************************************/

static size_t readAndExpandData( const int gid,
                                 const int id,
                                 const nc_type type,
                                 const size_t starts[],
                                 const size_t counts[],
                                 const double validMinimum,
                                 const double validMaximum,
                                 double data[] ) {

  size_t result = 0;
  int status = NC_NOERR - 1;
  float* fdata = (float*) data;
  int* idata = (int*) data;
  unsigned int* uidata = (unsigned int*) data;
  short* sdata = (short*) data;
  unsigned short* usdata = (unsigned short*) data;
  signed char* cdata = (signed char*) data;
  unsigned char* ucdata = (unsigned char*) data;
  long long* lldata = (long long*) data;
  unsigned long long* ulldata = (unsigned long long*) data;

  assert( gid >= 0 ); assert( id >= 0 ); assert( type > -1 );
  assert( starts );
  assert( counts );
  assert( validMinimum >= -DBL_MAX );
  assert( validMinimum <=  DBL_MAX );
  assert( validMaximum >= validMinimum );
  assert( validMaximum <= DBL_MAX );
  assert( data );

  if ( type == NC_DOUBLE ) {
    status = nc_get_vara_double( gid, id, starts, counts, data );
  } else if ( type == NC_FLOAT ) {
    status = nc_get_vara_float( gid, id, starts, counts, fdata );
  } else if ( type == NC_INT ) {
    status = nc_get_vara_int( gid, id, starts, counts, idata );
  } else if ( type == NC_UINT ) {
    status = nc_get_vara_uint( gid, id, starts, counts, uidata );
  } else if ( type == NC_SHORT ) {
    status = nc_get_vara_short( gid, id, starts, counts, sdata );
  } else if ( type == NC_USHORT ) {
    status = nc_get_vara_ushort( gid, id, starts, counts, usdata );
  } else if ( type == NC_CHAR ) {
    status = nc_get_vara_schar( gid, id, starts, counts, cdata );
  } else if ( type == NC_BYTE ) {
    status = nc_get_vara_schar( gid, id, starts, counts, cdata );
  } else if ( type == NC_UBYTE ) {
    status = nc_get_vara_uchar( gid, id, starts, counts, ucdata );
  } else if ( type == NC_INT64 ) {
    status = nc_get_vara_longlong( gid, id, starts, counts, lldata );
  } else if ( type == NC_UINT64 ) {
    status = nc_get_vara_ulonglong( gid, id, starts, counts, ulldata );
  }

  if ( status == NC_NOERR - 1 ) {
    fprintf( stderr, "Unsupported data type %d\n", type );
  } else if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read data because: %s\n", message );
  } else {
    size_t points = counts[ 0 ] * counts[ 1 ] * counts[ 2 ];

    /* Loop backwards to expand data values to 64-bit reals: */

    while ( points-- ) {
      double value = -DBL_MAX;

      switch ( type ) {
      case NC_DOUBLE: value =           data[   points ]; break;
      case NC_FLOAT:  value = (double) fdata[   points ]; break;
      case NC_INT:    value = (double) idata[   points ]; break;
      case NC_UINT:   value = (double) uidata[  points ]; break;
      case NC_SHORT:  value = (double) sdata[   points ]; break;
      case NC_USHORT: value = (double) usdata[  points ]; break;
      case NC_CHAR:   value = (double) cdata[   points ]; break;
      case NC_BYTE:   value = (double) cdata[   points ]; break;
      case NC_UBYTE:  value = (double) ucdata[  points ]; break;
      case NC_INT64:  value = (double) lldata[  points ]; break;
      case NC_UINT64: value = (double) ulldata[ points ]; break;
      default: break;
      }

      const int filtered = ! ( value >= validMinimum && value <= validMaximum );

      if ( filtered ) {
        data[ points ] = MISSING_VALUE;
      } else {
        data[ points ] = value;
        ++result;
      }
    }

    DEBUG( fprintf( stderr,
                    "Expanded/filtered %lu points of type %d [%g, %g] "
                    "result = %lu points.\n",
                    counts[ 0 ] * counts[ 1 ] * counts[ 2 ],
                    type, validMinimum, validMaximum, result ); )
  }

  return result;
}



/******************************************************************************
PURPOSE: sumDataAndScratch - Sum data with auxiliary variable.
INPUTS:  const int gid               Group of auxiliary variable.
         const int id                Id of auxiliary variable.
         const size_t starts[]       Subset start 0-based indices to read.
         const size_t counts[]       Subset counts to read.
         const double validMinimum   Valid minimum of data and auxiliary data.
         const double validMaximum   Valid maximum of data and auxiliary data.
         double data[]               Values read.
OUTPUTS: double scratch[]            Auxiliary variable to read into.
         double data[]               Summed data values.
RETURNS: size_t number of unfiltered points.
******************************************************************************/

static size_t sumDataAndScratch( const int gid,
                                 const int id,
                                 const size_t starts[],
                                 const size_t counts[],
                                 const double validMinimum,
                                 const double validMaximum,
                                 double scratch[],
                                 double data[] ) {

  size_t result = 0;
  int status = NC_NOERR - 1;

  assert( gid >= 0 ); assert( id >= 0 );
  assert( starts );
  assert( counts );
  assert( validMinimum >= -DBL_MAX );
  assert( validMinimum <=  DBL_MAX );
  assert( validMaximum >= validMinimum );
  assert( validMaximum <= DBL_MAX );
  assert( scratch );
  assert( data );
  assert( data != scratch );

  status = nc_get_vara_double( gid, id, starts, counts, scratch );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    fprintf( stderr, "Failed to read data to sum because: %s\n",
             message );
  } else {
    const size_t points = counts[ 0 ] * counts[ 1 ] * counts[ 2 ];
    size_t point = 0;

    for ( point = 0; point < points; ++point ) {
      const double value = data[ point ];
      int filtered = ! ( value >= validMinimum && value <= validMaximum );

      if ( filtered ) {
        data[ point ] = MISSING_VALUE;
      } else {
        const double value2 = scratch[ point ];
        filtered = ! ( value2 >= validMinimum && value2 <= validMaximum );

        if ( filtered ) {
          data[ point ] = MISSING_VALUE;
        } else {
          data[ point ] += value2;
          ++result;
        }
      }
    }

    DEBUG( fprintf( stderr,
                    "Summed %lu points filter range [%g, %g] "
                    "result = %lu points.\n",
                    points, validMinimum, validMaximum, result ); )
  }

  return result;
}



/******************************************************************************
PURPOSE: filterDataByQC - Filter data by qc variable.
INPUTS:  const int file              Input file.
         const char* const qcVariable  Name of QC variable to read & filter by.
         const size_t starts[]       Subset start 0-based indices to read.
         const size_t counts[]       Subset counts to read.
         const int qcMinimum         Minimum acceptable qc flag value.
         const int qcMaximum         Maximum acceptable qc flag value.
         const int shift             Right-shift bits for qc flag else 0.
         double data[]               Values read before filtering.
OUTPUTS: double scratch[]            Auxiliary storage to read qc values into.
         double data[]               Filtered data values.
RETURNS: size_t number of unfiltered points.
******************************************************************************/

static size_t filterDataByQC( const int file,
                              const char* const qcVariable,
                              const size_t starts[],
                              const size_t counts[],
                              const int qcMinimum,
                              const int qcMaximum,
                              const int shift,
                              double scratch[],
                              double data[] ) {

  size_t result = 0;
  int rank = 0;
  nc_type type = -1;
  int gid = -1;
  int id = -1;

  assert( file >= 0 );
  assert( qcVariable ); assert( *qcVariable );
  assert( starts ); assert( counts );
  assert( qcMinimum <= qcMaximum );
  assert( scratch ); assert( data ); assert( data != scratch );

  id = findVariable( file, qcVariable, &gid, &rank, &type );

  if ( id < 0 ) {
    fprintf( stderr, "Failed to find variable %s.\n", qcVariable );
  } else {
    int status = NC_NOERR - 1;

    switch ( type ) {
    case NC_SHORT:
    {
      short* qdata = (short*) scratch;
      status = nc_get_vara_short( gid, id, starts, counts, qdata );

      if ( status == NC_NOERR ) {
        const size_t points = counts[ 0 ] * counts[ 1 ] * counts[ 2 ];
        size_t point = 0;

        for ( point = 0; point < points; ++point ) {
          const short qValue = qdata[ point ];
          const int filtered = qValue < qcMinimum || qValue > qcMaximum;

          if ( filtered ) {
            data[ point ] = MISSING_VALUE;
          } else if ( data[ point ] > MISSING_VALUE ) {
            ++result;
          }
        }
      }
    }

    break;
    case NC_USHORT:
    {
      /*
       * TEMPO_Level-2-3_O3TOT_user_guide_V1.0.pdf, pages 13-14.
       * quality_flag in product group reports the detailed quality of the
       * retrieval.
       * Bits 0 to 3 together contain several output error flags:
       * 0: Good sample
       * 1: Glint contamination (corrected)
       * 2: SZA > 84°
       * 3: 360 residual > threshold
       * 4: Residual at unused ozone wavelength > 4 sigma
       * 5: SOI > 4 sigma (SO2 present)
       * 6: Non-convergence
       * 7: Abs (residual) > 16.0 (fatal)
       * 8: Row anomaly error (same as bit 6 in this field)
       * Bits 4 to 5 are reserved for future use (currently set to 0).
       * Bit 7 is set to 0 when TEMPO CLDO4 cloud pressure is used and set to 1 when
       * climatological cloud pressure is used.
       * Bits 8 to 15 are flags that are set to 0 for FALSE (good value), or 1 for TRUE
       * (bad value)
       * Bit 8: Geolocation error (anomalous FOV Earth location)
       * Bit 9: SZA > 88°
       * Bit 10: Missing input radiance
       * Bit 11: Error input radiance
       * Bit 12: Warning input Radiance
       * Bit 13: Missing input irradiance
       * Bit 14: Error input irradiance
       * Bit 15: Warning input irradiance
       */

      unsigned short* qdata = (unsigned short*) scratch;
      status = nc_get_vara_ushort( gid, id, starts, counts, qdata );

      if ( status == NC_NOERR ) {
        const size_t points = counts[ 0 ] * counts[ 1 ] * counts[ 2 ];
        size_t point = 0;

        for ( point = 0; point < points; ++point ) {
          const unsigned short qValue = qdata[ point ];
          const int filtered = qcMaximum == 0 ? qValue != 0 : qValue > 2;

          if ( filtered ) {
            data[ point ] = MISSING_VALUE;
          } else if ( data[ point ] > MISSING_VALUE ) {
            ++result;
          }
        }
      }
    }

    break;
    case NC_UBYTE:
    {
      unsigned char* qdata = (unsigned char*) scratch;
      status = nc_get_vara_uchar( gid, id, starts, counts, qdata );

      if ( status == NC_NOERR ) {
        const size_t points = counts[ 0 ] * counts[ 1 ] * counts[ 2 ];
        size_t point = 0;

        for ( point = 0; point < points; ++point ) {
          const unsigned char qValue = qdata[ point ];
          const int filtered = qValue < qcMinimum || qValue > qcMaximum;

          if ( filtered ) {
            data[ point ] = MISSING_VALUE;
          } else if ( data[ point ] > MISSING_VALUE ) {
            ++result;
          }
        }
      }
    }

    break;
    case NC_BYTE:
    {
      /*
       * TEMPO-ABI_ADP_Users_Guide_V1_20250107.pdf, Table 4 byte qc_flag
       * bits 2-3 = QC_SMOKE_CONFIDENCE
       * bits 4-5 = QC_DUST_CONFIDENCE
       * bits 6-7 = QC_NUC_CONFIDENCE
       * 00 = high, 01 = low, 10 (2) = medium, 11 (3) = bad/missing.
       * To extract the 2-bits, shift right then & 0x3.
       */

      signed char* qdata = (signed char*) scratch;
      status = nc_get_vara_schar( gid, id, starts, counts, qdata );

      if ( status == NC_NOERR ) {
        const size_t points = counts[ 0 ] * counts[ 1 ] * counts[ 2 ];
        size_t point = 0;

        for ( point = 0; point < points; ++point ) {
          const unsigned char qValue = (unsigned char) qdata[ point ];
          const unsigned char qBits = ( qValue >> shift ) & 0x03;
          const int filtered =
            qcMaximum == 0 ? qBits != 0
            : qcMaximum == 1 ? ( qBits != 0 && qBits != 2 )
            : 0;

          if ( filtered ) {
            data[ point ] = MISSING_VALUE;
          } else if ( data[ point ] > MISSING_VALUE ) {
            ++result;
          }
        }
      }
    }

    break;
    default:
    break;
    }

    if ( status == NC_NOERR - 1 ) {
      fprintf( stderr, "Unsupported qc data type %d\n", type );
    } else if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "Failed to read data because: %s\n", message );
    }

    DEBUG( fprintf( stderr,
                    "filterDataByQC %lu points of type %d "
                    "filters [%d, %d] shift = %d "
                    "result = %lu points.\n",
                    counts[ 0 ] * counts[ 1 ] * counts[ 2 ], type,
                    qcMinimum, qcMaximum, shift,
                    result ); )
  }

  return result;
}



/******************************************************************************
PURPOSE: filterDataByAuxiliaryVariable - Filter data by auxiliary variable.
INPUTS:  const int file              File to read.
         const char* const auxiliaryVariable  Variable to read and filter by.
         const size_t starts[]       Subset start 0-based indices to read.
         const size_t counts[]       Subset counts to read.
         const double maximum        Valid maximum of data and auxiliary data.
         double data[]               Values read.
OUTPUTS: double scratch[]            Auxiliary variable to read into.
         double data[]               Summed data values.
RETURNS: size_t number of unfiltered points.
******************************************************************************/

static size_t filterDataByAuxiliaryVariable( const int file,
                                            const char* const auxiliaryVariable,
                                             const size_t starts[],
                                             const size_t counts[],
                                             const double maximum,
                                             double scratch[],
                                             double data[] ) {

  size_t result = 0;
  int rank = 0;
  nc_type type = -1;
  int gid = -1;
  int id = -1;

  assert( file >= 0 );
  assert( auxiliaryVariable ); assert( *auxiliaryVariable );
  assert( starts ); assert( counts );
  assert( maximum >= 0.0 );
  assert( scratch ); assert( data ); assert( data != scratch );

  id = findVariable( file, auxiliaryVariable, &gid, &rank, &type );

  if ( id < 0 ) {
    fprintf( stderr, "Failed to find variable %s.\n", auxiliaryVariable );
  } else {
    double* ddata = 0;
    float* fdata = 0;
    int status = NC_NOERR - 1;

    if ( type == NC_DOUBLE ) {
      ddata = scratch;
      status = nc_get_vara_double( gid, id, starts, counts, ddata );
    } else if ( type == NC_FLOAT ) {
      fdata = (float*) scratch;
      status = nc_get_vara_float( gid, id, starts, counts, fdata );
    }

    if ( status == NC_NOERR - 1 ) {
      fprintf( stderr, "Unsupported aux filter data type %d\n", type );
    } else if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      fprintf( stderr, "Failed to read aux filter data because: %s\n", message);
    } else {
      const size_t points = counts[ 0 ] * counts[ 1 ] * counts[ 2 ];
      size_t point = 0;

      for ( point = 0; point < points; ++point ) {
        const double value = ddata ? ddata[ point ] : (double) fdata[ point ];
        const int filtered = ! ( value >= 0.0 && value <= maximum );

        if ( filtered ) {
          data[ point ] = MISSING_VALUE;
        } else if ( data[ point ] > MISSING_VALUE ) {
          ++result;
        }
      }

      DEBUG( fprintf( stderr,
                      "Filter by %s: %lu points filter range [0, %g] "
                       "result = %lu points.\n",
                       auxiliaryVariable, points, maximum, result ); )
    }
  }

  return result;
}




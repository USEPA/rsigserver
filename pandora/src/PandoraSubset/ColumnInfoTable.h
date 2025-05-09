/******************************************************************************
PURPOSE: ColumnInfoTable.h - Internal header that declares metadata table
         describing Pandora input file columns.

NOTES:   Included by PandoraSubset.c.

HISTORY: 2023-02-01 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/

/*================================ CONSTANTS ================================*/

/* Data points outside these valid ranges are filtered-out: */

#define MINIMUM_VALID_SURFACE_ELEVATION_METERS (-500.0)
#define MAXIMUM_VALID_SURFACE_ELEVATION_METERS 10000.0
#define MAXIMUM_VALID_ELEVATION_METERS 100000.0

/* Data points above surface can be aggregated to the following levels: */

#define AGGREGATION_LEVELS 100
#define AGGREGATION_LEVEL_THICKNESS_METERS 50.0
#define MAXIMUM_AGGREGATION_ELEVATION \
  (AGGREGATION_LEVELS * AGGREGATION_LEVEL_THICKNESS_METERS + \
   MINIMUM_VALID_SURFACE_ELEVATION_METERS)

/* https://en.wikipedia.org/wiki/Dobson_unit */

#define DU_TO_MOLECULES_PER_CM2 2.69e16

/* https://www.convertunits.com/from/mole/to/molecule */

#define MOLECULES_PER_MOL 6.0221415e23
#define M2_PER_CM2 1e-4
#define MOL_PER_M2_TO_MOLECULES_PER_CM2 (MOLECULES_PER_MOL * M2_PER_CM2)

#define KM_TO_M 1e3

/*============================= PRIVATE FUNCTIONS ============================*/

/* Convert bering angle to trig angle: */

static double beringToTrig(double bering) {
  double result = 450.0 - bering;
  while (result >= 360.0) result -= 360.0;
  return result;
}

/*================================== TYPES ==================================*/

/*
 *
 * List of known Pandora input file variables: variable, units and column.
 * variable is the name of the variable in pandoraserver
 * units is the output units of the variable in pandoraserver
 * column is the 1-based index of the variable on the data line.
 *
 * Note: the metadata below is based on Pandora files on maple as of 2023-01-01.
 * If new versions of Pandora files are added to maple then this table below
 * will need to be extended (as well as pandoraserver).
 * See NOTES file for excerpts of known Pandora file types handled by this code.
 */

typedef struct {
  const char* type;     /* Type of file. E.g., "_L2Tot_rnvs0p1-7.txt".    */
  const char* name;     /* Name of output variable. E.g., "temperature".  */
  const double offset;  /* Value to add to original value. E.g., -32.0.   */
  const double scale;   /* Value to multiply the offset value. E.g., 5.0/9/0 */
  double (*converter)(double); /* Convert function. new = converter(original)*/
  const char* units;    /* Units of output variable. E.g., "C".          */
  const double minimum; /* Minimum valid data value. E.g., -50.0.        */
  const double maximum; /* Minimum valid data value. E.g., 50.0.         */
  const int index;           /* 1-based column number of name. E.g., 20. */
  const int indexStride;     /* Stride to next level measure or 0.       */
  const int elevationIndex;  /* 1-based column number of elevation or 0. */
  const int elevationStride; /* Stride to next elevation or 0.           */
  const int filterIndex;     /* 1-based column number of elevation or 0. */
} ColumnInfo;

static const ColumnInfo columnInfoTable[] = {

  /* _L2Tot_rnvs0p1-7.txt ---------------------------------------------------*/

  { "_L2Tot_rnvs0p1-7.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "nitrogen_dioxide_total_vertical_column_amount",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    8, 0, 0, 0, 12 },

  { "_L2Tot_rnvs0p1-7.txt",
    "nitrogen_dioxide_total_vertical_column_amount_uncertainty",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    9, 0, 0, 0, 12 },

  { "_L2Tot_rnvs0p1-7.txt", "direct_nitrogen_dioxide_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 12 },

  { "_L2Tot_rnvs0p1-7.txt", "diffuse_correction_for_nitrogen_dioxide",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 22.0,
    12, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    13, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    14, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "fitting_result_index",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    15, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    16, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt",
    "expected_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    17, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt",
    "expected_normalized_rms_weighted_instrumental_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    18, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt",
    "mean_nitrogen_dioxide_total_vertical_column_amount",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    19, 0, 0, 0, 12 },

  { "_L2Tot_rnvs0p1-7.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    20, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    21, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    22, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    23, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "l2_fit_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    24, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "sum_of_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    25, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "sum_of_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    26, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    27, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    28, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    29, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    30, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    31, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    32, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "wavelength_shift_from_spectral_fitting",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    33, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    34, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "number_of_dark_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    35, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    36, 0, 0, 0, 0 },

  { "_L2Tot_rnvs0p1-7.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    37, 0, 0, 0, 0 },


  /* _L2Tot_rnvs1p1-7.txt ---------------------------------------------------*/

  { "_L2Tot_rnvs1p1-7.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "nitrogen_dioxide_total_vertical_column_amount",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    8, 0, 0, 0, 12 },

  { "_L2Tot_rnvs1p1-7.txt",
    "nitrogen_dioxide_total_vertical_column_amount_uncertainty",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    9, 0, 0, 0, 12 },

  { "_L2Tot_rnvs1p1-7.txt", "direct_nitrogen_dioxide_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 12 },

  { "_L2Tot_rnvs1p1-7.txt", "diffuse_correction_for_nitrogen_dioxide",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 22.0,
    12, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    13, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    14, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "fitting_result_index",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    15, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    16, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt",
    "expected_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    17, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt",
    "expected_normalized_rms_weighted_instrumental_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    18, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt",
    "mean_nitrogen_dioxide_total_vertical_column_amount",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    19, 0, 0, 0, 12 },

  { "_L2Tot_rnvs1p1-7.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    20, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    21, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    22, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    23, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "l2_fit_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    24, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "sum_of_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    25, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "sum_of_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    26, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    27, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    28, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    29, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    30, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    31, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    32, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "wavelength_shift_from_spectral_fitting",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    33, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    34, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "number_of_dark_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    35, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    36, 0, 0, 0, 0 },

  { "_L2Tot_rnvs1p1-7.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    37, 0, 0, 0, 0 },


  /* _L2Tot_rout0p1-7.txt ---------------------------------------------------*/

  { "_L2Tot_rout0p1-7.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "ozone_total_vertical_column_amount",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    8, 0, 0, 0, 12 },

  { "_L2Tot_rout0p1-7.txt",
    "ozone_total_vertical_column_amount_uncertainty",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    9, 0, 0, 0, 12 },

  { "_L2Tot_rout0p1-7.txt", "direct_ozone_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 12 },

  { "_L2Tot_rout0p1-7.txt", "diffuse_correction_for_ozone",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "ozone_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 22.0,
    12, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_ozone_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    13, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_ozone_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    14, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sulfur_dioxide_total_vertical_column_amount",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    15, 0, 0, 0, 19 },

  { "_L2Tot_rout0p1-7.txt",
    "sulfur_dioxide_total_vertical_column_amount_uncertainty",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    16, 0, 0, 0, 19 },

  { "_L2Tot_rout0p1-7.txt", "direct_sulfur_dioxide_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    17, 0, 0, 0, 19 },

  { "_L2Tot_rout0p1-7.txt", "diffuse_correction_for_sulfur_dioxide",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    18, 0, 0, 0, 19 },

  { "_L2Tot_rout0p1-7.txt", "sulfur_dioxide_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 22.0,
    19, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_sulfur_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    20, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_sulfur_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    21, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "fitting_result_index",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    22, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    23, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt",
    "expected_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    24, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt",
    "expected_normalized_rms_weighted_instrumental_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    25, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "mean_sulfur_dioxide_total_vertical_column_amount",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    26, 0, 0, 0, 19 },

  { "_L2Tot_rout0p1-7.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    27, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    28, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    29, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    30, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "l2_fit_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    31, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    32, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    33, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    34, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    35, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    36, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    37, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    38, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    39, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "wavelength_shift_from_spectral_fitting",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    40, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    41, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "number_of_dark_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    42, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    43, 0, 0, 0, 0 },

  { "_L2Tot_rout0p1-7.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    44, 0, 0, 0, 0 },


  /* _L2_rfuh5p1-8.txt ------------------------------------------------------*/

  { "_L2_rfuh5p1-8.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "pointing_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    8, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "pointing_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    9, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "normalized_rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "expected_rms_unweighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    12, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "expected_rms_weighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    13, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    14, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    15, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "o2_height",
    0.0, KM_TO_M, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    16, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "o2o2_height",
    0.0, KM_TO_M, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    17, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "surface_o2",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e3,
    18, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "surface_o2o2",
    0.0, 1.0, 0, "mol2/m6", 0.0, 1e3,
    19, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "total_o2_column",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    20, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "total_o2o2_column",
    0.0, 1.0, 0, "mol2/m5", 0.0, 1e7,
    21, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    22, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    23, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    24, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "mean_formaldehyde_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    25, 0, 0, 0, 42 },

  { "_L2_rfuh5p1-8.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    26, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    27, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    28, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "total_wavelength_shift",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    29, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "resolution_change",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    30, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    31, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "number_of_bright_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    32, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    33, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    34, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "atmospheric_variability",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    35, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    36, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    37, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    38, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "l2_fit_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    39, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    40, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    41, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "formaldehyde_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    42, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "sum_of_formaldehyde_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    43, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "sum_of_formaldehyde_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    44, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "surface_formaldehyde",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e6,
    45, 0, 0, 0, 42 },

  { "_L2_rfuh5p1-8.txt", "surface_formaldehyde_uncertainty",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e6,
    46, 0, 0, 0, 42 },

  { "_L2_rfuh5p1-8.txt", "surface_formaldehyde_index",
    0.0, 1.0, 0, "-", 0.0, 1e6,
    47, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "formaldehyde_heterogeneity_flag",
    0.0, 1.0, 0, "-", 0.0, 1.0,
    48, 0, 0, 0, 0 },

  { "_L2_rfuh5p1-8.txt", "formaldehyde_tropospheric_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    49, 0, 0, 0, 42 },

  { "_L2_rfuh5p1-8.txt",
    "formaldehyde_tropospheric_vertical_column_amount_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    50, 0, 0, 0, 42 },

  { "_L2_rfuh5p1-8.txt",
    "maximum_horizontal_distance_formaldehyde_tropospheric_column",
    0.0, KM_TO_M, 0, "m", 0.0, 1e8,
    51, 0, 0, 0, 42 },

  { "_L2_rfuh5p1-8.txt",
    "maximum_vertical_distance_formaldehyde_tropospheric_column",
    0.0, KM_TO_M, 0, "m", 0.0, 1e8,
    52, 0, 0, 0, 42 },

  /* Per-layer measures: */

  { "_L2_rfuh5p1-8.txt", "formaldehyde_layer_top_height",
    0.0, KM_TO_M, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    53, 2, 0, 0, 42 },

  { "_L2_rfuh5p1-8.txt", "formaldehyde_layer_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    54, 2, 53, 2, 42 },

  /*
   * Optional additional pairs of layer top height and layer amount (mol/m2)
   * may follow the above layer1 pair of values.
   */


  /* _L2_rfus5p1-8.txt ------------------------------------------------------*/

  { "_L2_rfus5p1-8.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    8, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "normalized_rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    9, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "expected_rms_unweighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "expected_normalized_rms_weighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    12, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    13, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    14, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    15, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt",
    "mean_formaldehyde_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    16, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    17, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    18, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    19, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "total_wavelength_shift",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    20, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "resolution_change",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    21, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    22, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "number_of_bright_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    23, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    24, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    25, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "atmospheric_variability",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    26, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "aerosol_optical_depth_start",
    0.0, 1.0, 0, "-", -1.0, 100.0,
    27, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "aerosol_optical_depth_center",
    0.0, 1.0, 0, "-", -1.0, 100.0,
    28, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "aerosol_optical_depth_end",
    0.0, 1.0, 0, "-", -1.0, 100.0,
    29, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    30, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    31, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    32, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "l2_fit_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    33, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    34, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    35, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "fomaldehyde_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 22.0,
    36, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "sum_of_formaldehyde_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    37, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "sum_of_formaldehyde_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    38, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "formaldehyde_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    39, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt",
    "formaldehyde_vertical_column_amount_independent_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    40, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt",
    "formaldehyde_vertical_column_amount_structured_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    41, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt",
    "formaldehyde_vertical_column_amount_common_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    42, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt", "formaldehyde_vertical_column_amount_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    43, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt", "formaldehyde_vertical_column_amount_rms_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    44, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    45, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "temperature_independent_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    46, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "temperature_structured_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    47, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "temperature_common_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    48, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "temperature_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    49, 0, 0, 0, 0 },

  { "_L2_rfus5p1-8.txt", "direct_formaldehyde_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    50, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt", "direct_formaldehyde_air_mass_factor_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    51, 0, 0, 0, 36 },

  { "_L2_rfus5p1-8.txt", "diffuse_correction_for_formaldehyde",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    52, 0, 0, 0, 36 },


  /* _L2_rnvh3p1-8.txt ------------------------------------------------------*/

  { "_L2_rnvh3p1-8.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "pointing_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    8, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "pointing_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    9, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "normalized_rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "expected_rms_unweighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    12, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "expected_normalized_rms_weighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    13, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    14, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    15, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "o2_height",
    0.0, KM_TO_M, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    16, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "o2o2_height",
    0.0, KM_TO_M, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    17, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "surface_o2",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e3,
    18, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "surface_o2o2",
    0.0, 1.0, 0, "mol2/m6", 0.0, 1e3,
    19, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "total_o2_column",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    20, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "total_o2o2_column",
    0.0, 1.0, 0, "mol2/m5", 0.0, 1e7,
    21, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    22, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    23, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    24, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "mean_formaldehyde_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    25, 0, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    26, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    27, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    28, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "total_wavelength_shift",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    29, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "resolution_change",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    30, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    31, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "number_of_bright_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    32, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    33, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    34, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "atmospheric_variability",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    35, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    36, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    37, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    38, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "l2_fit_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    39, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    40, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    41, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "water_vapor_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    42, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_water_vapor_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    43, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_water_vapor_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    44, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "surface_water_vapor",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e6,
    45, 0, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt", "surface_water_vapor_uncertainty",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e6,
    46, 0, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt", "surface_water_vapor_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    47, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "water_vapor_heterogeneity_flag",
    0.0, 1.0, 0, "-", 0.0, 1.0,
    48, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "water_vapor_tropospheric_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    49, 0, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt",
    "water_vapor_tropospheric_vertical_column_amount_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    50, 0, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt",
    "maximum_horizontal_distance_water_vapor_tropospheric_column",
    0.0, KM_TO_M, 0, "m", 0.0, 1e8,
    51, 0, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt",
    "maximum_vertical_distance_water_vapor_tropospheric_column",
    0.0, KM_TO_M, 0, "m", 0.0, 1e8,
    52, 0, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt", "nitrogen_dioxide_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    53, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    54, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    55, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "surface_nitrogen_dioxide",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e6,
    56, 0, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt", "surface_nitrogen_dioxide_uncertainty",
    0.0, 1.0, 0, "mol/m3", 0.0, 1e6,
    57, 0, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt", "surface_nitrogen_dioxide_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    58, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "nitrogen_dioxide_heterogeneity_flag",
    0.0, 1.0, 0, "-", 0.0, 1.0,
    59, 0, 0, 0, 0 },

  { "_L2_rnvh3p1-8.txt", "stratospheric_nitrogen_dioxide",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    60, 0, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt", "stratospheric_nitrogen_dioxide_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    61, 0, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt", "tropospheric_nitrogen_dioxide",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    62, 0, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt", "tropospheric_nitrogen_dioxide_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    63, 0, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt",
    "maximum_horizontal_distance_nitrogen_dioxide_tropospheric_column",
    0.0, KM_TO_M, 0, "m", 0.0, 1e8,
    64, 0, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt",
    "maximum_vertical_distance_nitrogen_dioxide_tropospheric_column",
    0.0, KM_TO_M, 0, "m", 0.0, 1e8,
    65, 0, 0, 0, 53 },

  /* Per-layer measures: */

  { "_L2_rnvh3p1-8.txt", "water_vapor_layer_top_height",
    0.0, KM_TO_M, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    66, 4, 0, 0, 42 },

  { "_L2_rnvh3p1-8.txt", "water_vapor_layer_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    67, 4, 66, 4, 42 },

  { "_L2_rnvh3p1-8.txt", "nitrogen_dioxide_layer_top_height",
    0.0, KM_TO_M, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    68, 4, 0, 0, 53 },

  { "_L2_rnvh3p1-8.txt", "nitrogen_dioxide_layer_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    69, 4, 68, 4, 53 },

  /*
   * Optional additional 4 columns per layer
   * top height and layer amount (mol/m2)
   * may follow the above layer1 pair of values.
   */


  /* _L2_rnvs3p1-8.txt ------------------------------------------------------*/

  { "_L2_rnvs3p1-8.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    8, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "normalized_rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    9, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "expected_rms_unweighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "expected_normalized_rms_weighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    12, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    13, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    14, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    15, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "mean_nitrogen_dioxide_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    16, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    17, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    18, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    19, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "total_wavelength_shift",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    20, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "resolution_change",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    21, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    22, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "number_of_bright_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    23, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    24, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    25, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "atmospheric_variability",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    26, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "aerosol_optical_depth_start",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    27, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "aerosol_optical_depth_center",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    28, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "aerosol_optical_depth_end",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    29, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    30, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    31, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    32, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "l2_fit_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    33, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    34, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    35, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "nitrogen_dioxide_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    36, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    37, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    38, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "nitrogen_dioxide_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    39, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt",
    "nitrogen_dioxide_vertical_column_independent_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    40, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt",
    "nitrogen_dioxide_vertical_column_structured_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    41, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt",
    "nitrogen_dioxide_vertical_column_common_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    42, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt", "nitrogen_dioxide_vertical_column_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    43, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt", "nitrogen_dioxide_vertical_rms_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    44, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    45, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "temperature_independent_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    46, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "temperature_structured_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    47, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "temperature_common_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    48, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "temperature_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    49, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "direct_nitrogen_dioxide_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    50, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt", "direct_nitrogen_dioxide_air_mass_factor_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    51, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt", "diffuse_correction",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    52, 0, 0, 0, 0 },

  { "_L2_rnvs3p1-8.txt", "stratospheric_nitrogen_dioxide",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    53, 0, 0, 0, 36 },

  { "_L2_rnvs3p1-8.txt", "stratospheric_nitrogen_dioxide_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    54, 0, 0, 0, 36 },


  /* _L2_rout2p1-8.txt ------------------------------------------------------*/

  { "_L2_rout2p1-8.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    8, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "normalized_rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    9, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "expected_rms_unweighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "expected_normalized_rms_weighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    12, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    13, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    14, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    15, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "mean_ozone_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    16, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    17, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    18, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    19, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "total_wavelength_shift",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    20, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "resolution_change",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    21, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    22, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "number_of_bright_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    23, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    24, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    25, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "atmospheric_variability",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    26, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "aerosol_optical_depth_start",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    27, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "aerosol_optical_depth_center",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    28, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "aerosol_optical_depth_end",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    29, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    30, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    31, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    32, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "l2_fit_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    33, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    34, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    35, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "ozone_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    36, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "sum_of_ozone_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    37, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "sum_of_ozone_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    38, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "ozone_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    39, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "ozone_vertical_column_independent_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    40, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "ozone_vertical_column_structured_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    41, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "ozone_vertical_column_common_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    42, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "ozone_vertical_column_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    43, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "ozone_vertical_column_rms_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    44, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    45, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "temperature_independent_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    46, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "temperature_structured_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    47, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "temperature_common_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    48, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "temperature_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    49, 0, 0, 0, 0 },

  { "_L2_rout2p1-8.txt", "direct_ozone_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    50, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "ozone_air_mass_factor_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    51, 0, 0, 0, 36 },

  { "_L2_rout2p1-8.txt", "diffuse_correction",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    52, 0, 0, 0, 0 },


  /* _L2_rsus1p1-8.txt ------------------------------------------------------*/

  { "_L2_rsus1p1-8.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    8, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "normalized_rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    9, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "expected_rms_unweighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt",
    "expected_normalized_rms_weighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    12, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    13, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    14, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    15, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "mean_sulfur_dioxide_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    16, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    17, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    18, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    19, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "total_wavelength_shift",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    20, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "resolution_change",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    21, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    22, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "number_of_bright_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    23, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    24, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    25, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "atmospheric_variability",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    26, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "aerosol_optical_depth_start",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    27, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "aerosol_optical_depth_center",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    28, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "aerosol_optical_depth_end",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    29, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    30, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    31, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    32, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "l2_fit_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    33, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    34, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sum_of_l2_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    35, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sulfur_dioxide_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    36, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sum_of_sulfur_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    37, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sum_of_sulfur_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    38, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "sulfur_dioxide_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    39, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt",
    "sulfur_dioxide_vertical_column_independent_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    40, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt",
    "sulfur_dioxide_vertical_column_structured_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    41, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt", "sulfur_dioxide_vertical_column_common_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    42, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt", "sulfur_dioxide_vertical_column_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    43, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt", "sulfur_dioxide_vertical_column_rms_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    44, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    45, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "temperature_independent_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    46, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "temperature_structured_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    47, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "temperature_common_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    48, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "temperature_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    49, 0, 0, 0, 0 },

  { "_L2_rsus1p1-8.txt", "direct_sulfur_dioxide_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    50, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt", "sulfur_dioxide_air_mass_factor_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    51, 0, 0, 0, 36 },

  { "_L2_rsus1p1-8.txt", "diffuse_correction",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    52, 0, 0, 0, 0 },


  /* _L2Trop_rnvh1p1-7.txt ------------------------------------------------------*/


  { "_L2Trop_rnvh1p1-7.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "pointing_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "pointing_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "water_vapor_surface",
    0.0, 1.0, 0, "ppb", 0.0, 1e9,
    8, 0, 0, 0, 14 },

  { "_L2Trop_rnvh1p1-7.txt", "water_vapor_surface_uncertainty",
    0.0, 1.0, 0, "ppb", 0.0, 1e9,
    9, 0, 0, 0, 14 },

  { "_L2Trop_rnvh1p1-7.txt", "water_vapor_tropospheric",
    0.0, 1.0, 0, "cm", 0.0, 1e9,
    10, 0, 0, 0, 14 },

  { "_L2Trop_rnvh1p1-7.txt", "water_vapor_tropospheric_uncertainty",
    0.0, 1.0, 0, "cm", 0.0, 1e9,
    11, 0, 0, 0, 14 },

  { "_L2Trop_rnvh1p1-7.txt", "water_vapor_surface_index",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    12, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "water_vapor_heterogeneity_flag",
    0.0, 1.0, 0, "-", 0.0, 1.0,
    13, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "water_vapor_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 22.0,
    14, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_water_vapor_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    15, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_water_vapor_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    16, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "direct_water_vapor_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    17, 0, 0, 0, 14 },

  { "_L2Trop_rnvh1p1-7.txt", "nitrogen_dioxide_surface",
    0.0, 1.0, 0, "ppb", 0.0, 1e3,
    18, 0, 0, 0, 24 },

  { "_L2Trop_rnvh1p1-7.txt", "nitrogen_dioxide_surface_uncertainty",
    0.0, 1.0, 0, "ppb", 0.0, 1e3,
    19, 0, 0, 0, 24 },

  { "_L2Trop_rnvh1p1-7.txt", "nitrogen_dioxide_tropospheric",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    20, 0, 0, 0, 24 },

  { "_L2Trop_rnvh1p1-7.txt", "nitrogen_dioxide_tropospheric_uncertainty",
    0.0, DU_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    21, 0, 0, 0, 24 },

  { "_L2Trop_rnvh1p1-7.txt", "nitrogen_dioxide_surface_index",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    22, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "nitrogen_dioxide_heterogeneity_flag",
    0.0, 1.0, 0, "-", 0.0, 1.0,
    23, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "nitrogen_dioxide_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 22.0,
    24, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    25, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    26, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "direct_nitrogen_dioxide_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    27, 0, 0, 0, 24 },

  { "_L2Trop_rnvh1p1-7.txt", "fitting_result_index",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    28, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "mean_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    29, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "minimum_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    30, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "maximum_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    31, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "mean_expected_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    32, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "minimum_expected_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    33, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "maximum_expected_normalized_rms_weighted_measured_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    34, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "mean_expected_normalized_rms_weighted_instrumental_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    35, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "minimum_expected_normalized_rms_weighted_instrumental_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    36, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt",
    "maximum_expected_normalized_rms_weighted_instrumental_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 1000.0,
    37, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    38, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    39, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "o2o2_height",
    0.0, 1.0, 0, "m", 0.0, MAXIMUM_VALID_ELEVATION_METERS,
    40, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    41, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    42, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    43, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_l2_fit_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    44, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_l2_fit_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    45, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    46, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    47, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    48, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    49, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    50, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    51, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    52, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "mean_residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    53, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    54, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    55, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "mean_wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    56, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    57, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    58, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "mean_wavelength_shift_from_spectral_fitting",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    59, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_wavelength_shift_from_spectral_fitting",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    60, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_wavelength_shift_from_spectral_fitting",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    61, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    62, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    63, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "mean_number_of_dark_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    64, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_number_of_dark_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    65, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_number_of_dark_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    66, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    67, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    68, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "minimum_filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    69, 0, 0, 0, 0 },

  { "_L2Trop_rnvh1p1-7.txt", "maximum_filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    70, 0, 0, 0, 0 },


  /* _L2_rnvssp1-8.txt ------------------------------------------------------*/

  { "_L2_rnvssp1-8.txt", "timestamp",
    0.0, 1.0, 0, "-", 20000101000000.0, 21001231235959.0,
    1, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "days_since_2000",
    0.0, 1.0, 0, "days", 0.0, 36600.0,
    2, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "measurement_duration",
    0.0, 1.0, 0, "s", 1e-3, 300.0,
    3, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "solar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    4, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "solar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    5, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "lunar_zenith_angle",
    0.0, 1.0, 0, "deg", 0.0, 180.0,
    6, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "lunar_azimuth_angle",
    0.0, 1.0, beringToTrig, "deg", -360.0, 360.0,
    7, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    8, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "normalized_rms_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    9, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "expected_rms_unweighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    10, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "expected_normalized_rms_weighted_fitting_residuals",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    11, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "pressure",
    0.0, 1.0, 0, "hPa", 500.0, 1500.0,
    12, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "data_processing_type_index",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    13, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "calibration_file_version",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    14, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "calibration_file_date",
    0.0, 1.0, 0, "yyyymmdd", 20000101.0, 21001231.0,
    15, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "mean_nitrogen_dioxide_total_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    16, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt", "wavelength_temperature",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    17, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "residual_stray_light_level",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    18, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "wavelength_shift_from_l1",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    19, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "total_wavelength_shift",
    0.0, 1.0, 0, "nm", -10.0, 10.0,
    20, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "resolution_change",
    0.0, 1.0, 0, "%", 0.0, 100.0,
    21, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "integration_time",
    0.0, 1.0, 0, "ms", 1e-3, 1e4,
    22, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "number_of_bright_count_cycles",
    0.0, 1.0, 0, "-", 0.0, 1e4,
    23, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "filterwheel1_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    24, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "filterwheel2_position",
    0.0, 1.0, 0, "-", 1.0, 9.0,
    25, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "atmospheric_variability",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    26, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "aerosol_optical_depth_start",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    27, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "aerosol_optical_depth_center",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    28, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "aerosol_optical_depth_end",
    0.0, 1.0, 0, "-", 0.0, 10.0,
    29, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "l1_data_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    30, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "sum_of_l1_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    31, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "sum_of_l1_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    32, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "l2_fit_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    33, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "sum_of_l2_fit_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    34, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "sum_of_l2_fit_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    35, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "nitrogen_dioxide_l2_quality_flag",
    0.0, 1.0, 0, "-", 0.0, 12.0,
    36, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq1",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    37, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "sum_of_nitrogen_dioxide_quality_exceeding_dq2",
    0.0, 1.0, 0, "-", 0.0, 256.0,
    38, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "nitrogen_dioxide_vertical_column_amount",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    39, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt",
    "nitrogen_dioxide_vertical_column_independent_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    40, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt",
    "nitrogen_dioxide_vertical_column_structured_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    41, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt",
    "nitrogen_dioxide_vertical_column_common_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    42, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt", "nitrogen_dioxide_vertical_column_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    43, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt", "nitrogen_dioxide_vertical_rms_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    44, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt", "temperature",
    -273.15, 1.0, 0, "C", -50.0, 50.0,
    45, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "temperature_independent_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    46, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "temperature_structured_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    47, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "temperature_common_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    48, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "temperature_uncertainty",
    0.0, 1.0, 0, "C", -50.0, 50.0,
    49, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "direct_nitrogen_dioxide_air_mass_factor",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    50, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt", "direct_nitrogen_dioxide_air_mass_factor_uncertainty",
    0.0, 1.0, 0, "-", 0.0, 100.0,
    51, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt", "diffuse_correction",
    0.0, 1.0, 0, "%", -100.0, 100.0,
    52, 0, 0, 0, 0 },

  { "_L2_rnvssp1-8.txt", "stratospheric_nitrogen_dioxide",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    53, 0, 0, 0, 36 },

  { "_L2_rnvssp1-8.txt", "stratospheric_nitrogen_dioxide_uncertainty",
    0.0, MOL_PER_M2_TO_MOLECULES_PER_CM2, 0, "molecules/cm2", 0.0, 1e30,
    54, 0, 0, 0, 36 },


  { "_", "_", 0.0, 0.0, 0, "_", 0.0, 0.0, 0, 0, 0, 0, 0 }, /* End of table. */

};



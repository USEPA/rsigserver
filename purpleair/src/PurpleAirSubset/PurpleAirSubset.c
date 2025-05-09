/******************************************************************************
PURPOSE: PurpleAirSubset.c - Extract a lon-lat subset of data from a list of
         PurpleAir files and write it to stdout as XDR binary format.

NOTES:   This program calls system( "/usr/bin/sort ..." ).

         Compile:
         gcc -Wall -g -o PurpleAirSubset PurpleAirSubset.c Utilities.c \
                   -L../../../lib/$platform \
                   -lm -lc

         Usage:
         PurpleAirSubset \
           -files <listfile> \
           -tmpdir <temp_directory> \
           -desc "description text" \
           -timerange <yyyymmddhhmmss> <yyyymmddhhmmss> \
           -variable <name> \
           -bounds <minimum_longitude> <minimum_latitude> \
                   <maximum_longitude> <maximum_latitude> \
           -format ascii | xdr \
           [-sensor sensor_id] \
           [-out_in_flag 0|1] (0 = outside, 1 = inside, default is either)
           [-default_humidity 0-100] (value for missing/invalid humidity,
                                      default none)
           [-maximum_difference >= 0 (ug/m3) (default 5)
           [-maximum_ratio 0.0-1.0 (default 0.7)
           [-aggregate hourly | daily | all] (default is none)
           [-minimum_aggregation_count_percentage 0-100] (default 75)

          Example:
          ../../../bin/$platform/PurpleAirSubset \
          -files testdata/file_list \
          -tmpdir testdata \
          -desc "https://api.purpleair.com,PurpleAirSubset" \
          -timerange 20201202000000 20201202235959 \
          -variable pm25_corrected \
          -bounds -124 49 -123 50 \
          -format xdr \
          -out_in_flag 0 \
          -aggregate hourly \
          > subset.xdr

          Outputs to stdout a data stream in the following format -
          a 10-line ASCII header followed by binary 64-bit big-endian arrays:

Point 1.0
https://api.purpleair.com,PurpleAirSubset
2020-12-02T00:00:00-0000 2020-12-02T23:59:59-0000
# Dimensions: variables points
6 459
# Variable names:
timestamp longitude latitude elevation id pm25
# Variable units:
yyyymmddhhmmss deg deg m - ug/m3
# char notes[points][80] and
# IEEE-754 64-bit reals data[variables][points]:
<big-endian binary format array>

HISTORY: 2020-12-02 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/


/*================================ INCLUDES =================================*/

#include <assert.h>    /* For assert(). */
#include <stdio.h>     /* For FILE, printf(), snprintf(). */
#include <string.h>    /* For memset(), strcmp(), strstr(), strtok_r(). */
#include <ctype.h>     /* For isalpha(). */
#include <stdlib.h>    /* For malloc(), free(), atoll(), atof(), system(). */
#include <limits.h>    /* For INT_MAX. */
#include <math.h>      /* For exp(). */
#include <unistd.h>    /* For unlink(), getpid() */

#include "Utilities.h" /* For LONGITUDE, Bounds, readFile(). */

/*================================= MACROS =================================*/

/*
 * If USE_SIGMOID is 0 use the piecewise-formula for pm25_atm_a/b
 * else use the sigmoid-formula for pm25_cf_1_a/b.
 */

#define USE_SIGMOID 0

/* Name of temp per-variable files created in -tmpdir with PID appended: */

#define TEMP_FILE_NAME "junk_PurpleAirSubset"

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

/*================================ CONSTANTS ================================*/

/*
 * On or before May 30, 2019 sensors report a set of values every 80 seconds.
 * After that date they report a set of values every 120 seconds.
 */

enum {
  YYYYMMDD_PREVIOUS = 20190530,
  PREVIOUS_SECONDS_PER_VALUE = 80,
  SECONDS_PER_VALUE = 120
};

/*
 * If aggregating, the following variable will be initialized just once to
 * seconds of YYYYMMDD_PREVIOUS.
 */

static double previousSeconds = 0.0;

#define MISSING_VALUE (-999.0)

/* Data points outside these valid ranges are filtered-out: */

#define MINIMUM_VALID_HUMIDITY 0.0
#define MAXIMUM_VALID_HUMIDITY 100.0

#define MINIMUM_VALID_PM 0.0
#define MAXIMUM_VALID_PM 1e6

#define MINIMUM_VALID_ELEVATION_METERS (-500.0)
#define MAXIMUM_VALID_ELEVATION_METERS 10000.0

#define FEET_TO_METERS 0.3048
#define METERS_TO_FEET 3.2808

/*
 * By default, when computing pm25_corrected, accept either
 * 1. absolute difference <= 5 ug/m3 of channel A & B measures or
 * 2. absolute ratio <= 0.7.
 */

static const double defaultMaximumDifference = 5.0;
static const double defaultMaximumRatio      = 0.7;

/*
 * Percentage of maximum number of available values per sensor,
 * per aggregation period. E.g., if aggregating hourly and a sensor reports a
 * value every 2 minutes then the maximum number of values from
 * that sensor would be 60 / 2 = 30. So 75% of 30 = 22.5. Rounding up = 23.
 * So omit that sensor from that hour if it reports less than 23 valid values.
 */

static const double defaultMinimumAggregationCountPercentage = 75.0;

/*
 * channel_state:
 * 0 = No channel is available
 * 1 = Only channel A is available
 * 2 = Only channel B is available
 * 3 = Both channel A and B are available
 *
 * Do we apply channel_state to channel-dependent variables?
 * Examples:
 * pm10 depends on both channel A and B so filter-out unless channel_state = 3.
 * pm25_a depends on channel A only so filter-out unless channel_state = 1 or 3.
 * humidity does not depend on either channel so channel_state is not applied.
 */

static const int applyChannelState = 0; /* 1 = apply, 0 = ignore. */

/*
 * channel_flag:
 * 0 = Normal
 * 1 = Channel A is degraded
 * 2 = Channel B is degraded
 * 3 = Both channel A and B are degraded
 *
 * Do we apply channel_flag to channel-dependent variables?
 * Examples:
 * pm10 depends on both channel A and B so filter-out unless channel_flag = 0.
 * pm25_a depends on channel A only so filter-out unless channel_flag = 0 or 2.
 * humidity does not depend on either channel so channel_flag is not applied.
 */

static const int applyChannelFlag = 0; /* 1 = apply, 0 = ignore. */


/*
 * After valid_minimum/maximum range filtering and possible channel filtering
 * do we output possible negative (presumably small) pm25_corrected values?
 */

static const int allowNegativePM25Corrected = 1; /* 1 = allow, 0 = filter-out*/


#if USE_SIGMOID

/******************************************************************************
PURPOSE: pm25_corrected_sigmoid - Compute pm25_corrected from pm25_cf_1_a/b
         channel measures and relative humidity with filtering by
         absolute difference and ratio in channel measures.
INPUTS:  const double pm25_cf_1_a   Channel A PM2.5 measure (ug/m3).
         const double pm25_cf_1_b   Channel B PM2.5 measure (ug/m3).
         const double humidity      Relative humidity measure (%).
         const double maximum_channel_difference
                                    Maximum acceptable absolute difference of
                                    channel measures (ug/m3).
         const double maximum_channel_ratio
                                    Maximum acceptable absolute ratio of
                                    channel measures [0.0, 1.0].
RETURNS: double pm25_corrected or returns MISSING_VALUE if invalid values are
         given or generated.
NOTES: This version is used for data with channel A & B PM25 both available.
******************************************************************************/

static double pm25_corrected_sigmoid(
  const double pm25_cf_1_a,
  const double pm25_cf_1_b,
  const double humidity,
  const double maximum_channel_difference,
  const double maximum_channel_ratio ) {

  double pm25_corrected = MISSING_VALUE; /* Default result. */

  /* Formula constants: */

#define SIGMOID_PM25_MEAN_OFFSET (-343.0)

#define LINEAR_TERM_PM25_COEFFICIENT       0.52
#define LINEAR_TERM_HUMIDITY_COEFFICIENT (-0.086)
#define LINEAR_TERM_OFFSET                 5.75

#define QUADRATIC_TERM_PM25_SQUARED_COEFFICIENT  3.93e-4
#define QUADRATIC_TERM_PM25_COEFFICIENT          0.46
#define QUADRATIC_TERM_OFFSET                    2.97

  /* Check if humidity is usable: */

  int usable =
    IN_RANGE( humidity, MINIMUM_VALID_HUMIDITY, MAXIMUM_VALID_HUMIDITY );

  assert( IN_RANGE( maximum_channel_difference, 0.0, 100.0 ) );
  assert( IN_RANGE( maximum_channel_ratio, 0.0, 1.0 ) );

  if ( usable ) {

    /* Check if a and b channel measures are usable: */

    usable =
      IN_RANGE( pm25_cf_1_a, MINIMUM_VALID_PM, MAXIMUM_VALID_PM ) &&
      IN_RANGE( pm25_cf_1_b, MINIMUM_VALID_PM, MAXIMUM_VALID_PM );

    if ( usable ) {
      const double sum = pm25_cf_1_a + pm25_cf_1_b;

      /* Check if a and b channel measures differ by <= given ug/m3: */

      const double absolute_difference =
        pm25_cf_1_a > pm25_cf_1_b ?
          pm25_cf_1_a - pm25_cf_1_b
        : pm25_cf_1_b - pm25_cf_1_a;

      usable = absolute_difference <= maximum_channel_difference;

      if ( ! usable ) {

        /* Else check if a and b are <= given relative ratio: */

        const double absolute_ratio =
          ( absolute_difference + absolute_difference ) / sum;

        usable = absolute_ratio <= maximum_channel_ratio;
      }

      if ( usable ) {

        /* Compute pm25_corrected formula: */

        const double mean_pm25 = 0.5 * sum;
        const double adjusted_mean_pm25 = mean_pm25 + SIGMOID_PM25_MEAN_OFFSET;
        const double sigmoid = 1.0 / ( 1.0 + exp( -adjusted_mean_pm25 ) );
        const double linear_term =
          LINEAR_TERM_PM25_COEFFICIENT * mean_pm25 +
          LINEAR_TERM_HUMIDITY_COEFFICIENT * humidity +
          LINEAR_TERM_OFFSET;
        const double quadratic_term =
          QUADRATIC_TERM_PM25_SQUARED_COEFFICIENT * mean_pm25 * mean_pm25 +
          QUADRATIC_TERM_PM25_COEFFICIENT * mean_pm25 +
          QUADRATIC_TERM_OFFSET;
        pm25_corrected =
          ( 1.0 - sigmoid ) * linear_term + sigmoid * quadratic_term;

        /* Check to filter invalid negative or NaN results: */

        if ( ! ( allowNegativePM25Corrected || pm25_corrected >= 0.0 ) ) {
          pm25_corrected = MISSING_VALUE;
        }
      }
    }
  }

#undef SIGMOID_PM25_MEAN_OFFSET

#undef LINEAR_TERM_PM25_COEFFICIENT
#undef LINEAR_TERM_HUMIDITY_COEFFICIENT
#undef LINEAR_TERM_OFFSET

#undef QUADRATIC_TERM_PM25_SQUARED_COEFFICIENT
#undef QUADRATIC_TERM_PM25_COEFFICIENT
#undef QUADRATIC_TERM_OFFSET

  return pm25_corrected;
}


#else

/******************************************************************************
PURPOSE: pm25_corrected_piecewise - Compute pm25_corrected from pm25_atm_a/b
         channel measures and relative humidity with filtering by
         absolute difference and ratio in channel measures.
INPUTS:  const double pm25_atm_a    Channel A PM2.5 measure (ug/m3).
         const double pm25_atm_b    Channel B PM2.5 measure (ug/m3).
         const double humidity      Relative humidity measure (%).
         const double maximum_channel_difference
                                    Maximum acceptable absolute difference of
                                    channel measures (ug/m3).
         const double maximum_channel_ratio
                                    Maximum acceptable absolute ratio of
                                    channel measures [0.0, 1.0].
RETURNS: double pm25_corrected or returns MISSING_VALUE if invalid values are
         given or generated.
NOTES: This version is used for data with channel A & B PM25 both available.
       Formula inputs (atm) revised on 2022-01-11 per Karoline Barkjohn email.
******************************************************************************/

static double pm25_corrected_piecewise(
  const double pm25_atm_a,
  const double pm25_atm_b,
  const double humidity,
  const double maximum_channel_difference,
  const double maximum_channel_ratio ) {

  double pm25_corrected = MISSING_VALUE; /* Default result. */

  /* Formula constants: */

#define PM25_LIMIT1  30.0
#define PM25_LIMIT2  50.0
#define PM25_LIMIT3 210.0
#define PM25_LIMIT4 260.0

#define ADJUSTED_PM25_SCALE1    0.05
#define ADJUSTED_PM25_SCALE2    0.02
#define ADJUSTED_PM25_OFFSET1 (-1.5)
#define ADJUSTED_PM25_OFFSET2 (-4.2)

#define PM25_COEFFICIENT1       0.524
#define PM25_COEFFICIENT2       0.786
#define PM25_COEFFICIENT3       0.69
#define PM25_COEFFICIENT4       0.000884

#define HUMIDITY_COEFFICIENT  (-0.0862)

#define OFFSET1                 5.75
#define OFFSET2                 2.966

  /* Check if humidity is usable: */

  int usable =
    IN_RANGE( humidity, MINIMUM_VALID_HUMIDITY, MAXIMUM_VALID_HUMIDITY );

  assert( IN_RANGE( maximum_channel_difference, 0.0, 100.0 ) );
  assert( IN_RANGE( maximum_channel_ratio, 0.0, 1.0 ) );

  if ( usable ) {

    /* Check if a and b channel measures are usable: */

    usable =
      IN_RANGE( pm25_atm_a, MINIMUM_VALID_PM, MAXIMUM_VALID_PM ) &&
      IN_RANGE( pm25_atm_b, MINIMUM_VALID_PM, MAXIMUM_VALID_PM );

    if ( usable ) {
      const double sum = pm25_atm_a + pm25_atm_b;

      /* Check if a and b channel measures differ by <= given ug/m3: */

      const double absolute_difference =
        pm25_atm_a > pm25_atm_b ?
          pm25_atm_a - pm25_atm_b
        : pm25_atm_b - pm25_atm_a;

      usable = absolute_difference <= maximum_channel_difference;

      if ( ! usable ) {

        /* Else check if a and b are <= given relative ratio: */

        const double absolute_ratio =
          ( absolute_difference + absolute_difference ) / sum;

        usable = absolute_ratio <= maximum_channel_ratio;
      }

      if ( usable ) {

        /* Compute pm25_corrected formula: */

        const double mean_pm25 = 0.5 * sum;

        if ( mean_pm25 < PM25_LIMIT1 ) {
          pm25_corrected =
            PM25_COEFFICIENT1 * mean_pm25 +
            HUMIDITY_COEFFICIENT * humidity +
            OFFSET1;
        } else if ( mean_pm25 < PM25_LIMIT2 ) {
          const double adjusted_mean_pm25 =
            mean_pm25 * ADJUSTED_PM25_SCALE1 + ADJUSTED_PM25_OFFSET1;
          const double one_minus_adjusted_mean_pm25 = 1.0 - adjusted_mean_pm25;
          pm25_corrected =
            ( PM25_COEFFICIENT2 * adjusted_mean_pm25 +
              PM25_COEFFICIENT1 * one_minus_adjusted_mean_pm25 ) * mean_pm25 +
            HUMIDITY_COEFFICIENT * humidity +
            OFFSET1;
        } else if ( mean_pm25 < PM25_LIMIT3 ) {
          pm25_corrected =
            PM25_COEFFICIENT2 * mean_pm25 +
            HUMIDITY_COEFFICIENT * humidity +
            OFFSET1;
        } else if ( mean_pm25 < PM25_LIMIT4 ) {
          const double adjusted_mean_pm25 =
            mean_pm25 * ADJUSTED_PM25_SCALE2 + ADJUSTED_PM25_OFFSET2;
          const double one_minus_adjusted_mean_pm25 = 1.0 - adjusted_mean_pm25;
          const double term1 =
            ( PM25_COEFFICIENT3 * adjusted_mean_pm25 +
              PM25_COEFFICIENT2 * one_minus_adjusted_mean_pm25 ) * mean_pm25;
          const double term2 =
            HUMIDITY_COEFFICIENT * humidity * one_minus_adjusted_mean_pm25;
          const double term3 = OFFSET2 * adjusted_mean_pm25;
          const double term4 = OFFSET1 * one_minus_adjusted_mean_pm25;
          const double term5 =
            PM25_COEFFICIENT4 * mean_pm25 * mean_pm25 * adjusted_mean_pm25;
          pm25_corrected = term1 + term2 + term3 + term4 + term5;
        } else {
          pm25_corrected =
            PM25_COEFFICIENT3 * mean_pm25 +
            PM25_COEFFICIENT4 * mean_pm25 * mean_pm25 +
            OFFSET2;
        }

        /* Check to filter invalid negative or NaN results: */

        if ( ! ( allowNegativePM25Corrected || pm25_corrected >= 0.0 ) ) {
          pm25_corrected = MISSING_VALUE;
        }
      }
    }
  }

#undef PM25_LIMIT1
#undef PM25_LIMIT2
#undef PM25_LIMIT3
#undef PM25_LIMIT4

#undef ADJUSTED_PM25_SCALE1
#undef ADJUSTED_PM25_SCALE2
#undef ADJUSTED_PM25_OFFSET1
#undef ADJUSTED_PM25_OFFSET1

#undef PM25_COEFFICIENT1
#undef PM25_COEFFICIENT2
#undef PM25_COEFFICIENT3
#undef PM25_COEFFICIENT4

#undef HUMIDITY_COEFFICIENT

#undef OFFSET1
#undef OFFSET2

  return pm25_corrected;
}

#endif


/*---------------------------------------------------------------------------*/

/* Output vars: timestamp, longitude, latitude, elevation, id, count, pm25: */

enum { VARIABLES = 7 };

enum { TEMP_FILE_1 = VARIABLES, TEMP_FILE_2, TEMP_AGGREGATED_FILE, TEMP_FILES};

static const char* const tempFileNames[ TEMP_FILES ] = {
  "timestamp", "longitude", "latitude", "elevation", "id", "count", "data",
  "temp1", "temp2", "aggregated"
};

enum { FORMAT_XDR, FORMAT_ASCII };
#define FORMAT_STRING "xdr ascii"

enum {
  AGGREGATE_NONE, AGGREGATE_ALL, AGGREGATE_HOURLY, AGGREGATE_DAILY,
  AGGREGATE_MONTHLY
};
#define AGGREGATE_STRING "none all hourly daily monthly"

enum {
  ID_INDEX, TIMESTAMP_INDEX, LONGITUDE_INDEX, LATITUDE_INDEX, ELEVATION_INDEX,
  INSIDE_INDEX, CHANNEL_STATE_INDEX, CHANNEL_FLAGS_INDEX,
  PM25_CF_1_A_INDEX, PM25_CF_1_B_INDEX, /* pm25_cf_1_a/b used for sigmoid */
  PM25_ATM_A_INDEX, PM25_ATM_B_INDEX,   /* pm25_atm_a/b  used for piecewise */
  HUMIDITY_INDEX, DESCRIPTION_INDEX, VARIABLE_INDEX, COLUMN_INDICES
};


/*================================== TYPES ==================================*/


/*
 * List of known variable names: variable,units,field
 * variable is the name of the variable in purpleairserver (or implicit)
 * units is the output units of the variable in purpleairserver
 * field is the name of the field/column in the .json input file.
 */

typedef struct {
  const char* name;          /* Name of output variable. "pm25_corrected". */
  const char* units;         /* Units of output variable. "ug/m3". */
  const char* field;         /* Column name in input file. "pm2.5". */
  const double validMinimum; /* Minimum valid value. */
  const double validMaximum; /* Maximum valid value. */
  int         index;         /* 0-based column number of variable. */
} ColumnInfo;

static ColumnInfo columnInfo[] = {
  { "id",              "-",               "sensor_index",  1, LLONG_MAX, -1 },
  { "description",     "-",               "name",          0, 0, -1 },
  { "timestamp",       "s",               "last_seen",     0, INT_MAX, -1 },
  { "longitude",       "deg",             "longitude",     -180.0, 180.0, -1 },
  { "latitude",        "deg",             "latitude",      -90.0, 90.0, -1 },
  { "elevation",       "m",               "altitude",
    MINIMUM_VALID_ELEVATION_METERS, MAXIMUM_VALID_ELEVATION_METERS, -1 },
  { "channel_state",   "-",               "channel_state",    0, 3, -1 },
  { "channel_flag",    "-",               "channel_flags",    0, 3, -1 },
  { "inside",          "-",               "location_type",    0, 1, -1 },

  { "pm1",             "ug/m3",           "pm1.0",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_a",           "ug/m3",           "pm1.0_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_b",           "ug/m3",           "pm1.0_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_atm",         "ug/m3",           "pm1.0_atm",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_atm_a",       "ug/m3",           "pm1.0_atm_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_atm_b",       "ug/m3",           "pm1.0_atm_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_cf_1",        "ug/m3",           "pm1.0_cf_1",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_cf_1_a",      "ug/m3",           "pm1.0_cf_1_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm1_cf_1_b",      "ug/m3",           "pm1.0_cf_1_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },

  { "pm25",            "ug/m3",           "pm2.5",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_a",          "ug/m3",           "pm2.5_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_b",          "ug/m3",           "pm2.5_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_atm",        "ug/m3",           "pm2.5_atm",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_atm_a",      "ug/m3",           "pm2.5_atm_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_atm_b",      "ug/m3",           "pm2.5_atm_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_cf_1",       "ug/m3",           "pm2.5_cf_1",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_cf_1_a",     "ug/m3",           "pm2.5_cf_1_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_cf_1_b",     "ug/m3",           "pm2.5_cf_1_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },

  { "pm25_10minute",   "ug/m3",           "pm2.5_10minute",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_10minute_a", "ug/m3",           "pm2.5_10minute_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_10minute_b", "ug/m3",           "pm2.5_10minute_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_60minute",   "ug/m3",           "pm2.5_60minute",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_60minute_a", "ug/m3",           "pm2.5_60minute_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm25_60minute_b", "ug/m3",           "pm2.5_60minute_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },

  { "pm10",            "ug/m3",           "pm10.0",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_a",          "ug/m3",           "pm10.0_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_b",          "ug/m3",           "pm10.0_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_atm",        "ug/m3",           "pm10.0_atm",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_atm_a",      "ug/m3",           "pm10.0_atm_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_atm_b",      "ug/m3",           "pm10.0_atm_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_cf_1",       "ug/m3",           "pm10.0_cf_1",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_cf_1_a",     "ug/m3",           "pm10.0_cf_1_a",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },
  { "pm10_cf_1_b",     "ug/m3",           "pm10.0_cf_1_b",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 },

  { "0_3_um_count",    "particles/100ml", "0.3_um_count",   0.0, 1e6, -1 },
  { "0_3_um_count_a",  "particles/100ml", "0.3_um_count_a", 0.0, 1e6, -1 },
  { "0_3_um_count_b",  "particles/100ml", "0.3_um_count_b", 0.0, 1e6, -1 },

  { "0_5_um_count",    "particles/100ml", "0.5_um_count",   0.0, 1e6, -1 },
  { "0_5_um_count_a",  "particles/100ml", "0.5_um_count_a", 0.0, 1e6, -1 },
  { "0_5_um_count_b",  "particles/100ml", "0.5_um_count_b", 0.0, 1e6, -1 },

  { "1_um_count",      "particles/100ml", "1.0_um_count",   0.0, 1e6, -1 },
  { "1_um_count_a",    "particles/100ml", "1.0_um_count_a", 0.0, 1e6, -1 },
  { "1_um_count_b",    "particles/100ml", "1.0_um_count_b", 0.0, 1e6, -1 },

  { "2_5_um_count",    "particles/100ml", "2.5_um_count",   0.0, 1e6, -1 },
  { "2_5_um_count_a",  "particles/100ml", "2.5_um_count_a", 0.0, 1e6, -1 },
  { "2_5_um_count_b",  "particles/100ml", "2.5_um_count_b", 0.0, 1e6, -1 },

  { "5_um_count",      "particles/100ml", "5.0_um_count",   0.0, 1e6, -1 },
  { "5_um_count_a",    "particles/100ml", "5.0_um_count_a", 0.0, 1e6, -1 },
  { "5_um_count_b",    "particles/100ml", "5.0_um_count_b", 0.0, 1e6, -1 },

  { "10_um_count",     "particles/100ml", "10.0_um_count",   0.0, 1e6, -1 },
  { "10_um_count_a",   "particles/100ml", "10.0_um_count_a", 0.0, 1e6, -1 },
  { "10_um_count_b",   "particles/100ml", "10.0_um_count_b", 0.0, 1e6, -1 },

  { "humidity",        "%",               "humidity",
    MINIMUM_VALID_HUMIDITY, MAXIMUM_VALID_HUMIDITY, -1 },
  { "temperature",     "C",               "temperature", -100.0, 100.0, -1 },
  { "pressure",        "hPa",             "pressure", 500.0, 1500.0, -1 },
  { "voc",             "IAQ",             "voc", 0, 1e6, -1 },
  { "ozone1",          "ppb",             "ozone1", 0.0, 1e6, -1 },

  { "pm25_corrected",  "ug/m3",           "pm2.5",
    MINIMUM_VALID_PM, MAXIMUM_VALID_PM, -1 }
};

enum { NAME_LENGTH = 256 };
typedef char FileName[ NAME_LENGTH ];

enum { NOTE_LENGTH = 79 };
typedef char Note[ NOTE_LENGTH + 1 ];

/* User-supplied command-line arguments: */

typedef struct {
  const char* listFile;       /* File containing list of input files to read.*/
  const char* tmpdir;         /* Name of directory to write temp files. */
  const char* description;    /* User-supplied description. */
  const char* variable;       /* Name of variable to read. */
  Bounds      bounds; /* Subset bounds[LONGITUDE,LATITUDE][MINIMUM,MAXIMUM]. */
  long long   yyyymmddhhmmss[ 2 ]; /* Beginning/ending timestamp of subset. */
  double      maximumDifference; /* Maximum acceptable channel A,B difference*/
  double      maximumRatio;      /* Maximum acceptable channel A,B ratio. */
  double      defaultHumidity; /* Value used if humidity is missing/invalid. */
  double      minimumAggregationCountPercentage; /* Default is 75%. */
  int         format;         /* FORMAT_XDR, FORMAT_ASCII. */
  int         sensor;          /* 0 = all. > 0 for specific sensor id. */
  int         out_in_flag;     /* 0 = outside, 1 = inside, 2=either (default)*/
  int         aggregate;       /* AGGREGATE_NONE/ALL/HOURLY/DAILY/MONTHLY. */
} Arguments;

/* Data type: */

typedef struct {
  Arguments   arguments;     /* User-supplied (command-line) arguments. */
  const char* units;         /* Units of output variable. E.g., "ug/m3". */
  double validMinimum;       /* Minimum valid value for variable. */
  double validMaximum;       /* Maximum valid value for variable. */
  FileName    tempFileNames[ TEMP_FILES ]; /* Pathed names of temp files. */
  FILE*       tempFiles[ TEMP_FILES ]; /* Temp files of output subset data. */
  size_t      bufferSize;    /* Byte size of input/output buffer. */
  size_t      inputLength;   /* Length of input buffer. */
  char*       inputBuffer;   /* Holds current input file content. */
  char*       outputBuffer;  /* Holds an output subset data content. */
  long long   seconds1;      /* Seconds from 1970 to yyyymmddhhmmss[ 0 ]. */
  long long   seconds2;      /* Seconds from 1970 to yyyymmddhhmmss[ 1 ]. */
  int         columnIndices[ COLUMN_INDICES ]; /*0-based index of needed cols*/
  double      columnValues[  COLUMN_INDICES ]; /* Value of required columns. */
  Note        note;          /* Column description. */
  size_t      points;        /* Number of valid data points in subset. */
  int         ok;            /* Did last command succeed? */
} Data;

/*========================== FORWARD DECLARATIONS ===========================*/

static void removeTempFiles( Data* const data );

static void createTempFiles( Data* const data );

static void closeTempFiles( Data* const data );

static void openVariableTempFiles( Data* const data );

static void allocateBuffers( Data* const data );

static void deallocateBuffers( Data* const data );

static void reallocateOutputBuffer( Data* const data );

static void printUsage( const char* const name );

static int parseArguments( int argc, char* argv[], Arguments* const arguments);

static int checkAggregate( const int aggregate,
                           const char* const fileName,
                           int* yyyymmddhh );

static long long computeMinimumCount( const Data* const data,
                                      const long long timestamp );

static long long dataFileTimestamp( const char* const fileName );

static void parseColumnIndices( Data* data );

static void readData( Data* const data );

static void extractSubset( Data* data );

static int parseColumnValues( char* line,
                              const long long seconds1,
                              const long long seconds2,
                              const double west,
                              const double east,
                              const double south,
                              const double north,
                              const int out_in_flag,
                              const int sensor,
                              const int isChannelDependent,
                              const int isChannelA,
                              const int isChannelB,
                              const int isPM25Corrected,
                              const int isTemperature,
                              const double maximumDifference,
                              const double maximumRatio,
                              const double defaultHumidity,
                              const char* const units,
                              const double validMinimum,
                              const double validMaximum,
                              const int columnIndices[],
                              double columnValues[],
                              Note note );

static void aggregateData( Data* data );

static void reformatData( Data* data );

static void writeTempVariableFiles( Data* data );

static void sortTempData( Data* data );

static void streamData( Data* const data );

static void streamHeader( const Data* const data );


/*================================ FUNCTIONS ================================*/


/******************************************************************************
PURPOSE: main - Extract a subset of data from a list of PurpleAir files and write
         it to stdout in XDR format.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Data data;
 memset( &data, 0, sizeof data );
  data.ok = parseArguments( argc, argv, &data.arguments );

  if ( ! data.ok ) {
    printUsage( argv[ 0 ] );
  } else {
    createTempFiles( &data );

    if ( data.ok ) {
      allocateBuffers( &data );

      if ( data.ok ) {
        readData( &data ); /* Read input data and write temp files. */

        if ( data.ok && data.points > 0 ) {
          reformatData( &data );

          if ( data.ok ) {
            streamData( &data ); /* Write header & temp files to stdout. */
          }
        }
      }
    }

    if ( ! data.ok ) {
      fprintf( stderr, "\n%s: No points were in the subset.\n", argv[ 0 ] );
    }
  }

#ifndef DEBUGGING
  removeTempFiles( &data );
#endif
  deallocateBuffers( &data );
  data.ok = data.ok && data.points > 0;
  return ! data.ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: removeTempFiles - Close and remove temp files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void removeTempFiles( Data* const data ) {
  const size_t count =
    data ? sizeof data->tempFiles / sizeof data->tempFiles[0] : 0;
  size_t index = 0;
  assert( data );

  for ( index = 0; index < count; ++index ) {

    if ( data->tempFiles[ index ] ) {
      fclose( data->tempFiles[ index ] );
      data->tempFiles[ index ] = 0;
    }

    if ( data->tempFileNames[ index ] && data->tempFileNames[ index ][ 0 ] ) {
      unlink( data->tempFileNames[ index ] );
      memset( data->tempFileNames[ index ],
              0, sizeof data->tempFileNames[ index ] );
    }
  }
}



/******************************************************************************
PURPOSE: openVariableTempFiles - Open temp variable output files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->ok, data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void openVariableTempFiles( Data* const data ) {
  size_t index = 0;
  assert( data );

  for ( index = 0; data->ok && index < VARIABLES; ++index ) {
    assert( data->tempFileNames[ index ] );
    assert( data->tempFileNames[ index ][ 0 ] );
    data->tempFiles[ index ] = fopen( data->tempFileNames[ index ], "wb" );

    if ( ! data->tempFiles[ index ] ) {
      fprintf( stderr, "\nCan't create temporary output file '%s'.\n",
               data->tempFileNames[ index ] );
      data->ok = 0;
    }
  }
}



/******************************************************************************
PURPOSE: createTempFiles - Create temp output files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->ok, data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void createTempFiles( Data* const data ) {
  assert( data );

  {
    const int pid = getpid();
    size_t index = 0;
    data->ok = 1;

    for ( index = 0; data->ok && index < TEMP_FILES; ++index ) {
      memset( data->tempFileNames[ index ], 0, sizeof (FileName) );
      assert( data->arguments.tmpdir ); assert( data->arguments.tmpdir[ 0 ] );
      assert( tempFileNames[ index ] ); assert( tempFileNames[ index ][ 0 ] );
      snprintf( data->tempFileNames[ index ],
                sizeof (FileName) / sizeof (char) - 1,
                "%s/%s_%s.%d",
                data->arguments.tmpdir, TEMP_FILE_NAME,
                tempFileNames[ index ], pid );
      data->tempFiles[ index ] = fopen( data->tempFileNames[ index ], "wb" );

      if ( ! data->tempFiles[ index ] ) {
        fprintf( stderr, "\nCan't create temporary output file '%s'.\n",
                 data->tempFileNames[ index ] );
        data->ok = 0;
      }
    }
  }
}



/******************************************************************************
PURPOSE: closeTempFiles - Close temp files.
INPUTS:  Data* const data  data->tempFileNames[], data->tempfiles[].
OUTPUTS: Data* const data  data->tempFileNames[], data->tempfiles[].
******************************************************************************/

static void closeTempFiles( Data* const data ) {
  const size_t count =
    data ? sizeof data->tempFiles / sizeof data->tempFiles[0] : 0;
  size_t index = 0;
  assert( data );


  for ( index = 0; index < count; ++index ) {

    if ( data->tempFiles[ index ] ) {
      fclose( data->tempFiles[ index ] );
      data->tempFiles[ index ] = 0;
    }
  }
}



/******************************************************************************
PURPOSE: allocateBuffers - Allocate input/output buffers.
OUTPUTS: Data* const data  data->inputLength  Allocated length of input buffer.
                           data->inputBuffer  Allocated input buffer.
                           data->outputBuffer Allocated output buffer.
                           data->ok 1 if successful, else 0.
******************************************************************************/

static void allocateBuffers( Data* const data ) {
  assert( data );
  assert( data->inputBuffer  == 0 );
  assert( data->outputBuffer == 0 );
  data->bufferSize = 3 * 1024 * 1024 * sizeof (char);
  data->inputBuffer = malloc( data->bufferSize );
  data->outputBuffer = data->inputBuffer ? malloc( data->bufferSize ) : 0;
  data->ok = data->outputBuffer != 0;

  if ( ! data->ok ) {
      fprintf( stderr,
               "\nCan't allocate %lu bytes "
               "to complete the requested action.\n", 2 * data->bufferSize );

    if ( data->inputBuffer ) {
      free( data->inputBuffer ), data->inputBuffer = 0;
    }

    data->bufferSize = 0;
  } else {
    memset( data->inputBuffer,  0, data->bufferSize );
    memset( data->outputBuffer, 0, data->bufferSize );
  }

  data->inputLength = 0;
  assert( data->inputLength == 0 );
}



/******************************************************************************
PURPOSE: deallocateBuffers - Deallocate input/output buffers.
OUTPUTS: Data* const data  data->inputLength  0.
                           data->inputBuffer  0.
                           data->outputBuffer 0.
******************************************************************************/

static void deallocateBuffers( Data* const data ) {
  assert( data );

  if ( data->inputBuffer ) {
    free( data->inputBuffer ), data->inputBuffer = 0;
  }

  if ( data->outputBuffer ) {
    free( data->outputBuffer ), data->outputBuffer = 0;
  }

  data->bufferSize = 0;
  data->inputLength = 0;

  assert( data->inputBuffer == 0 );
  assert( data->outputBuffer == 0 );
  assert( data->bufferSize  == 0 );
  assert( data->inputLength == 0 );
}



/******************************************************************************
PURPOSE: reallocateOutputBuffer - Reallocate output buffer.
INPUTS:  Data* const data  data->inputLength  Allocated length of input buffer.
                           data->outputBuffer Allocated output buffer.
OUTPUTS: Data* const data  data->bufferSize   Size of output buffer.
                           data->outputBuffer Allocated output buffer.
                           data->ok 1 if successful, else 0.
******************************************************************************/

static void reallocateOutputBuffer( Data* const data ) {
  assert( data );
  assert( data->inputBuffer );
  assert( data->inputLength * sizeof (char) > data->bufferSize );

  if ( data->outputBuffer ) {
    free( data->outputBuffer ), data->outputBuffer = 0;
  }

  data->bufferSize = data->inputLength * sizeof (char);
  data->outputBuffer = malloc( data->bufferSize );
  data->ok = data->outputBuffer != 0;

  if ( ! data->ok ) {
    fprintf( stderr,
             "\nCan't allocate %lu bytes "
             "to complete the requested action.\n", data->bufferSize );
    data->bufferSize = 0;
  } else {
    memset( data->outputBuffer, 0, data->bufferSize );
  }

  assert( (data->ok == 0 && data->outputBuffer == 0 && data->bufferSize == 0) ||
          (data->ok == 1 && data->outputBuffer != 0 &&
            data->bufferSize == data->inputLength * sizeof (char) ) );
}



/******************************************************************************
PURPOSE: printUsage - Print program usage instructions.
INPUTS:  const char* name  Name of program.
******************************************************************************/

static void printUsage( const char* name ) {
  assert( name ); assert( *name );
  fprintf( stderr,
           "\n%s - Extract a subset of data from a time-sorted list of\n"
          "PurpleAir files and write it to stdout in XDR binary format.\n",
           name );
  fprintf( stderr, "Data is subsetted by " );
  fprintf( stderr, "date-time range, lon-lat rectangle and variable.\n");
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "  -files <listfile> \\\n" );
  fprintf( stderr, "  -tmpdir <temp_directory> \\\n" );
  fprintf( stderr, "  -desc \"description text\" \\\n" );
  fprintf( stderr, "  -timerange <yyyymmddhhmmss> <yyyymmddhhmmss> \\\n" );
  fprintf( stderr, "  -variable <name> \\\n" );
  fprintf( stderr, " [-bounds <minimum_longitude> <minimum_latitude>]" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> \\\n" );
  fprintf( stderr, "  -format ascii | xdr \\\n" );
  fprintf( stderr, " [-sensor sensor_id] (subset to specific sensor id) \\\n" );
  fprintf( stderr, " [-out_in_flag 0|1] (0 = outside, 1 = inside, default = either) \\\n" );
  fprintf( stderr, " [-maximum_difference 0-100 ug/m3] (default = 5) \\\n" );
  fprintf( stderr, " [-maximum_ratio 0-1 ] (default = 0.7) \\\n" );
  fprintf( stderr, " [-default_humidity 0-100] " );
  fprintf( stderr, "(value used for missing/invalid humidity, default = none) \\\n" );
  fprintf( stderr, " [-aggregate hourly | daily | monthly | all] (default = none) \\\n" );
  fprintf( stderr, " [-minimum_aggregation_count_percentage 0-100 (default = 75)\n\n" );
  fprintf( stderr, "Note:\ntimes are in UTC (GMT)\n" );
  fprintf( stderr, "-tmpdir specifies a directory were temp files are " );
  fprintf ( stderr, "written.\nIt should have enough disk space (1TB).\n" );
  fprintf( stderr, "-minimum_aggregation_count_percentage specifies the \n" );
  fprintf( stderr, "minimum number of values allowed for aggregation\n");
  fprintf( stderr, "expressed as a percentage of the maximum number of\n");
  fprintf( stderr, "values for the aggregation time period.\n");
  fprintf( stderr, "E.g., if -agggregate hourly and a sensor can report at\n");
  fprintf( stderr, "most 30 times per hour and \n" );
  fprintf( stderr, "-minimum_aggregation_count_percentage is 75\n");
  fprintf( stderr, "Then 75%% of 30 = 22.5 so omit sensor if it reports less\n");
  fprintf( stderr, "Than 23 values for that hour.\n");
  fprintf( stderr, "\nExample:\n\n" );
  fprintf( stderr, "%s \\\n", name );
  fprintf( stderr, "-files file_list \\\n");
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc "
                   "\"https://api.purpleair.com,PurpleAirSubset\" \\\n");
  fprintf( stderr, "-timerange 20201202000000 20201202235959 \\\n" );
  fprintf( stderr, "-variable pm25_corrected \\\n" );
  fprintf( stderr, "-bounds -124 49 -123 50 \\\n" );
  fprintf( stderr, "-format xdr \\\n" );
  fprintf( stderr, "-out_in_flag 0 \\\n" );
  fprintf( stderr, "-aggregate hourly \\\n" );
  fprintf( stderr, "> subset.xdr\n\n" );
  fprintf( stderr, "Hourly corrected PM2.5 over BC on December 2, 2020.\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary arrays:\n\n" );
  fprintf( stderr, "Point 1.0\n" );
  fprintf( stderr, "https://api.purpleair.com,PurpleAirSubset\n" );
  fprintf( stderr, "2020-12-02T00:00:00-0000 2020-12-02T23:59:59-0000\n" );
  fprintf( stderr, "# Dimensions: variables points\n" );
  fprintf( stderr, "6 20\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "timestamp longitude latitude elevation id pm25_corrected\n");
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "yyyymmddhhmmss deg deg m - ug/m3\n" );
  fprintf( stderr, "# char notes[points][80] and\n" );
  fprintf( stderr, "# IEEE-754 64-bit reals data[variables][points]:\n" );
  fprintf( stderr, "<big-endian binary format array>\n" );
  fprintf( stderr, "\n\n\n");
}



/******************************************************************************
PURPOSE: parseArguments - Parse command-line arguments.
INPUTS:  Integer argc          Number of command-line arguments.
         char* argv[]          Command-line argument strings.
OUTPUTS: Arguments* arguments  Initialized arguments.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int parseArguments(int argc, char* argv[], Arguments* const arguments) {
  int result = argc >= 14;
  static const int zero_one[ 2 ] = { 0, 1 };
  static const int one_int_max[ 2 ] = { 1, INT_MAX };
  static const double zero_100[ 2 ] = { 0.0, 100.0 };
  static const double zero_1[ 2 ] = { 0.0, 1.0 };
  static Option options[] = {
    { "-files",              1, FILE_TYPE,           1, 0, 0, 0, 0 },
    { "-tmpdir",             0, DIRECTORY_TYPE,      1, 0, 0, 0, 0 },
    { "-desc",               0, STRING_TYPE,         1, 0, 0, 0, 0 },
    { "-variable",           1, ENUM_TYPE,           1, 0, 0, 0, 0 },
    { "-timerange",          1, YYYYMMDDHHMMSS_TYPE, 2, 0, 0, 0, 0 },
    { "-format",             0, ENUM_TYPE,           1, 0, FORMAT_STRING, 0, 0},
    { "-bounds",             0, BOUNDS_TYPE,         4, 0, 0, 0, 0 },
    { "-sensor",             0, INT_TYPE,            1, one_int_max, 0, 0, 0 },
    { "-out_in_flag",        0, INT_TYPE,            1, zero_one, 0, 0, 0 },
    { "-default_humidity",   0, REAL64_TYPE,         1, zero_100, 0, 0, 0 },
    { "-maximum_difference", 0, REAL64_TYPE,         1, zero_100, 0, 0,0 },
    { "-maximum_ratio",      0, REAL64_TYPE,         1, zero_1, 0, 0,0 },
    { "-aggregate",          0, ENUM_TYPE,           1, 0,AGGREGATE_STRING,0,0},
    { "-minimum_aggregation_count_percentage",
                             0, REAL64_TYPE,         1, zero_100, 0, 0, 0 }
  };
  char variableNames[ sizeof columnInfo / sizeof columnInfo[ 0 ] * 32 ] = "";
  int variableIndex = -1;

  assert( argc > 0 ); assert( argv ); assert( argv[ 0 ] );
  assert( argv[ argc - 1 ] ); assert( arguments );

  if ( result ) {

    /* Initialize arguments to defaults: */

    memset( arguments, 0, sizeof (Arguments) );
    arguments->tmpdir = ".";
    arguments->description = "https://api.purpleair.com,PurpleAirSubset";
    arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
    arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
    arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
    arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
    arguments->maximumDifference = defaultMaximumDifference;
    arguments->maximumRatio      = defaultMaximumRatio;
    arguments->defaultHumidity = MISSING_VALUE;
    arguments->minimumAggregationCountPercentage =
      defaultMinimumAggregationCountPercentage;
    arguments->out_in_flag = 2; /* By default, accept both outside and inside*/

    /* Init space-delimited string of valid variable names for parse check: */

    {
      const size_t variables = sizeof columnInfo / sizeof columnInfo[ 0 ];
      size_t variable = 0;
      memset( variableNames, 0, sizeof variableNames );

      for ( variable = 0; variable < variables; ++variable ) {
        assert( strlen( variableNames ) + 32
                < sizeof variableNames / sizeof *variableNames );
        strncat( variableNames, columnInfo[ variable ].name, 30 );
        strcat( variableNames, " " );
      }
    }

    /* Finish initializing non-compile-time-constant parts of options: */

    options[ 0  ].values = &arguments->listFile;
    options[ 1  ].values = &arguments->tmpdir;
    options[ 2  ].values = &arguments->description;
    options[ 3  ].valids = variableNames; /* Space-delimited var names. */
    options[ 3  ].values = &variableIndex; /* Temp int to get index of name. */
    options[ 4  ].values = &arguments->yyyymmddhhmmss[ 0 ];
    options[ 5  ].values = &arguments->format;
    options[ 6  ].values = &arguments->bounds[ 0 ][ 0 ];
    options[ 7  ].values = &arguments->sensor;
    options[ 8  ].values = &arguments->out_in_flag;
    options[ 9  ].values = &arguments->defaultHumidity;
    options[ 10 ].values = &arguments->maximumDifference;
    options[ 11 ].values = &arguments->maximumRatio;
    options[ 12 ].values = &arguments->aggregate;
    options[ 13 ].values = &arguments->minimumAggregationCountPercentage;

    result =
      parseOptions( argc, argv, sizeof options / sizeof *options, options );

    if ( result ) { /* Initialize argument->variable to static string name: */
      arguments->variable = columnInfo[ variableIndex ].name;
    }

  } else {
    fprintf( stderr, "\n%s: Invalid/insufficient command-line arguments.\n",
             argv[ 0 ] );
  }

  assert( result == 0 ||
         ( arguments->listFile && arguments->listFile[ 0 ] &&
           arguments->tmpdir && arguments->tmpdir[ 0 ] &&
           arguments->description && arguments->description[ 0 ] &&
           arguments->variable && arguments->variable[ 0 ] &&
           ( arguments->format == FORMAT_XDR ||
             arguments->format == FORMAT_ASCII ) &&
           isValidYYYYMMDDHHMMSS( arguments->yyyymmddhhmmss[ 0 ] ) &&
           isValidYYYYMMDDHHMMSS( arguments->yyyymmddhhmmss[ 1 ] ) &&
           arguments->yyyymmddhhmmss[ 0 ] <= arguments->yyyymmddhhmmss[ 1 ] &&
           isValidBounds( (const double (*)[2]) arguments->bounds ) &&
           IN_RANGE( arguments->maximumDifference, 0.0, 100.0 ) &&
           IN_RANGE( arguments->maximumRatio, 0.0, 1.0 ) &&
           IN_RANGE( arguments->minimumAggregationCountPercentage,0.0,100.0) &&
           IN_RANGE( arguments->out_in_flag, 0, 2 ) &&
           ( arguments->defaultHumidity == MISSING_VALUE ||
             IN_RANGE( arguments->defaultHumidity, 0.0, 100.0 ) ) &&
           arguments->sensor >= 0 ) );

  return result;
}



/******************************************************************************
PURPOSE: readData - Read data from each listed data file and
         write the subset of data to the temporary files.
INPUTS:  Data* data  Data to read.
OUTPUTS: Data* data  data->points, ok, tempFiles[] = 0 closed.
******************************************************************************/

static void readData( Data* const data ) {
  const Arguments* const arguments = &( data->arguments );
  const int aggregate = arguments->aggregate;
  int yyyymmddhh = 0;
  size_t length = 0;
  char* listFileContent = 0;
  data->ok = readFile( arguments->listFile, &length, &listFileContent );

  if ( data->ok ) {
    int wroteSomeData = 0;
    char* inputFileName = 0;
    char* end = 0;
    data->seconds1 = secondsSince1970( arguments->yyyymmddhhmmss[ 0 ] );
    data->seconds2 = secondsSince1970( arguments->yyyymmddhhmmss[ 1 ] );

    /* Get each line of list file. It is the data file to read: */

    for ( inputFileName = strtok_r( listFileContent, "\n", &end );
          inputFileName;
          inputFileName = strtok_r( 0, "\n", &end ) ) {

      const int aggregateNow = aggregate != AGGREGATE_NONE &&
        checkAggregate( aggregate, inputFileName, &yyyymmddhh );

      data->ok = 1; /* Reset to ok to carry on after encountering a bad file.*/

      if ( aggregateNow ) {

        /*
         * Before aggregating, we must sort temp file 1 into temp file 2
         * because aggregating reads the sorted temp 2 file.
         * After sorting, temp file 1 is emptied in preparation for subsequent
         * appends by extractSubset().
         */

        sortTempData( data );

        if ( data->ok ) {

          /* Aggregate sorted temp file 2 appending onto temp aggregate file:*/

          aggregateData( data );
        }
      }

      if ( data->ok ) {
        data->ok =
          readFile( inputFileName, &data->inputLength, &data->inputBuffer );

        if ( data->ok ) {

          if ( data->inputLength * sizeof (char) > data->bufferSize ) {
            reallocateOutputBuffer( data );
          }

          if ( data->ok ) {
            parseColumnIndices( data );

            if ( data->ok ) {
              extractSubset( data ); /* Appends to temp file 1. */

              if ( data->ok ) {
                wroteSomeData = 1;
              }
            }

            if ( ! data->ok ) {
              fprintf( stderr, "\nOmitting invalid file %s\n", inputFileName );
            }
          }
        } else {
          data->ok = 1; /* Ignore bad input files. */
        }
      }
    } /* End loop on listFile. */

    free( listFileContent ), listFileContent = 0;
    data->ok = wroteSomeData;

    if ( wroteSomeData ) {

      /* Sort temp file 1 into temp file 2 and empty temp file 1: */

      sortTempData( data );

      if ( data->ok && aggregate != AGGREGATE_NONE ) {

        /* Aggregate sorted temp file 2 appending onto temp aggregate file: */

        aggregateData( data );
        DEBUG( fprintf( stderr, "after final aggregateData() data->ok = %d\n",
                        data->ok ); )
      }

      closeTempFiles( data );
    }
  }

  DEBUG( fprintf( stderr,
                  "End of readData() data->points = %lu, data->ok = %d\n",
                  data->points, data->ok ); )
}



/******************************************************************************
PURPOSE: checkAggregate - Check if we need to aggregate now.
INPUTS:  const int aggregate         AGGREGATE_HOURLY, AGGREGATE_DAILY...
         const char* const fileName  Pathed name of PurpleAir file.
         int* yyyymmddhh             Timestamp to compare or 0 to initialize.
OUTPUT:  int* yyyymmddhh             Updated to file timestamp if result == 1.
RETURNS: int 1 if time to aggregate now, else 0.
NOTES:   Data file names are of the form:
           .../YYYYMMDD/purpleair_hhmmss.json
         For example,
           /data/PurpleAir/data/2020/20201202/purpleair_181802.json
******************************************************************************/

static int checkAggregate( const int aggregate,
                           const char* const fileName,
                           int* yyyymmddhh ) {
  int result = 0;
  assert( aggregate == AGGREGATE_NONE   ||
          aggregate == AGGREGATE_ALL    ||
          aggregate == AGGREGATE_HOURLY ||
          aggregate == AGGREGATE_DAILY  ||
          aggregate == AGGREGATE_MONTHLY );
  assert( fileName ); assert( *fileName );
  assert( yyyymmddhh );
  assert( *yyyymmddhh == 0 || isValidYYYYMMDDHH( *yyyymmddhh ) );

  /* Set/check yyyymmddhh for aggregation: */

  if ( aggregate != AGGREGATE_NONE &&
       aggregate != AGGREGATE_ALL ) {
    const long long fileYYYYMMDDHHMMSS = dataFileTimestamp( fileName );

    if ( fileYYYYMMDDHHMMSS ) {
      const int fileYYYYMMDDHH = fileYYYYMMDDHHMMSS / 10000;
      assert( isValidYYYYMMDDHH( fileYYYYMMDDHH ) );

      if ( *yyyymmddhh == 0 ) { /* Initialize: */
        *yyyymmddhh = fileYYYYMMDDHH;
      } else if ( aggregate == AGGREGATE_HOURLY ) {

        if ( fileYYYYMMDDHH > *yyyymmddhh ) {
          *yyyymmddhh = fileYYYYMMDDHH;
          result = 1;
        }
      } else if ( aggregate == AGGREGATE_DAILY ) {
        const int fileYYYYMMDD = fileYYYYMMDDHH / 100;
        const int yyyymmdd     = *yyyymmddhh    / 100;

        if ( fileYYYYMMDD > yyyymmdd ) {
          *yyyymmddhh = fileYYYYMMDDHH;
          result = 1;
        }
      } else if ( aggregate == AGGREGATE_MONTHLY ) {
        const int fileYYYYMM = fileYYYYMMDDHH / 10000;
        const int yyyymm     = *yyyymmddhh    / 10000;

        if ( fileYYYYMM > yyyymm ) {
          *yyyymmddhh = fileYYYYMMDDHH;
          result = 1;
        }
      }

      assert( isValidYYYYMMDDHH( *yyyymmddhh ) );
      assert( result == 0 || ( result == 1 && *yyyymmddhh == fileYYYYMMDDHH ) );
    }
  }

  assert( result == 0 || result == 1 );
  return result;
}



/******************************************************************************
PURPOSE: dataFileTimestamp - Timestamp of pathed data file.
INPUTS:  const char* const fileName   Name of PurpleAir file.
RETURNS: long long yyyymmddhh of file or 0 if failed and message on stderr.
NOTES:   Data file names are of the form:
           .../YYYYMMDD/purpleair_hhmmss.json
         For example,
          /data/PurpleAir/data/2020/20201202/purpleair_181802.json
         which yields 20201202181802.
******************************************************************************/

static long long dataFileTimestamp( const char* const fileName ) {
  long long result = 0;
  const char* lastUnderscore = strrchr( fileName, '_' );

  if ( lastUnderscore ) {
    const int hhmmss = atoi( lastUnderscore + 1 );
    const char* slash = strrchr( fileName, '/' );

    if ( slash ) {

      while ( slash > fileName && result == 0 ) {
        --slash;

        if ( *slash == '/' ) {
          result = atoi( slash + 1 );
        }
      }

      if ( result == 0 ) {
        result = atoi( fileName );
      }

      if ( result ) {
        result *= 1000000;
        result += hhmmss;
      }
    }
  }

  if ( ! isValidYYYYMMDDHHMMSS( result ) ) {
    fprintf( stderr, "\nInvalid file name timestamp '%s'.\n", fileName );
    result = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: parseColumnIndices - Parse data column indices from input string.
INPUTS:  Data* data  data->arguments->variable  Name of output variable.
                     data->inputBuffer Content of input data file.
OUTPUTS: Data* data  data->columnIndices[] 0-based indices.
                     data->ok = 1 if successful, else 0.
******************************************************************************/

static void parseColumnIndices( Data* data ) {
  const char* const tag = "\"fields\" : [";
  char* beginSection = 0; /* Beginning of section to parse. */
  size_t index = 0;
  assert( data );
  assert( data->arguments.variable ); assert( data->arguments.variable[ 0 ] );
  assert( data->inputBuffer );

  data->units = 0;
  data->ok = 0;

  /* Initialize to invalid value: */

  for ( index = 0; index < COLUMN_INDICES; ++index ) {
    data->columnIndices[ index ] = -1;
  }

  beginSection = strstr( data->inputBuffer, tag );

  if ( beginSection ) {
    char* const endSection = strchr( beginSection += strlen( tag ), ']' );

    if ( endSection ) {
      const char* const variable = data->arguments.variable;
      const char* const delimiters = ", \"\n";
      char* word = 0;
      char* end = 0;
      int wordIndex = 0;
      *endSection = '\0'; /* Terminate portion of input to parse. */

      for ( word = strtok_r( beginSection, delimiters, &end );
            word;
            word = strtok_r( 0, delimiters, &end ), ++wordIndex ) {
        const size_t count = sizeof columnInfo / sizeof columnInfo[ 0 ];

        for ( index = 0; index < count; ++index ) {
          const char* const field = columnInfo[ index ].field;


          if ( ! strcmp( word, field ) ) {
            const char* const name = columnInfo[ index ].name;
            columnInfo[ index ].index = wordIndex;

            if ( ! strcmp( name, "id" ) ) {
              data->columnIndices[ ID_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "description" ) ) {
              data->columnIndices[ DESCRIPTION_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "inside" ) ) {
              data->columnIndices[ INSIDE_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "latitude" ) ) {
              data->columnIndices[ LATITUDE_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "longitude" ) ) {
              data->columnIndices[ LONGITUDE_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "elevation" ) ) {
              data->columnIndices[ ELEVATION_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "timestamp" ) ) {
              data->columnIndices[ TIMESTAMP_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "channel_state" ) ) {
              data->columnIndices[ CHANNEL_STATE_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "channel_flag" ) ) {
              data->columnIndices[ CHANNEL_FLAGS_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "humidity" ) ) {
              data->columnIndices[ HUMIDITY_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "pm25_cf_1_a" ) ) {
              data->columnIndices[ PM25_CF_1_A_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "pm25_cf_1_b" ) ) {
              data->columnIndices[ PM25_CF_1_B_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "pm25_atm_a" ) ) {
              data->columnIndices[ PM25_ATM_A_INDEX ] = wordIndex;
            } else if ( ! strcmp( name, "pm25_atm_b" ) ) {
              data->columnIndices[ PM25_ATM_B_INDEX ] = wordIndex;
            }

            if ( ! strcmp( name, variable ) ) {
              data->columnIndices[ VARIABLE_INDEX ] = wordIndex;
              data->units = columnInfo[ index ].units;
              data->validMinimum = columnInfo[ index ].validMinimum;
              data->validMaximum = columnInfo[ index ].validMaximum;
            }
          }
        }
      } /* End of parsing loop. */

      /* Restore input. Change each '\0' to ' ': */

      *endSection = ']';

      {
        char* c = 0;

        for ( c = beginSection; c < endSection; ++c ) {

          if ( *c == '\0' ) {
            *c = ' ';
          }
        }
      }

      /* Check that all required columns have a valid index: */

      data->ok = 1;

      for ( index = 0; data->ok && index < COLUMN_INDICES; ++index ) {
        data->ok = data->columnIndices[ index ] >= 0;
      }
    }
  }

  if ( ! data->ok ) {
    fprintf( stderr, "\nFailed to parse column names.\n" );
  }

  assert( data->ok == 0 ||
          ( data->units != 0 && data->units[ 0 ] &&
            data->columnIndices[ 0 ] >= 0 ) );
}



/******************************************************************************
PURPOSE: extractSubset - Parse required column data from inputBuffer and write
         to outputBuffer then to temp file 1.
INPUTS:  Data* data  data->columnIndices[]  0-based indices of required columns
                     data->inputBuffer Content of input data file.
OUTPUTS: Data* data  data->outputBuffer  Content of extract.
                     data->tempFiles[ TEMP_FILE_1 ]  Write outputBuffer to it.
                     data->points  Increased if AGGREGSTE_NONE, else unchanged.
                     data->ok = 1 if successful, else 0.
NOTES:   For performance, both the input data file and the subset output are
         written to internal string buffers and the output buffer is written
         to temp file 1 at the end of extracting subset for the current file.
******************************************************************************/

static void extractSubset( Data* data ) {
  const char* const tag = "\"data\" : [";
  char* beginSection = 0; /* Beginning of section to parse. */
  assert( data );
  assert( data->arguments.variable ); assert( data->arguments.variable[ 0 ] );
  assert( data->inputBuffer );
  assert( data->tempFiles[ TEMP_FILE_1 ] );
  data->outputBuffer[ 0 ] = '\0';

  beginSection = strstr( data->inputBuffer, tag );

  if ( beginSection ) {
    char* const endSection = strrchr( beginSection += strlen( tag ), ']' );
    /* Note: above line logic assumes data is last [] in input file! */

    if ( endSection ) {
      const Arguments* const arguments = &data->arguments;
      const char* const variable = arguments->variable;
      const size_t variableLength = strlen( variable );
      const int isPM = variable[ 0 ] == 'p' && variable[ 1 ] == 'm';
      const int isChannelDependent = isPM || strstr( variable, "_count" ) != 0;
      const int isChannelA = isChannelDependent &&
        variableLength > 2 &&
        variable[ variableLength - 2 ] == '_' &&
        variable[ variableLength - 1 ] == 'a';
      const int isChannelB = isChannelDependent &&
        variableLength > 2 &&
        variable[ variableLength - 2 ] == '_' &&
        variable[ variableLength - 1 ] == 'b';
      const int isPM25Corrected = isPM && ! strcmp( variable, "pm25_corrected");
      const int isTemperature = ! isPM && ! strcmp( variable, "temperature" );
      const double west  = arguments->bounds[ LONGITUDE ][ MINIMUM ];
      const double east  = arguments->bounds[ LONGITUDE ][ MAXIMUM ];
      const double south = arguments->bounds[ LATITUDE  ][ MINIMUM ];
      const double north = arguments->bounds[ LATITUDE  ][ MAXIMUM ];
      const double maximumDifference = arguments->maximumDifference;
      const double maximumRatio      = arguments->maximumRatio;
      const double defaultHumidity   = arguments->defaultHumidity;
      const long long seconds1 = data->seconds1;
      const long long seconds2 = data->seconds2;
      const int sensor = arguments->sensor;
      const int out_in_flag = arguments->out_in_flag;
      double* const columnValues  = data->columnValues;
      int*    const columnIndices = data->columnIndices;
      char* output = data->outputBuffer;
      const char* const outputEnd = output + data->bufferSize / sizeof (char);
      char* line = beginSection;
      char* endLine = 0;

      *endSection = '\0'; /* Terminate portion of input to parse. */
      eraseQuotedCommasAndQuotesAndBrackets( beginSection ); /* perf hotspot */
      memset( data->outputBuffer, 0, data->bufferSize );

      for ( ; nextLine( &line, &endLine ); line = endLine + 1 ) {
        int inSubset = 0;
        assert( strchr( line, '[' ) == 0 );
        assert( strchr( line, ']' ) == 0 );
        assert( strchr( line, '"' ) == 0 );
        assert( strchr( line, ',' ) != 0 );

        inSubset =
          parseColumnValues( line,
                             seconds1, seconds2, west, east, south, north,
                             out_in_flag, sensor,
                             isChannelDependent,
                             isChannelA, isChannelB, isPM25Corrected,
                             isTemperature,
                             maximumDifference,
                             maximumRatio,
                             defaultHumidity,
                             data->units,
                             data->validMinimum,
                             data->validMaximum,
                             columnIndices,
                             columnValues, data->note );

        /* Append subset data to string buffer for speed: */

        if ( inSubset ) {
          const int id            = columnValues[ ID_INDEX ];
          const long long seconds = columnValues[ TIMESTAMP_INDEX ];
          const double longitude = columnValues[ LONGITUDE_INDEX ];
          const double latitude  = columnValues[ LATITUDE_INDEX ];
          const double elevation = columnValues[ ELEVATION_INDEX ];
          const double measure   = columnValues[ VARIABLE_INDEX ];
          const char* const description = data->note;
          enum { OUTPUT_LINE_LENGTH = 255 };
          char outputLine[ OUTPUT_LINE_LENGTH + 1 ] = "";
          const size_t length =
            snprintf( outputLine, sizeof outputLine / sizeof *outputLine,
                      "%d,%lld,%f,%f,%f,0,%f,%s\n",
                      id, seconds, longitude, latitude, elevation,
                      measure, description );

          /* If line is not truncated and line fits in outputBuffer, append: */

          if ( length < sizeof outputLine / sizeof *outputLine &&
               output + length < outputEnd ) {
            strcpy( output, outputLine );
            output += length;

            if ( data->arguments.aggregate == AGGREGATE_NONE ) {
              data->points += 1; /* Update number of output points. */
            }
          }
        }
      } /* End of parsing data lines. */

      if ( data->outputBuffer[ 0 ] ) { /* Write string buffer to temp file: */
        data->ok =
          fputs( data->outputBuffer, data->tempFiles[ TEMP_FILE_1 ] ) != EOF;
        data->outputBuffer[ 0 ] = '\0';
      }
    }
  }
}



/******************************************************************************
PURPOSE: parseColumnValues - Parse column values from a line.
INPUTS:  const char* line                Data line to parse.
         const long long seconds1        Beginning seconds of time subset.
         const long long seconds2        Ending    seconds of time subset.
         const double west               West bound of subset.
         const double east               East bound of subset.
         const double south              South bound of subset.
         const double north              North bound of subset.
         const int out_in_flag           0 = outside, 1 = inside, 2 = either.
         const int sensor                0 = all, > 0 = selected.
         const int isChannelDependent    Is selected variable a channel-
                                         dependent measure?
         const int isChannelA            Is variable channel A only?
         const int isChannelB            Is variable channel B only?
         const int isPM25Corrected       Is variable PM2.5 corrected?
         const int isTemperature         Is variable temperature?
         const double maximumDifference  Acceptable PM-channel difference.
         const double maximumRatio       Acceptable PM-channel ratio.
         const double defaultHumidity    Substitute value for humidity if it is
                                         missing or invalid.
         const char* const units         Units of measure.
         const double validMinimum       Minimum valid value of measure.
         const double validMaximum       Maximum valid value of measure.
         const int columnIndices[ COLUMN_INDICES] 0-based column index to parse.
OUTPUTS: double columnValues[  COLUMN_INDICES ]  Parsed column values.
         Note note                               Description of sensor location
RETURNS: int 1 if line values are in subset, else 0.
******************************************************************************/

static int parseColumnValues( char* line,
                              const long long seconds1,
                              const long long seconds2,
                              const double west,
                              const double east,
                              const double south,
                              const double north,
                              const int out_in_flag,
                              const int sensor,
                              const int isChannelDependent,
                              const int isChannelA,
                              const int isChannelB,
                              const int isPM25Corrected,
                              const int isTemperature,
                              const double maximumDifference,
                              const double maximumRatio,
                              const double defaultHumidity,
                              const char* const units,
                              const double validMinimum,
                              const double validMaximum,
                              const int columnIndices[],
                              double columnValues[],
                              Note note ) {
  int inSubset = 1;

  assert( line );
  assert( seconds1 > 0 );
  assert( seconds2 > seconds1 );
  assert( IN_RANGE( west, -180.0, 180.0 ) );
  assert( IN_RANGE( east, west, 180.0 ) );
  assert( IN_RANGE( south, -90.0, 90.0 ) );
  assert( IN_RANGE( north, south, 90.0 ) );
  assert( IN_RANGE( out_in_flag, 0, 2 ) );
  assert( isChannelDependent == 0 || isChannelDependent == 1 );
  assert( isChannelA == 0 || isChannelA == 1 );
  assert( isChannelB == 0 || isChannelB == 1 );
  assert( isChannelA + isChannelB <= 1 );
  assert( isPM25Corrected == 0 || isPM25Corrected == 1 );
  assert( isTemperature == 0 || isTemperature == 1 );
  assert( ! ( isChannelDependent && isTemperature ) );
  assert( sensor >= 0 );
  assert( IN_RANGE( maximumDifference, 0.0, 100.0 ) );
  assert( IN_RANGE( maximumRatio, 0.0, 1.0 ) );
  assert( defaultHumidity == MISSING_VALUE ||
          IN_RANGE( defaultHumidity, 0.0, 100.0 ) );
  assert( units ); assert( *units );
  assert( columnValues );
  assert( columnIndices );
  assert( note );

  {
    /* Hoist loop-invariant expressions for performance: */

    const int idColumn           = columnIndices[ ID_INDEX ];
    const int descriptionColumn  = columnIndices[ DESCRIPTION_INDEX ];
    const int insideColumn       = columnIndices[ INSIDE_INDEX ];
    const int latitudeColumn     = columnIndices[ LATITUDE_INDEX ];
    const int longitudeColumn    = columnIndices[ LONGITUDE_INDEX ];
    const int elevationColumn    = columnIndices[ ELEVATION_INDEX ];
    const int timestampColumn    = columnIndices[ TIMESTAMP_INDEX ];
    const int channelStateColumn = columnIndices[ CHANNEL_STATE_INDEX ];
    const int channelFlagsColumn = columnIndices[ CHANNEL_FLAGS_INDEX ];
    const int humidityColumn     = columnIndices[ HUMIDITY_INDEX ];
    const int pm25CF1AColumn     = columnIndices[ PM25_CF_1_A_INDEX ];
    const int pm25CF1BColumn     = columnIndices[ PM25_CF_1_B_INDEX ];
    const int pm25ATMAColumn     = columnIndices[ PM25_ATM_A_INDEX ];
    const int pm25ATMBColumn     = columnIndices[ PM25_ATM_B_INDEX ];
    const int variableColumn     = columnIndices[ VARIABLE_INDEX ];
    double* const idValue           = columnValues + ID_INDEX;
    double* const latitudeValue     = columnValues + LATITUDE_INDEX;
    double* const longitudeValue    = columnValues + LONGITUDE_INDEX;
    double* const elevationValue    = columnValues + ELEVATION_INDEX;
    double* const timestampValue    = columnValues + TIMESTAMP_INDEX;
    double* const humidityValue     = columnValues + HUMIDITY_INDEX;
    double* const pm25CF1AValue     = columnValues + PM25_CF_1_A_INDEX;
    double* const pm25CF1BValue     = columnValues + PM25_CF_1_B_INDEX;
    double* const pm25ATMAValue     = columnValues + PM25_ATM_A_INDEX;
    double* const pm25ATMBValue     = columnValues + PM25_ATM_B_INDEX;
    double* const variableValue     = columnValues + VARIABLE_INDEX;
    double timestamp  = MISSING_VALUE;
    double longitude  = MISSING_VALUE;
    double latitude   = MISSING_VALUE;
    double elevation  = MISSING_VALUE;
    double measure    = MISSING_VALUE;
    double humidity   = defaultHumidity;
    double pm25CF1A   = MISSING_VALUE;
    double pm25CF1B   = MISSING_VALUE;
    double pm25ATMA   = MISSING_VALUE;
    double pm25ATMB   = MISSING_VALUE;
    char* value = 0;
    char* endValue = 0;
    int column= 0;

    for ( column = 0; column < COLUMN_INDICES; ++column ) {
      columnValues[ column ] = MISSING_VALUE;
    }

    memset( note, 0, sizeof (Note) );
    note[ 0 ] = '-';

    for ( value = strtok_r( line, ",", &endValue ), column = 0;
          value && inSubset;
          value = strtok_r( 0, ",", &endValue ), ++column ) {
      const int isNull =
        value[ 0 ] == 'n' &&
        value[ 1 ] == 'u' &&
        value[ 2 ] == 'l' &&
        value[ 3 ] == 'l' &&
        value[ 4 ] == '\0';

      if ( ! isNull ) {

        if ( column == idColumn ) {
          long long id = 0;
          inSubset = parseLongLong( value, 1, LLONG_MAX, &id );

          if ( inSubset ) {
            *idValue = id;
            inSubset = sensor == 0 || *idValue == sensor;
          }

        } else if ( column == descriptionColumn ) {
          strncpy( note, value, sizeof (Note) / sizeof (char) - 1 );
        } else if ( column == insideColumn ) {
          inSubset = out_in_flag == 2 || *value - '0' == out_in_flag;
        } else if ( column == latitudeColumn ) {
          inSubset = parseDouble( value, south, north, latitudeValue );
          latitude = *latitudeValue;
        } else if ( column == longitudeColumn ) {
          inSubset = parseDouble( value, west, east, longitudeValue );
          longitude = *longitudeValue;
        } else if ( column == elevationColumn ) {
          inSubset =
            parseDouble( value,
                         METERS_TO_FEET * MINIMUM_VALID_ELEVATION_METERS,
                         METERS_TO_FEET * MAXIMUM_VALID_ELEVATION_METERS,
                         elevationValue );
          *elevationValue *= FEET_TO_METERS;
          elevation = *elevationValue;
        } else if ( column == timestampColumn ) {
          long long seconds = 0;
          inSubset = parseLongLong( value, seconds1, seconds2, &seconds );
          timestamp = seconds;
          *timestampValue = timestamp;
        } else if ( USE_SIGMOID && isPM25Corrected &&
                    column == pm25CF1AColumn ) {
          inSubset = parseDouble( value, MINIMUM_VALID_PM, MAXIMUM_VALID_PM,
                                  &pm25CF1A );
          *pm25CF1AValue = pm25CF1A;
        } else if ( USE_SIGMOID && isPM25Corrected &&
                    column == pm25CF1BColumn ) {
          inSubset = parseDouble( value, MINIMUM_VALID_PM, MAXIMUM_VALID_PM,
                                  &pm25CF1B );
          *pm25CF1BValue = pm25CF1B;
        } else if ( ! USE_SIGMOID && isPM25Corrected &&
                    column == pm25ATMAColumn ) {
           inSubset =
            parseDouble( value, MINIMUM_VALID_PM, MAXIMUM_VALID_PM,
                         &pm25ATMA );
          *pm25ATMAValue = pm25ATMA;
        } else if ( ! USE_SIGMOID && isPM25Corrected &&
                    column == pm25ATMBColumn ) {
           inSubset =
            parseDouble( value, MINIMUM_VALID_PM, MAXIMUM_VALID_PM,
                         &pm25ATMB );
          *pm25ATMBValue = pm25ATMB;
        } else if ( column == humidityColumn ) {
          inSubset =
            parseDouble( value, MINIMUM_VALID_HUMIDITY, MAXIMUM_VALID_HUMIDITY,
                         &humidity );

          /*
           * If computing pm25_corrected and humidity is missing/invalid then
           * just use defaultHumidity.
           */

          if ( ! inSubset && isPM25Corrected ) {
            inSubset = 1;
            humidity = defaultHumidity;
          }

          *humidityValue = humidity;
        }

        if ( isChannelDependent ) { /* Check channelState, channelFlags: */

          if ( applyChannelState && column == channelStateColumn ) {
            const char v = *value;

            /* 0 = No PM, 1 = PM_A, 2 = PM_B, 3 = PM_A + PM_B. */

            if ( isChannelA ) {
              inSubset = v == '3' || v == '1'; /* Both channels or just A ok */
            } else if ( isChannelB ) {
              inSubset = v == '3' || v == '2'; /* Both channels or just B ok */
            } else {
              inSubset = v == '3'; /* Both PM channels ok. */
            }
          } else if ( applyChannelFlag && column == channelFlagsColumn ) {
            const char v = *value;

            /*
             * 0 = Normal, 1 = PM_A downgraded, 2 = PM_B downgraded,
             * 3 = Both PM_A and PM_B downgraded.
             */

            if ( isChannelA ) {
              inSubset = v == '0' || v == '2'; /* None downgraded or just B */
            } else if ( isChannelB ) {
              inSubset = v == '0' || v == '1'; /* None downgraded or just A */
            } else {
              inSubset = v == '0'; /* Neither PM channels are downgraded. */
            }
          }
        }

        if ( inSubset && column == variableColumn ) {
          const double minimum = isTemperature ? -100.0 : validMinimum;
          const double maximum = isTemperature ?  212.0 : validMaximum;
          inSubset = parseDouble( value, minimum, maximum, &measure );

          if ( isTemperature ) {
            measure = ( measure - 32.0 ) * ( 5.0 / 9.0 ); /* F to C */
          }

          *variableValue = measure;
        }
      }
    } /* End of parsing columns. */

    /* Check if data line is in the subset and compute corrected PM25: */

    if ( inSubset ) {

      /* Check that required column values were parsed: */

      inSubset =
        *idValue > 0 &&
        IN_RANGE( longitude, west, east ) &&
        IN_RANGE( latitude, south, north ) &&
        IN_RANGE( timestamp, seconds1, seconds2 ) &&
        measure != MISSING_VALUE;

      if ( inSubset ) {

        /* Elevation is allowed to be missing. In this case set it to 0: */

        if ( elevation == MISSING_VALUE ) {
          *elevationValue = 0.0;
        }

        /* If specified then compute PM2.5 corrected value: */

        if ( isPM25Corrected ) {
#if USE_SIGMOID
          const double pm25Corrected =
              pm25_corrected_sigmoid( pm25CF1A, pm25CF1B, humidity,
                                      maximumDifference, maximumRatio );
#else
          const double pm25Corrected =
              pm25_corrected_piecewise( pm25ATMA, pm25ATMB, humidity,
                                        maximumDifference, maximumRatio );
#endif
          *variableValue = measure = pm25Corrected;
          inSubset = pm25Corrected != MISSING_VALUE;
        }

        /* Filter-out negative variable values depending on units: */

        if ( inSubset && measure < 0.0 ) {
          inSubset =
            ( isPM25Corrected && allowNegativePM25Corrected ) ||
            ( units[ 0 ] == 'd' &&
              units[ 1 ] == 'e' &&
              units[ 2 ] == 'g' &&
              units[ 3 ] == '\0' ) ||
            ( units[ 0 ] == 'm' && units[ 1 ] == '\0' ) ||
            ( units[ 0 ] == 'C' && units[ 1 ] == '\0' );
        }
      }
    }
  }

  return inSubset;
}



/******************************************************************************
PURPOSE: aggregateData - Compute per-id means of temp subset data extracted
         thus far.
INPUTS:  Data* data  data->tempFiles[ TEMP_FILE_2 ]  id-time sorted temp data
                                                     to aggregate.
                     data->points  Previous number of output points.
OUTPUTS: Data* data  data->tempFiles[ TEMP_AGGREGATED_FILE ]
                     Aggregated data appended.
                     data->points  Increased number of output points.
                     data->ok = 1 if successful, else 0.
******************************************************************************/

static void aggregateData( Data* data ) {
  assert( data );
  assert( data->tempFiles[ TEMP_FILE_2 ] == 0 );
  assert( data->tempFiles[ TEMP_AGGREGATED_FILE ] );
  data->tempFiles[ TEMP_FILE_2 ] =
    fopen( data->tempFileNames[ TEMP_FILE_2 ], "rb" );
  data->ok = data->tempFiles[ TEMP_FILE_2 ] != 0;

  if ( data->ok ) {
    FILE* const input = data->tempFiles[ TEMP_FILE_2 ];
    FILE* const output = data->tempFiles[ TEMP_AGGREGATED_FILE ];
    enum { LINE_LENGTH = 255 };
    char line[ LINE_LENGTH + 1 ] = "";
    int id           = -1;
    long long timestamp = 0;
    double longitude = 0.0;
    double latitude  = 0.0;
    double elevation = 0.0;
    double measure   = 0.0;
    double sum = 0.0;
    long long count = 0;
    long long minimumCount = 0;

    DEBUG( fprintf( stderr, "start aggregateData() data->points = %lu\n",
                    data->points ); )

    memset( line, 0, sizeof line );
    memset( data->note, 0, sizeof data->note );

    /* Read each line and parse column values comparing id: */

    while ( data->ok &&
            fgets( line, sizeof line, input ) != 0 ) {
      const int thisId = atoi( line ); /* atoi() is faster than sscanf(). */

      if ( thisId == id ) { /* Handle most frequent case first: aggregate sum*/

        /* Parse measure which appears after penultimate comma: */

        char* comma = strrchr( line, ',' );
        int commaCount = comma != 0;

        while ( comma > line && commaCount < 2 ) {
          --comma;

          if ( *comma == ',' ) {
            ++commaCount;
          }
        }

        data->ok = commaCount == 2 && comma > line;

        if ( data->ok ) {
          measure = atof( comma + 1 );
          sum += measure; /* Aggregate. */
          ++count;
        }
      } else { /* thisId != id */

        if ( id != -1 && count > 0 && count >= minimumCount ) {

          /* Output mean: */

          const double mean = sum / count;
          data->ok =
            fprintf( output, "%d,%lld,%lf,%lf,%lf,%d,%lf,%s\n",
                     id, timestamp, longitude, latitude, elevation,
                     (int) count, mean,
                     data->note ) > 8;
          data->points += data->ok;
        }

        /* Re-initialize: */

        data->ok =
          sscanf( line, "%d,%lld,%lf,%lf,%lf,%lld,%lf",
                  &id, &timestamp,
                  &longitude, &latitude, &elevation, &count, &measure ) == 7;
        sum = measure;
        count = 1;

        if ( minimumCount == 0 ) {
          minimumCount = computeMinimumCount( data, timestamp );
          DEBUG( fprintf( stderr, "  minimumCount = %lld\n", minimumCount ); )
        }

        /* Parse description (appearing after last comma) into data->note: */

        {
          char* comma = strrchr( line, ',' );
          data->ok = comma != 0;

          if ( data->ok ) {
            char* newline = 0;
            strncpy( data->note, comma + 1,
                     sizeof data->note / sizeof data->note[ 0 ] - 1 );
            newline = strchr( data->note, '\n' );

            if ( newline ) {
              *newline = '\0';
            }
          }
        }
      }
    }

    DEBUG( fprintf( stderr, "  count = %lld\n", count ); )

    /* Write final mean: */

    if ( data->ok && count > 0 && count >= minimumCount ) {
      const double mean = sum / count;
      data->ok =
        fprintf( output, "%d,%lld,%lf,%lf,%lf,%d,%lf,%s\n",
                 id, timestamp, longitude, latitude, elevation,
                 (int) count, mean,
                 data->note ) > 7;
      data->points += data->ok; /* Update number of output points. */
    }
  }

  DEBUG( fprintf( stderr,
                  "end aggregateData() data->points = %lu, data->ok = %d\n",
                  data->points, data->ok ); )
}



/******************************************************************************
PURPOSE: computeMinimumCount - Compute the minimum required count for the
         given aggregation period.
INPUTS:  const Data* const data->arguments.minimumAggregationCountPercentage.
                           data->arguments.aggregate.
                           data->arguments.yyyymmddhhmmss[].
         const long long timestamp  Data timestamp (seconds since 1970).
RETURNS: long long minimum required count.
NOTES:   Ensures global previousSeconds is initialized.
******************************************************************************/

static long long computeMinimumCount( const Data* const data,
                                      const long long timestamp ) {

  long long result = 0;
  assert( data );
  assert( data->arguments.aggregate );
  assert( timestamp > 0 );

  /* If not yet initialized then initialize just once now: */

  if ( previousSeconds == 0.0 ) {
    previousSeconds = secondsSince1970( YYYYMMDD_PREVIOUS * 1000000LL );
  }

  {
    const int secondsPerValue =
      timestamp <= previousSeconds ? PREVIOUS_SECONDS_PER_VALUE
      : SECONDS_PER_VALUE;
    const Arguments* const arguments = &data->arguments;
    const int aggregate = arguments->aggregate;
    enum {
      SECONDS_PER_MINUTE = 60,
      MINUTES_PER_HOUR   = 60,
      HOURS_PER_DAY      = 24,
      SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR,
      SECONDS_PER_DAY  = SECONDS_PER_HOUR * HOURS_PER_DAY
    };
    long long secondsPerAggregationPeriod = 0;

    /* Compute secondsPerAggregationPeriod: */

    if ( aggregate == AGGREGATE_HOURLY ) {
      secondsPerAggregationPeriod = SECONDS_PER_HOUR;
    } else if ( aggregate == AGGREGATE_DAILY ) {
      secondsPerAggregationPeriod = SECONDS_PER_DAY;
    } else if ( aggregate == AGGREGATE_MONTHLY ) {
      const long long yyyymmddhhmmss = secondsToYYYYMMDDHHMMSS( timestamp );
      const int yyyymm   = yyyymmddhhmmss / 100000000;
      const int yyyy = yyyymm / 100;
      const int mm   = yyyymm % 100;
      secondsPerAggregationPeriod = SECONDS_PER_DAY * daysInMonth( yyyy, mm );
    } else {
      const long long firstSeconds =
        secondsSince1970( arguments->yyyymmddhhmmss[ 0 ] );
      const long long lastSeconds =
        secondsSince1970( arguments->yyyymmddhhmmss[ 1 ] );
      secondsPerAggregationPeriod = 1 + lastSeconds - firstSeconds;
    }

    {
      const long long maximumValuesPerAggregationPeriod =
        secondsPerAggregationPeriod / secondsPerValue;
      const double fraction =
        arguments->minimumAggregationCountPercentage / 100.0;
      result = (long long) (maximumValuesPerAggregationPeriod * fraction + 0.5);
    }
  }

  assert( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: reformatData - Reformat final subset temp file to output format.
INPUTS:  Data* data  data->tempFiles[ TEMP_FILE_2 ]  id-time sorted temp data
                                                     to reformat or
                     data->tempFiles[ TEMP_AGGREGATED_FILE ]
                     Aggregated data to reformat.
                     data->points  Number of output points.
OUTPUTS: Data* data  data->tempFiles[ TEMP_FILE_1 ] if ascii else
                     data->tempFiles[ VARIABLES ] if xdr.
                     Reformatted data.
                     data->ok = 1 if successful, else 0.
NOTES: If data was aggregated then read from the temp aggregated file
       else read from temp file 2.
       Write intermediate unsorted reformatted file to temp file 1.
       Then sort temp file 1 into temp file 2.
       If format is ascii, the final result is temp file 2,
       else for xdr format the final results will be in per-variable temp files
******************************************************************************/

static void reformatData( Data* data ) {
  FILE* input = 0;
  int inputIndex = 0;
  assert( data );
  assert( data->ok );
  assert( data->points > 0 );
  assert( data->tempFiles[ TEMP_FILE_1 ] == 0 );
  assert( data->tempFiles[ TEMP_FILE_2 ] == 0 );
  assert( data->tempFiles[ TEMP_AGGREGATED_FILE ] == 0 );

  inputIndex =
    data->arguments.aggregate == AGGREGATE_NONE ? TEMP_FILE_2
    : TEMP_AGGREGATED_FILE;

  input =
    data->tempFiles[ inputIndex ] =
    fopen( data->tempFileNames[ inputIndex ], "rb" );

  data->ok = input != 0;

  if ( input ) {
    FILE* output =
      data->tempFiles[ TEMP_FILE_1 ] =
      fopen( data->tempFileNames[ TEMP_FILE_1 ], "wb" );
    data->ok = output != 0;

    if ( output ) {
      const int isASCII = data->arguments.format == FORMAT_ASCII;
      const int isAggregated = data->arguments.aggregate != AGGREGATE_NONE;
      enum { LINE_LENGTH = 255 };
      char line[ LINE_LENGTH + 1 ] = "";
      size_t points = 0;
      int id           = -1;
      long long timestamp = 0;
      double longitude = 0.0;
      double latitude  = 0.0;
      double elevation = 0.0;
      double measure   = 0.0;
      int count        = 0;
      long long yyyymmddhhmmss = 0;
      long long timestamp0 = 0; /* Memoize previous value for performance. */
      memset( line, 0, sizeof line );
      memset( data->note, 0, sizeof data->note );

      /* Read each line: */

      while ( data->ok &&
              fgets( line, sizeof line, input ) != 0 ) {

        data->ok =
          sscanf( line, "%d,%lld,%lf,%lf,%lf,%d,%lf",
                  &id, &timestamp, &longitude, &latitude, &elevation,
                  &count, &measure ) == VARIABLES;

        /* Parse description (appearing after last comma) into data->note: */

        if ( data->ok ) {
          char* comma = strrchr( line, ',' );
          data->ok = comma != 0;

          if ( data->ok ) {
            char* newline = 0;
            strncpy( data->note, comma + 2, /* Skip space replacing quote. */
                     sizeof data->note / sizeof data->note[ 0 ] - 1 );
            newline = strchr( data->note, '\n' );

            if ( newline ) {
              *newline = '\0';
            }

            /* Write reformatted line to temp file 1: */

            if ( timestamp != timestamp0 ) { /* Check to avoid expensive call*/
              yyyymmddhhmmss = secondsToYYYYMMDDHHMMSS( timestamp );
              timestamp0 = timestamp;
            }

            if ( isASCII ) {
              const int yyyy = yyyymmddhhmmss / 10000000000LL;
              const int mm   = yyyymmddhhmmss / 100000000 % 100;
              const int dd   = yyyymmddhhmmss / 1000000 % 100;
              const int hh   = yyyymmddhhmmss / 10000 % 100;
              const int mm2  = yyyymmddhhmmss / 100 % 100;
              const int ss   = yyyymmddhhmmss % 100;

              if ( isAggregated ) {
                data->ok =
                  fprintf( output,
                           "%d-%02d-%02dT%02d:%02d:%02d-0000"
                           "\t%f\t%f\t%f\t%d\t%d\t%f\t%s\n",
                           yyyy, mm, dd, hh, mm2, ss,
                           longitude, latitude, elevation, id, count, measure,
                           data->note ) > 12;
              } else {
                data->ok =
                  fprintf( output,
                           "%d-%02d-%02dT%02d:%02d:%02d-0000"
                           "\t%f\t%f\t%f\t%d\t%f\t%s\n",
                           yyyy, mm, dd, hh, mm2, ss,
                           longitude, latitude, elevation, id, measure,
                           data->note ) > 12;
              }

            } else {
              data->ok =
                fprintf( output,
                         "%lld,%f,%f,%f,%d,%d,%f,%s\n",
                         yyyymmddhhmmss,
                         longitude, latitude, elevation, id, count, measure,
                         data->note ) > 12;
            }

            points += data->ok;
          }
        }
      }

      data->ok = points == data->points;

      if ( data->ok ) {

        /* Sort temp file 1 into temp file 2 and empty temp file 1: */

        sortTempData( data );

        if ( data->ok && ! isASCII ) { /* Read temp 2. Write per-var temps: */
          writeTempVariableFiles( data );
        }
      }
    }
  }

  closeTempFiles( data );
}



/******************************************************************************
PURPOSE: writeTempVariableFiles - Read temp file 2 and write per-variable
         temp files in big-endian binary format.
INPUTS:  Data* data  data->tempFiles[ TEMP_FILE_2 ]  time sorted temp data
                                                     to read.
                     data->points  Number of output points.
OUTPUTS: Data* data  data->tempFiles[ VARIABLES ]  Temp files to write.
                     data->ok = 1 if successful, else 0.
******************************************************************************/

static void writeTempVariableFiles( Data* data ) {
  enum { LINE_LENGTH = 255 };
  char line[ LINE_LENGTH + 1 ] = "";
  double longitude = 0.0;
  double latitude  = 0.0;
  double elevation = 0.0;
  double measure   = 0.0;
  int count = 0;
  long long yyyymmddhhmmss = 0;
  int id = 0;
  size_t points = 0;

  assert( data );
  assert( data->ok );
  assert( data->points > 0 );
  assert( data->tempFiles[ TEMP_FILE_1 ] );

  memset( line, 0, sizeof line );
  memset( data->note, 0, sizeof data->note );

  /* Read each line: */

  openVariableTempFiles( data );

  if ( data->ok ) {
    FILE* notesOutput = data->tempFiles[ TEMP_FILE_1 ];
    FILE* input =
      data->tempFiles[ TEMP_FILE_2 ] =
      fopen( data->tempFileNames[ TEMP_FILE_2 ], "rb" );
      data->ok = input != 0;

    if ( input ) {
      FILE** tempFiles = data->tempFiles;

      /* Read each line: */

      while ( data->ok &&
              fgets( line, sizeof line, input ) != 0 ) {
        data->ok =
          sscanf( line, "%lld,%lf,%lf,%lf,%d,%d,%lf",
                  &yyyymmddhhmmss,
                  &longitude, &latitude, &elevation,
                  &id, &count, &measure ) == 7;

        /* Parse description (after last comma) into data->note: */

        if ( data->ok ) {
          char* comma = strrchr( line, ',' );
          data->ok = comma != 0;

          if ( data->ok ) {
            double words[ VARIABLES ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
            char* newline = 0;
            strncpy( data->note, comma + 1,
                     sizeof data->note / sizeof data->note[ 0 ] - 1 );
            newline = strchr( data->note, '\n' );

            if ( newline ) {
              *newline = '\0';
            }

            /* Write note to temp file 1: */

            data->ok = fprintf( notesOutput, "%-79s\n", data->note ) == 80;

            if ( data->ok ) { /* Write data to per-variable temp files: */
              words[ 0 ] = yyyymmddhhmmss;
              words[ 1 ] = longitude;
              words[ 2 ] = latitude;
              words[ 3 ] = elevation;
              words[ 4 ] = id;
              words[ 5 ] = count;
              words[ 6 ] = measure;
              rotate8ByteArrayIfLittleEndian( words, VARIABLES );
              data->ok =
                fwrite( &words[ 0 ], 8, 1, tempFiles[ 0 ] ) == 1 &&
                fwrite( &words[ 1 ], 8, 1, tempFiles[ 1 ] ) == 1 &&
                fwrite( &words[ 2 ], 8, 1, tempFiles[ 2 ] ) == 1 &&
                fwrite( &words[ 3 ], 8, 1, tempFiles[ 3 ] ) == 1 &&
                fwrite( &words[ 4 ], 8, 1, tempFiles[ 4 ] ) == 1 &&
                fwrite( &words[ 5 ], 8, 1, tempFiles[ 5 ] ) == 1 &&
                fwrite( &words[ 6 ], 8, 1, tempFiles[ 6 ] ) == 1;
              points += data->ok;
            }
          }
        }
      }
    }
  }

  data->ok = points == data->points;
}



/******************************************************************************
PURPOSE: sortTempData - Sort temp file 1 into temp file 2 then empty temp 1.
INPUTS:  Data* data  data->tempFileNames[ TEMP_FILE_1 ]  Name of file to sort.
                     data->tempFileNames[ TEMP_FILE_2 ]  Name of sorted file.
OUTPUTS: Data* data  data->tempFileNames[ TEMP_FILE_1 ]  Empty file.
                     data->tempFileNames[ TEMP_FILE_2 ]  Sorted closed file.
                     data->tempFiles[ TEMP_FILE_1 ]  Rewound.
                     data->tempFiles[ TEMP_FILE_2 ]  Closed.
                     data->ok = 1 if successful, else 0.
NOTES:   HACK: invokes UNIX /usr/bin/sort!
******************************************************************************/

static void sortTempData( Data* data ) {
  char command[ 256 ] = "";
  assert( data );
  assert( data->tempFileNames ); assert( data->tempFileNames[ 0 ] );

  /* Close temp files to flush content: */

  if ( data->tempFiles[ TEMP_FILE_1 ] ) {
    fclose( data->tempFiles[ TEMP_FILE_1 ] );
    data->tempFiles[ TEMP_FILE_1 ] = 0;
  }

  if ( data->tempFiles[ TEMP_FILE_2 ] ) {
    fclose( data->tempFiles[ TEMP_FILE_2 ] );
    data->tempFiles[ TEMP_FILE_2 ] = 0;
  }

  /* Sort temp file 1 into temp file 2 (overwriting it): */

  memset( command, 0, sizeof command );
  snprintf( command, sizeof command / sizeof *command - 1,
            "/usr/bin/sort -T %s -o %s %s",
            data->arguments.tmpdir,
            data->tempFileNames[ TEMP_FILE_2 ],
            data->tempFileNames[ TEMP_FILE_1 ] );
  data->ok = system( command ) == 0;

  if ( data->ok ) { /* Remove temp file 1 then reopen it for more appending: */
    data->ok = unlink( data->tempFileNames[ TEMP_FILE_1 ] ) == 0;
    data->tempFiles[ TEMP_FILE_1 ] =
      fopen( data->tempFileNames[ TEMP_FILE_1 ], "wb" );
    data->ok = data->tempFiles[ TEMP_FILE_1 ] != 0;
  } else {
    fprintf( stderr, "\nFailed to sort file %s into %s\n",
             data->tempFileNames[ TEMP_FILE_1 ],
             data->tempFileNames[ TEMP_FILE_2 ] );
  }
}



/******************************************************************************
PURPOSE: streamData - Write final content of temp files to stdout.
INPUTS:  Data* const data  Data to write.
OUTPUTS: Data* const data  data->ok, tempFiles[] = 0 (closed and removed).
******************************************************************************/

static void streamData( Data* const data ) {

  assert( data ); assert( data->ok ); assert( data->points > 0 );
  assert( data->outputBuffer ); assert( data->bufferSize );

  /* Temp file names are initialized but temp files are closed after writing:*/

  assert( data->tempFileNames[ 0 ][ 0 ] );
  assert( data->tempFiles[     0 ] == 0 );
  assert( data->tempFileNames[ 1 ][ 0 ] );
  assert( data->tempFiles[     1 ] == 0 );
  assert( data->tempFileNames[ 2 ][ 0 ] );
  assert( data->tempFiles[     2 ] == 0 );
  assert( data->tempFileNames[ 3 ][ 0 ] );
  assert( data->tempFiles[     3 ] == 0 );

  {
    const size_t bytes = data->bufferSize;
    void* buffer = data->outputBuffer;
    const int isASCII = data->arguments.format == FORMAT_ASCII;
    const int first = isASCII ? TEMP_FILE_2 : -1; /* -1 = notes first. */
    const int last  = isASCII ? first : VARIABLES - 1;
    const int isAggregated = data->arguments.aggregate != AGGREGATE_NONE;
    int fileIndex = 0;
    int index = 0;
    streamHeader( data );

    for ( index = first; data->ok && index <= last; ++index ) {
      fileIndex = index < 0 ? TEMP_FILE_1 : index;

      if ( isAggregated || index != VARIABLES - 2 ) {
        data->tempFiles[ fileIndex ] =
          fopen( data->tempFileNames[ fileIndex ], "rb" );
        data->ok = data->tempFiles[ fileIndex ] != 0;

        if ( ! data->ok ) {
          fprintf( stderr, "\nCan't open temp data file '%s' for reading.\n",
                   data->tempFileNames[ fileIndex ] );
        } else {

          while ( data->ok && ! feof( data->tempFiles[ fileIndex ] ) ) {
            const size_t bytesRead =
              fread( buffer, 1, bytes, data->tempFiles[ fileIndex ] );

            if ( bytesRead ) {
              const size_t bytesWritten = fwrite(buffer, 1, bytesRead, stdout);
              data->ok = bytesWritten == bytesRead;
            }
          }
        }
      }
    }

    if ( ! data->ok ) {
      fprintf( stderr, "\nFailed to stream subset data from temp file '%s'.\n",
               data->tempFileNames[ fileIndex ] );
    }
  }
}



/******************************************************************************
PURPOSE: streamHeader - Write ASCII header of subset to stdout.
INPUTS:  Data* data  Data to write.
******************************************************************************/

static void streamHeader( const Data* const data ) {
  const Arguments* const arguments = &data->arguments;

  /* Append to name to indicate aggregation: */

  const char* const aggregation =
    arguments->aggregate == AGGREGATE_HOURLY ? "_hourly"
    : arguments->aggregate == AGGREGATE_DAILY ? "_daily"
    : arguments->aggregate == AGGREGATE_MONTHLY ? "_monthly"
    : arguments->aggregate == AGGREGATE_ALL ? "_mean"
    : "";

  if ( arguments->format == FORMAT_ASCII ) {
    const char* const countColumn = *aggregation ? "\tCOUNT(-)" : "";
    printf( "Timestamp(UTC)\tLONGITUDE(deg)\tLATITUDE(deg)\tELEVATION(m)\t"
            "STATION(-)%s\t%s%s(%s)\tNOTE\n",
            countColumn, arguments->variable, aggregation, data->units );
  } else {
    const long long yyyymmddhhmmss1 = arguments->yyyymmddhhmmss[ 0 ];
    const long long yyyymmddhhmmss2 = arguments->yyyymmddhhmmss[ 1 ];
    const int yyyy1 = yyyymmddhhmmss1 / 10000000000;
    const int mm1   = yyyymmddhhmmss1 / 100000000 % 100;
    const int dd1   = yyyymmddhhmmss1 / 1000000 % 100;
    const int hh1   = yyyymmddhhmmss1 / 10000 % 100;
    const int mm21  = yyyymmddhhmmss1 / 100 % 100;
    const int ss1   = yyyymmddhhmmss1 % 100;
    const int yyyy2 = yyyymmddhhmmss2 / 10000000000;
    const int mm2   = yyyymmddhhmmss2 / 100000000 % 100;
    const int dd2   = yyyymmddhhmmss2 / 1000000 % 100;
    const int hh2   = yyyymmddhhmmss2 / 10000 % 100;
    const int mm22  = yyyymmddhhmmss2 / 100 % 100;
    const int ss2   = yyyymmddhhmmss2 % 100;
    const char* const countName  = *aggregation ? "count " : "";
    const char* const countUnits = *aggregation ? "- " : "";
    printf( "Point 1.0\n"
            "%s\n"
            "%04d-%02d-%02dT%02d:%02d:%02d-0000 "
            "%04d-%02d-%02dT%02d:%02d:%02d-0000\n",
            arguments->description,
            yyyy1, mm1, dd1, hh1, mm21, ss1,
            yyyy2, mm2, dd2, hh2, mm22, ss2 );
    printf( "# Dimensions: variables points:\n%d %lu\n",
            *aggregation ? VARIABLES : VARIABLES - 1, data->points );
    printf( "# Variable names:\n" );
    printf( "timestamp longitude latitude elevation id %s%s%s\n",
            countName, arguments->variable, aggregation );
    printf( "# Variable units:\nyyyymmddhhmmss deg deg m - %s%s\n",
            countUnits, data->units );
    printf( "# char notes[points][80] and\n" );
    printf( "# IEEE-754 64-bit reals data[variables][points]:\n" );
  }
}




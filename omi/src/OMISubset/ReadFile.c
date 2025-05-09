/******************************************************************************
PURPOSE: ReadFile.c - Simple to use wrapper routines to read data from
                      OMI-AURA HDF5 files.

NOTES:   Uses HDF5 libraries and libraries they depend on (curl, z, dl).

HISTORY: 2017-12-06 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <stdlib.h> /* For malloc(), free(). */
#include <string.h> /* For strncpy(), strcmp(), memset(). */
#include <time.h>   /* For time(), gmtime_r(). */
#include <float.h>  /* For DBL_MAX. */

#ifdef _WIN32
/* gmtime_r() is missing, just take a chance on using gmtime(). */
static struct tm* gmtime_r( const time_t* clock, struct tm* result ) {
  const struct tm* const copy = gmtime( clock ); /* HACK: not thread-safe! */
  memcpy( result, copy, sizeof (struct tm) );
  return result;
}
#endif

#include <hdf5.h> /* For H5*. */

#include "ReadFile.h"  /* For public interface. */

/*=========================== FORWARD DECLARATIONS ==========================*/

static int openDataset( const int file, const char* const path,
                        const char* const variable );

static void closeDataset( const int dataset );

static int readDatasetDimensions( const int dataset,
                                  size_t* rows, size_t* columns );

static int readFileData( const int dataset,
                         const size_t rows,
                         const size_t columns,
                         double data[] );

static size_t filterDataByVariable( const int file,
                                    const size_t rows,
                                    const size_t columns,
                                    const char* const product,
                                    const char* const variable,
                                    const double dataMinimum,
                                    const double dataMaximum,
                                    double data[],
                                    double temp[] );

static const char* variablePath( const char* const product,
                                 const char* const variable );

static double readAttribute( const int dataset, const char* const name );

static long long toSecondsUTC70( const double secondsTAI93 );

static long long toUTC( const long long seconds );

static int isValidTimestamp( const long long yyyydddhhmm );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: openFile - Open HDF5 file for reading.
INPUTS:  const char* const fileName   Name of HDF5 file to open for reading.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int openFile( const char* const fileName ) {
  int result = 0;
  assert( fileName ); assert( *fileName );
  result = H5Fopen( fileName, 0, 0 );

  if ( result == -1 ) {
    fprintf( stderr, "\a\nFailed to open HDF5 file for reading: %s.\n",
             fileName );
  }

  return result;
}



/******************************************************************************
PURPOSE: closeFile - Close HDF5 file.
INPUTS:  const int file HDF5 file to close.
******************************************************************************/

void closeFile( const int file ) {
  H5Fclose( file );
}



/******************************************************************************
PURPOSE: readDimensions - Read dimensions of swath in file.
INPUTS:  const int file    ID of file to query.
         const char* const product   Name of file type: OMNO2, OMCLDRR, etc.
OUTPUTS: size_t* rows      Number of rows in swath.
         size_t* columns   Number of columns in swath.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readDimensions( const int file, const char* const product,
                    size_t* const rows, size_t* const columns ) {
  int result = 0;

  assert( file > -1 ); assert( product ); assert( *product );
  assert( rows ); assert( columns );
  *rows = *columns = 0;

  {
    const char* const dataPath =
      ! strcmp( product, "OMNO2" ) ?
        "/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields"
      : ! strcmp( product, "OMTO3" ) ?
        "/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields"
      : ! strcmp( product, "OMCLDRR" ) ?
        "/HDFEOS/SWATHS/Cloud Product/Geolocation Fields"
      : "";
    const int dataset = openDataset( file, dataPath, "Latitude" );

    if ( dataset > -1 ) {
      result = readDatasetDimensions( dataset, rows, columns );
      closeDataset( dataset );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readDataset - Read appropriately-filtered variable data.
INPUTS:  const int file              ID of file to query.
         const size_t rows           Rows of dataset to match.
         const size_t columns        Columns of dataset to match.
         const char* const product   Name of file type: OMNO2, OMCLDRR, etc.
         const char* const variable  Name of dataset variable.
         const double maximumCloudFraction  Maximum allowed cloud fraction.
         const double maximumSolarZenithAngle  Max allowed solar zenith angle.
         const int allowNegativeCounts  Allow negative molecules/cm2? UGLY.
OUTPUTS: char units[ 80 ]            Units of variable.
         double data[ dimension0 * dimension1 ]  Filtered variable data values.
         double temp[ dimension0 * dimension1 ]  Temp buffer for QC, etc.
RETURNS: size_t number of unfiltered points, else 0.
******************************************************************************/

size_t readDataset( const int file,
                    const size_t rows,
                    const size_t columns,
                    const char* const product,
                    const char* const variable,
                    const double maximumCloudFraction,
                    const double maximumSolarZenithAngle,
                    const int allowNegativeCounts,
                    char units[ 80],
                    double data[],
                    double temp[] ) {

  typedef struct {
    const char* product;         /* File type OMNO2, OMCLDRR, etc. */
    const char* variable;        /* Name of data variable to filter. */
    const char* units;           /* Name of data units. */
    double dataMinimum;          /* Minimum valid value for this variable. */
    double dataMaximum;          /* Maximum valid value for this variable. */
    const char* filterVariable1; /* Name of 1st data variable to filter by. */
    const char* filterVariable2; /* Name of 2nd data variable to filter by. */
  } Entry;

  static const Entry table[] = {

    /* OMNO2 */

    /* These variables are in coordinatesPath: */

    { "OMNO2", "Latitude", "deg", -90.0, 90.0, 0, 0 },
    { "OMNO2", "Longitude", "deg", -180.0, 180.0, 0, 0 },
    { "OMNO2", "SolarAzimuthAngle", "deg", -180.0, 180.0, 0, 0 },
    { "OMNO2", "SolarZenithAngle", "deg", 0.0,  180.0, 0, 0 },
    { "OMNO2", "Time", "YYYYDDDHHMM", 0.0, 1e20, 0, 0 },
    { "OMNO2", "ViewingAzimuthAngle", "deg", -180.0, 180.0, 0, 0 },
    { "OMNO2", "ViewingZenithAngle", "deg", 0.0, 90.0, 0, 0 },

    /* These variables are in dataPath: */

    { "OMNO2", "AMFQualityFlags", "-", 0.0, 65535.0, 0, 0 },
    { "OMNO2", "AlgorithmFlags", "-", 0.0, 65535.0, 0, 0 },

    { "OMNO2", "AmfStrat", "-", 0.0, 1e38, "XTrackQualityFlags", "AMFQualityFlags" },
    { "OMNO2", "AmfStratClear", "-", 0.0, 1e38, "XTrackQualityFlags", "AMFQualityFlags"},
    { "OMNO2", "AmfStratCloudy", "-", 0.0, 1e38, "XTrackQualityFlags","AMFQualityFlags"},
    { "OMNO2", "AmfStratStd", "-", 0.0, 1e38, "XTrackQualityFlags", "AMFQualityFlags" },

    { "OMNO2", "AmfTrop", "-", 0.0, 1e38, "XTrackQualityFlags", "AMFQualityFlags" },
    { "OMNO2", "AmfTropClear", "-", 0.0, 1e38, "XTrackQualityFlags", "AMFQualityFlags" },
    { "OMNO2", "AmfTropCloudy", "-", 0.0, 1e38, "XTrackQualityFlags", "AMFQualityFlags"},
    { "OMNO2", "AmfTropStd", "-", 0.0, 1e38, "XTrackQualityFlags", "AMFQualityFlags" },

    { "OMNO2", "CloudFraction", "-", 0.0, 1.0, "XTrackQualityFlags", 0 },
    { "OMNO2", "CloudFractionStd", "-", 0.0, 1.0, "XTrackQualityFlags", 0 },
    { "OMNO2", "CloudPressure", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMNO2", "CloudPressureStd", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMNO2", "CloudRadianceFraction", "-", 0.0, 1.0, "XTrackQualityFlags", 0 },

    { "OMNO2", "ColumnAmountNO2", "molecules/cm2", 0.0, 1e38, "XTrackQualityFlags",
      "VcdQualityFlags" },
    { "OMNO2", "ColumnAmountNO2Std", "molecules/cm2", 0.0, 1e38, "XTrackQualityFlags",
      "VcdQualityFlags" },
    { "OMNO2", "ColumnAmountNO2Strat", "molecules/cm2", 0.0, 1e38, "XTrackQualityFlags",
      "VcdQualityFlags" },
    { "OMNO2", "ColumnAmountNO2StratStd", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "ColumnAmountNO2Trop", "molecules/cm2", 0.0, 1e38, "XTrackQualityFlags",
      "VcdQualityFlags" },
    { "OMNO2", "ColumnAmountNO2TropStd", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },

    { "OMNO2", "ScdApStrat", "molecules/cm2", 0.0, 1e38, "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "ScdApTrop", "molecules/cm2", 0.0, 1e38, "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "SceneLER", "-", 0.0, 1e38, "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "ScenePressure", "hPa", 0.0, 2000.0, "XTrackQualityFlags",
      "VcdQualityFlags" },
    { "OMNO2", "SlantColumnAmountCHOCHO", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "SlantColumnAmountCHOCHOStd", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "SlantColumnAmountH2O", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "SlantColumnAmountH2OStd", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "SlantColumnAmountNO2", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "SlantColumnAmountNO2Destriped", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "SlantColumnAmountNO2Std", "molecules/cm2", 0.0, 1e38,
      "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "TerrainHeight", "m", -500.0, 10000.0, "XTrackQualityFlags", 0 },
    { "OMNO2", "TerrainPressure", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMNO2", "TerrainReflectivity", "-", 0.0, 1e38, "XTrackQualityFlags", 0 },
    { "OMNO2", "TropopausePressure", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMNO2", "VcdApBelowCloud", "-", 0.0, 1e38, "XTrackQualityFlags",
      "VcdQualityFlags" },
    { "OMNO2", "VcdApStrat", "-", 0.0, 1e38, "XTrackQualityFlags", "VcdQualityFlags" },
    { "OMNO2", "VcdApTrop", "-", 0.0, 1e38, "XTrackQualityFlags", "VcdQualityFlags" },

    { "OMNO2", "VcdQualityFlags", "-", 0.0, 65535.0, 0, 0 },
    { "OMNO2", "XTrackQualityFlags", "-", 0.0, 255.0, 0, 0 },

    /* OMTO3 */

    { "OMTO3", "GroundPixelQualityFlags", "-", 0.0, 65535.0, 0, 0 },
    { "OMTO3", "LandWaterClassification", "-", 0.0, 255.0, 0, 0 },
    { "OMTO3", "Latitude", "deg", -90.0, 90.0, 0, 0 },
    { "OMTO3", "Longitude", "deg", -180.0, 180.0, 0, 0 },
    { "OMTO3", "RelativeAzimuthAngle", "deg", -180.0, 180.0, 0, 0 },
    { "OMTO3", "SolarAzimuthAngle", "deg", -180.0, 180.0, 0, 0 },
    { "OMTO3", "SolarZenithAngle", "deg", 0.0,  180.0, 0, 0 },
    { "OMTO3", "TerrainHeight", "m", -500.0, 10000.0, "XTrackQualityFlags", 0 },
    { "OMTO3", "Time", "YYYYDDDHHMM", 0.0, 1e20, 0, 0 },
    { "OMTO3", "ViewingAzimuthAngle", "deg", -180.0, 180.0, 0, 0 },
    { "OMTO3", "ViewingZenithAngle", "deg", 0.0, 90.0, 0, 0 },
    { "OMTO3", "WaterFraction", "%", 0.0, 100.0, 0, 0 },
    { "OMTO3", "XTrackQualityFlags", "-", 0.0, 255.0, 0, 0 },

    { "OMTO3", "AlgorithmFlags", "-", 0.0, 255.0, 0, 0 },
    { "OMTO3", "CloudPressure", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMTO3", "ColumnAmountO3", "DU", 0.0, 1000.0, "XTrackQualityFlags", 0 },
    { "OMTO3", "O3BelowCloud", "DU", 0.0, 1000.0, "XTrackQualityFlags", 0 },
    { "OMTO3", "QualityFlags", "-", 0.0, 65535.0, 0, 0 },
    { "OMTO3", "RadianceBadPixelFlagAccepted", "-", 0.0, 65535.0, 0, 0 },
    { "OMTO3", "RadiativeCloudFraction", "-", 0.0, 1.0, "XTrackQualityFlags", 0},
    { "OMTO3", "Reflectivity331", "%", -15.0, 115.0, "XTrackQualityFlags", 0},
    { "OMTO3", "Reflectivity360", "%", -15.0, 115.0, "XTrackQualityFlags", 0},
    { "OMTO3", "SO2index", "-", -300.0, 300.0, "XTrackQualityFlags", 0},
    { "OMTO3", "StepOneO3", "DU", 0.0, 1000.0, "XTrackQualityFlags", 0 },
    { "OMTO3", "StepTwoO3", "DU", 0.0, 1000.0, "XTrackQualityFlags", 0 },
    { "OMTO3", "TerrainPressure", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMTO3", "UVAerosolIndex", "-", -30.0, 30.0, "XTrackQualityFlags", 0},
    { "OMTO3", "fc", "-", 0.0, 1.0, "XTrackQualityFlags", 0},

    /* OMCLDRR */

    { "OMCLDRR", "GroundPixelQualityFlags", "-", 0.0, 65535.0, 0, 0 },
    { "OMCLDRR", "Latitude", "deg", -90.0, 90.0, 0, 0 },
    { "OMCLDRR", "Longitude", "deg", -180.0, 180.0, 0, 0 },
    { "OMCLDRR", "RelativeAzimuthAngle", "deg", -180.0, 180.0, 0, 0 },
    { "OMCLDRR", "SolarZenithAngle", "deg", 0.0,  180.0, 0, 0 },
    { "OMCLDRR", "TerrainHeight", "m", -500.0, 10000.0, "XTrackQualityFlags",0},
    { "OMCLDRR", "Time", "YYYYDDDHHMM", 0.0, 1e20, 0, 0 },
    { "OMCLDRR", "ViewingZenithAngle", "deg", 0.0, 90.0, 0, 0 },
    { "OMCLDRR", "XTrackQualityFlags", "-", 0.0, 255.0, 0, 0 },

    { "OMCLDRR", "Chlorophyll", "mg/m3", 0.0, 1e9, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "CloudFractionforO3", "-", 0.0, 1.0, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "CloudPressureforO3", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "CloudPressureforO3_uncorrected", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "Convergence_factor", "-", 0.0, 1e20, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "Filling-In", "-", 0.0, 1e20, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "ProcessingQualityFlagsforO3", "-", 0.0, 65535.0, 0, 0 },
    { "OMCLDRR", "RadiativeCloudFraction", "-", 0.0, 1.0, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "Reflectivity", "-", 0.0, 1.0, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "Residual_bias", "-", 0.0, 1e20, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "Residual_stddev", "-", 0.0, 1e20, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "SurfaceReflectivity", "-", 0.0, 1.0, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "TerrainPressure", "hPa", 0.0, 2000.0, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "WavelengthShift", "nm", 0.0, 1e20, "XTrackQualityFlags", 0 },
    { "OMCLDRR", "dIdR", "nm", 0.0, 1e20, "XTrackQualityFlags", 0 },

    { "", "", "", 0.0, 0.0, 0, 0 } /* End of table. */
  };

  const size_t entries = sizeof table / sizeof *table - 1;
  int result = 0;
  size_t index = 0;

  assert( file > -1 );
  assert( rows ); assert( columns );
  assert( product ); assert( *product );
  assert( variable ); assert( *variable );
  assert( maximumCloudFraction >= 0.0 ); assert( maximumCloudFraction <= 1.0 );
  assert( maximumSolarZenithAngle >= 0.0 );
  assert( maximumSolarZenithAngle <= 90.0 );
  assert( allowNegativeCounts == 0 || allowNegativeCounts == 1 );
  assert( units ); assert( data ); assert( temp );

  units[ 0 ] = '\0';

  /* Find variable in table. If not found it is an error: */

  while ( index < entries &&
          ( strcmp( table[ index ].product, product ) ||
            strcmp( table[ index ].variable, variable ) ) ) {
    ++index;
  }

  if ( index == entries ) {
    fprintf( stderr, "\nInvalid product variable '%s %s'.\n",
             product, variable );
  } else {
    const Entry* const entry = table + index;

    const char* const datasetPath = variablePath( product, variable );
    int dataset = openDataset( file, datasetPath, variable );

    if ( dataset > -1 ) {
      const int isTime = ! strcmp( variable, "Time" );
      result = readFileData( dataset, rows, columns, data );
      closeDataset( dataset ), dataset = -1;

      if ( result ) {
        const size_t points = rows * columns;
        const char* const filterVariable1 = entry->filterVariable1;
        const char* const filterVariable2 = entry->filterVariable2;
        strncpy( units, entry->units, 79 );
        assert( ! strchr( units, ' ' ) ); /* No spaces allowed in units. */

        if ( isTime ) { /* Convert from seconds TAI93 to yyyymmddhhmmss UTC: */

          for ( index = 0, result = 0; index < points; ++index ) {
            const double secondsTAI93 = temp[ index ];
            const long long secondsUTC70 = toSecondsUTC70( secondsTAI93 );
            const long long yyyydddhhmm = toUTC( secondsUTC70 );

            if ( isValidTimestamp( yyyydddhhmm ) ) {
              data[ index ] = yyyydddhhmm;
              ++result;
            } else {
              data[ index ] = MISSING_VALUE;
            }
          }
        }

        if ( result && filterVariable1 ) {
          result = filterDataByVariable( file, rows, columns,
                                         product, filterVariable1,
                                         0.0, 0.0, data, temp );
        }

        if ( result && filterVariable2 ) {
          result = filterDataByVariable( file, rows, columns,
                                         product, filterVariable2,
                                         0.0, 0.0, data, temp );
        }

        if ( result && maximumSolarZenithAngle != 90.0 ) {
          result = filterDataByVariable( file, rows, columns,
                                         product, "SolarZenithAngle",
                                         0.0, maximumSolarZenithAngle,
                                         data, temp );
        }

        if ( result && maximumCloudFraction != 1.0 ) {
          const char* const cloudFractionVariable =
            ! strcmp( product, "OMNO2" ) ? "CloudFraction"
            : "RadiativeCloudFraction";
          result = filterDataByVariable( file, rows, columns,
                                         product, cloudFractionVariable,
                                         0.0, maximumCloudFraction,
                                         data, temp );
        }

        /* Filter data by valid range: */

        if ( result ) {
          const double dataMinimum =
            allowNegativeCounts &&
            ! strcmp( entry->units, "molecules/cm2" ) ? -1e29
            : entry->dataMinimum;
          const double dataMaximum = entry->dataMaximum;

          for ( index = 0, result = 0; index < points; ++index ) {
            const double value = data[ index ];
            const int valid = value >= dataMinimum && value <= dataMaximum;

            if ( valid ) {
              ++result;
            } else {
              data[ index ] = MISSING_VALUE;
            }
          }
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: filterDataByVariable - Filter data by another variable.
INPUTS:  const int file              ID of file to query.
         const size_t rows           Rows of dataset to match.
         const size_t columns        Columns of dataset to match.
         const char* const product   Name of file type: OMNO2, OMCLDRR, etc.
         const char* const variable  Name of dataset variable.
         const double dataMinimum    Minimum valid filter variable value.
         const double dataMaximum    Maximum valid filter variable value.
OUTPUTS: double data[ dimension0 * dimension1 ]  Filtered variable data values.
         double temp[ dimension0 * dimension1 ]  Temp buffer for QC, etc.
RETURNS: size_t number of non-missing points, else 0.
******************************************************************************/

static size_t filterDataByVariable( const int file,
                                    const size_t rows,
                                    const size_t columns,
                                    const char* const product,
                                    const char* const variable,
                                    const double dataMinimum,
                                    const double dataMaximum,
                                    double data[],
                                    double temp[] ) {

  const char* const path = variablePath( product, variable );
  const int dataset = openDataset( file, path, variable );
  size_t result = 0;

  assert( file >= 0 );
  assert( product ); assert( *product );
  assert( variable ); assert( *variable );
  assert( dataMinimum >= -DBL_MAX );
  assert( dataMaximum <=  DBL_MAX );
  assert( dataMinimum <= dataMaximum );
  assert( rows );assert( columns );
  assert( data ); assert( temp );

  if ( dataset > -1 ) {
    result = readFileData( dataset, rows, columns, temp );
    closeDataset( dataset );
  }

  if ( result ) {
    const size_t points = rows * columns;
    size_t index = 0;

    for ( index = 0, result = 0; index < points; ++index ) {
      const double filterValue = temp[ index ];
      const int valid =
        filterValue >= dataMinimum && filterValue <= dataMaximum;

      if ( ! valid ) {
        data[ index ] = MISSING_VALUE;
      } else {
        result += data[ index ] > MISSING_VALUE;
      }
    }
  }

  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: variablePath - Get path to variable.
 INPUTS:  const char* const product  File type: OMNO2, etc.
         const char* const variable  Unpathed name of variable.
RETURNS: const char* Path to variable.
******************************************************************************/

static const char* variablePath( const char* const product,
                                 const char* const variable ) {
  const char* result = "";
  char spaceVariableSpace[ 80 ] = "";
  snprintf( spaceVariableSpace,
            sizeof spaceVariableSpace / sizeof *spaceVariableSpace,
            " %s ", variable );
  {
    const char* const coordinatesPath =
      ! strcmp( product, "OMNO2" ) ?
        "/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields"
      : ! strcmp( product, "OMTO3" ) ?
        "/HDFEOS/SWATHS/OMI Column Amount O3/Geolocation Fields"
      : ! strcmp( product, "OMCLDRR" ) ?
        "/HDFEOS/SWATHS/Cloud Product/Geolocation Fields"
      : "";
    const char* const dataPath =
       ! strcmp( product, "OMNO2" ) ?
         "/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields"
       : ! strcmp( product, "OMTO3" ) ?
        "/HDFEOS/SWATHS/OMI Column Amount O3/Data Fields"
       : ! strcmp( product, "OMCLDRR" ) ?
         "/HDFEOS/SWATHS/Cloud Product/Data Fields"
       : "";
    const char* const coordinatesPathVariables =
       ! strcmp( product, "OMNO2" ) ?
          " Latitude Longitude SolarAzimuthAngle SolarZenithAngle "
          " Time ViewingAzimuthAngle ViewingZenithAngle "
       : ! strcmp( product, "OMTO3" ) ?
          " GroundPixelQualityFlags LandWaterClassification Latitude Longitude "
          " RelativeAzimuthAngle SolarAzimuthAngle SolarZenithAngle "
          " TerrainHeight Time ViewingAzimuthAngle ViewingZenithAngle "
          " WaterFraction XTrackQualityFlags "
       : ! strcmp( product, "OMCLDRR" ) ?
          " GroundPixelQualityFlags Latitude Longitude "
          " RelativeAzimuthAngle "
          " SolarZenithAngle TerrainHeight Time "
          " ViewingZenithAngle XTrackQualityFlags "
       : "";
    result =
      strstr( coordinatesPathVariables, spaceVariableSpace ) ?
        coordinatesPath
      : dataPath;
  }

  return result;
}



/******************************************************************************
PURPOSE: openDataset - Open HDF5 dataset for reading.
INPUTS:  const int file                 ID of file to query.
         const char* const datasetPath  Path to variable to open.
         const char* const variable     Unpathed name of variable to open.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

static int openDataset( const int file, const char* const datasetPath,
                        const char* const variable ) {
  int result = -1;
  char pathedDatasetName[ 512 ] = "";

  assert( file > -1 );
  assert( datasetPath ); assert( *datasetPath == '/' );
  assert( variable ); assert( *variable );
  assert( ! strchr( variable, '/' ) ); /* variable name is unpathed. */

  memset( pathedDatasetName, 0, sizeof pathedDatasetName );
  snprintf( pathedDatasetName,
            sizeof pathedDatasetName / sizeof *pathedDatasetName,
            "%s/%s", datasetPath, variable );

  result = H5Dopen( file, pathedDatasetName, 0 );

  if ( result == -1 ) {
    fprintf( stderr, "\nFailed to open HDF5 dataset for reading: %s.\n",
             pathedDatasetName );
  }

  return result;
}



/******************************************************************************
PURPOSE: closeDataset - Close HDF5 dataset.
INPUTS:  const int dataset HDF5 dataset to close.
******************************************************************************/

static void closeDataset( const int dataset ) {
  H5Dclose( dataset );
}



/******************************************************************************
PURPOSE: readDatasetDimensions - Read dataset dimensions.
INPUTS:  const int dataset   ID of dataset to query.
OUTPUTS: size_t* rows     1st dimension of dataset.
         size_t* columns  2nd dimension of dataset or 0 if rank 1.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int readDatasetDimensions( const int dataset,
                                  size_t* const rows,
                                  size_t* const columns ) {
  int result = 0;

  assert( dataset > -1 ); assert( rows ); assert( columns );


  {
    const int dataspace = H5Dget_space( dataset );

    if ( dataspace > -1 ) {
      unsigned long long dims[ 32 ] = { 0, 0 };
      const int rank = H5Sget_simple_extent_dims( dataspace, dims, 0 );

      if ( rank == 1 && dims[ 0 ] > 0 ) {
        *rows = dims[ 0 ];
        *columns = 0;
        result = 1;
      } else if ( rank == 2 && dims[ 0 ] > 0 && dims[ 1 ] > 0 ) {
        *rows = dims[ 0 ];
        *columns = dims[ 1 ];
        result = 1;
      }

      H5Sclose( dataspace );
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read valid dimensions of dataset.\n" );
    *rows = *columns = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read, decode and expand file data.
INPUTS:  const int dataset              ID of dataset to read.
         const size_t rows              Rows data to read.
         const size_t columns           Columns of data to read.
OUTPUTS: double data[ rows * columns ]  Data read and decoded.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int readFileData( const int dataset,
                         const size_t rows,
                         const size_t columns,
                         double data[] ) {

  int result = 0;
  size_t dim0 = 0;
  size_t dim1 = 0;

  assert( dataset > -1 ); assert( rows ); assert( columns );
  assert( data );

  result = readDatasetDimensions( dataset, &dim0, &dim1 );

  if ( result ) {
    result = dim0 == rows && ( dim1 == columns || dim1 == 0 );

    if ( ! result ) {
      fprintf( stderr, "\nInvalid/mismatched dimensions in dataset: "
               "%lu %lu (expected %lu %lu)\n",
               dim0, dim1, rows, columns );
    }
  }

  if ( result ) {
    result = H5Dread( dataset, H5T_NATIVE_DOUBLE_g, 0,0,0, data ) >= 0;

    if ( ! result ) {
      fprintf( stderr, "\nFailed to read matched file data.\n" );
    } else { /* Decode values: */
      const double fillValue = readAttribute( dataset, "_FillValue" );
      const double scaleFactor = readAttribute( dataset, "ScaleFactor" );
      const size_t count = dim0 * ( dim1 > 0 ? dim1 : 1 );
      size_t index = 0;

      for ( index = 0; index < count; ++index ) {
        const double value = data[ index ];

        if ( fillValue != MISSING_VALUE && value == fillValue ) {
          data[ index ] = MISSING_VALUE;
        } else if ( scaleFactor != MISSING_VALUE ) {
          data[ index ] *= scaleFactor;
        }
      }

      if ( dim1 == 0 ) { /* Copy single-column row values to all columns: */
        const double* input = data + count;
        const size_t count2 = rows * columns;
        double* output = data + count2;

        while ( output != input ) {
          const double value = *--input;
          size_t counter = columns;

          while ( counter-- ) {
            *--output = value;
          }
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readAttribute - Read named dataset attribute.
INPUTS:  const int dataset       ID of dataset to read.
         const char* const name  Name of attribute.
RETURNS: double value read or MSSING_VALUE if named attribute is undefined.
******************************************************************************/

static double readAttribute( const int dataset, const char* const name ) {
  double result = MISSING_VALUE;
  assert( dataset >= 0 ); assert( name ); assert( *name );

  if ( H5Aexists( dataset, name ) > 0 ) {
    const int id = H5Aopen_name( dataset, name );

    if ( id >= 0 ) {
      const int status = H5Aread( id, H5T_NATIVE_DOUBLE_g, &result );

      if ( status < 0 ) {
        result = MISSING_VALUE;
      }

      H5Aclose( id );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: toSecondsUTC70 - Convert scan start time
         TAI seconds since 1993-01-01T00:00:00Z) to
         UTC seconds since January 1, 1970.
INPUTS:  const double secondsTAI93  Seconds TAI since January 1, 1993.
RETURNS: long long seconds UTC since January 1, 1970.
******************************************************************************/

static long long toSecondsUTC70( const double secondsTAI93 ) {

  /*
   * Compute approximate seconds from UNIX time: 1970-01-01T00:00:00Z
   * to TAI Scan Start Time base 1993-01-01T00:00:00Z
   * including 6 leap years: 1972, 1976, 1980, 1984, 1988, 1992 and
   * 17 leap seconds. See http://en.wikipedia.org/wiki/Leap_second
   */

  const int daysFrom1970To1993 = 8401; /* Includes 6 leap years. */
  const int leapSecondsFrom1970To1993 = 17;
  const int secondsDifferenceUTC2TAI = -10;
  const int offset = -13; /* HACK to match MODIS/CALIPSO UTC timestamp. */
  const int hoursPerDay = 24;
  const int minutesPerHour = 60;
  const int secondsPerMinute = 60;
  long long result =
    daysFrom1970To1993 * hoursPerDay * minutesPerHour * secondsPerMinute +
    leapSecondsFrom1970To1993 + secondsDifferenceUTC2TAI + offset;

  /*
   * Add TAI scan start time seconds after 1993-01-01T00:00:00Z
   * which includes leap seconds:
   */

  result += (long long) ( secondsTAI93 + 0.5 );
  return result;
}



/******************************************************************************
PURPOSE: toUTC - Convert UTC seconds since January 1, 1970 to UTC yyyydddhhmm.
INPUTS:  const long long seconds  Seconds since January 1, 1970.
RETURNS: long long yyyyddddhhmm UTC.
******************************************************************************/

static long long toUTC( const long long seconds ) {
  long long result = 19000010000LL;
  const time_t tSeconds = (const time_t) seconds;
  struct tm t;
  memset( &t, 0, sizeof t );

  if ( gmtime_r( &tSeconds, &t ) ) {
    const int yyyy = t.tm_year + 1900;
    const int ddd = t.tm_yday + 1;
    const int hh   = t.tm_hour;
    const int mm   = t.tm_min;
    result = yyyy;
    result *= 1000;
    result += ddd;
    result *= 100;
    result += hh;
    result *= 100;
    result += mm;
  }

  return result;
}



/******************************************************************************
PURPOSE: isValidTimestamp - Is the timestamp valid?
INPUTS:  const long long yyyydddhhmm The timestamp to examine.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidTimestamp( const long long yyyydddhhmm ) {
  const long long yyyy = yyyydddhhmm / 10000000;
  const long long ddd  = yyyydddhhmm / 10000 % 1000;
  const long long hh   = yyyydddhhmm / 100 % 100;
  const long long mm   = yyyydddhhmm % 100;
  const int isLeapYear = yyyy % 4 == 0 && (yyyy % 100 != 0 || yyyy % 400 == 0);
  const int result =
    yyyy >= 1900 && yyyy <= 9999 &&
    ddd >= 1 && ddd <= 365 + isLeapYear &&
    hh >= 0 && hh <= 23 &&
    mm >= 0 && mm <= 59;
  return result;
}


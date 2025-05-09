/******************************************************************************
PURPOSE: ReadFile.c - Simple to use wrapper routines to read data from
                      OMI-BEHR HDF5 files.

NOTES:   Uses HDF5 libraries and libraries they depend on (curl, z, dl).

HISTORY: 2017-12-06 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <stdlib.h> /* For malloc(), free(). */
#include <string.h> /* For strncpy(), strcmp(), memset(). */
#include <limits.h> /* For INT_MAX. */
#include <time.h>   /* For time(), gmtime_r(), timegm(). */

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

static int openDataset( const int file, const int swath,
                        const char* const variable );

static void closeDataset( const int dataset );

static int readDatasetDimensions( const int dataset,
                                  size_t* dimension0, size_t* dimension1 );

static int getDatasetName( const int file, const int swath,
                           const char* const variable,
                           char datasetName[ 256 ] );

static int readFileData( const int dataset,
                         const int isTime,
                         const size_t dimension0,
                         const size_t dimension1,
                         float data[] );

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
PURPOSE: swathsInFile - Query the number of swaths in file.
INPUTS:  const int file      ID of file to query.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int swathsInFile( const int file ) {
  int result = 0;
  int group = -1;

  assert( file > -1 );

  group = H5Gopen( file, "/Data", 0 );

  if ( group > -1 ) {
    unsigned long long count = 0;
    const int status = H5Gget_num_objs( group, &count );
    H5Gclose( group );
    group = -1;

    if ( status >= 0 && count > 0 && count < INT_MAX ) {
      result = count;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readDimensions - Read dimensions of swath in file.
INPUTS:  const int file    ID of file to query.
         const int swath   0-based index of swath within file.
OUTPUTS: size_t* rows      Number of rows in swath.
         size_t* columns   Number of columns in swath.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readDimensions( const int file, const int swath,
                    size_t* const rows, size_t* const columns ) {
  int result = 0;
  int dataset = -1;

  assert( file > -1 ); assert( swath >= 0 );
  assert( rows ); assert( columns );

  dataset = openDataset( file, swath, "Longitude" );

  if ( dataset > -1 ) {
    result = readDatasetDimensions( dataset, rows, columns );
    closeDataset( dataset );
    dataset = -1;
  }

  return result;
}



/******************************************************************************
PURPOSE: readDataset - Read appropriately-filtered variable data.
INPUTS:  const int file              ID of file to query.
         const int swath             0-based index of swath within file.
OUTPUTS: const size_t rows           Rows of dataset to match.
         const size_t columns        Columns of dataset to match.
         const char* const variable  Name of dataset variable.
         char units[ 80 ]            Units of variable.
         double data[ dimension0 * dimension1 ]  Filtered variable data values.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readDataset( const int file,
                 const int swath,
                 const size_t rows,
                 const size_t columns,
                 const char* const variable,
                 char units[ 80],
                 double data[] ) {

  typedef struct {
    const char* variable;        /* Name of data variable to filter. */
    const char* units;           /* Name of data units. */
    double dataMinimum;          /* Minimum valid value for this variable. */
    double dataMaximum;          /* Maximum valid value for this variable. */
    const char* filterVariable1; /* Name of 1st data variable to filter by. */
    const char* filterVariable2; /* Name of 2nd data variable to filter by. */
  } Entry;

  static const Entry table[] = {
    { "Longitude",                "deg", -180.0, 180.0, 0, 0 },
    { "Latitude",                 "deg", -90.0, 90.0, 0, 0 },
    { "Time",                     "YYYYDDDHHMM", 0.0, 1e20, 0, 0 },

    { "AMFStrat",                 "-",   0.0, 1e38, "XTrackQualityFlags", 0 },
    { "AMFTrop",                  "-",   0.0, 1e38, "XTrackQualityFlags", 0 },
    { "BEHRAMFTrop",              "-",   0.0, 1e38, "XTrackQualityFlags", 0 },
    { "BEHRAMFTropVisOnly",       "-",   0.0, 1e38, "XTrackQualityFlags", 0 },
    { "BEHRColumnAmountNO2Trop",  "molecules/cm2", 0.0, 1e38,
                                    "XTrackQualityFlags", "vcdQualityFlags" },
    { "BEHRColumnAmountNO2TropVisOnly",  "molecules/cm2", 0.0, 1e38,
                                    "XTrackQualityFlags", "vcdQualityFlags" },
    { "CloudFraction",            "-",   0.0, 1.0,  "XTrackQualityFlags", 0 },
    { "CloudPressure",            "hPa", 0.0, 2e4,  "XTrackQualityFlags", 0 },
    { "CloudRadianceFraction",    "-",   0.0, 1.0,  "XTrackQualityFlags", 0 },
    { "ColumnAmountNO2",          "molecules/cm2", 0.0, 1e38,
                                    "XTrackQualityFlags", "vcdQualityFlags" },
    { "ColumnAmountNO2Strat",     "molecules/cm2", 0.0, 1e38,
                                    "XTrackQualityFlags", "vcdQualityFlags" },
    { "ColumnAmountNO2Trop",      "molecules/cm2", 0.0, 1e38,
                                    "XTrackQualityFlags", "vcdQualityFlags" },
    { "ColumnAmountNO2TropStd",   "molecules/cm2", 0.0, 1e38,
                                    "XTrackQualityFlags", "vcdQualityFlags" },
    { "GLOBETerpres",             "hPa", 0.0, 2e4,    0, 0 },
    { "MODISAlbedo",              "-",   0.0, 1.0,    0, 0 },
    { "MODISCloud",               "-",   0.0, 1.0,    0, 0 },
    { "RelativeAzimuthAngle",     "deg", 0.0, 180.0,  0, 0 },
    { "Row",                      "-",   0.0, 1e5,    0, 0 },
    { "SlantColumnAmountNO2",     "molecules/cm2", 0.0, 1e38,
                                    "XTrackQualityFlags", 0 },
    { "SolarAzimuthAngle",        "deg", -180.0, 180.0,  0, 0 },
    { "SolarZenithAngle",         "deg",    0.0,  90.0,  0, 0 },
    { "Swath",                    "-",      0.0,  1e10,  0, 0 },
    { "TerrainHeight",            "m",   -500.0,  1e4,   0, 0 },
    { "TerrainPressure",          "hPa",   0.0,   2e4,   0, 0 },
    { "TerrainReflectivity",      "-",     0.0,   1.0,   0, 0 },
    { "ViewingAzimuthAngle",      "deg", -180.0, 180.0,  0, 0 },
    { "ViewingZenithAngle",       "deg",    0.0,  90.0,  0, 0 }
  };
  const size_t entries = sizeof table / sizeof *table;
  int result = 0;
  size_t index = 0;

  assert( file > -1 ); assert( swath >= 0 );
  assert( rows ); assert( columns );
  assert( variable ); assert( *variable );
  assert( units ); assert( data );
  assert( sizeof (double) == 2 * sizeof (float) );

  units[ 0 ] = '\0';

  /* Find variable in table. If not found it is an error: */

  while ( index < entries && strcmp( table[ index ].variable, variable ) ) {
    ++index;
  }

  if ( index == entries ) {
    fprintf( stderr, "\a\nInvalid variable '%s'.\n", variable );
  } else {
    const Entry* const entry = table + index;
    int dataset = openDataset( file, swath, variable );

    if ( dataset > -1 ) {
      float* fdata = (float*) data; /* Read file floats into 1st half of data*/
      const int isTime = ! strcmp( variable, "Time" );
      int ok = readFileData( dataset, isTime, rows, columns, fdata );
      closeDataset( dataset ), dataset = -1;

      if ( ok ) {
        const size_t points = rows * columns;
        const char* const filterVariable1 = entry->filterVariable1;
        const char* const filterVariable2 = entry->filterVariable2;
        strncpy( units, entry->units, 79 );

        if ( isTime ) { /* Convert from seconds TAI93 to yyyymmddhhmmss UTC: */
          double secondsTAI93 = fdata[ 0 ];
          const long long secondsUTC70 = toSecondsUTC70( secondsTAI93 );
          const long long yyyydddhhmm = toUTC( secondsUTC70 );
          result = isValidTimestamp( yyyydddhhmm );

          if ( result ) {

            /*
             * Since conversion is expensive and
             * entire swath scans complete in 512 seconds and
             * we only need hourly accuracy,
             * just replicate the first value converted:
             */

            for ( index = 0; index < points; ++index ) {
              data[ index ] = yyyydddhhmm;
            }
          }
        } else if ( filterVariable1 ) { /* Read and apply filterVariable1: */
          float* const filterData = fdata + points; /* data 2nd half is flags*/
          dataset = openDataset( file, swath, filterVariable1 );

          if ( dataset > -1 ) {
            ok = readFileData( dataset, 0, rows, columns, filterData );
            closeDataset( dataset ), dataset = -1;
          }

          if ( ok ) { /* Filter fdata by filterVariable1's filterData: */
            ok = 0;

            for ( index = 0; index < points; ++index ) {
              const float filterValue = filterData[ index ];
              const int valid = filterValue == 0.0;

              if ( ! valid ) {
                fdata[ index ] = MISSING_VALUE;
              } else {
                ok = 1;
              }
            }
          }

          if ( ok && filterVariable2 ) { /* Read and apply filterVariable2: */
            dataset = openDataset( file, swath, filterVariable2 );

            if ( dataset > -1 ) {
              ok = readFileData( dataset, 0, rows, columns, filterData );
              closeDataset( dataset ), dataset = -1;
            }

            if ( ok ) { /* Filter fdata by filterVariable2's filterData: */
              ok = 0;

              for ( index = 0; index < points; ++index ) {
                const float filterValue = filterData[ index ];
                const int valid = filterValue == 0.0;

                if ( ! valid ) {
                  fdata[ index ] = MISSING_VALUE;
                } else {
                  ok = 1;
                }
              }
            }
          }
        }

        /* Filter fdata by valid range and expand to double data: */

        if ( ok && ! isTime ) {
          const double dataMinimum = entry->dataMinimum;
          const double dataMaximum = entry->dataMaximum;
          index = points; /* Must loop backward to avoid overwrite. */

          while ( index-- ) {
            const float fvalue = fdata[ index ];
            const int valid = fvalue >= dataMinimum && fvalue <= dataMaximum;

            if ( valid ) {
              data[ index ] = fvalue; /* Expand float to double. */
              result = 1; /* At least 1 valid data value was read. */
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



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: openDataset - Open HDF5 dataset for reading.
INPUTS:  const int file               ID of file to query.
         const int swath              0-based index of swath within file.
         const char* const variable   Unpathed name of variable to open.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

static int openDataset( const int file, const int swath,
                        const char* const variable ) {
  int result = -1;
  char pathedDatasetName[ 256 ] = "";

  assert( file > -1 ); assert( swath >= 0 );
  assert( variable ); assert( *variable );
  assert( ! strchr( variable, '/' ) ); /* variable name is unpathed. */

  memset( pathedDatasetName, 0, sizeof pathedDatasetName );

  if ( getDatasetName( file, swath, variable, pathedDatasetName ) ) {
    result = H5Dopen( file, pathedDatasetName, 0 );

    if ( result == -1 ) {
      fprintf( stderr, "\a\nFailed to open HDF5 dataset for reading: %s.\n",
               pathedDatasetName );
    }
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
OUTPUTS: size_t* dimension0  1st dimension of dataset.
         size_t* dimension1  2nd dimension of dataset or 0 if rank 1.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int readDatasetDimensions( const int dataset,
                                  size_t* const dimension0,
                                  size_t* const dimension1 ) {
  int result = 0;

  assert( dataset > -1 ); assert( dimension0 ); assert( dimension1 );

  {
    const int dataspace = H5Dget_space( dataset );

    if ( dataspace > -1 ) {
      unsigned long long dims[ 32 ] = { 0, 0 };
      const int rank = H5Sget_simple_extent_dims( dataspace, dims, 0 );

      if ( rank == 1 && dims[ 0 ] > 0 ) {
        *dimension0 = dims[ 0 ];
        *dimension1 = 0;
        result = 1;
      } else if ( rank == 2 && dims[ 0 ] > 0 && dims[ 1 ] > 0 ) {
        *dimension0 = dims[ 0 ];
        *dimension1 = dims[ 1 ];
        result = 1;
      }

      H5Sclose( dataspace );
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read valid dimensions of dataset.\n" );
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file data.
INPUTS:  const int dataset              ID of dataset to read.
         const int isTime               Is dataset Time variable?
         const size_t dimension0        1st dimension of data to read.
         const size_t dimension1        2nd dimension of data to read or 0.
OUTPUTS: float data[ dimension0 * dimension1 ]  Data read.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int readFileData( const int dataset,
                         const int isTime,
                         const size_t dimension0,
                         const size_t dimension1,
                         float data[] ) {

  int result = 0;
  size_t dim0 = 0;
  size_t dim1 = 0;

  assert( dataset > -1 ); assert( dimension0 );
  assert( data );

  result = readDatasetDimensions( dataset, &dim0, &dim1 );

  if ( result ) {
    result =
      ( dim0 == dimension0 || ( isTime == 1 && dim0 == 1 ) ) &&
      dim1 == dimension1;

    if ( ! result ) {
      fprintf( stderr, "\a\nInvalid/mismatched dimensions in dataset: "
               "%lu %lu (expected %lu %lu)\n",
               dim0, dim1, dimension0, dimension1 );
    }
  }

  if ( result ) {
    result = H5Dread( dataset, H5T_NATIVE_FLOAT_g, 0,0,0, data ) >= 0;

    if ( ! result ) {
      fprintf( stderr, "\a\nFailed to read matched file data.\n" );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: getDatasetName - Get fully pathed name of index swath dataset.
INPUTS:  const int file              ID of file.
         const int swath             0-based index of swath within file.
         const char* const variable  Name of variable/dataset.
OUTPUTS: char datasetName[ 256 ]     Pathed name of indexed variable dataset.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int getDatasetName( const int file, const int swath,
                           const char* const variable,
                           char datasetName[ 256 ] ) {
  int result = 0;
  int group = -1;
  char temp[ 256 ] = "";

  assert( file > -1 ); assert( swath >= 0 ); assert( datasetName );

  memset( temp, 0, sizeof temp );
  memset( datasetName, 0, 256 * sizeof (char) );
  group = H5Gopen( file, "/Data", 0 );

  if ( group >= 0 ) {
    const int status = H5Gget_objname_by_idx( group, swath, temp, 128 );
  
    if ( status > 0 && *temp ) {
      strcpy( datasetName, "/Data/" );
      strcat( datasetName, temp );
      strcat( datasetName, "/" );
      strcat( datasetName, variable );
      assert( datasetName[ 255 ] == '\0' );
      result = 1;
    }

    H5Gclose( group );
    group = -1;
  }

  /*
  fprintf( stderr,
           "getDatasetName( file = %d, swath = %d, variable = %s ) = %s\n",
           file, swath, variable, datasetName );
  */

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
  const long long result =
    yyyy >= 1900 && yyyy <= 9999 &&
    ddd >= 1 && ddd <= 365 + isLeapYear &&
    hh >= 0 && hh <= 23 &&
    mm >= 0 && mm <= 59;
  return result;
}


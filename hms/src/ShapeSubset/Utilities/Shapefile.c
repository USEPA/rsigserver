
/******************************************************************************
PURPOSE: Shapefile.c - Write time-varying scalar grid and point data to ESRI
         Shapefiles (shp, shx, dbf) and ASCIIGrid files (asc, prj).
NOTES:   Uses Shapefile and GPC free, open-source libraries.
         http://shapelib.maptools.org/shp_api.html
         http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
         http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
         http://en.wikipedia.org/wiki/ESRI_grid
         http://en.wikipedia.org/wiki/Well-known_text
         When adding new shapefile files, edit table in writePolygonDBF().
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For FILE, stderr, fprintf(), fopen(), perror(),fputs()*/
#include <string.h> /* For strncpy(), memset(), strdup(). */
#include <ctype.h>  /* For tolower(), isprint(), isspace(). */
#include <limits.h> /* For INT_MAX. */
#include <stdlib.h> /* For malloc(), free(). */
#include <unistd.h> /* For unlink(). */

#include <Assertions.h>    /* For PRE*(), POST*(), DEBUG(), IMPLIES_ELSE(). */
#include <BasicNumerics.h> /* For CLAMPED_TO_RANGE(). */
#include <DateTime.h>      /* For isValidYearMonthDay(). */
#include <Utilities.h>     /* For daysInMonth(), copyFile(), streamFile(). */
#include <projections.h>   /* For is_valid_longitude_latitude(). */
#include <albers.h>        /* For initialize_albers(), project_albers(). */
#include <Shapefile.h>     /* For public interface. */

/* Memory allocation/deallocation: */

#define NEW( type, count ) (type*) allocate( sizeof (type), count )
#define FREE( p ) ( ( (p) ? free(p) : (void) 0 ), (p) = 0 )

/* The NEW() macro above calls allocate() which zeros the returned memory: */

static void* allocate( size_t bytesEach, size_t count );


static long long binarySearchAny( const long long value,
                                  const long long count,
                                  const long long values[] ) {
  long long result = -1;

  if ( count > 0 ) {
    long long left = 0;
    long long right = count;
    long long middle = 0;
    long long middleValue = 0;

    while ( left < right ) {
      middle = left + ( right - left ) / 2;
      middleValue = values[ middle ];

      if ( middleValue < value ) {
        left = middle + 1;
      } else if ( middleValue > value ) {
        right = middle;
      } else {
        left = middle;
        right = left;
      }
    }

    if ( middleValue == value ) {
      result = middle;
    }
  }

  return result;
}



/*================================== TYPES ==================================*/


static const int includeMissingValuesInCSVFile = 0; /* Write -9999.0 values? */

enum { BIG = 4321, LITTLE = 1234 };
enum { MAXIMUM_FILE_NAME_LENGTH = 255 };
enum { MAXIMUM_CSV_HEADER_LINE_LENGTH = 1023 };
typedef char CSVHeader[ MAXIMUM_CSV_HEADER_LINE_LENGTH + 1 ];

/*
 * The following table defines the translation of input DBF columns to
 * output DBF columns.
 * This is used by the ShapeSubset program when it reads PI-supplied DBF files
 * (stored on rtpmeta.epa.gov) and, after spatial subsetting, creates
 * subsetted SHP and DBF files for streaming back to users.
 */

#define ACRES_TO_HECTARES 0.404685642
#define FEET_TO_METERS 0.3048
#define CUBIC_FEET_TO_CUBIC_METERS \
  ( FEET_TO_METERS * FEET_TO_METERS * FEET_TO_METERS )
#define MILES_TO_KM 1.609344
#define GALLONS_TO_LITERS 3.785412
#define LITERS_PER_CUBIC_METER 0.001

/* Define columns of output DBF file: */

typedef struct {
  const char* fileName;   /* Unique part of input file name "temperature_". */
  int inputColumn;        /* 0-based column number of original input DBF. */
  const char* columnName; /* E.g., "TEMP_C". */
  int columnType;         /* FTString, FTInteger, FTDouble. */
  int fieldWidth;         /* E.g., 7 characters wide. */
  int decimals;           /* E.g., 3 digits to right of the decimal. */
  double offset;          /* Offset to map input to output. -32.0. */
  double scale;           /* Scale  to map input to output. 5.0/9.0. */
} ColumnEntry;

static const ColumnEntry table[] = {

  /* States: */

  { "states_",  0, "STATE_FIPS", FTInteger, 2, 0, 0.0, 1.0 },
  { "states_",  1, "STATE",      FTString,  2, 0, 0.0, 1.0 },
  { "states_",  2, "STATE_NAME", FTString, 24, 0, 0.0, 1.0 },
  { "states_",  3, "EPA_REGION", FTInteger, 2, 0, 0.0, 1.0 },
  { "states_",  4, "EPA_GEOID",  FTString,  2, 0, 0.0, 1.0 },
  { "states_",  5, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "states_",  6, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "provinces_",  0, "PRUID",     FTInteger, 3, 0, 0.0, 1.0 },
  { "provinces_",  1, "NAME",      FTString, 24, 0, 0.0, 1.0 },
  { "provinces_",  2, "PR_ABBR",   FTString,  8, 0, 0.0, 1.0 },
  { "provinces_",  3, "AREA_SQKM", FTDouble, 11, 3, 0.0, 1.0 },
  { "provinces_",  4, "HECTARES",  FTDouble, 20, 6, 0.0, 1.0 },

  /* Counties: */

  { "counties_",  0, "STATE_FIPS", FTInteger, 2, 0, 0.0, 1.0 },
  { "counties_",  1, "STATE",      FTString,  2, 0, 0.0, 1.0 },
  { "counties_",  2, "COUNTYNAME", FTString, 24, 0, 0.0, 1.0 },
  { "counties_",  3, "COUNTYFIPS", FTInteger, 4, 0, 0.0, 1.0 },
  { "counties_",  4, "EPA_GEOID",  FTString,  5, 0, 0.0, 1.0 },
  { "counties_",  5, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "counties_",  6, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Cities: */

  { "cities_",  0, "STATE_FIPS", FTInteger, 2, 0, 0.0, 1.0 },
  { "cities_",  1, "STATE",      FTString,  2, 0, 0.0, 1.0 },
  { "cities_",  2, "CITY_NAME",  FTString, 64, 0, 0.0, 1.0 },
  { "cities_",  3, "GEOID10",    FTInteger, 5, 0, 0.0, 1.0 },
  { "cities_",  4, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "cities_",  5, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Roads: */

  { "roads_",  0, "LINEARID",   FTString,  16, 0, 0.0, 1.0 },
  { "roads_",  1, "NAME",       FTString,  64, 0, 0.0, 1.0 },
  { "roads_",  2, "HWY_TYPE",   FTString,   1, 0, 0.0, 1.0 },

  /* Tributaries: */

  { "tributaries_great_lakes",  6, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  7, "FDATE",      FTString,  10, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  8, "GNIS_ID",    FTString,   8, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  9, "GNIS_NAME",  FTString,  40, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  4, "LENGTH_KM",  FTDouble,   8, 3, 0.0, 1e-3 },
  { "tributaries_great_lakes", 10, "REACH_CODE", FTString,  14, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 11, "WBAREACOMI", FTInteger,  9, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  0, "CD_LINK",    FTString,   8, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  1, "ORIGUNIT",   FTString,   2, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  2, "LAKE",       FTString,   2, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  3, "COUNTRY",    FTString,  10, 0, 0.0, 1.0 },
  { "tributaries_great_lakes",  5, "GLHDID",     FTInteger,  5, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 12, "GNIS_NBR",   FTInteger,  1, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 13, "TERTIARY",   FTString,   3, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 14, "CAN_NAME",   FTString,  32, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 15, "SHREVEDL",   FTInteger,  4, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 16, "STRAHLERDL", FTInteger,  4, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 17, "STATE_FIPS", FTInteger,  2, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 18, "STATE",      FTString,   2, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 19, "STATE_NAME", FTInteger, 16, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 20, "EPA_REGION", FTInteger,  2, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 21, "PRUID",      FTInteger,  2, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 22, "NAME",       FTString,   8, 0, 0.0, 1.0 },
  { "tributaries_great_lakes", 23, "PR_ABBR",    FTString,   4, 0, 0.0, 1.0 },

  /* All other tributaries NOT great_lakes: */

  { "tributaries_!great_lakes",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes",  1, "FDATE",      FTString,  10, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes",  2, "RESOLUTION", FTString,   8, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes",  3, "GNIS_ID",    FTString,   8, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes",  4, "GNIS_NAME",  FTString,  40, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes",  5, "LENGTH_KM",  FTDouble,   8, 3, 0.0, 1e-3},
  { "tributaries_!great_lakes",  6, "REACH_CODE", FTString,  14, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes",  8, "WBAREACOMI", FTInteger, 10, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes",  9, "FTYPE",      FTString,  16, 0, 0.0, 1.0 },
  { "tributaries_!great_lakes", 10, "FCODE",      FTInteger,  5, 0, 0.0, 1.0 },

  /* HUCs data: */

  { "hucs",  0, "HUC",        FTString,  8, 0, 0.0, 1.0 },
  { "hucs",  1, "NAME",       FTString, 48, 0, 0.0, 1.0 },
  { "hucs",  2, "STATES",     FTString, 24, 0, 0.0, 1.0 },
  { "hucs",  3, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "hucs",  4, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Watersheds data: */

  { "watersheds_great_lakes",  0, "GLHDID",     FTInteger,  5, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  8, "WATER_BODY", FTString,  48, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  1, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  2, "FDATE",      FTString,  10, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  3, "GNIS_ID",    FTString,   8, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  4, "GNIS_NAME",  FTString,  40, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  5, "REACH_CODE", FTString,  14, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  6, "WBAREACOMI", FTInteger, 10, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  7, "TERTIARY",   FTString,   3, 0, 0.0, 1.0 },
  { "watersheds_great_lakes",  9, "ORIGUNIT",   FTString,   2, 0, 0.0, 1.0 },
  { "watersheds_great_lakes", 10, "LAKE",       FTString,   2, 0, 0.0, 1.0 },
  { "watersheds_great_lakes", 11, "COUNTRY",    FTString,  10, 0, 0.0, 1.0 },
  { "watersheds_great_lakes", 12, "SHREVEDL",   FTInteger,  4, 0, 0.0, 1.0 },
  { "watersheds_great_lakes", 13, "STRAHLERDL", FTInteger,  4, 0, 0.0, 1.0 },
  { "watersheds_great_lakes", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "watersheds_great_lakes", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* All other watersheds NOT great_lakes: */

  { "watersheds_!great_lakes",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  2, "STATE",      FTString,  2, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  3, "STATE2",     FTString,  2, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  4, "STATE3",     FTString,  2, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  5, "STATE4",     FTString,  2, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  6, "STATE5",     FTString,  2, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  7, "STATE6",     FTString,  2, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes",  8, "STATE7",     FTString,  2, 0, 0.0, 1.0 },
  { "watersheds_!great_lakes", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "watersheds_!great_lakes", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Estuaries data: */

  { "estuaries",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "estuaries",  1, "STATE",      FTString,  2, 0, 0.0, 1.0 },
  { "estuaries",  2, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "estuaries",  3, "NCA_NAME",   FTString, 40, 0, 0.0, 1.0 },
  { "estuaries",  4, "SUBCODE",    FTString, 48, 0, 0.0, 1.0 },
  { "estuaries",  5, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },
  { "estuaries",  6, "PROVINCE",   FTString, 48, 0, 0.0, 1.0 },
  { "estuaries",  7, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "estuaries",  8, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Lakes data: */

  { "/lakes_",  1, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "/lakes_", 16, "STATE",      FTString,   2, 0, 0.0, 1.0 },
  { "/lakes_",  5, "GNIS_NAME",  FTString,  48, 0, 0.0, 1.0 },
  { "/lakes_",  4, "GNIS_ID",    FTString,   8, 0, 0.0, 1.0 },
  { "/lakes_",  8, "REACHCODE",  FTString,  16, 0, 0.0, 1.0 },
  { "/lakes_",  2, "FDATE",      FTString,  16, 0, 0.0, 1.0 },
  { "/lakes_",  7, "ELEVATIONM", FTDouble,   8, 1, 0.0, 1.0 },
  { "/lakes_",  3, "RESOLUTION", FTString,   8, 1, 0.0, 1.0 },
  { "/lakes_",  9, "FTYPE",      FTString,  16, 1, 0.0, 1.0 },
  { "/lakes_", 10, "FCODE",      FTInteger,  5, 0, 0.0, 1.0 },
  { "/lakes_", 11, "SHAPE_LENG", FTDouble,  20, 4, 0.0, 1.0 },
  { "/lakes_", 12, "SHAPE_AREA", FTDouble,  20, 4, 0.0, 1.0 },
  { "/lakes-", 13, "SHORE_DIST", FTDouble,  20, 4, 0.0, 1.0 },
  { "/lakes_", 14, "MAX_WINDOW", FTDouble,  20, 4, 0.0, 1.0 },
  { "/lakes-", 15, "AREA",       FTDouble,  20, 4, 0.0, 1.0 },
  { "/lakes_",  6, "AREA_SQKM",  FTDouble,  11, 3, 0.0, 1.0 },
  { "/lakes_", -1, "HECTARES",   FTDouble,  20, 6, 0.0, 1.0 },

  /* Coastal Zone Management Areas: */

  { "coastal_zone_man",  0, "STATE_FIPS", FTInteger, 2, 0, 0.0, 1.0 },
  { "coastal_zone_man",  1, "COUNTYNAME", FTString, 24, 0, 0.0, 1.0 },
  { "coastal_zone_man",  2, "COUNTYFIPS", FTInteger, 4, 0, 0.0, 1.0 },
  { "coastal_zone_man",  3, "COUNTYNS",   FTInteger, 8, 0, 0.0, 1.0 },
  { "coastal_zone_man",  4, "CNTYIDFP",   FTInteger, 5, 0, 0.0, 1.0 },
  { "coastal_zone_man",  5, "CLASSFP",    FTString,  2, 0, 0.0, 1.0 },
  { "coastal_zone_man",  6, "CBSAFP",     FTInteger, 5, 0, 0.0, 1.0 },
  { "coastal_zone_man", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "coastal_zone_man", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Legislative Districts: */

  { "legislative_districts_",  0, "STATE_FIPS", FTInteger, 2, 0, 0.0, 1.0},
  { "legislative_districts_",  1, "DISTR_FIPS", FTInteger, 4, 0, 0.0, 1.0},
  { "legislative_districts_",  2, "GEOID",      FTInteger, 4, 0, 0.0, 1.0},
  { "legislative_districts_",  3, "LSA",        FTString,  2, 0, 0.0, 1.0},
  { "legislative_districts_",  4, "CDSESSN",    FTInteger, 3, 0, 0.0, 1.0},
  { "legislative_districts_",  5, "MTFCC",      FTString,  5, 0, 0.0, 1.0},
  { "legislative_districts_",  6, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0},
  { "legislative_districts_",  7, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0},

#if 1

  /* Original chlorophyll dbf files before splitting into domains 2022-08-22:*/

  { "nventoryunique", 8, "SITE_ID",    FTInteger, 7,  0, 0.0, 1.0 },
  { "nventoryunique", 3, "ESTCODE",    FTString,  5,  0, 0.0, 1.0 },
  { "nventoryunique", 1, "STATE",      FTString,  2,  0, 0.0, 1.0 },
  { "nventoryunique", 6, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "nventoryunique", 7, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },
  { "nventoryunique", 0, "PROJECT",    FTString, 24,  0, 0.0, 1.0 },
  { "nventoryunique", 2, "ESTUARY",    FTString, 48,  0, 0.0, 1.0 },
  { "nventoryunique", 4, "STATION",    FTString, 64,  0, 0.0, 1.0 },

  /* Original chlorophyll dbf files before splitting into domains 2022-08-22:*/

  { "entinel2matched", 11, "SITE_ID",    FTInteger, 7,  0, 0.0, 1.0 },
  { "entinel2matched",  6, "ESTCODE",    FTString,  5,  0, 0.0, 1.0 },
  { "entinel2matched",  4, "STATE",      FTString,  2,  0, 0.0, 1.0 },
  { "entinel2matched",  1, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "entinel2matched",  0, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },
  { "entinel2matched",  2, "PROJECT",    FTString, 24,  0, 0.0, 1.0 },
  { "entinel2matched",  5, "ESTUARY",    FTString, 48,  0, 0.0, 1.0 },
  { "entinel2matched",  3, "STATION",    FTString, 32,  0, 0.0, 1.0 },
  { "entinel2matched", 10, "STATIONS",   FTInteger, 2,  0, 0.0, 1.0 },

#endif

  /* Processed chlorophyll dbf files after splitting into domains 2022-08-22:*/

  { "esat1_water_quality_", 0, "SITE_ID",    FTInteger, 7,  0, 0.0, 1.0 },
  { "esat1_water_quality_", 1, "ESTCODE",    FTString,  5,  0, 0.0, 1.0 },
  { "esat1_water_quality_", 2, "STATE",      FTString,  2,  0, 0.0, 1.0 },
  { "esat1_water_quality_", 3, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "esat1_water_quality_", 4, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },
  { "esat1_water_quality_", 5, "PROJECT",    FTString, 24,  0, 0.0, 1.0 },
  { "esat1_water_quality_", 6, "ESTUARY",    FTString, 48,  0, 0.0, 1.0 },
  { "esat1_water_quality_", 7, "STATION",    FTString, 48,  0, 0.0, 1.0 },

  /* Processed chlorophyll dbf files after splitting into domains 2022-08-22:*/

  { "esat2_water_quality_", 0, "SITE_ID",    FTInteger, 7,  0, 0.0, 1.0 },
  { "esat2_water_quality_", 1, "ESTCODE",    FTString,  5,  0, 0.0, 1.0 },
  { "esat2_water_quality_", 2, "STATE",      FTString,  2,  0, 0.0, 1.0 },
  { "esat2_water_quality_", 3, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "esat2_water_quality_", 4, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },
  { "esat2_water_quality_", 5, "PROJECT",    FTString, 24,  0, 0.0, 1.0 },
  { "esat2_water_quality_", 6, "ESTUARY",    FTString, 48,  0, 0.0, 1.0 },
  { "esat2_water_quality_", 7, "STATION",    FTString, 48,  0, 0.0, 1.0 },
  { "esat2_water_quality_", 8, "STATIONS",   FTInteger, 2,  0, 0.0, 1.0 },

  /* Soil data: */

  { "soil_",  2, "MUID",       FTString, 7, 0, 0.0, 1.0 },
  { "soil_",  3, "STATE",      FTString, 2, 0, 0.0, 1.0 },
  { "soil_",  4, "AVWATERCAP", FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_",  5, "CLAY_%",     FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_",  6, "KFFACT",     FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_",  7, "ORGANIC_%",  FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_",  8, "PERM_mmhr",  FTDouble, 8, 3, 0.0, 25.4 },
  { "soil_",  9, "THICK_mm",   FTDouble, 8, 3, 0.0, 25.4 },
  { "soil_", 10, "HYGRP",      FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_", 11, "DRAINAGE",   FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_", 12, "SLOPE",      FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_", 13, "LIQUID_LIM", FTDouble, 7, 3, 0.0, 1.0 },
  { "soil_", 14, "HYDRIC_%",   FTDouble, 7, 3, 0.0, 100.0 },
  { "soil_", 15, "ANN_FLOOD",  FTDouble, 7, 3, 0.0, 1.0 },
  
  /* Wetlands data: */

  { "wetlands_",  1, "CODE",       FTString, 16, 0, 0.0, 1.0 },
  { "wetlands_",  4, "WETLAND_TY", FTString, 40, 0, 0.0, 1.0 },
  { "wetlands_", -1, "ACRES",      FTDouble, 20, 5, 0.0, 1.0 },
  { "wetlands_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Seagrass point data pacific (updated 2023-11-01): */

  { "seagrass_point_pacific", 21, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 22, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 23, "SUBCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 24, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_pacific",  1, "ESTUARY_NA", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 19, "CUR_SOURCE", FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 11, "CUR_DATE",   FTString, 24, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 14, "CURRENT_AC",  FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_pacific", 14, "CURRENT_HA", FTDouble, 16, 6, 0.0,
                                                   ACRES_TO_HECTARES },
  { "seagrass_point_pacific", 18, "MAX_OBS_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_pacific", 18, "MAX_OBS_HA", FTDouble, 16, 6, 0.0,
                                                   ACRES_TO_HECTARES },
  { "seagrass_point_pacific", 15, "ZOSTERA_M",  FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 16, "ZOSTERA_J",  FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 17, "ZOSTERA_P",  FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_pacific",  7, "EELGSTATUS", FTString, 32, 0, 0.0, 1.0 },
  { "seagrass_point_pacific",  8, "EELGSOURCE", FTString, 180, 0, 0.0, 1.0 },
  { "seagrass_point_pacific",  3, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point_pacific",  4, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point_pacific",  0, "OBJECTID",   FTInteger, 4,  0, 0.0, 1.0 },
  { "seagrass_point_pacific",  2, "SYSTEM_ORD", FTInteger, 4,  0, 0.0, 1.0 },
  { "seagrass_point_pacific",  5, "LINK",       FTString, 48,  0, 0.0, 1.0 },
  { "seagrass_point_pacific",  6, "PMEPREGION", FTString, 48,  0, 0.0, 1.0 },
  { "seagrass_point_pacific", 20, "PMEP_CODE",  FTInteger, 4,  0, 0.0, 1.0 },
  { "seagrass_point_pacific",  9, "HABITAT_CH", FTString, 80,  0, 0.0, 1.0 },
  { "seagrass_point_pacific", 10, "NOTES",      FTString, 160, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 12, "OTHER_YEAR", FTString, 128, 0, 0.0, 1.0 },
  { "seagrass_point_pacific", 13, "DATA_COUNT", FTInteger, 4, 0, 0.0, 1.0 },

  /* Seagrass point data atlantic (updated 2024-04-18): */

  { "seagrass_point_atlantic",  4, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic",  3, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic",  9, "CUR_SOURCE", FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 20, "SOURCE2",    FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 21, "SOURCE3",    FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 12, "CUR_YEAR",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 15, "CURRENT_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_atlantic", 15, "CURRENT_HA", FTDouble, 16, 6, 0.0,
                                                   ACRES_TO_HECTARES },
  { "seagrass_point_atlantic", 19, "MAX_OBS_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_atlantic", 19, "MAX_OBS_HA", FTDouble, 16, 6, 0.0,
                                                   ACRES_TO_HECTARES },
  { "seagrass_point_atlantic", 14, "DATA_COUNT", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic",  8, "STATUS",     FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 13, "PREV_YEARS", FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic",  6, "LINK",       FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 27, "OTHERLINKS", FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 23, "LINK2",      FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 24, "SOURCE2_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_atlantic", 25, "LINK3",      FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 26, "SOURCE3_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_atlantic", 45, "EELGRASSUI", FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic",  2, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point_atlantic",  1, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point_atlantic", 31, "ECO_CODE",   FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 32, "PROV_CODE",  FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 33, "PROVINCE",   FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 34, "REALM_CODE", FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 35, "REALM",      FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 38, "LAT_ZONE",   FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 39, "FID",        FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 36, "ALT_CODE",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 37, "ECO_CODE_X", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 29, "SUBESTUARY", FTString, 32, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 46, "WATER_BODY", FTString, 80, 0, 0.0, 1.0 },
  { "seagrass_point_atlantic", 47, "REGION",     FTString, 40, 0, 0.0, 1.0 },

  /* Seagrass point data gulf (2024-04-18): */

  { "seagrass_point_gulf",  3, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point_gulf",  2, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_gulf",  8, "CUR_SOURCE", FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 19, "SOURCE2",    FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 20, "SOURCE3",    FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 11, "CUR_YEAR",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 14, "CURRENT_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_gulf", 14, "CURRENT_HA", FTDouble, 16, 6, 0.0,
                                                   ACRES_TO_HECTARES },
  { "seagrass_point_gulf", 18, "MAX_OBS_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_gulf", 18, "MAX_OBS_HA", FTDouble, 16, 6, 0.0,
                                                   ACRES_TO_HECTARES },
  { "seagrass_point_gulf", 13, "DATA_COUNT", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_gulf",  7, "STATUS",     FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 12, "PREV_YEARS", FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_gulf",  5, "LINK",       FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 26, "OTHERLINKS", FTString, 256, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 22, "LINK2",      FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 23, "SOURCE2_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_gulf", 24, "LINK3",      FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 25, "SOURCE3_AC", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point_gulf", 44, "EELGRASSUI", FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_gulf",  1, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point_gulf",  0, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point_gulf", 31, "ECO_CODE",   FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 32, "PROV_CODE",  FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 33, "PROVINCE",   FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 34, "REALM_CODE", FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 35, "REALM",      FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 38, "LAT_ZONE",   FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 39, "FID",        FTInteger, 8, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 36, "ALT_CODE",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 37, "ECO_CODE_X", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 28, "SUBESTUARY", FTString, 32, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 45, "WATER_BODY", FTString, 80, 0, 0.0, 1.0 },
  { "seagrass_point_gulf", 46, "REGION",     FTString, 40, 0, 0.0, 1.0 },

  /* Seagrass point2 data atlantic (updated 2024-01-25): */

  { "seagrass_point2_atlantic",  1, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  4, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  2, "SUBCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  5, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  3, "ESTUARY_NA", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  7, "CUR_YEAR",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  13, "ACRES",      FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  13, "HECTARES",   FTDouble, 16, 6, 0.0,
                                                     ACRES_TO_HECTARES },
  { "seagrass_point2_atlantic", 10, "DATA_COUNT", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic", 11, "CMECS_BIO",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic", 12, "AREA_TYPE",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  8, "EELGRASSUI", FTInteger, 8,  0, 0.0, 1.0 },
  { "seagrass_point2_atlantic",  9, "LINK",       FTString, 48,  0, 0.0, 1.0 },
  { "seagrass_point2_atlantic", 14, "Shape_Leng", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point2_atlantic", 17, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point2_atlantic", 16, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },

  /* Seagrass point2 data gulf (updated 2023-10-18): */

  { "seagrass_point2_gulf",  1, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  4, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  2, "SUBCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  5, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  3, "ESTUARY_NA", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  6, "CUR_YEAR",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  12, "ACRES",      FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point2_gulf",  12, "HECTARES",   FTDouble, 16, 6, 0.0,
                                                     ACRES_TO_HECTARES },
  { "seagrass_point2_gulf",  9, "DATA_COUNT", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf", 10, "CMECS_BIO",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf", 11, "AREA_TYPE",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  7, "EELGRASSUI", FTInteger, 8,  0, 0.0, 1.0 },
  { "seagrass_point2_gulf",  8, "LINK",       FTString, 48,  0, 0.0, 1.0 },
  { "seagrass_point2_gulf", 13, "Shape_Leng", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point2_gulf", 14, "Shape_Area", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_point2_gulf", 18, "LONGITUDE",  FTDouble, 16, 10, 0.0, 1.0 },
  { "seagrass_point2_gulf", 17, "LATITUDE",   FTDouble, 16, 10, 0.0, 1.0 },

  /* Seagrass polygon data pacific (updated 2023-11-14): */

  { "seagrass_polygon_pacific", 11, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  2, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific", 13, "SUBCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific", 14, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific", 12, "ESTUARY_NA", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  4, "CUR_YEAR",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  8, "ACRES",      FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  8, "HECTARES",   FTDouble, 16, 6, 0.0,
                                                     ACRES_TO_HECTARES },
  { "seagrass_polygon_pacific",  5, "DATA_COUNT", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  6, "CMECS_BIO",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  7, "AREA_TYPE",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  1, "EELGRASSUI", FTInteger, 8,  0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  3, "LINK",       FTString, 48,  0, 0.0, 1.0 },
  { "seagrass_polygon_pacific",  9, "Shape_Leng", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_polygon_pacific", 10, "Shape_Area", FTDouble, 16, 6, 0.0, 1.0 },

  /* Seagrass polygon data atlantic and gulf (last update 2023-10-18): */

  { "seagrass_polygon_!pacific",  1, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  4, "ESTUARY",    FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  2, "SUBCODE",    FTString,  5, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  5, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  3, "ESTUARY_NA", FTString, 48, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  7, "CUR_YEAR",   FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  12, "ACRES",      FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  12, "HECTARES",   FTDouble, 16, 6, 0.0,
                                                     ACRES_TO_HECTARES },
  { "seagrass_polygon_!pacific",  8, "DATA_COUNT", FTInteger, 4, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific", 10, "CMECS_BIO",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific", 11, "AREA_TYPE",  FTString, 16, 0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  9, "EELGRASSUI", FTInteger, 8,  0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific",  6, "LINK",       FTString, 48,  0, 0.0, 1.0 },
  { "seagrass_polygon_!pacific", 19, "Shape_Leng", FTDouble, 16, 6, 0.0, 1.0 },
  { "seagrass_polygon_!pacific", 20, "Shape_Area", FTDouble, 16, 6, 0.0, 1.0 },


  /*
   * Population data (single file per coast with all scenario columns)
   * Renamed 2017-02-08:
   */

  { "population_iclus1",  0, "COUNTYFIPS", FTInteger, 5, 0, 0.0, 1.0 },

  { "population_iclus1", -1, "A1_2010PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2020PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2030PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2040PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2050PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2060PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2070PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2080PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A1_2090PKM", FTDouble,  7, 1, 0.0, 1.0 },

  { "population_iclus1", -1, "A2_2010PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2020PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2030PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2040PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2050PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2060PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2070PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2080PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "A2_2090PKM", FTDouble,  7, 1, 0.0, 1.0 },

  { "population_iclus1", -1, "B1_2010PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2020PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2030PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2040PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2050PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2060PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2070PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2080PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B1_2090PKM", FTDouble,  7, 1, 0.0, 1.0 },

  { "population_iclus1", -1, "B2_2010PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2020PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2030PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2040PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2050PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2060PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2070PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2080PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "B2_2090PKM", FTDouble,  7, 1, 0.0, 1.0 },

  { "population_iclus1", -1, "BC_2010PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2020PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2030PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2040PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2050PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2060PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2070PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2080PKM", FTDouble,  7, 1, 0.0, 1.0 },
  { "population_iclus1", -1, "BC_2090PKM", FTDouble,  7, 1, 0.0, 1.0 },

  { "population_iclus1",  2, "A1_2010POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1",  3, "A1_2020POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1",  4, "A1_2030POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1",  5, "A1_2040POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1",  6, "A1_2050POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1",  7, "A1_2060POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1",  8, "A1_2070POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1",  9, "A1_2080POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 10, "A1_2090POP", FTInteger, 9, 0, 0.0, 1.0 },

  { "population_iclus1", 11, "A2_2010POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 12, "A2_2020POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 13, "A2_2030POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 14, "A2_2040POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 15, "A2_2050POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 16, "A2_2060POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 17, "A2_2070POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 18, "A2_2080POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 19, "A2_2090POP", FTInteger, 9, 0, 0.0, 1.0 },

  { "population_iclus1", 20, "B1_2010POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 21, "B1_2020POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 22, "B1_2030POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 23, "B1_2040POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 24, "B1_2050POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 25, "B1_2060POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 26, "B1_2070POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 27, "B1_2080POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 28, "B1_2090POP", FTInteger, 9, 0, 0.0, 1.0 },

  { "population_iclus1", 29, "B2_2010POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 30, "B2_2020POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 31, "B2_2030POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 32, "B2_2040POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 33, "B2_2050POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 34, "B2_2060POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 35, "B2_2070POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 36, "B2_2080POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 37, "B2_2090POP", FTInteger, 9, 0, 0.0, 1.0 },

  { "population_iclus1", 38, "BC_2010POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 39, "BC_2020POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 40, "BC_2030POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 41, "BC_2040POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 42, "BC_2050POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 43, "BC_2060POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 44, "BC_2070POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 45, "BC_2080POP", FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus1", 46, "BC_2090POP", FTInteger, 9, 0, 0.0, 1.0 },

  { "population_iclus1",  1, "COUNTYSQKM", FTDouble, 11, 3, 0.0, 1e-6 },
  { "population_iclus1", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },


  /*
   * Population data (single file per coast with all scenario columns)
   * Added 2017-02-08:
   */

  { "population_iclus2",  0, "ICLUSGEOID", FTInteger, 5, 0, 0.0, 1.0 },
  { "population_iclus2",  6, "POP_1990",   FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2",  7, "POP_2000",   FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2",  8, "POP_2010",   FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2",  9, "SSP2_2020",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 10, "SSP2_2030",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 11, "SSP2_2040",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 12, "SSP2_2050",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 13, "SSP2_2060",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 14, "SSP2_2070",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 15, "SSP2_2080",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 16, "SSP2_2090",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 17, "SSP2_2100",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 18, "SSP5_2020",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 19, "SSP5_2030",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 20, "SSP5_2040",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 21, "SSP5_2050",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 22, "SSP5_2060",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 23, "SSP5_2070",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 24, "SSP5_2080",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 25, "SSP5_2090",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2", 26, "SSP5_2100",  FTInteger, 9, 0, 0.0, 1.0 },
  { "population_iclus2",  2, "FIPS_SQKM",  FTDouble, 11, 3, 0.0, 1e-6 },
  { "population_iclus2", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "population_iclus2",  3, "FIPS_NAME",  FTString, 64, 0, 0.0, 1.0 },
  { "population_iclus2",  4, "STATE",      FTString,  2, 0, 0.0, 1.0 },
  { "population_iclus2",  5, "NCA_REGION", FTString, 16, 0, 0.0, 1.0 },


  /*
   * Land use data (single file per coast with all scenario columns)
   * Added 2017-08-22:
   */

  { "land_use_iclus_",   0, "ESTCODE",     FTString,  4, 0, 0.0, 1.0 },

  { "land_use_iclus_",   4, "LUV0_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",   5, "LUV1_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",   6, "LUV2_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",   7, "LUV3_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",   8, "LUV4_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",   9, "LUV5_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  10, "LUV6_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  11, "LUV7_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  12, "LUV8_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  13, "LUV9_2000",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  14, "LUV10_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  15, "LUV11_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  16, "LUV12_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  17, "LUV13_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  18, "LUV14_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  19, "LUV15_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  20, "LUV16_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  21, "LUV17_2000",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  22, "LUV18_2000",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_",  23, "LUV0_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  24, "LUV1_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  25, "LUV2_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  26, "LUV3_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  27, "LUV4_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  28, "LUV5_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  29, "LUV6_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  30, "LUV7_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  31, "LUV8_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  32, "LUV9_2010",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  33, "LUV10_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  34, "LUV11_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  35, "LUV12_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  36, "LUV13_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  37, "LUV14_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  38, "LUV15_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  39, "LUV16_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  40, "LUV17_2010",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  41, "LUV18_2010",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_",  42, "LUV0_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  43, "LUV1_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  44, "LUV2_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  45, "LUV3_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  46, "LUV4_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  47, "LUV5_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  48, "LUV6_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  49, "LUV7_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  50, "LUV8_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  51, "LUV9_2020",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  52, "LUV10_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  53, "LUV11_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  54, "LUV12_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  55, "LUV13_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  56, "LUV14_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  57, "LUV15_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  58, "LUV16_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  59, "LUV17_2020",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  60, "LUV18_2020",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_",  61, "LUV0_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  62, "LUV1_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  63, "LUV2_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  64, "LUV3_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  65, "LUV4_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  66, "LUV5_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  67, "LUV6_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  68, "LUV7_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  69, "LUV8_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  70, "LUV9_2030",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  71, "LUV10_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  72, "LUV11_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  73, "LUV12_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  74, "LUV13_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  75, "LUV14_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  76, "LUV15_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  77, "LUV16_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  78, "LUV17_2030",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  79, "LUV18_2030",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_",  80, "LUV0_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  81, "LUV1_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  82, "LUV2_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  83, "LUV3_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  84, "LUV4_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  85, "LUV5_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  86, "LUV6_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  87, "LUV7_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  88, "LUV8_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  89, "LUV9_2040",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  90, "LUV10_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  91, "LUV11_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  92, "LUV12_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  93, "LUV13_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  94, "LUV14_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  95, "LUV15_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  96, "LUV16_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  97, "LUV17_2040",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_",  98, "LUV18_2040",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_",  99, "LUV0_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 100, "LUV1_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 101, "LUV2_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 102, "LUV3_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 103, "LUV4_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 104, "LUV5_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 105, "LUV6_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 106, "LUV7_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 107, "LUV8_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 108, "LUV9_2050",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 109, "LUV10_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 110, "LUV11_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 111, "LUV12_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 112, "LUV13_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 113, "LUV14_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 114, "LUV15_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 115, "LUV16_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 116, "LUV17_2050",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 117, "LUV18_2050",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_", 118, "LUV0_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 119, "LUV1_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 120, "LUV2_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 121, "LUV3_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 122, "LUV4_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 123, "LUV5_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 124, "LUV6_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 125, "LUV7_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 126, "LUV8_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 127, "LUV9_2060",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 128, "LUV10_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 129, "LUV11_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 130, "LUV12_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 131, "LUV13_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 132, "LUV14_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 133, "LUV15_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 134, "LUV16_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 135, "LUV17_2060",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 136, "LUV18_2060",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_", 137, "LUV0_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 138, "LUV1_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 139, "LUV2_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 140, "LUV3_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 141, "LUV4_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 142, "LUV5_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 143, "LUV6_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 144, "LUV7_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 145, "LUV8_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 146, "LUV9_2070",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 147, "LUV10_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 148, "LUV11_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 149, "LUV12_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 150, "LUV13_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 151, "LUV14_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 152, "LUV15_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 153, "LUV16_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 154, "LUV17_2070",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 155, "LUV18_2070",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_", 156, "LUV0_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 157, "LUV1_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 158, "LUV2_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 159, "LUV3_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 160, "LUV4_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 161, "LUV5_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 162, "LUV6_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 163, "LUV7_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 164, "LUV8_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 165, "LUV9_2080",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 166, "LUV10_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 167, "LUV11_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 168, "LUV12_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 169, "LUV13_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 170, "LUV14_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 171, "LUV15_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 172, "LUV16_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 173, "LUV17_2080",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 174, "LUV18_2080",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_", 175, "LUV0_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 176, "LUV1_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 177, "LUV2_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 178, "LUV3_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 179, "LUV4_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 180, "LUV5_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 181, "LUV6_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 182, "LUV7_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 183, "LUV8_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 184, "LUV9_2090",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 185, "LUV10_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 186, "LUV11_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 187, "LUV12_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 188, "LUV13_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 189, "LUV14_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 190, "LUV15_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 191, "LUV16_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 192, "LUV17_2090",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 193, "LUV18_2090",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_", 194, "LUV0_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 195, "LUV1_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 196, "LUV2_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 197, "LUV3_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 198, "LUV4_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 199, "LUV5_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 200, "LUV6_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 201, "LUV7_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 202, "LUV8_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 203, "LUV9_2100",   FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 204, "LUV10_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 205, "LUV11_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 206, "LUV12_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 207, "LUV13_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 208, "LUV14_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 209, "LUV15_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 210, "LUV16_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 211, "LUV17_2100",  FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_iclus_", 212, "LUV18_2100",  FTDouble, 10, 6, 0.0, 1.0 },

  { "land_use_iclus_",   2, "FIPS_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "land_use_iclus_",  -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "land_use_iclus_",  -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "land_use_iclus_",   1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },


  /*
   * NLCD land use data 2011 and 2006 (single file per coast non-time-varying)
   * Added 2017-09-12:
   */

  { "land_use_nlcd_",   0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },

  { "land_use_nlcd_",   4, "LUV11", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",   5, "LUV12", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",   6, "LUV21", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",   7, "LUV22", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",   8, "LUV23", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",   9, "LUV24", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  10, "LUV31", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  11, "LUV41", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  12, "LUV42", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  13, "LUV43", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  14, "LUV52", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  15, "LUV71", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  16, "LUV81", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  17, "LUV82", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  18, "LUV90", FTDouble, 10, 6, 0.0, 1.0 },
  { "land_use_nlcd_",  19, "LUV95", FTDouble, 10, 6, 0.0, 1.0 },
 
  { "land_use_nlcd_",   2, "FIPS_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "land_use_nlcd_",  -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "land_use_nlcd_",  -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "land_use_nlcd_",   1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },


  /* Temperature data: */

  { "/temperature_",  0, "STATE_FIPS", FTInteger, 2, 0, 0.0, 1.0 },
  { "/temperature_",  1, "TEMP_C",     FTDouble,  5, 1, -32.0, 5.0 / 9.0 },
  { "/temperature_", -1, "ACRES",      FTDouble, 20, 5, 0.0, 1.0 },
  { "/temperature_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Precipitation data: */

  { "precipitation_",  1, "STATE_FIPS", FTInteger, 2, 0, 0.0, 1.0 },
  { "precipitation_",  0, "PRECIP_IN",  FTDouble,  8, 1, 0.0, 1.0 },
  { "precipitation_",  0, "PRECIP_MM",  FTDouble,  8, 1, 0.0, 25.4 },
  { "precipitation_", -1, "ACRES",      FTDouble, 20, 5, 0.0, 1.0 },
  { "precipitation_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Sediment data: */

  { "sediment_kriged",  3, "SEDIMENT",   FTString, 20, 0, 0.0, 1.0 },
  { "sediment_kriged",  2, "CLASS_CODE", FTInteger, 2, 0, 0.0, 1.0 },
  { "sediment_kriged",  4, "CLASS_NAME", FTString, 48, 0, 0.0, 1.0 },
  { "sediment_kriged",  5, "GLOBAL_ID",  FTString, 40, 0, 0.0, 1.0 },
  { "sediment_kriged", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "sediment_kriged", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* Nitrogen Deposition Estuary CMAQ data (original): */

  { "/nitrogen_estuary_cmaq",   0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",   3, "TOTN_2002",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  18, "TOTN_02_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  19, "TOTN_02_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  20, "TOTN_02_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  21, "TOTN_02_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  22, "TOTN_02_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  23, "TOTN_02_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  24, "TOTN_02_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  25, "TOTN_02_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  26, "TOTN_02_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  81, "TOTN_02_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  82, "TOTN_02_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  83, "TOTN_02_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",   8, "DRYN_2002",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 102, "DRYN_02_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 103, "DRYN_02_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 104, "DRYN_02_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 105, "DRYN_02_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 106, "DRYN_02_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 107, "DRYN_02_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 108, "DRYN_02_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 109, "DRYN_02_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 110, "DRYN_02_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 165, "DRYN_02_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 166, "DRYN_02_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 167, "DRYN_02_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  13, "WETN_2002",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 186, "WETN_02_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 187, "WETN_02_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 188, "WETN_02_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 189, "WETN_02_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 190, "WETN_02_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 191, "WETN_02_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 192, "WETN_02_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 193, "WETN_02_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 194, "WETN_02_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 249, "WETN_02_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 250, "WETN_02_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 251, "WETN_02_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",   4, "TOTN_2003",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  27, "TOTN_03_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  28, "TOTN_03_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  29, "TOTN_03_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  30, "TOTN_03_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  31, "TOTN_03_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  32, "TOTN_03_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  33, "TOTN_03_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  34, "TOTN_03_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  35, "TOTN_03_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  84, "TOTN_03_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  85, "TOTN_03_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  86, "TOTN_03_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",   9, "DRYN_2003",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 111, "DRYN_03_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 112, "DRYN_03_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 113, "DRYN_03_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 114, "DRYN_03_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 115, "DRYN_03_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 116, "DRYN_03_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 117, "DRYN_03_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 118, "DRYN_03_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 119, "DRYN_03_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 168, "DRYN_03_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 169, "DRYN_03_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 170, "DRYN_03_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  14, "WETN_2003",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 195, "WETN_03_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 196, "WETN_03_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 197, "WETN_03_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 198, "WETN_03_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 199, "WETN_03_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 200, "WETN_03_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 201, "WETN_03_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 202, "WETN_03_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 203, "WETN_03_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 252, "WETN_03_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 253, "WETN_03_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 254, "WETN_03_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",   5, "TOTN_2004",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  36, "TOTN_04_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  37, "TOTN_04_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  38, "TOTN_04_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  39, "TOTN_04_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  40, "TOTN_04_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  41, "TOTN_04_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  42, "TOTN_04_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  43, "TOTN_04_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  44, "TOTN_04_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  87, "TOTN_04_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  88, "TOTN_04_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  89, "TOTN_04_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  10, "DRYN_2004",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 120, "DRYN_04_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 121, "DRYN_04_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 122, "DRYN_04_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 123, "DRYN_04_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 124, "DRYN_04_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 125, "DRYN_04_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 126, "DRYN_04_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 127, "DRYN_04_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 128, "DRYN_04_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 171, "DRYN_04_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 172, "DRYN_04_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 173, "DRYN_04_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  15, "WETN_2004",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 204, "WETN_04_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 205, "WETN_04_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 206, "WETN_04_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 207, "WETN_04_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 208, "WETN_04_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 209, "WETN_04_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 210, "WETN_04_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 211, "WETN_04_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 212, "WETN_04_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 255, "WETN_04_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 256, "WETN_04_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 257, "WETN_04_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",   6, "TOTN_2005",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  45, "TOTN_05_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  46, "TOTN_05_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  47, "TOTN_05_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  48, "TOTN_05_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  49, "TOTN_05_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  50, "TOTN_05_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  51, "TOTN_05_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  52, "TOTN_05_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  53, "TOTN_05_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  90, "TOTN_05_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  91, "TOTN_05_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  92, "TOTN_05_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  11, "DRYN_2005",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 129, "DRYN_05_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 130, "DRYN_05_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 131, "DRYN_05_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 132, "DRYN_05_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 133, "DRYN_05_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 134, "DRYN_05_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 135, "DRYN_05_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 136, "DRYN_05_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 137, "DRYN_05_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 174, "DRYN_05_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 175, "DRYN_05_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 176, "DRYN_05_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  16, "WETN_2005",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 213, "WETN_05_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 214, "WETN_05_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 215, "WETN_05_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 216, "WETN_05_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 217, "WETN_05_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 218, "WETN_05_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 219, "WETN_05_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 220, "WETN_05_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 221, "WETN_05_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 258, "WETN_05_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 259, "WETN_05_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 260, "WETN_05_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",   7, "TOTN_2006",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  54, "TOTN_06_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  55, "TOTN_06_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  56, "TOTN_06_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  57, "TOTN_06_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  58, "TOTN_06_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  59, "TOTN_06_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  60, "TOTN_06_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  61, "TOTN_06_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  62, "TOTN_06_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  93, "TOTN_06_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  94, "TOTN_06_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  95, "TOTN_06_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  12, "DRYN_2006",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 138, "DRYN_06_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 139, "DRYN_06_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 140, "DRYN_06_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 141, "DRYN_06_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 142, "DRYN_06_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 143, "DRYN_06_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 144, "DRYN_06_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 145, "DRYN_06_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 146, "DRYN_06_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 177, "DRYN_06_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 178, "DRYN_06_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 179, "DRYN_06_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  17, "WETN_2006",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 222, "WETN_06_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 223, "WETN_06_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 224, "WETN_06_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 225, "WETN_06_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 226, "WETN_06_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 227, "WETN_06_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 228, "WETN_06_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 229, "WETN_06_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 230, "WETN_06_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 261, "WETN_06_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 262, "WETN_06_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 263, "WETN_06_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  -1, "TOTN_2007",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  63, "TOTN_07_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  64, "TOTN_07_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  65, "TOTN_07_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  66, "TOTN_07_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  67, "TOTN_07_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  68, "TOTN_07_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  69, "TOTN_07_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  70, "TOTN_07_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  71, "TOTN_07_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  96, "TOTN_07_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  97, "TOTN_07_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  98, "TOTN_07_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  -1, "DRYN_2007",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 147, "DRYN_07_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 148, "DRYN_07_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 149, "DRYN_07_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 150, "DRYN_07_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 151, "DRYN_07_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 152, "DRYN_07_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 153, "DRYN_07_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 154, "DRYN_07_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 155, "DRYN_07_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 180, "DRYN_07_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 181, "DRYN_07_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 182, "DRYN_07_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  -1, "WETN_2007",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 231, "WETN_07_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 232, "WETN_07_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 233, "WETN_07_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 234, "WETN_07_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 235, "WETN_07_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 236, "WETN_07_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 237, "WETN_07_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 238, "WETN_07_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 239, "WETN_07_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 264, "WETN_07_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 265, "WETN_07_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 266, "WETN_07_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  -1, "TOTN_2008",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  72, "TOTN_08_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  73, "TOTN_08_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  74, "TOTN_08_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  75, "TOTN_08_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  76, "TOTN_08_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  77, "TOTN_08_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  78, "TOTN_08_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  79, "TOTN_08_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  80, "TOTN_08_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  99, "TOTN_08_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 100, "TOTN_08_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 101, "TOTN_08_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  -1, "DRYN_2008",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 156, "DRYN_08_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 157, "DRYN_08_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 158, "DRYN_08_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 159, "DRYN_08_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 160, "DRYN_08_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 161, "DRYN_08_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 162, "DRYN_08_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 163, "DRYN_08_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 164, "DRYN_08_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 183, "DRYN_08_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 184, "DRYN_08_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 185, "DRYN_08_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  -1, "WETN_2008",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 240, "WETN_08_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 241, "WETN_08_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 242, "WETN_08_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 243, "WETN_08_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 244, "WETN_08_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 245, "WETN_08_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 246, "WETN_08_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 247, "WETN_08_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 248, "WETN_08_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 267, "WETN_08_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 268, "WETN_08_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq", 269, "WETN_08_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_estuary_cmaq",  -1, "UNITS",      FTString, 11, 0, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "/nitrogen_estuary_cmaq",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },


  /* CMAQ estuary data (new replacement version 2016-03-30): */

  { "chloride_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "chloride_estuary_cmaq_",  5, "CL_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_",  6, "CL_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_",  7, "CL_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_",  8, "CL_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_",  9, "CL_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 10, "CL_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 11, "CL_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 12, "CL_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 13, "CL_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 14, "CL_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 15, "CL_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 16, "CL_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 17, "CL_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 18, "CL_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 19, "CL_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 20, "CL_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 21, "CL_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 22, "CL_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 23, "CL_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 24, "CL_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 25, "CL_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 26, "CL_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 27, "CL_25_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", 28, "CL_25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "chloride_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "chloride_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nh3_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_",  5, "NH3_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_",  6, "NH302KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_",  7, "NH3_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_",  8, "NH303KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_",  9, "NH3_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 10, "NH304KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 11, "NH3_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 12, "NH305KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 13, "NH3_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 14, "NH306KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 15, "NH3_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 16, "NH307KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 17, "NH3_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 18, "NH308KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 19, "NH3_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 20, "NH309KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 21, "NH3_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 22, "NH310KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 23, "NH3_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 24, "NH311KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 25, "NH3_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 26, "NH312KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 27, "NH3_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", 28, "NH325KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nh3_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nhx_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_",  5, "NHX_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_",  6, "NHX02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_",  7, "NHX_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_",  8, "NHX03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_",  9, "NHX_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 10, "NHX04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 11, "NHX_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 12, "NHX05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 13, "NHX_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 14, "NHX06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 15, "NHX_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 16, "NHX07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 17, "NHX_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 18, "NHX08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 19, "NHX_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 20, "NHX09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 21, "NHX_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 22, "NHX10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 23, "NHX_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 24, "NHX11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 25, "NHX_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 26, "NHX12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 27, "NHX_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", 28, "NHX25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nhx_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_no3_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_no3_estuary_cmaq_",  5, "NO3_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_",  6, "NO302KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_",  7, "NO3_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_",  8, "NO303KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_",  9, "NO3_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 10, "NO304KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 11, "NO3_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 12, "NO305KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 13, "NO3_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 14, "NO306KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 15, "NO3_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 16, "NO307KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 17, "NO3_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 18, "NO308KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 19, "NO3_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 20, "NO309KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 21, "NO3_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 22, "NO310KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 23, "NO3_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 24, "NO311KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 25, "NO3_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 26, "NO312KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 27, "NO3_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", 28, "NO325KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_no3_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_no3_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nox_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nox_estuary_cmaq_",  5, "NOX_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_",  6, "NOX02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_",  7, "NOX_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_",  8, "NOX03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_",  9, "NOX_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 10, "NOX04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 11, "NOX_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 12, "NOX05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 13, "NOX_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 14, "NOX06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 15, "NOX_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 16, "NOX07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 17, "NOX_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 18, "NOX08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 19, "NOX_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 20, "NOX09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 21, "NOX_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 22, "NOX10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 23, "NOX_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 24, "NOX11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 25, "NOX_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 26, "NOX12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 27, "NOX_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", 28, "NOX25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_nox_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nox_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nitrogen_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_",  5, "N_02_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_",  6, "N_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_",  7, "N_03_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_",  8, "N_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_",  9, "N_04_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 10, "N_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 11, "N_05_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 12, "N_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 13, "N_06_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 14, "N_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 15, "N_07_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 16, "N_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 17, "N_08_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 18, "N_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 19, "N_09_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 20, "N_09_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 21, "N_10_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 22, "N_10_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 23, "N_11_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 24, "N_11_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 25, "N_12_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 26, "N_12_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 27, "N_25_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", 28, "N_25_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0},
  { "_nitrogen_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0},

  { "_so2_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_so2_estuary_cmaq_",  5, "SO2_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_",  6, "SO202KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_",  7, "SO2_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_",  8, "SO203KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_",  9, "SO2_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 10, "SO204KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 11, "SO2_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 12, "SO205KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 13, "SO2_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 14, "SO206KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 15, "SO2_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 16, "SO207KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 17, "SO2_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 18, "SO208KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 19, "SO2_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 20, "SO209KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 21, "SO2_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 22, "SO210KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 23, "SO2_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 24, "SO211KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 25, "SO2_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 26, "SO212KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 27, "SO2_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", 28, "SO225KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_so2_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_so2_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_so4_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_so4_estuary_cmaq_",  5, "SO4_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_",  6, "SO402KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_",  7, "SO4_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_",  8, "SO403KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_",  9, "SO4_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 10, "SO404KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 11, "SO4_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 12, "SO405KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 13, "SO4_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 14, "SO406KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 15, "SO4_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 16, "SO407KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 17, "SO4_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 18, "SO408KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 19, "SO4_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 20, "SO409KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 21, "SO4_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 22, "SO410KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 23, "SO4_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 24, "SO411KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 25, "SO4_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 26, "SO412KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 27, "SO4_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", 28, "SO425KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_so4_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_so4_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "sulfur_estuary_cmaq_",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_",  5, "S_02_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_",  6, "S_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_",  7, "S_03_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_",  8, "S_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_",  9, "S_04_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 10, "S_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 11, "S_05_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 12, "S_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 13, "S_06_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 14, "S_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 15, "S_07_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 16, "S_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 17, "S_08_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 18, "S_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 19, "S_09_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 20, "S_09_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 21, "S_10_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 22, "S_10_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 23, "S_11_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 24, "S_11_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 25, "S_12_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 26, "S_12_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 27, "S_25_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", 28, "S_25_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_",  4, "LATITUDE",   FTDouble, 11, 6, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "sulfur_estuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },


  /* Nitrogen Deposition Subestuary CMAQ data (original): */

  { "/nitrogen_subestuary_cmaq",   0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   1, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   3, "TOTN_2002",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  18, "TOTN_02_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  19, "TOTN_02_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  20, "TOTN_02_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  21, "TOTN_02_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  22, "TOTN_02_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  23, "TOTN_02_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  24, "TOTN_02_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  25, "TOTN_02_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  26, "TOTN_02_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  81, "TOTN_02_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  82, "TOTN_02_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  83, "TOTN_02_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   8, "DRYN_2002",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 102, "DRYN_02_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 103, "DRYN_02_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 104, "DRYN_02_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 105, "DRYN_02_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 106, "DRYN_02_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 107, "DRYN_02_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 108, "DRYN_02_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 109, "DRYN_02_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 110, "DRYN_02_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 165, "DRYN_02_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 166, "DRYN_02_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 167, "DRYN_02_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  13, "WETN_2002",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 186, "WETN_02_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 187, "WETN_02_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 188, "WETN_02_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 189, "WETN_02_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 190, "WETN_02_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 191, "WETN_02_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 192, "WETN_02_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 193, "WETN_02_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 194, "WETN_02_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 249, "WETN_02_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 250, "WETN_02_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 251, "WETN_02_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   4, "TOTN_2003",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  27, "TOTN_03_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  28, "TOTN_03_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  29, "TOTN_03_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  30, "TOTN_03_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  31, "TOTN_03_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  32, "TOTN_03_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  33, "TOTN_03_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  34, "TOTN_03_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  35, "TOTN_03_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  84, "TOTN_03_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  85, "TOTN_03_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  86, "TOTN_03_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   9, "DRYN_2003",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 111, "DRYN_03_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 112, "DRYN_03_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 113, "DRYN_03_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 114, "DRYN_03_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 115, "DRYN_03_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 116, "DRYN_03_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 117, "DRYN_03_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 118, "DRYN_03_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 119, "DRYN_03_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 168, "DRYN_03_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 169, "DRYN_03_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 170, "DRYN_03_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  14, "WETN_2003",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 195, "WETN_03_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 196, "WETN_03_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 197, "WETN_03_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 198, "WETN_03_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 199, "WETN_03_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 200, "WETN_03_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 201, "WETN_03_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 202, "WETN_03_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 203, "WETN_03_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 252, "WETN_03_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 253, "WETN_03_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 254, "WETN_03_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   5, "TOTN_2004",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  36, "TOTN_04_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  37, "TOTN_04_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  38, "TOTN_04_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  39, "TOTN_04_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  40, "TOTN_04_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  41, "TOTN_04_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  42, "TOTN_04_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  43, "TOTN_04_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  44, "TOTN_04_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  87, "TOTN_04_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  88, "TOTN_04_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  89, "TOTN_04_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  10, "DRYN_2004",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 120, "DRYN_04_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 121, "DRYN_04_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 122, "DRYN_04_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 123, "DRYN_04_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 124, "DRYN_04_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 125, "DRYN_04_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 126, "DRYN_04_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 127, "DRYN_04_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 128, "DRYN_04_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 171, "DRYN_04_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 172, "DRYN_04_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 173, "DRYN_04_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  15, "WETN_2004",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 204, "WETN_04_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 205, "WETN_04_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 206, "WETN_04_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 207, "WETN_04_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 208, "WETN_04_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 209, "WETN_04_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 210, "WETN_04_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 211, "WETN_04_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 212, "WETN_04_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 255, "WETN_04_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 256, "WETN_04_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 257, "WETN_04_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   6, "TOTN_2005",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  45, "TOTN_05_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  46, "TOTN_05_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  47, "TOTN_05_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  48, "TOTN_05_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  49, "TOTN_05_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  50, "TOTN_05_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  51, "TOTN_05_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  52, "TOTN_05_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  53, "TOTN_05_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  90, "TOTN_05_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  91, "TOTN_05_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  92, "TOTN_05_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  11, "DRYN_2005",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 129, "DRYN_05_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 130, "DRYN_05_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 131, "DRYN_05_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 132, "DRYN_05_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 133, "DRYN_05_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 134, "DRYN_05_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 135, "DRYN_05_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 136, "DRYN_05_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 137, "DRYN_05_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 174, "DRYN_05_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 175, "DRYN_05_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 176, "DRYN_05_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  16, "WETN_2005",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 213, "WETN_05_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 214, "WETN_05_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 215, "WETN_05_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 216, "WETN_05_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 217, "WETN_05_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 218, "WETN_05_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 219, "WETN_05_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 220, "WETN_05_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 221, "WETN_05_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 258, "WETN_05_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 259, "WETN_05_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 260, "WETN_05_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",   7, "TOTN_2006",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  54, "TOTN_06_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  55, "TOTN_06_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  56, "TOTN_06_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  57, "TOTN_06_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  58, "TOTN_06_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  59, "TOTN_06_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  60, "TOTN_06_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  61, "TOTN_06_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  62, "TOTN_06_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  93, "TOTN_06_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  94, "TOTN_06_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  95, "TOTN_06_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  12, "DRYN_2006",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 138, "DRYN_06_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 139, "DRYN_06_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 140, "DRYN_06_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 141, "DRYN_06_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 142, "DRYN_06_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 143, "DRYN_06_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 144, "DRYN_06_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 145, "DRYN_06_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 146, "DRYN_06_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 177, "DRYN_06_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 178, "DRYN_06_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 179, "DRYN_06_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  17, "WETN_2006",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 222, "WETN_06_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 223, "WETN_06_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 224, "WETN_06_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 225, "WETN_06_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 226, "WETN_06_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 227, "WETN_06_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 228, "WETN_06_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 229, "WETN_06_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 230, "WETN_06_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 261, "WETN_06_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 262, "WETN_06_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 263, "WETN_06_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  -1, "TOTN_2007",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  63, "TOTN_07_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  64, "TOTN_07_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  65, "TOTN_07_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  66, "TOTN_07_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  67, "TOTN_07_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  68, "TOTN_07_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  69, "TOTN_07_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  70, "TOTN_07_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  71, "TOTN_07_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  96, "TOTN_07_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  97, "TOTN_07_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  98, "TOTN_07_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  -1, "DRYN_2007",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 147, "DRYN_07_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 148, "DRYN_07_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 149, "DRYN_07_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 150, "DRYN_07_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 151, "DRYN_07_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 152, "DRYN_07_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 153, "DRYN_07_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 154, "DRYN_07_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 155, "DRYN_07_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 180, "DRYN_07_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 181, "DRYN_07_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 182, "DRYN_07_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  -1, "WETN_2007",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 231, "WETN_07_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 232, "WETN_07_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 233, "WETN_07_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 234, "WETN_07_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 235, "WETN_07_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 236, "WETN_07_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 237, "WETN_07_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 238, "WETN_07_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 239, "WETN_07_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 264, "WETN_07_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 265, "WETN_07_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 266, "WETN_07_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  -1, "TOTN_2008",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  72, "TOTN_08_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  73, "TOTN_08_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  74, "TOTN_08_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  75, "TOTN_08_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  76, "TOTN_08_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  77, "TOTN_08_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  78, "TOTN_08_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  79, "TOTN_08_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  80, "TOTN_08_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  99, "TOTN_08_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 100, "TOTN_08_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 101, "TOTN_08_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  -1, "DRYN_2008",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 156, "DRYN_08_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 157, "DRYN_08_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 158, "DRYN_08_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 159, "DRYN_08_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 160, "DRYN_08_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 161, "DRYN_08_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 162, "DRYN_08_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 163, "DRYN_08_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 164, "DRYN_08_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 183, "DRYN_08_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 184, "DRYN_08_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 185, "DRYN_08_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  -1, "WETN_2008",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 240, "WETN_08_01", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 241, "WETN_08_02", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 242, "WETN_08_03", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 243, "WETN_08_04", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 244, "WETN_08_05", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 245, "WETN_08_06", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 246, "WETN_08_07", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 247, "WETN_08_08", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 248, "WETN_08_09", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 267, "WETN_08_10", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 268, "WETN_08_11", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq", 269, "WETN_08_12", FTDouble, 16, 4, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",  -1, "UNITS",      FTString, 11, 0, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "/nitrogen_subestuary_cmaq",  -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "/nitrogen_subestuary_cmaq",2, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },


  /* CMAQ subestuary data (new replacement version 2016-03-30): */

  { "chloride_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_",  4, "CL_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_",  5, "CL_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_",  6, "CL_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_",  7, "CL_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_",  8, "CL_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_",  9, "CL_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 10, "CL_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 11, "CL_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 12, "CL_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 13, "CL_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 14, "CL_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 15, "CL_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 16, "CL_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 17, "CL_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 18, "CL_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 19, "CL_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 20, "CL_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 21, "CL_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 22, "CL_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 23, "CL_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 24, "CL_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 25, "CL_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 26, "CL_25_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", 27, "CL_25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "chloride_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nh3_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_",  4, "NH3_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_",  5, "NH302KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_",  6, "NH3_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_",  7, "NH303KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_",  8, "NH3_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_",  9, "NH304KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 10, "NH3_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 11, "NH305KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 12, "NH3_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 13, "NH306KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 14, "NH3_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 15, "NH307KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 16, "NH3_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 17, "NH308KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 18, "NH3_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 19, "NH309KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 20, "NH3_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 21, "NH310KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 22, "NH3_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 23, "NH311KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 24, "NH3_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 25, "NH312KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 26, "NH3_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", 27, "NH325KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nh3_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nhx_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_",  4, "NHX_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_",  5, "NHX02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_",  6, "NHX_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_",  7, "NHX03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_",  8, "NHX_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_",  9, "NHX04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 10, "NHX_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 11, "NHX05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 12, "NHX_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 13, "NHX06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 14, "NHX_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 15, "NHX07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 16, "NHX_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 17, "NHX08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 18, "NHX_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 19, "NHX09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 20, "NHX_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 21, "NHX10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 22, "NHX_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 23, "NHX11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 24, "NHX_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 25, "NHX12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 26, "NHX_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", 27, "NHX25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nhx_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_no3_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_",  4, "NO3_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_",  5, "NO302KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_",  6, "NO3_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_",  7, "NO303KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_",  8, "NO3_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_",  9, "NO304KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 10, "NO3_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 11, "NO305KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 12, "NO3_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 13, "NO306KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 14, "NO3_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 15, "NO307KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 16, "NO3_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 17, "NO308KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 18, "NO3_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 19, "NO309KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 20, "NO3_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 21, "NO310KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 22, "NO3_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 23, "NO311KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 24, "NO3_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 25, "NO312KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 26, "NO3_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", 27, "NO325KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_no3_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nox_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_",  4, "NOX_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_",  5, "NOX02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_",  6, "NOX_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_",  7, "NOX03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_",  8, "NOX_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_",  9, "NOX04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 10, "NOX_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 11, "NOX05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 12, "NOX_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 13, "NOX06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 14, "NOX_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 15, "NOX07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 16, "NOX_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 17, "NOX08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 18, "NOX_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 19, "NOX09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 20, "NOX_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 21, "NOX10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 22, "NOX_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 23, "NOX11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 24, "NOX_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 25, "NOX12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 26, "NOX_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", 27, "NOX25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nox_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_nitrogen_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_",  4, "N_02_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_",  5, "N_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_",  6, "N_03_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_",  7, "N_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_",  8, "N_04_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_",  9, "N_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 10, "N_05_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 11, "N_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 12, "N_06_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 13, "N_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 14, "N_07_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 15, "N_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 16, "N_08_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 17, "N_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 18, "N_09_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 19, "N_09_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 20, "N_10_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 21, "N_10_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 22, "N_11_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 23, "N_11_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 24, "N_12_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 25, "N_12_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 26, "N_25_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", 27, "N_25_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0},
  { "_nitrogen_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0},

  { "_so2_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_",  4, "SO2_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_",  5, "SO202KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_",  6, "SO2_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_",  7, "SO203KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_",  8, "SO2_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_",  9, "SO204KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 10, "SO2_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 11, "SO205KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 12, "SO2_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 13, "SO206KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 14, "SO2_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 15, "SO207KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 16, "SO2_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 17, "SO208KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 18, "SO2_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 19, "SO209KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 20, "SO2_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 21, "SO210KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 22, "SO2_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 23, "SO211KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 24, "SO2_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 25, "SO212KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 26, "SO2_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", 27, "SO225KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_so2_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "_so4_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_",  4, "SO4_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_",  5, "SO402KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_",  6, "SO4_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_",  7, "SO403KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_",  8, "SO4_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_",  9, "SO404KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 10, "SO4_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 11, "SO405KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 12, "SO4_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 13, "SO406KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 14, "SO4_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 15, "SO407KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 16, "SO4_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 17, "SO408KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 18, "SO4_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 19, "SO409KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 20, "SO4_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 21, "SO410KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 22, "SO4_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 23, "SO411KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 24, "SO4_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 25, "SO412KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 26, "SO4_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", 27, "SO425KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_so4_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  { "sulfur_subestuary_cmaq_",  0, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_",  1, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_",  4, "S_02_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_",  5, "S_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_",  6, "S_03_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_",  7, "S_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_",  8, "S_04_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_",  9, "S_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 10, "S_05_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 11, "S_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 12, "S_06_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 13, "S_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 14, "S_07_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 15, "S_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 16, "S_08_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 17, "S_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 18, "S_09_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 19, "S_09_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 20, "S_10_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 21, "S_10_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 22, "S_11_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 23, "S_11_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 24, "S_12_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 25, "S_12_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 26, "S_25_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", 27, "S_25_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "sulfur_subestuary_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },



  /*
   * Nitrogen Source Watershed NADP data (replaced 2016-11-21).
   * Note: wet_ files have more years (1985-2014) than rest (2000-2014).
   */

  { "total_nitrogen_source_watershed_nadp",  0, "ESTCODE",    FTString,  8, 0, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp",  4, "TN_00_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp",  5, "TN_00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp",  6, "TN_00_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp",  7, "TN_01_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp",  8, "TN_01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp",  9, "TN_01_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 10, "TN_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 11, "TN_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 12, "TN_02_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 13, "TN_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 14, "TN_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 15, "TN_03_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 16, "TN_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 17, "TN_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 18, "TN_04_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 19, "TN_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 20, "TN_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 21, "TN_05_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 22, "TN_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 23, "TN_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 24, "TN_06_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 25, "TN_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 26, "TN_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 27, "TN_07_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 28, "TN_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 29, "TN_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 30, "TN_08_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 31, "TN_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 32, "TN_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 33, "TN_09_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 34, "TN_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 35, "TN_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 36, "TN_10_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 37, "TN_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 38, "TN_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 39, "TN_11_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 40, "TN_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 41, "TN_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 42, "TN_12_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 43, "TN_13_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 44, "TN_13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 45, "TN_13_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 46, "TN_14_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 47, "TN_14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", 48, "TN_14_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "total_nitrogen_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "dry_nitrogen_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp",  4, "DN_00_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp",  5, "DN_00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp",  6, "DN_00_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp",  7, "DN_01_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp",  8, "DN_01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp",  9, "DN_01_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 10, "DN_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 11, "DN_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 12, "DN_02_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 13, "DN_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 14, "DN_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 15, "DN_03_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 16, "DN_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 17, "DN_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 18, "DN_04_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 19, "DN_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 20, "DN_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 21, "DN_05_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 22, "DN_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 23, "DN_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 24, "DN_06_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 25, "DN_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 26, "DN_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 27, "DN_07_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 28, "DN_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 29, "DN_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 30, "DN_08_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 31, "DN_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 32, "DN_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 33, "DN_09_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 34, "DN_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 35, "DN_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 36, "DN_10_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 37, "DN_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 38, "DN_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 39, "DN_11_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 40, "DN_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 41, "DN_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 42, "DN_12_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 43, "DN_13_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 44, "DN_13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 45, "DN_13_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 46, "DN_14_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 47, "DN_14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", 48, "DN_14_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "dry_nitrogen_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* / indicates match beginning of file name (to avoid conflict below). */

  { "/oxidized_nitrogen_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp",  4, "ON_00_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp",  5, "ON_00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp",  6, "ON_00_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp",  7, "ON_01_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp",  8, "ON_01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp",  9, "ON_01_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 10, "ON_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 11, "ON_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 12, "ON_02_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 13, "ON_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 14, "ON_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 15, "ON_03_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 16, "ON_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 17, "ON_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 18, "ON_04_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 19, "ON_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 20, "ON_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 21, "ON_05_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 22, "ON_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 23, "ON_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 24, "ON_06_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 25, "ON_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 26, "ON_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 27, "ON_07_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 28, "ON_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 29, "ON_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 30, "ON_08_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 31, "ON_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 32, "ON_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 33, "ON_09_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 34, "ON_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 35, "ON_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 36, "ON_10_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 37, "ON_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 38, "ON_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 39, "ON_11_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 40, "ON_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 41, "ON_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 42, "ON_12_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 43, "ON_13_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 44, "ON_13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 45, "ON_13_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 46, "ON_14_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 47, "ON_14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", 48, "ON_14_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "/oxidized_nitrogen_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "/reduced_nitrogen_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp",  4, "RN_00_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp",  5, "RN_00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp",  6, "RN_00_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp",  7, "RN_01_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp",  8, "RN_01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp",  9, "RN_01_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 10, "RN_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 11, "RN_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 12, "RN_02_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 13, "RN_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 14, "RN_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 15, "RN_03_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 16, "RN_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 17, "RN_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 18, "RN_04_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 19, "RN_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 20, "RN_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 21, "RN_05_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 22, "RN_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 23, "RN_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 24, "RN_06_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 25, "RN_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 26, "RN_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 27, "RN_07_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 28, "RN_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 29, "RN_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 30, "RN_08_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 31, "RN_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 32, "RN_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 33, "RN_09_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 34, "RN_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 35, "RN_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 36, "RN_10_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 37, "RN_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 38, "RN_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 39, "RN_11_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 40, "RN_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 41, "RN_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 42, "RN_12_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 43, "RN_13_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 44, "RN_13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 45, "RN_13_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 46, "RN_14_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 47, "RN_14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", 48, "RN_14_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "/reduced_nitrogen_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "dry_oxidized_nitrogen_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp",  4, "DON_00_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp",  5, "DON00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp",  6, "DON_00_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp",  7, "DON_01_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp",  8, "DON01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp",  9, "DON_01_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 10, "DON_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 11, "DON02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 12, "DON_02_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 13, "DON_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 14, "DON03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 15, "DON_03_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 16, "DON_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 17, "DON04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 18, "DON_04_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 19, "DON_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 20, "DON05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 21, "DON_05_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 22, "DON_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 23, "DON06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 24, "DON_06_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 25, "DON_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 26, "DON07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 27, "DON_07_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 28, "DON_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 29, "DON08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 30, "DON_08_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 31, "DON_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 32, "DON09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 33, "DON_09_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 34, "DON_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 35, "DON10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 36, "DON_10_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 37, "DON_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 38, "DON11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 39, "DON_11_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 40, "DON_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 41, "DON12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 42, "DON_12_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 43, "DON_13_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 44, "DON13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 45, "DON_13_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 46, "DON_14_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 47, "DON14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", 48, "DON_14_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "dry_oxidized_nitrogen_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "dry_reduced_nitrogen_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp",  4, "DRN_00_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp",  5, "DRN00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp",  6, "DRN_00_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp",  7, "DRN_01_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp",  8, "DRN01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp",  9, "DRN_01_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 10, "DRN_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 11, "DRN02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 12, "DRN_02_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 13, "DRN_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 14, "DRN03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 15, "DRN_03_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 16, "DRN_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 17, "DRN04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 18, "DRN_04_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 19, "DRN_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 20, "DRN05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 21, "DRN_05_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 22, "DRN_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 23, "DRN06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 24, "DRN_06_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 25, "DRN_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 26, "DRN07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 27, "DRN_07_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 28, "DRN_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 29, "DRN08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 30, "DRN_08_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 31, "DRN_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 32, "DRN09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 33, "DRN_09_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 34, "DRN_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 35, "DRN10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 36, "DRN_10_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 37, "DRN_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 38, "DRN11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 39, "DRN_11_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 40, "DRN_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 41, "DRN12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 42, "DRN_12_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 43, "DRN_13_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 44, "DRN13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 45, "DRN_13_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 46, "DRN_14_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 47, "DRN14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", 48, "DRN_14_%",   FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "dry_reduced_nitrogen_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "wet_nitrogen_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp",  4, "WN_85_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp",  5, "WN_85KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp",  6, "WN_85_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp",  7, "WN_86_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp",  8, "WN_86KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp",  9, "WN_86_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 10, "WN_87_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 11, "WN_87KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 12, "WN_87_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 13, "WN_88_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 14, "WN_88KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 15, "WN_88_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 16, "WN_89_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 17, "WN_89KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 18, "WN_89_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 19, "WN_90_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 20, "WN_90KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 21, "WN_90_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 22, "WN_91_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 23, "WN_91KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 24, "WN_91_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 25, "WN_92_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 26, "WN_92KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 27, "WN_92_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 28, "WN_93_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 29, "WN_93KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 30, "WN_93_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 31, "WN_94_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 32, "WN_94KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 33, "WN_94_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 34, "WN_95_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 35, "WN_95KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 36, "WN_95_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 37, "WN_96_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 38, "WN_96KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 39, "WN_96_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 40, "WN_97_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 41, "WN_97KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 42, "WN_97_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 43, "WN_98_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 44, "WN_98KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 45, "WN_98_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 46, "WN_99_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 47, "WN_99KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 48, "WN_99_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 49, "WN_00_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 50, "WN_00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 51, "WN_00_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 52, "WN_01_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 53, "WN_01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 54, "WN_01_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 55, "WN_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 56, "WN_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 57, "WN_02_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 58, "WN_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 59, "WN_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 60, "WN_03_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 61, "WN_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 62, "WN_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 63, "WN_04_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 64, "WN_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 65, "WN_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 66, "WN_05_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 67, "WN_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 68, "WN_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 69, "WN_06_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 70, "WN_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 71, "WN_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 72, "WN_07_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 73, "WN_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 74, "WN_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 75, "WN_08_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 76, "WN_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 77, "WN_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 78, "WN_09_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 79, "WN_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 80, "WN_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 81, "WN_10_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 82, "WN_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 83, "WN_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 84, "WN_11_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 85, "WN_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 86, "WN_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 87, "WN_12_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 88, "WN_13_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 89, "WN_13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 90, "WN_13_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 91, "WN_14_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 92, "WN_14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", 93, "WN_14_%",    FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "wet_nitrogen_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "dry_nh4_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp",  4, "DNH400_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp",  5, "DNH00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp",  6, "DNH4_00_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp",  7, "DNH401_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp",  8, "DNH01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp",  9, "DNH4_01_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 10, "DNH402_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 11, "DNH02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 12, "DNH4_02_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 13, "DNH403_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 14, "DNH03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 15, "DNH4_03_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 16, "DNH404_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 17, "DNH04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 18, "DNH4_04_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 19, "DNH405_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 20, "DNH05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 21, "DNH4_05_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 22, "DNH406_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 23, "DNH06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 24, "DNH4_06_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 25, "DNH407_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 26, "DNH07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 27, "DNH4_07_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 28, "DNH408_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 29, "DNH08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 30, "DNH4_08_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 31, "DNH409_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 32, "DNH09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 33, "DNH4_09_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 34, "DNH410_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 35, "DNH10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 36, "DNH4_10_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 37, "DNH411_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 38, "DNH11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 39, "DNH4_11_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 40, "DNH412_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 41, "DNH12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 42, "DNH4_12_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 43, "DNH413_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 44, "DNH13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 45, "DNH4_13_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 46, "DNH414_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 47, "DNH14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", 48, "DNH4_14_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "dry_nh4_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "dry_no3_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp",  4, "DNO300_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp",  5, "DNO00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp",  6, "DNO3_00_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp",  7, "DNO301_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp",  8, "DNO01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp",  9, "DNO3_01_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 10, "DNO302_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 11, "DNO02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 12, "DNO3_02_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 13, "DNO303_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 14, "DNO03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 15, "DNO3_03_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 16, "DNO304_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 17, "DNO04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 18, "DNO3_04_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 19, "DNO305_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 20, "DNO05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 21, "DNO3_05_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 22, "DNO306_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 23, "DNO06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 24, "DNO3_06_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 25, "DNO307_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 26, "DNO07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 27, "DNO3_07_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 28, "DNO308_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 29, "DNO08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 30, "DNO3_08_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 31, "DNO309_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 32, "DNO09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 33, "DNO3_09_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 34, "DNO310_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 35, "DNO10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 36, "DNO3_10_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 37, "DNO311_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 38, "DNO11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 39, "DNO3_11_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 40, "DNO312_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 41, "DNO12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 42, "DNO3_12_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 43, "DNO313_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 44, "DNO13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 45, "DNO3_13_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 46, "DNO314_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 47, "DNO14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", 48, "DNO3_14_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "dry_no3_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "wet_no3_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp",  4, "WNO385_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp",  5, "WNO85KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp",  6, "WNO3_85_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp",  7, "WNO386_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp",  8, "WNO86KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp",  9, "WNO3_86_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 10, "WNO387_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 11, "WNO87KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 12, "WNO3_87_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 13, "WNO388_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 14, "WNO88KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 15, "WNO3_88_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 16, "WNO389_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 17, "WNO89KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 18, "WNO3_89_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 19, "WNO390_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 20, "WNO90KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 21, "WNO3_90_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 22, "WNO391_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 23, "WNO91KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 24, "WNO3_91_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 25, "WNO392_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 26, "WNO92KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 27, "WNO3_92_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 28, "WNO393_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 29, "WNO93KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 30, "WNO3_93_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 31, "WNO394_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 32, "WNO94KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 33, "WNO3_94_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 34, "WNO395_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 35, "WNO95KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 36, "WNO3_95_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 37, "WNO396_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 38, "WNO96KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 39, "WNO3_96_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 40, "WNO397_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 41, "WNO97KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 42, "WNO3_97_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 43, "WNO398_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 44, "WNO98KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 45, "WNO3_98_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 46, "WNO399_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 47, "WNO99KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 48, "WNO3_99_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 49, "WNO300_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 50, "WNO00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 51, "WNO3_00_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 52, "WNO301_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 53, "WNO01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 54, "WNO3_01_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 55, "WNO302_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 56, "WNO02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 57, "WNO3_02_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 58, "WNO303_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 59, "WNO03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 60, "WNO3_03_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 61, "WNO304_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 62, "WNO04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 63, "WNO3_04_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 64, "WNO305_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 65, "WNO05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 66, "WNO3_05_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 67, "WNO306_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 68, "WNO06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 69, "WNO3_06_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 70, "WNO307_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 71, "WNO07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 72, "WNO3_07_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 73, "WNO308_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 74, "WNO08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 75, "WNO3_08_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 76, "WNO309_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 77, "WNO09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 78, "WNO3_09_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 79, "WNO310_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 80, "WNO10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 81, "WNO3_10_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 82, "WNO311_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 83, "WNO11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 84, "WNO3_11_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 85, "WNO312_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 86, "WNO12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 87, "WNO3_12_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 88, "WNO313_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 89, "WNO13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 90, "WNO3_13_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 91, "WNO314_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 92, "WNO14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", 93, "WNO3_14_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "wet_no3_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "wet_nh4_source_watershed_nadp",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp",  4, "WNH485_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp",  5, "WNH85KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp",  6, "WNH4_85_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp",  7, "WNH486_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp",  8, "WNH86KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp",  9, "WNH4_86_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 10, "WNH487_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 11, "WNH87KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 12, "WNH4_87_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 13, "WNH488_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 14, "WNH88KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 15, "WNH4_88_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 16, "WNH489_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 17, "WNH89KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 18, "WNH4_89_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 19, "WNH490_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 20, "WNH90KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 21, "WNH4_90_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 22, "WNH491_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 23, "WNH91KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 24, "WNH4_91_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 25, "WNH492_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 26, "WNH92KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 27, "WNH4_92_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 28, "WNH493_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 29, "WNH93KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 30, "WNH4_93_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 31, "WNH494_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 32, "WNH94KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 33, "WNH4_94_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 34, "WNH495_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 35, "WNH95KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 36, "WNH4_95_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 37, "WNH496_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 38, "WNH96KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 39, "WNH4_96_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 40, "WNH497_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 41, "WNH97KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 42, "WNH4_97_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 43, "WNH498_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 44, "WNH98KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 45, "WNH4_98_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 46, "WNH499_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 47, "WNH99KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 48, "WNH4_99_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 49, "WNH400_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 50, "WNH00KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 51, "WNH4_00_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 52, "WNH401_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 53, "WNH01KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 54, "WNH4_01_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 55, "WNH402_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 56, "WNH02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 57, "WNH4_02_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 58, "WNH403_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 59, "WNH03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 60, "WNH4_03_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 61, "WNH404_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 62, "WNH04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 63, "WNH4_04_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 64, "WNH405_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 65, "WNH05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 66, "WNH4_05_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 67, "WNH406_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 68, "WNH06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 69, "WNH4_06_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 70, "WNH407_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 71, "WNH07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 72, "WNH4_07_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 73, "WNH408_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 74, "WNH08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 75, "WNH4_08_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 76, "WNH409_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 77, "WNH09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 78, "WNH4_09_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 79, "WNH410_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 80, "WNH10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 81, "WNH4_10_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 82, "WNH411_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 83, "WNH11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 84, "WNH4_11_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 85, "WNH412_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 86, "WNH12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 87, "WNH4_12_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 88, "WNH413_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 89, "WNH13KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 90, "WNH4_13_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 91, "WNH414_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 92, "WNH14KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", 93, "WNH4_14_%",  FTDouble, 16, 4, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "wet_nh4_source_watershed_nadp",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Source Watershed CMAQ data (original): */

  { "/nitrogen_source_watershed_cmaq",  0, "ESTCODE",     FTString,  6, 0, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 23, "TN_02_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 24, "TN_03_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 25, "TN_04_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 26, "TN_05_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 27, "TN_06_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 28, "TN_07_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 29, "TN_08_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 48, "WN_02_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 49, "WN_03_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 50, "WN_04_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 51, "WN_05_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 52, "WN_06_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 53, "WN_07_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 54, "WN_08_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 30, "TN_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 31, "TN_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 32, "TN_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 33, "TN_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 34, "TN_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 35, "TN_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 36, "TN_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 55, "WN_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 56, "WN_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 57, "WN_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 58, "WN_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 59, "WN_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 60, "WN_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", 61, "WN_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "/nitrogen_source_watershed_cmaq",  3, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* CMAQ watershed data (new replacement version 2016-03-30): */

  { "chloride_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_",  4, "CL_02_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_",  5, "CL_02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_",  6, "CL_03_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_",  7, "CL_03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_",  8, "CL_04_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_",  9, "CL_04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 10, "CL_05_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 11, "CL_05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 12, "CL_06_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 13, "CL_06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 14, "CL_07_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 15, "CL_07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 16, "CL_08_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 17, "CL_08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 18, "CL_09_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 19, "CL_09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 20, "CL_10_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 21, "CL_10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 22, "CL_11_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 23, "CL_11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 24, "CL_12_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 25, "CL_12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 26, "CL_25_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", 27, "CL_25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "chloride_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "_nh3_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_",  4, "NH3_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_",  5, "NH302KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_",  6, "NH3_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_",  7, "NH303KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_",  8, "NH3_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_",  9, "NH304KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 10, "NH3_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 11, "NH305KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 12, "NH3_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 13, "NH306KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 14, "NH3_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 15, "NH307KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 16, "NH3_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 17, "NH308KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 18, "NH3_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 19, "NH309KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 20, "NH3_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 21, "NH310KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 22, "NH3_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 23, "NH311KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 24, "NH3_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 25, "NH312KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 26, "NH3_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", 27, "NH325KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "_nh3_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "_nhx_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_",  4, "NHX_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_",  5, "NHX02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_",  6, "NHX_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_",  7, "NHX03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_",  8, "NHX_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_",  9, "NHX04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 10, "NHX_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 11, "NHX05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 12, "NHX_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 13, "NHX06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 14, "NHX_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 15, "NHX07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 16, "NHX_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 17, "NHX08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 18, "NHX_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 19, "NHX09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 20, "NHX_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 21, "NHX10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 22, "NHX_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 23, "NHX11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 24, "NHX_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 25, "NHX12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 26, "NHX_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", 27, "NHX25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "_nhx_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "_no3_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_",  4, "NO3_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_",  5, "NO302KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_",  6, "NO3_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_",  7, "NO303KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_",  8, "NO3_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_",  9, "NO304KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 10, "NO3_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 11, "NO305KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 12, "NO3_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 13, "NO306KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 14, "NO3_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 15, "NO307KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 16, "NO3_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 17, "NO308KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 18, "NO3_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 19, "NO309KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 20, "NO3_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 21, "NO310KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 22, "NO3_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 23, "NO311KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 24, "NO3_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 25, "NO312KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 26, "NO3_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", 27, "NO325KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "_no3_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "_nox_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_",  4, "NOX_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_",  5, "NOX02KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_",  6, "NOX_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_",  7, "NOX03KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_",  8, "NOX_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_",  9, "NOX04KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 10, "NOX_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 11, "NOX05KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 12, "NOX_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 13, "NOX06KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 14, "NOX_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 15, "NOX07KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 16, "NOX_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 17, "NOX08KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 18, "NOX_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 19, "NOX09KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 20, "NOX_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 21, "NOX10KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 22, "NOX_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 23, "NOX11KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 24, "NOX_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 25, "NOX12KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 26, "NOX_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", 27, "NOX25KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "_nox_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "_nitrogen_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_",  4, "N_02_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_",  5, "N_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_",  6, "N_03_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_",  7, "N_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_",  8, "N_04_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_",  9, "N_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 10, "N_05_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 11, "N_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 12, "N_06_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 13, "N_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 14, "N_07_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 15, "N_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 16, "N_08_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 17, "N_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 18, "N_09_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 19, "N_09_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 20, "N_10_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 21, "N_10_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 22, "N_11_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 23, "N_11_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 24, "N_12_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 25, "N_12_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 26, "N_25_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", 27, "N_25_KGKMY", FTDouble, 16, 4, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0},
  { "_nitrogen_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "_so2_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_",  4, "SO2_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_",  5, "SO202KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_",  6, "SO2_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_",  7, "SO203KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_",  8, "SO2_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_",  9, "SO204KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 10, "SO2_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 11, "SO205KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 12, "SO2_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 13, "SO206KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 14, "SO2_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 15, "SO207KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 16, "SO2_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 17, "SO208KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 18, "SO2_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 19, "SO209KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 20, "SO2_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 21, "SO210KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 22, "SO2_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 23, "SO211KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 24, "SO2_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 25, "SO212KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 26, "SO2_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", 27, "SO225KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "_so2_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "_so4_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_",  4, "SO4_02_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_",  5, "SO402KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_",  6, "SO4_03_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_",  7, "SO403KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_",  8, "SO4_04_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_",  9, "SO404KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 10, "SO4_05_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 11, "SO405KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 12, "SO4_06_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 13, "SO406KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 14, "SO4_07_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 15, "SO407KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 16, "SO4_08_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 17, "SO408KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 18, "SO4_09_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 19, "SO409KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 20, "SO4_10_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 21, "SO410KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 22, "SO4_11_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 23, "SO411KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 24, "SO4_12_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 25, "SO412KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 26, "SO4_25_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", 27, "SO425KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "_so4_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  { "sulfur_source_watershed_cmaq_",  0, "ESTCODE",    FTString,  6, 0, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_",  4, "S_02_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_",  5, "S_02_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_",  6, "S_03_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_",  7, "S_03_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_",  8, "S_04_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_",  9, "S_04_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 10, "S_05_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 11, "S_05_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 12, "S_06_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 13, "S_06_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 14, "S_07_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 15, "S_07_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 16, "S_08_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 17, "S_08_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 18, "S_09_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 19, "S_09_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 20, "S_10_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 21, "S_10_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 22, "S_11_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 23, "S_11_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 24, "S_12_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 25, "S_12_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 26, "S_25_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", 27, "S_25_KGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "sulfur_source_watershed_cmaq_",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },


  /* Nitrogen Source Watershed point data (updated August 28, 2016): */

  { "n_source_watershed_point",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "n_source_watershed_point",  1, "TN_07_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  2, "TN_08_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  3, "TN_09_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  4, "TN_10_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  5, "TN_11_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 19, "TN_12_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 20, "TN_13_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 21, "TN_14_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 22, "TN_15_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  6, "TN_07_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  7, "TN_08_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  8, "TN_09_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point",  9, "TN_10_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 10, "TN_11_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 23, "TN_12_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 24, "TN_13_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 25, "TN_14_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 26, "TN_15_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 11, "TN_07_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 12, "TN_08_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 13, "TN_09_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 14, "TN_10_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 15, "TN_11_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 27, "TN_12_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 28, "TN_13_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 29, "TN_14_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", 30, "TN_15_%", FTDouble, 16, 4, 0.0, 1.0 },
  { "n_source_watershed_point", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "n_source_watershed_point", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "n_source_watershed_point", 18, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Phosphorus Source Watershed point data (delivered August 28, 2016): */

  { "s_source_watershed_point",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0},
  { "s_source_watershed_point",  1, "TP_07_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  2, "TP_08_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  3, "TP_09_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  4, "TP_10_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  5, "TP_11_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 19, "TP_12_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 20, "TP_13_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 21, "TP_14_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 22, "TP_15_KGY",   FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  6, "TP_07_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  7, "TP_08_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  8, "TP_09_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point",  9, "TP_10_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 10, "TP_11_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 23, "TP_12_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 24, "TP_13_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 25, "TP_14_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 26, "TP_15_KGHAY", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 11, "TP_07_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 12, "TP_08_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 13, "TP_09_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 14, "TP_10_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 15, "TP_11_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 27, "TP_12_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 28, "TP_13_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 29, "TP_14_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", 30, "TP_15_%", FTDouble, 16, 4, 0.0, 1.0},
  { "s_source_watershed_point", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0},
  { "s_source_watershed_point", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0},
  { "s_source_watershed_point", 18, "WATER_BODY", FTString, 48, 0, 0.0, 1.0},

  /* Nitrogen Source Watershed non-point data: */

  { "watershed_nonpoint",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "watershed_nonpoint", 20, "CROP_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "watershed_nonpoint", 21, "FERT_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "watershed_nonpoint", 22, "MANU_KGY",   FTDouble, 16, 4, 0.0, 1.0 },
  { "watershed_nonpoint", 23, "CROP_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "watershed_nonpoint", 24, "FERT_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "watershed_nonpoint", 25, "MANU_KGHAY", FTDouble, 16, 4, 0.0, 1.0 },
  { "watershed_nonpoint", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "watershed_nonpoint", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "watershed_nonpoint",  3, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 1992 Atlantic data: */

  { "load_estuary_sparrow_1992_atlantic",
    0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    29, "PREDOM_SRC", FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    -1, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    16, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    28, "TOTNAT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    17, "HUMPOP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    18, "WETDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    27, "TOTFER_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    19, "FERTCS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    20, "FERTAL_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    21, "FERTWT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    22, "FROTHF_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    23, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    24, "FOREST_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    25, "BARREN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    26, "SHRUB_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_atlantic",
     3, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 1992 non-Atlantic (Gulf & Pacific) data: */

  { "load_estuary_sparrow_1992_!atlantic",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    30, "PREDOM_SRC", FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    -1, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    17, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    29, "TOTNAT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    18, "HUMPOP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    19, "WETDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    28, "TOTFER_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    20, "FERTCS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    21, "FERTAL_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    22, "FERTWT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    23, "FROTHF_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    24, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    25, "FOREST_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    26, "BARREN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    27, "SHRUB_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_1992_!atlantic",
     3, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2002 MRB1 data: */

  { "load_estuary_sparrow_2002_mrb1",
    0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    -1, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    4, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    9, "FERTCS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    5, "FERTOT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    7, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    8, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    6, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    10, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb1",
    1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2002 MRB2 data: */

  { "load_estuary_sparrow_2002_mrb2",
    0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    -1, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    8, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    10, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    9, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    12, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    11, "URBANR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    13, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb2",
    1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2002 MRB5 data: */

  { "load_estuary_sparrow_2002_mrb5",
    0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    -1, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    4, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    5, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    7, "MANURF_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    6, "MANURP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    10, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    8, "URBANR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    9, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb5",
    1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2002 MRB7 data: */

  { "load_estuary_sparrow_2002_mrb7",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    -1, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    4, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    10, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    9, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    12, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    5, "FORALD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    7, "FORWES_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    6, "FOREAS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    8, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    11, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    13, "CANADA_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    14, "BOUNDS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2002_mrb7",
    1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2011 MRB1+MRB2 data: */

  { "load_estuary_sparrow_2011_mrb1mrb2",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
    14, "DOM_LOAD",   FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     2, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     4, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     5, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     6, "FERMNR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     7, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
    13, "ESATDP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     8, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     9, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
    12, "WSDAREAKM2", FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb1mrb2",
     1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2011 MRB5 data: */

  { "load_estuary_sparrow_2011_mrb5",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
    16, "DOM_LOAD",   FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     2, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     4, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     5, "MANURF_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     6, "MANURP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     7, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     8, "FERMNR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     9, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
    15, "ESATDP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
    10, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
    11, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
    14, "WSDAREAKM2", FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb5",
     1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2011 MRB7 data: */

  { "load_estuary_sparrow_2011_mrb7",
    0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    18, "DOM_LOAD",  FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    2, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    4, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    5, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    6, "FERMNR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    7, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    17, "ESATDP_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    8, "FORALD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    9, "FOR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    10, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    11, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    12, "NSEWER_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    13, "SPRPWR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    16, "WSDAREAKM2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb7",
    1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Estuary SPARROW 2011 MRB8 data: */

  { "load_estuary_sparrow_2011_mrb8",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
    15, "DOM_LOAD",   FTString,  4, 0, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     2, "TOT_YKGKMY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     4, "FERTCM_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     5, "MANURP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     6, "FERMNR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     7, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
    14, "ESATDP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     8, "FOR_KGY",    FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     9, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
    10, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
    13, "WSDAREAKM2", FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_estuary_sparrow_2011_mrb8",
     1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Subestuary SPARROW 1992 data: */

  { "load_subestuary_sparrow_1992_",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     1, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    16, "PREDOM_SRC", FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    15, "TOTNAT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     4, "HUMPOP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     5, "WETDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    14, "TOTFER_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     6, "FERTCS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     7, "FERTAL_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     8, "FERTWT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     9, "FROTHF_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    10, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    11, "FOREST_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    12, "BARREN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    13, "SHRUB_KGY",  FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     2, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_1992_",
     3, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Subestuary SPARROW 2002 MRB1 data: */

  { "load_subestuary_sparrow_2002_mrb1_",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     1, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     8, "FERTCS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     4, "FERTOT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     6, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     7, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     5, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     9, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb1_",
     2, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Subestuary SPARROW 2002 MRB2 data: */

  { "load_subestuary_sparrow_2002_mrb2_",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     1, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     5, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     4, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     7, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     6, "URBANR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     8, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb2_",
     2, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Subestuary SPARROW 2002 MRB5 data: */

  { "load_subestuary_sparrow_2002_mrb5_",
    0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    1, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    4, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    6, "MANURF_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    5, "MANURP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    9, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    7, "URBANR_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    8, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb5_",
    2, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },

  /* Nitrogen Load Subestuary SPARROW 2002 MRB7 data: */

  { "load_subestuary_sparrow_2002_mrb7_",
     0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     1, "SUBCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     3, "TOT_LD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     9, "FERTAG_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     8, "MANURE_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
    11, "ATMDEP_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     4, "FORALD_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     6, "FORWES_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     5, "FOREAS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     7, "DEVLAN_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
    10, "MUNIPT_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
    12, "CANADA_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
    13, "BOUNDS_KGY", FTDouble, 16, 4, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
    -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
    -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "load_subestuary_sparrow_2002_mrb7_",
     2, "SUBEMBAYMT", FTString, 48, 0, 0.0, 1.0 },

  /* Coastal Vulnerability Index data: */

  { "coastal_vulnerability_atl", -1, "LENGTH_KM",  FTDouble, 10,4,0.0,1.0},
  { "coastal_vulnerability_atl", 22, "CVI",        FTDouble, 10,4,0.0,1.0},
  { "coastal_vulnerability_atl", 23, "CVI_RANK",   FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_atl", 30, "CVI_RISK",   FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_atl", 20, "MEANWAVE_M", FTDouble, 10,2,0.0,1.0},
  { "coastal_vulnerability_atl", 21, "WAVE_RANK",  FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_atl", 25, "WAVE_RISK",  FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_atl", 11, "MEANTIDE_M", FTDouble, 10,2,0.0,1.0},
  { "coastal_vulnerability_atl", 12, "TIDE_RANK",  FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_atl", 24, "TIDE_RISK",  FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_atl", 17, "SLRISE_MMY", FTDouble, 10,1,0.0,1.0},
  { "coastal_vulnerability_atl", 18, "SL_RANK",    FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_atl", 27, "SL_RISK",    FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_atl", 13, "SLOPE_PCNT", FTDouble, 10,4,0.0,1.0},
  { "coastal_vulnerability_atl", 14, "SLOPE_RANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_atl", 29, "SLOPE_RISK", FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_atl", 15, "EROACC_MYR", FTDouble, 10,2,0.0,1.0},
  { "coastal_vulnerability_atl", 16, "EROACCRANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_atl", 26, "EROACCRISK", FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_atl", 19, "GEOMO_RANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_atl", 28, "GEOMO_RISK", FTString, 9,0,0.0,1.0},

  { "coastal_vulnerability_gulf", -1, "LENGTH_KM",  FTDouble,10,4,0.0,1.0},
  { "coastal_vulnerability_gulf", 22, "CVI",        FTDouble,10,4,0.0,1.0},
  { "coastal_vulnerability_gulf", 23, "CVI_RANK",   FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 30, "CVI_RISK",   FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 20, "MEANWAVE_M", FTDouble,10,2,0.0,1.0},
  { "coastal_vulnerability_gulf", 21, "WAVE_RANK",  FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 25, "WAVE_RISK",  FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 11, "MEANTIDE_M", FTDouble,10,2,0.0,1.0},
  { "coastal_vulnerability_gulf", 12, "TIDE_RANK",  FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 24, "TIDE_RISK",  FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 17, "SLRISE_MMY", FTDouble,10,1,0.0,1.0},
  { "coastal_vulnerability_gulf", 18, "SL_RANK",    FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 28, "SL_RISK",    FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 13, "SLOPE_PCNT", FTDouble,10,4,0.0,1.0},
  { "coastal_vulnerability_gulf", 14, "SLOPE_RANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 26, "SLOPE_RISK", FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 15, "EROACC_MYR", FTDouble,10,2,0.0,1.0},
  { "coastal_vulnerability_gulf", 16, "EROACCRANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 27, "EROACCRISK", FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 19, "GEOMO_RANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_gulf", 29, "GEOMO_RISK", FTString, 9,0,0.0,1.0},

  { "coastal_vulnerability_pac", -1, "LENGTH_KM",  FTDouble, 10,4,0.0,1.0},
  { "coastal_vulnerability_pac", 25, "CVI",        FTDouble, 10,4,0.0,1.0},
  { "coastal_vulnerability_pac", 26, "CVI_RANK",   FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_pac", 27, "CVI_RISK",   FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_pac", 17, "MEANWAVE_M", FTDouble, 10,2,0.0,1.0},
  { "coastal_vulnerability_pac", 18, "WAVE_RANK",  FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_pac", 19, "WAVE_RISK",  FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_pac", 14, "MEANTIDE_M", FTDouble, 10,2,0.0,1.0},
  { "coastal_vulnerability_pac", 15, "TIDE_RANK",  FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_pac", 16, "TIDE_RISK",  FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_pac",  8, "SLRISE_MMY", FTDouble, 10,1,0.0,1.0},
  { "coastal_vulnerability_pac",  9, "SL_RANK",    FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_pac", 10, "SL_RISK",    FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_pac", 11, "SLOPE_PCNT", FTDouble, 10,4,0.0,1.0},
  { "coastal_vulnerability_pac", 12, "SLOPE_RANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_pac", 13, "SLOPE_RISK", FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_pac", 20, "EROACC_MYR", FTDouble, 10,2,0.0,1.0},
  { "coastal_vulnerability_pac", 21, "EROACCRANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_pac", 22, "EROACCRISK", FTString, 9,0,0.0,1.0},
  { "coastal_vulnerability_pac", 23, "GEOMO_RANK", FTInteger,1,0,0.0,1.0},
  { "coastal_vulnerability_pac", 24, "GEOMO_RISK", FTString, 9,0,0.0,1.0},

  /* Estuary nutrient sensitivity: */

  { "sensitivity",  0, "ESTCODE",    FTString,  4, 0, 0.0, 1.0 },
  { "sensitivity",  5, "ESTRY_SQKM", FTDouble, 11, 3, 0.0, 1.0 },
  { "sensitivity",  6, "MIX_SQKM",   FTDouble, 11, 3, 0.0, 1.0 },
  { "sensitivity",  7, "SEA_SQKM",   FTDouble, 11, 3, 0.0, 1.0 },
  { "sensitivity",  8, "FRESH_SQKM", FTDouble, 11, 3, 0.0, 1.0 },
  { "sensitivity",  9, "AVFL_M3/DY", FTDouble, 14, 3, 0.0, 1.0 },
  { "sensitivity", 10, "MXFL_M3/DY", FTDouble, 14, 3, 0.0, 1.0 },
  { "sensitivity", 11, "ESTVOL_BM3", FTDouble,  6, 3, 0.0, 1.0 },
  { "sensitivity", 12, "TIDPRM_BM3", FTDouble,  6, 3, 0.0, 1.0 },
  { "sensitivity", 13, "TIDE_HT_M",  FTDouble,  6, 2, 0.0, 1.0 },
  { "sensitivity", 14, "BOTSAL_PPT", FTDouble,  6, 2, 0.0, 1.0 },
  { "sensitivity", 15, "TOPSAL_PPT", FTDouble,  6, 2, 0.0, 1.0 },
  { "sensitivity", 16, "DEPTH_M",    FTDouble,  6, 2, 0.0, 1.0 },
  { "sensitivity", 17, "DCP_MG/L",   FTDouble, 10, 2, 0.0, 1.0 },
  { "sensitivity", 18, "PRE_DAYS",   FTDouble, 10, 2, 0.0, 1.0 },
  { "sensitivity", 20, "AREA_M2",    FTDouble, 12, 1, 0.0, 1.0 },
  { "sensitivity", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "sensitivity", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },
  { "sensitivity",  3, "EDA",        FTString,  5, 0, 0.0, 1.0 },
  { "sensitivity",  2, "EDA_CDA",    FTString,  5, 0, 0.0, 1.0 },
  { "sensitivity", 19, "EDA_SUBEDA", FTString,  5, 0, 0.0, 1.0 },
  { "sensitivity",  4, "EDA_NAME",   FTString, 32, 0, 0.0, 1.0 },
  { "sensitivity",  1, "WATER_BODY", FTString, 48, 0, 0.0, 1.0 },

  /* Ground water contact in catchment & watershed (added 2017-09-29): */

  { "contact_catchment_",  2, "MEAN_DAYS",  FTDouble,  11, 3, 0.0, 1.0 },
  { "contact_catchment_",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "contact_catchment_",  1, "FIPS_SQKM",  FTDouble,  11, 3, 0.0, 1.0 },
  { "contact_catchment_", -1, "AREA_SQKM",  FTDouble,  11, 3, 0.0, 1.0 },
  { "contact_catchment_", -1, "HECTARES",   FTDouble,  20, 6, 0.0, 1.0 },

  { "contact_watershed_",  0, "ESTCODE",    FTString,   4, 0, 0.0, 1.0 },
  { "contact_watershed_",  3, "MEAN_DAYS",  FTDouble,  11, 3, 0.0, 1.0 },
  { "contact_watershed_",  2, "FIPS_SQKM",  FTDouble,  11, 3, 0.0, 1.0 },
  { "contact_watershed_", -1, "AREA_SQKM",  FTDouble,  11, 3, 0.0, 1.0 },
  { "contact_watershed_", -1, "HECTARES",   FTDouble,  20, 6, 0.0, 1.0 },
  { "contact_watershed_",  1, "WATER_BODY", FTString,  48, 0, 0.0, 1.0 },

  /* Estuary nutrient sensitivity volume 2 2020-11-17: */

  /* stream discharge: */

  { "stream_discharge_usgs_", 14, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "stream_discharge_usgs_",  3, "SOURCE_FEA", FTString,  16, 0, 0.0, 1.0 },
  { "stream_discharge_usgs_",  6, "FLCOMID",    FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_discharge_usgs_",  9, "STATE",      FTString,   2, 0, 0.0, 1.0 },
  { "stream_discharge_usgs_", 13, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "stream_discharge_usgs_", 12, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "stream_discharge_usgs_",  0, "REACH_CODE", FTString,  14, 0, 0.0, 1.0 },
  { "stream_discharge_usgs_", 11, "DRAIN_SQKM", FTDouble,  12, 4, 0.0, 1.0 },
  { "stream_discharge_usgs_", 10, "DRAIN_SQMI", FTDouble,  12, 4, 0.0, 1.0 },
  { "stream_discharge_usgs_", 15, "GAUGE_ID",   FTString,  16, 0, 0.0, 1.0 },
  { "stream_discharge_usgs_",  4, "GAUGE_URL",  FTString,  64, 0, 0.0, 1.0 },
  { "stream_discharge_usgs_",  8, "LOCATION",   FTString,  64, 0, 0.0, 1.0 },

  { "stream_discharge_nhd_",  73, "ESTCODE",    FTString,   5, 0, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",   0, "COMID",      FTInteger, 10, 0, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  33, "FLOW_MA",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  36, "FLOW_01",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  39, "FLOW_02",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  42, "FLOW_03",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  45, "FLOW_04",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  48, "FLOW_05",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  51, "FLOW_06",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  54, "FLOW_07",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  57, "FLOW_08",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  60, "FLOW_09",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  63, "FLOW_10",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  66, "FLOW_11",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  69, "FLOW_12",    FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  31, "RUNOFF_MA",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  34, "RUNOFF_01",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  37, "RUNOFF_02",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  40, "RUNOFF_03",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  43, "RUNOFF_04",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  46, "RUNOFF_05",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  49, "RUNOFF_06",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  52, "RUNOFF_07",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  55, "RUNOFF_08",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  58, "RUNOFF_09",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  61, "RUNOFF_10",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  64, "RUNOFF_11",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  67, "RUNOFF_12",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  32, "ADJUST_MA",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  35, "ADJUST_01",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  38, "ADJUST_02",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  41, "ADJUST_03",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  44, "ADJUST_04",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  47, "ADJUST_05",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  50, "ADJUST_06",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  53, "ADJUST_07",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  56, "ADJUST_08",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  59, "ADJUST_09",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  62, "ADJUST_10",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  65, "ADJUST_11",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",  68, "ADJUST_12",  FTDouble,  10, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "stream_discharge_nhd_",   5, "REACH_CODE", FTString,  14, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",   8, "FCODE",      FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  28, "TIDAL",      FTInteger,  1, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  10, "STREAM_LEV", FTInteger,  1, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  11, "STREAM_ORD", FTInteger,  1, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  26, "DRAIN_SQKM", FTDouble,  10, 4, 0.0, 1.0 },
  { "stream_discharge_nhd_",  27, "DIVDA_SQKM", FTDouble,  10, 4, 0.0, 1.0 },
  { "stream_discharge_nhd_",  13, "FROM_NODE",  FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  14, "TO_NODE",    FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  16, "PATH_ID",    FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  18, "TER_PATHID", FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",   1, "FDATE",      FTString,  10, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",  72, "GAUGE_ID",   FTString,  10, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",   2, "GNIS_ID",    FTString,   8, 0, 0.0, 1.0 },
  { "stream_discharge_nhd_",   3, "GNIS_NAME",  FTString,  64, 0, 0.0, 1.0 },

  /* tide_point: */

  { "tide_point_",  1, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "tide_point_",  0, "SITE",       FTString,   8, 0, 0.0, 1.0 },
  { "tide_point_",  7, "AVG_MLW",    FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_",  6, "AVG_MHW",    FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_",  9, "AVG_LLW",    FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_",  8, "AVG_HHW",    FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_", 11, "AVG_MINTR",  FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_", 10, "AVG_MAXTR",  FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_", 17, "AVGMEANTPV", FTDouble,  16, 2, 0.0, 1.0 },
  { "tide_point_", 19, "AVG_MINTPV", FTDouble,  16, 2, 0.0, 1.0 },
  { "tide_point_", 18, "AVG_MAXTPV", FTDouble,  16, 2, 0.0, 1.0 },
  { "tide_point_", 20, "AVG_MTP",    FTDouble,   8, 4, 0.0, 1.0 },
  { "tide_point_",  4, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "tide_point_",  3, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "tide_point_", 16, "LTIDE_CONV", FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_", 15, "HTIDE_CONV", FTDouble,  12, 4, 0.0, 1.0 },
  { "tide_point_", 12, "DEM_VDATUM", FTString,   8, 0, 0.0, 1.0 },
  { "tide_point_", 13, "TIDEVDATUM", FTString,   8, 0, 0.0, 1.0 },
  { "tide_point_",  5, "METHOD",     FTString,  16, 0, 0.0, 1.0 },
  { "tide_point_", 14, "TIDEADJUST", FTString,  32, 0, 0.0, 1.0 },
  { "tide_point_",  2, "LOCATION",   FTString,  64, 0, 0.0, 1.0 },

  /* tide_current: */

  { "tide_current_",  0, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "tide_current_",  1, "SITE",       FTString,   8, 0, 0.0, 1.0 },
  { "tide_current_",  5, "AVGCURRENT", FTDouble,  10, 5, 0.0, 1.0 },
  { "tide_current_",  2, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "tide_current_",  3, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "tide_current_",  4, "METHOD",     FTString,  16, 0, 0.0, 1.0 },

  /* longshore_current: */

  { "longshore_current_",  0, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "longshore_current_",  7, "AVGCURRENT", FTDouble,   8, 2, 0.0, 1.0 },
  { "longshore_current_",  3, "AVG_UMAX",   FTDouble,   8, 2, 0.0, 1.0 },
  { "longshore_current_",  4, "AVG_VMAX",   FTDouble,   8, 2, 0.0, 1.0 },
  { "longshore_current_",  5, "SHORE_DEGN", FTDouble,   8, 3, 0.0, 1.0 },
  { "longshore_current_",  6, "SHORE_RADN", FTDouble,   8, 3, 0.0, 1.0 },
  { "longshore_current_",  1, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "longshore_current_",  2, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },

  /* estuary_flushing: */

  { "estuary_flushing_",  0, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "estuary_flushing_", 15, "TPM_FT_AVG", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 18, "TPM_FT_MED", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 16, "TPM_FT_MIN", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 17, "TPM_FT_MAX", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 22, "FFM_FT_AVG", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 19, "FFM_FT_MED", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 20, "FFM_FT_MIN", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 21, "FFM_FT_MAX", FTDouble,  16, 5, 0.0, 1.0 },
  { "estuary_flushing_", 12, "AVG_TPV",    FTDouble,  16, 3, 0.0, 1.0 },
  { "estuary_flushing_",  5, "AVG_VOLUME", FTDouble,  16, 3, 0.0, 1.0 },
  { "estuary_flushing_",  3, "AVG_AREA",   FTDouble,  16, 2, 0.0, 1.0 },
  { "estuary_flushing_",  6, "AVG_DEPTH",  FTDouble,   8, 2, 0.0, 1.0 },
  { "estuary_flushing_",  4, "MOUTHWIDTH", FTDouble,   8, 2, 0.0, 1.0 },
  { "estuary_flushing_",  7, "AVG_MTP",    FTDouble,   8, 4, 0.0, 1.0 },
  { "estuary_flushing_",  8, "AVG_SAL",    FTDouble,   8, 3, 0.0, 1.0 },
  { "estuary_flushing_",  9, "AVG_FLOW",   FTDouble,  10, 3, 0.0, 1.0 },
  { "estuary_flushing_", 14, "AVG_PRECIP", FTDouble,  10, 3, 0.0, 1.0 },
  { "estuary_flushing_", 10, "STRAT_TYPE", FTString,  16, 0, 0.0, 1.0 },
  { "estuary_flushing_", 11, "STRAT_METH", FTString,   8, 0, 0.0, 1.0 },
  { "estuary_flushing_", 13, "ECO_REGION", FTString,   8, 0, 0.0, 1.0 },
  { "estuary_flushing_",  2, "STATE",      FTString,   2, 0, 0.0, 1.0 },
  { "estuary_flushing_", 23, "FLOW_TYPE",  FTString,  16, 0, 0.0, 1.0 },

  /* flowlines_: */

  /* 2021-08-05 RBEROST dbf file is preprocessed once into flowlines_ columns:*/

  { "rberost",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "rberost",  6, "TOTAL_N",    FTDouble,  16, 5, 0.0, 1.0 },
  { "rberost",  7, "INCREM_N",   FTDouble,  16, 5, 0.0, 1.0 },
  { "rberost",  8, "TOTAL_P",    FTDouble,  16, 5, 0.0, 1.0 },
  { "rberost",  9, "INCREM_P",   FTDouble,  16, 5, 0.0, 1.0 },
  { "rberost",  2, "LENGTH_KM",  FTDouble,   8, 3, 0.0, 1.0 },
  { "rberost",  3, "FROM_NODE",  FTInteger, 10, 0, 0.0, 1.0 },
  { "rberost",  4, "TO_NODE",    FTInteger, 10, 0, 0.0, 1.0 },
  { "rberost",  5, "HYDRO_SEQ",  FTInteger, 10, 0, 0.0, 1.0 },
  { "rberost",  1, "GNIS_NAME",  FTString,  48, 0, 0.0, 1.0 },

  /* These flowlines_ columns result from the above preprocessing: */

  { "flowlines_upper_ct",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_upper_ct",  1, "TOTAL_N",    FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_upper_ct",  2, "INCREM_N",   FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_upper_ct",  3, "TOTAL_P",    FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_upper_ct",  4, "INCREM_P",   FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_upper_ct",  5, "LENGTH_KM",  FTDouble,   8, 3, 0.0, 1.0 },
  { "flowlines_upper_ct",  6, "FROM_NODE",  FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_upper_ct",  7, "TO_NODE",    FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_upper_ct",  8, "HYDRO_SEQ",  FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_upper_ct",  9, "GNIS_NAME",  FTString,  48, 0, 0.0, 1.0 },

  /* Added 2024-07-02 includes watershed in last column: */

  { "flowlines_puget_sound",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_puget_sound",  1, "TOTAL_N",    FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_puget_sound",  2, "INCREM_N",   FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_puget_sound",  3, "TOTAL_P",    FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_puget_sound",  4, "INCREM_P",   FTDouble,  16, 5, 0.0, 1.0 },
  { "flowlines_puget_sound",  5, "LENGTH_KM",  FTDouble,   8, 3, 0.0, 1.0 },
  { "flowlines_puget_sound",  6, "FROM_NODE",  FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_puget_sound",  7, "TO_NODE",    FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_puget_sound",  8, "HYDRO_SEQ",  FTInteger, 10, 0, 0.0, 1.0 },
  { "flowlines_puget_sound",  9, "GNIS_NAME",  FTString,  48, 0, 0.0, 1.0 },
  { "flowlines_puget_sound", 10, "WATERSHED",  FTString,  48, 0, 0.0, 1.0 },

  /* discharge_points_federal_puget_sound_watershed.dbf (points) */

  { "_federal_puget_sound_", 14, "COMID",       FTInteger, 14, 0, 0.0, 1.0 },
  { "_federal_puget_sound_", 10, "YEAR",        FTInteger,  4, 0, 0.0, 1.0 },
  { "_federal_puget_sound_", 11, "MONTH",       FTInteger,  2, 0, 0.0, 1.0 },
  { "_federal_puget_sound_",  4, "LONGITUDE",   FTDouble,  11, 6, 0.0, 1.0 },
  { "_federal_puget_sound_",  3, "LATITUDE",    FTDouble,  10, 6, 0.0, 1.0 },
  { "_federal_puget_sound_",  5, "FLOW_MG_DY",  FTDouble,  18, 8, 0.0, 1.0 },
  { "_federal_puget_sound_",  5, "FLOW_ML_DY",  FTDouble,  18, 8, 0.0,
                                                           GALLONS_TO_LITERS },
  { "_federal_puget_sound_",  5, "FLOW_M3_DY",  FTDouble,  18, 8, 0.0,
                                  GALLONS_TO_LITERS * LITERS_PER_CUBIC_METER },
  { "_federal_puget_sound_",  6, "TN_MG_L",     FTDouble,  18, 8, 0.0, 1.0 },
  { "_federal_puget_sound_",  7, "TP_MG_L",     FTDouble,  18, 8, 0.0, 1.0 },
  { "_federal_puget_sound_",  8, "TN_KG_DAY",   FTDouble,  18, 8, 0.0, 1.0 },
  { "_federal_puget_sound_",  9, "TP_KG_DAY",   FTDouble,  18, 8, 0.0, 1.0 },
  { "_federal_puget_sound_", 12, "TN_KG_MONT",  FTDouble,  18, 8, 0.0, 1.0 },
  { "_federal_puget_sound_", 13, "TP_KG_MONT",  FTDouble,  18, 8, 0.0, 1.0 },
  { "_federal_puget_sound_",  0, "NPDES_ID",    FTString,  10, 0, 0.0, 1.0 },
  { "_federal_puget_sound_",  1, "PERMITTEE",   FTString,  70, 0, 0.0, 1.0 },
  { "_federal_puget_sound_",  2, "WATERBODY",   FTString,  50, 0, 0.0, 1.0 },

  /* discharge_points_state_puget_sound_watershed.dbf (points) */

  { "_state_puget_sound_", 17, "COMID",       FTInteger, 14, 0, 0.0, 1.0 },
  { "_state_puget_sound_", 12, "YEAR",        FTInteger,  4, 0, 0.0, 1.0 },
  { "_state_puget_sound_", 13, "MONTH",       FTInteger,  2, 0, 0.0, 1.0 },
  { "_state_puget_sound_", 19, "LONGITUDE",   FTDouble,  11, 6, 0.0, 1.0 },
  { "_state_puget_sound_", 18, "LATITUDE",    FTDouble,  10, 6, 0.0, 1.0 },
  { "_state_puget_sound_",  1, "FLOW_MG_DY",  FTDouble,  18, 8, 0.0, 1.0 },
  { "_state_puget_sound_",  1, "FLOW_ML_DY",  FTDouble,  18, 8, 0.0,
                                                           GALLONS_TO_LITERS },
  { "_state_puget_sound_",  1, "FLOW_M3_DY",  FTDouble,  18, 8, 0.0,
                                  GALLONS_TO_LITERS * LITERS_PER_CUBIC_METER },
  { "_state_puget_sound_",  2, "TN_MG_L",     FTDouble,  18, 8, 0.0, 1.0 },
  { "_state_puget_sound_",  3, "TP_MG_L",     FTDouble,  18, 8, 0.0, 1.0 },
  { "_state_puget_sound_", 10, "TN_KG_DAY",   FTDouble,  18, 8, 0.0, 1.0 },
  { "_state_puget_sound_", 11, "TP_KG_DAY",   FTDouble,  18, 8, 0.0, 1.0 },
  { "_state_puget_sound_", 15, "TN_KG_MONT",  FTDouble,  18, 8, 0.0, 1.0 },
  { "_state_puget_sound_", 16, "TP_KG_MONT",  FTDouble,  18, 8, 0.0, 1.0 },
  { "_state_puget_sound_",  0, "NPDES_ID",    FTString,  20, 0, 0.0, 1.0 },
  { "_state_puget_sound_",  4, "PERMITTEE",   FTString,  70, 0, 0.0, 1.0 },
  { "_state_puget_sound_",  6, "WATERBODY",   FTString,  50, 0, 0.0, 1.0 },
  { "_state_puget_sound_",  5, "PERMITTEEW",  FTString,  70, 0, 0.0, 1.0 },
  { "_state_puget_sound_",  7, "NOTE",        FTString,  80, 0, 0.0, 1.0 },

  /* salinity_point: */

  { "salinity_point_atlantic", 11, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "salinity_point_atlantic",  0, "SITE_ID",    FTString,  20, 0, 0.0, 1.0 },
  { "salinity_point_atlantic", 10, "SAL_AVG",    FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_atlantic",  2, "SAL_DIFF",   FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_atlantic", 12, "SAL_OCEAN",  FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_atlantic", 13, "FRESH_FRAC", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_atlantic",  1, "GAUGEDEPTH", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_atlantic",  3, "STRATIFY_P", FTDouble,  12, 9, 0.0, 1.0 },
  { "salinity_point_atlantic",  4, "YEAR",       FTInteger,  4, 0, 0.0, 1.0 },
  { "salinity_point_atlantic",  5, "MONTH",      FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_atlantic",  6, "DAY",        FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_atlantic",  8, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "salinity_point_atlantic",  7, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },

  { "salinity_point_lower_mi", 11, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "salinity_point_lower_mi",  0, "SITE_ID",    FTString,  20, 0, 0.0, 1.0 },
  { "salinity_point_lower_mi", 10, "SAL_AVG",    FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_lower_mi",  2, "SAL_DIFF",   FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_lower_mi", 12, "SAL_OCEAN",  FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_lower_mi", 13, "FRESH_FRAC", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_lower_mi",  1, "GAUGEDEPTH", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_lower_mi",  3, "STRATIFY_P", FTDouble,  12, 9, 0.0, 1.0 },
  { "salinity_point_lower_mi",  4, "YEAR",       FTInteger,  4, 0, 0.0, 1.0 },
  { "salinity_point_lower_mi",  5, "MONTH",      FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_lower_mi",  6, "DAY",        FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_lower_mi",  8, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "salinity_point_lower_mi",  7, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },

  { "salinity_point_gulf", 11, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "salinity_point_gulf",  0, "SITE_ID",    FTString,  20, 0, 0.0, 1.0 },
  { "salinity_point_gulf", 10, "SAL_AVG",    FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_gulf",  2, "SAL_DIFF",   FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_gulf", 12, "SAL_OCEAN",  FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_gulf", 16, "FRESH_FRAC", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_gulf", 13, "FRESH_FLOW", FTDouble,  12, 6, 0.0, 1.0 },
  { "salinity_point_gulf",  1, "GAUGEDEPTH", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_gulf",  3, "STRATIFY_P", FTDouble,  12, 9, 0.0, 1.0 },
  { "salinity_point_gulf",  4, "YEAR",       FTInteger,  4, 0, 0.0, 1.0 },
  { "salinity_point_gulf",  5, "MONTH",      FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_gulf",  6, "DAY",        FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_gulf",  8, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "salinity_point_gulf",  7, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "salinity_point_gulf", 14, "STRAT_TYPE", FTString,  10, 0, 0.0, 1.0 },
  { "salinity_point_gulf", 15, "STRAT_METH", FTString,   8, 0, 0.0, 1.0 },

  { "salinity_point_pacific", 11, "ESTCODE",    FTString,   5, 0, 0.0, 1.0 },
  { "salinity_point_pacific",  0, "SITE_ID",    FTString,  20, 0, 0.0, 1.0 },
  { "salinity_point_pacific", 10, "SAL_AVG",    FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_pacific",  2, "SAL_DIFF",   FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_pacific", 12, "SAL_OCEAN",  FTDouble,  12, 8, 0.0, 1.0 },
  { "salinity_point_pacific", 16, "FRESH_FRAC", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_pacific", 13, "FRESH_FLOW", FTDouble,  12, 6, 0.0, 1.0 },
  { "salinity_point_pacific",  1, "GAUGEDEPTH", FTDouble,   8, 3, 0.0, 1.0 },
  { "salinity_point_pacific",  3, "STRATIFY_P", FTDouble,  12, 9, 0.0, 1.0 },
  { "salinity_point_pacific",  4, "YEAR",       FTInteger,  4, 0, 0.0, 1.0 },
  { "salinity_point_pacific",  5, "MONTH",      FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_pacific",  6, "DAY",        FTInteger,  2, 0, 0.0, 1.0 },
  { "salinity_point_pacific",  8, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },
  { "salinity_point_pacific",  7, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "salinity_point_pacific", 14, "STRAT_TYPE", FTString,  10, 0, 0.0, 1.0 },
  { "salinity_point_pacific", 15, "STRAT_METH", FTString,   8, 0, 0.0, 1.0 },

  /* sediment_nca original files (delete after public deployment 2021-08-20) : */

  { "sediment_nca_atlantic", 23, "STATION",    FTString, 16, 0, 0.0, 1.0 },
  { "sediment_nca_atlantic", 24, "DATE",       FTInteger, 8, 0, 0.0, 1.0 },
  { "sediment_nca_atlantic",  5, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "sediment_nca_atlantic",  6, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "sediment_nca_atlantic", 21, "TOC_%",      FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 32, "CLAY_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 32, "SILTCLAY_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 29, "SILT_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 31, "SAND_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 36, "MOISTURE_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 27, "25th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 34, "50th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 33, "75th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 35, "DEVIATION",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic", 28, "SKEWNESS",   FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_atlantic",  0, "AGENCY",     FTString, 16, 0, 0.0, 1.0 },

  { "sediment_nca_gulf", 23, "STATION",    FTString, 16, 0, 0.0, 1.0 },
  { "sediment_nca_gulf", 24, "DATE",       FTInteger, 8, 0, 0.0, 1.0 },
  { "sediment_nca_gulf",  5, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "sediment_nca_gulf",  6, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "sediment_nca_gulf", 21, "TOC_%",      FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 32, "CLAY_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 32, "SILTCLAY_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 29, "SILT_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 31, "SAND_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 36, "MOISTURE_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 27, "25th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 34, "50th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 33, "75th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 35, "DEVIATION",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf", 28, "SKEWNESS",   FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_gulf",  0, "AGENCY",     FTString, 16, 0, 0.0, 1.0 },

  { "sediment_nca_pacific", 23, "STATION",    FTString, 16, 0, 0.0, 1.0 },
  { "sediment_nca_pacific", 24, "DATE",       FTInteger, 8, 0, 0.0, 1.0 },
  { "sediment_nca_pacific",  5, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "sediment_nca_pacific",  6, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "sediment_nca_pacific", 21, "TOC_%",      FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 32, "CLAY_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 32, "SILTCLAY_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 29, "SILT_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 31, "SAND_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 36, "MOISTURE_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 27, "25th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 34, "50th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 33, "75th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 35, "DEVIATION",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific", 28, "SKEWNESS",   FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_pacific",  0, "AGENCY",     FTString, 16, 0, 0.0, 1.0 },

  /* 2021-08-20 sediment_nca_1990-2006 renamed copy of original files: */

  { "sediment_nca_1990-2006", 23, "STATION",    FTString, 16, 0, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 24, "DATE",       FTInteger, 8, 0, 0.0, 1.0 },
  { "sediment_nca_1990-2006",  5, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "sediment_nca_1990-2006",  6, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 21, "TOC_%",      FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 32, "CLAY_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 32, "SILTCLAY_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 29, "SILT_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 31, "SAND_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 36, "MOISTURE_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 27, "25th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 34, "50th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 33, "75th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 35, "DEVIATION",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006", 28, "SKEWNESS",   FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_1990-2006",  0, "AGENCY",     FTString, 16, 0, 0.0, 1.0 },

  /* 2021-08-20 sediment_nca_2015: */

  { "sediment_nca_2015",  0, "STATION",    FTString, 16, 0, 0.0, 1.0 },
  { "sediment_nca_2015",  1, "DATE",       FTInteger, 8, 0, 0.0, 1.0 },
  { "sediment_nca_2015",  2, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "sediment_nca_2015",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "sediment_nca_2015",  4, "TOC_%",      FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015",  5, "CLAY_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015",  6, "SILTCLAY_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015",  7, "SILT_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015",  8, "SAND_%",     FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015",  9, "MOISTURE_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015", 10, "25th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015", 11, "50th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015", 12, "75th%_PHI",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015", 13, "DEVIATION",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015", 14, "SKEWNESS",   FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_nca_2015", 15, "AGENCY",     FTString, 16, 0, 0.0, 1.0 },
  { "sediment_nca_2015", 16, "ESTCODE",    FTString,  5, 0, 0.0, 1.0 },

  /* All other sediment NOT nca (calculated, extracted, kriged, parsed): */

  { "sediment_!nca",  7, "SITE_KEY",   FTInteger, 6, 0, 0.0, 1.0 },
  { "sediment_!nca",  8, "SAMPLE_KEY", FTInteger, 6, 0, 0.0, 1.0 },
  { "sediment_!nca",  1, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "sediment_!nca",  0, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "sediment_!nca",  2, "WATERDEPTH", FTDouble, 10, 2, 0.0, 1.0 },
  { "sediment_!nca",  3, "SAMPLE_TOP", FTDouble, 10, 2, 0.0, 1.0 },
  { "sediment_!nca",  4, "SAMPLEBASE", FTDouble, 10, 2, 0.0, 1.0 },
  { "sediment_!nca", 23, "CARBONATE%", FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 25, "ORGCARBON%", FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 14, "CLAY_%",     FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 13, "MUD_%",      FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 12, "SAND_%",     FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 11, "GRAVEL_%",   FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 27, "POROSITY_%", FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 26, "SS_LOGKPA",  FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_!nca", 30, "CSS_LOGKPA", FTDouble, 10, 3, 0.0, 1.0 },
  { "sediment_!nca", 15, "GRAINS_PHI", FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 16, "SORTING",    FTDouble, 10, 1, 0.0, 1.0 },
  { "sediment_!nca", 19, "FOLK_CODE",  FTString, 32, 0, 0.0, 1.0 },
  { "sediment_!nca", 20, "SHEPARD_CO", FTString, 32, 0, 0.0, 1.0 },
  { "sediment_!nca",  5, "SITE_NAME",  FTString, 48, 0, 0.0, 1.0 },
  { "sediment_!nca", 31, "SAMPLEPHAS", FTString, 64, 0, 0.0, 1.0 },
  { "sediment_!nca",  9, "SAMPLER",    FTString, 32, 0, 0.0, 1.0 },

  /* GI_BMP_Installations: */

  { "gi_bmp_inst", 4, "BMP_IC(ha)", FTDouble, 20, 6, 0.0, ACRES_TO_HECTARES },
  { "gi_bmp_inst", 1, "SOURCE",     FTString, 16, 0, 0.0, 1.0 },
  { "gi_bmp_inst", 3, "LONGITUDE",  FTDouble, 20, 6, 0.0, 1.0 },
  { "gi_bmp_inst", 2, "LATITUDE",   FTDouble, 20, 6, 0.0, 1.0 },
  { "gi_bmp_inst", 0, "STATE_NAME", FTString, 24, 0, 0.0, 1.0 },

  /* impervious_nhdplus: */

  { "impervious_nhdplus",  7, "IC_TREAT_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "impervious_nhdplus",  8, "IC_CATCH_%", FTDouble, 10, 3, 0.0, 1.0 },
  { "impervious_nhdplus",  9, "IC_TREA_ha", FTDouble, 20, 6, 0.0,
    ACRES_TO_HECTARES },
  { "impervious_nhdplus",  6, "IC_AREA_ha", FTDouble, 20, 6, 0.0,
    ACRES_TO_HECTARES },
  { "impervious_nhdplus",  0, "NHDFLOW_ID", FTInteger, 9, 0, 0.0, 1.0 },
  { "impervious_nhdplus",  1, "GRID_ID",    FTInteger, 9, 0, 0.0, 1.0 },
  { "impervious_nhdplus",  2, "GRID_COUNT", FTInteger,10, 0, 0.0, 1.0 },
  { "impervious_nhdplus",  3, "PROD_UNIT",  FTString,  3, 0, 0.0, 1.0 },
  { "impervious_nhdplus", -1, "AREA_SQKM",  FTDouble, 11, 3, 0.0, 1.0 },
  { "impervious_nhdplus", -1, "HECTARES",   FTDouble, 20, 6, 0.0, 1.0 },

  /* land_change_: */

  { "land_change_",  13, "AveAgricM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  49, "AveAgricPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  43, "AveAgricPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  19, "MedAgricM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  61, "MedAgricPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  55, "MedAgricPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  25, "MinAgricM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  73, "MinAgricPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  67, "MinAgricPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  31, "MaxAgricM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  85, "MaxAgricPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  79, "MaxAgricPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  37, "StdAgricM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  97, "StdAgricPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  91, "StdAgricPC", FTDouble, 5, 1, 0.0, 1.0 },

  { "land_change_",  16, "AveBarreM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  52, "AveBarrePT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  46, "AveBarrePC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  22, "MedBarreM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  64, "MedBarrePT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  58, "MedBarrePC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  28, "MinBarreM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  76, "MinBarrePT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  70, "MinBarrePC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  34, "MaxBarreM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  88, "MaxBarrePT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  82, "MaxBarrePC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  40, "StdBarreM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_", 100, "StdBarrePT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  94, "StdBarrePC", FTDouble, 5, 1, 0.0, 1.0 },

  { "land_change_",  11, "AveForesM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  47, "AveForesPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  41, "AveForesPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  17, "MedForesM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  59, "MedForesPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  53, "MedForesPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  23, "MinForesM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  71, "MinForesPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  65, "MinForesPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  29, "MaxForesM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  83, "MaxForesPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  77, "MaxForesPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  35, "StdForesM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  95, "StdForesPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  89, "StdForesPC", FTDouble, 5, 1, 0.0, 1.0 },

  { "land_change_",  12, "AveShrubM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  48, "AveShrubPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  42, "AveShrubPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  18, "MedShrubM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  60, "MedShrubPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  54, "MedShrubPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  24, "MinShrubM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  72, "MinShrubPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  66, "MinShrubPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  30, "MaxShrubM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  84, "MaxShrubPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  78, "MaxShrubPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  36, "StdShrubM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  96, "StdShrubPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  90, "StdShrubPC", FTDouble, 5, 1, 0.0, 1.0 },

  { "land_change_",  15, "AveWaterM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  51, "AveWaterPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  45, "AveWaterPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  21, "MedWaterM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  63, "MedWaterPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  57, "MedWaterPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  27, "MinWaterM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  75, "MinWaterPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  69, "MinWaterPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  33, "MaxWaterM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  87, "MaxWaterPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  81, "MaxWaterPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  39, "StdWaterM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  99, "StdWaterPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  93, "StdWaterPC", FTDouble, 5, 1, 0.0, 1.0 },

  { "land_change_",  14, "AveWetlaM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  50, "AveWetlaPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  44, "AveWetlaPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  20, "MedWetlaM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  62, "MedWetlaPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  56, "MedWetlaPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  26, "MinWetlaM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  74, "MinWetlaPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  68, "MinWetlaPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  32, "MaxWetlaM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  86, "MaxWetlaPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  80, "MaxWetlaPC", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  38, "StdWetlaM2", FTDouble, 9, 1, 0.0, 1.0 },
  { "land_change_",  98, "StdWetlaPT", FTDouble, 5, 1, 0.0, 1.0 },
  { "land_change_",  92, "StdWetlaPC", FTDouble, 5, 1, 0.0, 1.0 },

  { "land_change_",  1, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "land_change_",  2, "GRID_ID",    FTInteger, 7, 0, 0.0, 1.0 },
  { "land_change_",  3, "GRID_COUNT", FTInteger, 5, 0, 0.0, 1.0 },
  { "land_change_",  4, "PROD_UNIT",  FTString,  3, 0, 0.0, 1.0 },
  { "land_change_", -1, "AREA_SQKM",  FTDouble,  6, 3, 0.0, 1.0 },

  /* stream_temperature_point_median_07_new_england: */

  { "_point_median_07_new_",  7, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_07_new_",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "_point_median_07_new_", 13, "TEMP_OBS_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 32, "RESIDUAL_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 34, "RES_STUD_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 11, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 23, "WATER_%",    FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 12, "SLOPE",      FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 25, "FROM_LAKE",  FTInteger, 1, 0, 0.0, 1.0 },
  { "_point_median_07_new_", 26, "URBAN_HT_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 27, "AIR_TEMP_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 28, "SOLR_Wh/ha", FTDouble, 10, 3, 0.0, 1.0 },

  { "_point_median_07_new_",  3, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_point_median_07_new_",  2, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_07_new_", 10, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "_point_median_07_new_",  6, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 20, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "_point_median_07_new_", 19, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 18, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_07_new_", 24, "StIDYear",   FTString, 10, 0, 0.0, 1.0 },
  { "_point_median_07_new_",  9, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_07_new_",  0, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_07_new_", 17, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_07_new_",  1, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "_point_median_07_new_", 21, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_07_new_", 22, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_07_new_",  4, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_07_new_", 14, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_median_07_new_", 15, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_median_07_new_", 16, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_median_07_new_",  5, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_07_new_", 29, "SLOPE^2",    FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_07_new_", 30, "SLOPE^3",    FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_07_new_", 31, "TEMP_FIT_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 33, "RES_STD_C",  FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 37, "RES_CROV_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 38, "CROV_PRE_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 39, "CROV_ERR_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_07_new_", 35, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "_point_median_07_new_", 36, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_median_08_new_england: */

  { "_point_median_08_new_",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_08_new_",  9, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "_point_median_08_new_", 32, "TEMP_OBS_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 34, "RESIDUAL_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 36, "RES_STUD_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 11, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 24, "WATER_%",    FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 12, "SLOPE",      FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 13, "COARSE_SED", FTDouble, 15, 9, 0.0, 1.0 },
  { "_point_median_08_new_", 27, "FROM_LAKE",  FTInteger, 1, 0, 0.0, 1.0 },
  { "_point_median_08_new_", 26, "LN_DEPTH_m", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 28, "AIR_TEMP_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 29, "SOLR_Wh/ha", FTDouble, 10, 3, 0.0, 1.0 },

  { "_point_median_08_new_",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_point_median_08_new_",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_08_new_", 10, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "_point_median_08_new_",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 21, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "_point_median_08_new_", 20, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 19, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_08_new_", 25, "StIDYear",   FTString, 10, 0, 0.0, 1.0 },
  { "_point_median_08_new_",  8, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_08_new_",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_08_new_", 18, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_08_new_",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "_point_median_08_new_", 22, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_08_new_", 23, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_08_new_",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_median_08_new_", 15, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_median_08_new_", 16, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_median_08_new_", 17, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_median_08_new_",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_08_new_", 30, "SLOPE^2",    FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_08_new_", 31, "SLOPE^3",    FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_median_08_new_", 33, "TEMP_FIT_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 35, "RES_STD_C",  FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 39, "RES_CROV_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 40, "CROV_PRE_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 41, "CROV_ERR_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_median_08_new_", 37, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "_point_median_08_new_", 38, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_mean_daily_range_07_new_england: */

  { "t_mean_daily_range_07_new",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 27, "TEMP_OBS_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 29, "RESIDUAL_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 31, "RES_STUD_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 10, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 11, "IMPERVIOUS", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 12, "SLOPE",      FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 13, "COARSE_SED", FTDouble, 15, 9, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 14, "FLOW_m3/s",  FTDouble, 15, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "t_mean_daily_range_07_new", 15, "WATER_%",    FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 25, "WIDTHDEPTH", FTDouble, 15, 6, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 26, "W/DxATDR_C", FTDouble, 15, 6, 0.0, 1.0 },

  { "t_mean_daily_range_07_new",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 37, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 22, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 21, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 20, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  9, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 19, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 23, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 24, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 16, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 17, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 18, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "t_mean_daily_range_07_new",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 28, "TEMP_FIT_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 30, "RES_STD_C",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 34, "RES_CROV_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 35, "CROV_PRE_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 36, "CROV_ERR_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 32, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "t_mean_daily_range_07_new", 33, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_mean_daily_range_08_new_england: */

  { "t_mean_daily_range_08_new",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 28, "TEMP_OBS_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 30, "RESIDUAL_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 32, "RES_STUD_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 12, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 13, "SLOPE",      FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 15, "FLOW_m3/s",  FTDouble, 15, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },
  { "t_mean_daily_range_08_new", 16, "WATER_%",    FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  9, "AIR_TEMP_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 27, "SOLR_Wh/ha", FTDouble, 10, 3, 0.0, 1.0 },

  { "t_mean_daily_range_08_new",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 11, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 23, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 22, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 21, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 10, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 26, "StIDYear",   FTString, 10, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 20, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 24, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 25, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 17, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 18, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 19, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "t_mean_daily_range_08_new",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 29, "TEMP_FIT_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 31, "RES_STD_C",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 35, "RES_CROV_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 36, "CROV_PRE_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 37, "CROV_ERR_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 33, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "t_mean_daily_range_08_new", 34, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_maximum_daily_increase_07_new_england: */

  { "t_maximum_daily_increase_07_new",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 25, "T_OBS_BOX",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 35, "T_OBS_RAW",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 27, "RESIDUAL",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 29, "RES_STUD",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 11, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 12, "IMPERVIOUS", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 13, "COARSE_SED", FTDouble, 15, 9, 0.0, 1.0 },

  { "t_maximum_daily_increase_07_new",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 10, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 20, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 19, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 18, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  9, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 23, "StIDYrn",    FTString, 10, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 17, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 21, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 22, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 14, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 15, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 16, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 26, "TEMP_FIT",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 28, "RES_STD",    FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 32, "RES_CROV",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 33, "CROV_PRE",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 34, "CROV_ERR",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 30, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "t_maximum_daily_increase_07_new", 31, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_maximum_daily_increase_08_new_england: */

  { "t_maximum_daily_increase_08_new",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 25, "T_OBS_BOX",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 23, "T_OBS_RAW",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 27, "RESIDUAL",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 29, "RES_STUD",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 11, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 12, "IMPERVIOUS", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 13, "COARSE_SED", FTDouble, 15, 9, 0.0, 1.0 },

  { "t_maximum_daily_increase_08_new",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 10, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 20, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 19, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 18, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  9, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 24, "StIDYrn",    FTString, 10, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 17, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 21, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 22, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 14, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 15, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 16, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 26, "TEMP_FIT",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 28, "RES_STD",    FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 32, "RES_CROV",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 33, "CROV_PRE",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 34, "CROV_ERR",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 30, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "t_maximum_daily_increase_08_new", 31, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_maximum_daily_decrease_07_new_england: */

  { "t_maximum_daily_decrease_07_new",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 26, "T_OBS_BOX",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 22, "T_OBS_RAW",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 28, "RESIDUAL",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 30, "RES_STUD",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 11, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 12, "IMPERVIOUS", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 25, "FROM_LAKE",  FTInteger, 1, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 24, "LN_DEPTH_m", FTDouble, 10, 3, 0.0, 1.0 },

  { "t_maximum_daily_decrease_07_new", 27, "TEMP_FIT",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 10, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 19, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 18, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 17, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  9, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 23, "StIDYrn",    FTString, 10, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 16, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 20, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 21, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 13, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 14, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 15, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 29, "RES_STD",    FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 33, "RES_CROV",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 34, "CROV_PRE",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 35, "CROV_ERR",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 31, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "t_maximum_daily_decrease_07_new", 32, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_maximum_daily_decrease_08_new_england: */

  { "t_maximum_daily_decrease_08_new",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 25, "T_OBS_BOX",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 24, "T_OBS_RAW",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 27, "RESIDUAL",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 29, "RES_STUD",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 11, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 12, "IMPERVIOUS", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 13, "SLOPE",      FTDouble, 10, 3, 0.0, 1.0 },

  { "t_maximum_daily_decrease_08_new",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 10, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 20, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 19, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 18, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  9, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 23, "StIDYrn",    FTString, 10, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 17, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 21, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 22, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 14, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 15, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 16, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 26, "TEMP_FIT",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 28, "RES_STD",    FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 32, "RES_CROV",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 33, "CROV_PRE",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 34, "CROV_ERR",   FTDouble, 10, 3, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 30, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "t_maximum_daily_decrease_08_new", 31, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_maximum_new_england: */

  { "_point_maximum_new_eng",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng", 25, "TEMP_OBS_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 27, "RESIDUAL_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 29, "RES_STUD_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 12, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 13, "GRADIENT",   FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 14, "WATER_%",    FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 24, "WIDTHDEPTH", FTDouble, 15, 6, 0.0, 1.0 },
  { "_point_maximum_new_eng",  9, "AIR_TEMP_C", FTDouble, 10, 3, 0.0, 1.0 },

  { "_point_maximum_new_eng",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "_point_maximum_new_eng",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_maximum_new_eng", 11, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "_point_maximum_new_eng",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 21, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "_point_maximum_new_eng", 20, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 19, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_maximum_new_eng", 10, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng", 18, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng", 22, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng", 23, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "_point_maximum_new_eng", 15, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_maximum_new_eng", 16, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_maximum_new_eng", 17, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "_point_maximum_new_eng",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "_point_maximum_new_eng", 26, "TEMP_FIT_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 28, "RES_STD_C",  FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 32, "RES_CROV_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 33, "CROV_PRE_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 34, "CROV_ERR_C", FTDouble, 10, 3, 0.0, 1.0 },
  { "_point_maximum_new_eng", 30, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "_point_maximum_new_eng", 31, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },

  /* stream_temperature_point_day_of_maximum_new_england: */

  { "t_day_of_maximum_new_eng",  0, "Station_ID", FTInteger, 5, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng",  8, "YEAR",       FTInteger, 4, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 26, "DAYMAX_OBS", FTInteger, 4, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 28, "RESIDUAL_D", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 30, "RES_STUD_D", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 12, "DRAIN_km2",  FTDouble, 11, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 13, "DEN_km/km2", FTDouble, 10, 6, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 14, "COARSE_SED", FTDouble, 15, 9, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 15, "FLOW_m3/s",  FTDouble, 15, 5, 0.0,
    CUBIC_FEET_TO_CUBIC_METERS },

  { "t_day_of_maximum_new_eng",  4, "LONGITUDE",  FTDouble, 11, 6, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng",  3, "LATITUDE",   FTDouble, 10, 6, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 11, "Area_m2",    FTDouble, 20, 1, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng",  7, "DistKm",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 22, "afvArea",    FTDouble, 20, 8, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 21, "upDist",     FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 20, "ratio",      FTDouble, 10, 6, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 10, "STUDY_ID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng",  1, "pid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 19, "rid",        FTInteger, 5, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng",  2, "ReachCode",  FTString, 15, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 25, "StIDYrn",    FTString, 10, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 23, "locID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 24, "netID",      FTInteger, 5, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng",  5, "NEAR_FID",   FTInteger, 5, 0, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 16, "NEAR_X",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 17, "NEAR_Y",     FTDouble, 15, 5, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 18, "NEAR_ANGLE", FTDouble, 15, 5, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng",  6, "NEAR_DIST",  FTDouble, 10, 6, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 27, "DAYMAX_FIT", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 29, "RES_STD_D",  FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 33, "RES_CROV_D", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 34, "CROV_PRE_D", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 35, "CROV_ERR_D", FTDouble, 10, 3, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 31, "LEVERAGE",   FTDouble, 15, 8, 0.0, 1.0 },
  { "t_day_of_maximum_new_eng", 32, "CooksDist",  FTDouble, 15, 8, 0.0, 1.0 },


  /* stream_temperature_line_07_new_england: */

  { "stream_temperature_line_07_new_england",  7, "BEAUCLASS",  FTString,   4, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  1, "FDATE",      FTString,  10, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  2, "RESOLUTION", FTString,   6, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  3, "GNIS_ID",    FTString,   8, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  4, "GNIS_NAME",  FTString,  48, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  5, "REACHCODE",  FTString,  14, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  6, "FTYPE",      FTString,  14, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  8, "LTorGT",     FTString,   1, 0, 0.0, 1.0},
  { "stream_temperature_line_07_new_england",  9, "INPUTOOR",   FTInteger,  1, 0, 0.0, 1.0},

  /* stream_temperature_line_08_lower_columbia_river: */

  { "stream_temperature_line_08_lower_columbia_river",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river",  1, "SITE_ID",    FTInteger, 7, 0, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river",  4, "T_2000_CUR", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river",  5, "SE_2000_CU", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river",  6, "T_2000_ADD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river",  7, "SE_2000_AD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river",  8, "T_2000_REM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river",  9, "SE_2000_RM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 10, "T_2040_CUR", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 11, "SE_2040_CU", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 12, "T_2040_ADD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 13, "SE_2040_AD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 14, "T_2040_REM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 15, "SE_2040_RM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 16, "T_2080_CUR", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 17, "SE_2080_CU", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 18, "T_2080_ADD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 19, "SE_2080_AD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 20, "T_2080_REM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_lower_columbia_river", 21, "SE_2080_RM", FTDouble, 10, 3, 0.0, 1.0 },

  /* stream_temperature_line_08_middle_columbia_riverumbia_river: */

  { "stream_temperature_line_08_middle_columbia_river",  0, "COMID",      FTInteger, 10, 0, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river",  1, "SITE_ID",    FTInteger, 7, 0, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river",  4, "T_2000_CUR", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river",  5, "SE_2000_CU", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river",  6, "T_2000_ADD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river",  7, "SE_2000_AD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river",  8, "T_2000_REM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river",  9, "SE_2000_RM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 10, "T_2040_CUR", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 11, "SE_2040_CU", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 12, "T_2040_ADD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 13, "SE_2040_AD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 14, "T_2040_REM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 15, "SE_2040_RM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 16, "T_2080_CUR", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 17, "SE_2080_CU", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 18, "T_2080_ADD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 19, "SE_2080_AD", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 20, "T_2080_REM", FTDouble, 10, 3, 0.0, 1.0 },
  { "stream_temperature_line_08_middle_columbia_river", 21, "SE_2080_RM", FTDouble, 10, 3, 0.0, 1.0 },

  /* stream_temperature_line_meduxnekeag_river: */

  { "stream_temperature_line_meduxnekeag_river",  2, "REACH_ID",   FTInteger, 5, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  1, "NHD_CODE",   FTString,  40, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  3, "2010_07",    FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  5, "2010_08",    FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  7, "2010_MAX",   FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  9, "2010_07_R",  FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river", 11, "2010_MAX_R", FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  4, "2011_07",    FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  6, "2011_08",    FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river",  8, "2011_MAX",   FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river", 10, "2011_07_R",  FTString,  16, 0, 0.0, 1.0},
  { "stream_temperature_line_meduxnekeag_river", 12, "2011_MAX_R", FTString,  16, 0, 0.0, 1.0},

  /* stream_temperature_line_current_08_upper_rogue_river (updated 2024-09-25): */

  { "stream_temperature_line_current_08_upper_rogue_river",  3, "REACH_ID",   FTInteger, 5, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river",  0, "NHDPLUS_ID", FTString, 16, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 12, "SUBSHED_ID", FTInteger, 2, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river",  2, "LENGTH_M",   FTDouble, 16, 11, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river",  5, "CATCH_KM2",  FTDouble,  7, 4, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river",  6, "WSHED_KM2",  FTDouble,  9, 4, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 15, "T_2011_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 19, "T_2011_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 20, "T_2011_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 16, "T_2015_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 21, "T_2015_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 22, "T_2015_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 14, "T_9018_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 17, "T_9018_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 18, "T_9018_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 24, "C_2011_CUR", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 28, "C_2011_ANT", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 29, "C_2011_ORE", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 25, "C_2015_CUR", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 30, "C_2015_ANT", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 31, "C_2015_ORE", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 23, "C_9018_CUR", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 26, "C_9018_ANT", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_current_08_upper_rogue_river", 27, "C_9018_ORE", FTInteger, 3, 0, 0.0, 1.0},

  /* stream_temperature_line_future_08_upper_rogue_river (updated 2024-09-24): */

  { "stream_temperature_line_future_08_upper_rogue_river",  2, "REACH_ID",   FTInteger, 5, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  0, "NHDPLUS_ID", FTString, 16, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  5, "SUBSHED_ID", FTInteger, 2, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  1, "LENGTH_M",   FTDouble, 16, 11, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  3, "CATCH_KM2",  FTDouble,  7, 4, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  4, "WSHED_KM2",  FTDouble,  9, 4, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  6, "T_1992_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  7, "T_1992_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  8, "T_1992_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river",  9, "T_1993_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 10, "T_1993_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 11, "T_1993_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 12, "T_1995_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 13, "T_1995_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 14, "T_1995_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 15, "T_1997_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 16, "T_1997_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 17, "T_1997_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 18, "T_1998_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 19, "T_1998_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 20, "T_1998_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 21, "T_2002_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 22, "T_2002_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 23, "T_2002_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 24, "T_2005_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 25, "T_2005_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 26, "T_2005_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 27, "T_2007_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 28, "T_2007_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 29, "T_2007_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 30, "T_2011_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 31, "T_2011_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 32, "T_2011_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 33, "T_2012_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 34, "T_2012_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 35, "T_2012_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 36, "T_9018_CUR", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 37, "T_9018_ANT", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 38, "T_9018_ORE", FTDouble, 17, 10, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 39, "C_1997_CUR", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 40, "C_1997_ANT", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 41, "C_1997_ORE", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 42, "C_2002_CUR", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 43, "C_2002_ANT", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 44, "C_2002_ORE", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 45, "C_9018_CUR", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 46, "C_9018_ANT", FTInteger, 3, 0, 0.0, 1.0},
  { "stream_temperature_line_future_08_upper_rogue_river", 47, "C_9018_ORE", FTInteger, 3, 0, 0.0, 1.0},

  /* stream_temperature_point_08_lower_columbia_river, mid_columbia_river: */

  { "stream_temperature_point_08_", 14, "TRIB_ID",   FTDouble, 4, 0, 0.0, 1.0},
  { "stream_temperature_point_08_",  2, "SITE_ID",   FTInteger, 7, 0, 0.0, 1.0},
  { "stream_temperature_point_08_",  0, "TRIBUTARY", FTString, 84, 0, 0.0, 1.0},
  { "stream_temperature_point_08_",  1, "MI_TO_SEA", FTDouble, 8, 3, 0.0, 1.0 },
  { "stream_temperature_point_08_",  1, "KM_TO_SEA", FTDouble, 8, 3, 0.0,
    MILES_TO_KM
  },
  { "stream_temperature_point_08_",  3, "T_2000_CUR", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_",  4, "T_2000_ADD", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_",  5, "T_2000_REM", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_",  6, "T_2040_CUR", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_",  7, "T_2040_ADD", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_",  8, "T_2040_REM", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_",  9, "T_2080_CUR", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_", 10, "T_2080_ADD", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_", 11, "T_2080_REM", FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_", 12, "FLOW_CMS",   FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_", 13, "FLOW_CFS",   FTDouble, 10, 3, 0.0,1.0},
  { "stream_temperature_point_08_", 16, "LONGITUDE",  FTDouble, 11, 6, 0.0,1.0},
  { "stream_temperature_point_08_", 15, "LATITUDE",   FTDouble, 10, 6, 0.0,1.0},

  /* stream_shade_line_08_lower_columbia_river: */

  { "stream_shade_line_08_lower_columbia_river",  4, "SITE_ID",    FTInteger, 7, 0, 0.0, 1.0},
  { "stream_shade_line_08_lower_columbia_river",  0, "BANK_WIDTH", FTDouble, 7, 2, 0.0, 1.0 },
  { "stream_shade_line_08_lower_columbia_river",  2, "SHADE_CURV", FTDouble, 6, 2, 0.0, 1.0 },
  { "stream_shade_line_08_lower_columbia_river",  3, "SHADE_ADDV", FTDouble, 6, 2, 0.0, 1.0 },
  { "stream_shade_line_08_lower_columbia_river",  1, "SHADE_REMV", FTDouble, 6, 2, 0.0, 1.0 },

  /* stream_shade_line_08_middle_columbia_riverumbia_river: */

  { "stream_shade_line_08_middle_columbia_river",  4, "SITE_ID",    FTInteger, 7, 0, 0.0,1.0},
  { "stream_shade_line_08_middle_columbia_river",  3, "BANK_WIDTH", FTDouble, 7, 2, 0.0, 1.0},
  { "stream_shade_line_08_middle_columbia_river",  1, "SHADE_CURV", FTDouble, 6, 2, 0.0, 1.0},
  { "stream_shade_line_08_middle_columbia_river",  2, "SHADE_ADDV", FTDouble, 6, 2, 0.0, 1.0},
  { "stream_shade_line_08_middle_columbia_river",  0, "SHADE_REMV", FTDouble, 6, 2, 0.0, 1.0},



  /*
   * HSPF: hspf_charles3_huc10_polygon_atlantic, etc.
   * Note: don't change column names or order because the WMOST model
   * which will read these files requires
   * these column names in this order:
   */

  { "/hspf_",  0, "HRU_ID",     FTInteger,  4, 0, 0.0, 1.0 },
  { "/hspf_",  1, "KGW",        FTDouble,  12, 9, 0.0, 1.0 },
  { "/hspf_",  2, "EIA",        FTDouble,  10, 6, 0.0, 1.0 },
  { "/hspf_",  3, "INFILT",     FTDouble,  10, 6, 0.0, 1.0 },
  { "/hspf_",  4, "HRU_NAME",   FTString,  80, 0, 0.0, 1.0 },
  { "/hspf_",  5, "ACRES",      FTDouble,  10, 3, 0.0, 1.0 },
  { "/hspf_",  6, "START_DATE", FTString,  10, 0, 0.0, 1.0 },
  { "/hspf_",  7, "END_DATE",   FTString,  10, 0, 0.0, 1.0 },
  { "/hspf_",  8, "LAT_MODEL",  FTDouble,  10, 4, 0.0, 1.0 },
  { "/hspf_",  9, "MODEL_ID",   FTInteger,  4, 0, 0.0, 1.0 },
  { "/hspf_", 10, "SUB_ID",     FTInteger,  4, 0, 0.0, 1.0 },
  { "/hspf_", 11, "HUC_ID",     FTDouble,  14, 1, 0.0, 1.0 },
  { "/hspf_", 12, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "/hspf_", 13, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },

  { "/swmm_",  0, "HRU_ID",     FTInteger,  4, 0, 0.0, 1.0 },
  { "/swmm_",  1, "KGW",        FTDouble,  12, 9, 0.0, 1.0 },
  { "/swmm_",  2, "EIA",        FTDouble,  10, 6, 0.0, 1.0 },
  { "/swmm_",  3, "INFILT",     FTDouble,  10, 6, 0.0, 1.0 },
  { "/swmm_",  4, "HRU_NAME",   FTString,  80, 0, 0.0, 1.0 },
  { "/swmm_",  5, "ACRES",      FTDouble,  10, 3, 0.0, 1.0 },
  { "/swmm_",  6, "START_DATE", FTString,  10, 0, 0.0, 1.0 },
  { "/swmm_",  7, "END_DATE",   FTString,  10, 0, 0.0, 1.0 },
  { "/swmm_",  8, "LAT_MODEL",  FTDouble,  10, 4, 0.0, 1.0 },
  { "/swmm_",  9, "MODEL_ID",   FTInteger,  4, 0, 0.0, 1.0 },
  { "/swmm_", 10, "SUB_ID",     FTInteger,  4, 0, 0.0, 1.0 },
  { "/swmm_", 11, "HUC_ID",     FTDouble,  14, 1, 0.0, 1.0 },
  { "/swmm_", 12, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "/swmm_", 13, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },

  { "/swat_",  0, "HRU_ID",     FTInteger,  4, 0, 0.0, 1.0 },
  { "/swat_",  1, "KGW",        FTDouble,  12, 9, 0.0, 1.0 },
  { "/swat_",  2, "EIA",        FTDouble,  10, 6, 0.0, 1.0 },
  { "/swat_",  3, "INFILT",     FTDouble,  10, 6, 0.0, 1.0 },
  { "/swat_",  4, "HRU_NAME",   FTString,  80, 0, 0.0, 1.0 },
  { "/swat_",  5, "ACRES",      FTDouble,  10, 3, 0.0, 1.0 },
  { "/swat_",  6, "START_DATE", FTString,  10, 0, 0.0, 1.0 },
  { "/swat_",  7, "END_DATE",   FTString,  10, 0, 0.0, 1.0 },
  { "/swat_",  8, "LAT_MODEL",  FTDouble,  10, 4, 0.0, 1.0 },
  { "/swat_",  9, "MODEL_ID",   FTInteger,  4, 0, 0.0, 1.0 },
  { "/swat_", 10, "SUB_ID",     FTInteger,  4, 0, 0.0, 1.0 },
  { "/swat_", 11, "HUC_ID",     FTDouble,  14, 1, 0.0, 1.0 },
  { "/swat_", 12, "LATITUDE",   FTDouble,  10, 6, 0.0, 1.0 },
  { "/swat_", 13, "LONGITUDE",  FTDouble,  11, 6, 0.0, 1.0 },

  { "greenspace_housing", 13, "MHVKMmean",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 14, "MHVKLmean",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 15, "MHVKOmean",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 16, "MHVKHmean",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 21, "ptCanopyH",  FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 22, "ptCanopyM",  FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 23, "ptCanopyL",  FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 24, "ptCanopyO",  FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 33, "Can0_250m",  FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 34, "Can250_500", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 47, "hu10pha_21", FTDouble,  10, 6, 0.0, 1.0 },
  { "greenspace_housing", 49, "hu10pha_22", FTDouble,  10, 6, 0.0, 1.0 },
  { "greenspace_housing", 51, "hu10pha_23", FTDouble,  10, 6, 0.0, 1.0 },
  { "greenspace_housing", 53, "hu10pha_24", FTDouble,  10, 0, 0.0, 1.0 },
  { "greenspace_housing", 55, "m2_21_0250", FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 56, "m2_22_0250", FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 57, "m2_23_0250", FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 58, "m2_24_0250", FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 67, "m221250500", FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 68, "m222250500", FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 69, "m223250500", FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 70, "m224250500", FTInteger, 10, 0, 0.0, 1.0 },

  { "greenspace_housing",  4, "HUC12",      FTString,  12, 0, 0.0, 1.0 },
  { "greenspace_housing",  0, "GNIS_ID",    FTString,   8, 0, 0.0, 1.0 },
  { "greenspace_housing", 17, "PD10p900H",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 18, "PD10p900M",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 19, "PD10p900L",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 20, "PD10p900O",  FTDouble,  10, 4, 0.0, 1.0 },
  { "greenspace_housing", 25, "ptCan250mH", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 26, "ptCan250mM", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 27, "ptCan250mL", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 28, "ptCan250mO", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 29, "ptCan500mH", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 30, "ptCan500mM", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 31, "ptCan500mL", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 32, "ptCan500mO", FTDouble,   6, 1, 0.0, 1.0 },
  { "greenspace_housing", 35, "AD21_0250m", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 36, "AD22_0250m", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 37, "AD23_0250m", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 38, "AD24_0250m", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 39, "AD21250500", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 40, "AD22250500", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 41, "AD23250500", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 42, "AD24250500", FTInteger,  4, 0, 0.0, 1.0 },
  { "greenspace_housing", 43, "AREAM2_21",  FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 44, "AREAM2_22",  FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 45, "AREAM2_23",  FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 46, "AREAM2_24",  FTInteger, 10, 0, 0.0, 1.0 },
  { "greenspace_housing", 48, "haphu10_21", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 50, "haphu10_22", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 52, "haphu10_23", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 54, "haphu10_24", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 59, "hu_21_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 60, "hh_21_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 61, "hu_22_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 62, "hh_22_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 63, "hu_23_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 64, "hh_23_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 65, "hu_24_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 66, "hh_24_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 71, "hu21250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 72, "hh21250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 73, "hu22250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 74, "hh22250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 75, "hu23250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 76, "hh23250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 77, "hu24250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 78, "hh24250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 79, "hapbldg_21", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 80, "hapbldg_22", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 81, "hapbldg_23", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 82, "hapbldg_24", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 83, "hb_21_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 84, "hb_22_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 85, "hb_23_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 86, "hb_24_0250", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 87, "hb21250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 88, "hb22250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 89, "hb23250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing", 90, "hb24250500", FTDouble,  12, 6, 0.0, 1.0 },
  { "greenspace_housing",  1, "AREAACRES",  FTDouble,  12, 3, 0.0, 1.0 },
  { "greenspace_housing",  2, "AREASQKM",   FTDouble,  12, 3, 0.0, 1.0 },
  { "greenspace_housing",  3, "STATES",     FTString,  20, 0, 0.0, 1.0 },
  { "greenspace_housing",  6, "HUTYPE",     FTString,  12, 0, 0.0, 1.0 },
  { "greenspace_housing",  7, "HUMOD",      FTString,  32, 0, 0.0, 1.0 },
  { "greenspace_housing",  8, "TOHUC",      FTString,  12, 0, 0.0, 1.0 },
  { "greenspace_housing",  9, "NONCONTRIB", FTDouble,  12, 3, 0.0, 1.0 },
  { "greenspace_housing", 10, "NONCONTR_1", FTDouble,  12, 3, 0.0, 1.0 },
  { "greenspace_housing", 11, "SHAPE_Leng", FTDouble,  20, 6, 0.0, 1.0 },
  { "greenspace_housing", 12, "SHAPE_Area", FTDouble,  20, 6, 0.0, 1.0 },
  { "greenspace_housing",  5, "NAME",       FTString,  80, 0, 0.0, 1.0 },

  /* HMS Smoke files used by RSIG: */

  { "hms_smoke", -1, "YYYYDDD1",  FTInteger, 7, 0, 0.0, 1.0 },
  { "hms_smoke", -1, "HHMM1",     FTInteger, 4, 0, 0.0, 1.0 },
  { "hms_smoke", -1, "YYYYDDD2",  FTInteger, 7, 0, 0.0, 1.0 },
  { "hms_smoke", -1, "HHMM2",     FTInteger, 4, 0, 0.0, 1.0 },
  { "hms_smoke", -1, "DENS_UGM3", FTInteger, 2, 0, 0.0, 1.0 },

  /* End of table: */

  { 0, 0, 0, 0, 0, 0, 0.0, 0.0 } /* End of table. */
};


/******************************************************************************
PURPOSE: defineDBFColumns - Define columns of an output DBF file.
INPUTS:  const char* inputFileName  File name of input dbf.
         const int defineColumns    If 1 then define columns else just compute
                                    tableIndex and column count.
OUTPUTS: int* tableIndex            Index into table of first column.
         int* longitudeColumn       Index of column LONGITUDE, else -1 if none.
         int* latitudeColumn        Index of column LATITUDE,  else -1 if none.
         int* hucIdColumn           Index of column HUC_ID,    else -1 if none.
         int* estcodeColumn         Index of column ESTCODE,   else -1 if none.
         int* siteIdColumn          Index of column SITE_ID,   else -1 if none.
         DBFHandle outputFile       DBF file to write to.
RETURNS: int number of columns defined if successful,
         else 0 and a failure message is printed to stderr.
******************************************************************************/

static int defineDBFColumns( const char* const inputFileName,
                             const int defineColumns,
                             int* const tableIndex,
                             int* const longitudeColumn,
                             int* const latitudeColumn,
                             int* const hucIdColumn,
                             int* const estcodeColumn,
                             int* const siteIdColumn,
                             DBFHandle outputFile ) {
  PRE05( inputFileName, *inputFileName, IS_BOOL( defineColumns ), tableIndex,
         outputFile );
  int result = 0;
  int ok = 0;
  int entry = 0;
  char lowercaseInputFileName[ 256 ] = "";
  memset( lowercaseInputFileName, 0, sizeof lowercaseInputFileName );
  strncpy( lowercaseInputFileName, inputFileName,
           sizeof lowercaseInputFileName /
           sizeof *lowercaseInputFileName - 1 );
  lowercase( lowercaseInputFileName );
  *tableIndex = -1;

  for ( entry = 0, ok = 1; AND2( ok, table[ entry ].fileName ); ++entry ) {
    const ColumnEntry* const columnEntry = table + entry;
    int matchedColumn = 0;

    /* Handle possible filter such as "tributaries_!great_lakes": */

    {
      char columnEntryFileName[ 256 ] = "";
      memset( columnEntryFileName, 0, sizeof columnEntryFileName );
      strncpy( columnEntryFileName, columnEntry->fileName,
               sizeof columnEntryFileName / sizeof *columnEntryFileName - 1 );
      char* except = strchr( columnEntryFileName, '!' );

      if ( except ) {
        *except = '\0';
        ++except;
      }

      matchedColumn =
        AND2( strstr( lowercaseInputFileName, columnEntryFileName ),
              OR2( except == 0,
                   strstr( lowercaseInputFileName, except ) == 0 ) );
    }

    if ( matchedColumn ) {
      ++result;

      if ( *tableIndex == -1 ) {
        *tableIndex = entry;
      }

      if ( longitudeColumn ) {

        if ( ! strcmp( columnEntry->columnName, "LONGITUDE" ) ) {
          *longitudeColumn = columnEntry->inputColumn;
        }
      }

      if ( latitudeColumn ) {

        if ( ! strcmp( columnEntry->columnName, "LATITUDE" ) ) {
          *latitudeColumn = columnEntry->inputColumn;
        }
      }

      if ( hucIdColumn ) {

        if ( OR2( ! strcmp( columnEntry->columnName, "HUC_ID" ),
                  ! strcmp( columnEntry->columnName, "HUC12" ) ) ) {
          *hucIdColumn = columnEntry->inputColumn;
        }
      }

      if ( estcodeColumn ) {

        if ( ! strcmp( columnEntry->columnName, "ESTCODE" ) ) {
          *estcodeColumn = columnEntry->inputColumn;
        } else if ( ! strcmp( columnEntry->columnName, "ESTCODE_N" ) ) {
          *estcodeColumn = columnEntry->inputColumn;
        } else if ( ! strcmp( columnEntry->columnName, "WATERSHED" ) ) {
          *estcodeColumn = columnEntry->inputColumn;
        }
      }

      if ( siteIdColumn ) {

        if ( ! strcmp( columnEntry->columnName, "SITE_ID" ) ) {
          *siteIdColumn = columnEntry->inputColumn;
        }
      }

      if ( defineColumns ) {
        ok = DBFAddField( outputFile,
                          columnEntry->columnName,
                          columnEntry->columnType,
                          columnEntry->fieldWidth,
                          columnEntry->decimals ) != -1;
      }
    }
  }

  ok = AND3( ok,
             IN_RANGE( result,      1, sizeof table / sizeof *table - 1 ),
             IN_RANGE( *tableIndex, 0, sizeof table / sizeof *table - 2 ) );

  if ( ! ok ) {
    fprintf( stderr, "\nFailed to define DBF columns for file %s.\n",
             inputFileName );
    result = 0;
    *tableIndex = -1;

    if ( longitudeColumn ) {
      *longitudeColumn = -1;
    }

    if ( latitudeColumn ) {
      *latitudeColumn = -1;
    }

    if ( hucIdColumn ) {
      *hucIdColumn = -1;
    }

    if ( estcodeColumn ) {
      *estcodeColumn = -1;
    }
  }

  POST0( OR2( AND7( result == 0,
                    *tableIndex == -1,
                    IMPLIES( longitudeColumn, *longitudeColumn == -1 ),
                    IMPLIES( latitudeColumn,  *latitudeColumn  == -1 ),
                    IMPLIES( hucIdColumn,     *hucIdColumn     == -1 ),
                    IMPLIES( estcodeColumn,   *estcodeColumn   == -1 ),
                    IMPLIES( siteIdColumn,    *siteIdColumn    == -1 ) ),
              AND2( result > 0,
                    *tableIndex >= 0 ) ) );
                        
  return result;
}



/******************************************************************************
PURPOSE: matchesWithUnderscores - Does the first string match the second
         string after converting any spaces in the seconds string to
         underscores and ignoring case?
INPUTS:  const char* withUnderscores  String with possible underscores.
                                      E.g., "Liberty_Bay".
         const char* value            String with possible spaces.
                                      E.g., "Liberty Bay".
RETURNS: int 1 if strings match else 0.
******************************************************************************/

static int matchesWithUnderscores( const char* withUnderscores,
                                   const char* value ) {

  PRE02( withUnderscores, *withUnderscores );
  int result = 0;

  if ( AND3( value, *value, strlen( value ) == strlen( withUnderscores ) ) ) {
    const char* u = withUnderscores;
    const char* v = value;

    do {
      const char uc = tolower( *u );
      const char vc = tolower( *v );

      result = OR2( vc == uc, AND2( vc == ' ', uc == '_' ) );
      ++u;
      ++v;
    } while ( AND2( result, *u ) );
  }

  return result;
}



/*=========================== FORWARD DECLARATIONS ==========================*/

static int writePRJFile( const char* fileName, const int useASCIIGridForm );

static int writeDataToDBFFile( const char* fileName, const char* variable,
                               const char* units,
                               const int timesteps, const int yyyymmddhh[],
                               const int timestepType,
                               const int count, const int type,
                               const int components,
                               const void* data,
                               const float lonlats[], const float z[],
                               const char* const sids[],
                               const int ids[],
                               const char* const metadata[],
                               const int writeCSV );

static int writeDataToCSVFile( const char* fileName, const char* variable,
                               const char* units,
                               const int timesteps, const int yyyymmddhh[],
                               const int timestepType,
                               const int count, const int type,
                               const int components,
                               const void* data,
                               const float lonlats[], const float z[],
                               const char* const sids[],
                               const int ids[] );

static int appendCSVFile( const char* inputCSVDirectory, const char* code,
                          const int yyyymmdd1, const int yyyymmdd2,
                          const int timestepSize,
                          const int columns, CSVHeader header,
                          FILE* outputCSVFile );

static int writeCSVLine( char* line, const int columns, FILE* outputCSVFile );

static char* parseCSVHeader( char* line, const int columns, CSVHeader header );

static int parseCSVTimestamp( char* line, int* timestamp );

static int hasIntegerColumn( const ShapeData* const shapeData,
                             const char* const name );

static int flagUpstreamNodes( const int rows,
                              const int columns,
                              const Value values[],
                              const int fromNodeColumn,
                              const int toNodeColumn,
                              const int toNode,
                              char mask[] );

static size_t copyMatchedLines( FILE* input,
                                const size_t count,
                                const long long values[],
                                FILE* output );

static int compareDate( const int yyyymmdd1, const int yyyymmdd2,
                        const int timestepSize );

static const char* storeStringValue( const char* value, ShapeData* shapeData );

static int copyStringAttribute( DBFHandle inputFile,
                                int inputRow, int inputColumn,
                                DBFHandle outputFile,
                                int outputRow, int outputColumn );

static int copyIntegerAttribute( DBFHandle inputFile,
                                 int inputRow, int inputColumn,
                                 DBFHandle outputFile,
                                 int outputRow, int outputColumn );

static int copyDoubleAttribute( DBFHandle inputFile,
                                int inputRow, int inputColumn,
                                DBFHandle outputFile,
                                int outputRow, int outputColumn,
                                int invalidIfNegative, double missing,
                                double offset, double scale );

static int computeSparsedPartCount( const SHPObject* shape,
                                    const double minimumAdjacentVertexDistance,
                                    const int minimumSparsedVertices );

static int computeSparsedVertexCount( const int vertexCount,
                                      const double x[],
                                      const double y[],
                                      const double minimumAdjacentVertexDistance,
                                      const int isPolygon );

static void copySparseVertices( const int vertexCount,
                                const int sparsedVertices,
                                const double x[],
                                const double y[],
                                const double minimumAdjacentVertexDistance,
                                const int isPolygon,
                                int* const initializedBounds,
                                Bounds bounds,
                                gpc_vertex* const vertices );

static void computeGridBounds( const int rows,
                               const int columns,
                               const double westEdge,
                               const double southEdge,
                               const double cellWidth,
                               const double cellHeight,
                               Unproject unproject,
                               double range[ 4 ] );

static void updatePointMinmax( Unproject unproject,
                               const double x,
                               const double y,
                               double* minimumX,
                               double* maximumX,
                               double* minimumY,
                               double* maximumY );

static void computeGridCellCenters( const int rows,
                                    const int columns,
                                    const double westEdge,
                                    const double southEdge,
                                    const double cellWidth,
                                    const double cellHeight,
                                    Unproject unproject,
                                    float lonlats[] );

static void computePointBounds( const int points,
                                const float lonlats[],
                                double range[ 4 ] );

static void computeVertexBounds( const int count, const float xy[],
                                 double xyRange[] );

static void computeRange( const double array[], const int count,
                          const int stride, double range[ 2 ] );

static int convertTimestamp( const char* inputFileName,
                             const char* timestamp,
                             int* const yyyyddd, int* const hhmm );

static int writeShort( unsigned char bytes[], const int index,
                       const short value, const int endian );

static int writeInt( unsigned char bytes[], const int index,
                     const int value, const int endian );

static int writeDouble( unsigned char bytes[], const int index,
                        const double value, const int endian );

static int minimumInt( int count, const int array[] );

static int maximumInt( int count, const int array[] );



/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: writeASCIIGridFile - Write a single timestep-layer of grid cell scalar
         data to an ESRI ASCII Grid file.
INPUTS:  const char* const fileName  File to create. E.g.,"example.asc".
         int rows              Number of grid rows.
         int columns           Number of grid columns.
         double westEdge       Distance from origin to west edge of grid.
         double southEdge      Distance from origin to south edge of grid.
         double cellSize       Width/height of each grid cell.
         int type              Grid cell scalar data type: BYTE_TYPE, etc.
         const void* data      data[ rows * columns ] Scalar at cell centers.
NOTES:   Creates file fileName.
         http://en.wikipedia.org/wiki/ESRI_grid
******************************************************************************/

int writeASCIIGridFile( const char* fileName, int rows, int columns,
                        double westEdge, double southEdge, double cellSize,
                        int type, const void* data ) {
  PRE012( fileName,
          *fileName,
          rows > 0,
          columns > 0,
          rows < INT_MAX / columns,
          ! is_nan( westEdge ),
          ! is_nan( southEdge ),
          ! is_nan( cellSize ),
          cellSize > 0.0,
          IS_VALID_GRID_DATA_TYPE( type ),
          data,
          IMPLIES( type == FLOAT_TYPE,
                   AND2( ! is_nan( ((const float*) data)[ 0 ] ),
                         ! is_nan(((const float*)data)[rows * columns - 1]))));
  int result = 0;
  FILE* file = fopen( fileName, "w" );

  if ( file ) {
    const size_t last = ( rows - 1 ) * columns;
    const float* const fdata = data;
    const char* const cdata = data;
    const unsigned short* const sdata = data;
    const float* fDataRow = fdata + last;
    const char*  cDataRow = cdata + last;
    const unsigned short* sDataRow = sdata + last;
    const char* const noDataValue = type == UINT16_TYPE ? "0" : "-9999";
    int row = 0;
    result =
      fprintf( file, "ncols %d\nnrows %d\nxllcorner %lg\nyllcorner %lg\n"
                     "cellsize %lg\nNODATA_value %s\n",
               columns, rows, westEdge, southEdge, cellSize, noDataValue ) > 0;

    for ( row = rows - 1; AND2( result, row >= 0 ); --row,
          fDataRow -= columns, cDataRow -= columns, sDataRow -= columns ) {
      int column = 0;

      for ( column = 0; AND2( result, column < columns ); ++column ) {
        const float dataValue =
          type == FLOAT_TYPE ? fDataRow[ column ]
          : type == BYTE_TYPE ? (float) ( cDataRow[ column ] )
          : (float) ( sDataRow[ column ] );
        const float clampedValue = dataValue > -9999.0 ? dataValue : -9999.0;
        result = fprintf( file, "%g ", clampedValue ) > 0;
      }

      result = AND2( result, fprintf( file, "\n" ) > 0 );
    }

    fclose( file ), file = 0;
  }

  if ( ! result ) {
    perror( "\n\nFailed because" );
  }

  return result;
}



/******************************************************************************
PURPOSE: writeWGS84PRJFile - Write a WGS84 ESRI projection file.
INPUTS:  const char* const fileName  File to create. E.g.,"example.prj".
         const int         useASCIIGridForm  1 = use grid form, 0 = shape form.
NOTES:   Creates file fileName.
         http://en.wikipedia.org/wiki/Well-known_text
******************************************************************************/

int writeWGS84PRJFile( const char* fileName, const int useASCIIGridForm ) {
  PRE03( fileName, *fileName, IS_BOOL( useASCIIGridForm ) );
  const char* const content =
    useASCIIGridForm ?
      "Projection    GEOGRAPHIC\n"
      "Datum         NAD83\n"
      "Spheroid      GRS80\n"
      "Units         DD\n"
      "Zunits        NO\n"
    :
      "GEOGCS[\"GCS_North_American_1983\","
      "DATUM[\"D_North_American_1983\","
      "SPHEROID[\"GRS_1980\",6378137.0,298.257223563]],"
      "PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]]";
  int result = 0;
  FILE* file = fopen( fileName, "w" );

  if ( file ) {
    result = fprintf( file, "%s", content ) == strlen( content );
    fclose( file ), file = 0;
  }

  if ( ! result ) {
    perror( "\n\nFailed because" );
  }

  return result;
}



/******************************************************************************
PURPOSE: writeLambertPRJFile - Write a Lambert ESRI projection file.
INPUTS:  const char* const fileName  File to create. E.g.,"example.prj".
         double centralLongitude     Longitude of center of projection. -90.0.
         double centralLatitude      Latitude of center of projection. 40.0.
         double lowerLatitude        Latitude of lower tangent. 30.0.
         double upperLatitude        Latitude of upper tangent. 60.0.
         const int useASCIIGridForm  1 = use grid form, 0 = shape form.
NOTES:   Creates file fileName. Uses MM5 sphere of radius 6,370,000m.
         http://en.wikipedia.org/wiki/Well-known_text
******************************************************************************/

int writeLambertPRJFile( const char* fileName,
                         const double centralLongitude,
                         const double centralLatitude,
                         const double lowerLatitude,
                         const double upperLatitude,
                         const int useASCIIGridForm ) {
  PRE011( fileName,
          *fileName,
          is_valid_longitude_latitude( centralLongitude, centralLatitude ),
          IN_RANGE( centralLatitude, -89.0, 89.0 ),
          is_valid_latitude( lowerLatitude ),
          is_valid_latitude( upperLatitude ),
          lowerLatitude <= upperLatitude,
          SIGN( lowerLatitude ) == SIGN( upperLatitude ),
          IMPLIES_ELSE( lowerLatitude >= 0.0,
                        IN_RANGE( lowerLatitude, 1.0, 89.0 ),
                        IN_RANGE( lowerLatitude, -89.0, -1.0 ) ),
          IMPLIES_ELSE( upperLatitude >= 0.0,
                        IN_RANGE( upperLatitude, 1.0, 89.0 ),
                        IN_RANGE( upperLatitude, -89.0, -1.0 ) ),
          IS_BOOL( useASCIIGridForm ) );
  const char* const content =
    useASCIIGridForm ?
      "Projection    Lambert Conformal Conic\n"
      "Datum         NAD83\n"
      "Spheroid      GRS80\n"
      "Units         METERS\n"
      "Zunits        NO\n"
      "Xshift        0.0\n"
      "Yshift        0.0\n"
      "Parameters\n"
      "%lg /* 1st standard parallel */\n"
      "%lg /* 2nd standard parallel */\n"
      "%lg /* central meridian */\n"
      "%lg /* latitude of projection's origin */\n"
      "0.0 /* false easting (meters) */\n"
      "0.0 /* false northing (meters) */\n"
    :
      "PROJCS[\"Lambert Conformal Conic\","
      "GEOGCS[\"GCS_Sphere_EMEP\","
      "DATUM[\"D_Sphere_EMEP\","
      "SPHEROID[\"Sphere_EMEP\",6370000.0,0.0]],"
      "PRIMEM[\"Greenwich\",0.0],"
      "UNIT[\"Degree\",0.0174532925199433]],"
      "PROJECTION[\"Lambert_Conformal_Conic\"],"
      "PARAMETER[\"Standard_Parallel_1\",%lg],"
      "PARAMETER[\"Standard_Parallel_2\",%lg],"
      "PARAMETER[\"Latitude_Of_Origin\",%lg],"
      "PARAMETER[\"Central_Meridian\",%lg],"
      "PARAMETER[\"False_Easting\",0.0],"
      "PARAMETER[\"False_Northing\",0.0],"
      "UNIT[\"Meter\",1]]";
  int result = 0;
  FILE* file = fopen( fileName, "w" );

  if ( file ) {
    result = fprintf( file, content,
                      lowerLatitude, upperLatitude,
                      centralLatitude, centralLongitude )
             > strlen( content ) - 26;
    fclose( file ), file = 0;
  }

  if ( ! result ) {
    perror( "\n\nFailed because" );
  }

  return result;
}



/******************************************************************************
PURPOSE: writeGridToShapefile - Write a single layer of grid cells as a lon-lat
         Shapefile Polygon file set (shp, shx, dbf, prj) and a csv file
         containing time-varying data.
INPUTS:  const char* const fileName  Base name of file to create. "example".
         const int timesteps         Number of timesteps of data.
         const int yyyymmddhh[ timesteps ]  Timestamps of data.
         const int timestepType      HOURLY, DAILY, MONTHLY, YEARLY.
         const int rows              Number of grid rows.
         const int columns           Number of grid columns.
         const double westEdge       Distance from origin to west edge of grid.
         const double southEdge      Distance from origin to south edge of ".
         const double cellWidth      Width of each grid cell (e.g., 12000 m).
         const double cellWHeight    Height of each grid cell (e.g., 12000 m).
         const char* variable        Name of data variable.
         const char* units           Name of data units.
         int type                    Grid cell scalar data type: BYTE_TYPE...
         const void* data            data[ timesteps * rows * columns ]
                                     Scalar data at grid cell centers.
         Unproject unproject         Function to unproject (x,y) to (lon,lat).
RETURNS: int 1 if successful, else 0 and failure message is printed to stderr.
NOTES:   Creates files fileName.shp and fileName.shx
         then calls routines that create fileName.dbf and fileName.csv.
         See 1998 ESRI Shapefile Specification pages 2, 4, 5, 16, 23, 24.
         http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
         This routine does not use shapelib API since it was written before
         that library was used in this project and also because this routine
         is more straight-forward and much faster than using the shapelib API.
******************************************************************************/

int writeGridToShapefile( const char* fileName,
                          const int timesteps,
                          const int yyyymmddhh[],
                          const int timestepType,
                          const int rows,
                          const int columns,
                          const double westEdge,
                          const double southEdge,
                          const double cellWidth,
                          const double cellHeight,
                          const char* variable,
                          const char* units,
                          int type,
                          int components,
                          const void* data,
                          Unproject unproject ) {

  PRE016( fileName,
          *fileName,
          timesteps > 0,
          yyyymmddhh,
          IS_VALID_TIMESTEP_TYPE( timestepType ),
          rows > 0,
          columns > 0,
          rows * columns > 0,
          variable,
          *variable,
          units,
          *units,
          IS_VALID_GRID_DATA_TYPE( type ),
          data,
          IMPLIES( type == FLOAT_TYPE,
                   AND2( ! is_nan( ((const float*) data)[ 0 ] ),
             ! is_nan(((const float*)data)[timesteps * rows * columns - 1]))),
          IMPLIES( unproject == 0,
                   AND2 ( is_valid_longitude_latitude( westEdge, southEdge ),
                          is_valid_longitude_latitude(
                                      westEdge  + columns * cellWidth,
                                      southEdge + rows    * cellHeight ) ) ) );

  int result = 0;
  enum { /* Constants from the Shapefile Spec (see above link): */
    BYTES_PER_INT = 4,
    BYTES_PER_DOUBLE = 8,
    POLYGON = 5,                 /* ShapeType = POLYGON. */
    PARTS_PER_POLYGON = 1,
    VERTICES_PER_POLYGON = 5,    /* Quads with redundant last=first vertex.*/
    HEADER_BYTES = 100,
    RECORD_HEADER_BYTES = 8,
    RECORD_CONTENT_BYTES =
      1 * BYTES_PER_INT +        /* int ShapeType = POLYGON */
      4 * BYTES_PER_DOUBLE +     /* double Box[ 4 ] = xMin,yMin,xMax,yMax */
      1 * BYTES_PER_INT +        /* int NumParts = 1 */
      1 * BYTES_PER_INT +        /* int NumPoints = 5 */
      1 * BYTES_PER_INT +        /* int Parts[ NumParts = 1 ] = 0 */
      VERTICES_PER_POLYGON * 2 * BYTES_PER_DOUBLE /* double [NumPoints*2]. */
  };
  unsigned char header[ HEADER_BYTES ] = "";
  unsigned char recordHeader[ RECORD_HEADER_BYTES ] = "";
  unsigned char recordContents[ RECORD_CONTENT_BYTES ] = "";
  const int records = rows * columns;
  const int shxFileBytes = HEADER_BYTES + records * RECORD_HEADER_BYTES;
  const int shpFileBytes = shxFileBytes + records * RECORD_CONTENT_BYTES;
  int byteIndex = 0;
  double xyRange[ 2 * 2 ] = { 0.0, 0.0, 0.0, 0.0 };
  enum { FILE_NAME_LENGTH = 255 };
  char shxFileName[ FILE_NAME_LENGTH + 1 ] = "";
  char shpFileName[ FILE_NAME_LENGTH + 1 ] = "";
  FILE* file = 0;
  float* lonlats = NEW( float, rows * columns * 2 ); /* Of grid cell centers.*/

  if ( lonlats ) {

    /* Construct file names of resulting files from base file name: */

    memset( shxFileName, 0, sizeof shxFileName );
    memset( shpFileName, 0, sizeof shpFileName );
    strncpy( shxFileName, fileName, FILE_NAME_LENGTH - 8 );
    strncpy( shpFileName, fileName, FILE_NAME_LENGTH - 8 );
    strcat( shxFileName, ".shx" );
    strcat( shpFileName, ".shp" );

    computeGridBounds( rows, columns, westEdge, southEdge,
                       cellWidth, cellHeight, unproject, xyRange );

    computeGridCellCenters( rows, columns, westEdge, southEdge,
                            cellWidth, cellHeight, unproject, lonlats );

    /* Initialize shx file header and records: */

    memset( header, 0, sizeof header );
    memset( recordHeader, 0, sizeof recordHeader );
    memset( recordContents, 0, sizeof recordContents );

    writeInt( header, 0, 9994, BIG );
    byteIndex = writeInt( header, 24, shxFileBytes / 2, BIG );
    byteIndex = writeInt( header, byteIndex, 1000, LITTLE );
    byteIndex = writeInt( header, byteIndex, POLYGON, LITTLE );
    byteIndex = writeDouble( header, byteIndex, xyRange[ 0 ], LITTLE );
    byteIndex = writeDouble( header, byteIndex, xyRange[ 1 ], LITTLE );
    byteIndex = writeDouble( header, byteIndex, xyRange[ 2 ], LITTLE );
    byteIndex = writeDouble( header, byteIndex, xyRange[ 3 ], LITTLE );

    writeInt( recordHeader, 0, HEADER_BYTES / 2, BIG );
    writeInt( recordHeader, 4, RECORD_CONTENT_BYTES / 2, BIG );

    writeInt( recordContents, 0, POLYGON, LITTLE );
    writeInt( recordContents, 36, PARTS_PER_POLYGON, LITTLE );
    writeInt( recordContents, 40, VERTICES_PER_POLYGON, LITTLE );

    /* Write shx file: */

    file = fopen( shxFileName, "wb" );

    if ( file ) {
      int record = 0;
      result = fwrite( header, HEADER_BYTES, 1, file ) == 1;

      for ( record = 0; result && record < records; ++record ) {
        const int offsetBytes =
          HEADER_BYTES +
          record * ( RECORD_HEADER_BYTES + RECORD_CONTENT_BYTES );
        writeInt( recordHeader, 0, offsetBytes / 2, BIG );
        result = fwrite( recordHeader, RECORD_HEADER_BYTES, 1, file ) == 1;
      }

      fclose( file ), file = 0;
    }

    /* Write shp file: */

    if ( result ) {
      file = fopen( shpFileName, "wb" );
      result = 0;

      if ( file ) {
        int record = 0;
        writeInt( header, 24, shpFileBytes / 2, BIG );
        result = fwrite( header, HEADER_BYTES, 1, file ) == 1;

        for ( record = 0; result && record < records; ++record ) {
          writeInt( recordHeader, 0, record + 1, BIG );
          result = fwrite( recordHeader, RECORD_HEADER_BYTES, 1, file ) == 1;

          /* Compute and write POLYGON record contents: */

          if ( result ) {
            const int row    = record / columns;
            const int column = record % columns;
            double xy[ VERTICES_PER_POLYGON * 2 ] = {
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            };

            computePolygonVertices( row, column,
                                    westEdge, southEdge, cellWidth, cellHeight,
                                    unproject, xy, xyRange );

            /* Write polygon bounds: */

            byteIndex = writeDouble( recordContents, 4, xyRange[ 0 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xyRange[ 1 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xyRange[ 2 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xyRange[ 3 ], LITTLE);

            /* Write polygon vertices: */

            byteIndex =
              writeDouble( recordContents, 48, xy[ 0 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 1 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 2 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 3 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 4 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 5 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 6 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 7 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 8 ], LITTLE );
            byteIndex =
              writeDouble( recordContents, byteIndex, xy[ 9 ], LITTLE );
            result =
              fwrite( recordContents, RECORD_CONTENT_BYTES, 1, file ) == 1;
          }
        }

        fclose( file ), file = 0;
      }
    }

    if ( ! result ) {
      perror( "\n\nFailed because" );
    } else {
      result = writePRJFile( fileName, 0 );

      if ( AND3( result, variable, data ) ) {
        result =
          writeDataToDBFFile( fileName, variable, units,
                              timesteps, yyyymmddhh, timestepType,
                              records, type, components, data, lonlats,
                              0, 0, 0, 0, 0 );
      }
    }

    FREE( lonlats );
  }

  return result;
}



/******************************************************************************
PURPOSE: writePointsToShapefile - Write points as a lon-lat Shapefile Point
         file set (shp, shx, dbf, prj) and a csv file containing time-varying
         data.
INPUTS:  const char* const fileName      Base name of file to create. "storet".
         const char* const variableName  Name of variable. E.g., "salinity".
         const char* const units         Units of variable. E.g., "PSU".
         const int timesteps             Number of timesteps of data.
         const int hoursPerTimestep      Number of hours per timestep. 1, 24.
         const int yyyymmddhh[ timesteps ]  Timestamps of data or 0 if no data.
         const int count                    Number of points per timestep.
         const float lonlats[ count * 2 ]   Lon-lats.
         const float z[ count ]             Optional: z coordinate of points.
         const int components               1 = scalar, 2, 3 = vector.
         const float data[ components ][ timesteps * count ]
                                         Optional: component data per point.
         const char* const sids[ count ]   Station id strings or 0.
         const int ids[ count ]            Station ids or 0.
         const char* const metadata[ count ]   Station metadata or 0.
         const int writeCSV                Write csv file?
RETURNS: int 1 if successful, else 0 and failure message is printed to stderr.
NOTES:   Creates files fileName.shp and fileName.shx
         then calls routines that create fileName.dbf and fileName.csv.
         See 1998 ESRI Shapefile Specification pages 2, 4, 5, 15, 23, 24.
         http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
         This routine does not use shapelib API since it was written before
         that library was used in this project and also because this routine
         is more straight-forward and much faster than using the shapelib API.
******************************************************************************/

void writePointsToShapefile( const char* const fileName,
                             const char* const variableName,
                             const char* const units,
                             const int timesteps,
                             const int hoursPerTimestep,
                             const int yyyymmddhh[],
                             const int count,
                             const float lonlats[],
                             const float z[],
                             const int components,
                             const float data[],
                             const char* const sids[],
                             const int ids[],
                             const char* const metadata[],
                             const int writeCSV ) {

  PRE013( fileName,
          *fileName,
          IMPLIES( data, AND3( variableName, *variableName, units ) ),
          IMPLIES( data, timesteps > 0 ),
          IMPLIES( data, hoursPerTimestep > 0 ),
          IMPLIES( data, yyyymmddhh ),
          count > 0,
          lonlats,
          is_valid_longitude_latitude( lonlats[ 0 ], lonlats[ 1 ] ),
          is_valid_longitude_latitude( lonlats[  count * 2 - 2 ],
                                       lonlats[  count * 2 - 1 ] ),
          IMPLIES( z, AND2( ! is_nan( z[ 0 ] ), ! is_nan( z[ count - 1 ] ) ) ),
          IMPLIES( timesteps > 1, AND2( components > 0, data ) ),
          IMPLIES( data, AND2( ! is_nan( data[ 0 ] ),
                               ! is_nan( data[ timesteps * count - 1 ] ) ) ) );

  int result = 0;
  enum { /* Constants from the Shapefile Spec (see above link): */
    BYTES_PER_INT = 4,
    BYTES_PER_DOUBLE = 8,
    POINTZ = 11,
    HEADER_BYTES = 100,
    RECORD_HEADER_BYTES = 8,
    RECORD_CONTENT_BYTES =
      BYTES_PER_INT +     /* int ShapeType = POINTZ */
      4 * BYTES_PER_DOUBLE  /* x, y, z, m */
  };
  unsigned char header[ HEADER_BYTES ] = "";
  unsigned char recordHeader[ RECORD_HEADER_BYTES ] = "";
  unsigned char recordContents[ RECORD_CONTENT_BYTES ] = "";
  const int records = count;
  const int shxFileBytes = HEADER_BYTES + records * RECORD_HEADER_BYTES;
  const int shpFileBytes = shxFileBytes + records * RECORD_CONTENT_BYTES;
  int byteIndex = 0;
  double xyRange[ 2 * 2 ] = { 0.0, 0.0, 0.0, 0.0 };
  enum { FILE_NAME_LENGTH = 255 };
  char shxFileName[ FILE_NAME_LENGTH + 1 ] = "";
  char shpFileName[ FILE_NAME_LENGTH + 1 ] = "";
  FILE* file = 0;

  /* Construct file names of resulting files from base file name: */

  memset( shxFileName, 0, sizeof shxFileName );
  memset( shpFileName, 0, sizeof shpFileName );
  strncpy( shxFileName, fileName, FILE_NAME_LENGTH - 8 );
  strncpy( shpFileName, fileName, FILE_NAME_LENGTH - 8 );
  strcat( shxFileName, ".shx" );
  strcat( shpFileName, ".shp" );

  computePointBounds( count, lonlats, xyRange );

  /* Initialize shx file header and records: */

  memset( header, 0, sizeof header );
  memset( recordHeader, 0, sizeof recordHeader );
  memset( recordContents, 0, sizeof recordContents );

  writeInt( header, 0, 9994, BIG );
  byteIndex = writeInt( header, 24, shxFileBytes / 2, BIG );
  byteIndex = writeInt( header, byteIndex, 1000, LITTLE );
  byteIndex = writeInt( header, byteIndex, POINTZ, LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 0 ], LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 1 ], LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 2 ], LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 3 ], LITTLE );
  
  if ( z ) { /* Compute and store z range (else defaults to [0, 0]): */
    double minimum = z[ 0 ];
    double maximum = minimum;
    int index = 0;

    for ( index = 1; index < count; ++index ) {
      const double value = z[ index ];

      if ( value < minimum ) {
        minimum = value;
      } else if ( value > maximum ) {
        maximum = value;
      }
    }

    byteIndex = writeDouble( header, byteIndex, minimum, LITTLE );
    byteIndex = writeDouble( header, byteIndex, maximum, LITTLE );
  }

  writeInt( recordHeader, 0, HEADER_BYTES / 2, BIG );
  writeInt( recordHeader, 4, RECORD_CONTENT_BYTES / 2, BIG );

  writeInt( recordContents, 0, POINTZ, LITTLE );

  /* Write shx file: */

  file = fopen( shxFileName, "wb" );

  if ( file ) {
    int record = 0;
    result = fwrite( header, HEADER_BYTES, 1, file ) == 1;

    for ( record = 0; result && record < records; ++record ) {
      const int offsetBytes =
        HEADER_BYTES +
        record * ( RECORD_HEADER_BYTES + RECORD_CONTENT_BYTES );
      writeInt( recordHeader, 0, offsetBytes / 2, BIG );
      result = fwrite( recordHeader, RECORD_HEADER_BYTES, 1, file ) == 1;
    }

    fclose( file ), file = 0;
  }

  /* Write shp file: */

  if ( result ) {
    file = fopen( shpFileName, "wb" );

    if ( file ) {
      int record = 0;
      writeInt( header, 24, shpFileBytes / 2, BIG );
      result = fwrite( header, HEADER_BYTES, 1, file ) == 1;

      for ( record = 0; result && record < records; ++record ) {
        const int record2 = record + record;
        writeInt( recordHeader, 0, record + 1, BIG );
        result = fwrite( recordHeader, RECORD_HEADER_BYTES, 1, file ) == 1;

        /* Write POINTZ record contents: */

        if ( result ) {
          const double longitude = lonlats[ record2 ];
          const double latitude = lonlats[ record2 + 1 ];
          byteIndex = writeDouble( recordContents, 4, longitude, LITTLE );
          byteIndex = writeDouble( recordContents, 12, latitude, LITTLE );

          if ( z ) {
            byteIndex = writeDouble( recordContents, 20, z[ record ], LITTLE );
          }

          result =
            fwrite( recordContents, RECORD_CONTENT_BYTES, 1, file ) == 1;
        }
      }

      fclose( file ), file = 0;
    }
  }

  if ( ! result ) {
    perror( "\n\nFailed because" );
  } else {
    result = writePRJFile( fileName, 0 );

    if ( AND2( result, data ) ) {
      const int timestepType = hoursPerTimestep == 24 ? DAILY : HOURLY;
      result =
        writeDataToDBFFile( fileName, variableName, units,
                            timesteps, yyyymmddhh, timestepType,
                            records, FLOAT_TYPE,
                            components, data,
                            lonlats, z, sids, ids, metadata, writeCSV );
    }
  }
}



/******************************************************************************
PURPOSE: writePolylinesToShapefile - Write a polyline coordinates as a lon-lat
         Shapefile Polyline file set (shp, shx, dbf, prj).
INPUTS:  const char* baseFileName  Base name of file to create. "edm_bounds".
         const int polylineCount   Number of polylines.
         const int vertexCount     Number of vertices.
         const int counts[ polylineCount ]   Number of vertices per polyline.
         const float lonlats[ vertexCount * 2 ]  Lon-lat vertices.
RETURNS: int 1 if successful, else 0 and a message is printed to stderr.
NOTES:   Creates files fileName.shp and fileName.shx.
         See 1998 ESRI Shapefile Specification pages 2, 4, 5, 16, 23, 24.
         http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf
         This routine does not use shapelib API since it was written before
         that library was used in this project and also because this routine
         is more straight-forward and much faster than using the shapelib API.
******************************************************************************/

int writePolylinesToShapefile( const char* baseFileName,
                               const int polylineCount,
                               const int vertexCount,
                               const int counts[],
                               const float lonlats[] ) {
  int result = 0;
  enum { /* Constants from the Shapefile Spec (see above link): */
    BYTES_PER_INT = 4,
    BYTES_PER_DOUBLE = 8,
    POLYLINE = 3,                 /* ShapeType = POLYLINE. */
    HEADER_BYTES = 100,
    RECORD_HEADER_BYTES = 8,
    RECORD_CONTENT_BYTES =
      1 * BYTES_PER_INT +        /* int ShapeType = POLYLINE */
      4 * BYTES_PER_DOUBLE +     /* double Box[ 4 ] = xMin,yMin,xMax,yMax */
      1 * BYTES_PER_INT +        /* int NumParts */
      1 * BYTES_PER_INT          /* int NumPoints */
  };
  const int recordContentBytes =
    RECORD_CONTENT_BYTES +
    polylineCount * BYTES_PER_INT + vertexCount * 2 * BYTES_PER_DOUBLE;
  unsigned char header[ HEADER_BYTES ] = "";
  unsigned char recordHeader[ RECORD_HEADER_BYTES ] = "";
  unsigned char recordContents[ RECORD_CONTENT_BYTES ] = "";
  const int shxFileBytes = HEADER_BYTES + RECORD_HEADER_BYTES;
  const int shpFileBytes =
    shxFileBytes + RECORD_CONTENT_BYTES +
    polylineCount * BYTES_PER_INT +
    vertexCount * 2 * BYTES_PER_DOUBLE;
  int byteIndex = 0;
  double xyRange[ 2 * 2 ] = { 0.0, 0.0, 0.0, 0.0 };
  enum { FILE_NAME_LENGTH = 255 };
  char shxFileName[ FILE_NAME_LENGTH + 1 ] = "";
  char shpFileName[ FILE_NAME_LENGTH + 1 ] = "";
  FILE* file = 0;

  DEBUG( fprintf( stderr, "baseFileName = %s\n", baseFileName ); )

  /* Construct file names of resulting files from base file name: */

  memset( shxFileName, 0, sizeof shxFileName );
  memset( shpFileName, 0, sizeof shpFileName );
  strncpy( shxFileName, baseFileName, FILE_NAME_LENGTH - 4 );
  strncpy( shpFileName, baseFileName, FILE_NAME_LENGTH - 4 );
  strcat( shxFileName, ".shx" );
  strcat( shpFileName, ".shp" );

  computeVertexBounds( vertexCount, lonlats, xyRange );

  /* Initialize shx file header and records: */

  memset( header, 0, sizeof header );
  memset( recordHeader, 0, sizeof recordHeader );
  memset( recordContents, 0, sizeof recordContents );

  writeInt( header, 0, 9994, BIG );
  byteIndex = writeInt( header, 24, shxFileBytes / 2, BIG );
  byteIndex = writeInt( header, byteIndex, 1000, LITTLE );
  byteIndex = writeInt( header, byteIndex, POLYLINE, LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 0 ], LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 1 ], LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 2 ], LITTLE );
  byteIndex = writeDouble( header, byteIndex, xyRange[ 3 ], LITTLE );

  writeInt( recordHeader, 0, HEADER_BYTES / 2, BIG );
  writeInt( recordHeader, 4, recordContentBytes / 2, BIG );

  writeInt( recordContents, 0, POLYLINE, LITTLE );
  writeInt( recordContents, 36, polylineCount, LITTLE );
  writeInt( recordContents, 40, vertexCount, LITTLE );

  /* Write shx file: */

  file = fopen( shxFileName, "wb" );

  if ( file ) {
    result = fwrite( header, HEADER_BYTES, 1, file ) == 1;
    writeInt( recordHeader, 0, HEADER_BYTES / 2, BIG );
    result = fwrite( recordHeader, RECORD_HEADER_BYTES, 1, file ) == 1;
    fclose( file ), file = 0;
  }

  /* Write shp file: */

  if ( result ) {
    file = fopen( shpFileName, "wb" );

    if ( file ) {
      writeInt( header, 24, shpFileBytes / 2, BIG );
      result = fwrite( header, HEADER_BYTES, 1, file ) == 1;
      writeInt( recordHeader, 0, 1, BIG );
      result = fwrite( recordHeader, RECORD_HEADER_BYTES, 1, file ) == 1;

      /* Compute and write POLYLINE record contents: */

      if ( result ) {
        byteIndex = writeDouble( recordContents, 4, xyRange[ 0 ], LITTLE );
        byteIndex =
          writeDouble( recordContents, byteIndex, xyRange[ 1 ], LITTLE );
        byteIndex =
          writeDouble( recordContents, byteIndex, xyRange[ 2 ], LITTLE );
        byteIndex =
          writeDouble( recordContents, byteIndex, xyRange[ 3 ], LITTLE);
        byteIndex =
          writeInt( recordContents, byteIndex, polylineCount, LITTLE );
        byteIndex =
          writeInt( recordContents, byteIndex, vertexCount, LITTLE );
        result = fwrite( recordContents, RECORD_CONTENT_BYTES, 1, file ) == 1;

        /* Write Parts array: */

        if ( result ) {
          int index = 0;
          int polyline = 0;

          for ( polyline = 0; result && polyline < polylineCount; ++polyline) {
            unsigned char value[ BYTES_PER_INT ] = "";
            writeInt( value, 0, index, LITTLE );
            result = fwrite( value, BYTES_PER_INT, 1, file ) == 1;
            index += counts[ polyline ];
          }
        }

        /* Write Points array: */

        if ( result ) {
          const float* longitudes = lonlats;
          const float* latitudes  = lonlats + 1;
          int vertex = 0;

          for ( vertex = 0; result && vertex < vertexCount;
                ++vertex, longitudes += 2, latitudes += 2 ) {
            unsigned char value[ 2 * BYTES_PER_DOUBLE ] = "";
            byteIndex = writeDouble( value, 0, *longitudes, LITTLE );
            writeDouble( value, byteIndex, *latitudes, LITTLE );
            DEBUG2( fprintf( stderr, "(%g, %g)\n", *longitudes, *latitudes );)
            result = fwrite( value, 2 * BYTES_PER_DOUBLE, 1, file ) == 1;
          }
        }
      }

      fclose( file ), file = 0;
    }
  }

  if ( ! result ) {
    perror( "\n\nFailed because" );
  } else {
    result = writePRJFile( baseFileName, 0 );
  }

  return result;
}



/******************************************************************************
PURPOSE: computePolygonVertices - Compute vertices of grid cell as an
         explicitly closed 2D 5-vertex polygon ring in clockwise order.
INPUTS:  const int row             0-based row index of grid cell.
         const int column          0-based column index of grid cell.
         const double westEdge     Distance from origin to west edge of grid.
         const double southEdge    Distance from origin to south edge of ".
         const double cellWidth    Width of each grid cell (e.g., 12000 m).
         const double cellWHeight  Height of each grid cell (e.g., 12000 m).
OUTPUTS: double xy[ 5 * 2 ]        Sequence of polygon (x, y) vertices.
         double xyRange[ 2 * 2 ]   xMin, yMin, xMax, yMax.
******************************************************************************/

void computePolygonVertices( const int row, const int column,
                             const double westEdge,
                             const double southEdge,
                             const double cellWidth,
                             const double cellHeight,
                             Unproject unproject,
                             double xy[ 5 * 2 ],
                             double xyRange[ 2 * 2 ] ) {

  double x = westEdge + column * cellWidth;
  double y = southEdge + row * cellHeight;
  double px = x;
  double py = y;
  double swap = 0.0;

  if ( unproject ) {
    unproject( x, y, &px, &py );
  }

  xy[ 0 ] = px;
  xy[ 1 ] = py;

  y += cellHeight;

  if ( unproject ) {
    unproject( x, y, &px, &py );
  } else {
    px = x;
    py = y;
  }

  xy[ 2 ] = px;
  xy[ 3 ] = py;

  x += cellWidth;

  if ( unproject ) {
    unproject( x, y, &px, &py );
  } else {
    px = x;
    py = y;
  }

  xy[ 4 ] = px;
  xy[ 5 ] = py;

  y -= cellHeight;

  if ( unproject ) {
    unproject( x, y, &px, &py );
  } else {
    px = x;
    py = y;
  }

  xy[ 6 ] = px;
  xy[ 7 ] = py;

  xy[ 8 ] = xy[ 0 ];
  xy[ 9 ] = xy[ 1 ];

  xyRange[ 0 ] = xyRange[ 1 ] = xyRange[ 2 ] = xyRange[ 3 ] = 0.0;
  computeRange( xy, 10, 2, xyRange );
  computeRange( xy + 1, 10, 2, xyRange + 2 );
  swap = xyRange[ 1 ];
  xyRange[ 1 ] = xyRange[ 2 ];
  xyRange[ 2 ] = swap;
}



/******************************************************************************
PURPOSE: computeGridCellVertices - Compute lon-lat coordinates of rectangular
         grid cell corners.
INPUTS:  const int rows                            Number of grid row cells.
         const int columns                         Number of grid column cells.
         const float longitudes[ rows * columns ]  Longitude of cell center.
         const float latitudes[  rows * columns ]  Latitude of cell center.
OUTPUTS: float vertices[ ( rows + 1 ) * ( columns + 1 ) * 2 ]
         Corners in grid cell order with interleaved lon-lat coordinates.
NOTES:   Uses linear interpolation and extrapolation to the edges.
******************************************************************************/

void computeGridCellVertices( const int rows, const int columns,
                              const float longitudes[],
                              const float latitudes[], float vertices[] ) {

  PRE08( rows > 1, columns > 1, rows * columns >= 4,
         longitudes, latitudes, vertices,
         is_valid_longitude_latitude( longitudes[ 0 ], latitudes[ 0 ] ),
         is_valid_longitude_latitude( longitudes[ rows * columns - 1 ],
                                      latitudes[ rows * columns - 1 ] ) );

  const int rows_1    = rows - 1;
  const int columns_1 = columns - 1;
  const int columnsPlus1 = columns + 1;
  const int columnsPlus1Times2 = columnsPlus1 + columnsPlus1;
  const int count = ( rows + 1 ) * ( columns + 1 ) * 2;
  int row    = 0;
  int column = 0;
  int index = 0;
  int vIndex = 0;

  /*
   * First compute linearly interpolated corners of all interior cells:
   * Note: rows increase north to south and columns increase west to east.
   */

/* TEMP HACK #pragma omp parallel for private( column ) */

  for ( row = 0; row < rows_1; ++row ) {
    const int nextRow        = row + 1;
    const int rowOffset      = row     * columns;
    const int nextRowOffset  = nextRow * columns;
    const int verticesOffset = nextRow * columnsPlus1Times2 + 2;

    /* Interior row, interior columns: */

    for ( column = 0; column < columns_1; ++column ) {
      const int nextColumn             = column + 1;
      const int verticesIndex          = verticesOffset + column + column;
      const int dataIndex              = rowOffset + column;
      const int nextColumnIndex        = dataIndex + 1;
      const int nextRowIndex           = nextRowOffset + column;
      const int nextRowNextColumnIndex = nextRowOffset + nextColumn;

      const float longitude                  = longitudes[ dataIndex ];
      const float nextColumnLongitude        = longitudes[ nextColumnIndex ];
      const float nextRowLongitude           = longitudes[ nextRowIndex ];
      const float nextRowNextColumnLongitude =
        longitudes[ nextRowNextColumnIndex ];

      const float latitude                  = latitudes[ dataIndex ];
      const float nextColumnLatitude        = latitudes[ nextColumnIndex ];
      const float nextRowLatitude           = latitudes[ nextRowIndex ];
      const float nextRowNextColumnLatitude =
        latitudes[ nextRowNextColumnIndex ];
      const float interpolatedLongitude = 0.25 *
        ( longitude + nextColumnLongitude +
          nextRowLongitude + nextRowNextColumnLongitude );
      const float interpolatedLatitude = 0.25 *
        ( latitude + nextColumnLatitude +
          nextRowLatitude + nextRowNextColumnLatitude );

      vertices[ verticesIndex     ] = interpolatedLongitude;
      vertices[ verticesIndex + 1 ] = interpolatedLatitude;

    } /* End loop on interior columns. */

  } /* End parallel loop on interior rows. */

  /* Serial region (not worth parallelizing): */

  /* Last row, interior columns (extrapolated top edge, except corners): */

  for ( column = 1, index = rows_1 * columns + 1,
        vIndex = rows_1 * columnsPlus1Times2 + 2;
        column < columns; ++column, vIndex += 2, ++index ) {
    const int previousColumnIndex       = index - 1;
    const int extrapolatedIndex         = vIndex + columnsPlus1Times2;
    const float longitude               = longitudes[ index ];
    const float latitude                = latitudes[ index ];
    const float previousColumnLongitude = longitudes[ previousColumnIndex ];
    const float previousColumnLatitude  = latitudes[ previousColumnIndex ];
    const float midpointLongitude       =
      0.5 * ( longitude + previousColumnLongitude );
    const float midpointLatitude        =
      0.5 * ( latitude  + previousColumnLatitude );
    const float extrapolatedInputLongitude = vertices[ vIndex ];
    const float extrapolatedInputLatitude  = vertices[ vIndex + 1 ];
    const float longitudeDifference =
      midpointLongitude - extrapolatedInputLongitude;
    const float latitudeDifference  =
      midpointLatitude - extrapolatedInputLatitude;
    const float extrapolatedLongitude = midpointLongitude +longitudeDifference;
    const float extrapolatedLatitude = midpointLatitude + latitudeDifference;
    vertices[ extrapolatedIndex     ] = extrapolatedLongitude;
    vertices[ extrapolatedIndex + 1 ] = extrapolatedLatitude;
  }

  /* First row, interior columns (extrapolated bottom edge, except corners): */

  for ( column = 1, index = 1, vIndex = columnsPlus1Times2 + 2;
        column < columns; ++column, vIndex += 2, ++index ) {
    const int previousColumnIndex       = index - 1;
    const int extrapolatedIndex         = vIndex - columnsPlus1Times2;
    const float longitude               = longitudes[ index ];
    const float latitude                = latitudes[ index ];
    const float previousColumnLongitude = longitudes[ previousColumnIndex ];
    const float previousColumnLatitude  = latitudes[ previousColumnIndex ];
    const float midpointLongitude       =
      0.5 * ( longitude + previousColumnLongitude );
    const float midpointLatitude        =
      0.5 * ( latitude  + previousColumnLatitude );
    const float extrapolatedInputLongitude = vertices[ vIndex ];
    const float extrapolatedInputLatitude  = vertices[ vIndex + 1 ];
    const float longitudeDifference =
      midpointLongitude - extrapolatedInputLongitude;
    const float latitudeDifference  =
      midpointLatitude - extrapolatedInputLatitude;
    const float extrapolatedLongitude = midpointLongitude +longitudeDifference;
    const float extrapolatedLatitude = midpointLatitude + latitudeDifference;
    vertices[ extrapolatedIndex     ] = extrapolatedLongitude;
    vertices[ extrapolatedIndex + 1 ] = extrapolatedLatitude;
  }

  /* First column, interior rows (extrapolated left edge, except corners): */

  for ( row = 1, index = columns, vIndex = columnsPlus1Times2 + 2;
        row < rows; ++row, vIndex += columnsPlus1Times2, index += columns ) {
    const int previousRowIndex          = index - columns;
    const int extrapolatedIndex         = vIndex - 2;
    const float longitude               = longitudes[ index ];
    const float latitude                = latitudes[ index ];
    const float previousRowLongitude    = longitudes[ previousRowIndex ];
    const float previousRowLatitude     = latitudes[ previousRowIndex ];
    const float midpointLongitude       =
      0.5 * ( longitude + previousRowLongitude );
    const float midpointLatitude        =
      0.5 * ( latitude  + previousRowLatitude );
    const float extrapolatedInputLongitude = vertices[ vIndex ];
    const float extrapolatedInputLatitude  = vertices[ vIndex + 1 ];
    const float longitudeDifference =
      midpointLongitude - extrapolatedInputLongitude;
    const float latitudeDifference  =
      midpointLatitude - extrapolatedInputLatitude;
    const float extrapolatedLongitude = midpointLongitude +longitudeDifference;
    const float extrapolatedLatitude = midpointLatitude + latitudeDifference;
    vertices[ extrapolatedIndex     ] = extrapolatedLongitude;
    vertices[ extrapolatedIndex + 1 ] = extrapolatedLatitude;
  }

  /* Last column, interior rows (extrapolated right edge, except corners): */

  for ( row = 1, index = columns + columns - 1,
        vIndex = columnsPlus1Times2 + columnsPlus1Times2 - 4;
        row < rows; ++row, vIndex += columnsPlus1Times2, index += columns ) {
    const int previousRowIndex          = index - columns;
    const int extrapolatedIndex         = vIndex + 2;
    const float longitude               = longitudes[ index ];
    const float latitude                = latitudes[ index ];
    const float previousRowLongitude    = longitudes[ previousRowIndex ];
    const float previousRowLatitude     = latitudes[ previousRowIndex ];
    const float midpointLongitude       =
      0.5 * ( longitude + previousRowLongitude );
    const float midpointLatitude        =
      0.5 * ( latitude  + previousRowLatitude );
    const float extrapolatedInputLongitude = vertices[ vIndex ];
    const float extrapolatedInputLatitude  = vertices[ vIndex + 1 ];
    const float longitudeDifference =
      midpointLongitude - extrapolatedInputLongitude;
    const float latitudeDifference  =
      midpointLatitude - extrapolatedInputLatitude;
    const float extrapolatedLongitude = midpointLongitude +longitudeDifference;
    const float extrapolatedLatitude = midpointLatitude + latitudeDifference;
    vertices[ extrapolatedIndex     ] = extrapolatedLongitude;
    vertices[ extrapolatedIndex + 1 ] = extrapolatedLatitude;
  }

  /* First row, first column cell (extrapolated bottom-left corner): */

  vIndex = columnsPlus1Times2 + 2;

  {
    const float longitude             = longitudes[ 0 ];
    const float latitude              = latitudes[ 0 ];
    const float diagonalLongitude     = vertices[ vIndex ];
    const float diagonalLatitude      = vertices[ vIndex + 1 ];
    const float longitudeDifference   = longitude - diagonalLongitude;
    const float latitudeDifference    = latitude  - diagonalLatitude;
    const float extrapolatedLongitude = longitude + longitudeDifference;
    const float extrapolatedLatitude  = latitude  + latitudeDifference;
    vertices[ 0 ]                     = extrapolatedLongitude;
    vertices[ 1 ]                     = extrapolatedLatitude;
  }

  /* First row, last column cell (extrapolated bottom-right corner): */

  vIndex = columnsPlus1Times2 + columnsPlus1Times2 - 4;

  {
    const int extrapolatedIndex       = columnsPlus1Times2 - 2;
    const int dataIndex               = columns - 1;
    const float longitude             = longitudes[ dataIndex ];
    const float latitude              = latitudes[ dataIndex ];
    const float diagonalLongitude     = vertices[ vIndex ];
    const float diagonalLatitude      = vertices[ vIndex + 1 ];
    const float longitudeDifference   = longitude - diagonalLongitude;
    const float latitudeDifference    = latitude  - diagonalLatitude;
    const float extrapolatedLongitude = longitude + longitudeDifference;
    const float extrapolatedLatitude  = latitude  + latitudeDifference;
    vertices[ extrapolatedIndex ]     = extrapolatedLongitude;
    vertices[ extrapolatedIndex + 1 ] = extrapolatedLatitude;
  }

  /* Last row, first column cell (extrapolated top-left corner): */

  vIndex = rows_1 * columnsPlus1Times2 + 2;

  {
    const int extrapolatedIndex       = rows * columnsPlus1Times2;
    const int dataIndex               = rows_1 * columns;
    const float longitude             = longitudes[ dataIndex ];
    const float latitude              = latitudes[ dataIndex ];
    const float diagonalLongitude     = vertices[ vIndex ];
    const float diagonalLatitude      = vertices[ vIndex + 1 ];
    const float longitudeDifference   = longitude - diagonalLongitude;
    const float latitudeDifference    = latitude  - diagonalLatitude;
    const float extrapolatedLongitude = longitude + longitudeDifference;
    const float extrapolatedLatitude  = latitude  + latitudeDifference;
    vertices[ extrapolatedIndex     ] = extrapolatedLongitude;
    vertices[ extrapolatedIndex + 1 ] = extrapolatedLatitude;
  }

  /* Last row, last column cell (extrapolated top-right corner): */

  vIndex = rows * columnsPlus1Times2 - 4;

  {
    const int extrapolatedIndex       = vIndex + columnsPlus1Times2 + 2;
    const int dataIndex               = rows * columns - 1;
    const float longitude             = longitudes[ dataIndex ];
    const float latitude              = latitudes[ dataIndex ];
    const float diagonalLongitude     = vertices[ vIndex ];
    const float diagonalLatitude      = vertices[ vIndex + 1 ];
    const float longitudeDifference   = longitude - diagonalLongitude;
    const float latitudeDifference    = latitude  - diagonalLatitude;
    const float extrapolatedLongitude = longitude + longitudeDifference;
    const float extrapolatedLatitude  = latitude  + latitudeDifference;
    vertices[ extrapolatedIndex ]     = extrapolatedLongitude;
    vertices[ extrapolatedIndex + 1]  = extrapolatedLatitude;
  }

  /* Clamp any out-of-range values: */

/* TEMP HACK #pragma omp parallel for */

  for ( vIndex = 0; vIndex < count; vIndex += 2 ) {
    const int vIndex1 = vIndex + 1;
    vertices[ vIndex  ] = CLAMPED_TO_RANGE( vertices[ vIndex ], -180.0, 180.0);
    vertices[ vIndex1 ] = CLAMPED_TO_RANGE( vertices[ vIndex1 ], -90.0, 90.0 );
  }

  POST02( is_valid_longitude_latitude( vertices[ 0 ], vertices[ 1 ] ),
          is_valid_longitude_latitude( vertices[(rows+1)*(columns+1)* 2 - 2 ],
                                       vertices[(rows+1)*(columns+1)* 2 - 1]));
}




/******************************************************************************
PURPOSE: printShape - Print a ESRI polygon shape for debugging purposes.
INPUTS:  const SHPObject* shape  Shape to print.
NOTES:   http://shapelib.maptools.org/shp_api.html
******************************************************************************/

void printShape( const SHPObject* shape ) {
  PRE02( shape, shape->nVertices >= 2 );
  int index = 0;
  fprintf( stderr, "nSHPType = %d\n", shape->nSHPType );
  fprintf( stderr, "nShapeId = %d\n", shape->nShapeId );
  fprintf( stderr, "nParts = %d\n", shape->nParts );
  fprintf( stderr, "nVertices = %d\n", shape->nVertices );

  fprintf( stderr, "panPartStart:\n" );

  for ( index = 0; index < shape->nParts; ++index ) {
    fprintf( stderr, "%2d %d\n", index, shape->panPartStart[ index ] );
  }

  fprintf( stderr, "\n" );

  fprintf( stderr, "panPartType/Start/Count/X/Y:\n" );

  for ( index = 0; index < shape->nParts; ++index ) {
    const int start = shape->panPartStart[ index ];
    const int count =
      shape->nParts == 1 ? shape->nVertices
      : index < shape->nParts - 1 ?
        shape->panPartStart[ index + 1 ] - shape->panPartStart[ index ]
      : shape->nVertices - shape->panPartStart[ index ];
    const double x0 = shape->padfX[ start ];
    const double y0 = shape->padfY[ start ];
    const double x1 = shape->padfX[ start + 1 ];
    const double y1 = shape->padfY[ start + 1 ];
    const double xn_2 = shape->padfX[ start + count - 2 ];
    const double yn_2 = shape->padfY[ start + count - 2 ];
    const double xn_1 = shape->padfX[ start + count - 1 ];
    const double yn_1 = shape->padfY[ start + count - 1 ];
    fprintf( stderr,
             "# %3d: <%d, %4d, #%4d (%lg, %lg),(%lg, %lg)..."
             "(%lg, %lg), (%lg, %lg)>\n",
             index,
             shape->panPartType[ index ], shape->panPartStart[ index ],
             count,
             x0, y0, x1, y1, xn_2, yn_2, xn_1, yn_1 );
    CHECK( shape->panPartStart[ index ] >= 0 );
    CHECK( shape->panPartStart[ index ] < shape->nVertices );
    CHECK( shape->panPartType[ index ] == SHPP_RING );
    CHECK( IMPLIES( IN3( shape->nSHPType, SHPT_POLYGON, SHPT_POLYGONZ ),
                    AND3( count >= 4, x0 == xn_1, y0 == yn_1 ) ) );
  }
}



/******************************************************************************
PURPOSE: printPolygon - Print a polygon for tracing/debugging purposes.
INPUTS:  const gpc_polygon* polygon  Polygon to print.
NOTES:   http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
******************************************************************************/

void printPolygon( const gpc_polygon* polygon ) {
  PRE0 ( polygon );
  int index = 0;
  fprintf( stderr, "num_contours = %d\n", polygon->num_contours );
  fprintf( stderr, "  hole[] =" );

  for ( index = 0; index < polygon->num_contours; ++index ) {
    fprintf( stderr, " %d", polygon->hole[ index ] );
  }

  fprintf( stderr, "\n  contour[]:\n" );

  for ( index = 0; index < polygon->num_contours; ++index ) {
    int v = 0;
    fprintf( stderr, "    num_vertices = %d:\n",
             polygon->contour[ index ].num_vertices );

    for ( v = 0; v < polygon->contour[ index ].num_vertices; ++v ) {

      if ( v < 2 || v > polygon->contour[ index ].num_vertices - 3 ) {
        fprintf( stderr, " #%4d (%lf, %lf)",
                 v,
                 polygon->contour[ index ].vertex[ v ].x,
                 polygon->contour[ index ].vertex[ v ].y );
      }
    }

    fprintf( stderr, "\n" );
  }
}



/******************************************************************************
PURPOSE: printTriangles - Print a triangle strip for debugging purposes.
INPUTS:  const gpc_tristrip* tristrip  Tristrip to print.
NOTES:   http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
******************************************************************************/

void printTriangles( const gpc_tristrip* tristrip ) {
  PRE0 ( tristrip );
  int index = 0;
  fprintf( stderr, "num_strips = %d\n", tristrip->num_strips );

  fprintf( stderr, "\n  strip[]:\n" );

  for ( index = 0; index < tristrip->num_strips; ++index ) {
    const int vertexCount = tristrip->strip[ index ].num_vertices;
    int v = 0;
    fprintf( stderr, "    num_vertices = %d:\n", vertexCount );

    for ( v = 0; v < vertexCount; ++v ) {

      if ( v < 5 || v > vertexCount - 6 ) {
        fprintf( stderr, " #%4d (%lg, %lg)",
                 v,
                 tristrip->strip[ index ].vertex[ v ].x,
                 tristrip->strip[ index ].vertex[ v ].y );
      }
    }

    fprintf( stderr, "\n" );
  }
}



/******************************************************************************
PURPOSE: deallocatePolygons - Deallocate polygons.
OUTPUTS: int count                  Number of polygons.
         PolygonShape polygons[ count ]  Polygons to deallocate.
******************************************************************************/

void deallocatePolygons( int count, PolygonShape polygons[] ) {
  PRE02( count >= 0, IMPLIES( count, polygons ) );

  while ( count-- ) {

    if ( polygons[ count ].polygon.num_contours > 0 ) {
      gpc_free_polygon( &polygons[ count ].polygon );
      polygons[ count ].polygon.num_contours = 0;
      polygons[ count ].polygon.hole = 0;
      polygons[ count ].polygon.contour = 0;
    }

    if ( polygons[ count ].triangles.num_strips > 0 ) {
      gpc_free_tristrip( &polygons[ count ].triangles );
      polygons[ count ].triangles.num_strips = 0;
      polygons[ count ].triangles.strip = 0;
    }

    polygons[ count ].id = 0;
  }

  FREE( polygons );
  POST0( polygons == 0 );
}



/******************************************************************************
PURPOSE: shapefileType - Get type of shapefile: SHPT_POLYGON, SHPT_ARC, etc.
INPUTS:  const char* baseFileName   Base name (no ext) of Shapefile to read.
RETURNS: int type if successful, else 0 and a failure message is printed.
******************************************************************************/

int shapefileType( const char* baseFileName ) {
  PRE02( baseFileName, *baseFileName );
  int result = 0;
  SHPHandle handle = SHPOpen( baseFileName, "rb" );

  if ( ! handle ) {
    fprintf( stderr, "\nFailed to open Shapefile '%s'\n", baseFileName );
    perror( "because" );
  } else {
    int type = 0;
    int shapes = 0;
    SHPGetInfo( handle, &shapes, &type, 0, 0 );
    SHPClose( handle ), handle = 0;

    if ( AND2( type >= 0, shapes > 0 ) ) {
      result = type;
    }
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: computeShapeSubsetBoundsMask - Create a mask array indicating which
         shapes intersect bounds.
INPUTS:  const char* baseFileName   Base name (no ext) of Shapefile to read.
         const Bounds bounds        Clip bounds.
         const int count            Number of shapes to match.
OUTPUTS: char mask[ count ]         1 = intersects bounds, else 0.
RETURNS: int 1 if at least one shape was within bounds, else 0.
******************************************************************************/

int computeShapeSubsetBoundsMask( const char* baseFileName,
                                  const Bounds bounds,
                                  const int count,
                                  char mask[] ) {

  PRE05( baseFileName, *baseFileName, isValidBounds( bounds), count > 0, mask);
  int result = 0;
  int ok = 0;
  SHPHandle handle = SHPOpen( baseFileName, "rb" );

  if ( ! handle ) {
    fprintf( stderr, "\nFailed to open Shapefile '%s'\n", baseFileName );
    perror( "because" );
  } else {
    int shapes = 0;
    int type = 0;
    double minimums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    double maximums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    Bounds dataBounds = { { 0.0, 0.0 }, { 0.0, 0.0 } };
    SHPGetInfo( handle, &shapes, &type, minimums, maximums );
    dataBounds[ LONGITUDE ][ MINIMUM ] = minimums[ 0 ];
    dataBounds[ LATITUDE  ][ MINIMUM ] = minimums[ 1 ];
    dataBounds[ LONGITUDE ][ MAXIMUM ] = maximums[ 0 ];
    dataBounds[ LATITUDE  ][ MAXIMUM ] = maximums[ 1 ];
    ok = shapes == count;

    DEBUG( fprintf( stderr,
                    "computeShapeSubsetBoundsMask(): "
                    "%s shapes = %d, type = %d, mask = %p\n",
                    baseFileName, shapes, type, mask ); )

    if ( ! ok ) {
      fprintf( stderr, "\nInvalid number of shapes in Shapefile '%s': "
               "actual %d, expected %d\n",
               baseFileName, count, shapes );
    } else if ( overlap( (const double (*)[2]) dataBounds, bounds ) ) {
      int index = 0;

      /* Check each shape: */

      for ( index = 0; AND2( ok, index < shapes ); ++index ) {
        SHPObject* shape = 0;
        DEBUG( fprintf( stderr, "Calling SHPReadObject()...\n" ); )
        shape = SHPReadObject( handle, index );
        ok = shape != 0;
        DEBUG( fprintf( stderr, "done. ok = %d\n", ok ); )

        if ( shape ) {
          int in = 0;
          dataBounds[ LONGITUDE ][ MINIMUM ] = shape->dfXMin;
          dataBounds[ LATITUDE  ][ MINIMUM ] = shape->dfYMin;
          dataBounds[ LONGITUDE ][ MAXIMUM ] = shape->dfXMax;
          dataBounds[ LATITUDE  ][ MAXIMUM ] = shape->dfYMax;
          in = overlap( (const double (*)[2]) dataBounds, bounds );
          mask[ index ] = in;
          result += in;

          DEBUG2( fprintf( stderr,
                          "overlap( [%lg %lg][%lg %lg], [%lg %lg][%lg %lg])"
                          " = %d\n",
                          dataBounds[0][0], dataBounds[0][1],
                          dataBounds[1][0], dataBounds[1][1],
                          bounds[0][0], bounds[0][1],
                          bounds[1][0], bounds[1][1], in ); )

          SHPDestroyObject( shape ), shape = 0;
        }
      }

      DEBUG( fprintf( stderr,
                      "%d shapes intersect bounds [%lg %lg][%lg %lg]\n",
                      result,
                      bounds[0][0], bounds[0][1],
                      bounds[1][0], bounds[1][1] ); )

      if ( ! ok ) {
        result = 0;
      }
    }

    SHPClose( handle ), handle = 0;
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: subsetByCOMID - Initialize a mask array indicating which
         rows are upstream of row with given comid.
INPUTS:  const ShapeData* const shapeData  Rows with comid, fromnode, tonode.
         const int comid                   COMID to match and search upstream.
         char mask[ shapeData->rows ]      Allocated mask array.
OUTPUTS: char mask[ shapeData->rows ]      Set to 1 if upstream of comid.
RETURNS: int number of upstream rows.
******************************************************************************/

int subsetByCOMID( const ShapeData* const shapeData,
                   const int comid,
                   char mask[] ) {

  PRE06( isValidShapeData( shapeData ),
         hasIntegerColumn( shapeData, "COMID" ),
         hasIntegerColumn( shapeData, "FROM_NODE" ),
         hasIntegerColumn( shapeData, "TO_NODE" ),
         comid > 0,
         mask );

  int result = 0;
  const int rows      = shapeData->rows;
  const int columns   = shapeData->columns;
  const Value* values = shapeData->values;
  const char** const columnNames = (const char** const) shapeData->columnNames;
  const int comidColumn    = indexOfString( "COMID",    columnNames, columns );
  const int fromNodeColumn = indexOfString( "FROM_NODE", columnNames, columns );
  const int toNodeColumn   = indexOfString( "TO_NODE",   columnNames, columns );
  DEBUG( const int gnisNameColumn =
         indexOfString( "GNIS_NAME", columnNames, columns ); )
  int fromNode = 0; /* The fromNode of row with given comid: */
  int row = 0;

  memset( mask, 0, rows * sizeof *mask ); /* Clear mask. */

  /* Search for row with given comid then get fromNode of that row: */

  for ( row = 0; AND2( fromNode == 0, row < rows ); ++row, values += columns ) {
    const int thisRowComid = values[ comidColumn ].i;

    if ( thisRowComid == comid ) { /* Found. */
      fromNode = values[ fromNodeColumn ].i;
      mask[ row ] = 1; /* Include selected flowline. */
      DEBUG( fprintf( stderr,
                       "subsetByComId( comid = %d ): "
                       "found row %5d %s %d -> %d\n",
                       comid,
                       row,
                       shapeData->values[ row * columns + gnisNameColumn ].s,
                       shapeData->values[ row * columns + fromNodeColumn ].i,
                       shapeData->values[ row * columns + toNodeColumn   ].i);)
    }
  }

  /*
   * If found row with comid then mark all (recursively) upstream rows that
   * flow into the selected flowline's fromNode:
   */

  if ( fromNode > 0 ) {
    result = 1 +
      flagUpstreamNodes( rows, columns, shapeData->values,
                         fromNodeColumn, toNodeColumn, fromNode,
                         mask );

    DEBUG( fprintf( stderr,
                    "after flagUpstreamNodes( fromNode = %d ), result = %d\n",
                    fromNode, result ); )

#ifdef DEBUGGING
for ( row = 0; row < rows; ++row )
  if ( mask[ row ] )
    fprintf( stderr, "masked row %5d %s %d -> %d\n",
            row,
            shapeData->values[ row * columns + gnisNameColumn ].s,
            shapeData->values[ row * columns + fromNodeColumn ].i,
            shapeData->values[ row * columns + toNodeColumn   ].i );
#endif

  }

  POST0( IN_RANGE( result, 0, shapeData->rows ) );
  return result;
}



/******************************************************************************
PURPOSE: writeSubsetCSVFile - Create a csv file containing lines from input csv
         files matching time and either single estcode or mask[] of estcodes.
INPUTS:  const char* inputDBFFileName   Name of dbf file to read.
         const char* inputCSVDirectory  Name of directory containing csv files
                                        to read.
         const char* outputCSVFileName  Name of csv file to create and write.
         const int yyyymmdd1        Starting date of subset.
         const int yyyymmdd2        Ending   date of subset.
         const char* estcode        0 or single ESTCODE to match.
         const char mask[]          If ESTCODE is 0 and mask is non-0 then
                                    use mask[] to lookup ESTCODEs from dbf.
RETURNS: int number of values written.
******************************************************************************/

int writeSubsetCSVFile( const char* inputDBFFileName,
                        const char* inputCSVDirectory,
                        const char* outputCSVFileName,
                        const int yyyymmdd1,
                        const int yyyymmdd2,
                        const char* estcode,
                        const char mask[] ) {

 PRE09( inputDBFFileName, *inputDBFFileName,
        inputCSVDirectory, *inputCSVDirectory,
        outputCSVFileName, *outputCSVFileName,
        isValidYearMonthDay( yyyymmdd1 ),
        isValidYearMonthDay( yyyymmdd2 ),
        yyyymmdd1 <= yyyymmdd2 );

  int result = 0;
  DBFHandle inputDBFFile = DBFOpen( inputDBFFileName, "rb" );

  if ( inputDBFFile ) {
    const int rows = DBFGetRecordCount( inputDBFFile );
    int ok = rows > 0;

    if ( ok ) {
      const int columns = DBFGetFieldCount( inputDBFFile );
      ok = columns > 0;

      if ( ok ) {
        const int estcodeColumn0 = DBFGetFieldIndex( inputDBFFile, "ESTCODE_N" );
        const int estcodeColumn =
          estcodeColumn0 != -1 ? estcodeColumn0
          : DBFGetFieldIndex( inputDBFFile, "ESTCODE" );
        ok = estcodeColumn >= 0;
        DEBUG2( fprintf( stderr, "writeSubsetCSVFile( %s ), column = %d\n",
                         estcode, estcodeColumn ); )

        if ( ok ) {
          FILE* outputCSVFile = fopen( outputCSVFileName , "w" );
          ok = outputCSVFile != 0;

          if ( ok ) {
            enum { TAG_SIZE = 8, MEMO_SIZE = TAG_SIZE * 4000 };
            /* Memoize codes to avoid processing csv files more than once. */
            char codes[ MEMO_SIZE ] = " "; /* Init with 1 space delimiter. */
            char tag[ TAG_SIZE ] = " ";    /* Init with 1 space delimiter. */
            size_t memoLength = 1;
            CSVHeader header = "";
            const int timestepSize =
              strstr( outputCSVFileName, "yearly" ) ? YEARLY
              : strstr( outputCSVFileName, "monthly" ) ? MONTHLY
              : DAILY;
            const int isMultiColumn =
              OR2( strstr( inputDBFFileName, "tide_" ),
                   strstr( inputDBFFileName, "estuary_flushing" ) );
            const int components =
              strstr( inputCSVDirectory, "current" ) ? 2 : 1;
            const int columns =
              OR2( estcode, isMultiColumn ) ? 0 : 2 + components;
            /* Output columns: ESTCODE,YYYY-MM-DD,value1[,value2]. */
            /* columns = 0 means output all csv columns, else just 1st 3 or 4*/
            int row = 0;
            memset( header, 0, sizeof header );

            for ( row = 0; row < rows; ++row ) {
              const char* const code =
                DBFReadStringAttribute( inputDBFFile, row, estcodeColumn );

              if ( code ) {
                const int selected =
                  estcode ? ! strcmp( code, estcode )
                  : mask ? mask[ row ]
                  : 1;

                if ( selected ) {
                  const size_t codeLength = strlen( code );
                  int found = 0;
                  int validCount = 1;

                  /* If code is not too long, construct a space-delimited tag*/

                  if ( codeLength + 2 < sizeof tag / sizeof *tag ) {
                    strcpy( tag + 1, code ); /* Copy after leading space. */
                    tag[ codeLength + 1 ] = ' '; /* Add ending space. */
                    tag[ codeLength + 2 ] = '\0'; /* And terminate. */
                    found = strstr( codes, tag ) != 0; /* Code processed? */

                    if ( ! found ) { /* If not, then try to store it: */
                      const size_t tagLength = codeLength + 2;

                      if ( memoLength + tagLength <
                           sizeof codes / sizeof *codes ) {
                        strcpy( codes + memoLength - 1, tag ); /* last space*/
                        memoLength += tagLength - 1; /* Omit written space. */
                      }
                    }
                  }

                  if ( ! found ) { /* Not already memoized, process it: */
                    validCount =
                      appendCSVFile( inputCSVDirectory, code,
                                     yyyymmdd1, yyyymmdd2, timestepSize,
                                     columns, header,
                                     outputCSVFile );
                    result += validCount;
                  }
                }
              }
            }

            fclose( outputCSVFile ), outputCSVFile = 0;
            DEBUG( fprintf( stderr, "%d csv values written to %s.\n",
                            result, outputCSVFileName ); )
          }

          if ( result ) {
            result = sortUniqFile( outputCSVFileName, 1 );
          }

          if ( result == 0 ) {
            unlink( outputCSVFileName );
          }
        }
      }
    }

    DBFClose( inputDBFFile ), inputDBFFile = 0;
  }

  DEBUG2( fprintf( stderr, "writeSubsetCSVFile() returning %d\n", result ); )
  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: writeSubsetCSVFileById - Create a csv file containing lines from
         input csv files matching unmasked id of columnName from dbf file.
INPUTS:  const char* inputDBFFileName   Name of dbf file to read.
         const char* inputCSVFileName   Name of csv files to read.
         const char* outputCSVFileName  Name of csv file to create and write.
         const char* const columnName   Name of dbf file column with id to match
         const int allowEmptyOutputCSV  Allow output csv file w/ no data lines?
         const char mask[]              If 1 then include this column in output
                                        csv.
RETURNS: int number of values written.
******************************************************************************/

int writeSubsetCSVFileById( const char* inputDBFFileName,
                            const char* inputCSVFileName,
                            const char* outputCSVFileName,
                            const char* const columnName,
                            const int allowEmptyOutputCSV,
                            const char mask[] ) {

 PRE010( inputDBFFileName, *inputDBFFileName,
         inputCSVFileName, *inputCSVFileName,
         outputCSVFileName, *outputCSVFileName,
         columnName, *columnName,
         IS_BOOL( allowEmptyOutputCSV ),
         mask );

  int result = 0;
  DBFHandle inputDBFFile = DBFOpen( inputDBFFileName, "rb" );
  DEBUG2( fprintf( stderr, "writeSubsetCSVFileById( %s )\n", columnName ); )

  if ( inputDBFFile ) {
    const int rows = DBFGetRecordCount( inputDBFFile );
    int ok = rows > 0;

    if ( ok ) {

      /* Array of subset ids to output: */

      long long* ids = NEW( long long, rows );
      ok = ids != 0;

      if ( ids ) {
        const int column = DBFGetFieldIndex( inputDBFFile, columnName );
        ok = column >= 0;
        DEBUG2( fprintf( stderr, "column = %d\n", column ); )

        if ( ok ) {
          size_t count = 0;
          int row = 0;

          /* Copy selected ids: */

          for ( row = 0; row < rows; ++row ) {

            if ( mask[ row ] ) {
              const int id =
                DBFReadIntegerAttribute( inputDBFFile, row, column );

              if ( id > 0 ) {
                ids[ count ] = id;
                ++count;
              }
            }
          }

          if ( count ) {
            FILE* inputCSVFile = fopen( inputCSVFileName , "r" );
            ok = inputCSVFile != 0;

            if ( ok ) {
              FILE* outputCSVFile = fopen( outputCSVFileName , "w" );
              ok = outputCSVFile != 0;

              if ( ok ) {

                /* Copy header line: */

                ok = copyFileLine( inputCSVFile, outputCSVFile );

                if ( ok ) {
                  shellsortI( ids, count );
                  result =
                    copyMatchedLines( inputCSVFile, count, ids, outputCSVFile );
                }

                fclose( outputCSVFile ), outputCSVFile = 0;

                if ( IS_ZERO2( result, allowEmptyOutputCSV ) ) {
                  unlink( outputCSVFileName );
                }
              }

              fclose( inputCSVFile ), inputCSVFile = 0;
            }
          }
        }

        FREE( ids );
      }
    }

    DBFClose( inputDBFFile ), inputDBFFile = 0;
  }

  DEBUG2(fprintf(stderr, "writeSubsetCSVFileById() returning %d\n",result);)
  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: readAndClipShapes - Read and clip and return array of
         shapes (with ids) clipped to the given bounds (and mask, if mask != 0).
INPUTS:  const char* baseFileName   Base name (no ext) of Shapefile to read.
         const Bounds bounds        Clip bounds.
         double minimumAdjacentVertexDistance  Adjacent vertices closer than
                                               this (in either x or y) will
                                               be merged.
         char* const mask            0 or mask[] to filter shapes before clip.
OUTPUTS: char* const mask            If non-0, updated mask[] after clip.
         int* count                 Length of returned array.
         int* isPolyline            Is shape a polyline?
RETURNS: PolygonShape* if successful, else 0 and a failure message is printed.
******************************************************************************/

PolygonShape* readAndClipShapes( const char* baseFileName, const Bounds bounds,
                                 double minimumAdjacentVertexDistance,
                                 char* const mask,
                                 int* count, int* isPolyline ) {

  PRE04( baseFileName, bounds, count, isPolyline );
  PolygonShape* result = 0;
  int ok = 0;
  SHPHandle handle = SHPOpen( baseFileName, "rb" );
  *count = 0;

  if ( ! handle ) {
    fprintf( stderr, "\nFailed to open Shapefile '%s'\n", baseFileName );
    perror( "because" );
  } else {
    int shapes = 0;
    int type = 0;
    double minimums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    double maximums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    Bounds dataBounds = { { 0.0, 0.0 }, { 0.0, 0.0 } };
    SHPGetInfo( handle, &shapes, &type, minimums, maximums );
    dataBounds[ LONGITUDE ][ MINIMUM ] = minimums[ 0 ];
    dataBounds[ LATITUDE  ][ MINIMUM ] = minimums[ 1 ];
    dataBounds[ LONGITUDE ][ MAXIMUM ] = maximums[ 0 ];
    dataBounds[ LATITUDE  ][ MAXIMUM ] = maximums[ 1 ];
    ok = 1;
    *isPolyline = IN3( type, SHPT_ARC, SHPT_ARCZ );

    DEBUG( fprintf( stderr,
                    "readAndClipShapes(): "
                    "%s shapes = %d, type = %d\n",
                    baseFileName, shapes, type ); )

    if ( AND3( shapes > 0,
               IN5( type, SHPT_POLYGON, SHPT_POLYGONZ, SHPT_ARC, SHPT_ARCZ ),
               overlap( (const double (*)[2]) dataBounds, bounds ) ) ) {
      result = NEW( PolygonShape, shapes );
      ok = result != 0;

      if ( result ) {
        int index = 0;

        /* Make a static GPC clipping polygon, 'clip', from bounds: */

        gpc_vertex clipVertices[ 5 ] =
          { {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0} };
        int holes[ 5 ] = { 0, 0, 0, 0, 0 };
        gpc_vertex_list clipContour = { 5, 0 };
        gpc_polygon clip = { 1, 0, 0 };

        clipVertices[ 0 ].x = bounds[ LONGITUDE ][ MINIMUM ];
        clipVertices[ 0 ].y = bounds[ LATITUDE  ][ MINIMUM ];
        clipVertices[ 1 ].x = bounds[ LONGITUDE ][ MINIMUM ];
        clipVertices[ 1 ].y = bounds[ LATITUDE  ][ MAXIMUM ];
        clipVertices[ 2 ].x = bounds[ LONGITUDE ][ MAXIMUM ];
        clipVertices[ 2 ].y = bounds[ LATITUDE  ][ MAXIMUM ];
        clipVertices[ 3 ].x = bounds[ LONGITUDE ][ MAXIMUM ];
        clipVertices[ 3 ].y = bounds[ LATITUDE  ][ MINIMUM ];
        clipVertices[ 4 ].x = clipVertices[ 0 ].x;
        clipVertices[ 4 ].y = clipVertices[ 0 ].y;

        clipContour.vertex = clipVertices;

        clip.hole = holes;
        clip.contour = &clipContour;
        DEBUG( fprintf( stderr, "clip:\n" ); printPolygon( &clip ); )

        /* Clip each shape: */

        for ( index = 0; AND2( ok, index < shapes ); ++index ) {
          DEBUG( fprintf( stderr, "index %d: mask = %p [%d]\n",
                          index, mask, mask ? mask[ index ] : 1 ); )

          if ( OR2( mask == 0, mask[ index ] ) ) {
            const int oldCount = *count;
            SHPObject* shape = 0;
            DEBUG( fprintf( stderr, "Calling SHPReadObject()...\n" ); )
            shape = SHPReadObject( handle, index );
            ok = shape != 0;
            DEBUG( fprintf( stderr, "done. ok = %d\n", ok ); )

            if ( shape ) {
              int shapeOverlaps = 0;
              dataBounds[ LONGITUDE ][ MINIMUM ] = shape->dfXMin;
              dataBounds[ LATITUDE  ][ MINIMUM ] = shape->dfYMin;
              dataBounds[ LONGITUDE ][ MAXIMUM ] = shape->dfXMax;
              dataBounds[ LATITUDE  ][ MAXIMUM ] = shape->dfYMax;
              shapeOverlaps =
                overlap( (const double (*)[2]) dataBounds, bounds );

              if ( shapeOverlaps ) {
                DEBUG2(fprintf(stderr,
                            "overlap( [%lg %lg][%lg %lg], [%lg %lg][%lg %lg])"
                              " = 1\n",
                              dataBounds[0][0], dataBounds[0][1],
                              dataBounds[1][0], dataBounds[1][1],
                              bounds[0][0], bounds[0][1],
                              bounds[1][0], bounds[1][1] ); )
                gpc_polygon copy = { 0, 0, 0 };
                DEBUG2(fprintf(stderr,"calling makePolygon() #%d ...\n",index);)
                ok = makePolygon( shape,  minimumAdjacentVertexDistance, &copy,
                                  result[ *count ].bounds );
                DEBUG2( fprintf( stderr, "done. makePolygon = %d, num_contours = %d\n", ok, copy.num_contours ); )
                DEBUG2( fprintf( stderr, "copy:\n" ); printPolygon( &copy ); )

                if ( AND2( ok, copy.num_contours > 0 ) ) {
                  int contours = 0;

                  if ( *isPolyline ) {
                    DEBUG2( fprintf( stderr, "calling clipPolylines()...\n"); )
                    clipPolylines( &copy, bounds, &result[ *count ].polygon,
                                   result[ *count ].bounds );
                  } else {
                    DEBUG2( fprintf(stderr,"calling gpc_polygon_clip()...\n");)
                    gpc_polygon_clip( GPC_INT, &copy, &clip,
                                      &result[ *count ].polygon );
                  }

                  DEBUG2( fprintf( stderr, "done.\n"); )
                  gpc_free_polygon( &copy );
                  DEBUG2( fprintf( stderr, "result:\n" ); )
                  DEBUG2( fprintf( stderr, "bounds = [%lg, %lg] [%lg, %lg]\n",
                                  result[ *count ].bounds[ 0 ][ 0 ],
                                  result[ *count ].bounds[ 0 ][ 1 ],
                                  result[ *count ].bounds[ 1 ][ 0 ],
                                  result[ *count ].bounds[ 1 ][ 1 ] ); )
                  DEBUG2( printPolygon( &result[ *count ].polygon ); )
                  contours = result[ *count ].polygon.num_contours;

                  if ( AND3( contours > 0,
                             minimumInt( contours,
                                         result[ *count ].polygon.hole ) == 0,
                             OR2( *isPolyline,
                                  ensureCorrectVertexOrder(
                                    &result[ *count ].polygon ) ) ) ) {
                    DEBUG2(fprintf(stderr, "keep shape id %d (%d contours).\n",
                                     shape->nShapeId, contours ); )
                    result[ *count ].id = shape->nShapeId;
                    *count += 1;
                  } else {
                    gpc_free_polygon( &result[ *count ].polygon );
                  }
                }
              }

              SHPDestroyObject( shape ), shape = 0;
            }

            if ( AND2( mask, *count == oldCount ) ) { /* If not in subset. */
              mask[ index ] = 0; /* Mask-out the row. */
            }
          }
        }
      }
    }

    SHPClose( handle ), handle = 0;
  }

  if ( AND2( OR2( ! ok, *count == 0 ), result != 0 ) ) {
    deallocatePolygons( *count, result ), *count = 0, result = 0;
  }

  POST0( IMPLIES_ELSE( result == 0,
                        *count == 0,
                        AND6( *count > 0,
                              IS_BOOL( *isPolyline ),
                              result[ 0 ].id >= 0,
                              result[ 0 ].polygon.num_contours > 0,
                              result[ 0 ].polygon.hole != 0,
                              result[ 0 ].polygon.contour != 0 ) ) );
  return result;
}




/******************************************************************************
PURPOSE: readAndTriangulateShapes - Read and triangulate and return array of
         shapes (with ids).
INPUTS:  const char* baseFileName   Base name (no ext) of Shapefile to read.
         const double minimumAdjacentVertexDistance  0 or distance in degrees
                                                     to skip adjacent vertices.
OUTPUTS: int* count                 Length of returned array.
RETURNS: PolygonShape* if successful, else 0 and a failure message is printed.
******************************************************************************/

PolygonShape* readAndTriangulateShapes( const char* baseFileName,
                                    const double minimumAdjacentVertexDistance,
                                        int* count) {

  PRE03( baseFileName, minimumAdjacentVertexDistance >= 0.0, count );
  PolygonShape* result = 0;
  int ok = 0;
  SHPHandle handle = SHPOpen( baseFileName, "rb" );
  *count = 0;

  if ( ! handle ) {
    fprintf( stderr, "\nFailed to open Shapefile '%s'\n", baseFileName );
    perror( "because" );
  } else {
    int shapes = 0;
    int type = 0;
    SHPGetInfo( handle, &shapes, &type, 0, 0 );
    ok = AND2( shapes > 0, IN3( type, SHPT_POLYGON, SHPT_POLYGONZ ) );

    if ( ok ) {
      result = NEW( PolygonShape, shapes );
      ok = result != 0;

      if ( result ) {
        int index = 0;

        /* Read each shape: */

        for ( index = 0; AND2( ok, index < shapes ); ++index ) {
          SHPObject* shape = 0;
          DEBUG( fprintf( stderr, "Calling SHPReadObject()...\n" ); )
          shape = SHPReadObject( handle, index );
          DEBUG( fprintf( stderr, "done.\n" ); )
          ok = shape != 0;

          if ( shape ) {
            gpc_polygon copy = { 0, 0, 0 };
            DEBUG2( fprintf( stderr, "calling makePolygon()...\n"); )
            ok = makePolygon( shape,  minimumAdjacentVertexDistance,
                              &copy, result[ *count ].bounds );
            DEBUG2( fprintf( stderr, "done. makePolygon = %d\n", ok ); )
            DEBUG( fprintf( stderr, "copy:\n" ); printPolygon( &copy ); )

            if ( AND2( ok, copy.num_contours > 0 ) ) {
#if 0
              DEBUG2( fprintf( stderr, "calling ensureCorrectVertexOrder()...\n"); )
              ok = ensureCorrectVertexOrder( &copy );
              DEBUG2( fprintf( stderr, "done. ensureCorrectVertexOrder = %d\n", ok ); )
              DEBUG( fprintf( stderr, "copy:\n" ); printPolygon( &copy ); )
              ok = 1; /* Allow polygons with incorrect vertex order. */
#endif
              if ( ok ) {
                int strips = 0;
                DEBUG2(fprintf(stderr,"calling gpc_polygon_to_tristrip()...\n");)
                gpc_polygon_to_tristrip( &copy, &result[ *count ].triangles );
                DEBUG2( fprintf( stderr, "done.\n"); )
                gpc_free_polygon( &copy );
                memset( &copy, 0, sizeof copy );
                DEBUG( fprintf( stderr, "result:\n" ); )
                DEBUG( fprintf( stderr, "bounds = [%lg, %lg] [%lg, %lg]\n",
                                result[ *count ].bounds[ 0 ][ 0 ],
                                result[ *count ].bounds[ 0 ][ 1 ],
                                result[ *count ].bounds[ 1 ][ 0 ],
                                result[ *count ].bounds[ 1 ][ 1 ] ); )
                DEBUG( printTriangles( &result[ *count ].triangles ); )
                strips = result[ *count ].triangles.num_strips;

                if ( strips > 0 ) {
                  DEBUG2( fprintf( stderr, "keep.\n" ); )
                  result[ *count ].id = shape->nShapeId;
                  *count += 1;
                } else {
                  gpc_free_tristrip( &result[ *count ].triangles );
                }
              }
            }

            SHPDestroyObject( shape ), shape = 0;
          }
        }
      }
    }

    SHPClose( handle ), handle = 0;
  }

  if ( AND2( OR2( ! ok, *count == 0 ), result != 0 ) ) {
    deallocatePolygons( *count, result ), *count = 0, result = 0;
  }

  POST0( IMPLIES_ELSE( result == 0,
                        *count == 0,
                        AND4( *count > 0,
                              result[ 0 ].id >= 0,
                              result[ 0 ].triangles.num_strips > 0,
                              result[ 0 ].triangles.strip != 0 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: pointInTriangles - Is the specified point (x, y) in any of the set of
         triangles? If so return its index, else -1.
INPUTS:  double x          X-coordinate of point to test.
         double y          Y-coordinate of point to test.
         int count         Number of polygons in polygons[].
         const PolygonShape polygons[ count ]  Array of triangulated polygons.
RETURNS: int index [0, count - 1] if the point is inside the indexed triangles,
         else -1.
******************************************************************************/

int pointInTriangles( double x, double y,
                      int count, const PolygonShape polygons[] ) {

  PRE07( ! isNan( x ) , ! isNan( y ), count > 0, polygons,
         isValidBounds( (const double (*)[2]) polygons[ 0 ].bounds ),
         polygons[ 0 ].triangles.num_strips > 0,
         polygons[ count - 1 ].triangles.num_strips > 0 );

  int result = -1;
  int index = 0;

  for ( index = 0; index < count; ++index ) {
    const PolygonShape* const polygon = polygons + index;
    const double xMinimum = polygon->bounds[ LONGITUDE ][ MINIMUM ];
    const double xMaximum = polygon->bounds[ LONGITUDE ][ MAXIMUM ];
    const double yMinimum = polygon->bounds[ LATITUDE  ][ MINIMUM ];
    const double yMaximum = polygon->bounds[ LATITUDE  ][ MAXIMUM ];
    const int outsideBounds =
      OR4( x < xMinimum, x > xMaximum, y < yMinimum, y > yMaximum );

    if ( ! outsideBounds ) {
      const gpc_tristrip* const tristrip = &polygon->triangles;
      const int strips = tristrip->num_strips;
      int strip = 0;

      for ( strip = 0; strip < strips; ++strip ) {
        const gpc_vertex_list* const vertex_list = tristrip->strip + strip;
        const int vertexCount = vertex_list->num_vertices;
        const gpc_vertex* const vertices = vertex_list->vertex;
        int vertexIndex = 0;
        double x1 = vertices[ 0 ].x;
        double y1 = vertices[ 0 ].y;
        double x2 = vertices[ 1 ].x;
        double y2 = vertices[ 1 ].y;
        CHECK( vertexCount >= 3 );

        for ( vertexIndex = 2; vertexIndex < vertexCount; ++vertexIndex ) {
          const double x3 = vertices[ vertexIndex ].x;
          const double y3 = vertices[ vertexIndex ].y;
          const int insideTriangle =
            pointInsideTriangle( x, y, x1, y1, x2, y2, x3, y3 );

          if ( insideTriangle ) {
            result = index;
            vertexIndex = vertexCount; /* Stop looping. */
            strip = strips;
            index = count;
          } else {
            x1 = x2;
            y1 = y2;
            x2 = x3;
            y2 = y3;
          }
        }
      }
    }
  }

  POST0( IN_RANGE( result, -1, count - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: nearestPolyline - Is the specified point (x, y) on any of the set of
         polylines? If so return the index of the closest one, else -1.
INPUTS:  const double x          X-coordinate of point to test.
         const double y          Y-coordinate of point to test.
         const int count         Number of polygons in polygons[].
         const PolygonShape polylines[ count ]  Array of polylines.
RETURNS: int index [0, count - 1] if the point is on the indexed polylines,
         else -1.
******************************************************************************/

int nearestPolyline( const double x, const double y,
                     const int count, const PolygonShape polylines[] ) {

  PRE05( ! isNan( x ) , ! isNan( y ), count > 0, polylines,
         isValidBounds( (const double (*)[2]) polylines[ 0 ].bounds ) );

  /*
   * For each polyline:
   *   If point (x, y) is within the expanded bounds of the polyline then
   *     for each line segment of polyline:
   *       If point is within the expanded bounds of the line segment then
   *         compute point-to-line distance and compare/store closest one.
   * If closest point-to-line distance is within tolerance return index else -1.
   * Note bounds to check are expanded slightly to handle degenerate
   * width/height dimensions such as due to vertical or horizontal lines.
   */

  const double boundsMargin = 1e-3; /* Expand line segment bounds slightly. */
  const double tolerance    = 1e-3; /* Close enough distance to line to accept*/
  double nearestDistance = DBL_MAX;
  int result = -1;
  int index = 0;

  for ( index = 0; index < count; ++index ) {
    const PolygonShape* const polygonShape = polylines + index;
    double xMinimum =
      polygonShape->bounds[ LONGITUDE ][ MINIMUM ] - boundsMargin;

    if ( x >= xMinimum ) {
      double xMaximum =
        polygonShape->bounds[ LONGITUDE ][ MAXIMUM ] + boundsMargin;

      if ( x <= xMaximum ) {
        double yMinimum =
          polygonShape->bounds[ LATITUDE ][ MINIMUM ] - boundsMargin;

        if ( y >= yMinimum ) {
          double yMaximum =
            polygonShape->bounds[ LATITUDE ][ MAXIMUM ] + boundsMargin;

          if ( y <= yMaximum ) {
            const gpc_polygon* const polyline = &polygonShape->polygon;
            const int contours = polyline->num_contours;
            int contour = 0;

            DEBUG2( fprintf( stderr,
                             "Point (%f, %f) inside polyline %4d "
                             "expanded bounds[%f %f][%f %f]\n",
                             x, y, index,
                             xMinimum, xMaximum, yMinimum, yMaximum ); )

            for ( contour = 0; contour < contours; ++contour ) {
              const gpc_vertex_list* const vertexList =
                polyline->contour + contour;
              const gpc_vertex* const vertices = vertexList->vertex;
              const int vertexCount = vertexList->num_vertices;
              int vertexIndex = 1;
              double x1 = vertices[ 0 ].x;
              double y1 = vertices[ 0 ].y;
              CHECK( vertexCount >= 2 );

              for ( ; vertexIndex < vertexCount; ++vertexIndex ) {
                double x2 = vertices[ vertexIndex ].x;
                double y2 = vertices[ vertexIndex ].y;

                if ( x1 < x2 ) {
                  xMinimum = x1 - boundsMargin;
                  xMaximum = x2 + boundsMargin;
                } else {
                  xMinimum = x2 - boundsMargin;
                  xMaximum = x1 + boundsMargin;
                }

                if ( IN_RANGE( x, xMinimum, xMaximum ) ) {

                  if ( y1 < y2 ) {
                    yMinimum = y1 - boundsMargin;
                    yMaximum = y2 + boundsMargin;
                  } else {
                    yMinimum = y2 - boundsMargin;
                    yMaximum = y1 + boundsMargin;
                  }

                  if ( IN_RANGE( y, yMinimum, yMaximum ) ) {
                    const double distance =
                      pointLineDistance( x, y, x1, y1, x2, y2 );

                    DEBUG2( fprintf( stderr,
                                     "pointLineDistance (%f, %f) to "
                                     "(%f, %f)--(%f %f ) = %e\n",
                                     x, y, x1, y1, x2, y2, distance ); )

                    if ( distance < nearestDistance ) {
                      nearestDistance = distance;
                      result = index;
                    }
                  }
                }

                x1 = x2;
                y1 = y2;
              }
            }
          }
        }
      }
    }
  }

  DEBUG2( fprintf( stderr, "result = %d, nearestDistance = %e\n",
                   result, nearestDistance ); )

  if ( nearestDistance > tolerance ) {
    result = -1;
  }

  POST0( IN_RANGE( result, -1, count - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: nearestPoint - Is the specified point (x, y) on any of the set of
         points? If so return its index, else -1.
INPUTS:  const double x          X-coordinate of point to test.
         const double y          Y-coordinate of point to test.
         const ShapeData* const shapeData  ShapeData with LONGITUDE, LATITUDE.
RETURNS: int index [0, count - 1] if the point is on the indexed polylines,
         else -1.
******************************************************************************/

int nearestPoint( const double x, const double y,
                  const ShapeData* const shapeData ) {

  PRE03( ! isNan( x ) , ! isNan( y ), isValidShapeData( shapeData ) );

  int result = -1;
  const double tolerance = 1e-3; /* Close enough distance. */
  double nearestDistance = DBL_MAX;
  const int rows = shapeData->rows;
  const int columns = shapeData->columns;
  const int longitudeColumn =
    indexOfString( "LONGITUDE", (const char**) shapeData->columnNames,
                   columns );
  const int latitudeColumn =
    indexOfString( "LATITUDE", (const char**) shapeData->columnNames,
                   columns );
  int row = 0;

  if ( AND2( IN_RANGE( longitudeColumn, 0, columns - 1 ),
             IN_RANGE( latitudeColumn, 0, columns - 1 ) ) ) {
    int longitudeIndex = longitudeColumn;
    int latitudeIndex  = latitudeColumn;

    for ( row = 0; row < rows; ++row,
          longitudeIndex += columns, latitudeIndex += columns ) {

      CHECK2( IN_RANGE( longitudeIndex, 0, rows * columns - 1 ),
              IN_RANGE( latitudeIndex,  0, rows * columns - 1 ) );
              
      {
        const double longitude = shapeData->values[ longitudeIndex ].d;
        const double latitude  = shapeData->values[ latitudeIndex  ].d;
        const double longitudeDistance =
          x < longitude ? longitude - x : x - longitude;
        const double latitudeDistance =
          y < latitude ? latitude - y : y - latitude;
        const double distance = longitudeDistance + latitudeDistance;

        if ( distance < nearestDistance ) {
          nearestDistance = distance;
          result = row;
        }
      }
    }

    if ( nearestDistance > tolerance ) {
      result = -1;
    }
  }

  POST0( IN_RANGE( result, -1, shapeData->rows - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: makePolygon - Make a GPC-polygon from an ESRI Shape polygon/polyline.
INPUTS:  const SHPObject* shape  Shape to copy.
         double minimumAdjacentVertexDistance  Adjacent vertices closer than
                                               this (in either x or y) will
                                               be merged.
OUTPUTS: gpc_polygon* polygon  GPC polygon.
         Bounds bounds         Bounds of sparsed polygon vertices.
RETURNS: int 1 if successful (that is no allocation failures, but possibly no
         contours in polygon), else 0 and a failure message is printed.
NOTES:   Call gpc_free_polygon( polygon ) when finished with it.
         http://shapelib.maptools.org/shp_api.html
         http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
******************************************************************************/

int makePolygon( const SHPObject* shape, double minimumAdjacentVertexDistance,
                 gpc_polygon* polygon, Bounds bounds ) {

  PRE013( shape,
          IN5( shape->nSHPType,
               SHPT_POLYGON, SHPT_POLYGONZ, SHPT_ARC, SHPT_ARCZ ),
          shape->nShapeId >= 0,
          shape->nParts > 0,
          shape->panPartStart,
          shape->panPartType,
          shape->panPartType[ 0 ] == SHPP_RING,
          shape->nVertices > 0,
          shape->padfX,
          shape->padfY,
          minimumAdjacentVertexDistance >= 0.0,
          polygon,
          bounds );

  int result = 1;
  const int isPolygon = IN3( shape->nSHPType, SHPT_POLYGON, SHPT_POLYGONZ );
  const int minimumSparsedVertices = isPolygon ? 3 : 2;
  const int parts = shape->nParts;
  int sparsedParts = 0;
  DEBUG2( printShape( shape ); )
  memset( polygon, 0, sizeof *polygon );
  memset( bounds, 0, sizeof (Bounds) );

  /* Compute sparsed number of parts (with at least minimumSparsedVertices): */

  sparsedParts =
    computeSparsedPartCount( shape, minimumAdjacentVertexDistance,
                             minimumSparsedVertices );

  if ( sparsedParts > 0 ) {
    result = 0;
    polygon->hole = NEW( int, sparsedParts );

    if ( polygon->hole ) {
      polygon->contour = NEW( gpc_vertex_list, sparsedParts );

      if ( polygon->contour ) {
        int ok = 1;
        int part = 0;
        int sparsedPart = 0;
        int initializedBounds = 0;

        polygon->num_contours = sparsedParts;

        /* For each part: */

        for ( part = 0; AND2( ok, part < parts ); ++part ) {
          const int partVertexCount =
            parts == 1 ? shape->nVertices
            : part < parts - 1 ?
              shape->panPartStart[ part + 1 ] - shape->panPartStart[ part ]
            : shape->nVertices - shape->panPartStart[ part ];
          const int offset = shape->panPartStart[ part ];
          const double* const x = shape->padfX + offset;
          const double* const y = shape->padfY + offset;
          const int sparsedVertices =
            computeSparsedVertexCount( partVertexCount, x, y,
                                       minimumAdjacentVertexDistance,
                                       isPolygon );

          /* Allocate and copy sparse vertices: */

          if ( sparsedVertices >= ( minimumSparsedVertices + isPolygon ) ) {
            gpc_vertex* vertices = NEW( gpc_vertex, sparsedVertices );
            ok = vertices != 0;

            if ( vertices ) {
              CHECK( IN_RANGE( sparsedPart, 0, sparsedParts - 1 ) );
              polygon->contour[ sparsedPart ].num_vertices = sparsedVertices;
              polygon->contour[ sparsedPart ].vertex = vertices;
              copySparseVertices( partVertexCount, sparsedVertices, x, y,
                                  minimumAdjacentVertexDistance, isPolygon,
                                  &initializedBounds, bounds, vertices );
              CHECK( IN_RANGE( sparsedPart, 0, sparsedParts - 1 ) );

              /* Compute hole flag of sparse contour: */

              if ( isPolygon ) {
                const int counterClockwise =
                  signedAreaOfPolygon( &polygon->contour[sparsedPart] ) >= 0.0;
                polygon->hole[ sparsedPart ] = counterClockwise;
              }

              ++sparsedPart;
            }
          }
        }

        result = ok;
      }
    }
  }

  if ( ! result ) {
    gpc_free_polygon( polygon );
    memset( polygon, 0, sizeof *polygon );
    memset( bounds, 0, sizeof (Bounds) );
  }

  DEBUG2( printPolygon( polygon ); )

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( AND2( result, polygon->num_contours > 0 ),
                        AND4( IN_RANGE( polygon->num_contours, 1,
                                        shape->nParts ),
                              polygon->hole,
                              polygon->contour,
                              isValidBounds( (const double (*)[2]) bounds ) ),
                        IS_ZERO7( polygon->num_contours,
                                  polygon->hole,
                                  polygon->contour,
                                  bounds[ 0 ][ 0 ],
                                  bounds[ 0 ][ 1 ],
                                  bounds[ 1 ][ 0 ],
                                  bounds[ 1 ][ 1 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: clipPolylines - Clip polylines to a given bounds.
INPUTS:  const gpc_polygon* const polylines  GPC polylines to clip.
         const Bounds clipBounds             Clip bounds.
OUTPUTS: gpc_polygon* clippedPolylines       Clipped GPC polylines.
         Bounds clippedPolylinesBounds       Bounds of clippedPolylines.
RETURNS: int 1 if successful (that is no allocation failures, but possibly no
         contours in bounds), else 0 and a failure message is printed.
NOTES:   Call gpc_free_polygon( clippedPolylines ) when finished with it.
         http://www.cs.man.ac.uk/~toby/alan/software/gpc.html
******************************************************************************/

int clipPolylines( const gpc_polygon* const polylines,
                   const Bounds clipBounds,
                   gpc_polygon* clippedPolylines,
                   Bounds clippedPolylinesBounds ) {

  PRE04( polylines,
         isValidBounds( (const double (*)[2]) clipBounds ),
         polylines,
         isValidBounds( (const double (*)[2]) clippedPolylinesBounds ) );
  int result = 0;
  const int polylineCount = polylines->num_contours;
  memset( clippedPolylines, 0, sizeof *clippedPolylines );
  memset( clippedPolylinesBounds, 0, sizeof (Bounds) );

  if ( polylineCount > 0 ) { /* Count number of vertices for all polylines: */
    const int vertexCount = polygonVertexCount( polylines );

    if ( vertexCount > 0 ) { /* Copy vertices: */
      int* inputCounts = NEW( int, polylineCount );
      double* inputVertices = inputCounts ? NEW( double, vertexCount * 2 ) : 0;

      if ( inputVertices ) {
        int outputPolylineCount = 0;
        int outputVertexCount = 0;
        copyPolylineVertices( polylines, inputCounts, inputVertices );

        /* 1st call: get number of clipped polylines and total vertices: */

        subsetMapDouble( polylineCount, vertexCount,
                         inputCounts, inputVertices, 0.0, clipBounds,
                         &outputPolylineCount, &outputVertexCount,
                         0, 0 );

        /* Allocate storage for clipped results: */

        if ( AND2( outputPolylineCount > 0, outputVertexCount >= 2 ) ) {
          int* outputCounts = NEW( int, outputPolylineCount );
          double* outputVertices =
            outputCounts ? NEW( double, outputVertexCount * 2 ) : 0;

          if ( outputVertices ) {
            int outputPolylineCount2 = 0;
            int outputVertexCount2 = 0;

            /* 2nd call: get clipped counts and vertices: */

            subsetMapDouble( polylineCount, vertexCount,
                             inputCounts, inputVertices, 0.0, clipBounds,
                             &outputPolylineCount2, &outputVertexCount2,
                             outputCounts, outputVertices );

            CHECK2( outputPolylineCount2 == outputPolylineCount,
                    outputVertexCount2 == outputVertexCount );

            result = createPolyline( outputPolylineCount, outputCounts,
                                     outputVertices, clippedPolylines,
                                     clippedPolylinesBounds );
          }

          FREE( outputCounts );
          FREE( outputVertices );
        }

        FREE( inputCounts );
        FREE( inputVertices );
      }
    }
  }

  if ( ! result ) {
    gpc_free_polygon( clippedPolylines );
    memset( clippedPolylines, 0, sizeof *clippedPolylines );
    memset( clippedPolylinesBounds, 0, sizeof (Bounds) );
  }

  DEBUG2( printPolygon( clippedPolylines ); )

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( AND2( result, clippedPolylines->num_contours > 0 ),
                        AND4( minimumInt( clippedPolylines->num_contours,
                                          clippedPolylines->hole ) == 0,
                              maximumInt( clippedPolylines->num_contours,
                                          clippedPolylines->hole ) == 0,
                              clippedPolylines->contour,
                              isValidBounds( (const double (*)[2])
                                             clippedPolylinesBounds ) ),
                        IS_ZERO7( clippedPolylines->num_contours,
                                  clippedPolylines->hole,
                                  clippedPolylines->contour,
                                  clippedPolylinesBounds[ 0 ][ 0 ],
                                  clippedPolylinesBounds[ 0 ][ 1 ],
                                  clippedPolylinesBounds[ 1 ][ 0 ],
                                  clippedPolylinesBounds[ 1 ][ 1 ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: maximumPolygonContours - Maximum number of contours in set of polygons
INPUTS:  const int count                  Number of polygons.
         const PolygonShape polygons[ count ]  Polygons to check.
RETURNS: int maximum number of contours in set of polygons.
******************************************************************************/

int maximumPolygonContours( int count, const PolygonShape polygons[] ) {
  PRE02( count > 0, polygons );
  int result = 0;
  int polygon = 0;

  for ( polygon = 0; polygon < count; ++polygon ) {
    const int contours = polygons[ polygon ].polygon.num_contours;

    if ( contours > result ) {
      result = contours;
    }
  }

  POST0( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: maximumPolygonVertices - Maximum number of vertices in set of polygons
INPUTS:  const int count                  Number of polygons.
         const PolygonShape polygons[ count ]  Polygons to check.
RETURNS: int maximum number of vertices in set of polygons.
******************************************************************************/

int maximumPolygonVertices( int count, const PolygonShape polygons[] ) {
  PRE02( count > 0, polygons );
  int result = 0;
  int polygon = 0;

  for ( polygon = 0; polygon < count; ++polygon ) {
    result += polygonVertexCount( &polygons[ polygon ].polygon );
  }

  POST0( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: polygonVertexCount - Number of vertices in polygon.
INPUTS:  const gpc_polygon* polygon  Polygon to check.
RETURNS: int number of vertices in polygon.
******************************************************************************/

int polygonVertexCount( const gpc_polygon* polygon ) {
  PRE0( polygon );
  const int contours = polygon->num_contours;
  int contour = 0;
  int result = 0;

  for ( contour = 0; contour < contours; ++contour ) {
    const gpc_vertex_list vertexList = polygon->contour[ contour ];
    const int vertices = vertexList.num_vertices;
    result += vertices;
  }

  POST0( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: copyPolygonVertices - Copy vertices from polygon to x and y arrays.
INPUTS:  const gpc_polygon* polygon  Polygon to copy from.
         int closeRing               Copy first vertex to last to close ring?
OUTPUTS: int starts[]                0-based indices of vertex lists.
         double x[]                  X-coordinates of vertices.
         double y[]                  Y-coordinates of vertices.
******************************************************************************/

void copyPolygonVertices( const gpc_polygon* polygon, int closeRing,
                          int starts[], double x[], double y[] ) {

  PRE07( polygon, IS_BOOL( closeRing ), polygon->num_contours > 0,
         IMPLIES( polygon->num_contours > 1, starts ),
         x, y, x != y );

  double* xp = x;
  double* yp = y;
  const int contours = polygon->num_contours;
  int contour = 0;
  int offset = 0;

  if ( starts ) {
    starts[ 0 ] = 0;
  }

  for ( contour = 0; contour < contours; ++contour ) {
    const gpc_vertex_list vertexList = polygon->contour[ contour ];
    const int vertices = vertexList.num_vertices;
    const gpc_vertex* v = vertexList.vertex;
    const double x0 = v->x;
    const double y0 = v->y;
    int vertex = 0;

    if ( AND2( starts, contour < contours - 1 ) ) {
      offset += vertices + closeRing;
      starts[ contour + 1 ] = offset;
    }

    for ( vertex = 0; vertex < vertices; ++vertex ) {
      *xp++ = v->x;
      *yp++ = v->y;
      ++v;
    }

    if ( closeRing ) { /* Copy 1st vertex to last to ensure closed ring: */
      *xp++ = x0;
      *yp++ = y0;
    }
  }

  POST04( x[ 0 ] == polygon->contour->vertex->x,
          y[ 0 ] == polygon->contour->vertex->y,
          IMPLIES( starts, starts[ 0 ] == 0 ),
          IMPLIES( AND2( starts, polygon->num_contours > 1 ),
                   AND2( starts[ 1 ] > starts[ 0 ],
                         starts[ polygon->num_contours - 1 ] >
                           starts[ polygon->num_contours - 2 ] ) ) );
}



/******************************************************************************
PURPOSE: copyPolylineVertices - Copy vertices from polyline to xy array.
INPUTS:  const gpc_polygon* polygon  Polyline to copy from.
OUTPUTS: int starts[]                Number of vertices per polyline.
         double xy[]                 X/Y-coordinates of vertices.
******************************************************************************/

void copyPolylineVertices( const gpc_polygon* polygon,
                           int counts[], double xy[] ) {

  PRE04( polygon, polygon->num_contours > 0, counts, xy );

  int* count = counts;
  double* xyp = xy;
  const int contours = polygon->num_contours;
  int contour = 0;

  for ( contour = 0; contour < contours; ++contour ) {
    const gpc_vertex_list* const vertexList = polygon->contour + contour;
    const int vertices = vertexList->num_vertices;
    const gpc_vertex* v = vertexList->vertex;
    int vertex = 0;
    *count++ = vertices;

    for ( vertex = 0; vertex < vertices; ++vertex ) {
      *xyp++ = v->x;
      *xyp++ = v->y;
      ++v;
    }
  }

  POST03( counts[ 0 ] > 0,
          xy[ 0 ] == polygon->contour->vertex->x,
          xy[ 1 ] == polygon->contour->vertex->y );
}



/******************************************************************************
PURPOSE: createPolyline - Allocate and copy xy vertices to polyline.
INPUTS:  const int polylineCount                  Number of polylines.
         const int vertexCounts[ polylineCount ]  # of vertices per polyline.
         const double xy[]                     X/Y-coordinates of all vertices.
OUTPUTS: gpc_polygon* polygon  Polygon with allocated/copied hole and vertices.
         Bounds bounds         Bounds of vertices.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stdout
******************************************************************************/

int createPolyline( const int polylineCount, const int vertexCounts[],
                    const double xy[], gpc_polygon* polygon,
                    Bounds bounds ) {

  PRE08( polylineCount > 0, vertexCounts, xy,
         polygon, polygon->num_contours == 0,
         polygon->hole == 0, polygon->contour == 0, bounds );

  const double* xyp = xy;
  double xMinimum = xy[ 0 ];
  double xMaximum = xMinimum;
  double yMinimum = xy[ 1 ];
  double yMaximum = yMinimum;
  const int contours = polygon->num_contours = polylineCount;
  int contour = 0;
  int result = 0;
  polygon->hole = NEW( int, contours );
  polygon->contour =
    polygon->hole == 0 ? 0 : NEW( gpc_vertex_list, contours );
  result = polygon->contour != 0;

  for ( contour = 0; AND2( result, contour < contours ); ++contour ) {
    gpc_vertex_list* const vertexList = polygon->contour + contour;
    const int vertices = vertexList->num_vertices = vertexCounts[ contour ];
    gpc_vertex* v = vertexList->vertex = NEW( gpc_vertex, vertices );
    result = v != 0;

    if ( result ) {
      int vertex = 0;

      for ( vertex = 0; vertex < vertices; ++vertex ) {
        const double x = *xyp++;
        const double y = *xyp++;
        v->x = x;
        v->y = y;
        ++v;

        if ( x < xMinimum ) {
          xMinimum = x;
        } else if ( x > xMaximum ) {
          xMaximum = x;
        }

        if ( y < yMinimum ) {
          yMinimum = y;
        } else if ( y > yMaximum ) {
          yMaximum = y;
        }
      }
    }
  }

  bounds[ LONGITUDE ][ MINIMUM ] = xMinimum;
  bounds[ LONGITUDE ][ MAXIMUM ] = xMaximum;
  bounds[ LATITUDE  ][ MINIMUM ] = yMinimum;
  bounds[ LATITUDE  ][ MAXIMUM ] = yMaximum;

  if ( ! result ) {
    gpc_free_polygon( polygon );
    memset( polygon, 0, sizeof *polygon );
  }

  POST03( IS_BOOL( result ),
          isValidBounds( (const double (*)[2]) bounds ),
          IMPLIES_ELSE( result,
                        AND5( polygon->hole,
                              polygon->contour,
                              polygon->num_contours == polylineCount,
                              xy[ 0 ] == polygon->contour->vertex->x,
                              xy[ 1 ] == polygon->contour->vertex->y ),
                        IS_ZERO3( polygon->hole,
                                  polygon->contour,
                                  polygon->num_contours ) ) );

  return result;
}



/******************************************************************************
PURPOSE: ensureCorrectVertexOrder - Check and correct the vertex order to match
         hole designation to ESRI spec - i.e., hole vertices are in CCW order.
INPUTS:  gpc_polygon* polygon  Polygon to check.
OUTPUTS: gpc_polygon* polygon  Corrected polygon.
RETURNS: int 1 if corrected polygon now has strictly positive net area,
         else 0 which indicates that the area of the holes is >=
         the area of the non-holes.
NOTES:   http://mathworld.wolfram.com/PolygonArea.html
******************************************************************************/

int ensureCorrectVertexOrder( gpc_polygon* polygon ) {
  PRE0( polygon );
  int result = 0;
  double polygonArea = 0.0;
  const int contours = polygon->num_contours;
  int contour = 0;
  DEBUG2( fprintf( stderr, "ensureCorrectVertexOrder():\n" ); )

  for ( contour = 0; contour < contours; ++contour ) {
    const int isHole = polygon->hole[ contour ];
    const gpc_vertex_list vertex_list = polygon->contour[ contour ];
    const gpc_vertex* v = vertex_list.vertex;
    const int vertices = vertex_list.num_vertices;
    double contourArea = 0.0;
    int vertex = 0;

    for ( vertex = 0; vertex < vertices; ++vertex ) {
      const int vertexp1 = vertex + 1;
      const int vertex1 = vertexp1 < vertices ? vertexp1 : 0;
      const double x0 = v[ vertex  ].x;
      const double y0 = v[ vertex  ].y;
      const double x1 = v[ vertex1 ].x;
      const double y1 = v[ vertex1 ].y;
      const double triangleArea = x0 * y1 - x1 * y0;
      contourArea += triangleArea;
    }

    DEBUG2( fprintf( stderr, "hole = %d, contour area = %lf\n",
                     isHole, contourArea ); )

    /* Ensure correct vertex order (per ESRI spec, holes are CCW order): */

    if ( OR2( AND2(   isHole, contourArea < 0.0 ),
              AND2( ! isHole, contourArea > 0.0 ) ) ) {
      DEBUG2( fprintf( stderr, "reversingVertexList...\n" ); )
      reverseVertexList( polygon->contour + contour );
      contourArea = -contourArea;
    }

    CHECK( IMPLIES_ELSE( isHole, contourArea >= 0.0, contourArea <= 0.0 ) );
    polygonArea += -contourArea;
  }

  polygonArea *= 0.5;
  DEBUG2( fprintf( stderr, "polygon area (in square degrees) = %lf\n",
                   polygonArea ); )
  result = polygonArea > 0.0;
  POST0( IS_BOOL( result ) );
  DEBUG2( fprintf( stderr, "result = %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: reverseVertexList - Reverse the vertex order.
INPUTS:  gpc_vertex_list* vertex_list  Vertices to reverse.
OUTPUTS: gpc_vertex_list* vertex_list  Vertices in reversed order.
******************************************************************************/

void reverseVertexList( gpc_vertex_list* vertex_list ) {
  PRE03( vertex_list, vertex_list->vertex, vertex_list->num_vertices > 0 );
  int lower = 0;
  int upper = vertex_list->num_vertices - 1;
  gpc_vertex* const vertices = vertex_list->vertex;

  for ( lower = 0; lower < upper; ++lower, --upper ) {
    const double temp_x = vertices[ lower ].x;
    const double temp_y = vertices[ lower ].y;
    vertices[ lower ].x = vertices[ upper ].x;
    vertices[ lower ].y = vertices[ upper ].y;
    vertices[ upper ].x = temp_x;
    vertices[ upper ].y = temp_y;
  }
}



/******************************************************************************
PURPOSE: signedAreaOfPolygon - Signed area of a single contour of a polygon.
INPUTS:  const gpc_vertex_list* vertex_list  Vertex list to compute area of.
RETURNS: double signed area of polygon.
         Negative if vertices are in counter-clockwise order.
NOTES:   http://mathworld.wolfram.com/PolygonArea.html
******************************************************************************/

double signedAreaOfPolygon( const gpc_vertex_list* vertex_list ) {
  PRE03( vertex_list, vertex_list->num_vertices > 3, vertex_list->vertex );
  double result = 0.0;
  const int count = vertex_list->num_vertices;
  const gpc_vertex* const vertices = vertex_list->vertex;
  int index = 0;

  for ( index = 0; index < count; ++index ) {
    const int indexp1 = index + 1;
    const int index1 = indexp1 < count ? indexp1 : 0;
    const double triangleArea =
      vertices[ index  ].x * vertices[ index1 ].y -
      vertices[ index1 ].x * vertices[ index  ].y;
    result += triangleArea;
  }

  result *= 0.5;
  return result;
}



/******************************************************************************
PURPOSE: polygonArea - Compute positive area (in square meters) of polygon.
INPUTS:  const gpc_polygon* polygon  Polygon to check.
RETURNS: double area (in square meters) of polygon.
NOTES:   http://mathworld.wolfram.com/PolygonArea.html
******************************************************************************/

double polygonArea( const gpc_polygon* polygon ) {
  PRE0( polygon );
  static int initialized = 0;
  double result = 0.0;

  if ( ! initialized ) {
    const double majorSemiaxis = 6378137.0;
    const double minorSemiaxis = 6356752.3;
    const double lowerSecantLatitude = 30.0;
    const double upperSecantLatitude = 60.0;
    const double centerLongitude = -100.0;
    const double centerLatitude  = 40.0;
    initialize_albers( majorSemiaxis, minorSemiaxis,
                       lowerSecantLatitude, upperSecantLatitude,
                       centerLatitude, centerLongitude, 0.0, 0.0 );
    initialized = 1;
  }

  {
    const int contours = polygon->num_contours;
    int contour = 0;

    for ( contour = 0; contour < contours; ++contour ) {
      const int isHole = polygon->hole[ contour ];
      const gpc_vertex_list vertex_list = polygon->contour[ contour ];
      const gpc_vertex* v = vertex_list.vertex;
      const int vertices = vertex_list.num_vertices;
      double contourArea = 0.0;
      int vertex = 0;

      for ( vertex = 0; vertex < vertices; ++vertex ) {
        const int vertexp1 = vertex + 1;
        const int vertex1 = vertexp1 < vertices ? vertexp1 : 0;
        const double longitude0 = v[ vertex  ].x;
        const double latitude0  = v[ vertex  ].y;
        const double longitude1 = v[ vertex1 ].x;
        const double latitude1  = v[ vertex1 ].y;
        double x0 = 0;
        double y0 = 0.0;
        double x1 = 0.0;
        double y1 = 0.0;
        double triangleArea = 0.0;
        project_albers( longitude0, latitude0, &x0, &y0 );
        project_albers( longitude1, latitude1, &x1, &y1 );
        triangleArea = x0 * y1 - x1 * y0;
        contourArea += triangleArea;
      }

      DEBUG2( fprintf( stderr, "hole = %d, contour area = %lf\n",
                       isHole, contourArea ); )

      /* Ensure correct vertex order (per ESRI spec, holes are CCW order): */

      if ( OR2( AND2(   isHole, contourArea < 0.0 ),
                AND2( ! isHole, contourArea > 0.0 ) ) ) {
        DEBUG2( fprintf( stderr, "reversingVertexList...\n" ); )
        reverseVertexList( polygon->contour + contour );
        contourArea = -contourArea;
      }

      CHECK( IMPLIES_ELSE( isHole, contourArea >= 0.0, contourArea <= 0.0 ) );
      result += -contourArea;
    }

    result *= 0.5;
  }

  DEBUG2( fprintf( stderr, "polygonArea = %lf\n", result ); )
  POST0( result >= 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: polygonPerimeter - Compute perimeter (in meters) of polygon.
INPUTS:  const gpc_polygon* polygon  Polygon to check.
RETURNS: double perimeter (in meters) of polygon.
******************************************************************************/

double polygonPerimeter( const gpc_polygon* polygon ) {
  PRE0( polygon );
  static int initialized = 0;
  double result = 0.0;

  if ( ! initialized ) {
    const double majorSemiaxis = 6378137.0;
    const double minorSemiaxis = 6356752.3;
    const double lowerSecantLatitude = 30.0;
    const double upperSecantLatitude = 60.0;
    const double centerLongitude = -100.0;
    const double centerLatitude  = 40.0;
    initialize_albers( majorSemiaxis, minorSemiaxis,
                       lowerSecantLatitude, upperSecantLatitude,
                       centerLatitude, centerLongitude, 0.0, 0.0 );
    initialized = 1;
  }

  {
    const int contours = polygon->num_contours;
    int contour = 0;

    for ( contour = 0; contour < contours; ++contour ) {
      const gpc_vertex_list vertexList = polygon->contour[ contour ];
      const gpc_vertex* v = vertexList.vertex;
      const int vertices = vertexList.num_vertices;
      int vertex = 0;
      const double longitude0 = v[ 0 ].x;
      const double latitude0  = v[ 0 ].y;
      double x0 = 0;
      double y0 = 0.0;
      project_albers( longitude0, latitude0, &x0, &y0 );

      for ( vertex = 1; vertex < vertices; ++vertex ) {
        const double longitude = v[ vertex ].x;
        const double latitude  = v[ vertex ].y;
        double x = 0.0;
        double y = 0.0;
        project_albers( longitude, latitude, &x, &y );

        {
          const double deltaX = x - x0;
          const double deltaY = y - y0;
          const double segmentLength =
            sqrt( SQUARE( deltaX ) + SQUARE( deltaY ) );
          result += segmentLength;
          x0 = x;
          y0 = y;
        }
      }
    }
  }

  DEBUG2( fprintf( stderr, "polygonPerimeter = %le\n", result ); )
  POST0( result >= 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: writePolygonsToShapefile - Write clipped polygons to Shapefile
         (shx, shp).
INPUTS:  SHPHandle outputFile                  File to write to.
         int isPolyline                        Is shapefile polyline?
         int count                             Number of shapes.
         const PolygonShape polygons[ count ]  Shapes to write.
RETURNS: int 1 if successful, else 0 failure message is printed.
NOTES:   Uses external shapelib API and GPC API polygon data type.
******************************************************************************/

int writePolygonsToShapefile( SHPHandle outputFile, int isPolyline,
                              int count, const PolygonShape polygons[] ) {

  PRE03( outputFile, count > 0, polygons );

  int result = 0;
  const int maximumParts = maximumPolygonContours( count, polygons );
  const int maximumVertices =
    maximumPolygonVertices( count, polygons ) + maximumParts;
  double* x = NEW( double, maximumVertices * 2 );
  int* starts = maximumParts == 1 ? 0 : NEW( int, maximumParts );

  if ( AND2( x, IMPLIES( maximumParts > 1, starts ) ) ) {
    double* const y = x + maximumVertices;
    const int type = isPolyline ? SHPT_ARC : SHPT_POLYGON;
    int index = 0;

    for ( index = 0; index < count; ++index ) {
      const gpc_polygon polygon = polygons[ index ].polygon;
      const int parts = polygon.num_contours;
      SHPObject* shape = 0;
      CHECK( parts >= 1 );

      if ( parts == 1 ) {
        const int vertices = polygon.contour->num_vertices + ! isPolyline;
        CHECK( ! polygon.hole[ 0 ] ); /* Don't allow single holes. */
        copyPolygonVertices( &polygon, ! isPolyline, 0, x, y );
        shape = SHPCreateSimpleObject( type, vertices, x, y, 0 );
      } else {
        const int vertices =
          polygonVertexCount( &polygon ) + parts * ( ! isPolyline );
        copyPolygonVertices( &polygon, ! isPolyline, starts, x, y );
        shape = SHPCreateObject( type, index, parts, starts,
                                 0, vertices, x, y, 0, 0 );
      }

      if ( shape ) {
        DEBUG( const int entityNumber = )
        SHPWriteObject( outputFile, -1, shape );
        DEBUG( fprintf( stderr,
                        "wrote shape # %d with %d parts and %d vertices\n",
                        entityNumber, shape->nParts, shape->nVertices ); )
        DEBUG( printShape( shape ); )
        SHPDestroyObject( shape ), shape = 0;
      } else {
        fprintf( stderr, "\nFailed to create shape.\n" );
        index = count; /* Stop looping. */
      }
    }

    result = index == count;
  }

  FREE( x );
  FREE( starts );

  DEBUG( fprintf(stderr, "writePolygonsToShapefile() returning %d\n", result);)
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: copySubsetShapefile - Copy masked subset of shapes to Shapefile
         (shx, shp).
INPUTS:  SHPHandle inputFile       File to read from.
         SHPHandle outputFile      File to write to.
         const int count           Number of shapes.
         const char mask[ count ]  If mask[i] == 1 write shape i else skip.
RETURNS: int number of subset shapes written if successful,
         else 0 failure message is printed.
NOTES:   Uses external shapelib API.
******************************************************************************/

int copySubsetShapefile( SHPHandle inputFile, SHPHandle outputFile,
                         const int count, const char* mask ) {
  PRE04( inputFile, outputFile, count > 0, mask );

  int result = 0;
  int index = 0;
  DEBUG( fprintf( stderr, "copySubsetShapefile( count = %d, mask = %p )\n",
                 count, mask ); )

  for ( index = 0; index < count; ++index ) {
    const char m = mask[ index ];

    if ( m ) {
      SHPObject* shape = SHPReadObject( inputFile, index );

      if ( ! shape ) {
        fprintf( stderr, "\nFailed to read shape #%d.\n", result );
        index = count; /* Stop looping. */
        result = 0;
      } else {

        if ( SHPWriteObject( outputFile, result, shape ) >= 0 ) {
          ++result;
        } else {
          fprintf( stderr, "\nFailed to write shape #%d.\n", result );
          index = count; /* Stop looping. */
          result = 0;
        }

        SHPDestroyObject( shape );
        shape = 0;
      }
    }
  }

  DEBUG( fprintf( stderr, "copySubsetShapefile() returning %d\n", result ); )
  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: writePolygonDBF - Write subset of DBF file for clipped polygons.
INPUTS:  const char* inputFileName        File name of input dbf.
         DBFHandle outputFile             DBF file to write to.
         int offset                       Record offset to starting writing.
         int count                        Number of shapes.
         const char* mask                 0 or mask (0,1) to filter shapes by.
         const PolygonShape polygons[ count ]  Shapes to write.
RETURNS: int 1 if successful, else 0 failure message is printed.
NOTES:   Uses static file-global table[] at the top of this file and
         external shapelib API.
******************************************************************************/

int writePolygonDBF( const char* inputFileName, DBFHandle outputFile,
                     int offset, int count, const char* mask,
                     const PolygonShape polygons[] ) {

  PRE06( inputFileName, *inputFileName, outputFile, offset >= 0, count > 0,
         polygons );

  int result = 0;
  DBFHandle inputFile = DBFOpen( inputFileName, "rb" );

  if ( inputFile ) {
    const int isSoil = strstr( inputFileName, "soil_" ) != 0;
    const int isStreamTemperature =
      strstr( inputFileName, "stream_temperature_" ) != 0;
    const int isStreamTemperatureLine =
      strstr( inputFileName, "stream_temperature_line_" ) != 0;
    const int isTide = strstr( inputFileName, "tide_" ) != 0;
    const int isCoastalVulnerability =
      strstr( inputFileName, "coastal_vulnerability" ) != 0;
    const int isCMAQ = strstr( inputFileName, "_cmaq_" ) != 0;
    const int isGreenspaceHousing =
      strstr( inputFileName, "greenspace_housing" ) != 0;
    const double F_MISSING = isSoil ? -99.999 : -9999.0;
    int outputIndex = -1; /* Index into table of first column info. */
    int longitudeColumn = -1;
    const int outputColumns =
      defineDBFColumns( inputFileName, offset == 0,
                        &outputIndex, &longitudeColumn, 0, 0, 0, 0,
                        outputFile );
    int ok = outputColumns > 0;
    int index = 0;

    /* For each row/shape: */

    for ( index = 0; AND2( ok, index < count ); ++index ) {
      const int id = polygons[ index ].id; /* Input record/row. */

      if ( OR2( mask == 0, mask[ id ] ) ) {
        const int record = offset + index; /* Output record/row. */
        double areaInSquareMeters = 0.0;   /* Of current clipped polygon. */
        int outputColumn = 0;

        /* Write each column value in order shown in file global table: */

        for ( outputColumn = 0; AND2( ok, outputColumn < outputColumns );
              ++outputColumn ) {
          const ColumnEntry* const columnEntry =
            table + outputIndex + outputColumn;
          const int inputColumn = columnEntry->inputColumn;
          const char* const columnName = columnEntry->columnName;

          if ( inputColumn > -1 ) { /* Just map input value to output value: */
            const int columnType = columnEntry ->columnType;

            if ( columnType == FTDouble ) {
              const int filterNegatives =
                AND6( inputColumn != longitudeColumn,
                      ! isCoastalVulnerability ,
                      ! isCMAQ,
                      OR2( ! isStreamTemperature, isStreamTemperatureLine ),
                      ! isTide,
                      ! isGreenspaceHousing );
              ok = copyDoubleAttribute( inputFile, id, inputColumn,
                                        outputFile, record, outputColumn,
                                        filterNegatives, F_MISSING,
                                        columnEntry->offset,
                                        columnEntry->scale );
            } else if ( columnType == FTInteger ) {
              ok = copyIntegerAttribute( inputFile, id, inputColumn,
                                         outputFile, record, outputColumn );
            } else {
              CHECK( columnType == FTString );
              ok = copyStringAttribute( inputFile, id, inputColumn,
                                        outputFile, record, outputColumn );
            }
          } else if ( OR3( ! strcmp( columnName, "ACRES" ),
                           ! strcmp( columnName, "HECTARES" ),
                           ! strcmp( columnName, "AREA_SQKM" ) ) ) {
            double outputValue = 0.0; /* Compute and write subset area: */

            if ( areaInSquareMeters == 0.0 ) { /* Compute once for this row: */
              areaInSquareMeters = polygonArea( &polygons[ index ].polygon );
            }

            if ( ! strcmp( columnName, "ACRES" ) ) {
              const double SQUARE_METERS_TO_ACRES = 0.000247105381;
              const double areaInAcres =
                areaInSquareMeters * SQUARE_METERS_TO_ACRES;
              outputValue = areaInAcres;
            } else if ( ! strcmp( columnName, "HECTARES" ) ) {
              const double SQUARE_METERS_TO_HECTARES = 1e-4;
              const double areaInHectares =
                areaInSquareMeters * SQUARE_METERS_TO_HECTARES;
              outputValue = areaInHectares;
            } else {
              const double SQUARE_METERS_TO_SQUARE_KILOMETERS = 1e-6;
              const double areaInSquareKilometers =
                areaInSquareMeters * SQUARE_METERS_TO_SQUARE_KILOMETERS;
              outputValue = areaInSquareKilometers;
              CHECK( ! strcmp( columnName, "AREA_SQKM" ) );
            }

            ok = DBFWriteDoubleAttribute( outputFile, record, outputColumn,
                                          outputValue );
          } else if ( ! strcmp( columnName, "POP_SQKM" ) ) { /* Special case:*/
            const double shapeArea = DBFReadDoubleAttribute( inputFile, id, 3);
            const int population = DBFReadIntegerAttribute( inputFile, id, 4 );
            const double SQUARE_METERS_TO_SQUARE_KILOMETERS = 1e-6;
            const double populationPerSquareKilometer =
              population / ( shapeArea * SQUARE_METERS_TO_SQUARE_KILOMETERS );
            const double areaInSquareKilometers =
              areaInSquareMeters * SQUARE_METERS_TO_SQUARE_KILOMETERS;
            const int subsetPopulation = (int)
              ( populationPerSquareKilometer * areaInSquareKilometers + 0.5 );
            ok = DBFWriteDoubleAttribute( outputFile, record, 3,
                                          populationPerSquareKilometer );

            /* Also write SUBSET_POP: */

            ok = AND2( ok,
                       DBFWriteIntegerAttribute( outputFile, record, 4,
                                                 (int) subsetPopulation ) );
          } else if ( strstr( columnName, "0PKM" ) ) {

            /*
             * Find corresponding population and compute and write
             * population per square kilometer:
             */

            char populationName[ 16 ] = "";
            int populationIndex = 0;
            char* c = 0;
            memset( populationName, 0, sizeof populationName );
            strncpy( populationName, columnName,
                     sizeof populationName / sizeof *populationName - 1 );
            CHECK( strlen( populationName ) <
                   sizeof populationName / sizeof *populationName );
            c = strstr( populationName, "0PKM" );
            CHECK( c );
            strcpy( c, "0POP" );

            while ( AND2( table[ populationIndex ].columnName,
                          strcmp( table[ populationIndex ].columnName,
                                  populationName ) ) ) {
              ++populationIndex;
            }

            CHECK( ! strcmp( table[ populationIndex ].columnName,
                             populationName ) );
            {
              const int populationColumn = table[ populationIndex ].inputColumn;
              const int population =
                DBFReadIntegerAttribute( inputFile, id, populationColumn );
              const double shapeArea = DBFReadDoubleAttribute(inputFile, id, 1);
              const double SQUARE_METERS_TO_SQUARE_KILOMETERS = 1e-6;
              const double populationPerSquareKilometer =
                population / ( shapeArea * SQUARE_METERS_TO_SQUARE_KILOMETERS);
              ok = DBFWriteDoubleAttribute( outputFile, record, outputColumn,
                                            populationPerSquareKilometer );
            }

          } else if ( ! strcmp( columnName, "LENGTH_KM" ) ) {
            const double METERS_TO_KILOMETERS = 1e-3;
            const double outputValue =
              METERS_TO_KILOMETERS *
              polygonPerimeter( &polygons[index].polygon );
            ok = DBFWriteDoubleAttribute( outputFile, record, outputColumn,
                                          outputValue );
          } else if ( ! strcmp( columnName, "TOT_YKGKMY" ) ) {
            const int isMRB = strstr( inputFileName, "sparrow_2002_mrb" ) != 0;
            const int isMRB2 =
              strstr( inputFileName, "load_estuary_sparrow_2002_mrb2" ) != 0;
            const int isNonAtlantic =
              ! strstr( inputFileName, "load_estuary_sparrow_1992_atlantic" );
            const int totalLoadColumn =
              isMRB2 ? 8 : isMRB ? 4 : 16 + isNonAtlantic;
            const int areaKm2Column = OR2( isMRB, isMRB2 ) ? 3 : 11;
            const double totalLoad =
              DBFReadDoubleAttribute( inputFile, id, totalLoadColumn );
            const double areaKm2 =
              DBFReadDoubleAttribute( inputFile, id, areaKm2Column );
            const double totalYield = totalLoad / areaKm2;
            ok = DBFWriteDoubleAttribute( outputFile, record, outputColumn,
                                          totalYield );
          } else if ( AND2( strstr( inputFileName, "estuary_cmaq" ),
                            ! strcmp( columnName, "UNITS" ) ) ) {
            ok = DBFWriteStringAttribute( outputFile, record, outputColumn,
                                          "kgN/ha/year" );
          } else if ( AND2( strstr( inputFileName, "estuary_cmaq" ),
                            ! strcmp( columnName, "TOTN_2007" ) ) ) {
            const double sum =
              DBFReadDoubleAttribute( inputFile, id, 63 ) +
              DBFReadDoubleAttribute( inputFile, id, 64 ) +
              DBFReadDoubleAttribute( inputFile, id, 65 ) +
              DBFReadDoubleAttribute( inputFile, id, 66 ) +
              DBFReadDoubleAttribute( inputFile, id, 67 ) +
              DBFReadDoubleAttribute( inputFile, id, 68 ) +
              DBFReadDoubleAttribute( inputFile, id, 69 ) +
              DBFReadDoubleAttribute( inputFile, id, 70 ) +
              DBFReadDoubleAttribute( inputFile, id, 71 ) +
              DBFReadDoubleAttribute( inputFile, id, 96 ) +
              DBFReadDoubleAttribute( inputFile, id, 97 ) +
              DBFReadDoubleAttribute( inputFile, id, 98 );
            ok = DBFWriteDoubleAttribute(outputFile, record, outputColumn, sum);
          } else if ( AND2( strstr( inputFileName, "estuary_cmaq" ),
                            ! strcmp( columnName, "DRYN_2007" ) ) ) {
            const double sum =
              DBFReadDoubleAttribute( inputFile, id, 147 ) +
              DBFReadDoubleAttribute( inputFile, id, 148 ) +
              DBFReadDoubleAttribute( inputFile, id, 149 ) +
              DBFReadDoubleAttribute( inputFile, id, 150 ) +
              DBFReadDoubleAttribute( inputFile, id, 151 ) +
              DBFReadDoubleAttribute( inputFile, id, 152 ) +
              DBFReadDoubleAttribute( inputFile, id, 153 ) +
              DBFReadDoubleAttribute( inputFile, id, 154 ) +
              DBFReadDoubleAttribute( inputFile, id, 155 ) +
              DBFReadDoubleAttribute( inputFile, id, 180 ) +
              DBFReadDoubleAttribute( inputFile, id, 181 ) +
              DBFReadDoubleAttribute( inputFile, id, 182 );
            ok = DBFWriteDoubleAttribute(outputFile, record, outputColumn, sum);
          } else if ( AND2( strstr( inputFileName, "estuary_cmaq" ),
                            ! strcmp( columnName, "WETN_2007" ) ) ) {
            const double sum =
              DBFReadDoubleAttribute( inputFile, id, 231 ) +
              DBFReadDoubleAttribute( inputFile, id, 232 ) +
              DBFReadDoubleAttribute( inputFile, id, 233 ) +
              DBFReadDoubleAttribute( inputFile, id, 234 ) +
              DBFReadDoubleAttribute( inputFile, id, 235 ) +
              DBFReadDoubleAttribute( inputFile, id, 236 ) +
              DBFReadDoubleAttribute( inputFile, id, 237 ) +
              DBFReadDoubleAttribute( inputFile, id, 238 ) +
              DBFReadDoubleAttribute( inputFile, id, 239 ) +
              DBFReadDoubleAttribute( inputFile, id, 264 ) +
              DBFReadDoubleAttribute( inputFile, id, 265 ) +
              DBFReadDoubleAttribute( inputFile, id, 266 );
            ok = DBFWriteDoubleAttribute(outputFile, record, outputColumn, sum);
          } else if ( AND2( strstr( inputFileName, "estuary_cmaq" ),
                            ! strcmp( columnName, "TOTN_2008" ) ) ) {
            const double sum =
              DBFReadDoubleAttribute( inputFile, id, 72 ) +
              DBFReadDoubleAttribute( inputFile, id, 73 ) +
              DBFReadDoubleAttribute( inputFile, id, 74 ) +
              DBFReadDoubleAttribute( inputFile, id, 75 ) +
              DBFReadDoubleAttribute( inputFile, id, 76 ) +
              DBFReadDoubleAttribute( inputFile, id, 77 ) +
              DBFReadDoubleAttribute( inputFile, id, 78 ) +
              DBFReadDoubleAttribute( inputFile, id, 79 ) +
              DBFReadDoubleAttribute( inputFile, id, 80 ) +
              DBFReadDoubleAttribute( inputFile, id, 99 ) +
              DBFReadDoubleAttribute( inputFile, id, 100 ) +
              DBFReadDoubleAttribute( inputFile, id, 101 );
            ok = DBFWriteDoubleAttribute(outputFile, record, outputColumn, sum);
          } else if ( AND2( strstr( inputFileName, "estuary_cmaq" ),
                            ! strcmp( columnName, "DRYN_2008" ) ) ) {
            const double sum =
              DBFReadDoubleAttribute( inputFile, id, 156 ) +
              DBFReadDoubleAttribute( inputFile, id, 157 ) +
              DBFReadDoubleAttribute( inputFile, id, 158 ) +
              DBFReadDoubleAttribute( inputFile, id, 159 ) +
              DBFReadDoubleAttribute( inputFile, id, 160 ) +
              DBFReadDoubleAttribute( inputFile, id, 161 ) +
              DBFReadDoubleAttribute( inputFile, id, 162 ) +
              DBFReadDoubleAttribute( inputFile, id, 163 ) +
              DBFReadDoubleAttribute( inputFile, id, 164 ) +
              DBFReadDoubleAttribute( inputFile, id, 183 ) +
              DBFReadDoubleAttribute( inputFile, id, 184 ) +
              DBFReadDoubleAttribute( inputFile, id, 185 );
            ok = DBFWriteDoubleAttribute(outputFile, record, outputColumn, sum);
          } else if ( AND2( strstr( inputFileName, "estuary_cmaq" ),
                            ! strcmp( columnName, "WETN_2008" ) ) ) {
            const double sum =
              DBFReadDoubleAttribute( inputFile, id, 240 ) +
              DBFReadDoubleAttribute( inputFile, id, 241 ) +
              DBFReadDoubleAttribute( inputFile, id, 242 ) +
              DBFReadDoubleAttribute( inputFile, id, 243 ) +
              DBFReadDoubleAttribute( inputFile, id, 244 ) +
              DBFReadDoubleAttribute( inputFile, id, 245 ) +
              DBFReadDoubleAttribute( inputFile, id, 246 ) +
              DBFReadDoubleAttribute( inputFile, id, 247 ) +
              DBFReadDoubleAttribute( inputFile, id, 248 ) +
              DBFReadDoubleAttribute( inputFile, id, 267 ) +
              DBFReadDoubleAttribute( inputFile, id, 268 ) +
              DBFReadDoubleAttribute( inputFile, id, 269 );
            ok = DBFWriteDoubleAttribute(outputFile, record, outputColumn, sum);
          } else if ( ! strcmp( columnName, "YYYYDDD1" ) ) { /* HMS smoke: */

            /*
             * Start / End timestamp may be either "HHMM" or "YYYYDDD HHMM".
             * If only HHMM then parse YYYYMMDD from file name and
             * convert it to YYYYDDD. Then write column as an integer.
             */

            const char* const timestamp =
              DBFReadStringAttribute( inputFile, id, 1 );
            int yyyyddd = 0;
            int hhmm = 0;
            ok = convertTimestamp( inputFileName, timestamp, &yyyyddd, &hhmm );

            if ( ok ) {
              ok = DBFWriteIntegerAttribute( outputFile, record, outputColumn,
                                             yyyyddd );
            }
          } else if ( ! strcmp( columnName, "HHMM1" ) ) { /* HMS smoke: */
            const char* const timestamp =
              DBFReadStringAttribute( inputFile, id, 1 );
            int yyyyddd = 0;
            int hhmm = 0;
            ok = convertTimestamp( inputFileName, timestamp, &yyyyddd, &hhmm );

            if ( ok ) {
              ok = DBFWriteIntegerAttribute( outputFile, record, outputColumn,
                                             hhmm );
            }
          } else if ( ! strcmp( columnName, "YYYYDDD2" ) ) { /* HMS smoke: */
            const char* const timestamp =
              DBFReadStringAttribute( inputFile, id, 2 );
            int yyyyddd = 0;
            int hhmm = 0;
            ok = convertTimestamp( inputFileName, timestamp, &yyyyddd, &hhmm );

            if ( ok ) {
              ok = DBFWriteIntegerAttribute( outputFile, record, outputColumn,
                                             yyyyddd );
            }
          } else if ( ! strcmp( columnName, "HHMM2" ) ) { /* HMS smoke: */
            const char* const timestamp =
              DBFReadStringAttribute( inputFile, id, 2 );
            int yyyyddd = 0;
            int hhmm = 0;
            ok = convertTimestamp( inputFileName, timestamp, &yyyyddd, &hhmm );

            if ( ok ) {
              ok = DBFWriteIntegerAttribute( outputFile, record, outputColumn,
                                             hhmm );
            }
          } else if ( ! strcmp( columnName, "DENS_UGM3" ) ) { /* HMS smoke: */

            /*
             * After 2022-07-18, density is no longer numeric but a string:
             * Light, Medium, Heavy.
             * Detect this and convert to numeric value in ug/m3:
             * Light = 5
             * Medium = 10
             * Heavy = 30
             */

            double density = 0.0;
            const char* const tag = "hms_smoke";
            const char* const match = strstr( inputFileName, tag );
            const int fileYYYYMMDD =
              match ? atoi( match + strlen( tag ) ) : 0;

            if ( fileYYYYMMDD < 20220719 ) {
              density = DBFReadDoubleAttribute( inputFile, id, 3 );
            } else {
              const char* const densityString =
                DBFReadStringAttribute( inputFile, id, 3 );

              if ( densityString ) {

                if ( ! strcmp( densityString, "Light" ) ) {
                  density = 5.0;
                } else if ( ! strcmp( densityString, "Medium" ) ) {
                  density = 10.0;
                } else if ( ! strcmp( densityString, "Heavy" ) ) {
                  density = 30.0;
                }
              }
            }

            if ( density > 0.0 ) {
              ok = DBFWriteIntegerAttribute( outputFile, record, outputColumn,
                                             (int) density );
            }
          }

          if ( ! ok ) {
            fprintf( stderr,
                     "Failed to write row %d column %d (%s) to dbf file.\n",
                     record, outputColumn, columnName );
          }
        }
      }
    }

    result = ok;
    DBFClose( inputFile ), inputFile = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "Failed to write to dbf file.\n" );
  }

  DEBUG2( fprintf( stderr, "writePolygonDBF() returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeSubsetDBF - Write subset of DBF file for unmasked rows in bounds
         or with specified huc or estcode.
INPUTS:  const char* inputFileName  File name of input dbf.
         const Bounds          Bounds to subset rows of DBF containing
                               LONGITUDE, LATITUDE columns.
         const long long huc   HUC to filter by (if huc > 0).
         const char* const estcode  Or ESTCODE to filter by or "all" (if != 0).
         const int siteIdCount      0 or number of unique siteIds.
         const int siteIds[ siteIdCount ]  0 or siteIds[] to filter by.
         const int count       Number of file rows expected.
         char mask[ count ]    mask[ row ] == 1 means test row for
                               optional additional filtering by
                               huc/estcode/bounds, else omit row.
         DBFHandle outputFile  DBF file to write to.
OUTPUTS: char mask[ count ]    Filtered mask by huc/estcode/bounds.
RETURNS: int > 0 number of output rows if successful else 0 and
         failure message is printed to stderr.
NOTES:   Uses static file-global table[] at the top of this file and
         external shapelib API.
******************************************************************************/

int writeSubsetDBF( const char* inputFileName, const Bounds bounds,
                    const long long huc, const char* const estcode,
                    const int siteIdCount, const int siteIds[],
                    const int count,
                    char mask[], DBFHandle outputFile ) {

  PRE08( inputFileName,
         *inputFileName,
         isValidBounds( (const double (*)[2]) bounds ),
         huc >= 0,
         IMPLIES_ELSE( siteIdCount > 0,
                       AND4( siteIds, siteIds[ 0 ], siteIds[ siteIdCount - 1 ],
                             siteIds[ 0 ] <= siteIds[ siteIdCount - 1 ] ),
                       siteIds == 0 ),
         count > 0,
         mask,
         outputFile );

  int result = 0;
  DBFHandle inputFile = DBFOpen( inputFileName, "rb" );

  if ( inputFile ) {
    const int inputRecords = DBFGetRecordCount( inputFile );

    if ( inputRecords != count ) {
      fprintf( stderr, "\nUnmatched rows in dbf file: "
               "actual = %d, expected = %d.\n",
              inputRecords, count );
    } else {
      const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
      const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
      const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
      const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
      const int isNCA =
        OR2( strstr( inputFileName, "sediment_nca"  ),
             strstr( inputFileName, "SEDIMENT_NCA"  ) );
      const int isSoil = strstr( inputFileName, "soil_" ) != 0;
      const int isStreamTemperature =
        strstr( inputFileName, "stream_temperature_" ) != 0;
      const int isTide = strstr( inputFileName, "tide_" ) != 0;
      const int isCoastalVulnerability =
        strstr( inputFileName, "coastal_vulnerability" ) != 0;
      const int isCMAQ = strstr( inputFileName, "_cmaq_" ) != 0;
      const int isGreenspaceHousing =
        strstr( inputFileName, "greenspace_housing" ) != 0;
      const int isFlowlinesWatershed =
        AND2( estcode,
              strstr( inputFileName, "flowlines_puget_sound_watershed" ) );
      const double F_MISSING = isSoil ? -99.999 : -9999.0;
      const int inputRecords = DBFGetRecordCount( inputFile );
      int inputRecord = 0;
      int outputRecord = 0;
      int outputIndex = -1;
      int longitudeColumn = -1;
      int latitudeColumn = -1;
      int hucColumn = -1;
      int estcodeColumn = -1;
      int siteIdColumn = -1;
      const int outputColumns =
        defineDBFColumns( inputFileName, 1,
                          &outputIndex, &longitudeColumn, &latitudeColumn,
                          huc > 0 ? &hucColumn : 0,
                          estcode ? &estcodeColumn : 0,
                          siteIds ? &siteIdColumn : 0,
                          outputFile );
      int ok = outputColumns > 0;

      DEBUG( fprintf( stderr, "huc = %lld, estcode = '%s', mask = %p, "
                      "isStreamTemperature = %d, "
                      "table index = %d, outputColumns = %d, "
                      "longitudeColumnn = %d, latitudeColumn = %d, "
                      "hucColumn = %d, estcodeColumn = %d, "
                      "siteIdColumn = %d\n",
                      huc, estcode ? estcode : "", mask,
                      isStreamTemperature, outputIndex, outputColumns,
                      longitudeColumn, latitudeColumn, hucColumn,
                      estcodeColumn, siteIdColumn ); )

      /* Read rows filtered by mask/huc/estcode/bounds: */

      for ( inputRecord = 0; AND2( ok, inputRecord < inputRecords );
            ++inputRecord ) {
        const double longitude =
          longitudeColumn < 0 ? -9999.0
          : DBFReadDoubleAttribute( inputFile, inputRecord, longitudeColumn );
        const double latitude =
          latitudeColumn < 0 ? -9999.0
          : DBFReadDoubleAttribute( inputFile, inputRecord, latitudeColumn );
        const long long hucId =
          hucColumn < 0 ? 0
          : (long long) DBFReadDoubleAttribute(inputFile,inputRecord,hucColumn);
        const char* const estcodeValue =
          estcodeColumn < 0 ? 0
          : DBFReadStringAttribute( inputFile, inputRecord, estcodeColumn );
        const int siteId =
          siteIdColumn < 0 ? 0
          : DBFReadIntegerAttribute( inputFile, inputRecord, siteIdColumn );
        const int m = mask[ inputRecord ];
        int inSubset =
          AND2( m,
                huc > 0 ? hucId == huc
                : AND3( estcode, estcodeValue, longitudeColumn >= 0 ) ?
                    AND3( OR2( ! strcmp( estcode, "all" ),
                               ! strcmp( estcodeValue, estcode ) ),
                          IN_RANGE( longitude,
                                    longitudeMinimum, longitudeMaximum ),
                          IN_RANGE( latitude,
                                    latitudeMinimum,  latitudeMaximum  ) )
                : isFlowlinesWatershed ?
                  matchesWithUnderscores( estcode, estcodeValue )
                : AND2( estcode, estcodeValue ) ?
                    OR2( ! strcmp( estcode, "all" ),
                         ! strcmp( estcodeValue, estcode ) )
                : IMPLIES( longitudeColumn != -1,
                           AND2( IN_RANGE( longitude,
                                           longitudeMinimum, longitudeMaximum ),
                                 IN_RANGE( latitude,
                                           latitudeMinimum, latitudeMaximum))));

        if ( AND2( inSubset, siteIds ) ) { /* Search sorted siteIds[]: */
          int index = 0;

          while ( AND2( index < siteIdCount, siteId > siteIds[ index ] ) ) {
            ++index;
          }

          inSubset = AND2( index < siteIdCount, siteId == siteIds[ index ] );
        }

        mask[ inputRecord ] = inSubset;

        if ( inSubset ) {
          int outputColumn = 0;

          /* Write each column value in order shown in table at top of file: */

          for ( outputColumn = 0; AND2( ok, outputColumn < outputColumns );
                ++outputColumn ) {
            const ColumnEntry* const columnEntry =
              table + outputIndex + outputColumn;
            const int inputColumn = columnEntry->inputColumn;
            const char* const columnName = columnEntry->columnName;
            const int columnType = columnEntry ->columnType;
            CHECK( inputColumn > -1 );

            if ( columnType == FTDouble ) {
              const int filterNegatives =
                AND6( inputColumn != longitudeColumn,
                      ! isCoastalVulnerability ,
                      ! isCMAQ,
                      ! isStreamTemperature,
                      ! isTide,
                      ! isGreenspaceHousing );
              const double offset = columnEntry ->offset;
              const double scale  = columnEntry ->scale;
              const int isNCATOC =
                AND3( isNCA, ! strcmp( columnName, "TOC_%" ),
                     strstr( inputFileName, "sediment_nca_2015" ) == 0 );

              if ( OR3( offset != 0.0, scale != 1.0, isNCATOC ) ) {
                double value =
                  DBFReadDoubleAttribute( inputFile, inputRecord, inputColumn );

                if ( AND2( filterNegatives, value < 0.0 ) ) {
                  value = F_MISSING;
                } else if ( isNCATOC ) {
                  const char* const units =
                    DBFReadStringAttribute( inputFile, inputRecord,
                                            inputColumn + 1 ); /* TOC_UNITS */

                  if ( OR2( ! strcmp( units, "ppm" ),
                            ! strcmp( units, "ug/g" ) ) ) {
                    value *= 1e-4; /* Convert TOC to %. */
                  }
                } else {
                  value += offset;
                  value *= scale;
                }

                ok = DBFWriteDoubleAttribute( outputFile,
                                              outputRecord, outputColumn,
                                              value );
              } else {
                ok = copyDoubleAttribute( inputFile, inputRecord, inputColumn,
                                          outputFile, outputRecord, outputColumn,
                                          filterNegatives, F_MISSING,
                                          columnEntry->offset,
                                          columnEntry->scale );
              }
            } else if ( columnType == FTInteger ) {
              ok = copyIntegerAttribute( inputFile, inputRecord, inputColumn,
                                         outputFile, outputRecord, outputColumn);
            } else {
              CHECK( columnType == FTString );
              ok = copyStringAttribute( inputFile, inputRecord, inputColumn,
                                        outputFile, outputRecord, outputColumn );
            }

            DEBUG( if ( ! ok ) fprintf( stderr, "Table entry %d %s\n",
                                        outputIndex, columnName ); )

          } /* End loop on output columns. */

          ++outputRecord;
        } /* End if inSubset. */
      } /* End loop on input records. */

      if ( ! ok ) {
        fprintf( stderr, "Failed to write to dbf file.\n" );
      } else {
        result = outputRecord;
      }

    } /* If process records. */

    DBFClose( inputFile ), inputFile = 0;
  }

  DEBUG2( fprintf( stderr, "writeSubsetDBF() returning %d\n", result ); )
  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: getRowsDBF - Get number of rows in a DBF file.
INPUTS:  const char* baseFileName  Base name of DBF file.
RETURNS: int > 0 number of rows if successful else 0 and
         failure message is printed to stderr.
******************************************************************************/

int getRowsDBF( const char* baseFileName ) {

  PRE02( baseFileName, *baseFileName );
  int result = 0;
  DBFHandle inputFile = DBFOpen( baseFileName, "rb" );

  if ( inputFile ) {
    const int inputRecords = DBFGetRecordCount( inputFile );

    if ( inputRecords <= 0 ) {
      fprintf( stderr, "\nInvalid rows in dbf file: %d.\n",
               inputRecords );
    } else {
      result = inputRecords;
    }

    DBFClose( inputFile ), inputFile = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: writeShapeData - Write point multi-data ShapeData to DBF and SHP.
INPUTS:  const char* fileName        Pathed base file name to create.
         const ShapeData* shapeData  Data to write.
RETURNS: int 1 if successful, else 0 failure message is printed.
NOTES:   Uses external shapelib API.
******************************************************************************/

int writeShapeData( const char* fileName, const ShapeData* shapeData ) {

  PRE03( fileName, *fileName, isValidShapeData( shapeData ) );

  int result = 0;
  const int rows = shapeData->rows;
  const int columns = shapeData->columns;
  const Value* const values = shapeData->values;
  DBFHandle outputFile = DBFCreate( fileName );
  DEBUG( fprintf( stderr, "writeShapeData( %s, %p )\n", fileName, shapeData );)

  if ( outputFile ) {
    int outputIndex = -1;
    int ok = defineDBFColumns(fileName, 1, &outputIndex, 0,0,0,0,0, outputFile);
    int row = 0;

    /* Write rows: */

    for ( row = 0; AND2( ok, row < rows ); ++row ) {
      const int rowOffset = row * columns;
      const Value* const rowValues = values + rowOffset;
      int column = 0;

      for ( column = 0; AND2( ok, column < columns ); ++column ) {
        const int type = shapeData->columnTypes[ column ];

        if ( type == FTDouble ) {
          ok = AND2( ok,
                     DBFWriteDoubleAttribute( outputFile, row, column,
                                              rowValues[ column ].d ) );
        } else if ( type == FTInteger ) {
          ok = AND2( ok,
                     DBFWriteIntegerAttribute( outputFile, row, column,
                                               rowValues[ column ].i ) );
        } else {
          CHECK( type == FTString );
          ok = AND2( ok,
                     DBFWriteStringAttribute( outputFile, row, column,
                                              rowValues[ column ].s ) );
        }
      }
    }

    result = ok;
    DBFClose( outputFile ), outputFile = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "Failed to write to dbf file.\n" );
  } else { /* Write longitude-latitude coordinates to shp file: */
    SHPHandle outputFile = SHPCreate( fileName, SHPT_POINT );
    int ok = outputFile != 0;

    if ( outputFile ) {
      const int longitudeColumn =
        indexOfString( "LONGITUDE", (const char**) shapeData->columnNames,
                       columns );
      const int latitudeColumn =
        indexOfString( "LATITUDE", (const char**) shapeData->columnNames,
                        columns );
      int row = 0;

      /* Write longitude latitude coordinates for each row: */

      for ( row = 0; AND2( ok, row < rows ); ++row ) {
        const int rowOffset = row * columns;
        const Value* const rowValues = values + rowOffset;
        double longitude = rowValues[ longitudeColumn ].d;
        double latitude  = rowValues[ latitudeColumn  ].d;
        SHPObject* object =
          SHPCreateSimpleObject( SHPT_POINT, 1, &longitude, &latitude, 0 );
        ok = object != 0;

        if ( ! object ) {
          fprintf( stderr, "Failed to create SHP object.\n" );
        } else {
          ok = SHPWriteObject( outputFile, -1, object ) >= 0;

          if ( ! ok ) {
            fprintf( stderr, "Failed to write SHP object.\n" );
          }

          SHPDestroyObject( object ), object = 0;
        }
      }

      SHPClose( outputFile ), outputFile = 0;
    }

    result = ok;

    if ( ! result ) {
      fprintf( stderr, "Failed to write to shp file.\n" );
    }
  }

  if ( result ) {
    result = writePRJFile( fileName, 0 );
  }

  DEBUG2( fprintf( stderr, "writeShapeData() returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeBoundsDBF - Write DBF file for bounds.
INPUTS:  const char* fileName                 File name of output dbf.
         const double areaInSquareKilometers  Area to write to DBF.
RETURNS: int 1 if successful, else 0 failure message is printed.
NOTES:   Uses external shapelib API.
******************************************************************************/

int writeBoundsDBF(const char* fileName, const double areaInSquareKilometers) {

  PRE03( fileName, *fileName, IN_RANGE( areaInSquareKilometers, 0.001, 1e+10));

  int result = 0;
  DBFHandle file = DBFCreate( fileName );

  if ( file ) {
    result = DBFAddField( file, "AREA_SQKM", FTDouble, 11, 3 ) != -1;

    if ( result ) {
      result = DBFWriteDoubleAttribute( file, 0, 0, areaInSquareKilometers );
    }

    DBFClose( file ), file = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "Failed to write to dbf file.\n" );
  }

  DEBUG2( fprintf( stderr, "writeBoundsDBF() returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: streamShapefiles - Stream shx, shp, dbf files to stdout preceeded
         by a one line ASCII header listing the (unpathed) base file name and
         sizes in bytes of each of the files.
INPUTS:  const char* baseFileName  Base (no ext) name of shx/shp/dbf files.
         const char* name          Name to use in header.
         const int dbfOnly         Only process DBF file (not shx/p)?
         const int csv             Stream csv file too?
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int streamShapefiles( const char* baseFileName, const char* name,
                      const int dbfOnly, const int csv ) {
  PRE06( baseFileName, *baseFileName, name, *name,
         IS_BOOL( dbfOnly ), IS_BOOL( csv ) );
  int result = 0;
  char shxFileName[ MAXIMUM_FILE_NAME_LENGTH + 1 ] = "";
  char shpFileName[ MAXIMUM_FILE_NAME_LENGTH + 1 ] = "";
  char dbfFileName[ MAXIMUM_FILE_NAME_LENGTH + 1 ] = "";
  char csvFileName[ MAXIMUM_FILE_NAME_LENGTH + 1 ] = "";
  size_t shxBytes = 0;
  size_t shpBytes = 0;
  size_t dbfBytes = 0;
  size_t csvBytes = 0;
  memset( shxFileName, 0, sizeof shxFileName );
  memset( shpFileName, 0, sizeof shpFileName );
  memset( dbfFileName, 0, sizeof dbfFileName );
  memset( csvFileName, 0, sizeof csvFileName );
  snprintf( shxFileName, sizeof shxFileName / sizeof *shxFileName,
            "%s.shx", baseFileName );
  snprintf( shpFileName, sizeof shpFileName / sizeof *shpFileName,
            "%s.shp", baseFileName );
  snprintf( dbfFileName, sizeof dbfFileName / sizeof *dbfFileName,
            "%s.dbf", baseFileName );
  snprintf( csvFileName, sizeof csvFileName / sizeof *csvFileName,
            "%s.csv", baseFileName );
  CHECK( shxFileName[ MAXIMUM_FILE_NAME_LENGTH - 1 ] == '\0' );
  CHECK( shpFileName[ MAXIMUM_FILE_NAME_LENGTH - 1 ] == '\0' );
  CHECK( dbfFileName[ MAXIMUM_FILE_NAME_LENGTH - 1 ] == '\0' );
  CHECK( csvFileName[ MAXIMUM_FILE_NAME_LENGTH - 1 ] == '\0' );

  if ( ! dbfOnly ) {
    shxBytes = fileSize( shxFileName );
    shpBytes = fileSize( shpFileName );
  }

  dbfBytes = fileSize( dbfFileName );

  if ( csv ) {
    csvBytes = fileSize( csvFileName );
  }

  if ( dbfOnly ) {

    if ( csv ) {
      printf( "%s 0 0 %lu %lu\n", name, dbfBytes, csvBytes );
      result = streamFile( dbfFileName );
      result = AND2( result, streamFile( csvFileName ) );
    } else {
      printf( "%s 0 0 %lu\n", name, dbfBytes );
      result = streamFile( dbfFileName );
    }
  } else {

    if ( csv ) {
      printf( "%s %lu %lu %lu %lu\n",
              name, shxBytes, shpBytes, dbfBytes, csvBytes );
      result = streamFile( shxFileName );
      result = AND2( result, streamFile( shpFileName ) );
      result = AND2( result, streamFile( dbfFileName ) );
      result = AND2( result, streamFile( csvFileName ) );
    } else {
      printf( "%s %lu %lu %lu\n", name, shxBytes, shpBytes, dbfBytes );
      result = streamFile( shxFileName );
      result = AND2( result, streamFile( shpFileName ) );
      result = AND2( result, streamFile( dbfFileName ) );
    }
  }

  DEBUG( fprintf( stderr, "streamFiles() returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: removeShapefiles - Remove set of temporary shx, shp, dbf files.
INPUTS:  const char* baseFileName  Base (no ext) file name of files to remove.
******************************************************************************/

void removeShapefiles( const char* baseFileName ) {
  PRE02( baseFileName, *baseFileName );
  const char* extensions[] = { "dbf", "shx", "shp", "prj", "xml", "csv" };
  const size_t count = sizeof extensions / sizeof *extensions;
  size_t index = 0;
  char outputFileName[ MAXIMUM_FILE_NAME_LENGTH + 1 ] = "";
  memset( outputFileName, 0, sizeof outputFileName );

  for ( index = 0; index < count; ++index ) {
    const char* const extension = extensions[ index ];
    snprintf( outputFileName, sizeof outputFileName / sizeof *outputFileName,
              "%s.%s", baseFileName, extension );
    CHECK( outputFileName[ MAXIMUM_FILE_NAME_LENGTH - 1 ] == '\0' );

    if ( fileExists( outputFileName ) ) {
      unlink( outputFileName );
    }
  }
}



/******************************************************************************
PURPOSE: deallocateShapeData - Deallocate storage of a ShapeData.
INPUTS:  ShapeData* shapeData  Shape data to deallocate.
OUTPUTS: ShapeData* shapeData  Deallocated shapeData.
******************************************************************************/

void deallocateShapeData( ShapeData* shapeData ) {

  if ( shapeData ) {

    if ( shapeData->stringStorage ) {
      const int capacity = shapeData->capacity;
      int index = 0;

      for ( index = 0; index < capacity; ++index ) {
        FREE( shapeData->stringStorage[ index ] );
      }
    }

    FREE( shapeData->stringStorage );
    FREE( shapeData->columnTypes );
    FREE( shapeData->values );
    memset( shapeData, 0, sizeof *shapeData );
    FREE( shapeData );
  }

  POST0( shapeData == 0 );
}



/******************************************************************************
PURPOSE: printShapeData - Print a ShapeData (for tracing/debugging).
INPUTS:  const ShapeData* shapeData  Shape data to print.
******************************************************************************/

void printShapeData( const ShapeData* shapeData ) {

  if ( shapeData ) {
    const int rows = shapeData->rows;
    const int columns = shapeData->columns;
    int row = 0;
    int column = 0;

    fprintf( stderr, "shapeData: rows = %d, columns = %d, capacity = %d\n",
             shapeData->rows, shapeData->columns, shapeData->capacity );
    CHECK3( shapeData->columns <= shapeData->capacity,
            shapeData->columnNames, shapeData->columnNames[ 0 ] );

    for ( column = 0; column < columns; ++column ) {
      fprintf( stderr, "%s\t", shapeData->columnNames[ column ] );
    }

    fprintf( stderr, "\n" );

    for ( column = 0; column < columns; ++column ) {
      const char* const typeNames[ 3 ] =
        { "FTString", "FTInteger", "FTDouble" };
      const int type = shapeData->columnTypes[ column ];
      CHECK( IN4( type, FTString,  FTInteger , FTDouble ) );
      fprintf( stderr, "%s\t", typeNames[ type ] );
    }

    fprintf( stderr, "\n" );

    for ( row = 0; row < rows; ++row ) {
      const int rowOffset = row * columns;

      for ( column = 0; column < columns; ++column ) {
        const int index = rowOffset + column;

        switch ( shapeData->columnTypes[ column ] ) {
        case FTString:
          fprintf( stderr, "%s\t", shapeData->values[ index ].s );
          break;
        case FTInteger:
          fprintf( stderr, "%d\t", shapeData->values[ index ].i );
          break;
        default:
          CHECK( shapeData->columnTypes[ column ] == FTDouble );
          fprintf( stderr, "%lg\t", shapeData->values[ index ].d );
          break;
        }
      }

      fprintf( stderr, "\n" );
    }
  }
}



/******************************************************************************
PURPOSE: isValidShapeData - Validate a ShapeData.
INPUTS:  const ShapeData* shapeData  Shape data to validate.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidShapeData( const ShapeData* shapeData ) {
  int result =
    AND9( shapeData,
          shapeData->rows > 0,
          shapeData->columns > 0,
          shapeData->rows * shapeData->columns > 0,
          shapeData->stringStorage,
          shapeData->columnNames == shapeData->stringStorage,
          shapeData->columnTypes,
          shapeData->values,
          shapeData->capacity >= shapeData->columns );

if ( ! result ) fprintf( stderr, "\n\n====== bad initial result.\n" );

  if ( result ) {
    const int rows = shapeData->rows;
    const int columns = shapeData->columns;
    int row = 0;
    int column = 0;

    for ( column = 0; AND2( result, column < columns ); ++column ) {
      const char* const columnName = shapeData->columnNames[ column ];
      result =
        AND2( IN4( shapeData->columnTypes[ column ],
                   FTString,  FTInteger , FTDouble ),
              isValidColumnName( columnName ) );
    }

    /* Check that column names are unique: */

    for ( column = 1; AND2( result, column < columns ); ++column ) {
      const char* const columnName = shapeData->columnNames[ column - 1 ];
      result =
        indexOfString( columnName, 
                       (const char* const *) shapeData->columnNames + column,
                       columns - column ) == -1;
    }

    /* Check data values: */

    for ( row = 0; AND2( result, row < rows ); ++row ) {
      const int rowOffset = row * columns;

      for ( column = 0; AND2( result, column < columns ); ++column ) {
        const int index = rowOffset + column;

        switch ( shapeData->columnTypes[ column ] ) {
        case FTString:
          result = shapeData->values[ index ].s != 0;
          break;
        case FTDouble:
          result = ! isNan( shapeData->values[ index ].d );
          break;
        default:
          break;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: isValidColumnName - Is DBF column name valid?
INPUTS:  const char* const columnName  Name to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidColumnName( const char* columnName ) {
  const char* c = columnName;
  int result = AND2( c, *c );

  while ( AND2( result, *c ) ) {
    result = AND2( isprint( *c ), ! isspace( *c ) );
    ++c;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidValue - Is Value valid - i.e., non-missing-valued?
INPUTS:  int type           Type portion to check.
         const char* units  Units of value, e.g., %.
         Value value        Value to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidValue( int type, const char* units, Value value ) {
  PRE0( IN4( type, FTString, FTInteger, FTDouble ) );
  const double validDoubleMinimum  = *units == '%' ? 0.0 : -98.0;
  const double validIntegerMinimum = *units == '%' ? 0   : -98;
  const int result =
    OR3( AND2( type == FTDouble,  value.d >= validDoubleMinimum ),
         AND2( type == FTInteger, value.i >= validIntegerMinimum ),
         AND3( type == FTString,  value.s, value.s[ 0 ] ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: compareValues - Compare values for sorting order.
INPUTS:  Value value1  1st value to check.
         Value value2  2nd value to check.
RETURNS: int -1 if value1 < value2 else 1 if value1 > value2 else 0.
******************************************************************************/

int compareValues( int type, Value value1, Value value2 ) {
  PRE0( IN4( type, FTString, FTInteger, FTDouble ) );
  int result = 0; /* null == null. */

  if ( type == FTDouble ) {
    result = value1.d < value2.d ? -1 : value1.d > value2.d ? 1 : 0;
  } else if ( type == FTInteger ) {
    result = value1.i < value2.i ? -1 : value1.i > value2.i ? 1 : 0;
  } else {
    CHECK( type == FTString );

    if ( value1.s ) {

      if ( value2.s ) {
        result = strcmp( value1.s, value2.s );

        if ( result < 0 ) {
          result = -1;
        } else if ( result > 0 ) {
          result = 1;
        }
      } else {
        result = 1; /* non-null > null. */
      }
    } else if ( value2.s ) {
      result = -1; /* null < non-null. */
    }
  }

  POST0( IN4( result, -1, 0, 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: readDBF - Read a DBF file.
INPUTS:  const char* fileName  Name of input dbf file to read.
RETURNS: ShapeData* if successful, else 0 and failure message is printed.
******************************************************************************/

ShapeData* readDBF( const char* fileName ) {
  PRE02( fileName, *fileName );
  ShapeData* result = 0;
  DBFHandle inputFile = DBFOpen( fileName, "rb" );
  int ok = 0;

  if ( inputFile ) {
    const int rows = DBFGetRecordCount( inputFile );

    if ( rows > 0 ) {
      const int columns = DBFGetFieldCount( inputFile );

      if ( columns > 0 ) {
        result = NEW( ShapeData, 1 );

        if ( result ) {
          result->rows = rows;
          result->columns = columns;
          result->columnTypes = NEW( int, columns );

          if ( result->columnTypes ) {
            const int estimatedUniqueStringCount = 1000;
            result->capacity = columns + estimatedUniqueStringCount;
            result->stringStorage = NEW( char*, result->capacity );
            result->columnNames = result->stringStorage;
            CHECK2( result->capacity >= columns,
                    result->columnNames == result->stringStorage );

            if ( result->stringStorage ) {
              result->values = NEW( Value, rows * columns );

              if ( result->values ) {
                int row = 0;
                int column = 0;
                ok = 1;

                /* Read header column names/types: */

                for ( column = 0; AND2( ok, column < columns ); ++column ) {
                  char columnName[ 16 ] = "";
                  const DBFFieldType type =
                    DBFGetFieldInfo( inputFile, column, columnName, 0, 0 );
                  ok = AND2( IN4( type, FTString, FTInteger, FTDouble ),
                             isValidColumnName( columnName ) );

                  if ( ok ) {
                    result->columnTypes[ column ] = type;
                    result->stringStorage[ column ] = strdup( columnName );
                    ok = result->stringStorage[ column ] != 0;

                    if ( ok ) {
                      DEBUG2( fprintf( stderr,
                              "result->stringStorage[ column = %d ] = '%s'\n",
                              column, result->stringStorage[ column ] ); )
                      CHECK( ! strcmp( result->stringStorage[ column ],
                                       columnName ) );
                    }
                  }
                }

                /* Read row data: */

                for ( row = 0; AND2( ok, row < rows ); ++row ) {
                  const int rowOffset = row * columns;

                  for ( column = 0; AND2( ok, column < columns ); ++column ) {
                    const int index = rowOffset + column;
                    ok = 0;

                    switch ( result->columnTypes[ column ] ) {
                    case FTString:
                      {
                        const char* stringValue =
                          DBFReadStringAttribute( inputFile, row, column );

                        if ( OR2( stringValue == 0, *stringValue == '\0' ) ) {
                          stringValue = "NULL";
                        }

                        result->values[ index ].s =
                          storeStringValue( stringValue, result );
                        ok = result->values[ index ].s != 0;
                      }
                      break;
                    case FTInteger:
                      result->values[ index ].i =
                        DBFReadIntegerAttribute( inputFile, row, column );
                      ok = 1;
                      break;
                    default:
                      CHECK( result->columnTypes[ column ] == FTDouble );
                      result->values[ index ].d =
                        DBFReadDoubleAttribute( inputFile, row, column );
                      ok = 1;
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    DBFClose( inputFile ), inputFile = 0;
  }

  if ( ! ok ) {
    fprintf( stderr, "Failed to read dbf file %s.\n", fileName );
    deallocateShapeData( result ), result = 0;
  }

  DEBUG2( printShapeData( result ); )
  POST0( IMPLIES( result, isValidShapeData( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: subsetDBFByTime - Compute mask of dbf rows with Start/End column
         values within a time range.
INPUTS:  const char* fileName  Base file name of DBF file to read.
         const int yyyymmdd1   1st timestamp to filter by.
         const int yyyymmdd2   2nd timestamp to filter by.
         const int count       Expected row count to match.
OUTPUTS: char mask[ count ]    1 = within time subset, else 0.
RETURNS: int number of rows within the time subset.
******************************************************************************/

int subsetDBFByTime( const char* fileName,
                     const int yyyymmdd1,
                     const int yyyymmdd2,
                     const int count,
                     char mask[] ) {

  PRE07( fileName,
         *fileName,
        isValidYearMonthDay( yyyymmdd1 ),
        isValidYearMonthDay( yyyymmdd2 ),
        yyyymmdd1 <= yyyymmdd2,
        count > 0,
        mask );

  int result = 0;
  ShapeData* shapeData = readDBF( fileName );
  memset( mask, 0, count * sizeof (char) );

  if ( shapeData ) {
    const int columns = shapeData->columns;
    const int startTimestampColumn =
      indexOfString( "Start",
                     (const char**) shapeData->columnNames, columns );
    const int endTimestampColumn =
      indexOfString( "End",
                     (const char**) shapeData->columnNames, columns );

    if ( AND4( startTimestampColumn != -1,
               endTimestampColumn != -1,
               shapeData->columnTypes[ startTimestampColumn ] == FTString,
               shapeData->columnTypes[ endTimestampColumn ] == FTString ) ) {
      const int yyyyddd1 = (int) ( convertYearMonthDay( yyyymmdd1 ) );
      const int yyyyddd2 = (int) ( convertYearMonthDay( yyyymmdd2 ) );
      const int rows = shapeData->rows;
      int row = 0;
      int ok = 1;
      CHECK2( isValidDate( yyyyddd1 ), isValidDate( yyyyddd2 ) );

      for ( row = 0; AND2( ok, row < rows ); ++row ) {
        const int rowOffset = row * columns;
        const int index1 = rowOffset + startTimestampColumn;
        const int index2 = rowOffset + endTimestampColumn;
        const char* const startTimestamp = shapeData->values[ index1 ].s;
        const char* const endTimestamp   = shapeData->values[ index2 ].s;
        int yyyydddStart = 0;
        int yyyydddEnd = 0;
        int hhmm = 0;
        ok = convertTimestamp(fileName, startTimestamp, &yyyydddStart, &hhmm);

        if ( ok ) {
          ok = convertTimestamp( fileName, endTimestamp, &yyyydddEnd, &hhmm );

          if ( ok ) {
            const int outside =
              OR2( yyyydddStart > yyyyddd2, yyyydddEnd < yyyyddd1 );

            if ( ! outside ) {
              mask[ row ] = 1;
              ++result;
            }
          }
        }
      }

      if ( ! ok ) {
        fprintf( stderr,
                 "\nInvalid date-time data read from DBF file %s.\n\n",
                 fileName );
        result = 0;
        memset( mask, 0, count * sizeof (char) );
      }
    }

    deallocateShapeData( shapeData ), shapeData = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: subsetDBFByTimeAndBoundsOrEstcode - Compute mask of dbf rows with
         YEAR, MONTH, DAY column values within a time range and either with
         matching estcode (if not 0) or else by STA_LON/STA_LAT columns
         values within bounds.
INPUTS:  const char* fileName  Base file name of DBF file to read.
         const int yyyymmdd1   1st timestamp to filter by.
         const int yyyymmdd2   2nd timestamp to filter by.
         const Bounds bounds   Lon-lat bounds to filter by if estcode unused.
         const char* estcode   ESTCODE_N to match unless 0.
         const int count       Expected row count to match.
OUTPUTS: char mask[ count ]    1 = within time subset, else 0.
RETURNS: int number of rows within the time subset.
******************************************************************************/

int subsetDBFByTimeAndBoundsOrEstcode( const char* fileName,
                                       const int yyyymmdd1,
                                       const int yyyymmdd2,
                                       const Bounds bounds,
                                       const char* estcode,
                                       const int count,
                                       char mask[] ) {

  PRE08( fileName,
         *fileName,
         isValidYearMonthDay( yyyymmdd1 ),
         isValidYearMonthDay( yyyymmdd2 ),
         yyyymmdd1 <= yyyymmdd2,
         isValidBounds( (const double (*)[2]) bounds ),
         count > 0,
         mask );

  int result = 0;
  ShapeData* shapeData = readDBF( fileName );
  memset( mask, 0, count * sizeof (char) );

  if ( shapeData ) {
    const int columns = shapeData->columns;
    const int estcodeColumn = estcode ?
      indexOfString( "ESTCODE_N",
                    (const char**) shapeData->columnNames, columns ) : -1;
    const int vstDateColumn =
      indexOfString( "VST_DATE",
                     (const char**) shapeData->columnNames, columns );
    const int dateColumn =
      vstDateColumn != -1 ? vstDateColumn :
        indexOfString( "DATE",
                      (const char**) shapeData->columnNames, columns );

    const int yearColumn =
      dateColumn == -1 ?
        indexOfString( "YEAR",
                       (const char**) shapeData->columnNames, columns ) : -1;
    const int monthColumn =
      dateColumn == -1 ?
        indexOfString( "MONTH",
                       (const char**) shapeData->columnNames, columns ) : -1;
    const int dayColumn =
      dateColumn == -1 ?
        indexOfString( "DAY",
                      (const char**) shapeData->columnNames, columns ) : -1;

    int longitudeColumn = -1;
    int latitudeColumn = -1;

    if ( estcodeColumn == -1 ) {
      longitudeColumn =
        indexOfString( "LONGITUDE",
                       (const char**) shapeData->columnNames, columns );

      if ( longitudeColumn == -1 ) {
        longitudeColumn =
          indexOfString( "STA_LNG",
                        (const char**) shapeData->columnNames, columns );

        if ( longitudeColumn == -1 ) {
          longitudeColumn =
            indexOfString( "LONGITUDE_",
                           (const char**) shapeData->columnNames, columns );
        }
      }

      latitudeColumn =
        indexOfString( "LATITUDE",
                       (const char**) shapeData->columnNames, columns );

      if ( latitudeColumn == -1 ) {
        latitudeColumn =
          indexOfString( "STA_LAT",
                        (const char**) shapeData->columnNames, columns );

        if ( latitudeColumn == -1 ) {
          latitudeColumn =
            indexOfString( "LATITUDE_D",
                           (const char**) shapeData->columnNames, columns );
        }
      }
    }

    DEBUG( fprintf( stderr, "estcodeColumn = %d, dateColumn = %d, "
                   "yearColumn = %d, monthColumn = %d, dayColumn = %d, "
                   "longitudeColumn = %d, latitudeColumn = %d\n",
                   estcodeColumn, dateColumn,
                   yearColumn, monthColumn, dayColumn,
                   longitudeColumn, latitudeColumn ); )

    if ( AND2( OR2( dateColumn != -1,
                    AND6( yearColumn  != -1,
                          monthColumn != -1,
                          dayColumn   != -1,
                          shapeData->columnTypes[ yearColumn  ] == FTInteger,
                          shapeData->columnTypes[ monthColumn ] == FTInteger,
                          shapeData->columnTypes[ dayColumn   ] == FTInteger)),
               IMPLIES_ELSE( estcode,
                             estcodeColumn != -1,
                             AND2( longitudeColumn != -1,
                                   latitudeColumn != -1 ) ) ) ) {
      const int dateColumnType =
        dateColumn == -1 ? -1 : shapeData->columnTypes[ dateColumn ];
      const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
      const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
      const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
      const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
      const int rows = shapeData->rows;
      const Value* rowValues = shapeData->values;
      int row = 0;

      for ( row = 0; row < rows; ++row, rowValues += columns ) {
        int yyyymmdd = 0;

        if ( dateColumn != -1 ) {

          if ( dateColumnType == FTInteger ) {
            yyyymmdd = rowValues[ dateColumn ].i;
          } else if ( dateColumnType == FTString ) {
            yyyymmdd = atoi( rowValues[ dateColumn ].s );
          }
       } else {
          const int year  = rowValues[ yearColumn  ].i;
          const int month = rowValues[ monthColumn ].i;
          const int day   = rowValues[ dayColumn   ].i;
          yyyymmdd = year * 10000 + month * 100 + day;
        }

        if ( isValidYearMonthDay( yyyymmdd ) ) {
          const int withinTime = IN_RANGE( yyyymmdd, yyyymmdd1, yyyymmdd2 );

          if ( withinTime ) {
            int inDomain = 0;

            if ( estcode ) {
              const char* const thisEstcode = rowValues[ estcodeColumn ].s;

              if ( AND2( thisEstcode, ! strcmp( thisEstcode, estcode ) ) ) {
                inDomain = 1;
              }
            } else {
              const double longitude = rowValues[ longitudeColumn ].d;
              const double latitude  = rowValues[ latitudeColumn  ].d;
              inDomain =
                AND2( IN_RANGE( longitude, longitudeMinimum, longitudeMaximum),
                      IN_RANGE( latitude,  latitudeMinimum,  latitudeMaximum));
            }

            if ( inDomain ) {
              mask[ row ] = 1;
              ++result;
            }
          }
        }
      }
    } else {
      fprintf( stderr,
               "\nInvalid date-time/location/ESTCODE_N data read from DBF "
               "file %s.\n\n", fileName );
      result = 0;
      memset( mask, 0, count * sizeof (char) );
    }

    deallocateShapeData( shapeData ), shapeData = 0;
  }

  return result;
}



/******************************************************************************
PURPOSE: writeShapeDataToText - Write ShapeData to a tab-delimited text file.
INPUTS:  const ShapeData* shapeData  ShapeData to write.
         const char* outputFileName  Name of output txt file to create.
RETURNS: int 1 if successful, else 0 and failure message is printed.
******************************************************************************/

int writeShapeDataToText( const ShapeData* shapeData,
                          const char* outputFileName ) {

  PRE03( isValidShapeData( shapeData ), outputFileName, *outputFileName );

  int result = 0;
  FILE* outputFile = fopen( outputFileName, "w" );

  if ( ! outputFile ) {
    fprintf( stderr, "Failed to create file %s.\n", outputFileName );
  } else {
    const int rows = shapeData->rows;
    const int columns = shapeData->columns;
    const int columns_1 = columns - 1;
    int column = 0;
    int row = 0;

    for ( column = 0; column < columns_1; ++column ) {
      fprintf( outputFile, "%s\t", shapeData->columnNames[ column ] );
    }

    fprintf( outputFile, "%s\n", shapeData->columnNames[ columns_1 ] );

    for ( row = 0; row < rows; ++row ) {

      for ( column = 0; column < columns; ++column ) {
        const int type = shapeData->columnTypes[ column ];
        const char delimiter = column < columns - 1 ? '\t' : '\n';
        const int index = row * columns + column;

        if ( type == FTString ) {
          fprintf( outputFile, "%s%c",
                   shapeData->values[ index ].s, delimiter );
        } else if ( type == FTInteger ) {
          fprintf( outputFile, "%d%c",
                   shapeData->values[ index ].i, delimiter );
        } else {
          const double value = shapeData->values[ index ].d;
          const double absValue = fabs( value );
          const int useExp =
            AND2( value != 0.0, OR2( absValue < 1e-3, absValue > 1e6 ) );

          if ( useExp ) {
            fprintf( outputFile, "%le%c", value, delimiter );
          } else {
            fprintf( outputFile, "%lf%c", value, delimiter );
          }
        }
      }
    }

    fclose( outputFile ), outputFile = 0;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: copyMaskedShapes - Copy masked shapes.
INPUTS:  const char* baseInputFileName    Base name of shp file to read.
         const char* baseOutputFileName   Base name of shp file to write.
         SHPHandle* shpFile               shpFile to open if not already open.
         const count                      Number of shapes.
         const char mask[ count ]         Mask: 1 = copy, 0 = omit.
OUTPUTS: SHPHandle* shpFile               Opened shpFile.
RETURNS: int subset shape count written if successful,
         else 0 and a failure message is printed to stderr.
******************************************************************************/

int copyMaskedShapes( const char* baseInputFileName,
                      const char* baseOutputFileName,
                      SHPHandle* shpFile,
                      const int count, const char mask[] ) {

  PRE07( baseInputFileName, *baseInputFileName,
         baseOutputFileName, *baseOutputFileName,
         shpFile, count > 0, mask );
  int result = 0;
  SHPHandle inputShpFile = SHPOpen( baseInputFileName, "rb" );

  if ( ! inputShpFile ) {
    fprintf( stderr, "Failed to open input shp file '%s'.\n",
             baseInputFileName );
  } else {
    int shapes = 0;
    int type = 0;
    double minimums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    double maximums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    SHPGetInfo( inputShpFile, &shapes, &type, minimums, maximums );

    if ( shapes != count ) {
      fprintf( stderr, "Mismatched shape count (%d) in shp file '%s'"
               " - expecting %d shapes to match dbf file rows.\n",
               shapes, baseInputFileName, count );
    } else {

      if ( ! *shpFile ) {
        *shpFile = SHPCreate( baseOutputFileName, type );

        if ( ! *shpFile ) {
          fprintf( stderr,
                   "Failed to create temporary shp file '%s'.\n",
                   baseOutputFileName );
        }
      }

      if ( *shpFile ) {
        const int subsetCount =
          copySubsetShapefile( inputShpFile, *shpFile, count, mask );
        DEBUG( fprintf( stderr, "%d = copySubsetShapefile()\n",
                        subsetCount ); )
        result = subsetCount;
        SHPClose( inputShpFile );
        inputShpFile = 0;
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: isPointType - Is shape file points?
INPUTS:  const char* baseInputFileName  Base name of shp file to read.
RETURNS: int 1 if point, else 0.
******************************************************************************/

int isPointType( const char* baseInputFileName ) {

  PRE02( baseInputFileName, *baseInputFileName );
  int result = 0;
  SHPHandle inputShpFile = SHPOpen( baseInputFileName, "rb" );

  if ( ! inputShpFile ) {
    fprintf( stderr, "Failed to open input shp file '%s'.\n",
             baseInputFileName );
  } else {
    int shapes = 0;
    int type = 0;
    double minimums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    double maximums[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
    SHPGetInfo( inputShpFile, &shapes, &type, minimums, maximums );

    result = IN4( type, SHPT_POINT, SHPT_POINTZ, SHPT_POINTM );
    SHPClose( inputShpFile );
    inputShpFile = 0;
  }

  return result;
}



/*============================= PRIVATE FUNCTIONS ===========================*/



/******************************************************************************
PURPOSE: writePRJFile - Write a .prj file.
INPUTS:  const char* fileName  Base name of file to create. "example".
         const int   useASCIIGridForm  1 = Write grid form, 0 = shape form.
NOTES:   Creates files fileName.prj.
         http://en.wikipedia.org/wiki/Well-known_text
******************************************************************************/

static int writePRJFile( const char* fileName, const int useASCIIGridForm ) {
  int result = 0;
  enum { FILE_NAME_LENGTH = 255 };
  char prjFileName[ FILE_NAME_LENGTH + 1 ] = "";
  memset( prjFileName, 0, sizeof prjFileName );
  strncpy( prjFileName, fileName, FILE_NAME_LENGTH - 4 );
  strcat( prjFileName, ".prj" );
  result = writeWGS84PRJFile( prjFileName, useASCIIGridForm );
  return result;
}



/******************************************************************************
PURPOSE: writeDataToDBFFile - Write location info to a dbf file and
         time-varying data to a csv file.
INPUTS:  const char* fileName  Base name of file to create. "IOOS".
         const char* variable  Name of data variable. E.g., "salinity".
         const char* units     Name of data units. E.g., "PSU".
         const int timesteps   Number of timesteps of data.
         const int yyyymmddhh[ timesteps ]  Timestamps of data.
         const int timestepType  HOURLY, DAILY, MONTHLY, YEARLY.
         const int count        Number of data rows, e.g., grid rows x columns.
         const int type         Data type: BYTE_TYPE, FLOAT_TYPE, etc.
         const int components   1 = scalar, 2 or 3 = vector.
         const void* data      data[ components * timesteps * count ]
                                           Component data at cells.
         const float lonlats[ count * 2 ]  Lon-lat of each row data point.
         const float z[ count ]            Depth(m) of each row data point or 0
         const char* const sids[ count ]   Station id strings or 0.
         const int ids[ count ]            Station ids or 0.
         const char* const metadata[ count ]   Station metadata or 0.
         const int writeCSV                Write CSV file too?
NOTES:   Creates files fileName.dbf containing only static location info:
         ID (int)    [Station (string)]  Longitude Latitude [Depth(m)]
         [] denotes optional (if available).
         ID is int in range [0, points - 1] index matching shapes in shp file.
         Station is string-id for the station, if station data, else omitted.
         Depth(m) is omitted if unavailable.
         Calls writeDataToCSVFile() to create the csf file.
         http://www.clicketyclick.dk/databases/xbase/format/dbf.html#DBF_STRUCT
         This routine does not use shapelib API since it was written before
         that library was used in this project and also because this routine
         is more straight-forward and much faster than using the shapelib API.
******************************************************************************/

static int writeDataToDBFFile( const char* fileName,
                               const char* variable, const char* units,
                               const int timesteps, const int yyyymmddhh[],
                               const int timestepType,
                               const int count, const int type,
                               const int components,
                               const void* data,
                               const float lonlats[], const float z[],
                               const char* const sids[],
                               const int ids[],
                               const char* const metadata[],
                               const int writeCSV ) {

  PRE017( fileName,
          *fileName,
          variable,
          *variable,
          units,
          *units,
          timesteps > 0,
          yyyymmddhh,
          IS_VALID_TIMESTEP_TYPE( timestepType ),
          count > 0,
          IN_RANGE( type, 0, GRID_DATA_TYPES ),
          data,
          IN_RANGE( components, 1, 3 ),
          lonlats,
          is_valid_longitude_latitude( lonlats[ 0 ], lonlats[ 1 ] ),
          is_valid_longitude_latitude( lonlats[  count * 2 - 2 ],
                                       lonlats[  count * 2 - 1 ] ),
          IMPLIES( z, AND2( ! is_nan( z[ 0 ] ), ! is_nan( z[ count - 1 ] ))));

  int result = 0;
  enum { /* Constants from the Shapefile Spec (see above link): */
    HEADER_BYTES = 32,
    FIELD_DESCRIPTOR_BYTES = 32,
    INT_FIELD_WIDTH = 20,     /* For ID. */
    FIELD_WIDTH = 15,         /* For Longitude, Latitude, Depth. */
    FIELD_DECIMAL_DIGITS = 6, /* Number of digits to the right of the decimal*/
    STRING_FIELD_WIDTH = 31,  /* For Station from sids[]. */
    FILE_NAME_LENGTH = 255
  };
  const char* const stringFieldFormat = "%31s";
  const char* const metadataColumns[] = { "ESTCODE", "STATE", "PROJECT" };
  const size_t metadataColumnCount =
    sizeof metadataColumns / sizeof *metadataColumns;
  const int metadataCount = metadata ? countChar( metadata[0], ';' ) + 1 : 0;
  const int metadataParts =
    metadataCount > metadataColumnCount ? metadataColumnCount : metadataCount;
  const int stringCount = ( sids != 0 ) + metadataParts;
  const int fields = 1 + stringCount + 1 + 1 + ( z != 0 );
  /* ID (Station,ESTCODE,STATE,PROJECT) Longitude Latitude Depth. */
  const int records = count;
  const int headerBytes = HEADER_BYTES + fields * FIELD_DESCRIPTOR_BYTES + 1;
  unsigned char* header = NEW( unsigned char, headerBytes + 1 );

  DEBUG( fprintf( stderr, "writeDataToDBFFile( fileName = '%s', "
                  "variable = '%s', units = '%s', "
                  "timesteps = %d, yyyymmddhh[] = [%d ... %d], "
                  "timestepType = %d, count = %d, type = %d, "
                  "components = %d, data = %p, "
                  "lonlats = %p, z = %p, sids = %p, ids = %p, metadata = %p )\n",
                  fileName, variable, units,
                  timesteps, yyyymmddhh[ 0 ], yyyymmddhh[ timesteps - 1 ],
                  timestepType, count, type, components, data,
                  lonlats, z, sids, ids, metadata ); )

  if ( header ) {
    int metadataPart = 0;
    int headerIndex = 32;
    FILE* file = 0;
    char dbfFileName[ FILE_NAME_LENGTH + 1 ] = "";
    memset( dbfFileName, 0, sizeof dbfFileName );
    strncpy( dbfFileName, fileName, FILE_NAME_LENGTH - 4 );
    strcat( dbfFileName, ".dbf" );

    /* Initialize header: */

    memset( header, 0, headerBytes );

    header[ 0 ] = 3;  /* dBASE Level 5. */
    header[ 1 ] = 95; /* Update timestamp: YY - 1900. */
    header[ 2 ] = 7;  /* Update timestamp: MM. */
    header[ 3 ] = 26; /* Update timestamp: DD. */
    writeInt( header, 4, count, LITTLE ); /* Number of records. */
    writeShort( header, 8, headerBytes, LITTLE ); /* Number of header bytes. */

    /* Length of each record: */

    {
      const int recordLength =
        INT_FIELD_WIDTH                    /* ID. */
        + stringCount * STRING_FIELD_WIDTH /* (Station,ESTCODE,STATE,PROJECT)*/
        + 2 * FIELD_WIDTH                  /* Longitude, Latitude. */
        + ( z != 0 ) * FIELD_WIDTH;        /* Depth. */
       writeShort( header, 10, recordLength, LITTLE );      
    }

    /* Array of Field Descriptors: */

    snprintf( (char*) header + headerIndex, 11, "%-10s", "ID" );
    header[ headerIndex + 11 ] = 'N';  /* Formatted number. */
    header[ headerIndex + 16 ] = INT_FIELD_WIDTH;
    header[ headerIndex + 17 ] = 0; /* No digits to the right of decimal. */
    headerIndex += FIELD_DESCRIPTOR_BYTES;

    if ( sids ) {
      snprintf( (char*) header + headerIndex, 11, "%-10s", "Station" );
      header[ headerIndex + 11 ] = 'C';  /* Character string. */
      header[ headerIndex + 16 ] = STRING_FIELD_WIDTH;
      header[ headerIndex + 17 ] = 0;
      headerIndex += FIELD_DESCRIPTOR_BYTES;
    }

    for ( metadataPart = 0; metadataPart < metadataParts; ++metadataPart ) {
      CHECK( metadataPart < sizeof metadataColumns / sizeof *metadataColumns );
      snprintf( (char*) header + headerIndex, 11, "%-10s",
                metadataColumns[ metadataPart ] );
      header[ headerIndex + 11 ] = 'C';  /* Character string. */
      header[ headerIndex + 16 ] = STRING_FIELD_WIDTH;
      header[ headerIndex + 17 ] = 0;
      headerIndex += FIELD_DESCRIPTOR_BYTES;
    }

    snprintf( (char*) header + headerIndex, 11, "%-10s", "Longitude" );
    header[ headerIndex + 11 ] = 'N';  /* Formatted real number. */
    header[ headerIndex + 16 ] = FIELD_WIDTH;
    header[ headerIndex + 17 ] = FIELD_DECIMAL_DIGITS;
    headerIndex += FIELD_DESCRIPTOR_BYTES;

    snprintf( (char*) header + headerIndex, 11, "%-10s", "Latitude" );
    header[ headerIndex + 11 ] = 'N';  /* Formatted real number. */
    header[ headerIndex + 16 ] = FIELD_WIDTH;
    header[ headerIndex + 17 ] = FIELD_DECIMAL_DIGITS;
    headerIndex += FIELD_DESCRIPTOR_BYTES;

    if ( z ) {
      snprintf( (char*) header + headerIndex, 11, "%-10s", "Depth(m)" );
      header[ headerIndex + 11 ] = 'N';  /* Formatted real number. */
      header[ headerIndex + 16 ] = FIELD_WIDTH;
      header[ headerIndex + 17 ] = FIELD_DECIMAL_DIGITS;
      headerIndex += FIELD_DESCRIPTOR_BYTES;
    }

    header[ headerIndex ] = 0x0D; /* Terminator. */

    file = fopen( dbfFileName, "wb" );

    if ( file ) {
      int record = 0;

      /* Write header bytes to DBF file: */

      result = fwrite( header, headerBytes, 1, file ) == 1;

      /* Write data rows (records = shapes in shp file) to DBF file: */

      for ( record = 0; AND2( result, record < records ); ++record ) {
        const int record2 = record + record;
        const int id = ids ? ids[ record ] : record;
        const char* const station = sids ? sids[ record ] : 0;
        const char* theMetadata = metadata ? metadata[ record ] : 0;
        const float longitude = lonlats[ record2 ];
        const float latitude  = lonlats[ record2 + 1 ];
        const float depth = z ? z[ record ] : 0.0;
        char format[ 32 ] = "";
        char formattedIntValue[ INT_FIELD_WIDTH + 1 ] = "";
        char formattedStringValue[ STRING_FIELD_WIDTH + 1 ] = "";
        char formattedValue[ FIELD_WIDTH + 1 ] = "";

        /* Write Id: */

        memset( format, 0, sizeof format );
        memset( formattedIntValue, 0, sizeof formattedIntValue );
        snprintf( format, sizeof format / sizeof *format,
                  "%%%dd", INT_FIELD_WIDTH );
        snprintf( formattedIntValue,
                  sizeof formattedIntValue / sizeof *formattedIntValue,
                  format, id );
        result = fwrite( formattedIntValue, INT_FIELD_WIDTH, 1, file ) == 1;

        if ( result ) {

          /* Write Station: */

          if ( station ) {
            memset( formattedStringValue, 0, sizeof formattedStringValue );
            snprintf( formattedStringValue,
                      sizeof  formattedStringValue /
                      sizeof *formattedStringValue,
                      stringFieldFormat, station );
            result =
              fwrite( formattedStringValue, STRING_FIELD_WIDTH, 1, file ) == 1;
          }

          for (metadataPart = 0; metadataPart < metadataParts; ++metadataPart) {
            char metadataPartString[ STRING_FIELD_WIDTH + 1 ] = "";

            if ( OR2( theMetadata == 0, *theMetadata == '\0' ) ) {
              theMetadata = "NULL";
            }

            memset( metadataPartString, 0, sizeof metadataPartString );
            strncpy( metadataPartString, theMetadata,
                     sizeof metadataPartString / sizeof *metadataPartString -1);
            eraseChar( metadataPartString, ';' );

            memset( formattedStringValue, 0, sizeof formattedStringValue );
            snprintf( formattedStringValue,
                      sizeof  formattedStringValue /
                      sizeof *formattedStringValue,
                      stringFieldFormat, metadataPartString );
            result =
              fwrite( formattedStringValue, STRING_FIELD_WIDTH, 1, file ) == 1;

            theMetadata = strchr( theMetadata + 1, ';' );

            if ( theMetadata ) {
              ++theMetadata;
            }
          }

          /* Write Longitude: */

          if ( result ) {
            memset( format, 0, sizeof format );
            snprintf( format, sizeof format / sizeof *format,
                      "%%%d.%df", FIELD_WIDTH, FIELD_DECIMAL_DIGITS );
            memset( formattedValue, 0, sizeof formattedValue );
            snprintf( formattedValue,
                      sizeof formattedValue / sizeof *formattedValue,
                      format, longitude );
            result = fwrite( formattedValue, FIELD_WIDTH, 1, file ) == 1;

            /* Write Latitude: */

            if ( result ) {
              memset( formattedValue, 0, sizeof formattedValue );
              snprintf( formattedValue,
                        sizeof formattedValue / sizeof *formattedValue,
                        format, latitude );
              result = fwrite( formattedValue, FIELD_WIDTH, 1, file ) == 1;

              /* Write Depth: */

              if (result ) {

                if ( z ) {
                  memset( formattedValue, 0, sizeof formattedValue );
                  snprintf( formattedValue,
                            sizeof formattedValue / sizeof *formattedValue,
                            format, depth );
                  result = fwrite( formattedValue, FIELD_WIDTH, 1, file ) == 1;
                }
              }
            }
          }
        }
      }

      fclose( file ), file = 0;
    }

    if ( ! result ) {
      perror( "\n\nFailed to write DBF file because" );
    }

    memset( header, 0, headerBytes );
    FREE( header );
  }

  if ( AND2( result, writeCSV ) ) {
    result = writeDataToCSVFile( fileName, variable, units,
                                 timesteps, yyyymmddhh, timestepType,
                                 count, type, components, data, lonlats,
                                 z, sids, ids );
  }

  return result;
}



/******************************************************************************
PURPOSE: writeDataToCSVFile - Write location info and
         time-varying data to a csv file.
INPUTS:  const char* fileName  Base name of file to create. "IOOS".
         const char* variable  Name of data variable. E.g., "salinity".
         const char* units     Name of data units. E.g., "PSU".
         const int timesteps   Number of timesteps of data.
         const int yyyymmddhh[ timesteps ]  Timestamps of data.
         const int timestepType  HOURLY, DAILY, MONTHLY, YEARLY.
         const int count        Number of data rows, e.g., grid rows x columns.
         const int type         Data type: BYTE_TYPE, FLOAT_TYPE, etc.
         const int components   1 = scalar, 2 or 3 = vector.
         const void* data      data[ components * timesteps * count ]
                                           Component data at cells.
         const float lonlats[ count * 2 ]  Lon-lat of each row data point.
         const float z[ count ]            Depth(m) of each row data point or 0
         const char* const sids[ count ]   Station id strings or 0.
         const int ids[ count ]            Station ids or 0.
NOTES:   Creates file fileName.csv containing static location info:
         ID (int)    [Station (string)]  Longitude Latitude [Depth(m)]
         [] denotes optional (if available).
         ID is int in range [0, points - 1] index matching shapes in shp file.
         Station is string-id for the station, if station data, else omitted.
         Depth(m) is omitted if unavailable.
         And time-varying data:
         YYYYMMDDHH   Quantity    Units    Value
         or YYYYMMDD
         or YYYYMM
         or YYYY
         With rows ordered by: ID then timestamp then variable/component.
******************************************************************************/

static int writeDataToCSVFile( const char* fileName,
                               const char* variable, const char* units,
                               const int timesteps, const int yyyymmddhh[],
                               const int timestepType,
                               const int count, const int type,
                               const int components,
                               const void* data,
                               const float lonlats[], const float z[],
                               const char* const sids[],
                               const int ids[] ) {

  PRE017( fileName,
          *fileName,
          variable,
          *variable,
          units,
          *units,
          timesteps > 0,
          yyyymmddhh,
          IS_VALID_TIMESTEP_TYPE( timestepType ),
          count > 0,
          IN_RANGE( type, 0, GRID_DATA_TYPES ),
          data,
          IN_RANGE( components, 1, 3 ),
          lonlats,
          is_valid_longitude_latitude( lonlats[ 0 ], lonlats[ 1 ] ),
          is_valid_longitude_latitude( lonlats[  count * 2 - 2 ],
                                       lonlats[  count * 2 - 1 ] ),
          IMPLIES( z, AND2( ! is_nan( z[ 0 ] ), ! is_nan( z[ count - 1 ] ))));

  int result = 0;
  enum { MAXIMUM_LINE_LENGTH = 512 }; /* Of formatted data row. */
  size_t bufferSize = timesteps * MAXIMUM_LINE_LENGTH;
  char* buffer = NEW( char, bufferSize ); /* Buffer ASCII output for speed. */

  if ( buffer ) {
    enum { FILE_NAME_LENGTH = 255 };
    FILE* csvFile = 0;
    char csvFileName[ FILE_NAME_LENGTH + 1 ] = "";

    DEBUG( fprintf( stderr, "writeDataToCSVFile( fileName = '%s', "
                    "variable = '%s', units = '%s', "
                    "timesteps = %d, yyyymmddhh[] = [%d ... %d], "
                    "timestepType = %d, count = %d, type = %d, "
                    "components = %d, data = %p, "
                    "lonlats = %p, z = %p, sids = %p )\n",
                    fileName, variable, units,
                    timesteps, yyyymmddhh[ 0 ], yyyymmddhh[ timesteps - 1 ],
                    timestepType, count, type, components, data,
                    lonlats, z, sids ); )

    memset( csvFileName, 0, sizeof csvFileName );
    strncpy( csvFileName, fileName, FILE_NAME_LENGTH - 4 );
    strcat( csvFileName, ".csv" );
    csvFile = fopen( csvFileName, "w" );
    result = csvFile != 0;

    if ( csvFile == 0 ) {
      perror( "\n\nFailed because" );
    } else {
      const int componentOffset = timesteps * count;
      /* ID Station Longitude Latitude Depth YYYYMMDDHH Quantity Units. */
      const float* const fdata = type == FLOAT_TYPE ? data : 0;
      const char* const  cdata = type == BYTE_TYPE  ? data : 0;
      const unsigned short* const sdata = type == UINT16_TYPE ? data : 0;
      const float missingValue = cdata ? -99.0 : -9999.0;
      int point = 0;
      const char* const timestampColumnName =
        timestepType == YEARLY ? "YYYY"
        : timestepType == MONTHLY ? "YYYYMM"
        : timestepType == DAILY ? "YYYYMMDD"
        : "YYYYMMDDHH";
      const char* const valueColumnNames =
        components == 3 ? "U_Value,V_Value,W_Value"
        : components == 2 ? "U_Value,V_Value"
        : "Value";
      const char* const units0 = ! strcmp( units, "%s" ) ? "-" : units;
      /* If units = '%s' then print it as '-' instead. */

      /* Write one-line header column names: */

      fprintf( csvFile, "ID" );

      if ( sids ) {
        fprintf( csvFile, ",Station" );
      }

      fprintf( csvFile, ",Longitude,Latitude" );

      if ( z ) {
        fprintf( csvFile, ",Depth(m)" );
      }

      fprintf( csvFile, ",%s", timestampColumnName );
      fprintf( csvFile, ",Quantity,Units,%s\n", valueColumnNames );

      /* Write data rows: */

      for ( point = 0; AND2( result, point < count ); ++point ) {
        const int point2 = point + point;
        const float longitude = lonlats[ point2 ];
        const float latitude  = lonlats[ point2 + 1 ];
        const int id = ids ? ids[ point ] : point;
        const char* const station = sids ? sids[ point ] : 0;
        const float depth = z ? z[ point ] : 0.0;
        int timestep = 0;
        int bytes = 0;
        size_t bufferBytesStored = 0;
        char* outputBuffer = buffer;

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const int yyyymmddhh0 = yyyymmddhh[ timestep ];
          const int timestamp =
            timestepType == YEARLY ? yyyymmddhh0 / 1000000
            : timestepType == MONTHLY ? yyyymmddhh0 / 10000
            : timestepType == DAILY ? yyyymmddhh0 / 100
            : yyyymmddhh0;
          const int timesteppedPointOffset = timestep * count + point;
          int component = 0;
          int rowBytes = 0;
          int eraseRow = 0;

          /* Write row line: */

          bytes = snprintf( outputBuffer, MAXIMUM_LINE_LENGTH, "%d", id );
          outputBuffer += bytes;
          bufferBytesStored += bytes;
          rowBytes += bytes;

          if ( station ) {
            bytes =
              snprintf( outputBuffer, MAXIMUM_LINE_LENGTH, ",%s", station );
            outputBuffer += bytes;
            bufferBytesStored += bytes;
            rowBytes += bytes;
          }

          bytes = snprintf( outputBuffer, MAXIMUM_LINE_LENGTH,
                            ",%f,%f", longitude, latitude );
          outputBuffer += bytes;
          bufferBytesStored += bytes;
          rowBytes += bytes;

          if ( z ) {
            bytes = snprintf( outputBuffer, MAXIMUM_LINE_LENGTH, ",%f", depth );
            outputBuffer += bytes;
            bufferBytesStored += bytes;
            rowBytes += bytes;
          }

          bytes =
            snprintf( outputBuffer, MAXIMUM_LINE_LENGTH, ",%d,%s,%s",
                      timestamp, variable, units0 );
          outputBuffer += bytes;
          bufferBytesStored += bytes;
          rowBytes += bytes;

          for ( component = 0; component < components; ++component ) {
            const int index =
              component * componentOffset + timesteppedPointOffset;
            const float value =
              fdata ? fdata[ index ]
              : cdata ? (float) (cdata[ index ])
              : sdata ? (float) (sdata[ index ])
              : 0.0;

            if ( OR2( includeMissingValuesInCSVFile, value > missingValue ) ) {
              const float absValue = fabs( value );
              const int useExp =
                AND2( value != 0.0, OR2( absValue < 1e-3, absValue > 1e6 ) );
              const char* const format = useExp ? ",%e" : ",%f";
              CHECK( IN_RANGE( index, 0, components * timesteps * count - 1 ));

              DEBUG( fprintf( stderr,
                        "point = %d, timestep = %d, data[ index = %d ] = %f\n",
                        point, timestep, index, value ); )

              bytes =
                snprintf( outputBuffer, MAXIMUM_LINE_LENGTH, format, value );
              outputBuffer += bytes;
              bufferBytesStored += bytes;
              rowBytes += bytes;
            } else {
              eraseRow = 1;
            }
          } /* End loop on components. */

          bytes = snprintf( outputBuffer, MAXIMUM_LINE_LENGTH, "\n" ); /* End*/
          outputBuffer += bytes;
          bufferBytesStored += bytes;
          rowBytes += bytes;

          if ( eraseRow ) {
            outputBuffer -= rowBytes;
            bufferBytesStored -= rowBytes;
          }

          CHECK( bufferBytesStored / ( timestep + 1 ) < MAXIMUM_LINE_LENGTH );
        } /* End loop on timesteps. */

        /* Write buffer to CSV file: */

        result =
          fwrite( buffer, 1, bufferBytesStored, csvFile ) == bufferBytesStored;
        memset( buffer, 0, bufferSize );
      } /* End loop on points. */

      fclose( csvFile ), csvFile = 0;
    }

    FREE( buffer );
  }

  if ( ! result ) {
    perror( "\n\nFailed to write CSV file because" );
  }

  return result;
}



/******************************************************************************
PURPOSE: appendCSVFile - Read a time-matched subset of data from a time-sorted
         csv file and append it to an output csv file.
INPUTS:  const char* inputCSVDirectory  Name of directory containing input
                                        csv file for code.
         const char* code               Name of input csv file - code.csv.
         const int yyyymmdd1            Starting timestamp to subset by.
         const int yyyymmdd2            Ending   timestamp to subset by.
         const int timestepSize         YEARLY, MONTHLY, DAILY.
         const int columns              Number of header columns to output.
                                        0 means all.
OUTPUT:  CSVHeader header               Header to match or
                                        if "" to initialize to first columns.
         FILE* outputCSVFile            File to write subset of data lines to.
RETURNS: int number of time-subsetted data lines written.
******************************************************************************/

static int appendCSVFile( const char* inputCSVDirectory, const char* code,
                          const int yyyymmdd1, const int yyyymmdd2,
                          const int timestepSize,
                          const int columns, CSVHeader header,
                          FILE* outputCSVFile ) {

  PRE011( inputCSVDirectory, *inputCSVDirectory, code, *code,
          isValidYearMonthDay( yyyymmdd1 ),
          isValidYearMonthDay( yyyymmdd2 ),
          yyyymmdd1 <= yyyymmdd2,
          IN4( timestepSize, YEARLY, MONTHLY, DAILY ),
          columns >= 0,
          header,
          outputCSVFile );

  enum { FILE_NAME_LENGTH = 255 };
  char inputCSVFileName[ FILE_NAME_LENGTH + 1 ] = "";
  int result = 0;
  int ok = 0;
  char* buffer = 0;
  size_t length = 0;
  size_t lines = 0;
  memset( inputCSVFileName, 0, sizeof inputCSVFileName );
  snprintf( inputCSVFileName,
            sizeof inputCSVFileName / sizeof *inputCSVFileName,
            "%s/%s.csv", inputCSVDirectory, code );

  if ( fileExists( inputCSVFileName ) ) { /* Check 1st for estuary_flushing. */
    buffer = readFile( inputCSVFileName, &length, &lines );
    ok = buffer != 0;

    if ( ok ) {
      const int outputHeader = *header == '\0';
      char* line = parseCSVHeader( buffer, columns, header );
      int yyyymmdd = 0;
      ok = line != 0;

      if ( AND2( ok, outputHeader ) ) {
        ok = fputs( header, outputCSVFile ) >= 0;

        if ( ok ) {
          ok = fputc( '\n', outputCSVFile ) == '\n';
        }
      }

      /* Check timestamp of each line until greater than yyyymmdd2: */

      while ( AND4( ok, line, ok = parseCSVTimestamp( line, &yyyymmdd ),
                    compareDate( yyyymmdd, yyyymmdd2, timestepSize ) < 1 ) ) {
        char* end = strchr( line, '\n' );
        DEBUG( fprintf( stderr, "  ok = %d\n", ok ); )

        if ( compareDate( yyyymmdd, yyyymmdd1, timestepSize ) >= 0 ) {
          /* Output time-subset-matched line: */
          ok = writeCSVLine( line, columns, outputCSVFile );
          result += ok;
        }

        /* Advance to next line if there is one: */

        if ( end ) {
          line = end + 1;
        } else {
          line = 0;
        }
      }

      FREE( buffer );
    }
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: writeCSVLine - Write 1st few columns or all of data line to CSV file.
INPUTS:  char* line          Line to write.
         const int columns   Number of columns to match or 0 for all.
OUTPUTS: FILE outputCSVFile  File to append line to.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int writeCSVLine( char* line, const int columns, FILE* outputCSVFile ) {
  PRE03( line, columns >= 0, outputCSVFile );
  char* end = 0;
  int column = 0;
  int result = 0;

  if ( columns == 0 ) { /* Output full line: */
    end = strchr( line, '\n' );
  } else { /* Output only first columns of line: */
    char* c = line;

    while ( AND3( *c, *c != '\n', column < columns ) ) {

      if ( *c == ',' ) {
        ++column;
        end = c;
      }

      ++c;
    }

    if ( column < columns ) {
      end = c;
      ++column;
    }
  }

  if ( IMPLIES( columns, column == columns ) ) {
    char c0 = 0;

    if ( end ) {
      c0 = *end; /* Save delimiter. */
      *end = '\0'; /* Terminate full or partial line for output: */
    }

    result = fputs( line, outputCSVFile ) >= 0;

    if ( result ) {
      result = fputc( '\n', outputCSVFile ) == '\n';
    }

    if ( end ) {
      *end = c0; /* Restore delimiter. */
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: parseCSVHeader - Parse CSV file header line and match or initialize
         to 1st few or all columns.
INPUTS:  char* line         Line to parse.
         const int columns  Number of columns to match or 0 for all.
         CSVHeader header   Header to match 1st columns or if "" initialize.
OUTPUTS: CSVHeader header   Initialized header if it was "".
RETURNS: char* Next line if successful and another line follows, else 0.
******************************************************************************/

static char* parseCSVHeader( char* line, const int columns, CSVHeader header) {
  PRE03( line, columns >= 0, header );
  char* result = 0;
  char* c = line;
  int column = 0;
  int ok = 0;

  while ( AND2( *c != '\0', *c != '\n' ) ) {

    if ( *c == ',' ) {
      ++column;

      if ( column == columns ) { /* Check for partial match or initialize. */
        *c = '\0'; /* 0-terminate for copy/compare. */

        if ( *header == '\0' ) { /* Initialize: */
          strncpy( header, line, sizeof (CSVHeader) / sizeof (char) - 1 );
          header[ sizeof (CSVHeader) / sizeof (char) - 1 ] = '\0';
          ok = 1;
        } else { /* Compare: */
          ok = ! strcmp( line, header );
        }

        *c = ','; /* Restore. */
      }
    }

    ++c;
  }

  ++column;

  if ( OR2( columns == 0, column == columns ) ) { /* Init or match: */
    const char c0 = *c;
    *c = '\0'; /* 0-terminate for copy/compare. */

    if ( *header == '\0' ) { /* Initialize: */
      strncpy( header, line, sizeof (CSVHeader) / sizeof (char) - 1 );
      header[ sizeof (CSVHeader) / sizeof (char) - 1 ] = '\0';
      ok = 1;
    } else { /* Compare: */
      ok = ! strcmp( line, header );
    }

    *c = c0; /* Restore. */
  }

  if ( *c == '\n' ) {
    ++c;
  }

  ok = AND2( ok, *c );

  if ( ok ) {
    result = c;
  }

  POST0( IMPLIES( result, *result ) );
  return result;
}



/******************************************************************************
PURPOSE: parseCSVTimestamp - Parse timestamp of form yyyy, yyyy-mm, or
         yyyy-mm-dd possibly wrapped in double quotes.
INPUTS:  char* line     Address of string to parse.
         int* yyyymmdd  Parsed timestamp.
RETURNS: int 1 if successfully parsed a valid yyyymmdd else 0.
******************************************************************************/

static int parseCSVTimestamp( char* line, int* yyyymmdd ) {
  PRE02( line, yyyymmdd );
  int result = 0;
  char* c = line;

  *yyyymmdd = 0;

  /* Skip to after comma: */

  while ( AND2( *c != '\0', *c != ',' ) ) {
    ++c;
  }

  if ( *c == ',' ) {
    ++c;

    {
      char* end = 0;
      const int yyyy = (int) strtol( c, &end, 10 );
      int ok = IN_RANGE( yyyy, 1800, 2999 );

      if ( ok ) {
        int mm = 1;
        int dd = 1;

        if ( *end == '-' ) { /* Parse mm: */
          c = end + 1;
          mm = (int) strtol( c, &end, 10 );
          ok = IN_RANGE( mm, 1, 12 );

          if ( ok ) {

            if ( *end == '-' ) { /* Parse dd: */
              c = end + 1;
              dd = (int) strtol( c, &end, 10 );
              ok = IN_RANGE( dd, 1, 31 );
            }
          }
        }

        if ( ok ) {
          ok = IN3( *end, ',', 'T' );

          if ( ok ) {
            *yyyymmdd = yyyy * 10000 + mm * 100 + dd;
            ok = isValidYearMonthDay( *yyyymmdd );

            if ( ok ) {
              result = 1;
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    *yyyymmdd = 0;
  }

  DEBUG( fprintf(stderr, "%d = parseCSVTimestamp( %d )\n", result,*yyyymmdd );)
  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        isValidYearMonthDay( *yyyymmdd ), *yyyymmdd == 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: hasIntegerColumn - Does shapeData have named integer column?
INPUTS:  const ShapeData* shapeData ShapeData to check.
         const char* const name     Column name to check for.
RETURNS: int 1 if true, else 0.
******************************************************************************/

static int hasIntegerColumn( const ShapeData* const shapeData,
                             const char* const name ) {

  PRE06( shapeData, shapeData->columns > 0, shapeData->columnNames,
         shapeData->columnTypes, name, *name );

  const int column =
    indexOfString( name, (const char * const *) shapeData->columnNames,
                   shapeData->columns );
  const int result =
    AND2( column != -1, shapeData->columnTypes[ column ] == FTInteger );
  POST0( IS_BOOL( result ) );
  return result;
}


#if 0


/******************************************************************************
PURPOSE: flagUpstreamNodes - Mark a mask array indicating which
         rows are upstream of toNode.
INPUTS:  const int rows                       Number of rows of values.
         const int columns                    Number of columns per row.
         const Value values[ rows * columns ] Array containing from/toNode.
         const int fromNodeColumn             Column of fromNode.
         const int toNodeColumn               Column of toNode.
         const int toNode                     toNode value to match.
         char mask[ rows ]                    Set to 1 if upstream of toNode.
OUTPUTS: char mask[ rows ]                    Set to 1 if upstream of toNode.
RETURNS: int number of upstream rows.
NOTES:   This routine is cycle-safe depth-first non-tail recursive.
******************************************************************************/

static int flagUpstreamNodes( const int rows,
                              const int columns,
                              const Value values[],
                              const int fromNodeColumn,
                              const int toNodeColumn,
                              const int toNode,
                              char mask[] ) {

  PRE08( rows > 0,
         columns > 0,
         values,
         IN_RANGE( fromNodeColumn, 0, columns - 1 ),
         IN_RANGE( toNodeColumn,   0, columns - 1 ),
         fromNodeColumn != toNodeColumn,
         toNode > 0.0,
         mask );

  int result = 0;
  int row = 0;
  int rowOffset = 0;

  DEBUG( fprintf( stderr, "flagUpstreamNodes( toNode = %d )\n", toNode ); )

  /* Check each row's toNode for a match to given toNode: */

  for ( row = 0; row < rows; ++row, rowOffset += columns ) {

    /* Must compare to 0 to avoid infinite loop in case graph has cycles! */

    if ( mask[ row ] == 0 ) { /* If row is not already marked: */
      const int thisRowToNode = values[ rowOffset + toNodeColumn ].i;

      if ( thisRowToNode == toNode ) {
        const int thisRowFromNode = values[ rowOffset + fromNodeColumn ].i;
        mask[ row ] = 1;
        DEBUG( fprintf( stderr, " found row %5d %d -> %d\n",
                        row,
                        values[ rowOffset + fromNodeColumn ].i,
                        values[ rowOffset + toNodeColumn   ].i ); )

        /*
         * Depth-first (non-tail) recursively mark all upstream flowlines
         * of this flowline (toNode == thisRowFromNode):
         */

        result = 1 +
          flagUpstreamNodes( rows, columns, values,
                             fromNodeColumn, toNodeColumn,
                             thisRowFromNode, mask ); /* Non-tail recursion! */
      }
    }
  }

  POST0( IN_RANGE( result, 0, rows ) );
  return result;
}


#else


/******************************************************************************
PURPOSE: flagUpstreamNodes - Mark a mask array indicating which
         rows are upstream of toNode.
INPUTS:  const int rows                       Number of rows of values.
         const int columns                    Number of columns per row.
         const Value values[ rows * columns ] Array containing from/toNode.
         const int fromNodeColumn             Column of fromNode.
         const int toNodeColumn               Column of toNode.
         const int toNode                     toNode value to match.
         char mask[ rows ]                    Set to 1 if upstream of toNode.
OUTPUTS: char mask[ rows ]                    Set to 1 if upstream of toNode.
RETURNS: int number of upstream rows.
NOTES:   This routine is cycle-safe non-recursive.
******************************************************************************/

static int flagUpstreamNodes( const int rows,
                              const int columns,
                              const Value values[],
                              const int fromNodeColumn,
                              const int toNodeColumn,
                              const int toNode,
                              char mask[] ) {

  PRE08( rows > 0,
         columns > 0,
         values,
         IN_RANGE( fromNodeColumn, 0, columns - 1 ),
         IN_RANGE( toNodeColumn,   0, columns - 1 ),
         fromNodeColumn != toNodeColumn,
         toNode > 0.0,
         mask );

  typedef struct { /* Use an explicit stack in lieu of (non-tail) recursion. */
    int toNode; /* The toNode to search for. */
    int row;    /* Row index to continue searching from. */
  } Entry;
  Entry* stack = NEW( Entry, rows );
  int result = 0;

  DEBUG( fprintf( stderr, "flagUpstreamNodes( toNode = %d )\n", toNode ); )

  if ( stack ) {
    int searchToNode = toNode;
    int stackIndex = 0;
    int row = 0;

    do { /* Loop until stack is empty: */

      /* Check each row's toNode for a match to given searchToNode: */

      for ( ; row < rows; ++row ) {
        CHECK( IN_RANGE( row, 0, rows - 1 ) );

        /* Must compare to 0 to avoid infinite loop in case graph has cycles!*/

        if ( mask[ row ] == 0 ) { /* If row is not already marked: */
          const int rowOffset = row * columns;
          const int thisRowToNode = values[ rowOffset + toNodeColumn ].i;

          if ( thisRowToNode == searchToNode ) {
            const int thisRowFromNode = values[ rowOffset + fromNodeColumn ].i;
            mask[ row ] = 1;
            ++result;
            DEBUG( fprintf( stderr, " found row %5d %d -> %d\n",
                            row,
                            values[ rowOffset + fromNodeColumn ].i,
                            values[ rowOffset + toNodeColumn   ].i ); )

            /*
             * Depth-first mark all upstream flowlines
             * of this flowline (toNode == thisRowFromNode):
             */

            /* First push current search state onto stack: */

            CHECK( IN_RANGE( stackIndex, 0, rows - 1 ) );
            stack[ stackIndex ].toNode = searchToNode;
            stack[ stackIndex ].row    = row;
            ++stackIndex;

            /* Then do a depth-first search for new item: */

            searchToNode = thisRowFromNode; /* Search for this new toNode. */
            row = -1; /* Starting at the beginning again. */
          }
        }
      } /* End loop on rows. */

      --stackIndex;

      if ( stackIndex >= 0 ) {

        /* Pop previous search state off stack and continue searching: */

        CHECK( IN_RANGE( stackIndex, 0, rows - 1 ) );
        searchToNode = stack[ stackIndex ].toNode;
        row          = stack[ stackIndex ].row;
      }

    } while ( stackIndex >= 0 ); /* While stack is not empty. */

    FREE( stack );
  }

  POST0( IN_RANGE( result, 0, rows ) );
  return result;
}


#endif



/******************************************************************************
PURPOSE: copyMatchedLines - Copy lines from an input file to an output file
         if the first value on the line is found in an array of values.
INPUTS:  FILE* input               Input file to read one line from.
         const size_t count        Number of values to match.
         const long long values[]  Array of values to compare to.
OUTPUTS: FILE* output  Output file to write matched lines to.
RETURNS: size_t number of matched lines written to output.
NOTES:   Requires that the values array and the input file are numerically
         sorted (sort -n). Uses a temporary a memory buffer for faster I/O.
******************************************************************************/

static size_t copyMatchedLines( FILE* input,
                                const size_t count,
                                const long long values[],
                                FILE* output ) {
  size_t result = 0;
  const size_t bufferSize = 10 * 1024 * 1024;
  char* buffer = NEW( char, bufferSize + bufferSize );
  int ok = buffer != 0;

  if ( ok ) {
    char* inputBuffer = buffer;
    char* outputBuffer = buffer + bufferSize;
    char* out = outputBuffer;
    char* const end = outputBuffer + bufferSize - 1;
    const long long minimumId = values[ 0 ];
    const long long maximumId = values [ count - 1 ];
    int done = 0;
    *outputBuffer = '\0';
    *end = '\0';

    while ( ! done && ok && ! feof( input ) ) {
      const size_t bytesRead = fread( inputBuffer, 1, bufferSize - 1, input );
      inputBuffer[ bufferSize - 1 ] = '\0';

      if ( bytesRead ) {
        char* line = inputBuffer;
        char* newline = 0;
        inputBuffer[ bytesRead ] = '\0';

        for ( ; ! done && ok && ( newline = strchr( line, '\n' ) );
              line = newline + 1 ) {
          *newline = '\0';

          if ( newline > line ) {
            const long long value = atoll( line );
            int matched = 0;

            if ( value > maximumId ) {
              done = 1;
            } else if ( value >= minimumId ) {
              const long long index = binarySearchAny( value, count, values );
              matched = index >= 0;
            }

            if ( matched ) {
              const size_t lineLength = strlen( line );
              char* const before = newline - 1;

              if ( *before == '\r' ) {
                *before = ' ';
              }

              if ( lineLength ) {

                if ( lineLength >= end - out ) {
                  ok = fputs( outputBuffer, output ) != EOF;
                  *outputBuffer = '\0';
                  out = outputBuffer;
                }

                strcpy( out, line );
                out += lineLength;
                *out = '\n';
                ++out;
                *out = '\0';
                ++result;
              }
            }
          }
        }

        if ( ok && ! done && newline == 0 && *line ) {
          const long partialLineLength = strlen( line );

          ok = partialLineLength < bufferSize - 1 &&
               fseek( input, -partialLineLength, SEEK_CUR ) == 0;
        }
      }
    }

    if ( ok && out != outputBuffer ) {
      ok = fputs( outputBuffer, output ) != EOF;
    }

    if ( ! ok ) {
      result = 0;
    }

    FREE( buffer );
  }

  return result;
}



/******************************************************************************
PURPOSE: compareDate - Compare dates with respect to timestepSize.
INPUTS:  const int yyyymmdd1     1st date to compare.
         const int yyyymmdd2     2nd date to compare.
         const int timestepSize  YEARLY, MONTHLY, DAILY.
RETURNS: int negative if yyyymmdd1 < yyyymmdd2
         or positive if yyyymmdd1 > yyyymmdd2, else 0 if equal
         comparing only portion relevant to timestepSize.
******************************************************************************/

static int compareDate( const int yyyymmdd1, const int yyyymmdd2,
                        const int timestepSize ) {
  PRE03( isValidYearMonthDay( yyyymmdd1 ),
         isValidYearMonthDay( yyyymmdd2 ),
         IN4( timestepSize, DAILY, MONTHLY, YEARLY ) );
  int result = 0;

  if ( timestepSize == YEARLY ) {
    result = yyyymmdd1 / 10000 - yyyymmdd2 / 10000;
  } else if ( timestepSize == MONTHLY ) {
    result = yyyymmdd1 / 100 - yyyymmdd2 / 100;
  } else {
    result = yyyymmdd1 - yyyymmdd2;
  }

  return result;
}



/******************************************************************************
PURPOSE: storeStringValue - Store string value in ShapeData.
INPUTS:  const char* value     String value to duplicate and store.
         ShapeData* shapeData  ShapeData to store string value in.
OUTPUTS: ShapeData* shapeData  ShapeData with stored string value.
RETURNS: const char* pointer to stored string if successful, else 0 and
         a failure message is printed to stderr.
NOTES:   Flyweight Pattern: store only unique strings with shared pointers to
         duplicates to greatly reduce the amount of memory required to store
         many instances of relatively few unique strings.
******************************************************************************/

static const char* storeStringValue( const char* value, ShapeData* shapeData) {

  PRE04( value, shapeData, shapeData->capacity > 0, shapeData->stringStorage );

  const char* result = 0;
  int index = 0;
  int ok = 0;
  const int capacity = shapeData->capacity;
  char** stringStorage = shapeData->stringStorage;
  CHECK2( capacity > 0, stringStorage != 0 );

  DEBUG2( fprintf( stderr,
                   "storeStringValue( %s, capacity = %d, [%s, %s ...] )\n",
                   value, capacity,
                   stringStorage[ 0 ],
                   stringStorage[ 1 ] ); )

  /* Linear search for matching string or empty slot: */

  while ( AND3( index < capacity,
                stringStorage[ index ],
                strcmp( stringStorage[ index ], value ) ) ) {
    ++index;
  }

  DEBUG2( fprintf( stderr, "index = %d\n", index ); )

  if ( index < capacity ) { /* Found match or empty slot for new string: */

    if ( stringStorage[ index ] == 0 ) { /* Store copy into empty slot: */
      stringStorage[ index ] = strdup( value );
    }

    ok = stringStorage[ index ] != 0;

    if ( ! ok ) {
      fprintf(stderr, "Failed to allocate enough memory to copy text data.\n");
    }

    CHECK( IMPLIES( ok, ! strcmp( shapeData->stringStorage[ index ], value )));
  } else { /* Must grow stringStorage. */
    const int growthAmount = capacity / 2; /* Grow storage by 50%. */
    const int newCapacity = capacity + growthAmount;
    char** newStrings = newCapacity > 0 ? NEW( char*, newCapacity ) : 0;

    if ( newStrings ) {

      /* Transfer pointers to strings from old storage to new storage: */

      for ( index = 0; index < capacity; ++index ) {
        newStrings[ index ] = stringStorage[ index ];
        stringStorage[ index ] = 0;
      }

      FREE( shapeData->stringStorage );
      stringStorage = 0;
      shapeData->stringStorage = newStrings;
      newStrings = 0;
      shapeData->capacity = newCapacity;
      shapeData->columnNames = shapeData->stringStorage; /* Alias. */

      CHECK2( IN_RANGE( index, 0, shapeData->capacity - 1 ),
              shapeData->stringStorage[ index ] == 0 );

      shapeData->stringStorage[ index ] = strdup( value );
      ok = shapeData->stringStorage[ index ] != 0;
      CHECK( IMPLIES( ok, ! strcmp( shapeData->stringStorage[index], value )));

      if ( ! ok ) {
        fprintf( stderr,
                 "Failed to allocate enough memory to copy text data.\n" );
      }
    }
  }

  if ( ok ) {
    CHECK( IN_RANGE( index, 0, shapeData->capacity - 1 ) );
    result = shapeData->stringStorage[ index ];
  }

  POST0( IMPLIES( result,
                  AND3( IN_RANGE( index, 0, shapeData->capacity - 1 ),
                        ! strcmp( shapeData->stringStorage[ index ], value ),
                        result == shapeData->stringStorage[ index ] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: copyStringAttribute - Copy a string attribute from one DBF file to
         another.
INPUTS:  DBFHandle inputFile      File to read column attribute from.
         int inputRow             Row index to read (0-based).
         int inputColumn          Column index to read (0-based).
         DBFHandle outputFile     File to write column attribute to.
         int outputRow            Row index to write (0-based).
         int outputColumn         Column index to write (0-based).
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int copyStringAttribute( DBFHandle inputFile,
                                int inputRow, int inputColumn,
                                DBFHandle outputFile,
                                int outputRow, int outputColumn ) {
  const char* const value0 =
    DBFReadStringAttribute( inputFile, inputRow, inputColumn );
  const char* const value = value0 ? value0 : "NULL";
  const int result =
    DBFWriteStringAttribute( outputFile, outputRow, outputColumn, value );
  return result;
}



/******************************************************************************
PURPOSE: copyIntegerAttribute - Copy an integer attribute from one DBF file to
         another.
INPUTS:  DBFHandle inputFile      File to read column attribute from.
         int inputRow             Row index to read (0-based).
         int inputColumn          Column index to read (0-based).
         DBFHandle outputFile     File to write column attribute to.
         int outputRow            Row index to write (0-based).
         int outputColumn         Column index to write (0-based).
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int copyIntegerAttribute( DBFHandle inputFile,
                                 int inputRow, int inputColumn,
                                 DBFHandle outputFile,
                                 int outputRow, int outputColumn ) {
  const int value =
    DBFReadIntegerAttribute( inputFile, inputRow, inputColumn );
  const int result =
    DBFWriteIntegerAttribute( outputFile, outputRow, outputColumn, value );
  return result;
}



/******************************************************************************
PURPOSE: copyDoubleAttribute - Copy a double attribute from one DBF file to
         another.
INPUTS:  DBFHandle inputFile      File to read column attribute from.
         int inputRow             Row index to read (0-based).
         int inputColumn          Column index to read (0-based).
         DBFHandle outputFile     File to write column attribute to.
         int outputRow            Row index to write (0-based).
         int outputColumn         Column index to write (0-based).
         int invalidIfNegative    Map all negative values to missing?
         double missing           Value to map all negative values to.
         double offset            Value to add to value read.
         double scale             Value to multiply offset-value read.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int copyDoubleAttribute( DBFHandle inputFile,
                                int inputRow, int inputColumn,
                                DBFHandle outputFile,
                                int outputRow, int outputColumn,
                                int invalidIfNegative, double missing,
                                double offset, double scale ) {
  const double value0 =
    DBFReadDoubleAttribute( inputFile, inputRow, inputColumn );
  const double value =
    value0 < -9999.0 ? missing
    : AND2( invalidIfNegative, value0 < 0.0 ) ? missing
    : ( value0 + offset ) * scale;
  const int result =
    DBFWriteDoubleAttribute( outputFile, outputRow, outputColumn, value );
  return result;
}



/******************************************************************************
PURPOSE: computeSparsedPartCount - Compute sparsed polyline/polygon part count.
INPUTS:  const SHPObject* shape  Shape to examine.
         const double minimumAdjacentVertexDistance  Adjacent vertices closer
                                                     than this
                                                     (in either x or y) will
                                                     be merged.
         const int minimumSparsedVertices  2 if polyline, 3 if polygon.
RETURNS: int sparsed part count.
******************************************************************************/

static int computeSparsedPartCount( const SHPObject* shape,
                                    const double minimumAdjacentVertexDistance,
                                    const int minimumSparsedVertices ) {

  PRE012( shape,
          IN5( shape->nSHPType,
               SHPT_POLYGON, SHPT_POLYGONZ, SHPT_ARC, SHPT_ARCZ ),
          shape->nShapeId >= 0,
          shape->nParts > 0,
          shape->panPartStart,
          shape->panPartType,
          shape->panPartType[ 0 ] == SHPP_RING,
          shape->nVertices > 0,
          shape->padfX,
          shape->padfY,
          minimumAdjacentVertexDistance >= 0.0,
          IN3( minimumSparsedVertices, 2, 3 ) );

  int result = 0;
  const int parts = shape->nParts;
  int part = 0;

  /* Compute sparsed number of parts (with at least minimumSparsedVertices): */

  for ( part = 0; part < parts; ++part ) {
    const int partVertexCount =
      parts == 1 ? shape->nVertices
      : part < parts - 1 ?
        shape->panPartStart[ part + 1 ] - shape->panPartStart[ part ]
      : shape->nVertices - shape->panPartStart[ part ];
    const int offset = shape->panPartStart[ part ];
    const double* const x = shape->padfX + offset;
    const double* const y = shape->padfY + offset;
    double x0 = x[ 0 ];
    double y0 = y[ 0 ];
    int sparsedVertices = 1; /* Count first vertex. */
    int vertex = 1;

    for ( ; AND2( vertex < partVertexCount,
                  sparsedVertices < minimumSparsedVertices );
          ++vertex ) {
      const double x1 = x[ vertex ];
      const double y1 = y[ vertex ];
      const double deltaX = x1 > x0 ? x1 - x0 : x0 - x1;

      if ( deltaX >= minimumAdjacentVertexDistance ) {
        ++sparsedVertices;
        x0 = x1;
        y0 = y1;
      } else {
        const double deltaY = y1 > y0 ? y1 - y0 : y0 - y1;

        if ( deltaY >= minimumAdjacentVertexDistance ) {
          ++sparsedVertices;
          x0 = x1;
          y0 = y1;
        }
      }
    }

    result += sparsedVertices == minimumSparsedVertices;
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: computeSparsedVertexCount - Compute sparsed polyline/polygon vertex
                                     count.
INPUTS:  const int vertexCount          Number of vertices in part.
         const double x[ vertexCount ]  X-coordinate of part.
         const double y[ vertexCount ]  Y-coordinate of part.
         const double minimumAdjacentVertexDistance  Adjacent vertices closer
                                                     than this
                                                     (in either x or y) will
                                                     be merged.
         const int isPolygon            1 if is polygon, else 0.
RETURNS: int sparsed vertex count of part.
******************************************************************************/

static int computeSparsedVertexCount( const int vertexCount,
                                      const double x[],
                                      const double y[],
                                      const double minimumAdjacentVertexDistance,
                                      const int isPolygon ) {

  PRE05( vertexCount > 0, x, y, minimumAdjacentVertexDistance >= 0.0,
         IS_BOOL( isPolygon ) );

  int result = 1; /* Count first vertex. */
  double x0 = x[ 0 ];
  double y0 = y[ 0 ];
  int vertex = 1;

  for ( ; vertex < vertexCount; ++vertex ) {
    const double x1 = x[ vertex ];
    const double y1 = y[ vertex ];
    const double deltaX = x1 > x0 ? x1 - x0 : x0 - x1;

    if ( deltaX >= minimumAdjacentVertexDistance ) {
      ++result;
      x0 = x1;
      y0 = y1;
    } else {
      const double deltaY = y1 > y0 ? y1 - y0 : y0 - y1;

      if ( deltaY >= minimumAdjacentVertexDistance ) {
        ++result;
        x0 = x1;
        y0 = y1;
      }
    }
  }

  if ( AND2( isPolygon, ! AND2( x0 == x[ 0 ], y0 == y[ 0 ] ) ) ) {
    ++result; /* Replicate first vertex to close polygon. */
  }

  POST0( result >= 1 );
  return result;
}



/******************************************************************************
PURPOSE: copySparseVertices - Copy sparsed polyline/polygon vertices and also
         compute bounds of sparsed vertices.
INPUTS:  const int vertexCount          Number of vertices in part.
         const int sparsedVertices      Number of resulting sparsed vertices.
         const double x[ vertexCount ]  X-coordinates of part.
         const double y[ vertexCount ]  Y-coordinates of part.
         const double minimumAdjacentVertexDistance  Adjacent vertices closer
                                                     than this
                                                     (in either x or y) will
                                                     be merged.
         const int isPolygon            1 if is polygon, else 0.
         const int* initializedBounds    Have bounds been initialized?
         Bounds bounds                   Bounds of sparsed vertices.
OUTPUTS: const int* initializedBounds    Have bounds been initialized?
         Bounds bounds                   Updated bounds of sparsed vertices.
         gpc_vertex* const vertices      GPC vertices copied into.
******************************************************************************/

static void copySparseVertices( const int vertexCount,
                                const int sparsedVertices,
                                const double x[],
                                const double y[],
                                const double minimumAdjacentVertexDistance,
                                const int isPolygon,
                                int* const initializedBounds,
                                Bounds bounds,
                                gpc_vertex* const vertices ) {

  PRE010( vertexCount > 0, sparsedVertices > 0, x, y,
          minimumAdjacentVertexDistance >= 0.0,
          IS_BOOL( isPolygon ), initializedBounds,
          IS_BOOL( *initializedBounds ), bounds, vertices );

  double xMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  double xMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  double yMinimum = bounds[ LATITUDE ] [ MINIMUM ];
  double yMaximum = bounds[ LATITUDE  ][ MAXIMUM ];
  int vertex = 0;
  int sparsedVertex = 1; /* Count first vertex always copied. */
  double x0 = x[ 0 ];
  double y0 = y[ 0 ];
  vertices[ 0 ].x = x0; /* Copy first vertex. */
  vertices[ 0 ].y = y0;

  for ( vertex = 1; vertex < vertexCount; ++vertex ) {
    const double x1 = x[ vertex ];
    const double y1 = y[ vertex ];
    const double deltaX = x1 > x0 ? x1 - x0 : x0 - x1;
    int nonCoincident = deltaX >= minimumAdjacentVertexDistance;

    if ( ! nonCoincident ) {
      const double deltaY = y1 > y0 ? y1 - y0 : y0 - y1;
      nonCoincident = deltaY >= minimumAdjacentVertexDistance;
    }

    if ( nonCoincident ) {

      /* Copy vertex: */

      CHECK( IN_RANGE( sparsedVertex, 0, sparsedVertices - 1 ) );
      vertices[ sparsedVertex ].x = x1;
      vertices[ sparsedVertex ].y = y1;
      ++sparsedVertex;
      x0 = x1;
      y0 = y1;

      /* Update bounds: */

      if ( ! *initializedBounds ) {
        *initializedBounds = 1;
        xMinimum = xMaximum = x0;
        yMinimum = yMaximum = y0;
      }

      if ( x1 < xMinimum ) {
        xMinimum = x1;
      } else if ( x1 > xMaximum ) {
        xMaximum = x1;
      }

      if ( y1 < yMinimum ) {
        yMinimum = y1;
      } else if ( y1 > yMaximum ) {
        yMaximum = y1;
      }
    }
  }

  /* Update bounds: */

  bounds[ LONGITUDE ][ MINIMUM ] = xMinimum;
  bounds[ LONGITUDE ][ MAXIMUM ] = xMaximum;
  bounds[ LATITUDE  ][ MINIMUM ] = yMinimum;
  bounds[ LATITUDE  ][ MAXIMUM ] = yMaximum;

  /* If needed, replicate first vertex to close polygon. */

  if ( AND2( isPolygon, sparsedVertex < sparsedVertices ) ) {
    vertices[ sparsedVertex ].x = x[ 0 ];
    vertices[ sparsedVertex ].y = y[ 0 ];
    ++sparsedVertex;
  }

  CHECK( sparsedVertex == sparsedVertices );
  POST02( *initializedBounds, isValidBounds( (const double (*)[2]) bounds ) );
}



/******************************************************************************
PURPOSE: computeVertexBounds - Compute bounds of vertices.
INPUTS:  const int count          Number of vertices.
OUTPUTS: const float xy[          Sequence of (x, y) vertices.
         double xyRange[ 2 * 2 ]  xMin, yMin, xMax, yMax.
******************************************************************************/

static void computeVertexBounds( const int count, const float xy[],
                                 double xyRange[] ) {

  const float* x = xy;
  const float* y = xy + 1;
  float xMinimum = *x;
  float xMaximum = xMinimum;
  float yMinimum = *y;
  float yMaximum = yMinimum;
  int vertex = 0;

  for ( vertex = 0; vertex < count; ++vertex, x += 2, y += 2 ) {

    if ( *x < xMinimum ) {
      xMinimum = *x;
    } else if ( *x > xMaximum ) {
      xMaximum = *x;
    }

    if ( *y < yMinimum ) {
      yMinimum = *y;
    } else if ( *y > yMaximum ) {
      yMaximum = *y;
    }
  }

  xyRange[ 0 ] = xMinimum;
  xyRange[ 1 ] = yMinimum;
  xyRange[ 2 ] = xMaximum;
  xyRange[ 3 ] = yMaximum;
}



/******************************************************************************
PURPOSE: computeGridBounds - Compute bounds of grid coordinates.
INPUTS:  const int rows            Number of grid rows.
         const int columns         Number of grid columns.
         const double westEdge     Distance from origin to west edge of grid.
         const double southEdge    Distance from origin to south edge of ".
         const double cellWidth    Width of each grid cell (e.g., 12000 m).
         const double cellWHeight  Height of each grid cell (e.g., 12000 m).
         Unproject unproject       Function to unproject (x,y) to (lon,lat).
OUTPUTS: double range[ 4 ]         xMin, yMin, xMax, yMax.
******************************************************************************/

static void computeGridBounds( const int rows,
                               const int columns,
                               const double westEdge,
                               const double southEdge,
                               const double cellWidth,
                               const double cellHeight,
                               Unproject unproject,
                               double range[ 4 ] ) {

  PRE05( rows > 0,
         columns > 0,
         rows * columns > 0,
         IMPLIES( unproject == 0,
                   AND2 ( isValidLongitudeLatitude( westEdge, southEdge ),
                         isValidLongitudeLatitude(
                                     westEdge  + columns * cellWidth,
                                     southEdge + rows    * cellHeight ) ) ),
         range );

  const int rows_1 = rows + 1;
  const int columns_1 = columns + 1;
  double x = westEdge;
  double y = southEdge;
  double px = x;
  double py = y;
  double minimumX = px;
  double maximumX = minimumX;
  double minimumY = py;
  double maximumY = minimumY;
  int row = 0;
  int column = 0;

  /* Initialize to west/south point: */

  if ( unproject ) {
    unproject( x, y, &px, &py );
    minimumX = maximumX = px;
    minimumY = maximumY = py;
  }

  /* Check against rest of west edge of grid: */

  y += cellHeight;

  for ( row = 1; row < rows_1; ++row, y += cellHeight ) {
    updatePointMinmax( unproject, x, y,
                       &minimumX, &maximumX, &minimumY, &maximumY );
  }

  /* Check against east edge of grid: */

  x = westEdge + columns * cellWidth;
  y = southEdge;

  for ( row = 0; row < rows_1; ++row, y += cellHeight ) {
    updatePointMinmax( unproject, x, y,
                       &minimumX, &maximumX, &minimumY, &maximumY );
  }

  /* Check against rest of south edge of grid: */

  x = westEdge + cellWidth;
  y = southEdge;

  for ( column = 1; column < columns_1; ++column, x += cellWidth ) {
    updatePointMinmax( unproject, x, y,
                       &minimumX, &maximumX, &minimumY, &maximumY );
  }

  /* Check against north edge of grid: */

  x = westEdge;
  y = southEdge + rows * cellHeight;

  for ( column = 0; column < columns_1; ++column, x += cellWidth ) {
    updatePointMinmax( unproject, x, y,
                       &minimumX, &maximumX, &minimumY, &maximumY );
  }

  range[ 0 ] = minimumX;
  range[ 1 ] = minimumY;
  range[ 2 ] = maximumX;
  range[ 3 ] = maximumY;

  POST04( isValidLongitudeLatitude( range[ 0 ], range[ 1 ] ),
          isValidLongitudeLatitude( range[ 2 ], range[ 3 ] ),
          range[ 0 ] <= range[ 2 ],
          range[ 1 ] <= range[ 3 ] );
}



/******************************************************************************
PURPOSE: updatePointMinmax - Update minimum and maximum of grid coordinates.
INPUTS:  Unproject unproject       Function to unproject (x,y) to (lon,lat).
         const double x            X-coordinate to check.
         const double y            Y-coordinate to check.
         double* minimumX          Current minimum X-coordinate.
         double* maximumX          Current maximum X-coordinate.
         double* minimumY          Current minimum Y-coordinate.
         double* maximumY          Current maximum Y-coordinate.
OUTPUTS: double* minimumX          Updated minimum X-coordinate.
         double* maximumX          Updated maximum X-coordinate.
         double* minimumY          Updated minimum Y-coordinate.
         double* maximumY          Updated maximum Y-coordinate.
******************************************************************************/

static void updatePointMinmax( Unproject unproject,
                               const double x,
                               const double y,
                               double* minimumX,
                               double* maximumX,
                               double* minimumY,
                               double* maximumY ) {

  PRE05( minimumX, maximumX, minimumY, maximumY,
         IMPLIES( unproject == 0,
                   AND3( isValidLongitudeLatitude( x, y ),
                         isValidLongitudeLatitude(*minimumX, *minimumY ),
                         isValidLongitudeLatitude(*maximumX, *maximumY ) ) ) );

  double px = x;
  double py = y;

  if ( unproject ) {
    unproject( x, y, &px, &py );
  }

  if ( px < *minimumX ) {
    *minimumX = px;
  } else if ( px > *maximumX ) {
    *maximumX = px;
  }

  if ( py < *minimumY ) {
    *minimumY = py;
  } else if ( py > *maximumY ) {
    *maximumY = py;
  }
}



/******************************************************************************
PURPOSE: computeGridCellCenters - Compute longitude-latitude of each grid cell
         center.
INPUTS:  const int rows            Number of grid rows.
         const int columns         Number of grid columns.
         const double westEdge     Distance from origin to west edge of grid.
         const double southEdge    Distance from origin to south edge of ".
         const double cellWidth    Width of each grid cell (e.g., 12000 m).
         const double cellWHeight  Height of each grid cell (e.g., 12000 m).
         Unproject unproject       Function to unproject (x,y) to (lon,lat).
OUTPUTS: float lonlats[ rows * columns * 2 ]   Lon-lats of grid cell centers.
******************************************************************************/

static void computeGridCellCenters( const int rows,
                                    const int columns,
                                    const double westEdge,
                                    const double southEdge,
                                    const double cellWidth,
                                    const double cellHeight,
                                    Unproject unproject,
                                    float lonlats[] ) {

  PRE05( rows > 0,
         columns > 0,
         rows * columns > 0,
         IMPLIES( unproject == 0,
                   AND2 ( is_valid_longitude_latitude( westEdge, southEdge ),
                          is_valid_longitude_latitude(
                                     westEdge  + columns * cellWidth,
                                     southEdge + rows    * cellHeight ) ) ),
         lonlats );

  float* coordinates = lonlats;
  double y = southEdge;
  int row = 0;

  for ( row = 0; row < rows; ++row, y += cellHeight ) {
    const double centerY = ( y + y + cellHeight ) * 0.5;
    int column = 0;
    double x = westEdge;

    for ( column = 0; column < columns; ++column, x += cellWidth ) {
      const double centerX = ( x + x + cellWidth ) * 0.5;
      double px = x;
      double py = y;

      if ( unproject ) {
        unproject( centerX, centerY, &px, &py );
      }

      *coordinates++ = px;
      *coordinates++ = py;
    }
  }

  POST02( is_valid_longitude_latitude( lonlats[ 0 ], lonlats[ 1 ] ),
          is_valid_longitude_latitude( lonlats[ rows * columns * 2 - 2 ],
                                       lonlats[ rows * columns * 2 - 1 ] ) );
}



/******************************************************************************
PURPOSE: computePointBounds - Compute bounds of point coordinates.
INPUTS:  const int points                   Number of points.
         const float lonlats[ points * 2 ]  Longitude latitude coordinates.
OUTPUTS: double range[ 4 ]                  xMin, yMin, xMax, yMax.
******************************************************************************/

static void computePointBounds( const int points,
                                const float lonlats[],
                                double range[ 4 ] ) {
  double minimumX = lonlats[ 0 ];
  double maximumX = minimumX;
  double minimumY = lonlats[ 1 ];
  double maximumY = minimumY;
  int point = 0;

  /* Check against rest of points: */

  for ( point = 1; point < points; ++point ) {
    const int point2 = point + point;
    const float longitude = lonlats[ point2 ];
    const float latitude = lonlats[ point2 + 1 ];

    if ( longitude < minimumX ) {
      minimumX = longitude;
    } else if ( longitude > maximumX ) {
      maximumX = longitude;
    }

    if ( latitude < minimumY ) {
      minimumY = latitude;
    } else if ( latitude > maximumY ) {
      maximumY = latitude;
    }
  }

  range[ 0 ] = minimumX;
  range[ 1 ] = minimumY;
  range[ 2 ] = maximumX;
  range[ 3 ] = maximumY;
}



/******************************************************************************
PURPOSE: computeRange - Compute minimum and maximum value in array.
INPUTS:  const double array[ count ]  Items to examine.
         const int count              Number of items to compare.
         const int stride             Stride of index.
OUTPUTS: double range[ 2 ]            Minimum and maximum values in array.
******************************************************************************/

static void computeRange( const double array[], const int count,
                          const int stride, double range[ 2 ] ) {
  double minimum = array[ 0 ];
  double maximum = minimum;
  int index = 0;

  for ( index = stride; index < count; index += stride ) {
    const double item = array[ index ];

    if ( item < minimum ) {
      minimum = item;
    } else if ( item > maximum ) {
      maximum = item;
    }
  }

  range[ 0 ] = minimum;
  range[ 1 ] = maximum;
}



/******************************************************************************
PURPOSE: convertTimestamp - Convert HMS smoke timestamp to yyyyddd, hhmm.
INPUTS:  const char* inputFileName  File name of input: hms_smoke20161017.dbf
         const char* timestamp      Timestamp hhmm or yyyyddd hhmm.
OUTPUTS: int* const yyyyddd         Date.
         int* const hhmm            Time.
RETURNS: int 1 if successful, else 0 failure message is printed.
******************************************************************************/

static int convertTimestamp( const char* inputFileName,
                             const char* timestamp,
                             int* const yyyyddd, int* const hhmm ) {

  PRE05( inputFileName, *inputFileName, timestamp, yyyyddd, hhmm );
  int result = 0;

  if ( strlen( timestamp ) == 12 ) { /* YYYYDDD HHMM: */
    result = sscanf( timestamp, "%d %d", yyyyddd, hhmm ) == 2;

    if ( result ) {
      result = isValidDate( *yyyyddd );

      if ( result ) {
        result = isValidTime( *hhmm * 100 );
      }
    }
  } else { /* HHMM only. Parse YYYYMMDD from file name & convert to YYYYDDD: */
    const char* const tag = "hms_smoke";
    const char* name = strstr( inputFileName, tag );

    if ( name ) {
      name += strlen( tag );

      if ( isdigit( *name ) ) {
        const int yyyymmdd = atoi( name );
        result = isValidYearMonthDay( yyyymmdd );

        if ( result ) {
          *yyyyddd = convertYearMonthDay( yyyymmdd );
          *hhmm = atoi( timestamp );
          result = isValidTime( *hhmm * 100 );
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeShort - Write a 2-byte integer to a byte array with the given
                    endian format.
INPUTS:  const int index        Starting index into bytes[] to write.
         const short value      Value to write.
         const int endian       Endian format to write.
OUTPUTS: unsigned char bytes[]  bytes[ index ... index + 4 ].
RETURNS: int index + 4.
******************************************************************************/

static int writeShort( unsigned char bytes[], const int index,
                      const short value, const int endian ) {
  const int result = index + 4;
  const unsigned char byte1 =   value & 0x00ff;
  const unsigned char byte2 = ( value & 0xff00 ) >> 8;

  if ( endian == BIG ) {
    bytes[ index     ] = byte2;
    bytes[ index + 1 ] = byte1;
  } else {
    bytes[ index     ] = byte1;
    bytes[ index + 1 ] = byte2;
  }

  return result;
}



/******************************************************************************
PURPOSE: writeInt - Write a 4-byte integer to a byte array with the given
                    endian format.
INPUTS:  const int index        Starting index into bytes[] to write.
         const int value        Value to write.
         const int endian       Endian format to write.
OUTPUTS: unsigned char bytes[]  bytes[ index ... index + 4 ].
RETURNS: int index + 4.
******************************************************************************/

static int writeInt( unsigned char bytes[], const int index,
                     const int value, const int endian ) {
  const int result = index + 4;
  const unsigned char byte1 =   value & 0x000000ff;
  const unsigned char byte2 = ( value & 0x0000ff00 ) >> 8;
  const unsigned char byte3 = ( value & 0x00ff0000 ) >> 16;
  const unsigned char byte4 = ( value & 0xff000000 ) >> 24;

  if ( endian == BIG ) {
    bytes[ index     ] = byte4;
    bytes[ index + 1 ] = byte3;
    bytes[ index + 2 ] = byte2;
    bytes[ index + 3 ] = byte1;
  } else {
    bytes[ index     ] = byte1;
    bytes[ index + 1 ] = byte2;
    bytes[ index + 2 ] = byte3;
    bytes[ index + 3 ] = byte4;
  }

  return result;
}



/******************************************************************************
PURPOSE: writeDouble - Write a 8-byte real to a byte array with the given
                    endian format.
INPUTS:  const int index        Starting index into bytes[] to write.
         const double value     Value to write.
         const int endian       Endian format to write.
OUTPUTS: unsigned char bytes[]  bytes[ index ... index + 8 ].
RETURNS: int index + 8.
******************************************************************************/

static int writeDouble( unsigned char bytes[], const int index,
                        const double value, const int endian ) {
  const int result = index + 8;
  const void* const pvvalue      = &value;
  const long long* const pivalue = pvvalue;
  const long long ivalue         = *pivalue;
  const unsigned char byte1 =   ivalue & 0x00000000000000ffLL;
  const unsigned char byte2 = ( ivalue & 0x000000000000ff00LL ) >> 8;
  const unsigned char byte3 = ( ivalue & 0x0000000000ff0000LL ) >> 16;
  const unsigned char byte4 = ( ivalue & 0x00000000ff000000LL ) >> 24;
  const unsigned char byte5 = ( ivalue & 0x000000ff00000000LL ) >> 32;
  const unsigned char byte6 = ( ivalue & 0x0000ff0000000000LL ) >> 40;
  const unsigned char byte7 = ( ivalue & 0x00ff000000000000LL ) >> 48;
  const unsigned char byte8 = ( ivalue & 0xff00000000000000LL ) >> 56;

  if ( endian == BIG ) {
    bytes[ index     ] = byte8;
    bytes[ index + 1 ] = byte7;
    bytes[ index + 2 ] = byte6;
    bytes[ index + 3 ] = byte5;
    bytes[ index + 4 ] = byte4;
    bytes[ index + 5 ] = byte3;
    bytes[ index + 6 ] = byte2;
    bytes[ index + 7 ] = byte1;
  } else {
    bytes[ index     ] = byte1;
    bytes[ index + 1 ] = byte2;
    bytes[ index + 2 ] = byte3;
    bytes[ index + 3 ] = byte4;
    bytes[ index + 4 ] = byte5;
    bytes[ index + 5 ] = byte6;
    bytes[ index + 6 ] = byte7;
    bytes[ index + 7 ] = byte8;
  }

  return result;
}



/******************************************************************************
PURPOSE: minimumInt - Smallest value in array.
INPUTS:  int count                 Number of values in array.
         const int array[ count ]  Values to check.
RETURNS: int smallest value in array.
******************************************************************************/

static int minimumInt( int count, const int array[] ) {
  PRE02( count > 0, array );
  int index = 0;
  int result = array[ 0 ];

  for ( index = 1; index < count; ++index ) {
    const int value = array[ index ];
    result = result < value ? result : value;
  }

  POST02( result <= array[ 0 ], result <= array[ count - 1 ] );
  return result;
}



/******************************************************************************
PURPOSE: maximumInt - Largest value in array.
INPUTS:  int count                 Number of values in array.
         const int array[ count ]  Values to check.
RETURNS: int largest value in array.
******************************************************************************/

static int maximumInt( int count, const int array[] ) {
  PRE02( count > 0, array );
  int index = 0;
  int result = array[ 0 ];

  for ( index = 1; index < count; ++index ) {
    const int value = array[ index ];
    result = result > value ? result : value;
  }

  POST02( result >= array[ 0 ], result >= array[ count - 1 ] );
  return result;
}



/******************************************************************************
PURPOSE: allocate - Allocate memory by calling malloc() to allocate
         (then zero) an array of count items, each of size bytesEach.
         If allocation fails then a message is printed to stderr.
INPUTS:  size_t sizeEach  Number of bytes per item.
         size_t count     Number of items to allocate.
RETURNS: void*  The resulting address (zeroed), or 0 if unsuccessful.
******************************************************************************/

static void* allocate( size_t bytesEach, size_t count ) {
  PRE02( bytesEach > 0, count > 0 );
  void* result = 0;

  if ( bytesEach > ULONG_MAX / count || count > ULONG_MAX / bytesEach ) {
    fprintf( stderr, "\a\nCannot allocate %lu x %lu bytes ", count, bytesEach);
  } else {
    const size_t bytes = bytesEach * count;
    result = malloc( bytes );

    if ( ! result ) {
      fprintf(stderr, "\a\nCannot allocate %lu x %lu bytes ", count,bytesEach);
      perror( "because" );
    } else {
      memset( result, 0, bytes );
    }
  }

  return result;
}




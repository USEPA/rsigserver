
/******************************************************************************
PURPOSE: CMAQSubset.c - Read a subset of a CMAQ file and write it to stdout as
         XDR (IEEE-754)FORMAT binary, ASCII spreadsheet or NetCDF file.
NOTES:   https://cmascenter.org/ioapi/documentation/all_versions/html/GRIDS.html
HISTORY: 2005-08-01 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>     /* For stderr, fprintf(), tempnam(). */
#include <string.h>    /* For memset(), strcmp().  */
#include <stdlib.h>    /* For malloc(), free(), atoi(), atof().  */
#include <limits.h>    /* For INT_MAX. */
#include <math.h>      /* For exp(), log(). */
#include <unistd.h>    /* For unlink(), getpid().  */

#include <Assertions.h>      /* For PRE(), POST(), AND2(), DEBUG(), etc. */
#include <Utilities.h>       /* For LONGITUDE,MINIMUM,Bounds,parseOptions(). */
#include <Projector.h>       /* For Projector, isValidEllipsoid()*/
#include <Albers.h>          /* For Albers, newAlbers()*/
#include <Lambert.h>         /* For Lambert, newLambert()*/
#include <Mercator.h>        /* For Mercator, newMercator()*/
#include <Stereographic.h>   /* For Stereographic, newStereographic(). */
#include <NetCDFUtilities.h> /* For createNetCDFFile(), etc. */

/*================================== TYPES ==================================*/

/*
 * METCRO3D WWIND is always valid but CCTM_CONC W_VEL is sometimes BADVAL3 for
 * the first timestep of a file.
 * This is problematic for vis and particle tracing / back-trajectory, etc.
 * Setting the flag below to 1 will replace BADVAL3 vertical wind velocities
 * with 0 for more useful downstream applications.
 */

#define ZERO_BAD_WWIND 1

/*
 * When subsetting grid by bounds,
 * 0 = Just check that grid cell axis-aligned bounds overlap subset bounds and
 * include grid rows/columns with at least one such grid cell.
 * 1 = If subset bounds does not subsume the axis-aligned cell bounds,
 *     then also test clip such grid cells quadrilaterals to subset bounds and
 *     only include grid rows/columns with at least one such clipped grid cell.
 * 1 is slightly slower but may yield a slightly smaller subset if grid cells
 * are large.
 */

#define TEST_CLIP_GRID_CELLS 1

static const char* const temporaryFilePrefix = "CMAQSubset.";
static char temporaryFileName[ 256 ] = "";
static const char* const defaultNote = "https://www.epa.gov/cmaq,CMAQSubset";
static const char* const defaultTemporaryDirectory = ".";

enum { MAX_FILES = 512 };

/* From M3IO specification: Note CCTM_CONC file has 258 variables! */

enum { NAMLEN3 = 16, MXDLEN3 = 80, MXVARS3 = /* 120 */ 512, MXLAYS3 = 100 };
enum { LATGRD3 = 1, LAMGRD3 = 2, POLGRD3 = 6, EQMGRD3 = 7, ALBGRD3 = 9 };
enum {
  IMISS3  = -9999, /* None. */
  VGSGPH3 = 1,     /* Hydrostatic sigma-P. */
  VGSGPN3,         /* Non-hydrostatic sigma-P. */
  VGSIGZ3,         /* Sigma-Z. */
  VGPRES3,         /* Pressure (pascals). */
  VGZVAL3,         /* Z (m) (above sea level). */
  VGHVAL3,         /* Z (m) (above terrain). */
  VGWRFEM          /* Sigma-P WRF. */
};

#define ELEVATION_MINIMUM (-1e3)
#define ELEVATION_MAXIMUM 1e5


#define IS_VALID_VERTICAL_GRID_TYPE( vgtyp ) \
  IN8( vgtyp, VGSGPH3, VGSGPN3, VGSIGZ3, VGPRES3, VGZVAL3, VGHVAL3, VGWRFEM )

#define ELLIPSOID_MINIMUM 6e6
#define ELLIPSOID_MAXIMUM 7e6
#define DEFAULT_EARTH_RADIUS 6370000.0

/* Command-line arguments: */

typedef struct {
  int          fileCount;                   /* Count of input files.         */
  int          zfFileCount;                 /* Count of zf files.            */
  int          wwindFileCount;              /* Count of files with WWIND.    */
  int          format;                      /* FORMAT_XDR..FORMAT_IOAPI.     */
  int          auxMode;                     /* 0, INTEGRATE, WIND.           */
  int          aggregateMode;               /* 0, ... AGGREGATE_DAILY_MAX8.  */
  int          lonlat;                      /* Output LONGITUDE, LATITUDE?   */
  int          elevation;                   /* Output ELEVATION?             */
  int          variables;                   /* Count of variables to output. */
  int          subset[DIMENSIONS][2];  /* 1-based [COLUMN..TIME][MIN/MAXIMUM]*/
  const char*  fileNames[ MAX_FILES ];      /* Array of input files.         */
  const char*  htFileName;                  /* Name of file with LON,LAT,HT. */
  const char*  zfFileNames[ MAX_FILES ];    /* Array of files with ZH,ZF,DENS*/
  const char*  wwindFileNames[ MAX_FILES ]; /* Array of files with WWIND.    */
  const char*  note;                        /* File description text.        */
  const char*  tmpDir;                      /* Directory to write temp file. */
  const char*  outputFileName;              /* Optional. stdout if 0.        */
  const char*  lsDir;                       /* Directory to list.            */
  char         variableNames[ MXVARS3 ][ NAMLEN3 + 1 ]; /* Output variables. */
  double       ellipsoid[ 2 ];              /* Equitorial, polar radius (m). */
  Bounds       bounds;                      /* Subset lon-lat bounds.        */
} Arguments;

typedef struct {
  Arguments* arguments;  /* Command-line arguments. */
  /* *TimeRange[ file ][ first, last, timesteps, hoursPerTimestep ]: */
  int fileTimeRange[      MAX_FILES ][ 4 ]; /* Of each data file. */
  int zfFileTimeRange[    MAX_FILES ][ 4 ]; /* Of each zf METCRO3D file. */
  int wwindFileTimeRange[ MAX_FILES ][ 4 ]; /* Of each wwind METCRO3D file. */
  int yyyymmddhh;        /* First timestamp of available output data. */
  int readTimesteps;     /* Number of subset timesteps to read. */
  int skipFileCount;     /* Number of first input files to skip. */
  int isHourlyTimesteps; /* Are input file timesteps hourly? */
  int outputTimesteps;   /* Number of output timesteps. */
  int layers;            /* Layers  in input file(s). */
  int rows;              /* Rows    in input file(s). */
  int columns;           /* Columns in input file(s). */
  int isProjected;       /* Is grid in projected space? I.e., not lon-lat. */
  double* longitudes;    /* Grid cell edge longitudes[ rows + 1][columns + 1]*/
  double* latitudes;     /* Grid cell edge latitudes[  rows + 1][columns + 1]*/
  double* elevations;    /* " cell center elevations[layers][rows][columns] */
  float* heights;        /* " " " height in meters above mean sea level [r*c]*/
  char variableUnits[ MXVARS3 ][ NAMLEN3 + 1 ]; /* Output variable units. */
  char variableDescriptions[ MXVARS3 ][ MXDLEN3 + 1 ]; /* Output var desc. */
  char* wwindVariable;   /* Input either "WWIND" or "W_VEL". Output "WWIND". */
} Data;

/* Aggregation options: */

enum {
  AGGREGATE_NONE,
  AGGREGATE_DAILY_MEAN,
  AGGREGATE_DAILY_MAX,
  AGGREGATE_DAILY_MAX8,
  AGGREGATE_MEAN,
  AGGREGATE_SUM,
  AGGREGATE_MODES
};

#define AGGREGATE_STRING "none daily_mean daily_max daily_max8 mean sum"


/* Aux modes: */

enum {
  INTEGRATE = 1, WIND,
  VERSION, LIST, PRINT_WORKING_DIRECTORY, DIRECTORY_LISTING
};


/*
 * Implemented output formats are listed below.
 * To implement a new format extend each entry below:
 *   1. Add FORMAT_new in enum.
 *   2. Extend FORMAT_STRING.
 *   3. Add new forward declaration of writer routine.
 *   4. Add new writer routine to writers array.
 *   5. Add new writer routine in this file below.
 */

enum {
  FORMAT_XDR,
  FORMAT_ASCII,
  FORMAT_COARDS,
  FORMAT_IOAPI,
  FORMATS
};

#define FORMAT_STRING "xdr ascii coards ioapi"

static int writeXDR(    Data* const data );
static int writeASCII(  Data* const data );
static int writeCOARDS( Data* const data );
static int writeIOAPI(  Data* const data );

typedef int (*Writer)( Data* );

static const Writer writers[ FORMATS ] = {
  writeXDR,
  writeASCII,
  writeCOARDS,
  writeIOAPI
};


/* This table specifies changes to some variable names/units/var_desc output:*/

/*
 * The edit global variable is set from command-line option -edit
 * to edit output variable name/units/descriptions of CMAQ EQUATES files.
 * edit = 1 means search/edit, edit = 0 means do not search/edit.
 */

static int edit = 0;

typedef struct {
  const char* const name;        /* Variable name in input file. */
  const char* const newName;     /* Variable name in output file. */
  const char* const units;       /* Variable units. */
  const char* const description; /* Variable description. */
} VariableMetadata;

static const VariableMetadata variableMetadata[] = {

  /* 2021-06-16: COMBINE_DEP_var_names_units_desc.xlsx: */

  { "RT",            0,               "cm",       "Precipitation" },
  { "DD_OXN_NOX",    0,               "kgN ha-1", "Dry Deposition of NOX (NO, NO2)" },
  { "WD_OXN_NOX",    0,               "kgN ha-1", "Wet Deposition of NOX (NO, NO2)" },
  { "DD_OXN_TNO3",   0,               "kgN ha-1", "Dry Deposition of Total Nitrate (HNO3, NO3)" },
  { "WD_OXN_TNO3",   0,               "kgN ha-1", "Wet Deposition of Total Nitrate (HNO3, NO3)" },
  { "DD_OXN_PANT",   0,               "kgN ha-1", "Dry Depostion of PANs (PAN, PANX, OPAN)" },
  { "WD_OXN_PANT",   0,               "kgN ha-1", "Wet Depostion of PANs (PAN, PANX, OPAN)" },
  { "DD_OXN_ORGN",   0,               "kgN ha-1", "Dry Deposition of Organic N (NTR1, NTR2, INTR)" },
  { "WD_OXN_ORGN",   0,               "kgN ha-1", "Wet Deposition of Organic N (NTR1, NTR2, INTR)" },
  { "DD_OXN_OTHR",   0,               "kgN ha-1", "Dry Deposition of Other Oxidized N (N2O5, HONO)" },
  { "WD_OXN_OTHR",   0,               "kgN ha-1", "Dry Deposition of Other Oxidized N (N2O5, HONO, PNA)" },
  { "DD_OXN_TOT",    "DRYDEP_OXN",    "kgN ha-1", "Dry Deposition of Oxidized Nitrogen (NOX, TNO3, PANs, Org N, N205, HONO, PNA)" },
  { "WD_OXN_TOT",    "WETDEP_OXN",    "kgN ha-1", "Wet Deposition of Oxidized Nitrogen (NOX, TNO3, PANs, Org N, N205, HONO, PNA)" },
  { "TD_OXN_TOT",    "TOTDEP_OXN",    "kgN ha-1", "Total (Dry + Wet) Deposition of Oxidized Nitrogen" },
  { "DD_REDN_TOT",   "DRYDEP_REDN",   "kgN ha-1", "Dry Deposition of Reduced Nitorgen (NH4, NH3)" },
  { "WD_REDN_TOT",   "WETDEP_REDN",   "kgN ha-1", "Wet Deposition of Reduced Nitorgen (NH4, NH3)" },
  { "TD_REDN_TOT",   "TOTDEP_REDN",   "kgN ha-1", "Total (Dry + Wet) Deposition of Reduced Nitrogen" },
  { "DD_N_TOT",      "DRYDEP_N",      "kgN ha-1", "Dry Deposition of Nitrogen" },
  { "WD_N_TOT",      "WETDEP_N",      "kgN ha-1", "Wet Deposition of Nitrogen" },
  { "TD_N_TOT",      "TOTDEP_N",      "kgN ha-1", "Total (Dry + Wet) Deposition of Nitrogen" },
  { "DD_S_TOT",      "DRYDEP_S",      "kgS ha-1", "Dry Deposition of Sulfur" },
  { "WD_S_TOT",      "WETDEP_S",      "kgS ha-1", "Wet Deposition of Sulfur" },
  { "TD_S_TOT",      "TOTDEP_S",      "kgS ha-1", "Total (Dry + Wet) Deposition of Sulfur" },

  /* 2021-06-16: COMBINE_ACONC_var_names_units_desc.xlsx: */

  { "AALJ",          "PMF_AL",        "ug m-3",   "Fine Particle Aluminum" },
  { "ACAJ",          "PMF_CA",        "ug m-3",   "Fine Particle Calcium" },
  { "ACAK",          "PMC_CA",        "ug m-3",   "Coarse Particle Calcium" },
  { "ACLIJ",         "PMF_CL",        "ug m-3",   "Fine Particle Chloride" },
  { "AECIJ",         "PMF_EC",        "ug m-3",   "Fine Particle Elemental Carbon" },
  { "AFEJ",          0,               "ug m-3",   "Fine Particle Iron	PMF_FE" },
  { "AHPLUSIJ",      "PMF_HPLUS",     "umol m-3", "Fine Particle Hydronium Ion" },
  { "AIR_DENS",      0,               "kg m-3",   "Air Density" },
  { "AKJ",           "PMF_K",         "ug m-3",   "Fine Particle Potassium" },
  { "AKK",           "PMC_K",         "ug m-3",   "Coarse Particle Potassium" },
  { "ALD2",          0,               "ppbV",     "Acetaldehyde" },
  { "AMGK",          "PMC_MG",        "ug m-3",   "Coarse Particle Magnesium" },
  { "AMGJ",          "PMF_MG",        "ug m-3",   "Fine Particle Magnesium" },
  { "AMNJ",          "PMF_MN",        "ug m-3",   "Fine Particle Manganese" },
  { "ANAIJ",         "PMF_NA",        "ug m-3",   "Fine Particle Sodium" },
  { "ANAK",          "PMC_NA",        "ug m-3",   "Coarse Particle Sodium" },
  { "ANCOMIJ",       "PMF_NCOM",      "ug m-3",   "Fine Particle Non-Carbon Organic Mass (OM - OC)" },
  { "ANH4IJ",        "PMF_NH4",       "ug m-3",   "Fine Particle Ammonium" },
  { "ANH4K",         "PMC_NH4",       "ug m-3",   "Coarse Particle Ammonium" },
  { "ANO3K",         "PMC_NO3",       "ug m-3",   "Coarse Particle Nitrate" },
  { "ANO3IJ",        "PMF_NO3",       "ug m-3",   "Fine Particle Nitrate" },
  { "ANO3_PPB",      "PMF_NO3_PPB",   "ppbV",     "Fine Particle Nitrate (mixing ratio)" },
  { "AOCIJ",         "PMF_OC",        "ugC m-3",  "Fine Particle Organic Carbon (C only)" },
  { "AOMIJ",         "PMF_OM",        "ug m-3",   "Fine Particle Organic Matter (C,H,O,N, etc)" },
  { "AOMOCRAT_TOT",  "PMF_OMOC",      "ug ug-1",  "Fine Particle OM/OC Ratio" },
  { "AORGCJ",        "PMF_CLDGLY",    "ug m-3",   "Glyoxal and methylglyoxal SOA produced in cloud water" },
  { "APOCIJ",        "PMF_POC",       "ugC m-3",  "Fine Particle Primary Organic Carbon" },
  { "APOMIJ",        "PMF_POA",       "ug m-3",   "Fine Particle Primary Organic Matter" },
  { "ASIJ",          "PMF_SI",        "ug m-3",   "Fine Particle Silicon" },
  { "ASO4K",         "PMC_SO4",       "ug m-3",   "Coarse Particle Sulfate" },
  { "ASO4IJ",        "PMF_SO4",       "ug m-3",   "Fine Particle Sulfate" },
  { "ASOCIJ",        "PMF_SOC",       "ugC m-3",  "Fine Particle Secondary Organic Carbon" },
  { "ASOILJ",        "PMF_SOIL_IMPV", "ug m-3",   "Fine Particle Lumped Crustal Species calculated with IMPROVE method" },
  { "ASOMIJ",        "PMF_SOA",       "ug m-3",   "Fine Particle Secondary Organic Matter" },
  { "ATIJ",          "PMF_TI",        "ug m-3",   "Fine Particle Titanium" },
  { "ATOTIJK",       "PM_MASS",       "ug m-3",   "Total Particle Mass" },
  { "ATOTI",         "PMAIT_MASS",    "ug m-3",   "Aitken Particle Mass" },
  { "ATOTK",         "PMC_MASS",      "ug m-3",   "Coarse Particle Mass" },
  { "ATOTJ",         "PMACC_MASS",    "ug m-3",   "Accumulation Particle Mass" },
  { "ATOTIJ",        "PMF_MASS",      "ug m-3",   "Fine Particle Mass" },
  { "BENZENE",       0,               "ppbV",     "Benzene" },
  { "CO",            0,               "ppbV",     "Carbon Monoxide" },
  { "ETH",           0,               "ppbV",     "Ethene" },
  { "ETHA",          0,               "ppbV",     "Ethane" },
  { "FORM",          0,               "ppbV",     "Formaldehyde" },
  { "H2O2",          0,               "ppbV",     "Hydrogen Peroxide" },
  { "HNO3",          0,               "ppbV",     "Nitric Acid" },
  { "HNO3_UGM3",     0,               "ug m-3",   "Nitric Acid (concentration)" },
  { "HONO",          0,               "ppbV",     "Nitrous Acid" },
  { "HOX",           0,               "ppbV",     "Hydroxyl Radical (OH) + Hydroperoxy Radical (HO2)" },
  { "ISOP",          0,               "ppbV",     "Isoprene" },
  { "N2O5",          0,               "ppbV",     "Dinitrogen Pentoxide" },
  { "NH3",           0,               "ppbV ",    "Ammonia" },
  { "NH3_UGM3",      0,               "ug m-3",   "Ammonia (concentration)" },
  { "NHX",           0,               "ug m-3",   "Inorganic Nitrogen (ammonia gas plus particulate ammonium)" },
  { "NO",            0,               "ppbV",     "Nitric Oxide" },
  { "NO2",           0,               "ppbV",     "Nitrogen Dioxide" },
  { "NOX",           0,               "ppbV",     "Nitrogen Oxides (NO + NO2)" },
  { "NOY",           0,               "ppbV",     "Total Reative Nitrogen (NO + NO2 + HNO3 + PAN + other organic nitrates)" },
  { "NTR",           0,               "ppbV",     "Monofunctional Organic Nitrates (NTR1) + Multifunctional Organic Nitrates (NTR2) +Nitrate from Isoprene (INTR)" },
  { "O3",            0,               "ppbV",     "Ozone" },
  { "OH",            0,               "ppbV",     "Hydroxyl Radical" },
  { "PANS",          0,               "ppbV",     "Peroxyacylnitrate (PAN) + peroxyacylnitrates with 3 or morecarbons (PANX) + peroxyacylnitrate from OPO3 (OPAN)" },
  { "PBLH",          0,               "m",        "Planetary Boundary Layer Height" },
  { "PM1_TOT",       "PM1",           "ug m-3",   "PM1 (sharp 1 micrometer cutoff computed using modeled size distribution)" },
  { "PM10",          0,               "ug m-3",   "Particulate Matter up to 10 micrometers in Size" },
  { "PM25_CA",       0,               "ug m-3",   "PM2.5 Calcium  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_CL",       0,               "ug m-3",   "PM2.5 Chloride  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_EC",       0,               "ug m-3",   "PM2.5 Elemental Carbon  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_FRM",      0,               "ug m-3",   "FRM Equivalent PM2.5 (computed using modeled size distribution)" },
  { "PM25_HP",       0,               "ug m-3",   "Hydronium Ion (sharp cutoff computed using modeled size distribution)" },
  { "PM25_K",        0,               "ug m-3",   "PM2.5 Potassium  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_MG",       0,               "ug m-3",   "PM2.5 Magnesium  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_NA",       0,               "ug m-3",   "PM2.5 Sodium  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_NH4",      0,               "ug m-3",   "PM2.5 Ammonium  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_NO3",      0,               "ug m-3",   "PM2.5 Nitrate  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_OC",       0,               "ugC m-3",  "PM2.5 Organic Carbon  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_OM",       0,               "ug m-3",   "PM2.5 Organic Matter  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_SO4",      0,               "ug m-3",   "PM2.5 Sulfate  (sharp cutoff computed using modeled size distribution)" },
  { "PM25_SOIL",     0,               "ug m-3",   "PM2.5 Lumped Crustal Species (sharp cutoff computed using modeled size distribution)" },
  { "PM25_TOT",      "PM25",          "ug m-3",   "Total PM2.5 (sharp cutoff computed using modeled size distribution)" },
  { "PM25_UNSPEC1",  0,               "ug m-3",   "Other PM2.5 Species (Total - (CL+EC+NA+NH4+NO3+OC+SOIL+SO4))	" },
  { "PMC_CL",       "PM25to10_CL",    "ug m-3",   "Coarse Mode Chlorine (Total CL - PM25_CL)" },
  { "PMC_NA",       "PM25to10_NA",    "ug m-3",   "Coarse Mode Sodium (Total NA - PM25_NA)" },
  { "PMC_NH4",      "PM25to10_NH4",   "ug m-3",   "Coarse Mode Ammonium (Total NH4 - PM25_NH4)" },
  { "PMC_NO3",      "PM25to10_NO3",   "ug m-3",   "Coarse Mode Nitrate (Total Particle NO3 - PM25_NO3)" },
  { "PMC_SO4",      "PM25to10_SO4",   "ug m-3",   "Coarse Mode Sulfate (Total Particle SO4 - PM25_SO4)" },
  { "PMC_TOT",      "PM25to10",       "ug m-3",   "Coarse Mode Particulate Matter (Total PM - PM25_TOT)" },
  { "PMIJ_FRM",     "PMF_FRM",        "ug m-3",   "FRM Equivalent Particulate Matter (Fine Mode)" },
  { "precip",       0,                "cm",       "Precipitation" },
  { "RH",           0,                "%",        "Relative Humidity" },
  { "SFC_TMP",      0,                "C",        "Surface Temperature" },
  { "SO2",          0,                "ppbV",     "Sulfur Dioxide" },
  { "SO2_UGM3",     0,                "ug m-3",   "Sulfur Dioxide (concentration)" },
  { "SOL_RAD",      0,                "W m-2",    "Solar Radiation" },
  { "TERP",         0,                "ppbV",     "Monoterpenes" },
  { "TNO3",         0,                "ug m-3",   "Total Nitrate" },
  { "TOL",          0,                "ppbV",     "Toluene and Other Monoalkyl Aromatics" },
  { "WDIR10",       0,                "deg",      "10-m Wind Speed" },
  { "WSPD10",       0,                "m s-1",    "10-m Wind Direction" },
  { "XYL",          0,                "ppbV",     "Xylene and Other Polyalkyl Aromatics except Naphthalene" },

  /* Replace incorrect units in MET files: */

  { "QC",           0,                "kg kg-1",  "Cloud water mixing ratio" },
  { "QR",           0,                "kg kg-1",  "Rain water mixing ratio" },
  { "QI",           0,                "kg kg-1",  "Ice mixing ratio" },
  { "QS",           0,                "kg kg-1",  "Snow mixing ratio" },
  { "QG",           0,                "kg kg-1",  "Graupel mixing ratio" },

  /* Edit variable names of LST files : */

  { "O3_MDA8",      "DAILY_O3MAX8",   "ppbV",  "Local daily 8-hour maximum ozone." },
  { "O3_AVG",       "DAILY_O3",       "ppbV",  "Local daily average ozone." },
  { "CO_AVG",       "DAILY_CO",       "ppbV",  "Local daily average carbon monoxide." },
  { "NO_AVG",       "DAILY_NO",       "ppbV",  "Local daily average nitrogen oxide." },
  { "NO2_AVG",      "DAILY_NO2",      "ppbV",  "Local daily average nitrogen dioxide." },
  { "SO2_AVG",      "DAILY_SO2",      "ppbV",  "Local daily average sulfur dioxide." },
  { "CH2O_AVG",     "DAILY_CH2O",     "ppbV",  "Local daily average formaldehyde." },
  { "PM10_AVG",     "DAILY_PM10",     "ug/m3", "Local daily average particulate matter up to 10 micrometers in size." },
  { "PM25_AVG",     "DAILY_PM25",     "ug/m3", "Local daily average particulate matter up to 2.5 micrometers in size." },
  { "PM25_SO4_AVG", "DAILY_PM25_SO4", "ug/m3", "Local daily average sulfate particulate matter up to 2.5 micrometers in size." },
  { "PM25_NO3_AVG", "DAILY_PM25_NO3", "ug/m3", "Local daily average nitrate particulate matter up to 2.5 micrometers in size." },
  { "PM25_NH4_AVG", "DAILY_PM25_NH4", "ug/m3", "Local daily average ammonium particulate matter up to 2.5 micrometers in size." },
  { "PM25_OC_AVG",  "DAILY_PM25_CO",  "ug/m3", "Local daily average organic carbon particulate matter up to 2.5 micrometers in size." },
  { "PM25_EC_AVG",  "DAILY_PM25_EC",  "ug/m3", "Local daily average elemental carbon particulate matter up to 2.5 micrometers in size." },

  { 0, 0, 0, 0 } /* End of table. */
};

/*========================== FORWARD DECLARATIONS ===========================*/

static const char* outputVariableName( const char* const variableName );

static const char* outputVariableUnits( const char* const variableName,
                                        const char* const variableUnits,
                                        const int convertSpaces );

static const char* outputVariableDescription( const char* const variableName,
                                              const char* const
                                                variableDescription );

static void printUsage( const char* programName );

static int parseCommandLineOptions( int argc, char* argv[],
                                    Arguments* arguments );

static int isValidArguments( const Arguments* const arguments );
static int initializeData( Data* const data );
static int checkOrSetTimeSubset( Data* const data );
static int checkOrSetVariables( Data* const data );
static int checkOrSetGridSubset( Data* const data );
static int checkFilesAreCompatible( Data* const data );
static int checkFileVariables( const Data* const data );
static int checkHTFile( const Data* const data );
static int checkZFFiles( Data* const data );
static int checkWWINDFiles( Data* const data );
static int computeGridCellCoordinates( Data* const data );
static int computeGridCellCenterElevations( Data* const data );

static void checkAndFixVerticalGridParameters( const int layers,
                                               int* vgtyp,
                                               float* vgtop,
                                               float vglvls[] );

static int readHT( Data* const data );
static int boundsSubset( Data* const data );

static void computeCellBounds( const Data* const data,
                               const int row,
                               const int column,
                               Bounds cellBounds );

#if TEST_CLIP_GRID_CELLS
static void getCellVertices( const Data* const data,
                             const int row, const int column,
                             double longitudes[ 4 ], double latitudes[ 4 ] );
#endif

static void computeBounds( const size_t points,
                           const double longitudes[],
                           const double latitudes[],
                           Bounds bounds );

static int findTimestampedVariable( const Data* const data,
                                    const char* const variableName,
                                    const int yyyymmddhh,
                                    int* const variableId,
                                    int* const timestep );

static int writeXDRData( Data* const data, FILE* const output );
static int writeXDRHeader( const Data* const data, FILE* const output );
static int writeXDRProjection( const Data* const data, FILE* const output );
static int writeXDRGrid( const Data* const data, FILE* const output );
static int writeXDRVariableNamesAndUnits( const Data* const data,
                                          FILE* const output );

static int readSubset( Data* const data,
                       const char* const variable,
                       const int yyyymmddhh,
                       float subsetData[],
                       float subsetZF[],
                       float subsetDENS[] );

static int readSubsetVariable( Data* const data,
                               const char* const variableName,
                               const int yyyymmddhh,
                               float subsetData[] );

static int readSubsetZH( Data* const data, const int yyyymmddhh,
                         float subsetData[] );

static int readSubsetWWIND( Data* const data, const int yyyymmddhh,
                            float subsetData[] );

static int readSubsetZFAndDENS( Data* const data,
                                const int yyyymmddhh,
                                float subsetZF[],
                                float subsetDENS[] );

static void copySubsetCoordinates( const Data* const data,
                                   const char* const variableName,
                                   float subsetData[] );

static void copySubsetElevations( const Data* const data, float subsetData[] );

static int expandSubsetData( const int layers, const int rows,
                             const int columns, const int expandRow,
                             const int expandColumn, float subsetData[] );

static int createTemporaryXDRFile( Data* const data );

static int readXDRHeader( FILE* file,
                          int* timesteps,
                          int* variables,
                          int* layers,
                          int* rows,
                          int* columns,
                          int* yyyymmddhh,
                          char variableNames[ MXVARS3 ][ NAMLEN3 + 1 ],
                          char variableUnits[ MXVARS3 ][ NAMLEN3 + 1 ] );

static int createTemporaryCOARDSFileHeader( Data* const data );

static int writeCOARDSData( Data* const data, const int file );

static int createCOARDSStandardVariables( Data* const data,
                                          const int file,
                                          const int dimids[ 4 ] );

static int createTemporaryIOAPIFileHeader( Data* const data );

static int writeTFLAG( Data* const data, const int file );

static int writeIOAPIData( Data* const data, const int file );

static void aggregateData( const int mode,
                           const size_t timesteps, const size_t cells,
                           float data[], float means[], int counts[] ) ;

static void aggregateMean( const size_t timesteps, const size_t cells,
                           float data[] );

static void aggregateMax( const size_t timesteps, const size_t cells,
                          float data[] );

static void aggregateMax8( const size_t timesteps, const size_t cells,
                           float data[] );

static void aggregateAll( const size_t timesteps, const size_t cells,
                          const float data[], float means[], int counts[] );

static void aggregateSum( const size_t timesteps, const size_t cells,
                          const float data[], float means[] );

static void computeZ( const double g,
                      const double R,
                      const double A,
                      const double T0s,
                      const double P00,
                      const int layers,
                      const int type,
                      const double topPressure,
                      const double heightOfTerrainInMeters,
                      const float levels[],
                      double z[] );

static double pressureAtSigmaLevel( const double sigmaLevel,
                                    const double pressureAtTop );

static double heightAtPressure( double pressure );

static void elevationsAtSigmaPressures( const double g,
                                        const double R,
                                        const double A,
                                        const double T0s,
                                        const double P00,
                                        const double surfaceElevation,
                                        const int levels,
                                        const double topPressure,
                                        const float sigmaPressures[],
                                        double elevations[] );

static void integrateLayers( const size_t timesteps,
                             const size_t layers,
                             const size_t rows,
                             const size_t columns,
                             const float* const zf,
                             const float* const dens,
                             float* const data );

#ifdef DEBUGGING
static int printArgs( int argc, char* argv[] ) {
  int index = 0;

  while ( index < argc ) {
    fprintf( stderr, "%s ", argv[ index ] );
    ++index;
  }

  fputc( '\n', stderr );
  return argc;
}
#endif

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: main - Parse options and write subset data to stdout.
INPUTS:  int argc      Number of command-line arguments.
         char* argv[]  Command-line argument strings.
RETURNS: 0 if successful, else 1.
******************************************************************************/

int main( int argc, char* argv[] ) {
  Arguments arguments;
  DEBUG( int unused_ = printArgs( argc, argv ); )
  int ok = 0;
  memset( &arguments, 0, sizeof arguments );
  ok = parseCommandLineOptions( argc, argv, &arguments );

  if ( ! ok ) {
    const char* const programName = AND2( argv, *argv ) ? argv[ 0 ] : "";
    printUsage( programName );
  } else if ( arguments.auxMode == VERSION ) {
    ok = printf( "%d\n", fileDateUTC( argv[ 0 ] ) ) > 0;
  } else if ( arguments.auxMode == LIST ) {
    ok = printM3IOVariables( arguments.fileNames[ 0 ] );
  } else if ( arguments.auxMode == PRINT_WORKING_DIRECTORY ) {
    ok = printWorkingDirectory();
  } else if ( arguments.auxMode == DIRECTORY_LISTING ) {
    ok = printDirectoryListing( arguments.lsDir );
  } else {
    Data data;
    memset( &data, 0, sizeof data );
    data.arguments = &arguments;
    ok = initializeData( &data );

    if ( ok ) {
      Writer writer = 0;
      CHECK2( IN_RANGE( arguments.format, 0, FORMATS - 1 ),
              IN_RANGE( arguments.format, 0,
                        sizeof writers / sizeof *writers - 1 ) );
      writer = writers[ arguments.format ];
      CHECK( writer );
      ok = writer( &data );
      FREE( data.longitudes );
      FREE( data.latitudes );
      FREE( data.elevations );
      FREE( data.heights );
    }
  }

  DEBUG( fprintf( stderr, "main returning %d\n", ! ok ); )
  return ! ok;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: printUsage - Print usage instructions for the program.
INPUTS:  const char* programName  The name of the program.
******************************************************************************/

static void printUsage( const char* programName ) {
  PRE0( programName );
  fprintf( stderr, "\n\n\n%s - "
           "Read CMAQ (NetCDF M3IO/IOAPI FORMAT) grid files\n",
           programName );
  fprintf( stderr, "and write the specified subset of data variables\n" );
  fprintf( stderr, "to stdout in XDR, ASCII, COARDS or IOAPI format.\n" );
  fprintf( stderr, "\nUsage:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "[-pwd] (Print current working directory.)\\\n" );
  fprintf( stderr, "[-ls directory] (Print subdirectories & NetCDF files)\\\n");
  fprintf( stderr, "-files <file> [<file> ...] \\\n" );
  fprintf( stderr, "[-tmpdir directory] (Default is .)\\\n" );
  fprintf( stderr, "[-output file] (Default is stdout)\\\n" );
  fprintf( stderr, "[-desc 'description text'] \\\n" );
  fprintf( stderr, "[-format xdr | ascii | coards | ioapi] (Default is ioapi.)\\\n" );
  fprintf( stderr, "[-ht <gridcro2d> ] \\\n" );
  fprintf( stderr, "[-zf <metcro3d> [<metcro3d> ...]] \\\n" );
  fprintf( stderr, "[-wwind <metcro3d> [<metcro3d> ...]] \\\n" );
  fprintf( stderr, "[-integrate_layers]\\\n" );
  fprintf( stderr, "[-lonlat] \\\n" );
  fprintf( stderr, "[-elevation] \\\n" );
  fprintf( stderr, "[-edit] \\\n" );
  fprintf( stderr, "[-list] \\\n" );
  fprintf( stderr, "[-ellipsoid <major_semiaxis> <minor_semiaxis>] \\\n" );
  fprintf( stderr, "[-variable <name> ...] \\\n" );
  fprintf( stderr, "[-time   <yyyymmddhh1> [<yyyymmddhh2>]] \\\n" );
  fprintf( stderr, "[-layer  <first> [<last>]] \\\n" );
  fprintf( stderr, "[-row    <first> [<last>]] \\\n" );
  fprintf( stderr, "[-column <first> [<last>]] \\\n" );
  fprintf( stderr, "[-bounds <minimum_longitude> <minimum_latitude>" );
  fprintf( stderr, " <maximum_longitude> <maximum_latitude> ] \\\n" );
  fprintf( stderr, "[-aggregate daily_mean | daily_max | daily_max8 | mean | sum] \\\n");
  fprintf( stderr, "Note: -layer/row/column are 1-based.\n");
  fprintf( stderr, "If -bounds option is used then -row/-column options" );
  fprintf( stderr, " cannot be used.\n");
  fprintf( stderr, "-lonlat adds variables LONGITUDE LATITUDE in output.\n");
  fprintf( stderr, "-elevation adds variable ELEVATION in output.\n");
  fprintf( stderr, "-edit enables rename of variable names/units/descriptions\n");
  fprintf( stderr, "of CMAQ EQUATES files (using a built-in table).\n");
  fprintf( stderr, "-list lists variable names\n" );
  fprintf( stderr, "-integrate_layers option integrates the given variable ");
  fprintf( stderr, "(with units ppmV or ppbV) over the layers.\n");
  fprintf( stderr, "-wwind option specifies the METCRO3D files containing ");
  fprintf( stderr, "the variable WWIND or CCTM_CONC files containing W_VEL.\n");
  fprintf( stderr, "-aggregate daily_max8 is daily max of 17 8-hour means.\n");
  fprintf( stderr, "Default ellipsoid is a sphere of radius 6,370,000m.\n");
  fprintf( stderr, "Use ncdump -h to list variables in each file.\n");
  fprintf( stderr, "\n\n\n--------------------------------------------\n\n" );

  fprintf( stderr, "Example #1:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc \"https://www.epa.gov/cmaq,CMAQSubset\" \\\n" );
  fprintf( stderr, "-format xdr \\\n" );
  fprintf( stderr, "-ellipsoid 6370000 6370000 \\\n" );
  fprintf( stderr, "-files CCTM_2pCONC.20020423 CCTM_2pCONC.20020424 \\\n");
  fprintf( stderr, "-ht GRIDCRO2D_020423 \\\n");
  fprintf( stderr, "-zf METCRO3D_020423 METCRO3D_020424 \\\n");
  fprintf( stderr, "-variable NO2 O3 \\\n");
  fprintf( stderr, "-time 2002042300 2002042423 \\\n");
  fprintf( stderr, "-layer 1 5 \\\n");
  fprintf( stderr, "-bounds -123 24 -10 30 \\\n");
  fprintf( stderr, "> subset.xdr \n");
  fprintf( stderr, "\n" );
  fprintf( stderr, "Outputs an ASCII header followed by binary\n");
  fprintf( stderr,"array data[variables][timesteps][layers][rows][columns]\n");
  fprintf( stderr, "\n" );
  fprintf( stderr, "SUBSET 9.0 CMAQ\n" );
  fprintf( stderr, "M_02_99BRACE\n" );
  fprintf( stderr, "https://www.epa.gov/cmaq,CMAQSubset\n" );
  fprintf( stderr, "2000-04-23T00:00:00-0000\n" );
  fprintf( stderr, "# data dimensions: " );
  fprintf( stderr, "timesteps variables layers rows columns:\n" );
  fprintf( stderr, "48 4 3 65 83\n" );
  fprintf( stderr, "# subset indices (0-based time, " );
  fprintf( stderr, "1-based layer/row/column): " );
  fprintf( stderr, "first-timestep last-timestep first-layer last-layer");
  fprintf( stderr, " first-row last-row first-column last-column:\n" );
  fprintf( stderr, "0 47 1 5 1 65 1 83\n" );
  fprintf( stderr, "# Variable names:\n" );
  fprintf( stderr, "LONGITUDE LATITUDE ELEVATION NO2 O3\n" );
  fprintf( stderr, "# Variable units:\n" );
  fprintf( stderr, "deg deg m ppmV ppmV\n" );
  fprintf( stderr, "# lcc projection: lat_1 lat_2 lat_0 lon_0 " );
  fprintf( stderr, "major_semiaxis minor_semiaxis\n" );
  fprintf( stderr, "30 60 40 -100 6.36747e+06 6.36747e+06\n" );
  fprintf( stderr,
      "# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[22]:\n");
  fprintf( stderr,
          "268 259 1.578e+06 -1.27e+06 2000 2000 2 10000 1 0.995 0.99 " );
  fprintf( stderr, "0.985 0.98 0.97 0.96 0.945 0.93 0.91 0.89 0.865 0.84 " );
  fprintf( stderr, "0.78 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0\n" );
  fprintf( stderr, "# IEEE-754 32-bit doubles data[variables][timesteps]" );
  fprintf( stderr, "[layers][rows][columns]:\n<binary data array here>\n\n\n");

  fprintf( stderr, "Example #2:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc \"https://www.epa.gov/cmaq,CMAQSubset\" \\\n" );
  fprintf( stderr, "-format ioapi \\\n" );
  fprintf( stderr, "-ellipsoid 6370000 6370000 \\\n" );
  fprintf( stderr, "-files CCTM_2pCONC.20020423 CCTM_2pCONC.20020424 \\\n");
  fprintf( stderr, "-time 2002042300 2002042423 \\\n");
  fprintf( stderr, "-variable O3 \\\n");
  fprintf( stderr, "> subset.ncf\n");
  fprintf( stderr, "\nOutputs a subset in IOAPI format which is\n" );
  fprintf( stderr, "redirected to a local file 'subset.ncf'.\n");
  fprintf( stderr, "The file may be viewed with ncdump subset.ncf | more.\n\n");

  fprintf( stderr, "Example #3:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc \"https://www.epa.gov/cmaq,CMAQSubset\" \\\n" );
  fprintf( stderr, "-format ascii \\\n" );
  fprintf( stderr, "-lonlat \\\n" );
  fprintf( stderr, "-files CCTM_2pCONC.20020423 CCTM_2pCONC.20020424 \\\n");
  fprintf( stderr, "-ht GRIDCRO2D_020423 \\\n");
  fprintf( stderr, "-zf METCRO3D_020423 METCRO3D_020424 \\\n");
  fprintf( stderr, "-time 2002042300 2002042423 \\\n");
  fprintf( stderr, "-variable O3 \\\n");
  fprintf( stderr, "-aggregate daily_max8 \\\n");
  fprintf( stderr, "-layer 1 \\\n");
  fprintf( stderr, "-bounds -123 24 -120 30 \\\n");
  fprintf( stderr, "> subset.xdr\n");
  fprintf( stderr, "\n" );
  fprintf( stderr, "\nO3 daily 8-hour max data is written in a spreadsheet\n");
  fprintf( stderr, "importable format (tab-separated values).\n\n");

  fprintf( stderr, "Example #4:\n\n" );
  fprintf( stderr, "%s \\\n", programName );
  fprintf( stderr, "-tmpdir /data/tmp \\\n");
  fprintf( stderr, "-desc \"https://www.epa.gov/cmaq,CMAQSubset\" \\\n" );
  fprintf( stderr, "-format coards \\\n" );
  fprintf( stderr, "-ellipsoid 6370000 6370000 \\\n" );
  fprintf( stderr, "-files CCTM_2pCONC.20020423 CCTM_2pCONC.20020424 \\\n");
  fprintf( stderr, "-ht GRIDCRO2D_020423 \\\n");
  fprintf( stderr, "-zf METCRO3D_020423 METCRO3D_020424 \\\n");
  fprintf( stderr, "-integrate_layers \\\n");
  fprintf( stderr, "-variable CO \\\n");
  fprintf( stderr, "-time 2002042300 2002042423 \\\n");
  fprintf( stderr, "-bounds -90 30 -89 31 \\\n");
  fprintf( stderr, "> subset.nc\n");
  fprintf( stderr, "\n" );
  fprintf( stderr, "\nLayer-integrated CO data written in a spreadsheet\n");
  fprintf( stderr, "importable format (tab-separated values).\n");
  fprintf( stderr, "\n\n");
}



/******************************************************************************
PURPOSE: parseCommandLineOptions - Parse command line options.
INPUTS:  int argc              Number of command-line arguments.
         const char* argv[]    Command-line argument strings.
OUTPUTS: Arguments* arguments  User-specified arguments.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int parseCommandLineOptions( int argc, char* argv[],
                                    Arguments* arguments ) {

  PRE0( arguments );
  int result = 0;
  static const int one_int_max[ 2 ] = { 1, INT_MAX };
  static const double ellipsoidRange[ 2 ] =
    { ELLIPSOID_MINIMUM, ELLIPSOID_MAXIMUM };
  static Option options[] = {
    { "-tmpdir",             0, DIRECTORY_TYPE,      1, 0, 0, 0, 0 },
    { "-desc",               0, STRING_TYPE,         1, 0, 0, 0, 0 },
    { "-format",             0, ENUM_TYPE,           1, 0, FORMAT_STRING, 0,0},
    { "-ellipsoid",          0, REAL64_TYPE,         2, ellipsoidRange, 0,0,0},
    { "-files",              0, FILE_TYPE,           -MAX_FILES, 0, 0, 0, 0 },
    { "-ht",                 0, FILE_TYPE,           1, 0, 0, 0, 0 },
    { "-zf",                 0, FILE_TYPE,           -MAX_FILES, 0, 0, 0, 0 },
    { "-wwind",              0, FILE_TYPE,           -MAX_FILES, 0, 0, 0, 0 },
    { "-integrate_layers",   0, INT_TYPE,            0, 0, 0, 0, 0 },
    { "-variable",           0, STRING_TYPE,         -MXVARS3, 0, 0,0,0},
    { "-time",               0, YYYYMMDDHH_TYPE,     -2, 0, 0, 0, 0 },
    { "-layer",              0, INT_TYPE,            -2, one_int_max, 0, 0, 0 },
    { "-row",                0, INT_TYPE,            -2, one_int_max, 0, 0, 0 },
    { "-column",             0, INT_TYPE,            -2, one_int_max, 0, 0, 0 },
    { "-bounds",             0, BOUNDS_TYPE,         4, 0, 0, 0, 0 },
    { "-lonlat",             0, INT_TYPE,            0, 0, 0, 0, 0 },
    { "-elevation",          0, INT_TYPE,            0, 0, 0, 0, 0 },
    { "-aggregate",          0, ENUM_TYPE,           1, 0,AGGREGATE_STRING,0,0},
    { "-list",               0, INT_TYPE,            0, 0, 0, 0, 0 },
    { "-edit",               0, INT_TYPE,            0, 0, 0, 0, 0 },
    { "-output",             0, STRING_TYPE,         1, 0, 0, 0, 0 },

    /* These options are used to support remote file access: */

    { "-pwd",                0, INT_TYPE,            0, 0, 0, 0, 0 },
    { "-ls",                 0, DIRECTORY_TYPE,      1, 0, 0, 0, 0 },
    { "-version",            0, INT_TYPE,            0, 0, 0, 0, 0 }
  };
  memset( arguments, 0, sizeof (Arguments) );

  if ( AND6( argc > 0, argv, argv[ 0 ], argv[ 0 ][ 0 ],
             argv[ argc - 1 ], argv[ argc - 1 ][ 0 ] ) ) {
    const char* variableNames[ MXVARS3 ];
    memset( variableNames, 0, sizeof variableNames );

    /* Initialize arguments to defaults: */

    arguments->bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
    arguments->bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
    arguments->bounds[ LATITUDE  ][ MINIMUM ] =  -90.0;
    arguments->bounds[ LATITUDE  ][ MAXIMUM ] =   90.0;
    arguments->ellipsoid[ MINIMUM ] = arguments->ellipsoid[ MAXIMUM ] =
      DEFAULT_EARTH_RADIUS;

    /* Finish initializing non-compile-time-constant parts of options: */

    options[  0 ].values = &arguments->tmpDir;
    options[  1 ].values = &arguments->note;
    options[  2 ].values = &arguments->format;
    options[  3 ].values = &arguments->ellipsoid[ 0 ];
    options[  4 ].values = &arguments->fileNames[ 0 ];
    options[  5 ].values = &arguments->htFileName;
    options[  6 ].values = &arguments->zfFileNames[ 0 ];
    options[  7 ].values = &arguments->wwindFileNames[ 0 ];
    options[  8 ].values = 0;
    options[  9 ].values = &variableNames[ 0 ];
    options[ 10 ].values = &arguments->subset[ TIME ][ 0 ];
    options[ 11 ].values = &arguments->subset[ LAYER ][ 0 ];
    options[ 12 ].values = &arguments->subset[ ROW ][ 0 ];
    options[ 13 ].values = &arguments->subset[ COLUMN ][ 0 ];
    options[ 14 ].values = &arguments->bounds[ 0 ][ 0 ];
    options[ 15 ].values = 0;
    options[ 16 ].values = 0;
    options[ 17 ].values = &arguments->aggregateMode;
    options[ 18 ].values = 0;
    options[ 19 ].values = 0;
    options[ 20 ].values = &arguments->outputFileName;

    options[ 21 ].values = 0;
    options[ 22 ].values = &arguments->lsDir;
    options[ 23 ].values = 0;

    result =
      parseOptions( argc, argv, sizeof options / sizeof *options, options );

    {
      int hack = 0; /* -files is required unless -pwd, -ls, -version options.*/

      if ( options[ 21 ].parsed ) {
        arguments->auxMode = PRINT_WORKING_DIRECTORY;
        hack = 1;
      } else if ( options[ 22 ].parsed ) {
        arguments->auxMode = DIRECTORY_LISTING;
        hack = 1;
      } else if ( options[ 23 ].parsed ) {
        arguments->auxMode = VERSION;
        hack = 1;
      }

      if ( hack ) { /* Pretend -files hack was specified: */
        options[ 4 ].parsed = 2;
        arguments->fileCount = options[ 4 ].parsed - 1;
        arguments->fileNames[ 0 ] = "hack";
      } else {
        result = result && options[ 4 ].parsed >= 2; /* -files name required.*/
      }
    }

    if ( result ) {

      if ( arguments->tmpDir == 0 ) { /* If -tmpdir not specified. */
        arguments->tmpDir = defaultTemporaryDirectory; /* Set to default. */
      }

      if ( arguments->note == 0 ) { /* If -note not specified. */
        arguments->note = defaultNote; /* Set to defaultNote. */
      }

      if ( options[ 2 ].parsed == 0 ) { /* If -format is not specified. */
        arguments->format = FORMAT_IOAPI; /* Default to IOAPI. */
      }

      if ( options[ 18 ].parsed ) {
        arguments->auxMode = LIST;
      } else if ( options[ 8 ].parsed ) {
        arguments->auxMode = INTEGRATE;
      }

      edit = options[ 19 ].parsed;
      arguments->lonlat    = options[ 15 ].parsed;
      arguments->elevation = options[ 16 ].parsed;

      CHECK( options[ 4 ].parsed >= 2 ); /* -files file_name is required. */
      arguments->fileCount = options[ 4 ].parsed - 1;

      /* Get parsed string array counts: */

      arguments->zfFileCount =
        options[ 6 ].parsed > 1 ? options[ 6 ].parsed - 1 : 0;
      arguments->wwindFileCount =
        options[ 7 ].parsed > 1 ? options[ 7 ].parsed - 1 : 0;
      arguments->variables =
        options[ 9 ].parsed > 1 ? options[ 9 ].parsed - 1 : 0;

      if ( arguments->wwindFileCount > 0 ) {

        if ( arguments->auxMode ) {
          fputs( "\nCannot specify both -integrate_layers and -wwind.\n",
                 stderr );
        } else {
          arguments->auxMode = WIND;
        }
      }

      /* If only 1 (of 2 allowed) values was specified then copy to 2nd value:*/

      if ( options[ 10 ].parsed == 2 ) {
        arguments->subset[ TIME ][ MAXIMUM ] =
        arguments->subset[ TIME ][ MINIMUM ];
      }

      if ( options[ 11 ].parsed == 2 ) {
        arguments->subset[ LAYER ][ MAXIMUM ] =
        arguments->subset[ LAYER ][ MINIMUM ];
      }

      if ( options[ 12 ].parsed == 2 ) {
        arguments->subset[ ROW ][ MAXIMUM ] =
        arguments->subset[ ROW ][ MINIMUM ];
      }

      if ( options[ 13 ].parsed == 2 ) {
        arguments->subset[ COLUMN ][ MAXIMUM ] =
        arguments->subset[ COLUMN ][ MINIMUM ];
      }

      if ( arguments->variables ) { /* Copy variableNames to arguments: */
        const int variables = arguments->variables;
        int variable = 0;

        for ( ; variable < variables; ++variable ) {
          CHECK2( variableNames[ variable ], variableNames[ variable ][ 0 ] );
          strncpy( arguments->variableNames[ variable ],
                   variableNames[ variable ], NAMLEN3 );
        }
      }

      /* Integration requires zf files: */

      result =
        IMPLIES( arguments->auxMode == INTEGRATE,
                 AND3( IN_RANGE( arguments->zfFileCount, 1, MAX_FILES ),
                       arguments->zfFileNames[ 0 ],
                       arguments->zfFileNames[ arguments->zfFileCount - 1 ]));

      if ( ! result ) {
        fputs( "\n-integrate_layers requires -zf files.\n", stderr );
      } else {

        /* -bounds cannot be used with -row or -column: */

        result =
          ! AND2( OR2( arguments->subset[ ROW    ][ MAXIMUM ],
                       arguments->subset[ COLUMN ][ MAXIMUM ] ),
                  options[ 14 ].parsed );

        if ( ! result ) {
          fputs( "\nThe -bounds option cannot be used with "
                 "-row/-column options.\n", stderr );
        }
      }

      if ( AND2( result, arguments->outputFileName ) ) {

        /*
         * Redirect stdout to user-specified output file.
         * Note that the output file is opened for write-binary ("wb").
         * This prevents control-M characters from being inserted before
         * '\n' characters on Windows.
         */

        result = freopen( arguments->outputFileName, "wb", stdout ) != 0;
      }
    }
  }

  if ( ! result ) {
    fputs( "\nInvalid/insufficient command-line arguments.\n", stderr );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidArguments - Are arguments valid?
INPUTS:  const Arguments* const arguments  User-specified arguments to check.
RETURNS: 1 if valid, else 0.
******************************************************************************/

static int isValidArguments( const Arguments* const arguments ) {
  int result = arguments != 0;

  if ( result ) {
    result = AND2( result, IS_BOOL( arguments->lonlat ) );
    result = AND2( result, IS_BOOL( arguments->elevation ) );
    result = AND2( result, IN_RANGE( arguments->format, 0, FORMATS - 1 ) );
    result = AND2( result, IN4( arguments->auxMode, 0, INTEGRATE, WIND ) );
    result = AND2( result,
                   IN_RANGE( arguments->aggregateMode,
                             0, AGGREGATE_MODES - 1 ) );
    result = AND2( result, IN_RANGE( arguments->fileCount, 1, MAX_FILES ) );
    result = AND3( result, arguments->fileNames[ 0 ],
                   arguments->fileNames[ arguments->fileCount - 1 ] );
    result =
      AND3( result,
            IN_RANGE( arguments->zfFileCount,
                      arguments->auxMode == INTEGRATE, MAX_FILES ),
            IMPLIES( arguments->zfFileCount > 0,
                     AND2( arguments->zfFileNames[ 0 ],
                           arguments->zfFileNames[arguments->zfFileCount-1])));
    result =
      AND2( result,
           IMPLIES( arguments->auxMode == WIND,
                    AND3( IN_RANGE( arguments->wwindFileCount, 1, MAX_FILES ),
                          arguments->wwindFileNames[ 0 ],
                          arguments->wwindFileNames[
                                        arguments->wwindFileCount - 1 ] ) ) );
    result = AND4( result,
                   arguments->variables > 0,
                   arguments->variableNames[ 0 ][ 0 ],
                   arguments->variableNames[
                     arguments->variables - 1 ][ 0 ] );
    result =
      AND4( result,
           isValidYYYYMMDDHH( arguments->subset[ TIME ][ MINIMUM ] ),
           isValidYYYYMMDDHH( arguments->subset[ TIME ][ MAXIMUM ] ),
           arguments->subset[ TIME ][ MINIMUM ] <=
             arguments->subset[ TIME ][ MAXIMUM ] );
    result =
      AND3( result,
            IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ], 1, MXLAYS3 ),
            IN_RANGE( arguments->subset[ LAYER ][ MAXIMUM ],
                      arguments->subset[ LAYER ][ MINIMUM ], MXLAYS3 ) );
    result =
      AND2( result,
            IN_RANGE( arguments->subset[ ROW ][ MINIMUM ], 1,
                      arguments->subset[ ROW ][ MAXIMUM ] ) );
    result =
      AND2( result,
            IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ], 1,
                      arguments->subset[ COLUMN ][ MAXIMUM ] ) );

    result = AND3( result, arguments->note, arguments->note[ 0 ] );
    result = AND4( result, arguments->tmpDir, arguments->tmpDir[ 0 ],
                  isDirectory( arguments->tmpDir ) );
    result = AND2( result,
                   IN_RANGE( arguments->ellipsoid[ MINIMUM ],
                             ELLIPSOID_MINIMUM, ELLIPSOID_MAXIMUM ) );
    result = AND2( result,
                   IN_RANGE( arguments->ellipsoid[ MAXIMUM ],
                             arguments->ellipsoid[ MINIMUM ],
                             ELLIPSOID_MAXIMUM ) );
    result = AND2( result,
                   isValidBounds( (const double (*)[2]) arguments->bounds ) );
  }

  return result;
}



/******************************************************************************
PURPOSE: initializeData - Initialize data and check files/variables/units.
INPUTS:  Data* const data  Data to initialize/check.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int initializeData( Data* const data ) {
  PRE02( data, data->arguments );
  int result = checkOrSetTimeSubset( data );

  if ( result ) {
    result = checkOrSetVariables( data );

    if ( result ) {
      result = checkOrSetGridSubset( data );

      if ( result ) {
        result = isValidArguments( data->arguments );

        if ( ! result ) {
          fputs( "\nInvalid arguments.\n", stderr );
        } else {
          result = checkFilesAreCompatible( data );

          if ( result ) {
            result = computeGridCellCoordinates( data );

            if ( result ) {
              const Arguments* const arguments = data->arguments;

              /* If bounds are specified, reduce subset to bounds: */

              if ( ! AND4( arguments->bounds[ LONGITUDE ][ MINIMUM ] == -180.0,
                           arguments->bounds[ LONGITUDE ][ MAXIMUM ] ==  180.0,
                           arguments->bounds[ LATITUDE ][ MINIMUM] == -90.0,
                           arguments->bounds[ LATITUDE ][ MAXIMUM] ==  90.0)) {
                result = boundsSubset( data ); /* Reduce subset to bounds: */
              }

              if ( AND2( result, data->arguments->elevation ) ) {

                if ( data->arguments->htFileName ) {
                  result = readHT( data );
                }

                if ( result ) {
                  result = computeGridCellCenterElevations( data );
                }
              }
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: checkOrSetTimeSubset - Check time subset against the input or, if the
         time subset is unspecified then set to full range of input file.
INPUTS:  Data* const data  arguments->subset[TIME][MIN/MAXIMUM] to check/set.
OUTPUTS: Data* const data  subset[ TIME ][ MINIMUM, MAXIMUM,
                             2 = timesteps, 3 = hoursPerTimestep ] set/adjusted
                           data->arguments->fileCount  Possibly reduced.
                           data->yyyymmddhh      First timestamp to read.
                           data->readTimesteps   Number of timesteps to read
                                                 at once. 1 or 24 if aggregating
                           data->outputTimesteps  Total number of timesteps to
                                                  output.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkOrSetTimeSubset( Data* const data ) {
  PRE04( data, data->arguments, data->arguments->fileNames[ 0 ],
         data->arguments->fileNames[ 0 ][ 0 ] );
  Arguments* const arguments = data->arguments;
  const int count = arguments->fileCount;
  int yyyymmddhhMin = arguments->subset[ TIME ][ MINIMUM ];
  int yyyymmddhhMax = arguments->subset[ TIME ][ MAXIMUM ];
  int initialized = 0;
  int hoursPerTimestep0 = -1;
  int index = 0;
  int result = 1;
  data->isHourlyTimesteps = 1;
  data->skipFileCount = 0;

  for ( ; AND3( result, index < count, initialized != 2 ); ++index ) {
    const int file = openNetCDFFile( arguments->fileNames[ index ], 'r' );
    result = file >= 0;

    if ( result ) {
      int yyyymmddhh1 = 0;
      int yyyymmddhh2 = 0;
      int timesteps = 0;
      int hoursPerTimestep = 0;
      result = getM3IOFileTimeRange( file, &yyyymmddhh1, &yyyymmddhh2,
                                     &timesteps, &hoursPerTimestep );

      if ( result ) {
        data->fileTimeRange[ index ][ 0 ] = yyyymmddhh1;
        data->fileTimeRange[ index ][ 1 ] = yyyymmddhh2;
        data->fileTimeRange[ index ][ 2 ] = timesteps;
        data->fileTimeRange[ index ][ 3 ] = hoursPerTimestep;

        DEBUG( fprintf( stderr,
                        "%s: %d timesteps, %d hoursPerTimestep = [%d %d]\n",
                        arguments->fileNames[ index ],
                        data->fileTimeRange[ index ][ 2 ],
                        data->fileTimeRange[ index ][ 3 ],
                        data->fileTimeRange[ index ][ 0 ],
                        data->fileTimeRange[ index ][ 1 ] ); )

        if ( hoursPerTimestep > 24 /* 1 */ ) {
          data->isHourlyTimesteps = 0; /* Stride one file per timestep. */
        }

        if ( hoursPerTimestep0 == -1 ) {
          hoursPerTimestep0 = hoursPerTimestep;
        } else if ( AND2( data->isHourlyTimesteps,
                          hoursPerTimestep != hoursPerTimestep0 ) ) {
          fprintf( stderr, "\nMismatched TSTEP size (%d vs %d) in file %s.\n",
                   hoursPerTimestep, hoursPerTimestep0,
                   arguments->fileNames[ index ] );
          result = 0;
        }

        if ( result ) {

          if ( initialized == 0 ) {

            if ( yyyymmddhhMin <= yyyymmddhh2 ) { /* Within file. */

              if ( yyyymmddhhMin < yyyymmddhh1 ) { /* Before file. */
                yyyymmddhhMin = yyyymmddhh1;
              }

              initialized = 1;
            } else {
              data->skipFileCount += 1;
            }
          }

          if ( arguments->subset[ TIME ][ MAXIMUM ] == 0 ) { /* Unset. */
            yyyymmddhhMax = yyyymmddhh2; /* Set to last time of file. */
          } else if ( initialized == 1 ) {

            if ( yyyymmddhhMax <= yyyymmddhh2 ) { /* Within file. */
              initialized = 2;
            } else if ( AND2( count == 1, yyyymmddhh2 <= yyyymmddhhMax ) ) {

              /*
               * If a single file is specified and its last timestep is within
               * the specified time range then reduce the specified time range
               * to the single file's last timestep:
               */

              yyyymmddhhMax = yyyymmddhh2; /* Set to last time of file. */
              initialized = 2;
            }
          }
        }
      }

      closeNetCDFFile( file );
    }
  }

  if ( result ) {

    if ( arguments->subset[ TIME ][ MAXIMUM ] == 0 ) { /* Unset. */
      initialized = 2;
    }

    result = initialized == 2;

    if ( result ) {
      const int hours = hoursInRange( yyyymmddhhMin, yyyymmddhhMax );
      DEBUG( fprintf( stderr, "hours = %d\n", hours ); )

      /* Can't aggregrate daily if less than 24 hours of data per day: */

      if ( OR3( ! data->isHourlyTimesteps,
                data->fileTimeRange[ 0 ][ 3 ] != 1, hours % 24 ) ) {
        arguments->aggregateMode = 0;
      }

      data->yyyymmddhh = yyyymmddhhMin;
      data->readTimesteps = arguments->aggregateMode ? 24 : 1;

      if ( IN3( arguments->aggregateMode, AGGREGATE_MEAN, AGGREGATE_SUM ) ) {
        data->outputTimesteps = 1;
      } else if ( data->isHourlyTimesteps ) {
        data->outputTimesteps = hours / hoursPerTimestep0 / data->readTimesteps;
      } else {
        data->outputTimesteps = index - data->skipFileCount; /* File count. */
      }

      if ( data->outputTimesteps < 1 ) {
        data->outputTimesteps = 1;
      }

      arguments->subset[ TIME ][ MINIMUM ] = yyyymmddhhMin;
      arguments->subset[ TIME ][ MAXIMUM ] = yyyymmddhhMax;
      arguments->fileCount = index - data->skipFileCount; /* Possibly reduced*/
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nInput files do not contain subset time range.\n" );
  }

  DEBUG( fprintf( stderr,
                  "data yyyymmdd = %d, readTimesteps = %d, "
                  "outputTimesteps = %d, "
                  "isHourlyTimesteps = %d, skipFileCount = %d, "
                  "fileCount = %d, "
                  "subset[ TIME ] = [%d %d]\n",
                  data->yyyymmddhh,
                  data->readTimesteps,
                  data->outputTimesteps,
                  data->isHourlyTimesteps,
                  data->skipFileCount,
                  arguments->fileCount,
                  arguments->subset[ TIME ][ MINIMUM ],
                  arguments->subset[ TIME ][ MAXIMUM ] ); )

  POST0( IMPLIES( result,
           AND5( isValidYYYYMMDDHH( arguments->subset[ TIME ][ MINIMUM ] ),
                 isValidYYYYMMDDHH( arguments->subset[ TIME ][ MAXIMUM ] ),
                 arguments->subset[ TIME ][ MINIMUM ] <=
                   arguments->subset[ TIME ][ MAXIMUM ],
                 data->fileTimeRange[ 2 ] > 0,
                 data->fileTimeRange[ 3 ] > 0 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: checkOrSetVariables - Check that variable names are in input file or,
         if the variable names are unspecified then set them to all of those
         in the input file. Also read variable units.
INPUTS:  Data* const data  Data and user-specified arguments to check.
OUTPUTS: Data* const data  variableNames checked or copied,
                           data->layers
                           data->rows
                           data->columns
                           data->variableUnits
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkOrSetVariables( Data* const data ) {

  PRE02( data, data->arguments );
  int result = 0;
  Arguments* const arguments = data->arguments;
  int file = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );

  if ( file ) {
    const int timesteps = getNetCDFDimension( file, "TSTEP" );
    const int variables = getNetCDFDimension( file, "VAR" );
    const int layers    = getNetCDFDimension( file, "LAY" );
    const int rows      = getNetCDFDimension( file, "ROW" );
    const int columns   = getNetCDFDimension( file, "COL" );
    result = GT_ZERO5( timesteps, variables, layers, rows, columns );

    if ( result ) {
      data->layers  = layers;
      data->rows    = rows;
      data->columns = columns;

      if ( arguments->variables ) { /* Check specified variables & get units */
        const int count = arguments->variables;
        int variable = 0;

        for ( ; AND2( result, variable < count ); ++variable ) {
          const int id =
            getNetCDFVariableId( file, arguments->variableNames[ variable ] );
          result = id >= 0;

          if ( result ) {
            char units[ 256 ] = "";
            char var_desc[ 256 ] = "";
            int dimensions[ 32 ] = { 0, 0, 0, 0 };
            int type = 0;
            int rank = 0;
            memset( units, 0, sizeof units );
            memset( dimensions, 0, sizeof dimensions );
            result =
              getNetCDFVariableInfo( file, id,
                                     0, &type, &rank, dimensions, units,
                                     var_desc );

            if ( result ) {
              result =
                AND6( isNetCDFFloat( type ),
                      rank == 4,
                      dimensions[ 0 ] == timesteps,
                      dimensions[ 1 ] == layers,
                      dimensions[ 2 ] == rows,
                      dimensions[ 3 ] == columns );

              if ( ! result ) {
                fprintf( stderr, "\nInvalid data variable '%s'.\n",
                         arguments->variableNames[ variable ] );
              } else {
                strncpy( data->variableUnits[ variable ], units, NAMLEN3 );
                strncpy( data->variableDescriptions[ variable ], var_desc,
                         MXDLEN3 );
              }
            }
          }
        }

      } else { /* Copy all variable names and units: */
        const int count = MIN( variables, MXVARS3 );
        int variable = 0;
        int variableCount = 0;

        for ( ; AND2( result, variable < count ); ++variable ) {
          char name[  256 ] = "";
          char units[ 256 ] = "";
          char var_desc[ 256 ] = "";
          int dimensions[ 32 ] = { 0, 0, 0, 0 };
          int type = 0;
          int rank = 0;
          memset( name, 0, sizeof name );
          memset( units, 0, sizeof units );
          memset( dimensions, 0, sizeof dimensions );
          result =
            getNetCDFVariableInfo( file, variable + 1, /* + 1 to skip TFLAG. */
                                   name, &type, &rank, dimensions, units,
                                   var_desc );
          if ( result ) {

            if ( AND6( isNetCDFFloat( type ),
                       rank == 4,
                       dimensions[ 0 ] == timesteps,
                       dimensions[ 1 ] == layers,
                       dimensions[ 2 ] == rows,
                       dimensions[ 3 ] == columns ) ) {
              strncpy( arguments->variableNames[variableCount], name, NAMLEN3);
              strncpy( data->variableUnits[ variableCount ], units, NAMLEN3 );
              strncpy( data->variableDescriptions[ variableCount ], var_desc,
                       MXDLEN3 );
              ++variableCount;
            }
          }
        }

        arguments->variables = variableCount;
      }

      if ( AND2( result, arguments->auxMode == INTEGRATE ) ) {

        /* Check that variable units are ppmV or ppbV: */

        const int count = arguments->variables;
        int variable = 0;

        for ( ; AND2( result, variable < count ); ++variable ) {
          result =
            OR2( ! strcmp( data->variableUnits[ variable ], "ppmV" ),
                 ! strcmp( data->variableUnits[ variable ], "ppbV" ) );

          if ( ! result ) {
            fprintf( stderr,
                     "\nInvalid units '%s' (require ppmV or ppbV) "
                     "for integration variable '%s'.\n",
                     data->variableUnits[ variable ],
                     arguments->variableNames[ variable ] );
          }
        }
      }
    }

    closeNetCDFFile( file ), file = -1;
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND8( data->arguments->variables > 0,
                         data->arguments->variableNames[ 0 ],
                         data->arguments->variableNames[
                                              data->arguments->variables - 1 ],
                         data->variableUnits[ 0 ],
                         data->variableUnits[ data->arguments->variables - 1 ],
                         data->variableDescriptions[ 0 ],
                         data->variableDescriptions[
                                              data->arguments->variables - 1 ],
                         GT_ZERO3(data->layers, data->rows, data->columns))));
  return result;
}



/******************************************************************************
PURPOSE: checkOrSetGridSubset - Check subset subset against the input or,
         if the subset is unspecified then set them to full range of
         input file.
INPUTS:  Data* const data  data->layers
                           data->rows
                           data->columns
                           data->arguments->subset[][] subset to check/set.
OUTPUTS: Data* const data  data->arguments->subset[][] set/adjusted.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkOrSetGridSubset( Data* const data ) {
  PRE05( data,
         data->arguments,
         data->arguments->fileNames[ 0 ],
         data->arguments->fileNames[ 0 ][ 0 ],
         GT_ZERO3( data->layers, data->rows, data->columns ) );
  Arguments* const arguments = data->arguments;
  int result = 0;

  if ( arguments->subset[ LAYER ][ MINIMUM ] == 0 ) { /* Initialize: */
    arguments->subset[ LAYER  ][ MINIMUM ] = 1;
    arguments->subset[ LAYER  ][ MAXIMUM ] = data->layers;
    result = 1;
  } else { /* Check: */
    result =
      AND2( IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ], 1, data->layers ),
            IN_RANGE( arguments->subset[ LAYER ][ MAXIMUM ],
                      arguments->subset[ LAYER ][ MINIMUM ], data->layers ) );
  }

  if ( result ) {

    if ( arguments->subset[ ROW ][ MINIMUM ] == 0 ) { /* Initialize: */
      arguments->subset[ ROW  ][ MINIMUM ] = 1;
      arguments->subset[ ROW  ][ MAXIMUM ] = data->rows;
      result = 1;
    } else { /* Check: */
      result =
        AND2( IN_RANGE( arguments->subset[ ROW ][ MINIMUM ], 1, data->rows ),
              IN_RANGE( arguments->subset[ ROW ][ MAXIMUM ],
                        arguments->subset[ ROW ][ MINIMUM ], data->rows ) );
    }

    if ( result ) {

      if ( arguments->subset[ COLUMN ][ MINIMUM ] == 0 ) { /* Initialize: */
        arguments->subset[ COLUMN  ][ MINIMUM ] = 1;
        arguments->subset[ COLUMN  ][ MAXIMUM ] = data->columns;
        result = 1;
      } else { /* Check: */
        result =
          AND2( IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                          1, data->columns ),
                IN_RANGE( arguments->subset[ COLUMN ][ MAXIMUM ],
                          arguments->subset[ COLUMN ][ MINIMUM ],
                          data->columns ) );
      }
    }
  }

  POST0( IMPLIES( result,
           AND6( IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                           1, data->layers ),
                 IN_RANGE( arguments->subset[ LAYER ][ MAXIMUM ],
                           arguments->subset[ LAYER ][ MINIMUM ],
                           data->layers ),
                 IN_RANGE( arguments->subset[ ROW ][ MINIMUM ],
                           1, data->rows ),
                 IN_RANGE( arguments->subset[ ROW ][ MAXIMUM ],
                           arguments->subset[ ROW ][ MINIMUM ], data->rows ),
                 IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                           1, data->columns ),
                 IN_RANGE( arguments->subset[ COLUMN ][ MAXIMUM ],
                           arguments->subset[ COLUMN ][ MINIMUM ],
                           data->columns ) ) ) );

  return result;
}



/******************************************************************************
PURPOSE: checkFilesAreCompatible - Check that the set of input files are
         compatible.
INPUTS:  Data* const data  Data to check.
RETURNS: int 1 if compatible, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkFilesAreCompatible( Data* const data ) {

  PRE013( data, data->arguments, data->arguments->fileNames[ 0 ],
          data->arguments->fileNames[ 0 ][ 0 ],
          GT_ZERO3( data->layers, data->rows, data->columns ),
          isValidYYYYMMDDHH( data->fileTimeRange[ 0 ][ MINIMUM ] ),
          isValidYYYYMMDDHH( data->fileTimeRange[ 0 ][ MAXIMUM ] ),
          data->fileTimeRange[ 0 ][ 2 ] > 0,
          data->fileTimeRange[ 0 ][ 3 ] > 0,
          isValidYYYYMMDDHH( data->fileTimeRange[
                               data->arguments->fileCount - 1 ][ MINIMUM ] ),
          isValidYYYYMMDDHH( data->fileTimeRange[
                               data->arguments->fileCount - 1 ][ MAXIMUM ] ),
          data->fileTimeRange[ data->arguments->fileCount - 1 ][ 2 ]
            == data->fileTimeRange[ 0 ][ 2 ],
          IMPLIES( data->isHourlyTimesteps,
                   data->fileTimeRange[ data->arguments->fileCount - 1 ][ 3 ]
                     == data->fileTimeRange[ 0 ][ 3 ] ) );

  Arguments* const arguments = data->arguments;
  int count = arguments->fileCount;
  int result = 1;
  int index = 0;

  /* Check that each input file has the specified variables on the same grid:*/

  for ( ; AND2( result, index < count ); ++index ) {
    result = checkFileVariables( data );
  }

  /* Check that the HT (GRIDCRO2D) file has variable HT on compatible 2D grid*/

  if ( AND2( result, arguments->htFileName ) ) {
    result = checkHTFile( data );
  }

  /*
   * Check that ZF files (METCRO3D) have variable ZH on a compatible grid
   * and, if integrating layers, variables ZF, DENS on a matched grid:
   */

  if ( AND2( result, arguments->zfFileCount ) ) {
    result = checkZFFiles( data );
  }

  /* Check that the WWIND (METCRO3D) files contain variable WWIND: */

  if ( AND2( result, arguments->wwindFileCount > 0 ) ) {
    result = checkWWINDFiles( data );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: checkFileVariables - Check that each input file has all of the
         specified variables on the matched 3D grid.
INPUTS:  Data* const data  Data to check.
RETURNS: int 1 if compatible, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkFileVariables( const Data* const data ) {

  PRE011( data,
          data->arguments,
          data->arguments->fileCount,
          data->arguments->variables > 0,
          data->arguments->fileNames[ 0 ][ 0 ],
          data->arguments->fileNames[ data->arguments->fileCount - 1 ][ 0 ],
          data->arguments->variableNames[ 0 ][ 0 ],
          data->arguments->variableNames[ data->arguments->variables - 1 ][ 0],
          data->variableUnits[ 0 ][ 0 ],
          data->variableUnits[ data->arguments->variables - 1 ][ 0 ],
          GT_ZERO3( data->layers, data->rows, data->columns ) );

  Arguments* const arguments = data->arguments;
  const int fileCount = arguments->fileCount;
  const int variables = data->arguments->variables;
  int result = 1;
  int index = 0;

  for ( ; AND2( result, index < fileCount ); ++index ) {
    int file = openNetCDFFile( arguments->fileNames[ index ], 'r' );
    result = file != -1;

    if ( result ) {
      int variable = 0;

      for ( ; AND2( result, variable < variables ); ++variable ) {
        const int id =
          getNetCDFVariableId( file, arguments->variableNames[ variable ] );
        result = id != -1;

        if ( result ) {
          char units[ 256 ] = "";
          int dimensions[ 32 ] = { 0, 0, 0, 0 };
          int type = 0;
          int rank = 0;
          memset( units, 0, sizeof units );
          memset( dimensions, 0, sizeof dimensions );
          result =
            getNetCDFVariableInfo( file, id, 0, &type, &rank, dimensions,
                                   units, 0 );

          if ( result ) {
            result =
              AND6( isNetCDFFloat( type ),
                    rank == 4,
                    ! strcmp( units, data->variableUnits[ variable ] ),
                    dimensions[ 1 ] == data->layers,
                    dimensions[ 2 ] == data->rows,
                    dimensions[ 3 ] == data->columns );

            if ( ! result ) {
              fprintf( stderr, "\nInvalid input file '%s' "
                       "has incompatible variable '%s', units (%s).\n",
                       arguments->fileNames[ index ],
                       arguments->variableNames[ variable ], units );
            }
          }
        }
      }
    }

    closeNetCDFFile( file ), file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: checkHTFile - Check that the HT (GRIDCRO2D) file has variable HT on a
         compatible 2D grid
INPUTS:  Data* const data  Data to check.
RETURNS: int 1 if compatible, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkHTFile( const Data* const data ) {

  PRE05( data, data->arguments, data->arguments->htFileName,
         data->arguments->htFileName[ 0 ],
         GT_ZERO3( data->layers, data->rows, data->columns ) );

  Arguments* const arguments = data->arguments;
  int file = openNetCDFFile( arguments->htFileName, 'r' );
  int result = file != -1;

  if ( result ) {
    const int id = getNetCDFVariableId( file, "HT" );
    result = id != -1;

    if ( result ) {
      char units[ 256 ] = "";
      int dimensions[ 32 ] = { 0, 0, 0, 0 };
      int type = 0;
      int rank = 0;
      memset( units, 0, sizeof units );
      memset( dimensions, 0, sizeof dimensions );
      result =
        getNetCDFVariableInfo(file, id, 0, &type, &rank, dimensions, units, 0);

      if ( result ) {
        result =
          AND7( isNetCDFFloat( type ),
                rank == 4,
                OR2( units[ 0 ] == 'm', units[ 0 ] == 'M' ),
                units[ 1 ] == '\0',
                dimensions[ 0 ] == 1,
                dimensions[ 1 ] == 1,
                IMPLIES_ELSE( dimensions[ 2 ] == data->rows,
                              dimensions[ 3 ] == data->columns,
                              AND2( dimensions[ 2 ] == data->rows - 1,
                                    dimensions[ 3 ] == data->columns - 1 ) ) );

        if ( ! result ) {
          fprintf( stderr, "\nInvalid HT file specified '%s'.\n",
                   arguments->htFileName );
        }
      }
    }

    closeNetCDFFile( file ), file = -1;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: checkZFFiles - Check that the ZF (METCRO3D) files have variable
         ZH, ZF and DENS on a compatible 3D CRO grid.
INPUTS:  Data* const data  Data to check.
OUTPUTS: Data* const data  data->zfFileTimeRange[][2]
RETURNS: int 1 if compatible, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkZFFiles( Data* const data ) {

  PRE05( data, data->arguments, data->arguments->zfFileCount,
         data->arguments->zfFileNames[ 0 ][ 0 ],
         GT_ZERO3( data->layers, data->rows, data->columns ) );

  Arguments* const arguments = data->arguments;
  const int count = arguments->zfFileCount;
  int index = 0;
  int result = 1;

  for ( ; AND2( result, index < count ); ++index ) {
    int file = openNetCDFFile( arguments->zfFileNames[ index ], 'r' );
    result = file != -1;

    if ( result ) { /* Get file time range for later use during reading: */
      result =
        getM3IOFileTimeRange( file,
                             &data->zfFileTimeRange[ index ][ MINIMUM ],
                             &data->zfFileTimeRange[ index ][ MAXIMUM ],
                             &data->zfFileTimeRange[ index ][ 2 ],
                             &data->zfFileTimeRange[ index ][ 3 ] );

      if ( result ) {

        if ( data->zfFileTimeRange[ index ][ 3 ] !=
             data->fileTimeRange[ 0 ][ 3 ] ) {
          fprintf( stderr, "\nMismatched TSTEP size (%d vs %d) in file %s.\n",
                   data->zfFileTimeRange[ index ][ 3 ],
                   data->fileTimeRange[ 0 ][ 3 ],
                   arguments->zfFileNames[ index ] );
          result = 0;
        } else {
          int id = getNetCDFVariableId( file, "ZH" );
          result = id != -1;

          if ( result ) {

            /*
             * Check that ZH row/col dims match or are 1 less than DOT dims.
             * Note multi-layer ZH can be used with single-layer data (ACONC).
             */

            int dimensions[ 32 ] = { 0, 0, 0, 0 };
            int type = 0;
            int rank = 0;
            memset( dimensions, 0, sizeof dimensions );
            result =
              getNetCDFVariableInfo( file, id, 0, &type, &rank, dimensions,
                                     0, 0 );

            if ( result ) {
              result =
                AND4( isNetCDFFloat( type ),
                      rank == 4,
                      dimensions[ 1 ] >= data->layers,
                      IMPLIES_ELSE( dimensions[ 2 ] == data->rows,
                                    dimensions[ 3 ] == data->columns,
                                    AND2( dimensions[ 2 ] == data->rows - 1,
                                          dimensions[ 3 ] == data->columns - 1)));

              if ( result) {

                if ( arguments->auxMode == INTEGRATE ) {

                  /* Check that ZF, DENS 3D dims match CRO variable dims: */

                  id = getNetCDFVariableId( file, "ZF" );
                  result = id != -1;

                  if ( result ) {
                    result =
                      getNetCDFVariableInfo( file, id, 0, &type, &rank,
                                             dimensions, 0, 0 );

                    if ( result ) {
                      result =
                        AND5( isNetCDFFloat( type ),
                              rank == 4,
                              dimensions[ 1 ] == data->layers,
                              dimensions[ 2 ] == data->rows,
                              dimensions[ 3 ] == data->columns );

                      if ( result ) {
                        id = getNetCDFVariableId( file, "DENS" );
                        result = id != -1;

                        if ( result ) {
                          result =
                            getNetCDFVariableInfo( file, id, 0, &type, &rank,
                                                   dimensions, 0, 0 );

                          if ( result ) {
                            result =
                              AND5( isNetCDFFloat( type ),
                                   rank == 4,
                                   dimensions[ 1 ] == data->layers,
                                   dimensions[ 2 ] == data->rows,
                                   dimensions[ 3 ] == data->columns );
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

      closeNetCDFFile( file ), file = -1;
    }

    if ( ! result ) {
      fprintf( stderr, "\nInvalid ZF file specified '%s'.\n",
               arguments->zfFileNames[ index ] );
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND4( isValidYYYYMMDDHH( data->zfFileTimeRange[0][0] ),
                         isValidYYYYMMDDHH( data->zfFileTimeRange[0][1] ),
                         isValidYYYYMMDDHH( data->zfFileTimeRange[
                            data->arguments->zfFileCount - 1 ][ 0 ] ),
                         isValidYYYYMMDDHH( data->zfFileTimeRange[
                            data->arguments->zfFileCount - 1 ][ 1 ] ) ) ) );

  return result;
}



/******************************************************************************
PURPOSE: checkWWINDFiles - Check that the WWIND (METCRO3D) file has variable
         WWIND on a compatible 3D CRO grid.
INPUTS:  Data* const data  Data to check.
OUTPUTS: Data* const data  data->wwindFileTimeRange[][2],
                           wwindVariable is "WWIND" or "W_VEL".
RETURNS: int 1 if compatible, else 0 and a failure message is printed to stderr
******************************************************************************/

static int checkWWINDFiles( Data* const data ) {

  PRE05( data, data->arguments, data->arguments->wwindFileCount,
         data->arguments->wwindFileNames[ 0 ][ 0 ],
         GT_ZERO3( data->layers, data->rows, data->columns ) );

  Arguments* const arguments = data->arguments;
  const int count = arguments->wwindFileCount;
  int index = 0;
  int result = 1;
  data->wwindVariable = 0;

  for ( ; AND2( result, index < count ); ++index ) {
    int file = openNetCDFFile( arguments->wwindFileNames[ index ], 'r' );
    result = file != -1;

    if ( result ) { /* Get file time range for later use during reading: */
      result =
        getM3IOFileTimeRange( file,
                              &data->wwindFileTimeRange[ index ][ MINIMUM ],
                              &data->wwindFileTimeRange[ index ][ MAXIMUM ],
                              &data->wwindFileTimeRange[ index ][ 2 ],
                              &data->wwindFileTimeRange[ index ][ 3 ] );

      if ( result ) { /* Check that CRO dims are 1 less than DOT U/VWIND dims*/

        if ( data->wwindFileTimeRange[ index ][ 3 ] !=
             data->fileTimeRange[ 0 ][ 3 ] ) {
          fprintf( stderr, "\nMismatched TSTEP size (%d vs %d) in file %s.\n",
                   data->wwindFileTimeRange[ index ][ 3 ],
                   data->fileTimeRange[ 0 ][ 3 ],
                   arguments->wwindFileNames[ index ] );
          result = 0;
        } else {
          int id = -1;

          if ( ! data->wwindVariable ) { /* If uninitialized, initialize now.*/
            data->wwindVariable = "WWIND"; /* Might have WWIND (default). */
            id = checkNetCDFVariableId( file, data->wwindVariable );

            if ( id == -1 ) {
              data->wwindVariable = "W_VEL"; /* Otherwise must have W_VEL. */
              id = getNetCDFVariableId( file, data->wwindVariable );
            }
          } else { /* Must have same wwind variable as first file. */
            id = getNetCDFVariableId( file, data->wwindVariable );
          }

          result = id != -1;

          if ( result ) {
            int dimensions[ 32 ] = { 0, 0, 0, 0 };
            int type = 0;
            int rank = 0;
            memset( dimensions, 0, sizeof dimensions );
            result =
              getNetCDFVariableInfo( file, id, 0, &type, &rank, dimensions,
                                     0, 0 );

            /* Check METCRO3D or CCTM_CONC grid 'matches' METDOT3D grid. */

            if ( result ) {
              result =
                AND5( isNetCDFFloat( type ),
                      rank == 4,
                      dimensions[ 1 ] == data->layers,
                      dimensions[ 2 ] == data->rows - 1,
                      dimensions[ 3 ] == data->columns - 1 );
            }
          }
        }
      }

      closeNetCDFFile( file ), file = -1;
    }

    if ( ! result ) {
      fprintf( stderr, "\nInvalid WWIND file specified '%s'.\n",
               arguments->wwindFileNames[ index ] );
    }
  }

  result = AND2( result, data->wwindVariable != 0 );

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND5( isValidYYYYMMDDHH( data->wwindFileTimeRange[0][0] ),
                         isValidYYYYMMDDHH( data->wwindFileTimeRange[0][1] ),
                         isValidYYYYMMDDHH( data->wwindFileTimeRange[
                            data->arguments->wwindFileCount - 1 ][ 0 ] ),
                         isValidYYYYMMDDHH( data->wwindFileTimeRange[
                            data->arguments->wwindFileCount - 1 ][ 1 ] ),
                         OR2( ! strcmp( data->wwindVariable, "WWIND" ),
                              ! strcmp( data->wwindVariable, "W_VEL" ) ) ) ) );

  return result;
}



/******************************************************************************
PURPOSE: computeGridCellCoordinates - Compute grid cell edge lon-lats.
INPUTS:  Data* const data  Data to initialize.
OUTPUTS: Data* data        data->longitudes, data->latitudes.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int computeGridCellCoordinates( Data* const data ) {
  PRE03( data, data->arguments, IS_ZERO2( data->longitudes, data->latitudes ));
  const Arguments* const arguments = data->arguments;
  const double majorSemiaxis = arguments->ellipsoid[ MAXIMUM ];
  const double minorSemiaxis = arguments->ellipsoid[ MINIMUM ];
  int file = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );
  int result = file != -1;

  if ( result ) {
    int nrows = 0;
    int ncols = 0;
    double xorig = 0.0;
    double yorig = 0.0;
    double xcell = 0.0;
    double ycell = 0.0;
    int gdtyp = 0;
    result =
      AND11( getNetCDFIntAttribute( file, "NROWS", &nrows ),
             getNetCDFIntAttribute( file, "NCOLS", &ncols ),
             getNetCDFDoubleAttribute( file, "XORIG", &xorig ),
             getNetCDFDoubleAttribute( file, "YORIG", &yorig ),
             getNetCDFDoubleAttribute( file, "XCELL", &xcell ),
             getNetCDFDoubleAttribute( file, "YCELL", &ycell ),
             getNetCDFIntAttribute( file, "GDTYP", &gdtyp ),
             nrows > 0, ncols > 0, xcell > 0.0, ycell > 0.0 );

    if ( result ) {
      const size_t count = ( nrows + 1 ) * ( ncols + 1 );
      data->longitudes = NEW( double, count );
      data->latitudes = data->longitudes ? NEW( double, count ) : 0;
      result = data->latitudes != 0;

      if ( result ) {
        Projector* projector = 0;
        DEBUG( fprintf( stderr, "gdtyp = %d\n", gdtyp ); )

        /* Read CMAQ projection parameters: */

        switch ( gdtyp ) {
        case LATGRD3:
          break;
        case ALBGRD3:
          /* lint -fallthrough */
        case LAMGRD3:
          {
            double p_alp = 0.0;
            double p_bet = 0.0;
            double xcent = 0.0;
            double ycent = 0.0;
            result =
              AND14( getNetCDFDoubleAttribute( file, "P_ALP", &p_alp ),
                     getNetCDFDoubleAttribute( file, "P_BET", &p_bet ),
                     getNetCDFDoubleAttribute( file, "XCENT", &xcent ),
                     getNetCDFDoubleAttribute( file, "YCENT", &ycent ),
                     isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                     isValidLatitude( p_alp ),
                     isValidLatitude( p_bet ),
                     isValidLongitude( xcent ),
                     isValidLatitude( ycent ),
                     p_alp <= p_bet,
                     SIGN( p_alp ) == SIGN( p_bet ),
                     IMPLIES_ELSE( p_alp >= 0.0,
                                   IN_RANGE( p_alp, 1.0, 89.0 ),
                                   IN_RANGE( p_alp, -89.0, -1.0 ) ),
                     IMPLIES_ELSE( p_bet >= 0.0,
                                   IN_RANGE( p_bet, 1.0, 89.0 ),
                                   IN_RANGE( p_bet, -89.0, -1.0 ) ),
                     IN_RANGE( ycent, -89.0, 89.0 ) );

            if ( result ) {

              if ( gdtyp == LAMGRD3 ) {
                projector = (Projector*)
                  newLambert( majorSemiaxis, minorSemiaxis,
                              p_alp, p_bet, xcent, ycent, 0.0, 0.0 );
              } else {
                CHECK( gdtyp == ALBGRD3 );
                projector = (Projector*)
                  newAlbers( majorSemiaxis, minorSemiaxis,
                             p_alp, p_bet, xcent, ycent, 0.0, 0.0 );
              }

              result = projector != 0;
            }
          }
          break;
        case EQMGRD3:
          {
            double xcent = 0.0;
            result =
              AND3( getNetCDFDoubleAttribute( file, "XCENT", &xcent ),
                    isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                    isValidLongitude( xcent ) );

            if ( result ) {
              projector = (Projector*)
                newMercator( majorSemiaxis, minorSemiaxis, xcent, 0.0, 0.0 );
              result = projector != 0;
            }
          }
          break;
        case POLGRD3:
          {
            double p_bet = 0.0;
            double xcent = 0.0;
            double ycent = 0.0;
            result =
              AND7( getNetCDFDoubleAttribute( file, "P_BET", &p_bet ),
                    getNetCDFDoubleAttribute( file, "XCENT", &xcent ),
                    getNetCDFDoubleAttribute( file, "YCENT", &ycent ),
                    isValidEllipsoid( majorSemiaxis, minorSemiaxis ),
                    isValidLongitude( xcent ),
                    isValidLatitude( ycent ),
                    isValidLatitude( p_bet ) );

            if ( result ) {
              projector = (Projector*)
                newStereographic( majorSemiaxis, minorSemiaxis,
                                  xcent, ycent, p_bet, 0.0, 0.0 );
              result = projector != 0;
            }
          }
          break;
        default:
          fprintf(stderr, "\nUnsupported projection type GDTYP = %d.\n",gdtyp);
          result = 0;
          break;
        }

        if ( IS_ZERO2( result, projector ) ) {
          fputs( "\nRead invalid projection parameters.\n", stderr );
        } else {

          /*
           * Unproject (x, y) grid cell corner points to (longitude, latitude)
           * and store them in data->longitudes[], data->latitudes[]:
           */

          double* longitudes = data->longitudes;
          double* latitudes  = data->latitudes;
          double y = yorig;
          int row = 0;

          for ( ; row <= nrows; ++row, y += ycell ) {
            double x = xorig;
            int column = 0;

            for ( ; column <= ncols; ++column, x += xcell ) {
              double longitude = x;
              double latitude  = y;

              if ( projector ) {
                projector->unproject( projector, x, y, &longitude, &latitude );
              }

#ifdef DEBUGGING
              if ( OR2( AND2( row <= 1, column <= 1 ),
                        AND2( row >= nrows - 1, column >= ncols - 1 ) ) ) {
                fprintf( stderr, "[row %d column %d] = (%f %f) -> (%f %f)\n",
                         row, column, x, y, longitude, latitude );
              }
#endif

              CHECK( isValidLongitudeLatitude( longitude, latitude ) );
              *longitudes++ = longitude;
              *latitudes++  = latitude;
            }
          }

          CHECK2( isValidLongitudeLatitude( data->longitudes[ 0 ],
                                            data->latitudes[ 0 ] ),
                  isValidLongitudeLatitude( data->longitudes[( nrows + 1 ) *
                                                             ( ncols + 1 ) -1],
                                            data->latitudes[( nrows + 1 ) *
                                                            ( ncols + 1 ) -1]));
        }

        data->isProjected = projector != 0;
        FREE_OBJECT( projector );
      }
    }

    closeNetCDFFile( file ), file = -1;
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND4( data->longitudes, data->latitudes,
                         IN_RANGE( data->longitudes[ 0 ], -180.0, 180.0 ),
                         IN_RANGE( data->latitudes[  0 ],  -90.0,  90.0 ) )));
  return result;
}


/******************************************************************************
PURPOSE: computeGridCellCenterElevations - Compute grid cell center elevations.
INPUTS:  Data* const data  Data to initialize.
OUTPUTS: Data* data        data->elevations.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int computeGridCellCenterElevations( Data* const data ) {
  PRE03( data,
         data->arguments,
         GT_ZERO3( data->layers, data->rows, data->columns ) );
  const Arguments* const arguments = data->arguments;
  int file = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );
  int result = file != -1;

  if ( result ) {
    const int layers = data->layers;
    int vgtyp = 0;
    float vgtop = 0.0;
    float vglvls[ MXLAYS3 + 1 ];
    memset( vglvls, 0, sizeof vglvls );
    result =
      AND3( getNetCDFIntAttribute( file, "VGTYP", &vgtyp ),
            getNetCDFFloatAttribute( file, "VGTOP", &vgtop ),
            getNetCDFFloatArrayAttribute( file, "VGLVLS", layers + 1, vglvls ));
    closeNetCDFFile( file ), file = -1;
    checkAndFixVerticalGridParameters( layers, &vgtyp, &vgtop, vglvls );

    if ( result ) {
      const size_t count = layers * data->rows * data->columns;
      double* z = NEW( double, data->layers + 1 );
      data->elevations = z ? NEW( double, count ) : 0;
      result = data->elevations != 0;

      if ( ! result ) {
        FREE( z );
      } else {
        const double topPressure = 10000.0; /* Pascals. */
        const double g           = 9.81;
        const double R           = 287.04;
        const double A           = 50.0;
        const double T0s         = 290.0;
        const double P00         = 100000.0;

        DEBUG( fprintf( stderr,
                        "vgtyp = %d, vgtop = %f, vglvls = [%f ... %f]\n",
                         vgtyp, vgtop, vglvls[ 0 ], vglvls[ data->layers ] ); )


        if ( data->heights == 0 ) {

          /*
           * No HT is available so compute z[] for just one sea-level cell
           * and replicate it for all grid cells.
           */

          double* elevations = data->elevations;
          const size_t cells = data->rows * data->columns;
          int layer = 0;

          computeZ( g, R, A, T0s, P00, layers, vgtyp, topPressure,
                    0.0, vglvls, z );

          /* Copy z values to all grid cells: */

          for ( layer = 0; layer < layers; ++layer ) {
            const double zMidLayer = ( z[ layer ] + z[ layer + 1 ] ) * 0.5;
            const double cellCenterElevation =
              CLAMPED_TO_RANGE( zMidLayer,
                                ELEVATION_MINIMUM, ELEVATION_MAXIMUM );
            size_t counter = cells;

            DEBUG( fprintf( stderr, "zMidLayer at %d = %f\n",
                            layer, zMidLayer ); )

            while ( counter-- ) {
              *elevations++ = cellCenterElevation;
            }
          }

        } else { /* Use HT and compute z[] above each surface cell: */
          const float* heights = data->heights;
          double* elevations = data->elevations;
          double previousHeight = BADVAL3;
          const size_t cells = data->rows * data->columns;
          size_t cell = 0;
          int layer = 0;

          for ( cell = 0; cell < cells; ++cell ) {
            const double height = heights[ cell ];

            if ( height != previousHeight ) {
              computeZ( g, R, A, T0s, P00, layers, vgtyp, topPressure,
                        height, vglvls, z );
              previousHeight = height;
            }

            /* Copy z values to all layers above the cell: */

            for ( layer = 0; layer < layers; ++layer ) {
              const size_t index = layer * cells + cell;
              const double zMidLayer = ( z[ layer ] + z[ layer + 1 ] ) * 0.5;
              const double cellCenterElevation =
                CLAMPED_TO_RANGE( zMidLayer,
                                  ELEVATION_MINIMUM, ELEVATION_MAXIMUM );

              DEBUG( if ( cell < 3 || cell > cells - 3 ) fprintf( stderr,
                              "elevations[ %lu ] = zMidLayer at %d = %f\n",
                              index, layer, zMidLayer ); )

              CHECK( index < data->layers * data->rows * data->columns );
              elevations[ index ] = cellCenterElevation;
            }
          }
        }
      }
    }
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND2( data->elevations,
                         inRange( data->layers * data->rows * data->columns,
                                  data->elevations,
                                  ELEVATION_MINIMUM, ELEVATION_MAXIMUM ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: checkAndFixVerticalGridParameters - Check vertical grid parameters and
         edit them as needed to make them valid.
INPUTS:  const int layers  Number of layers.
         int* vgtyp        Input vertical grid type.
         float* vgtop      Pressure in Pascals of the top of the atmosphere.
         float vglvls[ layers + 1 ]  Vertical grid levels.
OUTPUTS: int* vgtyp        Valid vertical grid type.
         float* vgtop      Valid pressure in Pascals of the top of the atmos.
         float vglvls[ layers + 1 ]  Valid vertical grid levels.
******************************************************************************/

static void checkAndFixVerticalGridParameters( const int layers,
                                               int* vgtyp,
                                               float* vgtop,
                                               float vglvls[] ) {
  PRE04( layers > 0, vgtyp, vgtop, vglvls );

  DEBUG( fprintf( stderr, "before "
                  "check/fix: vgtyp = %d, vgtop = %f, vglvls = [%f ... %f]\n",
                  *vgtyp, *vgtop, vglvls[ 0 ], vglvls[ layers ] ) );

  if ( ! IS_VALID_VERTICAL_GRID_TYPE( *vgtyp ) ) {
    *vgtyp = VGSGPN3; /* Non-hydrostatic sigma-P. */
  }

  if ( *vgtop <= 0.0 ) {
    *vgtop = 5000.0;
  }

  if ( ! IN4( *vgtyp, VGPRES3, VGZVAL3, VGHVAL3 ) ) {
    float previousValue = vglvls[ 0 ];
    int level = 0;
    int ok = 1;

    for ( level = 1; AND2( ok, level <= layers ); ++level ) {
      const float value = vglvls[ level ];
      ok = AND2( value < previousValue, IN_RANGE( value, 0.0, 1.0 ) );
      previousValue = value;
    }

    if ( ! ok ) {
      const double delta = 1.0 / MXLAYS3;
      double sigma = 1.0;

      for ( level = 0; level <= layers; ++level, sigma -= delta ) {
        vglvls[ level ] = sigma;
      }
    }
  }

  DEBUG( fprintf( stderr, "after  "
                  "check/fix: vgtyp = %d, vgtop = %f, vglvls = [%f ... %f]\n",
                  *vgtyp, *vgtop, vglvls[ 0 ], vglvls[ layers ] ) );
}



/******************************************************************************
PURPOSE: readHT - Read variable HT from HT file into grid-matched heights.
INPUTS:  Data* const data  data->htFileName.
OUTPUTS: const Data* const data  data->height[ rows * columns ]
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readHT( Data* const data ) {

  PRE05( data, data->arguments, isValidArguments( data->arguments ),
         data->arguments->htFileName, data->heights == 0 );

  const Arguments* const arguments = data->arguments;
  int file = openNetCDFFile( arguments->htFileName, 'r' );
  int result = file != -1;

  if ( result ) {
    int dimsHT[ 4 ] = { 0, 0, 0, 0 }; /* TIME=1, LAYER=1, ROW, COLUMN. */
    result = getM3IOVariableDimensions( file, "HT", dimsHT );
    result = AND3( result, dimsHT[ TIME ] == 1, dimsHT[ LAYER ] == 1 );

    if ( result ) {
      const size_t countHT = dimsHT[ ROW ] * dimsHT[ COLUMN ];
      float* ht = NEW( float, countHT );
      result = ht != 0;

      if ( result ) {
        int id = getNetCDFVariableId( file, "HT" );
        result = id >= 0;

        if ( result ) {
          result = readM3IOVariable( file, id, 0, 0, 0, 0,
                                     0, dimsHT[ ROW    ] - 1,
                                     0, dimsHT[ COLUMN ] - 1,
                                     ht );
          if ( result ) {
            closeNetCDFFile( file );
            file = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );
            result = file != -1;

            if ( result ) {
              int dims[ 4 ] = { 0, 0, 0, 0 }; /* TIME, LAYER, ROW, COLUMN. */
              result = getM3IOVariableDimensions( file,
                                                  arguments->variableNames[0],
                                                  dims );

              /* If input file is a DOT grid then extend ht at edges: */

              if ( AND2( dims[ ROW    ] == dimsHT[ ROW    ] + 1,
                         dims[ COLUMN ] == dimsHT[ COLUMN ] + 1 ) ) {
                const int rows    = dims[ ROW    ];
                const int columns = dims[ COLUMN ];
                data->heights = NEW( float, rows * columns );
                result = data->heights != 0;

                if ( result ) {
                  float* input = ht;
                  float* output = data->heights;
                  const int rows_1    = rows - 1;
                  const int columns_1 = columns - 1;
                  int row = 0;

                  for ( ; row < rows_1; ++row ) {
                    int column = 0;

                    for ( ; column < columns_1; ++column ) {
                      *output++ = *input++;
                    }

                    *output++ = *( input - 1 ); /* Replicate last value. */
                  }

                  /* Replicate last row: */

                  input = output - columns;
                  memcpy( output, input, columns * sizeof *input );
                }
              } else { /* Input is also a CRO grid: */
                CHECK( AND2( dims[ ROW    ] == dimsHT[ ROW    ],
                             dims[ COLUMN ] == dimsHT[ COLUMN ] ) );
                data->heights = ht; /* Transfer ownership of pointer. */
                ht = 0;
              }
            }
          }
        }

        FREE( ht );
      }
    }

    if ( file != -1 ) {
      closeNetCDFFile( file ), file = -1;
    }
  }

  POST02( IS_BOOL( result ), IMPLIES( result, data->heights ) );
  return result;
}




/******************************************************************************
PURPOSE: boundsSubset - Apply bounds to possibly reduce subset.
INPUTS:  Data* const data  data->arguments->bounds, subset.
OUTPUTS: const Data* const data  data->arguments->subset.
RETURNS: int 1 if subset intersects, else 0 and a failure message is printed
         to stderr.
******************************************************************************/

static int boundsSubset( Data* const data ) {

  PRE07( data, data->arguments, isValidArguments( data->arguments ),
         IN_RANGE( data->arguments->subset[ ROW ][ MINIMUM ], 1, data->rows ),
         IN_RANGE( data->arguments->subset[ ROW ][ MAXIMUM ],
                   data->arguments->subset[ ROW ][ MINIMUM ], data->rows ),
         IN_RANGE( data->arguments->subset[ COLUMN ][ MINIMUM ],
                   1, data->columns ),
         IN_RANGE( data->arguments->subset[ COLUMN ][ MAXIMUM ],
                   data->arguments->subset[ COLUMN ][ MINIMUM ],
                   data->columns ) );

  Arguments* const arguments = data->arguments;
#if TEST_CLIP_GRID_CELLS
  double longitudes[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
  double latitudes[  4 ] = { 0.0, 0.0, 0.0, 0.0 };
  double clipLongitudes[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
  double clipLatitudes[  4 ] = { 0.0, 0.0, 0.0, 0.0 };
  const double west  = arguments->bounds[ LONGITUDE ][ MINIMUM ];
  const double east  = arguments->bounds[ LONGITUDE ][ MAXIMUM ];
  const double south = arguments->bounds[ LATITUDE  ][ MINIMUM ];
  const double north = arguments->bounds[ LATITUDE  ][ MAXIMUM ];
#endif
  int firstRow    = arguments->subset[ ROW    ][ MINIMUM ];
  int lastRow     = arguments->subset[ ROW    ][ MAXIMUM ];
  int firstColumn = arguments->subset[ COLUMN ][ MINIMUM ];
  int lastColumn  = arguments->subset[ COLUMN ][ MAXIMUM ];
  int row = 0;
  int column = 0;
  int intersects = 0;
  int result = 0;

  DEBUG( fprintf( stderr,"boundsSubset(): subset row indices [%d %d] "
                  "column indices [%d %d]\n",
                  firstRow, lastRow, firstColumn, lastColumn ); )

  /*
   * Check if bounds overlaps with each grid cell's bounds.
   * This is efficient and includes cells slightly outside the bounds
   * which is perhaps acceptable.
   * If tighter bounds checking is desired then #define TEST_CLIP_GRID_CELLS 1
   * which adds a more expensive quadrilateral clipping test
   * but only for cells that are not subsumed by (completely inside) bounds.
   */

  /* Find first subset row that has a grid cell that intersects bounds: */

  for ( row = firstRow, intersects = 0;
        AND2( ! intersects, row <= lastRow ); ++row ) {

    for ( column = firstColumn;
          AND2( ! intersects, column <= lastColumn ); ++column ) {
      Bounds cellBounds;
      computeCellBounds( data, row - 1, column - 1, cellBounds );
      intersects = boundsOverlap( (const double (*)[2]) arguments->bounds,
                                  (const double (*)[2]) cellBounds );

      DEBUG( if ( intersects )
          fprintf( stderr, "first    row: intersects cell = [%f %f][%f %f]\n",
                   cellBounds[ 0 ][ 0 ], cellBounds[ 0 ][ 1 ],
                   cellBounds[ 1 ][ 0 ], cellBounds[ 1 ][ 1 ] ); )

#if TEST_CLIP_GRID_CELLS
      if ( AND2( intersects,
                 ! boundsSubsumes( (const double (*)[2]) arguments->bounds,
                                   (const double (*)[2]) cellBounds ) ) ) {
        getCellVertices( data, row - 1, column - 1, longitudes, latitudes );
        intersects = clipPolygon( 0, west, south, east, north, 4,
                                  longitudes, latitudes,
                                  clipLongitudes, clipLatitudes );

#ifdef DEBUGGING
        fprintf( stderr, "bounds [%f %f][%f %f]\n", west, east, south, north );
        fprintf( stderr, "cell [%f %f %f %f], [%f %f %f %f]\n",
                 longitudes[ 0 ], longitudes[ 1 ],
                 longitudes[ 2 ], longitudes[ 3 ],
                 latitudes[ 0 ], latitudes[ 1 ],
                 latitudes[ 2 ], latitudes[ 3 ] );
        fprintf( stderr, "First subset row %d intersects = %d\n",
                 row, intersects );
#endif
      }
#endif
    }
  }

  if ( intersects ) {
    arguments->subset[ ROW ][ MINIMUM ] = firstRow = row - 1;
    result = 1; /* At least one grid cell is within bounds. */

    DEBUG( fprintf( stderr, "firstRow = %d\n", firstRow ); )

    /* Find last subset row that has a grid cell that intersects bounds: */

    for ( row = lastRow, intersects = 0;
          AND2( ! intersects, row > firstRow ); --row ) {

      for ( column = firstColumn;
            AND2( ! intersects, column <= lastColumn ); ++column ) {
        Bounds cellBounds;
        computeCellBounds( data, row - 1, column - 1, cellBounds );
        intersects = boundsOverlap( (const double (*)[2]) arguments->bounds,
                                    (const double (*)[2]) cellBounds );

      DEBUG( if ( intersects )
          fprintf( stderr, "last     row: intersects cell = [%f %f][%f %f]\n",
                   cellBounds[ 0 ][ 0 ], cellBounds[ 0 ][ 1 ],
                   cellBounds[ 1 ][ 0 ], cellBounds[ 1 ][ 1 ] ); )

#if TEST_CLIP_GRID_CELLS
        if ( AND2( intersects,
                   ! boundsSubsumes( (const double (*)[2]) arguments->bounds,
                                     (const double (*)[2]) cellBounds ) ) ) {
          getCellVertices( data, row - 1, column - 1, longitudes, latitudes );
          intersects = clipPolygon( 0, west, south, east, north, 4,
                                    longitudes, latitudes,
                                    clipLongitudes, clipLatitudes );

#ifdef DEBUGGING
          fprintf( stderr, "bounds [%f %f][%f %f]\n", west, east, south, north);
          fprintf( stderr, "cell [%f %f %f %f], [%f %f %f %f]\n",
                   longitudes[ 0 ], longitudes[ 1 ],
                   longitudes[ 2 ], longitudes[ 3 ],
                   latitudes[ 0 ], latitudes[ 1 ],
                   latitudes[ 2 ], latitudes[ 3 ] );
          fprintf( stderr, "Last  subset row %d intersects = %d\n",
                   row, intersects );
#endif
        }
#endif
      }
    }

    arguments->subset[ ROW ][ MAXIMUM ] = lastRow = row + 1;

    DEBUG( fprintf( stderr, "lastRow = %d\n", lastRow ); )

    /* Find first subset column that has a grid cell that intersects bounds*/

    for ( column = firstColumn, intersects = 0;
          AND2( ! intersects, column <= lastColumn ); ++column ) {

      for ( row = firstRow;
            AND2( ! intersects, row <= lastRow ); ++row ) {
        Bounds cellBounds;
        computeCellBounds( data, row - 1, column - 1, cellBounds );
        intersects = boundsOverlap( (const double (*)[2]) arguments->bounds,
                                    (const double (*)[2]) cellBounds );

      DEBUG( if ( intersects )
          fprintf( stderr, "first column: intersects cell = [%f %f][%f %f]\n",
                   cellBounds[ 0 ][ 0 ], cellBounds[ 0 ][ 1 ],
                   cellBounds[ 1 ][ 0 ], cellBounds[ 1 ][ 1 ] ); )

#if TEST_CLIP_GRID_CELLS
        if ( AND2( intersects,
                   ! boundsSubsumes( (const double (*)[2]) arguments->bounds,
                                     (const double (*)[2]) cellBounds ) )) {
          getCellVertices( data, row - 1, column - 1, longitudes, latitudes);
          intersects = clipPolygon( 0, west, south, east, north, 4,
                                    longitudes, latitudes,
                                    clipLongitudes, clipLatitudes );

#ifdef DEBUGGING
          fprintf( stderr, "bounds [%f %f][%f %f]\n", west, east, south, north);
          fprintf( stderr, "cell [%f %f %f %f], [%f %f %f %f]\n",
                   longitudes[ 0 ], longitudes[ 1 ],
                   longitudes[ 2 ], longitudes[ 3 ],
                   latitudes[ 0 ], latitudes[ 1 ],
                   latitudes[ 2 ], latitudes[ 3 ] );
          fprintf( stderr, "First subset column %d intersects = %d\n",
                   column, intersects );
#endif
        }
#endif
      }
    }

    arguments->subset[ COLUMN ][ MINIMUM ] = firstColumn = column - 1;

    DEBUG( fprintf( stderr, "firstColumn = %d\n", firstColumn ); )

    /* Find last subset column that has a grid cell intersecting bounds: */

    for ( column = lastColumn, intersects = 0;
          AND2( ! intersects, column > firstColumn ); --column ) {

      for ( row = firstRow;
            AND2( ! intersects, row <= lastRow ); ++row ) {
        Bounds cellBounds;
        computeCellBounds( data, row - 1, column - 1, cellBounds );
        intersects = boundsOverlap((const double (*)[2]) arguments->bounds,
                                   (const double (*)[2]) cellBounds );

      DEBUG( if ( intersects )
          fprintf( stderr, "last  column: intersects cell = [%f %f][%f %f]\n",
                   cellBounds[ 0 ][ 0 ],cellBounds[ 0 ][ 1 ],
                   cellBounds[ 1 ][ 0 ], cellBounds[ 1 ][ 1 ] ); )

#if TEST_CLIP_GRID_CELLS
        if ( AND2( intersects,
                   ! boundsSubsumes((const double (*)[2]) arguments->bounds,
                                    (const double (*)[2]) cellBounds ) )) {
          getCellVertices(data, row - 1, column - 1, longitudes,latitudes);
          intersects = clipPolygon( 0, west, south, east, north, 4,
                                    longitudes, latitudes,
                                    clipLongitudes, clipLatitudes );

#ifdef DEBUGGING
          fprintf( stderr, "bounds [%f %f][%f %f]\n", west, east, south, north);
          fprintf( stderr, "cell [%f %f %f %f], [%f %f %f %f]\n",
                   longitudes[ 0 ], longitudes[ 1 ],
                   longitudes[ 2 ], longitudes[ 3 ],
                   latitudes[ 0 ], latitudes[ 1 ],
                   latitudes[ 2 ], latitudes[ 3 ] );
          fprintf( stderr, "Last  subset column %d intersects = %d\n",
                   column, intersects );
#endif
        }
#endif
      }
    }

    arguments->subset[ COLUMN ][ MAXIMUM ] = lastColumn = column + 1;

    DEBUG( fprintf( stderr, "lastColumn = %d\n", lastColumn ); )
  }

#ifdef DEBUGGING
  fprintf( stderr, "subset: bounds = [%f %f][%f %f], "
           "rows = [%d %d], columns = [%d %d]\n",
           arguments->bounds[ LONGITUDE ][ MINIMUM ],
           arguments->bounds[ LONGITUDE ][ MAXIMUM ],
           arguments->bounds[ LATITUDE  ][ MINIMUM ],
           arguments->bounds[ LATITUDE  ][ MAXIMUM ],
           firstRow, lastRow, firstColumn, lastColumn );
#endif

  if ( ! result ) {
    fputs( "\nNo data is within the spatial subset.\n", stderr );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND4( IN_RANGE( data->arguments->subset[ ROW ][ MINIMUM ],
                                   1, data->rows ),
                         IN_RANGE( data->arguments->subset[ ROW ][ MAXIMUM ],
                                   data->arguments->subset[ ROW ][ MINIMUM ],
                                   data->rows ),
                         IN_RANGE( data->arguments->subset[ COLUMN ][ MINIMUM],
                                   1, data->columns ),
                         IN_RANGE( data->arguments->subset[ COLUMN ][ MAXIMUM],
                                   data->arguments->subset[ COLUMN ][ MINIMUM],
                                   data->columns ) ) ) );

  return result;
}



/******************************************************************************
PURPOSE: computeCellBounds - Compute grid cell lon-lat bounds.
INPUTS:  Data* const data  data->arguments->bounds, subset.
         const int row     0-based row index.
         const int column  0-based column index.
OUTPUTS: const Data* const data  data->arguments->subset.
RETURNS: int 1 if subset intersects, else 0 and a failure message is printed
         to stderr.
******************************************************************************/

static void computeCellBounds( const Data* const data,
                               const int row,
                               const int column,
                               Bounds bounds ) {

  PRE07( data, data->rows, data->columns, data->longitudes, data->latitudes,
         IN_RANGE( row,    0, data->rows    - 1 ),
         IN_RANGE( column, 0, data->columns - 1 ) );

  const int columns1 = data->columns + 1;
  const int rowOffset = row * columns1;
  const int index1 = rowOffset + column;
  const int index2 = index1 + 1;
  const int index3 = index2 + columns1;
  const int index4 = index1 + columns1;

  CHECK4( IN_RANGE( index1, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1),
          IN_RANGE( index2, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1),
          IN_RANGE( index3, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1),
          IN_RANGE( index4, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1));
  {
    const double longitude1 = data->longitudes[ index1 ];
    const double longitude2 = data->longitudes[ index2 ];
    const double longitude3 = data->longitudes[ index3 ];
    const double longitude4 = data->longitudes[ index4 ];
    const double latitude1  = data->latitudes[  index1 ];
    const double latitude2  = data->latitudes[  index2 ];
    const double latitude3  = data->latitudes[  index3 ];
    const double latitude4  = data->latitudes[  index4 ];
    double minimum = 0.0;
    double maximum = 0.0;

#ifdef DEBUGGING
    if ( ( row == 81 && ( column <= 1 || column >= data->columns - 2 ) ) ||
         OR2( AND2( row <= 1, column <= 1 ),
              AND2( row >= data->rows - 2, column >= data->columns - 2 ) ) )
      fprintf( stderr, "computeCellBounds( row = %d, column = %d ) = "
               "[%d %d %d %d], "
               " = [%f %f %f %f], [%f %f %f %f]\n",
               row, column, index1, index2, index3, index4,
               longitude1, longitude2, longitude3, longitude4,
               latitude1, latitude2, latitude3, latitude4 );
#endif

    /* Compute longitude range: */

    if ( longitude1 < longitude2 ) {
      minimum = longitude1;
      maximum = longitude2;
    } else {
      minimum = longitude2;
      maximum = longitude1;
    }

    if ( longitude3 < minimum ) {
      minimum = longitude3;
    } else if ( longitude3 > maximum ) {
      maximum = longitude3;
    }

    if ( longitude4 < minimum ) {
      minimum = longitude4;
    } else if ( longitude4 > maximum ) {
      maximum = longitude4;
    }

    /*
     * HACK to handle Stereographic case with cells that unproject across the
     * -180/180 line. Truncate the cell so it does not cross the -180/180 line:
     */

    CHECK( minimum <= maximum );

    if ( AND2( data->isProjected, maximum - minimum > 180.0 ) ) {
      maximum = minimum;
      minimum = -179.999;

      if ( minimum > maximum ) {
        const double temp = minimum;
        minimum = maximum;
        maximum = temp;
      }
    }

    CHECK( minimum <= maximum );
    bounds[ LONGITUDE ][ MINIMUM ] = minimum;
    bounds[ LONGITUDE ][ MAXIMUM ] = maximum;

    /* Compute latitude range: */

    if ( latitude1 < latitude2 ) {
      minimum = latitude1;
      maximum = latitude2;
    } else {
      minimum = latitude2;
      maximum = latitude1;
    }

    if ( latitude3 < minimum ) {
      minimum = latitude3;
    } else if ( latitude3 > maximum ) {
      maximum = latitude3;
    }

    if ( latitude4 < minimum ) {
      minimum = latitude4;
    } else if ( latitude4 > maximum ) {
      maximum = latitude4;
    }

    CHECK( minimum <= maximum );
    bounds[ LATITUDE ][ MINIMUM ] = minimum;
    bounds[ LATITUDE ][ MAXIMUM ] = maximum;
  }

#ifdef DEBUGGING
    if ( ( row == 81 && ( column <= 1 || column >= data->columns - 2 ) ) ||
         OR2( AND2( row <= 1, column <= 1 ),
              AND2( row >= data->rows - 2, column >= data->columns - 2 ) ) )
      fprintf( stderr, "  cell bounds = [%f %f][%f %f]\n",
               bounds[ LONGITUDE ][ MINIMUM ],
               bounds[ LONGITUDE ][ MAXIMUM ],
               bounds[ LATITUDE  ][ MINIMUM ],
               bounds[ LATITUDE  ][ MAXIMUM ] );
#endif

  POST0( isValidBounds( (const double (*)[2]) bounds ) );
}



#if TEST_CLIP_GRID_CELLS

/******************************************************************************
PURPOSE: getCellVertices - Get grid cell lon-lat vertices in counter-clockwise
         order.
INPUTS:  Data* const data  data->rows, columns, longitudes, latitudes.
         const int row     0-based row index.
         const int column  0-based column index.
OUTPUTS: const Data* const data  data->arguments->subset.
RETURNS: int 1 if subset intersects, else 0 and a failure message is printed
         to stderr.
******************************************************************************/

static void getCellVertices( const Data* const data,
                             const int row, const int column,
                             double longitudes[ 4 ], double latitudes[ 4 ] ) {

  PRE09( data, data->rows > 0, data->columns > 0,
         data->longitudes, data->latitudes,
         IN_RANGE( row,    0, data->rows   - 1 ),
         IN_RANGE( column, 0, data->columns - 1 ),
         longitudes, latitudes );

  const int columns1 = data->columns + 1;
  const int rowOffset = row * columns1;
  const int index1 = rowOffset + column; /* Counter-clockwise vertex order. */
  const int index2 = index1 + 1;
  const int index3 = index2 + columns1;
  const int index4 = index1 + columns1;

  CHECK4( IN_RANGE( index1, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1),
          IN_RANGE( index2, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1),
          IN_RANGE( index3, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1),
          IN_RANGE( index4, 0, ( data->rows + 1 ) * ( data->columns + 1 ) - 1));

  longitudes[ 0 ] = data->longitudes[ index1 ];
  longitudes[ 1 ] = data->longitudes[ index2 ];
  longitudes[ 2 ] = data->longitudes[ index3 ];
  longitudes[ 3 ] = data->longitudes[ index4 ];
  latitudes[  0 ] = data->latitudes[  index1 ];
  latitudes[  1 ] = data->latitudes[  index2 ];
  latitudes[  2 ] = data->latitudes[  index3 ];
  latitudes[  3 ] = data->latitudes[  index4 ];

  POST04( isValidLongitudeLatitude( longitudes[ 0 ], latitudes[ 0 ] ),
          isValidLongitudeLatitude( longitudes[ 1 ], latitudes[ 1 ] ),
          isValidLongitudeLatitude( longitudes[ 2 ], latitudes[ 2 ] ),
          isValidLongitudeLatitude( longitudes[ 3 ], latitudes[ 3 ] ) );
}

#endif



/******************************************************************************
PURPOSE: computeBounds - Compute overall bounds of grid.
INPUTS:  const size_t point     Number of points.
         const double longitudes[ points ]  Longitudes of points.
         const double latitudes[  points ]  Latitudes of points.
OUTPUTS: Bounds bounds                      Bounds of grid.
******************************************************************************/

static void computeBounds( const size_t points,
                           const double longitudes[],
                           const double latitudes[],
                           Bounds bounds ) {

  PRE05( points, longitudes, latitudes,
         validLongitudesAndLatitudes( points, longitudes, latitudes ),
         bounds );

  double west = longitudes[ 0 ];
  double east = west;
  double south = latitudes[ 0 ];
  double north = south;
  size_t index = 1;

  for ( ; index < points; ++index ) {
    const double longitude = longitudes[ index ];
    const double latitude  = latitudes[  index ];

    if ( longitude < west ) {
      west = longitude;
    } else if ( longitude > east ) {
      east = longitude;
    }

    if ( latitude < south ) {
      south = latitude;
    } else if ( latitude > north ) {
      north = latitude;
    }
  }

  bounds[ LONGITUDE ][ MINIMUM ] = west;
  bounds[ LONGITUDE ][ MAXIMUM ] = east;
  bounds[ LATITUDE  ][ MINIMUM ] = south;
  bounds[ LATITUDE  ][ MAXIMUM ] = north;

  POST0( isValidBounds( (const double (*)[2]) bounds ) );
}



/******************************************************************************
PURPOSE: findTimestampedVariable - Get NetCDF file and variable id and time
         index of named variable at timestamp.
INPUTS:  const Data* const data          Data structure.
         const char* const variableName  Name of variable.
         const int yyyymmddhh            Timestamp to find.
OUTPUTS: int* const variableId           NetCDF Id of variable.
         int* const timestep             0-based time index of variable data.
RETURNS: int NetCDF file id of file containing variable at timestep,
         else 0 and a failure message is printed to stderr.
******************************************************************************/

static int findTimestampedVariable( const Data* const data,
                                    const char* const variableName,
                                    const int yyyymmddhh,
                                    int* const variableId,
                                    int* const timestep ) {

  PRE08( data, data->arguments, isValidArguments( data->arguments ),
         variableName, *variableName, isValidYYYYMMDDHH( yyyymmddhh ),
         variableId, timestep );

  const Arguments* const arguments = data->arguments;
  int fileCount = arguments->fileCount;
  const char* const* fileNames = arguments->fileNames;
  int fileTimeRange[ MAX_FILES ][ 4 ]; /* First/last timestamp of files */
  int found = 0;
  int fileIndex = 0;
  int result  = -1;
  *variableId = -1;
  *timestep = -1;
  memcpy( fileTimeRange, data->fileTimeRange, sizeof fileTimeRange );

  DEBUG( fprintf( stderr, "findTimestampedVariable( %s, %d )\n",
                  variableName, yyyymmddhh ); )

  /* Determine which set of files to search: */

  if ( AND2( arguments->zfFileCount > 0,
             OR3( ! strcmp( variableName, "ZH" ),
                  ! strcmp( variableName, "ZF" ),
                  ! strcmp( variableName, "DENS" ) ) ) ) {
    fileCount = arguments->zfFileCount;
    fileNames = arguments->zfFileNames;
    memcpy( fileTimeRange, data->zfFileTimeRange, sizeof fileTimeRange );
  } else if ( AND2( arguments->wwindFileCount > 0,
                    ! strcmp( variableName, data->wwindVariable ) ) ) {
    fileCount = arguments->wwindFileCount;
    fileNames = arguments->wwindFileNames;
    memcpy( fileTimeRange, data->wwindFileTimeRange, sizeof fileTimeRange );
  }

  DEBUG( fprintf( stderr, "  fileCount = %d %s: [%d %d] ... %s: [%d %d]\n",
                  fileCount, fileNames[ 0 ],
                  fileTimeRange[ 0 ][ MINIMUM ],
                  fileTimeRange[ 0 ][ MAXIMUM ],
                  fileNames[ fileCount - 1 ],
                  fileTimeRange[ fileCount - 1 ][ MINIMUM ],
                  fileTimeRange[ fileCount - 1 ][ MAXIMUM ] ); )

  /* Find file encompassing yyyymmddhh: */

  do {
    const int fileFirstTimesamp = fileTimeRange[ fileIndex ][ MINIMUM ];
    const int fileLastTimesamp  = fileTimeRange[ fileIndex ][ MAXIMUM ];

    if ( IN_RANGE( yyyymmddhh, fileFirstTimesamp, fileLastTimesamp ) ) {
      found = 1;
    } else {
      ++fileIndex;
    }

  } while ( AND2( ! found, fileIndex < fileCount ) );

  /* If found then get timestep index and variable id: */

  if ( found ) {
    const int file = openNetCDFFile( fileNames[ fileIndex ], 'r' );

    if ( file >= 0 ) {
      const int hoursPerTimestep = fileTimeRange[ fileIndex ][ 3 ];
      *timestep = 0;

      if ( hoursPerTimestep > 0 ) {
        const int timesteps = fileTimeRange[ fileIndex ][ 2 ];
        int timestamp = fileTimeRange[ fileIndex ][ MINIMUM ];
        CHECK2( timesteps > 0, isValidYYYYMMDDHH( timestamp ) );

        while ( AND2( *timestep < timesteps, timestamp < yyyymmddhh ) ) {
          timestamp = incrementHours( timestamp, hoursPerTimestep );
          CHECK( isValidYYYYMMDDHH( timestamp ) );
          *timestep += 1;
        }
      }

      *variableId = getNetCDFVariableId( file, variableName );

      if ( *variableId >= 0 ) {
        result = file;
      } else {
        closeNetCDFFile( file );
        result      = -1;
        *variableId = -1;
        *timestep   = -1;
      }
    }
  }

  DEBUG( fprintf( stderr, "findTimestampedVariable() returning file = %d, "
                  "*variableId = %d, *timestep = %d\n",
                  result, *variableId, *timestep ) );

  POST0( IMPLIES_ELSE( result >= 0,
                       AND2( *variableId >= 0, *timestep >= 0 ),
                       AND3( result == -1, *variableId == -1,
                             *timestep == -1 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: writeXDRData - Read the subset of data and write it to
         output as XDR binary data.
INPUTS:  const Data* const data    Data structure to write.
OUTPUTS: FILE* const output        File to write to.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   The output binary data is of the form:
# IEEE-754 32-bit reals data[variables][timesteps][layers][rows][columns]:
******************************************************************************/

static int writeXDRData( Data* const data, FILE* const output ) {

  PRE04( data, data->arguments, isValidArguments( data->arguments ), output );

  const Arguments* const arguments = data->arguments;
  const int integrate = arguments->auxMode == INTEGRATE;
  const int yyyymmddhh1 = arguments->subset[ TIME ][ MINIMUM ];
  const int yyyymmddhh2 = arguments->subset[ TIME ][ MAXIMUM ];
  const int timestepHours =
    arguments->aggregateMode ? 24
    : data->fileTimeRange[ data->skipFileCount ][ 3 ];
  const int subsetLayers =
    COUNT_IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                    arguments->subset[ LAYER ][ MAXIMUM ] );
  const int subsetRows =
    COUNT_IN_RANGE( arguments->subset[ ROW ][ MINIMUM ],
                    arguments->subset[ ROW ][ MAXIMUM ] );
  const int subsetColumns =
    COUNT_IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                    arguments->subset[ COLUMN ][ MAXIMUM ] );
  const size_t subsetCells = subsetLayers * subsetRows * subsetColumns;
  const int writeSubsetCells =
    integrate ? subsetCells / subsetLayers : subsetCells;
  const size_t subsetHours = data->readTimesteps;
  const size_t variableSize = subsetHours * subsetCells;
  const size_t subsetVariables = 1 + integrate * 2; /* var, DENS, ZF. */
  const size_t subsetSize = subsetVariables * variableSize;
  const size_t aggregateAllSize =
    IN3( arguments->aggregateMode, AGGREGATE_MEAN, AGGREGATE_SUM ) ?
         writeSubsetCells : 0;
  float* subsetData = NEW( float, subsetSize + aggregateAllSize * 2 );
  int result = subsetData != 0;

  if ( result ) {
    float* const subsetZF   = integrate ? subsetData + variableSize : 0;
    float* const subsetDENS = integrate ? subsetZF   + variableSize : 0;
    float* const aggregateAllData =
      aggregateAllSize ? subsetData + subsetSize : 0;
    int* const aggregateAllCounts =
      aggregateAllSize ? (int*) aggregateAllData + aggregateAllSize : 0;
    const int variables = arguments->variables + (arguments->auxMode == WIND);
    const int coordinateVariables =
      2 * arguments->lonlat + arguments->elevation;
    int variable = -coordinateVariables;
    const int outputTimesteps = data->outputTimesteps;

    DEBUG( fprintf( stderr, "\n====writeXDRData(): "
                    "arguments->auxMode = %d, timestepHours = %d, "
                    "outputTimesteps = %d\n",
                    arguments->auxMode, timestepHours, outputTimesteps ); )
    do {
      const char* const variableName =
        variable == -3 ? "LONGITUDE" :
        variable == -2 ?
          ( coordinateVariables == 2 ? "LONGITUDE" : "LATITUDE" ) :
        variable == -1 ?
          ( coordinateVariables == 2 ? "LATITUDE" : "ELEVATION" ) :
        ( variable < arguments->variables ?
            arguments->variableNames[ variable ] : data->wwindVariable );
      int yyyymmddhh = yyyymmddhh1;
      int timestep = 0;

      do {

        if ( ! data->isHourlyTimesteps ) {
          yyyymmddhh =
            data->fileTimeRange[ timestep + data->skipFileCount ][ 0 ];
        }

        result =
          readSubset( data, variableName, yyyymmddhh,
                      subsetData, subsetZF, subsetDENS );

        if ( result ) {

          if ( AND3( arguments->aggregateMode,
                     strcmp( variableName, "LONGITUDE" ),
                     strcmp( variableName, "LATITUDE" ) ) ) {
            aggregateData( arguments->aggregateMode, subsetHours,
                           writeSubsetCells, subsetData,
                           aggregateAllData, aggregateAllCounts );
          }

          if ( ! aggregateAllData ) {
            result = writeFloats( writeSubsetCells, subsetData, output );
            ++timestep;
          }
        }

        if ( data->isHourlyTimesteps ) {
          yyyymmddhh = incrementHours( yyyymmddhh, timestepHours );
        }
      } while ( AND3( result, timestep < outputTimesteps,
                     yyyymmddhh <= yyyymmddhh2 ) );

      if ( aggregateAllData ) {
        result = writeFloats( writeSubsetCells, aggregateAllData, output );
      }

      ++variable;
    } while ( AND2( result, variable < variables ) );

    FREE( subsetData );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeXDRHeader - Write ASCII header of XDR format metadata to output.
INPUTS:  const Data* const data    Data structure to write.
OUTPUTS: FILE* const output        File to write to.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES: Header looks like the following:
SUBSET 9.0 CMAQ
M_02_99BRACE
http://www.epa.gov/ttn/scram/,CMAQSubset
2002-04-23T00:00:00-0000
# data dimensions: timesteps variables layers rows columns:
2 4 2 3 4
# subset indices (0-based time, 1-based layer/row/column): first-timestep last-timestep first-layer last-layer first-row last-row first-column last-column:
 0 1 1 2 1 3 1 4
# Variable names:
LONGITUDE LATITUDE ELEVATION DAILY_MAX8_O3
# Variable units:
deg deg m ppmV
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
30 60 40 -100 6.371e+06 6.371e+06
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[3]:
4 3 1.578e+06 -1.27e+06 2000 2000 2 10000 1 0.995 0.99
# IEEE-754 32-bit reals data[variables][timesteps][layers][rows][columns]:
******************************************************************************/

static int writeXDRHeader( const Data* const data, FILE* const output ) {

  PRE04( data, data->arguments, isValidArguments( data->arguments ), output );

  const Arguments* const arguments = data->arguments;
  int file = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );
  int result = file >= 0;

  if ( result ) {
    char gdnam[ NAMLEN3 + 1 ] = "";
    memset( gdnam, 0, sizeof gdnam );
    result = getNetCDFStringAttribute( file, -1, "GDNAM",
                                       sizeof gdnam / sizeof *gdnam, gdnam );
    closeNetCDFFile( file ), file = -1;

    if ( result ) {
      const int integrate = arguments->auxMode == INTEGRATE;
      const int layers =
        integrate ? 1
        : COUNT_IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                          arguments->subset[ LAYER ][ MAXIMUM ] );
      const int rows =
        COUNT_IN_RANGE( arguments->subset[ ROW ][ MINIMUM ],
                        arguments->subset[ ROW ][ MAXIMUM ] );
      const int columns =
        COUNT_IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                        arguments->subset[ COLUMN ][ MAXIMUM ] );
      const int outputVariables =
        arguments->lonlat * 2 + arguments->elevation + arguments->variables +
        ( arguments->auxMode == WIND );
      const int yyyymmddhh = data->yyyymmddhh;
      const int timesteps  = data->outputTimesteps;

      /* Determine index of first output timestep: */

      int firstTimeIndex = 0;
      int fileCount = arguments->fileCount;
      int index = 0;

      for ( ; index < fileCount; ++index ) {
        const int yyyymmddhhFirst = data->fileTimeRange[ index ][ MINIMUM ];
        const int yyyymmddhhLast  = data->fileTimeRange[ index ][ MAXIMUM ];

        if ( yyyymmddhh <= yyyymmddhhLast ) {
          firstTimeIndex = timestepsUntil( yyyymmddhhFirst, yyyymmddhh, 1 );
          index = fileCount; /* Stop looping. */
        }
      }

      result =
        fprintf( output,
                "SUBSET 9.0 CMAQ\n"
                "%s\n"
                "%s\n"
                "%04d-%02d-%02dT%02d:00:00-0000\n"
                "# data dimensions: timesteps variables layers rows columns:\n"
                "%d %d %d %d %d\n"
                "# subset indices (0-based time, 1-based layer/row/column): "
                  "first-timestep last-timestep first-layer last-layer "
                  "first-row last-row first-column last-column:\n"
                "%d %d %d %d %d %d %d %d\n",
                gdnam,
                arguments->note,
                yyyymmddhh / 1000000,
                yyyymmddhh / 10000 % 100,
                yyyymmddhh / 100 % 100,
                yyyymmddhh % 100,
                timesteps,
                outputVariables,
                layers,
                rows,
                columns,
                firstTimeIndex,
                firstTimeIndex + timesteps - 1,
                arguments->subset[ LAYER  ][ MINIMUM ],
                arguments->subset[ LAYER ][integrate ? MINIMUM : MAXIMUM],
                arguments->subset[ ROW    ][ MINIMUM ],
                arguments->subset[ ROW    ][ MAXIMUM ],
                arguments->subset[ COLUMN ][ MINIMUM ],
                arguments->subset[ COLUMN ][ MAXIMUM ] );

      if ( result ) {
        result = writeXDRVariableNamesAndUnits( data, output );

        if ( result ) {
          result = writeXDRProjection( data, output );

          if ( result ) {
            result = writeXDRGrid( data, output );

            if ( result ) {
              result =
                fputs( "# IEEE-754 32-bit reals "
                       "data[variables][timesteps][layers][rows][columns]:\n",
                       output )
                != EOF;
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeXDRProjection - Write part of XDR header for CMAQ
         projection and grid.
INPUTS:  const Data* const data    Data structure to write.
OUTPUTS: FILE* const output        File to write to.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES: CMAQ projection and grid lines look this this:
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
30 60 40 -100 6.371e+06 6.371e+06
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[3]:
4 3 1.578e+06 -1.27e+06 2000 2000 2 10000 1 0.995 0.99
******************************************************************************/

static int writeXDRProjection( const Data* const data, FILE* const output ) {
  PRE04( data, data->arguments, isValidArguments( data->arguments ), output );

  const Arguments* const arguments = data->arguments;
  const int file = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );
  int result = file != -1;

  if ( result ) {
    int gdtyp = 0;
    result = getNetCDFIntAttribute( file, "GDTYP", &gdtyp );

    if ( result ) {

      /* Read and write CMAQ projection parameters: */

      switch ( gdtyp ) {
      case LATGRD3:
        result =
          fprintf( output,
                   "# lonlat projection: major_semiaxis minor_semiaxis\n"
                   "%g %g\n",
                   arguments->ellipsoid[ MAXIMUM ],
                   arguments->ellipsoid[ MINIMUM ] ) > 50;
        break;
      case ALBGRD3:
        /* lint -fallthrough */
      case LAMGRD3:
        {
          double p_alp = 0.0;
          double p_bet = 0.0;
          double xcent = 0.0;
          double ycent = 0.0;
          result =
            AND4( getNetCDFDoubleAttribute( file, "P_ALP", &p_alp ),
                  getNetCDFDoubleAttribute( file, "P_BET", &p_bet ),
                  getNetCDFDoubleAttribute( file, "XCENT", &xcent ),
                  getNetCDFDoubleAttribute( file, "YCENT", &ycent ) );

          if ( result ) {
            result =
              fprintf( output,
                       "# %s projection: lat_1 lat_2 lat_0 lon_0 "
                       "major_semiaxis minor_semiaxis\n"
                       "%g %g %g %g %g %g\n",
                       gdtyp == LAMGRD3 ? "lcc" : "albers",
                       p_alp, p_bet, ycent, xcent,
                       arguments->ellipsoid[ MAXIMUM ],
                       arguments->ellipsoid[ MINIMUM ] ) > 70;
          }
        }
        break;
      case EQMGRD3:
        {
          double xcent = 0.0;
          result = getNetCDFDoubleAttribute( file, "XCENT", &xcent );

          if ( result ) {
            result =
              fprintf( output,
                       "# mercator projection: lon_0 "
                       "major_semiaxis minor_semiaxis\n"
                       "%g %g %g\n",
                       xcent,
                       arguments->ellipsoid[ MAXIMUM ],
                       arguments->ellipsoid[ MINIMUM ] ) > 60;
          }
        }
        break;
      case POLGRD3:
        {
          double p_bet = 0.0;
          double xcent = 0.0;
          double ycent = 0.0;
          result =
            AND3( getNetCDFDoubleAttribute( file, "P_BET", &p_bet ),
                  getNetCDFDoubleAttribute( file, "XCENT", &xcent ),
                  getNetCDFDoubleAttribute( file, "YCENT", &ycent ) );

          if ( result ) {
            result =
              fprintf( output,
                       "# stereographic projection: lat_0 lon_0 lat_sec "
                       "major_semiaxis minor_semiaxis\n"
                       "%g %g %g %g %g\n",
                       ycent, xcent, p_bet,
                       arguments->ellipsoid[ MAXIMUM ],
                       arguments->ellipsoid[ MINIMUM ] ) > 70;
          }
        }
        break;
      default:
        fprintf( stderr, "\nUnsupported projection type GDTYP = %d.\n", gdtyp);
        result = 0;
        break;
      }
    }

    closeNetCDFFile( file );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeXDRGrid - Write part of XDR header for CMAQ subset grid.
INPUTS:  const Data* const data    Data structure to write.
OUTPUTS: FILE* const output        File to write to.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES: CMAQ projection and grid lines look this this:
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[1+layers]: ...
******************************************************************************/

static int writeXDRGrid( const Data* const data, FILE* const output ) {
  PRE04( data, data->arguments, isValidArguments( data->arguments ), output );

  const Arguments* const arguments = data->arguments;
  const int file = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );
  int result = file != -1;

  if ( result ) {
    double xorig = 0.0;
    double yorig = 0.0;
    double xcell = 0.0;
    double ycell = 0.0;
    float vgtop = 0.0;
    int vgtyp = 0;
    int nlays = 0;
    float vglvls[ MXLAYS3 + 1 ];
    memset( vglvls, 0, sizeof vglvls );
    result =
      AND9( getNetCDFDoubleAttribute( file, "XORIG", &xorig ),
            getNetCDFDoubleAttribute( file, "YORIG", &yorig ),
            getNetCDFDoubleAttribute( file, "XCELL", &xcell ),
            getNetCDFDoubleAttribute( file, "YCELL", &ycell ),
            getNetCDFFloatAttribute( file, "VGTOP", &vgtop ),
            getNetCDFIntAttribute( file, "VGTYP", &vgtyp ),
            getNetCDFIntAttribute( file, "NLAYS", &nlays ),
            IN_RANGE( nlays, 1, MXLAYS3 ),
            getNetCDFFloatArrayAttribute( file, "VGLVLS", nlays + 1, vglvls ));

    if ( result ) {
      const int layers  = data->layers;
      const int rows    = data->rows;
      const int columns = data->columns;
      int level = 0;
      checkAndFixVerticalGridParameters( layers, &vgtyp, &vgtop, vglvls );
      result =
        fprintf( output, "# Grid: ncols nrows xorig yorig xcell ycell "
                         "vgtyp vgtop vglvls[%d]:\n"
                         "%d %d %g %g %g %g %d %g",
                         1 + layers,
                         columns, rows, xorig, yorig, xcell, ycell,
                         vgtyp, vgtop ) > 60;

      for ( ; AND2( result, level <= layers ); ++level ) {
        result = fprintf( output, " %g", vglvls[ level ] ) > 1;
      }

      if ( result ) {
        result = fputc( '\n', output ) != EOF;
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeXDRVariableNamesAndUnits - Write part of XDR header for variable
         names and units.
INPUTS:  const Data* const data    Data structure to write.
OUTPUTS: FILE* const output        File to write to.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int writeXDRVariableNamesAndUnits( const Data* const data,
                                          FILE* const output ) {

  PRE04( data, data->arguments, isValidArguments( data->arguments ), output );

  int result = fputs( "# Variable names:\n", output ) != EOF;

  if ( result ) {
    const Arguments* const arguments = data->arguments;

    if ( arguments->lonlat ) {
      result = fputs( "LONGITUDE LATITUDE ", output ) != EOF;
    }

    if ( result ) {
      const int aggregateMode = arguments->aggregateMode;
      const int auxMode       = arguments->auxMode;
      const int variables = arguments->variables;
      int variable = 0;

      if ( arguments->elevation ) {
        result = fputs( "ELEVATION ", output ) != EOF;
      }

      for ( ; AND2( result, variable < variables ); ++variable ) {
        const char* const prefix =
        aggregateMode   == AGGREGATE_DAILY_MEAN ? "DAILY_MEAN_"
        : aggregateMode == AGGREGATE_DAILY_MAX  ? "DAILY_MAX_"
        : aggregateMode == AGGREGATE_DAILY_MAX8 ? "DAILY_MAX8_"
        : "";
        const char* const variableName =
          outputVariableName( arguments->variableNames[ variable ] );

        if ( AND2( auxMode == WIND, ! strcmp( variableName, "VWIND" ) ) ) {
          result =
            fprintf( output, "%sVWIND %sWWIND", prefix, prefix ) > 10;
        } else {
          const char* const delimiter = variable + 1 < variables ? " " : "";
          result =
            fprintf( output, "%s%s%s", prefix, variableName, delimiter ) > 0;
        }
      }

      if ( result ) {
        result = fputs( "\n# Variable units:\n", output ) != EOF;

        if ( result ) {

          if ( arguments->lonlat ) {
            result = fputs( "deg deg ", output ) != EOF;
          }

          if ( result ) {

            if ( arguments->elevation ) {
              result = fputs( "m ", output ) != EOF;
            }

            for ( variable = 0; AND2( result, variable < variables );
                  ++variable ) {
              const char* const variableName =
                arguments->variableNames[ variable ];
              const char* const variableUnits = data->variableUnits[ variable];
              const char* const units =
                auxMode == INTEGRATE ? "molecules/cm2"
                : outputVariableUnits( variableName, variableUnits, 1 );
              const char* const delimiter =
                variable + 1 < variables ? " " : "";

              if ( AND2( auxMode == WIND, ! strcmp( variableName, "VWIND" ))) {
                result = /* Add units for WWIND: */
                  fprintf( output, "%s %s%s", units, units, delimiter ) > 0;
              } else {
                result = fprintf( output, "%s%s", units, delimiter ) > 0;
              }
            }

            if ( result ) {
              result = fputc( '\n', output ) != EOF;
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: outputVariableName - Return possibly changed output variable name.
INPUTS:  const char* const variableName  Name of variable to lookup.
RETURNS: const char* name of variable to output.
******************************************************************************/

static const char* outputVariableName( const char* const variableName ) {
  PRE02( variableName, *variableName );
  const size_t count =
    edit ? sizeof variableMetadata / sizeof *variableMetadata : 0;
  const char* result = variableName;
  size_t index = 0;

  for ( ; index < count; ++index ) {
    const VariableMetadata* const entry = variableMetadata + index;

    if ( AND2( entry->name, ! strcmp( entry->name, variableName ) ) ) {

      if ( entry->newName ) {
        result = entry->newName;
      }

      index = count;
    }
  }

  POST02( result, *result );
  return result;
}



/******************************************************************************
PURPOSE: outputVariableUnits - Return possibly changed output variable units.
INPUTS:  const char* const variableName   Name of variable to lookup.
         const char* const variableUnits  Units of variable to lookup.
         const int         convertSpaces  Convert spaces to _?
RETURNS: const char* units of variable to output.
NOTES:   Not thread-safe since it may return a pointer to modified static
         storage.
******************************************************************************/

static const char* outputVariableUnits( const char* const variableName,
                                        const char* const variableUnits,
                                        const int convertSpaces ) {
  PRE04( variableName, *variableName, variableUnits, IS_BOOL( convertSpaces ));
  const size_t count =
    edit ? sizeof variableMetadata / sizeof *variableMetadata : 0;
  const char* result = variableUnits;
  size_t index = 0;

  for ( ; index < count; ++index ) {
    const VariableMetadata* const entry = variableMetadata + index;

    if ( AND2( entry->name, ! strcmp( entry->name, variableName ) ) ) {

      if ( entry->units ) {
        result = entry->units;
      }

      index = count;
    }
  }

  if ( convertSpaces ) {
    static char staticResult[ NAMLEN3 + 1 ] = "";
    memset( staticResult, 0, sizeof staticResult );

    /* Convert ' ' to '_': */

    for ( index = 0; index < NAMLEN3; ++index ) {
      char c = result[ index ];

      if ( c == ' ' ) {
        c = '_';
      }

      staticResult[ index ] = c;
    }

    /* Trim trailing '_': */

    {
      int done = 0;

      while ( AND2( index--, done == 0 ) ) {

        if ( staticResult[ index ] == '_' ) {
          staticResult[ index ] = '\0';
        } else {
          done = 1;
        }
      }
    }

    result = staticResult;
  }

  POST03( result, *result, IMPLIES( convertSpaces, ! strchr( result, ' ' ) ) );
  return result;
}



/******************************************************************************
PURPOSE: outputVariableDescription - Return possibly changed output variable
         description.
INPUTS:  const char* const variableName  Name of variable to lookup.
         const char* const variableDescription  Description of variable.
RETURNS: const char* possibly modified description of variable to output.
******************************************************************************/

static const char* outputVariableDescription( const char* const variableName,
                                              const char* const
                                                variableDescription ) {
  PRE04( variableName, *variableName,
         variableDescription, *variableDescription );
  const size_t count =
    edit ? sizeof variableMetadata / sizeof *variableMetadata : 0;
  const char* result = variableDescription;
  size_t index = 0;

  for ( ; index < count; ++index ) {
    const VariableMetadata* const entry = variableMetadata + index;

    if ( AND2( entry->name, ! strcmp( entry->name, variableName ) ) ) {
      result = entry->description;
      index = count;
    }
  }

  POST02( result, *result );
  return result;
}



/******************************************************************************
PURPOSE: readSubset - Read given timestamp of layer/row/column subset of
         variable.
INPUTS:  Data* const data  Data that specify subset.
         const char* const variableName  Name of variable to read.
         const int yyyymmddhh            Timestamp to read.
OUTPUTS: float subsetData[ data->readTimesteps * subsetLayers * subsetRows *
                           subsetColumns ]   Subset variable data read.
         float subsetZF[]     0 or subset ZF data read.
         float subsetDENS[]   0 or subset DENS data read.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readSubset( Data* const data,
                       const char* const variableName,
                       const int yyyymmddhh,
                       float subsetData[],
                       float subsetZF[],
                       float subsetDENS[] ) {

  PRE011( data, data->arguments, isValidArguments( data->arguments ),
          variableName, *variableName, isValidYYYYMMDDHH( yyyymmddhh ),
          data->readTimesteps > 0,
          IMPLIES( AND2( data->arguments->lonlat,
                         OR4( ! strcmp( variableName, "LONGITUDE" ),
                              ! strcmp( variableName, "longitude" ),
                              ! strcmp( variableName, "LATITUDE" ),
                              ! strcmp( variableName, "latitude" ) ) ),
                   AND2( data->longitudes, data->latitudes ) ),
          IMPLIES( AND2( data->arguments->elevation,
                         OR2( ! strcmp( variableName, "ELEVATION" ),
                              ! strcmp( variableName, "elevation" ) ) ),
                   OR2( data->heights, data->elevations ) ),
          IMPLIES( data->arguments->auxMode == INTEGRATE,
                   AND2( subsetZF, subsetDENS ) ),
          subsetData );

  const Arguments* const arguments = data->arguments;
  int result = 0;

  DEBUG( fprintf( stderr, "readSubset( data, %s, %d, "
                  "subsetData = %p, subsetZF = %p, subsetDENS = %p )\n",
                  variableName, yyyymmddhh,
                  subsetData, subsetZF, subsetDENS ); )

  if ( OR4( ! strcmp( variableName, "LONGITUDE" ),
            ! strcmp( variableName, "longitude" ),
            ! strcmp( variableName, "LATITUDE" ),
            ! strcmp( variableName, "latitude" ) ) ) {
    copySubsetCoordinates( data, variableName, subsetData );
    result = 1;
  } else if ( OR2( ! strcmp( variableName, "ELEVATION" ),
                   ! strcmp( variableName, "elevation" ) ) ) {
    result = 1;

    if ( arguments->zfFileCount > 0 ) { /* Read/expand layer 1 subset of ZH: */
      result = readSubsetZH( data, yyyymmddhh, subsetData );
    }

    if ( result ) {
      copySubsetElevations( data, subsetData );
    }

  } else if ( AND3( arguments->auxMode == WIND,
                    arguments->wwindFileCount > 0,
                    ! strcmp( variableName, data->wwindVariable ) ) ) {

    result = readSubsetWWIND( data, yyyymmddhh, subsetData );
  } else {
    result = readSubsetVariable( data, variableName, yyyymmddhh, subsetData );

    if ( AND2( result, arguments->auxMode == INTEGRATE ) ) {
      result = readSubsetZFAndDENS( data, yyyymmddhh, subsetZF, subsetDENS );

      if ( result ) {
        const int layer0  = arguments->subset[ LAYER  ][ MINIMUM ] - 1;
        const int layer1  = arguments->subset[ LAYER  ][ MAXIMUM ] - 1;
        const int row0    = arguments->subset[ ROW    ][ MINIMUM ] - 1;
        const int row1    = arguments->subset[ ROW    ][ MAXIMUM ] - 1;
        const int column0 = arguments->subset[ COLUMN ][ MINIMUM ] - 1;
        const int column1 = arguments->subset[ COLUMN ][ MAXIMUM ] - 1;
        const int subsetLayers  = COUNT_IN_RANGE( layer0, layer1 );
        const int subsetRows    = COUNT_IN_RANGE( row0, row1 );
        const int subsetColumns = COUNT_IN_RANGE( column0, column1 );
        const int subsetTimesteps  = data->readTimesteps;

        integrateLayers( subsetTimesteps, subsetLayers, subsetRows,
                         subsetColumns, subsetZF, subsetDENS, subsetData );
      }
    }
  }

  DEBUG( fprintf( stderr, "  readSubset(): returning %d\n", result ); )

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readSubsetZH - Read given timestamp of layer/row/column subset
         of ZH data from ZF file.
INPUTS:  Data* const data  Data that specify subset.
         const char* const variableName  Name of variable to read.
         const int yyyymmddhh            Timestamp to read.
OUTPUTS: float subsetData[]   Single timestep of subset elevation data.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readSubsetZH( Data* const data, const int yyyymmddhh,
                         float subsetData[] ) {

  PRE08( data, data->arguments, isValidArguments( data->arguments ),
         data->arguments->zfFileCount > 0,
         isValidYYYYMMDDHH( yyyymmddhh ),
         data->arguments->elevation,
         OR2( data->heights, data->elevations ),
         subsetData );

  const Arguments* const arguments = data->arguments;
  const int layer0  = arguments->subset[ LAYER  ][ MINIMUM ] - 1;
  const int layer1  = arguments->subset[ LAYER  ][ MAXIMUM ] - 1;
  const int row0    = arguments->subset[ ROW    ][ MINIMUM ] - 1;
  const int row1    = arguments->subset[ ROW    ][ MAXIMUM ] - 1;
  const int column0 = arguments->subset[ COLUMN ][ MINIMUM ] - 1;
  const int column1 = arguments->subset[ COLUMN ][ MAXIMUM ] - 1;
  const size_t subsetLayers  = COUNT_IN_RANGE( layer0, layer1 );
  const size_t subsetRows    = COUNT_IN_RANGE( row0, row1 );
  const size_t subsetColumns = COUNT_IN_RANGE( column0, column1 );
  int variableId = -1;
  int timestep = -1;
  int file = 0;
  int result = 0;

  DEBUG( fprintf( stderr, "  readSubsetZH(): "
                  "timestep = %d, layer0 = %d, layer1 = %d, "
                  "row0 = %d, row1 = %d, column0 = %d, column1 = %d\n",
                  timestep, layer0, layer1, row0, row1, column0, column1 ); )

  file =
    findTimestampedVariable( data, "ZH", yyyymmddhh, &variableId, &timestep );
  result = file != -1;

  if ( result ) {
    int dims[ 4 ] = { 0, 0, 0, 0 };
    result =
      AND4( getM3IOVariableDimensions( file, "ZH", dims ),
            dims[ LAYER ] >= data->layers,
            OR2( dims[ ROW    ] == data->rows,
                 dims[ ROW    ] == data->rows - 1 ),
            OR2( dims[ COLUMN ] == data->columns,
                 dims[ COLUMN ] == data->columns - 1 ) );

    if ( result ) {
      int zhRow0          = row0;
      int zhRow1          = row1;
      int zhColumn0       = column0;
      int zhColumn1       = column1;
      int zhSubsetRows    = subsetRows;
      int zhSubsetColumns = subsetColumns;
      int expandRow = 0;
      int expandColumn = 0;

      if ( zhRow1 + 1 > dims[ ROW ] ) {
        --zhRow1;
        --zhRow0;

        if ( zhRow0 < 0 ) {
          ++zhRow0;
        }

        zhSubsetRows = COUNT_IN_RANGE( zhRow0, zhRow1 );
        expandRow = zhSubsetRows < subsetRows;
      }

      if ( zhColumn1 + 1 > dims[ COLUMN ] ) {
        --zhColumn1;
        --zhColumn0;

        if ( zhColumn0 < 0 ) {
          ++zhColumn0;
        }

        zhSubsetColumns = COUNT_IN_RANGE( zhColumn0, zhColumn1 );
        expandColumn = zhSubsetColumns < subsetColumns;
      }

      CHECK4( zhSubsetRows > 0,
              zhSubsetRows == subsetRows - expandRow,
              zhSubsetColumns > 0,
              zhSubsetColumns == subsetColumns - expandColumn );

      result =
        readM3IOVariable( file, variableId, timestep, timestep,
                          layer0, layer1,
                          zhRow0, zhRow1,
                          zhColumn0, zhColumn1,
                          subsetData );

      if ( result ) {

        if ( OR2( expandRow, expandColumn ) ) {
          result =
            expandSubsetData( subsetLayers, zhSubsetRows, zhSubsetColumns,
                              expandRow, expandColumn, subsetData );
        }
      }
    }

    closeNetCDFFile( file );
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read matched ZH data for timestamp %d.\n",
             yyyymmddhh );
  }

  DEBUG( fprintf( stderr, "  readSubsetZH(): returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readSubsetVariable - Read subset timesteps/layers/rows/columns
         data starting at the given timestamp.
INPUTS:  Data* const data  Data that specify subset.
         const char* const variableName  Name of variable to read.
         const int yyyymmddhh            Timestamp to start reading.
OUTPUTS: float subsetData[ data->readTimesteps * subsetLayers * subsetRows *
                           subsetColumns ]   Subset data read.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readSubsetVariable( Data* const data,
                               const char* const variableName,
                               const int yyyymmddhh,
                               float subsetData[] ) {

  PRE08( data, data->arguments, isValidArguments( data->arguments ),
         data->readTimesteps > 0, variableName, *variableName,
         isValidYYYYMMDDHH( yyyymmddhh ),
         subsetData );

  const Arguments* const arguments = data->arguments;
  const int layer0  = arguments->subset[ LAYER  ][ MINIMUM ] - 1;
  const int layer1  = arguments->subset[ LAYER  ][ MAXIMUM ] - 1;
  const int row0    = arguments->subset[ ROW    ][ MINIMUM ] - 1;
  const int row1    = arguments->subset[ ROW    ][ MAXIMUM ] - 1;
  const int column0 = arguments->subset[ COLUMN ][ MINIMUM ] - 1;
  const int column1 = arguments->subset[ COLUMN ][ MAXIMUM ] - 1;
  const int subsetLayers  = COUNT_IN_RANGE( layer0, layer1 );
  const int subsetRows    = COUNT_IN_RANGE( row0, row1 );
  const int subsetColumns = COUNT_IN_RANGE( column0, column1 );
  const int subsetTimesteps  = data->readTimesteps;
  const int hoursPerTimestep = data->fileTimeRange[ 0 ][ 3 ];
  int timestamp = yyyymmddhh;
  int subsetTimestep = 0;
  int result = 0;
  const size_t timestepSubsetSize = subsetLayers * subsetRows * subsetColumns;
  size_t timestepOffset = 0;

  do {
    int variableId = -1;
    int fileTimestep = -1;
    const int file =
      findTimestampedVariable( data, variableName, timestamp,
                               &variableId, &fileTimestep );
    result = file != -1;

    DEBUG( fprintf( stderr, "  readSubsetVariable( %s ): "
                    "fileTimestep = %d, layer0 = %d, layer1 = %d, "
                    "row0 = %d, row1 = %d, column0 = %d, column1 = %d, "
                    "timestepOffset = %lu\n",
                    variableName,
                    fileTimestep, layer0, layer1,
                    row0, row1, column0, column1,
                    timestepOffset ); )

    if ( result ) {
      result =
        readM3IOVariable( file, variableId, fileTimestep, fileTimestep,
                          layer0, layer1, row0, row1, column0, column1,
                          subsetData + timestepOffset );
    }

    if ( AND2( result, arguments->auxMode == INTEGRATE ) ) {

      /* Check if units need to be converted. */

      int type = 0;
      int rank = 0;
      int dims[ 32 ];
      char units[ 256 ] = "";
      memset( units, 0, sizeof units );
      result =
        getNetCDFVariableInfo( file, variableId, 0, &type, &rank,
                               dims, units, 0 );

      if ( result ) {

        if ( ! strcmp( units, "ppbV" ) ) { /* Convert to ppmV. */
          size_t count = timestepSubsetSize;
          float* p = subsetData + timestepOffset;

          while ( count-- ) {
            *p++ *= 1e-3;
          }
        } else if ( strcmp( units, "ppmV" ) ) {
          result = 0;
          fprintf( stderr, "\nInvalid units for integration: %s.\n",
                   units );
        }
      }
    }

    closeNetCDFFile( file );
    timestepOffset += timestepSubsetSize;
    timestamp = incrementHours( timestamp, hoursPerTimestep );
    ++subsetTimestep;
  } while ( AND2( result, subsetTimestep < subsetTimesteps ) );

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read %d hours of %s data starting at "
             "timestamp %d.\n",
             subsetTimesteps, variableName, yyyymmddhh );
  }

  DEBUG( fprintf( stderr, "  readSubsetVariable(): returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readSubsetWWIND - Read subset timesteps/layers/rows/columns WWIND
         data (from a METCRO3D file) starting at the given timestamp.
INPUTS:  Data* const data  Data that specify subset.
         const int yyyymmddhh            Timestamp to start reading.
OUTPUTS: float subsetData[ data->readTimesteps * subsetLayers * subsetRows *
                           subsetColumns ]   Subset WWIND data read.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readSubsetWWIND( Data* const data, const int yyyymmddhh,
                            float subsetData[] ) {

  PRE08( data, data->arguments, isValidArguments( data->arguments ),
         data->arguments->wwindFileCount > 0,
         isValidYYYYMMDDHH( yyyymmddhh ), data->readTimesteps > 0,
         data->arguments->auxMode == WIND,
         subsetData );

  const Arguments* const arguments = data->arguments;
  const int layer0  = arguments->subset[ LAYER  ][ MINIMUM ] - 1;
  const int layer1  = arguments->subset[ LAYER  ][ MAXIMUM ] - 1;
  const int row0    = arguments->subset[ ROW    ][ MINIMUM ] - 1;
  const int row1    = arguments->subset[ ROW    ][ MAXIMUM ] - 1;
  const int column0 = arguments->subset[ COLUMN ][ MINIMUM ] - 1;
  const int column1 = arguments->subset[ COLUMN ][ MAXIMUM ] - 1;
  const int subsetLayers  = COUNT_IN_RANGE( layer0, layer1 );
  const int subsetRows    = COUNT_IN_RANGE( row0, row1 );
  const int subsetColumns = COUNT_IN_RANGE( column0, column1 );
  const int subsetTimesteps  = data->readTimesteps;
  const int hoursPerTimestep = data->wwindFileTimeRange[ 0 ][ 3 ];
  int timestamp = yyyymmddhh;
  int subsetTimestep = 0;
  int result = 0;
  const size_t timestepSubsetSize = subsetLayers * subsetRows * subsetColumns;
  size_t timestepOffset = 0;

  do {
    int variableId = -1;
    int fileTimestep = -1;
    const int file =
      findTimestampedVariable( data, data->wwindVariable, timestamp,
                               &variableId, &fileTimestep );
    result = file != -1;

    DEBUG( fprintf( stderr, "  readSubsetWWIND(): "
                    "fileTimestep = %d, layer0 = %d, layer1 = %d, "
                    "row0 = %d, row1 = %d, column0 = %d, column1 = %d, "
                   "timestepOffset = %lu\n",
                    fileTimestep, layer0, layer1,
                    row0, row1, column0, column1,
                    timestepOffset ); )

    if ( result ) {

      /*
       * If WWIND is from METCRO3D file but data files are METDOT3D then
       * we must expand last row/column edges to match the input data since a
       * CRO grid has one less row and column than the corresponding DOT grid.
       */

      int dims[ 4 ] = { 0, 0, 0, 0 };
      result =
        AND4( getM3IOVariableDimensions( file, data->wwindVariable, dims ),
              dims[ LAYER ] == data->layers,
              OR2( dims[ ROW    ] == data->rows,
                   dims[ ROW    ] == data->rows - 1 ),
              OR2( dims[ COLUMN ] == data->columns,
                   dims[ COLUMN ] == data->columns - 1 ) );

      if ( result ) {
        int wwindRow0          = row0;
        int wwindRow1          = row1;
        int wwindColumn0       = column0;
        int wwindColumn1       = column1;
        int wwindSubsetRows    = subsetRows;
        int wwindSubsetColumns = subsetColumns;
        int expandRow = 0;
        int expandColumn = 0;

        if ( wwindRow1 + 1 > dims[ ROW ] ) {
          --wwindRow1;
          --wwindRow0;

          if ( wwindRow0 < 0 ) {
            ++wwindRow0;
          }

          wwindSubsetRows = COUNT_IN_RANGE( wwindRow0, wwindRow1 );
          expandRow = wwindSubsetRows < subsetRows;
        }

        if ( wwindColumn1 + 1 > dims[ COLUMN ] ) {
          --wwindColumn1;
          --wwindColumn0;

          if ( wwindColumn0 < 0 ) {
            ++wwindColumn0;
          }

          wwindSubsetColumns = COUNT_IN_RANGE( wwindColumn0, wwindColumn1 );
          expandColumn = wwindSubsetColumns < subsetColumns;
        }

        CHECK4( wwindSubsetRows > 0,
                wwindSubsetRows == subsetRows - expandRow,
                wwindSubsetColumns > 0,
                wwindSubsetColumns == subsetColumns - expandColumn );

        result =
          readM3IOVariable( file, variableId, fileTimestep, fileTimestep,
                            layer0, layer1,
                            wwindRow0, wwindRow1,
                            wwindColumn0, wwindColumn1,
                            subsetData + timestepOffset );

        if ( result ) {

          if ( OR2( expandRow, expandColumn ) ) {
            result =
              expandSubsetData( subsetLayers,
                                wwindSubsetRows, wwindSubsetColumns,
                                expandRow, expandColumn,
                                subsetData + timestepOffset );
          }
        }
      }
    }

    closeNetCDFFile( file );
    timestepOffset += timestepSubsetSize;
    timestamp = incrementHours( timestamp, hoursPerTimestep );
    ++subsetTimestep;
  } while ( AND2( result, subsetTimestep < subsetTimesteps ) );

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read %d hours of WWIND data starting at "
             "timestamp %d.\n",
             subsetTimesteps, yyyymmddhh );
  }

#if ZERO_BAD_WWIND

  if ( AND2( result, ! strcmp( data->wwindVariable, "W_VEL" ) ) ) {
    long long index = 0; /* OpenMP requires a signed integer loop index. */
    long long count = data->readTimesteps;
    count *= subsetLayers;
    count *= subsetRows;
    count *= subsetColumns;

#pragma omp parallel for

    for ( index = 0; index < count ; ++index ) {

      if ( ! IS_VALID_VALUE( subsetData[ index ] ) ) {
        subsetData[ index ] = 0.0;
      }
    }
  }

#endif

  DEBUG( fprintf( stderr, "  readSubsetWWIND(): returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readSubsetZFAndDENS - Read given timestamp of layer/row/column
         subset of integration variables - ZF and DENS from METCRO3D file.
INPUTS:  Data* const data  Data that specify subset.
         const int yyyymmddhh            Timestamp to read.
OUTPUTS: float subsetZF[ data->readTimesteps * subsetLayers * subsetRows *
                           subsetColumns ]   Subset ZF data read.
         float subsetDENS[]   Subset DENS data.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readSubsetZFAndDENS( Data* const data,
                                const int yyyymmddhh,
                                float subsetZF[],
                                float subsetDENS[] ) {

  PRE08( data, data->arguments, isValidArguments( data->arguments ),
         data->readTimesteps > 0, isValidYYYYMMDDHH( yyyymmddhh ),
         data->arguments->auxMode == INTEGRATE,
         subsetZF, subsetDENS );

  const Arguments* const arguments = data->arguments;
  const int layer0  = arguments->subset[ LAYER  ][ MINIMUM ] - 1;
  const int layer1  = arguments->subset[ LAYER  ][ MAXIMUM ] - 1;
  const int row0    = arguments->subset[ ROW    ][ MINIMUM ] - 1;
  const int row1    = arguments->subset[ ROW    ][ MAXIMUM ] - 1;
  const int column0 = arguments->subset[ COLUMN ][ MINIMUM ] - 1;
  const int column1 = arguments->subset[ COLUMN ][ MAXIMUM ] - 1;
  const int subsetLayers  = COUNT_IN_RANGE( layer0, layer1 );
  const int subsetRows    = COUNT_IN_RANGE( row0, row1 );
  const int subsetColumns = COUNT_IN_RANGE( column0, column1 );
  const int subsetTimesteps  = data->readTimesteps;
  const int hoursPerTimestep = data->zfFileTimeRange[ 0 ][ 3 ];
  int timestamp = yyyymmddhh;
  int subsetTimestep = 0;
  int result = 0;
  const size_t timestepSubsetSize = subsetLayers * subsetRows * subsetColumns;
  size_t timestepOffset = 0;

  do {
    int variableId = -1;
    int fileTimestep = -1;
    int file =
      findTimestampedVariable( data, "ZF", timestamp,
                               &variableId, &fileTimestep );
    result = file != -1;

    DEBUG( fprintf( stderr, "  readSubsetZFAndDENS( ZF ): "
                    "fileTimestep = %d, layer0 = %d, layer1 = %d, "
                    "row0 = %d, row1 = %d, column0 = %d, column1 = %d, "
                   "timestepOffset = %lu\n",
                    fileTimestep, layer0, layer1,
                    row0, row1, column0, column1,
                    timestepOffset); )

    if ( result ) {
      result =
        readM3IOVariable( file, variableId, fileTimestep, fileTimestep,
                          layer0, layer1, row0, row1, column0, column1,
                          subsetZF + timestepOffset );
      closeNetCDFFile( file ), file = -1;

      if ( result ) {
        file =
          findTimestampedVariable( data, "DENS", timestamp,
                                   &variableId, &fileTimestep );
        result = file != -1;

        DEBUG( fprintf( stderr, "  readSubsetZFAndDENS( DENS ): "
                        "fileTimestep = %d, layer0 = %d, layer1 = %d, "
                        "row0 = %d, row1 = %d, column0 = %d, column1 = %d, "
                       "timestepOffset = %lu\n",
                        fileTimestep, layer0, layer1,
                        row0, row1, column0, column1,
                        timestepOffset ); )

        if ( result ) {
          result =
            readM3IOVariable( file, variableId, fileTimestep, fileTimestep,
                              layer0, layer1, row0, row1, column0, column1,
                              subsetDENS + timestepOffset );
          closeNetCDFFile( file ), file = -1;
        }
      }
    }

    timestepOffset += timestepSubsetSize;
    timestamp = incrementHours( timestamp, hoursPerTimestep );
    ++subsetTimestep;
  } while ( AND2( result, subsetTimestep < subsetTimesteps ) );

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read %d hours of ZF/DENS data starting at "
             "timestamp %d.\n",
             subsetTimesteps, yyyymmddhh );
  }

  DEBUG( fprintf( stderr, "  readSubsetZFAndDENS(): returning %d\n", result);)
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: copySubsetCoordinates - Copy longitude or latitude coordinates to
         all subset timesteps/layers/rows/columns.
INPUTS:  Data* const data  Data that specify subset.
         const char* const variableName  Name of variable to read.
OUTPUTS: float subsetData[ data->readTimesteps * subsetLayers * subsetRows *
                           subsetColumns ]  Subset variable data.
******************************************************************************/

static void copySubsetCoordinates( const Data* const data,
                                   const char* const variableName,
                                   float subsetData[] ) {

  PRE010( data, data->arguments, isValidArguments( data->arguments ),
          data->readTimesteps > 0,
          variableName,
          OR4( ! strcmp( variableName, "LONGITUDE" ),
               ! strcmp( variableName, "longitude" ),
               ! strcmp( variableName, "LATITUDE" ),
               ! strcmp( variableName, "latitude" ) ),
          data->columns, data->longitudes, data->latitudes, subsetData );

  const Arguments* const arguments = data->arguments;
  const int integrate = arguments->auxMode == INTEGRATE;
  const int subsetLayers =
    integrate ? 1 :
    COUNT_IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                    arguments->subset[ LAYER ][ MAXIMUM ] );
  const int firstRow      = arguments->subset[ ROW    ][ MINIMUM ];
  const int lastRow       = arguments->subset[ ROW    ][ MAXIMUM ];
  const int firstColumn   = arguments->subset[ COLUMN ][ MINIMUM ];
  const int lastColumn    = arguments->subset[ COLUMN ][ MAXIMUM ];
  const int subsetRows    = COUNT_IN_RANGE( firstRow, lastRow );
  const int subsetColumns = COUNT_IN_RANGE( firstColumn, lastColumn );
  const size_t subsetCount = subsetRows * subsetColumns;
  const double* const input =
    OR2( ! strcmp( variableName, "LONGITUDE" ),
         ! strcmp( variableName, "longitude" ) ) ? data->longitudes
    : data->latitudes;
  float* output = subsetData;
  const int columns1 = data->columns + 1;
  int row = 0;

  for ( row = firstRow - 1; row < lastRow; ++row ) {
    const int rowOffset = row * columns1;
    int column = 0;

    for ( column = firstColumn - 1; column < lastColumn; ++column ) {
      const int index1 = rowOffset + column;
      const int index2 = index1 + 1;
      const int index3 = index2 + columns1;
      const int index4 = index1 + columns1;

      CHECK5( IN_RANGE( index1, 0,
                        ( data->rows + 1 ) * ( data->columns + 1 ) - 1 ),
              IN_RANGE( index2, 0,
                        ( data->rows + 1 ) * ( data->columns + 1 ) - 1 ),
              IN_RANGE( index3, 0,
                        ( data->rows + 1 ) * ( data->columns + 1 ) - 1 ),
              IN_RANGE( index4, 0,
                        ( data->rows + 1 ) * ( data->columns + 1 ) - 1 ),
              output < subsetData + data->rows * data->columns );

      {
        const double coordinate1 = input[ index1 ];
        const double coordinate2 = input[ index2 ];
        const double coordinate3 = input[ index3 ];
        const double coordinate4 = input[ index4 ];
        const double centerCoordinate =
          0.25 * ( coordinate1 + coordinate2 + coordinate3 + coordinate4 );
        *output++ = centerCoordinate;
      }
    }
  }

  {
    /* Replicate to other subset layers: */

    const int subsetTimesteps = data->readTimesteps;
    int timestep = 0;
    int layer = 0;
    size_t count = subsetCount;

    for ( layer = 1; layer < subsetLayers; ++layer ) {
      memcpy( output, subsetData, count * sizeof *subsetData );
      output += count;
    }

    /* Replicate to other subset timesteps: */

    count = subsetLayers * subsetCount;

    for ( timestep = 1; timestep < subsetTimesteps; ++timestep ) {
      memcpy( output, subsetData, count * sizeof *subsetData );
      output += count;
    }
  }

  POST0( IMPLIES_ELSE( OR2( ! strcmp( variableName, "LONGITUDE" ),
                            ! strcmp( variableName, "longitude" ) ),
                       AND2( isValidLongitude( subsetData[ 0 ] ),
                             isValidLongitude(subsetData[ subsetCount - 1 ])),
                       AND2( isValidLatitude( subsetData[ 0 ] ),
                             isValidLatitude(subsetData[ subsetCount - 1 ]))));
}



/******************************************************************************
PURPOSE: copySubsetElevations - Compute and copy elevations to all subset
         timesteps/layers/rows/columns.
INPUTS:  Data* const data    Data that specify subset.
OUTPUTS: float subsetData[ data->readTimesteps * subsetLayers * subsetRows *
                           subsetColumns ]  Subset variable data.
******************************************************************************/

static void copySubsetElevations( const Data* const data, float subsetData[]) {

  PRE07( data, data->arguments, isValidArguments( data->arguments ),
         data->columns > 0, data->readTimesteps > 0,
         OR2( data->heights, data->elevations ), subsetData );

  const Arguments* const arguments = data->arguments;
  const int subsetTimesteps = data->readTimesteps;
  const int integrate = arguments->auxMode == INTEGRATE;
  const int haveZF = AND2( arguments->zfFileCount > 0, data->heights );
  const int firstLayer = integrate ? 1 : arguments->subset[ LAYER ][ MINIMUM ];
  const int lastLayer  = integrate ? 1 : arguments->subset[ LAYER ][ MAXIMUM ];
  const int firstRow      = arguments->subset[ ROW    ][ MINIMUM ];
  const int lastRow       = arguments->subset[ ROW    ][ MAXIMUM ];
  const int firstColumn   = arguments->subset[ COLUMN ][ MINIMUM ];
  const int lastColumn    = arguments->subset[ COLUMN ][ MAXIMUM ];
  const float* const heights = data->heights;
  const double* const elevations = data->elevations;
  float* output = subsetData;
  const size_t columns = data->columns;
  const size_t rows    = data->rows;
  const size_t surfaceCells = columns * rows;
  int layer = 0;

  DEBUG( fprintf( stderr, "copySubsetElevatinons(): "
                  "data->layers = %d, data->rows = %d, data->columns = %d, "
                  "subset: timesteps = %d, "
                  "layers = [%d %d], rows = [%d %d], columns = [%d %d], "
                  "haveZF = %d\n",
                  data->layers, data->rows, data->columns,
                  subsetTimesteps,
                  firstLayer, lastLayer, firstRow, lastRow,
                  firstColumn, lastColumn,
                  haveZF ); )

  for ( layer = firstLayer - 1; layer < lastLayer; ++layer ) {
    const size_t layerOffset = layer * surfaceCells;
    int row = 0;

    for ( row = firstRow - 1; row < lastRow; ++row ) {
      const size_t rowOffset = row * columns;
      int column = 0;

      for ( column = firstColumn - 1; column < lastColumn; ++column ) {
        const size_t columnOffset = layerOffset + rowOffset;

        if ( haveZF ) { /* Add HT to ZH: */
          const size_t index = rowOffset + column;
          CHECK( index < data->rows * data->columns );
          const float height = heights[ index ];
          const float zh = *output;
          const float elevation = height + zh;
          *output++ = elevation;
        } else { /* Copy non-time-varying elevations: */
          const size_t index = columnOffset + column;
          CHECK( index < data->layers * data->rows * data->columns );
          const float elevation = elevations[ index ];
          CHECK( IN_RANGE( elevation, ELEVATION_MINIMUM, ELEVATION_MAXIMUM ) );
          *output++ = elevation;
        }
      }
    }
  }

  if ( subsetTimesteps > 1 ) { /* Replicate to other subset timesteps: */
    const size_t subsetLayers  = COUNT_IN_RANGE( firstLayer,  lastLayer );
    const size_t subsetRows    = COUNT_IN_RANGE( firstRow,    lastRow );
    const size_t subsetColumns = COUNT_IN_RANGE( firstColumn, lastColumn );
    const size_t count = subsetLayers * subsetRows * subsetColumns;
    int timestep = 0;

    for ( timestep = 1; timestep < subsetTimesteps; ++timestep ) {
      memcpy( output, subsetData, count * sizeof *subsetData );
      output += count;
    }
  }
}



/******************************************************************************
PURPOSE: expandSubsetData - Copy last row/column edge data to next row/column.
INPUTS:  const int layers        Number of layers.
         const int rows          Number of rows.
         const int columns       Number of columns.
         const int expandRow     Expand data by 1 row?
         const int expandColumn  Expand data by 1 column?
         float subsetData[]  subsetData[ layers ][ rows ][ columns ].
OUTPUTS: float subsetData[]  subsetData[ layers ][ rows + 1 ][ columns + 1 ].
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   subsetData must have been allocated large enough to hold the extra
         row and column of data. UGLY.
******************************************************************************/

static int expandSubsetData( const int layers, const int rows,
                             const int columns, const int expandRow,
                             const int expandColumn, float subsetData[] ) {

  PRE06( layers > 0, rows > 0, columns > 0, IS_BOOL( expandRow ),
         IS_BOOL( expandColumn ), subsetData );

  const size_t count = layers * rows * columns;
  float* copy = NEW( float, count );
  const int result = copy != 0;

  DEBUG( fprintf( stderr, "expandSubsetData( layers = %d, rows = %d, "
                  "columns = %d, expandRow = %d, expandColumn = %d, "
                  "subsetData = [%f ... %f] )\n",
                  layers, rows, columns, expandRow, expandColumn,
                  subsetData[ 0 ],
                  subsetData[ layers * rows * columns - 1 ] ); )

  if ( result ) {
    float* input = copy;
    float* output = subsetData;
    const int columns1 = columns + expandColumn;
    int layer = 0;
    memcpy( copy, subsetData, count * sizeof *copy );

    for ( layer = 0; layer < layers; ++layer ) {
      int row = 0;
      int column = 0;

      for ( row = 0; row < rows; ++row ) {

        for ( column = 0; column < columns; ++column ) {
          *output++ = *input++;
        }

        if ( expandColumn ) {
          *output = *( output - 1 ); /* Copy last column value. */
          ++output;
        }
      }

      if ( expandRow ) { /* Copy last row: */

        for ( column = 0; column < columns1; ++column ) {
          *output = *( output - columns1 );
          ++output;
        }
      }
    }

    DEBUG( fprintf( stderr, "  expanded subsetData = [%f ... %f] )\n",
                    subsetData[ 0 ],
                    subsetData[ layers *
                                ( rows + expandRow ) *
                                ( columns + expandColumn ) - 1 ] ); )

    FREE( copy );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeXDR - Read the subset of data from input files and write it to
         stdout in XDR format.
INPUTS:  Data* const    data    Data structure of subset.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int writeXDR( Data* const data ) {

  PRE03( data, data->arguments, isValidArguments( data->arguments ) );

  int result = writeXDRHeader( data, stdout );

  if ( result ) {
    result = writeXDRData( data, stdout );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createTemporaryXDRFile - Create a temporary XDR format file of
         the subset output.
INPUTS:  Data* const data  Data structure of subset.
OUTPUTS: Initializes file global varaible temporaryFileName.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int createTemporaryXDRFile( Data* const data ) {

  PRE03( data, data->arguments, isValidArguments( data->arguments ) );

  int result = 0;
  FILE* file = 0;
  const int pid = getpid();
  memset( temporaryFileName, 0, sizeof temporaryFileName );
  snprintf( temporaryFileName,
            sizeof temporaryFileName / sizeof *temporaryFileName,
            "%s/%s%d",
            data->arguments->tmpDir, temporaryFilePrefix, pid );
  file = fopen( temporaryFileName, "wb" );

  if ( file ) {
    result = writeXDRHeader( data, file );

    if ( result ) {
      result = writeXDRData( data, file );
    }

    result = AND2( fclose( file ) == 0, result );
    file = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to create temporary XDR file '%s'.\n",
             temporaryFileName );
  }

  POST02( temporaryFileName[ 0 ],
          IMPLIES( result, fileSize( temporaryFileName ) > 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: readXDRHeader - Read XDR format ASCII header.
INPUTS:  FILE* file  XDR file opened for reading.
OUTPUTS  FILE* file  XDR file pointing at binary data array.
         int* timesteps
         int* variables
         int* layers
         int* rows
         int* columns
         int* yyyymmddhh
         char variableNames[ MXVARS3 ][ NAMLEN3 + 1 ]
         char variableUnits[ MXVARS3 ][ NAMLEN3 + 1 ]
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES: Header looks like the following:
SUBSET 9.0 CMAQ
M_02_99BRACE
http://www.epa.gov/ttn/scram/,CMAQSubset
2002-04-23T00:00:00-0000
# data dimensions: timesteps variables layers rows columns:
2 4 2 3 4
# subset indices (0-based time, 1-based layer/row/column): first-timestep last-timestep first-layer last-layer first-row last-row first-column last-column:
 0 1 1 2 1 3 1 4
# Variable names:
LONGITUDE LATITUDE ELEVATION DAILY_MAX8_O3
# Variable units:
deg deg m ppmV
# lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
30 60 40 -100 6.371e+06 6.371e+06
# Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[3]:
4 3 1.578e+06 -1.27e+06 2000 2000 2 10000 1 0.995 0.99
# IEEE-754 32-bit reals data[variables][timesteps][layers][rows][columns]:
******************************************************************************/

static int readXDRHeader( FILE* file,
                          int* timesteps,
                          int* variables,
                          int* layers,
                          int* rows,
                          int* columns,
                          int* yyyymmddhh,
                          char variableNames[ MXVARS3 ][ NAMLEN3 + 1 ],
                          char variableUnits[ MXVARS3 ][ NAMLEN3 + 1 ] ) {

  PRE09( file, timesteps, variables, layers, rows, columns, yyyymmddhh,
         variableNames, variableUnits );

  int yyyy = 0;
  int mm   = 0;
  int dd   = 0;
  int hh = 0;
  int result =
    fscanf( file,
            "%*[^\n]\n"
            "%*[^\n]\n"
            "%*[^\n]\n"
            "%04d-%02d-%02dT%02d:00:00-0000\n"
            "%*[^\n]\n"
            "%d %d %d %d %d\n"
            "%*[^\n]\n"
            "%*[^\n]\n"
            "%*[^\n]\n",
            &yyyy, &mm, &dd, &hh,
            timesteps, variables, layers, rows, columns ) == 9;

  if ( result ) {
    int variable = 0;
    *yyyymmddhh = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
    result =
      AND2( isValidYYYYMMDDHH( *yyyymmddhh ),
            GT_ZERO5( *timesteps, *variables, *layers, *rows, *columns  ) );

    for ( variable = 0; AND2( result, variable < *variables ); ++variable ) {
      result = fscanf( file, "%16s", variableNames[ variable ] ) == 1;
    }

    if ( result ) {
      char line[ 256 ] = "";
      result = /* Skip newline and next line: */
        AND2( fgets( line, sizeof line, file ),
              fgets( line, sizeof line, file ) );

      for ( variable = 0; AND2( result, variable < *variables ); ++variable ) {
        result = fscanf( file, "%16s", variableUnits[ variable ] ) == 1;
      }

      /* Skip rest of header: */

      result =
        AND6( fgets( line, sizeof line, file ),
              fgets( line, sizeof line, file ),
              fgets( line, sizeof line, file ),
              fgets( line, sizeof line, file ),
              fgets( line, sizeof line, file ),
              fgets( line, sizeof line, file ) );
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to read header of temporary file '%s'\n",
                     temporaryFileName );
  }

  POST02( IS_BOOL( result ),
          IMPLIES( result,
                   AND5( GT_ZERO5( *timesteps, *variables, *layers, *rows,
                                   *columns  ),
                         variableNames[ 0 ][ 0 ],
                         variableUnits[ 0 ][ 0 ],
                         variableNames[ *variables - 1 ][ 0 ],
                         variableUnits[ *variables - 1 ][ 0 ] ) ) );

  return result;
}



/******************************************************************************
PURPOSE: writeASCII - Read the subset of data from input files and write it to
         stdout as a tab-delimited ASCII spreadsheet format.
INPUTS:  Data* const    data    Data structure of subset.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   First write temporary XDR file
         then read input_data[V][T][L][R][C] & write output_data[T][L][R][C][V]
******************************************************************************/

static int writeASCII(  Data* const data ) {

  PRE03( data, data->arguments, isValidArguments( data->arguments ) );

  int result = createTemporaryXDRFile( data );

  if ( result ) {
    FILE* temporaryFile = fopen( temporaryFileName, "rb" );
    int timesteps = 0;
    int variables = 0;
    int layers = 0;
    int rows = 0;
    int columns = 0;
    int yyyymmddhh = 0;
    char variableNames[ MXVARS3 ][ NAMLEN3 + 1 ];
    char variableUnits[ MXVARS3 ][ NAMLEN3 + 1 ];
    memset( variableNames, 0, sizeof variableNames );
    memset( variableUnits, 0, sizeof variableUnits );
    result =
      readXDRHeader( temporaryFile,
                     &timesteps,
                     &variables,
                     &layers,
                     &rows,
                     &columns,
                     &yyyymmddhh,
                     variableNames,
                     variableUnits );

    DEBUG( fprintf( stderr, "after readXDRHeader(): "
                    "timesteps = %d, variables = %d, "
                    "layers = %d, rows = %d, columns = %d, "
                    "yyyymmddhh %d, "
                    "variableNames = [%s ... %s], "
                    "variableUnits = [%s ... %s]\n",
                    timesteps, variables, layers, rows, columns, yyyymmddhh,
                     variableNames[ 0 ], variableNames[ variables - 1 ],
                     variableUnits[ 0 ], variableUnits[ variables - 1 ] ); )

    if ( result ) {
      const char* const headerStart = "Timestamp(UTC)";
      const char* const headerFormat = "\t%s(%s)";
      const char* const dataFormat = "\t%28.18e";
      int isDaily = 0;
      int variable = 0;

      /* Write spreadsheet header row: */

      printf( headerStart );

      for ( ;  variable < variables; ++variable ) {
        printf( headerFormat,
                variableNames[ variable ], variableUnits[ variable ] );

        if ( strstr( variableNames[ variable ], "DAILY_" ) ) {
          isDaily = 1;
        }
      }

      putchar( '\n' );

      /* Write data rows: */

      {
        const int cells = layers * rows * columns;
        const int timestepSize = variables * cells;
        const int variableSize = timesteps * cells;
        float* timestepData = NEW( float, timestepSize );
        result = timestepData != 0;

        if ( result ) {
          const size_t dataOffset = ftell( temporaryFile );
          const int hours = isDaily ? 24 : 1;
          int timestep = 0;

          for ( timestep = 0; AND2(result, timestep < timesteps); ++timestep) {

            if ( ! data->isHourlyTimesteps ) {
              yyyymmddhh =
                data->fileTimeRange[ timestep + data->skipFileCount ][ 0 ];
            }

            for ( variable = 0; AND2( result, variable < variables );
                  ++variable ) {
              const int variableOffset = variable * cells;
              const size_t offset =
                variable * variableSize + timestep * cells;
              const size_t offsetBytes = dataOffset + offset * 4;
              result = fseek( temporaryFile, offsetBytes, SEEK_SET ) == 0;

              if ( result ) {
                result =
                  readFloats( cells, timestepData + variableOffset,
                              temporaryFile );
              }
            }

            if ( result ) { /* Write data rows for this timestep: */
              int cell = 0;

              for ( cell = 0; cell < cells; ++cell ) {
                printf( "%04d-%02d-%02dT%02d:00:00-0000",
                        yyyymmddhh / 1000000,
                        yyyymmddhh / 10000 % 100,
                        yyyymmddhh / 100 % 100,
                        yyyymmddhh % 100 );

                for ( variable = 0; variable < variables; ++variable ) {
                  const int index = variable * cells + cell;
                  const float value = timestepData[ index ];
                  printf( dataFormat, value );
                }

                putchar( '\n' );
              }
            }

            if ( data->isHourlyTimesteps ) {
              yyyymmddhh = incrementHours( yyyymmddhh, hours );
            }
          }

          FREE( timestepData );
        }
      }
    }

    fclose( temporaryFile ), temporaryFile = 0;
  }

  if ( *temporaryFileName ) {
    unlink( temporaryFileName );
    *temporaryFileName = '\0';
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDS - Read the subset of data from input files and write it to
         stdout as NetCDF format using COARDS conventions.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   Must first write temporary NetCDF file then stream it to stdout.
******************************************************************************/

static int writeCOARDS( Data* const data ) {

  PRE03( data, data->arguments, isValidArguments( data->arguments ) );

  int file = createTemporaryCOARDSFileHeader( data );
  int result = file != -1;

  if ( result ) {
    result = writeCOARDSData( data, file );
    closeNetCDFFile( file ), file = -1;

    if ( result ) {
      result = streamFile( temporaryFileName );
    }
  }

  if ( *temporaryFileName ) {
    unlink( temporaryFileName );
    *temporaryFileName = '\0';
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createTemporaryCOARDSFileHeader - Create temporary NetCDF-COARDS file
         and write its header.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int file id if successful, else -1 and a failure message is printed to
         stderr
NOTES:   NetCDF-COARDS header looks like:
dimensions:
        time = 2 ;
        z = 2 ;
        y = 3 ;
        x = 4 ;
variables:
        float longitude(z, y, x) ;
                longitude:units = "degrees_east" ;
                longitude:range = "[-180, 180]" ;
        float latitude(z, y, x) ;
                latitude:units = "degrees_north" ;
                latitude:range = "[-90, 90]" ;
        float elevation(z, y, x) ;
                elevation:units = "meters" ;
                elevation:positive = "up" ;
                elevation:datum = "NAD83" ;
        int yyyymmdd(time) ;
                yyyymmdd:units = "yyyymmdd" ;
        int hhmmss(time) ;
                hhmmss:units = "hhmmss" ;
        float DAILY_MAX8_O3(time, z, y, x) ;
                DAILY_MAX8_O3:units = "ppmV" ;
                DAILY_MAX8_O3:missing_value = -9999.f ;
        float time(time) ;
                time:units = "hours since 2002-04-23 00:00:00.0 -00:00" ;

// global attributes:
                :grid = "M_02_99BRACE" ;
                :Conventions = "COARDS" ;
                :history = "http://www.epa.gov/ttn/scram/,CMAQSubset,RSIG3D" ;
                :west_bound = -84.18996f ;
                :east_bound = -84.10011f ;
                :south_bound = 27.0082f ;
                :north_bound = 27.07433f ;
******************************************************************************/

static int createTemporaryCOARDSFileHeader( Data* const data ) {

  PRE03( data, data->arguments, isValidArguments( data->arguments ) );

  const Arguments* const arguments = data->arguments;
  int result = 0;
  int file = 0;
  const int pid = getpid();
  memset( temporaryFileName, 0, sizeof temporaryFileName );
  snprintf( temporaryFileName,
            sizeof temporaryFileName / sizeof *temporaryFileName,
            "%s/%s%d",
            arguments->tmpDir, temporaryFilePrefix, pid );

  file = createNetCDFFile( temporaryFileName );

  if ( file != -1 ) {
    const int integrate = arguments->auxMode == INTEGRATE;
    const int timesteps = data->outputTimesteps;
    const int layers =
      integrate ? 1
      : COUNT_IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                        arguments->subset[ LAYER ][ MAXIMUM ] );
    const int rows =
      COUNT_IN_RANGE( arguments->subset[ ROW ][ MINIMUM ],
                      arguments->subset[ ROW ][ MAXIMUM ] );
    const int columns =
      COUNT_IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                      arguments->subset[ COLUMN ][ MAXIMUM ] );
    int dimids[ 4 ] = { 0, 0, 0, 0 };
    result =
      AND4( createNetCDFDimension( file, "time", timesteps, &dimids[ 0 ] ),
            createNetCDFDimension( file, "z", layers,       &dimids[ 1 ] ),
            createNetCDFDimension( file, "y", rows,         &dimids[ 2 ] ),
            createNetCDFDimension( file, "x", columns,      &dimids[ 3 ] ) );

    if ( result ) {
      const int variables = arguments->variables;
      int variable = 0;
      result = createCOARDSStandardVariables( data, file, dimids );

      for ( ; AND2( result, variable < variables ); ++variable ) {
        const char* const name =
          outputVariableName( arguments->variableNames[ variable ] );
        const char* const units =
          integrate ? "molecules/cm2"
          : outputVariableUnits( arguments->variableNames[ variable ],
                                 data->variableUnits[ variable ], 0 );
        result =
          createNetCDFVariable( file, name, units, 0, 0, 4, dimids ) != -1;
      }

      if ( result ) {

        if ( arguments->auxMode == WIND ) {
          CHECK( variable > 0 );
          result =
            createNetCDFVariable( file, "WWIND",
                                  data->variableUnits[ variable - 1 ],
                                  0, 0, 4, dimids ) != -1;
        }

        if ( result ) {
          int input = openNetCDFFile(arguments->fileNames[ 0 ], 'r');
          result = input >= 0;

          if ( result ) {
            char gdnam[ NAMLEN3 + 1 ] = "";
            memset( gdnam, 0, sizeof gdnam );
            result =
              getNetCDFStringAttribute( input, -1, "GDNAM",
                                        sizeof gdnam / sizeof *gdnam,
                                        gdnam );
            closeNetCDFFile( input ), input = -1;

            if ( result ) {
              const size_t points = ( data->rows + 1  ) * (data->columns + 1 );
              Bounds bounds = { { 0.0, 0.0 }, { 0.0, 0.0 } };
              computeBounds(points, data->longitudes, data->latitudes, bounds);
              result =
                AND7( createNetCDFStringAttribute( file, -1,
                                                   "Conventions", "COARDS" ),
                      createNetCDFStringAttribute( file, -1,
                                                   "history", arguments->note),
                      createNetCDFStringAttribute( file, -1, "grid", gdnam ),
                      createNetCDFDoubleAttribute( file, "west_bound",
                                              bounds[ LONGITUDE ][ MINIMUM ] ),
                      createNetCDFDoubleAttribute( file, "east_bound",
                                              bounds[ LONGITUDE ][ MAXIMUM ] ),
                      createNetCDFDoubleAttribute( file, "south_bound",
                                              bounds[ LATITUDE ][ MINIMUM ] ),
                      createNetCDFDoubleAttribute( file, "north_bound",
                                              bounds[ LATITUDE ][ MAXIMUM ] ) );

              if ( result ) {
                result = endNetCDFHeader( file );
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to create temporary NetCDF file '%s'.\n",
             temporaryFileName );
  } else {
    result = flushNetCDFFile( file );

    if ( ! result ) {
      closeNetCDFFile( file ), file = -1;
    }
  }

  result = file;
  POST0( IMPLIES( result > -1, fileSize( temporaryFileName ) > 0 ) );
  return result;
}




/******************************************************************************
PURPOSE: createCOARDSStandardVariables - Create COARDS standard coordinate and
         time variables.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   NetCDF-COARDS header looks like:
dimensions:
        time = 2 ;
        z = 2 ;
        y = 3 ;
        x = 4 ;
variables:
        float longitude(z, y, x) ;
                longitude:units = "degrees_east" ;
                longitude:range = "[-180, 180]" ;
        float latitude(z, y, x) ;
                latitude:units = "degrees_north" ;
                latitude:range = "[-90, 90]" ;
        float elevation(z, y, x) ;
                elevation:units = "meters" ;
                elevation:positive = "up" ;
                elevation:datum = "NAD83" ;
        int yyyymmdd(time) ;
                yyyymmdd:units = "yyyymmdd" ;
        int hhmmss(time) ;
                hhmmss:units = "hhmmss" ;
        float time(time) ;
                time:units = "hours since 2002-04-23 00:00:00.0 -00:00" ;
        float DAILY_MAX8_O3(time, z, y, x) ;
                DAILY_MAX8_O3:units = "ppmV" ;
                DAILY_MAX8_O3:missing_value = -9999.f ;

// global attributes:
                :grid = "M_02_99BRACE" ;
                :Conventions = "COARDS" ;
                :history = "http://www.epa.gov/ttn/scram/,CMAQSubset,RSIG3D" ;
                :west_bound = -84.18996f ;
                :east_bound = -84.10011f ;
                :south_bound = 27.0082f ;
                :north_bound = 27.07433f ;
******************************************************************************/

static int createCOARDSStandardVariables( Data* const data,
                                          const int file,
                                          const int dimids[ 4 ] ) {

  PRE06( data, data->arguments, isValidArguments( data->arguments ),
         file >= 0, dimids,
         GE_ZERO4( dimids[ 0 ], dimids[ 1 ], dimids[ 2 ], dimids[ 3 ] ) );

  const Arguments* const arguments = data->arguments;
  int result = 1;

  if ( arguments->lonlat ) {
    result =
      GE_ZERO2( createNetCDFVariable( file, "longitude", "degrees_east",
                                      "range", "[-180, 180]", 2, dimids + 2 ),
                createNetCDFVariable( file, "latitude",  "degrees_north",
                                      "range", "[-90, 90]", 2, dimids + 2 ) );
  }

  if ( result ) {

    if ( arguments->elevation ) {
      const int id =
        createNetCDFVariable( file, "elevation", "meters",
                              "positive", "up", 4, dimids );
      result = id != -1;

      if ( result ) {
        result = createNetCDFStringAttribute( file, id, "datum", "NAD83" );
      }
    }

    if ( result ) {
      const int yyyy = data->yyyymmddhh / 1000000;
      const int mm = data->yyyymmddhh / 10000 % 100;
      const int dd = data->yyyymmddhh / 100 % 100;
      const int hh = data->yyyymmddhh % 100;
      char timeUnits[ 64 ] = "";
      const char* const timestepSize =
        data->isHourlyTimesteps ? "hours" : "months";
      memset( timeUnits, 0, sizeof timeUnits );
      snprintf( timeUnits, sizeof timeUnits / sizeof *timeUnits,
                "%s since %4d-%02d-%02d %02d:00:00.0 -00:00",
                timestepSize, yyyy, mm, dd, hh );
      result =
        GE_ZERO3( createNetCDFVariable( file, "time", timeUnits, 0, 0, 1,
                                        dimids ),
                  createNetCDFVariable( file, "yyyymmdd", "yyyymmdd", 0, 0, 1,
                                        dimids ),
                  createNetCDFVariable( file, "hhmmss", "hhmmss", 0, 0, 1,
                                        dimids ) );
    }
  }

  if ( ! result ) {
    fprintf( stderr,
             "\nFailed to create COARDS standard variables in file '%s'.\n",
             temporaryFileName );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeCOARDSData - Write variable data to NetCDF-COARDS file.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int file id if successful, else -1 and a failure message is printed to
         stderr
NOTES:   NetCDF-COARDS variables looks like:
dimensions:
        time = 2 ;
        z = 2 ;
        y = 3 ;
        x = 4 ;
variables:
        float longitude(z, y, x) ;
                longitude:units = "degrees_east" ;
                longitude:range = "[-180, 180]" ;
        float latitude(z, y, x) ;
                latitude:units = "degrees_north" ;
                latitude:range = "[-90, 90]" ;
        float elevation(z, y, x) ;
                elevation:units = "meters" ;
                elevation:positive = "up" ;
                elevation:datum = "NAD83" ;
        int yyyymmdd(time) ;
                yyyymmdd:units = "yyyymmdd" ;
        int hhmmss(time) ;
                hhmmss:units = "hhmmss" ;
        float DAILY_MAX8_O3(time, z, y, x) ;
                DAILY_MAX8_O3:units = "ppmV" ;
                DAILY_MAX8_O3:missing_value = -9999.f ;
        float time(time) ;
                time:units = "hours since 2002-04-23 00:00:00.0 -00:00" ;

// global attributes:
                :grid = "M_02_99BRACE" ;
                :Conventions = "COARDS" ;
                :history = "http://www.epa.gov/ttn/scram/,CMAQSubset,RSIG3D" ;
                :west_bound = -84.18996f ;
                :east_bound = -84.10011f ;
                :south_bound = 27.0082f ;
                :north_bound = 27.07433f ;
******************************************************************************/

static int writeCOARDSData( Data* const data, const int file) {

  PRE04( data, data->arguments, isValidArguments( data->arguments), file > -1);

  const Arguments* const arguments = data->arguments;
  const int integrate = arguments->auxMode == INTEGRATE;
  const int yyyymmddhh1 = arguments->subset[ TIME ][ MINIMUM ];
  const int yyyymmddhh2 = arguments->subset[ TIME ][ MAXIMUM ];
  const int timestepHours =
    arguments->aggregateMode ? 24
    : data->fileTimeRange[ data->skipFileCount ][ 3 ];
  const int subsetLayers =
    COUNT_IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                    arguments->subset[ LAYER ][ MAXIMUM ] );
  const int subsetRows =
    COUNT_IN_RANGE( arguments->subset[ ROW ][ MINIMUM ],
                    arguments->subset[ ROW ][ MAXIMUM ] );
  const int subsetColumns =
    COUNT_IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                    arguments->subset[ COLUMN ][ MAXIMUM ] );
  const size_t subsetCells = subsetLayers * subsetRows * subsetColumns;
  const int writeSubsetCells =
    integrate ? subsetCells / subsetLayers : subsetCells;
  const size_t subsetHours = data->readTimesteps;
  const size_t variableSize = subsetHours * subsetCells;
  const size_t subsetVariables = 1 + integrate * 2; /* var, DENS, ZF. */
  const size_t subsetSize = subsetVariables * variableSize;
  const size_t aggregateAllSize =
    IN3( arguments->aggregateMode, AGGREGATE_MEAN, AGGREGATE_SUM ) ?
         writeSubsetCells : 0;
  float* subsetData = NEW( float, subsetSize + aggregateAllSize * 2 );
  int result = subsetData != 0;

  if ( result ) {
    float* const subsetZF   = integrate ? subsetData + variableSize : 0;
    float* const subsetDENS = integrate ? subsetZF   + variableSize : 0;
    float* const aggregateAllData =
      aggregateAllSize ? subsetData + subsetSize : 0;
    int* const aggregateAllCounts =
      aggregateAllSize ? (int*) aggregateAllData + aggregateAllSize : 0;
    const int variables = arguments->variables + (arguments->auxMode == WIND);
    const int coordinateVariables =
      2 * arguments->lonlat + arguments->elevation;
    int variable = -coordinateVariables;
    const int outputTimesteps = data->outputTimesteps;
    const int writeLayers = integrate ? 1 : subsetLayers;
    int wroteLongitudes = 0;
    int wroteLatitudes = 0;

    do {
      const char* const variableName =
        variable == -3 ? "longitude" :
        variable == -2 ?
          ( coordinateVariables == 2 ? "longitude" : "latitude" ) :
        variable == -1 ?
          ( coordinateVariables == 2 ? "latitude" : "elevation" ) :
        ( variable < arguments->variables ?
            arguments->variableNames[ variable ] : data->wwindVariable );
      const char* const writeVariableName =
        AND2( variable >= arguments->variables,
              ! strcmp( variableName, data->wwindVariable ) ) ? "WWIND"
        : variableName;
      const int isLongitude = ! strcmp( variableName, "longitude" );
      const int isLatitude  = ! strcmp( variableName, "latitude" );
      int yyyymmddhh = yyyymmddhh1;
      int timestep = 0;

      do {

        if ( ! data->isHourlyTimesteps ) {
          yyyymmddhh =
            data->fileTimeRange[ timestep + data->skipFileCount ][ 0 ];
        }

        result =
          readSubset( data, variableName, yyyymmddhh,
                      subsetData, subsetZF, subsetDENS );

        if ( result ) {

          if ( isLongitude ) {

            if ( ! wroteLongitudes ) {
              result =
                writeCOARDS2DVariable( file, variableName,
                                       subsetRows,
                                       subsetColumns,
                                       subsetData );
              wroteLongitudes = 1;
            }
          } else if ( isLatitude ) {

            if ( ! wroteLatitudes ) {
              result =
                writeCOARDS2DVariable( file, variableName,
                                       subsetRows,
                                       subsetColumns,
                                       subsetData );
              wroteLatitudes = 1;
            }
          } else {

            if ( arguments->aggregateMode ) {
              aggregateData( arguments->aggregateMode, subsetHours,
                             writeSubsetCells, subsetData,
                             aggregateAllData, aggregateAllCounts );
            }

            if ( ! aggregateAllData ) {
              const char* const theOutputVariableName =
                outputVariableName( writeVariableName );
              result =
                writeM3IOVariable( file, theOutputVariableName, timestep,
                                   writeLayers, subsetRows, subsetColumns,
                                   subsetData );

              if ( AND2( result, variable == 0 ) ) {
                result = writeCOARDSTimeVariables( file, timestep, yyyymmddhh );
              }

              ++timestep;
            }
          }
        }

        if ( data->isHourlyTimesteps ) {
          yyyymmddhh = incrementHours( yyyymmddhh, timestepHours );
        }

      } while ( AND3( result, timestep < outputTimesteps,
                      yyyymmddhh <= yyyymmddhh2) );

      if ( AND3( aggregateAllData,
                 strcmp( writeVariableName, "longitude" ),
                 strcmp( writeVariableName, "latitude" ) ) ) {
        const char* const theOutputVariableName =
          outputVariableName( writeVariableName );
        result =
          writeM3IOVariable( file, theOutputVariableName, 0,
                             writeLayers, subsetRows, subsetColumns,
                             aggregateAllData );

        if ( AND2( result, variable == 0 ) ) {
          result = writeCOARDSTimeVariables( file, 0, yyyymmddhh1 );
        }
      }

      ++variable;
    } while ( AND2( result, variable < variables ) );

    FREE( subsetData );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeIOAPI - Read the subset of data from input files and write it to
         stdout as NetCDF format using IOAPI conventions.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   Must first write temporary NetCDF file then stream it to stdout.
******************************************************************************/

static int writeIOAPI( Data* const data ) {

  PRE03( data, data->arguments, isValidArguments( data->arguments ) );

  int file = createTemporaryIOAPIFileHeader( data );
  int result = file != -1;

  if ( result ) {
    result = writeIOAPIData( data, file );
    closeNetCDFFile( file ), file = -1;

    if ( result ) {
      result = streamFile( temporaryFileName );
    }
  }

  if ( *temporaryFileName ) {
    unlink( temporaryFileName );
    *temporaryFileName = '\0';
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: createTemporaryIOAPIFileHeader - Create temporary NetCDF-IOAPI file
         and write its header and TFLAG variable data.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int file id if successful, else -1 and a failure message is printed to
         stderr
NOTES:   NetCDF-COARDS header looks like:
dimensions:
        TSTEP = UNLIMITED ; // (24 currently)
        DATE-TIME = 2 ;
        LAY = 1 ;
        VAR = 1 ;
        ROW = 299 ;
        COL = 459 ;
variables:
        int TFLAG(TSTEP, VAR, DATE-TIME) ;
                TFLAG:units = "<YYYYDDD,HHMMSS>" ;
                TFLAG:long_name = "TFLAG           " ;
                TFLAG:var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
        float O3(TSTEP, LAY, ROW, COL) ;
                O3:long_name = "O3              " ;
                O3:units = "ppmV            " ;
                O3:var_desc = "Variable O3                                                                     " ;
// global attributes:
                :IOAPI_VERSION = "$Id: @(#) ioapi library version 3.0 $                                           " ;
                :EXEC_ID = "????????????????                                                                " ;
                :FTYPE = 1 ;
                :CDATE = 2010139 ;
                :CTIME = 143719 ;
                :WDATE = 2010139 ;
                :WTIME = 143719 ;
                :SDATE = 2006184 ;
                :STIME = 0 ;
                :TSTEP = 10000 ;
                :NTHIK = 1 ;
                :NCOLS = 459 ;
                :NROWS = 299 ;
                :NLAYS = 1 ;
                :NVARS = 1 ;
                :GDTYP = 2 ;
                :P_ALP = 33. ;
                :P_BET = 45. ;
                :P_GAM = -97. ;
                :XCENT = -97. ;
                :YCENT = 40. ;
                :XORIG = -2556000. ;
                :YORIG = -1728000. ;
                :XCELL = 12000. ;
                :YCELL = 12000. ;
                :VGTYP = 7 ;
                :VGTOP = 5000.f ;
                :VGLVLS = 1.f, 0.f ;
                :GDNAM = "12US1           " ;
                :UPNAM = "WR_ACONC        " ;
                :VAR-LIST = "O3              " ;
                :FILEDESC = "Concentration file output..." ;
                :HISTORY = "" ;
variables:
        int TFLAG(TSTEP, VAR, DATE-TIME) ;
                TFLAG:units = "<YYYYDDD,HHMMSS>" ;
                TFLAG:long_name = "TFLAG           " ;
                TFLAG:var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
******************************************************************************/

static int createTemporaryIOAPIFileHeader( Data* const data ) {

  PRE03( data, data->arguments, isValidArguments( data->arguments ) );

  const Arguments* const arguments = data->arguments;
  int result = 0;
  int file = 0;
  const int pid = getpid();
  memset( temporaryFileName, 0, sizeof temporaryFileName );
  snprintf( temporaryFileName,
            sizeof temporaryFileName / sizeof *temporaryFileName,
            "%s/%s%d",
            arguments->tmpDir, temporaryFilePrefix, pid );

  file = createNetCDFFile( temporaryFileName );

  if ( file != -1 ) {
    const int integrate = arguments->auxMode == INTEGRATE;
    const int tstep = /* 0 = UNLIMITED. Actual = data->outputTimesteps; */
      IN3( arguments->aggregateMode, AGGREGATE_MEAN, AGGREGATE_SUM ) ? 1 : 0;
    const int firstLayer =
      integrate ? 1 : arguments->subset[ LAYER ][ MINIMUM ];
    const int layers =
      integrate ? 1
      : COUNT_IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                        arguments->subset[ LAYER ][ MAXIMUM ] );
    const int rows =
      COUNT_IN_RANGE( arguments->subset[ ROW ][ MINIMUM ],
                      arguments->subset[ ROW ][ MAXIMUM ] );
    const int columns =
      COUNT_IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                      arguments->subset[ COLUMN ][ MAXIMUM ] );
    const int coordinateVariables =
      arguments->lonlat * 2 + arguments->elevation;
    const int outputVariables =
      arguments->variables +
      coordinateVariables +
      ( arguments->auxMode == WIND );
    int dimids[ 4 ] = { 0, 0, 0, 0 };
    int tflagDimIds[ 3 ] = { 0, 0, 0 };

    result =
      AND6( createNetCDFDimension( file, "TSTEP", tstep, &dimids[0]),
            createNetCDFDimension( file, "DATE-TIME", 2, &tflagDimIds[ 2 ] ),
            createNetCDFDimension( file, "LAY", layers,  &dimids[ 1 ] ),
            createNetCDFDimension( file, "VAR", outputVariables,
                                   &tflagDimIds[ 1 ] ),
            createNetCDFDimension( file, "ROW", rows,    &dimids[ 2 ] ),
            createNetCDFDimension( file, "COL", columns, &dimids[ 3 ] ) );
    tflagDimIds[ 0 ] = dimids[ 0 ]; /* TSTEP. */

    if ( result ) {
      const char* const tflagDesc =
        "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS";
      char units[ NAMLEN3 + 1 ] = "";
      char desc[  MXDLEN3 + 1 ] = "";
      memset( units, 0, sizeof units );
      memset( desc,  0, sizeof desc  );
      result =
        createNetCDFVariable( file, "TFLAG",
                              paddedString( "<YYYYDDD,HHMMSS>", NAMLEN3,units),
                              "var_desc",
                              paddedString( tflagDesc, MXDLEN3, desc ),
                              3, tflagDimIds ) == 0;

      if ( result ) {

        if ( arguments->lonlat ) {
          result =
            GE_ZERO2( createNetCDFVariable( file, "LONGITUDE",
                                            paddedString("deg", NAMLEN3,units),
                                            "var_desc",
                                        paddedString("Longitude [-180, 180].",
                                                     MXDLEN3, desc ),
                                        4, dimids ),
                      createNetCDFVariable( file, "LATITUDE",
                                            paddedString("deg", NAMLEN3,units),
                                            "var_desc",
                                        paddedString( "Latitude [-90, 90].",
                                                      MXDLEN3, desc ),
                                        4, dimids ) );
        }

        if ( result ) {

          if ( arguments->elevation ) {
            result =
              createNetCDFVariable( file, "ELEVATION",
                                    paddedString( "m", NAMLEN3, units ),
                                    "var_desc",
                                   paddedString("Meters above mean sea level.",
                                                 MXDLEN3, desc ),
                                    4, dimids ) != -1;
          }

          if ( result ) {
            char name[ NAMLEN3 + 1 ] = "";
            char varlist[ MXVARS3 * ( NAMLEN3 + 1 ) ];
            const int variables = arguments->variables;
            int variable = 0;
            memset( name, 0, sizeof name );
            memset( varlist, 0, sizeof varlist );

            if ( arguments->lonlat ) {
              strncat( &varlist[0],
                       paddedString( "LONGITUDE", NAMLEN3, name ),
                       NAMLEN3 );
              strncat( &varlist[0],
                       paddedString( "LATITUDE", NAMLEN3, name ),
                       NAMLEN3 );
            }

            if ( arguments->elevation ) {
              strncat( &varlist[0],
                       paddedString( "ELEVATION", NAMLEN3, name ),
                       NAMLEN3 );
            }

            for ( ; AND2( result, variable < variables ); ++variable ) {
              const char* const variableName =
                arguments->variableNames[ variable ];
              const char* const outputName = outputVariableName( variableName );
              paddedString( integrate ? "molecules/cm2"
                              : outputVariableUnits( variableName,
                                              data->variableUnits[ variable ],
                                                     0 ),
                            NAMLEN3, units );
              underscoreToSpace( units );

              paddedString( outputVariableDescription( variableName,
                                          data->variableDescriptions[ variable ] ),
                            MXDLEN3, desc );
              result =
                createNetCDFVariable( file, outputName, units, "var_desc", desc,
                                      4, dimids ) != -1;
              strncat( &varlist[0],
                       paddedString( outputName, NAMLEN3, name ), NAMLEN3 );
            }

            if ( AND2( result, arguments->auxMode == WIND ) ) {
              const char* const variableName = "WWIND";
              const char* const outputName =
                outputVariableName( variableName );
              paddedString( outputVariableUnits(
                              variableName,
                              data->variableUnits[ variable - 1 ], 0 ),
                            NAMLEN3, units );
              underscoreToSpace( units );
              paddedString( outputVariableDescription( variableName,
                                                      "True W component of wind"),
                            MXDLEN3, desc );
              result =
                createNetCDFVariable( file, outputName, units, "var_desc", desc,
                                      4, dimids ) != -1;

              /* Also add WWIND to VAR-LIST: */

              strncat( &varlist[0],
                       paddedString( outputName, NAMLEN3, name ), NAMLEN3 );
            }

            CHECK( *( &varlist[0] + sizeof varlist / sizeof (char) - 1 )
                   == '\0' );

            if ( result ) {
              int input = openNetCDFFile( arguments->fileNames[ 0 ], 'r' );
              result = input >= 0;

              if ( result ) {
                double xorig = 0.0;
                double yorig = 0.0;
                double xcell = 0.0;
                double ycell = 0.0;
                float vgtop = 0.0;
                int vgtyp = 0;
                int nlays = 0;
                float vglvls[ MXLAYS3 + 1 ];
                memset( vglvls, 0, sizeof vglvls );
                result =
                  AND9( getNetCDFDoubleAttribute( input, "XORIG", &xorig ),
                        getNetCDFDoubleAttribute( input, "YORIG", &yorig ),
                        getNetCDFDoubleAttribute( input, "XCELL", &xcell ),
                        getNetCDFDoubleAttribute( input, "YCELL", &ycell ),
                        getNetCDFFloatAttribute(  input, "VGTOP", &vgtop ),
                        getNetCDFIntAttribute(    input, "VGTYP", &vgtyp ),
                        getNetCDFIntAttribute(    input, "NLAYS", &nlays ),
                        IN_RANGE( nlays, 1, MXLAYS3 ),
                        getNetCDFFloatArrayAttribute( input, "VGLVLS",
                                                      nlays + 1, vglvls ) );

                checkAndFixVerticalGridParameters( nlays, &vgtyp, &vgtop,
                                                   vglvls );

                if ( result ) {
                  const int firstColumn =
                    arguments->subset[ COLUMN ][ MINIMUM ];
                  const int firstRow =
                    arguments->subset[ ROW    ][ MINIMUM ];
                  const double xorigSubset = xorig + xcell * (firstColumn - 1);
                  const double yorigSubset = yorig + ycell * (firstRow    - 1);
                  const int yyyyddd = toYYYYDDD( data->yyyymmddhh / 100 );
                  const int hh0000 = data->yyyymmddhh % 100 * 10000;
                  const int outputTSTEP =
                    IN3( arguments->aggregateMode,
                         AGGREGATE_MEAN, AGGREGATE_SUM ) ?
                      hoursInRange( arguments->subset[ TIME ][ MINIMUM ],
                                    arguments->subset[ TIME ][ MAXIMUM ] )
                      * 10000
                    : arguments->aggregateMode ? 24 * 10000
                    : data->fileTimeRange[ data->skipFileCount ][ 3 ] * 10000;
                  int yyyy = 0;
                  int ddd = 0;
                  int hh = 0;
                  int mm = 0;
                  int ss = 0;
                  int hhmmss = 0;
                  int yyyyddd2 = 0;
                  nowUTC( &yyyy, &ddd, &hh, &mm, &ss );
                  yyyyddd2 = yyyy * 1000 + ddd;
                  hhmmss = ( hh * 100 + mm ) * 100 + ss;
                  result =
                    AND30( copyNetCDFAttribute( input, "IOAPI_VERSION", file ),
                           copyNetCDFAttribute( input, "EXEC_ID", file ),
                           copyNetCDFAttribute( input, "FTYPE", file ),
                           createNetCDFIntAttribute( file, "CDATE", yyyyddd2 ),
                           createNetCDFIntAttribute( file, "CTIME", hhmmss ),
                           createNetCDFIntAttribute( file, "WDATE", yyyyddd2 ),
                           createNetCDFIntAttribute( file, "WTIME", hhmmss ),
                           createNetCDFIntAttribute( file, "SDATE", yyyyddd ),
                           createNetCDFIntAttribute( file, "STIME", hh0000),
                           createNetCDFIntAttribute( file, "TSTEP",
                                                   outputTSTEP ),
                           createNetCDFIntAttribute( file, "NTHIK", 1 ),
                           createNetCDFIntAttribute( file, "NCOLS",columns ),
                           createNetCDFIntAttribute( file, "NROWS", rows ),
                           createNetCDFIntAttribute( file, "NLAYS", layers ),
                           createNetCDFIntAttribute( file, "NVARS",
                                                     outputVariables ),
                           copyNetCDFAttribute( input, "GDTYP", file ),
                           copyNetCDFAttribute( input, "P_ALP", file ),
                           copyNetCDFAttribute( input, "P_BET", file ),
                           copyNetCDFAttribute( input, "P_GAM", file ),
                           copyNetCDFAttribute( input, "XCENT", file ),
                           copyNetCDFAttribute( input, "YCENT", file ),
                           createNetCDFDoubleAttribute( file, "XORIG",
                                                        xorigSubset ),
                           createNetCDFDoubleAttribute( file, "YORIG",
                                                        yorigSubset ),
                           copyNetCDFAttribute( input, "XCELL", file ),
                           copyNetCDFAttribute( input, "YCELL", file ),
                           createNetCDFIntAttribute( file, "VGTYP", vgtyp ),
                           createNetCDFFloatAttribute( file, "VGTOP", vgtop ),
                           createNetCDFFloatArrayAttribute( file, "VGLVLS",
                                                            layers + 1,
                                                            vglvls +
                                                              firstLayer - 1 ),
                           copyNetCDFAttribute( input, "GDNAM", file ),
                           copyNetCDFAttribute( input, "UPNAM", file ) );

                  result =
                      AND4( result,
                            createNetCDFStringAttribute( file, -1, "VAR-LIST",
                                                         &varlist[0] ),
                            copyNetCDFAttribute( input, "FILEDESC", file ),
                            copyNetCDFAttribute( input, "HISTORY", file ) );

                  closeNetCDFFile( input ), input = -1;

                  if ( result ) {
                    result = endNetCDFHeader( file );

                    if ( result ) {
                      result = writeTFLAG( data, file );
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to create temporary NetCDF file '%s'.\n",
             temporaryFileName );

    if ( file >= 0 ) {
      closeNetCDFFile( file ), file = -1;
    }
  } else {
    flushNetCDFFile( file );
  }

  result = file;
  POST0( IMPLIES( result > -1, fileSize( temporaryFileName ) > 0 ) );
  return result;
}



/******************************************************************************
PURPOSE: writeTFLAG - Write TFLAG variable data to NetCDF-IOAPI file.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
NOTES:   NetCDF-COARDS TFLAG variable looks like:
int TFLAG(TSTEP, VAR, DATE-TIME) ;
    TFLAG:units = "<YYYYDDD,HHMMSS>" ;
    TFLAG:long_name = "TFLAG           " ;
    TFLAG:var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
******************************************************************************/


static int writeTFLAG( Data* const data, const int file ) {

  PRE04( data, data->arguments, isValidArguments( data->arguments ),
         file > -1 );

  const Arguments* const arguments = data->arguments;
  const int outputTimesteps = data->outputTimesteps;
  const int coordinateVariables =
    2 * arguments->lonlat + arguments->elevation;
  const int outputVariables =
    coordinateVariables +
    arguments->variables + ( arguments->auxMode == WIND );
  int* tflag = NEW( int, outputVariables * 2 );
  int result = tflag != 0;

  if ( result ) {
    const int timestepHours =
      IN3( arguments->aggregateMode, AGGREGATE_MEAN, AGGREGATE_SUM ) ? 1
      : arguments->aggregateMode ? 24
      : data->fileTimeRange[ data->skipFileCount ][ 3 ];

    int yyyymmddhh = data->yyyymmddhh;
    int timestep = 0;

    for ( ; AND2( result, timestep < outputTimesteps ); ++timestep ) {

      if ( ! data->isHourlyTimesteps ) {
        yyyymmddhh =
          data->fileTimeRange[ timestep + data->skipFileCount ][ 0 ];
      }

      {
        const int yyyyddd = toYYYYDDD( yyyymmddhh / 100 );
        const int hh0000 = yyyymmddhh % 100 * 10000;
        int* t = tflag;
        int variables = outputVariables;

        while ( variables-- ) {
          *t++ = yyyyddd;
          *t++ = hh0000;
        }
      }

      result =
        writeTFLAGVariable( file, timestep, outputVariables, 2, tflag );

      if ( data->isHourlyTimesteps ) {
        yyyymmddhh = incrementHours( yyyymmddhh, timestepHours );
      }
    }

    FREE( tflag );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeIOAPIData - Write variable data to NetCDF-IOAPI file.
INPUTS:  Data* const data  Data structure of subset.
RETURNS: int file id if successful, else -1 and a failure message is printed to
         stderr
NOTES:   NetCDF-IOAPI variables looks like:
        float O3(TSTEP, LAY, ROW, COL) ;
                O3:long_name = "O3              " ;
                O3:units = "ppmV            " ;
                O3:var_desc = "Variable O3                                                                     " ;
******************************************************************************/

static int writeIOAPIData( Data* const data, const int file) {

  PRE04( data, data->arguments, isValidArguments( data->arguments), file > -1);

  const Arguments* const arguments = data->arguments;
  const int integrate = arguments->auxMode == INTEGRATE;
  const int yyyymmddhh1 = arguments->subset[ TIME ][ MINIMUM ];
  const int yyyymmddhh2 = arguments->subset[ TIME ][ MAXIMUM ];
  const int timestepHours =
    arguments->aggregateMode ? 24
    : data->fileTimeRange[ data->skipFileCount ][ 3 ];
  const int subsetLayers =
    COUNT_IN_RANGE( arguments->subset[ LAYER ][ MINIMUM ],
                    arguments->subset[ LAYER ][ MAXIMUM ] );
  const int subsetRows =
    COUNT_IN_RANGE( arguments->subset[ ROW ][ MINIMUM ],
                    arguments->subset[ ROW ][ MAXIMUM ] );
  const int subsetColumns =
    COUNT_IN_RANGE( arguments->subset[ COLUMN ][ MINIMUM ],
                    arguments->subset[ COLUMN ][ MAXIMUM ] );
  const size_t subsetCells = subsetLayers * subsetRows * subsetColumns;
  const int writeSubsetCells =
    integrate ? subsetCells / subsetLayers : subsetCells;
  const size_t subsetHours = data->readTimesteps;
  const size_t variableSize = subsetHours * subsetCells;
  const size_t subsetVariables = 1 + integrate * 2; /* var, DENS, ZF. */
  const size_t subsetSize = subsetVariables * variableSize;
  const size_t aggregateAllSize =
    IN3( arguments->aggregateMode, AGGREGATE_MEAN, AGGREGATE_SUM ) ?
      writeSubsetCells
    : 0;
  float* subsetData = NEW( float, subsetSize + aggregateAllSize * 2 );
  int result = subsetData != 0;

  if ( result ) {
    float* const subsetZF   = integrate ? subsetData + variableSize : 0;
    float* const subsetDENS = integrate ? subsetZF   + variableSize : 0;
    float* const aggregateAllData =
      aggregateAllSize ? subsetData + subsetSize : 0;
    int* const aggregateAllCounts =
      aggregateAllSize ? (int*) aggregateAllData + aggregateAllSize : 0;
    const int variables = arguments->variables + (arguments->auxMode == WIND);
    const int coordinateVariables =
      2 * arguments->lonlat + arguments->elevation;
    int variable = -coordinateVariables;
    const int outputTimesteps  = data->outputTimesteps;
    const int outputLayers = integrate ? 1 : subsetLayers;

    do {
      const char* const variableName =
        variable == -3 ? "LONGITUDE" :
        variable == -2 ?
          ( coordinateVariables == 2 ? "LONGITUDE" : "LATITUDE" ) :
        variable == -1 ?
          ( coordinateVariables == 2 ? "LATITUDE" : "ELEVATION" ) :
        ( variable < arguments->variables ?
            arguments->variableNames[ variable ] : data->wwindVariable );
      const char* const writeVariableName =
        AND2( variable >= arguments->variables,
              ! strcmp( variableName, data->wwindVariable ) ) ? "WWIND"
        : variableName;
      int yyyymmddhh = yyyymmddhh1;
      int timestep = 0;

      do {

        if ( ! data->isHourlyTimesteps ) {
          yyyymmddhh =
            data->fileTimeRange[ timestep + data->skipFileCount ][ 0 ];
        }

        result = readSubset( data, variableName, yyyymmddhh,
                             subsetData, subsetZF, subsetDENS );

        if ( result ) {

          if ( AND3( arguments->aggregateMode,
                     strcmp( variableName, "LONGITUDE" ),
                     strcmp( variableName, "LATITUDE" ) ) ) {
            aggregateData( arguments->aggregateMode, subsetHours,
                           writeSubsetCells, subsetData,
                           aggregateAllData, aggregateAllCounts );
          }

          if ( ! aggregateAllData ) {
            const char* const theOutputVariableName =
              outputVariableName( writeVariableName );
            result =
              writeM3IOVariable( file, theOutputVariableName, timestep,
                                 outputLayers, subsetRows, subsetColumns,
                                 subsetData );
             ++timestep;
          }
        }

        if ( data->isHourlyTimesteps ) {
          yyyymmddhh = incrementHours( yyyymmddhh, timestepHours );
        }

      } while ( AND3( result, timestep < outputTimesteps,
                      yyyymmddhh <= yyyymmddhh2 ) );

      if ( aggregateAllData ) {
        const char* const theOutputVariableName =
          outputVariableName( writeVariableName );
        result =
          writeM3IOVariable( file, theOutputVariableName, 0,
                             outputLayers, subsetRows, subsetColumns,
                             aggregateAllData );
      }

      ++variable;
    } while ( AND2( result, variable < variables ) );

    FREE( subsetData );
  }

  DEBUG( fprintf( stderr, "writeIOAPIData(): returning %d\n", result ); )

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: aggregateData - Time-aggregate each cell over all timesteps
         into timestep 0.
INPUTS:  const int mode                   AGGREGATE_DAILY_MEAN, etc.
         const size_t timesteps           Number of timesteps.
         const size_t cells               Number of cells per timestep.
         float data[ timesteps * cells ]  Data to aggregate.
         float mean[ cells ]              0 or running mean.
         int counts[ cells ]              0 or running counts.
OUTPUTS: float data[ cells ]              Aggregated data into timestep 0.
******************************************************************************/

static void aggregateData( const int mode,
                           const size_t timesteps, const size_t cells,
                           float data[], float means[], int counts[] ) {

  PRE05( IN_RANGE( mode, 0, AGGREGATE_MODES - 1 ), timesteps, cells, data,
         IMPLIES_ELSE( IN3( mode, AGGREGATE_MEAN, AGGREGATE_SUM ),
                       NON_ZERO2( means, counts ),
                       IS_ZERO2( means, counts ) ) );

  switch ( mode ) {
  case AGGREGATE_DAILY_MEAN:
    aggregateMean( timesteps, cells, data );
    break;
  case AGGREGATE_DAILY_MAX:
    aggregateMax( timesteps, cells, data );
    break;
  case AGGREGATE_DAILY_MAX8:
    aggregateMax8( timesteps, cells, data );
    break;
  case AGGREGATE_MEAN:
    aggregateAll( timesteps, cells, data, means, counts );
    break;
  default:
    CHECK( mode == AGGREGATE_SUM );
    aggregateSum( timesteps, cells, data , means );
    break;
  }
}



/******************************************************************************
PURPOSE: aggregateMean - Time-aggregate mean of each cell over all timesteps
         into timestep 0.
INPUTS:  const size_t timesteps           Number of timesteps.
         const size_t cells               Number of cells per timestep.
         float data[ timesteps * cells ]  Data to aggregate.
OUTPUTS: float data[ cells ]              Aggregated data into timestep 0.
******************************************************************************/

static void aggregateMean( const size_t timesteps, const size_t cells,
                           float data[] ) {

  PRE03( timesteps, cells, data );

  long long cell = 0; /* OpenMP requires the loop index to be signed. */

#pragma omp parallel for

  for ( cell = 0; cell < cells; ++cell ) {
    size_t timestep = 0;
    size_t count = 0;
    size_t index = cell;
    double sum = 0.0;

    for ( timestep = 0; timestep < timesteps; ++timestep, index += cells ) {
      const double value = data[ index ];

      if ( IS_VALID_VALUE( value ) ) {
        sum += value;
        ++count;
      }
    }

    if ( count ) {
      data[ cell ] = sum / count;
    } else {
      data[ cell ] = BADVAL3;
    }
  }
}



/******************************************************************************
PURPOSE: aggregateMax - Time-aggregate max of data over all timesteps into
         timestep 0.
INPUTS:  const size_t timesteps           Number of timesteps.
         const size_t cells               Number of cells per timestep.
         float data[ timesteps * cells ]  Data to aggregate.
OUTPUTS: float data[ cells ]              Aggregated data into timestep 0.
******************************************************************************/

static void aggregateMax( const size_t timesteps, const size_t cells,
                          float data[] ) {

  PRE03( timesteps > 0, cells > 0, data );

  long long cell = 0; /* OpenMP requires the loop index to be signed. */

#pragma omp parallel for

  for ( cell = 0; cell < cells; ++cell ) {
    size_t timestep = 0;
    size_t index = cell;
    double maximum = BADVAL3;

    for ( timestep = 0; timestep < timesteps; ++timestep, index += cells ) {
      const double value = data[ index ];

      if ( value > maximum ) {
        maximum = value;
      }
    }

    data[ cell ] = maximum;
  }
}



/******************************************************************************
PURPOSE: aggregateMax8 - Time-aggregate cell-wise maximum of each 8-hour
         average of data over all timesteps into timestep 0.
INPUTS:  const size_t timesteps           Number of timesteps.
         const size_t cells               Number of cells per timestep.
         float data[ timesteps * cells ]  Data to aggregate.
OUTPUTS: float data[ cells ]              Aggregated data into timestep 0.
******************************************************************************/

static void aggregateMax8( const size_t timesteps, const size_t cells,
                           float data[] ) {

  PRE03( timesteps, cells, data );

  long long cell = 0; /* OpenMP requires the loop index to be signed. */

#pragma omp parallel for

  for ( cell = 0; cell < cells; ++cell ) {
    size_t timestep = 0;
    size_t index = cell;
    double maximum = BADVAL3;

    for ( timestep = 0; timestep < ( timesteps - 8 );
          ++timestep, index += cells ) {
      const size_t index1 = index  + cells;
      const size_t index2 = index1 + cells;
      const size_t index3 = index2 + cells;
      const size_t index4 = index3 + cells;
      const size_t index5 = index4 + cells;
      const size_t index6 = index5 + cells;
      const size_t index7 = index6 + cells;
      const double value0 = data[ index ];
      const double value1 = data[ index1 ];
      const double value2 = data[ index2 ];
      const double value3 = data[ index3 ];
      const double value4 = data[ index4 ];
      const double value5 = data[ index5 ];
      const double value6 = data[ index6 ];
      const double value7 = data[ index7 ];
#if 0
      /* This method only works if there are no BADVAL3 values (e.g., CMAQ): */

      const double sum =
        value0 + value1 + value2 + value3 + value4 + value5 + value6 + value7;
      const double average = sum * ( 1.0 / 8.0 );
#else
      /* This method skips BADVAL3 values (e.g., OMIBEHRIOAPI): */

      double sum = 0.0;
      double average = BADVAL3;
      size_t count = 0;

      if ( IS_VALID_VALUE( value0 ) ) {
        sum += value0;
        ++count;
      }

      if ( IS_VALID_VALUE( value1 ) ) {
        sum += value1;
        ++count;
      }

      if ( IS_VALID_VALUE( value2 ) ) {
        sum += value2;
        ++count;
      }

      if ( IS_VALID_VALUE( value3 ) ) {
        sum += value3;
        ++count;
      }

      if ( IS_VALID_VALUE( value4 ) ) {
        sum += value4;
        ++count;
      }

      if ( IS_VALID_VALUE( value5 ) ) {
        sum += value5;
        ++count;
      }

      if ( IS_VALID_VALUE( value6 ) ) {
        sum += value6;
        ++count;
      }

      if ( IS_VALID_VALUE( value7 ) ) {
        sum += value7;
        ++count;
      }

      if ( count ) {
        average = sum / count;
      }

#endif

      if ( average > maximum ) {
        maximum = average;
      }
    }

    data[ cell ] = maximum;
  }
}



/******************************************************************************
PURPOSE: aggregateAll - Time-aggregate mean of each cell over all timesteps.
INPUTS:  const size_t timesteps           Number of timesteps.
         const size_t cells               Number of cells per timestep.
         const float data[ timesteps * cells ]  Data to aggregate.
         float means[ cells ]  Running means.
         int counts[ cells ]   Running counts.
OUTPUTS: float means[ cells ]  Updated running means.
         int counts[ cells ]   Updated running counts.
******************************************************************************/

static void aggregateAll( const size_t timesteps, const size_t cells,
                          const float data[], float means[], int counts[] ) {

  PRE05( timesteps, cells, data, means, counts );

  long long cell = 0; /* OpenMP requires the loop index to be signed. */

#pragma omp parallel for

  for ( cell = 0; cell < cells; ++cell ) {
    size_t timestep = 0;
    size_t index = cell;

    for ( timestep = 0; timestep < timesteps; ++timestep, index += cells ) {
      const double value = data[ index ];

      if ( IS_VALID_VALUE( value ) ) {
        const int count = counts[ cell ];

        if ( count == 0 ) {
          means[ cell ] = value;
          counts[ cell ] = 1;
        } else {
          const int count1 = count + 1;
          float* const cellMean = means + cell;
          *cellMean *= count;
          *cellMean += value;
          *cellMean /= count1;
          counts[ cell ] = count1;
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: aggregateSum - Time-aggregate each cell over all timesteps.
INPUTS:  const size_t timesteps           Number of timesteps.
         const size_t cells               Number of cells per timestep.
         const float data[ timesteps * cells ]  Data to aggregate.
         float sums[ cells ]  Running sums.
OUTPUTS: float sums[ cells ]  Updated running sums.
******************************************************************************/

static void aggregateSum( const size_t timesteps, const size_t cells,
                          const float data[], float sums[] ) {

  PRE04( timesteps, cells, data, sums );

  long long cell = 0; /* OpenMP requires the loop index to be signed. */

#pragma omp parallel for

  for ( cell = 0; cell < cells; ++cell ) {
    size_t timestep = 0;
    size_t index = cell;

    for ( timestep = 0; timestep < timesteps; ++timestep, index += cells ) {
      const double value = data[ index ];

      if ( IS_VALID_VALUE( value ) ) {
        sums[ cell ] += value;
      }
    }
  }
}



/******************************************************************************
PURPOSE: computeZ - Compute elevation in meters above mean sea level from
         CMAQ/IOAPI vertical grid parameters.
INPUTS:  const double g         Gravitational force, e.g., 9.81 m/s^2.
         const double R         Gas constant e.g., 287.04 J/kg/K = m^3/s/K.
         const double A         Atmospheric lapse rate, e.g., 50.0 K/kg.
         const double T0s       Reference surface temperature, e.g., 290.0 K.
         const double P00       Reference surface pressure, e.g., 100000 P.
         const int layers       Number of elevation layers.
         const int type         Vertical grid type:
                                VGSGPH3  hydrostatic sigma-P
                                VGSGPN3  non-h sigma-P
                                VGSIGZ3  sigma-Z
                                VGPRES3  pressure (pascals)
                                VGZVAL3  Z (m) (above sea lvl)
                                VGHVAL3  H (m) (above ground)
                                VGWRFEM  WRF sigma-P
                                IMISS3   None (becomes a single layer at 0.0).
         const double topPressure  Pressure in pascals at the top of the model.
         const double heightOfTerrainInMeters  Height of terrain in meters.
         const float levels[]   Vertical grid levels. levels[ layers + 1 ].
OUTPUTS: double z[]             Elevation at each level.
******************************************************************************/

static void computeZ( const double g,
                      const double R,
                      const double A,
                      const double T0s,
                      const double P00,
                      const int layers,
                      const int type,
                      const double topPressure,
                      const double heightOfTerrainInMeters,
                      const float levels[],
                      double z[] ) {

  PRE07( GT_ZERO5( g, R, A, T0s, P00 ),
         layers > 0,
         IS_VALID_VERTICAL_GRID_TYPE( type ),
         topPressure > 0.0,
         heightOfTerrainInMeters >= ELEVATION_MINIMUM,
         levels,
         z );

  const int numberOfLevels = layers + 1;

  if ( IN4( type, VGSGPH3, VGSGPN3, VGWRFEM ) ) {
    /* Compute z using MM5 formula: */
    elevationsAtSigmaPressures( g, R, A, T0s, P00, heightOfTerrainInMeters,
                                numberOfLevels, topPressure, levels, z );
  } else { /* Compute z using other formulas: */
    int level;

    for ( level = 0; level < numberOfLevels; ++level ) {
      const double valueAtLevel = levels[ level ];
      double clampedValueAtLevel = valueAtLevel;
      double pressure = 0.0;

      switch ( type ) {
      case VGSGPH3: /* hydrostatic sigma-P */
      case VGSGPN3: /* non-h sigma-P */
      case VGWRFEM: /* WRF sigma-P */
        clampedValueAtLevel = CLAMPED_TO_RANGE( valueAtLevel, 0.0, 1.0 );
        pressure =
          pressureAtSigmaLevel( clampedValueAtLevel, topPressure / 100.0 );
       z[ level ] = heightAtPressure( pressure );
        break;
      case VGSIGZ3: /* sigma-Z */
        /* vgtop is in meters and valueAtLevel increases for each level. */
        clampedValueAtLevel = CLAMPED_TO_RANGE( valueAtLevel, 0.0, 1.0 );
        z[ level ] = heightOfTerrainInMeters +
                     clampedValueAtLevel *
                     ( topPressure - heightOfTerrainInMeters );
        break;
      case VGPRES3: /* pressure (Pascals) */
        clampedValueAtLevel = CLAMPED_TO_RANGE( valueAtLevel, 1.0, 1e6 );
        z[ level ] = heightAtPressure( clampedValueAtLevel / 100.0 );
        break;
      case VGZVAL3: /* Z (m) (above sea level) */
        clampedValueAtLevel = CLAMPED_TO_RANGE( valueAtLevel, -1e3, 1e5 );
        z[ level ] = clampedValueAtLevel;
        break;
      case VGHVAL3: /* H (m) (above ground)  */
        clampedValueAtLevel = CLAMPED_TO_RANGE( valueAtLevel, 0.0, 1e5 );
        z[ level ] = clampedValueAtLevel + heightOfTerrainInMeters;
        break;
      default:
        z[ level ] = level;
        break;
      }

      DEBUG( fprintf( stderr, "z[ %d ] = %lf\n", level, z[ level ] ); )
    }
  }

  POST02( inRange( layers + 1, z, ELEVATION_MINIMUM, ELEVATION_MAXIMUM ),
          z[ 0 ] <= z[ layers ] );
}



/******************************************************************************
PURPOSE: pressureAtSigmaLevel - Compute pressure (in millibars) at a given
         sigma level.
INPUTS:  const double sigmaLevel     Sigma level.
         const double pressureAtTop  Pressure (in millibars) at top of model.
RETURNS: double pressure in millibars at the given sigma level.
NOTES:   Based on formula in the documentation for Vis5d by Bill Hibbard.
******************************************************************************/

static double pressureAtSigmaLevel( const double sigmaLevel,
                                    const double pressureAtTop ) {
  PRE02( IN_RANGE( sigmaLevel, 0.0, 1.0 ), pressureAtTop > 1.0 );
#define SURFACE_PRESSURE_IN_MB 1012.5
  const double result =
    pressureAtTop + sigmaLevel * ( SURFACE_PRESSURE_IN_MB - pressureAtTop );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: heightAtPressure - Compute the height (in meters) at a given
         pressure (in millibars).
INPUTS:  double pressure   Pressure value in millibars.
RETURNS: double height in meters at the given pressure in millibars.
NOTES:   Based on formula in the documentation for Vis5d by Bill Hibbard.
******************************************************************************/

static double heightAtPressure( double pressure ) {
  PRE0( pressure > 1.0 );
  const double pressureToHeightScaleFactor = -7.2 * 1000.0;
  double result = 0.0;

  if ( pressure <= 0.0 ) {
    pressure = 1e-10; /* HACK: prevent crash on non-IEEE.*/
  }

  result =
    pressureToHeightScaleFactor * log( pressure / SURFACE_PRESSURE_IN_MB );

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: elevationsAtSigmaPressures - Compute elevations in meters above mean
         sea-level at sigma-pressures.
INPUTS:  const double g                Gravitational force, e.g., 9.81 m/s^2.
         const double R                Gas constant e.g., 287.04 J/kg/K = m^3/s/K.
         const double A                Atmospheric lapse rate, e.g., 50.0 K/kg.
         const double T0s              Reference surface temperature, e.g., 290.0 K.
         const double P00              Reference surface pressure, e.g., 100000 P.
         const double surfaceElevation Elevation of surface in meters AMSL.
         const int levels              Number of levels of sigmaPressures.
         const double topPressure      Pressure in Pascals at the top of the model.
         const float sigmaPressures[ levels ]  Sigma-pressures at levels.
OUTPUTS: double elevations[ levels ]  Elevation in meters above MSL at sigmas.
NOTES:   Based on formula used in MM5.
******************************************************************************/

static void elevationsAtSigmaPressures( const double g,
                                        const double R,
                                        const double A,
                                        const double T0s,
                                        const double P00,
                                        const double surfaceElevation,
                                        const int levels,
                                        const double topPressure,
                                        const float sigmaPressures[],
                                        double elevations[] ) {

  PRE07( GT_ZERO5( g , R, A, T0s, P00 ),
         IN_RANGE( surfaceElevation, -1e3, 1e4 ),
         levels > 0,
         IN_RANGE( topPressure, 1e4, 1e5 ),
         inRangeF( levels, sigmaPressures, 0.0, 1.0 ),
         valuesDecrease( levels, sigmaPressures ),
         elevations );

  /* Derived constants: */

  const double H0s            = R * T0s / g;
  const double one_over_H0s   = 1.0 / H0s;
  const double A_over_T0s     = A / T0s;
  const double A_over_two_T0s = A / ( T0s + T0s );
  const double Pt             = topPressure;
  const double Zs             = surfaceElevation;
  const double two_Zs         = Zs + Zs;
  const double sqrt_factor    = sqrt(1.0 - A_over_T0s * one_over_H0s * two_Zs);
  const double q_factor =
    ( Pt / P00 ) * exp( two_Zs * one_over_H0s / sqrt_factor );
  int level = 0;

  /* Compute elevations at sigma-pressures: */

  for ( level = 0; level < levels; ++level ) {
    const double sigma_p0   = sigmaPressures[ level ];
    const double q0_star    = sigma_p0 + ( 1.0 - sigma_p0 ) * q_factor;
    const double ln_q0_star = log( q0_star );
    const double z_level    =
      Zs - H0s * ln_q0_star * ( A_over_two_T0s * ln_q0_star + sqrt_factor );
    elevations[ level ] =
      CLAMPED_TO_RANGE( z_level, ELEVATION_MINIMUM, ELEVATION_MAXIMUM );
  }

  POST02( valuesIncrease( levels, elevations ),
          inRange( levels, elevations, ELEVATION_MINIMUM, ELEVATION_MAXIMUM ) );

}



/******************************************************************************
PURPOSE: integrateLayers - Integrate data (ppmV) over layers to molecules/cm2.
INPUTS:  const size_t timesteps  Number of timesteps.
         const size_t layers     Number of layers.
         const size_t rows       Number of rows.
         const size_t columns    Number of columns.
         const float* const zf   Layer thicknesses in meters:
                                 zf[ timesteps ][ layers ][ rows ][ columns ]
         const float* const dens  Density of dry air in kg/m3:
                                  dens[ timesteps ][ layers ][ rows ][columns]
         float* const data        Concentration data (ppmV) to integrate:
                                  data[ timesteps ][ layers ][ rows ][columns]
OUTPUTS: float* const data        Layer-integrated data in molecules/cm2:
                                  data[ timesteps ][ 1 ][ rows ][ columns ]
NOTES:   Column integrated concentration won't handle BADVAL3 values.
Based on 2017-08-30 Luke Valin EPA/ORD/NERL valin.lukas@epa.gov 919-541-2355

The CMAS M3 package ignores the impact of water vapor variations on the molar
mass (and density) of air, so I will do the same.
This makes the calculation very easy.
Final units should be molecules per cm^2.
Trace gas column =
  Vertical integral of Full Layer height * air mass density / molar mass of air
  * gas mixing ratio * Avogadro’s number

CONC and METCRO Variables:
Layer Height = ZF(i) – ZF (i-1) where  ZF(0) =0;
Air mass density = DENS;
Gas mixing ratio = CONC;

Constants:
REAL, PARAMETER :: AVO   = 6.0221367e23 ! Avogadro's Constant [ number/mol ]
REAL, PARAMETER :: MWAIR = 28.9628 ! mean molecular weight for dry air [g/mol]
                             ! FSB: 78.06% N2, 21% O2, and 0.943% A on a mole
                             ! fraction basis ( Source : Hobbs, 1995) pp. 69-70

REAL, PARAMETER :: DENS_CONV = ( 1.0E3 * AVO / MWAIR ) * 1.0E-6
 ! convert from kg/m**3 to #/cc
REAL, PARAMETER :: PPM_MCM3  = 1.0E-06
 ! convert from ppm to molecules / cc mol_Spec/mol_Air = ppm * 1E-06
REAL, PARAMETER :: M2CM      = 1.0E2   ! meters to centimeters
REAL, PARAMETER :: M2CM1     = 1.0E-6  ! 1/ m**3 to 1/ cm**3

******************************************************************************/

static void integrateLayers( const size_t timesteps,
                             const size_t layers,
                             const size_t rows,
                             const size_t columns,
                             const float* const zf,
                             const float* const dens,
                             float* const data ) {

  PRE07( timesteps > 0, layers > 0, rows > 0, columns > 0, zf, dens, data );

  const double cm_per_m = 100.0; /* Centimenters per meter. */

  /* Cubic meters per cubic centimeter: */

  const double m3_per_cm3 = 1.0 / ( cm_per_m * cm_per_m * cm_per_m );

  const double g_per_kg = 1000.0; /* 1000 g air / kg air. */

  const double Avogadro_molecules_per_mol = 6.022140857e23; /* number/mol. */

  /* grams of air / mole of air: */

  const double mean_molecular_weight_of_dry_air_g_per_mol = 28.9628;

  const double kg_per_m3_to_moles_per_cm3 =
    m3_per_cm3 * g_per_kg / mean_molecular_weight_of_dry_air_g_per_mol;

  const double ppm_to_molecules_per_mole = 1e-6 * Avogadro_molecules_per_mol;

  const size_t layer_cells = rows * columns;
  const size_t cells = layers * layer_cells;
  size_t timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep ) {
    const size_t timestep_offset = timestep * cells;
    const size_t output_offset   = timestep * layer_cells;
    const float* timestep_height        = zf   + timestep_offset;
    const float* timestep_density       = dens + timestep_offset;
    const float* timestep_concentration = data + timestep_offset;
    float* timestep_integration         = data + output_offset;
    size_t layer_cell = 0;

    for ( layer_cell = 0; layer_cell < layer_cells; ++layer_cell,
          ++timestep_height, ++timestep_density, ++timestep_concentration,
          ++timestep_integration ) {
      size_t layer = 0;
      double previous_height_m = 0.0;
      double sum = 0.0;

      for ( layer = 0; layer < layers; ++layer ) {
        const size_t layer_offset = layer * layer_cells;

        const double height_m = timestep_height[ layer_offset ];
        const double cell_layer_thickness_m = height_m - previous_height_m;
        const double cell_layer_thickness_cm =
          cell_layer_thickness_m * cm_per_m;

        const double density_kg_per_m3 = timestep_density[ layer_offset ];
        const double density_moles_per_cm3 =
          density_kg_per_m3 * kg_per_m3_to_moles_per_cm3;

        const double concentration_ppm = timestep_concentration[ layer_offset];
        const double concentration_molecules_per_mole =
          concentration_ppm * ppm_to_molecules_per_mole;

        const double term_molecules_per_cm2 =
          cell_layer_thickness_cm *
          density_moles_per_cm3 *
          concentration_molecules_per_mole;
        sum += term_molecules_per_cm2;
        previous_height_m = height_m;
      }

      *timestep_integration = sum; /*Write layer-integrated result to layer 1*/
    }
  }
}




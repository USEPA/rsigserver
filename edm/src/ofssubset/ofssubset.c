/******************************************************************************
PURPOSE: ofssubset.c - C program to read a sequence of DODS gridded NOAA OFS
         data files, subset it to a given lonlat box and output it in bin
         format.

NOTES:
To compile:
$ cc -o ofssubset ofssubset.c -lm
Usage:
ofssubset source variable units minimum maximum lonmin latmin lonmax latmax \
          yyyymmddhh hours input_files > output.bin
data outside the range [minimum, maximum] is mapped to -9999.
Usage examples:
 See test/get_files and testit script.

https://tidesandcurrents.noaa.gov/models.html
https://opendap.co-ops.nos.noaa.gov/thredds/catalog.html

@ LAYER = 0
@ LAYER1 = $LAYER + 1

# Current (u,v,w) LAYER = 0..19
echo 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/CBOFS/MODELS\
 /201306/nos.cbofs.fields.nowcast.20130617.t00z.nc.dods?\
 s_rho[LAYER:1:LAYER],\
 h[0:1:290][0:1:331],\
 lon_rho[0:1:290][0:1:331],\
 lat_rho[0:1:290][0:1:331],\
 angle[0:1:290][0:1:331],\
 mask_rho[0:1:290][0:1:331],\
 mask_u[0:1:290][0:1:330],\
 mask_v[0:1:289][0:1:331],\
 ocean_time[0:1:6],\
 zeta[0:1:6][0:1:290][0:1:331],\
 u[0:1:6][LAYER:1:LAYER][0:1:290][0:1:330],\
 v[0:1:6][LAYER:1:LAYER][0:1:289][0:1:331],\
 w[0:1:6][LAYER:1:LAYER1][0:1:290][0:1:331]'\
 | sed "s/LAYER1/$LAYER1/g" | sed "s/LAYER/$LAYER/g" \
 | awk -v q="'" -v b='!' \
 '{printf "#%c/bin/csh -f\nwget -q -T 0 -O - %c%s%c\n", b, q, $0, q}' \
 > x ; chmod +x x ; time x > test/current_2013061700.dods ; \rm x

echo 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/CBOFS/MODELS\
 /201306/nos.cbofs.fields.nowcast.20130617.t06z.nc.dods?\
 s_rho[LAYER:1:LAYER],\
 h[0:1:290][0:1:331],\
 lon_rho[0:1:290][0:1:331],\
 lat_rho[0:1:290][0:1:331],\
 angle[0:1:290][0:1:331],\
 mask_rho[0:1:290][0:1:331],\
 mask_u[0:1:290][0:1:330],\
 mask_v[0:1:289][0:1:331],\
 ocean_time[0:1:6],\
 zeta[0:1:6][0:1:290][0:1:331],\
 u[0:1:6][LAYER:1:LAYER][0:1:290][0:1:330],\
 v[0:1:6][LAYER:1:LAYER][0:1:289][0:1:331],\
 w[0:1:6][LAYER:1:LAYER1][0:1:290][0:1:331]'\
 | sed "s/LAYER1/$LAYER1/g" | sed "s/LAYER/$LAYER/g" \
 | awk -v q="'" -v b='!' \
 '{printf "#%c/bin/csh -f\nwget -q -T 0 -O - %c%s%c\n", b, q, $0, q}' \
 > x ; chmod +x x ; time x > test/current_2013061706.dods ; \rm x

echo 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/CBOFS/MODELS\
 /201306/nos.cbofs.fields.nowcast.20130617.t12z.nc.dods?\
 s_rho[LAYER:1:LAYER],\
 h[0:1:290][0:1:331],\
 lon_rho[0:1:290][0:1:331],\
 lat_rho[0:1:290][0:1:331],\
 angle[0:1:290][0:1:331],\
 mask_rho[0:1:290][0:1:331],\
 mask_u[0:1:290][0:1:330],\
 mask_v[0:1:289][0:1:331],\
 ocean_time[0:1:6],\
 zeta[0:1:6][0:1:290][0:1:331],\
 u[0:1:6][LAYER:1:LAYER][0:1:290][0:1:330],\
 v[0:1:6][LAYER:1:LAYER][0:1:289][0:1:331],\
 w[0:1:6][LAYER:1:LAYER1][0:1:290][0:1:331]'\
 | sed "s/LAYER1/$LAYER1/g" | sed "s/LAYER/$LAYER/g" \
 | awk -v q="'" -v b='!' \
 '{printf "#%c/bin/csh -f\nwget -q -T 0 -O - %c%s%c\n", b, q, $0, q}' \
 > x ; chmod +x x ; time x > test/current_2013061712.dods ; \rm x

echo 'http://opendap.co-ops.nos.noaa.gov/thredds/dodsC/NOAA/CBOFS/MODELS\
 /201306/nos.cbofs.fields.nowcast.20130617.t18z.nc.dods?\
 s_rho[LAYER:1:LAYER],\
 h[0:1:290][0:1:331],\
 lon_rho[0:1:290][0:1:331],\
 lat_rho[0:1:290][0:1:331],\
 angle[0:1:290][0:1:331],\
 mask_rho[0:1:290][0:1:331],\
 mask_u[0:1:290][0:1:330],\
 mask_v[0:1:289][0:1:331],\
 ocean_time[0:1:6],\
 zeta[0:1:6][0:1:290][0:1:331],\
 u[0:1:6][LAYER:1:LAYER][0:1:290][0:1:330],\
 v[0:1:6][LAYER:1:LAYER][0:1:289][0:1:331],\
 w[0:1:6][LAYER:1:LAYER1][0:1:290][0:1:331]'\
 | sed "s/LAYER1/$LAYER1/g" | sed "s/LAYER/$LAYER/g" \
 | awk -v q="'" -v b='!' \
 '{printf "#%c/bin/csh -f\nwget -q -T 0 -O - %c%s%c\n", b, q, $0, q}' \
 > x ; chmod +x x ; time x > test/current_2013061718.dods ; \rm x


$ head -17 test/current_2013061700.dods
Dataset {
    Float64 s_rho[s_rho = 1];
    Float64 h[eta_rho = 291][xi_rho = 332];
    Float64 lon_rho[eta_rho = 291][xi_rho = 332];
    Float64 lat_rho[eta_rho = 291][xi_rho = 332];
    Float64 angle[eta_rho = 291][xi_rho = 332];
    Float64 mask_rho[eta_rho = 291][xi_rho = 332];
    Float64 mask_u[eta_u = 291][xi_u = 331];
    Float64 mask_v[eta_v = 290][xi_v = 332];
    Float64 ocean_time[ocean_time = 7];
    Float32 zeta[ocean_time = 7][eta_rho = 291][xi_rho = 332];
    Float32 u[ocean_time = 7][s_rho = 1][eta_u = 291][xi_u = 331];
    Float32 v[ocean_time = 7][s_rho = 1][eta_v = 290][xi_v = 332];
    Float32 w[ocean_time = 7][s_w = 2][eta_rho = 291][xi_rho = 332];
} NOAA%2fCBOFS%2fMODELS%2f201306%2fnos%2ecbofs%2efields%2enowcast%2e20130617%2\
et00z%2enc;

Data:

# Note: DODS input arrays are MSB/big-endian IEEE-754 64-bit floating-point
# format and are preceeded by two 4-byte words (MSB integer) containing the
# count of array items that follows.

$ ls -1 test/current_*.dods > test/current_files

$ ofssubset cbofs current m/s -50 50 -78 35 -70 40 2013061700 24 \
  test/current_files > test/current.bin

$ head -9 test/current.bin
Content-type: application/octet-stream; charset=iso-8859-1
# dimensions: variables timesteps points
4        24       25
# variable names:
depth current_u current_v current_w
# variable units:
m m/s m/s m/s
# MSB 32-bit int yyyymmddhh[timesteps] and
# IEEE-754 32-bit float longitudes[points] and
# IEEE-754 32-bit float latitudes[points] and
# IEEE-754 32-bit float data[variables][timesteps][points]:

HISTORY: 2013-07-22 plessel.todd@epa.gov
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For FILE, stderr, scanf(), printf(), fwrite(). */
#include <stdlib.h> /* For atof(), malloc(), free(). */
#include <string.h> /* For strcmp(), memset(). */
#include <ctype.h>  /* For isspace(). */
#include <math.h>   /* For M_PI. */
#include <time.h>   /* For time_t, struct tm, mktime(), gmtime(). */
#include <assert.h> /* For assert(). */

/*================================= MACROS ==================================*/

#ifndef IS_LITTLE_ENDIAN
#define IS_LITTLE_ENDIAN \
  ( \
    ( defined(__BYTE_ORDER__) && \
      defined(__ORDER_LITTLE_ENDIAN__) && \
      __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ ) || \
    defined(__x86_64__) || \
    defined(__i386__)   || \
    defined(__i486__)   || \
    defined(__i586__)   || \
    defined(__i686__)   || \
    defined(__alpha)    || \
    defined(__ia64__)   || \
    defined(__ARMEL__)  || \
    defined(__MIPSEL)   || \
    defined(_WIN32)     || \
    defined(_WIN64)     || \
    defined(_M_IX86)    || \
    defined(_M_AMD64)      \
  )
#endif

#define IN_RANGE( x, lower, upper ) ( (x) >= (lower) && (x) <= (upper) )

#ifdef DEBUGGING
#define DEBUG(s) s
#else
#define DEBUG(unused)
#endif

#define FREE( p ) ( ( (p) ? free(p) : (void) 0 ), (p) = 0 )
#define NEW( type, count ) allocate( sizeof (type) * (count) )

static void* allocate( size_t bytes ) {
  void* result = 0;
  assert( bytes > 0 );
  result = malloc( bytes );

  if ( result ) {
    memset( result, 0, bytes );
  } else {
    fprintf( stderr, "\nFailed to allocate %lu bytes "
             "to complete the requested action.\n", bytes );
  }

  return result;
}

#define MISSING (-9999.0)

enum {
  HOURS_PER_DAY = 24,
  MINUTES_PER_HOUR = 60,
  SECONDS_PER_MINUTE = 60,
  SECONDS_PER_HOUR = SECONDS_PER_MINUTE * MINUTES_PER_HOUR,
  SECONDS_PER_HALF_HOUR = SECONDS_PER_HOUR / 2,
  SECONDS_PER_DAY = SECONDS_PER_HOUR * HOURS_PER_DAY
};

/*================================== TYPES ==================================*/


/* Command-line arguments: */

typedef struct {
  const char* source;   /* Input data type: cbofs, etc. */
  const char* variable; /* Main output data variable. E.g., "current". */
  const char* units;
  double minimum;       /* Min/max of valid data variable. */
  double maximum;
  double longitude_minimum;
  double latitude_minimum;
  double longitude_maximum;
  double latitude_maximum;
  int yyyymmddhh;       /* Starting timestamp. */
  int hours;            /* Number of hours. */
  const char* files;
} Arguments;


/* Output data: */

enum { MAXIMUM_VARIABLES = 32, MAXIMUM_LINE_LENGTH = 255 };

typedef char Name[ 31 + 1 ];
typedef char Line[ MAXIMUM_LINE_LENGTH + 1 ];

typedef struct {
  int timesteps;
  int variables;
  int points;
  Name names[ MAXIMUM_VARIABLES ]; /* names[ variables ]. */
  Name units[ MAXIMUM_VARIABLES ]; /* units[ variables ]. */
  int* yyyymmddhh;   /* yyyymmddhh[ timesteps ]. */
  float* longitudes; /* longitudes[ points ]. */
  float* latitudes;  /* latitudes[ points ]. */
  float* data;       /* data[ variables ][ timesteps ][ points ]. */
} OutputData;

static void deallocate_output_data( OutputData* data ) {

  if ( data ) {
    FREE( data->yyyymmddhh );
    FREE( data->longitudes );
    FREE( data->latitudes );
    FREE( data->data );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}


/*
 Input data from CBOFS DODS files whose header looks like:

 Dataset {
    Float64 s_rho[s_rho = 1];
    Float64 s_w[s_w = 21];
    Float64 h[eta_rho = 291][xi_rho = 332];
    Float64 lon_rho[eta_rho = 291][xi_rho = 332];
    Float64 lat_rho[eta_rho = 291][xi_rho = 332];
    Float64 lon_u[eta_u = 291][xi_u = 331];
    Float64 lat_u[eta_u = 291][xi_u = 331];
    Float64 lon_v[eta_v = 290][xi_v = 332];
    Float64 lat_v[eta_v = 290][xi_v = 332];
    Float64 angle[eta_rho = 291][xi_rho = 332];
    Float64 mask_rho[eta_rho = 291][xi_rho = 332];
    Float64 mask_u[eta_u = 291][xi_u = 331];
    Float64 mask_v[eta_v = 290][xi_v = 332];
    Float64 ocean_time[ocean_time = 1];
    Float32 zeta[ocean_time = 1][eta_rho = 291][xi_rho = 332];
    Float32 Pair[ocean_time = 1][eta_rho = 291][xi_rho = 332];
    Float32 Uwind[ocean_time = 1][eta_rho = 291][xi_rho = 332];
    Float32 Vwind[ocean_time = 1][eta_rho = 291][xi_rho = 332];
    Float32 u[ocean_time = 1][s_rho = 1][eta_u = 291][xi_u = 331];
    Float32 v[ocean_time = 1][s_rho = 1][eta_v = 290][xi_v = 332];
    Float32 w[ocean_time = 1][s_w = 2][eta_rho = 291][xi_rho = 332];
    Float32 temp[ocean_time = 1][s_rho = 1][eta_rho = 291][xi_rho = 332];
    Float32 salt[ocean_time = 1][s_rho = 1][eta_rho = 291][xi_rho = 332];
    Float32 oxygen[ocean_time = 1][s_rho = 1][eta_rho = 291][xi_rho = 332];
} NOAA%2fCBOFS%2fMODELS%2f201309%2fnos%2ecbofs%2efields%2en006%2e20130923...

Data:
*/

typedef struct {
  int timesteps;     /* Of time-varying data. */
  int rows;          /* rho points. */
  int columns;       /* rho points. */
  double  s;         /* Normalized depth coordinates [-1, 0]. */
  double* h;         /* h[rows][columns] bathymmetry depth in meters. */
  double* longitude; /* lon_rho[rows][columns] longitude in deg [-180, 180]. */
  double* latitude;  /* lat_rho[rows][columns] latitude in deg [-90, 90]. */
  double* angle;     /* angle[rows][columns] angle in radians to adjust u, v.*/
  double* mask;      /* mask_rho[rows][columns] 0 = missing, 1 = valid. */
  double* mask_u;    /* mask_u[rows][columns-1] 0 = missing, 1 = valid. */
  double* mask_v;    /* mask_v[rows-1][columns] 0 = missing, 1 = valid. */
  double* seconds;   /* ocean_time[timesteps] seconds since 2009-01-01 00Z. */
  float*  msl;       /* zeta[timesteps][rows][columns] m above/below datum. */
  float*  air_pressure;      /* Pair[timesteps][rows][columns] hPa. */
  float*  wind_u;            /* Uwind[timesteps][rows][columns] m/s. */
  float*  wind_v;            /* Vwind[timesteps][rows][columns] m/s. */
  float*  current_u;         /* u[timesteps][rows][columns-1] m/s. */
  float*  current_v;         /* v[timesteps][rows-1][columns] m/s. */
  float*  current_w;         /* w[timesteps][rows][columns] m/s. */
  float*  water_temperature; /* temp[timesteps][rows][columns] C. */
  float*  salinity;          /* salt[timesteps][rows][columns] PSU. */
  float*  oxygen;            /* oxygen[timesteps][rows][columns] PSU. */
} CBOFSData;

static void deallocate_cbofs_data( CBOFSData* data ) {

  if ( data ) {
    FREE( data->h );
    FREE( data->longitude );
    FREE( data->latitude );
    FREE( data->angle );
    FREE( data->mask );
    FREE( data->mask_u );
    FREE( data->mask_v );
    FREE( data->seconds );
    FREE( data->msl );
    FREE( data->air_pressure );
    FREE( data->wind_u );
    FREE( data->wind_v );
    FREE( data->current_u );
    FREE( data->current_v );
    FREE( data->current_w );
    FREE( data->water_temperature );
    FREE( data->salinity );
    FREE( data->oxygen );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}


/*
 Input data from GBOFS DODS files whose header looks like:

Dataset {
    Float32 time[time = 25];
    Float32 lon[ny = 134][nx = 73];
    Float32 lat[ny = 134][nx = 73];
    Float32 mask[ny = 134][nx = 73];
    Float32 depth[ny = 134][nx = 73];
    Float32 sigma[sigma = 1];
    Float32 zeta[time = 25][ny = 134][nx = 73];
    Float32 air_u[time = 25][ny = 134][nx = 73];
    Float32 air_v[time = 25][ny = 134][nx = 73];
    Float32 u[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 v[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 w[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 temp[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 salt[time = 25][sigma = 1][ny = 134][nx = 73];
} NOAA%2fGBOFS%2fMODELS%2f201309%2fnos%2egbofs%2efields%2enowcast%2e20130901...

Data:
*/

typedef struct {
  int timesteps;    /* Of time-varying data. */
  int file_base_yyyymmddhh; /* DODS data 'file' timestamp. */
  int rows;         /* ny points. */
  int columns;      /* nx points. */
  float s;          /* Normalized depth coordinate [-1, 0]. */
  float* h;         /* h[rows][columns] bathymmetry depth in meters. */
  float* longitude; /* lon[rows][columns] longitude in deg [-180, 180]. */
  float* latitude;  /* lat[rows][columns] latitude in deg [-90, 90]. */
  float* mask;      /* mask_rho[rows][columns] 0 = missing, 1 = valid. */
  float* seconds;   /* time[timesteps] seconds since 2008-01-01 00Z. */
  float* msl;       /* zeta[timesteps][rows][columns] m above/below datum. */
  float* wind_u;    /* Uwind[timesteps][rows][columns] m/s. */
  float* wind_v;    /* Vwind[timesteps][rows][columns] m/s. */
  float* current_u; /* u[timesteps][1][rows][columns-1] m/s. */
  float* current_v; /* v[timesteps][1][rows-1][columns] m/s. */
  float* current_w; /* w[timesteps][2][rows][columns] m/s. */
  float* water_temperature; /* temp[timesteps][1][rows][columns] C. */
  float* salinity;          /* salt[timesteps][1][rows][columns] PSU. */
} GBOFSData;

static void deallocate_gbofs_data( GBOFSData* data ) {

  if ( data ) {
    FREE( data->h );
    FREE( data->longitude );
    FREE( data->latitude );
    FREE( data->mask );
    FREE( data->seconds );
    FREE( data->msl );
    FREE( data->wind_u );
    FREE( data->wind_v );
    FREE( data->current_u );
    FREE( data->current_v );
    FREE( data->current_w );
    FREE( data->water_temperature );
    FREE( data->salinity );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}


/*
 Input data from NGOFS DODS files whose header looks like a subset of:

Dataset {
    Float32 lon[node = 90267];
    Float32 lat[node = 90267];
    Float32 lonc[nele = 174474];
    Float32 latc[nele = 174474];
    Float32 siglay[siglay = 1][node = 90267];
    Float32 h[node = 90267];
    Float32 nv[three = 3][nele = 174474];
    Float32 time[time = 7];
    Float32 zeta[time = 7][node = 90267];
    Float32 u[time = 7][siglay = 1][nele = 174474];
    Float32 v[time = 7][siglay = 1][nele = 174474];
    Float32 tauc[time = 7][nele = 174474];
    Float32 temp[time = 7][siglay = 1][node = 90267];
    Float32 salinity[time = 7][siglay = 1][node = 90267];
    Float32 short_wave[time = 7][node = 90267];
    Float32 net_heat_flux[time = 7][node = 90267];
    Float32 uwind_speed[time = 7][nele = 174474];
    Float32 vwind_speed[time = 7][nele = 174474];
    Int32 wet_nodes[time = 7][node = 90267];
    Int32 wet_cells[time = 7][nele = 174474];
} NOAA%2fNGOFS%2fMODELS%2f201210%2fnos%2engofs%2efields%2enowcast%2e\
20121001%2et03z%2enc;

Data:
*/

static const int NGOFS_nodes = 90267;

typedef struct {
  int timesteps;            /* Of time-varying data. 1 or 7. */
  int points;               /* Number of data points. */
  int nodes;                /* Number of nodes per point or 0 if not current.*/
  float* mask;              /* mask[points] 0 = missing, 1 = valid. */
  float* s;                 /* Normalized depth coordinates[points] [-1, 0]. */
  float* h;                 /* h[points] bathymmetry depth in meters. */
  float* node;              /* node[3][points] nodes per element. */
  float* longitude;         /* lon[points] longitude in deg [-180, 180]. */
  float* latitude;          /* lat[points] latitude in deg [-90, 90]. */
  float* seconds;           /* time[timesteps] seconds since 2012-01-01 00Z. */
  float* msl;               /* zeta[timesteps][points] m above/below datum. */
  float* wind_u;            /* uwind_speed[timesteps][points] m/s. */
  float* wind_v;            /* vwind_speed[timesteps][points] m/s. */
  float* current_u;         /* u[timesteps][1][points] m/s. */
  float* current_v;         /* v[timesteps][1][points] m/s. */
  float* water_temperature; /* temp[timesteps][1][points] C. */
  float* salinity;          /* salt[timesteps][1][points] PSU. */
  float* seabed_stress;     /* tauc[timesteps][points] PSU. */
  float* radiation;         /* short_wave[timesteps][points] PSU. */
  float* net_heat_flux;     /* net_heat_flux[timesteps][points] PSU. */
} NGOFSData;

static void deallocate_ngofs_data( NGOFSData* data ) {

  if ( data ) {
    FREE( data->mask );
    FREE( data->s );
    FREE( data->h );
    FREE( data->node );
    FREE( data->longitude );
    FREE( data->latitude );
    FREE( data->seconds );
    FREE( data->msl );
    FREE( data->wind_u );
    FREE( data->wind_v );
    FREE( data->current_u );
    FREE( data->current_v );
    FREE( data->water_temperature );
    FREE( data->salinity );
    FREE( data->seabed_stress );
    FREE( data->radiation );
    FREE( data->net_heat_flux );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}


/*
 Input data from NGOFS2 DODS files whose header looks like a subset of:

Dataset {
    Float32 lon[node = 303714];
    Float32 lat[node = 303714];
    Float32 lonc[nele = 569405];
    Float32 latc[nele = 569405];
    Float32 siglay[siglay = 1][node = 303714];
    Float32 h[node = 303714];
    Float32 nv[three = 3][nele = 569405];
    Float64 time[time = 7];
    Float32 zeta[time = 7][node = 303714];
    Float32 u[time = 7][siglay = 1][nele = 569405];
    Float32 v[time = 7][siglay = 1][nele = 569405];
    Float32 tauc[time = 7][nele = 569405];
    Float32 temp[time = 7][siglay = 1][node = 303714];
    Float32 salinity[time = 7][siglay = 1][node = 303714];
    Float32 short_wave[time = 7][node = 303714];
    Float32 net_heat_flux[time = 7][node = 303714];
    Float32 sensible_heat_flux[time = 7][node = 303714];
    Float32 latent_heat_flux[time = 7][node = 303714];
    Float32 atmos_press[time = 7][node = 303714];
    Float32 uwind_speed[time = 7][nele = 569405];
    Float32 vwind_speed[time = 7][nele = 569405];
    Int32 wet_nodes[time = 7][node = 303714];
    Int32 wet_cells[time = 7][nele = 569405];
} NOAA%2fNGOFS2%2fMODELS%2f201210%2fnos%2engofs2%2efields%2enowcast%2e\
20121001%2et03z%2enc;

Data:
*/

static const int NGOFS2_nodes = 303714;

typedef struct {
  int timesteps;            /* Of time-varying data. 1 or 7. */
  int points;               /* Number of data points. */
  int nodes;                /* Number of nodes per point or 0 if not current.*/
  float* mask;              /* mask[points] 0 = missing, 1 = valid. */
  float* s;                 /* Normalized depth coordinates[points] [-1, 0]. */
  float* h;                 /* h[points] bathymmetry depth in meters. */
  float* node;              /* node[3][points] nodes per element. */
  float* longitude;         /* lon[points] longitude in deg [-180, 180]. */
  float* latitude;          /* lat[points] latitude in deg [-90, 90]. */
  double* seconds;          /* time[timesteps] seconds since 2012-01-01 00Z. */
  float* msl;               /* zeta[timesteps][points] m above/below datum. */
  float* wind_u;            /* uwind_speed[timesteps][points] m/s. */
  float* wind_v;            /* vwind_speed[timesteps][points] m/s. */
  float* current_u;         /* u[timesteps][1][points] m/s. */
  float* current_v;         /* v[timesteps][1][points] m/s. */
  float* water_temperature; /* temp[timesteps][1][points] C. */
  float* salinity;          /* salt[timesteps][1][points] PSU. */
  float* seabed_stress;     /* tauc[timesteps][points] m2/s2. */
  float* radiation;         /* short_wave[timesteps][points] W/m2. */
  float* net_heat_flux;     /* net_heat_flux[timesteps][points] W/m2. */
  float* sensible_heat_flux; /* sensible_heat_flux[timesteps][points] W/m2. */
  float* latent_heat_flux;  /* latent_heat_flux[timesteps][points] W/m2. */
  float* air_pressure;      /* atmos_press[timesteps][points] hPa. */
} NGOFS2Data;

static void deallocate_ngofs2_data( NGOFS2Data* data ) {

  if ( data ) {
    FREE( data->mask );
    FREE( data->s );
    FREE( data->h );
    FREE( data->node );
    FREE( data->longitude );
    FREE( data->latitude );
    FREE( data->seconds );
    FREE( data->msl );
    FREE( data->wind_u );
    FREE( data->wind_v );
    FREE( data->current_u );
    FREE( data->current_v );
    FREE( data->water_temperature );
    FREE( data->salinity );
    FREE( data->seabed_stress );
    FREE( data->radiation );
    FREE( data->net_heat_flux );
    FREE( data->sensible_heat_flux );
    FREE( data->latent_heat_flux );
    FREE( data->air_pressure );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}

/*
 Input data from SSCOFS DODS files whose header looks like a subset of:

 Dataset {
    Float32 lon[node = 239734];
    Float32 lat[node = 239734];
    Float32 lonc[nele = 433410];
    Float32 latc[nele = 433410];
    Float32 siglay[siglay = 10][node = 239734];
    Float32 h[node = 239734];
    Int32 nv[three = 3][nele = 433410];
    Float64 time[time = 1];
    Float32 zeta[time = 1][node = 239734];
    Float32 u[time = 1][siglay = 10][nele = 433410];
    Float32 v[time = 1][siglay = 10][nele = 433410];
    Float32 tauc[time = 1][nele = 433410];
    Float32 temp[time = 1][siglay = 10][node = 239734];
    Float32 salinity[time = 1][siglay = 10][node = 239734];
    Float32 short_wave[time = 1][node = 239734];
    Float32 net_heat_flux[time = 1][node = 239734];
    Float32 sensible_heat_flux[time = 1][node = 239734];
    Float32 latent_heat_flux[time = 1][node = 239734];
    Float32 long_wave[time = 1][node = 239734];
    Float32 uwind_speed[time = 1][nele = 433410];
    Float32 vwind_speed[time = 1][nele = 433410];
    Float32 atmos_press[time = 1][node = 239734];
    Int32 wet_nodes[time = 1][node = 239734];
    Int32 wet_cells[time = 1][nele = 433410];
} NOAA/SSCOFS/MODELS/2024/09/16/sscofs.t21z.20240916.fields.n006.nc;

Data:
*/

static const int SSCOFS_nodes = 239734;

typedef struct {
  int timesteps;            /* Of time-varying data. 1 or 7. */
  int points;               /* Number of data points. */
  int nodes;                /* Number of nodes per point or 0 if not current.*/
  float* mask;              /* mask[points] 0 = missing, 1 = valid. */
  float* s;                 /* Normalized depth coordinates[points] [-1, 0]. */
  float* h;                 /* h[points] bathymmetry depth in meters. */
  float* node;              /* node[3][points] nodes per element. */
  float* longitude;         /* lon[points] longitude in deg [-180, 180]. */
  float* latitude;          /* lat[points] latitude in deg [-90, 90]. */
  double* seconds;          /* time[timesteps] seconds since 2018-01-01 00Z. */
  float* msl;               /* zeta[timesteps][points] m above/below datum. */
  float* wind_u;            /* uwind_speed[timesteps][points] m/s. */
  float* wind_v;            /* vwind_speed[timesteps][points] m/s. */
  float* current_u;         /* u[timesteps][1][points] m/s. */
  float* current_v;         /* v[timesteps][1][points] m/s. */
  float* water_temperature; /* temp[timesteps][1][points] C. */
  float* salinity;          /* salt[timesteps][1][points] PSU. */
  float* seabed_stress;     /* tauc[timesteps][points] m2/s2. */
  float* radiation;         /* short_wave[timesteps][points] W/m2. */
  float* net_heat_flux;     /* net_heat_flux[timesteps][points] W/m2. */
  float* sensible_heat_flux; /* sensible_heat_flux[timesteps][points] W/m2. */
  float* latent_heat_flux;  /* latent_heat_flux[timesteps][points] W/m2. */
  float* air_pressure;      /* atmos_press[timesteps][points] hPa. */
} SSCOFSData;

static void deallocate_sscofs_data( SSCOFSData* data ) {

  if ( data ) {
    FREE( data->mask );
    FREE( data->s );
    FREE( data->h );
    FREE( data->node );
    FREE( data->longitude );
    FREE( data->latitude );
    FREE( data->seconds );
    FREE( data->msl );
    FREE( data->wind_u );
    FREE( data->wind_v );
    FREE( data->current_u );
    FREE( data->current_v );
    FREE( data->water_temperature );
    FREE( data->salinity );
    FREE( data->seabed_stress );
    FREE( data->radiation );
    FREE( data->net_heat_flux );
    FREE( data->sensible_heat_flux );
    FREE( data->latent_heat_flux );
    FREE( data->air_pressure );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}


/*
 Input data from SFBOFS DODS files whose header looks like a subset of:

Dataset {
    Float32 lon[node = 54120];
    Float32 lat[node = 54120];
    Float32 lonc[nele = 102264];
    Float32 latc[nele = 102264];
    Float32 siglay[siglay = 20][node = 54120];
    Float32 h[node = 54120];
    Float64 time[time = 1];
    Float32 zeta[time = 1][node = 54120];
    Float32 u[time = 1][siglay = 20][nele = 102264];
    Float32 v[time = 1][siglay = 20][nele = 102264];
    Float32 tauc[time = 1][nele = 102264];
    Float32 temp[time = 1][siglay = 20][node = 54120];
    Float32 salinity[time = 1][siglay = 20][node = 54120];
    Float32 short_wave[time = 1][node = 54120];
    Float32 net_heat_flux[time = 1][node = 54120];
    Float32 sensible_heat_flux[time = 1][node = 54120];
    Float32 latent_heat_flux[time = 1][node = 54120];
    Float32 uwind_speed[time = 1][nele = 102264];
    Float32 vwind_speed[time = 1][nele = 102264];
    Float32 atmos_press[time = 1][node = 54120];
    Int32 wet_nodes[time = 1][node = 54120];
    Int32 wet_cells[time = 1][nele = 102264];
} NOAA/SFBOFS/MODELS/201507/nos.sfbofs.fields.n006.20150727.t15z.nc;

Data:
*/

static const int SFBOFS_nodes = 54120;

typedef struct {
  int timesteps;            /* Of time-varying data. 1 or 7. */
  int points;               /* Number of data points. */
  int nodes;                /* Number of nodes per point or 0 if not current.*/
  float* mask;              /* mask[points] 0 = missing, 1 = valid. */
  float* s;                 /* Normalized depth coordinates[points] [-1, 0]. */
  float* h;                 /* h[points] bathymmetry depth in meters. */
  float* node;              /* node[3][points] nodes per element. */
  float* longitude;         /* lon[points] longitude in deg [-180, 180]. */
  float* latitude;          /* lat[points] latitude in deg [-90, 90]. */
  double* seconds;          /* time[timesteps] seconds since 2012-01-01 00Z. */
  float* msl;               /* zeta[timesteps][points] m above/below datum. */
  float* air_pressure;      /* atmos_press[timesteps][points] hPa. */
  float* wind_u;            /* uwind_speed[timesteps][points] m/s. */
  float* wind_v;            /* vwind_speed[timesteps][points] m/s. */
  float* current_u;         /* u[timesteps][1][points] m/s. */
  float* current_v;         /* v[timesteps][1][points] m/s. */
  float* water_temperature; /* temp[timesteps][1][points] C. */
  float* salinity;          /* salinity[timesteps][1][points] PSU. */
  float* seabed_stress;     /* tauc[timesteps][points] m2/s2. */
  float* radiation;         /* short_wave[timesteps][points] W/m2. */
  float* net_heat_flux;     /* net_heat_flux[timesteps][points] W/m2. */
  float* sensible_heat_flux; /* senisble_heat_flux[timesteps][points] W/m2. */
  float* latent_heat_flux;  /* latent_heat_flux[timesteps][points] W/m2. */
} SFBOFSData;

static void deallocate_sfbofs_data( SFBOFSData* data ) {

  if ( data ) {
    FREE( data->mask );
    FREE( data->s );
    FREE( data->h );
    FREE( data->node );
    FREE( data->longitude );
    FREE( data->latitude );
    FREE( data->seconds );
    FREE( data->msl );
    FREE( data->air_pressure );
    FREE( data->wind_u );
    FREE( data->wind_v );
    FREE( data->current_u );
    FREE( data->current_v );
    FREE( data->water_temperature );
    FREE( data->salinity );
    FREE( data->seabed_stress );
    FREE( data->radiation );
    FREE( data->net_heat_flux );
    FREE( data->sensible_heat_flux );
    FREE( data->latent_heat_flux );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}


/*
 Input data from CREOFS DODS files whose header looks like a subset of:

Dataset {
    Float32 lon[node = 74061];
    Float32 lat[node = 74061];
    Float32 h[node = 74061];
    Float32 sigma[sigma = 1];
    Float32 zeta[time = 1][node = 74061];
    Float32 Pair[time = 1][node = 74061];
    Float32 uwind_speed[time = 1][node = 74061];
    Float32 vwind_speed[time = 1][node = 74061];
    Float32 temp[time = 1][nv = 1][node = 74061];
    Float32 salinity[time = 1][nv = 1][node = 74061];
    Float32 u[time = 1][nv = 1][node = 74061];
    Float32 v[time = 1][nv = 1][node = 74061];
} NOAA%2fCREOFS%2fMODELS%2f201210%2fnos%2ecreofs%2efields%2en00%2e\
20121001%2et21z%2enc;

Data:
*/

typedef struct {
  int file_base_yyyymmddhh; /* DODS data 'file' timestamp. */
  int points;               /* Number of data points. */
  float  s;                 /* Normalized depth coordinates[points] [-1, 0]. */
  float* mask;              /* mask[points] 0 = missing, 1 = valid. */
  float* h;                 /* h[points] bathymmetry depth in meters. */
  float* longitude;         /* lon[points] longitude in deg [-180, 180]. */
  float* latitude;          /* lat[points] latitude in deg [-90, 90]. */
  float* msl;               /* zeta[1][points] m above/below datum. */
  float* air_pressure;      /* Pair[1][points] m/s. */
  float* wind_u;            /* uwind_speed[1][points] m/s. */
  float* wind_v;            /* vwind_speed[1][points] m/s. */
  float* current_u;         /* u[1][1][points] m/s. */
  float* current_v;         /* v[1][1][points] m/s. */
  float* water_temperature; /* temp[1][1][points] C. */
  float* salinity;          /* salinity[1][1][points] PSU. */
} CREOFSData;

static void deallocate_creofs_data( CREOFSData* data ) {

  if ( data ) {
    FREE( data->mask );
    FREE( data->h );
    FREE( data->longitude );
    FREE( data->latitude );
    FREE( data->msl );
    FREE( data->air_pressure );
    FREE( data->wind_u );
    FREE( data->wind_v );
    FREE( data->current_u );
    FREE( data->current_v );
    FREE( data->water_temperature );
    FREE( data->salinity );
    memset( data, 0, sizeof *data );
    FREE( data );
  }
}


/* Dispatcher for readers: */

typedef OutputData* (*Reader)( const Arguments*, FILE* );

typedef struct {
  const char* const source;  /* Input data type. E.g., "cbofs". */
  const int base_yyyymmddhh; /* Base timestamp for time array. */
  const int timestep_scale;  /* time array * scale is seconds. */
  Reader reader;             /* Reads and converts input data to OutputData. */
} DispatchTableEntry;

/* Declare readers for new _OFSData here and add an entry to the table: */

static OutputData* read_CBOFS( const Arguments* arguments, FILE* file );
static OutputData* read_GBOFS( const Arguments* arguments, FILE* file );
static OutputData* read_NGOFS( const Arguments* arguments, FILE* file );
static OutputData* read_NGOFS2( const Arguments* arguments, FILE* file );
static OutputData* read_SFBOFS( const Arguments* arguments, FILE* file );
static OutputData* read_CREOFS( const Arguments* arguments, FILE* file );
static OutputData* read_SSCOFS( const Arguments* arguments, FILE* file );

static const DispatchTableEntry dispatch_table[] = {
  { "cbofs",  2016010100, 1, read_CBOFS }, /* Curvilinear. */
  { "dbofs",  2016010100, 1, read_CBOFS }, /* Same format as CBOFS. */
  { "tbofs",  2009010100, 1, read_CBOFS }, /* Same format as CBOFS. */
  { "gomofs", 2016010100, 1, read_CBOFS }, /* Same format as CBOFS. */
  { "ciofs",  2016010100, 1, read_CBOFS }, /* Same format as CBOFS. */
  { "wcofs",  2016010100, 1, read_CBOFS }, /* Same format as CBOFS. */
  { "gbofs",  2003010100, SECONDS_PER_DAY, read_GBOFS  }, /* Easy curvilinear*/
  { "nyofs",  2008010100, SECONDS_PER_DAY, read_GBOFS  }, /* Same as GBOFS. */
  { "sjrofs", 2016010100, SECONDS_PER_DAY, read_GBOFS  }, /* Same as GBOFS.*/
  { "ngofs",  2012010100, SECONDS_PER_DAY, read_NGOFS  }, /* Unstructured. */
  { "ngofs2", 2019010100, 1,               read_NGOFS2 }, /* Unstructured. */
  { "sfbofs", 2013010100, 1,               read_SFBOFS }, /* Unstructured. */
  { "nwgofs", 2013010100, 1,               read_SFBOFS }, /* Same as SFBOFS. */
  { "creofs", 0,          1,               read_CREOFS }, /* Unstructured. */
  { "sscofs", 2018010100, 1,               read_SSCOFS }, /* Unstructured. */
  { 0, 0, 0, 0 } /* Terminator. */
};

static int lookup( const char* source ) {
  const int count = sizeof dispatch_table / sizeof *dispatch_table;
  int result = -1;
  int index = 0;

  for ( index = 0; index < count; ++index ) {
    const DispatchTableEntry* const entry = dispatch_table + index;

    if ( entry->source && source && ! strcmp( entry->source, source ) ) {
      result = index;
      index = count;
    }
  }

  return result;
}

static Reader lookup_reader( const char* source ) {
  const int index = lookup( source );
  const Reader result = index == -1 ? 0 : dispatch_table[ index ].reader;
  return result;
}

static int base_yyyymmddhh0( const char* source ) {
  const int index = lookup( source );
  const int result = index == -1 ? 0 : dispatch_table[ index ].base_yyyymmddhh;
  return result;
}

static int base_yyyymmddhh( const Arguments* arguments ) {
  const int yyyymmddhh = base_yyyymmddhh0( arguments->source );
  const int result = yyyymmddhh ? yyyymmddhh : arguments->yyyymmddhh;
  return result;
}

static int timestep_scale( const char* source ) {
  const int index = lookup( source );
  const int result = index == -1 ? 0 : dispatch_table[ index ].timestep_scale;
  return result;
}


/* Linked list of output data: */

typedef struct Node { OutputData* data; struct Node* next; } Node;

static Node* append_list_data( Node* node, OutputData* data ) {
  Node* new_node = NEW( Node, 1 );
  assert( data );

  if ( new_node ) {
    new_node->data = data;
    new_node->next = 0;

    if ( node ) {

      while ( node->next ) {
        node = node->next;
      }

      node->next = new_node;
    }
  }

  return new_node;
}

static void deallocate_list( Node* node ) {

  while ( node ) {
    Node* const next = node->next;
    deallocate_output_data( node->data );
    node->data = 0;
    node->next = 0;
    FREE( node );
    node = next;
  }
}

/*=========================== FORWARD DECLARATIONS ==========================*/

static void usage( const char* program );

static int parse_options( const int argc, const char* const argv[],
                          Arguments* arguments );

static void filter_duplicate_timestamps( Node* list );

static int write_output( const Node* list );

static void count_output_points( const Node* list,
                                 int* total_timesteps, int* maximum_points );

static void write_header( const Node* list,
                          const int timesteps, const int points );

static int write_output_timestamps( const Node* list,
                                    const int total_timesteps );

static int write_output_coordinates( const Node* node,
                                     const int maximum_points,
                                     float* const buffer );

static int write_output_data( const Node* list,
                              const int maximum_points,
                              float* const buffer );

/* CBOFS: */

static CBOFSData* read_CBOFS_header( const Arguments* arguments, FILE* file );

static int read_CBOFS_data( const Arguments* arguments, FILE* file,
                            CBOFSData* data );

static OutputData* convert_CBOFS( const Arguments* arguments, CBOFSData* data);


/* GBOFS: */

static GBOFSData* read_GBOFS_header( const Arguments* arguments, FILE* file );

static int read_GBOFS_data( const Arguments* arguments, FILE* file,
                            GBOFSData* data );

static OutputData* convert_GBOFS( const Arguments* arguments, GBOFSData* data);


/* NGOFS: */

static NGOFSData* read_NGOFS_header( const Arguments* arguments, FILE* file );

static int read_NGOFS_data( const Arguments* arguments, FILE* file,
                            NGOFSData* data );

static OutputData* convert_NGOFS( const Arguments* arguments, NGOFSData* data);


/* NGOFS2: */

static NGOFS2Data* read_NGOFS2_header( const Arguments* arguments, FILE* file );

static int read_NGOFS2_data( const Arguments* arguments, FILE* file,
                             NGOFS2Data* data );

static OutputData* convert_NGOFS2(const Arguments* arguments, NGOFS2Data* data);


/* SFBOFS: */

static SFBOFSData* read_SFBOFS_header( const Arguments* arguments, FILE* file);

static int read_SFBOFS_data( const Arguments* arguments, FILE* file,
                             SFBOFSData* data );

static OutputData* convert_SFBOFS(const Arguments* arguments,SFBOFSData* data);


/* CREOFS: */

static CREOFSData* read_CREOFS_header( const Arguments* arguments, FILE* file);

static int read_CREOFS_data( const Arguments* arguments, FILE* file,
                             CREOFSData* data );

static OutputData* convert_CREOFS(const Arguments* arguments,CREOFSData* data);

/* SSCOFS: */

static SSCOFSData* read_SSCOFS_header( const Arguments* arguments, FILE* file );

static int read_SSCOFS_data( const Arguments* arguments, FILE* file,
                             SSCOFSData* data );

static OutputData* convert_SSCOFS(const Arguments* arguments, SSCOFSData* data);



/* Helpers: */

static int spatially_subset( const double longitude_minimum,
                             const double longitude_maximum,
                             const double latitude_minimum,
                             const double latitude_maximum,
                             const int count,
                             const int word_size,
                             const void* longitudes,
                             const void* latitudes,
                             void* mask );

static void copy_unmasked_values( const int count, const double mask[],
                                  const double values[], float output[] );

static void copy_unmasked_values2( const int count, const double mask[],
                                   const float values[], float output[] );

static void copy_unmasked_float_values( const int count, const float mask[],
                                        const float values[], float output[] );

static void copy_depth( const int count, const double mask[],
                        const double s, const double h[],
                        const float msl[], float output[] );

static void copy_float_depth( const int count, const float mask[],
                              const double s, const float s2[],
                              const float h[],
                              const float msl[], float output[] );

static void copy_float_depth3( const int elements,
                               const float node[],
                               const float mask[],
                               const float s[],
                               const float h[],
                               const float msl[],
                               float output[] );

static void copy_wind( const int count, const double mask[],
                       const double angle[],
                       const float u[], const float v[],
                       float u_east[], float v_north[] );

static void copy_current( const int timestep, const CBOFSData* data,
                          float u_east[], float v_north[], float w_up[] );

static int read_ignored_line( FILE* file );

static int read_matched_line( FILE* file, const char* const match );

static int read_dimension1_line( FILE* file, const char* const format,
                                 int* dimension );

static int read_dimension2_line( FILE* file, const char* const format,
                                 int* dimension1, int* dimension2 );

static int read_file_base_yyyymmddhh( FILE* file );

static double* read_new_double_data( FILE* file,
                                     const double scale,
                                     const double offset,
                                     const double minimum,
                                     const double maximum,
                                     const size_t count );

static float* read_new_float_data( FILE* file,
                                   const double scale,
                                   const double offset,
                                   const double minimum,
                                   const double maximum,
                                   const size_t count );

static float* read_new_int_data( FILE* file,
                                 const double scale,
                                 const double offset,
                                 const double minimum,
                                 const double maximum,
                                 const size_t count );

static int read_double_data( FILE* file,
                             const double scale, const double offset,
                             const double minimum, const double maximum,
                             const size_t count, double data[] );

static int read_float_data( FILE* file,
                            const double scale, const double offset,
                            const double minimum, const double maximum,
                            const size_t count, float data[] );

static int read_int_data( FILE* file,
                          const double scale, const double offset,
                          const double minimum, const double maximum,
                          const size_t count, float data[] );

static int skip_8_bytes( FILE* file );

static void reverse_8_byte_words_if_little_endian( const size_t count,
                                                   void* array );

static void reverse_4_byte_words_if_little_endian( const size_t count,
                                                   void* array );

static void transform_point2( const double angle,
                              const double x, const double y,
                              float* xp, float* yp );

static double seconds_difference(const int yyyymmddhh1, const int yyyymmddhh0);

static int convert_seconds( const int yyyymmddhh0, const int seconds );

static int is_valid_yyyymmddhh( const int yyyymmddhh );

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours );

static void decrement_yyyymmddhh( int* yyyymmddhh, const int hours );

static int days_in_month( const int year, const int month );



/* Read input options and DODS files and write bin data to stdout: */

int main( int argc, char* argv[] ) {
  int ok = 0;
  Arguments arguments;
  memset( &arguments, 0, sizeof arguments );

  if ( parse_options( argc, (const char**) argv, &arguments ) ) {
    Reader reader = lookup_reader( arguments.source );

    if ( reader ) {
      FILE* inputs = fopen( arguments.files, "r" );

      if ( inputs ) {
        Node* list = 0;
        Node* tail = 0;
        char file_name[ 256 ] = "";
        memset( file_name, 0, sizeof file_name );

        while ( fgets( file_name,
                       sizeof file_name / sizeof *file_name, inputs ) ) {
          char* const newline = strrchr( file_name, '\n' );

          if ( newline ) {
            *newline = '\0';
          }

          {
            FILE* input = fopen( file_name, "rb" );
            DEBUG( fprintf( stderr, "\n\nfile_name = '%s', input = %p\n",
                            file_name, input ); )

            if ( input ) {
              OutputData* output_data = reader( &arguments, input );

              if ( output_data ) {
                tail = append_list_data( tail, output_data );

                if ( list == 0 ) {
                  list = tail;
                }
              }

              fclose( input ), input = 0;
            }
          }
        }

        fclose( inputs ), inputs = 0;

        if ( list ) {
          filter_duplicate_timestamps( list );
          ok = write_output( list );
          deallocate_list( list ), list = tail = 0;
        }
      }
    }
  }

  return ! ok;
}



/* Print program usage instructions to stderr: */

static void usage( const char* program ) {
  fprintf( stderr, "\n%s - Read a sequence of DODS gridded NOAA OFS data "
                   "files,\nsubset it to a given lonlat box and output it "
                   "in bin format.\n", program );
  fprintf( stderr, "usage: %s source variable units minimum maximum "
                   "lonmin latmin lonmax latmax input_files "
                   "> output.bin\n", program );
  fprintf( stderr, "example: %s cbofs current m/s -50 50 -78 35 -70 40 "
                   "2013061700 24 test/current_files > current.bin\n",
                   program );
  fprintf( stderr, "head -8 current.bin\n\n" );
}



/* Read/check command-line options: */

static int parse_options( const int argc, const char* const argv[],
                          Arguments* arguments ) {

  int result = 0;
  memset( arguments, 0, sizeof *arguments );

  if ( argc == 13 ) {
    arguments->source = argv[ 1 ];

    if ( lookup_reader( arguments->source ) ) {
      arguments->variable = argv[ 2 ];

      if ( arguments->variable[ 0 ] ) {
        arguments->units = argv[ 3 ];

        if ( arguments->units[ 0 ] ) {
          arguments->minimum = atof( argv[ 4 ] );
          arguments->maximum = atof( argv[ 5 ] );

          if ( arguments->maximum > arguments->minimum ) {
            arguments->longitude_minimum = atof( argv[ 6 ] );

            if ( IN_RANGE( arguments->longitude_minimum, -180.0, 180.0 ) ) {
              arguments->latitude_minimum = atof( argv[ 7 ] );

              if ( IN_RANGE( arguments->latitude_minimum, -90.0, 90.0 ) ) {
                arguments->longitude_maximum = atof( argv[ 8 ] );

                if ( IN_RANGE( arguments->longitude_maximum,
                               arguments-> longitude_minimum, 180.0 ) ) {
                  arguments->latitude_maximum = atof( argv[ 9 ] );

                  if ( IN_RANGE( arguments->latitude_maximum,
                                 arguments->latitude_minimum, 90.0 ) ) {
                    arguments->yyyymmddhh = atoi( argv[ 10 ] );

                    if ( is_valid_yyyymmddhh( arguments->yyyymmddhh ) ) {
                      arguments->hours = atoi( argv[ 11 ] );

                      if ( arguments->hours > 0 ) {
                        arguments->files = argv[ 12 ];
                        result = arguments->files[ 0 ] != '\0';
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

  if ( ! result ) {
    fprintf( stderr, "\nInvalid command-line options.\n" );
    usage( argv[ 0 ] );
  }

  return result;
}



/*
 * Filter duplicate timestamps.
 * DODS data from sequential files have overlapping first/last timestamps
 * so filter-out these redundant timestamps.
 */

static void filter_duplicate_timestamps( Node* list ) {
  Node* node = 0;

  for ( node = list; node && node->next; node = node->next ) {

    if ( node->data ) {
      const int timesteps = node->data->timesteps;
      const int last_yyyymmddhh = node->data->yyyymmddhh[ timesteps - 1 ];

      while ( node->next &&
              node->next->data &&
              node->next->data->yyyymmddhh[ 0 ] <= last_yyyymmddhh ) {
        Node* next_node = node->next;
        OutputData* const next_node_data = next_node->data;

        /*
         * the first timestamp of the next node <=
         * the last timestamp of the current node
         * so remove the first timestamp of the next node
         * and if the next node has no more timesteps then remove it.
         */

        const int next_timesteps_1 = next_node_data->timesteps - 1;
        int* const yyyymmddhh = next_node_data->yyyymmddhh;
        int timestep = 0;
        DEBUG( fprintf( stderr,
                        "last_yyyymmddhh = %d, yyyymmddhh[ 0 ] = %d\n",
                        last_yyyymmddhh, yyyymmddhh[ 0 ] ); )

        /* Copy down the next node's timestamps: */

        for ( timestep = 0; timestep < next_timesteps_1; ++timestep ) {
          yyyymmddhh[ timestep ] = yyyymmddhh[ timestep + 1 ];
          DEBUG( fprintf( stderr, "yyyymmddhh[ %d ] = %d\n",
                          timestep, yyyymmddhh[ timestep ] ); )
        }

        yyyymmddhh[ next_timesteps_1 ] = (int) MISSING;

        if ( next_timesteps_1 == 0 ) { /* Remove the next node. */
          DEBUG( fprintf( stderr, "Removing next node.\n" ); )
          node->next = next_node->next;
          deallocate_output_data( next_node_data );
          next_node->data = 0;
          next_node->next = 0;
          FREE( next_node );
        } else if ( next_node_data->variables == 1 ) {
          DEBUG( fprintf( stderr, "Decrementing next node timesteps.\n" );)
          next_node_data->timesteps = next_timesteps_1;
        } else { /* Erase last timestep data[variable][timestep][point]: */
          const int variables = next_node_data->variables;
          const int points = next_node_data->points;
          const int timesteps_points = next_node_data->timesteps * points;
          const int variable_size = next_timesteps_1 * points;
          const float* read_data = next_node_data->data + timesteps_points;
          float* write_data = next_node_data->data + variable_size;
          const size_t copy_bytes = variable_size * sizeof *write_data;
          int variable = 0;

          for ( variable = 1; variable < variables; ++variable ) {
            DEBUG( fprintf( stderr, "Shifting data for variable %d.\n",
                           variable ); )
            assert( read_data > next_node_data->data );
            assert( read_data <=
                      next_node_data->data +
                      ( next_node_data->variables - 1 ) *
                      next_node_data->timesteps * next_node_data->points );
            assert( write_data > next_node_data->data );
            assert( write_data < read_data );
            memcpy( write_data, read_data, copy_bytes );
            write_data += variable_size;
            read_data += timesteps_points;
          }

          next_node_data->timesteps = next_timesteps_1;
        }
      }
    }
  }
}



/* Write bin-format ASCII header and IEEE-754 binary data to stdout: */

static int write_output( const Node* list ) {
  int result = 0;
  int total_timesteps = 0;
  int maximum_points = 0;
  count_output_points( list, &total_timesteps, &maximum_points );

  if ( total_timesteps > 0 && maximum_points > 0 ) {
    write_header( list, total_timesteps, maximum_points );

    if ( write_output_timestamps( list, total_timesteps ) ) {
      float* buffer = NEW( float, maximum_points );

      if ( buffer ) {

        if ( write_output_coordinates( list,  maximum_points, buffer ) ) {
          result = write_output_data( list, maximum_points, buffer );
        }

        FREE( buffer );
      }
    }
  }

  return result;
}



/* Count total number of timesteps and maximum number of points per timestep:*/

static void count_output_points( const Node* list,
                                 int* total_timesteps, int* maximum_points ) {
  const Node* node = list;
  int variables = 0;
  *total_timesteps = *maximum_points = 0;

  while ( node ) {

    if ( node->data ) {

      if ( variables == 0 ) {
        variables = node->data->variables;
      }

      assert( node->data->variables == variables );
      *total_timesteps += node->data->timesteps;
      *maximum_points =
        node->data->points > *maximum_points ? node->data->points
        : *maximum_points;
    }

    node = node->next;
  }
}



/* write ASCII header: */

static void write_header( const Node* list,
                          const int timesteps, const int points ) {
  const OutputData* data = list->data;
  const int variables = data->variables;
  int variable = 0;

  /* Write ASCII header: */

  printf( "Content-type: application/octet-stream; charset=iso-8859-1\n" );
  printf( "# dimensions: variables timesteps points\n" );
  printf( "%-3d %-5d %-10d \n", data->variables, timesteps, points );
  printf( "# variable names:\n" );

  for ( variable = 0; variable < variables; ++variable ) {
    const char delimiter = variable < variables - 1 ? ' ' : '\n';
    assert( variable < sizeof data->names / sizeof data->names[ 0 ] );
    printf( "%s%c", data->names[ variable ], delimiter );
  }

  printf( "# variable units:\n" );

  for ( variable = 0; variable < variables; ++variable ) {
    const char delimiter = variable < variables - 1 ? ' ' : '\n';
    assert( variable < sizeof data->units / sizeof data->units[ 0 ] );
    printf( "%s%c", data->units[ variable ], delimiter );
  }

  printf( "# MSB 32-bit int yyyymmddhh[timesteps] and\n" );
  printf( "# IEEE-754 32-bit float longitudes[points] and\n" );
  printf( "# IEEE-754 32-bit float latitudes[points] and\n" );
  printf( "# IEEE-754 32-bit float data[variables][timesteps][points]:\n" );
}



/* Write big-endian copy of yyyymmddhh timestamps to stdout: */

static int write_output_timestamps( const Node* list,
                                    const int total_timesteps ) {
  int result = 0;
  int* yyyymmddhh = NEW( int, total_timesteps );
  assert( list ); assert( total_timesteps > 0 );

  if ( yyyymmddhh ) {
    const Node* node = list;
    int index = 0;

    for ( node = list; node; node = node->next ) {
      const OutputData* const output_data = node->data;
      const int timesteps = output_data->timesteps;
      const int* const input = output_data->yyyymmddhh;
      int timestep = 0;
      assert( timesteps > 0 );
      assert( is_valid_yyyymmddhh( input[ 0 ] ) );

      for ( timestep = 0; timestep < timesteps; ++timestep, ++index ) {
        assert( IN_RANGE( index, 0, total_timesteps ) );
        assert( IN_RANGE( timestep, 0, output_data->timesteps ) );
        yyyymmddhh[ index ] = input[ timestep ];
      }
    }

    assert( is_valid_yyyymmddhh( yyyymmddhh[ 0 ] ) );
    assert( is_valid_yyyymmddhh( yyyymmddhh[ total_timesteps - 1 ] ) );
    reverse_4_byte_words_if_little_endian( total_timesteps, yyyymmddhh );
    result =
      fwrite( yyyymmddhh, sizeof (int), total_timesteps, stdout )
      == total_timesteps;
    FREE( yyyymmddhh );
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write timestamps.\n" );
  }

  return result;
}



/* Write big-endian copy of longitude-latitude coordinates to stdout: */

static int write_output_coordinates( const Node* node,
                                     const int maximum_points,
                                     float* const buffer ) {
  int result = 0;
  int variable = 0;

  for ( variable = 0, result = 1; result && variable < 2; ++variable ) {
    const OutputData* const output_data = node->data;
    const int points = output_data->points;
    const float* input =
      variable == 0 ? output_data->longitudes : output_data->latitudes;
    float* copy = buffer;
    int point = 0;

    DEBUG( fprintf( stderr,
                    "node points = %d, "
                    "maximum_points = %d, input[ 0 ] = %f\n",
                    points, maximum_points, input[ 0 ] ); )
    assert( points <= maximum_points );

    for ( point = 0; point < points; ++point ) {
      assert( copy < buffer + maximum_points );
      *copy++ = *input++;
    }

    for ( ; point < maximum_points; ++point ) {
      assert( copy < buffer + maximum_points );
      *copy++ = MISSING;
    }

    reverse_4_byte_words_if_little_endian( maximum_points, buffer );
    result =
      fwrite( buffer, sizeof (float), maximum_points, stdout )
      == maximum_points;
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write coordinates.\n" );
  }

  return result;
}



/* Write big-endian copy of data to stdout: */

static int write_output_data( const Node* list,
                              const int maximum_points,
                              float* const buffer ) {
  int result = 1;
  const int variables = list->data->variables;
  int variable = 0;

  for ( variable = 0; result && variable < variables; ++variable ) {
    const Node* node = 0;

    for ( node = list; result && node; node = node->next ) {
      const OutputData* const output_data = node->data;
      const int timesteps = output_data->timesteps;
      const int points = output_data->points;
      const int timesteps_points = timesteps * points;
      const float* const input = output_data->data;
      const int variable_offset = variable * timesteps_points;
      int timestep = 0;

      DEBUG( fprintf( stderr,
                      "node points = %d, timesteps = %d, variables = %d, "
                      "maximum_points = %d, input[ 0 ] = %f\n",
                      points, timesteps, variables, maximum_points,
                      input[ 0 ] ); )
      assert( points <= maximum_points );
      assert( output_data->variables == variables );

      for ( timestep = 0; result && timestep < timesteps; ++timestep ) {
        float* copy = buffer;
        const int timestep_offset = variable_offset + timestep * points;
        int point = 0;

        for ( point = 0; point < points; ++point ) {
          const int input_index = timestep_offset + point;
          assert( copy < buffer + maximum_points );
          assert( IN_RANGE( input_index, 0,
                  output_data->variables * output_data->timesteps *
                    output_data->points - 1 ) );

          DEBUG( if ( copy == buffer )
                   fprintf( stderr, "input[ %d ] = %f\n",
                            input_index, input[ input_index ] ); )

          *copy++ = input[ input_index ];
        }

        DEBUG( fprintf( stderr, "last copy = %f\n", *( copy - 1 ) ); )

        for ( ; point < maximum_points; ++point ) {
          assert( copy < buffer + maximum_points );
          *copy++ = MISSING;
        }

        DEBUG( fprintf( stderr, "buffer = [%f %f]\n",
                        buffer[ 0 ], buffer[ maximum_points - 1 ] );)

        reverse_4_byte_words_if_little_endian( maximum_points, buffer );
        result =
          fwrite( buffer, sizeof (float), maximum_points, stdout )
          == maximum_points;
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\nFailed to write all data.\n" );
  }

  return result;
}



/* Read and convert CBOFS data: */

static OutputData* read_CBOFS( const Arguments* arguments, FILE* file ) {
  OutputData* result = 0;
  CBOFSData* data = read_CBOFS_header( arguments, file );

  if ( data ) {

    if ( read_CBOFS_data( arguments, file, data ) ) {
      result = convert_CBOFS( arguments, data );
    }

    deallocate_cbofs_data( data ), data = 0;
  }

  return result;
}



/* Read/check DODS header for CBOFS file that looks like a subset of:
 Dataset {
    Float64 s_rho[s_rho = 1];
    Float64 h[eta_rho = 291][xi_rho = 332];
    Float64 lon_rho[eta_rho = 291][xi_rho = 332];
    Float64 lat_rho[eta_rho = 291][xi_rho = 332];
    Float64 angle[eta_rho = 291][xi_rho = 332];
    Float64 mask_rho[eta_rho = 291][xi_rho = 332];
    Float64 mask_u[eta_u = 291][xi_u = 331];
    Float64 mask_v[eta_v = 290][xi_v = 332];
    Float64 ocean_time[ocean_time = 7];
    Float32 zeta[ocean_time = 7][eta_rho = 291][xi_rho = 332];
    Float32 Pair[ocean_time = 7][eta_rho = 291][xi_rho = 332];
    Float32 Uwind[ocean_time = 7][eta_rho = 291][xi_rho = 332];
    Float32 Vwind[ocean_time = 7][eta_rho = 291][xi_rho = 332];
    Float32 u[ocean_time = 7][s_rho = 1][eta_u = 291][xi_u = 331];
    Float32 v[ocean_time = 7][s_rho = 1][eta_v = 290][xi_v = 332];
    Float32 w[ocean_time = 7][s_w = 2][eta_rho = 291][xi_rho = 332];
    Float32 temp[ocean_time = 7][s_rho = 1][eta_rho = 291][xi_rho = 332];
    Float32 salt[ocean_time = 7][s_rho = 1][eta_rho = 291][xi_rho = 332];
} NOAA%2fCBOFS%2fMODELS%2f201306%2fnos%2ecbofs%2efields%2enowcast%2e20130617...

Data:
*/

static CBOFSData* read_CBOFS_header( const Arguments* arguments, FILE* file ) {
  const char* const variable = arguments ? arguments->variable : "";
  const int has_depth =
    ! strcmp( variable, "water_temperature" ) ||
    ! strcmp( variable, "oxygen" ) ||
    ! strcmp( variable, "salinity" ) ||
    ! strcmp( variable, "current" );
  const int is_vector =
    ! strcmp( variable, "wind" ) || ! strcmp( variable, "current" );
  CBOFSData* result = 0;
  int timesteps = 0;
  int rows = 0;
  int columns = 0;
  int ok = 0;
  assert( arguments ); assert( file );

  ok = read_matched_line( file, "Dataset {\n" );

  if ( ok ) {
    int h_rows = 0;
    int h_columns = 0;
    Line line = "";
    memset( line, 0, sizeof line );

    /* If data has depth then read s_rho and h: */

    if ( has_depth ) {
      ok = read_matched_line( file, "Float64 s_rho[s_rho = 1];\n" );
      ok = ok && read_dimension2_line( file,
                                       "Float64 "
                                       "h[eta_rho = %d][xi_rho = %d];\n",
                                       &h_rows, &h_columns );
    }

    /* Read lon_rho and lat_rho: */

    if ( ok ) {
      ok = read_dimension2_line( file,
                                 "Float64 "
                                 "lon_rho[eta_rho = %d][xi_rho = %d];\n",
                                 &rows, &columns );
      ok = ok && ( ! has_depth || ( rows == h_rows && columns == h_columns ) );

      if ( ok ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float64 lat_rho[eta_rho = %d][xi_rho = %d];\n",
                  rows, columns );
        ok = read_matched_line( file, line );
      }
    }

    /* If data is a vector then read angle: */

    if ( ok && is_vector) {
      snprintf( line, sizeof line / sizeof *line,
                "Float64 angle[eta_rho = %d][xi_rho = %d];\n",
                rows, columns );
      ok = read_matched_line( file, line );
    }

    /* Read mask_rho: */

    if ( ok ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float64 mask_rho[eta_rho = %d][xi_rho = %d];\n",
                rows, columns );
      ok = read_matched_line( file, line );
    }

    /* If data is current then read mask_u, mask_v: */

    if ( ok && ! strcmp( variable, "current" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float64 mask_u[eta_u = %d][xi_u = %d];\n",
                rows, columns - 1 );
      ok = read_matched_line( file, line );

      if ( ok ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float64 mask_v[eta_v = %d][xi_v = %d];\n",
                  rows - 1, columns );
        ok = read_matched_line( file, line );
      }
    }

    /* Read ocean_time: */

    ok = ok && read_dimension1_line( file,
                                     "Float64 "
                                     "ocean_time[ocean_time = %d];\n",
                                     &timesteps ) && timesteps > 0;

    /* If data has depth then read zeta: */

    if ( ok && has_depth ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 "
                "zeta[ocean_time = %d][eta_rho = %d][xi_rho = %d];\n",
                timesteps, rows, columns );
      ok = read_matched_line( file, line );
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "air_pressure" ) ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float32 "
                  "Pair[ocean_time = %d][eta_rho = %d][xi_rho = %d];\n",
                  timesteps, rows, columns );
        ok = read_matched_line( file, line );
      } else if ( ! strcmp( variable, "wind" ) ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float32 "
                  "Uwind[ocean_time = %d][eta_rho = %d][xi_rho = %d];\n",
                  timesteps, rows, columns );
        ok = read_matched_line( file, line );

        if ( ok ) {
          snprintf( line, sizeof line / sizeof *line,
                    "Float32 "
                    "Vwind[ocean_time = %d][eta_rho = %d][xi_rho = %d];\n",
                    timesteps, rows, columns );
          ok = read_matched_line( file, line );
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float32 u[ocean_time = %d][s_rho = 1]"
                  "[eta_u = %d][xi_u = %d];\n",
                  timesteps, rows, columns - 1 );
        ok = read_matched_line( file, line );

        if ( ok ) {
          snprintf( line, sizeof line / sizeof *line,
                    "Float32 v[ocean_time = %d][s_rho = 1]"
                    "[eta_v = %d][xi_v = %d];\n",
                    timesteps, rows - 1, columns );
          ok = read_matched_line( file, line );
        }

        if ( ok ) {
          snprintf( line, sizeof line / sizeof *line,
                    "Float32 w[ocean_time = %d][s_w = 2]"
                    "[eta_rho = %d][xi_rho = %d];\n",
                    timesteps, rows, columns );
          ok = read_matched_line( file, line );
        }
      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float32 "
                  "temp[ocean_time = %d][s_rho = 1]"
                  "[eta_rho = %d][xi_rho = %d];\n",
                  timesteps, rows, columns );
        ok = read_matched_line( file, line );
      } else if ( ! strcmp( variable, "salinity" ) ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float32 "
                  "salt[ocean_time = %d][s_rho = 1]"
                  "[eta_rho = %d][xi_rho = %d];\n",
                  timesteps, rows, columns );
        ok = read_matched_line( file, line );
      } else if ( ! strcmp( variable, "oxygen" ) ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float32 "
                  "oxygen[ocean_time = %d][s_rho = 1]"
                  "[eta_rho = %d][xi_rho = %d];\n",
                  timesteps, rows, columns );
        ok = read_matched_line( file, line );
      }
    }

    /* Read end of header: */

    DEBUG( fprintf( stderr, "before end of header: ok = %d\n", ok ); )
#if 0
    ok = ok && fscanf( file, "} %*s\n\nData:\n" );
#else
    ok = ok &&
         read_ignored_line( file ) &&
         read_ignored_line( file ) &&
         read_ignored_line( file );
#endif
    DEBUG( fprintf( stderr, "after end of header: ok = %d\n", ok ); )
  }

  if ( ! ok ) {
    fprintf( stderr, "\nInvalid input DODS header.\n" );
  } else {
    result = NEW( CBOFSData, 1 );

    if ( result ) {
      result->timesteps = timesteps;
      result->rows = rows;
      result->columns = columns;
      DEBUG( fprintf( stderr,
                      "timesteps = %d, rows = %d, columns = %d\n",
                      result->timesteps, result->rows, result->columns ); )
    }
  }

  DEBUG( fprintf( stderr, "read_CBOFS_header result = %p\n", result ); )
  return result;
}



/* Read CBOFS data: */

static int read_CBOFS_data( const Arguments* arguments, FILE* file,
                            CBOFSData* data ) {
  int ok = 1;
  assert( arguments ); assert( data ); assert( file );
  assert( data->timesteps > 0 );
  assert( data->rows > 0 ); assert( data->columns > 0 );
  DEBUG( fprintf( stderr,
                  "read_CBOFS_data( data->timsteps = %d, rows = %d, "
                  "columns = %d )\n",
                  data->timesteps, data->rows, data->columns ); )
  {
    const char* const variable = arguments->variable;
    const int has_depth =
      ! strcmp( variable, "water_temperature" ) ||
      ! strcmp( variable, "oxygen" ) ||
      ! strcmp( variable, "salinity" ) ||
      ! strcmp( variable, "current" );
    const int is_vector =
      ! strcmp( variable, "wind" ) || ! strcmp( variable, "current" );
    const int timesteps = data->timesteps;
    const int rows = data->rows;
    const int columns = data->columns;
    const int rows_columns = rows * columns;
    const int rows_1_columns = ( rows - 1 ) * columns;
    const int rows_columns_1 = rows * ( columns - 1 );
    const int timesteps_rows_columns = timesteps * rows_columns;
    const int timesteps_rows_1_columns = timesteps * rows_1_columns;
    const int timesteps_rows_columns_1 = timesteps * rows_columns_1;

    DEBUG( fprintf( stderr, "has_depth = %d, is_vector = %d\n",
                    has_depth, is_vector ); )

    /* If data has depth then read s_rho and h: */

    if ( has_depth ) {
      ok = read_double_data( file, 1.0, 0.0, -1.0, 0.0, 1, &data->s );

      if ( ok ) {
        data->h =
          read_new_double_data( file, 1.0, 0.0, 0.0, 1.1e4, rows_columns );
        ok = data->h != 0;
      }
    }

    /* Read lon_rho and lat_rho: */

    if ( ok ) {
      data->longitude =
        read_new_double_data( file, 1.0, 0.0, -180.0, 180.0, rows_columns );
      ok = data->longitude != 0;
    }

    if ( ok ) {
      data->latitude =
        read_new_double_data( file, 1.0, 0.0, -90.0, 90.0, rows_columns );
      ok = data->latitude != 0;
    }

    /* If data is a vector then read angle: */

    if ( ok && is_vector) {
      data->angle =
        read_new_double_data( file, 1.0, 0.0, -2.0 * M_PI, 2.0 * M_PI,
                              rows_columns );
      ok = data->angle != 0;
    }

    /* Read mask_rho: */

    if ( ok ) {
      data->mask =
        read_new_double_data( file, 1.0, 0.0, 0.0, 1.0, rows_columns );
      ok = data->mask != 0;
    }

    /* If data is current then read mask_u, mask_v: */

    if ( ok && ! strcmp( variable, "current" ) ) {
      data->mask_u =
        read_new_double_data( file, 1.0, 0.0, 0.0, 1.0, rows_columns_1 );
      ok = data->mask_u != 0;

      if ( ok ) {
        data->mask_v =
          read_new_double_data( file, 1.0, 0.0, 0.0, 1.0, rows_1_columns );
        ok = data->mask_v != 0;
      }
    }

    /* Read ocean_time: */

    if ( ok ) {
      const int yyyymmddhh = base_yyyymmddhh( arguments );
      const double start_seconds =
        seconds_difference( arguments->yyyymmddhh, yyyymmddhh );
      const double end_seconds =
        start_seconds + arguments->hours * SECONDS_PER_HOUR - 1.0;
      const double scale = timestep_scale( arguments->source );
      data->seconds =
        read_new_double_data( file, scale, SECONDS_PER_HALF_HOUR,
                              start_seconds, end_seconds,
                              timesteps );
      ok = data->seconds != 0;
    }

    /* If data has depth then read zeta: */

    if ( ok && has_depth ) {
      data->msl =
        read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                             timesteps_rows_columns );
      ok = data->msl != 0;
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "air_pressure" ) ) {
        data->air_pressure =
          read_new_float_data( file, 1.0, 0.0, 500.0, 1500.0,
                               timesteps_rows_columns );
        ok = data->air_pressure != 0;
      } else if ( ! strcmp( variable, "wind" ) ) {
        data->wind_u =
          read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                               timesteps_rows_columns );
        ok = data->wind_u != 0;

        if ( ok ) {
          data->wind_v =
            read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                                 timesteps_rows_columns );
          ok = data->wind_v != 0;
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        data->current_u =
          read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                               timesteps_rows_columns_1 );
        ok = data->current_u != 0;

        if ( ok ) {
          data->current_v =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_rows_1_columns);
          ok = data->current_v != 0;
        }

        if ( ok ) {
          data->current_w =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps * 2 * rows_columns );
          ok = data->current_w != 0;
        }
      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        data->water_temperature =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_rows_columns );
        ok = data->water_temperature != 0;
      } else if ( ! strcmp( variable, "oxygen" ) ) {
        data->oxygen =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1e30,
                               timesteps_rows_columns );
        ok = data->oxygen != 0;
      } else {
        assert( ! strcmp( variable, "salinity" ) );
        data->salinity =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_rows_columns );
        ok = data->salinity != 0;
      }
    }
  }

  DEBUG( fprintf( stderr, "read_CBOFS_data result = %d\n", ok ); )
  return ok;
}



/* Convert CBOFS data: */

static OutputData* convert_CBOFS( const Arguments* arguments,
                                  CBOFSData* data ) {
  OutputData* result = 0;
  const int subset_points =
    spatially_subset( arguments->longitude_minimum,
                      arguments->longitude_maximum,
                      arguments->latitude_minimum,
                      arguments->latitude_maximum,
                      data->rows * data->columns,
                      sizeof data->mask[ 0 ],
                      data->longitude,
                      data->latitude,
                      data->mask );
  DEBUG( fprintf( stderr, "subset points = %d\n", subset_points ); )

  if ( subset_points ) {
    result = NEW( OutputData, 1 );

    if ( result ) {
      const size_t length = sizeof (Name) / sizeof (char);
      const int timesteps = data->timesteps;
      int timestep = 0;
      result->points = subset_points;

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        const double seconds = data->seconds[ timestep ];

        if ( seconds != MISSING ) {
          result->timesteps += 1;
        }
      }

      assert( result->timesteps > 0 );

      if ( data->air_pressure ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "air_pressure", length );
        strncpy( result->units[ 0 ], "hPa", length );
      } else if ( data->wind_u ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "wind_u", length );
        strncpy( result->names[ 1 ], "wind_v", length );
        strncpy( result->units[ 0 ], "m/s", length );
        strncpy( result->units[ 1 ], "m/s", length );
      } else if ( data->current_u ) {
        result->variables = 4;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "current_u", length );
        strncpy( result->names[ 2 ], "current_v", length );
        strncpy( result->names[ 3 ], "current_w", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "m/s", length );
        strncpy( result->units[ 2 ], "m/s", length );
        strncpy( result->units[ 3 ], "m/s", length );
      } else if ( data->water_temperature ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "water_temperature", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "C", length );
      } else if ( data->oxygen ){
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "oxygen", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "mmol/m3", length );
      } else {
        assert( data->salinity );
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "salinity", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "PSU", length );
      }

      result->yyyymmddhh = NEW( int, timesteps );
      result->longitudes = result->yyyymmddhh ? NEW( float, subset_points) : 0;
      result->latitudes  = result->longitudes ? NEW( float, subset_points) : 0;
      result->data = result->latitudes ?
        NEW( float, result->variables * timesteps * subset_points) : 0;

      if ( result->data ) {
        float* output_data = result->data;
        const int rows = data->rows;
        const int columns = data->columns;
        const int rows_columns = rows * columns;
        const double* const mask = data->mask;
        const double* const longitudes = data->longitude;
        const double* const latitudes = data->latitude;
        const float* input_data =
          data->air_pressure ? data->air_pressure
          : data->wind_u ? data->wind_u
          : data->water_temperature ? data->water_temperature
          : data->salinity ? data->salinity
          : data->oxygen ? data->oxygen
          : 0;
        const float* input_data_2 = data->wind_v;
        const int yyyymmddhh0 = base_yyyymmddhh( arguments );
        int variable_offset = 0;
        int output_timesteps = 0;

        assert( input_data_2 == 0 ||
                ( input_data == data->wind_u && input_data_2 == data->wind_v));

        /* Copy timestamps and count number of output timesteps: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            const int yyyymmddhh = convert_seconds(yyyymmddhh0, (int) seconds);
            assert( output_timesteps < 1 ||
                    yyyymmddhh > result->yyyymmddhh[ output_timesteps - 1 ] );
            result->yyyymmddhh[ output_timesteps ] = yyyymmddhh;
            ++output_timesteps;
            variable_offset += subset_points;
          }
        }

        /* Copy longitudes: */

        copy_unmasked_values( rows_columns, mask, longitudes,
                              result->longitudes );
        assert( IN_RANGE( result->longitudes[ 0 ], -180.0, 180.0 ) );
        assert( IN_RANGE( result->longitudes[ subset_points - 1],
                          -180.0, 180.0 ) );

        /* Copy latitudes: */

        copy_unmasked_values(rows_columns, mask, latitudes, result->latitudes);
        assert( IN_RANGE( result->latitudes[ 0 ], -90.0, 90.0 ) );
        assert( IN_RANGE( result->latitudes[ subset_points - 1 ],
                          -90.0, 90.0 ) );

        /* Copy depth if present: */

        if ( data->s ) {
          assert( data->h ); assert( data->msl );

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_depth( rows_columns, mask, data->s, data->h,
                          data->msl + timestep * rows_columns,
                          output_data );
              assert( IN_RANGE( output_data[ 0 ], -3e5, 0.0 ) );
              assert( IN_RANGE( output_data[ subset_points - 1 ],
                                -3e5, 0.0 ) );
              output_data += subset_points;
            }
          }
        }

        /* Copy data variable/components: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {

            if ( input_data_2 ) { /* Convert and copy wind: */
              copy_wind( rows_columns, mask, data->angle,
                         input_data, input_data_2,
                         output_data, output_data + variable_offset );
              output_data  += subset_points;
              input_data   += rows_columns;
              input_data_2 += rows_columns;
            } else if ( input_data ) { /* Copy scalar: */
              copy_unmasked_values2( rows_columns, mask, input_data,
                                     output_data );
              output_data += subset_points;
              input_data  += rows_columns;
            } else { /* Copy current_u/v/w for this timestep: */
              assert( data->current_u );
              copy_current( timestep, data,
                            output_data,
                            output_data + variable_offset,
                            output_data + 2 * variable_offset );
              output_data += subset_points;
            }
          }
        }
      }
    }
  }

  return result;
}



/* Read and convert GBOFS data: */

static OutputData* read_GBOFS( const Arguments* arguments, FILE* file ) {
  OutputData* result = 0;
  GBOFSData* data = read_GBOFS_header( arguments, file );

  if ( data ) {

    if ( read_GBOFS_data( arguments, file, data ) ) {
      result = convert_GBOFS( arguments, data );
    }

    deallocate_gbofs_data( data ), data = 0;
  }

  return result;
}


/* Read/check DODS header for GBOFS file that looks like a subset of:
Dataset {
    Float32 time[time = 25];
    Float32 lon[ny = 134][nx = 73];
    Float32 lat[ny = 134][nx = 73];
    Float32 mask[ny = 134][nx = 73];
    Float32 depth[ny = 134][nx = 73];
    Float32 sigma[sigma = 1];
    Float32 zeta[time = 25][ny = 134][nx = 73];
    Float32 air_u[time = 25][ny = 134][nx = 73];
    Float32 air_v[time = 25][ny = 134][nx = 73];
    Float32 u[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 v[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 w[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 temp[time = 25][sigma = 1][ny = 134][nx = 73];
    Float32 salt[time = 25][sigma = 1][ny = 134][nx = 73];
} NOAA%2fGBOFS%2fMODELS%2f201309%2fnos%2egbofs%2efields%2enowcast%2e20130901...

Data:
*/

static GBOFSData* read_GBOFS_header( const Arguments* arguments, FILE* file ) {
  const char* const variable = arguments ? arguments->variable : "";
  const int has_w = ! strcmp( arguments->source, "nyofs" );
  const int has_depth = strcmp( variable, "wind" ) != 0;
  GBOFSData* result = 0;
  int timesteps = 0;
  int file_base_yyyymmddhh = 0;
  int rows = 0;
  int columns = 0;
  int ok = 0;
  Line line = "";
  memset( line, 0, sizeof line );
  assert( arguments ); assert( file );

  ok = read_matched_line( file, "Dataset {\n" );

  /* Read time: */

  ok = ok && read_dimension1_line( file, "Float32 time[time = %d];\n",
                                   &timesteps );

  /* Read lon and its dimensions: */

  ok = ok &&
    read_dimension2_line( file, "Float32 lon[ny = %d][nx = %d];\n",
                          &rows, &columns );

  /* Read lat: */

  snprintf( line, sizeof line / sizeof *line,
            "Float32 lat[ny = %d][nx = %d];\n", rows, columns );
  ok = ok && read_matched_line( file, line );

  /* Read mask: */

  snprintf( line, sizeof line / sizeof *line,
            "Float32 mask[ny = %d][nx = %d];\n", rows, columns );
  ok = ok && read_matched_line( file, line );

  /* If data has depth then read depth, sigma and zeta: */

  if ( has_depth ) {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 depth[ny = %d][nx = %d];\n", rows, columns );
    ok = ok && read_matched_line( file, line );
    ok = ok && read_matched_line( file, "Float32 sigma[sigma = 1];\n" );
    snprintf( line, sizeof line / sizeof *line,
              "Float32 zeta[time = %d][ny = %d][nx = %d];\n",
              timesteps, rows, columns );
    ok = ok && read_matched_line( file, line );
  }

  /* Read main variable: */

  if ( ok ) {

    if ( ! strcmp( variable, "wind" ) ) { /* Read air_u, air_v: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 air_u[time = %d][ny = %d][nx = %d];\n",
                timesteps, rows, columns );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 air_v[time = %d][ny = %d][nx = %d];\n",
                timesteps, rows, columns );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "current" ) ) { /* Read u, v, w: */
      assert ( ! strcmp( variable, "current" ) );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 u[time = %d][sigma = 1][ny = %d][nx = %d];\n",
                timesteps, rows, columns );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 v[time = %d][sigma = 1][ny = %d][nx = %d];\n",
                timesteps, rows, columns );
      ok = ok && read_matched_line( file, line );

      if ( ok && has_w ) {
        snprintf( line, sizeof line / sizeof *line,
                  "Float32 w[time = %d][sigma = 1][ny = %d][nx = %d];\n",
                  timesteps, rows, columns );
        ok = ok && read_matched_line( file, line );
      }
    } else if ( ! strcmp( variable, "water_temperature" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 temp[time = %d][sigma = 1][ny = %d][nx = %d];\n",
                timesteps, rows, columns );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "salinity" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 salt[time = %d][sigma = 1][ny = %d][nx = %d];\n",
                timesteps, rows, columns );
      ok = read_matched_line( file, line );
    } else {
      fprintf( stderr, "Invalid variable '%s'.\n", variable );
      ok = 0;
    }

    /* Read end of header: */

    DEBUG( fprintf( stderr, "before end of header: ok = %d\n", ok ); )

    if ( ok ) {

      if ( base_yyyymmddhh0( arguments->source ) == 0 ) {
        int end_yyyymmddhh = arguments->yyyymmddhh;
        increment_yyyymmddhh( &end_yyyymmddhh, arguments->hours );
        file_base_yyyymmddhh = read_file_base_yyyymmddhh( file );
        ok = is_valid_yyyymmddhh( file_base_yyyymmddhh ) &&
             file_base_yyyymmddhh < end_yyyymmddhh &&
             read_ignored_line( file ) &&
             read_ignored_line( file );
      } else {
        ok = read_ignored_line( file ) &&
             read_ignored_line( file ) &&
             read_ignored_line( file );
      }
    }

    DEBUG( fprintf( stderr, "after end of header: ok = %d\n", ok ); )
  }

  if ( ! ok ) {

    if ( file_base_yyyymmddhh == 0 ) {
      fprintf( stderr, "\nInvalid input DODS header.\n" );
    }

  } else {
    result = NEW( GBOFSData, 1 );

    if ( result ) {
      result->timesteps = timesteps;
      result->file_base_yyyymmddhh = file_base_yyyymmddhh;
      result->rows = rows;
      result->columns = columns;
      DEBUG( fprintf( stderr, "timesteps = %d, file_base_yyyymmddhh = %d, "
                      "rows = %d, columns = %d\n",
                      result->timesteps, result->file_base_yyyymmddhh,
                      result->rows, result->columns ); )
    }
  }

  DEBUG( fprintf( stderr, "read_GBOFS_header result = %p\n", result ); )
  return result;
}



/* Read GBOFS data: */

static int read_GBOFS_data( const Arguments* arguments, FILE* file,
                            GBOFSData* data ) {
  int ok = 1;
  assert( arguments ); assert( data ); assert( file );
  assert( data->timesteps > 0 );
  assert( data->rows > 0 ); assert( data->columns > 0 );
  DEBUG( fprintf( stderr,
                  "read_GBOFS_data: timesteps = %d, rows = %d, columns = %d\n",
                  data->timesteps, data->rows, data->columns ); )
  {
    const char* const variable = arguments->variable;
    const int has_w = ! strcmp( arguments->source, "nyofs" );
    const int has_depth = strcmp( variable, "wind" ) != 0;
    const int timesteps = data->timesteps;
    const int rows = data->rows;
    const int columns = data->columns;
    const int rows_columns = rows * columns;
    const int timesteps_rows_columns = timesteps * rows_columns;

    DEBUG( fprintf( stderr, "has_depth = %d\n", has_depth ); )
    assert( ! data->file_base_yyyymmddhh || data->timesteps == 1 );

    /* Read time: */

    {
      const int yyyymmddhh = base_yyyymmddhh( arguments );
      const double start_seconds =
        data->file_base_yyyymmddhh ? 0.0
        : seconds_difference( arguments->yyyymmddhh, yyyymmddhh );
      const double end_seconds =
        start_seconds + arguments->hours * SECONDS_PER_HOUR - 1.0;
      const int scale = timestep_scale( arguments->source );
      data->seconds =
        read_new_float_data( file, scale, SECONDS_PER_HALF_HOUR,
                             start_seconds, end_seconds, timesteps );
      ok = data->seconds != 0;
    }

    /* Read lon_rho and lat_rho: */

    if ( ok ) {
      data->longitude =
        read_new_float_data( file, 1.0, 0.0, -180.0, 180.0, rows_columns );
      ok = data->longitude != 0;
    }

    if ( ok ) {
      data->latitude =
        read_new_float_data( file, 1.0, 0.0, -90.0, 90.0, rows_columns );
      ok = data->latitude != 0;
    }

    /* Read mask: */

    if ( ok ) {
      data->mask = read_new_float_data(file, 1.0, 0.0, 0.0, 1.0, rows_columns);
      ok = data->mask != 0;
    }

    /* If data has depth then read depth, sigma and zeta: */

    if ( ok && has_depth ) {
      data->h = read_new_float_data( file, 1.0, 0.0, 0.0, 1.1e4, rows_columns);
      ok = data->h != 0;

      if ( ok ) {
        ok = read_float_data( file, 1.0, -1.0, -1.0, 0.0, 1, &data->s );

        if ( ok ) {
          data->msl =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_rows_columns );
          ok = data->msl != 0;
        }
      }
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "wind" ) ) {
        data->wind_u =
          read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                               timesteps_rows_columns );
        ok = data->wind_u != 0;

        if ( ok ) {
          data->wind_v =
            read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                                 timesteps_rows_columns );
          ok = data->wind_v != 0;
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        data->current_u =
          read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                               timesteps_rows_columns );
        ok = data->current_u != 0;

        if ( ok ) {
          data->current_v =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_rows_columns);
          ok = data->current_v != 0;
        }

        if ( ok && has_w ) {
          data->current_w =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_rows_columns );
          ok = data->current_w != 0;
        }
      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        data->water_temperature =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_rows_columns );
        ok = data->water_temperature != 0;
      } else {
        assert( ! strcmp( variable, "salinity" ) );
        data->salinity =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_rows_columns );
        ok = data->salinity != 0;
      }
    }
  }

  DEBUG( fprintf( stderr, "read_GBOFS_data result = %d\n", ok ); )
  return ok;
}



/* Convert GBOFS data: */

static OutputData* convert_GBOFS( const Arguments* arguments,
                                  GBOFSData* data ) {
  OutputData* result = 0;
  const int subset_points =
    spatially_subset( arguments->longitude_minimum,
                      arguments->longitude_maximum,
                      arguments->latitude_minimum,
                      arguments->latitude_maximum,
                      data->rows * data->columns,
                      sizeof data->mask[ 0 ],
                      data->longitude,
                      data->latitude,
                      data->mask );
  const int has_w = ! strcmp( arguments->source, "nyofs" );
  DEBUG( fprintf( stderr, "subset points = %d\n", subset_points ); )

  if ( subset_points ) {
    result = NEW( OutputData, 1 );

    if ( result ) {
      const size_t length = sizeof (Name) / sizeof (char);
      const int timesteps = data->timesteps;
      int timestep = 0;
      result->points = subset_points;

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        const double seconds = data->seconds[ timestep ];

        if ( seconds != MISSING ) {
          result->timesteps += 1;
        }
      }

      assert( result->timesteps > 0 );

      if ( data->wind_u ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "wind_u", length );
        strncpy( result->names[ 1 ], "wind_v", length );
        strncpy( result->units[ 0 ], "m/s", length );
        strncpy( result->units[ 1 ], "m/s", length );
      } else if ( data->current_u ) {
        result->variables = 3;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "current_u", length );
        strncpy( result->names[ 2 ], "current_v", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "m/s", length );
        strncpy( result->units[ 2 ], "m/s", length );

        if ( has_w ) {
          result->variables = 4;
          strncpy( result->names[ 3 ], "current_w", length );
          strncpy( result->units[ 3 ], "m/s", length );
        }
      } else if ( data->water_temperature ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "water_temperature", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "C", length );
      } else {
        assert( data->salinity );
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "salinity", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "PSU", length );
      }

      result->yyyymmddhh = NEW( int, timesteps );
      result->longitudes = result->yyyymmddhh ? NEW( float, subset_points) : 0;
      result->latitudes  = result->longitudes ? NEW( float, subset_points) : 0;
      result->data = result->latitudes ?
        NEW( float, result->variables * timesteps * subset_points ) : 0;

      if ( result->data ) {
        float* output_data = result->data;
        const int rows = data->rows;
        const int columns = data->columns;
        const int rows_columns = rows * columns;
        const float* const mask = data->mask;
        const float* const longitudes = data->longitude;
        const float* const latitudes = data->latitude;
        const float* input_data =
          data->wind_u ? data->wind_u
          : data->current_u ? data->current_u
          : data->water_temperature ? data->water_temperature
          : data->salinity ? data->salinity
          : 0;
        const float* input_data_2 =
          data->wind_v ? data->wind_v : data->current_v;
        const float* input_data_3 = data->current_w;
        const int yyyymmddhh0 = base_yyyymmddhh( arguments );
        int yyyymmddhh = data->file_base_yyyymmddhh;
        int output_timesteps = 0;

        /* Copy timestamps and count number of output timesteps: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {

          if ( data->file_base_yyyymmddhh ) {
            assert( output_timesteps < 1 ||
                    yyyymmddhh > result->yyyymmddhh[ output_timesteps - 1 ] );
            result->yyyymmddhh[ output_timesteps ] = yyyymmddhh;
            increment_yyyymmddhh( &yyyymmddhh, 1 );
            ++output_timesteps;
          } else {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              yyyymmddhh = convert_seconds( yyyymmddhh0, (int) seconds );
              assert( output_timesteps < 1 ||
                      yyyymmddhh > result->yyyymmddhh[ output_timesteps - 1 ]);
              result->yyyymmddhh[ output_timesteps ] = yyyymmddhh;
              ++output_timesteps;
            }
          }
        }

        /* Copy longitudes: */

        copy_unmasked_float_values( rows_columns, mask, longitudes,
                                    result->longitudes );
        assert( IN_RANGE( result->longitudes[ 0 ], -180.0, 180.0 ) );
        assert( IN_RANGE( result->longitudes[ subset_points - 1],
                          -180.0, 180.0 ) );

        /* Copy latitudes: */

        copy_unmasked_float_values( rows_columns, mask, latitudes,
                                    result->latitudes );
        assert( IN_RANGE( result->latitudes[ 0 ], -90.0, 90.0 ) );
        assert( IN_RANGE( result->latitudes[ subset_points - 1 ],
                          -90.0, 90.0 ) );

        /* Copy depth if present: */

        if ( data->s ) {
          assert( data->h ); assert( data->msl );

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_float_depth( rows_columns, mask, data->s, 0, data->h,
                                data->msl + timestep * rows_columns,
                                output_data );
              assert( IN_RANGE( output_data[ 0 ], -3e5, 0.0 ) );
              assert( IN_RANGE( output_data[ subset_points - 1 ], -3e5, 0.0 ));
              output_data += subset_points;
            }
          }
        }

        /* Copy data variable: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            copy_unmasked_float_values( rows_columns, mask, input_data,
                                        output_data );
            output_data += subset_points;
            input_data  += rows_columns;
          }
        }

        /* Copy additional component of data variable: */

        if ( input_data_2 ) {

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_unmasked_float_values( rows_columns, mask, input_data_2,
                                          output_data );
              output_data += subset_points;
              input_data_2 += rows_columns;
            }
          }
        }

        /* Copy additional component of data variable: */

        if ( input_data_3 ) {

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_unmasked_float_values( rows_columns, mask, input_data_3,
                                          output_data );
              output_data += subset_points;
              input_data_3 += rows_columns;
            }
          }
        }
      }
    }
  }

  return result;
}



/* Read and convert NGOFS data: */

static OutputData* read_NGOFS( const Arguments* arguments, FILE* file ) {
  OutputData* result = 0;
  NGOFSData* data = read_NGOFS_header( arguments, file );

  if ( data ) {

    if ( read_NGOFS_data( arguments, file, data ) ) {
      result = convert_NGOFS( arguments, data );
    }

    deallocate_ngofs_data( data ), data = 0;
  }

  return result;
}


/* Read/check DODS header for NGOFS file that looks like a subset of:
Dataset {
    Float32 lon[node = 90267];
    Float32 lat[node = 90267];
    Float32 lonc[nele = 174474];
    Float32 latc[nele = 174474];
    Float32 siglay[siglay = 1][node = 90267];
    Float32 h[node = 90267];
    Float32 nv[three = 3][nele = 174474];
    Float32 time[time = 7];
    Float32 zeta[time = 7][node = 90267];
    Float32 u[time = 7][siglay = 1][nele = 174474];
    Float32 v[time = 7][siglay = 1][nele = 174474];
    Float32 tauc[time = 7][nele = 174474];
    Float32 temp[time = 7][siglay = 1][node = 90267];
    Float32 salinity[time = 7][siglay = 1][node = 90267];
    Float32 short_wave[time = 7][node = 90267];
    Float32 net_heat_flux[time = 7][node = 90267];
    Float32 uwind_speed[time = 7][nele = 174474];
    Float32 vwind_speed[time = 7][nele = 174474];
    Int32 wet_nodes[time = 7][node = 90267];
    Int32 wet_cells[time = 7][nele = 174474];
} NOAA%2fNGOFS%2fMODELS%2f201210%2fnos%2engofs%2efields%2enowcast%2e\
20121001%2et03z%2enc;

Data:
*/

static NGOFSData* read_NGOFS_header( const Arguments* arguments, FILE* file ) {
  const char* const variable = arguments ? arguments->variable : "";
  const int is_current = strcmp( variable, "current" ) == 0;
  const int has_depth =
    is_current ||
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0;
  const int is_nodal =
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0 ||
    strcmp( variable, "radiation" ) == 0 ||
    strcmp( variable, "net_heat_flux" ) == 0;
  NGOFSData* result = 0;
  int timesteps = 0;
  int points = 0;
  const int nodes = is_current ? NGOFS_nodes : 0;
  int ok = 0;
  Line line = "";
  memset( line, 0, sizeof line );
  assert( arguments ); assert( file );

  ok = read_matched_line( file, "Dataset {\n" );

  /* Read lon and its dimensions: */

  if ( is_nodal ) {
    ok = ok &&
      read_dimension1_line( file, "Float32 lon[node = %d];\n", &points );
  } else {
    ok = ok &&
      read_dimension1_line( file, "Float32 lonc[nele = %d];\n", &points );
  }

  /* Read lat: */

  if ( is_nodal ) {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 lat[node = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  } else {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 latc[nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  if ( has_depth ) { /* Read sigma, depth: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 siglay[siglay = 1][node = %d];\n", count );
    ok = ok && read_matched_line( file, line );

    snprintf( line, sizeof line / sizeof *line,
              "Float32 h[node = %d];\n", count );
    ok = ok && read_matched_line( file, line );
  }

  if ( is_current ) { /* Read nv: */
    snprintf( line, sizeof line / sizeof *line,
              "Int32 nv[three = 3][nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  /* Read time: */

  ok = ok && read_dimension1_line( file, "Float32 time[time = %d];\n",
                                   &timesteps );

  if ( has_depth ) { /* Read zeta: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 zeta[time = %d][node = %d];\n", timesteps, count );
    ok = ok && read_matched_line( file, line );
  }

  /* Read main variable: */

  if ( ok ) {

    if ( ! strcmp( variable, "wind" ) ) { /* Read uwind_speed, vwind_speed: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 uwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 vwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "current" ) ) { /* Read u, v: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 u[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 v[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "water_temperature" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 temp[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "salinity" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 salinity[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "seabed_stress" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 tauc[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "radiation" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 short_wave[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "net_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 net_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else {
      fprintf( stderr, "Invalid variable '%s'.\n", variable );
      ok = 0;
    }

    /* Read mask: */

    if ( is_nodal ) {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_nodes[time = %d][node = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    } else {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_cells[time = %d][nele = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    }

    /* Read end of header: */

    DEBUG( fprintf( stderr, "before end of header: ok = %d\n", ok ); )

    if ( ok ) {
      ok = read_ignored_line( file ) &&
           read_ignored_line( file ) &&
           read_ignored_line( file );
    }

    DEBUG( fprintf( stderr, "after end of header: ok = %d\n", ok ); )
  }

  if ( ! ok ) {
    fprintf( stderr, "\nInvalid input DODS header.\n" );
  } else {
    result = NEW( NGOFSData, 1 );

    if ( result ) {
      result->timesteps = timesteps;
      result->points = points;
      result->nodes = nodes;
      DEBUG( fprintf( stderr, "timesteps = %d, points = %d, nodes = %d\n",
                      result->timesteps, result->points, result->nodes ); )
    }
  }

  DEBUG( fprintf( stderr, "read_NGOFS_header result = %p\n", result ); )
  return result;
}



/* Read NGOFS data: */

static int read_NGOFS_data( const Arguments* arguments, FILE* file,
                            NGOFSData* data ) {
  int ok = 1;
  assert( arguments ); assert( data ); assert( file );
  assert( data->timesteps > 0 );
  assert( data->points > 0 );
  DEBUG( fprintf( stderr,
                  "read_NGOFS_data: timesteps = %d, points = %d\n",
                  data->timesteps, data->points ); )
  {
    const char* const variable = arguments->variable;
    const int is_current = strcmp( variable, "current" ) == 0;
    const int has_depth =
      is_current ||
      strcmp( variable, "water_temperature" ) == 0 ||
      strcmp( variable, "salinity" ) == 0;
    const int timesteps = data->timesteps;
    const int points = data->points;
    const int nodes = data->nodes;
    const int timesteps_points = timesteps * points;

    DEBUG( fprintf( stderr, "has_depth = %d\n", has_depth ); )

    /* Read and convert longitudes: */

    data->longitude =
      read_new_float_data( file, 1.0, -360.0, -180.0, 180.0, points );
    ok = data->longitude != 0;

    /* Read latitudes: */

    if ( ok ) {
      data->latitude =
        read_new_float_data( file, 1.0, 0.0, -90.0, 90.0, points );
      ok = data->latitude != 0;
    }

    /* If data has depth then read siglay and h: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->s = read_new_float_data( file, 1.0, 0.0, -1.0, 0.0, count );
      ok = data->s != 0;

      if ( ok ) {
        data->h = read_new_float_data( file, 1.0, 0.0, 0.0, 1.1e4, count );
        ok = data->h != 0;
      }
    }

    if ( ok && is_current ) { /* Read nv: */
      data->node =
        read_new_int_data( file, 1.0, 0.0, 0.0, NGOFS_nodes - 1, 3 * points );
      ok = data->node != 0;
    }

    /* Read time: */

    if ( ok ) {
      const int yyyymmddhh = base_yyyymmddhh( arguments );
      const double start_seconds =
        seconds_difference( arguments->yyyymmddhh, yyyymmddhh );
      const double end_seconds =
        start_seconds + arguments->hours * SECONDS_PER_HOUR - 1.0;
      const int scale = timestep_scale( arguments->source );
      data->seconds =
        read_new_float_data( file, scale, SECONDS_PER_HALF_HOUR,
                             start_seconds, end_seconds, timesteps );
      ok = data->seconds != 0;
    }

    /* If data has depth then read zeta: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->msl =
        read_new_float_data( file, 1.0, 0.0, -50.0, 50.0, timesteps * count );
      ok = data->msl != 0;
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "wind" ) ) {
        data->wind_u =
          read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                               timesteps_points );
        ok = data->wind_u != 0;

        if ( ok ) {
          data->wind_v =
            read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                                 timesteps_points );
          ok = data->wind_v != 0;
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        data->current_u =
          read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                               timesteps_points );
        ok = data->current_u != 0;

        if ( ok ) {
          data->current_v =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_points);
          ok = data->current_v != 0;
        }

      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        data->water_temperature =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->water_temperature != 0;
      } else if ( ! strcmp( variable, "salinity" ) ) {
        data->salinity =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->salinity != 0;
      } else if ( ! strcmp( variable, "seabed_stress" ) ) {
        data->seabed_stress =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1.0,
                               timesteps_points );
        ok = data->seabed_stress != 0;
      } else if ( ! strcmp( variable, "radiation" ) ) {
        data->radiation =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1500.0,
                               timesteps_points );
        ok = data->radiation != 0;
      } else {
        assert( ! strcmp( variable, "net_heat_flux" ) );
        data->net_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->net_heat_flux != 0;
      }
    }

    /* Read mask: */

    if ( ok ) {
      data->mask =
        read_new_int_data( file, 1.0, 0.0, 0.0, 1.0, timesteps_points );
      ok = data->mask != 0;
    }
  }

  DEBUG( fprintf( stderr, "read_NGOFS_data result = %d\n", ok ); )
  return ok;
}



/* Convert NGOFS data: */

static OutputData* convert_NGOFS( const Arguments* arguments,
                                  NGOFSData* data ) {
  OutputData* result = 0;
  const int subset_points =
    spatially_subset( arguments->longitude_minimum,
                      arguments->longitude_maximum,
                      arguments->latitude_minimum,
                      arguments->latitude_maximum,
                      data->points,
                      sizeof data->mask[ 0 ],
                      data->longitude,
                      data->latitude,
                      data->mask );
  DEBUG( fprintf( stderr, "subset points = %d\n", subset_points ); )

  if ( subset_points ) {
    result = NEW( OutputData, 1 );

    if ( result ) {
      const size_t length = sizeof (Name) / sizeof (char);
      const int timesteps = data->timesteps;
      int timestep = 0;
      result->points = subset_points;

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        const double seconds = data->seconds[ timestep ];

        if ( seconds != MISSING ) {
          result->timesteps += 1;
        }
      }

      assert( result->timesteps > 0 );

      if ( data->wind_u ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "wind_u", length );
        strncpy( result->names[ 1 ], "wind_v", length );
        strncpy( result->units[ 0 ], "m/s", length );
        strncpy( result->units[ 1 ], "m/s", length );
      } else if ( data->current_u ) {
        result->variables = 3;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "current_u", length );
        strncpy( result->names[ 2 ], "current_v", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "m/s", length );
        strncpy( result->units[ 2 ], "m/s", length );
      } else if ( data->water_temperature ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "water_temperature", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "C", length );
      } else if ( data->salinity ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "salinity", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "PSU", length );
      } else if ( data->seabed_stress ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "seabed_stress", length );
        strncpy( result->units[ 0 ], "m2/s2", length );
      } else if ( data->radiation ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "radiation", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else {
        assert( data->net_heat_flux );
        result->variables = 1;
        strncpy( result->names[ 0 ], "net_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      }

      result->yyyymmddhh = NEW( int, timesteps );
      result->longitudes = result->yyyymmddhh ? NEW( float, subset_points) : 0;
      result->latitudes  = result->longitudes ? NEW( float, subset_points) : 0;
      result->data = result->latitudes ?
        NEW( float, result->variables * timesteps * subset_points ) : 0;

      if ( result->data ) {
        float* output_data = result->data;
        const int points = data->points;
        const int nodes = data->nodes;
        const float* const mask = data->mask;
        const float* const longitudes = data->longitude;
        const float* const latitudes = data->latitude;
        const float* input_data =
          data->wind_u ? data->wind_u
          : data->current_u ? data->current_u
          : data->water_temperature ? data->water_temperature
          : data->salinity ? data->salinity
          : data->seabed_stress ? data->seabed_stress
          : data->radiation ? data->radiation
          : data->net_heat_flux ? data->net_heat_flux
          : 0;
        const float* input_data_2 =
          data->wind_v ? data->wind_v : data->current_v;
        const int yyyymmddhh0 = base_yyyymmddhh( arguments );
        int output_timesteps = 0;

        /* Copy timestamps and count number of output timesteps: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            const int yyyymmddhh =
              convert_seconds( yyyymmddhh0, (int) seconds );
            assert( output_timesteps < 1 ||
                    yyyymmddhh > result->yyyymmddhh[ output_timesteps - 1 ] );
            result->yyyymmddhh[ output_timesteps ] = yyyymmddhh;
            ++output_timesteps;
          }
        }

        /* Copy longitudes: */

        copy_unmasked_float_values( points, mask, longitudes,
                                    result->longitudes );
        assert( IN_RANGE( result->longitudes[ 0 ], -180.0, 180.0 ) );
        assert( IN_RANGE( result->longitudes[ subset_points - 1],
                          -180.0, 180.0 ) );

        /* Copy latitudes: */

        copy_unmasked_float_values( points, mask, latitudes,
                                    result->latitudes );
        assert( IN_RANGE( result->latitudes[ 0 ], -90.0, 90.0 ) );
        assert( IN_RANGE( result->latitudes[ subset_points - 1 ],
                          -90.0, 90.0 ) );

        /* Copy depth if present: */

        if ( data->s ) {
          assert( data->h ); assert( data->msl );

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {

              if ( nodes ) {
                copy_float_depth3( points, data->node,
                                   mask, data->s, data->h,
                                   data->msl + timestep * nodes,
                                   output_data );
              } else {
                copy_float_depth( points, mask, 0.0, data->s, data->h,
                                  data->msl + timestep * points,
                                  output_data );
              }
              assert( IN_RANGE( output_data[ 0 ], -3e5, 0.0 ) );
              assert( IN_RANGE( output_data[ subset_points - 1 ], -3e5, 0.0 ));
              output_data += subset_points;
            }
          }
        }

        /* Copy data variable: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            copy_unmasked_float_values( points, mask, input_data, output_data);
            output_data += subset_points;
            input_data  += points;
          }
        }

        /* Copy additional component of data variable: */

        if ( input_data_2 ) {

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_unmasked_float_values( points, mask, input_data_2,
                                          output_data );
              output_data += subset_points;
              input_data_2 += points;
            }
          }
        }
      }
    }
  }

  return result;
}



/* Read and convert NGOFS2 data: */

static OutputData* read_NGOFS2( const Arguments* arguments, FILE* file ) {
  OutputData* result = 0;
  NGOFS2Data* data = read_NGOFS2_header( arguments, file );

  if ( data ) {

    if ( read_NGOFS2_data( arguments, file, data ) ) {
      result = convert_NGOFS2( arguments, data );
    }

    deallocate_ngofs2_data( data ), data = 0;
  }

  return result;
}


/*
 Input data from NGOFS2 DODS files whose header looks like a subset of:

Dataset {
    Float32 lon[node = 303714];
    Float32 lat[node = 303714];
    Float32 lonc[nele = 569405];
    Float32 latc[nele = 569405];
    Float32 siglay[siglay = 1][node = 303714];
    Float32 h[node = 303714];
    Float32 nv[three = 3][nele = 569405];
    Float64 time[time = 7];
    Float32 zeta[time = 7][node = 303714];
    Float32 u[time = 7][siglay = 1][nele = 569405];
    Float32 v[time = 7][siglay = 1][nele = 569405];
    Float32 tauc[time = 7][nele = 569405];
    Float32 temp[time = 7][siglay = 1][node = 303714];
    Float32 salinity[time = 7][siglay = 1][node = 303714];
    Float32 short_wave[time = 7][node = 303714];
    Float32 net_heat_flux[time = 7][node = 303714];
    Float32 sensible_heat_flux[time = 7][node = 303714];
    Float32 latent_heat_flux[time = 7][node = 303714];
    Float32 atmos_press[time = 7][node = 303714];
    Float32 uwind_speed[time = 7][nele = 569405];
    Float32 vwind_speed[time = 7][nele = 569405];
    Int32 wet_nodes[time = 7][node = 303714];
    Int32 wet_cells[time = 7][nele = 569405];
} NOAA%2fNGOFS2%2fMODELS%2f201210%2fnos%2engofs%2efields%2enowcast%2e\
20121001%2et03z%2enc;

Data:
*/

static NGOFS2Data* read_NGOFS2_header( const Arguments* arguments, FILE* file ) {
  const char* const variable = arguments ? arguments->variable : "";
  const int is_current = strcmp( variable, "current" ) == 0;
  const int has_depth =
    is_current ||
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0;
  const int is_nodal =
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0 ||
    strcmp( variable, "radiation" ) == 0 ||
    strcmp( variable, "net_heat_flux" ) == 0 ||
    strcmp( variable, "sensible_heat_flux" ) == 0 ||
    strcmp( variable, "latent_heat_flux" ) == 0 ||
    strcmp( variable, "air_pressure" ) == 0;
  NGOFS2Data* result = 0;
  int timesteps = 0;
  int points = 0;
  const int nodes = is_current ? NGOFS2_nodes : 0;
  int ok = 0;
  Line line = "";
  memset( line, 0, sizeof line );
  assert( arguments ); assert( file );

  ok = read_matched_line( file, "Dataset {\n" );

  /* Read lon and its dimensions: */

  if ( is_nodal ) {
    ok = ok &&
      read_dimension1_line( file, "Float32 lon[node = %d];\n", &points );
  } else {
    ok = ok &&
      read_dimension1_line( file, "Float32 lonc[nele = %d];\n", &points );
  }

  /* Read lat: */

  if ( is_nodal ) {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 lat[node = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  } else {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 latc[nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  if ( has_depth ) { /* Read sigma, depth: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 siglay[siglay = 1][node = %d];\n", count );
    ok = ok && read_matched_line( file, line );

    snprintf( line, sizeof line / sizeof *line,
              "Float32 h[node = %d];\n", count );
    ok = ok && read_matched_line( file, line );
  }

  if ( is_current ) { /* Read nv: */
    snprintf( line, sizeof line / sizeof *line,
              "Int32 nv[three = 3][nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  /* Read time: */

  ok = ok && read_dimension1_line( file, "Float64 time[time = %d];\n",
                                   &timesteps );

  if ( has_depth ) { /* Read zeta: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 zeta[time = %d][node = %d];\n", timesteps, count );
    ok = ok && read_matched_line( file, line );
  }

  /* Read main variable: */

  if ( ok ) {

    if ( ! strcmp( variable, "wind" ) ) { /* Read uwind_speed, vwind_speed: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 uwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 vwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "current" ) ) { /* Read u, v: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 u[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 v[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "water_temperature" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 temp[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "salinity" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 salinity[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "seabed_stress" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 tauc[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "radiation" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 short_wave[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "net_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 net_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "sensible_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 sensible_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "latent_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 latent_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "air_pressure" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 atmos_press[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else {
      fprintf( stderr, "Invalid variable '%s'.\n", variable );
      ok = 0;
    }

    /* Read mask: */

    if ( is_nodal ) {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_nodes[time = %d][node = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    } else {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_cells[time = %d][nele = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    }

    /* Read end of header: */

    DEBUG( fprintf( stderr, "before end of header: ok = %d\n", ok ); )

    if ( ok ) {
      ok = read_ignored_line( file ) &&
           read_ignored_line( file ) &&
           read_ignored_line( file );
    }

    DEBUG( fprintf( stderr, "after end of header: ok = %d\n", ok ); )
  }

  if ( ! ok ) {
    fprintf( stderr, "\nInvalid input DODS header.\n" );
  } else {
    result = NEW( NGOFS2Data, 1 );

    if ( result ) {
      result->timesteps = timesteps;
      result->points = points;
      result->nodes = nodes;
      DEBUG( fprintf( stderr, "timesteps = %d, points = %d, nodes = %d\n",
                      result->timesteps, result->points, result->nodes ); )
    }
  }

  DEBUG( fprintf( stderr, "read_NGOFS2_header result = %p\n", result ); )
  return result;
}



/* Read NGOFS2 data: */

static int read_NGOFS2_data( const Arguments* arguments, FILE* file,
                             NGOFS2Data* data ) {
  int ok = 1;
  assert( arguments ); assert( data ); assert( file );
  assert( data->timesteps > 0 );
  assert( data->points > 0 );
  DEBUG( fprintf( stderr,
                  "read_NGOFS2_data: timesteps = %d, points = %d\n",
                  data->timesteps, data->points ); )
  {
    const char* const variable = arguments->variable;
    const int is_current = strcmp( variable, "current" ) == 0;
    const int has_depth =
      is_current ||
      strcmp( variable, "water_temperature" ) == 0 ||
      strcmp( variable, "salinity" ) == 0;
    const int timesteps = data->timesteps;
    const int points = data->points;
    const int nodes = data->nodes;
    const int timesteps_points = timesteps * points;

    DEBUG( fprintf( stderr, "has_depth = %d\n", has_depth ); )

    /* Read and convert longitudes: */

    data->longitude =
      read_new_float_data( file, 1.0, -360.0, -180.0, 180.0, points );
    ok = data->longitude != 0;

    /* Read latitudes: */

    if ( ok ) {
      data->latitude =
        read_new_float_data( file, 1.0, 0.0, -90.0, 90.0, points );
      ok = data->latitude != 0;
    }

    /* If data has depth then read siglay and h: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->s = read_new_float_data( file, 1.0, 0.0, -1.0, 0.0, count );
      ok = data->s != 0;

      if ( ok ) {
        data->h = read_new_float_data( file, 1.0, 0.0, 0.0, 1.1e4, count );
        ok = data->h != 0;
      }
    }

    if ( ok && is_current ) { /* Read nv: */
      data->node =
        read_new_int_data( file, 1.0, 0.0, 0.0, NGOFS_nodes - 1, 3 * points );
      ok = data->node != 0;
    }

    /* Read time: */

    if ( ok ) {
      const int yyyymmddhh = base_yyyymmddhh( arguments );
      const double start_seconds =
        seconds_difference( arguments->yyyymmddhh, yyyymmddhh );
      const double end_seconds =
        start_seconds + arguments->hours * SECONDS_PER_HOUR - 1.0;
      const int scale = timestep_scale( arguments->source );
      data->seconds =
        read_new_double_data( file, scale, SECONDS_PER_HALF_HOUR,
                              start_seconds, end_seconds, timesteps );
      ok = data->seconds != 0;
    }

    /* If data has depth then read zeta: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->msl =
        read_new_float_data( file, 1.0, 0.0, -50.0, 50.0, timesteps * count );
      ok = data->msl != 0;
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "wind" ) ) {
        data->wind_u =
          read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                               timesteps_points );
        ok = data->wind_u != 0;

        if ( ok ) {
          data->wind_v =
            read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                                 timesteps_points );
          ok = data->wind_v != 0;
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        data->current_u =
          read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                               timesteps_points );
        ok = data->current_u != 0;

        if ( ok ) {
          data->current_v =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_points);
          ok = data->current_v != 0;
        }

      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        data->water_temperature =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->water_temperature != 0;
      } else if ( ! strcmp( variable, "salinity" ) ) {
        data->salinity =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->salinity != 0;
      } else if ( ! strcmp( variable, "seabed_stress" ) ) {
        data->seabed_stress =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1.0,
                               timesteps_points );
        ok = data->seabed_stress != 0;
      } else if ( ! strcmp( variable, "radiation" ) ) {
        data->radiation =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1500.0,
                               timesteps_points );
        ok = data->radiation != 0;
      } else if ( ! strcmp( variable, "net_heat_flux" ) ) {
        data->net_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->net_heat_flux != 0;
      } else if ( ! strcmp( variable, "sensible_heat_flux" ) ) {
        data->sensible_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->sensible_heat_flux != 0;
      } else if ( ! strcmp( variable, "latent_heat_flux" ) ) {
        data->latent_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->latent_heat_flux != 0;
      } else {
        assert( ! strcmp( variable, "air_pressure" ) );
        data->air_pressure =
          read_new_float_data( file, 0.01, 0.0, 500.0, 1500.0,
                               timesteps_points );
        ok = data->air_pressure != 0;
      }
    }

    /* Read mask: */

    if ( ok ) {
      data->mask =
        read_new_int_data( file, 1.0, 0.0, 0.0, 1.0, timesteps_points );
      ok = data->mask != 0;
    }
  }

  DEBUG( fprintf( stderr, "read_NGOFS2_data result = %d\n", ok ); )
  return ok;
}



/* Convert NGOFS data: */

static OutputData* convert_NGOFS2( const Arguments* arguments,
                                   NGOFS2Data* data ) {
  OutputData* result = 0;
  const int subset_points =
    spatially_subset( arguments->longitude_minimum,
                      arguments->longitude_maximum,
                      arguments->latitude_minimum,
                      arguments->latitude_maximum,
                      data->points,
                      sizeof data->mask[ 0 ],
                      data->longitude,
                      data->latitude,
                      data->mask );
  DEBUG( fprintf( stderr, "subset points = %d\n", subset_points ); )

  if ( subset_points ) {
    result = NEW( OutputData, 1 );

    if ( result ) {
      const size_t length = sizeof (Name) / sizeof (char);
      const int timesteps = data->timesteps;
      int timestep = 0;
      result->points = subset_points;

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        const double seconds = data->seconds[ timestep ];

        if ( seconds != MISSING ) {
          result->timesteps += 1;
        }
      }

      assert( result->timesteps > 0 );

      if ( data->wind_u ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "wind_u", length );
        strncpy( result->names[ 1 ], "wind_v", length );
        strncpy( result->units[ 0 ], "m/s", length );
        strncpy( result->units[ 1 ], "m/s", length );
      } else if ( data->current_u ) {
        result->variables = 3;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "current_u", length );
        strncpy( result->names[ 2 ], "current_v", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "m/s", length );
        strncpy( result->units[ 2 ], "m/s", length );
      } else if ( data->water_temperature ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "water_temperature", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "C", length );
      } else if ( data->salinity ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "salinity", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "PSU", length );
      } else if ( data->seabed_stress ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "seabed_stress", length );
        strncpy( result->units[ 0 ], "m2/s2", length );
      } else if ( data->radiation ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "radiation", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else if ( data->air_pressure ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "air_pressure", length );
        strncpy( result->units[ 0 ], "hPa", length );
      } else if ( data->net_heat_flux ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "net_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else if ( data->sensible_heat_flux ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "sensible_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else {
        assert( data->latent_heat_flux );
        result->variables = 1;
        strncpy( result->names[ 0 ], "latent_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      }

      result->yyyymmddhh = NEW( int, timesteps );
      result->longitudes = result->yyyymmddhh ? NEW( float, subset_points) : 0;
      result->latitudes  = result->longitudes ? NEW( float, subset_points) : 0;
      result->data = result->latitudes ?
        NEW( float, result->variables * timesteps * subset_points ) : 0;

      if ( result->data ) {
        float* output_data = result->data;
        const int points = data->points;
        const int nodes = data->nodes;
        const float* const mask = data->mask;
        const float* const longitudes = data->longitude;
        const float* const latitudes = data->latitude;
        const float* input_data =
          data->wind_u ? data->wind_u
          : data->current_u ? data->current_u
          : data->water_temperature ? data->water_temperature
          : data->salinity ? data->salinity
          : data->seabed_stress ? data->seabed_stress
          : data->radiation ? data->radiation
          : data->net_heat_flux ? data->net_heat_flux
          : data->sensible_heat_flux ? data->sensible_heat_flux
          : data->latent_heat_flux ? data->latent_heat_flux
          : data->air_pressure ? data->air_pressure
          : 0;
        const float* input_data_2 =
          data->wind_v ? data->wind_v : data->current_v;
        const int yyyymmddhh0 = base_yyyymmddhh( arguments );
        int output_timesteps = 0;

        /* Copy timestamps and count number of output timesteps: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            const int yyyymmddhh =
              convert_seconds( yyyymmddhh0, (int) seconds );
            assert( output_timesteps < 1 ||
                    yyyymmddhh > result->yyyymmddhh[ output_timesteps - 1 ] );
            result->yyyymmddhh[ output_timesteps ] = yyyymmddhh;
            ++output_timesteps;
          }
        }

        /* Copy longitudes: */

        copy_unmasked_float_values( points, mask, longitudes,
                                    result->longitudes );
        assert( IN_RANGE( result->longitudes[ 0 ], -180.0, 180.0 ) );
        assert( IN_RANGE( result->longitudes[ subset_points - 1],
                          -180.0, 180.0 ) );

        /* Copy latitudes: */

        copy_unmasked_float_values( points, mask, latitudes,
                                    result->latitudes );
        assert( IN_RANGE( result->latitudes[ 0 ], -90.0, 90.0 ) );
        assert( IN_RANGE( result->latitudes[ subset_points - 1 ],
                          -90.0, 90.0 ) );

        /* Copy depth if present: */

        if ( data->s ) {
          assert( data->h ); assert( data->msl );

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {

              if ( nodes ) {
                copy_float_depth3( points, data->node,
                                   mask, data->s, data->h,
                                   data->msl + timestep * nodes,
                                   output_data );
              } else {
                copy_float_depth( points, mask, 0.0, data->s, data->h,
                                  data->msl + timestep * points,
                                  output_data );
              }
              assert( IN_RANGE( output_data[ 0 ], -3e5, 0.0 ) );
              assert( IN_RANGE( output_data[ subset_points - 1 ], -3e5, 0.0 ));
              output_data += subset_points;
            }
          }
        }

        /* Copy data variable: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            copy_unmasked_float_values( points, mask, input_data, output_data);
            output_data += subset_points;
            input_data  += points;
          }
        }

        /* Copy additional component of data variable: */

        if ( input_data_2 ) {

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_unmasked_float_values( points, mask, input_data_2,
                                          output_data );
              output_data += subset_points;
              input_data_2 += points;
            }
          }
        }
      }
    }
  }

  return result;
}



/* Read and convert SSCOFS data: */

static OutputData* read_SSCOFS( const Arguments* arguments, FILE* file ) {
  OutputData* result = 0;
  SSCOFSData* data = read_SSCOFS_header( arguments, file );

  if ( data ) {

    if ( read_SSCOFS_data( arguments, file, data ) ) {
      result = convert_SSCOFS( arguments, data );
    }

    deallocate_sscofs_data( data ), data = 0;
  }

  return result;
}

/*
 Input data from SSCOFS DODS files whose header looks like a subset of:

 Dataset {
    Float32 lon[node = 239734];
    Float32 lat[node = 239734];
    Float32 lonc[nele = 433410];
    Float32 latc[nele = 433410];
    Float32 siglay[siglay = 10][node = 239734];
    Float32 h[node = 239734];
    Int32 nv[three = 3][nele = 433410];
    Float64 time[time = 1];
    Float32 zeta[time = 1][node = 239734];
    Float32 u[time = 1][siglay = 10][nele = 433410];
    Float32 v[time = 1][siglay = 10][nele = 433410];
    Float32 tauc[time = 1][nele = 433410];
    Float32 temp[time = 1][siglay = 10][node = 239734];
    Float32 salinity[time = 1][siglay = 10][node = 239734];
    Float32 short_wave[time = 1][node = 239734];
    Float32 net_heat_flux[time = 1][node = 239734];
    Float32 sensible_heat_flux[time = 1][node = 239734];
    Float32 latent_heat_flux[time = 1][node = 239734];
    Float32 long_wave[time = 1][node = 239734];
    Float32 uwind_speed[time = 1][nele = 433410];
    Float32 vwind_speed[time = 1][nele = 433410];
    Float32 atmos_press[time = 1][node = 239734];
    Int32 wet_nodes[time = 1][node = 239734];
    Int32 wet_cells[time = 1][nele = 433410];
} NOAA/SSCOFS/MODELS/2024/09/16/sscofs.t21z.20240916.fields.n006.nc;

Data:
*/

static SSCOFSData* read_SSCOFS_header( const Arguments* arguments, FILE* file ) {
  const char* const variable = arguments ? arguments->variable : "";
  const int is_current = strcmp( variable, "current" ) == 0;
  const int has_depth =
    is_current ||
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0;
  const int is_nodal =
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0 ||
    strcmp( variable, "radiation" ) == 0 ||
    strcmp( variable, "net_heat_flux" ) == 0 ||
    strcmp( variable, "sensible_heat_flux" ) == 0 ||
    strcmp( variable, "latent_heat_flux" ) == 0 ||
    strcmp( variable, "air_pressure" ) == 0;
  SSCOFSData* result = 0;
  int timesteps = 0;
  int points = 0;
  const int nodes = is_current ? SSCOFS_nodes : 0;
  int ok = 0;
  Line line = "";
  memset( line, 0, sizeof line );
  assert( arguments ); assert( file );

  ok = read_matched_line( file, "Dataset {\n" );

  /* Read lon and its dimensions: */

  if ( is_nodal ) {
    ok = ok &&
      read_dimension1_line( file, "Float32 lon[node = %d];\n", &points );
  } else {
    ok = ok &&
      read_dimension1_line( file, "Float32 lonc[nele = %d];\n", &points );
  }

  /* Read lat: */

  if ( is_nodal ) {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 lat[node = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  } else {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 latc[nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  if ( has_depth ) { /* Read sigma, depth: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 siglay[siglay = 1][node = %d];\n", count );
    ok = ok && read_matched_line( file, line );

    snprintf( line, sizeof line / sizeof *line,
              "Float32 h[node = %d];\n", count );
    ok = ok && read_matched_line( file, line );
  }

  if ( is_current ) { /* Read nv: */
    snprintf( line, sizeof line / sizeof *line,
              "Int32 nv[three = 3][nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  /* Read time: */

  ok = ok && read_dimension1_line( file, "Float64 time[time = %d];\n",
                                   &timesteps );

  if ( has_depth ) { /* Read zeta: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 zeta[time = %d][node = %d];\n", timesteps, count );
    ok = ok && read_matched_line( file, line );
  }

  /* Read main variable: */

  if ( ok ) {

    if ( ! strcmp( variable, "wind" ) ) { /* Read uwind_speed, vwind_speed: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 uwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 vwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "current" ) ) { /* Read u, v: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 u[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 v[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "water_temperature" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 temp[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "salinity" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 salinity[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "seabed_stress" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 tauc[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "radiation" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 short_wave[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "net_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 net_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "sensible_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 sensible_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "latent_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 latent_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "air_pressure" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 atmos_press[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else {
      fprintf( stderr, "Invalid variable '%s'.\n", variable );
      ok = 0;
    }

    /* Read mask: */

    if ( is_nodal ) {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_nodes[time = %d][node = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    } else {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_cells[time = %d][nele = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    }

    /* Read end of header: */

    DEBUG( fprintf( stderr, "before end of header: ok = %d\n", ok ); )

    if ( ok ) {
      ok = read_ignored_line( file ) &&
           read_ignored_line( file ) &&
           read_ignored_line( file );
    }

    DEBUG( fprintf( stderr, "after end of header: ok = %d\n", ok ); )
  }

  if ( ! ok ) {
    fprintf( stderr, "\nInvalid input DODS header.\n" );
  } else {
    result = NEW( SSCOFSData, 1 );

    if ( result ) {
      result->timesteps = timesteps;
      result->points = points;
      result->nodes = nodes;
      DEBUG( fprintf( stderr, "timesteps = %d, points = %d, nodes = %d\n",
                      result->timesteps, result->points, result->nodes ); )
    }
  }

  DEBUG( fprintf( stderr, "read_SSCOFS_header result = %p\n", result ); )
  return result;
}



/* Read SSCOFS data: */

static int read_SSCOFS_data( const Arguments* arguments, FILE* file,
                             SSCOFSData* data ) {
  int ok = 1;
  assert( arguments ); assert( data ); assert( file );
  assert( data->timesteps > 0 );
  assert( data->points > 0 );
  DEBUG( fprintf( stderr,
                  "read_SSCOFS_data: timesteps = %d, points = %d\n",
                  data->timesteps, data->points ); )
  {
    const char* const variable = arguments->variable;
    const int is_current = strcmp( variable, "current" ) == 0;
    const int has_depth =
      is_current ||
      strcmp( variable, "water_temperature" ) == 0 ||
      strcmp( variable, "salinity" ) == 0;
    const int timesteps = data->timesteps;
    const int points = data->points;
    const int nodes = data->nodes;
    const int timesteps_points = timesteps * points;

    DEBUG( fprintf( stderr, "has_depth = %d\n", has_depth ); )

    /* Read and convert longitudes: */

    data->longitude =
      read_new_float_data( file, 1.0, -360.0, -180.0, 180.0, points );
    ok = data->longitude != 0;

    /* Read latitudes: */

    if ( ok ) {
      data->latitude =
        read_new_float_data( file, 1.0, 0.0, -90.0, 90.0, points );
      ok = data->latitude != 0;
    }

    /* If data has depth then read siglay and h: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->s = read_new_float_data( file, 1.0, 0.0, -1.0, 0.0, count );
      ok = data->s != 0;

      if ( ok ) {
        data->h = read_new_float_data( file, 1.0, 0.0, 0.0, 1.1e4, count );
        ok = data->h != 0;
      }
    }

    if ( ok && is_current ) { /* Read nv: */
      data->node =
        read_new_int_data( file, 1.0, 0.0, 0.0, NGOFS_nodes - 1, 3 * points );
      ok = data->node != 0;
    }

    /* Read time: */

    if ( ok ) {
      const int yyyymmddhh = base_yyyymmddhh( arguments );
      const double start_seconds =
        seconds_difference( arguments->yyyymmddhh, yyyymmddhh );
      const double end_seconds =
        start_seconds + arguments->hours * SECONDS_PER_HOUR - 1.0;
      const int scale = timestep_scale( arguments->source );
      data->seconds =
        read_new_double_data( file, scale, SECONDS_PER_HALF_HOUR,
                              start_seconds, end_seconds, timesteps );
      ok = data->seconds != 0;
    }

    /* If data has depth then read zeta: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->msl =
        read_new_float_data( file, 1.0, 0.0, -50.0, 50.0, timesteps * count );
      ok = data->msl != 0;
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "wind" ) ) {
        data->wind_u =
          read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                               timesteps_points );
        ok = data->wind_u != 0;

        if ( ok ) {
          data->wind_v =
            read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                                 timesteps_points );
          ok = data->wind_v != 0;
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        data->current_u =
          read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                               timesteps_points );
        ok = data->current_u != 0;

        if ( ok ) {
          data->current_v =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_points);
          ok = data->current_v != 0;
        }

      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        data->water_temperature =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->water_temperature != 0;
      } else if ( ! strcmp( variable, "salinity" ) ) {
        data->salinity =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->salinity != 0;
      } else if ( ! strcmp( variable, "seabed_stress" ) ) {
        data->seabed_stress =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1.0,
                               timesteps_points );
        ok = data->seabed_stress != 0;
      } else if ( ! strcmp( variable, "radiation" ) ) {
        data->radiation =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1500.0,
                               timesteps_points );
        ok = data->radiation != 0;
      } else if ( ! strcmp( variable, "net_heat_flux" ) ) {
        data->net_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->net_heat_flux != 0;
      } else if ( ! strcmp( variable, "sensible_heat_flux" ) ) {
        data->sensible_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->sensible_heat_flux != 0;
      } else if ( ! strcmp( variable, "latent_heat_flux" ) ) {
        data->latent_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->latent_heat_flux != 0;
      } else {
        assert( ! strcmp( variable, "air_pressure" ) );
        data->air_pressure =
          read_new_float_data( file, 0.01, 0.0, 500.0, 1500.0,
                               timesteps_points );
        ok = data->air_pressure != 0;
      }
    }

    /* Read mask: */

    if ( ok ) {
      data->mask =
        read_new_int_data( file, 1.0, 0.0, 0.0, 1.0, timesteps_points );
      ok = data->mask != 0;
    }
  }

  DEBUG( fprintf( stderr, "read_SSCOFS_data result = %d\n", ok ); )
  return ok;
}



/* Convert SSCOFS data: */

static OutputData* convert_SSCOFS( const Arguments* arguments,
                                   SSCOFSData* data ) {
  OutputData* result = 0;
  const int subset_points =
    spatially_subset( arguments->longitude_minimum,
                      arguments->longitude_maximum,
                      arguments->latitude_minimum,
                      arguments->latitude_maximum,
                      data->points,
                      sizeof data->mask[ 0 ],
                      data->longitude,
                      data->latitude,
                      data->mask );
  DEBUG( fprintf( stderr, "subset points = %d\n", subset_points ); )

  if ( subset_points ) {
    result = NEW( OutputData, 1 );

    if ( result ) {
      const size_t length = sizeof (Name) / sizeof (char);
      const int timesteps = data->timesteps;
      int timestep = 0;
      result->points = subset_points;

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        const double seconds = data->seconds[ timestep ];

        if ( seconds != MISSING ) {
          result->timesteps += 1;
        }
      }

      assert( result->timesteps > 0 );

      if ( data->wind_u ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "wind_u", length );
        strncpy( result->names[ 1 ], "wind_v", length );
        strncpy( result->units[ 0 ], "m/s", length );
        strncpy( result->units[ 1 ], "m/s", length );
      } else if ( data->current_u ) {
        result->variables = 3;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "current_u", length );
        strncpy( result->names[ 2 ], "current_v", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "m/s", length );
        strncpy( result->units[ 2 ], "m/s", length );
      } else if ( data->water_temperature ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "water_temperature", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "C", length );
      } else if ( data->salinity ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "salinity", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "PSU", length );
      } else if ( data->seabed_stress ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "seabed_stress", length );
        strncpy( result->units[ 0 ], "m2/s2", length );
      } else if ( data->radiation ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "radiation", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else if ( data->air_pressure ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "air_pressure", length );
        strncpy( result->units[ 0 ], "hPa", length );
      } else if ( data->net_heat_flux ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "net_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else if ( data->sensible_heat_flux ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "sensible_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else {
        assert( data->latent_heat_flux );
        result->variables = 1;
        strncpy( result->names[ 0 ], "latent_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      }

      result->yyyymmddhh = NEW( int, timesteps );
      result->longitudes = result->yyyymmddhh ? NEW( float, subset_points) : 0;
      result->latitudes  = result->longitudes ? NEW( float, subset_points) : 0;
      result->data = result->latitudes ?
        NEW( float, result->variables * timesteps * subset_points ) : 0;

      if ( result->data ) {
        float* output_data = result->data;
        const int points = data->points;
        const int nodes = data->nodes;
        const float* const mask = data->mask;
        const float* const longitudes = data->longitude;
        const float* const latitudes = data->latitude;
        const float* input_data =
          data->wind_u ? data->wind_u
          : data->current_u ? data->current_u
          : data->water_temperature ? data->water_temperature
          : data->salinity ? data->salinity
          : data->seabed_stress ? data->seabed_stress
          : data->radiation ? data->radiation
          : data->net_heat_flux ? data->net_heat_flux
          : data->sensible_heat_flux ? data->sensible_heat_flux
          : data->latent_heat_flux ? data->latent_heat_flux
          : data->air_pressure ? data->air_pressure
          : 0;
        const float* input_data_2 =
          data->wind_v ? data->wind_v : data->current_v;
        const int yyyymmddhh0 = base_yyyymmddhh( arguments );
        int output_timesteps = 0;

        /* Copy timestamps and count number of output timesteps: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            const int yyyymmddhh =
              convert_seconds( yyyymmddhh0, (int) seconds );
            assert( output_timesteps < 1 ||
                    yyyymmddhh > result->yyyymmddhh[ output_timesteps - 1 ] );
            result->yyyymmddhh[ output_timesteps ] = yyyymmddhh;
            ++output_timesteps;
          }
        }

        /* Copy longitudes: */

        copy_unmasked_float_values( points, mask, longitudes,
                                    result->longitudes );
        assert( IN_RANGE( result->longitudes[ 0 ], -180.0, 180.0 ) );
        assert( IN_RANGE( result->longitudes[ subset_points - 1],
                          -180.0, 180.0 ) );

        /* Copy latitudes: */

        copy_unmasked_float_values( points, mask, latitudes,
                                    result->latitudes );
        assert( IN_RANGE( result->latitudes[ 0 ], -90.0, 90.0 ) );
        assert( IN_RANGE( result->latitudes[ subset_points - 1 ],
                          -90.0, 90.0 ) );

        /* Copy depth if present: */

        if ( data->s ) {
          assert( data->h ); assert( data->msl );

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {

              if ( nodes ) {
                copy_float_depth3( points, data->node,
                                   mask, data->s, data->h,
                                   data->msl + timestep * nodes,
                                   output_data );
              } else {
                copy_float_depth( points, mask, 0.0, data->s, data->h,
                                  data->msl + timestep * points,
                                  output_data );
              }
              assert( IN_RANGE( output_data[ 0 ], -3e5, 0.0 ) );
              assert( IN_RANGE( output_data[ subset_points - 1 ], -3e5, 0.0 ));
              output_data += subset_points;
            }
          }
        }

        /* Copy data variable: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            copy_unmasked_float_values( points, mask, input_data, output_data);
            output_data += subset_points;
            input_data  += points;
          }
        }

        /* Copy additional component of data variable: */

        if ( input_data_2 ) {

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_unmasked_float_values( points, mask, input_data_2,
                                          output_data );
              output_data += subset_points;
              input_data_2 += points;
            }
          }
        }
      }
    }
  }

  return result;
}



/* Read and convert SFBOFS data: */

static OutputData* read_SFBOFS( const Arguments* arguments, FILE* file ) {
  OutputData* result = 0;
  SFBOFSData* data = read_SFBOFS_header( arguments, file );

  if ( data ) {

    if ( read_SFBOFS_data( arguments, file, data ) ) {
      result = convert_SFBOFS( arguments, data );
    }

    deallocate_sfbofs_data( data ), data = 0;
  }

  return result;
}


/* Read/check DODS header for SFBOFS file that looks like a subset of:
Dataset {
    Float32 lon[node = 54120];
    Float32 lat[node = 54120];
    Float32 lonc[nele = 102264];
    Float32 latc[nele = 102264];
    Float32 siglay[siglay = 20][node = 54120];
    Float32 h[node = 54120];
    Float64 time[time = 1];
    Float32 zeta[time = 1][node = 54120];
    Float32 u[time = 1][siglay = 20][nele = 102264];
    Float32 v[time = 1][siglay = 20][nele = 102264];
    Float32 tauc[time = 1][nele = 102264];
    Float32 temp[time = 1][siglay = 20][node = 54120];
    Float32 salinity[time = 1][siglay = 20][node = 54120];
    Float32 short_wave[time = 1][node = 54120];
    Float32 net_heat_flux[time = 1][node = 54120];
    Float32 sensible_heat_flux[time = 1][node = 54120];
    Float32 latent_heat_flux[time = 1][node = 54120];
    Float32 uwind_speed[time = 1][nele = 102264];
    Float32 vwind_speed[time = 1][nele = 102264];
    Float32 atmos_press[time = 1][node = 54120];
    Int32 wet_nodes[time = 1][node = 54120];
    Int32 wet_cells[time = 1][nele = 102264];
} NOAA/SFBOFS/MODELS/201507/nos.sfbofs.fields.n006.20150727.t15z.nc;

Data:
*/

static SFBOFSData* read_SFBOFS_header( const Arguments* arguments,
                                       FILE* file ) {
  const char* const variable = arguments ? arguments->variable : "";
  const int is_current = strcmp( variable, "current" ) == 0;
  const int has_depth =
    is_current ||
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0;
  const int is_nodal =
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0 ||
    strcmp( variable, "radiation" ) == 0 ||
    strcmp( variable, "net_heat_flux" ) == 0 ||
    strcmp( variable, "sensible_heat_flux" ) == 0 ||
    strcmp( variable, "latent_heat_flux" ) == 0 ||
    strcmp( variable, "air_pressure" ) == 0;
  SFBOFSData* result = 0;
  int timesteps = 0;
  int points = 0;
  const int nodes =
    ! is_current ? 0
    : ! strcmp( arguments->source, "sfbofs" ) ? SFBOFS_nodes
    : 0;
  int ok = 0;
  Line line = "";
  memset( line, 0, sizeof line );
  assert( arguments ); assert( file );

  ok = read_matched_line( file, "Dataset {\n" );

  /* Read lon and its dimensions: */

  if ( is_nodal ) {
    ok = ok &&
      read_dimension1_line( file, "Float32 lon[node = %d];\n", &points );
  } else {
    ok = ok &&
      read_dimension1_line( file, "Float32 lonc[nele = %d];\n", &points );
  }

  /* Read lat: */

  if ( is_nodal ) {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 lat[node = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  } else {
    snprintf( line, sizeof line / sizeof *line,
              "Float32 latc[nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  if ( has_depth ) { /* Read sigma, depth: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 siglay[siglay = 1][node = %d];\n", count );
    ok = ok && read_matched_line( file, line );

    snprintf( line, sizeof line / sizeof *line,
              "Float32 h[node = %d];\n", count );
    ok = ok && read_matched_line( file, line );
  }

  if ( is_current ) { /* Read nv: */
    snprintf( line, sizeof line / sizeof *line,
              "Int32 nv[three = 3][nele = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  /* Read time: */

  ok = ok && read_dimension1_line( file, "Float64 time[time = %d];\n",
                                   &timesteps );

  if ( has_depth ) { /* Read zeta: */
    const int count = is_current ? nodes : points;
    snprintf( line, sizeof line / sizeof *line,
              "Float32 zeta[time = %d][node = %d];\n", timesteps, count );
    ok = ok && read_matched_line( file, line );
  }

  /* Read main variable: */

  if ( ok ) {

    if ( ! strcmp( variable, "wind" ) ) { /* Read uwind_speed, vwind_speed: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 uwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 vwind_speed[time = %d][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "current" ) ) { /* Read u, v: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 u[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 v[time = %d][siglay = 1][nele = %d];\n",
                timesteps, points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "water_temperature" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 temp[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "salinity" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 salinity[time = %d][siglay = 1][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "seabed_stress" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 tauc[time = %d][nele = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "radiation" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 short_wave[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "net_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 net_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "sensible_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 sensible_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "latent_heat_flux" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 latent_heat_flux[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "air_pressure" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 atmos_press[time = %d][node = %d];\n",
                timesteps, points );
      ok = read_matched_line( file, line );
    } else {
      fprintf( stderr, "Invalid variable '%s'.\n", variable );
      ok = 0;
    }

    /* Read mask: */

    if ( is_nodal ) {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_nodes[time = %d][node = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    } else {
      snprintf( line, sizeof line / sizeof *line,
                "Int32 wet_cells[time = %d][nele = %d];\n", timesteps, points);
      ok = ok && read_matched_line( file, line );
    }

    /* Read end of header: */

    DEBUG( fprintf( stderr, "before end of header: ok = %d\n", ok ); )

    if ( ok ) {
      ok = read_ignored_line( file ) &&
           read_ignored_line( file ) &&
           read_ignored_line( file );
    }

    DEBUG( fprintf( stderr, "after end of header: ok = %d\n", ok ); )
  }

  if ( ! ok ) {
    fprintf( stderr, "\nInvalid input DODS header.\n" );
  } else {
    result = NEW( SFBOFSData, 1 );

    if ( result ) {
      result->timesteps = timesteps;
      result->points = points;
      result->nodes = nodes;
      DEBUG( fprintf( stderr, "timesteps = %d, points = %d, nodes = %d\n",
                      result->timesteps, result->points, result->nodes ); )
    }
  }

  DEBUG( fprintf( stderr, "read_SFBOFS_header result = %p\n", result ); )
  return result;
}



/* Read SFBOFS data: */

static int read_SFBOFS_data( const Arguments* arguments, FILE* file,
                             SFBOFSData* data ) {
  int ok = 1;
  assert( arguments ); assert( data ); assert( file );
  assert( data->timesteps > 0 );
  assert( data->points > 0 );
  DEBUG( fprintf( stderr,
                  "read_SFBOFS_data: timesteps = %d, points = %d\n",
                  data->timesteps, data->points ); )
  {
    const char* const variable = arguments->variable;
    const int is_current = strcmp( variable, "current" ) == 0;
    const int has_depth =
      is_current ||
      strcmp( variable, "water_temperature" ) == 0 ||
      strcmp( variable, "salinity" ) == 0;
    const int timesteps = data->timesteps;
    const int points = data->points;
    const int nodes = data->nodes;
    const int timesteps_points = timesteps * points;

    DEBUG( fprintf( stderr, "has_depth = %d\n", has_depth ); )

    /* Read and convert longitudes: */

    data->longitude =
      read_new_float_data( file, 1.0, -360.0, -180.0, 180.0, points );
    ok = data->longitude != 0;

    /* Read latitudes: */

    if ( ok ) {
      data->latitude =
        read_new_float_data( file, 1.0, 0.0, -90.0, 90.0, points );
      ok = data->latitude != 0;
    }

    /* If data has depth then read siglay and h: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->s = read_new_float_data( file, 1.0, 0.0, -1.0, 0.0, count );
      ok = data->s != 0;

      if ( ok ) {
        data->h = read_new_float_data( file, 1.0, 0.0, 0.0, 1.1e4, count );
        ok = data->h != 0;
      }
    }

    if ( ok && is_current ) { /* Read nv: */
      data->node =
        read_new_int_data( file, 1.0, 0.0, 0.0, SFBOFS_nodes - 1, 3 * points );
      ok = data->node != 0;
    }

    /* Read time: */

    if ( ok ) {
      const int yyyymmddhh = base_yyyymmddhh( arguments );
      const double start_seconds =
        seconds_difference( arguments->yyyymmddhh, yyyymmddhh );
      const double end_seconds =
        start_seconds + arguments->hours * SECONDS_PER_HOUR - 1.0;
      const int scale = timestep_scale( arguments->source );
      data->seconds =
        read_new_double_data( file, scale, SECONDS_PER_HALF_HOUR,
                              start_seconds, end_seconds, timesteps );
      ok = data->seconds != 0;
    }

    /* If data has depth then read zeta: */

    if ( ok && has_depth ) {
      const int count = nodes ? nodes : points;
      data->msl =
        read_new_float_data( file, 1.0, 0.0, -50.0, 50.0, timesteps * count );
      ok = data->msl != 0;
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "wind" ) ) {
        data->wind_u =
          read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                               timesteps_points );
        ok = data->wind_u != 0;

        if ( ok ) {
          data->wind_v =
            read_new_float_data( file, 1.0, 0.0, -100.0, 100.0,
                                 timesteps_points );
          ok = data->wind_v != 0;
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        data->current_u =
          read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                               timesteps_points );
        ok = data->current_u != 0;

        if ( ok ) {
          data->current_v =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0,
                                 timesteps_points);
          ok = data->current_v != 0;
        }

      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        data->water_temperature =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->water_temperature != 0;
      } else if ( ! strcmp( variable, "salinity" ) ) {
        data->salinity =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0,
                               timesteps_points );
        ok = data->salinity != 0;
      } else if ( ! strcmp( variable, "seabed_stress" ) ) {
        data->seabed_stress =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1.0,
                               timesteps_points );
        ok = data->seabed_stress != 0;
      } else if ( ! strcmp( variable, "radiation" ) ) {
        data->radiation =
          read_new_float_data( file, 1.0, 0.0, 0.0, 1500.0,
                               timesteps_points );
        ok = data->radiation != 0;
      } else if ( ! strcmp( variable, "net_heat_flux" ) ) {
        data->net_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->net_heat_flux != 0;
      } else if ( ! strcmp( variable, "sensible_heat_flux" ) ) {
        data->sensible_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->sensible_heat_flux != 0;
      } else if ( ! strcmp( variable, "latent_heat_flux" ) ) {
        data->latent_heat_flux =
          read_new_float_data( file, 1.0, 0.0, -1500.0, 1500.0,
                               timesteps_points );
        ok = data->latent_heat_flux != 0;
      } else {
        const double Pa_to_hPa = 1e-2;
        assert( ! strcmp( variable, "air_pressure" ) );
        data->air_pressure =
          read_new_float_data( file, Pa_to_hPa, 0.0, 500.0, 1500.0,
                               timesteps_points );
        ok = data->air_pressure != 0;
      }
    }

    /* Read mask: */

    if ( ok ) {
      data->mask =
        read_new_int_data( file, 1.0, 0.0, 0.0, 1.0, timesteps_points );
      ok = data->mask != 0;
    }
  }

  DEBUG( fprintf( stderr, "read_SFBOFS_data result = %d\n", ok ); )
  return ok;
}



/* Convert SFBOFS data: */

static OutputData* convert_SFBOFS( const Arguments* arguments,
                                   SFBOFSData* data ) {
  OutputData* result = 0;
  const int subset_points =
    spatially_subset( arguments->longitude_minimum,
                      arguments->longitude_maximum,
                      arguments->latitude_minimum,
                      arguments->latitude_maximum,
                      data->points,
                      sizeof data->mask[ 0 ],
                      data->longitude,
                      data->latitude,
                      data->mask );
  DEBUG( fprintf( stderr, "subset points = %d\n", subset_points ); )

  if ( subset_points ) {
    result = NEW( OutputData, 1 );

    if ( result ) {
      const size_t length = sizeof (Name) / sizeof (char);
      const int timesteps = data->timesteps;
      int timestep = 0;
      result->points = subset_points;

      for ( timestep = 0; timestep < timesteps; ++timestep ) {
        const double seconds = data->seconds[ timestep ];

        if ( seconds != MISSING ) {
          result->timesteps += 1;
        }
      }

      assert( result->timesteps > 0 );

      if ( data->wind_u ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "wind_u", length );
        strncpy( result->names[ 1 ], "wind_v", length );
        strncpy( result->units[ 0 ], "m/s", length );
        strncpy( result->units[ 1 ], "m/s", length );
      } else if ( data->current_u ) {
        result->variables = 3;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "current_u", length );
        strncpy( result->names[ 2 ], "current_v", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "m/s", length );
        strncpy( result->units[ 2 ], "m/s", length );
      } else if ( data->water_temperature ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "water_temperature", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "C", length );
      } else if ( data->salinity ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "salinity", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "PSU", length );
      } else if ( data->seabed_stress ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "seabed_stress", length );
        strncpy( result->units[ 0 ], "m2/s2", length );
      } else if ( data->radiation ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "radiation", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else if ( data->net_heat_flux ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "net_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else if ( data->sensible_heat_flux ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "sensible_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else if ( data->latent_heat_flux ) {
        result->variables = 1;
        strncpy( result->names[ 0 ], "latent_heat_flux", length );
        strncpy( result->units[ 0 ], "W/m2", length );
      } else {
        assert( data->air_pressure );
        result->variables = 1;
        strncpy( result->names[ 0 ], "air_pressure", length );
        strncpy( result->units[ 0 ], "hPa", length );
      }

      result->yyyymmddhh = NEW( int, timesteps );
      result->longitudes = result->yyyymmddhh ? NEW( float, subset_points) : 0;
      result->latitudes  = result->longitudes ? NEW( float, subset_points) : 0;
      result->data = result->latitudes ?
        NEW( float, result->variables * timesteps * subset_points ) : 0;

      if ( result->data ) {
        float* output_data = result->data;
        const int points = data->points;
        const int nodes = data->nodes;
        const float* const mask = data->mask;
        const float* const longitudes = data->longitude;
        const float* const latitudes = data->latitude;
        const float* input_data =
          data->wind_u ? data->wind_u
          : data->current_u ? data->current_u
          : data->water_temperature ? data->water_temperature
          : data->salinity ? data->salinity
          : data->seabed_stress ? data->seabed_stress
          : data->radiation ? data->radiation
          : data->net_heat_flux ? data->net_heat_flux
          : data->sensible_heat_flux ? data->sensible_heat_flux
          : data->latent_heat_flux ? data->latent_heat_flux
          : data->air_pressure ? data->air_pressure
          : 0;
        const float* input_data_2 =
          data->wind_v ? data->wind_v : data->current_v;
        const int yyyymmddhh0 = base_yyyymmddhh( arguments );
        int output_timesteps = 0;

        /* Copy timestamps and count number of output timesteps: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            const int yyyymmddhh =
              convert_seconds( yyyymmddhh0, (int) seconds );
            assert( output_timesteps < 1 ||
                    yyyymmddhh > result->yyyymmddhh[ output_timesteps - 1 ] );
            result->yyyymmddhh[ output_timesteps ] = yyyymmddhh;
            ++output_timesteps;
          }
        }

        /* Copy longitudes: */

        copy_unmasked_float_values( points, mask, longitudes,
                                    result->longitudes );
        assert( IN_RANGE( result->longitudes[ 0 ], -180.0, 180.0 ) );
        assert( IN_RANGE( result->longitudes[ subset_points - 1],
                          -180.0, 180.0 ) );

        /* Copy latitudes: */

        copy_unmasked_float_values( points, mask, latitudes,
                                    result->latitudes );
        assert( IN_RANGE( result->latitudes[ 0 ], -90.0, 90.0 ) );
        assert( IN_RANGE( result->latitudes[ subset_points - 1 ],
                          -90.0, 90.0 ) );

        /* Copy depth if present: */

        if ( data->s ) {
          assert( data->h ); assert( data->msl );

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {

              if ( nodes ) {
                copy_float_depth3( points, data->node,
                                   mask, data->s, data->h,
                                   data->msl + timestep * nodes,
                                   output_data );
              } else {
                copy_float_depth( points, mask, 0.0, data->s, data->h,
                                  data->msl + timestep * points,
                                  output_data );
              }
              assert( IN_RANGE( output_data[ 0 ], -3e5, 0.0 ) );
              assert( IN_RANGE( output_data[ subset_points - 1 ], -3e5, 0.0 ));
              output_data += subset_points;
            }
          }
        }

        /* Copy data variable: */

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          const double seconds = data->seconds[ timestep ];

          if ( seconds != MISSING ) {
            copy_unmasked_float_values( points, mask, input_data, output_data);
            output_data += subset_points;
            input_data  += points;
          }
        }

        /* Copy additional component of data variable: */

        if ( input_data_2 ) {

          for ( timestep = 0; timestep < timesteps; ++timestep ) {
            const double seconds = data->seconds[ timestep ];

            if ( seconds != MISSING ) {
              copy_unmasked_float_values( points, mask, input_data_2,
                                          output_data );
              output_data += subset_points;
              input_data_2 += points;
            }
          }
        }
      }
    }
  }

  return result;
}



/* Read and convert CREOFS data: */

static OutputData* read_CREOFS( const Arguments* arguments, FILE* file ) {
  OutputData* result = 0;
  CREOFSData* data = read_CREOFS_header( arguments, file );

  if ( data ) {

    if ( read_CREOFS_data( arguments, file, data ) ) {
      result = convert_CREOFS( arguments, data );
    }

    deallocate_creofs_data( data ), data = 0;
  }

  return result;
}



/* Read/check DODS header for CREOFS file that looks like a subset of:
Dataset {
    Float32 lon[node = 74061];
    Float32 lat[node = 74061];
    Float32 h[node = 74061];
    Float32 sigma[sigma = 1];
    Float32 zeta[time = 1][node = 74061];
    Float32 Pair[time = 1][node = 74061];
    Float32 uwind_speed[time = 1][node = 74061];
    Float32 vwind_speed[time = 1][node = 74061];
    Float32 temp[time = 1][nv = 1][node = 74061];
    Float32 salinity[time = 1][nv = 1][node = 74061];
    Float32 u[time = 1][nv = 1][node = 74061];
    Float32 v[time = 1][nv = 1][node = 74061];
} NOAA%2fCREOFS%2fMODELS%2f201210%2fnos%2ecreofs%2efields%2en00%2e\
20121001%2et21z%2enc;

Data:
*/

static CREOFSData* read_CREOFS_header( const Arguments* arguments,
                                       FILE* file ) {
  const char* const variable = arguments ? arguments->variable : "";
  const int has_depth =
    strcmp( variable, "current" ) == 0 ||
    strcmp( variable, "water_temperature" ) == 0 ||
    strcmp( variable, "salinity" ) == 0;
  CREOFSData* result = 0;
  int file_base_yyyymmddhh = 0;
  int points = 0;
  int ok = 0;
  Line line = "";
  memset( line, 0, sizeof line );
  assert( arguments ); assert( file );

  ok = read_matched_line( file, "Dataset {\n" );

  /* Read lon and its dimensions: */

  ok = ok && read_dimension1_line( file, "Float32 lon[node = %d];\n", &points);

  /* Read lat: */

  snprintf( line, sizeof line / sizeof *line,
            "Float32 lat[node = %d];\n", points );
  ok = ok && read_matched_line( file, line );

  if ( has_depth ) { /* Read h, sigma, zeta: */
    snprintf( line, sizeof line / sizeof *line,
              "Float32 h[node = %d];\n", points );
    ok = ok && read_matched_line( file, line );
    ok = ok && read_matched_line( file, "Float32 sigma[sigma = 1];\n" );
    snprintf( line, sizeof line / sizeof *line,
              "Float32 zeta[time = 1][node = %d];\n", points );
    ok = ok && read_matched_line( file, line );
  }

  /* Read main variable: */

  if ( ok ) {

    if ( ! strcmp( variable, "air_pressure" ) ) { /* Read Pair: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 Pair[time = 1][node = %d];\n", points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "wind" ) ) { /* Read uwind_speed, vwind: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 uwind_speed[time = 1][node = %d];\n", points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 vwind_speed[time = 1][node = %d];\n", points );
      ok = ok && read_matched_line( file, line );
    } else if ( ! strcmp( variable, "water_temperature" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 temp[time = 1][nv = 1][node = %d];\n", points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "salinity" ) ) {
      snprintf( line, sizeof line / sizeof *line,
                "Float32 salinity[time = 1][nv = 1][node = %d];\n", points );
      ok = read_matched_line( file, line );
    } else if ( ! strcmp( variable, "current" ) ) { /* Read u, v: */
      snprintf( line, sizeof line / sizeof *line,
                "Float32 u[time = 1][nv = 1][node = %d];\n", points );
      ok = read_matched_line( file, line );
      snprintf( line, sizeof line / sizeof *line,
                "Float32 v[time = 1][nv = 1][node = %d];\n", points );
      ok = ok && read_matched_line( file, line );
    } else {
      fprintf( stderr, "Invalid variable '%s'.\n", variable );
      ok = 0;
    }

    /* Read end of header: */

    DEBUG( fprintf( stderr, "before end of header: ok = %d\n", ok ); )

    if ( ok ) {
      int end_yyyymmddhh = arguments->yyyymmddhh;
      increment_yyyymmddhh( &end_yyyymmddhh, arguments->hours );
      file_base_yyyymmddhh = read_file_base_yyyymmddhh( file );
      ok = file_base_yyyymmddhh >= arguments->yyyymmddhh &&
           file_base_yyyymmddhh < end_yyyymmddhh &&
           read_ignored_line( file ) &&
           read_ignored_line( file );
    }

    DEBUG( fprintf( stderr, "after end of header: file_base_yyyymmddhh = %d, "
                    "ok = %d\n", file_base_yyyymmddhh, ok ); )
  }

  if ( ! ok ) {

    if ( file_base_yyyymmddhh == 0 ) {
      fprintf( stderr, "\nInvalid input DODS header.\n" );
    }

  } else {
    result = NEW( CREOFSData, 1 );

    if ( result ) {
      result->file_base_yyyymmddhh = file_base_yyyymmddhh;
      result->points = points;
      DEBUG( fprintf( stderr, "points = %d\n", result->points ); )
    }
  }

  DEBUG( fprintf( stderr, "read_CREOFS_header result = %p\n", result ); )
  return result;
}



/* Read CREOFS data: */

static int read_CREOFS_data( const Arguments* arguments, FILE* file,
                             CREOFSData* data ) {
  int ok = 1;
  assert( arguments ); assert( data ); assert( file );
  assert( data->points > 0 );
  DEBUG( fprintf( stderr, "read_CREOFS_data: points = %d\n", data->points ); )

  {
    const char* const variable = arguments->variable;
    const int has_depth =
      strcmp( variable, "current" ) == 0 ||
      strcmp( variable, "water_temperature" ) == 0 ||
      strcmp( variable, "salinity" ) == 0;
    const int points = data->points;

    DEBUG( fprintf( stderr, "has_depth = %d\n", has_depth ); )

    /* Read longitudes: */

    data->longitude =
      read_new_float_data( file, 1.0, 0.0, -180.0, 180.0, points );
    ok = data->longitude != 0;

    /* Read latitudes: */

    if ( ok ) {
      data->latitude =
        read_new_float_data( file, 1.0, 0.0, -90.0, 90.0, points );
      ok = data->latitude != 0;
    }

    /* If data has depth then read h, sigma, and zeta: */

    if ( ok && has_depth ) {
      data->h = read_new_float_data( file, 1.0, 0.0, 0.0, 1.1e4, points );
      ok = data->h != 0;

      if ( ok ) {
        ok = read_float_data( file, 1.0, 0.0, -1.0, 0.0, 1, &data->s );
        ok = data->s != 0;

        if ( ok ) {
          data->msl = read_new_float_data(file, 1.0, 0.0, -50.0, 50.0, points);
          ok = data->msl != 0;
        }
      }
    }

    /* Read main variable: */

    if ( ok ) {

      if ( ! strcmp( variable, "wind" ) ) {
        data->wind_u =
          read_new_float_data( file, 1.0, 0.0, -100.0, 100.0, points );
        ok = data->wind_u != 0;

        if ( ok ) {
          data->wind_v =
            read_new_float_data( file, 1.0, 0.0, -100.0, 100.0, points );
          ok = data->wind_v != 0;
        }
      } else if ( ! strcmp( variable, "current" ) ) {
        data->current_u =
          read_new_float_data( file, 1.0, 0.0, -50.0, 50.0, points );
        ok = data->current_u != 0;

        if ( ok ) {
          data->current_v =
            read_new_float_data( file, 1.0, 0.0, -50.0, 50.0, points );
          ok = data->current_v != 0;
        }

      } else if ( ! strcmp( variable, "water_temperature" ) ) {
        data->water_temperature =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0, points );
        ok = data->water_temperature != 0;
      } else if ( ! strcmp( variable, "salinity" ) ) {
        data->salinity =
          read_new_float_data( file, 1.0, 0.0, 0.0, 50.0, points );
        ok = data->salinity != 0;
      } else {
        const double Pa_to_hPa = 1e-2;
        assert( ! strcmp( variable, "air_pressure" ) );
        data->air_pressure =
          read_new_float_data( file, Pa_to_hPa, 0.0, 500.0, 1500.0, points );
        ok = data->air_pressure != 0;
      }
    }

    /* Create mask: */

    if ( ok ) {
      data->mask = malloc( points * sizeof (float) );
      ok = data->mask != 0;

      if ( ok ) {
        int point = 0;

        for ( point = 0; point < points; ++point ) {
          data->mask[ point ] = 1.0;
        }
      }
    }
  }

  DEBUG( fprintf( stderr, "read_CREOFS_data result = %d\n", ok ); )
  return ok;
}



/* Convert CREOFS data: */

static OutputData* convert_CREOFS( const Arguments* arguments,
                                   CREOFSData* data ) {
  OutputData* result = 0;
  const int subset_points =
    spatially_subset( arguments->longitude_minimum,
                      arguments->longitude_maximum,
                      arguments->latitude_minimum,
                      arguments->latitude_maximum,
                      data->points,
                      sizeof data->mask[ 0 ],
                      data->longitude,
                      data->latitude,
                      data->mask );
  DEBUG( fprintf( stderr, "subset points = %d\n", subset_points ); )

  if ( subset_points ) {
    result = NEW( OutputData, 1 );

    if ( result ) {
      const size_t length = sizeof (Name) / sizeof (char);
      result->points = subset_points;
      result->timesteps = 1;

      if ( data->wind_u ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "wind_u", length );
        strncpy( result->names[ 1 ], "wind_v", length );
        strncpy( result->units[ 0 ], "m/s", length );
        strncpy( result->units[ 1 ], "m/s", length );
      } else if ( data->current_u ) {
        result->variables = 3;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "current_u", length );
        strncpy( result->names[ 2 ], "current_v", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "m/s", length );
        strncpy( result->units[ 2 ], "m/s", length );
      } else if ( data->water_temperature ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "water_temperature", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "C", length );
      } else if ( data->salinity ) {
        result->variables = 2;
        strncpy( result->names[ 0 ], "depth", length );
        strncpy( result->names[ 1 ], "salinity", length );
        strncpy( result->units[ 0 ], "m", length );
        strncpy( result->units[ 1 ], "PSU", length );
      } else {
        assert( data->air_pressure );
        result->variables = 1;
        strncpy( result->names[ 0 ], "air_pressure", length );
        strncpy( result->units[ 0 ], "hPa", length );
      }

      result->yyyymmddhh = NEW( int, 1 );
      result->longitudes = result->yyyymmddhh ? NEW( float, subset_points) : 0;
      result->latitudes  = result->longitudes ? NEW( float, subset_points) : 0;
      result->data = result->latitudes ?
        NEW( float, result->variables * subset_points ) : 0;

      if ( result->data ) {
        float* output_data = result->data;
        const int points = data->points;
        const float* const mask = data->mask;
        const float* const longitudes = data->longitude;
        const float* const latitudes = data->latitude;
        const float* input_data =
          data->wind_u ? data->wind_u
          : data->current_u ? data->current_u
          : data->water_temperature ? data->water_temperature
          : data->salinity ? data->salinity
          : data->air_pressure;
        const float* input_data_2 =
          data->wind_v ? data->wind_v : data->current_v;
        assert( input_data );
        assert( is_valid_yyyymmddhh( data->file_base_yyyymmddhh ) );

        /* Copy timestamps: */

        result->yyyymmddhh[ 0 ] = data->file_base_yyyymmddhh;

        /* Copy longitudes: */

        copy_unmasked_float_values( points, mask, longitudes,
                                    result->longitudes );
        assert( IN_RANGE( result->longitudes[ 0 ], -180.0, 180.0 ) );
        assert( IN_RANGE( result->longitudes[ subset_points - 1],
                          -180.0, 180.0 ) );

        /* Copy latitudes: */

        copy_unmasked_float_values( points, mask, latitudes,
                                    result->latitudes );
        assert( IN_RANGE( result->latitudes[ 0 ], -90.0, 90.0 ) );
        assert( IN_RANGE( result->latitudes[ subset_points - 1 ],
                          -90.0, 90.0 ) );

        /* Copy depth if present: */

        if ( data->s ) {
          assert( data->h ); assert( data->msl );
          copy_float_depth( points, mask, data->s, 0, data->h, data->msl,
                            output_data );
          assert( IN_RANGE( output_data[ 0 ], -3e5, 0.0 ) );
          assert( IN_RANGE( output_data[ subset_points - 1 ], -3e5, 0.0 ) );
          output_data += subset_points;
        }

        /* Copy data variable: */

        copy_unmasked_float_values( points, mask, input_data, output_data );
        output_data += subset_points;
        input_data  += points;

        /* Copy additional component of data variable: */

        if ( input_data_2 ) {
          copy_unmasked_float_values( points, mask, input_data_2, output_data);
          output_data += subset_points;
          input_data_2 += points;
        }
      }
    }
  }

  return result;
}



/* Reduce mask by user-specified lon-lat bounds: */

static int spatially_subset( const double longitude_minimum,
                             const double longitude_maximum,
                             const double latitude_minimum,
                             const double latitude_maximum,
                             const int count,
                             const int word_size,
                             const void* longitudes,
                             const void* latitudes,
                             void* mask ) {
  int result = 0;
  float* const fmask = mask;
  const float* const flongitudes = longitudes;
  const float* const flatitudes = latitudes;
  double* const dmask = mask;
  const double* const dlongitudes = longitudes;
  const double* const dlatitudes = latitudes;
  int index = 0;
  assert( IN_RANGE( longitude_minimum, -180.0, longitude_maximum ) );
  assert( IN_RANGE( longitude_maximum, longitude_minimum, 180.0 ) );
  assert( word_size == sizeof (float) || word_size == sizeof (double) );
  assert( mask ); assert( longitudes ); assert( latitudes );
  assert( count > 0 );
  DEBUG( fprintf( stderr,
                  "spatially_subset: word_size = %d, "
                  "[%lf %lf][%lf %lf]\n",
                  word_size,
                  longitude_minimum,
                  longitude_maximum,
                  latitude_minimum,
                  latitude_maximum ); )

  for ( index = 0; index < count; ++index ) {
    int m =
      word_size == sizeof (double) ? dmask[ index ] != 0.0
      : fmask[ index ] != 0.0;

    if ( m ) {
      const double longitude =
        word_size == sizeof (double) ? dlongitudes[ index ]
      : (double) ( flongitudes[ index ] );

      if ( ! IN_RANGE( longitude, longitude_minimum, longitude_maximum ) ) {
        m = 0;
      } else {
        const double latitude =
          word_size == sizeof (double) ? dlatitudes[ index ]
        : (double) ( flatitudes[ index ] );

        if ( ! IN_RANGE( latitude, latitude_minimum, latitude_maximum ) ) {
          m = 0;
        }
      }

      if ( word_size == sizeof (double) ) {
        dmask[ index ] = m;
      } else {
        fmask[ index ] = m;
      }

      result += m;
    }
  }

  DEBUG( fprintf( stderr, "spatially_subset = %d (of %d)\n",
                 result, count ); )
  return result;
}



/* Copy unmasked values to output[]: */

static void copy_unmasked_values( const int count, const double mask[],
                                  const double values[], float output[] ) {
  int index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( mask[ index ] ) {
      *output++ = values[ index ];
    }
  }
}



/* Copy unmasked values to output[]: */

static void copy_unmasked_values2( const int count, const double mask[],
                                   const float values[], float output[] ) {

  int index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( mask[ index ] ) {
      *output++ = values[ index ];
    }
  }
}



/* Copy unmasked values to output[]: */

static void copy_unmasked_float_values( const int count, const float mask[],
                                        const float values[], float output[]) {
  int index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( mask[ index ] ) {
      *output++ = values[ index ];
    }
  }
}



/*
 * Copy unmasked computed depth (m) from
 * s = normalized s-coordinate. s = -0.975 near surface, -0.025 near bottom.
 * h = depth in meters of bottom. h > 0.
 * msl = + or - meters from mean sea level of surface due to tide.
 */

static void copy_depth( const int count, const double mask[],
                        const double s, const double h[],
                        const float msl[], float output[] ) {
  assert( s > -1.0 ); assert( s < 0.0 );
  int index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( mask[ index ] ) {
      const double bottom = h[ index ];
      const double top = msl[ index ];
      const double depth = s * ( bottom - top );

      if ( bottom > 0.0 && bottom - top > 0.0 && depth < 0.0 ) {
        *output++ = depth;
      } else {
        *output++ = MISSING;
      }
    }
  }
}



/*
 * Copy unmasked computed depth (m) from
 * s = normalized s-coordinate. s = -0.975 near surface, -0.025 near bottom.
 * h = depth in meters of bottom. h > 0.
 * msl = + or - meters from mean sea level of surface due to tide.
 */

static void copy_float_depth( const int count, const float mask[],
                              const double s, const float s2[],
                              const float h[],
                              const float msl[], float output[] ) {
  int index = 0;
  assert( s >= -1.0 ); assert( s <= 0.0 );

  for ( index = 0; index < count; ++index ) {

    if ( mask[ index ] ) {
      const double bottom = h[ index ];
      const double top = msl[ index ];
      const double s0 = s2 ? s2[ index ] : s;
      const double depth0 = s0 * ( bottom - top );
      const double depth = depth0 <= 0.0 ? depth0 : -depth0;
      assert( depth <= 0.0 );
      *output++ = depth;
    }
  }
}



/*
 * Copy unmasked computed depth (m) from
 * s = normalized s-coordinate. s = -0.975 near surface, -0.025 near bottom.
 * h = depth in meters of bottom. h > 0.
 * msl = + or - meters from mean sea level of surface due to tide.
 * node[ 3 * elements ] - 3 nodes for a given element.
 * mask[ elements ] - 1 if element is not missing.
 * s[ nodes ] - sigma per node.
 * h[ nodes ] - h per node.
 * msl[ nodes ] - msl per node.
 */

static void copy_float_depth3( const int elements,
                               const float node[],
                               const float mask[],
                               const float s[],
                               const float h[],
                               const float msl[],
                               float output[] ) {
  const double one_third = 1.0 / 3.0;
  double previous_depth = 0.0;
  const int elements2 = elements + elements;
  int element = 0;
  assert( s[0] >= -1.0 ); assert( s[0] <= 0.0 );

  for ( element = 0; element < elements; ++element ) {

    if ( mask[ element ]  ) {
      const int node1 = node[ element ];
      const int node2 = node[ elements + element ];
      const int node3 = node[ elements2 + element ];

      const double bottom1 = h[ node1 ];
      const double top1 = msl[ node1 ];
      const double s1 = s[ node1 ];
      const double depth1 = s1 * ( bottom1 - top1 );

      const double bottom2 = h[ node2 ];
      const double top2 = msl[ node2 ];
      const double s2 = s[ node2 ];
      const double depth2 = s2 * ( bottom2 - top2 );

      const double bottom3 = h[ node3 ];
      const double top3 = msl[ node3 ];
      const double s3 = s[ node3 ];
      const double depth3 = s3 * ( bottom3 - top3 );

      const double depth = ( depth1 + depth2 + depth3 ) * one_third;

/*
      assert( bottom1 > 0.0 );
      assert( bottom1 - top1 > 0.0 );
      assert( depth1 <= 0.0 );

      assert( bottom2 > 0.0 );
      assert( bottom2 - top2 > 0.0 );
      assert( depth2 <= 0.0 );

      assert( bottom3 > 0.0 );
      assert( bottom3 - top3 > 0.0 );
      assert( depth3 <= 0.0 );

      assert( depth <= 0.0 );

      *output++ = depth;
*/

      if ( ! ( bottom1 > 0.0 &&
               bottom1 - top1 > 0.0 &&
               depth1 <= 0.0 &&
               bottom2 > 0.0 &&
               bottom2 - top2 > 0.0 &&
               depth2 <= 0.0 &&
               bottom3 > 0.0 &&
               bottom3 - top3 > 0.0 &&
               depth3 <= 0.0 ) ) {
        DEBUG( fprintf( stderr, "Invalid depth: element %d, depth = %lf "
                        "[%lf %lf] = %lf, [%lf %lf] = %lf, [%lf %lf] = %lf\n",
                        element, depth,
                        bottom1, top1, depth1,
                        bottom2, top2, depth2,
                        bottom3, top3, depth3 ); )
        *output++ = previous_depth;
      } else {
        *output++ = depth;
        previous_depth = depth;
      }
    }
  }
}



/* Convert and copy unmasked wind: */

static void copy_wind( const int count, const double mask[],
                       const double angle[],
                       const float u[], const float v[],
                       float u_east[], float v_north[] ) {
  int index = 0;

  for ( index = 0; index < count; ++index ) {

    if ( mask[ index ] ) {
      transform_point2( angle[ index ], u[ index ], v[index], u_east, v_north);
      u_east++;
      v_north++;
    }
  }
}



/*
 Convert and copy unmasked current:
    Float64 s_rho[s_rho = 1];
    Float64 h[eta_rho = 291][xi_rho = 332];
    Float64 lon_rho[eta_rho = 291][xi_rho = 332];
    Float64 lat_rho[eta_rho = 291][xi_rho = 332];
    Float64 angle[eta_rho = 291][xi_rho = 332];
    Float64 mask_rho[eta_rho = 291][xi_rho = 332];
    Float64 mask_u[eta_u = 291][xi_u = 331];
    Float64 mask_v[eta_v = 290][xi_v = 332];
    Float64 ocean_time[ocean_time = 7];
    Float64 zeta[ocean_time = 7][eta_rho = 291][xi_rho = 332];
    Float64 u[ocean_time = 7][s_rho = 1][eta_u = 291][xi_u = 331];
    Float64 v[ocean_time = 7][s_rho = 1][eta_v = 290][xi_v = 332];
    Float64 w[ocean_time = 7][s_w = 2][eta_rho = 291][xi_rho = 332];
 */

static void copy_current( const int timestep, const CBOFSData* data,
                          float u_east[], float v_north[], float w_up[] ) {
  const int rows = data->rows;
  const int columns = data->columns;
  const int rows_1 = rows - 1;
  const int columns_1 = columns - 1;
  const int rows_columns = rows * columns;
  const double* const mask = data->mask;
  const double* const angle = data->angle;
  const float* const u = data->current_u + timestep * rows * columns_1;
  const float* const v = data->current_v + timestep * rows_1 * columns;
  const float* const w = data->current_w + timestep * 2 * rows_columns;
  int index = 0;
  int row = 0;

  for ( row = 0; row < rows; ++row ) {
    const int v_row = row < rows_1 ? row : row - 1;
    int column = 0;

    for ( column = 0; column < columns; ++column, ++index ) {

      if ( mask[ index ] ) {

        /* Compute mean u from current u and previous u: */

        const int u_column = column < columns_1 ? column : column - 1;
        const int u_index = row * columns_1 + u_column;
        const double u_grid1 = u[ u_index ];
        const double u_grid2 =
        ( column == 0 || column == columns_1 ) ? u_grid1 : u[ u_index - 1 ];
        const double u_grid =
          ( u_grid1 >= -50.0 && u_grid2 >= -50.0 ) ? 0.5 * ( u_grid1 + u_grid2)
          : MISSING;

        /* Compute mean v from current v and previous v: */

        const int v_index = v_row * columns + column;
        const double v_grid1 = v[ v_index ];
        const double v_grid2 =
          (row == 0 || row == rows_1) ? v_grid1 : v[ v_index - columns ];
        const double v_grid =
          ( v_grid1 >= -50.0 && v_grid2 >= -50.0 ) ? 0.5 * ( v_grid1 + v_grid2)
          : MISSING;

        /* Compute mean w from current w and next w: */

        const double w_grid1 = w[ index ];
        const double w_grid2 = w[ index + rows_columns ];
        const double w_grid  =
          ( w_grid1 >= -50.0 && w_grid2 >= -50.0 ) ? 0.5 * ( w_grid1 + w_grid2)
          : MISSING;

        if ( u_grid >= -50.0 && v_grid >= -50.0 && w_grid >= -50.0 ) {
          transform_point2( angle[ index ], u_grid, v_grid, u_east, v_north );
          ++u_east;
          ++v_north;
          *w_up++ = w_grid;
        } else {
          *u_east++  = MISSING;
          *v_north++ = MISSING;
          *w_up++    = MISSING;
        }
      }
    }
  }
}



/* Read and match a line: */

static int read_ignored_line( FILE* file ) {
  int result = 0;
  Line line = "";
  memset( line, 0, sizeof line );

  result = fgets( line, sizeof line / sizeof *line, file ) != 0;
  DEBUG( fprintf( stderr, "read line = '%s'\n", line ); )

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read an input line.\n" );
  }

  DEBUG( fprintf( stderr, "read_ignored_line result = %d\n", result ); )
  return result;
}



/* Read and match a line: */

static int read_matched_line( FILE* file, const char* const match ) {
  int result = 0;
  Line line = "";
  memset( line, 0, sizeof line );

  if ( fgets( line, sizeof line / sizeof *line, file ) ) {
    char* s = &line[ 0 ];

    while ( isspace( *s ) ) {
      ++s;
    }

    result = ! strcmp( s, match );
  }

  DEBUG( fprintf( stderr, "match = '%s', read line = '%s'\n", match, line ); )

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read valid input line '%s'\n", match );
  }

  DEBUG( fprintf( stderr, "read_matched_line result = %d\n", result ); )
  return result;
}



/* Read line with 1-dimensioned array: */

static int read_dimension1_line( FILE* file, const char* const format,
                                 int* dimension ) {
#if 0
  const int result = fscanf( file, format, dimension ) == 1 && *dimension > 0;
#else
  int result = 0;
  Line line = "";
  memset( line, 0, sizeof line );

  if ( fgets( line, sizeof line / sizeof *line, file ) ) {
    char* s = &line[ 0 ];

    while ( isspace( *s ) ) {
      ++s;
    }

    result = sscanf( s, format, dimension ) == 1 && *dimension > 0;
  }
#endif

  DEBUG( fprintf( stderr, "format = '%s', *dimension = %d\n",
                  format, *dimension ); )

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read valid input line '%s'\n", format );
  }

  DEBUG( fprintf( stderr, "read_dimension1_line result = %d\n", result ); )
  return result;
}



/* Read line with 2-dimensioned array: */

static int read_dimension2_line( FILE* file, const char* const format,
                                 int* dimension1, int* dimension2 ) {
#if 0
  const int result =
    fscanf( file, format, dimension1, dimension2 ) == 2 &&
    *dimension1 > 0 && *dimension2 > 0;
#else
  int result = 0;
  Line line = "";
  memset( line, 0, sizeof line );

  if ( fgets( line, sizeof line / sizeof *line, file ) ) {
    char* s = &line[ 0 ];

    while ( isspace( *s ) ) {
      ++s;
    }

    result = sscanf( s, format, dimension1, dimension2 ) == 2 &&
             *dimension1 > 0 && *dimension2 > 0;
  }
#endif

  DEBUG( fprintf( stderr,
                  "format = '%s', *dimension1 = %d, *dimension2 = %d\n",
                  format, *dimension1, *dimension2 ); )

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read valid input line '%s'\n", format );
  }

  DEBUG( fprintf( stderr, "read_dimension2_line result = %d\n", result ); )
  return result;
}



/*
 * Read line and parse the file timestamp:
"} NOAA%2fSJROFS%2fMODELS%2f201210%2fnos%2esjrofs%2efields%2enowcast\
%2e20121022%2et00z%2enc;\n"
result = 2012102200
*/

static int read_file_base_yyyymmddhh( FILE* file ) {
  int result = 0;
  Line line = "";
  memset( line, 0, sizeof line );

  if ( fgets( line, sizeof line / sizeof *line, file ) ) {
    const char* const t = strrchr( line, 't' );
    const int isPercent = strrchr( line, '%' ) != 0;
    DEBUG( fprintf( stderr, "read_file_base_yyyymmddhh line = '%s'\n", line );)

    if ( t ) {
      const int hh = atoi( t + 1 );
      const int dateOffset = isPercent ? 11 : 9;
      const int yyyymmdd = atoi( t - dateOffset );
      const int is_creofs = strstr( line, "CREOFS" ) != 0;
      result = yyyymmdd * 100 + hh;

      if ( ! is_valid_yyyymmddhh( result ) ) {
        result = 0;
      } else if ( is_creofs ) {
        const int nnOffset = isPercent ? 16 : 13;
        const int nn = 6 - atoi( t - nnOffset );

        if ( nn > 0 ) {
          decrement_yyyymmddhh( &result, nn );
          assert( is_valid_yyyymmddhh( result ) );          
        }
      }
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read valid input line '%s'\n", line );
  }

  DEBUG( fprintf( stderr, "read_file_base_yyyymmddhh result = %d\n", result );)
  assert( result == 0 || is_valid_yyyymmddhh( result ) );
  return result;
}



/* Allocate and read input file double data: */

static double* read_new_double_data( FILE* file,
                                     const double scale,
                                     const double offset,
                                     const double minimum,
                                     const double maximum,
                                     const size_t count ) {
  double* result = 0;
  assert( file ); assert( minimum <= maximum ); assert( count );
  result = NEW( double, count );

  if ( result ) {
    const int ok =
      read_double_data( file, scale, offset, minimum, maximum, count, result );

    if ( ! ok ) {
      FREE( result );
    }
  }

  return result;
}



/* Allocate and read input file float data: */

static float* read_new_float_data( FILE* file,
                                   const double scale, const double offset,
                                   const double minimum, const double maximum,
                                   const size_t count ) {
  float* result = 0;
  assert( file ); assert( minimum <= maximum ); assert( count );
  result = NEW( float, count );

  if ( result ) {
    const int ok =
      read_float_data( file, scale, offset, minimum, maximum, count, result );

    if ( ! ok ) {
      FREE( result );
    }
  }

  return result;
}



/* Allocate and read input file int data: */

static float* read_new_int_data( FILE* file,
                                 const double scale, const double offset,
                                 const double minimum, const double maximum,
                                 const size_t count ) {
  float* result = 0;
  assert( file ); assert( minimum <= maximum ); assert( count );
  result = NEW( float, count );

  if ( result ) {
    const int ok =
      read_int_data( file, scale, offset, minimum, maximum, count, result );

    if ( ! ok ) {
      FREE( result );
    }
  }

  return result;
}



/* Read input file double data: */

static int read_double_data( FILE* file,
                             const double scale, const double offset,
                             const double minimum, const double maximum,
                             const size_t count, double data[] ) {
  int result = skip_8_bytes( file );

  if ( result ) {
    result = fread( data, sizeof (double), count, file ) == count;
    DEBUG( fprintf( stderr, "Read %lu doubles = %d [%lf %lf].\n",
                    count, result, minimum, maximum ); )

    if ( result ) {
      size_t index = 0;
      result = 0;
      reverse_8_byte_words_if_little_endian( count, data );
      DEBUG( fprintf( stderr, "data[] = [%lf, ... %lf, ... %lf]\n",
                      *data, data[ count / 2 ], data[ count - 1 ] ); )

      for ( index = 0; index < count; ++index ) {
        double value = data[ index ];
        DEBUG( if ( count < 50 ) fprintf( stderr, "%lf -> ", value ); )

        if ( value != MISSING ) {
          value *= scale;
          value += offset;

          if ( ! IN_RANGE( value, minimum, maximum ) ) {
            value = MISSING;
          } else {
            ++result;
          }

          data[ index ] = value;
        }

        DEBUG( if ( count < 50 ) fprintf( stderr, "%lf\n", data[ index ] ); )
      }
    }
  }

#ifdef DEBUGGING
  if ( ! result ) {
    fprintf( stderr, "\nFailed to read any DODS file data "
             "in range [%lf, %lf].\n", minimum, maximum );
  }
#endif

  DEBUG( fprintf( stderr, "read_double_data result = %d\n", result ); )
  return result;
}



/* Read input file float data: */

static int read_float_data( FILE* file,
                            const double scale, const double offset,
                            const double minimum, const double maximum,
                            const size_t count, float data[] ) {
  int result = skip_8_bytes( file );

  if ( result ) {
    result = fread( data, sizeof (float), count, file ) == count;
    DEBUG( fprintf( stderr, "Read %lu floats = %d [%lf %lf].\n",
                    count, result, minimum, maximum ); )

    if ( result ) {
      size_t index = 0;
      result = 0;
      reverse_4_byte_words_if_little_endian( count, data );
      DEBUG( fprintf( stderr, "data[] = [%f, ... %f, ... %f]\n",
                      *data, data[ count / 2 ], data[ count - 1 ] ); )

      for ( index = 0; index < count; ++index ) {
        double value = data[ index ];
        DEBUG( if ( count < 50 ) fprintf( stderr, "%lf -> ", value ); )

        if ( value != MISSING ) {
          value *= scale;
          value += offset;

          if ( ! IN_RANGE( value, minimum, maximum ) ) {
            value = MISSING;
          } else {
            ++result;
          }

          data[ index ] = value;
        }

        DEBUG( if ( count < 50 ) fprintf( stderr, "%f\n", data[ index ] ); )
      }
    }
  }

#ifdef DEBUGGING
  if ( ! result ) {
    fprintf( stderr, "\nFailed to read any DODS file data "
             "in range [%lf, %lf].\n", minimum, maximum );
  }
#endif

  DEBUG( fprintf( stderr, "read_float_data result = %d\n", result ); )
  return result;
}



/* Read input file float data: */

static int read_int_data( FILE* file,
                          const double scale, const double offset,
                          const double minimum, const double maximum,
                          const size_t count, float data[] ) {
  int result = skip_8_bytes( file );

  if ( result ) {
    int* const idata = (int*) data;
    result = fread( idata, sizeof (int), count, file ) == count;
    DEBUG( fprintf( stderr, "Read %lu ints = %d [%lf %lf].\n",
                    count, result, minimum, maximum ); )

    if ( result ) {
      size_t index = 0;
      result = 0;
      reverse_4_byte_words_if_little_endian( count, idata );
      DEBUG( fprintf( stderr, "idata[] = [%d, ... %d, ... %d]\n",
                      idata[ 0 ], idata[ count / 2 ], idata[ count - 1 ] ); )

      for ( index = 0; index < count; ++index ) {
        int ivalue = idata[ index ];
        float value = ivalue;
        DEBUG( if ( count < 50 ) fprintf( stderr, "%f -> ", value ); )

        if ( value != MISSING ) {
          value *= scale;
          value += offset;

          if ( ! IN_RANGE( value, minimum, maximum ) ) {
            value = MISSING;
          } else {
            ++result;
          }

          data[ index ] = value;
        }

        DEBUG( if ( count < 50 ) fprintf( stderr, "%f\n", data[ index ] ); )
      }
    }
  }

#ifdef DEBUGGING
  if ( ! result ) {
    fprintf( stderr, "\nFailed to read any DODS file data "
             "in range [%lf, %lf].\n", minimum, maximum );
  }
#endif

  DEBUG( fprintf( stderr, "read_int_data result = %d\n", result ); )
  return result;
}



/* Read and ignore 8 bytes from file: */

static int skip_8_bytes( FILE* file ) {
  double unused = 0.0;
  const int result = fread( &unused, 8, 1, file ) == 1;
  return result;
}



/* On little-endian platforms, change endianess of an array of 4-byte words: */

static void reverse_4_byte_words_if_little_endian( const size_t count,
                                                   void* array ) {

#if IS_LITTLE_ENDIAN

  unsigned int* const array4 = array;
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const unsigned int value = array4[ index ];
    const unsigned int new_value =
      ( value & 0xff000000 ) >> 24 |
      ( value & 0x00ff0000 ) >>  8 |
      ( value & 0x0000ff00 ) <<  8 |
      ( value & 0x000000ff ) << 24;
    array4[ index ] = new_value;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/* On little-endian platforms, change endianess of an array of 8-byte words: */

static void reverse_8_byte_words_if_little_endian( const size_t count,
                                                   void* array ) {

#if IS_LITTLE_ENDIAN

  unsigned long long* const array8 = array;
  size_t index = 0;

  for ( index = 0; index < count; ++index ) {
    const unsigned long long value = array8[ index ];
    const unsigned long long new_value =
      ( value & 0xff00000000000000LL ) >> 56 |
      ( value & 0x00ff000000000000LL ) >> 40 |
      ( value & 0x0000ff0000000000LL ) >> 24 |
      ( value & 0x000000ff00000000LL ) >>  8 |
      ( value & 0x00000000ff000000LL ) <<  8 |
      ( value & 0x0000000000ff0000LL ) << 24 |
      ( value & 0x000000000000ff00LL ) << 40 |
      ( value & 0x00000000000000ffLL ) << 56;
    array8[ index ] = new_value;
  }

#endif /* IS_LITTLE_ENDIAN */

}



/*
 * Given a 2D point (x, y) on coordinate-axes rotated by angle,
 * compute the point's coordinates (xp, yp) on non-rotated axes:
 */

static void transform_point2( const double angle,
                              const double x, const double y,
                              float* xp, float* yp ) {
  const double angle2 = atan2( y, x );
  const double r = hypot( x, y );
  const double angle_sum = angle + angle2;
  *xp = r * cos( angle_sum );
  *yp = r * sin( angle_sum );
}



/* Seconds difference between two timetamps: */

static double seconds_difference( const int yyyymmddhh1,
                                  const int yyyymmddhh0 ) {
  double result = 0.0;
  int yyyymmddhh = yyyymmddhh0;
  assert( yyyymmddhh1 >= yyyymmddhh0 );

  while ( yyyymmddhh < yyyymmddhh1 ) {
    increment_yyyymmddhh( &yyyymmddhh, 1 );
    result += SECONDS_PER_HOUR;
  }

  DEBUG( fprintf( stderr, "seconds_difference( %d, %d ) = %lf\n",
                  yyyymmddhh1, yyyymmddhh0, result ); )
  return result;
}



/*
 * The non-standard routine timegm(), when available, is not either thread-safe
 * or is thread-safe but employs mutex-locking that is too slow or even hangs!
 */

static time_t my_timegm( struct tm *t ) {
  time_t result = mktime( t );

  /* Compute seconds from 1970-01-01T00:00:00-0000 to get timezone offset: */

  time_t timezone_offset_seconds = 0;
  struct tm t0;
  memset( &t0, 0, sizeof t0 );
  t0.tm_mday = 1;
  t0.tm_year = 70;
  t0.tm_wday = 4; /* January 1, 1970 was Thu. Sun = 0, Mon = 1, ... Thu = 4. */
  timezone_offset_seconds = mktime( &t0 );

  result -= timezone_offset_seconds;
  return result;
}



/* Return yyyymmddhh = yyyymmddhh0 + seconds: */

static int convert_seconds( const int yyyymmddhh0, const int seconds ) {
  int yyyy = yyyymmddhh0 / 1000000;
  int mm   = yyyymmddhh0 / 10000 % 100;
  int dd   = yyyymmddhh0 / 100 % 100;
  int hh   = yyyymmddhh0 % 100;
  int result = 0;
  struct tm s;
  memset( &s, 0, sizeof s );
  s.tm_year = yyyy - 1900;
  s.tm_mon  = mm - 1;
  s.tm_mday = dd;
  s.tm_hour = hh;
  s.tm_sec  = seconds;
  s.tm_zone = "UTC";

  {
    const time_t t = my_timegm( &s );
    const struct tm* const s2 = gmtime( &t );
    yyyy = s2->tm_year + 1900;
    mm = s2->tm_mon + 1;
    dd = s2->tm_mday;
    hh = s2->tm_hour;
    result = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
  }

  DEBUG( fprintf( stderr, "%d + %d = %d\n", yyyymmddhh0, seconds, result ); )
  return result;
}



/* Is yyyymmddhh valid? */

static int is_valid_yyyymmddhh( const int yyyymmddhh ) {
  const int yyyy = yyyymmddhh / 1000000;
  const int mm   = yyyymmddhh / 10000 % 100;
  const int dd   = yyyymmddhh / 100 % 100;
  const int hh   = yyyymmddhh % 100;
  const int result =
    IN_RANGE( yyyy, 1900, 3000 ) &&
    IN_RANGE( mm, 1, 12 ) &&
    IN_RANGE( dd, 1, days_in_month( yyyy, mm ) ) &&
    IN_RANGE( hh, 0, 23 );
  return result;
}



/* Increment yyyymmddhh by hours: */

static void increment_yyyymmddhh( int* yyyymmddhh, const int hours ) {
  int yyyy = *yyyymmddhh / 1000000;
  int mm   = *yyyymmddhh / 10000 % 100;
  int dd   = *yyyymmddhh / 100 % 100;
  int hh   = *yyyymmddhh % 100;
  int hour = 0;
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) ); assert( hours > 0 );

  for ( hour = 0; hour < hours; ++hour ) {
    ++hh;

    if ( hh > 23 ) {
      hh = 0;
      ++dd;

      if ( dd > 28 && dd > days_in_month( yyyy, mm ) ) {
        dd = 1;
        ++mm;

        if ( mm > 12 ) {
          mm = 1;
          ++yyyy;
        }
      }
    }
  }

  *yyyymmddhh = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) );
}



/* Decrement yyyymmddhh by hours: */

static void decrement_yyyymmddhh( int* yyyymmddhh, const int hours ) {
  int yyyy = *yyyymmddhh / 1000000;
  int mm   = *yyyymmddhh / 10000 % 100;
  int dd   = *yyyymmddhh / 100 % 100;
  int hh   = *yyyymmddhh % 100;
  int hour = 0;
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) ); assert( hours > 0 );

  for ( hour = 0; hour < hours; ++hour ) {
    --hh;

    if ( hh < 0 ) {
      hh = 23;
      --dd;

      if ( dd < 1 ) {
        --mm;

        if ( mm < 1 ) {
          mm = 12;
          --yyyy;
        }

        dd = days_in_month( yyyy, mm );
      }
    }
  }

  *yyyymmddhh = yyyy * 1000000 + mm * 10000 + dd * 100 + hh;
  assert( is_valid_yyyymmddhh( *yyyymmddhh ) );
}



/* Number of days in year/month: */

static int days_in_month( const int year, const int month ) {
  static const int days_per_month[ 2 ][ 12 ] = {
    { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }, /* Non-leap year. */
    { 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 }  /* Leap year. */
  };
  const int leap =
    month != 2 ? 0 : year % 4 == 0 && ! ( year % 100 == 0 && year % 400 != 0 );
  const int result = days_per_month[ leap ][ month - 1 ];
  return result;
}






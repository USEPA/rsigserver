#ifndef PARAMETERS_H
#define PARAMETERS_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Parameters.h - Declare parameters for translator routines.

NOTES:

HISTORY: 2007/12, plessel.todd@epa.gov, Created.

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <Utilities.h> /* For Integer, Stream, Grid. */
#include <Helpers.h>   /* For CompareFunction. */

/*================================== TYPES ==================================*/

/* Output formats: */

enum { FORMAT_XDR, FORMAT_ASCII, FORMAT_COARDS, FORMAT_IOAPI };
  
#define IS_VALID_FORMAT( format ) IN_RANGE( (format), FORMAT_XDR, FORMAT_IOAPI )

typedef struct {
  Integer format;         /* Output format: FORMAT_XDR, FORMAT_ASCII, ... */
  Stream* input;          /* Opened, readable stream to read. */
  const char* temporaryDirectory; /* Name of writable directory for temp files*/
  char regridFileName[ 256 ]; /* Name of temporary regrid file to create.*/
  char netcdfFileName[ 256 ]; /* Name of temporary NetCDF file to create. */
  Grid*   grid;           /* Projects lon-lat-elv points onto grid. */
  Integer regrid;         /* Regrid: 0 or AGGREGATE_MEAN, etc,... */
  Integer aggregationTimesteps; /* 0, 24 or timesteps timesteps to mean. */
  Integer ok;             /* Did last command succeed? */
  CompareFunction compareFunction; /* Or 0 if not comparing. */
  ConvertFunction convertFunction; /* Or 0 if not converting. */
  UTCTimestamp timestamp; /* First CMAQ timestamp to compare to. */
  Name variable;          /* CMAQ data variable name. */
  Name units;             /* CMAQ data variable units. */
  Integer timesteps;      /* Number of CMAQ data timesteps. */
  Integer firstLayer;     /* 1-based subset layer number. */
  Integer lastLayer;      /* 1-based subset layer number. */
  Integer firstRow;       /* 1-based subset row number. */
  Integer lastRow;        /* 1-based subset row number. */
  Integer firstColumn;    /* 1-based subset column number. */
  Integer lastColumn;     /* 1-based subset column number. */
  Real* data;             /* data[ timesteps ][ rows ][ columns ]. */
  Real* data2;            /* data2[ timesteps ][ rows ][ columns ]. */
} Parameters;

/*================================ FUNCTIONS ================================*/

extern Integer isValidParameters( const Parameters* parameters );


#ifdef __cplusplus
}
#endif

#endif /* PARAMETERS_H */




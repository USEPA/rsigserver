
/******************************************************************************
PURPOSE: M3IO.c - Define convenience routines for M3IO files.

NOTES:   For a description of M3IO Conventions see:

HISTORY: 2008/02 plessel.todd@epa.gov, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h>  /* For stderr, fprintf(). */
#endif
#include <string.h> /* For memset(). */

#include <netcdf.h> /* For nc_*. */

#include <Utilities.h> /* For PRE*(), NEW_ZERO(), Integer, Real, Stream. */

#include <Helpers.h>         /* For Name. */
#include <NetCDFUtilities.h> /* For createDimension(), etc. */
#include <M3IO.h>            /* For public interface. */

/*================================== TYPES ==================================*/

/* From M3IO Library: */

enum {
  LATGRD3 = 1, LAMGRD3 = 2, POLGRD3 = 6, EQMGRD3 = 7,
  NAMLEN3 = 16, MXDLEN3 = 80, MXDESC3 = 60, MXVARS3 = 120
};

enum { TSTEP, DATE_TIME, LAY, VAR, ROW_DIM, COL, M3IO_DIMS };

/*========================== FORWARD DECLARATIONS ===========================*/

static Integer writeM3IODimensions( Integer file,
                                    Integer timesteps,
                                    Integer variables,
                                    Integer layers,
                                    Integer rows,
                                    Integer columns,
                                    Integer dimensions[ M3IO_DIMS ] );

static Integer writeM3IOVariables( Integer file, Integer variables,
                                   const Name variableNames[],
                                   const Name variableUnits[],
                                   const Integer dimensions[ M3IO_DIMS ] );

static Integer writeM3IOAttributes( Integer file,
                                   Integer hoursPerTimestep,
                                   Integer firstTimestamp,
                                    Integer variables,
                                    Integer layers,
                                    const Name variableNames[],
                                    const char* description,
                                    const Grid* grid );

static Integer writeVarListAttribute( Integer file,
                                      const Name variableNames[],
                                      Integer variables );

static Integer writeM3IOTFLAG( Integer file,
                               Integer variables,
                               Integer timesteps,
                               Integer hoursPerTimestep,
                               Integer firstTimestamp );

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: writeM3IOHeader - Write M3IO header (plus TFLAG variable).
INPUTS:  Integer file                NetCDF file to write to.
         Integer hoursPerTimestep    Number of hours per timesteps: 1, 24, etc.
         Integer timesteps           Number of timesteps.
         Integer firstTimestamp      First timestamp.
         Integer variables           Number of variables in variableNames array.
         Integer layers              Number of layers in grid.
         const Name variableNames[]  Names of variables.
         const Name variableUnits[]  Units of variables.
         const char* description     Description of the data.
         const Grid* grid            Grid description object.
RETURNS: Integer 1 if successful, else 0.
NOTES:   If unsuccessful then failureMessage() is called.
         Resulting NetCDF file will also contain TFLAG.
******************************************************************************/

Integer writeM3IOHeader( Integer file,
                         Integer timesteps,
                         Integer hoursPerTimestep,
                         Integer firstTimestamp,
                         Integer variables,
                         Integer layers,
                         const Name variableNames[],
                         const Name variableUnits[],
                         const char* description,
                         const Grid* grid ) {

  PRE010( file > -1,
          GT_ZERO4( timesteps, hoursPerTimestep, variables, layers ),
          variableNames, *variableNames, **variableNames,
          variableUnits, *variableUnits, **variableUnits,
          grid, grid->invariant( grid ) );

  const Integer rows    = grid->rows( grid );
  const Integer columns = grid->columns( grid );
  Integer dimensionIds[ M3IO_DIMS ] = { -1, -1, -1, -1, -1 };
  Integer result =
    writeM3IODimensions( file, timesteps, variables,
                         layers, rows, columns, dimensionIds );

  result = AND2( result,
                 writeM3IOVariables( file, variables,
                                     variableNames, variableUnits,
                                     dimensionIds ) );

  result = AND2( result,
                 writeM3IOAttributes( file, hoursPerTimestep, firstTimestamp,
                                      variables, layers, variableNames,
                                      description, grid ) );

  DEBUG( if ( result ) fprintf( stderr, "calling nc_enddef()...\n" ); )

  if ( result ) {
    const int status = nc_enddef( file ); /* SLOW! */
    DEBUG( fprintf( stderr, "...done\n" ); )
    result = status == NC_NOERR;

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't end definition because %s.", message );
    } else {
      result =
        writeM3IOTFLAG( file, variables, timesteps, hoursPerTimestep,
                        firstTimestamp );
    }
  }

  if ( ! result ) {
    nc_close( file );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeM3IOGrid - Write grid longitudes and latitudes to an M3IO file.
INPUTS:  const Grid* grid          The grid to write.
         Integer timesteps         Number of timesteps to write.
         Integer layers            Number of layers  to write.
         Integer file              Handle of NetCDF file to write to.
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeM3IOGrid( const Grid* grid, Integer timesteps, Integer layers,
                       Integer file ) {

  PRE05( grid, grid->invariant( grid ), timesteps > 0, layers > 0, file > -1);

  Integer result = 0;
  const Integer rows    = grid->rows( grid );
  const Integer columns = grid->columns( grid );
  const Integer cells   = rows * columns;
  float* data = NEW_ZERO( float, cells ); /* Copy to 32-bits... */

  if ( data ) {
    const char* const variableNames[] = { "LONGITUDE", "LATITUDE", "ELEVATION"};
    typedef Real (*CoordinateFunction)( const Grid* self,
                                        Integer row, Integer column );
    CoordinateFunction coordinateFunctions[] = {
      grid->longitude, grid->latitude, 0
    };
    const Integer variableCount = 2 + ( layers > 1 );
    Integer variableIndex = 0;
    int status = 0;
    size_t start[ 4 ] = { 0, 0, 0, 0 };
    size_t count[ 4 ] = { 1, 1, 0, 0 };
    count[ 2 ] = rows;
    count[ 3 ] = columns;

    for ( variableIndex = 0; variableIndex < variableCount; ++variableIndex ) {
      const char* const variableName = variableNames[ variableIndex ];
      int variableId = -1;
      status = nc_inq_varid( file, variableName, &variableId );

      if ( status != NC_NOERR ) {
        const char* const message = nc_strerror( status );
        failureMessage( "Can't determine id of variable %s because %s.",
                        variableName, message );
        variableIndex = variableCount; /* Stop looping. */
      } else {
        Integer cell = 0;
        Integer timestep = 0;
        CoordinateFunction coordinateFunction =
          coordinateFunctions[ variableIndex ];

        if ( coordinateFunction ) {

#pragma omp parallel for

          for ( cell = 0; cell < cells; ++cell ) {
            const Integer row    = cell / columns;
            const Integer column = cell % columns;
            const Real coordinate = coordinateFunction( grid, row, column );
            data[ cell ] = coordinate; /* 64-bits to 32-bits. */
          }
        }

        for ( timestep = 0; timestep < timesteps; ++timestep ) {
          Integer layer = 0;
          start[ 0 ] = timestep;

          for ( layer = 0; layer < layers; ++layer ) {
            start[ 1 ] = layer;

            if ( ! coordinateFunction ) {
              const Real elevation = grid->elevation( grid, layer );

              for ( cell = 0; cell < cells; ++cell ) {
                data[ cell ] = elevation; /* 64-bits to 32-bits. */
              }
            }

            DEBUG( fprintf( stderr, "calling nc_put_vara_float()...\n" ); )
            status = nc_put_vara_float(file, variableId, start, count, data);
            DEBUG( fprintf( stderr, "...done\n" ); )

            if ( status != NC_NOERR ) {
              const char* const message = nc_strerror( status );
              failureMessage( "Can't write variable %s because %s.",
                              variableName, message );
              layer = layers; /* Stop looping. */
              timestep = timesteps;
            }
          }
        }
      }
    }

    result = status == NC_NOERR;
    FREE( data );
  }

  return result;
}



/******************************************************************************
PURPOSE: writeM3IOData - Write data to the M3IO file.
INPUTS:  Integer file              NetCDF file handle.
         const Name variableName   Name of variable to write.
         Integer timestep          Timestep to write.
         Integer layers            Number of layers  to write.
         Integer rows              Number of rows    to write.
         Integer columns           Number of columns to write.
         const void* data          Data to write.
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeM3IOData( Integer file, const Name variableName,
                       Integer timestep, Integer layers,
                       Integer rows, Integer columns,
                       const void* data ) {

  PRE05( file > -1, variableName, *variableName,
         GT_ZERO3( layers, rows, columns ), data );

  Integer result = 0;
  int id = -1;
  int status = nc_inq_varid( file, variableName, &id );

  if ( status == NC_NOERR ) {
    nc_type type = 0;
    status = nc_inq_vartype( file, id, &type );

    if ( AND2( status == NC_NOERR, ! IN3( type, NC_FLOAT, NC_INT ) ) ) {
      status = NC_EBADTYPE;
    }

    if ( status == NC_NOERR ) {
      const Integer size = layers * rows * columns;
      float* fdata = type == NC_FLOAT ? NEW_ZERO( float, size ) : 0;
      int*   idata = type == NC_INT   ? NEW_ZERO( int,   size ) : 0;

      if ( OR2( fdata, idata ) ) {
        Integer index = 0;

        if ( type == NC_FLOAT ) {
          const Real* const input = data;

#pragma omp parallel for

          for ( index = 0; index < size; ++index ) {
            const Real value = input[ index ];
            const float fvalue = CLAMPED_TO_RANGE( value, -FLT_MAX, FLT_MAX );
            fdata[ index ] = fvalue;
          }

        } else {
          const Integer* const input = data;

#pragma omp parallel for

          for ( index = 0; index < size; ++index ) {
            const Integer value = input[ index ];
            const int ivalue = CLAMPED_TO_RANGE( value, 0, INT_MAX );
            idata[ index ] = ivalue;
          }
        }

        {
          size_t start[ 4 ] = { 0, 0, 0, 0 };
          size_t count[ 4 ] = { 1, 0, 0, 0 };
          start[ 0 ] = timestep;
          count[ 1 ] = layers;
          count[ 2 ] = rows;
          count[ 3 ] = columns;

          DEBUG( fprintf( stderr,
                          "calling nc_put_vara_%s( "
                          "start = [%lu, %lu, %lu, %lu] )\n"
                          "count = [%lu, %lu, %lu, %lu] )\n",
                          type == NC_FLOAT ? "float" : "int",
                          start[ 0 ], start[ 1 ], start[ 2 ], start[ 3 ],
                          count[ 0 ], count[ 1 ], count[ 2 ], count[ 3 ] ); )

          if ( type == NC_FLOAT ) {
            status = nc_put_vara_float( file, id, start, count, fdata );
          } else {
            status = nc_put_vara_int( file, id, start, count, idata );
          }
        }

        DEBUG( fprintf( stderr, "...done\n" ); )

        result = status == NC_NOERR;
        FREE( fdata );
        FREE( idata );
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    variableName, message );
  }

  return result;
}



/******************************************************************************
PURPOSE: copyDataToGrid - Copy sparse data onto a 2D grid.
INPUTS:  Integer points                   Number of sparse data points.
         const Integer rows[ points ]     1-based grid row number per point.
         const Integer columns[ points ]  1-based grid column number per point.
         const Real pointData[ points * gridLayers ]   Data value per point.
         Real scale                       Scale factor for pointData.
         Integer gridLayers               Number of grid layers.
         Integer gridRows                 Number of grid rows.
         Integer gridColumns              Number of grid columns.
OUTPUTS: Real gridData[ gridLayers * gridRows * gridColumns ]
                                          Initialized grid data.
******************************************************************************/

void copyDataToGrid( Integer points, const Integer rows[],
                     const Integer columns[], const Real pointData[],
                     Real scale,
                     Integer gridLayers, Integer gridRows, Integer gridColumns,
                     Real gridData[] ) {

  PRE05( points >= 0,
         GT_ZERO3( gridLayers, gridRows, gridColumns ),
         IMPLIES( points > 0,
                  AND5( minimumItemI( rows, points ) >= 1,
                        maximumItemI( rows, points ) <= gridRows,
                        minimumItemI( columns, points ) >= 1,
                        maximumItemI( columns, points ) <= gridColumns,
                        isNanFree( pointData, points * gridLayers ) ) ),
         ! isNan( scale ),
         gridData );

  /* Initialize all grid cells to 'missing value': */

  const Integer rowsTimesColumns = gridRows * gridColumns;
  const Integer gridPoints = gridLayers * rowsTimesColumns;
  Integer point = 0;

#pragma omp parallel
  { /* Begin parallel region. */

#pragma omp for

  for ( point = 0; point < gridPoints; ++point ) {
    gridData[ point ] = BADVAL3;
  }

  /* Over-write non-missing points: */

#pragma omp for

  for ( point = 0; point < points; ++point ) {
    const Integer row    = rows[    point ] - 1;
    const Integer column = columns[ point ] - 1;
    Integer gridIndex    = row * gridColumns + column;
    Integer pointIndex   = point * gridLayers;
    Integer layer        = 0;

    for ( layer = 0; layer < gridLayers;
          ++layer, ++pointIndex, gridIndex += rowsTimesColumns ) {
      CHECK2( IN_RANGE( gridIndex, 0, gridPoints - 1 ),
              IN_RANGE( pointIndex, 0, points * gridLayers - 1 ) );
      gridData[ gridIndex ] = pointData[ pointIndex ] * scale;
    }
  }

  } /* End parallel region. */

  POST0( isNanFree( gridData, gridLayers * gridRows * gridColumns ) );
}



/******************************************************************************
PURPOSE: copyIntDataToGrid - Copy sparse integer data onto a 2D grid.
INPUTS:  Integer points                   Number of sparse data points.
         const Integer rows[ points ]     1-based grid row number per point.
         const Integer columns[ points ]  1-based grid column number per point.
         const Integer pointData[ points * gridLayers ]   Data value per point.
         Integer gridLayers               Number of grid layers.
         Integer gridRows                 Number of grid rows.
         Integer gridColumns              Number of grid columns.
OUTPUTS: Integer gridData[ gridLayers * gridRows * gridColumns ]
                                          Initialized grid data.
******************************************************************************/

void copyIntDataToGrid( Integer points, const Integer rows[],
                        const Integer columns[], const Integer pointData[],
                        Integer gridLayers, Integer gridRows,
                        Integer gridColumns, Integer gridData[] ) {

  PRE04( points >= 0,
         GT_ZERO3( gridLayers, gridRows, gridColumns ),
         IMPLIES( points > 0,
                  AND4( minimumItemI( rows, points ) >= 1,
                        maximumItemI( rows, points ) <= gridRows,
                        minimumItemI( columns, points ) >= 1,
                        maximumItemI( columns, points ) <= gridColumns ) ),
         gridData );

  /* Initialize all grid cells to 'missing value': */

  const Integer rowsTimesColumns = gridRows * gridColumns;
  const Integer gridPoints = gridLayers * rowsTimesColumns;
  Integer point = 0;

#pragma omp parallel
  { /* Begin parallel region. */

#pragma omp for

    for ( point = 0; point < gridPoints; ++point ) {
      gridData[ point ] = IMISS3;
    }

    /* Over-write non-missing points: */

#pragma omp for

    for ( point = 0; point < points; ++point ) {
      const Integer row    = rows[    point ] - 1;
      const Integer column = columns[ point ] - 1;
      Integer gridIndex    = row * gridColumns + column;
      Integer pointIndex   = point * gridLayers;
      Integer layer        = 0;

      for ( layer = 0; layer < gridLayers;
            ++layer, ++pointIndex, gridIndex += rowsTimesColumns ) {
        CHECK2( IN_RANGE( gridIndex, 0, gridPoints - 1 ),
                IN_RANGE( pointIndex, 0, points * gridLayers - 1 ) );
        gridData[ gridIndex ] = pointData[ pointIndex ];
      }
    }

  } /* End parallel region. */

}



/******************************************************************************
PURPOSE: copyDataToGrid3 - Copy sparse data onto a 3D grid.
INPUTS:  Integer points                   Number of sparse data points.
         const Integer layers[ points ]   1-based grid layer number per point.
         const Integer rows[ points ]     1-based grid row number per point.
         const Integer columns[ points ]  1-based grid column number per point.
         const Real pointData[ points ]   Data value per point.
         Real scale                       Scale factor for pointData.
         Integer gridLayers               Number of grid layers.
         Integer gridRows                 Number of grid rows.
         Integer gridColumns              Number of grid columns.
OUTPUTS: Real gridData[ gridLayers * gridRows * gridColumns ]
                                          Initialized grid data.
******************************************************************************/

void copyDataToGrid3( Integer points,
                      const Integer layers[],
                      const Integer rows[],
                      const Integer columns[],
                      const Real pointData[],
                      Real scale,
                      Integer gridLayers,
                      Integer gridRows,
                      Integer gridColumns,
                      Real gridData[] ) {

  PRE05( points >= 0,
         GT_ZERO3( gridLayers, gridRows, gridColumns ),
         IMPLIES( points > 0,
                  AND7( minimumItemI( layers, points ) >= 1,
                        maximumItemI( layers, points ) <= gridLayers,
                        minimumItemI( rows, points ) >= 1,
                        maximumItemI( rows, points ) <= gridRows,
                        minimumItemI( columns, points ) >= 1,
                        maximumItemI( columns, points ) <= gridColumns,
                        isNanFree( pointData, points ) ) ),
         ! isNan( scale ),
         gridData );

  /* Initialize all grid cells to 'missing value': */

  const Integer rowsTimesColumns = gridRows * gridColumns;
  const Integer gridPoints = gridLayers * rowsTimesColumns;
  Integer point = 0;

#pragma omp parallel
  { /* Begin parallel region. */

#pragma omp for

  for ( point = 0; point < gridPoints; ++point ) {
    gridData[ point ] = BADVAL3;
  }

  /* Over-write non-missing points: */

#pragma omp for

  for ( point = 0; point < points; ++point ) {
    const Integer layer  = layers[  point ] - 1;
    const Integer row    = rows[    point ] - 1;
    const Integer column = columns[ point ] - 1;
    Integer gridIndex    =
      layer * rowsTimesColumns + row * gridColumns + column;
    Integer pointIndex   = point;
    CHECK2( IN_RANGE( gridIndex, 0, gridPoints - 1 ),
            IN_RANGE( pointIndex, 0, points - 1 ) );
    gridData[ gridIndex ] = pointData[ pointIndex ] * scale;
  }

  } /* End parallel region. */

  POST0( isNanFree( gridData, gridLayers * gridRows * gridColumns ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: writeM3IODimensions - Write the dimensions for an M3IO file.
INPUTS:  Integer file          NetCDF file to write to.
         Integer timesteps      Number of timesteps.
         Integer variables      Number of variables in variableNames array.
         Integer layers         Number of layers.
         Integer rows           Number of rows.
         Integer columns        Number of columns.
OUTPUTS: Integer dimensions[ M3IO_DIMS ] Initialized NetCDF dimension handles.
                 dimensions[ TSTEP, DATE_TIME, LAY, VAR, ROW_DIM, COL ]
RETURNS: Integer 1 if successful, else 0.
******************************************************************************/

static Integer writeM3IODimensions( Integer file,
                                    Integer timesteps,
                                    Integer variables,
                                    Integer layers,
                                    Integer rows,
                                    Integer columns,
                                    Integer dimensions[ M3IO_DIMS ] ) {

  PRE03( file > -1,
         GT_ZERO5( timesteps, variables, layers, rows, columns ), dimensions );

  const char* const names[ M3IO_DIMS ] = {
    "TSTEP", "DATE-TIME", "LAY", "VAR", "ROW", "COL"
  };

  Integer result = 0;

  Integer values[ M3IO_DIMS ] = { 0, 2, 0, 0, 0, 0 };
  values[ TSTEP ]     = timesteps;
  values[ LAY ]       = layers;
  values[ VAR ]       = variables;
  values[ ROW_DIM ]   = rows;
  values[ COL ]       = columns;

  result = createDimensions( file, M3IO_DIMS, names, values, dimensions );

  if ( ! result ) {
    dimensions[ TSTEP ]     = -1;
    dimensions[ DATE_TIME ] = -1;
    dimensions[ LAY ]       = -1;
    dimensions[ VAR ]       = -1;
    dimensions[ ROW_DIM ]   = -1;
    dimensions[ COL ]       = -1;
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        GE_ZERO6( dimensions[ TSTEP ], dimensions[ DATE_TIME ],
                                  dimensions[ LAY ], dimensions[ VAR ],
                                  dimensions[ ROW_DIM ], dimensions[ COL ] ),
                        AND6( dimensions[ TSTEP ]     == -1,
                              dimensions[ DATE_TIME ] == -1,
                              dimensions[ LAY ]       == -1,
                              dimensions[ VAR ]       == -1,
                              dimensions[ ROW_DIM ]   == -1,
                              dimensions[ COL ]       == -1 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: writeM3IOVariables - Write variables to a M3IO file.
INPUTS:  Integer file          M3IO file to write to.
         Integer variables     Number of variables in variableNames array.
         const Name variableNames[]  Names of variables.
         const Name variableUnits[]  Units of variables.
         const Integer dimensions[ M3IO_DIMS ]
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer writeM3IOVariables( Integer file, Integer variables,
                                   const Name variableNames[],
                                   const Name variableUnits[],
                                   const Integer dimensions[ M3IO_DIMS ] ) {

  PRE09( file > -1, variables > 0,
         variableNames, *variableNames, **variableNames,
         variableUnits, *variableUnits, **variableUnits,
         dimensions );

  Integer result = 0;
  Integer variable = 0;
  Integer variableId = -1;
  Integer dimensionIds[ M3IO_DIMS ] = { -1, -1, -1, -1, -1, -1 };

  dimensionIds[ 0 ] = dimensions[ TSTEP ];
  dimensionIds[ 1 ] = dimensions[ VAR ];
  dimensionIds[ 2 ] = dimensions[ DATE_TIME ];
  variableId = createVariable( file, "TFLAG", "<YYYYDDD,HHMMSS>", NC_INT,
                               0, 3, dimensionIds );

  result = variableId > -1;
  result = AND2( result, writeTextAttribute( file, variableId, "long_name",
                                             "TFLAG           " ) );
  result = AND2( result, writeTextAttribute( file, variableId, "var_desc",
                                             "Timestep-valid flags:  "
                                             "(1) YYYYDDD or (2) HHMMSS"
                                             "                         "
                                             "       " ) );

  dimensionIds[ 1 ] = dimensions[ LAY ];
  dimensionIds[ 2 ] = dimensions[ ROW_DIM ];
  dimensionIds[ 3 ] = dimensions[ COL ];

  for ( variable = 0; AND2( result, variable < variables ); ++variable ) {
    const nc_type type =
      strcmp( variableNames[ variable ], "COUNT" ) == 0 ? NC_INT : NC_FLOAT;
    char desc[ MXDLEN3 + 1 ] = "";
    char name[ NAMLEN3 + 1 ] = "";
    char unit[ NAMLEN3 + 1 ] = "";
    memset( desc, 0, sizeof desc );
    memset( name, 0, sizeof name );
    memset( unit, 0, sizeof unit );
    strncpy( name, variableNames[ variable ], NAMLEN3 );
    strncpy( unit, variableUnits[ variable ], NAMLEN3 );
    strncpy( desc, variableNames[ variable ], MXDLEN3 );

    variableId = createVariable( file, name, unit, type, 0, 4, dimensionIds );

    expandString( desc, desc, MXDLEN3 );
    expandString( name, name, NAMLEN3 );
    expandString( unit, unit, NAMLEN3 );

    result = variableId > -1;
    result = AND2( result,
                   writeTextAttribute( file, variableId, "long_name", name ) );

    if ( strcmp( name, "LONGITUDE       " ) == 0 ) {
      result = AND2( result, writeTextAttribute( file, variableId,
                                                 "var_desc",
                                                 "Longitude at the center "
                                                 "of each grid cell       "
                                                 "                        "
                                                 "        " ) );
    } else if ( strcmp( name, "LATITUDE        " ) == 0 ) {
      result = AND2( result, writeTextAttribute( file, variableId,
                                                 "var_desc",
                                                 "Latitude at the center "
                                                 "of each grid cell       "
                                                 "                        "
                                                 "         " ) );
    } else if ( strcmp( name, "COUNT           " ) == 0 ) {
      result = AND2( result, writeTextAttribute( file, variableId,
                                                 "var_desc",
                                                 "Number of data points "
                                                 "regridded into grid cell "
                                                 "                        "
                                                 "         " ) );
    } else {
      result =
        AND2( result, writeTextAttribute( file, variableId, "var_desc", desc));
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeM3IOAttributes - Write the global attributes of an M3IO file.
INPUTS:  Integer file                NetCDF file to write to.
         Integer hoursPerTimestep    Number of hours per timestep: 1, 24, etc.
         Integer firstTimestamp      First timestamp. YYYYDDDHHMMSS.
         Integer variables           Number of variables in variableNames array.
         Integer layers              Number of grid layers to use.
         const Name variableNames[]  Names of variables.
         const char* description     Description of the data.
         const Grid* grid            Grid description object.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer writeM3IOAttributes( Integer file,
                                    Integer hoursPerTimestep,
                                    Integer firstTimestamp,
                                    Integer variables,
                                    Integer layers,
                                    const Name variableNames[],
                                    const char* description,
                                    const Grid* grid ) {

  PRE010( file > -1,
          variables > 0,
          layers > 0,
          hoursPerTimestep > 0,
          isValidTimestamp( firstTimestamp ),
          variableNames, *variableNames, **variableNames,
          grid, grid->invariant( grid ) );

  const Integer rows    = grid->rows( grid );
  const Integer columns = grid->columns( grid );
  const Projector* const projector = grid->projector( grid );
  const Real xcent = projector ? projector->centralLongitude( projector) : 0.0;
  const Real ycent = projector ? projector->centralLatitude( projector ) : 0.0;
  const char* const name = projector ? projector->name( projector ) : "LonLat";
  const Integer gdtyp =
    ! strcmp( name, "Stereographic" ) ? POLGRD3
    : ! strcmp( name, "Lambert" ) ? LAMGRD3
    : ! strcmp( name, "Mercator" ) ? EQMGRD3
    : LATGRD3;

  const Lambert* const lambert = (const Lambert*)
    ( ! strcmp( name, "Lambert" ) ? grid->projector( grid ) : 0 );

  const Stereographic* const stereographic = (const Stereographic*)
    ( ! strcmp( name, "Stereographic" ) ? grid->projector( grid ) : 0 );

  const Mercator* const mercator = (const Mercator*)
    ( ! strcmp( name, "Mercator" ) ? grid->projector( grid ) : 0 );

  const Real p_alp =
    lambert ? lambert->lowerLatitude( lambert )
    : stereographic ? SIGN( ycent )
    : mercator ? mercator->centralLongitude( mercator )
    : 0.0;

  const Real p_bet =
    lambert ? lambert->upperLatitude( lambert )
    : stereographic ? stereographic->secantLatitude( stereographic )
    : 0.0;
  
  const Real p_gam = xcent;

  const Real xorig = grid->westEdge( grid );
  const Real yorig = grid->southEdge( grid );
  const Real xcell = grid->cellWidth( grid );
  const Real ycell = grid->cellHeight( grid );
  Integer result = 0;
  const Integer yyyydddhhmm = nowUTC();
  const Integer cdate = yyyydddhhmm  / 10000; /* YYYYDDD. */
  const Integer ctime = yyyydddhhmm  % 10000; /* HHMMSS.  */
  const Integer sdate = firstTimestamp / 10000; /* YYYYDDD. */
  const Integer stime = firstTimestamp % 10000; /* HHMMSS.  */
  const Integer tstep = 10000 * hoursPerTimestep; /* HHMMSS.  */
  const Real levels[ 2 ] = { 1.0, 0.995 };
  const char* const version = "1.0 1997349 (Dec. 15, 1997)";
  const char* const exec_id = "????????????????"
    "                                                                ";
  enum { FILE_DESCRIPTION_LENGTH = MXDLEN3 * MXDESC3 };
  char fileDescription[ FILE_DESCRIPTION_LENGTH + 1 ] = "";
  memset( fileDescription, ' ', sizeof fileDescription );
  strncpy( fileDescription, description, FILE_DESCRIPTION_LENGTH );
  expandString( fileDescription, fileDescription, FILE_DESCRIPTION_LENGTH );
  CHECK( strlen( fileDescription ) == FILE_DESCRIPTION_LENGTH );
  CHECK( XOR4( lambert, stereographic, mercator, projector == 0 ) );

  result = writeTextAttribute( file, NC_GLOBAL, "IOAPI_VERSION", version );

  result = AND2( result,
                 writeTextAttribute( file, NC_GLOBAL, "EXEC_ID", exec_id ) );

  result = AND2( result, writeIntegerAttribute( file, "FTYPE", 1 ) );
  result = AND2( result, writeIntegerAttribute( file, "CDATE", cdate ) );
  result = AND2( result, writeIntegerAttribute( file, "CTIME", ctime ) );
  result = AND2( result, writeIntegerAttribute( file, "WDATE", cdate ) );
  result = AND2( result, writeIntegerAttribute( file, "WTIME", ctime ) );
  result = AND2( result, writeIntegerAttribute( file, "SDATE", sdate ) );
  result = AND2( result, writeIntegerAttribute( file, "STIME", stime ) );
  result = AND2( result, writeIntegerAttribute( file, "TSTEP", tstep ) );
  result = AND2( result, writeIntegerAttribute( file, "NTHIK", 1 ) );
  result = AND2( result, writeIntegerAttribute( file, "NCOLS", columns ) );
  result = AND2( result, writeIntegerAttribute( file, "NROWS", rows ) );
  result = AND2( result, writeIntegerAttribute( file, "NLAYS", layers ) );
  result = AND2( result, writeIntegerAttribute( file, "NVARS", variables ) );
  result = AND2( result, writeIntegerAttribute( file, "GDTYP", gdtyp ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "P_ALP", p_alp ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "P_BET", p_bet ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "P_GAM", p_gam ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "XCENT", xcent ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "YCENT", ycent ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "XORIG", xorig ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "YORIG", yorig ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "XCELL", xcell ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_DOUBLE,
                                             "YCELL", ycell ) );
  result = AND2( result, writeIntegerAttribute( file, "VGTYP", VGSGPN3 ) );
  result = AND2( result, writeRealAttribute( file, NC_GLOBAL, NC_FLOAT,
                                             "VGTOP",10000.0 ) );

  if ( layers > 1 ) {
    const Integer gridLevels = layers + 1;
    Real copyLevels[ MXLAYS3 ];
    Integer level = 0;

    do {
      copyLevels[ level ] = grid->level( grid, level );
      ++level;
    } while ( level < gridLevels );

    result = AND2( result,
                   writeRealArrayAttribute( file, NC_FLOAT, "VGLVLS",
                                            copyLevels, layers + 1 ) );
  } else {
    result = AND2( result, writeRealArrayAttribute( file, NC_FLOAT, "VGLVLS",
                                                    levels, 2 ) );
  }

  result = AND2( result, writeTextAttribute( file, NC_GLOBAL, "GDNAM",
                                             "M_02_99BRACE    " ) );

  result = AND2( result, writeTextAttribute( file, NC_GLOBAL, "UPNAM",
                                             "XDRConvert      " ) );

  result = AND2( result, writeVarListAttribute( file, variableNames,
                                                variables ) );

  result = AND2( result, writeTextAttribute( file, NC_GLOBAL, "FILEDESC",
                                             fileDescription ) );

  result = AND2( result, writeTextAttribute( file, NC_GLOBAL, "HISTORY",
                                             "XDRConvert" ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeVarListAttribute - Write the value of the VAR-LIST attribute.
INPUTS:  Integer file                NetCDF file ID.
         const Name variableNames[]  Names of variables.
         Integer variables           Number of subset variables.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

static Integer writeVarListAttribute( Integer file,
                                      const Name variableNames[],
                                      Integer variables ) {

  PRE04( file >= 0, variableNames, *variableNames, variables > 0 );

  Integer result = 0;
  int status = 0;
  char attribute[ 10 * MXVARS3 * ( NC_MAX_NAME + 1 ) ] = "";
  Integer variable = 0;
  Integer length = 0;

  memset( attribute, 0, sizeof attribute );

  do {
    const char* const variableName = variableNames[ variable ];
    CHECK( length < sizeof attribute / sizeof *attribute );
    expandString( attribute + length, variableName, NAMLEN3 );
    length += NAMLEN3;
    ++variable;
  } while ( variable < variables );

  status = nc_put_att_text( file, NC_GLOBAL, "VAR-LIST",
                            strlen( attribute ), attribute );

  result = status == NC_NOERR;

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage("Can't write text attribute VAR-LIST because %s.", message);
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeM3IOTFLAG - Write the TFLAG variable to a subset file.
INPUTS:  Integer file           NetCDF file to write to.
         Integer variables      Number of variables in variableNames array.
         Integer timesteps      Number of timesteps.
         Integer hoursPerTimestep  Number of hours per timestep: 1, 24, etc.
         Integer firstTimestep  First timestep.
RETURNS: Integer 1 if successful, else 0 and failureMessage() is called.
******************************************************************************/

static Integer writeM3IOTFLAG( Integer file,
                               Integer variables,
                               Integer timesteps,
                               Integer hoursPerTimestep,
                               Integer firstTimestamp ) {

  PRE03( file > -1,
         GT_ZERO3( variables, timesteps, hoursPerTimestep ),
         isValidTimestamp( firstTimestamp ) );

  Integer result = 0;
  int id = -1;
  int status = nc_inq_varid( file, "TFLAG", &id );

  if ( OR2( status != NC_NOERR, id < 0 ) ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable TFLAG because %s.",
                    message );
  } else {
    const Integer count = timesteps * variables * 2;
    int* data = NEW_ZERO( int, count );

    if ( data ) {
      const Integer tstep = 10000 * hoursPerTimestep;
      Integer yyyyddd     = firstTimestamp / 10000;
      Integer hhmmss      = firstTimestamp % 10000 * 100;
      int* dateTime = data;
      Integer timestep = 0;

      /*
       * For each timestep, advance the timestamp
       * and replicate it for each subset variable:
       */

      do {
        Integer variable = 0;
        CHECK2( isSignedInt( yyyyddd ), isSignedInt( hhmmss ) );

        do {
          *dateTime++ = yyyyddd;
          *dateTime++ = hhmmss;
          ++variable;
        } while ( variable < variables );

        incrementTime( &yyyyddd, &hhmmss, tstep );
        ++timestep;
      } while ( timestep < timesteps );

      {
        size_t starts[ 3 ] = { 0, 0, 0 };
        size_t counts[ 3 ] = { 0, 0, 2 };
        counts[ 0 ] = timesteps;
        counts[ 1 ] = variables;
        status = nc_put_vara_int( file, id, starts, counts, data );
        result = status == NC_NOERR;
      }

      if ( status != NC_NOERR ) {
        const char* const message = nc_strerror( status );
        failureMessage( "Can't write TFLAG variable because %s.", message );
      }

      FREE_ZERO( data );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}




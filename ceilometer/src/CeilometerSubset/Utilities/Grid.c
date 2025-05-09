/******************************************************************************
PURPOSE: Grid.c - Define ADT for grids supporting projection.

NOTES:   See Grid.h.

HISTORY: 2007/12 plessel.todd@epa.gov, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <string.h> /* For memset(), memcpy(). */
#include <math.h>   /* For log(), exp(), fabs(). */

#ifdef DEBUGGING
#include <stdio.h> /* For fprintf(). */
#endif

#if defined( _OPENMP ) && ! defined( SERIAL_REGRID )
#include <omp.h> /* For omp_*_lock(). */
#else
typedef long long omp_lock_t;
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#define omp_init_lock(unused)
#define omp_set_lock(unused)
#define omp_unset_lock(unused)
#define omp_destroy_lock(unused)
#endif

#include <Assertions.h>    /* For PRE*(), etc. */
#include <BasicNumerics.h> /* For Integer, Real, isNan(). */
#include <Failure.h>       /* For failureMessage(). */
#include <Memory.h>        /* For NEW_ZERO(), FREE_ZERO(). */
#include <Projector.h>     /* For validLongitudesAndLatitudes(). */
#include <Lambert.h>       /* For Lambert. */
#include <Mercator.h>      /* For Mercator. */
#include <Stereographic.h> /* For Stereographic. */
#include <Utilities.h>     /* For appendNote(). */
#include <RegridQuadrilaterals.h> /* For regridQuadrilaterals(). */
#include <Grid.h>          /* For public interface. */

#ifndef MIN
#define MIN( a, b ) ( (a) < (b) ? (a) : (b) )
#endif

#ifndef MAX
#define MAX( a, b ) ( (a) > (b) ? (a) : (b) )
#endif

/*================================== TYPES ==================================*/

/* Grid cell used for aggregate: */

typedef struct {
  Integer count;             /* Number of points in this grid cell. */
  Integer column;            /* 1-based column number of grid cell. */
  Integer row;               /* 1-based row    number of grid cell. */
  Integer layer;             /* 1-based layer  number of grid cell. */
  Real    longitude;         /* Longituyde of grid cell center. */
  Real    latitude;          /* Latitude   of grid cell center. */
  Real    elevation;         /* Elevation  of grid cell center in MAMSL. */
  Real    surfaceElevation;  /* Elevation  of surface in MAMSL. */
  Real    radius;            /* Normalized squared-distance from cell center*/
  Real    data;              /* Aggregated data. */
  Real    data2;             /* Aggregated data of 2nd vector component. */
  Real    weights;           /* Data weight sum. */
  Real    minimumValidValue; /* Minimum value of valid data. */
  RegriddedNote regriddedNote; /* Optional: appended notes. */
  omp_lock_t lock;           /* Lock to prevent race-condition on updates. */
} Cell;

static void lockCell(   Cell* cell ) { omp_set_lock(   &cell->lock ); }
static void unlockCell( Cell* cell ) { omp_unset_lock( &cell->lock ); }

struct GridPrivate {
  Projector* projector;  /* Projector for (lon, lat) -> (x, y).    */

  /* Input 2D grid parameters: */

  Integer columns;     /* Number of columns of grid cells.         */
  Integer rows;        /* Number of rows    of grid cells.         */
  Real    xMinimum;    /* Projected west  edge of grid, in meters. */
  Real    yMinimum;    /* Projected south edge of grid, in meters. */
  Real    cellWidth;   /* Width  of each grid cell,     in meters. */
  Real    cellHeight;  /* Height of each grid cell,     in meters. */

  /* Derived 2D grid parameters: */

  Real xMaximum;       /* xMinimum + cellWidth  * columns.         */
  Real yMaximum;       /* yMinimum + cellHeight * rows.            */
  Real oneOverWidth;   /* 1 / cellWidth.                           */
  Real oneOverHeight;  /* 1 / cellHeight.                          */
  Real* longitudes;    /* longitudes[ rows ][ columns ] of grid cell centers*/
  Real* latitudes;     /* latitudes[  rows ][ columns ] of grid cell centers*/

  /* Input 3D grid parameters: */

  Integer layers;      /* Or 1 if 2D grid only.                         */
  Integer type;        /* Vertical grid type: e.g., VGSGPN3.            */
  Real    topPressure; /* Pressure, in pascals, at top of grid.         */
  Real*   levels;      /* levels[ layers + 1 ]. E.g., sigma-pressures   */
  Real    g;           /* Gravitational force, e.g., 9.81 m/s^2.        */
  Real    R;           /* Gas constant e.g., 287.04 J/kg/K = m^3/s/K.   */
  Real    A;           /* Atmospheric lapse rate, e.g., 50.0 K/kg.      */
  Real    T0s;         /* Reference surface temperature, e.g., 290.0 K. */
  Real    P00;         /* Reference surface pressure, e.g., 100000 P.   */

  /* Derived 3D grid parameters: */

  Real*   z;           /* z[thread][layers + 1] meters above mean sea-level.*/
  Cell*   cells;       /* cells[ rows ][ columns ][ layers ].           */
};

typedef void (*PreAggregator)( Integer column, Integer row,
                               Real gridLongitude, Real gridLatitude,
                               Real xOffset, Real yOffset, Real zOffset,
                               Real inputData, Real inputData2, Cell* cell );

typedef void (*Aggregator)( Real xOffset, Real yOffset, Real zOffset,
                            Real inputData, Real inputData2, Cell* cell );

typedef void (*PostAggregator)( Cell* cell );

typedef struct {
  PreAggregator  preAggregator;  /* Called to initilize cell. */
  Aggregator     aggregator;     /* Called to aggregate cell. */
  PostAggregator postAggregator; /* Called after aggregating matched cells. */
} AggregatorEntry;


/*========================== FORWARD DECLARATIONS ===========================*/

/* Commands: */

static void free__( Grid* self );

static void projectXY( Grid* self, Integer count,
                       const Real longitudes[], const Real latitudes[],
                       Integer* griddedPoints,
                       Integer columns[], Integer rows[],
                       Real xCenterOffsets[], Real yCenterOffsets[],
                       Real gridLongitudes[], Real gridLatitudes[] );

static void projectZ( Grid* self, Integer count,
                      const Real elevations[],
                      Integer* griddedPoints,
                      Integer layers[], Real centerOffsets[],
                      Real gridElevations[] );

static void aggregate( Grid* grid,
                       Integer method, Real minimumValidValue,
                       Integer inputPoints,
                       Integer columns[], Integer rows[],
                       const Real xOffsets[], const Real yOffsets[],
                       Real gridLongitudes[], Real gridLatitudes[],
                       Integer layers, const Real elevations[],
                       const Real inputData[], const Real inputData2[],
                       Integer* outputPoints,
                       Real outputData[], Real outputData2[],
                       Real gridElevations[] );

static void regrid( Grid* grid,
                    Integer method, Real minimumValidValue,
                    Integer points, Integer levels,
                    const Real longitudes[], const Real latitudes[],
                    const Real elevations[],
                    const Real inputData[], const Real inputData2[],
                    const Note notes[],
                    Integer* outputPoints,
                    Integer columns[], Integer rows[], Integer layers[],
                    Real gridLongitudes[], Real gridLatitudes[],
                    Real gridElevations[],
                    Real outputData[], Real outputData2[],
                    RegriddedNote regriddedNotes[] );

static void regridSwath( Grid* self,
                         const Integer method,
                         const Real minimumValidValue,
                         const Integer points,
                         const Real longitudesSW[],
                         const Real longitudesSE[],
                         const Real longitudesNW[],
                         const Real longitudesNE[],
                         const Real latitudesSW[],
                         const Real latitudesSE[],
                         const Real latitudesNW[],
                         const Real latitudesNE[],
                         const Real inputData[],
                         Integer* regriddedPoints,
                         Integer gridColumns[],
                         Integer gridRows[],
                         Real gridLongitudes[],
                         Real gridLatitudes[],
                         Real outputData[] );

/* Queries: */

static Integer invariant( const Grid* self );

static Integer equal( const Grid* self, const Grid* other );

static Grid* clone( const Grid* self );

static Projector* projector( const Grid* self );

static Integer layers( const Grid* self );

static Integer rows( const Grid* self );

static Integer columns( const Grid* self );

static Real longitude( const Grid* self, Integer row, Integer column );

static Real latitude( const Grid* self, Integer row, Integer column );

static Real elevation( const Grid* self, Integer layer );

static Real level( const Grid* self, Integer level );

static Real westEdge( const Grid* self );

static Real southEdge( const Grid* self );

static Real cellWidth( const Grid* self );

static Real cellHeight( const Grid* self );

/* Helpers: */

static void assignMembers( Grid* self );

static Integer hasMembers( const Grid* self );

static void computeLongitudesAndLatitudes( GridPrivate* data );

static void computeZ( Real g, Real R, Real A, Real T0s, Real P00,
                      Integer layers, Integer type, Real topPressure,
                      const Real levels[], Real z[] );

static Real pressureAtSigmaLevel( Real sigmaLevel, Real pressureAtTop );

static Real heightAtPressure( Real pressure );

static void parseLambert( int argc, char* argv[],
                          Integer* argument,
                          Real* lowerLatitude,    Real* upperLatitude,
                          Real* centralLongitude, Real* centralLatitude,
                          Integer* ok );

static void parseMercator( int argc, char* argv[],
                           Integer* argument,
                           Real* centralLongitude,
                           Integer* ok );

static void parseStereographic( int argc, char* argv[],
                                Integer* argument,
                                Real* centralLongitude, Real* centralLatitude,
                                Real* secantLatitude,
                                Integer* ok );

static void parseCentralLongitudeAndLatitude( int argc, char* argv[],
                                              Integer* argument,
                                              Real* centralLongitude,
                                              Real* centralLatitude,
                                              Integer* ok );

static Real* parseLayers( int argc, char* argv[], Integer* argument,
                          Real* g, Real* R, Real* A, Real* T0s, Real* P00,
                          Integer* layers, Integer* type, Real* topPressure );

static void initializeCells(Integer count, Real minimumValidValue, Cell* cell);

static void finalizeCells( Integer count, Cell cells[] );

static void zeroUnused( Integer index, Integer count,
                        Integer columns[], Integer rows[],
                        Real longitudes[], Real latitudes[] );


static Real radiusSquared( Real x, Real y, Real z );

static void commonPreAggregator( Integer column, Integer row,
                                 Real gridLongitude, Real gridLatitude,
                                 Real unused_xOffset, Real unused_yOffset,
                                 Real unused_zOffset,
                                 Real inputData, Real inputData2, Cell* cell );

static void nearestPreAggregator( Integer column, Integer row,
                                  Real gridLongitude, Real gridLatitude,
                                  Real xOffset, Real yOffset, Real zOffset,
                                  Real inputData, Real inputData2, Cell* cell );

static void nearestAggregator( Real xOffset, Real yOffset, Real zOffset,
                               Real inputData, Real inputData2, Cell* cell );

static void meanAggregator( Real unused_xOffset, Real unused_yOffset,
                            Real unused_zOffset, Real inputData,
                            Real inputData2, Cell* cell );

static void weightedPreAggregator( Integer column, Integer row,
                                   Real gridLongitude, Real gridLatitude,
                                   Real xOffset, Real yOffset, Real zOffset,
                                   Real inputData, Real inputData2, Cell* cell);


static void weightedAggregator( Real xOffset, Real yOffset, Real zOffset,
                                Real inputData, Real inputData2, Cell* cell );

static void weightedPostAggregator( Cell* cell );

static void postNull( Cell* unused_cell ) { }

static void aggregateCellData( Grid* self,
                               PreAggregator preAggregator,
                               Aggregator aggregator,
                               Integer column, Integer row,
                               Real gridLongitude, Real gridLatitude,
                               Real xOffset, Real yOffset,
                               Integer point,
                               const Real data[], const Real data2[],
                               Integer layers, const Real elevations[],
                               const Note note,
                               Cell cells[] );

static Real surfacePointData( Integer point, const Real data[],
                              const Real data2[],
                              Integer layers, const Real elevations[],
                              Real* dataValue2 );

static Integer surfaceIndex( Integer layers, const Real elevations[] );

static Integer binElevation( Real dataElevation, Integer gridLevels,
                             const Real cellElevations[],
                             Integer startingLayer, Real* zOffset );

static void elevationsAtSigmaPressures( Real g, Real R, Real A, Real T0s,
                                        Real P00,
                                        Real surfaceElevation,
                                        Integer levels, Real topPressure,
                                        const Real sigmaPressures[],
                                        Real elevations[] );

/*============================= GLOBAL VARIABLES ============================*/

static const AggregatorEntry aggregators[] = {
  { nearestPreAggregator,  nearestAggregator,  postNull },
  { commonPreAggregator,   meanAggregator,     postNull },
  { weightedPreAggregator, weightedAggregator, weightedPostAggregator }
};

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: newGrid - Construct a Grid.
INPUTS:  Projector* projector  Owned projector for (lon, lat) -> (x, y) or 0.
                               Object ownership is transfered by this call.
         Integer columns       Number of grid cell columns.
         Integer rows          Number of grid cell rows.
         Real    westEdge      Projected west  edge of grid, in meters.
         Real    southEdge     Projected south edge of grid, in meters.
         Real    cellWidth     Width  of each grid cell,     in meters.
         Real    cellHeight    Height of each grid cell,     in meters.
         Integer layers        Number of grid cell layers, or 0 if 2D.
         Integer type          Vertical grid type: e.g., VGSGPN3
         Real    topPressure   Pressure, in pascals, at top of grid.
         const Real levels[]   E.g., sigma-pressure levels at each grid level
                               that bounds a grid cell center.
                               levels[ layers + 1 ] = 0.999,0.995,0.99...
         Real g                Gravitational force, e.g., 9.81 m/s^2.
         Real R                Gas constant e.g., 287.04 J/kg/K = m^3/s/K.
         Real A                Atmospheric lapse rate, e.g., 50.0 K/kg.
         Real T0s              Reference surface temperature, e.g., 290.0 K.
         Real P00              Reference surface pressure, e.g., 100000 P.

RETURNS: Grid* else 0 if failed and failureMessage() called.
******************************************************************************/

Grid* newGrid( Projector* projector,
               Integer columns, Integer rows,
               Real westEdge, Real southEdge,
               Real cellWidth, Real cellHeight,
               Integer layers,
               Integer type, Real topPressure, const Real levels[],
               Real g, Real R, Real A, Real T0s, Real P00 ) {

  PRE017( columns > 0,
          rows > 0,
          rows < INTEGER_MAX / columns,
          ! isNan( westEdge ),
          ! isNan( southEdge ),
          ! isNan( cellWidth ),
          ! isNan( cellHeight ),
          cellWidth > 0.0,
          cellHeight > 0.0,
          IMPLIES_ELSE( projector,
                        projector->invariant( projector ),
                        AND2( isValidLongitudeLatitude( westEdge, southEdge ),
                              isValidLongitudeLatitude(
                                westEdge + columns * cellWidth,
                                southEdge + rows * cellHeight ) ) ),
          ! isNan( g ),
          ! isNan( R ),
          ! isNan( A ),
          ! isNan( T0s ),
          ! isNan( P00 ),
          layers >= 0,
          IMPLIES_ELSE( layers == 0,
                        IS_ZERO8( type, topPressure, levels, g, R, A,T0s,P00),
            AND7( IS_VALID_VERTICAL_GRID_TYPE( type ),
                  ! isNan( topPressure ),
                  topPressure > 0.0,
                  levels,
                  isNanFree( levels, layers + 1 ),
                  IMPLIES( ! IN6( type, 0, IMISS3, VGPRES3, VGZVAL3, VGHVAL3 ),
                           AND2( minimumItem( levels, layers + 1 ) >= 0.0,
                                 maximumItem( levels, layers + 1 ) <= 1.0 ) ),
                  GT_ZERO5( g, R, A, T0s, P00 ) ) ) );

  const Integer threads = omp_get_max_threads();
  const Integer layerPoints = columns * rows;
  Real* longitudes  = NEW_ZERO( Real, layerPoints * 2 );
  Real* latitudes   = longitudes ? longitudes + layerPoints : 0;
  const Integer zCount = ( threads + 1 ) * ( MAX( 1, layers ) + 1 );
  Real* z = longitudes ? NEW_ZERO( Real, zCount ) : 0;
  Cell* cells = z ? NEW_ZERO( Cell, rows * columns * MAX( 1, layers ) ) : 0;
  GridPrivate* data = cells ? NEW_ZERO( GridPrivate, 1 ) : 0;
  Grid* result      = data ? NEW_ZERO( Grid, 1 ) : 0;

  if ( result ) {

    /* Initialize input 2D grid attributes: */

    data->columns    = columns;
    data->rows       = rows;
    data->xMinimum   = westEdge;
    data->yMinimum   = southEdge;
    data->cellWidth  = cellWidth;
    data->cellHeight = cellHeight;
    data->projector  = projector;
    projector = 0; /* Transfered ownership. */

    /* Initialize derived 2D grid attributes: */

    data->xMaximum      = westEdge  + cellWidth  * columns;
    data->yMaximum      = southEdge + cellHeight * rows;
    data->oneOverWidth  = 1.0 / cellWidth;
    data->oneOverHeight = 1.0 / cellHeight;
    data->longitudes    = longitudes;
    data->latitudes     = latitudes;
    longitudes = latitudes = 0; /* Transfered ownership. */
    computeLongitudesAndLatitudes( data );

    /* Initialize input 3D grid attributes: */

    if ( layers == 0 ) {
      data->layers      = 1;
      data->type        = VGSGPN3;
      data->topPressure = 10000.0; /* Pascals. */
      data->g           = 9.81;
      data->R           = 287.04;
      data->A           = 50.0;
      data->T0s         = 290.0;
      data->P00         = 100000.0;
    } else {
      data->type        = type;
      data->layers      = layers;
      data->topPressure = topPressure;
      data->g           = g;
      data->R           = R;
      data->A           = A;
      data->T0s         = T0s;
      data->P00         = P00;
    }

    /* Initialize 3D grid attributes: */

    data->z = z;
    z = 0; /* Transfered ownership. */
    data->cells = cells;
    cells = 0; /* Transfered ownership. */
    data->levels = data->z + threads * ( data->layers + 1 );

    if ( levels ) {
      memcpy( data->levels, levels,
              ( data->layers + 1 ) * sizeof data->levels[ 0 ] );
    } else {
      data->type = VGSGPN3;
      data->levels[ 0 ] = 1.0;
      data->levels[ 1 ] = 0.995;
    }

    computeZ( data->g, data->R, data->A, data->T0s, data->P00,
              data->layers, data->type, data->topPressure, data->levels,
              data->z );

    result->data = data;
    data = 0; /* Transfered ownership. */
    assignMembers( result );
  }

  if ( ! result ) {
    FREE_OBJECT( projector );
  }

  FREE( longitudes );
  latitudes = 0;
  FREE( z );
  FREE( cells );
  FREE( data );

  POST02( projector == 0, IMPLIES( result, result->invariant( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: newSubsetGrid - Construct a subset Grid.
INPUTS:  const Grid* grid  Object to subset.
         const Integer firstLayer   First layer  index (0-based).
         const Integer lastLayer    Last  layer  index (0-based).
         const Integer firstRow     First row    index (0-based).
         const Integer lastRow      Last  row    index (0-based).
         const Integer firstColumn  First column index (0-based).
         const Integer lastColumn   Last  column index (0-based).
RETURNS: Grid* else 0 if failed and failureMessage() called.
******************************************************************************/

Grid* newSubsetGrid( const Grid* grid,
                     const Integer firstLayer, const Integer lastLayer,
                     const Integer firstRow, const Integer lastRow,
                     const Integer firstColumn, const Integer lastColumn ) {

  PRE08( grid, grid->invariant( grid ),
         IN_RANGE( firstLayer, 0, grid->layers( grid ) - 1 ),
         IN_RANGE( lastLayer, firstLayer, grid->layers( grid ) - 1 ),
         IN_RANGE( firstRow, 0, grid->rows( grid ) - 1 ),
         IN_RANGE( lastRow, firstRow, grid->rows( grid ) - 1 ),
         IN_RANGE( firstColumn, 0, grid->columns( grid ) - 1 ),
         IN_RANGE( lastColumn, firstRow, grid->columns( grid ) - 1 ) );

  const Integer layers  = 1 + lastLayer  - firstLayer;
  const Integer rows    = 1 + lastRow    - firstRow;
  const Integer columns = 1 + lastColumn - firstColumn;

  const GridPrivate* const data = grid->data;
  const Projector* const projector = grid->projector( grid );
  Grid* const result =
    newGrid( projector ? projector->clone( projector ) : 0,
             columns, rows,
             westEdge( grid ) + cellWidth( grid ) * ( columns - 1 ),
             southEdge( grid ) + cellHeight( grid ) * ( rows - 1 ),
             cellWidth( grid ), cellHeight( grid ),
             layers, data->type, data->topPressure, data->levels + firstLayer,
             data->g, data->R, data->A, data->T0s, data->P00 );

  POST0( IMPLIES( result, result->invariant( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: projectXY - Project longitude-latitude points onto a regular 2D grid
         in a projected cartographic space and obtain the 1-based grid cell
         column and row numbers where the points lie and also the x and y
         normalized [-1, 1] offsets from the cell centers.
INPUTS:  Grid* self              The Grid.
         Integer count           Number of longitudes/latitudes to test.
         const Real longitudes[] Longitudes to project and test.
         const Real latitudes[]  Latitudes  to project and test.
                                 These input latitudes are assumed to be on a
                                 GRS80/WGS84/NAD83 spheroid/datum.
OUTPUTS: Integer* griddedPoints  Number of points projected onto grid.
         Integer columns[]       1-based grid column number or 0 if outside.
         Integer rows[]          1-based grid row    number or 0 if outside.
         Real xCenterOffsets[]   [-1, 1] offset from center of grid cell.
         Real yCenterOffsets[]   [-1, 1] offset from center of grid cell.
         Real gridLongitudes[]   Optional longitudes of regridded points or 0.
         Real gridLatitudes[]    Optional latitudes  of regridded points or 0.
******************************************************************************/

static void projectXY( Grid* self, Integer count,
                       const Real longitudes[], const Real latitudes[],
                       Integer* griddedPoints,
                       Integer columns[], Integer rows[],
                       Real xCenterOffsets[], Real yCenterOffsets[],
                       Real gridLongitudes[], Real gridLatitudes[] ) {

  PRE10( count > 0, longitudes, latitudes, griddedPoints, columns, rows,
         xCenterOffsets, yCenterOffsets,
         validLongitudesAndLatitudes( count, longitudes, latitudes ),
         IMPLIES_ELSE( gridLongitudes, gridLatitudes, ! gridLatitudes ) );

  const GridPrivate* const data = self->data;
  const Real xMinimum        = data->xMinimum;
  const Real xMaximum        = data->xMaximum;
  const Real yMinimum        = data->yMinimum;
  const Real yMaximum        = data->yMaximum;
  const Real oneOverWidth    = data->oneOverWidth;
  const Real oneOverHeight   = data->oneOverHeight;
  const Integer gridColumns  = data->columns;
  const Integer gridRows     = data->rows;
  const Real* const gridCellCenterLongitudes = data->longitudes;
  const Real* const gridCellCenterLatitudes  = data->latitudes;
  Projector* const projector = data->projector;
  Real (*latitudeAdjuster)( Real ) = 0;
  Integer index = 0;
  Integer griddedPointCount = 0;

  /*
   * Input latitudes are on a WGS84 spheroid.
   * If the grid uses a sphere (like CMAQ) then the input latitudes must be
   * adjusted to be on the sphere before projecting for regridding:
   */

  if ( projector ) {
    Real majorSemiaxis = 0.0;
    Real minorSemiaxis = 0.0;
    projector->ellipsoid( projector, &majorSemiaxis, &minorSemiaxis );
    latitudeAdjuster = majorSemiaxis == minorSemiaxis ? latitudeSphere : 0;
  } else {
    latitudeAdjuster = latitudeSphere;
  }

  *griddedPoints = 0;
  memset( columns,        0, count * sizeof *columns );
  memset( rows,           0, count * sizeof *rows );
  memset( xCenterOffsets, 0, count * sizeof *xCenterOffsets );
  memset( yCenterOffsets, 0, count * sizeof *yCenterOffsets );

  if ( gridLongitudes ) {
    memset( gridLongitudes, 0, count * sizeof *gridLongitudes );
    memset( gridLatitudes,  0, count * sizeof *gridLatitudes );
  }

#pragma omp parallel for reduction( + : griddedPointCount )

  for ( index = 0; index < count; ++index ) {
    const Real longitude = longitudes[ index ];
    const Real latitudeOnWGS84Spheroid = latitudes[ index ];
    const Real latitude =
      latitudeAdjuster ? latitudeAdjuster( latitudeOnWGS84Spheroid )
      : latitudeOnWGS84Spheroid;
    Real x = longitude;
    Real y = latitude;

    if ( projector ) {
      projector->project( projector, longitude, latitude, &x, &y );
    }

    DEBUG( fprintf( stderr, "(%lf %lf)->(%lf, %lf)\n",
                    longitude, latitude, x, y ); )

    if ( AND2( IN_RANGE( x, xMinimum, xMaximum ),
               IN_RANGE( y, yMinimum, yMaximum ) ) ) {
      const Real fractionalColumn = (x - xMinimum) * oneOverWidth  + 1.0;
      const Real fractionalRow    = (y - yMinimum) * oneOverHeight + 1.0;
      Integer column     = fractionalColumn; /* Truncate fraction. */
      Integer row        = fractionalRow;
      Real xCenterOffset = fractionalColumn - column - 0.5;
      Real yCenterOffset = fractionalRow    - row    - 0.5;
      xCenterOffset += xCenterOffset;
      yCenterOffset += yCenterOffset;

      if ( column > gridColumns ) {
        column = gridColumns;
        xCenterOffset = 1.0;
      }

      if ( row > gridRows ) {
        row = gridRows;
        yCenterOffset = 1.0;
      }

      CHECK4( IN_RANGE( column, 1, data->columns ),
              IN_RANGE( row,    1, data->rows ),
              IN_RANGE( xCenterOffset, -1.0, 1.0 ),
              IN_RANGE( yCenterOffset, -1.0, 1.0 ) );
      DEBUG( fprintf(stderr, "  %lf %lf\n", fractionalColumn, fractionalRow);)
      DEBUG( fprintf( stderr, "  %"INTEGER_FORMAT" %"INTEGER_FORMAT
                      " %25.16lf %25.16lf\n",
              column, row, xCenterOffset, yCenterOffset ); )

      columns[ index ]        = column;
      rows[ index ]           = row;
      xCenterOffsets[ index ] = xCenterOffset;
      yCenterOffsets[ index ] = yCenterOffset;

      if ( gridLongitudes ) { /* Store unprojected grid cell center: */
        const Integer offset = ( row - 1 ) * gridColumns + ( column - 1 );
        CHECK( IN_RANGE( offset, 0, data->rows * data->columns - 1 ) );
        gridLongitudes[ index ] = gridCellCenterLongitudes[ offset ];
        gridLatitudes[  index ] = gridCellCenterLatitudes[  offset ];

        DEBUG( fprintf( stderr, " @ (%lf %lf)\n",
                        gridLongitudes[ index ], gridLatitudes[ index ] ); )
      }

      ++griddedPointCount;
    }
  }

  *griddedPoints = griddedPointCount;

  POST12( IN_RANGE( *griddedPoints, 0, count ),
          minimumItemI( columns, count ) >= 0,
          maximumItemI( columns, count ) <= gridColumns,
          minimumItemI( rows,    count ) >= 0,
          maximumItemI( rows,    count ) <= gridRows,
          IMPLIES_ELSE( columns[ 0 ] == 0, rows[ 0 ] == 0, rows[ 0 ] != 0 ),
          IMPLIES_ELSE( columns[ count - 1 ] == 0,
                        rows[ count - 1 ] == 0,
                        rows[ count - 1 ] != 0 ),
          minimumItem( xCenterOffsets, count ) >= -1.0,
          maximumItem( xCenterOffsets, count ) <=  1.0,
          minimumItem( yCenterOffsets, count ) >= -1.0,
          maximumItem( yCenterOffsets, count ) <=  1.0,
          IMPLIES( gridLongitudes,
                   AND4( isValidLongitude( minimumItem(gridLongitudes, count)),
                         isValidLongitude( maximumItem(gridLongitudes, count)),
                         isValidLatitude( minimumItem(gridLatitudes, count)),
                         isValidLatitude( maximumItem(gridLatitudes,count)))));


}



/******************************************************************************
PURPOSE: projectZ - Project elevation points (in meters above mean sea level)
         onto a vertical grid and obtain the 1-based grid cell layer numbers
         where the points lie and also the normalized [-1, 1] offsets from the
         cell centers.
INPUTS:  Grid* self              The Grid.
         Integer count           Number of elevation points.
         const Real elevations[] Elevations to project.
OUTPUTS: Integer* griddedPoints  Number of points projected onto grid.
         Integer layers[]        1-based grid layer number or 0 if outside.
         Real centerOffsets[]    [-1, 1] offset from center of grid cell.
         Real gridElevations[]   Optional elevations of regridded points or 0.
******************************************************************************/

static void projectZ( Grid* self, Integer count,
                      const Real elevations[], Integer* griddedPoints,
                      Integer layers[], Real centerOffsets[],
                      Real gridElevations[] ) {

  PRE6( self->layers( self ) > 0,
        count > 0, elevations, griddedPoints, layers, centerOffsets );

  const GridPrivate* const data = self->data;
  const Integer gridLayers = data->layers;
  const Real minimum       = data->z[ 0 ];
  const Real maximum       = data->z[ gridLayers ];
  Integer index = 0;
  Integer griddedPointCount = 0;

  *griddedPoints = 0;
  memset( layers,        0, count * sizeof *layers );
  memset( centerOffsets, 0, count * sizeof *centerOffsets );

  if ( gridElevations ) {
    memset( gridElevations, 0, count * sizeof *gridElevations );
  }

  DEBUG( fprintf( stderr, "[%lf %lf]\n", minimum, maximum ); )

#pragma omp parallel for reduction( + : griddedPointCount )

  for ( index = 0; index < count; ++index ) {
    const Real z = elevations[ index ];

    DEBUG( fprintf( stderr, "z = %lf\n", z); )

    if ( IN_RANGE( z, minimum, maximum ) ) {
      Integer layer = 0;

      for ( layer = 0; layer < gridLayers; ++layer ) {
        const Integer layerPlus1 = layer + 1;
        const Real layerMinimum = data->z[ layer ];
        const Real layerMaximum = data->z[ layer + 1 ];

        if ( IN_RANGE( z, layerMinimum, layerMaximum ) ) {
          const Real layerHeight = layerMaximum - layerMinimum;
          const Real fractionalLayer =
            (z - layerMinimum) / layerHeight + layerPlus1;
          Integer integerLayer = fractionalLayer; /* Truncate fraction. */
          Real centerOffset    = fractionalLayer - integerLayer - 0.5;
          centerOffset += centerOffset;

          if ( integerLayer > gridLayers ) {
            integerLayer = gridLayers;
            centerOffset = 1.0;
          }

          CHECK2( IN_RANGE( integerLayer, 1, data->layers ),
                  IN_RANGE( centerOffset, -1.0, 1.0 ) );

          DEBUG( fprintf(stderr, "  %lf %"INTEGER_FORMAT" %lf\n",
                         fractionalLayer, integerLayer, centerOffset ); )

          layers[ index ]        = integerLayer;
          centerOffsets[ index ] = centerOffset;

          if ( gridElevations ) {
            const Real gridZ = layerMinimum + 0.5 * layerHeight;
            gridElevations[ index ] = gridZ;
          }

          ++griddedPointCount;
        }
      }
    }
  }

  *griddedPoints = griddedPointCount;

  POST6( IN_RANGE( *griddedPoints, 0, count ),
         minimumItemI( layers, count ) >= 0,
         maximumItemI( layers, count ) <= gridLayers,
         minimumItem( centerOffsets, count ) >= -1.0,
         maximumItem( centerOffsets, count ) <=  1.0,
         IMPLIES( gridElevations,
                  AND2( minimumItem( gridElevations, count ) >= 0.0,
                        maximumItem( gridElevations, count ) < maximum ) ) );
}



/******************************************************************************
PURPOSE: aggregate - Aggregate grid-projected surface points in compact arrays.
INPUTS:  Grid* self                         The Grid.
         Integer method                     E.g., AGGREGATE_MEAN.
         Real    minimumValidValue          Minimum value of valid data.
         Integer inputPoints                Number of points to regrid.
         Integer columns[ inputPoints ]     1-based point grid columns.
         Integer rows[    inputPoints ]     1-based point grid rows.
         const Real xOffsets[ inputPoints ] Normalized cell center offsets.
         const Real yOffsets[ inputPoints ] Normalized cell center offsets.
         Real gridLongitudes[ inputPoints ] Longitudes of point cells.
         Real gridLatitudes[ inputPoints ]  Latitudes  of point cells.
         Integer layers                     Number of vertical data points.
         const Real elevations[ inputPoints * layers ]  Optional: m above MSL.
         const Real inputData[ inputPoints * layers ] Data to regrid.
         const Real inputData2[ inputPoints * layers ] 0 or vector 2nd component
                                                       to regrid.
         Integer*  outputPoints             Number of points that map onto grid
OUTPUTS: Integer*  outputPoints             Reduced (aggregated) point count.
         Integer columns[ outputPoints ]    Compact 1-based point grid columns.
         Integer rows[    outputPoints ]    Compact 1-based point grid rows.
         Real gridLongitudes[ outputPoints] Compact longitudes of point cells.
         Real gridLatitudes[ outputPoints ] Compact latitudes  of point cells.
         Real outputData[ outputPoints *
                          MAX( 1, ( elevations != 0 ) * self->layers( self ))]
                                            Regridded data.
         Real outputData2[ outputPoints *
                           MAX( 1, ( elevations != 0 ) * self->layers( self ))]
                                            0 or regridded vector 2nd component
         Real gridElevations[ outputPoints *
                              ( elevations != 0 ) * self->layers( self ) ]
                                            Optional regridded z.
******************************************************************************/

static void aggregate( Grid* self,
                       Integer method, Real minimumValidValue,
                       Integer inputPoints,
                       Integer columns[], Integer rows[],
                       const Real xOffsets[], const Real yOffsets[],
                       Real gridLongitudes[], Real gridLatitudes[],
                       const Integer layers, const Real elevations[],
                       const Real inputData[], const Real inputData2[],
                       Integer* outputPoints,
                       Real outputData[], Real outputData2[],
                       Real gridElevations[] ) {

  PRE30( IS_VALID_AGGREGATE_METHOD( method ),
         ! isNan( minimumValidValue ),
         inputPoints > 0, columns, rows,
         xOffsets, yOffsets, gridLongitudes, gridLatitudes,
         layers > 0,
         IMPLIES( layers > 1, elevations ),
         IMPLIES( gridElevations, elevations ),
         inputData,
         outputPoints, IN_RANGE( *outputPoints, 1, inputPoints ),
         outputData,
         IN_RANGE( minimumItemI( columns, inputPoints ),
                   0, self->columns( self ) ),
         IN_RANGE( maximumItemI( columns, inputPoints ),
                   1, self->columns( self ) ),
         IN_RANGE( minimumItemI( rows, inputPoints ),
                   0, self->rows( self ) ),
         IN_RANGE( maximumItemI( rows, inputPoints ),
                   1, self->rows( self ) ),
         IN_RANGE( minimumItem( xOffsets, inputPoints ), -1.0, 1.0 ),
         IN_RANGE( maximumItem( xOffsets, inputPoints ), -1.0, 1.0 ),
         IN_RANGE( minimumItem( yOffsets, inputPoints ), -1.0, 1.0 ),
         IN_RANGE( maximumItem( yOffsets, inputPoints ), -1.0, 1.0 ),
         isValidLongitude( minimumItem( gridLongitudes, inputPoints ) ),
         isValidLongitude( maximumItem( gridLongitudes, inputPoints ) ),
         isValidLatitude( minimumItem( gridLatitudes, inputPoints ) ),
         isValidLatitude( maximumItem( gridLatitudes, inputPoints ) ),
         isNanFree( inputData, inputPoints * layers ),
         IMPLIES( elevations, isNanFree( elevations, inputPoints * layers )));

  CHECKING( const Integer OLD( outputPoints ) = *outputPoints; )

  AggregatorEntry aggregatorEntry = aggregators[ method - 1 ];
  Aggregator      aggregator      = aggregatorEntry.aggregator;
  PreAggregator   preAggregator   = aggregatorEntry.preAggregator;
  PostAggregator  postAggregator  = aggregatorEntry.postAggregator;
  Integer aggregatedOutputPoints = 0;
  GridPrivate* const data = self->data;
  const Integer gridLayers  = elevations ? data->layers : 1;
  const Integer gridRows    = data->rows;
  const Integer gridColumns = data->columns;
  const Integer gridRowsTimesGridColumns = gridRows * gridColumns;
  const Integer gridColumnsTimesGridLayers = gridColumns * gridLayers;
  Cell* const gridCells = data->cells;
  Integer input = 0;
  Integer gridCell = 0;

  CHECK4( gridLayers <= MXLAYS3, aggregator, preAggregator, postAggregator );

  DEBUG( fprintf( stderr, "aggregate(): \n" ); )

  memset( outputData, 0, *outputPoints * gridLayers * sizeof *outputData );

  if ( gridElevations ) {
    memset( gridElevations, 0,
            *outputPoints * gridLayers * sizeof *gridElevations );
  }

  /* Aggregate each projected input into its grid cell: */

  initializeCells( gridRows * gridColumns * gridLayers,
                   minimumValidValue, gridCells );

  DEBUG( fprintf( stderr, "  loop on %lld inputPoints\n", inputPoints ); )

  for ( input = 0; input < inputPoints; ++input ) {
    const Integer row = rows[ input ];

    if ( row ) {
      const Integer column     = columns[ input ];
      const Real gridLongitude = gridLongitudes[ input ];
      const Real gridLatitude  = gridLatitudes[ input ];
      const Real xOffset       = xOffsets[ input ];
      const Real yOffset       = yOffsets[ input ];
      Cell* const cells =
        gridCells +
        ( row - 1 ) * gridColumnsTimesGridLayers + ( column - 1 ) * gridLayers;

      CHECK( IMPLIES_ELSE( cells[ 0 ].count > 0,
                           GT_ZERO2( cells[ 0 ].column, cells[ 0 ].row ),
                           IS_ZERO2( cells[ 0 ].column, cells[ 0 ].row ) ) );

      CHECK2( IN_RANGE( xOffset, -1.0, 1.0 ), IN_RANGE( yOffset, -1.0, 1.0 ) );

      aggregateCellData( self, preAggregator, aggregator,
                         column, row, gridLongitude, gridLatitude,
                         xOffset, yOffset,
                         input, inputData, inputData2, layers, elevations,
                         0, cells );
    }
  } /* End loop on input. */

  DEBUG( fprintf( stderr, "  loop on %lld 2D grid cells\n",
                  gridRowsTimesGridColumns ); )

  /* For each grid cell, post-aggregate and copy result to outputs: */

  for ( gridCell = 0; gridCell < gridRowsTimesGridColumns; ++gridCell ) {
    Cell* const cells = gridCells + gridCell * gridLayers;

    if ( AND2( cells[ 0 ].count, cells[ 0 ].data >= minimumValidValue ) ) {

      /* Store cell-aggregated output point: */

      const Integer offset = aggregatedOutputPoints * gridLayers;
      Integer layer = 0;

      CHECK6( isValidLongitude( cells[ 0 ].longitude ),
              isValidLatitude(  cells[ 0 ].latitude  ),
              IN_RANGE( cells[ 0 ].column, 1, self->columns( self ) ),
              IN_RANGE( cells[ 0 ].row,    1, self->rows( self ) ),
              ! isNan( cells[ 0 ].data ),
              cells[ 0 ].data >= minimumValidValue );

      gridLongitudes[ aggregatedOutputPoints ] = cells[ 0 ].longitude;
      gridLatitudes[  aggregatedOutputPoints ] = cells[ 0 ].latitude;
      columns[        aggregatedOutputPoints ] = cells[ 0 ].column;
      rows[           aggregatedOutputPoints ] = cells[ 0 ].row;

      for ( layer = 0; layer < gridLayers; ++layer ) {
        const Integer index = offset + layer;
        Cell* const cell = cells + layer;
        postAggregator( cell );
        CHECK( IN_RANGE( index, 0, *outputPoints * gridLayers - 1 ) );

        if ( cell->count ) {
          outputData[ index ] = cell->data;

          if ( outputData2 ) {
            outputData2[ index ] = cell->data2;
          }
        } else {
          outputData[ index ] = BADVAL3;

          if ( outputData2 ) {
            outputData2[ index ] = BADVAL3;
          }
        }

        if ( gridElevations ) {
          gridElevations[ index ] = cell->elevation;
        }
      }

      ++aggregatedOutputPoints;
    } /* End if aggregated cell. */
  } /* End loop on grid cells. */

  /* Clear any unused portion of output arrays: */

  zeroUnused( aggregatedOutputPoints, inputPoints,
              columns, rows, gridLongitudes, gridLatitudes );

  *outputPoints = aggregatedOutputPoints;

  DEBUG( fprintf( stderr, "Aggregated %lld points.\n", *outputPoints ); )

  POST2( IN_RANGE( *outputPoints, 0, OLD( outputPoints ) ),
         IMPLIES_ELSE( *outputPoints > 0,
                       AND10( IN_RANGE( minimumItemI( columns, *outputPoints ),
                                        1, self->columns( self ) ),
                              IN_RANGE( maximumItemI( columns, *outputPoints ),
                                        1, self->columns( self ) ),
                              IN_RANGE( minimumItemI( rows, *outputPoints ),
                                        1, self->rows( self ) ),
                              IN_RANGE( maximumItemI( rows, *outputPoints ),
                                        1, self->rows( self ) ),
                              isValidLongitude( minimumItem( gridLongitudes,
                                                *outputPoints ) ),
                              isValidLongitude( maximumItem( gridLongitudes,
                                                *outputPoints ) ),
                              isValidLatitude(  minimumItem( gridLatitudes,
                                                *outputPoints ) ),
                              isValidLatitude(  maximumItem( gridLatitudes,
                                                *outputPoints ) ),
                              isNanFree(outputData, *outputPoints *gridLayers),
                              IMPLIES( gridElevations,
                                       isNanFree( gridElevations,
                                                 *outputPoints * gridLayers))),
                       AND5( allZeroI( columns, inputPoints ),
                             allZeroI( rows, inputPoints ),
                             allZero( gridLongitudes, inputPoints ),
                             allZero( gridLatitudes, inputPoints ),
                             allZero( outputData, inputPoints ) ) ) );
}



/******************************************************************************
PURPOSE: regrid - Project and aggregate 2D/3D data points onto a 2D/3D grid.
INPUTS:  Grid* self                         Grid to project/aggregate onto.
         Integer method                     E.g., AGGREGATE_MEAN.
         Real minimumValidValue             Minimum value of valid data.
         Integer points                     Number of surface points to regrid.
         Integer levels                     Number of vertical data points.
         const Real longitudes[ points ]    Longitudes of surface point cells.
         const Real latitudes[ points ]     Latitudes  of surface point cells.
         const Real elevations[ points * levels ]  Optional: m above MSL.
         const Real inputData[ points * levels ]   Input data to regrid.
         const Real inputData2[ points * levels ]  0 or vector 2nd component
                                                   data to regrid.
         const Note notes[ points * levels ]       Optional: notes.
OUTPUTS: Integer*  regriddedPoints          Reduced (aggregated) point count.
         Integer columns[ regriddedPoints ]     1-based point grid columns.
         Integer rows[    regriddedPoints ]     1-based point grid rows.
         Integer layers[  regriddedPoints ] Optional 1-based point grid layers.
         Real gridLongitudes[ regriddedPoints]  Longitudes of point cells.
         Real gridLatitudes[ regriddedPoints ]  Latitudes  of point cells.
         Real gridElevations[ regriddedPoints *
                              ( MAX( 1, ( levels > 1 ) * self->layers( self ) ]
                                                Optional regridded z.
         Real outputData[ regriddedPoints * ... ] Regridded data.
         Real outputData2[ regriddedPoints * ... ] 0 or regridded vector 2nd
                                                   component data.
         RegriddedNote regriddedNotes[ regriddedPoints  * ... ] Optional.
******************************************************************************/

static void regrid( Grid* self,
                    Integer method, Real minimumValidValue,
                    Integer points, Integer levels,
                    const Real longitudes[], const Real latitudes[],
                    const Real elevations[],
                    const Real inputData[], const Real inputData2[],
                    const Note notes[],
                    Integer* regriddedPoints,
                    Integer columns[], Integer rows[], Integer layers[],
                    Real gridLongitudes[], Real gridLatitudes[],
                    Real gridElevations[],
                    Real outputData[], Real outputData2[],
                    RegriddedNote regriddedNotes[] ) {

  PRE12( IS_VALID_AGGREGATE_METHOD( method ),
         ! isNan( minimumValidValue ),
         points > 0,
         levels > 0,
         NON_ZERO9( longitudes, latitudes, inputData, regriddedPoints,
                    columns, rows, gridLongitudes, gridLatitudes, outputData),
         isValidLongitude( minimumItem( longitudes, points ) ),
         isValidLongitude( maximumItem( longitudes, points ) ),
         isValidLatitude( minimumItem( latitudes, points ) ),
         isValidLatitude( maximumItem( latitudes, points ) ),
         IMPLIES( elevations, isNanFree( elevations, points * levels ) ),
         isNanFree( inputData, points * levels ),
         ( elevations != 0 ) == ( gridElevations != 0 ) );

  AggregatorEntry aggregatorEntry = aggregators[ method - 1 ];
  Aggregator      aggregator      = aggregatorEntry.aggregator;
  PreAggregator   preAggregator   = aggregatorEntry.preAggregator;
  PostAggregator  postAggregator  = aggregatorEntry.postAggregator;
  const GridPrivate* const data = self->data;
  Projector* const projector = data->projector;
  const Real xMinimum        = data->xMinimum;
  const Real xMaximum        = data->xMaximum;
  const Real yMinimum        = data->yMinimum;
  const Real yMaximum        = data->yMaximum;
  const Real oneOverWidth    = data->oneOverWidth;
  const Real oneOverHeight   = data->oneOverHeight;
  const Integer gridColumns  = data->columns;
  const Integer gridRows     = data->rows;
  const Integer gridLayers   = elevations ? data->layers : 1;
  const Integer gridRowsTimesGridColumns = gridRows * gridColumns;
  const Integer gridColumnsTimesGridLayers = gridColumns * gridLayers;
  Cell* const gridCells = data->cells;
  const Real* const gridCellCenterLongitudes = data->longitudes;
  const Real* const gridCellCenterLatitudes  = data->latitudes;
  Integer index = 0;
  Integer gridCell = 0;
  Integer outputPoints = 0;

  DEBUG( fprintf( stderr, "regrid(): levels = %lld elevation = %lf, "
                  "gridLayers = %lld, [%f %f][%f %f]\n",
                  levels, elevations ? elevations[ 0 ] : 0.0, gridLayers,
                 xMinimum, xMaximum, yMinimum, yMaximum ); )

  *regriddedPoints = 0;

  initializeCells( gridRows * gridColumns * gridLayers,
                   minimumValidValue, gridCells );

  DEBUG( fprintf( stderr, "  loop on %lld points\n", points ); )

/* #pragma omp parallel for */

  for ( index = 0; index < points; ++index ) {
    const Real longitude = longitudes[ index ];
    const Real latitude  = latitudes[ index ];
    Real x = longitude;
    Real y = latitude;

    if ( projector ) {
      projector->project( projector, longitude, latitude, &x, &y );
    }

    DEBUG( fprintf( stderr, "(%lf %lf)->(%lf, %lf)\n",
                    longitude, latitude, x, y ); )

    if ( AND2( IN_RANGE( x, xMinimum, xMaximum ),
               IN_RANGE( y, yMinimum, yMaximum ) ) ) {
      const Real fractionalColumn = ( x - xMinimum ) * oneOverWidth  + 1.0;
      const Real fractionalRow    = ( y - yMinimum ) * oneOverHeight + 1.0;
      Integer column     = fractionalColumn; /* Truncate fraction. */
      Integer row        = fractionalRow;
      Real xCenterOffset = fractionalColumn - column - 0.5;
      Real yCenterOffset = fractionalRow    - row    - 0.5;
      xCenterOffset += xCenterOffset;
      yCenterOffset += yCenterOffset;

      if ( column > gridColumns ) {
        column = gridColumns;
        xCenterOffset = 1.0;
      }

      if ( row > gridRows ) {
        row = gridRows;
        yCenterOffset = 1.0;
      }

      CHECK4( IN_RANGE( column, 1, data->columns ),
              IN_RANGE( row,    1, data->rows ),
              IN_RANGE( xCenterOffset, -1.0, 1.0 ),
              IN_RANGE( yCenterOffset, -1.0, 1.0 ) );
      DEBUG( fprintf(stderr, "  %lf %lf\n", fractionalColumn, fractionalRow);)
      DEBUG( fprintf( stderr, "  %"INTEGER_FORMAT" %"INTEGER_FORMAT
                      " %25.16lf %25.16lf\n",
                      column, row, xCenterOffset, yCenterOffset ); )

      {
        const Integer row1    = row    - 1;
        const Integer column1 = column - 1;
        const Integer offset2 = row1 * gridColumns + column1;
        const Integer offset3 =
          row1 * gridColumnsTimesGridLayers + column1 * gridLayers;
        const Real gridLongitude = gridCellCenterLongitudes[ offset2 ];
        const Real gridLatitude  = gridCellCenterLatitudes[  offset2 ];
        Cell* const cells = gridCells + offset3;
        CHECK( IN_RANGE( offset2, 0, data->rows * data->columns - 1 ) );
        CHECK( IN_RANGE( offset3, 0,
                         data->rows * data->columns * data->layers - 1 ) );
        lockCell( cells );
        CHECK( IMPLIES_ELSE( cells[ 0 ].count > 0,
                             GT_ZERO2( cells[ 0 ].column, cells[ 0 ].row ),
                             IS_ZERO2( cells[ 0 ].column, cells[ 0 ].row ) ) );

        aggregateCellData( self, preAggregator, aggregator,
                           column, row, gridLongitude, gridLatitude,
                           xCenterOffset, yCenterOffset,
                           index, inputData, inputData2, levels, elevations,
                           notes ? notes[ index ] : 0, cells );
        unlockCell( cells );
      }
    }
  } /* End loop on input points. */

  DEBUG( fprintf( stderr, "  loop on %lld 2D grid cells\n",
                  gridRowsTimesGridColumns ); )

  /* For each grid cell, post-aggregate and copy result to outputs: */

  DEBUG( fprintf( stderr, "levels = %lld, elevations = %p, layers = %p, "
                  "gridElevations = %p, gridLayers = %lld\n",
                  levels, elevations, layers, gridElevations, gridLayers ); )

  for ( gridCell = 0; gridCell < gridRowsTimesGridColumns; ++gridCell ) {
    Cell* const cells = gridCells + gridCell * gridLayers;

    if ( AND4( levels == 1, elevations, layers, gridElevations ) ) {

      /* Store a single layer cell-aggregated output point: */

      Integer layer = 0;

      for ( layer = 0; layer < gridLayers; ++layer ) {
        Cell* const cell = cells + layer;
        postAggregator( cell );

        if ( cell->count ) {
          CHECK7( isValidLongitude( cell->longitude ),
                  isValidLatitude(  cell->latitude  ),
                  IN_RANGE( cell->column, 1, self->columns( self ) ),
                  IN_RANGE( cell->row,    1, self->rows( self ) ),
                  IN_RANGE( layer + 1,    1, self->layers( self ) ),
                  ! isNan( cell->data ),
                  cell->data >= minimumValidValue );
          gridLongitudes[ outputPoints ] = cell->longitude;
          gridLatitudes[  outputPoints ] = cell->latitude;
          columns[        outputPoints ] = cell->column;
          rows[           outputPoints ] = cell->row;
          layers[         outputPoints ] = layer + 1;
          outputData[     outputPoints ] = cell->data;
          gridElevations[ outputPoints ] = cell->elevation;

          if ( regriddedNotes ) {
            CHECK( strlen( cell->regriddedNote ) < 256 );
            strcpy( regriddedNotes[ outputPoints ], cell->regriddedNote );
          }

          if ( outputData2 ) {
            outputData2[ outputPoints ] = cell->data2;
          }

          ++outputPoints;
        }
      }
    } else { /* Store each layer cell-aggregated output point: */
      Integer layer = 0;
      const Integer offset = outputPoints * gridLayers;

      for ( layer = 0; layer < gridLayers; ++layer ) {
        const Integer index = offset + layer;
        Cell* const cell = cells + layer;
        postAggregator( cell );


        if ( cell->count ) {
          CHECK( IMPLIES( gridLayers == 1, index < points * levels ) );
          outputData[ index ] = cell->data;
          gridLongitudes[ outputPoints ] = cell->longitude;
          gridLatitudes[  outputPoints ] = cell->latitude;
          columns[        outputPoints ] = cell->column;
          rows[           outputPoints ] = cell->row;

          if ( outputData2 ) {
            outputData2[ index ] = cell->data2;
          }

          ++outputPoints;
        } else if ( gridLayers > 1 ) {
          outputData[ index ] = BADVAL3;

          if ( outputData2 ) {
            outputData2[ index ] = BADVAL3;
          }
        }

        if ( gridElevations ) {
          gridElevations[ index ] = cell->elevation;
        }
      }
    } /* End if aggregated cell. */
  } /* End loop on grid cells. */

  finalizeCells( gridRows * gridColumns * gridLayers, gridCells );
  *regriddedPoints = outputPoints;

  DEBUG( fprintf( stderr, "Regridded to %lld points.\n", *regriddedPoints ); )

  POST2( IN_RANGE( *regriddedPoints, 0, points ),
         IMPLIES( *regriddedPoints > 0,
           AND11( IN_RANGE( minimumItemI( columns, *regriddedPoints ),
                            1, self->columns( self ) ),
                  IN_RANGE( maximumItemI( columns, *regriddedPoints ),
                            1, self->columns( self ) ),
                  IN_RANGE( minimumItemI( rows, *regriddedPoints ),
                            1, self->rows( self ) ),
                  IN_RANGE( maximumItemI( rows, *regriddedPoints ),
                            1, self->rows( self ) ),
                  IMPLIES( layers,
                           IN_RANGE( maximumItemI( layers, *regriddedPoints ),
                            1, self->layers( self ) ) ),
                  isValidLongitude( minimumItem( gridLongitudes,
                                    *regriddedPoints ) ),
                  isValidLongitude( maximumItem( gridLongitudes,
                                    *regriddedPoints ) ),
                  isValidLatitude(  minimumItem( gridLatitudes,
                                    *regriddedPoints ) ),
                  isValidLatitude(  maximumItem( gridLatitudes,
                                    *regriddedPoints ) ),
                  isNanFree( outputData,
                             *regriddedPoints *
                            ( levels > 1 ? self->layers(self) : 1 ) ),
                  IMPLIES( gridElevations,
                           isNanFree( gridElevations,
                                      *regriddedPoints *
                                     (levels > 1 ? self->layers(self) :1))))));
}




/******************************************************************************
PURPOSE: regridSwath - Project and aggregate 2D swath (quadrilateral) scalar
         data onto 2D grid cells covered by the swath quadrilaterals.
INPUTS:  Grid* self                        Grid to project/aggregate onto.
         const Integer method              AGGREGATE_MEAN, AGGREGATE_WEIGHTED.
         const const Real minimumValidValue  Minimum valid mean data value.
         Integer points                     Number of surface points to regrid.
         const Real longitudesSW[ points ]  Longitudes of SW corner of quads.
         const Real longitudesSE[ points ]  Longitudes of SE corner of quads.
         const Real longitudesNW[ points ]  Longitudes of NW corner of quads.
         const Real longitudesNE[ points ]  Longitudes of NE corner of quads.
         const Real latitudesSW[  points ]  Latitudes  of SW corner of quads.
         const Real latitudesSE[  points ]  Latitudes  of SE corner of quads.
         const Real latitudesNW[  points ]  Latitudes  of NW corner of quads.
         const Real latitudesNE[  points ]  Latitudes  of NE corner of quads.
         const Real inputData[ points * levels ]   Input data to regrid.
OUTPUTS: Integer*  regriddedPoints          Number of grid cells with data.
                                            This could be as high as the number
                                            of surface grid cells!
                                            grid->rows( grid ) *
                                            grid->columns( grid ).
         Integer gridColumns[ regriddedPoints ]  1-based cell grid columns.
         Integer gridRows[    regriddedPoints ]  1-based cell grid rows.
         Real gridLongitudes[ regriddedPoints]   Longitudes of cells.
         Real gridLatitudes[ regriddedPoints ]   Latitudes  of cells.
         Real outputData[ regriddedPoints ]      Regridded data.
******************************************************************************/

static void regridSwath( Grid* self,
                         const Integer method,
                         const Real minimumValidValue,
                         const Integer points,
                         const Real longitudesSW[],
                         const Real longitudesSE[],
                         const Real longitudesNW[],
                         const Real longitudesNE[],
                         const Real latitudesSW[],
                         const Real latitudesSE[],
                         const Real latitudesNW[],
                         const Real latitudesNE[],
                         const Real inputData[],
                         Integer* regriddedPoints,
                         Integer gridColumns[],
                         Integer gridRows[],
                         Real gridLongitudes[],
                         Real gridLatitudes[],
                         Real outputData[] ) {

  PRE21( IS_VALID_AGGREGATE_METHOD( method ),
         ! isNan( minimumValidValue ),
         points > 0,
         NON_ZERO15( longitudesSW, longitudesSE, longitudesNW, longitudesNE,
                     latitudesSW, latitudesSE, latitudesNW, latitudesNE,
                     inputData, regriddedPoints,
                     gridColumns, gridRows,
                     gridLongitudes, gridLatitudes, outputData),
         isValidLongitude( minimumItem( longitudesSW, points ) ),
         isValidLongitude( minimumItem( longitudesSE, points ) ),
         isValidLongitude( minimumItem( longitudesNW, points ) ),
         isValidLongitude( minimumItem( longitudesNE, points ) ),
         isValidLongitude( maximumItem( longitudesSW, points ) ),
         isValidLongitude( maximumItem( longitudesSE, points ) ),
         isValidLongitude( maximumItem( longitudesNW, points ) ),
         isValidLongitude( maximumItem( longitudesNE, points ) ),
         isValidLatitude( minimumItem( latitudesSW, points ) ),
         isValidLatitude( minimumItem( latitudesSE, points ) ),
         isValidLatitude( minimumItem( latitudesNW, points ) ),
         isValidLatitude( minimumItem( latitudesNE, points ) ),
         isValidLatitude( maximumItem( latitudesSW, points ) ),
         isValidLatitude( maximumItem( latitudesSE, points ) ),
         isValidLatitude( maximumItem( latitudesNW, points ) ),
         isValidLatitude( maximumItem( latitudesNE, points ) ),
         isNanFree( inputData, points ) );

  /* Project corners in counter-clockwise order to contiguous x, y arrays: */

  double* vx = NEW_ZERO( Real, points * 8 );

  *regriddedPoints = 0;

  if ( vx ) {
    double* vy = vx + 4 * points;
    const GridPrivate* const gridPrivate = self->data;
    Projector* const projector = gridPrivate->projector;
    projectAndOrReorderQuadrilateralVertices( (size_t) points,
                                              longitudesSW, longitudesSE,
                                              longitudesNW, longitudesNE,
                                              latitudesSW,  latitudesSE,
                                              latitudesNW,  latitudesNE,
                                              projector,
                                              projector ?
                                                (ProjectFunction)
                                                projector->project
                                                : 0,
                                              vx, vy );

    const size_t rows    = gridPrivate->rows;
    const size_t columns = gridPrivate->columns;
    const size_t gridCells = rows * columns;
    size_t* cellCounts = NEW_ZERO( size_t, gridCells );
    double* cellMeans  = cellCounts ? NEW_ZERO( double, gridCells ) : 0;
    double* cellWeights =
      AND2( method == AGGREGATE_WEIGHTED, cellMeans ) ?
        NEW_ZERO( double, gridCells )
      : 0;

    if ( AND2( cellMeans,
               IMPLIES( method == AGGREGATE_WEIGHTED, cellWeights ) ) ) {
      const double gridXMinimum   = gridPrivate->xMinimum;
      const double gridYMinimum   = gridPrivate->yMinimum;
      const double cellWidth  = gridPrivate->cellWidth;
      const double cellHeight = gridPrivate->cellHeight;

      /* Bin the swath data into grid cells: */

      const size_t binnedPoints =
        binQuadrilateralData( points, inputData, vx, vy,
                              rows, columns,
                              gridXMinimum, gridYMinimum, cellWidth, cellHeight,
                              cellCounts, cellWeights, cellMeans );

      if ( binnedPoints ) { /* At least 1 point was within the grid: */

        /*
         * Compute mean value in each non-empty cell
         * but filter-out any cells whose mean is less than minimumValidValue:
         */

        *regriddedPoints =
          computeCellMeans( minimumValidValue, points,
                            cellCounts, cellWeights, cellMeans );

        CHECK( *regriddedPoints <= gridCells ); // Count of cells with data.

        if ( *regriddedPoints > 0 ) { /* At least 1 cell mean is valid: */

          /*
           * Compute compact arrays of
           * data->gridLongitudes, data->gridLatitudes,
           * data->columns, data->rows and
           * data->counts, data->gridData
           * so afterwards, all these arrays have length outputCells - i.e.,
           * they contain no missing data from index [0, outputCells - 1]:
           */

          memcpy( outputData, cellMeans, gridCells * sizeof (double) );
          compactCells( projector,
                       (UnprojectFunction)
                         ( projector ? projector->unproject : 0 ),
                       columns, rows,
                       gridXMinimum, gridYMinimum, cellWidth, cellHeight,
                       *regriddedPoints,
                       (size_t*) cellCounts, outputData,
                       gridLongitudes, gridLatitudes,
                       (size_t*) gridColumns, (size_t*) gridRows );
        }
      }
    }

    FREE( cellCounts );
    FREE( cellWeights );
    FREE( cellMeans );
  }

  FREE( vx );

  DEBUG( fprintf( stderr, "Regridded to %lld points.\n", *regriddedPoints ); )

  POST2( IN_RANGE( *regriddedPoints, 0,
                   self->rows( self ) * self->columns( self ) ),
         IMPLIES( *regriddedPoints > 0,
           AND9( IN_RANGE( minimumItemI( gridColumns, *regriddedPoints ),
                           1, self->columns( self ) ),
                 IN_RANGE( maximumItemI( gridColumns, *regriddedPoints ),
                           1, self->columns( self ) ),
                 IN_RANGE( minimumItemI( gridRows, *regriddedPoints ),
                           1, self->rows( self ) ),
                 IN_RANGE( maximumItemI( gridRows, *regriddedPoints ),
                           1, self->rows( self ) ),
                 isValidLongitude( minimumItem( gridLongitudes,
                                                *regriddedPoints ) ),
                 isValidLongitude( maximumItem( gridLongitudes,
                                                *regriddedPoints ) ),
                 isValidLatitude(  minimumItem( gridLatitudes,
                                                *regriddedPoints ) ),
                 isValidLatitude(  maximumItem( gridLatitudes,
                                                *regriddedPoints ) ),
                 isNanFree( outputData, *regriddedPoints ) ) ) );
}



/******************************************************************************
PURPOSE: free - Destruct a Grid.
INPUTS:  Grid* self  Object to destruct.
NOTE:    Use FREE_OBJECT( object ) instead since it zeros argument.
******************************************************************************/

static void free__( Grid* self ) {
  PRE( self );

  if ( self->data ) {
    FREE_OBJECT( self->data->projector );
    FREE_ZERO( self->data->longitudes );
    FREE_ZERO( self->data->z );
    FREE_ZERO( self->data->cells );
    FREE_ZERO( self->data );
  }

  FREE_ZERO( self );
  POST0( self == 0 );
}



/******************************************************************************
PURPOSE: invariant - Class invariant.
INPUTS:  const Grid* self  Object to check.
RETURNS: Integer 1 if valid, else 0.
NOTE:    If this query ever returns 0 then there is a defect in the code.
******************************************************************************/

static Integer invariant( const Grid* self ) {

  const GridPrivate* const data = self ? self->data : 0;
  const Integer result =
    AND22( self,
           hasMembers( self ),
           self->data,
           data == self->data,
           data->columns > 0,
           data->rows > 0,
           data->rows < INTEGER_MAX / data->columns,
           ! isNan( data->xMinimum ),
           ! isNan( data->yMinimum ),
           ! isNan( data->cellWidth ),
           ! isNan( data->cellHeight ),
           data->cellWidth > 0.0,
           data->cellHeight > 0.0,
           IMPLIES_ELSE( data->projector,
                         data->projector->invariant( data->projector ),
                         AND2( isValidLongitudeLatitude( data->xMinimum,
                                                         data->yMinimum ),
                               isValidLongitudeLatitude(
                                 data->xMinimum +
                                   data->columns * data->cellWidth,
                                 data->yMinimum +
                                   data->rows * data->cellHeight ) ) ),
           ! isNan( data->g ),
           ! isNan( data->R ),
           ! isNan( data->A ),
           ! isNan( data->T0s ),
           ! isNan( data->P00 ),
           data->cells,
           data->layers >= 0,
           IMPLIES_ELSE( data->layers == 0,
             IS_ZERO9( data->type, data->topPressure, data->z, data->levels,
                       data->g, data->R, data->A, data->T0s, data->P00 ),
             AND11( IS_VALID_VERTICAL_GRID_TYPE( data->type ),
                    data->type != IMISS3,
                    ! isNan( data->topPressure ),
                    data->topPressure > 0.0,
                    data->levels,
                    data->z,
                    isNanFree( data->z, data->layers + 1 ),
                    increasing( data->z, data->layers + 1 ),
                    isNanFree( data->levels, data->layers + 1 ),
                    IMPLIES_ELSE( IN4( data->type, VGSIGZ3, VGZVAL3, VGHVAL3 ),
                                  AND3( increasing( data->levels,
                                                    data->layers + 1 ),
                                        minimumItem( data->levels,
                                                     data->layers + 1)
                                                    >= -1000.0,
                                        maximumItem( data->levels,
                                                     data->layers + 1)
                                                    <= 1000000.0 ),
                                  AND3( decreasing( data->levels,
                                                    data->layers + 1 ),
                                        minimumItem( data->levels,
                                                     data->layers + 1) >= 0.0,
                                        maximumItem( data->levels,
                                                     data->layers + 1) <=1.0)),
                    GT_ZERO5( data->g, data->R, data->A, data->T0s,
                              data->P00 ) ) ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: equal - Is self functionally equivalent to other?
INPUTS:  const Grid* self   Object to compare.
         const Grid* other  Object to compare.
RETURNS: Integer 1 if equal, else 0.
******************************************************************************/

static Integer equal( const Grid* self, const Grid* other ) {

  PRE2( other, other->invariant( other ) );

  const GridPrivate* const data = self->data;
  const GridPrivate* const otherData = other->data;
  const Integer result =
    AND18( data->columns == otherData->columns,
           data->rows    == otherData->rows,
           aboutEqual( data->xMinimum,   otherData->xMinimum ),
           aboutEqual( data->yMinimum,   otherData->yMinimum ),
           aboutEqual( data->xMaximum,   otherData->xMaximum ),
           aboutEqual( data->yMaximum,   otherData->yMaximum ),
           aboutEqual( data->cellWidth,  otherData->cellWidth ),
           aboutEqual( data->cellHeight, otherData->cellHeight ),
           OR2( IS_ZERO2( data->projector, otherData->projector ),
                data->projector->equal(data->projector, otherData->projector)),
           data->layers  == otherData->layers,
           data->type    == otherData->type,
           data->topPressure == otherData->topPressure,
           data->g == otherData->g,
           data->R == otherData->R,
           data->A == otherData->A,
           data->T0s == otherData->T0s,
           data->P00 == otherData->P00,
           IMPLIES_ELSE( data->layers == 0,
                         IS_ZERO4( data->z, otherData->z,
                                   data->levels, otherData->levels ),
                         AND2( aboutEqual(
                                 minimumItem( data->levels, data->layers + 1 ),
                                 minimumItem( otherData->levels,
                                              data->layers + 1 ) ),
                               aboutEqual(
                                 maximumItem( data->levels, data->layers + 1 ),
                                 maximumItem( otherData->levels,
                                              data->layers + 1 ) ) ) ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: clone - Yield a clone - a NEW functionally equivalent Grid.
INPUTS:  const Grid* self  Object to clone.
RETURNS: Grid* result clone of self.
******************************************************************************/

static Grid* clone( const Grid* self ) {

  PRE( self );

  const GridPrivate* const data = self->data;
  Grid* result =
  newGrid( data->projector ? data->projector->clone( data->projector ) : 0,
           data->columns, data->rows,
           data->xMinimum, data->yMinimum,
           data->cellWidth, data->cellHeight,
           0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0 );

  if ( AND2( result, data->layers ) ) {
    const Integer threads = omp_get_max_threads();
    const Integer layers = data->layers;
    const Integer zCount = ( threads + 1 ) * ( layers + 1 );
    result->data->z = NEW_ZERO( Real, zCount );
    result->data->cells =
      result->data->z ?
        NEW_ZERO( Cell, data->rows * data->columns * data->layers )
      : 0;

    if ( ! result->data->z ) {
      FREE_OBJECT( result );
    } else {
      result->data->layers = layers;
      memcpy( result->data->z, data->z, zCount * sizeof result->data->z[ 0 ] );
    }
  }

  POST( IMPLIES( result,
                 AND2( result->invariant( result ),
                       result->equal( result, self ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: projector - The cartographic projector.
INPUTS:  const Grid* self  Object to query.
RETURNS: Projector* The cartographic projector.
******************************************************************************/

static Projector* projector( const Grid* self ) {
  PRE( self );
  Projector* const result = self->data->projector;
  POST( IMPLIES( result, result->invariant( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: layers - Number of layers of grid cells or 0 if 2D.
INPUTS:  const Grid* self  Object to query.
RETURNS: Integer number of layers of grid cells.
******************************************************************************/

static Integer layers( const Grid* self ) {
  PRE( self );
  const Integer result = self->data->layers;
  POST( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: rows - Number of rows of grid cells.
INPUTS:  const Grid* self  Object to query.
RETURNS: Integer number of rows of grid cells.
******************************************************************************/

static Integer rows( const Grid* self ) {
  PRE( self );
  const Integer result = self->data->rows;
  POST( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: columns - Number of columns of grid cells.
INPUTS:  const Grid* self  Object to query.
RETURNS: Integer number of columns of grid cells.
******************************************************************************/

static Integer columns( const Grid* self ) {
  PRE( self );
  const Integer result = self->data->columns;
  POST( result > 0 );
  return result;
}



/******************************************************************************
PURPOSE: longitude - Longitude at the center of a grid cell.
INPUTS:  const Grid* self  Object to query.
         Integer row        0-Based row    number.
         Integer column     0-Based column number.
RETURNS: Real longitude at the center of the grid cell.
******************************************************************************/

static Real longitude( const Grid* self, Integer row, Integer column ) {

  PRE2( IN_RANGE( row,    0, self->rows( self ) - 1 ),
        IN_RANGE( column, 0, self->columns( self ) - 1 ) );

  const GridPrivate* const data = self->data;
  const Integer columns = data->columns;
  const Integer index = row * columns + column;
  const Real result = data->longitudes[ index ];

  POST( isValidLongitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: latitude - Latitude at the center of a grid cell.
INPUTS:  const Grid* self  Object to query.
         Integer row        0-Based row    number.
         Integer column     0-Based column number.
RETURNS: Real latitude at the center of the grid cell.
******************************************************************************/

static Real latitude( const Grid* self, Integer row, Integer column ) {

  PRE2( IN_RANGE( row,    0, self->rows( self ) - 1 ),
        IN_RANGE( column, 0, self->columns( self ) - 1 ) );

  const GridPrivate* const data = self->data;
  const Integer columns = data->columns;
  const Integer index = row * columns + column;
  const Real result = data->latitudes[ index ];

  POST( isValidLatitude( result ) );
  return result;
}



/******************************************************************************
PURPOSE: elevation - Elevation at the center of a grid cell.
INPUTS:  const Grid* self  Object to query.
         Integer layer  0-Based layer number.
RETURNS: Real elevation, in meters above mean sea level, at the center of the
         grid cell.
******************************************************************************/

static Real elevation( const Grid* self, Integer layer ) {

  PRE( IN_RANGE( layer, 0, self->layers( self ) - 1 ) );

  const GridPrivate* const data = self->data;
  const Integer layers = data->layers;
  const Real result =
    layers ? 0.5 * ( data->z[ layer ] + data->z[ layer + 1 ] ) : 0.0;

  POST( IMPLIES_ELSE( self->layers( self ), result >= -1000.0, result == 0.0));
  return result;
}



/******************************************************************************
PURPOSE: level - Level at the edges of the grid cell.
INPUTS:  const Grid* self  Object to query.
         Integer level      0-Based level number.
RETURNS: Real level, in originally-specified units (e.g., sigma-pressure).
******************************************************************************/

static Real level( const Grid* self, Integer level ) {

  PRE( IN_RANGE( level, 0, self->layers( self ) ) );

  const GridPrivate* const data = self->data;
  const Real result = data->levels[ level ];

  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: westEdge - West edge of the grid in meters from the projection center.
INPUTS:  const Grid* self  Object to query.
RETURNS: Real West edge of the grid in meters from the projection center.
******************************************************************************/

static Real westEdge( const Grid* self ) {
  PRE( self );
  const Real result = self->data->xMinimum;
  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: southEdge - South edge of the grid in meters from the projection
         center.
INPUTS:  const Grid* self  Object to query.
RETURNS: Real South edge of the grid in meters from the projection center.
******************************************************************************/

static Real southEdge( const Grid* self ) {
  PRE( self );
  const Real result = self->data->yMinimum;
  POST( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: cellWidth - Width of grid cells in meters.
INPUTS:  const Grid* self  Object to query.
RETURNS: Real Width of grid cells in meters.
******************************************************************************/

static Real cellWidth( const Grid* self ) {
  PRE( self );
  const Real result = self->data->cellWidth;
  POST2( ! isNan( result ), result > 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: cellHeight - Height of grid cells in meters.
INPUTS:  const Grid* self  Object to query.
RETURNS: Real Height of grid cells in meters.
******************************************************************************/

static Real cellHeight( const Grid* self ) {
  PRE( self );
  const Real result = self->data->cellHeight;
  POST2( ! isNan( result ), result > 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: parseGrid - Parse grid command-line arguments.
INPUTS:  int argc              Number of command-line arguments.
         char* argv[ argc ]    Command-line argument strings.
         Integer* argument     Index of argument to begin parsing.
         Projector* projector  Projector to use (take ownership of).
OUTPUTS: Integer* argument     Updated index of argument to parse next.
         Projector* projector  Zero.
RETURNS: Grid* grid descibed by arguments.
******************************************************************************/

Grid* parseGrid( int argc, char* argv[], Integer* argument,
                 Projector* projector ) {

  PRE07( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], argument,
         IN_RANGE( *argument, 1, argc ),
         IMPLIES( projector, projector->invariant( projector ) ) );

  Grid* result = 0;
  Integer columns     = 0;
  Integer rows        = 0;
  Integer layers      = 0;
  Real    westEdge    = 0.0;
  Real    southEdge   = 0.0;
  Real    cellWidth   = 0.0;
  Real    cellHeight  = 0.0;
  Integer type        = 0;
  Real    topPressure = 0.0;
  Real    g           = 0.0;
  Real    R           = 0.0;
  Real    A           = 0.0;
  Real    T0s         = 0.0;
  Real    P00         = 0.0;
  Real*   levels = 0;
  Integer ok = 0;

  if ( ! AND3( *argument + 6 < argc,
               argv[ *argument ], ! strcmp( argv[ *argument ], "-grid" ) ) ) {
    failureMessage( "Invalid -grid command-line argument '%s'.",
                    argv[ *argument ] );
  } else {
    const Integer failures = failureCount();
    ++*argument; /* Skip "-grid". */
    columns = toInteger( argv[ *argument ], 1, INTEGER_MAX, &ok );

    if ( ok ) {
      ++*argument;
      rows = toInteger( argv[ *argument ], 1, INTEGER_MAX / columns, &ok );

      if ( ok ) {
        ++*argument;
        westEdge = toReal( argv[ *argument ], -1e8, 1e8, &ok );
        ok = AND2( ok,
                   IMPLIES( projector == 0, isValidLongitude( westEdge ) ) );

        if ( ok ) {
          ++*argument;
          southEdge = toReal( argv[ *argument ], -1e8, 1e8, &ok );
          ok = AND2( ok,
                     IMPLIES( projector == 0, isValidLatitude( southEdge ) ) );

          if ( ok ) {
            ++*argument;
            cellWidth = toReal( argv[ *argument ], 0.001, 1000000.0, &ok );

            if ( ok ) {
              ++*argument;
              cellHeight = toReal( argv[ *argument ], 0.001, 1000000.0, &ok );

              if ( ok ) {
                ++*argument;
                ok = IMPLIES( projector == 0,
                  isValidLongitudeLatitude( westEdge + columns * cellWidth,
                                            southEdge + rows * cellHeight ) );

                if ( ! ok ) {
                  failureMessage( "Invalid -grid." );
                } else if ( AND2( *argument < argc,
                           ! strcmp( argv[ *argument ], "-layers" ) ) ) {
                  levels = parseLayers( argc, argv, argument,
                                        &g, &R, &A, &T0s, &P00,
                                        &layers, &type, &topPressure );
                  ok = levels != 0;
                }
              }
            }
          }
        }
      }
    }

    if ( AND2( ! ok, failures == failureCount() ) ) {
      failureMessage( "Invalid argument '%s'.", argv[ *argument ] );
    }
  }

  if ( ok ) {
    result = newGrid( projector,
                      columns, rows,
                      westEdge, southEdge,
                      cellWidth, cellHeight,
                      layers, type, topPressure, levels,
                      g, R, A, T0s, P00 );
  }

  FREE( levels );
  POST0( IMPLIES( result, result->invariant( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseProjection - Parse projection command-line arguments.
INPUTS:  int argc       Number of command-line argument strings.
         char* argv[]   Command-line argument strings.
         int* argument  Index of current argument to parse.
OUTPUTS: int* argument  Index of next argument to parse.
RETURNS: Projector* cartographic projector specified by arguments.
******************************************************************************/

Projector* parseProjection( int argc, char* argv[], Integer* argument ) {

  PRE06( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], argument,
         IN_RANGE( *argument, 1, argc ) );

  Projector* result = 0;
  Integer ok = 0;

  if ( AND3( *argument + 2 < argc, argv[ *argument ], argv[ *argument ][0])) {
    Real majorSemiaxis    = 0.0;
    Real minorSemiaxis    = 0.0;
    Real centralLongitude = 0.0;
    Real centralLatitude  = 0.0;
    Real lowerLatitude    = 0.0;
    Real upperLatitude    = 0.0;
    Real secantLatitude   = 0.0;
    Integer itemsToParse    = 2;
    Integer isLambert       = 0;
    Integer isMercator      = 0;
    Integer isStereographic = 0;

    do {
      const char* const option = argv[ *argument ];

      if ( ! strcmp( option, "-ellipsoid" ) ) {
        parseEllipsoid( argc, argv, argument,
                        &majorSemiaxis, &minorSemiaxis, &ok );
      } else if ( ! strcmp( option, "-lambert" ) ) {
        parseLambert( argc, argv, argument,
                      &lowerLatitude,    &upperLatitude,
                      &centralLongitude, &centralLatitude, &ok );
        isLambert = 1;
      } else if ( ! strcmp( option, "-mercator" ) ) {
        parseMercator( argc, argv, argument, &centralLongitude, &ok );
        isMercator = 1;
      } else if ( ! strcmp( option, "-stereographic" ) ) {
        parseStereographic( argc, argv, argument,
                            &centralLongitude, &centralLatitude,
                            &secantLatitude, &ok );
        isStereographic = 1;
      }

    } while ( AND2( ok, --itemsToParse ) );

    if ( AND2( ok, majorSemiaxis > 0.0 ) ) {

      if ( isLambert ) {
        result = (Projector*)
          newLambert( majorSemiaxis, minorSemiaxis,
                      lowerLatitude, upperLatitude,
                      centralLongitude, centralLatitude, 0.0, 0.0 );
      } else if ( isMercator ) {
        result = (Projector*)
          newMercator( majorSemiaxis, minorSemiaxis, centralLongitude,
                       0.0, 0.0 );
      } else if ( isStereographic ) {
        result = (Projector*)
          newStereographic( majorSemiaxis, minorSemiaxis,
                            centralLongitude, centralLatitude,
                            secantLatitude, 0.0, 0.0 );
      }
    }
  }

  if ( ! result ) {
    failureMessage("Invalid/insufficient projection command-line arguments.");
  }

  POST0( IMPLIES( result, result->invariant( result ) ) );
  return result;
}



/******************************************************************************
PURPOSE: parseEllipsoid - Parse ellipsoid command-line arguments.
INPUTS:  int argc             Number of command-line argument strings.
         char* argv[]         Command-line argument strings.
         int* argument        Index of current argument to parse.
OUTPUTS: int* argument        Index of next argument to parse.
         Real* majorSemiaxis  Planet mean equitorial radius in meters.
         Real* minorSemiaxis  Planet mean polar      radius in meters.
         Integer* ok          1 if successful, else 0.
******************************************************************************/

void parseEllipsoid( int argc, char* argv[],
                     Integer* argument,
                     Real* majorSemiaxis, Real* minorSemiaxis,
                     Integer* ok ) {

  PRE09( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], argument,
         IN_RANGE( *argument, 1, argc ), majorSemiaxis, minorSemiaxis, ok );

  const Integer failures = failureCount();
  *ok = 0;
  *majorSemiaxis = *minorSemiaxis = 0.0;

  if ( ! AND2( *argument + 2 < argc,
               ! strcmp( argv[ *argument ], "-ellipsoid" ) ) ) {
    failureMessage( "Invalid -ellipsoid command-line argument '%s'.",
                    argv[ *argument ] );
  } else {
    ++*argument;
    *majorSemiaxis = toReal( argv[ *argument ], 1.0, REAL_MAX, ok );

    if ( *ok ) {
      ++*argument;
      *minorSemiaxis = toReal( argv[ *argument ], 1.0, *majorSemiaxis, ok );
      *argument += *ok;
    }
  }

  if ( ! *ok ) {

    if ( failures == failureCount() ) {
      failureMessage( "Invalid argument '%s'.", argv[ *argument ] );
    }

    *majorSemiaxis = *minorSemiaxis = 0.0;
  }

  POST03( IS_BOOL( *ok ),
          IN_RANGE( *argument, 0, argc ),
          IMPLIES_ELSE( *ok,
                        AND3( *majorSemiaxis >= 1.0, *minorSemiaxis >= 1.0,
                              *minorSemiaxis <= *majorSemiaxis ),
                        IS_ZERO2( *majorSemiaxis, *minorSemiaxis ) ) );
}



/******************************************************************************
PURPOSE: writeProjectionAndGrid - Write projection and grid info to stdout.
INPUTS:  const Grid* grid  Grid to write.
OUTPUTS: Stream* stream    Updated stream.
******************************************************************************/

void writeProjectionAndGrid( const Grid* grid, Stream* stream ) {
  PRE03( grid, grid->invariant( grid ), stream );

  /* Write portion of XDR header that looks like this:

    # lcc projection: lat_1 lat_2 lat_0 lon_0 major_semiaxis minor_semiaxis
    33 45 40 -97 6.37e+06 6.37e+06

  */

  const Projector* const projector = grid->projector( grid );

  if ( ! projector ) {
    const double sphereRadius = 6370000.0;
    stream->writeString( stream,
                         "# lonlat projection: "
                         "major_semiaxis minor_semiaxis\n%lg %lg\n",
                         sphereRadius, sphereRadius );
  } else {
    const char* const name = projector->name( projector );
    const Real centralLongitude = projector->centralLongitude( projector );
    const Real centralLatitude  = projector->centralLatitude(  projector );
    Real majorSemiaxis = 0.0;
    Real minorSemiaxis = 0.0;
    projector->ellipsoid( projector, &majorSemiaxis, &minorSemiaxis );

    if ( ! strcmp( name, "Lambert" ) ) {
      const Lambert* const lambert = (const Lambert*) projector;
      const Real lowerLatitude = lambert->lowerLatitude( lambert );
      const Real upperLatitude = lambert->upperLatitude( lambert );
      stream->writeString( stream,
                           "# lcc projection: "
                           "lat_1 lat_2 lat_0 lon_0 "
                           "major_semiaxis minor_semiaxis\n"
                           "%lg %lg %lg %lg %lf %lf\n",
                           lowerLatitude, upperLatitude,
                           centralLatitude, centralLongitude,
                           majorSemiaxis, minorSemiaxis );
    } else if ( ! strcmp( name, "Albers" ) ) { 
#if 0
      const Albers* const albers = (const Albers*) projector;
      const Real lowerLatitude = albers->lowerLatitude( albers );
      const Real upperLatitude = albers->upperLatitude( albers );
      stream->writeString( stream,
                           "# albers projection: "
                           "lat_1 lat_2 lat_0 lon_0 "
                           "major_semiaxis minor_semiaxis\n"
                           "%lg %lg %lg %lg %lf %lf\n",
                           lowerLatitude, upperLatitude,
                           centralLatitude, centralLongitude,
                           majorSemiaxis, minorSemiaxis );
#endif
    } else if ( ! strcmp( name, "Mercator" ) ) { 
      stream->writeString( stream,
                           "# mercator projection: lon_0 "
                           "major_semiaxis minor_semiaxis\n"
                           "%lg %lf %lf\n",
                           centralLongitude, majorSemiaxis, minorSemiaxis );
    } else if ( ! strcmp( name, "Stereographic" ) ) { 
      const Stereographic* const stereographic =
        (const Stereographic*) projector;
      const Real secantLatitude = stereographic->secantLatitude(stereographic);
      stream->writeString( stream,
                           "# stereographic projection:"
                           " lat_0 lon_0 lat_sec "
                           "major_semiaxis minor_semiaxis\n"
                           "%lg %lg %lg %lf %lf\n",
                           centralLatitude, centralLongitude, secantLatitude,
                           majorSemiaxis, minorSemiaxis );
    }
  }

  /* Now write portion of XDR header that looks like this:

    # Grid: ncols nrows xorig yorig xcell ycell vgtyp vgtop vglvls[25]:
    279 240 -1.008e+06 -1.62e+06 12000 12000 2 10000 1 0.995 0.99 0.98 0.97 \
      0.96 0.94 0.92 0.9 0.88 0.86 0.84 0.82 0.8 0.77 0.74 0.7 0.65 0.6 0.5 \
      0.4 0.3 0.2 0.1 0

  */

  if ( stream->ok( stream ) ) {
    const int levels = grid->layers( grid ) + 1;
    const GridPrivate* const gridPrivate = grid->data;
    const int vgtyp = gridPrivate->type;
    const double vgtop = gridPrivate->topPressure;
    const Real* const vglvls = gridPrivate->levels;
    int level = 0;

    stream->writeString( stream,
                         "# Grid: ncols nrows xorig yorig xcell ycell vgtyp "
                         "vgtop vglvls[%d]:\n"
                         "%d %d %lf %lf %lf %lf %d %lf",
                         levels,
                         (int) grid->columns( grid ), (int) grid->rows( grid ),
                         grid->westEdge( grid ), grid->southEdge( grid ),
                         grid->cellWidth( grid ), grid->cellHeight( grid ),
                         vgtyp, vgtop );

    for ( level = 0; AND2( stream->ok( stream ), level < levels ); ++level ) {
      const Real value = vglvls[ level ];
      stream->writeString( stream, " %lg", value );
    }

    stream->writeString( stream, "\n" );
  }
}



/*============================= PRIVATE FUNCTIONS ===========================*/



/******************************************************************************
PURPOSE: assignMembers - Assign pointers to member functions.
OUTPUTS: Grid* self  Grid to initialize.
******************************************************************************/

static void assignMembers( Grid* self ) {

  PRE0( self );

  self->free       = free__;
  self->projectXY  = projectXY;
  self->projectZ   = projectZ;
  self->aggregate  = aggregate;
  self->regrid     = regrid;
  self->regridSwath = regridSwath;
  self->invariant  = invariant;
  self->equal      = equal;
  self->clone      = clone;
  self->projector  = projector;
  self->layers     = layers;
  self->rows       = rows;
  self->columns    = columns;
  self->longitude  = longitude;
  self->latitude   = latitude;
  self->elevation  = elevation;
  self->level      = level;
  self->westEdge   = westEdge;
  self->southEdge  = southEdge;
  self->cellWidth  = cellWidth;
  self->cellHeight = cellHeight;

  POST0( hasMembers( self ) );
}



/******************************************************************************
PURPOSE: hasMembers - Does self have all pointers to member functions?
INPUTS:  const Grid* self  Grid to check.
RETURNS: Integer 1 if valid, else 0.
******************************************************************************/

static Integer hasMembers( const Grid* self ) {

  PRE0( self );

  const Integer result =
    AND21( self->free       == free__,
           self->projectXY  == projectXY,
           self->projectZ   == projectZ,
           self->aggregate  == aggregate,
           self->regrid     == regrid,
           self->regridSwath == regridSwath,
           self->invariant  == invariant,
           self->equal      == equal,
           self->clone      == clone,
           self->projector  == projector,
           self->rows       == rows,
           self->layers     == layers,
           self->columns    == columns,
           self->longitude  == longitude,
           self->latitude   == latitude,
           self->elevation  == elevation,
           self->level      == level,
           self->westEdge   == westEdge,
           self->southEdge  == southEdge,
           self->cellWidth  == cellWidth,
           self->cellHeight == cellHeight );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: computeLongitudesAndLatitudes - Compute grid cell center lon-lats.
INPUTS:  GridPrivate* data  Relevant parameters.
OUTPUTS: GridPrivate* data->longitudes, data->latitudes Initialized arrays.
NOTES:   If the grid uses a sphere (like CMAQ) then the computed latitudes are
         adjusted to be on the WGS84 spheroid for compatibility with other
         data (e.g., satellites).
******************************************************************************/

static void computeLongitudesAndLatitudes( GridPrivate* data ) {

  PRE013( data,
          data->columns > 0,
          data->rows > 0,
          data->rows < INTEGER_MAX / data->columns,
          ! isNan( data->xMinimum ),
          ! isNan( data->yMinimum ),
          ! isNan( data->cellWidth ),
          ! isNan( data->cellHeight ),
          data->cellWidth > 0.0,
          data->cellHeight > 0.0,
          IMPLIES( data->projector,
                   data->projector->invariant( data->projector ) ),
          data->longitudes,
          data->latitudes );

  Projector* const projector = data->projector;
  Real* const longitudes     = data->longitudes;
  Real* const latitudes      = data->latitudes;
  const Integer columns      = data->columns;
  const Integer rows         = data->rows;
  const Real cellWidth       = data->cellWidth;
  const Real cellHeight      = data->cellHeight;
  const Real xStart          = data->xMinimum - 0.5 * cellWidth;
  const Real yStart          = data->yMinimum - 0.5 * cellHeight;
  const Integer count = columns * rows;
  Real (*latitudeAdjuster)( Real ) = 0;
  Integer index = 0;

  if ( projector ) {
    Real majorSemiaxis = 0.0;
    Real minorSemiaxis = 0.0;
    projector->ellipsoid( projector, &majorSemiaxis, &minorSemiaxis );
    latitudeAdjuster = majorSemiaxis == minorSemiaxis ? latitudeWGS84 : 0;    
  } else {
    latitudeAdjuster = latitudeWGS84;    
  }

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    const Integer row    = index / columns;
    const Integer column = index % columns;
    const Integer offset = row * columns + column;
    const Real x = xStart + ( column + 1 ) * cellWidth;
    const Real y = yStart + ( row    + 1 ) * cellHeight;

    if ( projector ) {
      projector->unproject( projector, x, y,
                            longitudes + offset, latitudes + offset );
    } else {
      longitudes[ offset ] = x;
      latitudes[  offset ] = y;
    }

    if ( latitudeAdjuster ) {
      latitudes[ offset ] = latitudeAdjuster( latitudes[ offset ] );
    }
  }

  POST04( isValidLongitude(
            minimumItem( data->longitudes, data->columns * data->rows ) ),
          isValidLongitude(
            maximumItem( data->longitudes, data->columns * data->rows ) ),
          isValidLatitude(
            minimumItem( data->latitudes,  data->columns * data->rows ) ),
          isValidLongitude(
            maximumItem( data->latitudes,  data->columns * data->rows ) ) );
}




/******************************************************************************
PURPOSE: computeZ - Compute elevation in meters above mean sea level from
         CMAQ/IOAPI vertical grid parameters.
INPUTS:  Real g               Gravitational force, e.g., 9.81 m/s^2.
         Real R               Gas constant e.g., 287.04 J/kg/K = m^3/s/K.
         Real A               Atmospheric lapse rate, e.g., 50.0 K/kg.
         Real T0s             Reference surface temperature, e.g., 290.0 K.
         Real P00             Reference surface pressure, e.g., 100000 P.
         Integer layers       Number of elevation layers.
         Integer type         Vertical grid type:
                                VGSGPH3  hydrostatic sigma-P
                                VGSGPN3  non-h sigma-P
                                VGSIGZ3  sigma-Z
                                VGPRES3  pressure (pascals)
                                VGZVAL3  Z (m) (above sea lvl)
                                VGHVAL3  H (m) (above ground)
                                VGWRFEM  WRF sigma-P
                                IMISS3   None (becomes a single layer at 0.0).
         Real topPressure     Pressure in pascals at the top of the model.
         const Real levels[]  Vertical grid levels. levels[ layers + 1 ].
OUTPUTS: Real z[]             Elevation at each level.
RETURNS: Integer number of points that projected onto the grid.
******************************************************************************/

static void computeZ( Real g, Real R, Real A, Real T0s, Real P00,
                      Integer layers, Integer type, Real topPressure,
                      const Real levels[], Real z[] ) {

  PRE012( ! isNan( g ),
          ! isNan( R ),
          ! isNan( A ),
          ! isNan( T0s ),
          ! isNan( P00 ),
          GT_ZERO5( g, R, A, T0s, P00 ),
          layers > 0,
          IS_VALID_VERTICAL_GRID_TYPE( type ),
          topPressure > 0.0,
          levels,
          z,
          isNanFree( levels, layers + 1 ) );

  const Real HEIGHT_OF_TERRAIN_IN_METERS = 0.0;
  const Integer numberOfLevels = layers + 1;

  if ( IN4( type, VGSGPN3, VGSGPN3, VGWRFEM ) ) {
    /* Compute z using MM5 formula: */
    elevationsAtSigmaPressures( g, R, A, T0s, P00, HEIGHT_OF_TERRAIN_IN_METERS,
                                numberOfLevels, topPressure, levels, z );
  } else { /* Compute z using other formulas: */
    Integer level;

    for ( level = 0; level < numberOfLevels; ++level ) {
      const Real valueAtLevel = levels[ level ];
      Real pressure = 0.0;

      switch ( type ) {
        case VGSGPH3: /* hydrostatic sigma-P */
        case VGSGPN3: /* non-h sigma-P */
        case VGWRFEM: /* WRF sigma-P */
          pressure = pressureAtSigmaLevel( valueAtLevel, topPressure / 100.0 );
          z[ level ] = heightAtPressure( pressure );
          break;
        case VGSIGZ3: /* sigma-Z */
          /* vgtop is in meters and valueAtLevel increases for each level. */
          z[ level ] = HEIGHT_OF_TERRAIN_IN_METERS +
                       valueAtLevel *
                       ( topPressure - HEIGHT_OF_TERRAIN_IN_METERS );
          break;
        case VGPRES3: /* pressure (Pascals) */
          z[ level ] = heightAtPressure( valueAtLevel / 100.0 );
          break;
        case VGZVAL3: /* Z (m) (above sea level) */
          z[ level ] = valueAtLevel;
          break;
        case VGHVAL3: /* H (m) (above ground)  */
          z[ level ] = valueAtLevel + HEIGHT_OF_TERRAIN_IN_METERS;
          break;
        default:
          z[ level ] = level;
          break;
      }

      DEBUG( fprintf( stderr, "grid z[ %"INTEGER_FORMAT" ] = %lf\n",
                      level, z[ level ] ); )
    }
  }

  POST02( isNanFree( z, layers + 1 ), z[ 0 ] <= z[ layers ] );
}



/******************************************************************************
PURPOSE: pressureAtSigmaLevel - Compute pressure (in millibars) at a given
         sigma level.
INPUTS:  Real sigmaLevel     Sigma level.
         Real pressureAtTop  Pressure (in millibars) at top of the model.
OUTPUTS: None
RETURNS: Real pressure in millibars at the given sigma level.
NOTES:   Based on formula in the documentation for Vis5d by Bill Hibbard.
******************************************************************************/

static Real pressureAtSigmaLevel( Real sigmaLevel, Real pressureAtTop ) {
#define SURFACE_PRESSURE_IN_MB 1012.5
  const Real result =
    pressureAtTop + sigmaLevel * ( SURFACE_PRESSURE_IN_MB - pressureAtTop );
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: heightAtPressure - Compute the height (in meters) at a given
         pressure (in millibars).
INPUTS:  Real pressure   Pressure value in millibars.
OUTPUTS: None
RETURNS: Real height in meters at the given pressure in millibars.
NOTES:   Based on formula in the documentation for Vis5d by Bill Hibbard.
******************************************************************************/

static Real heightAtPressure( Real pressure ) {
  const Real pressureToHeightScaleFactor = -7.2 * 1000.0;
  Real result = 0.0;

  if ( pressure == 0.0 ) {
    pressure = 1e-10; /* HACK: prevent core on non-IEEE.*/
  }

  result =
    pressureToHeightScaleFactor * log( pressure / SURFACE_PRESSURE_IN_MB );

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: parseLambert - Parse Lambert projection command-line arguments.
INPUTS:  int argc       Number of command-line argument strings.
         char* argv[]   Command-line argument strings.
         int* argument  Index of current argument to parse.
OUTPUTS: int* argument  Index of next argument to parse.
         Real* lowerLatitude    Lower tangent in degrees, e.g., 30.0.
         Real* upperLatitude    Upper tangent in degrees, e.g., 60.0.
         Real* centralLongitude Projects to zero, e.g., -100.0 degrees.
         Real* centralLatitude  Projects to zero, e.g., 40.0 degrees.
         Integer* ok             1 if successful, else 0 and failureMessage().
******************************************************************************/

static void parseLambert( int argc, char* argv[],
                          Integer* argument,
                          Real* lowerLatitude,    Real* upperLatitude,
                          Real* centralLongitude, Real* centralLatitude,
                          Integer* ok ) {

  PRE011( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], argument,
          IN_RANGE( *argument, 1, argc ),
          lowerLatitude, upperLatitude, centralLongitude, centralLatitude,
          ok );

  *ok = 0;
  *lowerLatitude    = *upperLatitude   = 0.0;
  *centralLongitude = *centralLatitude = 0.0;

  if ( ! AND3( *argument + 4 < argc, argv[ *argument ],
               ! strcmp( argv[ *argument ], "-lambert" ) ) ) {
    failureMessage( "Invalid -lambert command-line argument '%s'.",
                    argv[ *argument ] );
  } else {
    ++*argument; /* Skip "-lambert". */
    *lowerLatitude = toReal( argv[ *argument ], -89.0, 89.0, ok );

    if ( ! *ok ) {
      failureMessage( "Invalid lowerLatitude '%s'.", argv[ *argument ] );
    } else {
      *ok = fabs( *lowerLatitude ) >= 1.0;

      if ( ! *ok ) {
        failureMessage( "Invalid lowerLatitude '%s'.", argv[ *argument ] );
      } else {
        const Real maximum = *lowerLatitude < 0.0 ? -1.0 : 89.0;
        ++*argument;
        *upperLatitude =
          toReal( argv[ *argument ], *lowerLatitude, maximum, ok );

        if ( ! *ok ) {
          failureMessage( "Invalid upperLatitude '%s'.", argv[ *argument ] );
        } else {
          ++*argument;
          parseCentralLongitudeAndLatitude( argc, argv, argument,
                                            centralLongitude,
                                            centralLatitude, ok );

          if ( *ok ) {
            *ok = IN_RANGE( *centralLatitude, -89.0, 89.0 );

            if ( ! *ok ) {
              failureMessage( "Invalid central-latitude." );
            }
          }
        }
      }
    }
  }

  if ( ! *ok ) {
    *lowerLatitude    = *upperLatitude   = 0.0;
    *centralLongitude = *centralLatitude = 0.0;
  }

  POST03( IS_BOOL( *ok ),
          IN_RANGE( *argument, 0, argc ),
          IMPLIES_ELSE( *ok,
                        AND4( IMPLIES_ELSE( *lowerLatitude >= 0.0,
                                  IN_RANGE( *lowerLatitude, 1.0, 89.0 ),
                                  IN_RANGE( *lowerLatitude, -89.0, -1.0 )),
                              IMPLIES_ELSE( *upperLatitude >= 0.0,
                                  IN_RANGE( *upperLatitude, 1.0, 89.0 ),
                                  IN_RANGE( *upperLatitude, -89.0, -1.0 )),
                              *lowerLatitude <= *upperLatitude,
                              isValidLongitudeLatitude( *centralLongitude,
                                                        *centralLatitude ) ),
                        IS_ZERO4( *lowerLatitude,    *upperLatitude,
                                  *centralLongitude, *centralLatitude ) ) );
}



/******************************************************************************
PURPOSE: parseMercator - Parse Mercator projection command-line arguments.
INPUTS:  int argc       Number of command-line argument strings.
         char* argv[]   Command-line argument strings.
         int* argument  Index of current argument to parse.
OUTPUTS: int* argument  Index of next argument to parse.
         Real* centralLongitude Projects to zero, e.g., -100.0 degrees.
         Integer* ok             1 if successful, else 0 and failureMessage().
******************************************************************************/

static void parseMercator( int argc, char* argv[],
                           Integer* argument,
                           Real* centralLongitude,
                           Integer* ok ) {

  PRE08( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], argument,
         IN_RANGE( *argument, 1, argc ),
         centralLongitude, ok );

  *ok = 0;
  *centralLongitude = 0.0;

  if ( ! AND3( *argument + 1 < argc, argv[ *argument ],
               ! strcmp( argv[ *argument ], "-mercator" ) ) ) {
    failureMessage( "Invalid -mercator command-line argument '%s'.",
                    argv[ *argument ] );
  } else {
    ++*argument; /* Skip "-mercator". */
    *centralLongitude = toReal( argv[ *argument ], -180.0, 180.0, ok );

    if ( ! *ok ) {
      failureMessage( "Invalid central-longitude '%s'.", argv[ *argument ] );
    } else {
      ++*argument;
    }
  }

  if ( ! *ok ) {
    *centralLongitude = 0.0;
  }

  POST03( IS_BOOL( *ok ),
          IN_RANGE( *argument, 0, argc ),
          IMPLIES_ELSE( *ok,
                        IN_RANGE( *centralLongitude, -180.0, 180.0 ),
                        *centralLongitude == 0 ) );
}



/******************************************************************************
PURPOSE: parseStereographic - Parse Stereographic command-line arguments.
INPUTS:  int argc                Number of command-line argument strings.
         char* argv[]            Command-line argument strings.
         int* argument           Index of current argument to parse.
OUTPUTS: int* argument           Index of next argument to parse.
         Real* centralLongitude  Projects to zero, e.g., -98.0 degrees.
         Real* centralLatitude   Projects to zero, e.g., 90.0 degrees.
         Real* secantLatitude    Secant latitude in degrees, e.g., 45.0.
         Integer* ok             1 if successful, else 0 and failureMessage().
******************************************************************************/

static void parseStereographic( int argc, char* argv[],
                                Integer* argument,
                                Real* centralLongitude, Real* centralLatitude,
                                Real* secantLatitude,
                                Integer* ok ) {

  PRE010( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], argument,
          IN_RANGE( *argument, 1, argc ),
          centralLongitude, centralLatitude, secantLatitude, ok );

  *ok = 0;
  *centralLongitude = *centralLatitude = *secantLatitude = 0;

  if ( ! AND3( *argument + 3 < argc, argv[ *argument ],
               ! strcmp( argv[ *argument ], "-stereographic" ) ) ) {
    failureMessage( "Invalid -stereographic command-line argument '%s'.",
                    argv[ *argument ] );
  } else {
    ++*argument; /* Skip "-stereographic". */
    parseCentralLongitudeAndLatitude( argc, argv, argument,
                                      centralLongitude, centralLatitude, ok );

    if ( *ok ) {
      *secantLatitude = toReal( argv[ *argument ], -90.0, 90.0, ok );
      *argument += *ok;
    }
  }

  if ( ! *ok ) {
    *centralLongitude = *centralLatitude = *secantLatitude = 0;
  }

  POST03( IS_BOOL( *ok ),
          IN_RANGE( *argument, 0, argc ),
          IMPLIES_ELSE( *ok,
                        AND2( isValidLongitudeLatitude( *centralLongitude,
                                                        *centralLatitude ),
                              isValidLatitude( *secantLatitude ) ),
                        IS_ZERO3( *centralLongitude, *centralLatitude,
                                  *secantLatitude ) ) );
}



/******************************************************************************
PURPOSE: parseCentralLongitudeAndLatitude - Parse central longitude and
         latitude command-line arguments.
INPUTS:  int argc                Number of command-line argument strings.
         char* argv[]            Command-line argument strings.
         int* argument           Index of current argument to parse.
OUTPUTS: int* argument           Index of next argument to parse.
         Real* centralLongitude  Central longitude of projection.
         Real* centralLatitude   Central latitude  of projection.
         Integer* ok             1 if successful, else 0 and failureMessage().
******************************************************************************/

static void parseCentralLongitudeAndLatitude( int argc, char* argv[],
                                              Integer* argument,
                                              Real* centralLongitude,
                                              Real* centralLatitude,
                                              Integer* ok ) {

  PRE09( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ], argument,
         IN_RANGE( *argument, 1, argc ),
         centralLongitude, centralLatitude, ok );

  const Integer failures = failureCount();
  *ok = 0;
  *centralLongitude = *centralLatitude = 0.0;

  if ( ! AND2( *argument + 2 <= argc, argv[ *argument ] ) ) {
    failureMessage( "Invalid central-longitude-latitude command-line"
                    " argument '%s'.",
                    argv[ *argument ] );
  } else {
    *centralLongitude = toReal( argv[ *argument ], -180.0, 180.0, ok );

    if ( *ok ) {
      ++*argument;
      *centralLatitude = toReal( argv[ *argument ], -90.0, 90.0, ok );
      *argument += *ok;
    }
  }

  if ( ! *ok ) {

    if ( failures == failureCount() ) {
      failureMessage( "Invalid argument '%s'.", argv[ *argument ] );
    }

    *centralLongitude = *centralLatitude = 0.0;
  }

  POST03( IS_BOOL( *ok ),
          IN_RANGE( *argument, 0, argc ),
          IMPLIES_ELSE( *ok,
                        isValidLongitudeLatitude( *centralLongitude,
                                                  *centralLatitude ),
                        IS_ZERO2( *centralLongitude, *centralLatitude ) ) );
}


/******************************************************************************
PURPOSE: parseLayers - Parse layers command-line arguments into
         parameters.
INPUTS:  int argc           Number of command-line arguments.
         char* argv[]       Command-line argument strings.
         Integer* argument  Index of argument to parse.
OUTPUTS: Integer* argument  Updated index of argument to parse.
         Real* g            Gravitational force, e.g., 9.81 m/s^2.
         Real* R            Gas constant e.g., 287.04 J/kg/K = m^3/s/K.
         Real* A            Atmospheric lapse rate, e.g., 50.0 K/kg.
         Real* T0s          Reference surface temperature, e.g., 290.0 K.
         Real* P00          Reference surface pressure, e.g., 100000 P.
         Integer* layers    Number of layers of grid cells.
         Integer* type      Vertical grid type: e.g., VGSGPN3.
         Real* topPressure  Pressure, in Pascals, at the top of the model.
RETURNS: Real* levels       Allocated array of levels around layers.
******************************************************************************/

static Real* parseLayers( int argc, char* argv[], Integer* argument,
                          Real* g, Real* R, Real* A, Real* T0s, Real* P00,
                          Integer* layers, Integer* type, Real* topPressure) {

  PRE015( argc > 0, argv, argv[ 0 ], argv[ argc - 1 ],
          argument, *argument > 0, *argument < argc,
          g, R, A, T0s, P00, layers, type, topPressure );

  const Integer failures = failureCount();
  Real* result = 0;
  Integer ok = 0;

  if ( ! AND2( *argument + 4 < argc,
               ! strcmp( argv[ *argument ], "-layers" ) ) ) {
    failureMessage( "Invalid -layers command-line argument '%s'.",
                    argv[ *argument ] );
  } else {
    ++*argument; /* Skip "-layers". */
    *layers = toInteger( argv[ *argument ], 1, argc - 5, &ok );

    if ( ok ) {
      ++*argument;
      *type = toInteger( argv[ *argument ], 1, VGWRFEM, &ok );

      if ( ok ) {
        CHECK( IS_VALID_VERTICAL_GRID_TYPE( *type ) );
        ++*argument;
        *topPressure = toReal( argv[ *argument ], 0.01, 1e8, &ok );
        ok = AND2( ok, *argument + *layers + 1 + 5 < argc );

        if ( ok ) {
          const Integer levels = *layers + 1;
          ++*argument;
          result = NEW_ZERO( Real, levels );
          ok = result != 0;

          if ( result ) {
            const Real minimumLevelDifference = 1e-6;
            Integer level = 0;

            do {
              Real minimum = 0.0;
              Real maximum =
                level > 0 ? result[ level - 1 ] - minimumLevelDifference : 1.0;

              if ( IN3( *type, VGZVAL3, VGHVAL3 ) ) {
                minimum =
                  level > 0 ? result[ level - 1 ] + minimumLevelDifference
                  : -1000.0;
                maximum = 100000.0;
              } else if ( *type == VGPRES3 ) {
                minimum =
                  level > 0 ? result[ level - 1 ] - minimumLevelDifference
                  : 100.0;
                maximum = 10000.0;
              }

              ok = minimum < maximum;

              if ( ok ) {
                CHECK2( *argument < argc, IN_RANGE( level, 0, levels - 1 ) );
                result[ level ] =
                  toReal( argv[ *argument ], minimum, maximum, &ok );                
              }

              if ( ! ok ) {
                level = levels + 1;
              }

              *argument += ok;
              ++level;
            } while ( level < levels );

            if ( ok ) {
              *g = toReal( argv[ *argument ], 0.01, 1e2, &ok );
              *argument += ok;

              if ( ok ) {
                *R = toReal( argv[ *argument ], 0.01, 1e4, &ok );
                *argument += ok;

                if ( ok ) {
                  *A = toReal( argv[ *argument ], 0.01, 1e4, &ok );
                  *argument += ok;

                  if ( ok ) {
                    *T0s = toReal( argv[ *argument ], 0.01, 1e4, &ok );
                    *argument += ok;

                    if ( ok ) {
                      *P00 = toReal( argv[ *argument ], 0.01, 1e6, &ok );
                      *argument += ok;
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

  if ( ! ok ) {

    if ( failures == failureCount() ) {
      failureMessage( "Invalid argument '%s'.", argv[ *argument ] );
    }

    *layers = *type = 0;
    *topPressure = *g = *R = *A = *T0s = *P00 = 0.0;
    FREE( result );
  }

  POST08( *argument <= argc,
          ! isNan( *topPressure ),
          ! isNan( *g ),
          ! isNan( *R ),
          ! isNan( *A ),
          ! isNan( *T0s ),
          ! isNan( *P00 ),
          IMPLIES_ELSE( result,
                        AND9( *layers > 0,
                              IS_VALID_VERTICAL_GRID_TYPE( *type ),
                              IN_RANGE( *topPressure, 0.01, 1e8 ),
                              IN_RANGE( *g, 0.01, 1e2 ),
                              IN_RANGE( *R, 0.01, 1e4 ),
                              IN_RANGE( *A, 0.01, 1e4 ),
                              IN_RANGE( *T0s, 0.01, 1e4 ),
                              IN_RANGE( *P00, 0.01, 1e6 ),
                              isNanFree( result, *layers + 1 ) ),
                        IS_ZERO8( *layers, *type, *topPressure,
                                  *g, *R, *A, *T0s, *P00 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: initializeCells - Zero-out and reinitialize array of cells.
INPUTS:  Integer count           Number of cells.
         Real minimumValidValue  The minimum valid value of data.
OUTPUTS: Cell cells[]            The zeroed/initialized cells.
******************************************************************************/

static void initializeCells( Integer count, Real minimumValidValue,
                             Cell cells[] ) {

  PRE03( count > 0, ! isNan( minimumValidValue ), cells );
  Integer index = 0;

  memset( cells, 0, count * sizeof (Cell) );

#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
    Cell* cell = cells + index;
    cell->minimumValidValue = minimumValidValue;
    omp_init_lock( &cell->lock );
  }

  POST02( IS_ZERO3( cells[ 0 ].column, cells[ 0 ].row, cells[ 0 ].count ),
          cells[ 0 ].minimumValidValue == minimumValidValue );
}



/******************************************************************************
PURPOSE: finalizeCells - Finalize array of cells.
INPUTS:  Integer count           Number of cells.
OUTPUTS: Cell cells[]            The zeroed/initialized cells.
******************************************************************************/

static void finalizeCells( Integer count, Cell cells[] ) {

  PRE02( count > 0, cells );
  Integer index = 0;


#pragma omp parallel for

  for ( index = 0; index < count; ++index ) {
#if defined( _OPENMP ) && ! defined( SERIAL_REGRID )
    Cell* cell = cells + index;
#endif
    omp_destroy_lock( &cell->lock );
  }

  memset( cells, 0, count * sizeof (Cell) );

  POST0( IS_ZERO3( cells[ 0 ].column, cells[ 0 ].row, cells[ 0 ].count ) );
}



/******************************************************************************
PURPOSE: zeroUnused - Zero-out portion of array beyond compact result.
INPUTS:  Integer start             0-based index into array to start clearing.
         Integer count             Size of arrays.
OUTPUTS: Integer columns[ count ]  1-based grid column numbers.
         Integer rows[ count ]     1-based grid row    numbers.
         Real longitudes[ count ]  Grid cell center longitudes.
         Real latitudes[ count ]   Grid cell center latitudes.
******************************************************************************/

static void zeroUnused( Integer start, Integer count,
                        Integer columns[], Integer rows[],
                        Real longitudes[], Real latitudes[] ) {

  PRE06( IN_RANGE( start, 0, count ), count > 0,
         columns, rows, longitudes, latitudes );

  Integer index = 0;

#pragma omp parallel for

  for ( index = start; index < count; ++index ) {
    columns[ index ] = rows[ index ] = 0;
    longitudes[ index ] = latitudes[ index ] = 0.0;
  }

  POST0( IMPLIES( start < count,
                  IS_ZERO4( columns[ count - 1 ], rows[ count - 1 ],
                            longitudes[ count - 1 ], latitudes[ count - 1 ])));
}



/******************************************************************************
PURPOSE: radiusSquared - Sum of squares of x, y and z.
INPUTS:  Real x  Normalized x-component to square.
         Real y  Normalized y-component to square.
         Real z  Normalized z-component to square.
RETURNS: Real sum of squares of x and y and z clamped to at least TOLERANCE.
******************************************************************************/

static Real radiusSquared( Real x, Real y, Real z ) {

  PRE03( ! isNan( x ), ! isNan( y ), ! isNan( z ) );

  Real result = SQUARE( x ) + SQUARE( y ) + SQUARE( z );

  if ( result < TOLERANCE ) {
    result = TOLERANCE;
  }

  POST03( ! isNan( result ), result >= TOLERANCE,
          OR2( aboutEqual( result, x * x + y * y + z * z ),
               result == TOLERANCE ) );
  return result;
}



/******************************************************************************
PURPOSE: commonPreAggregator - Initialize cell with data.
INPUTS:  Integer column      1-based grid column number.
         Integer row         1-based grid column number.
         Real gridLongitude  Grid cell center longitude.
         Real gridLatitude   Grid cell center latitude.
         Real unused_xOffset
         Real unused_yOffset
         Real unused_zOffset
         Real inputData      Initial data value for the cell.
         Real inputData2     Initial vector 2nd component value for the cell.
OUTPUTS: Cell* cell          Initialized cell.
******************************************************************************/

static void commonPreAggregator( Integer column, Integer row,
                                 Real gridLongitude, Real gridLatitude,
                                 Real unused_xOffset, Real unused_yOffset,
                                 Real unused_zOffset,
                                 Real inputData, Real inputData2, Cell* cell ) {

  PRE08( column > 0,
         row > 0,
         isValidLongitude( gridLongitude ),
         isValidLatitude( gridLatitude ),
         ! isNan( inputData ),
         cell,
         inputData >= cell->minimumValidValue,
         IS_ZERO3( cell->count, cell->column, cell->row ) );

  cell->count     = 1;
  cell->column    = column;
  cell->row       = row;
  cell->longitude = gridLongitude;
  cell->latitude  = gridLatitude;
  cell->data      = inputData;
  cell->data2     = inputData2;

  POST07( cell->count == 1,
          cell->column == column,
          cell->row == row,
          cell->longitude == gridLongitude,
          cell->latitude == gridLatitude,
          cell->data == inputData,
          cell->data2 == inputData2 );
}



/******************************************************************************
PURPOSE: nearestPreAggregator - Initialize cell with data for nearest method.
INPUTS:  Integer column      1-based grid column number.
         Integer row         1-based grid column number.
         Real gridLongitude  Grid cell center longitude.
         Real gridLatitude   Grid cell center latitude.
         Real xOffset        Normalized x-offset from grid cell center.
         Real yOffset        Normalized y-offset from grid cell center.
         Real zOffset        Normalized z-offset from grid cell center.
         Real inputData      Initial data value for the cell.
         Real inputData2     Initial vector 2nd component value for the cell.
OUTPUTS: Cell* cell          Initialized cell.
******************************************************************************/

static void nearestPreAggregator( Integer column, Integer row,
                                  Real gridLongitude, Real gridLatitude,
                                  Real xOffset, Real yOffset, Real zOffset,
                                  Real inputData, Real inputData2, Cell* cell) {

  PRE011( column > 0,
          row > 0,
          isValidLongitude( gridLongitude ),
          isValidLatitude( gridLatitude ),
          ! isNan( inputData ),
          cell,
          inputData >= cell->minimumValidValue,
          IS_ZERO3( cell->count, cell->column, cell->row ),
          IN_RANGE( xOffset, -1.0, 1.0 ),
          IN_RANGE( yOffset, -1.0, 1.0 ),
          IN_RANGE( zOffset, -1.0, 1.0 ) );

  commonPreAggregator( column, row, gridLongitude, gridLatitude,
                       xOffset, yOffset, zOffset, inputData, inputData2, cell );

  cell->radius = radiusSquared( xOffset, yOffset, zOffset );

  POST08( cell->count == 1,
          cell->column == column,
          cell->row == row,
          cell->longitude == gridLongitude,
          cell->latitude == gridLatitude,
          cell->data == inputData,
          cell->data2 == inputData2,
          cell->radius >= TOLERANCE );
}



/******************************************************************************
PURPOSE: nearestAggregator - Aggregate new data into cell for nearest method.
INPUTS:  Real xOffset    Normalized x-offset from grid cell center.
         Real yOffset    Normalized y-offset from grid cell center.
         Real zOffset    Normalized z-offset from grid cell center.
         Real inputData  Additional data value to aggregate into the cell.
         Real inputData2 Additional vector 2nd component value for the cell.
OUTPUTS: Cell* cell      Aggregated cell.
******************************************************************************/

static void nearestAggregator( Real xOffset, Real yOffset, Real zOffset,
                               Real inputData, Real inputData2, Cell* cell ) {

  PRE09( IN_RANGE( xOffset, -1.0, 1.0 ),
         IN_RANGE( yOffset, -1.0, 1.0 ),
         IN_RANGE( zOffset, -1.0, 1.0 ),
         ! isNan( inputData ),
         cell,
         inputData >= cell->minimumValidValue,
         GT_ZERO3( cell->count, cell->column, cell->row ),
         isValidLongitude( cell->longitude ),
         isValidLatitude( cell->latitude ) );

  const Real minimumValidValue = cell->minimumValidValue;

  if ( inputData >= minimumValidValue ) {
    const Real radius = radiusSquared( xOffset, yOffset, zOffset );
    const Real cellData = cell->data;
    const Real cellRadius = cell->radius;

    if ( OR2( cellData < minimumValidValue, radius < cellRadius ) ) {
      cell->radius = radius;
      cell->data = inputData;
      cell->data2 = inputData2;
    }
  }

  POST0( ! isNan( cell->data ) );
}



/******************************************************************************
PURPOSE: meanAggregator - Aggregate new data into cell for mean method.
INPUTS:  Real unused_xOffset
         Real unused_yOffset
         Real unused_zOffset
         Real inputData  Additional data value to aggregate into the cell.
         Real inputData2 Additional vector 2nd component value for the cell.
OUTPUTS: Cell* cell      Aggregated cell.
******************************************************************************/

static void meanAggregator( Real unused_xOffset, Real unused_yOffset,
                            Real unused_zOffset, Real inputData,
                            Real inputData2, Cell* cell ) {

  PRE06( ! isNan( inputData ),
         cell,
         inputData >= cell->minimumValidValue,
         GT_ZERO3( cell->count, cell->column, cell->row ),
         isValidLongitude( cell->longitude ),
         isValidLatitude( cell->latitude ) );

  const Real minimumValidValue = cell->minimumValidValue;

  if ( inputData >= minimumValidValue ) {
    Real cellData = cell->data;
    Real cellData2 = cell->data2;

    if ( cellData < minimumValidValue ) {
      cell->data = inputData; /* Replace invalid value. */
      cell->data2 = inputData2; /* Replace invalid value. */
    } else { /* Update mean: */
      const Integer cellCount = cell->count;
      cellData *= cellCount;
      cellData += inputData;
      cellData /= cellCount + 1;
      cell->data = cellData;
      cellData2 *= cellCount;
      cellData2 += inputData2;
      cellData2 /= cellCount + 1;
      cell->data2 = cellData2;
      ++cell->count;
    }
  }

  POST0( ! isNan( cell->data ) );
}



/******************************************************************************
PURPOSE: weightedPreAggregator - Initialize cell with data for weighted method.
INPUTS:  Integer column      1-based grid column number.
         Integer row         1-based grid column number.
         Real gridLongitude  Grid cell center longitude.
         Real gridLatitude   Grid cell center latitude.
         Real xOffset        Normalized x-offset from grid cell center.
         Real yOffset        Normalized y-offset from grid cell center.
         Real zOffset        Normalized z-offset from grid cell center.
         Real inputData      Initial data value for the cell.
         Real inputData2     Initial vector 2nd component value for the cell.
OUTPUTS: Cell* cell          Initialized cell.
******************************************************************************/

static void weightedPreAggregator( Integer column, Integer row,
                                   Real gridLongitude, Real gridLatitude,
                                   Real xOffset, Real yOffset, Real zOffset,
                                   Real inputData, Real inputData2,
                                   Cell* cell ) {

  PRE011( column > 0,
          row > 0,
          isValidLongitude( gridLongitude ),
          isValidLatitude( gridLatitude ),
          IN_RANGE( xOffset, -1.0, 1.0 ),
          IN_RANGE( yOffset, -1.0, 1.0 ),
          IN_RANGE( zOffset, -1.0, 1.0 ),
          ! isNan( inputData ),
          cell,
          inputData >= cell->minimumValidValue,
          IS_ZERO3( cell->count, cell->column, cell->row ) );

  const Real radius = radiusSquared( xOffset, yOffset, zOffset );
  const Real weight = 1.0 / radius;

  commonPreAggregator( column, row, gridLongitude, gridLatitude,
                       xOffset, yOffset, zOffset, inputData, inputData2, cell );

  cell->radius = radius;

  if ( cell->data >= cell->minimumValidValue ) {
    cell->data *= weight;
    cell->data2 *= weight;
    cell->weights = weight;
  }

  POST08( cell->count == 1,
          cell->column == column,
          cell->row == row,
          cell->longitude == gridLongitude,
          cell->latitude == gridLatitude,
          cell->radius >= TOLERANCE,
          cell->weights >= 0.0,
          ! isNan( cell->data ) );
}



/******************************************************************************
PURPOSE: weightedAggregator - Aggregate new data into cell for weighted method.
INPUTS:  Real xOffset    Normalized x-offset from grid cell center.
         Real yOffset    Normalized y-offset from grid cell center.
         Real zOffset    Normalized z-offset from grid cell center.
         Real inputData  Additional data value to aggregate into the cell.
         Real inputData2 Additional vector 2nd component value for the cell.
OUTPUTS: Cell* cell      Aggregated cell.
******************************************************************************/

static void weightedAggregator( Real xOffset, Real yOffset, Real zOffset,
                                Real inputData, Real inputData2, Cell* cell ) {

  PRE09( IN_RANGE( xOffset, -1.0, 1.0 ),
         IN_RANGE( yOffset, -1.0, 1.0 ),
         IN_RANGE( zOffset, -1.0, 1.0 ),
         ! isNan( inputData ),
         cell,
         inputData >= cell->minimumValidValue,
         GT_ZERO3( cell->count, cell->column, cell->row ),
         isValidLongitude( cell->longitude ),
         isValidLatitude( cell->latitude ) );

  const Real minimumValidValue = cell->minimumValidValue;

  if ( inputData >= minimumValidValue ) {
    const Real radius = radiusSquared( xOffset, yOffset, zOffset );
    const Real weight = 1.0 / radius;

    if ( cell->data < minimumValidValue ) { /* Replace invalid data: */
      cell->weights = weight;
      cell->data = inputData * weight;
      cell->data2 = inputData2 * weight;
    } else { /* Accumulate valid data: */
      cell->weights += weight;
      cell->data += inputData * weight;
      cell->data2 += inputData2 * weight;
      ++cell->count;
    }
  }

  POST02( ! isNan( cell->data ), cell->weights >= 0.0 );
}



/******************************************************************************
PURPOSE: weightedPostAggregator - Unscale aggregated cell data.
INPUTS:  Cell* cell  Scaled   aggregated cell.
OUTPUTS: Cell* cell  Unscaled aggregated cell.
******************************************************************************/

static void weightedPostAggregator( Cell* cell ) {

  PRE03( cell, cell->weights >= 0.0, ! isNan( cell->data ) );

  if ( cell->weights > 0.0 ) {
    cell->data /= cell->weights;
    cell->data2 /= cell->weights;
  }

  POST0( ! isNan( cell->data ) );
}



/******************************************************************************
PURPOSE: aggregateCellData - Initialize and/or aggregate cell data.
INPUTS:  Grid* self                         Grid object.
         PreAggregator preAggregator        Pre-aggregator routine.
         Aggregator aggregator              Aggregator routine.
         Integer column                     1-based index of grid cell.
         Integer row                        1-based index of grid cell.
         Real gridLongitude                 Longitude of center of grid cell.
         Real gridLatitude                  Latitude  of center of grid cell.
         Real xOffset                       Normalized offset from cell center.
         Real yOffset                       Normalized offset from cell center.
         Integer point                      0-based index of point data.
         const Real data[ points * levels ] Data of the point.
         const Real data2[ points * levels ] 0 or vector 2nd component data
                                             of the point.
         Integer levels                     Number of levels of data.
         const Real elevations[ points * levels ] Elevations of the point data.
         const Note note                    Optional note or 0.
OUTPUTS: Cell cells[ grid->layers( grid ) ] Initialized cell data.
******************************************************************************/

static void aggregateCellData( Grid* self,
                               PreAggregator preAggregator,
                               Aggregator aggregator,
                               Integer column, Integer row,
                               Real gridLongitude, Real gridLatitude,
                               Real xOffset, Real yOffset,
                               Integer point, const Real data[],
                               const Real data2[],
                               Integer levels, const Real elevations[],
                               const Note note,
                               Cell cells[] ) {

  PRE14( preAggregator,
         aggregator,
         IN_RANGE( column, 1, self->columns( self ) ),
         IN_RANGE( row,    1, self->rows( self ) ),
         isValidLongitude( gridLongitude ),
         isValidLatitude( gridLatitude ),
         IN_RANGE( xOffset, -1.0, 1.0 ),
         IN_RANGE( yOffset, -1.0, 1.0 ),
         point >= 0,
         data,
         levels > 0,
         IMPLIES( levels > 1, elevations ),
         IMPLIES( elevations, isNanFree( elevations, levels ) ),
         cells );

  const GridPrivate* const selfData = self->data;
  const Real minimumValidValue = cells[ 0 ].minimumValidValue;

  if ( ! elevations ) {
    const Real zOffset = 0.0;
    Real dataValue2 = 0.0;
    const Real dataValue =
      surfacePointData( point, data, data2, levels, elevations, &dataValue2 );

    if ( dataValue >= minimumValidValue ) {

      if ( cells[ 0 ].count == 0 ) {
        preAggregator( column, row, gridLongitude, gridLatitude,
                       xOffset, yOffset, zOffset, dataValue, dataValue2,
                       &cells[ 0 ] );
      } else {
        aggregator( xOffset, yOffset, zOffset, dataValue, dataValue2,
                    &cells[ 0 ] );
      }

      if ( note ) {
        appendNote( cells[ 0 ].regriddedNote, note );
      }
    }
  } else {

    /*
     * Compute cellElevations[ gridLevels ] to contain the elevation in meters
     * above mean sea level of the boundaries between the vertical stack of
     * grid cells at this row-column containing the data point.
     * Note that it is based on the point's elevation (e.g., CALIPSO LIDAR
     * measured) rather than the 12km averaged grid cell terrain height.
     */

    const Integer gridLayers          = selfData->layers;
    const Integer gridLevels          = gridLayers + 1;
    const Integer thread              = omp_get_thread_num();
    Real* const cellElevations        = selfData->z + thread * gridLevels;
    const Real* const sigmaPressures  = selfData->levels;
    const Integer pointOffset         = point * levels;
    const Real* const pointElevations = elevations + pointOffset;
    const Real* const pointData       = data       + pointOffset;
    const Real* const pointData2      = data2 ? data2 + pointOffset : 0;
    const Integer multiLevel          = levels > 1;
    const Integer surfaceLayer =
      multiLevel ? surfaceIndex( levels, pointElevations) : 0;
    extern float elevationAt( float, float ); /* See file ./elevation.c */
    const Real surfaceElevation0 =
      multiLevel ? pointElevations[ surfaceLayer ]
      : elevationAt( gridLongitude, gridLatitude );
    const Real surfaceElevation =
      AND2( ! multiLevel, surfaceElevation0 < 0.0 ) ? 0.0 : surfaceElevation0;
    const Integer dataLayers          = levels;
    Integer dataLayer                 = surfaceLayer;
    Integer previousCellIndex         = 0;
    Real previousDataElevation        = surfaceElevation;
    const Real minimumElevationDifference = 40.0; /* If > then recompute Z. */

    /* (Re)initialize cell center elevations and layers: */

    if ( OR2( cells[ 0 ].layer == 0,
              fabs( surfaceElevation - cells[ 0 ].surfaceElevation )
              > minimumElevationDifference ) ) {
      Integer gridLayer = 0;
      cells[ 0 ].surfaceElevation = surfaceElevation;

      if ( IN4( selfData->type, VGSGPH3, VGSGPN3, VGWRFEM ) ) {
        DEBUG( fprintf( stderr, "calling elevationAtSigmaPressues()...\n" ); )
        elevationsAtSigmaPressures( selfData->g, selfData->R, selfData->A,
                                    selfData->T0s, selfData->P00,
                                    surfaceElevation, gridLevels,
                                    selfData->topPressure,
                                    sigmaPressures, cellElevations );
        DEBUG( fprintf( stderr, "selfData->g = %lf, selfData->R = %lf, "
                                "selfData->A = %lf, selfData->T0s = %lf, "
                                "selfData->P00 = %lf, "
                                "surfaceElevation = %lf, gridLevels = %lld, "
                                "selfData->topPressure = %lf, "
                                "sigmaPressures = [%lf %lf ...], "
                                "cellElevations = [%lf %lf ...]\n",
                                selfData->g, selfData->R, selfData->A,
                                selfData->T0s, selfData->P00,
                                surfaceElevation, gridLevels,
                                selfData->topPressure,
                                sigmaPressures[ 0 ],
                                sigmaPressures[ 1 ],
                                cellElevations[ 0 ],
                                cellElevations[ 1 ] ); )
      } else if ( selfData->type == VGZVAL3 ) {
        Integer level = 0;

        for ( level = 0; level < gridLevels; ++level ) {
          cellElevations[ level ] = selfData->levels[ level ];
        }

      } else if ( selfData->type == VGHVAL3 ) {
        Integer level = 0;

        for ( level = 0; level < gridLevels; ++level ) {
          cellElevations[ level ] = surfaceElevation + selfData->levels[level];
        }

      } else { /* HACK: Other vertical grid schemes are UNIMPLEMENTED! */
        Integer level = 0;

        for ( level = 0; level < gridLevels; ++level ) {
          cellElevations[ level ] = selfData->levels[ level ];
        }

      }

      /* Compute cell center as average of cell lower and upper edges: */

      for ( gridLayer = 0; gridLayer < gridLayers; ++gridLayer ) {
        const Real gridCellCenterElevation =
          0.5 * ( cellElevations[ gridLayer ] + cellElevations[gridLayer + 1]);
        cells[ gridLayer ].elevation = gridCellCenterElevation;
        cells[ gridLayer ].layer = gridLayer + 1;
      }
    }

    /* Initialize or aggregate each above-ground value into grid cells: */

    for ( dataLayer = surfaceLayer; dataLayer < dataLayers; ++dataLayer ) {
      const Real dataElevation = pointElevations[ dataLayer ];
      Real zOffset = 0.0;
      Integer cellIndex  =
        binElevation( dataElevation, gridLevels, cellElevations,
                      previousCellIndex, &zOffset );

      DEBUG( fprintf( stderr, "dataElevation = %lf, cellIndex = %lld\n",
                      dataElevation, cellIndex ); )

      previousCellIndex =
        AND2( cellIndex != -1, dataElevation >= previousDataElevation ) ?
          cellIndex
        : 0;
      previousDataElevation =
        previousCellIndex ? dataElevation : surfaceElevation;

      /* HACK: force data points below grid layer 1 into layer 1: */

      if ( AND4( ! multiLevel, cellIndex == -1, dataElevation >= 0.0,
                 dataElevation < cellElevations[ 0 ] ) ) {
        cellIndex = 0;  /* Pretend the data point is in layer 1. */
        zOffset = -1.0; /* At the bottom edge of layer 1. */
      }

      if ( cellIndex != -1 ) {
        const Real dataValue = pointData[ dataLayer ];
        const Real dataValue2 = pointData2 ? pointData2[ dataLayer ] : 0.0;
        DEBUG( fprintf( stderr, "dataValue = %lf, dataValue2 = %lf\n",
                        dataValue, dataValue2 ); )

        if ( dataValue >= minimumValidValue ) {

          if ( cells[ cellIndex ].count == 0 ) {
            preAggregator( column, row, gridLongitude, gridLatitude,
                           xOffset, yOffset, zOffset, dataValue, dataValue2,
                           &cells[ cellIndex ] );
          } else {
            aggregator( xOffset, yOffset, zOffset, dataValue, dataValue2,
                        &cells[ cellIndex ] );
          }

          DEBUG( fprintf( stderr, "cells[ cellIndex %lld ].count = %lld\n",
                          cellIndex, cells[ cellIndex ].count ); )

          if ( note ) {
            appendNote( cells[ cellIndex ].regriddedNote, note );
          }
        }
      }
    }
  }

  POST( ! isNan( cells[ 0 ].data ) );
}



/******************************************************************************
PURPOSE: surfacePointData - Obtain the point data at the surface.
INPUTS:  Integer point                       0-based index of point data.
         const Real data[ points * layers ]  Data of the point.
         const Real data2[ points * layers ] 0 or vector 2nd component data of
                                             the point.
         Integer layers                      Number of layers of data.
         const Real elevations[ points * layers ] Elevations of the point data.
OUTPUTS: Real* dataValue2                   0 or vector 2nd component at point.
RETURNS: Real the data at the surface (after skipping any possibly collapsed
         elevations which indicate sub-surface points, e.g., CALIPSO).
******************************************************************************/

static Real surfacePointData( Integer point, const Real data[],
                              const Real data2[],
                              Integer layers, const Real elevations[],
                              Real* dataValue2 ) {

  PRE04( point >= 0, data, layers > 0, IMPLIES( layers > 1, elevations ) );

  Real result = data[ point ];

  if ( AND2( data2, dataValue2 ) ) {
    *dataValue2 = data2[ point ];
  }

  if ( layers > 1 ) { /* Skip over points with collapsed elevations: */
    const Integer pointOffset  = point * layers;
    const Real* const pointElevations = elevations + pointOffset;
    const Integer surfaceLayer = surfaceIndex( layers, pointElevations );
    result = data[ pointOffset + surfaceLayer ];

    if ( AND2( data2, dataValue2 ) ) {
      *dataValue2 = data2[ pointOffset + surfaceLayer ];
    }
  }

  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: surfaceIndex - 0-based index of surface elevation.
INPUTS:  Integer layers                      Number of layers of data.
         const Real elevations[ layers ]     Elevations of the point data.
RETURNS: Integer 0-based index of the surface (after skipping any collapsed
         elevations which indicate sub-surface points, e.g., CALIPSO).
******************************************************************************/

static Integer surfaceIndex( Integer layers, const Real elevations[] ) {

  PRE02( layers > 0, IMPLIES( layers > 1, elevations ) );

  Integer result = 0;

  if ( layers > 1 ) {
    Real surfaceZ = elevations[ 0 ];
    Integer index = 1;

    do {
      const Real z = elevations[ index ];

      if ( ! aboutEqual( surfaceZ, z ) ) {
        result = index - 1;
        index = layers;
      } else {
        surfaceZ = z;
      }

      ++index;
    } while ( index < layers );
  }

  POST0( IN_RANGE( result, 0, layers - 1 ) );
  return result;
}



/******************************************************************************
PURPOSE: binElevation - 0-based index of grid cell containing elevation.
INPUTS:  Real dataElevation      Elevation of a data point to bin.
         Integer gridLevels      Number of levels of grid cell boundaries.
         const Real cellElevations[ gridLevels ]  Elevations at gridLevels.
         Integer startingLayer   Starting grid layer index to search.
OUTPUTS: Real* zOffset           Normalized offset from grid cell center.
RETURNS: Integer 0-based grid layer index of grid cell containing elevation
         or -1 if outside the grid.
******************************************************************************/

static Integer binElevation( Real dataElevation, Integer gridLevels,
                             const Real cellElevations[],
                             Integer startingLayer, Real* zOffset ) {

  PRE07( dataElevation > -1000.0,
         gridLevels > 1,
         cellElevations,
         isNanFree( cellElevations, gridLevels ),
         increasing( cellElevations, gridLevels ),
         IN_RANGE( startingLayer, 0, gridLevels - 2 ),
         zOffset );

  const Integer gridLayers = gridLevels - 1;
  Integer result = -1;
  Integer layer = startingLayer;
  Real lower = cellElevations[ layer ];
  *zOffset = 0.0;

  for ( layer = startingLayer; layer < gridLayers; ++layer ) {
    const Real upper = cellElevations[ layer + 1 ];
    CHECK( upper > lower );

    if ( IN_RANGE( dataElevation, lower, upper ) ) {

      /* Compute normalized offset from geometric grid cell center: */

      Real z = ( dataElevation - lower ) / ( upper - lower ); /* [0, 1] */
      z += z;   /* [0, 2] */
      z -= 1.0; /* [-1, 1] */
      *zOffset = z;
      result = layer;
      layer  = gridLayers; /* Stop looping. */
    } else {
      lower = upper;
    }
  }

  POST0( IMPLIES_ELSE( result == -1,
                       *zOffset == 0.0,
                       AND2( IN_RANGE( result, 0, gridLevels - 1 ),
                             IN_RANGE( *zOffset, -1.0, 1.0 ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: elevationsAtSigmaPressures - Compute elevations in meters above mean
         sea-level at sigma-pressures.
INPUTS:  Real g                Gravitational force, e.g., 9.81 m/s^2.
         Real R                Gas constant e.g., 287.04 J/kg/K = m^3/s/K.
         Real A                Atmospheric lapse rate, e.g., 50.0 K/kg.
         Real T0s              Reference surface temperature, e.g., 290.0 K.
         Real P00              Reference surface pressure, e.g., 100000 P.
         Real surfaceElevation     Elevation of surface in meters AMSL.
         Integer levels            Number of levels of sigmaPressures.
         Real topPressure          Pressure in Pascals at the top of the model.
         const Real sigmaPressures[ levels ]  Sigma-pressures at levels.
OUTPUTS: Real elevations[ levels ]  Elevation in meters above MSL at sigmas.
NOTES:   Based on formula used in MM5.
******************************************************************************/

static void elevationsAtSigmaPressures( Real g, Real R, Real A, Real T0s,
                                        Real P00,
                                        Real surfaceElevation,
                                        Integer levels, Real topPressure,
                                        const Real sigmaPressures[],
                                        Real elevations[] ) {

  PRE015( ! isNan( g ),
          ! isNan( R ),
          ! isNan( A ),
          ! isNan( T0s ),
          ! isNan( P00 ),
          ! isNan( surfaceElevation ),
          surfaceElevation > -1000.0,
          levels > 0,
          ! isNan( topPressure ),
          GT_ZERO6( topPressure, g, R, A, T0s, P00 ),
          isNanFree( sigmaPressures, levels ),
          decreasing( sigmaPressures, levels ),
          minimumItem( sigmaPressures, levels ) >= 0.0,
          maximumItem( sigmaPressures, levels ) <= 1.0,
          elevations );

  /* Derived constants: */

  const Real H0s            = R * T0s / g;
  const Real one_over_H0s   = 1.0 / H0s;
  const Real A_over_T0s     = A / T0s;
  const Real A_over_two_T0s = A / ( T0s + T0s );
  const Real Pt             = topPressure;
  const Real Zs             = surfaceElevation;
  const Real two_Zs         = Zs + Zs;
  const Real sqrt_factor    = sqrt( 1.0 - A_over_T0s * one_over_H0s * two_Zs );
  const Real q_factor =
    ( Pt / P00 ) * exp( two_Zs * one_over_H0s / sqrt_factor );
  Integer level = 0;

  /* Compute elevations at sigma-pressures: */

  for ( level = 0; level < levels; ++level ) {
    const Real sigma_p0   = sigmaPressures[ level ];
    const Real q0_star    = sigma_p0 + ( 1.0 - sigma_p0 ) * q_factor;
    const Real ln_q0_star = log( q0_star );
    const Real z_level    =
      Zs - H0s * ln_q0_star * ( A_over_two_T0s * ln_q0_star + sqrt_factor );
    elevations[ level ]   = z_level;
  }

  POST04( isNanFree( elevations, levels ),
          increasing( elevations, levels ),
          minimumItem( elevations, levels ) >= -1000.0,
          maximumItem( elevations, levels ) <= 1e6 );

}






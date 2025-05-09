/******************************************************************************
PURPOSE: RegridQuadrilaterals.h - Declares routine to regrid 2D quadrilaterals
         onto a regular 2D grid with or without area weighting
         (i.e., scaling the data by the fractional area of the clipped polygon).

HISTORY: 2019-12-23 plessel.todd@epa.gov
STATUS:  unreviewed tested
NOTES:
         Uses the Liang-Barsky polygon clipping algorithm. (Fastest known.)
         "An Analysis and Algorithm for Polygon Clipping",
         You-Dong Liang and Brian Barsky, UC Berkeley,
         CACM Vol 26 No. 11, November 1983.
         https://www.longsteve.com/fixmybugs/?page_id=210
         Also:
http://geomalgorithms.com/a01-_area.html
http://geomalgorithms.com/a03-_inclusion.html
http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
https://demonstrations.wolfram.com/PointInTriangle/
https://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
https://en.wikipedia.org/wiki/Point_in_polygon
https://stackoverflow.com/questions/15490795/determine-if-a-2d-point-is-within-a-quadrilateral
******************************************************************************/


/*================================ INCLUDES =================================*/

#include <float.h>  /* For DBL_MAX */
#include <assert.h> /* For macro assert(). */
#ifndef NDEBUG
#include <stdio.h>  /* For stderr, fprintf(). */
#endif

#include <RegridQuadrilaterals.h>  /* For public interface. */

/*================================= MACROS ==================================*/

#ifdef NDEBUG
#define DEBUG(unused)
#else
#define DEBUG(s) s
#endif

#define IN_RANGE(x,low,high) ((low)<=(x)&&(x)<=(high))

#define IS_LONGITUDE( value ) ( IN_RANGE( (value), -180.0, 180.0 ) )
#define IS_LATITUDE(  value ) ( IN_RANGE( (value),  -90.0,  90.0 ) )
#define IS_LONGITUDE_LATITUDE( longitude, latitude ) \
( IS_LONGITUDE( longitude ) && IS_LATITUDE( latitude ) )

#define CLAMPED_TO_RANGE( value, low, high ) \
((value) < (low) ? (low) : (value) > (high) ? (high) : (value))

/*
 * Compute 2D vector cross-product:
 * (from1X, from1Y)->(to1X, to1Y) X (from2X, from2Y)->(to2X, to2Y)
 */

#define CROSS2D( from1X, from1Y, to1X, to1Y, from2X, from2Y, to2X, to2Y ) \
  ( ( (to1X) - (from1X) ) * ( (to2Y) - (from2Y) ) - \
    ( (to2X) - (from2X) ) * ( (to1Y) - (from1Y) ) )

/*================================== TYPES ==================================*/

enum { X, Y };
enum { MINIMUM, MAXIMUM };

typedef double Bounds[ 2 ][ 2 ]; /* bounds[ X, Y ][ MINIMUM, MAXIMUM ]. */

/*========================== FORWARD DECLARATIONS ===========================*/

static int regridQuadrilateral( const double data,
                                const double x[],
                                const double y[],
                                const size_t rows,
                                const size_t columns,
                                const double gridXMinimum,
                                const double gridYMinimum,
                                const double cellWidth,
                                const double cellHeight,
                                size_t cellCounts[],
                                double cellWeights[],
                                double cellSums[] );

static double areaOfClippedQuadrilateral( const double cellXMinimum,
                                          const double cellYMinimum,
                                          const double cellXMaximum,
                                          const double cellYMaximum,
                                          const double x[],
                                          const double y[] );

static void computePolygonBounds( const size_t count,
                                  const double x[],
                                  const double y[],
                                  Bounds bounds );

static int boundsOverlap( const Bounds bounds1, const Bounds bounds2 );

#if 0
/* Unused unless weighting by grid cell is used. */
static int pointIsInsidePolygon( const double x, const double y,
                                 const size_t vertexCount,
                                 const double vx[], const double vy[] );

static double cross2d( const double from1X, const double from1Y,
                       const double to1X,   const double to1Y,
                       const double from2X, const double from2Y,
                       const double to2X,   const double to2Y );
#endif

static double areaOfQuadrilateral( const double x1, const double y1,
                                   const double x2, const double y2,
                                   const double x3, const double y3,
                                   const double x4, const double y4 );

static double areaOfTriangle( const double x1, const double y1,
                              const double x2, const double y2,
                              const double x3, const double y3 );

static double signedAreaOfPolygon( const size_t count,
                                   const double x[], const double y[] );

static int clipPolygon( const int discardDegenerates,
                        const double clipXMin,
                        const double clipYMin,
                        const double clipXMax,
                        const double clipYMax,
                        const int count,
                        const double x[],
                        const double y[],
                        double clippedX[],
                        double clippedY[] );


/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: projectAndOrReorderQuadrilateralVertices - Project quadrilateral
         vertices (if project function provded) and
         reorder as contiguous counter-clockwise arrays.
INPUTS:  const size_t count  Number of quadrilaterals.
         const double longitudesSW[ count ]  SW corner longitudes.
         const double longitudesSE[ count ]  SE corner longitudes.
         const double longitudesNW[ count ]  NW corner longitudes.
         const double longitudesNE[ count ]  NE corner longitudes.
         const double latitudesSW[  count ]  SW corner latitudes.
         const double latitudesSE[  count ]  SE corner latitudes.
         const double latitudesNW[  count ]  NW corner latitudes.
         const double latitudesNE[  count ]  NE corner latitudes.
         void* projector                     Projector object or 0 if lonlat.
         ProjectFunction project             Project function or 0 if lonlat.
OUTPUTS: double vx[ count * 4 ]              Projected ccw x-vertices per quad.
OUTPUTS: double vy[ count * 4 ]              Projected ccw y-vertices per quad.
******************************************************************************/

void projectAndOrReorderQuadrilateralVertices( const size_t count,
                                               const double longitudesSW[],
                                               const double longitudesSE[],
                                               const double longitudesNW[],
                                               const double longitudesNE[],
                                               const double latitudesSW[],
                                               const double latitudesSE[],
                                               const double latitudesNW[],
                                               const double latitudesNE[],
                                               void* projector,
                                               ProjectFunction project,
                                               double vx[],
                                               double vy[] ) {

  const long long count0 = (long long) count; /* Must use signed type for omp*/

  assert( count );
  assert( longitudesSW );
  assert( longitudesSE );
  assert( longitudesNW );
  assert( longitudesNE );
  assert( latitudesSW );
  assert( latitudesSE );
  assert( latitudesNW );
  assert( latitudesNE );
  assert( vx );
  assert( vy );
  assert( ! projector || project ); /* projector implies project function. */

  if ( ! project ) { /* Just copy and reorder ccw SW,SE,NE,NW per quad: */
    long long index = 0;

#pragma omp parallel
    { /* Start of parallel region: */

#pragma omp for

    for ( index = 0; index < count0; ++index ) {
      const long long index2 = index + index;
      const long long index4 = index2 + index2;
      assert( IS_LONGITUDE( longitudesSW[ index ] ) );
      assert( IS_LONGITUDE( longitudesSE[ index ] ) );
      assert( IS_LONGITUDE( longitudesNE[ index ] ) );
      assert( IS_LONGITUDE( longitudesNW[ index ] ) );
      vx[ index4     ] = longitudesSW[ index ];
      vx[ index4 + 1 ] = longitudesSE[ index ];
      vx[ index4 + 2 ] = longitudesNE[ index ];
      vx[ index4 + 3 ] = longitudesNW[ index ];
    }

#pragma omp for

    for ( index = 0; index < count0; ++index ) {
      const long long index2 = index + index;
      const long long index4 = index2 + index2;
      assert( IS_LATITUDE( latitudesSW[ index ] ) );
      assert( IS_LATITUDE( latitudesSE[ index ] ) );
      assert( IS_LATITUDE( latitudesNE[ index ] ) );
      assert( IS_LATITUDE( latitudesNW[ index ] ) );
      vy[ index4     ] = latitudesSW[  index ];
      vy[ index4 + 1 ] = latitudesSE[  index ];
      vy[ index4 + 2 ] = latitudesNE[  index ];
      vy[ index4 + 3 ] = latitudesNW[  index ];
    }

    } /* End of parallel region. */

  } else { /* Project (longitude, latitude) to (x, y) and store ccw sequence:*/
    long long index = 0;

#pragma omp parallel for

    for ( index = 0; index < count0; ++index ) {
      const long long index2 = index + index;
      const long long index4_0 = index2 + index2;
      const long long index4_1 = index4_0 + 1;
      const long long index4_2 = index4_0 + 2;
      const long long index4_3 = index4_0 + 3;
      assert( IS_LONGITUDE_LATITUDE( longitudesSW[index], latitudesSW[index]));
      assert( IS_LONGITUDE_LATITUDE( longitudesSE[index], latitudesSE[index]));
      assert( IS_LONGITUDE_LATITUDE( longitudesNE[index], latitudesNE[index]));
      assert( IS_LONGITUDE_LATITUDE( longitudesNW[index], latitudesNW[index]));
      project( projector, longitudesSW[ index ], latitudesSW[ index ],
               vx + index4_0, vy + index4_0 );
      project( projector, longitudesSE[ index ], latitudesSE[ index ],
               vx + index4_1, vy + index4_1 );
      project( projector, longitudesNE[ index ], latitudesNE[ index ],
               vx + index4_2, vy + index4_2 );
      project( projector, longitudesNW[ index ], latitudesNW[ index ],
               vx + index4_3, vy + index4_3 );
      assert( IN_RANGE( vx[ index4_0 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_0 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vx[ index4_1 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_1 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vx[ index4_2 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_2 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vx[ index4_3 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_3 ], -DBL_MAX, DBL_MAX ) );
    }
  }
}



/******************************************************************************
PURPOSE: projectAndOrReorderQuadrilateralVertices2 - Project quadrilateral
         vertices (if project function provded) and
         reorder as contiguous counter-clockwise arrays.
INPUTS:  const size_t count  Number of quadrilaterals.
         const double longitudes[ count ]    Cell center longitudes.
         const double latitudes[  count ]    Cell center latitudes.
         const double cellWidth              Cell width in degrees.
         const double cellHeight             Cell height in degrees.
         void* projector                     Projector object or 0 if lonlat.
         ProjectFunction project             Project function or 0 if lonlat.
OUTPUTS: double vx[ count * 4 ]              Projected ccw x-vertices per quad.
OUTPUTS: double vy[ count * 4 ]              Projected ccw y-vertices per quad.
******************************************************************************/

void projectAndOrReorderQuadrilateralVertices2( const size_t count,
                                                const double longitudes[],
                                                const double latitudes[],
                                                const double cellWidth,
                                                const double cellHeight,
                                                void* projector,
                                                ProjectFunction project,
                                                double vx[],
                                                double vy[] ) {

  const long long count0 = (long long) count; /* Must use signed type for omp*/

  assert( count );
  assert( longitudes );
  assert( latitudes );
  assert( cellWidth > 0.0 );
  assert( cellWidth <= 180.0 );
  assert( cellHeight > 0.0 );
  assert( cellHeight <= 180.0 );
  assert( vx );
  assert( vy );
  assert( ! projector || project ); /* projector implies project function. */

  if ( ! project ) { /* Just copy and reorder ccw SW,SE,NE,NW per quad: */
    long long index = 0;

#pragma omp parallel
    { /* Start of parallel region: */

#pragma omp for

    for ( index = 0; index < count0; ++index ) {
      const long long index2 = index + index;
      const long long index4 = index2 + index2;
      const double longitude = longitudes [ index ];
      const double longitudeW = longitude - cellWidth;
      const double longitudeE = longitude + cellWidth;
      assert( IS_LONGITUDE( longitude ) );
      assert( IS_LONGITUDE( longitudeW ) );
      assert( IS_LONGITUDE( longitudeE ) );
      vx[ index4     ] = longitudeW;
      vx[ index4 + 1 ] = longitudeE;
      vx[ index4 + 2 ] = longitudeE;
      vx[ index4 + 3 ] = longitudeW;
    }

#pragma omp for

    for ( index = 0; index < count0; ++index ) {
      const long long index2 = index + index;
      const long long index4 = index2 + index2;
      const double latitude = latitudes [ index ];
      const double latitudeS = latitude - cellHeight;
      const double latitudeN = latitude + cellHeight;
      assert( IS_LATITUDE( latitude ) );
      assert( IS_LATITUDE( latitudeS ) );
      assert( IS_LATITUDE( latitudeN ) );
      vy[ index4     ] = latitudeS;
      vy[ index4 + 1 ] = latitudeS;
      vy[ index4 + 2 ] = latitudeN;
      vy[ index4 + 3 ] = latitudeN;
    }

    } /* End of parallel region. */

  } else { /* Project (longitude, latitude) to (x, y) and store ccw sequence:*/
    long long index = 0;

#pragma omp parallel for

    for ( index = 0; index < count0; ++index ) {
      const long long index2 = index + index;
      const long long index4_0 = index2 + index2;
      const long long index4_1 = index4_0 + 1;
      const long long index4_2 = index4_0 + 2;
      const long long index4_3 = index4_0 + 3;
      const double longitude = longitudes [ index ];
      const double longitudeW = longitude - cellWidth;
      const double longitudeE = longitude + cellWidth;
      const double latitude = latitudes [ index ];
      const double latitudeS = latitude - cellHeight;
      const double latitudeN = latitude + cellHeight;
      assert( IS_LONGITUDE( longitude ) );
      assert( IS_LONGITUDE( longitudeW ) );
      assert( IS_LONGITUDE( longitudeE ) );
      assert( IS_LATITUDE( latitude ) );
      assert( IS_LATITUDE( latitudeS ) );
      assert( IS_LATITUDE( latitudeN ) );
      project( projector, longitudeW, latitudeS, vx + index4_0, vy + index4_0 );
      project( projector, longitudeE, latitudeS, vx + index4_1, vy + index4_1 );
      project( projector, longitudeE, latitudeN, vx + index4_2, vy + index4_2 );
      project( projector, longitudeW, latitudeN, vx + index4_3, vy + index4_3 );
      assert( IN_RANGE( vx[ index4_0 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_0 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vx[ index4_1 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_1 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vx[ index4_2 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_2 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vx[ index4_3 ], -DBL_MAX, DBL_MAX ) );
      assert( IN_RANGE( vy[ index4_3 ], -DBL_MAX, DBL_MAX ) );
    }
  }
}



/******************************************************************************
PURPOSE: binQuadrilateralData - Bin 2D quadrilaterals onto a regular 2D grid
         and return the number of cells updated with data.
INPUTS:  const size_t count           Number of quadrilaterals.
         const double data[  count ]  Data values of quadrilaterals.
         const double x[ 4 * count ]  X-coordinates of counter-clockwise
                                      vertices of quadrilaterals.
                                      [ quad1_X1, quad1_X2, quad1_X3, quad1_X4,
                                        quad2_X1, quad2_X2, quad2_X3, quad2_X4,
                                        ...
                                        quadN_X1, quadN_X2, quadN_X3, quadN_X4]
         const double y[ 4 * count ]  Y-coordinates of counter-clockwise
                                      vertices of quadrilaterals.
         const size_t rows            Number of grid rows    of cells.
         const size_t columns         Number of grid columns of cells.
         const double gridXMinimum    X-coordinate of lower-left  corner of grid
         const double gridYMinimum    Y-coordinate of lower-left  corner of grid
         const double cellWidth       Width  of each grid cell.
         const double cellHeight      Height of each grid cell.
         size_t cellCounts[ rows * columns ]  Allocated and initialized array
                                              of count of values in each cell.
         double cellSums[ rows * columns ]    Allocated and initialized array
                                              of data values in each cell.
         double cellWeights[ rows * columns ] Allocated and initialized array
                                              of data weights in each cell or
                                              0 if weighting is not to be done.
OUTPUTS: size_t cellCounts[ rows * columns ]  Updated array of counts.
         double cellSums[  rows * columns ]   Updated array of sums.
RETURNS: size_t number of quadrilaterals that were binned.
******************************************************************************/

size_t binQuadrilateralData( const size_t count,
                             const double data[],
                             const double x[],
                             const double y[],
                             const size_t rows,
                             const size_t columns,
                             const double gridXMinimum,
                             const double gridYMinimum,
                             const double cellWidth,
                             const double cellHeight,
                             size_t cellCounts[],
                             double cellWeights[],
                             double cellSums[] ) {

  const long long count0 = count;
  long long index  = 0;
  long long result0 = 0; /* Must use signed type for OpenMP. */
  size_t result = 0;

  assert( count ); assert( data ); assert( x ); assert( y );
  assert( rows ); assert( columns );
  assert( cellWidth > 0.0 ); assert( cellHeight > 0.0 );
  assert( cellCounts );
  assert( cellSums );

  /*
   * Regrid each quadrilateral into grid cells.
   * Parallel execution uses critical section around updates to
   * cellCounts[], cellSums[], cellWeights[] in regridQuadrilateral().
   */

#pragma omp parallel for reduction ( + : result0 )

  for ( index = 0; index < count0; ++index ) {
    const long long index2 = index + index;
    const long long index4 = index2 + index2;
    result0 +=
      regridQuadrilateral( data[ index ], x + index4, y + index4,
                           rows, columns, gridXMinimum, gridYMinimum,
                           cellWidth, cellHeight,
                           cellCounts, cellWeights, cellSums );
  }

  result = result0;
  assert( result <= count );
  return result;
}



/******************************************************************************
PURPOSE: computeCellMeans - Compute mean of each cell's data.
INPUTS:  const double minimumValidValue       If a cell's mean is less than
                                              this then zero that cell count.
         const size_t count                   Number of cells.
         size_t cellCounts[ cellCount ]       Number of values in each cell.
         double cellWeights[ rows * columns ] Sum of weights in each cell or
                                              0 if weighting was not done.
         double cellSums[ rows * columns ]    Sum of values in each cell.
OUTPUTS: size_t cellCounts[ cellCount ]       Possibly reduced # of values/cell
         double cellWeights[  rows * columns ] Final weight value per cell.
         double cellSums[  rows * columns ]    Final mean value per cell.
RETURNS: size_t number of cells that contain data.
******************************************************************************/

size_t computeCellMeans( const double minimumValidValue,
                         const size_t count,
                         size_t cellCounts[],
                         double cellWeights[],
                         double cellSums[] ) {

  size_t result = 0;
  const long long count0 = (const long long) count;
  long long result0 = 0; /* Must use signed type for OpenMP. */
  long long cell = 0;

  assert( count >= 1 );
  assert( cellCounts );
  assert( cellSums );

#pragma omp parallel for reduction ( + : result0 )

  for ( cell = 0; cell < count0; ++cell ) {
    const size_t cellCount = cellCounts[ cell ];

    if ( cellCount ) {

      /* Convert sum to mean by dividing by weight or count: */

      if ( cellWeights ) {
        const double cellWeightSum = cellWeights[ cell ];
        assert( cellWeightSum > 0.0 );
        cellSums[ cell ] /= cellWeightSum;
      } else if ( cellCount > 1 ) {
        cellSums[ cell ] /= cellCount;
      }

      /* Filter-out cells whose final mean is below the mininumValidValue: */

      if ( cellSums[ cell ] >= minimumValidValue ) {
        ++result0;
      } else {
        cellCounts[ cell ] = 0;
      }
    }
  }

  result = result0;
  assert( result <= count );
  return result;
}



/******************************************************************************
PURPOSE: compactCells - Compute compact arrays of non-empty cells.
INPUTS:  void* projector                 Projector object or 0 if lonlat.
         UnprojectFunction unproject     Unproject function or 0 if lonlat.
         const size_t columns            Number of grid columns.
         const size_t rows               Number of grid rows.
         const double gridXMinimum       West  edge coordinate of grid.
         const double gridYMinimum       South edge coordinate of grid.
         const double cellWidth          Width  of each grid cell.
         const double cellHeight         Height of each grid cell.
         const size_t outputCells        Number of output grid cells with data.
         size_t cellCounts[ rows * columns ]  Count in each grid cell.
         double cellMeans[  rows * columns ]  Mean data in each grid cell.
OUTPUTS: size_t cellCounts[ outputCells ]  Compact count in each non-empty cell
         size_t cellMeans[ outputCells ]   Compact mean  in each non-empty cell
         double cellCenterLongitudes[  outputCells ]  Lon of output cells
         double cellCenterLatitudes[   outputCells ]  Lat of output cells
         size_t cellColumns[  outputCells ]           Col of output cells
         size_t cellRows[  outputCells ]              Row of output cells
******************************************************************************/

void compactCells( void* projector,
                   UnprojectFunction unproject,
                   const size_t columns,
                   const size_t rows,
                   const double gridXMinimum,
                   const double gridYMinimum,
                   const double cellWidth,
                   const double cellHeight,
                   const size_t outputCells,
                   size_t cellCounts[],
                   double cellMeans[],
                   double cellCenterLongitudes[],
                   double cellCenterLatitudes[],
                   size_t cellColumns[],
                   size_t cellRows[] ) {

  assert( ! projector || unproject ); /* projector implies unproject function*/
  assert(  columns ); assert( rows );
  assert( cellWidth > 0.0 ); assert( cellHeight > 0.0 );
  assert( outputCells >= 1 ); assert( outputCells <= rows * columns );
  assert( cellCounts ); assert( cellMeans );
  assert( cellCenterLongitudes ); assert( cellCenterLatitudes );
  assert( cellColumns ); assert( cellRows );

  size_t row = 0;
  size_t index = 0;
  size_t output = 0;
  double cellCenterY = gridYMinimum + 0.5 * cellHeight;

  for ( row = 0; row < rows; ++row, cellCenterY += cellHeight ) {
    double cellCenterX = gridXMinimum + 0.5 * cellWidth;
    size_t column = 0;

    for ( column = 0; column < columns; ++column, ++index,
          cellCenterX += cellWidth ) {
      assert( index < rows * columns );

      if ( cellCounts[ index ] ) {
        double cellCenterLongitude = cellCenterX;
        double cellCenterLatitude  = cellCenterY;
        const double mean = cellMeans[ index ];
        const size_t counts = cellCounts[ index ];

        if ( projector ) {
          unproject( projector, cellCenterX, cellCenterY,
                     &cellCenterLongitude, &cellCenterLatitude );
          DEBUG( fprintf( stderr, "(%lf %lf)->(%lf, %lf)\n",
                          cellCenterX, cellCenterY,
                          cellCenterLongitude, cellCenterLatitude ); )
        }

        assert( output < outputCells );
        cellCenterLongitudes[ output ] = cellCenterLongitude;
        cellCenterLatitudes[  output ] = cellCenterLatitude;
        cellRows[             output ] = row + 1;
        cellColumns[          output ] = column + 1;

        /* Compact cellCounts[], cellMeans[] to length outputCells: */

        assert( output <= index ); /* Not overwriting data to read. */
        cellMeans[  output ] = mean;
        cellCounts[ output ] = counts;
        ++output;
      }
    }
  }

  assert( output == outputCells );
  assert( IN_RANGE( cellCenterLongitudes[ 0 ], -180.0, 180.0 ) );
  assert( IN_RANGE( cellCenterLongitudes[ outputCells - 1 ], -180.0, 180.0 ) );
  assert( IN_RANGE( cellCenterLatitudes[ 0 ], -90.0, 90.0 ) );
  assert( IN_RANGE( cellCenterLatitudes[ outputCells - 1 ], -90.0, 90.0 ) );
  assert( cellCounts[ 0 ] );
  assert( cellCounts[ outputCells - 1 ] );
  assert( IN_RANGE( cellColumns[ 0 ], 1, rows * columns ) );
  assert( IN_RANGE( cellColumns[ outputCells - 1 ], 1, rows * columns ) );
  assert( IN_RANGE( cellRows[ 0 ], 1, rows * columns ) );
  assert( IN_RANGE( cellRows[ outputCells - 1 ], 1, rows * columns ) );
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: regridQuadrilateral - Regrid a 2D quadrilateral onto a regular 2D grid.
INPUTS:  const double data    Data value of quadrilateral.
         const double x[ 4 ]  X-coordinates of counter-clockwise
                              vertices of quadrilateral.
         const double y[ 4 ]  Y-coordinates of counter-clockwise
                                      vertices of quadrilateral.
         const size_t rows            Number of grid rows    of cells.
         const size_t columns         Number of grid columns of cells.
         const double gridXMinimum    X-coordinate of lower-left  corner of grid
         const double gridYMinimum    Y-coordinate of lower-left  corner of grid
         const double cellWidth       Width  of each grid cell.
         const double cellHeight      Height of each grid cell.
         size_t cellCounts[ rows * columns ]  Allocated array to hold number of
                                              values aggregated in each cell.
         double cellWeights[ rows * columns ] Allocated array to hold weights
                                              of values aggregated in each cell
                                              or 0 if weighting is not done.
         double cellSums[ rows * columns ]    Allocated array to hold means of
                                              values aggregated in each cell.
OUTPUTS: size_t cellCounts[ rows * columns ]  Initialized array:
                                              1 = cell has data,
                                              0 = cell has no data.
         double cellSums[  rows * columns ]   Initialized array:
                                              sum of values per cell.
RETURNS: int 1 if quadrilateral was regridded (i.e., overlaps the grid).
******************************************************************************/

static int regridQuadrilateral( const double data,
                                const double x[],
                                const double y[],
                                const size_t rows,
                                const size_t columns,
                                const double gridXMinimum,
                                const double gridYMinimum,
                                const double cellWidth,
                                const double cellHeight,
                                size_t cellCounts[],
                                double cellWeights[],
                                double cellSums[] ) {

  /*
   * If cellWeights is non-0 then
   *   For each grid cell intersected by the swath quadrilateral:
   *     Clip the swath quadrilateral to the rectangular grid cell yielding
   *     a polygon. clipped_quadrilateral_area = area of this clipped polygon
   *     and data_scale_factor = clipped_quadrilateral_area / quadrilateral_area
   */

  int result = 0;

  /* Compute grid bounds: */

  const double gridXMaximum = gridXMinimum + columns * cellWidth;
  const double gridYMaximum = gridYMinimum + rows    * cellHeight;

  Bounds gridBounds          = { { 0.0, 0.0 }, { 0.0, 0.0 } };
  Bounds quadrilateralBounds = { { 0.0, 0.0 }, { 0.0, 0.0 } };

  assert( x ); assert( y );
  assert( rows > 0 ); assert( columns > 0 );
  assert( cellWidth > 0.0 ); assert( cellHeight > 0.0 );
  assert( cellCounts );
  assert( cellSums );

  gridBounds[ X ][ MINIMUM ] = gridXMinimum;
  gridBounds[ X ][ MAXIMUM ] = gridXMaximum;
  gridBounds[ Y ][ MINIMUM ] = gridYMinimum;
  gridBounds[ Y ][ MAXIMUM ] = gridYMaximum;
  computePolygonBounds( 4, x, y, quadrilateralBounds );

  result =
    boundsOverlap( (const double (*)[2]) gridBounds,
                   (const double (*)[2]) quadrilateralBounds );

  DEBUG( fprintf( stderr,
                  "\n\nquad: %f [%f %f][%f %f] grid: [%f %f][%f %f] "
                  "overlap = %d\n",
                  data,
                  quadrilateralBounds[ X ][ MINIMUM ],
                  quadrilateralBounds[ X ][ MAXIMUM ],
                  quadrilateralBounds[ Y ][ MINIMUM ],
                  quadrilateralBounds[ Y ][ MAXIMUM ],
                  gridBounds[ X ][ MINIMUM ], gridBounds[ X ][ MAXIMUM ],
                  gridBounds[ Y ][ MINIMUM ], gridBounds[ Y ][ MAXIMUM ],
                  result ); )

  if ( result ) { /* Quadrilateral overlaps the grid: */
    const double oneOverCellWidth  = 1.0 / cellWidth;
    const double oneOverCellHeight =
      cellWidth == cellHeight ? oneOverCellWidth : 1.0 / cellHeight;
    const size_t rows_1    = rows    - 1;
    const size_t columns_1 = columns - 1;

    /* Get 0-based row/column ranges of grid cells that overlap the quad: */

    size_t firstRow    = 0;
    size_t lastRow     = 0;
    size_t firstColumn = 0;
    size_t lastColumn  = 0;
    const double deltaYMinimum =
      quadrilateralBounds[ Y ][ MINIMUM ] - gridYMinimum;
    const double deltaYMaximum =
      quadrilateralBounds[ Y ][ MAXIMUM ] - gridYMinimum;
    const double deltaXMinimum =
      quadrilateralBounds[ X ][ MINIMUM ] - gridXMinimum;
    const double deltaXMaximum =
      quadrilateralBounds[ X ][ MAXIMUM ] - gridXMinimum;
    int singleCell = 0;

    if ( deltaYMinimum > 0.0 ) {
      const double rowMinimum = deltaYMinimum * oneOverCellHeight + 1.0;
      DEBUG( fprintf( stderr, "rowMinimum = %f\n", rowMinimum ); )
      firstRow = rowMinimum;
      --firstRow;

      if ( firstRow > rows_1 ) {
        firstRow = rows_1;
      }
    }

    lastRow = firstRow;

    if ( deltaYMaximum > 0.0 ) {
      const double rowMaximum = deltaYMaximum * oneOverCellHeight + 1.0;
      DEBUG( fprintf( stderr, "rowMaximum = %f\n", rowMaximum ); )
      lastRow = rowMaximum;
      --lastRow;
      lastRow = CLAMPED_TO_RANGE( lastRow, firstRow, rows_1 );
    }

    if ( deltaXMinimum > 0.0 ) {
      const double columnMinimum = deltaXMinimum * oneOverCellWidth + 1.0;
      DEBUG( fprintf( stderr, "columnMinimum = %f\n", columnMinimum ); )
      firstColumn = columnMinimum;
      --firstColumn;

      if ( firstColumn > columns_1 ) {
        firstColumn = columns_1;
      }
    }

    lastColumn = firstColumn;

    if ( deltaXMaximum > 0.0 ) {
      const double columnMaximum = deltaXMaximum * oneOverCellWidth + 1.0;
      DEBUG( fprintf( stderr, "columnMaximum = %f\n", columnMaximum ); )
      lastColumn = columnMaximum;
      --lastColumn;
      lastColumn = CLAMPED_TO_RANGE( lastColumn, firstColumn, columns_1 );
    }

    assert( IN_RANGE( firstRow, 0, rows - 1 ) );
    assert( IN_RANGE( lastRow, firstRow, rows - 1 ) );
    assert( IN_RANGE( firstColumn, 0, columns - 1 ) );
    assert( IN_RANGE( lastColumn, firstColumn, columns - 1 ) );

    DEBUG( fprintf( stderr, "rows = [%lu %lu] columns = [%lu %lu]\n",
                    firstRow, lastRow, firstColumn, lastColumn ); )

    if ( cellWeights ) {

      /*
       * Check if quadrilateral is completely inside a single grid cell.
       * If true then no expensive clipping calculation is required.
       */

      singleCell = lastRow == firstRow && lastColumn == firstColumn;

      if ( singleCell ) {
        const int isInteriorCell =
          firstRow    > 0 && firstRow    < rows_1 &&
          firstColumn > 0 && firstColumn < columns_1;

        if ( ! isInteriorCell ) {

          /*
           * If the quad is not contained within an interior cell then
           * we must also check that quad does not extend outside the grid
           * (otherwise we must still perform clipping).
           * This is done by checking that the quad bounds are completely
           * within the single cell:
           */

          const double cellYMinimum = gridYMinimum + firstRow * cellHeight;
          const double cellYMaximum = cellYMinimum + cellHeight;
          int inside = IN_RANGE( quadrilateralBounds[ Y ][ MINIMUM ],
                                 cellYMinimum, cellYMaximum );

          if ( inside ) {
            inside = IN_RANGE( quadrilateralBounds[ Y ][ MAXIMUM ],
                               cellYMinimum, cellYMaximum );

            if ( inside ) {
              const double cellXMinimum =
                gridXMinimum + firstColumn * cellWidth;
              const double cellXMaximum = cellXMinimum + cellWidth;
              inside = IN_RANGE( quadrilateralBounds[ X ][ MINIMUM ],
                                 cellXMinimum, cellXMaximum );

              if ( inside ) {
                inside = IN_RANGE( quadrilateralBounds[ X ][ MAXIMUM ],
                                   cellXMinimum, cellXMaximum );
              }
            }
          }

          singleCell = inside;
        }
      }
    }

    if ( singleCell ) {

      /*
       * The quadrilateral is completely contained within a single grid cell
       * so just add the quadrlateral data value to the cell.
       */

      const size_t cellIndex = firstRow * columns + firstColumn;
      assert( cellIndex < rows * columns );
#pragma omp critical
      { /* Start of critical section: */

        if ( cellCounts[ cellIndex ] == 0 ) {
          cellCounts[ cellIndex ] = 1;
          cellSums[   cellIndex ] = data;

          if ( cellWeights ) {
            cellWeights[ cellIndex ] = 1.0;
          }

        } else {
          cellCounts[ cellIndex ] += 1;
          cellSums[   cellIndex ] += data;

          if ( cellWeights ) {
            cellWeights[ cellIndex ] += 1.0;
          }
        }

        DEBUG( fprintf( stderr, "single cell: data = %f ==>  "
                        "cell #%lu: count = %lu, sum = %f, weight = %f\n",
                        data, cellIndex,
                        cellCounts[ cellIndex ],
                        cellSums[ cellIndex ],
                        cellWeights ? cellWeights[ cellIndex ] : 0.0 ); )

      } /* End of critical section. */

    } else { /* For each overlapping cell, add the quad data. */

      if ( ! cellWeights ) {
        size_t row = 0;

        for ( row = firstRow; row <= lastRow; ++row ) {
          const size_t rowOffset = row * columns;
          size_t column = 0;

          for ( column = firstColumn; column <= lastColumn; ++column ) {
            const size_t cellIndex = rowOffset + column;
            assert( cellIndex < rows * columns );
#pragma omp critical
            { /* Start of critical section: */

              if ( cellCounts[ cellIndex ] == 0 ) {
                cellCounts[ cellIndex ] = 1;
                cellSums[   cellIndex ] = data;
              } else {
                cellCounts[ cellIndex ] += 1;
                cellSums[   cellIndex ] += data;
              }

              DEBUG( fprintf( stderr, "cell: data = %f ==> "
                              "cell #%lu: count = %lu, sum = %f\n",
                              data, cellIndex,
                              cellCounts[ cellIndex ], cellSums[ cellIndex ]);)

            } /* End of critical section. */
          }
        }

      } else { /* Use cellWeights[]: */

        /*
         * For each overlapping cell:
         *   add the quad data scaled by
         *   the fraction of the (assumed convex) quadrilateral in the cell
         *   fraction = clipped_quadrilateral_area / quadrilateral_area
         */

        double cellYMinimum = gridYMinimum + firstRow * cellHeight;
        size_t row = 0;

        for ( row = firstRow; row <= lastRow; ++row,
              cellYMinimum += cellHeight ) {
          const double cellYMaximum = cellYMinimum + cellHeight;
          double cellXMinimum = gridXMinimum + firstColumn * cellWidth;
          const size_t rowOffset = row * columns;
          size_t column = 0;

          for ( column = firstColumn; column <= lastColumn;
                ++column, cellXMinimum += cellWidth ) {
            const double cellXMaximum = cellXMinimum + cellWidth;
            const double quadrilateralArea =
              areaOfQuadrilateral( x[ 0 ], y[ 0 ],
                                   x[ 1 ], y[ 1 ],
                                   x[ 2 ], y[ 2 ],
                                   x[ 3 ], y[ 3 ] );

            DEBUG( fprintf( stderr, "quadrilateralArea = %f\n",
                            quadrilateralArea ); )
            assert( IN_RANGE( cellXMinimum, gridXMinimum, gridXMaximum ) );
            assert( IN_RANGE( cellXMaximum, cellXMinimum, gridXMaximum ) );
            assert( IN_RANGE( cellYMinimum, gridYMinimum, gridYMaximum ) );
            assert( IN_RANGE( cellYMaximum, cellYMinimum, gridYMaximum ) );

            if ( quadrilateralArea > 0.0 ) { /* Filter degenerates like TEMPO*/
              const double clippedPolygonArea =
                areaOfClippedQuadrilateral( cellXMinimum, cellYMinimum,
                                            cellXMaximum, cellYMaximum, x, y );
              const double fraction = clippedPolygonArea / quadrilateralArea;
              DEBUG( fprintf(stderr, "clippedPolygonArea = %f, fraction = %f\n",
                              clippedPolygonArea, fraction ); )
              assert( clippedPolygonArea >= 0.0 );

              if ( fraction > 0.0 ) {
                const double scaledData = fraction * data;
                const size_t cellIndex = rowOffset + column;
                assert( cellIndex < rows * columns );
#pragma omp critical
                { /* Start of critical section: */

                  if ( cellCounts[ cellIndex ] == 0 ) {
                    cellCounts[  cellIndex ] = 1;
                    cellSums[    cellIndex ] = scaledData;
                    cellWeights[ cellIndex ] = fraction;
                  } else {
                    cellCounts[  cellIndex ] += 1;
                    cellSums[    cellIndex ] += scaledData;
                    cellWeights[ cellIndex ] += fraction;
                  }

                  DEBUG( fprintf( stderr, "cell: data = %f ==>  "
                            "cell #%lu: count = %lu, sum = %f, weight = %f\n",
                                 data, cellIndex,
                                 cellCounts[  cellIndex ],
                                 cellSums[    cellIndex ],
                                 cellWeights[ cellIndex ] ); )

                } /* End of critical section. */
              }
            }
          }
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: areaOfClippedQuadrilateral - Clip quadrilateral to rectangular cell
         and return the area of the resulting polygon.
INPUTS:  const double cellXMinimum  X-coordinate of cell minimum.
         const double cellYMinimum  Y-coordinate of cell minimum.
         const double cellXMaximum  X-coordinate of cell maximum.
         const double cellYMaximum  Y-coordinate of cell maximum.
         const double x[ 4 ]        X-coordinates of counter-clockwise vertices
                                    of quadrilateral.
         const double y[ 4 ]        Y-coordinates of counter-clockwise vertices
                                    of quadrilateral.
RETURNS: double area of polygon (portion of quadrilateral inside cell).
******************************************************************************/

static double areaOfClippedQuadrilateral( const double cellXMinimum,
                                          const double cellYMinimum,
                                          const double cellXMaximum,
                                          const double cellYMaximum,
                                          const double x[],
                                          const double y[] ) {
  double clippedX[ 8 + 2 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0};
  double clippedY[ 8 + 2 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0};
  const int discardDegenerates = 0; /* Faster if 0. But use 1 if rendering. */
  const int clippedVertexCount =
    clipPolygon( discardDegenerates,
                 cellXMinimum, cellYMinimum, cellXMaximum, cellYMaximum,
                 4, x, y, clippedX, clippedY );
  double result = 0.0;
  assert( clippedVertexCount >= 0 );
  assert( clippedVertexCount < sizeof clippedX / sizeof *clippedX );

  if ( clippedVertexCount == 3 ) {
    result =
      areaOfTriangle( clippedX[ 0 ], clippedY[ 0 ],
                      clippedX[ 1 ], clippedY[ 1 ],
                      clippedX[ 2 ], clippedY[ 2 ] );
  } else if ( clippedVertexCount == 4 ) {
    result =
      areaOfQuadrilateral( clippedX[ 0 ], clippedY[ 0 ],
                           clippedX[ 1 ], clippedY[ 1 ],
                           clippedX[ 2 ], clippedY[ 2 ],
                           clippedX[ 3 ], clippedY[ 3 ] );
  } else {
    result = signedAreaOfPolygon( clippedVertexCount, clippedX, clippedY );

    if ( result < 0.0 ) {
      result = -result;
    }
  }

  DEBUG( fprintf( stderr,
                  "clippedVertexCount = %d, clippedPolygonArea = %f\n",
                  clippedVertexCount, result ); )
  assert( result >= 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: computePolygonBounds - Compute (axis-aligned)) bounds of a polygon.
         const size_t count   Number of polygon vertices.
         const double x[ count ]  X-coordinates of counter-clockwise
                                  vertices of polygon.
         const double y[ count ]  Y-coordinates of counter-clockwise
                                  vertices of polygon.
OUTPUTS: Bounds bounds            Bounds of polygon.
******************************************************************************/

static void computePolygonBounds( const size_t count,
                                  const double x[],
                                  const double y[],
                                  Bounds bounds ) {
  size_t index = 0;
  double xMinimum = x[ 0 ];
  double yMinimum = y[ 0 ];
  double xMaximum = xMinimum;
  double yMaximum = yMinimum;

  for ( index = 1; index < count; ++index ) {
    const double vx = x[ index ];
    const double vy = y[ index ];

    if ( vx < xMinimum ) {
      xMinimum = vx;
    } else if ( vx > xMaximum ) {
      xMaximum = vx;
    }

    if ( vy < yMinimum ) {
      yMinimum = vy;
    } else if ( vy > yMaximum ) {
      yMaximum = vy;
    }
  }

  bounds[ X ][ MINIMUM ] = xMinimum;
  bounds[ X ][ MAXIMUM ] = xMaximum;
  bounds[ Y ][ MINIMUM ] = yMinimum;
  bounds[ Y ][ MAXIMUM ] = yMaximum;
}



/******************************************************************************
PURPOSE: boundsOverlap - Do bounds overlap/intersect?
INPUTS:  Bounds a  1st bounds to check.
         Bounds b  2nd bounds to check.
RETURNS: int 1 if a overlaps/intersects b, else 0.
******************************************************************************/

static int boundsOverlap( const Bounds a, const Bounds b ) {
  const int outside =
    a[ Y ][ MINIMUM ] > b[ Y ][ MAXIMUM ] ||
    a[ Y ][ MAXIMUM ] < b[ Y ][ MINIMUM ] ||
    a[ X ][ MINIMUM ] > b[ X ][ MAXIMUM ] ||
    a[ X ][ MAXIMUM ] < b[ X ][ MINIMUM ];
  const int result = ! outside;
  return result;
}



#if 0
/* Unused unless weighting by grid cell is used. */

/******************************************************************************
PURPOSE: pointIsInsidePolygon - Determine if point (x, y) is inside polygon
         with vertices (vx[], vy[]).
INPUTS:  const double x                  X-Coordinate of point to test.
         const double y                  Y-Coordinate of point to test.
         const size_t vertexCount        Number of vertices in polygon.
         const double vx[ vertexCount ]  X-Coordinates of vertices of polygon.
         const double vy[ vertexCount ]  Y-Coordinates of vertices of polygon.
RETURNS: int 1 if inside, else 0.
NOTES:   Uses Winding Number algorithm:
         http://geomalgorithms.com/a03-_inclusion.html
******************************************************************************/

static int pointIsInsidePolygon( const double x, const double y,
                                 const size_t vertexCount,
                                 const double vx[], const double vy[] ) {
  int result = 1;
  long long windingNumber = 0; /* 0 = outside, else inside. */
  size_t vertex = 0;

  for ( vertex = 0; vertex < vertexCount; ++vertex ) {
    const size_t vertexp1 = vertex + 1;
    const size_t vertex1 = vertexp1 < vertexCount ? vertexp1 : 0;
    const double y1 = vy[ vertex  ];
    const double y2 = vy[ vertex1 ];

    if ( y1 <= y ) {

      if ( y2 > y ) {
        const double x1 = vx[ vertex  ];
        const double x2 = vx[ vertex1 ];
        const double cross = cross2d( x, y, x1, y1, x, y, x2, y2 );

        if ( cross > 0.0 ) {
          ++windingNumber;
        }
      }
    } else { /* y1 > y */

      if ( y2 <= y ) {
        const double x1 = vx[ vertex  ];
        const double x2 = vx[ vertex1 ];
        const double cross = cross2d( x, y, x1, y1, x, y, x2, y2 );

        if ( cross < 0.0 ) {
          --windingNumber;
        }
      }
    }
  }

  result = windingNumber != 0;
  return result;
}



/******************************************************************************
PURPOSE: cross2d - Computes 2D vector cross-product.
INPUTS:  const double from1X  X-coordinate of start of 1st vector.
         const double from1Y  Y-coordinate of start of 1st vector.
         const double to1X    X-coordinate of end   of 1st vector.
         const double to1Y    Y-coordinate of end   of 1st vector.
         const double from2X  X-coordinate of start of 2nd vector.
         const double from2Y  Y-coordinate of start of 2nd vector.
         const double to2X    X-coordinate of end   of 2nd vector.
         const double to2Y    Y-coordinate of end   of 2nd vector.
RETURNS: double signed 2D vector cross-product.
NOTES:   http://mathworld.wolfram.com/CrossProduct.html
******************************************************************************/

static double cross2d( const double from1X, const double from1Y,
                       const double to1X,   const double to1Y,
                       const double from2X, const double from2Y,
                       const double to2X,   const double to2Y ) {
  const double px = to1X - from1X;
  const double py = to1Y - from1Y;
  const double qx = to2X - from2X;
  const double qy = to2Y - from2Y;
  const double result = px * qy - qx * py;
  return result;
}

#endif



/******************************************************************************
PURPOSE: areaOfTriangle - Absolute/unsigned area of triangle with vertices
         (x1, y1), (x2, y2), (x3, y3).
INPUTS:  const double x1  X-Coordinate of 1st vertex of triangle.
         const double y1  Y-Coordinate of 1st vertex of triangle.
         const double x2  X-Coordinate of 2nd vertex of triangle.
         const double y2  Y-Coordinate of 2nd vertex of triangle.
         const double x3  X-Coordinate of 3rd vertex of triangle.
         const double y3  Y-Coordinate of 3rd vertex of triangle.
RETURNS: double area.
NOTES:   Area = 0.5 * |p X q|
         Where p and q are vectors of the triangle:
         p being from 1st vertex to 2nd vertex and
         q being from 1st vertex to 3rd vertex.
         and X is the 2D vector cross-product binary operator:
         p X q = px * qy - qx * py.
         http://geomalgorithms.com/a01-_area.html
******************************************************************************/

static double areaOfTriangle( const double x1, const double y1,
                              const double x2, const double y2,
                              const double x3, const double y3 ) {
  const double px = x2 - x1;
  const double py = y2 - y1;
  const double qx = x3 - x1;
  const double qy = y3 - y1;
  const double cross = px * qy - qx * py;
  const double abscross = cross < 0.0 ? -cross : cross;
  const double result = 0.5 * abscross;
  assert( result >= 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: areaOfQuadrilateral - Absolute/unsigned area of quadrilateral with
         vertices (x1, y1), (x2, y2), (x3, y3), (x4, y4).
INPUTS:  const double x1  X-Coordinate of 1st vertex of quadrilateral.
         const double y1  Y-Coordinate of 1st vertex of quadrilateral.
         const double x2  X-Coordinate of 2nd vertex of quadrilateral.
         const double y2  Y-Coordinate of 2nd vertex of quadrilateral.
         const double x3  X-Coordinate of 3rd vertex of quadrilateral.
         const double y3  Y-Coordinate of 3rd vertex of quadrilateral.
         const double x4  X-Coordinate of 4th vertex of quadrilateral.
         const double y4  Y-Coordinate of 4th vertex of quadrilateral.
RETURNS: double area.
NOTES:   Area = 0.5 * |p X q|
         Where p and q are diagonal vectors of the quadrilateral:
         p being from 1st vertex to 3rd vertex and
         q being from 4th vertex to 2nd vertex.
         and X is the 2D vector cross-product binary operator:
         p X q = px * qy - qx * py.
         http://mathworld.wolfram.com/Quadrilateral.html
         http://mathworld.wolfram.com/CrossProduct.html
******************************************************************************/

static double areaOfQuadrilateral( const double x1, const double y1,
                                   const double x2, const double y2,
                                   const double x3, const double y3,
                                   const double x4, const double y4 ) {
  const double px = x3 - x1;
  const double py = y3 - y1;
  const double qx = x2 - x4;
  const double qy = y2 - y4;
  const double cross = px * qy - qx * py;
  const double abscross = cross < 0.0 ? -cross : cross;
  const double result = 0.5 * abscross;
  assert( result >= 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: signedAreaOfPolygon - Signed area of a single contour of a polygon.
INPUTS:  const size_t count  Number of vertices in polygon.
         const double x[]    X-coordinates of vertices.
         const double y[]    Y-coordinates of vertices.
RETURNS: double signed area of polygon.
         Negative if vertices are in clockwise order.
NOTES:   http://mathworld.wolfram.com/PolygonArea.html
******************************************************************************/

static double signedAreaOfPolygon( const size_t count,
                                   const double x[], const double y[] ) {
  double result = 0.0;
  size_t index = 0;

  DEBUG( putc( '\n', stderr ); )

  for ( index = 0; index < count; ++index ) {
    const size_t indexp1 = index + 1;
    const size_t index1 = indexp1 < count ? indexp1 : 0;
    const double triangleArea =
      x[ index ] * y[ index1 ] - x[ index1 ] * y[ index ];
    result += triangleArea;
    DEBUG( fprintf( stderr, "  (%f, %f)", x[ index ], y[ index ] ); )
  }

  DEBUG( putc( '\n', stderr ); putc( '\n', stderr ); )

  result *= 0.5;
  return result;
}



/******************************************************************************
PURPOSE: clipPolygon - Clip polygon to an axis-aligned rectangle
         and return the number of vertices in clipped polygon.
INPUTS:  const int discardDegenerates  Discard degenerate triangles?
                                       Slower, but useful if rendering.
         const double clipXMin  X-coordinate of lower-left  corner of clip rect
         const double clipYMin  Y-coordinate of lower-left  corner of clip rect
         const double clipXMax  X-coordinate of upper-right corner of clip rect
         const double clipYMax  Y-coordinate of upper-right corner of clip rect
         const int count          Number of input vertices.
         const double x[ count ]  X-coordinates of input polygon to clip.
         const double y[ count ]  Y-coordinates of input polygon to clip.
         double cx[ 2 * count + 2 ] Storage for 2 * count + 2 X-coordinates.
         double cy[ 2 * count + 2 ] Storage for 2 * count + 2 Y-coordinates.
OUTPUTS: double cx[ result ]      X-coordinates of vertices of clipped poly.
         double cy[ result ]      Y-coordinates of vertices of clipped poly.
RETURNS: int number of vertices in clipped polygon.
NOTES:   Uses the Liang-Barsky polygon clipping algorithm. (Fastest known.)
         "An Analysis and Algorithm for Polygon Clipping",
         You-Dong Liang and Brian Barsky, UC Berkeley,
         CACM Vol 26 No. 11, November 1983.
         https://www.longsteve.com/fixmybugs/?page_id=210
******************************************************************************/

static int clipPolygon( const int discardDegenerates,
                        const double clipXMin,
                        const double clipYMin,
                        const double clipXMax,
                        const double clipYMax,
                        const int count,
                        const double x[],
                        const double y[],
                        double cx[],
                        double cy[] ) {
  int result = 0;
  const double inf = DBL_MAX;
  double xIn   = 0.0; /* X-coordinate of entry point. */
  double yIn   = 0.0; /* Y-coordinate of entry point. */
  double xOut  = 0.0; /* X-coordinate of exit point. */
  double yOut  = 0.0; /* Y-coordinate of exit point. */
  double tInX  = 0.0; /* Parameterized X-coordinate of entry intersection. */
  double tInY  = 0.0; /* Parameterized Y-coordinate of entry intersection. */
  double tOutX = 0.0; /* Parameterized X-coordinate of exit intersection. */
  double tOutY = 0.0; /* Parameterized Y-coordinate of exit intersection. */
  int vertex = 0;

  for ( vertex = 0; vertex < count; ++vertex ) {
    const int vertexp1 = vertex + 1;
    const int vertex1 = vertexp1 < count ? vertexp1 : 0;
    const double vx = x[ vertex ];
    const double vy = y[ vertex ];
    const double deltaX = x[ vertex1 ] - vx; /* Edge direction. */
    const double deltaY = y[ vertex1 ] - vy;
    const double oneOverDeltaX = deltaX ? 1.0 / deltaX : 0.0;
    const double oneOverDeltaY = deltaY ? 1.0 / deltaY : 0.0;
    double tOut1 = 0.0;
    double tOut2 = 0.0;
    double tIn2 = 0.0;
    assert( result < count + count + 2 );

    /*
     * Determine which bounding lines for the clip window the containing line
     * hits first:
     */

    if ( deltaX > 0.0 || ( deltaX == 0.0 && vx > clipXMax ) ) {
      xIn  = clipXMin;
      xOut = clipXMax;
    } else {
      xIn  = clipXMax;
      xOut = clipXMin;
    }

    if ( deltaY > 0.0 || ( deltaY == 0.0 && vy > clipYMax ) ) {
      yIn  = clipYMin;
      yOut = clipYMax;
    } else {
      yIn  = clipYMax;
      yOut = clipYMin;
    }

    /* Find the t values for the x and y exit points: */

    if ( deltaX != 0.0 ) {
      tOutX = ( xOut - vx ) * oneOverDeltaX;
    } else if ( vx <= clipXMax && clipXMin <= vx ) {
      tOutX = inf;
    } else {
      tOutX = -inf;
    }

    if ( deltaY != 0.0 ) {
      tOutY = ( yOut - vy ) * oneOverDeltaY;
    } else if ( vy <= clipYMax && clipYMin <= vy ) {
      tOutY = inf;
    } else {
      tOutY = -inf;
    }

    /* Set tOut1 = min( tOutX, tOutY ) and tOut2 = max( tOutX, tOutY ): */

    if ( tOutX < tOutY ) {
      tOut1 = tOutX;
      tOut2 = tOutY;
    } else {
      tOut1 = tOutY;
      tOut2 = tOutX;
    }

    if ( tOut2 > 0.0 ) {

      if ( deltaX != 0.0 ) {
        tInX = ( xIn - vx ) * oneOverDeltaX;
      } else {
        tInX = -inf;
      }

      if ( deltaY != 0.0 ) {
        tInY = ( yIn - vy ) * oneOverDeltaY;
      } else {
        tInY = -inf;
      }

      /* Set tIn2 = max( tInX, tInY ): */

      if ( tInX < tInY ) {
        tIn2 = tInY;
      } else {
        tIn2 = tInX;
      }

      if ( tOut1 < tIn2 ) { /* No visible segment. */

        if ( 0.0 < tOut1 && tOut1 <= 1.0 ) {
          assert( result < count + count );

          /* Line crosses over intermediate corner region. */

          if ( tInX < tInY ) {
            cx[ result ] = xOut;
            cy[ result ] = yIn;
          } else {
            cx[ result ] = xIn;
            cy[ result ] = yOut;
          }

          ++result;
        }
      } else { /* Line crosses through window: */

        if ( 0.0 < tOut1 && tIn2 <= 1.0 ) {

          if ( 0.0 <= tIn2 ) { /* Visible segment: */
            assert( result < count + count );

            if ( tInX > tInY ) {
              cx[ result ] = xIn;
              cy[ result ] = vy + ( tInX * deltaY );
            } else {
              cx[ result ] = vx + ( tInY * deltaX );
              cy[ result ] = yIn;
            }

            ++result;
          }

          assert( result < count + count );

          if ( 1.0 >= tOut1 ) {

            if ( tOutX < tOutY ) {
              cx[ result ] = xOut;
              cy[ result ] = vy + ( tOutX * deltaY );
            } else {
              cx[ result ] = vx + ( tOutY * deltaX );
              cy[ result ] = yOut;
            }

            ++result;
          } else {
            cx[ result ] = x[ vertex1 ];
            cy[ result ] = y[ vertex1 ];
            ++result;
          }
        }
      }

      if ( 0.0 < tOut2 && tOut2 <= 1.0 ) {
        assert( result < count + count );
        cx[ result ] = xOut;
        cy[ result ] = yOut;
        ++result;
      }
    }
  }

  /*
   * The above algorithm can generate 5-vertex 'line' or 'hat' polygons: _/\_
   * where the last 3 vertices are colinear
   * which yields a degenerate 'triangle' (i.e., with 0 area).
   * Here we discard the last 2 verticies in such cases.
   */

  if ( discardDegenerates && result == 5 ) {
    int twice = 2; /* Check twice in case of 5-vertex 'line'. */

    do {

      if ( result >= 3 ) {
        const size_t count_3 = result - 3;
        const size_t count_2 = result - 2;
        const size_t count_1 = result - 1;

        const double lastTriangleArea =
          areaOfTriangle( cx[ count_3 ], cy[ count_3 ],
                          cx[ count_2 ], cy[ count_2 ],
                          cx[ count_1 ], cy[ count_1 ] );

        DEBUG( fprintf( stderr, "lastTriangleArea = %f\n", lastTriangleArea );)

        if ( lastTriangleArea == 0.0 ) {
          result -= 2;
        }
      }
    } while (--twice );
  }

  /* Always discard any result less than a triangle. */

  if ( result < 3 ) {
    result = 0;
  }

  assert( result == 0 || IN_RANGE( result, 3, count + count + 2 ) );
  return result;
}



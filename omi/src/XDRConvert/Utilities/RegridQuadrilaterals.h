#ifndef REGRIDQUADRILATERALS_H
#define REGRIDQUADRILATERALS_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: RegridQuadrilaterals.h - Declares routine to regrid 2D quadrilaterals
         onto a regular 2D grid.

HISTORY: 2019-12-23 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

typedef void (*ProjectFunction)( void* object,
                                 double longitude, double latitude,
                                 double* x, double* y );

typedef void (*UnprojectFunction)( void* object,
                                   double x, double y,
                                   double* longitude, double* latitude );

/*================================ FUNCTIONS ================================*/

extern void
projectAndOrReorderQuadrilateralVertices( const size_t count,
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
                                          double vy[] );

extern void
projectAndOrReorderQuadrilateralVertices2( const size_t count,
                                           const double longitudes[],
                                           const double latitudes[],
                                           const double cellWidth,
                                           const double cellHeight,
                                           void* projector,
                                           ProjectFunction project,
                                           double vx[],
                                           double vy[] );

extern size_t binQuadrilateralData( const size_t count,
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
                                    double cellSums[] );

extern size_t computeCellMeans( const double minimumValidValue,
                                const size_t count,
                                size_t cellCounts[],
                                double cellWeights[],
                                double cellSums[] );

extern void compactCells( void* projector,
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
                          size_t cellRows[] );



#ifdef __cplusplus
}
#endif

#endif /* REGRIDQUADRILATERALS_H */




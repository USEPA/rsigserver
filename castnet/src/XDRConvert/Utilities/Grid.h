
#ifndef GRID_H
#define GRID_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Grid.h - Declare Grid ADT for projecting longitude-
         latitude-elevation points onto a 2D/3D regular grid in projected
         cartographic space, etc.

NOTES:   Commands:

           Grid* newGrid( Projector* projector,
                          Integer columns,   Integer rows,
                          Real    westEdge,  Real    southEdge,
                          Real    cellWidth, Real    cellHeight,
                          Integer layers, Integer type, Real topPressure,
                          const Real levels[],
                          Real g, Real R, Real A, Real T0s, Real P00 );

           Grid* newSubsetGrid( const Grid* grid,
                                const Integer firstLayer,
                                const Integer lastLayer,
                                const Integer firstRow,
                                const Integer lastRow,
                                const Integer firstColumn,
                                const Integer lastColumn );

           Grid* parseGrid( int argc, char* argv[], Integer* argument,
                            Projector* projector );

           Projector* parseProjection( int argc, char* argv[],
                                       Integer* argument );

           void free( Grid* self );

           Integer projectXY( Grid* self, Integer count,
                              const Real longitudes[], const Real latitudes[],
                              Integer columns[], Integer rows[],
                              Real xCenterOffsets[], Real yCenterOffsets[],
                              Real gridLongitudes[], Real gridLatitudes[] );


           Integer projectZ( Grid* self, Integer count,
                             const Real elevations[],
                             Integer layers[], Real centerOffsets[],
                             Real gridElevations[] );

           void aggregate( Grid* grid,
                           Integer method, Real minimumValidValue,
                           Integer inputPoints,
                           Integer columns[], Integer rows[],
                           const Real xOffsets[], const Real yOffsets[],
                           Real gridLongitudes[], Real gridLatitudes[],
                           const Integer layers, const Real elevations[],
                           const Real inputData[],
                           Integer* outputPoints, Real outputData[] );

           void regrid( Grid* self,
                        Integer method, Real minimumValidValue,
                        Integer points, Integer levels,
                        const Real longitudes[], const Real latitudes[],
                        const Real elevations[], const Real inputData[],
                        Integer* outputPoints,
                        Integer columns[], Integer rows[], Integer layers[],
                        Real gridLongitudes[], Real gridLatitudes[],
                        Real gridElevations[], Real outputData[] );

           void regridSwath( Grid* self,
                             Integer method,
                             Integer points,
                             const Real longitudesSW[],
                             const Real longitudesSE[],
                             const Real longitudesNW[],
                             const Real longitudesNE[],
                             const Real latitudesSW[],
                             const Real latitudesSE[],
                             const Real latitudesNW[],
                             const Real latitudesNE[],
                             const Real inputData[],
                             Integer* outputPoints,
                             Integer columns[],
                             Integer rows[],
                             Real gridLongitudes[],
                             Real gridLatitudes[],
                             Real outputData[] );

         Queries:

           Integer invariant( const Grid* self );

           Integer equal( const Grid* self,
                          const Grid* other );

           Grid* clone( const Grid* self );

           Projector* projector( const Grid* self );

           Integer layers(  const Grid* self );

           Integer rows(    const Grid* self );

           Integer columns( const Grid* self );

           Real longitude(  const Grid* self, Integer row, Integer column );

           Real latitude(   const Grid* self, Integer row, Integer column );

           Real elevation(  const Grid* self, Integer layer );

           Real level(      const Grid* self, Integer level );

           Real westEdge(   const Grid* self );

           Real southEdge(  const Grid* self );

           Real cellWidth(  const Grid* self );

           Real cellHeight( const Grid* self );

         Example:

#include <stdio.h>
#include <Lambert.h>
#include <Grid.h>

int main( void ) {
  const Integer method = AGGREGATE_WEIGHTED;
  const Real minimumValidValue = -9000.0;
  Projector* projector = (Projector*)
    newLambert( 6370997.0, 6370997.0, 33.0, 45.0, -97.0, 40.0, 0.0, 0.0 );

  if ( projector ) {
    enum { LAYERS = 22 };
    const Real levels[ LAYERS + 1 ] = {
      1.0, 0.995, 0.988, 0.979, 0.97, 0.96, 0.938, 0.914, 0.889, 0.862,
      0.834, 0.804, 0.774, 0.743, 0.694, 0.644, 0.592, 0.502, 0.408, 0.311,
      0.21, 0.106, 0.0
    };
    Grid* grid =
      newGrid( projector, 268, 259, -420000.0, -1716000.0,
               12000.0, 12000.0, LAYERS, VGSGPN3, 10000.0, levels,
               9.81, 287.04, 50.0, 290.0, 100000.0 );

    if ( grid ) {
      enum { POINTS = 7, ELEVATIONS = 8 };
      const Real longitudes[ POINTS] =
        { -75.0, -77.0, 78.0, -77.01, -74.99, -74.98, -70.0 };
      const Real latitudes[ POINTS ] =
        { 35.0, 36.0, 37.0, 36.01, 34.0, 33.9, 30.0 };
      const Real elevations[ ELEVATIONS ] =
        { 100.0, 100.0, 120.0, 121.0, 1000.0, 1001.0, 5000.0, 20000.0 };
      const Real data[ POINTS * ELEVATIONS ] = {
        11.1, 11.5, 22.2, 22.3, -9999.0, 44.4, -9999.0, 1.0,
        111.1, 111.5, 222.2, 222.3, -9999.0, 444.4, -9999.0, 1.0,
        1111.1, 1111.5, 2222.2, 2222.3, -9999.0, 4444.4, -9999.0, 1.0,
        11111.1, 11111.5, 22222.2, 22222.3, -9999.0, 44444.4, -9999.0, 1.0,
        111111.1, 111111.5, 222222.2, 222222.3, -9999.0, 444444.4, -9999.0,1.0,
        1111111.1,1111111.5,2222222.2,2222222.3,-9999.0,4444444.4,-9999.0,1.0,
        11111111.1,11111111.5,22222222.2,22222222.3,-9999.0,44444444.4,-9999.0,
        1.0,
      };
      Integer columns[ POINTS ]  = { 0, 0, 0, 0, 0, 0, 0 };
      Integer rows[ POINTS ]     = { 0, 0, 0, 0, 0, 0, 0 };
      Real    xOffsets[ POINTS ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      Real    yOffsets[ POINTS ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      Real    gridLongitudes[ POINTS ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      Real    gridLatitudes[ POINTS ]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      Real    gridData[ POINTS ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
      Integer pointsXY = 0;
      projector = 0;

      printf( "calling projectXY:\n" );
      grid->projectXY( grid, POINTS, longitudes, latitudes,
                       &pointsXY, columns, rows, xOffsets, yOffsets,
                       gridLongitudes, gridLatitudes );
      printf( "pointsXY = %"INTEGER_FORMAT"\n", pointsXY );

      if ( pointsXY ) {
        Integer point = 0;

        do {
          printf( "(%lf, %lf) -> [%"
                  INTEGER_FORMAT", %"INTEGER_FORMAT"] + [%lf, %lf] @ (%lf, %lf)\n",
                  longitudes[ point ], latitudes[ point ],
                  columns[ point ], rows[ point ],
                  xOffsets[ point ], yOffsets[ point ],
                  gridLongitudes[ point ], gridLatitudes[ point ] );
          ++point;
        } while ( point < POINTS );

        printf( "calling aggregate:\n" );

        grid->aggregate( grid, method, minimumValidValue, POINTS,
                         columns, rows, xOffsets, yOffsets,
                         gridLongitudes, gridLatitudes,
                         ELEVATIONS, elevations, data, 0,
                         &pointsXY, gridData, 0, 0 );

        printf( "pointsXY = %"INTEGER_FORMAT"\n", pointsXY );

        if ( pointsXY ) {
          point = 0;

          do {
            printf( "[%"INTEGER_FORMAT", %"INTEGER_FORMAT
                    "] @ (%lf, %lf): %lf\n",
                    columns[ point ], rows[ point ],
                    gridLongitudes[ point ], gridLatitudes[ point ],
                    gridData[ point ] );
            ++point;
          } while ( point < pointsXY );
        }
      }

      FREE_OBJECT( grid );
    }

    FREE_OBJECT( projector );
  }

  return 0;
}


HISTORY: 2007/12, plessel.todd@epa.gov, Created.

STATUS:  unreviewed, untested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <BasicNumerics.h> /* For Real, Integer. */
#include <Memory.h>        /* For macro FREE_OBJECT() used by clients. */
#include <Stream.h>        /* For Stream. */
#include <Projector.h>     /* For Projector. */

/*================================== TYPES ==================================*/

/* Note and RegriddedNote are used to hold Aircraft flight info: */

/*
 * For non-regrided files, there is only one flight per track:
 * cruising, decending, ascending, for example:
 * M20060703014:FRANKFURT->ATLANTA
 * MD20060703014:FRANKFURT->ATLANTA
 * MA20060703054:ATLANTA->FRANKFURT
 */

enum { NOTE_LENGTH = 79 };
typedef char Note[ NOTE_LENGTH + 1 ];

/*
 * For regrided files, there are no tracks, rather points from multiple tracks
 * (e.g., cruising, decending, ascending - even independet flights) are merged
 * into a 1-hour cell value. For large CMAQ cells, (e.g., 108km) there could
 * be multiple airports and multiple sets of tracks per cell.
 * How many is unknowable a priori.
 * Instead a reasonable length limit is used while concatenating above Notes.
 * Caution: increasing this estimate significantly increases the size
 * (memory, transmission time) of the regridded data since it is per
 * regridded cell value.
 */

enum { REGRIDDED_NOTE_LENGTH = 255 };
typedef char RegriddedNote[ REGRIDDED_NOTE_LENGTH + 1 ];


typedef struct Grid Grid;

/* Constructor: */

extern Grid* newGrid( Projector* projector,
                      Integer columns, Integer rows,
                      Real westEdge,   Real southEdge,
                      Real cellWidth,  Real cellHeight,
                      Integer layers,
                      Integer type, Real topPressure, const Real levels[],
                      Real g, Real R, Real A, Real T0s, Real P00 );

extern Grid* newSubsetGrid( const Grid* grid,
                            const Integer firstLayer, const Integer lastLayer,
                            const Integer firstRow, const Integer lastRow,
                          const Integer firstColumn, const Integer lastColumn);

extern Grid* parseGrid( int argc, char* argv[], Integer* argument,
                        Projector* projector );

extern Projector* parseProjection( int argc, char* argv[], Integer* argument );

extern void parseEllipsoid( int argc, char* argv[],
                            Integer* argument,
                            Real* majorSemiaxis, Real* minorSemiaxis,
                            Integer* ok );

extern void writeProjectionAndGrid( const Grid* grid, Stream* stream );

/* Destructor: use macro FREE_OBJECT( object ); */

/* Missing data: */

#define AMISS3 (-9.000E36)
#define BADVAL3 (-9.999E36)

/* Vertical grid types: */

enum {
  IMISS3  = -9999, /* None (becomes a single layer at 0.0). */
  VGSGPH3 = 1,     /* Hydrostatic sigma-P */
  VGSGPN3,         /* Non-h sigma-P */
  VGSIGZ3,         /* Sigma-Z */
  VGPRES3,         /* Pressure (pascals) */
  VGZVAL3,         /* Z (m) (above sea level) */
  VGHVAL3,         /* Z (m) (above terrain) */
  VGWRFEM          /* Sigma-P WRF. */
};

enum { MXLAYS3 = 100 };

#define IS_VALID_VERTICAL_GRID_TYPE( type ) \
          OR2( IN_RANGE( (type), VGSGPH3, VGWRFEM ), (type) == IMISS3 )

enum { AGGREGATE_NEAREST = 1, AGGREGATE_MEAN, AGGREGATE_WEIGHTED };

#define IS_VALID_AGGREGATE_METHOD( method ) \
  IN_RANGE( (method), AGGREGATE_NEAREST, AGGREGATE_WEIGHTED )

/*= PRIVATE =*/

typedef struct GridPrivate GridPrivate;

struct Grid {

  /* Commands: */

  void (*free)( Grid* self ); /* Called by FREE_OBJECT(). */

  void (*projectXY)( Grid* self, Integer count,
                     const Real longitudes[], const Real latitudes[],
                     Integer* griddedPoints,
                     Integer columns[], Integer rows[],
                     Real xCenterOffsets[], Real yCenterOffsets[],
                     Real gridLongitudes[], Real gridLatitudes[] );

  void (*projectZ)( Grid* self, Integer count,
                    const Real elevations[],
                    Integer* griddedPoints,
                    Integer layers[], Real centerOffsets[],
                    Real gridElevations[] );

  void (*aggregate)( Grid* grid,
                     Integer method, Real minimumValidValue,
                     Integer inputPoints,
                     Integer columns[], Integer rows[],
                     const Real xOffsets[], const Real yOffsets[],
                     Real gridLongitudes[], Real gridLatitudes[],
                     const Integer layers, const Real elevations[],
                     const Real inputData[], const Real inputData2[],
                     Integer* outputPoints,
                     Real outputData[], Real outputData2[],
                     Real gridElevations[] );

  void (*regrid)( Grid* self,
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
                  RegriddedNote regriddedNotes[] );

  void (*regridSwath)( Grid* self,
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

  Integer (*invariant)( const Grid* self );
  Integer (*equal)( const Grid* self, const Grid* other );
  Grid* (*clone)( const Grid* self );
  Projector* (*projector)( const Grid* self );
  Integer (*layers)(  const Grid* self );
  Integer (*rows)(    const Grid* self );
  Integer (*columns)( const Grid* self );
  Real (*longitude)(  const Grid* self, Integer row, Integer column );
  Real (*latitude)(   const Grid* self, Integer row, Integer column );
  Real (*elevation)(  const Grid* self, Integer layer );
  Real (*level)(      const Grid* self, Integer level );
  Real (*westEdge)(   const Grid* self );
  Real (*southEdge)(  const Grid* self );
  Real (*cellWidth)(  const Grid* self );
  Real (*cellHeight)( const Grid* self );

  /* Private: */

  GridPrivate* data;
};


#ifdef __cplusplus
}
#endif

#endif /* GRID_H */




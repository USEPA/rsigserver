
#ifndef SHAPEFILE_H
#define SHAPEFILE_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Shapefile.h - Declare routines for ESRI Shapefile creation.

NOTES:   Uses Shapefile and GPC free, open-source libraries.
         http://shapelib.maptools.org/shp_api.html
         http://www.cs.man.ac.uk/~toby/alan/software/gpc.html

HISTORY: 2010/11-23 Todd Plessel EPA/LM Created.

STATUS:  unreviewed, untested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <shapefil.h> /* For SHP*, DBF*. */
#include <gpc.h> /* For gpc_**. */

#include <Utilities.h> /* For Bounds. */

/*================================== TYPES ==================================*/

typedef struct {
  int id;                 /* From original input Shapefile. */
  Bounds bounds;          /* Bounds of polygon or triangles. */
  gpc_polygon polygon;    /* After clipping to bounds. */
  gpc_tristrip triangles; /* Triangle strip for rendering. */
} PolygonShape;

/* Data read from a DBF file: */

typedef union {
  const char* s; /* String value. */
  int         i; /* int value. */
  double      d; /* double value. */
} Value;

typedef struct {
  int    rows;          /* Number of data rows. */
  int    columns;       /* Number of data columns. */
  char** columnNames;   /* columnNames[ columns ] name of each column. */
  int*   columnTypes;   /* columnTypes[columns] FTString, FTInteger, FTDouble*/
  Value* values;        /* values[ rows * columns ]. */
  char** stringStorage; /* Array of allocated unique strings Flyweight.*/
  int    capacity;      /* Number of strings in stringStorage[]. */
} ShapeData;

typedef Color (*TextColor)( const char* );

/*================================ FUNCTIONS ================================*/

/* Point, Polyline and Grid routines: */

extern int writeASCIIGridFile( const char* fileName,
                               int rows, int columns,
                               double westEdge, double southEdge,
                               double cellSize,
                               int type,
                               const void* data );

extern int writeLambertPRJFile( const char* fileName,
                                const double centralLongitude,
                                const double centralLatitude,
                                const double lowerLatitude,
                                const double upperLatitude,
                                const int useASCIIGridForm );

extern int writeWGS84PRJFile( const char* fileName,
                              const int useASCIIGridForm );

extern int writeGridToShapefile( const char* baseFileName,
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
                                 Unproject unproject );

extern void writePointsToShapefile( const char* const fileName,
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
                                    const int writeCSV );

extern int writePolylinesToShapefile( const char* baseFileName,
                                      const int polyLineCount,
                                      const int vertexCount,
                                      const int counts[],
                                      const float lonlats[] );

extern void computePolygonVertices( const int row, const int column,
                                    const double westEdge,
                                    const double southEdge,
                                    const double cellWidth,
                                    const double cellHeight,
                                    Unproject unproject,
                                    double xy[ 5 * 2 ],
                                    double xyRange[ 2 * 2 ] );

extern void computeGridCellVertices( const int rows,
                                     const int columns,
                                     const float longitudes[],
                                     const float latitudes[],
                                     float vertices[] );

/* SHP, GPC PolygonShape routines: */

extern void printShape( const SHPObject* shape );

extern void printPolygon( const gpc_polygon* polygon );
  
extern void printTriangles( const gpc_tristrip* tristrip );

extern void deallocatePolygons( int count, PolygonShape polygons[] );

extern int shapefileType( const char* baseFileName );

extern int computeShapeSubsetBoundsMask( const char* baseFileName,
                                         const Bounds bounds,
                                         const int count,
                                         char mask[] );

extern int writeSubsetCSVFile( const char* inputDBFFileName,
                               const char* inputCSVDirectory,
                               const char* outputCSVFileName,
                               const int yyyymmdd1,
                               const int yyyymmdd2,
                               const char* estcode,
                               const char mask[] );

extern int writeSubsetCSVFileById( const char* const inputDBFFileName,
                                   const char* const inputCSVFileName,
                                   const char* const outputCSVFileName,
                                   const char* const columnName,
                                   const int allowEmptyOutputCSV,
                                   const char mask[] );

extern PolygonShape* readAndClipShapes( const char* baseFileName,
                                        const Bounds bounds,
                                        double minimumAdjacentVertexDistance,
                                        char* const mask,
                                        int* count, int* isPolyline );

extern PolygonShape* readAndTriangulateShapes( const char* baseFileName,
                                           const double minimumAdjacentVertexDistance,
                                               int* count );

extern int pointInTriangles( double x, double y,
                             int count, const PolygonShape polygons[] );

extern int makePolygon( const SHPObject* shape,
                        double minimumAdjacentVertexDistance,
                        gpc_polygon* polygon,
                        Bounds bounds );

extern int clipPolylines( const gpc_polygon* const polylines,
                          const Bounds clipBounds,
                          gpc_polygon* clippedPolylines,
                          Bounds polylinesBounds );

extern int nearestPolyline( const double x, const double y,
                            const int count, const PolygonShape polylines[] );

extern int nearestPoint( const double x, const double y,
                         const ShapeData* const shapeData );

extern int ensureCorrectVertexOrder( gpc_polygon* polygon );

extern void reverseVertexList( gpc_vertex_list* vertex_list );

extern double signedAreaOfPolygon( const gpc_vertex_list* vertex_list );

extern double polygonArea( const gpc_polygon* polygon );

extern double polygonPerimeter( const gpc_polygon* polygon );

extern int maximumPolygonContours( int count, const PolygonShape polygons[] );

extern int maximumPolygonertices( int count, const PolygonShape polygons[] );

extern int polygonVertexCount( const gpc_polygon* polygon );

extern void copyPolygonVertices( const gpc_polygon* polygon, int closeRing,
                                 int starts[], double x[], double y[] );

extern void copyPolylineVertices( const gpc_polygon* polygon,
                                  int counts[], double xy[] );

extern int createPolyline( const int polylineCount, const int vertexCounts[],
                           const double xy[], gpc_polygon* polygon,
                           Bounds bounds );

extern int writePolygonsToShapefile( SHPHandle outputFile, int isPolyline,
                                     int count, const PolygonShape polygons[]);

extern int copySubsetShapefile( SHPHandle inputFile, SHPHandle outputFile,
                                const int count, const char* mask );

extern int writePolygonDBF( const char* inputFileName, DBFHandle outputFile,
                            int offset, int count, const char* mask,
                            const PolygonShape polygons[] );

extern int writePointSubsetDBF( const char* inputFileName,
                                const Bounds,
                                const long long huc,
                                const char* const estcode,
                                char** mask,
                                DBFHandle outputFile );

extern int writeSubsetDBF( const char* inputFileName,
                           const Bounds bounds,
                           const long long huc,
                           const char* const estcode,
                           const int siteIdCount,
                           const int siteIds[],
                           const int count,
                           char mask[],
                           DBFHandle outputFile );

extern int getRowsDBF( const char* baseFileName );

extern int writeBoundsDBF( const char* fileName,
                           const double areaInSquareKilometers );

extern int streamShapefiles( const char* baseFileName, const char* name,
                             const int dbfOnly, const int csv );

extern void removeShapefiles( const char* baseFileName );

/* ShapeData routines: */

extern ShapeData* readDBF( const char* fileName );

extern int subsetDBFByTime( const char* fileName,
                            const int yyyymmdd1,
                            const int yyyymmdd2,
                            const int count,
                            char mask[] );

extern int subsetDBFByTimeAndBoundsOrEstcode( const char* fileName,
                                              const int yyyymmdd1,
                                              const int yyyymmdd2,
                                              const Bounds bounds,
                                              const char* estcode,
                                              const int count,
                                              char mask[] );

extern int subsetByCOMID( const ShapeData* const shapeData,
                          const int comid,
                          char mask[] );

extern int writeShapeDataToText( const ShapeData* shapeData,
                                 const char* outputFileName );

extern int writeShapeData( const char* baseFileName,
                           const ShapeData* shapeData );

extern void deallocateShapeData( ShapeData* shapeData );

extern void printShapeData( const ShapeData* shapeData );

extern int isValidShapeData( const ShapeData* shapeData );

extern int copyMaskedShapes( const char* baseInputFileName,
                             const char* baseOutputFileName,
                             SHPHandle* shpFile,
                             const int count, const char mask[] );
  
extern int isPointType( const char* baseFileName );

/* Value routines: */

extern int isValidColumnName( const char* columnName );

extern int isValidValue( int type, const char* units, Value value );

extern int compareValues( int type, Value value1, Value value2 );

#ifdef __cplusplus
}
#endif

#endif /* SHAPEFILE_H */



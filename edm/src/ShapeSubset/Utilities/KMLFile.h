
#ifndef KMLFILE_H
#define KMLFILE_H


#ifdef __cplusplus
extern "C" {
#endif
  
/******************************************************************************
PURPOSE: KMLFile.h - Declare routines for KML file creation.

NOTES:   

HISTORY: 2010/12/22 Todd Plessel EPA/LM Created.

STATUS:  unreviewed, untested.
******************************************************************************/

/*================================= INCLUDES ================================*/

#include <stdio.h> /* For FILE. */

#include <Utilities.h> /* For Bounds, Color, BGR, RGBColormap, PointData. */
#include <Shapefile.h> /* For PolygonShape, ShapeData. */

/*================================ FUNCTIONS ================================*/

extern void writeStartKML( FILE* file,
                           const char* name, const char* description,
                           const char* boundsName,
                           const char* boundsDescription,
                           const Bounds bounds );
  
extern void writeEndKML( FILE* file );

extern void writeMapSubsetToKML( FILE* file,
                                 int polylineCount,
                                 const int counts[],
                                 const float vertices[],
                                 const char* name,
                                 Color color );

extern void writePointsToKML( FILE* file,
                              const char* const source,
                              const char* const variableName,
                              const char* const units,
                              const int timesteps,
                              const int hoursPerTimestep,
                              const int timestamps[],
                              const int stations,
                              const long long id,
                              const char* sids[],
                              const float lonlats[],
                              const float z[],
                              const float zAll,
                              const int components,
                              const float dataMinimum,
                              const float dataMaximum,
                              const float data[] );

extern void writeShapeDataToKML( FILE* file, const char* name,
                                 const char* units,
                                 const char* description,
                                 int column,
                                 double minimumValue, double maximumValue,
                                 DataColor dataColor, TextColor textColor,
                                 const ShapeData* shapeData );

extern void writeGridToKML( FILE* file,
                            const int timestamps[],
                            int yyyymmdd,
                            int components,
                            int timesteps,
                            int rows,
                            int columns,
                            double westEdge,
                            double southEdge,
                            double cellWidth,
                            double cellHeight,
                            const LongName source,
                            const Name name,
                            const Units units, 
                            const float corners[],
                            int type,
                            const void* data,
                            const float dataRange[ 2 ],
                            const RGBColormap colormap,
                            DataColor dataColor );

extern void writePolygonsToKML( FILE* file,
                                const int yyyymmddhhStart,
                                const int yyyymmddhhEnd,
                                const int count,
                                const PolygonShape polygons[],
                                const ShapeData* shapeData,
                                const char* stringValues[],
                                const double csvValues[],
                                const LongName source,
                                const Name name,
                                const Units units,
                                const double dataRange[ 2 ],
                                const RGBColormap* colormap,
                                DataColor dataColor,
                                TextColor textColor );

extern void writeMapImagesToKML( FILE* file,
                                 const char* const inputDirectory,
                                 const char* const outputDirectory,
                                 const Bounds bounds );
  
#ifdef __cplusplus
}
#endif
    
#endif /* KMLFILE_H */



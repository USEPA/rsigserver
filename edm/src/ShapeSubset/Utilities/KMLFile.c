
/******************************************************************************
PURPOSE: KMLFile.c - Write KML files.
NOTES:
HISTORY: 2011-06-17 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For FILE, stderr, fprintf(), fopen(), perror(). */
#include <string.h> /* For strncpy(), memset(). */
#include <limits.h> /* For INT_MAX. */
#include <math.h>   /* For sqrt(). */

#include <Assertions.h>  /* For PRE(), IMPLIES_ELSE(). */
#include <BasicNumerics.h> /* For isNan().*/
#include <Utilities.h>  /* For Bounds, Color, RGBColormap, forEachFile(). */
#include <Shapefile.h>  /* For isValidValue(). */
#include <DateTime.h>   /* For nowUTC(), monthAndDay(). */
#include <KMLFile.h>    /* For public interface. */

/*================================== TYPES ==================================*/

enum { KML_COLOR_LEVELS = 8 };
assert_static( 256 % KML_COLOR_LEVELS == 0 );

typedef char BGR[ 7 ]; /* KML blue-green-red hexedecimal color string. */

typedef struct {
  FILE* file;
  const char* directory;
  Bounds bounds;
} MapImageFileInfo;

/*=========================== FORWARD DECLARATIONS ==========================*/

static int isValidMapImageFileInfo( const MapImageFileInfo* const info );

static void writeMapImagesToKMLHelper( const char* const fileName, void* data );

static void writeKMLColors( FILE* file, int levels );

static void writeKMLColormap( FILE* file, const RGBColormap colormap,
                              unsigned char opacity );

static void colorKML( Color color, int levels, BGR bgr );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: writeStartKML - Write start of KML file.
INPUTS:  FILE* file                     File to write to.
         const char* name               Name of document.
         const char* description        Description of document.
         const char* boundsName         Name of bounds.
         const char* boundsDescription  Description of bounds.
         const Boubnds bounds           Bounds of KML data.
******************************************************************************/

void writeStartKML( FILE* file,
                    const char* name, const char* description,
                    const char* boundsName, const char* boundsDescription,
                    const Bounds bounds ) {

  PRE010( file, name, *name, description, *description,
          boundsName, *boundsName, boundsDescription, *boundsDescription,
          isValidBounds( bounds ) );

  const char* const fileStartingContent =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<kml xmlns=\"http://earth.google.com/kml/2.1\">\n"
    "  <Document>\n"
    "    <name>EDM</name>\n"
    "    <description>Estuary Data Mappper static estuarine boundary "
    "lines "
    "and hourly measured site and/or modeled gridded data."
    "</description>\n";

  fprintf( file, "%s", fileStartingContent );
  writeKMLColors( file, KML_COLOR_LEVELS );

  /* Save Overall Bounds: */

  fprintf( file,
           "    <Placemark>\n"
           "      <name>Overall_Bounds</name>\n"
           "      <description>Boundary of estuarine region of study."
           "</description>\n"
           "      <styleUrl>#ffffffff</styleUrl>\n"
           "      <LineString>\n"
           "        <coordinates>\n"
           "          %lf,%lf\n"
           "          %lf,%lf\n"
           "          %lf,%lf\n"
           "          %lf,%lf\n"
           "          %lf,%lf\n"
           "        </coordinates>\n"
           "      </LineString>\n"
           "    </Placemark>\n",
           bounds[ LONGITUDE ][ MINIMUM ],
           bounds[ LATITUDE  ][ MINIMUM ],
           bounds[ LONGITUDE ][ MAXIMUM ],
           bounds[ LATITUDE  ][ MINIMUM ],
           bounds[ LONGITUDE ][ MAXIMUM ],
           bounds[ LATITUDE  ][ MAXIMUM ],
           bounds[ LONGITUDE ][ MINIMUM ],
           bounds[ LATITUDE  ][ MAXIMUM ],
           bounds[ LONGITUDE ][ MINIMUM ],
           bounds[ LATITUDE  ][ MINIMUM ] );
}



/******************************************************************************
PURPOSE: writeEndKML - Write end of KML file.
INPUTS:  FILE* file  File to write to.
******************************************************************************/

void writeEndKML( FILE* file ) {
  PRE0( file );
  const char* const fileEndingContent =
    "  </Document>\n"
    "</kml>\n";
  fprintf( file, "%s", fileEndingContent );
}



/******************************************************************************
PURPOSE: writeMapSubsetToKML - Write map as polylines to KML file.
INPUTS:  FILE* file                         KML file to write to.
         int polylineCOunt                  Number of polylines.
         const int counts[ polylineCount ]  Number of vertices per polyline.
         const float vertices[]             Lon-lat vertices.
         const char* name                   Name of polyline set placemark.
         Color color                        Color of lines.
******************************************************************************/

void writeMapSubsetToKML( FILE* file,
                          int polylineCount,
                          const int counts[],
                          const float vertices[],
                          const char* name,
                          Color color ) {

  PRE010( file,
          polylineCount > 0,
          counts,
          counts[ 0 ] > 1,
          counts[ polylineCount - 1 ] > 1,
          vertices,
          isValidLongitudeLatitude( vertices[ 0 ], vertices[ 1 ] ),
          name,
          *name,
          isValidColor( color ) );

  BGR bgr = "";
  colorKML( color, KML_COLOR_LEVELS, bgr );

  DEBUG( fprintf( stderr, "writeMapSubsetToKML: polylineCount = %d...\n",
                  polylineCount ); )

  fprintf( file,
           "    <Placemark>\n"
           "      <name>%s</name>\n"
           "      <description>%s of estuarine region of study."
           "</description>\n"
           "      <styleUrl>#ff%s</styleUrl>\n"
           "      <MultiGeometry>\n",
           name, name, bgr );
  {
    const float* longitudes = vertices;
    const float* latitudes  = vertices + 1;
    size_t polyline = 0;

    for ( polyline = 0; polyline < polylineCount; ++polyline ) {
      const size_t count = counts[ polyline ];
      size_t vertex = 0;
      CHECK( count > 1 );
      fprintf( file,
               "        <LineString>\n"
               "          <coordinates>\n" );

      for ( vertex = 0; vertex < count; ++vertex,
            longitudes += 2, latitudes += 2 ) {
        CHECK( isValidLongitudeLatitude( *longitudes, *latitudes ) );
        fprintf( file,
                 "            %f,%f\n", *longitudes, *latitudes );

      }

      fprintf( file,
               "          </coordinates>\n"
               "        </LineString>\n" );
    }

    fprintf( file,
             "      </MultiGeometry>\n"
             "    </Placemark>\n" );
  }
}



/******************************************************************************
PURPOSE: writePointsToKML - Write point data to KML.
INPUTS:  FILE* file                           KML file to write to.
         const char* const source             Data source. E.g., "STORET".
         const char* const variableName       Name of scalar variable.
         const char* const units              Units of scalar variable.
         const int timesteps                  Number of timesteps of data.
         const int hoursPerTimestep           Number of hours per timestep.
         const int timestamps[ timesteps ]    yyyymmddhh timestamps of data.
         const int stations                   Number of lon-lat stations.
         const long long id                   Optional: station ID or -1.
         const char* sids[]                   String station ids or 0.
         const float lonlats[ stations * 2 ]  Lon-lat coordinates of stations.
         const float z[ stations ]        Optional: z coordinates of stations.
         const float zAll                 Default z-coordinate for all points.
         const int components             1 = scalar, 2 or 3 is vector.
         const float dataMinimum          Minimum data value (for coloring).
         const float dataMaximum          Maximum data value (for coloring).
         const float data[ components ][ timesteps * stations ]  Variable data.
******************************************************************************/

void writePointsToKML( FILE* file,
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
                       const float data[] ) {

  PRE020( file,
          source,
          *source,
          variableName,
          *variableName,
          units,
          timesteps > 0,
          hoursPerTimestep >= 1,
          timestamps,
          stations > 0,
          lonlats,
          isValidLongitudeLatitude( lonlats[ 0 ], lonlats[ 1 ] ),
          isValidLongitudeLatitude( lonlats[ stations * 2 - 2 ],
                                    lonlats[ stations * 2 - 1 ] ),
          IMPLIES( z, AND2( ! isNan( z[ 0 ] ), ! isNan( z[ stations - 1 ] ))),
          ! isNan( zAll ),
          IN_RANGE( components, 1, 3 ),
          dataMinimum <= dataMaximum,
          data,
          ! isNan( data[ 0 ] ),
          ! isNan( data[ (long long) timesteps * (long long) stations - 1 ] ) );

  const float missing = -9999.0;
  const size_t count = timesteps * stations;
  const size_t maximumPointsToDraw = 2000; /* Per timestep. */
  const size_t stride0 = stations / maximumPointsToDraw;
  const size_t stride = stride0 > 0 ? stride0 : 1;
  const float* data1 = data;
  const float* data2 = components > 1 ? data + count : 0;
  const float* data3 = components > 2 ? data2 + count : 0;
  const char* const atOrOn = hoursPerTimestep == 24 ? "on" : "at";
  char hhmm[ 6 ] = "";
  size_t index = 0;
  size_t timestep = 0;

  if ( sids ) {
    fprintf( file,
             "    <Folder>\n"
             "      <name>Point_Data:%s_%s</name>\n"
             "      <description>Surface stations.</description>\n",
             source, variableName );
  } else if ( id >= 0 ) {
    fprintf( file,
             "    <Folder>\n"
             "      <name>Point_Data:%s_%s(%s)</name>\n"
             "      <description>Hourly point data from surface stations."
             "</description>\n",
             source, variableName, units );
  } else {
    fprintf( file,
             "    <Folder>\n"
             "      <name>Point_Data:%s_%s(%s)</name>\n"
             "      <description>Hourly point data."
             "</description>\n",
             source, variableName, units );
  }

  for ( timestep = 0, index = 0; timestep < timesteps; ++timestep ) {
    const float* longitudes = lonlats;
    const float* latitudes  = lonlats + 1;
    const int yyyymmddhh = timestamps[ timestep ];
    const int hh   = yyyymmddhh % 100;
    const int yyyymmdd = yyyymmddhh / 100;
    const int yyyy = yyyymmdd / 10000;
    const int mm   = yyyymmdd / 100 % 100;
    const int dd   = yyyymmdd % 100;
    size_t station = 0;

    if ( hoursPerTimestep != 24 ) {
      memset( hhmm, 0, sizeof hhmm );
      snprintf( hhmm, sizeof hhmm / sizeof *hhmm, "%02d:00", hh );
    }

    for ( station = 0; station < stations;
          ++station, longitudes += 2, latitudes += 2, ++index ) {
      float value = missing;
      CHECK( index < count );

      if ( station % stride == 0 ) {
        value = data1[ index ];

        if ( components == 2 ) {
          const float value2 = data2[ index ];

          if ( value2 > missing ) {
            value = sqrt( value * value + value2 * value2 );
          } else {
            value = missing;
          }
        } else if ( components == 3 ) {
          const float value2 = data2[ index ];

          if ( value2 > missing ) {
            const float value3 = data3[ index ];

            if ( value3 > missing ) {
              value =
                sqrt( value * value + value2 * value2 + value3 * value3 );
            } else {
              value = missing;
            }
          } else {
            value = missing;
          }
        }
      }

      if ( value > -999.0 ) { /* Don't output points with missing data. */
        const float longitude = *longitudes;
        const float latitude = *latitudes;
        const float pointZ = z ? z[ station ] : zAll;
        char stationId[ 64 ] = "";
        const Color color = dataColor( value, dataMinimum, dataMaximum );
        BGR bgr = "";
        colorKML( color, KML_COLOR_LEVELS, bgr );
        memset( stationId, 0, sizeof stationId );

        if ( sids ) {
          CHECK2( sids[ station ], sids[ station ][ 0 ] );
          strncpy( stationId, sids[ station ],
                   sizeof stationId / sizeof *stationId - 1 );
        } else if ( id >= 0 ) {
          snprintf( stationId, sizeof stationId / sizeof *stationId,
                    "%"INTEGER_FORMAT, id );
        }

        CHECK( strlen( stationId ) < sizeof stationId / sizeof *stationId );

        fprintf( file, "      <Placemark>\n" );

        if ( units[ 0 ] == '\0' ) {
          fprintf( file, "        <name>%d</name>\n", (int) value );

          fprintf( file,
                   "        <description>%s %s station %s at (%g, %g, %g) "
                   "%s %04d-%02d-%02d %s UTC.</description>\n",
                   source, variableName, stationId,
                   longitude, latitude, pointZ, atOrOn, yyyy, mm, dd, hhmm );
        } else if ( ! strcmp( units, "%s" ) ) {
          fprintf( file, "        <name>%s</name>\n", stationId );
          fprintf( file,
                   "        <description>%s %s station %s at (%g, %g, %g) "
                   "%s %04d-%02d-%02d %s UTC.</description>\n",
                   source, variableName, stationId,
                   longitude, latitude, pointZ, atOrOn, yyyy, mm, dd, hhmm );
        } else if ( *stationId ) {
          fprintf( file, "        <name>%5.3g(%s)</name>\n", value, units );
          fprintf( file,
                   "        <description> "
                   "%s %s (%s) station %s at (%g, %g, %g) "
                   "%s %04d-%02d-%02d %s UTC.</description>\n",
                   source, variableName, units,
                   stationId,
                   longitude, latitude, pointZ, atOrOn, yyyy, mm, dd, hhmm );
        } else {
          fprintf( file, "        <name>%5.3g(%s)</name>\n", value, units );
          fprintf( file,
                   "        <description> "
                   "%s %s (%s) at (%g, %g, %g) "
                   "%s %04d-%02d-%02d %s UTC.</description>\n",
                   source, variableName, units,
                   longitude, latitude, pointZ, atOrOn, yyyy, mm, dd, hhmm );
        }

        fprintf( file,
                 "        <TimeSpan>\n"
                 "          <begin>%04d-%02d-%02dT%02d:00:00Z</begin>\n"
                 "          <end>%04d-%02d-%02dT%02d:59:59Z</end>\n"
                 "        </TimeSpan>\n"
                 "        <styleUrl>#ff%s</styleUrl>\n",
                 yyyy, mm, dd, hh,
                 yyyy, mm, dd,
                 hoursPerTimestep == 1 ? hh : hoursPerTimestep - 1,
                 bgr );

        if ( components == 1 ) { /* Draw a ruler icon: */
          fprintf( file,
                   "        <Point>\n"
                   "          <coordinates>%f,%f,%f</coordinates>\n"
                   "        </Point>\n"
                   "      </Placemark>\n",
                   longitude, latitude, pointZ );
        } else { /* Draw an arrow vector: */
          const double degreesPerPixel = 0.001;
          const double pixelsPerUnitLength =
            ! strcmp( variableName, "wind" ) ? 5.0 : 20.0;
          double point0X = 0.0;
          double point0Y = 0.0;
          double point1X = 0.0;
          double point1Y = 0.0;
          double point2X = 0.0;
          double point2Y = 0.0;
          double point3X = 0.0;
          double point3Y = 0.0;
          const double z0 = components > 2 ? data3[ index ] : pointZ;
          computeArrowVectorCoordinates( longitude, latitude,
                                         data1[ index ], data2[ index ],
                                         degreesPerPixel,
                                         pixelsPerUnitLength,
                                         &point0X,
                                         &point0Y,
                                         &point1X,
                                         &point1Y,
                                         &point2X,
                                         &point2Y,
                                         &point3X,
                                         &point3Y );
          fprintf( file,
                   "        <LineString>\n"
                   "          <coordinates>\n"
                   "            %lf,%lf,%lf\n"
                   "            %lf,%lf,%lf\n"
                   "            %lf,%lf,%lf\n"
                   "            %lf,%lf,%lf\n"
                   "            %lf,%lf,%lf\n"
                   "          </coordinates>\n"
                   "        </LineString>\n"
                   "      </Placemark>\n",
                   point0X, point0Y, z0,
                   point1X, point1Y, z0,
                   point2X, point2Y, z0,
                   point3X, point3Y, z0,
                   point1X, point1Y, z0 );
        }
      }
    }
  }

  fprintf( file,
           "    </Folder>\n" );

}



/******************************************************************************
PURPOSE: writeShapeDataToKML - Write point data to KML.
INPUTS:  FILE* file                  KML file to write to.
         const char* name            Name of data variable.
         const char* units           Units of data variable.
         const char* description     Description of data.
         int column                  Column index of data to write.
         double minimumValue         Minimum value for coloring.
         double maximumValue         Maximum value for coloring.
         DataColor dataColor         Function for coloring.
         TextColor textColor         Function for coloring.
         const ShapeData* shapeData  Point data info.
******************************************************************************/

void writeShapeDataToKML( FILE* file, const char* name, const char* units,
                          const char* description,
                          int column, double minimumValue, double maximumValue,
                          DataColor dataColor, TextColor textColor,
                          const ShapeData* shapeData ) {

  PRE011( file, name, *name, units, *units, column >= 0, shapeData,
          shapeData->columns >= 4,
          column < shapeData->columns,
          shapeData->columnTypes,
          IMPLIES_ELSE( shapeData->columnTypes[ column ] == FTString,
                        textColor, dataColor ) );

  const size_t rows    = shapeData->rows;
  const size_t columns = shapeData->columns;
  const int type = shapeData->columnTypes[ column ];
  const Value* const values = shapeData->values;
  const int longitudeColumn =
    indexOfString("LONGITUDE", (const char**) shapeData->columnNames, columns);
  const int latitudeColumn =
    indexOfString( "LATITUDE", (const char**) shapeData->columnNames, columns);
  const int dateColumn =
    indexOfString( "DATE", (const char**) shapeData->columnNames, columns );
  const int depthColumn =
    indexOfString("WATERDEPTH", (const char**) shapeData->columnNames,columns);
  size_t row = 0;

  CHECK2( IN_RANGE( longitudeColumn, 2, shapeData->columns - 1 ),
          IN_RANGE( latitudeColumn,  2, shapeData->columns - 1 ) );

  fprintf( file,
           "    <Folder>\n"
           "     <name>Point_Data:%s(%s)</name>\n"
           "        <description>%s</description>\n",
           name, units, description );

  for ( row = 0; row < rows; ++row ) {
    const size_t rowOffset = row * columns;
    const Value* const rowValues = values + rowOffset;
    const Value value = rowValues[ column ];

    if ( isValidValue( type, units, value ) ) { /* Don't output invalid. */
      const double longitude = rowValues[ longitudeColumn ].d;
      const double latitude  = rowValues[ latitudeColumn ].d;
      char siteLabel[ 64 ] = "";
      char valueLabel[ 128 ] = "";
      char dateLabel[ 32 ] = "";
      char depthLabel[ 32 ] = "";
      const Color color =
        type == FTString ? textColor( value.s )
        : type == FTDouble ? dataColor( value.d, minimumValue, maximumValue )
        : dataColor( (double) value.i, minimumValue, maximumValue );
      BGR bgr = "";
      colorKML( color, KML_COLOR_LEVELS, bgr );


      /* Site label: */

      if ( AND3( columns >= 2,
                 shapeData->columnTypes[ 0 ] == FTString,
                 shapeData->columnTypes[ 1 ] == FTInteger ) ) {
        CHECK( rowValues[ 0 ].s );
        snprintf( siteLabel, sizeof siteLabel / sizeof *siteLabel, "%s.%d",
                 rowValues[ 0 ].s, rowValues[ 1 ].i );
      } else if ( AND3( columns >= 2,
                        shapeData->columnTypes[ 0 ] == FTInteger,
                        shapeData->columnTypes[ 1 ] == FTInteger ) ) {
        snprintf( siteLabel, sizeof siteLabel / sizeof *siteLabel, "%d.%d",
                 rowValues[ 0 ].i, rowValues[ 1 ].i );
      } else if ( AND2( columns >= 2,
                        shapeData->columnTypes[ 1 ] == FTString ) ) {
        CHECK( rowValues[ 1 ].s );
        snprintf( siteLabel, sizeof siteLabel / sizeof *siteLabel, "%s",
                 rowValues[ 1 ].s );
      }

      CHECK( strlen( siteLabel ) < sizeof siteLabel / sizeof *siteLabel );

      /* Value label: */

      if ( type == FTString ) {
        CHECK( value.s );
        snprintf( valueLabel, sizeof valueLabel / sizeof *valueLabel,
                  "%s(%s)", value.s, units );
      } else if ( type == FTInteger ) {
        snprintf( valueLabel, sizeof valueLabel / sizeof *valueLabel,
                  "%d(%s)", value.i, units );
      } else {
        CHECK( type == FTDouble );
        snprintf( valueLabel, sizeof valueLabel / sizeof *valueLabel,
                  "%5.3g(%s)", value.d, units );
      }

      CHECK( strlen( valueLabel ) < sizeof valueLabel / sizeof *valueLabel );

      /* Date label: */

      if ( IN_RANGE( dateColumn, 0, columns - 1 ) ) {
        const int yyyymmdd = rowValues[ dateColumn ].i;
        const int yyyy = yyyymmdd / 10000;
        const int mm   = yyyymmdd / 100 % 100;
        const int dd   = yyyymmdd % 100;
        snprintf( dateLabel, sizeof dateLabel / sizeof *dateLabel,
                  " on %04d-%02d-%02d.", yyyy, mm, dd );
      }

      CHECK( strlen( dateLabel ) < sizeof dateLabel / sizeof *dateLabel );

      /* Depth label: */

      if ( IN_RANGE( depthColumn, 0, columns - 1 ) ) {
        CHECK( shapeData->columnTypes[ depthColumn ] == FTDouble );
        snprintf( depthLabel, sizeof depthLabel / sizeof *depthLabel,
                  ",-%lf", rowValues[ depthColumn ].d );
      }

      CHECK( strlen( depthLabel ) < sizeof depthLabel / sizeof *depthLabel );

      fprintf( file,
               "      <Placemark>\n"
               "        <name>%s</name>\n"
               "        <description>%s (%s) at site %s%s."
               "</description>\n"
               "        <styleUrl>#ff%s</styleUrl>\n"
               "        <Point>\n"
               "          <coordinates>%lf,%lf%s</coordinates>\n"
               "        </Point>\n"
               "      </Placemark>\n",
               valueLabel,
               name, units, siteLabel, dateLabel,
               bgr,
               longitude, latitude, depthLabel );
    }
  }

  fprintf( file, "    </Folder>\n" );
}



/******************************************************************************
PURPOSE: writeGridToKML - Write grid data to KML.
INPUTS:  FILE* file                  File to write to.
         const int timestamps[ timesteps ] Optional: timestamps (yyyymmddhh)
                                     to use instead of yyyymmdd.
         int yyyymmdd                Starting timestamp (at hour 0)
                                     if yyyymmddhh is not 0.
         int components              1 = scalar, 2 or 3 is vector.
         int timesteps               Number of hourly timesteps of data.
         int rows                    Number of rows of grid cells.
         int columns                 Number of columns of grid cells.
         double westEdge             Meters from center of projection to
                                     west edge of grid.
         double southEdge            Meters from center of projection to
                                     south edge of grid.
         double cellWidth            Width of grid cell in meters.
         double cellHeight           Height of grid cell in meters.
         const LongName source       E.g., "CMAQ".
         const Name name             E.g., "NOx".
         const Units units           E.g., "ug/ha".
         const float corners[]       corners[ (rows + 1 ) * (columns + 1) ].
                                     or 0 if regular lon-lat grid.
         int type                    Grid cell data type: FLOAT_TYPE, etc.
         const void* data      data[ components * timesteps * rows * columns ]
         const float dataRange[ 2 ]  Minimum and maximum data values.
         const RGBColormap colormap  Colormap to use, else 0.
         const DataColor dataColor   DataColor function to use, else 0.
******************************************************************************/

void writeGridToKML( FILE* file,
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
                     DataColor dataColor ) {

  PRE019( file, timesteps > 0, rows > 0, columns > 0,
          IN_RANGE( components, 1, 3 ),
          cellWidth > 0.0, cellHeight > 0.0,
          IMPLIES( corners == 0,
                   AND2( isValidLongitudeLatitude( westEdge, southEdge ),
                         isValidLongitudeLatitude( westEdge +
                                                   columns * cellWidth,
                                                   southEdge +
                                                   rows * cellHeight ) ) ),
          source, *source, name, *name, units, *units,
          IS_VALID_GRID_DATA_TYPE( type ),
          data, dataRange, dataRange[ MINIMUM ] <= dataRange[ MAXIMUM ],
          XOR2( colormap, dataColor ) );

  const float missing = -9999.0;
  const double degreesPerPixel = 0.001;
  const double pixelsPerUnitLength = ! strcmp( name, "wind" ) ? 5.0 : 20.0;
  const unsigned char opacity = 127; /* 50% transparent so ground is seen. */
  const float* const fdata = type == FLOAT_TYPE ? data : 0;
  const char* const cdata = type == BYTE_TYPE ? data : 0;
  const unsigned short* const sdata = type == UINT16_TYPE ? data : 0;
  const float dataMinimum = dataRange[ MINIMUM ];
  const float dataMaximum = dataRange[ MAXIMUM ];
  const size_t columnsPlus1 = columns + 1;
  const size_t nextRowOffset = columnsPlus1 + columnsPlus1;
  int yyyymmddhh = timestamps ? timestamps[ 0 ] : yyyymmdd * 100;
  size_t timestep = 0;
  const size_t componentSize = timesteps * rows * columns;
  CHECK( IMPLIES( components > 1, fdata ) );

  fprintf( file,
           "    <Folder>\n"
           "      <name>Grid_Data:%s_%s_(%s)</name>\n"
           "      <description>Gridded data.</description>\n",
           source, name, units );

  if ( colormap ) {
    writeKMLColormap( file, colormap, opacity );
  }

  for ( timestep = 0; timestep < timesteps; ++timestep,
        yyyymmddhh = AND2( timestamps, timestep < timesteps ) ?
          timestamps[ timestep ]
        : incrementDateTime( yyyymmddhh, 1 ) ) {
    const int hh   = yyyymmddhh % 100;
    const int yyyymmdd = yyyymmddhh / 100;
    const int yyyy = yyyymmdd / 10000;
    const int mm   = yyyymmdd / 100 % 100;
    const int dd   = yyyymmdd % 100;
    const size_t timestepOffset = timestep * rows * columns;
    size_t row = 0;
    size_t index = 0;
    size_t indexSW = 0;

    /* Draw colored label value (units) in center of cell: */

    for ( row = 0, index = timestepOffset, indexSW = 0;
          row < rows;
          ++row, indexSW += 2 ) {
      size_t column = 0;

      for ( column = 0; column < columns;
           ++column, ++index, indexSW += 2 ) {
        float value =
          fdata ? fdata[ index ]
          : cdata ? (float) ( cdata[ index ] )
          : sdata ? (float) ( sdata[ index ] )
          : -9999.0;
        float u = 0.0;
        float v = 0.0;
        int drawThisCell = 1;
        Color color = { 0.0, 0.0, 0.0 };

        if ( colormap ) {
          drawThisCell = IN_RANGE( value, 0.0, 127.0 );

          if ( drawThisCell ) {
            const unsigned char ivalue = (unsigned char) value;
            const RGB rgb = colormap[ ivalue ];
            drawThisCell = ! IS_ZERO3( rgb.r, rgb.g, rgb.b );

            if ( drawThisCell ) {
              color.r = rgb.r / 255.0;
              color.g = rgb.g / 255.0;
              color.b = rgb.b / 255.0;
              units = rgb.s; /* Land use category, e.g., "Evergreen Forest".*/
            }
          }
        } else {

          if ( components > 1 ) {
            const size_t cellOffset = row * columns + column;
            size_t component = 0;

            for ( component = 0; component < components; ++component ) {
              const size_t index =
                component * componentSize + timestepOffset + cellOffset;
              const float componentValue = fdata[ index ];

              if ( component == 0 ) {
                u = componentValue;
              } else if ( component == 1 ) {
                v = componentValue;
              }

              if ( componentValue > missing ) {

                if ( value == missing ) {
                  value = componentValue * componentValue;
                } else {
                  value += componentValue * componentValue;
                }
              }
            }

            if ( value >= 0.0 ) {
              value = sqrt( value );
            }
          }

          color = dataColor( value, dataMinimum, dataMaximum );
          drawThisCell = ! IS_ZERO3( color.r, color.g, color.b );
        }

        if ( drawThisCell ) {
          BGR bgr = "";
          colorKML( color, KML_COLOR_LEVELS, bgr );
          float longitudeSW = 0.0;
          float latitudeSW  = 0.0;
          float longitudeSE = 0.0;
          float latitudeSE  = 0.0;
          float longitudeNE = 0.0;
          float latitudeNE  = 0.0;
          float longitudeNW = 0.0;
          float latitudeNW  = 0.0;

          if ( corners ) {
            const size_t indexSE = indexSW + 2;
            const size_t indexNE = indexSE + nextRowOffset;
            const size_t indexNW = indexNE - 2;
            longitudeSW = corners[ indexSW ];
            latitudeSW  = corners[ indexSW + 1 ];
            longitudeSE = corners[ indexSE ];
            latitudeSE  = corners[ indexSE + 1 ];
            longitudeNE = corners[ indexNE ];
            latitudeNE  = corners[ indexNE + 1 ];
            longitudeNW = corners[ indexNW ];
            latitudeNW  = corners[ indexNW + 1 ];
          } else {
            longitudeSW = longitudeNW = westEdge + column * cellWidth;
            longitudeSE = longitudeNE = longitudeSW + cellWidth;
            latitudeSE  = latitudeSW = southEdge + row * cellHeight;
            latitudeNW  = latitudeNE = latitudeSE + cellHeight;
          }

          fprintf( file,
                   "      <Placemark>\n"
                   "        <name>%5.3g(%s)</name>\n"
                   "        <description>%s %s (%s) at "
                   "%04d-%02d-%02d %02d:00 UTC.</description>\n"
                   "        <TimeSpan>\n"
                   "          <begin>%04d-%02d-%02dT%02d:00:00Z</begin>\n"
                   "          <end>%04d-%02d-%02dT%02d:59:59Z</end>\n"
                   "        </TimeSpan>\n"
                   "        <styleUrl>#%02x%s</styleUrl>\n",
                  value, units,
                  source, name, units,
                  yyyy, mm, dd, hh,
                  yyyy, mm, dd, hh,
                  yyyy, mm, dd, hh,
                  opacity, bgr );

          if ( components == 1 ) { /* Draw colored quad grid cell: */
            fprintf( file,
                     "        <Style>\n"
                     "          <PolyStyle>\n"
                     "            <outline>0</outline>\n"
                     "          </PolyStyle>\n"
                     "        </Style>\n"
                     "        <Polygon>\n"
                     "          <outerBoundaryIs>\n"
                     "            <LinearRing>\n"
                     "              <coordinates>\n"
                     "                %f,%f\n"
                     "                %f,%f\n"
                     "                %f,%f\n"
                     "                %f,%f\n"
                     "                %f,%f\n"
                     "              </coordinates>\n"
                     "            </LinearRing>\n"
                     "          </outerBoundaryIs>\n"
                     "        </Polygon>\n"
                     "      </Placemark>\n",
                    longitudeSW, latitudeSW,
                    longitudeSE, latitudeSE,
                    longitudeNE, latitudeNE,
                    longitudeNW, latitudeNW,
                    longitudeSW, latitudeSW );
          } else { /* Draw colored arrow vector: */
            const float longitude =
              westEdge + column * cellWidth + cellWidth * 0.5;
            const float latitude =
              southEdge + row * cellHeight + cellHeight * 0.5;
            double point0X = 0.0;
            double point0Y = 0.0;
            double point1X = 0.0;
            double point1Y = 0.0;
            double point2X = 0.0;
            double point2Y = 0.0;
            double point3X = 0.0;
            double point3Y = 0.0;
            const size_t cellOffset = row * columns + column;
            const size_t index3 =
              components == 3 ? 2 * componentSize + cellOffset : 0;
            const float value3 = components == 3 ? fdata[ index3 ] : 0.0;
            const float z0 = value3 > missing ? value3 : 0.0;
            computeArrowVectorCoordinates( longitude, latitude,
                                           u, v,
                                           degreesPerPixel,
                                           pixelsPerUnitLength,
                                           &point0X,
                                           &point0Y,
                                           &point1X,
                                           &point1Y,
                                           &point2X,
                                           &point2Y,
                                           &point3X,
                                           &point3Y );
            fprintf( file,
                     "        <LineString>\n"
                     "          <coordinates>\n"
                     "            %lf,%lf,%lf\n"
                     "            %lf,%lf,%lf\n"
                     "            %lf,%lf,%lf\n"
                     "            %lf,%lf,%lf\n"
                     "            %lf,%lf,%lf\n"
                     "          </coordinates>\n"
                     "        </LineString>\n"
                     "      </Placemark>\n",
                     point0X, point0Y, z0,
                     point1X, point1Y, z0,
                     point2X, point2Y, z0,
                     point3X, point3Y, z0,
                     point1X, point1Y, z0 );
          }
        }
      }
    }
  }

  fprintf( file,
           "    </Folder>\n" );
}


/******************************************************************************
PURPOSE: writePolygonsToKML - Write polygon data to KML.
INPUTS:  FILE* file                       File to write to.
         const int yyyymmddhhStart        Starting timestamp or 0.
         const int yyyymmddhhEnd          Ending timestamp or 0.
         const int count                        Number of polygons.
         const PolygonShape polygons[ count ]  Polygons to write.
         const ShapeData* shapeData       ShapeData for polygons.
         const char* stringValues[count]  If not 0 then use instead of
                                          csvValues or shapeData->values.
         const double csvValues[count]    If not 0 then use instead of
                                          shapeData->values.
         const LongName source            E.g., "PRISM".
         const Name name                  E.g., "TEMP_C".
         const Units units                E.g., "degC".
         const double dataRange[ 2 ]      Minimum and maximum data values.
         const RGBColormap* colormap      Colormap to use, else 0.
         const DataColor dataColor        DataColor function to use, else 0.
         const TextColor textColor        TextColor function to use, else 0.
******************************************************************************/

void writePolygonsToKML( FILE* file,
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
                         TextColor textColor ) {

  PRE015( file,
          OR2( IS_ZERO2( yyyymmddhhStart, yyyymmddhhEnd ),
               AND3( isValidYYYYMMDDHH( yyyymmddhhStart ),
                     isValidYYYYMMDDHH( yyyymmddhhEnd ),
                     yyyymmddhhStart <= yyyymmddhhEnd ) ),
          count > 0, polygons, shapeData, source, *source, name, *name,
          OR3( stringValues, csvValues,
               indexOfString( name, (const char**) shapeData->columnNames,
                              shapeData->columns ) > -1 ),
          units, *units, ! isNan( dataRange[ 0 ] ), ! isNan( dataRange[ 1 ] ),
          dataRange[ MINIMUM ] <= dataRange[ MAXIMUM ] );

  const unsigned char opacity = 127; /* 50% transparent so ground is seen. */
  const double dataMinimum = dataRange[ MINIMUM ];
  const double dataMaximum = dataRange[ MAXIMUM ];
  const int columns = shapeData->columns;
  const int column =
    OR2( stringValues, csvValues ) ? 0
    : indexOfString( name, (const char**) shapeData->columnNames, columns );
  const int dataType = csvValues ? FTDouble : shapeData->columnTypes[ column ];
  const Value* const values = shapeData->values;
  const int isPolyline = polygons[ 0 ].triangles.num_strips == 0;
  const char* const geoTypeName = isPolyline ? "Polyline" : "Polygon";
  const int yyyy1 = yyyymmddhhStart / 1000000;
  const int mm1   = yyyymmddhhStart / 10000 % 100;
  const int dd1   = yyyymmddhhStart / 100 % 100;
  const int hh1   = yyyymmddhhStart % 100;
  const int yyyy2 = yyyymmddhhEnd / 1000000;
  const int mm2   = yyyymmddhhEnd / 10000 % 100;
  const int dd2   = yyyymmddhhEnd / 100 % 100;
  const int hh2   = yyyymmddhhEnd % 100;
  int set = 0;
  int index = 0;

  CHECK( IMPLIES_ELSE( csvValues,
                       dataType == FTDouble,
                       OR2( stringValues,
                            IN_RANGE( column, 0, columns - 1 ) ) ) );

  fprintf( file,
           "    <Folder>\n"
           "      <name>%s_Data:%s_%s_(%s)</name>\n"
           "      <description>%s data.</description>\n",
           geoTypeName, source, name, units, geoTypeName );

  if ( colormap ) {
    writeKMLColormap( file, *colormap, opacity );
  }

  for ( set = 0, index = column; set < count; ++set, index += columns ) {
    const PolygonShape* const polygonShape = polygons + set;
    const gpc_vertex_list* const vertexList =
      isPolyline ? polygonShape->polygon.contour
      : polygonShape->triangles.strip;
    const int strips =
      isPolyline ? polygonShape->polygon.num_contours
      : polygonShape->triangles.num_strips;
    const char* const stringValue =
      stringValues ? stringValues[ set ] : &name[ 0 ];
    const double value =
      csvValues ? csvValues[ set ]
      : dataType == FTDouble ? values[ index ].d
      : dataType == FTInteger ? (double) values[ index ].i
      : 0.0;
    Color color = { 0.0, 0.0, 0.0 };
    BGR bgr = "";
    int strip = 0;

    if ( dataType == FTString ) {
      color = textColor( values[ index ].s );
    } else if ( colormap ) {
      const unsigned char ivalue = (unsigned char) value;
      const RGB rgb = (*colormap)[ ivalue ];
      color.r = rgb.r / 255.0;
      color.g = rgb.g / 255.0;
      color.b = rgb.b / 255.0;
    } else {

      if ( textColor ) {
        color = textColor( name );
      } else if ( dataColor ) {
        color = dataColor( value, dataMinimum, dataMaximum );
      }
    }

    if ( ! IS_ZERO3( color.r, color.g, color.b ) ) { /* Non-missing data: */
      colorKML( color, KML_COLOR_LEVELS, bgr );
      fprintf( file,
              "      <Placemark>\n" );

      if ( dataType == FTString ) {
        const char* const textValue =
          stringValues ? stringValue : values[ index ].s;
        fprintf( file,
                 "        <name>%s</name>\n", textValue );
      } else if ( dataType == FTInteger ) {
        fprintf( file,
                 "        <name>%d(%s)</name>\n", values[ index ].i, units );
      } else {
        fprintf( file,
                 "        <name>%5.3lg(%s)</name>\n", value, units );
      }

      if ( yyyymmddhhStart ) {
        fprintf( file,
                 "        <description>%s %s (%s) at "
                 "%04d-%02d-%02d %02d:00 UTC.</description>\n"
                 "        <TimeSpan>\n"
                 "          <begin>%04d-%02d-%02dT%02d:00:00Z</begin>\n"
                 "          <end>%04d-%02d-%02dT%02d:59:59Z</end>\n"
                 "        </TimeSpan>\n",
                 source, name, units,
                 yyyy1, mm1, dd1, hh1,
                 yyyy1, mm1, dd1, hh1,
                 yyyy2, mm2, dd2, hh2 );
      } else {
        fprintf( file,
                 "        <description>%s %s (%s).</description>\n",
                 source, name, units );
      }

      fprintf( file,
               "        <styleUrl>#%02x%s</styleUrl>\n"
               "        <Style>\n"
               "          <PolyStyle>\n"
               "            <outline>0</outline>\n"
               "          </PolyStyle>\n"
               "        </Style>\n"
               "        <MultiGeometry>\n",
               opacity, bgr );

      /*
       * For each triangle strip, draw it as a concave polygon by
       * first drawing one edge (even vertices) then drawing the
       * odd vertices starting at the end back to the first vertex.
       */

      for ( strip = 0; strip < strips; ++strip ) {
        const int vertexCount = vertexList[ strip ].num_vertices;
        const gpc_vertex* const v = vertexList[ strip ].vertex;
        int vertex = 0;
        const int lastOddVertex =
          vertexCount % 2 == 0 ? vertexCount - 1 : vertexCount - 2;

        if ( isPolyline ) {
          fprintf( file,
                   "          <LineString>\n"
                   "            <coordinates>\n" );
        } else {
          fprintf( file,
                   "          <Polygon>\n"
                   "            <outerBoundaryIs>\n"
                   "              <LinearRing>\n"
                   "                <coordinates>\n" );
        }

        for ( vertex = 0; vertex < vertexCount; vertex += 2 ) {
          const double x = v[ vertex ].x;
          const double y = v[ vertex ].y;
          fprintf( file, "                %lf,%lf\n", x, y );
        }

        for ( vertex = lastOddVertex; vertex >= 0; vertex -= 2 ) {
          const double x = v[ vertex ].x;
          const double y = v[ vertex ].y;
          fprintf( file, "                %lf,%lf\n", x, y );
        }

        if ( vertex != 0 ) { /* End with first vertex to close polygon: */
          const double x = v[ 0 ].x;
          const double y = v[ 0 ].y;
          fprintf( file, "                %lf,%lf\n", x, y );
        }

        if ( isPolyline ) {
          fprintf( file,
                   "            </coordinates>\n"
                   "          </LineString>\n" );          
        } else {
          fprintf( file,
                   "                </coordinates>\n"
                   "              </LinearRing>\n"
                   "            </outerBoundaryIs>\n"
                   "          </Polygon>\n" );          
        }
      }

      fprintf( file,
               "        </MultiGeometry>\n"
               "      </Placemark>\n" );
    } /* Non-black (missing data) color. */
  }

  fprintf( file, "    </Folder>\n" );
}



/******************************************************************************
PURPOSE: writeMapImagesToKML - Write map_image_yyyymmdd.png files to KML.
INPUTS:  FILE* file                   File to write to.
         const char* const inputDirectory   Directory path to check for
                                            map_image_yyyymmdd.png files.
         const char* const outputDirectory  Directory path for outputing files.
         const Bounds bounds                Lon-lat bounds of images.
******************************************************************************/

void writeMapImagesToKML( FILE* file,
                          const char* const inputDirectory,
                          const char* const outputDirectory,
                          const Bounds bounds ) {

  PRE07( file, inputDirectory, *inputDirectory,
         outputDirectory, *outputDirectory,
         strcmp( inputDirectory, outputDirectory ),
         isValidBounds( bounds ) );

  const char* const startsWith = "map_image_";
  const char* const endsWith = ".png";
  MapImageFileInfo info;
  memset( &info, 0, sizeof info );
  info.file = file;
  info.directory = outputDirectory;
  info.bounds[ LONGITUDE ][ MINIMUM ] = bounds[ LONGITUDE ][ MINIMUM ];
  info.bounds[ LONGITUDE ][ MAXIMUM ] = bounds[ LONGITUDE ][ MAXIMUM ];
  info.bounds[ LATITUDE  ][ MINIMUM ] = bounds[ LATITUDE  ][ MINIMUM ];
  info.bounds[ LATITUDE  ][ MAXIMUM ] = bounds[ LATITUDE  ][ MAXIMUM ];

  fprintf( file, "    <Folder>\n" );
  fprintf( file, "      <name>Satellite images</name>\n" );
  fprintf( file, "      <description>Satellite images."
                       "</description>\n" );

  /* Calls writeMapImagesToKMLHelper() with each matched file name: */

  forEachFile( inputDirectory, startsWith, endsWith,
               writeMapImagesToKMLHelper, (void*) &info );

  fprintf( file, "    </Folder>\n" );
}



/*============================= PRIVATE FUNCTIONS ===========================*/



/******************************************************************************
PURPOSE: isValidMapImageFileInfo - Check MapImageFileInfo object.
INPUTS:  const MapImageFileInfo* const info  Object to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

static int isValidMapImageFileInfo( const MapImageFileInfo* const info ) {
  const int result = 
    AND5( info, info->file, info->directory, info->directory[ 0 ],
          isValidBounds( info->bounds ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeMapImagesToKMLHelper - Write map_image_yyyymmdd.png file to KML.
INPUTS:  const char* const fileName   File name map_image_yyyymmdd.png.
         void* data                   MapImageFileInfo.
******************************************************************************/

static void writeMapImagesToKMLHelper( const char* const fileName,
                                       void* data ) {

  PRE04( fileName, *fileName, data,
         isValidMapImageFileInfo( (const MapImageFileInfo*) data ) );

  const MapImageFileInfo* const info = (const MapImageFileInfo*) data;
  FILE* const file = info->file;
  const char* const directory = info->directory;
  const double west  = info->bounds[ LONGITUDE ][ MINIMUM ];
  const double east  = info->bounds[ LONGITUDE ][ MAXIMUM ];
  const double south = info->bounds[ LATITUDE  ][ MINIMUM ];
  const double north = info->bounds[ LATITUDE  ][ MAXIMUM ];
  const char* const startsWith = "map_image_";
  const size_t startLength = strlen( startsWith );
  const size_t fileNameLength = strlen( fileName );

  if ( fileNameLength > startLength ) {
    const char* const yyyymmddString = fileName + startLength;
    const int yyyymmdd = (int) atoi( yyyymmddString );
    const int yyyy = yyyymmdd / 10000;
    const int mm   = yyyymmdd / 100 % 100;
    const int dd   = yyyymmdd % 100;

    if ( isValidYearMonthDay( yyyymmdd ) ) {
      fprintf( file, "      <GroundOverlay>\n" );
      fprintf( file, "        <name>image_%d</name>\n", yyyymmdd );
      fprintf( file, "        <description>Satellite image on "
               "%d-%02d-%02d UTC.</description>\n", yyyy, mm, dd );
      fprintf( file, "        <TimeSpan>\n" );
      fprintf( file, "          <begin>%d-%02d-%02dT00:00:00Z</begin>\n",
               yyyy, mm, dd );
      fprintf( file, "          <end>%d-%02d-%02dT23:59:59Z</end>\n",
               yyyy, mm, dd );
      fprintf( file, "        </TimeSpan>\n" );
      fprintf( file, "        <Icon>\n" );
      fprintf( file, "          <href>%s/map_image_%d.png</href>\n",
               directory, yyyymmdd );
      fprintf( file, "        </Icon>\n" );
      fprintf( file, "        <LatLonBox>\n" );
      fprintf( file, "          <north>%lf</north>\n", north );
      fprintf( file, "          <south>%lf</south>\n", south );
      fprintf( file, "          <east>%lf</east>\n", east );
      fprintf( file, "          <west>%lf</west>\n", west );
      fprintf( file, "        </LatLonBox>\n" );
      fprintf( file, "      </GroundOverlay>\n" );
    }
  }
}



/******************************************************************************
PURPOSE: writeKMLColors - Write color styles to KML file.
INPUTS:  FILE* file  KML file to write to.
int levels  Number of levels per color component.
******************************************************************************/

static void writeKMLColors( FILE* file, int levels ) {

  PRE02( file, 256 % levels == 0 );

  const int increment = 256 / levels;
  int alpha = 0;

  for ( alpha = 127; alpha < 256; alpha += 128 ) {
    int blue = 0;

    for ( blue = 0; blue <= 256; blue += increment ) {
      int green = 0;
      blue = blue == 256 ? 255 : blue;

      for ( green = 0; green <= 256; green += increment ) {
        int red = 0;
        green = green == 256 ? 255 : green;

        for ( red = 0; red <= 256; red += increment ) {
          red = red == 256 ? 255 : red;

          fprintf( file,
                  "    <Style id=\"%02x%02x%02x%02x\">\n"
                  "      <IconStyle>\n"
                  "        <color>%02x%02x%02x%02x</color>\n"
                  "        <Icon>\n"
                  "          <href>"
                  "http://maps.google.com/mapfiles/kml/pal5/icon5.png"
                  "</href>\n"
                  "        </Icon>\n"
                  "      </IconStyle>\n"
                  "      <LabelStyle>\n"
                  "        <color>%02x%02x%02x%02x</color>\n"
                  "      </LabelStyle>\n"
                  "      <LineStyle>\n"
                  "        <color>%02x%02x%02x%02x</color>\n"
                  "        <width>2</width>\n"
                  "      </LineStyle>\n"
                  "      <PolyStyle>\n"
                  "        <color>%02x%02x%02x%02x</color>\n"
                  "      </PolyStyle>\n"
                  "    </Style>\n",
                  alpha, blue, green, red,
                  alpha, blue, green, red,
                  alpha, blue, green, red,
                  alpha, blue, green, red,
                  alpha, blue, green, red );
        }
      }
    }
  }
}



/******************************************************************************
PURPOSE: writeKMLColormap - Write colormap styles to KML file.
INPUTS:  FILE* file                  KML file to write to.
         const RGBColormap colormap  Colormap to write.
         unsigned char opacity       Opacity for all colormap entries.
******************************************************************************/

static void writeKMLColormap( FILE* file, const RGBColormap colormap,
                              unsigned char opacity ) {

  PRE02( file, colormap );

  int index = 0;

#if 0
  <ScreenOverlay>
    <name>Legend</name>
    <color>bfffffff</color>
    <Icon>
      <href>rsig2dviz_18751_legend.png</href>
    </Icon>
    <overlayXY  x="0" y="0" xunits="fraction" yunits="fraction"/>
    <screenXY   x="0" y="0" xunits="fraction" yunits="fraction"/>
    <rotationXY x="0" y="0" xunits="fraction" yunits="fraction"/>
    <size       x="0" y="0" xunits="fraction" yunits="fraction"/>
  </ScreenOverlay>
#endif

  for ( index = 0; index < 255; ++index ) {
    const RGB rgb = colormap[ index ];

    if ( ! IS_ZERO3( rgb.r, rgb.g, rgb.b ) ) {
      fprintf( file,
               "    <Style id=\"%02x%02x%02x%02x\">\n"
               "      <IconStyle>\n"
               "        <color>%02x%02x%02x%02x</color>\n"
               "        <Icon>\n"
               "          <href>"
               "http://maps.google.com/mapfiles/kml/pal5/icon5.png"
               "</href>\n"
               "        </Icon>\n"
               "      </IconStyle>\n"
               "      <LabelStyle>\n"
               "        <color>%02x%02x%02x%02x</color>\n"
               "      </LabelStyle>\n"
               "      <LineStyle>\n"
               "        <color>%02x%02x%02x%02x</color>\n"
               "        <width>2</width>\n"
               "      </LineStyle>\n"
               "      <PolyStyle>\n"
               "        <color>%02x%02x%02x%02x</color>\n"
               "      </PolyStyle>\n"
               "    </Style>\n",
               opacity, rgb.b, rgb.g, rgb.r,
               opacity, rgb.b, rgb.g, rgb.r,
               opacity, rgb.b, rgb.g, rgb.r,
               opacity, rgb.b, rgb.g, rgb.r,
               opacity, rgb.b, rgb.g, rgb.r );
    }
  }
}



/******************************************************************************
PURPOSE: colorKML - Compute KML color string from Color.
INPUTS:  Color color  Color to convert to KML.
         int levels   Number of levels per color component.
OUTPUTS: BGR bgr      KML bblue-green-red color string.
******************************************************************************/

static void colorKML( Color color, int levels, BGR result ) {
  PRE02( isValidColor( color ), result );
  const int increment = 256 / levels;
  int b = (int) ( color.b * 255 );
  int g = (int) ( color.g * 255 );
  int r = (int) ( color.r * 255 );
  b = (int) ( (float) b / increment + 0.5 );
  b *= increment;
  b = b < 256 ? b : 255;
  g = (int) ( (float) g / increment + 0.5 );
  g *= increment;
  g = g < 256 ? g : 255;
  r = (int) ( (float) r / increment + 0.5 );
  r *= increment;
  r = r < 256 ? r : 255;
  snprintf( result, sizeof (BGR) / sizeof (char), "%02x%02x%02x", b, g, r );
  POST03( strlen( result ) == 6,
          strcmp( result, "000000" ) >= 0,
          strcmp( "ffffff", result ) >= 0 );
}




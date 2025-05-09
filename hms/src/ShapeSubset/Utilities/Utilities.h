
#ifndef UTILITIES_H
#define UTILITIES_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: Utilities.h - Declare utility routines.
HISTORY: 2009/05/26 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h> /* For FILE. */

#ifdef _WIN32
extern double drand48( void );
extern char* strtok_r( char*, const char*, char** );
#endif

/*================================= TYPES ==================================*/

enum { MINIMUM, MAXIMUM };
enum { LONGITUDE, LATITUDE };

enum { BYTE_TYPE, UINT16_TYPE, FLOAT_TYPE, GRID_DATA_TYPES };

#define IS_VALID_GRID_DATA_TYPE(type) ((type) >= 0 && (type) < GRID_DATA_TYPES)

enum { HOURLY, DAILY, MONTHLY, YEARLY }; /* Timestep of data. */

#define IS_VALID_TIMESTEP_TYPE(type) IN5((type), HOURLY, DAILY, MONTHLY,YEARLY)

typedef double Bounds[ 2 ][ 2 ]; /* [ LONGITUDE LATITUDE ][MINIMUM MAXIMUM].*/

typedef struct { float r, g, b; } Color;

typedef Color (*DataColor)( double, double, double ); /* Color function. */

#ifdef RGB
#undef RGB
#endif

typedef struct { unsigned char r, g, b; const char* s; } RGB;

typedef RGB RGBColormap[ 256 ];

typedef char StateCode[ 2 + 1 ]; /* E.g., "ME". */
typedef char LongName[ 63 + 1 ];
typedef char Name[  31 + 1 ];
typedef char Units[ 15 + 1 ];
typedef char Line[ 79 + 1 ];

typedef struct {
  long long id;     /* Station ID, e.g., 01010000, given or else hash of sid.*/
  int yyyymmdd;     /* Datestamp, e.g., 20020701. */
  int hhmmss;       /* timestamp, e.g., -1 if daily. */
  int aggregate;    /* Aggregate points to hour? */
  double longitude;  /* E.g., -69.71542500. */
  double latitude;   /* E.g., 46.70045278. */
  double measure;    /* E.g., 1.23. */
  double measure2;   /* -9999.0 if unused else v-component of wind/current. */
  double z;          /* Elevation or, if negative, depth in meters. */
  LongName source;  /* E.g., "Water_Quality_NWIS". */
  Name sid;         /* String id, e.g., "MDEDAT07-0103031". */
  Name name;        /* E.g., "Discharge". */
  Units units;      /* E.g., "m3/s". */
  Line metadata;    /* E.g., "PUGE;WA". */
  Color color;      /* E.g., color = dataColor( data, 0.0, 10.0 ). */
} PointData;

typedef void (*Project)( double, double, double*, double* );
typedef void (*Unproject)( double, double, double*, double* );
typedef void (*ForEachFileCallback)( const char* const fileName, void* data );

/*=============================== FUNCTIONS ================================*/

/* Geometric functions: */

extern int isValidLongitudeLatitude( double longitude, double latitude );

extern int isValidBounds( const Bounds bounds );

extern int overlap( const Bounds a, const Bounds b );

extern int clampedRangesOverlap( const int minimum, const int maximum,
                                 int* lower, int* upper );

extern double degreesPerPixel( const Bounds bounds, int width, int height );

extern int uniquePoints( double longitude1, double latitude1,
                         double longitude2, double latitude2,
                         double tolerance );

extern int clipLine( double wxl, double wyl, double wxu, double wyu,
                     double* x1, double* y1, double* x2, double* y2 );

extern int clipCoordinate( double p, double q, double* t1, double* t2 );

extern int pointInsideTriangle( double x, double y,
                                double x1, double y1,
                                double x2, double y2,
                                double x3, double y3 );

extern double areaOfTriangle( double x1, double y1,
                              double x2, double y2,
                              double x3, double y3 );

extern double pointLineDistance( double x, double y,
                                 double x1, double y1,
                                 double x2, double y2 );

extern int colinear( const double x1, const double y1,
                     const double x2, const double y2,
                     const double x3, const double y3 );

extern int readMapFileHeader( FILE* file,
                              int* polylineCount, int* vertexCount );

extern int readMapFileData( FILE* file,
                            int polylineCount, int vertexCount,
                            int counts[], float vertices[] );

extern void subsetMap( int inputPolylineCount, int inputVertexCount,
                       const int inputCounts[], const float inputVertices[],
                       double resolution, const double bounds[ 2 ][ 2 ],
                       int* outputPolylineCount, int* outputVertexCount,
                       int outputCounts[], float outputVertices[] );

extern void subsetMapDouble( int inputPolylineCount, int inputVertexCount,
                             const int inputCounts[],
                             const double inputVertices[],
                             double resolution, const double bounds[ 2 ][ 2 ],
                             int* outputPolylineCount, int* outputVertexCount,
                             int outputCounts[], double outputVertices[] );

/* Date-time functions: */

extern int daysInMonth( int year, int month );
extern int isValidYYYYMMDD( int yyyymmdd );
extern int incrementDate( int yyyymmdd, int days );
extern int isValidYYYYMMDDHH( int yyyymmddhh );
extern int incrementDateTime( int yyyymmddhh, int hours );
extern int endOfDay( int yyyymmddhh );
extern int endOfMonth( int yyyymmddhh );
extern int endOfYear( int yyyymmddhh );
extern int isValidColor( Color color );
extern double pointValue( const PointData* const pointData );

extern void pointDataRange( const int count, const PointData pointData[],
                            int* indexOfMinimum, int* indexOfMaximum );

extern int pointMatches( const PointData* const pointData,
                         const int yyyymmddhh, const int timestepType );

extern RGB makeRGB( unsigned char r, unsigned char g, unsigned char b,
                    const char* s );

extern Color dataColor( double value, double minimum, double maximum );
extern Color grayColor( double value, double minimum, double maximum );
extern Color greenColor( double value, double minimum, double maximum );
extern Color greenGrayColor( double value, double minimum, double maximum );
extern Color darkGreenGrayColor( double value, double minimum, double maximum );
extern Color cyanGrayColor( double value, double minimum, double maximum );
extern Color blueGrayColor( double value, double minimum, double maximum );
extern Color brownGrayColor( double value, double minimum, double maximum );
extern Color tanGrayColor( double value, double minimum, double maximum );
extern Color modulo6Color( double value, double minimum, double maximum );
extern Color elevationColor( double value, double minimum, double maximum );
extern Color seagrassColor( double value, double minimum, double maximum );
extern Color riskColor( double value, double minimum, double maximum );
extern Color riskColor4( double value, double minimum, double maximum );
extern Color riskColorText( const char* text );
extern Color soilColorText( const char* text );
extern Color soilColorDarkLight( double value, double minimum, double maximum );
extern Color soilColorLightDark7( double value, double minimum, double maximum );
extern Color soilColorLightDark( double value, double minimum, double maximum );
extern Color soilColorDarkLight4( double value, double minimum, double maximum );
extern Color precipitationColor( double value, double minimum, double maximum );
extern Color populationDensityColor(double value, double minimum, double maximum);
extern Color populationColor( double value, double minimum, double maximum );
extern Color lexicographicTextColor( const char* text );
extern Color sedimentColorText( const char* text );
extern Color sedimentColorMud( double value, double minimum, double maximum );
extern Color sedimentColorMudu( double value, double minimum, double maximum );
extern Color sedimentColorSand( double value, double minimum, double maximum );
extern Color sedimentColorGravel( double value, double minimum, double maximum );
extern Color sedimentColorRock( double value, double minimum, double maximum );
extern Color nitrogenColorText( const char* text );
extern Color wordColor( const char* word );
extern Color blackstoneHRUColor( double hru, double unused1, double unused2 );
extern Color charlesHRUColor( double hru, double unused1, double unused2 );
extern Color farmingtonHRUColor( double hru, double unused1, double unused2 );
extern Color ipswichHRUColor( double hru, double unused1, double unused2 );
extern Color pawcatuckHRUColor( double hru, double unused1, double unused2 );
extern Color sudburyHRUColor( double hru, double unused1, double unused2 );
extern Color tauntonHRUColor( double hru, double unused1, double unused2 );

extern Color iclusNaturalWaterColor(double percent,double unused1,double maximum);
extern Color iclusReservoirsCanalsColor( double percent,
                                         double unused1, double maximum );
extern Color iclusWetlandsColor(double percent,double unused1, double maximum);
extern Color iclusRecreationalConservationColor( double percent,
                                                 double unused1, double maximum);
extern Color iclusTimberColor(double percent,double unused1, double maximum);
extern Color iclusGrazingColor(double percent,double unused1, double maximum);
extern Color iclusPastureColor(double percent,double unused1, double maximum);
extern Color iclusCroplandColor(double percent,double unused1, double maximum);
extern Color iclusMiningBarrenColor(double percent,double unused1,double maximum);
extern Color iclusParksGolfCoursesColor( double percent,
                                         double unused1, double maximum );
extern Color iclusExurbanLowDensityColor( double percent,
                                          double unused1, double maximum);
extern Color iclusExurbanHighDensityColor( double percent,
                                           double unused1, double maximum );
extern Color iclusSuburbanColor(double percent,double unused1, double maximum);
extern Color iclusUrbanLowDensityColor( double percent,
                                        double unused1, double maximum );
extern Color iclusUrbanHighDensityColor( double percent,
                                         double unused1, double maximum );
extern Color iclusCommercialColor(double percent,double unused1, double maximum);
extern Color iclusIndustrialColor(double percent,double unused1, double maximum);
extern Color iclusInstitutionalColor( double percent,
                                      double unused1,double maximum);
extern Color iclusTransportationColor( double percent,
                                       double unused1, double maximum );

extern Color nlcdOpenWaterColor( double percent,
                                 double unused1, double maximum );
extern Color nlcdSnowIceColor( double percent, double unused1, double maximum );
extern Color nlcdWDevelopedOpenSpaceColor( double percent,
                                           double unused1, double maximum );
extern Color nlcdDevelopedLowIntensityColor( double percent,
                                             double unused1, double maximum );
extern Color nlcdDevelopedMediumIntensityColor( double percent,
                                                double unused1, double maximum );
extern Color nlcdDevelopedHighIntensityColor( double percent,
                                              double unused1, double maximum );
extern Color nlcdBarrenColor( double percent, double unused1, double maximum );
extern Color nlcdDeciduousForestColor( double percent,
                                       double unused1, double maximum );
extern Color nlcdEvergreenForestColor( double percent,
                                       double unused1, double maximum );
extern Color nlcdMixedForestColor(double percent, double unused1, double maximum);
extern Color nlcdShrubScrubColor( double percent, double unused1, double maximum);
extern Color nlcdGrasslandColor( double percent, double unused1, double maximum );
extern Color nlcdPastureHayColor( double percent, double unused1, double maximum);
extern Color nlcdCropsColor( double percent, double unused1, double maximum );
extern Color nlcdWetlandForestColor( double percent,
                                     double unused1, double maximum );
extern Color nlcdWetlandEmergentColor( double percent,
                                       double unused1, double maximum );

extern Color stratTypeTextColor( const char* text );
extern Color stratMethTextColor( const char* text );
extern Color ecoRegionTextColor( const char* text );

extern Color nhdCodeTextColor( const char* text );
extern Color nhdPlusIdTextColor( const char* text );
extern Color temperatureRegimeTextColor( const char* text );
extern int temperatureFlagLineStipple( const char* text );

extern Color streamTemperatureColor( double value,
                                     double unused_1, double unused_2 );

extern Color streamTemperatureCategoryColor( double value,
                                             double unused_1, double unused_2);

extern int hruStippling( const char* const name, const int hruId );

extern void daltonize( const int dichromacy, const size_t length, float rgb[] );

extern int indexOfString( const char* string, const char* const strings[],
                          int count );

extern int matchesWord( const char* word, const char* const words );

extern int matchesPattern(const char* const string, const char* const pattern);

extern void accumulate( const float* input, float* output, int count );

extern double scaledMaximum( float array[], int count, double scale,
                            double threshold );

extern int unsignedShortMaximum( const unsigned short array[],
                                 const size_t count );

extern int charMaximum( const char array[], const size_t count );

extern void swapCharDataRows( char array[],
                              const int timesteps,
                              const int rows, const int columns );

extern void swapUnsignedShortDataRows( unsigned short array[],
                                       const int timesteps,
                                       const int rows, const int columns );

extern void swapFloatDataRows( float array[],
                               const int timesteps,
                               const int rows, const int columns );

extern void expandBytesToFloats( float array[], int count );

extern double widthInMeters( const Bounds bounds );

extern void computeArrowVectorCoordinates( const double longitude,
                                           const double latitude,
                                           const double x,
                                           const double y,
                                           const double degreesPerPixel,
                                           const double pixelsPerUnitLength,
                                           double* const point0X,
                                           double* const point0Y,
                                           double* const point1X,
                                           double* const point1Y,
                                           double* const point2X,
                                           double* const point2Y,
                                           double* const point3X,
                                           double* const point3Y );

extern int wordCount( const char* string );
extern int lineCount( const char* string );
extern void lowercase( char string[] );
extern void uppercase( char string[] );
extern void changeChar( char string[], const char from, const char to );
extern size_t countChar( const char string[], const char ch );
extern void eraseChar( char string[], const char ch );
extern void shortenName( const LongName name, LongName shortenedName );
extern void eraseTrailingWhitespace( char string[] );
extern void eraseLeadingWhitespace( char string[] );
extern void eraseDoubleQuotedCommas( char string[] );

extern void substituteWord( const char* input,
                            const char* oldWord,
                            const char* newWord,
                            char* output );

  

extern int parseInts( char** string,
                      const char* const delimiters,
                      const size_t count,
                      const int range[][ 2 ],
                      const int clamp,
                      int values[] );

extern int parseDoubles( char** string,
                         const char* const delimiters,
                         const size_t count,
                         const double range[][ 2 ],
                         const int clamp,
                         double values[] );

extern int parseWords( char** string,
                       const char* const delimiters,
                       const size_t valueCount,
                       const size_t wordCount,
                       const char* const words[],
                       int values[] );

extern int fileExists( const char* fileName );
extern size_t fileSize( const char* fileName );
extern char* readFile( const char* name, size_t* length, size_t* lines );
extern size_t controlMToSpace( char* string );
extern size_t copyFileLine( FILE* input, FILE* output );
extern size_t skipFileLine( FILE* file );
extern int copyFile( const char* inputFileName, const char* outputFileName );

extern int copyFileBytes( FILE* inputFile, const char* outputFileName,
                          size_t bytes );

extern int streamBytesToFile( FILE* stream, const char* outputFileName );

extern int streamFile( const char* fileName );

extern int directoryExists( const char* name );

extern void removeAllFiles( const char* name );

extern void removeMatchedFiles( const char* const directoryName,
                                const char* const startsWith,
                                const char* const endsWith );

extern void copyFiles( const char* const fromDirectory,
                       const char* const toDirectory,
                       const char* const startsWith,
                       const char* const endsWith );

extern void forEachFile( const char* const directory,
                         const char* const startsWith,
                         const char* const endsWith,
                         ForEachFileCallback callback,
                         void* callbackData );

extern void directoryListing( const char* directory, const char* extensions,
                              int size, char buffer[] );

extern const char* homeDirectory( void );

extern size_t sortUniqFile( const char* name, const size_t headerLines );

#ifdef __cplusplus
}
#endif

#endif /* UTILITIES_H */


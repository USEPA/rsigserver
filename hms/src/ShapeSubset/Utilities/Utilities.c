
/******************************************************************************
PURPOSE: Utilities.c - Utility routines.
HISTORY: 2009/05/26 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>     /* For FILE, fopen(), fread(), snprintf(). fprintf(). */
#include <string.h>    /* For strcmp(), strncpy(), strncat(), strtok_r(). */
#include <stdlib.h>    /* For malloc(), free(), strtol(), strtod(). */
#include <ctype.h>     /* For isspace(), isprint(). */
#include <math.h>      /* For sqrt(). */
#include <unistd.h>    /* For unlink(). */
#include <time.h>      /* For struct tm, struct timespec, localtime(). */
#include <dirent.h>    /* For DIR, struct dirent, opendir(), readdir(), etc.*/
#include <sys/types.h> /* For struct stat. */
#include <sys/stat.h>  /* For stat(). */


#ifdef _WIN32

/* HACKS to fix deficiencies in stupid Windoze! */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h> /* For rand(). */

double drand48( void ) { return (double) rand() / RAND_MAX; }

char* strtok_r( char* s, const char* delim, char** lasts ) {

  char* spanp = 0;
  int c = 0, sc = 0;
  char* tok = 0;

  if ( s == 0 && ( s = *lasts ) == 0 ) {
    return 0;
  }

cont:
  c = *s++;

  for ( spanp = (char*) delim; ( sc = *spanp++ ) != 0; ) {

    if ( c == sc ) {
      goto cont;
    }
  }

  if ( c == 0 ) {
    *lasts = 0;
    return 0;
  }

  tok = s - 1;

  while ( 1 ) {
    c = *s++;
    spanp = (char*) delim;

    do {

      if ( ( sc = *spanp++ ) == c ) {

        if ( c == 0 ) {
          s = 0;
        } else {
          s[ -1 ] = 0;
        }

        *lasts = s;
        return tok;
      }

    } while ( sc != 0 );
  }

  return 0;
}


#ifdef __cplusplus
} /* extern "C". */
#endif


#endif /* Winhell. */



#include <Assertions.h>    /* For PRE*(), POST*(), CHECK*(), AND*(). */
#include <BasicNumerics.h> /* For rotate4ByteArrayIfLittleEndian(), isNan().*/
#include <Failure.h>       /* For failureMessage(). */
#include <Utilities.h>     /* For public interface. */

#ifdef MIN
#undef MIN
#endif
#define MIN( a, b ) ((a) < (b) ? (a) : (b))

/* Minimum rgb value for ICLUS and NLCD category data: */

#define CATEGORY_MINIMUM ( 32.0 / 255.0 )

/* If USE_SCALE_FACTOR is 1 use non-linear scaling to boost category hue. */

#define USE_CATEGORY_SCALE_FACTOR 1

#if USE_CATEGORY_SCALE_FACTOR
/* #define CATEGORY_SCALE_FACTOR( x ) ( exp( x ) * ( 1.0 / M_E ) ) */
/* #define CATEGORY_SCALE_FACTOR( x ) sqrt( x ) */
#define CATEGORY_SCALE_FACTOR( x ) ( ( (x) + sqrt( x ) * 0.5 ) * ( 2.0 / 3.0 ) )
#else
#define CATEGORY_SCALE_FACTOR( x ) ( x )
#endif

/*=========================== FORWARD DECLARATIONS ==========================*/

static void removeFile0( const char* const name, void* unused ) {
  PRE02( name, *name );
  unlink( name );
}

static int stringCompare( const void* va, const void* vb ) {
  const char** a = (const char**) va;
  const char** b = (const char**) vb;
  return strcmp( *a, *b );
}

static int minimumItemi( const int array[], int count );

static int sum( const int array[], int count );

static int inBounds( const float vertices[], int count,
                     const double bounds[ 2 ][ 2 ] );

static int inBoundsDouble( const double vertices[], int count,
                           const double bounds[ 2 ][ 2 ] );

static Color soilColor( double t );

static Color soilColor4( double t );

static Color soilColor7( double t );

static Color categoryColor( const double percent, const double maximum,
                            const int red, const int green, const int blue );

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: isValidLongitudeLatitude - Is longitude, latitude a valid point?
INPUTS:  double longitude Longitude to check.
         double latitude  Latitude to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidLongitudeLatitude( double longitude, double latitude ) {
  const int result =
    AND2( IN_RANGE( longitude, -180.0, 180.0 ),
          IN_RANGE( latitude,   -90.0,  90.0 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidBounds - Is bounds a valid lon-lat rectangle?
INPUTS:  const Bounds bounds  Bounds to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidBounds( const Bounds bounds ) {
  const int result =
    AND5( bounds,
          IN_RANGE( bounds[ LONGITUDE ][ MINIMUM ], -180.0, 180.0 ),
          IN_RANGE( bounds[ LONGITUDE ][ MAXIMUM ],
                    bounds[ LONGITUDE ][ MINIMUM ], 180.0 ),
          IN_RANGE( bounds[ LATITUDE ][ MINIMUM ], -90.0, 90.0 ),
          IN_RANGE( bounds[ LATITUDE ][ MAXIMUM ],
                    bounds[ LATITUDE ][ MINIMUM ], 90.0 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: overlap - Do rectangles overlap?
INPUTS:  const Bounds a  a[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
         const Bounds b  b[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
RETURNS: int 1 if overlap, else 0.
******************************************************************************/

int overlap( const Bounds a, const Bounds b ) {
  PRE02( isValidBounds( a ), isValidBounds( b ) );
  const int outside =
    OR4( a[ LATITUDE  ][ MINIMUM ] > b[ LATITUDE  ][ MAXIMUM ],
         a[ LATITUDE  ][ MAXIMUM ] < b[ LATITUDE  ][ MINIMUM ],
         a[ LONGITUDE ][ MINIMUM ] > b[ LONGITUDE ][ MAXIMUM ],
         a[ LONGITUDE ][ MAXIMUM ] < b[ LONGITUDE ][ MINIMUM ] );
  const int result = ! outside;
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: clampedRangesOverlap - If the range [lower, upper] overlaps
         the range [minimum, maximum] return 1 and clamp [lower, upper] to
         within [minimum, maximum],
         else return 0 and [lower, upper] are unchanged.
INPUTS:  const int minimum   Minimum of clamp-to range.
         const int maximum   Maximum of clamp-to range.
         int* lower          Minimum of range to clamp.
         int* upper          Maximum of range to clamp.
OUTPUTS: int* lower          If overlap, clamped minimum of range to clamp.
         int* upper          If overlap, clamped maximum of range to clamp.
RETURNS: int 1 if overlap, else 0.
******************************************************************************/

int clampedRangesOverlap( const int minimum, const int maximum,
                          int* lower, int* upper ) {
  PRE04( minimum <= maximum, lower, upper, *lower <= *upper );
  const int outside = OR2( *upper < minimum, *lower > maximum );
  const int result = ! outside;

  if ( result ) {

    if ( *lower < minimum ) {
      *lower = minimum;
    }

    if ( *upper < *lower ) {
      *upper = *lower;
    } else if ( *upper > maximum ) {
      *upper = maximum;
    }
  }

  POST02( IS_BOOL( result ),
                   IMPLIES( result,
                            AND2( IN_RANGE( *lower, minimum, maximum ),
                                  IN_RANGE( *upper, *lower, maximum ) ) ) );
  return result;
}



/******************************************************************************
PURPOSE: degreesPerPixel - Degrees lon-lat per screen pixel.
INPUTS:  const Bounds bounds  Bounds of lon-lat window.
         int width            Width of window in screen pixels.
         int height           Height of window in screen pixels.
RETURNS double degrees of lon-lat per screen pixel.
******************************************************************************/

double degreesPerPixel( const Bounds bounds, int width, int height ) {
  PRE03( isValidBounds( bounds ), width > 0, height > 0);
  const double longitudeRange =
    bounds[ LONGITUDE ][ MAXIMUM ] - bounds[ LONGITUDE ][ MINIMUM ];
  const double latitudeRange =
    bounds[ LATITUDE ][ MAXIMUM ] - bounds[ LATITUDE ][ MINIMUM ];
  const double longitudeDegreesPerPixel = longitudeRange / width;
  const double latitudeDegreesPerPixel  = latitudeRange / height;
  double result =
    longitudeDegreesPerPixel > latitudeDegreesPerPixel ?
  longitudeDegreesPerPixel : latitudeDegreesPerPixel;

  if ( result <= 0.0 ) {
    result = 1e-4;;
  }

  DEBUG( fprintf( stderr, "degreesPerPixel = %lf\n", result ); )
  POST0( result > 0.0 );
/*POST0( IN_RANGE( result, MINIMUM_DEGREES_PER_PIXEL, 360.0 ) ); */
  return result;
}



/******************************************************************************
PURPOSE: pointInsideTriangle - Determine if point (x, y) is inside triangle
         with vertices (x1, y1), (x2, y2), (x3, y3).
INPUTS:  double x    X-Coordinate of point to test.
         double y    Y-Coordinate of point to test.
         double x1   X-Coordinate of 1st vertex of triangle.
         double y1   Y-Coordinate of 1st vertex of triangle.
         double x2   X-Coordinate of 2nd vertex of triangle.
         double y2   Y-Coordinate of 2nd vertex of triangle.
         double x3   X-Coordinate of 3rd vertex of triangle.
         double y3   Y-Coordinate of 3rd vertex of triangle.
RETURNS: int 1 if inside, else 0.
******************************************************************************/

int pointInsideTriangle( double x, double y,
                         double x1, double y1,
                         double x2, double y2,
                         double x3, double y3 ) {
  const double scale = 1.01; /* Allow within 1% larger for round-off error. */
  const double triangleArea = scale * areaOfTriangle( x1, y1, x2, y2, x3, y3 );
  double area = areaOfTriangle( x, y, x2, y2, x3, y3 );
  int result = area <= triangleArea;

  if ( result ) {
    area += areaOfTriangle( x, y, x1, y1, x2, y2 );
    result = area <= triangleArea;

    if ( result ) {
      area += areaOfTriangle( x, y, x1, y1, x3, y3 );
      result = area <= triangleArea;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: areaOfTriangle - Area of triangle with vertices
         (x1, y1), (x2, y2), (x3, y3).
INPUTS:  double x1   X-Coordinate of 1st vertex of triangle.
         double y1   Y-Coordinate of 1st vertex of triangle.
         double x2   X-Coordinate of 2nd vertex of triangle.
         double y2   Y-Coordinate of 2nd vertex of triangle.
         double x3   X-Coordinate of 3rd vertex of triangle.
         double y3   Y-Coordinate of 3rd vertex of triangle.
RETURNS: double area.
******************************************************************************/

double areaOfTriangle( double x1, double y1,
                       double x2, double y2,
                       double x3, double y3 ) {
  const double a = x1 - x3;
  const double b = y1 - y3;
  const double c = x2 - x3;
  const double d = y2 - y3;
  const double triangleArea = 0.5 * ( a * d - b * c );
  const double result = triangleArea < 0.0 ? -triangleArea : triangleArea;
  return result;
}



/******************************************************************************
PURPOSE: pointLineDistance - Nearest distance from point (x, y) to infinite
         undirected line containing points (x1, y1)--(x2, y2).
INPUTS:  const double x   X-Coordinate of point to check.
         const double y   Y-Coordinate of point to check.
         const double x1  X-Coordinate of 1st endpoint of line.
         const double y1  Y-Coordinate of 1st endpoint of line.
         const double x2  X-Coordinate of 2nd endpoint of line.
         const double y2  Y-Coordinate of 2nd endpoint of line.
RETURNS: double distance  Nearest distance from point to line.
NOTES:   Swokowski, Earl, "Calculus With Analytic Geometry, Second Edition",
         Prindle, Weber & Schmidt, Boston, MA, 1979, page 696.
******************************************************************************/

double pointLineDistance( const double x, const double y,
                          const double x1, const double y1,
                          const double x2, const double y2 ) {

  PRE06( ! isNan( x ) , ! isNan( y ),
         ! isNan( x1 ) , ! isNan( y1 ),
         ! isNan( x2 ) , ! isNan( y2 ) );
  double result = 0.0;
  const double dx = x2 - x1;
  const double dy = y2 - y1;
  const double lineLength = sqrt( dx * dx + dy * dy );

  if ( lineLength == 0.0 ) { /* Degenerate line so compute point distance: */
    const double dx0 = x - x1;
    const double dy0 = y - y1;
    result = sqrt( dx0 * dx0 + dy0 * dy0 );
  } else {
    const double recipricolLineLength = 1.0 / lineLength;

    if ( recipricolLineLength <= 1e-12 ) {

      if ( ! colinear( x, y, x1, y1, x2, y2 ) ) {
        result = DBL_MAX;
      }

    } else {
      const double px = x - x1;
      const double py = y - y1;
      const double z0 = dx * py - px * dy;
      const double z = z0 < 0.0 ? -z0 : z0;
      result = recipricolLineLength * z;
    }
  }

  POST0( result >= 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: colinear - Do the given three points lie along a line (or are they
         coincident)?
INPUTS:  const double x1  X-Coordinate of 1st point to check.
         const double y1  Y-Coordinate of 1st point to check.
         const double x2  X-Coordinate of 2nd point to check.
         const double y2  Y-Coordinate of 2nd point to check.
         const double x3  X-Coordinate of 3rd point to check.
         const double y3  Y-Coordinate of 3rd point to check.
RETURNS: int 1 if colinear, else 0.
NOTES:   Swokowski, Earl, "Calculus With Analytic Geometry, Second Edition",
         Prindle, Weber & Schmidt, Boston, MA, 1979, page 670.
******************************************************************************/

int colinear( const double x1, const double y1,
              const double x2, const double y2,
              const double x3, const double y3 ) {

  PRE06( ! isNan( x1 ) , ! isNan( y1 ),
         ! isNan( x2 ) , ! isNan( y2 ),
         ! isNan( x3 ) , ! isNan( y3 ) );
  /*
   * First check for coincidence:
   *   if any two (or all three) points are coincident then consider this
   *   degenerate case colinear
   * Else check for reflectance:
   *   if one of the points is the origin and the other two are negatives
   *   of each other then the three points are colinear.
   * Otherwise compute the result as follows:
   * 1. Form three vectors: P1->P2 and P1->P3 and P2->P3.
   * 2. Normalize the vectors and compute the
   *    dot-products: P1->P2 * P1->P3 and P1->P2 * P2->P3
   * 3. If both dot-products are approximately +/-1 then
   *    the points are colinear since the angle between the vectors is
   *    arccosine (+/-1) = 0 (or pi).
   * Swokowski, Earl, "Calculus With Analytic Geometry, Second Edition",
   * Prindle, Weber & Schmidt, Boston, MA, 1979, page 670.
   * Note: the coincidence and reflectance checks are necessary because
   * the angle-test is not numerically robust enough to handle these cases
   * correctly for the full range of non-NaN values.
   */

  const double tolerance = 1e-6;
  int result = 1;

  if ( ! AND2( x1 == x2, y1 == y2 )  ) { /* non-coincident. */

    if ( ! AND2( x1 == x3, y1 == y3 )  ) { /* non-coincident. */

      if ( ! AND2( x2 == x3, y2 == y3 )  ) { /* non-coincident. */

        if ( ! AND4( x1 == 0.0, y1 == 0.0,
                     x2 == -x3, y2 == -y3 ) ) { /* non-reflected */

          if ( ! AND4( x2 == 0.0, y2 == 0.0,
                       x1 == -x3, y1 == -y3 ) ) { /* non-reflected */

            if ( ! AND4( x3 == 0.0, y3 == 0.0,
                         x1 == -x2, y1 == -y2 ) ) { /* non-reflected */
              const double v1x_ = x2 - x1;
              const double v1y_ = y2 - y1;
              const double v1m = 1.0 / sqrt( v1x_ * v1x_ + v1y_ * v1y_ );
              const double v1x = v1x_ * v1m;
              const double v1y = v1y_ * v1m;
              const double v2x_ = x3 - x1;
              const double v2y_ = y3 - y1;
              const double v2m = 1.0 / sqrt( v2x_ * v2x_ + v2y_ * v2y_ );
              const double v2x = v2x_ * v2m;
              const double v2y = v2y_ * v2m;
              const double v1DotV2    = v1x * v2x + v1y * v2y;
              const double magnitude1 = fabs( v1DotV2 );
              const double oneMinus   = 1.0 - tolerance;
              const double onePlus    = 1.0 + tolerance;
              result = IN_RANGE( magnitude1, oneMinus, onePlus );

              if ( result ) {
                const double v3x_ = x3 - x2;
                const double v3y_ = y3 - y2;
                const double v3m = 1.0 / sqrt( v3x_ * v3x_ + v3y_ * v3y_ );
                const double v3x = v3x_ * v3m;
                const double v3y = v3y_ * v3m;
                const double v1DotV3    = v1x * v3x + v1y * v3y;
                const double magnitude2 = fabs( v1DotV3 );
                result = IN_RANGE( magnitude2, oneMinus, onePlus );
              }
            }
          }
        }
      }
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: clipLine - Check and clip a line segment to a rectangle window.
INPUTS:  double wxl  Window x-lower coordinate.
         double wyl  Window y-lower coordinate.
         double wxu  Window x-upper coordinate.
         double wyu  Window y-upper coordinate.
         double* x1  (x1, y1)-(x2, y2) is the line segment to clip.
         double* y1
         double* x2
         double* y2
OUTPUTS: double* x1  (x1, y1)-(x2, y2) is the clipped line segment.
         double* y1
         double* x2
         double* y2
RETURNS: int 1 if line segment intersects the rectangle.
NOTES:   Uses the Liang-Barsky 2D Line Clipping Algorithm (fastest known).
         ACM Transactions on Graphics Volume 3, Issue 1, Jan. 1984 pp 1–22
         https://doi.org/10.1145/357332.357333
******************************************************************************/

int clipLine( double wxl, double wyl, double wxu, double wyu,
              double* x1, double* y1, double* x2, double* y2 ) {

  int result = 0;
  const double dx = *x2 - *x1;
  double t1 = 0.0; /* t holds new start point. */
  double t2 = 1.0;

  /* Check boundaries: left, right, bottom, top: */

  if ( clipCoordinate( -dx, *x1 - wxl, &t1, &t2 ) ) { /* left. */

    if ( clipCoordinate( dx, wxu - *x1, &t1, &t2 ) ) { /* right. */
      const double dy = *y2 - *y1;

      if ( clipCoordinate( -dy, *y1 - wyl, &t1, &t2 ) ) { /* bottom. */

        if ( clipCoordinate( dy, wyu - *y1, &t1, &t2 ) ) { /* top. */

          /*
           * At least some of the line is within the window so
           * calculate the new end and start points (in that order).
           */

          if ( t2 < 1.0 ) { /* Calculate new end point first. */
            *x2 = *x1 + t2 * dx;
            *y2 = *y1 + t2 * dy;
          }

          if ( t1 > 0.0 ) { /* Calculate new start point. */
            *x1 += t1 * dx;
            *y1 += t1 * dy;
          }

          result = 1; /* Successfully clipped. */
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: clipCoordinate - Check and clip a line segment to a boundary.
INPUTS:  double p
         double q
         double* t1
         double* t2
OUTPUTS: double* t1
         double* t2
RETURNS: int 1 if point crosses the boundary.
NOTES:   Part of the Liang-Barsky 2D Line Clipping Algorithm.
         ACM Transactions on Graphics Volume 3, Issue 1, Jan. 1984 pp 1–22
         https://doi.org/10.1145/357332.357333
******************************************************************************/

int clipCoordinate( double p, double q, double* t1, double* t2 ) {
  int result = 1; /* Clipped. */

  if ( p < 0.0 ) { /* Line from outside to inside of that boundary. */
    const double r = q / p; /* Intersection coordinate. */

    if ( r > *t2 ) {
      result = 0; /* Intersection past segment end point. */
    } else if ( r > *t1 ) {
      *t1 = r; /* Intersection is past start point. */
    }
  } else if ( p > 0.0 ) { /* Line from inside to outside of that boundary. */
    const double r = q / p;

    if ( r < *t1 ) {
      result = 0; /* Intersection is before start point. */
    } else if ( r < *t2 ) {
      *t2 = r; /* Intersection is before end point. */
    }
  } else if ( q < 0.0 ) { /* p == 0.0. */
    result = 0; /* Line is parallel to that boundary. */
  }

  return result;
}



/******************************************************************************
PURPOSE: uniquePoints - Is |(longitude1,lattitude1)-(longitude2, latitude2)|max
         > tolerance?
INPUTS:  double longitude1  Longitude of 1st point to compare.
         double latitude1   Latitude  of 1st point to compare.
         double longitude2  Longitude of 2nd point to compare.
         double latitude2   Latitude  of 2nd point to compare.
         double tolerance   Minimum distance between points considered unique.
RETURNS: int 1 if both coordinates of points are more than tolerance apart.
******************************************************************************/

int uniquePoints( double longitude1, double latitude1,
                  double longitude2, double latitude2,
                  double tolerance ) {
  const double deltaLongitude = longitude1 - longitude2;
  const double deltaLatitude  = latitude1  - latitude2;
  const int result =
    deltaLongitude > tolerance || deltaLongitude < -tolerance ||
    deltaLatitude  > tolerance || deltaLatitude  < -tolerance;
  return result;
}



/******************************************************************************
PURPOSE: readMapFileHeader - Read header dimensions of map_*.bin file.
INPUTS:  FILE* file  Opened readable file.
OUTPUTS: int* polylineCount  Number of polylines in file.
         int* vertexCount    Number of lon-lat vertices in file.
RETURNS: int 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

int readMapFileHeader( FILE* file, int* polylineCount, int* vertexCount ) {
  PRE03( file, polylineCount, vertexCount );
  const int result =
    fscanf( file, "%*[^\n]\n%*[^\n]\n%d %d\n%*[^\n]%*c",
            polylineCount, vertexCount ) == 2 &&
    *polylineCount > 0 && *vertexCount > 0;

  if ( ! result ) {
    failureMessage( "Invalid map file header." );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: readMapFileData - Read data of map_*.bin file.
INPUTS:  FILE* file  Opened readable file, seeked to data.
         int polylineCount  Number of polylines in file.
         int vertexCount    Number of lon-lat vertex pairs in file.
OUTPUTS: int counts[ polylineCount ]  Array of vertex counts per polyline.
         float vertices[ vertexCount * 2 ] Array of lon-lat vertices.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int readMapFileData( FILE* file,
                     int polylineCount, int vertexCount,
                     int counts[], float vertices[] ) {
  PRE05( file, polylineCount > 0, vertexCount > 0, counts, vertices );
  assert_static( sizeof (int) == 4 );
  assert_static( sizeof (float) == 4 );
  const int vertexValues = vertexCount + vertexCount;
  int result = fread( counts, polylineCount * sizeof (int), 1, file );

  if ( result == 1 ) {
    result = fread( vertices, vertexValues * sizeof (float), 1, file );

    if ( result == 1 ) {
      rotate4ByteArrayIfLittleEndian( counts, polylineCount );
      rotate4ByteArrayIfLittleEndian( vertices, vertexValues );
    }
  }

  result = result == 1 ? 1 : 0;

  if ( ! result ) {
    failureMessage( "Invalid map file data." );
  }

#ifdef DEBUGGING
  {
    int index = 0;

    for ( index = 0; index < polylineCount; ++index ) {
      CHECK( counts[ index ] > 0 );
    }

    for ( index = 0; index < vertexCount; ++index ) {
      const int index2 = index + index;
      const float longitude = vertices[ index2 ];
      const float latitude = vertices[ index2 + 1 ];
      CHECK2( IN_RANGE( longitude, -180.0, 180.0 ),
              IN_RANGE( latitude,   -90.0,  90.0 ) );
    }
  }
#endif

  return result;
}



/******************************************************************************
PURPOSE: subsetMap - Count and clip polylines to bounds.
INPUTS:  int inputPolylineCount  Number of polylines in input data.
         int inputVertexCount    Number of vertices in input data.
         const int inputCounts[inputPolylines] Vertex count per input polyline.
         const float inputVertices[ inputVertexCount * 2 ] Input vertices.
OUTPUTS: double resolution  Minimum distance between included points.
         const double bounds[ 2 ][ 2 ]  Longitude-latitude bounds to include.
                           bounds[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
         int* outputPolylineCount  Number of polylines in output data.
         int* outputVertexCount    Number of vertices in output data.
         int outputCounts[outputPolylines] Vertex count per output polyline.
         float outputVertices[ outputVertexCount * 2 ] Output vertices if non-0
******************************************************************************/

void subsetMap( int inputPolylineCount, int inputVertexCount,
                const int inputCounts[], const float inputVertices[],
                double resolution, const double bounds[ 2 ][ 2 ],
                int* outputPolylineCount, int* outputVertexCount,
                int outputCounts[], float outputVertices[] ) {

  PRE012( inputPolylineCount > 0,
          inputVertexCount > 0,
          inputCounts,
          inputCounts[ 0 ] > 0,
          inputCounts[ inputPolylineCount - 1 ] > 0,
          inputVertices,
          IN_RANGE( inputVertices[ 0 ], -180.0, 180.0 ),
          IN_RANGE( inputVertices[ 1 ],  -90.0,  90.0 ),
          resolution >= 0.0,
          isValidBounds( (const double (*)[2]) bounds ),
          outputPolylineCount,
          outputVertexCount );

  const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const float* inputLongitudes  = inputVertices;
  const float* inputLatitudes   = inputVertices + 1;
  float* outputLongitudes = outputVertices;
  float* outputLatitudes  = outputVertices ? outputVertices + 1 : 0;
  int outputPolyline = 0;
  int polyline = 0;
  double lastStoredLongitude = -999.0; /* Sentinel value. */
  double lastStoredLatitude  = -999.0;

  *outputPolylineCount = *outputVertexCount = 0; /* Initialize outputs. */

  if ( outputCounts ) {
    outputCounts[ 0 ] = 0;
  }

  if ( outputVertices ) {
    outputVertices[ 0 ] = outputVertices[ 1 ] = 0.0;
  }

  DEBUG( fprintf( stderr, "\nsubset( %d polylines %d vertices ):\n",
                  inputPolylineCount, inputVertexCount ); )

  for ( polyline = 0; polyline < inputPolylineCount; ++polyline ) {
    const int thisPolylineVertexCount = inputCounts[ polyline ];
    int vertex = 0;
    double longitude1 = *inputLongitudes;
    double latitude1  = *inputLatitudes;
    inputLongitudes += 2;
    inputLatitudes  += 2;

    DEBUG2( fprintf( stderr, "polyline # %6d, vertex count = %6d\n",
                    polyline, thisPolylineVertexCount ); )

    for ( vertex = 1; vertex < thisPolylineVertexCount; ++vertex ) {
      const double longitude2 = *inputLongitudes;
      const double latitude2  = *inputLatitudes;
      double clippedLongitude1 = longitude1;
      double clippedLatitude1  = latitude1;
      double clippedLongitude2 = longitude2;
      double clippedLatitude2  = latitude2;
      inputLongitudes += 2;
      inputLatitudes  += 2;

      DEBUG2( fprintf( stderr, "before (%g, %g)---(%g, %g)\n",
                       clippedLongitude1, clippedLatitude1,
                       clippedLongitude2, clippedLatitude2 ); )

      if ( ( resolution == 0.0 ||
             uniquePoints( longitude1, latitude1, longitude2, latitude2,
                           resolution ) ) &&
           clipLine( longitudeMinimum, latitudeMinimum,
                     longitudeMaximum, latitudeMaximum,
                     &clippedLongitude1, &clippedLatitude1,
                     &clippedLongitude2, &clippedLatitude2 ) ) {

        const int discontiguous =
          clippedLongitude1 != lastStoredLongitude ||
          clippedLatitude1  != lastStoredLatitude;
        const int addedVertices = 1 + discontiguous;

        DEBUG2( fprintf( stderr, "after  (%g, %g)---(%g, %g)\n",
                         clippedLongitude1, clippedLatitude1,
                         clippedLongitude2, clippedLatitude2 ); )
        DEBUG2( fprintf( stderr, "discontiguous = %d\n", discontiguous );)
        DEBUG2( fprintf( stderr, "addedVertices = %d\n", addedVertices );)
        CHECK( IN_RANGE( clippedLongitude1,
                         bounds[ LONGITUDE ][ MINIMUM ],
                         bounds[ LONGITUDE ][ MAXIMUM ] ) );
        CHECK( IN_RANGE( clippedLongitude2,
                         bounds[ LONGITUDE ][ MINIMUM ],
                         bounds[ LONGITUDE ][ MAXIMUM ] ) );
        CHECK( IN_RANGE( clippedLatitude1,
                         bounds[ LATITUDE ][ MINIMUM ],
                         bounds[ LATITUDE ][ MAXIMUM ] ) );
        CHECK( IN_RANGE( clippedLatitude2,
                         bounds[ LATITUDE ][ MINIMUM ],
                         bounds[ LATITUDE ][ MAXIMUM ] ) );

        if ( outputVertices ) {

          if ( addedVertices == 2 ) {
            *outputLongitudes = clippedLongitude1;
            outputLongitudes += 2;
            *outputLatitudes  = clippedLatitude1;
            outputLatitudes  += 2;
            *outputLongitudes = clippedLongitude2;
            outputLongitudes += 2;
            *outputLatitudes  = clippedLatitude2;
            outputLatitudes  += 2;
          } else {
            *outputLongitudes = clippedLongitude2;
            outputLongitudes += 2;
            *outputLatitudes  = clippedLatitude2;
            outputLatitudes  += 2;
          }
        }

        *outputVertexCount += addedVertices;

        if ( discontiguous && lastStoredLongitude > -900.0 ) {
          ++outputPolyline;
          DEBUG2(fprintf(stderr,"outputPolyline index = %d\n",outputPolyline);)

          if ( outputCounts ) {
            outputCounts[ outputPolyline ] = 0;
          }
        }

        if ( outputCounts ) {
          outputCounts[ outputPolyline ] += addedVertices;
          DEBUG2( fprintf( stderr, "outputCounts[ %6d ] = %6d\n",
                           outputPolyline, outputCounts[ outputPolyline ] ); )
          CHECK( outputCounts[ outputPolyline ] >= 2 );
        }

        DEBUG2( fprintf( stderr, "total = %6d\n", *outputVertexCount ); )

        lastStoredLongitude = clippedLongitude2;
        lastStoredLatitude  = clippedLatitude2;

        DEBUG2( fprintf( stderr, "lastStored = (%g, %g)\n",
                         lastStoredLongitude, lastStoredLatitude ); )
      }

      longitude1 = longitude2;
      latitude1  = latitude2;

    } /* End loop on this polyline's vertices. */
  } /* End loop on polylines. */

  *outputPolylineCount = outputPolyline + ( *outputVertexCount != 0 );

  if ( *outputPolylineCount == 0 ) {
    *outputVertexCount = 0;

    if ( outputCounts ) {
      outputCounts[ 0 ] = 0;
    }

    if ( outputVertices ) {
      outputVertices[ 0 ] = outputVertices[ 1 ] = 0.0;
    }
  }

  DEBUG( fprintf( stderr, "*outputPolylineCount = %d\n",*outputPolylineCount);)
  DEBUG( fprintf( stderr, "*outputVertexCount = %d\n",*outputVertexCount);)
  CHECK( IMPLIES( AND2( *outputPolylineCount > 0, outputVertices ),
                  inBounds( outputVertices, *outputVertexCount,
                            (const double (*)[2]) bounds ) ) );

  POST0( IMPLIES_ELSE( *outputPolylineCount == 0, *outputVertexCount == 0,
           AND3( *outputVertexCount >= 2,
                 IMPLIES( outputCounts,
                          AND2( minimumItemi( outputCounts,
                                              *outputPolylineCount ) >= 2,
                                sum( outputCounts, *outputPolylineCount )
                                  == *outputVertexCount ) ),
                 IMPLIES( outputVertices,
                          inBounds( outputVertices, *outputVertexCount,
                                   (const double (*)[2]) bounds ) ) ) ) );
}



/******************************************************************************
PURPOSE: subsetMapDouble - Count and clip polylines to bounds.
INPUTS:  int inputPolylineCount  Number of polylines in input data.
         int inputVertexCount    Number of vertices in input data.
         const int inputCounts[inputPolylines] Vertex count per input polyline.
         const double inputVertices[ inputVertexCount * 2 ] Input vertices.
OUTPUTS: double resolution  Minimum distance between included points.
         const double bounds[ 2 ][ 2 ]  Longitude-latitude bounds to include.
                           bounds[ LONGITUDE, LATITUDE ][ MINIMUM, MAXIMUM ].
         int* outputPolylineCount  Number of polylines in output data.
         int* outputVertexCount    Number of vertices in output data.
         int outputCounts[outputPolylines] Vertex count per output polyline.
         double outputVertices[ outputVertexCount * 2] Output vertices if non-0
******************************************************************************/

void subsetMapDouble( int inputPolylineCount, int inputVertexCount,
                      const int inputCounts[], const double inputVertices[],
                      double resolution, const double bounds[ 2 ][ 2 ],
                      int* outputPolylineCount, int* outputVertexCount,
                      int outputCounts[], double outputVertices[] ) {

  PRE012( inputPolylineCount > 0,
          inputVertexCount > 0,
          inputCounts,
          inputCounts[ 0 ] > 0,
          inputCounts[ inputPolylineCount - 1 ] > 0,
          inputVertices,
          IN_RANGE( inputVertices[ 0 ], -180.0, 180.0 ),
          IN_RANGE( inputVertices[ 1 ],  -90.0,  90.0 ),
          resolution >= 0.0,
          isValidBounds( (const double (*)[2]) bounds ),
          outputPolylineCount,
          outputVertexCount );

  const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const double* inputLongitudes  = inputVertices;
  const double* inputLatitudes   = inputVertices + 1;
  double* outputLongitudes = outputVertices;
  double* outputLatitudes  = outputVertices ? outputVertices + 1 : 0;
  int outputPolyline = 0;
  int polyline = 0;
  double lastStoredLongitude = -999.0; /* Sentinel value. */
  double lastStoredLatitude  = -999.0;

  *outputPolylineCount = *outputVertexCount = 0; /* Initialize outputs. */

  if ( outputCounts ) {
    outputCounts[ 0 ] = 0;
  }

  DEBUG( fprintf( stderr, "\nsubset( %d polylines %d vertices ):\n",
                  inputPolylineCount, inputVertexCount ); )

  for ( polyline = 0; polyline < inputPolylineCount; ++polyline ) {
    const int thisPolylineVertexCount = inputCounts[ polyline ];
    int vertex = 0;
    double longitude1 = *inputLongitudes;
    double latitude1  = *inputLatitudes;
    inputLongitudes += 2;
    inputLatitudes  += 2;

    DEBUG2( fprintf( stderr, "polyline # %6d, vertex count = %6d\n",
                    polyline, thisPolylineVertexCount ); )

    for ( vertex = 1; vertex < thisPolylineVertexCount; ++vertex ) {
      const double longitude2 = *inputLongitudes;
      const double latitude2  = *inputLatitudes;
      double clippedLongitude1 = longitude1;
      double clippedLatitude1  = latitude1;
      double clippedLongitude2 = longitude2;
      double clippedLatitude2  = latitude2;
      inputLongitudes += 2;
      inputLatitudes  += 2;

      DEBUG2( fprintf( stderr, "before (%lg, %lg)---(%lg, %lg)\n",
                       clippedLongitude1, clippedLatitude1,
                       clippedLongitude2, clippedLatitude2 ); )

      if ( ( resolution == 0.0 ||
             uniquePoints( longitude1, latitude1, longitude2, latitude2,
                           resolution ) ) &&
           clipLine( longitudeMinimum, latitudeMinimum,
                     longitudeMaximum, latitudeMaximum,
                     &clippedLongitude1, &clippedLatitude1,
                     &clippedLongitude2, &clippedLatitude2 ) ) {

        const int discontiguous =
          clippedLongitude1 != lastStoredLongitude ||
          clippedLatitude1  != lastStoredLatitude;
        const int addedVertices = 1 + discontiguous;

        DEBUG2( fprintf( stderr, "after  (%lg, %lg)---(%lg, %lg)\n",
                         clippedLongitude1, clippedLatitude1,
                         clippedLongitude2, clippedLatitude2 ); )
        DEBUG2( fprintf( stderr, "discontiguous = %d\n", discontiguous );)
        DEBUG2( fprintf( stderr, "addedVertices = %d\n", addedVertices );)
        CHECK( IN_RANGE( clippedLongitude1,
                         bounds[ LONGITUDE ][ MINIMUM ],
                         bounds[ LONGITUDE ][ MAXIMUM ] ) );
        CHECK( IN_RANGE( clippedLongitude2,
                         bounds[ LONGITUDE ][ MINIMUM ],
                         bounds[ LONGITUDE ][ MAXIMUM ] ) );
        CHECK( IN_RANGE( clippedLatitude1,
                         bounds[ LATITUDE ][ MINIMUM ],
                         bounds[ LATITUDE ][ MAXIMUM ] ) );
        CHECK( IN_RANGE( clippedLatitude2,
                         bounds[ LATITUDE ][ MINIMUM ],
                         bounds[ LATITUDE ][ MAXIMUM ] ) );

        if ( outputVertices ) {

          if ( addedVertices == 2 ) {
            *outputLongitudes = clippedLongitude1;
            outputLongitudes += 2;
            *outputLatitudes  = clippedLatitude1;
            outputLatitudes  += 2;
            *outputLongitudes = clippedLongitude2;
            outputLongitudes += 2;
            *outputLatitudes  = clippedLatitude2;
            outputLatitudes  += 2;
          } else {
            *outputLongitudes = clippedLongitude2;
            outputLongitudes += 2;
            *outputLatitudes  = clippedLatitude2;
            outputLatitudes  += 2;
          }
        }

        *outputVertexCount += addedVertices;

        if ( discontiguous && lastStoredLongitude > -900.0 ) {
          ++outputPolyline;
          DEBUG2(fprintf(stderr,"outputPolyline index = %d\n",outputPolyline);)

          if ( outputCounts ) {
            outputCounts[ outputPolyline ] = 0;
          }
        }

        if ( outputCounts ) {
          outputCounts[ outputPolyline ] += addedVertices;
          DEBUG2( fprintf( stderr, "outputCounts[ %6d ] = %6d\n",
                           outputPolyline, outputCounts[ outputPolyline ] ); )
          CHECK( outputCounts[ outputPolyline ] >= 2 );
        }

        DEBUG2( fprintf( stderr, "total = %6d\n", *outputVertexCount ); )

        lastStoredLongitude = clippedLongitude2;
        lastStoredLatitude  = clippedLatitude2;

        DEBUG2( fprintf( stderr, "lastStored = (%lg, %lg)\n",
                         lastStoredLongitude, lastStoredLatitude ); )
      }

      longitude1 = longitude2;
      latitude1  = latitude2;

    } /* End loop on this polyline's vertices. */
  } /* End loop on polylines. */

  *outputPolylineCount = outputPolyline + ( *outputVertexCount != 0 );

  if ( *outputPolylineCount == 0 ) {
    *outputVertexCount = 0;

    if ( outputCounts ) {
      outputCounts[ 0 ] = 0;
    }

    if ( outputVertices ) {
      outputVertices[ 0 ] = outputVertices[ 1 ] = 0.0;
    }
  }

  DEBUG( fprintf( stderr, "*outputPolylineCount = %d\n",*outputPolylineCount);)
  DEBUG( fprintf( stderr, "*outputVertexCount = %d\n",*outputVertexCount);)
  CHECK( IMPLIES( AND2( *outputPolylineCount > 0, outputVertices ),
                  inBoundsDouble( outputVertices, *outputVertexCount,
                            (const double (*)[2]) bounds ) ) );

  POST0( IMPLIES_ELSE( *outputPolylineCount == 0, *outputVertexCount == 0,
           AND3( *outputVertexCount >= 2,
                 IMPLIES( outputCounts,
                          AND2( minimumItemi( outputCounts,
                                              *outputPolylineCount ) >= 2,
                                sum( outputCounts, *outputPolylineCount )
                                  == *outputVertexCount ) ),
                 IMPLIES( outputVertices,
                          inBoundsDouble( outputVertices, *outputVertexCount,
                                          (const double (*)[2]) bounds )))));
}



/******************************************************************************
PURPOSE: daysInMonth - The number of days in the given [year, month].
INPUTS:  int year  4-digit year.
         int month [1, 12].
RETURNS: int [1, 31].
******************************************************************************/

int daysInMonth( int year, int month ) {
  PRE02( IN_RANGE( year, 1800, 2147 ), IN_RANGE( month, 1, 12 ) );
  int result = 31;
  int isLeapYear = 0;

  switch ( month ) {
  case 2:
    isLeapYear = year % 4 == 0 && ! ( year % 100 == 0 && year % 400 != 0 );
    result = isLeapYear ? 29 : 28;
    break;
  case 4:
  case 6:
  case 9:
  case 11:
    result = 30;
    break;
  default:
    result = 31;
    break;
  }

  POST0( IN_RANGE( result, 1, 31 ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidYYYYMMDD - Is the date yyyymmdd valid?
INPUTS:  int yyyymmdd  E.g., 20090501 is May 1, 2009.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidYYYYMMDD( int yyyymmdd ) {
  int yyyy   = yyyymmdd / 10000;
  int mm     = yyyymmdd / 100 % 100;
  int dd     = yyyymmdd % 100;
  const int validYearMonth =
    AND2( IN_RANGE( yyyy, 1800, 2147 ), IN_RANGE( mm, 1, 12 ) );
  const int daysInThisMonth = validYearMonth ? daysInMonth( yyyy, mm ) : 0;
  const int result = AND2( validYearMonth, IN_RANGE( dd, 1, daysInThisMonth ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: incrementDate - The date (yyyymmdd) incremented by the given of days.
INPUTS:  int yyyymmdd  E.g., 20090501 is May 1, 2009.
         int days      Number of days to increment from yyyymmdd.
RETURNS: int incremented yyyymmdd or maximum value upon overflow.
******************************************************************************/

int incrementDate( int yyyymmdd, int days ) {
  PRE02( isValidYYYYMMDD( yyyymmdd ), days >= 0 );
  int result = yyyymmdd;
  int yyyy   = yyyymmdd / 10000;
  int mm     = yyyymmdd / 100 % 100;
  int dd     = yyyymmdd % 100;
  int daysInThisMonth = daysInMonth( yyyy, mm );
  int d = 0;

  for ( d = 0; d < days; ++d ) {
    ++dd;

    if ( dd > daysInThisMonth ) {
      dd = 1;
      ++mm;

      if ( mm > 12 ) {
        mm = 1;
        ++yyyy;

        if ( yyyy > 2147 ) {
          yyyy = 2147;
          mm = 12;
          dd = 31;
          d = days - 1;
        }
      }

      daysInThisMonth = daysInMonth( yyyy, mm );
    }
  }

  result = yyyy * 10000 + mm * 100 + dd;
  POST0( isValidYYYYMMDD( result ) );
  return result;
}



/******************************************************************************
PURPOSE: isValidYYYYMMDDHH - Is the date yyyymmddhh valid?
INPUTS:  int yyyymmddhh  E.g., 2009050112 is May 1, 2009 at noon.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidYYYYMMDDHH( int yyyymmddhh ) {
  const int yyyymmdd = yyyymmddhh / 100;
  const int hh = yyyymmddhh % 100;
  const int result = AND2( isValidYYYYMMDD( yyyymmdd ), IN_RANGE( hh, 0, 23));
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: incrementDateTime - The date (yyyymmddhh) incremented by hours.
INPUTS:  int yyyymmddhh  E.g., 2009050112 is May 1, 2009 noon.
         int hours      Number of hours to increment from yyyymmddhh.
RETURNS: int incremented yyyymmddhh or maximum value upon overflow.
******************************************************************************/

int incrementDateTime( int yyyymmddhh, int hours ) {
  PRE02( isValidYYYYMMDDHH( yyyymmddhh ), hours >= 0 );
  const int days = hours / 24;
  int yyyymmdd = yyyymmddhh / 100;
  const int remainingHours = hours % 24;
  int hh = yyyymmddhh % 100 + remainingHours;
  int result = 0;
  yyyymmdd = incrementDate( yyyymmdd, days );

  if ( hh > 23 ) {
    yyyymmdd = incrementDate( yyyymmdd, 1 );
    hh -= 24;
  }

  result = yyyymmdd * 100 + hh;
  POST0( isValidYYYYMMDDHH( result ) );
  return result;
}



/******************************************************************************
PURPOSE: endOfDay - Compute timestamp at of end of day.
INPUTS: const int yyyymmddhh   Starting timestamp.
RETURNS: int yyyymmddhh at end of day.
******************************************************************************/

int endOfDay( int yyyymmddhh ) {
  PRE0( isValidYYYYMMDDHH( yyyymmddhh ) );
  const int result = yyyymmddhh / 100 * 100 + 23;
  POST03( isValidYYYYMMDDHH( result ), result / 100 == yyyymmddhh / 100,
          result % 100 == 23 );
  return result;
}



/******************************************************************************
PURPOSE: endOfMonth - Compute timestamp at of end of month.
INPUTS: const int yyyymmddhh   Starting timestamp.
RETURNS: int yyyymmddhh at end of month.
******************************************************************************/

int endOfMonth( int yyyymmddhh ) {
  PRE0( isValidYYYYMMDDHH( yyyymmddhh ) );
  const int yyyymm = yyyymmddhh / 10000;
  const int yyyy = yyyymm / 100;
  const int mm   = yyyymm % 100;
  const int dd = daysInMonth( yyyy, mm );
  const int result = ( yyyymm * 100 + dd ) * 100 + 23;
  POST04( isValidYYYYMMDDHH( result ), result / 10000 == yyyymmddhh / 10000,
          result / 100 % 100
            == daysInMonth( result / 1000000, result / 10000 % 100 ),
          result % 100 == 23 );
  return result;
}



/******************************************************************************
PURPOSE: endOfMYear - Compute timestamp at of end of year.
INPUTS: const int yyyymmddhh   Starting timestamp.
RETURNS: int yyyymmddhh at end of month.
******************************************************************************/

int endOfYear( int yyyymmddhh ) {
  PRE0( isValidYYYYMMDDHH( yyyymmddhh ) );
  const int result =
            ( ( yyyymmddhh / 1000000 * 100 + 12 ) * 100 + 31 ) * 100 + 23;
  POST02( isValidYYYYMMDDHH( result ),
          result ==
            ( ( yyyymmddhh / 1000000 * 100 + 12 ) * 100 + 31 ) * 100 + 23 );
  return result;
}



/******************************************************************************
PURPOSE: indexOfString - Index of string in strings[] or -1 if not present.
INPUTS:  const char* string           String to search for.
         const char* const strings[]  Strings to search.
         int count                    Size of strings[].
RETURNS: int index of string in strings[], else -1 if not present.
******************************************************************************/

int indexOfString(const char* string, const char* const strings[], int count) {

  PRE08( string, *string, strings, strings[ 0 ], *strings[ 0 ],
         count > 0, strings[ count - 1 ], *strings[ count - 1 ] );

  int result = -1;
  int index = 0;

  do {

    if ( strcmp( string, strings[ index ] ) == 0 ) {
      result = index;
      index = count;
    }

    ++index;
  } while ( index < count );

  POST0( OR2( result == -1,
              AND2( result >= 0,
                    strcmp( string, strings[ result ] ) == 0 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: matchesWord - Does word appear in space-delimited words?
INPUTS:  const char* word   String to search for.
         const char* words  String of space-delimited words to search.
RETURNS: int 1 if found, else 0.
******************************************************************************/

int matchesWord( const char* word, const char* const words ) {
  PRE06( word, *word, ! strchr( word, ' ' ),
         words, *words, strchr( words, ' ' ) );
  int result = 0;
  const char* nextWord = words;
  const size_t length = strlen( word );

  do {
    nextWord = strstr( nextWord, word );

    if ( nextWord ) {
      const int atStart = nextWord == words;
      const int nothingBefore = OR2( atStart, *( nextWord - 1 ) == ' ' );

      if ( nothingBefore ) {
        const char after = *( nextWord + length );
        result = IN3( after, '\0', ' ' );
      }

      ++nextWord;
    }
  } while ( AND2( nextWord, ! result ) );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: matchesPattern - Does word match pattern with integer format(s)?
INPUTS:  const char* string   String  to check, e.g., "WN_85_%"".
         const char* pattern  Pattern to match, e.g., "WN_%02d_%%".
RETURNS: int 1 if matched, else 0.
******************************************************************************/

int matchesPattern( const char* const string, const char* const pattern ) {
  PRE05( string, *string, pattern, *pattern,
         OR3( strstr( pattern, "%02d" ),
              strstr( pattern, "%04d" ),
              strstr( pattern, "%d" ) ) );
  int result = 0;
  const char* s = string;
  const char* p = pattern;

  do {
    int sSkip = 1; /* Number of characters to advance s after matched compare*/
    int pSkip = 1; /* Number of characters to advance p after matched compare*/

    if ( ! strncmp( p, "%d", 3 ) ) {
      int i = 1;
      result = isdigit( *s );

      while ( isdigit( *s + i ) ) {
        ++sSkip; /* Skip possible digits. */
        ++i;
      }

      pSkip = 2; /* Skip format. */
    } else if ( ! strncmp( p, "%02d", 4 ) ) {
      result = AND2( isdigit( *s ), isdigit( s[ 1 ] ) );
      sSkip = 2; /* Skip possible 2 digits. */
      pSkip = 4; /* Skip format. */
    } else if ( ! strncmp( p, "%04d", 4 ) ) {
      result =
        AND7( isdigit( *s ), s[ 1 ] != '\0', s[ 2 ] != '\0', s[ 3 ] != '\0',
              isdigit( s[ 1 ] ), isdigit( s[ 2 ] ), isdigit( s[ 3 ] ) );

      sSkip = 4; /* Skip possible 4 digits. */
      pSkip = 4; /* Skip format. */
    } else if ( AND2( *p == '%', p[ 1 ] == '%' ) ) {
      result = *s == '%';
      sSkip = 1; /* Skip %. */
      pSkip = 2; /* Skip %%. */
    } else {
      result = *s == *p;
    }

    if ( result ) {
      s += sSkip;
      p += pSkip;
    }

  } while ( AND3( result, *s, *p ) );

  result = AND3( result, *s == '\0', *p == '\0' );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: makeRGB - Initialize and return an RGB.
INPUTS:  unsigned char r  Red component.
         unsigned char g  Green component.
         unsigned char b  Blue component.
         const char* s    Name component.
RETURNS: RGB initialized to the given components.
******************************************************************************/

RGB makeRGB(unsigned char r, unsigned char g, unsigned char b, const char* s) {
  const RGB result = { r, g, b, s };
  return result;
}



/******************************************************************************
PURPOSE: isValidColor - Is Color in the unit RGB cube?
INPUTS:  Color color  Color to check.
RETURNS: int 1 if valid, else 0.
******************************************************************************/

int isValidColor( Color color ) {
  const int result =
    AND3( IN_RANGE( color.r, 0.0, 1.0 ),
          IN_RANGE( color.g, 0.0, 1.0 ),
          IN_RANGE( color.b, 0.0, 1.0 ) );
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: pointValue - Measure of point data.
INPUTS:  const PointData* const pointData  Point to check.
RETURNS: double measure or magnitude if measure2 != -9999.0.
******************************************************************************/

double pointValue( const PointData* const pointData ) {
  PRE0( pointData );
  const double measure = pointData->measure;
  const double measure2 = pointData->measure2;
  const double result =
    measure2 > -9999.0 ?
      sqrt( measure * measure + measure2 * measure2 )
    : measure;
  POST0( ! isNan( result ) );
  return result;
}



/******************************************************************************
PURPOSE: pointDataRange - Range of measure of point data.
INPUTS:  const int count                     Number of points in array.
         const PointData[ count ] pointData  Points to check.
OUTPUTS: int* indexOfMinimum  Index of minimum measure.
         int* indexOfMaximum  Index of maximum measure.
******************************************************************************/

void pointDataRange( const int count, const PointData pointData[],
                     int* indexOfMinimum, int* indexOfMaximum ) {
  PRE04( count > 0, pointData, indexOfMinimum, indexOfMaximum );
  int index = 0;
  double minimum = pointValue( pointData );
  double maximum = minimum;

  *indexOfMinimum = *indexOfMaximum = 0;

  for ( index = 1; index < count; ++index ) {
    const double value = pointValue( pointData + index );

    if ( value < minimum ) {
      minimum = value;
      *indexOfMinimum = index;
    } else if ( value > maximum ) {
      maximum = value;
      *indexOfMaximum = index;
    }
  }

  POST03( IN_RANGE( *indexOfMinimum, 0, count - 1 ),
          IN_RANGE( *indexOfMaximum, 0, count - 1 ),
          pointValue( pointData + *indexOfMinimum ) <=
            pointValue( pointData + *indexOfMaximum ) );
}



/******************************************************************************
PURPOSE: pointMatches - Does point match timestamp?
INPUTS:  const PointData* const pointData  Point to check.
         const int yyyymmddhh              Timestamp to match.
         const int timestepType            HOURLY, DAILY, MONTHLY, YEARLY.
RETURNS: int 1 if matched, else 0.
******************************************************************************/

int pointMatches( const PointData* const pointData,
                  const int yyyymmddhh, const int timestepType ) {

  PRE03( pointData,
         isValidYYYYMMDDHH( yyyymmddhh ),
         IS_VALID_TIMESTEP_TYPE( timestepType ) );
  const int yyyymmdd = pointData->yyyymmdd;
  const char* const name = pointData->name;
  /* Non-timestamped data, e.g., station ids match. */
  int result =
    OR3( yyyymmdd == -1,
        ! strcmp( name, "station" ),
        ! strcmp( name, "uid" ) );

  if ( ! result ) {

    switch ( timestepType ) {
    case YEARLY:
      result = yyyymmddhh / 1000000 == yyyymmdd / 10000;
      break;
    case MONTHLY:
      result = yyyymmddhh / 10000 == yyyymmdd / 100;
      break;
    case DAILY:
      result = yyyymmddhh / 100 == yyyymmdd;
      break;
    default:
      CHECK( timestepType == HOURLY );
      {
        const int hhmmss = pointData->hhmmss; /* hhmmss == -1 is daily data. */
        const int hour = hhmmss == -1 ? -1 : hhmmss / 10000;
        result = OR2( hour == -1, yyyymmddhh == yyyymmdd * 100 + hour );
      }
      break;
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



#if 0

/******************************************************************************
PURPOSE: dataColor - Map value in range [minimum, maximum] to
         RGB color cube range: blue-cyan-green-yellow-red.
INPUTS:  double value    Data value to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color dataColor( double value, double minimum, double maximum ) {
  PRE03( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ) );
  Color result = { 0.0, 0.0, 0.0 };

  if ( value != -9999.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    if ( t <= 0.25 ) {
      result.g = t * 4.0;
      result.b = 1.0;
    } else if ( t <= 0.5 ) {
      t = ( t - 0.25 ) * 4.0;
      result.g = 1.0;
      result.b = 1.0 - t;
    } else if ( t <= 0.75 ) {
      t = ( t - 0.5 ) * 4.0;
      result.g = 1.0;
      result.r = t;
    } else {
      t = ( t - 0.75 ) * 4.0;
      result.r = 1.0;
      result.g = 1.0 - t;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}

#else

/******************************************************************************
PURPOSE: dataColor - Map value in range [minimum, maximum] to
         one of the discrete colors: blue-green-yellow-orange-red.
INPUTS:  double value    Data value to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color dataColor( double value, double minimum, double maximum ) {
  PRE03( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ) );
  Color result = { 0.0, 0.0, 0.0 };

  if ( value != -9999.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    if ( t < 0.2 ) {
      result.b = 1.0;
    } else if ( t < 0.4 ) {
      result.g = 1.0;
    } else if ( t < 0.6 ) {
      result.g = 1.0;
      result.r = 1.0;
    } else if ( t < 0.8 ) {
      result.g = 0.5;
      result.r = 1.0;
    } else {
      result.r = 1.0;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}

#endif



/******************************************************************************
PURPOSE: modulo6Color - Map value (modulo 6) to
         RGB color cube range: blue-cyan-green-yellow-red-magenta.
INPUTS:  double value    Data value to map.
         double unused1
         double unused2
RETURNS: Color mapped value.
******************************************************************************/

Color modulo6Color( double value, double unused1, double unused2 ) {
  PRE03( ! isNan( value ), ! isNan( unused1 ), ! isNan( unused2 ) );
  Color result = { 0.0, 0.0, 0.0 };

  if ( value != -9999.0 ) {
    const int ivalue = (int) value;
    const int modValue = ivalue % 7;

    switch ( modValue ) {
    case 0:
      result.b = 1.0;
      break;
    case 1:
      result.b = 1.0;
      result.g = 1.0;
      break;
    case 2:
      result.g = 1.0;
      break;
    case 3:
      result.r = 1.0;
      result.g = 1.0;
      break;
    case 4:
      result.r = 1.0;
      break;
    default:
      result.r = 1.0;
      result.b = 1.0;
      break;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: grayColor - Map value in range [minimum, maximum] to gray.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color grayColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    {
      const double colorMinimum = 1.0 / 8.0;
      const double colorMaximum = 1.0; /* 1.0 - colorMinimum; */
      const double colorRange = colorMaximum - colorMinimum;
      result.r = result.g = result.b = colorMinimum + t * colorRange;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: greenColor - Map value in range [minimum, maximum] to green.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color greenColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = grayColor( value, minimum, maximum );
  result.r = result.b = 0.0;
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: greenGrayColor - Map value in range [minimum, maximum] to green-gray.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color greenGrayColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    result.r = result.b = t;
    result.g = 0.75;

    if ( t > 0.5 ) {
      result.g += ( t - 0.5 ) * 0.5;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: darkGreenGrayColor - Map value in range [minimum, maximum] to
         dark_green-gray.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color darkGreenGrayColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    if ( t <= 0.5 ) {
      result.r = result.b = t * 0.5;
      result.g = 0.125 + 0.375 * ( t + t );
    } else {
      const double t2 = t - 0.5;
      result.r = result.b = 0.25 + t2 * 1.5;
      result.g = 0.5 + t2;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: cyanGrayColor - Map value in range [minimum, maximum] to cyan-gray.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color cyanGrayColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    if ( t <= 0.5 ) {
      result.r = t * 0.5;
      result.g = result.b = 0.25 + t * 0.5;
    } else {
      const double t2 = t - 0.5;
      result.r = 0.25 + t2 * 1.5;
      result.g = result.b = 0.5 + t2;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: blueGrayColor - Map value in range [minimum, maximum] to blue-gray.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color blueGrayColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    if ( t <= 0.5 ) {
      result.r = result.g = t * 0.5;
      result.b = 0.25 + t * 0.5;
    } else {
      const double t2 = t - 0.5;
      result.r = result.g = 0.25 + t2 * 1.5;
      result.b = 0.5 + t2;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: brownGrayColor - Map value in range [minimum, maximum] to brown-gray.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color brownGrayColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    if ( t <= 0.5 ) {
      const double t2 = t + t;
      result.r = 0.125 + 0.375 * t2;
      result.g = 0.125 + 0.125 * t2;
      result.b = 0.125 * t2;
    } else {
      const double t2 = t - 0.5;
      result.r = t;
      result.g = 0.25 + t2 * 1.5;
      result.b = 0.125 + t2 * 1.75;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: tanGrayColor - Map value in range [minimum, maximum] to tan-gray.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color tanGrayColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double range = maximum - minimum;
    double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    if ( t <= 0.5 ) {
      result.r = result.g = 0.25 + t * 0.5;
      result.b = t * 0.5;
    } else {
      const double t2 = t - 0.5;
      result.r = result.g = t;
      result.b = 0.25 + t2 * 1.5;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: elevationColor - Map value in range [minimum, maximum] to
         RGB color for elevation.
INPUTS:  double value    Data value to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color elevationColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value != -9999.0 ) {

    if ( value == 0.0 ) { /* Water edge: light cyan. */
      result.r = 0.5;
      result.g = 0.85;
      result.b = 0.99;
    } else {

      if ( value < 0.0 ) {
        maximum = 0.0;
      } else {
        minimum = 0.0;
      }

      {
        const double range = maximum - minimum;
        double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

        if ( t < 0.0 || isNan( t ) ) {
          t = 0.0;
        } else if ( t > 1.0 ) {
          t = 1.0;
        }

        if ( value < 0.0 ) { /* Water: dark blue to light blue-green. */
          double blue = 0.5;
          double green = 0.0;

          if ( value <= -1000.0 ) {
            blue = 0.15; /* Very deep water: very dark blue. */
          } else if ( value <= -500.0 ) {
            blue = 0.25; /* Deep water: dark blue. */
          } else if ( value >= -5 ) {
            green = 0.25; /* Swimming-pool-depth: swimming pool green. */
          }

          t *= t; /* Non-linear scaling enhances surface variation. */
          result.g = green + t * 0.5;
          result.b = blue + ( 1.0 - blue ) * t;
        } else { /* Land: lighter brown to grayish white. */
          t = sqrt( t ); /* Non-linear scaling enhances surface variation. */
          result.r = 0.5 * 0.58 + 0.39 * t;
          result.g = 0.5 * 0.39 + 0.58 * t;
          result.b = 0.5 * 0.19 + 0.78 * t;
        }
      }
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: seagrassColor - Map value in range [minimum, maximum] to
         RGB color for seagrass.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color seagrassColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value > 0.0 ) {
    const int ivalue = value;
    const int density = ivalue / 10000; /* Extract just 1st of 5 digits. */

    if ( density >= 1 && density <= 4 ) {
      result.g = density * 0.25;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: riskColor - Map value in range [minimum, maximum] to
         RGB color for risk.
INPUTS:  double value    Data value to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color riskColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */
  const double range = maximum - minimum;
  double t = range <= 0.0 ? 0.0 : ( value - minimum ) / range;

  if ( t < 0.0 || isNan( t ) ) {
    t = 0.0;
  } else if ( t > 1.0 ) {
    t = 1.0;
  }

  if ( t <= 0.2 ) {
    result.b = 1.0;
  } else if ( t <= 0.4 ) {
    result.g = 1.0;
  } else if ( t <= 0.6 ) {
    result.r = 1.0;
    result.g = 1.0;
  } else if ( t <= 0.8 ) {
    result.r = 1.0;
    result.g = 0.5;
  } else {
    result.r = 1.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: riskColor4 - Map value [1-4] in range [minimum, maximum] to
         RGB color for risk.
INPUTS:  double value    Data value to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color riskColor4( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
        minimum <= maximum );
  const Color result = riskColor( value + 1.0, minimum, maximum + 1.0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: riskColorText - Map risk text value to color for vulnerability data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color riskColorText( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( ! strcmp( text, "Very High" ) ) {
    result.r = 1.0;
  } else if ( ! strcmp( text, "High" ) ) {
    result.r = 1.0;
    result.g = 0.5;
  } else if ( ! strcmp( text, "Moderate" ) ) {
    result.r = 1.0;
    result.g = 1.0;
  } else if ( ! strcmp( text, "Low" ) ) {
    result.g = 1.0;
  } else if ( ! strcmp( text, "Very Low" ) ) {
    result.b = 1.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: soilColorText - Map text value to color for soil data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color soilColorText( const char* text ) {
  PRE02( text, *text );
  const char c = tolower( *text );
  const double t = c - 'a';
  const Color result = soilColor( t );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: soilColorDarkLight - Map value in range [minimum, maximum] to
         RGB color for soil.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color soilColorDarkLight( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double t = ( value - minimum ) / ( maximum - minimum );
    result = soilColor( t );
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: soilColorLightDark7 - Map value in range [minimum, maximum] to
         RGB color for soil.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color soilColorLightDark7( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double t = ( value - minimum ) / ( maximum - minimum );
    result = soilColor7( 1.0 - t );
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: soilColorLightDark - Map value in range [minimum, maximum] to
         RGB color for soil.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color soilColorLightDark( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double t = ( value - minimum ) / ( maximum - minimum );
    result = soilColor( 1.0 - t );
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: soilColorDarkLight4 - Map value in range [minimum, maximum] to
         RGB color for soil.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color soilColorDarkLight4( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value >= 0.0 ) {
    const double t = ( value - minimum ) / ( maximum - minimum );
    result = soilColor4( t );
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: precipitationColor - Map value in range [minimum, maximum] to
         RGB color for precipitation.
INPUTS:  double value    Data value (in mm) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color precipitationColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value != -9999.0 ) {
    const double middle = ( minimum + maximum ) * 0.5;

    if ( value <= middle ) {
      double t = value / maximum;

      if ( t < 0.0 || isNan( t ) ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      {
        const double oneMinusT = 1.0 - t;
        const double oneOver255 = 1.0 / 255.0;
        result.r = 219.0 * oneMinusT;
        result.g = 187.0 + t * ( 191.0 - 187.0 );
        result.b = 127.0 * oneMinusT;
        result.r *= oneOver255;
        result.g *= oneOver255;
        result.b *= oneOver255;
      }
    } else {
      double t = ( value - middle ) / ( maximum - middle );

      if ( t < 0.0 || isNan( t ) ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      result.r = 0.0;
      result.g = 191.0 - t * ( 191.0 - 61.0 );
      result.g /= 255.0;
      result.b = 0.0;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: populationDensityColor - Map value in range [minimum, maximum] to
         RGB color for population density.
INPUTS:  double value    Data value (in people / km2) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
NOTE:    Uses a log-normal scale so very high population counties won't
         dominate color range:
         black->gray  for value in [0, 10] (people per square km)
         gray->yellow for value in [10, 100]
         yellow->red  for value in [100, 1000+].
******************************************************************************/

Color populationDensityColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value > 0.0 ) {

    if ( value < 10.0 ) { /* Ramp black rgb=10 to gray rgb=127: */
      double t = value / 10.0;

      if ( t < 0.0 || isNan( t ) ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      result.r = 10.0 + t * ( 127.0 - 10.0 );
      result.r /= 255.0;
      result.g = result.r;
      result.b = result.r;

    } else if ( value < 100.0 ) { /* Ramp gray rgb=127 to yellow: */
      double t = value / ( 100.0 - 10.0 );

      if ( t < 0.0 || isNan( t ) ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      {
        const double oneMinusT = 1.0 - t;
        const double oneOver255 = 1.0 / 255.0;
        result.r = 127.0 + t * ( 255.0 - 127.0 );
        result.r *= oneOver255;
        result.g = result.r;
        result.b = 127.0 * oneMinusT;
        result.b *= oneOver255;
      }
    } else { /* Ramp yellow to red: */
      double t = value / ( 1000.0 - 100.0 );

      if ( t < 0.0 || isNan( t ) ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      {
        const double oneMinusT = 1.0 - t;
        result.r = 1.0;
        result.g = oneMinusT;
      }
    }
  }

  DEBUG2( fprintf( stderr, "pop-color: %f -> RGB<%f %f %f>\n",
                  value, result.r, result.g, result.b ); )
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: populationColor - Map value in range [minimum, maximum] to
         RGB color for population.
INPUTS:  double value    Data value (in people / km2) to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
NOTE:    Uses a log-normal scale so very high population counties won't
         dominate color range. For t (normalized value):
         black->gray  for t in [0, 0.00025]
         gray->yellow for t in [0.00025, 0.0025]
         yellow->red  for t in [0.0025, 1].
******************************************************************************/

Color populationColor( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  Color result = { 0.0, 0.0, 0.0 }; /* Missing data is black. */

  if ( value > 0.0 ) {
    double t = ( value - minimum ) / ( maximum - minimum );

    if ( t < 0.00025 ) { /* Ramp black rgb=10 to gray rgb=127: */
      t *= 4000.0; /* Normalize t from range [0, 0.00025] to [0, 1]: */

      if ( t < 0.0 ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      result.r = 10.0 + t * ( 127.0 - 10.0 );
      result.r /= 255.0;
      result.g = result.r;
      result.b = result.r;

    } else if ( t < 0.0025 ) { /* Ramp gray rgb=127 to yellow: */

      /* Normalize t from range [0.00025, 0.0025] to [0, 1]: */

      t -= 0.00025;       /* t in [0, 0.00225]. */
      t *= 1.0 / 0.00225; /* t in [0, 1]. */

      if ( t < 0.0 || isNan( t ) ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      {
        const double oneMinusT = 1.0 - t;
        const double oneOver255 = 1.0 / 255.0;
        result.r = 127.0 + t * ( 255.0 - 127.0 );
        result.r *= oneOver255;
        result.g = result.r;
        result.b = 127.0 * oneMinusT;
        result.b *= oneOver255;
      }
    } else { /* Ramp yellow to red: */

      /* Normalize t from range [0.0025, 1] to [0, 1]: */

      t -= 0.0025;       /* t in [0, 0.9975]. */
      t *= 1.0 / 0.9975; /* t in [0, 1]. */

      if ( t < 0.0 || isNan( t ) ) {
        t = 0.0;
      } else if ( t > 1.0 ) {
        t = 1.0;
      }

      {
        const double oneMinusT = 1.0 - t;
        result.r = 1.0;
        result.g = oneMinusT;
      }
    }
  }

  DEBUG2( fprintf( stderr, "pop-color: %f -> RGB<%f %f %f>\n",
                  value, result.r, result.g, result.b ); )
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: lexicographicTextColor - Map text to Color based on first letter
         compared to 'a'..'z' or '0'..'9'.
INPUTS:  const char* text  Text to map.
RETURNS: Color mapped value.
******************************************************************************/

Color lexicographicTextColor( const char* text ) {
  PRE0( text );
  const char c = tolower( *text );
  const int isNumber = isdigit( c );
  const double t = isNumber ? c - '0' : c - 'a';
  const Color result = dataColor( t, 0.0, isNumber ? 9.0 : 26.0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sedimentColorText - Map text value to color for sediment data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color sedimentColorText( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( strstr( text, "- 0.03" ) ) {
    result.r = 64.0 / 255.0;
    result.g = 31.0 / 255.0;
  } else if ( strstr( text, "- 0.17" ) ) {
    result.r = 96.0 / 255.0;
    result.g = 47.0 / 255.0;
    result.b = 31.0 / 255.0;
  } else if ( strstr( text, "- 0.35" ) ) {
    result.r = 127.0 / 255.0;
    result.g =  64.0 / 255.0;
    result.b =  31.0 / 255.0;
  } else if ( strstr( text, "- 0.36" ) ) {
    result.r = 191.0 / 255.0;
    result.g = 127.0 / 255.0;
    result.b =  64.0 / 255.0;
  } else if ( strstr( text, "- 0.48" ) ) {
    result.r = 222.0 / 255.0;
    result.g = 174.0 / 255.0;
    result.b = 116.0 / 255.0;
  } else if ( strstr( text, "0.48+" ) ) {
    result.r = 255.0 / 255.0;
    result.g = 221.0 / 255.0;
    result.b = 167.0 / 255.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sedimentColor - Map value in range [minimum, maximum] to
         RGB color for sediment data.
INPUTS:  double value    Data value percent to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
         int r          Red color maximum.
         int g          Green color maximum.
         int b          Blue color maximum.
RETURNS: Color mapped value.
******************************************************************************/

Color sedimentColor( double value, double minimum, double maximum,
                     int r, int g, int b ) {
  PRE07( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum,
         IN_RANGE( r, 0, 255 ), IN_RANGE( g, 0, 255 ), IN_RANGE( b, 0, 255 ) );
  Color result = { 0.0, 0.0, 0.0 };

  if ( AND2( value > 0.0, maximum > minimum ) ) {
    double t = ( value - minimum ) / ( maximum - minimum );

    if ( t < 0.0 || isNan( t ) ) {
      t = 0.0;
    } else if ( t > 1.0 ) {
      t = 1.0;
    }

    {
      const double tOver255 = t / 255.0;
      result.r = r * tOver255;
      result.g = g * tOver255;
      result.b = b * tOver255;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sedimentColorMud - Map value in range [minimum, maximum] to
         RGB color for sediment mud data.
INPUTS:  double value    Data value percent to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color sedimentColorMud( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  const Color result = sedimentColor( value, minimum, maximum, 64, 31, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sedimentColorMudu - Map value in range [minimum, maximum] to
         RGB color for sediment mudu data.
INPUTS:  double value    Data value percent to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color sedimentColorMudu( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  const Color result = sedimentColor( value, minimum, maximum, 64, 31, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sedimentColorSand - Map value in range [minimum, maximum] to
         RGB color for sediment sand data.
INPUTS:  double value    Data value percent to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color sedimentColorSand( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  const Color result = sedimentColor( value, minimum, maximum, 127, 64, 31 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sedimentColorGravel - Map value in range [minimum, maximum] to
         RGB color for sediment gravel data.
INPUTS:  double value    Data value percent to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color sedimentColorGravel( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  const Color result = sedimentColor( value, minimum, maximum, 255, 221, 167 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sedimentColorRock - Map value in range [minimum, maximum] to
         RGB color for sediment rock data.
INPUTS:  double value    Data value percent to map.
         double minimum  Minimum value for data.
         double maximum  Maximum value for data.
RETURNS: Color mapped value.
******************************************************************************/

Color sedimentColorRock( double value, double minimum, double maximum ) {
  PRE04( ! isNan( value ), ! isNan( minimum ), ! isNan( maximum ),
         minimum <= maximum );
  const Color result = sedimentColor( value, minimum, maximum, 255, 238, 211 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nitrogenColorText - Map text value to color for nitrogen source data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color nitrogenColorText( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( ! strcmp( text, "BKGD" ) ) {
    result.r = 0.5;
    result.g = 1.0;
    result.b = 0.5;
  } else if ( ! strcmp( text, "ATDP" ) ) {
    result.g = 1.0;
    result.b = 1.0;
  } else if ( ! strcmp( text, "FERT" ) ) {
    result.r = 1.0;
    result.g = 0.5;
  } else if ( ! strcmp( text, "MANR" ) ) {
    result.r = 0.5;
    result.g = 0.25;
  } else if ( ! strcmp( text, "WW" ) ) {
    result.g = 0.5;
    result.b = 1.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: wordColor - Map word to color .
INPUTS:  const char* word  Word to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color wordColor( const char* word ) {
  PRE02( word, *word );
  Color result = { 0.0, 0.0, 0.0 };

  if ( ! strcmp( word, "COLD" ) ) {
    result.b = 1.0;
  } else if ( ! strcmp( word, "COOL" ) ) {
    result.g = 1.0;
    result.b = 1.0;
  } else if ( ! strcmp( word, "WARM" ) ) {
    result.r = 1.0;
  } else if ( ! strcmp( word, "TIDL" ) ) {
    result.g = 1.0;
  } else if ( ! strcmp( word, "NTDL" ) ) {
    result.r = 1.0;
    result.g = 1.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: blackstoneHRUColor - Map Blackstone watershed HRU_ID to color .
INPUTS:  const int hru  HRU_ID to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/WSM/blackstone_hru_legend.png.
******************************************************************************/

Color blackstoneHRUColor( double value, double unused1_, double unused2_ ) {
  PRE0( value >= 0.0 );
  Color result = { 0.5, 0.5, 0.5 }; /* Unknown/Other is gray. */
  const int hru = (int) ( value + 0.5 );

  if ( hru == 1 ) { /* Commercial-Industrial: */
    result.r = 255.0 / 255.0;
    result.g = 127.0 / 255.0;
    result.b =   0.0 / 255.0;
  } else if ( IN3( hru, 2, 10 ) ) { /* Open: */
    result.r = 200.0 / 255.0;
    result.g = 175.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 3, 11 ) ) { /* Forest: */
    result.r =  36.0 / 255.0;
    result.g = 103.0 / 255.0;
    result.b =  51.0 / 255.0;
  } else if ( IN3( hru, 4, 12 ) ) { /* Medium-densty residential: */
    result.r = 215.0 / 255.0;
    result.g = 150.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 7, 15 ) ) { /* High-densty residential: */
    result.r = 240.0 / 255.0;
    result.g =   0.0 / 255.0;
    result.b =  24.0 / 255.0;
  } else if ( hru == 100 ) { /* Water: */
    result.r =  79.0 / 255.0;
    result.g = 107.0 / 255.0;
    result.b = 161.0 / 255.0;
  } else if ( hru == 100 ) { /* Wetland: */
    result.r = 120.0 / 255.0;
    result.g = 160.0 / 255.0;
    result.b = 160.0 / 255.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: charlesHRUColor - Map Charles watershed HRU_ID to color .
INPUTS:  const int hru  HRU_ID to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/WSM/charles_hru_legend.png.
******************************************************************************/

Color charlesHRUColor( double value, double unused1_, double unused2_ ) {
  PRE0( value >= 0.0 );
  Color result = { 0.5, 0.5, 0.5 }; /* Unknown/Other is gray. */
  const int hru = (int) ( value + 0.5 );
  const int firstDigit = hru / 10;

  if ( firstDigit <= 9 ) {

    switch ( firstDigit ) {
    case 0: /* Water: */
      result.r =  79.0 / 255.0;
      result.g = 107.0 / 255.0;
      result.b = 161.0 / 255.0;
      break;
    case 1: /* Open: */
      result.r = 200.0 / 255.0;
      result.g = 175.0 / 255.0;
      result.b = 127.0 / 255.0;
      break;
    case 2: /* Forest: */
      result.r =  36.0 / 255.0;
      result.g = 103.0 / 255.0;
      result.b =  51.0 / 255.0;
      break;
    case 3: /* Forested wetland: */
      result.r = 190.0 / 255.0;
      result.g = 215.0 / 255.0;
      result.b = 127.0 / 255.0;
      break;
    case 4: /* Non-forested wetland: */
      result.r = 120.0 / 255.0;
      result.g = 160.0 / 255.0;
      result.b = 160.0 / 255.0;
      break;
    case 5: /* Low-densty residential: */
      result.r = 240.0 / 255.0;
      result.g = 210.0 / 255.0;
      result.b = 200.0 / 255.0;
      break;
    case 6: /* Medium-densty residential: */
      result.r = 215.0 / 255.0;
      result.g = 150.0 / 255.0;
      result.b = 127.0 / 255.0;
      break;
    case 7: /* High-densty residential: */
      result.r = 240.0 / 255.0;
      result.g =   0.0 / 255.0;
      result.b =  24.0 / 255.0;
      break;
    case 8: /* Multi-family residential: */
      result.r = 164.0 / 255.0;
      result.g =   0.0 / 255.0;
      result.b =  16.0 / 255.0;
      break;
    default: /* Commercial-Industrial: */
      CHECK( firstDigit == 9 );
      result.r = 255.0 / 255.0;
      result.g = 127.0 / 255.0;
      result.b =   0.0 / 255.0;
      break;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: farmingtonHRUColor - Map Farmington watershed HRU_ID to color .
INPUTS:  const int hru  HRU_ID to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/WSM/farmington_hru_legend.png.
******************************************************************************/

Color farmingtonHRUColor( double value, double unused1_, double unused2_ ) {
  PRE0( value >= 0.0 );
  Color result = { 0.5, 0.5, 0.5 }; /* Unknown/Other is gray. */
  const int hru = (int) ( value + 0.5 );
  const int lastTwoDigits = hru % 100;

  if ( lastTwoDigits < 9 ) {

    switch ( lastTwoDigits ) {
    case 1: /* Forest: */
      result.r =  36.0 / 255.0;
      result.g = 103.0 / 255.0;
      result.b =  51.0 / 255.0;
      break;
    case 2: /* Agriculture/Grass: */
      result.r = 170.0 / 255.0;
      result.g = 114.0 / 255.0;
      result.b =  46.0 / 255.0;
      break;
    case 3: /* Impervious: */
      result.r = 175.0 / 255.0;
      result.g = 175.0 / 255.0;
      result.b = 175.0 / 255.0;
      break;
    case 4: /* Wetland: */
      result.r = 120.0 / 255.0;
      result.g = 160.0 / 255.0;
      result.b = 160.0 / 255.0;
      break;
    case 5: /* High-densty residential: */
      result.r = 240.0 / 255.0;
      result.g =   0.0 / 255.0;
      result.b =  24.0 / 255.0;
      break;
    case 6: /* Medium-densty residential: */
      result.r = 215.0 / 255.0;
      result.g = 150.0 / 255.0;
      result.b = 127.0 / 255.0;
      break;
    case 7: /* Barren: */
      result.r = 200.0 / 255.0;
      result.g = 175.0 / 255.0;
      result.b = 127.0 / 255.0;
      break;
    case 8: /* Major road: */
      result.r = 200.0 / 255.0;
      result.g = 200.0 / 255.0;
      result.b = 200.0 / 255.0;
      break;
    default:

      if ( hru == 1000 ) { /* Water: */
        result.r =  79.0 / 255.0;
        result.g = 107.0 / 255.0;
        result.b = 161.0 / 255.0;
      }

      break;
    }
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: ipswichHRUColor - Map Ipswich watershed HRU_ID to color .
INPUTS:  const int hru  HRU_ID to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/WSM/ipswich_hru_legend.png.
******************************************************************************/

Color ipswichHRUColor( double value, double unused1_, double unused2_ ) {
  PRE0( value >= 0.0 );
  Color result = { 0.5, 0.5, 0.5 }; /* Unknown/Other is gray. */
  const int hru = (int) ( value + 0.5 );

  if ( IN4( hru, 1, 8, 14 ) ) { /* Forest: */
    result.r =  36.0 / 255.0;
    result.g = 103.0 / 255.0;
    result.b =  51.0 / 255.0;
  } else if ( IN4( hru, 2, 9, 15 ) ) { /* Open: */
    result.r = 200.0 / 255.0;
    result.g = 175.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 3, 10 ) ) { /* Low-densty residential: */
    result.r = 240.0 / 255.0;
    result.g = 210.0 / 255.0;
    result.b = 200.0 / 255.0;
  } else if ( IN3( hru, 5, 12 ) ) { /* High-densty residential: */
    result.r = 240.0 / 255.0;
    result.g =   0.0 / 255.0;
    result.b =  24.0 / 255.0;
  } else if ( hru == 7 ) { /* Commercial-Industrial: */
    result.r = 255.0 / 255.0;
    result.g = 127.0 / 255.0;
    result.b =   0.0 / 255.0;
  } else if ( hru == 100 ) { /* Water: */
    result.r =  79.0 / 255.0;
    result.g = 107.0 / 255.0;
    result.b = 161.0 / 255.0;
  } else if ( hru == 201 ) { /* Forest Wetland: */
    result.r = 190.0 / 255.0;
    result.g = 215.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( hru == 202 ) { /* Nonforest Wetland: */
    result.r = 120.0 / 255.0;
    result.g = 160.0 / 255.0;
    result.b = 160.0 / 255.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: pawcatuckHRUColor - Map Pawcatuck watershed HRU_ID to color .
INPUTS:  const int hru  HRU_ID to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/WSM/pawcatuck_hru_legend.png.
******************************************************************************/

Color pawcatuckHRUColor( double value, double unused1_, double unused2_ ) {
  PRE0( value >= 0.0 );
  Color result = { 0.5, 0.5, 0.5 }; /* Unknown/Other is gray. */
  const int hru = (int) ( value + 0.5 );

  if ( IN3( hru, 1, 11 ) ) { /* Commercial-Industrial: */
    result.r = 255.0 / 255.0;
    result.g = 127.0 / 255.0;
    result.b =   0.0 / 255.0;
  } else if ( IN3( hru, 2, 12 ) ) { /* High-densty residential: */
    result.r = 240.0 / 255.0;
    result.g =   0.0 / 255.0;
    result.b =  24.0 / 255.0;
  } else if ( IN3( hru, 3, 13 ) ) { /* Low-densty residential: */
    result.r = 240.0 / 255.0;
    result.g = 210.0 / 255.0;
    result.b = 200.0 / 255.0;
  } else if ( IN3( hru, 4, 14 ) ) { /* Open: */
    result.r = 200.0 / 255.0;
    result.g = 175.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 5, 15 ) ) { /* Forest: */
    result.r =  36.0 / 255.0;
    result.g = 103.0 / 255.0;
    result.b =  51.0 / 255.0;
  } else if ( hru == 6 ) { /* Agriculture: */
    result.r = 170.0 / 255.0;
    result.g = 114.0 / 255.0;
    result.b =  46.0 / 255.0;
  } else if ( hru == 16 ) { /* Rocky: */
    result.r = 140.0 / 255.0;
    result.g = 140.0 / 255.0;
    result.b = 130.0 / 255.0;
  } else if ( hru == 18 ) { /* Nonforest Wetland: */
    result.r = 120.0 / 255.0;
    result.g = 160.0 / 255.0;
    result.b = 160.0 / 255.0;
  } else if ( hru == 19 ) { /* Forest Wetland: */
    result.r = 190.0 / 255.0;
    result.g = 215.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( hru == 100 ) { /* Water: */
    result.r =  79.0 / 255.0;
    result.g = 107.0 / 255.0;
    result.b = 161.0 / 255.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: sudburyHRUColor - Map Sudbury watershed HRU_ID to color .
INPUTS:  const int hru  HRU_ID to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/WSM/sudbury_hru_legend.png.
******************************************************************************/

Color sudburyHRUColor( double value, double unused1_, double unused2_ ) {
  PRE0( value >= 0.0 );
  Color result = { 0.5, 0.5, 0.5 }; /* Unknown/Other is gray. */
  const int hru = (int) ( value + 0.5 );

  if ( IN3( hru, 3, 33 ) ) { /* Commercial-Industrial: */
    result.r = 255.0 / 255.0;
    result.g = 127.0 / 255.0;
    result.b =   0.0 / 255.0;
  } else if ( IN3( hru, 6, 36 ) ) { /* High-densty residential: */
    result.r = 240.0 / 255.0;
    result.g =   0.0 / 255.0;
    result.b =  24.0 / 255.0;
  } else if ( IN3( hru, 9, 39 ) ) { /* Medium-densty residential: */
    result.r = 215.0 / 255.0;
    result.g = 150.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 12, 42 ) ) { /* Low-densty residential: */
    result.r = 240.0 / 255.0;
    result.g = 210.0 / 255.0;
    result.b = 200.0 / 255.0;
  } else if ( IN3( hru, 13, 43 ) ) { /* Open: */
    result.r = 200.0 / 255.0;
    result.g = 175.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 14, 44 ) ) { /* Forest: */
    result.r =  36.0 / 255.0;
    result.g = 103.0 / 255.0;
    result.b =  51.0 / 255.0;
  } else if ( hru == 51 ) { /* Shrub Wetland: */
    result.r = 120.0 / 255.0;
    result.g = 160.0 / 255.0;
    result.b = 160.0 / 255.0;
  } else if ( hru == 52 ) { /* Forest Wetland: */
    result.r = 190.0 / 255.0;
    result.g = 215.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( hru == 100 ) { /* Water: */
    result.r =  79.0 / 255.0;
    result.g = 107.0 / 255.0;
    result.b = 161.0 / 255.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: tauntonHRUColor - Map Taunton watershed HRU_ID to color .
INPUTS:  const int hru  HRU_ID to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/WSM/taunton_hru_legend.png.
******************************************************************************/

Color tauntonHRUColor( double value, double unused1_, double unused2_ ) {
  PRE0( value >= 0.0 );
  Color result = { 0.5, 0.5, 0.5 }; /* Unknown/Other is gray. */
  const int hru = (int) ( value + 0.5 );

  if ( IN3( hru, 1, 9 ) ) { /* Forest: */
    result.r =  36.0 / 255.0;
    result.g = 103.0 / 255.0;
    result.b =  51.0 / 255.0;
  } else if ( IN3( hru, 2, 10 ) ) { /* Open: */
    result.r = 200.0 / 255.0;
    result.g = 175.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 3, 11 ) ) { /* Medium-densty residential: */
    result.r = 215.0 / 255.0;
    result.g = 150.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( IN3( hru, 6, 14 ) ) { /* High-densty residential: */
    result.r = 240.0 / 255.0;
    result.g =   0.0 / 255.0;
    result.b =  24.0 / 255.0;
  } else if ( IN3( hru, 7, 15 ) ) { /* Commercial-Industrial: */
    result.r = 255.0 / 255.0;
    result.g = 127.0 / 255.0;
    result.b =   0.0 / 255.0;
  } else if ( IN3( hru, 8, 16 ) ) { /* Agriculture: */
    result.r = 170.0 / 255.0;
    result.g = 114.0 / 255.0;
    result.b =  46.0 / 255.0;
  } else if ( hru == 17 ) { /* Cranberry bogs: */
    result.r = 134.0 / 255.0;
    result.g = 188.0 / 255.0;
    result.b = 157.0 / 255.0;
  } else if ( hru == 18 ) { /* Forest Wetland: */
    result.r = 190.0 / 255.0;
    result.g = 215.0 / 255.0;
    result.b = 127.0 / 255.0;
  } else if ( hru == 19 ) { /* Nonforest Wetland: */
    result.r = 120.0 / 255.0;
    result.g = 160.0 / 255.0;
    result.b = 160.0 / 255.0;
  } else if ( hru == 100 ) { /* Water: */
    result.r =  79.0 / 255.0;
    result.g = 107.0 / 255.0;
    result.b = 161.0 / 255.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusNaturalWaterColor - Map ICLUS V0 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusNaturalWaterColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 0, 0, 255 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusReservoirsCanalsColor - Map ICLUS V1 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusReservoirsCanalsColor( double percent, double unused1, double maximum) {
  const Color result = categoryColor( percent, maximum, 0, 127, 255 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusWetlandsColor - Map ICLUS V2 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusWetlandsColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 93, 255, 255 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusRecreationalConservationColor - Map ICLUS V3 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusRecreationalConservationColor( double percent,
                                          double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 0, 255, 80 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusTimberColor - Map ICLUS V4 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusTimberColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 0, 78, 14 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusGrazingColor - Map ICLUS V5 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusGrazingColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 127, 255, 159 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusPastureColor - Map ICLUS V6 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusPastureColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 185, 255, 144 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusCroplandColor - Map ICLUS V7 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusCroplandColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 99, 127, 7 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusMiningBarrenColor - Map ICLUS V8 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusMiningBarrenColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 168, 168, 168 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusParksGolfCoursesColor - Map ICLUS V9 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusParksGolfCoursesColor( double percent, double unused1,double maximum) {
  const Color result = categoryColor( percent, maximum, 0, 129, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusExurbanLowDensityColor - Map ICLUS V10 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusExurbanLowDensityColor( double percent, double unused1,double maximum) {
  const Color result = categoryColor( percent, maximum, 255, 127, 255 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusExurbanHighDensityColor - Map ICLUS V11 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusExurbanHighDensityColor(double percent, double unused1,double maximum) {
  const Color result = categoryColor( percent, maximum, 255, 0, 255 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusSuburbanColor - Map ICLUS V12 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusSuburbanColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 255, 186, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusUrbanLowDensityColor - Map ICLUS V13 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusUrbanLowDensityColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 255, 127, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusUrbanHighDensityColor - Map ICLUS V14 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusUrbanHighDensityColor( double percent, double unused1, double maximum) {
  const Color result = categoryColor( percent, maximum, 255, 0, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusCommercialColor - Map ICLUS V15 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusCommercialColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 255, 255, 129 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusIndustrialColor - Map ICLUS V16 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusIndustrialColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 255, 255, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusInstitutionalColor - Map ICLUS V17 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusInstitutionalColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 255, 205, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: iclusTransportationColor - Map ICLUS V18 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/ICLUS_legend.png.
******************************************************************************/

Color iclusTransportationColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 255, 64, 0 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdOpenWaterColor - Map NLCD V11 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdOpenWaterColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 79, 107, 161 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdSnowIceColor - Map NLCD V12 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdSnowIceColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 232, 239, 252 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdWDevelopedOpenSpaceColor - Map NLCD V21 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdWDevelopedOpenSpaceColor(double percent, double unused1,double maximum) {
  const Color result = categoryColor( percent, maximum, 224, 205, 206 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdDevelopedLowIntensityColor - Map NLCD V22 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdDevelopedLowIntensityColor( double percent,
                                      double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 216, 151, 130 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdDevelopedMediumIntensityColor - Map NLCD V23 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdDevelopedMediumIntensityColor( double percent,
                                         double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 233, 0, 23 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdDevelopedHighIntensityColor - Map NLCD V24 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdDevelopedHighIntensityColor( double percent,
                                       double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 165, 0, 13 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdBarrenColor - Map NLCD V31 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdBarrenColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 178, 175, 164 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdDeciduousForestColor - Map NLCD V41 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdDeciduousForestColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 110, 172, 103 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdEvergreenForestColor - Map NLCD V42 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdEvergreenForestColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 36, 103, 51 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdMixedForestColor - Map NLCD V43 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdMixedForestColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 189, 206, 148 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdShrubScrubColor - Map NLCD V52 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdShrubScrubColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 206, 188, 132 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdGrasslandColor - Map NLCD V71 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdGrasslandColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 214, 216, 159 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdPastureHayColor - Map NLCD V81 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdPastureHayColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 218, 219, 68 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdCropsColor - Map NLCD V82 to color.
INPUTS:  double percent  Value to map to a color.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdCropsColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 170, 114, 46 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdWetlandForestColor - Map NLCD V90 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdWetlandForestColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 190, 215, 236 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nlcdWetlandEmergentColor - Map NLCD V95 to color.
INPUTS:  double percent  Value to map to a color.
         double maximum  Maximum percentage.
RETURNS: Color mapped value.
NOTES:   Must match data/land_use/NLCD_2006_legend.png.
******************************************************************************/

Color nlcdWetlandEmergentColor( double percent, double unused1, double maximum ) {
  const Color result = categoryColor( percent, maximum, 119, 163, 192 );
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: stratTypeTextColor - Map text value to color for STRAT_TYPE data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color stratTypeTextColor( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( ! strcmp( text, "well-mixed" ) ) {
    result.b = 1.0;
  } else if ( ! strcmp( text, "partial" ) ) {
    result.g = 1.0;
  } else if ( ! strcmp( text, "strong" ) ) {
    result.g = 1.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: stratMethTextColor - Map text value to color for STRAT_METH data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color stratMethTextColor( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( ! strcmp( text, "Froude" ) ) {
    result.b = 1.0;
  } else if ( ! strcmp( text, "QV" ) ) {
    result.r = 1.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: ecoRegionTextColor - Map text value to color for ECO_REGION data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color ecoRegionTextColor( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( ! strcmp( text, "CRN" ) ) {
    result.b = 1.0;
  } else if ( ! strcmp( text, "FLN" ) ) {
    result.g = 0.5;
    result.b = 1.0;
  } else if ( ! strcmp( text, "GM" ) ) {
    result.r = 0.5;
    result.g = 0.5;
  } else if ( ! strcmp( text, "GOM" ) ) {
    result.r = 1.0;
    result.g = 0.5;
  } else if ( ! strcmp( text, "NCA" ) ) {
    result.g = 1.0;
    result.b = 1.0;
  } else if ( ! strcmp( text, "PT" ) ) {
    result.g = 1.0;
  } else if ( ! strcmp( text, "SCB" ) ) {
    result.g = 1.0;
    result.r = 1.0;
  } else if ( ! strcmp( text, "VCF" ) ) {
    result.r = 1.0;
  } else if ( ! strcmp( text, "VGN" ) ) {
    result.r = 1.0;
    result.b = 1.0;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nhdCodeTextColor - Map text NHD_CODE value to color for stream water
         temperature data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color nhdCodeTextColor( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( *text == '{' ) { /* Hexadecimal {0078DCE9-0495-870D-FB59EA815619}. */
    const unsigned long minimum = 7920873; /* Only 1st part: 0078DCE9. */
    const unsigned long maximum = 4290652823U; /* FFBE2A97. */
    unsigned long value = 0;
    sscanf( text + 1, "%lx", &value ); /* Skip opening brace, parse hex. */
    result = dataColor( value, minimum, maximum );
  } else { /* Integer value: */
    const int minimum = 145056585;
    const int maximum = 145061063;
    const int value = atoi( text );
    result = dataColor( value, minimum, maximum );
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: nhdPlusIdTextColor - Map text NHDPLUS_ID value to color for stream
         water temperature data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color nhdPlusIdTextColor( const char* text ) {
  PRE02( text, *text );
  const double minimum = 55000100000009.0;
  const double maximum = 55000100973830.0;
  const double value = atof( text );
  const Color result = dataColor( value, minimum, maximum );

fprintf( stderr, "------ %s %f in [%f, %f]\n", text, value, minimum, maximum );

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: temperatureRegimeTextColor - Map text value to color for stream water
         temperature regime data.
INPUTS:  const char* text  Text to map to a color.
RETURNS: Color mapped value.
******************************************************************************/

Color temperatureRegimeTextColor( const char* text ) {
  PRE02( text, *text );
  Color result = { 0.0, 0.0, 0.0 };

  if ( strstr( text, "cold" ) == text ) {
    result.b = 1.0;
  } else if ( strstr( text, "cool" ) == text ) { /* light blue. */
    result.g = 0.75;
    result.b = 1.0;
  } else if ( strstr( text, "warm" ) == text ) {
    result.r = 1.0;
  } else if ( strstr( text, "no s" ) == text ) {

    /* "no s" (means no survival of salmon) = brown (dead fish). */

    result.r = 0.5;
    result.g = 0.25;
  }

  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: temperatureFlagLineStipple - Get line style mode for a stream water
         temperature regime_flag value.
INPUTS:  const char* text  Text to map to a stipple pattern.
RETURNS: int 0 = no draw, 1 = solid line (no stipple), else 0xAAAA = dotted.
******************************************************************************/

int temperatureFlagLineStipple( const char* text ) {
  PRE02( text, *text );
  int result = 1; /* Draw. */

  if ( strstr( text, "no flag") ) {
    result = 1; /* Solid line, no stippling. */
  } else if ( strstr( text, "flag" ) ) {
    result = 0xAAAA;
  } else if ( strstr( text, "drop" ) ) {
    result = 0;
  }

  POST0( IN4( result, 0, 1, 0xAAAA ) );
  return result;
}



/******************************************************************************
PURPOSE: streamTemperatureColor - Map stream temperature value to color
         for fish status.
INPUTS:  double value  Value (degrees C) to map to color.
RETURNS: Color mapped value.
******************************************************************************/

Color streamTemperatureColor( double value, double unused_1, double unused_2) {
  Color result = { 0.0, 0.0, 0.0 };

  if ( value >= 0.0 ) {

    if ( value < 12.0 ) { /* Too cold (ice) = gray-blue. */
      result.r = 0.5;
      result.g = 0.5;
      result.b = 0.75;
    } else if ( value <= 16.0 ) { /* cold = blue. */
      result.b = 1.0;
    } else if ( value <= 18.0 ) { /* cool = light blue. */
      result.g = 0.75;
      result.b = 1.0;
    } else { /* too warm = red. */
      result.r = 1.0;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: streamTemperatureCategoryColor - Map stream temperature category
         value to color for fish status.
INPUTS:  double value  Value to map to color. -99 = missing,
                                             1 = < 12 degrees C,
                                             2 = [12, 16)C,
                                             3 = [16, 18]C,
                                             4 = > 18C
RETURNS: Color mapped value.
******************************************************************************/

Color streamTemperatureCategoryColor( double value,
                                      double unused_1, double unused_2 ) {
  Color result = { 0.0, 0.0, 0.0 };
  const int ivalue = (int) value;

  if ( ivalue == 1 ) { /* 1 means < 12 degees C = Too cold = gray-blue. */
    result.r = 0.5;
    result.g = 0.5;
    result.b = 0.75;
  } else if ( ivalue == 2 ) { /* 2 = [12, 16) = cold = blue. */
    result.b = 1.0;
  } else if ( ivalue == 3 ) { /* 3 = [16, 18] = cool = light blue. */
    result.g = 0.75;
    result.b = 1.0;
  } else if ( ivalue == 4 ) { /* 4 = > 18C = too warm = red. */
    result.r = 1.0;
  }

  return result;
}



/******************************************************************************
PURPOSE: hruStippling - Get stippling mode for given WSM dataset and hruId
INPUTS:  const char* const name  Dataset name.
         const int hruId         HRU_ID value.
RETURNS: int 1 for striped, 2 for zigzag, else 0.
******************************************************************************/

int hruStippling( const char* const name, const int hruId ) {
  PRE02( name, *name );
  int result = 0;

  if ( hruId > 0 ) {

    if ( strstr( name, "_blackstone" ) ) {
      result = IN_RANGE( hruId, 1, 9 );
    } else if ( strstr( name, "_charles" ) ) {

      if ( IN_RANGE( hruId, 10, 99 ) ) {
        const int lastDigit = hruId % 10;

        if ( lastDigit == 2 ) {
          result = 1;
        } else if ( lastDigit == 3 ) {
          result = 2;
        }
      }
    } else if ( strstr( name, "_farmington" ) ) {

      if ( IN_RANGE( hruId, 100, 999 ) ) {
        const int firstDigit = hruId / 100;

        if ( firstDigit <= 2 ) {
          result = firstDigit;
        }
      }
    } else if ( strstr( name, "_ipswich" ) ) {

      if ( IN_RANGE( hruId, 8, 12 ) ) {
        result = 1;
      } else if ( IN_RANGE( hruId, 14, 15 ) ) {
        result = 2;
      }
    } else if ( strstr( name, "_pawcatuck" ) ) {
      result = IN_RANGE( hruId, 11, 15 );
    } else if ( strstr( name, "_sudbury" ) ) {
      result = IN_RANGE( hruId, 30, 49 );
    } else if ( strstr( name, "_taunton" ) ) {
      result = IN_RANGE( hruId, 10, 16 );
    }
  }

  POST0( IN_RANGE( result, 0, 2 ) );
  return result;
}



/******************************************************************************
PURPOSE: daltonize - Apply Dalton Algorithm to adjust an image for better
         differentiation by people with one of three forms of dichromacy.
 INPUTS:  const int dichromacy  Type of 2-cone dificiency to correct for.
                1 = Protanope:   greatly reduced red   in 1% of males.
                2 = Deuteranope: greatly reduced green in 1% of males.
                3 = Tritanope:   greatly reduced blue  in 0.003% of people.
         const size_t length  Length of rgb[].
         float rgb[ length ]  Normalized RGB values [0.0, 1.0] to adjust.
                              Stored as triplets rgb, rgb, rgb, ...
OUTPUTS: float rgb[ length ]  dichromacy-adjusted RGB values to better
                              differentiate hues.
NOTES:   See: http://www.daltonize.org
         This routine is based on the (less efficient) JavaScript routine:
https://galactic.ink/labs/Color-Vision/Javascript/Color.Vision.Daltonize.js
******************************************************************************/

void daltonize( const int dichromacy, const size_t length, float rgb[] ) {

  /* size_t */ long long index = 0; /* OpenMP requires signed loop index. */

  PRE07( IN4( dichromacy, 1, 2 , 3 ),
         length,
         length % 3 == 0,
         rgb,
         IN_RANGE( rgb[ 0 ], 0.0, 1.0 ),
         IN_RANGE( rgb[ length / 2 ], 0.0, 1.0 ),
         IN_RANGE( rgb[ length - 1 ], 0.0, 1.0 ) );

#pragma omp parallel for

  for ( index = 0; index < length; index += 3 ) {
    const size_t index_1 = index + 1;
    const size_t index_2 = index + 2;

    /* Get original image r, g, b: */

    const double r = rgb[ index ];
    const double g = rgb[ index_1 ];
    const double b = rgb[ index_2 ];

    /* Check to avoid useless computation on unhued pixels: */

    if ( r != g || r != b || g != b ) {

      /* Compute lms (long/medium/short-wave) from rgb: */

      const double l = 17.8824   * r + 43.5161  * g + 4.11935 * b;
      const double m = 3.45565   * r + 27.1554  * g + 3.86714 * b;
      const double s = 0.0299566 * r + 0.184309 * g + 1.46709 * b;

      /* Simulate color-blind/impaired lms: */

      double l_i = l;
      double m_i = m;
      double s_i = s;

      switch ( dichromacy ) {
      case 1: /* Protanope: Reds are greatly reduced. */
        l_i = 2.02344 * m + -2.52581 * s;
        break;
      case 2: /* Deuteranope: Greens are greatly reduced. */
        m_i = 0.494207 * l + 1.24827 * s;
        break;
      default: /* Tritanope: Blues are greatly reduced.*/
        s_i = -0.395913 * l + 0.801109 * m;
        break;
      }

      /* Compute impaired rgb from impaired lms: */

      const double r_i =
        0.0809444479 * l_i + -0.130504409 * m_i + 0.116721066 * s_i;
      const double g_i =
        -0.0102485335 * l_i + 0.0540193266 * m_i + -0.113614708 * s_i;
      const double b_i =
        -0.000365296938 * l_i + -0.00412161469 * m_i + 0.693511405 * s_i;

      /* Compute difference between original rgb and imparied rgb: */

      const double r_delta = r - r_i;
      const double g_delta = g - g_i;
      const double b_delta = b - b_i;

      /* Shift color towards visible spectrum (note no change to red): */

      const double scaled_r_delta = 0.7 * r_delta;
      const double g_shifted = scaled_r_delta + g_delta;
      const double b_shifted = scaled_r_delta + b_delta;

      /* Add shift to original rgb: */

      double g_adjusted = g + g_shifted;
      double b_adjusted = b + b_shifted;

      /* Clamp to valid range: */

      if ( g_adjusted < 0.0 ) {
        g_adjusted = 0.0;
      } else if ( g_adjusted > 1.0 ) {
        g_adjusted = 1.0;
      }

      if ( b_adjusted < 0.0 ) {
        b_adjusted = 0.0;
      } else if ( b_adjusted > 1.0 ) {
        b_adjusted = 1.0;
      }

      /* Store modified color: */

      rgb[ index_1 ] = g_adjusted;
      rgb[ index_2 ] = b_adjusted;
    }
  }
}



/******************************************************************************
PURPOSE: accumulate - Sum input into output.
INPUTS:  const float* input  Data to accumulate.
         int count           Number of values in input/output.
OUTPUTS: float* output       Accumulated data.
******************************************************************************/

void accumulate( const float* input, float* output, int count ) {
  PRE03( input, output, count > 0 );
  int index = 0;

/* TEMP HACK until pthreads/OpenMP BUG is fixed. #pragma omp parallel for */

  for ( index = 0; index < count; ++index ) {
    output[ index ] += input[ index ];
  }
}



/******************************************************************************
PURPOSE: scaledMaximum - Scale array values and return maximum.
INPUTS:  float array[ count ]  Values to scale.
         int   count           Number of values.
         double scale           Scale factor.
         double threshold       If value < threshold then it is not scaled.
OUTPUTS: float array[ count ]  Scaled values.
RETURNS: double maximum scaled value.
******************************************************************************/

double scaledMaximum( float array[], int count, double scale, double threshold ) {
  double result = 0.0;
  int index = 0;

  if ( array[ 0 ] > threshold ) {
    array[ 0 ] *= scale;
  }

  result = array[ 0 ];

  for ( index = 1; index < count; ++index ) {

    if ( array[ index ] > threshold ) {
      array[ index ] *= scale;
    }

    if ( array[ index ] > result ) {
      result = array[ index ];
    }
  }

  DEBUG( fprintf( stderr, "scaledMaximum = %f\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: unsignedShortMaximum - Maximum of unsigned short values.
INPUTS:  const short array[ count ]  Values to check.
         const size_t count          Number of values.
RETURNS: int maximum value.
******************************************************************************/

int unsignedShortMaximum( const unsigned short array[], const size_t count ) {
  PRE02( array, count );
  int index = 0;
  int result = array[ 0 ];

  for ( index = 1; index < count; ++index ) {

    if ( array[ index ] > result ) {
      result = array[ index ];
    }
  }

  DEBUG( fprintf( stderr, "shortMaximum = %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: charMaximum - Maximum of char values.
INPUTS:  const char array[ count ]  Values to check.
         const size_t count          Number of values.
RETURNS: int maximum value.
******************************************************************************/

int charMaximum( const char array[], const size_t count ) {
  PRE02( array, count );
  int index = 0;
  int result = array[ 0 ];

  for ( index = 1; index < count; ++index ) {

    if ( array[ index ] > result ) {
      result = array[ index ];
    }
  }

  DEBUG( fprintf( stderr, "charMaximum = %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: swapCharDataRows - Swap char data rows in-place.
INPUTS:  char array[ timesteps * rows * columns ] Array of 1-byte values.
         const int timesteps                      Number of timesteps in array.
         const int rows                           Number of rows in array.
         const int columns                        Number of columns in array.
OUTPUTS: char array[ timesteps * rows * columns ] Array with data rows swapped.
******************************************************************************/

void swapCharDataRows( char array[], const int timesteps,
                       const int rows, const int columns ) {
  PRE05( array, timesteps > 0, rows > 0, columns > 0,
         timesteps * rows * columns > 0 );
  char* a = array;
  const int timestepOffset = rows * columns;
  int timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep, a += timestepOffset ) {
    int lowerRowIndex = 0;
    int upperRowIndex = rows - 1;

    for ( ; lowerRowIndex < upperRowIndex; ++lowerRowIndex, --upperRowIndex ) {
      char* const lowerRowData = a + lowerRowIndex * columns;
      char* const upperRowData = a + upperRowIndex * columns;
      int column = 0;

      for ( column = 0; column < columns; ++column ) {
        char temp = lowerRowData[ column ];
        lowerRowData[ column ] = upperRowData[ column ];
        upperRowData[ column ] = temp;
      }
    }
  }
}



/******************************************************************************
PURPOSE: swapUnsignedShortDataRows - Swap unsigned short data rows in-place.
INPUTS:  unsigned short array[ timesteps * rows * columns ] 2-byte values.
         const int timesteps                      Number of timesteps in array.
         const int rows                           Number of rows in array.
         const int columns                        Number of columns in array.
OUTPUTS: unsigned short array[ timesteps * rows * columns ] Array rows swapped.
******************************************************************************/

void swapUnsignedShortDataRows( unsigned short array[], const int timesteps,
                                const int rows, const int columns ) {
  PRE05( array, timesteps > 0, rows > 0, columns > 0,
         timesteps * rows * columns > 0 );
  unsigned short* a = array;
  const int timestepOffset = rows * columns;
  int timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep, a += timestepOffset ) {
    int lowerRowIndex = 0;
    int upperRowIndex = rows - 1;

    for ( ; lowerRowIndex < upperRowIndex; ++lowerRowIndex, --upperRowIndex ) {
      unsigned short* const lowerRowData = a + lowerRowIndex * columns;
      unsigned short* const upperRowData = a + upperRowIndex * columns;
      int column = 0;

      for ( column = 0; column < columns; ++column ) {
        unsigned short temp = lowerRowData[ column ];
        lowerRowData[ column ] = upperRowData[ column ];
        upperRowData[ column ] = temp;
      }
    }
  }
}



/******************************************************************************
PURPOSE: swapFloatDataRows - Swap float data rows in-place.
INPUTS:  float array[ timesteps * rows * columns ] Array of 4-byte values.
         const int timesteps                      Number of timesteps in array.
         const int rows                           Number of rows in array.
         const int columns                        Number of columns in array.
OUTPUTS: float array[ timesteps * rows * columns ] Array w/ data rows swapped.
******************************************************************************/

void swapFloatDataRows( float array[], const int timesteps,
                        const int rows, const int columns ) {
  PRE05( array, timesteps > 0, rows > 0, columns > 0,
         timesteps * rows * columns > 0 );
  float* a = array;
  const int timestepOffset = rows * columns;
  int timestep = 0;

  for ( timestep = 0; timestep < timesteps; ++timestep, a += timestepOffset ) {
    int lowerRowIndex = 0;
    int upperRowIndex = rows - 1;

    for ( ; lowerRowIndex < upperRowIndex; ++lowerRowIndex, --upperRowIndex ) {
      float* const lowerRowData = a + lowerRowIndex * columns;
      float* const upperRowData = a + upperRowIndex * columns;
      int column = 0;

      for ( column = 0; column < columns; ++column ) {
        float temp = lowerRowData[ column ];
        lowerRowData[ column ] = upperRowData[ column ];
        upperRowData[ column ] = temp;
      }
    }
  }
}



/******************************************************************************
PURPOSE: expandBytesToFloats - Expand byte array into float array in-place.
INPUTS:  float array[ count ]  Array of 1-byte values.
         int count             Number of values in array.
OUTPUTS: float array[ count ]  Array of 4-byte values.
******************************************************************************/

void expandBytesToFloats( float array[], int count ) {
  PRE02( array, count > 0 );
  const signed char* input = (const signed char*) array;
  float* output = array + count;
  int counter = count;
  input += count;

  while ( counter-- ) {
    *--output = *--input;
  }
}



/******************************************************************************
PURPOSE: widthInMeters - Width of lon-lat window in meters.
INPUTS:  const Bounds bounds  Lon-lat bounds of window.
RETURNS: double width of (center line of) window in meters.
******************************************************************************/

double widthInMeters( const Bounds bounds ) {
  PRE0( isValidBounds( bounds ) );
  const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const double latitudeMean     = 0.5 * ( latitudeMinimum + latitudeMaximum );
  const double longitudeRange   = longitudeMaximum - longitudeMinimum;
  const double meanEarthRadius = 6371009.0; /* meters. */
  const double degreesPerMeter =
    1.0 / ( 2.0 * M_PI * meanEarthRadius / 360.0 * cos(radians(latitudeMean)));
  const double result = longitudeRange / degreesPerMeter;
  POST0( result > 0.0 );
  return result;
}



/******************************************************************************
PURPOSE: computeArrowVectorCoordinates - Compute coordinates of arrow vector.
INPUTS:  const double longitude           Longitude of base of arrow.
         const double latitude            Latitude of base of arrow.
         const double x                   X-component of arrow.
         const double y                   Y-component of arrow.
         const double degreesPerPixel     Degrees longitude per viewport pixel.
         const double pixelsPerUnitLength Viewport pixels per unit vector.
OUTPUTS: double* const point0X            X-coordinate of arrow base vertex.
         double* const point0Y            Y-coordinate of arrow base vertex.
         double* const point1X            X-coordinate of arrow tip vertex.
         double* const point1Y            Y-coordinate of arrow tip vertex.
         double* const point2X            X-coordinate of arrow head top vertex
         double* const point2Y            Y-coordinate of arrow head top vertex
         double* const point3X            X-coordinate of arrow head bot vertex
         double* const point3Y            Y-coordinate of arrow head bot vertex
******************************************************************************/

void computeArrowVectorCoordinates( const double longitude,
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
                                    double* const point3Y ) {

  PRE013( isValidLongitudeLatitude( longitude, latitude ),
          ! isNan( x ),
          ! isNan( y ),
          IN_RANGE( degreesPerPixel, 1e-8, 1.0 ),
          IN_RANGE( pixelsPerUnitLength, 1.0, 100.0 ),
          point0X, point0Y, point1X, point1Y, point2X, point2Y,
          point3X, point3Y );

  const double arrowX = 0.8; /* Normalized X-coord of barb of unit vector. */
  const double arrowY = 0.2; /* Normalized Y-coord of barb of unit vector. */
  const double angle = atan2( y, x );
  const double r = hypot( x, y );
  const double scale = r * pixelsPerUnitLength * degreesPerPixel;
  const double scaledCosineAngle = scale * cos( angle );
  const double scaledSineAngle   = scale * sin( angle );
  const double arrowXScaledCosineAngle = arrowX * scaledCosineAngle;
  const double arrowXScaledSineAngle   = arrowX * scaledSineAngle;
  const double arrowYScaledCosineAngle = arrowY * scaledCosineAngle;
  const double arrowYScaledSineAngle   = arrowY * scaledSineAngle;
  const double arrowXScaledCosineAngleMinusArrowYScaledSineAngle =
    arrowXScaledCosineAngle - arrowYScaledSineAngle;
  const double arrowYScaledCosineAnglePlusArrowXScaledSineAngle =
    arrowYScaledCosineAngle + arrowXScaledSineAngle;
  const double arrowXScaledCosineAnglePlusArrowYScaledSineAngle =
    arrowXScaledCosineAngle + arrowYScaledSineAngle;
  *point0X = longitude;
  *point0Y = latitude;
  *point1X = longitude + scaledCosineAngle;
  *point1Y = latitude  + scaledSineAngle;
  *point2X = longitude + arrowXScaledCosineAngleMinusArrowYScaledSineAngle;
  *point2Y = latitude  + arrowYScaledCosineAnglePlusArrowXScaledSineAngle;
  *point3X = longitude + arrowXScaledCosineAnglePlusArrowYScaledSineAngle;
  *point3Y = latitude  - arrowYScaledCosineAngle + arrowXScaledSineAngle;

  POST05( *point0X == longitude, *point0Y == latitude,
          isValidLongitudeLatitude( *point1X, *point1Y ),
          isValidLongitudeLatitude( *point2X, *point2Y ),
          isValidLongitudeLatitude( *point3X, *point3Y ) );
}



/******************************************************************************
PURPOSE: wordCount - Count number of whitespace-separated printable words
         in a string.
INPUTS:  const char* string  String to count.
RETURNS: int number of words in string.
******************************************************************************/

int wordCount( const char* string ) {
  PRE0( string );
  const char* s = string;
  int result = 0;

  while ( *s ) {

    while ( isspace( *s ) ) {
      ++s;
    }

    result += isprint( *s ) != 0;

    while ( *s && ! isspace( *s ) ) {
      ++s;
    }
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: wordCount - Count lines in a string.
INPUTS:  const char* string  String to count.
RETURNS: int number of lines in string.
******************************************************************************/

int lineCount( const char* string ) {
  PRE0( string );
  int result = 0;
  const char* s = string;

  for ( s = string; *s; ++s ) {
    result += *s == '\n';
  }

  POST0( result >= 0 );
  return result;
}



/******************************************************************************
PURPOSE: lowercase - Convert a string to lowercase.
INPUTS:  char string[]  String to convert.
OUTPUTS: char string[]  Lowercase string.
******************************************************************************/

void lowercase( char string[] ) {
  PRE0( string );
  char* s = string;

  for ( s = string; *s; ++s ) {
    *s = tolower( *s );
  }
}



/******************************************************************************
PURPOSE: uppercase - Convert a string to uppercase.
INPUTS:  char string[]  String to convert.
OUTPUTS: char string[]  Uppercase string.
******************************************************************************/

void uppercase( char string[] ) {
  PRE0( string );
  char* s = string;

  for ( s = string; *s; ++s ) {
    *s = toupper( *s );
  }
}



/******************************************************************************
PURPOSE: changeChar - Change each instance of a character to another.
INPUTS:  char string[]  String to convert.
OUTPUTS: char string[]  Changed string.
******************************************************************************/

void changeChar( char string[], const char from, const char to ) {
  PRE0( string );
  char* s = string;

  while ( *s ) {

    if ( *s == from ) {
      *s = to;
    }

    ++s;
  }

  POST0( ! strchr( string, from ) );
}



/******************************************************************************
PURPOSE: countChar - Count each instance of a character in a string.
INPUTS:  const char string[]  String to check.
         const char ch        Character to check for.
RETURNS: size_t number of instances of ch in string.
******************************************************************************/

size_t countChar( const char string[], const char ch ) {
  PRE0( string );
  size_t result = 0;
  const char* s = string;

  while ( *s ) {

    if ( *s == ch ) {
      ++result;
    }

    ++s;
  }

  POST0( IMPLIES( result, strchr( string, ch ) ) );
  return result;
}



/******************************************************************************
PURPOSE: eraseChar - Erase all characters of string starting with the first
         occurrence of ch.
INPUTS:  char string[]  String to erase.
         const char ch  Character to locate.
OUTPUTS: char string[]  Erased string.
******************************************************************************/

void eraseChar( char string[], const char ch ) {
  PRE0( string );
  char* s = string;

  while ( ! IN3( *s, ch, '\0' ) ) {
    ++s;
  }

  while ( *s ) {
    *s = '\0';
    ++s;
  }

  POST0( ! strchr( string, ch ) );
}



/******************************************************************************
PURPOSE: shortenName - Convert a string to a shorter variable name.
INPUTS:  const LongName name[]  Name to shorten.
OUTPUTS: LongName shortenedName[]  Shortened name.
******************************************************************************/

void shortenName( const LongName name, LongName shortenedName ) {
  PRE03( name, shortenedName, name != shortenedName );
  LongName temp = "";
  memset( temp, 0, sizeof temp );
  strncpy( temp, name, sizeof temp / sizeof *temp - 1 );
  uppercase( temp );

  memset( shortenedName, 0, sizeof (LongName) / sizeof (char) );
  substituteWord( temp, "OXIDIZED", "OX", shortenedName );

  strncpy( temp, shortenedName, sizeof temp / sizeof *temp - 1 );
  memset( shortenedName, 0, sizeof (LongName) / sizeof (char) );
  substituteWord( temp, "REDUCED", "RD", shortenedName );

  strncpy( temp, shortenedName, sizeof temp / sizeof *temp - 1 );
  memset( shortenedName, 0, sizeof (LongName) / sizeof (char) );
  substituteWord( temp, "DRY_", "D_", shortenedName );

  strncpy( temp, shortenedName, sizeof temp / sizeof *temp - 1 );
  memset( shortenedName, 0, sizeof (LongName) / sizeof (char) );
  substituteWord( temp, "WET_", "W_", shortenedName );
}



/******************************************************************************
PURPOSE: eraseTrailingWhitespace - Erase any trailing whitespace in string.
INPUTS:  char string[]  String to erase.
OUTPUTS: char string[]  String without any trailing whitespace.
******************************************************************************/

void eraseTrailingWhitespace( char string[] ) {
  PRE0( string );
  size_t index = strlen( string );

  if ( index ) {
    int done = 0;

    do {
      --index;

      if ( isspace( string[ index ] ) ) {
        string[ index ] = '\0';
      } else {
        done = 1;
      }

    } while ( AND2( index, ! done ) );
  }
}



/******************************************************************************
PURPOSE: eraseLeadingWhitespace - Erase any leading whitespace in string.
INPUTS:  char string[]  String to erase.
OUTPUTS: char string[]  String without any leading whitespace.
******************************************************************************/

void eraseLeadingWhitespace( char string[] ) {
  PRE0( string );
  char* input  = string;
  char* output = string;

  while ( AND2( *input, isspace( *input ) ) ) {
    ++input;
  }

  while ( *output ) {
    *output = *input;
    ++output;

    if ( *input ) {
      ++input;
    }
  }
}



/******************************************************************************
PURPOSE: eraseDoubleQuotedCommas - Change commas to space in each double-quoted
         string.
INPUTS:  char string[]  String to convert.
OUTPUTS: char string[]  Converted string.
******************************************************************************/

void eraseDoubleQuotedCommas( char string[] ) {
  PRE0( string );
  char* s = string;
  int quoted = 0;

  for ( s = string; *s; ++s ) {

    if ( *s == '"' ) {
      quoted = ! quoted;
    } else if ( *s == ',' ) {

      if ( quoted ) {
        *s = ' ';
      }
    }
  }
}



/******************************************************************************
PURPOSE: substituteWord - Change occurrences of oldWord to newWord in input.
INPUTS:  const char* input    String to read.
         const char* oldWord  Word to search for.
         const char* newWord  Word to substitute.
OUTPUTS: char* output         Resulting string.
******************************************************************************/

void substituteWord( const char* input,
                     const char* oldWord,
                     const char* newWord,
                     char* output ) {
  PRE06( input, oldWord, *oldWord, newWord, output, input != output );
  const size_t length = strlen( oldWord );
  const char* i = input;
  char* o = output;

  do {
    const char* const p = strstr( i, oldWord ); /* Point to oldWord in input:*/

    if ( p ) {
      const char* n = newWord;

      while ( i < p ) { /* Copy input up to oldWord: */
        *o++ = *i++;
      }

      i += length; /* Skip oldWord. */

      while ( *n ) { /* Copy newWord: */
        *o++ = *n++;
      }
    } else { /* Copy remainder of input: */

      while ( *i ) {
        *o++ = *i++;
      }
    }
  } while ( *i );

  *o = '\0';
  POST0( ! strstr( output, oldWord ) );
}



/******************************************************************************
PURPOSE: parseInts - Parse a string into a sequence of values, each within
         a given range (possibly clamped to range) and return the values.
INPUTS:  char** string           Address of string to parse.
         const char* const delimiters  String of delimiter characters.
         const size_t count      Number of values expected to parse.
         const int range[][ 2 ]  MINIMUM and MAXIMUM allowed for each value.
         const int clamp         Clamp values to within range?
OUTPUTS: char** string           Updated pointer past parsed values.
         int values[ count ]     Parsed values within range.
RETURNS: int 1 if successful else 0.
******************************************************************************/

int parseInts( char** string,
               const char* const delimiters,
               const size_t count,
               const int range[][ 2 ],
               const int clamp,
               int values[] ) {

  int result = 1;
  size_t index = 0;

  for ( ; AND2( result, index < count ); ++index ) {
    const char* word = strtok_r( 0, delimiters, string );
    result = AND2( word, *word );

    if ( result ) {
      char* end2 = 0;
      int value = strtol( word, &end2, 10 );
      result = end2 != word;

      if ( result ) {
        const int minimum = range[ index ][ MINIMUM ];
        const int maximum = range[ index ][ MAXIMUM ];

        if ( clamp ) {

          if ( value < minimum ) {
            value = minimum;
          } else if ( value > maximum ) {
            value = maximum;
          }
        }

        result = IN_RANGE( value, minimum, maximum );

        if ( result ) {
          values[ index ] = value;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: parseDoubles - Parse a string into a sequence of values, each within
         a given range (possibly clamped to range) and return the values.
INPUTS:  char** string              Address of string to parse.
         const char* const delimiters  String of delimiter characters.
         const size_t count         Number of values expected to parse.
         const double range[][ 2 ]  MINIMUM and MAXIMUM allowed for each value.
         const int clamp            Clamp values to within range?
OUTPUTS: char** string              Updated pointer past parsed values.
         double values[ count ]     Parsed values within range.
RETURNS: int 1 if successful else 0.
******************************************************************************/

int parseDoubles( char** string,
                  const char* const delimiters,
                  const size_t count,
                  const double range[][ 2 ],
                  const int clamp,
                  double values[] ) {

  int result = 1;
  size_t index = 0;

  for ( ; AND2( result, index < count ); ++index ) {
    const char* word = strtok_r( 0, delimiters, string );
    result = AND2( word, *word );

    if ( result ) {
      char* end2 = 0;
      double value = strtod( word, &end2 );
      result = end2 != word;

      if ( result ) {
        const double minimum = range[ index ][ MINIMUM ];
        const double maximum = range[ index ][ MAXIMUM ];

        if ( clamp ) {

          if ( value < minimum ) {
            value = minimum;
          } else if ( value > maximum ) {
            value = maximum;
          }
        }

        result = IN_RANGE( value, minimum, maximum );

        if ( result ) {
          values[ index ] = value;
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: parseWords - Parse a sequence of words matching a given set of words
                      and store their 0-based indices.
INPUTS:  char** string              Address of string to parse.
         const char* const delimiters  String of delimiter characters.
         const size_t valueCount    Number of values expected to parse.
         const size_t wordCount     Number of words in matching string.
         const char* const words[ wordCount ]  Words to match.
OUTPUTS: char** string              Updated pointer past parsed values.
         int values[ count ]   0-based index into words[] of each parsed word.
RETURNS: int 1 if successful else 0.
******************************************************************************/

int parseWords( char** string,
                const char* const delimiters,
                const size_t valueCount,
                const size_t wordCount,
                const char* const words[],
                int values[] ) {

  int result = 1;
  size_t index = 0;

  for ( ; AND2( result, index < valueCount ); ++index ) {
    const char* word = strtok_r( 0, delimiters, string );
    result = AND2( word, *word );

    if ( result ) {
      const int wordIndex = indexOfString( word, words, wordCount );
      result = wordIndex != -1;

      if ( result ) {
        values[ index ] = wordIndex;
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: copyFile - Copy a file.
INPUTS:  const char* inputFileName   Pathed name of file to copy.
         const char* outputFileName  Created file.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int copyFile( const char* inputFileName, const char* outputFileName ) {
  PRE03( inputFileName, outputFileName, *outputFileName );
  int result = 1;

  if ( strcmp( inputFileName, outputFileName ) ) {
    FILE* inputFile = fopen( inputFileName, "rb" );
    FILE* outputFile = inputFile ? fopen( outputFileName, "wb" ) : 0;
    result = 0;

    if ( outputFile ) {
      enum { BUFFER_SIZE = 16384 };
      char buffer[ BUFFER_SIZE ] = "";
      int done = 0;

      do {
        const size_t bytes = fread( buffer, 1, BUFFER_SIZE, inputFile );
        result = 1;

        if ( bytes ) {
          result = fwrite( buffer, 1, bytes, outputFile ) == bytes;
          done = ! result;
        } else {
          done = 1;
        }

      } while ( ! done );
    }

    if ( inputFile ) {
      fclose( inputFile ), inputFile = 0;
    }

    if ( outputFile ) {
      fclose( outputFile ), outputFile = 0;
    }

    if ( ! result ) {
      failureMessage( "Failed to copy file %s\n", inputFileName );
    }
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: copyFileBytes - Copy a file's bytes to another file.
INPUTS:  FILE* inputFile             File to copy.
         const char* outputFileName  Created file.
         size_t bytes                Bytes to read and write.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int copyFileBytes( FILE* inputFile, const char* outputFileName, size_t bytes) {

  PRE04( inputFile, outputFileName, *outputFileName, bytes > 0 );

  int result = 0;
  FILE* outputFile = fopen( outputFileName, "wb" );

  if ( outputFile ) {
    enum { BUFFER_SIZE = 16384 };
    char buffer[ BUFFER_SIZE ] = "";
    size_t bytesRemaining = bytes;
    DEBUG2(fprintf(stderr, "copyFileBytes(%s, %lu)\n", outputFileName, bytes);)

    do {
      const size_t bytesToRead = MIN( bytesRemaining, BUFFER_SIZE );
      const size_t bytesRead = fread( buffer, 1, bytesToRead, inputFile );

      if ( bytesRead ) {
        const size_t bytesWritten = fwrite( buffer, 1, bytesRead, outputFile );
        result = bytesWritten == bytesRead;
        bytesRemaining -= bytesWritten;
        DEBUG2( fprintf( stderr, "  bytesWritten = %lu\n", bytesWritten ); )
      } else {
        result = 0;
      }

    } while ( AND2( result, bytesRemaining ) );

    fclose( outputFile ), outputFile = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "Failed to copy %lu bytes to output file %s ",
             bytes, outputFileName );
    perror( "because" );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: streamBytesToFile - Copy all of a stream's bytes to a file.
INPUTS:  FILE* stream                Stream to read.
         const char* outputFileName  Created file.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int streamBytesToFile( FILE* stream, const char* outputFileName ) {
  PRE03( stream, outputFileName, *outputFileName );
  int result = 0;
  FILE* outputFile = fopen( outputFileName, "wb" );

  if ( outputFile ) {
    enum { BUFFER_SIZE = 16384 };
    char buffer[ BUFFER_SIZE ] = "";
    size_t bytesRead = 0;

    do {
      bytesRead = fread( buffer, 1, BUFFER_SIZE, stream );

      if ( bytesRead ) {
        const size_t bytesWritten = fwrite( buffer, 1, bytesRead, outputFile );
        result = bytesWritten == bytesRead;
      }

    } while ( AND2( bytesRead, result ) );

    fclose( outputFile ), outputFile = 0;
  }

  if ( ! result ) {
    fprintf( stderr, "Failed to stream all bytes to output file %s ",
             outputFileName );
    perror( "because" );
  }

  return result;
}



/******************************************************************************
PURPOSE: streamFile - Stream bytes of a file to stdout.
INPUTS:  const char* fileName  Name of file to stream.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int streamFile( const char* fileName ) {
  PRE02( fileName, *fileName );
  int result = 0;
  FILE* file = fopen( fileName, "rb" );

  if ( file ) {
    enum { BUFFER_SIZE = 16384 };
    char buffer[ BUFFER_SIZE ] = "";
    size_t bytesRemaining = fileSize( fileName );

    do {
      const size_t bytesToStreamNow = MIN( bytesRemaining, BUFFER_SIZE );
      result = fread( buffer, bytesToStreamNow, 1, file ) == 1;
      result = AND2( result,
                     fwrite( buffer, bytesToStreamNow, 1, stdout ) == 1 );
      bytesRemaining -= bytesToStreamNow;
    } while ( AND2( result, bytesRemaining ) );

    fclose( file ), file = 0;
  }

  DEBUG( fprintf( stderr, "streamFile() returning %d\n", result ); )
  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: fileExists - Does the named file exist?
INPUTS:  const char* name  Name of file to examine.
RETURNS: int 1 if the file exists, else 0.
******************************************************************************/

int fileExists( const char* name ) {
  PRE02( name, *name );
  struct stat unused;
  const int result = stat( name, &unused ) != -1;
  return result;
}



/******************************************************************************
PURPOSE: fileSize - Determine size of named file.
INPUTS:  const char* name  Name of file to examine.
RETURNS: size_t Size, in bytes, of named file, else 0 if failed.
******************************************************************************/

size_t fileSize( const char* name ) {
  PRE02( name, *name );
  size_t result = 0;
  struct stat buf;

  if ( stat( name, &buf ) == -1 ) {
    fprintf( stderr, "\nFailed to determine size of file '%s'.\n", name );
  } else {
    result = buf.st_size;
  }

  return result;
}



/******************************************************************************
PURPOSE: readFile - Read named file into memory and return it as an allocated
         string.
INPUTS:  const char* name  Name of file to examine.
OUTPUTS: size_t* length    Length of string.
         size_t* lines     Number of lines (\n) in string.
RETURNS: char* string contents of file
         (with any '\r' characters converted to ' '),
         else 0 if failed and failureMessage() is called.
******************************************************************************/

char* readFile( const char* name, size_t* length, size_t* lines ) {
  PRE03( name, length, lines );
  char* result = 0;
  *length = fileSize( name ) / sizeof (char);
  *lines = 0;

  if ( *length > 0 ) {
    result = (char*) malloc( sizeof (char) * ( *length + 1 ) );

    if ( ! result ) {
      fprintf( stderr, "\nFailed to allocate %lu bytes to read file %s.\n",
               *length + 1, name );
      *length = 0;
    } else {
      FILE* file = fopen( name, "rb" );

      if ( file ) {
        const size_t itemsRead =
          fread( result, *length * sizeof (char), 1, file );

        if ( itemsRead != 1 ) {
          fprintf( stderr, "\nFailed to read entire file '%s'.\n", name );
          free( result ), result = 0;
          *length = 0;
        } else {
          result[ *length ] = '\0'; /* Terminate string. */
          *lines = controlMToSpace( result );
        }

        fclose( file );
        file = 0;
      }
    }
  }

  POST0( IMPLIES_ELSE( result,
                       AND2( *length > 0, *lines > 0 ),
                       AND2( *length == 0, *lines == 0 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: controlMToSpace - Convert any '\r' characters to ' '.
INPUTS:  char* string  String to filter.
OUTPUTS: char* string  Filtered string.
RETURNS: size_t number of lines (\n) in string.
******************************************************************************/

size_t controlMToSpace( char* string ) {
  PRE0( string );
  size_t result = 0;
  char* s = 0;

  for ( s = string; *s; ++s ) {

    if ( *s == '\r' ) {
      *s = ' ';
    } else if ( *s == '\n' ) {
      ++result;
    }
  }

  POST0( ! strchr( string, '\r' ) );
  return result;
}



/******************************************************************************
PURPOSE: copyFileLine - Copy a line from an input file to an output file.
INPUTS:  FILE* input   Input file to read one line from.
OUTPUTS: FILE* output  Output file to write one line to.
RETURNS: size_t number of characters written.
******************************************************************************/

size_t copyFileLine( FILE* input, FILE* output ) {
  PRE02( input, output );
  size_t result = 0;
  char c = 0;

  do {
    c = fgetc( input );

    if ( c != EOF ) {
      c = fputc( c, output );
      result += c != EOF;
    }
  } while ( AND2( c != '\n', c != EOF ) );

  return result;
}



/******************************************************************************
PURPOSE: skipFileLine - Read and discard a line from a file.
INPUTS:  FILE* input   Input file to read one line from.
RETURNS: size_t number of characters read.
******************************************************************************/

size_t skipFileLine( FILE* input ) {
  PRE0( input );
  size_t result = 0;
  char c = 0;

  do {
    c = fgetc( input );
    result += c != EOF;
  } while ( AND2( c != '\n', c != EOF ) );

  return result;
}



/******************************************************************************
PURPOSE: removeAllFiles - Remove all files in a directory.
INPUTS:  const char* name  Name of directory to remove files from.
******************************************************************************/

void removeAllFiles( const char* name ) {
  PRE02( name, *name );
  DIR* directory = opendir( name );

  if ( directory ) {
    struct dirent* entry = 0;

    do {
      entry = readdir( directory );

      if ( AND3( entry, entry->d_name, entry->d_name[ 0 ] != '.' ) ) {
        enum { LENGTH = 255 };
        char pathedName[ LENGTH + 1 ] = "";
        memset( pathedName, 0, sizeof pathedName );
        snprintf( pathedName, sizeof pathedName / sizeof *pathedName,
                  "%s/%s", name, entry->d_name );
        unlink( pathedName );
      }

    } while ( entry );

    closedir( directory ), directory = 0;
  }
}



/******************************************************************************
PURPOSE: removeMatchedFiles - Remove matched files in a directory.
INPUTS:  const char* const directoryName Name of directory to remove files from
         const char* const startsWith Optional matching file name start.
         const char* const endsWith   Optional matching file name end.
******************************************************************************/

void removeMatchedFiles( const char* const directoryName,
                         const char* const startsWith,
                         const char* const endsWith ) {

  PRE04( directoryName,
         *directoryName,
         IMPLIES( startsWith, *startsWith ),
         IMPLIES( endsWith, *endsWith ) );

  forEachFile( directoryName, startsWith, endsWith, removeFile0, 0 );
}



/******************************************************************************
PURPOSE: copyFiles - Copies matched files in a directory to another directory.
INPUTS:  const char* const fromDirectory  Name of directory to copy files from.
         const char* const toDirectory    Name of directory to copy files to.
         const char* const startsWith     Optional matching file name start.
         const char* const endsWith       Optional matching file name end.
******************************************************************************/

void copyFiles( const char* const fromDirectory,
                const char* const toDirectory,
                const char* const startsWith,
                const char* const endsWith ) {

  PRE07( fromDirectory, *fromDirectory, toDirectory, *toDirectory,
         strcmp( fromDirectory, toDirectory ),
         IMPLIES( startsWith, *startsWith ),
         IMPLIES( endsWith, *endsWith ) );

  DIR* directory = opendir( fromDirectory );

  if ( directory ) {
    struct dirent* entry = 0;

    do {
      entry = readdir( directory );

      if ( AND3( entry, entry->d_name, entry->d_name[ 0 ] != '.' ) ) {
        const char* const fileName = entry->d_name;
        int matched = 1;

        if ( startsWith ) {
          const char* const found = strstr( fileName, startsWith );
          matched = found == fileName;
        }

        if ( matched ) {

          if ( endsWith ) {
            const size_t endLength = strlen( endsWith );
            const size_t fileNameLength = strlen( fileName );

            if ( fileNameLength >= endLength ) {
              const char* const endOfFileName =
                fileName + fileNameLength - endLength;
              matched = strcmp( endOfFileName, endsWith ) == 0;
            }
          }
        }

        if ( matched ) {
          enum { LENGTH = 255 };
          char fromFileName[ LENGTH + 1 ] = "";
          char toFileName[ LENGTH + 1 ] = "";
          memset( fromFileName, 0, sizeof fromFileName );
          memset( toFileName, 0, sizeof toFileName );
          snprintf( fromFileName, sizeof fromFileName / sizeof *fromFileName,
                    "%s/%s", fromDirectory, fileName );
          snprintf( toFileName, sizeof toFileName / sizeof *toFileName,
                    "%s/%s", toDirectory, fileName );
          copyFile( fromFileName, toFileName );
        }
      }

    } while ( entry );

    closedir( directory ), directory = 0;
  }
}



/******************************************************************************
PURPOSE: forEachFile - Calls callback for each matched file in a directory.
INPUTS:  const char* const directory  Name of directory to check for files.
         const char* const startsWith     Optional matching file name start.
         const char* const endsWith       Optional matching file name end.
         ForEachFileCallback callback     Called with each matched file name.
         void* callbackData               Second argument to callback.
******************************************************************************/

void forEachFile( const char* const directory,
                  const char* const startsWith,
                  const char* const endsWith,
                  ForEachFileCallback callback,
                  void* callbackData ) {

  PRE05( directory,
         *directory,
         IMPLIES( startsWith, *startsWith ),
         IMPLIES( endsWith, *endsWith ),
         callback );

  DIR* dir = opendir( directory );

  if ( dir ) {
    struct dirent* entry = 0;

    do {
      entry = readdir( dir );

      if ( AND3( entry, entry->d_name, entry->d_name[ 0 ] != '.' ) ) {
        const char* const fileName = entry->d_name;
        int matched = 1;

        if ( startsWith ) {
          const char* const found = strstr( fileName, startsWith );
          matched = found == fileName;
        }

        if ( matched ) {

          if ( endsWith ) {
            const size_t endLength = strlen( endsWith );
            const size_t fileNameLength = strlen( fileName );

            if ( fileNameLength >= endLength ) {
              const char* const endOfFileName =
                fileName + fileNameLength - endLength;
              matched = strcmp( endOfFileName, endsWith ) == 0;
            }
          }
        }

        if ( matched ) {
          callback( fileName, callbackData );
        }
      }

    } while ( entry );

    closedir( dir ), dir = 0;
  }
}



/******************************************************************************
PURPOSE: directoryExists - Does directory exist?
INPUTS:  const char* name  Name of directory to remove files from.
******************************************************************************/

int directoryExists( const char* name ) {
  PRE02( name, *name );
  DIR* directory = opendir( name );
  const int result = directory != 0;

  if ( directory ) {
    closedir( directory ), directory = 0;
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: directoryListing - List files in a directory that end with one of
         the given extensions and are dated today.
INPUTS:  const char* directory   Name of directory.
         const char* extensions  Optional: file extensions to filter by.
                                 E.g., "shx shp txt".
         int size                Size of buffer[].
OUTPUTS: char buffer[ size ]     Result containing formatted columns:
                                 bytes hh:mm name
                                 sorted by file name.
NOTES:   Omits file names that begin with dot (.).
******************************************************************************/

void directoryListing( const char* directory, const char* extensions, int size,
                       char buffer[] ) {
  PRE04( directory, *directory, size > 0, buffer );
  int remaining = size - 1;
  enum { NAME_LENGTH = 1024 };
  char pathedFileName[ NAME_LENGTH + 1 ] = "";
  DIR* dir = opendir( directory );
  memset( pathedFileName, 0, sizeof pathedFileName );
  memset( buffer, 0, size );

  if ( dir ) {
    const time_t clock = time( 0 );
    const struct tm* const timeInfo = localtime( &clock );

    if ( timeInfo ) {
      const int todayYear  = timeInfo->tm_year + 1900;
      const int todayMonth = timeInfo->tm_mon + 1;
      const int todayDay   = timeInfo->tm_mday;
      const struct dirent* entry = 0;

      while ( ( entry = readdir( dir ) ) != 0 ) {
        const char* const fileName = entry->d_name;

        if ( AND2( fileName, *fileName != '.' ) ) {
          int matched = extensions == 0; /* If no filter then all match. */

          if ( extensions ) { /* Check that file extension matches filter: */
            const char* const dot = strrchr( fileName, '.' );

            if ( dot ) {
              const char* const extension = dot + 1;

              if ( *extension ) { /* File name did not end in '.': */
                matched = matchesWord( extension, extensions );
              }
            }
          }

          if ( matched ) {
            struct stat info;
            memset( &info, 0, sizeof info );
            snprintf( pathedFileName, NAME_LENGTH, "%s/%s",
                      directory, fileName );

            if ( ! stat( pathedFileName, &info ) ) {
              const time_t seconds = info.st_mtime;
              const struct tm* const fileTimeInfo = localtime( &seconds );
              const int year   = fileTimeInfo->tm_year + 1900;
              const int month  = fileTimeInfo->tm_mon + 1;
              const int day    = fileTimeInfo->tm_mday;

              if ( AND3( year  == todayYear,
                         month == todayMonth,
                         day   == todayDay ) ) {
                const int hh = fileTimeInfo->tm_hour;
                const int mm = fileTimeInfo->tm_min;
                const size_t fileSize = info.st_size;
                enum { LINE_LENGTH = 255 };
                char line[ LINE_LENGTH + 1 ] = "";
                int lineLength = 0;
                memset( line, 0, sizeof line );
                snprintf( line, LINE_LENGTH, "%9lu %02d:%02d %s\n",
                          fileSize, hh, mm, fileName );
                line[ LINE_LENGTH ] = '\0';
                lineLength = strlen( line );

                if ( remaining > lineLength ) {
                  strncat( buffer, line, remaining );
                  remaining -= lineLength;
                }
              }
            }
          }
        }
      }
    }

    closedir( dir ), dir = 0;
  }
}



/******************************************************************************
PURPOSE: homeDirectory - Name of user's home directory.
RETURNS: const char* Name of user's home directory. Examples:
         "/Users/plessel"
         "C:\Documents and Settings\tplessel"
         "."
******************************************************************************/

const char* homeDirectory( void ) {
  static char result[ 256 ] = "";

  if ( result[ 0 ] == '\0' ) {
#if defined( _WIN32 ) || defined( __CYGWIN__ )
    const char* const homeDrive0 = getenv( "HOMEDRIVE" );
    const char* const homePath0  = getenv( "HOMEPATH" );
    const char* const homeDrive = homeDrive0 ? homeDrive0 : "";
    const char* const homePath = homePath0 ? homePath0 : ".";
    memset( result, 0, sizeof result );
    snprintf( result, sizeof result / sizeof *result,
              "%s%s", homeDrive, homePath );
#else
    const char* const home0 = getenv( "HOME" );
    const char* const home = home0 ? home0 : ".";
    memset( result, 0, sizeof result );
    strncpy( result, home, sizeof result / sizeof *result - 1 );
#endif
  }

  DEBUG( fprintf( stderr, "homeDirectory = '%s'\n", result ); )
  POST02( *result, result[ strlen( result ) - 1 ] != '/' );
  return result;
}



/******************************************************************************
PURPOSE: sortUniqFile - Sort lines of a text file and filter-out consecutive
         matched lines.
INPUTS:  const char* name          Name of file to sort.
         const size_t headerLines  Number of lines at the beginning to keep.
RETURNS: size_t number of lines rewritten to file.
******************************************************************************/

size_t sortUniqFile( const char* name, const size_t headerLines ) {
  PRE02( name, *name );
  size_t result = 0;
  size_t length = 0;
  size_t lines = 0;
  char* buffer = readFile( name, &length, &lines );

  if ( buffer ) {
    char** fileLines = malloc( lines * sizeof (char*) );

    if ( fileLines ) {
      char* s = buffer;
      size_t line = 1;
      fileLines[ 0 ] = buffer;

      for ( s = buffer; *s; ++s ) {

        if ( *s == '\n' ) {
          *s = '\0';

          if ( line < lines ) {
            fileLines[ line ] = s + 1;
            ++line;
          }
        }
      }

      if ( line == lines ) {
        size_t uniqLines = headerLines;
        qsort( fileLines + headerLines, lines - headerLines, sizeof (char*),
               stringCompare );

        /* uniq the data lines: */

        for ( line = uniqLines + 1; line < lines; ++line ) {

          if ( strcmp( fileLines[ line ], fileLines[ uniqLines ] ) ) {
            ++uniqLines;
            fileLines[ uniqLines ] = fileLines[ line ];
          }
        }

        if ( uniqLines < lines ) {
          ++uniqLines;
        }

        /* Write sorted unique lines to named file: */

        {
          FILE* file = fopen( name, "w" );

          if ( file ) {
            int ok = 1;

            for ( line = 0; AND2( ok, line < uniqLines ); ++line ) {
              ok = AND2( fputs( fileLines[ line ], file ) != EOF,
                         fputc( '\n', file ) == '\n' );
              result += ok;
            }

            fclose( file ), file = 0;
          }
        }
      }

      free( fileLines ), fileLines = 0;
    }

    free( buffer ), buffer = 0;
  }

  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



static int minimumItemi( const int array[], int count ) {
  int result = array[ 0 ];
  int index = 0;

  for ( index = 1; index < count; ++index ) {
    const int value = array[ index ];

    if ( value < result ) {
      result = value;
    }
  }

  DEBUG( fprintf( stderr, "minimumItemi = %d\n", result ); )
  return result;
}



static int sum( const int array[], int count ) {
  int result = 0;
  int index = 0;

  for ( index = 0; index < count; ++index ) {
    result += array[ index ];
  }

  DEBUG( fprintf( stderr, "sum = %d\n", result ); )
  return result;
}



static int inBounds( const float vertices[], int count,
                     const double bounds[ 2 ][ 2 ] ) {
  PRE03( vertices, count > 0, isValidBounds( bounds ) );
  const double tolerance = 1e-3;
  const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const float* longitudes = vertices;
  const float* latitudes  = vertices + 1;
  int vertex = 0;
  int result = 1;

  for ( vertex = 0; vertex < count && result; ++vertex,
        longitudes += 2, latitudes += 2 ) {
    result =
      IN_RANGE( *longitudes,
                longitudeMinimum - tolerance,
                longitudeMaximum + tolerance ) &&
      IN_RANGE( *latitudes,
                latitudeMinimum - tolerance,
                latitudeMaximum + tolerance );

#ifdef DEBUGGING
    if ( ! result ) {
      fprintf( stderr, "Invalid vertex: index %d (%f, %f)\n",
               vertex, *longitudes, *latitudes ) ;
    }
#endif
  }

  DEBUG( fprintf( stderr, "inBounds = %d\n", result ); )
  return result;
}



static int inBoundsDouble( const double vertices[], int count,
                           const double bounds[ 2 ][ 2 ] ) {
  PRE03( vertices, count > 0, isValidBounds( bounds ) );
  const double tolerance = 1e-3;
  const double longitudeMinimum = bounds[ LONGITUDE ][ MINIMUM ];
  const double longitudeMaximum = bounds[ LONGITUDE ][ MAXIMUM ];
  const double latitudeMinimum  = bounds[ LATITUDE  ][ MINIMUM ];
  const double latitudeMaximum  = bounds[ LATITUDE  ][ MAXIMUM ];
  const double* longitudes = vertices;
  const double* latitudes  = vertices + 1;
  int vertex = 0;
  int result = 1;

  for ( vertex = 0; vertex < count && result; ++vertex,
        longitudes += 2, latitudes += 2 ) {
    result =
      IN_RANGE( *longitudes,
                longitudeMinimum - tolerance,
                longitudeMaximum + tolerance ) &&
      IN_RANGE( *latitudes,
                latitudeMinimum - tolerance,
                latitudeMaximum + tolerance );

#ifdef DEBUGGING
    if ( ! result ) {
      fprintf( stderr, "Invalid vertex: index %d (%lf, %lf)\n",
               vertex, *longitudes, *latitudes ) ;
    }
#endif
  }

  DEBUG( fprintf( stderr, "inBoundsDouble = %d\n", result ); )
  return result;
}



/******************************************************************************
PURPOSE: soilColor - Map normalized t value to color for soil data.
INPUTS:  double t  Normalized t value.
RETURNS: Color mapped value.
******************************************************************************/

static Color soilColor( double t ) {
  Color result = { 0.0, 0.0, 0.0 };
  const double oneOver255 = 1.0 / 255.0;

  if ( t < 0.0 || isNan( t ) ) {
    t = 0.0;
  } else if ( t > 1.0 ) {
    t = 1.0;
  }

  result.r = 64.0 + ( 255.0 - 64.0 ) * t;
  result.g = 31.0 + ( 221.0 - 31.0 ) * t;
  result.b = 167.0 * t;
  result.r *= oneOver255;
  result.g *= oneOver255;
  result.b *= oneOver255;
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: soilColor4 - Map normalized t value to discrete 4 color for soil data.
INPUTS:  double t  Normalized t value.
RETURNS: Color mapped value.
******************************************************************************/

static Color soilColor4( double t ) {
  Color result = { 0.0, 0.0, 0.0 };
  const double oneOver255 = 1.0 / 255.0;

  if ( t < 0.0 || isNan( t ) ) {
    t = 0.0;
  } else if ( t > 1.0 ) {
    t = 1.0;
  }

  if ( t < 1.0 / 4.0 ) {
    result.r = 64.0;
    result.g = 31.0;
  } else if ( t < 2.0 / 4.0 ) {
    result.r = 127.0;
    result.g = 64.0;
    result.b = 31.0;
  } else if ( t < 3.0 / 4.0 ) {
    result.r = 191.0;
    result.g = 127.0;
    result.b = 64.0;
  } else {
    result.r = 255.0;
    result.g = 221.0;
    result.b = 167.0;
  }

  result.r *= oneOver255;
  result.g *= oneOver255;
  result.b *= oneOver255;
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: soilColor7 - Map normalized t value to discrete 7 color for soil data.
INPUTS:  double t  Normalized t value.
RETURNS: Color mapped value.
******************************************************************************/

static Color soilColor7( double t ) {
  Color result = { 0.0, 0.0, 0.0 };
  const double oneOver255 = 1.0 / 255.0;

  if ( t < 0.0 || isNan( t ) ) {
    t = 0.0;
  } else if ( t > 1.0 ) {
    t = 1.0;
  }

  if ( t < 1.0 / 7.0 ) {
    result.r = 64.0;
    result.g = 31.0;
  } else if ( t < 2.0 / 7.0 ) {
    result.r = 96.0;
    result.g = 48.0;
    result.b = 16.0;
  } else if ( t < 3.0 / 7.0 ) {
    result.r = 127.0;
    result.g = 64.0;
    result.b = 31.0;
  } else if ( t < 4.0 / 7.0 ) {
    result.r = 159.0;
    result.g = 96.0;
    result.b = 48.0;
  } else if ( t < 5.0 / 7.0 ) {
    result.r = 191.0;
    result.g = 127.0;
    result.b = 64.0;
  } else if ( t < 6.0 / 7.0 ) {
    result.r = 233.0;
    result.g = 174.0;
    result.b = 116.0;
  } else {
    result.r = 255.0;
    result.g = 221.0;
    result.b = 167.0;
  }

  result.r *= oneOver255;
  result.g *= oneOver255;
  result.b *= oneOver255;
  POST0( isValidColor( result ) );
  return result;
}



/******************************************************************************
PURPOSE: categoryColor - Map percent of category to color.
INPUTS:  const double percent  Percent of category to map to a color.
         const double maximum  Maximum percent of category to map.
         const int    red      Red   component of color. [0, 255].
         const int    green    Green component of color. [0, 255].
         const int    blue     Blue  component of color. [0, 255].
RETURNS: Color mapped value.
******************************************************************************/

static Color categoryColor( const double percent, const double maximum,
                            const int red, const int green, const int blue ) {
  PRE03( IN_RANGE( red,   0, 255 ),
         IN_RANGE( green, 0, 255 ),
         IN_RANGE( blue,  0, 255 ) );
  Color result = { CATEGORY_MINIMUM, CATEGORY_MINIMUM, CATEGORY_MINIMUM };

  if ( AND2( IN_RANGE( percent, 0.0, 100.0 ),
             IN_RANGE( maximum, 0.0, 100.0 ) ) ) {
    const double t = percent / maximum;
    const double normalizedValue =
      IN_RANGE( t, 0.0, 1.0 ) ? t
      : t > 1.0 ? 1.0
      : t < 0.0 ? 0.0
      : 0.0; /* Could be NaN from division. If so, map to 0. */
    const double factor = CATEGORY_SCALE_FACTOR( normalizedValue );

    if ( factor > CATEGORY_MINIMUM ) {
      const double oneOver255 = 1.0 / 255.0;
      const double redValue   = red   * oneOver255 * factor;
      const double greenValue = green * oneOver255 * factor;
      const double blueValue  = blue  * oneOver255 * factor;

      if ( redValue > CATEGORY_MINIMUM ) {
        result.r = redValue;
      }

      if ( greenValue > CATEGORY_MINIMUM ) {
        result.g = greenValue;
      }

      if ( blueValue > CATEGORY_MINIMUM ) {
        result.b = blueValue;
      }
    }
  }

  POST0( isValidColor( result ) );
  return result;
}


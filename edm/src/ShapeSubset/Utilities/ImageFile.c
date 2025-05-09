
/******************************************************************************
PURPOSE: ImageFile.c - Read subsets of *.bin image files.
NOTES:
HISTORY: 2010-12-12 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For FILE, stderr, fprintf(), fopen(), perror(). */
#include <string.h> /* For  memset(). */
#include <limits.h> /* For INT_MAX. */
#include <stdlib.h> /* For malloc(), free(). */

#include <Assertions.h> /* For PRE(), IMPLIES_ELSE(). */
#include <Utilities.h>  /* For Bounds, MINIMUM LONGITUDE, isValidBounds(). */
#include <ImageFile.h>  /* For public interface. */

/* Value clamped to range [low, high]. */

#define CLAMPED_TO_RANGE( value, low, high ) \
((value) < (low) ? (low) : (value) > (high) ? (high) : (value))

/*=========================== FORWARD DECLARATIONS ==========================*/

enum { BYTES_PER_PIXEL = 3 }; /* RGB. */

static int readImageFileHeader( FILE* file,
                                int* width, int* height, Bounds corners );

static void subsetIndices( int isWidth, const double domain[ 2 ],
                           double corners[ 2 ],
                           int count, int indices[ 2 ] );

/*============================= PUBLIC FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: readImageFile - Read a subset of an image.
INPUTS:  const char* fileName  Name of image file (.bin).
         const Bounds clip     Optional Lon-lat bounds to clip image to or 0.
OUTPUTS: int* width            Subsetted image width in pixels.
         int* height           Subsetted image height in pixels.
         Bounds domain         Lon-lat bounds of unsubsetted image.
         Bounds corners        Lon-lat bounds of subsetted image.
RETURNS: unsigned char* array of RGB pixels of subsetted image.
NOTES:   free() result when no longer needed, e.g., when used as a texture.
******************************************************************************/

unsigned char* readImageFile( const char* fileName, const Bounds clip,
                              int* width, int* height,
                              Bounds domain, Bounds corners ) {

  PRE07( fileName, *fileName,
         IMPLIES( clip, isValidBounds( clip ) ),
         width, height, domain, corners );

  unsigned char* result = 0;
  int ok = 0;
  FILE* file = fopen( fileName, "rb" );
  *width = *height = 0;
  memset( domain, 0, sizeof (Bounds) );
  memset( corners, 0, sizeof (Bounds) );

  DEBUG( fprintf( stderr, "readImageFile %s [%lf, %lf] [%lf, %lf]\n",
                  fileName,
                  clip ? clip[ LONGITUDE ][ MINIMUM ] : -180.0,
                  clip ? clip[ LONGITUDE ][ MAXIMUM ] : 180.0,
                  clip ? clip[ LATITUDE  ][ MINIMUM ] : -90.0,
                  clip ? clip[ LATITUDE  ][ MAXIMUM ] : 90.0 ); )

  if ( ! file ) {
    fprintf( stderr, "\nFailed to open image file %s.\n", fileName );
  } else if ( ! readImageFileHeader( file, width, height, domain ) ) {
    fprintf( stderr, "\nFailed to read valid header of image file %s.\n",
             fileName );
  } else {
    enum { BYTES_PER_PIXEL = 3 };
    int subsetRows[ 2 ] = { 0, *height - 1 };
    int subsetColumns[ 2 ] = { 0, *width - 1 };
    memcpy( corners, domain, sizeof (Bounds) );

    if ( clip ) {
      subsetIndices( 0, clip[ LATITUDE ],
                     corners[ LATITUDE ], *height, subsetRows );

      if ( subsetRows[ MINIMUM ] != -1 ) { /* Adjust for upside-down rows. */
        const int swapTemp = subsetRows[ MINIMUM ];
        subsetRows[ MINIMUM ] = subsetRows[ MAXIMUM ];
        subsetRows[ MAXIMUM ] = swapTemp;
        subsetRows[ MINIMUM ] = *height - 1 - subsetRows[ MINIMUM ];
        subsetRows[ MAXIMUM ] = *height - 1 - subsetRows[ MAXIMUM ];
        subsetIndices( 1, clip[ LONGITUDE ], corners[ LONGITUDE ],
                       *width, subsetColumns );
      }
    }

    DEBUG( fprintf( stderr,
                    "subsetRows = [%d %d] of %d "
                    "subsetColumns = [%d %d] of %d\n",
                    subsetRows[ MINIMUM ], subsetRows[ MAXIMUM ], *height,
                    subsetColumns[ MINIMUM ], subsetColumns[MAXIMUM], *width);)

    if ( AND2( subsetRows[ MINIMUM ] != -1, subsetColumns[ MINIMUM ] != -1 )) {
      const int rows = subsetRows[ MAXIMUM ] - subsetRows[ MINIMUM ] + 1;
      const int columns = subsetColumns[MAXIMUM] - subsetColumns[MINIMUM] + 1;
      const int bytes = rows * columns * BYTES_PER_PIXEL;
      result = malloc( bytes );

      if ( ! result ) {
        fprintf( stderr,
                "\nFailed to allocate enough memory (%d bytes) "
                "to read image file %s.\n", bytes, fileName );
      } else {
        const int startColumnBytes = subsetColumns[MINIMUM] * BYTES_PER_PIXEL;
        const int seekBytes =
        subsetRows[ MINIMUM ] * *width * BYTES_PER_PIXEL + startColumnBytes;
        memset( result, 0, bytes );

        if ( AND2( seekBytes, fseek( file, seekBytes, SEEK_CUR ) ) ) {
          fprintf( stderr,
                  "\nFailed to seek to %d bytes into image file %s.\n",
                  seekBytes, fileName );
        } else {
          const int rowBytes = columns * BYTES_PER_PIXEL;
          const int readRows = columns == *width ? rows : 1;
          const int readBytes = readRows * rowBytes;
          const int skipBytes = ( *width - columns ) * BYTES_PER_PIXEL;
          int bytesRead = 0;
          unsigned char* input = result;
          ok = 1;

          DEBUG( fprintf( stderr,
                          "rows = %d, columns = %d, bytes = %d, "
                          "seekBytes = %d, rowBytes = %d, readRows = %d, "
                          "readBytes = %d, skipBytes = %d\n",
                          rows, columns, bytes,
                          seekBytes, rowBytes, readRows,
                          readBytes, skipBytes ); )

          while ( AND2( ok, bytesRead != bytes ) ) {
            ok = fread( input, readBytes, 1, file ) == 1;
            ok = AND2( ok,
                       IMPLIES( skipBytes,
                                fseek( file, skipBytes, SEEK_CUR ) == 0 ) );
            input += readBytes;
            bytesRead += readBytes;
          }

          if ( ! ok ) {
            fprintf( stderr,
                    "\nFailed to read %d bytes of image data from file %s.\n",
                    bytes, fileName );
          } else {
            *width = columns;
            *height = rows;
          }
        }
      }
    }
  }

  if ( file ) {
    fclose( file );
    file = 0;
  }

  if ( ! ok ) {

    if ( result ) {
      free( result );
      result = 0;
    }

    *width = *height = 0;
    memset( domain, 0, sizeof (Bounds) );
    memset( corners, 0, sizeof (Bounds) );
  }

  DEBUG( fprintf( stderr,
                  "readImageFile output: %d x %d [%lf %lf] [%lf, %lf]\n",
                  *width, *height,
                  corners[ LONGITUDE ][ MINIMUM ],
                  corners[ LONGITUDE ][ MAXIMUM ],
                  corners[ LATITUDE  ][ MINIMUM ],
                  corners[ LATITUDE  ][ MAXIMUM ] ); )

  DEBUG( if ( result ) fprintf( stderr,
                  "result %p = [%u %u %u] ... [%u %u %u]\n",
                  result,
                  result[ 0 ], result[ 1 ], result[ 2 ],
                  result[ *width * *height - 3 ],
                  result[ *width * *height - 2 ],
                  result[ *width * *height - 1 ] ); )

  POST0( IMPLIES_ELSE( result,
                       AND4( *width > 0, *height > 0,
                             isValidBounds( (const double (*)[2]) domain ),
                             isValidBounds( (const double (*)[2]) corners ) ),
                       IS_ZERO10( *width, *height,
                                  domain[0][0], domain[0][1],
                                  domain[1][0], domain[1][1],
                                  corners[0][0], corners[0][1],
                                  corners[1][0], corners[1][1] ) ) );

  return result;
}



/*============================= PRIVATE FUNCTIONS ===========================*/



/******************************************************************************
PURPOSE: readImageFileHeader - Read image file (.bin) header dimensions/range.
INPUTS:  FILE* file      Image file to read.
OUTPUTS: int* width      Width of image in pixels.
         int* height     Height of image in pixels.
         Bounds corners  Lon-lat extent of image.
RETURNS: int 1 if successful, else 0 and a failure message is printed to stderr
******************************************************************************/

static int readImageFileHeader( FILE* file,
                                int* width, int* height, Bounds corners ) {
  int result = 0;
  *width = *height = 0;
  memset( corners, 0, sizeof (Bounds) );

  result =
  AND6( fscanf( file,
               "%*[^\n]\n%*[^\n]\n%d %d %lf %lf %lf %lf\n%*[^\n]%*c",
               width, height,
               &corners[ LONGITUDE ][ MINIMUM ],
               &corners[ LONGITUDE ][ MAXIMUM ],
               &corners[ LATITUDE  ][ MINIMUM ],
               &corners[ LATITUDE  ][ MAXIMUM ] ) == 6,
       *width > 0, *height > 0, *width * *height > 0,
       *width % 2 == 0,
       isValidBounds( (const double (*)[2]) corners ) );

  DEBUG( fprintf( stderr,
                 "readImageFileHeader result = %d: "
                 "%d x %d [%lf %lf] [%lf, %lf]\n",
                 result, *width, *height,
                 corners[ LONGITUDE ][ MINIMUM ],
                 corners[ LONGITUDE ][ MAXIMUM ],
                 corners[ LATITUDE  ][ MINIMUM ],
                 corners[ LATITUDE  ][ MAXIMUM ] ); )

  if ( ! result ) {
    fprintf( stderr,
            "\nFailed to read valid header of image file: "
            "%d x %d [%lf %lf] [%lf, %lf]\n",
            *width, *height,
            corners[ LONGITUDE ][ MINIMUM ],
            corners[ LONGITUDE ][ MAXIMUM ],
            corners[ LATITUDE  ][ MINIMUM ],
            corners[ LATITUDE  ][ MAXIMUM ] );
    *width = *height = 0;
    memset( corners, 0, sizeof (Bounds) );
  }

  POST02( IS_BOOL( result ),
         IMPLIES_ELSE( result,
                       AND4( *width > 0, *height > 0, *width * *height > 0,
                             isValidBounds( (const double (*)[2]) corners ) ),
                       IS_ZERO6( *width, *height,
                                 corners[0][0], corners[0][1],
                                 corners[1][0], corners[1][1] ) ) );
  return result;
}



/******************************************************************************
PURPOSE: subsetIndices - Subset image range based on clip and compute indices.
 INPUTS: int isWidth              1 if subsetting width, else 0.
         const double clip[ 2 ]   clip [ MINIMUM MAXIMUM ] lon-lat bounds.
         double range[ 2 ]        range[ MINIMUM MAXIMUM ] of unclipped image.
         int count                Pixels in image along subset dimension.
OUTPUTS: double range[ 2 ]        range[ MINIMUM MAXIMUM ] of clipped image.
         int indices[ 2 ]         indices[ MINIMUM MAXIMUM ] of clipped image.
******************************************************************************/

static void subsetIndices( int isWidth, const double clip[ 2 ],
                           double range[ 2 ], int count, int indices[ 2 ] ) {

  PRE07( IS_BOOL( isWidth ),
         clip,
         range,
         IMPLIES_ELSE( isWidth,
                       AND4( IN_RANGE( clip[ MINIMUM ], -180.0, 180.0 ),
                             IN_RANGE( clip[ MAXIMUM ], clip[MINIMUM], 180.0),
                             IN_RANGE( range[ MINIMUM ], -180.0, 180.0 ),
                             IN_RANGE( range[MAXIMUM], range[MINIMUM], 180.0)),
                       AND4( IN_RANGE( clip[ MINIMUM ], -90.0, 90.0 ),
                             IN_RANGE( clip[ MAXIMUM ], clip[MINIMUM], 90.0),
                             IN_RANGE( range[ MINIMUM ], -90.0, 90.0 ),
                             IN_RANGE( range[MAXIMUM], range[MINIMUM], 90.0))),
          count > 0,
          count % 2 == 0,
          indices );

  const double clipMinimum = clip[ MINIMUM ];
  const double clipMaximum = clip[ MAXIMUM ];
  const double rangeMinimum = range[ MINIMUM ];
  const double rangeMaximum = range[ MAXIMUM ];
  const double clipRange = clipMaximum - clipMinimum;
  const double rangeRange = rangeMaximum - rangeMinimum;
  const double TOO_SMALL = 1e-6;

  indices[ MINIMUM ] = indices[ MAXIMUM ] = -1;

  DEBUG( fprintf( stderr,
                  "subsetIndices input:  "
                  "range[%lf %lf] clip(%lf %lf) count = %d, isWidth = %d\n",
                  range[ MINIMUM ], range[ MAXIMUM ],
                  clip[ MINIMUM ], clip[ MAXIMUM ], count, isWidth ); )

  if ( ! OR4( clipRange < TOO_SMALL, rangeRange < TOO_SMALL,
              rangeMaximum < clipMinimum, rangeMinimum > clipMaximum ) ) {
    const double scale = 1.0 / rangeRange;
    const double rangeIncrement = rangeRange / count;

    DEBUG( fprintf( stderr, "scale = %lf, rangeIncrement = %lf\n",
                    scale, rangeIncrement ); )

    if ( rangeMinimum > clipMinimum ) {
      indices[ MINIMUM ] = 0;
      DEBUG( fprintf( stderr, "set indices[ MINIMUM ] = 0\n" ); )
    } else {
      const double interpolation = ( clipMinimum - rangeMinimum ) * scale;
      DEBUG( fprintf( stderr, "interpolation = %lf\n", interpolation ); )
      indices[ MINIMUM ] = (int) ( interpolation * count - 0.5 );
      DEBUG( fprintf( stderr, "indices[ MINIMUM ] = %d\n", indices[ MINIMUM]);)
      indices[ MINIMUM ] = CLAMPED_TO_RANGE( indices[ MINIMUM ], 0, count - 1);
      DEBUG( fprintf( stderr, "indices[ MINIMUM ] = %d\n", indices[ MINIMUM]);)
      range[ MINIMUM ] += indices[ MINIMUM ] * rangeIncrement;
      DEBUG( fprintf( stderr,
                      "indices[ MINIMUM ] = %d, range[ MINIMUM ] = %lf\n",
                      indices[ MINIMUM ], range[ MINIMUM ] ); )
    }

    if ( rangeMaximum < clipMaximum ) {
      indices[ MAXIMUM ] = count - 1;
      DEBUG( fprintf( stderr, "set indices[ MAXIMUM ] = count - 1\n" ); )
    } else {
      const double interpolation = ( clipMaximum - rangeMinimum ) * scale;
      DEBUG( fprintf( stderr, "interpolation = %lf\n", interpolation ); )
      indices[ MAXIMUM ] = (int) ( interpolation * count + 0.5 );
      DEBUG( fprintf( stderr, "indices[ MAXIMUM ] = %d\n", indices[ MAXIMUM]);)
      indices[ MAXIMUM ] =
        CLAMPED_TO_RANGE( indices[ MAXIMUM ], indices[ MINIMUM ], count - 1 );
      DEBUG( fprintf( stderr,
                      "indices[ MAXIMUM ] = %d\n", indices[ MAXIMUM ] ); )
    }

    /* Ensure subset width is a multiple of 4: */

    if ( isWidth ) {
      int subsetCount = indices[ MAXIMUM ] - indices[ MINIMUM ] + 1;

      if ( subsetCount < 4 ) {
        indices[ 0 ] = indices[ 1 ] = -1;
      } else {
        const int remainder = ( subsetCount * BYTES_PER_PIXEL ) % 4;
        DEBUG( fprintf( stderr, "subsetCount = %d\n", subsetCount ); )
        DEBUG( fprintf( stderr, "remainder = %d\n", remainder ); )

        if ( remainder ) {
          subsetCount -= 4 - remainder;
          DEBUG( fprintf( stderr, "subsetCount = %d\n", subsetCount ); )
          indices[ MAXIMUM ] = indices[ MINIMUM ] + subsetCount - 1;
        }
      }

      DEBUG( fprintf( stderr, "indices[ MAXIMUM ] = %d\n", indices[MAXIMUM]);)
    }

    CHECK( IMPLIES( AND2( isWidth, indices[ 0 ] != -1 ),
                    ( indices[ MAXIMUM ] - indices[ MINIMUM] + 1 ) % 4 == 0 ));

    range[ MAXIMUM ] = range[ MINIMUM ] +
      ( indices[ MAXIMUM ] - indices[ MINIMUM ] ) * rangeIncrement;

    DEBUG( fprintf( stderr, "range[ MAXIMUM ] = %lf\n", range[ MAXIMUM ] ); )

    if ( ! AND5( IN_RANGE( range[ MINIMUM ],
                           clip[ MINIMUM ] - 0.1, clip[ MAXIMUM ] + 0.1 ),
                 IN_RANGE( range[ MAXIMUM ],
                           range[ MINIMUM ], clip[ MAXIMUM ] + 0.1 ),
                 IN_RANGE( indices[ MINIMUM ], 0, count - 1 ),
                 IN_RANGE( indices[ MAXIMUM ],
                           indices[ MINIMUM ], count - 1 ),
                 IMPLIES( isWidth,
                          (indices[MAXIMUM] - indices[MINIMUM] + 1) % 4==0))) {

      DEBUG( fprintf( stderr, "set indices[] = -1\n" ); )
      DEBUG( fprintf( stderr, "because:\n"
                      "IN_RANGE( range[ MINIMUM ],"
                      "clip[ MINIMUM ] - 0.1, clip[ MAXIMUM ] + 0.1 ) = %d, "
                      "IN_RANGE( range[ MAXIMUM ],"
                      "range[ MINIMUM ], clip[ MAXIMUM ] + 0.1 ) = %d, "
                      "IN_RANGE( indices[ MINIMUM ], 0, count - 1 ) = %d, "
                      "IN_RANGE( indices[ MAXIMUM ],"
                      "indices[ MINIMUM ], count - 1 ) = %d, "
                      "IMPLIES( isWidth,"
                         " (indices[MAXIMUM] - indices[MINIMUM] + 1) %% 4==0) "
                      "= %d\n",
                     IN_RANGE( range[ MINIMUM ],
                               clip[ MINIMUM ] - 0.1, clip[ MAXIMUM ] + 0.1 ),
                     IN_RANGE( range[ MAXIMUM ],
                               range[ MINIMUM ], clip[ MAXIMUM ] + 0.1 ),
                     IN_RANGE( indices[ MINIMUM ], 0, count - 1 ),
                     IN_RANGE( indices[ MAXIMUM ],
                               indices[ MINIMUM ], count - 1 ),
                     IMPLIES( isWidth,
                             (indices[MAXIMUM] - indices[MINIMUM] +1)%4==0));)

      indices[ MINIMUM ] = indices[ MAXIMUM ] = -1;
    }
  }

  DEBUG( fprintf( stderr,
                  "subsetIndices output: "
                  "range[%lf %lf] indices[%d %d]\n",
                  range[ MINIMUM ], range[ MAXIMUM ],
                  indices[ MINIMUM ], indices[ MAXIMUM ] ); )

  POST0( IMPLIES( ! AND2( indices[ MINIMUM ] == -1, indices[ MAXIMUM ] == -1),
                  AND5( IN_RANGE( range[ MINIMUM ],
                                  clip[ MINIMUM ] - 0.1, clip[ MAXIMUM] + 0.1),
                        IN_RANGE( range[ MAXIMUM ],
                                  range[ MINIMUM ], clip[ MAXIMUM ] + 0.1 ),
                        IN_RANGE( indices[ MINIMUM ], 0, count - 1 ),
                        IN_RANGE( indices[ MAXIMUM ],
                                  indices[ MINIMUM ], count - 1 ),
                        IMPLIES( isWidth,
                                 (indices[MAXIMUM] - indices[MINIMUM] + 1) % 4
                                   == 0 ) ) ) );
}



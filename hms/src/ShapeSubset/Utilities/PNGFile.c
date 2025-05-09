/******************************************************************************
PURPOSE: PNGFile.c - Simple PNG image file read/write routines.
HISTORY: 2009/05/26 plessel.todd@epa.gov Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdio.h>  /* For FILE, stderr, fprintf(), fread(). */
#include <stdlib.h> /* For malloc(), free(). */
#include <string.h> /* For memset(). */

#include <png.h> /* For png_*. */

#include "PNGFile.h" /* To match declarations. */

/*============================ PUBLIC FUNCTIONS =============================*/



/******************************************************************************
PURPOSE: readPNGFile - Read a PNG image stream.
INPUTS:  FILE* input  Stream to read.
         int width    Image width in pixels to match.
         int height   Image height in pixels to match.
         int reverseImageRows Reverse image rows in rgb output?
OUTPUTS: unsigned char rgb[ height * width * 3 ]  RGB pixel values.
RETURNS: int 1 if successful, else 0 and a message is printed to stderr.
******************************************************************************/

int readPNGFile( FILE* input, int width, int height, int reverseImageRows,
                 unsigned char rgb[] ) {
  int result = 0;

  if ( ! ( width > 0 && height > 0 ) ) {
    fprintf( stderr, "\nUnsupported PNG file dimensions.\n" );
  } else {
    png_byte header[ 8 ];

    if ( fread( header, 8, 1, input ) != 1 || png_sig_cmp( header, 0, 8) != 0) {
      fprintf( stderr, "\nInvalid PNG header.\n" );
    } else {
      png_structp png = png_create_read_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0);

      if ( ! png ) {
        fprintf( stderr, "\nFailed to create PNG read structure.\n" );
      } else {
        png_infop info = png_create_info_struct( png );

        if ( ! info ) {
          fprintf( stderr, "\nFailed to create PNG info structure.\n" );
        } else {
          png_uint_32 imageWidth = 0;
          png_uint_32 imageHeight = 0;
          int bitDepth = 0;
          int colorType = 0;
          int interlacing = 0;
          int unused_ = 0;

          png_init_io( png, input );
          png_set_sig_bytes( png, 8 );
          png_read_info( png, info );

          png_get_IHDR( png, info, &imageWidth, &imageHeight,
                        &bitDepth, &colorType, &interlacing,
                        &unused_, &unused_ );

          if ( ! ( imageWidth == (png_uint_32) width &&
                   imageHeight == (png_uint_32) height &&
                   bitDepth == 8 &&
                   ( colorType == PNG_COLOR_TYPE_PALETTE ||
                     colorType == PNG_COLOR_TYPE_RGB ||
                     colorType == PNG_COLOR_TYPE_RGBA ) &&
                   interlacing == PNG_INTERLACE_NONE ) ) {
            fprintf( stderr, "\nUnsupported PNG file type.\n" );
          } else {
            const int components =
              colorType == PNG_COLOR_TYPE_PALETTE ? 1
              : colorType == PNG_COLOR_TYPE_RGBA ? 4
              : 3;
            unsigned char* output = rgb;
            size_t bytes = height * sizeof (png_bytep);
            png_bytep* rowPointers = malloc( bytes );

            if ( ! rowPointers ) {
              fprintf( stderr,
                       "\nFailed to allocate %lu bytes for rowPointers.\n",
                       bytes );
            } else {
              memset( rowPointers, 0, bytes );
              bytes = height * width * components * sizeof (png_byte);

              {
                png_byte* pixels = malloc( bytes );

                if ( ! pixels ) {
                  fprintf( stderr,
                           "\nFailed to allocate %lu bytes for pixels.\n",
                           bytes );
                } else {
                  const size_t rowSize = width * components;
                  png_byte* pixelRow = pixels;
                  png_colorp palette = 0;
                  int row = 0;
                  memset( pixels, 0, bytes );

                  for ( row = 0; row < height; ++row, pixelRow += rowSize ) {
                    rowPointers[ row ] = pixelRow;
                  }

                  png_read_image( png, rowPointers );

                  if ( colorType == PNG_COLOR_TYPE_PALETTE ) {
                    png_get_PLTE( png, info, &palette, &unused_ );
                  }

                  /* Decode palette, rgb or rgba values to rgb: */

                  for ( row = 0; row < height; ++row ) {
                    const int rowIndex =
                      reverseImageRows ? height - row - 1 : row;
                    const unsigned char* const rowPixels =
                      rowPointers[ rowIndex ];
                    int column = 0;

                    for ( column = 0; column < width; ++column, output += 3 ) {

                      if ( colorType == PNG_COLOR_TYPE_PALETTE ) {
                        const int index = rowPixels[ column ];
                        output[ 0 ] = palette[ index ].red;
                        output[ 1 ] = palette[ index ].green;
                        output[ 2 ] = palette[ index ].blue;
                      } else {
                        const png_byte* const pixel =
                          rowPixels + column * components;
                        output[ 0 ] = pixel[ 0 ];
                        output[ 1 ] = pixel[ 1 ];
                        output[ 2 ] = pixel[ 2 ];
                      }
                    }
                  }

                  free( pixels ), pixels = 0;
                  result = 1;
                }
              }

              free( rowPointers ), rowPointers = 0;
            }
          }

          png_destroy_info_struct( png, &info );
        }

        png_destroy_read_struct( &png, 0, 0 );
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: writePNGFile - Write a PNG file.
INPUTS:  FILE* output  File to write to.
         int width     Image width in pixels.
         int height    Image height in pixels.
         unsigned char rgba[ height * width * 4 ]  RGBA pixel values in
                                                   OpenGL bottom-up row order.
RETURNS: int 1 if successful, else 0 and a message is printed to stderr.
******************************************************************************/

int writePNGFile( FILE* output, int width, int height, unsigned char rgba[] ) {
  int result = 0;

  if ( ! ( width > 0 && height > 0 ) ) {
    fprintf( stderr, "\nUnsupported PNG file dimensions.\n" );
  } else {
    png_structp png = png_create_write_struct( PNG_LIBPNG_VER_STRING, 0, 0, 0);

    if ( ! png ) {
      fprintf( stderr, "\nFailed to create PNG write structure.\n" );
    } else {
      png_infop info = png_create_info_struct( png );

      if ( ! info ) {
        fprintf( stderr, "\nFailed to create PNG info structure.\n" );
      } else {
        const int components = 4; /* R, G, B, A */
        const size_t rowSize = width * components;
        unsigned char* rowPointer = rgba + height * rowSize;
        const size_t bytes = height * sizeof (png_byte*);
        png_byte** rowPointers = malloc( bytes );

        if ( ! rowPointers ) {
          fprintf( stderr, "\nFailed to allocate %lu bytes for rowPointers.\n",
                   bytes );
        } else {
          int row = 0;
          memset( rowPointers, 0, bytes );

          /* Reverse row order from bottom-up OpenGL to top-down PNG: */

          for ( row = 0; row < height; ++row ) {
            rowPointer -= rowSize;
            rowPointers[ row ] = rowPointer;
          }

          png_init_io( png, output );
          png_set_IHDR( png, info, width, height,
                        8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
                        PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT );
          png_set_rows( png, info, rowPointers );
          png_write_png( png, info, PNG_TRANSFORM_IDENTITY, 0 );
          free( rowPointers ), rowPointers = 0;
          result = 1;
        }

        png_destroy_info_struct( png, &info );
      }

      png_destroy_write_struct( &png, 0 );
    }
  }

  return result;
}



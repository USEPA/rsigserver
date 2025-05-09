
#ifndef IMAGEFILE_H
#define IMAGEFILE_H

#ifdef __cplusplus
extern "C" {
#endif

  
/******************************************************************************
PURPOSE: ImageFile.c - Read subsets of *.bin image files.
NOTES:   
HISTORY: 2010-12-12 plessel.todd@epa.gov Created.
STATUS:  unreviewed, untested.
******************************************************************************/


/*================================ INCLUDES =================================*/

#include <Utilities.h> /* For Bounds. */

/*================================ FUNCTIONS ================================*/

unsigned char* readImageFile( const char* fileName, const Bounds clip,
                              int* width, int* height,
                              Bounds domain, Bounds corners );

#ifdef __cplusplus
}
#endif

#endif /* IMAGEFILE_H */



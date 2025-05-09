
#ifndef READDATA_H
#define READDATA_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: ReadData.h - Declare simple to use wrapper routines to read data from
         HRRR .grib2 files.

NOTES:   Uses GRIB2 library.

HISTORY: 2020-02-21 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <stdlib.h> /* For size_t. */

/*================================ FUNCTIONS ================================*/

extern int readData( const char* const fileName,
                     const int isVector2,
                     const size_t count,
                     double data[] );

#ifdef __cplusplus
}
#endif

#endif /* READDATA_H */




#ifndef M3IO_H
#define M3IO_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: M3IO.h - Declare convenience routines for writing M3IO files.

NOTES:

HISTORY: 2008/02, plessel.todd@epa.gov, Created.

STATUS:  unreviewed, tested.
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <BasicNumerics.h> /* For Integer, Real. */
#include <Grid.h>          /* For Grid. */
#include <Helpers.h>       /* For Name. */

/*================================ FUNCTIONS ================================*/

extern Integer writeM3IOHeader( Integer file,
                                Integer timesteps,
                                Integer hoursPerTimestep,
                                Integer firstTimestamp,
                                Integer variables,
                                Integer layers,
                                const Name variableNames[],
                                const Name variableUnits[],
                                const char* description,
                                const Grid* grid );

extern Integer writeM3IOGrid( const Grid* grid,
                              Integer timesteps, Integer layers, Integer file);

extern Integer writeM3IOData( Integer file, const Name variableName,
                              Integer timestep, Integer layers,
                              Integer rows, Integer columns,
                              const void* data );

extern void copyDataToGrid( Integer points,
                            const Integer rows[],
                            const Integer columns[],
                            const Real pointData[],
                            Real scale,
                            Integer gridLayers,
                            Integer gridRows,
                            Integer gridColumns,
                            Real gridData[] );

extern void copyIntDataToGrid( Integer points,
                               const Integer rows[],
                               const Integer columns[],
                               const Integer pointData[],
                               Integer gridLayers,
                               Integer gridRows,
                               Integer gridColumns,
                               Integer gridData[] );

extern void copyDataToGrid3( Integer points,
                             const Integer layers[],
                             const Integer rows[],
                             const Integer columns[],
                             const Real pointData[],
                             Real scale,
                             Integer gridLayers,
                             Integer gridRows,
                             Integer gridColumns,
                             Real gridData[] );

#ifdef __cplusplus
}
#endif

#endif /* M3IO_H */




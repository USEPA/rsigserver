#ifndef NETCDFUTILITIES_H
#define NETCDFUTILITIES_H

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************
PURPOSE: NetCDFUtilities.h - Declare convenience routines for NetCDF files.
NOTES:
HISTORY: 2007-12-01 plessel.todd@epa.gov
STATUS:  unreviewed tested
******************************************************************************/


/*================================ FUNCTIONS ================================*/

extern int printM3IOVariables( const char* const fileName );

extern int createNetCDFFile( const char* const fileName );

extern int openNetCDFFile( const char* const fileName, const char rw );

extern int closeNetCDFFile( const int file );

extern int flushNetCDFFile( const int file );

extern int endNetCDFHeader( const int file );

extern int getNetCDFDimension( const int file, const char* const name );

extern int checkNetCDFVariableId( const int file, const char* const name );

extern int getNetCDFVariableId( const int file, const char* const name );


extern int isNetCDFFloat( const int type );

extern int getM3IOVariableDimensions( const int file,
                                      const char* const variable,
                                      int dims[ 4 ] );

extern int getNetCDFVariableInfo( const int file, const int id,
                                  char name[], int* type, int* rank,
                                  int dimensions[], char units[],
                                  char description[] );

extern int getNetCDFStringAttribute( const int file, const int id,
                                     const char* const name,
                                     const size_t size, char value[] );

extern int getNetCDFDoubleAttribute( const int file,
                                     const char* const name,
                                     double* const value );

extern int getNetCDFFloatAttribute( const int file,
                                    const char* const name,
                                    float* const value );

extern int getNetCDFIntAttribute( const int file,
                                  const char* const name,
                                  int* const value );

extern int getNetCDFFloatArrayAttribute( const int file,
                                         const char* const name,
                                         const int count,
                                         float values[] );

extern int getM3IOFileTimeRange( const int file,
                                 int* const yyyymmddhh1,
                                 int* const yyyymmddhh2,
                                 int* timesteps,
                                 int* hoursPerTimestep );

extern int readM3IOVariable( const int file,
                             const int id,
                             const int time0,   const int time1,
                             const int layer0,  const int layer1,
                             const int row0,    const int row1,
                             const int column0, const int column1,
                             float array[] );

extern int createNetCDFDimension( const int file, const char* const name,
                                  const int size, int* const id );

extern int createNetCDFVariable( const int file,
                                 const char* const name,
                                 const char* const units,
                                 const char* const attributeName,
                                 const char* const attributeValue,
                                 const int rank,
                                 const int dimids[] );

extern int copyNetCDFAttribute( const int file,
                                const char* const name,
                                const int output );

extern int createNetCDFStringAttribute( const int file,
                                        const int id,
                                        const char* const name,
                                        const char* const value );

extern int createNetCDFIntAttribute( const int file,
                                     const char* const name,
                                     const int value );

extern int createNetCDFDoubleAttribute( const int file,
                                        const char* const name,
                                        const double value );

extern int createNetCDFFloatAttribute( const int file,
                                       const char* const name,
                                       const float value );

extern int createNetCDFFloatArrayAttribute( const int file,
                                            const char* const name,
                                            const int count,
                                            const float values[] );

extern int writeCOARDS2DVariable( const int file,
                                  const char* const name,
                                  const int rows,
                                  const int columns,
                                  const float array[] );

extern int writeCOARDSTimeVariables( const int file,
                                     const int timestep,
                                     const int yyyymmddhh );

extern int writeM3IOVariable( const int file,
                              const char* const name,
                              const int timestep,
                              const int layers,
                              const int rows,
                              const int columns,
                              const float array[] );

extern int writeTFLAGVariable( const int file,
                               const int timestep,
                               const int variables,
                               const int datetime,
                               const int array[] );


#ifdef __cplusplus
}
#endif

#endif /* NETCDFUTILITIES_H */




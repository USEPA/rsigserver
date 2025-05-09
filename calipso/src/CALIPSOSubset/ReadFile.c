/******************************************************************************
PURPOSE: ReadFile.c - Simple to use wrapper routines to read data from CALIPSO
                      HDF files.

NOTES:   Uses HDF-EOS and HDF4 libraries and libs they depend on (z, jpeg).

HISTORY: 2017-01-02 plessel.todd@epa.gov
STATUS: inchoate
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For strcmp(), strstr(), memset(). */

/* Declare simplified, const-correct prototypes of HDF-EOS routines called: */

extern int EHchkfid( int, const char*, int*, int*, unsigned char* );
extern int SWopen( const char*, int );
extern int SWclose( int );

/* Declare simplified, const-correct prototypes of HDF4 routines called: */

extern int SDnametoindex( int, const char* );
extern int SDselect( int, int );
extern int SDfindattr( int, const char* );
extern int SDattrinfo( int, int, char*, int*, int* );
extern int SDreadattr( int, int, void* );
extern int SDgetinfo( int, const char*, int*, int*, int*, int* );
extern int SDreaddata( int, const int*, const int*, const int*, void* );

extern int Vinitialize( int );
extern int Vfinish( int );
#define Vstart Vinitialize
#define Vend Vfinish
extern int VSattach( int, int, const char* );
extern int VSdetach( int );
extern int VSfind( int, const char* );
extern int VSsetfields( int, const char* );
extern int VSread( int, void*, int, int );

#include "Utilities.h" /* For MISSING_VALUE, LONGITUDE, IN_RANGE(), Bounds. */
#include "ReadFile.h"  /* For public interface. */

/*================================== TYPES ==================================*/

/* Type of file variable/dataset data: */

enum {
  CHAR   = 4,
  INT8   = 20, UINT8  = 21,
  INT16  = 22, UINT16 = 23,
  INT32  = 24, UINT32 = 25,
  REAL32 = 5,  REAL64 = 6
};

#define IS_VALID_TYPE( type ) \
  IN10( type, CHAR, INT8, UINT8, INT16, UINT16, INT32, UINT32, REAL32, REAL64 )

/*=========================== FORWARD DECLARATIONS ==========================*/

static size_t dimensionsProduct( const int rank, const int dimensions[] );

static int dimsMatch( const int rank,
                     const int dimensions1[], const int dimensions2[] );

static int readScaleFactorAttribute( const int file, const char* variable,
                                     double* value );

static int readTextAttribute( const int file,
                              const char* const variable,
                              const char* const attribute,
                              char value[] );

static int lookupAttribute( const int file, const char* const variable,
                            const char* const attribute,
                            int* const variableId, int* const attributeIndex,
                            int* const type, int* const count );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: openFile - Open HDF file for reading.
INPUTS:  const char* const fileName  Name of file to open for reading.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int openFile( const char* const fileName ) {
  int result = 0;
  assert( fileName ); assert( *fileName );
  DEBUG( fprintf( stderr, "\n%s\n", fileName ); )
  result = SWopen( fileName, 1 );

  if ( result == -1 ) {
    fprintf( stderr, "\a\nFailed to open HDF file for reading: %s.\n",
             fileName );
  }

  DEBUG( fprintf( stderr, "file = %d\n", result ); )
  assert( result >= -1 );
  return result;
}



/******************************************************************************
 PURPOSE: closeFile - Close HDF file.
 INPUTS:  int file  HDF file to close.
 ******************************************************************************/

void closeFile( const int file ) {
  assert( file > -1 );
  SWclose( file );
}



/******************************************************************************
PURPOSE: readFileBounds - Read file lon-lat bounds from metadata header.
INPUTS:  const int file  HDF file ID to read.
OUTPUTS: Bounds bounds   Bounds of file.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileBounds( const int file, Bounds bounds ) {
  int result = 0;
  const char* const attributeName = "coremetadata";
  int sdId = -1;
  int unused1 = 0;
  int type = 0;
  int count = 0;
  unsigned char unused2 = 0;
  char buffer[ 50000 ] = ""; /* HACK: avoid dynamic allocation delays. */

  assert( file > -1 ); assert( bounds );

  memset( buffer, 0, sizeof buffer );
  memset( bounds, 0, 2 * 2 * sizeof bounds[ 0 ][ 0 ] );

  if ( EHchkfid( file, attributeName, &unused1, &sdId, &unused2 ) != -1 ) {
    const int attributeIndex = SDfindattr( sdId, attributeName );

    if ( attributeIndex != -1 ) {

      if ( AND4( SDattrinfo(sdId, attributeIndex, buffer, &type, &count) != -1,
                 type == CHAR, count > 0, count < sizeof buffer ) ) {

        if ( SDreadattr( sdId, attributeIndex, buffer ) != -1 ) {
          const char* const matchString = "= MINLAT";
          const char* start = strstr( buffer, matchString );

          if ( start ) {
            start = strstr( start, "VALUE" );

            if ( start ) {
              const char* const format =
                "%*s%*s%lf"
                "%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lf"
                "%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lf"
                "%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lf";

              result =
                sscanf( start, format,
                        &bounds[ LATITUDE  ][ MINIMUM ],
                        &bounds[ LONGITUDE ][ MINIMUM ],
                        &bounds[ LATITUDE  ][ MAXIMUM ],
                        &bounds[ LONGITUDE ][ MAXIMUM ] ) == 4;
            }
          }
        }
      }
    }
  }

  if ( result ) {

    /* HACK: expand longitude when crossing the +/-180 line: */

    if ( bounds[ LONGITUDE ][ MINIMUM ] >
         bounds[ LONGITUDE ][ MAXIMUM ] ) {
         bounds[ LONGITUDE ][ MINIMUM ] = -180.0;
         bounds[ LONGITUDE ][ MAXIMUM ] =  180.0;
    }

    result = isValidBounds( (const double (*)[2]) bounds );
  }

  if ( ! result ) {
    fprintf( stderr, "\a\n\nInvalid file metadata for lon-lat bounds.\n" );
    memset( bounds, 0, 2 * 2 * sizeof bounds[ 0 ][ 0 ] );
  }

  DEBUG( fprintf( stderr, "swath bounds = [%lf %lf][%lf %lf]\n",
                 bounds[ LONGITUDE ][ MINIMUM ],
                 bounds[ LONGITUDE ][ MAXIMUM ],
                 bounds[ LATITUDE  ][ MINIMUM ],
                 bounds[ LATITUDE  ][ MAXIMUM ] ); )
  assert( IS_BOOL( result ) );
  assert( IMPLIES( result, isValidBounds( (const double (*)[2]) bounds ) ) );
  return result;
}



/******************************************************************************
PURPOSE: fileVariableExists - Does the named variable exist in file?
INPUTS:  const int file              HDF file ID to read.
         const char* const variable  Name of variable to query.
RETURNS: int 1 if exists, else 0.
******************************************************************************/

int fileVariableExists( const int file, const char* const variable ) {
  int result = 0;

  assert( file > -1 );
  assert( variable ); assert( *variable );

  {
    int sdId = -1;
    int unused1 = 0;
    unsigned char unused2 = 0;
    int status = EHchkfid( file, variable, &unused1, &sdId, &unused2 );

    if ( status != -1 ) {
      const int variableIndex = SDnametoindex( sdId, variable );

      if ( variableIndex != -1 ) {
        const int variableId = SDselect( sdId, variableIndex );

        if ( variableId != -1 ) {
          int type = 0;
          int rank = 0;
          int dims[ 32 ]; /* UGLY Big enough? */
          int unused3 = 0;
          memset( dims, 0, sizeof dims );
          status = SDgetinfo( variableId, 0, &rank, dims, &type, &unused3 );
          result =
            AND5( status != -1, rank == 2, IN3( type, 5, 22 ),
                  dims[ 0 ] > 0, dims[ 1 ] > 0 );
        }
      }
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readVariableDimensions - Read variable dimensions.
INPUTS:  const int file              HDF file ID to read.
         const char* const variable  Name of variable to query.
OUTPUTS: int* const rank             Number of dimensions. (2 or 3)
         int dimensions[ 32 ]        Dimensions.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readVariableDimensions( const int file, const char* const variable,
                            int* const rank, int dimensions[ 32 ] ) {
  int result = 0;

  assert( file > -1 );
  assert( variable ); assert( *variable );
  assert( rank ); assert( dimensions );

  {
    int sdId = -1;
    int unused1 = 0;
    unsigned char unused2 = 0;
    int status = EHchkfid( file, variable, &unused1, &sdId, &unused2 );

    if ( status == -1 ) {
      fprintf( stderr, "\a\n\nFailed to find SD-ID for %s.\n", variable );
    } else {
      const int variableIndex = SDnametoindex( sdId, variable );

      if ( variableIndex == -1 ) {
        fprintf( stderr, "\a\n\nFailed to find SDS index for %s.\n", variable);
      } else {
        const int variableId = SDselect( sdId, variableIndex );

        if ( variableId == -1 ) {
          fprintf( stderr, "\a\n\nFailed to select SDS ID for %s.\n",variable);
        } else {
          int type = 0;
          int unused3 = 0;
          status =
            SDgetinfo( variableId, 0, rank, dimensions, &type, &unused3 );

          if ( OR5( status == -1, ! IN3( *rank, 2, 3), ! IS_VALID_TYPE( type ),
                    dimensions[ 0 ] < 1, dimensions[ *rank - 1 ] < 1 ) ) {
            fprintf( stderr, "\a\n\nFailed to get valid info on %s.\n",
                     variable );
          } else {
            result = 1;
          }
        }
      }
    }
  }

  if ( ! result ) {
    fprintf(stderr, "\a\nFailed to read valid dimensions of %s.\n",variable);
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file variable data.
INPUTS:  const int file                  HDF file ID to read.
         const char* const variable      Name of variable to read.
         const int rank                  Number of variable dimensions.
         const int dimensions[ rank ]    Dimensions of variable to read.
OUTPUTS: char units[ 80 ]                Units of variable.
         double data[ dimensions[0] * ... * dimensions[ rank-1 ] ]
                                         Data for variable.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileData( const int file, const char* const variable,
                  const int rank, const int dimensions[],
                  char units[ 80 ], double data[] ) {

  int result = 0;
  int status = -1;
  int sdId = -1;
  int unused1 = 0;
  unsigned char unused2 = 0;

  assert( file >= 0 );
  assert( variable ); assert( *variable );
  assert( rank > 0 ); assert( dimensions );
  assert( dimensions[ 0 ] > 0 ); assert( dimensions[ rank - 1 ] > 0 );
  assert( units ); assert( data );

  DEBUG( fprintf( stderr,
                  "readFileData( %s, rank = %d, dimensions = [%d .. %d]\n",
                 variable, rank, dimensions[ 0 ], dimensions[ rank - 1 ] ) );

  *units = '\0';
  status = EHchkfid( file, variable, &unused1, &sdId, &unused2 );

  if ( status == -1 ) {
    fprintf( stderr, "\a\n\nFailed to find SD-ID for %s.\n", variable );
  } else {
    const int variableIndex = SDnametoindex( sdId, variable );

    if ( variableIndex == -1 ) {
      fprintf( stderr, "\a\n\nFailed to find SDS index for %s.\n", variable);
    } else {
      const int variableId = SDselect( sdId, variableIndex );

      if ( variableId == -1 ) {
        fprintf( stderr, "\a\n\nFailed to select SDS ID for %s.\n",
                 variable );
      } else {
        int type = 0;
        int rank0 = 0;
        int dims[ 32 ]; /* UGLY Big enough? */
        int unused3 = 0;
        memset( dims, 0, sizeof dims );
        status = SDgetinfo( variableId, 0, &rank0, dims, &type, &unused3 );

        if ( OR4( status == -1, rank0 != rank, ! IS_VALID_TYPE( type ),
                  ! dimsMatch( rank, dims, dimensions ) ) ) {
          fprintf( stderr, "\a\n\nFailed to get valid/matching info on %s.\n",
                   variable );
        } else {
          const int starts[ 3 ] = { 0, 0, 0 };

          if ( SDreaddata( variableId, starts, 0, dims, data ) == -1 ) {
            fprintf( stderr, "\a\n\nFailed to read '%s'.\n", variable );
          } else {
            const size_t count = dimensionsProduct( rank, dimensions );
            double scaleFactor = 1.0; /* Default value if not found. */
            readScaleFactorAttribute( file, variable, &scaleFactor );
            readTextAttribute( file, variable, "units", units );
            result = 1;
            DEBUG( fprintf( stderr, "type = %d, rank = %d, "
                            "dims = [%d ... %d], count = %lu, units = %s, "
                            "scaleFactor = %lf\n",
                            type, rank, dims[ 0 ], dims[ rank - 1 ], count,
                            units, scaleFactor ); )

            switch ( type ) {
            case INT8:
              expandInt8( count, data );
              break;
            case UINT8:
              expandUint8( count, data );
              break;
            case INT16:
              expandInt16( count, data );
              break;
            case UINT16:
              expandUint16( count, data );
              break;
            case INT32:
              expandInt32( count, data );
              break;
            case UINT32:
              expandUint32( count, data );
              break;
            case REAL32:
              expandReals( count, data );
              break;
            default:
              assert( type == REAL64 );
              break;
            }

            if ( scaleFactor != 1.0 ) {
              scaleValues( scaleFactor, count, data );
            }
          }
        }
      }
    }
  }

  DEBUG( fprintf( stderr,
                  "%s (%s) data = [%lf ... %lf], result = %d\n",
                  variable, units, data[ 0 ],
                  data[ dimensionsProduct( rank, dimensions ) - 1 ],
                  result ); )
  return result;
}



/******************************************************************************
PURPOSE: readFileVData - Read Vdata of variable.
INPUTS:  const int file              HDF file ID to read.
         const char* const variable  Name of variable to read.
         const int count             Number of data values to read.
OUTPUTS: double data[ count ]        Data for variable.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

int readFileVData( const int file, const char* const variable,
                   const int count, double data[] ) {

  int result = 0;
  int unused1;
  int HDFId = -1;
  unsigned char unused2 = 0;

  assert( file >= 0 ); assert( variable ); assert( *variable );
  assert( count > 0 ); assert( data );

  DEBUG( fprintf( stderr, "readFileVData( %s, %d )\n", variable, count ) );

  if ( EHchkfid( file, "unused", &HDFId, &unused1, &unused2 ) != 0 ) {
    fprintf( stderr, "\a\n\nFailed to find valid HDF ID for '%s'.\n",
            variable );
  } else if ( Vstart( HDFId ) != 0 ) {
    fprintf( stderr, "\a\n\nFailed to start Vdata interface for '%s'.\n",
            variable );
  } else {
    const int vdataRef = VSfind( HDFId, "metadata" );

    if ( vdataRef == 0 ) {
      fprintf( stderr, "\a\n\nFailed to find valid Vdata ref for '%s'.\n",
              variable );
    } else {
      const int vdataId = VSattach( HDFId, vdataRef, "r" );

      if ( vdataId == 0 ) {
        fprintf( stderr, "\a\n\nFailed to attach to Vdata ref for '%s'.\n",
                variable );
      } else if ( VSsetfields( vdataId, variable ) != 0 ) {
        fprintf( stderr, "\a\n\nFailed to set field to '%s'.\n", variable );
      } else if ( VSread( vdataId, data, 1, 0 ) != 1 ) {
        fprintf( stderr, "\a\n\nFailed to read Vdata for '%s'.\n", variable );
      } else {
        expandReals( count, data ); /* 32-bit to 64-bit. */
        result = 1;
      }

      VSdetach( vdataId );
    }

    Vend( HDFId );
  }

  DEBUG( fprintf( stderr, "result = %d\n", result ) );
  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/




/******************************************************************************
 PURPOSE: dimensionsProduct - Product of dimensions.
 INPUTS:  const int rank                Number of dimensions.
          const int dimensions[ rank ]  Set of dimensions to multiply.
 RETURNS: size_t product of dimensions.
 ******************************************************************************/

static size_t dimensionsProduct( const int rank, const int dimensions[] ) {
  size_t result = 1;
  int index = 0;
  assert( rank > 0 ); assert( dimensions );

  for ( index = 0; index < rank; ++index ) {
    assert( dimensions[ index ] > 0 );
    result *= dimensions[ index ];
  }

  return result;
}



/******************************************************************************
 PURPOSE: dimsMatch - Do dimensions match?
 INPUTS:  const int rank                 Number of dimensions.
          const int dimensions1[ rank ]  1st set of dimensions to compare.
          const int dimensions2[ rank ]  2nd set of dimensions to compare.
 RETURNS: int 1 if matched, else 0.
 ******************************************************************************/

static int dimsMatch( const int rank,
                      const int dimensions1[], const int dimensions2[] ) {
  int result = 1;
  int index = 0;
  assert( rank > 0 ); assert( dimensions1 ); assert( dimensions2 );

  for ( index = 0; AND2( result, index < rank ); ++index ) {
    result = dimensions1[ index ] == dimensions2[ index ];
  }

  return result;
}



/******************************************************************************
 PURPOSE: readScaleFactorAttribute - Read possible scale_factor attribute of
          variable.
 INPUTS:  const int file         Id of file or dataset.
          const char* variable   E.g., "Aerosol_Layer_Fraction".
 OUTPUTS: double* value
 RETURNS: int 1 if found and successfully read, else 0.
 ******************************************************************************/

static int readScaleFactorAttribute( const int file, const char* variable,
                                     double* value ) {

  const char* const attribute = "scale_factor";
  int result = 0;
  int type = -1;
  int count = 0;
  int variableId = -1;
  int attributeIndex = -1;

  assert( file > -1 ); assert( variable ); assert( *variable );
  assert( value );

  *value = 1.0; /* Default value if attribute is not found. */

  if ( lookupAttribute( file, variable, attribute,
                        &variableId, &attributeIndex, &type, &count ) ) {

    if ( ! AND2( IS_VALID_TYPE( type ), count == 1 ) ) {
      fprintf( stderr,
              "\a\n\nUnsupported type (%d) count (%d) of attribute '%s'.\n",
              type, count, attribute );
    } else if ( SDreadattr( variableId, attributeIndex, value ) == -1 ) {
      fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
              attribute );
    } else {

      switch ( type ) {
      case INT8:
        expandInt8( 1, value );
        break;
      case UINT8:
        expandUint8( 1, value );
        break;
      case INT16:
        expandInt16( 1, value );
        break;
      case UINT16:
        expandUint16( 1, value );
        break;
      case INT32:
        expandInt32( 1, value );
        break;
      case UINT32:
        expandUint32( 1, value );
        break;
      case REAL32:
        expandReals( 1, value );
        break;
      default:
        assert( type == REAL64 );
        break;
      }

      result = 1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readTextAttribute - Read named attribute of variable.
INPUTS:  const int file         Id of file or dataset.
         const char* variable   E.g., "Optical_Depth_Land_And_Ocean".
         const char* attribute  E.g., "units".
OUTPUTS: char value[ 80 ]       Value read if successful else unwritten.
 RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int readTextAttribute( const int file, const char* variable,
                              const char* attribute, char value[ 80 ] ) {

  int result = 0;
  int type = -1;
  int count = 0;
  int variableId = -1;
  int attributeIndex = -1;

  assert( file > -1 ); assert( variable ); assert( *variable );
  assert( attribute ); assert( *attribute ); assert( value );

  if ( lookupAttribute( file, variable, attribute,
                        &variableId, &attributeIndex, &type, &count ) ) {

    if ( OR3( type != CHAR, count < 1, count > 79 ) ) {
      fprintf( stderr,
               "\a\n\nUnsupported type (%d) count (%d) of attribute '%s'.\n",
               type, count, attribute );
    } else if ( SDreadattr( variableId, attributeIndex, value ) == -1 ) {
      fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
               attribute );
    } else {
      result = 1;
    }
  } else if ( strcmp( attribute, "units" ) == 0 ) {
    value[ 0 ] = '-';
    value[ 1 ] = '\0';
    result = 1;
  }

  if ( strcmp( attribute, "units" ) == 0 ) { /* Convert problematic units: */
    spacesToUnderscores( value );

    if ( OR3( strcmp( value, "mb" ) == 0,
              strcmp( value, "millibars" ) == 0,
              strcmp( value, "hPA" ) == 0 ) ) {
      strcpy( value, "hPa" ); /* Go metric every inch of the way! */
    } else if ( ! strcmp( variable, "Profile_Time" ) ) {
      strcpy( value, "seconds_since_1993-01-01" );
    } else if ( ! strcmp( variable, "Profile_UTC_Time" ) ) {
      strcpy( value, "yyyymmdd.f" );
    } else if ( strcmp( value, "NoUnits" ) == 0 ) {
      value[ 0 ] = '-';
      value[ 1 ] = '\0';
    } else if ( strcmp( value, "None" ) == 0 ) {
      value[ 0 ] = '-';
      value[ 1 ] = '\0';
    } else if ( strcmp( value, "none" ) == 0 ) {
      value[ 0 ] = '-';
      value[ 1 ] = '\0';
    } else if ( strstr( value, "egrees" ) ) {
      value[ 0 ] = 'd';
      value[ 1 ] = 'e';
      value[ 2 ] = 'g';
      value[ 3 ] = '\0';
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: lookupAttribute - Lookup attribute of named variable.
INPUTS:  const int file         File id.
         const char* variable   E.g., "Optical_Depth_Land_And_Ocean".
         const char* attribute  E.g., "units".
OUTPUTS: int* const variableId        Id number of variable in file.
         int* const attributeIndex    Index of attribute of variable.
         int* const type              Type code of attribute. E.g., 4 = char.
         int* const count             Length of attribute. E.g., 16.
******************************************************************************/

static int lookupAttribute( const int file, const char* const variable,
                            const char* const attribute,
                            int* const variableId, int* const attributeIndex,
                            int* const type, int* const count ) {

  int result = 0;
  int sdId = 0;
  int unused1 = 0;
  unsigned char unused2 = 0;

  assert( file > -1 ); assert( variable ); assert( *variable );
  assert( attribute ); assert( *attribute ); assert( variableId );
  assert( attributeIndex ); assert( type ); assert( count );

  if ( EHchkfid( file, variable, &unused1, &sdId, &unused2 ) == -1 ) {
    fprintf( stderr, "\a\n\nInvalid variable '%s'.\n", variable );
  } else {
    const int variableIndex = SDnametoindex( sdId, variable );

    if ( variableIndex == -1 ) {
      fprintf( stderr, "\a\n\nInvalid variable '%s'.\n", variable );
    } else {
      *variableId = SDselect( sdId, variableIndex );

      if ( *variableId == -1 ) {
        fprintf( stderr, "\a\n\nInvalid variable '%s'.\n", variable );
      } else {
        *attributeIndex = SDfindattr( *variableId, attribute );

        if ( *attributeIndex == -1 ) {

          if ( AND2( strcmp( attribute, "units" ),
                     strcmp( attribute, "scale_factor" ) ) ) {
            fprintf( stderr, "\a\n\nInvalid attribute '%s'.\n", attribute );
          }
        } else {
          char unused[ 64 ] = "";

          if ( OR3( SDattrinfo( *variableId, *attributeIndex, unused,
                                type, count ) == -1,
                    ! IS_VALID_TYPE( *type ), *count < 1 ) ) {
            fprintf( stderr, "\a\n\nInvalid attribute '%s'.\n", attribute );
          } else {
            result = 1;
          }
        }
      }
    }
  }

  if ( ! result ) {
    *variableId = *attributeIndex = *type = -1;
    *count = 0;
  }

  return result;
}




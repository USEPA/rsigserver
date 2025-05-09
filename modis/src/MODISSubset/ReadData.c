/******************************************************************************
PURPOSE: ReadData.c - Simple to use wrapper routines to read data from MODIS
                      HDF-EOS files.

NOTES:   Uses HDF-EOS and HDF libraries and libs they depend on (z, jpeg).

HISTORY: 2017-12-19 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */
#include <string.h> /* For strcmp(), strstr(), strncpy(). */
#include <ctype.h>  /* For tolower(). */
#include <float.h>  /* For DBL_MAX. */

/* Declare simplified, const-correct prototypes of HDF-EOS routines called: */

extern int EHchkfid( int, const char*, int*, int*, unsigned char* );
extern int SWopen( const char*, int );
extern int SWclose( int );
extern int SWattach( int, const char* );
extern int SWdetach( int );
extern int SWfieldinfo( int, const char*, int*, int*, int*, char* );
extern int SWreadfield( int, const char*, int*, int*, int*, void* );

/* Declare simplified, const-correct prototypes of HDF4 routines called: */

extern int SDnametoindex( int, const char* );
extern int SDselect( int, int );
extern int SDfindattr( int, const char* );
extern int SDattrinfo( int, int, char*, int*, int* );
extern int SDreadattr( int, int, void* );

#include "Utilities.h" /* For MISSING_VALUE, LONGITUDE, IN_RANGE(), Bounds. */
#include "ReadData.h"  /* For public interface. */

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

/*========================= PRIVATE GLOBAL VARIABLES ========================*/

static int swathId = -1; /* Attached/detached when open/closeFile() called. */

/*=========================== FORWARD DECLARATIONS ==========================*/

static void parseSwathNameFromFileName( const char* const fileName,
                                        char swathName[ 6 ] );

static size_t decodeAndFilterData( const size_t points,
                                   const double scale,
                                   const double offset,
                                   const double validMinimum,
                                   const double validMaximum,
                                   const int type,
                                   void* data );

static int readAttribute( const int file,
                          const char* const variable,
                          const char* const attribute,
                          double* value );

static int readAttribute2( const int file,
                           const char* const variable,
                           const char* const attribute,
                           double values[] );

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
PURPOSE: openFile - Open HDF-EOS file for reading.
INPUTS:  const char* const fileName  Name of file to open for reading.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int openFile( const char* const fileName ) {
  int result = 0;
  assert( fileName ); assert( *fileName );
  DEBUG( fprintf( stderr, "%s\n", fileName ); )

  if ( swathId != -1 ) {
    SWdetach( swathId ), swathId = -1;
  }

  result = SWopen( fileName, 1 );

  if ( result == -1 ) {
    fprintf( stderr, "\a\nFailed to open HDF-EOS file for reading: %s.\n",
             fileName );
  } else {
    char swathName[ 6 ] = ""; /* Parse from file name: mod04 or mod06. */
    memset( swathName, 0, sizeof swathName );
    parseSwathNameFromFileName( fileName, swathName );
    swathId = SWattach( result, swathName );

    if ( swathId == -1 ) {
      fprintf( stderr, "\a\nFailed to attach to swath: %s.\n", swathName );
      SWclose( result ), result = -1;
    }
  }

  DEBUG( fprintf( stderr, "file = %d, swathId = %d\n", result, swathId ); )
  assert( result >= -1 );
  assert( OR2( AND2( result == -1, swathId == -1 ),
               AND2( result >  -1, swathId >  -1 ) ) );
  return result;
}



/******************************************************************************
PURPOSE: closeFile - Close HDF-EOS file.
INPUTS:  int file  HDF-EOS file to close.
******************************************************************************/

void closeFile( const int file ) {
  assert( file > -1 );

  if ( swathId > -1 ) {
    SWdetach( swathId ), swathId = -1;
  }

  SWclose( file );
  assert( swathId == -1 );
}



/******************************************************************************
PURPOSE: readFileBounds - Read file lon-lat bounds from metadata header.
INPUTS:  const int file  HDF-EOS file ID to read.
OUTPUTS: Bounds bounds   Bounds of file.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileBounds( const int file, Bounds bounds ) {
  int result = 0;
  const char* const attributeName = "CoreMetadata.0";
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
          const char* const matchString = "= EASTBOUNDINGCOORDINATE";
          const char* start = strstr( buffer, matchString );

          if ( start ) {
            start = strstr( start, "VALUE" );

            if ( start ) {
              const char* const format =
                "%*s%*s%lf"
                "%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lf"
                "%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lf"
                "%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%lf";

              result = sscanf( start, format,
                               &bounds[ LONGITUDE ][ MAXIMUM ],
                               &bounds[ LONGITUDE ][ MINIMUM ],
                               &bounds[ LATITUDE  ][ MINIMUM ],
                               &bounds[ LATITUDE  ][ MAXIMUM ] ) == 4;

              if ( ! result ) {
                start = strstr( buffer, "= WESTBOUNDINGCOORDINATE" );

                if ( start ) {
                  start = strstr( start, "VALUE" );

                  if ( start ) {
                    result = sscanf( start, format,
                                     &bounds[ LONGITUDE ][ MINIMUM ],
                                     &bounds[ LATITUDE  ][ MAXIMUM ],
                                     &bounds[ LONGITUDE ][ MAXIMUM ],
                                     &bounds[ LATITUDE  ][ MINIMUM ] ) == 4;
                  }
                }
              }
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
PURPOSE: readFileDimensions - Read file swath row/column dimensions.
INPUTS:  const int file   HDF-EOS file ID to read.
OUTPUTS: size_t* rows     Rows in swath.
         size_t* columns  Columns in swath.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileDimensions( const int file,
                        size_t* const rows, size_t* const columns ) {
  int result = 0;

  assert( file > -1 ); assert( swathId > -1 );
  assert( rows ); assert( columns );

  {
    const char* const variable = "Longitude";
    int type = 0;
    int rank = 0;
    int dims[ 32 ];
    int status = 0;
    dims[ 0 ] = dims[ 1 ] = 0;
    status = SWfieldinfo( swathId, variable, &rank, dims, &type, 0 );

    if ( OR4( status == -1, rank != 2, dims[ 0 ] < 1, dims[ 1 ] < 1 ) ) {
      fprintf( stderr, "\a\n\nFailed to get valid info on Longitude.\n");
    } else {
      *rows    = dims[ 0 ];
      *columns = dims[ 1 ];
      result = 1;
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read valid dimensions of Longitude.\n" );
  }

  DEBUG( fprintf( stderr, "swath dimensions = %lu x %lu\n", *rows, *columns );)
  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file swath data.
INPUTS:  const int file                 HDF-EOS file ID to read.
         const char* const variable     Name of variable to read.
         const size_t rows              Rows of variable to read.
         const size_t columns           Columns of variable to read.
OUTPUTS: char units[ 80 ]               Units of variable.
         double data[ rows * columns ]  Swath data for variable.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
NOTES:   Data is filtered by QC flag - i.e., set to MISSING_VALUE if QC != 2
         (good) or 3 (very good).
******************************************************************************/

int readFileData( const int file, const char* const variable,
                  const size_t rows, const size_t columns,
                  char units[ 80 ], double data[] ) {

  int result = 0;

  assert( file >= 0 ); assert( swathId > -1 );
  assert( variable ); assert( *variable );
  assert( rows != 0 ); assert( columns != 0 ); assert( data );

  *units = '\0';

  {
    int type = 0;
    int rank = 0;
    int dims[ 32 ];
    int status = 0;
    dims[ 0 ] = dims[ 1 ] = 0;
    status = SWfieldinfo( swathId, variable, &rank, dims, &type, 0 );

    if ( OR5( status == -1, rank != 2, ! IS_VALID_TYPE( type ),
              dims[ 0 ] != rows, dims[ 1 ] != columns ) ) {
      fprintf( stderr, "\a\n\nFailed to get valid/matching info on %s.\n",
               variable );
      DEBUG( fprintf( stderr, "status = %d, rank = %d, type = %d, "
                      "dims[ 0 ] = %d, dims[ 1 ] = %d\n",
                      status, rank, type, dims[ 0 ], dims[ 1 ] ); )
    } else {
      const size_t points = rows * columns;

      if ( SWreadfield( swathId, variable, 0, 0, 0, data ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read %s.\n", variable );
      } else {
        double scale = 1.0;
        double offset = 0.0;
        double validRange[ 2 ] = { -DBL_MAX, DBL_MAX };
        units[ 0 ] = '-';
        units[ 1 ] = '\0';

        readAttribute( file, variable, "scale_factor", &scale );
        readAttribute( file, variable, "add_offset", &offset );
        readAttribute2( file, variable, "valid_range", validRange );
        readTextAttribute( file, variable, "units", units );
        result =
          decodeAndFilterData( points, scale, offset,
                               validRange[ MINIMUM ], validRange[ MAXIMUM ],
                               type, data ) != 0;
      }
    }
  }

  DEBUG( fprintf( stderr,
                  "%s (%s) data = [%lf ... %lf], result = %d\n",
                  variable, units, data[ 0 ], data[ rows * columns - 1 ],
                  result ); )
  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: parseSwathNameFromFileName - Parse swath name from MODIS file name.
INPUTS:  const char* fileName  Name of MODIS file to parse.
OUTPUTS: char swathName[ 1 ]   Name of swath, e.g., "mod04".
******************************************************************************/

static void parseSwathNameFromFileName( const char* fileName,
                                        char swathName[ 6 ] ) {

  const char* prefix = 0;
  int index = 0;

  assert( fileName ); assert( *fileName ); assert( swathName );

  prefix = strrchr( fileName, '/' );

  if ( prefix == 0 ) {
    prefix = fileName;
  } else {
    ++prefix;
  }

  for ( ; AND3( index < 5, *prefix, *prefix != '_' ); ++index, ++prefix ) {
    swathName[ index ] = tolower( *prefix );
  }

  swathName[ 1 ] = 'o'; /* MYD04 files use swath name 'mod04'! */
  swathName[ 5 ] = '\0';
  assert( strlen( swathName ) == 5 );

  DEBUG( fprintf( stderr, "swathName = %s\n", swathName ); )
}



/******************************************************************************
PURPOSE: decodeAndFilterData - Decode and filter data based on quality flags.
INPUTS:  const size_t points        Number of data values read.
         const double scale         Multiplier for decoding data.
         const double offset        Offset for decoding data.
         const double validMinimum  Minimum valid value of encoded data.
         const double validMaximum  Minimum valid value of encoded data.
         void* const data           Encoded data[ points ].
OUTPUTS: void* const data           Decoded, filtered double data[ points ].
RETURNS: size_t number of valid/unfiltered decoded values.
******************************************************************************/

static size_t decodeAndFilterData( const size_t points,
                                   const double scale,
                                   const double offset,
                                   const double validMinimum,
                                   const double validMaximum,
                                   const int type,
                                   void* const data ) {

  size_t result = 0;
  int applyValidRange = 0;
  size_t point = points;
  double* const ddata = (double*) data;

  assert( points );
  assert( IN_RANGE( scale, -DBL_MAX, DBL_MAX ) ); /* Finite. */
  assert( scale != 0.0 );
  assert( IN_RANGE( offset, -DBL_MAX, DBL_MAX ) ); /* Finite. */
  assert( IN_RANGE( validMinimum, -DBL_MAX, DBL_MAX ) ); /* Finite. */
  assert( IN_RANGE( validMaximum, -DBL_MAX, DBL_MAX ) ); /* Finite. */
  assert( IS_VALID_TYPE( type ) );
  assert( data );
  assert( sizeof (char)   == 1 );
  assert( sizeof (short)  == 2 );
  assert( sizeof (int)    == 4 );
  assert( sizeof (float)  == 4 );
  assert( sizeof (double) == 8 );

  applyValidRange = validMinimum <= validMaximum; /* False for 8-bit QA flags */

  DEBUG( fprintf( stderr,
                  "decodeAndFilterData: scale = %lf, offset = %lf, "
                  "validMinimum = %lf, validMaximum = %lf, type = %d\n",
                  scale, offset, validMinimum, validMaximum, type ); )

  switch ( type ) {
  case INT8:
    {
      const signed char* const encodedData = (signed char*) data;

      while ( point-- ) {
        const signed char encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  case UINT8:
    {
      const unsigned char* const encodedData = (unsigned char*) data;

      while ( point-- ) {
        const unsigned char encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  case INT16:
    {
      const signed short* const encodedData = (signed short*) data;

      while ( point-- ) {
        const signed short encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  case UINT16:
    {
      const unsigned short* const encodedData = (unsigned short*) data;

      while ( point-- ) {
        const unsigned short encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  case INT32:
    {
      const signed int* const encodedData = (signed int*) data;

      while ( point-- ) {
        const signed int encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  case UINT32:
    {
      const unsigned int* const encodedData = (unsigned int*) data;

      while ( point-- ) {
        const unsigned int encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  case REAL32:
    {
      const float* const encodedData = (float*) data;

      while ( point-- ) {
        const float encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  default:
    {
      const double* const encodedData = (double*) data;
      assert( type == REAL64 );

      while ( point-- ) {
        const double encodedValue = encodedData[ point ];
        const int valid =
          IMPLIES( applyValidRange,
                   IN_RANGE( encodedValue, validMinimum, validMaximum ) );

        if ( valid ) {
          const double decodedValue = ( encodedValue - offset ) * scale;
          ddata[ point ] = decodedValue;
          ++result;
        } else {
          ddata[ point ] = MISSING_VALUE;
        }
      }
    }
    break;
  } /* End of switch on data type. */

  assert( result <= points );
  return result;
}



/******************************************************************************
PURPOSE: readAttribute - Read a real dataset attribute of an HDF-EOS file.
INPUTS:  const int file          Id of file.
         const char* const variable   Name of variable/dataset with attribute.
         const char* const attribute  Attribute name to read.
OUTPUTS: double* value  Value read if successful else unwritten.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int readAttribute( const int file, const char* variable,
                          const char* attribute, double* value ) {

  int result = 0;
  int type = -1;
  int count = 0;
  int variableId = -1;
  int attributeIndex = -1;

  assert( file > -1 ); assert( variable ); assert( *variable );
  assert( attribute ); assert( *attribute ); assert( value );

  if ( lookupAttribute( file, variable, attribute,
                        &variableId, &attributeIndex, &type, &count ) ) {

    if ( OR2( type != REAL64, count != 1 ) ) {
      fprintf( stderr,
               "\a\n\nUnsupported type (%d) count (%d) of attribute '%s'.\n",
               type, count, attribute );
    } else if ( SDreadattr( variableId, attributeIndex, value ) == -1 ) {
      fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
               attribute );
    } else {
      result = 1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readAttribute2 - Read a real-2 dataset attribute of an HDF-EOS file.
INPUTS:  const int file          Id of file.
         const char* const variable   Name of variable/dataset with attribute.
         const char* const attribute  Attribute name to read. "valid_range".
OUTPUTS: double values[ 2 ]           Values read if successful else unwritten.
RETURNS: int 1 if successful, else 0.
******************************************************************************/

static int readAttribute2( const int file, const char* variable,
                           const char* attribute, double values[ 2 ] ) {

  int result = 0;
  int type = -1;
  int count = 0;
  int variableId = -1;
  int attributeIndex = -1;

  assert( file > -1 ); assert( variable ); assert( *variable );
  assert( attribute ); assert( *attribute ); assert( values );

  if ( lookupAttribute( file, variable, attribute,
                        &variableId, &attributeIndex, &type, &count ) ) {

    if ( OR2( ! IS_VALID_TYPE( type ), count != 2 ) ) {
      fprintf( stderr,
               "\a\n\nUnsupported type (%d) count (%d) of attribute '%s'.\n",
               type, count, attribute );
    } else if ( type == INT8 ) {
      signed char ivalues[ 2 ] = { 0, 0 };
   
      if ( SDreadattr( variableId, attributeIndex, ivalues ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else {
        values[ MINIMUM ] = ivalues[ MINIMUM ];
        values[ MAXIMUM ] = ivalues[ MAXIMUM ];
        result = 1;
      }
    } else if ( type == UINT8 ) {
      unsigned char ivalues[ 2 ] = { 0, 0 };
   
      if ( SDreadattr( variableId, attributeIndex, ivalues ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else {
        values[ MINIMUM ] = ivalues[ MINIMUM ];
        values[ MAXIMUM ] = ivalues[ MAXIMUM ];
        result = 1;
      }
    } else if ( type == INT16 ) {
      signed short ivalues[ 2 ] = { 0, 0 };
   
      if ( SDreadattr( variableId, attributeIndex, ivalues ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else {
        values[ MINIMUM ] = ivalues[ MINIMUM ];
        values[ MAXIMUM ] = ivalues[ MAXIMUM ];
        result = 1;
      }
    } else if ( type == UINT16 ) {
      unsigned short ivalues[ 2 ] = { 0, 0 };
   
      if ( SDreadattr( variableId, attributeIndex, ivalues ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else {
        values[ MINIMUM ] = ivalues[ MINIMUM ];
        values[ MAXIMUM ] = ivalues[ MAXIMUM ];
        result = 1;
      }
    } else if ( type == INT32 ) {
      signed int ivalues[ 2 ] = { 0, 0 };
   
      if ( SDreadattr( variableId, attributeIndex, ivalues ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else {
        values[ MINIMUM ] = ivalues[ MINIMUM ];
        values[ MAXIMUM ] = ivalues[ MAXIMUM ];
        result = 1;
      }
    } else if ( type == UINT32 ) {
      unsigned int ivalues[ 2 ] = { 0, 0 };
   
      if ( SDreadattr( variableId, attributeIndex, ivalues ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else {
        values[ MINIMUM ] = ivalues[ MINIMUM ];
        values[ MAXIMUM ] = ivalues[ MAXIMUM ];
        result = 1;
      }
    } else if ( type == REAL32 ) {
      float fvalues[ 2 ] = { 0.0, 0.0 };

      if ( SDreadattr( variableId, attributeIndex, fvalues ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else if ( AND2( IN_RANGE( fvalues[ MINIMUM ], -DBL_MAX, DBL_MAX ),
                        IN_RANGE( fvalues[ MAXIMUM ], -DBL_MAX, DBL_MAX ) ) ) {
        values[ MINIMUM ] = fvalues[ MINIMUM ];
        values[ MAXIMUM ] = fvalues[ MAXIMUM ];
        result = 1;
      }
    } else {
      assert( type == REAL64 );

      if ( SDreadattr( variableId, attributeIndex, values ) == -1 ) {
        fprintf( stderr, "\a\n\nFailed to read value of attribute '%s'.\n",
                 attribute );
      } else if ( AND2( IN_RANGE( values[ MINIMUM ], -DBL_MAX, DBL_MAX ),
                        IN_RANGE( values[ MAXIMUM ], -DBL_MAX, DBL_MAX ) ) ) {
        result = 1;
      }
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

          if ( strcmp( attribute, "units" ) ) {
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




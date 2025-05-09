/******************************************************************************
PURPOSE: ReadFile.c - Simple to use wrapper routines to read data from
                      ceilometer HDF5 files.

NOTES:   Uses HDF5 libraries and libs they depend on (curl, z, dl).

HISTORY: 2017-11-29 plessel.todd@epa.gov
STATUS: unreviewed tested
******************************************************************************/

/*================================ INCLUDES =================================*/

#include <assert.h> /* For assert(). */
#include <stdio.h>  /* For stderr, fprintf(). */

#include <hdf5.h> /* For H5*. */

#include "ReadFile.h"  /* For public interface. */

/*=========================== FORWARD DECLARATIONS ==========================*/

static int matchedDatasetDimensions( const int dataset,
                                     const size_t dimension0,
                                     const size_t dimension1 );

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: openFile - Open HDF5 file for reading.
INPUTS:  const char* const fileName   Name of HDF5 file to open for reading.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int openFile( const char* const fileName ) {
  int result = 0;
  assert( fileName ); assert( *fileName );
  result = H5Fopen( fileName, 0, 0 );

  if ( result == -1 ) {
    fprintf( stderr, "\a\nFailed to open HDF5 file for reading: %s.\n",
             fileName );
  }

  return result;
}



/******************************************************************************
PURPOSE: closeFile - Close HDF5 file.
INPUTS:  const int file HDF5 file to close.
******************************************************************************/

void closeFile( const int file ) {
  H5Fclose( file );
}



/******************************************************************************
PURPOSE: openDataset - Open HDF5 dataset for reading.
INPUTS:  const char* const fileName   Name of HDF5 file to open for reading.
RETURNS: int id if ok, else -1 and a failure message is printed to stderr.
******************************************************************************/

int openDataset( const int file, const char* const datasetName ) {
  int result = 0;
  assert( file > -1 ); assert( datasetName ); assert( *datasetName );
  result = H5Dopen( file, datasetName, 0 );

  if ( result == -1 ) {
    fprintf( stderr, "\a\nFailed to open HDF5 dataset for reading: %s.\n",
             datasetName );
  }

  return result;
}



/******************************************************************************
PURPOSE: closeDataset - Close HDF5 dataset.
INPUTS:  const int dataset HDF5 dataset to close.
******************************************************************************/

void closeDataset( const int dataset ) {
  H5Dclose( dataset );
}



/******************************************************************************
PURPOSE: readFileAttribute - Read global dataset attribute.
INPUTS:  const int file               Id of file to read.
         const char* const attribute  Name of global attribute.
OUTPUTS: double* value                Value of attribute read.
RETURNS: int 1 if successfule, else 0 and failure message on stderr.
******************************************************************************/

int readFileAttribute( const int file, const char* const attribute,
                       double* value ) {
  int result = 0;
  int dataset = 0;

  assert( file > -1 ); assert( attribute ); assert( *attribute );
  assert( value );

  dataset = openDataset( file, attribute );
  result = dataset > -1;

  if ( result ) {
    const int status = H5Dread( dataset, H5T_NATIVE_DOUBLE_g, 0, 0, 0, value );
    result = status > -1;
    closeDataset( dataset ), dataset = -1;

    if ( ! result ) {
      fprintf( stderr, "\a\nFailed to read attribute '%s'.\n", attribute );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readDatasetDimensions - Read dataset dimensions.
INPUTS:  const int dataset   ID of dataset to query.
OUTPUTS: size_t* dimension0  1st dimension of dataset.
         size_t* dimension1  2nd dimension of dataset or 0 if rank 1.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readDatasetDimensions( const int dataset,
                           size_t* const dimension0,
                           size_t* const dimension1 ) {
  int result = 0;

  assert( dataset > -1 ); assert( dimension0 ); assert( dimension1 );

  {
    const int dataspace = H5Dget_space( dataset );

    if ( dataspace > -1 ) {
      unsigned long long dims[ 32 ] = { 0, 0 };
      const int rank = H5Sget_simple_extent_dims( dataspace, dims, 0 );

      if ( rank == 1 && dims[ 0 ] > 0 ) {
        *dimension0 = dims[ 0 ];
        *dimension1 = 0;
        result = 1;
      } else if ( rank == 2 && dims[ 0 ] > 0 && dims[ 1 ] > 0 ) {
        *dimension0 = dims[ 0 ];
        *dimension1 = dims[ 1 ];
        result = 1;
      }

      H5Sclose( dataspace );
    }
  }

  if ( ! result ) {
    fprintf( stderr, "\a\nFailed to read valid dimensions of dataset.\n" );
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileData - Read file data.
INPUTS:  const int dataset              ID of dataset to read.
         const size_t dimension0        1st dimension of data to read.
         const size_t dimension1        2nd dimension of data to read or 0.
OUTPUTS: double data[ dimension0 * dimension1 ]  Data read.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileData( const int dataset,
                  const size_t dimension0,
                  const size_t dimension1,
                  double data[] ) {

  int result = 0;

  assert( dataset > -1 ); assert( dimension0 );
  assert( data );

  result = matchedDatasetDimensions( dataset, dimension0, dimension1 );

  if ( result ) {
    result = H5Dread( dataset, H5T_NATIVE_DOUBLE_g, 0,0,0, data ) >= 0;

    if ( ! result ) {
      fprintf( stderr, "\a\nFailed to read matched file data.\n" );
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: readFileDataIntegers - Read file data.
INPUTS:  const int dataset              ID of dataset to read.
         const size_t dimension0        1st dimension of data to read.
         const size_t dimension1        2nd dimension of data to read or 0.
OUTPUTS: long long data[ dimension0 * dimension1 ]  Data read.
RETURNS: int 1 if ok, else 0 and a failure message is printed to stderr.
******************************************************************************/

int readFileDataIntegers( const int dataset,
                          const size_t dimension0,
                          const size_t dimension1,
                          long long data[] ) {

  int result = 0;

  assert( dataset > -1 ); assert( dimension0 );
  assert( data );

  result = matchedDatasetDimensions( dataset, dimension0, dimension1 );

  if ( result ) {
    result = H5Dread( dataset, H5T_NATIVE_LLONG_g, 0,0,0, data ) >= 0;

    if ( ! result ) {
      fprintf( stderr, "\a\nFailed to read matched file data.\n" );
    }
  }

  return result;
}



/*============================ PRIVATE FUNCTIONS ============================*/



/******************************************************************************
PURPOSE: matchedDatasetDimensions - Does the dataset have the expected
         dimensions?
INPUTS:  const int dataset       ID of dataset to check.
         const int dimension0    1st dimension.
         const int dimension1    2nd dimension or 0 if rank 1.
RETURNS: int 1 if matched, else 0 and a failure message is printed to stderr.
******************************************************************************/

static int matchedDatasetDimensions( const int dataset,
                                     const size_t dimension0,
                                     const size_t dimension1 ) {

  int result = 0;

  assert( dataset > -1 );
  assert( dimension0 > 0 || dimension1 > 0 );

  {
    const int dataspace = H5Dget_space( dataset );

    if ( dataspace < 0 ) {
      fprintf( stderr, "\a\n\nFailed to get dataspace for dataset %d\n",
               dataset );
    } else {
      const int rank = H5Sget_simple_extent_dims( dataspace, 0, 0 );
      const int expectedRank = ( dimension0 > 0 ) + ( dimension1 > 0 );

      if ( rank != expectedRank ) {
        fprintf( stderr, "\a\n\nInvalid/mismatched rank in dataset: "
                 "%d (expected %d)\n", rank, expectedRank );
      } else {
        unsigned long long dims[ 2 ] = { 0, 0 };
        const int status = H5Sget_simple_extent_dims( dataspace, dims, 0 );

        if ( status < 0 ) {
          fprintf( stderr, "\a\n\nFailed to get rank for dataset %d\n",
                   dataset );
        } else if ( dimension0 > 0 && dims[ 0 ] != dimension0 ) {
          fprintf( stderr, "\a\n\nInvalid/mismatched dims[0] in dataset: "
                   "%llu (expected %lu)\n", dims[ 0 ], dimension0 );
        } else if ( dimension1 > 0 && dims[ 1 ] != dimension1 ) {
          fprintf( stderr, "\a\n\nInvalid/mismatched dims[1] in dataset: "
                   "%llu (expected %lu)\n", dims[ 1 ], dimension1 );
        } else {
          result = 1;
        }
      }

      H5Sclose( dataspace );
    }
  }

  return result;
}



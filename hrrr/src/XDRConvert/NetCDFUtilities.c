
/******************************************************************************
PURPOSE: NetCDFUtilities.c - Define convenience routines for NetCDF files.

NOTES:   For a description of NetCDF COARDS Conventions see:
         http://ferret.wrc.noaa.gov/noaa_coop/coop_cdf_profile.html

HISTORY: 2007/12 plessel.todd@epa.gov, Created.
******************************************************************************/

/*================================ INCLUDES =================================*/

#ifdef DEBUGGING
#include <stdio.h>  /* For stderr, fprintf(). */
#endif
#include <string.h> /* For memset().  */

#include <netcdf.h> /* For nc_*. */

#include <Utilities.h> /* For PRE*(), NEW_ZERO(), Integer, Real, Stream. */

#include <Helpers.h>         /* For Name. */
#include <NetCDFUtilities.h> /* For public interface. */

/*================================ FUNCTIONS ================================*/



/******************************************************************************
PURPOSE: createLongitudeAndLatitude - Create longitude and latitude variables.
INPUTS:  Integer file                  NetCDF file to write to.
         Integer dimensionality        Number of dimensions of longitude.
         const Integer dimensionIds[]  NetCDF ids of dimensions.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

Integer createLongitudeAndLatitude( Integer file,
                                    Integer dimensionality,
                                    const Integer dimensionIds[] ) {

  PRE05( file != -1, dimensionality > 0, dimensionIds,
         dimensionIds[ 0 ] > -1, dimensionIds[ dimensionality - 1 ] > -1 );

  const Integer result =
    AND2( createVariable( file, "longitude", "degrees_east", NC_FLOAT, 1,
                          dimensionality, dimensionIds ) != -1,
          createVariable( file, "latitude", "degrees_north", NC_FLOAT, 1,
                          dimensionality, dimensionIds ) != -1 );

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeExtraAttributes - Write extra header variables and attributes.
INPUTS:  Integer file             NetCDF file to write to.
         const Real domain[2][2]  domain[LONGITUDE LATITUDE][MIMIMUM MAXIMUM].
         Integer dimensionId      NetCDF dimension id for variables.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

Integer writeExtraAttributes( Integer file, const Real domain[2][2],
                              Integer dimensionId ) {

  PRE07( file != -1, domain, dimensionId != -1,
         isValidLongitudeLatitude( domain[ LONGITUDE ][ MINIMUM ],
                                   domain[ LATITUDE  ][ MINIMUM ] ),
         isValidLongitudeLatitude( domain[ LONGITUDE ][ MAXIMUM ],
                                   domain[ LATITUDE  ][ MAXIMUM ] ),
         domain[ LONGITUDE ][ MINIMUM ] <= domain[ LONGITUDE ][ MAXIMUM ],
         domain[ LATITUDE  ][ MINIMUM ] <= domain[ LATITUDE  ][ MAXIMUM ] );

  Integer result = 0;

  if ( createVariable( file, "yyyyddd", "date",
                       NC_INT, 0, 1, &dimensionId ) != -1 ) {

    if ( createVariable( file, "hhmmss", "time",
                         NC_INT, 0, 1, &dimensionId ) != -1 ) {

      if ( writeRealAttribute( file, NC_GLOBAL, NC_FLOAT, "west_bound",
                               domain[ LONGITUDE ][ MINIMUM ] ) ) {

        if ( writeRealAttribute( file, NC_GLOBAL, NC_FLOAT, "east_bound",
                                 domain[ LONGITUDE ][ MAXIMUM ] ) ) {

          if ( writeRealAttribute( file, NC_GLOBAL, NC_FLOAT, "south_bound",
                                   domain[ LATITUDE ][ MINIMUM ] ) ) {

            if ( writeRealAttribute( file, NC_GLOBAL, NC_FLOAT, "north_bound",
                                    domain[ LATITUDE ][ MAXIMUM ] ) ) {
              result = 1;
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
PURPOSE: writeTimeData - Write time variables yyyyddd, hhmmss, time to file.
INPUTS:  Integer file               NetCDF file to write to.
         Integer count              Number of timestamps / dimensions.
         Integer dims               Number of dimensions in dimensions[][dims].
         Integer useBothDims        1 to multiply dimensions, 0 to skip 2nd.
         const Integer* timestamps  timestamps[ count ] yyyydddhhmm.
         const Integer* dimensions  dimensions[ count ][ dim ].
OUTPUTS: Real* buffer               Allocated storage for writing data.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

Integer writeTimeData( Integer file, Integer count, Integer dims,
                       Integer useBothDims,
                       const Integer* timestamps,
                       const Integer* dimensions,
                       Real* buffer ) {

  PRE013( file != -1, count > 0, IN3( dims, 1, 2 ),
          IS_BOOL( useBothDims ), IMPLIES( useBothDims, dims == 2 ),
          IMPLIES( dims == 1, useBothDims == 0 ),
          timestamps, dimensions, buffer,
          isValidTimestamp( timestamps[ 0 ] ),
          isValidTimestamp( timestamps[ count - 1 ] ),
          dimensions[ 0 ] > 0,
          dimensions[ ( count - 1 ) * dims + ( dims - 1 ) ] > 0 );

  Integer result = 0;
  const Integer yyyydddhhmmStart = timestamps[ 0 ];
  Integer dimension = 0;
  Integer start = 0;

  do {
    const Integer offset = dimension * dims;
    const Integer dataCount =
      useBothDims ? dimensions[ offset ] * dimensions[ offset + 1 ]
      : dimensions[ offset ];
    const Integer yyyydddhhmm = timestamps[ dimension ];
    const Integer yyyyddd = yyyydddhhmm / 10000;
    const Integer hhmmss  = yyyydddhhmm % 10000 * 100;
    int* const dateData = (int*) buffer; /* UGLY: Reuse storage. */
    int* const timeData = dateData + dataCount;
    assert_static( sizeof (int) == 4 ); assert_static( sizeof (Real) == 8 );
    replicateIntValue( dateData, dataCount, yyyyddd );
    replicateIntValue( timeData, dataCount, hhmmss );

    if ( ! AND2( writeSomeIntData( file, "yyyyddd", start,
                                   dataCount, 1, 1, 1, dateData ),
                 writeSomeIntData( file, "hhmmss", start,
                                   dataCount, 1, 1, 1, timeData ) ) ) {
      dimension = count;
    } else { /* Write float time array (hours from start time): */
      const Integer yyyydddhhmmNow = timestamps[ dimension ];
      const Real fractionalTime =
        fractionalHours( yyyydddhhmmStart, yyyydddhhmmNow );
      replicateRealValue( buffer, dataCount, fractionalTime );

      if ( ! writeSomeData( file, "time", start, dataCount, 1,1,1, buffer )) {
        dimension = count;
      }
    }

    start += dataCount;
    ++dimension;
  } while ( dimension < count );

  result = dimension == count;

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeTimeData1 - Write time variables yyyyddd, hhmmss, time to file.
INPUTS:  const Integer file          NetCDF file to write to.
         const Integer count         Number of timestamps / dimensions.
         const int yyyyddd[ count ]  Dates to write.
         const int hhmmss[  count ]  Times to write.
         const float fhour[ count ]  Fractional hours to write.
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
******************************************************************************/

Integer writeTimeData1( const Integer file, const Integer count,
                        const int yyyyddd[], const int hhmmss[],
                        const float fhour[] ) {

  PRE014( file != -1,
          count > 0,
          yyyyddd,
          isValidDate( yyyyddd[ 0 ] ),
          isValidDate( yyyyddd[ count - 1 ] ),
          yyyyddd[ count - 1 ] >= yyyyddd[ 0 ],
          hhmmss,
          isValidTime( hhmmss[ 0 ] ),
          isValidTime( hhmmss[ count - 1 ]),
          hhmmss[ count - 1 ] >= hhmmss[ 0 ],
          fhour,
          fhour[ 0 ] >= 0.0,
          fhour[ count - 1 ] >= fhour[ 0 ],
          fhour[ count - 1 ] >= fhour[ 0 ] );

  Integer result = 0;
  const char* variableName = "yyyyddd";
  int id = 0;
  int status = nc_inq_varid( file, variableName, &id );

  if ( status == NC_NOERR ) {
    const size_t start = 0;
    const size_t counts = count;
    status = nc_put_vara_int( file, id, &start, &counts, yyyyddd );

    if ( status == NC_NOERR ) {
      variableName = "hhmmss";
      status = nc_inq_varid( file, variableName, &id );

      if ( status == NC_NOERR ) {
        status = nc_put_vara_int( file, id, &start, &counts, hhmmss );

        if ( status == NC_NOERR ) {
          variableName = "time";
          status = nc_inq_varid( file, variableName, &id );

          if ( status == NC_NOERR ) {
            status = nc_put_vara_float( file, id, &start, &counts, fhour );
            result = status == NC_NOERR;
          }
        }
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't write variable '%s' because %s.",
                    variableName, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeStandardContents - Write standard COARDS contents.
INPUTS:  Integer file                  NetCDF file to write to.
         const char* history           Source/history of data (e.g., URL).
         const UTCTimestamp timestamp  YYYY-MM-DDTHH:MM:SS-0000.
         Integer timeDimension         Dimension of time.
         Integer timesteps             Number of timesteps.
         Integer writeTime             Write sequential time data values?
RETURNS: Integer 1 if successful, else 0 and failureMessage is called.
NOTES:   nc_enddef() is called so only data may be written to file now.
******************************************************************************/

Integer writeStandardContents( Integer file, const char* history,
                               const UTCTimestamp timestamp,
                               Integer timeDimension,
                               Integer timesteps, Integer writeTime ) {

  PRE09( file != -1, history, *history, strlen( history ) < 128,
         timestamp, *timestamp, timeDimension != -1,
         IS_BOOL( writeTime ),
         IMPLIES( writeTime, timesteps > 0 ) );

  Integer result = 0;

  if ( writeTextAttribute( file, NC_GLOBAL, "Conventions", "COARDS" ) ) {

    if ( writeTextAttribute( file, NC_GLOBAL, "crs", "latitude_longitude" ) ) {

      if ( writeTextAttribute( file, NC_GLOBAL, "history", history ) ) {
        Integer time = 0;
        char timeUnits[ 50 ] = "hours since ";
        CHECK( strlen( timeUnits ) + strlen( timestamp ) < 50 );
        strcat( timeUnits, timestamp );
        timeUnits[ 22 ] = ' ';
        strcpy( timeUnits + 31, ".0 -00:00" );
        CHECK( strlen( timeUnits ) < COUNT( timeUnits ) );
        time = createVariable( file, "time", timeUnits, NC_FLOAT, 0,
                               1, &timeDimension );

        if ( time != -1 ) {
          int status = nc_enddef( file );

          if ( status != NC_NOERR ) {
            const char* const message = nc_strerror( status );
            failureMessage( "Can't create file definition because %s.",
                            message );
          } else {
            int crsId = -1;
            int status = nc_inq_varid( file, "crs", &crsId );

            if ( status != NC_NOERR ) {
              const char* const message = nc_strerror( status );
              failureMessage( "Can't get int crs variable because %s.",
                              message );
            } else {
              const int crsValue = -9999;
              status = nc_put_var_int( file, crsId, &crsValue );

              if ( status != NC_NOERR ) {
                const char* const message = nc_strerror( status );
                failureMessage( "Can't write int crs variable because %s.",
                                message );
              } else {

                if ( ! writeTime ) {
                  result = 1;
                } else {
                  float* hours = NEW_ZERO( float, timesteps );

                  if ( hours ) {
                    const size_t start = 0;
                    const size_t count = timesteps;
                    Integer hour = 0;

                    do {
                      hours[ hour ] = hour;
                      ++hour;
                    } while ( hour < timesteps );

                    status = nc_put_vara_float( file, time, &start, &count, hours );
                    result = status == NC_NOERR;

                    if ( ! result ) {
                      const char* const message = nc_strerror( status );
                      failureMessage( "Can't write time variable because %s.",
                                      message );
                    }

                    FREE( hours );
                  }
                }
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
PURPOSE: createNetCDFFile - Create a NetCDFFile for writing.
INPUTS:  const char* fileName     Name of file to create.
         Integer create64BitFile  Create 64-bit NetCDF file?
RETURNS: Integer NetCDF file ID if successful, else -1 if failed and
         failureMessage() called.
******************************************************************************/

Integer createNetCDFFile( const char* fileName, Integer create64BitFile ) {

  PRE03( fileName, *fileName, IS_BOOL( create64BitFile ) );

  Integer result = -1;
  const int mode = create64BitFile ? NC_CLOBBER | NC_64BIT_OFFSET : NC_CLOBBER;
  int ncid = -1;
  const int status = nc_create( fileName, mode, &ncid );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't create file '%s' because %s.", fileName, message );
  } else if ( ncid < 0 ) {
    nc_close( ncid );
    failureMessage( "Invalid id for file '%s'.", fileName );
  } else {
    result = ncid;
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: createDimensions - Create the given named dimensions of a NetCDF file.
INPUTS:  Integer file              NetCDF file ID.
         Integer count             Number of dimensions.
         const char* const names[] Name of each dimension.
         const Integer values[]    Value of each dimension.
OUTPUTS: Integer ids[] >           Ids of each dimension or -1 if unsuccessful.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer createDimensions( Integer file, Integer count,
                          const char* const names[],
                          const Integer values[],
                          Integer ids[] ) {

  PRE09( file >= 0, count > 0, names, *names, names[ count - 1 ],
         values, *values > 0, values[ count - 1 ] > 0, ids );

  Integer result = 0;
  Integer index = 0;

  do {
    int id = -1;
    const int status =
      nc_def_dim( file, names[ index ], values[ index ], &id );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't create dimension %s because %s.",
                      names[ index ], message );
      index = count;
    } else {
      ids[ index ] = id;
    }

    ++index;
  } while ( index < count );

  result = index == count;

  if ( ! result ) {
    memset( ids, -1, count * sizeof *ids );
  }

  POST02( IS_BOOL( result ),
          IMPLIES_ELSE( result,
                        AND2( ids[ 0 ] >=  0, ids[ count - 1 ] >=  0 ),
                        AND2( ids[ 0 ] == -1, ids[ count - 1 ] == -1 ) ) );

  return result;
}



/******************************************************************************
PURPOSE: createCRSVariable - Create the crs variable of a NetCDF file
         which is needed so the resulting file can be displayed in
         visualization programs such as VERDI, QGIS, etc.
INPUTS:  const Integer file  NetCDF file ID.
RETURNS: Integer variable id > -1 if successful,
         else -1 and a failure message is printed to stderr.
NOTES:  The created variable will look like the following:
 int crs ;
		crs:spatial_ref = "GEOGCS[\"GCS_WGS_1984\",DATUM[\"WGS_1984\",\
                      SPHEROID[\"WGS_84\",6378137.0,298.257223563]],\
                      PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",\
                      0.017453292519943295]]" ;
		crs:grid_mapping_name = "latitude_longitude" ;
		crs:longitude_of_prime_meridian = 0. ;
		crs:semi_major_axis = 6378137. ;
		crs:inverse_flattening = 298.257223563 ;
See: https://gis.stackexchange.com/questions/218441/\
displaying-netcdf-data-with-correct-crs
******************************************************************************/

Integer createCRSVariable( const Integer file ) {

  PRE0( file >= 0 );

  Integer result = -1;
  int id = -1;
  int status = 0;

  DEBUG( fprintf( stderr, "createCRSVariable( %lld )\n", file ); )

  status = nc_def_var( file, "crs", NC_INT, 0, 0, &id );

  if ( status == NC_NOERR ) {
    const char* const spatial_ref =
      "GEOGCS[\"GCS_WGS_1984\",DATUM[\"WGS_1984\",SPHEROID[\"WGS_84\","
      "6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\","
      "0.017453292519943295]]";

    if ( writeTextAttribute( file, id, "spatial_ref", spatial_ref ) ) {

      if ( writeTextAttribute( file, id, "grid_mapping_name",
                               "latitude_longitude" ) ) {


        if ( writeRealAttribute( file, id, NC_DOUBLE,
                                 "longitude_of_prime_meridian", 0.0 ) ) {

          if ( writeRealAttribute( file, id, NC_DOUBLE,
                                   "semi_major_axis", 6378137.0 ) ) {

            if ( writeRealAttribute( file, id, NC_DOUBLE,
                                     "inverse_flattening", 298.257223563 ) ) {
              result = id;
            }
          }
        }
      }
    }
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't create variable createLongitudeAndLatitude because %s.", message );
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: createVariable - Create the given named variable of a NetCDF file.
INPUTS:  Integer file                 NetCDF file ID.
         const char* name             Name of variable.
         const char* units            Units of variable.
         Integer type                 Type of variable.
         Integer hasMissingValues     Can there be missing data values?
         Integer dimensionality       Number of dimensions.
         const Integer* dimensionIds  Ids of dimensions.
RETURNS: Integer variable id > -1 if successful, else -1 and failureMessage().
******************************************************************************/

Integer createVariable( Integer file, const char* name,
                        const char* units,
                        Integer type, Integer hasMissingValues,
                        Integer dimensionality,
                        const Integer* dimensionIds ) {

  PRE010( file >= 0, name, *name, units, *units,
          IN5( type, NC_INT, NC_FLOAT, NC_DOUBLE, NC_CHAR ),
          IS_BOOL( hasMissingValues ),
          dimensionality >= 1, dimensionIds,
          OR2( dimensionIds[ 0 ] == NC_UNLIMITED, dimensionIds[ 0 ] >= 0 ) );

  Integer result = -1;
  int id = -1;
  int ids[ NC_MAX_DIMS ];
  int status = 0;
  Integer dimension = 0;

  DEBUG( fprintf( stderr, "createVariable( %lld, '%s', '%s', "
                  "%lld, %lld, %lld, %lld )\n",
                  file, name, units, type, hasMissingValues, dimensionality,
                  dimensionIds[ 0 ] ) );

  do {
    ids[ dimension ] = dimensionIds[ dimension ];
    ++dimension;
  } while ( dimension < dimensionality );

  status = nc_def_var( file, name, type, dimensionality, ids, &id );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't create variable %s because %s.", name, message );
  } else {
    const char* const renamedUnits = strcmp( units, "-" ) ? units : "none";

    if ( writeTextAttribute( file, id, "units", renamedUnits ) ) {

      if ( hasMissingValues ) {

        if ( writeRealAttribute( file, id, NC_FLOAT, "missing_value",
                                 -9999.0 ) ) {
          if  ( writeTextAttribute( file, id, "grid_mapping", "crs" ) ) {
            result = id;
          }
        }
      } else {
        result = id;
      }
    }
  }

  POST0( result >= -1 );
  return result;
}



/******************************************************************************
PURPOSE: writeIntegerAttribute - Write the value of a global integer attribute
         to a NetCDF file.
INPUTS:  Integer file      NetCDF file ID.
         const char* name  Name of attribute.
         Integer value     Value for attribute.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer writeIntegerAttribute( Integer file, const char* name,
                               Integer value ) {

  PRE03( file >= 0, name, *name );

  const int attribute = CLAMPED_TO_RANGE( value, INT_MIN, INT_MAX );
  int status = nc_put_att_int( file, NC_GLOBAL, name, NC_INT, 1, &attribute );
  const Integer result = status == NC_NOERR;

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't write value of attribute %s because %s.",
                    name, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRealAttribute - Write the value of a real attribute to a
         NetCDF file.
INPUTS:  Integer file      NetCDF file ID.
         Integer id        NetCDF variable ID or NC_GLOBAL.
         Integer type      NC_FLOAT or NC_DOUBLE.
         const char* name  Name of attribute.
         Real value        Value for attribute.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer writeRealAttribute( Integer file, Integer id, Integer type,
                            const char* name, Real value ) {

  PRE06( file >= 0, OR2( id == NC_GLOBAL, id >= 0 ),
         IN3( type, NC_FLOAT, NC_DOUBLE ),
         name, *name, ! isNan( value ) );

  Integer result = 0;
  int status = 0;

  if ( type == NC_FLOAT ) {
    const float attribute = CLAMPED_TO_RANGE( value, -FLT_MAX, FLT_MAX );
    status = nc_put_att_float( file, id, name, NC_FLOAT, 1, &attribute );
  } else {
    const double attribute = value;
    status = nc_put_att_double( file, id, name, NC_DOUBLE, 1, &attribute );
  }

  result = status == NC_NOERR;

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't write value of attribute %s because %s.",
                    name, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeTextAttribute - Write the value of a text attribute.
INPUTS:  Integer file      NetCDF file ID.
         Integer id        NetCDF variable ID or NC_GLOBAL.
         const char* name  Name of attribute.
         const char* value Value of named attribute.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer writeTextAttribute( Integer file, Integer id, const char* name,
                            const char* value ) {

  PRE06( file >= 0, OR2( id == NC_GLOBAL, id >= 0 ),
         name, *name, value, *value );

  const int status = nc_put_att_text( file, id, name, strlen( value ), value );
  const Integer result = status == NC_NOERR;

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage("Can't write text attribute %s because %s.", name,message);
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeRealArrayAttribute - Write a real array attribute.
INPUTS:  Integer file       NetCDF file ID to write to.
         Integer type       NC_FLOAT or NC_DOUBLE.
         const char* name   Name of attribute.
         Integer count      Number of elements in values array.
         const Real* values Value of attribute.
RETURNS: Integer 1 if successful, else 0 and failureMessage() was called.
******************************************************************************/

Integer writeRealArrayAttribute( Integer file, Integer type,
                                 const char* name,
                                 const Real* values, Integer count ) {

  PRE07( file >= 0, IN3( type, NC_FLOAT, NC_DOUBLE ),
         name, *name, values, count > 0, isNanFree( values, count ) );

  Integer result = 0;
  int status = 0;

  if ( type == NC_FLOAT ) {
    float* attributes = NEW_ZERO( float, count );

    if ( attributes ) {
      Integer index = 0;

      for ( index = 0; index < count; ++index ) {
        const Real value = values[ index ];
        const float attribute = CLAMPED_TO_RANGE( value, -FLT_MAX, FLT_MAX );
        attributes[ index ] = attribute;
      }

      status =
        nc_put_att_float( file, NC_GLOBAL, name, NC_FLOAT, count, attributes);

      result = status == NC_NOERR;
      FREE_ZERO( attributes );
    }
  } else {
    status =
      nc_put_att_double( file, NC_GLOBAL, name, NC_DOUBLE, count, values );
  }

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't write value of attribute %s because %s.",
                    name, message );
  }

  POST0( IS_BOOL( result ) );
  return result;
}



/******************************************************************************
PURPOSE: writeAllData - Write all data to the NetCDF file.
INPUTS:  Integer file              NetCDF file handle.
         const char* variableName  Name of variable to write.
         Integer dimension1        E.g., timesteps.
         Integer dimension2        E.g., layers or 1.
         Integer dimension3        E.g., rows or 1.
         Integer dimension4        E.g., columns or 1.
         Real* data                Data to write.
OUTPUTS: Real* data                Ruined data...
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeAllData( Integer file, const char* variableName,
                      Integer dimension1, Integer dimension2,
                      Integer dimension3, Integer dimension4,
                      Real* data ) {

  PRE06( file > -1, variableName, *variableName, dimension1 > 0,
          GT_ZERO3( dimension2, dimension3, dimension4 ), data );

  Integer result = 0;
  int id = -1;
  int status = nc_inq_varid( file, variableName, &id );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    variableName, message );
  } else {
    const size_t start[ 4 ] = { 0, 0, 0, 0 };
    size_t count[ 4 ] = { 0, 0, 0, 0 };
    float* const fdata = (float*) data; /* UGLY const cast away, avoid copy!*/
    count[ 0 ] = dimension1;
    count[ 1 ] = dimension2;
    count[ 2 ] = dimension3;
    count[ 3 ] = dimension4;
    compress64BitValues( data,
                         dimension4 * dimension3 * dimension2 * dimension1 );
    status = nc_put_vara_float( file, id, start, count, fdata );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't write subset of variable '%s' because %s.",
                      variableName, message );
    } else {
      result = 1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: writeAllIntData - Write all data to the NetCDF file.
INPUTS:  Integer file              NetCDF file handle.
         const char* variableName  Name of variable to write.
         Integer dimension1        E.g., timesteps.
         Integer dimension2        E.g., layers or 1.
         Integer dimension3        E.g., rows or 1.
         Integer dimension4        E.g., columns or 1.
         Integer* data             Data to write.
OUTPUTS: Integer* data             Ruined data...
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeAllIntData( Integer file, const char* variableName,
                         Integer dimension1, Integer dimension2,
                         Integer dimension3, Integer dimension4,
                         Integer* data ) {

  PRE06( file > -1, variableName, *variableName, dimension1 > 0,
          GT_ZERO3( dimension2, dimension3, dimension4 ), data );

  Integer result = 0;
  int id = -1;
  int status = nc_inq_varid( file, variableName, &id );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    variableName, message );
  } else {
    const size_t start[ 4 ] = { 0, 0, 0, 0 };
    size_t count[ 4 ] = { 0, 0, 0, 0 };
    int* const idata = (int*) data; /* UGLY const cast away, avoid copy!*/
    count[ 0 ] = dimension1;
    count[ 1 ] = dimension2;
    count[ 2 ] = dimension3;
    count[ 3 ] = dimension4;
    compress64BitIntegerValues( data,
                                dimension4 * dimension3 * dimension2 *
                                dimension1 );
    status = nc_put_vara_int( file, id, start, count, idata );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't write subset of variable '%s' because %s.",
                      variableName, message );
    } else {
      result = 1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: writeAllCharData - Write all character data to the NetCDF file.
INPUTS:  Integer file              NetCDF file handle.
         const char* variableName  Name of variable to write.
         Integer count             Number of items in data[].
         Integer length            Length of each item in data[].
         const char data[ count ] Data to write: data[ count * length ].
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeAllCharData( Integer file, const char* variableName,
                          Integer count, Integer length,
                          const char data[] ) {

  PRE07( file > -1, variableName, *variableName,
         count > 0, length > 0, length < 1000, data );

  Integer result = 0;
  int id = -1;
  int status = nc_inq_varid( file, variableName, &id );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    variableName, message );
  } else {
    const size_t starts[ 2 ] = { 0, 0 };
    size_t counts[ 2 ] = { 0, 0 };
    counts[ 0 ] = count;
    counts[ 1 ] = length;
    status = nc_put_vara_text( file, id, starts, counts, data );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't write subset of variable '%s' because %s.",
                      variableName, message );
    } else {
      result = 1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: writeSomeData - Write some data to the NetCDF file.
INPUTS:  Integer file              NetCDF file handle.
         const char* variableName  Name of variable to write.
         Integer timestep          Timestep to write [ 0, dimension1 - 1 ].
         Integer dimension1        E.g., timesteps or 1.
         Integer dimension2        E.g., layers or 1.
         Integer dimension3        E.g., rows or 1.
         Integer dimension4        E.g., columns or 1.
         Real* data                Data to write.
OUTPUTS: Real* data                Ruined data...
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeSomeData( Integer file, const char* variableName,
                       Integer timestep,
                       Integer dimension1, Integer dimension2,
                       Integer dimension3, Integer dimension4,
                       Real* data ) {

  PRE06( file > -1, variableName, *variableName, timestep >= 0,
          GT_ZERO4( dimension1, dimension2, dimension3, dimension4 ), data );

  Integer result = 0;
  int id = -1;
  int status = nc_inq_varid( file, variableName, &id );

  DEBUG( fprintf( stderr, "writeSomeData( %lld, '%s', "
                  "%lld, %lld, %lld, %lld, %lld, [%lf ...] )\n",
                  file, variableName, timestep,
                  dimension1, dimension2, dimension3, dimension4,
                  data[ 0 ] ) );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    variableName, message );
  } else {
    size_t start[ 4 ] = { 0, 0, 0, 0 };
    size_t count[ 4 ] = { 0, 0, 0, 0 };
    float* const fdata = (float*) data; /* UGLY const cast away, avoid copy!*/
    const Integer size = dimension4 * dimension3 * dimension2 * dimension1;
    start[ 0 ] = timestep;
    count[ 0 ] = dimension1;
    count[ 1 ] = dimension2;
    count[ 2 ] = dimension3;
    count[ 3 ] = dimension4;
    compress64BitValues( data, size );
    status = nc_put_vara_float( file, id, start, count, fdata );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't write subset of variable '%s' because %s.",
                      variableName, message );
    } else {
      result = 1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: writeSomeIntData - Write some int data to the NetCDF file.
INPUTS:  Integer file              NetCDF file handle.
         const char* variableName  Name of variable to write.
         Integer timestep          Timestep to write [ 0, dimension1 - 1 ].
         Integer dimension1        E.g., timesteps or 1.
         Integer dimension2        E.g., layers or 1.
         Integer dimension3        E.g., rows or 1.
         Integer dimension4        E.g., columns or 1.
         const int* data                Data to write.
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeSomeIntData( Integer file, const char* variableName,
                          Integer timestep,
                          Integer dimension1, Integer dimension2,
                          Integer dimension3, Integer dimension4,
                          const int* data ) {

  PRE06( file > -1, variableName, *variableName, timestep >= 0,
          GT_ZERO4( dimension1, dimension2, dimension3, dimension4 ), data );

  Integer result = 0;
  int id = -1;
  int status = nc_inq_varid( file, variableName, &id );

  DEBUG( fprintf( stderr, "writeSomeIntData( %lld, '%s', "
                  "%lld, %lld, %lld, %lld, %lld, [%d ...] )\n",
                  file, variableName, timestep,
                  dimension1, dimension2, dimension3, dimension4,
                  data[ 0 ] ) );

  if ( status != NC_NOERR ) {
    const char* const message = nc_strerror( status );
    failureMessage( "Can't determine id of variable '%s' because %s.",
                    variableName, message );
  } else {
    size_t start[ 4 ] = { 0, 0, 0, 0 };
    size_t count[ 4 ] = { 0, 0, 0, 0 };
    start[ 0 ] = timestep;
    count[ 0 ] = dimension1;
    count[ 1 ] = dimension2;
    count[ 2 ] = dimension3;
    count[ 3 ] = dimension4;
    status = nc_put_vara_int( file, id, start, count, data );

    if ( status != NC_NOERR ) {
      const char* const message = nc_strerror( status );
      failureMessage( "Can't write subset of variable '%s' because %s.",
                      variableName, message );
    } else {
      result = 1;
    }
  }

  return result;
}



/******************************************************************************
PURPOSE: writeSomeIntegerData - Write some integer data to the NetCDF file.
INPUTS:  Integer file              NetCDF file handle.
         const char* variableName  Name of variable to write.
         Integer timestep          Timestep to write [ 0, dimension1 - 1 ].
         Integer dimension1        E.g., timesteps or 1.
         Integer dimension2        E.g., layers or 1.
         Integer dimension3        E.g., rows or 1.
         Integer dimension4        E.g., columns or 1.
         Integer* data             Data to write.
OUTPUTS: Real* data                Ruined data...
RETURNS: Integer 1 if unsuccessful, else 0 and failureMessage() is called.
******************************************************************************/

Integer writeSomeIntegerData( Integer file, const char* variableName,
                              Integer timestep,
                              Integer dimension1, Integer dimension2,
                              Integer dimension3, Integer dimension4,
                              Integer* data ) {

  PRE06( file > -1, variableName, *variableName, timestep >= 0,
          GT_ZERO4( dimension1, dimension2, dimension3, dimension4 ), data );

  Integer result = 0;

  compress64BitIntegerValues( data,
                              dimension4 * dimension3 * dimension2 *
                              dimension1 );

  result =
    writeSomeIntData( file, variableName, timestep,
                      dimension1,  dimension2, dimension3,  dimension4,
                      (const int*) data );

  return result;
}





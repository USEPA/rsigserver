#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: erddap.pm - NOAA ERDDAP coverage metadata used by erddapserver.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use erddap; # Use my $coverage_metadata_reference = erddap::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-05-06 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package erddap;

##############################################################################


# Metadata:
# name         - The COVERAGE= name used by this script.
# units        - Units for the parameter, e.g., mg/m3.
# description  - Parameter description.
# start_date   - Earliest date of available data.
# end_date     - Latest date of available data.
#
# min          - Minimum valid value.
# max          - Maximum valid value.
# previous_day - 1 to query previous day, else 0 for current day.
# query        - template for remote query.
# negative_lon - 1 if query longitude range is [-180, 180], 0 if [0, 360].
# swap_rows    - 1 if query latitude  range is [MAX_LAT, MIN_LAT], else 0.
# file_format       - Format of data returned by the remote web service: grid or csv
# header_lines - Number of result header lines to skip.
# fields       - Number of result data fields separated by commas.
# timestamp    - 0-based index of timestamp in result line.
# longitude    - 0-based index of longitude in result line.
# latitude     - 0-based index of latitude  in result line.
# depth        - 0-based index of depth     in result line or -1 if none/ignore
# id           - 0-based index of site-id   in result line or -1 if none/ignore
# sid          - 0-based index of site-name in result line or -1 if none/ignore
# value        - 0-based index of measure   in result line.
#
# For a list of all potentially available servers/datasets:
# https://coastwatch.pfeg.noaa.gov/erddap/info/index.html?page=1&itemsPerPage=5000

my @coverage_metadata  = (

  # Gridded:

  # MODIS 12-line DODS header looks like this:
  # Dataset {
  #  GRID {
  #    ARRAY:
  #      Float32 chlorophyll[time = 1][altitude = 1][latitude = 241][longitude = 241];
  #    MAPS:
  #      Float64 time[time = 1];
  #      Float64 altitude[altitude = 1];
  #      Float64 latitude[latitude = 241];
  #      Float64 longitude[longitude = 241];
  #  } chlorophyll;
  #} erdMHchla1day;
  #
  {
    name => 'MODIS.daytime_sea_surface_temperature',
    units => 'C',
    description =>
      'MODIS Aqua satellite-measured daytime daily composite 0.04166 degree' .
      'sea surface temperature.',
    begin_date => '2019-01-01',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdMH1sstd1day_R2022NRTMasked.dods?' .
               'sstMasked' .
               '[(YYYY1-MM1-DD1T00:00:00Z):1:(YYYY2-MM2-DD2T23:59:59Z)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    previous_day => 1,
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },
  {
    name => 'MODIS.attenuation_coefficient',
    units => '/m', # range: [0, 5].
    description =>
      'MODIS Aqua satellite-measured daily composite 4km' .
      'diffuse attenuation coefficient K measured at 490nm wavelength.',
    begin_date => '2003-01-01',
    end_date   => '',
    min => 0.0,
    max => 10.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdMH1kd4901day_R2022NRT.dods?' .
               'Kd_490' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               #'[(0.0):1:(0.0)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },
  {
    name => 'MODIS.atlantic.fluorescence',
    units => 'W/m2/s2/um/sr',
    description =>
      'MODIS Aqua satellite-measured daily composite 4km fluorescence.',
    begin_date => '2003-01-01',
    end_date   => '',
    min => 0.0,
    max => 5.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdMH1cflh1day_R2022NRT.dods?' .
               'nflh' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 0,
    file_format  => 'dods'
  },
  {
    name => 'MODIS.gulf.fluorescence',
    units => 'W/m2/s2/um/sr',
    description =>
      'MODIS Aqua satellite-measured daily composite 4km fluorescence.',
    begin_date => '2003-01-01',
    end_date   => '',
    min => 0.0,
    max => 5.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdMH1cflh1day_R2022NRT.dods?' .
               'nflh' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 0,
    file_format  => 'dods'
  },
  {
    name => 'MODIS.pacific.fluorescence',
    units => 'W/m2/s2/um/sr',
    description =>
      'MODIS Aqua satellite-measured daily composite 4km fluorescence.',
    begin_date => '2003-01-01',
    end_date   => '',
    min => 0.0,
    max => 5.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdMH1cflh1day_R2022NRT.dods?' .
               'nflh' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 0,
    file_format  => 'dods'
  },
  {
    name => 'MODIS.chlorophyll',
    units => 'mg/m3',
    description =>
      'MODIS Aqua satellite-measured daily composite chlorophyll-a.',
    begin_date => '2003-09-29',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdMH1chla1day_R2022NRT.dods?' .
               'chlorophyll' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },
  # CDOM is only available for West domain. No Atlantic or Gulf versions.
  {
    name => 'MODIS.cdom',
    units => '/m', # Light absorbtion coefficient range: [0, 5].
    description =>
      'MODIS Aqua satellite-measured daily composite 0.0125-degree ' .
      'chromophoric dissolved organic material.',
    begin_date => '2012-03-11',
    end_date   => '2016-07-18',
    min => 0.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               #'griddap/erdMEcdommday_Lon0360.dods?' .
               'griddap/erdMWcdom1day.dods?' .
               'cdom' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(0.0):1:(0.0)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 0,
    swap_rows    => 0,
    file_format  => 'dods'
  },

  # SeaWiFS 12-line DODS header looks like this:
  # Dataset {
  #  GRID {
  #    ARRAY:
  #      Float32 chlorophyll[time = 1][altitude = 1][latitude = 241][longitude = 241];
  #    MAPS:
  #      Float64 time[time = 1];
  #      Float64 altitude[altitude = 1];
  #      Float64 latitude[latitude = 241];
  #      Float64 longitude[longitude = 241];
  #  } chlorophyll;
  #} erdSWchla1day;
  #
  # erdSW1chla1day deprecated and superceded by erdSW2018chla1day.
  {
    name => 'SEAWIFS.chlorophyll',
    units => 'mg/m3',
    description =>
      'SeaWiFS Orbview-2 satellite-measured daily composite 0.1-degree chlorophyll-a.',
    begin_date => '1997-01-11',
    end_date   => '2010-12-31',
    min => 0.0,
    max => 1000.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdSW2018chla1day.dods?' .
               'chlorophyll' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },

  # VIIRS Pacific 750m 13-line DODS header looks like this:
  # Dataset {
  #   GRID {
  #     ARRAY:
  #       Float32 chla[time = 1][altitude = 1][latitude = 11985][longitude = 9338];
  #     MAPS:
  #       Float64 time[time = 1];
  #       Float64 altitude[altitude = 1];
  #       Float64 latitude[latitude = 11985];
  #       Float64 longitude[longitude = 9338];
  #   } chla;
  # } erdVHNchla1day;
  #
  # Data:
  #
  {
    name => 'VIIRS.pacific.chlorophyll',
    units => 'mg/m3',
    description =>
      'VIIRS North Pacific 750m satellite-measured daily composite chlorophyll-a.',
    begin_date => '2015-02-25',
    end_date   => '',
    min => 0.0,
    max => 1000.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVHNchla1day.dods?' .
               'chla' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(0.0):1:(0.0)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },

  {
    name => 'VIIRS.pacific.chlorophyll_bloom',
    units => 'mg/m3',
    description =>
      'VIIRS North Pacific 750m satellite-measured daily composite chlorophyll-a.',
    begin_date => '2015-02-25',
    end_date   => '',
    min => 0.0,
    max => 1000.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVHNchla1day.dods?' .
               'chla' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(0.0):1:(0.0)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },


  # VIIRS 4km Chlorophyll 12-line DODS header looks like this:
  # Dataset {
  #   GRID {
  #     ARRAY:
  #       Float32 chla[time = 1][latitude = 4320][longitude = 8640];
  #     MAPS:
  #       Float64 time[time = 1];
  #       Float32 latitude[latitude = 4320];
  #       Float32 longitude[longitude = 8640];
  #   } chla;
  # } erdVH2chla1day;
  #
  # Data:
  #
  # erdVH2chla1day deprecated and superceded by erdVH2018chla1day.
  {
    name => 'VIIRS.chlorophyll',
    units => 'mg/m3',
    description =>
      'VIIRS 4km satellite-measured daily composite chlorophyll-a.',
    begin_date => '2012-01-02',
    end_date   => '2022-07-26',
    min => 0.0,
    max => 1000.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVH2018chla1day.dods?' .
               'chla' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },


  # VIIRS 4km K490 12-line DODS header looks like this:
  # Dataset {
  #   GRID {
  #     ARRAY:
  #       Float32 k490[time = 1][latitude = 4320][longitude = 8640];
  #     MAPS:
  #       Float64 time[time = 1];
  #       Float32 latitude[latitude = 4320];
  #       Float32 longitude[longitude = 8640];
  #   } k490;
  # } erdVH2k4901day;
  #
  # Data:
  #
  {
    name => 'VIIRS.k490',
    units => '/m',
    description =>
      'VIIRS 4km satellite-measured daily composite ' .
      'diffuse attenuation coefficient in downwelling radiative flux ' .
      'in sea water.',
    begin_date => '2012-01-02',
    end_date   => '2022-07-26',
    min => 0.0,
    max => 0.5,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVH2018k4901day.dods?' .
               'k490' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },


  # VIIRS 4km PAR 12-line DODS header looks like this:
  # Dataset {
  #   GRID {
  #     ARRAY:
  #       Float32 par[time = 1][latitude = 4320][longitude = 8640];
  #     MAPS:
  #       Float64 time[time = 1];
  #       Float32 latitude[latitude = 4320];
  #       Float32 longitude[longitude = 8640];
  #   } par;
  # } erdVH2par1day;
  #
  # Data:
  #
  # erdVH2par1day deprecated and superceded by erdVH2018par1day.
  {
    name => 'VIIRS.par',
    units => 'Einsteins/m2/dy',
    description =>
      'VIIRS 4km satellite-measured daily composite ' .
      'downwelling photosynthetic photon radiance in sea water.',
    begin_date => '2012-01-02',
    end_date   => '2022-07-26',
    min => 0.0,
    max => 130.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVH2018par1day.dods?' .
               'par' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },


  # VIIRS 4km PIC 12-line DODS header looks like this:
  # Dataset {
  #   GRID {
  #     ARRAY:
  #       Float32 pic[time = 1][latitude = 4320][longitude = 8640];
  #     MAPS:
  #       Float64 time[time = 1];
  #       Float32 latitude[latitude = 4320];
  #       Float32 longitude[longitude = 8640];
  #   } pic;
  # } erdVH2pic1day;
  #
  # Data:
  #
  # erdVH2pic1day deprecated and superceded by erdVH2018pic1day.
  {
    name => 'VIIRS.pic',
    units => 'mol/m3',
    description =>
      'VIIRS 4km satellite-measured daily composite ' .
      'particulate inorganic carbon.',
    begin_date => '2012-01-02',
    end_date   => '2022-07-26',
    min => 0.0,
    max => 0.125,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVH2018pic1day.dods?' .
               'pic' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },


  # VIIRS 4km POC 12-line DODS header looks like this:
  # Dataset {
  #   GRID {
  #     ARRAY:
  #       Float32 poc[time = 1][latitude = 4320][longitude = 8640];
  #     MAPS:
  #       Float64 time[time = 1];
  #       Float32 latitude[latitude = 4320];
  #       Float32 longitude[longitude = 8640];
  #   } poc;
  # } erdVH2poc1day;
  #
  # Data:
  #
  # erdVH2poc1day deprecated and superceded by erdVH2018poc1day.
  {
    name => 'VIIRS.poc',
    units => 'mg/m3',
    description =>
      'VIIRS 4km satellite-measured daily composite ' .
      'particulate organic carbon.',
    begin_date => '2012-01-02',
    end_date   => '2022-07-26',
    min => 0.0,
    max => 1000.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVH2018poc1day.dods?' .
               'poc' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },


  # VIIRS 4km reflectance 12-line DODS header looks like this:
  # Dataset {
  #   GRID {
  #     ARRAY:
  #       Float32 r671[time = 1][latitude = 4320][longitude = 8640];
  #     MAPS:
  #       Float64 time[time = 1];
  #       Float32 latitude[latitude = 4320];
  #       Float32 longitude[longitude = 8640];
  #   } r671;
  # } erdVH2r6711day;
  #
  # Data:
  #
  # erdVH2r6711day deprecated and superceded by erdVH2018r6711day.
  {
    name => 'VIIRS.reflectance',
    units => '/sr',
    description =>
      'VIIRS 4km satellite-measured daily composite ' .
      'surface ratio of upwelling radiance emerging from sea water to ' .
      'downwelling radiative flux in air at 671nm wavelength.',
    begin_date => '2012-01-02',
    end_date   => '2022-07-26',
    min => 0.0,
    max => 0.004,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/erdVH2018r6711day.dods?' .
               'r671' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MAX_LAT):1:(MIN_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 1,
    file_format  => 'dods'
  },


  # MUR SST 12-line DODS header looks like this:
  #Dataset {
  #  GRID {
  #    ARRAY:
  #      Float64 analysed_sst[time = 3][latitude = 92][longitude = 92];
  #    MAPS:
  #      Float64 time[time = 3];
  #      Float32 latitude[latitude = 92];
  #      Float32 longitude[longitude = 92];
  #  } analysed_sst;
  #} jplMURSST;
  #
  # jplMURSST deprecated and superceded by jplMURSST41.
  {
    # DODS binary file_format, Missing data value = -7.768
    name => 'MUR.water_temperature',
    units => 'C',
    description =>
      'JPL Multi-scale Ultra-high Resolution 1km daily ' .
      'sea surface temperature analysis.',
    begin_date => '2002-01-01',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'griddap/jplMURSST41.dods?' .
               'analysed_sst' .
               '[(YYYY1-MM1-DD1T12:00:00Z):1:(YYYY2-MM2-DD2T12:00:00Z)]' .
               '[(MIN_LAT):1:(MAX_LAT)]' .
               '[(MIN_LON):1:(MAX_LON)]',
    negative_lon => 1,
    swap_rows    => 0,
    file_format  => 'dods'
  },



  ########################################################################
  # Non-gridded: buoy stations ###########################################
  # Note as of 2022 the erddap webservice now requires encoding of < and >
  # as %3C and %3E.
  # https://coastwatch.pfeg.noaa.gov/erddap/tabledap/documentation.html#PercentEncoded
  #
  # GTSPP data looks like this:
  #time,station_id,longitude,latitude,depth,temperature,salinity
  #UTC,,degrees_east,degrees_north,m,degree_C,PSU
  #2006-07-18T00:00:00Z,3259901,-80.862,24.843,1.0,31.02,35.0
  #2006-07-18T01:00:00Z,3259902,-80.862,24.843,1.0,30.97,34.86
  #...
  #2006-07-20T22:00:00Z,3262749,-80.862,24.843,1.0,29.75,34.8
  #2006-07-20T23:00:00Z,3262750,-80.862,24.843,1.0,29.76,34.9
  #2006-07-21T00:00:00Z,3262751,-80.862,24.843,1.0,29.71,34.97

  {
    name => 'GTSPP.station',
    units => '-',
    description => 'GTSPP station id.',
    begin_date => '1985-02-15',
    end_date   => '',
    min => '',
    max => '',
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/erdGtsppBest.csv?' .
      'station_id,longitude,latitude' .
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 3,
    timestamp    => -1,
    longitude    => 1,
    latitude     => 2,
    depth        => -1,
    id           => 0,
    sid          => -1,
    value        => -1
  },
  {
    name => 'GTSPP.water_temperature',
    units => 'C',
    description => 'GTSPP measured hourly water temperature.',
    begin_date => '1985-02-15',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/erdGtsppBest.csv?' .
               'time,station_id,longitude,latitude,depth,temperature'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => 1,
    sid          => -1,
    value        => 5
  },
  {
    name => 'GTSPP.salinity',
    units => 'PSU',
    description => 'GTSPP measured hourly salinity.',
    begin_date => '1985-02-15',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/erdGtsppBest.csv?' .
               'time,station_id,longitude,latitude,depth,salinity'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => 1,
    sid          => -1,
    value        => 5
  },
  {
    name => 'GTSPP.all',
    units => '',
    description => 'GTSPP measured hourly water temperature and salinity.',
    begin_date => '1985-02-15',
    end_date   => '',
    minimums => [
      0.0, 0.0, 0.0, 0.0, # Unused: timestamp, station, longitude, latitude.
      -500.0,
      -50.0,
      0.0
    ],
    maximums => [
      0.0, 0.0, 0.0, 0.0, # Unused: timestamp, station, longitude, latitude.
      0.0,
      50.0,
      50.0
    ],
    factors => [
      0.0, 0.0, 0.0, 0.0, # Unused: timestamp, station, longitude, latitude.
      1.0,
      1.0,
      1.0
    ],
    output_header =>
      "\twater_temperature(C)" .
      "\tsalinity(PSU)\n",
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/erdGtsppBest.csv?' .
               'time,station_id,longitude,latitude,depth,temperature,salinity'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 7,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => 1,
    sid          => -1,
    value        => 5
  },


  # NDBC buoy data looks like this:
  #time,station,longitude,latitude,wd,wspd,gst,wvht,dpd,apd,mwd,bar,atmp,wtmp,dewp,vis,ptdy,tide,wspu,wspv
  #UTC,,degrees_east,degrees_north,degrees_true,m s-1,m s-1,m,s,s,degrees_true,hPa,degree_C,degree_C,degree_C,km,hPa,m,m s-1,m s-1
  #2006-07-18T00:00:00Z,41001,-72.698,34.675,91,1.7,2.7,1.09,8.33,5.51,NaN,1015.1,26.6,26.6,NaN,NaN,NaN,NaN,-1.7,0.0
  #2006-07-18T01:00:00Z,41001,-72.698,34.675,102,1.9,2.4,1.09,6.25,5.29,NaN,1015.7,26.6,26.6,NaN,NaN,NaN,NaN,-1.9,0.4
  #...
  #2006-07-18T23:00:00Z,YKTV2,-76.478,37.227,100,3.6,4.0,NaN,NaN,NaN,NaN,1014.5,31.3,27.8,NaN,NaN,NaN,NaN,-3.5,0.6
  #2006-07-19T00:00:00Z,YKTV2,-76.478,37.227,132,3.5,4.8,NaN,NaN,NaN,NaN,1014.7,31.3,27.8,NaN,NaN,NaN,NaN,-2.6,2.3
  {
    name => 'NDBC.station',
    units => '-',
    description => 'NDBC buoy id.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => '',
    max => '',
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'station,longitude,latitude' .
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 3,
    timestamp    => -1,
    longitude    => 1,
    latitude     => 2,
    depth        => -1,
    id           => -1,
    sid          => 0,
    value        => -1
  },
  {
    name => 'NDBC.air_temperature',
    units => 'C',
    description => 'NDBC buoy-measured hourly air temperature.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => -50.0,
    max => 50.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,atmp'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.dew_point',
    units => 'C',
    description => 'NDBC buoy-measured hourly dew point temperature.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => -50.0,
    max => 50.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,dewp'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.pressure',
    units => 'hPa',
    description => 'NDBC buoy-measured hourly air pressure.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 500.0,
    max => 1500.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,bar'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.pressure_tendency',
    units => 'hPa',
    description =>
      'NDBC buoy-measured hourly air pressure tendency (plus or minus) ' .
      'for a three hour period ending at the time of measurement.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => -500.0,
    max => 500.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,ptdy'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.air_visibility',
    units => 'm',
    description => 'NDBC buoy-measured hourly atmospheric visibility.',
    begin_date => '1970-02-26',
    end_date   => '',
    factor => 1000.0, # km to meters.
    min => 0.0,
    max => 100000.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,vis'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind_speed',
    units => 'm/s',
    description => 'NDBC buoy-measured hourly average wind speed.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wspd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind_gust_speed',
    units => 'm/s',
    description => 'NDBC buoy-measured hourly peak wind gust speed.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,gst'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind_direction',
    units => 'bearing',
    description =>
       'NDBC buoy-measured hourly average wind bearing direction ' .
       '(bearing angle).',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 360.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind_u',
    units => 'm/s',
    description =>
       'NDBC buoy-measured hourly average wind speed eastward component.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => -100.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wspu'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind_v',
    units => 'm/s',
    description =>
       'NDBC buoy-measured hourly average wind speed northward component.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => -100.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wspv'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind',
    units => 'm/s',
    description =>
       'NDBC buoy-measured hourly average wind vector u and v components.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => -100.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wspu,wspv'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind_speed_and_wind_direction',
    units => 'm/s',
    description =>
      'NDBC buoy-measured hourly average wind speed and direction.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wspd,wd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wind_gust_speed_and_wind_direction',
    units => 'm/s',
    description =>
      'NDBC buoy-measured hourly average wind speed and direction.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,gst,wd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.water_temperature',
    units => 'C',
    description => 'NDBC buoy-measured hourly water temperature.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wtmp'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.tide',
    units => 'm',
    description =>
      'NDBC buoy-measured hourly tide water level above or below ' .
      'Mean Lower Low Water (MLLW).',
    begin_date => '1970-02-26',
    end_date   => '',
    min => -30.0,
    max => 30.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,tide'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wave_height',
    units => 'm',
    description =>
      'NDBC buoy-measured hourly wave height average of the highest '.
      'one-third of all of the wave heights during the 20-minute sampling ' .
      'period.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,wvht'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wave_average_period',
    units => 's',
    description =>
      'NDBC buoy-measured hourly average wave period of all waves during ' .
      'the 20-minute period.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 600.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,apd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wave_dominant_period',
    units => 's',
    description =>
      'NDBC buoy-measured hourly dominant wave period with the maximum ' .
      'wave energy.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 600.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,dpd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.wave_direction',
    units => 'degrees',
    description => 'NDBC buoy-measured hourly dominant wave direction.',
    begin_date => '1970-02-26',
    end_date   => '',
    min => 0.0,
    max => 360.0,
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,mwd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NDBC.all', # Get all measured data.
    units => '',
    description => 'NDBC buoy-measured hourly met data.',
    begin_date => '1970-02-26',
    end_date   => '',
    minimums => [
      0.0, 0.0, 0.0, 0.0, # Unused: timestamp, station, longitude, latitude.
      -50.0,
      -50.0,
      500.0,
      -500.0,
      0.0,
      0.0,
      0.0,
      0.0,
      -100.0,
      -100.0,
      0.0,
      -30.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    maximums => [
      0.0, 0.0, 0.0, 0.0, # Unused: timestamp, station, longitude, latitude.
      50.0,
      50.0,
      1500.0,
      500.0,
      100000.0,
      360.0,
      100.0,
      100.0,
      100.0,
      100.0,
      50.0,
      30.0,
      100.0,
      600.0,
      600.0,
      360.0
    ],
    factors => [
      0.0, 0.0, 0.0, 0.0, # Unused: timestamp, station, longitude, latitude.
      1.0,
      1.0,
      1.0,
      1.0,
      1000.0, # kilometers to meters.
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0
    ],
    output_header =>
      "\tair_temperature(C)" .
      "\tdew_point(C)" .
      "\tpressure(hPa)" .
      "\tpressure_tendency(hPa)" .
      "\tair_visibility(m)" .
      "\twind_direction(bearing)" .
      "\twind_speed(m/s)" .
      "\twind_gust_speed(m/s)" .
      "\twind_u(m/s)" .
      "\twind_v(m/s)" .
      "\twater_temperature(C)" .
      "\ttide(m)" .
      "\twave_height(m)" .
      "\twave_average_period(s)" .
      "\twave_dominant_period(s)" .
      "\twave_direction(degrees)\n",
    query => 'https://coastwatch.pfeg.noaa.gov/erddap/' .
               'tabledap/cwwcNDBCMet.csv?' .
               'time,station,longitude,latitude,' .
               'atmp,dewp,bar,ptdy,vis,' .
               'wd,wspd,gst,wspu,wspv,' .
               'wtmp,tide,wvht,apd,dpd,mwd'.
               '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
               '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
               '&time%3E=START_TIME&time%3C=END_TIME' .
               '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 20,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },

  # FIX: as of 2024-05-23 neracoosSos is not found at
  # https://www.neracoos.org/erddap/tabledap/index.html?page=1&itemsPerPage=1000
  # HACK for now use the air temperature webservice.
  # NERACOOS SOS buoy data looks like this:
  #time,station_id,longitude,latitude,altitude,air_temperature,sea_level_pressure,visibility_in_air,wind_from_direction,wind_speed,wind_gust,sea_water_temperature,sea_water_density,sea_water_salinity,sea_water_electrical_conductivity,direction_of_sea_water_velocity,sea_water_speed,wave_height,dominant_wave_period,chlorophyll
  #UTC,,degrees_east,degrees_north,m,degree_C,mbar,m,degrees_true,m s-1,m s-1,degree_C,kg m-3,PSU,S m-1,degrees_true,cm s-1,m,s,mg m-3
  #UTC,,degrees_east,degrees_north,degree_C,mbar,m,degrees_true,m s-1,m s-1,degree_C,kg m-3,PSU,S m-1,degrees_true,cm s-1,m,s,mg m-3
  #2010-08-17T00:00:00Z,44022,-73.73,40.88,24.7,1013.5,NaN,240.0,1.0,2.0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN
  #2010-08-17T00:00:00Z,44022,-73.73,40.88,NaN,NaN,NaN,NaN,NaN,NaN,23.3,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN
  #2010-08-17T00:00:00Z,44039,-72.655,41.138,24.7,1014.1,NaN,210.0,4.0,5.0,NaN,NaN,NaN,NaN,NaN,NaN,0.4,2.0,NaN
  #...
  #2010-08-19T23:50:00Z,I01,-68.1087,44.1058,NaN,NaN,NaN,162.199996948242,2.28099989891052,2.65799999237061,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN
  #2010-08-19T23:50:00Z,M01,-67.8716659545898,43.490140914917,19.1599998474121,1009.65179443359,2956.12133789062,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN
  #2010-08-19T23:50:00Z,M01,-67.8716659545898,43.490140914917,NaN,NaN,NaN,171.5,2.18099999427795,2.43400001525879,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN
  {
    name => 'NERACOOSSOS.station',
    units => '-',
    description => 'NERACOOS-SOS buoy station id.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => '',
    max => '',
    query =>
      #'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/neracoosSos.csv?' .
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosATemp.csv?' .
      'station_id,longitude,latitude' .
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 3,
    timestamp    => -1,
    longitude    => 1,
    latitude     => 2,
    depth        => -1,
    id           => -1,
    sid          => 0,
    value        => -1
  },
  {
    name => 'NERACOOSSOS.air_temperature',
    units => 'C',
    description =>
      'NERACOOS-SOS buoy station-measured hourly air temperature.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => -50.0,
    max => 50.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosATemp.csv?' .
      'time,station_id,longitude,latitude,altitude,air_temperature' .
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => -1,
    sid          => 1,
    value        => 5
  },
  {
    name => 'NERACOOSSOS.air_pressure',
    units => 'hPa',
    description =>
      'NERACOOS-SOS buoy station-measured hourly atmospheric pressure.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 500.0,
    max => 1500.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosBPres.csv?' .
      'time,station_id,longitude,latitude,altitude,air_pressure' .
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => -1,
    sid          => 1,
    value        => 5
  },
  {
    name => 'NERACOOSSOS.precipitation',
    units => 'mm',
    description =>
      'NERACOOS-SOS buoy station-measured 6-minute precipitation.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsMRF.csv?' .
      'time,stationID,longitude,latitude,RF' .
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.humidity',
    units => '%',
    description =>
      'NERACOOS-SOS buoy station-measured 6-minute relative humidity.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsMRH.csv?' .
      'time,stationID,longitude,latitude,RH' .
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.air_visibility',
    units => 'm',
    description =>
      'NERACOOS-SOS buoy station-measured hourly visibility in the air.',
    begin_date => '1970-01-01',
    end_date   => '',
    factor => 1852.0, # Nautical miles to meters.
    min => 0.0,
    max => 100000.0,
    query =>
      #'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/neracoosSos.csv?' .
      #'time,station_id,longitude,latitude,altitude,visibility_in_air'.
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsMV.csv?' .
      'time,stationID,longitude,latitude,Vis'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.wind_direction',
    units => 'bearing',
    description =>
      'NERACOOS-SOS buoy station-measured hourly wind from direction ' .
      '(bearing angle).',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 360.0,
    query =>
      #'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/neracoosSos.csv?' .
      #'time,station_id,longitude,latitude,altitude,wind_from_direction'.
      # FIX: 2024-05-28 The new webservice below returns no data.
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsMW.csv?' .
      'time,stationID,longitude,latitude,WD'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.wind_speed',
    units => 'm/s',
    description => 'NERACOOS-SOS buoy station-measured hourly wind speed.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query =>
       #'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/neracoosSos.csv?' .
       #'time,station_id,longitude,latitude,altitude,wind_speed'.
       # FIX: 2024-05-28 The new webservice below returns no data.
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsMW.csv?' .
      'time,stationID,longitude,latitude,WS'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.wind_gust_speed',
    units => 'm/s',
    description =>
      'NERACOOS-SOS buoy station-measured hourly wind gust speed.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query =>
      #'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/neracoosSos.csv?' .
      #'time,station_id,longitude,latitude,altitude,wind_gust'.
      # FIX: 2024-05-28 The new webservice below returns no data.
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsMW.csv?' .
      'time,stationID,longitude,latitude,WG'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.wind',
    units => 'm/s',
    description => 'NERACOOS-SOS buoy station-measured hourly wind vector.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query =>
      #'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/neracoosSos.csv?' .
      #'time,station_id,longitude,latitude,altitude,' .
      #'wind_speed,wind_from_direction'.
      # FIX: 2024-05-28 The new webservice below returns no data.
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsMW.csv?' .
      'time,stationID,longitude,latitude,WS,WD'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.current_direction',
    units => 'degrees',
    description =>
      'NERACOOS-SOS buoy station-measured hourly ' .
      'water current flow direction (degrees).',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 360.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsCA.csv?' .
      'time,stationID,longitude,latitude,CD'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.current_speed',
    units => 'm/s',
    description => 'NERACOOS-SOS buoy station-measured hourly current speed.',
    begin_date => '1970-01-01',
    end_date   => '',
    factor => 0.51444, # knots to meters/second
    min => 0.0,
    max => 20.0, # After multiplying by factor.
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsCA.csv?' .
      'time,stationID,longitude,latitude,CS'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 5,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.current',
    units => 'm/s',
    description =>
      'NERACOOS-SOS buoy station-measured hourly water current flow vector.',
    begin_date => '1970-01-01',
    end_date   => '',
    factor => 0.51444, # knots to meters/second
    min => 0.0,
    max => 20.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosCoopsCA.csv?' .
      'time,stationID,longitude,latitude,CS,CD'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => -1,
    id           => -1,
    sid          => 1,
    value        => 4
  },
  {
    name => 'NERACOOSSOS.water_temperature',
    units => 'C',
    description =>
      'NERACOOS-SOS buoy station-measured hourly water temperature.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query =>
      #'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/neracoosSos.csv?' .
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosWTemp.csv?' .
      'time,station_id,longitude,latitude,altitude,sea_water_temperature'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => -1,
    sid          => 1,
    value        => 5
  },
  {
    name => 'NERACOOSSOS.water_level',
    units => 'm',
    description =>
      'NERACOOS-SOS buoy station-measured hourly water surface height above refernece datum.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosWLevel.csv?' .
      'time,station_id,longitude,latitude,altitude,water_level'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => -1,
    sid          => 1,
    value        => 5
  },
  {
    name => 'NERACOOSSOS.salinity',
    units => 'PSU',
    description =>
      'NERACOOS-SOS buoy station-measured hourly water salinity.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 50.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosSalinity.csv?' .
      'time,station_id,longitude,latitude,altitude,sea_water_salinity'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => -1,
    sid          => 1,
    value        => 5
  },
  {
    name => 'NERACOOSSOS.conductivity',
    units => 'S/m',
    description =>
      'NERACOOS-SOS buoy station-measured hourly water electrical ' .
      'conductivity.',
    begin_date => '1970-01-01',
    end_date   => '',
    min => 0.0,
    max => 100.0,
    query =>
      'https://coastwatch.pfeg.noaa.gov/erddap/tabledap/nosSosCond.csv?' .
      'time,station_id,longitude,latitude,altitude,' .
      'conductivity'.
      '&longitude%3E=MIN_LON&longitude%3C=MAX_LON' .
      '&latitude%3E=MIN_LAT&latitude%3C=MAX_LAT' .
      '&time%3E=START_TIME&time%3C=END_TIME' .
      '&distinct()',
    negative_lon => 1,
    file_format  => 'csv',
    header_lines => 2,
    fields       => 6,
    timestamp    => 0,
    longitude    => 2,
    latitude     => 3,
    depth        => 4,
    id           => -1,
    sid          => 1,
    value        => 5
  }
);

###############################################################################

sub get_coverage_metadata_reference { return \@coverage_metadata; }

1; # Required last executed statement of module. Do not add code below here!

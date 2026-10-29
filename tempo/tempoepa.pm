#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: tempoepa.pm - TEMPO coverage metadata used by tempoepaserver.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use tempoepa; # Use my $coverage_metadata_reference = tempoepa::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-03-23 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package tempoepa;

##############################################################################

my @coverage_metadata  = (

  # AODALH_L2 -----------------------------------------------------------------:

  {
    name =>        'l2.aodalh.aod550',
    units =>       '-',
    description => 'Aerosol Optical Depth at 550 nanometer wavelength.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.aodalh.alh',
    units =>       'km',
    description => 'Aerosol layer height.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.aodalh.aermodel',
    units =>       '-',
    description => 'Aerosol model number: 0 = smoke model 1, 1 = smoke model 2, 2 = dust, 3 = generic.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.aodalh.lwmask',
    units =>       '-',
    description => 'Land-water mask: 0 = water, 1 = land, 2 = coastal or shallow water.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.aodalh.dqf',
    units =>       '-',
    description => 'Data quality flag: 0 = high quality, 1 = medium quality, 2 = low quality, 3 = no retrieval.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },



  # PM25_L3 -----------------------------------------------------------------:

  {
    name =>        'l3.pm25.pm25sat_ge',
    units =>       'ug/m3',
    description => 'Hourly PM2.5 estimates from ABI AOD, TEMPO AOD ALH over CONUS on GOES-East grid (east of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.pm25sat_gw',
    units =>       'ug/m3',
    description => 'Hourly PM2.5 estimates from ABI AOD, TEMPO AOD ALH over CONUS on GOES-West grid (west of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.aod_ge',
    units =>       '-',
    description => 'Hourly composite ABI AOD and TEMPO AOD over CONUS on GOES-East grid (east of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.aod_gw',
    units =>       '-',
    description => 'Hourly composite ABI AOD and TEMPO AOD over CONUS on GOES-West grid (west of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.alh_ge',
    units =>       'km',
    description => 'Hourly TEMPO ALH over CONUS on GOES-East grid (east of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.alh_gw',
    units =>       'km',
    description => 'Hourly TEMPO ALH over CONUS on GOES-West grid (west of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.slope_aod_ge',
    units =>       'ug/m3',
    description => 'Slopes that convert AOD ALH to PM2.5 on GOES-East grid (for AOD, east of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.slope_aod_gw',
    units =>       'ug/m3',
    description => 'Slopes that convert AOD ALH to PM2.5 on GOES-East grid on GOES-West grid (for AOD, west of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.slope_alh_ge',
    units =>       'ug/m3/km',
    description => 'Slopes that convert AOD ALH to PM2.5 on GOES-East grid (for ALH, east of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.slope_alh_gw',
    units =>       'ug/m3/km',
    description => 'Slopes that convert AOD ALH to PM2.5 on GOES-West grid (for ALH, west of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.intercept_ge',
    units =>       'ug/m3',
    description => 'Intercepts that convert AOD ALH to PM2.5 on GOES-East grid (east of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.intercept_gw',
    units =>       'ug/m3',
    description => 'Intercepts that convert AOD ALH to PM2.5 on GOES-West grid (for ALH, west of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.count_aod_ge',
    units =>       '-',
    description => 'Number of AOD retrievals counts within an hour on GOES-East grid (east of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l3.pm25.count_aod_gw',
    units =>       '-',
    description => 'Number of AOD retrievals counts within an hour on GOES-West grid (west of 106W).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },



  # ADP_L2 -----------------------------------------------------------------:

  {
    name =>        'l2.adp.smoke',
    units =>       '-',
    description => 'smoke detection flag: 1 yes, 0 no.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.dust',
    units =>       '-',
    description => 'dust detection flag: 1 yes, 0 no.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.cloud',
    units =>       '-',
    description => 'cloud detection flag: 1 yes, 0 no.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.nuc',
    units =>       '-',
    description => 'none/unknown/clear-sky detection flag: 1 yes, 0 no".',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.snowice',
    units =>       '-',
    description => 'snow/ice detection flag: 1 yes, 0 no.',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.saai',
    units =>       '-',
    description => 'scaled absorbing aerosol index [0, 30].',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.dsdi',
    units =>       '-',
    description => 'dust/smoke discrimination index [-50, 50].',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.deepblue_aai',
    units =>       '-',
    description => 'deep-blue (410/440 nm) absorbing aerosol index [-30, 30].',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.uv_aai',
    units =>       '-',
    description => 'uv (354/384 nm) absorbing aerosol index [-50, 50].',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.qc_flag',
    units =>       '-',
    description => 'quality flag for smoke, dust, nuc, and cloud detections (see doc table 4).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.pqi1',
    units =>       '-',
    description => 'product quality information flag 1 (see doc table 5).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.pqi2',
    units =>       '-',
    description => 'product quality information flag 2 (see doc table 6).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.pqi3',
    units =>       '-',
    description => 'product quality information flag 3 (see doc table 7).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.pqi4',
    units =>       '-',
    description => 'product quality information flag 4 (see doc table 8).',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.std_dev_410nm',
    units =>       '-',
    description => 'standard deviation of the reflectance at 410 nm in a 3 by 3 box [0, 10].',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.std_dev_865nm',
    units =>       '-',
    description => 'standard deviation of the reflectance at 865 nm in a 3 by 3 box [0, 10].',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },
  {
    name =>        'l2.adp.std_dev_2210nm',
    units =>       '-',
    description => 'standard deviation of the reflectance at 2210 nm in a 3 by 3 box [0, 10].',
    begin_date =>  '2023-08-02',
    end_date   =>  '2023-08-31'
  },



  # NO2_L2 -------------------------------------------------------------------:
  # (subset of variables)

  {
    name =>        'l2.no2.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.vertical_column_sum',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column sum of ' .
    'vertical_column_troposphere and vertical_column_stratosphere .',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.vertical_column_troposphere',
    units =>       'molecules/cm2',
    description => 'Troposphere nitrogen dioxide vertical column.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.vertical_column_stratosphere',
    units =>       'molecules/cm2',
    description => 'Stratosphere nitrogen dioxide vertical column.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide fitted slant column.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.amf_stratosphere',
    units =>       '-',
    description => 'Nitrogen dioxide stratospheric air mass factor.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.no2.amf_troposphere',
    units =>       '-',
    description => 'Nitrogen dioxide tropospheric air mass factor.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },



  # HCHO_L2 -------------------------------------------------------------------:
  # (subset of variables)

  {
    name =>        'l2.hcho.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.vertical_column_sum',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column sum of ' .
    'vertical_column_troposphere and vertical_column_stratosphere .',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.vertical_column_troposphere',
    units =>       'molecules/cm2',
    description => 'Troposphere formaldehyde vertical column.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.vertical_column_stratosphere',
    units =>       'molecules/cm2',
    description => 'Stratosphere formaldehyde vertical column.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.amf_stratosphere',
    units =>       '-',
    description => 'Formaldehyde stratospheric air mass factor.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  },
  {
    name =>        'l2.hcho.amf_troposphere',
    units =>       '-',
    description => 'Formaldehyde tropospheric air mass factor.',
    begin_date =>  '2023-08-02',
    end_date   =>  ''
  }
);

###############################################################################

sub get_coverage_metadata_reference { return \@coverage_metadata; }

1; # Required last executed statement of module. Do not add code below here!


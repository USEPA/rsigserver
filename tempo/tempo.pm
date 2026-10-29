#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: tempo.pm - TEMPO coverage metadata used by temposerver.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use tempo; # Use my $coverage_metadata_reference = tempo::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-03-23 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package tempo;

##############################################################################

my @coverage_metadata  = (

  # NO2_L2 -------------------------------------------------------------------:

  {
    name =>        'l2.no2.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.no2.solar_azimuth_angle',
    units =>       'deg',
    description => 'Solar azimuth angle at pixel center.'
  },
  {
    name =>        'l2.no2.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.no2.viewing_azimuth_angle',
    units =>       'deg',
    description => 'Viewing azimuth angle at pixel center.'
  },
  {
    name =>        'l2.no2.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.no2.vertical_column_sum',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column sum of ' .
    'vertical_column_troposphere and vertical_column_stratosphere .'
  },
  {
    name =>        'l2.no2.vertical_column_total',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column determined from ' .
    'fitted slant column and total AMF calculated from surface to top ' .
    'of atmosphere.'
  },
  {
    name =>        'l2.no2.vertical_column_total_uncertainty',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l2.no2.vertical_column_troposphere',
    units =>       'molecules/cm2',
    description => 'Troposphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l2.no2.vertical_column_troposphere_uncertainty',
    units =>       'molecules/cm2',
    description => 'Troposphere nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l2.no2.vertical_column_stratosphere',
    units =>       'molecules/cm2',
    description => 'Stratosphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l2.no2.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.'
  },
  {
    name =>        'l2.no2.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide fitted slant column.'
  },
  {
    name =>        'l2.no2.fitted_slant_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide fitted slant column uncertainty.'
  },
  {
    name =>        'l2.no2.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l2.no2.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l2.no2.ground_pixel_quality_flag',
    units =>       '-',
    description => 'shallow_ocean land shallow_inland_water shoreline ' .
    'intermittent_water deep_inland_water continental_shelf_ocean ' .
    'deep_ocean land_water_error sun_glint_possibility ' .
    'solar_eclipse_possibility water evergreen_needleleaf_forest ' .
    'evergreen_broadleaf_forest deciduous_needleleaf_forest ' .
    'deciduous_broadleaf_forest mixed_forest closed_shrublands ' .
    'open_shrublands woody_savannas savannas grasslands permanent_wetlands ' .
    'croplands urban_and_built_up cropland_natural_vegetation_mosaic ' .
    'snow_and_ice barren_or_sparsely_vegetated unclassified fill_value.'
  },
  {
    name =>        'l2.no2.tropopause_pressure',
    units =>       'hPa',
    description => 'Tropopause pressure.'
  },
  {
    name =>        'l2.no2.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l2.no2.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l2.no2.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l2.no2.amf_total',
    units =>       '-',
    description => 'Total nitrogen dioxide air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere.'
  },
  {
    name =>        'l2.no2.amf_total_clear_sky',
    units =>       '-',
    description => 'Total nitrogen dioxide air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere assuming 0 cloud fraction.'
  },
  {
    name =>        'l2.no2.amf_diagnostic_flag',
    units =>       '-',
    description => 'geometric_AMF glint snow_correction no_cloud_pressure ' .
    'adjusted_surface_pressure adjusted_cloud_pressure no_albedo ' .
    'no_cloud_fraction no_gas_profile no_scattering_weights AMF_disabled.'
  },
  {
    name =>        'l2.no2.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l2.no2.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l2.no2.amf_cloud_pressure',
    units =>       '-',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l2.no2.amf_troposphere',
    units =>       '-',
    description => 'Nitrogen dioxide tropospheric air mass factor.'
  },
  {
    name =>        'l2.no2.amf_stratosphere',
    units =>       '-',
    description => 'Nitrogen dioxide stratospheric air mass factor.'
  },
  {
    name =>        'l2.no2.amf_troposphere_clear_sky',
    units =>       '-',
    description => 'Nitrogen dioxide clear sky tropospheric air mass factor.'
  },
  {
    name =>        'l2.no2.amf_stratosphere_clear_sky',
    units =>       '-',
    description => 'Nitrogen dioxide clear sky stratospheric air mass factor.'
  },
  {
    name =>        'l2.no2.fit_rms_residual',
    units =>       '-',
    description => 'Radiance fit RMS residual.'
  },
  {
    name =>        'l2.no2.fit_convergence_flag',
    units =>       '-',
    description => 'failed maxiter_exceeded suspect good.'
  },
  {
    name =>        'l2.no2.wind_speed',
    units =>       'm/s',
    description => '2-meter wind speed determined from the Eastward and ' .
    'Northward components.'
  },



  # NO2_L2_NRT ---------------------------------------------------------------:

  {
    name =>        'l2.nrt.no2.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.nrt.no2.solar_azimuth_angle',
    units =>       'deg',
    description => 'Solar azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.no2.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.nrt.no2.viewing_azimuth_angle',
    units =>       'deg',
    description => 'Viewing azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.no2.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.no2.vertical_column_sum',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column sum of ' .
    'vertical_column_troposphere and vertical_column_stratosphere .'
  },
  {
    name =>        'l2.nrt.no2.vertical_column_total',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column determined from ' .
    'fitted slant column and total AMF calculated from surface to top ' .
    'of atmosphere.'
  },
  {
    name =>        'l2.nrt.no2.vertical_column_total_uncertainty',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l2.nrt.no2.vertical_column_troposphere',
    units =>       'molecules/cm2',
    description => 'Troposphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l2.nrt.no2.vertical_column_troposphere_uncertainty',
    units =>       'molecules/cm2',
    description => 'Troposphere nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l2.nrt.no2.vertical_column_stratosphere',
    units =>       'molecules/cm2',
    description => 'Stratosphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l2.nrt.no2.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.'
  },
  {
    name =>        'l2.nrt.no2.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide fitted slant column.'
  },
  {
    name =>        'l2.nrt.no2.fitted_slant_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide fitted slant column uncertainty.'
  },
  {
    name =>        'l2.nrt.no2.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l2.nrt.no2.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l2.nrt.no2.ground_pixel_quality_flag',
    units =>       '-',
    description => 'shallow_ocean land shallow_inland_water shoreline ' .
    'intermittent_water deep_inland_water continental_shelf_ocean ' .
    'deep_ocean land_water_error sun_glint_possibility ' .
    'solar_eclipse_possibility water evergreen_needleleaf_forest ' .
    'evergreen_broadleaf_forest deciduous_needleleaf_forest ' .
    'deciduous_broadleaf_forest mixed_forest closed_shrublands ' .
    'open_shrublands woody_savannas savannas grasslands permanent_wetlands ' .
    'croplands urban_and_built_up cropland_natural_vegetation_mosaic ' .
    'snow_and_ice barren_or_sparsely_vegetated unclassified fill_value.'
  },
  {
    name =>        'l2.nrt.no2.tropopause_pressure',
    units =>       'hPa',
    description => 'Tropopause pressure.'
  },
  {
    name =>        'l2.nrt.no2.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l2.nrt.no2.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l2.nrt.no2.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l2.nrt.no2.amf_total',
    units =>       '-',
    description => 'Total nitrogen dioxide air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere.'
  },
  {
    name =>        'l2.nrt.no2.amf_total_clear_sky',
    units =>       '-',
    description => 'Total nitrogen dioxide air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere assuming 0 cloud fraction.'
  },
  {
    name =>        'l2.nrt.no2.amf_diagnostic_flag',
    units =>       '-',
    description => 'geometric_AMF glint snow_correction no_cloud_pressure ' .
    'adjusted_surface_pressure adjusted_cloud_pressure no_albedo ' .
    'no_cloud_fraction no_gas_profile no_scattering_weights AMF_disabled.'
  },
  {
    name =>        'l2.nrt.no2.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l2.nrt.no2.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l2.nrt.no2.amf_cloud_pressure',
    units =>       '-',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l2.nrt.no2.amf_troposphere',
    units =>       '-',
    description => 'Nitrogen dioxide tropospheric air mass factor.'
  },
  {
    name =>        'l2.nrt.no2.amf_stratosphere',
    units =>       '-',
    description => 'Nitrogen dioxide stratospheric air mass factor.'
  },
  {
    name =>        'l2.nrt.no2.amf_troposphere_clear_sky',
    units =>       '-',
    description => 'Nitrogen dioxide clear sky tropospheric air mass factor.'
  },
  {
    name =>        'l2.nrt.no2.amf_stratosphere_clear_sky',
    units =>       '-',
    description => 'Nitrogen dioxide clear sky stratospheric air mass factor.'
  },
  {
    name =>        'l2.nrt.no2.fit_rms_residual',
    units =>       '-',
    description => 'Radiance fit RMS residual.'
  },
  {
    name =>        'l2.nrt.no2.fit_convergence_flag',
    units =>       '-',
    description => 'failed maxiter_exceeded suspect good.'
  },
  {
    name =>        'l2.nrt.no2.wind_speed',
    units =>       'm/s',
    description => '2-meter wind speed determined from the Eastward and ' .
    'Northward components.'
  },



  # HCHO_L2 ------------------------------------------------------------------:

  {
    name =>        'l2.hcho.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.hcho.solar_azimuth_angle',
    units =>       'deg',
    description => 'Solar azimuth angle at pixel center.'
  },
  {
    name =>        'l2.hcho.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.hcho.viewing_azimuth_angle',
    units =>       'deg',
    description => 'Viewing azimuth angle at pixel center.'
  },
  {
    name =>        'l2.hcho.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.hcho.vertical_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column determined from fitted ' .
    'slant column and total AMF calculated from surface to top of atmosphere.'
  },
  {
    name =>        'l2.hcho.vertical_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column uncertainty.'
  },
  {
    name =>        'l2.hcho.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.'
  },
  {
    name =>        'l2.hcho.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column.'
  },
  {
    name =>        'l2.hcho.fitted_slant_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column uncertainty.'
  },
  {
    name =>        'l2.hcho.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l2.hcho.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l2.hcho.ground_pixel_quality_flag',
    units =>       '-',
    description => 'shallow_ocean land shallow_inland_water shoreline ' .
    'intermittent_water deep_inland_water continental_shelf_ocean deep_ocean ' .
    'land_water_error sun_glint_possibility solar_eclipse_possibility water ' .
    'evergreen_needleleaf_forest evergreen_broadleaf_forest ' .
    'deciduous_needleleaf_forest deciduous_broadleaf_forest mixed_forest ' .
    'closed_shrublands open_shrublands woody_savannas savannas grasslands ' .
    'permanent_wetlands croplands urban_and_built_up ' .
    'cropland_natural_vegetation_mosaic snow_and_ice ' .
    'barren_or_sparsely_vegetated unclassified fill_value.'
  },
  {
    name =>        'l2.hcho.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l2.hcho.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l2.hcho.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l2.hcho.amf',
    units =>       '-',
    description => 'Total formaldehyde air mass factor (AMF) calculated from ' .
    'surface to top of atmosphere.'
  },
  {
    name =>        'l2.hcho.amf_diagnostic_flag',
    units =>       '-',
    description => 'Geometric_AMF glint snow_correction no_cloud_pressure ' .
    'adjusted_surface_pressure adjusted_cloud_pressure no_albedo ' .
    'no_cloud_fraction no_gas_profile no_scattering_weights AMF_disabled.'
  },
  {
    name =>        'l2.hcho.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l2.hcho.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l2.hcho.amf_cloud_pressure',
    units =>       '-',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l2.hcho.amf_clear_sky',
    units =>       '-',
    description => 'Total formaldehyde air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere assuming 0 cloud fraction.'
  },
  {
    name =>        'l2.hcho.fit_rms_residual',
    units =>       '-',
    description => 'Radiance fit RMS residual.'
  },
  {
    name =>        'l2.hcho.fit_convergence_flag',
    units =>       '-',
    description => 'failed maxiter_exceeded suspect good.'
  },
  {
    name =>        'l2.hcho.wind_speed',
    units =>       'm/s',
    description => '2-meter wind speed determined from the Eastward and ' .
    'Northward components.'
  },



  # HCHO_L2_NRT ------------------------------------------------------------------:

  {
    name =>        'l2.nrt.hcho.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.nrt.hcho.solar_azimuth_angle',
    units =>       'deg',
    description => 'Solar azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.hcho.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.nrt.hcho.viewing_azimuth_angle',
    units =>       'deg',
    description => 'Viewing azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.hcho.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.hcho.vertical_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column determined from fitted ' .
    'slant column and total AMF calculated from surface to top of atmosphere.'
  },
  {
    name =>        'l2.nrt.hcho.vertical_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column uncertainty.'
  },
  {
    name =>        'l2.nrt.hcho.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.'
  },
  {
    name =>        'l2.nrt.hcho.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column.'
  },
  {
    name =>        'l2.nrt.hcho.fitted_slant_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column uncertainty.'
  },
  {
    name =>        'l2.nrt.hcho.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l2.nrt.hcho.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l2.nrt.hcho.ground_pixel_quality_flag',
    units =>       '-',
    description => 'shallow_ocean land shallow_inland_water shoreline ' .
    'intermittent_water deep_inland_water continental_shelf_ocean deep_ocean ' .
    'land_water_error sun_glint_possibility solar_eclipse_possibility water ' .
    'evergreen_needleleaf_forest evergreen_broadleaf_forest ' .
    'deciduous_needleleaf_forest deciduous_broadleaf_forest mixed_forest ' .
    'closed_shrublands open_shrublands woody_savannas savannas grasslands ' .
    'permanent_wetlands croplands urban_and_built_up ' .
    'cropland_natural_vegetation_mosaic snow_and_ice ' .
    'barren_or_sparsely_vegetated unclassified fill_value.'
  },
  {
    name =>        'l2.nrt.hcho.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l2.nrt.hcho.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l2.nrt.hcho.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l2.nrt.hcho.amf',
    units =>       '-',
    description => 'Total formaldehyde air mass factor (AMF) calculated from ' .
    'surface to top of atmosphere.'
  },
  {
    name =>        'l2.nrt.hcho.amf_diagnostic_flag',
    units =>       '-',
    description => 'Geometric_AMF glint snow_correction no_cloud_pressure ' .
    'adjusted_surface_pressure adjusted_cloud_pressure no_albedo ' .
    'no_cloud_fraction no_gas_profile no_scattering_weights AMF_disabled.'
  },
  {
    name =>        'l2.nrt.hcho.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l2.nrt.hcho.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l2.nrt.hcho.amf_cloud_pressure',
    units =>       '-',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l2.nrt.hcho.amf_clear_sky',
    units =>       '-',
    description => 'Total formaldehyde air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere assuming 0 cloud fraction.'
  },
  {
    name =>        'l2.nrt.hcho.fit_rms_residual',
    units =>       '-',
    description => 'Radiance fit RMS residual.'
  },
  {
    name =>        'l2.nrt.hcho.fit_convergence_flag',
    units =>       '-',
    description => 'failed maxiter_exceeded suspect good.'
  },
  {
    name =>        'l2.nrt.hcho.wind_speed',
    units =>       'm/s',
    description => '2-meter wind speed determined from the Eastward and ' .
    'Northward components.'
  },



  # O3TOT_L2 -----------------------------------------------------------------:

  {
    name =>        'l2.o3tot.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.o3tot.solar_azimuth_angle',
    units =>       'deg',
    description => 'Solar azimuth angle at pixel center.'
  },
  {
    name =>        'l2.o3tot.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.o3tot.viewing_azimuth_angle',
    units =>       'deg',
    description => 'Viewing azimuth angle at pixel center.'
  },
  {
    name =>        'l2.o3tot.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.o3tot.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l2.o3tot.column_amount_o3',
    units =>       'DU',
    description => 'Best total ozone solution.'
  },
  {
    name =>        'l2.o3tot.radiative_cloud_frac',
    units =>       '-',
    description => 'Cloud radiance fraction = fc*Ic331/Im331.'
  },
  {
    name =>        'l2.o3tot.fc',
    units =>       '-',
    description => 'Effective cloud fraction (mixed LER model).'
  },
  {
    name =>        'l2.o3tot.o3_below_cloud',
    units =>       '-',
    description => 'Ozone below fractional cloud".'
  },
  {
    name =>        'l2.o3tot.uv_aerosol_index',
    units =>       '-',
    description => 'UV aerosol index.'
  },
  {
    name =>        'l2.o3tot.cloud_pressure',
    units =>       'hPa',
    description => 'Effective cloud pressure.'
  },
  {
    name =>        'l2.o3tot.terrain_pressure',
    units =>       'hPa',
    description => 'Terrain pressure.'
  },
  {
    name =>        'l2.o3tot.surface_reflectivity_at_331nm',
    units =>       '-',
    description => 'Effective surface reflectivity at 331 nm.'
  },
  {
    name =>        'l2.o3tot.surface_reflectivity_at_360nm',
    units =>       '-',
    description => 'Effective surface reflectivity at 360 nm.'
  },
  {
    name =>        'l2.o3tot.step1_o3',
    units =>       'DU',
    description => 'Step 1 ozone solution.'
  },
  {
    name =>        'l2.o3tot.step2_o3',
    units =>       'DU',
    description => 'Step 2 ozone solution.'
  },

  # These variables are not yielding data from the initial files 2023-10-17:

  {
    name =>        'l2.o3tot.quality_flag',
    units =>       '-',
    description => 'Quality_flags.'
  },
  {
    name =>        'l2.o3tot.so2_index',
    units =>       '-',
    description => 'SO2 index.'
  },
  {
    name =>        'l2.o3tot.ground_pixel_quality_flag',
    units =>       '-',
    description => 'Ground pixel quality flag.'
  },
  {
    name =>        'l2.o3tot.lut_wavelength',
    units =>       'nm',
    description => 'Lookup table wavelength.'
  },
  {
    name =>        'l2.o3tot.algorithm_flags',
    units =>       '-',
    description => 'Algorithm_flags.'
  },
  {
    name =>        'l2.o3tot.radiance_bpix_flag_accepted',
    units =>       '-',
    description => 'Radiance bad pixel flag accepted.'
  },
  {
    name =>        'l2.o3tot.cal_adjustment',
    units =>       '-',
    description => 'Calibration adjustment.'
  },



  # CLDO4_L2 -----------------------------------------------------------------:

  {
    name =>        'l2.cldo4.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.cldo4.solar_azimuth_angle',
    units =>       'deg',
    description => 'Solar azimuth angle at pixel center.'
  },
  {
    name =>        'l2.cldo4.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.cldo4.viewing_azimuth_angle',
    units =>       'deg',
    description => 'Viewing azimuth angle at pixel center.'
  },
  {
    name =>        'l2.cldo4.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.cldo4.CloudRadianceFraction440',
    units =>       '-',
    description => 'Cloud radiance fraction at 440nm.'
  },
  {
    name =>        'l2.cldo4.CloudRadianceFraction466',
    units =>       '-',
    description => 'Cloud radiance fraction at 466nm.'
  },
  {
    name =>        'l2.cldo4.cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction at 466nm.'
  },
  {
    name =>        'l2.cldo4.cloud_pressure',
    units =>       'hPa',
    description => 'Optical centroid pressure for cloud.'
  },
  {
    name =>        'l2.cldo4.processing_quality_flag',
    units =>       '-',
    description => '0_error_geoloc_or_angles 1_warn_466nm_crf_invalid ' .
    '2_warn_ocp_repaced for_small_ecf 3_error_input_psurf_albedo ' .
    '4_warn_ocp_replaced_for_snowice 5_warn_ocp_scd_iteration_exceeds_max ' .
    '6_error_scd_negative_or_bad 7_info_440nm_rad_or_irr_error ' .
    '8_error_466nm_rad_or_irr_error 9_warn_ecf_beyond_normal_range ' .
    '10_info_pscene_SurfaceLER_TerrainP_invalid '.
    '11_info_pscene_SceneLER_SceneP_invalid '.
    '12_error_ecf_calc_skipped ' .
    '13_error_ocp_calc_skipped 14_warning_ocp_beyond_normal_range ' .
    '15_info_pscene_calc_skipped 16-29_unused_set_to_zero (4-byte only) ' .
    '30-31_reserved_set_to_zero (4-byte only).'
  },
  {
    name =>        'l2.cldo4.fit_convergence_flag',
    units =>       '-',
    description => 'failed maxiter_exceeded suspect good.'
  },
  {
    name =>        'l2.cldo4.fit_rms_residual',
    units =>       '-',
    description => 'Radiance fit RMS residual.'
  },
  {
    name =>        'l2.cldo4.GLER440',
    units =>       '-',
    description => '440nm surface reflectivity used in calculation.'
  },
  {
    name =>        'l2.cldo4.GLER466',
    units =>       '-',
    description => '466nm surface reflectivity used in calculation.'
  },
  {
    name =>        'l2.cldo4.SCD_MainDataQualityFlags',
    units =>       '-',
    description => 'Main data quality flags for fitted_slant_column: ' .
    '0=normal, 1=suspicious, 2=bad.'
  },
  {
    name =>        'l2.cldo4.SceneLER440',
    units =>       'hPa',
    description => 'LER at 440nm calculated at ScenePressure.'
  },
  {
    name =>        'l2.cldo4.SceneLER466',
    units =>       'hPa',
    description => 'LER at 466nm calculated at ScenePressure.'
  },
  {
    name =>        'l2.cldo4.ScenePressure',
    units =>       'hPa',
    description => 'Scene pressure.'
  },
  {
    name =>        'l2.cldo4.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l2.cldo4.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l2.cldo4.scddes',
    units =>       'molecules2/cm5',
    description => 'Destriped column amount.'
  },
  {
    name =>        'l2.cldo4.ground_pixel_quality_flag',
    units =>       '-',
    description => 'shallow_ocean land shallow_inland_water shoreline ' .
    'intermittent_water deep_inland_water continental_shelf_ocean ' .
    'deep_ocean land_water_error sun_glint_possibility ' .
    'solar_eclipse_possibility water evergreen_needleleaf_forest ' .
    'evergreen_broadleaf_forest deciduous_needleleaf_forest ' .
    'deciduous_broadleaf_forest mixed_forest closed_shrublands ' .
    'open_shrublands woody_savannas savannas grasslands permanent_wetlands ' .
    'croplands urban_and_built_up cropland_natural_vegetation_mosaic ' .
    'snow_and_ice barren_or_sparsely_vegetated unclassified fill_value.'
  },
  {
    name =>        'l2.cldo4.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l2.cldo4.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l2.cldo4.fitted_slant_column',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column.'
  },
  {
    name =>        'l2.cldo4.fitted_slant_column_uncertainty',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column ' .
    'uncertainty.'
  },



  # CLDO4_L2_NRT -------------------------------------------------------------:

  {
    name =>        'l2.nrt.cldo4.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.nrt.cldo4.solar_azimuth_angle',
    units =>       'deg',
    description => 'Solar azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.cldo4.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.nrt.cldo4.viewing_azimuth_angle',
    units =>       'deg',
    description => 'Viewing azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.cldo4.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.nrt.cldo4.CloudRadianceFraction440',
    units =>       '-',
    description => 'Cloud radiance fraction at 440nm.'
  },
  {
    name =>        'l2.nrt.cldo4.CloudRadianceFraction466',
    units =>       '-',
    description => 'Cloud radiance fraction at 466nm.'
  },
  {
    name =>        'l2.nrt.cldo4.cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction at 466nm.'
  },
  {
    name =>        'l2.nrt.cldo4.cloud_pressure',
    units =>       'hPa',
    description => 'Optical centroid pressure for cloud.'
  },
  {
    name =>        'l2.nrt.cldo4.processing_quality_flag',
    units =>       '-',
    description => '0_error_geoloc_or_angles 1_warn_466nm_crf_invalid ' .
    '2_warn_ocp_repaced for_small_ecf 3_error_input_psurf_albedo ' .
    '4_warn_ocp_replaced_for_snowice 5_warn_ocp_scd_iteration_exceeds_max ' .
    '6_error_scd_negative_or_bad 7_info_440nm_rad_or_irr_error ' .
    '8_error_466nm_rad_or_irr_error 9_warn_ecf_beyond_normal_range ' .
    '10_info_pscene_SurfaceLER_TerrainP_invalid '.
    '11_info_pscene_SceneLER_SceneP_invalid '.
    '12_error_ecf_calc_skipped ' .
    '13_error_ocp_calc_skipped 14_warning_ocp_beyond_normal_range ' .
    '15_info_pscene_calc_skipped 16-29_unused_set_to_zero (4-byte only) ' .
    '30-31_reserved_set_to_zero (4-byte only).'
  },
  {
    name =>        'l2.nrt.cldo4.fit_convergence_flag',
    units =>       '-',
    description => 'failed maxiter_exceeded suspect good.'
  },
  {
    name =>        'l2.nrt.cldo4.fit_rms_residual',
    units =>       '-',
    description => 'Radiance fit RMS residual.'
  },
  {
    name =>        'l2.nrt.cldo4.GLER440',
    units =>       '-',
    description => '440nm surface reflectivity used in calculation.'
  },
  {
    name =>        'l2.nrt.cldo4.GLER466',
    units =>       '-',
    description => '466nm surface reflectivity used in calculation.'
  },
  {
    name =>        'l2.nrt.cldo4.SCD_MainDataQualityFlags',
    units =>       '-',
    description => 'Main data quality flags for fitted_slant_column: ' .
    '0=normal, 1=suspicious, 2=bad.'
  },
  {
    name =>        'l2.nrt.cldo4.SceneLER440',
    units =>       'hPa',
    description => 'LER at 440nm calculated at ScenePressure.'
  },
  {
    name =>        'l2.nrt.cldo4.SceneLER466',
    units =>       'hPa',
    description => 'LER at 466nm calculated at ScenePressure.'
  },
  {
    name =>        'l2.nrt.cldo4.ScenePressure',
    units =>       'hPa',
    description => 'Scene pressure.'
  },
  {
    name =>        'l2.nrt.cldo4.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l2.nrt.cldo4.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l2.nrt.cldo4.scddes',
    units =>       'molecules2/cm5',
    description => 'Destriped column amount.'
  },
  {
    name =>        'l2.nrt.cldo4.ground_pixel_quality_flag',
    units =>       '-',
    description => 'shallow_ocean land shallow_inland_water shoreline ' .
    'intermittent_water deep_inland_water continental_shelf_ocean ' .
    'deep_ocean land_water_error sun_glint_possibility ' .
    'solar_eclipse_possibility water evergreen_needleleaf_forest ' .
    'evergreen_broadleaf_forest deciduous_needleleaf_forest ' .
    'deciduous_broadleaf_forest mixed_forest closed_shrublands ' .
    'open_shrublands woody_savannas savannas grasslands permanent_wetlands ' .
    'croplands urban_and_built_up cropland_natural_vegetation_mosaic ' .
    'snow_and_ice barren_or_sparsely_vegetated unclassified fill_value.'
  },
  {
    name =>        'l2.nrt.cldo4.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l2.nrt.cldo4.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l2.nrt.cldo4.fitted_slant_column',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column.'
  },
  {
    name =>        'l2.nrt.cldo4.fitted_slant_column_uncertainty',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column ' .
    'uncertainty.'
  },



  # O3PROF_L2 -------------------------------------------------------------------:

  {
    name =>        'l2.o3prof.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l2.o3prof.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l2.o3prof.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l2.o3prof.ozone_profile',
    units =>       'DU',
    description => 'Ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.ozone_profile_precision',
    units =>       'DU',
    description => 'Precision of ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.ozone_profile_error',
    units =>       'DU',
    description => 'Error of ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.total_ozone_column',
    units =>       'DU',
    description => 'Integrated total ozone from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.total_ozone_precision',
    units =>       'DU',
    description => 'Precision of integrated total ozone from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.total_ozone_error',
    units =>       'DU',
    description => 'Error of integrated total ozone from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.stratosphere_ozone_column',
    units =>       'DU',
    description => 'Integrated ozone in the stratosphere.'
  },
  {
    name =>        'l2.o3prof.stratosphere_ozone_precision',
    units =>       'DU',
    description => 'Precision of integrated ozone in the stratosphere.'
  },
  {
    name =>        'l2.o3prof.stratosphere_ozone_error',
    units =>       'DU',
    description => 'Error of integrated ozone in the stratosphere.'
  },
  {
    name =>        'l2.o3prof.troposphere_ozone_column',
    units =>       'DU',
    description => 'Integrated ozone in the troposphere.'
  },
  {
    name =>        'l2.o3prof.troposphere_ozone_precision',
    units =>       'DU',
    description => 'Precision of integrated ozone in the troposphere.'
  },
  {
    name =>        'l2.o3prof.troposphere_ozone_error',
    units =>       'DU',
    description => 'Error of integrated ozone in the troposphere.'
  },
  {
    name =>        'l2.o3prof.ground_pixel_quality_flag',
    units =>       '-',
    description => 'Ground pixel quality flag.'
  },
  {
    name =>        'l2.o3prof.aerosol_index',
    units =>       '-',
    description => 'Aerosol index.'
  },
  {
    name =>        'l2.o3prof.cloud_flag',
    units =>       '-',
    description => 'Cloud flag.'
  },
  {
    name =>        'l2.o3prof.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction.'
  },
  {
    name =>        'l2.o3prof.eff_cloud_pressure',
    units =>       'hPa',
    description => 'Effective cloud pressure.'
  },
  {
    name =>        'l2.o3prof.glint_probability',
    units =>       '-',
    description => 'Probability of sun glint.'
  },
  {
    name =>        'l2.o3prof.ozone_apriori_profile',
    units =>       'DU',
    description => 'A priori ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.ozone_apriori_profile_error',
    units =>       'DU',
    description => 'Error of a priori ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.ozone_profile_altitude',
    units =>       'km',
    description => 'Altitude of each of the 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.ozone_profile_temperature',
    units =>       'K',
    description => 'Temperature at each of the 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.ozone_profile_pressure',
    units =>       'hPa',
    description => 'Pressure at each of the 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l2.o3prof.tropopause_index',
    units =>       '-',
    description => 'Tropopause index.'
  },
  {
    name =>        'l2.o3prof.surface_albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l2.o3prof.num_wavelengths_in_fit',
    units =>       '-',
    description => 'Number of wavelengths used in fitting.'
  },
  {
    name =>        'l2.o3prof.ozone_information_content',
    units =>       '-',
    description => 'Ozone information content.'
  },
  {
    name =>        'l2.o3prof.ozone_degree_of_freedom',
    units =>       '-',
    description => 'Ozone degree of freedom for signal.'
  },
  {
    name =>        'l2.o3prof.exit_status',
    units =>       '-',
    description => 'retrieval exit status (0=not converged, 1-99 converged, 100-110 converged but with negative ozone value).'
  },
  {
    name =>        'l2.o3prof.num_iterations',
    units =>       '-',
    description => 'Number of iterations.'
  },



  # O3PROF_L3 -------------------------------------------------------------------:

  {
    name =>        'l3.o3prof.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.o3prof.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.o3prof.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.o3prof.ozone_profile',
    units =>       'DU',
    description => 'Ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.ozone_profile_precision',
    units =>       'DU',
    description => 'Precision of ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.ozone_profile_error',
    units =>       'DU',
    description => 'Error of ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.total_ozone_column',
    units =>       'DU',
    description => 'Integrated total ozone from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.total_ozone_precision',
    units =>       'DU',
    description => 'Precision of integrated total ozone from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.total_ozone_error',
    units =>       'DU',
    description => 'Error of integrated total ozone from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.stratosphere_ozone_column',
    units =>       'DU',
    description => 'Integrated ozone in the stratosphere.'
  },
  {
    name =>        'l3.o3prof.stratosphere_ozone_precision',
    units =>       'DU',
    description => 'Precision of integrated ozone in the stratosphere.'
  },
  {
    name =>        'l3.o3prof.stratosphere_ozone_error',
    units =>       'DU',
    description => 'Error of integrated ozone in the stratosphere.'
  },
  {
    name =>        'l3.o3prof.troposphere_ozone_column',
    units =>       'DU',
    description => 'Integrated ozone in the troposphere.'
  },
  {
    name =>        'l3.o3prof.troposphere_ozone_precision',
    units =>       'DU',
    description => 'Precision of integrated ozone in the troposphere.'
  },
  {
    name =>        'l3.o3prof.troposphere_ozone_error',
    units =>       'DU',
    description => 'Error of integrated ozone in the troposphere.'
  },
  {
    name =>        'l3.o3prof.weight',
    units =>       '/km2',
    description => 'Sum of Level 2 pixel overlap areas.'
  },
  {
    name =>        'l3.o3prof.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction.'
  },
  {
    name =>        'l3.o3prof.eff_cloud_pressure',
    units =>       'hPa',
    description => 'Effective cloud pressure.'
  },
  {
    name =>        'l3.o3prof.ozone_apriori_profile',
    units =>       'DU',
    description => 'A priori ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.ozone_apriori_profile_error',
    units =>       'DU',
    description => 'Error of a priori ozone at 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.ozone_profile_altitude',
    units =>       'km',
    description => 'Altitude of each of the 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.ozone_profile_pressure',
    units =>       'hPa',
    description => 'Pressure at each of the 24 layers from the surface up to 100km.'
  },
  {
    name =>        'l3.o3prof.tropopause_index',
    units =>       '-',
    description => 'Tropopause index.'
  },



  # NO2_L3 -------------------------------------------------------------------:

  {
    name =>        'l3.no2.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.no2.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.no2.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.no2.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l3.no2.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l3.no2.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l3.no2.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l3.no2.vertical_column_sum',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column sum of ' .
    'vertical_column_troposphere and vertical_column_stratosphere .'
  },
  {
    name =>        'l3.no2.vertical_column_total',
    units =>       'molecules2/cm2',
    description => 'Total nitrogen dioxide vertical column.'
  },
  {
    name =>        'l3.no2.vertical_column_total_uncertainty',
    units =>       'molecules2/cm2',
    description => 'Total nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l3.no2.vertical_column_troposphere',
    units =>       'molecules2/cm2',
    description => 'Troposphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l3.no2.vertical_column_troposphere_uncertainty',
    units =>       'molecules2/cm2',
    description => 'Troposphere nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l3.no2.vertical_column_stratosphere',
    units =>       'molecules2/cm2',
    description => 'Stratosphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l3.no2.fitted_slant_column',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column.'
  },
  {
    name =>        'l3.no2.fitted_slant_column_uncertainty',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column ' .
    'uncertainty.'
  },
  {
    name =>        'l3.no2.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l3.no2.tropopause_pressure',
    units =>       'hPa',
    description => 'Tropopause pressure.'
  },
  {
    name =>        'l3.no2.amf_total',
    units =>       '-',
    description => 'Total nitrogen dioxide air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere.'
  },
  {
    name =>        'l3.no2.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l3.no2.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l3.no2.amf_cloud_pressure',
    units =>       'hPa',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l3.no2.amf_troposphere',
    units =>       '-',
    description => 'Nitrogen dioxide tropospheric air mass factor.'
  },
  {
    name =>        'l3.no2.amf_stratosphere',
    units =>       '-',
    description => 'Nitrogen dioxide stratospheric air mass factor.'
  },
  {
    name =>        'l3.no2.weight',
    units =>       'km2',
    description => 'Sum of area weights.'
  },
  {
    name =>        'l3.no2.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.'
  },
  {
    name =>        'l3.no2.num_vertical_column_total_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.min_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.max_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.num_vertical_column_troposphere_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.min_vertical_column_troposphere_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.max_vertical_column_troposphere_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.num_vertical_column_troposphere_uncertainty_samples',
    units =>       '-',
    description => 'Uncertainty of number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.min_vertical_column_troposphere_uncertainty_sample',
    units =>       'molecules/cm2',
    description => 'Uncertainty of smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.max_vertical_column_troposphere_uncertainty_sample',
    units =>       'molecules/cm2',
    description => 'Uncertainty of largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.num_vertical_column_stratosphere_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.min_vertical_column_stratosphere_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.no2.max_vertical_column_stratosphere_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },



  # NO2_L3_NRT ---------------------------------------------------------------:

  {
    name =>        'l3.nrt.no2.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.nrt.no2.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.nrt.no2.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.nrt.no2.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l3.nrt.no2.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l3.nrt.no2.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l3.nrt.no2.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l3.nrt.no2.vertical_column_sum',
    units =>       'molecules/cm2',
    description => 'Nitrogen dioxide vertical column sum of ' .
    'vertical_column_troposphere and vertical_column_stratosphere .'
  },
  {
    name =>        'l3.nrt.no2.vertical_column_total',
    units =>       'molecules2/cm2',
    description => 'Total nitrogen dioxide vertical column.'
  },
  {
    name =>        'l3.nrt.no2.vertical_column_total_uncertainty',
    units =>       'molecules2/cm2',
    description => 'Total nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l3.nrt.no2.vertical_column_troposphere',
    units =>       'molecules2/cm2',
    description => 'Troposphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l3.nrt.no2.vertical_column_troposphere_uncertainty',
    units =>       'molecules2/cm2',
    description => 'Troposphere nitrogen dioxide vertical column uncertainty.'
  },
  {
    name =>        'l3.nrt.no2.vertical_column_stratosphere',
    units =>       'molecules2/cm2',
    description => 'Stratosphere nitrogen dioxide vertical column.'
  },
  {
    name =>        'l3.nrt.no2.fitted_slant_column',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column.'
  },
  {
    name =>        'l3.nrt.no2.fitted_slant_column_uncertainty',
    units =>       'molecules2/cm5',
    description => 'Collision induced oxygen complex fitted slant column ' .
    'uncertainty.'
  },
  {
    name =>        'l3.nrt.no2.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l3.nrt.no2.tropopause_pressure',
    units =>       'hPa',
    description => 'Tropopause pressure.'
  },
  {
    name =>        'l3.nrt.no2.amf_total',
    units =>       '-',
    description => 'Total nitrogen dioxide air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere.'
  },
  {
    name =>        'l3.nrt.no2.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l3.nrt.no2.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l3.nrt.no2.amf_cloud_pressure',
    units =>       'hPa',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l3.nrt.no2.amf_troposphere',
    units =>       '-',
    description => 'Nitrogen dioxide tropospheric air mass factor.'
  },
  {
    name =>        'l3.nrt.no2.amf_stratosphere',
    units =>       '-',
    description => 'Nitrogen dioxide stratospheric air mass factor.'
  },
  {
    name =>        'l3.nrt.no2.weight',
    units =>       'km2',
    description => 'Sum of area weights.'
  },
  {
    name =>        'l3.nrt.no2.main_data_quality_flag',
    units =>       '-',
    description => 'Main data quality flag: 0 = normal, 1 = suspect, 2 = bad.'
  },
  {
    name =>        'l3.nrt.no2.num_vertical_column_total_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.min_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.max_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.num_vertical_column_troposphere_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.min_vertical_column_troposphere_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.max_vertical_column_troposphere_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.num_vertical_column_troposphere_uncertainty_samples',
    units =>       '-',
    description => 'Uncertainty of number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.min_vertical_column_troposphere_uncertainty_sample',
    units =>       'molecules/cm2',
    description => 'Uncertainty of smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.max_vertical_column_troposphere_uncertainty_sample',
    units =>       'molecules/cm2',
    description => 'Uncertainty of largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.num_vertical_column_stratosphere_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.min_vertical_column_stratosphere_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.no2.max_vertical_column_stratosphere_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the area-weighted Level 3 pixel value.'
  },



  # HCHO_L3 ------------------------------------------------------------------:

  {
    name =>        'l3.hcho.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.hcho.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.hcho.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.hcho.vertical_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column determined from fitted ' .
    'slant column and total AMF calculated from surface to top of atmosphere.'
  },
  {
    name =>        'l3.hcho.vertical_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column determined from fitted ' .
    'slant column and total AMF calculated from surface to top of atmosphere ' .
    'uncertainty.'
  },
  {
    name =>        'l3.hcho.main_data_quality_flag',
    units =>       '-',
    description => 'normal suspicious bad.'
  },
  {
    name =>        'l3.hcho.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l3.hcho.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l3.hcho.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l3.hcho.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l3.hcho.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column.'
  },
  {
    name =>        'l3.hcho.fitted_slant_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column uncertainty.'
  },
  {
    name =>        'l3.hcho.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l3.hcho.amf',
    units =>       '-',
    description => 'Total formaldehyde air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere.'
  },
  {
    name =>        'l3.hcho.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l3.hcho.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l3.hcho.amf_cloud_pressure',
    units =>       'hPa',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l3.hcho.weights',
    units =>       'km2',
    description => 'Sum of Level 2 pixel overlap areas.'
  },
  {
    name =>        'l3.hcho.num_vertical_column_total_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.hcho.min_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.hcho.max_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },



  # HCHO_L3_NRT --------------------------------------------------------------:

  {
    name =>        'l3.nrt.hcho.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.nrt.hcho.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.nrt.hcho.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.nrt.hcho.vertical_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column determined from fitted ' .
    'slant column and total AMF calculated from surface to top of atmosphere.'
  },
  {
    name =>        'l3.nrt.hcho.vertical_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde vertical column determined from fitted ' .
    'slant column and total AMF calculated from surface to top of atmosphere ' .
    'uncertainty.'
  },
  {
    name =>        'l3.nrt.hcho.main_data_quality_flag',
    units =>       '-',
    description => 'normal suspicious bad.'
  },
  {
    name =>        'l3.nrt.hcho.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l3.nrt.hcho.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l3.nrt.hcho.pbl_height',
    units =>       'm',
    description => 'Planetary boundary layer height.'
  },
  {
    name =>        'l3.nrt.hcho.snow_ice_fraction',
    units =>       '-',
    description => 'Fraction of pixel area covered by snow and/or ice.'
  },
  {
    name =>        'l3.nrt.hcho.fitted_slant_column',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column.'
  },
  {
    name =>        'l3.nrt.hcho.fitted_slant_column_uncertainty',
    units =>       'molecules/cm2',
    description => 'Formaldehyde fitted slant column uncertainty.'
  },
  {
    name =>        'l3.nrt.hcho.albedo',
    units =>       '-',
    description => 'Surface albedo.'
  },
  {
    name =>        'l3.nrt.hcho.amf',
    units =>       '-',
    description => 'Total formaldehyde air mass factor (AMF) calculated ' .
    'from surface to top of atmosphere.'
  },
  {
    name =>        'l3.nrt.hcho.eff_cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction from cloud retrieval.'
  },
  {
    name =>        'l3.nrt.hcho.amf_cloud_fraction',
    units =>       '-',
    description => 'Cloud radiance fraction for AMF computation.'
  },
  {
    name =>        'l3.nrt.hcho.amf_cloud_pressure',
    units =>       'hPa',
    description => 'Cloud pressure for AMF computation.'
  },
  {
    name =>        'l3.nrt.hcho.weights',
    units =>       'km2',
    description => 'Sum of Level 2 pixel overlap areas.'
  },
  {
    name =>        'l3.nrt.hcho.num_vertical_column_total_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.hcho.min_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Smallest Level 2 pixel value contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.nrt.hcho.max_vertical_column_total_sample',
    units =>       'molecules/cm2',
    description => 'Largest Level 2 pixel value contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },



  # O3TOT_L3 -----------------------------------------------------------------:

  {
    name =>        'l3.o3tot.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.o3tot.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.o3tot.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.o3tot.terrain_height',
    units =>       'm',
    description => 'Terrain height.'
  },
  {
    name =>        'l3.o3tot.column_amount_o3',
    units =>       'DU',
    description => 'Best total ozone solution.'
  },
  {
    name =>        'l3.o3tot.radiative_cloud_frac',
    units =>       '-',
    description => 'Cloud radiance fraction = fc*Ic331/Im331.'
  },
  {
    name =>        'l3.o3tot.fc',
    units =>       '-',
    description => 'Effective cloud fraction (mixed LER model).'
  },
  {
    name =>        'l3.o3tot.o3_below_cloud',
    units =>       'DU',
    description => 'Ozone below fractional cloud.'
  },
  {
    name =>        'l3.o3tot.so2l_index',
    units =>       '-',
    description => 'SO2 index.'
  },
  {
    name =>        'l3.o3tot.uv_aerosol_index',
    units =>       '-',
    description => 'UV aerosol index.'
  },
  {
    name =>        'l3.o3tot.so2_index',
    units =>       '-',
    description => 'SO2 index.'
  },
  {
    name =>        'l3.o3tot.uv_aerosol_index',
    units =>       '-',
    description => 'UV aerosol index.'
  },
  {
    name =>        'l3.o3tot.num_column_samples',
    units =>       '-',
    description => 'Number of Level 2 pixel values contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.o3tot.min_column_samples',
    units =>       'DU',
    description => 'Smallest Level 2 pixel value contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.o3tot.max_column_samples',
    units =>       'DU',
    description => 'Largest Level 2 pixel value contributing to the ' .
    'area-weighted Level 3 pixel value.'
  },
  {
    name =>        'l3.o3tot.cloud_pressures',
    units =>       'hPa',
    description => 'Effective cloud pressure.'
  },
  {
    name =>        'l3.o3tot.terrain_pressures',
    units =>       'hPa',
    description => 'Terrain pressure.'
  },
  {
    name =>        'l3.o3tot.weight',
    units =>       'km2',
    description => 'Sum of Level 2 pixel overlap areas.'
  },



  # CLDO4_L3 -----------------------------------------------------------------:

  {
    name =>        'l3.cldo4.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.cldo4.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.cldo4.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.cldo4.cloud_pressure',
    units =>       'hPa',
    description => 'Optical centroid pressure for cloud.'
  },
  {
    name =>        'l3.cldo4.cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction at 466nm.'
  },
  {
    name =>        'l3.cldo4.CloudRadianceFraction440',
    units =>       '-',
    description => 'Cloud radiance fraction at 440nm.'
  },
  {
    name =>        'l3.cldo4.CloudRadianceFraction466',
    units =>       '-',
    description => 'Cloud radiance fraction at 466nm.'
  },
  {
    name =>        'l3.cldo4.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l3.cldo4.weight',
    units =>       'km2',
    description => 'Sum of Level 2 pixel overlap areas.'
  },
  {
    name =>        'l3.cldo4.GLER440',
    units =>       '-',
    description => '440nm surface reflectivity used in calculation.'
  },
  {
    name =>        'l3.cldo4.GLER466',
    units =>       '-',
    description => '466nm surface reflectivity used in calculation.'
  },


  # CLDO4_L3_NRT -------------------------------------------------------------:

  {
    name =>        'l3.nrt.cldo4.solar_zenith_angle',
    units =>       'deg',
    description => 'Solar zenith angle at pixel center.'
  },
  {
    name =>        'l3.nrt.cldo4.viewing_zenith_angle',
    units =>       'deg',
    description => 'Viewing zenith angle at pixel center.'
  },
  {
    name =>        'l3.nrt.cldo4.relative_azimuth_angle',
    units =>       'deg',
    description => 'Relative azimuth angle at pixel center.'
  },
  {
    name =>        'l3.nrt.cldo4.cloud_pressure',
    units =>       'hPa',
    description => 'Optical centroid pressure for cloud.'
  },
  {
    name =>        'l3.nrt.cldo4.cloud_fraction',
    units =>       '-',
    description => 'Effective cloud fraction at 466nm.'
  },
  {
    name =>        'l3.nrt.cldo4.CloudRadianceFraction440',
    units =>       '-',
    description => 'Cloud radiance fraction at 440nm.'
  },
  {
    name =>        'l3.nrt.cldo4.CloudRadianceFraction466',
    units =>       '-',
    description => 'Cloud radiance fraction at 466nm.'
  },
  {
    name =>        'l3.nrt.cldo4.surface_pressure',
    units =>       'hPa',
    description => 'Surface pressure.'
  },
  {
    name =>        'l3.nrt.cldo4.weight',
    units =>       'km2',
    description => 'Sum of Level 2 pixel overlap areas.'
  },
  {
    name =>        'l3.nrt.cldo4.GLER440',
    units =>       '-',
    description => '440nm surface reflectivity used in calculation.'
  },
  {
    name =>        'l3.nrt.cldo4.GLER466',
    units =>       '-',
    description => '466nm surface reflectivity used in calculation.'
  }

);

###############################################################################

sub get_coverage_metadata_reference { return \@coverage_metadata; }

1; # Required last executed statement of module. Do not add code below here!


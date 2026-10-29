#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: modis.pm - MODIS coverage metadata used by modisserver.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use modis; # Use my $coverage_metadata_reference = modis::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-04-03 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package modis;

##############################################################################

my @coverage_metadata  = (

  # MOD4 ----------------------------------------------------------------------

  {
    name        => 'MOD4.AOD_550_Dark_Target_Deep_Blue_Combined',
    units       => '-',
    description => 'Combined Dark Target, Deep Blue AOT at 0.55 micron for land and ocean.'
  },
  {
    name        => 'MOD4.AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag',
    units       => '-',
    description => 'Combined Dark Target, Deep Blue AOT at 0.55 micron Algorithm Flag (0=Dark Target, 1=Deep Blue, 2=Mixed)'
  },
  {
    name        => 'MOD4.AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag',
    units       => '-',
    description => 'Combined Dark Target, Deep Blue Aerosol Confidence Flag (0= No Confidence (or fill),1= Marginal, 2= Good, 3= Very Good)'
  },
  {
    name        => 'MOD4.Aerosol_Cldmask_Land_Ocean',
    units       => '-',
    description => 'Aerosol Cloud Mask 500 meter resolution 0 = cloud 1 = clear'
  },
  {
    name        => 'MOD4.Aerosol_Cloud_Fraction_Land',
    units       => '-',
    description => 'Cloud fraction from Land aerosol cloud mask from retrieved and overcast pixels not including cirrus mask'
  },
  {
    name        => 'MOD4.Aerosol_Cloud_Fraction_Ocean',
    units       => '-',
    description => 'Cloud fraction from Land aerosol cloud mask from retrieved and overcast pixels not including cirrus mask'
  },
  {
    name        => 'MOD4.Aerosol_Type_Land',
    units       => '-',
    description => 'Aerosol Type: 1 = Continental, 2 = Moderate Absorption Fine, 3 = Strong Absorption Fine,4 = Weak Absorption Fine, 5 = Dust Coarse'
  },
  {
    name        => 'MOD4.Angstrom_Exponent_1_Ocean',
    units       => '-',
    description => 'Calculated Angstrom Exponent for 0.55 vs 0.86 micron for Average Solution'
  },
  {
    name        => 'MOD4.Angstrom_Exponent_2_Ocean',
    units       => '-',
    description => 'Calculated Angstrom Exponent for 0.86 vs 2.13 micron for Average Solution'
  },
  {
    name        => 'MOD4.Average_Cloud_Pixel_Distance_Land_Ocean',
    units       => 'pixels',
    description => 'Average Distance (number of pixels) to nearest pixel identified as cloudy from each clear pixel used for Aerosol Retrieval in 10 km retrieval box'
  },
  {
    name        => 'MOD4.Cloud_Pixel_Distance_Land_Ocean',
    units       => 'pixels',
    description => 'Distance (number of pixels) to nearest pixel identified as cloudy (500 m resolution)'
  },
  {
    name        => 'MOD4.Corrected_Optical_Depth_Land_wav2p1',
    units       => '-',
    description => 'RetrievedAOT at 2.13 micron'
  },
  {
    name        => 'MOD4.Deep_Blue_Aerosol_Optical_Depth_550_Land',
    units       => '-',
    description => 'AOT at 0.55 micron for landwith all quality data (Quality flag=1,2,3)'
  },
  {
    name        => 'MOD4.Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate',
    units       => '-',
    description => 'Deep Blue AOT at 0.55 micron for land with higher quality data (Quality flag=2,3)'
  },
  {
    name        => 'MOD4.Deep_Blue_Aerosol_Optical_Depth_550_Land_Estimated_Uncertainty',
    units       => '-',
    description => 'Estimated uncertainty (one-sigma confidence envelope) of Deep Blue AOT at 0.55 micron for land for all quality data (Quality flag=1,2,3)'
  },
  {
    name        => 'MOD4.Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag',
    units       => '-',
    description => 'Deep Blue Aerosol Confidence Flag (0= No Confidence (or fill), 1= Marginal, 2= Good, 3= Very Good)'
  },
  {
    name        => 'MOD4.Deep_Blue_Aerosol_Optical_Depth_550_Land_STD',
    units       => '-',
    description => 'Standard deviation of Deep Blue AOT at 0.55 micron for land with all quality data (Quality flag=1,2,3)'
  },
  {
    name        => 'MOD4.Deep_Blue_Algorithm_Flag_Land',
    units       => '-',
    description => 'Deep Blue Aerosol Algorithm Flag (0=DeepBlue, 1=Vegetated, 2=Mixed)'
  },
  {
    name        => 'MOD4.Deep_Blue_Angstrom_Exponent_Land',
    units       => '-',
    description => 'Deep Blue Angstrom Exponent for land with all quality data (Quality flag=1,2,3)'
  },
  {
    name        => 'MOD4.Deep_Blue_Cloud_Fraction_Land',
    units       => '-',
    description => 'Cloud fraction from Deep Blue Aerosol cloud mask over land'
  },
  {
    name        => 'MOD4.Deep_Blue_Number_Pixels_Used_550_Land',
    units       => 'pixels',
    description => 'Number of pixels used for AOT retrieval 0.55 micron for land'
  },
  {
    name        => 'MOD4.Effective_Optical_Depth_0p55um_Ocean',
    units       => '-',
    description => 'Retrieved AOT for average solution at 0.55um For easy L3 processing'
  },
  {
    name        => 'MOD4.Fitting_Error_Land',
    units       => '-',
    description => 'Spectral Fitting error for inversion over land'
  },
  {
    name        => 'MOD4.Glint_Angle',
    units       => 'deg',
    description => 'Glint angle'
  },
  {
    name        => 'MOD4.Image_Optical_Depth_Land_And_Ocean',
    units       => '-',
    description => 'AOT at 0.55 micron for both ocean (Average) and land (corrected) with all quality data (Quality flag=0,1,2,3)'
  },
  {
    name        => 'MOD4.Land_Ocean_Quality_Flag',
    units       => '-',
    description => 'Quality Flag for Land and ocean Aerosol retreivals 0=bad 1= Marginal 2=Good 3=Very Good)'
  },
  {
    name        => 'MOD4.Land_sea_Flag',
    units       => '-',
    description => 'Land_sea_Flag(based on MOD03 Landsea mask 0 = Ocean, 1 = Land and Ephemeral water 2 =Coastal)'
  },
  {
    name        => 'MOD4.Mass_Concentration_Land',
    units       => '1.0e-6g/cm^2',
    description => 'Estimated Column Mass(per area) using assumed mass extinction efficiency'
  },
  {
    name        => 'MOD4.Optical_Depth_Land_And_Ocean',
    units       => '-',
    description => 'AOT at 0.55 micron for both ocean (Average) (Quality flag=1,2,3) and land (corrected) (Quality flag=3)'
  },
  {
    name        => 'MOD4.Optical_Depth_Ratio_Small_Land',
    units       => '-',
    description => 'Fraction of AOT contributed by fine dominated model'
  },
  {
    name        => 'MOD4.Optical_Depth_Ratio_Small_Ocean_0p55micron',
    units       => '-',
    description => 'Fraction of AOT (at 0.55 micron) contributed by fine mode for average solution'
  },
  {
    name        => 'MOD4.Scan_Start_Time',
    units       => 'Seconds_since_1993-01-01',
    description => 'TAI Time at Start of Scan replicated across the swath'
  },
  {
    name        => 'MOD4.Scattering_Angle',
    units       => 'deg',
    description => 'Scattering angle'
  },
  {
    name        => 'MOD4.Sensor_Azimuth',
    units       => 'deg',
    description => 'Sensor azimuth angle, cell to sensor'
  },
  {
    name        => 'MOD4.Sensor_Zenith',
    units       => 'deg',
    description => 'Sensor zenith angle, cell to sensor'
  },
  {
    name        => 'MOD4.Solar_Azimuth',
    units       => 'deg',
    description => 'Solar azimuth angle, cell to sun'
  },
  {
    name        => 'MOD4.Solar_Zenith',
    units       => 'deg',
    description => 'Solar zenith angle, cell to sun'
  },
  {
    name        => 'MOD4.Topographic_Altitude_Land',
    units       => 'km',
    description => 'Averaged topographic altitude (in km) for Land'
  },
  {
    name        => 'MOD4.Wind_Speed_Ncep_Ocean',
    units       => 'm/s',
    description => 'Wind Speed based on NCEP reanalysis for Ocean'
  },


  # MOD43k --------------------------------------------------------------------

  {
    name        => 'MOD43K.Aerosol_Cloud_Fraction_Land',
    units       => '-',
    description => 'Cloud fraction from Land aerosol cloud mask from retrieved and overcast pixels not including cirrus mask'
  },
  {
    name        => 'MOD43K.Aerosol_Cloud_Fraction_Ocean',
    units       => '-',
    description => 'Cloud fraction from land aerosol cloud mask from retrieved and overcast pixels not including cirrus mask'
  },
  {
    name        => 'MOD43K.Aerosol_Type_Land',
    units       => '-',
    description => 'Aerosol Type: 1 = continental, 2 = moderate absorption fine, 3 = strong bbsorption fine, 4 = weak absorption fine, 5 = dust coarse'
  },
  {
    name        => 'MOD43K.Angstrom_Exponent_1_Ocean',
    units       => '-',
    description => 'Calculated angstrom exponent for 0.55 vs 0.86 microns for average solution'
  },
  {
    name        => 'MOD43K.Angstrom_Exponent_2_Ocean',
    units       => '-',
    description => 'Calculated angstrom exponent for 0.86 vs 2.13 microns for average solution'
  },
  {
    name        => 'MOD43K.Corrected_Optical_Depth_Land_wav2p1',
    units       => '-',
    description => 'Retrieved AOT at 2.13 microns'
  },
  {
    name        => 'MOD43K.Fitting_Error_Land',
    units       => '-',
    description => 'Spectral fitting error for inversion over land'
  },
  {
    name        => 'MOD43K.Glint_Angle',
    units       => 'deg',
    description => 'Glint angle'
  },
  {
    name        => 'MOD43K.Image_Optical_Depth_Land_And_Ocean',
    units       => '-',
    description => 'AOT at 0.55 micron for both ocean (average) and land (corrected) with all quality data (quality flag=0,1,2,3)'
  },
  {
    name        => 'MOD43K.Land_Ocean_Quality_Flag',
    units       => '-',
    description => 'Quality flag for land and ocean aerosol retreivals: 0 = bad, 1 = marginal, 2 = good, 3 = very good'
  },
  {
    name        => 'MOD43K.Land_sea_Flag',
    units       => '-',
    description => 'Land_sea_Flag(based on MOD03 Landsea mask: 0 = Ocean, 1 = Land and Ephemeral water, 2 = Coastal'
  },
  {
    name        => 'MOD43K.Mass_Concentration_Land',
    units       => '1.0e-6g/cm^2',
    description => 'Estimated Column Mass(per area) using assumed mass extinction efficiency'
  },
  {
    name        => 'MOD43K.Number_Pixels_Used_Ocean',
    units       => '-',
    description => 'Number of pixels used for ocean retrieval at 865 nm'
  },
  {
    name        => 'MOD43K.Optical_Depth_Land_And_Ocean',
    units       => '-',
    description => 'AOT at 0.55 micron for both ocean (average) (quality flag=1,2,3) and land (corrected) (quality flag=3)'
  },
  {
    name        => 'MOD43K.Optical_Depth_Ratio_Small_Land',
    units       => '-',
    description => 'Fraction of AOT contributed by fine dominated model'
  },
  {
    name        => 'MOD43K.Optical_Depth_Ratio_Small_Ocean_0p55micron',
    units       => '-',
    description => 'Fraction of AOT (at 0.55 microns) contributed by fine mode for average solution'
  },
  {
    name        => 'MOD43K.Scan_Start_Time',
    units       => 'Seconds_since_1993-01-01',
    description => 'TAI Time at Start of Scan replicated across the swath'
  },
  {
    name        => 'MOD43K.Scattering_Angle',
    units       => 'deg',
    description => 'Scattering angle'
  },
  {
    name        => 'MOD43K.Sensor_Azimuth',
    units       => 'deg',
    description => 'Sensor azimuth angle, cell to sensor'
  },
  {
    name        => 'MOD43K.Sensor_Zenith',
    units       => 'deg',
    description => 'Sensor zenith angle, cell to sensor'
  },
  {
    name        => 'MOD43K.Solar_Azimuth',
    units       => 'deg',
    description => 'Sensor azimuth angle, cell to Sun'
  },
  {
    name        => 'MOD43K.Solar_Zenith',
    units       => 'deg',
    description => 'Sensor zenith angle, cell to Sun'
  },
  {
    name        => 'MOD43K.Topographic_Altitude_Land',
    units       => 'km',
    description => 'Terrain height above mean sea level'
  },
  {
    name        => 'MOD43K.Wind_Speed_Ncep_Ocean',
    units       => 'm/s',
    description => 'Wind speed based on NCEP reanalysis for ocean'
  },


  # MOD6 ----------------------------------------------------------------------

  {
    name        => 'MOD6.Cloud_Effective_Emissivity',
    units       => '-',
    description => 'Cloud Effective Emissivity from Cloud Top Pressure Retrieval'
  },
  {
    name        => 'MOD6.Cloud_Effective_Emissivity_Day',
    units       => '-',
    description => 'Cloud Effective Emissivity from Cloud Top Pressure Retrieval, Day Only'
  },
  {
    name        => 'MOD6.Cloud_Effective_Emissivity_Nadir',
    units       => '-',
    description => 'Cloud Effective Emissivity from Cloud Top Pressure Retrieval for Sensor Zenith (View) Angles <= 32 Degrees'
  },
  {
    name        => 'MOD6.Cloud_Effective_Emissivity_Nadir_Day',
    units       => '-',
    description => 'Cloud Effective Emissivity from Cloud Top Pressure Retrieval for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only'
  },
  {
    name        => 'MOD6.Cloud_Effective_Emissivity_Nadir_Night',
    units       => '-',
    description => 'Cloud Effective Emissivity from Cloud Top Pressure Retrieval for Sensor Zenith (View) Angles <= 32 Degrees, Night Data Only'
  },
  {
    name        => 'MOD6.Cloud_Effective_Emissivity_Night',
    units       => '-',
    description => 'Cloud Effective Emissivity from Cloud Top Pressure Retrieval, Night Only'
  },
  {
    name        => 'MOD6.Cloud_Fraction',
    units       => '-',
    description => 'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask'
  },
  {
    name        => 'MOD6.Cloud_Fraction_Day',
    units       => '-',
    description => 'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask, Day Only'
  },
  {
    name        => 'MOD6.Cloud_Fraction_Nadir',
    units       => '-',
    description => 'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask for Sensor Zenith (View) Angles <= 32 Degrees'
  },
  {
    name        => 'MOD6.Cloud_Fraction_Nadir_Day',
    units       => '-',
    description => 'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only'
  },
  {
    name        => 'MOD6.Cloud_Fraction_Nadir_Night',
    units       => '-',
    description => 'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask for Sensor Zenith (View) Angles <= 32 Degrees, Night Data Only'
  },
  {
    name        => 'MOD6.Cloud_Fraction_Night',
    units       => '-',
    description => 'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask, Night Only'
  },
  {
    name        => 'MOD6.Cloud_Height_Method',
    units       => '-',
    description => 'Index Indicating MODIS Bands Used for Cloud Top Pressure Retrieval'
  },
  {
    name        => 'MOD6.Cloud_Phase_Infrared',
    units       => '-',
    description => 'Cloud Phase from 8.5 and 11 um Bands'
  },
  {
    name        => 'MOD6.Cloud_Phase_Infrared_Day',
    units       => '-',
    description => 'Cloud Phase from 8.5 and 11 um Bands, Day Only'
  },
  {
    name        => 'MOD6.Cloud_Phase_Infrared_Night',
    units       => '-',
    description => 'Cloud Phase from 8.5 and 11 um Bands, Night Only'
  },
  {
    name        => 'MOD6.Cloud_Top_Height',
    units       => 'm',
    description => 'Geopotential Height at Retrieved Cloud Top Pressure Level (rounded to nearest 50m)'
  },
  {
    name        => 'MOD6.Cloud_Top_Height_Nadir',
    units       => 'm',
    description => 'Geopotential Height at Retrieved Cloud Top Pressure Level for Sensor Zenith (View) Angles <=32 Degrees (rounded to nearest 50m)'
  },
  {
    name        => 'MOD6.Cloud_Top_Height_Nadir_Day',
    units       => 'm',
    description => 'Geopotential Height at Retrieved Cloud Top Pressure Level for Sensor Zenith (View) Angles <=32 Degrees, Day Data Only (rounded to nearest 50m)'
  },
  {
    name        => 'MOD6.Cloud_Top_Height_Nadir_Night',
    units       => 'm',
    description => 'Geopotential Height at Retrieved Cloud Top Pressure Level for Sensor Zenith (View) Angles <=32 Degrees, Night Data Only (rounded to nearest 50m)'
  },
  {
    name        => 'MOD6.Cloud_Top_Pressure',
    units       => 'hPa',
    description => 'Cloud Top Pressure Level (rounded to nearest 5 hPa)'
  },
  {
    name        => 'MOD6.Cloud_Top_Pressure_Day',
    units       => 'hPa',
    description => 'Cloud Top Pressure Level, Day Data Only (rounded to nearest 5 hPa)'
  },
  {
    name        => 'MOD6.Cloud_Top_Pressure_Infrared',
    units       => 'hPa',
    description => 'Cloud Top Pressure from IR Window Retrieval'
  },
  {
    name        => 'MOD6.Cloud_Top_Pressure_Nadir',
    units       => 'hPa',
    description => 'Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees (rounded to nearest 5 hPa)'
  },
  {
    name        => 'MOD6.Cloud_Top_Pressure_Nadir_Day',
    units       => 'hPa',
    description => 'Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees (rounded to nearest 5 hPa), Day Data Only'
  },
  {
    name        => 'MOD6.Cloud_Top_Pressure_Nadir_Night',
    units       => 'hPa',
    description => 'Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees (rounded to nearest 5 hPa), Night Data Only'
  },
  {
    name        => 'MOD6.Cloud_Top_Pressure_Night',
    units       => 'hPa',
    description => 'Cloud Top Pressure Level, Night Data Only (rounded to nearest 5 hPa)'
  },
  {
    name        => 'MOD6.Cloud_Top_Temperature',
    units       => 'K',
    description => 'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level'
  },
  {
    name        => 'MOD6.Cloud_Top_Temperature_Day',
    units       => 'K',
    description => 'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level, Day Only'
  },
  {
    name        => 'MOD6.Cloud_Top_Temperature_Nadir',
    units       => 'K',
    description => 'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees'
  },
  {
    name        => 'MOD6.Cloud_Top_Temperature_Nadir_Day',
    units       => 'K',
    description => 'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only'
  },
  {
    name        => 'MOD6.Cloud_Top_Temperature_Nadir_Night',
    units       => 'K',
    description => 'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level for Sensor Z enith (View) Angles <= 32 Degrees, NightData Only'
  },
  {
    name        => 'MOD6.Cloud_Top_Temperature_Night',
    units       => 'K',
    description => 'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level, Night Only'
  },
  {
    name        => 'MOD6.Radiance_Variance',
    units       => 'W/m2/sr/micron',
    description => 'Band 31 Radiance Standard Deviation'
  },
  {
    name        => 'MOD6.Scan_Start_Time',
    units       => 'Seconds_since_1993-01-01',
    description => 'TAI Time at Start of Scan replicated across the swath'
  },
  {
    name        => 'MOD6.Sensor_Azimuth',
    units       => 'deg',
    description => 'Sensor Azimuth Angle, Cell to Sensor'
  },
  {
    name        => 'MOD6.Sensor_Azimuth_Day',
    units       => 'deg',
    description => 'Sensor Azimuth Angle, Cell to Sensor, Day Data Only'
  },
  {
    name        => 'MOD6.Sensor_Azimuth_Night',
    units       => 'deg',
    description => 'Sensor Azimuth Angle, Cell to Sensor, Night Data Only'
  },
  {
    name        => 'MOD6.Sensor_Zenith',
    units       => 'deg',
    description => 'Solar Zenith Angle, Cell to Sun'
  },
  {
    name        => 'MOD6.Sensor_Zenith_Day',
    units       => 'deg',
    description => 'Solar Zenith Angle, Cell to Sun, Day Data Only'
  },
  {
    name        => 'MOD6.Sensor_Zenith_Night',
    units       => 'deg',
    description => 'Solar Zenith Angle, Cell to Sun, Night Data Only'
  },
  {
    name        => 'MOD6.Solar_Azimuth',
    units       => 'deg',
    description => 'Solar Azimuth Angle, Cell to Sun'
  },
  {
    name        => 'MOD6.Solar_Azimuth_Day',
    units       => 'deg',
    description => 'Solar Azimuth Angle, Cell to Sun, Day Data Only'
  },
  {
    name        => 'MOD6.Solar_Azimuth_Night',
    units       => 'deg',
    description => 'Solar Azimuth Angle, Cell to Sun, Night Data Only'
  },
  {
    name        => 'MOD6.Solar_Zenith',
    units       => 'deg',
    description => 'Sensor Zenith Angle, Cell to Sensor'
  },
  {
    name        => 'MOD6.Solar_Zenith_Day',
    units       => 'deg',
    description => 'Sensor Zenith Angle, Cell to Sensor, Day Data Only'
  },
  {
    name        => 'MOD6.Solar_Zenith_Night',
    units       => 'deg',
    description => 'Sensor Zenith Angle, Cell to Sensor, Night Data Only'
  },
  {
    name        => 'MOD6.Surface_Pressure',
    units       => 'hPa',
    description => 'Surface Pressure from Ancillary Data'
  },
  {
    name        => 'MOD6.Surface_Temperature',
    units       => 'K',
    description => 'Surface Temperature from Ancillary Data'
  },
  {
    name        => 'MOD6.Tropopause_Height',
    units       => 'hPa',
    description => 'Tropopause Height from Ancillary Data'
  },


  # MOD7 ----------------------------------------------------------------------

  {
    name        => 'MOD7.Cloud_Mask',
    units       => '-',
    description => 'MODIS Cloud Mask, First Byte'
  },
  {
    name        => 'MOD7.K_Index',
    units       => 'K',
    description => 'K Index'
  },
  {
    name        => 'MOD7.Lifted_Index',
    units       => 'K',
    description => 'Lifted Index'
  },
  {
    name        => 'MOD7.Scan_Start_Time',
    units       => 'Seconds_since_1993-01-01',
    description => 'TAI Time at Start of Scan replicated across the swath'
  },
  {
    name        => 'MOD7.Sensor_Azimuth',
    units       => 'deg',
    description => 'Sensor Azimuth Angle, Cell to Sensor'
  },
  {
    name        => 'MOD7.Sensor_Zenith',
    units       => 'deg',
    description => 'Sensor Zenith Angle, Cell to Sensor'
  },
  {
    name        => 'MOD7.Skin_Temperature',
    units       => 'K',
    description => 'Skin Temperature'
  },
  {
    name        => 'MOD7.Solar_Azimuth',
    units       => 'deg',
    description => 'Solar Azimuth Angle, Cell to Sun'
  },
  {
    name        => 'MOD7.Solar_Zenith',
    units       => 'deg',
    description => 'Solar Zenith Angle, Cell to Sun'
  },
  {
    name        => 'MOD7.Surface_Elevation',
    units       => 'm',
    description => 'Surface Elevation'
  },
  {
    name        => 'MOD7.Surface_Pressure',
    units       => 'hPa',
    description => 'Surface Pressure'
  },
  {
    name        => 'MOD7.Total_Ozone',
    units       => 'Dobson',
    description => 'Total Ozone Burden'
  },
  {
    name        => 'MOD7.Total_Totals',
    units       => 'K',
    description => 'Total Totals'
  },
  {
    name        => 'MOD7.Tropopause_Height',
    units       => 'hPa',
    description => 'Tropopause Height'
  },
  {
    name        => 'MOD7.Water_Vapor',
    units       => 'cm',
    description => 'Total Column Precipitable Water Vapor - IR Retrieval'
  },
  {
    name        => 'MOD7.Water_Vapor_Direct',
    units       => 'cm',
    description => 'Total Column Precipitable Water Vapor - Direct IR Retrieval'
  },
  {
    name        => 'MOD7.Water_Vapor_High',
    units       => 'cm',
    description => 'Precipitable Water Vapor High - IR Retrieval'
  },
  {
    name        => 'MOD7.Water_Vapor_Low',
    units       => 'cm',
    description => 'Precipitable Water Vapor Low - IR Retrieval'
  }

);

###############################################################################

sub get_coverage_metadata_reference { return \@coverage_metadata; }

1; # Required last executed statement of module. Do not add code below here!


#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: calipso.pm - CALIPSO coverage metadata used by calipsoserver.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use calipso; # Use my $coverage_metadata_reference = calipso::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-04-05 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package calipso;

##############################################################################

my @coverage_metadata  = (

  # L1 -----------------------------------------------------------------------

  {
    name        => 'L1.Amplifier_Gain_1064',
    units       => 'V/V',
    description => 'Gain of the variable gain amplifier for the 1064nm parallel channel.'
  },
  {
    name        => 'L1.Attenuated_Backscatter_1064',
    units       => '/km/sr',
    description => 'Attenuated backscatter at 1064nm, derived from the calibrated, range-corrected, gain and laser energy normalized, background and baseline subtracted lidar return signal. The first 33 bins of each profile contain fill values (-9999), because no 1064nm data is downlinked from the 30.1-40km altitude range.'
  },
  {
    name        => 'L1.Calibration_Constant_1064',
    units       => 'km3*sr*counts',
    description => 'Calibration constant at 1064 nm reported for each single shot profile. CALIOP 10 64nm channel is calibrated relative to the 532nm channel using the backscatter from a carefully selected subset of cirrus cloud measurements. Unlike the 532nm calibrations, at 1064nm the same procedure is used for both daytime and nighttime calibrations.'
  },
  {
    name        => 'L1.Calibration_Constant_532',
    units       => 'km3*sr*count/J',
    description => 'Calibration constant at 532 nm reported for each single shot profile using a molecular normalization technique. Separate techniques are applied for day and night measurements.'
  },
  {
    name        => 'L1.Calibration_Constant_Uncertainty_1064',
    units       => 'km3*sr*counts',
    description => 'Estimated uncertainty in the 1064nm calibration constant, computed and reported for each single shot profile. Unlike the 532nm calibration constant uncertainties, the uncertainties reported for the 1064nm calibration constants include contributions from both random errors, which can be reduced by averaging, and systematic biases, which cannot be reduced by averaging.'
  },
  {
    name        => 'L1.Calibration_Constant_Uncertainty_532',
    units       => 'km3*sr*count/J',
    description => 'Estimated calibration uncertainty due to random noise in the 532 nm calibration constant, computed and reported for each single shot profile. Random uncertainties for calibration constants subsequently derived using a multi-orbit sliding window averaging are reduced from the single, individual sample estimate by the square root of the number of samples included in the average.'
  },
  {
    name        => 'L1.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L1.Depolarization_Gain_Ratio_532',
    units       => '-',
    description => 'Ratio of the opto-electric gains between the 532nm perpendicular and 532nm parallel channels, i.e. the polarization gain ratio (PGR). The 532 nm nighttime PGR is calculated from routine PGR calibration operations, in which a pseudo-depolarizer is inserted into the optical path to generate equal backscatter intensities in both the 532nm channels. The 532nm daytime PGR is calculated based on the Perpendicular and Parallel RMS Baseline measurements above optically thick ice clouds.'
  },
  {
    name        => 'L1.Depolarization_Gain_Ratio_Uncertainty_532',
    units       => '-',
    description => 'Uncertainty in the 532 nm depolarization gain ratio (PGR) due to random noise computed based on the 532nm parallel and perpendicular noise scale factors. Systematic errors are not included but are estimated to be small.'
  },
  {
    name        => 'L1.Earth-Sun_Distance',
    units       => 'AU',
    description => 'The distance from the Earth surface to the Sun.'
  },
  {
    name        => 'L1.Frame_Number',
    units       => '-',
    description => 'The number of the current frame (5km horizonal segment) within the sequence of 11 frames within a Payload Data Acquisition Cycle (PDAC).'
  },
  {
    name        => 'L1.GMAO_Surface_Elevation',
    units       => '-',
    description => 'Surface elevation at the laser footprint, above local mean sealeve, obtained ancillary GMAQ data. This is at a lower resolution as compared with Surface_Elevation.'
  },
  {
    name        => 'L1.IGBP_Surface_Type',
    units       => '-',
    description => 'IGBP surface type at the laser footprint provided by the CERES/SARB surface map: 1 = Evergreen-Needleleaf-Forest, 2 = Evergreen-Broadleaf-Forest, 3 = Deciduous-Needleleaf-Forest, 4 = Deciduous-Broadleaf-Forest, 5 = Mixed-Forest, 6 = Closed-Shrubland, 7 = Open-Shrubland/Desert, 8 = Woody-Savanna, 9 = Savanna, 10 = Grassland, 11 = Wetland, 12 = Cropland, 13 = Urban, 14 = Crop-Mosaic, 15 = Permanent-Snow, 16 = Barren/Desert, 17 = Water, 18 = Tundra.'
  },
  {
    name        => 'L1.Land_Water_Mask',
    units       => '-',
    description => 'Surface type at the laser footprint provided by the CloudSat science team digital elevation model (DEM): 0 = shallow ocean, 1 land, 2 = coastlines, 3 = shallow inland water, 4 = intermittent water, 5 = deep inland water, 6 = continental ocean, 7 = deep ocean.'
  },
  {
    name        => 'L1.Laser_Energy_1064',
    units       => 'J',
    description => 'The laser energy at 1064nm measured by the laser energy monitor for each shot.'
  },
  {
    name        => 'L1.Laser_Energy_532',
    units       => 'J',
    description => 'The laser energy at 532 nm measured by the laser energy monitor for each shot.'
  },
  {
    name        => 'L1.Lidar_Mode',
    units       => '-',
    description => 'Lidar operations mode: 1 = safe, 2 = standby, 3 = autonomous data collection, 4 = built-in test system (BITS), 5 = special operations, 6 = diagnostic.'
  },
  {
    name        => 'L1.Lidar_Submode',
    units       => '-',
    description => 'Lidar operations submode (value - submode): 1 - safe, 2 - warm-up, 3 - standby, 4 - normal, 5 - BITS profile #1, no background, 6 - BITS profile #2, no background, 7 - BITS profile #3, no background, 8 - BITS profile #4, no background, 9 - BITS profile #1, with background,  10 - BITS profile #2, with background, 11 - BITS profile #3, with background, 12 - BITS profile #4, with background, 13 - no baseline removal, 14 - baseline measurement, 15 - polarization calibration, 16 - extended background measurement,  17 - boresight alignment, 18 - manual boresight operation, 19 - boresight search, 20 - build on-board baseline profiles, 21 - baseline corrected high gain channel, 22 - baseline corrected low gain channel, 23 - uncorrected high gain channel, 24 - uncorrected low gain channel, 25 - diagnostic.'
  },
  {
    name        => 'L1.Molecular_Number_Density',
    units       => 'molecules/m3',
    description => 'Mean molecular number density reported for the midpoint of each range bin in the profile; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L1.Negative_Signal_Anomaly_Index_1064',
    units       => '-',
    description => 'Reports the zero-based altitude index of the lowest range bin where the negative signal anomaly is detected in the 1064nm channel. For profiles where the negative signal anomaly is not detected, as value of -1 is reported.'
  },
  {
    name        => 'L1.Negative_Signal_Anomaly_Index_532Par',
    units       => '-',
    description => 'Reports the zero-based altitude index of the lowest range bin where the negative signal anomaly is detected in the 532nm parallel channel. For profiles where the negative signal anomaly is not detected, as value of -1 is reported.'
  },
  {
    name        => 'L1.Negative_Signal_Anomaly_Index_532Perp',
    units       => '-',
    description => 'Reports the zero-based altitude index of the lowest range bin where the negative signal anomaly is detected in the 532nm perpendicular channel. For profiles where the negative signal anomaly is not detected, as value of -1 is reported.'
  },
  {
    name        => 'L1.Noise_Scale_Factor_1064',
    units       => 'square_root_of_counts',
    description => 'Noise scale factor (NSF) for the 1064nm channel. CALIOP does not measure the background signal level at 1064 nm, because the APD detector dark noise is dominant during both nighttime and daytime measurement. For this reason, the procedure to estimate the NSF for the 532 nm channels cannot be used for the 1064 nm channel. The 1064nm NSF is therefore set to 0, which causes negligible error because, as above, the APD detector dark noise is the dominant error source.'
  },
  {
    name        => 'L1.Noise_Scale_Factor_532_Parallel',
    units       => 'square_root_of_counts',
    description => 'Noise scale factor (NSF) for each shot of the 532nm parallel channel, computed from daytime measurement of the Parallel RMS Baseline 532 and the Parallel Background Monitor 532.'
  },
  {
    name        => 'L1.Noise_Scale_Factor_532_Perpendicular',
    units       => 'square_root_of_counts',
    description => 'Noise scale factor (NSF) for each shot of the 532nm perpendicular channel, computed from daytime measurement of the Perpendicular RMS Baseline 532 and the Perpendicular Background Monitor 532.'
  },
  {
    name        => 'L1.Number_Bins_Shift',
    units       => '-',
    description => 'Number of 30 meter bins shifted by the altitude registration algorithm to match the range bin at mean sea level in the measured profile to the range bin at mean sealevel in the fixed altitude grid; positive = shift up and negative = shift down.'
  },
  {
    name        => 'L1.Off_Nadir_Angle',
    units       => 'degrees',
    description => 'Angle of the lidar viewing vector; 0.3 degrees off-nadir from launch/first-light to 3.0 degrees after November 28, 2007.'
  },
  {
    name        => 'L1.Ozone_Number_Density',
    units       => 'molecules/m3',
    description => 'Ozone number density, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L1.Parallel_Amplifier_Gain_532',
    units       => 'V/V',
    description => 'Gain of the variable gain amplifier for the 532 nm parallel channel.'
  },
  {
    name        => 'L1.Parallel_Background_Monitor_532',
    units       => 'counts',
    description => 'Amplitude of the background signal measured using the 532 nm parallel channel background monitor and rescaled to science digitizer units.'
  },
  {
    name        => 'L1.Parallel_Column_Reflectance_532',
    units       => '-',
    description => 'Bi-directional column reflectance derived from the root-mean-square variation of the 532 nm parallel channel baseline measurements.'
  },
  {
    name        => 'L1.Parallel_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Column reflectance uncertainty for the 532nm parallel channel.'
  },
  {
    name        => 'L1.Parallel_RMS_Baseline_532',
    units       => 'counts',
    description => 'Root mean square (RMS) noise, in science digitizer counts, of the high-altitude baseline signal measured for each laser pulse over 1000 15m samples acquired in the 64-80km altitude range in the 532nm parallel channel.'
  },
  {
    name        => 'L1.Perpendicular_Amplifier_Gain_532',
    units       => 'V/V',
    description => 'Gain of the variable gain amplifier for the 532 nm perpendicular channel.'
  },
  {
    name        => 'L1.Perpendicular_Attenuated_Backscatter_532',
    units       => '/km/sr',
    description => 'Perpendicular channel component of the 532nm total attenuated backscatter coefficients. Profiles of the parallel channel component of the 532nm total attenuated backscatter coefficients can be obtained by subtraction of the perpendicular from the total.'
  },
  {
    name        => 'L1.Perpendicular_Background_Monitor_532',
    units       => 'counts',
    description => ' Amplitude of the background signal measured using the 532nm perpendicular channel background monitor and rescaled to science digitizer units.'
  },
  {
    name        => 'L1.Perpendicular_Column_Reflectance_532',
    units       => '-',
    description => 'Bi-directional column reflectance derived from the root-mean-square variation of the 532 nm perpendicular channel baseline measurements.'
  },
  {
    name        => 'L1.Perpendicular_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Column reflectance uncertainty for the 532nm perpendicular channel.'
  },
  {
    name        => 'L1.Perpendicular_RMS_Baseline_532',
    units       => 'counts',
    description => 'Root mean square (RMS) noise, in science digitizer counts, of the high-altitude baseline signal measured for each laser pulse over 1000 15m samples acquired in the 64-80km altitude range in the 532nm perpendicular channel.'
  },
  {
    name        => 'L1.Pressure',
    units       => 'hPa',
    description => 'Pressure, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L1.Profile_ID',
    units       => '-',
    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
  },
  {
    name        => 'L1.Profile_Time',
    units       => 'seconds_since_1993-01-01',
    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
  },
  {
    name        => 'L1.QC_Flag',
    units       => '-',
    description => 'Quality Control Flag #1 (bit - interpretation): bit 1 - 532nm parallel channel missing bit, 2 - 532nm perpendicular channel missing, bit 3 - 1064nm channel missing, bit 4 - Not geolocated, bit 5 - Single shot 532 laser energy below calibration threshold (below 0.01 J; near zero energy), bit 6 - Single shot 1064 laser energy below calibration threshold (below 0.0408 J; near zero energy), bit 7 - Historical value used for the depolarization gain ratio, bit 8 - Historical calibration constant used, 532nm parallel channel, bit 9 - Historical calibration constant used, 532nm perpendicular channel, bit 10 - Historical calibration constant used, 1064nm channel, bit 11 - Averaged calibration constant used, 532nm parallel channel, bit 1 2 - Averaged calibration constant used, 532nm perpendicular channel, bit 1 3 - Single shot 532 laser energy below data quality threshold (below 0.05 J; low energ y), bit 14 - Single shot 1064 laser energy below data quality threshold (below 0.05J; low energy), bit 15 - Near zero 532nm laser energy profile included in region 3 average, bit 16 - Near zero 1064nm laser energy profile inclu ded in region 3 average, bit 17 - Near zero 532nm laser energy profile included in region 4 average, bit 18 - Near zero 1064nm laser energy profile included in region 4 average, bit 19 - Near zero 532nm laser energy profile included in region 5 average, bit 20 - Low 532 nm laser energy profile included in region 3 average, bit 21 - Low 1064nm laser energy profile included in region 3 average, bit 22 - Low 532nm laser energy profile included in region 4 average, bit 23 - Low 532nm laser energy profile included in region 4 average, bit 24 - Low 532nm laser energy profile included in region 5 average, bits 25...32 - Spare.'
  },
  {
    name        => 'L1.QC_Flag_2',
    units       => '-',
    description => 'Quality Control Flag #2 (bit - interpretation):, bit 1 - Reserve, bit 2 - Excessive digitizer underflows, 532nm parallel channel in region 6, bit 3 - Excessive digitizer underflows, 532nm perpendicular parallel channel, region 6, bit 4 - Excessive digitizer underflows, 1064nm channel, region 6, bit 5 - Excessive digitizer overflows, 532nm parallel channel, region 6, bit 6 - Excessive digitizer overflows, 532nm perpendicular parallel channel, region 6, bit 7 - Excessive digitizer overflows, 1064nm channel, region 6, bit 8 - Excessive digitizer overflows, 532nm parallel channel, region 2, bit 9 - Excessive digitizer overflows, 532nm perpendicular parallel channel, region 2, bit 10 - Excessive digitizer overflows, 1064nm channel, region 2, bit 11 - LRE Flags in SAD packet indicate bad data, 532nm parallel channel, bit 12 - LRE Flags in SAD packet indicate bad data, 532nm perpendicular channel, bit 13 - LRE Flags in SAD packet indicate bad data, 1064nm channel, bit 14 - Quality Flags in SAD packet indicate bad data, 532nm parallel channel, bit 15 - Quality Flags in SAD packet indicate bad data, 532nm perpendicular channel, bit 16 - Quality Flags in SAD packet indicate bad data, 1064nm channel, bit 17 - Suspicious offset calculation, 532nm parallel channel, bit 18 - Suspicious offset calculation, 532nm perpendicular channel, bit 19 - Suspicious offset calculation, 1064nm channel, bit 20 - Suspicious mean signal value, 532nm parallel channel (any/all regions), bit 21 - Suspicious mean signal value, 532nm perpendicular channel (any/all regions), bit 22 - Suspicious mean signal value, 1064nm channel (any/all regions), bit 23 - RMS noise out of range, 532nm parallel channel, bit 24 - RMS noise out of range, 532nm perpendicular parallel channel, bit 25 - RMS noise out of range, 1064nm channel, bit 26 - Near surface meteorological parameters were remapped to DEM surface altitude, bit 27 - 1064nm calibration coefficients suspect due to low temperatures associated with power up, bits 28 ...32 Spare.'
  },
  {
    name        => 'L1.Reflectance_Calibration_Coefficient_532',
    units       => 'W/m^2/sr/counts',
    description => 'Calibration coefficient that converts the sum of the squares of the RMS baseline measurements (units = CALIOP science digitizer counts) to radiances (units = W * 1/ m^2 * 1/sr).'
  },
  {
    name        => 'L1.Reflectance_Calibration_Coefficient_Uncertainty_532',
    units       => 'W/m^2/sr/counts',
    description => 'Absolute uncertainty in the calibration coefficients that convert the sum of the squares of the RMS baseline measurements to radiances.'
  },
  {
    name        => 'L1.Relative_Humidity',
    units       => '%',
    description => 'Relative humidity, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L1.RMS_Baseline_1064',
    units       => 'counts',
    description => 'Root mean square (RMS) noise, in science digitizer counts, of the high-altitude baseline signal measured for each laser pulse over 1000 15m samples acquired in the 65-80km altitude range in the 1064nm channel.'
  },
  {
    name        => 'L1.Scattering_Angle',
    units       => 'degrees',
    description => 'Angle between the lidar viewing angle and the line of sight of the sun.'
  },
  {
    name        => 'L1.Snow_Ice_Surface_Type',
    units       => '-',
    description => 'Snow and ice coverage for the surface at the laser footprint provided by the U.S. Airforce Weather Agency (AFWA). Values between 0 and 100 represent the percentage of ice, while 101 = permanent ice, 103 = snow, and 255 = mixed pixels at coastline.'
  },
  {
    name        => 'L1.Solar_Azimuth_Angle',
    units       => 'degrees',
    description => 'Azimuth angle from north of the line of sight to the sun.'
  },
  {
    name        => 'L1.Solar_Zenith_Angle',
    units       => 'degrees',
    description => 'Angle between the zenith at the lidar footprint on the surface and the line of sight to the sun.'
  },
  {
    name        => 'L1.Spacecraft_Altitude',
    units       => 'km',
    description => 'Altitude above mean sea level of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Attitude_X',
    units       => 'degrees',
    description => 'X-coordinate of attitude (roll, pitch, yaw) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Attitude_Y',
    units       => 'degrees',
    description => 'Y-coordinate of attitude (roll, pitch, yaw) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Attitude_Z',
    units       => 'degrees',
    description => 'Z-coordinate of attitude (roll, pitch, yaw) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Attitude_Rate_X',
    units       => 'deg/s',
    description => 'X-coordinate of change of attitude of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Attitude_Rate_Y',
    units       => 'deg/s',
    description => 'Y-coordinate of change of attitude of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Attitude_Rate_Z',
    units       => 'deg/s',
    description => 'Z-coordinate of change of attitude of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Position_X',
    units       => 'km',
    description => 'X-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Position_Y',
    units       => 'km',
    description => 'Y-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Position_Z',
    units       => 'km',
    description => 'Z-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Velocity_X',
    units       => 'km/s',
    description => 'X-coordinate of velocity (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Velocity_Y',
    units       => 'km/s',
    description => 'Y-coordinate of velocity (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L1.Spacecraft_Velocity_Z',
    units       => 'km/s',
    description => 'Z-coordinate of velocity (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L1.Subsatellite_Latitude',
    units       => 'degrees',
    description => 'Latitude of the geodetic subsolar point, which is the point on the surface where the geodetic zenith vector points toward the satellite.'
  },
  {
    name        => 'L1.Subsatellite_Longitude',
    units       => 'degrees',
    description => 'Longitude of the geodetic subsolar point, which is the point on the surface where the geodetic zenith vector points toward the satellite.'
  },
  {
    name        => 'L1.Subsolar_Latitude',
    units       => 'degrees',
    description => 'Latitude of the geodetic subsolar point, which is the point on the surface where the geodetic zenith vector points toward the Sun.'
  },
  {
    name        => 'L1.Subsolar_Longitude',
    units       => 'degrees',
    description => 'Longitude of the geodetic subsolar point, which is the point on the surface where the geodetic zenith vector points toward the Sun.'
  },
  {
    name        => 'L1.Surface_Altitude_Shift',
    units       => 'km',
    description => 'Altitude difference between the profile specific 30 meter altitude array and the fixed 30 meter altitude array at the array element that includes mean sea level. Profile specific altitude arrays are computed as a function of the actual spacecraft off-nadir angle, which varies slightly from the commanded spacecraft off-nadir angle. The difference is calculated as: Surface_Altitude_Shift = altitude (profile specific 30 meter mean sea level bin) - altitude (fixed 30 meter mean sea level bin); positive = shift up and negative = shift down.'
  },
  {
    name        => 'L1.Surface_Elevation',
    units       => 'km',
    description => 'Surface elevation at the laser footprint, above local mean sealeve, obtained from the CloudSat science team digital elevation model (DEM).'
  },
  {
    name        => 'L1.Surface_Saturation_Flag_1064',
    units       => '-',
    description => 'Reports the likelihood that the surface signal in the 1064nm channel is saturated: 0 - not saturated, 1 - possibly saturated, 2 - certainly saturated.'
  },
  {
    name        => 'L1.Surface_Saturation_Flag_532Par',
    units       => '-',
    description => 'Reports the likelihood that the surface signal in the 532 nm parallel channel is saturated: 0 - not saturated, 1 - possibly saturated, 2 - certainly saturated.'
  },
  {
    name        => 'L1.Surface_Saturation_Flag_532Perp',
    units       => '-',
    description => 'Reports the likelihood that the surface signal in the 532nm perpendicular channel is saturated: 0 - not saturated, 1 - possibly saturated 2 - certainly saturated.'
  },
  {
    name        => 'L1.Surface_Saturation_Index_1064',
    units       => '-',
    description => 'Reports the zero-based altitude index of the lowest range bin where the surface saturation flag is either possibly saturated or certainly saturated for the 1064nm channel. For profiles classified as not saturated, a value of -1 is reported.'
  },
  {
    name        => 'L1.Surface_Saturation_Index_532Par',
    units       => '-',
    description => 'Reports the zero-based altitude index of the lowest range bin where the surface saturation flag is either possibly saturated or certainly saturated for the 532nm parallel channel. For profiles classified as not saturated, a value of -1 is reported.'
  },
  {
    name        => 'L1.Surface_Saturation_Index_532Perp',
    units       => '-',
    description => 'Reports the zero-based altitude index of the lowest range bin where the surface saturation flag is either possibly saturated of certainly saturated for the 532nm perpendictular channel. For profiles classified as not saturated, a value of -1 is reported.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L1.Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Surface wind speeds measured 10m above the Earth surface. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v. Provided by MERRA-2.'
#  },
  {
    name        => 'L1.Temperature',
    units       => 'degC',
    description => 'Temperature, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L1.Total_Attenuated_Backscatter_532',
    units       => '/km/sr',
    description => 'Total attenuated backscatter coefficientsat 532nm, equal to the sum of the parallel and perpendicular channel 532nm attenuated backscatter coefficients. Attenuated backscatter coefficients are derived from the calibrated, range-corrected, gain and laser energy normalized, background and baseline subtracted lidar return signal.'
  },
  {
    name        => 'L1.Tropopause_Height',
    units       => 'km',
    description => 'Tropopause height. Provided by MERRA-2.'
  },
  {
    name        => 'L1.Tropopause_Temperature',
    units       => 'degC',
    description => 'Tropopause temperature. Provided by MERRA-2.'
  },
  {
    name        => 'L1.Viewing_Azimuth_Angle',
    units       => 'degrees',
    description => 'Azimuth angle from north of the lidar viewing vector.'
  },
  {
    name        => 'L1.Viewing_Zenith_Angle',
    units       => 'degrees',
    description => 'Angle between the lidar viewing vector and the zenith at the lidar footprint on the surface.'
  },



  # L2_01kmCLay ---------------------------------------------------------------

  {
    name        => 'L2_01kmCLay.Column_IAB_Cumulative_Probability',
    units       => '-',
    description => 'The cumulative probability of measuring a total column integrated attenuated backscatter value equal to the value computed for the current profile.'
  },
  {
    name        => 'L2_01kmCLay.Column_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'The integral with respect to altitude of the 532nm total attenuated backscatter coefficients. The limits of integration are from the onset of the backscatter signal at ~40 km, down to the range bin immediately prior to the surface elevation specified by the digital elevation map (DEM).'
  },
  {
    name        => 'L2_01kmCLay.Column_Particulate_Optical_Depth_Above_Opaque_Water_Cloud_532',
    units       => '-',
    description => '532nm particulate optical depths in the vertical column above an opaque water cloud, derived according to the method detailed in Hu et al., 2007 (DOI: 10.1109/LGRS.2007.901085).'
  },
  {
    name        => 'L2_01kmCLay.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L2_01kmCLay.Feature_Classification_Flags',
    units       => '-',
    description => 'For each layer detected in the CALIPSO backscatter data a feature classification flags are derived. The complete set of flags are stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): 011 bits 1-3; Feature Type, 0110 = invalid (bad or missing data), 0111 = clear air, 0112 = cloud, 0113 = tropospheric aerosol, 0114 = stratospheric aerosol, 0115 = surface, 0116 = subsurface, 0117 = no signal, bits 4-5; Feature Type QA, 011 = 0 = none, 0111 = low, 0112 = medium, 0113 = high 011 bits 6-7; Ice/Water Phase (cloud), 0110 = unknown/not determined, 0111 = ice, 0112 = water, 0113 = oriented ice crystals, bits 8-9; Ice/Water Phase QA (cloud), 0110 = none, 0111 = low, 0112 = medium, 0113 = high, 011bits 10-12; Feature Subtype (clear-air), 0110 = n/a, 0111 = not searched for features at 80km resolution, 0112 = not searched for features at 20 and 8km resolution, bits 10-12; Feature Subtype (Cloud), 0110 = low overcast, transparent, 0111 = low overcast, opaque, 0112 = transition stratocumulus, 0113 = low, broken cloud.'
  },
  {
    name        => 'L2_01kmCLay.IGBP_Surface_Type',
    units       => '-',
    description => 'IGBP surface type at the laser footprint provided by the CERES/SARB surface map: 1 = Evergreen-Needleleaf-Forest, 2 = Evergreen-Broadleaf-Forest, 3 = Deciduous-Needleleaf-Forest, 4 = Deciduous-Broadleaf-Forest, 5 = Mixed-Forest, 6 = Closed-Shrubland, 7 = Open-Shrubland/Desert, 8 = Woody-Savanna, 9 = Savanna, 10 = Grassland, 11 = Wetland, 12 = Cropland, 13 = Urban, 14 = Crop-Mosaic, 15 = Permanent-Snow, 16 = Barren/Desert, 17 = Water, 18 = Tundra.'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Attenuated_Backscatter_1064',
    units       => '1/sr',
    description => 'Integrated_Attenuated_Backscatter_1064'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'Integrated_Attenuated_Backscatter_532'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Attenuated_Backscatter_Uncertainty_1064',
    units       => '1/sr',
    description => 'Integrated_Attenuated_Backscatter_Uncertainty_1064'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Attenuated_Backscatter_Uncertainty_532',
    units       => '1/sr',
    description => 'Integrated_Attenuated_Backscatter_Uncertainty_532'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Attenuated_Total_Color_Ratio',
    units       => '-',
    description => 'Integrated_Attenuated_Total_Color_Ratio'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Attenuated_Total_Color_Ratio_Uncertainty',
    units       => '-',
    description => 'Integrated_Attenuated_Total_Color_Ratio_Uncertainty'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Volume_Depolarization_Ratio',
    units       => '-',
    description => 'Integrated_Volume_Depolarization_Ratio'
  },
  {
    name        => 'L2_01kmCLay.Integrated_Volume_Depolarization_Ratio_Uncertainty',
    units       => '-',
    description => 'Integrated_Volume_Depolarization_Ratio_Uncertainty'
  },
  {
    name        => 'L2_01kmCLay.Layer_Base_Altitude',
    units       => 'km',
    description => 'Base of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_01kmCLay.Layer_Base_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer base altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_01kmCLay.Layer_Centroid_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer attenuated backscatter centroid altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_01kmCLay.Layer_IAB_QA_Factor',
    units       => '-',
    description => 'Single layer analog of the column integrated attenuated backscatter cumulative probability.'
  },
  {
    name        => 'L2_01kmCLay.Layer_Top_Altitude',
    units       => 'km',
    description => 'Top of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_01kmCLay.Layer_Top_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_01kmCLay.Layer_Top_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_01kmCLay.Low_Energy_Mitigation_Column_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each column; least significant bit = bit 1 (bit; interpretation): bit 1; Column contains LEM affected data (data has been rejected or low energy shots that LEM accepts, bit 2; Column belongs to a 5km resolution frame that has been rejected due to too many unusable profiles, bit 3; Column belongs to a 5 km resolution frame that has been rejected due to too many rejected subregions in altitude region 3, bit 4; Column belongs to a 5km resolution frame that has been rejected due to too many rejected subregions in altitude region 4, bit 5; Feature detection at 20km resolution not performed in this column due to insufficient fraction of usable 5km resolution frames, bit 6; Feature detection at 80 km resolution not performed in this column due to insufficient fraction of usable 5 km resolution frames, bits 7-16; Unused.'
  },
  {
    name        => 'L2_01kmCLay.Low_Energy_Mitigation_Feature_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each layer detected (flag - interpretation): 0 - No low energy shots affecting the layer, 1 - Some portion of the layer contains low energy shots and/or data rejected by the LEM algorithm'
  },
  {
    name        => 'L2_01kmCLay.Midlayer_Pressure',
    units       => 'hPa',
    description => 'Pressure at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_01kmCLay.Midlayer_Temperature',
    units       => 'degC',
    description => 'Temperature at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_01kmCLay.Minimum_Laser_Energy_532',
    units       => 'J',
    description => 'Minimum laser energy at 532 nm measured within each 80 km along-track data segment (80km = 240 single-shot laser pulses). The 80km distance matches the largest horizontal extent considered in CALIOP standard level 2 data analyses.'
  },
  {
    name        => 'L2_01kmCLay.Number_Layers_Found',
    units       => '-',
    description => 'The number of cloud layers found in a column. The merged 5km layer product reports layers detected at multiple averaging resolutions (5km, 20km, and 80km). Because the reporting resolution (5-km) is not always identical to the detection resolution, layers may appear to overlap in the vertical dimension.'
  },
  {
    name        => 'L2_01kmCLay.ODCOD_Effective_Optical_Depth_532',
    units       => '-',
    description => 'Column estimate of apparent (multiple scattering effects not accounted for) optical depth at 532nm wavelength for all particulate (clouds and aerosols) within the illuminated volume from just above the lidar detected ocean surface to the CALIPSO top of atmosphere.'
  },
  {
    name        => 'L2_01kmCLay.ODCOD_Effective_Optical_Depth_532_Uncertainty',
    units       => '-',
    description => 'An analytical estimate of the absolute standard deviation uncertainty in the ODCOD Effective Optical Depth estimate at 532nm.'
  },
  {
    name        => 'L2_01kmCLay.ODCOD_QC_Flag_532',
    units       => '-',
    description => 'Ocean Derived Column Optical Depth (ODCOD) quality control flags; least significant bit = bit 1 (bit; interpretation): bit 1; The fundamental time delay was shifted from the first point, bit 2; Lidar surface data too many data points for a valid surface spike, bit 3; Range bin above recorded lidar surface data added, bit 4; Range bin below recorded lidar surface data added, bit 5; The surface point ratio indicated that the first point of the surface spike is missing, bit 6; To find the time delay, the surface spike index range had to be shifted down, bit 7; The retrieval inputs fall outside of confidence ranges indicating a confident retrieval, bits 8-10; Unused, none, bit 11; No Lidar surface was found, bit 12; The surface is not ocean, bit 13; The surface is ice, bit 14; The wind speed is outside the valid range; 0.025 to 43 m/s, bit 15; The fundamental time delay found was too low and there were not enough points to find a new one, bit 16; Too few data points to solve, bit 17; The low pass filter impulse response model curve area was too large, bit 18; The low pass filter impulse response model curve scale factor could no t be found, bit 19; Lidar surface saturated, bit 20; Lidar surface has negative signal anomaly, bit 21; Lidar surface data is not valid, bit 22; Flags bad input/ancillary data, bit 23; No Lidar surface was detected via derivative method at greater than single shot resolutions, bits 24-32; Unused, none.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_01kmCLay.ODCOD_Surface_Wind_Speed_Correction',
#    units       => 'm/s',
#    description => 'An additive correction derived from AMSR analysis applied to the MERRA-2 wind speed magnitude calculated from the zonal winds reported in ODCOD Surface Wind Speeds 10m.'
#  },
#  {
#    name        => 'L2_01kmCLay.ODCOD_Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Mean zonal and meridional component of the surface wind speed 10 meters above the surface provided by MERRA-2. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v.'
#  },
  {
    name        => 'L2_01kmCLay.Off_Nadir_Angle',
    units       => 'deg',
    description => 'Angle of the lidar viewing vector; 0.3 degrees off-nadir from launch/first-light to 3.0 degrees after November 28, 2007.'
  },
  {
    name        => 'L2_01kmCLay.Opacity_Flag',
    units       => '-',
    description => 'Flag that identifies a layer in which the backscatter signal becomes completely attenuated. An opacity flag of 0 indicates a transparent layer, while a value of 1 indicates an opaque layer. A layer is considered opaque if (a) it is the lowest feature detected in a column, and (b) it is not subsequently classified as a surface return. For those layers identified as opaque, the reported base altitude must be considered as an apparent base, rather than a true base.'
  },
  {
    name        => 'L2_01kmCLay.Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'Overlying_Integrated_Attenuated_Backscatter_532'
  },
  {
    name        => 'L2_01kmCLay.Parallel_Column_Reflectance_532',
    units       => '-',
    description => 'Bi-directional column reflectance derived from the root-mean-square variation of the 532 nm parallel channel baseline measurements.'
  },
  {
    name        => 'L2_01kmCLay.Parallel_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Column reflectance uncertainty for the 532nm parallel channel.'
  },
  {
    name        => 'L2_01kmCLay.Perpendicular_Column_Reflectance_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_532'
  },
  {
    name        => 'L2_01kmCLay.Perpendicular_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_Uncertainty_532'
  },
  {
    name        => 'L2_01kmCLay.Profile_ID',
    units       => '-',
    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
  },
  {
    name        => 'L2_01kmCLay.Profile_Time',
    units       => 'seconds_since_1993-01-01',
    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
  },
  {
    name        => 'L2_01kmCLay.Scattering_Angle',
    units       => 'deg',
    description => 'Angle between the lidar viewing angle and the line of sight of the sun.'
  },
  {
    name        => 'L2_01kmCLay.Scene_Flag',
    units       => '-',
    description => 'Bit mapped flag where the following bits indicate the atmospheric feature type was detected at some resolution in the colum n; least significant bit = bit 1 (bit; interpretation): bit 1; tropospheric aerosol subtype clean marine, bit 2; tropospheric aerosol subtype dust, bit 3; tropospheric aerosol subtype polluted continental/smoke, bit 4; tropospheric aerosol subtype clean continental, bit 5; tropospheric aerosol subtype polluted dust, bit 6; tropospheric aerosol subtype elevated smoke, bit 7; tropospheric aerosol subtype dusty marine, bit 8; stratospheric aerosol subtype polar stratospheric aerosol, bit 9; stratospheric aerosol subtype volcanic ash, bit 10; stratospheric aerosol subtype sulfate, bit 11; stratospheric aerosol subtype elevated smoke bit 12; stratospheric aerosol unclassified. bit 13; On if column has ROI clouds at any resolution. bit 14; On if column has HOI clouds at any resolution. bit 15; On if column has water clouds at any resolution. bit 16; On if column has unknown phase clouds at any resolution.'
  },
  {
    name        => 'L2_01kmCLay.Snow_Ice_Surface_Type',
    units       => '-',
    description => 'Snow and ice coverage for the surface at the laser footprint provided by the U.S. Airforce Weather Agency (AFWA). Values between 0 and 100 represent the percentage of ice, while 101 = permanent ice, 103 = snow, and 255 = mixed pixels at coastline.'
  },
  {
    name        => 'L2_01kmCLay.Solar_Azimuth_Angle',
    units       => 'degrees',
    description => 'Azimuth angle from north of the line of sight to the sun.'
  },
  {
    name        => 'L2_01kmCLay.Solar_Zenith_Angle',
    units       => 'degrees',
    description => 'Angle between the zenith at the lidar footprint on the surface and the line of sight to the sun.'
  },
  {
    name        => 'L2_01kmCLay.Spacecraft_Position_X',
    units       => 'km',
    description => 'X-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_01kmCLay.Spacecraft_Position_X',
    units       => 'km',
    description => 'X-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_01kmCLay.Spacecraft_Position_Y',
    units       => 'km',
    description => 'Y-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_01kmCLay.Spacecraft_Position_Z',
    units       => 'km',
    description => 'Z-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_01kmCLay.Surface_1064_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 1064 to Surface Base 1064, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_01kmCLay.Surface_1064_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 1064 to Surface Base Altitude 1064, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_01kmCLay.Surface_532_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_01kmCLay.Surface_532_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Base_Altitude_1064',
    units       => 'km',
    description => 'Base altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Base_Altitude_532',
    units       => 'km',
    description => 'Base altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Detection_Flags_1064',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Detection_Flags_532',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Detections_333m_1064',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 1064nm channel.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Detections_333m_532',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 532nm channel.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 1064nm of the surface return from Surface Top Altitude 1064 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 532nm of the surface return from Surface Top Altitude 532 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Overlying_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 1064 from 40km to one range bin above Surface Top Altitude 1064.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 532 from 40km to one range bin above Surface Top Altitude 532.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Peak_Signal_1064',
    units       => '/km/sr',
    description => 'Maximum 1064nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Peak_Signal_532',
    units       => '/km/sr',
    description => 'Maximum 532nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Scaled_RMS_Background_1064',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 1064nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Scaled_RMS_Background_532',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 532nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Top_Altitude_1064',
    units       => 'km',
    description => 'Top altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_01kmCLay.Surface_Top_Altitude_532',
    units       => 'km',
    description => 'Top altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_01kmCLay.Tropopause_Height',
    units       => 'km',
    description => 'Tropopause height. Provided by MERRA-2.'
  },
  {
    name        => 'L2_01kmCLay.Tropopause_Temperature',
    units       => 'degC',
    description => 'Tropopause temperature. Provided by MERRA-2.'
  },



  # L2_05kmALay ---------------------------------------------------------------

  {
    name        => 'L2_05kmALay.Column_Feature_Fraction',
    units       => '-',
    description => 'The fraction of the 5-km horizontally averaged profile, between 30-km and the DEM surface elevation, which has been identified as containing a feature (i.e., either a cloud, an aerosol, or a stratospheric layer).'
  },
  {
    name        => 'L2_05kmALay.Column_IAB_Cumulative_Probability',
    units       => '-',
    description => 'The cumulative probability of measuring a total column integrated attenuated backscatter value equal to the value computed for the current profile.'
  },
  {
    name        => 'L2_05kmALay.Column_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'The integral with respect to altitude of the 532nm total attenuated backscatter coefficients. The limits of integration are from the onset of the backscatter signal at ~40km, down to the range bin immediately prior to the surface elevation specified by the digital elevation map (DEM).'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Cloud_532',
    units       => '-',
    description => 'Optical depth of all clouds detected with in a 5km averaged profile, obtained by integrating the 532nm cloud extinction profile reported in the CALIPSO 5km Cloud Profile Products.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Cloud_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Cloud 532 parameter.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Stratospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5-km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 1064nm in the CALIPSO 5 km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Stratospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5-km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 532nm in the CALIPSO 5 km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 1064 parameter.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 532 parameter.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Tropospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols detected within a 5-km averaged profile, obtained by integrating the 1064nm aerosol extinction profile reported in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Tropospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols detected within a 5-km averaged profile, obtained by integrating the 532nm aerosol extinction profile reported in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric Aerosol 1064 parameter.'
  },
  {
    name        => 'L2_05kmALay.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric Aerosol 532 parameter.'
  },
  {
    name        => 'L2_05kmALay.Column_Particulate_Optical_Depth_Above_Opaque_Water_Cloud_532',
    units       => '-',
    description => '532nm particulate optical depths in the vertical column above an opaque water cloud, derived according to the method detailed in Hu et al., 2007 (DOI: 10.1109/LGRS.2007.901085).'
  },
  {
    name        => 'L2_05kmALay.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L2_05kmALay.Extinction_QC_Flag_1064',
    units       => '-',
    description => 'The 1064nm extinction Quality Control flags are bit-mapped 16-bit integers, reported for each layer for which an extinction retrieval was attempted; least significant bit = bit 1 (bit; interpretation): bit 1; Retrieval, 0 = Unconstrained retrieval, initial lidar ratio unchanged during solution process, 1 = Constrained retrieval, bit 2; Initial lidar ratio reduced to achieve successful full-column solutions, bit 3; Suspicious retrieval due to layer or overlying integrated attenuated backscatter being too high or excessive lidar ratio reductions, bit 4; Initial lidar ratio has been reduced and has converged, but backscatter uncertainty solution does not exist, bit 5; Layer being analyzed has been identified by the feature finder as being totally attenuating (i.e., opaque), bit 6; Estimated optical depth error exceeds the maximum allowable value, bit 7; Negative signal anomaly detected, bit 8; Retrieval terminated at maximum iterations for a constrained retrieval, bit 9; No solution possible within allowable lidar ratio bounds, bit 10; Two-way particulate transmittance has converged but constrained retrieval still not achieved, bit 11; Backscatter coefficients not converging and maximum lidar ratio correction iterations reached, bit 12; Uncertainties not converging and maximum lidar ratio correction iterations achieved, bit 13; Lidar ratio converged but retrieval still not converging, bit 14; Range bin contains two different feature types that may have different lidar ratios, bit 15; Complex retrieval failure, bit 16; Fill value or no solution attempted.'
  },
  {
    name        => 'L2_05kmALay.Extinction_QC_Flag_532',
    units       => '-',
    description => 'The 532nm extinction Quality Control flags are bit-mapped 16-bit integers, reported for each layer for which an extinction retrieval was attempted; least significant bit = bit 1 (bit; interpretation): bit 1; Retrieval, 0 = Unconstrained retrieval, initial lidar ratio unchanged during solution process, 1 = Constrained retrieval, bit 2; Initial lidar ratio reduced to achieve successful full-column solutions, bit 3; Suspicious retrieval due to layer or overlying integrated attenuated backscatter being too high or excessive lidar ratio reductions, bit 4; Initial lidar ratio has been reduced and has converged, but backscatter uncertainty solution does not exist, bit 5; Layer being analyzed has been identified by the feature finder as being totally attenuating (i.e., opaque), bit 6; Estimated optical depth error exceeds the maximum allowable value, bit 7; Negative signal anomaly detected, bit 8; Retrieval terminated at maximum iterations for a constrained retrieval, bit 9; No solution possible within allowable lidar ratio bounds, bit 10; Two-way particulate transmittance has converged but constrained retrieval still not achieved, bit 11; Backscatter coefficients not converging and maximum lidar ratio correction iterations reached, bit 12; Uncertainties not converging and maximum lidar ratio correction iterations achieved, bit 13; Lidar ratio converged but retrieval still not converging, bit 14; Range bin contains two different feature types that may have different lidar ratios, bit 15; Complex retrieval failure, bit 16; Fill value or no solution attempted.'
  },
  {
    name        => 'L2_05kmALay.Feature_Classification_Flags',
    units       => '-',
    description => 'For each layer detected in the CALIPSO backscatter data a feature classification flags are derived. The complete set of flags are stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): 011 bits 1-3; Feature Type, 0110 = invalid (bad or missing data), 0111 = clear air, 0112 = cloud, 0113 = tropospheric aerosol, 0114 = stratospheric aerosol, 0115 = surface, 0116 = subsurface, 0117 = no signal, bits 4-5; Feature Type QA, 011 = 0 = none, 0111 = low, 0112 = medium, 0113 = high 011 bits 6-7; Ice/Water Phase (cloud), 0110 = unknown/not determined, 0111 = ice, 0112 = water, 0113 = oriented ice crystals, bits 8-9; Ice/Water Phase QA (cloud), 0110 = none, 0111 = low, 0112 = medium, 0113 = high, 011bits 10-12; Feature Subtype (clear-air), 0110 = n/a, 0111 = not searched for features at 80km resolution, 0112 = not searched for features at 20 and 8km resolution, bits 10-12; Feature Subtype (Cloud), 0110 = low overcast, transparent, 0111 = low overcast, opaque, 0112 = transition stratocumulus, 0113 = low, broken cloud.'
  },
  {
    name        => 'L2_05kmALay.FeatureFinderQC',
    units       => '-',
    description => 'To generate data at a nominal 5-km horizontal resolution requires averaging 15 consecutive single-shot laser pulses. For each 5km average a set of feature finder QC flags are reported, which are a set of 15 Boolean values that report whether or not a feature (cloud, aerosol, or surface echo) was detected in each of the 15 laser pulses. The flags are implemented as a 16-bit integer, with the most significant bit unused.'
  },
  {
    name        => 'L2_05kmALay.Feature_Optical_Depth_1064',
    units       => '-',
    description => 'Reports estimates of 1064nm layer optical depth.'
  },
  {
    name        => 'L2_05kmALay.Feature_Optical_Depth_532',
    units       => '-',
    description => 'Reports estimates of 532nm layer optical depth.'
  },
  {
    name        => 'L2_05kmALay.Feature_Optical_Depth_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty, reported as an absolute error, in the 1064nm layer optical  depth.'
  },
  {
    name        => 'L2_05kmALay.Feature_Optical_Depth_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty, reported as an absolute error, in the 532nm layer optical  depth.'
  },
  {
    name        => 'L2_05kmALay.Final_1064_Lidar_Ratio',
    units       => 'sr',
    description => 'Reports the lidar ratio in use at the conclusion of the extinction processing for each layer. The final lidar ratio may be (1) the initial lidar ratio supplied by the Scene Classification Algorithms, (2) the result of modifications to this initial lidar ratio to avoid a non-physical solution, or (3) a lidar ratio determined from a measured layer transmittance.'
  },
  {
    name        => 'L2_05kmALay.Final_1064_Lidar_Ratio_Uncertainty',
    units       => 'sr',
    description => 'Reports the lidar ratio uncertainty at the conclusion of the extinction processing for each layer.'
  },
  {
    name        => 'L2_05kmALay.Final_532_Lidar_Ratio',
    units       => 'sr',
    description => 'Reports the lidar ratio in use at the conclusion of the extinction processing for each layer. The final lidar ratio may be (1) the initial lidar ratio supplied by the Scene Classification Algorithms, (2) the result of modifications to this initial lidar ratio to avoid a non-physical solution, or (3) a lidar ratio determined from a measured layer transmittance.'
  },
  {
    name        => 'L2_05kmALay.Final_532_Lidar_Ratio_Uncertainty',
    units       => 'sr',
    description => 'Reports the lidar ratio uncertainty at the conclusion of the extinction processing for each layer.'
  },
  {
    name        => 'L2_05kmALay.High_Resolution_Layers_Cleared',
    units       => '-',
    description => 'Bit-mapped integer identifying which single-shot resolution columns in a 5km average contain clouds detected and removed by the boundary layer cloud-clearing algorithm. The flags are implemented as a 16-bit integer, with the most significant bit unused.'
  },
  {
    name        => 'L2_05kmALay.Horizontal_Averaging',
    units       => 'km',
    description => 'The amount of horizontal averaging required for a feature to be detected.'
  },
  {
    name        => 'L2_05kmALay.IGBP_Surface_Type',
    units       => '-',
    description => 'IGBP surface type at the laser footprint provided by the CERES/SARB surface map: 1 = Evergreen-Needleleaf-Forest, 2 = Evergreen-Broadleaf-Forest, 3 = Deciduous-Needleleaf-Forest, 4 = Deciduous-Broadleaf-Forest, 5 = Mixed-Forest, 6 = Closed-Shrubland, 7 = Open-Shrubland/Desert, 8 = Woody-Savanna, 9 = Savanna, 10 = Grassland, 11 = Wetland, 12 = Cropland, 13 = Urban, 14 = Crop-Mosaic, 15 = Permanent-Snow, 16 = Barren/Desert, 17 = Water, 18 = Tundra.'
  },
  {
    name        => 'L2_05kmALay.Initial_1064_Lidar_Ratio',
    units       => 'sr',
    description => 'Retrieving optical depth and profiles of extinction and backscatter coefficients from the CALIOP measurements requires an estimate of the particulate extinction-to-backscatter ratio, which in the lidar community is commonly known as the lidar ratio. These initial estimates are selected based on the type and subtype of the layer being analyzed.'
  },
  {
    name        => 'L2_05kmALay.Final_1064_Lidar_Ratio',
    units       => 'sr',
    description => 'Reports the lidar ratio in use at the conclusion of the extinction processing for each layer. The final lidar ratio may be (1) the initial lidar ratio supplied by the Scene Classification Algorithms, (2) the result of modifications to this initial lidar ratio to avoid a non-physical solution, or (3) a lidar ratio determined from a measured layer transmittance.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'The 532 nm integrated attenuated backscatter (IAB) computed for each layer after correcting the backscatter signal for 532nm attenuation by any overlying layers.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Attenuated_Backscatter_1064',
    units       => '1/sr',
    description => 'The 1064nm integrated attenuated backscatter (IAB) computed for each layer.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Attenuated_Backscatter_Uncertainty_1064',
    units       => '1/sr',
    description => 'The uncertainties reported for the 1064nm integrated attenuated backscatters, which provide an estimate of the random error in the backscatter signal.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Attenuated_Backscatter_Uncertainty_532',
    units       => '1/sr',
    description => 'Estimated uncertainties the 532 nm integrated attenuated backscatter values. These uncertainties quantify the random error in the integral of the attenuation-corrected 532nm backscatter signal.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Attenuated_Total_Color_Ratio',
    units       => '-',
    description => 'Integrated attenuated total color ratio computed for each layer.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Attenuated_Total_Color_Ratio_Uncertainty',
    units       => '-',
    description => 'Integrated_Attenuated_Total_Color_Ratio_Uncertainty'
  },
  {
    name        => 'L2_05kmALay.Integrated_Particulate_Color_Ratio',
    units       => '-',
    description => 'The integrated particulate color ratio derived from ratio of the particulate backscatter coefficients at 532nm and 1064nm, each summed over the vertical extent of the layer.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Particulate_Color_Ratio_Uncertainty',
    units       => '-',
    description => 'Uncertainty associated with the estimate of the layer-integrated particulate backscatter color ratio, calculated by first integrating the uncertainties for the particulate backscatter coefficients measured at each wavelength, and then applying standard techniques for error propagation in ratioed quantities.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Particulate_Depolarization_Ratio',
    units       => '-',
    description => 'Layer-integrated particulate depolarization ratio. The particulate depolarization ratio represents the contribution to the volume depolarization ratio that is due only to the cloud and/or aerosol particles within the layer.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Particulate_Depolarization_Ratio_Uncertainty',
    units       => '-',
    description => 'Uncertainty associated with the estimate of the layer-integrated particulate depolarization ratio, calculated by integrating the uncertainties for the particulate backscatter coefficients measured in the 532nm parallel and perpendicular channels, and then applying standard techniques for error propagation in ratioed quantities.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Volume_Depolarization_Ratio',
    units       => '-',
    description => 'The layer integrated 532 nm volume depolarization ratio (IVDR) computed for each layer. The quality of the estimate for IVDR is determined by the accuracy of the top and base identification, the reliability of the polarization gain ratio calibration, and by the signal-to-noise ratio (SNR) of the backscatter data within the layer.'
  },
  {
    name        => 'L2_05kmALay.Integrated_Volume_Depolarization_Ratio_Uncertainty',
    units       => '-',
    description => 'The uncertainties reported for the 532nm layer-integrated volume depolarization ratios provide an estimate of the total random error in the combined backscatter signals.'
  },
  {
    name        => 'L2_05kmALay.Layer_Base_Altitude',
    units       => 'km',
    description => 'Base of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_05kmALay.Layer_Base_Extended',
    units       => '-',
    description => 'A non-zero value indicates that the base of the layer has been extended from the initial altitude assigned by the layer detection algorithm to a new, lower altitude lying three range bins (90 m) above the Earth surface. Non-zero values represent the layer feature classification flags prior to the base being extended.'
  },
  {
    name        => 'L2_05kmALay.Layer_Base_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer base altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmALay.Layer_Base_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer base altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmALay.Layer_Effective_1064_Multiple_Scattering_Factor',
    units       => '-',
    description => 'The layer effective 1064 nm multiple scattering factor corresponds to the limit of single scattering only, with smaller values indicating increasing contributions to the backscatter signal from multiple scattering.'
  },
  {
    name        => 'L2_05kmALay.Layer_Effective_532_Multiple_Scattering_Factor',
    units       => '-',
    description => 'The layer effective 532 nm multiple scattering factor corresponds to the limit of single scattering only, with smaller values indicating increasing contributions to the backscatter signal from multiple scattering.'
  },
  {
    name        => 'L2_01kmALay.Layer_IAB_QA_Factor',
    units       => '-',
    description => 'Single layer analog of the column integrated attenuated backscatter cumulative probability.'
  },
  {
    name        => 'L2_01kmALay.Layer_Top_Altitude',
    units       => 'km',
    description => 'Top of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_05kmALay.Layer_Top_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmALay.Layer_Top_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmALay.Lidar_Ratio_1064_Selection_Method',
    units       => '-',
    description => 'Specifies the internal procedure used to select the initial lidar ratio for each layer:, - not determined, - Transmittance Method, constrained retrieval, - Opaque Layer Method, based on opaque clouds or aerosols, - Cirrus Temperature Method, unconstrained ice clouds, - Unknown Cloud Phase, unconstrained unknown phase clouds with centroid temperature between -37 and 0 degrees C, - Fixed Layer Method, unconstrained water clouds and aerosols, - Solution at 532nm Method, unconstrained 1064 nm cloud only when 532 nm extinction solution is available, - MAC Table Interpolation, unconstrained aerosol using MAC interpolated lookup table.'
  },
  {
    name        => 'L2_05kmALay.Lidar_Ratio_532_Selection_Method',
    units       => '-',
    description => 'Specifies the internal procedure used to select the initial lidar ratio for each layer: - not determined, - Transmittance Method, constrained retrieval, - Opaque Layer Method, based on opaque clouds or aerosols, - Cirrus Temperature Method, unconstrained ice clouds, - Unknown Cloud Phase, unconstrained unknown phase clouds with centroid temperature between -37 and 0 degrees C, - Fixed Layer Method, unconstrained water clouds and aerosols, - Solution at 532nm Method, unconstrained 1064 nm cloud only when 532 nm extinction solution is available\012\0117 - MAC Table Interpolation, unconstrained aerosol using MAC interpolated lookup table.'
  },
  {
    name        => 'L2_05kmALay.Low_Energy_Mitigation_Column_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each column; least significant bit = bit 1 (bit; interpretation): bit 1; Column contains LEM affected data (data has been rejected or low energy shots that LEM accepts, bit 2; Column belongs to a 5km resolution frame that has been rejected due to too many unusable profiles, bit 3; Column belongs to a 5 km resolution frame that has been rejected due to too many rejected subregions in altitude region 3, bit 4; Column belongs to a 5km resolution frame that has been rejected due to too many rejected subregions in altitude region 4, bit 5; Feature detection at 20km resolution not performed in this column due to insufficient fraction of usable 5km resolution frames, bit 6; Feature detection at 80 km resolution not performed in this column due to insufficient fraction of usable 5 km resolution frames, bits 7-16; Unused.'
  },
  {
    name        => 'L2_05kmALay.Low_Energy_Mitigation_Feature_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each layer detected (flag - interpretation): 0 - No low energy shots affecting the layer, 1 - Some portion of the layer contains low energy shots and/or data rejected by the LEM algorithm'
  },
  {
    name        => 'L2_05kmALay.Measured_Two_Way_Transmittance_532',
    units       => '-',
    description => 'Provides the measured value of the layer two-way transmittance, derived by computing the ratio of the mean attenuated scattering ratios in the "clear air" regions immediately below and above the layer, for isolated transparent layers.'
  },
  {
    name        => 'L2_05kmALay.Measured_Two_Way_Transmittance_Uncertainty_532',
    units       => '-',
    description => 'The relative error in the two-way transmittance measurement, calculated using standard techniques for error propagation in ratioed quantities.'
  },
  {
    name        => 'L2_05kmALay.Midlayer_Pressure',
    units       => 'hPa',
    description => 'Pressure at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmALay.Midlayer_Temperature',
    units       => 'degC',
    description => 'Temperature at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmALay.Minimum_Laser_Energy_532',
    units       => 'J',
    description => 'Minimum laser energy at 532 nm measured within each 80 km along-track data segment (80km = 240 single-shot laser pulses). The 80km distance matches the largest horizontal extent considered in CALIOP standard level 2 data analyses.'
  },
# unmatched dimensions
#  {
#    name        => 'L2_05kmALay.Normalization_Constant_Uncertainty',
#    units       => '-',
#    description => 'Uncertainty, reported as a relative error, in the 532nm and 1064nm calibration constants due solely to random error in the backscatter measurements in the calibration region.'
#  },
  {
    name        => 'L2_05kmALay.Number_Layers_Found',
    units       => '-',
    description => 'The number of cloud layers found in a column. The merged 5km layer product reports layers detected at multiple averaging resolutions (5km, 20km, and 80km). Because the reporting resolution (5-km) is not always identical to the detection resolution, layers may appear to overlap in the vertical dimension.'
  },
  {
    name        => 'L2_05kmALay.ODCOD_Effective_Optical_Depth_532',
    units       => '-',
    description => 'Column estimate of apparent (multiple scattering effects not accounted for) optical depth at 532nm wavelength for all particulate (clouds and aerosols) within the illuminated volume from just above the lidar detected ocean surface to the CALIPSO top of atmosphere.'
  },
  {
    name        => 'L2_05kmALay.ODCOD_Effective_Optical_Depth_532_Uncertainty',
    units       => '-',
    description => 'An analytical estimate of the absolute standard deviation uncertainty in the ODCOD Effective Optical Depth estimate at 532nm.'
  },
  {
    name        => 'L2_05kmALay.ODCOD_QC_Flag_532',
    units       => '-',
    description => 'Ocean Derived Column Optical Depth (ODCOD) quality control flags; least significant bit = bit 1 (bit; interpretation): bit 1; The fundamental time delay was shifted from the first point, bit 2; Lidar surface data too many data points for a valid surface spike, bit 3; Range bin above recorded lidar surface data added, bit 4; Range bin below recorded lidar surface data added, bit 5; The surface point ratio indicated that the first point of the surface spike is missing, bit 6; To find the time delay, the surface spike index range had to be shifted down, bit 7; The retrieval inputs fall outside of confidence ranges indicating a confident retrieval, bits 8-10; Unused, none, bit 11; No Lidar surface was found, bit 12; The surface is not ocean, bit 13; The surface is ice, bit 14; The wind speed is outside the valid range; 0.025 to 43 m/s, bit 15; The fundamental time delay found was too low and there were not enough points to find a new one, bit 16; Too few data points to solve, bit 17; The low pass filter impulse response model curve area was too large, bit 18; The low pass filter impulse response model curve scale factor could no t be found, bit 19; Lidar surface saturated, bit 20; Lidar surface has negative signal anomaly, bit 21; Lidar surface data is not valid, bit 22; Flags bad input/ancillary data, bit 23; No Lidar surface was detected via derivative method at greater than single shot resolutions, bits 24-32; Unused, none.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_05kmALay.ODCOD_Surface_Wind_Speed_Correction',
#    units       => 'm/s',
#    description => 'An additive correction derived from AMSR analysis applied to the MERRA-2 wind speed magnitude calculated from the zonal winds reported in ODCOD Surface Wind Speeds 10m.'
#  },
#  {
#    name        => 'L2_05kmALay.ODCOD_Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Mean zonal and meridional component of the surface wind speed 10 meters above the surface provided by MERRA-2. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v.'
#  },
  {
    name        => 'L2_05kmALay.Off_Nadir_Angle',
    units       => 'deg',
    description => 'Angle of the lidar viewing vector; 0.3 degrees off-nadir from launch/first-light to 3.0 degrees after November 28, 2007.'
  },
  {
    name        => 'L2_05kmALay.Opacity_Flag',
    units       => '-',
    description => 'Flag that identifies a layer in which the backscatter signal becomes completely attenuated. An opacity flag of 0 indicates a transparent layer, while a value of 1 indicates an opaque layer. A layer is considered opaque if (a) it is the lowest feature detected in a column, and (b) it is not subsequently classified as a surface return. For those layers identified as opaque, the reported base altitude must be considered as an apparent base, rather than a true base.'
  },
  {
    name        => 'L2_05kmALay.Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'Overlying_Integrated_Attenuated_Backscatter_532'
  },
  {
    name        => 'L2_05kmALay.Parallel_Column_Reflectance_532',
    units       => '-',
    description => 'Bi-directional column reflectance derived from the root-mean-square variation of the 532 nm parallel channel baseline measurements.'
  },
  {
    name        => 'L2_05kmALay.Parallel_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Column reflectance uncertainty for the 532nm parallel channel.'
  },
  {
    name        => 'L2_05kmALay.Perpendicular_Column_Reflectance_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_532'
  },
  {
    name        => 'L2_05kmALay.Perpendicular_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_Uncertainty_532'
  },
# unmatched dimensions
#  {
#    name        => 'L2_05kmALay.Profile_ID',
#    units       => '-',
#    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
#  },
#  {
#    name        => 'L2_05kmALay.Profile_Time',
#    units       => 'seconds_since_1993-01-01',
#    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
#  },
  {
    name        => 'L2_05kmALay.Relative_Humidity',
    units       => '%',
    description => 'Relative humidity, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmALay.Scattering_Angle',
    units       => 'deg',
    description => 'Angle between the lidar viewing angle and the line of sight of the sun.'
  },
  {
    name        => 'L2_05kmALay.Scene_Flag',
    units       => '-',
    description => 'Bit mapped flag where the following bits indicate the atmospheric feature type was detected at some resolution in the colum n; least significant bit = bit 1 (bit; interpretation): bit 1; tropospheric aerosol subtype clean marine, bit 2; tropospheric aerosol subtype dust, bit 3; tropospheric aerosol subtype polluted continental/smoke, bit 4; tropospheric aerosol subtype clean continental, bit 5; tropospheric aerosol subtype polluted dust, bit 6; tropospheric aerosol subtype elevated smoke, bit 7; tropospheric aerosol subtype dusty marine, bit 8; stratospheric aerosol subtype polar stratospheric aerosol, bit 9; stratospheric aerosol subtype volcanic ash, bit 10; stratospheric aerosol subtype sulfate, bit 11; stratospheric aerosol subtype elevated smoke bit 12; stratospheric aerosol unclassified. bit 13; On if column has ROI clouds at any resolution. bit 14; On if column has HOI clouds at any resolution. bit 15; On if column has water clouds at any resolution. bit 16; On if column has unknown phase clouds at any resolution.'
  },
  {
    name        => 'L2_05kmALay.Single_Shot_Cloud_Cleared_Fraction',
    units       => '-',
    description => 'Layers detected in the planetary boundary layer (PBL) are subjected to an additional cloud clearing procedure to separate small-scale boundary layer clouds from any surrounding PBL aerosols at the highest possible spatial resolution (i.e., single shot data). The single shot cloud cleared fraction reports the fraction of the nominal layer area (i.e., horizontal averaging distance times layer height) that was removed by the cloud clearing process.'
  },
  {
    name        => 'L2_05kmALay.Smoke_Layer_Base_Extended_QC_Flag',
    units       => '-',
    description => 'Flag used to identify aerosol layers which are candidates for smoke layer base extension, and, if selected, if the extension was successful; least significant bit = bit 1 (bit; interpretation): 1; Candidate for smoke layer extension, 2; Success, smoke layer extended, 3; Failed, smoke layer base already in contact with a cloud layer, 4; Failed, no cloud layer is found below the smoke layer, 5; Failed, gap distance between smoke ad cloud layer exceeds limit, 6; Failed, 1064nm scattering ratio test is below threshold, 7; Failed, cloud color ratio test is below threshold.'

  },
  {
    name        => 'L2_05kmALay.Snow_Ice_Surface_Type',
    units       => '-',
    description => 'Snow and ice coverage for the surface at the laser footprint provided by the U.S. Airforce Weather Agency (AFWA). Values between 0 and 100 represent the percentage of ice, while 101 = permanent ice, 103 = snow, and 255 = mixed pixels at coastline.'
  },
  {
    name        => 'L2_05kmALay.Solar_Azimuth_Angle',
    units       => 'degrees',
    description => 'Azimuth angle from north of the line of sight to the sun.'
  },
  {
    name        => 'L2_05kmALay.Solar_Zenith_Angle',
    units       => 'degrees',
    description => 'Angle between the zenith at the lidar footprint on the surface and the line of sight to the sun.'
  },
  {
    name        => 'L2_05kmALay.Spacecraft_Position_X',
    units       => 'km',
    description => 'X-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_05kmALay.Spacecraft_Position_Y',
    units       => 'km',
    description => 'Y-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_05kmALay.Spacecraft_Position_Z',
    units       => 'km',
    description => 'Z-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_05kmALay.Surface_1064_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 1064 to Surface Base 1064, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmALay.Surface_1064_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 1064 to Surface Base Altitude 1064, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmALay.Surface_532_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmALay.Surface_532_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmALay.Surface_Base_Altitude_1064',
    units       => 'km',
    description => 'Base altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmALay.Surface_Base_Altitude_532',
    units       => 'km',
    description => 'Base altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmALay.Surface_Detection_Flags_1064',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmALay.Surface_Detection_Flags_532',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmALay.Surface_Detections_333m_1064',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 1064nm channel.'
  },
  {
    name        => 'L2_05kmALay.Surface_Detections_333m_532',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 532nm channel.'
  },
  {
    name        => 'L2_05kmALay.Surface_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 1064nm of the surface return from Surface Top Altitude 1064 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmALay.Surface_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 532nm of the surface return from Surface Top Altitude 532 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmALay.Surface_Overlying_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 1064 from 40km to one range bin above Surface Top Altitude 1064.'
  },
  {
    name        => 'L2_05kmALay.Surface_Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 532 from 40km to one range bin above Surface Top Altitude 532.'
  },
  {
    name        => 'L2_05kmALay.Surface_Peak_Signal_1064',
    units       => '/km/sr',
    description => 'Maximum 1064nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmALay.Surface_Peak_Signal_532',
    units       => '/km/sr',
    description => 'Maximum 532nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmALay.Surface_Scaled_RMS_Background_1064',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 1064nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmALay.Surface_Scaled_RMS_Background_532',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 532nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmALay.Surface_Top_Altitude_1064',
    units       => 'km',
    description => 'Top altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmALay.Surface_Top_Altitude_532',
    units       => 'km',
    description => 'Top altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmALay.Tropopause_Height',
    units       => 'km',
    description => 'Tropopause height. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmALay.Tropopause_Temperature',
    units       => 'degC',
    description => 'Tropopause temperature. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmALay.Unique_Layer_ID',
    units       => '-',
    description => 'Unique number for each layer detected at
                         5-km or coarser resolution.'
  },



  # L2_05kmAPro ---------------------------------------------------------------

  {
    name        => 'L2_05kmAPro.Aerosol_Layer_Fraction',
    units       => '-',
    description => 'Aerosol_Layer_Fraction'
  },
  {
    name        => 'L2_05kmAPro.Aerosol_Multiple_Scattering_Profile_1064',
    units       => '-',
    description => 'Aerosol_Multiple_Scattering_Profile_1064'
  },
  {
    name        => 'L2_05kmAPro.Aerosol_Multiple_Scattering_Profile_532',
    units       => '-',
    description => 'Aerosol_Multiple_Scattering_Profile_532'
  },
  {
    name        => 'L2_05kmAPro.Atmospheric_Volume_Description',
    units       => '-',
    description => 'Atmospheric_Volume_Description'
  },
  {
    name        => 'L2_05kmAPro.Backscatter_Coefficient_1064',
    units       => 'sr^-1km^-1',
    description => 'Backscatter_Coefficient_1064'
  },
  {
    name        => 'L2_05kmAPro.Backscatter_Coefficient_Uncertainty_1064',
    units       => 'sr^-1km^-1',
    description => 'Backscatter_Coefficient_Uncertainty_1064'
  },
  {
    name        => 'L2_05kmAPro.Cloud_Layer_Fraction',
    units       => '-',
    description => 'Cloud_Layer_Fraction'
  },
  {
    name        => 'L2_05kmAPro.Column_Feature_Fraction',
    units       => '-',
    description => 'The fraction of the 5-km horizontally averaged profile, between 30km and the DEM surface elevation, which has been identified as containing a feature (i.e., either a cloud, an aerosol, or a stratospheric layer).'
  },
  {
    name        => 'L2_05kmAPro.Column_IAB_Cumulative_Probability',
    units       => '-',
    description => 'The cumulative probability of measuring a total column integrated attenuated backscatter value equal to the value computed for the current profile.'
  },
  {
    name        => 'L2_05kmAPro.Column_Integrated_Attenuated_Backscatter_532',
    units       => 'sr^-1',
    description => 'The integral with respect to altitude of the 532nm total attenuated backscatter coefficients. The limits of integration are from the onset of the backscatter signal at ~40 km, down to the range bin immediately prior to the surface elevation specified by the digital elevation map (DEM).'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Cloud_532',
    units       => '-',
    description => 'Column_Optical_Depth_Cloud_532'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Cloud_Uncertainty_532',
    units       => '-',
    description => 'Column_Optical_Depth_Cloud_Uncertainty_532'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Stratospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5-km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 1064nm in the CALIPSO 5 km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Stratospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5-km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 532nm in the CALIPSO 5 km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 1064 parameter.'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 532 parameter.'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Tropospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols detected within a 5-km averaged profile, obtained by integrating the 1064nm aerosol extinction profile reported in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Tropospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols detected within a 5-km averaged profile, obtained by integrating the 532nm aerosol extinction profile reported in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric Aerosol 1064 parameter.'
  },
  {
    name        => 'L2_05kmAPro.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric Aerosol 532 parameter.'
  },
  {
    name        => 'L2_05kmAPro.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L2_05kmAPro.Extinction_Coefficient_1064',
    units       => 'km^-1',
    description => 'Particulate 1064nm extinction coefficients reported for each profile range bin in which the appropriate particulates were detected.'
  },
  {
    name        => 'L2_05kmAPro.Extinction_Coefficient_532',
    units       => 'km^-1',
    description => 'Particulate 532nm extinction coefficients reported for each profile range bin in which the appropriate particulates were detected.'
  },
  {
    name        => 'L2_05kmAPro.Extinction_Coefficient_Uncertainty_1064',
    units       => 'km^-1',
    description => 'Uncertainties reported for the particulate 1064nm extinction coefficients for each profile range bin.'
  },
  {
    name        => 'L2_05kmAPro.Extinction_Coefficient_Uncertainty_532',
    units       => 'km^-1',
    description => 'Uncertainties reported for the particulate 532nm extinction coefficients for each profile range bin.'
  },
  {
    name        => 'L2_05kmAPro.Extinction_QC_Flag_1064',
    units       => '-',
    description => 'The 1064nm extinction Quality Control flags are bit-mapped 16-bit integers, reported for each layer for which an extinction retrieval was attempted; least significant bit = bit 1 (bit; interpretation): bit 1; Retrieval, 0 = Unconstrained retrieval, initial lidar ratio unchanged during solution process, 1 = Constrained retrieval, bit 2; Initial lidar ratio reduced to achieve successful full-column solutions, bit 3; Suspicious retrieval due to layer or overlying integrated attenuated backscatter being too high or excessive lidar ratio reductions, bit 4; Initial lidar ratio has been reduced and has converged, but backscatter uncertainty solution does not exist, bit 5; Layer being analyzed has been identified by the feature finder as being totally attenuating (i.e., opaque), bit 6; Estimated optical depth error exceeds the maximum allowable value, bit 7; Negative signal anomaly detected, bit 8; Retrieval terminated at maximum iterations for a constrained retrieval, bit 9; No solution possible within allowable lidar ratio bounds, bit 10; Two-way particulate transmittance has converged but constrained retrieval still not achieved, bit 11; Backscatter coefficients not converging and maximum lidar ratio correction iterations reached, bit 12; Uncertainties not converging and maximum lidar ratio correction iterations achieved, bit 13; Lidar ratio converged but retrieval still not converging, bit 14; Range bin contains two different feature types that may have different lidar ratios, bit 15; Complex retrieval failure, bit 16; Fill value or no solution attempted.'
  },
  {
    name        => 'L2_05kmAPro.Extinction_QC_Flag_532',
    units       => '-',
    description => 'The 532nm extinction Quality Control flags are bit-mapped 16-bit integers, reported for each layer for which an extinction retrieval was attempted; least significant bit = bit 1 (bit; interpretation): bit 1; Retrieval, 0 = Unconstrained retrieval, initial lidar ratio unchanged during solution process, 1 = Constrained retrieval, bit 2; Initial lidar ratio reduced to achieve successful full-column solutions, bit 3; Suspicious retrieval due to layer or overlying integrated attenuated backscatter being too high or excessive lidar ratio reductions, bit 4; Initial lidar ratio has been reduced and has converged, but backscatter uncertainty solution does not exist, bit 5; Layer being analyzed has been identified by the feature finder as being totally attenuating (i.e., opaque), bit 6; Estimated optical depth error exceeds the maximum allowable value, bit 7; Negative signal anomaly detected, bit 8; Retrieval terminated at maximum iterations for a constrained retrieval, bit 9; No solution possible within allowable lidar ratio bounds, bit 10; Two-way particulate transmittance has converged but constrained retrieval still not achieved, bit 11; Backscatter coefficients not converging and maximum lidar ratio correction iterations reached, bit 12; Uncertainties not converging and maximum lidar ratio correction iterations achieved, bit 13; Lidar ratio converged but retrieval still not converging, bit 14; Range bin contains two different feature types that may have different lidar ratios, bit 15; Complex retrieval failure, bit 16; Fill value or no solution attempted.'
  },
  {
    name        => 'L2_05kmAPro.IGBP_Surface_Type',
    units       => '-',
    description => 'IGBP surface type at the laser footprint provided by the CERES/SARB surface map: 1 = Evergreen-Needleleaf-Forest, 2 = Evergreen-Broadleaf-Forest, 3 = Deciduous-Needleleaf-Forest, 4 = Deciduous-Broadleaf-Forest, 5 = Mixed-Forest, 6 = Closed-Shrubland, 7 = Open-Shrubland/Desert, 8 = Woody-Savanna, 9 = Savanna, 10 = Grassland, 11 = Wetland, 12 = Cropland, 13 = Urban, 14 = Crop-Mosaic, 15 = Permanent-Snow, 16 = Barren/Desert, 17 = Water, 18 = Tundra.'
  },
  {
    name        => 'L2_05kmAProy.Low_Energy_Mitigation_Column_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each column; least significant bit = bit 1 (bit; interpretation): bit 1; Column contains LEM affected data (data has been rejected or low energy shots that LEM accepts, bit 2; Column belongs to a 5km resolution frame that has been rejected due to too many unusable profiles, bit 3; Column belongs to a 5 km resolution frame that has been rejected due to too many rejected subregions in altitude region 3, bit 4; Column belongs to a 5km resolution frame that has been rejected due to too many rejected subregions in altitude region 4, bit 5; Feature detection at 20km resolution not performed in this column due to insufficient fraction of usable 5km resolution frames, bit 6; Feature detection at 80 km resolution not performed in this column due to insufficient fraction of usable 5 km resolution frames, bits 7-16; Unused.'
  },
  {
    name        => 'L2_05kmAPro.Low_Energy_Mitigation_Feature_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each layer detected (flag - interpretation): 0 - No low energy shots affecting the layer, 1 - Some portion of the layer contains low energy shots and/or data rejected by the LEM algorithm'
  },
  {
    name        => 'L2_05kmAPro.Minimum_Laser_Energy_532',
    units       => 'J',
    description => 'Minimum laser energy at 532 nm measured within each 80 km along-track data segment (80km = 240 single-shot laser pulses). The 80km distance matches the largest horizontal extent considered in CALIOP standard level 2 data analyses.'
  },
  {
    name        => 'L2_05kmAPro.Molecular_Number_Density',
    units       => 'm^-3',
    description => 'Mean molecular number density reported for the midpoint of each range bin in the profile; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmAPro.ODCOD_Effective_Optical_Depth_532',
    units       => '-',
    description => 'Column estimate of apparent (multiple scattering effects not accounted for) optical depth at 532nm wavelength for all particulate (clouds and aerosols) within the illuminated volume from just above the lidar detected ocean surface to the CALIPSO top of atmosphere.'
  },
  {
    name        => 'L2_05kmAPro.ODCOD_Effective_Optical_Depth_532_Uncertainty',
    units       => '-',
    description => 'An analytical estimate of the absolute standard deviation uncertainty in the ODCOD Effective Optical Depth estimate at 532nm.'
  },
  {
    name        => 'L2_05kmAPro.ODCOD_QC_Flag_532',
    units       => '-',
    description => 'Ocean Derived Column Optical Depth (ODCOD) quality control flags; least significant bit = bit 1 (bit; interpretation): bit 1; The fundamental time delay was shifted from the first point, bit 2; Lidar surface data too many data points for a valid surface spike, bit 3; Range bin above recorded lidar surface data added, bit 4; Range bin below recorded lidar surface data added, bit 5; The surface point ratio indicated that the first point of the surface spike is missing, bit 6; To find the time delay, the surface spike index range had to be shifted down, bit 7; The retrieval inputs fall outside of confidence ranges indicating a confident retrieval, bits 8-10; Unused, none, bit 11; No Lidar surface was found, bit 12; The surface is not ocean, bit 13; The surface is ice, bit 14; The wind speed is outside the valid range; 0.025 to 43 m/s, bit 15; The fundamental time delay found was too low and there were not enough points to find a new one, bit 16; Too few data points to solve, bit 17; The low pass filter impulse response model curve area was too large, bit 18; The low pass filter impulse response model curve scale factor could no t be found, bit 19; Lidar surface saturated, bit 20; Lidar surface has negative signal anomaly, bit 21; Lidar surface data is not valid, bit 22; Flags bad input/ancillary data, bit 23; No Lidar surface was detected via derivative method at greater than single shot resolutions, bits 24-32; Unused, none.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_05kmAPro.ODCOD_Surface_Wind_Speed_Correction',
#    units       => 'm/s',
#    description => 'An additive correction derived from AMSR analysis applied to the MERRA-2 wind speed magnitude calculated from the zonal winds reported in ODCOD Surface Wind Speeds 10m.'
#  },
# {
#    name        => 'L2_05kmAPro.ODCOD_Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Mean zonal and meridional component of the surface wind speed 10 meters above the surface provided by MERRA-2. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v.'
#  },
  {
    name        => 'L2_05kmAPro.Ozone_Number_Density',
    units       => 'molecules/m3',
    description => 'Ozone number density, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmAPro.Particulate_Depolarization_Ratio_Profile_532',
    units       => '-',
    description => 'Particulate_Depolarization_Ratio_Profile_532'
  },
  {
    name        => 'L2_05kmAPro.Particulate_Depolarization_Ratio_Uncertainty_532',
    units       => '-',
    description => 'Particulate_Depolarization_Ratio_Uncertainty_532'
  },
  {
    name        => 'L2_05kmAPro.Perpendicular_Backscatter_Coefficient_532',
    units       => 'sr^-1km^-1',
    description => 'Perpendicular_Backscatter_Coefficient_532'
  },
  {
    name        => 'L2_05kmAPro.Perpendicular_Backscatter_Coefficient_Uncertainty_532',
    units       => 'sr^-1km^-1',
    description => 'Perpendicular_Backscatter_Coefficient_Uncertainty_532'
  },
  {
    name        => 'L2_05kmAPro.Pressure',
    units       => 'hPa',
    description => 'Mean pressure reported for the midpoint of each range bin in the profile; derived from the MERRA-2 data product.'
  },
# unmatched dimensions
#  {
#    name        => 'L2_05kmAPro.Profile_ID',
#    units       => '-',
#    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
#  },
#  {
#    name        => 'L2_05kmAPro.Profile_Time',
#    units       => 'seconds_since_1993-01-01',
#    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
#  },
  {
    name        => 'L2_05kmAPro.Relative_Humidity',
    units       => '%',
    description => 'Relative humidity, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmAPro.Surface_1064_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 1064 to Surface Base 1064, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmAPro.Surface_1064_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 1064 to Surface Base Altitude 1064, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmAPro.Surface_532_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmAPro.Surface_532_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Base_Altitude_1064',
    units       => 'km',
    description => 'Base altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Base_Altitude_532',
    units       => 'km',
    description => 'Base altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Detection_Flags_1064',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Detection_Flags_532',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Detections_1km_1064',
    units       => '-',
    description => 'Number of 1km resolution profiles within each 5 km resolution profile where surface was detected for the 1064nm channel.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Detections_1km_532',
    units       => '-',
    description => 'Number of 1km resolution profiles within each 5 km resolution profile where surface was detected for the 532nm channel.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Detections_333m_1064',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 1064nm channel.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Detections_333m_532',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 532nm channel.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 1064nm of the surface return from Surface Top Altitude 1064 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 532nm of the surface return from Surface Top Altitude 532 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Overlying_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 1064 from 40km to one range bin above Surface Top Altitude 1064.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 532 from 40km to one range bin above Surface Top Altitude 532.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Peak_Signal_1064',
    units       => '/km/sr',
    description => 'Maximum 1064nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Peak_Signal_532',
    units       => '/km/sr',
    description => 'Maximum 532nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Scaled_RMS_Background_1064',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 1064nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Scaled_RMS_Background_532',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 532nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Top_Altitude_1064',
    units       => 'km',
    description => 'Top altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmAPro.Surface_Top_Altitude_532',
    units       => 'km',
    description => 'Top altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_05kmAPro.Surface_Wind_Speeds_02m',
#    units       => 'm/s',
#    description => 'Surface wind speeds measured 2m above the Earth surface. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v. Provided by MERRA-2.'
#  },
#  {
#    name        => 'L2_05kmAPro.Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Surface wind speeds measured 10m above the Earth surface. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v. Provided by MERRA-2.'
#  },
  {
    name        => 'L2_05kmAPro.Temperature',
    units       => 'degC',
    description => 'Temperature, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmAPro.Total_Backscatter_Coefficient_532',
    units       => 'sr^-1km^-1',
    description => 'Total_Backscatter_Coefficient_532'
  },
  {
    name        => 'L2_05kmAPro.Total_Backscatter_Coefficient_Uncertainty_532',
    units       => 'sr^-1km^-1',
    description => 'Total_Backscatter_Coefficient_Uncertainty_532'
  },
  {
    name        => 'L2_05kmAPro.Tropopause_Height',
    units       => 'km',
    description => 'Tropopause height. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmAPro.Tropopause_Temperature',
    units       => 'degC',
    description => 'Tropopause temperature. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmAPro.Unique_Layer_ID',
    units       => '-',
    description => 'Unique number for each layer detected at 5-km or coarser resolution.'
  },



  # L2_05kmCLay ---------------------------------------------------------------

  {
    name        => 'L2_05kmCLay.Column_Feature_Fraction',
    units       => '-',
    description => 'The fraction of the 5-km horizontally averaged profile, between 30-km and the DEM surface elevation, which has been identified as containing a feature (i.e., either a cloud, an aerosol, or a stratospheric layer).'
  },
  {
    name        => 'L2_05kmCLay.Column_IAB_Cumulative_Probability',
    units       => '-',
    description => 'The cumulative probability of measuring a total column integrated attenuated backscatter value equal to the value computed for the current profile.'
  },
  {
    name        => 'L2_05kmCLay.Column_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'The integral with respect to altitude of the 532nm total attenuated backscatter coefficients. The limits of integration are from the onset of the backscatter signal at ~40km, down to the range bin immediately prior to the surface elevation specified by the digital elevation map (DEM).'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Cloud_532',
    units       => '-',
    description => 'Optical depth of all clouds detected with in a 5km averaged profile, obtained by integrating the 532nm cloud extinction profile reported in the CALIPSO 5km Cloud Profile Products.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Cloud_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Cloud 532 parameter.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Stratospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5-km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 1064nm in the CALIPSO 5 km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Stratospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5-km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 532nm in the CALIPSO 5 km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 1064 parameter.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 532 parameter.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Tropospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols detected within a 5-km averaged profile, obtained by integrating the 1064nm aerosol extinction profile reported in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Tropospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols detected within a 5-km averaged profile, obtained by integrating the 532nm aerosol extinction profile reported in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric Aerosol 1064 parameter.'
  },
  {
    name        => 'L2_05kmCLay.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric Aerosol 532 parameter.'
  },
  {
    name        => 'L2_05kmCLay.Column_Particulate_Optical_Depth_Above_Opaque_Water_Cloud_532',
    units       => '-',
    description => '532nm particulate optical depths in the vertical column above an opaque water cloud, derived according to the method detailed in Hu et al., 2007 (DOI: 10.1109/LGRS.2007.901085).'
  },
  {
    name        => 'L2_05kmCLay.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L2_05kmCLay.Extinction_QC_Flag_532',
    units       => '-',
    description => 'The 532nm extinction Quality Control flags are bit-mapped 16-bit integers, reported for each layer for which an extinction retrieval was attempted; least significant bit = bit 1 (bit; interpretation): bit 1; Retrieval, 0 = Unconstrained retrieval, initial lidar ratio unchanged during solution process, 1 = Constrained retrieval, bit 2; Initial lidar ratio reduced to achieve successful full-column solutions, bit 3; Suspicious retrieval due to layer or overlying integrated attenuated backscatter being too high or excessive lidar ratio reductions, bit 4; Initial lidar ratio has been reduced and has converged, but backscatter uncertainty solution does not exist, bit 5; Layer being analyzed has been identified by the feature finder as being totally attenuating (i.e., opaque), bit 6; Estimated optical depth error exceeds the maximum allowable value, bit 7; Negative signal anomaly detected, bit 8; Retrieval terminated at maximum iterations for a constrained retrieval, bit 9; No solution possible within allowable lidar ratio bounds, bit 10; Two-way particulate transmittance has converged but constrained retrieval still not achieved, bit 11; Backscatter coefficients not converging and maximum lidar ratio correction iterations reached, bit 12; Uncertainties not converging and maximum lidar ratio correction iterations achieved, bit 13; Lidar ratio converged but retrieval still not converging, bit 14; Range bin contains two different feature types that may have different lidar ratios, bit 15; Complex retrieval failure, bit 16; Fill value or no solution attempted.'
  },
  {
    name        => 'L2_05kmCLay.Feature_Classification_Flags',
    units       => '-',
    description => 'For each layer detected in the CALIPSO backscatter data a feature classification flags are derived. The complete set of flags are stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): 011 bits 1-3; Feature Type, 0110 = invalid (bad or missing data), 0111 = clear air, 0112 = cloud, 0113 = tropospheric aerosol, 0114 = stratospheric aerosol, 0115 = surface, 0116 = subsurface, 0117 = no signal, bits 4-5; Feature Type QA, 011 = 0 = none, 0111 = low, 0112 = medium, 0113 = high 011 bits 6-7; Ice/Water Phase (cloud), 0110 = unknown/not determined, 0111 = ice, 0112 = water, 0113 = oriented ice crystals, bits 8-9; Ice/Water Phase QA (cloud), 0110 = none, 0111 = low, 0112 = medium, 0113 = high, 011bits 10-12; Feature Subtype (clear-air), 0110 = n/a, 0111 = not searched for features at 80km resolution, 0112 = not searched for features at 20 and 8km resolution, bits 10-12; Feature Subtype (Cloud), 0110 = low overcast, transparent, 0111 = low overcast, opaque, 0112 = transition stratocumulus, 0113 = low, broken cloud.'
  },
  {
    name        => 'L2_05kmCLay.FeatureFinderQC',
    units       => '-',
    description => 'To generate data at a nominal 5-km horizontal resolution requires averaging 15 consecutive single-shot laser pulses. For each 5km average a set of feature finder QC flags are reported, which are a set of 15 Boolean values that report whether or not a feature (cloud, aerosol, or surface echo) was detected in each of the 15 laser pulses. The flags are implemented as a 16-bit integer, with the most significant bit unused.'
  },
  {
    name        => 'L2_05kmCLay.Feature_Optical_Depth_532',
    units       => '-',
    description => 'Reports estimates of 532nm layer optical depth.'
  },
  {
    name        => 'L2_05kmCLay.Feature_Optical_Depth_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty, reported as an absolute error, in the 532nm layer optical  depth.'
  },
  {
    name        => 'L2_05kmCLay.Final_532_Lidar_Ratio',
    units       => 'sr',
    description => 'Reports the lidar ratio in use at the conclusion of the extinction processing for each layer. The final lidar ratio may be (1) the initial lidar ratio supplied by the Scene Classification Algorithms, (2) the result of modifications to this initial lidar ratio to avoid a non-physical solution, or (3) a lidar ratio determined from a measured layer transmittance.'
  },
  {
    name        => 'L2_05kmCLay.Final_532_Lidar_Ratio_Uncertainty',
    units       => 'sr',
    description => 'Reports the lidar ratio uncertainty at the conclusion of the extinction processing for each layer.'
  },
  {
    name        => 'L2_05kmCLay.High_Resolution_Layers_Cleared',
    units       => '-',
    description => 'Bit-mapped integer identifying which single-shot resolution columns in a 5km average contain clouds detected and removed by the boundary layer cloud-clearing algorithm. The flags are implemented as a 16-bit integer, with the most significant bit unused.'
  },
  {
    name        => 'L2_05kmCLay.Horizontal_Averaging',
    units       => 'km',
    description => 'The amount of horizontal averaging required for a feature to be detected.'
  },
  {
    name        => 'L2_05kmCLay.Ice_Water_Path',
    units       => 'g/m2',
    description => 'The integral, from layer top to layer bas, of the ice-water content profile within any ice cloud layer. Ice water content is derived from the cloud particulate 532nm extinction coefficient using a temperature-dependent parameterization derived from in-situ measurements.'
  },
  {
    name        => 'L2_05kmCLay.Ice_Water_Path_Uncertainty',
    units       => 'g/m2',
    description => 'Uncertainty associated with the estimate of ice water path, computed by integrating the uncertainties computed for the ice water content within the layer.'
  },
  {
    name        => 'L2_05kmCLay.IGBP_Surface_Type',
    units       => '-',
    description => 'IGBP surface type at the laser footprint provided by the CERES/SARB surface map: 1 = Evergreen-Needleleaf-Forest, 2 = Evergreen-Broadleaf-Forest, 3 = Deciduous-Needleleaf-Forest, 4 = Deciduous-Broadleaf-Forest, 5 = Mixed-Forest, 6 = Closed-Shrubland, 7 = Open-Shrubland/Desert, 8 = Woody-Savanna, 9 = Savanna, 10 = Grassland, 11 = Wetland, 12 = Cropland, 13 = Urban, 14 = Crop-Mosaic, 15 = Permanent-Snow, 16 = Barren/Desert, 17 = Water, 18 = Tundra.'
  },
  {
    name        => 'L2_05kmCLay.Final_1064_Lidar_Ratio',
    units       => 'sr',
    description => 'Reports the lidar ratio in use at the conclusion of the extinction processing for each layer. The final lidar ratio may be (1) the initial lidar ratio supplied by the Scene Classification Algorithms, (2) the result of modifications to this initial lidar ratio to avoid a non-physical solution, or (3) a lidar ratio determined from a measured layer transmittance.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'The 532 nm integrated attenuated backscatter (IAB) computed for each layer after correcting the backscatter signal for 532nm attenuation by any overlying layers.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Attenuated_Backscatter_1064',
    units       => '1/sr',
    description => 'The 1064nm integrated attenuated backscatter (IAB) computed for each layer.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Attenuated_Backscatter_Uncertainty_1064',
    units       => '1/sr',
    description => 'The uncertainties reported for the 1064nm integrated attenuated backscatters, which provide an estimate of the random error in the backscatter signal.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Attenuated_Backscatter_Uncertainty_532',
    units       => '1/sr',
    description => 'Estimated uncertainties the 532 nm integrated attenuated backscatter values. These uncertainties quantify the random error in the integral of the attenuation-corrected 532nm backscatter signal.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Attenuated_Total_Color_Ratio',
    units       => '-',
    description => 'Integrated attenuated total color ratio computed for each layer.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Attenuated_Total_Color_Ratio_Uncertainty',
    units       => '-',
    description => 'Integrated_Attenuated_Total_Color_Ratio_Uncertainty'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Particulate_Depolarization_Ratio',
    units       => '-',
    description => 'Layer-integrated particulate depolarization ratio. The particulate depolarization ratio represents the contribution to the volume depolarization ratio that is due only to the cloud and/or aerosol particles within the layer.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Particulate_Depolarization_Ratio_Uncertainty',
    units       => '-',
    description => 'Uncertainty associated with the estimate of the layer-integrated particulate depolarization ratio, calculated by integrating the uncertainties for the particulate backscatter coefficients measured in the 532nm parallel and perpendicular channels, and then applying standard techniques for error propagation in ratioed quantities.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Volume_Depolarization_Ratio',
    units       => '-',
    description => 'The layer integrated 532 nm volume depolarization ratio (IVDR) computed for each layer. The quality of the estimate for IVDR is determined by the accuracy of the top and base identification, the reliability of the polarization gain ratio calibration, and by the signal-to-noise ratio (SNR) of the backscatter data within the layer.'
  },
  {
    name        => 'L2_05kmCLay.Integrated_Volume_Depolarization_Ratio_Uncertainty',
    units       => '-',
    description => 'The uncertainties reported for the 532nm layer-integrated volume depolarization ratios provide an estimate of the total random error in the combined backscatter signals.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Base_Altitude',
    units       => 'km',
    description => 'Base of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Base_Extended',
    units       => '-',
    description => 'A non-zero value indicates that the base of the layer has been extended from the initial altitude assigned by the layer detection algorithm to a new, lower altitude lying three range bins (90 m) above the Earth surface. Non-zero values represent the layer feature classification flags prior to the base being extended.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Base_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer base altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Base_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer attenuated backscatter centroid altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Centroid_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer base altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Effective_532_Multiple_Scattering_Factor',
    units       => '-',
    description => 'The layer effective 532 nm multiple scattering factor corresponds to the limit of single scattering only, with smaller values indicating increasing contributions to the backscatter signal from multiple scattering.'
  },
  {
    name        => 'L2_05kmCLay.Layer_IAB_QA_Factor',
    units       => '-',
    description => 'Single layer analog of the column integrated attenuated backscatter cumulative probability.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Top_Altitude',
    units       => 'km',
    description => 'Top of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Top_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCLay.Layer_Top_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCLay.Lidar_Ratio_532_Selection_Method',
    units       => '-',
    description => 'Specifies the internal procedure used to select the initial lidar ratio for each layer: - not determined, - Transmittance Method, constrained retrieval, - Opaque Layer Method, based on opaque clouds or aerosols, - Cirrus Temperature Method, unconstrained ice clouds, - Unknown Cloud Phase, unconstrained unknown phase clouds with centroid temperature between -37 and 0 degrees C, - Fixed Layer Method, unconstrained water clouds and aerosols, - Solution at 532nm Method, unconstrained 1064 nm cloud only when 532 nm extinction solution is available\012\0117 - MAC Table Interpolation, unconstrained aerosol using MAC interpolated lookup table.'
  },
  {
    name        => 'L2_05kmCLay.Low_Energy_Mitigation_Column_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each column; least significant bit = bit 1 (bit; interpretation): bit 1; Column contains LEM affected data (data has been rejected or low energy shots that LEM accepts, bit 2; Column belongs to a 5km resolution frame that has been rejected due to too many unusable profiles, bit 3; Column belongs to a 5 km resolution frame that has been rejected due to too many rejected subregions in altitude region 3, bit 4; Column belongs to a 5km resolution frame that has been rejected due to too many rejected subregions in altitude region 4, bit 5; Feature detection at 20km resolution not performed in this column due to insufficient fraction of usable 5km resolution frames, bit 6; Feature detection at 80 km resolution not performed in this column due to insufficient fraction of usable 5 km resolution frames, bits 7-16; Unused.'
  },
  {
    name        => 'L2_05kmCLay.Low_Energy_Mitigation_Feature_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each layer detected (flag - interpretation): 0 - No low energy shots affecting the layer, 1 - Some portion of the layer contains low energy shots and/or data rejected by the LEM algorithm'
  },
  {
    name        => 'L2_05kmCLay.Measured_Two_Way_Transmittance_532',
    units       => '-',
    description => 'Provides the measured value of the layer two-way transmittance, derived by computing the ratio of the mean attenuated scattering ratios in the "clear air" regions immediately below and above the layer, for isolated transparent layers.'
  },
  {
    name        => 'L2_05kmCLay.Measured_Two_Way_Transmittance_Uncertainty_532',
    units       => '-',
    description => 'The relative error in the two-way transmittance measurement, calculated using standard techniques for error propagation in ratioed quantities.'
  },
  {
    name        => 'L2_05kmCLay.Midlayer_Pressure',
    units       => 'hPa',
    description => 'Pressure at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCLay.Midlayer_Temperature',
    units       => 'degC',
    description => 'Temperature at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCLay.Minimum_Laser_Energy_532',
    units       => 'J',
    description => 'Minimum laser energy at 532 nm measured within each 80 km along-track data segment (80km = 240 single-shot laser pulses). The 80km distance matches the largest horizontal extent considered in CALIOP standard level 2 data analyses.'
  },
# unmatched dimensionw
#  {
#    name        => 'L2_05kmCLay.Normalization_Constant_Uncertainty',
#    units       => '-',
#    description => 'Uncertainty, reported as a relative error, in the 532nm and 1064nm calibration constants due solely to random error in the backscatter measurements in the calibration region.'
#  },
  {
    name        => 'L2_05kmCLay.Number_Layers_Found',
    units       => '-',
    description => 'The number of cloud layers found in a column. The merged 5km layer product reports layers detected at multiple averaging resolutions (5km, 20km, and 80km). Because the reporting resolution (5-km) is not always identical to the detection resolution, layers may appear to overlap in the vertical dimension.'
  },
  {
    name        => 'L2_05kmCLay.ODCOD_Effective_Optical_Depth_532',
    units       => '-',
    description => 'Column estimate of apparent (multiple scattering effects not accounted for) optical depth at 532nm wavelength for all particulate (clouds and aerosols) within the illuminated volume from just above the lidar detected ocean surface to the CALIPSO top of atmosphere.'
  },
  {
    name        => 'L2_05kmCLay.ODCOD_Effective_Optical_Depth_532_Uncertainty',
    units       => '-',
    description => 'An analytical estimate of the absolute standard deviation uncertainty in the ODCOD Effective Optical Depth estimate at 532nm.'
  },
  {
    name        => 'L2_05kmCLay.ODCOD_QC_Flag_532',
    units       => '-',
    description => 'Ocean Derived Column Optical Depth (ODCOD) quality control flags; least significant bit = bit 1 (bit; interpretation): bit 1; The fundamental time delay was shifted from the first point, bit 2; Lidar surface data too many data points for a valid surface spike, bit 3; Range bin above recorded lidar surface data added, bit 4; Range bin below recorded lidar surface data added, bit 5; The surface point ratio indicated that the first point of the surface spike is missing, bit 6; To find the time delay, the surface spike index range had to be shifted down, bit 7; The retrieval inputs fall outside of confidence ranges indicating a confident retrieval, bits 8-10; Unused, none, bit 11; No Lidar surface was found, bit 12; The surface is not ocean, bit 13; The surface is ice, bit 14; The wind speed is outside the valid range; 0.025 to 43 m/s, bit 15; The fundamental time delay found was too low and there were not enough points to find a new one, bit 16; Too few data points to solve, bit 17; The low pass filter impulse response model curve area was too large, bit 18; The low pass filter impulse response model curve scale factor could no t be found, bit 19; Lidar surface saturated, bit 20; Lidar surface has negative signal anomaly, bit 21; Lidar surface data is not valid, bit 22; Flags bad input/ancillary data, bit 23; No Lidar surface was detected via derivative method at greater than single shot resolutions, bits 24-32; Unused, none.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_05kmCLay.ODCOD_Surface_Wind_Speed_Correction',
#    units       => 'm/s',
#    description => 'An additive correction derived from AMSR analysis applied to the MERRA-2 wind speed magnitude calculated from the zonal winds reported in ODCOD Surface Wind Speeds 10m.'
#  },
#  {
#    name        => 'L2_05kmCLay.ODCOD_Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Mean zonal and meridional component of the surface wind speed 10 meters above the surface provided by MERRA-2. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v.'
#  },
  {
    name        => 'L2_05kmCLay.Off_Nadir_Angle',
    units       => 'deg',
    description => 'Angle of the lidar viewing vector; 0.3 degrees off-nadir from launch/first-light to 3.0 degrees after November 28, 2007.'
  },
  {
    name        => 'L2_05kmCLay.Opacity_Flag',
    units       => '-',
    description => 'Flag that identifies a layer in which the backscatter signal becomes completely attenuated. An opacity flag of 0 indicates a transparent layer, while a value of 1 indicates an opaque layer. A layer is considered opaque if (a) it is the lowest feature detected in a column, and (b) it is not subsequently classified as a surface return. For those layers identified as opaque, the reported base altitude must be considered as an apparent base, rather than a true base.'
  },
  {
    name        => 'L2_05kmCLay.Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'Overlying_Integrated_Attenuated_Backscatter_532'
  },
  {
    name        => 'L2_05kmCLay.Parallel_Column_Reflectance_532',
    units       => '-',
    description => 'Bi-directional column reflectance derived from the root-mean-square variation of the 532 nm parallel channel baseline measurements.'
  },
  {
    name        => 'L2_05kmCLay.Parallel_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Column reflectance uncertainty for the 532nm parallel channel.'
  },
  {
    name        => 'L2_05kmCLay.Perpendicular_Column_Reflectance_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_532'
  },
  {
    name        => 'L2_05kmCLay.Perpendicular_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_Uncertainty_532'
  },
# unmatched dimensions
#  {
#    name        => 'L2_05kmCLay.Profile_ID',
#    units       => '-',
#    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
#  },
#  {
#    name        => 'L2_05kmCLay.Profile_Time',
#    units       => 'seconds_since_1993-01-01',
#    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
#  },
  {
    name        => 'L2_05kmCLay.Scattering_Angle',
    units       => 'deg',
    description => 'Angle between the lidar viewing angle and the line of sight of the sun.'
  },
  {
    name        => 'L2_05kmCLay.Scene_Flag',
    units       => '-',
    description => 'Bit mapped flag where the following bits indicate the atmospheric feature type was detected at some resolution in the colum n; least significant bit = bit 1 (bit; interpretation): bit 1; tropospheric aerosol subtype clean marine, bit 2; tropospheric aerosol subtype dust, bit 3; tropospheric aerosol subtype polluted continental/smoke, bit 4; tropospheric aerosol subtype clean continental, bit 5; tropospheric aerosol subtype polluted dust, bit 6; tropospheric aerosol subtype elevated smoke, bit 7; tropospheric aerosol subtype dusty marine, bit 8; stratospheric aerosol subtype polar stratospheric aerosol, bit 9; stratospheric aerosol subtype volcanic ash, bit 10; stratospheric aerosol subtype sulfate, bit 11; stratospheric aerosol subtype elevated smoke bit 12; stratospheric aerosol unclassified. bit 13; On if column has ROI clouds at any resolution. bit 14; On if column has HOI clouds at any resolution. bit 15; On if column has water clouds at any resolution. bit 16; On if column has unknown phase clouds at any resolution.'
  },
  {
    name        => 'L2_05kmCLay.Single_Shot_Cloud_Cleared_Fraction',
    units       => '-',
    description => 'Layers detected in the planetary boundary layer (PBL) are subjected to an additional cloud clearing procedure to separate small-scale boundary layer clouds from any surrounding PBL aerosols at the highest possible spatial resolution (i.e., single shot data). The single shot cloud cleared fraction reports the fraction of the nominal layer area (i.e., horizontal averaging distance times layer height) that was removed by the cloud clearing process.'
  },
  {
    name        => 'L2_05kmCLay.Snow_Ice_Surface_Type',
    units       => '-',
    description => 'Snow and ice coverage for the surface at the laser footprint provided by the U.S. Airforce Weather Agency (AFWA). Values between 0 and 100 represent the percentage of ice, while 101 = permanent ice, 103 = snow, and 255 = mixed pixels at coastline.'
  },
  {
    name        => 'L2_05kmCLay.Solar_Azimuth_Angle',
    units       => 'degrees',
    description => 'Azimuth angle from north of the line of sight to the sun.'
  },
  {
    name        => 'L2_05kmCLay.Solar_Zenith_Angle',
    units       => 'degrees',
    description => 'Angle between the zenith at the lidar footprint on the surface and the line of sight to the sun.'
  },
  {
    name        => 'L2_05kmCLay.Spacecraft_Position_X',
    units       => 'km',
    description => 'X-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_05kmCLay.Spacecraft_Position_Y',
    units       => 'km',
    description => 'Y-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_05kmCLay.Spacecraft_Position_Z',
    units       => 'km',
    description => 'Z-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_05kmCLay.Surface_1064_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 1064 to Surface Base 1064, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCLay.Surface_1064_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 1064 to Surface Base Altitude 1064, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCLay.Surface_532_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCLay.Surface_532_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Base_Altitude_1064',
    units       => 'km',
    description => 'Base altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Base_Altitude_532',
    units       => 'km',
    description => 'Base altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Detection_Flags_1064',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Detection_Flags_532',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Detections_333m_1064',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 1064nm channel.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Detections_333m_532',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 532nm channel.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 1064nm of the surface return from Surface Top Altitude 1064 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 532nm of the surface return from Surface Top Altitude 532 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Overlying_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 1064 from 40km to one range bin above Surface Top Altitude 1064.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 532 from 40km to one range bin above Surface Top Altitude 532.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Peak_Signal_1064',
    units       => '/km/sr',
    description => 'Maximum 1064nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Peak_Signal_532',
    units       => '/km/sr',
    description => 'Maximum 532nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Scaled_RMS_Background_1064',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 1064nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Scaled_RMS_Background_532',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 532nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Top_Altitude_1064',
    units       => 'km',
    description => 'Top altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmCLay.Surface_Top_Altitude_532',
    units       => 'km',
    description => 'Top altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmCLay.Tropopause_Height',
    units       => 'km',
    description => 'Tropopause height. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmCLay.Tropopause_Temperature',
    units       => 'degC',
    description => 'Tropopause temperature. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmCLay.Unique_Layer_ID',
    units       => '-',
    description => 'Unique number for each layer detected at 5-km or coarser resolution.'
  },



  # L2_05kmCPro ---------------------------------------------------------------

  {
    name        => 'L2_05kmCPro.Aerosol_Layer_Fraction',
    units       => '-',
    description => 'Aerosol_Layer_Fraction'
  },
  {
    name        => 'L2_05kmCPro.Atmospheric_Volume_Description',
    units       => '-',
    description => 'Atmospheric_Volume_Description'
  },
  {
    name        => 'L2_05kmCPro.Cloud_Layer_Fraction',
    units       => '-',
    description => 'Cloud_Layer_Fraction'
  },
  {
    name        => 'L2_05kmCPro.Cloud_Multiple_Scattering_Profile_532',
    units       => '-',
    description => 'Cloud_Multiple_Scattering_Profile_532'
  },
  {
    name        => 'L2_05kmCPro.Column_Feature_Fraction',
    units       => '-',
    description => 'The fraction of the 5-km horizontally averaged profile, between 30km and the DEM surface elevation, which has been identified as containing a feature (i.e., either a cloud, an aerosol, or a stratospheric layer).'
  },
  {
    name        => 'L2_05kmCPro.Column_IAB_Cumulative_Probability',
    units       => '-',
    description => 'The cumulative probability of measuring a total column integrated attenuated backscatter value equal to the value computed for the current profile.'
  },
  {
    name        => 'L2_05kmCPro.Column_Integrated_Attenuated_Backscatter_532',
    units       => 'sr^-1',
    description => 'The integral with respect to altitude of the 532nm total attenuated backscatter coefficients. The limits of integration are from the onset of the backscatter signal at ~40 km, down to the range bin immediately prior to the surface elevation specified by the digital elevation map (DEM).'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Cloud_532',
    units       => '-',
    description => 'Column_Optical_Depth_Cloud_532'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Cloud_Uncertainty_532',
    units       => '-',
    description => 'Column_Optical_Depth_Cloud_Uncertainty_532'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Stratospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 1064nm in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Stratospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all stratospheric aerosols within a 5km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 532nm in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 1064 parameter.'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Stratospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Stratospheric 532 parameter.'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Tropospheric_Aerosols_1064',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols within a 5km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 1064nm in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Tropospheric_Aerosols_532',
    units       => '-',
    description => 'Optical depth of all tropospheric aerosols within a 5km averaged profile, obtained by integrating the stratospheric particulate extinction coefficients reported at 532nm in the CALIPSO 5km Aerosol Profile Products.'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_1064',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric 1064 parameter.'
  },
  {
    name        => 'L2_05kmCPro.Column_Optical_Depth_Tropospheric_Aerosols_Uncertainty_532',
    units       => '-',
    description => 'Estimated uncertainty in the Column Optical Depth Tropospheric 532 parameter.'
  },
  {
    name        => 'L2_05kmCPro.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L2_05kmCPro.Extinction_Coefficient_532',
    units       => 'km^-1',
    description => 'Particulate 532nm extinction coefficients reported for each profile range bin in which the appropriate particulates were detected.'
  },
  {
    name        => 'L2_05kmCPro.Extinction_Coefficient_532_Uncertainty',
    units       => 'km^-1',
    description => 'Uncertainties reported for the particulate 532nm extinction coefficients for each profile range bin.'
  },
  {
    name        => 'L2_05kmCPro.Extinction_QC_Flag_532',
    units       => '-',
    description => 'The 532nm extinction Quality Control flags are bit-mapped 16-bit integers, reported for each layer for which an extinction retrieval was attempted; least significant bit = bit 1 (bit; interpretation): bit 1; Retrieval, 0 = Unconstrained retrieval, initial lidar ratio unchanged during solution process, 1 = Constrained retrieval, bit 2; Initial lidar ratio reduced to achieve successful full-column solutions, bit 3; Suspicious retrieval due to layer or overlying integrated attenuated backscatter being too high or excessive lidar ratio reductions, bit 4; Initial lidar ratio has been reduced and has converged, but backscatter uncertainty solution does not exist, bit 5; Layer being analyzed has been identified by the feature finder as being totally attenuating (i.e., opaque), bit 6; Estimated optical depth error exceeds the maximum allowable value, bit 7; Negative signal anomaly detected, bit 8; Retrieval terminated at maximum iterations for a constrained retrieval, bit 9; No solution possible within allowable lidar ratio bounds, bit 10; Two-way particulate transmittance has converged but constrained retrieval still not achieved, bit 11; Backscatter coefficients not converging and maximum lidar ratio correction iterations reached, bit 12; Uncertainties not converging and maximum lidar ratio correction iterations achieved, bit 13; Lidar ratio converged but retrieval still not converging, bit 14; Range bin contains two different feature types that may have different lidar ratios, bit 15; Complex retrieval failure, bit 16; Fill value or no solution attempted.'
  },
  {
    name        => 'L2_05kmCPro.Ice_Water_Content_Profile',
    units       => 'gm^-3',
    description => 'Ice_Water_Content_Profile'
  },
  {
    name        => 'L2_05kmCPro.Ice_Water_Content_Profile_Uncertainty',
    units       => 'gm^-3',
    description => 'Ice_Water_Content_Profile_Uncertainty'
  },
  {
    name        => 'L2_05kmCPro.IGBP_Surface_Type',
    units       => '-',
    description => 'IGBP surface type at the laser footprint provided by the CERES/SARB surface map: 1 = Evergreen-Needleleaf-Forest, 2 = Evergreen-Broadleaf-Forest, 3 = Deciduous-Needleleaf-Forest, 4 = Deciduous-Broadleaf-Forest, 5 = Mixed-Forest, 6 = Closed-Shrubland, 7 = Open-Shrubland/Desert, 8 = Woody-Savanna, 9 = Savanna, 10 = Grassland, 11 = Wetland, 12 = Cropland, 13 = Urban, 14 = Crop-Mosaic, 15 = Permanent-Snow, 16 = Barren/Desert, 17 = Water, 18 = Tundra.'
  },
  {
    name        => 'L2_05kmCPro.Low_Energy_Mitigation_Column_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each column; least significant bit = bit 1 (bit; interpretation): bit 1; Column contains LEM affected data (data has been rejected or low energy shots that LEM accepts, bit 2; Column belongs to a 5km resolution frame that has been rejected due to too many unusable profiles, bit 3; Column belongs to a 5 km resolution frame that has been rejected due to too many rejected subregions in altitude region 3, bit 4; Column belongs to a 5km resolution frame that has been rejected due to too many rejected subregions in altitude region 4, bit 5; Feature detection at 20km resolution not performed in this column due to insufficient fraction of usable 5km resolution frames, bit 6; Feature detection at 80 km resolution not performed in this column due to insufficient fraction of usable 5 km resolution frames, bits 7-16; Unused.'
  },
  {
    name        => 'L2_05kmCPro.Low_Energy_Mitigation_Feature_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each layer detected (flag - interpretation): 0 - No low energy shots affecting the layer, 1 - Some portion of the layer contains low energy shots and/or data rejected by the LEM algorithm'
  },
  {
    name        => 'L2_05kmCPro.Minimum_Laser_Energy_532',
    units       => 'J',
    description => 'Minimum laser energy at 532 nm measured within each 80 km along-track data segment (80km = 240 single-shot laser pulses). The 80km distance matches the largest horizontal extent considered in CALIOP standard level 2 data analyses.'
  },
  {
    name        => 'L2_05kmCPro.Molecular_Number_Density',
    units       => 'm^-3',
    description => 'Mean molecular number density reported for the midpoint of each range bin in the profile; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_05kmCPro.ODCOD_Effective_Optical_Depth_532',
    units       => '-',
    description => 'Column estimate of apparent (multiple scattering effects not accounted for) optical depth at 532nm wavelength for all particulate (clouds and aerosols) within the illuminated volume from just above the lidar detected ocean surface to the CALIPSO top of atmosphere.'
  },
  {
    name        => 'L2_05kmCPro.ODCOD_Effective_Optical_Depth_532_Uncertainty',
    units       => '-',
    description => 'An analytical estimate of the absolute standard deviation uncertainty in the ODCOD Effective Optical Depth estimate at 532nm.'
  },
  {
    name        => 'L2_05kmCPro.ODCOD_QC_Flag_532',
    units       => '-',
    description => 'Ocean Derived Column Optical Depth (ODCOD) quality control flags; least significant bit = bit 1 (bit; interpretation): bit 1; The fundamental time delay was shifted from the first point, bit 2; Lidar surface data too many data points for a valid surface spike, bit 3; Range bin above recorded lidar surface data added, bit 4; Range bin below recorded lidar surface data added, bit 5; The surface point ratio indicated that the first point of the surface spike is missing, bit 6; To find the time delay, the surface spike index range had to be shifted down, bit 7; The retrieval inputs fall outside of confidence ranges indicating a confident retrieval, bits 8-10; Unused, none, bit 11; No Lidar surface was found, bit 12; The surface is not ocean, bit 13; The surface is ice, bit 14; The wind speed is outside the valid range; 0.025 to 43 m/s, bit 15; The fundamental time delay found was too low and there were not enough points to find a new one, bit 16; Too few data points to solve, bit 17; The low pass filter impulse response model curve area was too large, bit 18; The low pass filter impulse response model curve scale factor could no t be found, bit 19; Lidar surface saturated, bit 20; Lidar surface has negative signal anomaly, bit 21; Lidar surface data is not valid, bit 22; Flags bad input/ancillary data, bit 23; No Lidar surface was detected via derivative method at greater than single shot resolutions, bits 24-32; Unused, none.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_05kmCPro.ODCOD_Surface_Wind_Speed_Correction',
#    units       => 'm/s',
#    description => 'An additive correction derived from AMSR analysis applied to the MERRA-2 wind speed magnitude calculated from the zonal winds reported in ODCOD Surface Wind Speeds 10m.'
#  },
#  {
#    name        => 'L2_05kmCPro.ODCOD_Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Mean zonal and meridional component of the surface wind speed 10 meters above the surface provided by MERRA-2. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v.'
#  },
  {
    name        => 'L2_05kmCPro.Ozone_Number_Density',
    units       => 'molecules/m3',
    description => 'Ozone number density, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmCPro.Particulate_Depolarization_Ratio_Profile_532',
    units       => '-',
    description => 'Particulate_Depolarization_Ratio_Profile_532'
  },
  {
    name        => 'L2_05kmCPro.Particulate_Depolarization_Ratio_Uncertainty_532',
    units       => '-',
    description => 'Particulate_Depolarization_Ratio_Uncertainty_532'
  },
  {
    name        => 'L2_05kmCPro.Perpendicular_Backscatter_Coefficient_532',
    units       => 'sr^-1km^-1',
    description => 'Perpendicular_Backscatter_Coefficient_532'
  },
  {
    name        => 'L2_05kmCPro.Perpendicular_Backscatter_Coefficient_Uncertainty_532',
    units       => 'sr^-1km^-1',
    description => 'Perpendicular_Backscatter_Coefficient_Uncertainty_532'
  },
  {
    name        => 'L2_05kmCPro.Pressure',
    units       => 'hPa',
    description => 'Mean pressure reported for the midpoint of each range bin in the profile; derived from the MERRA-2 data product.'
  },
# unmatched dimensions
#  {
#    name        => 'L2_05kmCPro.Profile_ID',
#    units       => '-',
#    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
#  },
#  {
#    name        => 'L2_05kmCPro.Profile_Time',
#    units       => 'seconds_since_1993-01-01',
#    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
#  },
  {
    name        => 'L2_05kmCPro.Relative_Humidity',
    units       => '%',
    description => 'Relative humidity, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmCPro.Surface_1064_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 1064 to Surface Base 1064, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCPro.Surface_1064_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 1064 to Surface Base Altitude 1064, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCPro.Surface_532_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCPro.Surface_532_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Base_Altitude_1064',
    units       => 'km',
    description => 'Base altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Base_Altitude_532',
    units       => 'km',
    description => 'Base altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Detection_Flags_1064',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Detection_Flags_532',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Detections_333m_1064',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 1064nm channel.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Detections_333m_532',
    units       => '-',
    description => 'Number of single shot profiles (1/3 km resolution) within each 5 km resolution profile where surface was detected for the 532nm channel.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 1064nm of the surface return from Surface Top Altitude 1064 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 532nm of the surface return from Surface Top Altitude 532 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Overlying_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 1064 from 40km to one range bin above Surface Top Altitude 1064.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 532 from 40km to one range bin above Surface Top Altitude 532.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Peak_Signal_1064',
    units       => '/km/sr',
    description => 'Maximum 1064nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Peak_Signal_532',
    units       => '/km/sr',
    description => 'Maximum 532nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Scaled_RMS_Background_1064',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 1064nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Scaled_RMS_Background_532',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 532nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Top_Altitude_1064',
    units       => 'km',
    description => 'Top altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_05kmCPro.Surface_Top_Altitude_532',
    units       => 'km',
    description => 'Top altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_05kmCPro.Surface_Wind_Speeds_02m',
#    units       => 'm/s',
#    description => 'Surface wind speeds measured 2m above the Earth surface. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v. Provided by MERRA-2.'
#  },
  {
    name        => 'L2_05kmCPro.Temperature',
    units       => 'degC',
    description => 'Temperature, for each profile, reported at the Met_Data_Altitudes resolution. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmCPro.Total_Backscatter_Coefficient_532',
    units       => 'sr^-1km^-1',
    description => 'Total_Backscatter_Coefficient_532'
  },
  {
    name        => 'L2_05kmCPro.Total_Backscatter_Coefficient_Uncertainty_532',
    units       => 'sr^-1km^-1',
    description => 'Total_Backscatter_Coefficient_Uncertainty_532'
  },
  {
    name        => 'L2_05kmCPro.Tropopause_Height',
    units       => 'km',
    description => 'Tropopause height. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmCPro.Tropopause_Temperature',
    units       => 'degC',
    description => 'Tropopause temperature. Provided by MERRA-2.'
  },
  {
    name        => 'L2_05kmCPro.Unique_Layer_ID',
    units       => '-',
    description => 'Unique number for each layer detected at
                         5-km or coarser resolution.'
  },



  # L2_333mMLay ---------------------------------------------------------------

  {
    name        => 'L2_333mMLay.Column_Feature_Fraction',
    units       => '-',
    description => 'The fraction of the 5-km horizontally averaged profile, between 30-km and the DEM surface elevation, which has been identified as containing a feature (i.e., either a cloud, an aerosol, or a stratospheric layer).'
  },
  {
    name        => 'L2_333mMLay.Column_IAB_Cumulative_Probability',
    units       => '-',
    description => 'The cumulative probability of measuring a total column integrated attenuated backscatter value equal to the value computed for the current profile.'
  },
  {
    name        => 'L2_333mMLay.Column_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'The integral with respect to altitude of the 532nm total attenuated backscatter coefficients. The limits of integration are from the onset of the backscatter signal at ~40km, down to the range bin immediately prior to the surface elevation specified by the digital elevation map (DEM).'
  },
  {
    name        => 'L2_333mMLay.Column_Particulate_Optical_Depth_Above_Opaque_Water_Cloud_532',
    units       => '-',
    description => '532nm particulate optical depths in the vertical column above an opaque water cloud, derived according to the method detailed in Hu et al., 2007 (DOI: 10.1109/LGRS.2007.901085).'
  },
  {
    name        => 'L2_333mMLay.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L2_333mMLay.Feature_Classification_Flags',
    units       => '-',
    description => 'For each layer detected in the CALIPSO backscatter data a feature classification flags are derived. The complete set of flags are stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): 011 bits 1-3; Feature Type, 0110 = invalid (bad or missing data), 0111 = clear air, 0112 = cloud, 0113 = tropospheric aerosol, 0114 = stratospheric aerosol, 0115 = surface, 0116 = subsurface, 0117 = no signal, bits 4-5; Feature Type QA, 011 = 0 = none, 0111 = low, 0112 = medium, 0113 = high 011 bits 6-7; Ice/Water Phase (cloud), 0110 = unknown/not determined, 0111 = ice, 0112 = water, 0113 = oriented ice crystals, bits 8-9; Ice/Water Phase QA (cloud), 0110 = none, 0111 = low, 0112 = medium, 0113 = high, 011bits 10-12; Feature Subtype (clear-air), 0110 = n/a, 0111 = not searched for features at 80km resolution, 0112 = not searched for features at 20 and 8km resolution, bits 10-12; Feature Subtype (Cloud), 0110 = low overcast, transparent, 0111 = low overcast, opaque, 0112 = transition stratocumulus, 0113 = low, broken cloud.'
  },
  {
    name        => 'L2_333mMLay.IGBP_Surface_Type',
    units       => '-',
    description => 'IGBP surface type at the laser footprint provided by the CERES/SARB surface map: 1 = Evergreen-Needleleaf-Forest, 2 = Evergreen-Broadleaf-Forest, 3 = Deciduous-Needleleaf-Forest, 4 = Deciduous-Broadleaf-Forest, 5 = Mixed-Forest, 6 = Closed-Shrubland, 7 = Open-Shrubland/Desert, 8 = Woody-Savanna, 9 = Savanna, 10 = Grassland, 11 = Wetland, 12 = Cropland, 13 = Urban, 14 = Crop-Mosaic, 15 = Permanent-Snow, 16 = Barren/Desert, 17 = Water, 18 = Tundra.'
  },
  {
    name        => 'L2_333mMLay.Final_1064_Lidar_Ratio',
    units       => 'sr',
    description => 'Reports the lidar ratio in use at the conclusion of the extinction processing for each layer. The final lidar ratio may be (1) the initial lidar ratio supplied by the Scene Classification Algorithms, (2) the result of modifications to this initial lidar ratio to avoid a non-physical solution, or (3) a lidar ratio determined from a measured layer transmittance.'
  },
  {
    name        => 'L2_333mMLay.Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'The 532 nm integrated attenuated backscatter (IAB) computed for each layer after correcting the backscatter signal for 532nm attenuation by any overlying layers.'
  },
  {
    name        => 'L2_333mMLay.Integrated_Attenuated_Backscatter_1064',
    units       => '1/sr',
    description => 'The 1064nm integrated attenuated backscatter (IAB) computed for each layer.'
  },
  {
    name        => 'L2_333mMLay.Integrated_Attenuated_Backscatter_Uncertainty_1064',
    units       => '1/sr',
    description => 'The uncertainties reported for the 1064nm integrated attenuated backscatters, which provide an estimate of the random error in the backscatter signal.'
  },
  {
    name        => 'L2_333mMLay.Integrated_Attenuated_Backscatter_Uncertainty_532',
    units       => '1/sr',
    description => 'Estimated uncertainties the 532 nm integrated attenuated backscatter values. These uncertainties quantify the random error in the integral of the attenuation-corrected 532nm backscatter signal.'
  },
  {
    name        => 'L2_333mMLay.Integrated_Attenuated_Total_Color_Ratio',
    units       => '-',
    description => 'Integrated attenuated total color ratio computed for each layer.'
  },
  {
    name        => 'L2_333mMLay.Integrated_Attenuated_Total_Color_Ratio_Uncertainty',
    units       => '-',
    description => 'Integrated_Attenuated_Total_Color_Ratio_Uncertainty'
  },
  {
    name        => 'L2_333mMLay.Integrated_Volume_Depolarization_Ratio',
    units       => '-',
    description => 'The layer integrated 532 nm volume depolarization ratio (IVDR) computed for each layer. The quality of the estimate for IVDR is determined by the accuracy of the top and base identification, the reliability of the polarization gain ratio calibration, and by the signal-to-noise ratio (SNR) of the backscatter data within the layer.'
  },
  {
    name        => 'L2_333mMLay.Integrated_Volume_Depolarization_Ratio_Uncertainty',
    units       => '-',
    description => 'The uncertainties reported for the 532nm layer-integrated volume depolarization ratios provide an estimate of the total random error in the combined backscatter signals.'
  },
  {
    name        => 'L2_333mMLay.Laser_Energy_532',
    units       => 'J',
    description => 'The laser energy at 532nm measured by the laser energy monitor for each shot.'
  },
  {
    name        => 'L2_333mMLay.Layer_Base_Altitude',
    units       => 'km',
    description => 'Base of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_333mMLay.Layer_Base_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer base altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_333mMLay.Layer_Base_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer attenuated backscatter centroid altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_333mMLay.Layer_Centroid_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer base altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_333mMLay.Layer_IAB_QA_Factor',
    units       => '-',
    description => 'Single layer analog of the column integrated attenuated backscatter cumulative probability.'
  },
  {
    name        => 'L2_333mMLay.Layer_Top_Altitude',
    units       => 'km',
    description => 'Top of the layer above mean sea level. Due to the on-board data averaging scheme, the precision with which CALIPSO can make this measurement is itself a function of altitude.'
  },
  {
    name        => 'L2_333mMLay.Layer_Top_Pressure',
    units       => 'hPa',
    description => 'Pressure at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_333mMLay.Layer_Top_Temperature',
    units       => 'degC',
    description => 'Temperature at the layer top altitude; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_333mMLay.Low_Energy_Mitigation_Column_QC_Flag',
    units       => '-',
    description => 'Low energy mitigation (LEM) quality control flag reported for each column; least significant bit = bit 1 (bit; interpretation): bit 1; Column contains LEM affected data (data has been rejected or low energy shots that LEM accepts, bit 2; Column belongs to a 5km resolution frame that has been rejected due to too many unusable profiles, bit 3; Column belongs to a 5 km resolution frame that has been rejected due to too many rejected subregions in altitude region 3, bit 4; Column belongs to a 5km resolution frame that has been rejected due to too many rejected subregions in altitude region 4, bit 5; Feature detection at 20km resolution not performed in this column due to insufficient fraction of usable 5km resolution frames, bit 6; Feature detection at 80 km resolution not performed in this column due to insufficient fraction of usable 5 km resolution frames, bits 7-16; Unused.'
  },
  {
    name        => 'L2_333mMLay.Midlayer_Pressure',
    units       => 'hPa',
    description => 'Pressure at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_333mMLay.Midlayer_Temperature',
    units       => 'degC',
    description => 'Temperature at the geometric midpoint of the layer in the vertical dimension; derived from the MERRA-2 data product.'
  },
  {
    name        => 'L2_333mMLay.Number_Bins_Shift',
    units       => '-',
    description => 'Number of 30 meter bins shifted by the altitude registration algorithm to match the range bin at mean sea level in the measured profile to the range bin at mean sealevel in the fixed altitude grid; positive = shift up and negative = shift down.'
  },
  {
    name        => 'L2_333mMLay.Number_Layers_Found',
    units       => '-',
    description => 'The number of cloud layers found in a column. The merged 5km layer product reports layers detected at multiple averaging resolutions (5km, 20km, and 80km). Because the reporting resolution (5-km) is not always identical to the detection resolution, layers may appear to overlap in the vertical dimension.'
  },
  {
    name        => 'L2_333mMLay.ODCOD_Effective_Optical_Depth_532',
    units       => '-',
    description => 'Column estimate of apparent (multiple scattering effects not accounted for) optical depth at 532nm wavelength for all particulate (clouds and aerosols) within the illuminated volume from just above the lidar detected ocean surface to the CALIPSO top of atmosphere.'
  },
  {
    name        => 'L2_333mMLay.ODCOD_Effective_Optical_Depth_532_Uncertainty',
    units       => '-',
    description => 'An analytical estimate of the absolute standard deviation uncertainty in the ODCOD Effective Optical Depth estimate at 532nm.'
  },
  {
    name        => 'L2_333mMLay.ODCOD_QC_Flag_532',
    units       => '-',
    description => 'Ocean Derived Column Optical Depth (ODCOD) quality control flags; least significant bit = bit 1 (bit; interpretation): bit 1; The fundamental time delay was shifted from the first point, bit 2; Lidar surface data too many data points for a valid surface spike, bit 3; Range bin above recorded lidar surface data added, bit 4; Range bin below recorded lidar surface data added, bit 5; The surface point ratio indicated that the first point of the surface spike is missing, bit 6; To find the time delay, the surface spike index range had to be shifted down, bit 7; The retrieval inputs fall outside of confidence ranges indicating a confident retrieval, bits 8-10; Unused, none, bit 11; No Lidar surface was found, bit 12; The surface is not ocean, bit 13; The surface is ice, bit 14; The wind speed is outside the valid range; 0.025 to 43 m/s, bit 15; The fundamental time delay found was too low and there were not enough points to find a new one, bit 16; Too few data points to solve, bit 17; The low pass filter impulse response model curve area was too large, bit 18; The low pass filter impulse response model curve scale factor could no t be found, bit 19; Lidar surface saturated, bit 20; Lidar surface has negative signal anomaly, bit 21; Lidar surface data is not valid, bit 22; Flags bad input/ancillary data, bit 23; No Lidar surface was detected via derivative method at greater than single shot resolutions, bits 24-32; Unused, none.'
  },
# Vectors are UNIMPLEMENTED.
#  {
#    name        => 'L2_333mMLay.ODCOD_Surface_Wind_Speed_Correction',
#    units       => 'm/s',
#    description => 'An additive correction derived from AMSR analysis applied to the MERRA-2 wind speed magnitude calculated from the zonal winds reported in ODCOD Surface Wind Speeds 10m.'
#  },
#  {
#    name        => 'L2_333mMLay.ODCOD_Surface_Wind_Speeds_10m',
#    units       => 'm/s',
#    description => 'Mean zonal and meridional component of the surface wind speed 10 meters above the surface provided by MERRA-2. First Dimension = eastward, zonal, u; Second Dimension = northward, meridional, v.'
#  },
  {
    name        => 'L2_333mMLay.Off_Nadir_Angle',
    units       => 'deg',
    description => 'Angle of the lidar viewing vector; 0.3 degrees off-nadir from launch/first-light to 3.0 degrees after November 28, 2007.'
  },
  {
    name        => 'L2_333mMLay.Opacity_Flag',
    units       => '-',
    description => 'Flag that identifies a layer in which the backscatter signal becomes completely attenuated. An opacity flag of 0 indicates a transparent layer, while a value of 1 indicates an opaque layer. A layer is considered opaque if (a) it is the lowest feature detected in a column, and (b) it is not subsequently classified as a surface return. For those layers identified as opaque, the reported base altitude must be considered as an apparent base, rather than a true base.'
  },
  {
    name        => 'L2_333mMLay.Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '1/sr',
    description => 'Overlying_Integrated_Attenuated_Backscatter_532'
  },
  {
    name        => 'L2_333mMLay.Parallel_Column_Reflectance_532',
    units       => '-',
    description => 'Bi-directional column reflectance derived from the root-mean-square variation of the 532 nm parallel channel baseline measurements.'
  },
  {
    name        => 'L2_333mMLay.Parallel_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Column reflectance uncertainty for the 532nm parallel channel.'
  },
  {
    name        => 'L2_333mMLay.Perpendicular_Column_Reflectance_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_532'
  },
  {
    name        => 'L2_333mMLay.Perpendicular_Column_Reflectance_Uncertainty_532',
    units       => '-',
    description => 'Perpendicular_Column_Reflectance_Uncertainty_532'
  },
# unmatched dimensions
#  {
#    name        => 'L2_333mMLay.Profile_ID',
#    units       => '-',
#    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
#  },
#  {
#    name        => 'L2_333mMLay.Profile_Time',
#    units       => 'seconds_since_1993-01-01',
#    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
#  },
  {
    name        => 'L2_333mMLay.Scattering_Angle',
    units       => 'deg',
    description => 'Angle between the lidar viewing angle and the line of sight of the sun.'
  },
  {
    name        => 'L2_333mMLay.Scene_Flag',
    units       => '-',
    description => 'Bit mapped flag where the following bits indicate the atmospheric feature type was detected at some resolution in the colum n; least significant bit = bit 1 (bit; interpretation): bit 1; tropospheric aerosol subtype clean marine, bit 2; tropospheric aerosol subtype dust, bit 3; tropospheric aerosol subtype polluted continental/smoke, bit 4; tropospheric aerosol subtype clean continental, bit 5; tropospheric aerosol subtype polluted dust, bit 6; tropospheric aerosol subtype elevated smoke, bit 7; tropospheric aerosol subtype dusty marine, bit 8; stratospheric aerosol subtype polar stratospheric aerosol, bit 9; stratospheric aerosol subtype volcanic ash, bit 10; stratospheric aerosol subtype sulfate, bit 11; stratospheric aerosol subtype elevated smoke bit 12; stratospheric aerosol unclassified. bit 13; On if column has ROI clouds at any resolution. bit 14; On if column has HOI clouds at any resolution. bit 15; On if column has water clouds at any resolution. bit 16; On if column has unknown phase clouds at any resolution.'
  },
  {
    name        => 'L2_333mMLay.Snow_Ice_Surface_Type',
    units       => '-',
    description => 'Snow and ice coverage for the surface at the laser footprint provided by the U.S. Airforce Weather Agency (AFWA). Values between 0 and 100 represent the percentage of ice, while 101 = permanent ice, 103 = snow, and 255 = mixed pixels at coastline.'
  },
  {
    name        => 'L2_333mMLay.Solar_Azimuth_Angle',
    units       => 'degrees',
    description => 'Azimuth angle from north of the line of sight to the sun.'
  },
  {
    name        => 'L2_333mMLay.Solar_Zenith_Angle',
    units       => 'degrees',
    description => 'Angle between the zenith at the lidar footprint on the surface and the line of sight to the sun.'
  },
  {
    name        => 'L2_333mMLay.Spacecraft_Position_X',
    units       => 'km',
    description => 'X-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_333mMLay.Spacecraft_Position_Y',
    units       => 'km',
    description => 'Y-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_333mMLay.Spacecraft_Position_Z',
    units       => 'km',
    description => 'Z-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_333mMLay.Surface_1064_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 1064 to Surface Base 1064, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_333mMLay.Surface_1064_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 1064 to Surface Base Altitude 1064, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_333mMLay.Surface_532_Integrated_Attenuated_Color_Ratio',
    units       => '-',
    description => 'Attenuated color ratio of surface return from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Attenuated Backscatter 1064 to vertically integrated Total Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_333mMLay.Surface_532_Integrated_Depolarization_Ratio',
    units       => '-',
    description => 'Depolarization ratio of surface return at 532nm from Surface Top Altitude 532 to Surface Base Altitude 532, computed as the ratio of vertically integrated Perpendicular Attenuated Backscatter 532 to vertically integrated Parallel Attenuated Backscatter 532.'
  },
  {
    name        => 'L2_333mMLay.Surface_Base_Altitude_1064',
    units       => 'km',
    description => 'Base altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_333mMLay.Surface_Base_Altitude_532',
    units       => 'km',
    description => 'Base altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_333mMLay.Surface_Detection_Flags_1064',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_333mMLay.Surface_Detection_Flags_532',
    units       => '-',
    description => 'Surface detection flag, stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): bit 1; Surface detected (0=no, 1=yes) bits 2-3; Surface detection method, 0 = derivative test, 1 = multi-shot averaged data test, 2 = single shot surface detection fraction test, 3 = unused, bit 4; 532nm parallel saturated, bit 5; 532nm perpendicular saturated, bit 6; 1064nm saturated, bit 7; 532nm parallel negative signal anomaly, bit 8; 532nm perpendicular negative signal anomaly, bit 9; 1064 nm negative signal anomaly, bit 10; Derivative method failure: Z(min gradient) below Z(max gradient), bit 11; Derivative method failure: vertical extent exceeds limit, bit 12; Derivative method failure: peak signal below threshold, bit 13; Failure: vertical separation between 532 and 1064 surface top altitudes exceeds allowable limit, bit 14; Failure: surface detected at 1064nm only, but overlying color ratio below threshold, bits 15-16; Unused.'
  },
  {
    name        => 'L2_333mMLay.Surface_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 1064nm of the surface return from Surface Top Altitude 1064 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_333mMLay.Surface_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Total Attenuated Backscatter 532nm of the surface return from Surface Top Altitude 532 to Surface Base Altitude 1064.'
  },
  {
    name        => 'L2_333mMLay.Surface_Overlying_Integrated_Attenuated_Backscatter_1064',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 1064 from 40km to one range bin above Surface Top Altitude 1064.'
  },
  {
    name        => 'L2_333mMLay.Surface_Overlying_Integrated_Attenuated_Backscatter_532',
    units       => '-',
    description => 'Vertically integrated Attenuated Backscatter 532 from 40km to one range bin above Surface Top Altitude 532.'
  },
  {
    name        => 'L2_333mMLay.Surface_Peak_Signal_1064',
    units       => '/km/sr',
    description => 'Maximum 1064nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_333mMLay.Surface_Peak_Signal_532',
    units       => '/km/sr',
    description => 'Maximum 532nm attenuated backscatter value of the surface signal.'
  },
  {
    name        => 'L2_333mMLay.Surface_Scaled_RMS_Background_1064',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 1064nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_333mMLay.Surface_Scaled_RMS_Background_532',
    units       => '/km/sr',
    description => 'Background noise estimate computed from RMS baseline noise measurements between 65 and 80km above mean sea level, rescaled to create a pseudo 532nm attenuated backscatter coefficient.'
  },
  {
    name        => 'L2_333mMLay.Surface_Top_Altitude_1064',
    units       => 'km',
    description => 'Top altitude of the surface return at 1064nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_333mMLay.Surface_Top_Altitude_532',
    units       => 'km',
    description => 'Top altitude of the surface return at 532nm at the lidar footprint above local mean sea level.'
  },
  {
    name        => 'L2_333mMLay.Tropopause_Height',
    units       => 'km',
    description => 'Tropopause height. Provided by MERRA-2.'
  },
  {
    name        => 'L2_333mMLay.Tropopause_Temperature',
    units       => 'degC',
    description => 'Tropopause temperature. Provided by MERRA-2.'
  },
  {
    name        => 'L2_333mMLay.Was_Cleared',
    units       => '-',
    description => 'Layers which were cleared by the boundary layer cloud-clearing algorithm; 0=not cleared, 1=cleared.'
  },



  # L2_VFM --------------------------------------------------------------------

  {
    name        => 'L2_VFM.Day_Night_Flag',
    units       => '-',
    description => 'Lighting conditions at an altitude of about 24km above mean sea level: 0 = daytime, 1 = nighttime.'
  },
  {
    name        => 'L2_VFM.Feature_Classification_Flags',
    units       => '-',
    description => 'For each layer detected in the CALIPSO backscatter data a feature classification flags are derived. The complete set of flags are stored as a single 16-bit integer; least significant bit = bit 1 (bits; interpretation): 011 bits 1-3; Feature Type, 0110 = invalid (bad or missing data), 0111 = clear air, 0112 = cloud, 0113 = tropospheric aerosol, 0114 = stratospheric aerosol, 0115 = surface, 0116 = subsurface, 0117 = no signal, bits 4-5; Feature Type QA, 011 = 0 = none, 0111 = low, 0112 = medium, 0113 = high 011 bits 6-7; Ice/Water Phase (cloud), 0110 = unknown/not determined, 0111 = ice, 0112 = water, 0113 = oriented ice crystals, bits 8-9; Ice/Water Phase QA (cloud), 0110 = none, 0111 = low, 0112 = medium, 0113 = high, 011bits 10-12; Feature Subtype (clear-air), 0110 = n/a, 0111 = not searched for features at 80km resolution, 0112 = not searched for features at 20 and 8km resolution, bits 10-12; Feature Subtype (Cloud), 0110 = low overcast, transparent, 0111 = low overcast, opaque, 0112 = transition stratocumulus, 0113 = low, broken cloud.'
  },
  {
    name        => 'L2_VFM.Land_Water_Mask',
    units       => '-',
    description => 'Surface type at the laser footprint provided by the CloudSat science team digital elevation model (DEM): 0 = shallow ocean, 1 land, 2 = coastlines, 3 = shallow inland water, 4 = intermittent water, 5 = deep inland water, 6 = continental ocean, 7 = deep ocean.'
  },
  {
    name        => 'L2_VFM.Minimum_Laser_Energy_532',
    units       => 'J',
    description => 'Minimum laser energy at 532 nm measured within each 80 km along-track data segment (80km = 240 single-shot laser pulses). The 80km distance matches the largest horizontal extent considered in CALIOP standard level 2 data analyses.'
  },
  {
    name        => 'L2_VFM.Profile_ID',
    units       => '-',
    description => 'Unique profile identifier generated sequentially for each single-shot profile record. Profile IDs are unique within each granule but not unique over multiple granules.'
  },
  {
    name        => 'L2_VFM.Profile_Time',
    units       => 'seconds_since_1993-01-01',
    description => 'Seconds since January 1, 1983 International Atomic Time (TAI)'
  },
  {
    name        => 'L2_VFM.Spacecraft_Position_X',
    units       => 'km',
    description => 'X-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_VFM.Spacecraft_Position_Y',
    units       => 'km',
    description => 'Y-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_VFM.Spacecraft_Position_Z',
    units       => 'km',
    description => 'Z-coordinate of position (ECR coordinates) of the satellite.'
  },
  {
    name        => 'L2_VFM.VFM_Feature_Detection_Quality_Flag',
    units       => '-',
    description => 'Quality control flag, used to directly map the actions of the low energy mitigation (LEM) algorithm for each range bin in the VFM (flag - interpretation): 0 - First single-shot profile has low energy or missing data, 1 - Second single-shot profile has low energy or missing data, 2 - Third single-shot profile has low energy or missing data, 3 - Fourth single-shot profile has low energy or missing data, 4 - Fifth single-shot profile has low energy or missing data, 5 - Bin in first single-shot profile is rejected by LEM, 6 - Bin in second single-shot profile is rejected by LEM, 7 - Bin in third single-shot profile is rejected by LEM, 8 - Bin in fourth single-shot profile is rejected by LEM, 9 - Bin in fifth single-shot profile is rejected by LEM, 10 - Contributed to feature detection at resolution 1 (1/3km), 11 - Contributed to feature detection at resolution 2 (1-km), 12 - Contributed to feature detection at resolution 3 (5km), 13 - Contributed to feature detection at resolution 4 (20km), 14 - Contributed to feature detection at resolution 5 (80-km)'
  }


);

###############################################################################

sub get_coverage_metadata_reference { return \@coverage_metadata; }

1; # Required last executed statement of module. Do not add code below here!

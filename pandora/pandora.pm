#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: pandora.pm - Pandora coverage metadata used by pandoraserver.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use pandora; # Use my $coverage_metadata_reference = pandora::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-03-23 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package pandora;

##############################################################################

my @coverage_metadata  = (
  {
    name => 'id',
    units => '-',
    description =>
    'Ground-based Pandora Spectrometer instrument globally unique integer id.'
  },
  {
    name => 'locations.id',
    units => '-',
    description =>
    'Ground-based Pandora Spectrometer instrument globally unique integer id.'
  },
  #'{
  #' name' => 'L2Tot_rnvs0p1_7.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2Tot_rnvs0p1_7.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2Tot_rnvs0p1_7.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2Tot_rnvs0p1_7.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs0p1_7.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs0p1_7.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs0p1_7.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs0p1_7.nitrogen_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Nitrogen dioxide total vertical column amount'
  },
  { name => 'L2Tot_rnvs0p1_7.nitrogen_dioxide_total_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of nitrogen dioxide total vertical column amount'
  },
  { name => 'L2Tot_rnvs0p1_7.direct_nitrogen_dioxide_air_mass_factor',
    units => '-',
    description => 'Direct nitrogen dioxide air mass factor'
  },
  { name => 'L2Tot_rnvs0p1_7.diffuse_correction_for_nitrogen_dioxide',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for nitrogen dioxide, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },
  { name => 'L2Tot_rnvs0p1_7.l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for nitrogen dioxide: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2Tot_rnvs0p1_7.sum_of_nitrogen_dioxide_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Uncertainty too high, 2=Signal to noise ratio too low, 3=Air mass factor too large'
  },
  { name => 'L2Tot_rnvs0p1_7.sum_of_nitrogen_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rnvs0p1_7.fitting_result_index',
    units => '-',
    description => 'Fitting result index: 0,1,2=no error or warning, 3,4=warning, >4=error'
  },
  { name => 'L2Tot_rnvs0p1_7.normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Normalized rms of spectral fitting residuals weighted with measured uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rnvs0p1_7.expected_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on measured uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rnvs0p1_7.expected_normalized_rms_weighted_instrumental_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on instrumental uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rnvs0p1_7.mean_nitrogen_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2Tot_rnvs0p1_7.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2Tot_rnvs0p1_7.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2Tot_rnvs0p1_7.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2Tot_rnvs0p1_7.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2Tot_rnvs0p1_7.l2_fit_data_quality_flag',
    units => '-',
    description => 'Level 2 Fit data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Tot_rnvs0p1_7.sum_of_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Weighted rms of spectral fitting too large, 3=Wavelength shift too large'
  },
  { name => 'L2Tot_rnvs0p1_7.sum_of_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rnvs0p1_7.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Tot_rnvs0p1_7.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=Absolute value of estimated average residual stray light level too high, 6=Although attempted, no wavelength change could be retrieved, 7=Absolute value of retrieved wavelength shift too large, 8=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2Tot_rnvs0p1_7.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rnvs0p1_7.temperature',
    units => 'C',
    description => 'Wavelength effective temperature, 999=no effective temperature given'
  },
  { name => 'L2Tot_rnvs0p1_7.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2Tot_rnvs0p1_7.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2Tot_rnvs0p1_7.wavelength_shift_from_spectral_fitting',
    units => 'nm',
    description => 'Retrieved wavelength shift from spectral fitting, -9=no wavelength change fitting'
  },
  { name => 'L2Tot_rnvs0p1_7.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2Tot_rnvs0p1_7.number_of_dark_count_cycles',
    units => '-',
    description => 'Number of dark count cycles'
  },
  { name => 'L2Tot_rnvs0p1_7.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2Tot_rnvs0p1_7.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },

  #'{
  #' name' => 'L2Tot_rnvs1p1_7.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2Tot_rnvs1p1_7.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2Tot_rnvs1p1_7.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2Tot_rnvs1p1_7.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs1p1_7.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs1p1_7.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs1p1_7.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Tot_rnvs1p1_7.nitrogen_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Nitrogen dioxide total vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2Tot_rnvs1p1_7.nitrogen_dioxide_total_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of nitrogen dioxide total vertical column amount based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no uncertainty could be retrieved'
  },
  { name => 'L2Tot_rnvs1p1_7.direct_nitrogen_dioxide_air_mass_factor',
    units => '-',
    description => 'Direct nitrogen dioxide air mass factor'
  },
  { name => 'L2Tot_rnvs1p1_7.diffuse_correction_for_nitrogen_dioxide',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for nitrogen dioxide, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },
  { name => 'L2Tot_rnvs1p1_7.l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for nitrogen dioxide: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2Tot_rnvs1p1_7.sum_of_nitrogen_dioxide_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Uncertainty too high, 2=Signal to noise ratio too low, 3=Air mass factor too large'
  },
  { name => 'L2Tot_rnvs1p1_7.sum_of_nitrogen_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rnvs1p1_7.fitting_result_index',
    units => '-',
    description => 'Fitting result index: 0,1,2=no error or warning, 3,4=warning, >4=error'
  },
  { name => 'L2Tot_rnvs1p1_7.normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Normalized rms of spectral fitting residuals weighted with measured uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rnvs1p1_7.expected_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on measured uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rnvs1p1_7.expected_normalized_rms_weighted_instrumental_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on instrumental uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rnvs1p1_7.mean_nitrogen_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2Tot_rnvs1p1_7.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2Tot_rnvs1p1_7.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2Tot_rnvs1p1_7.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2Tot_rnvs1p1_7.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2Tot_rnvs1p1_7.l2_fit_data_quality_flag',
    units => '-',
    description => 'Level 2 Fit data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Tot_rnvs1p1_7.sum_of_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Weighted rms of spectral fitting too large, 3=Wavelength shift too large'
  },
  { name => 'L2Tot_rnvs1p1_7.sum_of_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rnvs1p1_7.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Tot_rnvs1p1_7.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=Absolute value of estimated average residual stray light level too high, 6=Although attempted, no wavelength change could be retrieved, 7=Absolute value of retrieved wavelength shift too large, 8=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2Tot_rnvs1p1_7.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rnvs1p1_7.temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2Tot_rnvs1p1_7.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2Tot_rnvs1p1_7.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2Tot_rnvs1p1_7.wavelength_shift_from_spectral_fitting',
    units => 'nm',
    description => 'Retrieved wavelength shift from spectral fitting, -9=no wavelength change fitting'
  },
  { name => 'L2Tot_rnvs1p1_7.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2Tot_rnvs1p1_7.number_of_dark_count_cycles',
    units => '-',
    description => 'Number of dark count cycles'
  },
  { name => 'L2Tot_rnvs1p1_7.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2Tot_rnvs1p1_7.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },

  #'{
  #' name' => 'L2Tot_rout0p1_7.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2Tot_rout0p1_7.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2Tot_rout0p1_7.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2Tot_rout0p1_7.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2Tot_rout0p1_7.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Tot_rout0p1_7.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2Tot_rout0p1_7.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Tot_rout0p1_7.ozone_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Ozone total vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2Tot_rout0p1_7.ozone_total_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of ozone total vertical column amount based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no uncertainty could be retrieved'
  },
  { name => 'L2Tot_rout0p1_7.direct_ozone_air_mass_factor',
    units => '-',
    description => 'Direct ozone air mass factor'
  },
  { name => 'L2Tot_rout0p1_7.diffuse_correction_for_ozone',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for ozone, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },
  { name => 'L2Tot_rout0p1_7.ozone_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for ozone: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_ozone_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for ozone exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Uncertainty too high, 2=Signal to noise ratio too low, 3=Air mass factor too large'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_ozone_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for ozone exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rout0p1_7.sulfur_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Sulfur dioxide total vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2Tot_rout0p1_7.sulfur_dioxide_total_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of sulfur dioxide total vertical column amount based on measured uncertainty, -8=retrieval not successful, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no uncertainty could be retrieved'
  },
  { name => 'L2Tot_rout0p1_7.direct_sulfur_dioxide_air_mass_factor',
    units => '-',
    description => 'Direct sulfur dioxide air mass factor'
  },
  { name => 'L2Tot_rout0p1_7.diffuse_correction_for_sulfur_dioxide',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for sulfur dioxide, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },
  { name => 'L2Tot_rout0p1_7.sulfur_dioxide_l2_quality_flag',
    units => '-',
    description => ' L2 data quality flag for sulfur dioxide: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_sulfur_dioxide_quality_exceeding_dq1',
    units => '-',
    description => ' Sum over 2^i using those i, for which the corresponding L2 data quality parameter for sulfur dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Uncertainty too high, 2=Signal to noise ratio too low, 3=Air mass factor too large'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_sulfur_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for sulfur dioxide exceeds the DQ2 limit (same parameters as for DQ1'
  },
  { name => 'L2Tot_rout0p1_7.fitting_result_index',
    units => '-',
    description => 'Fitting result index: 0,1,2=no error or warning, 3,4=warning, >4=error'
  },
  { name => 'L2Tot_rout0p1_7.normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Normalized rms of spectral fitting residuals weighted with measured uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rout0p1_7.expected_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on measured uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rout0p1_7.expected_normalized_rms_weighted_instrumental_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on instrumental uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2Tot_rout0p1_7.mean_sulfur_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2Tot_rout0p1_7.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2Tot_rout0p1_7.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2Tot_rout0p1_7.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2Tot_rout0p1_7.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2Tot_rout0p1_7.l2_fit_data_quality_flag',
    units => '-',
    description => 'Level 2 Fit data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Weighted rms of spectral fitting too large, 3=Wavelength shift too large'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rout0p1_7.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=Absolute value of estimated average residual stray light level too high, 6=Although attempted, no wavelength change could be retrieved, 7=Absolute value of retrieved wavelength shift too large, 8=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2Tot_rout0p1_7.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Tot_rout0p1_7.temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2Tot_rout0p1_7.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2Tot_rout0p1_7.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2Tot_rout0p1_7.wavelength_shift_from_spectral_fitting',
    units => 'nm',
    description => 'Retrieved wavelength shift from spectral fitting, -9=no wavelength change fitting'
  },
  { name => 'L2Tot_rout0p1_7.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2Tot_rout0p1_7.number_of_dark_count_cycles',
    units => '-',
    description => 'Number of dark count cycles'
  },
  { name => 'L2Tot_rout0p1_7.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2Tot_rout0p1_7.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },

  #'{
  #' name' => 'L2_rfuh5p1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rfuh5p1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rfuh5p1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rfuh5p1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rfuh5p1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rfuh5p1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rfuh5p1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rfuh5p1_8.pointing_zenith_angle',
    units => 'deg',
    description => 'Pointing zenith angle for measurement center'
  },
  { name => 'L2_rfuh5p1_8.pointing_azimuth_angle',
    units => 'deg',
    description => 'Pointing azimuth trig angle for measurement center'
  },
  { name => 'L2_rfuh5p1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rfuh5p1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rfuh5p1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rfuh5p1_8.expected_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rfuh5p1_8.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2_rfuh5p1_8.temperature',
    units => 'C',
    description => 'Climatological station temperature'
  },
  { name => 'L2_rfuh5p1_8.o2_height',
    units => 'm',
    description => 'Climatological effective O2 height'
  },
  { name => 'L2_rfuh5p1_8.o2o2_height',
    units => 'm',
    description => 'Climatological effective O2O2 height'
  },
  { name => 'L2_rfuh5p1_8.surface_o2',
    units => 'mol/m3',
    description => 'Climatological surface O2 concentration'
  },
  { name => 'L2_rfuh5p1_8.surface_o2o2',
    units => 'mol2/m6',
    description => 'Climatological surface O2O2 concentration'
  },
  { name => 'L2_rfuh5p1_8.total_o2_column',
    units => 'molecules/cm2',
    description => 'Climatological total O2 column'
  },
  { name => 'L2_rfuh5p1_8.total_o2o2_column',
    units => 'mol2/m5',
    description => 'Climatological total O2O2 colum'
  },
  { name => 'L2_rfuh5p1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rfuh5p1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rfuh5p1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rfuh5p1_8.mean_formaldehyde_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rfuh5p1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2_rfuh5p1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rfuh5p1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rfuh5p1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift'
  },
  { name => 'L2_rfuh5p1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change'
  },
  { name => 'L2_rfuh5p1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rfuh5p1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rfuh5p1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rfuh5p1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rfuh5p1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rfuh5p1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rfuh5p1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rfuh5p1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rfuh5p1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2 fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rfuh5p1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rfuh5p1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rfuh5p1_8.formaldehyde_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for formaldehyde, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rfuh5p1_8.sum_of_formaldehyde_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for formaldehyde exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rfuh5p1_8.sum_of_formaldehyde_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for formaldehyde exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rfuh5p1_8.surface_formaldehyde',
    units => 'mol/m3',
    description => 'Formaldehyde surface concentration, -9e99=retrieval not successful'
  },
  { name => 'L2_rfuh5p1_8.surface_formaldehyde_uncertainty',
    units => 'mol/m3',
    description => 'Independent uncertainty of formaldehyde surface concentration, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rfuh5p1_8.surface_formaldehyde_index',
    units => '-',
    description => 'Formaldehyde surface concentration index, 1=Fully mixed case from extrapolation to horizon, 2=Fully mixed case from largest pointing zenith angle, 3=Heterogeneous case from extrapolation to horizon, 4=Heterogeneous case from largest pointing zenith angle, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rfuh5p1_8.formaldehyde_heterogeneity_flag',
    units => '-',
    description => 'Formaldehyde heterogeneity flag, 0=well mixed conditions, 1=heterogeneous conditions, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rfuh5p1_8.formaldehyde_tropospheric_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Formaldehyde tropospheric vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2_rfuh5p1_8.formaldehyde_tropospheric_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of formaldehyde tropospheric vertical column amount, -4=tropospheric column was estimated from direct sun measurements using effective temperature, -5=tropospheric column was estimated from direct sun measurements using stratospheric climatology, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rfuh5p1_8.maximum_horizontal_distance_formaldehyde_tropospheric_column',
    units => 'm',
    description => 'Maximum horizontal distance for formaldehyde tropospheric column'
  },
  { name => 'L2_rfuh5p1_8.maximum_vertical_distance_formaldehyde_tropospheric_column',
    units => 'm',
    description => 'Maximum vertical distance for formaldehyde tropospheric column'
  },
  { name => 'L2_rfuh5p1_8.formaldehyde_layer_top_height',
    units => 'm',
    description => 'Top height of formaldehyde layer 1'
  },
  { name => 'L2_rfuh5p1_8.formaldehyde_layer_amount',
    units => 'molecules/cm2',
    description => 'Partial formaldehyde vertical column amount in layer 1'
  },

  #'{
  #' name' => 'L2_rfus5p1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rfus5p1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rfus5p1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rfus5p1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rfus5p1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rfus5p1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rfus5p1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rfus5p1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rfus5p1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rfus5p1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rfus5p1_8.expected_normalized_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rfus5p1_8.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2_rfus5p1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rfus5p1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rfus5p1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rfus5p1_8.mean_formaldehyde_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rfus5p1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2_rfus5p1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rfus5p1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rfus5p1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift, -9=no wavelength change fitting'
  },
  { name => 'L2_rfus5p1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change, -999=no resolution change fitting'
  },
  { name => 'L2_rfus5p1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rfus5p1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rfus5p1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rfus5p1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rfus5p1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rfus5p1_8.aerosol_optical_depth_start',
    units => '-',
    description => 'Estimated aerosol optical depth at starting wavelength of fitting window'
  },
  { name => 'L2_rfus5p1_8.aerosol_optical_depth_center',
    units => '-',
    description => 'Estimated aerosol optical depth at center wavelength of fitting window'
  },
  { name => 'L2_rfus5p1_8.aerosol_optical_depth_end',
    units => '-',
    description => 'Estimated aerosol optical depth at ending wavelength of fitting window'
  },
  { name => 'L2_rfus5p1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rfus5p1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rfus5p1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rfus5p1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2Fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rfus5p1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rfus5p1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rfus5p1_8.fomaldehyde_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for formaldehyde, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rfus5p1_8.sum_of_formaldehyde_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for formaldehyde exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rfus5p1_8.sum_of_formaldehyde_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for formaldehyde exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rfus5p1_8.formaldehyde_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Formaldehyde total vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2_rfus5p1_8.formaldehyde_vertical_column_amount_independent_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of formaldehyde total vertical column amount, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.formaldehyde_vertical_column_amount_structured_uncertainty',
    units => 'molecules/cm2',
    description => 'Structured uncertainty of formaldehyde total vertical column amount, -1=cross section is zero in this wavelength range, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.formaldehyde_vertical_column_amount_common_uncertainty',
    units => 'molecules/cm2',
    description => 'Common uncertainty of formaldehyde total vertical column amount , -1=cross section is zero in this wavelength range, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.formaldehyde_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of formaldehyde total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rfus5p1_8.formaldehyde_vertical_column_amount_rms_uncertainty',
    units => 'molecules/cm2',
    description => 'rms-based uncertainty of formaldehyde total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.temperature',
    units => 'C',
    description => 'Formaldehyde effective temperature'
  },
  { name => 'L2_rfus5p1_8.temperature_independent_uncertainty',
    units => 'C',
    description => 'Independent uncertainty of formaldehyde effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.temperature_structured_uncertainty',
    units => 'C',
    description => 'Structured uncertainty of formaldehyde effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.temperature_common_uncertainty',
    units => 'C',
    description => 'Common uncertainty of formaldehyde effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.temperature_uncertainty',
    units => 'C',
    description => 'Total uncertainty of formaldehyde effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty or common uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rfus5p1_8.direct_formaldehyde_air_mass_factor',
    units => '-',
    description => 'Direct formaldehyde air mass factor'
  },
  { name => 'L2_rfus5p1_8.direct_formaldehyde_air_mass_factor_uncertainty',
    units => '-',
    description => 'Uncertainty of direct formaldehyde air mass factor, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rfus5p1_8.diffuse_correction_for_formaldehyde',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for formaldehyde, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },

  #'{
  #' name' => 'L2_rnvh3p1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rnvh3p1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rnvh3p1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rnvh3p1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rnvh3p1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rnvh3p1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rnvh3p1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rnvh3p1_8.pointing_zenith_angle',
    units => 'deg',
    description => 'Pointing zenith angle for measurement center'
  },
  { name => 'L2_rnvh3p1_8.pointing_azimuth_angle',
    units => 'deg',
    description => 'Pointing azimuth trig angle for measurement center'
  },
  { name => 'L2_rnvh3p1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rnvh3p1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvh3p1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvh3p1_8.expected_normalized_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvh3p1_8.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2_rnvh3p1_8.temperature',
    units => 'C',
    description => 'Climatological station temperature'
  },
  { name => 'L2_rnvh3p1_8.o2_height',
    units => 'm',
    description => 'Climatological effective O2 height'
  },
  { name => 'L2_rnvh3p1_8.o2o2_height',
    units => 'm',
    description => 'Climatological effective O2O2 height'
  },
  { name => 'L2_rnvh3p1_8.surface_o2',
    units => 'mol/m3',
    description => 'Climatological surface O2 concentration'
  },
  { name => 'L2_rnvh3p1_8.surface_o2o2',
    units => 'mol2/m6',
    description => 'Climatological surface O2O2 concentration'
  },
  { name => 'L2_rnvh3p1_8.total_o2_column',
    units => 'molecules/cm2',
    description => 'Climatological total O2 column'
  },
  { name => 'L2_rnvh3p1_8.total_o2o2_column',
    units => 'mol2/m5',
    description => 'Climatological total O2O2 column'
  },
  { name => 'L2_rnvh3p1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rnvh3p1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rnvh3p1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rnvh3p1_8.mean_formaldehyde_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rnvh3p1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2_rnvh3p1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rnvh3p1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rnvh3p1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift, -9=no wavelength change fitting'
  },
  { name => 'L2_rnvh3p1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change, -999=no resolution change fitting'
  },
  { name => 'L2_rnvh3p1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rnvh3p1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rnvh3p1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rnvh3p1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rnvh3p1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rnvh3p1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rnvh3p1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rnvh3p1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvh3p1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2Fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rnvh3p1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rnvh3p1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvh3p1_8.water_vapor_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for water vapor, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rnvh3p1_8.sum_of_water_vapor_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for water vapor exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rnvh3p1_8.sum_of_water_vapor_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for water vapor exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvh3p1_8.surface_water_vapor',
    units => 'mol/m3',
    description => 'Water vapor surface concentration, -9e99=retrieval not successful'
  },
  { name => 'L2_rnvh3p1_8.surface_water_vapor_uncertainty',
    units => 'mol/m3',
    description => 'Independent uncertainty of water vapor surface concentration, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rnvh3p1_8.surface_water_vapor_index',
    units => '-',
    description => 'Water vapor surface concentration index, 1=Fully mixed case from extrapolation to horizon, 2=Fully mixed case from largest pointing zenith angle, 3=Heterogeneous case from extrapolation to horizon, 4=Heterogeneous case from largest pointing zenith angle, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rnvh3p1_8.water_vapor_heterogeneity_flag',
    units => '-',
    description => 'Water vapor heterogeneity flag, 0=well mixed conditions, 1=heterogeneous conditions, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rnvh3p1_8.water_vapor_tropospheric_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Water vapor tropospheric vertical column amount , -9e99=retrieval not successful'
  },
  { name => 'L2_rnvh3p1_8.water_vapor_tropospheric_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of water vapor tropospheric vertical column amount , -4=tropospheric column was estimated from direct sun measurements using effective temperature, -5=tropospheric column was estimated from direct sun measurements using stratospheric climatology, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rnvh3p1_8.maximum_horizontal_distance_water_vapor_tropospheric_column',
    units => 'm',
    description => 'Maximum horizontal distance for water vapor tropospheric column'
  },
  { name => 'L2_rnvh3p1_8.maximum_vertical_distance_water_vapor_tropospheric_column',
    units => 'm',
    description => 'Maximum vertical distance for water vapor tropospheric column'
  },
  { name => 'L2_rnvh3p1_8.nitrogen_dioxide_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for nitrogen dioxide, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rnvh3p1_8.sum_of_nitrogen_dioxide_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rnvh3p1_8.sum_of_nitrogen_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvh3p1_8.surface_nitrogen_dioxide',
    units => 'mol/m3',
    description => 'Nitrogen dioxide surface concentration'
  },
  { name => 'L2_rnvh3p1_8.surface_nitrogen_dioxide_uncertainty',
    units => 'mol/m3',
    description => 'Independent uncertainty of nitrogen dioxide surface concentration, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rnvh3p1_8.surface_nitrogen_dioxide_index',
    units => '-',
    description => 'Nitrogen dioxide surface concentration index, 1=Fully mixed case from extrapolation to horizon, 2=Fully mixed case from largest pointing zenith angle, 3=Heterogeneous case from extrapolation to horizon, 4=Heterogeneous case from largest pointing zenith angle, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rnvh3p1_8.nitrogen_dioxide_heterogeneity_flag',
    units => '-',
    description => 'Nitrogen dioxide heterogeneity flag, 0=well mixed conditions, 1=heterogeneous conditions, -6=no surface concentration was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rnvh3p1_8.stratospheric_nitrogen_dioxide',
    units => 'molecules/cm2',
    description => 'Climatological nitrogen dioxide stratospheric column amount'
  },
  { name => 'L2_rnvh3p1_8.stratospheric_nitrogen_dioxide_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of climatological nitrogen dioxide stratospheric column amount'
  },
  { name => 'L2_rnvh3p1_8.tropospheric_nitrogen_dioxide',
    units => 'molecules/cm2',
    description => 'Nitrogen dioxide tropospheric vertical column amount'
  },
  { name => 'L2_rnvh3p1_8.tropospheric_nitrogen_dioxide_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of nitrogen dioxide tropospheric vertical column amount , -4=tropospheric column was estimated from direct sun measurements using effective temperature, -5=tropospheric column was estimated from direct sun measurements using stratospheric climatology, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rnvh3p1_8.maximum_horizontal_distance_nitrogen_dioxide_tropospheric_column',
    units => 'm',
    description => 'Maximum horizontal distance for nitrogen dioxide tropospheric column'
  },
  { name => 'L2_rnvh3p1_8.maximum_vertical_distance_nitrogen_dioxide_tropospheric_column',
    units => 'm',
    description => 'Maximum vertical distance for nitrogen dioxide tropospheric column'
  },
  { name => 'L2_rnvh3p1_8.water_vapor_layer_top_height',
    units => 'm',
    description => 'Top height of water vapor layer 1, -6=no profile was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rnvh3p1_8.water_vapor_layer_amount',
    units => 'molecules/cm2',
    description => 'Partial water vapor vertical column amount in layer 1 , -9e99=retrieval not successful'
  },
  { name => 'L2_rnvh3p1_8.nitrogen_dioxide_layer_top_height',
    units => 'm',
    description => 'Top height of nitrogen dioxide layer 1, -6=no profile was retrieved since the maximum viewing zenith angle was below 87deg'
  },
  { name => 'L2_rnvh3p1_8.nitrogen_dioxide_layer_amount',
    units => 'molecules/cm2',
    description => 'Partial nitrogen dioxide vertical column amount in layer 1 , -9e99=retrieval not successful'
  },

  #'{
  #' name' => 'L2_rnvs3p1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rnvs3p1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rnvs3p1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rnvs3p1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rnvs3p1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rnvs3p1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rnvs3p1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rnvs3p1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rnvs3p1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvs3p1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvs3p1_8.expected_normalized_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvs3p1_8.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2_rnvs3p1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rnvs3p1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rnvs3p1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rnvs3p1_8.mean_nitrogen_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rnvs3p1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature, 999=no effective temperature given'
  },
  { name => 'L2_rnvs3p1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rnvs3p1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rnvs3p1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift, -9=no wavelength change fitting'
  },
  { name => 'L2_rnvs3p1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change, -999=no resolution change fitting'
  },
  { name => 'L2_rnvs3p1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rnvs3p1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rnvs3p1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rnvs3p1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rnvs3p1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rnvs3p1_8.aerosol_optical_depth_start',
    units => '-',
    description => 'Estimated aerosol optical depth at starting wavelength of fitting window'
  },
  { name => 'L2_rnvs3p1_8.aerosol_optical_depth_center',
    units => '-',
    description => 'Estimated aerosol optical depth at center wavelength of fitting window'
  },
  { name => 'L2_rnvs3p1_8.aerosol_optical_depth_end',
    units => '-',
    description => 'Estimated aerosol optical depth at ending wavelength of fitting window'
  },
  { name => 'L2_rnvs3p1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rnvs3p1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rnvs3p1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvs3p1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2Fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rnvs3p1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rnvs3p1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvs3p1_8.nitrogen_dioxide_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for nitrogen dioxide, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rnvs3p1_8.sum_of_nitrogen_dioxide_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rnvs3p1_8.sum_of_nitrogen_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvs3p1_8.nitrogen_dioxide_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Nitrogen dioxide total vertical column amount'
  },
  { name => 'L2_rnvs3p1_8.nitrogen_dioxide_vertical_column_independent_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvs3p1_8.nitrogen_dioxide_vertical_column_structured_uncertainty',
    units => 'molecules/cm2',
    description => 'Structured uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvs3p1_8.nitrogen_dioxide_vertical_column_common_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rnvs3p1_8.nitrogen_dioxide_vertical_column_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rnvs3p1_8.nitrogen_dioxide_vertical_rms_uncertainty',
    units => 'molecules/cm2',
    description => 'rms-based uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvs3p1_8.temperature',
    units => 'C',
    description => 'Nitrogen dioxide effective temperature'
  },
  { name => 'L2_rnvs3p1_8.temperature_independent_uncertainty',
    units => 'C',
    description => 'Independent uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvs3p1_8.temperature_structured_uncertainty',
    units => 'C',
    description => 'Structured uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvs3p1_8.temperature_common_uncertainty',
    units => 'C',
    description => 'Common uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvs3p1_8.temperature_uncertainty',
    units => 'C',
    description => 'Total uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty or common uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvs3p1_8.direct_nitrogen_dioxide_air_mass_factor',
    units => '-',
    description => 'Direct nitrogen dioxide air mass factor'
  },
  { name => 'L2_rnvs3p1_8.direct_nitrogen_dioxide_air_mass_factor_uncertainty',
    units => '-',
    description => ' Uncertainty of direct nitrogen dioxide air mass factor, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rnvs3p1_8.diffuse_correction',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for nitrogen dioxide, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },
  { name => 'L2_rnvs3p1_8.stratospheric_nitrogen_dioxide',
    units => 'molecules/cm2',
    description => 'Climatological nitrogen dioxide stratospheric column amount'
  },
  { name => 'L2_rnvs3p1_8.stratospheric_nitrogen_dioxide_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of climatological nitrogen dioxide stratospheric column amount'
  },

  #'{
  #' name' => 'L2_rout2p1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rout2p1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rout2p1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rout2p1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rout2p1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rout2p1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rout2p1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rout2p1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rout2p1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rout2p1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rout2p1_8.expected_normalized_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rout2p1_8.pressure',
    units => 'hPa',
    description => ' Climatological station pressure'
  },
  { name => 'L2_rout2p1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rout2p1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rout2p1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rout2p1_8.mean_ozone_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rout2p1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2_rout2p1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rout2p1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rout2p1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift, -9=no wavelength change fitting'
  },
  { name => 'L2_rout2p1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change, -999=no resolution change fitting'
  },
  { name => 'L2_rout2p1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rout2p1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rout2p1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rout2p1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rout2p1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rout2p1_8.aerosol_optical_depth_start',
    units => '-',
    description => 'Estimated aerosol optical depth at starting wavelength of fitting window'
  },
  { name => 'L2_rout2p1_8.aerosol_optical_depth_center',
    units => '-',
    description => 'stimated aerosol optical depth at center wavelength of fitting window'
  },
  { name => 'L2_rout2p1_8.aerosol_optical_depth_end',
    units => '-',
    description => 'Estimated aerosol optical depth at ending wavelength of fitting window'
  },
  { name => 'L2_rout2p1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rout2p1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rout2p1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rout2p1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2Fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rout2p1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rout2p1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rout2p1_8.ozone_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for ozone, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rout2p1_8.sum_of_ozone_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for ozone exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rout2p1_8.sum_of_ozone_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for ozone exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rout2p1_8.ozone_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Ozone total vertical column amount , -9e99=retrieval not successful'
  },
  { name => 'L2_rout2p1_8.ozone_vertical_column_independent_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of ozone total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rout2p1_8.ozone_vertical_column_structured_uncertainty',
    units => 'molecules/cm2',
    description => 'Structured uncertainty of ozone total vertical column amount , -1=cross section is zero in this wavelength range, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rout2p1_8.ozone_vertical_column_common_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of ozone total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rout2p1_8.ozone_vertical_column_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of ozone total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rout2p1_8.ozone_vertical_column_rms_uncertainty',
    units => 'molecules/cm2',
    description => 'rms-based uncertainty of ozone total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -9=spectral fitting was not successful'
  },
  { name => 'L2_rout2p1_8.temperature',
    units => 'C',
    description => 'Ozone effective temperature'
  },
  { name => 'L2_rout2p1_8.temperature_independent_uncertainty',
    units => 'C',
    description => 'Independent uncertainty of ozone effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rout2p1_8.temperature_structured_uncertainty',
    units => 'C',
    description => 'Structured uncertainty of ozone effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -9=spectral fitting was not successful'
  },
  { name => 'L2_rout2p1_8.temperature_common_uncertainty',
    units => 'C',
    description => 'Common uncertainty of ozone effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rout2p1_8.temperature_uncertainty',
    units => 'C',
    description => 'Total uncertainty of ozone effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty or common uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rout2p1_8.direct_ozone_air_mass_factor',
    units => '-',
    description => 'Direct ozone air mass factor'
  },
  { name => 'L2_rout2p1_8.ozone_air_mass_factor_uncertainty',
    units => '-',
    description => 'Uncertainty of direct ozone air mass factor, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rout2p1_8.diffuse_correction',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for ozone, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },

  #'{
  #' name' => 'L2_rsus1p1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rsus1p1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rsus1p1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rsus1p1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rsus1p1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rsus1p1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rsus1p1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rsus1p1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rsus1p1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rsus1p1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rsus1p1_8.expected_normalized_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rsus1p1_8.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2_rsus1p1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rsus1p1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rsus1p1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rsus1p1_8.mean_sulfur_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rsus1p1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature, 999=no effective temperature given'
  },
  { name => 'L2_rsus1p1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rsus1p1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rsus1p1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift, -9=no wavelength change fitting'
  },
  { name => 'L2_rsus1p1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change, -999=no resolution change fitting'
  },
  { name => 'L2_rsus1p1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rsus1p1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rsus1p1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rsus1p1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rsus1p1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rsus1p1_8.aerosol_optical_depth_start',
    units => '-',
    description => 'Estimated aerosol optical depth at starting wavelength of fitting window'
  },
  { name => 'L2_rsus1p1_8.aerosol_optical_depth_center',
    units => '-',
    description => 'Estimated aerosol optical depth at center wavelength of fitting window'
  },
  { name => 'L2_rsus1p1_8.aerosol_optical_depth_end',
    units => '-',
    description => 'Estimated aerosol optical depth at ending wavelength of fitting window'
  },
  { name => 'L2_rsus1p1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rsus1p1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rsus1p1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rsus1p1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2Fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rsus1p1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rsus1p1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for sulfur dioxide, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rsus1p1_8.sum_of_sulfur_dioxide_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for sulfur dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rsus1p1_8.sum_of_sulfur_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for sulfur dioxide exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Sulfur dioxide total vertical column amount , -9e99=retrieval not successful'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_vertical_column_independent_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of sulfur dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_vertical_column_structured_uncertainty',
    units => 'molecules/cm2',
    description => 'Structured uncertainty of sulfur dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_vertical_column_common_uncertainty',
    units => 'molecules/cm2',
    description => 'Common uncertainty of sulfur dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_vertical_column_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of sulfur dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_vertical_column_rms_uncertainty',
    units => 'molecules/cm2',
    description => 'rms-based uncertainty of sulfur dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.temperature',
    units => 'C',
    description => 'Sulfur dioxide effective temperature'
  },
  { name => 'L2_rsus1p1_8.temperature_independent_uncertainty',
    units => 'C',
    description => 'Independent uncertainty of sulfur dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.temperature_structured_uncertainty',
    units => 'C',
    description => 'Structured uncertainty of sulfur dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.temperature_common_uncertainty',
    units => 'C',
    description => 'Common uncertainty of sulfur dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.temperature_uncertainty',
    units => 'C',
    description => 'Total uncertainty of sulfur dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty or common uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rsus1p1_8.direct_sulfur_dioxide_air_mass_factor',
    units => '-',
    description => 'Direct sulfur dioxide air mass factor'
  },
  { name => 'L2_rsus1p1_8.sulfur_dioxide_air_mass_factor_uncertainty',
    units => '-',
    description => 'Uncertainty of direct sulfur dioxide air mass factor, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rsus1p1_8.diffuse_correction',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for sulfur dioxide, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },

  #'{
  #' name' => 'L2Trop_rnvh1p1_7.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2Trop_rnvh1p1_7.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2Trop_rnvh1p1_7.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2Trop_rnvh1p1_7.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2Trop_rnvh1p1_7.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Trop_rnvh1p1_7.pointing_zenith_angle',
    units => 'deg',
    description => 'Pointing zenith angle for center of measurement in degree'
  },
  { name => 'L2Trop_rnvh1p1_7.pointing_azimuth_angle',
    units => 'deg',
    description => 'Pointing azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2Trop_rnvh1p1_7.water_vapor_surface',
    units => 'ppb',
    description => 'Water vapor surface concentration, -9e99=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.water_vapor_surface_uncertainty',
    units => 'ppb',
    description => 'Uncertainty of water vapor surface concentration, -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.water_vapor_tropospheric',
    units => 'cm',
    description => 'Water vapor tropospheric vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.water_vapor_tropospheric_uncertainty',
    units => 'cm',
    description => 'Uncertainty of water vapor tropospheric vertical column amount, -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.water_vapor_surface_index',
    units => '-',
    description => 'Water vapor surface concentration index, integer >0 (see manual), -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.water_vapor_heterogeneity_flag',
    units => '-',
    description => 'Water vapor heterogeneity flag, 0=well mixed conditions, 1=heterogeneous conditions, -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.water_vapor_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for water vapor: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_water_vapor_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for water vapor exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Uncertainty too high, 2=Signal to noise ratio too low, 3=Air mass factor too large'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_water_vapor_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for water vapor exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Trop_rnvh1p1_7.direct_water_vapor_air_mass_factor',
    units => '-',
    description => 'Direct water vapor air mass factor'
  },
  { name => 'L2Trop_rnvh1p1_7.nitrogen_dioxide_surface',
    units => 'ppb',
    description => 'Nitrogen dioxide surface concentration, -9e99=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.nitrogen_dioxide_surface_uncertainty',
    units => 'ppb',
    description => 'Uncertainty of nitrogen dioxide surface concentration, -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.nitrogen_dioxide_tropospheric',
    units => 'molecules/cm2',
    description => 'Nitrogen dioxide tropospheric vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.nitrogen_dioxide_tropospheric_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of nitrogen dioxide tropospheric vertical column amount, -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.nitrogen_dioxide_surface_index',
    units => '-',
    description => 'Nitrogen dioxide surface concentration index, integer >0 (see manual), -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.nitrogen_dioxide_heterogeneity_flag',
    units => '-',
    description => 'Nitrogen dioxide heterogeneity flag, 0=well mixed conditions, 1=heterogeneous conditions, -8=retrieval not successful'
  },
  { name => 'L2Trop_rnvh1p1_7.nitrogen_dioxide_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for nitrogen dioxide: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_nitrogen_dioxide_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Uncertainty too high, 2=Signal to noise ratio too low, 3=Air mass factor too large'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_nitrogen_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2Trop_rnvh1p1_7.direct_nitrogen_dioxide_air_mass_factor',
    units => '-',
    description => 'Direct nitrogen dioxide air mass factor'
  },
  { name => 'L2Trop_rnvh1p1_7.fitting_result_index',
    units => '-',
    description => 'Maximum fitting result index: 0,1,2=no error or warning, 3,4=warning, >4=error'
  },
  { name => 'L2Trop_rnvh1p1_7.mean_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Normalized rms of spectral fitting residuals weighted with measured uncertainty, mean over all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Normalized rms of spectral fitting residuals weighted with measured uncertainty, minimum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Normalized rms of spectral fitting residuals weighted with measured uncertainty, maximum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.mean_expected_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on measured uncertainty, mean over all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_expected_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on measured uncertainty, minimum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_expected_normalized_rms_weighted_measured_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on measured uncertainty, maximum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.mean_expected_normalized_rms_weighted_instrumental_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on instrumental uncertainty, mean over all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_expected_normalized_rms_weighted_instrumental_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on instrumental uncertainty, minimum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_expected_normalized_rms_weighted_instrumental_uncertainty',
    units => '-',
    description => 'Expected normalized rms of weighted spectral fitting residuals based on instrumental uncertainty, maximum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2Trop_rnvh1p1_7.temperature',
    units => 'C',
    description => 'Climatological station temperature'
  },
  { name => 'L2Trop_rnvh1p1_7.o2o2_height',
    units => 'm',
    description => 'Climatological effective O2O2 height'
  },
  { name => 'L2Trop_rnvh1p1_7.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2Trop_rnvh1p1_7.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2Trop_rnvh1p1_7.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_l2_fit_data_quality_flag',
    units => '-',
    description => 'Level 2 Fit data quality flag, minimum of all sets: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_l2_fit_data_quality_flag',
    units => '-',
    description => 'Level 2 Fit data quality flag, maximum of all sets: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, maximum of all sets, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Weighted rms of spectral fitting too large, 3=Wavelength shift too large'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit, maximum of all sets (same parameters as for DQ1)'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag, minimum of all sets: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag, maximum of all sets: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, maximum of all sets, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=Absolute value of estimated average residual stray light level too high, 6=Although attempted, no wavelength change could be retrieved, 7=Absolute value of retrieved wavelength shift too large, 8=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2Trop_rnvh1p1_7.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit, maximum of all sets (same parameters as for DQ1)'
  },
  { name => 'L2Trop_rnvh1p1_7.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature, mean over all sets, 999=no effective temperature given'
  },
  { name => 'L2Trop_rnvh1p1_7.mean_residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level, mean over all sets (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level, minimum of all sets (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level, maximum of all sets (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2Trop_rnvh1p1_7.mean_wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, mean over all sets, -9=no wavelength change determination'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, minimum of all sets, -9=no wavelength change determination'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, maximum of all sets, -9=no wavelength change determination'
  },
  { name => 'L2Trop_rnvh1p1_7.mean_wavelength_shift_from_spectral_fitting',
    units => 'nm',
    description => 'Retrieved wavelength shift from spectral fitting, mean over all sets, -9=no wavelength change fitting'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_wavelength_shift_from_spectral_fitting',
    units => 'nm',
    description => 'Retrieved wavelength shift from spectral fitting, minimum of all sets, -9=no wavelength change fitting'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_wavelength_shift_from_spectral_fitting',
    units => 'nm',
    description => 'Retrieved wavelength shift from spectral fitting, maximum of all sets, -9=no wavelength change fitting'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_integration_time',
    units => 'ms',
    description => 'Integration time, minimum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_integration_time',
    units => 'ms',
    description => 'Integration time, maximum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.mean_number_of_dark_count_cycles',
    units => '-',
    description => 'Number of dark count cycles, mean over all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_number_of_dark_count_cycles',
    units => '-',
    description => 'Number of dark count cycles, minimum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_number_of_dark_count_cycles',
    units => '-',
    description => 'Number of dark count cycles, maximum of all sets'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, minimum of all sets, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, maximum of all sets, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2Trop_rnvh1p1_7.minimum_filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, minimum of all sets, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2Trop_rnvh1p1_7.maximum_filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, maximum of all sets, 0=filterwheel not used, 1-9 are valid positions'
  },

  #'{
  #' name' => 'L2_rnvssp1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rnvssp1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rnvssp1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rnvssp1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rnvssp1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rnvssp1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rnvssp1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rnvssp1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rnvssp1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvssp1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvssp1_8.expected_normalized_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rnvssp1_8.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2_rnvssp1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rnvssp1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rnvssp1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rnvssp1_8.mean_nitrogen_dioxide_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rnvssp1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2_rnvssp1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rnvssp1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rnvssp1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift, -9=no wavelength change fitting'
  },
  { name => 'L2_rnvssp1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change, -999=no resolution change fitting'
  },
  { name => 'L2_rnvssp1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rnvssp1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rnvssp1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rnvssp1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rnvssp1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rnvssp1_8.aerosol_optical_depth_start',
    units => '-',
    description => 'Estimated aerosol optical depth at starting wavelength of fitting window'
  },
  { name => 'L2_rnvssp1_8.aerosol_optical_depth_center',
    units => '-',
    description => 'Estimated aerosol optical depth at center wavelength of fitting window'
  },
  { name => 'L2_rnvssp1_8.aerosol_optical_depth_end',
    units => '-',
    description => 'Estimated aerosol optical depth at ending wavelength of fitting window'
  },
  { name => 'L2_rnvssp1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rnvssp1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => ' Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rnvssp1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvssp1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2Fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rnvssp1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rnvssp1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvssp1_8.nitrogen_dioxide_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for nitrogen dioxide, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rnvssp1_8.sum_of_nitrogen_dioxide_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rnvssp1_8.sum_of_nitrogen_dioxide_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for nitrogen dioxide exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rnvssp1_8.nitrogen_dioxide_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Nitrogen dioxide total vertical column amount , -9e99=retrieval not successful'
  },
  { name => 'L2_rnvssp1_8.nitrogen_dioxide_vertical_column_independent_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.nitrogen_dioxide_vertical_column_structured_uncertainty',
    units => 'molecules/cm2',
    description => 'Structured uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.nitrogen_dioxide_vertical_column_common_uncertainty',
    units => 'molecules/cm2',
    description => ' Common uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.nitrogen_dioxide_vertical_column_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rnvssp1_8.nitrogen_dioxide_vertical_rms_uncertainty',
    units => 'molecules/cm2',
    description => 'rms-based uncertainty of nitrogen dioxide total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.temperature',
    units => 'C',
    description => 'Nitrogen dioxide effective temperature'
  },
  { name => 'L2_rnvssp1_8.temperature_independent_uncertainty',
    units => 'C',
    description => 'Independent uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.temperature_structured_uncertainty',
    units => 'C',
    description => 'Structured uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.temperature_common_uncertainty',
    units => 'C',
    description => 'Common uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.temperature_uncertainty',
    units => 'C',
    description => 'Total uncertainty of nitrogen dioxide effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty or common uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rnvssp1_8.direct_nitrogen_dioxide_air_mass_factor',
    units => '-',
    description => 'Direct nitrogen dioxide air mass factor'
  },
  { name => 'L2_rnvssp1_8.direct_nitrogen_dioxide_air_mass_factor_uncertainty',
    units => '-',
    description => 'Uncertainty of direct nitrogen dioxide air mass factor, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rnvssp1_8.diffuse_correction',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for nitrogen dioxide, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  },
  { name => 'L2_rnvssp1_8.stratospheric_nitrogen_dioxide',
    units => 'molecules/cm2',
    description => 'Climatological nitrogen dioxide stratospheric column amount'
  },
  { name => 'L2_rnvssp1_8.stratospheric_nitrogen_dioxide_uncertainty',
    units => 'molecules/cm2',
    description => 'Uncertainty of climatological nitrogen dioxide stratospheric column amount'
  },

  #'{
  #' name' => 'L2_rwvt1p1_8.timestamp',
  # units => 'yyyymmddhhmmss',
  # description => 'UT date and time for center of measurement, yyyymmddThhmmssZ (ISO 8601)'
  #},
  { name => 'L2_rwvt1p1_8.days_since_2000',
    units => 'days',
    description => 'Fractional days since 1-Jan-2000 UT midnight for center of measurement'
  },
  { name => 'L2_rwvt1p1_8.measurement_duration',
    units => 's',
    description => 'Effective duration of measurement in seconds'
  },
  { name => 'L2_rwvt1p1_8.solar_zenith_angle',
    units => 'deg',
    description => 'Solar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rwvt1p1_8.solar_azimuth_angle',
    units => 'deg',
    description => 'Solar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rwvt1p1_8.lunar_zenith_angle',
    units => 'deg',
    description => 'Lunar zenith angle for center of measurement in degree'
  },
  { name => 'L2_rwvt1p1_8.lunar_azimuth_angle',
    units => 'deg',
    description => 'Lunar azimuth trig angle for center of measurement in degree'
  },
  { name => 'L2_rwvt1p1_8.rms_fitting_residuals',
    units => '-',
    description => 'rms of unweighted fitting residuals, -9=fitting not successful'
  },
  { name => 'L2_rwvt1p1_8.normalized_rms_fitting_residuals',
    units => '-',
    description => 'Normalized rms of fitting residuals weighted with independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rwvt1p1_8.expected_rms_unweighted_fitting_residuals',
    units => '-',
    description => 'Expected rms of unweighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rwvt1p1_8.expected_normalized_rms_weighted_fitting_residuals',
    units => '-',
    description => 'Expected normalized rms of weighted fitting residuals based on independent uncertainty, -9=fitting not successful or no uncertainty given'
  },
  { name => 'L2_rwvt1p1_8.pressure',
    units => 'hPa',
    description => 'Climatological station pressure'
  },
  { name => 'L2_rwvt1p1_8.data_processing_type_index',
    units => '-',
    description => 'Data processing type index'
  },
  { name => 'L2_rwvt1p1_8.calibration_file_version',
    units => '-',
    description => 'Calibration file version'
  },
  { name => 'L2_rwvt1p1_8.calibration_file_date',
    units => 'yyyymmdd',
    description => 'Calibration file validity starting date'
  },
  { name => 'L2_rwvt1p1_8.mean_water_vapor_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Mean value of measured data inside fitting window'
  },
  { name => 'L2_rwvt1p1_8.wavelength_temperature',
    units => 'C',
    description => 'Wavelength effective temperature'
  },
  { name => 'L2_rwvt1p1_8.residual_stray_light_level',
    units => '%',
    description => 'Estimated average residual stray light level (only valid for stray light correction methods 2 and higher)'
  },
  { name => 'L2_rwvt1p1_8.wavelength_shift_from_l1',
    units => 'nm',
    description => 'Retrieved wavelength shift from level 1 data, -9=no wavelength change determination'
  },
  { name => 'L2_rwvt1p1_8.total_wavelength_shift',
    units => 'nm',
    description => 'Retrieved total wavelength shift, -9=no wavelength change fitting'
  },
  { name => 'L2_rwvt1p1_8.resolution_change',
    units => '%',
    description => 'Retrieved resolution change, -999=no resolution change fitting'
  },
  { name => 'L2_rwvt1p1_8.integration_time',
    units => 'ms',
    description => 'Integration time'
  },
  { name => 'L2_rwvt1p1_8.number_of_bright_count_cycles',
    units => '-',
    description => 'Number of bright count cycles'
  },
  { name => 'L2_rwvt1p1_8.filterwheel1_position',
    units => '-',
    description => 'Effective position of filterwheel #1, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rwvt1p1_8.filterwheel2_position',
    units => '-',
    description => 'Effective position of filterwheel #2, 0=filterwheel not used, 1-9 are valid positions'
  },
  { name => 'L2_rwvt1p1_8.atmospheric_variability',
    units => '%',
    description => 'Atmospheric variability, 999=no atmospheric variability was determined'
  },
  { name => 'L2_rwvt1p1_8.aerosol_optical_depth_start',
    units => '-',
    description => 'Estimated aerosol optical depth at starting wavelength of fitting window'
  },
  { name => 'L2_rwvt1p1_8.aerosol_optical_depth_center',
    units => '-',
    description => 'Estimated aerosol optical depth at center wavelength of fitting window'
  },
  { name => 'L2_rwvt1p1_8.aerosol_optical_depth_end',
    units => '-',
    description => 'Estimated aerosol optical depth at ending wavelength of fitting window'
  },
  { name => 'L2_rwvt1p1_8.l1_data_quality_flag',
    units => '-',
    description => 'Level 1 data quality flag: 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rwvt1p1_8.sum_of_l1_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ1 limit, 0=Saturated data, 1=Too few dark counts measurements, 2=No temperature given or effective temperature too different from the reference temperature, 3=Dark count too high, 4=Unsuccessful dark background fitting, 5=The dark count differs significantly from the dark map for too many pixels, 6=Absolute value of estimated average residual stray light level too high, 7=Although attempted, no wavelength change could be retrieved, 8=Absolute value of retrieved wavelength shift too large, 9=Retrieved wavelength shift differs too much from the shift predicted by the effective temperature'
  },
  { name => 'L2_rwvt1p1_8.sum_of_l1_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L1 data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rwvt1p1_8.l2_fit_quality_flag',
    units => '-',
    description => 'L2Fit data quality flag, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality'
  },
  { name => 'L2_rwvt1p1_8.sum_of_l2_fit_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ1 limit, 0=L1 data quality above 0, 1=Spectral fitting was not successful, 2=Wavelength shift too large, 3=Normalized rms of fitting residuals weighted with independent uncertainty too large'
  },
  { name => 'L2_rwvt1p1_8.sum_of_l2_fit_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2Fit data quality parameter exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rwvt1p1_8.fomaldehyde_l2_quality_flag',
    units => '-',
    description => 'L2 data quality flag for water_vapor, 0=assured high quality, 1=assured medium quality, 2=assured low quality, 10=not-assured high quality, 11=not-assured medium quality, 12=not-assured low quality, 20=unusable high quality, 21=unusable medium quality, 22=unusable low quality'
  },
  { name => 'L2_rwvt1p1_8.sum_of_water_vapor_quality_exceeding_dq1',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for water vapor exceeds the DQ1 limit, 0=L2Fit data quality above 0, 1=Retrieval error, 2=Air mass factor too large, 3=Atmospheric variability too large'
  },
  { name => 'L2_rwvt1p1_8.sum_of_water_vapor_quality_exceeding_dq2',
    units => '-',
    description => 'Sum over 2^i using those i, for which the corresponding L2 data quality parameter for water vapor exceeds the DQ2 limit (same parameters as for DQ1)'
  },
  { name => 'L2_rwvt1p1_8.water_vapor_total_vertical_column_amount',
    units => 'molecules/cm2',
    description => 'Water vapor total vertical column amount, -9e99=retrieval not successful'
  },
  { name => 'L2_rwvt1p1_8.water_vapor_vertical_column_amount_independent_uncertainty',
    units => 'molecules/cm2',
    description => 'Independent uncertainty of water vapor total vertical column amount, -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.water_vapor_vertical_column_amount_structured_uncertainty',
    units => 'molecules/cm2',
    description => 'Structured uncertainty of water vapor total vertical column amount, -1=cross section is zero in this wavelength range, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.water_vapor_vertical_column_amount_common_uncertainty',
    units => 'molecules/cm2',
    description => 'Common uncertainty of water vapor total vertical column amount , -1=cross section is zero in this wavelength range, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.water_vapor_vertical_column_amount_uncertainty',
    units => 'molecules/cm2',
    description => 'Total uncertainty of water vapor total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -5=no independent uncertainty input was given, -6=no common uncertainty input was given, -7=not given since method "MEAS" was chosen, -8=not given, since not all components are given'
  },
  { name => 'L2_rwvt1p1_8.water_vapor_vertical_column_amount_rms_uncertainty',
    units => 'molecules/cm2',
    description => 'rms-based uncertainty of water vapor total vertical column amount , -1=cross section is zero in this wavelength range, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.temperature',
    units => 'C',
    description => 'Water vapor effective temperature'
  },
  { name => 'L2_rwvt1p1_8.temperature_independent_uncertainty',
    units => 'C',
    description => 'Independent uncertainty of water vapor effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.temperature_structured_uncertainty',
    units => 'C',
    description => 'Structured uncertainty of water vapor effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.temperature_common_uncertainty',
    units => 'C',
    description => 'Common uncertainty of water vapor effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -6=no independent uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.temperature_uncertainty',
    units => 'C',
    description => 'Total uncertainty of water vapor effective temperature , -1=temperature fitting was requested, but cross section is zero in this wavelength range, -2=no temperature fitting was requested and output for effective temperature and common uncertainty of it is based on f-code, -3=spectral fitting was done, but no independent uncertainty could be retrieved, -4=temperature fitting was requested, but differential optical depth is too small to retrieve temperature, -5=spectral fitting was done without using the independent uncertainty input, -6=no independent uncertainty or common uncertainty input was given, -9=spectral fitting was not successful'
  },
  { name => 'L2_rwvt1p1_8.direct_water_vapor_air_mass_factor',
    units => '-',
    description => 'Direct water vapor air mass factor'
  },
  { name => 'L2_rwvt1p1_8.direct_water_vapor_air_mass_factor_uncertainty',
    units => '-',
    description => 'Uncertainty of direct water vapor air mass factor, -7=uncertainty could not be retrieved since slant column uncertainties were missing'
  },
  { name => 'L2_rwvt1p1_8.diffuse_correction_for_water_vapor',
    units => '%',
    description => 'Diffuse correction applied before fitting at effective fitting wavelength for water_vapor, 0=no diffuse correction applied or fitting not successful, >0=measured diffuse correction, <0=(negative value of) calculated diffuse correction'
  }
);


###############################################################################

sub get_coverage_metadata_reference { return \@coverage_metadata; }

1; # Required last executed statement of module. Do not add code below here!

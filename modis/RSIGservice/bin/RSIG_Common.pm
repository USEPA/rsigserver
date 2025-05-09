package RSIG_Common;
use strict;

=head1 SYNOPSIS

use RSIG_Common;
RSIG_Common::check_cgi(new CGI);
  or
RSIG_Common::check_args(@ARGV);

=head1 DESCRIPTION

Define subroutines that are common to the EPA front and back-ends for
processing EPA WCS queries.

=head1 CONFIG

The following keys must be defined when configuring an instance:

FileList                a reference to an array of hash of at least
                        {FileId,FileType,FileName,DiskId}. This is the list
                        of files to post-process.

=head1 DESIGN NOTES

This module is used by both the front-end and the back-end of the EPA WCS
service hosted by MODAPS Web services. The front-end taint checks and
validates args coming from CGI, and makes a call to SLURM, which schedules
and runs the back-end to actually produce the results. Results get streamed
back to the requestor via stdout/http.

The WCS queries that the EPA runs are augmented with custom features. These
queries do NOT work with standard WCS services as defined by OGC.

=head1 VERSION

 1.0

=head1 AUTHORS AND MAINTAINERS

Greg Ederer

based on the file "modisserver"
from Todd Plessel, Plessel.Todd@epa.gov, 919-541-5500

=head1 ACKNOWLEDGEMENTS

This software is developed by the MODAPS Team for the National
Aeronautics and Space Administrationn, Goddard Space Flight
Center, under contract NAS5-32373.

=cut

###########################################################################

use FileHandle;

# Start figuring out the runtime environment.

# SLURM job id and host we are running on...
$RSIG_Common::JOB_ID = $ENV{SLURM_JOB_ID};
$RSIG_Common::JOB_ID = time() unless $RSIG_Common::JOB_ID;
$RSIG_Common::JOB_ID = $1 if "$RSIG_Common::JOB_ID" =~ /^(\d+)$/;
die "bad JOB_ID" unless $RSIG_Common::JOB_ID;

$ENV{PATH} = undef;  # clean for taint mode
$RSIG_Common::HOST = $1
  if `/usr/bin/hostname` =~ /^([\w\d_-]+)$/;          # taint check
die "HOST not defined in environment" unless $RSIG_Common::HOST;

# Figure out what to set MODAPSROOT to.
# Try environment variable...
$RSIG_Common::MODAPSROOT = $ENV{MODAPSROOT};
if ("$RSIG_Common::MODAPSROOT" =~ /^([\w\d\/]+)$/)
{
  $RSIG_Common::MODAPSROOT = $1;      # taint check
}
die "MODAPSROOT not defined in environment" unless $RSIG_Common::MODAPSROOT;

$RSIG_Common::INSTANCE = $ENV{MODAPS_INSTANCE};
if ("$ENV{MODAPS_INSTANCE}" =~ /^([\w\d]+)$/)
{
  $RSIG_Common::INSTANCE = $1;      # taint check
}
die "MODAPS_INSTANCE not defined in environment" unless $RSIG_Common::INSTANCE;

#####################################################################
# RSIG uses lots of file scope variables.
# Declare them up here so that the subs can find them
#####################################################################
# Server where this program is installed:
my $GRANULE_LIMIT = 10000;
my $debugging = 0;
my $server_path = 'http://modwebsrv.modaps.eosdis.nasa.gov/cgi-bin';
my $service_origin = "$server_path/RSIGservice";

my $subsetter  = "$RSIG_Common::MODAPSROOT/RSIG/bin/MODISSubset";
my $xdrconvert = "$RSIG_Common::MODAPSROOT/RSIG/bin/XDRConvert";
my $compressor = "$RSIG_Common::MODAPSROOT/RSIG/bin/gzip -c";
my $getfiles   = "$RSIG_Common::MODAPSROOT/RSIG/bin/get_modis_files";
# changed in modwebsrv (see above), but might still be used by EPA
#my $subsetter  = "$RSIG_Common::MODAPSROOT/RSIG/bin/Linux.x86_64/MODISSubset";
#my $xdrconvert = "$RSIG_Common::MODAPSROOT/RSIG/bin/Linux.x86_64/XDRConvert";
#my $compressor = "$RSIG_Common::MODAPSROOT/RSIG/bin/Linux.x86_64/gzip -c";
#my $getfiles   = "$RSIG_Common::MODAPSROOT/RSIG/bin/Linux.x86_64/get_modis_files";

# These are output when REQUEST=GetMetadata. See routine print_metadata.
# Always print this first:

my $metadata_content = '
NASA MODIS satellite (Aqua+Terra) measured data accessed using RSIG.
MODIS:
https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/modis/
https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD04_L2/
https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD04_3K/
https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD06_L2/
https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MOD07_L2/
RSIG: https://www.epa.gov/rsig
';

# If FORMAT=original then print this second:

my $metadata_content_original = '
The original data files are shown below.
';

# Else FORMAT!=original then print this second:

my $metadata_content_processed = '
Data processing was done using the RSIG programs MODISSubset and XDRConvert.
MODISSubset is used to subset (by variable, lon-lat box and time range)
and reorganize/reformat and filter the data.
XDRConvert is optionally used to regrid, aggregate and reformat the data
to other file formats.
The original data files and RSIG command used to process them are shown below.
';


# Query string parsing routine dispatch table:

my %parsers = (
 'service'       => \&parse_service_option,
 'version'       => \&parse_version_option,
 'request'       => \&parse_request_option,
 'coverage'      => \&parse_coverage_option,
 'time'          => \&parse_time_option,
 'bbox'          => \&parse_bbox_option,
 'format'        => \&parse_format_option,
 'stride'        => \&parse_stride_option,
 'sparse'        => \&parse_sparse_option,
 'corners'       => \&parse_corners_option,
 'compress'      => \&parse_compress_option,
 'regrid'        => \&parse_regrid_option,
 'regrid_aggregate' => \&parse_regrid_aggregate_option,
 'lambert'       => \&parse_lambert_option,
 'stereographic' => \&parse_stereographic_option,
 'mercator'      => \&parse_mercator_option,
 'lonlat'        => \&parse_lonlat_option,
 'ellipsoid'     => \&parse_ellipsoid_option,
 'grid'          => \&parse_grid_option
);

# Webserver content types for each output format:

my %content_types = (
 'xdr'           => 'application/octet-stream',
 'ascii'         => 'text/plain',
 'netcdf-coards' => 'application/netcdf',
 'netcdf-ioapi'  => 'application/netcdf',
 'original'      => 'application/octet-stream' # Compressed tar file (tgz).
);

# Variable and units lists for each mod file type.
#
# Note: WCS will still ignore case, but use capitalized variable names here
# (as they appear in the MODIS files).
# Except change any '.' in the variable name to 'p' so it works with
# COVERAGE naming scheme. For example:
# Optical_Depth_Ratio_Small_Ocean_0.55micron
# is referred to below as:
# Optical_Depth_Ratio_Small_Ocean_0p55micron
# until passed to MODISSubset (when the dot is put back).

my @mod4_variables = (
  'AOD_550_Dark_Target_Deep_Blue_Combined',
  'AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag',
  'AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag',
  'Aerosol_Cldmask_Land_Ocean',
  'Aerosol_Cloud_Fraction_Land',
  'Aerosol_Cloud_Fraction_Ocean',
  'Aerosol_Type_Land',
  'Angstrom_Exponent_1_Ocean',
  'Angstrom_Exponent_2_Ocean',
  'Average_Cloud_Pixel_Distance_Land_Ocean',
  'Cloud_Pixel_Distance_Land_Ocean',
  'Corrected_Optical_Depth_Land_wav2p1',
  'Deep_Blue_Aerosol_Optical_Depth_550_Land',
  'Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate',
  'Deep_Blue_Aerosol_Optical_Depth_550_Land_Estimated_Uncertainty',
  'Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag',
  'Deep_Blue_Aerosol_Optical_Depth_550_Land_STD',
  'Deep_Blue_Algorithm_Flag_Land',
  'Deep_Blue_Angstrom_Exponent_Land',
  'Deep_Blue_Cloud_Fraction_Land',
  'Deep_Blue_Number_Pixels_Used_550_Land',
  'Effective_Optical_Depth_0p55um_Ocean',
  'Fitting_Error_Land',
  'Glint_Angle',
  'Image_Optical_Depth_Land_And_Ocean',
  'Land_Ocean_Quality_Flag',
  'Land_sea_Flag',
  'Mass_Concentration_Land',
  'Optical_Depth_Land_And_Ocean',
  'Optical_Depth_Ratio_Small_Land',
  'Optical_Depth_Ratio_Small_Ocean_0p55micron',
  'Scan_Start_Time',
  'Scattering_Angle',
  'Sensor_Azimuth',
  'Sensor_Zenith',
  'Solar_Azimuth',
  'Solar_Zenith',
  'Topographic_Altitude_Land',
  'Wind_Speed_Ncep_Ocean'
);

my @mod4_descriptions = (
  'Combined Dark Target, Deep Blue AOT at 0.55 micron for land and ocean.',
  'Combined Dark Target, Deep Blue AOT at 0.55 micron Algorithm Flag' .
    ' (0=Dark Target, 1=Deep Blue, 2=Mixed)',
  'Combined Dark Target, Deep Blue Aerosol Confidence Flag' .
    ' (0= No Confidence (or fill),  1= Marginal, 2= Good, 3= Very Good)',
  'Aerosol Cloud Mask 500 meter resolution 0 = cloud 1 = clear',
  'Cloud fraction from Land aerosol cloud mask' .
    ' from retrieved and overcast pixels not including cirrus mask',
  'Cloud fraction from Land aerosol cloud mask' .
    ' from retrieved and overcast pixels not including cirrus mask',
  'Aerosol Type: 1 = Continental, 2 = Moderate Absorption Fine,' .
    ' 3 = Strong Absorption Fine,4 = Weak Absorption Fine,' .
    ' 5 = Dust Coarse',
  'Calculated Angstrom Exponent for 0.55 vs 0.86 micron for Average Solution',
  'Calculated Angstrom Exponent for 0.86 vs 2.13 micron for Average Solution',
  'Average Distance (number of pixels) to nearest pixel identified as cloudy' .
    ' from each clear pixel used for Aerosol Retrieval in 10 km retrieval box',
  'Distance (number of pixels) to nearest pixel identified as cloudy' .
    ' (500 m resolution)',
  'Retrieved  AOT at 2.13 micron',
  'AOT at 0.55 micron for land  with all quality data (Quality flag=1,2,3)',
  'Deep Blue AOT at 0.55 micron for land with higher quality data' .
    ' (Quality flag=2,3)',
  'Estimated uncertainty (one-sigma confidence envelope) of Deep Blue AOT' .
    ' at 0.55 micron for land for all quality data (Quality flag=1,2,3)',
  'Deep Blue Aerosol Confidence Flag (0= No Confidence (or fill),' .
    ' 1= Marginal, 2= Good, 3= Very Good)',
  'Standard deviation of Deep Blue AOT at 0.55 micron for land with' .
    ' all quality data (Quality flag=1,2,3)',
  'Deep Blue Aerosol Algorithm Flag (0=DeepBlue, 1=Vegetated, 2=Mixed)',
  'Deep Blue Angstrom Exponent for land with all quality data' .
    ' (Quality flag=1,2,3)',
  'Cloud fraction from Deep Blue Aerosol cloud mask over land',
  'Number of pixels used for AOT retrieval 0.55 micron for land',
  'Retrieved AOT for average solution at 0.55um For easy L3 processing',
  'Spectral Fitting error for inversion over land',
  'Glint angle',
  'AOT at 0.55 micron for both ocean (Average) and land (corrected) with' .
    ' all quality data (Quality flag=0,1,2,3)',
  'Quality Flag for Land and ocean Aerosol retreivals 0=bad 1= Marginal' .
    ' 2=Good 3=Very Good)',
  'Land_sea_Flag(based on MOD03 Landsea mask 0 = Ocean, 1 = Land and' .
    ' Ephemeral water 2 =Coastal)',
  'Estimated Column Mass(per area) using assumed mass extinction efficiency',
  'AOT at 0.55 micron for both ocean (Average) (Quality flag=1,2,3) and' .
    ' land (corrected) (Quality flag=3)',
  'Fraction of AOT contributed by fine dominated model',
  'Fraction of AOT (at 0.55 micron) contributed by fine mode for average' .
    ' solution',
  'TAI Time at Start of Scan replicated across the swath',
  'Scattering angle',
  'Sensor azimuth angle, cell to sensor',
  'Sensor zenith angle, cell to sensor',
  'Solar azimuth angle, cell to sun',
  'Solar zenith angle, cell to sun',
  'Averaged topographic altitude (in km) for Land',
  'Wind Speed based on NCEP reanalysis for Ocean'
);

my @mod4_units = (
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  'pixels',
  'pixels',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  'pixels',
  '-',
  '-',
  'deg',
  '-',
  '-',
  '-',
  '1.0e-6g/cm^2',
  '-',
  '-',
  '-',
  'Seconds_since_1993-01-01',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'km',
  'm/s'
);



my @mod43k_variables = (
  'Aerosol_Cloud_Fraction_Land',
  'Aerosol_Cloud_Fraction_Ocean',
  'Aerosol_Type_Land',
  'Angstrom_Exponent_1_Ocean',
  'Angstrom_Exponent_2_Ocean',
  'Corrected_Optical_Depth_Land_wav2p1',
  'Fitting_Error_Land',
  'Glint_Angle',
  'Image_Optical_Depth_Land_And_Ocean',
  'Land_Ocean_Quality_Flag',
  'Land_sea_Flag',
  'Mass_Concentration_Land',
  'Number_Pixels_Used_Ocean',
  'Optical_Depth_Land_And_Ocean',
  'Optical_Depth_Ratio_Small_Land',
  'Optical_Depth_Ratio_Small_Ocean_0p55micron',
  'Scan_Start_Time',
  'Scattering_Angle',
  'Sensor_Azimuth',
  'Sensor_Zenith',
  'Solar_Azimuth',
  'Solar_Zenith',
  'Topographic_Altitude_Land',
  'Wind_Speed_Ncep_Ocean'
);

my @mod43k_descriptions = (
  'Cloud fraction from Land aerosol cloud mask from retrieved and ' .
    'overcast pixels not including cirrus mask',
  'Cloud fraction from land aerosol cloud mask from retrieved and ' .
    'overcast pixels not including cirrus mask',
  'Aerosol Type: 1 = continental, 2 = moderate absorption fine, ' .
    '3 = strong bbsorption fine, 4 = weak absorption fine, ' .
    '5 = dust coarse',
  'Calculated angstrom exponent for 0.55 vs 0.86 microns for average solution',
  'Calculated angstrom exponent for 0.86 vs 2.13 microns for average solution',
  'Retrieved AOT at 2.13 microns',
  'Spectral fitting error for inversion over land',
  'Glint angle',
  'AOT at 0.55 micron for both ocean (average) and ' .
    'land (corrected) with all quality data (quality flag=0,1,2,3)',
  'Quality flag for land and ocean aerosol retreivals: ' .
    '0 = bad, 1 = marginal, 2 = good, 3 = very good',
  'Land_sea_Flag(based on MOD03 Landsea mask: ' .
    '0 = Ocean, 1 = Land and Ephemeral water, 2 = Coastal',
  'Estimated Column Mass(per area) using assumed mass extinction efficiency',
  'Number of pixels used for ocean retrieval at 865 nm',
  'AOT at 0.55 micron for both ocean (average) (quality flag=1,2,3) and ' .
    'land (corrected) (quality flag=3)',
  'Fraction of AOT contributed by fine dominated model',
  'Fraction of AOT (at 0.55 microns) contributed by fine mode for ' .
    'average solution',
  'TAI Time at Start of Scan replicated across the swath',
  'Scattering angle',
  'Sensor azimuth angle, cell to sensor',
  'Sensor zenith angle, cell to sensor',
  'Sensor azimuth angle, cell to Sun',
  'Sensor zenith angle, cell to Sun',
  'Terrain height above mean sea level',
  'Wind speed based on NCEP reanalysis for ocean'
);

my @mod43k_units = (
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  'deg',
  '-',
  '-',
  '-',
  '1.0e-6g/cm^2',
  '-',
  '-',
  '-',
  '-',
  'Seconds_since_1993-01-01',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'km',
  'm/s'
);



my @mod6_variables = (
  'Cloud_Effective_Emissivity',
  'Cloud_Effective_Emissivity_Day',
  'Cloud_Effective_Emissivity_Nadir',
  'Cloud_Effective_Emissivity_Nadir_Day',
  'Cloud_Effective_Emissivity_Nadir_Night',
  'Cloud_Effective_Emissivity_Night',
  'Cloud_Fraction',
  'Cloud_Fraction_Day',
  'Cloud_Fraction_Nadir',
  'Cloud_Fraction_Nadir_Day',
  'Cloud_Fraction_Nadir_Night',
  'Cloud_Fraction_Night',
  'Cloud_Height_Method',
  'Cloud_Phase_Infrared',
  'Cloud_Phase_Infrared_Day',
  'Cloud_Phase_Infrared_Night',
  'Cloud_Top_Height',
  'Cloud_Top_Height_Nadir',
  'Cloud_Top_Height_Nadir_Day',
  'Cloud_Top_Height_Nadir_Night',
  'Cloud_Top_Pressure',
  'Cloud_Top_Pressure_Day',
  'Cloud_Top_Pressure_Infrared',
  'Cloud_Top_Pressure_Nadir',
  'Cloud_Top_Pressure_Nadir_Day',
  'Cloud_Top_Pressure_Nadir_Night',
  'Cloud_Top_Pressure_Night',
  'Cloud_Top_Temperature',
  'Cloud_Top_Temperature_Day',
  'Cloud_Top_Temperature_Nadir',
  'Cloud_Top_Temperature_Nadir_Day',
  'Cloud_Top_Temperature_Nadir_Night',
  'Cloud_Top_Temperature_Night',
  'Radiance_Variance',
  'Scan_Start_Time',
  'Sensor_Azimuth',
  'Sensor_Azimuth_Day',
  'Sensor_Azimuth_Night',
  'Sensor_Zenith',
  'Sensor_Zenith_Day',
  'Sensor_Zenith_Night',
  'Solar_Azimuth',
  'Solar_Azimuth_Day',
  'Solar_Azimuth_Night',
  'Solar_Zenith',
  'Solar_Zenith_Day',
  'Solar_Zenith_Night',
  'Surface_Pressure',
  'Surface_Temperature',
  'Tropopause_Height'
);

my @mod6_descriptions = (
  'Cloud Effective Emissivity from Cloud Top Pressure Retrieval',
  'Cloud Effective Emissivity from Cloud Top Pressure Retrieval, Day Only',
  'Cloud Effective Emissivity from Cloud Top Pressure Retrieval ' .
    'for Sensor Zenith (View) Angles <= 32 Degrees',
  'Cloud Effective Emissivity from Cloud Top Pressure Retrieval for ' .
    'Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only',
  'Cloud Effective Emissivity from Cloud Top Pressure Retrieval for ' .
    'Sensor Zenith (View) Angles <= 32 Degrees, Night Data Only',
  'Cloud Effective Emissivity from Cloud Top Pressure Retrieval, Night Only',
  'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask',
  'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask,'.
    ' Day Only',
  'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask '.
    'for Sensor Zenith (View) Angles <= 32 Degrees',
  'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask '.
    'for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only',
  'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask '.
    'for Sensor Zenith (View) Angles <= 32 Degrees, Night Data Only',
  'Cloud Fraction in Retrieval Region (5x5 1-km Pixels) from 1-km Cloud Mask, '.
    'Night Only',
  'Index Indicating MODIS Bands Used for Cloud Top Pressure Retrieval',
  'Cloud Phase from 8.5 and 11 um Bands',
  'Cloud Phase from 8.5 and 11 um Bands, Day Only',
  'Cloud Phase from 8.5 and 11 um Bands, Night Only',
  'Geopotential Height at Retrieved Cloud Top Pressure Level ' .
    '(rounded to nearest 50m)',
  'Geopotential Height at Retrieved Cloud Top Pressure Level ' .
    'for Sensor Zenith (View) Angles <=32 Degrees (rounded to nearest 50m)',
  'Geopotential Height at Retrieved Cloud Top Pressure Level ' .
    'for Sensor Zenith (View) Angles <=32 Degrees, Day Data Only ' .
    '(rounded to nearest 50m)',
  'Geopotential Height at Retrieved Cloud Top Pressure Level ' .
    'for Sensor Zenith (View) Angles <=32 Degrees, Night Data Only ' .
    '(rounded to nearest 50m)',
  'Cloud Top Pressure Level (rounded to nearest 5 hPa)',
  'Cloud Top Pressure Level, Day Data Only (rounded to nearest 5 hPa)',
  'Cloud Top Pressure from IR Window Retrieval',
  'Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees ' .
    '(rounded to nearest 5 hPa)',
  'Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees ' .
    '(rounded to nearest 5 hPa), Day Data Only',
  'Cloud Top Pressure Level for Sensor Zenith (View) Angles <= 32 Degrees ' .
    '(rounded to nearest 5 hPa), Night Data Only',
  'Cloud Top Pressure Level, Night Data Only (rounded to nearest 5 hPa)',
  'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level',
  'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level, ' .
    'Day Only',
  'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level ' .
    'for Sensor Zenith (View) Angles <= 32 Degrees',
  'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level ' .
    'for Sensor Zenith (View) Angles <= 32 Degrees, Day Data Only',
  'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level ' .
    'for Sensor Z enith (View) Angles <= 32 Degrees, NightData Only',
  'Temperature from Ancillary Data at Retrieved Cloud Top Pressure Level, ' .
    'Night Only',
  'Band 31 Radiance Standard Deviation',
  'TAI Time at Start of Scan replicated across the swath',
  'Sensor Azimuth Angle, Cell to Sensor',
  'Sensor Azimuth Angle, Cell to Sensor, Day Data Only',
  'Sensor Azimuth Angle, Cell to Sensor, Night Data Only',
  'Solar Zenith Angle, Cell to Sun',
  'Solar Zenith Angle, Cell to Sun, Day Data Only',
  'Solar Zenith Angle, Cell to Sun, Night Data Only',
  'Solar Azimuth Angle, Cell to Sun',
  'Solar Azimuth Angle, Cell to Sun, Day Data Only',
  'Solar Azimuth Angle, Cell to Sun, Night Data Only',
  'Sensor Zenith Angle, Cell to Sensor',
  'Sensor Zenith Angle, Cell to Sensor, Day Data Only',
  'Sensor Zenith Angle, Cell to Sensor, Night Data Only',
  'Surface Pressure from Ancillary Data',
  'Surface Temperature from Ancillary Data',
  'Tropopause Height from Ancillary Data'
);

my @mod6_units = (
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  '-',
  'm',
  'm',
  'm',
  'm',
  'hPa',
  'hPa',
  'hPa',
  'hPa',
  'hPa',
  'hPa',
  'hPa',
  'K',
  'K',
  'K',
  'K',
  'K',
  'K',
  'W/m2/sr/micron',
  'Seconds_since_1993-01-01',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'deg',
  'hPa',
  'K',
  'hPa'
);




my @mod7_variables = (
  'Cloud_Mask',
  'K_Index',
  'Lifted_Index',
  'Scan_Start_Time',
  'Sensor_Azimuth',
  'Sensor_Zenith',
  'Skin_Temperature',
  'Solar_Azimuth',
  'Solar_Zenith',
  'Surface_Elevation',
  'Surface_Pressure',
  'Total_Ozone',
  'Total_Totals',
  'Tropopause_Height',
  'Water_Vapor',
  'Water_Vapor_Direct',
  'Water_Vapor_High',
  'Water_Vapor_Low'
);

my @mod7_descriptions = (
  'MODIS Cloud Mask, First Byte',
  'K Index',
  'Lifted Index',
  'TAI Time at Start of Scan replicated across the swath',
  'Sensor Azimuth Angle, Cell to Sensor',
  'Sensor Zenith Angle, Cell to Sensor',
  'Skin Temperature',
  'Solar Azimuth Angle, Cell to Sun',
  'Solar Zenith Angle, Cell to Sun',
  'Surface Elevation',
  'Surface Pressure',
  'Total Ozone Burden',
  'Total Totals',
  'Tropopause Height',
  'Total Column Precipitable Water Vapor - IR Retrieval',
  'Total Column Precipitable Water Vapor - Direct IR Retrieval',
  'Precipitable Water Vapor High - IR Retrieval',
  'Precipitable Water Vapor Low - IR Retrieval'
);

my @mod7_units = (
  '-',
  'K',
  'K',
  'Seconds_since_1993-01-01',
  'deg',
  'deg',
  'K',
  'deg',
  'deg',
  'm',
  'hPa',
  'Dobson',
  'K',
  'hPa',
  'cm',
  'cm',
  'cm',
  'cm'
);



my %mod_metadata = (
 'mod4' => {
   'variables'    => \@mod4_variables,
   'units'        => \@mod4_units,
   'descriptions' => \@mod4_descriptions,
   'begin_date'   => '2001-01-01T00:00:00Z',
   'end_date'     => '2017-12-31T23:59:59Z',
   'directory'    => '/rsig/MODIS'
  },
 'mod6' => {
   'variables'    => \@mod6_variables,
   'units'        => \@mod6_units,
   'descriptions' => \@mod6_descriptions,
   'begin_date'   => '2001-08-01T00:00:00Z',
   'end_date'     => '2017-08-31T23:59:59Z',
   'directory'    => '/rsig2/MODIS'
 },
 'mod7' => {
   'variables'    => \@mod7_variables,
   'units'        => \@mod7_units,
   'descriptions' => \@mod7_descriptions,
   'begin_date'   => '2001-08-01T00:00:00Z',
   'end_date'     => '2017-08-31T23:59:59Z',
   'directory'    => '/rsig2/MODIS'
 },
 'mod43k' => {
   'variables'    => \@mod43k_variables,
   'units'        => \@mod43k_units,
   'descriptions' => \@mod43k_descriptions,
   'begin_date'   => '2001-01-01T00:00:00Z',
   'end_date'     => '2017-12-31T23:59:59Z',
   'directory'    => '/rsig/MODIS'
  },
);

my @mods = keys %mod_metadata; # mod4, mod6, mod7, mod43k.


################################## VARIABLES #################################


# Parsed from the URL query string:

my $service       = ''; # wcs.
my $version       = ''; # 1.0.0.
my $request       = ''; # getcapabilities or describecoverage or getcoverage.
my $variable      = ''; # mod4.optical_depth_land_and_ocean, ...
my $format        = ''; # xdr, ascii, netcdf, original.
my $compress      = ''; # 1 = | gzip -c otherwise don't compress (default).
my $time          = ''; # E.g., 2001-08-29t00:00:00z/2001-08-30t23:59:59z.
my $bbox          = ''; # E.g., -90,28,-80,32,0,0.
my $stride        = ''; # Stride between points in domain. E.g., 10.
my $sparse        = ''; # Sparse target count of points in domain. E.g., 100.
my $corners       = ''; # 1 = compute corners otherwise don't (default).
my $regrid        = ''; # E.g., nearest, mean, weighted.
my $regrid_aggregate = ''; # E.g., none, all, daily.
my $lambert       = ''; # E.g., 33,45,-97,40.
my $stereographic = ''; # E.g., -98,90,45.
my $mercator      = ''; # E.g., -98.
my $lonlat        = '';
my $ellipsoid     = ''; # E.g., 6370997,6370997.
my $grid          = ''; # E.g., 268,259,-420000,-1716000,12000,12000.

# Derived from the above parsed values:

my $starting_timestamp = 0;  # yyyymmddhh, e.g., 2005082600.
my $hours              = 0;  # E.g., 5 days = 5 x 24 = 120.
my $command            = ''; # Subset command to run on the remote host.
my $mod                = ''; # mod4 or mod6 or mod7.

my $list_file = ''; # Pathed file containing list of MODIS files to read.
my $temp_dir  = ''; # Dir where temp files are created/removed by MODISSubset.

################################## ROUTINES ##################################

# routine to sort by PGE_StartTime
sub by_pge_starttime
{
  my $asplit = index($a,".A");
  my $asat = substr($a,0,$asplit);
  my $at = substr($a,$asplit);
  my $bsplit = index($b,".A");
  my $bsat = substr($b,0,$bsplit);
  my $bt = substr($b,$bsplit);
  return "$at$asat" cmp "$bt$bsat";
}

# http header
my $has_printed_hdr = 0;
sub http_header {
  return if $has_printed_hdr;
  my $hdr = "Content-type: text/plain; charset=iso-8859-1\n\n";
  if ($format && $content_types{$format})
  {
    $hdr = "Content-type: $content_types{$format}; charset=iso-8859-1\n\n";
  }

  print $hdr;
  $has_printed_hdr = -1;
}

# usage
my $sent_header = 0;
sub WARN_USER
{
  $format = '';  # reset the format, since all warnings are text
  http_header() unless $sent_header;
  $sent_header = 'sent';
  print join("\n", @_), "\n";
}

sub LOG
{
  print STDERR join("\n", @_), "\n";
}

sub check_cgi
{
  my ($q) = @_;
  die "no q" unless $q;
  my $result = 1;
  my $qstr = '';
  my @names = $q->param;
  foreach my $option (@names) {
    my $value = $q->param( $option );
    my $str = parse_arg("$option=$value");
    if ($str) {
      $qstr = join(" ", $qstr, $str);
    }
    else
    {
      $result = 0;
    }
  }
  # check for existence of required options
  $result = $result && required_options_specified();

  #  print header for httpd to strip-off:
  http_header();

  return $qstr;
}

sub check_args
{
  my $result = 1;
  foreach my $item (@_)
  {
    my $str = parse_arg($item);
    if ($str) {
      $result = $result && 1;
    }
    else
    {
      $result = 0;
    }
  }
  $result = $result && required_options_specified();
  return $request if $result;
  return 0;
}


# validate parameters, and assign into variables.

sub parse_arg {
 my ($item) = @_;
 return undef unless $item;

 my ($option,$value) = split(/\s*=\s*/, $item);
 $option =~ s/[^A-Za-z]/_/go;                        # Untaint.
 $value =~ s/[^\w\-.,:\/]/_/go;                      # Untaint.
 my $lowercase_option = lc( $option );
 my $lowercase_value  = lc( $value );
 debug( "$lowercase_option $lowercase_value" );

 if ( $parsers{ $lowercase_option } ) {
   my $result = $parsers{ $lowercase_option }->( $lowercase_value );
   $item = undef unless $result;
 } else {
   WARN_USER( "unknown option");
   $item = undef;
 }
 return $item;
}

# print examples
sub examples {
 WARN_USER( "\n\nUsage Examples:\n");
 WARN_USER( "wget -q -c -T 0 -O - '$service_origin?SERVICE=wcs&VERSION=1.0.0");
 WARN_USER( "&REQUEST=GetCapabilities'\n");
 WARN_USER( "wget -q -c -T 0 -O - '$service_origin?SERVICE=wcs&VERSION=1.0.0");
 WARN_USER( "&REQUEST=DescribeCoverage'\n");
 WARN_USER( "wget -q -c -T 0 -O - '$service_origin?SERVICE=wcs&VERSION=1.0.0");
 WARN_USER( "&REQUEST=GetCoverage");
 WARN_USER( "&COVERAGE=mod4.optical_depth_land_and_ocean");
 WARN_USER( "&TIME=2005-08-26T00:00:00Z/2005-08-30T23:59:59Z");
 WARN_USER( "&BBOX=-90,28,-80,32,0,0&FORMAT=ascii'\n");
 WARN_USER( "wget -q -c -T 0 -O - '$service_origin?SERVICE=wcs&VERSION=1.0.0");
 WARN_USER( "&REQUEST=GetCoverage");
 WARN_USER( "&COVERAGE=mod6.cloud_optical_thickness");
 WARN_USER( "&TIME=2005-08-26T00:00:00Z/2005-08-30T23:59:59Z");
 WARN_USER( "&BBOX=-90,28,-80,32,0,0&STRIDE=10&FORMAT=xdr'");
 WARN_USER( "  > modis_cot.xdr ; head -12 modis_cot.xdr\n");
 WARN_USER( "wget -q -c -T 0 -O - '$service_origin?SERVICE=wcs&VERSION=1.0.0");
 WARN_USER( "&REQUEST=GetCoverage");
 WARN_USER( "&COVERAGE=mod7.total_ozone");
 WARN_USER( "&TIME=2005-08-26T00:00:00Z/2005-08-30T23:59:59Z");
 WARN_USER( "&BBOX=-90,28,-80,32,0,0&FORMAT=netcdf-coards'");
 WARN_USER( "  > modis_ozone.nc ; ncdump modis_ozone.nc | more\n");
 WARN_USER( "wget -q -c -T 0 -O - '$service_origin?SERVICE=wcs&VERSION=1.0.0");
 WARN_USER( "&REQUEST=GetCoverage");
 WARN_USER( "&COVERAGE=mod7.total_ozone");
 WARN_USER( "&TIME=2005-08-26T00:00:00Z/2005-08-30T23:59:59Z");
 WARN_USER( "&BBOX=-90,28,-80,32,0,0&FORMAT=original'");
 WARN_USER( "  > modis_files.tgz ; gtar -zxf modis_files.tgz\n");
 WARN_USER( "wget -q -c -T 0 -O - '$service_origin?SERVICE=wcs&VERSION=1.0.0");
 WARN_USER( "&REQUEST=GetCoverage&COVERAGE=mod7.total_ozone");
 WARN_USER( "&TIME=2005-08-26T00:00:00Z/2005-08-26T23:59:59Z");
 WARN_USER( "&BBOX=-130,20,-60,55,0,0&FORMAT=netcdf-ioapi");
 WARN_USER( "&REGRID=weighted&LAMBERT=33,45,-97,40");
 WARN_USER( "&ELLIPSOID=6370997,6370997");
 WARN_USER( "&GRID=268,259,-420000,-1716000,12000,12000'");
 WARN_USER( "  > modis_ozone.ncf ; ncdump modis_ozone.ncf | more\n");
 WARN_USER( "Note: these variables are implicitly included:");
 WARN_USER( "  longitude(degrees_east)");
 WARN_USER( "  latitude(degrees_north)");
 WARN_USER( "  scan_start_time(seconds_since_1993-1-1_00:00:00.0_0)");
 WARN_USER( );
 return 0;
}

# Print web server capabilities metadata.

sub print_capabilities {
 print '<WCS_Capabilities version="1.0.0" ';
 print 'xmlns="http://www.opengeospatial.org/standards/wcs" ';
 print 'xmlns:gml="http://www.opengis.net//gml" ';
 print 'xmlns:xlink="http://www.w3.org/1999/xlink">';
 print <<"END_CAP";
   <Service>
       <metadataLink xlink:type="simple" xlink:href="http://www.epa.gov/airdata" metadataType="other" />
       <description>EPA MODIS Web Server 1.0.0</description>
       <name>EPA_MODIS_OGC_WCS_1.0.0</name>
       <label>EPA MODIS Web Server 1.0.0</label>
       <keywords>
           <keyword>EPA</keyword>
           <keyword>MODIS</keyword>
           <keyword>interoperability</keyword>
       </keywords>
       <responsibleParty>
           <individualName>Todd Plessel</individualName>
           <organisationName>EPA Vislab</organisationName>
           <contactInfo>
               <onlineResource xlink:type="simple" xlink:href="mailto:plessel.todd\@epa.gov" />
           </contactInfo>
       </responsibleParty>
       <fees>NONE</fees>
       <accessConstraints>NONE</accessConstraints>
   </Service>
   <Capability>
       <Request>
           <GetCapabilities>
               <DCPType>
                   <HTTP>
                       <Get>
                           <OnlineResource xlink:type="simple" xlink:href="$service_origin?" />
                       </Get>
                       <Post>
                           <OnlineResource xlink:type="simple" xlink:href="$service_origin" />
                       </Post>
                   </HTTP>
               </DCPType>
           </GetCapabilities>
           <DescribeCoverage>
               <DCPType>
                   <HTTP>
                       <Get>
                           <OnlineResource xlink:type="simple" xlink:href="$service_origin?" />
                       </Get>
                       <Post>
                           <OnlineResource xlink:type="simple" xlink:href="$service_origin" />
                       </Post>
                   </HTTP>
               </DCPType>
           </DescribeCoverage>
           <GetCoverage>
               <DCPType>
                   <HTTP>
                       <Get>
                           <OnlineResource xlink:type="simple" xlink:href="$service_origin?" />
                       </Get>
                       <Post>
                           <OnlineResource xlink:type="simple" xlink:href="$service_origin" />
                       </Post>
                   </HTTP>
               </DCPType>
           </GetCoverage>
       </Request>
       <Exception>
           <Format>text/xml; charset="utf-8"</Format>
       </Exception>
   </Capability>
   <ContentMetadata version="1.0.0">

END_CAP

 foreach my $this_mod ( @mods ) {
   my $metadata         = $mod_metadata{ $this_mod };
   my $mod_variables    = $metadata->{ variables };
   my $mod_units        = $metadata->{ units };
   my $mod_descriptions = $metadata->{ descriptions };
   my $variable_count   = @$mod_variables;

   for ( my $variable_index = 0; $variable_index < $variable_count;
         ++$variable_index ) {
     my $this_variable     = @$mod_variables[ $variable_index ];
     my $this_units        = @$mod_units[ $variable_index ];
     my $this_description  = @$mod_descriptions[ $variable_index ];
     print "        <CoverageOfferingBrief>\n";
     print "            <name>$this_mod.$this_variable</name>\n";
     print "            <label>$this_variable($this_units)</label>\n";
     print "            <description>$this_description</description>";
     print '
           <lonLatEnvelope srsName="WGS84(DD)">
               <gml:pos>-180 -90</gml:pos>
               <gml:pos>180 90</gml:pos>
           </lonLatEnvelope>
       </CoverageOfferingBrief>';
     print "\n";
   }
 }

 print "    </ContentMetadata>\n";
 print "</WCS_Capabilities>\n";
}



# Print web server coverage description metadata.

sub print_coverage_description {
 print '<CoverageDescription version="1.0.0" ';
 print 'xmlns="http://www.opengeospatial.org/standards/wcs" ';
 print 'xmlns:gml="http://www.opengis.net/gml" ';
 print 'xmlns:xlink="http://www.w3.org/1999/xlink">';
 print "\n";

 foreach my $this_mod ( @mods ) {
   my $metadata         = $mod_metadata{ $this_mod };
   my $mod_variables    = $metadata->{ variables };
   my $mod_units        = $metadata->{ units };
   my $mod_descriptions = $metadata->{ descriptions };
   my $begin_date       = $metadata->{ begin_date };
   my $end_date         = $metadata->{ end_date };
   my $variable_count   = @$mod_variables;

   for ( my $variable_index = 0; $variable_index < $variable_count;
         ++$variable_index ) {
     my $this_variable     = @$mod_variables[ $variable_index ];
     my $this_units        = @$mod_units[ $variable_index ];
     my $this_description  = @$mod_descriptions[ $variable_index ];

     if ( $variable eq '' || $variable eq lc( "$this_mod.$this_variable" ) ) {
       print "    <CoverageOffering>\n";
       print "        <name>$this_mod.$this_variable</name>\n";
       print "        <label>$this_variable($this_units)</label>\n";
       print "        <description>$this_description</description>";
       print '
       <domainSet>
           <spatialDomain>
               <gml:Envelope srsName="WGS84(DD)">
                   <gml:pos>-180 -90</gml:pos>
                   <gml:pos>180 90</gml:pos>
               </gml:Envelope>
           </spatialDomain>
           <temporalDomain>
               <timePeriod>';
       print "\n";
       print "                    <beginPosition>$begin_date</beginPosition>\n";
       print "                    <endPosition>$end_date</endPosition>\n";
       print '
                   <timeResolution>PT1H</timeResolution>
               </timePeriod>
           </temporalDomain>
       </domainSet>
       <rangeSet>
           <RangeSet>';
       print "\n";
       print "                <name>$this_mod.$this_variable</name>\n";
       print "                <label>$this_variable($this_units)</label>\n";
       print "                <description>$this_description</description>";
       print '
               <nullValues>
                   <singleValue>-9999.0</singleValue>
               </nullValues>
           </RangeSet>
       </rangeSet>
       <supportedCRSs>
           <requestResponseCRSs>CRS:84</requestResponseCRSs>
           <nativeCRSs>CRS:84</nativeCRSs>
       </supportedCRSs>
       <supportedFormats>
           <formats>ASCII XDR NetCDF-COARDS NetCDF-IOAPI HDF</formats>
       </supportedFormats>
       <supportedInterpolations>
           <interpolationMethod>none</interpolationMethod>
       </supportedInterpolations>
   </CoverageOffering>
';
     }
   }
 }

 print "</CoverageDescription>\n";
}

# Parse service option.

sub parse_service_option {
 my $value = shift;
 return parse_option( $service, $value, 'SERVICE', 'wcs' );
}



# Parse version option.

sub parse_version_option {
 my $value = shift;
 my $result = parse_option( $version, $value, 'VERSION', '1.0.0' );
 return $result;
}



# Parse request option.

sub parse_request_option {
 my $value = shift;
 my $result = parse_option(
  $request,
  $value,
  'REQUEST',
  'getcoverage getcapabilities describecoverage getmetadata'
 );
 return $result;
}



# Parse format option:

sub parse_format_option {
 my $value = shift;
 my $result = parse_option( $format, $value, 'FORMAT',
                            'xdr ascii netcdf-coards netcdf-ioapi original' );
 return $result;
}



# Parse coverage option.

sub parse_coverage_option {
 my $value = shift;
 my $result = 0;

 if ( $variable ne '' ) {
  WARN_USER( "\nRedundant COVERAGE option.");
 } else {
   my @components = split( /\./, $value );
   my $count = @components;

   if ( $count == 2 ) {
     $mod = $components[ 0 ];
     my $metadata = $mod_metadata{ $mod };

     if ( $metadata ) {
       my $variables = $metadata->{ variables };
       my @selected_variables = split( /\,/, $components[ 1 ] );

       foreach my $selected_variable ( @selected_variables ) {
         my $variable_count = @$variables;
         my $found = 0;

         # Compare lowercase selected_variable to capitalized file variables:

         for ( my $variable_index = 0;
               $variable_index < $variable_count && $found == 0;
               ++$variable_index ) {
           my $this_variable = @$variables[ $variable_index ];

           if ( $selected_variable eq lc( $this_variable ) ) {
             $found = 1;
             # HACK for problematic MODIS file variable names containing dot:
             $this_variable =~ s/_0p55micron/_0.55micron/; # Use dot name now.
             $variable .= ' ' . $this_variable; # Use capitalized name.
           }
         }
 
         if ( $found == 0 ) {
           WARN_USER( "\nInvalid COVERAGE variable." );
           $result = -1;
         }

         debug( "$mod $variable" );
       }
     }
   }

   $result += 1;
 }

 if ( $result == 0 ) {
   WARN_USER( "\nInvalid COVERAGE option.");
 }

 return $result;
}



# Parse time option.

sub parse_time_option {
 my $value = shift;
 my $result = 0;

 if ( $time ne '' ) {
   WARN_USER( "\nRedundant TIME option.");
 } else {
   $time = $value;
   my $is_valid = is_valid_time( $time );

   if ( ! $is_valid ) {
     WARN_USER( "\nInvalid TIME option.");
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Parse bbox option.

sub parse_bbox_option {
 my $value = shift;
 my $result = 0;

 if ( $bbox ne '' ) {
   WARN_USER( "\nRedundant BBOX option.");
 } else {
   $bbox = $value;
   my @bounds = split( /,/, $bbox );
   my $bounds_count = @bounds;

   if ( $bounds_count != 4 && $bounds_count != 6 ) {
     WARN_USER( "\nInvalid bbox option.");
   } elsif ( ! in_range( $bounds[ 0 ], -180.0, 180.0 ) ) {
     WARN_USER( "\nInvalid bbox option.");
   } elsif ( ! in_range( $bounds[ 1 ], -90.0, 90.0 ) ) {
     WARN_USER( "\nInvalid bbox option.");
   } elsif ( ! in_range( $bounds[ 2 ], $bounds[ 0 ], 180.0 ) ) {
     WARN_USER( "\nInvalid bbox option.");
   } elsif ( ! in_range( $bounds[ 3 ], $bounds[ 1 ], 90.0 ) ) {
     WARN_USER( "\nInvalid bbox option.");
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Parse stride option.

sub parse_stride_option {
 my $value = shift;
 my $result = 0;

 if ( $stride ne '' || $sparse ne '' ) {
   WARN_USER( "\nRedundant STRIDE option.");
 } else {
   $stride = $value;
   my $is_valid = $stride > 0;

   if ( ! $is_valid ) {
     WARN_USER( "\nInvalid STRIDE option.");
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Parse sparse option.

sub parse_sparse_option {
 my $value = shift;
 my $result = 0;

 if ( $sparse ne '' || $stride ne '' ) {
   WARN_USER( "\nRedundant SPARSE option.");
 } else {
   $sparse = $value;
   my $is_valid = $sparse >= 0;

   if ( ! $is_valid ) {
     WARN_USER( "\nInvalid SPARSE option.");
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Parse corners option:

sub parse_corners_option {
  my $value = shift;
  my $result = 0;

  if ( $corners ne '' ) {
    WARN_USER( "\nRedundant CORNERS option.");
  } else {
    $corners = $value;
    my $is_valid = $corners == 0 || $corners == 1;

    if ( ! $is_valid ) {
      WARN_USER( "\nInvalid CORNERS option.");
    } else {
      $result = 1;
    }
  }

  return $result;
}



# Parse compress option:

sub parse_compress_option {
 my $value = shift;
 my $result = 0;

 if ( $compress ne '' ) {
   WARN_USER( "\nRedundant COMPRESS option.");
 } else {
   $compress = $value;
   my $is_valid = $compress == 0 || $compress == 1;

   if ( ! $is_valid ) {
     WARN_USER( "\nInvalid COMPRESS option.");
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Parse regrid option:

sub parse_regrid_option {
 my $value = shift;
 my $result = parse_option( $regrid, $value, 'REGRID',
                            'nearest mean weighted' );
 return $result;
}



# Parse regrid_aggregate option:

sub parse_regrid_aggregate_option {
  my $value = shift;
  my $result = parse_option( $regrid_aggregate, $value, 'REGRID_AGGREGATE',
                             'none all daily' );
  return $result;
}



# Parse lambert option.

sub parse_lambert_option {
 my $value = shift;
 my $result = 0;

 if ( $lambert ne '' ) {
   WARN_USER( "\nRedundant LAMBERT option.");
 } else {
   $lambert = $value;
   my @values = split( /,/, $value );
   my $count = @values;

   if ( $count != 4 ) {
     WARN_USER( "\nInvalid LAMBERT option.");
   } elsif ( ! in_range( $values[ 0 ], -89.0, 89.0 ) ) {
     WARN_USER( "\nInvalid LAMBERT option 1.'");
   } elsif ( ! in_range( $values[ 1 ], -89.0, 89.0 ) ) {
     WARN_USER( "\nInvalid LAMBERT option 2.'");
   } elsif ( ! in_range( $values[ 2 ], -180.0, 180.0 ) ) {
     WARN_USER( "\nInvalid LAMBERT option 3.");
   } elsif ( ! in_range( $values[ 3 ], -89.0, 89.0 ) ) {
     WARN_USER( "\nInvalid LAMBERT option 4.");
   } elsif ( $values[ 0 ] > $values[ 1 ] ) {
     WARN_USER( "\nInvalid LAMBERT option 1.");
   } elsif ( $values[ 0 ] > 0.0 && $values[ 1 ] < 0.0 ) {
     WARN_USER( "\nInvalid LAMBERT option 1.");
   } elsif ( $values[ 0 ] < 0.0 && $values[ 1 ] > 0.0 ) {
     WARN_USER( "\nInvalid LAMBERT option 1.");
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Parse stereographic option.

sub parse_stereographic_option {
  my $value = shift;
  my $result = 0;

  if ( $stereographic ne '' ) {
    WARN_USER( "\nRedundant STEREOGRAPHIC option.");
  } else {
    $stereographic = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 3 ) {
      WARN_USER( "\nInvalid STEREOGRAPHIC option.");
    } elsif ( ! in_range( $values[ 0 ], -180.0, 180.0 ) ) {
      WARN_USER( "\nInvalid STEREOGRAPHIC option 1.");
    } elsif ( ! in_range( $values[ 1 ], -90.0, 90.0 ) ) {
      WARN_USER( "\nInvalid STEREOGRAPHIC option 2.");
    } elsif ( ! in_range( $values[ 2 ], -90.0, 90.0 ) ) {
      WARN_USER( "\nInvalid STEREOGRAPHIC option 3.");
    } else {
      $result = 1;
    }
  }

  return $result;
}



# Parse mercator option.

sub parse_mercator_option {
  my $value = shift;
  my $result = 0;

  if ( $mercator ne '' ) {
    WARN_USER( "\nRedundant MERCATOR option.");
  } else {
    $mercator = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 1 ) {
      WARN_USER( "\nInvalid MERCATOR option.");
    } elsif ( ! in_range( $values[ 0 ], -180.0, 180.0 ) ) {
      WARN_USER( "\nInvalid MERCATOR option 1.");
    } else {
      $result = 1;
    }
  }

  return $result;
}



# Parse lonlat option.

sub parse_lonlat_option {
  my $value = shift;
  my $result = 0;

  if ( $lonlat ne '' ) {
    WARN_USER( "\nRedundant LONLAT option.");
  } else {
    $lonlat = 1;
    $result = 1;
  }

  return $result;
}



# Parse ellipsoid option.

sub parse_ellipsoid_option {
 my $value = shift;
 my $result = 0;

 if ( $ellipsoid ne '' ) {
   WARN_USER( "\nRedundant ELLIPSOID option.");
 } else {
   $ellipsoid = $value;
   my @values = split( /,/, $value );
   my $count = @values;

   if ( $count != 1 && $count != 2 ) {
     WARN_USER( "\nInvalid ELLIPSOID option.");
   } elsif ( ! in_range( $values[ 0 ], 1.0, 1e10 ) ) {
     WARN_USER( "\nInvalid ELLIPSOID option 1.");
   } elsif ( $count == 2 ) {

     if ( ! in_range( $values[ 1 ], $values[ 0 ], 1e10 ) ) {
       WARN_USER( "\nInvalid ELLIPSOID option 2.");
     } else {
       $result = 1;
     }
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Parse grid option.

sub parse_grid_option {
 my $value = shift;
 my $result = 0;

 if ( $grid ne '' ) {
   WARN_USER( "\nRedundant GRID option.");
 } else {
   $grid = $value;
   my @values = split( /,/, $value );
   my $count = @values;

   if ( $count != 6 ) {
     WARN_USER( "\nInvalid GRID option.");
   } elsif ( $values[ 0 ] < 1 ) {
     WARN_USER( "\nInvalid GRID option 1.");
   } elsif ( $values[ 1 ] < 1 ) {
     WARN_USER( "\nInvalid GRID option 2.");
   } elsif ( ! in_range( $values[ 2 ], -1e10, 1e10 ) ) {
     WARN_USER( "\nInvalid GRID option 3.");
   } elsif ( ! in_range( $values[ 3 ], -1e10, 1e10 ) ) {
     WARN_USER( "\nInvalid GRID option 4.");
   } elsif ( ! in_range( $values[ 4 ], 0.01, 1e10 ) ) {
     WARN_USER( "\nInvalid GRID option 5.");
   } elsif ( ! in_range( $values[ 5 ], 0.01, 1e10 ) ) {
     WARN_USER( "\nInvalid GRID option 6.");
   } else {
     $result = 1;
   }
 }

 return $result;
}



# Check that all required 'options' have been specified:

sub required_options_specified {
  my $result = 0;

  if ( $service eq '' ) {
    WARN_USER( "\nMissing option: 'SERVICE'");
  } elsif ( $version eq '' ) {
    WARN_USER( "\nMissing option: 'VERSION'");
  } elsif ( $request eq '' ) {
    WARN_USER( "\nMissing option: 'REQUEST'");
  } elsif ( $request eq 'getcoverage' ) {

    if ( $variable eq '' ) {
      WARN_USER( "\nMissing option: 'COVERAGE'");
    } elsif ( $format eq '' ) {
      WARN_USER( "\nMissing option: 'FORMAT'");
    } elsif ( $time eq '' ) {
      WARN_USER( "\nMissing option: 'TIME'");
    } elsif ( $bbox eq '' ) {
      WARN_USER( "\nMissing option: 'BBOX'");
    } else {
      my $regrid_count = 0;
      $regrid_count += $regrid ne '';
      $regrid_count += $ellipsoid ne '';
      $regrid_count += $grid ne '';
      my $projection_count = 0;
      $projection_count += $lambert ne '';
      $projection_count += $stereographic ne '';
      $projection_count += $mercator ne '';
      $projection_count += $lonlat ne '';

      $result =
        $regrid_count == 0 && $projection_count == 0 ||
        $regrid_count == 3 && $projection_count == 1;

      if ( ! $result ) {
        WARN_USER( "\nInvalid options: 'REGRID/");
        WARN_USER( "LAMBERT/STEREOGRAPHIC/MERCATOR/LONLAT/");
        WARN_USER( "GRID/ELLIPSOID'");
      }
    }
  } else {
    $result = 1;
  }

  return $result;
}



# Compute starting_timestamp and hours.
# inputs:  $time = '2001-08-26t00:00:00z/2001-08-30t23:59:59z'
# outputs: $starting_timestamp = 2005082600
#          $hours = 120

sub compute_time_range {
 my $yyyy1 = substr( $time, 0, 4 );
 my $mm1   = substr( $time, 5, 2 );
 my $dd1   = substr( $time, 8, 2 );
 my $hh1   = substr( $time, 11, 2 );
 my $i     = index( $time, '/' );

 $starting_timestamp = integer_timestamp( $yyyy1, $mm1, $dd1, $hh1 );
 $hours = 1;

 if ( $i != -1 ) {
   ++$i;
   my $yyyy2 = substr( $time, $i + 0, 4 );
   my $mm2   = substr( $time, $i + 5, 2 );
   my $dd2   = substr( $time, $i + 8, 2 );
   my $hh2   = substr( $time, $i + 11, 2 );
   my $yyyy  = $yyyy1;
   my $mm    = $mm1;
   my $dd    = $dd1;
   my $hh    = $hh1;

   while ( integer_timestamp( $yyyy, $mm, $dd, $hh ) !=
           integer_timestamp( $yyyy2, $mm2, $dd2, $hh2 ) ) {
     increment_timestamp( $yyyy, $mm, $dd, $hh );
     ++$hours;
   }
 }
}



# Construct command.

sub construct_command {
  my @bounds = split( /,/, $bbox );
  my $domain = " -domain $bounds[ 0 ] $bounds[ 1 ] $bounds[ 2 ] $bounds[ 3 ] ";

  if ( $format eq 'original' ) {
    $command = "$getfiles $list_file";
  } else {
    my $my_xdrconvert = '';
    my $my_compressor = '';
    my $corners_option = '';
    my $subsetter_format = $format;
    my @format_parts = split( /-/, $format );
    my $format_parts_count = @format_parts;

    if ( $format_parts_count == 2 ) {
      $subsetter_format = $format_parts[ 1 ];
    }

    if ( $format ne 'xdr' || $regrid ne '' || $subsetter_format eq 'coards' ) {
      my $regrid_args = '';
      my $xdrconvert_format = $subsetter_format;

      if ( $regrid ne '' ) {
        my $projection_args =
          $lambert ne '' ? "-lambert $lambert "
          : $stereographic ne '' ? "-stereographic $stereographic "
          : $mercator ne '' ? "-mercator $mercator "
          : "-lonlat ";

        $projection_args =~ tr/,/ /;
        my @ellipsoid_args = split( /,/, $ellipsoid );
        my $ellipsoid_args_count = @ellipsoid_args;
        my $major_semiaxis = $ellipsoid_args[ 0 ];
        my $minor_semiaxis =
          $ellipsoid_args_count == 1 ? $major_semiaxis : $ellipsoid_args[ 1 ];
        my $grid_args = "-grid $grid ";
        $grid_args =~ tr/,/ /;

        my $regrid_aggregate_option = '';

        if ( $regrid_aggregate eq 'daily' ) {
          $regrid_aggregate_option = '-aggregate 24';
        } elsif ( $regrid_aggregate eq 'all' ) {
          $regrid_aggregate_option = "-aggregate $hours";
        }

        $regrid_args =
          "-regrid $regrid " .
          $projection_args .
          "-ellipsoid $major_semiaxis $minor_semiaxis " .
          $grid_args .
          $regrid_aggregate_option;
      }

      $my_xdrconvert = " | $xdrconvert -tmpdir $temp_dir $regrid_args -$xdrconvert_format";
    }

    if ( $compress ne '' && $compress == 1 ) {
      $my_compressor = " | $compressor";
    }

    if ( $corners ne '' && $corners == 1 && $regrid ne 'nearest' ) {
      $corners_option = " -corners";
    }

    $command =
      "$subsetter" .
      " -files $list_file" .
      " -tmpdir $temp_dir" .
      " -desc https://modis.gsfc.nasa.gov/data/,MODISSubset" .
      " -timestamp $starting_timestamp -hours $hours" .
      " -variable $variable" .
      $domain .
      "$corners_option" .
      "$my_xdrconvert$my_compressor";
  }
  return $command;
}



# For REQUEST=GetMetadata, first call construct_command() to get the text of
# the command that will create the actual coverage, then call this
# function to generate a command that will print the metadata, the command
# string from construct_command, and the list of files that would be processed
# by that command.
#
# see run_RSIG_Backend::run_RSIG() for the logic for making this all go.

sub print_metadata {
  my ($c) = @_;
  # clean up $c to remove local path info
  my $remove = "$RSIG_Common::MODAPSROOT/";
  $c =~ s|$remove||g;
  my $result = '';

  print $metadata_content;

  if ( $format eq 'original' ) {
    print $metadata_content_original;
    $result = "/bin/cat $list_file";
  } else {
    print $metadata_content_processed;
    $result = "/bin/cat $list_file ; /bin/echo ; /bin/echo '$c' ; /bin/echo";
    print $metadata_content_original;
  }

  return $result;
}



############################### HELPER ROUTINES ##############################



# debug( message );

sub debug {
 my $message = shift;

 if ( $debugging ) {
   LOG( "\n$message\n");
####print STDOUT "\n$message\n";
 }
}

# my $result = parse_option( $option, $value, $option_name, $valid_values );
# my $result = parse_option( $variable, $value, 'COVERAGE', 'ozone pm25' );

sub parse_option {
 my ( $option, $value, $option_name, $valid_values ) = @_;
 my $result = 0;

 if ( $option ne '' ) {
   WARN_USER( "\nRedundant $option_name option.");
 } else {
   $result = index( " $valid_values ", " $value " ) != -1;

   if ( $result ) {
     $_[ 0 ] = $value;
   } else {
     WARN_USER( "\nInvalid $option_name option. Permitted values are: $valid_values");
   }
 }

 return $result;
}



# my $ok = in_range( $value, $minimum, $maximum );

sub in_range {
 my ( $value, $minimum, $maximum ) = @_;
 my $result = $value >= $minimum && $value <= $maximum;
 return $result;
}



# my $is_valid = is_valid_time( '2001-08-26t20:00:00z/2001-08-27t23:59:59z' );

sub is_valid_time {
 my $time = shift;
 my $result = 0;
 my $length = length( $time );

 if ( $length == 41 ) {
   $result = is_valid_time( substr( $time, 0, 20 ) );
   $result = $result && substr( $time, 20, 1 ) eq '/';
   $result = $result && is_valid_time( substr( $time, 21, 20 ) );
   $result = $result && substr( $time, 0, 20 ) le substr( $time, 21, 20 );
 } elsif ( $length == 20 ) {
   my $year   = substr( $time, 0, 4 );
   my $month  = substr( $time, 5, 2 );
   my $day    = substr( $time, 8, 2 );
   my $hour   = substr( $time, 11, 2 );
   my $minute = substr( $time, 14, 2 );
   my $second = substr( $time, 17, 2 );
   $result = in_range( $year, 1900, 3000 );
   $result = $result && in_range( $month, 1, 12 );
   $result = $result && in_range( $day, 1, days_in_month( $year, $month ) );
   $result = $result && in_range( $hour, 0, 23 );
   $result = $result && in_range( $minute, 0, 59 );
   $result = $result && in_range( $second, 0, 59 );
   $result = $result && substr( $time, 4, 1 ) eq '-';
   $result = $result && substr( $time, 7, 1 ) eq '-';
   $result = $result && substr( $time, 10, 1 ) eq 't';
   $result = $result && substr( $time, 13, 1 ) eq ':';
   $result = $result && substr( $time, 16, 1 ) eq ':';
   $result = $result && substr( $time, 19, 1 ) eq 'z';
 }

 return $result;
}



# increment_timestamp( $yyyy, $mm, $dd, $hh );

sub increment_timestamp {
 my ( $yyyy, $mm, $dd, $hh ) = @_;
 my $hours_per_day   = 23;
 my $months_per_year = 12;
 my $days_this_month = days_in_month( $yyyy, $mm );
 ++$hh;

 if ( $hh > $hours_per_day ) {
   $hh = 0;
   ++$dd;

   if ( $dd > $days_this_month ) {
     $dd = 1;
     ++$mm;

     if ( $mm > $months_per_year ) {
       $mm = 1;
       ++$yyyy;
     }
   }
 }

 $_[ 0 ] = $yyyy;
 $_[ 1 ] = $mm;
 $_[ 2 ] = $dd;
 $_[ 3 ] = $hh;
}



# my $yyyymmddhh = integer_timestamp( $yyyy, $mm, $dd, $hh );

sub integer_timestamp {
 my ( $yyyy, $mm, $dd, $hh ) = @_;
 my $result = int( $yyyy * 1000000 + $mm * 10000 + $dd * 100 + $hh );
 return $result;
}



# my $leap = is_leap_year( $year );

sub is_leap_year {
 my $year = shift;
 my $result = $year % 4 == 0 && ! ( $year % 100 == 0 && ! $year % 400 == 0 );
 return $result;
}



# my $days = days_in_month( $year, $month );

sub days_in_month {
 my $year  = shift;
 my $month = shift;

 # 30 days hath September, April, June and November...

 my @days_per_month = ( 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 );
 my $result = $days_per_month[ $month - 1 ];

 if ( $month == 2 ) {

   if ( is_leap_year( $year ) ) {
     ++$result;
   }
 }

 return $result;
}



# SQL Query and get data locations
sub getLogin
{
    my ($name, $passwd, $separator) = @_;
    return undef unless defined $name;

    $passwd ||= $RSIG_Common::MODAPSROOT.'/etc/passwd';
    $separator ||= '\s+';

    my $pwfile = new FileHandle($passwd) or do {
        warn "E Failed open configuration file $passwd";
	return undef;
    };

    while (<$pwfile>)
    {
	if (/^$name\s+/) {
	    my @param = split /$separator/, $_;
	    shift @param;
	    return wantarray ? @param : \@param;
	}
    }
    return undef;
}

sub get_modis_data
{
  my ($work_dir, $alldata_dir) = @_;
  die "no work_dir" unless $work_dir;
  die "no alldata_dir" unless $alldata_dir;

  use DBI;
  use POSIX qw(floor ceil);

  # set so we can use rest of MTVS stuff...
  $ENV{MODAPSROOT} = $RSIG_Common::MODAPSROOT;

  # setup parameters for sql query...
  my ($west,$south,$east,$north) = split( /,/, $bbox );

  my $esdt_number = substr( $mod, 3 );
  my $esdt_mod   =  "MOD0${esdt_number}_L2";
  my $esdt_myd   =  "MYD0${esdt_number}_L2";
  if ($esdt_number =~ /43k/)
  {
    $esdt_mod   =  "MOD04_3K";
    $esdt_myd   =  "MYD04_3K";
  }
  my $collection =  61;

  #....TIME=2005-08-26T00:00:00Z/2005-08-30T23:59:59Z
  my $user_input_time = length( $time ) == 20 ? "$time/$time" : $time;

  #  Create the start and end time strings for the SQL command
  my ($starttime,$endtime) = split(/\//, $user_input_time);
  $starttime =~ s/T/ /i;
  $starttime =~ s/Z//i;

  $endtime =~ s/T/ /i;
  $endtime =~ s/Z//i;

  # ready for database connection...
  my ($driver, $server, $port, $dbname, $user, $password) =
    getLogin("lads_www");
  die "unknown driver $driver" unless $driver eq "Pg";

  my $dbh = DBI->connect( 
    "DBI:Pg(Taint=>1):dbname=$dbname;host=$server;port=$port", 
    $user, $password, { AutoCommit => 1, RaiseError => 1 } 
  ) or die("Failed to open $dbname database!");

  if ( $dbh->{'pg_server_version'} !~ /^8/ )
  {
    $dbh->do( qq{set application_name = 'RSIG'} );
  }

  # search for files...
  my $esdt_list = [$esdt_mod, $esdt_myd];
  my $constraints = {
    ARCHIVE_SET => $collection,
    START_TIME  => $starttime,
    END_TIME    => $endtime,
    NORTH => $north,
    SOUTH => $south,
    EAST => $east,
    WEST => $west,
  };
  my $results = search_metadata($dbh, $esdt_list, $constraints);

  # make list_file name unique so multiple users don't write over each other
  $list_file = "$work_dir/list.$RSIG_Common::JOB_ID";
  $temp_dir = $work_dir; # MODISSubset will create and remove temp files here.

  open (FILEHANDLE, ">$list_file") or die ("cannot open $list_file ");
  # process rows...
  foreach my $filename (@$results)
  {
    my ($esdt,$yd,$gtime,$as) = split(/\./, $filename);
    my ($yr,$doy) = ($1,$2) if $yd =~ /A(\d{4})(\d{3})/;
    print FILEHANDLE "$alldata_dir/$collection/$esdt/$yr/$doy/$filename\n";
  }
  close FILEHANDLE;
}

#====================================================================
my %types = (
    GEO  => "FileMeta_Geolocation",
    L1L2 => "FileMeta_L1L2",
    TILE => "FileMeta_L3L4Tiled",
);

sub get_meta_tables
{
  my ($dbh, $args) = @_;
  die "no args" unless $args;

  my $products = $args->{Products};

  my ($sql, @bind, %tables);

  #
  # Query the ESDT_Def table for the satellite/instrument and table type
  #

  if ( @$products )
  {
      $sql = join("\n",
        qq{select "SatelliteInstrument", "FileSubTableType", "ESDT"},
        qq{from "ESDT_Def"},
        qq{where "ESDT" = ?},
      );
      @bind = @$products;
  }
  else
  {
      return [];
  }

  # Assign the satellite/instrument and table to the hash
  my $sth = $dbh->prepare( $sql );
  foreach my $item ( @bind )
  {
      $sth->execute( $item );
      while ( my ($si, $type, $value) = $sth->fetchrow_array )
      {
          push @{ $tables{$si}{$types{$type}} }, $value;
      }
  }
  return %tables;
}
#====================================================================
# search for files as specified in the $constraints hash. Only
# return files that meet the specified constraints.
#
my %coverageTables = (
    "FileMeta_L1L2"        => 1,
    "FileMeta_Geolocation" => 1,
);

my %boxTables = (
    "FileMeta_L1L2"        => 1,
    "FileMeta_Geolocation" => 1,
);

my %overlapTables = (
    "FileMeta_L0"          => 1,
    "FileMeta_L1L2"        => 1,
);

my %spatialTables = (
    "FileMeta_Geolocation" => 1,
);

my %tileTables = (
    "FileMeta_L3L4Tiled"   => 1,
);

#====================================================================
# NOTE: $instance is necessary because the sql query depends
# on the instance. Not all instances use tables the same way.
sub search_metadata
{
  my ($dbh, $esdt_list, $constraints) = @_;

  # constraints...
  #??? Should we really default START_TIME/END_TIME if supplied values invalid?
  my $startTime  = $constraints->{START_TIME};
  die "invalid date/time : $constraints->{START_TIME}" unless $startTime;
  my $endTime    = $constraints->{END_TIME};
  die "invalid date/time : $constraints->{END_TIME}" unless $endTime;
  my $north = $constraints->{NORTH};
  my $east = $constraints->{EAST};
  my $south = $constraints->{SOUTH};
  my $west = $constraints->{WEST};
  my $polar_tiles = 0;
  my $archiveSet     = $constraints->{ARCHIVE_SET};
  my @coverageOptions  = ();
  my $granuleLimit   = $constraints->{MAX_FILES_PER_ORDER};
  $granuleLimit = $GRANULE_LIMIT unless $granuleLimit;

  my $count = 0;
  my $fileNames = [];

  #
  # Create a hash of dates
  #

  my %allDates = ();
  $allDates{$startTime}{START} = $startTime;
  $allDates{$startTime}{END}   = $endTime;

  # Determine if the spatial coordinates fit in a single box tile
  # NOTE: THIS IS ONLY APPLICABLE FOR LAADS INSTANCE, NOT OTHER INSTANCES.
  my $boxTile;
  if ( floor($south/10) == floor($north/10)
  and  floor($west/10)  == floor($east/10) )
  {
    $boxTile = sprintf("%02d0%02d", floor($south/10)+10, floor($west/10)+19);
  }

#
# Determine the metadata tables for all products
#

  $esdt_list = [$esdt_list] unless ref($esdt_list) eq 'ARRAY';
  my %tables = get_meta_tables( $dbh, {Products => $esdt_list } );

#
# For each satellite/instrument and table, search the table for all granules
# that match the defined search parameters
#

foreach my $si ( keys %tables ) {

    foreach my $table ( keys %{ $tables{$si} } ) {

        my @bind;

        my $geoESDT =
            ( $si eq "AM1M" ) ? "MOD03" :
            ( $si eq "PM1M" ) ? "MYD03" :
            ( $si eq "NPP"  ) ? "NPP_VMAE_L1" : undef;

        #
        # Use the spatial coordinates if the coordinates are defined and
        # do not include the whole earth and if no tiles are defined or
        # the products belonging to this table are not tiled products
        #

        my $spatial = 0;

        $spatial = 1
            if ( defined $north and defined $south
                 and defined $east and defined $west
                 and not ( $north ==  90 and $south ==  -90 and
		       $east  == 180 and $west  == -180
                 )
            );

        my $sql = qq{select "FileName", fma."FileId", fma."PGE_StartTime", fma."ESDT"};

	$sql .= qq{, fmg."GRingLongitude1", fmg."GRingLongitude2"}.
                qq{, fmg."GRingLongitude3", fmg."GRingLongitude4"}.
                qq{, fmg."GRingLatitude1",  fmg."GRingLatitude2"}.
                qq{, fmg."GRingLatitude3",  fmg."GRingLatitude4"}
            if ( $spatial and $overlapTables{$table} );

	$sql .= qq{, fma."GRingLongitude1", fma."GRingLongitude2"}.
                qq{, fma."GRingLongitude3", fma."GRingLongitude4"}.
                qq{, fma."GRingLatitude1",  fma."GRingLatitude2"}.
                qq{, fma."GRingLatitude3",  fma."GRingLatitude4"}
            if ( $spatial and $table eq "FileMeta_Geolocation" );

        $sql .= qq{ from "File" f, "$table" fma};

	$sql .= qq{, "Tile_GeoStartTime" tgst}
	    if ( $boxTile and $boxTables{$table} and $si ne "NPP" and $si ne "MERIS" );

	$sql .= qq{, "FileMeta_Geolocation" fmg}
	    if ( $spatial and $overlapTables{$table} );

	$sql .= qq{, "Tile_Def" td}
	    if ( $spatial and $tileTables{$table} );

	$sql .= qq{ where f."FileId"=fma."FileId" and};
	
	    if ( $boxTile and $boxTables{$table} and $si ne "NPP" and $si ne "MERIS" )
            {
              $sql .= qq{ tgst."SatelliteInstrument" = ?}.
	        qq{ and tgst."BoxTile" = ?}.
		qq{ and tgst."PGE_StartTime" = fma."PGE_StartTime" and};
              push @bind, $si, $boxTile;
            }

            if ( $spatial and $overlapTables{$table} )
            {
              $sql .= qq{ fmg."ESDT" = ?}.
                qq{ and fmg."PGE_StartTime" = fma."PGE_StartTime"}.
                qq{ and fmg."ArchiveSet" = fma."ArchiveSet" and};

              push @bind, $geoESDT;
            }

        $sql .= qq{ fma."Tile" = td."Tile" and}
            if ( $spatial and $tileTables{$table} );

      $sql = join("\n", $sql,
        qq{not exists (select 1 from "FilesToPhysSync" where "FileId" = fma."FileId")},
        qq{and fma."ESDT" in (?} . ", ?" x $#{$tables{$si}{$table}} . qq{)},
      );
      push @bind, ( @{ $tables{$si}{$table} } );

      if ($archiveSet)
      {
        $sql = join("\n", $sql,
          qq{and fma."ArchiveSet" = ?},
        );
        push @bind, ( $archiveSet );
      }

        if ( $spatial and ( $overlapTables{$table} or
	                    $spatialTables{$table} or
                            $tileTables{$table} ) )
        {

            my ($N, $S, $E, $W);

            if ( $tileTables{$table} ) {

                ($N, $S, $E, $W) = (
                    "NorthBoundCoord",
                    "SouthBoundCoord",
                    "EastBoundCoord",
                    "WestBoundCoord",
                );

            } else {

                ($N, $S, $E, $W) = (
                    "NorthBoundingCoord",
                    "SouthBoundingCoord",
                    "EastBoundingCoord",
                    "WestBoundingCoord",
                );

            } #endif

            $sql .= qq{ and "$N" >= ? and "$S" <= ?};

            push @bind, $south, $north;

            if ( $east >= $west ) {

                $sql .= qq{ and ( ( ( "$E" > "$W" and "$W" < ? and "$E" > ? ) or}
                               .qq{ ( "$E" < "$W" and ( "$W" < ? or "$E" > ? ) ) )}
                             .qq{ or "$E" = "$W" )};

                push @bind, $east, $west, $east, $west;

            } elsif ( $east < $west ) {

                $sql .= qq{ and ( ( ( "$E" > "$W" and ( "$W" < ? or "$E" > ? ) ) or}
                               .qq{ "$E" < "$W" ) or "$E" = "$W" )};

                push @bind, $east, $west;

            }#endif

        } #endif

        if ( @coverageOptions and $coverageTables{$table} ) {

	    $sql .= qq{ and fma."DayNightFlag" in (?} . ", ?" x $#coverageOptions . qq{)};

	    push @bind, ( @coverageOptions );

	} #endif


        if ( $boxTile and $boxTables{$table} and $si ne "NPP" and $si ne "MERIS" ) {

            $sql .= qq{ and tgst."PGE_StartTime" between ? and ?};

        } else {

            $sql .= qq{ and fma."PGE_StartTime" between ? and ?};

	} #endif

        my $sth = $dbh->prepare( $sql );

        foreach my $date ( sort keys %allDates ) {

            $sth->execute( @bind, $allDates{$date}{START}, $allDates{$date}{END} );

            while ( my ($file_name, $fileID, $time, $esdt, @gring) = $sth->fetchrow_array ) {

                push @$fileNames, $file_name;
                ++$count;

                #
	        # Stop searching if the granule limit is reached
	        #

                if ( defined $granuleLimit and $count >= $granuleLimit ) {

                  print "  Search results exceeded the $granuleLimit granule",
		        "  limit.\n  Only the first $granuleLimit granules",
		        "  are displayed.\n";

                   $sth->finish;

	           last;

                } #endif

            } #endwhile

            last if ( defined $granuleLimit and $count >= $granuleLimit );

        } #endforeach

       last if ( defined $granuleLimit and $count >= $granuleLimit );

    } #endforeach

    last if ( defined $granuleLimit and $count >= $granuleLimit );

} #endforeach

# EPA wants files sorted by PGE_StartTime
return [sort by_pge_starttime @$fileNames];
}

#------------------------------------------------------------------------------
1;  # end of package RSIG_Common

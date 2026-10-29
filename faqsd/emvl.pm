#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: emvl.pm - Utilities to parse WCS URL used by RSIG *server files.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use emvl; # Local utility functions.
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-03-09 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package emvl;
$| = 1; # Turn off output buffering so messages appear in correct order.
# Restrict PATH and ENV. (But before popen, ENV is cleared!)
delete @ENV{ 'PATH', 'IFS', 'CDPATH', 'ENV', 'BASH_ENV' };

##############################################################################

my $debugging = 0;
sub get_debugging { return $debugging; }
sub set_debugging { $debugging = shift; $debugging = int( $debugging ); }
sub debug { if ( $debugging ) { my @args = @_; print STDERR "\n@args\n"; } }

my $webservice_url = 'https://ofmpub.epa.gov/rsig/rsigserver?';
my $organization_url = 'https://www.epa.gov/rsig';
my $email_address ='rsig@epa.gov';


################################## VARIABLES #################################

# Parsed from the URL query string:

my $parsed_query_string = '';

my $service       = ''; # wcs.
my $version       = ''; # 1.0.0.
my $request       = ''; # getcapabilities or describecoverage or getcoverage.
my $coverage      = ''; # E.g., aqs.pm25.
my $time          = ''; # E.g., 2008-06-21t00:00:00z/2008-06-22t23:59:59z.
my $bbox          = ''; # E.g., -76,34,-74,36,0,0.
my $format        = ''; # xdr, ascii, netcdf-coards, netcdf-ioapi, original.
my $compress      = ''; # 1 = | gzip -c otherwise don't compress (default).
my $regrid        = ''; # E.g., nearest, mean, weighted.
my $regrid_aggregate = ''; # E.g., none, all, daily.
my $lambert       = ''; # E.g., 33,45,-97,40.
my $stereographic = ''; # E.g., -98,90,45.
my $mercator      = ''; # E.g., -98.
my $lonlat        = '';
my $ellipsoid     = ''; # E.g., 6370000,6370000.
my $grid          = ''; # E.g., 279,240,-1008000,-1620000,12000,12000.
my $levels        = ''; # E.g,22,2,10000,1.0,...,0.0,9.81,287.04,50,290,100000.

# Derived from the above parsed values and initialized during parsing:

my $coverage_metadata_reference = undef; # Optionally set_coverage_metadata().
my $coverage_index = -1; # Index into coverage_metadata array of COVERAGE=value.
my $yyyymmddhhmmss1 = '';
my $yyyymmddhhmmss2 = '';
my $yyyymmdd1 = 0;
my $yyyymmdd2 = 0;
my $days = 0;
my $hours = 0;
my @bounds = ();

# These variables can be set as defaults used in print_capabilities() and
# print_coverage_description() if they are not included in coverage_metadata
# entries.

my $metadata_begin_date = '';
my $metadata_end_date = '';
my $metadata_west_bound = '';
my $metadata_east_bound = '';
my $metadata_south_bound = '';
my $metadata_north_bound = '';
sub set_metadata_begin_date { $metadata_begin_date = shift; }
sub set_metadata_end_date { $metadata_end_date = shift; }
sub set_metadata_west_bound { $metadata_west_bound = shift; }
sub set_metadata_east_bound { $metadata_east_bound = shift; }
sub set_metadata_south_bound { $metadata_south_bound = shift; }
sub set_metadata_north_bound { $metadata_north_bound = shift; }

################################## ROUTINES ##################################

# Get the global variables:

sub get_parsed_query_string { return $parsed_query_string; }
sub get_service { return $service; }
sub get_version { return $version; }
sub get_request { return $request; }
sub get_coverage { return $coverage; }
sub get_coverage_index { return $coverage_index; }
sub get_bbox { return $bbox; }
sub get_bounds { return @bounds; }
sub get_west {  return $bounds[ 0 ]; }
sub get_south { return $bounds[ 1 ]; }
sub get_east {  return $bounds[ 2 ]; }
sub get_north { return $bounds[ 3 ]; }
sub get_west_south_east_north {
  return "$bounds[ 0 ] $bounds[ 1 ] $bounds[ 2 ] $bounds[ 3 ]";
}
sub get_lower { my $count = @bounds; return $count == 6 ? $bounds[ 4 ] : ''; }
sub get_upper { my $count = @bounds; return $count == 6 ? $bounds[ 5 ] : ''; }
sub get_format { return $format; }
sub get_time { return $time; }
sub get_yyyymmddhhmmss1 { return $yyyymmddhhmmss1; }
sub get_yyyymmddhhmmss2 { return $yyyymmddhhmmss2; }
sub get_yyyymmddhh1 { return substr( $yyyymmddhhmmss1, 0, 10 ); }
sub get_yyyymmddhh2 { return substr( $yyyymmddhhmmss2, 0, 10 ); }
sub get_yyyymmdd1 { return $yyyymmdd1; }
sub get_yyyymmdd2 { return $yyyymmdd2; }
sub get_days { return $days; }
sub get_hours { return $hours; }
sub get_compress { return $compress eq '1' ? 1 : 0; }
sub get_regrid { return $regrid; }
sub get_regrid_aggregate { return $regrid_aggregate; }
sub get_lambert { return $lambert; }
sub get_stereographic { return $stereographic; }
sub get_mercator { return $mercator; }
sub get_lonlat { return $lonlat; }
sub get_ellipsoid { return $ellipsoid; }
sub get_grid { return $grid; }
sub get_levels { return $levels; }

# Query string parsing routine dispatch table:

my %parsers = (
  'service'          => \&parse_service_option,
  'version'          => \&parse_version_option,
  'request'          => \&parse_request_option,
  'coverage'         => \&parse_coverage_option,
  'time'             => \&parse_time_option,
  'bbox'             => \&parse_bbox_option,
  'format'           => \&parse_format_option,
  'compress'         => \&parse_compress_option,
  'regrid'           => \&parse_regrid_option,
  'regrid_aggregate' => \&parse_regrid_aggregate_option,
  'lambert'          => \&parse_lambert_option,
  'stereographic'    => \&parse_stereographic_option,
  'mercator'         => \&parse_mercator_option,
  'lonlat'           => \&parse_lonlat_option,
  'ellipsoid'        => \&parse_ellipsoid_option,
  'grid'             => \&parse_grid_option,
  'levels'           => \&parse_levels_option,
  'rsig_id'          => \&parse_rsigid_option
);

# Before calling emvl::parse_url(),
# Callers can append extra parsers (or replace default parsers). For example,
# emvl::append_parser( 'aggregate', \&parse_aggregate_option );

sub append_parser {
  my ($option_name, $parser) = @_;
  $parsers{ $option_name } = $parser;
}


# Webserver content types for each output format:

my %content_types = (
  'ascii'         => 'text/plain',
  'xdr'           => 'application/octet-stream',
  'netcdf-coards' => 'application/netcdf',
  'netcdf-ioapi'  => 'application/netcdf',
  'original'      => 'application/octet-stream'
);

# Before calling emvl::parse_url(),
# callers can append extra format and content types, for example,
# emvl::append_content_type( 'shapefile', 'application/zip' );
# emvl::append_content_type( 'geojson', 'application/geo_json' );

sub append_content_type {
  my ($format_name, $content_type) = @_;
  $content_types{ $format_name } = $content_type;
}

sub get_formats { return join( " ", sort keys %content_types ); }
sub uniq { my %seen; return grep { !$seen{$_}++ } @_; }
sub get_content_types { return join( " ",  uniq( sort values %content_types ) ); }
sub get_content_type_header {
  my $content_type = $content_types{ $format };
  return "Content-type: $content_type; charset=iso-8859-1";
}


# Before calling parse_urls(),
# caller can pass a reference to an array of coverage metadata structures
# with required key-value pairs shown in the example below.
# Use end_date => '', for data that exists through the present day.
# Other key-value pairs can be included (and are ignored by emvl.pm code).
#
#my @coverage_metadata  = (
#  {
#    name => 'NO2',
#    units => 'ppb',
#    description =>
#      'UTC sub-hourly DISCOVER-AQ ACAM aircraft-measured ' .
#      'nitrogen dioxide concentration.',
#    begin_date => '2011-07-01',
#    end_date   => '2011-07-29',
#    west_bound => -79.6,
#    east_bound => -72.3,
#    south_bound => 36.4,
#    north_bound => 40.7
#  },
#  {
#    name => 'VCD',
#    units => 'mol/cm2',
#    description =>
#      'UTC sub-hourly DISCOVER-AQ ACAM aircraft-measured ' .
#      'vertical column density.',
#    begin_date => '2011-07-01',
#    end_date   => '2011-07-29', # Or use '' for no end date for recurring data.
#    west_bound => -79.6,
#    east_bound => -72.3,
#    south_bound => 36.4,
#    north_bound => 40.7
#  },
#);
#
# emvl::set_coverage_metadata( \@coverage_metadata );
# Then the default parse_coverage(), print_capabilities() and
# print_describe_coverage() will work with this metadata information.

sub set_coverage_metadata { ($coverage_metadata_reference) = @_; }

# Before calling emvl::parse_url()
# callers can either specify a custom parser callback (Template Design Pattern)
# (that checks and returns coverage or '' if failed)
# like this: emvl::set_parse_coverage_option( \&parse_coverage_option );
# Or else specify a coverages array to search like this:
# emvl::set_coverages( \@variables );
# Or else specify a coverage_metadata array as shown above.
# Otherwise $coverage will just be initialized to whatever was specified in
# COVERAGE= parameter and emvl::get_coverage_index() will return -1;

my $custom_parse_coverage_option = undef;
sub get_parse_coverage_option { return $custom_parse_coverage_option; }
sub set_parse_coverage_option { $custom_parse_coverage_option = shift; }
my $coverages_reference = undef;
sub set_coverages { $coverages_reference = @_; }


# Before calling emvl::parse_url()
# callers can set the name of the source (e.g., CMAQ, AQS, etc.) to be used in
# pprint_capabilities().

my $source = '';
sub get_source { return $source; }
sub set_source { $source = shift; }

# Before calling emvl::parse_url()
# callers can specify a custom required_options_specified() function that will
# be called at the end of parsing instead of the default one.
# emvl::set_required_options_specified( \&required_options_specified );

my $custom_required_options_specified = undef;
sub set_required_options_specified {$custom_required_options_specified = shift;}


# Parse OGC-WCS query string URL parameters into global variables.

sub parse_url {
  use CGI qw/ untaint /;     # Load Safe untainting CGI routines.
  $CGI::DISABLE_UPLOADS = 1; # Disable uploads.
  $CGI::POST_MAX = 1024;     # Limit posts to 1024 bytes.
  my $query = new CGI;       # Parse QUERY_STRING.
  %ENV = ();                 # Unset all environment variables.
  my @names = $query->param; # Extract names of parameters.
  my $count = @names;
  my $result = 0;
  $parsed_query_string = '';
  debug( "@names\n$count" );

  for ( my $parameter = 0; $parameter < $count; ++$parameter) {
    $_ = $names[ $parameter ];     # $_ is tainted.
    s/[^A-Za-z1]/_/go;             # Untaint $_.
    my $option = $_;               # $option is untainted.
    $_ = $query->param( $option ); # $_ is tainted.
    s/[^\w\-.,:\/]/_/go;           # Untaint $_.
    my $value = $_;                # $value is untainted.
    my $lowercase_option = lc( $option );
    my $lowercase_value  = lc( $value );
    debug( "$lowercase_option $lowercase_value" );

    if ( $parsers{ $lowercase_option } ) {
      $result = $parsers{ $lowercase_option }->( $lowercase_value );

      if ( $result ) {

        if ( $parsed_query_string ne '' ) {
          $parsed_query_string .= '&';
        }

        $parsed_query_string .= $option . '=' . $value;
      }
    } else {
      print STDERR "\nInvalid option #$parameter\n";
      $result = 0;
    }

    if ( ! $result ) {
      $parameter = $count; # Finish looping.
    }
  }

  if ( $result ) {

    if ( defined( $custom_required_options_specified ) ) {
      $result = $custom_required_options_specified->();
    } else {
      $result = required_options_specified();
    }
  }

  return $result;
}



# Remove named option from parsed_query_string:
# my $edited_parsed_query_string =
#  emvl::remove_option_from_parsed_query_string( 'RSIG_ID' );


sub remove_option_from_parsed_query_string {
  my $option = shift;
  my $option_length = length( $option );
  my $result = '';
  my @parts = split( /&/, $parsed_query_string );

  # Append all option=value pairs except $option=value:

  foreach my $option_value_pair ( @parts ) {

    if ( substr( uc( $option_value_pair ), 0, $option_length ) ne $option ) {

      if ( $result ne '' ) {
        $result .= '&';
      }

      $result .= $option_value_pair;
    }
  }

  return $result;
}



# Print web server capabilities metadata.
# Caller may optionally pass a callback function that prints coverage offering
# brief for each variable as shown below.
# print_capabilities( \&print_coverage_offering_brief );
# Template Design Pattern.

sub print_capabilities {
  my $callback = shift;

  print '<?xml version="1.0" encoding="UTF-8" ?>' . "\n";
  print '<WCS_Capabilities version="1.0.0" ';
  print 'xmlns="https://www.opengeospatial.org/standards/wcs" ';
  print 'xmlns:gml="https://www.opengis.net//gml" ';
  print 'xmlns:xlink="https://www.w3.org/1999/xlink">';
  print "\n";
  print "    <Service>\n";
  print "        <metadataLink xlink:type=\"simple\" xlink:href=\"$organization_url\" metadataType=\"other\"/>\n";
  print "        <description>EPA RSIG Web Server 1.0.0</description>\n";
  print "        <name>EPA_RSIG_OGC_WCS_1.0.0</name>\n";
  print "        <label>EPA RSIG Web Server 1.0.0</label>\n";
  print "        <keywords>\n";
  print "            <keyword>EPA</keyword>\n";
  print "            <keyword>RSIG</keyword>\n";
  print "            <keyword>interoperability</keyword>\n";
  print "        </keywords>\n";
  print "        <responsibleParty>\n";
  print "            <individualName>Todd Plessel</individualName>\n";
  print "            <organisationName>EPA Vislab</organisationName>\n";
  print "            <contactInfo>\n";
  print "                <onlineResource xlink:type=\"simple\" xlink:href=\"$email_address\" />\n";
  print "            </contactInfo>\n";
  print "        </responsibleParty>\n";
  print "        <fees>NONE</fees>\n";
  print "        <accessConstraints>NONE</accessConstraints>\n";
  print "    </Service>\n";
  print "    <Capability>\n";
  print "        <Request>\n";
  print "            <GetCapabilities>\n";
  print "                <DCPType>\n";
  print "                    <HTTP>\n";
  print "                        <Get>\n";
  print "                            <OnlineResource xlink:type=\"simple\" xlink:href=\"$webservice_url\"/>\n";
  print "                        </Get>\n";
  print "                        <Post>\n";
  print "                            <OnlineResource xlink:type=\"simple\" xlink:href=\"$webservice_url\"/>\n";
  print "                        </Post>\n";
  print "                    </HTTP>\n";
  print "                </DCPType>\n";
  print "            </GetCapabilities>\n";
  print "            <DescribeCoverage>\n";
  print "                <DCPType>\n";
  print "                    <HTTP>\n";
  print "                        <Get>\n";
  print "                            <OnlineResource xlink:type=\"simple\" xlink:href=\"$webservice_url\"/>\n";
  print "                        </Get>\n";
  print "                        <Post>\n";
  print "                            <OnlineResource xlink:type=\"simple\" xlink:href=\"$webservice_url\"/>\n";
  print "                        </Post>\n";
  print "                    </HTTP>\n";
  print "                </DCPType>\n";
  print "            </DescribeCoverage>\n";
  print "            <GetCoverage>\n";
  print "                <DCPType>\n";
  print "                    <HTTP>\n";
  print "                        <Get>\n";
  print "                            <OnlineResource xlink:type=\"simple\" xlink:href=\"$webservice_url\"/>\n";
  print "                        </Get>\n";
  print "                        <Post>\n";
  print "                            <OnlineResource xlink:type=\"simple\" xlink:href=\"$webservice_url\"/>\n";
  print "                        </Post>\n";
  print "                    </HTTP>\n";
  print "                </DCPType>\n";
  print "            </GetCoverage>\n";
  print "        </Request>\n";
  print "        <Exception>\n";
  print "            <Format>text/xml; charset=\"utf-8\"</Format>\n";
  print "        </Exception>\n";
  print "    </Capability>\n";
  print "    <ContentMetadata version=\"1.0.0\">\n";

  if ( defined( $callback ) ) {
    debug( "print_capabilities callback = $callback" );
    $callback->();
  } elsif ( $source ne '' && defined( $coverage_metadata_reference ) ) {
    my $variable_count = scalar @$coverage_metadata_reference;
    debug( "print_capabilities coverage_metadata_reference variable_count = $variable_count" );

    for ( my $variable_index = 0; $variable_index < $variable_count;
          ++$variable_index ) {
      my $this_entry = $coverage_metadata_reference->[ $variable_index ];
      my $this_variable = $this_entry->{ name };
      my $this_units = $this_entry->{ units };
      my $this_description = $this_entry->{ description };
      my $this_west_bound =
        exists( $this_entry->{west_bound} ) ? $this_entry->{west_bound} :
        $metadata_west_bound;
      my $this_east_bound =
        exists( $this_entry->{east_bound} ) ? $this_entry->{east_bound} :
        $metadata_east_bound;
      my $this_south_bound =
        exists( $this_entry->{south_bound} ) ? $this_entry->{south_bound} :
        $metadata_south_bound;
      my $this_north_bound =
        exists( $this_entry->{north_bound} ) ? $this_entry->{north_bound} :
        $metadata_north_bound;

      my @parts = split( /\./, $this_variable );
      my $parts_count = @parts;
      my $this_label = $parts[ $parts_count - 1 ];

      print "        <CoverageOfferingBrief>\n";
      print "            <name>$source.$this_variable</name>\n";
      print "            <label>$this_label($this_units)</label>\n";
      print "            <description>$this_description</description>\n";
      print '            <lonLatEnvelope srsName="WGS84(DD)">' . "\n";
      print "                <gml:pos>$this_west_bound $this_south_bound</gml:pos>\n";
      print "                <gml:pos>$this_east_bound $this_north_bound</gml:pos>\n";
      print "            </lonLatEnvelope>\n";
      print "        </CoverageOfferingBrief>\n";
    }
  }

  print "    </ContentMetadata>\n";
  print "</WCS_Capabilities>\n";
}



# Print web server capabilities metadata header.
# Caller may optionally pass a callback function that prints coverage offering
# for each variable as shown below.
# print_coverage_description( \&print_coverage_offering );
# Template Design Pattern.

sub print_coverage_description {
  my $callback = shift;

  print '<?xml version="1.0" encoding="UTF-8" ?>' . "\n";
  print '<CoverageDescription version="1.0.0" ';
  print 'xmlns="https://www.opengeospatial.org/standards/wcs" ';
  print 'xmlns:gml="https://www.opengis.net/gml" ';
  print 'xmlns:xlink="https://www.w3.org/1999/xlink">';
  print "\n";

  if ( defined( $callback ) ) {
    debug( "print_coverage_description callback" );
    $callback->();
  } elsif ( $source ne '' && defined( $coverage_metadata_reference ) ) {
    my $coverage_lc = lc( $coverage );
    my $variable_count = scalar @$coverage_metadata_reference;
    my $formats = get_formats();
    debug( "print_coverage_description coverage = '$coverage' " .
           "coverage_metadata_reference variable_count = $variable_count" );

    for ( my $variable_index = 0; $variable_index < $variable_count;
          ++$variable_index ) {
      my $this_entry = $coverage_metadata_reference->[ $variable_index ];
      my $this_variable = $this_entry->{ name };

      if ( $coverage eq '' || $coverage_lc eq lc( $this_variable ) ) {
        my $this_units = $this_entry->{ units };
        my $this_description = $this_entry->{ description };
        my $this_west_bound =
          exists( $this_entry->{west_bound} ) ? $this_entry->{west_bound} :
          $metadata_west_bound;
        my $this_east_bound =
          exists( $this_entry->{east_bound} ) ? $this_entry->{east_bound} :
          $metadata_east_bound;
        my $this_south_bound =
          exists( $this_entry->{south_bound} ) ? $this_entry->{south_bound} :
          $metadata_south_bound;
        my $this_north_bound =
          exists( $this_entry->{north_bound} ) ? $this_entry->{north_bound} :
          $metadata_north_bound;
        my $this_begin_date =
          exists( $this_entry->{ begin_date } ) ? $this_entry->{ begin_date } :
          $metadata_begin_date;
        my $this_end_date =
          exists( $this_entry->{ end_date } ) ? $this_entry->{ end_date } :
          $metadata_end_date;

        my @parts = split( /\./, $this_variable );
        my $parts_count = @parts;
        my $this_label = $parts[ $parts_count - 1 ];

        print "    <CoverageOffering>\n";
        print "        <name>$source.$this_variable</name>\n";
        print "        <label>$this_label($this_units)</label>\n";
        print "        <description>$this_description</description>\n";
        print "        <domainSet>\n";
        print "            <spatialDomain>\n";
        print '                <gml:Envelope srsName="WGS84(DD)">' . "\n";
        print "                    <gml:pos>$this_west_bound $this_south_bound</gml:pos>\n";
        print "                    <gml:pos>$this_east_bound $this_north_bound</gml:pos>\n";
        print "                </gml:Envelope>\n";
        print "            </spatialDomain>\n";
        print "            <temporalDomain>\n";
        print "                <timePeriod>\n";
        print "                    <beginPosition>$this_begin_date";
        print "T00:00:00Z</beginPosition>\n";

        if ( $this_end_date ne '' ) {
          print "                    <endPosition>$this_end_date";
          print "T23:59:59Z</endPosition>\n";
        }

        print "                    <timeResolution>PT1H</timeResolution>\n";
        print "                </timePeriod>\n";
        print "            </temporalDomain>\n";
        print "        </domainSet>\n";
        print "        <rangeSet>\n";
        print "            <RangeSet>\n";
        print "                <name>$source.$this_variable</name>\n";
        print "                <label>$this_label($this_units)</label>\n";
        print "                <description>$this_description</description>\n";
        print "                <nullValues>\n";
        print "                    <singleValue>-9999.0</singleValue>\n";
        print "                </nullValues>\n";
        print "            </RangeSet>\n";
        print "        </rangeSet>\n";
        print "        <supportedCRSs>\n";
        print "            <requestResponseCRSs>CRS:84</requestResponseCRSs>\n";
        print "            <nativeCRSs>CRS:84</nativeCRSs>\n";
        print "        </supportedCRSs>\n";
        print "        <supportedFormats>\n";
        print "            <formats>$formats</formats>\n";
        print "        </supportedFormats>\n";
        print "        <supportedInterpolations>\n";
        print "            <interpolationMethod>none</interpolationMethod>\n";
        print "        </supportedInterpolations>\n";
        print "    </CoverageOffering>\n";
      }
    }
  }

  print "</CoverageDescription>\n";
}



# Parse and ignore RSIG_ID option, if present.

sub parse_rsigid_option { return 1; }



# Parse service option.

sub parse_service_option {
  my $value = shift;
  my $result = parse_option( $service, $value, 'SERVICE', 'wcs' );
  return $result;
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
  my $result =
    parse_option( $request, $value, 'REQUEST',
                  'getcoverage getcapabilities describecoverage getmetadata' );
  return $result;
}



# Parse coverage option.
# Default implementation handles case of multiple variables such as
# COVERAGE=uwind,vwind,wwind
# and sets $coverage to the lookup value(s) (which may not be all lowercase
# as parsed $value is). E.g., $coverage might be 'UWIND,VWIND,WWIND' if that is
# how it is specified in the caller-supplied coverage array.
# $coverage_index is set to the index of the (last) looked-up value, e.g.,
# index of WWIND.

sub parse_coverage_option {
  my $value = shift;
  my $result = 0;
  debug( "parse_coverage_option value = '$value'\n" );

  if ( $coverage ne '' ) {
    print STDERR "\nRedundant COVERAGE option\n";
  } elsif ( defined( $custom_parse_coverage_option ) ) {
    debug( "custom_parse_coverage_option\n" );
    $coverage = $custom_parse_coverage_option->( $value );
  } else {
    my $variable_count = 0;

    if ( defined( $coverages_reference ) ) {
      $variable_count = scalar @$coverages_reference;
      debug( "coverages_reference variable_count = $variable_count\n" );
    } elsif ( defined( $coverage_metadata_reference ) ) {
      $variable_count = scalar @$coverage_metadata_reference;
      debug( "coverage_metadata_reference variable_count = $variable_count\n" );
    }

    if ( $variable_count > 0 ) {
      my @values = split( /,/, $value );
      my $value_count = @values;
      my $found_count = 0;

      for ( my $value_index = 0; $value_index < $value_count; ++$value_index ) {
        my $this_value = $values[ $value_index ];
        my $found_this_value = 0;

        for ( my $variable_index = 0; $variable_index < $variable_count;
              ++$variable_index ) {
          my $this_coverage = '';

          if ( defined( $coverages_reference ) ) {
            $this_coverage = $coverages_reference->[ $variable_index ];
          } else {
            my $this_entry = $coverage_metadata_reference->[ $variable_index ];
            $this_coverage = $this_entry->{ name };
          }

          if ( lc( $this_coverage ) eq $this_value ) {

            if ( $coverage eq '' ) {
              $coverage = $this_coverage;
            } else {
              $coverage .= ',' . $this_coverage;
            }

            $coverage_index = $variable_index;
            $variable_index = $variable_count;
            $found_this_value = 1;
            ++$found_count;
          }
        }

        if ( ! $found_this_value ) {
          $value_index = $value_count;
        }
      }

      if ( $found_count != $value_count ) {
        $coverage = '';
      }
    }
  }

  if ( $coverage eq '' ) {
    print STDERR "\nInvalid COVERAGE option\n";
  } else {
    $result = 1;
  }

  return $result;
}



# Parse time option.
# Initializes variables:
# $time, $yyyymmdd1, $yyyymmdd2, $yyyymmddhhmmss1, $yyyymmddhhmmss2,
# $days, $hours.

sub parse_time_option {
  my $value = shift;
  my $result = 0;

  if ( $time ne '' ) {
    print STDERR "\nRedundant TIME option\n";
  } else {
    my $is_valid = is_valid_time( $value );

    if ( ! $is_valid ) {
      print STDERR "\nInvalid TIME option\n";
    } else {
      $time = $value;
      $yyyymmdd1 =
        int( substr( $time, 0, 4 ) .
             substr( $time, 5, 2 ) .
             substr( $time, 8, 2 ) );
      $yyyymmdd2 = $yyyymmdd1;

      $yyyymmddhhmmss1 = $yyyymmdd1 .
        substr( $time, 11, 2 ) .
        substr( $time, 14, 2 ) .
        substr( $time, 17, 2 );
      $yyyymmddhhmmss2 = $yyyymmddhhmmss1;

      if ( length( $time ) == 41 ) {
        $yyyymmdd2 = int( substr( $time, 21, 4 ) .
                          substr( $time, 26, 2 ) .
                          substr( $time, 29, 2 ) );
        $yyyymmddhhmmss2 = $yyyymmdd2 .
          substr( $time, 32, 2 ) .
          substr( $time, 35, 2 ) .
          substr( $time, 38, 2 );
      }

      my $yyyymmdd = $yyyymmdd1;
      $days = 1;

      while ( $yyyymmdd < $yyyymmdd2 ) {
        $yyyymmdd = increment_date( $yyyymmdd );
        ++$days;
      }

      my $hh = int( substr( $yyyymmddhhmmss1, 8, 2 ) );
      my $yyyymmddhh = int( $yyyymmdd1 * 100 + $hh );
      my $hh2 = int( substr( $yyyymmddhhmmss2, 8, 2 ) );
      my $yyyymmddhh2 = int( $yyyymmdd2 * 100 + $hh2 );
      $hours = 1;

      while ( $yyyymmddhh < $yyyymmddhh2 ) {
        $yyyymmddhh = increment_hour( $yyyymmddhh );
        ++$hours;
      }

      debug( "$yyyymmddhhmmss1, $yyyymmddhhmmss2, $yyyymmdd1, $yyyymmdd2, " .
             "$days, $hours" );
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
    print STDERR "\nRedundant BBOX option\n";
  } else {
    $bbox = $value;
    @bounds = split( /,/, $bbox );
    my $bounds_count = @bounds;

    if ( $bounds_count != 4 && $bounds_count != 6 ) {
      print STDERR "\nInvalid BBOX option (requires 4 or 6 numbers).\n";
    } elsif ( ! in_range( $bounds[ 0 ], -180.0, 180.0 ) ) {
      print STDERR "\nInvalid BBOX option (1st number)\n";
    } elsif ( ! in_range( $bounds[ 1 ], -90.0, 90.0 ) ) {
      print STDERR "\nInvalid BBOX option (2nd number)\n";
    } elsif ( ! in_range( $bounds[ 2 ], $bounds[ 0 ], 180.0 ) ) {
      print STDERR "\nInvalid BBOX option (3rd number)\n";
    } elsif ( ! in_range( $bounds[ 3 ], $bounds[ 1 ], 90.0 ) ) {
      print STDERR "\nInvalid BBOX option (4th number)\n";
    } else {
      $result = 1;
    }
  }

  return $result;
}



# Parse format option:

sub parse_format_option {
  my $value = shift;
  my $result = 0;

  if ( $format ne '' ) {
    print STDERR "\nRedundant FORMAT option\n";
  } else {

    if ( ! exists( $content_types{ $value } ) ) {
      print STDERR "\nInvalid FORMAT option\n";
    } else {
      $format = $value;
      $result = 1;
    }
  }

  return $result;
}



# Parse compress option.

sub parse_compress_option {
  my $value = shift;
  my $result = 0;

  if ( $compress ne '' ) {
    print STDERR "\nRedundant COMPRESS option\n";
  } else {
    my $is_valid = $value eq '0' || $value eq '1';

    if ( ! $is_valid ) {
      print STDERR "\nInvalid COMPRESS option (requires 0 or 1)\n";
    } else {
      $compress = int( $value );
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
    print STDERR "\nRedundant LAMBERT option\n";
  } else {
    $lambert = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 4 ) {
      print STDERR "\nInvalid LAMBERT option (requires 4 numbers)\n";
    } elsif ( ! in_range( $values[ 0 ], -89.0, 89.0 ) ) {
      print STDERR "\nInvalid LAMBERT option (1st number)\n";
    } elsif ( ! in_range( $values[ 1 ], -89.0, 89.0 ) ) {
      print STDERR "\nInvalid LAMBERT option (2nd number)\n";
    } elsif ( ! in_range( $values[ 2 ], -180.0, 180.0 ) ) {
      print STDERR "\nInvalid LAMBERT option: (3rd number)'\n";
    } elsif ( ! in_range( $values[ 3 ], -89.0, 89.0 ) ) {
      print STDERR "\nInvalid LAMBERT option:  (4th number)\n";
    } elsif ( $values[ 0 ] > $values[ 1 ] ) {
      print STDERR "\nInvalid LAMBERT option: (require 1st > 2nd number)\n";
    } elsif ( $values[ 0 ] > 0.0 && $values[ 1 ] < 0.0 ) {
      print STDERR "\nInvalid LAMBERT option: (require 1st and 2nd same sign)\n";
    } elsif ( $values[ 0 ] < 0.0 && $values[ 1 ] > 0.0 ) {
      print STDERR "\nInvalid LAMBERT option: (require 1st and 2nd same sign)\n";
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
    print STDERR "\nRedundant STEREOGRAPHIC option\n";
  } else {
    $stereographic = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 3 ) {
      print STDERR "\nInvalid STEREOGRAPHIC option (require 3 numbers)\n";
    } elsif ( ! in_range( $values[ 0 ], -180.0, 180.0 ) ) {
      print STDERR "\nInvalid STEREOGRAPHIC option: (require 1st in [-180, 180])\n";
    } elsif ( ! in_range( $values[ 1 ], -90.0, 90.0 ) ) {
      print STDERR "\nInvalid STEREOGRAPHIC option: (require 2nd in [-90, 90])\n";
    } elsif ( ! in_range( $values[ 2 ], -90.0, 90.0 ) ) {
      print STDERR "\nInvalid STEREOGRAPHIC option: (require 3rd in [-90, 90])\n";
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
    print STDERR "\nRedundant MERCATOR option\n";
  } else {
    $mercator = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 1 ) {
      print STDERR "\nInvalid MERCATOR option (require 1 number)\n";
    } elsif ( ! in_range( $values[ 0 ], -180.0, 180.0 ) ) {
      print STDERR "\nInvalid MERCATOR option (require number in [-180, 180])\n";
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
    print STDERR "\nRedundant LONLAT option\n";
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
    print STDERR "\nRedundant ELLIPSOID option\n";
  } else {
    $ellipsoid = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 1 && $count != 2 ) {
      print STDERR "\nInvalid ELLIPSOID option (require 1 or 2 numbers)\n";
    } elsif ( ! in_range( $values[ 0 ], 1.0, 1e10 ) ) {
      print STDERR "\nInvalid ELLIPSOID option (1st number)\n";
    } elsif ( $count == 2 ) {

      if ( ! in_range( $values[ 1 ], $values[ 0 ], 1e10 ) ) {
        print STDERR "\nInvalid ELLIPSOID option: (2nd number)\n";
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
    print STDERR "\nRedundant GRID option\n";
  } else {
    $grid = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 6 ) {
      print STDERR "\nInvalid GRID option (require 6 numbers)\n";
    } elsif ( $values[ 0 ] < 1 ) {
      print STDERR "\nInvalid GRID option (1st number)\n";
    } elsif ( $values[ 1 ] < 1 ) {
      print STDERR "\nInvalid GRID option (2nd number)\n";
    } elsif ( ! in_range( $values[ 2 ], -1e10, 1e10 ) ) {
      print STDERR "\nInvalid GRID option (3rd number)\n";
    } elsif ( ! in_range( $values[ 3 ], -1e10, 1e10 ) ) {
      print STDERR "\nInvalid GRID option (4th number)\n";
    } elsif ( ! in_range( $values[ 4 ], 0.01, 1e10 ) ) {
      print STDERR "\nInvalid GRID option (5th number)\n";
    } elsif ( ! in_range( $values[ 5 ], 0.01, 1e10 ) ) {
      print STDERR "\nInvalid GRID option (6th number)\n";
    } else {
      $result = 1;
    }
  }

  return $result;
}



# Parse levels option.

sub parse_levels_option {
  my $value = shift;
  my $result = 0;

  if ( $levels ne '' ) {
    print STDERR "\nRedundant LEVELS option\n";
  } else {
    $levels = $value;
    my @values = split( /,/, $value );
    my $count = @values;

    if ( $count != 3 + $values[ 0 ] + 1 + 5 ) {
      print STDERR "\nInvalid LEVELS option (count)\n";
    } elsif ( ! in_range( $values[ 0 ], 1, 100 ) ) {
      print STDERR "\nInvalid LEVELS option (1st number)\n";
    } elsif ( ! in_range( $values[ 1 ], 1, 7 ) ) {
      print STDERR "\nInvalid LEVELS option (2nd number)\n";
    } elsif ( ! in_range( $values[ 2 ], 0.0, 1e8 ) ) {
      print STDERR "\nInvalid LEVELS option (3rd number)\n";
    } else {

      # Ensure level values are ordered decreasing if sigma-p else increasing:

      my $end = $count - 5;
      my $minimum = $lonlat ne '' ?   -500.0 : 0.0;
      my $maximum = $lonlat ne '' ? 100000.0 : 1.0;
      my $ok = 0;

      for ( my $index = 3; $index < $end; ++$index ) {
        my $level_value = $values[ $index ];
        $ok = in_range( $level_value, $minimum, $maximum );

        if ( ! $ok ) {
          print STDERR "\nInvalid LEVELS option (number $index)\n";
          $index = $count;
        } elsif ( $lonlat ne '' ) {
          $minimum = $level_value;
        } else {
          $maximum = $level_value;
        }
      }

      # Check last parameter constants:

      if ( $ok == 1 ) {
        my $index = $count - 5;

        if ( ! in_range( $values[ $index ], 0.01, 1e2 ) ) {
          print STDERR "\nInvalid LEVELS option\n";
        } elsif ( ! in_range( $values[ $index + 1 ], 0.01, 1e4 ) ) {
          print STDERR "\nInvalid LEVELS option\n";
        } elsif ( ! in_range( $values[ $index + 2 ], 0.01, 1e4 ) ) {
          print STDERR "\nInvalid LEVELS option\n";
        } elsif ( ! in_range( $values[ $index + 3 ], 0.01, 1e4 ) ) {
          print STDERR "\nInvalid LEVELS option\n";
        } elsif ( ! in_range( $values[ $index + 4 ], 0.01, 1e6 ) ) {
          print STDERR "\nInvalid LEVELS option\n";
        } else {
          $result = 1;
        }
      }
    }
  }

  return $result;
}



# Check that all required 'options' have been specified:

sub required_options_specified {
  my $result = 0;

  if ( $service eq '' ) {
    print STDERR "\nMissing option: 'SERVICE'\n";
  } elsif ( $version eq '' ) {
    print STDERR "\nMissing option: 'VERSION'\n";
  } elsif ( $request eq '' ) {
    print STDERR "\nMissing option: 'REQUEST'\n";
  } elsif ( $request eq 'getcoverage' || $request eq 'getmetadata' ) {

    if ( $coverage eq '' ) {
      print STDERR "\nMissing option: 'COVERAGE'\n";
    } elsif ( $format eq '' ) {
      print STDERR "\nMissing option: 'FORMAT'\n";
    } elsif ( $time eq '' ) {
      print STDERR "\nMissing option: 'TIME'\n";
    } elsif ( $bbox eq '' ) {
      print STDERR "\nMissing option: 'BBOX'\n";
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
        print STDERR "\nInvalid options: 'REGRID/";
        print STDERR "LAMBERT/STEREOGRAPHIC/MERCATOR/LONLAT/";
        print STDERR "GRID/ELLIPSOID'\n";
      }
    }
  } else {
    $result = 1;
  }

  return $result;
}



# my $result = parse_option( $option, $value, $option_name, $valid_values );
# my $result = parse_option( $variable, $value, 'COVERAGE', 'ozone pm25' );

sub parse_option {
  my ( $option, $value, $option_name, $valid_values ) = @_;
  my $result = 0;

  if ( $option ne '' ) {
    print STDERR "\nRedundant $option_name option\n";
  } else {
    my $start = index( lc( " $valid_values " ), lc( " $value " ) );

    if ( $start > -1 ) {
      my $case_matched_value =
        substr( " $valid_values ", $start + 1, length( $value ) );
      $_[ 0 ] = $case_matched_value;
      $result = 1;
    } else {
      print STDERR "\nInvalid $option_name option\n";
    }
  }

  return $result;
}



# my $result = parse_range_option( $option, $value, $option_name ) ;
# my $result = parse_range_option( $aod_range, $value, 'AOD_RANGE' );

sub parse_range_option {
  my ( $option, $value, $option_name ) = @_;
  my $result = 0;
  my @values = split( /,/, $value );
  my $count = @values;

  if ( $option ne '' ) {
    print STDERR "\nRedundant $option_name option\n";
  } elsif ( $count != 2 ) {
    print STDERR "\nRedundant $option_name option (requires 2 numbers)\n";
  } else {
    $_[ 0 ] = $value;
    $result = 1;
  }

  return $result;
}



# If $regrid and $formats is not xdr then this function returns a string
# containing the various command-line arguments to XDRConvert pertaining to
# regridding (including regrid aggregating) based on the parsed WCS options.
# my $regrid_args = xdrconvert_regrid_args( 1 );
# 1 = include LEVELS, 0 = do not include LEVELS.

sub xdrconvert_regrid_args {
  my $with_levels = shift;
  my $result = '';

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

    my $layer_args = '';

    if ( $with_levels ) {

      if ( $levels ne '' ) {
        $layer_args = "-layers $levels ";
        $layer_args =~ tr/,/ /;
      }
    }

    $result =
      "-regrid $regrid " .
      $projection_args .
      "-ellipsoid $major_semiaxis $minor_semiaxis " .
      $grid_args .
      $layer_args .
      $regrid_aggregate_option;
  }

  return $result;
}



# my $my_xdrconvert_format = xdrconvert_format();
# E.g., If FORMAT=netcdf-coads the result is 'coards'.

sub xdrconvert_format {
  my $result = '';

  if ( $regrid ne '' || $format ne 'xdr' ) {
    $result = $format;
    my @format_parts = split( /-/, $format );
    my $format_parts_count = @format_parts;

    if ( $format_parts_count == 2 ) {
      $result = $format_parts[ 1 ];
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



# my $leap = is_leap_year( 2016 );
# leap is 1

sub is_leap_year {
  my $year = shift;
  my $result = $year % 4 == 0 && ( $year % 100 != 0 || $year % 400 == 0 );
  return $result;
}



# my $days = days_in_month( 2016, 2 );
# days is 29.

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



# my $is_valid = is_valid_time( '2001-08-26t20:00:00z/2001-08-27t23:59:59z' );

sub is_valid_time {
  my $time_string = shift;
  my $result = 0;
  my $length = length( $time_string );

  if ( $length == 41 ) {
    $result = is_valid_time( substr( $time_string, 0, 20 ) );
    $result = $result && substr( $time_string, 20, 1 ) eq '/';
    $result = $result && is_valid_time( substr( $time_string, 21, 20 ) );
    $result = $result && substr( $time_string, 0, 20 ) le substr( $time_string, 21, 20 );
  } elsif ( $length == 20 ) {
    my $year   = substr( $time_string, 0, 4 );
    my $month  = substr( $time_string, 5, 2 );
    my $day    = substr( $time_string, 8, 2 );
    my $hour   = substr( $time_string, 11, 2 );
    my $minute = substr( $time_string, 14, 2 );
    my $second = substr( $time_string, 17, 2 );
    $result = in_range( $year, 1900, 3000 );
    $result = $result && in_range( $month, 1, 12 );
    $result = $result && in_range( $day, 1, days_in_month( $year, $month ) );
    $result = $result && in_range( $hour, 0, 23 );
    $result = $result && in_range( $minute, 0, 59 );
    $result = $result && in_range( $second, 0, 59 );
    $result = $result && substr( $time_string, 4, 1 ) eq '-';
    $result = $result && substr( $time_string, 7, 1 ) eq '-';
    $result = $result && substr( $time_string, 10, 1 ) eq 't';
    $result = $result && substr( $time_string, 13, 1 ) eq ':';
    $result = $result && substr( $time_string, 16, 1 ) eq ':';
    $result = $result && substr( $time_string, 19, 1 ) eq 'z';
  }

  return $result;
}



# my $yyyymmdd = increment_date( 20160229 );
# $yyyymmdd is 20160301

sub increment_date {
  my $yyyymmdd  = shift;
  my $yyyy = int( $yyyymmdd / 10000 );
  my $mm   = int( $yyyymmdd / 100 ) % 100;
  my $dd   = int( $yyyymmdd ) % 100;

  ++$dd;

  if ( $dd > 28 ) {
    my $days_this_month = days_in_month( $yyyy, $mm );

    if ( $dd > $days_this_month ) {
      $dd = 1;
      ++$mm;

      if ( $mm > 12 ) {
        $mm = 1;
        ++$yyyy;
      }
    }
  }

  $yyyymmdd = int( $yyyy * 10000 + $mm * 100 + $dd );
  return $yyyymmdd;
}



# my $yyyymmdd = decrement_date( 20160301 );
# $yyyymmdd is 20160229

sub decrement_date {
  my $yyyymmdd  = shift;
  my $yyyy = int( $yyyymmdd / 10000 );
  my $mm   = int( $yyyymmdd / 100 ) % 100;
  my $dd   = int( $yyyymmdd ) % 100;

  --$dd;

  if ( $dd < 1 ) {
    --$mm;

    if ( $mm < 1 ) {
      $mm = 12;
      --$yyyy;
    }

    $dd = days_in_month( $yyyy, $mm );
  }

  $yyyymmdd = int( $yyyy * 10000 + $mm * 100 + $dd );
  return $yyyymmdd;
}



# my $yyyymmdd = increment_days( $yyyymmdd, $days );

sub increment_days {
  my ( $yyyymmdd, $days ) = @_;

  while ( $days-- ) {
    $yyyymmdd = increment_date( $yyyymmdd );
  }

  return $yyyymmdd;
}



# increment_month( $yyyy, $mm );

sub increment_month {
  my ( $yyyy, $mm ) = @_;
  my $months_per_year = 12;

  ++$mm;

  if ( $mm > $months_per_year ) {
    $mm = 1;
    ++$yyyy;
  }

  $_[ 0 ] = $yyyy;
  $_[ 1 ] = $mm;
}



# my $yyyymmddhh = increment_hour( 2016022923 );
# $yyyymmddhh is 2016030100

sub increment_hour {
  my $yyyymmddhh  = shift;
  my $yyyy = int( $yyyymmddhh / 1000000 );
  my $mm   = int( $yyyymmddhh / 10000 ) % 100;
  my $dd   = int( $yyyymmddhh / 100 ) % 100;
  my $hh   = int( $yyyymmddhh ) % 100;

  ++$hh;

  if ( $hh > 23 ) {
    $hh = 0;
    ++$dd;

    if ( $dd > 28 ) {
      my $days_this_month = days_in_month( $yyyy, $mm );

      if ( $dd > $days_this_month ) {
        $dd = 1;
        ++$mm;

        if ( $mm > 12 ) {
          $mm = 1;
          ++$yyyy;
        }
      }
    }
  }

  $yyyymmddhh = int( $yyyy * 1000000 + $mm * 10000 + $dd * 100 + $hh );
  return $yyyymmddhh;
}



# my $yyyymmddhh = decrement_hour( 2016030100 );
# $yyyymmddhh is 2016022923

sub decrement_hour {
  my $yyyymmddhh  = shift;
  my $yyyy = int( $yyyymmddhh / 1000000 );
  my $mm   = int( $yyyymmddhh / 10000 ) % 100;
  my $dd   = int( $yyyymmddhh / 100 ) % 100;
  my $hh   = int( $yyyymmddhh ) % 100;

  --$hh;

  if ( $hh < 0 ) {
    $hh = 23;
    --$dd;

    if ( $dd < 1 ) {
      --$mm;

      if ( $mm < 1 ) {
        $mm = 12;
        --$yyyy;
      }

      $dd = days_in_month( $yyyy, $mm );
    }
  }

  $yyyymmddhh = int( $yyyy * 1000000 + $mm * 10000 + $dd * 100 + $hh );
  return $yyyymmddhh;
}



# my $yyyymmddhhmmss = increment_second( '20160229235959' );
# $yyyymmddhhmmss is 20160301000000

sub increment_second {
  my $yyyymmddhhmmss = shift;
  my $months_per_year    = 12;
  my $hours_per_day      = 24;
  my $minutes_per_hour   = 60;
  my $seconds_per_minute = 60;
  my $yyyymmdd = int( substr( $yyyymmddhhmmss, 0, 8 ) );
  my $yyyy = int( $yyyymmdd / 10000 );
  my $mm   = int( $yyyymmdd / 100 ) % 100;
  my $dd   = int( $yyyymmdd % 100 );
  my $hhmmss = substr( $yyyymmddhhmmss, 8, 6 );
  my $hh   = int( $hhmmss / 10000 );
  my $min  = int( $hhmmss / 100 ) % 100;
  my $ss   = int( $hhmmss % 100 );

  ++$ss;

  if ( $ss > $seconds_per_minute ) {
    $ss = 0;
    ++$min;

    if ( $min > $minutes_per_hour ) {
      $min = 0;
      ++$hh;

      if ( $hh > $hours_per_day ) {
        $hh = 0;
        ++$dd;

        if ( $dd > 28 ) {
          my $days_this_month = days_in_month( $yyyy, $mm );

          if ( $dd > $days_this_month ) {
            $dd = 1;
            ++$mm;

            if ( $mm > $months_per_year ) {
              $mm = 1;
              ++$yyyy;
            }
          }
        }

        $yyyymmdd = int( $yyyy * 10000 + $mm * 100 + $dd );
      }
    }
  }

  $hhmmss = sprintf( "%02d%02d%02d", $hh, $min, $ss );
  $yyyymmddhhmmss = $yyyymmdd . $hhmmss;

  return $yyyymmddhhmmss;
}



# Split yyyymmdd into yyyy mm0 dd0 (0-prefixed):
# my ($yyyy, $mm0, $dd0) = split_yyyymmdd(20160229);
# yyyy is 2016
# mm0 is "02"
# dd0 is "29"

sub split_yyyymmdd {
  my $yyyymmdd = shift;
  my $yyyy = substr( $yyyymmdd, 0, 4 );
  my $mm0  = substr( $yyyymmdd, 4, 2 );
  my $dd0  = substr( $yyyymmdd, 6, 2 );
  return $yyyy, $mm0, $dd0;
}



# Split yyyymmddhh into yyyy mm0 dd0 hh0 (0-prefixed):
# my ($yyyy, $mm0, $dd0, $hh0) = split_yyyymmdd(2016022900);
# yyyy is 2016
# mm0 is "02"
# dd0 is "29"
# hh0 is "00"

sub split_yyyymmddhh {
  my $yyyymmddhh = shift;
  my ($yyyy, $mm0, $dd0) = split_yyyymmdd( int( $yyyymmddhh / 100 ) );
  my $hh0 = substr( $yyyymmddhh, 8, 2 );
  return $yyyy, $mm0, $dd0, $hh0;
}



# my $yyyymmdd = convert_yyyyddd( 2009186 );
# $yyyymmdd is 20090826

sub convert_yyyyddd {
  my $yyyyddd = shift;
  my $yyyy = int( $yyyyddd / 1000 );
  my $ddd = int( $yyyyddd % 1000 );
  my $month = 0;
  my $day = 0;

  for ( $month = 1; $month <= 12 && $ddd > 0; ++$month ) {
    $day = $ddd;
    $ddd -= days_in_month( $yyyy, $month );
  }

  --$month;
  my $yyyymmdd = $yyyy * 10000 + $month * 100 + $day;
  return $yyyymmdd;
}



# my $yyyyddd = to_yyyyddd( 20090826 );
# $yyyyddd is 2009186

sub to_yyyyddd {
  my $yyyymmdd = shift;
  my $yyyy = int( $yyyymmdd / 10000 );
  my $mm   = int( $yyyymmdd / 100 ) % 100;
  my $dd   =  $yyyymmdd % 100;
  my $ddd = 0;

  for ( my $month = 1; $month < $mm; ++$month ) {
    $ddd += days_in_month( $yyyy, $month );
  }

  $ddd += $dd;
  my $yyyyddd = $yyyy * 1000 + $ddd;
  return $yyyyddd;
}



# convert_timestamp( '20090826T1500', $yyyy, $ddd, $hour );
# $yyyy is 2009
# $ddd is '186'
# $hour is '15.

sub convert_timestamp {
  my ( $timestamp, $yyyy, $ddd, $hour ) = @_;
  my $timestamp_length = length( $timestamp );
  my $month = 0;
  my $day   = 0;
  $yyyy = 0;
  $ddd  = 0;
  $hour = 0;

  if ( $timestamp_length == 13 ) {
    $yyyy  = substr( $timestamp, 0, 4 );
    $month = substr( $timestamp, 4, 2 );
    $day   = substr( $timestamp, 6, 2 );
    $hour  = substr( $timestamp, 9, 2 );

    my $t = sprintf( "%04d-%02d-%02dt%02d:00:00z", $yyyy, $month, $day, $hour);

    if ( is_valid_time( $t ) ) {

      for ( my $mo = 1; $mo < $month; ++$mo ) {
        $ddd += days_in_month( $yyyy, $mo );
      }

      $ddd += $day;
    }
  }

  $_[ 1 ] = $yyyy;
  $_[ 2 ] = $ddd;
  $_[ 3 ] = $hour;
}



# my $ddd = day_of_year( $year, $month, $day );

sub day_of_year {
  my ( $year, $month, $day ) = @_;
  my $result = 0;
  my $mm = 1;

  while ( $mm < $month ) {
    $result += days_in_month( $year, $mm );
    ++$mm;
  }

  $result += $day;
  return $result;
}



# increment_date_day( $yyyy, $mm, $dd );

sub increment_date_day {
  my ( $yyyy, $mm, $dd ) = @_;
  my $months_per_year = 12;
  my $days_this_month = days_in_month( $yyyy, $mm );

  ++$dd;

  if ( $dd > $days_this_month ) {
    $dd = 1;
    ++$mm;

    if ( $mm > $months_per_year ) {
      $mm = 1;
      ++$yyyy;
    }
  }

  $_[ 0 ] = $yyyy;
  $_[ 1 ] = $mm;
  $_[ 2 ] = $dd;
}



# increment_date_days( $yyyy, $mm, $dd, $days );

sub increment_date_days {
  my ( $yyyy, $mm, $dd, $days ) = @_;

  while ( $days-- ) {
    increment_date_day( $yyyy, $mm, $dd );
  }

  $_[ 0 ] = $yyyy;
  $_[ 1 ] = $mm;
  $_[ 2 ] = $dd;
}



# decrement_date_day( $yyyy, $mm, $dd );

sub decrement_date_day {
  my ( $yyyy, $mm, $dd ) = @_;
  --$dd;

  if ( $dd < 1 ) {
    --$mm;

    if ( $mm < 1 ) {
      --$yyyy;
      $mm = 12;
    }

    $dd = days_in_month( $yyyy, $mm );
  }

  $_[ 0 ] = $yyyy;
  $_[ 1 ] = $mm;
  $_[ 2 ] = $dd;
}



# seven_days_of_year( $yyyy, $mm, $dd, $ddd, $yyyy2, $ddd2 );

sub seven_days_of_year {
  my ( $yyyy, $mm, $dd, $ddd, $yyyy2, $ddd2 ) = @_;

  # Back-up to beginning of week:

  $ddd = day_of_year( $yyyy, $mm, $dd );

  while ( ( $ddd - 1 ) % 7 ) {
    decrement_date_day( $yyyy, $mm, $dd );
    $ddd = day_of_year( $yyyy, $mm, $dd );
  }

  $ddd = day_of_year( $yyyy, $mm, $dd );
  debug( "  $yyyy $mm $dd = $ddd" );

  # Compute end of week:

  $yyyy2 = $yyyy;
  my $mo   = $mm;
  my $d    = $dd;
  increment_date_days( $yyyy2, $mo, $d, 6 );
  $ddd2 = day_of_year( $yyyy2, $mo, $d );
  debug( "  $yyyy2 $mo $d = $ddd2" );

  $_[ 0 ] = $yyyy;
  $_[ 1 ] = $mm;
  $_[ 2 ] = $dd;
  $_[ 3 ] = $ddd;
  $_[ 4 ] = $yyyy2;
  $_[ 5 ] = $ddd2;
}



# Check if given timestamp is within the range of $time.
# my $ok = within_time_range( $timestamp );

sub within_time_range {
  my $timestamp = shift;
  my $result = is_valid_time( $timestamp );

  if ( $result ) {
    my $length = length( $time );

    if ( $length == 20 ) {
      $result = ( $timestamp eq substr( $time, 0, 20 ) );
    } elsif ( $length == 41 ) {
      my $first_timestamp = substr( $time, 0, 20 );
      my $last_timestamp  = substr( $time, 21, 20 );
      $result =
        ( $timestamp ge $first_timestamp && $timestamp le $last_timestamp );
    }
  }

  return $result;
}



# Hours offset from UTC:

my $timezone_offset_table = {
  AST => -4,
  EST => -5,
  EDT => -4,
  CST => -6,
  CDT => -5,
  MST => -7,
  MDT => -6,
  PST => -8,
  PDT => -7,
  AKST => -9,
  AKDT => -8,
  HAST => -10,
  HASDT => -9
};

# Compute UTC timestamp from local one.
# my $timestamp = UTC_timestamp( $yyyy_mm_dd, $hh_mm_ss, $timezone );

sub UTC_timestamp {
  my ( $yyyy_mm_dd, $hh_mm_ss, $timezone ) = @_;
  my $yyyy = substr( $yyyy_mm_dd, 0, 4 );
  my $mm   = substr( $yyyy_mm_dd, 5, 2 );
  my $dd   = substr( $yyyy_mm_dd, 8, 2 );
  my $hh_mm_ss_length = length( $hh_mm_ss );
  my $hh   = substr( $hh_mm_ss, 0, 2 );
  my $mn   = $hh_mm_ss_length >= 5 ? substr( $hh_mm_ss, 3, 2 ) : '00';
  my $ss   = $hh_mm_ss_length >= 8 ? substr( $hh_mm_ss, 6, 2 ) : '00';
  my $mm_ss = "$mn:$ss";

  if ( exists( $timezone_offset_table->{ $timezone } ) ) {
    my $hours_offset_from_utc = $timezone_offset_table->{ $timezone };
    $hh += -$hours_offset_from_utc;
  }

  increment_timestamp( $yyyy, $mm, $dd, $hh );
  my $result =
    sprintf( "%04d-%02d-%02dt%02d:%sz",
             $yyyy, $mm, $dd, $hh, $mm_ss );
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



# increment_timestamp_hours( $yyyy, $mm, $dd, $hh, $hours );

sub increment_timestamp_hours {
  my ( $yyyy, $mm, $dd, $hh, $hours ) = @_;

  if ( $hours > 0 ) {

    while ( $hours-- ) {
      increment_timestamp( $yyyy, $mm, $dd, $hh );
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



# my $string = truncate_string( 44, $string );

sub truncate_string {
  my $len = shift;
  my $string = shift;
  my $string_len = length( $string );
  my $i = 0;

  # Skip leading spaces:

  while ( $i < $string_len && substr( $string, $i, 1 ) eq ' ' ) {
    ++$i;
  }

  my $result = substr( $string, $i, $len );
  return $result;
}



# Change commas to spaces in double-quoted portions of a string and
# remove double quotes:
# my $result = erase_double_quoted_commas_and_quotes( $string );

sub erase_double_quoted_commas_and_quotes {
  my $string = shift;
  my $count = length( $string );
  my $quoted = 0;
  my $result = '';

  for ( my $i = 0; $i < $count; ++$i ) {
    my $c = substr( $string, $i, 1 );

    if ( $c eq '"' ) {
      $quoted = ! $quoted;
    } elsif ( $c eq ',' ) {

      if ( $quoted ) {
        $result .= ' ';
      } else {
        $result .= ',';
      }
    } elsif ( $c ne "\r" ) {
      $result .= $c;
    }
  }

  return $result;
}



# my $trig_angle = bearing_to_trig_angle( $bearing );

sub bearing_to_trig_angle {
  my $bearing = shift;

  # Convert bearing angle to standard trig angle:

  my $angle = 90.0 - $bearing;

  if ( $angle < 0.0 ) {
    $angle += 360.0;
  }

  # Convert 'from_direction' to 'to_direction':

  $angle -= 180.0;

  if ( $angle < 0.0 ) {
    $angle += 360.0;
  }

  return $angle;
}



# my ( $u, $v ) = speed_bearing_to_u_v( $speed, $bearing );

sub speed_bearing_to_u_v {
  my $speed = shift;
  my $bearing = shift;
  my $angle = bearing_to_trig_angle( $bearing );
  my ( $u, $v ) = speed_angle_to_u_v( $speed, $angle );
  return $u, $v;
}



# my ( $u, $v ) = speed_angle_to_u_v( $speed, $angle );

sub speed_angle_to_u_v {
  my $speed = shift;
  my $angle = shift;
  my $to_radians = 0.017453292519943; # Pi / 180.
  my $theta = $angle * $to_radians;
  my $u = $speed * cos( $theta );
  my $v = $speed * sin( $theta );
  return $u, $v;
}



# Match or convert units and measure.
# my $converted_measure = matches_units( $from_units, $to_units, $measure );

sub matches_units {
  my ( $from_units, $to_units, $measure ) = @_;
  my $result = -9999; # Missing/invalid value.

  # Various conversion scale factors:

  my $feet_to_meters = 0.3048;
  my $to_micro = 1e6;
  my $milli_to_micro = 1e3;
  my $micro_to_milli = 1e-3;

  if ( $from_units eq $to_units ) {
    $result = $measure;
  } elsif ( lc( $from_units ) eq lc( $to_units ) ) {
    $result = $measure;
  } elsif ( $to_units eq '%' && index( $from_units, '%' ) == 0 ) {
    $result = $measure;
  } elsif ( ( $from_units eq 'degF' || $from_units eq 'F' )
            && $to_units eq 'C' ) {
    $result = ( $measure - 32.0 ) * ( 5.0 / 9.0 );
  } elsif ( $from_units eq 'degC' && $to_units eq 'C' ) {
    $result = $measure;
  } elsif ( $from_units eq 'perm' && $to_units eq '/m' ) {
    $result = $measure;
  } elsif ( $from_units eq 'ft' && $to_units eq 'm' ) {
    $result = $measure * $feet_to_meters;
  } elsif ( ( $from_units eq 'cfs' || $from_units eq 'ft3/s' )
            && $to_units eq 'm3/s' ) {
    $result = $measure * $feet_to_meters * $feet_to_meters * $feet_to_meters;
  } elsif ( $from_units eq 'mg/m3' && $to_units eq 'mg/L' ) {
    $result = $measure * 1e-3;
  } elsif ( $from_units eq 'mg/lasN' && $to_units eq 'mg/L' ) {
    $result = $measure;
  } elsif ( $from_units eq 'ug/l' && $to_units eq 'mg/L' ) {
    $result = $measure * $micro_to_milli;
  } elsif ( $from_units eq 'mg/l' && $to_units eq 'ug/L' ) {
    $result = $measure * $milli_to_micro;
  } elsif ( $from_units eq 'mho/cm' && $to_units eq 'uS/cm' ) {
    $result = $measure * $to_micro;
  } elsif ( $from_units eq 'umho/cm' && $to_units eq 'uS/cm' ) {
    $result = $measure;
  } elsif ( $from_units eq 'mS/cm' && $to_units eq 'uS/cm' ) {
    $result = $measure * $milli_to_micro;
  } elsif ( $from_units eq 'PSS' && $to_units eq 'PSU' ) {
    $result = $measure;
  } elsif ( $from_units eq 'ppt' && $to_units eq 'ppm' ) {
    $result = $measure * 1e3;
  } elsif ( $from_units eq 'ppth' && $to_units eq 'ppm' ) {
    $result = $measure * 1e3;
  } elsif ( $from_units eq 'uE/m2/sec' && $to_units eq 'uE/m2/s' ) {
    $result = $measure;
  } elsif ( $to_units eq 'pH' ) {
    $result = $measure;
  } elsif ( $to_units eq '-' ) {
    $result = $measure;
  }

  if ( $result != -9999 && $to_units eq 'm' ) {

    if ( $result > 0.0 ) {
      $result = -$result; # Express depth as negative meters.
    }
  }

  return $result;
}



###############################################################################
1; # Required last executed statement of module. Do not add code below here!

#!/usr/bin/perl -wT

##############################################################################
# PURPOSE: castnet.pm - CASTNET coverage metadata used by castnetserver and
#          /data/CASTNET/reformat_quantaq and data/CASTNET/reformat_purpleair.
# NOTES:   Usage:
#          use strict;
#          BEGIN { push @INC, '/data/Apache/cgi-bin'; }
#          use castnet; # Use my $coverage_metadata_reference = castnet::get_coverage_metadata_reference();
#          pop @INC; # Remove local directory from path.
#          package main;
#          ...
# HISTORY: 2026-07-04 plessel.todd@epa.gov
# STATUS:  unreviewed tested
##############################################################################

use strict;
package castnet;

##############################################################################

# column is the 1-based column of the measure in purpleair/quantaq_web.csv.

my @coverage_metadata  = (
  {
    name => 'purpleair.hour_valid',
    units => '%',
    description => 'Percentage of hour with valid measurements.',
    begin_date => '2025-07-01',
    end_date   => '',
    column => 4,
    min => 0.0,
    max => 100.0
  },
  {
    name => 'purpleair.humidity',
    units => '%',
    description => 'Relative humidity.',
    begin_date => '2025-07-01',
    end_date   => '',
    column => 6,
    min => 0.0,
    max => 100
  },
  {
    name => 'purpleair.pm25_corrected',
    units => 'ug/m3',
    description => 'Scaled humidity-corrected particulate matter not more than 2.5 microns in diameter.',
    begin_date => '2025-07-01',
    end_date   => '',
    column => 5,
    min => 0.0,
    max => 1000.0
  },
  {
    name => 'purpleair.temperature',
    units => 'C',
    description => 'Temperature.',
    begin_date => '2025-07-01',
    end_date   => '',
    column => 7,
    min => -50.0,
    max => 50.0
  },
  {
    name => 'quantaq.co',
    units => 'ppm',
    description => 'Carbon monoxide concentration.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 9,
    min => 0.0,
    max => 1000.0,
    scale => 0.001 # ppb to ppm
  },
  {
    name => 'quantaq.hour_valid',
    units => '%',
    description => 'Percentage of hour with valid measurements.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 4,
    min => 0.0,
    max => 100.0
  },
  {
    name => 'quantaq.humidity',
    units => '%',
    description => 'Relative humidity.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 17,
    min => 0.0,
    max => 100
  },
  {
    name => 'quantaq.no',
    units => 'ppb',
    description => 'Nitric oxide concentration.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 11,
    min => 0.0,
    max => 1e6
  },
  {
    name => 'quantaq.no2',
    units => 'ppb',
    description => 'Nitrogen dioxide concentration.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 13,
    min => 0.0,
    max => 1e6
  },
  {
    name => 'quantaq.o3',
    units => 'ppb',
    description => 'Ozone concentration.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 15,
    min => 0.0,
    max => 1e6
  },
  {
    name => 'quantaq.pm1',
    units => 'ug/m3',
    description => 'Particulate matter not more than 1 micron in diameter.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 6,
    min => 0.0,
    max => 1000.0
  },
  {
    name => 'quantaq.pm10',
    units => 'ug/m3',
    description => 'Particulate matter not more than 10 microns in diameter.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 7,
    min => 0.0,
    max => 1000.0
  },
  {
    name => 'quantaq.pm25',
    units => 'ug/m3',
    description => 'Particulate matter not more than 2.5 microns in diameter.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 6,
    min => 0.0,
    max => 1000.0
  },
  {
    name => 'quantaq.temperature',
    units => 'C',
    description => 'Temperature.',
    begin_date => '2025-09-17',
    end_date   => '',
    column => 19,
    min => -50.0,
    max => 50.0
  }
);

###############################################################################

sub get_coverage_metadata_reference { return \@coverage_metadata; }

1; # Required last executed statement of module. Do not add code below here!


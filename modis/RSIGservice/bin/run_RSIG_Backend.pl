#!/usr/bin/env perl
# This script runs the RSIG backend used for EPA Web Coverage Server.
# Output from RSIG appears on stdout.
#
# If we get a USR1 signal from sge, then we took too long, and it's time
# to stop. This means killing child processes and cleaning up.
#
# Design Notes:
# Be very careful how you use "die", because the script should clean up
# before exiting.
#
# Must not use MODAPS/MTVS perl modules, since the MODAPS environment has
# not necessarily been set.

use strict;
use warnings;

use File::Path qw(remove_tree make_path);
use POSIX;
use File::Spec;

use FindBin;
use lib $FindBin::Bin;
use RSIG_Common;

my $MODAPSROOT = $RSIG_Common::MODAPSROOT;
my $HOST = $RSIG_Common::HOST;
my $JOB_ID = $RSIG_Common::JOB_ID;

RSIG_Common::LOG( "running $HOST\:\:$JOB_ID");

# setup signal handlers
my $SIGNAME = '';
my $SLEEP_TIME = 2;  # seconds to sleep while waiting for forked child
my $TIMED_OUT = 0;
$SIG{TERM} = $SIG{INT} = $SIG{USR1} = sub { $SIGNAME = shift; ++$TIMED_OUT; };

# working directory
my $WORKING_DIR = $ENV{SCRATCH};
chdir($WORKING_DIR);

my $ALLDATA = "/tis/modaps/allData";

#--------------------------------------------------------------
# need to fork so we can keep watching for signals while
# waiting for child to complete.
#--------------------------------------------------------------
my $result = -1;
if (! $TIMED_OUT)
{
  my $child = fork;
  if ($child)                      # proud parent
  {
    # check for signals, and check if child is done...
    while (1)
    {
      # if we're being shut down, let child know first...
      if ($TIMED_OUT && $child)
      {
        kill 9, $child;
      }

      # wait a bit so child can make progress
      sleep($SLEEP_TIME);

      # is child done yet?
      my $wpid = POSIX::waitpid(-1, POSIX::WNOHANG);
      if (defined $wpid)
      {
        if ($wpid == $child)
        {
          $result = $?;
          $child = undef;
          last;
        }
	elsif ($wpid == -1)
        {
          # no child
          $child = undef;
          last;
	}
      }
    }
  }
  elsif (defined $child)           # child here 
  {
    my $result = run_RSIG(@ARGV);
    $result = 1 if $result; # make sure status <= 255, since it is shifted <<8
    exit $result;
    # IMPORTANT. THE CHILD HAS COMPLETED ITS TASK AND
    # MUST NEVER GET BEYOND THIS POINT.
  }
  else                             # really bad karma 
  {
    LOG("Failed fork : $!\n");
    $result = -1;
  }
}
else                    # error
{
  if ($TIMED_OUT)
  {
    LOG("Caught signal $SIGNAME -- exiting\n");
  }
}

# clean up done automatically by SLURM

# exit with whatever RSIG exited with, so that errors get propogated
# to the front end...
exit($result);

################################## ROUTINES ##################################
#------------------------------------------------------------------------------
# Run the RSIG applications after parsing command line params and
# searching for input files...
#------------------------------------------------------------------------------
sub run_RSIG
{
 my $result = 0;

 my $request = RSIG_Common::check_args(@_);

 if ( ! $request ) {
   $result = 0;
 } elsif ( $request eq 'getcapabilities' ) {
   RSIG_Common::print_capabilities();
   $result = 1;
 } elsif ( $request eq 'describecoverage' ) {
   RSIG_Common::print_coverage_description();
   $result = 1;
 } elsif ( $request eq 'getcoverage' ) {
   # get the text of the getcoverage command
   RSIG_Common::compute_time_range();
   RSIG_Common::get_modis_data($WORKING_DIR, $ALLDATA);
   my $c = RSIG_Common::construct_command();
   $result = system( $c ) ? 0 : 1; # Run either the original or altered command.
 } elsif ( $request eq 'getmetadata' ) {
   # get the text of the getcoverage command
   RSIG_Common::compute_time_range();
   RSIG_Common::get_modis_data($WORKING_DIR, $ALLDATA);
   my $c = RSIG_Common::construct_command();
   # get a new command to print $c from above with other metadata.
   $c = RSIG_Common::print_metadata($c);
   $result = system( $c ) ? 0 : 1; # Run either the original or altered command.
 }

 $result =  1 - $result; # UNIX: zero is success, non-zero is failure.
 RSIG_Common::LOG("system result=$result\n");
 return $result;
}

0;

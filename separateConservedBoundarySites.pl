#!/usr/bin/env perl
#==============================================================================
#     FILENAME:  separateConservedBoundarySites.pl
#      VERSION:  1.0
#      CREATED:  2018-02-02 11:15:47
#     REVISION:  ---
#	    AUTHOR:  Zhijun Han, hangeneral@126.com
#==============================================================================

use strict;
use warnings;

die "Usage: perl $0 diff.boundary overlap.boundary output\n" unless @ARGV == 3;

my ($fdiff, $foverlap, $output) = @ARGV;
my (%boundary, @t);
open my $in, '<', "$fdiff" or die "Can't open file $fdiff. $!\n";
my $head = <$in>;
while (<$in>) {
	@t = split;
	$t[0] =~ s/-.*$//;
	$boundary{$t[0]} = $t[-1];
}
close $in;

open $in, '<', "$foverlap" or die "Can't open file $foverlap. $!\n";
open my $to, '>', "$output" or die "Can't create file $output. $!\n";
while (<$in>) {
	@t = split;
	if ($t[0] =~ /,/) { # conserved sites
		foreach my $x (split /,|\|/, $t[1]) {
			$x =~ s/:\w+$//;
			next unless exists $boundary{$x};
			print {$to} join("\t", $x, $boundary{$x}, 'conserved'), "\n";
			delete $boundary{$x};
		}
	}
	else { # non-conserved sites
		foreach my $x (split /,/, $t[1]) {
			$x =~ s/:\w+$//;
			next unless exists $boundary{$x};
			print {$to} join("\t", $x, $boundary{$x}, 'unique'), "\n";
			delete $boundary{$x};
		}
	}
}
close $in;
close $to;

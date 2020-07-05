#!/usr/bin/env perl
#==============================================================================
#     FILENAME:  getBiasedUniqueSites.pl
#      VERSION:  1.0
#      CREATED:  2018-05-01 16:40:10
#     REVISION:  ---
#	    AUTHOR:  Zhijun Han, hangeneral@126.com
#==============================================================================

use strict;
use warnings;


die "perl $0 diff.boundary diff.boundary.separated bias_cutoff region_size sites.bed regions.bed\n" unless @ARGV == 6;

my ($fbias, $ftype, $cutoff, $size, $fsites, $fregions) = @ARGV;

my (%unique, @t);
open my $in, '<', "$ftype" or die "Can't open file $ftype. $!\n";
while (<$in>) {
	@t = split;
	$unique{$t[0]} = undef if $t[2] eq 'unique';
}
close $in;

open my $to1, '>', "$fsites" or die "Can't create file $fsites. $!\n";
open my $to2, '>', "$fregions" or die "Can't create file $fregions. $!\n";

my (@s, $mean);
my $count = 0;
open $in, '<', "$fbias" or die "Can't open file $fbias. $!\n";
my $head = <$in>;
while (<$in>) {
	@t = split;
	next if $t[-1] < $cutoff;
	@s = split /:|-/, $t[0];
	next unless exists $unique{"$s[0]:$s[1]"};
	print {$to1} join("\t", @s, $count++, $t[3], '.'), "\n"; # BED6 format
	$mean = ($s[1] + $s[2]) / 2;
	print {$to2} join("\t", $s[0], $mean - $size, $mean + $size), "\n";
}
close $in;
close $to1;
close $to2;

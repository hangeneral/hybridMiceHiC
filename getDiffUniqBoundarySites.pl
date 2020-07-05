#!/usr/bin/env perl
#==============================================================================
#     FILENAME:  getDiffUniqBoundarySites.pl
#      VERSION:  1.0
#      CREATED:  2018-03-02 17:44:30
#     REVISION:  ---
#	    AUTHOR:  Zhijun Han, hangeneral@126.com
#==============================================================================

use strict;
use warnings;


my ($foverlap, $output, $bedGraphs, $labels) = @ARGV;

my @bedGraphs = split /,/, $bedGraphs;
my @labels = split /,/, $labels;
die "Number of bedGraphs differ with labels!\n" if $#bedGraphs != $#labels;

my (%unique, @t);
open my $in, '<', "$foverlap" or die "Can't open file $foverlap. $!\n";
while (<$in>) {
	@t = split;
	next if $t[0] =~ /,/; # overlapped
	for (split /,/, $t[1]) {
		@t = split /:/, $_;
		$unique{$t[0]}{$t[1]} = undef;
	}
}
close $in;

for my $i (0 .. $#bedGraphs) {
	open $in, '<', "$bedGraphs[$i]" or die "Can't open file $bedGraphs[$i]. $!\n";
	while (<$in>) {
		@t = split;
		if (exists $unique{$t[0]}{$t[1]}) {
			$unique{$t[0]}{$t[1]}[$i] = $t[3];
		}
	}
	close $in;
}
my $binsize = $t[2] - $t[1];

my @array;
while (my ($chr, $ref1) = each %unique) {
	while (my ($pos, $ref2) = each %$ref1) {
		push @array, ["$chr:$pos-" . ($pos + $binsize), @$ref2, abs( $ref2->[0] - $ref2->[1] )]
	}
}
@array = sort { $b->[-1] <=> $a->[-1] } @array;

open my $to, '>', "$output" or die "Can't create file $output. $!\n";
print {$to} join( "\t", "position", @labels, "Bias" ), "\n";
foreach my $x (@array) {
    print {$to} join( "\t", @$x ), "\n";
}
close $to;

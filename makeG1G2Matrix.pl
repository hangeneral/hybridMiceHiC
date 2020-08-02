#!/usr/bin/env perl
#==============================================================================
#     FILENAME:  makeG1G2Matrix.pl
#      VERSION:  1.0
#      CREATED:  2018-09-19 15:28:49
#     REVISION:  ---
#	    AUTHOR:  Zhijun Han, hangeneral@126.com
#==============================================================================

use strict;
use warnings;


my $finput = '../paired.pup.txt';
my $output = 'G1G2.count.txt';

my ($in, $to, @t, %count);

my %type = (
	G1 => 'p',
	G2 => 'm',
);

open $in, '<', "$finput" or die "Can't open file $finput. $!\n";
while (<$in>) {
	@t = split;
	next unless exists $type{$t[2]} and exists $type{$t[5]};
	$t[0] =~ s/^chr/$type{$t[2]}/;
	$t[3] =~ s/^chr/$type{$t[5]}/;
	$count{$t[0]}{$t[3]}++;
	$count{$t[3]}{$t[0]}++;
}
close $in;

open $to, '>', "$output" or die "Can't create file $output. $!\n";
my @chrs = sort keys %count;
print {$to} join("\t", 'chr', @chrs, 'Genotype'), "\n";
for my $i (@chrs) {
	@t = ($i, @{$count{$i}}{@chrs});
	if ($i =~ /^p/) {
		push @t, 'Paternal';
	}
	else {
		push @t, 'Maternal';
	}
	print {$to} join("\t", @t), "\n";
}
close $to;


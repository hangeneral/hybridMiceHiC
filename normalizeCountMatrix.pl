#!/usr/bin/env perl 
#===============================================================================
#
#      COMPANY:  Group of Epigenome Biology, PICB
#      VERSION:  1.0
#      CREATED:  09/15/2017 10:31:24 AM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use List::Util qw(sum);


my ($input, $output) = @ARGV;

my (@count, @t);


my	$in_file_name = "$input";		# input file name
open  my $in, '<', $in_file_name
    or die  "$0 : failed to open  input file '$in_file_name' : $!\n";

my $head = <$in>;
while ( <$in> ) {
    @t = split /,/;
    $count[$. - 1] = sum(@t[1 .. $#t-1]); # the last column is genotype
}

close  $in
    or warn "$0 : failed to close input file '$in_file_name' : $!\n";

$in_file_name = "$input";		# input file name
open $in, '<', $in_file_name
    or die  "$0 : failed to open  input file '$in_file_name' : $!\n";

my	$to_file_name = "$output";		# output file name
open  my $to, '>', $to_file_name
    or die  "$0 : failed to open  output file '$to_file_name' : $!\n";

$head = <$in>;
print {$to} $head;
while ( <$in> ) {
    @t = split /,/;
    for (1 .. $#t-1) {
        $t[$_] = ($count[$.-1] + $count[$_]) * $t[$_] / $count[$.-1] / $count[$_];
    }
    print {$to} join(",", @t);
}

close  $to
    or warn "$0 : failed to close output file '$to_file_name' : $!\n";

close  $in
    or warn "$0 : failed to close input file '$in_file_name' : $!\n";




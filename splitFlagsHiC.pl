#!/usr/bin/env perl 
#===============================================================================
#
#      COMPANY:  Group of Epigenome Biology, PICB
#      VERSION:  1.0
#      CREATED:  08/07/2017 03:09:32 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

my $prefix = $ARGV[0];
$prefix =~ s/.allele_flagged.txt$//;

my	$G1_file_name = "$prefix.G1.txt";		# output file name
open  my $G1, '>', $G1_file_name
    or die  "$0 : failed to open  output file '$G1_file_name' : $!\n";
my	$G2_file_name = "$prefix.G2.txt";		# output file name
open  my $G2, '>', $G2_file_name
    or die  "$0 : failed to open  output file '$G2_file_name' : $!\n";
my	$UA_file_name = "$prefix.UA.txt";		# output file name
open  my $UA, '>', $UA_file_name
    or die  "$0 : failed to open  output file '$UA_file_name' : $!\n";
my	$MIX_file_name = "$prefix.MIX.txt";		# output file name
open  my $MIX, '>', $MIX_file_name
    or die  "$0 : failed to open  output file '$MIX_file_name' : $!\n";

my %count = (
    UA => 0, 
    G1 => 0, 
    G2 => 0, 
    MIX => 0, 
);
my %value = (
    UA => 0, 
    CF => 0, 
    G1 => 1, 
    G2 => 3, 
);
my (@t);
while ( <> ) {
    @t = split;
    if ($value{$t[2]} + $value{$t[5]} == 0) {
        $count{UA}++;
        print {$UA} join("\t", @t), "\n";
    }
    elsif ($value{$t[2]} + $value{$t[5]} == 4) {
        $count{MIX}++;
        print {$MIX} join("\t", @t), "\n";
    }
    elsif ($value{$t[2]} + $value{$t[5]} > 2) {
        $count{G2}++;
        print {$G2} join("\t", @t), "\n";
    }
    else {
        $count{G1}++;
        print {$G1} join("\t", @t), "\n";
    }
}

print "Total\tG1\tG2\tUA\tMIX\n";
print "$.\t$count{G1}\t$count{G2}\t$count{UA}\t$count{MIX}\n";

close  $UA
    or warn "$0 : failed to close output file '$UA_file_name' : $!\n";
close  $G2
    or warn "$0 : failed to close output file '$G2_file_name' : $!\n";
close  $G1
    or warn "$0 : failed to close output file '$G1_file_name' : $!\n";
close  $MIX
    or warn "$0 : failed to close output file '$MIX_file_name' : $!\n";

#!/usr/bin/env perl 
#===============================================================================
#
#      COMPANY:  Group of Epigenome Biology, PICB
#      VERSION:  1.0
#      CREATED:  09/16/2017 02:25:21 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;



open  my $toG1, "| samtools view -bt $ENV{mm10}/mm10.chromSizes -o G1.reads.bam";
open  my $toG2, "| samtools view -bt $ENV{mm10}/mm10.chromSizes -o G2.reads.bam";
open  my $toUA, "| samtools view -bt $ENV{mm10}/mm10.chromSizes -o UA.reads.bam";

my (@t);
while ( <> ) {
    @t = split;
    if ($t[-1] =~ /G2$/) {
        print {$toG2} $_;
    }
    elsif ($t[-1] =~ /G1$/) {
        print {$toG1} $_;
    }
    else {
        print {$toUA} $_;
    }
}

close  $toUA;
close  $toG2;
close  $toG1;




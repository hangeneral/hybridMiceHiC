#!/usr/bin/env perl 
#===============================================================================
#
#      COMPANY:  Group of Epigenome Biology, PICB
#      VERSION:  1.0
#      CREATED:  12/29/2017 07:23:03 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

# call peaks from a bedGraph file
# input: TAD boundary bedGraph file
# output: TAD boundaries in BED format
# delta: value to define a peak

die "Usage: perl $0 boundary.bedGraph output.bed delta\n" unless @ARGV == 3;
peakDetect(@ARGV);

sub peakDetect {
    my ( $input, $output, $delta ) = @_;
    open my $in, '<', $input or die "Can't open $input. $!\n";
    open my $to, '>', $output;
    my ( $flag, @t );
    my $chr = '';
    my $max;
    my $pos;
    my $min;
    my $count = 1;
    my $binsize;

    while (<$in>) {
        @t = split;
        if ( $t[0] ne $chr ) {
            $chr     = $t[0];
            $binsize = $t[2] - $t[1];
            $max     = $t[3];
            $pos     = $t[1];
            $min     = $t[3];
            $flag    = 0;               # 1 for search peaks
        }
        else {
            if ( $max < $t[3] ) {
                $max = $t[3];
                $pos = $t[1];
            }
            if ( $min > $t[3] ) {
                $min = $t[3];
            }
            if ($flag) {
                if ( $max - $t[3] > $delta ) {
                    print {$to} join( "\t",
                        $chr, $pos, $pos + $binsize,
                        $count, $max, '.' ),
                      "\n";
                    $count++;
                    $min  = $t[3];
                    $flag = 0;
                }
            }
            else {
                if ( $t[3] - $min > $delta ) {
                    $max  = $t[3];
                    $pos  = $t[1];
                    $flag = 1;
                }
            }
        }
    }
    close $in;
    close $to;
}


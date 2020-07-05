#!/usr/bin/env perl 
#===============================================================================
#
#      COMPANY:  Group of Epigenome Biology, PICB
#      VERSION:  1.0
#      CREATED:  01/08/2018 09:11:15 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use List::MoreUtils 'pairwise';


my ($file1, $file2, $output) = @ARGV;

my (%score1, @t, @score1, @score2, $median1, $median2);
open my $in, '<', $file1 or die "Can't open $file1. $!\n";
my ($chr, @output, %start1, $start1, $start2) = ('');
while (<$in>) {
    @t = split;
    push @score1, $t[3];
    push @{$score1{$t[0]}}, $t[3];
    if ($t[0] ne $chr) {
        $start1{$t[0]} = $t[1];
        $chr = $t[0];
    }
}
close $in;
my $binsize1 = $t[2] - $t[1];
@score1 = sort {$a <=> $b} @score1;
$median1 = $score1[$#score1 / 2];
open $in, '<', $file2 or die "Can't open $file2. $!\n";
while (<$in>) {
    @t = split;
    push @score2, $t[3];
}
close $in;
my $binsize = $t[2] - $t[1];
die "Binsize of $file1 and $file2 differ!\n" unless $binsize1 == $binsize;
@score2 = sort {$a <=> $b} @score2;
$median2 = $score2[$#score2 / 2];
my $diff = $median1 - $median2;
$chr = '';
open $in, '<', $file2 or die "Can't open $file2. $!\n";
open my $to, '>', $output;
while (<$in>) {
    @t = split;
    if ($t[0] ne $chr) {
        if ($chr) {
            @score1 = @{$score1{$chr}};
            $start1 = $start1{$chr};
            if ($start1 < $start2) {
                while ($start1 < $start2) {
                    unshift @score2, 0;
                    $start2 -= $binsize;
                }
            }
            if ($start1 > $start2) {
                while ($start1 < $start2) {
                    unshift @score1, 0;
                    $start1 -= $binsize;
                }
            }
            if ($#score1 < $#score2) {
                push @score1, 0 until $#score1 == $#score2;
            }
            if ($#score1 > $#score2) {
                push @score2, 0 until $#score1 == $#score2;
            }
            @output = pairwise {$a - $b - $diff} @score1, @score2;
            for (@output) {
                print {$to} join("\t", $chr, $start1, $start1 + $binsize, $_), "\n";
                $start1 += $binsize;
            }
        }
        $chr = $t[0];
        $start2 = $t[1];
        @score2 = ($t[3]);
    }
    else {
        push @score2, $t[3];
    }
}
close $in;
@score1 = @{$score1{$chr}};
$start1 = $start1{$chr};
if ($start1 < $start2) {
    while ($start1 < $start2) {
        unshift @score2, 0;
        $start2 -= $binsize;
    }
}
if ($start1 > $start2) {
    while ($start1 < $start2) {
        unshift @score1, 0;
        $start1 -= $binsize;
    }
}
if ($#score1 < $#score2) {
    push @score1, 0 until $#score1 == $#score2;
}
if ($#score1 > $#score2) {
    push @score2, 0 until $#score1 == $#score2;
}
@output = pairwise {$a - $b - $diff} @score1, @score2;
for (@output) {
    print {$to} join("\t", $chr, $start1, $start1 + $binsize, $_), "\n";
    $start1 += $binsize;
}
close $to;


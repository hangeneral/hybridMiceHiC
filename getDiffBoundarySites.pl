#!/usr/bin/env perl
#==============================================================================
#     FILENAME:  getDiffBoundarySites.pl
#      VERSION:  1.0
#      CREATED:  2018-02-01 15:51:26
#     REVISION:  ---
#	    AUTHOR:  Zhijun Han, hangeneral@126.com
#==============================================================================

use strict;
use warnings;

die
    "Usage: perl $0 boundaryFile1,boundaryFile2 bedGraphFile1,bedGraphFile2 label1,label2 output\n"
    unless @ARGV == 4;

my ( $fboundary, $fbedGraph, $labels, $output ) = @ARGV;

my @fboundary = split /,/, $fboundary;
my @fbedGraph = split /,/, $fbedGraph;
my @labels    = split /,/, $labels;

die "Number of boundary files, bedGraph files and labels must equal to 2!\n"
    unless $#fboundary == 1
    and $#fbedGraph == 1
    and $#fboundary == 1;

my ( %boundary, @t );
foreach my $i ( 0 .. $#fboundary ) {
    open my $in, '<', "$fboundary[$i]"
        or die "Can't open file $fboundary[$i]. $!\n";
    while (<$in>) {
        @t = split;
        $boundary{ $t[0] }{ $t[1] } = undef;
    }
    close $in;
}
foreach my $i ( 0 .. $#fbedGraph ) {
    open my $in, '<', "$fbedGraph[$i]"
        or die "Can't open file $fbedGraph[$i]. $!\n";
    while (<$in>) {
        @t = split;
        next unless exists $boundary{ $t[0] }{ $t[1] };
        push @{ $boundary{ $t[0] }{ $t[1] } }, $t[3];
    }
    close $in;
}
my $size = $t[2] - $t[1];    # bin size
my @array;
while ( my ( $chr, $ref ) = each %boundary ) {
    while ( my ( $pos, $ref2 ) = each %$ref ) {
        push @array,
            [
            "$chr:$pos-" . ($pos + $size),
            @$ref2, abs( $ref2->[0] - $ref2->[1] )
            ];
    }
}
@array = sort { $b->[-1] <=> $a->[-1] } @array;

open my $to, '>', "$output" or die "Can't create file $output. $!\n";
print {$to} join( "\t", "position", @labels, "Bias" ), "\n";
foreach my $x (@array) {
    print {$to} join( "\t", @$x ), "\n";
}
close $to;

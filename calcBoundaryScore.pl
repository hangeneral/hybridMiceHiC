#!/usr/bin/env perl
#===============================================================================
#
#      COMPANY:  Group of Epigenome Biology, PICB
#      VERSION:  1.0
#      CREATED:  12/26/2017 02:13:15 PM
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use List::MoreUtils qw(pairwise);
use List::Util qw(sum);
use Parallel::ForkManager;

die "Usage: perl $0 genome output binsize short_cutoff maximal_cutoff\n"
    unless @ARGV == 5;
my ( $genome, $output, $binsize, $short, $maximal ) = @ARGV;

my $cpu = 20;
my ( $chromSizesFile, @chrs, %chrs );
_getGenomeInfo();

my $pm = Parallel::ForkManager->new($cpu);
$pm->run_on_finish(
    sub {
        my ( $pid, $exit_code, $indent, $exit_signal ) = @_;
        die "Failed to calculate boundrary score for '$indent'. $exit_signal\n"
          unless $exit_code == 0;
    }
);
my $pid = $$;
for my $chr (@chrs) {
    $pm->start($chr) and next;
    _calcBoundaryScore( "$chr.realintra", "$chr.$pid", $binsize, $short,
        $maximal );
    $pm->finish(0);
}
$pm->wait_all_children;

unlink $output if -e $output;
for my $chr (@chrs) {
    system("cat $chr.$pid >> $output") == 0
      or die "CMD: cat $chr.$pid >> $output failed. $!\n";
    unlink "$chr.$pid";
}

sub _calcBoundaryScore {
    my ( $input, $output, $binsize, $short, $maximal ) = @_;
    open my $in, '<', $input or die "Can't open $input. $!\n";
    my ( @pets, $within, @t, $cross, $valid, @boundary, @within, @cross );

    # read in the first line to initialize
    chomp( my $line = <$in> );
    @t = split /\s+/, $line;    # site1 chr site2
    my $chr      = $t[1];
    my $halfbin  = int( $binsize / 2 );
    my $firstBin = $t[0] - $t[0] % $binsize;   # start of the first non-zero bin
    my $locus    = $firstBin + $halfbin;       # exact locus position
    my $start = $locus - int( $maximal / 2 ); # domain start for each locus
    my $end   = $start + $maximal;            # domain end: locus + $maximal / 2
    push @pets, [ @t[ 0, 2 ] ] if $t[2] - $t[0] > $short and $t[2] - $t[0] < $maximal;

    while (<$in>) {
        @t = split;
        next if $t[2] - $t[0] < $short;        # remove short distance PETs
        next if $t[2] - $t[0] > $maximal;      # remove long distance PETs
        if ( $t[0] > $end ) {                  # process @pets for this locus
            if (@pets) {    # calculate valid and cross PETs number
                ( $cross, $valid ) = ( 0, 0 );
                for (@pets) {
                    if ( $_->[1] < $end ) {
                        $valid++;
                        $cross++ if $_->[0] < $locus and $_->[1] > $locus;
                    }
                }
                $within = $valid - $cross;
                push @within, $within;
                push @cross, $cross;

                # boundary score = log2(within / new-cross)
#                push @boundary,
#                  log( ( $within + 1 ) / ( $cross + 1 ) / log(2) );
            }
            else {    # empty @pets
                push @within, 0;
                push @cross, 0;
#                push @boundary, 0;
            }

            # upadte values
            $start += $binsize;
            $end   += $binsize;
            $locus += $binsize;
            while (@pets) {    # @pets are sorted by first elements
                if ( $pets[0][0] < $start ) {
                    shift @pets;
                }
                else {
                    last;
                }
            }
            redo;
        }
        else {                 # just record it
            push @pets, [ @t[ 0, 2 ] ];
        }
    }
    close $in;
    while (@pets) {            # $in reached the end
        ( $cross, $valid ) = ( 0, 0 );
        for (@pets) {
            if ( $_->[1] < $end ) {
                $valid++;
                $cross++ if $_->[0] < $locus and $_->[1] > $locus;
            }
        }
        $within = $valid - $cross;
        push @within, $within;
        push @cross, $cross;

        # boundary score = log2(within / new-cross)
#        push @boundary, log( ( $within + 1 ) / ( $cross + 1 ) / log(2) );

        # upadte values
        $start += $binsize;
        $end   += $binsize;
        $locus += $binsize;
        while (@pets) {    # @pets are sorted by first elements
            if ( $pets[0][0] < $start ) {
                shift @pets;
            }
            else {
                last;
            }
        }
    }

    my $meanWithin = sum(@within) / scalar(@within);
    my $meanCross = sum(@cross) / scalar(@cross);
    @within = map {$_ + $meanWithin} @within;
    @cross = map {$_ + $meanCross} @cross;
    @boundary = pairwise {log($a / $b) / log(2)} @within, @cross;
    # filter negative values since they largely located in unmappable regions
#    @boundary = map { $_ < 0 ? 0 : $_ } @boundary;
    open my $to, '>>', $output;
    for (@boundary) {
        print {$to} join( "\t",
            $chr, $firstBin,
            $firstBin + $binsize,
            sprintf( "%.2f", $_ ) ),
          "\n";
        $firstBin += $binsize;
        last if $firstBin + $binsize > $chrs{$chr};
    }
    close $to;
}

sub _getGenomeInfo {

    # search and parse chromSizes file to get chromosomes
    die "-genome is NEEDED!\n" unless $genome;
    if ( -e "$ENV{BOWTIE_INDEXES}/$genome.chromSizes" ) {
        $chromSizesFile = "$ENV{BOWTIE_INDEXES}/$genome.chromSizes";
    }
    elsif ( -e "$ENV{$genome}/$genome.chromSizes" ) {
        $chromSizesFile = "$ENV{$genome}/$genome.chromSizes";
    }
    else {
        system("fetchChromSizes $genome > $genome.chromSizes") == 0
          or die "Failed to fetchChromSizes for '$genome'.\n";
        $chromSizesFile = "$genome.chromSizes";
    }
    my $in_file_name = $chromSizesFile;    # input file name
    open my $in, '<', $in_file_name
      or die "$0 : failed to open  input file '$in_file_name' : $!\n";
    %chrs = map { (split)[ 0, 1 ] } <$in>;
    for ( keys %chrs ) {
        delete $chrs{$_} if /chr.*_/ or /chrM/i or /chrY/i;
    }
    close $in
      or warn "$0 : failed to close input file '$in_file_name' : $!\n";
    unlink $chromSizesFile if $chromSizesFile eq "$genome.chromSizes";
    @chrs = sort { $a cmp $b } keys %chrs;
}

#!/usr/bin/env perl

use warnings;
use strict;

# debug with stacktrace:
use Carp;
$SIG{__DIE__} = \&Carp::confess;


######################################################
#
# simple library functions:
#
my $LOGFH =\*STDERR;
my $OUTFH =\*STDOUT;

sub errorX($) {
    confess "ERROR: " . $_[0] . "\n";
}

sub logX($) {
    print $LOGFH "INFO: " . $_[0] . "\n";
}

sub checkFile($;$) {
    my ($file,$label) = @_;
    return if(-f $file);
    errorX("Can't find" . (defined($label) ? " $label" : "") . " file: '$file'");
}


#####################################

# optionally point to a samtools binary:
#
my $samtoolsBin=$ENV{'SAMTOOLS'};

if(! defined($samtoolsBin)) {
   $samtoolsBin = "samtools"; 
}

if(scalar(@ARGV)!=1) {
  die("usage: $0 bamfile");
}

my $bamfile=$ARGV[0];

checkFile($bamfile,"bam");
checkFile($bamfile.".bai","bam index");



my $cmd1="$samtoolsBin idxstats $bamfile";
open(my $FP1,"$cmd1 |");

my %chrom;
my @chroms;

while(<$FP1>) {
    my @X =split(/\t/);
    next if($X[0] eq "*");
    $chrom{$X[0]} = [ $X[1] , $X[2] ];
    push @chroms, $X[0];
}

close($FP1);

#my $cmd2="$samtoolsBin view -F 4 -s 0.1 $bamfile";

# pass 0 is a subsampled approximation of read length,
# if that fails, then pass 1 runs the exact computation
#
my $length = 0;
my $count = 0;
for my $pass ((0,1)) {
    my $cmd2;
    if($pass == 0) {
        $cmd2="$samtoolsBin view -F 4 -s 0.1 $bamfile";
    } else {
        $cmd2="$samtoolsBin view -F 4 $bamfile";
    }

    my $sid=open(my $FP2,"$cmd2 |");

    $length = 0;
    $count = 0;
    while(<$FP2>) {
        my @F = split(/\t/,$_,7);
        next unless($F[5] =~ /[0-9]+M/);
        $length += $1 while($F[5] =~ /([0-9]+)M/g);
        $count++;
        last if($count >= 200000);
    }
    kill('INT',$sid);
    close($FP2);

    last if($count > 100000);
    logX("Poor read length approximation results. Count: '$count' Rerunning exact estimate"); # rerun pass 1 for exact read length
}


my $avg_length = ($length/$count);

for my $c (@chroms) {
    next if($chrom{$c}->[0] < $avg_length);
    my $depth = $chrom{$c}->[1]*$avg_length / $chrom{$c}->[0];
    printf "%s\t%.3f\t%s\t%.3f\n",$c,$depth,$count,$avg_length;
}



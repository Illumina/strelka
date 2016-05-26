#!/usr/bin/perl

# Open feature file for input:
open FEATFILE, "scoring_features.txt" or die $!;

# Open all output files:
open SNVTRAIN, ">snv_training_data.csv" or die $!;
open INDELTRAIN, ">indel_training_data.csv" or die $!;

# Get headers and print to output files:
while (<FEATFILE>) {
    @content = split(/=/);
    chomp($content[1]);
    if ($content[0] =~ /snv_scoring_features/) {
	$SNVHEADER = sprintf("CHROM,POS,TYPE,%s,tag\n",$content[1]);
    } else {
	$INDELHEADER = sprintf("CHROM,POS,TYPE,%s,tag\n",$content[1]);
    }
}
print SNVTRAIN $SNVHEADER;
print INDELTRAIN $INDELHEADER;

# Work through input file:
while (<>) {
    next if (/^#/);
    @items = split;
    $gotEVSF=0;
    $starlingtype="SNV";
    foreach $infoentry (split(/;/, $items[7])) {
	@content = split(/=/, $infoentry);
	if ($content[0] eq "EVSF") {
	    $EVSF = $content[1];
	    $gotEVSF = 1;
	}
        if ($content[0] eq "CIGAR") { # if CIGAR is present it means Starling called the variant as an indel
            $starlingtype = "INDEL";
        }
    }
    next unless $gotEVSF; # skip entries with no EVSF field in INFO

    @filter = split(/;/, $items[6]);
    next if ("OffTarget" ~~ @filter); # Skip entries matching OffTarget in the filter field (for WES data).

    @truth = split(/:/, $items[9]);
    @query = split(/:/, $items[10]);
    next if ("NOCALL" ~~ @query); # skip entries with NOCALL in QUERY                     
    $label = $query[1]; # should be FP, TP or UNK
    ($label eq "TP" || $label eq "FP" || $label eq "UNK") or die "Found a line that is neither FP nor TP nor UNK:\n$_\n";

    # Get variant type:
    $type = $query[4];
    next unless ($type eq $starlingtype);

    # Print to appropriate output files according to type and geno:
    $out = join(",",$items[0],$items[1],$type,$EVSF,$label)."\n";

    if ($type eq "SNV") {
	print SNVTRAIN $out;
    }
    if ($type eq "INDEL") {
	print INDELTRAIN $out;
    }
}

# Clean up:
close SNVTRAIN;
close INDELTRAIN;

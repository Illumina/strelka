#!/usr/bin/perl

# Filters out some gvcf entries before passing it to hap.py. This script makes hard assumptions about the order of entries in the vcf.

open FEATFILE, ">scoring_features.txt" or die $!;
while (<>) {
    if (/^#/) {
	print;
	if (/scoring_features/) {
	    print FEATFILE
	}
	next;
    }

    @items = split;

    @filter = split(/;/, $items[6]);
    next if ("OffTarget" ~~ @filter); # Skip entries matching OffTarget in the filter field (for WES data).
    next if (/.*Conflict/ ~~ @filter); # Skip entries matching any type of conflict

    @format = split(/:/, $items[9]);
    next if ($format[0] !~ /.\/./); # Skip entries the do not have diploid genotype (hemizygotes)

    print;
}

close FEATFILE

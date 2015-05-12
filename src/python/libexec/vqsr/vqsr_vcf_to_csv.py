#!/usr/bin/env python
# coding=utf-8
#
# 20/11/2014
#
# Convert a Strelka VCF to CSV format, annotate TP and FP given a
# truth VCF and FP / AMBI bed files.
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import os
import sys
import argparse

scriptDir = os.path.abspath(os.path.dirname(__file__))
scriptName = os.path.basename(__file__)
workflowDir = os.path.abspath(os.path.join(scriptDir, "@THIS_RELATIVE_PYTHON_LIBDIR@"))
templateConfigDir = os.path.abspath(os.path.join(scriptDir, '@THIS_RELATIVE_CONFIGDIR@'))

sys.path.append(workflowDir)

import vqsr
import vqsr.features
from vqsr.tools.bamstats import bamStats
from vqsr.tools.bedintervaltree import BedIntervalTree

import pandas


def main():
    parser = argparse.ArgumentParser("vqsr precision/recall script")
    parser.add_argument("input", help="Strelka VCF file", nargs=1)

    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file CSV name")

    parser.add_argument("-t", "--truth", dest="truth", required=False, default="",
                        help="Truth VCF file to use.")

    parser.add_argument("-f", "--fp-regions", dest="fp_regions", required=False, default="",
                        help="Use specific regions for FPs only.")

    parser.add_argument("-a", "--ambiguous", dest="ambi", action='append',
                        help="Ambiguous region bed file(s) (places we don't want to label as FP).")

    parser.add_argument("--feature-table", dest="features", default="posandalleles",
                        choices=vqsr.features.FeatureSet.sets.keys(),
                        help="Select a feature table to output.")

    parser.add_argument("--bam", dest="bams", default=[], action="append",
                        help="pass one or more BAM files for feature table extraction")

    args = parser.parse_args()

    bams = []
    md = None
    for x in args.bams:
        bams.append(bamStats(x))

    if bams:
        bres = pandas.concat(bams).groupby("CHROM").mean()
        md = {}
        for x in bres.index:
            print >>sys.stderr, "Mean coverage on %s is %f" % (x, bres.loc[x]["COVERAGE"])
            md[x] = float(bres.loc[x]["COVERAGE"])*3.0

    fset = vqsr.features.FeatureSet.make(args.features)
    fset.setChrDepths(md)

    featuretable = fset.collect(args.input[0])
    featuretable["tag"] = ""

    if args.truth:
        fset2 = vqsr.features.FeatureSet.make("posandalleles")
        truth_alleles = fset2.collect(args.truth)
        truth_alleles["tag"] = "TP"
        featuretable["tag"] = "FP"
        featuretable = pandas.merge(featuretable, truth_alleles, how="outer", on=["CHROM", "POS", "REF", "ALT"],
                                    suffixes=(".query", ".truth"))

        featuretable["tag.truth"].fillna("", inplace=True)
        featuretable["tag.query"].fillna("", inplace=True)
        featuretable.loc[(featuretable["tag.query"] == "FP") & (featuretable["tag.truth"] == "TP"), "tag"] = "TP"
        featuretable.loc[(featuretable["tag.query"] == "") & (featuretable["tag.truth"] == "TP"), "tag"] = "FN"
        featuretable.loc[(featuretable["tag.query"] == "FP") & (featuretable["tag.truth"] == ""), "tag"] = "FP"

        to_keep = [x for x in list(featuretable) if not x.endswith(".query") and not x.endswith(".truth")]
        featuretable = featuretable[to_keep]

    fp = BedIntervalTree()
    if args.fp_regions:
        fp.addFromBed(args.fp_regions, "FP")

    if args.ambi:
        # can have multiple ambiguous BED files
        for aBED in args.ambi:
            fp.addFromBed(aBED, lambda xe: xe[4])

    if args.ambi or args.fp_regions:
        has_fp = (fp.count("FP") > 0) or (fp.count("fp") > 0 and args.ambi_fp)

        def relabeller(xx):
            if xx["tag"] == "TP" or xx["tag"] == "FN":
                return xx
            chrom = xx["CHROM"]
            start = xx["POS"]
            stop = xx["POS"] + len(xx["REF"])
            overlap = fp.intersect(chrom, start, stop)

            is_fp = False
            is_ambi = False

            classes_this_pos = set()
            for o in overlap:
                reason = o.value[0].upper()
                classes_this_pos.add(reason)
                if reason == "FP":
                    is_fp = True
                else:
                    is_ambi = True

            if is_fp:
                xx["tag"] = "FP"
            elif is_ambi:
                xx["tag"] = ",".join(list(classes_this_pos))
            elif not has_fp:
                # when we don't have FP regions, unk stuff becomes FP
                xx["tag"] = "FP"
            else:
                xx["tag"] = "UNK"

            return xx

        featuretable = featuretable.apply(relabeller, axis=1)

    featuretable.to_csv(args.output)


if __name__ == '__main__':
    main()

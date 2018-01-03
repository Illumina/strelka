#!/usr/bin/env python2
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""
Convert a Strelka somatic VCF to CSV format, annotate TP and FP given a
truth VCF and FP / ambiguous region bed files.
"""

__author__ = "Peter Krusche <pkrusche@illumina.com>"

import os
import sys

import pandas

scriptDir = os.path.abspath(os.path.dirname(__file__))
scriptName = os.path.basename(__file__)
workflowDir = os.path.abspath(os.path.join(scriptDir, "../lib"))

sys.path.append(workflowDir)

import evs
import evs.features
from evs.tools.bedintervaltree import BedIntervalTree


def parseArgs():

    import argparse
    parser = argparse.ArgumentParser(description="Converts somatic VCF to annotated CSV",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", help="Strelka VCF file", nargs=1)
    parser.add_argument("-o", "--output", required=True,
                        help="Output CSV filename for training data")
    parser.add_argument("--testSet", action='append', help="Chromosome (e.g. chr20) to hold out as test data (may be specified more than once; if omitted, all data will be used for training)")
    parser.add_argument("--testOutput", help="Output CSV filename for test data")
    parser.add_argument("--truth", help="Truth VCF file")
    parser.add_argument("--fp-regions", dest="fpRegionsFile",
                        help="Bed file indicating regions where variants that are not true can be labeled as false positives. Outside of these regions variants will be labeled as unknown.")
    parser.add_argument("--ambiguous", dest="ambiguousRegionsFiles", action='append',
                        help="Bed file conforming to the curium ambiguous region file format"
                             " (may be specified more than once)")
    parser.add_argument("--features", required=True,
                        choices=evs.features.FeatureSet.sets.keys(),
                        help="Select a feature table to output.")

    args = parser.parse_args()

    def checkFile(filename, label) :
        if not os.path.isfile(filename) :
            raise Exception("Can't find input %s file: '%s'" % (label,filename))

    def checkOptionalFile(filename, label) :
        if filename is None : return
        checkFile(filename, label)

    checkOptionalFile(args.truth,"truth")
    checkOptionalFile(args.fpRegionsFile,"false positive regions")
    if args.ambiguousRegionsFiles is not None :
        for ambiguousRegionsFile in args.ambiguousRegionsFiles :
            checkFile(ambiguousRegionsFile,"ambiguous regions")

    return args



def main():
    args = parseArgs()

    fset = evs.features.FeatureSet.make(args.features)

    featuretable = fset.collect(args.input[0])
    featuretable["tag"] = "FP" # If no truth set is specified, label all variants as FP. Useful for normal-normal.

    if args.truth:
        fset2 = evs.features.FeatureSet.make("posandalleles")
        truth_alleles = fset2.collect(args.truth)
        truth_alleles["tag"] = "TP"
        featuretable = pandas.merge(featuretable, truth_alleles, how="outer", on=["CHROM", "POS", "REF", "ALT"],
                                    suffixes=(".query", ".truth"))

        featuretable["tag.truth"].fillna("", inplace=True)
        featuretable["tag.query"].fillna("", inplace=True)
        featuretable.loc[(featuretable["tag.query"] == "FP") & (featuretable["tag.truth"] == "TP"), "tag"] = "TP"
        featuretable.loc[(featuretable["tag.query"] == "") & (featuretable["tag.truth"] == "TP"), "tag"] = "FN"
        featuretable.loc[(featuretable["tag.query"] == "FP") & (featuretable["tag.truth"] == ""), "tag"] = "FP"

        to_keep = [x for x in list(featuretable) if not x.endswith(".query") and not x.endswith(".truth")]
        featuretable = featuretable[to_keep]

    if args.ambiguousRegionsFiles or args.fpRegionsFile:
        #
        # 1. Load all false positive and ambiguous region information into labeledIntervals
        #
        labeledIntervals = BedIntervalTree()
        if args.fpRegionsFile:
            labeledIntervals.addFromBed(args.fpRegionsFile, "FP")

        if args.ambiguousRegionsFiles:
            # can have multiple ambiguous BED files
            for ambiguousRegionsFile in args.ambiguousRegionsFiles:
                labeledIntervals.addFromBed(ambiguousRegionsFile, lambda xe: xe[4])

        #
        # 2. Resolve all interaction rules between truth sets, fp and amiguous regions to produce a final labeling
        #
        areFPRegionsProvided = (labeledIntervals.count("FP") > 0) or (labeledIntervals.count("fp") > 0 and args.ambiguousRegionsFiles)

        def relabeller(xx):
            """
            Resolve various rules regarding how variants should interact with the fp and ambiguous regions they
            intersect.

            Rules:
            - All TP and FN calls are untouched -- even if they fall in a false positive or ambiguous region
            - Otherwise...
                - Any call intersecting an FP region is labeled as "FP", regardless of ambiguous region input
                - Any call intersecting an ambiguous region gets a comma separated list of all ambiguous region labels
                - Any call falling outside of an ambiguous or fp region will be labeled as:
                  - FP if no fp regions are given the ambiguous region file contains no false positive regions
                  - UNK otherwise.
            """
            if xx["tag"] == "TP" or xx["tag"] == "FN":
                return xx
            chrom = xx["CHROM"]
            start = xx["POS"]
            stop = xx["POS"] + len(xx["REF"])
            overlap = labeledIntervals.intersect(chrom, start, stop)

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
            elif not areFPRegionsProvided:
                # when we don't have FP regions, unk stuff becomes FP
                xx["tag"] = "FP"
            else:
                xx["tag"] = "UNK"

            return xx

        featuretable = featuretable.apply(relabeller, axis=1)

    if args.testSet is not None:
        if args.testOutput is not None:
            featuretable[featuretable["CHROM"].isin(args.testSet)].to_csv(args.testOutput)
        featuretable = featuretable[~featuretable["CHROM"].isin(args.testSet)]
    featuretable.to_csv(args.output)


if __name__ == '__main__':
    main()

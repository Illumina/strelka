#!/usr/bin/env python
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2016 Illumina, Inc.
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
workflowDir = os.path.abspath(os.path.join(scriptDir, "@THIS_RELATIVE_PYTHON_LIBDIR@"))

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
                        help="Output file CSV name")
    parser.add_argument("-t", "--truth", help="Truth VCF file")
    parser.add_argument("-f", "--fp-regions", dest="fp_regions",
                        help="Bed file indicating regions for FPs only.")
    parser.add_argument("-a", "--ambiguous", dest="ambi", action='append',
                        help="Bed file indicating ambiguous regions ie. places we don't want to label as FP"
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
    checkOptionalFile(args.fp_regions,"false positive regions")
    if args.ambi is not None :
        for ambiFile in args.ambi :
            checkFile(ambiFile,"ambiguous regions")

    return args



def main():
    args = parseArgs()

    fset = evs.features.FeatureSet.make(args.features)

    featuretable = fset.collect(args.input[0])
    featuretable["tag"] = ""

    if args.truth:
        fset2 = evs.features.FeatureSet.make("posandalleles")
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

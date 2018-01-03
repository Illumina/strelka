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
Given an EVS model evaluation, provide precision/recall
"""

__author__ = "Peter Krusche <pkrusche@illumina.com>"


import os
import sys
import random

import pandas
import numpy


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser("evs precision/recall script")
    parser.add_argument("inputs", help="Quality/tag CSV files", nargs="+")

    parser.add_argument("-q", "--quality-fields", dest="q", default="qual",
                        help="Select fields which give quality values, separated by comma.")

    parser.add_argument("-o", "--output", required=True,
                        help="Output file name")

    parser.add_argument("-N", "--max-qual-values", dest="max_qual", type=int, default=1000,
                        help="Maximum number of distinct quality values to calculate P/R for.")

    args = parser.parse_args()

    def checkFile(filename, label) :
        if not os.path.isfile(filename) :
            raise Exception("Can't find input %s file: '%s'" % (label,filename))

    if len(args.inputs) == 0 :
        raise Exception("No input file(s) given")

    for inputFile in args.inputs :
        checkFile(inputFile,"features CSV")

    return args


def main():
    args = parseArgs()

    datasets = []
    for i in args.inputs:
        i = os.path.abspath(i)
        print "Reading %s" % i
        df = pandas.read_csv(i)
        datasets.append(df)

    if len(datasets) > 1:
        dataset = pandas.concat(datasets)
    else:
        dataset = datasets[0]

    result = []

    fns = dataset[dataset["tag"] == "FN"].shape[0]
    # total_tp = dataset[dataset["tag"] == "TP"].shape[0]
    # total_fp = dataset[dataset["tag"] == "FP"].shape[0]
    dataset = pandas.DataFrame(dataset[dataset["tag"] != "FN"])

    for f in args.q.split(","):
        allquals = pandas.DataFrame(dataset.dropna(subset=[f]))

        qual_vals = sorted(set(allquals[f].values))
        if len(qual_vals) > args.max_qual:
            qual_vals = numpy.linspace(qual_vals[0], qual_vals[-1], args.max_qual+1)

        data_remaining = dataset

        if f == "SomaticSNVQualityAndHomRefGermlineGenotype":
            strelka_f = (dataset["NT"] == "ref") & \
                        (dataset["N_FDP_RATE"] < 0.4) & \
                        (dataset["T_FDP_RATE"] < 0.4) & \
                        (dataset["N_SDP_RATE"] < 0.75) & \
                        (dataset["T_SDP_RATE"] < 0.75) & \
                        (dataset[
                         "NormalSampleRelativeTotalLocusDepth"] < 0.3333)  # -> NormalSampleRelativeTotalLocusDepth = N_DP / chr_depth < 1 -> N_DP < chr_depth
            strelka_def_f = dataset[True != strelka_f]
            strelka_f_tps = strelka_def_f[
                strelka_def_f["tag"] == "TP"].shape[0]
            strelka_f_fps = strelka_def_f[
                strelka_def_f["tag"] == "FP"].shape[0]
            data_remaining = dataset[True == strelka_f]
        elif f == "SomaticIndelQualityAndHomRefGermlineGenotype":
            strelka_f_ref = (dataset["NT"] == "ref")
            strelka_f_ihpol = (dataset["InterruptedHomopolymerLength"] <= 14)
            strelka_f_bcnoise = (dataset["bcn"] < 0.3)
            strelka_f_repeat = (dataset["RefRepeatCount"] <= 8)
            strelka_f = strelka_f_ref & strelka_f_ihpol & strelka_f_bcnoise & strelka_f_repeat

            strelka_def_f = dataset[True != strelka_f]
            strelka_f_tps = strelka_def_f[
                strelka_def_f["tag"] == "TP"].shape[0]
            strelka_f_fps = strelka_def_f[
                strelka_def_f["tag"] == "FP"].shape[0]

            print "Strelka default filtering:"
            print "NT != ref: %i" % dataset[True != strelka_f_ref].shape[0]
            print "iHpol: %i" % dataset[True != strelka_f_ihpol].shape[0]
            print "BCNoise: %i" % dataset[True != strelka_f_bcnoise].shape[0]
            print "Repeat: %i" % dataset[True != strelka_f_repeat].shape[0]
            print "Total: %i (%i TP and %i FP)" % (dataset[True != strelka_f].shape[0], strelka_f_tps, strelka_f_fps)

            data_remaining = dataset[True == strelka_f]
        else:
            strelka_f_tps = 0
            strelka_f_fps = 0

        counter = 0
        for q in qual_vals[:-1]:
            rec = {"field": f, "qual": q}
            rec["tp"] = data_remaining[
                (data_remaining["tag"] == "TP") & (data_remaining[f] > q)].shape[0]
            rec["fp"] = data_remaining[
                (data_remaining["tag"] == "FP") & (data_remaining[f] > q)].shape[0]
            rec["tp_filtered"] = data_remaining[
                (data_remaining["tag"] == "TP") & (data_remaining[f] <= q)].shape[0]
            rec["fn"] = fns + rec["tp_filtered"] + strelka_f_tps
            rec["fp_filtered"] = data_remaining[(data_remaining["tag"] == "FP") & (
                data_remaining[f] <= q)].shape[0] + strelka_f_fps

            try:
                rec["recall"] = rec["tp"] / float(rec["tp"] + rec["fn"])
            except:
                rec["recall"] = None
            try:
                rec["precision"] = rec["tp"] / float(rec["tp"] + rec["fp"])
            except:
                rec["precision"] = None
            result.append(rec)

            counter += 1
            if counter % 10 == 0:
                print "Processed %i / %i qual values for %s" % (counter, len(qual_vals), f)

    pandas.DataFrame(result, columns=[
        "field", "qual", "tp", "fp", "fn",
        "tp_filtered", "fp_filtered", "precision", "recall"]
    ).to_csv(args.output)


if __name__ == '__main__':
    main()

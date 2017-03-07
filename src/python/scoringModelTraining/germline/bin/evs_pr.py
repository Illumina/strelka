#!/usr/bin/env python
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2017 Illumina, Inc.
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
    dataset = pandas.DataFrame(dataset[dataset["tag"] != "FN"])

    for f in args.q.split(","):
        allquals = pandas.DataFrame(dataset.dropna(subset=[f]))

        qual_vals = sorted(set(allquals[f].values))
        if len(qual_vals) > args.max_qual:
            qual_vals = random.sample(qual_vals, args.max_qual)

        data_remaining = dataset

        strelka_f_tps = 0
        strelka_f_fps = 0

        counter = 0
        for q in qual_vals:
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
            rec["na"] = data_remaining[(data_remaining["tag"] == "UNK") & (data_remaining[f] > q)].shape[0]
            rec["na_filtered"] = data_remaining[(data_remaining["tag"] == "UNK") & (data_remaining[f] <= q)].shape[0]

            try:
                rec["recall"] = rec["tp"] / float(rec["tp"] + rec["fn"])
            except:
                rec["recall"] = None
            try:
                rec["precision"] = rec["tp"] / float(rec["tp"] + rec["fp"])
            except:
                rec["precision"] = None
            try:
                rec["fracNA"] = rec["na"] / float(rec["tp"] + rec["fp"] + rec["na"])
            except:
                rec["fracNA"] = None
            result.append(rec)

            counter += 1
            if counter % 10 == 0:
                print "Processed %i / %i qual values for %s" % (counter, len(qual_vals), f)

    pandas.DataFrame(result, columns=[
        "field", "qual", "tp", "fp", "fn",
        "tp_filtered", "fp_filtered", "na", "na_filtered", "precision", "recall", "fracNA"]
    ).to_csv(args.output)


if __name__ == '__main__':
    main()

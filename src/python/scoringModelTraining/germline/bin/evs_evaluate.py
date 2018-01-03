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
Given a EVS model, evaluate a set of variant calls
"""

__author__ = "Peter Krusche <pkrusche@illumina.com>"


import os
import sys

import pandas

scriptDir = os.path.abspath(os.path.dirname(__file__))
workflowDir = os.path.abspath(
    os.path.join(scriptDir, "../lib"))

sys.path.append(workflowDir)


import evs
import evs.tools
import evs.features


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser("evs evaluation script")

    parser.add_argument("inputs", help="Feature CSV files", nargs="+")
    parser.add_argument("-c", "--classifier", required=True,
                        help="Classifier pickle file name")
    parser.add_argument("-f", "--features",
                        choices=evs.features.FeatureSet.sets.keys(),
                        required=True,
                        help="Which feature set to use (or a comma-separated list of feature names,"
                             " e.g. -f QSS_NT,T_DP_RATE")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file name")

    args = parser.parse_args()

    def checkFile(filename, label) :
        if not os.path.isfile(filename) :
            raise Exception("Can't find input %s file: '%s'" % (label,filename))

    def checkOptionalFile(filename, label) :
        if filename is None : return
        checkFile(filename,label)

    if len(args.inputs) == 0 :
        raise Exception("No input file(s) given")

    for inputFile in args.inputs :
        checkFile(inputFile,"features CSV")

    checkFile(args.classifier, "classifier")

    return args


def main():
    args = parseArgs()

    datasets = []
    for i in args.inputs:
        i = os.path.abspath(i)
        print "Reading %s" % i
        df = pandas.read_csv(i, na_values=".")
        df.fillna("0", inplace=True)
        datasets.append(df)

    if len(datasets) > 1:
        dataset = pandas.concat(datasets)
    else:
        dataset = datasets[0]

    fns = pandas.DataFrame(dataset[dataset["tag"] == "FN"])
    fns["qual"] = 1.0
    fns["ptag"] = "FN"

    dataset = pandas.DataFrame(dataset[dataset["tag"] != "FN"])

    try:
        fset = evs.features.FeatureSet.make(args.features)
        features = fset.trainingfeatures()
    except:
        features = args.features.split(",")

    model = evs.EVSModel.createFromFile(args.classifier)
    preds = model.classify(
        dataset[dataset["tag"].isin(["TP", "FP", "UNK"])], features)

    rlist = [fns, preds]

    result = pandas.concat(rlist)

    print pandas.crosstab(result['tag'], result['ptag'])

    result.to_csv(args.output)


if __name__ == '__main__':
    main()

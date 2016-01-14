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

# coding=utf-8
#
# 20/11/2014
#
# Given a EVS model, evaluate a set of variant calls
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
workflowDir = os.path.abspath(
    os.path.join(scriptDir, "@THIS_RELATIVE_PYTHON_LIBDIR@"))
templateConfigDir = os.path.abspath(
    os.path.join(scriptDir, '@THIS_RELATIVE_CONFIGDIR@'))

sys.path.append(workflowDir)

import pandas

import evs
import evs.tools
import evs.features


def main():
    parser = argparse.ArgumentParser("evs evaluation/prediction script")

    parser.add_argument("inputs", help="Feature CSV files", nargs="+")

    parser.add_argument("-c", "--classifier", dest="clf", required=True,
                        help="Classifier pickle file name")

    modelNames=evs.EVSModel.names()
    parser.add_argument("-m", "--model", dest="model", choices=modelNames, required=True,
                        help="Which model to use (options are: %s)" % str(modelNames))

    parser.add_argument("-f", "--featuresets", dest="features",
                        choices=evs.features.FeatureSet.sets.keys(),
                        required=True,
                        help="Which feature set to use (or a comma-separated list of feature names,"
                             " e.g. -f QSS_NT,T_DP_RATE")

    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file name")

    args = parser.parse_args()

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

    fns = pandas.DataFrame(dataset[dataset["tag"] == "FN"])
    fns["qual"] = 1.0
    fns["ptag"] = "FN"

    dataset = pandas.DataFrame(dataset[dataset["tag"] != "FN"])

    try:
        fset = evs.features.FeatureSet.make(args.features)
        features = fset.trainingfeatures()
    except:
        features = args.features.split(",")

    model = evs.EVSModel.create(args.model)
    model.load(args.clf)
    preds = model.classify(
        dataset[dataset["tag"].isin(["TP", "FP", "UNK"])], features)

    rlist = [fns, preds]

    result = pandas.concat(rlist)

    print pandas.crosstab(result['tag'], result['ptag'])

    result.to_csv(args.output)

if __name__ == '__main__':
    main()

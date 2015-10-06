#!/usr/bin/env python
#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

# coding=utf-8
#
# 20/11/2014
#
# Given a VQSR classifier, evaluate a set of variant calls
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

import vqsr
import vqsr.tools
import vqsr.features


def main():
    parser = argparse.ArgumentParser("vqsr evaluation/prediction script")

    parser.add_argument("inputs", help="Feature CSV files", nargs="+")

    parser.add_argument("-c", "--classifier", dest="clf", required=True,
                        help="Classifier pickle file name")

    parser.add_argument("-m", "--model", dest="model", choices=vqsr.VQSRModel.names(), required=True,
                        help="Which model to use (options are: %s)" % str(vqsr.VQSRModel.names()))

    parser.add_argument("-f", "--featuresets", dest="features",
                        choices=vqsr.features.FeatureSet.sets.keys(),
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
        fset = vqsr.features.FeatureSet.make(args.features)
        features = fset.trainingfeatures()
    except:
        features = args.features.split(",")

    model = vqsr.VQSRModel.create(args.model)
    model.load(args.clf)
    preds = model.classify(
        dataset[dataset["tag"].isin(["TP", "FP", "UNK"])], features)

    rlist = [fns, preds]

    result = pandas.concat(rlist)

    print pandas.crosstab(result['tag'], result['ptag'])

    result.to_csv(args.output)

if __name__ == '__main__':
    main()

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
# Learn VQSR Classifiers from (set of) feature files
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
import json

scriptDir = os.path.abspath(os.path.dirname(__file__))
scriptName = os.path.basename(__file__)
workflowDir = os.path.abspath(
    os.path.join(scriptDir, "@THIS_RELATIVE_PYTHON_LIBDIR@"))
templateConfigDir = os.path.abspath(
    os.path.join(scriptDir, '@THIS_RELATIVE_CONFIGDIR@'))

sys.path.append(workflowDir)

import pandas
import random

import vqsr
import vqsr.tools
import vqsr.features


def main():
    parser = argparse.ArgumentParser("vqsr learning script")

    parser.add_argument("inputs", help="Feature CSV files", nargs="+")

    parser.add_argument("-m", "--model", dest="model", choices=vqsr.VQSRModel.names(), required=True,
                        help="Which model to use (options are: %s)" % str(vqsr.VQSRModel.names()))

    parser.add_argument("-f", "--featuresets", dest="features",
                        choices=vqsr.features.FeatureSet.sets.keys(),
                        required=True,
                        help="Which feature set to use (or a comma-separated list of feature names,"
                             " e.g. -f QSS_NT,T_DP_RATE")

    parser.add_argument("-p", "--parameter-file", dest="parameters", default=None,
                        help="Specify additional parameters for the learning algorithm (in a JSON file)")

    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file name")

    parser.add_argument("--balance", dest="balance", default=False, action="store_true",
                        help="Balance the number of samples from each label.")

    parser.add_argument("--sample-input", dest="sample_input", default=0, type=int,
                        help="Number of rows to subsample from each input data file")

    parser.add_argument("--plots", default=False, action="store_true",
                        help="Make plots.")

    args = parser.parse_args()

    datasets = []
    for i in args.inputs:
        i = os.path.abspath(i)
        print "Reading %s" % i
        df = pandas.read_csv(i)

        if args.sample_input:
            p_rows = min(df.shape[0], args.sample_input)
            p_rows_selected = random.sample(df.index, p_rows)
            df = pandas.DataFrame(df.ix[p_rows_selected])

        datasets.append(df)

    if len(datasets) > 1:
        dataset = pandas.concat(datasets)
    else:
        dataset = datasets[0]

    dataset = pandas.DataFrame(dataset[dataset["tag"] != "FN"])

    try:
        fset = vqsr.features.FeatureSet.make(args.features)
        features = fset.trainingfeatures()
    except:
        features = args.features.split(",")

    try:
        pars = json.load(open(args.parameters))
        print "Using learning parameters: %s" % str(pars)
    except:
        print "Using default parameters."
        pars = {}

    model = vqsr.VQSRModel.create(args.model)

    if not args.balance:
        tpdata = dataset[dataset["tag"] == "TP"]
        fpdata = dataset[dataset["tag"] == "FP"]
    else:
        tpdata2 = dataset[dataset["tag"] == "TP"]
        fpdata2 = dataset[dataset["tag"] == "FP"]
        nrows = min(tpdata2.shape[0], fpdata2.shape[0])
        tp_rows_selected = random.sample(tpdata2.index, nrows)
        fp_rows_selected = random.sample(fpdata2.index, nrows)
        tpdata = tpdata2.ix[tp_rows_selected]
        fpdata = fpdata2.ix[fp_rows_selected]

    model.train(tpdata, fpdata, features, **pars)
    model.save(args.output)

    if args.plots and hasattr(model, "plots"):
        model.plots(args.output + ".plots", features)
    elif args.plots:
        print "No plots created, this is not supported by %s" % args.model


if __name__ == '__main__':
    main()

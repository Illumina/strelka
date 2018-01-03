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
Given an EVS model evaluation, provide binned scores and empirical precision
"""

__author__ = "Konrad Scheffler <kscheffler@illumina.com>"

import os
import sys
import random
import numpy
import pandas
import json
from sklearn import linear_model


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser("evs QQ script")
    parser.add_argument("inputs", help="Quality/tag CSV files", nargs="+")
    parser.add_argument("-o", "--output", required=True,
                        help="QQ output file name (csv)")
    parser.add_argument("-c", "--calibration", required=True,
                        help="Calibration output file name (json)")
    parser.add_argument("-N", "--bins", dest="nbins", type=int, default=20,
                        help="Number of bins for which to calculate empirical precision.")
    parser.add_argument("-G", "--geno", dest="geno", type=int, default="-1",
                        help="If specified, restrict analysis to variants with this GENO value.")

    args = parser.parse_args()

    def checkFile(filename, label) :
        if not os.path.isfile(filename) :
            raise Exception("Can't find input %s file: '%s'" % (label,filename))

    if len(args.inputs) == 0 :
        raise Exception("No input file(s) given")

    for inputFile in args.inputs :
        checkFile(inputFile,"features CSV")

    return args


def phred(corrprob):
    return numpy.minimum(-10*numpy.log10(1-corrprob),numpy.full(corrprob.shape,50))


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

    if args.geno >= 0:
        dataset = dataset[dataset["GENO"] == args.geno]

    result = []
    regression = []

    print "Total: %d" % len(dataset)
    fns = dataset[dataset["tag"] == "FN"].shape[0]
    print "fn: %d" % fns
    print "tp: %d" % len(dataset[dataset["tag"] == "TP"])
    print "fp: %d" % len(dataset[dataset["tag"] == "FP"])
    print "unk: %d" % len(dataset[dataset["tag"] == "UNK"])
    dataset = pandas.DataFrame(dataset[dataset["tag"].isin(["TP", "FP"])])
    data = pandas.DataFrame(dataset.dropna(subset=["qual"]))
    data = data.iloc[numpy.random.permutation(len(data))] # shuffle before sorting: this ensures random order within tied groups
    data = data.sort("qual")

    # note this call fails on numpy 1.9.x (default package for Ubuntu 14.04) TODO: workaround?
    bins = numpy.array_split(data,args.nbins)

    counter = 0
    for bin in bins:
        rec = {"qual": numpy.median(bin["qual"])}
        rec["tp"] = bin[bin["tag"] == "TP"].shape[0]
        rec["fp"] = bin[bin["tag"] == "FP"].shape[0]
        rec["na"] = bin[bin["tag"] == "UNK"].shape[0]

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
            print "Processed %i / %i values" % (counter, len(bins))

    df = pandas.DataFrame(result)
    df["Q"] = phred(df["qual"])
    df["empiricalQ"] = phred(df["precision"])
    regr = linear_model.LinearRegression()
    regr.fit(phred(df["qual"]).reshape(-1,1), phred(df["precision"]).reshape(-1,1))
    df["calibratedQ"] = regr.predict(phred(df["qual"]).reshape(-1,1))
    print df
    print("Coefficient: {}, Intercept: {}".format(regr.coef_[0,0], regr.intercept_[0]))
    df.to_csv(args.output)
    json.dump({"Coefficient" : regr.coef_[0,0], "Intercept" : regr.intercept_[0]}, open(args.calibration, 'w'))


if __name__ == '__main__':
    main()

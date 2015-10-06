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
# Randomly sample subsets of TPs and FPs that have a given size.
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
import pandas
import random


def main():
    parser = argparse.ArgumentParser("Split a set of CSV files")

    parser.add_argument("inputs", help="Feature CSV files", nargs="+")

    parser.add_argument("-n", "--sample-size", dest="samplesize", help="Number of TPs/FPs to sample", type=int, default=1000)

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

    tpdata2 = dataset[dataset["tag"] == "TP"]
    fpdata2 = dataset[dataset["tag"] == "FP"]

    nrows = min([tpdata2.shape[0], fpdata2.shape[0]])

    if nrows < args.samplesize:
        print >> sys.stderr, "Warning: sampling will use some records multiple times because there is not enough data!"

    tp_rows_selected = random.sample(tpdata2.index, args.samplesize)
    fp_rows_selected = random.sample(fpdata2.index, args.samplesize)
    tpdata = tpdata2.ix[tp_rows_selected]
    fpdata = fpdata2.ix[fp_rows_selected]

    dataset = pandas.concat([tpdata, fpdata])
    dataset.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()

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

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
# Split a (set of) CSV files into two disjoint random subsets (cross-validation).
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
import argparse
import pandas
import random


def main():
    parser = argparse.ArgumentParser("Split a set of CSV files")

    parser.add_argument("inputs", help="Feature CSV files", nargs="+")

    parser.add_argument("-r", "--ratio", help="p/q ratio", type=float, default=0.5)

    parser.add_argument("-p", "--output-1", dest="output_1", required=True,
                        help="Output file name 1")

    parser.add_argument("-q", "--output-2", dest="output_2", required=False,
                        help="Output file name 2")

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

    p_rows = int(args.ratio*dataset.shape[0])

    p_rows_selected = random.sample(dataset.index, p_rows)

    try:
        dataset.ix[p_rows_selected].sort(["CHROM", "POS"]).to_csv(args.output_1, index=False)
    except:
        dataset.ix[p_rows_selected].to_csv(args.output_1, index=False)

    if args.output_2:
        dataset.drop(p_rows_selected).to_csv(args.output_2, index=False)


if __name__ == '__main__':
    main()

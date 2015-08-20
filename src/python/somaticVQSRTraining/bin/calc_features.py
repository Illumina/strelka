#!/illumina/development/somaticVQSR/python-virtualenv/bin/python
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
# Calculate additional features and write updated CSV
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
libDir = os.path.abspath(
    os.path.join(scriptDir, "@THIS_RELATIVE_PYTHON_LIBDIR@"))
sys.path.append(libDir)


import pandas

from vqsr.features.ref import getReference
from vqsr.features.entropy import getEntropy
from vqsr.features.repeat import getRepeats

def main():
    parser = argparse.ArgumentParser("Split a set of CSV files")

    parser.add_argument("inputs", help="Feature CSV files", nargs="+")

    parser.add_argument("--ref-fasta", help="Fasta reference sequence file name", dest="ref",
                        default="/illumina/development/iSAAC/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")

    parser.add_argument("--window", dest="window", help="Window length for entropy calculation",
                        type=int, default=200)

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

    dataset = getReference(dataset, args.ref, args.window)
    dataset = getEntropy(dataset)
    dataset = getRepeats(dataset)

    dataset.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()

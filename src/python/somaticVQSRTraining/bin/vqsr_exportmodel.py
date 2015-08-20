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
# Export VQSR model to JSON format
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

import vqsr
import vqsr.tools


def main():
    parser = argparse.ArgumentParser("vqsr learning script")

    parser.add_argument("-c", "--classifier", dest="clf", required=True,
                        help="Classifier pickle file name")

    parser.add_argument("-m", "--model", dest="model", choices=vqsr.VQSRModel.names(), required=True,
                        help="Which model to use (options are: %s)" % str(vqsr.VQSRModel.names()))

    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file name")

    args = parser.parse_args()

    if not args.output.endswith(".json"):
        args.output += ".json"

    model = vqsr.VQSRModel.create(args.model)
    model.load(args.clf)
    model.save(args.output)


if __name__ == '__main__':
    main()

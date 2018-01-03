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
Export EVS model to JSON format
"""

__author__ = "Peter Krusche <pkrusche@illumina.com>"


import os
import sys
import pandas
import json

scriptDir = os.path.abspath(os.path.dirname(__file__))
workflowDir = os.path.abspath(
#    os.path.join(scriptDir, "@THIS_RELATIVE_PYTHON_LIBDIR@"))
    os.path.join(scriptDir, "../lib"))

sys.path.append(workflowDir)

import evs


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser("evs learning script")

    parser.add_argument("-c", "--classifier", required=True,
                        help="Classifier pickle file name")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file name")
    parser.add_argument("--calltype", default="Germline",
                        help="Call type (Germline/RNAseq)")
    parser.add_argument("-v", "--varianttype", default="SNV",
                        help="Variant type (SNV/INDEL)")
    parser.add_argument("-C", "--calibration",
                        help="Calibration file name (if omitted, model is assumed to be already calibrated)")
    parser.add_argument("-t", "--threshold", type = int, default = 3,
                        help="Q-score threshold")

    args = parser.parse_args()

    def checkFile(filename, label) :
        if not os.path.isfile(filename) :
            raise Exception("Can't find input %s file: '%s'" % (label,filename))

    checkFile(args.classifier, "classifier")
    if args.calibration != None :
        checkFile(args.calibration, "calibration")

    return args


def main() :
    args = parseArgs()
    if args.calibration == None :
        caldict = {"Coefficient" : 1, "Intercept" : 1}
    else:
        caldict = json.load(open(args.calibration))

    # Convert to linear domain for Strelka:
    power = caldict["Coefficient"]
    scale = 10**(-0.1*caldict["Intercept"])

    model = evs.EVSModel.createFromFile(args.classifier)
    model.save_json_strelka_format(args.output, args.calltype, args.varianttype, power, scale, args.threshold)

if __name__ == '__main__':
    main()

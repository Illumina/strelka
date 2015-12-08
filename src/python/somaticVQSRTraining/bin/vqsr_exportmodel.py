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

    model = vqsr.VQSRModel.create(args.model)
    model.load(args.clf)
    
    model.save_json_strelka_format(args.output)


if __name__ == '__main__':
    main()

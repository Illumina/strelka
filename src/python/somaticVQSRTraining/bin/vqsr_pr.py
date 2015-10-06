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
import random


def main():
    parser = argparse.ArgumentParser("vqsr precision/recall script")
    parser.add_argument("inputs", help="Quality/tag CSV files", nargs="+")

    parser.add_argument("-q", "--quality-fields", dest="q", default="QSS_NT,qual",
                        help="Select fields which give quality values, separated by comma. "
                             "If QSS_NT is selected, standard Strelka filtering will be output.")

    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file name")

    parser.add_argument("-N", "--max-qual-values", dest="max_qual", type=int, default=1000,
                        help="Maximum number of distinct quality values to calculate P/R for.")

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

    result = []

    fns = dataset[dataset["tag"] == "FN"].shape[0]
    # total_tp = dataset[dataset["tag"] == "TP"].shape[0]
    # total_fp = dataset[dataset["tag"] == "FP"].shape[0]
    dataset = pandas.DataFrame(dataset[dataset["tag"] != "FN"])

    for f in args.q.split(","):
        allquals = pandas.DataFrame(dataset.dropna(subset=[f]))

        qual_vals = sorted(set(allquals[f].values))
        if len(qual_vals) > args.max_qual:
            qual_vals = random.sample(qual_vals, args.max_qual)

        data_remaining = dataset

        if f == "QSS_NT":
            strelka_f = (dataset["NT"] == "ref") & \
                        (dataset["N_FDP_RATE"] < 0.4) & \
                        (dataset["T_FDP_RATE"] < 0.4) & \
                        (dataset["N_SDP_RATE"] < 0.75) & \
                        (dataset["T_SDP_RATE"] < 0.75) & \
                        (dataset[
                         "N_DP_RATE"] < 1.0)  # -> N_DP_RATE = N_DP / chr_depth < 1 -> N_DP < chr_depth
            strelka_def_f = dataset[True != strelka_f]
            strelka_f_tps = strelka_def_f[
                strelka_def_f["tag"] == "TP"].shape[0]
            strelka_f_fps = strelka_def_f[
                strelka_def_f["tag"] == "FP"].shape[0]
            data_remaining = dataset[True == strelka_f]
        elif f == "QSI_NT":
            strelka_f_ref = (dataset["NT"] == "ref")
            strelka_f_ihpol = (dataset["IHP"] <= 14)
            strelka_f_bcnoise = (dataset["bcn"] < 0.3)
            strelka_f_repeat = (dataset["RC"] <= 8)
            strelka_f = strelka_f_ref & strelka_f_ihpol & strelka_f_bcnoise & strelka_f_repeat

            strelka_def_f = dataset[True != strelka_f]
            strelka_f_tps = strelka_def_f[
                strelka_def_f["tag"] == "TP"].shape[0]
            strelka_f_fps = strelka_def_f[
                strelka_def_f["tag"] == "FP"].shape[0]

            print "Strelka default filtering:"
            print "NT != ref: %i" % dataset[True != strelka_f_ref].shape[0]
            print "iHpol: %i" % dataset[True != strelka_f_ihpol].shape[0]
            print "BCNoise: %i" % dataset[True != strelka_f_bcnoise].shape[0]
            print "Repeat: %i" % dataset[True != strelka_f_repeat].shape[0]
            print "Total: %i (%i TP and %i FP)" % (dataset[True != strelka_f].shape[0], strelka_f_tps, strelka_f_fps)

            data_remaining = dataset[True == strelka_f]
        else:
            strelka_f_tps = 0
            strelka_f_fps = 0

        counter = 0
        for q in qual_vals:
            rec = {"field": f, "qual": q}
            rec["tp"] = data_remaining[
                (data_remaining["tag"] == "TP") & (data_remaining[f] > q)].shape[0]
            rec["fp"] = data_remaining[
                (data_remaining["tag"] == "FP") & (data_remaining[f] > q)].shape[0]
            rec["tp_filtered"] = data_remaining[
                (data_remaining["tag"] == "TP") & (data_remaining[f] <= q)].shape[0]
            rec["fn"] = fns + rec["tp_filtered"] + strelka_f_tps
            rec["fp_filtered"] = data_remaining[(data_remaining["tag"] == "FP") & (
                data_remaining[f] <= q)].shape[0] + strelka_f_fps

            try:
                rec["recall"] = rec["tp"] / float(rec["tp"] + rec["fn"])
            except:
                rec["recall"] = None
            try:
                rec["precision"] = rec["tp"] / float(rec["tp"] + rec["fp"])
            except:
                rec["precision"] = None
            result.append(rec)

            counter += 1
            if counter % 10 == 0:
                print "Processed %i / %i qual values for %s" % (counter, len(qual_vals), f)

    pandas.DataFrame(result, columns=[
        "field", "qual", "tp", "fp", "fn",
        "tp_filtered", "fp_filtered", "precision", "recall"]
    ).to_csv(args.output)


if __name__ == '__main__':
    main()

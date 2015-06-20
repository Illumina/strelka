#!/illumina/development/haplocompare/hc-virtualenv/bin/python
# coding=utf-8
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/sequencing/licenses/blob/master/Simplified-BSD-License.txt

import pandas

from vqsr.tools.vcfextract import vcfExtract
from . import FeatureSet


@FeatureSet.register("posandalleles")
class PosAndAlleles(FeatureSet):
    """ Collect position and alleles from VCF """

    def collect(self, vcfname):
        """ Return a data frame with features collected from the
            given VCF, tagged by given type """
        features = ["CHROM", "POS", "REF", "ALT"]
        records = []
        for vr in vcfExtract(vcfname, features):
            rec = {}
            for i, ff in enumerate(features):
                rec[ff] = vr[i]

            # Gather the computed data into a dict
            qrec = {
                "CHROM": rec["CHROM"],
                "POS": int(rec["POS"]),
                "REF": rec["REF"],
                "ALT": ",".join(rec["ALT"]),
            }
            records.append(qrec)

        cols = ["CHROM", "POS", "REF", "ALT"]

        if records:
            df = pandas.DataFrame(records, columns=cols)
        else:
            df = pandas.DataFrame(columns=cols)

        return df

    def trainingfeatures(self):
        """ Return a list of columns that are features to use for VQSR training """
        return []

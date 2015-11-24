#!/illumina/development/haplocompare/hc-virtualenv/bin/python
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

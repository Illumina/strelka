#!/illumina/development/haplocompare/hc-virtualenv/bin/python
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

import pandas

from evs.tools.vcfextract import vcfExtract
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
        """ Return a list of columns that are features to use for EVS model training """
        return []

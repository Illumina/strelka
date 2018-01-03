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

import pandas

from evs.tools.vcf import openMaybeGzip, VCFID
from . import FeatureSet


class VcfFeatureSet(FeatureSet):

    def collectCore(self, vcfname, headerKey = None):
        """
        Return a data frame with features collected from the given VCF

        If headerKey is provided, then use this header value to extract labels
        for the INFO EVSF feature tag
        """
        feature_labels = ["CHROM", "POS", "REF", "ALT"]
        header_feature_labels = None

        records = []

        isHeader = True
        isHeaderKey = (headerKey is not None)

        for line in openMaybeGzip(vcfname):
            if isHeader :
                if line[0] == "#" :
                    if isHeaderKey and line.startswith("##") :
                        word = line[2:].strip().split("=")
                        if word[0] == headerKey :
                            assert(header_feature_labels is None)
                            header_feature_labels = word[1].split(",")
#                            print header_feature_labels
                            assert(len(header_feature_labels) > 0)
                    continue
                else :
                    if isHeaderKey :
                        assert(header_feature_labels is not None)
                    isHeader = False

            word = line.strip().split('\t')

            qrec = {
                "CHROM": word[VCFID.CHROM],
                "POS": int(word[VCFID.POS]),
                "REF": word[VCFID.REF],
                "ALT": word[VCFID.ALT]
            }

            if isHeaderKey :
                for ikv in word[VCFID.INFO].split(';') :
                    iword = ikv.split("=",1)
                    if iword[0] != "EVSF" : continue
                    assert(len(iword) == 2)
                    if len(word[VCFID.ALT]) > 1: continue # skip indels
                    features = [float(f) for f in iword[1].split(',')]
#                    print features
#                    assert(len(features) == len(header_feature_labels))
                    for i in range(len(features)) :
                        qrec[header_feature_labels[i]] = features[i]

            records.append(qrec)

        cols = feature_labels
        if isHeaderKey :
            cols += header_feature_labels

        return pandas.DataFrame(records, columns=cols)

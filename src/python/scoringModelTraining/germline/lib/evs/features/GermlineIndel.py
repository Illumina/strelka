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

from VcfFeatureSet import VcfFeatureSet


@VcfFeatureSet.register("germline.indel")
class GermlineSNVFeatures(VcfFeatureSet):
    """Collect germline INDEL features from VCF"""

    def collect(self, vcfname):
        """ Return a data frame with features collected from the
            given VCF"""

        return self.collectCore(vcfname,"germline_indel_scoring_features")


    def trainingfeatures(self):
        """ Return a list of columns that are features to
            use for EVS training

            Any change here must be done together with changing
            src/c++/lib/applications/strelka/strelkaScoringFeatures.hh
        """
        return [
            "GENO",
#            "REFREP1",
            "IDREP1",
            "RULEN1",
            "ABlower",
            "AB",
            "F_MQ",
#            "mapqZeroFraction",
#            "F_DPI_NORM",
            "TDP_NORM",
#            "QUAL_NORM",
#            "F_GQX_NORM",
#            "F_GQ_NORM",
#            "AD0_NORM",
            "AD1_NORM",
#            "AD2_NORM",
#            "QUAL_EXACT",
            "F_GQX_EXACT"]
#            "F_GQ_EXACT"]



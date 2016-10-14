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

"""
parses a hap.py-annotated Strelka VCF with EVS features into csv format used for EVS model training
"""

import os.path, sys


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < annotated.vcf"
    parser = OptionParser(usage=usage)

    parser.add_option("--scoringFeatures", type="string", dest="scoringFeaturesPath",metavar="FILE",
                      help="Scoring feature lists extracted from the EVS feature variant VCF")
    parser.add_option("--snvOutput", type="string", dest="snvOutputPath",metavar="FILE",
                      help="Write labeled SNV feature output in csv format to this file (required)")
    parser.add_option("--indelOutput", type="string", dest="indelOutputPath",metavar="FILE",
                      help="Write labeled indel feature output in csv format to this file (required)")

    (options,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    if options.scoringFeaturesPath is None :
        parser.error("Scoring features filename is required")
    if not os.path.isfile(options.scoringFeaturesPath) :
        parser.error("Can't find scoring features file: '%s'" % (options.scoringFeaturesPath))
    if options.snvOutputPath is None :
        parser.error("SNV output filename is required")
    if options.indelOutputPath is None :
        parser.error("Indel output filename is required")

    return (options,args)


class VCFID :
    CHROM = 0
    POS = 1
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7
    FORMAT = 8
    SAMPLE = 9


def getKeyVal(string,key) :
    import re
    match=re.search("%s=([^;\t]*);?" % (key) ,string)
    if match is None : return None
    return match.group(1)


def main() :

    import re

    (options,args) = getOptions()

    infp = sys.stdin

    snv_outfp = open(options.snvOutputPath,"w")
    indel_outfp = open(options.indelOutputPath,"w")

    class HeaderData :
        isOutputHeaderInitialized = False
        snvFeatures = None
        indelFeatures = None

    def processHeaderLine(line) :
        if line.find("snv_scoring_features=") != -1 :
            assert(HeaderData.snvFeatures is None)
            tword = line.strip().split("snv_scoring_features=")
            assert(len(tword) == 2)
            HeaderData.snvFeatures = tword[1]
        elif line.find("indel_scoring_features=") != -1 :
            assert(HeaderData.indelFeatures is None)
            tword = line.strip().split("indel_scoring_features=")
            assert(len(tword) == 2)
            HeaderData.indelFeatures = tword[1]

    def finalizeHeader() :
        HeaderData.isOutputHeaderInitialized = True

        def writeCsvHeader(fp, features, label) :
            if features is None :
                raise Exception("Can't find '%s_scoring_features' in input VCF header" % (label))
            fp.write("CHROM,POS,TYPE,%s,tag\n" % (features))

        writeCsvHeader(snv_outfp, HeaderData.snvFeatures, "snv")
        writeCsvHeader(indel_outfp, HeaderData.indelFeatures, "indel")

    for line in open(options.scoringFeaturesPath) :
        if line[0] == "#" :
            processHeaderLine(line)

    finalizeHeader()

    for line in infp :
        if line[0] == "#" :
            continue

        word = line.strip().split('\t')

        # expecting standard happy annotation with truth/query samples:
        assert(len(word) == (VCFID.SAMPLE+2))

        filterVals = word[VCFID.FILTER].split(';')

        # Skip entries matching OffTarget in the filter field (for WES data)
        if "OffTarget" in filterVals : continue

        evsf = getKeyVal(word[VCFID.INFO],"EVSF")
        if evsf is None : continue

        isSNV = (word[VCFID.INFO].find("CIGAR=") == -1)

        formatVals = word[VCFID.FORMAT].split(":")
        queryVals = word[VCFID.SAMPLE+1].split(":")

        if "NOCALL" in queryVals : continue

        sampleBDIndex = formatVals.index("BD")
        qlabel = queryVals[sampleBDIndex]
        if qlabel not in ("TP","FP","UNK") :
            raise Exception("Query value is not TP|FP|UNK as expected:\n%s" % (line))

        def typeLabel(isSNV) :
            if isSNV : return "SNP"
            else : return "INDEL"

        sampleBVTIndex = formatVals.index("BVT")
        qtype = queryVals[sampleBVTIndex]
        if qtype != typeLabel(isSNV) : continue

        def typeStream(isSNV) :
            if isSNV : return snv_outfp
            else :     return indel_outfp

        typeStream(isSNV).write(",".join([word[VCFID.CHROM], word[VCFID.POS], qtype, evsf, qlabel]) +"\n")


if __name__ == "__main__" :
    main()

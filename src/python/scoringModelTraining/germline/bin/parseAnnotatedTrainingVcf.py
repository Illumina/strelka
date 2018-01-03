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
parses a hap.py-annotated Strelka VCF with EVS features into csv format used for EVS model training
"""

import os.path, sys


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < annotated.vcf"
    parser = OptionParser(usage=usage)

    parser.add_option("--snvOutput", type="string", dest="snvOutputPath", metavar="FILE",
                      help="Write labeled SNV feature output for training data (default: all of it) in csv format to this file (required)")
    parser.add_option("--indelOutput", type="string", dest="indelOutputPath", metavar="FILE",
                      help="Write labeled indel feature output for training data (default: all of it) in csv format to this file (required)")
    parser.add_option("--testSet", type="string", action='append',
                      help="Chromosome (e.g. chr20) to hold out as test data (may be specified more than once; if omitted, all data will be used for training)")
    parser.add_option("--snvTestOutput", type="string", dest="snvTestOutputPath",metavar="FILE",
                      help="Write labeled SNV feature output for test data in csv format to this file (optional)")
    parser.add_option("--indelTestOutput", type="string", dest="indelTestOutputPath",metavar="FILE",
                      help="Write labeled indel feature output for test data in csv format to this file (optional)")
    parser.add_option("--suppressGTMismatch", dest="suppressGTMismatch", default=False, action="store_true",
                      help="When variant sequence matches but GT does not, treat as a true positive (recommended for RNA EVS)")
    parser.add_option("--discardFNs", dest="discardFNs", default=False, action="store_true",
                      help="Do not output FN variants (recommended for RNA EVS)")
    parser.add_option("--removeRNAEditing", dest="removeRNAEditing", default=False, action="store_true",
                      help="Label potential RNA editing sites (variants with A->G or T->C) as unknown (recommended for RNA EVS)")

    (options,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    if options.snvOutputPath is None :
        parser.error("SNV output filename is required")
    if options.indelOutputPath is None :
        parser.error("Indel output filename is required")
    if options.testSet and (options.snvTestOutputPath is None) :
        parser.error("SNV test output filename is required when specifying a test set")
    if options.testSet and (options.indelTestOutputPath is None) :
        parser.error("Indel test output filename is required when specifying a test set")

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
    if options.testSet :
        snv_test_outfp = open(options.snvTestOutputPath,"w")
        indel_test_outfp = open(options.indelTestOutputPath,"w")

    class HeaderData :
        isOutputHeaderInitialized = False
        snvFeatures = None
        indelFeatures = None
        snvEmptyFeatureString = None
        indelEmptyFeatureString = None

    def processHeaderLine(line) :
        if line.find("snv_scoring_features=") != -1 :
            assert(HeaderData.snvFeatures is None)
            tword = line.strip().split("snv_scoring_features=")
            assert(len(tword) == 2)
            HeaderData.snvFeatures = tword[1]
            HeaderData.snvEmptyFeatureString = ',' * HeaderData.snvFeatures.count(',')
        elif line.find("indel_scoring_features=") != -1 :
            assert(HeaderData.indelFeatures is None)
            tword = line.strip().split("indel_scoring_features=")
            assert(len(tword) == 2)
            HeaderData.indelFeatures = tword[1]
            HeaderData.indelEmptyFeatureString = ',' * HeaderData.indelFeatures.count(',')

    def finalizeHeader() :
        HeaderData.isOutputHeaderInitialized = True

        def writeCsvHeader(fp, features, label) :
            if features is None :
                raise Exception("Can't find '%s_scoring_features' in input VCF header" % (label))
            fp.write("CHROM,POS,TYPE,%s,tag\n" % (features))

        writeCsvHeader(snv_outfp, HeaderData.snvFeatures, "snv")
        writeCsvHeader(indel_outfp, HeaderData.indelFeatures, "indel")
        if options.testSet :
            writeCsvHeader(snv_test_outfp, HeaderData.snvFeatures, "snv")
            writeCsvHeader(indel_test_outfp, HeaderData.indelFeatures, "indel")


    for line in infp :
        if line[0] == "#" :
            processHeaderLine(line)
        if (HeaderData.snvFeatures is not None) and (HeaderData.indelFeatures is not None) :
            finalizeHeader()
            break

    for line in infp :
        if line[0] == "#" :
            continue

        word = line.strip().split('\t')

        # expecting standard happy annotation with truth/query samples:
        assert(len(word) == (VCFID.SAMPLE+2))

        isTrain = (not options.testSet) or (word[VCFID.CHROM] not in options.testSet)

        filterVals = word[VCFID.FILTER].split(';')

        # Skip entries matching OffTarget in the filter field (for legacy WES data)
        if "OffTarget" in filterVals : continue

        def typeLabel(isSNV) :
            if isSNV : return "SNP"
            else : return "INDEL"

        formatVals = word[VCFID.FORMAT].split(":")
        truthVals = word[VCFID.SAMPLE].split(":")
        queryVals = word[VCFID.SAMPLE+1].split(":")
        sampleBDIndex = formatVals.index("BD")
        sampleBVTIndex = formatVals.index("BVT")
        qtype = queryVals[sampleBVTIndex]
        evsf = getKeyVal(word[VCFID.INFO],"EVSF")
        if evsf is None :
            # hap.py adds FN entries with no EVSF. Let hap.py decide from the truth variant whether these are indel or snv:
            isSNV = (truthVals[sampleBVTIndex]==typeLabel(True))
            # populate the empty evsf field with dummy information of the right length:
            if isSNV :
                evsf = HeaderData.snvEmptyFeatureString
            else :
                evsf = HeaderData.indelEmptyFeatureString
        else:
            # For entries with EVSF, use presence/absence of CIGAR to decide if indel or snv. Discard variant if this is not consistent with the type assigned by hap.py:
            isSNV = (word[VCFID.INFO].find("CIGAR=") == -1)
            if (qtype != typeLabel(isSNV)) : continue

        if "NOCALL" in queryVals :
            if truthVals[sampleBDIndex] != "FN" : continue
            label = "FN"
        else :
            label = queryVals[sampleBDIndex]

        if options.suppressGTMismatch and label == "FP":
            if getKeyVal(word[VCFID.INFO],"kind") == "gtmismatch":
                label = "TP"

        if options.removeRNAEditing and isSNV:
            if ((word[VCFID.REF] == "A" and word[VCFID.ALT] == "G") or
                (word[VCFID.REF] == "T" and word[VCFID.ALT] == "C")):
                label = "UNK"

        if label not in ("TP","FP","FN","UNK") :
            raise Exception("Variant label is not TP|FP|FN|UNK as expected:\n%s" % (line))

        if qtype is None :
            raise Exception("No valid type in input line:\n%s" % (line))

        def outputStream(isSNV, isTrain) :
            if isTrain :
                if isSNV : return snv_outfp
                else :     return indel_outfp
            else :
                if isSNV : return snv_test_outfp
                else :     return indel_test_outfp

        def writeVariant(outputlabel) :
            if not(options.discardFNs and outputlabel == "FN"):
                outputStream(isSNV, isTrain).write(",".join([word[VCFID.CHROM], word[VCFID.POS], qtype, evsf, outputlabel]) +"\n")

        writeVariant(label)
        # Add an extra FN entry if a variant is truth FN,query FP:
        if (label == "FP" and truthVals[sampleBDIndex] == "FN") :
            writeVariant("FN")


if __name__ == "__main__" :
    main()

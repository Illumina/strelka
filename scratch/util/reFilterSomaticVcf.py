#!/usr/bin/env python
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2017 Illumina, Inc.
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

import sys
import re


def getKeyVal(string,key) :
    match=re.search("%s=([^;\t]*);?" % (key) ,string)
    if match is None : return None
    return match.group(1);


class VCFID :
    CHROM = 0
    POS = 1
    REF = 3
    ALT = 4
    QUAL = 5
    FILTER = 6
    INFO = 7


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [options] < in.vcf > out.vcf"
    description="""
Reset LowEVS filter in Strelka2 somatic VCFs to reflect a new filtration theshold. 
"""
    parser = OptionParser(usage=usage)

    parser.add_option("--minSomaticEVS", dest="minSomaticEVS", type="float",
                      help="minimum somatic EVS score, no default (required)")

    (opt,args) = parser.parse_args()

    if (opt.minSomaticEVS is None) or (len(args) != 0) :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if opt.minSomaticEVS < 0:
        raise Exception("Invalid minSomaticEVS value: %f" % (opt.minSomaticEVS))

    return (opt,args)


def main() :

    (opt,args) = getOptions()

    filterTag = "LowEVS"

    infp = sys.stdin
    outfp = sys.stdout

    for line in infp :
        if line[0] == "#" :
            outfp.write(line)
            continue

        w=line.strip().split('\t')
        assert(len(w) > VCFID.INFO)

        
        val = getKeyVal(w[VCFID.INFO],"SomaticEVS")
        if val is None :    
            outfp.write(line)
            continue

        isPass = (float(val) >= opt.minSomaticEVS)

        filters = set(w[VCFID.FILTER].split(';'))

        if len(filters) == 1 and ("PASS" in filters or "." in filters) : filters = set() 

        if isPass :
            if filterTag in filters : filters.remove(filterTag)
        else      : filters.add(filterTag)
        
        if len(filters) == 0 :
            w[VCFID.FILTER] = "PASS"
        else :
            w[VCFID.FILTER] = ';'.join(filters) 

        sys.stdout.write('\t'.join(w) + '\n')


main()

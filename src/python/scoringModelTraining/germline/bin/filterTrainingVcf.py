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
Filters out some gvcf entries before passing it to hap.py.
"""

import sys


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog < input.vcf > filtered.vcf"
    parser = OptionParser(usage=usage)

    (options,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

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


def main() :

    import re

    (options,args) = getOptions()

    infp = sys.stdin
    outfp = sys.stdout

    for line in infp :
        if line[0] == "#" :
            outfp.write(line)
            continue

        word = line.strip().split('\t')

        # EVS training procedure requires exactly one sample
        assert(len(word) == (VCFID.SAMPLE+1))

        filterVals = word[VCFID.FILTER].split(';')

        # Skip entries matching OffTarget in the filter field (for legacy WES data)
        if "OffTarget" in filterVals : continue

        # Skip entries matching any type of conflict
        foundConflict = False
        for filterVal in filterVals :
            if filterVal.endswith("Conflict") :
                foundConflict = True
                continue
        if foundConflict: continue

        formatVals = word[VCFID.FORMAT].split(':')
        sampleVals = word[VCFID.SAMPLE].split(':')

        outfp.write(line)

if __name__ == "__main__" :
    main()

#!/usr/bin/env bash
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

#
# facilitate comparison of strelka vcf output from two runs by stripping out the header and uninteresting field differences
#

set -o nounset
set -o pipefail

scriptName=$(basename $0)

#
# unzip or cat file as required:
#
optionalUngzip() {
    infile=$1
    if file --dereference -b $infile | grep -q gzip; then
        gzip -dc $infile
    else
        cat $infile
    fi
}

filterHeader() {
    awk '!/^#/'
}


#
# Optional add-on filters
#
stripQSX() {
    sed "s/QSS=[0-9]*//g" | sed "s/QSS_NT=[0-9]*//g" |\
    sed "s/QSI=[0-9]*//g" | sed "s/QSI_NT=[0-9]*//g"
}

stripEVS() {
    sed "s/EVS=[0-9]*//g"
}

stripEVSF() {
    sed "s/EVSF=[0-9.,e\-]*//g"
}

stripQual() {
    awk 'BEGIN {FS="\t"; OFS="\t";} {$6=""; print}'
}

stripSample() {
    awk 'BEGIN {FS="\t"; OFS="\t";} {for(i=10;i<=NF;++i) { $i=""; } print;}'
}



#
# optionally ungzip, then remove header and other fields from each vcf
#
stripVcfCore() {
    infile=$1

    optionalUngzip $infile |\
    filterHeader |\
    awk '$6!="."' |\
    stripEVSF |\
    cat
}


#
# strip vcf and add additional info to make diff more informative:
#
stripVcf() {
    infile=$1

    # print input filename so that it's easy to figure out the diff polarity and see global lineCount diff:
    echo "$scriptName filteredVcf: $infile"
    lineCount=$(stripVcfCore $infile | wc -l)
    echo "$scriptName variantLineCount: $lineCount"

    # extra spaces keep the filename/lineCount diff above from attaching to a change on line 1 of the file:
    echo -e "\n\n\n"

    stripVcfCore $infile
}



#
# parse cmdline and diff:
#
if [ $# != 2 ]; then
    cat <<END
usage: $0 file1.vcf[.gz] file2.vcf[.gz]

This script helps to compare two strelka germline vcf files. The vcfs can be
bgziped, gziped or uncompressed. Each file will be uncompressed and the
header + various nuisance felds will be stripped out.
END
    exit 2
fi

file1=$1
file2=$2

diff -U 0 <(stripVcf $file1) <(stripVcf $file2)


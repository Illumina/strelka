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
replace vcf header cmdline field

usage $0 "new cmdline" < in > out
"""

import os, sys

prefix="##cmdline="


def main() :
    infp = sys.stdin
    outfp = sys.stdout
    
    for line in infp :
        if line.startswith(prefix):
            line = prefix + " ".join(sys.argv[1:]) + "\n"
        outfp.write(line)


main()

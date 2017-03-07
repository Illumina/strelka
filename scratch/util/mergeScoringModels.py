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

import collections
import json
import sys


def update(d, u):
    """
    recursive merge of u into d
    """
    for k, v in u.iteritems():
        if isinstance(v, collections.Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            assert(k not in d)
            d[k] = u[k]
    return d

def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog [file1.json [file2.json [...]]] > out.json"
    description="""
merge json files
"""
    parser = OptionParser(usage=usage)

    (opt,args) = parser.parse_args()

    if len(args) == 0 :
        parser.print_help()
        sys.exit(2)

    return (opt,args)

def main() :

    (opt,args) = getOptions()

    data = {}
    for infile in args :
        sys.stderr.write("starting file: %s\n" % (infile))
        indata = json.load(open(infile))
        update(data, indata)

    json.dump(data,sys.stdout)

main()


#!/usr/bin/env python

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


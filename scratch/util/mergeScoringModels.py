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

data = {} 
for infile in sys.argv[1:] :
    sys.stderr.write("starting file: %s\n" % (infile))
    indata = json.load(open(infile))
    update(data, indata)

json.dump(data,sys.stdout)


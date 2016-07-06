#!/usr/bin/env python

import sys


def getOptions() :

    from optparse import OptionParser

    usage = "usage: %prog -schema FILE [options] < in.json"
    description="""
validate json file against specified schema
"""
    parser = OptionParser(usage=usage)

    parser.add_option("--schema", type="string",
                      help="JSON schema file (required)")

    (opt,args) = parser.parse_args()

    if len(args) != 0 :
        parser.print_help()
        sys.exit(2)

    # validate input:
    if opt.schema is None:
        raise Exception("JSON schema file is required")

    return (opt,args)



def main() :
    import json,jsonschema

    (opt,args) = getOptions()

    schemaData = json.load(open(opt.schema))
    jsonschema.Draft4Validator.check_schema(schemaData)
    jsonschema.validate(json.load(sys.stdin), schemaData)


main()


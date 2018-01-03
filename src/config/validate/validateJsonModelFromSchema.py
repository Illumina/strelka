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


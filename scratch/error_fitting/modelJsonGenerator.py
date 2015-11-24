#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

# 12/11/2014
#
# Formatting of model json files for Strelka/Starling variant callers
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Morten Kallberg <mkallberg@illumina.com>
#

import json
import argparse
import logging


def add_sections_indel(json):
    # 'Xten', 'Hiseq', 'HiseqPCR', u'Nextseq'
    models = {'Hiseq':hiseqBase(),'HiseqPCR':pcrBase(),'Nextseq':novaBase(),'Xten':xtenBase()}
    for model in models.keys():
        add_section(json,model,models[model])
        

def add_sections_filtering():
    e=1

def add():
    e=1

def add_section(json,name,modelSection,type='IndelModels'):
    del json['IndelModels'][name]
    json['IndelModels'][name] = modelSection
    

#parser = OptionParser(usage='%prog <OPTIONS>')
#parser.add_option('--force-interactive', dest='force_interactive', help='.')
#parser.add_option('--ini-file', dest='ini_file', help='')

#opt, args = parser.parse_args()

parser.add_argument("truth", help="Truth VCF file")
parser.add_argument("query", help="Query VCF file")

parser.add_argument("-o", "--output", dest="output", required=True,
                    help="Output file prefix for statistics and feature table (when selected)")

parser.add_argument("-l", "--location", dest="location", default="",
                    help="Location for bcftools view (e.g. chr1)")

parser.add_argument("-f", "--false-positives", dest="FP", 
                    help="False-positive region bed file to distinguish UNK from FP")


inModel  = '/home/mkallberg/workspace/starka/src/config/indel_models.json' 
outModel = '/home/mkallberg/workspace/starka/src/config/indel_models2.json' 
f = open(inModel)
j = json.load(f)
f.close()

add_sections_indel(j)
for k in  j['IndelModels'].keys():
    print k
    print j['IndelModels'][k][0]


f = open(outModel,'w')
json.dump(j,f)
f.close()

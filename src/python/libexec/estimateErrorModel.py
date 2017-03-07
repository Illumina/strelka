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
import os
from math import *
import json

def error(hpol_len,indel_length=1):
    insert_A = 1.49133831e-03
    insert_B = 1.03348683e+01
    insert_C = 1.13646811e+00
    insert_D = 1.18488282e-05

    delete_A = 1.09573511e-03
    delete_B = 9.82226042e+00
    delete_C = 1.03579658e+00
    delete_D = 8.31843836e-06

    insert_g = insert_A/(1 + exp((insert_B-hpol_len)/insert_C))+insert_D
    insert_error_prob= (1.- exp(-insert_g/hpol_len))

    delete_g = delete_A/ (1 + exp((delete_B-hpol_len)/delete_C))+delete_D
    delete_error_prob=(1.-exp(-delete_g/hpol_len))

    return [insert_error_prob,delete_error_prob]




def writeModelFile(path,models={}):
    f = open(os.path.join(path,'model_test.json'),'w')
    s = {'IndelModels':{'Hiseq':[],'HiseqPCR':[],'Nextseq':[],'Xten':[],},}
    s['CalibrationModels'] = {'qscorehpol':[]}
    for i in xrange(1,41):
        s['IndelModels']['Hiseq'].append(error(i))
        s['IndelModels']['HiseqPCR'].append(error(i))
        s['IndelModels']['Nextseq'].append(error(i))
        s['IndelModels']['Xten'].append(error(i))
        s['CalibrationModels']['qscorehpol'].append(0.1)
    json.dump(s,f)
    f.close()

writeModelFile('/home/mkallberg/workspace/starka_debug/data/')
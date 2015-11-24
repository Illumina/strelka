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
import os.path
import shutil
import sys
import subprocess
from dynamicCaseCounter import *
import pickle
import numpy as np
import matplotlib.pyplot as plt
import json

base = '/home/mkallberg/new_indel_models/'
models  = ['Xten', 'Hiseq', 'HiseqPCR', 'Nextseq']
pars = {}
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
colors = ['Red','Blue','Green','Gold',"Cyan"]

for i,m in enumerate(models):
    case = pickle.load(open(os.path.join(base,m+'.pickle')))
    case.calcPars('del')
    case.calcPars('ins')
    model = 'simple'
    caseInd = 'ins'
    pars[m] =  {'ins':case.parsHpol['ins'][model]['props'],'del':case.parsHpol['del'][model]['props']}
    x = np.linspace(1, len(pars[m][caseInd]),len(pars[m][caseInd]))
    ax.scatter(x, pars[m][caseInd],label=m,color=colors[i])

plt.ylim([1e-6,5e-2])
ax.set_yscale('log')
ax.legend()
plt.show() 

inModel  = '/home/mkallberg/workspace/starka/src/config/indel_models.json' 
outModel = os.path.join(base,'indel_models.json') 
f = open(inModel)
j = json.load(f)
f.close()

for model in pars.keys():
    m = []
    for i in xrange(0,34):
        m.append([pars[model]['ins'][i],pars[model]['del'][i]])
    j['IndelModels'][model] = m
    print j['IndelModels'][model]


f = open(outModel,'w')
json.dump(j,f)
f.close()
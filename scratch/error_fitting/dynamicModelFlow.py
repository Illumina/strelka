#!/illumina/thirdparty/python/python-2.7.5/bin/python
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

########################################
# Pyflow for running error-model estimator in parallele
#
# Author: Morten Kallberg
# Date: September 16th, 2014
#
########################################

DEBUG=0

import os.path
import shutil
import sys
import subprocess
import random

# add script path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(scriptDir,"pyflow"))
from pyflow import WorkflowRunner
#
def massageTaskName(s):
    return s.replace(':',"_")


class indelErrorWorkflow(WorkflowRunner) :
    def __init__(self, bamFile, tempPath, outModel, reference, samtools,
                 depth = None, range = None, inModel = None) :
        self.bamFile    = bamFile
        self.tempPath   = tempPath
        self.inModel    = inModel
        self.outModel   = outModel
        self.reference  = reference
        self.depth      = depth
        self.samtools   = samtools
        self.range      = range     # maximum walking range
#        self.workflowSetup()

    def generateRegionCandidates(self,randomSelection=0,cases=100,regionSize=500000,totalMB=10):
        '''
        Generate a list of candidate regions to estimate error profile from by sampling random reads from bam file.
        '''
        # list that tries to steer clear of centromers replace with a mean depth/ cleanliness criteria
#        regions = ['chr1:10000-120000000','chr1:120000000-240000000','chr2:10000-90000000','chr2:100000000-240000000',\
#        'chr3:10000-88000000','chr3:94000000-200000000','chr4:10000-43000000','chr4:50000000-182000000','chr5:10000-45000000',\
#        'chr5:50000000-178000000','chr6:10000-58000000','chr6:62000000-165000000','chr7:10000-55000000','chr7:65000000-158000000',\
#        'chr8:10000-42000000','chr8:48000000-140000000']
        regions = ['chr20:8700000-9100000','chr20:9100000-9500000','chr20:9500000-9900000','chr20:9900000-10200000',]
        if randomSelection:
            regions = []
            cmd = [self.samtools, 'view', '-H', self.bamFile]
            process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
            process.wait()
            ranges = []
            for l in process.stdout.read().split('@'):
                s = l.split()
                if len(s) and s[0]=='SQ':
                    ranges.append([s[1].replace('SN:',''),int(s[2].replace('LN:',''))])
    
            for i in xrange(0,cases):
                randomChr = random.choice(ranges)
                myChr = randomChr[0]
                myStart = int(random.random()*randomChr[1])
                myEnd = myStart + regionSize
                if not myChr in ['chrX','chrY','chrM']:
                    regions.append(myChr + ":"  + str(myStart) + "-" + str(myEnd))
        return regions

    def workflow(self) :
        callPreReqs = set()
        for c in self.generateRegionCandidates(randomSelection=1):
            cmd = self.generateJob(c)
            callPreReqs.add(self.addTask("Estimate_" + massageTaskName(c)," ".join(cmd),dependencies=None) )
        cmd = ['/illumina/thirdparty/python/python-2.7.5/bin/python',os.path.join(scriptDir,'dynamicModelAggregator.py'), '-i', self.tempPath,'-o',self.outModel]
        if self.inModel is not None:
            cmd.append('-m', self.inModel)
        print " ".join(cmd)
        self.addTask("JoinIndelModel"," ".join(cmd),dependencies=callPreReqs)

    def generateJob(self,region):
        cmd = ['/illumina/thirdparty/python/python-2.7.5/bin/python',os.path.join(scriptDir,'dynamicModel.py')]
        cmd += ['-b',self.bamFile]
        cmd += ['-o',self.tempPath]
        cmd += ['-r',self.reference]
        cmd += ['-g',region]
        return cmd


if __name__ == "__main__":
    # w = indelErrorWorkflow('/home/mkallberg/workspace/starka_debug/bams/error_model/NA12878_hiseq_nopcr.bam',
    #                        '/home/mkallberg/workspace/starka_debug/bams/error_model/test/',
    #                        '/home/mkallberg/workspace/starka/src/config/indel_models.json',
    #                        '/home/mkallberg/workspace/starka_debug/bams/error_model/test/indel_models_test.json',
    #                        '/home/mkallberg/common/hg_19/genome.fa',
    #                        '/home/mkallberg/workspace/starka_debug/bams/error_model/test/',
    #                        'samtools','samtools','')
    w = indelErrorWorkflow('/home/mbekritsky/workspace/starka/error_model/bams/NA12878PG4L.bam',
                           '/home/mbekritsky/workspace/starka/error_model/BinomTest/',
                           '/home/mbekritsky/workspace/starka/error_model/indel_models.json',
                           '/illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
                           '/illumina/thirdparty/samtools/samtools-1.2/bin/samtools')
    w.run(mode='local', nCores = 3)
    

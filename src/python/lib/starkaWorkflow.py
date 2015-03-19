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

"""
Shared small varaint calling workflow components
"""

import os.path
import sys

# add this path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(scriptDir)

# add pyflow path:
sys.path.append(os.path.join(scriptDir,"pyflow"))


from pyflow import WorkflowRunner
from workflowUtil import checkFile, ensureDir, getFastaChromOrderSize, cleanPyEnv, preJoin


class StarkaCallWorkflow(WorkflowRunner) :

    def __init__(self,params) :
        self.params = params
        
        
    def concatIndexBgzipFile(self, taskPrefix, dependencies, inputList, output, label, fileType) :
        """
        Internal helper function
        @param inputList files to be concatenated (in order), already bgzipped
        @param output output filename
        @param label used for error task id
        @param fileType provided to tabix
        """
        assert(len(inputList) > 0)
    
        if len(inputList) > 1 :
            catCmd = [self.params.bgcatBin,"-o",output]
            catCmd.extend(inputList)
        else :
            catCmd = "mv -f %s %s" % (inputList[0],output)
    
        indexCmd = "%s -p %s %s" % (self.params.tabixBin, fileType, output)
        catTask = self.addTask(preJoin(taskPrefix,label+"_concat_"+fileType), catCmd,
                               dependencies=dependencies, isForceLocal=True)
        return self.addTask(preJoin(taskPrefix,label+"_index_"+fileType), indexCmd,
                                    dependencies=catTask, isForceLocal=True)



    def concatIndexVcf(self, taskPrefix, dependencies, inputList, output, label) :
        """
        Concatenate bgzipped vcf segments
        @param inputList files to be concatenated (in order), already bgzipped
        @param output output filename
        @param label used for error task id
        """
        assert(len(inputList) > 0)
        return self.concatIndexBgzipFile(taskPrefix, dependencies, inputList, output, label, "vcf")



    def concatIndexBed(self, taskPrefix, dependencies, inputList, output, label) :
        """
        Concatenate bgzipped bed segments
        @param inputList files to be concatenated (in order), already bgzipped
        @param output output filename
        @param label used for error task id
        """
        assert(len(inputList) > 0)
        return self.concatIndexBgzipFile(taskPrefix, dependencies, inputList, output, label, "bed")



class StarkaWorkflow(WorkflowRunner) :
    """
    small variant calling workflow base
    """

    def __init__(self,params,iniSections) :

        cleanPyEnv()

        self.params=params
        self.iniSections=iniSections

        # make sure run directory is setup:
        self.params.runDir=os.path.abspath(self.params.runDir)
        ensureDir(self.params.runDir)

        # everything that's not intended to be a final result should dump directories/files in workDir
        self.params.workDir=os.path.join(self.params.runDir,"workspace")
        ensureDir(self.params.workDir)

        # all finalized pretty results get transfered to resultsDir
        self.params.resultsDir=os.path.join(self.params.runDir,"results")
        ensureDir(self.params.resultsDir)
        self.params.variantsDir=os.path.join(self.params.resultsDir,"variants")
        ensureDir(self.params.variantsDir)

        indexRefFasta=self.params.referenceFasta+".fai"

        if self.params.referenceFasta is None:
            raise Exception("No reference fasta defined.")
        else:
            checkFile(self.params.referenceFasta,"reference fasta")
            checkFile(indexRefFasta,"reference fasta index")

        # read fasta index
        (self.params.chromOrder,self.params.chromSizes) = getFastaChromOrderSize(indexRefFasta)

        self.params.isHighDepthFilter = (not self.params.isExome)


    def setCallMemMb(self) :
        # set default mem requirements:
        if self.params.callMemMbOverride is not None :
            self.params.callMemMb = self.params.callMemMbOverride
        else :
            if self.getRunMode() == "sge" :
                self.params.callMemMb = self.params.callSGEMemMb
            else :
                self.params.callMemMb = self.params.callLocalMemMb


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

"""
Strelka noise estimate workflow
"""


import os.path
import sys

# add script path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(scriptDir))

# add pyflow path:
# TODO: get a more robust link to the pyflow dir at config time:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from configBuildTimeInfo import workflowVersion
from pyflow import WorkflowRunner
from strelkaSharedWorkflow import SharedPathInfo, StrelkaSharedWorkflow
from workflowUtil import preJoin, getNextGenomeSegment



__version__ = workflowVersion



class TempSegmentFiles :
    def __init__(self) :
        self.gvcf = []
        self.bamRealign = []



def callGenomeSegment(self, gseg, segFiles, taskPrefix="", dependencies=None) :

    isFirstSegment = (len(segFiles.gvcf) == 0)

    segStr = str(gseg.id)

    segCmd = [ self.params.snoiseBin ]

    segCmd.extend(["--region", gseg.chromLabel + ":" + str(gseg.beginPos) + "-" + str(gseg.endPos)])
    segCmd.extend(["--min-mapping-quality",self.params.minMapq])
    segCmd.extend(["--ref", self.params.referenceFasta ])
    segCmd.extend(["--max-indel-size", self.params.maxIndelSize])

    for bamPath in self.params.bamList :
        segCmd.extend(["--align-file",bamPath])

    if not isFirstSegment :
        segCmd.append("--skip-vcf-header")

    if self.params.indelCandidates is not None :
        segCmd.extend(['--candidate-indel-input-vcf', self.params.indelCandidates])

     # vcf is written to stdout so we need shell features:
    segCmd = " ".join(segCmd)

    segFiles.gvcf.append(self.paths.getTmpSegmentGvcfPath(segStr))
    segCmd += " | %s -c >| %s" % (self.params.bgzip9Bin, segFiles.gvcf[-1])

    nextStepWait = set()

    setTaskLabel=preJoin(taskPrefix,"callGenomeSegment_"+gseg.id)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)
    nextStepWait.add(setTaskLabel)

    return nextStepWait



def callGenome(self,taskPrefix="",dependencies=None):
    """
    run variant caller on all genome segments
    """

    tmpGraphDir=self.paths.getTmpSegmentDir()
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), "mkdir -p "+tmpGraphDir, dependencies=dependencies, isForceLocal=True)

    graphTasks = set()

    segFiles = TempSegmentFiles()
    for gseg in getNextGenomeSegment(self.params) :

        graphTasks |= callGenomeSegment(self, gseg, segFiles, dependencies=dirTask)

    # create a checkpoint for all segments:
    completeSegmentsTask = self.addTask(preJoin(taskPrefix,"completedAllGenomeSegments"),dependencies=graphTasks)

    finishTasks = set()

    def finishVcf(tmpList, output, label) :
        assert(len(tmpList) > 0)

        if len(tmpList) > 1 :
            catCmd=[self.params.bgcatBin,"-o",output]
            catCmd.extend(tmpList)
            catCmd = " ".join(catCmd)
        else :
            catCmd="mv -f %s %s" % (tmpList[0],output)

        catCmd += " && %s -p vcf %s" % (self.params.tabixBin, output)
        finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeVCF"), catCmd, dependencies=completeSegmentsTask))

    finishVcf(segFiles.gvcf, self.paths.getGvcfOutputPath(),"gVCF")

    cleanTask=self.addTask(preJoin(taskPrefix,"cleanTmpDir"), "rm -rf "+tmpGraphDir, dependencies=finishTasks, isForceLocal=True)

    nextStepWait = finishTasks

    return nextStepWait



class CallWorkflow(WorkflowRunner) :
    """
    A separate call workflow is setup so that we can delay the workflow execution until
    the ref count file exists
    """

    def __init__(self,params,paths) :
        self.params = params
        self.paths = paths

    def workflow(self) :
        callGenome(self)



class PathInfo(SharedPathInfo):
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        super(PathInfo,self).__init__(params)

    def getTmpSegmentDir(self) :
        return os.path.join(self.params.workDir, "genomeSegment.tmpdir")

    def getTmpSegmentGvcfPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "noise.%s.vcf.gz" % (segStr))

    def getGvcfOutputPath(self) :
        return os.path.join( self.params.variantsDir, "noise.vcf.gz")



class snoiseWorkflow(StrelkaSharedWorkflow) :
    """
    germline small variant calling workflow
    """

    def __init__(self,params) :
        global PathInfo
        super(snoiseWorkflow,self).__init__(params, PathInfo)

        # format bam lists:
        if self.params.bamList is None : self.params.bamList = []

        # sanity check some parameter typing:
        MEGABASE = 1000000
        self.params.scanSize = int(self.params.scanSizeMb) * MEGABASE


    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Strelka noise estimation workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Strelka noise estimation workflow version: %s" % (__version__))

        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths))

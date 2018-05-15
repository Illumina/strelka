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
Sequence error counts workflow
"""


import os.path
import sys

# add this path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(scriptDir)

# add pyflow path:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from configBuildTimeInfo import workflowVersion
from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getRmdirCmd, getDepthFromAlignments
from strelkaSharedWorkflow import SharedPathInfo, \
                           StrelkaSharedCallWorkflow, StrelkaSharedWorkflow
from workflowUtil import ensureDir, preJoin, getNextGenomeSegment


__version__ = workflowVersion



def strelkaGermlineRunDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    bamList=[]
    if len(self.params.bamList) :
        bamList.append(self.params.bamList[0])
    else :
        return set()

    outputPath=self.paths.getChromDepth()
    return getDepthFromAlignments(self, bamList, outputPath, taskPrefix, dependencies)



class TempSegmentFiles :
    def __init__(self) :
        self.counts = []
        self.observedIndelBed = []



def callGenomeSegment(self, gseg, segFiles, taskPrefix="", dependencies=None) :

    segStr = str(gseg.id)

    segCmd = [ self.params.getCountsBin ]

    segCmd.extend(["--region", gseg.chromLabel + ":" + str(gseg.beginPos) + "-" + str(gseg.endPos)])
    segCmd.extend(["--ref", self.params.referenceFasta ])
    segCmd.extend(["--max-indel-size", self.params.maxIndelSize])

    segFiles.counts.append(self.paths.getTmpSegmentAlleleCountsPath(segStr))
    segCmd.extend(["--counts-file", segFiles.counts[-1]])

    for bamPath in self.params.bamList :
        segCmd.extend(["--align-file",bamPath])

    if self.params.isHighDepthFilter :
        segCmd.extend(["--chrom-depth-file", self.paths.getChromDepth()])

    def addListCmdOption(optList,arg) :
        if optList is None : return
        for val in optList :
            segCmd.extend([arg, val])

    addListCmdOption(self.params.indelCandidatesList, '--candidate-indel-input-vcf')
    addListCmdOption(self.params.forcedGTList, '--force-output-vcf')

    addListCmdOption(self.params.excludedRegions,"--excluded-regions-bed-file")
    if self.params.knownVariants is not None :
        segCmd.extend(["--known-variants-vcf-file",self.params.knownVariants])

    if self.params.isReportObservedIndels :
        tmpObservedIndelBedPath = self.paths.getTmpObservedIndelBedPath(segStr)
        segFiles.observedIndelBed.append(tmpObservedIndelBedPath + ".gz")
        segCmd.extend(['--observation-bed-file', tmpObservedIndelBedPath])

    if self.params.extraCountsArguments is not None :
        for arg in self.params.extraCountsArguments.strip().split() :
            segCmd.append(arg)

    nextStepWait = set()

    setTaskLabel=preJoin(taskPrefix,"countGenomeSegment_"+gseg.id)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)
    nextStepWait.add(setTaskLabel)

    if self.params.isReportObservedIndels :
        compressTask=preJoin(taskPrefix,"compressSegmentOutput_"+gseg.id)
        compressCmd = "\"%s\" \"%s\"" % (self.params.bgzipBin, tmpObservedIndelBedPath)
        self.addTask(compressTask, compressCmd, dependencies=setTaskLabel, isForceLocal=True)
        nextStepWait.add(compressTask)

    return nextStepWait


def mergeSequenceAlleleCounts(self, taskPrefix, dependencies, runStatsLogPaths) :

    runMergeLabel=preJoin(taskPrefix,"mergeCounts")
    runMergeCmd=[self.params.mergeCountsBin]
    for statsFile in runStatsLogPaths :
        runMergeCmd.extend(["--counts-file",statsFile])
    runMergeCmd.extend(["--output-file", self.paths.getAlleleCountsOutputPath()])
    return self.addTask(runMergeLabel, runMergeCmd, dependencies=dependencies, isForceLocal=True)



def callGenome(self,taskPrefix="",dependencies=None):
    """
    run counter on all genome segments
    """

    tmpSegmentDir=self.paths.getTmpSegmentDir()
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), getMkdirCmd() + [tmpSegmentDir], dependencies=dependencies, isForceLocal=True)

    segmentTasks = set()

    segFiles = TempSegmentFiles()
    for gseg in getNextGenomeSegment(self.params) :

        segmentTasks |= callGenomeSegment(self, gseg, segFiles, dependencies=dirTask)

    if len(segmentTasks) == 0 :
        raise Exception("No genome regions to analyze. Possible target region parse error.")

    # create a checkpoint for all segments:
    completeSegmentsTask = self.addTask(preJoin(taskPrefix,"completedAllGenomeSegments"),dependencies=segmentTasks)

    finishTasks = set()

    # merge segment stats:
    finishTasks.add(mergeSequenceAlleleCounts(self, taskPrefix, completeSegmentsTask, segFiles.counts))

    if self.params.isReportObservedIndels :
        finishTasks.add(self.concatIndexBed(taskPrefix, completeSegmentsTask, segFiles.observedIndelBed,
                                            self.paths.getObservedIndelBedPath(), "observedIndels"))

    if not self.params.isRetainTempFiles :
        rmTmpCmd = getRmdirCmd() + [tmpSegmentDir]
        rmTask=self.addTask(preJoin(taskPrefix,"rmTmpDir"),rmTmpCmd,dependencies=finishTasks, isForceLocal=True)

    nextStepWait = finishTasks

    return nextStepWait



class CallWorkflow(StrelkaSharedCallWorkflow) :
    """
    A separate call workflow is setup so that we can delay the workflow execution until
    the ref count file exists
    """

    def __init__(self,params,paths) :
        super(CallWorkflow,self).__init__(params)
        self.paths = paths

    def workflow(self) :
        callGenome(self)



class PathInfo(SharedPathInfo):
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        super(PathInfo,self).__init__(params)

    def getTmpSegmentAlleleCountsPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "strelkaAlleleCounts.%s.bin" % (segStr))

    def getTmpObservedIndelBedPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "strelkaObservedIndel.%s.bed" % segStr)

    def getAlleleCountsOutputPath(self) :
        return os.path.join( self.params.variantsDir, "strelkaAlleleCounts.bin")

    def getObservedIndelBedPath(self) :
        return os.path.join( self.params.debugDir, "strelkaObservedIndel.bed.gz")



class SequenceAlleleCountsWorkflow(StrelkaSharedWorkflow) :
    """
    Sequence allele counts workflow
    """

    def __init__(self, params) :
        global PathInfo
        super(SequenceAlleleCountsWorkflow, self).__init__(params, PathInfo)

        # if debugging output is going to be produced, add a results/debug dir
        if self.params.isReportObservedIndels:
            self.params.debugDir=os.path.join(self.params.resultsDir,"debug")
            ensureDir(self.params.debugDir)

        # format bam lists:
        if self.params.bamList is None : self.params.bamList = []


    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Strelka sequence error counts workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg


    def workflow(self) :
        self.flowLog("Initiating Strelka sequence error counts workflow version: %s" % (__version__))
        self.setCallMemMb()

        callPreReqs = set()
        if self.params.isHighDepthFilter :
            callPreReqs |= strelkaGermlineRunDepthFromAlignments(self)
        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)

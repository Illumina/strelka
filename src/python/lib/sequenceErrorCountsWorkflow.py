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
from configureUtil import safeSetBool, getIniSections, dumpIniSections
from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getRmdirCmd, runDepthFromAlignments
from starkaWorkflow import runCount, SharedPathInfo, \
                           StarkaCallWorkflow, StarkaWorkflow
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getNextGenomeSegment, bamListCatCmd


__version__ = workflowVersion



def starlingRunDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    bamList=[]
    if len(self.params.bamList) :
        bamList.append(self.params.bamList[0])
    else :
        return set()

    outputPath=self.paths.getChromDepth()
    return runDepthFromAlignments(self, bamList, outputPath, taskPrefix, dependencies)



class TempSegmentFiles :
    def __init__(self) :
        self.counts = []



def callGenomeSegment(self, gseg, segFiles, taskPrefix="", dependencies=None) :

    segStr = str(gseg.id)

    segCmd = [ self.params.getCountsBin ]

    segCmd.append("-clobber")
    segCmd.extend(["-bam-seq-name", gseg.chromLabel] )
    segCmd.extend(["-report-range-begin", str(gseg.beginPos) ])
    segCmd.extend(["-report-range-end", str(gseg.endPos) ])
    segCmd.extend(["-samtools-reference", self.params.referenceFasta ])
    segCmd.extend(["-genome-size", str(self.params.knownSize)] )
    segCmd.extend(["-max-indel-size", "50"] )

    segFiles.counts.append(self.paths.getTmpSegmentCountsPath(segStr))
    segCmd.extend(["--counts-file", segFiles.counts[-1]])

    for bamPath in self.params.bamList :
        segCmd.extend(["-bam-file", bamPath])

    if self.params.isHighDepthFilter :
        segCmd.extend(["--chrom-depth-file", self.paths.getChromDepth()])

    def addListCmdOption(optList,arg) :
        if optList is None : return
        for val in optList :
            segCmd.extend([arg, val])

    addListCmdOption(self.params.indelCandidatesList, '--candidate-indel-input-vcf')
    addListCmdOption(self.params.forcedGTList, '--force-output-vcf')

    if self.params.targetRegionsBed is not None :
        segCmd.extend(['--targeted-regions-bed', self.params.targetRegionsBed])

    if self.params.ploidyBed is not None :
        segCmd.extend(['--ploidy-region-bed', self.params.ploidyBed])

    if self.params.extraCountsArguments is not None :
        for arg in self.params.extraCountsArguments.strip().split() :
            segCmd.append(arg)

    nextStepWait = set()

    setTaskLabel=preJoin(taskPrefix,"countGenomeSegment_"+gseg.id)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)
    nextStepWait.add(setTaskLabel)

    return nextStepWait



def mergeSequenceErrorCounts(self, taskPrefix, dependencies, runStatsLogPaths) :

    runMergeLabel=preJoin(taskPrefix,"mergeCounts")
    runMergeCmd=[self.params.mergeCountsBin]
    for statsFile in runStatsLogPaths :
        runMergeCmd.extend(["--counts-file",statsFile])
    runMergeCmd.extend(["--output-file",self.paths.getCountsOutputPath()])
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
    finishTasks.add(mergeSequenceErrorCounts(self,taskPrefix,completeSegmentsTask, segFiles.counts))

    if not self.params.isRetainTempFiles :
        rmTmpCmd = getRmdirCmd() + [tmpSegmentDir]
        rmTask=self.addTask(preJoin(taskPrefix,"rmTmpDir"),rmTmpCmd,dependencies=finishTasks, isForceLocal=True)

    nextStepWait = finishTasks

    return nextStepWait



"""
A separate call workflow is setup so that we can delay the workflow execution until
the ref count file exists
"""
class CallWorkflow(StarkaCallWorkflow) :

    def __init__(self,params,paths) :
        super(CallWorkflow,self).__init__(params)
        self.paths = paths

    def workflow(self) :

        if True :
            knownSize = 0
            for line in open(self.paths.getRefCountFile()) :
                word = line.strip().split('\t')
                if len(word) != 4 :
                    raise Exception("Unexpected format in ref count file: '%s'" % (self.paths.getRefCountFile()))
                knownSize += int(word[2])

            self.params.knownSize = knownSize

        callGenome(self)



class PathInfo(SharedPathInfo):
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        super(PathInfo,self).__init__(params)

    def getTmpSegmentCountsPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "strelkaErrorCounts.%s.bin" % (segStr))

    def getCountsOutputPath(self) :
        return os.path.join( self.params.variantsDir, "strelkaErrorCounts.bin")



class SequenceErrorCountsWorkflow(StarkaWorkflow) :
    """
    sequence error counts workflow
    """

    def __init__(self,params,iniSections) :
        global PathInfo
        super(SequenceErrorCountsWorkflow,self).__init__(params,iniSections,PathInfo)

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
        callPreReqs |= runCount(self)
        if self.params.isHighDepthFilter :
            callPreReqs |= starlingRunDepthFromAlignments(self)
        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)


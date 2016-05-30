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
Starling germline small variant calling workflow
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



def runIndelModel(self,taskPrefix="",dependencies=None) :
    """
    estimate indel error parameters and write back a modified model file
    """

    bamFile=""
    if len(self.params.bamList) :
        bamFile = self.params.bamList[0]
    else :
        return set()
    tempPath  = self.params.getIndelSegmentDir
    inModel   = self.params.inputIndelErrorModelsFile
    outModel  = self.params.getRunSpecificModel
    reference = self.params.referenceFasta
    scriptDir = os.path.abspath(scriptDir)
    depth     = self.params.getChromDepth

    nextStepWait = set()
    nextStepWait.add(self.addWorkflowTask("GenerateIndelModel", indelErrorWorkflow(bamFile,tempPath,inModel,outModel,reference,scriptDir,depth), dependencies=dependencies))

    # edit the scoring model used to reflect the restated model
    self.params.dynamicIndelErrorModelsFile = outModel

    return nextStepWait


class TempSegmentFiles :
    def __init__(self) :
        self.gvcf = []
        self.bamRealign = []
        self.stats = []


# we need extra quoting for files with spaces in this workflow because command is stringified below to enable gVCF pipe:
def quote(instr):
    return "\"%s\"" % (instr)



def callGenomeSegment(self, gseg, segFiles, taskPrefix="", dependencies=None) :

    isFirstSegment = (len(segFiles.gvcf) == 0)

    segStr = str(gseg.id)

    segCmd = [ quote(self.params.starlingBin) ]

    segCmd.append("-clobber")
    segCmd.extend(["-min-mapping-quality",self.params.minMapq])
    segCmd.extend(["-bam-seq-name", gseg.chromLabel] )
    segCmd.extend(["-report-range-begin", str(gseg.beginPos) ])
    segCmd.extend(["-report-range-end", str(gseg.endPos) ])
    segCmd.extend(["-samtools-reference", quote(self.params.referenceFasta) ])
    segCmd.extend(["-max-window-mismatch", "2", "20" ])
    segCmd.extend(["-genome-size", str(self.params.knownSize)] )
    segCmd.extend(["-max-indel-size", "50"] )

    segCmd.extend(["--gvcf-file","-"])
    segCmd.extend(['--gvcf-min-gqx','15'])
    segCmd.extend(['--gvcf-max-snv-strand-bias','10'])
    segCmd.extend(['--gvcf-max-indel-ref-repeat', '-1'])
    segCmd.extend(['-min-qscore','17'])
    segCmd.extend(['-bsnp-ssd-no-mismatch', '0.35'])
    segCmd.extend(['-bsnp-ssd-one-mismatch', '0.6'])
    segCmd.extend(['-min-vexp', '0.25'])
    segCmd.extend(['--do-short-range-phasing'])
    # currently short-range phasing is not enabled with continuous variant calling. This ensures
    # the header value for the relevant phasing filters is still emitted
    if len(self.params.callContinuousVf) > 0 :
        segCmd.extend(["--gvcf-include-header", "Phasing"])
    segCmd.extend(["--report-file", quote(self.paths.getTmpSegmentReportPath(gseg.id))])

    segFiles.stats.append(self.paths.getTmpRunStatsPath(segStr))
    segCmd.extend(["--stats-file", quote(segFiles.stats[-1])])

    # Empirical Variant Scoring(EVS):
    if self.params.isEVS :
        segCmd.extend(['--variant-scoring-models-file',quote(self.params.evsModelFile)])
        segCmd.extend(['--variant-scoring-model-name',self.params.evsModelName])

    if self.params.indelErrorModelName is not None :
        segCmd.extend(['--indel-error-model-name',self.params.indelErrorModelName])
    if self.params.inputIndelErrorModelsFile is not None :
        segCmd.extend(['--indel-error-models-file', quote(self.params.inputIndelErrorModelsFile)])

    if self.params.isReportEVSFeatures :
        segCmd.append("--report-evs-features")

    for bamPath in self.params.bamList :
        segCmd.extend(["-bam-file",quote(bamPath)])

    if not isFirstSegment :
        segCmd.append("--gvcf-skip-header")
    elif len(self.params.callContinuousVf) > 0 :
        segCmd.extend(["--gvcf-include-header", "VF"])

    if self.params.isHighDepthFilter :
        segCmd.extend(["--chrom-depth-file", quote(self.paths.getChromDepth())])

    if self.params.isWriteRealignedBam :
        segCmd.extend(["-realigned-read-file", quote(self.paths.getTmpUnsortRealignBamPath(segStr))])

    def addListCmdOption(optList,arg) :
        if optList is None : return
        for val in optList :
            segCmd.extend([arg, val])

    addListCmdOption(self.params.indelCandidatesList, '--candidate-indel-input-vcf')
    addListCmdOption(self.params.forcedGTList, '--force-output-vcf')

    if self.params.noCompressBed is not None :
        segCmd.extend(['--nocompress-bed', quote(self.params.noCompressBed)])

    if self.params.targetRegionsBed is not None :
        segCmd.extend(['--targeted-regions-bed', quote(self.params.targetRegionsBed)])

    if self.params.ploidyBed is not None :
        segCmd.extend(['--ploidy-region-bed', quote(self.params.ploidyBed)])

    if self.params.callContinuousVf is not None and gseg.chromLabel in self.params.callContinuousVf :
        segCmd.append('--call-continuous-vf')

    if self.params.extraStarlingArguments is not None :
        for arg in self.params.extraStarlingArguments.strip().split() :
            segCmd.append(arg)

     # gvcf is written to stdout so we need shell features:
    segCmd = " ".join(segCmd)

    # swap parent pyflow command-line into vcf header
    if isFirstSegment :
        def getHeaderFixCmd() :
            cmd  = "\"%s\" -E \"%s\"" % (sys.executable, self.params.vcfCmdlineSwapper)
            cmd += ' "' + " ".join(self.params.configCommandLine) + '"'
            return cmd

        segCmd += " | " + getHeaderFixCmd()

    segFiles.gvcf.append(self.paths.getTmpSegmentGvcfPath(segStr))
    segCmd += " | \"%s\" -c >| \"%s\"" % (self.params.bgzip9Bin, segFiles.gvcf[-1])

    nextStepWait = set()

    setTaskLabel=preJoin(taskPrefix,"callGenomeSegment_"+gseg.id)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)
    nextStepWait.add(setTaskLabel)

    if self.params.isWriteRealignedBam :
        def sortRealignBam(sortList) :
            unsorted = self.paths.getTmpUnsortRealignBamPath(segStr)
            sorted   = self.paths.getTmpRealignBamPath(segStr)
            sortList.append(sorted)

            # adjust sorted to remove the ".bam" suffix
            sorted = sorted[:-4]
            sortCmd="\"%s\" sort \"%s\" \"%s\" && rm -f \"%s\"" % (self.params.samtoolsBin,unsorted,sorted,unsorted)

            sortTaskLabel=preJoin(taskPrefix,"sortRealignedSegment_"+gseg.id)
            self.addTask(sortTaskLabel,sortCmd,dependencies=setTaskLabel,memMb=self.params.callMemMb)
            nextStepWait.add(sortTaskLabel)

        sortRealignBam(segFiles.bamRealign)

    return nextStepWait



def callGenome(self,taskPrefix="",dependencies=None):
    """
    run variant caller on all genome segments
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

    # merge gVCF
    finishTasks.add(self.concatIndexVcf(taskPrefix, completeSegmentsTask, segFiles.gvcf, self.paths.getGvcfOutputPath(),"gVCF"))

    # merge segment stats:
    finishTasks.add(self.mergeRunStats(taskPrefix,completeSegmentsTask, segFiles.stats))

    if self.params.isWriteRealignedBam :
        def finishBam(tmpList, output, label) :
            cmd = bamListCatCmd(self.params.samtoolsBin, tmpList, output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))

        finishBam(segFiles.bamRealign, self.paths.getRealignedBamPath(), "realigned")

    if not self.params.isRetainTempFiles :
        rmStatsTmpCmd = getRmdirCmd() + [tmpSegmentDir]
        rmTask=self.addTask(preJoin(taskPrefix,"rmTmpDir"),rmStatsTmpCmd,dependencies=finishTasks, isForceLocal=True)

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

    def getRunSpecificModel(self) :
        return os.path.join(self.params.workDir,"Indel_model_run.json")

    def getIndelSegmentDir(self) :
        return os.path.join(self.params.workDir, "indelSegment.tmpdir")

    def getTmpSegmentGvcfPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "genome.%s.vcf.gz" % (segStr))

    def getTmpUnsortRealignBamPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "%s.unsorted.realigned.bam" % (segStr))

    def getTmpRealignBamPath(self, segStr,) :
        return os.path.join( self.getTmpSegmentDir(), "%s.realigned.bam" % (segStr))

    def getGvcfOutputPath(self) :
        return os.path.join( self.params.variantsDir, "genome.vcf.gz")

    def getRealignedBamPath(self) :
        return os.path.join( self.params.realignedDir, 'realigned.bam');



class StarlingWorkflow(StarkaWorkflow) :
    """
    germline small variant calling workflow
    """

    def __init__(self,params,iniSections) :
        global PathInfo
        super(StarlingWorkflow,self).__init__(params,iniSections,PathInfo)

        # format bam lists:
        if self.params.bamList is None : self.params.bamList = []

        # format other:
        safeSetBool(self.params,"isWriteRealignedBam")

        if self.params.isWriteRealignedBam :
            self.params.realignedDir=os.path.join(self.params.resultsDir,"realigned")
            ensureDir(self.params.realignedDir)

        if self.params.isExome :
            self.params.isEVS = False


    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Starling workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Starling workflow version: %s" % (__version__))
        self.setCallMemMb()

        callPreReqs = set()
        callPreReqs |= runCount(self)
        if self.params.isHighDepthFilter :
            callPreReqs |= starlingRunDepthFromAlignments(self)
        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)


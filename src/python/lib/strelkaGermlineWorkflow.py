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
Strelka germline small variant calling workflow
"""


import os.path
import sys

# add this path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(scriptDir)

# add pyflow path:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from pyflow import WorkflowRunner, LogState
from configBuildTimeInfo import workflowVersion
from configureUtil import safeSetBool
from sharedWorkflow import getMkdirCmd, getRmdirCmd, getDepthFromAlignments
from strelkaSequenceErrorEstimation import EstimateSequenceErrorWorkflow
from strelkaSharedWorkflow import SharedPathInfo, \
                           StrelkaSharedCallWorkflow, StrelkaSharedWorkflow
from workflowUtil import ensureDir, preJoin, bamListCatCmd

__version__ = workflowVersion



def strelkaGermlineGetDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    if len(self.params.bamList) == 0 :
        return set()

    bamList = self.params.bamList
    outputPath=self.paths.getChromDepth()
    return getDepthFromAlignments(self, bamList, outputPath, taskPrefix, dependencies)



class TempVariantCallingSegmentFilesPerSample :
    def __init__(self) :
        self.gvcf = []
        self.bamRealign = []


class TempVariantCallingSegmentFiles :
    def __init__(self, sampleCount) :
        self.variants = []
        self.stats = []
        self.sample = [TempVariantCallingSegmentFilesPerSample() for _ in range(sampleCount)]


# we need extra quoting for files with spaces in this workflow because some commands are stringified as shell calls:
def quote(instr):
    return "\"%s\"" % (instr)


def gvcfSampleLabel(sampleIndex) :
    return "gVCF_S%i" % (sampleIndex+1)


def callGenomeSegment(self, gsegGroup, segFiles, taskPrefix="", dependencies=None) :

    assert(len(gsegGroup) != 0)
    genomeSegmentLabel = gsegGroup[0].id
    if len(gsegGroup) > 1 :
        genomeSegmentLabel += "_to_"+gsegGroup[-1].id

    isFirstSegment = (len(segFiles.variants) == 0)

    segCmd = [ self.params.strelkaGermlineBin ]

    self.appendCommonGenomeSegmentCommandOptions(gsegGroup, genomeSegmentLabel, segCmd)

    segCmd.extend(["--min-mapping-quality",self.params.minMapq])

    segCmd.extend(["--gvcf-output-prefix", self.paths.getTmpSegmentGvcfPrefix(genomeSegmentLabel)])
    segCmd.extend(['--gvcf-min-gqx','15'])
    segCmd.extend(['--gvcf-min-homref-gqx','15'])
    segCmd.extend(['--gvcf-max-snv-strand-bias','10'])
    segCmd.extend(['--enable-read-backed-phasing'])

    segFiles.stats.append(self.paths.getTmpRunStatsPath(genomeSegmentLabel))
    segCmd.extend(["--stats-file", segFiles.stats[-1]])

    if self.params.isRNA:
        segCmd.extend(['--het-variant-frequency-extension', '0.45'])
        segCmd.extend(['--use-rna-scoring'])
        segCmd.extend(['--retain-optimal-soft-clipping'])

    # Empirical Variant Scoring(EVS):
    if self.params.isEVS :
        if self.params.snvScoringModelFile is not None :
            segCmd.extend(['--snv-scoring-model-file', self.params.snvScoringModelFile])
        if self.params.indelScoringModelFile is not None :
            segCmd.extend(['--indel-scoring-model-file', self.params.indelScoringModelFile])

    for bamPath in self.params.bamList :
        segCmd.extend(["--align-file",bamPath])

    if not isFirstSegment :
        segCmd.append("--gvcf-skip-header")
    elif len(self.params.callContinuousVf) > 0 :
        segCmd.extend(["--gvcf-include-header", "VF"])

    if self.params.isHighDepthFilter :
        segCmd.extend(["--chrom-depth-file", self.paths.getChromDepth()])

    if self.params.noCompressBed is not None :
        segCmd.extend(['--nocompress-bed', self.params.noCompressBed])

    if self.params.ploidyFilename is not None :
        segCmd.extend(['--ploidy-region-vcf', self.params.ploidyFilename])

    for gseg in gsegGroup :
        # we have special logic to prevent the continuousVF targets from being grouped, the assertion here
        # verifies that this is working as expected:
        if self.params.callContinuousVf is not None and gseg.chromLabel in self.params.callContinuousVf :
            assert(len(gsegGroup) == 1)
            segCmd.append('--call-continuous-vf')

    if self.params.isEstimateSequenceError :
        for bamIndex in range(len(self.params.bamList)) :
            segCmd.extend(['--indel-error-models-file', self.paths.getIndelErrorModelPath(bamIndex)])
    else :
        segCmd.extend(['--indel-error-models-file', self.params.indelErrorRateDefault])

    segCmd.extend(['--theta-file', self.params.thetaParamFile])

    segTaskLabel=preJoin(taskPrefix,"callGenomeSegment_"+genomeSegmentLabel)
    self.addTask(segTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)

    # clean up and compress genome segment files:
    nextStepWait = set()

    def compressRawVcf(rawVcfFilename, label) :
        """
        Process each raw vcf file with header modifications and bgzip compression
        """

        compressedVariantsPath = rawVcfFilename +".gz"
        compressCmd = "cat "+quote(rawVcfFilename)

        if isFirstSegment :
            def getHeaderFixCmd() :
                cmd  = "\"%s\" -E \"%s\"" % (sys.executable, self.params.vcfCmdlineSwapper)
                cmd += ' "' + " ".join(self.params.configCommandLine) + '"'
                return cmd
            compressCmd += " | " + getHeaderFixCmd()

        compressCmd += " | \"%s\" -c >| \"%s\"" % (self.params.bgzip9Bin, compressedVariantsPath)

        compressTaskLabel=preJoin(taskPrefix,"compressGenomeSegment_"+genomeSegmentLabel+"_"+label)
        self.addTask(compressTaskLabel, compressCmd, dependencies=segTaskLabel, memMb=self.params.callMemMb)
        nextStepWait.add(compressTaskLabel)
        return compressedVariantsPath

    rawVariantsPath = self.paths.getTmpSegmentVariantsPath(genomeSegmentLabel)
    compressedVariantsPath = compressRawVcf(rawVariantsPath, "variants")
    segFiles.variants.append(compressedVariantsPath)

    sampleCount = len(self.params.bamList)
    for sampleIndex in range(sampleCount) :
        rawVariantsPath = self.paths.getTmpSegmentGvcfPath(genomeSegmentLabel, sampleIndex)
        compressedVariantsPath = compressRawVcf(rawVariantsPath, gvcfSampleLabel(sampleIndex))
        segFiles.sample[sampleIndex].gvcf.append(compressedVariantsPath)


    if self.params.isWriteRealignedBam :
        def sortRealignBam(sampleIndex) :
            """
            Sort each realigned bam output
            """
            sortList = segFiles.sample[sampleIndex].bamRealign

            unsorted = self.paths.getTmpUnsortRealignBamPath(genomeSegmentLabel, sampleIndex)
            sorted   = self.paths.getTmpSortRealignBamPath(genomeSegmentLabel, sampleIndex)
            sortList.append(sorted)

            sortCmd="\"%s\" sort \"%s\" -o \"%s\" && rm -f \"%s\"" %\
                    (self.params.samtoolsBin, unsorted, sorted, unsorted)

            sortTaskLabel=preJoin(taskPrefix,"sortRealignedSegment_"+ genomeSegmentLabel + "_" + self.paths.sampleLabel(sampleIndex))
            self.addTask(sortTaskLabel, sortCmd, dependencies=segTaskLabel, memMb=self.params.callMemMb)
            nextStepWait.add(sortTaskLabel)

        for sampleIndex in range(sampleCount) :
            sortRealignBam(sampleIndex)

    return nextStepWait



def callGenome(self,taskPrefix="",dependencies=None):
    """
    run variant caller on all genome segments
    """

    tmpSegmentDir=self.paths.getTmpSegmentDir()
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), getMkdirCmd() + [tmpSegmentDir],
                         dependencies=dependencies, isForceLocal=True)

    segmentTasks = set()
    sampleCount = len(self.params.bamList)

    segFiles = TempVariantCallingSegmentFiles(sampleCount)

    for gsegGroup in self.getStrelkaGenomeSegmentGroupIterator(contigsExcludedFromGrouping = self.params.callContinuousVf) :
        segmentTasks |= callGenomeSegment(self, gsegGroup, segFiles, dependencies=dirTask)

    if len(segmentTasks) == 0 :
        raise Exception("No genome regions to analyze. Possible target region parse error.")

    # create a checkpoint for all segments:
    completeSegmentsTask = self.addTask(preJoin(taskPrefix,"completedAllGenomeSegments"),dependencies=segmentTasks)

    finishTasks = set()

    # merge various VCF outputs
    finishTasks.add(self.concatIndexVcf(taskPrefix, completeSegmentsTask, segFiles.variants,
                                        self.paths.getVariantsOutputPath(), "variants"))
    for sampleIndex in range(sampleCount) :
        concatTask = self.concatIndexVcf(taskPrefix, completeSegmentsTask, segFiles.sample[sampleIndex].gvcf,
                                         self.paths.getGvcfOutputPath(sampleIndex), gvcfSampleLabel(sampleIndex))
        finishTasks.add(concatTask)
        if sampleIndex == 0 :
            outputPath = self.paths.getGvcfOutputPath(sampleIndex)
            outputDirname=os.path.dirname(outputPath)
            outputBasename=os.path.basename(outputPath)
            def linkLegacy(extension) :
                return "ln -s " + quote(outputBasename + extension) + " " + quote(self.paths.getGvcfLegacyFilename() + extension)
            linkCmd = linkLegacy("") + " && " + linkLegacy(".tbi")
            self.addTask(preJoin(taskPrefix, "addLegacyOutputLink"), linkCmd, dependencies=concatTask,
                         isForceLocal=True, cwd=outputDirname)

    # merge segment stats:
    finishTasks.add(self.mergeRunStats(taskPrefix,completeSegmentsTask, segFiles.stats))

    if self.params.isWriteRealignedBam :
        def catRealignedBam(sampleIndex) :
            segmentList = segFiles.sample[sampleIndex].bamRealign
            output = self.paths.getRealignedBamPath(sampleIndex)

            bamCatCmd = bamListCatCmd(self.params.samtoolsBin, segmentList, output)
            bamCatTaskLabel = preJoin(taskPrefix, "realignedBamCat_" + self.paths.sampleLabel(sampleIndex))
            finishTasks.add(self.addTask(bamCatTaskLabel, bamCatCmd, dependencies=completeSegmentsTask))

        for sampleIndex in range(sampleCount) :
            catRealignedBam(sampleIndex)

    if not self.params.isRetainTempFiles :
        rmTmpCmd = getRmdirCmd() + [tmpSegmentDir]
        self.addTask(preJoin(taskPrefix,"removeTmpDir"), rmTmpCmd, dependencies=finishTasks, isForceLocal=True)

    nextStepWait = finishTasks

    return nextStepWait



def validateEstimatedParameters(self, sampleIndex) :
    """
    Check whether the indel error rates are static or estimated values
    """

    jsonFileName = self.paths.getIndelErrorModelPath(sampleIndex)
    if os.path.isfile(jsonFileName) :
        try :
            import json
            jsonFile = open(jsonFileName, 'r')
            indelModel = json.load(jsonFile)
            jsonFile.close()
            if indelModel['sample'][0]['isStatic'] :
                self.flowLog("Adaptive indel error rate estimation for sample '" + indelModel['sample'][0]['sampleName'] + "' did not succeed, using static indel error model instead.", logState=LogState.WARNING)

        except ImportError :
            self.flowLog("Python version does not support json parsing. Skipping json validation.", logState=pyflowDir.LogState.WARNING)

        except :
            self.flowLog("Error while parsing the indel rate model file '"+jsonFileName+"'. The file format may be invalid", logState=LogState.ERROR)
    else :
        self.flowLog("Error while parsing the indel rate model file '"+jsonFileName+"'. File does not exist.", logState=LogState.ERROR)


class ValidateEstimatedParametersWorkflow(WorkflowRunner) :
    """
    A separate workflow is setup around validateEstimatedParameters() so that the workflow execution can be
    delayed until the json files exist
    """

    def __init__(self, params, paths) :
        self.paths = paths
        self.params = params

    def workflow(self) :
        sampleCount = len(self.params.bamList)

        for sampleIndex in range(sampleCount) :
            validateEstimatedParameters(self,sampleIndex)


class CallWorkflow(StrelkaSharedCallWorkflow) :
    """
    TODO:
    callGenome() is wrapped by a workflow for legacy reasons, this could be reduced to a standard task tree in the
    parent workflow when there is time to impl/test.
    """

    def __init__(self, params, paths) :
        super(CallWorkflow,self).__init__(params)
        self.paths = paths

    def workflow(self) :
        callGenome(self)



class PathInfo(SharedPathInfo):
    """
    Object to centralize shared workflow path names
    """

    def __init__(self, params) :
        super(PathInfo,self).__init__(params)

    def getTmpSegmentGvcfPrefix(self, genomeSegmentLabel) :
        return os.path.join(self.getTmpSegmentDir(), "segment.%s." % (genomeSegmentLabel))

    def getTmpSegmentVariantsPath(self, genomeSegmentLabel) :
        return self.getTmpSegmentGvcfPrefix(genomeSegmentLabel) + "variants.vcf"

    def sampleLabel(self, sampleIndex):
        return "S%i" % (sampleIndex+1)

    def getTmpSegmentGvcfPath(self, genomeSegmentLabel, sampleIndex) :
        return self.getTmpSegmentGvcfPrefix(genomeSegmentLabel) + "genome.%s.vcf" % (self.sampleLabel(sampleIndex))

    def getTmpUnsortRealignBamPath(self, genomeSegmentLabel, sampleIndex) :
        return self.getTmpUnsortRealignBamPrefix(genomeSegmentLabel) + "%s.bam" % (self.sampleLabel(sampleIndex))

    def getTmpSortRealignBamPath(self, genomeSegmentLabel, sampleIndex) :
        return os.path.join(self.getTmpSegmentDir(), "segment.%s.sorted.realigned.%s.bam" % (genomeSegmentLabel, self.sampleLabel(sampleIndex)))

    def getVariantsOutputPath(self) :
        return os.path.join( self.params.variantsDir, "variants.vcf.gz")

    def getGvcfOutputPath(self, sampleIndex) :
        return os.path.join( self.params.variantsDir, "genome.%s.vcf.gz" % (self.sampleLabel(sampleIndex)))

    def getGvcfLegacyFilename(self) :
        return "genome.vcf.gz"

    def getRealignedBamPath(self, sampleIndex) :
        return os.path.join( self.params.realignedDir, 'realigned.%s.bam' % (self.sampleLabel(sampleIndex)))

    def getTmpSegmentNonemptySiteCountsPath(self, sampleIndex, genomeSegmentLabel) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join(self.getTmpErrorEstimationDir(), "nonEmptySiteCounts.Sample%s.%s.tsv" % (sampleIndexStr, genomeSegmentLabel))

    def getTmpSegmentAlleleCountsPath(self, sampleIndex, genomeSegmentLabel) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join(self.getTmpErrorEstimationDir(), "sequenceAlleleCounts.Sample%s.%s.bin" % (sampleIndexStr, genomeSegmentLabel))

    def getAlleleCountsOutputPath(self, sampleIndex) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join( self.params.workDir, "sequenceAlleleCounts.Sample%s.bin" % (sampleIndexStr))

    def getIndelErrorModelPath(self, sampleIndex) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join( self.params.workDir, "sequenceErrorModel.Sample%s.json" % (sampleIndexStr))



class StrelkaGermlineWorkflow(StrelkaSharedWorkflow) :
    """
    Germline small variant calling workflow
    """

    def __init__(self,params) :
        global PathInfo
        super(StrelkaGermlineWorkflow,self).__init__(params,PathInfo)

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

        msg  = "Strelka germline workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg


    def workflow(self) :
        self.flowLog("Initiating Strelka germline workflow version: %s" % (__version__))
        self.setCallMemMb()

        callPreReqs = set()
        estimatePreReqs = set()
        if self.params.isHighDepthFilter :
            estimatePreReqs |= strelkaGermlineGetDepthFromAlignments(self)

        if self.params.isEstimateSequenceError :
            validatePreReq = self.addWorkflowTask("EstimateSeqErrorParams", EstimateSequenceErrorWorkflow(self.params, self.paths), dependencies=estimatePreReqs)
            callPreReqs.add(self.addWorkflowTask("ValidateSeqErrorParams", ValidateEstimatedParametersWorkflow(self.params, self.paths), dependencies=validatePreReq))
        else :
            callPreReqs = estimatePreReqs

        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)

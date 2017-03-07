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

from configBuildTimeInfo import workflowVersion
from configureUtil import safeSetBool, joinFile
from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getRmdirCmd, runDepthFromAlignments
from strelkaSharedWorkflow import runCount, SharedPathInfo, \
                           StrelkaSharedCallWorkflow, StrelkaSharedWorkflow
from workflowUtil import ensureDir, preJoin, bamListCatCmd


__version__ = workflowVersion



def strelkaGermlineRunDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    bamList=[]
    if len(self.params.bamList) :
        bamList = self.params.bamList
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

class TempSegmentFilesPerSample :
    def __init__(self) :
        self.gvcf = []


class TempSegmentFiles :
    def __init__(self, sampleCount) :
        self.variants = []
        self.bamRealign = []
        self.stats = []
        self.sample = [TempSegmentFilesPerSample() for _ in range(sampleCount)]


# we need extra quoting for files with spaces in this workflow because some commands are stringified as shell calls:
def quote(instr):
    return "\"%s\"" % (instr)


def gvcfSampleLabel(sampleIndex) :
    return "gVCF_S%i" % (sampleIndex+1)


def callGenomeSegment(self, gsegGroup, segFiles, taskPrefix="", dependencies=None) :

    assert(len(gsegGroup) != 0)
    gid=gsegGroup[0].id
    if len(gsegGroup) > 1 :
        gid += "_to_"+gsegGroup[-1].id

    isFirstSegment = (len(segFiles.variants) == 0)

    segCmd = [ self.params.strelkaGermlineBin ]

    self.appendCommonGenomeSegmentCommandOptions(gsegGroup, segCmd)

    segCmd.extend(["-min-mapping-quality",self.params.minMapq])
    segCmd.extend(["-max-window-mismatch", "2", "20" ])

    segCmd.extend(["--gvcf-output-prefix", self.paths.getTmpSegmentGvcfPrefix(gid)])
    segCmd.extend(['--gvcf-min-gqx','15'])
    segCmd.extend(['--gvcf-min-homref-gqx','15'])
    segCmd.extend(['--gvcf-max-snv-strand-bias','10'])
    segCmd.extend(['-min-qscore','17'])
    segCmd.extend(['-bsnp-ssd-no-mismatch', '0.35'])
    segCmd.extend(['-bsnp-ssd-one-mismatch', '0.6'])
    segCmd.extend(['-min-vexp', '0.25'])
    segCmd.extend(['--enable-read-backed-phasing'])

    segFiles.stats.append(self.paths.getTmpRunStatsPath(gid))
    segCmd.extend(["--stats-file", segFiles.stats[-1]])

    if self.params.isRNA:
        # RNA-Seq het calls are considered over a wider frequency range:
        segCmd.extend(['-bsnp-diploid-het-bias', '0.40'])
        segCmd.extend(['--use-rna-scoring'])

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

    # TODO STREL-125 come up with new solution for outbams
    if self.params.isWriteRealignedBam :
        segCmd.extend(["-realigned-read-file", self.paths.getTmpUnsortRealignBamPath(gid)])

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

    segTaskLabel=preJoin(taskPrefix,"callGenomeSegment_"+gid)
    self.addTask(segTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)

    # clean up and compress genome segment files:
    nextStepWait = set()

    def compressRawVcf(rawVcfFilename, label) :
        """
        process each raw vcf file with header modifications and bgzip compression
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

        compressTaskLabel=preJoin(taskPrefix,"compressGenomeSegment_"+gid+"_"+label)
        self.addTask(compressTaskLabel, compressCmd, dependencies=segTaskLabel, memMb=self.params.callMemMb)
        nextStepWait.add(compressTaskLabel)
        return compressedVariantsPath

    rawVariantsPath = self.paths.getTmpSegmentVariantsPath(gid)
    compressedVariantsPath = compressRawVcf(rawVariantsPath, "variants")
    segFiles.variants.append(compressedVariantsPath)

    sampleCount = len(self.params.bamList)
    for sampleIndex in range(sampleCount) :
        rawVariantsPath = self.paths.getTmpSegmentGvcfPath(gid, sampleIndex)
        compressedVariantsPath = compressRawVcf(rawVariantsPath, gvcfSampleLabel(sampleIndex))
        segFiles.sample[sampleIndex].gvcf.append(compressedVariantsPath)


    if self.params.isWriteRealignedBam :
        def sortRealignBam(sortList) :
            unsorted = self.paths.getTmpUnsortRealignBamPath(gid)
            sorted   = self.paths.getTmpRealignBamPath(gid)
            sortList.append(sorted)

            # adjust sorted to remove the ".bam" suffix
            sorted = sorted[:-4]
            sortCmd="\"%s\" sort \"%s\" \"%s\" && rm -f \"%s\"" % (self.params.samtoolsBin,unsorted,sorted,unsorted)

            sortTaskLabel=preJoin(taskPrefix,"sortRealignedSegment_"+gid)
            self.addTask(sortTaskLabel,sortCmd,dependencies=segTaskLabel,memMb=self.params.callMemMb)
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

    sampleCount = len(self.params.bamList)

    segFiles = TempSegmentFiles(sampleCount)

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
        def finishBam(tmpList, output, label) :
            cmd = bamListCatCmd(self.params.samtoolsBin, tmpList, output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))

        finishBam(segFiles.bamRealign, self.paths.getRealignedBamPath(), "realigned")

    if not self.params.isRetainTempFiles :
        rmStatsTmpCmd = getRmdirCmd() + [tmpSegmentDir]
        rmTask=self.addTask(preJoin(taskPrefix,"rmTmpDir"),rmStatsTmpCmd,dependencies=finishTasks, isForceLocal=True)

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

    def getTmpSegmentGvcfPrefix(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "segment.%s." % (segStr))

    def getTmpSegmentVariantsPath(self, segStr) :
        return self.getTmpSegmentGvcfPrefix(segStr) + "variants.vcf"

    def getTmpSegmentGvcfPath(self, segStr, sampleIndex) :
        return self.getTmpSegmentGvcfPrefix(segStr) + "genome.S%i.vcf" % (sampleIndex+1)

    def getTmpUnsortRealignBamPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "%s.unsorted.realigned.bam" % (segStr))

    def getTmpRealignBamPath(self, segStr,) :
        return os.path.join( self.getTmpSegmentDir(), "%s.realigned.bam" % (segStr))

    def getVariantsOutputPath(self) :
        return os.path.join( self.params.variantsDir, "variants.vcf.gz")

    def getGvcfOutputPath(self, sampleIndex) :
        return os.path.join( self.params.variantsDir, "genome.S%i.vcf.gz" % (sampleIndex+1))

    def getGvcfLegacyFilename(self) :
        return "genome.vcf.gz"

    def getRealignedBamPath(self) :
        return os.path.join( self.params.realignedDir, 'realigned.bam')



class StrelkaGermlineWorkflow(StrelkaSharedWorkflow) :
    """
    germline small variant calling workflow
    """

    def __init__(self,params,iniSections) :
        global PathInfo
        super(StrelkaGermlineWorkflow,self).__init__(params,iniSections,PathInfo)

        # format bam lists:
        if self.params.bamList is None : self.params.bamList = []

        # format other:
        safeSetBool(self.params,"isWriteRealignedBam")

        if self.params.isWriteRealignedBam :
            self.params.realignedDir=os.path.join(self.params.resultsDir,"realigned")
            ensureDir(self.params.realignedDir)

        if self.params.isExome :
            self.params.isEVS = False

        if self.params.isRNA :
            self.params.snvScoringModelFile = joinFile(self.params.configDir,'RNAVariantScoringModels.json')
            self.params.indelScoringModelFile = None


    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Strelka germline workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg


    def workflow(self) :
        self.flowLog("Initiating Strelka germline workflow version: %s" % (__version__))
        self.setCallMemMb()

        callPreReqs = set()
        callPreReqs |= runCount(self)
        if self.params.isHighDepthFilter :
            callPreReqs |= strelkaGermlineRunDepthFromAlignments(self)
        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)

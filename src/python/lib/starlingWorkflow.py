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
Starling germline small variant calling workflow
"""


import os.path
import shutil
import sys

# add script path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(scriptDir))

# add pyflow path:
# TODO: get a more robust link to the pyflow dir at config time:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from pyflow import WorkflowRunner
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getNextGenomeSegment, getFastaChromOrderSize, bamListCatCmd

from configureUtil import argToBool, getIniSections, dumpIniSections



def getVersion() :
    return "@STARKA_VERSION@"


__version__ = getVersion()



def runCount(self, taskPrefix="", dependencies=None) :
    """
    count size of fasta chromosomes
    """
    cmd  = "%s %s > %s"  % (self.params.countFastaBin, self.params.referenceFasta, self.paths.getRefCountFile())

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"RefCount"), cmd, dependencies=dependencies))

    return nextStepWait


def runDepth(self,taskPrefix="",dependencies=None) :
    """
    estimate chrom depth
    """

    bamFile=""
    if len(self.params.bamList) :
        bamFile = self.params.bamList[0]
    else :
        return set()


    cmd  = "%s -E %s" % (sys.executable, self.params.getChromDepth)
    cmd += " --bam '%s'" % (bamFile)
    cmd += " > %s" % (self.paths.getChromDepth())

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"estimateChromDepth"),cmd,dependencies=dependencies))

    return nextStepWait

def runIndelModel(self,taskPrefix="",dependencies=None) :
    """
    estimate indel error paramters
    """

    bamFile=""
    if len(self.params.bamList) :
        bamFile = self.params.bamList[0]
    else :
        return set()

    nextStepWait = set()
    nextStepWait.add(self.addWorkflowTask("GenerateIndelModel", indelErrorWorkflow(self.params, self.paths), dependencies=dependencies))

    return nextStepWait

class TempSegmentFiles :
    def __init__(self) :
        self.gvcf = []
        self.bamRealign = []



def callGenomeSegment(self, gseg, segFiles, taskPrefix="", dependencies=None) :

    isFirstSegment = (len(segFiles.gvcf) == 0)

    segStr = str(gseg.id)

    segCmd = [ self.params.starlingBin ]

    segCmd.append("-clobber")
    segCmd.extend(["-min-paired-align-score",self.params.minMapq])
    segCmd.extend(["-min-single-align-score",self.params.minMapq])
    segCmd.extend(["-bam-seq-name", gseg.chromLabel] )
    segCmd.extend(["-report-range-begin", str(gseg.beginPos) ])
    segCmd.extend(["-report-range-end", str(gseg.endPos) ])
    segCmd.extend(["-samtools-reference", self.params.referenceFasta ])
    segCmd.extend(["-max-window-mismatch", "2", "20" ])
    segCmd.extend(["-genome-size", str(self.params.knownSize)] )
    segCmd.extend(["-max-indel-size", "50"] )

    segCmd.extend(["--gvcf-file","-"])
    segCmd.extend(['--gvcf-min-gqx','30'])
    segCmd.extend(['--gvcf-max-snv-strand-bias','10'])
    segCmd.extend(['--gvcf-max-indel-ref-repeat', '-1'])
    segCmd.extend(['-min-qscore','17'])
    segCmd.extend(['-bsnp-ssd-no-mismatch', '0.35'])
    segCmd.extend(['-bsnp-ssd-one-mismatch', '0.6'])
    segCmd.extend(['-min-vexp', '0.25'])
    segCmd.extend(['--calibration-model-file',self.params.vqsrModel])
#    segCmd.extend(['--scoring-models', self.params.scoringModelFile])

    for bamPath in self.params.bamList :
        segCmd.extend(["-bam-file",bamPath])

    segCmd.extend(["--report-file", self.paths.getTmpSegmentReportPath(gseg.pyflowId)])

    if not isFirstSegment :
        segCmd.append("--gvcf-skip-header")

    if not self.params.isSkipDepthFilters :
        segCmd.extend(["--chrom-depth-file", self.paths.getChromDepth()])

    if self.params.isWriteRealignedBam :
        segCmd.extend(["-realigned-read-file", self.paths.getTmpUnsortRealignBamPath(segStr)])

    if self.params.indelCandidates is not None :
        segCmd.extend(['--candidate-indel-input-vcf', self.params.indelCandidates])

    if self.params.forcedGTIndels is not None :
        segCmd.extend(['--force-output-vcf', self.params.forcedGTIndels])

    if self.params.minorAllele is not None :
        segCmd.extend(['--minor-allele-bed-file', self.params.minorAllele])

    if self.params.extraStarlingArguments is not None :
        for arg in self.params.extraStarlingArguments.strip().split() :
            segCmd.append(arg)

     # gvcf is written to stdout so we need shell features:
    segCmd = " ".join(segCmd)

    segFiles.gvcf.append(self.paths.getTmpSegmentGvcfPath(segStr))
    segCmd += " | %s -c >| %s" % (self.params.bgzip9Bin, segFiles.gvcf[-1])

    nextStepWait = set()

    setTaskLabel=preJoin(taskPrefix,"callGenomeSegment_"+gseg.pyflowId)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)
    nextStepWait.add(setTaskLabel)

    if self.params.isWriteRealignedBam :
        def sortRealignBam(sortList) :
            unsorted = self.paths.getTmpUnsortRealignBamPath(segStr)
            sorted   = self.paths.getTmpRealignBamPath(segStr)
            sortList.append(sorted)

            # adjust sorted to remove the ".bam" suffix
            sorted = sorted[:-4]
            sortCmd="%s sort %s %s && rm -f %s" % (self.params.samtoolsBin,unsorted,sorted,unsorted)

            sortTaskLabel=preJoin(taskPrefix,"sortRealignedSegment_"+gseg.pyflowId)
            self.addTask(sortTaskLabel,sortCmd,dependencies=setTaskLabel,memMb=self.params.callMemMb)
            nextStepWait.add(sortTaskLabel)

        sortRealignBam(segFiles.bamRealign)

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
        finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeGVCF"), catCmd, dependencies=completeSegmentsTask))

    finishVcf(segFiles.gvcf, self.paths.getGvcfOutputPath(),"gVCF")

    if self.params.isWriteRealignedBam :
        def finishBam(tmpList, output, label) :
            cmd = bamListCatCmd(self.params.samtoolsBin, tmpList, output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))

        finishBam(segFiles.bamRealign, self.paths.getRealignedBamPath(), "realigned")

    cleanTask=self.addTask(preJoin(taskPrefix,"cleanTmpDir"), "rm -rf "+tmpGraphDir, dependencies=finishTasks, isForceLocal=True)

    nextStepWait = finishTasks

    return nextStepWait



"""
A separate call workflow is setup so that we can delay the workflow execution until
the ref count file exists
"""
class CallWorkflow(WorkflowRunner) :

    def __init__(self,params,paths) :
        self.params = params
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



class PathInfo:
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        self.params = params

    def getChromDepth(self) :
        return os.path.join(self.params.workDir,"chromDepth.txt")

    def getTmpSegmentDir(self) :
        return os.path.join(self.params.workDir, "genomeSegment.tmpdir")

    def getTmpSegmentGvcfPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "genome.%s.vcf.gz" % (segStr))

    def getTmpUnsortRealignBamPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "%s.unsorted.realigned.bam" % (segStr))

    def getTmpRealignBamPath(self, segStr,) :
        return os.path.join( self.getTmpSegmentDir(), "%s.realigned.bam" % (segStr))

    def getTmpSegmentReportPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "stats.%s.txt" % (segStr))

    def getVariantsDir(self) :
        return self.params.variantsDir

    def getGvcfOutputPath(self) :
        return os.path.join( self.getVariantsDir(), "genome.vcf.gz")

    def getRealignedBamPath(self) :
        return os.path.join( self.params.realignedDir, 'realigned.bam');

    def getRefCountFile(self) :
        return os.path.join( self.params.workDir, "refCount.txt")



class StarlingWorkflow(WorkflowRunner) :
    """
    germline small variant calling workflow
    """

    def __init__(self,params,iniSections) :

        # clear out some potentially destabilizing env variables:
        clearList = [ "PYTHONPATH", "PYTHONHOME"]
        for key in clearList :
            if key in os.environ :
                del os.environ[key]

        self.params=params
        self.iniSections=iniSections

        # format bam lists:
        if self.params.bamList is None : self.params.bamList = []

        # format other:
        self.params.isWriteRealignedBam = argToBool(self.params.isWriteRealignedBam)
        self.params.isSkipDepthFilters = argToBool(self.params.isSkipDepthFilters)
        self.params.isSkipIndelErrorModel = argToBool(self.params.isSkipIndelErrorModel)

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

        if self.params.isWriteRealignedBam :
            self.params.realignedDir=os.path.join(self.params.resultsDir,"realigned")
            ensureDir(self.params.realignedDir)

        indexRefFasta=self.params.referenceFasta+".fai"

        if self.params.referenceFasta is None:
            raise Exception("No reference fasta defined.")
        else:
            checkFile(self.params.referenceFasta,"reference fasta")
            checkFile(indexRefFasta,"reference fasta index")

        # read fasta index
        (self.params.chromOrder,self.params.chromSizes) = getFastaChromOrderSize(indexRefFasta)

        # sanity check some parameter typing:
        MEGABASE = 1000000
        self.params.scanSize = int(self.params.scanSizeMb) * MEGABASE

        self.paths = PathInfo(self.params)


    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Starling workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Starling workflow version: %s" % (__version__))

        callPreReqs = set()
        callPreReqs |= runCount(self)
        if not self.params.isSkipDepthFilters :
            callPreReqs |= runDepth(self)
        if not self.params.isSkipIndelErrorModel :
            callPreReqs |= runIndelModel(self)      
        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)


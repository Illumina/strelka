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

# add this path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(scriptDir)

# add pyflow path:
sys.path.append(os.path.join(scriptDir,"pyflow"))


from pyflow import WorkflowRunner
from starkaWorkflow import StarkaWorkflow
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getNextGenomeSegment, bamListCatCmd

from configureUtil import safeSetBool, getIniSections, dumpIniSections



def getVersion() :
    return "@STARKA_VERSION@"


__version__ = getVersion()



def runCount(self, taskPrefix="", dependencies=None) :
    """
    count size of fasta chromosomes
    """
    cmd  = "%s '%s' > %s"  % (self.params.countFastaBin, self.params.referenceFasta, self.paths.getRefCountFile())

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

    # we need extra quoting for files with spaces in this workflow because command is stringified below to enable gVCF pipe:
    def quote(instr):
        return "'%s'" % (instr)

    segCmd = [ self.params.starlingBin ]

    segCmd.append("-clobber")
    segCmd.extend(["-min-paired-align-score",self.params.minMapq])
    segCmd.extend(["-min-single-align-score",self.params.minMapq])
    segCmd.extend(["-bam-seq-name", gseg.chromLabel] )
    segCmd.extend(["-report-range-begin", str(gseg.beginPos) ])
    segCmd.extend(["-report-range-end", str(gseg.endPos) ])
    segCmd.extend(["-samtools-reference", quote(self.params.referenceFasta) ])
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
    segCmd.extend(['--calibration-model-file',self.params.vqsrModelFile])
    segCmd.extend(['--scoring-model',self.params.vqsrModel ])
    segCmd.extend(['--indel-ref-error-factor',self.params.indelRefErrorFactor])

    if self.params.isReportVQSRMetrics :
        segCmd.append("--gvcf-report-VQSRmetrics")

    segCmd.extend(['--do-short-range-phasing'])

    for bamPath in self.params.bamList :
        segCmd.extend(["-bam-file",quote(bamPath)])

    segCmd.extend(["--report-file", self.paths.getTmpSegmentReportPath(gseg.pyflowId)])

    if not isFirstSegment :
        segCmd.append("--gvcf-skip-header")

    if self.params.isHighDepthFilter :
        segCmd.extend(["--chrom-depth-file", self.paths.getChromDepth()])

    if self.params.isWriteRealignedBam :
        segCmd.extend(["-realigned-read-file", self.paths.getTmpUnsortRealignBamPath(segStr)])

    def addListCmdOption(optList,arg) :
        if optList is None : return
        for val in optList :
            segCmd.extend([arg, val])

    addListCmdOption(self.params.indelCandidatesList, '--candidate-indel-input-vcf')
    addListCmdOption(self.params.forcedGTList, '--force-output-vcf')

    if self.params.noCompressBed is not None :
        segCmd.extend(['--nocompress-bed', self.params.noCompressBed])

    if self.params.ploidyBed is not None :
        segCmd.extend(['--ploidy-region-bed', self.params.ploidyBed])

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

    if len(graphTasks) == 0 :
        raise Exception("No genome regions to analyze. Possible target region parse error.")

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



class StarlingWorkflow(StarkaWorkflow) :
    """
    germline small variant calling workflow
    """

    def __init__(self,params,iniSections) :

        super(StarlingWorkflow,self).__init__(params,iniSections)

        # format bam lists:
        if self.params.bamList is None : self.params.bamList = []

        # format other:
        safeSetBool(self.params,"isWriteRealignedBam")
        safeSetBool(self.params,"isSkipIndelErrorModel")

        if self.params.isWriteRealignedBam :
            self.params.realignedDir=os.path.join(self.params.resultsDir,"realigned")
            ensureDir(self.params.realignedDir)

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
        if self.params.isHighDepthFilter :
            callPreReqs |= runDepth(self)
        if not self.params.isSkipIndelErrorModel :
            callPreReqs |= runIndelModel(self)
        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)


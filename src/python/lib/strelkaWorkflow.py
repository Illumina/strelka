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
Strelka somatic small variant calling workflow
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
    cmd  = "%s '%s' > %s"  % (self.params.countFastaBin, self.params.referenceFasta, self.paths.getRefCountFile())

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"RefCount"), cmd, dependencies=dependencies))

    return nextStepWait



def runDepth(self,taskPrefix="",dependencies=None) :
    """
    estimate chrom depth
    """

    bamFile=""
    if len(self.params.normalBamList) :
        bamFile = self.params.normalBamList[0]
    elif len(self.params.tumorBamList) :
        bamFile = self.params.tumorBamList[0]
    else :
        return set()


    cmd  = "%s -E %s" % (sys.executable, self.params.getChromDepth)
    cmd += " --bam '%s'" % (bamFile)
    cmd += " > %s" % (self.paths.getChromDepth())

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"estimateChromDepth"),cmd,dependencies=dependencies))

    return nextStepWait



class TempSegmentFiles :
    def __init__(self) :
        self.snv = []
        self.indel = []
        self.callable = []
        self.normalRealign = []
        self.tumorRealign = []



def callGenomeSegment(self, gseg, segFiles, taskPrefix="", dependencies=None) :

    isFirstSegment = (len(segFiles.snv) == 0)

    segStr = str(gseg.id)

    segCmd = [ self.params.strelkaBin ]
    segCmd.append("-clobber")
    segCmd.append("-filter-unanchored")
    segCmd.extend(["-min-paired-align-score",str(self.params.minTier1Mapq)])
    segCmd.extend(["-min-single-align-score","10"])
    segCmd.extend(["-min-qscore","0"])
    segCmd.extend(["-report-range-begin", str(gseg.beginPos) ])
    segCmd.extend(["-report-range-end", str(gseg.endPos) ])
    segCmd.extend(["-samtools-reference", self.params.referenceFasta ])
    segCmd.extend(["-max-window-mismatch", "3", "20" ])
    segCmd.append("-print-used-allele-counts")
    segCmd.extend(["-bam-seq-name", gseg.chromLabel] )
    segCmd.extend(["-genome-size", str(self.params.knownSize)] )
    segCmd.extend(["-max-indel-size", "50"] )
    segCmd.extend(["-indel-nonsite-match-prob", "0.5"] )
    segCmd.extend(["--somatic-snv-rate", str(self.params.ssnvPrior) ] )
    segCmd.extend(["--shared-site-error-rate", str(self.params.ssnvNoise) ] )
    segCmd.extend(["--shared-site-error-strand-bias-fraction", str(self.params.ssnvNoiseStrandBiasFrac) ] )
    segCmd.extend(["--somatic-indel-rate", str(self.params.sindelPrior) ] )
    segCmd.extend(["--shared-indel-error-factor", str(self.params.sindelNoiseFactor)])
    segCmd.extend(["--tier2-min-single-align-score", str(self.params.minTier2Mapq) ] )
    segCmd.extend(["--tier2-min-paired-align-score", str(self.params.minTier2Mapq) ] )
    segCmd.append("--tier2-single-align-score-rescue-mode")
    segCmd.extend(["--tier2-mismatch-density-filter-count", "10"] )
    segCmd.append("--tier2-no-filter-unanchored")
    segCmd.extend(["--tier2-indel-nonsite-match-prob", "0.25"] )
    segCmd.append("--tier2-include-singleton")
    segCmd.append("--tier2-include-anomalous")

    segCmd.extend(["--strelka-snv-max-filtered-basecall-frac", str(self.params.snvMaxFilteredBasecallFrac)])
    segCmd.extend(["--strelka-snv-max-spanning-deletion-frac", str(self.params.snvMaxSpanningDeletionFrac)])
    segCmd.extend(["--strelka-snv-min-qss-ref", str(self.params.ssnvQuality_LowerBound)])

    segCmd.extend(["--strelka-indel-max-ref-repeat", str(self.params.indelMaxRefRepeat)])
    segCmd.extend(["--strelka-indel-max-int-hpol-length", str(self.params.indelMaxIntHpolLength)])
    segCmd.extend(["--strelka-indel-max-window-filtered-basecall-frac", str(self.params.indelMaxWindowFilteredBasecallFrac)])
    segCmd.extend(["--strelka-indel-min-qsi-ref", str(self.params.sindelQuality_LowerBound)])

    # do not apply VQSR in exome case
    if not self.params.isExome :
        segCmd.extend(['--indel-scoring-models', self.params.scoringModelFile])

    for bamPath in self.params.normalBamList :
        segCmd.extend(["-bam-file", bamPath])
    for bamPath in self.params.tumorBamList :
        segCmd.extend(["--tumor-bam-file", bamPath])

    tmpSnvPath = self.paths.getTmpSegmentSnvPath(segStr)
    segFiles.snv.append(tmpSnvPath)
    segCmd.extend(["--somatic-snv-file ", tmpSnvPath ] )

    tmpIndelPath = self.paths.getTmpSegmentIndelPath(segStr)
    segFiles.indel.append(tmpIndelPath)
    segCmd.extend(["--somatic-indel-file", tmpIndelPath ] )

    if self.params.isWriteCallableRegion :
        tmpCallablePath = self.paths.getTmpSegmentRegionPath(segStr)
        segFiles.callable.append(tmpCallablePath)
        segCmd.extend(["--somatic-callable-region-file", tmpCallablePath ])

    if self.params.isWriteRealignedBam :
        segCmd.extend(["-realigned-read-file", self.paths.getTmpUnsortRealignBamPath(segStr, "normal")])
        segCmd.extend(["--tumor-realigned-read-file",self.paths.getTmpUnsortRealignBamPath(segStr, "tumor")])

    def addListCmdOption(optList,arg) :
        if optList is None : return
        for val in optList :
            segCmd.extend([arg, val])

    addListCmdOption(self.params.indelCandidatesList, '--candidate-indel-input-vcf')
    addListCmdOption(self.params.forcedGTList, '--force-output-vcf')
    addListCmdOption(self.params.noiseVcfList, '--noise-vcf')

    if self.params.extraStrelkaArguments is not None :
        for arg in self.params.extraStrelkaArguments.strip().split() :
            segCmd.append(arg)

    segCmd.extend(["--report-file", self.paths.getTmpSegmentReportPath(gseg.pyflowId)])

    if not isFirstSegment :
        segCmd.append("--strelka-skip-header")

    if self.params.isHighDepthFilter :
        segCmd.extend(["--strelka-chrom-depth-file", self.paths.getChromDepth()])
        segCmd.extend(["--strelka-max-depth-factor", self.params.depthFilterMultiple])

    nextStepWait = set()

    callTask=preJoin(taskPrefix,"callGenomeSegment_"+gseg.pyflowId)
    self.addTask(callTask,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)

    compressLabel=preJoin(taskPrefix,"compressSegmentOutput_"+gseg.pyflowId)
    compressCmd="%s %s && %s %s" % (self.params.bgzipBin, tmpSnvPath, self.params.bgzipBin, tmpIndelPath)
    if self.params.isWriteCallableRegion :
        compressCmd += " && %s %s" % (self.params.bgzipBin, self.paths.getTmpSegmentRegionPath(segStr))

    self.addTask(compressLabel, compressCmd, dependencies=callTask, isForceLocal=True)
    nextStepWait.add(compressLabel)

    if self.params.isWriteRealignedBam :
        def sortRealignBam(label, sortList) :
            unsorted = self.paths.getTmpUnsortRealignBamPath(segStr, label)
            sorted   = self.paths.getTmpRealignBamPath(segStr, label)
            sortList.append(sorted)

            # adjust sorted to remove the ".bam" suffix
            sorted = sorted[:-4]
            sortCmd="%s sort %s %s && rm -f %s" % (self.params.samtoolsBin,unsorted,sorted,unsorted)

            sortTaskLabel=preJoin(taskPrefix,"sortRealignedSegment_"+label+"_"+gseg.pyflowId)
            self.addTask(sortTaskLabel,sortCmd,dependencies=callTask,memMb=self.params.callMemMb)
            nextStepWait.add(sortTaskLabel)

        sortRealignBam("normal", segFiles.normalRealign)
        sortRealignBam("tumor", segFiles.tumorRealign)

    return nextStepWait



def callGenome(self,taskPrefix="",dependencies=None):
    """
    run strelka on all genome segments
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

    def finishTabixIndexedFile(tmpList, output, label, fileType) :
        assert(len(tmpList) > 0)

        if len(tmpList) > 1 :
            catCmd = [self.params.bgcatBin,"-o",output]
            catCmd.extend(tmpList)
        else :
            catCmd = "mv -f %s %s" % (tmpList[0],output)

        indexCmd = "%s -p %s %s" % (self.params.tabixBin, fileType, output)
        catTask = self.addTask(preJoin(taskPrefix,label+"_concat_"+fileType), catCmd,
                               dependencies=completeSegmentsTask, isForceLocal=True)
        finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_index_"+fileType), indexCmd,
                                     dependencies=catTask, isForceLocal=True))

    def finishVcf(inputList, output, label) :
        assert(len(inputList) > 0)
        tmpList = [file + ".gz" for file in inputList]
        finishTabixIndexedFile(tmpList, output, label, "vcf")

    def finishBed(inputList, output, label) :
        assert(len(inputList) > 0)
        tmpList = [file + ".gz" for file in inputList]
        finishTabixIndexedFile(tmpList, output, label, "bed")


    finishVcf(segFiles.snv, self.paths.getSnvOutputPath(),"SNV")
    finishVcf(segFiles.indel, self.paths.getIndelOutputPath(),"Indel")

    if self.params.isWriteCallableRegion :
        finishBed(segFiles.callable, self.paths.getRegionOutputPath(), "callableRegions")

    if self.params.isWriteRealignedBam :
        def finishBam(tmpList, output, label) :
            cmd = bamListCatCmd(self.params.samtoolsBin,tmpList,output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))

        finishBam(segFiles.normalRealign, self.paths.getRealignedBamPath("normal"), "realignedNormal")
        finishBam(segFiles.tumorRealign, self.paths.getRealignedBamPath("tumor"), "realignedTumor")

    # add a tmp folder rm step here....

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

    def getTmpSegmentSnvPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "somatic.snvs.unfiltered.%s.vcf" % (segStr))

    def getTmpSegmentIndelPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "somatic.indels.unfiltered.%s.vcf" % (segStr))

    def getTmpSegmentRegionPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "somatic.callable.region.%s.bed" % (segStr))

    def getTmpUnsortRealignBamPath(self, segStr, label) :
        return os.path.join( self.getTmpSegmentDir(), "%s.%s.unsorted.realigned.bam" % (label, segStr))

    def getTmpRealignBamPath(self, segStr, label) :
        return os.path.join( self.getTmpSegmentDir(), "%s.%s.realigned.bam" % (label, segStr))

    def getTmpSegmentReportPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "stats.%s.txt" % (segStr))

    def getVariantsDir(self) :
        return self.params.variantsDir

    def getSnvOutputPath(self) :
        return os.path.join( self.getVariantsDir(), "somatic.snvs.vcf.gz")

    def getIndelOutputPath(self) :
        return os.path.join( self.getVariantsDir(), "somatic.indels.vcf.gz")

    def getRegionOutputPath(self) :
        return os.path.join( self.params.regionsDir, 'somatic.callable.region.bed.gz');

    def getRealignedBamPath(self, label) :
        return os.path.join( self.params.realignedDir, '%s.realigned.bam' % (label));

    def getRefCountFile(self) :
        return os.path.join( self.params.workDir, "refCount.txt")



class StrelkaWorkflow(StarkaWorkflow) :
    """
    Strelka somatic small variant calling workflow
    """

    def __init__(self,params,iniSections) :

        super(StrelkaWorkflow,self).__init__(params,iniSections)

        # format bam lists:
        if self.params.normalBamList is None : self.params.normalBamList = []
        if self.params.tumorBamList is None : self.params.tumorBamList = []

        # bools coming from the ini file need to be cleaned up:
        safeSetBool(self.params,"isWriteRealignedBam")

        if self.params.isWriteCallableRegion :
            self.params.regionsDir=os.path.join(self.params.resultsDir,"regions")
            ensureDir(self.params.regionsDir)

        if self.params.isWriteRealignedBam :
            self.params.realignedDir=os.path.join(self.params.resultsDir,"realigned")
            ensureDir(self.params.realignedDir)

        self.paths = PathInfo(self.params)



    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Strelka workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Strelka workflow version: %s" % (__version__))
        self.setCallMemMb()

        callPreReqs = set()
        callPreReqs |= runCount(self)
        if self.params.isHighDepthFilter :
            callPreReqs |= runDepth(self)

        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)

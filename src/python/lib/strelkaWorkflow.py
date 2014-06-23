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

# add script path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(scriptDir))

# add pyflow path:
# TODO: get a more robust link to the pyflow dir at config time:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from pyflow import WorkflowRunner
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getChromIntervals, getFastaChromOrderSize

from configureUtil import argToBool, getIniSections, dumpIniSections



def getVersion() :
    return "@STARKA_VERSION@"


__version__ = getVersion()



def cleanId(input_id) :
    """
    filter id so that it's safe to use as a pyflow indentifier
    """
    import re
    return re.sub(r'([^a-zA-Z0-9_\-])', "_", input_id)



class GenomeSegment(object) :
    """
    organizes all variables which can change
    with each genomic segment.

    The genomic segment is defined by:

    1. chromosome
    2. begin position (1-indexed closed)
    3. end position (1-indexed closed)
    4. chromosome segment (ie. bin) number (0-indexed)
    """

    def __init__(self,chromIndex,chromLabel,beginPos,endPos,binId,genomeRegion) :
        """
        arguments are the 4 genomic interval descriptors detailed in class documentation
        """
        self.chromLabel = chromLabel
        self.beginPos = beginPos
        self.endPos = endPos
        self.bamRegion = chromLabel + ':' + str(beginPos) + '-' + str(endPos)
        self.binId = binId
        self.binStr = str(binId).zfill(4)
        self.id = chromLabel + "_" + self.binStr

        regionId=cleanId(chromLabel)
        if genomeRegion is not None :
            if genomeRegion['start'] is not None :
                regionId += "-"+str(genomeRegion['start'])
                if genomeRegion['end'] is not None :
                    regionId += "-"+str(genomeRegion['end'])
        self.pyflowId = "chromId_%s_%s_%s" % (str(chromIndex).zfill(3), regionId, self.binStr)



def getNextGenomeSegment(params) :
    """
    generator which iterates through all genomic segments and
    returns a segmentValues object for each one.
    """
    if params.genomeRegionList is None :
        for segval in getChromIntervals(params.chromOrder,params.chromSizes,params.scanSize) :
            yield GenomeSegment(*segval)
    else :
        for genomeRegion in params.genomeRegionList :
            for segval in getChromIntervals(params.chromOrder,params.chromSizes,params.scanSize, genomeRegion) :
                yield GenomeSegment(*segval)



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

    # TMP
    genomeSize=1000000000

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
    segCmd.extend(["-genome-size", str(genomeSize)] )
    segCmd.extend(["-max-indel-size", "50"] )
    segCmd.extend(["-indel-nonsite-match-prob", "0.5"] )
    segCmd.extend(["--somatic-snv-rate", str(self.params.ssnvPrior) ] )
    segCmd.extend(["--shared-site-error-rate", str(self.params.ssnvNoise) ] )
    segCmd.extend(["--shared-site-error-strand-bias-fraction", str(self.params.ssnvNoiseStrandBiasFrac) ] )
    segCmd.extend(["--somatic-indel-rate", str(self.params.sindelPrior) ] )
    segCmd.extend(["--shared-indel-error-rate", str(self.params.sindelNoise) ] )
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

    for bamPath in self.params.normalBamList :
        segCmd.extend(["-bam-file",bamPath])
    for bamPath in self.params.tumorBamList :
        segCmd.extend(["--tumor-bam-file",bamPath])

    tmpSnvPath = self.paths.getTmpSegmentSnvPath(segStr)
    segFiles.snv.append(tmpSnvPath)
    segCmd.extend(["--somatic-snv-file ", tmpSnvPath ] )

    tmpIndelPath = self.paths.getTmpSegmentIndelPath(segStr)
    segFiles.indel.append(tmpIndelPath)
    segCmd.extend(["--somatic-indel-file", tmpIndelPath ] )
    segCmd.extend(["--variant-window-flank-file", "50", self.paths.getTmpSegmentIndelWinPath(segStr) ] )

    if (self.params.maxInputDepth is not None) and (self.params.maxInputDepth > 0) :
        segCmd.extend(["--max-input-depth", str(self.params.maxInputDepth)])

    if self.params.isWriteCallableRegion :
        tmpCallablePath = self.paths.getTmpSegmentRegionPath(segStr)
        segFiles.callable.append(tmpCallablePath)
        segCmd.extend(["--somatic-callable-region-file", tmpCallablePath ])

    if self.params.isWriteRealignedBam :
        segCmd.extend(["-realigned-read-file", self.paths.getTmpUnsortRealignBamPath(segStr, "normal")])
        segCmd.extend(["--tumor-realigned-read-file",self.paths.getTmpUnsortRealignBamPath(segStr, "tumor")])

    if self.params.extraStrelkaArguments is not None :
        for arg in self.params.extraStrelkaArguments.strip().split() :
            segCmd.append(arg)

    segCmd.extend(["--report-file", self.paths.getTmpSegmentReportPath(gseg.pyflowId)])

    if not isFirstSegment :
        segCmd.append("--strelka-skip-header")

    if not self.params.isSkipDepthFilters :
        segCmd.extend(["--strelka-chrom-depth-file", self.paths.getChromDepth()])
        segCmd.extend(["--strelka-max-depth-factor", self.params.depthFilterMultiple])

    nextStepWait = set()

    setTaskLabel=preJoin(taskPrefix,"callGenomeSegment_"+gseg.pyflowId)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)
    nextStepWait.add(setTaskLabel)

    if self.params.isWriteRealignedBam :
        def sortRealignBam(label, sortList) :
            unsorted = self.paths.getTmpUnsortRealignBamPath(segStr, label)
            sorted   = self.paths.getTmpRealignBamPath(segStr, label)
            sortList.append(sorted)
            sortCmd="%s sort %s %s && rm -f %s" % (self.params.samtoolsBin,unsorted,sorted,unsorted)

            sortTaskLabel=preJoin(taskPrefix,"sortRealignedSegment_"+label+"_"+gseg.pyflowId)
            self.addTask(sortTaskLabel,sortCmd,dependencies=setTaskLabel,memMb=self.params.callMemMb)
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

    # create a checkpoint for all segments:
    completeSegmentsTask = self.addTask(preJoin(taskPrefix,"completedAllGenomeSegmetns"),dependencies=graphTasks)

    finishTasks = set()

    def finishVcf(tmpList, output, label) :
        cmd  = "cat " + " ".join(tmpList)
        cmd += " | %s -c >| %s" % (self.params.bgzipBin, output)
        cmd += " && %s -p vcf %s" % (self.params.tabixBin, output)
        finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeVCF"), cmd, dependencies=completeSegmentsTask))

    finishVcf(segFiles.snv, self.paths.getSnvOutputPath(),"SNV")
    finishVcf(segFiles.indel, self.paths.getIndelOutputPath(),"Indel")

    if self.params.isWriteCallableRegion :
        def finishBed(tmpList, output, label):
            cmd  = "cat " + " ".join(tmpList)
            cmd += " | %s -c >| %s" % (self.params.bgzipBin, output)
            cmd += " && %s -p bed %s" % (self.params.tabixBin, output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBED"), cmd, dependencies=completeSegmentsTask))
            
        finishBed(segFiles.callable, self.paths.getRegionOutputPath(), "callableRegions")
    
    if self.params.isWriteRealignedBam :
        def finishBam(tmpList, output, label) :
            assert(len(tmpList) > 0)
            headerTmp = tmpList[0] + "header"
            cmd  = "%s view -H %s >| %s" % (self.params.samtoolsBin, tmpList[0], headerTmp)
            cmd += " && %s merge  -h %s %s " % (self.params.samtoolsBin, headerTmp, output)
            cmd += " ".join(tmpList)
            cmd += " && %s index %s" % (self.params.samtoolsBin, output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))
            
        finishBam(segFiles.normalRealign, self.paths.getRealignedBamPath("normal"), "realignedNormal")
        finishBam(segFiles.tumorRealign, self.paths.getRealignedBamPath("tumor"), "realignedTumor")
    
    # add a tmp folder rm step here....

    nextStepWait = finishTasks

    return nextStepWait



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

    def getTmpSegmentIndelWinPath(self, segStr) :
        return self.getTmpSegmentIndelPath(segStr) + ".window"

    def getTmpSegmentRegionPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "somatic.callable.region.%s.bed" % (segStr))

    def getTmpUnsortRealignBamPath(self, segStr, label) :
        return os.path.join( self.getTmpSegmentDir(), "%s.%s.unsorted.realigned.bam" % (label, segStr))

    def getTmpRealignBamPath(self, segStr, label) :
        return os.path.join( self.getTmpSegmentDir(), "%s.%s.realigned" % (label, segStr))

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


RealignedNormalPath
class StrelkaWorkflow(WorkflowRunner) :
    """
    Strelka somatic small variant calling workflow
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
        if self.params.normalBamList is None : self.params.normalBamList = []
        if self.params.tumorBamList is None : self.params.tumorBamList = []

        # format other:
        self.params.isWriteRealignedBam = argToBool(self.params.isWriteRealignedBam)
        self.params.isWriteCallableRegion = argToBool(self.params.isWriteCallableRegion)
        self.params.isSkipDepthFilters = argToBool(self.params.isSkipDepthFilters)

        # make sure run directory is setup:
        self.params.runDir=os.path.abspath(self.params.runDir)
        ensureDir(self.params.runDir)

        # everything that's not intended to be a final result should dump directories/files in workDir
        self.params.workDir=os.path.join(self.params.runDir,"workspace")
        ensureDir(self.params.workDir)

        # all finalized pretty results get transfered to resultsDir
        self.params.resultsDir=os.path.join(self.params.runDir,"results")
        ensureDir(self.params.resultsDir)
        self.params.statsDir=os.path.join(self.params.resultsDir,"stats")
        ensureDir(self.params.statsDir)
        self.params.variantsDir=os.path.join(self.params.resultsDir,"variants")
        ensureDir(self.params.variantsDir)
        if self.params.isWriteCallableRegion :
            self.params.regionsDir=os.path.join(self.params.resultsDir,"regions")
            ensureDir(self.params.regionsDir)

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

        msg  = "Strelka workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Strelka workflow version: %s" % (__version__))

        callPreReqs = set()
        if not self.params.isSkipDepthFilters :
            depthTasks = runDepth(self)
            callPreReqs |= depthTasks

        callGenome(self, dependencies=callPreReqs)

        """
        graphTaskDependencies = set()

        if not self.params.useExistingAlignStats :
            statsTasks = runStats(self)
            graphTaskDependencies |= statsTasks

        if not ((not self.params.isHighDepthFilter) or self.params.useExistingChromDepths) :
            depthTasks = runDepth(self)
            graphTaskDependencies |= depthTasks

        graphTasks = runLocusGraph(self,dependencies=graphTaskDependencies)

        hygenTasks = runHyGen(self,dependencies=graphTasks)
        """

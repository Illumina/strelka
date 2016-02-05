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
Strelka somatic small variant calling workflow
"""


import os.path
import sys

# add this path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(scriptDir)

# add pyflow path:
sys.path.append(os.path.join(scriptDir,"pyflow"))


from configBuildTimeInfo import workflowVersion
from configureUtil import safeSetBool, getIniSections, dumpIniSections
from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getRmdirCmd, runDepthFromAlignments
from starkaWorkflow import runCount, StarkaCallWorkflow, StarkaWorkflow
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getNextGenomeSegment, bamListCatCmd


__version__ = workflowVersion



def strelkaRunDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    bamList=[]
    if len(self.params.normalBamList) :
        bamList.append(self.params.normalBamList[0])
    elif len(self.params.tumorBamList) :
        bamList.append(self.params.tumorBamList[0])
    else :
        return set()

    outputPath=self.paths.getChromDepth()
    return runDepthFromAlignments(self, bamList, outputPath, taskPrefix, dependencies)



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
    segCmd.extend(["--tier2-mismatch-density-filter-count", "10"] )
    segCmd.append("--tier2-no-filter-unanchored")
    segCmd.extend(["--tier2-indel-nonsite-match-prob", "0.25"] )
    segCmd.append("--tier2-include-singleton")
    segCmd.append("--tier2-include-anomalous")

    segCmd.extend(["--strelka-snv-max-filtered-basecall-frac", str(self.params.snvMaxFilteredBasecallFrac)])
    segCmd.extend(["--strelka-snv-max-spanning-deletion-frac", str(self.params.snvMaxSpanningDeletionFrac)])
    segCmd.extend(["--strelka-snv-min-qss-ref", str(self.params.ssnvQuality_LowerBound)])

    segCmd.extend(["--strelka-indel-max-window-filtered-basecall-frac", str(self.params.indelMaxWindowFilteredBasecallFrac)])
    segCmd.extend(["--strelka-indel-min-qsi-ref", str(self.params.sindelQuality_LowerBound)])

    segCmd.extend(['--indel-error-models-file', self.params.indelErrorModelsFile])
    segCmd.extend(['--indel-error-model-name', self.params.indelErrorModelName])

    if self.params.isEVS :
        segCmd.extend(['--somatic-snv-scoring-model-file', self.params.somaticSnvScoringModelFile])
        if self.params.isSomaticIndelEmpiricalScoring:
            segCmd.extend(['--somatic-indel-scoring-model-file', self.params.somaticIndelScoringModelFile])

    if self.params.isReportEVSFeatures :
        segCmd.append("--report-evs-features")

    for bamPath in self.params.normalBamList :
        segCmd.extend(["-bam-file", bamPath])
    for bamPath in self.params.tumorBamList :
        segCmd.extend(["--tumor-bam-file", bamPath])

    tmpSnvPath = self.paths.getTmpSegmentSnvPath(segStr)
    segFiles.snv.append(tmpSnvPath+".gz")
    segCmd.extend(["--somatic-snv-file ", tmpSnvPath ] )

    tmpIndelPath = self.paths.getTmpSegmentIndelPath(segStr)
    segFiles.indel.append(tmpIndelPath+".gz")
    segCmd.extend(["--somatic-indel-file", tmpIndelPath ] )

    if self.params.isWriteCallableRegion :
        tmpCallablePath = self.paths.getTmpSegmentRegionPath(segStr)
        segFiles.callable.append(tmpCallablePath+".gz")
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

    segCmd.extend(["--report-file", self.paths.getTmpSegmentReportPath(gseg.pyflowId)])

    if not isFirstSegment :
        segCmd.append("--strelka-skip-header")

    if self.params.isHighDepthFilter :
        segCmd.extend(["--strelka-chrom-depth-file", self.paths.getChromDepth()])
        segCmd.extend(["--strelka-max-depth-factor", self.params.depthFilterMultiple])

    if self.params.extraStrelkaArguments is not None :
        for arg in self.params.extraStrelkaArguments.strip().split() :
            segCmd.append(arg)


    nextStepWait = set()

    callTask=preJoin(taskPrefix,"callGenomeSegment_"+gseg.pyflowId)
    self.addTask(callTask,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)

    # fix vcf header to use parent pyflow cmdline instead of random segment command:
    compressWaitFor=callTask
    if isFirstSegment :
        headerFixTask=preJoin(taskPrefix,"fixVcfHeader_"+gseg.pyflowId)
        def getHeaderFixCmd(fileName) :
            tmpName=fileName+".reheader.tmp"
            cmd  = "\"%s\" -E \"%s\"" % (sys.executable, self.params.vcfCmdlineSwapper)
            cmd += ' "' + " ".join(self.params.configCommandLine) + '"'
            cmd += " < \"%s\" > \"%s\" && mv \"%s\" \"%s\"" % (fileName,tmpName,
                                                               tmpName, fileName)
            return cmd

        headerFixCmd  = getHeaderFixCmd(tmpSnvPath)
        headerFixCmd += " && "
        headerFixCmd += getHeaderFixCmd(tmpIndelPath)

        self.addTask(headerFixTask, headerFixCmd, dependencies=callTask, isForceLocal=True)
        compressWaitFor=headerFixTask

    compressTask=preJoin(taskPrefix,"compressSegmentOutput_"+gseg.pyflowId)
    compressCmd="\"%s\" \"%s\" && \"%s\" \"%s\"" % (self.params.bgzipBin, tmpSnvPath, self.params.bgzipBin, tmpIndelPath)
    if self.params.isWriteCallableRegion :
        compressCmd += " && \"%s\" \"%s\"" % (self.params.bgzipBin, self.paths.getTmpSegmentRegionPath(segStr))

    self.addTask(compressTask, compressCmd, dependencies=compressWaitFor, isForceLocal=True)
    nextStepWait.add(compressTask)

    if self.params.isWriteRealignedBam :
        def sortRealignBam(label, sortList) :
            unsorted = self.paths.getTmpUnsortRealignBamPath(segStr, label)
            sorted   = self.paths.getTmpRealignBamPath(segStr, label)
            sortList.append(sorted)

            # adjust sorted to remove the ".bam" suffix
            sorted = sorted[:-4]
            sortCmd="\"%s\" sort \"%s\" \"%s\" && rm -f \"%s\"" % (self.params.samtoolsBin,unsorted,sorted,unsorted)

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

    finishTasks.add(self.concatIndexVcf(taskPrefix, completeSegmentsTask, segFiles.snv,
                                        self.paths.getSnvOutputPath(),"SNV"))
    finishTasks.add(self.concatIndexVcf(taskPrefix, completeSegmentsTask, segFiles.indel,
                                        self.paths.getIndelOutputPath(),"Indel"))

    if self.params.isWriteCallableRegion :
        finishTasks.add(self.concatIndexBed(taskPrefix, completeSegmentsTask, segFiles.callable,
                                            self.paths.getRegionOutputPath(), "callableRegions"))

    if self.params.isWriteRealignedBam :
        def finishBam(tmpList, output, label) :
            cmd = bamListCatCmd(self.params.samtoolsBin,tmpList,output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))

        finishBam(segFiles.normalRealign, self.paths.getRealignedBamPath("normal"), "realignedNormal")
        finishBam(segFiles.tumorRealign, self.paths.getRealignedBamPath("tumor"), "realignedTumor")

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
            callPreReqs |= strelkaRunDepthFromAlignments(self)

        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)

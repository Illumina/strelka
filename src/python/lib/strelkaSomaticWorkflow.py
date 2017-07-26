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
from configureUtil import safeSetBool
from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getRmdirCmd, getDepthFromAlignments
from strelkaSharedWorkflow import getTotalKnownReferenceSize, runCount, SharedPathInfo, \
                           StrelkaSharedCallWorkflow, StrelkaSharedWorkflow
from workflowUtil import ensureDir, preJoin, bamListCatCmd


__version__ = workflowVersion



def strelkaSomaticGetDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    bamList=[]
    if len(self.params.normalBamList) :
        bamList.append(self.params.normalBamList[0])
    elif len(self.params.tumorBamList) :
        bamList.append(self.params.tumorBamList[0])
    else :
        return set()

    outputPath=self.paths.getChromDepth()
    return getDepthFromAlignments(self, bamList, outputPath, taskPrefix, dependencies)



class TempVariantCallingSegmentFiles :
    def __init__(self) :
        self.snv = []
        self.indel = []
        self.callable = []
        self.normalRealign = []
        self.tumorRealign = []
        self.stats = []



def callGenomeSegment(self, gsegGroup, segFiles, taskPrefix="", dependencies=None) :

    assert(len(gsegGroup) != 0)
    gid=gsegGroup[0].id
    if len(gsegGroup) > 1 :
        gid += "_to_"+gsegGroup[-1].id

    isFirstSegment = (len(segFiles.snv) == 0)

    segCmd = [ self.params.strelkaSomaticBin ]

    self.appendCommonGenomeSegmentCommandOptions(gsegGroup, segCmd)

    segCmd.extend(["-min-mapping-quality",str(self.params.minTier1Mapq)])
    segCmd.extend(["-min-qscore","0"])
    segCmd.extend(["-max-window-mismatch", "3", "20" ])
    segCmd.extend(["-indel-nonsite-match-prob", "0.5"] )
    segCmd.extend(["--somatic-snv-rate", str(self.params.ssnvPrior) ] )
    segCmd.extend(["--shared-site-error-rate", str(self.params.ssnvNoise) ] )
    segCmd.extend(["--shared-site-error-strand-bias-fraction", str(self.params.ssnvNoiseStrandBiasFrac) ] )
    segCmd.extend(["--somatic-indel-rate", str(self.params.sindelPrior) ] )
    segCmd.extend(["--shared-indel-error-factor", str(self.params.sindelNoiseFactor)])
    segCmd.extend(["--tier2-min-mapping-quality", str(self.params.minTier2Mapq) ] )
    segCmd.extend(["--tier2-mismatch-density-filter-count", "10"] )
    segCmd.extend(["--tier2-indel-nonsite-match-prob", "0.25"] )
    segCmd.append("--tier2-include-singleton")
    segCmd.append("--tier2-include-anomalous")

    segCmd.extend(["--strelka-snv-max-filtered-basecall-frac", str(self.params.snvMaxFilteredBasecallFrac)])
    segCmd.extend(["--strelka-snv-max-spanning-deletion-frac", str(self.params.snvMaxSpanningDeletionFrac)])
    segCmd.extend(["--strelka-snv-min-qss-ref", str(self.params.ssnvQuality_LowerBound)])

    segCmd.extend(["--strelka-indel-max-window-filtered-basecall-frac", str(self.params.indelMaxWindowFilteredBasecallFrac)])
    segCmd.extend(["--strelka-indel-min-qsi-ref", str(self.params.sindelQuality_LowerBound)])

    segCmd.extend(["--ssnv-contam-tolerance", str(self.params.ssnvContamTolerance) ] )
    segCmd.extend(["--indel-contam-tolerance", str(self.params.indelContamTolerance) ] )

    if self.params.isEVS :
        if self.params.snvScoringModelFile is not None :
            segCmd.extend(['--somatic-snv-scoring-model-file', self.params.snvScoringModelFile])
        if self.params.indelScoringModelFile is not None :
            segCmd.extend(['--somatic-indel-scoring-model-file', self.params.indelScoringModelFile])

    for bamPath in self.params.normalBamList :
        segCmd.extend(["--normal-align-file", bamPath])
    for bamPath in self.params.tumorBamList :
        segCmd.extend(["--tumor-align-file", bamPath])

    tmpSnvPath = self.paths.getTmpSegmentSnvPath(gid)
    segFiles.snv.append(tmpSnvPath+".gz")
    segCmd.extend(["--somatic-snv-file ", tmpSnvPath ] )

    tmpIndelPath = self.paths.getTmpSegmentIndelPath(gid)
    segFiles.indel.append(tmpIndelPath+".gz")
    segCmd.extend(["--somatic-indel-file", tmpIndelPath ] )

    if self.params.isOutputCallableRegions :
        tmpCallablePath = self.paths.getTmpSegmentRegionPath(gid)
        segFiles.callable.append(tmpCallablePath+".gz")
        segCmd.extend(["--somatic-callable-regions-file", tmpCallablePath ])

    if self.params.isWriteRealignedBam :
        segCmd.extend(["-realigned-read-file", self.paths.getTmpUnsortRealignBamPath(gid, "normal")])
        segCmd.extend(["--tumor-realigned-read-file",self.paths.getTmpUnsortRealignBamPath(gid, "tumor")])

    def addListCmdOption(optList,arg) :
        if optList is None : return
        for val in optList :
            segCmd.extend([arg, val])

    addListCmdOption(self.params.noiseVcfList, '--noise-vcf')

    segFiles.stats.append(self.paths.getTmpRunStatsPath(gid))
    segCmd.extend(["--stats-file", segFiles.stats[-1]])

    if not isFirstSegment :
        segCmd.append("--strelka-skip-header")

    if self.params.isHighDepthFilter :
        segCmd.extend(["--strelka-chrom-depth-file", self.paths.getChromDepth()])
        segCmd.extend(["--strelka-max-depth-factor", self.params.depthFilterMultiple])


    nextStepWait = set()

    callTask=preJoin(taskPrefix,"callGenomeSegment_"+gid)
    self.addTask(callTask,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)

    # fix vcf header to use parent pyflow cmdline instead of random segment command:
    compressWaitFor=callTask
    if isFirstSegment :
        headerFixTask=preJoin(taskPrefix,"fixVcfHeader_"+gid)
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

    compressTask=preJoin(taskPrefix,"compressSegmentOutput_"+gid)
    compressCmd="\"%s\" \"%s\" && \"%s\" \"%s\"" % (self.params.bgzipBin, tmpSnvPath, self.params.bgzipBin, tmpIndelPath)
    if self.params.isOutputCallableRegions :
        compressCmd += " && \"%s\" \"%s\"" % (self.params.bgzipBin, self.paths.getTmpSegmentRegionPath(gid))

    self.addTask(compressTask, compressCmd, dependencies=compressWaitFor, isForceLocal=True)
    nextStepWait.add(compressTask)

    if self.params.isWriteRealignedBam :
        def sortRealignBam(label, sortList) :
            unsorted = self.paths.getTmpUnsortRealignBamPath(gid, label)
            sorted   = self.paths.getTmpRealignBamPath(gid, label)
            sortList.append(sorted)

            # adjust sorted to remove the ".bam" suffix
            sorted = sorted[:-4]
            sortCmd="\"%s\" sort \"%s\" \"%s\" && rm -f \"%s\"" % (self.params.samtoolsBin,unsorted,sorted,unsorted)

            sortTaskLabel=preJoin(taskPrefix,"sortRealignedSegment_"+label+"_"+gid)
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
    dirTask=self.addTask(preJoin(taskPrefix,"makeTmpDir"), getMkdirCmd() + [tmpSegmentDir],
                         dependencies=dependencies, isForceLocal=True)

    segmentTasks = set()

    segFiles = TempVariantCallingSegmentFiles()

    for gsegGroup in self.getStrelkaGenomeSegmentGroupIterator() :
        segmentTasks |= callGenomeSegment(self, gsegGroup, segFiles, dependencies=dirTask)

    if len(segmentTasks) == 0 :
        raise Exception("No genome regions to analyze. Possible target region parse error.")

    # create a checkpoint for all segments:
    completeSegmentsTask = self.addTask(preJoin(taskPrefix,"completedAllGenomeSegments"),dependencies=segmentTasks)

    finishTasks = set()

    finishTasks.add(self.concatIndexVcf(taskPrefix, completeSegmentsTask, segFiles.snv,
                                        self.paths.getSnvOutputPath(),"SNV"))
    finishTasks.add(self.concatIndexVcf(taskPrefix, completeSegmentsTask, segFiles.indel,
                                        self.paths.getIndelOutputPath(),"Indel"))

    # merge segment stats:
    finishTasks.add(self.mergeRunStats(taskPrefix,completeSegmentsTask, segFiles.stats))

    if self.params.isOutputCallableRegions :
        finishTasks.add(self.concatIndexBed(taskPrefix, completeSegmentsTask, segFiles.callable,
                                            self.paths.getRegionOutputPath(), "callableRegions"))

    if self.params.isWriteRealignedBam :
        def finishBam(tmpList, output, label) :
            cmd = bamListCatCmd(self.params.samtoolsBin,tmpList,output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))

        finishBam(segFiles.normalRealign, self.paths.getRealignedBamPath("normal"), "realignedNormal")
        finishBam(segFiles.tumorRealign, self.paths.getRealignedBamPath("tumor"), "realignedTumor")

    if not self.params.isRetainTempFiles :
        rmTmpCmd = getRmdirCmd() + [tmpSegmentDir]
        self.addTask(preJoin(taskPrefix,"removeTmpDir"), rmTmpCmd, dependencies=finishTasks, isForceLocal=True)

    nextStepWait = finishTasks

    return nextStepWait



class CallWorkflow(StrelkaSharedCallWorkflow) :
    """
    A separate workflow is setup around callGenome() so that the workflow execution can be
    delayed until the ref count exists
    """

    def __init__(self, params, paths) :
        super(CallWorkflow,self).__init__(params)
        self.paths = paths

    def workflow(self) :
        self.params.totalKnownReferenceSize = getTotalKnownReferenceSize(self.paths.getReferenceSizePath())
        callGenome(self)



class PathInfo(SharedPathInfo):
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        super(PathInfo,self).__init__(params)

    def getTmpSegmentSnvPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "somatic.snvs.unfiltered.%s.vcf" % (segStr))

    def getTmpSegmentIndelPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "somatic.indels.unfiltered.%s.vcf" % (segStr))

    def getTmpSegmentRegionPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "somatic.callable.regions.%s.bed" % (segStr))

    def getTmpUnsortRealignBamPath(self, segStr, label) :
        return os.path.join( self.getTmpSegmentDir(), "%s.%s.unsorted.realigned.bam" % (label, segStr))

    def getTmpRealignBamPath(self, segStr, label) :
        return os.path.join( self.getTmpSegmentDir(), "%s.%s.realigned.bam" % (label, segStr))

    def getSnvOutputPath(self) :
        return os.path.join( self.params.variantsDir, "somatic.snvs.vcf.gz")

    def getIndelOutputPath(self) :
        return os.path.join( self.params.variantsDir, "somatic.indels.vcf.gz")

    def getRegionOutputPath(self) :
        return os.path.join( self.params.regionsDir, 'somatic.callable.regions.bed.gz')

    def getRealignedBamPath(self, label) :
        return os.path.join( self.params.realignedDir, '%s.realigned.bam' % (label))



class StrelkaSomaticWorkflow(StrelkaSharedWorkflow) :
    """
    Strelka somatic small variant calling workflow
    """

    def __init__(self,params) :
        global PathInfo
        super(StrelkaSomaticWorkflow,self).__init__(params, PathInfo)

        # format bam lists:
        if self.params.normalBamList is None : self.params.normalBamList = []
        if self.params.tumorBamList is None : self.params.tumorBamList = []

        # bools coming from the ini file need to be cleaned up:
        safeSetBool(self.params,"isWriteRealignedBam")

        if self.params.isOutputCallableRegions :
            self.params.regionsDir=os.path.join(self.params.resultsDir,"regions")
            ensureDir(self.params.regionsDir)

        if self.params.isWriteRealignedBam :
            self.params.realignedDir=os.path.join(self.params.resultsDir,"realigned")
            ensureDir(self.params.realignedDir)


    def getSuccessMessage(self) :
        "Message to be included in email for successful runs"

        msg  = "Strelka somatic workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Strelka somatic workflow version: %s" % (__version__))
        self.setCallMemMb()

        callPreReqs = set()
        callPreReqs.add(runCount(self))
        if self.params.isHighDepthFilter :
            callPreReqs |= strelkaSomaticGetDepthFromAlignments(self)

        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)

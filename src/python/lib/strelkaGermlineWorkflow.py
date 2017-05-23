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
from configureUtil import safeSetBool
from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getRmdirCmd, getDepthFromAlignments
from strelkaSharedWorkflow import DeepCopyProtector, runCount, SharedPathInfo, \
                           StrelkaSharedCallWorkflow, StrelkaSharedWorkflow
from workflowUtil import ensureDir, preJoin, bamListCatCmd, GenomeSegment, getChromIntervals

__version__ = workflowVersion



def strelkaGermlineGetDepthFromAlignments(self,taskPrefix="getChromDepth",dependencies=None):
    bamList=[]
    if len(self.params.bamList) :
        bamList = self.params.bamList
    else :
        return set()

    outputPath=self.paths.getChromDepth()
    return getDepthFromAlignments(self, bamList, outputPath, taskPrefix, dependencies)



class TempVariantCallingSegmentFilesPerSample :
    def __init__(self) :
        self.gvcf = []


class TempVariantCallingSegmentFiles :
    def __init__(self, sampleCount) :
        self.variants = []
        self.bamRealign = []
        self.stats = []
        self.sample = [TempVariantCallingSegmentFilesPerSample() for _ in range(sampleCount)]


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
        segCmd.extend(['-bsnp-diploid-het-bias', '0.45'])
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

    if self.params.isIndelErrorRateEstimated :
        for bamIndex in range(len(self.params.bamList)) :
            segCmd.extend(['--indel-error-models-file', self.paths.getIndelErrorModelPath(bamIndex)])
    else :
        segCmd.extend(['--indel-error-models-file', self.params.indelErrorRateDefault])

    segCmd.extend(['--theta-file', self.params.thetaParamFile])

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
        def finishBam(tmpList, output, label) :
            cmd = bamListCatCmd(self.params.samtoolsBin, tmpList, output)
            finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeBAM"), cmd, dependencies=completeSegmentsTask))

        finishBam(segFiles.bamRealign, self.paths.getRealignedBamPath(), "realigned")

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

    def __init__(self, params, paths, dynamicParams) :
        super(CallWorkflow,self).__init__(params, dynamicParams)
        self.paths = paths

    def workflow(self) :
        callGenome(self)



def countGenomeSegment(self, sampleIndex, gseg, segFiles, taskPrefix="", dependencies=None) :
    """
    Extract sequencing error count data from the genome segment specified by gseg.bamRegion
    """

    segStr = str(gseg.id)

    segCmd = [ self.params.getCountsBin ]

    segCmd.extend(["--region", gseg.bamRegion])
    segCmd.extend(["--ref", self.params.referenceFasta ])
    segCmd.extend(["-genome-size", str(self.dynamicParams.totalKnownSize)] )
    segCmd.extend(["-max-indel-size", "50"] )

    segFiles.counts.append(self.paths.getTmpSegmentErrorCountsPath(sampleIndex, segStr))
    segCmd.extend(["--counts-file", segFiles.counts[-1]])

    segFiles.nonEmptySiteCounts.append(self.paths.getTmpSegmentNonemptySiteCountsPath(sampleIndex, segStr))
    segCmd.extend(["--nonempty-site-count-file", segFiles.nonEmptySiteCounts[-1]])

    bamPath = self.params.bamList[sampleIndex]
    segCmd.extend(["--align-file", bamPath])

    if self.params.isHighDepthFilter :
        segCmd.extend(["--chrom-depth-file", self.paths.getChromDepth()])

    def addListCmdOption(optList,arg) :
        if optList is None : return
        for val in optList :
            segCmd.extend([arg, val])

    addListCmdOption(self.params.indelCandidatesList, '--candidate-indel-input-vcf')
    addListCmdOption(self.params.forcedGTList, '--force-output-vcf')

    setTaskLabel=preJoin(taskPrefix,"countErrors_"+gseg.id)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)

    return setTaskLabel



def lockMethod(f):
    """
    Method decorator acquires/releases object's lock
    """

    def wrapped(self, *args, **kw):
        import threading

        if not hasattr(self,"lock") :
            self.lock = threading.RLock()

        self.lock.acquire()
        try:
            return f(self, *args, **kw)
        finally:
            self.lock.release()
    return wrapped



class SyncronizedAccumulator(DeepCopyProtector) :

    def __init__(self) :
        self._values = []

    @lockMethod
    def addOrderedValue(self, index, value):
        while index+1 > len(self._values) :
            self._values.append(None)
        self._values[index] = value

    @lockMethod
    def totalValue(self):
        count = 0
        sum = 0
        for v in self._values :
            count += 1
            if v is None : continue
            sum += v
        return (count, sum)

    @lockMethod
    def totalContinuousValue(self):
        count = 0
        sum = 0
        for v in self._values :
            if v is None : break
            count += 1
            sum += v
        return (count, sum)

    @lockMethod
    def countTasksRequiredToReachTarget(self, targetVal):
        """
        Return the tuple (taskCount, isContinuous), where:

        taskCount is the smallest index N such that sum of values [1:N] is >= targetVal, or None if no such value exists
        isContinuous is True if all tasks in range [1:N] are present in the value
        """
        taskCount = 0
        sum = 0
        isContinuous = True

        # Handle the edge case, targetVal <= 0
        if sum >= targetVal :
            return (taskCount, isContinuous)

        for v in self._values :
            taskCount += 1
            if v is None :
                isContinuous = False
                continue
            sum += v
            if sum >= targetVal :
                return (taskCount, isContinuous)

        return (None, isContinuous)



class UpdateCompletedTaskTrackerWorkflow(WorkflowRunner) :
    """
    This workflow reads the number of non-empty sites returned by each task and
    updates pyflow nonempty site tracking data-structures.
    """

    def __init__(self, taskIndex, inputFile, completedTaskTracker) :
        self.taskIndex=taskIndex
        self.inputFile=inputFile
        self.completedTaskTracker = completedTaskTracker

    def workflow(self) :
        infp = open(self.inputFile, "rb")
        nonEmptySiteCount = int(infp.read().strip().split()[1])
        infp.close()
        self.completedTaskTracker.addOrderedValue(self.taskIndex, nonEmptySiteCount)



def countEmTillTheCountinsDone(self, estimationIntervals, sampleIndex, segFiles, taskPrefix="", dependencies=None) :
    """
    This routine organizes the process of launching sequence error count jobs until a specific total
    evidence count has been gathered from the genome (or not genome segments are left)

    Note that this function will not return until its tasks are completed, so it will block conventional task parallelization
    """

    class Constants :
        Megabase = 1000000
        totalContinuousNonEmptySiteTarget = 50 * Megabase

    maxTaskCount = self.getNCores()
    assert(maxTaskCount > 0)

    taskByIndex = []
    allTasks = set()
    completedTasks = set()
    completedTaskTracker = SyncronizedAccumulator()

    class Shared :
        def __init__(self):
            self.lowestCanceledTaskIndex = None

    shared = Shared()

    def launchNextTask() :
        """
        Launch the next task in queue for this sample

        Return false if there are no more jobs to launch
        """
        taskIndex = len(allTasks)
        if taskIndex >= len(estimationIntervals) : return False

        gseg = estimationIntervals[taskIndex]
        countTask = countGenomeSegment(self, sampleIndex, gseg, segFiles,
                                       taskPrefix=taskPrefix, dependencies=dependencies)

        allTasks.add(countTask)
        taskByIndex.append(countTask)

        updateTaskLabel=preJoin(taskPrefix,"trackCounts_"+gseg.id)
        updateWorkflow =  UpdateCompletedTaskTrackerWorkflow(taskIndex, segFiles.nonEmptySiteCounts[-1], completedTaskTracker)
        self.addWorkflowTask(updateTaskLabel, updateWorkflow, dependencies=countTask)

        return True


    def updateCompletedTasks() :
        for task in allTasks :
            if task in completedTasks : continue
            (isDone, isError) = self.isTaskDone(task)
            if not isDone : continue
            if isError :
                raise Exception("Task %s failed." % (task))
            completedTasks.add(task)


    def stopRunningExtraTasks(nTargetTasks) :
        """
        Cancel tasks we don't need anymore
        """

        if shared.lowestCanceledTaskIndex is None :
            shared.lowestCanceledTaskIndex = len(taskByIndex)

        assert(nTargetTasks <= shared.lowestCanceledTaskIndex)

        for task in taskByIndex[nTargetTasks:shared.lowestCanceledTaskIndex] :
            self.cancelTaskTree(task)

        shared.lowestCanceledTaskIndex = nTargetTasks


    import time

    while True :
        # Loop until the total required counts have been found in this sample, or there are no more counting
        # tasks to launch.

        (nTargetTasks,isContinuous) = completedTaskTracker.countTasksRequiredToReachTarget(Constants.totalContinuousNonEmptySiteTarget)
        if nTargetTasks is not None :
            # nTargetTasks always provides an upper bound on the total number of tasks required to reach the goal,
            # so we can start cancelling tasks we know won't be needed. This upper bound becomes exact when
            # isContinuous is True
            stopRunningExtraTasks(nTargetTasks)

            # If the target tasks were found in one continuous block, then we're done!:
            if isContinuous :
                # Reduce completed file lists to just the content we're going to keep:
                del segFiles.counts[nTargetTasks:]
                del segFiles.nonEmptySiteCounts[nTargetTasks:]
                break

        updateCompletedTasks()
        runningTaskCount = len(allTasks)-len(completedTasks)
        assert(runningTaskCount >= 0)

        #self.flowLog("Sample%i Completed/Running/All/Input tasks: %i %i %i %i" %
        #             (sampleIndex, len(completedTasks), runningTaskCount, len(allTasks), len(estimationIntervals)))

        if len(allTasks) < len(estimationIntervals) :
            # launch new tasks unless it is already clear the total nonempty site threshold will be met
            if nTargetTasks is None :
                numberOfTasksToLaunch = max(maxTaskCount-runningTaskCount,0)
                for _ in range(numberOfTasksToLaunch) :
                    launchStatus = launchNextTask()
                    if not launchStatus : break
        elif runningTaskCount == 0 :
            # All task intervals have been run but total count has not been reached
            break

        time.sleep(1)

    # This function doesn't return until its tasks are complete or cancelled, so nothing to wait for:
    waitForTasks = set()
    return waitForTasks



def mergeSequenceErrorCounts(self, sampleIndex, segmentErrorCountFiles, taskPrefix="", dependencies=None) :
    """
    Given sequencing error counts generated from multiple genome regions, merge these into a single error count set
    """

    runMergeLabel=preJoin(taskPrefix,"mergeCounts")
    runMergeCmd=[self.params.mergeCountsBin]
    runMergeCmd.extend(["--output-file", self.paths.getErrorCountsOutputPath(sampleIndex)])
    for segmentErrorCountFile in segmentErrorCountFiles :
        runMergeCmd.extend(["--counts-file",segmentErrorCountFile])
    return self.addTask(runMergeLabel, runMergeCmd, dependencies=dependencies, isForceLocal=True)



def estimateParametersFromErrorCounts(self, sampleIndex, taskPrefix="", dependencies=None) :
    """
    Estimate variant error parameters from sequencing error count data
    """

    runEstimateLabel=preJoin(taskPrefix,"estimateVariantErrorRatesBin")
    runEstimateCmd=[self.params.estimateVariantErrorRatesBin]
    runEstimateCmd.extend(["--counts-file", self.paths.getErrorCountsOutputPath(sampleIndex)])
    runEstimateCmd.extend(["--theta-file",self.params.thetaParamFile])
    runEstimateCmd.extend(["--output-file", self.paths.getIndelErrorModelPath(sampleIndex)])
    runEstimateCmd.extend(["--fallback-file",self.params.indelErrorRateDefault])
    return self.addTask(runEstimateLabel, runEstimateCmd, dependencies=dependencies, isForceLocal=True)


class TempSequenceErrorCountSegmentFiles :
    def __init__(self) :
        self.counts = []
        self.nonEmptySiteCounts = []


def getSequenceErrorEstimatesForSample(self, estimationIntervals, sampleIndex, taskPrefix="", dependencies=None):
    """
    Count sequencing errors in one sample and use these to estimate sample error parameters
    """

    segmentTasks = set()

    segFiles = TempSequenceErrorCountSegmentFiles()

    # Launch tasks until the required counts are found
    segmentTasks |= countEmTillTheCountinsDone(self, estimationIntervals, sampleIndex, segFiles, taskPrefix, dependencies)

    # create a checkpoint for all segments:
    completeSegmentsTask = self.addTask(preJoin(taskPrefix,"completedAllGenomeSegments"),dependencies=segmentTasks)

    # merge segment stats:
    mergeCountsTask = mergeSequenceErrorCounts(self, sampleIndex, segFiles.counts,
                                               taskPrefix=taskPrefix, dependencies=completeSegmentsTask)

    # get error parameters:
    estimateTask = estimateParametersFromErrorCounts(self, sampleIndex,
                                                     taskPrefix=taskPrefix, dependencies=mergeCountsTask)

    nextStepWait = set()
    nextStepWait.add(estimateTask)
    return nextStepWait


class EstimateSequenceErrorWorkflowForSample(WorkflowRunner) :
    """
    A separate workflow is setup around the error estimation process for each sample so that:
    (1) sequence error estimation can run in parallel for all samples
    (2) workflow execution can be delayed until the ref count exists
    """

    def __init__(self, params, paths, dynamicParams, estimationIntervals, sampleIndex) :
        self.params = params
        self.paths = paths
        self.dynamicParams = dynamicParams
        self.estimationIntervals = estimationIntervals
        self.sampleIndex = sampleIndex

    def workflow(self) :
        getSequenceErrorEstimatesForSample(self, self.estimationIntervals, self.sampleIndex)



def getErrorEstimationIntervals(params) :
    """
    Return the list of genome intervals to use for error estimation
    """

    # Setup the estimation region sampling list for all samples, estimation uses smaller genome
    # segments than variant calling, and shuffles these with a constant seed to make results deterministic
    class Constants :
        Megabase = 1000000
        errorEstimationMinChromSize = params.errorEstimationMinChromMb * Megabase
        errorEstimationChunkSize = 2 * Megabase

    largeChroms = [chrom for chrom in params.chromOrder if params.chromSizes[chrom] >= Constants.errorEstimationMinChromSize]
    estimationIntervals = [GenomeSegment(*interval) for interval in getChromIntervals(largeChroms, params.chromSizes, Constants.errorEstimationChunkSize)]

    # Create a local random instance such that the fixed-seed shuffle procedure works as expected
    # in a multi-threaded context.
    import random
    localRandom = random.Random()
    localRandom.seed(4000)

    localRandom.shuffle(estimationIntervals)

    return estimationIntervals



def getSequenceErrorEstimates(self, taskPrefix="", dependencies=None):
    """
    Count sequence errors and use these to estimate error parameters
    """

    mkDirTask = preJoin(taskPrefix,"makeTmpDir")
    tmpErrorEstimationDir = self.paths.getTmpErrorEstimationDir()
    mkDirCmd = getMkdirCmd() + [tmpErrorEstimationDir]
    self.addTask(mkDirTask, mkDirCmd, dependencies=dependencies, isForceLocal=True)

    estimationIntervals = getErrorEstimationIntervals(self.params)
    assert(len(estimationIntervals) != 0)

    # The count and estimation processes are currently independent for each sample
    sampleTasks = set()
    for sampleIndex in range(len(self.params.bamList)) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        sampleTask=preJoin(taskPrefix,"Sample"+sampleIndexStr)
        workflow=EstimateSequenceErrorWorkflowForSample(self.params, self.paths, self.dynamicParams, estimationIntervals, sampleIndex)
        sampleTasks.add(self.addWorkflowTask(sampleTask, workflow, dependencies=mkDirTask))

    #if not self.params.isRetainTempFiles :
    #    rmTmpCmd = getRmdirCmd() + [tmpErrorEstimationDir]
    #    self.addTask(preJoin(taskPrefix,"removeTmpDir"), rmTmpCmd, dependencies=sampleTasks, isForceLocal=True)

    nextStepWait = sampleTasks
    return nextStepWait



class PathInfo(SharedPathInfo):
    """
    Object to centralize shared workflow path names
    """

    def __init__(self, params) :
        super(PathInfo,self).__init__(params)

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

    def getTmpSegmentNonemptySiteCountsPath(self, sampleIndex, segStr) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join( self.getTmpErrorEstimationDir(), "nonEmptySiteCounts.Sample%s.%s.tsv" % (sampleIndexStr,segStr))

    def getTmpSegmentErrorCountsPath(self, sampleIndex, segStr) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join( self.getTmpErrorEstimationDir(), "sequenceErrorCounts.Sample%s.%s.bin" % (sampleIndexStr,segStr))

    def getErrorCountsOutputPath(self, sampleIndex) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join( self.params.workDir, "sequenceErrorCounts.Sample%s.bin" % (sampleIndexStr))

    def getIndelErrorModelPath(self, sampleIndex) :
        sampleIndexStr = str(sampleIndex).zfill(3)
        return os.path.join( self.params.statsDir, "sequenceErrorModel.Sample%s.json" % (sampleIndexStr))



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
        estimatePreReqs |= runCount(self)
        if self.params.isHighDepthFilter :
            estimatePreReqs |= strelkaGermlineGetDepthFromAlignments(self)

        if self.params.isIndelErrorRateEstimated :
            callPreReqs |= getSequenceErrorEstimates(self, taskPrefix="EstimateSequenceError", dependencies=estimatePreReqs)
        else :
            callPreReqs = estimatePreReqs

        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths, self.dynamicParams), dependencies=callPreReqs)

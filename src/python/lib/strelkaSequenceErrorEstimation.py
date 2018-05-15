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
Strelka workflow components providing dynamic sequence error estimation capabilities
"""


import os.path
import sys

# add this path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(scriptDir)

# add pyflow path:
pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from pyflow import WorkflowRunner
from sharedWorkflow import getMkdirCmd, getRmdirCmd
from workflowUtil import preJoin, GenomeSegment, getChromIntervals



def countGenomeSegment(self, sampleIndex, gseg, segFiles, taskPrefix="", dependencies=None) :
    """
    Extract sequencing error count data from the genome segment specified by gseg.bamRegion
    """

    genomeSegmentLabel = gseg.id

    segCmd = [ self.params.getCountsBin ]

    segCmd.extend(["--region", gseg.bamRegion])
    segCmd.extend(["--ref", self.params.referenceFasta ])
    segCmd.extend(["--max-indel-size", self.params.maxIndelSize])

    segFiles.counts.append(self.paths.getTmpSegmentAlleleCountsPath(sampleIndex, genomeSegmentLabel))
    segCmd.extend(["--counts-file", segFiles.counts[-1]])

    segFiles.nonEmptySiteCounts.append(self.paths.getTmpSegmentNonemptySiteCountsPath(sampleIndex, genomeSegmentLabel))
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



def countAllEligibleSequenceEvidence(self, estimationIntervals, sampleIndex, segFiles, taskPrefix="", dependencies=None) :
    """
    This routine launches error counting jobs for a sample to cover all available input data.
    """

    segmentTasks = set()

    for gseg in estimationIntervals :
        segmentTasks.add(countGenomeSegment(self, sampleIndex, gseg, segFiles,
                                            taskPrefix=taskPrefix, dependencies=dependencies))

    return segmentTasks



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



class DeepCopyProtector(object) :
    """
    Any data attached to this object will remain aliased through a deepcopy operation

    Overloading __copy__ is provided here just to ensure that deep/shallow copy semantics are identical.
    """
    def __copy__(self) :
        return self

    def __deepcopy__(self, dict) :
        return self



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



def countSequenceEvidenceUntilTargetIsReached(self, estimationIntervals, sampleIndex, segFiles,
                                              taskPrefix="", dependencies=None) :
    """
    This routine organizes the process of launching sequence error count jobs until a specific total
    evidence count has been gathered from the genome (or no genome segments are left)

    Note that this function will not return until its tasks are completed, so it will block conventional task
    parallelization
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

        #self.flowLog("ZZZ Sample%i launching taskIndex/task %i %s" % (sampleIndex, taskIndex, countTask))

        allTasks.add(countTask)
        taskByIndex.append(countTask)

        updateTaskLabel=preJoin(taskPrefix,"trackCounts_"+gseg.id)
        updateWorkflow = UpdateCompletedTaskTrackerWorkflow(taskIndex, segFiles.nonEmptySiteCounts[-1],
                                                            completedTaskTracker)
        self.addWorkflowTask(updateTaskLabel, updateWorkflow, dependencies=countTask, isEphemeral=True)

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

        #self.flowLog("ZZZ Sample%i cancelQuery allTasks nTargetTasks lowestCanceledTaskIndex %i %i %i" %
        #             (sampleIndex, len(taskByIndex), nTargetTasks, shared.lowestCanceledTaskIndex))
        for task in taskByIndex[nTargetTasks:shared.lowestCanceledTaskIndex] :
            #self.flowLog("ZZZ Sample%i canceling task %s" % (sampleIndex, task))
            self.cancelTaskTree(task)

        shared.lowestCanceledTaskIndex = nTargetTasks


    import time

    while (not self.isWorkflowStopping()) :
        # Loop until the total required counts have been found in this sample, or there are no more counting
        # tasks to launch.

        #(count,sum) = completedTaskTracker.totalValue()
        #self.flowLog("ZZZ Sample%i total sum (count) %i (%i)" % (sampleIndex, count, sum))

        #(count,sum) = completedTaskTracker.totalContinuousValue()
        #self.flowLog("ZZZ Sample%i continuous sum (count) %i (%i)" % (sampleIndex, count, sum))

        (nTargetTasks,isContinuous) = completedTaskTracker.countTasksRequiredToReachTarget(Constants.totalContinuousNonEmptySiteTarget)
        #self.flowLog("ZZZ Sample%i nTargetTasks isContinuous %s %s" % (sampleIndex, str(nTargetTasks), str(isContinuous)))

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

        #self.flowLog("ZZZ Sample%i Completed/Running/All/Input tasks: %i %i %i %i" %
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



def mergeSequenceAlleleCounts(self, sampleIndex, segmentAlleleCountsFiles, taskPrefix="", dependencies=None) :
    """
    Given sequencing error counts generated from multiple genome regions, merge these into a single error count set
    """

    runMergeLabel=preJoin(taskPrefix,"mergeCounts")
    runMergeCmd=[self.params.mergeCountsBin]
    runMergeCmd.extend(["--output-file", self.paths.getAlleleCountsOutputPath(sampleIndex)])
    for segmentAlleleCountsFile in segmentAlleleCountsFiles :
        runMergeCmd.extend(["--counts-file",segmentAlleleCountsFile])
    return self.addTask(runMergeLabel, runMergeCmd, dependencies=dependencies, isForceLocal=True)



def estimateParametersFromAlleleCounts(self, sampleIndex, taskPrefix="", dependencies=None) :
    """
    Estimate variant error parameters from sequencing error count data
    """

    runEstimateLabel=preJoin(taskPrefix,"estimateVariantErrorRates")
    runEstimateCmd=[self.params.estimateVariantErrorRatesBin]
    runEstimateCmd.extend(["--counts-file", self.paths.getAlleleCountsOutputPath(sampleIndex)])
    runEstimateCmd.extend(["--theta-file",self.params.thetaParamFile])
    runEstimateCmd.extend(["--output-file", self.paths.getIndelErrorModelPath(sampleIndex)])
    runEstimateCmd.extend(["--fallback-file",self.params.indelErrorRateDefault])
    return self.addTask(runEstimateLabel, runEstimateCmd, dependencies=dependencies, isForceLocal=True)


class TempSequenceAlleleCountsSegmentFiles :
    def __init__(self) :
        self.counts = []
        self.nonEmptySiteCounts = []


def getSequenceErrorEstimatesForSample(self, estimationIntervals, sampleIndex, taskPrefix="", dependencies=None):
    """
    Count sequencing errors in one sample and use these to estimate sample error parameters
    """

    segmentTasks = set()

    segFiles = TempSequenceAlleleCountsSegmentFiles()

    if self.params.isErrorEstimationFromAllData :
        # get error counts from full data set:
        segmentTasks |= countAllEligibleSequenceEvidence(self, estimationIntervals, sampleIndex, segFiles, taskPrefix, dependencies)
    else :
        # Launch tasks until the required counts are found
        segmentTasks |= countSequenceEvidenceUntilTargetIsReached(self, estimationIntervals, sampleIndex, segFiles, taskPrefix, dependencies)

    # create a checkpoint for all segments:
    completeSegmentsTask = self.addTask(preJoin(taskPrefix,"completedAllGenomeSegments"),dependencies=segmentTasks)

    # merge segment stats:
    mergeCountsTask = mergeSequenceAlleleCounts(self, sampleIndex, segFiles.counts,
                                                taskPrefix=taskPrefix, dependencies=completeSegmentsTask)

    # get error parameters:
    estimateTask = estimateParametersFromAlleleCounts(self, sampleIndex,
                                                      taskPrefix=taskPrefix, dependencies=mergeCountsTask)

    nextStepWait = set()
    nextStepWait.add(estimateTask)
    return nextStepWait


class EstimateSequenceErrorWorkflowForSample(WorkflowRunner) :
    """
    A separate workflow is setup around the error estimation process for each sample so that:
    (1) sequence error estimation can run in parallel for all samples
    """

    def __init__(self, params, paths, estimationIntervals, sampleIndex) :
        self.params = params
        self.paths = paths
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
        sampledErrorEstimationChunkSize = 2 * Megabase
        allDataErrorEstimationChunkSize = 12 * Megabase


    largeChroms = [chrom for chrom in params.chromOrder if params.chromSizes[chrom] >= Constants.errorEstimationMinChromSize]

    if params.isErrorEstimationFromAllData :
        estimationIntervals = [GenomeSegment(*interval) for interval in getChromIntervals(largeChroms, params.chromSizes, Constants.allDataErrorEstimationChunkSize)]

    else :
        estimationIntervals = [GenomeSegment(*interval) for interval in getChromIntervals(largeChroms, params.chromSizes, Constants.sampledErrorEstimationChunkSize)]

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
        workflow=EstimateSequenceErrorWorkflowForSample(self.params, self.paths, estimationIntervals, sampleIndex)
        sampleTasks.add(self.addWorkflowTask(sampleTask, workflow, dependencies=mkDirTask))

    if not self.params.isRetainTempFiles :
        rmTmpCmd = getRmdirCmd() + [tmpErrorEstimationDir]
        self.addTask(preJoin(taskPrefix,"removeTmpDir"), rmTmpCmd, dependencies=sampleTasks, isForceLocal=True)

    nextStepWait = sampleTasks
    return nextStepWait


class EstimateSequenceErrorWorkflow(WorkflowRunner) :
    """
    A separate workflow is setup around the error estimation process for all samples so that:
    (1) execution can be delayed until refcount is complete
    """

    def __init__(self, params, paths) :
        self.params = params
        self.paths = paths

    def workflow(self) :
        getSequenceErrorEstimates(self)

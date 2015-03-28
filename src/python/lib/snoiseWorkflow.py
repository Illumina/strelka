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
Strelka noise estimate workflow
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

from configBuildTimeInfo import workflowVersion
from pyflow import WorkflowRunner
from workflowUtil import checkFile, ensureDir, preJoin, which, \
                         getNextGenomeSegment, getFastaChromOrderSize, bamListCatCmd

from configureUtil import argToBool, getIniSections, dumpIniSections


__version__ = workflowVersion



def runCount(self, taskPrefix="", dependencies=None) :
    """
    count size of fasta chromosomes
    """
    cmd  = "%s %s > %s"  % (self.params.countFastaBin, self.params.referenceFasta, self.paths.getRefCountFile())

    nextStepWait = set()
    nextStepWait.add(self.addTask(preJoin(taskPrefix,"RefCount"), cmd, dependencies=dependencies))

    return nextStepWait



class TempSegmentFiles :
    def __init__(self) :
        self.gvcf = []
        self.bamRealign = []



def callGenomeSegment(self, gseg, segFiles, taskPrefix="", dependencies=None) :

    isFirstSegment = (len(segFiles.gvcf) == 0)

    segStr = str(gseg.id)

    segCmd = [ self.params.snoiseBin ]

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

    segCmd.extend(['-min-qscore','17'])
    segCmd.extend(['-bsnp-ssd-no-mismatch', '0.35'])
    segCmd.extend(['-bsnp-ssd-one-mismatch', '0.6'])
    segCmd.extend(['-min-vexp', '0.25'])

    for bamPath in self.params.bamList :
        segCmd.extend(["-bam-file",bamPath])

    segCmd.extend(["--report-file", self.paths.getTmpSegmentReportPath(gseg.pyflowId)])

    if not isFirstSegment :
        segCmd.append("--skip-vcf-header")

    if self.params.indelCandidates is not None :
        segCmd.extend(['--candidate-indel-input-vcf', self.params.indelCandidates])

     # vcf is written to stdout so we need shell features:
    segCmd = " ".join(segCmd)

    segFiles.gvcf.append(self.paths.getTmpSegmentGvcfPath(segStr))
    segCmd += " | %s -c >| %s" % (self.params.bgzip9Bin, segFiles.gvcf[-1])

    nextStepWait = set()

    setTaskLabel=preJoin(taskPrefix,"callGenomeSegment_"+gseg.pyflowId)
    self.addTask(setTaskLabel,segCmd,dependencies=dependencies,memMb=self.params.callMemMb)
    nextStepWait.add(setTaskLabel)

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
        finishTasks.add(self.addTask(preJoin(taskPrefix,label+"_finalizeVCF"), catCmd, dependencies=completeSegmentsTask))

    finishVcf(segFiles.gvcf, self.paths.getGvcfOutputPath(),"gVCF")

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

    def getTmpSegmentDir(self) :
        return os.path.join(self.params.workDir, "genomeSegment.tmpdir")

    def getTmpSegmentGvcfPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "noise.%s.vcf.gz" % (segStr))

    def getTmpSegmentReportPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "stats.%s.txt" % (segStr))

    def getVariantsDir(self) :
        return self.params.variantsDir

    def getGvcfOutputPath(self) :
        return os.path.join( self.getVariantsDir(), "noise.vcf.gz")

    def getRefCountFile(self) :
        return os.path.join( self.params.workDir, "refCount.txt")



class snoiseWorkflow(WorkflowRunner) :
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

        msg  = "Strelka noise estimation workflow successfully completed.\n\n"
        msg += "\tworkflow version: %s\n" % (__version__)
        return msg



    def workflow(self) :
        self.flowLog("Initiating Strelka noise estimation workflow version: %s" % (__version__))

        callPreReqs = set()
        callPreReqs |= runCount(self)

        self.addWorkflowTask("CallGenome", CallWorkflow(self.params, self.paths), dependencies=callPreReqs)

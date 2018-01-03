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
Shared small variant calling workflow components
"""

import os.path
import sys

# add this path to pull in utils in same directory:
scriptDir=os.path.abspath(os.path.dirname(__file__))
sys.path.append(scriptDir)

# add pyflow path:
sys.path.append(os.path.join(scriptDir,"pyflow"))


from checkChromSet import getTabixChromSet
from pyflow import WorkflowRunner
from sharedWorkflow import getMvCmd
from workflowUtil import checkFile, cleanPyEnv, ensureDir, getFastaChromOrderSize, \
                  getGenomeSegmentGroups, getNextGenomeSegment, preJoin



def getChromIsSkipped(self) :
    """
    Determine subset of chroms from chromOrder which are completely skipped over

    here "skipped" means that not a single base on the chrom is requested for calling or error estimation

    \return The set of chromLabels which are skipped
    """

    chromIsSkipped = set()

    # return empty set when no region selections have been made:
    if ((self.params.genomeRegionList is None) and
        (self.params.callRegionsBed is None)) :
       return chromIsSkipped

    def allChromosomes() :
        """
        Return a set of all chromosomes from the reference/alignments in this analysis
        """
        return set(self.params.chromOrder)

    # first check chromosome coverage of "regions" arguments
    if self.params.genomeRegionList is not None :
        chromIsSkipped = allChromosomes()
        for genomeRegion in self.params.genomeRegionList :
            if genomeRegion["chrom"] in chromIsSkipped :
                chromIsSkipped.remove(genomeRegion["chrom"])

    # further refine coverage based on callRegions BED file
    if self.params.callRegionsBed is not None :
        callRegionsChroms = getTabixChromSet(self.params.tabixBin, self.params.callRegionsBed)
        chromsNotInCallRegions = allChromosomes() - callRegionsChroms

        # Skip the union of:
        # 1. chromosomes skipped already due to region arguments
        # 2. chromosomes skipped due to callRegions bed track
        chromIsSkipped = chromIsSkipped | chromsNotInCallRegions

    # if sequencing error estimation is turned on, make sure estimation targets are not skipped:
    if self.params.isEstimateSequenceError :
        class Constants :
            Megabase = 1000000
            errorEstimationMinChromSize = self.params.errorEstimationMinChromMb * Megabase

        for chrom in self.params.chromSizes :
            if self.params.chromSizes[chrom] < Constants.errorEstimationMinChromSize : continue
            if chrom in chromIsSkipped :
                chromIsSkipped.remove(chrom)

    return chromIsSkipped



class StrelkaSharedCallWorkflow(WorkflowRunner) :

    def __init__(self, params) :
        self.params = params

    def concatIndexBgzipFile(self, taskPrefix, dependencies, inputList, output, label, fileType) :
        """
        Internal helper function
        @param inputList files to be concatenated (in order), already bgzipped
        @param output output filename
        @param label used for error task id
        @param fileType provided to tabix
        """
        assert(len(inputList) > 0)

        if len(inputList) > 1 :
            catCmd = [self.params.bgcatBin,"-o",output]
            catCmd.extend(inputList)
        else :
            catCmd = getMvCmd() + [inputList[0], output]

        indexCmd = [self.params.tabixBin,"-p", fileType, output]
        catTask = self.addTask(preJoin(taskPrefix,label+"_concat_"+fileType), catCmd,
                               dependencies=dependencies, isForceLocal=True)
        return self.addTask(preJoin(taskPrefix,label+"_index_"+fileType), indexCmd,
                                    dependencies=catTask, isForceLocal=True)



    def concatIndexVcf(self, taskPrefix, dependencies, inputList, output, label) :
        """
        Concatenate bgzipped vcf segments
        @param inputList files to be concatenated (in order), already bgzipped
        @param output output filename
        @param label used for error task id
        """
        assert(len(inputList) > 0)
        return self.concatIndexBgzipFile(taskPrefix, dependencies, inputList, output, label, "vcf")



    def concatIndexBed(self, taskPrefix, dependencies, inputList, output, label) :
        """
        Concatenate bgzipped bed segments
        @param inputList files to be concatenated (in order), already bgzipped
        @param output output filename
        @param label used for error task id
        """
        assert(len(inputList) > 0)
        return self.concatIndexBgzipFile(taskPrefix, dependencies, inputList, output, label, "bed")



    def mergeRunStats(self, taskPrefix, dependencies, runStatsLogPaths) :
        """
        merge run stats:
        """
        runStatsMergeLabel=preJoin(taskPrefix,"mergeRunStats")
        runStatsMergeCmd=[self.params.statsMergeBin]
        for statsFile in runStatsLogPaths :
            runStatsMergeCmd.extend(["--stats-file",statsFile])
        runStatsMergeCmd.extend(["--output-file",self.paths.getRunStatsPath()])
        runStatsMergeCmd.extend(["--report-file",self.paths.getRunStatsReportPath()])
        return self.addTask(runStatsMergeLabel, runStatsMergeCmd, dependencies=dependencies, isForceLocal=True)



    def appendCommonGenomeSegmentCommandOptions(self, gsegGroup, genomeSegmentLabel, segCmd) :
        """
        Append cmdline components to the cmd arg list that are common to multiple strelka workflows

        @param genomeSegmentLabel Genome segment identifier string
        @param segCmd Command argument list
        """

        for gseg in gsegGroup :
            segCmd.extend(["--region", gseg.bamRegion])

        segCmd.extend(["--ref", self.params.referenceFasta ])
        segCmd.extend(["--max-indel-size", self.params.maxIndelSize] )

        if self.params.indelErrorModelName is not None :
            segCmd.extend(['--indel-error-model-name',self.params.indelErrorModelName])

        if self.params.isReportEVSFeatures :
            segCmd.append("--report-evs-features")

        def addListCmdOption(optList,arg) :
            if optList is None : return
            for val in optList :
                segCmd.extend([arg, val])

        addListCmdOption(self.params.indelCandidatesList, '--candidate-indel-input-vcf')
        addListCmdOption(self.params.forcedGTList, '--force-output-vcf')

        if self.params.callRegionsBed is not None :
            segCmd.extend(['--call-regions-bed', self.params.callRegionsBed])

        if self.params.isWriteRealignedBam :
            segCmd.extend(["--realigned-output-prefix", self.paths.getTmpUnsortRealignBamPrefix(genomeSegmentLabel)])

        if self.params.extraVariantCallerArguments is not None :
            for arg in self.params.extraVariantCallerArguments.strip().split() :
                segCmd.append(arg)


    def filterUncalledChromosomeSegments(self, genomeSegmentIterator):
        """
        Given a genome segment iterator, filter the genome segments to remove any chromosomes which are
        absent from callRegions bed file, if such a file has been specified, otherwise just pass the
        iterator through.
        """

        isCallRegionBed = (self.params.callRegionsBed is not None)

        for gseg in genomeSegmentIterator :
            if isCallRegionBed :
                if gseg.chromLabel in self.params.chromIsSkipped : continue

            yield gseg


    def getStrelkaGenomeSegmentGroupIterator(self, contigsExcludedFromGrouping = None) :
        """
        setup genome segment iteration for germline and somatic calling,
        including  clumping together small segments into groups.
        """
        genomeSegmentIterator = self.filterUncalledChromosomeSegments(getNextGenomeSegment(self.params))
        return getGenomeSegmentGroups(genomeSegmentIterator, contigsExcludedFromGrouping)



class SharedPathInfo(object):
    """
    object to centralize shared workflow path names
    """

    def __init__(self, params) :
        self.params = params

    def getChromDepth(self) :
        return os.path.join(self.params.workDir,"chromDepth.tsv")

    def getTmpSegmentDir(self) :
        return os.path.join(self.params.workDir, "genomeSegment.tmpdir")

    def getTmpErrorEstimationDir(self) :
        return os.path.join(self.params.workDir, "errorEstimation.tmpdir")

    def getTmpRunStatsPath(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "runStats.%s.xml" % (segStr))

    def getRunStatsPath(self) :
        return os.path.join(self.params.statsDir,"runStats.xml")

    def getRunStatsReportPath(self) :
        return os.path.join(self.params.statsDir,"runStats.tsv")

    def getTmpUnsortRealignBamPrefix(self, segStr) :
        return os.path.join( self.getTmpSegmentDir(), "segment.%s.unsorted.realigned." % (segStr))



class StrelkaSharedWorkflow(WorkflowRunner) :
    """
    small variant calling workflow base
    """

    def __init__(self,params,PathInfoType) :

        cleanPyEnv()

        self.params=params

        # make sure run directory is setup:
        self.params.runDir=os.path.abspath(self.params.runDir)
        ensureDir(self.params.runDir)

        # everything that's not intended to be a final result should dump directories/files in workDir
        self.params.workDir=os.path.join(self.params.runDir,"workspace")
        ensureDir(self.params.workDir)

        # all finalized pretty results get transferred to resultsDir
        self.params.resultsDir=os.path.join(self.params.runDir,"results")
        ensureDir(self.params.resultsDir)
        self.params.variantsDir=os.path.join(self.params.resultsDir,"variants")
        ensureDir(self.params.variantsDir)

        # timings and other stats go into statsDir
        self.params.statsDir=os.path.join(self.params.resultsDir,"stats")
        ensureDir(self.params.statsDir)

        self.paths = PathInfoType(self.params)

        referenceFastaIndex=self.params.referenceFasta+".fai"

        if self.params.referenceFasta is None:
            raise Exception("No reference fasta defined.")
        else:
            checkFile(self.params.referenceFasta,"reference fasta")
            checkFile(referenceFastaIndex,"reference fasta index")

        # read fasta index
        (self.params.chromOrder,self.params.chromSizes) = getFastaChromOrderSize(referenceFastaIndex)

        # determine subset of chroms where we can skip calling entirely
        self.params.chromIsSkipped = getChromIsSkipped(self)

        self.params.isHighDepthFilter = (not (self.params.isExome or self.params.isRNA))


    def setCallMemMb(self) :
        # set default mem requirements:
        if self.params.callMemMbOverride is not None :
            self.params.callMemMb = self.params.callMemMbOverride
        else :
            if self.getRunMode() == "sge" :
                self.params.callMemMb = self.params.callSGEMemMb
            else :
                self.params.callMemMb = self.params.callLocalMemMb

#!/usr/bin/env python
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
This script configures the Strelka somatic small variant calling workflow
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)
workflowDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))

version="@STARKA_FULL_VERSION@"


sys.path.append(workflowDir)

from strelkaOptions import StrelkaWorkflowOptionsBase
from configureUtil import assertOptionExists, OptParseException, validateFixExistingDirArg, validateFixExistingFileArg
from makeRunScript import makeRunScript
from strelkaWorkflow import StrelkaWorkflow
from workflowUtil import ensureDir, isValidSampleId, parseGenomeRegion
from checkChromSet import checkChromSet



class StrelkaWorkflowOptions(StrelkaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Strelka somatic small variant calling pipeline.
You must specify BAM file(s) for a pair of samples.
""" % (version)

    validAlignerModes = ["bwa","isaac"]


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--normalBam", type="string",dest="normalBamList",metavar="FILE", action="append",
                         help="Normal sample BAM file. [required] (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                          help="Tumor sample BAM file. [required] (no default)")
#         group.add_option("--aligner", type="string",dest="alignerMode",metavar="ALIGNER",
#                          help="Aligner type. Accepted option are {%s} [required] (no default)" % (",".join(['%s' % (x) for x in self.validAlignerModes])))
        #group.add_option("--exome", dest="isExome", action="store_true",
        #                 help="Set options for WES input: turn off depth filters")
        #group.add_option("--rna", dest="isRNA", action="store_true",
        #                 help="Set options for RNA-Seq input: turn off depth filters and don't treat "
        #                      "anomalous reads as SV evidence when the proper-pair bit is set.")
        group.add_option("--referenceFasta",type="string",dest="referenceFasta",metavar="FILE",
                         help="samtools-indexed reference fasta file [required] (default: %default)")

        MantaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def addExtendedGroupOptions(self,group) :
        group.add_option("--scanSizeMb",type="int",dest="scanSizeMb",metavar="scanSizeMb",
                         help="Maximum sequence region size (in Mb) scanned by each task during "
                         "SV locus graph generation. (default: %default)")
        group.add_option("--region",type="string",dest="regionStrList",metavar="samtoolsRegion", action="append",
                         help="Limit the SV analysis to a region of the genome for debugging purposes. "
                              "If this argument is provided multiple times all specified regions will "
                              "be analyzed together. All regions must be non-overlapping to get a "
                              "meaningful result. Examples: '--region chr20' (whole chromosome), "
                              "'--region chr2:100-2000 --region chr3:2500-3000' (two translocation regions)'")

        StrelkaWorkflowOptionsBase.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StrelkaWorkflowOptionsBase.getOptionDefaults(self)
        defaults.update({
            'alignerMode' : "isaac",
            'runDir' : 'StrelkaWorkflow',
            "minTier2Mapq" : 5,
            'scanSizeMb' : 12
            })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        def checkForBamIndex(bamFile):
            baiFile=bamFile + ".bai"
            if not os.path.isfile(baiFile) :
                raise OptParseException("Can't find expected BAM index file: '%s'" % (baiFile))

        def groomBamList(bamList, sampleLabel):
            if bamList is None : return
            for (index,bamFile) in enumerate(bamList) :
                bamList[index]=validateFixExistingFileArg(bamFile,"%s BAM file" % (sampleLabel))
                checkForBamIndex(bamList[index])

        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        # check alignerMode:
        if options.alignerMode is not None :
            options.alignerMode = options.alignerMode.lower()
            if options.alignerMode not in self.validAlignerModes :
                raise OptParseException("Invalid aligner mode: '%s'" % options.alignerMode)

        options.referenceFasta=validateFixExistingFileArg(options.referenceFasta,"reference")

        # check for reference fasta index file:
        if options.referenceFasta is not None :
            faiFile=options.referenceFasta + ".fai"
            if not os.path.isfile(faiFile) :
                raise OptParseException("Can't find expected fasta index file: '%s'" % (faiFile))

        if (options.regionStrList is None) or (len(options.regionStrList) == 0) :
            options.genomeRegionList = None
        else :
            options.genomeRegionList = [parseGenomeRegion(r) for r in options.regionStrList]

        MantaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)



    def validateOptionExistence(self,options) :

        # note that we inherit a multi-bam capable infrastructure from manta, but then restrict usage
        # to one bam from each sample (hopefully temporarily)
        #
        def checkBamList(bamList, label) :
            if (bamList is None) or (len(options.bamList) == 0) :
                raise OptParseException("No %s sample BAM files specified" % (label))

            if len(options.bamList) > 1 :
                raise OptParseException("More than one %s sample BAM files specified" % (label))

        checkBamList(options.normalBamList, "normal")
        checkBamList(options.tumorBamList, "tumor")

        assertOptionExists(options.alignerMode,"aligner mode")
        assertOptionExists(options.referenceFasta,"reference fasta file")

        StrelkaWorkflowOptionsBase.validateOptionExistence(self,options)

        # check that the reference and all bams are using the same
        # set of chromosomes:
        bamList=[]
        bamLabels=[]

        def appendBams(inputBamList,inputLabel) :
            if inputBamList is None : return
            for inputBamFile in inputBamList :
                bamList.append(inputBamFile)
                bamLabels.append(inputLabel)

        appendBams(options.normalBamList,"Normal")
        appendBams(options.tumorBamList,"Tumor")

        checkChromSet(options.samtoolsBin,
                      options.referenceFasta,
                      bamList,
                      bamLabels,
                      isReferenceLocked=True)

        # check for repeated bam entries:
        #
        bamSet=set()
        for bamFile in bamList :
            if bamFile in bamSet :
                raise OptParseException("Repeated input BAM file: %s" % (bamFile))
            bamSet.add(bamFile)



def main() :

    primarySectionName="strelka"
    options,iniSections=StrelkaWorkflowOptions().getRunOptions(primarySectionName, version=version)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    StrelkaWorkflow(options,iniSections)

    # generate runscript:
    #
    ensureDir(options.runDir)
    scriptFile=os.path.join(options.runDir,"runWorkflow.py")

    makeRunScript(scriptFile,os.path.join(workflowDir,"strelkaWorkflow.py"),"StrelkaWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (scriptFile))


if __name__ == "__main__" :
    main()


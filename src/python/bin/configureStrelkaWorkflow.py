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

from starkaOptions import StarkaWorkflowOptionsBase
from configureUtil import assertOptionExists, groomBamList, OptParseException
from makeRunScript import makeRunScript
from strelkaWorkflow import StrelkaWorkflow
from workflowUtil import ensureDir, isValidSampleId, parseGenomeRegion
from checkChromSet import checkChromSet



class StrelkaWorkflowOptions(StarkaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Strelka somatic small variant calling pipeline.
You must specify BAM file(s) for a pair of samples.
""" % (version)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--normalBam", type="string",dest="normalBamList",metavar="FILE", action="append",
                         help="Normal sample BAM file. [required] (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                          help="Tumor sample BAM file. [required] (no default)")
        #group.add_option("--exome", dest="isExome", action="store_true",
        #                 help="Set options for WES input: turn off depth filters")

        StarkaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StarkaWorkflowOptionsBase.getOptionDefaults(self)
        defaults.update({
            'runDir' : 'StrelkaWorkflow',
            "minTier2Mapq" : 5,
            })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        StarkaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)



    def validateOptionExistence(self,options) :

        # note that we inherit a multi-bam capable infrastructure from manta, but then restrict usage
        # to one bam from each sample (hopefully temporarily)
        #
        def checkBamList(bamList, label) :
            if (bamList is None) or (len(bamList) == 0) :
                raise OptParseException("No %s sample BAM files specified" % (label))

            if len(bamList) > 1 :
                raise OptParseException("More than one %s sample BAM files specified" % (label))

        checkBamList(options.normalBamList, "normal")
        checkBamList(options.tumorBamList, "tumor")

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

        StarkaWorkflowOptionsBase.validateOptionExistence(self,options)




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


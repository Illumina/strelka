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
templateConfigDir=os.path.abspath(os.path.join(scriptDir,'@THIS_RELATIVE_CONFIGDIR@'))

sys.path.append(workflowDir)

from configBuildTimeInfo import workflowVersion
from starkaOptions import StarkaWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, OptParseException, joinFile, checkTabixListOption
from makeRunScript import makeRunScript
from strelkaWorkflow import StrelkaWorkflow
from workflowUtil import ensureDir



class StrelkaWorkflowOptions(StarkaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Strelka somatic small variant calling pipeline.
You must specify BAM/CRAM file(s) for a pair of samples.
""" % (workflowVersion)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--normalBam", type="string",dest="normalBamList",metavar="FILE", action="append",
                         help="Normal sample BAM or CRAM file. [required] (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                          help="Tumor sample BAM or CRAM file. [required] (no default)")
        group.add_option("--noiseVcf", type="string",dest="noiseVcfList",metavar="FILE", action="append",
                          help="Noise vcf file (submit argument multiple times for more than one file)")
        group.add_option("--isWriteCallableRegion", action="store_true",
                         help="Write out a bed file describing somatic callable regions of thedupliates genome")
        group.add_option("--allowDuplicateReadIDs", action="store_true",
                         help="Allow duplicate read IDs in the same segment. This works around a bug in BWA-mem 0.7.12"
                              " which doesn't correctly mark secondary alignments.")

        StarkaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StarkaWorkflowOptionsBase.getOptionDefaults(self)

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'StrelkaWorkflow',
            "minTier2Mapq" : 0,
            'variantScoringModelFile' : joinFile(configDir,'somaticVariantScoringModels.json'),
            'indelErrorModelsFile' : joinFile(configDir,'indelErrorModels.json'),
            'indelErrorModelName': 'new',
            'isWriteCallableRegion' : False
            })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        StarkaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)
        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        checkTabixListOption(options.noiseVcfList,"noise vcf")



    def validateOptionExistence(self,options) :

        StarkaWorkflowOptionsBase.validateOptionExistence(self,options)

        bcheck = BamSetChecker()

        def singleAppender(bamList,label):
            if len(bamList) > 1 :
                raise OptParseException("More than one %s sample BAM/CRAM files specified" % (label))
            bcheck.appendBams(bamList,label)

        singleAppender(options.normalBamList,"Normal")
        singleAppender(options.tumorBamList,"Tumor")
        bcheck.check(options.samtoolsBin,
                     options.referenceFasta)

def main() :

    primarySectionName="strelka"
    options,iniSections=StrelkaWorkflowOptions().getRunOptions(primarySectionName,
                                                               version=workflowVersion)

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

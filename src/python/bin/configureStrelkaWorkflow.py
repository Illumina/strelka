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

version="@STARKA_FULL_VERSION@"


sys.path.append(workflowDir)

from starkaOptions import StarkaWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, OptParseException, joinFile
from makeRunScript import makeRunScript
from strelkaWorkflow import StrelkaWorkflow
from workflowUtil import ensureDir



class StrelkaWorkflowOptions(StarkaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Strelka somatic small variant calling pipeline.
You must specify BAM/CRAM file(s) for a pair of samples.
""" % (version)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--normalBam", type="string",dest="normalBamList",metavar="FILE", action="append",
                         help="Normal sample BAM or CRAM file. [required] (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                          help="Tumor sample BAM or CRAM file. [required] (no default)")
        group.add_option("--noiseVcf", type="string",dest="noiseVcfList",metavar="FILE", action="append",
                          help="Noise vcf file (submit argument multiple times for more than one file)")
#        group.add_option("--scoringModel", type="string",dest="scoringModel",metavar="FILE", action="append",
#                          help="Json file specifying filtering and indel model parameters")
        #group.add_option("--exome", dest="isExome", action="store_true",
        #                 help="Set options for WES input: turn off depth filters")

        StarkaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StarkaWorkflowOptionsBase.getOptionDefaults(self)

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'StrelkaWorkflow',
            "minTier2Mapq" : 5,
            'scoringModelFile' : joinFile(configDir,'indel_models.json')
            })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        StarkaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)
        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        if options.noiseVcfList is not None :
            for vcfname in options.noiseVcfList :
                tabixFile = vcfname + ".tbi"
                if os.path.isfile(tabixFile) : return
                raise OptParseException("Can't find expected noise vcf index file: '%s'" % (tabixFile))



    def validateOptionExistence(self,options) :

        StarkaWorkflowOptionsBase.validateOptionExistence(self,options)

        if options.userConfigPath is None :
            raise OptParseException("A config file must be specified for strelka workflow")

        bcheck = BamSetChecker()
        bcheck.appendBams(options.normalBamList,"Normal")
        bcheck.appendBams(options.tumorBamList,"Tumor")
        bcheck.check(options.samtoolsBin,
                     options.referenceFasta)



def main() :

    primarySectionName="strelka"

    # change the config help text for strelka:
    configHelp="Provide strelka workflow configuration file, see template configuration files in: '%s' (required)" % (templateConfigDir)

    options,iniSections=StrelkaWorkflowOptions().getRunOptions(primarySectionName,
                                                               version=version,
                                                               configHelp=configHelp)

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


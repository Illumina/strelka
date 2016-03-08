#!/usr/bin/env python
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
from configureUtil import BamSetChecker, groomBamList, OptParseException, joinFile, \
                            checkFixTabixListOption, validateFixExistingFileArg
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
                         help="Normal sample BAM or CRAM file. (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                         help="Tumor sample BAM or CRAM file. [required] (no default)")
        group.add_option("--isWriteCallableRegion", action="store_true",
                         help="Write out a bed file describing somatic callable regions of thedupliates genome")

        StarkaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)

    def addExtendedGroupOptions(self,group) :
        group.add_option("--somaticSnvScoringModelFile", type="string", dest="somaticSnvScoringModelFile", metavar="FILE",
                         help="Provide a custom EVS model file for somatic SNVs (default: %default)")
        group.add_option("--enableSomaticIndelScoring", action="store_true", dest="isSomaticIndelEmpiricalScoring",
                         help="Enable empirical variant scoring for somatic indels")
        group.add_option("--somaticIndelScoringModelFile", type="string", dest="somaticIndelScoringModelFile", metavar="FILE",
                         help="Provide a custom EVS model file for somatic Indels (default: %default)")
        group.add_option("--noiseVcf", type="string",dest="noiseVcfList",metavar="FILE", action="append",
                         help="Noise vcf file (submit argument multiple times for more than one file)")


        StarkaWorkflowOptionsBase.addExtendedGroupOptions(self,group)

    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StarkaWorkflowOptionsBase.getOptionDefaults(self)

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'StrelkaWorkflow',
            "minTier2Mapq" : 0,
            "isSomaticIndelEmpiricalScoring" : False,
            'somaticSnvScoringModelFile' : joinFile(configDir,'somaticVariantScoringModels.json'),
            'somaticIndelScoringModelFile' : joinFile(configDir,'somaticVariantScoringModels.json'),
            'indelErrorModelsFile' : joinFile(configDir,'indelErrorModels.json'),
            'indelErrorModelName': 'Binom',
            'isWriteCallableRegion' : False,
            'noiseVcfList' : None
            })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        StarkaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)
        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        checkFixTabixListOption(options.noiseVcfList,"noise vcf")

        options.somaticSnvScoringModelFile=validateFixExistingFileArg(options.somaticSnvScoringModelFile,"Somatic SNV empirical scoring file")
        options.somaticIndelScoringModelFile=validateFixExistingFileArg(options.somaticIndelScoringModelFile,"Somatic indel empirical scoring file")


    def validateOptionExistence(self,options) :

        StarkaWorkflowOptionsBase.validateOptionExistence(self,options)

        def checkRequired(bamList,label):
            if (bamList is None) or (len(bamList) == 0) :
                raise OptParseException("No %s sample BAM/CRAM files specified" % (label))

        checkRequired(options.tumorBamList,"tumor")

        bcheck = BamSetChecker()

        def singleAppender(bamList,label):
            if bamList is None : return
            if len(bamList) > 1 :
                raise OptParseException("More than one %s sample BAM/CRAM files specified" % (label))
            bcheck.appendBams(bamList,label)

        singleAppender(options.normalBamList,"normal")
        singleAppender(options.tumorBamList,"tumor")
        bcheck.check(options.htsfileBin,
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

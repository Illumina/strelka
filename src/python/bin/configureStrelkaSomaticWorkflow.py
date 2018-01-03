#!/usr/bin/env python2
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
This script configures the Strelka somatic small variant calling workflow
"""

import os,sys

if sys.version_info >= (3,0):
    import platform
    raise Exception("Strelka does not currently support python3 (version %s detected)" % (platform.python_version()))

if sys.version_info < (2,6):
    import platform
    raise Exception("Strelka requires python2 version 2.6+ (version %s detected)" % (platform.python_version()))


scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)
workflowDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))
templateConfigDir=os.path.abspath(os.path.join(scriptDir,'@THIS_RELATIVE_CONFIGDIR@'))

sys.path.append(workflowDir)

from configBuildTimeInfo import workflowVersion
from strelkaSharedOptions import StrelkaSharedWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, OptParseException, joinFile, \
                            checkFixTabixListOption, validateFixExistingFileArg
from makeRunScript import makeRunScript
from strelkaSomaticWorkflow import StrelkaSomaticWorkflow
from workflowUtil import ensureDir, exeFile



class StrelkaSomaticWorkflowOptions(StrelkaSharedWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures Strelka somatic small variant calling.
You must specify an alignment file (BAM or CRAM) for each sample of a matched tumor-normal pair.
""" % (workflowVersion)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--normalBam", type="string",dest="normalBamList",metavar="FILE", action="append",
                         help="Normal sample BAM or CRAM file. (no default)")
        group.add_option("--tumorBam","--tumourBam", type="string",dest="tumorBamList",metavar="FILE", action="append",
                         help="Tumor sample BAM or CRAM file. [required] (no default)")
        group.add_option("--outputCallableRegions", dest="isOutputCallableRegions", action="store_true",
                         help="Output a bed file describing somatic callable regions of the genome")

        StrelkaSharedWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def addExtendedGroupOptions(self,group) :
        group.add_option("--noiseVcf", type="string",dest="noiseVcfList",metavar="FILE", action="append",
                         help="Noise vcf file (submit argument multiple times for more than one file)")

        StrelkaSharedWorkflowOptionsBase.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StrelkaSharedWorkflowOptionsBase.getOptionDefaults(self)

        libexecDir=defaults["libexecDir"]

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'StrelkaSomaticWorkflow',
            'strelkaSomaticBin' : joinFile(libexecDir,exeFile("strelka2")),
            'minTier2Mapq' : 0,
            'snvScoringModelFile' : joinFile(configDir,'somaticSNVScoringModels.json'),
            'indelScoringModelFile' : joinFile(configDir,'somaticIndelScoringModels.json'),
            'isOutputCallableRegions' : False,
            'noiseVcfList' : None
            })
        return defaults


    def validateAndSanitizeOptions(self,options) :

        StrelkaSharedWorkflowOptionsBase.validateAndSanitizeOptions(self,options)

        checkFixTabixListOption(options.noiseVcfList,"noise vcf")

        groomBamList(options.normalBamList,"normal sample")
        groomBamList(options.tumorBamList, "tumor sample")

        def checkRequired(bamList,label):
            if (bamList is None) or (len(bamList) == 0) :
                raise OptParseException("No %s sample BAM/CRAM files specified" % (label))

        checkRequired(options.tumorBamList,"tumor")

        bamSetChecker = BamSetChecker()

        def singleAppender(bamList,label):
            if bamList is None : return
            if len(bamList) > 1 :
                raise OptParseException("More than one %s sample BAM/CRAM files specified" % (label))
            bamSetChecker.appendBams(bamList,label)

        singleAppender(options.normalBamList,"normal")
        singleAppender(options.tumorBamList,"tumor")
        bamSetChecker.check(options.htsfileBin,
                     options.referenceFasta)



def main() :

    primarySectionName="StrelkaSomatic"
    options,iniSections=StrelkaSomaticWorkflowOptions().getRunOptions(primarySectionName,
                                                               version=workflowVersion)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    StrelkaSomaticWorkflow(options)

    # generate runscript:
    #
    ensureDir(options.runDir)
    workflowScriptPath = os.path.join(options.runDir, options.workflowScriptName)

    makeRunScript(workflowScriptPath,os.path.join(workflowDir,"strelkaSomaticWorkflow.py"),"StrelkaSomaticWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (workflowScriptPath))


if __name__ == "__main__" :
    main()

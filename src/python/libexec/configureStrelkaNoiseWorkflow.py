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
This script configures the strelka noise estimation workflow
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)
workflowDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))

version="@STARKA_FULL_VERSION@"


sys.path.append(workflowDir)

from starkaOptions import StarkaWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, joinFile, OptParseException
from makeRunScript import makeRunScript
from snoiseWorkflow import snoiseWorkflow
from workflowUtil import ensureDir



class snoiseWorkflowOptions(StarkaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Strelka noise estimation pipeline.
You must specify a BAM or CRAM file.
""" % (version)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--bam", type="string",dest="bamList",metavar="FILE", action="append",
                         help="Sample BAM or CRAM file. [required] (no default)")

        StarkaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)



    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StarkaWorkflowOptionsBase.getOptionDefaults(self)

        libexecDir=defaults["libexecDir"]

        defaults.update({
            'runDir' : 'StrelkaNoiseWorkflow',
            'bgcatBin' : joinFile(libexecDir,"bgzf_cat"),
            'bgzip9Bin' : joinFile(libexecDir,"bgzip9"),
            'snoiseBin' : joinFile(libexecDir,"strelkaNoiseExtractor")
            })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        StarkaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)
        groomBamList(options.bamList,"input")



    def validateOptionExistence(self,options) :

        StarkaWorkflowOptionsBase.validateOptionExistence(self,options)
        bcheck = BamSetChecker()
        bcheck.appendBams(options.bamList,"Input")
        bcheck.check(options.samtoolsBin,
                     options.referenceFasta)



def main() :

    primarySectionName="snoise"
    options,iniSections=snoiseWorkflowOptions().getRunOptions(primarySectionName, version=version)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    snoiseWorkflow(options,iniSections)

    # generate runscript:
    #
    ensureDir(options.runDir)
    scriptFile=os.path.join(options.runDir,"runWorkflow.py")

    makeRunScript(scriptFile,os.path.join(workflowDir,"snoiseWorkflow.py"),"snoiseWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (scriptFile))


if __name__ == "__main__" :
    main()


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
This script configures the starling small variant calling workflow
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
from starlingWorkflow import StarlingWorkflow
from workflowUtil import ensureDir



class StarlingWorkflowOptions(StarkaWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Starling small variant calling pipeline.
You must specify a BAM or CRAM file.
""" % (version)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--bam", type="string",dest="bamList",metavar="FILE", action="append",
                         help="Sample BAM or CRAM file. [required] (no default)")
        group.add_option("--minorAllele", type="string", metavar="FILE",
                         help="Provide minor allele bed file. Must be tabix indexed. (no default)")

        StarkaWorkflowOptionsBase.addWorkflowGroupOptions(self,group)



    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StarkaWorkflowOptionsBase.getOptionDefaults(self)

        libexecDir=defaults["libexecDir"]

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'StarlingWorkflow',
            'bgcatBin' : joinFile(libexecDir,"bgzf_cat"),
            'bgzip9Bin' : joinFile(libexecDir,"bgzip9"),
            'vqsrModel' : "QScoreHpolmodel",
            'vqsrModelFile' : joinFile(configDir,'model.json'),
            'scoringModelFile' : joinFile(configDir,'indel_models.json')
            })
        return defaults



    def validateAndSanitizeExistingOptions(self,options) :

        StarkaWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)
        groomBamList(options.bamList,"input")



    def validateOptionExistence(self,options) :

        StarkaWorkflowOptionsBase.validateOptionExistence(self,options)

        if options.minorAllele is not None :
            alleleTabixFile = options.minorAllele + ".tbi"
            if not os.path.isfile(alleleTabixFile) :
                raise OptParseException("Can't find expected minor allele index file: '%s'" % (alleleTabixFile))

        bcheck = BamSetChecker()
        bcheck.appendBams(options.bamList,"Input")
        bcheck.check(options.samtoolsBin,
                     options.referenceFasta)



def main() :

    primarySectionName="starling"
    options,iniSections=StarlingWorkflowOptions().getRunOptions(primarySectionName, version=version)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    StarlingWorkflow(options,iniSections)

    # generate runscript:
    #
    ensureDir(options.runDir)
    scriptFile=os.path.join(options.runDir,"runWorkflow.py")

    makeRunScript(scriptFile,os.path.join(workflowDir,"starlingWorkflow.py"),"StarlingWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (scriptFile))


if __name__ == "__main__" :
    main()


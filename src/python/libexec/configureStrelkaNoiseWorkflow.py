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
This script configures the strelka somatic variant noise estimation workflow
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)
workflowDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))


sys.path.append(workflowDir)

from configBuildTimeInfo import workflowVersion
from strelkaSharedOptions import StrelkaSharedWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, joinFile, OptParseException
from makeRunScript import makeRunScript
from snoiseWorkflow import snoiseWorkflow
from workflowUtil import ensureDir



class snoiseWorkflowOptions(StrelkaSharedWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures Strelka noise estimation.
You must specify a BAM or CRAM file.
""" % (workflowVersion)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--bam", type="string",dest="bamList",metavar="FILE", action="append",
                         help="Sample BAM or CRAM file. [required] (no default)")

        StrelkaSharedWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StrelkaSharedWorkflowOptionsBase.getOptionDefaults(self)

        libexecDir=defaults["libexecDir"]

        defaults.update({
            'runDir' : 'StrelkaNoiseWorkflow',
            'workflowScriptName' : 'runWorkflow.py',
            'bgcatBin' : joinFile(libexecDir,"bgzf_cat"),
            'bgzip9Bin' : joinFile(libexecDir,"bgzip9"),
            'snoiseBin' : joinFile(libexecDir,"strelkaNoiseExtractor")
            })
        return defaults


    def validateAndSanitizeOptions(self,options) :

        StrelkaSharedWorkflowOptionsBase.validateAndSanitizeOptions(self,options)

        groomBamList(options.bamList,"input")

        bamSetChecker = BamSetChecker()
        bamSetChecker.appendBams(options.bamList,"Input")
        bamSetChecker.check(options.samtoolsBin,
                            options.referenceFasta)



def main() :

    if (sys.version_info[0] != 2):
        notefp=sys.stdout
        notefp.write("""Failed to create workflow run script.\nPyflow only supports python2. Detected python %s on the system.\n""" % sys.version_info[0])
        return

    primarySectionName="snoise"
    options,iniSections=snoiseWorkflowOptions().getRunOptions(primarySectionName, version=workflowVersion)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    snoiseWorkflow(options)

    # generate runscript:
    #
    ensureDir(options.runDir)
    workflowScriptPath = os.path.join(options.runDir, options.workflowScriptName)

    makeRunScript(workflowScriptPath,os.path.join(workflowDir,"snoiseWorkflow.py"),"snoiseWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (workflowScriptPath))


if __name__ == "__main__" :
    main()

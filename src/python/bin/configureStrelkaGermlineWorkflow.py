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
This script configures the strelka germline small variant calling workflow
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)
workflowDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))

sys.path.append(workflowDir)

from configBuildTimeInfo import workflowVersion
from strelkaSharedOptions import StrelkaSharedWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, joinFile, OptParseException, checkFixTabixIndexedFileOption
from makeRunScript import makeRunScript
from strelkaGermlineWorkflow import StrelkaGermlineWorkflow
from workflowUtil import ensureDir



class StrelkaGermlineWorkflowOptions(StrelkaSharedWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures Strelka germline small variant calling.
You must specify an alignment file (BAM or CRAM) for at least one sample.
""" % (workflowVersion)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--bam", type="string",dest="bamList",metavar="FILE", action="append",
                         help="Sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [required] (no default)")
        group.add_option("--ploidy", type="string", dest="ploidyFilename", metavar="FILE",
                         help="Provide ploidy file in VCF. The VCF should include one sample column per input sample labeled with the same sample names found in the input BAM/CRAM RG header sections."
                              " Ploidy should be provided in records using the FORMAT/CN field, which are interpreted to span the range [POS+1, INFO/END]. Any CN value besides 1 or 0 will be treated as 2."
                              " File must be tabix indexed. (no default)")
        group.add_option("--noCompress", type="string", dest="noCompressBed", metavar="FILE",
                         help="Provide BED file of regions where gVCF block compress is disallowed. File must be bgzip-compressed/tabix-indexed. (no default)")
        group.add_option("--targetRegions", type="string", dest="targetRegionsBed", metavar="FILE",
                         help="Provide BED file of regions to allow variant calls. Calls outside these ares are filtered "
                         "as OffTarget. File must be bgzip-compressed/tabix-indexed. (no default)")
        group.add_option("--callContinuousVf", type="string", dest="callContinuousVf", metavar="CHROM", action="append",
                         help="Call variants on CHROM without a ploidy prior assumption, issuing calls with continuous variant frequencies (no default)")
        group.add_option("--rna", dest="isRNA", action="store_true",
                         help="Set options for RNA-Seq input.")

        StrelkaSharedWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StrelkaSharedWorkflowOptionsBase.getOptionDefaults(self)

        libexecDir=defaults["libexecDir"]

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'StrelkaGermlineWorkflow',
            'bgzip9Bin' : joinFile(libexecDir,"bgzip9"),
            'configDir' : configDir,
            'snvScoringModelFile' : joinFile(configDir,'germlineSNVScoringModels.json'),
            'indelScoringModelFile' : joinFile(configDir,'germlineIndelScoringModels.json'),
            'callContinuousVf' : []
            })
        return defaults


    def validateAndSanitizeExistingOptions(self,options) :

        StrelkaSharedWorkflowOptionsBase.validateAndSanitizeExistingOptions(self,options)
        groomBamList(options.bamList,"input")

        options.ploidyFilename = checkFixTabixIndexedFileOption(options.ploidyFilename,"ploidy file")
        options.noCompressBed = checkFixTabixIndexedFileOption(options.noCompressBed,"no-compress bed")
        options.targetRegionsBed = checkFixTabixIndexedFileOption(options.targetRegionsBed,"targeted-regions bed")


    def validateOptionExistence(self,options) :

        StrelkaSharedWorkflowOptionsBase.validateOptionExistence(self,options)

        def safeLen(x) :
            if x is None : return 0
            return len(x)

        if safeLen(options.bamList) == 0 :
            raise OptParseException("No input sample alignment files specified")

        bcheck = BamSetChecker()
        bcheck.appendBams(options.bamList,"Input")
        bcheck.check(options.htsfileBin,
                     options.referenceFasta)



def main() :

    primarySectionName="StrelkaGermline"
    options,iniSections=StrelkaGermlineWorkflowOptions().getRunOptions(primarySectionName, version=workflowVersion)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    StrelkaGermlineWorkflow(options,iniSections)

    # generate runscript:
    #
    ensureDir(options.runDir)
    scriptFile=os.path.join(options.runDir,"runWorkflow.py")

    makeRunScript(scriptFile,os.path.join(workflowDir,"strelkaGermlineWorkflow.py"),"StrelkaGermlineWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (scriptFile))


if __name__ == "__main__" :
    main()

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
This script configures the strelka germline small variant calling workflow
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

sys.path.append(workflowDir)

from configBuildTimeInfo import workflowVersion
from strelkaSharedOptions import StrelkaSharedWorkflowOptionsBase
from configureUtil import BamSetChecker, groomBamList, joinFile, OptParseException, checkFixTabixIndexedFileOption
from makeRunScript import makeRunScript
from strelkaGermlineWorkflow import StrelkaGermlineWorkflow
from workflowUtil import ensureDir, exeFile



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
                         help="Provide BED file of regions where gVCF block compression is not allowed. File must be bgzip-compressed/tabix-indexed. (no default)")
        group.add_option("--callContinuousVf", type="string", dest="callContinuousVf", metavar="CHROM", action="append",
                         help="Call variants on CHROM without a ploidy prior assumption, issuing calls with continuous variant frequencies (no default)")
        group.add_option("--rna", dest="isRNA", action="store_true",
                         help="Set options for RNA-Seq input.")

        StrelkaSharedWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def addExtendedGroupOptions(self,group) :
        # note undocumented library behavior: "dest" is optional, but not including it here will
        # cause the hidden option to always print
        group.add_option("--disableSequenceErrorEstimation", dest="isEstimateSequenceError", action="store_false",
                         help="Disable estimation of sequence error rates from data.")
        group.add_option("--useAllDataForSequenceErrorEstimation", dest="isErrorEstimationFromAllData", action="store_true",
                         help="Instead of sampling a subset of data for error estimation, use all data from sufficiently large chromosomes."
                              " This could greatly increase the workflow's runtime.")

        StrelkaSharedWorkflowOptionsBase.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StrelkaSharedWorkflowOptionsBase.getOptionDefaults(self)

        libexecDir=defaults["libexecDir"]

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'StrelkaGermlineWorkflow',
            'strelkaGermlineBin' : joinFile(libexecDir,exeFile("starling2")),
            'bgzip9Bin' : joinFile(libexecDir, exeFile("bgzip9")),
            'configDir' : configDir,
            'germlineSnvScoringModelFile' : joinFile(configDir,'germlineSNVScoringModels.json'),
            'germlineIndelScoringModelFile' : joinFile(configDir,'germlineIndelScoringModels.json'),
            'rnaSnvScoringModelFile' : joinFile(configDir,'RNASNVScoringModels.json'),
            'rnaIndelScoringModelFile' : joinFile(configDir,'RNAIndelScoringModels.json'),
            'callContinuousVf' : [],
            'getCountsBin' : joinFile(libexecDir,exeFile("GetSequenceAlleleCounts")),
            'mergeCountsBin' : joinFile(libexecDir,exeFile("MergeSequenceAlleleCounts")),
            'estimateVariantErrorRatesBin' : joinFile(libexecDir,exeFile("EstimateVariantErrorRates")),
            'thetaParamFile' : joinFile(configDir,'theta.json'),
            'indelErrorRateDefault' : joinFile(configDir,'indelErrorModel.json'),
            'isEstimateSequenceError' : True,
            'isErrorEstimationFromAllData' : False
            })
        return defaults


    def validateAndSanitizeOptions(self,options) :

        StrelkaSharedWorkflowOptionsBase.validateAndSanitizeOptions(self,options)

        options.ploidyFilename = checkFixTabixIndexedFileOption(options.ploidyFilename,"ploidy file")
        options.noCompressBed = checkFixTabixIndexedFileOption(options.noCompressBed,"no-compress bed")
        if options.snvScoringModelFile is None :
            if options.isRNA :
                options.snvScoringModelFile = options.rnaSnvScoringModelFile
            else :
                options.snvScoringModelFile = options.germlineSnvScoringModelFile

        if options.indelScoringModelFile is None :
            if options.isRNA :
                options.indelScoringModelFile = options.rnaIndelScoringModelFile
            else :
                options.indelScoringModelFile = options.germlineIndelScoringModelFile

        # Disable dynamic error estimation for Exome
        if options.isExome :
            options.isEstimateSequenceError = False

        # Disable dynamic error estimation for RNA
        if options.isRNA :
            options.isEstimateSequenceError = False

        groomBamList(options.bamList,"input")

        def safeLen(x) :
            if x is None : return 0
            return len(x)

        if safeLen(options.bamList) == 0 :
            raise OptParseException("No input sample alignment files specified")

        bamSetChecker = BamSetChecker()
        bamSetChecker.appendBams(options.bamList,"Input")
        bamSetChecker.check(options.htsfileBin, options.referenceFasta)



def main() :

    primarySectionName="StrelkaGermline"
    options,iniSections=StrelkaGermlineWorkflowOptions().getRunOptions(primarySectionName, version=workflowVersion)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    StrelkaGermlineWorkflow(options)

    # generate runscript:
    #
    ensureDir(options.runDir)
    workflowScriptPath = os.path.join(options.runDir, options.workflowScriptName)

    makeRunScript(workflowScriptPath,os.path.join(workflowDir,"strelkaGermlineWorkflow.py"),"StrelkaGermlineWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (workflowScriptPath))


if __name__ == "__main__" :
    main()

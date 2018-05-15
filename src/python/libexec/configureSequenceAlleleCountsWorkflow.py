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
This script configures the sequence error counts workflow
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
from sequenceAlleleCountsWorkflow import SequenceAlleleCountsWorkflow
from workflowUtil import ensureDir, exeFile



class SequenceAlleleCountsWorkflowOptions(StrelkaSharedWorkflowOptionsBase) :

    def workflowDescription(self) :
        return """Version: %s

This script configures the Strelka sequence allele counts workflow.
""" % (workflowVersion)


    def addWorkflowGroupOptions(self,group) :
        group.add_option("--bam", type="string",dest="bamList",metavar="FILE", action="append",
                         help="Sample BAM or CRAM file. [required] (no default)")
        group.add_option("--excludedRegions", type="string", metavar="FILE", action="append",
                         help="Provide BED file of regions to be excluded from allele count analysis. BED file must be tabix indexed."
                         "argument may be specified multiple times to provide multiple exclusion regions (no default)")
        group.add_option("--knownVariants", type="string", metavar="FILE",
                         help="Provide VCF file of indels and SNVs with known genotype assignments. VCF file must be tabix indexed."
                         "Matching alt alleles in the allele counts file will be marked with a known copy number. There is no"
                         "handling of hom ref assertions, whether remaining unlabeled loci are treated as known hom ref is left to the"
                         "downstream estimation model. Note this option does not promote the known variants to candidate or forced GT"
                         "status, to do so the same VCF file can be resubmitted to the appropriate additional argument. "
                         " Input VCF must be tabixed and normalized.")
        group.add_option("--reportObservedIndels", dest="isReportObservedIndels", action="store_true", default = False,
                         help="Report all observed indels by location in a separate BED file in addition to the"
                         "summary counts")

        StrelkaSharedWorkflowOptionsBase.addWorkflowGroupOptions(self,group)


    def getOptionDefaults(self) :

        self.configScriptDir=scriptDir
        defaults=StrelkaSharedWorkflowOptionsBase.getOptionDefaults(self)

        libexecDir=defaults["libexecDir"]

        configDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_CONFIGDIR@"))
        assert os.path.isdir(configDir)

        defaults.update({
            'runDir' : 'SequenceAlleleCountsWorkflow',
            'workflowScriptName' : 'runWorkflow.py',
            'getCountsBin' : joinFile(libexecDir,exeFile("GetSequenceAlleleCounts")),
            'mergeCountsBin' : joinFile(libexecDir,exeFile("MergeSequenceAlleleCounts")),
            'extraCountsArguments' : None
            })
        return defaults


    def validateAndSanitizeOptions(self,options) :

        StrelkaSharedWorkflowOptionsBase.validateAndSanitizeOptions(self,options)

        if options.excludedRegions is not None :
            for excludeIndex in range(len(options.excludedRegions)) :
                options.excludedRegions[excludeIndex] = \
                    checkFixTabixIndexedFileOption(options.excludedRegions[excludeIndex],"excluded-regions bed")

        if options.knownVariants is not None :
            options.knownVariants = \
                checkFixTabixIndexedFileOption(options.knownVariants,"known-variants vcf")

        groomBamList(options.bamList,"input")

        bamSetChecker = BamSetChecker()

        def singleAppender(bamList,label):
            if len(bamList) > 1 :
                raise OptParseException("More than one %s sample BAM/CRAM files specified" % (label))
            bamSetChecker.appendBams(bamList,label)

        singleAppender(options.bamList,"Input")
        bamSetChecker.check(options.htsfileBin, options.referenceFasta)



def main() :

    primarySectionName="counts"
    options,iniSections=SequenceAlleleCountsWorkflowOptions().getRunOptions(primarySectionName, version=workflowVersion)

    # we don't need to instantiate the workflow object during configuration,
    # but this is done here to trigger additional parameter validation:
    #
    SequenceAlleleCountsWorkflow(options)

    # generate runscript:
    #
    ensureDir(options.runDir)
    workflowScriptPath = os.path.join(options.runDir, options.workflowScriptName)

    makeRunScript(workflowScriptPath,os.path.join(workflowDir,"sequenceAlleleCountsWorkflow.py"),"SequenceAlleleCountsWorkflow",primarySectionName,iniSections)

    notefp=sys.stdout
    notefp.write("""
Successfully created workflow run script.
To execute the workflow, run the following script and set appropriate options:

%s
""" % (workflowScriptPath))


if __name__ == "__main__" :
    main()

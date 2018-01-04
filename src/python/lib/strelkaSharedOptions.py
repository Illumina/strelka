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
Workflow configuration options shared by multiple
configuration scripts.
"""

import os,sys

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)

sys.path.append(scriptDir)

from checkChromSet import getFastaInfo, getTabixChromSet
from configureOptions import ConfigureWorkflowOptions
from configureUtil import assertOptionExists, joinFile, OptParseException, \
                          validateFixExistingDirArg, validateFixExistingFileArg, \
                          checkFixTabixListOption, checkFixTabixIndexedFileOption
from workflowUtil import exeFile, getFastaChromOrderSize, parseGenomeRegion


def cleanLocals(locals_dict) :
    """
    When passed a locals() dictionary, clean out all of the hidden keys and return
    """

    return dict((k,v) for (k,v) in locals_dict.items() if not k.startswith("__") and k != "self")



class StrelkaSharedWorkflowOptionsBase(ConfigureWorkflowOptions) :

    def addWorkflowGroupOptions(self,group) :
        group.add_option("--referenceFasta",type="string",metavar="FILE",
                         help="samtools-indexed reference fasta file [required]")
        group.add_option("--indelCandidates", type="string", dest="indelCandidatesList", metavar="FILE", action="append",
                         help="Specify a VCF of candidate indel alleles. These alleles are always"
                              " evaluated but only reported in the output when they are inferred to exist in the sample."
                              " The VCF must be tabix indexed."
                              " All indel alleles must be left-shifted/normalized, any unnormalized alleles will be ignored."
                              " This option may be specified more than once, multiple input VCFs will be merged."
                              " (default: None)")
        group.add_option("--forcedGT", type="string", dest="forcedGTList", metavar="FILE", action="append",
                         help="Specify a VCF of candidate alleles. These alleles are always"
                              " evaluated and reported even if they are unlikely to exist in the sample."
                              " The VCF must be tabix indexed."
                              " All indel alleles must be left-shifted/normalized, any unnormalized allele will trigger"
                              " a runtime error."
                              " This option may be specified more than once, multiple input VCFs will be merged."
                              " Note that for any SNVs provided in the VCF, the SNV site will be reported (and for gVCF,"
                              " excluded from block compression), but the specific SNV alleles are ignored."
                              " (default: None)")
        group.add_option("--exome", "--targeted", dest="isExome", action="store_true",
                         help="Set options for exome or other targeted input: note in particular that this flag turns off high-depth filters")
        group.add_option("--callRegions", dest="callRegionsBed", metavar="FILE",
                         help="Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. "
                              "No VCF output will be provided outside of these regions. The full genome will still be used "
                              "to estimate statistics from the input (such as expected depth per chromosome). "
                              "Only one BED file may be specified. (default: call the entire genome)")
        group.add_option("--runDir", type="string", metavar="DIR",
                         help="Name of directory to be created where all workflow scripts and output will be written. "
                              "Each analysis requires a separate directory. (default: %default)")


    def addExtendedGroupOptions(self,group) :
        # note undocumented library behavior: "dest" is optional, but not including it here will
        # cause the hidden option to always print
        group.add_option("--scanSizeMb", type="int", dest="scanSizeMb", metavar="INT",
                         help="Maximum sequence region size (in megabases) scanned by each task during "
                              "genome variant calling. (default: %default)")
        group.add_option("--region", type="string",dest="regionStrList",metavar="REGION", action="append",
                         help="Limit the analysis to one or more genome region(s) for debugging purposes. "
                              "If this argument is provided multiple times the union of all specified regions will "
                              "be analyzed. All regions must be non-overlapping to get a meaningful result. "
                              "Examples: '--region chr20' (whole chromosome), "
                              "'--region chr2:100-2000 --region chr3:2500-3000' (two regions)'. If this "
                              "option is specified (one or more times) together with the --callRegions BED file, then "
                              "all region arguments will be intersected with the callRegions BED track.")
        group.add_option("--callMemMb",dest="callMemMbOverride",type="int",metavar="INT",
                         help="Set variant calling task memory limit (in megabytes). It is not "
                              "recommended to change the default in most cases, but this might be required "
                              "for a sample of unusual depth.")
        group.add_option("--retainTempFiles", dest="isRetainTempFiles", action="store_true",
                         help="Keep all temporary files (for workflow debugging)")
        group.add_option("--disableEVS", dest="isEVS", action="store_false",
                         help="Disable empirical variant scoring (EVS).")
        group.add_option("--reportEVSFeatures", dest="isReportEVSFeatures", action="store_true",
                         help="Report all empirical variant scoring features in VCF output.")
        group.add_option("--snvScoringModelFile", type="string", dest="snvScoringModelFile", metavar="FILE",
                         help="Provide a custom empirical scoring model file for SNVs (default: %default)")
        group.add_option("--indelScoringModelFile", type="string", dest="indelScoringModelFile", metavar="FILE",
                         help="Provide a custom empirical scoring model file for indels (default: %default)")

        ConfigureWorkflowOptions.addExtendedGroupOptions(self,group)


    def getOptionDefaults(self) :
        """
        Set option defaults.

        Every local variable in this method becomes part of the default hash
        """

        configCommandLine=sys.argv

        libexecDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_LIBEXECDIR@"))
        assert os.path.isdir(libexecDir)

        bgzipBin=joinFile(libexecDir,exeFile("bgzip"))
        htsfileBin=joinFile(libexecDir,exeFile("htsfile"))
        samtoolsBin=joinFile(libexecDir,exeFile("samtools"))
        tabixBin=joinFile(libexecDir,exeFile("tabix"))
        bgcatBin=joinFile(libexecDir,exeFile("bgzf_cat"))

        getChromDepthBin=joinFile(libexecDir,exeFile("GetChromDepth"))

        mergeChromDepth=joinFile(libexecDir,"mergeChromDepth.py")
        catScript=joinFile(libexecDir,"cat.py")
        vcfCmdlineSwapper=joinFile(libexecDir,"vcfCmdlineSwapper.py")

        statsMergeBin=joinFile(libexecDir,exeFile("MergeRunStats"))

        workflowScriptName = "runWorkflow.py"

        # default memory request per process-type
        #
        # where different values are provided for SGE and local runs note:
        #  1. for SGE the memory limits must be greater than the highest memory use ever
        #      expected in a production run. The consequence of exceeding this limit is a failed
        #      run.
        #   2. for localhost the memory usage should be at least above the highest mean memory
        #       use ever expected in a production run. The consequence of exceeding the mean is
        #       a slow run due to swapping.
        #
        callSGEMemMb=4*1024
        callLocalMemMb=1.5*1024


        runDir = "variantCallWorkflow"

        # extended options
        maxIndelSize = 49
        scanSizeMb = 12
        regionStrList = None
        callMemMbOverride = None

        isExome = False

        # isRNA is shared by all Strelka workflows, but only can be set for the germline calling case:
        isRNA = False

        isRetainTempFiles = False

        # Empirical Variant Scoring:
        isEVS = True
        isReportEVSFeatures = False

        indelErrorModelName = None
        inputIndelErrorModelsFile = None

        snvScoringModelFile = None
        indelScoringModelFile = None

        # error estimation is planned for all workflows, but can only be set true in germline at present:
        isEstimateSequenceError = False

        errorEstimationMinChromMb = 5
        errorEstimationMinTotalMb = 50

        return cleanLocals(locals())


    def validateAndSanitizeOptions(self,options) :

        assertOptionExists(options.runDir,"run directory")
        options.runDir=os.path.abspath(options.runDir)

        workflowScriptPath = os.path.join(options.runDir, options.workflowScriptName)
        if os.path.exists(workflowScriptPath):
            raise OptParseException("Run directory already contains workflow script file '%s'. Each analysis must be configured in a separate directory." % (workflowScriptPath))

        assertOptionExists(options.referenceFasta,"reference fasta file")
        options.referenceFasta=validateFixExistingFileArg(options.referenceFasta,"reference fasta file")

        # check for reference fasta index file:
        referenceFastaIndex=options.referenceFasta + ".fai"
        if not os.path.isfile(referenceFastaIndex) :
            raise OptParseException("Can't find expected fasta index file: '%s'" % (referenceFastaIndex))

        if options.isEstimateSequenceError :
            # Determine if dynamic error estimation is feasible based on the reference size
            # - Given reference contig set (S) with sequence length of at least 5 Mb
            # - The total sequence length from S must be at least 50 Mb

            class Constants :
                Megabase = 1000000
                minChromSize = options.errorEstimationMinChromMb*Megabase
                minTotalSize = options.errorEstimationMinTotalMb*Megabase

            # read fasta index
            (_, chromSizes) = getFastaChromOrderSize(referenceFastaIndex)

            totalEstimationSize=0
            for chromSize in chromSizes.values() :
                if chromSize < Constants.minChromSize : continue
                totalEstimationSize += chromSize

            if totalEstimationSize < Constants.minTotalSize :
                sys.stderr.write("WARNING: Cannot estimate sequence errors from data due to small or overly fragmented reference sequence. Sequence error estimation disabled.\n")
                options.isEstimateSequenceError = False

        checkFixTabixListOption(options.indelCandidatesList,"candidate indel vcf")
        checkFixTabixListOption(options.forcedGTList,"forced genotype vcf")
        options.callRegionsBed = checkFixTabixIndexedFileOption(options.callRegionsBed,"call-regions bed")

        def extendedRegionStrList() :
            """
            A generator on the regionStrList which parses the (intentionally undocumented/possibly deprecated) '+' entry format
            to specify multiple regions in a single argument.
            """
            for r in options.regionStrList :
                for rr in r.split("+") :
                    yield rr

        if (options.regionStrList is None) or (len(options.regionStrList) == 0) :
            options.genomeRegionList = None
        else :
            options.genomeRegionList = [parseGenomeRegion(r) for r in extendedRegionStrList()]

        # validate chromosome names appearing in region tags and callRegions bed file
        if (options.callRegionsBed is not None) or (options.genomeRegionList is not None) :
            refChromInfo = getFastaInfo(options.referenceFasta)
            if options.callRegionsBed is not None :
                for chrom in getTabixChromSet(options.tabixBin, options.callRegionsBed) :
                    if chrom not in refChromInfo :
                        raise OptParseException("Chromosome label '%s', in call regions bed file '%s', not found in reference genome." %
                                                (chrom, options.callRegionsBed))

            if options.genomeRegionList is not None :
                for (genomeRegionIndex, genomeRegion) in enumerate(options.genomeRegionList) :
                    chrom = genomeRegion["chrom"]
                    if chrom not in refChromInfo :
                        raise OptParseException("Chromosome label '%s', parsed from region argument '%s', not found in reference genome." %
                                                (chrom, list(extendedRegionStrList())[genomeRegionIndex]))

        options.snvScoringModelFile=validateFixExistingFileArg(options.snvScoringModelFile,"SNV empirical scoring model file")
        options.indelScoringModelFile=validateFixExistingFileArg(options.indelScoringModelFile,"Indel empirical scoring model file")

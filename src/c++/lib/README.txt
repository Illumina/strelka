applications/X:
code specific to command-line application X

blt_common:
shared snv-calling logic (from old 'blt' snp-caller)

blt_util:
general utility functions from manta/starling/strelka/gvcftools

common:
general utility functions from CASAVA/Grouper/Isaac

errorAnalysis:
shared libraries for utilities that estimate sequence artifact
rates from data, used for various forms of model parameter estimation

htsapi:
various c++ wrapper objects built on top of samtools/htslib
and other utilities for standard genomic indexed file formats
like bam/cram,bed,vcf, etc...

options:
command-line options objects which are shared between applications

starling_common:
shared small variant calling logic

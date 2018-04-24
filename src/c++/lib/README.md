# Libraries

- alignment
  - sequence alignment utilities

- applications/X:
  - code specific to command-line application X

- appstats
  - shared performance tracking code between applications

- assembly
  - local sequence assembly utilities

- blt\_common
  - shared snv-calling logic (from old 'blt' snp-caller)

- blt\_util
  - general utility functions from manta/starling/strelka/gvcftools

- calibration
  - empirical variant scoring logic

- common
  - general utility functions from CASAVA/Grouper/Isaac

- errorAnalysis
  - shared libraries for utilities that estimate sequence artifact rates from data, used for various forms of model parameter estimation

- htsapi
  - various c++ wrapper objects built on top of samtools/htslib and other utilities for standard genomic indexed file formats like bam/cram,bed,vcf, etc...

- options
  - command-line options objects which are shared between applications

- starling\_common
  - shared small variant calling logic

- strelka\_common
  - shared comparitive small variant calling logic

- test
  - Logic used only for unit testing other libraries. These are linked into the unit tests but not production binaries

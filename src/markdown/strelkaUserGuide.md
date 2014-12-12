<link rel='stylesheet' href='userGuide.css' />

Strelka User Guide
==================

Version: @STARKA_VERSION@

<script src="tableOfContents.js"></script>

## Introduction

Strelka calls somatic SNVs and small indels from short sequencing reads corresponding to
a tumor and matched normal sample. It is designed to handle impurity in the tumor sample.

## Method Overview

The strelka calling algorithm fully is described in
[Strelka: Accurate somatic small-variant calling from sequenced tumor-normal sample pairs.][3]

In summary strelka scans throught the tumor and normal sample alignments, discovering SNV
and indel candidates. Each candidate variant is evaluated as potentially germline, somatic or
sequencer artifact, with a quality score reflecting the final probability of being somatic.

## Capabilities

xxxx

### Known Limitations

Strelka requires a matched normal sample to make somatic calls. The matched
normal is used to distinguish both germline variation and sequencing artifact from
somatic variation. The general depth guideline for the normal sample is either
one half the tumor depth or ~30x, whichever is higher.

## Input requirements

The sequencing reads provided as input to Strelka are expected to be from a
paired-end sequencing assay.

Strelka requires input sequencing reads to be mapped by an external tool and
provided as input in BAM or CRAM format.

The following limitations exist on the input BAM/CRAM alignment records provided
to Starka:

* Alignments cannot contain the "=" character in the SEQ field.
* Alignments cannot use the sequence match/mismatch ("="/"X") CIGAR notation
* RG (read group) tags in the BAMs are ignored -- each BAM must represent one
  sample.
* Alignments with basecall quality values greater than 70 are rejected (these
  are not supported on the assumption that this indicates an offset error)

## Outputs

### Somatic variant predictions

The primary strelka outputs are a set of [VCF 4.1][1] files, found in
`${RUNFOLDER}/results/variants`. Currently there are at least two vcf files
created for any run. These files are:

* __somatic.snvs.vcf.gz__
    * snvs
* __somatic.indels.vcf.gz__
    * indels
    
Strelka can also optionally produce a somatic callability track in:

'${RUNFOLDER}/results/regions/somatic.callable.region.bed.gz'

## Run configuration and Execution

Strelka is run in a two step procedure: (1) configuration and (2) workflow
execution. The configuration step is used to specify the input data and any
options pertaining to the variant calling methods themselves. The execution
step is used to specify any parameters pertaining to _how_ strelka is executed
(such as the total number of cores or SGE nodes over which the jobs should be
parallelized). The second execution step can also be interrupted and restarted
without changing the final result of the workflow.

### Configuration

The workflow is configured with the script: `${INSTALL_DIR}/bin/configureStrelkaWorkflow.py`
. Running this script with no arguments will display all standard configuration
options to specify input BAM files, the reference sequence and the output run folder.
Note that all input BAMs and reference sequence must contain the same chromosome names
in the same order. Strelka's default settings assume a whole genome DNA-Seq analysis.

Example Configuration:

    ${INSTALL_DIR}/bin/configureStrelkaWorkflow.py \
    --config ${INSTALL_DIR}/share/config/strelka_config_isaac_default.ini \
    --normalBam HCC1187BL.bam \
    --tumorBam HCC1187C.bam \
    --referenceFasta hg19.fa \
    --runDir ${ANALYSIS_RUN_DIR}

On completion, the configuration script will create the workflow run script `${ANALYSIS_RUN_DIR}/runWorkflow.py`
. This can be used to run the workflow in various parallel compute modes per the
instructions in the [Execution] section below.

#### Advanced configuration options

* Advanced options listed in: `${INSTALL_DIR}/bin/configureStrelkaWorkflow.py -- allHelp`
    * These options are indented primarily for workflow development and
      debugging, but could be useful for runtime optimization in some specialized
      cases.

### Execution

The configuration step creates a new workflow run script in the requested run directory:

`{ANALYSIS_RUN_DIR}/runWorkflow.py`

This script is used to control parallel execution of the workflow via the [pyFlow][2]
task engine. It can be used to parallelize structural variant analysis via one
of two modes:

1. Parallelized across multiple cores on a single node.
2. Parallelized across multiple nodes on an SGE cluster.

A running workflow can be interrupted at any time and resumed where it left
off. If desired, the resumed analysis can use a different running mode or total
core count.

For a full list of execution options, see:

`{ANALYSIS_RUN_DIR}/runWorkflow.py -h`

Example execution on a single node:

`${ANALYSIS_RUN_DIR}/runWorkflow.py -m local -j 8`

Example execution on an SGE cluster:

`${ANALYSIS_RUN_DIR}/runWorkflow.py -m sge -j 36`

#### Advanced execution options

These options are useful for workflow development and debugging:

* Stderr logging can be disabled with `--quiet` argument. Note this log is
  replicated to `${ANALYSIS_RUN_DIR}/workspace/pyflow.data/logs/pyflow_log.txt`
  so there is no loss of log information.

[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[2]: http://ctsa.github.io/pyflow/
[3]: http://bioinformatics.oxfordjournals.org/content/28/14/1811

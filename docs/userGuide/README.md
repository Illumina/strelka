Strelka User Guide
==================

## Table of Contents
[] (BEGIN automated TOC section, any edits will be overwritten on next source refresh)
* [Introduction](#introduction)
* [Installation](#installation)
* [Method Overview](#method-overview)
* [Capabilities](#capabilities)
  * [Known Limitations](#known-limitations)
* [Input requirements](#input-requirements)
* [Outputs](#outputs)
  * [Somatic variant predictions](#somatic-variant-predictions)
* [Run configuration and Execution](#run-configuration-and-execution)
  * [Configuration](#configuration)
    * [Advanced configuration options](#advanced-configuration-options)
  * [Execution](#execution)
    * [Advanced execution options](#advanced-execution-options)
* [Special Topics](#special-topics)
[] (END automated TOC section, any edits will be overwritten on next source refresh)

## Introduction

Strelka calls somatic SNVs and small indels from short sequencing reads corresponding to
a tumor and matched normal sample. It is designed to handle impurity in the tumor sample.

## Installation

Please see the [Strelka installation instructions](installation.md)

## Method Overview

The strelka somatic variant calling algorithm fully is described in
[Strelka: Accurate somatic small-variant calling from sequenced tumor-normal sample pairs.][3]

In summary strelka scans through the tumor and normal sample alignments, discovering SNV
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

The following limitations apply to the BAM/CRAM alignment records which can be used as input:

* Alignments cannot contain the "=" character in the SEQ field.
* RG (read group) tags are ignored -- each alignment file must represent one
  sample.
* Alignments with basecall quality values greater than 70 are rejected (these
  are not supported on the assumption that this indicates an offset error)

## Outputs

### Somatic variant predictions

The primary strelka outputs are a set of [VCF 4.1][1] files, found in
`${STRELKA_ANALYSIS_PATH}/results/variants`. Currently there are at least two vcf files
created for any run. These files are:

* __somatic.snvs.vcf.gz__
    * snvs
* __somatic.indels.vcf.gz__
    * indels

The somatic variant caller can also optionally produce a callability tract,
see the [somatic callability](#somatic-callability) section below for details.

## Run configuration and Execution

Strelka is run in a two step procedure: (1) configuration and (2) workflow
execution. The configuration step is used to specify the input data and any
options pertaining to the variant calling methods themselves. The execution
step is used to specify any parameters pertaining to _how_ strelka is executed
(such as the total number of cores or SGE nodes over which the jobs should be
parallelized). The second execution step can also be interrupted and restarted
without changing the final result of the workflow.

### Configuration

The workflow is configured with the script: `${STRELKA_INSTALL_PATH}/bin/configureStrelkaWorkflow.py`
. Running this script with no arguments will display all standard configuration
options to specify input alignment files, the reference sequence and the output run folder.
Note that all input alignment and reference sequence files must contain the same chromosome names
in the same order. Strelka's default settings assume a whole genome DNA-Seq analysis.

Example Configuration:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaWorkflow.py \
    --config ${STRELKA_INSTALL_PATH}/share/config/strelka_config_isaac_default.ini \
    --normalBam HCC1187BL.bam \
    --tumorBam HCC1187C.bam \
    --referenceFasta hg19.fa \
    --runDir ${STRELKA_ANALYSIS_PATH}

On completion, the configuration script will create the workflow run script `${STRELKA_ANALYSIS_PATH}/runWorkflow.py`
. This can be used to run the workflow in various parallel compute modes per the
instructions in the [Execution] section below.

#### Advanced configuration options

* Advanced options are listed in: `${STRELKA_INSTALL_PATH}/bin/configureStrelkaWorkflow.py -- allHelp`
    * These options are indented primarily for workflow development and
      debugging, but could be useful for runtime optimization in some specialized
      cases.

### Execution

The configuration step creates a new workflow run script in the requested run directory:

`${STRELKA_ANALYSIS_PATH}/runWorkflow.py`

This script is used to control parallel execution of the workflow via the [pyFlow][2]
task engine. It can be used to parallelize structural variant analysis via one
of two modes:

1. Parallelized across multiple cores on a single node.
2. Parallelized across multiple nodes on an SGE cluster.

A running workflow can be interrupted at any time and resumed where it left
off. If desired, the resumed analysis can use a different running mode or total
core count.

For a full list of execution options, see:

`${STRELKA_ANALYSIS_PATH}/runWorkflow.py -h`

Example execution on a single node:

`${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m local -j 8`

Example execution on an SGE cluster:

`${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m sge -j 36`

#### Advanced execution options

These options are useful for workflow development and debugging:

* Stderr logging can be disabled with `--quiet` argument. Note this log is
  replicated to `${STRELKA_ANALYSIS_PATH}/workspace/pyflow.data/logs/pyflow_log.txt`
  so there is no loss of log information.

### Extended use cases

#### Exome/Targeted

Supplying the '--exome' flag at configuration time will provide
appropriate settings for WES and other regional enrichment
analyses. At present this flag disables all high depth filters, which
are designed to exclude pericentromeric reference compressions in the
WGS case but cannot be applied correctly to a targeted analysis.

#### Somatic callability

The somatic variant caller can be configured with the option `--writeCallableRegion`, which
will extend the somatic SNV quality model calculation to be applied as a test of
somatic SNV callability at all positions in the genome.

The outcome of this callibility calculation will be summarized in a BED-formatted callability track
found in:`${STRELKA_ANALYSIS_PATH}/results/regions/somatic.callable.region.bed.gz`. This BED track
contains regions which are determined to be callable, indicating that there is sufficient evidence to
either call a somatic SNV or assert the absence of a somatic SNV with a variant frequency of 10% or greater.
Both somatic and non-somatic sites are determined to be 'callable' if the somatic or non-somatic quality
threshold is at least 15. See methods for details of the underlying quality scores.

This is still an experimental feature, which will considerably increase runtime cost of the analysis
(by approximately 2x).


## Special Topics

The following items provide an in-depth focus on a special topic or procedure

* [Training Procedure for Somatic Empirical Score](trainingSomaticEmpiricalScore.md)



[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[2]: http://Illumina.github.io/pyflow/
[3]: http://bioinformatics.oxfordjournals.org/content/28/14/1811

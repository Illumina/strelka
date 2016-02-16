Strelka User Guide - Germline gVCF workflow
===========================================

[User Guide Home](README.md)

## Introduction

Starling is a small variant caller for short sequencing reads. It is capable of calling snps and small indels. Starling will call simple and complex indels up to a specified maximum size. The output of all variant calls is combined with reference block calls to produce a single gVCF output file for the genome.

## Method Overview

For any genomic regions, starling works in multiple phases. Its calling process starts with
candidate indel discovery, followed by read realignment, SNP and indel genotyping, and finally
filtration and scoring. These steps are described in more detail below:

1. **Candidate indel discovery**

2. **Read Realignment** 

3. **SNP and indel genotyping**

4. **Filtration and scoring**

## Capabilities

WGS. WES. RNA-Seq. Amplicon

## Input requirements

The sequencing reads provided as input to starling are expected to be from a
paired-end sequencing assay with an "innie" orientation between the two reads of
each DNA fragment, each presenting a read from the outer edge of the fragment
insert inward.

Starling can tolerate non-paired reads in the input, but these will be ignored for
during variant calling by default.

Starling requires input sequencing reads to be mapped by an external tool and
provided as input in BAM or CRAM format.

At configuration time, a bam file must be provided for the query sample.

The following limitations exist on the input BAMCRAM files provided to Starling:

* Alignments cannot contain the "=" character in the SEQ field.
* RG (read group) tags are ignored -- each alignment file must represent one
  sample.
* Alignments with basecall quality values greater than 70 are rejected (these
  are not supported on the assumption that this indicates an offset error)

## Outputs

### Genome VCF (gVCF)

The primary starling output is a [VCF 4.1][1] file found in
`${STRELKA_ANALYSIS_PATH}/results/variants`:

* __genome.vcf.gz__
    * This file represents the genotype for all positions in the genome, in addition to all discovered indel variants


### Statistics

Additional secondary output is provided in ${RUNFOLDER}/results/stats

* __genomeCallStats.tsv__
    * runtime information accumulated for each genome segment. Excludes auxiliary steps such as BAM indexing and vcf merging.

* __genomeCallStats.xml__
    * xml data backing the genomeCallStats.tsv report


## Run configuration and Execution

Starling is run in a two step procedure: (1) configuration and (2) workflow
execution. The configuration step is used to specify the input data and any
options pertaining to the variant calling methods themselves. The execution
step is used to specify any parameters pertaining to _how_ starling is executed
(such as the total number of cores or SGE nodes over which the jobs should be
parallelized). The second execution step can also be interrupted and restarted
without changing the final result of the workflow.

Note in the guidelines below that starling is made available as part of the STARKA
package, which also includes the strelka somatic small variant caller. For this reason
the root installation referenced below is `${STRELKA_INSTALL_PATH}`.

### Configuration

The workflow is configured with the script: `${STRELKA_INSTALL_PATH}/bin/configureStarlingWorkflow.py`
. Running this script with no arguments will display all standard configuration
options to specify input alignment files, the reference sequence and the output run folder.
Note that all input alignment and reference sequence files must contain the same chromosome names
in the same order. Starka's default settings assume a whole genome DNA-Seq analysis.

Simple WGS Analysis -- Example Configuration:

    ${STRELKA_INSTALL_PATH}/bin/configureStarlingWorkflow.py \
    --bam NA12878_S1.bam \
    --referenceFasta hg19.fa \
    --runDir ${STRELKA_ANALYSIS_PATH}

On completion, the configuration script will create the workflow run script `${STRELKA_ANALYSIS_PATH}/runWorkflow.py`
. This can be used to run the workflow in various parallel compute modes per the
instructions in the [Execution] section below.

#### Advanced configuration options

There are two sources of advanced configuration options:

* Options listed in the file: `${STRELKA_INSTALL_PATH}/bin/configureStarlingWorkflow.py.ini`
    * These parameters are not expected to change frequently. Changing the file
  listed above will re-configure all starling runs for the installation. To change
  parameters for a single run, copy the configureStarlingWorkflow.py.ini file to another location,
  change the desired parameter values and supply the new file using the configuration
  script's `--config FILE` option.
* Advanced options listed in: `${STRELKA_INSTALL_PATH}/bin/configureStarlingWorkflow.py --allHelp`
    * These options are indented primarily for workflow development and
  debugging, but could be useful for runtime optimization in some specialized
  cases.

### Execution

The configuration step creates a new workflow run script in the requested run directory:

`{STRELKA_ANALYSIS_PATH}/runWorkflow.py`

This script is used to control parallel execution of the workflow via the [pyFlow][2]
task engine. It can be used to parallelize structural variant analysis via one
of two modes:

1. Parallelized across multiple cores on a single node.
2. Parallelized across multiple nodes on an SGE cluster.

A running workflow can be interrupted at any time and resumed where it left
off. If desired, the resumed analysis can use a different running mode or total
core count.

For a full list of execution options, see:

`{STRELKA_ANALYSIS_PATH}/runWorkflow.py -h`

Example execution on a single node:

`${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m local -j 8`

Example execution on an SGE cluster:

`${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m sge -j 36`

### Extended use cases

#### Exome/Targeted

Supplying the '--exome' flag at configuration time will provide
appropriate settings for WES and other regional enrichment
analyses. At present this flag disables all high depth filters, which
are designed to exclude pericentromeric reference compressions in the
WGS case but cannot be applied correctly to a targeted analysis. It also
changes the indel error model and quality scoring to a Qrule, which is
less tuned for WGS.

#### Heteroplasmic/non-ploidy calling

Supplying the '--callContinuousVf chrName' argument at configuration time will
change variant calling to not use any ploidy prior assumption on the named chromosome.
Instead, variants will be called with continuous frequencies and scored using a poisson
model. This is useful for detecting variants in mixed populations. This
argument may be passed multiple times to specify multiple chromosomes to
call in this manner.

#### Advanced execution options

These options are useful for workflow development and debugging:

* Stderr logging can be disabled with `--quiet` argument. Note this log is
  replicated to `${STRELKA_ANALYSIS_PATH}/workspace/pyflow.data/logs/pyflow_log.txt`
  so there is no loss of log information.

[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[2]: http://ctsa.github.io/pyflow/

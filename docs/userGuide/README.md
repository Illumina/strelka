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
  * [Variant prediction](#variant-prediction)
    * [Germline](#germline)
    * [Somatic](#somatic)
  * [Statistics](#statistics)
* [Run configuration and Execution](#run-configuration-and-execution)
  * [Configuration](#configuration)
    * [General configuration options](#general-configuration-options)
    * [Germline configuration options](#germline-configuration-options)
    * [Somatic configuration options](#somatic-configuration-options)
    * [Advanced configuration options](#advanced-configuration-options)
  * [Execution](#execution)
    * [Advanced execution options](#advanced-execution-options)
  * [Extended use cases](#extended-use-cases)
    * [Exome/Targeted](#exometargeted)
    * [RNA-Seq](#rna-seq)
    * [Heteroplasmic/pooled calling](#heteroplasmicpooled-calling)
    * [Somatic callability](#somatic-callability)
* [Special Topics](#special-topics)
[] (END automated TOC section, any edits will be overwritten on next source refresh)

## Introduction

Strelka calls small variants from mapped sequencing reads. It is optimized for rapid clinical analysis of germline variation in small cohorts and somatic
variation in tumor/normal sample pairs.

## Installation

Please see the [Strelka installation instructions](installation.md)

## Method Overview

The strelka workflow comprises of a number of common sequence analysis
steps followed by application-specific variant modeling and empirical rescoring
methods specific to the analysis of germline or somatic variation.

In all cases a preliminary step is run to estimate various genomic or regional
statistics from the input alignments, including the sequence depth distribution
and sequence error patterns. This is followed by a division of genome into segments
for parallel processing, where within each chunk the input samples are jointly analyzed
to identify candidate alleles, realign all input reads, analyze reads to make model
specific variant inferences, then compute properties of each variant used to apply filters
or empirically recalibrate confidence that each variant represents a germline or somatic
variant in the input sample(s). Finally, all parallel segment results are joined to produce
Strelka's final variant output.

## Capabilities

Strelka is capable of detecting SNVs and indels up to a predefined maximum size, currently
defaulting to 50 or less. Indels are detected from several sources, including alignments
inserted into read alignments by the mapper, candidate indel VCFs provided as input from an
SV caller (such as Manta), or under specific conditions the indel may be detected from the
assembly of an active region.

All methods are optimized for WGS of DNA, but are routinely tested for WES and amplicon inputs.
RNA-Seq germline analysis is still in development and not fully supported. It
can be configured with the `--rna` flag. This will adjust filtration
levels, heterozygous allele ratio expectations and empirical scoring

Strelka includes a short-range read-backed phasing capability for germline calls to facilitate
the correct annotation of codon changes induced by proximal SNVs.

While Strelka is capable of performing joint germline analysis on a family scale (10s of samples),
it is not optimized for population analysis and may become unstable or fail to leverage
population variant constraints to improve calls at higher sample counts.

Strelka's somatic calling capability is known to provide good results down to about 5-10% tumor purity given
sufficient normal and tumor sequencing depth. Strelka also accounts for minor contamination of the normal
sample with tumor cells (up to 10%) to better support liquid and late-stage solid tumor analysis.


### Known Limitations

Strelka requires a matched normal sample to make somatic calls. The matched
normal is used to distinguish both germline variation and sequencing artifact from
somatic variation. The general depth guideline for the normal sample is either
one half the tumor depth or ~30x, whichever is higher.

## Input requirements

The sequencing reads provided as input to Strelka are expected to be from a
paired-end sequencing assay. Any non-paired reads in the input are ignored
by default during variant calling.

Strelka requires input sequencing reads to be mapped by an external tool and
provided as input in BAM or CRAM format.

The following limitations apply to the BAM/CRAM alignment records which can be used as input:

* Alignments cannot contain the "=" character in the SEQ field.
* RG (read group) tags are ignored -- each alignment file must represent one
  sample.
* Alignments with basecall quality values greater than 70 are rejected (these
  are not supported on the assumption that this indicates an offset error)

## Outputs

### Variant prediction

Primary variant inferences are provided as a series of [VCF 4.1][1] files in
`${STRELKA_ANALYSIS_PATH}/results/variants`.

#### Germline

Germline analysis is reported to the following variant files:

* __variants.vcf.gz__
    * This describes all potential variant loci across all samples. Note this file includes non-variant loci if they have a non-trivial level of variant evidence or contain
      one or more alleles for which genotyping has been forced.
* __genome.S${N}.vcf.gz__
    * This is the genome VCF output for sample ${N}, which includes both variant records and compressed non-variant blocks. The sample index, ${N} is 1-indexed and corresponds to the input order of alignment files on the configuration command-line.

##### Germline VCF Sample Names

Sample names printed into the VCF output are extracted from each input
alignment file from the first read group ('@RG') record found in the
header. Any spaces found in the name will be replaced with
underscores. If no sample name is found a default SAMPLE1, SAMPLE2,
etc.. label will be used instead.


#### Somatic

Somatic analysis provides somatic variants in the following two files:

* __somatic.snvs.vcf.gz__
    * All somatic SNVs inferred in the tumor sample.
* __somatic.indels.vcf.gz__
    * All somatic indels inferred in the tumor sample.

The somatic variant caller can also optionally produce a callability track,
see the [somatic callability](#somatic-callability) section below for details.

### Statistics

Additional secondary output is provided in ${STRELKA_ANALYSIS_PATH}/results/stats

* __genomeCallStats.tsv__
    * This file provides runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging.

* __genomeCallStats.xml__
    * xml data backing the genomeCallStats.tsv report

## Run configuration and Execution

Strelka is run in a two step procedure: (1) configuration and (2) workflow
execution. The configuration step is used to specify the input data and any
options pertaining to the variant calling methods themselves. The execution
step is used to specify any parameters pertaining to _how_ strelka is executed
(such as the total number of cores or SGE nodes over which the jobs should be
parallelized). The second execution step can also be interrupted and restarted
without changing the final result of the workflow.

### Configuration

All workflows are configured with a set of workflow scripts following the pattern: `${STRELKA_INSTALL_PATH}/bin/configureStrelka${type}Workflow.py`
. Running any of these scripts with no arguments will display all standard configuration
options to specify input alignment files, the reference sequence and the output run folder.
Note that all input alignment and reference sequence files must contain the same chromosome names
in the same order. In all workflows, Strelka's default settings assume a whole genome DNA-Seq analysis,
but there are configuration options for exome/targeted analysis, in addition to RNA-Seq options for the
germline workflow.

On completion, the configuration script will create the workflow run script `${STRELKA_ANALYSIS_PATH}/runWorkflow.py`. This can be used to run the workflow in various
parallel compute modes per the instructions in the [Execution] section below.

#### General configuration options

##### Candidate indels

One or more candidate indel VCFs can be provided to any Strelka workflow with the `--indelCandidates` flag. Any indels provided in this way will be given candidate status and considered during realignment and genotyping steps, but will not be output unless the indel allele is found in the input samples. This is a useful mechanism to supply Strelka with larger indels (ie. larger than can be found by the read mapper), or to provide population variants.

Any candidate indel record which is not left-normalized will be skipped with a warning.

Multiple candidate indel VCFs may be submitted to the workflow (e.g. `--indelCandidates cand1.vcf.gz --indelCandidates cand2.vcf.gz ...`). All input VCFs must be bgzip compressed and tabix-indexed.

##### Forced genotypes

One or more forced genotype VCFs can be provided to any Strelka workflow with the `--forcedGT` flag.  Any indel allele provided in this way will be treated as a candidate (per the `--indelCandidates` option above), and additionally must appear in the output VCF, even when there is no support for the allele in the input samples. Be aware that in certain cases where the forced allele overlaps another called allele, the forced allele may appear in the output on a VCF record with a different position and/or an additional prefix/suffix added to the REF and ALT fields compared to the allele description in the VCF input. Any SNV listed in the forced genotype VCF will prevent the corresponding site form being compressed into a homozygous reference block and ensure that a VCF site record is output for the given position, but will not provide any special treatment of the alternate base(s) listed in the VCF.

Any forced genotype variant record which is not left-normalized will __trigger a runtime error__.

Multiple forced genotype VCFs may be submitted to the workflow (e.g. `--forcedGT fgt1.vcf.gz --forcedGT fgt2.vcf.gz ...`). All input VCFs must be bgzip compressed and tabix-indexed.

#### Germline configuration options

Germline analysis is configured with the script: `${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py`

Single Diploid Sample Analysis -- Example Configuration:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
    --bam NA12878.bam \
    --referenceFasta hg19.fa \
    --runDir ${STRELKA_ANALYSIS_PATH}

Joint Diploid Sample Analysis -- Example Configuration:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
    --bam NA12878.cram \
    --bam NA12891.cram \
    --bam NA12892.cram \
    --referenceFasta hg19.fa \
    --runDir ${STRELKA_ANALYSIS_PATH}

#### Somatic configuration options

Somatic analysis is configured with the script: `${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py`

Example Configuration:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam HCC1187BL.bam \
    --tumorBam HCC1187C.bam \
    --referenceFasta hg19.fa \
    --runDir ${STRELKA_ANALYSIS_PATH}

#### Advanced configuration options

There are two sources of advanced configuration options:

* Advanced options listed in: `${STRELKA_INSTALL_PATH}/bin/configureStrelka${type}Workflow.py --allHelp`
    * These options are indented primarily for workflow development and
  debugging, but could be useful for runtime optimization in some specialized
  cases.
* Options listed in the file: `${STRELKA_INSTALL_PATH}/bin/configureStrelka${type}Workflow.py.ini`
    * These parameters are not expected to change frequently. Changing the file
  listed above will re-configure all runs of the corresponding workflow from this installation.
  To change parameters for a single run, copy the `ini` file to another location,
  change the desired parameter values and supply the new file using the configuration
  script's `--config FILE` option.

##### Advanced configuration options for germline calling

###### Ploidy
Strelka includes an option to specify regions of a diploid genome which should be treated as
haploid (given ploidy of 1) or are expected to be absent (given ploidy of 0). In any region specified
as absent -- all variants will be called under the default (diploid) model, but filtered with the `PloidyConflict` label.

Ploidy is provided to the call using the `--ploidy` configuration option and supplying the copy number
information for all input samples in VCF format using the `FORMAT/CN` tag to records where the `ALT` tag is set to `<CNV>`. An example ploidy input file is:

    ##fileformat=VCFv4.1
    ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
    ##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
    #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  NA12882 NA12878 NA12877
    chrX    0   .   N   <CNV>   .   PASS    END=10000   CN  1   2   1
    chrX    2781479 .   N   <CNV>   .   PASS    END=155701382   CN  1   2   1
    chrX    156030895   .   N   <CNV>   .   PASS    END=156040895   CN  1   2   1
    chrY    0   .   N   <CNV>   .   PASS    END=57227415    CN  1   0   1

The span over which the copy number from each VCF record is applied is:  `[POS+1, INFO/END]`.

Strelka does not read any fields besides `CHROM`, `POS`, `ALT`, `INFO/END` and `FORMAT/CN`, so a ploidy specific record could be further simplified if desired, e.g:

    chrY    0   .   .   <CNV>   .   .    END=57227415    CN  1   0   1

...would be a valid input record for this option.

Note this feature is primarily intended to delineate the sex chromosome copy number but can be used to call small variants in the context of CNV calls as well.


### Execution

The configuration step creates a new workflow run script in the requested run directory:

`${STRELKA_ANALYSIS_PATH}/runWorkflow.py`

This script is used to control parallel execution of the workflow via the [pyFlow][2]
task engine. It can be used to parallelize small variant analysis via one
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

Supplying the `--exome` flag at configuration time will provide
appropriate settings for WES and other regional enrichment
analyses. At present this flag disables all high depth filters, which
are designed to exclude pericentromeric reference compressions in the
WGS case but cannot be applied correctly to a targeted analysis.

For germline analysis, this mode also disables the empirical variant scoring (EVS)
model, falling back to a set of simple threshold based filters instead. The somatic EVS model
remains in use for exome and targeted data.

#### RNA-Seq

The germline workflow can be configured with the `--rna` flag. This will provide
experimental settings for RNA-Seq variant calling. At present this flag
disables all high depth filters which are designed to exclude
pericentromeric reference compressions in the WGS case but cannot be
applied correctly to RNA-Seq analysis. In addition the expected allele frequency
of heterozygous variants is expanded to account for allele specific expression and
a custom RNA-Seq empirical scoring model is used.

#### Heteroplasmic/pooled calling

The germline workflow can be configured with the `--callContinuousVf ${CHROM}` argument. This will
change variant calling to treat `${CHROM}` as a pooled sample: variants will be called with continuous
frequencies and scored using a simple Poisson noise model.

#### Somatic callability

The somatic variant caller can be configured with the option `--outputCallableRegions`, which
will extend the somatic SNV quality model calculation to be applied as a test of
somatic SNV callability at all positions in the genome.

The outcome of this callability calculation will be summarized in a BED-formatted callability track
found in:`${STRELKA_ANALYSIS_PATH}/results/regions/somatic.callable.regions.bed.gz`. This BED track
contains regions which are determined to be callable, indicating that there is sufficient evidence to
either call a somatic SNV or assert the absence of a somatic SNV with a variant frequency of 10% or greater.
Both somatic and non-somatic sites are determined to be 'callable' if the somatic or non-somatic quality
threshold is at least 15. See methods for details of the underlying quality scores.

This is still an experimental feature, which will considerably increase runtime cost of the analysis
(by approximately 2x).


## Special Topics

The following items provide an in-depth focus on a special topic or procedure

* [Training Procedure for Germline Empirical Score](trainingGermlineEmpiricalScore.md)
* [Training Procedure for Somatic Empirical Score](trainingSomaticEmpiricalScore.md)



[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[2]: http://Illumina.github.io/pyflow/
[3]: http://bioinformatics.oxfordjournals.org/content/28/14/1811

Error Pattern Analyzer User Guide
=================================

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Introduction](#introduction)
* [Installation](#installation)
* [Method overview](#method-overview)
* [Input requirements](#input-requirements)
* [Outputs](#outputs)
  * [Counts files](#counts-files)
  * [Error model output](#error-model-output)
* [Allele counting workflow configuration and execution](#allele-counting-workflow-configuration-and-execution)
  * [Configuration](#configuration)
    * [Configuration: Excluding regions](#configuration-excluding-regions)
    * [Configuration: Annotating known variants](#configuration-annotating-known-variants)
    * [Configuration: Report observed indels](#configuration-report-observed-indels)
  * [Execution](#execution)
  * [Advanced execution options](#advanced-execution-options)
    * [`--quiet`](#--quiet)
* [Viewing allele counting workflow ouput](#viewing-allele-counting-workflow-ouput)
  * [Summary output](#summary-output)
  * [Excluding basecalls/indels](#excluding-basecallsindels)
  * [Extended output (for model development)](#extended-output-for-model-development)
* [Error model parameter estimation](#error-model-parameter-estimation)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


## Introduction

The error pattern analyzer is used to evaluate various models of spurious basecall and indel errors in sequencing data. This is an internal Strelka development tool and is not well supported for external users. This module is part of the Strelka code distribution because it (1) reuses many of Strelka's sequence handling libraries (2) potentially could be used to make error model decisions for the variant caller. Note that any usage of the pattern analyzer that directly impacts the variant calling model will be documented as part of the small variant caller itself.

## Installation

This module builds and installs with the Strelka variant caller by default. No special build steps are required.

## Method overview

The error pattern analyzer is comprised of two primary steps.

(1) An allele counting workflow analyzes BAMs to produce per-locus allele distributions over the genome for various allele types and context segmentations. Allele counts are found for segments of the genome in parallel and merged to produce a single counts file for the sample.

(2) Given a counts file from the first step, various error models can be run on the data to evaluate model fit or to parameterize a model to specific sequencing conditions.

More detailed methods documentation for both steps can be found [here](../methods/errorAnalysis/).

## Input requirements

The error counting process requires aligned sequencing reads in BAM or CRAM forma. The BAM/CRAM file restrictions are identical to those for the Strelka small variant caller. The workflow also optionally excepts region and variant input formated as indexed BED and VCF files respectively.

## Outputs

### Counts files

The primary output of the allele counting workflow is a binary allele counts file, currently written
to `${COUNTS_ANALYSIS_PATH}/results/variants/strelkaAlleleCounts.bin`. This is a binary format. To view the counts data without running a model, see [Viewing allele counting workflow ouput](#viewing-allele-counting-workflow-ouput).

### Error model output

All error models applied to the counts file currently write to stdout, typically in csv format.

## Allele counting workflow configuration and execution

Allele counting is run in a two step procedure: (1) configuration and (2) workflow
execution. The configuration step is used to specify the input data and any
options pertaining to the allele counting methods. The execution
step is used to specify any parameters pertaining to _how_ the workflow is executed
(such as the total number of cores or SGE nodes over which the jobs should be
parallelized). The second execution step can also be interrupted and resumed
without changing the final result of the workflow.

### Configuration

The workflow is configured with the script: `${STRELKA_INSTALL_PATH}/libexec/configureSequenceAlleleCountsWorkflow.py`.
Running this script with no arguments will display all standard configuration
options to specify input alignment files, the reference sequence and the output run folder.
Note that all input alignment and reference sequence files must contain the same chromosome names
in the same order. Note that many downstream error modeling routines assume the counts are gathered from diploid regions, for now non-autosomes need to be filtered out crudely by explicitily listing all autosomes as region targets.

Example configuration for human autosomes:

    ${STRELKA_INSTALL_PATH}/libexec/configureSequenceAlleleCountsWorkflow.py \
    --bam=sample.bam \
    --referenceFasta=hg19.fa \
    --runDir ${COUNTS_ANALYSIS_PATH}
    --region chr1 \
    --region chr2 \
    --region chr3 \
    --region chr4 \
    --region chr5 \
    --region chr6 \
    --region chr7 \
    --region chr8 \
    --region chr9 \
    --region chr10 \
    --region chr11 \
    --region chr12 \
    --region chr13 \
    --region chr14 \
    --region chr15 \
    --region chr16 \
    --region chr17 \
    --region chr18 \
    --region chr19 \
    --region chr20 \
    --region chr21 \
    --region chr22


On completion, the configuration script will create the workflow run script `${COUNTS_ANALYSIS_PATH}/runWorkflow.py`. This can be used to run the workflow in various parallel compute modes per the
instructions in the [Execution](#execution) section below.

#### Configuration: Excluding regions

Regions of the genome specified in a BED file may be marked for exclusion from the counting process. Multiple BED files may be given, in which case the union of all excluded regions are skipped in the counting process

    ${STRELKA_INSTALL_PATH}/libexec/configureSequenceAlleleCountsWorkflow.py \
    --bam=sample.bam \
    --referenceFasta=hg19.fa \
    --runDir ${COUNTS_ANALYSIS_PATH} \
    --excludedRegions=nonconfidentRegions.bed.gz \
    --excludedRegions=nonAutosomes.bed.gz

#### Configuration: Annotating known variants

A _single_ known variants VCF file can be supplied to distinguish pattern counts of known and unknown variants.  Pattern analyzer output will then be labeled according to the number of times it overlaps a known variant present in the VCF.  The labels are currently "Unknown" for patterns that do not overlap a known variant call, and the matching genotype from the known variant VCF for those that do (Homref, Het, Hetalt, or Homalt).  If a known variant VCF is not provided, all patterns are marked as "Unknown".  This is currently limited to indels.

    ${STRELKA_INSTALL_PATH}/libexec/configureSequenceAlleleCountsWorkflow.py \
    --bam=sample.bam \
    --referenceFasta=hg19.fa \
    --runDir ${COUNTS_ANALYSIS_PATH} \
    --knownVariants=sample.knownVariant.vcf.gz

One strategy would be to combine the `--excludedRegion` argument to exclude non-platinum regions and the `--knownVariant` argument to label platinum variants in a Platinum Genomes sample.  This would provide "Unknown" patterns that either are in confident homref regions or do not match a platinum variant at a variant site, and known variant patterns that overlap high-confidence indel calls.

#### Configuration: Report observed indels

Extended context details for each observed indel can optionally be reported in a BED file.


    ${STRELKA_INSTALL_PATH}/libexec/configureSequenceAlleleCountsWorkflow.py \
    --bam=sample.bam \
    --referenceFasta=hg19.fa \
    --runDir ${COUNTS_ANALYSIS_PATH} \
    --reportObservedIndels

This argument will create a `debug` directory in `${COUNTS_ANALYSIS_PATH}/results` with a BED file titled `strelkaObservedIndel.bed.gz` where each observed indel as a single record.  In addition to the standard BED fields, it will also output the following additional fields:

| Field # |         Description                                               |
|--------:|:------------------------------------------------------------------|
|       4 | Indel type (INSERT/DELETE)                                        |
|       5 | Repeat unit                                                       |
|       6 | Reference repeat count (in # of repeat units)                     |
|       7 | Variant status (Unknown/Homref/Het/Hetalt/Homalt)                 |
|       8 | context repeat count (i.e. homopolymer length)                    |
|       9 | indel category (i.e. one of I\_3+, I\_2, I\_1, D\_1, D\_2, D\_3+) |
|      10 | Indel magnitude/length (in bp)                                    |
|      11 | Indel index at site (both 1-indexed)                              |
|      12 | Indel allele coverage                                             |
|      13 | Reference allele coverage                                         |
|      14 | Total locus coverage                                              |


### Execution

The configuration step creates a new workflow run script in the requested run directory:

`${COUNTS_ANALYSIS_PATH}/runWorkflow.py`

This script is used to control parallel execution of the workflow via the [pyFlow][2]
task engine. It can be used to parallelize analysis via one
of two modes:

1. Parallelized across multiple cores on a single node.
2. Parallelized across multiple nodes on an SGE cluster.

A running workflow can be interrupted at any time and resumed where it left
off. If desired, the resumed analysis can use a different running mode or total
core count.

For a full list of execution options, see:

`${COUNTS_ANALYSIS_PATH}/runWorkflow.py -h`

Example execution on a single node:

`${COUNTS_ANALYSIS_PATH}/runWorkflow.py -m local -j 8`

Example execution on an SGE cluster:

`${COUNTS_ANALYSIS_PATH}/runWorkflow.py -m sge -j 36`

### Advanced execution options

These options are useful for workflow development and debugging:

#### `--quiet`
Stderr logging can be disabled with `--quiet` argument. Note this log is replicated to `${COUNTS_ANALYSIS_PATH}/workspace/pyflow.data/logs/pyflow_log.txt` so there is no loss of log information.

## Viewing allele counting workflow ouput

### Summary output

The following command allows you to view a summary of the output from the allele counting workflow:

    ${STRELKA_INSTALL_PATH}/libexec/DumpSequenceAlleleCounts --counts-file ${COUNTS_ANALYSIS_PATH}/results/variants/strelkaAlleleCounts.bin

### Excluding basecalls/indels

DumpSequenceAlleleCounts reports both basecalls (i.e. SNPs) and indels.  If you would like to limit the pattern analyzer output to one of the two pattern types, you can do so with the `--exclude-basecalls` and `--exclude-indels` to exclude SNP and indel patterns, respectively.

### Extended output (for model development)

If you would like to get the complete contents of the file, which could be useful for model development in external software (e.g. numpy or R), you can get a tab-delimited file by running `DumpSequenceAlleleCounts` with the `--extended` argument.  This will write a headered tab-delimited file to stdout (N.B. output should __definitely__ be redirected to a file) with the following information:

| Field # | Field name     | Description                                      |
|--------:|---------------:|:-------------------------------------------------|
|       1 | context        | Currently, this is the homopolymer tract length  |
|       2 | variant_status | Is this a known variant? (UNKNOWN/VARIANT)       |
|       3 | alts           | Alt alleles observed (0 if no alts observed)     |
|       4 | alt_counts     | Counts for each alt allele (0 if no alts)        |
|       5 | ref_count      | Total reference allele coverage                  |
|       6 | total_alt      | Total alternate allele coverage                  |
|       7 | times_observed | Number of times this configuration is observed   |

Both alts and alt_counts are comma-delimited lists.  The first 6 fields define the specific configuration that was observed (e.g. 1 bp homopolymer with no observed alternate alleles, no known variants, and 30x coverage), while the final field provides the number of times it was observed.


## Error model parameter estimation

Given a counts file, error estimation under a particular model can be run via:

        ${STRELKA_INSTALL_PATH}/libexec/EstimateParametersFromAlleleCounts --counts-file myCounts.bin --model-type indel --model 2

Note that the estimation process currently offers very little runtime flexibility or documentation, any user of this tool is assumed to be adding or modifying models by changes the estimator source code itself. As a minimal convenience to developers, the model type (indel or snv) and the model index can be choosen at runtime, per the above example.

[2]: http://Illumina.github.io/pyflow/

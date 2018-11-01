Strelka User Guide
==================

## Table of Contents

[//]: # (BEGIN automated TOC section, any edits will be overwritten on next source refresh)

* [Introduction](#introduction)
* [Installation](#installation)
* [Method Overview](#method-overview)
* [Capabilities](#capabilities)
  * [Known Limitations](#known-limitations)
* [Input requirements](#input-requirements)
  * [Sequencing Data](#sequencing-data)
  * [Alignment Files](#alignment-files)
  * [VCF Files](#vcf-files)
* [Outputs](#outputs)
  * [Variant prediction](#variant-prediction)
    * [Germline](#germline)
    * [Somatic](#somatic)
  * [Statistics](#statistics)
* [Run configuration and Execution](#run-configuration-and-execution)
  * [Configuration](#configuration)
    * [Somatic configuration example](#somatic-configuration-example)
    * [Germline configuration example](#germline-configuration-example)
    * [General configuration options](#general-configuration-options)
    * [Advanced configuration options](#advanced-configuration-options)
  * [Execution](#execution)
    * [Advanced execution options](#advanced-execution-options)
  * [Extended use cases](#extended-use-cases)
    * [Improving runtime for references with many short contigs, such as GRCh38](#improving-runtime-for-references-with-many-short-contigs-such-as-grch38)
    * [Exome/Targeted](#exometargeted)
    * [De novo variant calling](#de-novo-variant-calling)
    * [RNA-Seq](#rna-seq)
    * [Heteroplasmic/pooled calling](#heteroplasmicpooled-calling)
    * [Somatic callability](#somatic-callability)
* [Special Topics](#special-topics)

[//]: # (END automated TOC section, any edits will be overwritten on next source refresh)


## Introduction

Strelka calls germline and somatic small variants from mapped sequencing reads. It is optimized for rapid clinical analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. Strelka's germline caller employs a haplotype model to improve call quality and provide short-range read-backed phasing in addition to a probabilistic variant calling model using indel error rates adaptively estimated from each input sample's sequencing data. Both germline and somatic callers include a final empirical variant rescoring step using a random forest model to reflect numerous features indicative of call reliability which may not be represented in the core variant calling probability model.

Strelka accepts input read mappings from BAM or CRAM files, and optionally candidate and/or forced-call alleles from VCF. It reports all small variant predictions in VCF 4.1 format. Germline variant reporting uses the [gVCF conventions][gvcfPage] to represent both variant and reference
call confidence. For best somatic indel performance, Strelka is designed to be run with the [Manta structural variant and indel caller][manta], which provides additional indel candidates up to a given maxiumum indel size (by default this is 49). By design, Manta and Strelka run together with default settings provide complete coverage over all indel sizes (in additional to all SVs and SNVs) for clinical somatic and germline analysis scenarios.

[manta]:https://github.com/Illumina/manta
[gvcfPage]:https://sites.google.com/site/gvcftools/home/about-gvcf


## Installation

Please see the [Strelka installation instructions](installation.md)

## Method Overview

The strelka workflow comprises a number of common sequence analysis
steps followed by application-specific variant modeling and empirical re-scoring
methods specific to the analysis of germline or somatic variation.

In all cases, strelka starts with a preliminary step to estimate various genomic or regional
statistics from the input alignments, including the sequence depth distribution
and, for germline analysis, a detailed analysis of indel error rates for each input sample.
This is followed by segmenting of the genome for parallel processing, where within each segment
all input samples are jointly analyzed to identify candidate alleles, realign all input reads,
analyze reads to make model specific variant inferences, then compute properties of each variant
used to apply filters or empirically recalibrate confidence that each variant represents a germline or somatic
variant in the input sample(s). Finally, all parallel segment results are joined to produce
Strelka's final variant output.

Strelka's germline calling model employs a haplotype representation to improve variant call quality
and provide short-range read backed phasing of all variants. The haplotype idenfiication method uses both
a fast k-mer ranking approach for simple loci and local assembly for more complex or repetitive regions.

## Capabilities

Strelka is capable of detecting SNVs and indels up to a predefined maximum size, currently
defaulting to 49 bases or less. Indels are detected from several sources, including indels
present in the input read alignments, indels detected by Strelka from the assembly of an
active region, and candidate indel VCFs provided as input from an external SV/indel caller or
population database. For somatic variant calling, it is a recommended best practice
to provide indel candidates from the [Manta SV and indel caller][manta], as outlined
in the suggested configuration steps below.

All methods are optimized by default for whole genome DNA-Seq, but are routinely tested for exome and amplicon inputs
(note additional flags for this case described in [Exome/Targeted](#exometargeted)).
RNA-Seq germline analysis is still in development and not fully supported. It
can be configured with the `--rna` flag, more details of this mode are described below in [RNA-Seq](#rna-seq).

Strelka's somatic calling capability is known to provide good results down to about 5-10% tumor purity given
sufficient normal and tumor sequencing depth. Strelka also accounts for minor contamination of the normal
sample with tumor cells (up to 10%) to better support liquid and late-stage solid tumor analysis.

Strelka is capable of performing joint germline analysis on a family scale (10s of samples). This is primarily intended
to facilitate de-novo variant analysis in families. Strelka's germline analysis capabilities are not currently
optimized for population analysis and may become unstable or fail to leverage population variant constraints to improve
calls at higher sample counts.

Strelka includes a short-range read-backed phasing capability for germline calls to facilitate
the correct inference of haplotypes induced by proximal SNVs and indels.


### Known Limitations

Strelka requires a matched normal sample to make somatic calls. The matched
normal is used to distinguish both germline variation and sequencing artifact from
somatic variation. The general depth guideline for the normal sample is either
one half the tumor depth or ~30x, whichever is higher.

As described above, Strelka's joint germline analysis is limited to family scale (10s of samples), and not intended
to support case/control or population analysis.

## Input requirements

### Sequencing Data

The input sequencing reads are expected to come from a paired-end sequencing assay.
Any input other than paired-end reads are ignored by default except to double-check
for putative somatic variant evidence in the normal sample during somatic variant analysis.
Read lengths above ~400 bases are not tested.

### Alignment Files

All input sequencing reads should be mapped by an external tool and provided as input in
BAM or CRAM format.

The following limitations apply to the input BAM/CRAM alignment records:

* Alignments cannot contain the "=" character in the SEQ field.
* RG (read group) tags are ignored -- each alignment file must represent one
  sample.
* Alignments with basecall quality values greater than 70 will trigger a runtime error (these
  are not supported on the assumption that the high basecall quality indicates an offset error)

### VCF Files

Input VCFs files are accepted for a number of roles as described below. All input VCF records are checked for
compatibility with the given reference genome, in additional to role-specific checks described below. If any
VCF record's REF field is not compatible with the reference genome a runtime error will be triggered.
"Compatible with the reference genome" means that each VCF record's REF base either (1) matches the corresponding
reference genome base or the VCF record's REF base is 'N' or the reference genome base is any ambiguous IUPAC base code
(all ambiguous base codes are converted to 'N' while importing the reference).

## Outputs

### Variant prediction

Primary variant inferences are provided as a series of [VCF 4.1][1] files in
`${STRELKA_ANALYSIS_PATH}/results/variants`.

#### Germline

Germline analysis is reported to the following variant files:

* __variants.vcf.gz__
    * This describes all potential variant loci across all samples. Note this file includes non-variant loci if they
    have a non-trivial level of variant evidence or contain one or more alleles for which genotyping has been forced.
    Please see the [multi-sample variants VCF](#interpreting-the-germline-multi-sample-variants-vcf) section
    below for additional details on interpreting this file.
* __genome.S${N}.vcf.gz__
    * This is the genome VCF output for sample ${N}, which includes both variant records and compressed non-variant blocks. The sample index, ${N} is 1-indexed and corresponds to the input order of alignment files on the configuration command-line.

##### Germline VCF Sample Names

Sample names printed into the VCF output are extracted from each input
alignment file from the first read group ('@RG') record found in the
header. Any spaces found in the name will be replaced with
underscores. If no sample name is found a default SAMPLE1, SAMPLE2,
etc.. label will be used instead.

##### Interpreting the germline multi-sample variants VCF

The germline multi-sample variants VCF, `variants.vcf.gz`,  describes all potential variant loci across all analyzed
samples. This VCF includes both high-confidence variant loci and lower-confidence potential variant loci, where
'high-confidence variant loci' refer to those that include a variant genotype passing all filters in at least one
sample. To ease interpretation of this file, an additional filter `NoPassedVariantGTs` is appended to the VCF `FILTER`
field at loci lacking any samples with a variant genotype passing all sample-level filters. This allows the
high-confidence variant loci to be queried by simply requiring that the VCF `FILTER` field is set to `PASS`. For
instance, ts/tv for high-confidence loci could be computed as follows:

```bash
bcftools stats -f PASS variants.vcf.gz | grep TSTV
```

To ease further analysis, the `NoPassedVariantGTs` filter can be updated using the script
`updateNoPassedVariantGTsFilter.py` provided in the strelka distribution. This can be useful to refresh the
`NoPassedVariantGTs` filter when the sample composition of the original file has changed, such as when extracting a
subset of samples. For example, updating this filter to apply only to samples NA12878 and NA12877 could be accomplished
as shown below:

```bash
bcftools view -s NA12878,NA12877 variants.vcf.gz |\
python ${STRELKA_INSTALL_PATH}/libexec/updateNoPassedVariantGTsFilter.py |\
bgzip -c >|\
sampleSubset.variants.vcf.gz
```

#### Somatic

Somatic analysis provides somatic variants in the following two files:

* __somatic.snvs.vcf.gz__
    * All somatic SNVs inferred in the tumor sample.
* __somatic.indels.vcf.gz__
    * All somatic indels inferred in the tumor sample.

The somatic variant caller can also optionally produce a callability track,
see the [somatic callability](#somatic-callability) section below for details.

##### Somatic variant allele frequencies

The somatic allele frequency estimate in the tumor sample is not directly available in the VCF output. A recommend way to extract such a value from the strelka VCF record is:

* Somatic SNVs:
```
refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FOMRAT/AU)
altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FOMRAT/TU)
tier1RefCounts = First comma-delimited value from $refCounts
tier1AltCounts = First comma-delimited value from $altCounts
Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
```

* Somatic indels:
```
tier1RefCounts = First comma-delimited value from FORMAT/TAR
tier1AltCounts = First comma-delimited value from FORMAT/TIR
Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)
```

### Statistics

Additional diagnostic output is provided in `${STRELKA_ANALYSIS_PATH}/results/stats`

* __genomeCallStats.tsv__
    * A tab-delimited report of various internal statistics from the variant calling process:
        * Runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging.
        * Indel candidacy statistics

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

All workflows are configured with a set of workflow scripts following the pattern:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelka${type}Workflow.py

Running one of these scripts with no arguments will display all standard configuration
options to specify input alignment files, the reference sequence and the output run folder.
Note that all input alignment and reference sequence files must contain the same chromosome names
in the same order. The default settings in all workflows assume a whole genome DNA-Seq analysis,
but there are configuration options for exome/targeted analysis, in addition to RNA-Seq options for the
germline workflow.

On completion, the configuration script will create the workflow run script `${STRELKA_ANALYSIS_PATH}/runWorkflow.py`.
This can be used to run the workflow in various parallel compute modes per the instructions in the [Execution](#execution) section
below.

For the somatic workflow, the best-practice recommendation is to run the [Manta SV and indel caller][manta] on the same set of
samples first, then supply Manta's candidate indels as input to Strelka. Examples of this procedure are shown below.

#### Somatic configuration example

Somatic analysis is configured with the script: `${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py`

Example Configuration:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam HCC1187BL.bam \
    --tumorBam HCC1187C.bam \
    --referenceFasta hg19.fa \
    --indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz \
    --runDir ${STRELKA_ANALYSIS_PATH}

In the above example, including the candidate indel file `${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz`
is a recommended best practice but not required. To generate these candidate indels the corresponding configuration for
Manta is:

    ${MANTA_INSTALL_PATH}/bin/configManta.py \
    --normalBam HCC1187BL.bam \
    --tumorBam HCC1187C.bam \
    --referenceFasta hg19.fa \
    --runDir ${MANTA_ANALYSIS_PATH}

...followed by execution of the manta workflow. Note that full installation and usage instructions for Manta can be found in the [Manta User Guide][mantaUserGuide].

[mantaUserGuide]:https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md

#### Germline configuration example

Germline analysis is configured with the script: `${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py`. The
germline caller already includes its own assembly step, so indel candidate input from a tool like Manta is not
recommended in this case (in testing, Manta candidate indels have been found to provide no benefit and in some cases
may complicate the accurate resolution of a locus).

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

#### General configuration options

##### Candidate indels

As demonstrated in the above examples, one or more candidate indel VCFs can be provided to any Strelka workflow with
the `--indelCandidates` flag. Any indels provided in this way will be given candidate status and considered during
realignment and genotyping steps, but will not be output unless the indel allele is found in the input samples.
This is a useful mechanism to supply Strelka with larger indels (ie. larger than can be found by the read mapper),
or to provide population variants.

In addition to [standard input VCF checks](#vcf-files), any candidate indel record which is not left-normalized will be skipped with a warning.

Multiple candidate indel VCFs may be submitted to the workflow (e.g. `--indelCandidates cand1.vcf.gz --indelCandidates cand2.vcf.gz ...`). All input VCFs must be bgzip compressed and tabix-indexed.

##### Forced genotypes

One or more forced genotype VCFs can be provided to any Strelka workflow with the `--forcedGT` configuration option.  Any indel allele provided in this way will be treated as a candidate (per the `--indelCandidates` option above), and additionally must appear in the output VCF, even when there is no support for the allele in the input samples. Be aware that in certain cases where a forced allele conflicts with Strelka's internal haplotype model, the forced variant will not be genotyped, but it will still appear in the VCF output with the "NotGenotyped" filter. Any SNV listed in the forced genotype VCF will prevent the corresponding site form being compressed into a homozygous reference block and ensure that a VCF site record is output for the given position, but will not provide any special treatment of the alternate base(s) listed in the VCF.

In addition to [standard input VCF checks](#vcf-files), any forced genotype variant record which is not left-normalized will trigger a runtime error.

Multiple forced genotype VCFs may be submitted to the workflow (e.g. `--forcedGT fgt1.vcf.gz --forcedGT fgt2.vcf.gz ...`). All input VCFs must be bgzip compressed and tabix-indexed.

##### Call regions

Strelka calls the entire genome by default, however variant calling may be restricted to an arbitrary subset of the genome by providing a region file in BED format with the `--callRegions` configuration option. The BED file must be bgzip-compressed and tabix-indexed, and only one such BED file may be specified. When specified, all VCF output is restricted to the provided call regions only, however statistics derived from the input data (such as expected chromosome depth) will not be restricted to the call regions. Note in particular that even when `--callRegions` is specified, the `--exome` flag is still required for exome or targeted data to get appropriate depth filtration behavior for non-WGS cases.

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
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12882	NA12878	NA12877
    chrX	0	.	N	<CNV>	.	PASS	END=10000	CN	1	2	1
    chrX	2781479	.	N	<CNV>	.	PASS	END=155701382	CN	1	2	1
    chrX	156030895	.	N	<CNV>	.	PASS	END=156040895	CN	1 	2	1
    chrY	0	.	N	<CNV>	.	PASS	END=57227415	CN	1	0	1

The span over which the copy number from each VCF record is applied is:  `[POS+1, INFO/END]`.

Strelka does not require any fields besides `CHROM`, `POS`, `ALT`, `INFO/END` and `FORMAT/CN`, so a ploidy specific record could be further simplified if desired, e.g:

    chrY	0	.	.	<CNV>	.	.	END=57227415	CN	1	0	1

...would be a valid input record for this option.

Note this feature is primarily intended to delineate the sex chromosome copy number but can be used to call small variants in the context of CNV calls as well.

###### Controlling gVCF homozygous reference block compression

The germline configuration option `--noCompress` can be used to specify a BED file listing sites/regions which will
be excluded from gVCF homozygous reference block compression. This can be used for various use cases such as excluding
an entire chromosome or all known polymorphic sites from homozygous reference block compression.

Only one BED file can be submitted to the `--noCompress` option. The BED file must be bgzip compressed and tabix-indexed.

Note that even if such an option was not selected, a tool such as gvcftools `break_blocks` can be used to uncompress
all gVCF regions specified in a BED file after germline calling has already been run.

###### Adaptive sequence error estimation

By default the germline calling workflow includes an adaptive sequencing error estimation step, in which indel errors are estimated from a subsample of the input sequencing data and used to estimate indel error rates for each sample. This step can be disabled with the configuration option `--disableSequenceErrorEstimation` so that the workflow reverts to precomputed indel error rates reflecting an intermediate point between different sequencing assays.

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

#### Improving runtime for references with many short contigs, such as GRCh38

For both germline and somatic analysis, Strelka may have runtime issues while attempting to process the large number of
small decoys and unplaced/unlocalized contigs found in GRCh38 and other reference genomes. This is due to a known issue
with read realignment sporadically experiencing substantial slowdowns on very short, high-depth contigs. Until this
issue can be resolved, runtime can be improved for such cases by excluding  smaller contigs from analysis. This can be
done in Strelka by creating a bed file of all the chromosomes that should be included in the analysis, and providing it
as an argument to the [call regions configuration option](#call-regions). For instance, the following bed file could be
provided for GRCh38 to exclude all decoys and small contigs:

```
chr1	0	248956422
chr2	0	242193529
chr3	0	198295559
chr4	0	190214555
chr5	0	181538259
chr6	0	170805979
chr7	0	159345973
chr8	0	145138636
chr9	0	138394717
chr10	0	133797422
chr11	0	135086622
chr12	0	133275309
chr13	0	114364328
chr14	0	107043718
chr15	0	101991189
chr16	0	90338345
chr17	0	83257441
chr18	0	80373285
chr19	0	58617616
chr20	0	64444167
chr21	0	46709983
chr22	0	50818468
chrX	0	156040895
chrY	0	57227415
chrM	0	16569
```

#### Exome/Targeted

Supplying the `--exome` flag at configuration time will provide appropriate settings for WES and other
regional enrichment analyses. At present this flag disables all high depth filters, which are designed
to exclude pericentromeric reference compressions in the WGS case but cannot be applied correctly to a
targeted analysis.

For germline analysis, this mode also disables the empirical variant scoring (EVS) model, falling back
to a set of simple threshold based filters instead. The somatic EVS model remains in use for exome and
targeted data.

In an exome or targeted analysis it may be desirable to restrict calling to the targeted regions by
providing a BED file to the `--callRegions` option. Note that this option acts independently of `--exome`.

#### De novo variant calling

The recommended workflow for de novo variant calling with strelka is to first run a multi-sample germline analysis on the sample pedigree of interest, and then analyze the resulting multi-sample VCF output from strelka with a separate de novo variant calling script, `denovo.py`. This script infers the de novo event probability for each variant in the proband. The details of this procedure are described in [de novo variant calling using Strelka](deNovoQualityScore.md).

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

* [De novo variant calling using Strelka](deNovoQualityScore.md)
* [Training Procedure for Somatic Empirical Score](trainingSomaticEmpiricalScore.md)
* [Training Procedure for Germline Empirical Score](trainingGermlineEmpiricalScore.md)


[1]: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
[2]: http://Illumina.github.io/pyflow/

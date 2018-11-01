# De Novo Quality Score for Small Pedigrees

[User Guide Home](README.md)

## Introduction

This document outlines the recommended process for calling de novo variants from Strelka multi-sample VCFs.
This process uses a separate script, `denovo.py` implementing a DeNovoGear-like model. This script can be
used to assign de novo quality scores to each variant in the VCF given a specified sample pedigree. Note
that this procedure is provided to document the recommended process for Strelka, but it is not fully
integrated into Strelka and therefore requires additional package requirements as specified below.

## Requirements

The de novo quality score calculation script, `denovo.py` has the following requirements:

* python 2.7
* numpy >= 1.11.1
* cyvcf2 >= 0.6.2

Except for python, these are not part of the standard Strelka package requirements and are also not
provided by the Strelka build system or binary package distribution. A template for the installation of
all required packages to run `denovo.py` is provided in the following [Dockerfile](../../src/python/deNovoQualityScore/Dockerfile),
which can be used to either setup a docker image or as a guideline to install all dependencies on another system.

## Example

A small example multi-sample VCF file from Strelka is provided as an example to demonstrate the de novo
variant calling method. This VCF contains a short segment from chromosome 20 for a sample trio from CEPH
pedigree 1463. To artificially create de novo variant calls for demonstration purposes,
it can be run through `denovo.py` with an incorrect pedigree as follows:

```
${STRELKA_INSTALL_PATH}/share/deNovoQualityScore/denovo.py \
  --proband NA12877 --mother NA12878 --father NA12882 \
  ${STRELKA_INSTALL_PATH}/share/deNovoQualityScore/example/strelka_inverted_trio_segment.vcf.gz |\
  ${STRELKA_INSTALL_PATH}/libexec/bgzip -c > strelka_inverted_trio_segment-denovo-scores.vcf.gz
```

Here, the `proband`, `mother` and `father` arguments specify the sample names in the VCF. For a conventional de novo
 calling application these should be set to match the respective family members. For the purpose of creating a small
 demonstration, the sample assignments of the father and proband are inverted in the above command-line.

The returned VCF file contains an additional de novo quality score format field, `DQ`, for the proband sample. High
`DQ` scores indicate evidence of a de novo event in the proband. De novo variants are identified by `DQ` exceeding a
user-defined threshold. The best threshold choice will be application dependent, the recommended default threshold
is 7. Using the processed data from the example above, a command-line to select de novo variants using this default
threshold is:

```
bcftools view -i 'FORMAT/DQ >= 7' strelka_inverted_trio_segment-denovo-scores.vcf.gz
```

...which is expected to show 2 de novo SNV records.


## Interface, Data and Parameters

Abstract example for invoking denovo.py:

```
${STRELKA_INSTALL_PATH}/share/deNovoQualityScore/denovo.py [...] \
  --proband <probandID> \
  --mother <motherID> \
  --father <fatherID> \
  <input.vcf>
```

### Input VCF

All variants are read from an input VCF file which contains the trio
encoded as different samples. Minimal requirements for this file are:
- Follows the VCF v4.1 or v4.2 specification
- Contains a PL genotype field for each samples

Further remarks:
* Input file can be uncompressed, gzipped or bgzipped; all input formats
  supported by HTSLIB are covered. Instead of a path to a file, STDIN as stream
  input is also supported.
* Additional samples in the file will be silently ignored in the analysis.

*Command line argument*: Path to file or STDIN (denoted by `-`), required


### Output VCF

The primary output of the analysis is a VCF file which contains the same
information of as the input VCF, and a `DQ` FORMAT field if a meaningful score
could be computed. The proband sample contains the `DQ` score, in the other
samples it is marked as missing. The `DQ` field is defined in the header as

```
##FORMAT=<ID=DQ,Number=1,Type=Float,Description="Denovo quality">
```

and represents a score of the posterior probability of the variant being denovo
in the proband. While the underlying posterior probability p falls naturally
in the interval [0, 1], it is transformed to the final DQ score with
-log_{10} p. Thus, large scores are indicative of denovo variants. For
example, DQ scores of 13 and 20 would correspond to a posterior probability of a
denovo variant of 0.95 and 0.99, respectively. The reported score is rounded to
a defined number of decimals, by default to one decimal.

The VCF also receives a header line that specifies the denovo program version
used for the analysis:

```
##denovo_program=denovo.py 1.2.3
```

The VCF is written to stdout and can be piped into other tools, such as `bgzip` or `bcftools`.


### Pedigree information

The pedigree defines the role of a each sample in the trio of interest. This
information is supplied as individual samples IDs:

```
${STRELKA_INSTALL_PATH}/share/deNovoQualityScore/denovo.py --proband NA12882 --mother NA12878 --father NA12877 <input.vcf> <output.vcf>
```

In case the pedigree information is incomplete or ambiguous, an exception is
raised.

### Parameters

#### Model

`${STRELKA_INSTALL_PATH}/share/deNovoQualityScore/denovo.py` uses the DeNovoGear
model, as described in [Ramu et al.,
2013](http://www.nature.com/nmeth/journal/v10/n10/full/nmeth.2611.html) and
[DeNovoGear](https://github.com/denovogear/denovogear/). This estimates the
posterior probability of the variant being denovo given the genotype likelihoods
of the trio. The implementation follows the DeNovoGear codebase on github.
Scores derived from the denovo posterior probabilities are reported in the VCF.
Variants for which a concrete and meaningful DQ score cannot be calculated get
no DQ field.


### Priors

The analysis makes use of allele-specific priors. Currently, the precomputed
priors are stored in `${STRELKA_INSTALL_PATH}/share/deNovoQualityScore/data/prior` and need to be distributed with the script. A detailed description of the computation of the prior and its
parameters can be found in the denovogear publication.

Column descriptions for SNV priors (`snv_lookup_*.txt`):

| n_u_alleles              | case            | transmission      | gt                                     | mrate         | denovo_flag              | normal_flag                  | A / C / G / T prior (4 columns)      |
|--------------------------|:---------------:|:-----------------:|:--------------------------------------:|:-------------:|:------------------------:|:----------------------------:|:------------------------------------:|
| Number of unique alleles | Case identifier | Transmission rate | Trio genotypes (proband/mother/father) | Mutation rate | Flag for denovo variants | Flag for non-denovo variants | Transition probabilities for alleles |

Column descriptions for SNV priors (=indel_lookup_*.txt=):

| n_u_alleles              | case            | transmission      | gt                                                           | mrate         | denovo_flag              | normal_flag                  | prior                                        |
|--------------------------|:---------------:|:-----------------:|:------------------------------------------------------------:|:-------------:|:------------------------:|:----------------------------:|:--------------------------------------------:|
| Number of unique alleles | Case identifier | Transmission rate | Trio genotypes (R: reference, V: alt; proband/mother/father) | Mutation rate | Flag for denovo variants | Flag for non-denovo variants | Transition probabilities for variant alleles |


## Parallel Execution

Parallel execution can be achieved with [GNU
parallel](https://www.gnu.org/software/parallel/) by setting the `--parallel`
flag. In addition, the `--ncores` argument allows to specify the number of cpus
used for the calculation.

```
python ${STRELKA_INSTALL_PATH}/share/deNovoQualityScore/denovo.py \
  --proband PROBAND \
  --mother PARENT1 \
  --father PARENT2 \
  --parallel \
  input.vcf.gz \
  | bgzip -c > output.vcf.gz
```

# Installation and Testing

## Using docker

The [`Dockerfile`](../../src/python/deNovoQualityScore/Dockerfile) creates a container with `denovo.py` installed and runs the tests. This should be the preferred way of running the application.


## Using python

### Requirements

#### Install Requirements

* `python 2.7`
* `numpy`
* `cyvcf2`

Install/run requirements are listed in [requirements.txt](../../src/python/deNovoQualityScore/requirements.txt), and call be installed with

```
pip install -r requirements.txt
```

#### Test Requirements

* `pytest`
* `pytest-cov`
* `PyVCF`

Test requirements are listed in [test/test_requirements.txt](../../src/python/deNovoQualityScore/test/test_requirements.txt), and call be installed with

```
pip install -r test/test_requirements.txt
```

### Running the script

```
./denovo.py
```

or equivalently

```
python denovo.py
```

### Running the test suite

```
python -m pytest test_denovo.py
```

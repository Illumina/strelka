Strelka2 Small Variant Caller
============================

Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. The germline caller employs an efficient tiered haplotype model to improve accuracy and provide read-backed phasing, adaptively selecting between assembly and a faster alignment-based haplotyping approach at each variant locus. The germline caller also analyzes input sequencing data using a mixture-model indel error estimation method to improve robustness to indel noise. The somatic calling model improves on the original Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell contamination in the normal sample. A final empirical variant re-scoring step using random forest models trained on various call quality features has been added to both callers to further improve precision.

Compared with submissions to the recent PrecisonFDA Consistency and Truth challenges, the average indel F-score for Strelka2 running in its default configuration is 3.1% and 0.08% higher, respectively, than the best challenge submissions. Runtime on a 28-core server is ~40 minutes for 40x WGS germline analysis and ~3 hours for a 110x/40x WGS tumor-normal somatic analysis. More details on Strelka2 methods and benchmarking for both germline and somatic calling are described in the following open-access pre-print:

Kim, S., Scheffler, K. *et al.* (2017) Strelka2: Fast and accurate variant calling for clinical sequencing applications. *bioRxiv* [doi: 10.1101/192872][preprint]

Strelka accepts input read mappings from BAM or CRAM files, and optionally candidate and/or forced-call alleles from VCF. It reports all small variant predictions in VCF 4.1 format. Germline variant reporting uses the [gVCF conventions][gvcfPage] to represent both variant and reference
call confidence. For best somatic indel performance, Strelka is designed to be run with the [Manta structural variant and indel caller][manta], which provides additional indel candidates up to a given maxiumum indel size (49 by default). By design, Manta and Strelka run together with default settings provide complete coverage over all indel sizes (in additional to SVs and SNVs). See the [user guide][UserGuide] for a full description of capabilities and limitations.

[preprint]:http://dx.doi.org/10.1101/192872
[gvcfPage]:https://sites.google.com/site/gvcftools/home/about-gvcf
[manta]:https://github.com/Illumina/manta
[UserGuide]:docs/userGuide/README.md

Quick Start
---------------

### Installation
For Linux users, it is recommended to start from the most recent
[binary distribution on the Strelka releases page][releases] (replace x.y.z with a Strelka version):
```bash
wget https://github.com/Illumina/strelka/releases/download/v2.9.1/strelka-x.y.z.centos6_x86_64.tar.bz2
tar xvjf strelka-x.y.z.centos6_x86_64.tar.bz2
```
This distribution can be unpacked, moved to any convenient directory and
tested by [running a small demo](docs/userGuide/installation.md#demo)
included with the release distribution. Strelka can also be installed
from source code. Please see the [installation instructions](docs/userGuide/installation.md)
for full build and installation details.

[releases]:https://github.com/Illumina/strelka/releases


### Run configuration and execution

Strelka is run in two steps: (1) configuration (specifying input data and options) and 
(2) workflow execution (specifying parameters on how strelka is executed). The second execution step can also be interrupted and restarted without changing the final result of the workflow. 

#### Configuration
Example configuration for germline calling:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${ALIGNMENT_FILE1}.bam \
    --bam ${ALIGNMENT_FILE2}.bam \
    --bam ${ALIGNMENT_FILE3}.bam \
    --ref ${REFERENCE}.fa \
    --runDir ${STRELKA_ANALYSIS_PATH}

Example configuration for somatic calling:

    ${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${NORMAL}.bam \
    --tumorBam ${TUMOR}.bam \
    --ref ${REFERENCE}.fa \
    --runDir ${STRELKA_ANALYSIS_PATH}

[excludeContigs]:https://git.illumina.com/Bioinformatics/strelkadev/blob/develop/docs/userGuide/README.md#improving-runtime-for-references-with-many-short-contigs-such-as-grch38

For references with many short contigs, it is strongly recommended to 
[provide callable regions to avoid possible runtime issues][excludeContigs]:

    --callRegions ${CALLABLE_REGION_FILE}.bed.gz 

[mantaCandidates]: https://git.illumina.com/Bioinformatics/strelkadev/blob/develop/docs/userGuide/README.md#somatic-configuration-example

For somatic calling, it is recommended to [provide indel candidates from the Manta SV and indel caller][mantaCandidates]
to improve sensitivity to call indels of size larger than 20: 

    --indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz

For a full list of execution options, see:

    ${STRELKA_INSTALL_PATH}/configureStrelkaGermlineWorkflow.py -h

#### Workflow execution

Example execution on a single node:

    ${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m local -j $NUM_JOBS

Example execution on an SGE cluster:

    ${STRELKA_ANALYSIS_PATH}/runWorkflow.py -m sge -j $NUM_JOBS

[runWorkflowUserGuide]: https://git.illumina.com/Bioinformatics/strelkadev/blob/develop/docs/userGuide/README.md#execution

For more details see [the Execution section in the user guide][runWorkflowUserGuide].

Refer to the [Strelka user guide][UserGuide] for full instructions on how to run Strelka, 
interpret results and estimate hardware requirements/compute cost, 
in addition to a high-level methods overview.

License
-------

Strelka source code is provided under the [GPLv3 license](LICENSE.txt).
Strelka includes several third party packages provided under other
open source licenses, please see [COPYRIGHT.txt](COPYRIGHT.txt)
for additional details.



Strelka Code Development
------------------------

For strelka code development and debugging details, see the
[Strelka developer guide][DeveloperGuide]. This includes details
on Strelka's development protocols, special build instructions,
recommended workflows for investigating
calls, and internal documentation details.

[DeveloperGuide]:docs/developerGuide/README.md

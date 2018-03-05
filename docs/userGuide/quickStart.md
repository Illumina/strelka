Strelka2 Quick Start
====================

[releases]:https://github.com/Illumina/strelka/releases

### Installation
For Linux users, it is recommended to start from the most recent
[binary distribution on the Strelka releases page][releases] (replace x.y.z with a Strelka version):
```bash
wget https://github.com/Illumina/strelka/releases/download/v2.9.1/strelka-x.y.z.centos6_x86_64.tar.bz2
tar xvjf strelka-x.y.z.centos6_x86_64.tar.bz2
```
Strelka can also be installed from source code. Please see the [installation instructions](docs/userGuide/installation.md)
for full build and installation details.

### Configuration and execution

Strelka is run in two steps: (1) configuration (specifying input data and options) and 
(2) workflow execution (specifying parameters on how strelka is executed). The second execution step can also be interrupted and restarted without changing the final result of the workflow. 

Example for germline calling:

    strelka-x.y.z.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
        --bam strelka-x.y.z.centos6_x86_64/share/demo/strelka/data/NA12891_demo20.bam \
        --bam strelka-x.y.z.centos6_x86_64/share/demo/strelka/data/NA12892_demo20.bam \
        --ref strelka-x.y.z.centos6_x86_64/share/demo/strelka/data/demo20.fa \
        --runDir demo_germline
    demo_germline/runWorkflow.py -m local -j 1

Example for somatic calling:

    strelka-x.y.z.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
        --normalBam strelka-x.y.z.centos6_x86_64/share/demo/strelka/data/NA12891_demo20.bam \
        --tumorBam strelka-x.y.z.centos6_x86_64/share/demo/strelka/data/NA12892_demo20.bam \
        --ref strelka-x.y.z.centos6_x86_64/share/demo/strelka/data/demo20.fa \
        --runDir demo_somatic
    demo_somatic/runWorkflow.py -m local -j 1

[excludeContigs]:README.md#improving-runtime-for-references-with-many-short-contigs-such-as-grch38
[mantaCandidates]: README.md#somatic-configuration-example

### Tips

For references with many short contigs, it is strongly recommended to 
[provide callable regions to avoid possible runtime issues][excludeContigs]:

    --callRegions callable.bed.gz 


For somatic calling, it is recommended to [provide indel candidates from the Manta SV and indel caller][mantaCandidates]
to improve sensitivity to call indels of size 20 or larger: 

    --indelCandidates candidateSmallIndels.vcf.gz

For exome and amplicon inputs, add:

    --exome

### User guide

[UserGuide]: README.md
Refer to the [Strelka user guide][UserGuide] for full instructions on how to run Strelka, 
interpret results and estimate hardware requirements/compute cost, 
in addition to a high-level methods overview.

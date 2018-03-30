Strelka Quick Start
====================

[releases]:https://github.com/Illumina/strelka/releases

### Installation
For Linux users, it is recommended to start from the most recent
[binary distribution on the Strelka releases page][releases], this
distribution can be unpacked, moved to any convenient directory and
tested by [running a small demo](installation.md#demo)
included with the release distribution. For example, unpacking and running
the demo on the strelka 2.9.2 binary distribution could be accomplished as follows
(optionally replace 2.9.2 with a different Strelka version):
```bash
# download strelka binary
wget https://github.com/Illumina/strelka/releases/download/v2.9.2/strelka-2.9.2.centos6_x86_64.tar.bz2
# decompress
tar xvjf strelka-2.9.2.centos6_x86_64.tar.bz2
# run demo to check successful installation
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaSomaticWorkflowDemo.bash
bash strelka-2.9.2.centos6_x86_64/bin/runStrelkaGermlineWorkflowDemo.bash
```

Strelka can also be installed from source code. Please see the [installation instructions](installation.md)
for full build and installation details.

### Configuration and execution

Strelka is run in two steps: (1) configuration (specifying input data and options) and
(2) workflow execution (specifying parameters on how strelka is executed). The second execution step can also be interrupted and restarted without changing the final result of the workflow.

Example for germline calling:
```bash
# configuration
${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
    --bam sample1.bam \
    --bam sample2.bam \
    --ref hg38.fa \
    --runDir demo_germline
# execution on a single local machine with 20 parallel jobs
demo_germline/runWorkflow.py -m local -j 20
```

Example for somatic calling:
```bash
# configuration
${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam normal.bam \
    --tumorBam tumor.bam \
    --ref hg38.fa \
    --runDir demo_somatic
# execution on a single local machine with 20 parallel jobs
demo_somatic/runWorkflow.py -m local -j 20
```

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

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
[QuickStart]:docs/userGuide/quickStart.md
[UserGuide]:docs/userGuide/README.md

Getting Started
---------------
To get started installing and using Strelka, please consult the [quick start guide][QuickStart].

Data Analysis and Interpretation
---------------
After completing installation and reviewing the quick start guide, see the [Strelka user guide][UserGuide] for full instructions on how to run Strelka, interpret results and estimate hardware requirements/compute cost, in addition to a high-level methods overview.

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

## v2.9.10 - 2018-11-07

This is a minor feature update from v2.9.9 which adds a method that computes de novo variant quality scores.

### Added
- Add recommended de novo variant calling method and documentation (STREL-983)
  - The de novo variant quality score computation method, `denovo.py`, has been previously developed and maintained as a separate project. It is now added to the primary strelka distribution to make this easier for users to access.

## v2.9.9 - 2018-09-13

This is a minor update from v2.9.8.

### Changed
- Changed EVS threshold for RNA-seq inputs (STREL-975)
  - The EVS threshold was changed from 15 to 5 for both indels and SNVs.

## v2.9.8 - 2018-09-06

This is a minor bugfix update from v2.9.7.

### Fixed
- Fix depth filters for germline continuous frequency indel calls (STREL-978)
  - The LowDepth filter was being spuriously applied to all indels called through the continuous frequency model. This model is typically applied only to the mitochondrial chromosome to find heteroplasmic calls.
  - This filter is now fixed to work as documented for all variant types.

## v2.9.7 - 2018-08-14

This is a minor bugfix update from v2.9.6.

### Fixed
- Preserve softlinks in input reference fasta file paths (STREL-976)

## v2.9.6 - 2018-07-17

This is a minor bugfix update from v2.9.5.

### Fixed
- Fix realignment error when processing certain RNA-seq inputs (STREL-966)
  - Recent changes to read realignment/scoring introduced an edge case which fails on RNA-Seq reads from certain aligners - specifically when an insertion immediately follows an intron gap. The issue is now fixed.

## v2.9.5 - 2018-06-25

This is a minor bugfix update from v2.9.4.

### Fixed
- Fix germline joint sample failure case (STREL-955)
  - This fixes an uncommon issue where a multi-sample germline analysis qualifies for indel error rate estimation, but the estimator fails to converge for at least one sample. When this happens, the method is supposed to fallback to a set of default error parameters but was failing to do so. This is now fixed. The issue does not impact single-sample analysis or the somatic workflow.

## v2.9.4 - 2018-06-08

This is a minor bugfix update from v2.9.3.

### Fixed
- Fix rare instance where FORMAT/DP50, FORMAT/FDP50 and FORMAT/SUBDP50 are '-nan' (STREL-936)

## v2.9.3 - 2018-05-07

### Changed
- Add strict checks and improve error message for BED regions of size less than one (STREL-865/[#36])

### Fixed
- Prevent an assertion for certain rarely encountered combinations of forced genotype alleles in the germline workflow (STREL-896)

## v2.9.2 - 2018-03-02

This is a minor bugfix update from v2.9.1.

### Changed
- Turn on automated task resubmission for all workflow run modes (STREL-852)
  - Failed tasks have always been automatically resubmitted in SGE mode, this is now extended to localhost mode as well.
  - This change is intended to work around sporadic I/O issues on network filesystems.

### Fixed
- Update to pyflow v1.1.20 to close infrequent race condition in task resolution (STREL-853)
  - Under the race condition error, a non-failing task would be logged as failing with the message "TASKNAME has stopped without a traceable cause"

## v2.9.1 - 2018-02-28

This is a minor bugfix update from v2.9.0.

### Fixed
- Standardize header description of somatic snv and indel LowDepth filters (STREL-849)
  - This allows merging of the snv and indel vcf files.
- Fix build for gcc7+ and boost-1.64+ (STREL-824)

## v2.9.0 - 2018-02-08

This is a major update from v2.8.4. The most important change in this release is indirect: haplotype modeling and realignment have been improved such that, given the strelka2 germline VCF output of a trio at typical cWGS depth, false-positive de novo variant calls have been roughly cut in half. This is due to fixes for realignment artifacts that were too rare to noticeably impact germline call quality, but frequent relative to de novo variant rates. These changes should also accelerate the future transition to haplotype modeling for somatic variants. Many additional improvements to stability, error diagnostics, ease of use and accuracy are also included in this release, as enumerated below.

### Added
- Add strand bias feature to germline indel scoring model (STREL-676)
  - Improves filtration for a small number of false positive indels in typical WGS analysis.
- Add haplotyping constraints to the read alignment (STREL-743)
  - Phasing information from haplotyping is used to constrain combinations of variants within read alignments
  - Removes rare artifact which could trigger false de novo calls from multi-sample germline variant output, baseline false positive SNVs and indels reduced to approx half of previous count.
- Add new filter to make multi-sample germline variant output easier to interpret (STREL-819)
  - Locus filter 'NoPassedVariantGTs' added when no sample has a passing variant genotype.
  - This allows passing variants to be easily extracted with the FILTER field, without querying FORMAT/GT and FORMAT/FT.
- Add new filter to prevent interference between forced indels and other indels (STREL-607)
  - Locus filter 'NotGenotyped' added to ForcedGT indels if they can possibly interfere with indels discovered by Strelka.
  - All complex alleles are also not genotyped and appear in the VCF output with this NotGenotyped filter.

### Changed
- Change default maximum indel size from 50 to 49 (STREL-811)
  - This change is made as part of an effort to better align manta with GIAB SV size range conventions, such that strelka and manta together provide complete, non-overlapping coverage over the full indel spectrum using default settings.
- Remove preliminary step which counts the 'mappable' (non-N) size of the genome (STREL-772)
  - This has a legacy use in identifying noisy alignments. Now replaced with a simplified scheme.
- Lower default local task memory requirement from 2 to 1.5 Gb (STREL-802)
  - This enables all cores on an AWS c4.8xlarge with default configuration, use `--callMemMb` option to override for unusual cases.
- Update LowDepth filter for somatic calls to include cases where the normal sample depth is below 2 (STREL-745)
- Update htslib to incorporate CRAM file query fix (STREL-839/MANTA-1336)
  - This is expected to resolve possible issues with error parameter estimation from alignments in CRAM format.

### Fixed
- Fix empirical variant scoring (EVS) of complex somatic indels (STREL-774)
  - Previously the EVS model was overly pessimistic against complex somatic indels. This is now fixed by changing how EVS input features are computed for complex indels.
- Fix default sample name used in the VCF output for germline analysis (STREL-737)
  - Default is used when sample name cannot be parsed from the BAM header. Now fixed to insert SAMPLE1, SAMPLE2, etc. as documented.
- Fix rare instance where strand bias (FORMAT/SB) is 'inf' (STREL-741)
- Provide clear error message when attempting to configure/run with python3 (STREL-762)
- Fix python configure scripts to make maximum reported indel size configurable (STREL-763)
  - This can be done by configuring the maxIndelSize value inside the .ini file.
- Fix realignment slow-down issue that occurs when reads overlap too many candidate SNVs (STREL-805)
- Stop automatically clearing python environment variables (STREL-810)
  - This should allow python from certain module systems to be used, but may (rarely) cause instability due to conflicting content in a user's PYTHONPATH.
- Standardize germline FORMAT/GQ VCF tag to Integer type (STREL-812)
- Fix the issue that low depth filter is not applied to continuous variant frequency (e.g. mitochondrial) calls (STREL-803)

## v2.8.4 - 2017-10-23

This is a major bugfix update from v2.8.3. The two most notable changes are (1) a nearly 10-fold
reduction in memory usage during germline analysis and (2) fixing the default variant recall level
for male non-PAR chrX, or any region given a ploidy of 1 in the ploidy VCF file.

### Changed
- Switch to RapidJSON library for all json parsing (STREL-696)
  - Reduces germline calling memory usage ~10-fold due to improved parse of random forest rescoring models.
- Change active region detection method to create active regions shared by all samples (STREL-710)
- Verify region/callRegion values at configuration time (STREL-724)
  - Chromosome labels in BED records and region arguments must be found in the reference.
- Verify run directory has not already been configured (MANTA-1252/STREL-734)
- Verify alignment file extension at configuration time (MANTA-886)
- Update minimum supported linux OS from Centos 5 to Centos 6 (STREL-720)
- Move changelog to markdown format (STREL-571)

### Fixed
- Fix germline empirical variant scoring (EVS) for haploid regions (STREL-678)
  - Previously, EVS resulted in reduced recall for haploid regions such as non-PAR regions of chrX in male samples.
    After adding haploid training examples from NA12877 chrX, EVS performance for haploid regions is comparable to diploid.
- Fix debug option to provide realigned reads in bam output (STREL-721/[#15])

## v2.8.3 - 2017-09-22

This is a bugfix update from v2.8.2.

### Fixed
- Make minor correction to the non-error term used during adaptive indel error estimation (STREL-705)
- Make minor correction to somatic joint allele-frequency prior (STREL-632)
- Improve somatic EVS feature consistency (STREL-652)
- Improve CRAM reference handling (STREL-647)
  - The reference provided as input during workflow configuration is now prioritized over the URI in the CRAM header. This makes it easier to work with any CRAM file which contains a local file path in the header.

## v2.8.2 - 2017-08-03

This is a minor bugfix update from v2.8.1.

### Fixed
- Fix haplotype model issue occurring when contigs have no sequence coverage (STREL-653)

## v2.8.1 - 2017-08-02

This is minor bugfix update from v2.8.0.

### Fixed
- Fix allele noise filtration to synchronize across multiple samples (STREL-650)
- Fix minor inconsistency in sequence error counting genome segment bounds (STREL-651)
- Update htslib/samtools to v1.5 for improved error detection/messages and CRAM support (STREL-633)

## v2.8.0 - 2017-07-14

This is a major feature update from v2.7.1.

- STREL-608 Fix hang after error during adaptive estimation
- STREL-610 Fix workflow resumption after interrupt
- STREL-557 Retrain somatic EVS on updated alignments and truth sets
- STREL-609 Fix genotyping error on forced non-variant indels
- STREL-612 Fix sites overlapping forced non-variant indels
- STREL-602 Retrain germline SNV EVS on updated alignments
- STREL-597 Add low depth filter for both somatic and germline calls
- STREL-596 Add missing data to pooled indel calls
- STREL-586 Improve germline multi-sample calling runtime
- STREL-576 Reduce demo installation size
- STREL-579 Add filtered depth rates as somatic SNV EVS features
- STREL-580 Retrain germline EVS
- STREL-260 Add validation on all input vcf reference fields
- STREL-566 Use indel error estimation by default for germline calling
- STREL-577 Fix bug creating negative size active regions
- STREL-553 Standardize somatic EVS features
- STREL-564 Add filter preventing low depth PASS calls
- STREL-567 Change readConfidentSupportThreshold to 0.51
- STREL-524 Retrain germline EVS
- STREL-519 Fix callRegions option thread utilization
- STREL-459 Retain optimal soft-clipping for RNA analysis
- STREL-451 Enable somatic indel EVS and retrain somatic SNV EVS
- STREL-478 Change gVCF non-variant blocks to use mean depth
- STREL-469 Add RNAseq EVS models
- STREL-479 External candidate indels create active regions
- STREL-462 Do not penalize candidate SNVs in alignment scoring
- STREL-450 Add dinucs to indel error stats module
- STREL-465 Filter germline haplotypes with phasing error signature
- STREL-460 Do not assess indel candidacy in active regions if assembly fails
- STREL-454 Penalize for non-candidate indels in alignment scoring
- STREL-443 Remove old germline caller target BED file option
- STREL-178 Replace codon phaser with variant phaser
- STREL-356 Allow BED file to restrict variant calling regions
- STREL-392 Fix overlap del handling across segments
- STREL-405 Fix overcounting issue in error stats module
- STREL-401 Remove read edge events from error stats module
- STREL-342 Use assembly to generate haplotypes in long active regions
- STREL-251 Enable automatic germline EVS calibration
- STREL-357 Add separate threshold for homref calls
## v2.7.1
- STREL-336 Fix incorrect indel normalizations
- STREL-332 Revert to core somatic indel scoring with adjusted threshold
- STREL-331 Left-shift indels inside of active regions
- STREL-248 Update somatic EVS to include isaac4 data in training
- STREL-321 Fix rare shared variant prefix issue in multi-sample analysis
- STREL-322 Fix VCF namespace conflict on INFO/EVS
- STREL-200 Adjust mitochondrial strand-bias for consistency with autosomes
- STREL-319 Sync RNA-seq EVS feature names with DNA changes
- STREL-318 Accept VCFs with unknown ALT values as forced-call sites
- STREL-275 REF and ALTs do not share a common prefix of size over 1
- STREL-274 Simplify variant filter logic in multi-sample germline gVCF(s)
- STREL-269 Prevent indels larger than the max indel size from becoming candidate
- STREL-267 Forced indel call does not appear in somatic output
## v2.7.0
- STREL-184 Enable somatic indel EVS
- STREL-228 Improve runtime for references with many small contigs
- STREL-210 Retrain germline EVS
- STREL-164 Document short haplotyping
- STREL-254 Fix haploid indel call exception
- STREL-227 Shift away from legacy naming schemes
- STREL-240 Make active regions non-overlapping
- STREL-239 Fix haploid nonvariant indels in multi-sample analysis
- STREL-238 Add multi-sample ploidy specification using vcf
- STREL-217 Integrate short haplotying with multi-sample support
- STREL-158 Support multi-sample in standard germline analysis
- STREL-218 Tune indel parameters for improved het/hom ratio
- STREL-216 Relax normalization requirement for indel candidates to a warning
- STREL-198 Update germline EVS using mixed training data (variety of platforms/chemistries/depths)
- STREL-154 Adjust active window size depending on neighboring sequences
- STREL-188 Haplotype generation includes soft-clipped reads
## v2.6.0
- STREL-163 Update germline EVS model training with ambiguous call handling
- STREL-122/STARKA-454 Updated somatic scoring model, integrate somatic indel EVS model
- STREL-175 restore somatic callability tract
- STREL-155 determine indel candidacy in active regions
- STREL-75 fix sample limit in pedicure
- STREL-118 relax MMDF using short haplotyping
## v2.5.0
- STREL-138 adjust indel theta based on hpol context in germline model
- STREL-124 integrate new RF scoring model for germline variants
- STREL-128 correctly call overlapping same pos/same length insertions
- STARKA-371 enforce indel normalization on vcf input
- STARKA-372 prototype RNA-Seq support
- STREL-50 add active region detection and short haplotype enumeration
- STREL-103 add new diplotype model for germline indels
- STREL-111 Updated somatic EVS model to support liquid tumor and bwamem bams
- STREL-116 fix germline MQ/MQ0, EVS feature norm, add new EVS dev features
- STREL-104 provide corrected germline MQ/MQ0, EVS feature norm and exome support
- STREL-101 update to simplified log-linear indel error rates
- STREL-97 use EVSF field to report all variant scoring features for germline model
- STREL-47 add new basecall error estimator
- STREL-54 fix realignment tree pruning
- STREL-33 add method docs for the new somatic scoring model
- STREL-49 Liquid tumor EVS cut-off retuning
- PED-42 Filtering model for SNV and Indels, De novo filter tuning, Better indel overlap handling
- STREL-18 add new indel error estimator
- STREL-46 fix ALT matches REF error in somatic SNV output
- STREL-43 add new hts file stream merger
- STREL 32 left-shift/standardize input alignment records
- STARKA-403 change somatic SNV/indel scoring model for liquid tumor support
- STREL-36 fix indel breakend support test
- STREL-29 fix edge-case handling of overlapping indel candidate/forcedGT input
- STREL-27 refactor indel candidacy to uniformly test all indel types
- STREL-24 fix off-by-one error in test leading to extra candidate indel noise
- STARKA-388 add ADF/ADR tags to gVCF output
- STARKA-393 stabilize forced long del genotyping in gVCF aggregator
- STARKA-351 reorg all strelka documentation
- PED-33 Forced outputs of SNVs for pedicure workflow
- STARKA-310 train somatic EVS model from features in VCF output
- STARKA-369 add option to keep all temp files to support workflow debug
- STARKA-350 more accurate runtime instrumentation for diploid gVCF workflow
- PED-48 bcftools compatibility in VCF, better sample naming
- STARKA-335 restore breakend filters for germline gVCF output
- PED Multi-sample calling introduced for SNVs and Indel + PL field added (tickets PED10-19)
- STARKA-317 Enable CRAM input, improve per-chrom depth estimation
- STARKA-306 fix rare chunk size boundary defect in RangeMap
- STARKA-323 refactor STARKA-293 to allow reuse of general accelerated binom
test
- STARKA-297 adjust somatic EVS model for exome data
- STARKA-293 remove count and frequency cutoffs for non-STR indels, replace with statistical test
- STARKA-305 add simple germline SNV calling site simulator
- STARKA-302 add feature verification for file-based scoring models
- STARKA-298 enable fully automated empirical scoring model training and test
cycle for SNVs
- STARKA-296 add somatic indel empirical scoring model
- STARKA-260 add additional somatic indel features
- STARKA-275 prevent phaser from causing memory exaustion for amplicon input
- STARKA-279 prevent rare vcf syntax error in gVCF when ALT is repeated
- STARKA-284 update license to GPLv3
- WGSW-765 Failed SNV phasing if any of the original SNVs are not represented
  in the most common alleles
- STARKA-241 turned on binom indel error model in starling and strelka
## v2.4.3
- STAR-66 Correct GT reported for ALT alleles at forced-output sites in
  continuous vf mode from "1/1" to "0/1" as appropriate.
- WGSW-724 Do not print PLs for homref sites
## v2.4.2
- STAR-65 Take the minimum GQ/GQX when overlapping indels
- STAR-64 Output forced indels & SNPs in continuous calling mode. Modify
  output of GQX to be per-call, not per-site. Fix flushing of compressed blocks
  at end of processed region
- WGSW-714 fix PLs not printed at nonref site
- STAR-60 Update codon phaser to correctly handle interaction between GT and
  alleles that are dropped due to low phasing support.
- WGSW-711 Make inclusion of phasing-related headers conditional
## v2.4.1
- STAR-62 fix setting vqsrModelName on --exome
## v2.4.0
- STARKA-257 remap somatic SVN VQSR scores
- STARKA-254 Add PL values for germline SNV and indel calls
- STAR-54 Do not output homref call in continuous mode if the site has
  alt alleles
- STAR-55 Make GT of block compressed areas with overlapping deletions "0"
## v2.3.14
- STAR-16 Continuous variant frequency calling (a.k.a. somatic)
- TNW-371 Set MQ to 0 at forced calls with no coverage
- STARKA-250 Remove hpol and ihpol indel filters
- STAR-46 Don't block-compress forced GT SNVs in regions of no coverage
- STARKA-249 Fix minor issue with somatic indel prior
- STARKA-248 tolerate reference allele insert/delete alignments in read input
- STARKA-237 new format for scoring and indel models
- STARKA-239 update strelka VQSR and parameters for higher recall
## v2.3.13
- STAR-34: When SNVs are not phased due to low read depth (<10), add
Unphased to INFO field
- STAR-42: Support forced-output SNVs in the forced-output VCF.
## v2.3.12
- STAR-32: Trigger phasing of SNVs if an indel is encountered within the phasing interval.
- STAR-14: Add --targeted-regions-bed parameter to Starling for flagging OffTarget positions. Also
  added --targetRegions to the python workflow wrapper.
- Apply SiteConflict filter to block-compressed areas that overlap a filtered indel. This was a
  side-effect of a restructuring of the code, but was deemed to be more correct than the previous behavior.
## v2.3.11
- Modify workflow generation to remove use of reflection
## v2.3.10
- STARKA-234: make --exome flag affect indel-error-model, scoring-model and indel-ref-error-factor arguments to starling
- WGSW-357 ensure that ploidy conflict filter is consistently applied for all
records with no coverage/unknown genotype
## v2.3.9
- STARKA-231 correct codon phasing for insertions/small-fragments within the
phasing range
## v2.3.8
- STARKA-227 Wrong Qscore value reported as GQX in overlapping indels
- STARKA-226 codon phaser does not account for read N-trimming
- WGSW-370 - force overlapping indels to use the same scoring model. If either is a complex indel (i.e. contains both an insertion and a deletion), then both uise the default model. Otherwise, both use the VQSR model, if enabled.
- WGSW-358 - Make homref calls (GT=0/0) use the default model, not VQSR. This includes clinical indels/forced output.
## v2.3.7
- Fix VCF version labels
- STARKA-221 remove non-vqsr depth filter labels in strelka somatic snv vcf
## v2.3.6
- STARKA-222 enable partial win32 build/visual studio development
- WGSW-293, WGSW-297 re-update germline VQSR cutoffs to reduce trio conflicts
- STARKA-220 improve build system versioning
- STARKA-217 adjust codon phaser to match read and basecall filtration of
non-phased variants, and align phased candidates near indels correctly
- STARKA-191 add denovo calling model for parent/child trios
## v2.3.5
- WGSW-293, WGSW-297 update germline VQSR cutoffs to decrease excessive het/hom ratios, and hopefully reduce trio conflicts
- STARKA-211 filter low-confidence het snp candidates from codon phaser input,
changing phased block composition.
- STARKA-214 fix vcf and bed concat for very high contig counts, including human ref w/ decoys
- STARKA-215 roll back STARKA-198 (due to trio conflict elevation and anomolous gVCF records)
- STARKA-210 fix very low frequency realigner assertion on contig edges
## v2.3.4
- STARKA-204 remove inconsistent filters at homref sites
- STARKA-196 fix nocompress sites in empty regions
## v2.3.3
- STARKA-202 limit total read buffer size to prevent ultra-high depth memory exhuastion
- STARKA-201 set memory requirements based on run context
- STARKA-198 SNP records with FILTER LowGQX and GQX>30 (when running VQSR) fixed
## v2.3.2
- STARKA-190 apply VQSR to haploid regions
- Fix somatic callability track tabix index generation
## v2.3.1
- STARKA-186 support CIGAR sequence match/mismatch (=/X) in input BAM
- STARKA-185 turn off somatic VQSR when exome configuration is selected for
## v2.3.0
- STARKA-181 Use indel error estimates to improve indel quality computation
- Update germline indel error settings to reflect MIB defaults
- STARKA-163 Use indel error estimates to improve indel candidate selection
- STARKA-179 Provide config option for somatic callability track
- STARKA-178 Consolidate single strelka config for all aligners
- Expand random forest tree count in somatic SNV VQSR
## v2.2.2
- STARKA-175 Update somatic SNV VQSR model to include MAPQ0 with improved
normalization
- STARKA-176 Prevent large deletions from locking indels outside of the
realignment buffer.
- STARKA-174 Prevent segfault due to libstdc++ bug in certain gcc versions
## v2.2.1
- STARKA-170 Add "--exome" config option for starling/strelka
## v2.2.0
- STARKA-166 add binomial-distribution-based allele bias predictors for Starling VQSR models
- STARKA-165 patch LOH artifact in somatic VQSR output
- STARKA-162 revise depth normalization in Starling VQSR
- STARKA-161 refine somatic VQSR settings
- STARKA-158 generalize gvcf nocompress to support regions/compression/tabix index
- STARKA-155 add ploidy specification via bedfile (supporting haploid and deleted states only_
- STARKA-151 add VQSR model trained using RefError*100 and median of chromosome means for depth normalization
- STARKA-152 add ref indel error multiplier to reduce undercalling of homozygous indels
- STARKA-148 add support for hetalt ins and del VQSR models to starling
- STARKA-149 fix Qrule/default filtering being applied even when a VQSR model is specified; also remove some special-case scoring
- STARKA-150 allow tabs as separators in scoring model (model.json) file
- Extend address sanitizer build support, improve build option checks
- STARKA-144 add initial site-noise strelka feature
- STARKA-73 add strelka noise extractor workflow
- STARKA-143 Improve runtime for high-depth centromere regions, especially for strelka
- Add starling option to always genotype indels provided in vcf
- STARKA-135 add strand bias feature
- STARKA-138 add read position features to strelka SNV output
- STARKA-136 add mapping info values to strelka SNV output
## v2.1.5
- Improve VQSR depth normalization
- Add CRAM input support (still experimental pending samtools idxstats
capability for CRAM)
## v2.1.4
- Turned ON DP filter for rule filters as default
- Made scoring-model a configurable option for starling pyflow
- Sync many blt_util changes with manta
- Changed scoring model names to be case insensitive, added assertion for unknown model
- Change minimum support gcc version to 4.7
- Runtime opt: reduce excessive syscalls from position based std::map's
## v2.1.3
- STARKA-127 fix strelka workflow config file interface
- STARKA-128 fix build with gcc-4.9.0, remove solexa q-scores
- STARKA-125 add starling workflow and demo
## v2.1.2
- STARKA-124 fix strelka overlapping indel issue
- STARKA-118 fix GSNAP alignment example: leading deletion in exon
- Rolled back elimination of min-vexp and min-mismatch-window default settings, these are again required as command-line
- Set refitted hpol indel model as default
- STARKA-113 Set parameter values according to Isis defaults
- STARKA-130 Read buffer not being cleared in certain phasing context
- STARKA-131 Codon-phasing crash when run with external indel candidates
## v2.1.1
- STARKA-111 transfer full strelka workflow v1 logic into starka
- STARKA-58 Refit indel homopolymer model
- STARKA-89 LowGQX filters for no-calls in reference/depth records around indels is not se
- STARKA-91 Investigate low-depth passing hom alt calls in VQSR model
- STARKA-112 bwamem + starling 2.1 crash
## v2.1
- STARKA-53 remove grouper legacy contig logic
- STARKA-29 short-range SNP phasing and arbitrary phasing window
- STARKA-57 HighRefRep should not be default rule-based filter for indels
## v2.0.21
- Minor VQSR fixes
- STARKA-52 gVCF block compression filters not cleared on single record
## v2.0.20
- Added VQSR for indels
- Updated VQSR model parameters
- Updated homopolymer error-model
- Modified block-compression parameters for better NextSeq compression
- Added option for providing bed-file with sites that should not be block-compressed
- STARKA-50 option to output somatic-callable bed file
## v2.0.17
- STARKA-48 Fixed formatting bug for high GQX
- Header fix for adjusted Nova filters
## v2.0.16
- Adjusted filters for Nova release
## v2.0.15
- STARKA-47 Accept GATK-style bam indices
- STARKA-43 Accept edge indel pattern produced by freeBayes/BamLeftAlign
- STARKA-45 Filter indels greater than max indel size from candidate indel vcf
## v2.0.14
- STARKA-41 Fix consensus open-break-end error
## v2.0.13
- Skip BWA-mem supplementary reads
- STARKA-37 Handle Skip-Delete-Skip pattern in tophat output
- STARKA-23 Accept an input VCF file for which each alternate allele must be genotyped
- STARKA-35 Fixed 255 q-score bug for RNAseq workflow
## v2.0.12
- STARKA-34 Properly handle insertions adjacent to introns
- Filter out all open-breakends from vcf output
## v2.0.11
- STARKA-14 Add option to provide sample name in output VCF SAMPLE column
- VQSR features
## v2.0.10
- STARKA-30 fix stability issue encounted with large indels in 2x400 reads
## v2.0.9
- STARKA-32 fix handling of another complex indel on read edge
## v2.0.8
- STARKA-27 correctly handle complex insert/delete indels for read edges
## v2.0.7
- STARKA-12 tolerate all edge insertions/deletions
- STARKA-17 Add option to output gVCF with no block compression
- STARKA-24 fix gVCF site records to correctly inherit spanning deletion filters
## v2.0.6
- Change makefile to build when cwd is not in PATH
## v2.0.5
- Associated with new parent strelka workflow release, no major changes from v2.0.4
## v2.0.4
- STARKA-21 Add command-line control for snv hpol filter. Set snv and indel hpol filters
off by default.
- STARKA-22 Reorganize build system around libraries, add unit test framework as
part of every build and seed framework with a few tests for each library
## v2.0.3
- STARKA-16 fix gVCF output so that GT does not contain allele numbers which are
not in the ALT tag
- Add win32 compat fixes from Eric Roller
- STARKA-13 Groom CIGAR alignments on input to remove zero-length and PAD
segments
- STARKA-11 Samtools upgraded to 0.1.18 to resolve issues reported for
strelka run with long-line version of hg19 reference.
## v2.0.2
- Add haplotype score option to command-line
## v2.0.1
- Added command-line controls for for "R8" indel filter and strand-bias
## v2.0.0
- First RC. No changes from v2.0a3
## v2.0a3
- Completed gVCF output to pass vcf-validator
- added AD tags to snps and indels
- added haplotypescore but left this turned off
- added command-line controls for min-gqx,max-depth-factor and other filters/blocking thresholds
- cleaned up other final details
## v2.0a2
- Bugfix: gVCF output was being written +1 past end of requested range
## v2.0a1
- initial version of starling with direct gVCF output
## v1.1.0
- Import all updates from starling/strelka maintained on the strelka
standalone 0.4.10 tag.
## v1.0.0
- initial transfer of v1 starling/strelka from the public strelka release
branch

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file
///
/// \author Chris Saunders
///

#include "blt_common/blt_arg_validate.hh"
#include "blt_util/compat_util.hh"
#include "starling_common/starling_option_parser.hh"

#include "boost/format.hpp"

#include <iostream>



const unsigned max_flank_size(1000);



void validate(boost::any& v,
              const std::vector<std::string>& values,
              std::vector<avg_window_data>*, int)
{
    if (v.empty()) {
        v = boost::any(std::vector<avg_window_data>());
    }
    std::vector<avg_window_data>* tv = boost::any_cast< std::vector<avg_window_data> >(&v);
    assert(NULL != tv);

    avg_window_data awd;

    // Extract tokens from values string vector and populate avg_window_data struct.
    if (values.size() != 2) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }

    try {
        awd.flank_size = boost::lexical_cast<unsigned>(values[0]);
    } catch(const boost::bad_lexical_cast&) {
        throw po::validation_error(po::validation_error::invalid_option_value);
    }
    awd.filename = values[1];

    tv->push_back(awd);
}



po::options_description
get_starling_shared_option_parser(starling_options& opt) {


    po::options_description geno_opt("genotyping options");
    geno_opt.add_options()
        ("snp-theta", po::value<double>(&opt.bsnp_diploid_theta)->default_value(opt.bsnp_diploid_theta),
         "Set snp theta.")
        ("indel-theta", po::value<double>(&opt.bindel_diploid_theta)->default_value(opt.bindel_diploid_theta),
         "Set indel theta");

    po::options_description gvcf_opt("gVCF options");
    gvcf_opt.add_options()
        ("gvcf-file",po::value<std::string>(&opt.gvcf.out_file),"gVCF output file, enter '-' to write to stdout. An argument must be specified to activeate gVCF mode")
        ("chrom-depth-file",po::value<std::string>(&opt.gvcf.chrom_depth_file),"If provided, the mean depth for each chromosome will be read from file, and these values will be used for high depth filtration. File should contain one line per chromosome, where each line begins with: \"chrom_name<TAB>depth\" (default: no chrom depth filtration)");

    po::options_description hap_opt("haplotype-options");
    hap_opt.add_options()
        ("hap-model",
         "Turn on haplotype-based variant calling");

    po::options_description blt_nonref_opt("nonref-model-options");
    blt_nonref_opt.add_options()
        ("nonref-test-file", po::value<std::string>(&opt.nonref_test_filename),
         "Test for non-reference alleles at any frequency, write results to specified filename")
        ("nonref-sites-file", po::value<std::string>(&opt.nonref_sites_filename),
         "Print results of non-reference allele test at every site to file")
        ("nonref-variant-rate", po::value<double>(&opt.nonref_variant_rate)->default_value(opt.nonref_variant_rate),
         "The expected non-reference variant frequency used with nonref-test")
        ("min-nonref-freq", po::value<double>(&opt.min_nonref_freq)->default_value(opt.min_nonref_freq),
         "The minimum non-reference allele frequency considered in nonref-test")
        ("nonref-site-error-rate", po::value<double>(&opt.nonref_site_error_rate)->default_value(opt.nonref_site_error_rate),
         "The expected rate of erroneous non-reference allele sites applied to the nonref model. At error sites a nonref allele is expected in the frequency range [0,decay_freq], with a probability that linearly decays to zero at decay_freq.")
        ("nonref-site-error-decay-freq",
         po::value<double>(&opt.nonref_site_error_decay_freq)->default_value(opt.nonref_site_error_decay_freq),
         "The decay_freq used for the site-error state as described above.");

    po::options_description contig_opt("contig-options");
    contig_opt.add_options()
        ("min-contig-open-end-support",po::value<unsigned>(&opt.min_contig_open_end_support)->default_value(opt.min_contig_open_end_support),
         "Filter out any open-ended contig with an unaligned breakpoint sequence length of less than N.")
        ("min-contig-edge-alignment",po::value<unsigned>(&opt.min_contig_edge_alignment)->default_value(opt.min_contig_edge_alignment),
         "Filter out any contig with an edge match segment shorter than N.")
        ("min-contig-contiguous-match",po::value<unsigned>(&opt.min_contig_contiguous_match)->default_value(opt.min_contig_contiguous_match),
         "Filter out any contig without a match segment of length at least N.");

    po::options_description realign_opt("realignment-options");
    realign_opt.add_options()
        ("max-indel-toggle-depth", po::value<int>(&opt.max_read_indel_toggle)->default_value(opt.max_read_indel_toggle),
         "Controls the realignment stringency. Lowering this value will increase the realignment speed at the expense of indel-call quality")
        ("skip-realignment",
         "Turns off read realignment. Only accepted when there are no indel calling options turned on and no input grouper contigs");

    po::options_description indel_opt("indel-options");
    indel_opt.add_options()
        ("max-candidate-indel-depth", 
         po::value<unsigned>(&opt.max_candidate_indel_depth)->default_value(opt.max_candidate_indel_depth),
"Maximum estimated read depth for an indel to reach candidacy. If any one sample exceeds this depth at the indel, the indel will not reach candidacy in all indel-synchronized samples")
        ("min-candidate-open-length",
         po::value<unsigned>(&opt.min_candidate_indel_open_length)->default_value(opt.min_candidate_indel_open_length),
         "Minimum open-ended breakpoint sequence length required to become a breakpoint candidate")
         ("candidate-indel-input-vcf",
          po::value<std::vector<std::string> >(&opt.input_candidate_indel_vcf)->multitoken(),
          "Add candidate indels from the specified vcf file. Option can be provided multiple times to combine evidence from multiple vcf files.");

    po::options_description window_opt("window-options");
    window_opt.add_options()
        ("variant-window-flank-file",po::value<std::vector<avg_window_data> >(&opt.variant_windows)->multitoken(),
         "Print out regional average basecall statistics at variant sites within a window of the variant call. Must provide arguments for window flank size and output file. Option can be specified multiple times. (example: '--variant-window-flank-file 10 window10.txt')");

    po::options_description compat_opt("compatibility-options");
    compat_opt.add_options()
        ("eland-compatability",
         "When argument is provided the input reads are checked for an optional AS field corresponding to the ELAND PE map score.");

    po::options_description input_opt("input-options");
    input_opt.add_options()
        ("max-input-depth", po::value<unsigned>(&opt.max_input_depth),
         "Maximum allowed read depth per sample (prior to realignment). Input reads which would exceed this depth are filtered out.  (default: no limit)")
        ("ignore-conflicting-read-names",
         "Do not report an error if two input reads share the same QNAME and read number");

    po::options_description other_opt("other-options");
    other_opt.add_options()
        ("report-file", po::value<std::string>(&opt.report_filename),
         "Report non-error run info and statistics to file");

    po::options_description new_opt("New options");

    new_opt.add(geno_opt).add(gvcf_opt).add(hap_opt).add(blt_nonref_opt).add(contig_opt).add(realign_opt).add(indel_opt).add(window_opt).add(compat_opt).add(input_opt).add(other_opt);

    return new_opt;
}



po::options_description
get_starling_option_parser(starling_options& opt) {

    po::options_description starling_parse_opt(get_starling_shared_option_parser(opt));

    po::options_description help_parse_opt("Help");
    help_parse_opt.add_options()
        ("help,h","print this message");  
    
    po::options_description visible("Options");
    visible.add(starling_parse_opt).add(help_parse_opt);

    return visible;
}



void
write_starling_legacy_options(std::ostream& os) {

    static const starling_options default_opt;
        
    os << 
        " -bam-file file     - Analyze reads from 'file' in sorted BAM format (required) \n" // (use \"" << STDIN_FILENAME << "\" for stdin)\n"
        " -bam-seq-name name - Analyze reads aligned to chromosome 'name' in the bam file (required)\n"
        " -samtools-reference file\n"
        "                    - Get the reference sequence from the multi-sequence fasta 'file' following samtools reference conventions (single-seq or samtools reference required)\n"
        "\n"
        " -bsnp-diploid-het-bias x\n"
        "                    - Set bias term for the heterozygous state in the bsnp model, such that\n"
        "                      hets are expected at allele ratios in the range [0.5-x,0.5+x] (default: 0)\n"
#if 0
        " -bsnp-monoploid x  - Use Bayesian monoploid genotype snp caller with theta=x\n"
        " -bsnp-nploid n x   - Use Bayesian nploid genotype snp caller with ploidy=n and prior(snp)=x\n"
        " -lsnp-alpha x      - Use likelihood ratio test snp caller with alpha=x\n"
#endif
        " -bsnp-diploid-file file\n"
        "                    - Run bayesian diploid genotype model, write results to 'file'\n"
        " -bsnp-diploid-allele-file file\n"
        "                    - Write the most probable genotype at every position to file\n"
#if 0
        " -anom-distro-table-alpha x\n"
        "                    - Test whether strands were sampled from different distributions at snp\n"
        "                      call sites. The test has a false positive rate of x over all snp calls.\n"
        "                      Implemented as a contingency table test.\n"
        " -anom-distro-lrt-alpha x\n"
        "                    - Test whether strands were sampled from different distributions at snp\n"
        "                      call sites. The test has a false positive rate of x over all snp calls.\n"
        "                      Implemented as a likelihood ratio test.\n"
        " -anom-cov-alpha x  - Detect strand coverage anomaly with alpha=x\n"
        " -filter-anom-calls - Don't write any variants at positions where an anomaly is detected\n"
#endif
        " -min-qscore n      - Don't use base if qscore<n (default: " << default_opt.min_qscore << ")\n"
        " -max-window-mismatch n m\n"
        "                    - Don't use base if mismatch count>n within a window of m flanking bases\n"
        " -min-single-align-score n\n"
        "                    - Reads with single align score<n are marked as SE-failed. By default such reads are excluded\n"
        "                      from consideration unless a paired score is present. This behavior can be modified by\n"
        "                      setting the exclude or rescue modes below. (default: " << default_opt.min_single_align_score << ")\n"
        " -single-align-score-exclude-mode\n"
        "                    - Exclude SE-failed reads even when a paired score is present and the read is not PE-failed.\n"
        " -single-align-score-rescue-mode\n"
        "                    - Include non SE-failed reads even when a paired score is present and the read is PE-failed.\n"
        " -min-paired-align-score n\n"
        "                    - Reads with paired align score<n are marked as PE-failed if a paired score is present. By\n"
        "                      default such reads are excluded from consideration, but may still be used if the single score\n"
        "                      rescue mode is enabled. (default score: " << default_opt.min_paired_align_score << ")\n"
        " -filter-unanchored - Don't use unanchored read pairs during variant calling. Unanchored read pairs have a single-read\n"
        "                      mapping score of zero in both reads of the pair\n"
        " -include-singleton - Include paired reads with unmapped mates\n"
        " -include-anomalous - Include paired reads which are not part of a 'proper pair' (anomalous orientation or insert size)\n"
        " -counts file       - Write observation counts for every position to 'file'\n"
        " -clobber           - Overwrite pre-existing output files\n"
        " -print-evidence    - Print the observed data at single site events (does not include indels)\n"
        " -print-all-site-evidence\n"
        "                    - Print the observed data for all sites (does not include indels)\n"
        " -bindel-diploid-het-bias x\n"
        "                    - Set bias term for the heterozygous state in the bindel model, such that\n"
        "                      hets are expected at allele ratios in the range [0.5-x,0.5+x] (default: 0)\n"
        " -bindel-diploid-file file\n"
        "                    - Run Bayesian diploid genotype caller, write results to 'file'\n"
        " -indel-contigs file\n"
        "                    - Contig file produced by GROUPER indel-finder (must be specified together with -indel-contig-reads)\n"
        " -indel-contig-reads file\n"
        "                    - Contig reads file produced by GROUPER indel-finder (must be specified together with -indel-contigs)\n"
        " -indel-error-rate x\n"
        "                    - If calling indels, set the indel error rate to a constant value of x (0<=x<=1).\n"
        "                      The default indel error rate is taken from an empirical function accounting for\n"
        "                      homopolymer length and indel type (i.e. insertion or deletion). This option\n"
        "                      overrides the default behavior.\n"
        " -indel-nonsite-match-prob x\n"
        "                    - The probability of a base matching the reference in an 'average' mismapped read. This\n"
        "                      value is used by the indel-caller only. (default: " << default_opt.indel_nonsite_match_prob << ")\n"
        " -report-range-begin n\n"
        "                    - Event reports and coverage begin at base n\n"
        "                      (default: 0)\n"
        " -report-range-end n\n"
        "                    - Event reports and coverage end after base n or min(n,ref_size) if reference\n"
        "                      specified.\n"
        "                      (default: ref_size)\n"
        " -report-range-reference\n"
        "                    - Event reports and coverage span the entire reference sequence.\n"
        "                      A reference sequence is required to use this flag. This sets begin=1 and\n"
        "                      end=ref_size. This flag cannot be combined with -report-range-begin/-end.\n"
        "                      (NOTE: this behaviour is now default, but the flag is still accepted.)\n"
        " -genome-size n     - Specify the total number of non-ambiguous bases in the genome to which the input reads\n"
        "                      have been aligned for use in indel calling.\n"
        " -min-candidate-indel-reads n\n"
        "                    - Unless an indel is supported by at least this many reads, it cannot become a candidate\n"
        "                      for realignment and indel calling. A read counts if it's either a contig read or a genomic\n"
        "                      read which passes the mapping score threshold. (default: " << default_opt.default_min_candidate_indel_reads << ")\n"

        " -min-candidate-indel-read-frac x\n"
        "                    - Unless at least this fraction of intersecting reads contain the indel, it cannot become\n"
        "                      a candidate for realignment and indel calling. A read counts if it's either a contig read\n"
        "                      or a genomic read which passes the mapping score threshold. Only genomic reads which pass\n"
        "                      the mapping score threshold are used for the denominator of this metric. (default: " << default_opt.min_candidate_indel_read_frac << ")\n"
        " -min-small-candidate-indel-read-frac x\n"
        "                    - For small indels (no more than " << default_opt.max_small_candidate_indel_size << " bases), an additional indel candidacy filter is applied:\n"
        "                      Unless at least this fraction of intersecting reads contain the small indel, it cannot become\n"
        "                      a candidate for realignment and indel calling. A read counts if it's either a contig read\n"
        "                      or a genomic read which passes the mapping score threshold. Only genomic reads which pass\n"
        "                      the mapping score threshold are used for the denominator of this metric. (default: " << default_opt.default_min_small_candidate_indel_read_frac << ")\n"
        " -max-candidate-indel-density x\n"
        "                    - If there are more than x candidate indels per base intersecting a read, then realignment\n"
        "                      is truncated to only allow individual indel toggles of the starting alignments for that read.\n"
        "                      x must be greater than 0  (default: " << default_opt.max_candidate_indel_density << ")\n"
        " -candidate-indel-file file\n"
        "                      write out all candidate indels before realignment and genotyping\n"
        " -write-candidate-indels-only\n"
        "                      Skip all analysis steps besides writing candidate indels (only valid when -candidate-indel-file is given)\n"
        " -realigned-read-file file\n"
        "                    - Write reads which have had their alignments altered during realignemnt to a BAM file.\n"
        " -realign-submapped-reads\n"
        "                    - When this argument is provided, even reads which fail the variant calling mapping thresholds\n"
        "                      are realigned using the same procedure as the variant calling reads.\n"
        " -snp-max-basecall-filter-fraction x\n"
        "                    - Do not call snps at sites where the fraction of filtered basecalls exceeds x. (default: " << default_opt.max_basecall_filter_fraction << ")\n"
        " -no-ambiguous-path-clip\n"
        "                    - Turn off ambiguous read trimming after realignment.\n"
        " -max-indel-size    - Sets the maximum size for indels processed for indel genotype calling and realignment.\n"
        "                      Increasing this value should lead to an approx linear increase in memory consumption.\n"
        "                      (default: " << default_opt.max_indel_size << ")\n"
        "\n"
        " -print-all-poly-gt - Print all polymorphic-site genotype probabilties in the diploid sites and snps files\n"
        " -print-used-allele-counts\n"
        "                    - Print used base counts for each allele in the diploid sites and snps files\n"
        " -used-allele-count-min-qscore x\n"
        "                    - If printing used allele counts, filter them for qscore >= x\n"
        "\n"
        " -all-warnings      - print all warnings (default: errors and low-frequency warnings only)\n"
        "\n"
        " -skip-variable-metadata\n"
        "                    - do not print commmand-line or time stamp in data file metadata\n"
        "\n"
        " -h                 - Display usage (this page)\n";
}



static
void
finalize_legacy_starling_options(const prog_info& pinfo,
                                 starling_options& opt) {


    if(! opt.is_ref_set()){
        pinfo.usage("a reference sequence must be specified");
    }

    // canonicalize the reference seqeunce path:
    if(opt.is_samtools_ref_set) {
        if(! compat_realpath(opt.samtools_ref_seq_file)){
            std::ostringstream oss;
            oss << "can't resolve samtools reference path: " << opt.samtools_ref_seq_file << "\n";
            pinfo.usage(oss.str().c_str());
        }
    } else {
        assert(0);
    }

    if(! opt.is_user_genome_size) {
        // this requirment is not what we want, but it's the only way to make things reliable for now:
        pinfo.usage("must specify genome-size");
    } else {
        if(opt.user_genome_size<1) {
            pinfo.usage("genome-size must be greater than 0");
        }
    }


    if(! opt.is_call_indels()) {
        if(opt.is_simple_indel_error) {
            pinfo.usage("--indel-error-rate has no effect when not calling indels");
        }
    }

    {
        const bool is_contigs(! opt.indel_contig_filename.empty());
        const bool is_reads(! opt.indel_contig_read_filename.empty());

        if(is_contigs ^ is_reads){ //should work when they're both bools...
            if(is_contigs) {
                pinfo.usage("Contigs specified without corresponding contig reads.");
            } else {
                pinfo.usage("Contig reads specifed without corresponding contigs.");
            }
        }

        if(is_contigs || is_reads) {
            if(opt.is_ignore_read_names) {
                pinfo.usage("Cannot use 'ignore-conflicting-read-names' when indel contigs are specified.");
            }
        }
    }

    if(opt.is_write_candidate_indels_only &&
       opt.candidate_indel_filename.empty()) {
        pinfo.usage("Cannot specify -write-candidate-indels-only without providing candidate indel filename.");
    }
}



void
finalize_starling_options(const prog_info& pinfo,
                          const po::variables_map& vm,
                          starling_options& opt) {

    // blt section:
    check_option_arg_range(pinfo,opt.nonref_variant_rate,"nonref-variant-rate",0.,1.);
    check_option_arg_range(pinfo,opt.min_nonref_freq,"min-nonref-freq",0.,1.);

    // starling section:

    // max_theta for indels is actually 2./3., but because we don't
    // allow non-reference hets, we stick with the lower value
    // used for snps:
    //
    if(opt.bindel_diploid_theta>MAX_DIPLOID_THETA){
        std::ostringstream oss;
        oss << "indel diploid heterozygosity exceeds maximum value of: " << MAX_DIPLOID_THETA;
        pinfo.usage(oss.str().c_str());
    }

    if(vm.count("hap-model")) {
        opt.is_htype_calling=true;
    }

    if(vm.count("eland-compatability")) {
        opt.is_eland_compat=true;
    }

    if(vm.count("max-input-depth")) {
        opt.is_max_input_depth=true;
    }

    if(vm.count("ignore-conflicting-read-names")) {
        opt.is_ignore_read_names=true;
    }

    if(vm.count("skip-realignment")) {
        if(opt.is_call_indels()) {
            pinfo.usage("Cannot disable realignment when indel-calling is selected.");
        }

        const bool is_contigs(! opt.indel_contig_filename.empty());
        const bool is_reads(! opt.indel_contig_read_filename.empty());

        if(is_contigs || is_reads) {
            pinfo.usage("Cannot disable realignment when reading grouper contigs.");
        }

        opt.is_skip_realignment=true;
    }

    std::sort(opt.variant_windows.begin(),opt.variant_windows.end());
    const unsigned vs(opt.variant_windows.size());
    unsigned last_fs(0);
    std::set<std::string> filenames;
    for(unsigned i(0);i<vs;++i) {
        const std::string& filename(opt.variant_windows[i].filename);
        if(filenames.count(filename)) {
            pinfo.usage((boost::format("variant-window-flank-file options contain a repeated filename: '%s'") % filename).str().c_str());
        }
        filenames.insert(filename);
        const unsigned fs(opt.variant_windows[i].flank_size);
        if(fs>max_flank_size) {
            pinfo.usage((boost::format("variant-window-flank-file flank size %u exceeds maximum flank size of %u") % fs % max_flank_size).str().c_str());
        }
        if(fs<1) {
            pinfo.usage((boost::format("variant-window-flank-file flank size %u is less than minimum flank size of 1") % fs).str().c_str());
        }
        if((i!=0) && (fs==last_fs)) {
            pinfo.usage((boost::format("Repeated flank size of %u provided to variant-window-flank-file") % fs).str().c_str());
        }
        last_fs=fs;
    }

    finalize_legacy_starling_options(pinfo,opt);
}




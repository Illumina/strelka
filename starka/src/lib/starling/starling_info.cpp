// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright 2009 Illumina, Inc.
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

#include "starling/starling_info.hh"
#include "starling_common/starling_option_parser.hh"
#include "starling_common/starling_shared.hh"

#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>



void
starling_info::
usage(const char* xmessage) const {

    std::ostream& os(log_os);

    os <<
        "\n" << name() << " - joint snp/small-indel caller\n"
        "\tversion: " << version() << "\n"
        "\n"
        "usage: " << name() << " -bam-file file [options] > event_report\n"
        "\n";

    static starling_options default_opt;
    static const po::options_description visible(get_starling_option_parser(default_opt));
    os << "\n\n[ ***** new options ***** ]\n\n";
    os << visible
       << "\n\n\n[ ***** legacy options ***** ]\n\n";
    write_starling_legacy_options(os);
    os << "\n";

    if(xmessage) {
        os << "\n"
           << "******** COMMAND-LINE ERROR:: " << xmessage << " ********\n"
           << "\n";
        exit(EXIT_FAILURE);
    }

    exit(EXIT_SUCCESS);
}


#if 0
void
starling_info::
doc() const {

    std::ostream& os(log_os);

    os <<
        "\n"
        "\n\nDOCUMENTATION PAGE OUT OF DATE\n\n\n"
        "This is the " << name() << " documentation page, for command-line\n"
        "usage enter:\n"
        "'" << name() <<" -h'\n"
        "\n"
        "To call snps, the program reads lines from an export file sorted by aligned 5'\n"
        "end on the forward strand. All input reads are assumed to be from the same\n"
        "chromosome. To call indels, a reference sequence must be specified together\n"
        "with contig and contig read files produced by the CASAVA indel finder (GROUPER).\n"
        "\n"
        "Output is a sequence of single-line event reports. Events are snp calls,\n"
        "indels, anomalies and summary statistics. Events are reported with a tag\n"
        "as the first word, followed by event specific information (detailed below).\n"
        "Additional output can be provided in separate files for position counts,\n"
        "BaCON allele calls and BaCON snp calls.\n"
        "\n"
       << name() << " has two primary snp-calling methods:\n"
        "1) Bayesian genotype caller:\n"
        "    This method is designed for samples from individuals of known ploidy.\n"
        "  It predicts the probability of each possible genotype given the observed\n"
        "  basecalls at each site. A snp is called whenever the probability of any\n"
        "  genotype exceeds that of the reference genotype.\n"
        "    For the monoploid and diploid versions of this snp caller, a heterozygosity\n"
        "  parameter is required. This parameter corresponds to the expected\n"
        "  probability that any site in an individual is heterozygous (diploid\n"
        "  case) or the expected fraction of single base differences between any two\n"
        "  homologous chromosomes in the population (monoploid case).\n"
        "    The N-ploid caller requires two parameters, one specifying the ploidy\n"
        "  of the sample and another specifying the prior probability of a snp. Note\n"
        "  that the prior probability of all non-reference genotypes is uniform,\n"
        "  such that this caller at n=2 will not produce the same results as the\n"
        "  diploid snp caller. For example, priors for '-bsnp-nploid 2 0.01', given\n"
        "  reference base 'A' are:\n"
        "      P(AC)=P(AG)=P(AT)=P(CC)=P(GG)=P(TT)=P(CG)=P(CT)=P(GT)= 0.01/9\n"
        "      P(AA)= 1-0.01\n"
        "  ...while priors for '-bsnp-diploid 0.01', given reference base 'A' are:\n"
        "      P(AC)=P(AG)=P(AT)= 0.01/3\n"
        "      P(CC)=P(GG)=P(TT)= 0.01/6\n"
        "      P(CG)=P(CT)=P(GT)= 0.0001/9\n"
        "      P(AA)= 1-0.01*3/2-0.0001/3\n\n"
        "2) Likelihood ratio test (LRT) caller:\n"
        "    This method is designed for populations or samples of unknown/unusual ploidy.\n"
        "  It tests whether the expected reference allele frequency is significantly less\n"
        "  than one given the observed basecalls at each site.\n"
        "    A single parameter (alpha) is required. This corresponds to the expected false\n"
        "  positive rate under ideal model assumptions (such as perfectly calibrated quality\n"
        "  scores and accurate read alignments).\n"
        "\n"
        "There are also two anomaly detection methods:\n"
        "1) Strand-dependent base distribution at snp calls:\n"
        "    Test whether the distribution of basecalls observed at snp call sites\n"
        "  is significantly different on forward and reverse strands.\n"
        "    A single parameter (alpha) is required. This corresponds to the expected false\n"
        "  positive rate of the anomaly test over all snp calls.\n"
        "2) Strand-dependent coverage:\n"
        "    Test whether the number of basecalls observed on the forward and reverse strands\n"
        "  of each site significantly deviates from a binomial distribution with p=0.5\n"
        "    A single parameter (alpha) is required. This corresponds to the expected false\n"
        "  positive rate of the anomaly test over all sites.\n"
        "\n"
        "Possible events:\n"
        "LSNP - snp called by the LRT snp caller\n"
        "BSNP2 - snp called by the Bayesian diploid genotype caller\n"
        "BSNP1 - snp called by the Bayesian monoploid genotype caller\n"
        "BSNPN - snp called by the Bayesian nploid genotype caller\n"
        "BINDEL2 - indel called by the Bayesian diploid genotype caller\n"
        "ANOM_COV - anomalous strand coverage at site\n"
        "ANOM_DIS - anomalous dependency of base distribution on strand\n"
        "ALLSITES_COVERAGE - coverage statistics after alignment score filtration (see below)\n"
        "ALLSITES_COVERAGE_USED - coverage statistics after all basecall filtration (see below)\n"
        "NO_REF_N_COVERAGE - ALLSITES_COVERAGE without reference 'N' positions (see below)\n"
        "NO_REF_N_COVERAGE_USED - ALLSITES_COVERAGE_USED without reference 'N' positions (see below)\n"
        "CMDLINE - echo the snp_caller command line\n"
        "\n"
        "Note that when the '-print-evidence' flag is given, a multi-line event with the\n"
        "tag: SITE_EVIDENCE will be printed after all site events at that position. Note\n"
        "that indels are not site events.\n"
        "\n"
        "Coverage reports:\n"
        "  The coverage information reported with tags ALLSITES_COVERAGE and NO_REF_N_COVERAGE\n"
        "includes all basecalls from reads that pass the alignment score filters. 'N's are included\n"
        "except where a continuous sequence of 'N's is found at the end of a read. These rules are\n"
        "also applied to the 'bcalls' column reported in the BaCON allele-caller file.\n"
        "  The coverage reported with tags ALLSITES_COVERAGE_USED and NO_REF_N_COVERAGE_USED\n"
        "includes only basecalls used for snp-calling. This basecall count is also reported for\n"
        "individual sites in the 'bcalls_used' column of several output files.\n"
        "  The range of the coverage calculation is between the first and last base covered by\n"
        "any read, unless -report-range-begin or -report-range-end are set. Note that the default\n"
        "range infered from the input reads could reduce the accuracy of the coverage estimate,\n"
        "so use of the -report-range-{begin,end} flags is recomended.\n"
        "  ALLSITES_COVERAGE* results are generated using every position in the begin-end range\n"
        "described above. NO_REF_N_COVERAGE* results are generated using every position in this\n"
        "range where the reference base is not 'N'. NO_REF_N_COVERAGE* is only reported when a\n"
        "reference sequence is specified.\n"
        "\n"
        "Event information:\n"
        "pos:           position in reference genome\n"
        "ref:           reference base\n"
        "P(snp):        probability of...\n"
        "               ...reference allele frequency being less than one (LRT model)\n"
        "               ...any non-reference genotype (Bayesian genotype models)\n"
        "Q(snp):        Qphred(reference allele posterior probability)\n"
        "P(indel):      probability of any non-reference indel genotype\n"
        "Q(snp):        Qphred(reference sequence posterior probability)\n"
        "freq(X):       maximum likelihood freqeuncy of allele X\n"
        "max_gtype:     genotype with highest posterior probability\n"
        "P(max_gtype):  highest posterior probability\n"
        "Q(max_gtype):  Qphred(1-highest posterior probability)\n"
        "max2_gtype:    genotype with second-highest posterior probability\n"
        "P(max2_gtype): second-highest posterior probability\n"
        "\n"
        "Note that P(snp) >= P(max_gtype).\n"
        "\n"
        "Output files can also be produced for BaCON allele calls, BaCON snp calls or counts\n"
        "\n"
        "BaCON allele call file:\n"
        "  TBD -- see CASAVA documentation.\n"
        "\n"
        "BaCON snp call file:\n"
        "  TBD -- see CASAVA documentation, except note that filtration of snps at positions\n"
        "with depth >= 3 times the chromosomal mean is not done.\n"
        "\n"
        "Counts file:\n"
        "The program can also write a summary file of observation counts for every position when\n"
        "'-counts filename' is specified on the command-line. The file can start with a\n"
        "a series of comment lines indicated by a leading '#' character, followed by lines\n"
        "with the following tab-delimited fields:\n"
        "1. reference sequence position number\n"
        "2. no of A basecalls used\n"
        "3. no of C basecalls used\n"
        "4. no of G basecalls used\n"
        "5. no of T basecalls used\n"
        "5. no of unused basecalls\n"
        "\n"
        "  Reads that failed alignment score filters and trailing 'N's are not included in\n"
        "the counts file. Among the counts, unambiguous bases that passed all snp calling\n"
        "filters are summarized in the 'used' counts, all others are summarized in\n"
        "the 'unused' counts.\n"
        "\n"
        "Caveats:\n"
        "- No circular genome support, negative alignment positions in reads are ignored.\n"
        "- 'N' basecalls are ignored (except in coverage tests).\n"
        "- Only coverage tests are conducted at sites with 'N' in the reference sequence.\n"
        "\n";

    exit(EXIT_SUCCESS);
}
#endif

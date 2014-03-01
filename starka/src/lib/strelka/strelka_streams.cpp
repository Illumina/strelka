// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

/// \file

/// \author Chris Saunders
///

#include "blt_util/bam_dumper.hh"
#include "blt_util/vcf_util.hh"
#include "strelka_streams.hh"

#include <cassert>

#include <fstream>
#include <iostream>


const strelka_sample_info ssi;



strelka_streams::
strelka_streams(
    const strelka_options& opt,
    const prog_info& pinfo,
    const bam_header_t* const header)
    : base_t(opt,pinfo,ssi)
{
    using namespace STRELKA_SAMPLE_TYPE;

    if (opt.is_bindel_diploid_file) {
        _bindel_diploid_osptr[NORMAL].reset(initialize_bindel_file(opt,pinfo,opt.bindel_diploid_filename,"normal-sample"));
    }
    if (opt.is_tumor_bindel_diploid()) {
        _bindel_diploid_osptr[TUMOR].reset(initialize_bindel_file(opt,pinfo,opt.tumor_bindel_diploid_filename,"tumor-sample"));
    }

    if (opt.is_realigned_read_file) {
        _realign_bam_ptr[NORMAL].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.realigned_read_filename,"normal sample realigned-read BAM",header));
    }

    if (opt.is_tumor_realigned_read()) {
        _realign_bam_ptr[TUMOR].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.tumor_realigned_read_filename,"tumor sample realigned-read BAM",header));
    }

    if (opt.is_somatic_snv()) {
        const char* const cmdline(opt.cmdline.c_str());

        std::ofstream* fosptr(new std::ofstream);
        _somatic_snv_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.somatic_snv_filename,"somatic-snv",opt.is_clobber,fos);

        write_vcf_audit(opt,pinfo,cmdline,header,fos);
        fos << "##content=strelka somatic snv calls\n"
            << "##germlineSnvTheta=" << opt.bsnp_diploid_theta << "\n"
            << "##priorSomaticSnvRate=" << opt.somatic_snv_rate << "\n";

        // INFO:
        fos << "##INFO=<ID=QSS,Number=1,Type=Integer,Description=\"Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal\">\n";
        fos << "##INFO=<ID=TQSS,Number=1,Type=Integer,Description=\"Data tier used to compute QSS\">\n";
        fos << "##INFO=<ID=NT,Number=1,Type=String,Description=\"Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.\">\n";
        fos << "##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description=\"Quality score reflecting the joint probability of a somatic variant and NT\">\n";
        fos << "##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description=\"Data tier used to compute QSS_NT\">\n";
        fos << "##INFO=<ID=SGT,Number=1,Type=String,Description=\"Most likely somatic genotype excluding normal noise states\">\n";
        fos << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";

        // FORMAT:
        fos << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tier1 (used+filtered)\">\n";
        fos << "##FORMAT=<ID=FDP,Number=1,Type=Integer,Description=\"Number of basecalls filtered from original read depth for tier1\">\n";
        fos << "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Number of reads with deletions spanning this site at tier1\">\n";
        fos << "##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description=\"Number of reads below tier1 mapping quality threshold aligned across this site\">\n";
        fos << "##FORMAT=<ID=AU,Number=2,Type=Integer,Description=\"Number of 'A' alleles used in tiers 1,2\">\n";
        fos << "##FORMAT=<ID=CU,Number=2,Type=Integer,Description=\"Number of 'C' alleles used in tiers 1,2\">\n";
        fos << "##FORMAT=<ID=GU,Number=2,Type=Integer,Description=\"Number of 'G' alleles used in tiers 1,2\">\n";
        fos << "##FORMAT=<ID=TU,Number=2,Type=Integer,Description=\"Number of 'T' alleles used in tiers 1,2\">\n";

        fos << vcf_col_label() << "\tFORMAT";
        for (unsigned s(0); s<STRELKA_SAMPLE_TYPE::SIZE; ++s) {
            fos << "\t" << STRELKA_SAMPLE_TYPE::get_label(s);
        }
        fos << "\n";
    }

    if (opt.is_somatic_indel()) {
        const char* const cmdline(opt.cmdline.c_str());

        std::ofstream* fosptr(new std::ofstream);
        _somatic_indel_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);

        open_ofstream(pinfo,opt.somatic_indel_filename,"somatic-indel",opt.is_clobber,fos);

        write_vcf_audit(opt,pinfo,cmdline,header,fos);
        fos << "##content=strelka somatic indel calls\n"
            << "##germlineIndelTheta=" << opt.bindel_diploid_theta << "\n"
            << "##priorSomaticIndelRate=" << opt.somatic_indel_rate << "\n";

        // INFO:
        fos << "##INFO=<ID=QSI,Number=1,Type=Integer,Description=\"Quality score for any somatic variant, ie. for the ALT haplotype to be present at a significantly different frequency in the tumor and normal\">\n";
        fos << "##INFO=<ID=TQSI,Number=1,Type=Integer,Description=\"Data tier used to compute QSI\">\n";
        fos << "##INFO=<ID=NT,Number=1,Type=String,Description=\"Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.\">\n";
        fos << "##INFO=<ID=QSI_NT,Number=1,Type=Integer,Description=\"Quality score reflecting the joint probability of a somatic variant and NT\">\n";
        fos << "##INFO=<ID=TQSI_NT,Number=1,Type=Integer,Description=\"Data tier used to compute QSI_NT\">\n";
        fos << "##INFO=<ID=SGT,Number=1,Type=String,Description=\"Most likely somatic genotype excluding normal noise states\">\n";
        fos << "##INFO=<ID=RU,Number=1,Type=String,Description=\"Smallest repeating sequence unit in inserted or deleted sequence\">\n";
        fos << "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Number of times RU repeats in the reference allele\">\n";
        fos << "##INFO=<ID=IC,Number=1,Type=Integer,Description=\"Number of times RU repeats in the indel allele\">\n";
        fos << "##INFO=<ID=IHP,Number=1,Type=Integer,Description=\"Largest reference interrupted homopolymer length intersecting with the indel\">\n";
        fos << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
        fos << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";
        fos << "##INFO=<ID=OVERLAP,Number=0,Type=Flag,Description=\"Somatic indel possibly overlaps a second indel.\">\n";

        // FORMAT:
        fos << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tier1\">\n";
        fos << "##FORMAT=<ID=DP2,Number=1,Type=Integer,Description=\"Read depth for tier2\">\n";
        fos << "##FORMAT=<ID=TAR,Number=2,Type=Integer,Description=\"Reads strongly supporting alternate allele for tiers 1,2\">\n";
        fos << "##FORMAT=<ID=TIR,Number=2,Type=Integer,Description=\"Reads strongly supporting indel allele for tiers 1,2\">\n";
        fos << "##FORMAT=<ID=TOR,Number=2,Type=Integer,Description=\"Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2\">\n";

        fos << vcf_col_label() << "\tFORMAT";
        for (unsigned s(0); s<STRELKA_SAMPLE_TYPE::SIZE; ++s) {
            fos << "\t" << STRELKA_SAMPLE_TYPE::get_label(s);
        }
        fos << "\n";
    }

}



// dtor here to make auto_ptr work correctly:
strelka_streams::
~strelka_streams() {}

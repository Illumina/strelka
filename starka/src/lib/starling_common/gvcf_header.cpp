// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file

/// \author Chris Saunders
///

#include "blt_util/vcf_util.hh"
#include "starling_common/gvcf_header.hh"
#include "starling_common/gvcf_locus_info.hh"

#include <cassert>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>



static
void
write_vcf_filter(std::ostream& os,
                 const char* id,
                 const char* desc) {
    os << "##FILTER=<ID=" << id << ",Description=\"" << desc << "\">\n";
}



static
void
add_gvcf_filters(const gvcf_options& opt,
                 const cdmap_t& chrom_depth,
                 std::ostream& os) {

    using namespace VCF_FILTERS;

    write_vcf_filter(os,get_label(IndelConflict),"Locus is in region with conflicting indel calls");
    write_vcf_filter(os,get_label(SiteConflict),"Site genotype conflicts with proximal indel call. This is typically a heterozygous SNV call made inside of a heterozygous deletion");


    if(opt.is_min_gqx) {
        std::ostringstream oss;
        oss << "Locus GQX is less than " << opt.min_gqx << " or not present";
        write_vcf_filter(os,get_label(LowGQX),oss.str().c_str());
    }

    if(opt.is_max_base_filt) {
        std::ostringstream oss;
        oss << "The fraction of basecalls filtered out at a site is greater than " << opt.max_base_filt;
        write_vcf_filter(os,get_label(HighBaseFilt),oss.str().c_str());
    }

    if(opt.is_max_snv_sb) {
        std::ostringstream oss;
        oss << "SNV strand bias value (SNVSB) exceeds " << opt.max_snv_sb;
        write_vcf_filter(os,get_label(HighSNVSB),oss.str().c_str());
    }
    if(opt.is_max_snv_hpol) {
        std::ostringstream oss;
        oss << "SNV contextual homopolymer length (SNVHPOL) exceeds " << opt.max_snv_hpol;
        write_vcf_filter(os,get_label(HighSNVHPOL),oss.str().c_str());
    }

    if(opt.is_max_ref_rep) {
        std::ostringstream oss;
        oss << "Locus contains an indel allele occurring in a homopolymer or dinucleotide track with a reference repeat greater than " << opt.max_ref_rep;
        write_vcf_filter(os,get_label(HighRefRep),oss.str().c_str());
    }

    if(opt.is_max_depth_factor && (! chrom_depth.empty())) {
        std::ostringstream oss;
        oss << "Locus depth is greater than " << opt.max_depth_factor << "x the mean chromosome depth";
        write_vcf_filter(os,get_label(HighDepth),oss.str().c_str());

        std::ofstream tmp_os;
        tmp_os.copyfmt(os);
        os << std::fixed << std::setprecision(2);

        cdmap_t::const_iterator i(chrom_depth.begin()), i_end(chrom_depth.end());
        for(; i!=i_end; ++i) {
            const std::string& chrom(i->first);
            const double max_depth(opt.max_depth_factor*i->second);
            os << "##MaxDepth_" << chrom << '=' << max_depth << "\n";
        }

        os.copyfmt(tmp_os);
    }
}



void
finish_gvcf_header(const gvcf_options& opt,
                   const cdmap_t& chrom_depth,
                   std::ostream& os) {

    //INFO:
    os << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the region described in this record\">\n";
    os << "##INFO=<ID=BLOCKAVG_min30p3a,Number=0,Type=Flag,Description=\"Non-variant site block. All sites in a block are constrained to be non-variant, have the same filter value, and have all sample values in range [x,y], y <= max(x+3,(x*1.3)). All printed site block sample values are the minimum observed in the region spanned by the block\">\n";

    // site specific:
    os << "##INFO=<ID=SNVSB,Number=1,Type=Float,Description=\"SNV site strand bias\">\n";
    os << "##INFO=<ID=SNVHPOL,Number=1,Type=Integer,Description=\"SNV contextual homopolymer length\">\n";

    // indel specific:
    os << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">\n";
    os << "##INFO=<ID=RU,Number=A,Type=String,Description=\"Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs are not reported if longer than 20 bases.\">\n";
    os << "##INFO=<ID=REFREP,Number=A,Type=Integer,Description=\"Number of times RU is repeated in reference.\">\n";
    os << "##INFO=<ID=IDREP,Number=A,Type=Integer,Description=\"Number of times RU is repeated in indel allele.\">\n";

    //FORMAT:
    os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    os << "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n";
    os << "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">\n";
    os << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">\n";
    os << "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Basecalls filtered from input prior to site genotyping\">\n";

    os << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">\n";

    os << "##FORMAT=<ID=DPI,Number=1,Type=Integer,Description=\"Read depth associated with indel, taken from the site preceding the indel.\">\n";

    // FILTER:
    add_gvcf_filters(opt,chrom_depth,os);

    os << vcf_col_label() << "\tFORMAT\tSAMPLE\n";
}

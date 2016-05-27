// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

///
/// \author Chris Saunders
///

#include "gvcf_header.hh"
#include "germlineVariantEmpiricalScoringFeatures.hh"
#include "gvcf_locus_info.hh"
#include "blt_util/io_util.hh"
#include "htsapi/bam_header_util.hh"
#include "htsapi/vcf_util.hh"

#include <cassert>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>



static
void
add_gvcf_filters(
    const starling_options& sopt,
    const cdmap_t& chrom_depth,
    std::ostream& os,
    const ScoringModelManager& CM)
{
    const gvcf_options& opt(sopt.gvcf);

    using namespace GERMLINE_VARIANT_VCF_FILTERS;
    os << "##VariantQualityScoringModel=" << sopt.germline_variant_scoring_model_name << "\n";
    write_vcf_filter(os,get_label(IndelConflict),"Locus is in region with conflicting indel calls");
    write_vcf_filter(os,get_label(SiteConflict),"Site genotype conflicts with proximal indel call. This is typically a heterozygous SNV call made inside of a heterozygous deletion");

    if (opt.is_min_gqx)
    {
        std::ostringstream oss;
        oss << "Locus GQX is less than " << opt.min_gqx << " or not present";
        write_vcf_filter(os,get_label(LowGQX),oss.str().c_str());
    }

    if ( opt.is_max_base_filt)
    {
        std::ostringstream oss;
        oss << "The fraction of basecalls filtered out at a site is greater than " << opt.max_base_filt;
        write_vcf_filter(os,get_label(HighBaseFilt),oss.str().c_str());
    }

    if (opt.is_max_snv_sb && (!CM.is_current_logistic()))
    {
        std::ostringstream oss;
        oss << "SNV strand bias value (SNVSB) exceeds " << opt.max_snv_sb;
        write_vcf_filter(os,get_label(HighSNVSB),oss.str().c_str());
    }
    if (opt.is_max_snv_hpol)
    {
        std::ostringstream oss;
        oss << "SNV contextual homopolymer length (SNVHPOL) exceeds " << opt.max_snv_hpol;
        write_vcf_filter(os,get_label(HighSNVHPOL),oss.str().c_str());
    }

    if (opt.is_max_ref_rep && !CM.is_current_logistic())
    {
        std::ostringstream oss;
        oss << "Locus contains an indel allele occurring in a homopolymer or dinucleotide track with a reference repeat greater than " << opt.max_ref_rep;
        write_vcf_filter(os,get_label(HighRefRep),oss.str().c_str());
    }

    if (CM.is_current_logistic())
    {
        for (unsigned varType(0); varType<CALIBRATION_MODEL::SIZE; ++varType)
        {
            using namespace CALIBRATION_MODEL;
            if (varType == HetAltSNP) continue;

            const var_case vti(static_cast<var_case>(varType));
            std::ostringstream oss;
            oss << "Locus GQX is less than " << CM.get_case_cutoff(vti)  << " for " << get_label_header(vti);
            write_vcf_filter(os,GERMLINE_VARIANT_VCF_FILTERS::get_label(get_Qscore_filter(vti)),oss.str().c_str());
        }
    }

    // Inconsistent phasing, meaning we cannot confidently identify haplotypes in windows
    if (sopt.do_codon_phasing ||
        (opt.include_headers.end() != std::find(opt.include_headers.begin(), opt.include_headers.end(), "Phasing")))
    {
        std::ostringstream oss;
        oss << "Locus read evidence displays unbalanced phasing patterns";
        write_vcf_filter(os,get_label(PhasingConflict),oss.str().c_str());
    }

    if ((! chrom_depth.empty()))
    {
        std::ostringstream oss;
        oss << "Locus depth is greater than " << opt.max_depth_factor << "x the mean chromosome depth";
        write_vcf_filter(os,get_label(HighDepth),oss.str().c_str());

        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        for (const auto& val : chrom_depth)
        {
            const std::string& chrom(val.first);
            const double max_depth(opt.max_depth_factor*val.second);
            os << "##MaxDepth_" << chrom << '=' << max_depth << "\n";
        }
    }

    // even if no ploidy bed file is provided, this filter should still exist in principal, so I don't
    // see any reason to leave it in the header for all cases:
    if (true)
    {
        write_vcf_filter(os,get_label(PloidyConflict),"Genotype call from variant caller not consistent with chromosome ploidy");
    }
    if (!sopt.gvcf.targeted_regions_bedfile.empty())
    {
        write_vcf_filter(os,get_label(OffTarget),"Variant was found in a non-targeted region");
    }
}



void
finish_gvcf_header(
        const starling_options& opt,
        const gvcf_deriv_options& dopt,
        const cdmap_t& chrom_depth,
        const std::string& sample_name,
        std::ostream& os,
        const ScoringModelManager& CM)
{
    //INFO:
    os << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the region described in this record\">\n";

    {
        const StreamScoper ss(os);// scope precision changes for output to not affect high DP/GQX values in variant records
        os << "##INFO=<ID=" << dopt.block_label << ",Number=0,Type=Flag,Description=\"Non-variant site block. All sites in a block are constrained to be non-variant, have the same filter value, and have all sample values in range [x,y], y <= max(x+" << opt.gvcf.block_abs_tol << ",(x*" << std::setprecision(2) << static_cast<double>(100+opt.gvcf.block_percent_tol)/100. << ")). All printed site block sample values are the minimum observed in the region spanned by the block\">\n";
    }

    os << "##INFO=<ID=SNVSB,Number=1,Type=Float,Description=\"SNV site strand bias\">\n";
    os << "##INFO=<ID=SNVHPOL,Number=1,Type=Integer,Description=\"SNV contextual homopolymer length\">\n";

    // indel specific:
    os << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">\n";
    os << "##INFO=<ID=RU,Number=A,Type=String,Description=\"Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs are not reported if longer than 20 bases.\">\n";
    os << "##INFO=<ID=REFREP,Number=A,Type=Integer,Description=\"Number of times RU is repeated in reference.\">\n";
    os << "##INFO=<ID=IDREP,Number=A,Type=Integer,Description=\"Number of times RU is repeated in indel allele.\">\n";

#ifdef SUPPORT_LEGACY_EVS_TRAINING_SCRIPTS
    if (opt.isReportEVSFeatures)
    {
        os << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS of mapping quality\">\n";
        os << "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">\n";
        os << "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref mapping qualities\">\n";
        os << "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base-call qualities\">\n";
        os << "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position\">\n";
        /* not currently used */
        //os << "##INFO=<ID=MapQ0Count,Number=1,Type=Integer,Description=\"Number of overlapping reads with MAPQ=0\">\n";
        os << "##INFO=<ID=AvgBaseQ,Number=1,Type=Float,Description=\"Mean base Qscore\">\n";
        os << "##INFO=<ID=AvgPos,Number=1,Type=Float,Description=\"Mean position in aligned reads\">\n";
    }
#endif

    // Unphased, flag if a site that is within a phasing window hasn't been phased
    if (opt.do_codon_phasing ||
        opt.gvcf.include_headers.end() != std::find(opt.gvcf.include_headers.begin(), opt.gvcf.include_headers.end(), "Phasing"))
    {
        os << "##INFO=<ID=Unphased,Number=0,Type=Flag,Description=\"Indicates a record that is within the specified phasing window of another variant but could not be phased due to lack of minimum read support.\">\n";
    }

    if (opt.isReportEVSFeatures)
    {
        os << "##INFO=<ID=EVSF,Number=.,Type=Float,Description=\"Empirical variant scoring features.\">\n";
    }

    //FORMAT:
    os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    os << "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n";
    os << "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Empirically calibrated variant quality score for variant sites, otherwise Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">\n";
//    os << "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">\n";
    os << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">\n";
    os << "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Basecalls filtered from input prior to site genotyping\">\n";

    os << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">\n";
    os << "##FORMAT=<ID=ADF,Number=.,Type=Integer,Description=\"Allelic depths on the forward strand\">\n";
    os << "##FORMAT=<ID=ADR,Number=.,Type=Integer,Description=\"Allelic depths on the reverse strand\">\n";
    if (!opt.is_ploidy_prior ||
        opt.gvcf.include_headers.end() != std::find(opt.gvcf.include_headers.begin(), opt.gvcf.include_headers.end(), "VF"))
    {
        os << "##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant frequency\">\n";
    }


    os << "##FORMAT=<ID=DPI,Number=1,Type=Integer,Description=\"Read depth associated with indel, taken from the site preceding the indel.\">\n";
    os << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification.\">\n";

    // FILTER:

    add_gvcf_filters(opt,chrom_depth,os,CM);


    if (opt.isReportEVSFeatures)
    {
        os << "##snv_scoring_features=";
        for (unsigned featureIndex = 0; featureIndex < GERMLINE_SNV_SCORING_FEATURES::SIZE; ++featureIndex)
        {
            if (featureIndex > 0)
            {
                os << ",";
            }
            os << GERMLINE_SNV_SCORING_FEATURES::get_feature_label(featureIndex);
        }
        for (unsigned featureIndex = 0; featureIndex < GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::SIZE; ++featureIndex)
        {
            os << ',' << GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::get_feature_label(featureIndex);
        }
        os << "\n";

        os << "##indel_scoring_features=";
        for (unsigned featureIndex = 0; featureIndex < GERMLINE_INDEL_SCORING_FEATURES::SIZE; ++featureIndex)
        {
            if (featureIndex > 0)
            {
                os << ",";
            }
            os << GERMLINE_INDEL_SCORING_FEATURES::get_feature_label(featureIndex);
        }
        for (unsigned featureIndex = 0; featureIndex < GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::SIZE; ++featureIndex)
        {
            os << ',' << GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::get_feature_label(featureIndex);
        }
        os << "\n";
    }

    os << vcf_col_label() << "\tFORMAT\t" << sample_name << "\n";
}

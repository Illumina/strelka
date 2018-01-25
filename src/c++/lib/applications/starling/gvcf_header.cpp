//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#include "gvcf_header.hh"
#include "gvcf_locus_info.hh"
#include "blt_util/io_util.hh"
#include "htsapi/bam_header_util.hh"
#include "htsapi/vcf_util.hh"

#include <cassert>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>


/// \param isGenomeVCF If true, the header output is formatted for gVCF. Otherwise the header output is formatted
///                    for a conventional variants VCF file.
static
void
addFiltersToGermlineVCFHeader(
    const starling_options& sopt,
    const cdmap_t& chrom_depth,
    const bool isGenomeVCF,
    std::ostream& os)
{
    const gvcf_options& opt(sopt.gvcf);

    using namespace GERMLINE_VARIANT_VCF_FILTERS;
    write_vcf_filter(os,get_label(IndelConflict),"Indel genotypes from two or more loci conflict in at least one sample");
    write_vcf_filter(os,get_label(SiteConflict),"Site is filtered due to an overlapping indel call filter");

    if (opt.is_min_gqx)
    {
        std::ostringstream oss;
        oss << "Locus GQX is below threshold or not present";
        write_vcf_filter(os,get_label(LowGQX),oss.str().c_str());
    }

    if ( opt.is_max_base_filt)
    {
        std::ostringstream oss;
        oss << "The fraction of basecalls filtered out at a site is greater than " << opt.max_base_filt;
        write_vcf_filter(os,get_label(HighBaseFilt),oss.str().c_str());
    }

    if (opt.is_max_snv_sb)
    {
        std::ostringstream oss;
        oss << "Sample SNV strand bias value (SB) exceeds " << opt.max_snv_sb;
        write_vcf_filter(os,get_label(HighSNVSB),oss.str().c_str());
    }
    if (opt.is_max_snv_hpol)
    {
        std::ostringstream oss;
        oss << "SNV contextual homopolymer length (SNVHPOL) exceeds " << opt.max_snv_hpol;
        write_vcf_filter(os,get_label(HighSNVHPOL),oss.str().c_str());
    }

    if (opt.is_max_ref_rep())
    {
        std::ostringstream oss;
        oss << "Locus contains an indel allele occurring in a homopolymer or dinucleotide track with a reference repeat greater than " << opt.max_ref_rep;
        write_vcf_filter(os,get_label(HighRefRep),oss.str().c_str());
    }

    if (! chrom_depth.empty())
    {
        std::ostringstream oss;
        oss << "Locus depth is greater than " << opt.max_depth_factor << "x the mean chromosome depth";
        write_vcf_filter(os,get_label(HighDepth),oss.str().c_str());

        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        for (const auto& val : chrom_depth)
        {
            const std::string& chrom(val.first);
            const double expectedDepth(val.second);
            os << "##Depth_" << chrom << '=' << expectedDepth << "\n";
        }
    }

    // Add low depth filter
    {
        std::ostringstream oss;
        oss << "Locus depth is below " << opt.minPassedCallDepth;
        write_vcf_filter(os,get_label(LowDepth),oss.str().c_str());
    }

    // Add NotGenotyped filter
    {
        std::ostringstream oss;
        oss << "Locus contains forcedGT input alleles which could not be genotyped";
        write_vcf_filter(os,get_label(NotGenotyped),oss.str().c_str());
    }

    // even if no ploidy bed file is provided, this filter should still exist, so I don't
    // see any reason to leave it in the header for all cases:
    write_vcf_filter(os,get_label(PloidyConflict),"Genotype call from variant caller not consistent with chromosome ploidy");

    if (! isGenomeVCF)
    {
        write_vcf_filter(os, get_label(NoPassedVariantGTs), "No samples at this locus pass all sample filters and have a variant genotype");
    }
}



void
finishGermlineVCFheader(
    const starling_options& opt,
    const gvcf_deriv_options& dopt,
    const cdmap_t& chrom_depth,
    const std::vector<std::string>& sampleNames,
    const bool isGenomeVCF,
    std::ostream& os)
{
    //INFO:
    os << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the region described in this record\">\n";

    {
        const StreamScoper ss(os);// scope precision changes for output to not affect high DP/GQX values in variant records
        os << "##INFO=<ID=" << dopt.block_label << ",Number=0,Type=Flag,Description=\"Non-variant multi-site block. Non-variant blocks are defined independently for each sample. All sites in such a block are constrained to be non-variant, have the same filter value, and have sample values {GQX,DP,DPF} in range [x,y], y <= max(x+" << opt.gvcf.block_abs_tol << ",(x*" << std::setprecision(2) << static_cast<double>(100+opt.gvcf.block_percent_tol)/100. << ")).\">\n";
    }

    os << "##INFO=<ID=SNVHPOL,Number=1,Type=Integer,Description=\"SNV contextual homopolymer length\">\n";

    // indel specific:
    os << "##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">\n";
    os << "##INFO=<ID=RU,Number=A,Type=String,Description=\"Smallest repeating sequence unit extended or contracted in the indel allele relative to the reference. RUs are not reported if longer than 20 bases\">\n";
    os << "##INFO=<ID=REFREP,Number=A,Type=Integer,Description=\"Number of times RU is repeated in reference\">\n";
    os << "##INFO=<ID=IDREP,Number=A,Type=Integer,Description=\"Number of times RU is repeated in indel allele\">\n";
    os << "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"RMS of mapping quality\">\n";

#if 0
    os << "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">\n";
    os << "##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref mapping qualities\">\n";
    os << "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base-call qualities\">\n";
    os << "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position\">\n";
    os << "##INFO=<ID=AvgBaseQ,Number=1,Type=Float,Description=\"Mean base Qscore\">\n";
    os << "##INFO=<ID=AvgPos,Number=1,Type=Float,Description=\"Mean position in aligned reads\">\n";
#endif

    if (opt.isReportEVSFeatures)
    {
        os << "##INFO=<ID=EVSF,Number=.,Type=Float,Description=\"Empirical variant scoring features\">\n";
    }

    //FORMAT:
    os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    os << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    os << "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Empirically calibrated genotype quality score for variant sites, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">\n";
    os << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping. In a non-variant multi-site block this value represents the average of all sites in the block.\">\n";
    os << "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Basecalls filtered from input prior to site genotyping. In a non-variant multi-site block this value represents the average of all sites in the block.\">\n";
    os << "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum filtered basecall depth used for site genotyping within a non-variant multi-site block\">\n";

    os << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob " << opt.readConfidentSupportThreshold.strval() << " or higher that read contains indicated allele vs all other intersecting indel alleles)\">\n";
    os << "##FORMAT=<ID=ADF,Number=.,Type=Integer,Description=\"Allelic depths on the forward strand\">\n";
    os << "##FORMAT=<ID=ADR,Number=.,Type=Integer,Description=\"Allelic depths on the reverse strand\">\n";
    os << "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample filter, 'PASS' indicates that all filters have passed for this sample\">\n";
    if (!opt.is_ploidy_prior ||
        opt.gvcf.include_headers.end() != std::find(opt.gvcf.include_headers.begin(), opt.gvcf.include_headers.end(), "VF"))
    {
        os << "##FORMAT=<ID=VF,Number=1,Type=Float,Description=\"Variant frequency\">\n";
    }


    os << "##FORMAT=<ID=DPI,Number=1,Type=Integer,Description=\"Read depth associated with indel, taken from the site preceding the indel\">\n";
    os << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n";

    if (opt.isUseVariantPhaser)
    {
        os << "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set identifier\">\n";
    }

    os << "##FORMAT=<ID=SB,Number=1,Type=Float,Description=\"Sample site strand bias\">\n";

    // FILTER:
    addFiltersToGermlineVCFHeader(opt, chrom_depth, isGenomeVCF, os);


    if (opt.isReportEVSFeatures)
    {
        os << "##snv_scoring_features=";
        writeExtendedFeatureSet(dopt.snvFeatureSet,
                                dopt.snvDevelopmentFeatureSet,
                                "SNV", os);
        os << "\n";

        os << "##indel_scoring_features=";
        writeExtendedFeatureSet(dopt.indelFeatureSet,
                                dopt.indelDevelopmentFeatureSet,
                                "indel", os);
        os << "\n";
    }

    // Write the column header line:
    os << vcf_col_label() << "\tFORMAT";
    for (const auto& sampleName : sampleNames)
    {
        os << '\t' << sampleName;
    }
    os << '\n';
}

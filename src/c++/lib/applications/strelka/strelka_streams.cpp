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

///
/// \author Chris Saunders
///

#include "strelka_vcf_locus_info.hh"
#include "strelka_streams.hh"

#include "blt_util/chrom_depth_map.hh"
#include "blt_util/io_util.hh"
#include "common/Exceptions.hh"
#include "htsapi/vcf_util.hh"
#include "htsapi/bam_dumper.hh"

#include <cassert>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>


//#define DEBUG_HEADER

#ifdef DEBUG_HEADER
#include "blt_util/log.hh"
#endif



/// add vcf filter tags shared by all vcf types:
static
void
write_shared_vcf_header_info(
    const somatic_filter_options& opt,
    const somatic_filter_deriv_options& dopt,
    const bool isPrintRuleFilters,
    std::ostream& os)
{
    if (dopt.is_max_depth())
    {
        using namespace SOMATIC_VARIANT_VCF_FILTERS;

        if (isPrintRuleFilters)
        {
            std::ostringstream oss;
            oss << "Locus depth is greater than " << opt.max_depth_factor << "x the mean chromosome depth in the normal sample";
            write_vcf_filter(os,get_label(HighDepth),oss.str().c_str());
        }

        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        for (const auto& val : dopt.chrom_depth)
        {
            const std::string& chrom(val.first);
            const double expectedDepth(val.second);
            os << "##Depth_" << chrom << '=' << expectedDepth << "\n";
        }
    }
}



static
void
writeLowEVSFilter(
    std::ofstream& fos,
    const strelka_options& opt,
    const char* label)
{
    std::ostringstream oss;
    oss << "Somatic Empirical Variant Score (" << opt.SomaticEVSVcfInfoTag << ") is below threshold";
    write_vcf_filter(fos, label, oss.str().c_str());
}



strelka_streams::
strelka_streams(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const prog_info& pinfo,
    const bam_hdr_t& header,
    const StrelkaSampleSetSummary& ssi)
    : base_t(ssi.size())
{
    {
        using namespace STRELKA_SAMPLE_TYPE;
        if (opt.isWriteRealignedReads())
        {
            auto getBamPath = [&](const std::string& label)
            {
                std::ostringstream rfile;
                rfile << opt.realignedReadFilenamePrefix << label << ".bam";
                return rfile.str();
            };

            _realign_bam_ptr[NORMAL] = initialize_realign_bam(getBamPath("normal"),header);
            _realign_bam_ptr[TUMOR] = initialize_realign_bam(getBamPath("tumor"),header);
        }
    }

    if (opt.is_somatic_snv())
    {
        const char* const cmdline(opt.cmdline.c_str());

        std::ofstream* fosptr(new std::ofstream);
        _somatic_snv_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.somatic_snv_filename,"somatic-snv",fos);

        if (! opt.sfilter.is_skip_header)
        {
            write_vcf_audit(opt,pinfo,cmdline,header,fos);
            fos << "##content=strelka somatic snv calls\n"
                << "##priorSomaticSnvRate=" << opt.somatic_snv_rate << "\n";

            // this is already captured in commandline call to strelka written to the vcf header:
            //scoring_models::Instance().writeVcfHeader(fos);

            // INFO:
            fos << "##INFO=<ID=QSS,Number=1,Type=Integer,Description=\"Quality score for any somatic snv, ie. for the ALT allele to be present at a significantly different frequency in the tumor and normal\">\n";
            fos << "##INFO=<ID=TQSS,Number=1,Type=Integer,Description=\"Data tier used to compute QSS\">\n";
            fos << "##INFO=<ID=NT,Number=1,Type=String,Description=\"Genotype of the normal in all data tiers, as used to classify somatic variants. One of {ref,het,hom,conflict}.\">\n";
            fos << "##INFO=<ID=QSS_NT,Number=1,Type=Integer,Description=\"Quality score reflecting the joint probability of a somatic variant and NT\">\n";
            fos << "##INFO=<ID=TQSS_NT,Number=1,Type=Integer,Description=\"Data tier used to compute QSS_NT\">\n";
            fos << "##INFO=<ID=SGT,Number=1,Type=String,Description=\"Most likely somatic genotype excluding normal noise states\">\n";
            fos << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";
            fos << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Combined depth across samples\">\n";
            fos << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n";
            fos << "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">\n";
//            fos << "##INFO=<ID=ALTPOS,Number=1,Type=Integer,Description=\"Tumor alternate allele read position median\">\n";
//            fos << "##INFO=<ID=ALTMAP,Number=1,Type=Integer,Description=\"Tumor alternate allele read position MAP\">\n";
            fos << "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position in the tumor\">\n";
            fos << "##INFO=<ID=SNVSB,Number=1,Type=Float,Description=\"Somatic SNV site strand bias\">\n";
            fos << "##INFO=<ID=PNOISE,Number=1,Type=Float,Description=\"Fraction of panel containing non-reference noise at this site\">\n";
            fos << "##INFO=<ID=PNOISE2,Number=1,Type=Float,Description=\"Fraction of panel containing more than one non-reference noise obs at this site\">\n";

            const bool isUseEVS(dopt.somaticSnvScoringModel);
            if (isUseEVS)
            {
                fos << "##INFO=<ID=" << opt.SomaticEVSVcfInfoTag
                    << ",Number=1,Type=Float,Description=\"Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.\">\n";
            }

            if (opt.isReportEVSFeatures)
            {
                fos << "##INFO=<ID=EVSF,Number=.,Type=Float,Description=\"Empirical variant scoring features.\">\n";
            }

            // FORMAT:
            fos << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tier1 (used+filtered)\">\n";
            fos << "##FORMAT=<ID=FDP,Number=1,Type=Integer,Description=\"Number of basecalls filtered from original read depth for tier1\">\n";
            fos << "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Number of reads with deletions spanning this site at tier1\">\n";
            fos << "##FORMAT=<ID=SUBDP,Number=1,Type=Integer,Description=\"Number of reads below tier1 mapping quality threshold aligned across this site\">\n";
            fos << "##FORMAT=<ID=AU,Number=2,Type=Integer,Description=\"Number of 'A' alleles used in tiers 1,2\">\n";
            fos << "##FORMAT=<ID=CU,Number=2,Type=Integer,Description=\"Number of 'C' alleles used in tiers 1,2\">\n";
            fos << "##FORMAT=<ID=GU,Number=2,Type=Integer,Description=\"Number of 'G' alleles used in tiers 1,2\">\n";
            fos << "##FORMAT=<ID=TU,Number=2,Type=Integer,Description=\"Number of 'T' alleles used in tiers 1,2\">\n";

            // FILTERS:
            {
                using namespace SOMATIC_VARIANT_VCF_FILTERS;
                if (isUseEVS)
                {
                    assert(dopt.somaticSnvScoringModel);
                    writeLowEVSFilter(fos, opt, get_label(LowEVSsnv));
                }
                else
                {
                    {
                        std::ostringstream oss;
                        oss << "Fraction of basecalls filtered at this site in either sample is at or above " << opt.sfilter.snv_max_filtered_basecall_frac;
                        write_vcf_filter(fos, get_label(BCNoise), oss.str().c_str());
                    }
                    {
                        std::ostringstream oss;
                        oss << "Fraction of reads crossing site with spanning deletions in either sample exceeds " << opt.sfilter.snv_max_spanning_deletion_frac;
                        write_vcf_filter(fos, get_label(SpanDel), oss.str().c_str());
                    }
                    {
                        std::ostringstream oss;
                        oss << "Normal sample is not homozygous ref or ssnv Q-score < " << opt.sfilter.snv_min_qss_ref << ", ie calls with NT!=ref or QSS_NT < " << opt.sfilter.snv_min_qss_ref;
                        write_vcf_filter(fos, get_label(QSS_ref), oss.str().c_str());
                    }
                }
                {
                    std::ostringstream oss;
                    oss << "Tumor or normal sample read depth at this locus is below " << opt.sfilter.minPassedCallDepth;
                    write_vcf_filter(fos, get_label(LowDepth), oss.str().c_str());
                }
            }

            write_shared_vcf_header_info(opt.sfilter, dopt.sfilter, (! isUseEVS), fos);

            if (opt.isReportEVSFeatures)
            {
                fos << "##snv_scoring_features=";
                writeExtendedFeatureSet(SOMATIC_SNV_SCORING_FEATURES::getInstance(),
                                        SOMATIC_SNV_SCORING_DEVELOPMENT_FEATURES::getInstance(),
                                        "SNV", fos);
                fos << "\n";
            }

            fos << vcf_col_label() << "\tFORMAT";
            for (unsigned s(0); s<STRELKA_SAMPLE_TYPE::SIZE; ++s)
            {
                fos << "\t" << STRELKA_SAMPLE_TYPE::get_label(s);
            }
            fos << "\n";
        }
    }

    if (opt.is_somatic_indel())
    {
        const char* const cmdline(opt.cmdline.c_str());

        std::ofstream* fosptr(new std::ofstream);
        _somatic_indel_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);

        open_ofstream(pinfo,opt.somatic_indel_filename,"somatic-indel",fos);

        if (! opt.sfilter.is_skip_header)
        {
            write_vcf_audit(opt,pinfo,cmdline,header,fos);
            fos << "##content=strelka somatic indel calls\n"
                << "##priorSomaticIndelRate=" << opt.somatic_indel_rate << "\n";

            // this is already captured in commandline call to strelka written to the vcf header:
            //scoring_models::Instance().writeVcfHeader(fos);

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
            fos << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n";
            fos << "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">\n";
            fos << "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n";
            fos << "##INFO=<ID=OVERLAP,Number=0,Type=Flag,Description=\"Somatic indel possibly overlaps a second indel.\">\n";

            const bool isUseEVS(opt.isUseSomaticIndelScoring());

            if (isUseEVS)
            {
                fos << "##INFO=<ID=" << opt.SomaticEVSVcfInfoTag
                    << ",Number=1,Type=Float,Description=\"Somatic Empirical Variant Score (EVS) expressing the phred-scaled probability of the call being a false positive observation.\">\n";
            }

            if (opt.isReportEVSFeatures)
            {
                fos << "##INFO=<ID=EVSF,Number=.,Type=Float,Description=\"Empirical variant scoring features.\">\n";
            }

            // FORMAT:
            fos << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth for tier1\">\n";
            fos << "##FORMAT=<ID=DP2,Number=1,Type=Integer,Description=\"Read depth for tier2\">\n";
            fos << "##FORMAT=<ID=TAR,Number=2,Type=Integer,Description=\"Reads strongly supporting alternate allele for tiers 1,2\">\n";
            fos << "##FORMAT=<ID=TIR,Number=2,Type=Integer,Description=\"Reads strongly supporting indel allele for tiers 1,2\">\n";
            fos << "##FORMAT=<ID=TOR,Number=2,Type=Integer,Description=\"Other reads (weak support or insufficient indel breakpoint overlap) for tiers 1,2\">\n";

            fos << "##FORMAT=<ID=DP" << opt.sfilter.indelRegionFlankSize << ",Number=1,Type=Float,Description=\"Average tier1 read depth within " << opt.sfilter.indelRegionFlankSize << " bases\">\n";
            fos << "##FORMAT=<ID=FDP" << opt.sfilter.indelRegionFlankSize << ",Number=1,Type=Float,Description=\"Average tier1 number of basecalls filtered from original read depth within " << opt.sfilter.indelRegionFlankSize << " bases\">\n";
            fos << "##FORMAT=<ID=SUBDP" << opt.sfilter.indelRegionFlankSize << ",Number=1,Type=Float,Description=\"Average number of reads below tier1 mapping quality threshold aligned across sites within " << opt.sfilter.indelRegionFlankSize << " bases\">\n";

#if 0
            fos << "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Estimated Indel AF in tier1\">\n";
            fos << "##FORMAT=<ID=OF,Number=1,Type=Float,Description=\"Estimated frequency of supported alleles different from ALT in tier1\">\n";
            fos << "##FORMAT=<ID=SOR,Number=1,Type=Float,Description=\"Strand odds ratio, capped at [+/-]2 for tier1\">\n";
            fos << "##FORMAT=<ID=FS,Number=1,Type=Float,Description=\"Log p-value using Fisher's exact test to detect strand bias, based on tier1\">\n";
            fos << "##FORMAT=<ID=BSA,Number=1,Type=Float,Description=\"Binomial test log-pvalue for ALT allele in tier1\">\n";
            fos << "##FORMAT=<ID=RR,Number=1,Type=Float,Description=\"Read position ranksum for ALT allele in tier1 reads (U-statistic)\">\n";
#endif
            fos << "##FORMAT=<ID=BCN" << opt.sfilter.indelRegionFlankSize <<  ",Number=1,Type=Float,Description=\"Fraction of filtered reads within " << opt.sfilter.indelRegionFlankSize << " bases of the indel.\">\n";

            // FILTERS:
            {
                using namespace SOMATIC_VARIANT_VCF_FILTERS;
                if (isUseEVS)
                {
                    assert(dopt.somaticIndelScoringModel);
                    writeLowEVSFilter(fos, opt, get_label(LowEVSindel));
                }
                else
                {
                    {
                        std::ostringstream oss;
                        oss << "Average fraction of filtered basecalls within " << opt.sfilter.indelRegionFlankSize << " bases of the indel exceeds " << opt.sfilter.indelMaxWindowFilteredBasecallFrac;
                        write_vcf_filter(fos, get_label(IndelBCNoise), oss.str().c_str());
                    }
                    {
                        std::ostringstream oss;
                        oss << "Normal sample is not homozygous ref or sindel Q-score < " << opt.sfilter.sindelQuality_LowerBound << ", ie calls with NT!=ref or QSI_NT < " << opt.sfilter.sindelQuality_LowerBound;
                        write_vcf_filter(fos, get_label(QSI_ref), oss.str().c_str());
                    }
                }
                {
                    std::ostringstream oss;
                    oss << "Tumor or normal sample read depth at this locus is below " << opt.sfilter.minPassedCallDepth;
                    write_vcf_filter(fos, get_label(LowDepth), oss.str().c_str());
                }
            }

            // for indels only, we keep using the highdepth filter while EVS is on, so we need
            // to add this into the header too:
            const bool isPrintRuleFilters(true);
            write_shared_vcf_header_info(opt.sfilter, dopt.sfilter, isPrintRuleFilters, fos);

            if (opt.isReportEVSFeatures)
            {
                fos << "##indel_scoring_features=";
                writeExtendedFeatureSet(SOMATIC_INDEL_SCORING_FEATURES::getInstance(),
                                        SOMATIC_INDEL_SCORING_DEVELOPMENT_FEATURES::getInstance(),
                                        "indel", fos);
                fos << "\n";
            }

            fos << vcf_col_label() << "\tFORMAT";
            for (unsigned s(0); s<STRELKA_SAMPLE_TYPE::SIZE; ++s)
            {
                fos << "\t" << STRELKA_SAMPLE_TYPE::get_label(s);
            }
            fos << "\n";
        }
    }

    if (opt.is_somatic_callable())
    {
        std::ofstream* fosptr(new std::ofstream);
        _somatic_callable_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);

        open_ofstream(pinfo,opt.somatic_callable_filename,"somatic-callable-regions",fos);

        // post samtools 1.0 tabix doesn't handle header information anymore, so take this out entirely:
#if 0
        if (! opt.sfilter.is_skip_header)
        {
            fos << "track name=\"StrelkaCallableSites\"\t"
                << "description=\"Sites with sufficient information to call somatic alleles at 10% frequency or greater.\"\n";
        }
#endif
    }
}

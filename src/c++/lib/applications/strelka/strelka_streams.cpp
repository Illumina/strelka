// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Chris Saunders
///

#include "strelka_vcf_locus_info.hh"
#include "strelka_streams.hh"

#include "blt_util/chrom_depth_map.hh"
#include "blt_util/io_util.hh"
#include "htsapi/vcf_util.hh"
#include "htsapi/bam_dumper.hh"

#include <cassert>

#include <fstream>
#include <iomanip>
#include <iostream>
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
    std::ostream& os)
{
    if (dopt.is_max_depth())
    {
        using namespace STRELKA_VCF_FILTERS;

        std::ostringstream oss;
        oss << "Locus depth is greater than " << opt.max_depth_factor << "x the mean chromosome depth in the normal sample";
        //oss << "Greater than " << opt.max_depth_factor << "x chromosomal mean depth in Normal sample
        write_vcf_filter(os,get_label(HighDepth),oss.str().c_str());

        const StreamScoper ss(os);
        os << std::fixed << std::setprecision(2);

        for (const auto& val : dopt.chrom_depth)
        {
            const std::string& chrom(val.first);
            const double max_depth(opt.max_depth_factor*val.second);
            os << "##MaxDepth_" << chrom << '=' << max_depth << "\n";
        }
    }
}



strelka_streams::
strelka_streams(
    const strelka_options& opt,
    const strelka_deriv_options& dopt,
    const prog_info& pinfo,
    const bam_header_t* const header,
    const strelka_sample_info& ssi)
    : base_t(opt,pinfo,ssi)
{
    {
        using namespace STRELKA_SAMPLE_TYPE;
        if (opt.is_realigned_read_file)
        {
            _realign_bam_ptr[NORMAL].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.realigned_read_filename,"normal sample realigned-read BAM",header));
        }

        if (opt.is_tumor_realigned_read())
        {
            _realign_bam_ptr[TUMOR].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.tumor_realigned_read_filename,"tumor sample realigned-read BAM",header));
        }
    }

    if (opt.is_somatic_snv())
    {
        const char* const cmdline(opt.cmdline.c_str());

        std::ofstream* fosptr(new std::ofstream);
        _somatic_snv_osptr.reset(fosptr);
#ifdef DEBUG_HEADER
        std::ostream& fos = std::cout;
#else
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.somatic_snv_filename,"somatic-snv",opt.is_clobber,fos);
#endif

        if (! opt.sfilter.is_skip_header)
        {
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
            fos << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Combined depth across samples\">\n";
            fos << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n";
            fos << "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">\n";
            fos << "##INFO=<ID=ALTPOS,Number=1,Type=Integer,Description=\"Tumor alternate allele read position median\">\n";
            fos << "##INFO=<ID=ALTMAP,Number=1,Type=Integer,Description=\"Tumor alternate allele read position MAP\">\n";
            fos << "##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref read-position in the tumor\">\n";
            fos << "##INFO=<ID=SNVSB,Number=1,Type=Float,Description=\"Somatic SNV site strand bias\">\n";
            fos << "##INFO=<ID=PNOISE,Number=1,Type=Float,Description=\"Fraction of panel containing non-reference noise at this site\">\n";
            fos << "##INFO=<ID=PNOISE2,Number=1,Type=Float,Description=\"Fraction of panel containing more than one non-reference noise obs at this site\">\n";
            fos << "##INFO=<ID=VQSR,Number=1,Type=Float,Description=\"Recalibrated quality score expressing the phred scaled probability of the somatic call being a FP observation.\">\n";

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
                using namespace STRELKA_VCF_FILTERS;
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
                {
                    std::ostringstream oss;
                    oss << "The empirically fitted VQSR score is less than " << opt.sfilter.minimumQscore;
                    write_vcf_filter(fos, get_label(LowQscore), oss.str().c_str());
                }
            }

            write_shared_vcf_header_info(opt.sfilter,dopt.sfilter,fos);

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

        open_ofstream(pinfo,opt.somatic_indel_filename,"somatic-indel",opt.is_clobber,fos);

        if (! opt.sfilter.is_skip_header)
        {
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

            fos << "##FORMAT=<ID=DP" << opt.sfilter.indelRegionFlankSize << ",Number=1,Type=Float,Description=\"Average tier1 read depth within " << opt.sfilter.indelRegionFlankSize << " bases\">\n";
            fos << "##FORMAT=<ID=FDP" << opt.sfilter.indelRegionFlankSize << ",Number=1,Type=Float,Description=\"Average tier1 number of basecalls filtered from original read depth within " << opt.sfilter.indelRegionFlankSize << " bases\">\n";
            fos << "##FORMAT=<ID=SUBDP" << opt.sfilter.indelRegionFlankSize << ",Number=1,Type=Float,Description=\"Average number of reads below tier1 mapping quality threshold aligned across sites within " << opt.sfilter.indelRegionFlankSize << " bases\">\n";

            // FILTERS:
            {
                using namespace STRELKA_VCF_FILTERS;
                {
                    std::ostringstream oss;
                    oss << "Sequence repeat of more than " << opt.sfilter.indelMaxRefRepeat << "x in the reference sequence";
                    write_vcf_filter(fos, get_label(Repeat), oss.str().c_str());
                }
                {
                    std::ostringstream oss;
                    oss << "Indel overlaps an interrupted homopolymer longer than " << opt.sfilter.indelMaxIntHpolLength << "x in the reference sequence";
                    write_vcf_filter(fos, get_label(iHpol), oss.str().c_str());
                }
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
            write_shared_vcf_header_info(opt.sfilter,dopt.sfilter,fos);

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

        open_ofstream(pinfo,opt.somatic_callable_filename,"somatic-callable-regions",opt.is_clobber,fos);

        if (! opt.sfilter.is_skip_header)
        {
            fos << "track name=\"StrelkaCallableSites\"\t"
                << "description=\"Sites with sufficient information to call somatic alleles at 10% frequency or greater.\"\n";
        }
    }
}

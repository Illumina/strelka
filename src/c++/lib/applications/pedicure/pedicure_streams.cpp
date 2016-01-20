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

#include "pedicure_vcf_locus_info.hh"
#include "pedicure_streams.hh"

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
    const denovo_filter_options& opt,
    const denovo_filter_deriv_options& dopt,
    std::ostream& os)
{
    if (dopt.is_max_depth())
    {
        using namespace PEDICURE_VCF_FILTERS;

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



pedicure_streams::
pedicure_streams(
    const pedicure_options& opt,
    const pedicure_deriv_options& dopt,
    const prog_info& pinfo,
    const bam_hdr_t* const header,
    const PedicureSampleSetSummary& ssi)
    : base_t(opt,pinfo,ssi)
{
    assert(! opt.is_realigned_read_file);


    const char* const cmdline(opt.cmdline.c_str());

    {
        std::ofstream* fosptr(new std::ofstream);
        _denovo_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.denovo_filename,"denovo-small-variants",opt.is_clobber,fos);
    }

    if (! opt.dfilter.is_skip_header)
    {
        std::ostream& os(*_denovo_osptr);
//    	std::ostream& os(std::cout);

        write_vcf_audit(opt,pinfo,cmdline,header,os);
        os << "##content=pedicure snv calls\n"
           << "##germlineSnvTheta=" << opt.bsnp_diploid_theta << "\n";

        // INFO:
        os << "##INFO=<ID=DQ,Number=1,Type=Integer,Description=\"De novo quality score\">\n";
        os << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Combined depth across samples\">\n";
        os << "##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n";
        os << "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">\n";

        // Indel specific infos
        os << "##INFO=<ID=IHP,Number=1,Type=Integer,Description=\"Test\">\n";
        os << "##INFO=<ID=RU,Number=.,Type=String,Description=\"Test\">\n";
        os << "##INFO=<ID=RC,Number=1,Type=Integer,Description=\"Test\">\n";
        os << "##INFO=<ID=IC,Number=1,Type=Integer,Description=\"Test\">\n";

        // FORMAT:
        os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        os << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n";
        os << "##FORMAT=<ID=GQX,Number=1,Type=Integer,Description=\"Minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}\">\n";
        os << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n";
        os << "##FORMAT=<ID=DP2,Number=1,Type=Integer,Description=\"Read depth\">\n";
        os << "##FORMAT=<ID=FDP,Number=1,Type=Integer,Description=\"Filtered read depth\">\n";
        os << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">\n";
        os << "##FORMAT=<ID=PL,Number=.,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n";

        // FILTERS:
        {
            using namespace PEDICURE_VCF_FILTERS;
            {
                std::ostringstream oss;

                oss << "Calls that have GQX below " << 30;
                write_vcf_filter(os, get_label(lowGQX), oss.str().c_str());
            }
        }

        write_shared_vcf_header_info(opt.dfilter,dopt.dfilter,os);

        os << vcf_col_label() << "\tFORMAT";
        {
            const SampleInfoManager& si(opt.alignFileOpt.alignmentSampleInfo);
            unsigned parentCount = 1;
            for (unsigned sampleIndex(0); sampleIndex<si.size(); ++sampleIndex)
            {
                os << "\t" << PEDICURE_SAMPLETYPE::get_label(si.getSampleInfo(sampleIndex).stype);
                if (si.getSampleInfo(sampleIndex).stype==PEDICURE_SAMPLETYPE::PARENT)
                {
                    os << parentCount;
                    parentCount++;
                }

            }
        }
        os << "\n";
    }

    if (opt.is_denovo_callable())
    {
        std::ofstream* fosptr(new std::ofstream);
        _denovo_callable_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);

        open_ofstream(pinfo,opt.denovo_callable_filename,"denovo-callable-regions",opt.is_clobber,fos);

        // post samtools 1.0 tabix doesn't handle header information anymore, so take this out entirely:
#if 0
        if (! opt.dfilter.is_skip_header)
        {
            fos << "track name=\"DenovoCallableSites\"\t"
                << "description=\"Sites with sufficient information to call denovo variants.\"\n";
        }
#endif
    }
}

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

/// \author Chris Saunders
///

#include "starling_streams.hh"
#include "htsapi/bam_header_util.hh"

#include <cassert>

#include <fstream>
#include <iostream>


std::ostream*
starling_streams::
initialize_gvcf_file(
    const starling_options& opt,
    const prog_info& pinfo,
    const std::string& filename,
    const char* label,
    const bam_hdr_t& header)
{
    std::ofstream* fosPtr(new std::ofstream);
    open_ofstream(pinfo, filename, label, *fosPtr);

    if (not opt.gvcf.is_skip_header)
    {
        std::ostream& os(*fosPtr);
        const char* const cmdline(opt.cmdline.c_str());

        write_vcf_audit(opt,pinfo,cmdline,header,os);

        os << "##content=" << pinfo.name() << " germline small-variant calls\n";
    }
    return fosPtr;
}



starling_streams::
starling_streams(
    const starling_options& opt,
    const prog_info& pinfo,
    const std::vector<std::reference_wrapper<const bam_hdr_t>>& bamHeaders,
    const unsigned sampleCount)
    : base_t(opt, pinfo, sampleCount)
{
    for (const bam_hdr_t& bamHeader : bamHeaders)
    {
        _sampleNames.push_back(get_bam_header_sample_name(bamHeader));
    }

    assert(not bamHeaders.empty());
    const bam_hdr_t& referenceHeader(bamHeaders.front());

    if (opt.gvcf.is_gvcf_output())
    {
        const std::string gvcfVariantsPath(opt.gvcf.outputPrefix+"variants.vcf");
        _gvcfVariantsStreamPtr.reset(initialize_gvcf_file(opt, pinfo, gvcfVariantsPath, "variants", referenceHeader));
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            std::ostringstream sampleTag;
            sampleTag << "S" << (sampleIndex+1);
            const std::string gvcfSamplePath(opt.gvcf.outputPrefix+"genome." + sampleTag.str() + ".vcf");
            _gvcfSampleStreamPtr.emplace_back(
                initialize_gvcf_file(opt, pinfo, gvcfSamplePath, sampleTag.str().c_str(), referenceHeader));
        }
    }

    if (opt.is_realigned_read_file())
    {
        const unsigned inputAlignFileCount(bamHeaders.size());
        for (unsigned alignFileIndex(0); alignFileIndex < inputAlignFileCount; alignFileIndex++)
        {
            std::ostringstream rfile;
            rfile << opt.realignedReadFilenamePrefix << ".S" << alignFileIndex << ".bam";
            _realign_bam_ptr[alignFileIndex] = initialize_realign_bam(rfile.str(), bamHeaders[alignFileIndex]);
        }
    }
}

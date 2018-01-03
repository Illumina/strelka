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

/// \file
/// \author Chris Saunders
///

#include "starling_common/starling_streams_base.hh"
#include "blt_util/digt.hh"
#include "htsapi/vcf_util.hh"

#include <cassert>
#include <ctime>

#include <fstream>
#include <iostream>
#include <sstream>



static
void
open_ofstream(const prog_info& pinfo,
              const std::string& filename,
              const char* label,
              std::ofstream& fos)
{
    fos.open(filename.c_str());
    if (!fos)
    {
        std::ostringstream oss;
        oss << label << " file can't be opened: " << filename;
        pinfo.usage(oss.str().c_str());
    }
}



static
void
write_vcf_audit(
    const starling_base_options& opt,
    const prog_info& pinfo,
    const char* const cmdline,
    const bam_hdr_t& header,
    std::ostream& os)
{
    static const unsigned bufferSize(256);
    char timeBuffer[bufferSize];

    const time_t result(time(nullptr));
    strftime(timeBuffer, bufferSize, "%c", localtime(&result));

    os << "##fileformat=VCFv4.1\n";
    os << "##fileDate=" << vcf_fileDate << "\n";
    os << "##source=" << pinfo.name() << "\n";
    os << "##source_version=" << pinfo.version() << "\n";
    os << "##startTime=" << timeBuffer << "\n";
    os << "##cmdline=" << cmdline << "\n";
    if (not opt.referenceFilename.empty())
    {
        os << "##reference=file://" << opt.referenceFilename << "\n";
    }
    for (int32_t i(0); i<header.n_targets; ++i)
    {
        os << "##contig=<ID=" << header.target_name[i]
           << ",length=" << header.target_len[i] << ">\n";
    }
}



void
starling_streams_base::
write_vcf_audit(const starling_base_options& opt,
                const prog_info& pinfo,
                const char* const cmdline,
                const bam_hdr_t& header,
                std::ostream& os)
{
    ::write_vcf_audit(opt,pinfo,cmdline,header,os);
}



void
starling_streams_base::
open_ofstream(const prog_info& pinfo,
              const std::string& filename,
              const char* label,
              std::ofstream& fos)
{
    ::open_ofstream(pinfo,filename,label,fos);
}



std::unique_ptr<bam_dumper>
starling_streams_base::
initialize_realign_bam(
    const std::string& filename,
    const bam_hdr_t& header)
{
    // \TODO consider putting extra info into BAM header:
    //
    //fp->header = bam_header_dup((const bam_header_t*)aux);
    //fos << "@PG\tID:" << pinfo.name() << "\tVN:" << pinfo.version() << "\tCL:" << cmdline << "\n";

    return std::unique_ptr<bam_dumper>(new bam_dumper(filename.c_str(),header));
}



starling_streams_base::
starling_streams_base(
    const unsigned sampleCount)
    : _realign_bam_ptr(sampleCount),
      _sampleCount(sampleCount)
{
    assert(_sampleCount > 0);
}

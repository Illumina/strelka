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

/// \file

/// \author Chris Saunders
///

#include "htsapi/bam_dumper.hh"
#include "starling_common/starling_streams_base.hh"

#include <cassert>

#include <fstream>
#include <iostream>



bam_dumper*
starling_streams_base::
initialize_realign_bam(
    const bool is_clobber,
    const prog_info& pinfo,
    const std::string& filename,
    const char* label,
    const bam_hdr_t& header)
{
    // \TODO consider putting extra info into BAM header:
    //
    //fp->header = bam_header_dup((const bam_header_t*)aux);
    //fos << "@PG\tID:" << pinfo.name() << "\tVN:" << pinfo.version() << "\tCL:" << cmdline << "\n";

    if (! is_clobber)   // weak clobber test:
    {
        std::ofstream fos;
        open_ofstream(pinfo,filename,label,is_clobber,fos);
    }
    return new bam_dumper(filename.c_str(),header);
}



std::ostream*
starling_streams_base::
initialize_candidate_indel_file(const starling_base_options& opt,
                                const prog_info& pinfo,
                                const std::string& filename)
{
    const char* const cmdline(opt.cmdline.c_str());

    std::ofstream* fosptr(new std::ofstream);
    std::ofstream& fos(*fosptr);
    open_ofstream(pinfo,filename,"candidate-indel",opt.is_clobber,fos);

    fos << "# ** " << pinfo.name();
    fos << " candidate-indel file **\n";
    write_file_audit(opt,pinfo,cmdline,fos);
    fos << "#\n";
    fos << "#$ COLUMNS seq_name pos type length seq\n";

    return fosptr;
}



starling_streams_base::
starling_streams_base(const starling_base_options& opt,
                      const prog_info& pinfo,
                      const SampleSetSummary& si)
    : base_t(opt,pinfo)
    , _n_samples(si.size())
{
    assert((_n_samples>0) && (_n_samples<=MAX_SAMPLE));

    if (opt.is_write_candidate_indels())
    {
        _candidate_indel_osptr.reset(initialize_candidate_indel_file(opt,pinfo,opt.candidate_indel_filename));
    }
}

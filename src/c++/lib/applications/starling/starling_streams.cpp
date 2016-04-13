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
    const bam_hdr_t& header,
    std::unique_ptr<std::ostream>& os_ptr_auto)
{
    std::ostream* osptr(&std::cout);
    if (filename != "-")
    {
        std::ofstream* fos_ptr(new std::ofstream);
        open_ofstream(pinfo,filename,"gvcf",opt.is_clobber,*fos_ptr);
        os_ptr_auto.reset(fos_ptr);
        osptr=os_ptr_auto.get();
    }
    std::ostream& os(*osptr);

    if (! opt.gvcf.is_skip_header)
    {
        const char* const cmdline(opt.cmdline.c_str());

        write_vcf_audit(opt,pinfo,cmdline,header,os);

        os << "##content=" << pinfo.name() << " small-variant calls\n"
           << "##SnvTheta=" << opt.bsnp_diploid_theta << "\n"
           << "##IndelTheta=" << opt.bindel_diploid_theta << "\n";
    }
    return osptr;
}



starling_streams::
starling_streams(
    const starling_options& opt,
    const prog_info& pinfo,
    const bam_hdr_t& header,
    const SampleSetSummary& ssi)
    : base_t(opt,pinfo,ssi)
{
    assert(_n_samples == 1);
    _gvcf_osptr = nullptr;
    _sampleName = get_bam_header_sample_name(header);

    if (opt.gvcf.is_gvcf_output())
    {
        _gvcf_osptr = initialize_gvcf_file(opt,pinfo,opt.gvcf.out_file,header,_gvcf_osptr_auto);
    }

    if (opt.is_realigned_read_file)
    {
        _realign_bam_ptr[0].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.realigned_read_filename,"realigned-read BAM/CRAM",header));
    }
}

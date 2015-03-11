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

/// \author Chris Saunders
///

#include "starling_streams.hh"

#include <fstream>
#include <iostream>


std::ostream*
starling_streams::
initialize_gvcf_file(
    const starling_options& opt,
    const prog_info& pinfo,
    const std::string& filename,
    const bam_header_t* const header,
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
    const bam_header_t* const header,
    const sample_info& ssi)
    : base_t(opt,pinfo,ssi)
{
    for (unsigned i(0); i<_n_samples; ++i) _gvcf_osptr[i] = nullptr;

    if (opt.gvcf.is_gvcf_output())
    {
        _gvcf_osptr[0] = initialize_gvcf_file(opt,pinfo,opt.gvcf.out_file,header,_gvcf_osptr_auto[0]);
    }

    if (opt.is_realigned_read_file)
    {
        _realign_bam_ptr[0].reset(initialize_realign_bam(opt.is_clobber,pinfo,opt.realigned_read_filename,"realigned-read BAM",header));
    }
}

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

#include "snoise_streams.hh"

#include <iostream>


const sample_info starling_sample_info;



snoise_streams::
snoise_streams(
    const snoise_options& opt,
    const prog_info& pinfo,
    const bam_header_t* const header)
    : base_t(opt,pinfo,starling_sample_info)
{
    if (true)
    {
        const char* const cmdline(opt.cmdline.c_str());

        std::ofstream* fosptr(new std::ofstream);
        _snoise_osptr.reset(fosptr);
        std::ofstream& fos(*fosptr);
        open_ofstream(pinfo,opt.snoise_filename,"snoise",opt.is_clobber,fos);

        static const bool is_skip_header(false);

        if (! is_skip_header)
        {
            write_vcf_audit(opt,pinfo,cmdline,header,fos);
            fos << "##content=strelka sequencing noise extraction\n";

            // INFO:

            // FORMAT:

            // FILTERS:


            fos << vcf_col_label() << "\tFORMAT";
            for (unsigned s(0); s<STRELKA_SAMPLE_TYPE::SIZE; ++s)
            {
                fos << "\t" << STRELKA_SAMPLE_TYPE::get_label(s);
            }
            fos << "\n";
        }
    }
}

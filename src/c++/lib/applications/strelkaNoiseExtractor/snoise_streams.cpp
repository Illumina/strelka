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
#include "htsapi/bam_header_util.hh"
#include "htsapi/vcf_util.hh"

#include <fstream>
#include <iostream>



snoise_streams::
snoise_streams(
    const snoise_options& opt,
    const prog_info& pinfo,
    const bam_header_t* const header,
    const SampleSetSummary& ssi)
    : base_t(opt,pinfo,ssi),
      _snoise_osptr(&std::cout)
{
    const char* const cmdline(opt.cmdline.c_str());

    std::ostream& fos(*_snoise_osptr);

    if (! opt.is_skip_header)
    {
        write_vcf_audit(opt,pinfo,cmdline,header,fos);
        fos << "##content=strelka sequencing noise extraction\n";

        // INFO:

        // FORMAT:
        fos << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered read depth\">\n";
        fos << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">\n";

        // FILTERS:

        const std::string sample_name = get_bam_header_sample_name(header->text);
        fos << vcf_col_label() << "\tFORMAT\t" << sample_name << "\n";
    }
}

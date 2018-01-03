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
    const bam_hdr_t& header,
    const unsigned sampleCount)
    : base_t(sampleCount),
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
        fos << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob " << opt.readConfidentSupportThreshold.strval() << " or higher that read contains indicated allele vs all other intersecting indel alleles)\">\n";

        // FILTERS:

        const std::string sample_name = get_bam_header_sample_name(header);
        fos << vcf_col_label() << "\tFORMAT\t" << sample_name << "\n";
    }
}

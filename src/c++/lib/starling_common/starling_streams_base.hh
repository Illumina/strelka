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

#pragma once

#include "blt_util/prog_info.hh"
#include "htsapi/bam_util.hh"
#include "htsapi/bam_dumper.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_types.hh"

#include <iosfwd>
#include <memory>
#include <string>
#include <vector>


struct starling_streams_base
{
    explicit
    starling_streams_base(
        const unsigned sampleCount);

    bam_dumper*
    realign_bam_ptr(const unsigned sampleIndex) const
    {
        return _realign_bam_ptr[sampleIndex].get();
    }

    unsigned
    getSampleCount() const
    {
        return _sampleCount;
    }

protected:
    std::unique_ptr<bam_dumper>
    initialize_realign_bam(
        const std::string& filename,
        const bam_hdr_t& header);

    static
    void
    open_ofstream(const prog_info& pinfo,
                  const std::string& filename,
                  const char* label,
                  std::ofstream& fos);

    /// write the first few meta-data lines for a vcf file:
    ///
    static
    void
    write_vcf_audit(const starling_base_options& opt,
                    const prog_info& pinfo,
                    const char* const cmdline,
                    const bam_hdr_t& header,
                    std::ostream& os);

    std::vector<std::unique_ptr<bam_dumper>> _realign_bam_ptr;
private:
    unsigned _sampleCount;
};

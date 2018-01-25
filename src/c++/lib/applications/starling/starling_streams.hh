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

#pragma once

#include "starling_shared.hh"
#include "starling_common/starling_streams_base.hh"


struct starling_streams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    starling_streams(
        const starling_options& opt,
        const prog_info& pinfo,
        const std::vector<std::reference_wrapper<const bam_hdr_t>>& bamHeaders,
        const std::vector<std::string>& sampleNames);

    std::ostream&
    gvcfSampleStream(const unsigned sampleIndex) const
    {
        assert(sampleIndex < getSampleCount());
        return *(_gvcfSampleStreamPtr[sampleIndex]);
    }

    std::ostream&
    variantsVCFStream() const
    {
        return *(_variantsVCFStreamPtr);
    }

    const std::vector<std::string>&
    getSampleNames() const
    {
        return _sampleNames;
    }

private:
    static
    std::unique_ptr<std::ostream>
    initializeGermlineVCFStream(
        const starling_options& opt,
        const prog_info& pinfo,
        const std::string& filename,
        const char* label,
        const bam_hdr_t& header);

    std::unique_ptr<std::ostream> _variantsVCFStreamPtr;
    std::vector<std::unique_ptr<std::ostream>> _gvcfSampleStreamPtr;
    std::vector<std::string> _sampleNames;
};

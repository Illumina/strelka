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

#pragma once

#include "strelka_shared.hh"

#include "starling_common/starling_streams_base.hh"
#include "strelka_common/StrelkaSampleSetSummary.hh"



struct strelka_streams : public starling_streams_base
{
    typedef starling_streams_base base_t;

    strelka_streams(
        const strelka_options& opt,
        const strelka_deriv_options& dopt,
        const prog_info& pinfo,
        const bam_hdr_t& bam_header,
        const StrelkaSampleSetSummary& ssi);

    std::ostream*
    somatic_snv_osptr() const
    {
        return _somatic_snv_osptr.get();
    }

    std::ostream*
    somatic_indel_osptr() const
    {
        return _somatic_indel_osptr.get();
    }

    std::ostream*
    somatic_callable_osptr() const
    {
        return _somatic_callable_osptr.get();
    }

private:
    std::unique_ptr<std::ostream> _somatic_snv_osptr;
    std::unique_ptr<std::ostream> _somatic_indel_osptr;
    std::unique_ptr<std::ostream> _somatic_callable_osptr;
};

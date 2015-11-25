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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_common/blt_streams.hh"
#include "htsapi/bam_dumper.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_types.hh"

#include <vector>
#include "SampleSetSummary.hh"


struct starling_streams_base : public blt_streams
{
    typedef blt_streams base_t;

    starling_streams_base(
        const starling_base_options& opt,
        const prog_info& pinfo,
        const SampleSetSummary& si);

    bam_dumper*
    realign_bam_ptr(const unsigned sample_no) const
    {
        return _realign_bam_ptr[sample_no].get();
    }

    std::ostream*
    candidate_indel_osptr() const
    {
        return _candidate_indel_osptr.get();
    }

    std::ostream*
    variant_window_osptr(const unsigned window_no) const
    {
        return _window_osptr.at(window_no).get();
    }

protected:
    bam_dumper*
    initialize_realign_bam(const bool is_clobber,
                           const prog_info& pinfo,
                           const std::string& filename,
                           const char* label,
                           const bam_header_t* const header);

    static
    std::ostream*
    initialize_candidate_indel_file(const starling_base_options& client_opt,
                                    const prog_info& pinfo,
                                    const std::string& filename);

    static
    std::ostream*
    initialize_window_file(const starling_base_options& opt,
                           const prog_info& pinfo,
                           const avg_window_data& awd,
                           const SampleSetSummary& si);

    std::unique_ptr<bam_dumper> _realign_bam_ptr[MAX_SAMPLE];
private:
    std::unique_ptr<std::ostream> _candidate_indel_osptr;
protected:
    unsigned _n_samples;
private:
    std::vector<std::shared_ptr<std::ostream>> _window_osptr;
};

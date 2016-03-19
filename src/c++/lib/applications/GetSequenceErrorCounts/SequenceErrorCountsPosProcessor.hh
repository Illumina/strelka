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

#include "SequenceErrorCountsOptions.hh"
#include "SequenceErrorCountsStreams.hh"
#include "errorAnalysis/SequenceErrorCounts.hh"
#include "starling_common/starling_pos_processor_base.hh"


///
///
struct SequenceErrorCountsPosProcessor : public starling_pos_processor_base
{
    typedef starling_pos_processor_base base_t;

    SequenceErrorCountsPosProcessor(
        const SequenceErrorCountsOptions& opt,
        const SequenceErrorCountsDerivOptions& dopt,
        const reference_contig_segment& ref,
        const SequenceErrorCountsStreams& streams);

    virtual
    ~SequenceErrorCountsPosProcessor();

    void reset();

private:

    void
    process_pos_variants_impl(const pos_t pos) override
    {
        process_pos_indel_single_sample_digt(pos,0);
    }

    void process_pos_indel_single_sample_digt(
        const pos_t pos,
        const unsigned sample_no);

    void
    write_counts(const pos_range&) const override {}

    const SequenceErrorCountsOptions& _opt;
    const SequenceErrorCountsDerivOptions& _dopt;
    const SequenceErrorCountsStreams& _streams;

    SequenceErrorCounts _counts;

    double _max_candidate_normal_sample_depth = -1;

};

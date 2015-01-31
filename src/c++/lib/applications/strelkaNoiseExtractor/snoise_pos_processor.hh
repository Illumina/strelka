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

///
/// \author Chris Saunders
///

#pragma once

#include "snoise_streams.hh"
#include "starling_common/starling_pos_processor_base.hh"


///
///
struct snoise_pos_processor : public starling_pos_processor_base
{
    typedef starling_pos_processor_base base_t;

    snoise_pos_processor(
        const starling_options& opt,
        const starling_deriv_options& dopt,
        const reference_contig_segment& ref,
        const snoise_streams& client_io);

private:
    void
    process_pos_variants_impl(const pos_t pos) override
    {
        process_pos_snp_snoise(pos);
    }

    void
    process_pos_snp_snoise(const pos_t pos);

    void
    write_counts(const pos_range& output_report_range) const;

    const snoise_streams& _client_io;
};

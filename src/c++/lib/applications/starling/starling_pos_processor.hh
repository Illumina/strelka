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

#include "starling_common/starling_pos_processor_base.hh"



///
///
struct starling_pos_processor : public starling_pos_processor_base {

    typedef starling_pos_processor_base base_t;

    starling_pos_processor(const starling_options& client_opt,
                           const starling_deriv_options& client_dopt,
                           const reference_contig_segment& ref,
                           const starling_streams_base& client_io)
        : base_t(client_opt,client_dopt,ref,client_io,1) {}

private:
    void
    process_pos_variants(const pos_t pos) {
        process_pos_indel_single_sample(pos,0);
        process_pos_snp_single_sample(pos,0);
    }

    void
    write_counts(const pos_range& output_report_range) const;
};

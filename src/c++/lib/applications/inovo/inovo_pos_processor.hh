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

///
/// note coding convention for all ranges '_pos fields' is:
/// XXX_begin_pos is zero-indexed position at the beginning of the range
/// XXX_end_pos is zero-index position 1 step after the end of the range
///

#pragma once


#include "DenovoCallableProcessor.hh"
#include "inovo_shared.hh"
#include "inovo_streams.hh"

#include "starling_common/starling_pos_processor_base.hh"


///
///
struct inovo_pos_processor : public starling_pos_processor_base
{
    typedef starling_pos_processor_base base_t;

    inovo_pos_processor(
        const inovo_options& opt,
        const inovo_deriv_options& dopt,
        const reference_contig_segment& ref,
        const inovo_streams& streams);

private:

    void
    process_pos_variants_impl(const pos_t pos) override;

    void
    process_pos_snp_denovo(const pos_t pos);

    void
    process_pos_indel_denovo(const pos_t pos);

    void
    write_counts(const pos_range& output_report_range) const override;

    /////////////////////////////

    // keep some of the original inovo classes handy so we don't
    // have to down-cast:
    const inovo_options& _opt;
    const inovo_deriv_options& _dopt;
    const inovo_streams& _streams;

    DenovoCallableProcessor _icallProcessor;

//    extra_position_data _tier2_epd[MAX_SAMPLE];

//    unsigned _indelRegionIndexNormal;
//    unsigned _indelRegionIndexTumor;
};

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///

#include "blt_util/math_util.hh"
#include "starling_common/starling_indel_call_pprob_digt.hh"
#include "starling_common/starling_shared.hh"

#include <cmath>

#include <iostream>



std::ostream&
operator<<(std::ostream& os, const avg_window_data& awd) {

    os << "flank_size: " << awd.flank_size << " file: " << awd.filename << "\n";
    return os;
}



starling_deriv_options::
starling_deriv_options(const starling_options& opt,
                       const reference_contig_segment& ref)
    : base_t(opt,ref.end())
    , sal(opt.max_realignment_candidates)
    , _incaller(new indel_digt_caller(opt.bindel_diploid_theta))
{

    indel_nonsite_match_lnp=std::log(opt.indel_nonsite_match_prob);
    if(opt.is_tier2_indel_nonsite_match_prob) {
        tier2_indel_nonsite_match_lnp=std::log(opt.tier2_indel_nonsite_match_prob);
    } else {
        tier2_indel_nonsite_match_lnp=indel_nonsite_match_lnp;
    }

    {
        // set genome_size for indel model:
        uint32_t genome_size;
        if(opt.is_user_genome_size) {
            genome_size = opt.user_genome_size;
        } else {
            assert(0);
            //            genome_size = get_ref_seq_known_size(ref.seq());
        }

        // get read path posterior probs:
        const double site_prior(1./(2.*static_cast<double>(genome_size)));
        site_lnprior=std::log(site_prior);
        nonsite_lnprior=log1p_switch(-site_prior);
    }
}



// dtor must be here so that auto_ptr works correctly:
starling_deriv_options::
~starling_deriv_options() {}



void
starling_read_counts::
report(std::ostream& os) const {
    blt_read_counts::report(os);
    os << "STARLING_READ_COUNTS"
       << " normal_indel_used: " << normal_indel_used
       << " normal_indel_intersect: " << normal_indel_intersect
       << " grouper_indel_used: " << grouper_indel_used
       << " grouper_indel_intersect: " << grouper_indel_intersect
       << " grouper_unused: " << grouper_unused << "\n";
}

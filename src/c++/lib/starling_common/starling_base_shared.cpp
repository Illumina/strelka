//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

#include "starling_common/starling_base_shared.hh"

#include "blt_util/math_util.hh"
#include "calibration/IndelErrorModel.hh"
#include "htsapi/bam_streamer.hh"
#include "starling_common/AlleleGroupGenotype.hh"

#include <cmath>

#include <iostream>



starling_base_deriv_options::
starling_base_deriv_options(const starling_base_options& opt)
    : base_t(opt)
    , sal(opt.max_realignment_candidates)
    , variant_window_first_stage(0)
    , variant_window_last_stage(0)
    , countCache(opt.indel_candidate_signal_test_alpha)
    , logIndelRefErrorFactor(std::log(opt.indelRefErrorFactor))
    , _indelErrorModel(new IndelErrorModel(opt.getAlignmentFileOptions().alignmentFilenames, opt.indel_error_model_name,opt.indelErrorModelFilenames))
    , _indelGenotypePriors(new GenotypePriorSet(opt.thetaFilename))
{
    indel_nonsite_match_lnp=std::log(opt.indel_nonsite_match_prob);
    if (opt.tier2.is_tier2_indel_nonsite_match_prob)
    {
        tier2_indel_nonsite_match_lnp=std::log(opt.tier2.tier2_indel_nonsite_match_prob);
    }
    else
    {
        tier2_indel_nonsite_match_lnp=indel_nonsite_match_lnp;
    }

    {
        // set genome_size for indel model:
        uint32_t genome_size;
        if (opt.is_user_genome_size)
        {
            genome_size = opt.user_genome_size;
        }
        else
        {
            assert(0);
            //            genome_size = get_ref_seq_known_size(ref.seq());
        }

        // get read path posterior probs:
        const double site_prior(1./(2.*static_cast<double>(genome_size)));
        site_lnprior=std::log(site_prior);
        nonsite_lnprior=log1p_switch(-site_prior);
    }

    // register post-call stages:
    //

    // no post-call stages remain in latest design!
}



// dtor is required here for unique_ptr
starling_base_deriv_options::
~starling_base_deriv_options() {}



void
starling_read_counts::
report(std::ostream& os) const
{
    blt_read_counts::report(os);
    os << "STARLING_READ_COUNTS"
       << " normal_indel_used: " << normal_indel_used
       << " normal_indel_intersect: " << normal_indel_intersect
       << "\n";
}

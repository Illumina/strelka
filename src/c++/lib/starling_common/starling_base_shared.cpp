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
    randomBaseMatchLogProb=std::log(opt.randomBaseMatchProb);
    if (opt.tier2.isRandomBaseMatchProb)
    {
        tier2RandomBaseMatchLogProb=std::log(opt.tier2.randomBaseMatchProb);
    }
    else
    {
        tier2RandomBaseMatchLogProb=randomBaseMatchLogProb;
    }

    {
        // Used to evaluate the chance that a read is mismapped, independent of any read mapper assertions
        // - thus accounting for things like population specific sequence, etc. which may not be part of
        // the mapper's MAPQ model
        //
        // The correct mapping prior was originally the chance that a read will be correctly mapped at random, so
        // was set to 1/(2*genome_size). The usage/concept has evolved to the point where this is now an arbitrary
        // constant rather than something based on the real genome size
        //
        // TODO Note this is legacy logic that is targeted for replacement
        correctMappingLogPrior=std::log(1.7e-10);
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

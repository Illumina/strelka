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

#include "snoise_pos_processor.hh"

#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>



snoise_pos_processor::
snoise_pos_processor(
    const snoise_options& opt,
    const starling_base_deriv_options& dopt,
    const reference_contig_segment& ref,
    const snoise_streams& fileStreams,
    RunStatsManager& statsManager)
    : base_t(opt, dopt, ref, fileStreams, opt.alignFileOpt.alignmentFilenames.size(), statsManager),
      _fileStreams(fileStreams)
{
    // setup indel buffer:
    {
        double maxIndelCandidateDepthSumOverNormalSamples(-1.);

        if (opt.max_candidate_indel_depth > 0.)
        {
            if (maxIndelCandidateDepthSumOverNormalSamples > 0.)
            {
                maxIndelCandidateDepthSumOverNormalSamples = std::min(maxIndelCandidateDepthSumOverNormalSamples,static_cast<double>(opt.max_candidate_indel_depth));
            }
            else
            {
                maxIndelCandidateDepthSumOverNormalSamples = opt.max_candidate_indel_depth;
            }
        }

        getIndelBuffer().setMaxCandidateDepth(maxIndelCandidateDepthSumOverNormalSamples);

        const unsigned sampleCount(getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            sample_info& sif(sample(sampleIndex));
            getIndelBuffer().registerSample(sif.estdepth_buff, sif.estdepth_buff_tier2, true);
        }

        getIndelBuffer().finalizeSamples();
    }
}



void
snoise_pos_processor::
process_pos_snp_snoise(
    const pos_t pos)
{
    const unsigned sample_no(0);
    const sample_info& sif(sample(sample_no));

    const snp_pos_info& pi(sif.cleanedPileup.rawPileup());
    const snp_pos_info& good_pi(sif.cleanedPileup.cleanedPileup());
    const pos_t output_pos(pos+1);


    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no != 0) return;

    if (pi.calls.empty()) return;

    // make early filtration decision -- then get the allele distribution:
    static const unsigned min_used_calls(12);

    const unsigned n_used_calls(good_pi.calls.size());
    if (n_used_calls < min_used_calls) return;

    const auto ref_id(base_to_id(good_pi.get_ref_base()));

    // don't process sites where the reference base is ambiguous
    if (ref_id >= BASE_ID::ANY) return;

    std::array<unsigned,BASE_ID::ANY> base_count;
    std::fill(base_count.begin(),base_count.end(),0);

    for (const auto& bc : good_pi.calls)
    {
        assert(bc.base_id<BASE_ID::ANY);
        base_count[bc.base_id]++;
    }

    const unsigned ref_count(base_count[ref_id]);
    if (ref_count == n_used_calls) return;

    unsigned alt_id( (ref_id==0) ? 1 : 0);
    for (unsigned i(1); i<BASE_ID::ANY; ++i)
    {
        if (i == ref_id) continue;
        if (i == alt_id) continue;
        if (base_count[i] > base_count[alt_id]) alt_id = i;
    }
    const unsigned alt_count(base_count[alt_id]);

#if 0
    const double ref_ratio(static_cast<double>(ref_count)/_site_info.usedBasecallCount);

    if (ref_ratio > 0.2) return;
#endif

    const double alt_ratio(static_cast<double>(alt_count)/n_used_calls);

    static const double max_alt_ratio(0.2);

    /// too likely to be germline:
    if (alt_ratio > max_alt_ratio) return;

    {
        std::ostream& os(*_fileStreams.snoise_osptr());

        // CHROM POS ID:
        os << _chromName << '\t'
           << output_pos << '\t'
           << ".";

        //REF:
        os << '\t' << good_pi.get_ref_base();
        //ALT:
        os << '\t' << id_to_base(alt_id);

        //QUAL:
        os << "\t.";

        //FILTER:
        os << "\t.";

        //INFO:
        os << "\t.";


        //FORMAT:
        os << '\t'
           << "DP:AD";

        // SAMPLE:
        os << "\t";
        os << n_used_calls << ':' << ref_count << ',' << alt_count;
        os << "\n";
    }
}

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
    const snoise_streams& streams)
    : base_t(opt,dopt,ref,streams,1),
      _streams(streams)
{
    // setup indel syncronizers:
    {
        sample_info& normal_sif(sample(0));

        double max_candidate_normal_sample_depth(-1.);
#if 0
        if (dopt.gvcf.is_max_depth())
        {
            if (opt.max_candidate_indel_depth_factor > 0.)
            {
                max_candidate_normal_sample_depth = (opt.max_candidate_indel_depth_factor * dopt.gvcf.max_depth);
            }
        }
#endif

        if (opt.max_candidate_indel_depth > 0.)
        {
            if (max_candidate_normal_sample_depth > 0.)
            {
                max_candidate_normal_sample_depth = std::min(max_candidate_normal_sample_depth,static_cast<double>(opt.max_candidate_indel_depth));
            }
            else
            {
                max_candidate_normal_sample_depth = opt.max_candidate_indel_depth;
            }
        }

        indel_sync_data isdata;
        isdata.register_sample(normal_sif.indel_buff,normal_sif.estdepth_buff,normal_sif.estdepth_buff_tier2,
                               normal_sif.sample_opt, max_candidate_normal_sample_depth, 0);
        normal_sif.indel_sync_ptr.reset(new indel_synchronizer(opt, ref, isdata, 0));
    }
}



void
snoise_pos_processor::
write_counts(
    const pos_range& output_report_range) const
{
    std::ostream* report_osptr(get_report_osptr());
    if (nullptr==report_osptr) return;
    std::ostream& report_os(*report_osptr);

    const sample_info& sif(sample());

    report_os << std::setprecision(8);
    report_stream_stat(sif.ss,"ALLSITES_COVERAGE",output_report_range,report_os);
    report_stream_stat(sif.used_ss,"ALLSITES_COVERAGE_USED",output_report_range,report_os);

    if (_opt.is_ref_set())
    {
        report_stream_stat(sif.ssn,"NO_REF_N_COVERAGE",output_report_range,report_os);
        report_stream_stat(sif.used_ssn,"NO_REF_N_COVERAGE_USED",output_report_range,report_os);
    }
}



void
snoise_pos_processor::
process_pos_snp_snoise(
    const pos_t pos)
{
    const unsigned sample_no(0);
    const sample_info& sif(sample(sample_no));

    const snp_pos_info& pi(sif.cpi.rawPileup());
    const snp_pos_info& good_pi(sif.cpi.cleanedPileup());
    const pos_t output_pos(pos+1);


    // note multi-sample status -- can still be called only for one sample
    // and only for sample 0. working on generalization:
    //
    if (sample_no != 0) return;

    if (pi.calls.empty()) return;

    // make early filtration decision -- then get the allele distribution:
    constexpr unsigned min_used_calls(12);

    const unsigned n_used_calls(good_pi.calls.size());
    if (n_used_calls < min_used_calls) return;

    std::array<unsigned,N_BASE> base_count;
    std::fill(base_count.begin(),base_count.end(),0);

    for (const auto& bc : good_pi.calls)
    {
        assert(bc.base_id!=BASE_ID::ANY);
        base_count[bc.base_id]++;
    }

    const auto ref_id(base_to_id(good_pi.get_ref_base()));
    const unsigned ref_count(base_count[ref_id]);

    if (ref_count == n_used_calls) return;

    unsigned alt_id( (ref_id==0) ? 1 : 0);
    for (unsigned i(1); i<N_BASE; ++i)
    {
        if (i == ref_id) continue;
        if (i == alt_id) continue;
        if (base_count[i] > base_count[alt_id]) alt_id = i;
    }
    const unsigned alt_count(base_count[alt_id]);

#if 0
    const double ref_ratio(static_cast<double>(ref_count)/_site_info.n_used_calls);

    if (ref_ratio > 0.2) return;
#endif

    const double alt_ratio(static_cast<double>(alt_count)/n_used_calls);

    constexpr double max_alt_ratio(0.2);

    /// too likely to be germline:
    if (alt_ratio > max_alt_ratio) return;

    {
        std::ostream& os(*_streams.snoise_osptr());

        // CHROM POS ID:
        os << _chrom_name << '\t'
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

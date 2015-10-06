#pragma once
#include "starling_common/starling_base_shared.hh"
#include "blt_common/snp_pos_info.hh"

#include "blt_util/qscore.hh"
#include "gvcf_locus_info.hh"
#include <boost/utility.hpp>


class starling_continuous_variant_caller : private boost::noncopyable
{
    public:

    static void position_snp_call_continuous(
        const starling_base_options& opt,
        const snp_pos_info& good_pi,
        continuous_site_info& info);
    static void add_indel_call(
            const starling_base_options& opt,
            const indel_key& ik,
            const indel_data& id,
            const starling_indel_report_info& iri,
            const starling_indel_sample_report_info& isri,
            continuous_indel_info& info);
    static unsigned poisson_qscore(unsigned callCount, unsigned coverage, unsigned estimatedBaseCallQuality, unsigned maxQScore);

    static double strand_bias(unsigned fwdAlt, unsigned revAlt, unsigned fwdOther, unsigned revOther, double noise);
};

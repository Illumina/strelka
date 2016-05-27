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


#include "gvcf_locus_info.hh"
#include "blt_util/math_util.hh"
#include "common/Exceptions.hh"

#include "boost/math/distributions/binomial.hpp"
#include "boost/math/distributions.hpp"

#include <iostream>
#include <map>
#include <sstream>
#include <typeinfo>



void
GermlineVariantSimpleGenotypeInfo::
write_filters(std::ostream& os) const
{
    if (filters.none())
    {
        os << "PASS";
        return;
    }

    bool is_sep(false);
    for (unsigned i(0); i<GERMLINE_VARIANT_VCF_FILTERS::SIZE; ++i)
    {
        if (! filters.test(i)) continue;

        if (is_sep)
        {
            os << ";";
        }
        else
        {
            is_sep=true;
        }
        os << GERMLINE_VARIANT_VCF_FILTERS::get_label(i);
    }
}


pos_t
GermlineDiploidIndelCallInfo::
end() const
{
    pos_t result = 0;
    for (auto& x : _calls)
        result = std::max(result, x._ik.right_pos());
    return result;
}

void GermlineIndelSimpleGenotypeInfo::set_hap_cigar(
    const unsigned lead,
    const unsigned trail)
{
    using namespace ALIGNPATH;

    cigar.clear();
    if (lead)
    {
        cigar.push_back(path_segment(MATCH,lead));
    }
    if (_ik.delete_length())
    {
        cigar.push_back(path_segment(DELETE,_ik.delete_length()));
    }
    if (_ik.insert_length())
    {
        cigar.push_back(path_segment(INSERT,_ik.insert_length()));
    }
    if (trail)
    {
        cigar.push_back(path_segment(MATCH,trail));
    }
}


void
GermlineDiploidIndelSimpleGenotypeInfo::
computeEmpiricalScoringFeatures(
    const bool isComputeDevelopmentFeatures,
    const double chromDepth,
    const bool isHetalt)
{
    const double filteredLocusDepth(_isri.tier1Depth);
    const double locusDepth(_isri.mapqTracker.count);
    const double q30Depth(_isri.total_q30_reads());

    const double chromDepthFactor(safeFrac(1,chromDepth));
    const double filteredLocusDepthFactor(safeFrac(1,filteredLocusDepth));
    const double locusDepthFactor(safeFrac(1,locusDepth));
    const double q30DepthFactor(safeFrac(1,q30Depth));

    features.set(GERMLINE_INDEL_SCORING_FEATURES::QUAL, (_dindel.indel_qphred * chromDepthFactor));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::F_GQX, (gqx * chromDepthFactor));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::REFREP1, (_iri.ref_repeat_count));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::IDREP1, (_iri.indel_repeat_count));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::RULEN1, (_iri.repeat_unit.length()));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::AD0, (_isri.n_q30_ref_reads * chromDepthFactor));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::AD1, (_isri.n_q30_indel_reads * chromDepthFactor));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::AD2, (_isri.n_q30_alt_reads * chromDepthFactor));
    features.set(GERMLINE_INDEL_SCORING_FEATURES::F_DPI, (_isri.tier1Depth * chromDepthFactor));

    {
        // allele bias metrics
        const double r0(_isri.n_q30_ref_reads);
        const double r1(_isri.n_q30_indel_reads);
        const double r2(_isri.n_q30_alt_reads);

        // cdf of binomial prob of seeing no more than the number of 'allele A' reads out of A reads + B reads, given p=0.5
        // cdf of binomial prob of seeing no more than the number of 'allele B' reads out of A reads + B reads, given p=0.5
        double allelebiaslower;
        double allelebiasupper;
        if (isHetalt)
        {
            allelebiaslower = cdf(boost::math::binomial(r2+r1,0.5),r1);
            allelebiasupper = cdf(boost::math::binomial(r2+r1,0.5),r2);
        }
        else
        {
            allelebiaslower = cdf(boost::math::binomial(r0+r1,0.5),r0);
            allelebiasupper = cdf(boost::math::binomial(r0+r1,0.5),r1);
        }

        // +1e-30 to avoid log(0) in extreme cases
        features.set(GERMLINE_INDEL_SCORING_FEATURES::ABlower, (-std::log(allelebiaslower+1.e-30)));
        features.set(GERMLINE_INDEL_SCORING_FEATURES::AB, (-std::log(std::min(1.,2.*std::min(allelebiaslower,allelebiasupper))+1.e-30)));
    }

    // compute any experimental features not currently used in production
    if (isComputeDevelopmentFeatures)
    {
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQ, (gqx * chromDepthFactor));

        // ************ notes for Konrad ***************:

        // these two are anticipated to be used as two possibly new features for the indel scoring model:
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_MQ, (_isri.mapqTracker.getRMS()));
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, (_isri.mapqTracker.getZeroFrac()));

        // this should replace F_DPI
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_DPI_NORM, (filteredLocusDepth * locusDepthFactor));

        // this is a new feature that currently isn't captured in the model -- this expresses how surprising the depth
        // is relative to chromosome mean.
        //
        // This is the only feature we want to use chrom depth after everything is fixed up. It will be set to 1
        // for exome builds.
        //
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::TDP_NORM, (locusDepth * chromDepthFactor));

        // all of the features below are simply renormalized reqplacements of the current produciton feature set
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::QUAL_NORM, (_dindel.indel_qphred * filteredLocusDepthFactor));
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::F_GQX_NORM, (gqx * filteredLocusDepthFactor));

        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD0_NORM, (_isri.n_q30_ref_reads * q30DepthFactor));
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD1_NORM, (_isri.n_q30_indel_reads * q30DepthFactor));
        developmentFeatures.set(GERMLINE_INDEL_SCORING_DEVELOPMENT_FEATURES::AD2_NORM, (_isri.n_q30_alt_reads * q30DepthFactor));

    }
}


static void add_cigar_to_ploidy(
    const ALIGNPATH::path_t& apath,
    std::vector<unsigned>& ploidy)
{
    using namespace ALIGNPATH;
    int offset(-1);
    for (const auto& ps : apath)
    {
        if (is_segment_align_match(ps.type))
        {
            for (unsigned j(0); j<ps.length; ++j)
            {
                if (offset>=0) ploidy[offset]++;
                offset++;
            }
        }
        else if (ps.type==DELETE)
        {
            offset+=ps.length;
        }
    }
}


void
GermlineDiploidIndelCallInfo::
add_overlap(
    const reference_contig_segment& ref,
    GermlineDiploidIndelCallInfo& overlap)
{
    assert(_calls.size() == 1);
    assert(overlap._calls.size() == 1);

    auto& call(first());
    auto& overlap_call(overlap.first());

    // there's going to be 1 (possibly empty) fill range in front of one haplotype
    // and one possibly empty fill range on the back of one haplotype
    std::string leading_seq,trailing_seq;
    auto indel_end_pos=std::max(overlap_call._ik.right_pos(),call._ik.right_pos());

    const pos_t indel_begin_pos(pos-1);

    // add shared information (to the first indel only)
    // TODO: Reevaluate this. Since we're merging, we can do this on output
    // make extended vcf ref seq:
    std::string tmp;
    ref.get_substring(indel_begin_pos,(indel_end_pos-indel_begin_pos),tmp);
    call._iri.vcf_ref_seq = tmp;

    ploidy.resize(indel_end_pos-pos,0);

    auto munge_indel = [&] (GermlineDiploidIndelCallInfo& ii)
    {
        auto& this_call(ii.first());
        // extend leading sequence start back 1 for vcf compat, and end back 1 to concat with vcf_indel_seq
        ref.get_substring(indel_begin_pos,(ii.pos-indel_begin_pos)-1,leading_seq);
        const unsigned trail_len(indel_end_pos-this_call._ik.right_pos());
        ref.get_substring(indel_end_pos-trail_len,trail_len,trailing_seq);

        this_call._iri.vcf_indel_seq = leading_seq + this_call._iri.vcf_indel_seq + trailing_seq;

        this_call.set_hap_cigar(leading_seq.size()+1,
                                trailing_seq.size());

        // add to the ploidy object:
        add_cigar_to_ploidy(this_call.cigar,ploidy);
    };
    munge_indel(*this);
    munge_indel(overlap);

    // we only combine pairs of simple het indels on different haplotpyes, so this assertion must hold:
    for (const unsigned pl : ploidy)
    {
        assert(pl<2);
    }

    //reduce qual and gt to the lowest of the set:
    call._dindel.indel_qphred = std::min(call._dindel.indel_qphred, overlap.first()._dindel.indel_qphred);
    call._dindel.max_gt_qphred = std::min(call._dindel.max_gt_qphred, overlap.first()._dindel.max_gt_qphred);


    // combine filter flags
    call.filters |= overlap.first().filters;
    // combine QScores. Since the "unset" value is -1, this complex logic is necessary
    if (call.empiricalVariantScore <0)
        call.empiricalVariantScore = overlap.first().empiricalVariantScore;
    else if (overlap.first().empiricalVariantScore >= 0)
        call.empiricalVariantScore = std::min(call.empiricalVariantScore, overlap.first().empiricalVariantScore);
    call.gqx = std::min(call.gqx, overlap.first().gqx);
    call.gq = std::min(call.gq, overlap.first().gq);

    _calls.push_back(overlap.first());
}



void
GermlineDiploidIndelCallInfo::
getPloidyError(
    const unsigned offset) const
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: get_ploidy offset '" << offset << "' exceeds ploidy region size '" << ploidy.size() << "'\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



void
GermlineDiploidSiteCallInfo::
computeEmpiricalScoringFeatures(
    const bool isComputeDevelopmentFeatures,
    const double chromDepth,
    GermlineDiploidSiteSimpleGenotypeInfo& smod2) const
{
    const double chromDepthFactor(1./chromDepth);

    // get the alt base id (choose second in case of an alt het....)
    unsigned altBase(N_BASE);
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==dgt.ref_gt) continue;
        if (DIGT::expect2(b,smod.max_gt))
        {
            altBase=b;
        }
    }
    assert(altBase!=N_BASE);

    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::QUAL, (dgt.genome.snp_qphred * chromDepthFactor));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::F_GQX, (smod.gqx * chromDepthFactor));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::F_GQ, (smod.gq * chromDepthFactor));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVSB, (smod.strand_bias));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::I_SNVHPOL, (hpol));

    //we need to handle the scaling of DP better for high depth cases
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::F_DP, (n_used_calls * chromDepthFactor));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::F_DPF, (n_unused_calls * chromDepthFactor));

    const double r0 = alleleObservationCounts(dgt.ref_gt);
    const double r1 = alleleObservationCounts(altBase);
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::AD0, (r0 * chromDepthFactor));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::AD1, (r1 * chromDepthFactor));


    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::I_MQ, (mapqRMS));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::I_ReadPosRankSum, (ReadPosRankSum));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::I_BaseQRankSum, (BaseQRankSum));
    smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::I_MQRankSum, (MQRankSum));

    // allele bias metrics
    {
        const double allelebiaslower  = cdf(boost::math::binomial(r0+r1,0.5),r0);
        const double allelebiasupper  = cdf(boost::math::binomial(r0+r1,0.5),r1);

        // +1e-30 to avoid log(0) in extreme cases
        smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::ABlower, (-log(allelebiaslower+1.e-30)));
        smod2.features.set(GERMLINE_SNV_SCORING_FEATURES::AB, (-log(std::min(1.,2.*std::min(allelebiaslower,allelebiasupper))+1.e-30)));
    }

    //
    // compute any experimental features not currently used in production
    //
    if (isComputeDevelopmentFeatures)
    {
        //The average baseQ of the position of alt allele
        smod2.developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawBaseQ, (avgBaseQ));

        //the average position value within a read of alt allele
        smod2.developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::I_RawPos, (rawPos));

        smod2.developmentFeatures.set(GERMLINE_SNV_SCORING_DEVELOPMENT_FEATURES::mapqZeroFraction, (safeFrac(mapqZeroCount,mapqCount)));
    }
}



std::ostream&
operator<<(std::ostream& os,
           const GermlineVariantSimpleGenotypeInfo& shmod)
{
    os << "gqx: " << shmod.gqx
       << " gq: " << shmod.gq;
    if (typeid(shmod) == typeid(GermlineDiploidIndelSimpleGenotypeInfo))
    {
        auto imod = dynamic_cast<const GermlineDiploidIndelSimpleGenotypeInfo&>(shmod);

        os << " max_gt: " << DIGT::label(imod.max_gt);
    }

    os << " filters: ";
    shmod.write_filters(os);

    return os;
}

std::ostream&
operator<<(std::ostream& os,
           const GermlineDiploidSiteSimpleGenotypeInfo& smod)
{
    os << static_cast<GermlineVariantSimpleGenotypeInfo>(smod) << '\n';

    os << "is_unknown: " << smod.is_unknown;
    os << " is_covered: " << smod.is_covered;
    os << " is_used_coverage: " << smod.is_used_covered;
    os << " is_zero_ploidy: " << smod.is_zero_ploidy;

    if (smod.modified_gt != MODIFIED_SITE_GT::NONE)
    {
        os << " modgt: " << MODIFIED_SITE_GT::get_label(smod.modified_gt);
    }

    return os;
}



std::ostream&
operator<<(std::ostream& os,
           const GermlineDiploidSiteCallInfo& si)
{
    os << "pos: " << (si.pos+1) << " " << si.get_gt();
    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineIndelSimpleGenotypeInfo& shi)
{
    os << static_cast<GermlineVariantSimpleGenotypeInfo>(shi) << '\n';

    os << "indel_key: " << shi._ik << "\n";
    //os << "indel_data: " << shi._id << "\n";
    os << "indel_report_info: " << shi._iri << "\n";
    os << "indel_sample_info: " << shi._isri << "\n";
    os << "cigar: " << shi.cigar << "\n";

    return os;
}



std::ostream&
operator<<(
    std::ostream& os,
    const GermlineDiploidIndelSimpleGenotypeInfo& dic)
{
    os << static_cast<GermlineIndelSimpleGenotypeInfo>(dic) << '\n';

    dic._dindel.dump(os);

    os << "EVS: " << dic.empiricalVariantScore << " max_gt: " << dic.max_gt << "\n";

    return os;
}


void
GermlineDiploidIndelCallInfo::
dump(std::ostream& os) const
{
    os << "digt_indel_info\n";
    os << "nCalls: " << _calls.size() << " isOverlap: " << _is_overlap << "\n";
    os << "ploidy: ";
    for (const unsigned pl : ploidy)
    {
        os << " " << pl;
    }
    os << "\n";
    os << "Calls:\n";
    for (const auto& cl : _calls)
    {
        os << cl << "\n";
    }
}


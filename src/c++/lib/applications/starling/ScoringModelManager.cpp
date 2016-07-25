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

/*
 * \author Morten Kallberg
 */

#include "ScoringModelManager.hh"


//#define DEBUG_CAL

#ifdef DEBUG_CAL
#include "blt_util/log.hh"
#endif



ScoringModelManager::
ScoringModelManager(
    const starling_options& opt,
    const gvcf_deriv_options& gvcfDerivedOptions)
    : _opt(opt.gvcf),
      _dopt(gvcfDerivedOptions),
      _isReportEVSFeatures(opt.isReportEVSFeatures),
      _isRNA(opt.isRNA)
{
    const SCORING_CALL_TYPE::index_t callType(opt.isRNA ? SCORING_CALL_TYPE::RNA : SCORING_CALL_TYPE::GERMLINE);

    if (not opt.snv_scoring_model_filename.empty())
    {
        _snvScoringModelPtr.reset(
            new VariantScoringModelServer(
                _dopt.snvFeatureSet.getFeatureMap(),
                opt.snv_scoring_model_filename,
                callType,
                SCORING_VARIANT_TYPE::SNV)
        );
    }

    if (not opt.indel_scoring_model_filename.empty())
    {
        _indelScoringModelPtr.reset(
            new VariantScoringModelServer(
                _dopt.indelFeatureSet.getFeatureMap(),
                opt.indel_scoring_model_filename,
                callType,
                SCORING_VARIANT_TYPE::INDEL)
        );
    }
}



void
ScoringModelManager::
classify_site(
    GermlineDiploidSiteLocusInfo& si,
    GermlineDiploidSiteAlleleInfo& smod) const
{
    if (si.dgt.is_snp && _isReportEVSFeatures)
    {
        // when reporting is turned on, we need to compute EVS features
        // for any usable variant regardless of EVS model type:
        const bool isUniformDepthExpected(_dopt.is_max_depth());
        si.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, _isReportEVSFeatures, _dopt.norm_depth);
    }

    //si.smod.filters.reset(); // make sure no filters have been applied prior
    if (si.dgt.is_snp && isEVSSiteModel())
    {
        if (si.EVSFeatures.empty())
        {
            static const bool isComputeDevelopmentFeatures(false);
            const bool isUniformDepthExpected(_dopt.is_max_depth());
            si.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, isComputeDevelopmentFeatures, _dopt.norm_depth);
        }
        si.empiricalVariantScore = error_prob_to_qphred(_snvScoringModelPtr->scoreVariant(si.EVSFeatures.getAll()));

        static const int maxEmpiricalVariantScore(60);
        si.empiricalVariantScore = std::min(si.empiricalVariantScore, maxEmpiricalVariantScore);

        if (si.empiricalVariantScore < snvEVSThreshold())
        {
            si.filters.set(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
        }
    }
    else
    {
        // don't know what to do with this site, throw it to the old default filters
        default_classify_site(si, smod);
    }
}



bool
ScoringModelManager::
checkIsVariantUsableInEVSModel(const GermlineDiploidIndelLocusInfo& ii) const
{
    const auto& call(ii.getFirstAltAllele());
    return ((call._indelReportInfo.it == SimplifiedIndelReportType::INSERT ||
             call._indelReportInfo.it == SimplifiedIndelReportType::DELETE ||
             call._indelReportInfo.it == SimplifiedIndelReportType::SWAP) &&
            (call._dindel.max_gt != STAR_DIINDEL::NOINDEL) ); // empirical scoring does not handle homref sites
}



void
ScoringModelManager::
set_indel_modifiers(
    const GermlineDiploidIndelLocusInfo& ii,
    GermlineDiploidIndelAlleleInfo& call) const
{
    const auto& dindel(ii.getFirstAltAllele()._dindel);
    /// max_gt != max_gt_poly indicates we're in a boundary zone between variant and hom-ref call
    if (dindel.max_gt != dindel.max_gt_poly)
    {
        call.gqx=0;
    }
    else
    {
        call.gqx=std::min(dindel.max_gt_poly_qphred,dindel.max_gt_qphred);
    }
    call.max_gt=dindel.max_gt_poly;
    call.gq=dindel.max_gt_poly_qphred;
}



void
ScoringModelManager::
classify_indel_impl(
    const bool isVariantUsableInEVSModel,
    GermlineDiploidIndelLocusInfo& ii,
    GermlineDiploidIndelAlleleInfo& call) const
{
    set_indel_modifiers(ii, call);

    if (isVariantUsableInEVSModel && _isReportEVSFeatures)
    {
        // when reporting is turned on, we need to compute EVS features
        // for any usable variant regardless of EVS model type:
        const bool isUniformDepthExpected(_dopt.is_max_depth());
        ii.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, _isReportEVSFeatures, _dopt.norm_depth);
    }

    if (isEVSIndelModel() && isVariantUsableInEVSModel)
    {
        if (ii.features.empty())
        {
            static const bool isComputeDevelopmentFeatures(false);
            const bool isUniformDepthExpected(_dopt.is_max_depth());
            ii.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, isComputeDevelopmentFeatures, _dopt.norm_depth);
        }
        ii.empiricalVariantScore = error_prob_to_qphred(_indelScoringModelPtr->scoreVariant(ii.features.getAll()));

        static const int maxEmpiricalVariantScore(60);
        ii.empiricalVariantScore = std::min(ii.empiricalVariantScore, maxEmpiricalVariantScore);

        if (ii.empiricalVariantScore < indelEVSThreshold())
        {
            ii.filters.set(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
        }
    }
    else
    {
        default_classify_indel(ii, call);
    }
}



void
ScoringModelManager::
classify_indel(
    GermlineDiploidIndelLocusInfo& ii,
    GermlineDiploidIndelAlleleInfo& call) const
{
    classify_indel_impl(checkIsVariantUsableInEVSModel(ii),ii,call);
}



void
ScoringModelManager::
classify_indels(
    std::vector<std::unique_ptr<GermlineDiploidIndelLocusInfo>>& indels) const
{
    bool isVariantUsableInEVSModel = true;
    for (const auto& indel : indels)
    {
        if (! isVariantUsableInEVSModel) break;
        isVariantUsableInEVSModel = checkIsVariantUsableInEVSModel(*indel);
    }

    for (auto& indel : indels)
    {
        GermlineDiploidIndelLocusInfo& ii(*indel);
        classify_indel_impl(isVariantUsableInEVSModel,ii, ii.getFirstAltAllele());
    }
}



void
ScoringModelManager::
default_classify_site(
    GermlineSiteLocusInfo& si,
    const GermlineVariantAlleleInfo& call) const
{
    if (_opt.is_min_gqx)
    {
        if (call.gqx<_opt.min_gqx) si.filters.set(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    }
    if (_dopt.is_max_depth())
    {
        if ((si.n_used_calls+si.n_unused_calls) > _dopt.max_depth)
            si.filters.set(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }
    // high DPFratio filter
    if (_opt.is_max_base_filt)
    {
        const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
        if (total_calls>0)
        {
            const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
            if (filt>_opt.max_base_filt) si.filters.set(GERMLINE_VARIANT_VCF_FILTERS::HighBaseFilt);
        }
    }
    if (si.is_snp())
    {
        if (_opt.is_max_snv_sb)
        {
            if (call.strand_bias>_opt.max_snv_sb) si.filters.set(GERMLINE_VARIANT_VCF_FILTERS::HighSNVSB);
        }
        if (_opt.is_max_snv_hpol)
        {
            if (static_cast<int>(si.hpol)>_opt.max_snv_hpol) si.filters.set(GERMLINE_VARIANT_VCF_FILTERS::HighSNVHPOL);
        }
    }
}



void
ScoringModelManager::
default_classify_indel(
    GermlineIndelLocusInfo& ii,
    const GermlineIndelAlleleInfo& call) const
{
    if (this->_opt.is_min_gqx)
    {
        if (call.gqx<_opt.min_gqx) ii.filters.set(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    }

    if (this->_dopt.is_max_depth())
    {
        if (call._indelSampleReportInfo.tier1Depth > this->_dopt.max_depth)
            ii.filters.set(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }

    if (_opt.is_max_ref_rep())
    {
        const auto& iri(call._indelReportInfo);
        if (iri.is_repeat_unit())
        {
            if ((iri.repeat_unit.size() <= 2) &&
                (static_cast<int>(iri.ref_repeat_count) > _opt.max_ref_rep))
            {
                ii.filters.set(GERMLINE_VARIANT_VCF_FILTERS::HighRefRep);
            }
        }
    }
}

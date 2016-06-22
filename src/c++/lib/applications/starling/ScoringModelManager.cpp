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
    if (not opt.germline_variant_scoring_models_filename.empty())
    {
        /// this if is TEMPORARY until legacy support is taken out:
        if (opt.isRNA)
        {
            const SCORING_CALL_TYPE::index_t callType(opt.isRNA ? SCORING_CALL_TYPE::RNA : SCORING_CALL_TYPE::GERMLINE);

            assert (opt.germline_variant_scoring_model_name.empty());

            _snvScoringModelPtr.reset(
                new VariantScoringModelServer(
                    _dopt.snvFeatureSet.getFeatureMap(),
                    opt.germline_variant_scoring_models_filename,
                    callType,
                    SCORING_VARIANT_TYPE::SNV)
            );

            /// TEMPORARY suspend indel RF file load until a valid model file exists
#if 0
            _indelScoringModelPtr.reset(
                new VariantScoringModelServer(
                    _dopt.indelFeatureSet.getFeatureMap(),
                    opt.germline_variant_scoring_models_filename,
                    callType,
                    SCORING_VARIANT_TYPE::INDEL)
            );
#endif
        }
        else
        {
            assert (not opt.germline_variant_scoring_model_name.empty());
            _legacyModelPtr.reset(
                new LogisticAndRuleScoringModels(_dopt, opt.germline_variant_scoring_models_filename,
                                                 opt.germline_variant_scoring_model_name));
        }
    }
}



int
ScoringModelManager::
get_case_cutoff(
    const LEGACY_CALIBRATION_MODEL::var_case my_case) const
{
    if (not isLegacyModel()) return 0;
    return getLegacyModel().get_var_threshold(my_case);
}



bool ScoringModelManager::isLegacyLogisticEVSModel() const
{
    if (not isLegacyModel()) return false;
    return getLegacyModel().is_logistic_model();
}


void
ScoringModelManager::
classify_site(
    const GermlineDiploidSiteCallInfo& si,
    GermlineDiploidSiteSimpleGenotypeInfo& smod) const
{
    if (si.dgt.is_snp && _isReportEVSFeatures)
    {
        // when reporting is turned on, we need to compute EVS features
        // for any usable variant regardless of EVS model type:
        const bool isUniformDepthExpected(_dopt.is_max_depth());
        si.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, _isReportEVSFeatures, _dopt.norm_depth, smod);
    }

    //si.smod.filters.reset(); // make sure no filters have been applied prior
    if (si.dgt.is_snp && isEVSSiteModel())
    {
        if (isLegacyModel())
        {
            getLegacyModel().score_site_instance(si, smod);
        }
        else
        {
            if (smod.features.empty())
            {
                static const bool isComputeDevelopmentFeatures(false);
                const bool isUniformDepthExpected(_dopt.is_max_depth());
                si.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, isComputeDevelopmentFeatures, _dopt.norm_depth, smod);
            }
            smod.empiricalVariantScore = static_cast<int>(_snvScoringModelPtr->scoreVariant(smod.features.getAll()));

            static const int maxEmpiricalVariantScore(60);
            smod.empiricalVariantScore = std::min(smod.empiricalVariantScore, maxEmpiricalVariantScore);

            if (smod.empiricalVariantScore < snvEVSThreshold())
            {
                smod.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
            }
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
checkIsVariantUsableInEVSModel(const GermlineDiploidIndelCallInfo& ii) const
{
    const auto& call(ii.first());
    return ((call._indelReportInfo.it == INDEL::INSERT || call._indelReportInfo.it == INDEL::DELETE) &&
            (call._dindel.max_gt != STAR_DIINDEL::NOINDEL) ); // empirical scoring does not handle homref sites properly
}

void
ScoringModelManager::
set_indel_modifiers(
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    const auto& dindel(ii.first()._dindel);
    if ((dindel.max_gt != dindel.max_gt_poly) || dindel.is_zero_coverage)
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
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    set_indel_modifiers(ii, call);

    if (isVariantUsableInEVSModel && _isReportEVSFeatures)
    {
        // when reporting is turned on, we need to compute EVS features
        // for any usable variant regardless of EVS model type:
        const bool isUniformDepthExpected(_dopt.is_max_depth());
        call.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, _isReportEVSFeatures, _dopt.norm_depth, ii.is_hetalt());
    }

    if (isEVSIndelModel() && isVariantUsableInEVSModel)
    {
        if (isLegacyModel())
        {
            getLegacyModel().score_indel_instance(ii, call);
        }
        else if (_isRNA)
        {
            // TEMPORARY HACK TO DISABLE RNA RF for indels
            default_classify_indel(call);
        }
        else
        {
            if (call.features.empty())
            {
                static const bool isComputeDevelopmentFeatures(false);
                const bool isUniformDepthExpected(_dopt.is_max_depth());
                call.computeEmpiricalScoringFeatures(_isRNA, isUniformDepthExpected, isComputeDevelopmentFeatures, _dopt.norm_depth, ii.is_hetalt());
            }
            call.empiricalVariantScore = static_cast<int>(_indelScoringModelPtr->scoreVariant(call.features.getAll()));

            static const int maxEmpiricalVariantScore(60);
            call.empiricalVariantScore = std::min(call.empiricalVariantScore, maxEmpiricalVariantScore);

            if (call.empiricalVariantScore < indelEVSThreshold())
            {
                call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
            }
        }
    }
    else
    {
        default_classify_indel(call);
    }
}


void
ScoringModelManager::
classify_indel(
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    classify_indel_impl(checkIsVariantUsableInEVSModel(ii),ii,call);
}



void
ScoringModelManager::
classify_indels(
    std::vector<std::unique_ptr<GermlineDiploidIndelCallInfo>>& indels) const
{
    bool isVariantUsableInEVSModel = true;
    for (const auto& indel : indels)
    {
        if (! isVariantUsableInEVSModel) break;
        isVariantUsableInEVSModel = checkIsVariantUsableInEVSModel(*indel);
    }

    for (auto& indel : indels)
    {
        GermlineDiploidIndelCallInfo& ii(*indel);
        classify_indel_impl(isVariantUsableInEVSModel,ii,ii.first());
    }
}




void
ScoringModelManager::
default_classify_site(
    const GermlineSiteCallInfo& si,
    GermlineVariantSimpleGenotypeInfo& call) const
{
    if (_opt.is_min_gqx)
    {
        if (call.gqx<_opt.min_gqx) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    }
    if (_dopt.is_max_depth())
    {
        if ((si.n_used_calls+si.n_unused_calls) > _dopt.max_depth) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }
    // high DPFratio filter
    if (_opt.is_max_base_filt)
    {
        const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
        if (total_calls>0)
        {
            const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
            if (filt>_opt.max_base_filt) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighBaseFilt);
        }
    }
    if (si.is_snp())
    {
        if (_opt.is_max_snv_sb)
        {
            if (call.strand_bias>_opt.max_snv_sb) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighSNVSB);
        }
        if (_opt.is_max_snv_hpol)
        {
            if (static_cast<int>(si.hpol)>_opt.max_snv_hpol) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighSNVHPOL);
        }
    }
}



// default rules based indel model
void
ScoringModelManager::
default_classify_indel(
    GermlineIndelSimpleGenotypeInfo& call) const
{
    if (this->_opt.is_min_gqx)
    {
        if (call.gqx<_opt.min_gqx) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    }

    if (this->_dopt.is_max_depth())
    {
        if (call._indelSampleReportInfo.tier1Depth > this->_dopt.max_depth) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }

    if (this->_opt.is_max_ref_rep)
    {
        const auto& iri(call._indelReportInfo);
        if (iri.is_repeat_unit())
        {
            if ((iri.repeat_unit.size() <= 2) &&
                (static_cast<int>(iri.ref_repeat_count) > this->_opt.max_ref_rep))
            {
                call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighRefRep);
            }
        }
    }
}

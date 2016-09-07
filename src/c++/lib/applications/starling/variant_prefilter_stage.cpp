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
 *  Created on: Jun 4, 2015
 *      Author: jduddy
 */

#include "variant_prefilter_stage.hh"

#include "ScoringModelManager.hh"

variant_prefilter_stage::variant_prefilter_stage(const ScoringModelManager& model, std::shared_ptr<variant_pipe_stage_base> destination)
    : variant_pipe_stage_base(destination)
    , _model(model)
{
}

static
void
set_site_gt(
    const diploid_genotype::result_set& rs,
    LocusSampleInfo& sampleInfo,
    GermlineSiteSampleInfo& siteSampleInfo)
{
    siteSampleInfo.max_gt=rs.max_gt;
    sampleInfo.gqx=rs.max_gt_qphred;
}


void
variant_prefilter_stage::
add_site_modifiers(
    GermlineDiploidSiteLocusInfo& locus,
    const ScoringModelManager& model)
{
    auto& allele(locus.allele);
    allele.clear();

    /// TODO STREL-125 generalize to multi-sample
    const unsigned sampleIndex(0);
    LocusSampleInfo& sampleInfo(locus.getSample(sampleIndex));
    auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

    /// TODO STREL-125 find a way to restore strand bias feature
    // allele.strandBias=locus.dgt.strand_bias;

    if     (locus.isRefUnknown())
    {
        sampleInfo.gqx=0;
        sampleInfo.genotypeQualityPolymorphic=0;
        siteSampleInfo.max_gt=0;
    }
    else if (siteSampleInfo.dgt.genome.max_gt != siteSampleInfo.dgt.poly.max_gt)
    {
        sampleInfo.gqx=0;
        sampleInfo.genotypeQualityPolymorphic=siteSampleInfo.dgt.poly.max_gt_qphred;
        siteSampleInfo.max_gt=siteSampleInfo.dgt.poly.max_gt;
    }
    else
    {
        if (siteSampleInfo.dgt.genome.max_gt_qphred<siteSampleInfo.dgt.poly.max_gt_qphred)
        {
            set_site_gt(siteSampleInfo.dgt.genome, sampleInfo, siteSampleInfo);
        }
        else
        {
            set_site_gt(siteSampleInfo.dgt.poly, sampleInfo, siteSampleInfo);
        }
        sampleInfo.genotypeQualityPolymorphic=siteSampleInfo.dgt.poly.max_gt_qphred;
    }

    model.classify_site(locus);
}



void
variant_prefilter_stage::
process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    /// TODO STREL-125 generalize to N samples
    static const unsigned sampleIndex(0);

    const uint8_t refBaseId(base_to_id(locusPtr->ref));

    if (dynamic_cast<GermlineDiploidSiteLocusInfo*>(locusPtr.get()) != nullptr)
    {
        auto diploidLocusPtr(downcast<GermlineDiploidSiteLocusInfo>(std::move(locusPtr)));
        const auto& sampleInfo(diploidLocusPtr->getSample(sampleIndex));
        auto& siteSampleInfo(diploidLocusPtr->getSiteSample(sampleIndex));

        add_site_modifiers(*diploidLocusPtr, _model);
        if (sampleInfo.getPloidy().isHaploid())
        {
            if (siteSampleInfo.max_gt == refBaseId)
            {
                siteSampleInfo.modified_gt=MODIFIED_SITE_GT::ZERO;
            }
            else
            {
                siteSampleInfo.modified_gt=MODIFIED_SITE_GT::ONE;
            }
        }
        else if (sampleInfo.getPloidy().isNoploid())
        {
            if (! diploidLocusPtr->is_print_unknowngt(sampleIndex))
            {
                diploidLocusPtr->filters.set(GERMLINE_VARIANT_VCF_FILTERS::PloidyConflict);
            }
        }

        _sink->process(std::move(diploidLocusPtr));
    }
    else
    {
        for (auto& altAllele : locusPtr->getSiteAlleles())
        {
            _model.default_classify_site_locus(*locusPtr, altAllele);
        }

        _sink->process(std::move(locusPtr));
    }
}



void
variant_prefilter_stage::
process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
{
    // we can't handle breakends at all right now:
    for (const auto& altAllele : locusPtr->getIndelAlleles())
    {
        if (altAllele.indelKey.is_breakpoint()) return;
    }

    // add filter for all indels in 'no-ploid' regions:
    const unsigned sampleCount(locusPtr->getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        LocusSampleInfo& sampleInfo(locusPtr->getSample(sampleIndex));
        if (sampleInfo.isPloidyConflict())
        {
            sampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::PloidyConflict);
        }
    }

    // apply filtration/EVS model:
    if (dynamic_cast<GermlineContinuousIndelLocusInfo*>(locusPtr.get()) != nullptr)
    {
        _model.default_classify_indel_locus(*locusPtr);
    }
    else
    {
        _model.classify_indel(dynamic_cast<GermlineDiploidIndelLocusInfo&>(*locusPtr));
    }

    _sink->process(std::move(locusPtr));
}

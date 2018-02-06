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
/*
 *  Created on: Jun 4, 2015
 *      Author: jduddy
 */

#include "variant_prefilter_stage.hh"

#include "ScoringModelManager.hh"



void
variant_prefilter_stage::
applySharedLocusFilters(
    LocusInfo& locus) const
{
    // add filter for all sites in 'no-ploid' regions:
    const unsigned sampleCount(locus.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        LocusSampleInfo& sampleInfo(locus.getSample(sampleIndex));
        if (sampleInfo.isPloidyConflict())
        {
            sampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::PloidyConflict);
        }
    }
}



void
variant_prefilter_stage::
process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    applySharedLocusFilters(*locusPtr);

    // apply depth filter
    _model.applyDepthFilter(*locusPtr);

    // apply filtration/EVS model:
    if (dynamic_cast<GermlineContinuousSiteLocusInfo*>(locusPtr.get()) != nullptr)
    {
        _model.default_classify_site_locus(*locusPtr);
    }
    else
    {
        _model.classify_site(dynamic_cast<GermlineDiploidSiteLocusInfo&>(*locusPtr));
    }

    _sink->process(std::move(locusPtr));
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

    if (locusPtr->isNotGenotyped())
    {
        const unsigned sampleCount(locusPtr->getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            auto& sampleInfo(locusPtr->getSample(sampleIndex));
            sampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::NotGenotyped);
        }
    }
    else
    {
        applySharedLocusFilters(*locusPtr);

        // apply depth filter
        _model.applyDepthFilter(*locusPtr);

        // apply filtration/EVS model:
        if (dynamic_cast<GermlineContinuousIndelLocusInfo*>(locusPtr.get()) != nullptr)
        {
            _model.default_classify_indel_locus(*locusPtr);
        }
        else
        {
            _model.classify_indel(dynamic_cast<GermlineDiploidIndelLocusInfo&>(*locusPtr));
        }
    }

    _sink->process(std::move(locusPtr));
}

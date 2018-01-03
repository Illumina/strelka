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

/// \file
/// \author John Duddy

#include "VariantOverlapResolver.hh"
#include "ScoringModelManager.hh"
#include "blt_util/log.hh"

#include <iostream>

//#define DEBUG_GVCF



void
VariantOverlapResolver::
process(std::unique_ptr<GermlineSiteLocusInfo> siteLocusPtr)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: Begin " << __FUNCTION__ << " siteLocus\n";
    log_os << "CHIRP: " << __FUNCTION__ << " pos/bufferedVariantIndelRange: " << siteLocusPtr->pos << "/" << _bufferedVariantIndelRange << "\n";
#endif

    if (! _variantIndelBuffer.empty())
    {
        if (siteLocusPtr->pos >= _bufferedVariantIndelRange.end_pos())
        {
            // Resolve any current or previous indels before forwarding this site down the pipeline:
            processOverlappingVariants();
        }
        else
        {
            // Buffer this site in the context of an overlapping variant indel record
            //
            // the site must be buffered even if it occurs before the variant indel range so that nonvariant
            // indels and sites can be ordered correctly
            //
            // \TODO STREL-617
            // Currently need to downcast to diploid type in this object so that filters can be
            // rerun on sites overlapped by deletions. Would be better to hide this detail from the overlapper, and
            // have a general "refilter" method shared by this component and variant_prefilter stage, which hides this
            // detail from both
            //
            _siteBuffer.push_back(downcast<GermlineDiploidSiteLocusInfo>(std::move(siteLocusPtr)));
            return;
        }
    }

    _sink->process(std::move(siteLocusPtr));
}



void
VariantOverlapResolver::
process(std::unique_ptr<GermlineIndelLocusInfo> indelLocusPtr)
{
    const bool isVariantLocus(indelLocusPtr->isVariantLocus());
    const bool isForcedOutputLocus(indelLocusPtr->isAnyForcedOutputAtLocus());

    // only handle variant or forced indel loci
    if (! (isVariantLocus || isForcedOutputLocus)) return;

#ifdef DEBUG_GVCF
    log_os << "CHIRP: Begin " << __FUNCTION__ << " indelLocus\n";
    log_os << "CHIRP: " << __FUNCTION__ << " pos/indel_range: " << indelLocusPtr->pos << "/" << _bufferedVariantIndelRange << "\n";
    log_os << "CHIRP: " << __FUNCTION__ << " isVariant/isForcedOutput: " << isVariantLocus << "/" << isForcedOutputLocus << "\n";
#endif

    if (! _variantIndelBuffer.empty())
    {
        // the test is greater than, and not greater than or equals, because two adjacent indels are considered
        // as interacting/interfering
        if (indelLocusPtr->pos > _bufferedVariantIndelRange.end_pos())
        {
            processOverlappingVariants();
        }
    }

    if (isVariantLocus)
    {
        if (_variantIndelBuffer.empty())
        {
            _bufferedVariantIndelRange = known_pos_range2(indelLocusPtr->range());
        }
        else
        {
            _bufferedVariantIndelRange.merge_range(indelLocusPtr->range());
        }
        _variantIndelBuffer.push_back(std::move(indelLocusPtr));
    }
    else
    {
        if (_variantIndelBuffer.empty())
        {
            _sink->process(std::move(indelLocusPtr));
        }
        else
        {
            // If variant indels are buffered, we already know that this indel intersects or is upstream
            // of the buffered variant, and therefore must be buffered to resolve possible indel/site
            // reordering requirements.
            //
            _nonvariantIndelBuffer.push_back(std::move(indelLocusPtr));
        }
    }
}



template <typename T>
void
dumpLocusBuffer(
    const char* locusTypeLabel,
    const std::vector<std::unique_ptr<T>>& locusBuffer,
    std::ostream& os)
{
    // This function may need to deal with object data in an intermediate state when some locus
    // pointers have already been sunk (especially if called while building an exception report)
    //
    const unsigned locusCount(locusBuffer.size());
    os << locusTypeLabel << " count: (" << locusCount << ")\n";
    for (unsigned locusIndex(0); locusIndex < locusCount; ++locusIndex)
    {
        os << locusTypeLabel << locusIndex << " ";
        const auto& locus(locusBuffer[locusIndex]);
        if (locus)
        {
            os << *locus;
        }
        else
        {
            os << "ALREADY RELEASED";
        }
        os << '\n';
    }
}



void
VariantOverlapResolver::
dump(std::ostream& os) const
{
    os << "VariantOverlapResolver:"
       << " VariantIndelRange: " << _bufferedVariantIndelRange << "\n";
    dumpLocusBuffer("Site", _siteBuffer, os);
    dumpLocusBuffer("VariantIndel", _variantIndelBuffer, os);
    dumpLocusBuffer("NonVariantIndel", _nonvariantIndelBuffer, os);
}



void
VariantOverlapResolver::
processOverlappingVariants()
{
    try
    {
        processOverlappingVariantsImplementation();
    }
    catch (...)
    {
        log_os << "ERROR: exception caught in " << __FUNCTION__ << "\n";
        dump(log_os);

        // The buffers need to be cleared in case this object is in an unstable state, otherwise
        // the flush() call into indel_buffer could trigger another exception or segfault
        clearBuffers();

        throw;
    }
}


namespace VARQUEUE
{
enum index_t
{
    NONE,
    VARIANT_INDEL,
    NONVARIANT_INDEL,
    SITE
};
}


/// Implements expected ordering for all buffered variants in the object
///
/// Note this doesn't really generalize or tidy up the (implicit) indel/site priority queue, but
/// it at least dumps the ugliness into one place.
static
VARQUEUE::index_t
nextVariantType(
    const std::vector<std::unique_ptr<GermlineIndelLocusInfo>>& variantIndelBuffer,
    const std::vector<std::unique_ptr<GermlineIndelLocusInfo>>& nonvariantIndelBuffer,
    const std::vector<std::unique_ptr<GermlineDiploidSiteLocusInfo>>& siteBuffer,
    const unsigned variantIndelIndex,
    const unsigned nonvariantIndelIndex,
    const unsigned siteIndex)
{
    const bool areVariantIndelsAvailable(variantIndelIndex < variantIndelBuffer.size());
    const bool areNonvariantIndelsAvailable(nonvariantIndelIndex < nonvariantIndelBuffer.size());
    const bool areSitesAvailable(siteIndex < siteBuffer.size());

    if ((!areVariantIndelsAvailable) && (!areNonvariantIndelsAvailable) && (!areSitesAvailable))
    {
        return VARQUEUE::NONE;
    }

    // For A=variantIndels, B=nonvariantIndels, C=sites...
    const bool AlessB(areVariantIndelsAvailable && ((! areNonvariantIndelsAvailable) || (variantIndelBuffer[variantIndelIndex]->pos <= nonvariantIndelBuffer[nonvariantIndelIndex]->pos)));
    const bool AlessC(areVariantIndelsAvailable && ((! areSitesAvailable) || (variantIndelBuffer[variantIndelIndex]->pos <= siteBuffer[siteIndex]->pos)));
    const bool BlessC(areNonvariantIndelsAvailable && ((! areSitesAvailable) || (nonvariantIndelBuffer[nonvariantIndelIndex]->pos <= siteBuffer[siteIndex]->pos)));

    if (AlessB)
    {
        if (AlessC)
        {
            return VARQUEUE::VARIANT_INDEL;
        }
        else
        {
            return VARQUEUE::SITE;
        }
    }
    else
    {
        if (BlessC)
        {
            return VARQUEUE::NONVARIANT_INDEL;
        }
        else
        {
            return VARQUEUE::SITE;
        }
    }
}



void
VariantOverlapResolver::
processOverlappingVariantsImplementation()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: Begin " << __FUNCTION__ << "\n";
    dump(log_os);
#endif

    // If no indels are buffered, then there is nothing to process
    if (_variantIndelBuffer.empty() && _nonvariantIndelBuffer.empty())
    {
        assert(_siteBuffer.empty());
        return;
    }

    bool isVariantConflictFound(false);

    // Process conflicting indel loci (these should be rare)
    if (_variantIndelBuffer.size() > 1)
    {
        // sort the indel buffer by begin pos:
        std::stable_sort(std::begin(_variantIndelBuffer), std::end(_variantIndelBuffer),
                         [](const indel_ptr_t& lhs, const indel_ptr_t& rhs) -> bool
        {
            return (lhs->range().begin_pos() < rhs->range().begin_pos());
        });

        // mark the whole region as conflicting
        annotateVariantIndelRecordsAsConflicting();
        isVariantConflictFound = true;
    }

    // process sites to be consistent with the first overlapping variant indel:
    //

    // check that if anything is in the site buffer, we have at least one variant indel:
    // (this (1) must be true (2) guards the _variantIndelBuffer.front() access below)
    assert(_siteBuffer.empty() || (! _variantIndelBuffer.empty()));

    for (auto& siteLocusPtr : _siteBuffer)
    {
#ifdef DEBUG_GVCF
        log_os << "CHIRP: " << __FUNCTION__ << " indel overlapping site at position: " << siteLocusPtr->pos << "\n";
#endif
        if (siteLocusPtr->pos < _bufferedVariantIndelRange.begin_pos()) continue;
        modifySiteOverlappingVariantIndel(*(_variantIndelBuffer.front()), *siteLocusPtr, _scoringModels);
    }

    //
    // Order all buffered indel and site record output according to VCF formatting rules and strelka gvcf conventions
    //

    unsigned variantIndelIndex(0);
    unsigned nonvariantIndelIndex(0);
    unsigned siteIndex(0);

    while (true)
    {
        const VARQUEUE::index_t nextvar = nextVariantType(_variantIndelBuffer, _nonvariantIndelBuffer, _siteBuffer,
                                                          variantIndelIndex, nonvariantIndelIndex, siteIndex);

        if      (nextvar == VARQUEUE::NONE)
        {
            break;
        }
        else if (nextvar == VARQUEUE::VARIANT_INDEL)
        {
            _sink->process(std::move(_variantIndelBuffer[variantIndelIndex]));
            if (isVariantConflictFound)
            {
                // emit each conflict record
                variantIndelIndex++;
            }
            else
            {
                // just emit the overlapped or single non-conflict record
                variantIndelIndex=_variantIndelBuffer.size();
            }
        }
        else if (nextvar == VARQUEUE::NONVARIANT_INDEL)
        {
            _sink->process(std::move(_nonvariantIndelBuffer[nonvariantIndelIndex]));
            nonvariantIndelIndex++;
        }
        else if (nextvar == VARQUEUE::SITE)
        {
            _sink->process(std::move(_siteBuffer[siteIndex]));
            siteIndex++;
        }
        else
        {
            assert(false && "unexpected varqueue type");
        }
    }

    clearBuffers();
}



/// Annotate locus with the indel conflict filter
static
void
annotateVariantAsIndelConflict(LocusInfo& locus)
{
    locus.filters.set(GERMLINE_VARIANT_VCF_FILTERS::IndelConflict);
}



void
VariantOverlapResolver::
modifySiteOverlappingVariantIndel(
    const GermlineIndelLocusInfo& indelLocus,
    GermlineDiploidSiteLocusInfo& siteLocus,
    const ScoringModelManager& model)
{
    if (indelLocus.filters.test(GERMLINE_VARIANT_VCF_FILTERS::IndelConflict))
    {
        annotateVariantAsIndelConflict(siteLocus);
    }
    else
    {
        modifySiteOverlappingNonconflictingVariantIndel(indelLocus, siteLocus, model);
    }
}



void
VariantOverlapResolver::
modifySiteOverlappingNonconflictingVariantIndel(
    const GermlineIndelLocusInfo& indelLocus,
    GermlineDiploidSiteLocusInfo& siteLocus,
    const ScoringModelManager& model)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: Begin " << __FUNCTION__ << " siteLocus before: " << siteLocus << "\n";
    log_os << "CHIRP: Begin " << __FUNCTION__ << " indelLocus: " << indelLocus << "\n";
#endif

    // site must be positioned within the indel range
    assert(indelLocus.range().is_pos_intersect(siteLocus.pos));

    // If the overlapping indel has any filters, annotate the site with a site conflict
    // (note that we formerly had the site inherit indel filters, but this interacts
    // poorly with empirical scoring)
    //
    // Apply this site filtration at both the locus and sample level
    //
    if (indelLocus.filters.any())
    {
        siteLocus.filters.set(GERMLINE_VARIANT_VCF_FILTERS::SiteConflict);
    }

    const unsigned sampleCount(siteLocus.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& indelSampleInfo(indelLocus.getSample(sampleIndex));

        if (indelSampleInfo.filters.any())
        {
            auto& siteSampleInfo(siteLocus.getSample(sampleIndex));
            siteSampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::SiteConflict);
        }
    }

    // limit qual and gqx values to those of the indel:
    siteLocus.anyVariantAlleleQuality = std::min(siteLocus.anyVariantAlleleQuality, indelLocus.anyVariantAlleleQuality);

    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& indelSampleInfo(indelLocus.getSample(sampleIndex));
        auto& sampleInfo(siteLocus.getSample(sampleIndex));

        sampleInfo.gqx = std::min(sampleInfo.gqx, indelSampleInfo.gqx);
    }

    // after these changes we need to rerun the site filters:
    siteLocus.clearEVSFeatures();
    model.classify_site(siteLocus);
}



void
VariantOverlapResolver::
annotateVariantIndelRecordsAsConflicting()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
#endif

    assert(_variantIndelBuffer.size()>1);

    for (auto& indelLocusPtr : _variantIndelBuffer)
    {
        annotateVariantAsIndelConflict(*indelLocusPtr);
    }
}

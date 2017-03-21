//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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
 *
 *  Created on: Jun 3, 2015
 *      Author: jduddy
 */

#include "indel_overlapper.hh"
#include "blt_util/log.hh"
#include "ScoringModelManager.hh"

//#define DEBUG_GVCF



void
indel_overlapper::
process(std::unique_ptr<GermlineSiteLocusInfo> siteLocusPtr)
{
    std::unique_ptr<GermlineDiploidSiteLocusInfo> si(downcast<GermlineDiploidSiteLocusInfo>(std::move(siteLocusPtr)));
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " SITE START\n";
    log_os << "CHIRP: " << __FUNCTION__ << " pos/indel_range: " << si->pos << "/" << _indel_range << "\n";
#endif

    if (_indel_range.is_pos_intersect(si->pos))
    {
        _site_buffer.push_back(std::move(si));
        return;
    }
    else if (si->pos >= _indel_range.end_pos())
    {
        // resolve any current or previous indels before queuing site:
        process_overlaps();
    }

    assert(_nonvariant_indel_buffer.empty());

    _sink->process(std::move(si));
}



void
indel_overlapper::
process(std::unique_ptr<GermlineIndelLocusInfo> indelLocusPtr)
{
    const bool isVariantLocus(indelLocusPtr->isVariantLocus());
    const bool isForcedOutputLocus(indelLocusPtr->isAnyForcedOutputAtLocus());

    // only handle variant or forced calls
    if (not (isVariantLocus or isForcedOutputLocus)) return;

#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " INDEL START\n";
    log_os << "CHIRP: " << __FUNCTION__ << " pos/indel_range: " << indelLocusPtr->pos << "/" << _indel_range << "\n";
    log_os << "CHIRP: " << __FUNCTION__ << " isVariant/isForcedOutput: " << isVariantLocus << "/" << isForcedOutputLocus << "\n";
#endif

    if (indelLocusPtr->pos > _indel_range.end_pos())
    {
        process_overlaps();
    }

    if (isVariantLocus)
    {
        if (_indel_buffer.empty())
        {
            _indel_range = known_pos_range2(indelLocusPtr->range());
        }
        else
        {
            _indel_range.merge_range(indelLocusPtr->range());
        }
        _indel_buffer.push_back(std::move(indelLocusPtr));
    }
    else
    {
        // below is similar to normal test of pos intersecting with indel_range, except that
        // indel_range's 'end' is extended by one position:
        if ((indelLocusPtr->pos >= _indel_range.begin_pos()) && (indelLocusPtr->pos <= _indel_range.end_pos()))
        {
            _nonvariant_indel_buffer.push_back(std::move(indelLocusPtr));
        }
        else
        {
            _sink->process(std::move(indelLocusPtr));
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
    // dump function may need to deal with data structure in an intermediate state when certain site
    // and indel pointers are already released (principally if called while building an exception
    // report)
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
indel_overlapper::
dump(std::ostream& os) const
{
    os << "indel_overlapper:"
       << " indel_range: " << _indel_range << "\n";
    dumpLocusBuffer("Site", _site_buffer, os);
    dumpLocusBuffer("VariantIndel", _indel_buffer, os);
    dumpLocusBuffer("NonVariantIndel", _nonvariant_indel_buffer, os);
}



void
indel_overlapper::
process_overlaps()
{
    try
    {
        process_overlaps_impl();
    }
    catch (...)
    {
        log_os << "ERROR: exception caught in process_overlaps()\n";
        dump(log_os);

        // need to clear buffer in case indel_overlapper is in an unstable state, otherwise the flush()
        // call into indel_buffer could trigger another exception/segfault:
        clearBuffers();

        throw;
    }
}


namespace VARQUEUE
{
enum index_t
{
    NONE,
    INDEL,
    NONVARIANT_INDEL,
    SITE
};
}


// this doesn't really generalize or tidy up the (implicit) indel/site priority queue, but
// just dumps the ugliness into one place:
static
VARQUEUE::index_t
nextVariantType(
    const std::vector<std::unique_ptr<GermlineIndelLocusInfo>>& indel_buffer,
    const std::vector<std::unique_ptr<GermlineIndelLocusInfo>>& nonvariant_indel_buffer,
    const std::vector<std::unique_ptr<GermlineDiploidSiteLocusInfo>>& site_buffer,
    const unsigned indel_index,
    const unsigned nonvariant_indel_index,
    const unsigned site_index)
{
    const bool is_indel(indel_index<indel_buffer.size());
    const bool is_nonvariant_indel(nonvariant_indel_index<nonvariant_indel_buffer.size());
    const bool is_site(site_index<site_buffer.size());

    if ((!is_indel) && (!is_nonvariant_indel) && (!is_site))
    {
        return VARQUEUE::NONE;
    }

    const bool AlessB(is_indel && ((! is_nonvariant_indel) || (indel_buffer[indel_index]->pos <= nonvariant_indel_buffer[nonvariant_indel_index]->pos)));
    const bool AlessC(is_indel && ((! is_site) || (indel_buffer[indel_index]->pos <= site_buffer[site_index]->pos)));
    const bool BlessC(is_nonvariant_indel && ((! is_site) || (nonvariant_indel_buffer[nonvariant_indel_index]->pos <= site_buffer[site_index]->pos)));

    if (AlessB)
    {
        if (AlessC)
        {
            return VARQUEUE::INDEL;
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
indel_overlapper::
process_overlaps_impl()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
    dump(log_os);
#endif

    if (_indel_buffer.empty() && _nonvariant_indel_buffer.empty()) return;

    bool is_conflict(false);

    // process conflicting loci (these should be rare)
    if (_indel_buffer.size() > 1)
    {
        // sort the indel buffer by pos:
        std::stable_sort(std::begin(_indel_buffer), std::end(_indel_buffer),
                         [](const indel_ptr_t& lhs, const indel_ptr_t& rhs) -> bool
        {
            return (lhs->range().begin_pos() < rhs->range().begin_pos());
        });

        // mark the whole region as conflicting
        modify_conflict_indel_record();
        is_conflict = true;
    }

    // process sites to be consistent with overlapping indels:
    //

    // check that if anything is in the site buffer, we have at least one variant indel:
    // (this guards the _indel_buffer.front() access below)
    assert(_site_buffer.empty() || (not _indel_buffer.empty()));

    for (auto& siteLocusPtr : _site_buffer)
    {
#ifdef DEBUG_GVCF
        log_os << "CHIRP: indel overlapping site: " << siteLocusPtr->pos << "\n";
#endif
        modify_overlapping_site(*(_indel_buffer.front()), *siteLocusPtr, _scoringModels);
    }

    unsigned indel_index(0);
    unsigned nonvariant_indel_index(0);
    unsigned site_index(0);

    // order all buffered indel and site record output according to VCF formatting rules:
    while (true)
    {
        const VARQUEUE::index_t nextvar = nextVariantType(_indel_buffer,_nonvariant_indel_buffer,_site_buffer,indel_index,nonvariant_indel_index,site_index);

        if      (nextvar == VARQUEUE::NONE)
        {
            break;
        }
        else if (nextvar == VARQUEUE::INDEL)
        {
            _sink->process(std::move(_indel_buffer[indel_index]));
            if (is_conflict)
            {
                // emit each conflict record
                indel_index++;
            }
            else
            {
                // just emit the overlapped or single non-conflict record
                indel_index=_indel_buffer.size();
            }
        }
        else if (nextvar == VARQUEUE::NONVARIANT_INDEL)
        {
            _sink->process(std::move(_nonvariant_indel_buffer[nonvariant_indel_index]));
            nonvariant_indel_index++;
        }
        else if (nextvar == VARQUEUE::SITE)
        {
            _sink->process(std::move(_site_buffer[site_index]));
            site_index++;
        }
        else
        {
            assert(false && "unexpected varqueue type");
        }
    }

    clearBuffers();
}



void
indel_overlapper::
modify_overlapping_site(
    const GermlineIndelLocusInfo& indelLocus,
    GermlineDiploidSiteLocusInfo& siteLocus,
    const ScoringModelManager& model)
{
    if (indelLocus.filters.test(GERMLINE_VARIANT_VCF_FILTERS::IndelConflict))
    {
        modify_indel_conflict_site(siteLocus);
    }
    else
    {
        modify_indel_overlap_site(indelLocus, siteLocus, model);
    }
}



void
indel_overlapper::
modify_indel_overlap_site(
    const GermlineIndelLocusInfo& indelLocus,
    GermlineDiploidSiteLocusInfo& siteLocus,
    const ScoringModelManager& model)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site before: " << siteLocus << "\n";
#endif

    // if overlapping indel has any filters, mark as site conflict
    // (note that we formerly had the site inherit indel filters, but
    // this interacts poorly with empirical scoring)

    // apply at both locus level and sample level:
    if (not indelLocus.filters.none())
    {
        siteLocus.filters.set(GERMLINE_VARIANT_VCF_FILTERS::SiteConflict);
    }

    const unsigned sampleCount(siteLocus.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& indelSampleInfo(indelLocus.getSample(sampleIndex));
        auto& siteSampleInfo(siteLocus.getSample(sampleIndex));

        if (not indelSampleInfo.filters.none())
        {
            siteSampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::SiteConflict);
        }
    }

    const pos_t offset(siteLocus.pos-indelLocus.pos);
    assert(offset>=0);

    // limit qual and gq values to those of the indel, and modify site ploidy:
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
indel_overlapper::
modify_indel_conflict_site(GermlineSiteLocusInfo& siteLocus)
{
    siteLocus.filters.set(GERMLINE_VARIANT_VCF_FILTERS::IndelConflict);
}



void
indel_overlapper::
modify_conflict_indel_record()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
#endif

    assert(_indel_buffer.size()>1);

    for (auto& indelLocusPtr : _indel_buffer)
    {
        indelLocusPtr->filters.set(GERMLINE_VARIANT_VCF_FILTERS::IndelConflict);
    }
}

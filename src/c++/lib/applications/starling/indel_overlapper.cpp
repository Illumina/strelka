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
 *
 *  Created on: Jun 3, 2015
 *      Author: jduddy
 */

#include "indel_overlapper.hh"
#include "blt_util/log.hh"
#include "ScoringModelManager.hh"

//#define DEBUG_GVCF



indel_overlapper::indel_overlapper(const ScoringModelManager& model, const reference_contig_segment& ref, std::shared_ptr<variant_pipe_stage_base> destination)
    : variant_pipe_stage_base(destination)
    , _CM(model)
    , _ref(ref)
    , _indel_end_pos(0)
{
    // this component doesn't make any sense without a destination:
    assert(destination);
}



void
indel_overlapper::
process(std::unique_ptr<GermlineSiteLocusInfo> siteLocusPtr)
{
    std::unique_ptr<GermlineDiploidSiteLocusInfo> si(downcast<GermlineDiploidSiteLocusInfo>(std::move(siteLocusPtr)));

    // resolve any current or previous indels before queuing site:
    if (si->pos>=_indel_end_pos)
    {
        process_overlaps();
    }
    else
    {
        _site_buffer.push_back(std::move(si));
        return;
    }

    assert(si->pos>=_indel_end_pos);
    assert(_nonvariant_indel_buffer.empty());

    _sink->process(std::move(si));
}



void
indel_overlapper::
process(std::unique_ptr<GermlineIndelLocusInfo> indelLocusPtr)
{
    auto ii(downcast<GermlineDiploidIndelLocusInfo>(std::move(indelLocusPtr)));

    const bool isNonVariantLocus(not ii->isVariantLocus());

    // don't handle homozygous reference calls unless genotyping is forced
    if (isNonVariantLocus and (not ii->isAnyForcedOutputAtLocus())) return;

    if (ii->pos>_indel_end_pos)
    {
        process_overlaps();
    }

    if (isNonVariantLocus)
    {
        _nonvariant_indel_buffer.push_back(std::move(ii));
    }
    else
    {
        _indel_end_pos=std::max(_indel_end_pos, ii->end());
        _indel_buffer.push_back(std::move(ii));
    }
}



void
indel_overlapper::
dump(std::ostream& os) const
{
    os << "indel_overlapper:"
       << " nSites: " << _site_buffer.size()
       << " nIndels: " << _indel_buffer.size()
       << " indel_end_pos: " << _indel_end_pos << "\n";
    os << "buffered sites:\n";
    for (const auto& site : _site_buffer)
    {
        os << *site << "\n";
    }

    os << "buffered indels:\n";
    for (const auto& indel : _indel_buffer)
    {
        indel->dump(os);
        os << "\n";
    }
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
    const std::vector<std::unique_ptr<GermlineDiploidIndelLocusInfo>>& indel_buffer,
    const std::vector<std::unique_ptr<GermlineDiploidIndelLocusInfo>>& nonvariant_indel_buffer,
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



void indel_overlapper::process_overlaps_impl()
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: " << __FUNCTION__ << " START\n";
#endif

    if (_indel_buffer.empty() && _nonvariant_indel_buffer.empty()) return;

    bool is_conflict(false);

    // process conflicting loci (these should be rare)
    if (_indel_buffer.size() > 1)
    {
        // mark the whole region as conflicting
        modify_conflict_indel_record();
        is_conflict=true;
    }

    // process sites to be consistent with overlapping indels:
    for (auto& si : _site_buffer)
    {
#ifdef DEBUG_GVCF
        log_os << "CHIRP: indel overlapping site: " << si->pos << "\n";
#endif
        modify_overlapping_site(*_indel_buffer[0], *si, _CM);
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

    _indel_buffer.clear();
    _nonvariant_indel_buffer.clear();
    _site_buffer.clear();
}



void
indel_overlapper::
modify_overlapping_site(
    const GermlineDiploidIndelLocusInfo& indelLocus,
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
    const GermlineDiploidIndelLocusInfo& indelLocus,
    GermlineDiploidSiteLocusInfo& siteLocus,
    const ScoringModelManager& model)
{
#ifdef DEBUG_GVCF
    log_os << "CHIRP: indel_overlap_site smod before: " << si.smod << "\n";
#endif

    // if overlapping indel has any filters, mark as site conflict
    // (note that we formerly had the site inherit indel filters, but
    // this interacts poorly with empirical scoring)

    // apply at both locus level and sample level:
    if (! indelLocus.filters.none())
    {
        siteLocus.filters.set(GERMLINE_VARIANT_VCF_FILTERS::SiteConflict);
    }

    const unsigned sampleCount(siteLocus.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& indelSampleInfo(indelLocus.getSample(sampleIndex));
        auto& siteSampleInfo(siteLocus.getSample(sampleIndex));

        if (! indelSampleInfo.filters.none())
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
        auto& siteSampleInfo(siteLocus.getSiteSample(sampleIndex));

        sampleInfo.gqx = std::min(sampleInfo.gqx, indelSampleInfo.gqx);

        const unsigned indelInducedSitePloidy(indelLocus.getSitePloidy(sampleIndex, offset));

        // change ploidy:
        if (indelInducedSitePloidy == 1)
        {
            /// TODO STREL-125 change site GT to per-sample:
            if (DIGT::is_het(siteSampleInfo.max_gt))
            {
                sampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::SiteConflict);
            }
            else
            {
                if (siteSampleInfo.max_gt == siteLocus.refBaseIndex)
                {
                    siteSampleInfo.modified_gt = MODIFIED_SITE_GT::ZERO;
                }
                else
                {
                    siteSampleInfo.modified_gt = MODIFIED_SITE_GT::ONE;
                }
            }
        }
        else if (indelInducedSitePloidy == 0)
        {
            if (siteSampleInfo.max_gt == siteLocus.refBaseIndex)
            {
                siteSampleInfo.modified_gt = MODIFIED_SITE_GT::UNKNOWN;
                siteSampleInfo.isOverlappingHomAltDeletion = true;
                if (sampleInfo.getPloidy().isNoploid())
                {
                    sampleInfo.filters.unset(GERMLINE_VARIANT_VCF_FILTERS::PloidyConflict);
                }
            }
            else
            {
                sampleInfo.filters.set(GERMLINE_VARIANT_VCF_FILTERS::SiteConflict);
            }
        }
        else if (indelInducedSitePloidy != 2)
        {
            assert(false && "Unexpected ploidy value");
        }
    }

    // after all those changes we need to rerun the site filters:
    siteLocus.clearEVSFeatures();
    model.classify_site(siteLocus);
}



void
indel_overlapper::
modify_indel_conflict_site(GermlineDiploidSiteLocusInfo& siteLocus)
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

    for (auto& ii : _indel_buffer)
    {
        ii->filters.set(GERMLINE_VARIANT_VCF_FILTERS::IndelConflict);
    }
}



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

#include "gvcf_writer.hh"

#include "gvcf_header.hh"
#include "indel_overlapper.hh"
#include "LocusReportInfoUtil.hh"
#include "variant_prefilter_stage.hh"

#include "blt_util/io_util.hh"
#include "blt_util/log.hh"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <blt_common/ref_context.hh>



//#define DEBUG_GVCF


#ifdef DEBUG_GVCF
#include "blt_util/log.hh"
#endif



gvcf_writer::
gvcf_writer(
    const starling_options& opt,
    const starling_deriv_options& dopt,
    const starling_streams& streams,
    const reference_contig_segment& ref,
    const RegionTracker& nocompress_regions,
    const ScoringModelManager& scoringModels)
    : _opt(opt)
    , _streams(streams)
    , _ref(ref)
    , _dopt(dopt.gvcf)
    , _empty_site(_dopt, streams.getSampleCount())
    , _headPos(0)
    , _gvcf_comp(opt.gvcf,nocompress_regions)
    , _scoringModels(scoringModels)
{
    if (! opt.gvcf.is_gvcf_output())
        throw std::invalid_argument("gvcf_writer cannot be constructed with nothing to do.");

    const unsigned sampleCount(_streams.getSampleCount());
    const auto& sampleNames(_streams.getSampleNames());

    if (! _opt.gvcf.is_skip_header)
    {
        finish_gvcf_header(_opt, _dopt, _dopt.chrom_depth, sampleNames, _streams.gvcfVariantsStream());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const std::string& sampleName(sampleNames[sampleIndex]);
            finish_gvcf_header(_opt, _dopt, _dopt.chrom_depth, {sampleName}, _streams.gvcfSampleStream(sampleIndex));
        }
    }

    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        _blockPerSample.emplace_back(_opt.gvcf);
    }

    // add appropriate filters to empty_site, as if it had gone through the standard pipeline:
    _scoringModels.classify_site(_empty_site);
}



void
gvcf_writer::
writeSampleNonVariantBlockRecord(
    const unsigned sampleIndex)
{
    auto& block(_blockPerSample[sampleIndex]);
    if (block.count<=0) return;

    std::ostream& os(_streams.gvcfSampleStream(sampleIndex));
    write_site_record(block, os);
    block.reset();
}



void
gvcf_writer::
filter_site_by_last_indel_overlap(
    GermlineDiploidSiteLocusInfo& locus)
{
    if (_last_indel)
    {
        if (locus.pos >= _last_indel->end())
        {
            _last_indel.reset(nullptr);
        }
        else
        {
            indel_overlapper::modify_overlapping_site(*_last_indel, locus, _scoringModels);
        }
    }
}



void
gvcf_writer::
skip_to_pos(
    const pos_t target_pos)
{
    // advance through any indel region by adding individual sites
    while (_headPos<target_pos)
    {
        GermlineDiploidSiteLocusInfo si = get_empty_site(_headPos);

        add_site_internal(si);
        // Don't do compressed ranges if there is an overlapping indel
        // filters are being applied to the overlapping positions
        if (_last_indel) continue;

        if (_gvcf_comp.is_range_compressible(known_pos_range2(si.pos, target_pos)))
        {
            const int deltapos(target_pos - _headPos);
            for (auto& block : _blockPerSample)
            {
                assert(block.count != 0);
                block.count += deltapos;
            }
            _headPos= target_pos;
        }
    }
}



template <typename T>
void
safeLocusDump(
    const std::unique_ptr<T>& locusPtr)
{
    if (locusPtr)
    {
        log_os << *locusPtr;
    }
    else
    {
        log_os << "ALREADY RELEASED";
    }
    log_os << "\n";

}


void
gvcf_writer::
process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    try
    {
        assert(locusPtr->getSampleCount() == getSampleCount());

        skip_to_pos(locusPtr->pos);
        add_site_internal(*locusPtr);
    }
    catch (...)
    {
        log_os << "ERROR: Exception caught in gvcf_writer while processing site:\n";
        safeLocusDump(locusPtr);
        throw;
    }
}



void
gvcf_writer::
process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
{
    try
    {
        skip_to_pos(locusPtr->pos);

        // flush any non-variant block before starting:
        writeAllNonVariantBlockRecords();

        write_indel_record(*locusPtr);
        if (dynamic_cast<GermlineDiploidIndelLocusInfo*>(locusPtr.get()) != nullptr)
        {
            _last_indel = std::move(locusPtr);
        }
    }
    catch (...)
    {
        log_os << "ERROR: Exception caught in gvcf_writer while processing indel:\n";
        safeLocusDump(locusPtr);
        throw;
    }
}



void
gvcf_writer::
flush_impl()
{
    skip_to_pos(_reportRange.end_pos());
    writeAllNonVariantBlockRecords();
}



void
gvcf_writer::
add_site_internal(
    GermlineSiteLocusInfo& locus)
{
    GermlineDiploidSiteLocusInfo* diploidLocusPtr(dynamic_cast<GermlineDiploidSiteLocusInfo*>(&locus));
    if (diploidLocusPtr != nullptr)
    {
        filter_site_by_last_indel_overlap(*diploidLocusPtr);
    }

    _headPos=locus.pos+1;

    // write_site
    queue_site_record(locus);
}



/// queue site record for writing, after
/// possibly joining it into a compressed non-variant block
///
void
gvcf_writer::
queue_site_record(
    const GermlineSiteLocusInfo& locus)
{
    //test for basic blocking criteria
    if (!_gvcf_comp.is_site_compressible(locus))
    {
        writeAllNonVariantBlockRecords();
        write_site_record(locus);
        return;
    }

    const unsigned sampleCount(getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        gvcf_block_site_record& block(_blockPerSample[sampleIndex]);
        if (! block.testCanSiteJoinSampleBlock(locus, sampleIndex))
        {
            writeSampleNonVariantBlockRecord(sampleIndex);
        }
        block.joinSiteToSampleBlock(locus, sampleIndex);
    }
}



/// extend the locus filter set such that any sample filter applied to all samples is added to the locus filters
static
GermlineFilterKeeper
getExtendedLocusFilters(
    const LocusInfo& locus)
{
    GermlineFilterKeeper locusFilters = locus.filters;
    GermlineFilterKeeper sampleFilterUnion;
    const unsigned sampleCount(locus.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
    {
        const auto& sampleInfo(locus.getSample(sampleIndex));
        if (sampleIndex == 0)
        {
            sampleFilterUnion = sampleInfo.filters;
        }
        else
        {
            sampleFilterUnion.unionMerge(sampleInfo.filters);
        }
    }
    locusFilters.merge(sampleFilterUnion);
    return locusFilters;
}



static
void
writeSiteVcfAltField(
    const std::vector<GermlineSiteAlleleInfo>& siteAlleles,
    std::ostream& os)
{
    if (siteAlleles.empty())
    {
        os << '.';
    }
    else
    {
        const unsigned altAlleleCount(siteAlleles.size());
        for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; altAlleleIndex++)
        {
            if (altAlleleIndex != 0) os <<',';
            os << id_to_base(siteAlleles[altAlleleIndex].baseIndex);
        }
    }

}



/// print sample AD/ADF/ADR
static
void
printSampleAD(
    const LocusSupportingReadStats& counts,
    const unsigned expectedAltAlleleCount,
    std::ostream& os)
{
    // verify locus and sample allele counts are in sync:
    assert(counts.getAltCount() == expectedAltAlleleCount);
    const unsigned fullAlleleCount(expectedAltAlleleCount+1);

    // AD
    os << ':';
    for (unsigned alleleIndex(0); alleleIndex < fullAlleleCount; ++alleleIndex)
    {
        if (alleleIndex>0) os << ',';
        os << (counts.getCounts(true).confidentAlleleCount(alleleIndex) + counts.getCounts(false).confidentAlleleCount(alleleIndex));
    }

    // ADF/ADR
    static const unsigned strandCount(2);
    for (unsigned strandIndex(0); strandIndex < strandCount; ++strandIndex)
    {
        const bool isFwdStrand(strandIndex==0);
        const auto& strandCounts(counts.getCounts(isFwdStrand));

        os << ':';
        for (unsigned alleleIndex(0); alleleIndex < fullAlleleCount; ++alleleIndex)
        {
            if (alleleIndex>0) os << ',';
            os << strandCounts.confidentAlleleCount(alleleIndex);
        }
    }
}



void
gvcf_writer::
write_site_record_instance(
    const GermlineSiteLocusInfo& locus,
    std::ostream& os,
    const int targetSampleIndex) const
{
    const auto& siteAlleles(locus.getSiteAlleles());
    const unsigned altAlleleCount(siteAlleles.size());
    const bool isAltAlleles(altAlleleCount > 0);
    const unsigned sampleCount(locus.getSampleCount());

    os << getChromName() << '\t'  // CHROM
       << (locus.pos + 1) << '\t'  // POS
       << ".\t";           // ID

    os << id_to_base(locus.refBaseIndex) << '\t'; // REF

    // ALT
    writeSiteVcfAltField(locus.getSiteAlleles(), os);
    os << '\t';

    // QUAL:
    if (locus.isQual())
    {
        os << locus.anyVariantAlleleQuality;
    }
    else
    {
        os << '.';
    }
    os << '\t';

    // FILTER:
    getExtendedLocusFilters(locus).write(os);
    os << '\t';

    // INFO:
    // SNVHPOL
    {
        unsigned hpol(locus.hpol);
        if (not locus.isVariantLocus())
        {
            // we never computed this upfront for non-variants (b/c not running EVS and saves time)
            hpol = get_snp_hpol_size(locus.pos, _ref);
        }
        os << "SNVHPOL=" << hpol;
    }
    os << ';';

    // MQ
    {
        // compute global MQ over all samples
        MapqTracker mapqTracker;
        {
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));
                mapqTracker.merge(siteSampleInfo.mapqTracker);
            }
        }
        os << "MQ=" << std::lround(mapqTracker.getRMS());
    }

    const GermlineDiploidSiteLocusInfo* diploidLocusPtr(dynamic_cast<const GermlineDiploidSiteLocusInfo*>(&locus));
    if (diploidLocusPtr != nullptr)
    {
        const GermlineDiploidSiteLocusInfo& diploidLocus(*diploidLocusPtr);

        if (locus.isVariantLocus())
        {
            if (_opt.isReportEVSFeatures)
            {
                // EVS features may not be computed for certain records, so check first:
                if (not diploidLocus.evsFeatures.empty())
                {
                    const StreamScoper ss(os);
                    os << std::setprecision(5);
                    os << ";EVSF=";
                    diploidLocus.evsFeatures.writeValues(os);
                    os << ",";
                    diploidLocus.evsDevelopmentFeatures.writeValues(os);
                }
            }
        }
        os << '\t';

        //FORMAT
        os << "GT:GQ:GQX:DP:DPF";
        if (isAltAlleles)
        {
            os << ":AD:ADF:ADR";
        }
        if (isAltAlleles)
        {
            os << ":SB";
        }
        os << ":FT";
        if (isAltAlleles)
        {
            os << ":PL";
        }

        bool isAnyPhased(false);
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            const auto& sampleInfo(locus.getSample(sampleIndex));
            if (sampleInfo.phaseSetId >= 0)
            {
                isAnyPhased = true;
                break;
            }
        }

        if (isAnyPhased)
        {
            os << ":PS";
        }

        //SAMPLE
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            if (targetSampleIndex >= 0)
            {
                if (static_cast<int>(sampleIndex) != targetSampleIndex) continue;
            }

            const auto& sampleInfo(locus.getSample(sampleIndex));
            const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

            os << '\t';

            os << sampleInfo.max_gt() << ':';
            if (locus.is_gqx(sampleIndex))
            {
                os << sampleInfo.genotypeQualityPolymorphic << ':';
                os << ((sampleInfo.empiricalVariantScore >= 0) ? sampleInfo.empiricalVariantScore : sampleInfo.gqx);
            }
            else
            {
                os << ".:.";
            }
            os << ':';
            //print DP:DPF
            os << siteSampleInfo.n_used_calls << ':'
               << siteSampleInfo.n_unused_calls;

            // AD/ADF/ADR
            if (isAltAlleles)
            {
                printSampleAD(sampleInfo.supportCounts, altAlleleCount, os);
            }

            // SB
            if (isAltAlleles)
            {
                os << ":";
                const StreamScoper ss(os);
                os << std::fixed << std::setprecision(1) << siteSampleInfo.strandBias;
            }

            // FT
            os << ':';
            sampleInfo.filters.write(os);

            // PL
            if (isAltAlleles)
            {
                const bool isUnknonwnGT(sampleInfo.max_gt().isUnknown());

                os << ':';
                if (isUnknonwnGT)
                {
                    os << '.';
                }
                else
                {
                    bool isFirst(true);
                    for (const auto pls : sampleInfo.genotypePhredLoghood)
                    {
                        if (isFirst)
                        {
                            isFirst = false;
                        }
                        else
                        {
                            os << ',';
                        }
                        os << std::min(pls, maxPL);
                    }
                }
            }

            // PS
            if (isAnyPhased)
            {
                os << ':';
                if (sampleInfo.phaseSetId < 0)
                {
                    os << '.';
                }
                else
                {
                    os << sampleInfo.phaseSetId;
                }
            }
        }
    }
    else
    {
        // special constraint on continuous allele reporting right now:
        assert(altAlleleCount == 1);

        const GermlineContinuousSiteLocusInfo* contLocusPtr(dynamic_cast<const GermlineContinuousSiteLocusInfo*>(&locus));
        assert(contLocusPtr != nullptr);
        const GermlineContinuousSiteLocusInfo& contLocus(*contLocusPtr);

        os << '\t';

        //FORMAT
        os << "GT";
        os << ":GQ";
        os << ":GQX";
        os << ":DP:DPF";
        if (isAltAlleles)
        {
            os << ":AD:ADF:ADR";
        }
        if (isAltAlleles)
        {
            os << ":SB";
        }
        os << ":FT:VF";

        //SAMPLE
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            if (targetSampleIndex >= 0)
            {
                if (static_cast<int>(sampleIndex) != targetSampleIndex) continue;
            }

            const auto& sampleInfo(locus.getSample(sampleIndex));
            const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

            os << '\t';

            VcfGenotypeUtil::writeGenotype(sampleInfo.max_gt(),os);

            //SAMPLE
            os << ':' << sampleInfo.genotypeQualityPolymorphic
               << ':' << sampleInfo.gqx;

            // DP:DPF
            os << ':' << siteSampleInfo.n_used_calls << ':' << siteSampleInfo.n_unused_calls;

            // AD/ADF/ADR
            if (isAltAlleles)
            {
                printSampleAD(sampleInfo.supportCounts, altAlleleCount, os);
            }

            // SB
            if (isAltAlleles)
            {
                const StreamScoper ss(os);
                os << ':';
                os << std::fixed << std::setprecision(1) << siteSampleInfo.strandBias;
            }

            // FT
            os << ':';
            sampleInfo.filters.write(os);

            // VF
            {
                const auto& continuousSiteSampleInfo(contLocus.getContinuousSiteSample(sampleIndex));
                const StreamScoper ss(os);
                os << ':' << std::fixed << std::setprecision(3) << continuousSiteSampleInfo.getContinuousAlleleFrequency();
            }

        }
    }

    os << '\n';
}



void
gvcf_writer::
write_site_record(
    const GermlineSiteLocusInfo& locus) const
{
    const unsigned sampleCount(locus.getSampleCount());

    write_site_record_instance(locus, _streams.gvcfVariantsStream());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        write_site_record_instance(locus, _streams.gvcfSampleStream(sampleIndex), sampleIndex);
    }
}



void
gvcf_writer::
write_site_record(
    const gvcf_block_site_record& locus,
    std::ostream& os) const
{
    os << getChromName() << '\t'  // CHROM
       << (locus.pos+1) << '\t'  // POS
       << ".\t";           // ID

    os  << id_to_base(locus.refBaseIndex) << '\t'; // REF

    // ALT
    os << '.';
    os << '\t';

    // QUAL:
    os << '.';
    os << '\t';

    // FILTER:
    getExtendedLocusFilters(locus).write(os);
    os << '\t';

    // INFO:
    if (locus.count>1)
    {
        os << "END=" << (locus.pos+locus.count) << ';';
        os << _dopt.block_label;
    }
    else
    {
        os << '.';
    }
    os << '\t';

    //FORMAT
    os << "GT";
    os << ":GQX:DP:DPF";
    os << '\t';

    //SAMPLE
    // there should be exactly one sample:
    assert(1 == locus.getSampleCount());
    const auto& sampleInfo(locus.getSample(0));
    os << sampleInfo.max_gt() << ':';
    if (locus.isBlockGqxDefined)
    {
        os << locus.block_gqx.min();
    }
    else
    {
        os << '.';
    }
    os << ':';
    //print DP:DPF
    os << locus.block_dpu.min() << ':'
       << locus.block_dpf.min();
    os << '\n';
}



void
gvcf_writer::
write_indel_record_instance(
    const GermlineIndelLocusInfo& locus,
    std::ostream& os,
    const int targetSampleIndex) const
{
    const unsigned sampleCount(locus.getSampleCount());

    // create VCF specific transformation of the alt allele list
    const auto& indelAlleles(locus.getIndelAlleles());
    OrthogonalAlleleSetLocusReportInfo locusReportInfo;
    getLocusReportInfoFromAlleles(_ref, indelAlleles, locusReportInfo);

    os << getChromName() << '\t'   // CHROM
       << locusReportInfo.vcfPos << '\t'   // POS
       << ".\t"            // ID
       << locusReportInfo.vcfRefSeq << '\t'; // REF

    // ALT
    const unsigned altAlleleCount(locus.getAltAlleleCount());

    for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
    {
        if (altAlleleIndex > 0) os << ',';
        os << locusReportInfo.altAlleles[altAlleleIndex].vcfAltSeq;
    }
    os << '\t';

    os << locus.anyVariantAlleleQuality << '\t'; //QUAL

    // FILTER:
    getExtendedLocusFilters(locus).write(os);
    os << '\t';

    // INFO
    os << "CIGAR=";
    for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
    {
        if (altAlleleIndex > 0) os << ',';
        os << locusReportInfo.altAlleles[altAlleleIndex].vcfCigar;
    }
    os << ';';
    os << "RU=";
    for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
    {
        if (altAlleleIndex > 0) os << ',';
        const auto& iri(indelAlleles[altAlleleIndex].indelReportInfo);
        if (iri.is_repeat_unit() && iri.repeat_unit.size() <= 20)
        {
            os << iri.repeat_unit;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "REFREP=";
    for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
    {
        if (altAlleleIndex > 0) os << ',';
        const auto& iri(indelAlleles[altAlleleIndex].indelReportInfo);
        if (iri.is_repeat_unit())
        {
            os << iri.ref_repeat_count;
        }
        else
        {
            os << '.';
        }
    }
    os << ';';
    os << "IDREP=";
    for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
    {
        if (altAlleleIndex > 0) os << ',';
        const auto& iri(indelAlleles[altAlleleIndex].indelReportInfo);
        if (iri.is_repeat_unit())
        {
            os << iri.indel_repeat_count;
        }
        else
        {
            os << '.';
        }
    }

    // compute global MQ over all samples
    MapqTracker mapqTracker;
    {
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            const auto& indelSampleInfo(locus.getIndelSample(sampleIndex));
            mapqTracker.merge(indelSampleInfo.mapqTracker);
        }
    }
    os << ';';
    os << "MQ=" << std::lround(mapqTracker.getRMS());

    const GermlineDiploidIndelLocusInfo* diploidLocusPtr(dynamic_cast<const GermlineDiploidIndelLocusInfo*>(&locus));
    if (diploidLocusPtr != nullptr)
    {
        const GermlineDiploidIndelLocusInfo& diploidLocus(*diploidLocusPtr);

        //FORMAT
        if (_opt.isReportEVSFeatures)
        {
            // EVS features may not be computed for certain records, so check first:
            if (! diploidLocus.evsFeatures.empty())
            {
                const StreamScoper ss(os);
                os << std::setprecision(5);
                os << ";EVSF=";
                diploidLocus.evsFeatures.writeValues(os);
                os << ",";
                diploidLocus.evsDevelopmentFeatures.writeValues(os);
            }
        }

        os << '\t';

        //FORMAT
        os << "GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL";

        //SAMPLE
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            if (targetSampleIndex >= 0)
            {
                if (static_cast<int>(sampleIndex) != targetSampleIndex) continue;
            }

            const auto& sampleInfo(locus.getSample(sampleIndex));
            const auto& indelSampleInfo(locus.getIndelSample(sampleIndex));

            os << '\t';

            os << sampleInfo.max_gt();
            os << ':' << sampleInfo.genotypeQualityPolymorphic;

            os << ':' << ((sampleInfo.empiricalVariantScore >= 0) ? sampleInfo.empiricalVariantScore : sampleInfo.gqx);

            os << ':' << indelSampleInfo.tier1Depth;

            printSampleAD(sampleInfo.supportCounts, altAlleleCount, os);

            // FT
            os << ':';
            sampleInfo.filters.write(os);

            // PL
            os << ":";
            bool isFirst(true);
            for (const auto pls : sampleInfo.genotypePhredLoghood)
            {
                if (isFirst)
                {
                    isFirst = false;
                }
                else
                {
                    os << ',';
                }
                os << std::min(pls, maxPL);
            }
        }
    }
    else
    {
        // special constraint on continuous allele reporting right now:
        assert(altAlleleCount == 1);

        os << '\t';

        //FORMAT
        os << "GT:GQ:GQX:DPI:AD:ADF:ADR:FT:VF";

        //SAMPLE
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            if (targetSampleIndex >= 0)
            {
                if (static_cast<int>(sampleIndex) != targetSampleIndex) continue;
            }

            const auto& sampleInfo(locus.getSample(sampleIndex));
            const auto& indelSampleInfo(locus.getIndelSample(sampleIndex));

            os << '\t';

            //SAMPLE
            os << sampleInfo.max_gt();
            os << ':' << sampleInfo.genotypeQualityPolymorphic;

            os << ':' << sampleInfo.gqx;

            os << ':' << indelSampleInfo.tier1Depth;

            {
                const auto& sampleReportInfo(indelSampleInfo.legacyReportInfo);

                // AD:
                os << ':' << sampleReportInfo.n_confident_ref_reads
                   << ',' << sampleReportInfo.n_confident_indel_reads;

                // ADF
                os << ':' << sampleReportInfo.n_confident_ref_reads_fwd
                   << ',' << sampleReportInfo.n_confident_indel_reads_fwd;

                // ADR
                os << ':' << sampleReportInfo.n_confident_ref_reads_rev
                   << ',' << sampleReportInfo.n_confident_indel_reads_rev;
            }

            // FT
            os << ':';
            sampleInfo.filters.write(os);

            // VF
            {
                const StreamScoper ss(os);
                os << ':' << std::setprecision(3) << indelSampleInfo.alleleFrequency();
            }
        }
    }

    os << '\n';
}



void
gvcf_writer::
write_indel_record(
    const GermlineIndelLocusInfo& locus) const
{
    const unsigned sampleCount(locus.getSampleCount());

    write_indel_record_instance(locus, _streams.gvcfVariantsStream());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        write_indel_record_instance(locus, _streams.gvcfSampleStream(sampleIndex), sampleIndex);
    }
}

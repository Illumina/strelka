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

#include "gvcf_writer.hh"

#include "gvcf_header.hh"
#include "VariantOverlapResolver.hh"
#include "LocusReportInfoUtil.hh"
#include "variant_prefilter_stage.hh"

#include "blt_common/ref_context.hh"
#include "blt_util/io_util.hh"
#include "blt_util/log.hh"

#include <iomanip>
#include <iostream>
#include <sstream>



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
    const RegionTracker& nocompressRegions,
    const RegionTracker& callRegions,
    const ScoringModelManager& scoringModels)
    : _opt(opt)
    , _streams(streams)
    , _ref(ref)
    , _dopt(dopt.gvcf)
    , _empty_site(_dopt, streams.getSampleCount())
    , _callRegions(callRegions)
    , _headPos(0)
    , _gvcf_comp(opt.gvcf,nocompressRegions)
    , _scoringModels(scoringModels)
{
    if (! opt.gvcf.is_gvcf_output())
        throw std::invalid_argument("gvcf_writer cannot be constructed with nothing to do.");

    const unsigned sampleCount(_streams.getSampleCount());
    const auto& sampleNames(_streams.getSampleNames());

    if (! _opt.gvcf.is_skip_header)
    {
        // Create the header for the variants VCF
        bool isGenomeVCF(false);
        finishGermlineVCFheader(_opt, _dopt, _dopt.chrom_depth, sampleNames, isGenomeVCF, _streams.variantsVCFStream());

        // Create the header for each sample's gVCF
        isGenomeVCF = true;
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const std::string& sampleName(sampleNames[sampleIndex]);
            finishGermlineVCFheader(_opt, _dopt, _dopt.chrom_depth, {sampleName}, isGenomeVCF,
                                    _streams.gvcfSampleStream(sampleIndex));
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
modifySiteForConsistencyWithUpstreamIndels(
    GermlineDiploidSiteLocusInfo& locus)
{
    if (_lastVariantIndelWritten)
    {
        if (locus.pos >= _lastVariantIndelWritten->end())
        {
            _lastVariantIndelWritten.reset(nullptr);
        }
        else
        {
            VariantOverlapResolver::modifySiteOverlappingVariantIndel(*_lastVariantIndelWritten, locus, _scoringModels);
        }
    }
}



void
gvcf_writer::
skip_to_pos(
    const pos_t target_pos)
{
    // advance through any indel or uncovered region by adding individual empty sites
    while (_headPos<target_pos)
    {
        if (_opt.isUseCallRegions())
        {
            if (not _callRegions.isIntersectRegion(_headPos))
            {
                _headPos++;
                continue;
            }
        }

        GermlineDiploidSiteLocusInfo si = get_empty_site(_headPos);

        add_site_internal(si);

        // Don't do compressed ranges if there is an overlapping indel
        // because filters are being applied to the overlapping positions
        if (_lastVariantIndelWritten) continue;

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
        log_os << "Exception caught in gvcf_writer while processing site:\n";
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

        // flush any non-variant blocks across all samples before handling the indel:
        writeAllNonVariantBlockRecords();

        write_indel_record(*locusPtr);
        if (locusPtr->isVariantLocus())
        {
            if (dynamic_cast<GermlineDiploidIndelLocusInfo*>(locusPtr.get()) != nullptr)
            {
                _lastVariantIndelWritten = std::move(locusPtr);
            }
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
    if (not _chromName.empty())
    {
        skip_to_pos(_reportRange.end_pos());
        writeAllNonVariantBlockRecords();
    }

    _chromName.clear();
    _headPos = 0;
    _lastVariantIndelWritten.reset(nullptr);
}



void
gvcf_writer::
add_site_internal(
    GermlineSiteLocusInfo& locus)
{
    GermlineDiploidSiteLocusInfo* diploidLocusPtr(dynamic_cast<GermlineDiploidSiteLocusInfo*>(&locus));
    if (diploidLocusPtr != nullptr)
    {
        modifySiteForConsistencyWithUpstreamIndels(*diploidLocusPtr);
    }

    _headPos=locus.pos+1;

    // write_site
    queue_site_record(locus);
}



void
gvcf_writer::
queue_site_record(
    const GermlineSiteLocusInfo& locus)
{
    //test for basic blocking criteria
    if (! _gvcf_comp.is_site_compressible(locus))
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



/// \brief Apply final VCF record FILTER customizations based the type of output VCF
///
/// Currently the customizations are:
/// - Extend the locus filter set to add those sample filters which are present in all samples of the VCF
///   (Note that for single-sample gVCFs, this means that all sample filters will be copied to the locus filter field.)
/// - Add a filter to the multi-sample variant VCF output whenever there are no samples with a variant
///    genotype that pass filters, this makes it easy to naively interpret the variant variants.vcf file by just
///    looking at whether FILTER is 'PASS'
///
/// \param targetSampleIndex The sample index. This indicates the index of the sample-specific gVCF to target, or
///                          if the value is less than 0, this signifies processing for the variants VCF.
/// \return Filters to be used for the input locus
///
/// TODO: default sampleIndex of 0 is confusing and needs cleanup or docs, why is there a default at all?
static
GermlineFilterKeeper
getExtendedLocusFilters(
    const LocusInfo& locus,
    const int targetSampleIndex = 0)
{
    const unsigned sampleCount(locus.getSampleCount());

    GermlineFilterKeeper locusFilters = locus.filters;

    // 1. Extend locus filters with sample filters applied to every sample
    {
        bool isFirstSample(true);
        GermlineFilterKeeper sampleFilterIntersection;
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            if (targetSampleIndex >= 0)
            {
                if (static_cast<int>(sampleIndex) != targetSampleIndex) continue;
            }

            const auto& sampleInfo(locus.getSample(sampleIndex));
            if (isFirstSample)
            {
                sampleFilterIntersection = sampleInfo.filters;
                isFirstSample = false;
            }
            else
            {
                sampleFilterIntersection.intersectWith(sampleInfo.filters);
            }
        }

        locusFilters.merge(sampleFilterIntersection);
    }

    // 2. Add a filter in the multi-sample variant VCF to flag records without any passing variant genotypes
    if (targetSampleIndex < 0)
    {
        bool isFilter(true);
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            const auto& sampleInfo(locus.getSample(sampleIndex));
            if (sampleInfo.isVariant())
            {
                if (sampleInfo.filters.none())
                {
                    isFilter = false;
                    break;
                }
            }
        }

        if (isFilter)
        {
            locusFilters.set(GERMLINE_VARIANT_VCF_FILTERS::NoPassedVariantGTs);
        }
    }

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
    getExtendedLocusFilters(locus, targetSampleIndex).write(os);
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
            os << siteSampleInfo.usedBasecallCount << ':'
               << siteSampleInfo.unusedBasecallCount;

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
            os << ':' << siteSampleInfo.usedBasecallCount << ':' << siteSampleInfo.unusedBasecallCount;

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

    write_site_record_instance(locus, _streams.variantsVCFStream());
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
    os << ":GQX:DP:DPF:MIN_DP";
    os << '\t';

    //SAMPLE
    // there should be exactly one sample:
    assert(1 == locus.getSampleCount());
    const auto& sampleInfo(locus.getSample(0));
    os << sampleInfo.max_gt() << ':';
    if (locus.isBlockGqxDefined)
    {
        // The value we want here is the genotype confidence of the entire block. Instead we only have the GT
        // probabilities from each site in the block, so an approximation is used: The minimum site genotype
        // confidence from all sites in the block provides a reasonable upper-bound on our confidence that
        // the block as a whole is 0/0 (in all cases, completely neglecting possible indel evidence in the block)
        os << locus.block_gqx.min();
    }
    else
    {
        os << '.';
    }
    os << ':';
    //print DP:DPF:MIN_DP
    os << std::round(locus.block_dpu.mean()) << ':'
       << std::round(locus.block_dpf.mean()) << ':'
       << locus.block_dpu.min();
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
    getLocusReportInfoFromAlleles(_ref, indelAlleles, locus.getCommonPrefixLength(), locusReportInfo);

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
    getExtendedLocusFilters(locus, targetSampleIndex).write(os);
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
        if (iri.isRepeatUnit() && iri.repeatUnit.size() <= 20)
        {
            os << iri.repeatUnit;
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
        if (iri.isRepeatUnit())
        {
            os << iri.refRepeatCount;
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
        if (iri.isRepeatUnit())
        {
            os << iri.indelRepeatCount;
        }
        else
        {
            os << '.';
        }
    }

    // compute global MQ over all samples
    os << ';';
    if (! locus.isNotGenotyped())
    {
        MapqTracker mapqTracker;
        {
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const auto& indelSampleInfo(locus.getIndelSample(sampleIndex));
                mapqTracker.merge(indelSampleInfo.mapqTracker);
            }
        }
        os << "MQ=" << std::lround(mapqTracker.getRMS());
    }
    else
    {
        os << "MQ=.";
    }


    const GermlineDiploidIndelLocusInfo* diploidLocusPtr(dynamic_cast<const GermlineDiploidIndelLocusInfo*>(&locus));
    if (diploidLocusPtr != nullptr)
    {
        if (!diploidLocusPtr->isNotGenotyped())
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
            // not genotyped
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

                os << '\t';
                os << ".:.:.:.:.:.:.";

                // FT
                os << ':';
                sampleInfo.filters.write(os);

                // PL
                os << ":.";
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

            printSampleAD(sampleInfo.supportCounts, altAlleleCount, os);

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
    if (! _reportRange.is_pos_intersect(locus.pos)) return;
    if (_opt.isUseCallRegions())
    {
        if (! _callRegions.isIntersectRegion(locus.pos)) return;
    }

    const unsigned sampleCount(locus.getSampleCount());

    write_indel_record_instance(locus, _streams.variantsVCFStream());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        write_indel_record_instance(locus, _streams.gvcfSampleStream(sampleIndex), sampleIndex);
    }
}

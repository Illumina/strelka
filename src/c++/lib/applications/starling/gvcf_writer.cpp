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
#include "ScoringModelManager.hh"
#include "variant_prefilter_stage.hh"

#include "blt_util/io_util.hh"

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
    const reference_contig_segment& ref,
    const RegionTracker& nocompress_regions,
    const std::vector<std::string>& sampleNames,
    std::ostream* osptr,
    const ScoringModelManager& cm)
    : _opt(opt)
    , _report_range(dopt.report_range.begin_pos,dopt.report_range.end_pos)
    , _ref(ref)
    , _osptr(osptr)
    , _chrom(opt.bam_seq_name.c_str())
    , _dopt(dopt.gvcf)
    , _head_pos(dopt.report_range.begin_pos)
    , _empty_site(_dopt, sampleNames.size())
    , _gvcf_comp(opt.gvcf,nocompress_regions)
    , _CM(cm)
{
    assert(_report_range.is_begin_pos);
    assert(_report_range.is_end_pos);

    if (! opt.gvcf.is_gvcf_output())
        throw std::invalid_argument("gvcf_writer cannot be constructed with nothing to do.");

    assert(nullptr != _osptr);
    assert((nullptr !=_chrom) && (strlen(_chrom)>0));
    assert(not sampleNames.empty());

    if (! _opt.gvcf.is_skip_header)
    {
        /// TODO STREL-125 initialize output file array
        const std::string& sampleName(sampleNames.front());
        finish_gvcf_header(_opt, _dopt, _dopt.chrom_depth, sampleName, *_osptr);
    }

    const unsigned sampleCount(sampleNames.size());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        _blockPerSample.emplace_back(_opt.gvcf);
    }
}



void
gvcf_writer::
writeSampleNonVariantBlockRecord(
    const unsigned sampleIndex)
{
    auto& block(_blockPerSample[sampleIndex]);
    if (block.count<=0) return;
    write_site_record(block);
    block.reset();
}



void gvcf_writer::filter_site_by_last_indel_overlap(GermlineDiploidSiteLocusInfo& locus)
{
    if (_last_indel)
    {
        if (locus.pos >= _last_indel->end())
        {
            _last_indel.reset(nullptr);
        }
        else
        {
            indel_overlapper::modify_overlapping_site(*_last_indel, locus, _CM);
        }
    }
}



void
gvcf_writer::
skip_to_pos(const pos_t target_pos)
{
    // advance through any indel region by adding individual sites
    while (_head_pos<target_pos)
    {
        GermlineDiploidSiteLocusInfo si = get_empty_site(_head_pos);

        add_site_internal(si);
        // Don't do compressed ranges if there is an overlapping indel
        // filters are being applied to the overlapping positions
        if (_last_indel) continue;

        if (_gvcf_comp.is_range_compressable(known_pos_range2(si.pos,target_pos)))
        {
            const int deltapos(target_pos - _head_pos);
            for (auto& block : _blockPerSample)
            {
                assert(block.count != 0);
                block.count += deltapos;
            }
            _head_pos= target_pos;
        }
    }
}


void gvcf_writer::process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr)
{
    assert(locusPtr->getSampleCount() == getSampleCount());
    
    skip_to_pos(locusPtr->pos);

    if (dynamic_cast<GermlineDiploidSiteLocusInfo*>(locusPtr.get()) != nullptr)
    {
        add_site_internal(*downcast<GermlineDiploidSiteLocusInfo>(std::move(locusPtr)));
    }
    else
    {
        add_site_internal(*downcast<GermlineContinuousSiteLocusInfo>(std::move(locusPtr)));
    }

}



void
gvcf_writer::
process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr)
{
    skip_to_pos(locusPtr->pos);

    // flush any non-variant block before starting:
    writeAllNonVariantBlockRecords();

    write_indel_record(*locusPtr);
    if (dynamic_cast<GermlineDiploidIndelLocusInfo*>(locusPtr.get()) != nullptr)
    {
        _last_indel = downcast<GermlineDiploidIndelLocusInfo>(std::move(locusPtr));
    }
}



void
gvcf_writer::
flush_impl()
{
    skip_to_pos(_report_range.end_pos);
    writeAllNonVariantBlockRecords();
}


//Add sites to queue for writing to gVCF
void
gvcf_writer::
add_site_internal(
    GermlineDiploidSiteLocusInfo& locus)
{
    filter_site_by_last_indel_overlap(locus);
    if (locus.isPhasedRegion)
    {
        _head_pos=locus.pos+locus.phased_ref.length();
    }
    else
    {
        _head_pos=locus.pos+1;
    }
    // write_site
    queue_site_record(locus);
}

void
gvcf_writer::
add_site_internal(
    GermlineContinuousSiteLocusInfo& locus)
{
    // TODO: phasing
    _head_pos=locus.pos+1;
    // write_site
    queue_site_record(locus);
}




static
void
get_visible_alt_order(
    const GermlineDiploidSiteLocusInfo& locus,
    std::vector<uint8_t>& altOrder)
{
    altOrder.clear();

    /// TODO STREL-125 generalize to samples
    static const unsigned sampleIndex(0);
    const auto& siteSample(locus.getSiteSample(sampleIndex));

    // list max_gt alts first:
    const uint8_t refBaseId(base_to_id(locus.ref));
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==refBaseId) continue;
        if (! DIGT::expect2(b, siteSample.max_gt)) continue;
        altOrder.push_back(b);
    }

#if 0
    // include other alts based on known count:
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==si.dgt.ref_gt) continue;
        if (DIGT::expect2(b,si.smod.max_gt)) continue;
        if (si.known_counts[b] > 0) altOrder.push_back(b);
    }
#endif
}



static
void
print_vcf_alt(
    const std::vector<uint8_t>& altOrder,
    std::ostream& os)
{
    bool is_print(false);
    for (const auto& b : altOrder)
    {
        if (is_print) os << ',';
        os << id_to_base(b);
        is_print=true;
    }
    if (! is_print) os << '.';
}



static
void
print_site_ad(
    const uint8_t refBaseId,
    const GermlineSiteSampleInfo& siteSampleInfo,
    const std::vector<uint8_t>& altOrder,
    std::ostream& os)
{
    os << siteSampleInfo.alleleObservationCounts(refBaseId);

    for (const auto& b : altOrder)
    {
        os << ',' << siteSampleInfo.alleleObservationCounts(b);
    }
}



/// extend the locus filter set such that any sample filter applied to all samples is added to the locus filters
static
GermlineFilterKeeper
getExtendedLocusFilters(const LocusInfo& locus)
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
print_site_ad_strand(
    const uint8_t refBaseId,
    const GermlineSiteSampleInfo& siteSampleInfo,
    const std::vector<uint8_t>& altOrder,
    const bool is_fwd_strand,
    std::ostream& os)
{
    os << siteSampleInfo.alleleObservationCountsByStrand(is_fwd_strand, refBaseId);

    for (const auto& b : altOrder)
    {
        os << ',' << siteSampleInfo.alleleObservationCountsByStrand(is_fwd_strand,b);
    }
}



//writes out a SNP or block record
void
gvcf_writer::
write_site_record(
    const GermlineDiploidSiteLocusInfo& locus) const
{
    std::ostream& os(*_osptr);

    const auto& siteAlleles(locus.getSiteAlleles());
    const unsigned altAlleleCount(siteAlleles.size());
    const bool isAltAlleles(altAlleleCount > 0);
    const unsigned sampleCount(locus.getSampleCount());

    const uint8_t refBaseId(base_to_id(locus.ref));

    os << _chrom << '\t'  // CHROM
       << (locus.pos+1) << '\t'  // POS
       << ".\t";           // ID

    if (locus.isPhasedRegion)
    {
        os  << locus.phased_ref << '\t'; // REF
    }
    else
    {
        os  << locus.ref << '\t'; // REF
    }

    // ALT
    std::vector<uint8_t> altOrder;
    if (not isAltAlleles)
    {
        os << '.';
    }
    else if (locus.isPhasedRegion)
    {
        os << locus.phased_alt;
    }
    else
    {
        get_visible_alt_order(locus,altOrder);
        print_vcf_alt(altOrder,os);
    }
    os << '\t';

    // QUAL:
    if (locus.is_qual())
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
    if (locus.isVariantLocus())
    {
        os << "SNVSB=";
        {
            /// TODO STREL-125 generalize to multiple alts
            const StreamScoper ss(os);
            os << std::fixed << std::setprecision(1) << locus.allele.strandBias;
        }
        os << ';';
        os << "SNVHPOL=" << locus.hpol;

        // compute global MQ over all samples
        MapqTracker mapqTracker;
        {
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));
                mapqTracker.merge(siteSampleInfo.mapqTracker);
            }
        }
        os << ';';
        os << "MQ=" << mapqTracker.getRMS();

        if (_opt.isReportEVSFeatures)
        {
            // EVS features may not be computed for certain records, so check first:
            if (! locus.evsFeatures.empty())
            {
                const StreamScoper ss(os);
                os << std::setprecision(5);
                os << ";EVSF=";
                locus.evsFeatures.writeValues(os);
                os << ",";
                locus.evsDevelopmentFeatures.writeValues(os);
            }
        }
    }
    else
    {
        os << '.';
    }
    os << '\t';

    //FORMAT
    os << "GT";
    if (locus.isVariantLocus())
    {
        os << ":GQ";
    }
    os << ":GQX:DP:DPF";
    if (isAltAlleles)
    {
        os << ":AD:ADF:ADR";
    }
    os << ":FT";
    if (isAltAlleles)
    {
        os << ":PL";
    }

    //SAMPLE
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& sampleInfo(locus.getSample(sampleIndex));
        const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

        os << '\t';

        os << locus.get_gt(sampleIndex) << ':';
        if (locus.isVariantLocus())
        {
            os << sampleInfo.genotypeQualityPolymorphic << ':';
        }

        if (locus.is_gqx(sampleIndex))
        {
            os << ((sampleInfo.empiricalVariantScore >= 0) ? sampleInfo.empiricalVariantScore : sampleInfo.gqx);
        }
        else
        {
            os << '.';
        }
        os << ':';
        //print DP:DPF
        os << siteSampleInfo.n_used_calls << ':'
           << siteSampleInfo.n_unused_calls;

        if (isAltAlleles)
        {
            if (locus.isPhasedRegion)
            {
                os << ':' << locus.phased_AD
                   << ':' << locus.phased_ADF
                   << ':' << locus.phased_ADR;
            }
            else
            {
                os << ':';
                print_site_ad(refBaseId, siteSampleInfo, altOrder, os);
                os << ':';
                print_site_ad_strand(refBaseId, siteSampleInfo, altOrder, true, os);
                os << ':';
                print_site_ad_strand(refBaseId, siteSampleInfo, altOrder, false, os);
            }
        }

        // FT
        os << ':';
        sampleInfo.filters.write(os);

        if (isAltAlleles)
        {
            // print PL values
            os << ':';
            if (locus.is_hetalt(sampleIndex))
            {
                const unsigned print_gt(siteSampleInfo.max_gt);
                const uint8_t a0(DIGT::get_allele(print_gt, 0));
                const uint8_t a1(DIGT::get_allele(print_gt, 1));
                os << siteSampleInfo.dgt.phredLoghood[refBaseId] << ','
                   << siteSampleInfo.dgt.phredLoghood[DIGT::get_gt_with_alleles(refBaseId, a0)] << ','
                   << siteSampleInfo.dgt.phredLoghood[DIGT::get_gt_with_alleles(a0, a0)] << ','
                   << siteSampleInfo.dgt.phredLoghood[DIGT::get_gt_with_alleles(refBaseId, a1)] << ','
                   << siteSampleInfo.dgt.phredLoghood[DIGT::get_gt_with_alleles(a0, a1)] << ','
                   << siteSampleInfo.dgt.phredLoghood[DIGT::get_gt_with_alleles(a1, a1)];
            }
            else if (sampleInfo.getPloidy().isHaploid() || (siteSampleInfo.modified_gt == MODIFIED_SITE_GT::ONE))
            {
                os << siteSampleInfo.dgt.phredLoghood[refBaseId] << ','
                   << siteSampleInfo.dgt.phredLoghood[siteSampleInfo.max_gt];
            }
            else
            {
                const unsigned print_gt(siteSampleInfo.max_gt);
                const uint8_t a0(DIGT::get_allele(print_gt, 0));
                const uint8_t a1(DIGT::get_allele(print_gt, 1));
                uint8_t alt(a0);
                if (refBaseId == a0)
                {
                    alt = a1;
                }
                os << siteSampleInfo.dgt.phredLoghood[refBaseId] << ','
                   << siteSampleInfo.dgt.phredLoghood[DIGT::get_gt_with_alleles(refBaseId, alt)] << ','
                   << siteSampleInfo.dgt.phredLoghood[DIGT::get_gt_with_alleles(alt, alt)];
            }
        }
    }

    os << '\n';
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
            os << id_to_base(siteAlleles[altAlleleIndex].baseId);
        }
    }

}



void
gvcf_writer::
write_site_record(
    const GermlineContinuousSiteLocusInfo& locus) const
{
    const auto refBaseId = base_to_id(locus.ref);

    const unsigned sampleCount(locus.getSampleCount());
    const auto& siteAlleles(locus.getSiteAlleles());

    /// TODO STREL-125 tmp transitional structures:
    const bool is_no_alt(siteAlleles.empty());
    std::vector<uint8_t> altOrder;
    for (const auto& allele : siteAlleles)
    {
        assert(allele.baseId != refBaseId);
        altOrder.push_back(allele.baseId);
    }

    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (locus.pos+1) << '\t'  // POS
       << ".\t";           // ID

    os  << locus.ref << '\t'; // REF

    // ALT
    writeSiteVcfAltField(locus.getSiteAlleles(), os);
    os << '\t';

    // QUAL:
    os << locus.anyVariantAlleleQuality << '\t';

    // FILTER:
    getExtendedLocusFilters(locus).write(os);
    os << '\t';

    // INFO
    std::ostringstream info;

    if (locus.isVariantLocus())
    {
        {
            assert(not siteAlleles.empty());

            /// TODO STREL-125 generalize to multiple alts
            const auto& allele(siteAlleles.front());
            info << "SNVSB=";
            {
                const StreamScoper ss(info);
                info << std::fixed << std::setprecision(1) << allele.strandBias;
            }
            info << ';';
        }
        info << "SNVHPOL=" << locus.hpol;

        // compute global MQ over all samples
        MapqTracker mapqTracker;
        {
            for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
            {
                const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));
                mapqTracker.merge(siteSampleInfo.mapqTracker);
            }
        }
        info << ';';
        info << "MQ=" << mapqTracker.getRMS();

    }

    os << (info.str().empty() ? "." : info.str()) << "\t";

    //FORMAT
    os << "GT";
    os << ":GQ";
    os << ":GQX";
    os << ":DP:DPF";
    if (!is_no_alt)
    {
        os << ":AD:ADF:ADR";
    }
    os << "FT:VF";

    //SAMPLE
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& sampleInfo(locus.getSample(sampleIndex));
        const auto& siteSampleInfo(locus.getSiteSample(sampleIndex));

        os << '\t';

        VcfGenotypeUtil::writeGenotype(sampleInfo.getPloidy().getPloidy(),sampleInfo.max_gt(),os);

        //SAMPLE
        os << ':' << sampleInfo.genotypeQualityPolymorphic
           << ':' << sampleInfo.gqx;

        // DP:DPF
        os << ':' << siteSampleInfo.n_used_calls << ':' << siteSampleInfo.n_unused_calls;

        if (!is_no_alt)
        {
            os << ':';
            print_site_ad(refBaseId, siteSampleInfo, altOrder, os);
            os << ':';
            print_site_ad_strand(refBaseId, siteSampleInfo, altOrder, true, os);
            os << ':';
            print_site_ad_strand(refBaseId, siteSampleInfo, altOrder, false, os);
        }

        // FT
        os << ':';
        sampleInfo.filters.write(os);

        // VF
        {
            const auto& continuousSiteSampleInfo(locus.getContinuousSiteSample(sampleIndex));
            const StreamScoper ss(os);
            os << ':' << std::fixed << std::setprecision(3) << continuousSiteSampleInfo.getContinuousAlleleFrequency();
        }
    }
    os << '\n';
}



void
gvcf_writer::
write_site_record(
    const gvcf_block_site_record& locus) const
{
    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (locus.pos+1) << '\t'  // POS
       << ".\t";           // ID

    os  << locus.ref << '\t'; // REF

    // ALT
    os << '.';
    os << '\t';

    // QUAL:
    os << '.';
    os << '\t';

    // FILTER:
    locus.filters.write(os);
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
    os << locus.get_gt() << ':';
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
    for (unsigned strandIndex(0); strandIndex<2; ++strandIndex)
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
write_indel_record(
    const GermlineIndelLocusInfo& locus) const
{
    std::ostream& os(*_osptr);

    const unsigned sampleCount(locus.getSampleCount());

    // create VCF specific transformation of the alt allele list
    const auto& indelAlleles(locus.getIndelAlleles());
    OrthogonalAlleleSetLocusReportInfo locusReportInfo;
    getLocusReportInfoFromAlleles(_ref, indelAlleles, locusReportInfo);

    os << _chrom << '\t'   // CHROM
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
    os << "MQ=" << mapqTracker.getRMS();

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
            const auto& sampleInfo(locus.getSample(sampleIndex));
            const auto& indelSampleInfo(locus.getIndelSample(sampleIndex));

            os << '\t';

            VcfGenotypeUtil::writeGenotype(sampleInfo.getPloidy().getPloidy(),sampleInfo.max_gt(),os);
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
            const auto& sampleInfo(locus.getSample(sampleIndex));
            const auto& indelSampleInfo(locus.getIndelSample(sampleIndex));

            os << '\t';

            //SAMPLE
            // print GTs using a fake ploidy of 2, real ploidy is continuous...
            static const int printGTPloidy(2);
            VcfGenotypeUtil::writeGenotype(printGTPloidy,sampleInfo.max_gt(),os);
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

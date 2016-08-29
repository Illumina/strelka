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

    variant_prefilter_stage::add_site_modifiers(_empty_site, cm);
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



void gvcf_writer::filter_site_by_last_indel_overlap(GermlineDiploidSiteLocusInfo& si)
{
    if (_last_indel)
    {
        if (si.pos >= _last_indel->end())
        {
            _last_indel.reset(nullptr);
        }
        else
        {
            indel_overlapper::modify_overlapping_site(*_last_indel, si, _CM);
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


void gvcf_writer::process(std::unique_ptr<GermlineSiteLocusInfo> si)
{
    assert(si->getSampleCount() == getSampleCount());
    
    skip_to_pos(si->pos);

    if (dynamic_cast<GermlineDiploidSiteLocusInfo*>(si.get()) != nullptr)
    {
        add_site_internal(*downcast<GermlineDiploidSiteLocusInfo>(std::move(si)));
    }
    else
    {
        add_site_internal(*downcast<GermlineContinuousSiteLocusInfo>(std::move(si)));
    }

}



void
gvcf_writer::
process(std::unique_ptr<GermlineIndelLocusInfo> ii)
{
    skip_to_pos(ii->pos);

    // flush any non-variant block before starting:
    writeAllNonVariantBlockRecords();

    write_indel_record(*ii);
    if (dynamic_cast<GermlineDiploidIndelLocusInfo*>(ii.get()) != nullptr)
    {
        _last_indel = downcast<GermlineDiploidIndelLocusInfo>(std::move(ii));
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
    GermlineDiploidSiteLocusInfo& si)
{
    filter_site_by_last_indel_overlap(si);
    if (si.allele.is_phased_region)
    {
        _head_pos=si.pos+si.phased_ref.length();
    }
    else
    {
        _head_pos=si.pos+1;
    }
    // write_site
    queue_site_record(si);
}

void
gvcf_writer::
add_site_internal(
    GermlineContinuousSiteLocusInfo& si)
{
    // TODO: phasing
    _head_pos=si.pos+1;
    // write_site
    queue_site_record(si);
}




static
void
get_visible_alt_order(
    const GermlineDiploidSiteLocusInfo& si,
    std::vector<uint8_t>& altOrder)
{
    altOrder.clear();

    // list max_gt alts first:
    for (unsigned b(0); b<N_BASE; ++b)
    {
        if (b==si.dgt.ref_gt) continue;
        if (! DIGT::expect2(b,si.allele.max_gt)) continue;
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
    const GermlineSiteLocusInfo& si,
    const std::vector<uint8_t>& altOrder,
    std::ostream& os)
{
    os << si.alleleObservationCounts(base_to_id(si.ref));

    for (const auto& b : altOrder)
    {
        os << ',' << si.alleleObservationCounts(b);
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
    const GermlineSiteLocusInfo& si,
    const std::vector<uint8_t>& altOrder,
    const bool is_fwd_strand,
    std::ostream& os)
{
    os << si.alleleObservationCountsByStrand(is_fwd_strand, base_to_id(si.ref));

    for (const auto& b : altOrder)
    {
        os << ',' << si.alleleObservationCountsByStrand(is_fwd_strand,b);
    }
}



//writes out a SNP or block record
void
gvcf_writer::
write_site_record(
    const GermlineDiploidSiteLocusInfo& si) const
{
    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (si.pos+1) << '\t'  // POS
       << ".\t";           // ID

    if (si.allele.is_phased_region)
    {
        os  << si.phased_ref << '\t'; // REF
    }
    else
    {
        os  << si.ref << '\t'; // REF
    }

    // ALT
    std::vector<uint8_t> altOrder;
    const bool isNoAlt(si.allele.is_unknown);
    if (isNoAlt)
    {
        os << '.';
    }
    else if (si.allele.is_phased_region)
    {
        os << si.phased_alt;
    }
    else
    {
        get_visible_alt_order(si,altOrder);
        print_vcf_alt(altOrder,os);
    }
    os << '\t';

    // QUAL:
    if (si.is_qual())
    {
        os << si.dgt.genome.snp_qphred;
    }
    else
    {
        os << '.';
    }
    os << '\t';

    // FILTER:
    getExtendedLocusFilters(si).write(os);
    os << '\t';

    // INFO:
    if (si.dgt.is_snp)
    {
        os << "SNVSB=";
        {
            const StreamScoper ss(os);
            os << std::fixed << std::setprecision(1) << si.allele.strand_bias;
        }
        os << ';';
        os << "SNVHPOL=" << si.hpol;
        if (_opt.is_compute_hapscore)
        {
            os << ';';
            os << "HaplotypeScore=" << si.hapscore;
        }

        if (_opt.isReportEVSFeatures)
        {
#ifdef SUPPORT_LEGACY_EVS_TRAINING_SCRIPTS
            os << ';';
            os << "MQ=" << si.mapqRMS;
            os << ';';
            os << "MQ0=" << si.mapqZeroCount;
            os << ';';
            os << "MQRankSum=" << si.MQRankSum;
            os << ';';
            os << "BaseQRankSum=" << si.BaseQRankSum;
            os << ';';
            os << "ReadPosRankSum=" << si.ReadPosRankSum;
            os << ';';
            os << "AvgBaseQ=" << si.avgBaseQ;
            os << ';';
            os << "AvgPos=" << si.rawPos;
            // if you uncomment the following, make sure you also uncomment the matching INFO header entry in gvcf_header.cpp
            //                os << ';';
            //                os << "MapQ0Count=" << si.mapq_zero;

            // N.B. DP is in FORMAT already, and that seems to be where Nondas's code expects to find it, so suppress it here:
            //                os << ';';
            //                os << "DP=" << (si.n_used_calls+si.n_unused_calls);
#endif

            // EVS features may not be computed for certain records, so check first:
            if (! si.evsFeatures.empty())
            {
                const StreamScoper ss(os);
                os << std::setprecision(5);
                os << ";EVSF=";
                si.evsFeatures.writeValues(os);
                os << ",";
                si.evsDevelopmentFeatures.writeValues(os);
            }
        }

        if (si.allele.is_phasing_insufficient_depth)
        {
            os << ";Unphased";
        }
    }
    else
    {
        os << '.';
    }
    os << '\t';

    const bool is_nonref_gt(si.allele.max_gt != si.dgt.ref_gt);
    const bool is_print_pl(is_nonref_gt || si.dgt.is_snp);

    //FORMAT
    os << "GT";
    if (si.dgt.is_snp)
    {
        os << ":GQ";
    }
    os << ":GQX:DP:DPF";
    if (! isNoAlt)
    {
        os << ":AD:ADF:ADR";
    }
    os << ":FT";
    if (is_print_pl)
    {
        os << ":PL";
    }

    //SAMPLE
    const unsigned sampleCount(si.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& sampleInfo(si.getSample(sampleIndex));

        os << '\t';

        os << si.get_gt() << ':';
        if (si.dgt.is_snp)
        {
            os << sampleInfo.genotypeQualityPolymorphic << ':';
        }
        if (si.allele.is_gqx())
        {
            os << ((sampleInfo.empiricalVariantScore >= 0) ? sampleInfo.empiricalVariantScore : sampleInfo.gqx);
        }
        else
        {
            os << '.';
        }
        os << ':';
        //print DP:DPF
        os << si.n_used_calls << ':'
           << si.n_unused_calls;

        if (isNoAlt)
        {
            // pass
        }
        else if (si.allele.is_phased_region)
        {
            os << ':' << si.phased_AD
               << ':' << si.phased_ADF
               << ':' << si.phased_ADR;
        }
        else
        {
            os << ':';
            print_site_ad(si, altOrder, os);
            os << ':';
            print_site_ad_strand(si, altOrder, true, os);
            os << ':';
            print_site_ad_strand(si, altOrder, false, os);
        }

        // FT
        os << ':';
        sampleInfo.filters.write(os);

        if (is_print_pl)
        {
            // print PL values
            os << ':';
            if (si.is_hetalt())
            {
                const unsigned print_gt(si.allele.max_gt);
                const uint8_t a0(DIGT::get_allele(print_gt, 0));
                const uint8_t a1(DIGT::get_allele(print_gt, 1));
                os << si.dgt.phredLoghood[si.dgt.ref_gt] << ','
                   << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(si.dgt.ref_gt, a0)] << ','
                   << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(a0, a0)] << ','
                   << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(si.dgt.ref_gt, a1)] << ','
                   << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(a0, a1)] << ','
                   << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(a1, a1)];
            }
            else if (si.dgt.is_haploid() || (si.allele.modified_gt == MODIFIED_SITE_GT::ONE))
            {
                os << si.dgt.phredLoghood[si.dgt.ref_gt] << ','
                   << si.dgt.phredLoghood[si.allele.max_gt];
            }
            else
            {
                const unsigned print_gt(si.allele.max_gt);
                const uint8_t a0(DIGT::get_allele(print_gt, 0));
                const uint8_t a1(DIGT::get_allele(print_gt, 1));
                uint8_t alt(a0);
                if (si.dgt.ref_gt == a0)
                {
                    alt = a1;
                }
                os << si.dgt.phredLoghood[si.dgt.ref_gt] << ','
                   << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(si.dgt.ref_gt, alt)] << ','
                   << si.dgt.phredLoghood[DIGT::get_gt_with_alleles(alt, alt)];
            }
        }
    }

    os << '\n';
}



void
gvcf_writer::
write_site_record(
    const GermlineContinuousSiteLocusInfo& si) const
{
    auto ref_base_id = base_to_id(si.ref);

    assert(si.altAlleles.size() == 1);
    const auto& allele(si.altAlleles.front());

    std::vector<uint8_t> altOrder;
    const bool is_no_alt(allele._base == ref_base_id);
    if (! is_no_alt)
    {
        altOrder.push_back(allele._base);
    }

    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (si.pos+1) << '\t'  // POS
       << ".\t";           // ID

    os  << si.ref << '\t'; // REF

    std::string gt(si.get_gt(allele));

    // ALT
    if (is_no_alt)
        os << ".";
    else
        os << id_to_base(allele._base);
    os << '\t';

    // QUAL:
    os << si.anyVariantAlleleQuality << '\t';

    // FILTER:
    getExtendedLocusFilters(si).write(os);
    os << '\t';

    // INFO
    std::ostringstream info;

    if (si._is_snp)
    {
        info << "SNVSB=";
        {
            const StreamScoper ss(info);
            info << std::fixed << std::setprecision(1) << allele.strand_bias;
        }
        info << ';';
        info << "SNVHPOL=" << si.hpol;
    }

    if (!is_no_alt)
    {
        if (_opt.do_codon_phasing)
        {
            if (!info.str().empty())
                info << ";";
            info << "Unphased"; // TODO: placeholder until we do phasing on continuous variants
        }
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
    const unsigned sampleCount(si.getSampleCount());
    for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
    {
        const auto& sampleInfo(si.getSample(sampleIndex));

        os << '\t';

        //SAMPLE
        os << gt
           << ':' << sampleInfo.genotypeQualityPolymorphic
           << ':' << sampleInfo.gqx;

        // DP:DPF
        os << ':' << si.n_used_calls << ':' << si.n_unused_calls;

        if (!is_no_alt)
        {
            os << ':';
            print_site_ad(si, altOrder, os);
            os << ':';
            print_site_ad_strand(si, altOrder, true, os);
            os << ':';
            print_site_ad_strand(si, altOrder, false, os);
        }

        // FT
        os << ':';
        sampleInfo.filters.write(os);

        // VF
        {
            const StreamScoper ss(os);
            os << ':' << std::fixed << std::setprecision(3) << allele.variant_frequency();
        }
    }
    os << '\n';
}



void
gvcf_writer::
write_site_record(
    const gvcf_block_site_record& si) const
{
    std::ostream& os(*_osptr);

    os << _chrom << '\t'  // CHROM
       << (si.pos+1) << '\t'  // POS
       << ".\t";           // ID

    os  << si.ref << '\t'; // REF

    // ALT
    os << '.';
    os << '\t';

    // QUAL:
    os << '.';
    os << '\t';

    // FILTER:
    si.filters.write(os);
    os << '\t';

    // INFO:
    if (si.count>1)
    {
        os << "END=" << (si.pos+si.count) << ';';
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
    os << si.get_gt() << ':';
    if (si.has_call)
    {
        os << si.block_gqx.min();
    }
    else
    {
        os << '.';
    }
    os << ':';
    //print DP:DPF
    os << si.block_dpu.min() << ':'
       << si.block_dpf.min();
    os << '\n';
}



void
gvcf_writer::
write_indel_record(
    const GermlineIndelLocusInfo& ii) const
{
    std::ostream& os(*_osptr);

    // create VCF specific transformation of the alt allele list
    const auto& indelAlleles(ii.getIndelAlleles());
    OrthogonalAlleleSetLocusReportInfo locusReportInfo;
    getLocusReportInfoFromAlleles(_ref, indelAlleles, locusReportInfo);

    os << _chrom << '\t'   // CHROM
       << locusReportInfo.vcfPos << '\t'   // POS
       << ".\t"            // ID
       << locusReportInfo.vcfRefSeq << '\t'; // REF

    // ALT
    const unsigned altAlleleCount(ii.getAltAlleleCount());
    const unsigned fullAlleleCount(altAlleleCount+1);

    for (unsigned altAlleleIndex(0); altAlleleIndex < altAlleleCount; ++altAlleleIndex)
    {
        if (altAlleleIndex > 0) os << ',';
        os << locusReportInfo.altAlleles[altAlleleIndex].vcfAltSeq;
    }
    os << '\t';

    os << ii.anyVariantAlleleQuality << '\t'; //QUAL

    // FILTER:
    getExtendedLocusFilters(ii).write(os);
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

    const GermlineDiploidIndelLocusInfo* iiDiploidPtr(dynamic_cast<const GermlineDiploidIndelLocusInfo*>(&ii));
    if (iiDiploidPtr != nullptr)
    {
        const GermlineDiploidIndelLocusInfo& iiDiploid(*iiDiploidPtr);

        //FORMAT
        if (_opt.isReportEVSFeatures)
        {
            // EVS features may not be computed for certain records, so check first:
            if (! iiDiploid.evsFeatures.empty())
            {
                const StreamScoper ss(os);
                os << std::setprecision(5);
                os << ";EVSF=";
                iiDiploid.evsFeatures.writeValues(os);
                os << ",";
                iiDiploid.evsDevelopmentFeatures.writeValues(os);
            }
        }

        os << '\t';

        //FORMAT
        os << "GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL";

        //SAMPLE
        const unsigned sampleCount(ii.getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const auto& sampleInfo(ii.getSample(sampleIndex));
            const auto& indelSampleInfo(ii.getIndelSample(sampleIndex));
            const auto& sampleReportInfo(indelSampleInfo.reportInfo);

            os << '\t';

            VcfGenotypeUtil::writeGenotype(sampleInfo.getPloidy().getPloidy(),sampleInfo.max_gt(),os);
            os << ':' << sampleInfo.genotypeQualityPolymorphic;

            os << ':' << ((sampleInfo.empiricalVariantScore >= 0) ? sampleInfo.empiricalVariantScore : sampleInfo.gqx);

            os << ':' << sampleReportInfo.tier1Depth;

            // SAMPLE AD/ADF/ADR:
            {
                const auto& counts(sampleInfo.supportCounts);

                // verify locus and sample allele counts are in sync:
                assert(counts.getAltCount() == altAlleleCount);

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
        os << "GT:GQ:GQX:DPI:AD:ADF:ADR:VF";

        //SAMPLE
        const unsigned sampleCount(ii.getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            const auto& sampleInfo(ii.getSample(sampleIndex));
            const auto& indelSampleInfo(ii.getIndelSample(sampleIndex));
            const auto& sampleReportInfo(indelSampleInfo.reportInfo);

            os << '\t';

            //SAMPLE
            // print GTs using a fake ploidy of 2, real ploidy is continuous...
            static const int printGTPloidy(2);
            VcfGenotypeUtil::writeGenotype(printGTPloidy,sampleInfo.max_gt(),os);
            os << ':' << sampleInfo.genotypeQualityPolymorphic;

            os << ':' << sampleInfo.gqx;

            os << ':' << sampleReportInfo.tier1Depth;

            // AD:
            os << ':' << sampleReportInfo.n_confident_ref_reads
               << ',' << sampleReportInfo.n_confident_indel_reads;

            // ADF
            os << ':' << sampleReportInfo.n_confident_ref_reads_fwd
               << ',' << sampleReportInfo.n_confident_indel_reads_fwd;

            // ADR
            os << ':' << sampleReportInfo.n_confident_ref_reads_rev
               << ',' << sampleReportInfo.n_confident_indel_reads_rev;

            // VF
            {
                const StreamScoper ss(os);
                os << ':' << std::setprecision(3) << indelSampleInfo.alleleFrequency();
            }
        }
    }

    os << '\n';
}

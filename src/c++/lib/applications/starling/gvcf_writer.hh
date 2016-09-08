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

#pragma once

#include "gvcf_block_site_record.hh"
#include "gvcf_compressor.hh"
#include "starling_shared.hh"
#include "variant_pipe_stage_base.hh"

#include "blt_util/RegionTracker.hh"

#include <iosfwd>

struct ScoringModelManager;


/// Assembles all site and indel call information into a consistent set, blocks output
/// and writes to a VCF stream
///
struct gvcf_writer : public variant_pipe_stage_base
{
    gvcf_writer(
        const starling_options& opt,
        const starling_deriv_options& dopt,
        const reference_contig_segment& ref,
        const RegionTracker& nocompress_regions,
        const std::vector<std::string>& sampleNames,
        std::ostream* os,
        const ScoringModelManager& cm);

    void process(std::unique_ptr<GermlineSiteLocusInfo>) override;
    void process(std::unique_ptr<GermlineIndelLocusInfo>) override;

private:
    unsigned
    getSampleCount() const
    {
        return _blockPerSample.size();
    }
    void flush_impl() override;

    void add_site_internal(GermlineDiploidSiteLocusInfo& locus);
    void add_site_internal(GermlineContinuousSiteLocusInfo& locus);

    /// write out compressed non-variant block for one sample
    ///
    /// this function is used because non-variant block boundaries
    /// are independent for each sample and written to separate gVCF files
    ///
    void
    writeSampleNonVariantBlockRecord(const unsigned sampleIndex);

    /// write out compressed non-variant block for all samples
    ///
    /// frequenty the non-variant block ending criteria is shared by
    /// all samples (eg. introduction of a variant at the next position)
    /// so we frequently need to syncronize block end across all samples
    ///
    void
    writeAllNonVariantBlockRecords()
    {
        const unsigned sampleCount(getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            writeSampleNonVariantBlockRecord(sampleIndex);
        }
    }

    // queue site record for writing, after
    // possibly joining it into a compressed non-variant block
    //
    template<class TSiteInfo>
    void queue_site_record(const TSiteInfo& si)
    {
        //test for basic blocking criteria
        if (! _gvcf_comp.is_site_compressable(si))
        {
            writeAllNonVariantBlockRecords();
            write_site_record(si);
            return;
        }

        const unsigned sampleCount(getSampleCount());
        for (unsigned sampleIndex(0); sampleIndex<sampleCount; ++sampleIndex)
        {
            gvcf_block_site_record& block(_blockPerSample[sampleIndex]);
            if (! block.testCanSiteJoinSampleBlock(si, sampleIndex))
            {
                writeSampleNonVariantBlockRecord(sampleIndex);
            }
            block.joinSiteToSampleBlock(si, sampleIndex);
        }
    }

    void write_site_record(const GermlineDiploidSiteLocusInfo& locus) const;
    void write_site_record(const GermlineContinuousSiteLocusInfo& locus) const;
    void write_site_record(const gvcf_block_site_record& locus) const;

    void write_indel_record(const GermlineIndelLocusInfo& locus) const;

    /// fill in missing sites
    void skip_to_pos(const pos_t target_pos);

    const GermlineDiploidSiteLocusInfo&
    get_empty_site(const pos_t pos)
    {
        _empty_site.pos = pos;
        _empty_site.refBaseIndex = base_to_id(_ref.get_base(pos));
        _empty_site.isSiteUnphasable = true;
        return _empty_site;
    }

    const starling_options& _opt;
    const known_pos_range _report_range;
    const reference_contig_segment& _ref;
    std::ostream* _osptr;
    const char* _chrom;
    const gvcf_deriv_options _dopt;
    std::vector<gvcf_block_site_record> _blockPerSample;
    pos_t _head_pos;
    GermlineDiploidSiteLocusInfo _empty_site;

    std::unique_ptr<GermlineDiploidIndelLocusInfo> _last_indel;

    void filter_site_by_last_indel_overlap(GermlineDiploidSiteLocusInfo& locus);

    gvcf_compressor _gvcf_comp;
    const ScoringModelManager& _CM;

    /// print output limits:
    const unsigned maxPL = 999;
};

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

#pragma once

#include "gvcf_block_site_record.hh"
#include "gvcf_compressor.hh"
#include "ScoringModelManager.hh"
#include "starling_shared.hh"
#include "starling_streams.hh"
#include "variant_pipe_stage_base.hh"

#include "blt_util/RegionTracker.hh"

#include <iosfwd>


/// Assembles all site and indel call information into a consistent set, blocks output
/// and writes to a VCF stream
///
struct gvcf_writer : public variant_pipe_stage_base
{
    gvcf_writer(
        const starling_options& opt,
        const starling_deriv_options& dopt,
        const starling_streams& streams,
        const reference_contig_segment& ref,
        const RegionTracker& nocompressRegions,
        const RegionTracker& callRegions,
        const ScoringModelManager& scoringModels);

    void process(std::unique_ptr<GermlineSiteLocusInfo>) override;
    void process(std::unique_ptr<GermlineIndelLocusInfo>) override;

    void
    resetRegion(
        const std::string& chromName,
        const known_pos_range2& reportRange)
    {
        _chromName = chromName;
        _reportRange = reportRange;
        _headPos = _reportRange.begin_pos();
    }

private:
    unsigned
    getSampleCount() const
    {
        return _blockPerSample.size();
    }

    const std::string&
    getChromName() const
    {
        assert(not _chromName.empty());
        return _chromName;
    }

    void flush_impl() override;

    /// Add sites to queue for writing to gVCF
    void add_site_internal(GermlineSiteLocusInfo& locus);

    /// write out compressed non-variant block for one sample
    ///
    /// this function is used because non-variant block boundaries
    /// are independent for each sample and written to separate gVCF files
    ///
    void
    writeSampleNonVariantBlockRecord(const unsigned sampleIndex);

    /// Write out compressed non-variant block for all samples
    ///
    /// Although compressed non-variant blocks can start and end independently in each sample,
    /// the block ending criteria is often synchronized across all samples (due to, eg. introduction
    /// of a variant at the next position). This method assists with the latter case.
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

    /// queue site record for writing, after
    /// possibly joining it into a compressed non-variant block
    ///
    void
    queue_site_record(
        const GermlineSiteLocusInfo& locus);

    /// Write site record out to a single VCF stream
    ///
    /// \param targetSampleIndex The sample index. This indicates the index of the sample-specific gVCF to write to, or
    ///                          if the value is less than 0, this signifies writing to the variants VCF.
    void
    write_site_record_instance(
        const GermlineSiteLocusInfo& locus,
        std::ostream& os,
        const int targetSampleIndex = -1) const;

    /// write site record out to all VCF streams
    void
    write_site_record(
        const GermlineSiteLocusInfo& locus) const;

    /// special write function for gvcf compressed non-reference site blocks
    void
    write_site_record(
        const gvcf_block_site_record& locus,
        std::ostream& os) const;

    /// \brief Write indel record out to a single VCF stream
    ///
    /// \param targetSampleIndex The sample index. This indicates the index of the sample-specific gVCF to write to, or
    ///                          if the value is less than 0, this signifies writing to the variants VCF.
    void
    write_indel_record_instance(
        const GermlineIndelLocusInfo& locus,
        std::ostream& os,
        const int targetSampleIndex = -1) const;

    /// write indel record out to all VCF streams
    void
    write_indel_record(
        const GermlineIndelLocusInfo& locus) const;

    /// fill in missing sites
    void skip_to_pos(const pos_t target_pos);

    const GermlineDiploidSiteLocusInfo&
    get_empty_site(const pos_t pos)
    {
        _empty_site.pos = pos;
        _empty_site.refBaseIndex = base_to_id(_ref.get_base(pos));
        return _empty_site;
    }

    /// Check if the site overlaps with the last variant indel written
    ///
    /// Note this should only impact empty sites, because non-empty sites have already been
    /// handled in VariantOverlapResolver
    void modifySiteForConsistencyWithUpstreamIndels(GermlineDiploidSiteLocusInfo& locus);

    const starling_options& _opt;
    const starling_streams& _streams;
    const reference_contig_segment& _ref;
    const gvcf_deriv_options _dopt;
    std::vector<gvcf_block_site_record> _blockPerSample;
    GermlineDiploidSiteLocusInfo _empty_site;

    const RegionTracker& _callRegions;

    std::string _chromName;
    known_pos_range2 _reportRange;
    pos_t _headPos;

    /// Hold a reference to the last variant indel output
    ///
    /// This is held so that empty sites overlapping the indel can be modified if required. Note that
    /// non-empty sites have already been modified for indel overlap upstream of this object in the
    /// gvcf writer pipeline, but empty sites are created by this object, and these empty sites may
    /// require modification.
    ///
    std::unique_ptr<GermlineIndelLocusInfo> _lastVariantIndelWritten;

    gvcf_compressor _gvcf_comp;
    const ScoringModelManager& _scoringModels;

    /// print output limits:
    const unsigned maxPL = 999;
};

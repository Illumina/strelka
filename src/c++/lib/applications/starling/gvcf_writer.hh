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
#include "starling_streams.hh"
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
        const starling_streams& streams,
        const reference_contig_segment& ref,
        const RegionTracker& nocompress_regions,
        const ScoringModelManager& scoringModels);

    void process(std::unique_ptr<GermlineSiteLocusInfo>) override;
    void process(std::unique_ptr<GermlineIndelLocusInfo>) override;

private:
    unsigned
    getSampleCount() const
    {
        return _blockPerSample.size();
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

    /// queue site record for writing, after
    /// possibly joining it into a compressed non-variant block
    ///
    void
    queue_site_record(
        const GermlineSiteLocusInfo& locus);

    /// write site record out to a single VCF stream
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

    /// write indel record out to a single VCF stream
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

    const starling_options& _opt;
    const starling_streams& _streams;
    const reference_contig_segment& _ref;
    const known_pos_range _report_range;
    const char* _chrom;
    const gvcf_deriv_options _dopt;
    std::vector<gvcf_block_site_record> _blockPerSample;
    pos_t _head_pos;
    GermlineDiploidSiteLocusInfo _empty_site;

    std::unique_ptr<GermlineIndelLocusInfo> _last_indel;

    /// TODO STREL-125 why can't we get rid of this? Indel overlapper should already be doing the same thing!
    void filter_site_by_last_indel_overlap(GermlineDiploidSiteLocusInfo& locus);

    gvcf_compressor _gvcf_comp;
    const ScoringModelManager& _scoringModels;

    /// print output limits:
    const unsigned maxPL = 999;
};

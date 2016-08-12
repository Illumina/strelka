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


///
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
    void flush_impl() override;

    void add_site_internal(GermlineDiploidSiteLocusInfo& si);
    void add_site_internal(GermlineContinuousSiteLocusInfo& si);
    void write_block_site_record();

    // queue site record for writing, after
    // possibly joining it into a compressed non-variant block
    //
    template<class TSiteInfo>
    void queue_site_record(const TSiteInfo& si)
    {
        //test for basic blocking criteria
        if (! _gvcf_comp.is_site_compressable(si))
        {
            write_block_site_record();
            write_site_record(si);
            return;
        }

        if (! _block.test(si))
        {
            write_block_site_record();
        }
        _block.join(si);
    }

    void write_site_record(const GermlineDiploidSiteLocusInfo& si) const;
    void write_site_record(const GermlineContinuousSiteLocusInfo& si) const;
    void write_site_record(const gvcf_block_site_record& si) const;

    void write_indel_record(const GermlineDiploidIndelLocusInfo& ii) const;
    void write_indel_record(const GermlineContinuousIndelLocusInfo& ii) const;

    /// fill in missing sites
    void skip_to_pos(const pos_t target_pos);

    const GermlineDiploidSiteLocusInfo&
    get_empty_site(const pos_t pos)
    {
        _empty_site.pos = pos;
        _empty_site.ref = _ref.get_base(pos);
        _empty_site.Unphasable = true;
        return _empty_site;
    }

    const starling_options& _opt;
    const known_pos_range _report_range;
    const reference_contig_segment& _ref;
    std::ostream* _osptr;
    const char* _chrom;
    const gvcf_deriv_options _dopt;
    gvcf_block_site_record _block;
    pos_t _head_pos;
    GermlineDiploidSiteLocusInfo _empty_site;

    std::unique_ptr<GermlineDiploidIndelLocusInfo> _last_indel;

    void filter_site_by_last_indel_overlap(GermlineDiploidSiteLocusInfo& si);

    gvcf_compressor _gvcf_comp;
    const ScoringModelManager& _CM;
};

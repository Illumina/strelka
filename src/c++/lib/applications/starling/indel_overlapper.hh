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
 *  Created on: Jun 3, 2015
 *      Author: jduddy
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "variant_pipe_stage_base.hh"

#include <iosfwd>

struct ScoringModelManager;


struct indel_overlapper : public variant_pipe_stage_base
{
    indel_overlapper(
        const ScoringModelManager& scoringModels,
        const reference_contig_segment& ref,
        std::shared_ptr<variant_pipe_stage_base> destination)
        : variant_pipe_stage_base(destination)
        , _scoringModels(scoringModels)
        , _ref(ref)
        , _indel_range(-1, -1)
    {
        // this component doesn't make any sense without a destination:
        assert(destination);
    }

    void process(std::unique_ptr<GermlineSiteLocusInfo> siteLocusPtr) override;
    void process(std::unique_ptr<GermlineIndelLocusInfo> indelLocusPtr) override;

    static
    void
    modify_overlapping_site(
        const GermlineIndelLocusInfo& indelLocus,
        GermlineDiploidSiteLocusInfo& siteLocus,
        const ScoringModelManager& model);

private:
    void flush_impl() override
    {
        // flush out accumulated sites & indels
        process_overlaps();
    }

    static void modify_indel_conflict_site(GermlineSiteLocusInfo& siteLocus);
    static void modify_indel_overlap_site(
        const GermlineIndelLocusInfo& indelLocus,
        GermlineDiploidSiteLocusInfo& siteLocus,
        const ScoringModelManager& model);

    void process_overlaps();
    void process_overlaps_impl();
    void modify_single_indel_record(GermlineIndelLocusInfo& indelLocus);


    void modify_conflict_indel_record();

    void dump(std::ostream& os) const;

    void
    clearBuffers()
    {
        _indel_buffer.clear();
        _nonvariant_indel_buffer.clear();
        _site_buffer.clear();
    }

    const ScoringModelManager& _scoringModels;
    const reference_contig_segment& _ref;
    known_pos_range2 _indel_range;

    typedef std::unique_ptr<GermlineIndelLocusInfo> indel_ptr_t;

    std::vector<indel_ptr_t> _indel_buffer;
    std::vector<indel_ptr_t> _nonvariant_indel_buffer;
    std::vector<std::unique_ptr<GermlineDiploidSiteLocusInfo>> _site_buffer;
};


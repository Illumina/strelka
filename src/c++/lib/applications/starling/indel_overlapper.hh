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

class ScoringModelManager;


class indel_overlapper : public variant_pipe_stage_base
{
public:
    indel_overlapper(
        const ScoringModelManager& model,
        const reference_contig_segment& ref,
        std::shared_ptr<variant_pipe_stage_base> destination);

    void process(std::unique_ptr<site_info> si) override;
    void process(std::unique_ptr<indel_info> ii) override;

    static void modify_overlapping_site(const digt_indel_info& ii, digt_site_info& si, const ScoringModelManager& model);

private:
    void flush_impl() override
    {
        // flush out accumulated sites & indels
        process_overlaps();
    }

    static void modify_indel_conflict_site(digt_site_info& si);
    static void modify_indel_overlap_site(const digt_indel_info& ii,
                                          const unsigned ploidy,
                                          digt_site_info& si,
                                          const ScoringModelManager& model);

    void process_overlaps();
    void process_overlaps_impl();
    void modify_single_indel_record(digt_indel_info& ii);
    void modify_conflict_indel_record();
    void modify_overlap_indel_record();

    void dump(std::ostream& os) const;

    const ScoringModelManager& _CM;
    const reference_contig_segment& _ref;
    pos_t _indel_end_pos;

    std::vector<std::unique_ptr<digt_indel_info>> _indel_buffer;
    std::vector<std::unique_ptr<digt_indel_info>> _nonvariant_indel_buffer;
    std::vector<std::unique_ptr<digt_site_info>> _site_buffer;
};


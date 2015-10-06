// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/*
 *  Created on: Jun 3, 2015
 *      Author: jduddy
 */

#pragma once
#include "gvcf_locus_info.hh"
#include "variant_pipe_stage_base.hh"

class calibration_models;


class indel_overlapper : public variant_pipe_stage_base
{
public:
    indel_overlapper(const calibration_models& model, const reference_contig_segment& ref, std::shared_ptr<variant_pipe_stage_base> destination);

    void process(std::unique_ptr<site_info> si) override;
    void process(std::unique_ptr<indel_info> ii) override;

    void flush() override;

    static void modify_overlapping_site(const digt_indel_info& ii, digt_site_info& si, const calibration_models& model);

private:
    const calibration_models& _CM;
    const reference_contig_segment& _ref;
    pos_t _indel_end_pos;

    static void modify_indel_conflict_site(digt_site_info& si);
    static void modify_indel_overlap_site(const digt_indel_info& ii,
            const unsigned ploidy,
            digt_site_info& si,
            const calibration_models& model);



    std::vector<std::unique_ptr<digt_indel_info>> _indel_buffer;
    std::vector<std::unique_ptr<digt_site_info>> _site_buffer;


    void process_overlaps();
    void modify_single_indel_record();
    void modify_conflict_indel_record();
    void modify_overlap_indel_record();
};


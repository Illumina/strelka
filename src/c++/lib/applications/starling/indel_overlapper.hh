// -*- mode: c++; indent-tabs-mode: nil; -*-
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
    indel_overlapper(const calibration_models& model, const reference_contig_segment& ref, variant_pipe_stage_base& destination);

    void process(site_info& si) override;
    void process(indel_info& ii) override;

    void flush() override;

    bool has_buffered_indels() { return !_indel_buffer.empty(); }

private:
    const calibration_models& _CM;
    const reference_contig_segment& _ref;
    pos_t _indel_end_pos;

    std::vector<indel_info> _indel_buffer;
    std::vector<site_info> _site_buffer;


    void process_overlaps();
    void modify_single_indel_record();
    void modify_conflict_indel_record();
    void modify_overlap_indel_record();
    static void modify_indel_conflict_site(site_info& si);
    void modify_indel_overlap_site(const indel_info& ii,
        const unsigned ploidy,
        site_info& si);


};


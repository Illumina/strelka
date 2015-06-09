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

///
/// \author Chris Saunders
///

#pragma once


#include "calibration_models.hh"
#include "codon_phaser.hh"
#include "gvcf_block_site_record.hh"
#include "gvcf_locus_info.hh"
#include "gvcf_compressor.hh"
#include <iosfwd>
#include "bedstreamprocessor.hh"
#include "variant_prefilter_stage.hh"
#include "indel_overlapper.hh"


///
/// Assembles all site and indel call information into a consistent set, blocks output
/// and writes to a VCF stream
///
struct gvcf_aggregator : public variant_pipe_stage_base
{
    gvcf_aggregator(
        const starling_options& opt,
        const starling_deriv_options& dopt,
        const reference_contig_segment& ref,
        const RegionTracker& nocompress_regions,
        std::ostream* os,
        const pos_basecall_buffer& bc_buff);

    ~gvcf_aggregator();

    /// provide the phasing block status, if true, pileup buffer is
    /// preserved until the block is completed
    bool is_phasing_block() const
    {
        return _codon_phaser.is_in_block();
    }

    void add_site(site_info& si);

    void add_indel(const pos_t pos, const indel_key ik,
                   const starling_diploid_indel_core& dindel,
                   const starling_indel_report_info& iri,
                   const starling_indel_sample_report_info& isri);
    void reset();

    void process(site_info&) override;
    void process(indel_info&) override;
    void flush();

    pos_t
    headPos() const
    {
        return _head_pos;
    }

private:
    void add_site_internal(const site_info& si);
    void write_block_site_record();
    void write_site_record(const site_info& si) const;
    void queue_site_record(const site_info& si);
    void write_indel_record(const indel_info& ii);

    /// fill in missing sites
    void skip_to_pos(const pos_t target_pos);

    const site_info& get_empty_site(const pos_t pos)
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
    site_info _empty_site;

    pos_t _last_hetalt_indel_end_pos;

    calibration_models _CM;
    gvcf_compressor _gvcf_comp;

    indel_overlapper _overlapper;
    Codon_phaser _codon_phaser;
    bed_stream_processor _targeted_region_processor;
    variant_prefilter_stage _head;


};


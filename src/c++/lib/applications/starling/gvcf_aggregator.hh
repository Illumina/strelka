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
#include "gvcf_writer.hh"


///
/// Assembles all site and indel call information into a consistent set, blocks output
/// and writes to a VCF stream
///
class gvcf_aggregator
{
public:
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
    pos_t headPos() const
    {
        return _writer.headPos();
    }

private:
    calibration_models _CM;

    gvcf_writer _writer;
    indel_overlapper _overlapper;
    Codon_phaser _codon_phaser;
    bed_stream_processor _targeted_region_processor;
    variant_prefilter_stage _head;


};


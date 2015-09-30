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
        return _codon_phaser && _codon_phaser->is_in_block();
    }

    void add_site(std::unique_ptr<site_info> si);

    void add_indel(std::unique_ptr<indel_info> info);
    void reset();

private:
    calibration_models _CM;

    std::shared_ptr<gvcf_writer> _writer;
    std::shared_ptr<Codon_phaser> _codon_phaser;
    std::shared_ptr<variant_pipe_stage_base> _head;


};


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

#pragma once


#include "codon_phaser.hh"
#include "gvcf_block_site_record.hh"
#include "gvcf_locus_info.hh"
#include "gvcf_compressor.hh"

#include <iosfwd>

#include "ScoringModelManager.hh"


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
        const std::string& sampleName,
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
    ScoringModelManager _CM;

    std::shared_ptr<Codon_phaser> _codon_phaser;
    std::shared_ptr<variant_pipe_stage_base> _head;
};


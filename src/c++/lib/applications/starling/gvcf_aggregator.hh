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
#include "ScoringModelManager.hh"
#include "starling_streams.hh"

#include <iosfwd>


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
        const starling_streams& streams,
        const reference_contig_segment& ref,
        const RegionTracker& nocompress_regions,
        const std::vector<std::reference_wrapper<const pos_basecall_buffer>>& basecallBuffers);

    ~gvcf_aggregator();

    /// provide the phasing block status, if true, pileup buffer is
    /// preserved until the block is completed
    bool is_phasing_block() const
    {
        return (_codon_phaser && _codon_phaser->isBuffer());
    }

    void add_site(std::unique_ptr<GermlineSiteLocusInfo> si);

    void add_indel(std::unique_ptr<GermlineIndelLocusInfo> info);
    void reset();

private:
    ScoringModelManager _scoringModels;

    std::shared_ptr<Codon_phaser> _codon_phaser;
    std::shared_ptr<variant_pipe_stage_base> _head;
};


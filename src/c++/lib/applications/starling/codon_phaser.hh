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
 *  Test for codon-phasing.
 *
 *  Created on: Aug 10, 2013
 *  Author: Morten Kallberg
 */

#pragma once

#include "gvcf_locus_info.hh"
#include "starling_shared.hh"
#include "variant_pipe_stage_base.hh"

#include "starling_common/pos_basecall_buffer.hh"

#include <iosfwd>


/// short-range phasing utility for het-snps
///
/// requires extended preservation of the pileup buffer so that it
/// can go back and recover phase information form a candidate phasing block
///
/// \TODO generally recognized development direction is to record some kind of
///       read id in SNP pileups and indel support info so that we can go back
///       and phase from the hets without having to keep the whole read buffer (and so
///       read filtration, etc. is an exact match to the pileup).
///       Will this be worth doing before we transition to a haplotype assembly model
///       for short-range phasing?
///
struct Codon_phaser : public variant_pipe_stage_base
{
    Codon_phaser(
        const starling_options& opt,
        const reference_contig_segment& ref,
        const std::vector<std::reference_wrapper<const pos_basecall_buffer>>& initBasecallBuffers,
        std::shared_ptr<variant_pipe_stage_base> destination)
        : variant_pipe_stage_base(destination),
          _opt(opt),
          _ref(ref),
          _basecallBuffers(initBasecallBuffers),
          _sampleBlocks(initBasecallBuffers.size())
    {}

    void process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr) override;

    void process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr) override;

    bool
    isBuffer() const
    {
        return (not _buffer.empty());
    }

private:
    void flush_impl() override;

    static
    bool
    is_phasable_locus(
        const LocusInfo& locus,
        const unsigned sampleIndex)
    {
        if (not locus.isVariantLocus()) return false;
        return locus.getSample(sampleIndex).max_gt().isHet();
    }

    pos_t
    bufferStartPos() const
    {
        if (_buffer.empty()) return -1;
        return _buffer.front()->pos;
    }

    /// trigger phasing cycle for a single sample+block
    void phaseRecords(const unsigned sampleIndex);

    /// dump buffer contents to sink and clear object
    void output_buffer();

    struct allele_observations
    {
        int
        count() const
        {
            return (fwd+rev);
        }

        int fwd=0;
        int rev=0;
    };

    typedef std::map<std::string, allele_observations> allele_map_t;

    /// holds short-term per-sample per-block read information used to conclude phasing status with block
    struct BlockReadObservations
    {
        /// total used and unused reads spanning phasing region
        int total_reads = 0;
        int total_reads_unused = 0;

        allele_map_t observations;
    };

    /// \returns false if block fails consistency checks for phasing
    static
    bool
    evaluatePhasingConsistency(
        const BlockReadObservations& blockObs,
        const std::array<std::pair<std::string,allele_observations>, 2>& max_alleles);

    /// fill in allele count info for once sample+block
    void
    collect_pileup_evidence(
        const unsigned sampleIndex,
        BlockReadObservations& blockObs);

    /// update sample records with phasing information for one sample+block
    void
    create_phased_record(
        const unsigned sampleIndex,
        const BlockReadObservations& blockObs);

    struct SampleBlockInfo
    {
        /// Are we currently in a phasing block?
        bool is_in_block() const
        {
            return start != -1;
        }

        unsigned get_length() const
        {
            return (end-start+1);
        }

        void
        clear()
        {
            start = -1;
            end = -1;
            het_count = 0;
        }

        /// position of the first and last added het site to block
        pos_t start = -1;
        pos_t end = -1;

        /// total hets observed in block
        unsigned het_count = 0;
    };

    const starling_options& _opt;
    const reference_contig_segment& _ref;

    /// pass along the relevant pileup buffer
    const std::vector<std::reference_wrapper<const pos_basecall_buffer>> _basecallBuffers;

    std::vector<SampleBlockInfo> _sampleBlocks;
    std::vector<std::unique_ptr<GermlineSiteLocusInfo>> _buffer;
};

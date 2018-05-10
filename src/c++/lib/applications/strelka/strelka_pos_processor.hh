//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

/// \file
/// \author Chris Saunders
///

#pragma once


#include "NoiseBuffer.hh"
#include "strelka_shared.hh"
#include "SomaticIndelVcfWriter.hh"
#include "strelka_streams.hh"

#include "starling_common/PileupCleaner.hh"
#include "starling_common/starling_pos_processor_base.hh"
#include "SomaticCallableProcessor.hh"
#include "strelka_common/StrelkaSampleSetSummary.hh"


///
///
struct strelka_pos_processor : public starling_pos_processor_base
{
    typedef starling_pos_processor_base base_t;

    strelka_pos_processor(
        const strelka_options& opt,
        const strelka_deriv_options& dopt,
        const reference_contig_segment& ref,
        const strelka_streams& fileStreams,
        RunStatsManager& statsManager);

    void reset() override;

    void
    resetRegion(
        const std::string& chromName,
        const known_pos_range2& regionRange);

    void
    insert_noise_pos(
        const pos_t pos,
        const SiteNoise& sn);

private:

    void
    process_pos_variants_impl(
        const pos_t pos,
        const bool isPosPrecedingReportableRange) override;

    void
    process_pos_snp_somatic(const pos_t pos);

    void
    process_pos_indel_somatic(const pos_t pos);

    void
    post_align_clear_pos(const pos_t pos) override
    {
        _noisePos.clear_pos(pos);
    }

    void
    run_post_call_step(
        const int stage_no,
        const pos_t pos) override;

    bool
    derived_empty() const override
    {
        return _indelWriter.empty();
    }

    /////////////////////////////

    // keep some of the original strelka classes handy so we don't
    // have to down-cast:
    const strelka_options& _opt;
    const strelka_deriv_options& _dopt;
    const strelka_streams& _streams;

    double _normChromDepth = 0.;
    double _maxChromDepth = 0.;

    CleanedPileup _tier2_cpi[STRELKA_SAMPLE_TYPE::SIZE];

    SomaticCallableProcessor _scallProcessor;

    // enables delayed indel write:
    SomaticIndelVcfWriter _indelWriter;

    unsigned _indelRegionIndexNormal;
    unsigned _indelRegionIndexTumor;

    NoiseBuffer _noisePos;
};

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
/// \author John Duddy

#pragma once

#include "gvcf_locus_info.hh"
#include "variant_pipe_stage_base.hh"

#include <iosfwd>

struct ScoringModelManager;


/// Component of the gvcf output pipeline which handles overlapping variant issues
///
/// This object has several roles:
/// 1) Handle overlapping indels: mark all overlapping indels with "IndelConflict" filter (rare)
/// 2) Change the ploidy of all site record which are overlapped by an indel (common)
/// 3) Reorder sets of overlapping sites and indels according to gvcf output ordering conventions (common)
///
struct VariantOverlapResolver : public variant_pipe_stage_base
{
    VariantOverlapResolver(
        const ScoringModelManager& scoringModels,
        std::shared_ptr<variant_pipe_stage_base> destination)
        : variant_pipe_stage_base(destination)
        , _scoringModels(scoringModels)
        , _bufferedVariantIndelRange(-1, -1)
    {
        // this component doesn't make any sense without a destination:
        assert(destination);
    }

    void process(std::unique_ptr<GermlineSiteLocusInfo> siteLocusPtr) override;
    void process(std::unique_ptr<GermlineIndelLocusInfo> indelLocusPtr) override;

    /// Adjust site record details for greater consistency with the overlapping indel
    ///
    /// The indel must be a variant and overlap the site.
    static
    void
    modifySiteOverlappingVariantIndel(
        const GermlineIndelLocusInfo& indelLocus,
        GermlineDiploidSiteLocusInfo& siteLocus,
        const ScoringModelManager& model);

private:
    void flush_impl() override
    {
        // flush out accumulated sites & indels
        processOverlappingVariants();
        clearBuffers();
    }

    /// Adjust site record details for greater consistency with the overlapping indel
    ///
    /// The indel must be a variant and overlap the site.
    static
    void
    modifySiteOverlappingNonconflictingVariantIndel(
        const GermlineIndelLocusInfo& indelLocus,
        GermlineDiploidSiteLocusInfo& siteLocus,
        const ScoringModelManager& model);

    /// Process all buffered variant information, resolving any overlaps and/or reordering the output as required
    ///
    /// This version of process_overlaps will add helpful information if an exception is thrown and prevent the class
    /// from becoming unstable.
    void processOverlappingVariants();

    /// This is a helper method for processOverlappingVariants, it should never be called directly by any other method
    void processOverlappingVariantsImplementation();

    /// Change all buffered variant indel records to reflect an indel conflict
    void annotateVariantIndelRecordsAsConflicting();

    void dump(std::ostream& os) const;

    void
    clearBuffers()
    {
        _bufferedVariantIndelRange.set_range(-1, -1);
        _variantIndelBuffer.clear();
        _nonvariantIndelBuffer.clear();
        _siteBuffer.clear();
    }

    const ScoringModelManager& _scoringModels;

    /// Union of buffered variant indel ranges
    known_pos_range2 _bufferedVariantIndelRange;

    typedef std::unique_ptr<GermlineIndelLocusInfo> indel_ptr_t;

    /// All buffered (overlapping) variant indels
    std::vector<indel_ptr_t> _variantIndelBuffer;

    /// All buffered non-variant indels
    ///
    /// Non-variants have a separate buffer because their overlap/conflict rules are substantially different
    /// than those for variants
    std::vector<indel_ptr_t> _nonvariantIndelBuffer;

    /// All buffered sites
    ///
    /// Sites are only buffered when they overlap at least one variant indel
    std::vector<std::unique_ptr<GermlineDiploidSiteLocusInfo>> _siteBuffer;
};


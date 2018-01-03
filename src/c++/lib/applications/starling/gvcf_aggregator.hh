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

#pragma once


#include "VariantPhaser.hh"
#include "gvcf_block_site_record.hh"
#include "gvcf_locus_info.hh"
#include "gvcf_compressor.hh"
#include "gvcf_writer.hh"
#include "ScoringModelManager.hh"
#include "starling_streams.hh"

#include <iosfwd>


/// Assembles all site and indel call information into a consistent set,
/// compresses sites into hom-ref blocks, and writes output to several VCF streams
///
class gvcf_aggregator
{
public:
    gvcf_aggregator(
        const starling_options& opt,
        const starling_deriv_options& dopt,
        const starling_streams& streams,
        const reference_contig_segment& ref,
        const RegionTracker& nocompressRegions,
        const RegionTracker& callRegions,
        const unsigned sampleCount);

    ~gvcf_aggregator();

    void add_site(std::unique_ptr<GermlineSiteLocusInfo> si);

    void add_indel(std::unique_ptr<GermlineIndelLocusInfo> info);
    void reset();

    void
    resetRegion(
        const std::string& chromName,
        const known_pos_range2& reportRegion)
    {
        _scoringModels.resetChrom(chromName);
        assert(_gvcfWriterPtr);
        _gvcfWriterPtr->resetRegion(chromName, reportRegion);
    }

    double
    getMaxDepth() const
    {
        return _scoringModels.getMaxDepth();
    }

private:
    ScoringModelManager _scoringModels;

    std::shared_ptr<VariantPhaser> _variantPhaserPtr;
    std::shared_ptr<gvcf_writer> _gvcfWriterPtr;
    std::shared_ptr<variant_pipe_stage_base> _head;
};

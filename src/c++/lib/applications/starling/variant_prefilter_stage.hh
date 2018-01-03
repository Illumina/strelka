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
/*
 *      Author: jduddy
 */

#pragma once

#include "variant_pipe_stage_base.hh"

struct RegionTracker;
struct ScoringModelManager;

struct variant_prefilter_stage : public variant_pipe_stage_base
{
    variant_prefilter_stage(
        const ScoringModelManager& model,
        const std::shared_ptr<variant_pipe_stage_base>& destination)
        : variant_pipe_stage_base(destination)
        , _model(model)
    {}


    void process(std::unique_ptr<GermlineSiteLocusInfo> locusPtr) override;
    void process(std::unique_ptr<GermlineIndelLocusInfo> locusPtr) override;

private:
    void
    applySharedLocusFilters(
        LocusInfo& locus) const;

    const ScoringModelManager& _model;
};

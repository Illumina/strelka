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
 *      Author: Morten Kallberg
 */

#pragma once

#include "VariantScoringModelBase.hh"
#include "VariantScoringModelMetadata.hh"
#include "VariantScoringModelTypes.hh"

#include <cmath>

#include <algorithm>
#include <memory>


/// \brief Client interface to variant scoring models specified by file at runtime
///
struct VariantScoringModelServer
{
    /// \param[in] featureMap Names of features supported in the client code, each feature
    ///                       name should be mapped to a feature index number.
    VariantScoringModelServer(
        const VariantScoringModelMetadata::featureMap_t& featureMap,
        const std::string& modelFile,
        const SCORING_CALL_TYPE::index_t callType,
        const SCORING_VARIANT_TYPE::index_t variantType);

    /// \return Probability that the variant call is false
    double
    scoreVariant(
        const VariantScoringModelBase::featureInput_t& features) const
    {
        return std::max(0.,std::min(1.,(_meta.probScale * std::pow(_model->getProb(features), _meta.probPow))));
    }

    double scoreFilterThreshold() const
    {
        return _meta.filterCutoff;
    }

private:
    VariantScoringModelMetadata _meta;
    std::unique_ptr<VariantScoringModelBase> _model;
};

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
 *      Author: Morten Kallberg
 */

#pragma once

#include "RandomForestModel.hh"
#include "VariantScoringModelMetadata.hh"
#include "VariantScoringModelTypes.hh"



/// client interface to variant scoring models specified by file at runtime
///
struct VariantScoringModel
{
    /// \param featureMap Names of features supported in the client code, each feature
    ///                   name should be mapped to a feature index number.
    VariantScoringModel(
        const featureMap_t& featureMap,
        const std::string& model_file,
        const SCORING_CALL_TYPE::index_t ctype,
        const SCORING_VARIANT_TYPE::index_t vtype);

    double
    scoreVariant(
        const feature_type& features) const
    {
        return _model.getProb(features);
    }

    double scoreFilterThreshold() const
    {
        return _meta.FilterCutoff;
    }

private:
    VariantScoringModelMetadata _meta;
    RandomForestModel _model;
};

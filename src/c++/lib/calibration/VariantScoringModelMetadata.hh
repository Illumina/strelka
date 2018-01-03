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

#include "rapidjson/document.h"

#include <map>
#include <string>


/// Parse common meta-data format shared for all variant scoring models
///
struct VariantScoringModelMetadata
{
    typedef std::map<std::string,unsigned> featureMap_t;

    VariantScoringModelMetadata() {}

    void Deserialize(
        const featureMap_t& featureMap,
        const rapidjson::Value& root);

    std::string date;
    std::string modelType;

    /// as part of an optional calibration component to all models, raise the raw prob by this power before reporting:
    double probPow = 1.;

    /// as part of an optional calibration component to all models, scale the prob using this factor (after power transform) before reporting:
    double probScale = 1.;

    /// Phred-scale threshold: PASS variants will be >= filterCutoff
    double filterCutoff;
};

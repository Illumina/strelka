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
 * SerializedModel.hh
 *
 *  Created on: Jun 23, 2015
 *      Author: Morten Kallberg
 */

#pragma once

#include "json/json.h"

#include <cassert>

#include <map>



namespace VARIATION_NODE_TYPE
{
enum index_t
{
    SNP,
    INDEL,
    SIZE
};

inline
const std::string
get_label(const index_t i)
{
    switch (i)
    {
    case SNP:
        return "SNP";
    case INDEL:
        return "INDEL";
    default:
        assert(false && "Unknown variation node type in serializedModel.");
        return nullptr;
    }
}
}


class serialized_model
{
public:
    serialized_model() {}

    /** methods for serializing */
    void Deserialize( const Json::Value& root);

    const std::string&
    get_model_string() const
    {
        return name;
    }

protected:
    std::string name;
    std::string version;
    std::string date;
};


typedef std::map<int, double> feature_type;
class serialized_calibration_model : public serialized_model
{
public:
    serialized_calibration_model() {}

    /** methods for serializing */
    void Deserialize( const Json::Value& root);
    double getProb(const feature_type& features) const;
    bool doFilter() const;

protected:
    std::string ModelType;
    std::string Type;
    double FilterCutoff;
    // add feature sequence
};


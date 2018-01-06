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

#include <iostream>
#include "common/Exceptions.hh"
#include "rapidjson/document.h"

class RapidJsonHelper
{
public:
    static const
    rapidjson::Value& getNodeMember(const rapidjson::Value& node, const char* label)
    {
        const rapidjson::Value::ConstMemberIterator iter(node.FindMember(label));
        return iter->value;
    }
    static
    void
    wrongValueTypeError(
        const char* key,
        const char* keyType)
    {
        std::ostringstream oss;
        oss << "Node '" << key << "' does not have expected type '" << keyType << "'.";
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }

    static
    void
    missingNodeError(
        const std::string& fileName,
        const char* key)
    {
        std::ostringstream oss;
        oss << "Can't find expected node '" << key << "' in json file '" << fileName << "'.";
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }

};


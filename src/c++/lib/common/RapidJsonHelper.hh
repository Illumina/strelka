//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

/**
 ** \file
 ** \brief Declaration of the common exception mechanism.
 **
 ** All exceptions must carry the same data (independently of the
 ** exception type) to homogenize the reporting and processing of
 ** errors.
 **
 ** \author Come Raczy
 **/

#pragma once

#include <iostream>
#include "common/Exceptions.hh"
#include "rapidjson/document.h"

namespace illumina
{
namespace common
{

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
        const std::string& fileName,
        const char* key,
        const char* keyType)
    {
        std::ostringstream oss;
        oss << "ERROR: Node '" << key << "' does not have expected type '" << keyType
            << "' in json file '" << fileName << "'.";
        BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
    }

    static
    void
    missingNodeError(
        const std::string& fileName,
        const char* key)
    {
        std::ostringstream oss;
        oss << "ERROR: Can't find expected node '" << key << "' in json file '" << fileName << "'.";
        BOOST_THROW_EXCEPTION(illumina::common::LogicException(oss.str()));
    }

};

}
}

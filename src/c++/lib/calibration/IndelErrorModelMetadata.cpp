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
/*
 *
 *      Author: mkallberg
 */

#include "IndelErrorModelMetadata.hh"

#include <cassert>


namespace IMODEL_ENTRY_TYPE
{
enum index_t
{
    NAME,
    VERSION,
    DATE,
    SIZE
};

static
const char*
get_label(const index_t i)
{
    switch (i)
    {
    case NAME:
        return "Name";
    case VERSION:
        return "Version";
    case DATE:
        return "Date";
    default:
        assert(false && "Unknown indel error model entry type");
        return nullptr;
    }
}
}



/// remove newline and quotes from string
static
std::string
Clean_string(const std::string& str)
{
    std::string temp = str;
    // strip any trailing newline
    std::size_t lastNewline = temp.find_last_of("\n\r");
    temp = temp.substr(0, lastNewline);

    // strip flanking double quotes if present
    std::size_t doubleQuote = temp.find_last_of("\"");
    temp = temp.substr(0, doubleQuote);

    doubleQuote = temp.find_first_of("\"");
    if (doubleQuote != std::string::npos)
    {
        temp = temp.substr(doubleQuote);
    }

    return temp;
}



void
IndelErrorModelMetadata::
deserialize(const Json::Value& root)
{
    using namespace IMODEL_ENTRY_TYPE;
    this->name      = Clean_string(root[get_label(NAME)].asString());
    this->version   = Clean_string(root[get_label(VERSION)].asString());
    this->date      = Clean_string(root[get_label(DATE)].asString());
}

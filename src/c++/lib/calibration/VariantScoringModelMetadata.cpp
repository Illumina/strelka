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
 *      Author: mkallberg
 */

#include "VariantScoringModelMetadata.hh"
#include "common/Exceptions.hh"

#include <cassert>


namespace SMODEL_ENTRY_TYPE
{
enum index_t
{
    NAME,
    VERSION,
    DATE,
    MODELTYPE,
    FILTERCUTOFF,
    FEATURES,
    CALIBRATION,
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
    case MODELTYPE:
        return "ModelType";
    case FILTERCUTOFF:
        return "FilterCutoff";
    case FEATURES:
        return "Features";
    case CALIBRATION:
        return "Calibration";
    default:
        assert(false && "Unknown serialized calibration model entry type");
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



static
void
featureMapError(
    const VariantScoringModelMetadata::featureMap_t& featureMap,
    const std::string& fname)
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: Can't find requested scoring model feature '" << fname << "' in candidate feature map containing:\n";
    for (const auto& val : featureMap)
    {
        oss << val.first << "\n";
    }
    oss << "\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



static
void
featureMapOrderError(
    const std::string& fname,
    const unsigned expectedIdx,
    const unsigned foundIdx)
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "ERROR: Scoring model feature '" << fname << "' is at position " << foundIdx << " in {calltype}VariantEmpiricalScoringFeatures.hh but at position " << expectedIdx << " in json model file. Check that feature order in the two sources correspond.\n";
    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
}



void
VariantScoringModelMetadata::
Deserialize(
    const featureMap_t& featureMap,
    const Json::Value& root)
{
    using namespace SMODEL_ENTRY_TYPE;
    date  = Clean_string(root[get_label(DATE)].asString());
    ModelType = root[get_label(MODELTYPE)].asString();
    filterCutoff = root[get_label(FILTERCUTOFF)].asDouble();

    // read optional calibration items:
    const Json::Value caliRoot = root[get_label(CALIBRATION)];
    if (not caliRoot.isNull())
    {
        probPow = caliRoot.get("Power", probPow).asDouble();
        probScale = caliRoot.get("Scale", probScale).asDouble();
    }

    // read and validate features:
    const Json::Value featureRoot = root[get_label(FEATURES)];
    assert(!featureRoot.isNull());

    const auto fend(featureMap.end());

    unsigned expectedIndex=0;
    for (const auto& val : featureRoot)
    {
        const std::string fname(val.asString());
        const auto fiter(featureMap.find(fname));
        if (fiter == fend)
        {
            featureMapError(featureMap,fname);
        }
        if (expectedIndex != fiter->second)
        {
            featureMapOrderError(fname,expectedIndex,fiter->second);
        }
        expectedIndex++;
    }

    if (expectedIndex != featureMap.size())
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: scoring feature count specified in modelfile (" << expectedIndex << ")"
            << " does not match expected feature count (" << featureMap.size() << ")\n";
        {
            bool isFirst(true);
            oss << "\tModelfile features: {";
            for (const auto& val : featureRoot)
            {
                if (not isFirst) oss << ",";
                oss << val.asString();
                isFirst=false;
            }
            oss << "}\n";
        }
        {
            bool isFirst(true);
            oss << "\tExpected features: {";
            for (const auto& val : featureMap)
            {
                if (not isFirst) oss << ",";
                oss << val.first;
                isFirst=false;
            }
            oss << "}\n";
        }
        oss << "\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
}


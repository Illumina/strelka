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



static
void
missingNodeError(
    const char* key)
{
    std::ostringstream oss;
    oss << "Can't find expected node '" << key << "' in  json scoring model file.";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
}



static
void
wrongValueTypeError(
    const char* key,
    const char* keyType)
{
    std::ostringstream oss;
    oss << "Node '" << key << "' in json scoring model file does not have expected type '" << keyType << "'.";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
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
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
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
    oss << "Scoring model feature '" << fname << "' is at position " << foundIdx
        << " in {calltype}VariantEmpiricalScoringFeatures.hh but at position " << expectedIdx
        << " in json model file. Check that feature order in the two sources correspond.";
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
}



void
VariantScoringModelMetadata::
Deserialize(
    const featureMap_t& featureMap,
    const rapidjson::Value& root)
{
    using namespace SMODEL_ENTRY_TYPE;

    auto getNodeMember = [](const rapidjson::Value& node, const char* label) -> const rapidjson::Value&
    {
        const rapidjson::Value::ConstMemberIterator iter(node.FindMember(label));
        if (iter == node.MemberEnd()) missingNodeError(label);
        return iter->value;
    };

    const rapidjson::Value& dateValue(getNodeMember(root, get_label(DATE)));
    if (! dateValue.IsString()) wrongValueTypeError(get_label(DATE), "string");
    date = Clean_string(dateValue.GetString());

    const rapidjson::Value& modelTypeValue(getNodeMember(root, get_label(MODELTYPE)));
    if (! modelTypeValue.IsString()) wrongValueTypeError(get_label(MODELTYPE), "string");
    modelType = modelTypeValue.GetString();

    const rapidjson::Value& filterCutoffValue(getNodeMember(root, get_label(FILTERCUTOFF)));
    if (! filterCutoffValue.IsNumber()) wrongValueTypeError(get_label(FILTERCUTOFF), "number");
    filterCutoff = filterCutoffValue.GetDouble();

    // read optional calibration items:
    const rapidjson::Value::ConstMemberIterator calibrationRootIter(root.FindMember(get_label(CALIBRATION)));
    if (calibrationRootIter != root.MemberEnd())
    {
        const rapidjson::Value& calibrationRoot(calibrationRootIter->value);
        if (! calibrationRoot.IsObject())
        {
            wrongValueTypeError(get_label(CALIBRATION), "object");
        }

        auto optionalDoubleKeyUpdate = [&](const char* label, double& val)
        {
            const rapidjson::Value::ConstMemberIterator iter(calibrationRoot.FindMember(label));
            if (iter != calibrationRoot.MemberEnd())
            {
                if (! iter->value.IsNumber()) wrongValueTypeError(label, "number");
                val = iter->value.GetDouble();
            }
        };

        optionalDoubleKeyUpdate("Power", probPow);
        optionalDoubleKeyUpdate("Scale", probScale);
    }

    // read and validate features:
    const rapidjson::Value& featureRoot(getNodeMember(root, get_label(FEATURES)));
    if (! featureRoot.IsArray()) wrongValueTypeError(get_label(FEATURES), "array");

    const auto fend(featureMap.end());

    unsigned expectedIndex=0;
    for (const auto& val : featureRoot.GetArray())
    {
        const std::string featureName(val.GetString());
        const auto fiter(featureMap.find(featureName));
        if (fiter == fend)
        {
            featureMapError(featureMap, featureName);
        }
        if (expectedIndex != fiter->second)
        {
            featureMapOrderError(featureName, expectedIndex, fiter->second);
        }
        expectedIndex++;
    }

    if (expectedIndex != featureMap.size())
    {
        std::ostringstream oss;
        oss << "ERROR: scoring feature count specified in modelfile (" << expectedIndex << ")"
            << " does not match expected feature count (" << featureMap.size() << ")\n";
        {
            bool isFirst(true);
            oss << "\tModelfile features: {";
            for (const auto& val : featureRoot.GetArray())
            {
                if (not isFirst) oss << ",";
                oss << val.GetString();
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
        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
    }
}


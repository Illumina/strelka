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

#include "common/RapidJsonHelper.hh"
#include "ThetaJson.hh"

ThetaJson::ThetaJson(
    size_t repeatPatternSize,
    const std::vector<double>& theta):
    _repeatPatternSize(repeatPatternSize),
    _theta(theta)
{
}

ThetaJson
ThetaJson::deserialize(
    const rapidjson::Value& root)
{
    using namespace illumina::common;

    static const char* repeatPatternSizeLabel = "repeatPatternSize";
    const rapidjson::Value& repeatPatternSizeValue(RapidJsonHelper::getNodeMember(root, repeatPatternSizeLabel));
    if (!repeatPatternSizeValue.IsUint()) RapidJsonHelper::wrongValueTypeError(repeatPatternSizeLabel, "unsigned");
    size_t repeatPatternSize(repeatPatternSizeValue.GetUint());

    std::vector<double> theta;
    static const char* thetaLabel = "theta";
    const rapidjson::Value& thetaArray(RapidJsonHelper::getNodeMember(root, thetaLabel));
    if (!thetaArray.IsArray()) RapidJsonHelper::wrongValueTypeError(thetaLabel, "array");
    for (const auto& thetaValue : thetaArray.GetArray())
    {
        if (!thetaValue.IsNumber()) RapidJsonHelper::wrongValueTypeError(thetaLabel, "number");
        theta.push_back(thetaValue.GetDouble());
    }

    return ThetaJson(repeatPatternSize, theta);
}

ThetasJson::ThetasJson(std::map<unsigned, std::vector<double>>& thetasMap)
    : _thetasMap(thetasMap)
{
}

ThetasJson
ThetasJson::deserialize(
    const rapidjson::Value& root)
{
    using namespace illumina::common;
    static const char* thetasLabel = "thetas";
    std::map<unsigned, std::vector<double>> thetasMap;
    const rapidjson::Value& thetasArray(RapidJsonHelper::getNodeMember(root, thetasLabel));
    if (!thetasArray.IsArray()) RapidJsonHelper::wrongValueTypeError(thetasLabel, "array");
    for (const auto& thetasByPatternSize : thetasArray.GetArray())
    {
        const ThetaJson theta(ThetaJson::deserialize(thetasByPatternSize));
        thetasMap[theta.getRepeatPatternSize()] = theta.getTheta();
    }

    return ThetasJson(thetasMap);
}

const std::map<unsigned, std::vector<double>>&
                                           ThetasJson::getThetasMap() const
{
    return _thetasMap;
}

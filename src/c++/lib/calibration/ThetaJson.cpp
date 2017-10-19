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

#include "ThetaJson.hh"


static const rapidjson::Value &getNodeMember(const rapidjson::Value& node, const char* label)
{
    const rapidjson::Value::ConstMemberIterator iter(node.FindMember(label));
    return iter->value;
}

ThetaJson::ThetaJson(size_t repeatPatternSize, const std::vector<double> theta):
_repeatPatternSize(repeatPatternSize),
_theta(theta)
{
}

ThetaJson ThetaJson::Deserialize(const rapidjson::Value &root)
{
    static const char* repeatPatternSizeLabel = "repeatPatternSize";
    const rapidjson::Value& repeatPatternSizeValue(getNodeMember(root, repeatPatternSizeLabel));
    size_t repeatPatternSize(repeatPatternSizeValue.GetUint());

    std::vector<double> theta;
    static const char* thetaLabel = "theta";
    const rapidjson::Value& thetaArray(getNodeMember(root, thetaLabel));

    for (const auto& thetaValue : thetaArray.GetArray())
    {
        theta.push_back(thetaValue.GetDouble());
    }

    return ThetaJson(repeatPatternSize, theta);
}

ThetasJson::ThetasJson(std::map<unsigned, std::vector<double>> thetasMap)
: _thetasMap(thetasMap)
{
}

ThetasJson ThetasJson::Deserialize(const rapidjson::Value &root)
{
    static const char* thetasLabel = "thetas";
    std::map<unsigned, std::vector<double>> thetasMap;
    const rapidjson::Value& thetasArray(getNodeMember(root, thetasLabel));
    for (const auto& thetasByPatternSize : thetasArray.GetArray())
    {
        const ThetaJson theta(ThetaJson::Deserialize(thetasByPatternSize));
        thetasMap[theta.RepeatPatternSize()] = theta.Theta();
    }

    return ThetasJson(thetasMap);
}

const std::map<unsigned, std::vector<double>> &ThetasJson::ThetasMap() const
{
    return _thetasMap;
}

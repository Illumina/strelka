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

#pragma once

#include <vector>
#include <map>
#include "rapidjson/document.h"

/// \brief Data structure to store theta.json
///
class ThetaJson
{
public:
    ThetaJson() {}
    ThetaJson(
            size_t repeatPatternSize,
            const std::vector<double>& theta);
private:
    size_t _repeatPatternSize;
    std::vector<double> _theta;

public:
    static ThetaJson deserialize(const rapidjson::Value& root);

    const std::vector<double>& getTheta() const
    {
        return _theta;
    }

    size_t getRepeatPatternSize() const
    {
        return _repeatPatternSize;
    }

};

class ThetasJson
{
public:
    ThetasJson() {}
    explicit
    ThetasJson(
            std::map<unsigned, std::vector<double>>& thetasMap);
private:
    std::map<unsigned, std::vector<double>> _thetasMap;

public:
    static ThetasJson deserialize(const rapidjson::Value& root);
    const std::map<unsigned, std::vector<double>>& getThetasMap() const;
};

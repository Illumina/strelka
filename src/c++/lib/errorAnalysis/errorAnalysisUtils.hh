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

#include <map>


/// Iterate \p inputMap value for \p key by \p valueIncrement, assuming value is 0 if key is not present
///
/// \tparam V map value type, assumed to be integral
template <typename K, typename V>
void
iterateMapValue(
    std::map<K, V>& inputMap,
    const K& key,
    const unsigned valueIncrement = 1)
{
    const auto mapIterator(inputMap.find(key));
    if (mapIterator == inputMap.end())
    {
        inputMap[key] = valueIncrement;
    }
    else
    {
        mapIterator->second += valueIncrement;
    }
}

/// Merge map1 into map2 to produce map2', such that:
/// (1) map2' to will contain the union of map1 and map2 keys
/// (2) for each shared key, map2'[key] = map1[key]*value + map2[key]
template <typename K, typename V1, typename V2>
void
mergeMapKeys(
    const std::map<K,V1>& map1,
    std::map<K,V2>& map2,
    const unsigned scale = 1)
{
    for (const auto& map1Element : map1)
    {
        const auto iter(map2.find(map1Element.first));
        if (iter == map2.end())
        {
            map2.insert(map1Element);
        }
        else
        {
            iter->second += (map1Element.second * scale);
        }
    }
}

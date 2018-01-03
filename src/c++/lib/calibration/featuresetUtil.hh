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

#include "calibration/VariantScoringModelBase.hh"
#include "calibration/VariantScoringModelMetadata.hh"
#include "blt_util/PolymorphicObject.hh"

#include "boost/dynamic_bitset.hpp"

#include <cassert>

#include <sstream>
#include <string>


struct FeatureSet : public PolymorphicObject
{
    virtual
    unsigned
    size() const = 0;

    virtual
    const char*
    getName() const = 0;

    virtual
    const char*
    getFeatureLabel(const unsigned idx) const = 0;

    /// feature map is typically used to verify a feature set match to a particular model
    VariantScoringModelMetadata::featureMap_t
    getFeatureMap() const
    {
        VariantScoringModelMetadata::featureMap_t result;

        const unsigned featureSize(size());
        for (unsigned featureIndex(0); featureIndex<featureSize; ++featureIndex)
        {
            result.insert(std::make_pair(std::string(getFeatureLabel(featureIndex)),featureIndex));
        }

        return result;
    }

protected:
    FeatureSet() = default;

    explicit
    FeatureSet(const FeatureSet&) = default;

private:
    FeatureSet& operator=(const FeatureSet&) = delete;
};


/// output the primary and development feature set labels, and check for dup labels
void
writeExtendedFeatureSet(
    const FeatureSet& primaryFeatureSet,
    const FeatureSet& developmentFeatureSet,
    const char* featureTypeLabel,
    std::ostream& os);


/// simple feature organizer
///
/// (1) doesn't mix up features with other tracking info
/// (2) generates no system calls after initialization
///
struct VariantScoringFeatureKeeper
{
    explicit
    VariantScoringFeatureKeeper(
        const FeatureSet& featureSet)
        : _featureSet(featureSet),
          _isFeatureSet(featureSet.size()),
          _featureVal(featureSet.size())
    {
        clear();
    }

    void
    set(
        const unsigned featureIndex,
        const double featureValue)
    {
        assert(featureIndex < _featureSet.size());
        if (test(featureIndex))
        {
            featureError(featureIndex, "attempted to set scoring feature twice");
        }
        _featureVal[featureIndex] = featureValue;
        _isFeatureSet.set(featureIndex);
    }

    double
    get(const unsigned featureIndex) const
    {
        assert(featureIndex < _featureSet.size());
        if (! test(featureIndex))
        {
            featureError(featureIndex, "attempted to retrieve scoring feature before it was set");
        }
        return _featureVal[featureIndex];
    }

    const VariantScoringModelBase::featureInput_t&
    getAll() const
    {
        return _featureVal;
    }

    bool
    test(const unsigned featureIndex) const
    {
        assert(featureIndex < _isFeatureSet.size());
        return _isFeatureSet.test(featureIndex);
    }

    void
    writeValues(
        std::ostream& os) const
    {
        const unsigned featureSize(_featureSet.size());
        for (unsigned featureIndex(0); featureIndex<featureSize; ++featureIndex)
        {
            if (featureIndex > 0) os << ',';
            os << get(featureIndex);
        }
    }

    void
    dump(
        std::ostream& os) const
    {
        const unsigned featureSize(_featureSet.size());
        for (unsigned featureIndex(0); featureIndex<featureSize; ++featureIndex)
        {
            if (featureIndex > 0) os << ',';
            os << _featureSet.getFeatureLabel(featureIndex) << ":" << get(featureIndex);
        }
    }

    void
    clear()
    {
        _isFeatureSet.reset();
    }

    bool
    empty() const
    {
        return _isFeatureSet.none();
    }

    const FeatureSet&
    getFeatureSet() const
    {
        return _featureSet;
    }

private:

    void
    featureError(
        const unsigned featureIndex,
        const char* msg) const;

    const FeatureSet& _featureSet;
    boost::dynamic_bitset<> _isFeatureSet;
    VariantScoringModelBase::featureInput_t _featureVal;
};


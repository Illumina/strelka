// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
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

#include <cassert>

#include <bitset>


template <typename FEATURESET>
struct FeatureMapMaker
{
    FeatureMapMaker()
    {
        for (unsigned i(0); i<FEATURESET::SIZE; ++i)
        {
            result.insert(std::make_pair(std::string(FEATURESET::get_feature_label(i)),i));
        }
    }

    VariantScoringModelMetadata::featureMap_t result;
};


/// simple feature organizer
///
/// (1) doesn't mix up features with other tracking info
/// (2) generates no system calls after initialization
///
template <typename FEATURESET>
struct VariantScoringFeatureKeeper
{
    VariantScoringFeatureKeeper()
    {
        clear();
        _featureVal.resize(FEATURESET::SIZE);
    }

    void
    set(const typename FEATURESET::index_t i,double val)
    {
        if (test(i))
        {
            assert(false && "Set scoring feature twice");
        }
        _featureVal[i] = val;
        _isFeatureSet.set(i);
    }

    double
    get(const typename FEATURESET::index_t i) const
    {
        if (! test(i))
        {
            assert(false && "Requesting undefined feature");
        }
        return _featureVal[i];
    }

    const VariantScoringModelBase::featureInput_t&
    getAll() const
    {
        return _featureVal;
    }

    bool
    test(const typename FEATURESET::index_t i) const
    {
        return _isFeatureSet.test(i);
    }

    void
    write(
        std::ostream& os) const
    {
        const unsigned featureSize(FEATURESET::SIZE);
        for (unsigned featureIndex(0); featureIndex<featureSize; ++featureIndex)
        {
            if (featureIndex > 0) os << ',';
            os << FEATURESET::get_feature_label(featureIndex) << ":" << get(featureIndex);
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

private:
    std::bitset<FEATURESET::SIZE> _isFeatureSet;
    VariantScoringModelBase::featureInput_t _featureVal;
};


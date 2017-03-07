// -*- mode: c++; indent-tabs-mode: nil; -*-
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

#include "starling_common/IndelKey.hh"

#include <cassert>

#include <algorithm>
#include <vector>


namespace IndelErrorRateType
{
/// the indel error model current reduces all alternate alleles to the following
/// simplified states
///
enum index_t
{
    INSERT,
    DELETE,
    OTHER
};

inline
index_t
getRateType(
    const IndelKey& indelKey)
{
    if     (indelKey.isPrimitiveDeletionAllele())
    {
        return DELETE;
    }
    else if (indelKey.isPrimitiveInsertionAllele())
    {
        return INSERT;
    }
    else
    {
        return OTHER;
    }
}
}


/// helper object used by IndelErrorModel
///
/// For ref="ATACACACAT" alt="ATACACAT",
/// repeatingPatternSize is 2 and
/// patternRepeatCount is 3
///
struct IndelErrorRateSet
{
    double
    getRate(
        unsigned repeatingPatternSize,
        unsigned patternRepeatCount,
        const IndelErrorRateType::index_t simpleIndelType) const
    {
        assert(_isFinalized);
        assert(repeatingPatternSize>0);
        assert(patternRepeatCount>0);

        // revert back to baseline if our motif size isn't covered
        if (repeatingPatternSize > _errorRates.size())
        {
            repeatingPatternSize=1;
            patternRepeatCount=1;
        }

        const unsigned maxPatternSizeIndex(_errorRates.size()-1);
        repeatingPatternSize = std::min((repeatingPatternSize-1), maxPatternSizeIndex);

        const auto& patternRates(_errorRates[repeatingPatternSize]);
        const unsigned maxRepeatCountIndex(patternRates.size()-1);

        patternRepeatCount = std::min((patternRepeatCount-1), maxRepeatCountIndex);

        const auto& theRates(patternRates[patternRepeatCount]);

        return theRates.getRate(simpleIndelType);
    }

    void
    addRate(
        const unsigned repeatingPatternSize,
        const unsigned patternRepeatCount,
        const double insertionErrorRate,
        const double deletionErrorRate)
    {
        assert(repeatingPatternSize>0);
        assert(patternRepeatCount>0);
        assert(repeatingPatternSize <= maxRepeatingPatternSize);
        assert(patternRepeatCount <= maxPatternRepeatCount);
        if (repeatingPatternSize > _errorRates.size())
        {
            _errorRates.resize(repeatingPatternSize);
        }

        auto& patternRates(_errorRates[repeatingPatternSize-1]);

        if (patternRepeatCount > patternRates.size())
        {
            patternRates.resize(patternRepeatCount);
        }

        IndelErrorRates& theRates(patternRates[patternRepeatCount-1]);

        assert(! theRates.isInit);
        theRates.isInit=true;
        theRates.insertionErrorRate=insertionErrorRate;
        theRates.deletionErrorRate=deletionErrorRate;
    }

    //// this must be called before calling getRates
    void
    finalizeRates()
    {
        // ensure that the (possibly jagged) matrix is complete:
        assert(_errorRates.size() > 0);
        for (const auto& patternRates : _errorRates)
        {
            assert(patternRates.size() > 0);
            for (const auto& theRates : patternRates)
            {
                assert(theRates.isInit);
            }
        }

        _isFinalized = true;
    }

private:
    // these limits are for QC purposes only
    static const unsigned maxRepeatingPatternSize = 4;
    static const unsigned maxPatternRepeatCount = 40;

    struct IndelErrorRates
    {
        bool
        isRate() const
        {
            return isInit;
        }

        double
        getRate(const IndelErrorRateType::index_t simpleIndelType) const
        {
            using namespace IndelErrorRateType;
            switch (simpleIndelType)
            {
            case DELETE:
                return deletionErrorRate;
            case INSERT:
                return insertionErrorRate;
            default:
                assert(false && "Unexpected indel type");
                return 0.;
            }
        }

        bool isInit = false;
        double insertionErrorRate = 0;
        double deletionErrorRate = 0;
    };

    bool _isFinalized = false;
    std::vector<std::vector<IndelErrorRates>> _errorRates;
};

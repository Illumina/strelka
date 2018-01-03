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

#include "starling_common/IndelKey.hh"

#include <cassert>

#include <algorithm>
#include <vector>


namespace IndelErrorRateType
{
/// the indel error model reduces all alternate alleles to the following
/// simplified states
///
enum index_t
{
    INSERT,
    DELETE,
    NOISYLOCUS,
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


/// \brief Helper object used by IndelErrorModel to store a set of error rates
///
/// The rate storage grows dynamically as entries are given with different repeatingPatternSize and patternRepeatCount,
/// with the following constraints:
/// - At least one repeatingPatternSize must be defined
/// - If a repeatingPatternSize is defined, all values from 1..repeatingPatternSize must also be defined.
/// - For any repeatingPatternSize which is defined, at least one patternRepeatCount must be defined.
/// - For any repeatingPatternSize, if a patternRepeatCount value is defined, all values from 1..patternRepeatCount must
///   also be defined before finalizing the rate set.
///
/// Terminology:
///
/// For REF="TACACAC" ALT="TACAC",
/// repeatingPatternSize is 2 (referring to the "AC" repeat unit)
/// patternRepeatCount is 3 (taken from the reference)
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

        // revert back to baseline if the repeating pattern size isn't represented
        if (repeatingPatternSize > _errorRates.size())
        {
            repeatingPatternSize=1;
            patternRepeatCount=1;
        }

        const unsigned maxRepeatingPatternSizeIndex(_errorRates.size()-1);
        const unsigned repeatingPatternSizeIndex = std::min((repeatingPatternSize-1), maxRepeatingPatternSizeIndex);

        const auto& repeatingPatternSizeRates(_errorRates[repeatingPatternSizeIndex]);

        const unsigned maxPatternRepeatCountIndex(repeatingPatternSizeRates.size()-1);
        const unsigned patternRepeatCountIndex = std::min((patternRepeatCount-1), maxPatternRepeatCountIndex);

        const auto& indelRates(repeatingPatternSizeRates[patternRepeatCountIndex]);

        return indelRates.getRate(simpleIndelType);
    }

    void
    addRate(
        const unsigned repeatingPatternSize,
        const unsigned patternRepeatCount,
        const double insertionErrorRate,
        const double deletionErrorRate,
        const double noisyLocusRate = 0)
    {
        assert(repeatingPatternSize>0);
        assert(patternRepeatCount>0);
        assert(repeatingPatternSize <= maxRepeatingPatternSize);
        assert(patternRepeatCount <= maxPatternRepeatCount);
        if (repeatingPatternSize > _errorRates.size())
        {
            _errorRates.resize(repeatingPatternSize);
        }

        auto& repeatingPatternSizeRates(_errorRates[repeatingPatternSize-1]);

        if (patternRepeatCount > repeatingPatternSizeRates.size())
        {
            repeatingPatternSizeRates.resize(patternRepeatCount);
        }

        IndelErrorRates& indelRates(repeatingPatternSizeRates[patternRepeatCount-1]);

        assert(! indelRates.isInit);
        indelRates.isInit=true;
        indelRates.insertionErrorRate=insertionErrorRate;
        indelRates.deletionErrorRate=deletionErrorRate;
        indelRates.noisyLocusRate=noisyLocusRate;
    }

    /// \brief Check for a valid rate initialization pattern
    ///
    /// This must be called before calling getRates
    void
    finalizeRates()
    {
        // ensure that the (possibly jagged) matrix is complete:
        assert(_errorRates.size() > 0);
        for (const auto& repeatingPatternSizeRates : _errorRates)
        {
            assert(repeatingPatternSizeRates.size() > 0);
            for (const auto& indelRates : repeatingPatternSizeRates)
            {
                assert(indelRates.isInit);
            }
        }

        _isFinalized = true;
    }

private:
    // these limits are not used to define any data structure size bounds, rather
    // they check against implausible input values for QC purposes
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
            case NOISYLOCUS:
                return noisyLocusRate;
            default:
                assert(false && "Unexpected indel type");
                return 0.;
            }
        }

        bool isInit = false;
        double insertionErrorRate = 0;
        double deletionErrorRate = 0;
        double noisyLocusRate = 0;
    };

    bool _isFinalized = false;
    std::vector<std::vector<IndelErrorRates>> _errorRates;
};

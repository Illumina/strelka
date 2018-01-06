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

#include "featuresetUtil.hh"
#include "common/Exceptions.hh"

#include <iostream>
#include <set>
#include <sstream>



static
void
repeatedFeatureLabelError(
    const char* label,
    const std::string& featureLabel)
{
    using namespace illumina::common;

    std::ostringstream oss;
    oss << "Repeated " << label << " EVS training feature label '" << featureLabel << "'";
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
}



void
writeExtendedFeatureSet(
    const FeatureSet& primaryFeatureSet,
    const FeatureSet& developmentFeatureSet,
    const char* featureTypeLabel,
    std::ostream& os)
{
    std::set<std::string> featureLabels;
    bool isFirst(true);
    auto printFeatureSet = [&](const FeatureSet& featureSet)
    {
        const unsigned featureSize(featureSet.size());
        for (unsigned featureIndex = 0; featureIndex < featureSize; ++featureIndex)
        {
            if (isFirst)
            {
                isFirst=false;
            }
            else
            {
                os << ",";
            }
            const std::string featureLabel(featureSet.getFeatureLabel(featureIndex));
            os << featureLabel;

            const auto retVal = featureLabels.insert(featureLabel);
            if (not retVal.second)
            {
                repeatedFeatureLabelError(featureTypeLabel, featureLabel);
            }
        }
    };

    printFeatureSet(primaryFeatureSet);
    printFeatureSet(developmentFeatureSet);
}



void
VariantScoringFeatureKeeper::
featureError(
    const unsigned featureIndex,
    const char* msg) const
{
    using namespace illumina::common;

    assert(msg);

    std::ostringstream oss;
    oss << "ERROR: " << msg << "."
        << " Feature: '" << _featureSet.getFeatureLabel(featureIndex) << "'"
        << " from set: '" << _featureSet.getName() << "'\n";
    BOOST_THROW_EXCEPTION(GeneralException(oss.str()));
}

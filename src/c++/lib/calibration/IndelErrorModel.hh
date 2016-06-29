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
/*
 *      Author: mkallberg
 */

#pragma once

#include "IndelErrorModelMetadata.hh"
#include "IndelErrorRateSet.hh"

#include "starling_common/starling_indel_report_info.hh"


/// organizes indel error rate information
///
struct IndelErrorModel
{
    /// Initialize indel error model to one of the hard-coded variants compiled into
    /// Strelka (if modelFilename is empty), or from a json parameter file (otherwise)
    ///
    IndelErrorModel(
        const std::string& modelName,
        const std::string& modelFilename);

    /// Retrieve indel error rates for a specific indel type
    void
    getIndelErrorRate(
        const IndelKey& indelKey,
        const starling_indel_report_info& indelReportInfo,
        double& refToIndelErrorProb,
        double& indelToRefErrorProb) const;

private:
    IndelErrorModelMetadata _meta;
    IndelErrorRateSet _errorRates;

#if 0
    const std::string&
    getName() const
    {
        return _meta.name;
    }


    unsigned get_max_motif_length() const
    {
        return MaxMotifLength;
    }
#endif
};

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

/// \author Chris Saunders
///

#pragma once

#include "denovo_call.hh"
#include "denovo_indel_call.hh"
#include "pedicure_shared.hh"
#include "starling_common/starling_indel_report_info.hh"

#include <array>
#include <iosfwd>


typedef std::array<starling_indel_sample_report_info,PEDICURE_TIERS::SIZE> isriTiers_t;


void
denovo_indel_call_vcf(
    const pedicure_options& opt,
    const pedicure_deriv_options& dopt,
    const SampleInfoManager& sinfo,
    const denovo_indel_call& dinc,
    const starling_indel_report_info& indelReportInfo,
    const std::vector<isriTiers_t>& isri,
    std::ostream& os);

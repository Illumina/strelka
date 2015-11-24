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
#include "starling_common/starling_base_shared.hh"
#include "blt_common/snp_pos_info.hh"

#include "blt_util/qscore.hh"
#include "gvcf_locus_info.hh"
#include <boost/utility.hpp>


class starling_continuous_variant_caller : private boost::noncopyable
{
public:

    static void position_snp_call_continuous(
        const starling_base_options& opt,
        const snp_pos_info& good_pi,
        continuous_site_info& info);
    static void add_indel_call(
        const starling_base_options& opt,
        const indel_key& ik,
        const indel_data& id,
        const starling_indel_report_info& iri,
        const starling_indel_sample_report_info& isri,
        continuous_indel_info& info);
    static unsigned poisson_qscore(unsigned callCount, unsigned coverage, unsigned estimatedBaseCallQuality, unsigned maxQScore);

    static double strand_bias(unsigned fwdAlt, unsigned revAlt, unsigned fwdOther, unsigned revOther, double noise);
};

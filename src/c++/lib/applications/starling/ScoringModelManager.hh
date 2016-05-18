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
 *  Created on: Oct 10, 2013
 *  Author: Morten Kallberg
 */

#pragma once

#include "starling_shared.hh"

#include "boost/utility.hpp"
#include "LogisticAndRuleScoringModels.hh"


/// handles site and indel filter labeling and also possibly EVS scoring
/// for (1) a logisitic model read from a file (2) a default hard-threshold
/// model and (3) and file-based hard-threshold model.
///
class ScoringModelManager : private boost::noncopyable
{
public:
    ScoringModelManager(
        const starling_options& init_opt,
        const gvcf_deriv_options& init_dopt)
        : opt(init_opt.gvcf),
          dopt(init_dopt)
    {
        load_models(
                init_opt.germline_variant_scoring_models_filename,
                init_opt.germline_variant_scoring_model_name);
    }

    void
    classify_site(
        const digt_site_info& si,
        digt_call_info& smod) const;

    void
    classify_indel(
        const digt_indel_info& ii,
        digt_indel_call& call) const;

    void
    classify_indels(
        std::vector<std::unique_ptr<digt_indel_info>>& indels) const;

    // mimics behavior of previous hard filters
    void default_classify_site(const site_info& si,
                               shared_call_info& call) const;

    void default_classify_indel(shared_indel_call_info& call) const;


    bool is_current_logistic() const;

    int
    get_case_cutoff(const CALIBRATION_MODEL::var_case my_case) const;

private:
    LogisticAndRuleScoringModels& get_model() { return *modelPtr; }
    const LogisticAndRuleScoringModels& get_model() const { return *modelPtr; }

    bool check_is_model_usable(const digt_indel_info& ii) const;
    void set_indel_modifiers(const digt_indel_info& ii, digt_indel_call& call) const;

    void
    classify_indel_impl(
        const bool is_model_usable,
        const digt_indel_info& ii,
        digt_indel_call& call) const;

    // set options
    void load_models(
        const std::string& model_file,
        const std::string& name);

    bool
    is_default_model() const
    {
        return (!modelPtr);
    }

    // for setting the vcf header filters
    const gvcf_options& opt;
    const gvcf_deriv_options& dopt;

    std::unique_ptr<LogisticAndRuleScoringModels> modelPtr;
};


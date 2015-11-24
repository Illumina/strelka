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

#include "cmodel.hh"
#include "starling_shared.hh"

#include "boost/utility.hpp"



class calibration_models : private boost::noncopyable
{
public:
    calibration_models(
        const starling_options& init_opt,
        const gvcf_deriv_options& init_dopt)
        : opt(init_opt.gvcf),
          dopt(init_dopt)
    {
        load_models(init_opt.germline_variant_scoring_models_filename);
        set_model(init_opt.germline_variant_scoring_model_name);
    }

    void
    clasify_site(
        const digt_site_info& si,
        digt_call_info& smod) const;

    void
    clasify_indel(
        const digt_indel_info& ii,
        digt_indel_call& call) const;

    void
    clasify_indels(
        std::vector<std::unique_ptr<digt_indel_info>>& indels) const;

    // mimics behavior of previous hard filters
    void  default_clasify_site(const site_info& si,
                               shared_call_info& call) const;

    void default_clasify_indel(shared_indel_call_info& call) const;


    bool is_current_logistic() const;

    int
    get_case_cutoff(const CALIBRATION_MODEL::var_case my_case) const;

    const char* get_model_name() const
    {
        return model_name.c_str();
    }

private:
    c_model& get_model(const std::string& name);
    const c_model& get_model(const std::string& name) const;

    bool check_is_model_usable(const digt_indel_info& ii) const;
    void set_indel_modifiers(const digt_indel_info& ii, digt_indel_call& call) const;

    void
    clasify_indel_impl(
        const bool is_model_usable,
        const digt_indel_info& ii,
        digt_indel_call& call) const;

    // set options
    void set_model(const std::string& name);  // set the calibration model to use
    void load_models(std::string model_file); // read in model parameters

    void load_chr_depth_stats();
    void add_model_pars(std::string& name,parmap& my_pars);


    // for setting the vcf header filters
    const gvcf_options& opt;
    const gvcf_deriv_options& dopt;
    std::string model_name="DEFAULT";
    bool is_default_model=true;

    cdmap_t chrom_depth;
    bool has_depth=false;
    int chr_avg=30, chr_median=30;
    typedef std::map<std::string,c_model> modelmap;
    typedef std::map<std::string, double> featuremap;
    modelmap models;
};


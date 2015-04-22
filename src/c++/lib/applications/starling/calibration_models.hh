// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
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
        load_models(init_opt.calibration_models_filename);
        set_model(init_opt.calibration_model);
    }

    void
    clasify_site(
        const site_info& si,
        site_modifiers& smod) const;

    void
    clasify_indel(
        const indel_info& ii,
        indel_modifiers& imod) const;

    void
    clasify_indels(
        std::vector<indel_info>& indels) const;


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

    bool can_use_model(const indel_info& ii) const;
    void set_indel_modifiers(const indel_info& ii, indel_modifiers& imod) const;

    // set options
    void set_model(const std::string& name);  // set the calibration model to use
    void load_models(std::string model_file); // read in model parameters

    void load_chr_depth_stats();
    void add_model_pars(std::string& name,parmap& my_pars);

    // mimics behavior of previous hard filters
    void
    default_clasify_site(
        const site_info& si,
        site_modifiers& smod) const;

    void
    default_clasify_indel(
        const indel_info& ii,
        indel_modifiers& imod) const;

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


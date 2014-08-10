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
 * calbrationmodels.h
 *
 *  Created on: Oct 10, 2013
 *  Author: Morten Kallberg
 */

#pragma once

#include "blt_util/blt_exception.hh"
#include "starling_common/gvcf_block_site_record.hh"
#include "starling_common/gvcf_locus_info.hh"
#include "starling_common/cmodel.hh"


//forward declaration
struct gvcf_options;
struct gvcf_deriv_options;

class calibration_models
{
public:
    calibration_models();
    virtual ~calibration_models();

    // set options
    void set_model(const std::string& name);  // set the calibration model to use
    void load_models(std::string model_file); // read in model parameters

    void clasify_site(site_info& si);
    void clasify_site(indel_info& ii);

    c_model& get_model(std::string& name);

    // mimics behavior of previous hard filters
    void default_clasify_site(site_info& si);
    void default_clasify_site(indel_info& ii);
    bool is_current_logistic();

    // for setting the vcf header filters
    int get_case_cutoff(CALIBRATION_MODEL::var_case my_case);

    void add_model_pars(std::string& name,parmap& my_pars);
    const gvcf_deriv_options* dopt;
    const gvcf_options* opt;
    std::string model_name;
    bool is_default_model=true;
private:
    typedef std::map<std::string,c_model> modelmap;
    typedef std::map<std::string, double> featuremap;
    modelmap models;
};


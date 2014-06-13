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

    void set_model(const std::string& name);  // set the calibration model to use
    void load_models(std::string model_file);                       // read in model parameters

    void clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, site_info& si);
    void clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, indel_info& ii);

    c_model& get_model(std::string& name);

    // mimics behavior of previous hard filters
    void default_clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, site_info& si);
    void default_clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, indel_info& ii);

    void add_model_pars(std::string& name,parmap& my_pars);
private:
    typedef std::map<std::string,c_model> modelmap;
    typedef std::map<std::string, double> featuremap;
    std::string model_name;
    modelmap models;
};


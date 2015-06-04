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
 * Author: Morten Kallberg
 */

#include "calibration_models.hh"

#include <cassert>
#include <exception>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <stdlib.h>     /* atof */


//#define DEBUG_CAL

#ifdef DEBUG_CAL
#include "blt_util/log.hh"
#endif



int
calibration_models::
get_case_cutoff(
    const CALIBRATION_MODEL::var_case my_case) const
{
    return get_model(model_name).get_var_threshold(my_case);
}



bool calibration_models::is_current_logistic() const
{
    if (is_default_model) return false;
    return get_model(model_name).is_logistic_model();
}


void
calibration_models::
clasify_site(
    const site_info& si,
    site_modifiers& smod) const
{
    //si.smod.filters.reset(); // make sure no filters have been applied prior
    if ((si.dgt.is_snp) && (!is_default_model))
    {
        get_model(model_name).score_site_instance(si, smod);
    }
    else
    {
        // don't know what to do with this site, throw it to the old default filters
        default_clasify_site(si, smod);
    }
}

bool 
calibration_models::can_use_model(const indel_info& ii) const
{
    return ((ii.iri().it == INDEL::INSERT || ii.iri().it == INDEL::DELETE) &&
		(!is_default_model) && 
		(ii.dindel.max_gt != STAR_DIINDEL::NOINDEL) ); // VQSR does not handle homref sites properly
}

void
calibration_models::
set_indel_modifiers(const indel_info& ii, indel_modifiers& imod) const
{
    if ((ii.dindel.max_gt != ii.dindel.max_gt_poly) || ii.dindel.is_zero_coverage)
    {
        imod.gqx=0;
    }
    else
    {
        imod.gqx=std::min(ii.dindel.max_gt_poly_qphred,ii.dindel.max_gt_qphred);
    }
    imod.max_gt=ii.dindel.max_gt_poly;
    imod.gq=ii.dindel.max_gt_poly_qphred;
}


void
calibration_models::
clasify_indel(
    const indel_info& ii,
    indel_modifiers& imod) const
{
    set_indel_modifiers(ii, imod);
    
    if ( can_use_model(ii) )
    {
        get_model(model_name).score_indel_instance(ii, imod);
    }
    else
    {
        default_clasify_indel(ii, imod);
    }
}

void
calibration_models::
clasify_indels(
        std::vector<indel_info>& indels) const
{
    bool use_model = true;
    for (auto it = indels.begin(); it != indels.end() && use_model; ++it)
    {
        use_model = can_use_model(*it);
    }
    for (auto it = indels.begin(); it != indels.end(); ++it)
    {
        indel_info& ii(*it);
        
        set_indel_modifiers(ii, ii.imod());
        if ( use_model )
        {
            get_model(model_name).score_indel_instance(ii, ii.imod());
        }
        else
        {
            default_clasify_indel(ii, ii.imod());
        }
    }
}




void
calibration_models::
default_clasify_site(
    const site_info& si,
    site_modifiers& smod) const
{
    if (opt.is_min_gqx)
    {
        if (si.smod.gqx<opt.min_gqx) smod.set_filter(VCF_FILTERS::LowGQX);
    }
    if (dopt.is_max_depth())
    {
        if ((si.n_used_calls+si.n_unused_calls) > dopt.max_depth) smod.set_filter(VCF_FILTERS::HighDepth);
    }
    // high DPFratio filter
    if (opt.is_max_base_filt)
    {
        const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
        if (total_calls>0)
        {
            const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
            if (filt>opt.max_base_filt) smod.set_filter(VCF_FILTERS::HighBaseFilt);
        }
    }
    if (si.dgt.is_snp)
    {
        if (opt.is_max_snv_sb)
        {
            if (si.dgt.sb>opt.max_snv_sb) smod.set_filter(VCF_FILTERS::HighSNVSB);
        }
        if (opt.is_max_snv_hpol)
        {
            if (static_cast<int>(si.hpol)>opt.max_snv_hpol) smod.set_filter(VCF_FILTERS::HighSNVHPOL);
        }
    }
}



// default rules based indel model
void
calibration_models::
default_clasify_indel(
    const indel_info& ii,
    indel_modifiers& imod) const
{
    if (this->opt.is_min_gqx)
    {
        if (ii.imod().gqx<opt.min_gqx) imod.set_filter(VCF_FILTERS::LowGQX);
    }

    if (this->dopt.is_max_depth())
    {
        if (ii.isri().depth > this->dopt.max_depth) imod.set_filter(VCF_FILTERS::HighDepth);
    }

    if (this->opt.is_max_ref_rep)
    {
        if (ii.iri().is_repeat_unit())
        {
            if ((ii.iri().repeat_unit.size() <= 2) &&
                (static_cast<int>(ii.iri().ref_repeat_count) > this->opt.max_ref_rep))
            {
                imod.set_filter(VCF_FILTERS::HighRefRep);
            }
        }
    }
}



void
calibration_models::load_chr_depth_stats()
{
    if (! this->opt.chrom_depth_file.empty())
    {
        parse_chrom_depth(this->opt.chrom_depth_file,chrom_depth);
        std::vector<double> depths;
        for (cdmap_t::const_iterator iter = this->chrom_depth.begin(); iter != this->chrom_depth.end() ; ++iter)
        {
#if 0
            /// NO HARD CODED CHROMOSOME NAMES!!!!!!!!!!!!
            if (!(iter->first=="chrM" ||iter->first=="chrY"||iter->first=="chrX"))
#endif
            {
//                log_os << iter->first << " is " << iter->second << "\n";
                depths.push_back(iter->second);
            }
        }
        std::sort(depths.begin(), depths.end());
        int sum = 0;
        for (unsigned i=0; i<depths.size(); i++)
        {
//            log_os << depths.at(i) << "\n";
            sum += depths.at(i);
        }
        this->chr_median    = depths.at(depths.size()/2);
        this->chr_avg       = sum/depths.size();
//        log_os << "median " << this->chr_median << "\n";
//        log_os << "mean " << this->chr_avg << "\n";
        this->has_depth = true;
    }
}

void calibration_models::set_model(const std::string& name)
{
    if (name.empty()) return;

    modelmap::iterator it = this->models.find(boost::to_upper_copy(name));
    assert("Unrecognized calibration model given using --scoring-model option in set_model" && it != this->models.end());
    this->load_chr_depth_stats();

    if (it != this->models.end())
    {
//        if (this->has_depth && (this->chr_median>70 || this->chr_median<10))
//        {
//            this->model_name = "QRULE";     //TODO hacky fix for defaulting to qrule if we have a high median chromosome depth (VQSR not trained to handle these cases yet)
//        }
//        else
        this->model_name = boost::to_upper_copy(name);
        this->is_default_model = false;
    }
//    log_os << "Calibration model set to '" << this->model_name << "'\n";
//    log_os << "Calibration model is default '" << this->is_default_model << "'\n";
#ifdef DEBUG_CAL
    log_os << "Calibration model set to '" << this->model_name << "'\n";
#endif
}

c_model&
calibration_models::get_model(const std::string& name)
{
    auto it = this->models.find(boost::to_upper_copy(name));
    assert("Unrecognized calibration model given using --scoring-model option" && it != this->models.end());
    return it->second;
}

const c_model&
calibration_models::get_model(const std::string& name) const
{
    auto it = this->models.find(boost::to_upper_copy(name));
    assert("Unrecognized calibration model given using --scoring-model option" && it != this->models.end());
    return it->second;
}

void calibration_models::add_model_pars(std::string& name,parmap& my_pars)
{
#ifdef DEBUG_CAL
    log_os << "Adding pars for model " << name << "\n";
    log_os << "Adding pars " << my_pars.size() << "\n";
#endif
    this->get_model(name).add_parameters(my_pars);
    my_pars.clear();
}


void calibration_models::load_models(std::string model_file)
{
    using namespace boost::algorithm;

#ifdef DEBUG_CAL
    log_os << "Loading models from file: " << model_file << "\n";
#endif
    std::ifstream myReadFile;
    myReadFile.open(model_file.c_str());
    std::string output;
    std::string parspace;
    std::string submodel;
    std::string current_name;
    parmap pars;
    if (myReadFile.is_open())
    {
        while (!myReadFile.eof())
        {
            std::getline (myReadFile,output);
            std::vector<std::string> tokens;
            split(tokens, output, is_any_of(" \t")); // tokenize string
            //case new model
            if (tokens.at(0).substr(0,3)=="###")
            {
                if (! pars.empty())
                {
                    this->add_model_pars(current_name,pars);
                }
                current_name = boost::to_upper_copy(tokens.at(1));
                c_model current_model(current_name,tokens.at(2),dopt);

                this->models.insert(modelmap::value_type(current_name, current_model));
#ifdef DEBUG_CAL
                log_os << "Loading model: " << tokens.at(1) << " Type: " << tokens.at(2) << "\n";
#endif
            }
            //load submodel
            else if (tokens.at(0)=="#")
            {
#ifdef DEBUG_CAL
                log_os << "submodel: " << tokens.at(1) << " parspace: " << tokens.at(2) << "\n";
#endif
                submodel = tokens.at(1);
                parspace = tokens.at(2);
            }
            //case load parameters
            else
            {
                if (tokens.size()>1)
                {
#ifdef DEBUG_CAL
                    log_os << " setting " << tokens.at(0) << " = " << tokens.at(1) << "\n";
#endif
                    pars[submodel][parspace][tokens.at(0)] = atof(tokens.at(1).c_str());
                }
            }
        }
        this->add_model_pars(current_name,pars);
    }
#ifdef DEBUG_CAL
    log_os << "Done loading models" << "\n";
#endif
    myReadFile.close();
}

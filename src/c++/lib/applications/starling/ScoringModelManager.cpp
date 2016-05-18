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
 * Author: Morten Kallberg
 */

#include "ScoringModelManager.hh"

#include "common/Exceptions.hh"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/classification.hpp"

#include <cassert>
#include <cstdlib>     /* atof */

#include <exception>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <iterator>
#include <map>


//#define DEBUG_CAL

#ifdef DEBUG_CAL
#include "blt_util/log.hh"
#endif



int
ScoringModelManager::
get_case_cutoff(
    const CALIBRATION_MODEL::var_case my_case) const
{
    if (is_default_model()) return 0;
    return get_model().get_var_threshold(my_case);
}



bool ScoringModelManager::is_current_logistic() const
{
    if (is_default_model()) return false;
    return get_model().is_logistic_model();
}


void
ScoringModelManager::
classify_site(
    const digt_site_info& si,
    digt_call_info& smod) const
{
    //si.smod.filters.reset(); // make sure no filters have been applied prior
    if ((si.dgt.is_snp) && (!is_default_model()))
    {
        get_model().score_site_instance(si, smod);
    }
    else
    {
        // don't know what to do with this site, throw it to the old default filters
        default_classify_site(si, smod);
    }
}

bool
ScoringModelManager::check_is_model_usable(const digt_indel_info& ii) const
{
    const auto& call(ii.first());
    return ((call._iri.it == INDEL::INSERT || call._iri.it == INDEL::DELETE) &&
            (!is_default_model()) &&
            (call._dindel.max_gt != STAR_DIINDEL::NOINDEL) ); // empirical scoring does not handle homref sites properly
}

void
ScoringModelManager::
set_indel_modifiers(const digt_indel_info& ii, digt_indel_call& call) const
{
    const auto& dindel(ii.first()._dindel);
    if ((dindel.max_gt != dindel.max_gt_poly) || dindel.is_zero_coverage)
    {
        call.gqx=0;
    }
    else
    {
        call.gqx=std::min(dindel.max_gt_poly_qphred,dindel.max_gt_qphred);
    }
    call.max_gt=dindel.max_gt_poly;
    call.gq=dindel.max_gt_poly_qphred;
}



void
ScoringModelManager::
classify_indel_impl(
    const bool is_model_usable,
    const digt_indel_info& ii,
    digt_indel_call& call) const
{
    set_indel_modifiers(ii, call);

    if ( is_model_usable )
    {
        get_model().score_indel_instance(ii, call);
    }
    else
    {
        default_classify_indel(call);
    }
}


void
ScoringModelManager::
classify_indel(
    const digt_indel_info& ii,
    digt_indel_call& call) const
{
    classify_indel_impl(check_is_model_usable(ii),ii,call);
}

void
ScoringModelManager::
classify_indels(
    std::vector<std::unique_ptr<digt_indel_info>>& indels) const
{
    bool is_model_usable = true;
    for (const auto& indel : indels)
    {
        if (! is_model_usable) break;
        is_model_usable = check_is_model_usable(*indel);
    }

    for (auto& indel : indels)
    {
        digt_indel_info& ii(*indel);
        classify_indel_impl(is_model_usable,ii,ii.first());
    }
}




void
ScoringModelManager::
default_classify_site(const site_info& si,
                      shared_call_info& call) const
{
    if (opt.is_min_gqx)
    {
        if (call.gqx<opt.min_gqx) call.set_filter(VCF_FILTERS::LowGQX);
    }
    if (dopt.is_max_depth())
    {
        if ((si.n_used_calls+si.n_unused_calls) > dopt.max_depth) call.set_filter(VCF_FILTERS::HighDepth);
    }
    // high DPFratio filter
    if (opt.is_max_base_filt)
    {
        const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
        if (total_calls>0)
        {
            const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
            if (filt>opt.max_base_filt) call.set_filter(VCF_FILTERS::HighBaseFilt);
        }
    }
    if (si.is_snp())
    {
        if (opt.is_max_snv_sb)
        {
            if (call.strand_bias>opt.max_snv_sb) call.set_filter(VCF_FILTERS::HighSNVSB);
        }
        if (opt.is_max_snv_hpol)
        {
            if (static_cast<int>(si.hpol)>opt.max_snv_hpol) call.set_filter(VCF_FILTERS::HighSNVHPOL);
        }
    }
}



// default rules based indel model
void
ScoringModelManager::
default_classify_indel(shared_indel_call_info& call) const
{
    if (this->opt.is_min_gqx)
    {
        if (call.gqx<opt.min_gqx) call.set_filter(VCF_FILTERS::LowGQX);
    }

    if (this->dopt.is_max_depth())
    {
        if (call._isri.depth > this->dopt.max_depth) call.set_filter(VCF_FILTERS::HighDepth);
    }

    if (this->opt.is_max_ref_rep)
    {
        const auto& iri(call._iri);
        if (iri.is_repeat_unit())
        {
            if ((iri.repeat_unit.size() <= 2) &&
                (static_cast<int>(iri.ref_repeat_count) > this->opt.max_ref_rep))
            {
                call.set_filter(VCF_FILTERS::HighRefRep);
            }
        }
    }
}



void
ScoringModelManager::
load_models(
    const std::string& model_file,
    const std::string& name)
{
    using namespace boost::algorithm;

    if (model_file.empty()) return;
    if (name.empty()) return;

#ifdef DEBUG_CAL
    log_os << "Loading models from file: " << model_file << "\n";
#endif
    std::ifstream myReadFile;
    myReadFile.open(model_file.c_str());

    if (! myReadFile)
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: Failed to load germline variant scoring file '" << model_file << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }

    std::string parspace;
    std::string submodel;
    std::string modelType;
    parmap pars;

    bool isInCurrentModel(false);

    while (!myReadFile.eof())
    {
        std::vector<std::string> tokens;
        {
            std::string output;
            std::getline(myReadFile,output);
            split(tokens, output, is_any_of(" \t")); // tokenize string
        }

        //case new model
        if (tokens.at(0).substr(0,3)=="###")
        {
            if (isInCurrentModel) break;
            isInCurrentModel=(tokens.at(1)==name);

            if (! isInCurrentModel) continue;
            modelType=tokens.at(2);
        }

        if (! isInCurrentModel) continue;
        //load submodel
        if (tokens.at(0)=="#")
        {
            submodel = tokens.at(1);
            parspace = tokens.at(2);
        }
        //case load parameters
        else
        {
            if (tokens.size()>1)
            {
                pars[submodel][parspace][tokens.at(0)] = atof(tokens.at(1).c_str());
            }
        }
    }

    if (! isInCurrentModel)
    {
        using namespace illumina::common;

        std::ostringstream oss;
        oss << "ERROR: unrecognized variant scoring model name: '" << name << "'\n";
        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
    }
    assert(! pars.empty());
    modelPtr.reset(new LogisticAndRuleScoringModels(modelType,dopt,pars));
}

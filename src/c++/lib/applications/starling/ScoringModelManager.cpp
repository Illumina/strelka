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
    const GermlineDiploidSiteCallInfo& si,
    GermlineDiploidSiteSimpleGenotypeInfo& smod) const
{
    if (si.dgt.is_snp && _isReportEVSFeatures)
    {
        // when reporting is turned on, we need to compute EVS features
        // for any usable variant regardless of EVS model type:
        si.computeEmpiricalScoringFeatures(_isReportEVSFeatures, _dopt.norm_depth, smod);
    }

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
ScoringModelManager::
checkIsVariantUsableInEVSModel(const GermlineDiploidIndelCallInfo& ii) const
{
    const auto& call(ii.first());
    return ((call._iri.it == INDEL::INSERT || call._iri.it == INDEL::DELETE) &&
            (call._dindel.max_gt != STAR_DIINDEL::NOINDEL) ); // empirical scoring does not handle homref sites properly
}

void
ScoringModelManager::
set_indel_modifiers(
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
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
    const bool isVariantUsableInEVSModel,
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    set_indel_modifiers(ii, call);

    if (isVariantUsableInEVSModel && _isReportEVSFeatures)
    {
        // when reporting is turned on, we need to compute EVS features
        // for any usable variant regardless of EVS model type:
        call.computeEmpiricalScoringFeatures(_isReportEVSFeatures, _dopt.norm_depth, ii.is_hetalt());
    }

    if ( (! is_default_model()) && isVariantUsableInEVSModel )
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
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    classify_indel_impl(checkIsVariantUsableInEVSModel(ii),ii,call);
}

void
ScoringModelManager::
classify_indels(
    std::vector<std::unique_ptr<GermlineDiploidIndelCallInfo>>& indels) const
{
    bool isVariantUsableInEVSModel = true;
    for (const auto& indel : indels)
    {
        if (! isVariantUsableInEVSModel) break;
        isVariantUsableInEVSModel = checkIsVariantUsableInEVSModel(*indel);
    }

    for (auto& indel : indels)
    {
        GermlineDiploidIndelCallInfo& ii(*indel);
        classify_indel_impl(isVariantUsableInEVSModel,ii,ii.first());
    }
}




void
ScoringModelManager::
default_classify_site(
    const GermlineSiteCallInfo& si,
    GermlineVariantSimpleGenotypeInfo& call) const
{
    if (_opt.is_min_gqx)
    {
        if (call.gqx<_opt.min_gqx) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    }
    if (_dopt.is_max_depth())
    {
        if ((si.n_used_calls+si.n_unused_calls) > _dopt.max_depth) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }
    // high DPFratio filter
    if (_opt.is_max_base_filt)
    {
        const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
        if (total_calls>0)
        {
            const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
            if (filt>_opt.max_base_filt) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighBaseFilt);
        }
    }
    if (si.is_snp())
    {
        if (_opt.is_max_snv_sb)
        {
            if (call.strand_bias>_opt.max_snv_sb) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighSNVSB);
        }
        if (_opt.is_max_snv_hpol)
        {
            if (static_cast<int>(si.hpol)>_opt.max_snv_hpol) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighSNVHPOL);
        }
    }
}



// default rules based indel model
void
ScoringModelManager::
default_classify_indel(
    GermlineIndelSimpleGenotypeInfo& call) const
{
    if (this->_opt.is_min_gqx)
    {
        if (call.gqx<_opt.min_gqx) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    }

    if (this->_dopt.is_max_depth())
    {
        if (call._isri.depth > this->_dopt.max_depth) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }

    if (this->_opt.is_max_ref_rep)
    {
        const auto& iri(call._iri);
        if (iri.is_repeat_unit())
        {
            if ((iri.repeat_unit.size() <= 2) &&
                (static_cast<int>(iri.ref_repeat_count) > this->_opt.max_ref_rep))
            {
                call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighRefRep);
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
    modelPtr.reset(new LogisticAndRuleScoringModels(modelType,_dopt,pars));
}

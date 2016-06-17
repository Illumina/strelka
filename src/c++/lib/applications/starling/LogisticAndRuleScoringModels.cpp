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
 *  \author Morten Kallberg
 */

#include "LogisticAndRuleScoringModels.hh"

#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"

#include <iostream>
#include <fstream>



/// the logistic regression model still requires an actual string map -- something we want
/// to move away from, but this provides a bridge for now:
///
template <typename FEATURESET>
featuremap
getFeatureMap(
    const VariantScoringFeatureKeeper<FEATURESET>& features)
{
    typedef typename FEATURESET::index_t findex_t;
    featuremap ret;
    for (unsigned index(0); index<FEATURESET::SIZE; ++index)
    {
        const findex_t findex(static_cast<findex_t>(index));
        ret.insert(std::make_pair(FEATURESET::get_feature_label(findex),features.get(findex)));
    }
    return ret;
}



LogisticAndRuleScoringModels::
LogisticAndRuleScoringModels(
    const gvcf_deriv_options& dopt,
    const std::string& model_file,
    const std::string& name) :
    _dopt(dopt)
{

    using namespace boost::algorithm;

    assert (not model_file.empty());
    assert (not name.empty());

#ifdef DEBUG_CAL
    log_os << "Loading logisitic and qrule models from file: " << model_file << "\n";
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
            _modelType=tokens.at(2);
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
                _pars[submodel][parspace][tokens.at(0)] = atof(tokens.at(1).c_str());
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
    assert(not _pars.empty());
}



void
LogisticAndRuleScoringModels::
do_site_rule_model(
    const featuremap& cutoffs,
    const GermlineDiploidSiteCallInfo& si,
    GermlineDiploidSiteSimpleGenotypeInfo& smod) const
{
    if (si.smod.gqx<cutoffs.at("GQX")) smod.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    if (cutoffs.at("DP")>0 && _dopt.is_max_depth())
    {
        if ((si.n_used_calls+si.n_unused_calls) > _dopt.max_depth) smod.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }
    // high DPFratio filter
    const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
    if (total_calls>0)
    {
        const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
        if (filt>cutoffs.at("DPFratio")) smod.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighBaseFilt);
    }
    if (si.dgt.is_snp)
    {
        if (smod.strand_bias>cutoffs.at("HighSNVSB")) smod.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighSNVSB);
    }
}



void
LogisticAndRuleScoringModels::
do_indel_rule_model(
    const featuremap& cutoffs,
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    if (cutoffs.at("GQX")>0)
    {
        if (ii.first().gqx<cutoffs.at("GQX")) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::LowGQX);
    }

    if (cutoffs.at("DP")>0 && _dopt.is_max_depth())
    {
        if (ii.first()._indelSampleReportInfo.tier1Depth > _dopt.max_depth) call.set_filter(GERMLINE_VARIANT_VCF_FILTERS::HighDepth);
    }
}



featuremap
LogisticAndRuleScoringModels::
normalize(
    const featuremap& features,
    const featuremap& adjust_factor,
    const featuremap& norm_factor) const
{
    featuremap result(features);
    for (const auto& val : norm_factor)
    {
        const auto& key(val.first);
        result[key] = (result.at(key) - adjust_factor.at(key))/norm_factor.at(key);
    }
    return result;
}



double
LogisticAndRuleScoringModels::
log_odds(
    const featuremap& features,
    const featuremap& coeffs) const
{
    //Model: ln(p(TP|x)/p(FP|x))=w_1*x_1 ... w_n*x_n + w_1*w_2*x_1 ... w_{n-1}*w_{n}*x_1
    // calculate sum from feature map

    using namespace boost::algorithm;
    std::vector<std::string> tokens;
    double sum = coeffs.at("Intercept");
    for (const auto& val : coeffs)
    {
        // check that our coefficient is different from 0
        if ((val.first == "Intercept") || (val.second ==0)) continue;
        split(tokens, val.first, is_any_of(":"));
        double term = val.second;
        for (const std::string& tok : tokens)
        {
            const auto fiter(features.find(tok));
            if (fiter == features.end()) continue;
            term *= fiter->second;
        }
        sum += term;
    }
    //TO-DO sort map to find most predictive feature for setting filter
    return sum;
}

//- From the identities logOddsRatio = ln(p(TP|x)/p(FP|x)) and p(TP|x) + p(FP|x) = 1,
//  solve for p(TP): p(TP) = 1/(1+exp(-logOddsRatio))
//- Rescale the probabilities p(TP), p(FP) with the specified class priors,
//  as they were counted on the training set before balancing to 50/50.
//- Rescale with the class priors; hard-coded 0.5, because we have 1:1 ratio in the balanced data
//  p(FP).pos.res = p(FP).pos*(minorityPrior/0.5)
//  p(TP).pos.res = p(TP).pos*(minorityPrior/0.5)
//- Renormalize to sum to one
//  p(FP).pos.ret = p(FP).pos.res/(p(TP).pos.res + p(FP).pos.res)
//  p(TP).pos.ret = p(TP).pos.res/(p(TP).pos.res + p(FP).pos.res)
//- Above simplifies to:
//  p(FP).pos.ret = p(FP).pos*minorityPrior/(1+2*minorityPrior*p(FP).pos-minorityPrior-p(FP).pos)
//- Convert the rescaled probability p(FP) into a
//  Q-score: Q-score = round( -10*log10(p(FP).pos.ret) )
static
int
prior_adjustment(
    const double raw_score,
    const double minorityPrior)
{
    const double pFP          = 1.0/(1+std::exp(raw_score)); // this calculation can likely be simplified
    const double pFPrescale   = pFP*minorityPrior/(1+2*minorityPrior*pFP-minorityPrior-pFP);
    const int qscore          = error_prob_to_qphred(pFPrescale);

    return qscore;
}



void
LogisticAndRuleScoringModels::
apply_site_qscore_filters(
    const LEGACY_CALIBRATION_MODEL::var_case var_case,
    const GermlineDiploidSiteCallInfo& si,
    GermlineDiploidSiteSimpleGenotypeInfo& smod) const
{
    // do extreme case handling better
    smod.empiricalVariantScore = std::min(smod.empiricalVariantScore,60);

    const double dpfExtreme(0.85);
    const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
    if (total_calls>0)
    {
        const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
        if (filt>dpfExtreme) smod.empiricalVariantScore=3;
    }

    if (smod.empiricalVariantScore<0)
    {
        const auto orig_filters(smod.filters);
        do_site_rule_model(_pars.at("snp").at("cutoff"), si, smod);
        if (smod.filters.count()>0)
        {
            smod.empiricalVariantScore = 1;
            smod.filters = orig_filters;
        }
        else
        {
            smod.empiricalVariantScore = 35;
        }
    }

    if (smod.empiricalVariantScore < get_var_threshold(var_case))
    {
        smod.set_filter(LEGACY_CALIBRATION_MODEL::get_Qscore_filter(var_case)); // more sophisticated filter setting here
    }
}



void
LogisticAndRuleScoringModels::
apply_indel_qscore_filters(
    const LEGACY_CALIBRATION_MODEL::var_case var_case,
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    call.empiricalVariantScore = std::min(call.empiricalVariantScore,60);

    if (call.empiricalVariantScore<0)
    {
        const auto orig_filters(call.filters);
        do_indel_rule_model(_pars.at("indel").at("cutoff"), ii, call);
        if (call.filters.count()>0)
        {
            call.empiricalVariantScore = 1;
            call.filters = orig_filters;
        }
        else
        {
            call.empiricalVariantScore = 12;
        }
    }

    if (call.empiricalVariantScore < get_var_threshold(var_case))
    {
        call.set_filter(LEGACY_CALIBRATION_MODEL::get_Qscore_filter(var_case));
    }
}

// joint logistic regression for both SNPs and INDELs
int
LogisticAndRuleScoringModels::
logistic_score(
    const LEGACY_CALIBRATION_MODEL::var_case var_case,
    const featuremap& features) const
{
    const std::string var_case_label(LEGACY_CALIBRATION_MODEL::get_label(var_case));

    // normalize
    const auto case_pars(this->_pars.at(var_case_label));
    const featuremap norm_features = normalize(features,case_pars.at("CenterVal"),case_pars.at("ScaleVal"));

    //calculates log-odds ratio
    const double raw_score = log_odds(norm_features, case_pars.at("Coefs"));

    // adjust by prior and calculate q-score
    const int Qscore = prior_adjustment(raw_score, case_pars.at("Priors").at("fp.prior"));

    return Qscore;
}

int
LogisticAndRuleScoringModels::
get_var_threshold(
    const LEGACY_CALIBRATION_MODEL::var_case& my_case) const
{
    return this->_pars.at(LEGACY_CALIBRATION_MODEL::get_label(my_case)).at("PassThreshold").at("Q");
}

bool
LogisticAndRuleScoringModels::
is_logistic_model() const
{
    return (_modelType=="LOGISTIC");
}



void
LogisticAndRuleScoringModels::
score_site_instance(
    const GermlineDiploidSiteCallInfo& si,
    GermlineDiploidSiteSimpleGenotypeInfo& smod) const
{
    if (is_logistic_model())   //case we are using a logistic regression mode
    {
        LEGACY_CALIBRATION_MODEL::var_case var_case(LEGACY_CALIBRATION_MODEL::HomSNP);
        if (si.is_het())
            var_case = LEGACY_CALIBRATION_MODEL::HetSNP;

#undef HETALTSNPMODEL
#ifdef HETALTSNPMODEL // future-proofing: do not remove unless you are sure we will not be adding hetalt SNP model to scoring
        if (si.is_hetalt())
            var_case = LEGACY_CALIBRATION_MODEL::HetAltSNP;
        //else
#endif


        // compute scoring features if this hasn't already been done:
        if (smod.features.empty())
        {
            static const bool isComputeDevelopmentFeatures(false);
            const bool isUniformDepthExpected(_dopt.is_max_depth());
            si.computeEmpiricalScoringFeatures(isUniformDepthExpected, isComputeDevelopmentFeatures, normal_depth(), smod);
        }
        smod.empiricalVariantScore = logistic_score(var_case, getFeatureMap(smod.features));
        apply_site_qscore_filters(var_case, si, smod);
    }
//    else if(this->_modelType=="RFtree"){ // place-holder, put random forest here
//        si.Qscore = rf_score(var_case, features);
//    }
    else if (this->_modelType=="RULE")   //case we are using a rule based model
    {
        do_site_rule_model(_pars.at("snp").at("cutoff"), si, smod);
    }
    else
    {
        assert(false && "Unknown model type");
    }
}



void
LogisticAndRuleScoringModels::
score_indel_instance(
    const GermlineDiploidIndelCallInfo& ii,
    GermlineDiploidIndelSimpleGenotypeInfo& call) const
{
    if (is_logistic_model())   //case we are using a logistic regression mode
    {
        LEGACY_CALIBRATION_MODEL::var_case var_case(LEGACY_CALIBRATION_MODEL::HetDel);
        switch (ii.first()._indelReportInfo.it)
        {
        case INDEL::DELETE:
        {
            if (ii.is_hetalt())
                var_case = LEGACY_CALIBRATION_MODEL::HetAltDel;
            else if (! ii.is_het())
                var_case = LEGACY_CALIBRATION_MODEL::HomDel;
        }
        break;
        case INDEL::INSERT:
        {
            if (ii.is_hetalt())
                var_case = LEGACY_CALIBRATION_MODEL::HetAltIns;
            else if (ii.is_het())
                var_case = LEGACY_CALIBRATION_MODEL::HetIns;
            else
                var_case = LEGACY_CALIBRATION_MODEL::HomIns;
        }
        break;
        default:
            // block substitutions???
            do_indel_rule_model(_pars.at("indel").at("cutoff"), ii, call);
            return;
        }

        // compute scoring features if this hasn't already been done:
        if (call.features.empty())
        {
            static const bool isComputeDevelopmentFeatures(false);
            const bool isUniformDepthExpected(_dopt.is_max_depth());
            call.computeEmpiricalScoringFeatures(isUniformDepthExpected, isComputeDevelopmentFeatures, normal_depth(), ii.is_hetalt());
        }
        call.empiricalVariantScore = logistic_score(var_case, getFeatureMap(call.features));
        apply_indel_qscore_filters(var_case, ii, call);
    }
    else if (this->_modelType=="RULE")   //case we are using a rule based model
    {
        do_indel_rule_model(this->_pars.at("indel").at("cutoff"), ii, call);
    }
    else
    {
        assert(false && "Unknown model type");
    }
}

// give the normal depth for this model; used to depth-normalize various features in logistic model
double LogisticAndRuleScoringModels::normal_depth() const
{
    return this->_dopt.norm_depth;
}

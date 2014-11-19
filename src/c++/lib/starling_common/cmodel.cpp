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
 * cmodel.cpp
 *
 *      Created on: Dec 1, 2013
 *      Author: Morten Kallberg
 */

#include "cmodel.hh"
#include <stdlib.h>     /* atof */
#include <iostream>
#include <sstream>
#include <cassert>
#include <string>
#include <fstream>
#include <iterator>
//#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "blt_util/qscore.hh"

//#define DEBUG_MODEL

#ifdef DEBUG_MODEL
#include "blt_util/log.hh"
#endif

// add model paramaters
void c_model::add_parameters(const parmap& myPars)
{
    this->pars = myPars;
}

// rule-based filtering for SNPs
void c_model::do_rule_model(featuremap& cutoffs, site_info& si)
{
    if (si.smod.gqx<cutoffs["GQX"]) si.smod.set_filter(VCF_FILTERS::LowGQX);
    if (cutoffs["DP"]>0 && dopt.is_max_depth())
    {
        if ((si.n_used_calls+si.n_unused_calls) > dopt.max_depth) si.smod.set_filter(VCF_FILTERS::HighDepth);
    }
    // high DPFratio filter
    const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
    if (total_calls>0)
    {
        const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
        if (filt>cutoffs["DPFratio"]) si.smod.set_filter(VCF_FILTERS::HighBaseFilt);
    }
    if (si.dgt.is_snp)
    {
        if (si.dgt.sb>cutoffs["HighSNVSB"]) si.smod.set_filter(VCF_FILTERS::HighSNVSB);
    }
}

// rule-based filtering for INDELs
void c_model::do_rule_model(featuremap& cutoffs, indel_info& ii)
{
    ii.imod.max_gt=ii.dindel.max_gt_poly;
    ii.imod.gq=ii.dindel.max_gt_poly_qphred;

    if (cutoffs["GQX"]>0)
    {
        if (ii.imod.gqx<cutoffs["GQX"]) ii.imod.set_filter(VCF_FILTERS::LowGQX);
    }

    if (cutoffs["DP"]>0 && dopt.is_max_depth())
    {
        if (ii.isri.depth > dopt.max_depth) ii.imod.set_filter(VCF_FILTERS::HighDepth);
    }
}

//Transform the features with the specified scaling parameters that were used to standardize
//the dataset to zero mean and unit variance: newVal = (oldVal-centerVal)/scaleVal.
featuremap c_model::normalize(featuremap features, featuremap& adjust_factor, featuremap& norm_factor)
{
    for (featuremap::const_iterator it = norm_factor.begin(); it != norm_factor.end(); ++it)   // only normalize the features that are needed
    {
//        log_os << it->first << "=" << features[it->first] << " ";
        features[it->first] = (features[it->first]-adjust_factor[it->first])/norm_factor[it->first];
//        log_os << it->first << "=" << features[it->first] << " ";
    }
//    log_os << "\n";
    return features;
}

//Model: ln(p(TP|x)/p(FP|x))=w_1*x_1 ... w_n*x_n + w_1*w_2*x_1 ... w_{n-1}*w_{n}*x_1
// calculate sum from feature map
double c_model::log_odds(featuremap features, featuremap& coeffs)
{
    using namespace boost::algorithm;
    std::vector<std::string> tokens;
    double sum = coeffs["Intercept"];
//    log_os << "sum" << "=" << sum << "\n";
    for (featuremap::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it)
    {
        if (it->first !="Intercept" && it->second !=0)   // check that our coefficient is different from 0
        {
            split(tokens, it->first, is_any_of(":"));
//            log_os << it->first << "=" << it->second << "\n";
            double term = it->second;
            for (unsigned int i=0; i < tokens.size(); i++)
            {
                if (features.find( tokens[i] ) != features.end())
                {
                    term = term*features[tokens[i]];
//                    log_os << tokens[i] << "=" << features[tokens[i]] << "\n";
                }
                //should not get here, if we havent loaded the feature we are in trouble...
//                else{
//                    log_os << "I dont know feature " << tokens[i] << "\n";
//                }
            }
//            log_os << "\n";
//            log_os << "term" << "=" << term << "\n";
            // use term to determine the most predictive parameter
            sum += term;
//            log_os << "sum " << "=" << sum << "\n";
//            log_os << tokens.size() << "\n";
        }
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
int prior_adjustment(
    const double raw_score,
    const double minorityPrior)
{

    double pFP          = 1.0/(1+std::exp(raw_score)); // this calculation can likely be simplified
    double pFPrescale   = pFP*minorityPrior/(1+2*minorityPrior*pFP-minorityPrior-pFP);
    int qscore          = error_prob_to_qphred(pFPrescale);
#ifdef DEBUG_MODEL
//        log_os << "minorityPrior " << minorityPrior << "\n";
//        log_os << "raw_score=" << raw_score << "\n";
//        log_os << "rescale=" << pFPrescale << "\n";
//        log_os << "experimental=" << qscore << "\n";
#endif

    return qscore;
}
void c_model::apply_qscore_filters(site_info& si, const int qscore_cut, const CALIBRATION_MODEL::var_case my_case)   //, featuremap& most_predictive) {
{
//    most_predictive.size();
#ifdef DEBUG_MODEL
    log_os << "Qscore "<< si.Qscore << std::endl;
#endif

    // do extreme case handeling better
    if (si.Qscore>60)
        si.Qscore = 60;

    const double dpfExtreme(0.85);
    const unsigned total_calls(si.n_used_calls+si.n_unused_calls);
    if (total_calls>0)
    {
        const double filt(static_cast<double>(si.n_unused_calls)/static_cast<double>(total_calls));
        if (filt>dpfExtreme) si.Qscore=3;
    }

    if (si.Qscore<0)
    {
        featuremap cutoffs = {{"GQX", 30}, {"DP", 1}, {"DPFratio", 0.4}, {"HighSNVSB", 10}};
        this->do_rule_model(cutoffs,si);
        if (si.smod.filters.count()>0)
        {
            si.Qscore = 1;
            si.smod.filters.reset();
        }
        else
        {
            si.Qscore = 35;
        }
    }

    if (si.Qscore < qscore_cut)
    {
//        log_os << CALIBRATION_MODEL::get_label(my_case) << "\n";
        si.smod.set_filter(CALIBRATION_MODEL::get_Qscore_filter(my_case)); // more sophisticated filter setting here
    }
}

void c_model::apply_qscore_filters(indel_info& ii, const int qscore_cut, const CALIBRATION_MODEL::var_case my_case)   //, featuremap& most_predictive) {
{
    // do extreme case handling better
    if (ii.Qscore>60)
        ii.Qscore = 60;

    if (ii.Qscore<0)
    {
        featuremap cutoffs = {{"GQX", 30}, {"DP", 1},{"DPFratio", 0.2}};
        this->do_rule_model(cutoffs,ii);
        if (ii.imod.filters.count()>0)
        {
            ii.Qscore = 1;
            ii.imod.filters.reset();
        }
        else
        {
            ii.Qscore = 12;
        }
    }

    if (ii.Qscore < qscore_cut)
    {
//        log_os << CALIBRATION_MODEL::get_label(my_case) << "\n";
        ii.imod.set_filter(CALIBRATION_MODEL::get_Qscore_filter(my_case));
    }
}

// joint logistic regression for both SNPs and INDELs
int c_model::logistic_score(const CALIBRATION_MODEL::var_case var_case, featuremap features)
{
    std::string var_case_label(CALIBRATION_MODEL::get_label(var_case));
    // normalize
    featuremap norm_features = this->normalize(features,this->pars[var_case_label]["CenterVal"],this->pars[var_case_label]["ScaleVal"]);

    //calculates log-odds ratio
    double raw_score = this->log_odds(norm_features,this->pars[var_case_label]["Coefs"]);

    // adjust by prior and calculate q-score
    int Qscore = prior_adjustment(raw_score,this->pars[var_case_label]["Priors"]["fp.prior"]);

    // not active, but a good sanity-check:
    // assert(Qscore>=0);

    // not active, but possibly the safe thing to do, to avoid special-case treatment elsewhere that we no longer want:
    // if(Qscore == 0) Qscore = 1;

    return Qscore;
}

int c_model::get_var_threshold(CALIBRATION_MODEL::var_case& my_case)
{
    return this->pars[CALIBRATION_MODEL::get_label(my_case)]["PassThreshold"]["Q"];
}

bool c_model::is_logistic_model() const
{
    return this->model_type=="LOGISTIC";
}

//score snp case
void c_model::score_instance(featuremap features, site_info& si)
{
    if (this->model_type=="LOGISTIC")   //case we are using a logistic regression mode
    {
        CALIBRATION_MODEL::var_case var_case(CALIBRATION_MODEL::HomSNP);
#undef HETALTSNPMODEL
#ifdef HETALTSNPMODEL // future-proofing: do not remove unless you are sure we will not be adding hetalt SNP model to VQSR
        if (si.is_hetalt())
            var_case = CALIBRATION_MODEL::HetAltSNP;
        else
#endif
        if (si.is_het())
            var_case = CALIBRATION_MODEL::HetSNP;

#ifdef DEBUG_MODEL
        //log_os << "Im doing a logistic model varcase: " << var_case <<  "\n";
#endif
        si.Qscore = logistic_score(var_case, features);
        this->apply_qscore_filters(si,static_cast<int>(this->pars[CALIBRATION_MODEL::get_label(var_case)]["PassThreshold"]["Q"]), var_case); // set filters according to q-scores
    }
//    else if(this->model_type=="RFtree"){ // place-holder, put random forest here
//        si.Qscore = rf_score(var_case, features);
//    }
    else if (this->model_type=="RULE")   //case we are using a rule based model
    {
        this->do_rule_model(this->pars["snp"]["cutoff"],si);
    }
}

// score indel case
void c_model::score_instance(featuremap features, indel_info& ii)
{
    if (this->model_type=="LOGISTIC")   //case we are using a logistic regression mode
    {
        //TODO put into enum context
        CALIBRATION_MODEL::var_case var_case(CALIBRATION_MODEL::HetDel);
        if (ii.iri.it==INDEL::DELETE)
        {
          if (ii.is_hetalt())
            var_case = CALIBRATION_MODEL::HetAltDel;
          else if (! ii.is_het())
            var_case = CALIBRATION_MODEL::HomDel;
        }
        else if (ii.iri.it==INDEL::INSERT)
        {
            if (ii.is_hetalt())
                var_case = CALIBRATION_MODEL::HetAltIns;
            else if (ii.is_het())
                var_case = CALIBRATION_MODEL::HetIns;
            else
                var_case = CALIBRATION_MODEL::HomIns;
        }

        else
        {
          // block substitutions???
          this->do_rule_model(this->pars["indel"]["cutoff"],ii);
          return;
        }
        ii.Qscore = logistic_score(var_case, features);
        this->apply_qscore_filters(ii,static_cast<int>(this->pars[CALIBRATION_MODEL::get_label(var_case)]["PassThreshold"]["Q"]),var_case);
    }
//    else if(this->model_type=="RFtree"){
//        si.Qscore = rf_score(var_case, features);
//    }
    else if (this->model_type=="RULE")   //case we are using a rule based model
    {
        this->do_rule_model(this->pars["indel"]["cutoff"],ii);
    }
}

// give the normal depth for this model; used to depth-normalize various features in logistic model
double c_model::normal_depth() const { return this->dopt.norm_depth; }

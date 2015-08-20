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
 * Indelmodel.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: Morten Kallberg
 */

#include "calibration/Indelmodel.hh"
#include "blt_util/log.hh"

#include <cmath>
#include <cassert>


// Calculate p(error) of
//    CASE: del
//    FIT pars: [  1.49133831e-03   1.03348683e+01   1.13646811e+00   1.18488282e-05]
//    Function prob(error)=0.00149133830825/ (1 + exp((10.3348683003-x)/1.13646810558))+1.18488281756e-05
//    --------------------
//    CASE: ins
//    FIT pars: [  1.09573511e-03   9.82226042e+00   1.03579658e+00   8.31843836e-06]
//    Function prob(error)=0.00109573511176/ (1 + exp((9.82226041538-x)/1.03579658224))+8.31843836296e-06
//    --------------------
Indel_model
generate_new()
{
    Indel_model res("new","v1","",1,40);

    const double insert_A(1.49133831e-03);
    const double insert_B(1.03348683e+01);
    const double insert_C(1.13646811e+00);
    const double insert_D(1.18488282e-05);

    const double delete_A(1.09573511e-03);
    const double delete_B(9.82226042e+00);
    const double delete_C(1.03579658e+00);
    const double delete_D(8.31843836e-06);

    for (unsigned hpol_len=1; hpol_len <= res.MaxTractLength; ++hpol_len)
    {
        double insert_error_prob, delete_error_prob;

        const double insert_g(insert_A/ (1 + std::exp((insert_B-hpol_len)/insert_C))+insert_D);
        insert_error_prob=(1.-std::exp(-insert_g/hpol_len));
        const double delete_g(delete_A/ (1 + std::exp((delete_B-hpol_len)/delete_C))+delete_D);
        delete_error_prob=(1.-std::exp(-delete_g/hpol_len));

        std::pair<double,double> pair(insert_error_prob,delete_error_prob);
        res.add_prop(0,hpol_len-1,pair);
    }
    return res;
}


Indel_model
generate_old()
{
    Indel_model res("old","v1","",1,40);

    const double insert_A(5.03824e-7);
    const double insert_B(3.30572e-10);
    const double insert_C(6.99777);

    const double delete_hpol1_err(5.00057e-5);
    const double delete_A(1.09814e-5);
    const double delete_B(5.19742e-10);
    const double delete_C(6.99256);
    for (unsigned hpol_len=1; hpol_len <= res.MaxTractLength; ++hpol_len)
    {
        double insert_error_prob, delete_error_prob;

        const double insert_g(insert_A*hpol_len+insert_B*std::pow(hpol_len,insert_C));
        insert_error_prob=(1.-std::exp(-insert_g));

        double delete_g(delete_hpol1_err);
        if (hpol_len>1)
        {
            delete_g = delete_A*hpol_len+delete_B*std::pow(hpol_len,delete_C);
        }
        delete_error_prob=(1.-std::exp(-delete_g));

        std::pair<double,double> pair(insert_error_prob,delete_error_prob);
        res.add_prop(0,hpol_len-1,pair);
    }
    return res;
}



Indel_model::Indel_model()
{
    this->MaxMotifLength = 0;
    this->MaxTractLength = 0;
}

// read in the json
void Indel_model::Deserialize( const Json::Value& root)
{
    serialized_model::Deserialize(root);

    this->MaxMotifLength = root["MaxMotifLength"].asInt();
    this->MaxTractLength = root["MaxTractLength"].asInt();
    Json::Value jmodels = root["Model"];

    // Reads in model parameter matrix with entries as error-pair [ins_error,del_error]
    // in the following format:
    // unit length 1: [[[ins_hpol1,del_hpol1],[ins_hpol2,del_hpol2],...,[ins_hpol_m,del_hpol_m]],
    // unit length 2:  [ins_dinuc1,del_dinuc1],[ins_dinuc2,del_dinuc2],...,[ins_dinuc_m,del_dinuc_m]],
    //  ....
    // unit length N:  [[ins_repeatN,del_repeatN],[ins_repeatN2,del_repeatN2],...,]]
    unsigned tract(0);
    unsigned unit(0);
    for ( unit = 0; unit < jmodels.size(); ++unit)
        for ( tract = 0; tract < jmodels[unit].size(); ++tract)
        {
            std::pair<double,double> pair(jmodels[unit][tract][1].asDouble(),jmodels[unit][tract][0].asDouble());
            // assert((unit == 0 || ((tract + 1) < (2 * (unit + 1)) && pair.first == 0 && pair.second == 0) ||
            // 	((tract + 1) >= (2 * (unit + 1)))) &&
            // 	"Nonzero error probability for tract length below twice repeat unit size");
            this->add_prop(unit,tract,pair);
        }

    //Make sure the model is self-consistent -- contains the data
    // specified in Max{Motif,Tract}
    assert(unit==this->MaxMotifLength && "Unexpected motif length in indel model");
    assert(tract>=this->MaxTractLength && "Unexpected tract length in indel model");
}

void Indel_model::calc_prop(const starling_base_options& client_opt,
                            const starling_indel_report_info& iri,
                            double& indel_error_prob,
                            double& ref_error_prob) const
{
    bool use_length_dependence = false;
    calc_prop(client_opt,iri,indel_error_prob,ref_error_prob,use_length_dependence);
}

std::pair<double,double>
Indel_model::calc_prop(const starling_base_options& client_opt,
                       const starling_indel_report_info& iri)
{
    if (client_opt.acov_alpha>0 && iri.ihpol>0) {}
    return model[0][0];
}

void Indel_model::calc_prop(const starling_base_options& client_opt,
                            const starling_indel_report_info& iri,
                            double& indel_error_prob,
                            double& ref_error_prob,
                            bool use_length_dependence) const
{

    // determine simple case
    const bool is_simple_indel(iri.it==INDEL::INSERT || iri.it==INDEL::DELETE);

    // determine the tract length to use
    static const unsigned one(1);
    const unsigned repeat_unit    = std::min(std::max(iri.repeat_unit_length,one), this->MaxMotifLength);
    const unsigned ref_hpol_len   = std::min(repeat_unit*std::max(iri.ref_repeat_count,one),this->MaxTractLength);
    const unsigned indel_hpol_len = std::min(repeat_unit*std::max(iri.indel_repeat_count,one),this->MaxTractLength);

    // determine indel size
    int indel_size(1);
    if (use_length_dependence)
    {
        indel_size = std::abs(static_cast<long>(iri.ref_repeat_count)-static_cast<long>(iri.indel_repeat_count));
    }

    // Fall-back probs
    unsigned min_tract_length = get_min_tract_length(iri);
    if (repeat_unit == 1)
    {
        min_tract_length = 1;
    }

    if (! is_simple_indel)
    {
        // breakpoints and swaps -- // use zero repeat error for now.
        // TODO - provide estimates for complex indels NOTE: likely never utilized
        double baseline_ins_prob(this->model[0][0].first);
        double baseline_del_prob(this->model[0][0].second);
        // double baseline_ins_prob(this->model[repeat_unit - 1][min_tract_length - 1].first);
        // double baseline_del_prob(this->model[repeat_unit - 1][min_tract_length - 1].second);

        indel_error_prob=std::max(baseline_ins_prob,baseline_del_prob);
        ref_error_prob=indel_error_prob;
        return;
    }

    // if tract length is too short for repeat unit, set to shortest indel error rate for
    // that repeat unit length
    const unsigned ref_query_len = std::max(min_tract_length, ref_hpol_len);
    const unsigned indel_query_len = std::max(min_tract_length, indel_hpol_len);

    if (iri.repeat_unit_length <= MaxMotifLength)
    {
        if (iri.it == INDEL::INSERT)
        {
            indel_error_prob=std::max(model[0][0].first,
                                      std::pow(model[repeat_unit - 1][ref_query_len - 1].first,indel_size));

            // Reverse prob that true allele has been masked as reference by chance,
            // may want to leave this term for now.
            ref_error_prob=client_opt.indel_ref_error_factor
                           * std::max(model[0][0].second,
                                      std::pow(model[repeat_unit - 1][indel_query_len - 1].second,indel_size));
        }
        else if (iri.it == INDEL::DELETE)
        {
            indel_error_prob=std::max(model[0][0].second,
                                      std::pow(model[repeat_unit - 1][ref_query_len - 1].second,indel_size));

            ref_error_prob=client_opt.indel_ref_error_factor
                           * std::max(model[0][0].first,
                                      std::pow(model[repeat_unit - 1][indel_query_len - 1].first,indel_size));
        }
        else
        {
            // this should never happen, but just for completeness
            log_os << "ERROR: Unknown indel type: " << iri.desc << "\n";
            throw blt_exception("Unknown indel type.");
        }
    }
    else
    {
        // if there is no model for the repeat unit length observed, and the indel is in
        // non-repeat sequence (i.e. RC=0/IC=1 or RC=1/IC=0)
        if (iri.it == INDEL::INSERT)
        {
            // current model is too aggressive for hpol length 1, which will be fixed in
            // the model calculation in the very near future.  For now, use error prob
            // from the previous model
            indel_error_prob = model[0][0].first;
            ref_error_prob   = model[0][0].second;
        }
        else if (iri.it == INDEL::DELETE)
        {
            indel_error_prob = model[0][0].second;
            ref_error_prob   = model[0][0].first;
        }
    }
    // else
    // {
    //     log_os << "ERROR: Unknown indel type: " << iri.desc << "\n";
    //     throw blt_exception("Unknown indel type.");
    // }

}

unsigned Indel_model::get_min_tract_length(const starling_indel_report_info& iri) const
{
    return iri.repeat_unit_length * 2;
}

bool Indel_model::is_simple_tandem_repeat(const starling_indel_report_info& iri) const
{
    // an STR only has insertions or deletions, has a repeat unit length present in the model,
    // and has a tract length present in the model
    unsigned min_tract_length = get_min_tract_length(iri);
    if (iri.repeat_unit_length <= MaxMotifLength &&
        (iri.it == INDEL::DELETE || iri.it == INDEL::INSERT) &&
        (iri.ref_repeat_count >= min_tract_length ||
         iri.indel_repeat_count >= min_tract_length))
    {
        return true;
    }
    return false;
}

void Indel_model::add_prop(const unsigned& unit, const unsigned& tract, const std::pair<double,double>& myProps)
{
    this->model[unit][tract] = myProps;
}

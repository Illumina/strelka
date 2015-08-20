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
 * Indelmodel.hh
 *
 *  Created on: Jun 23, 2015
 *      Author: mkallberg
 */

#include "calibration/SerializedModel.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_indel_report_info.hh"
#include "blt_util/blt_exception.hh"

#ifndef C___LIB_CALIBRATION_INDELMODEL_HH_
#define C___LIB_CALIBRATION_INDELMODEL_HH_

static const unsigned max_unit_len(10);
static const unsigned max_tract_len(40);
typedef std::pair<double,double> error_model[max_unit_len][max_tract_len];


struct Indel_model : public serialized_model
{
    Indel_model();
    Indel_model(std::string n, std::string v, std::string d,int motif, int tract):
        MaxMotifLength (motif),
        MaxTractLength (tract)
    {
        this->name 		= n;
        this->version 	= v;
        this->date 		= d;
    }

    void Deserialize(const Json::Value& root);
    std::pair<double,double> calc_prop(const starling_base_options& client_opt, const starling_indel_report_info& iri);
    void calc_prop(const starling_base_options& client_opt,
                   const starling_indel_report_info& iri,
                   double& indel_error_prob,
                   double& ref_error_prob) const;
    void calc_prop(const starling_base_options& client_opt,
                   const starling_indel_report_info& iri,
                   double& indel_error_prob,
                   double& ref_error_prob,
                   bool use_length_dependence) const;
    std::pair<double,double> get_prop(const unsigned& unit, const unsigned& tract)
    {
        return model[unit][tract];
    }

    unsigned get_max_motif_length() const
    {
        return MaxMotifLength;
    }
    unsigned get_min_tract_length(const starling_indel_report_info& iri) const;
    bool is_simple_tandem_repeat(const starling_indel_report_info& iri) const;
    void add_prop(const unsigned& unit, const unsigned& tract, const std::pair<double,double>& myProps);
    error_model model;
    unsigned MaxMotifLength, MaxTractLength;
};


// Calculate p(error) of
//    CASE: del
//    FIT pars: [  1.49133831e-03   1.03348683e+01   1.13646811e+00   1.18488282e-05]
//    Function prob(error)=0.00149133830825/ (1 + exp((10.3348683003-x)/1.13646810558))+1.18488281756e-05
//    --------------------
//    CASE: ins
//    FIT pars: [  1.09573511e-03   9.82226042e+00   1.03579658e+00   8.31843836e-06]
//    Function prob(error)=0.00109573511176/ (1 + exp((9.82226041538-x)/1.03579658224))+8.31843836296e-06
//    --------------------

static Indel_model generate_new()
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

static Indel_model generate_old()
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


#endif /* C___LIB_CALIBRATION_INDELMODEL_HH_ */

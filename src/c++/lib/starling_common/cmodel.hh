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
 * cmodel.hh
 *
 *  Created on: Jan 15, 2014
 *      Author: Morten Kallberg
 */

#ifndef CMODEL_HH_
#define CMODEL_HH_
#include <vector>
#include <map>
#include "starling_common/gvcf_locus_info.hh"

namespace CALIBRATION_MODEL
{

enum var_case
{
    HetSNP,
    HomSNP,
    HetIns,
    HomIns,
    HetDel,
    HomDel
};

inline
const char*
get_label(const unsigned idx)
{
    switch (idx)
    {
    case HetSNP:
        return "snphet";
    case HomSNP:
        return "snphom";
    case HetIns:
        return "inshet";
    case HomIns:
        return "inshom";
    case HetDel:
        return "delhet";
    case HomDel:
        return "delhom";
    default:
        assert(0);
        return NULL;
    }
}

inline
const char*
get_label_header(const unsigned idx)
{
    switch (idx)
    {
    case HetSNP:
        return "het SNP";
    case HomSNP:
        return "hom SNP";
    case HetIns:
        return "het insertion";
    case HomIns:
        return "hom insertion";
    case HetDel:
        return "het deletion";
    case HomDel:
        return "hom deletion";
    default:
        assert(0);
        return NULL;
    }
}

inline
VCF_FILTERS::index_t
get_Qscore_filter(const unsigned var_case)
{
    switch (var_case)
    {
    case HetSNP:
        return VCF_FILTERS::LowQscoreHetSNP;
    case HomSNP:
        return VCF_FILTERS::LowQscoreHomSNP;
    case HetIns:
        return VCF_FILTERS::LowQscoreHetIns;
    case HomIns:
        return VCF_FILTERS::LowQscoreHomIns;
    case HetDel:
        return VCF_FILTERS::LowQscoreHetDel;
    case HomDel:
        return VCF_FILTERS::LowQscoreHomDel;
    default:
        assert(0);
        return VCF_FILTERS::LowGQX;
    }
}
}

typedef std::map<std::string, double> featuremap;
typedef std::map<std::string, std::map<std::string, featuremap > > parmap;
class c_model
{
public:
    c_model(
        const std::string& name,
        const std::string& type) :
        model_name(name),
        model_type(type)
    {}
    // add parameters to the model
    void add_parameters(const parmap& myPars);
    void score_instance(featuremap features, site_info& si);
    void score_instance(featuremap features, indel_info& ii);
private:
    int logistic_score(const CALIBRATION_MODEL::var_case var_case, featuremap features);
    void do_rule_model(featuremap& cutoffs, site_info& si);  //snp case
    void do_rule_model(featuremap& cutoffs, indel_info& ii); //indel case
//    void do_rule_model(featuremap& cutoffs, phased_info& ii); //phased record TODO
    featuremap normalize(featuremap features, featuremap& adjust_factor, featuremap& norm_factor);
    double log_odds(featuremap features, featuremap& coeffs);
    void apply_qscore_filters(site_info& si, const int qscore_cut, const CALIBRATION_MODEL::var_case my_case);
    void apply_qscore_filters(indel_info& ii, const int qscore_cut,const CALIBRATION_MODEL::var_case my_case);
    void sanity_check();
    std::string model_name;
    std::string model_type;
    parmap pars;
};
#endif /* CMODEL_HH_ */

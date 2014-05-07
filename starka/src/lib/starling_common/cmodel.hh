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

typedef std::map<std::string, double> featuremap;
typedef std::map<std::string, std::map<std::string, featuremap > > parmap;
class c_model {
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
    int logistic_score(std::string var_case, featuremap features);
    void do_rule_model(featuremap& cutoffs, site_info& si);  //snp case
    void do_rule_model(featuremap& cutoffs, indel_info& ii); //indel case
//    void do_rule_model(featuremap& cutoffs, phased_info& ii); //phased record TODO
    featuremap normalize(featuremap features, featuremap& adjust_factor, featuremap& norm_factor);
    double log_odds(featuremap features, featuremap& coeffs);
    void apply_qscore_filters(site_info& si, const int qscore_cut);
    void apply_qscore_filters(indel_info& ii, const int qscore_cut);
    void sanity_check();
    std::string model_name;
    std::string model_type;
    parmap pars;
};
#endif /* CMODEL_HH_ */

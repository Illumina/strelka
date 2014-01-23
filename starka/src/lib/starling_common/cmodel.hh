/*
 * cmodel.hh
 *
 *  Created on: Jan 15, 2014
 *      Author: mkallberg
 */

#ifndef CMODEL_HH_
#define CMODEL_HH_
#include "blt_util/blt_exception.hh"
#include <vector>
#include <map>
#include "starling_common/gvcf_locus_info.hh"

typedef std::map<std::string, double> featuremap;
typedef std::map<std::string, std::map<std::string, featuremap > > parmap;
class c_model {
public:
    c_model(std::string name, std::string type);
    virtual ~c_model();

    // add parameters to the model
    void add_parameters(parmap myPars);
    void score_instance(featuremap features, site_info& si);
private:
    void do_rule_model(featuremap& cutoffs, site_info& si);
    featuremap normalize(featuremap features, featuremap& adjust_factor, featuremap& norm_factor);
    double log_odds(featuremap features, featuremap& coeffs);
    double prior_adjustment(const double raw_score,featuremap& priors);
    void set_qscore_filters(const double raw_score,featuremap& priors);
    std::string model_name;
    std::string model_type;
    parmap pars;
};

#endif /* CMODEL_HH_ */

/*
 * calbrationmodels.h
 *
 *  Created on: Oct 10, 2013
 *  Author: Morten Kallberg
 */

#ifndef CALBRATIONMODELS_H_
#define CALBRATIONMODELS_H_
#include "blt_util/blt_exception.hh"
#include "starling_common/gvcf_block_site_record.hh"
#include "starling_common/gvcf_locus_info.hh"
#include "starling_common/cmodel.hh"
//#include "starling_common/gvcf_aggregator.hh"


struct gvcf_deriv_options {
    gvcf_deriv_options()
        : is_max_depth(false)
        , max_depth(0)
    {}
    bool is_max_depth;
    double max_depth;
};
//forward declaration
struct gvcf_options;

class calibration_models {
public:
    calibration_models();
    virtual ~calibration_models();

    void set_model(const std::string& name);  // set the calibration model to use
    void load_models(std::string model_file);                       // read in model parameters

    void clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, site_info& si);

    c_model get_model(std::string name);

    // mimics behaviour of previous hard filters
    void default_clasify_site(const gvcf_options& opt, const gvcf_deriv_options& dopt, site_info& si);
private:
    typedef std::map<std::string,c_model> modelmap;
    typedef std::map<std::string, double> featuremap;
    std::string model_name;
    modelmap models;
};

#endif /* CALIBRATIONMODELS_H_ */

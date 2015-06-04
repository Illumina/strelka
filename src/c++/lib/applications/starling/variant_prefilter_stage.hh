/*
 *
 *  Created on: Jun 4, 2015
 *      Author: jduddy
 */
#pragma once
#include "variant_pipe_stage.hh"

class calibration_models;

class variant_prefilter_stage : public variant_pipe_stage
{
public:
    variant_prefilter_stage(const calibration_models& model, variant_pipe_stage& destination);
    void process(site_info& si) override;
    void process(indel_info& ii) override;

    static void  add_site_modifiers(
        const site_info& si,
        site_modifiers& smod,
        const calibration_models& model);

private:
    const calibration_models& _model;
};



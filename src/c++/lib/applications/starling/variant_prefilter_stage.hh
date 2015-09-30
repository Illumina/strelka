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
 *
 *  Created on: Jun 4, 2015
 *      Author: jduddy
 */
#pragma once
#include "variant_pipe_stage_base.hh"

class calibration_models;

class variant_prefilter_stage : public variant_pipe_stage_base
{
public:
    variant_prefilter_stage(const calibration_models& model, std::shared_ptr<variant_pipe_stage_base> destination);
    void process(std::unique_ptr<site_info> si) override;
    void process(std::unique_ptr<indel_info> ii) override;

    static void  add_site_modifiers(
        const digt_site_info& si,
        digt_call_info& smod,
        const calibration_models& model);

private:
    const calibration_models& _model;
};



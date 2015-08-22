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

#pragma once

#include "calibration/SerializedModel.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_indel_report_info.hh"
#include "blt_util/blt_exception.hh"

static const unsigned max_unit_len(10);
static const unsigned max_tract_len(40);
typedef std::pair<double,double> error_model[max_unit_len][max_tract_len];


struct IndelErrorModel : public serialized_model
{
    IndelErrorModel();
    IndelErrorModel(std::string n, std::string v, std::string d,int motif, int tract):
        MaxMotifLength (motif),
        MaxTractLength (tract)
    {
        this->name 		= n;
        this->version 	= v;
        this->date 		= d;
    }

    void Deserialize(const Json::Value& root);

    void calc_prop(const starling_base_options& client_opt,
                   const starling_indel_report_info& iri,
                   double& indel_error_prob,
                   double& ref_error_prob,
                   bool use_length_dependence = false) const;

    bool is_simple_tandem_repeat(const starling_indel_report_info& iri) const;

    unsigned get_max_motif_length() const
    {
        return MaxMotifLength;
    }

    void add_prop(const unsigned& unit, const unsigned& tract, const std::pair<double,double>& myProps);

private:
    static
    unsigned
    get_min_tract_length(const starling_indel_report_info& iri)
    {
        return iri.repeat_unit_length * 2;
    }

public:
    error_model model;
    unsigned MaxMotifLength, MaxTractLength;
};


IndelErrorModel generate_new_indel_error_model();

IndelErrorModel generate_old_indel_error_model();

// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2016 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
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

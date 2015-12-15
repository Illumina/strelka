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

#include "IndelModelMetadata.hh"
#include "starling_common/starling_base_shared.hh"
#include "starling_common/starling_indel_report_info.hh"
#include "blt_util/blt_exception.hh"

static const unsigned max_unit_len(10);
static const unsigned max_tract_len(40);

struct indel_error_rates{
    double insert_rate;
    double delete_rate;

    indel_error_rates()
    {
        insert_rate = 0;
        delete_rate = 0;
    }
    indel_error_rates(double insert_error_rate, double delete_error_rate)
    {
        insert_rate = insert_error_rate;
        delete_rate = delete_error_rate;
    }
};

typedef indel_error_rates error_model[max_unit_len][max_tract_len];

struct IndelErrorModel
{
    IndelErrorModel();
    IndelErrorModel(const std::string& n, const std::string& v, const std::string& d,int motif, int tract):
        MaxMotifLength (motif),
        MaxTractLength (tract)
    {
        _meta.name 		= n;
        _meta.version 	= v;
        _meta.date 		= d;
    }

    const std::string&
    getName() const
    {
        return _meta.name;
    }

    void Deserialize(const Json::Value& root);

    void calc_prop(const starling_base_options& client_opt,
                   const starling_indel_report_info& iri,
                   double& indel_error_prob,
                   double& ref_error_prob,
                   bool use_length_dependence = false) const;

    void calc_abstract_prop(unsigned repeat_unit_length,
                            unsigned tract_length,
                            unsigned indel_size,
                            indel_error_rates& error_rates,
                            bool use_length_dependence) const;

    indel_error_rates calc_abstract_prop(unsigned repeat_unit_length,
                            unsigned tract_length,
                            unsigned indel_size,
                            bool use_length_dependence) const;

    bool is_simple_tandem_repeat(const starling_indel_report_info& iri) const;

    unsigned get_max_motif_length() const
    {
        return MaxMotifLength;
    }

    double adjusted_rate(unsigned repeat_unit_length,
                         unsigned tract_length,
                         unsigned indel_size,
                         INDEL::index_t it) const;

    void add_prop(const unsigned& unit, const unsigned& tract, const indel_error_rates& myProps);

private:
    static
    unsigned
    get_min_tract_length(const starling_indel_report_info& iri)
    {
        return iri.repeat_unit_length * 2;
    }

    IndelModelMetadata _meta;

public:
    error_model model;
    unsigned MaxMotifLength, MaxTractLength;
};


IndelErrorModel generate_new_indel_error_model();

IndelErrorModel generate_old_indel_error_model();

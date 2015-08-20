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
 * SerializedModel.hh
 *
 *  Created on: Jun 23, 2015
 *      Author: Morten Kallberg
 */

#include "json/json.h"
#include <cassert>
#include <map>

#ifndef C___LIB_CALIBRATION_SERIALIZEDMODEL_HH_
#define C___LIB_CALIBRATION_SERIALIZEDMODEL_HH_
#pragma GCC diagnostic ignored "-Wunused-function"


namespace VARIATION_NODE_TYPE
{
enum index_t
{
    SNP,
    INDEL,
    SIZE
};

static
const std::string
get_label(const index_t i)
{
    switch (i)
    {
    case SNP:
        return "SNP";
    case INDEL:
        return "INDEL";
    default:
        assert(false && "Unknown variation node type in serializedModel.");
        return nullptr;
    }
}
}


class serialized_model
{
public:
    serialized_model() {}

    /** methods for serializing */
//    virtual void Serialize( Json::Value& root );
    void Deserialize( const Json::Value& root);
    std::string get_model_string();
    std::string Clean_string(const std::string& str) const;

protected:
    std::string name;
    std::string version;
    std::string date;
};


typedef std::map<int, double> feature_type;
class serialized_calibration_model : public serialized_model
{
public:
    serialized_calibration_model() {}

    /** methods for serializing */
    void Deserialize( const Json::Value& root);
    double getProb(const feature_type& features) const;
    bool doFilter() const;

protected:
    std::string ModelType;
    std::string Type;
    double FilterCutoff;
    // add feature sequence
};

#endif /* C___LIB_CALIBRATION_SERIALIZEDMODEL_HH_ */

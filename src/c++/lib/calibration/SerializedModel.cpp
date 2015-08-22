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
 * SerializedModel.cpp
 *
 *  Created on: Jun 23, 2015
 *      Author: mkallberg
 */

#include "calibration/SerializedModel.hh"

#include <cassert>

namespace SMODEL_ENTRY_TYPE
{
enum index_t
{
    NAME,
    VERSION,
    DATE,
    MODELTYPE,
    TYPE,		//SNV,INDEL ect.
    FILTERCUTOFF,
    SIZE
};

static
const char*
get_label(const index_t i)
{
    switch (i)
    {
    case NAME:
        return "Name";
    case VERSION:
        return "Version";
    case DATE:
        return "Date";
    case MODELTYPE:
        return "ModelType";
    case TYPE:
        return "Type";
    case FILTERCUTOFF:
        return "FilterCutoff";
    default:
        assert(false && "Unknown serialized calibration model entry type");
        return nullptr;
    }
}
}



/// remove newline and quotes from string
static
std::string
Clean_string(const std::string& str)
{
    std::string temp = str;
    // strip any trailing newline
    std::size_t lastNewline = temp.find_last_of("\n\r");
    temp = temp.substr(0, lastNewline);

    // strip flanking double quotes if present
    std::size_t doubleQuote = temp.find_last_of("\"");
    temp = temp.substr(0, doubleQuote);

    doubleQuote = temp.find_first_of("\"");
    if (doubleQuote != std::string::npos)
    {
        temp = temp.substr(doubleQuote);
    }

    return std::move(temp);
}



void serialized_model::Deserialize(const Json::Value& root)
{
    using namespace SMODEL_ENTRY_TYPE;
    this->name 		= Clean_string(root[get_label(NAME)].asString());
    this->version 	= Clean_string(root[get_label(VERSION)].asString());
    this->date 		= Clean_string(root[get_label(DATE)].asString());
}

void serialized_calibration_model::Deserialize( const Json::Value& root)
{
    serialized_model::Deserialize(root);

    //switch these to enums
    using namespace SMODEL_ENTRY_TYPE;
    this->ModelType 		= root[get_label(MODELTYPE)].asString();
    this->Type 				= root[get_label(TYPE)].asString();
    this->FilterCutoff		= root[get_label(FILTERCUTOFF)].asDouble();

    //TODO load feature list here


}


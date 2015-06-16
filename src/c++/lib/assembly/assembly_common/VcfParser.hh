// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Rumovsky
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/// \file

#pragma once

#include <assembly_common/VcfRecord.hh>
#include <iostream>
#include <sstream>
#include <map>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>


//struct VcfRecord; // forward declaration

//typedef std::vector<VcfRecord>             vcfVectorT;  // vector of VCF entries
//typedef std::map<int,VcfRecord>            vcfMapT;     // hashes VCF entries by position
//typedef std::map<std::string, std::string> vcfInfoMapT; // stores key,value pairs from the VCF info field.


class VcfParser {
public:
    VcfParser() {};

    bool buildVcfRecordFromString(const std::string& vcfLine, VcfRecord& rec);
    void _tokenizeLine(const std::string& vcfLine, std::vector<std::string>& vcfTokens);
    bool _convertPosAndQual(const std::string& posStr, const std::string& qualStr, VcfRecord& rec);
};



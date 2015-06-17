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

#include <iostream>
#include <sstream>
#include <map>

//#include <boost/lexical_cast.hpp>
//#include <boost/tokenizer.hpp>
//#include <boost/foreach.hpp>
//#include <boost/algorithm/string.hpp>


//struct VcfRecord; // forward declaration

//typedef std::vector<VcfRecord>             vcfVectorT;  // vector of VCF entries
//typedef std::map<int,VcfRecord>            vcfMapT;     // hashes VCF entries by position
typedef std::map<std::string, std::string> vcfInfoMapT; // stores key,value pairs from the VCF info field.

struct VcfRecord {

    VcfRecord() : pos(0),len(0),qual(0) {}

    VcfRecord(std::string c, int p, int l, std::string i, std::string r, std::string a,
              int q, std::string fi, std::string in, std::string fo, std::string g);

    std::string chr;
    int pos;
    int len; // length can be negative for backward loops in a breakpoint
    std::string id;
    std::string ref;
    std::string alt;
    float qual;
    std::string filter;
    std::string format;
    std::string gt;
    vcfInfoMapT infoFields;

    void clear();
    bool isDeletion() const;
    bool isOpenEnded() const;
    bool isLeftOpenEnded() const;
    bool isRightOpenEnded() const;
    bool isTranslocation() const;
    void updateVcfFilter(const std::string& newFilterTag);
    void updateVcfGT (const std::string& , const std::string& );
    bool parseInfoStringAndAssignLength(const std::string& infoString, const bool reverse);
};

struct SortVcfRecordByPos {
    bool operator()(const VcfRecord& r1, const VcfRecord& r2) {
        return (r1.pos < r2.pos);
    }
};

std::ostream& operator<<(std::ostream& os, const VcfRecord& r);


// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2013 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/downloads/sequencing/licenses/>.
//

/// \file
///
/// \author Chris Saunders
///

#ifndef __BLT_ARG_PARSE_UTIL_HH
#define __BLT_ARG_PARSE_UTIL_HH

#include "blt_util/blt_types.hh"
#include "blt_util/prog_info.hh"

#include <boost/lexical_cast.hpp>

#include <sstream>
#include <string>
#include <vector>



struct arg_data {

    arg_data(int init_argc,
             char* init_argv[],
             const prog_info& init_pinfo);

    // This ctor is designed to continue parsing after another program
    // has already parsed a subset of the arguments. In this case the
    // program name won't be in $0 and cmdline needs to be provided.
    //
    arg_data(const std::vector<std::string>& arg,
             const prog_info& init_pinfo,
             const std::string& init_cmdline);

    void finalize_args();

    unsigned
    size() { return argstr.size(); }

    std::vector<bool> argmark;
    std::vector<std::string> argstr;
    std::string cmdline;
    const prog_info& pinfo;
};


template <typename T>
void
set_val(const prog_info& pinfo,
        const char* arg_label,
        const char* arg,
        T& val) {
    try {
        val=boost::lexical_cast<T>(arg);
    } catch(boost::bad_lexical_cast& e) {
        std::ostringstream oss;
        oss << "argument after flag " << arg_label << " (" << arg << ") cannot be parsed to expected type: " << e.what();
        pinfo.usage(oss.str().c_str());
    }
}


template <typename T>
void
set_arg(unsigned& argi,
        arg_data& ad,
        bool& is_val_set,
        T& val) {

    const char* arg_label(ad.argstr[argi].c_str());

    ad.argmark[argi] = true;
    if(++argi>=ad.argstr.size()) {
        ad.pinfo.usage((std::string("no value following ")+arg_label).c_str());
    }
    if(is_val_set) { ad.pinfo.usage((std::string("multiple ")+arg_label+" arguments").c_str()); }
    set_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val);
    is_val_set=true;
}


void
set_xrange_arg(unsigned& argi,
               arg_data& ad,
               bool& is_val_set,
               double& val,
               bool is_allow_zero=false,
               bool is_no_max_check=false);

void
set_filename_arg(unsigned& argi,
                 arg_data& ad,
                 bool& is_val_set,
                 std::string& file);

void
set_win_arg(unsigned& argi,
            arg_data& ad,
            bool& is_val_set,
            int& val1,
            unsigned& val2);

void
set_nploid_arg(unsigned& argi,
               arg_data& ad,
               bool& is_val_set,
               int& val1,
               double& val2);

void
set_xrange_win_arg(unsigned& argi,
                   arg_data& ad,
                   bool& is_val_set,
                   double& val1,
                   unsigned& val2);

#endif

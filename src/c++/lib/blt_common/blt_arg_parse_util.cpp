//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
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

/// \file
///
/// \author Chris Saunders
///

#include "blt_common/blt_arg_parse_util.hh"



arg_data::
arg_data(int init_argc,
         char* init_argv[],
         const prog_info& init_pinfo)
    : argmark(init_argc,false)
    , argstr(init_argc)
    , pinfo(init_pinfo)
{

    for (int i(0); i<init_argc; ++i) argstr[i]=init_argv[i];

    if (argmark.size()) argmark[0]=true;
    for (unsigned i(0); i<argstr.size(); ++i)
    {
        if (i) cmdline += ' ';
        cmdline += argstr[i];
    }
}

arg_data::
arg_data(const std::vector<std::string>& arg,
         const prog_info& init_pinfo,
         const std::string& init_cmdline)
    : argmark(arg.size(),false)
    , argstr(arg)
    , cmdline(init_cmdline)
    , pinfo(init_pinfo) {}


void
arg_data::
finalize_args()
{
    const unsigned as(size());
    for (unsigned i(0); i<as; ++i)
    {
        if (! argmark[i]) pinfo.usage((std::string("Invalid argument: ")+argstr[i]).c_str());
    }



}



static
bool
is_valid_xrange(const double val,
                bool is_allow_zero,
                bool is_no_max_check)
{
    return (((val > 0.) || (is_allow_zero && (val >= 0.))) && (is_no_max_check || (val <= 1.)));
}



static
void
set_xrange_val(const prog_info& pinfo,
               const char* arg_label,
               const char* arg,
               double& val,
               bool is_allow_zero,
               bool is_no_max_check)
{
    bool is_val_invalid(false);
    try
    {
        val=boost::lexical_cast<double>(arg);
    }
    catch (boost::bad_lexical_cast&)
    {
        is_val_invalid=true;
    }
    is_val_invalid = (is_val_invalid || (! is_valid_xrange(val,is_allow_zero,is_no_max_check)));
    if (is_val_invalid)
    {
        std::ostringstream oss;
        oss << "argument after flag " << arg_label << " (" << arg << ") is not valid";
        pinfo.usage(oss.str().c_str());
    }
}



void
set_xrange_arg(unsigned& argi,
               arg_data& ad,
               bool& is_val_set,
               double& val,
               bool is_allow_zero,
               bool is_no_max_check)
{

    const char* arg_label(ad.argstr[argi].c_str());

    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("no value following ")+arg_label).c_str());
    }
    if (is_val_set)
    {
        ad.pinfo.usage((std::string("multiple ")+arg_label+" arguments").c_str());
    }
    set_xrange_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val,is_allow_zero,is_no_max_check);
    is_val_set=true;
}



void
set_filename_arg(unsigned& argi,
                 arg_data& ad,
                 bool& is_val_set,
                 std::string& file)
{

    const char* arg_label(ad.argstr[argi].c_str());

    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("no value following: ")+arg_label).c_str());
    }
    if (is_val_set)
    {
        ad.pinfo.usage((std::string("multiple ")+arg_label+" arguments").c_str());
    }

    file=ad.argstr[argi];

    if (file.empty())
    {
        ad.pinfo.usage((std::string("empty filename following: ")+arg_label).c_str());
    }

    is_val_set=true;
}



void
set_win_arg(unsigned& argi,
            arg_data& ad,
            bool& is_val_set,
            int& val1,
            unsigned& val2)
{

    const char* arg_label(ad.argstr[argi].c_str());

    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("two arguments must follow ")+arg_label).c_str());
    }
    if (is_val_set)
    {
        ad.pinfo.usage((std::string("multiple ")+arg_label+" arguments").c_str());
    }

    set_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val1);
    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("two arguments must follow ")+arg_label).c_str());
    }

    int val2_tmp;
    set_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val2_tmp);
    if (val2_tmp<0)
    {
        ad.pinfo.usage((std::string("second argument following: ")+arg_label+" must be non-negative\n").c_str());
    }
    val2=val2_tmp;

    is_val_set=true;
}



void
set_nploid_arg(unsigned& argi,
               arg_data& ad,
               bool& is_val_set,
               int& val1,
               double& val2)
{

    const char* arg_label(ad.argstr[argi].c_str());

    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("two arguments must follow ")+arg_label).c_str());
    }
    if (is_val_set)
    {
        ad.pinfo.usage((std::string("multiple ")+arg_label+" arguments").c_str());
    }

    set_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val1);
    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("two arguments must follow ")+arg_label).c_str());
    }

    const bool is_allow_zero(false);
    const bool is_no_max_check(false);
    set_xrange_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val2,is_allow_zero,is_no_max_check);
    is_val_set=true;
}



void
set_xrange_win_arg(unsigned& argi,
                   arg_data& ad,
                   bool& is_val_set,
                   double& val1,
                   unsigned& val2)
{

    const char* arg_label(ad.argstr[argi].c_str());

    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("two arguments must follow ")+arg_label).c_str());
    }
    if (is_val_set)
    {
        ad.pinfo.usage((std::string("multiple ")+arg_label+" arguments").c_str());
    }

    const bool is_allow_zero(false);
    const bool is_no_max_check(false);
    set_xrange_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val1,is_allow_zero,is_no_max_check);

    ad.argmark[argi] = true;
    if (++argi>=ad.argstr.size())
    {
        ad.pinfo.usage((std::string("two arguments must follow ")+arg_label).c_str());
    }

    int val2_tmp;
    set_val(ad.pinfo,arg_label,ad.argstr[argi].c_str(),val2_tmp);
    if (val2_tmp<0)
    {
        ad.pinfo.usage((std::string("second argument following: ")+arg_label+" must be non-negative\n").c_str());
    }
    val2=val2_tmp;

    is_val_set=true;
}

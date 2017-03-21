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

///
/// \author Chris Saunders
///

#include "blt_common/blt_arg_validate.hh"
#include "blt_util/log.hh"

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <sstream>


const double MAX_DIPLOID_THETA(0.38068); // solution of: 0 = 1/3*(theta**2) + (5/2)*theta - 1



void
check_option_arg_range(const double val,
                       const char* label,
                       const double min,
                       const double max,
                       std::string& errorMsg)
{
    errorMsg.clear();
    if ((val >= min) && (val <= max)) return;

    std::ostringstream oss;
    oss << std::setprecision(10);
    oss << "Value provided for '" << label << "': '" << val
        << "', is not in expected range: [ " << min << " , " << max << " ]";
    errorMsg = oss.str();
}



void
check_option_arg_range(const prog_info& pinfo,
                       const double val,
                       const char* label,
                       const double min,
                       const double max)
{
    std::string errorMsg;
    check_option_arg_range(val,label,min,max,errorMsg);
    if (errorMsg.empty()) return;
    pinfo.usage(errorMsg.c_str());
}



void
validate_blt_opt(
    const prog_info& pinfo,
    const blt_options& opt)
{
    if (opt.is_max_win_mismatch && opt.max_win_mismatch_flank_size > MAX_FLANK_SIZE)
    {
        std::ostringstream oss;
        oss << "max-window-mismatch flank size exceeds max value of: " << MAX_FLANK_SIZE;
        pinfo.usage(oss.str().c_str());
    }

    if (opt.bsnp_diploid_theta>MAX_DIPLOID_THETA)
    {
        std::ostringstream oss;
        oss << "diploid heterozygosity exceeds maximum value of: " << MAX_DIPLOID_THETA;
        pinfo.usage(oss.str().c_str());
    }

    if (opt.is_bsnp_nploid && opt.bsnp_nploid_ploidy <= 0)
    {
        pinfo.usage("ERROR:: ploidy argument to -bsnp-nploid must be greater than 0\n");
    }

    if (opt.is_bsnp_diploid_het_bias)
    {
        if (! opt.is_bsnp_diploid())
        {
            pinfo.usage("-bsnp-diploid-het-bias does not make sense when bsnp-diploid model is not in use\n");
        }
        if ((opt.bsnp_diploid_het_bias < 0.) || (opt.bsnp_diploid_het_bias >= 0.5))
        {
            pinfo.usage("-bsnp-diploid-het-bias argument must be in range [0,0.5)\n");
        }
    }
}

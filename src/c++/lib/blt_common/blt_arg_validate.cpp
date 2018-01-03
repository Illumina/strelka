//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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

#include "blt_common/blt_arg_validate.hh"

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

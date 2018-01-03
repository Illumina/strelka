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

#pragma once


#include "blt_common/blt_shared.hh"
#include "blt_util/prog_info.hh"


extern const double MAX_DIPLOID_THETA;


void
check_option_arg_range(const double val,
                       const char* label,
                       const double min,
                       const double max,
                       std::string& errorMsg);

void
check_option_arg_range(const prog_info& pinfo,
                       const double val,
                       const char* label,
                       const double min,
                       const double max);

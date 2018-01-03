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

///
/// \author Chris Saunders
///

#pragma once

#include "blt_util/seq_util.hh"

#include <boost/utility.hpp>


struct nonref_test_call : private boost::noncopyable
{
    nonref_test_call()
        : is_snp(false),
          snp_qphred(0),
          max_gt_qphred(0),
          max_gt(0),
          nonref_id(BASE_ID::ANY) {}

    bool is_snp;
    int snp_qphred;
    int max_gt_qphred;
    unsigned  max_gt;
    unsigned nonref_id;
};


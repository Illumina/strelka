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

#include "starling_shared.hh"

#include "blt_common/snp_pos_info.hh"

#include <iosfwd>
#include <memory>


/// a simple pileup caller used in simulations
///
struct starling_pile_caller
{
    starling_pile_caller(
        starling_options& opt,
        std::ostream& os);

    void
    call(
        const unsigned pos,
        const snp_pos_info& pi);

private:
    starling_options& _opt;
    std::unique_ptr<starling_deriv_options> _dopt_ptr;
    std::ostream& _os;
};


void
starling_pile_test_run(
    starling_options& opt);

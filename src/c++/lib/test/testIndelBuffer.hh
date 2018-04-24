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

#include "starling_base_options_test.hh"
#include "blt_util/reference_contig_segment.hh"
#include "starling_common/IndelBuffer.hh"

#include <memory>


/// Manages a mock indel buffer approriate for unit testing
///
// TODO Refactor this to a factory function returning a unique_ptr to an IndelBuffer instead
struct TestIndelBuffer
{
    explicit
    TestIndelBuffer(
        const reference_contig_segment& ref);

    IndelBuffer&
    getIndelBuffer()
    {
        return *_IndelBufferPtr;
    }

private:
    starling_base_options_test _opt;
    std::unique_ptr<starling_base_deriv_options> _doptPtr;
    std::unique_ptr<IndelBuffer> _IndelBufferPtr;
};

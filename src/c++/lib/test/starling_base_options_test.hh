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

#include "options/AlignmentFileOptions.hh"
#include "starling_common/starling_base_shared.hh"


/// \brief This version of starling_base_options provides null implementations of required virtuals
///        so that unit tests are easier to setup
struct starling_base_options_test final : public starling_base_options
{
    const AlignmentFileOptions&
    getAlignmentFileOptions() const override
    {
        static const AlignmentFileOptions alignFileOpt(getTestAlignmentFileOptions());
        return alignFileOpt;
    }

private:
    static
    AlignmentFileOptions
    getTestAlignmentFileOptions()
    {
        AlignmentFileOptions alignFileOpt;
        alignFileOpt.alignmentFilenames.push_back("sample.bam");
        return alignFileOpt;
    }
};

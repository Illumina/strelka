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

#include "testIndelBuffer.hh"



TestIndelBuffer::
TestIndelBuffer(
    const reference_contig_segment& ref)
{
    _opt.is_candidate_indel_signal_test = false;

    const double maxDepth = 100.0;
    _doptPtr.reset(new starling_base_deriv_options(_opt));

    _IndelBufferPtr.reset(new IndelBuffer(_opt, *_doptPtr, ref));

    _IndelBufferPtr->registerSample(depth_buffer(), depth_buffer(), maxDepth);
    _IndelBufferPtr->finalizeSamples();
}

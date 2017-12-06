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

/// \file
/// \author Chris Saunders
///

#pragma once

#include <cassert>

#include <limits>


/// Compress integer input so that higher resolution is preserved for
/// values near zero, while reducing systematic bias between input and
/// compressed output
///
/// The motivating use case is to reduce the total number of integers used for
/// a map key, where the integer represents count data, such that high resolution
/// should be preserved near zero
///
/// The compression scheme is: The highest "bitCount" bits from the input integer
/// are preserved, followed by a suffix. The suffix is fixed for all input with the
/// same prefix, and is set in such a way that systematic bias between the input
/// and compressed output is reduced.
///
/// This function deterministically switches between two suffix values, one
/// has only the high bit set, the other is the complement of the first. The suffix
/// is chosen based on the low bit of the prefix (or shift value if bitCount is one).
///
/// Example: For bitCount=3 and input=67:
///
/// input  is 67 or 0b1000011
/// output is 72 or 0b1001000
///
/// Here the top three bits are preserved "0b100XXXX", and the suffix "0bXXX1000" is
/// replaced on the lower bits to create the compressed output.
///
/// T must be of an unsigned integral type
///
template <typename T>
T
compressInt(
    const T input,
    const unsigned bitCount)
{
    using Tinfo = std::numeric_limits<T>;
    static_assert((! Tinfo::is_signed) && Tinfo::is_integer, "T must have unsigned integral type");

    assert(bitCount>0);

    // find last bit (should match POSIX fls() function)
    unsigned highBitIndex(0);
    {
        T input2(input);
        while (input2 >> highBitIndex)
        {
            highBitIndex++;
        }
    }

    if (highBitIndex <= bitCount) return input;

    const unsigned shift(highBitIndex - bitCount);
    const T prefix(input >> shift);

    // switch off between two different suffix schemes to reduce bias
    // scheme 1: suffix is 0b10000...
    // scheme 2: suffix is 0b01111...
    static const T one (1);
    T suffix(one << (shift - 1));
    if ((bitCount == 1) ? (shift & 0x1) : (prefix & 0x1))
    {
        suffix -= one;
    }
    return ((prefix << shift) | suffix);
}


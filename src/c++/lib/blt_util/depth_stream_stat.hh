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

#include <cmath>

#include <iosfwd>
#include <limits>


/// \brief Simple on-line statistics for unsigned values
///
/// derived From Tony Cox's IndelFinder code
///
struct depth_stream_stat
{
    // Accumulate mean and standard dev using a single pass formula
    // Uses the cancellation-friendly formulae on p.26 of
    // Higham, Accuracy & Stability of Numerical Algorithms
    // Variable names follow his
    depth_stream_stat() : M_(0),Q_(0),max_(0),min_(0),k_(0),n_(0) {}

    void update (const unsigned d)
    {
        k_++;
        if (d!=0) n_++;
        if ((k_==1) || (d>max_)) max_=d;
        if ((k_==1) || (d<min_)) min_=d;

        // important to do M before Q as Q uses previous iterate of M
        const double delta(static_cast<double>(d)-M_);
        M_+=delta/static_cast<double>(k_);
        Q_+=delta*(static_cast<double>(d)-M_);
    }

    unsigned sample_size() const
    {
        return k_;
    }
    unsigned nonzero() const
    {
        return n_;
    }
    double min() const
    {
        return ((k_<1) ? std::numeric_limits<double>::quiet_NaN() : min_);
    }
    double max() const
    {
        return ((k_<1) ? std::numeric_limits<double>::quiet_NaN() : max_);
    }
    double mean() const
    {
        return ((k_<1) ? std::numeric_limits<double>::quiet_NaN() : M_);
    }
    double variance() const
    {
        return ((k_<2) ? std::numeric_limits<double>::quiet_NaN() : Q_/(static_cast<double>(k_-1)));
    }
    double sd() const
    {
        return std::sqrt(variance());
    }
    double stderror() const
    {
        return sd()/std::sqrt(static_cast<double>(k_));
    }

private:

    double M_;
    double Q_;
    unsigned max_;
    unsigned min_;
    unsigned k_;
    unsigned n_;
};


std::ostream& operator<<(std::ostream& os,const depth_stream_stat& ss);


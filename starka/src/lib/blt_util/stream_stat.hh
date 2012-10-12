// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

/// \file

/// \author Chris Saunders
///
#ifndef __STREAM_STAT_HH
#define __STREAM_STAT_HH

#include <cmath>

#include <iosfwd>



/// \brief from the indel finder's single_pass mean/sd calculator:
///
struct stream_stat {

    // Accumulate mean and standard dev using a single pass formula
    // Uses the cancellation-friendly formulae on p.26 of
    // Higham, Accuracy & Stability of Numerical Algorithms
    // Variable names follow his
    stream_stat() : M_(0),Q_(0),max_(0),min_(0),k_(0) {}

    void reset() {
        M_ = 0;
        Q_ = 0;
        max_ = 0;
        min_ = 0;
        k_ = 0;
    }

    void add(const double x) {
        k_++;
        if(k_==1 || x>max_) max_=x;
        if(k_==1 || x<min_) min_=x;

        // important to do M before Q as Q uses previous iterate of M
        const double delta(x-M_);
        M_+=delta/static_cast<double>(k_);
        Q_+=delta*(x-M_);
    }

    int size() const { return k_; }
    bool empty() const { return (k_==0); } 

    double min() const { return ((k_<1) ? nan() : min_); }
    double max() const { return ((k_<1) ? nan() : max_); }
    double mean() const { return ((k_<1) ? nan() : M_); }
    double variance() const { return ((k_<2) ? nan() : Q_/(static_cast<double>(k_-1))); }
    double sd() const { return std::sqrt(variance()); }
    double stderror() const { return sd()/std::sqrt(static_cast<double>(k_)); }

private:
    static
    double nan() { double a(0.); return 0./a; }  // 'a' supresses compiler warnings

    double M_;
    double Q_;
    double max_;
    double min_;
    unsigned k_;
};


std::ostream& operator<<(std::ostream& os,const stream_stat& ss);


#endif

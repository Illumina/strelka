// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2009-2012 Illumina, Inc.
//
// This software is covered by the "Illumina Genome Analyzer Software
// License Agreement" and the "Illumina Source Code License Agreement",
// and certain third party copyright/licenses, and any user of this
// source file is bound by the terms therein (see accompanying files
// Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
// Illumina_Source_Code_License_Agreement.pdf and third party
// copyright/license notices).
//
//

/// \file

/// \author Chris Saunders
///
#ifndef __STREAM_STAT_HH
#define __STREAM_STAT_HH

#include <ciso646>
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

    void update (const double x) {
        k_++;
        if(k_==1 || x>max_) max_=x;
        if(k_==1 || x<min_) min_=x;

        // important to do M before Q as Q uses previous iterate of M
        const double delta(x-M_);
        M_+=delta/static_cast<double>(k_);
        Q_+=delta*(x-M_);
    }

    int sample_size() const { return k_; }
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

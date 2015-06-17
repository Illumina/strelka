// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Rumovsky
// Copyright (c) 2013-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//

///
/// \author Ole Schulz-Trieglaff
///

#pragma once

//#include <iosfwd>
#include <string>
#include <vector>


/// \brief data pertaining to a de-novo assembly contig
///
/// stores for each contig the sequence and the number of reads
/// containing its seeding k-mer
///
struct Contig
{

    Contig() :
        seq(""),
        avgCoverage(0), // avg kmer depth
        numKmers(0),	// total kmer count in contig
        seedReadCount(0), //
        readSupport(0), //
        hasSinkKmer(false) // contig contains the final ref kmer
    {}

    Contig(const Contig& c) :
        seq(c.seq),
        avgCoverage(c.avgCoverage),
        numKmers(c.numKmers),
        seedReadCount(c.seedReadCount),
        readSupport(c.readSupport),
        hasSinkKmer(c.hasSinkKmer)
    {}

    void clear()
    {
        seq="";
        avgCoverage = 0.0;
        numKmers = 0;
        seedReadCount = 0;
        readSupport = 0;
        hasSinkKmer = false;
    }

    std::string seq;
    float avgCoverage;
    int numKmers;
    unsigned seedReadCount;
    unsigned readSupport;
    bool hasSinkKmer;
};

typedef std::vector<Contig> Assembly;

std::ostream& operator<<(std::ostream& os, const Contig& contig);


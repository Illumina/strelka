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


struct GreedyAssemblerOptions
{
    // sets reasonable default values for 30x DNA-seq, 100bp reads
    GreedyAssemblerOptions() :
        minWordLength(41),
        maxWordLength(76),
        wordStepSize(5),
        minContigLength(15),
        kmerPruningThreshold(1),
        minSeedReads(3),
        maxAssemblyIterations(3),
        doPruning(false),
        verbose(false)
    {}

    std::string outFullGraphPrefix; // write the full graph to this file
    std::string outGraphPrefix; // write the simplified graph to this file
    unsigned minWordLength; // initial word (kmer) length
    unsigned maxWordLength; // max word length
    unsigned wordStepSize;
    unsigned minContigLength; // min contig size
    unsigned kmerPruningThreshold; // cutoff for graph pruning
    unsigned minSeedReads; // min. number of reads required to start assembly
    unsigned maxAssemblyIterations; // Max. number of assembly iterations for a cluster before we give up
    bool doPruning;
    bool verbose;
};


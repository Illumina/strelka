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

#include <assembly/assembly_common/AssemblyReadInfo.hh>
#include <assembly/assembly_common/BamAddOns.hh>
#include <assembly/assembly_common/Contig.hh>
#include "assembly/kmer/DfsAssemblerOptions.hh"

#include "assembly/graph/GraphBase.hh"

#include "boost/foreach.hpp"
#include "boost/unordered_map.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/utility.hpp"

#include <cassert>
#include <vector>
#include <fstream>
#include <iostream>

//#define DEBUG_DB_GRAPH

namespace rumovsky {

///
/// implementation of a de Bruijn graph
class DeBruijnGraph {

public:
    DeBruijnGraph() : wordLength_(0) {}

    DeBruijnGraph(const unsigned wl) : wordLength_(wl) {}

    unsigned getWordLength() const { return wordLength_; }
    void setWordLength(const unsigned wl) { wordLength_ = wl; }

    void setRefKmer(const std::string& vertex) { kmerRefTable_[vertex] = true; }
    bool isRefKmer(const std::string& vertex) const { return (kmerRefTable_.find(vertex) == kmerRefTable_.end()); }

    ///
    /// sets all reference k-mers from reference sequence
    void setRefKmers(const std::string& ref);

    ///
    /// Removes k-mers below min freq
    void prune(const unsigned minFreq);

    ///
    /// build de-bruijn graph from a set of reads
    bool build(const AssemblyReadInput& reads, const AssemblyReadOutput& readInfo, const bool addreverse);



    ///
    /// Writes doxygen and json output
    void dumpToFile(const std::string& outFilePrefix, const std::string& refSeq, const std::string& src) const;

    ///
    /// Gives all vertices connected with \param vertex
    void getNeighbours (const std::string& vertex, AdjList& adj) const;

    ///
    /// Tests if \param vertex is contained in graph
    bool isInGraph(const std::string& vertex) const { return (kmerFreqTable_.find(vertex) != kmerFreqTable_.end()); }

    ///
    /// Returns the number of vertices
    size_t size() const { return kmerFreqTable_.size(); }

    const str_vec_map_t&  getAdjList() const { return adjList_;}

    const str_uint_map_t&  getKmerFreqTable() const { return kmerFreqTable_;}

    void contractPaths(const std::string& startVertex);


private:

    void buildAdjList_();


    std::string doPathContract_(Contig& contigSoFar,
                                str_bool_map_t& seenVertices,
                                str_uint_map_t& newKmerFreqTable,
                                str_vec_map_t& newAdjList);

    str_uint_map_t kmerFreqTable_;
    str_uint_vec_map_t kmerReadPosTable_;
    str_vec_map_t adjList_;
    std::vector<str_uint_map_t> readWordOffsets_;
    str_bool_map_t kmerRefTable_;

    unsigned wordLength_;
};

} // end of namespace rumovsky


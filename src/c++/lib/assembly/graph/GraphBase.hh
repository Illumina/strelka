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
#include <string>

enum AFFIX_TYPE { PREFIX, SUFFIX };
// keeps record of k-mer frequencies
typedef boost::unordered_map<std::string,unsigned> str_uint_map_t;
// maps kmers to positions in read
typedef boost::unordered_map<std::string,std::vector<unsigned> > str_uint_vec_map_t;
//
typedef std::pair<std::string,AFFIX_TYPE> Neighbour;
// adjacency list for k-mers
typedef std::vector<Neighbour > AdjList;
//
typedef boost::unordered_map<std::string, AdjList> str_vec_map_t;
// keeps track of visited k-mers
typedef boost::unordered_map<std::string,bool> str_bool_map_t;


namespace rumovsky {

///
/// abstract base class of an assembly graph
class GraphBase {
public:
    GraphBase() {}
    virtual ~GraphBase() {}

    ///
    /// Build assembly graph from reads
    /// \return true if successful
    virtual bool build(const AssemblyReadInput& reads, const AssemblyReadOutput& readInfo, const bool addReverse=false) = 0;

    virtual void dumpToFile(const std::string& outFilePrefix, const std::string& refSeq, const std::string& src) const = 0;

    virtual void getNeighbours(const std::string& v, AdjList& adj) const = 0;

    //virtual bool isInGraph(const std::string& kmer) const = 0;
};

} // end of namespace

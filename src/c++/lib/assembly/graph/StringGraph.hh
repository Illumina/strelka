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

#include "GraphBase.hh"

#include <string>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>

#include <boost/lexical_cast.hpp>

//#define DEBUG_WORDS
//#define DEBUG_OVERLAPS

namespace rumovsky {

enum NodeState
{
    Vacant=-1000,
    //  InPlay=1,
    Eliminated=-2000
};

///
/// Stores the location of a kmer and the number of the containing read
struct WordLocation
{
    WordLocation(unsigned read, unsigned pos) : read_(read), pos_(pos) {}
    unsigned read_;
    unsigned pos_;
    bool operator<(const WordLocation& rhs) const
    {
        return ((read_<rhs.read_)||
                ((read_==rhs.read_)&&(pos_<rhs.pos_)));
    } // ~operator <
}; // ~struct WordLocation

struct Edge
{
    Edge( unsigned node, std::string label, int weight=1) :
        label_(label), weight_(weight), node_(node), visited_(false)
    {
#ifdef DEBUG_EDGE
        cout << "constructing edge " << this << endl;
        printf("node=%d label=%s weight=%d\n",node_,label_.c_str(),weight);
#endif
    }

    Edge( const Edge& rhs ):
        label_(rhs.label_),
        weight_(rhs.weight_),
        node_(rhs.node_),
        visited_(rhs.visited_)
    {
#ifdef DEBUG_EDGE
        cout << "copy constructing edge " << this << endl;
#endif
    }

    ~Edge()
    {
#ifdef DEBUG_EDGE
        cout << "destructing edge " << this << endl;
#endif
    }

    //  state_(Vacant) {}
    std::string label_;
    int weight_;
    //  EdgeState state_;
    unsigned node_;
    bool visited_;
    bool operator<(const Edge& rhs) const
    {
        return (label_.size()<rhs.label_.size());
    } // ~operator <
}; // ~struct WordLocation};

struct Node
{
    Node() : state_(Vacant) {}
    std::list<Edge> edges_;
    int state_;
    //  NodeState state_;

    bool isInPlay() const { return state_>0; }
};

struct Walk
{
    Walk() : isNeeded_(true), start_(0), end_(0) {}
    std::string label_;
    bool isNeeded_;
    int start_;
    int end_;
    void print( void )
    {
        printf("%d %s %d\n",start_,label_.c_str(),end_);
    } // ~print
}; // ~struct Walk

struct WalkEnd
{
    WalkEnd( bool isEnd, int walkNum ): isEnd_(isEnd), walkNum_(walkNum) {}

    bool isEnd_; // true then end of walk else start
    int walkNum_;
    bool operator<(const WalkEnd& rhs) const
    {
        return ((isEnd_==true)&&(rhs.isEnd_==false));
    }
};

typedef std::map<std::string, std::vector<WordLocation> > WordStoreT;
typedef std::map< int, std::map<int, int > > HitStoreT; // idea is hits[read][offset]
typedef std::map<std::string, std::vector<int> > ReadStoreT;
typedef std::vector<Node> NodeStoreT;
typedef std::vector<unsigned> PathT;
typedef std::map< WordLocation, std::map< char, int > > AllMatchToN;

struct MatchToN
{
    MatchToN(int read, int pos, char base) :
        read_(read), pos_(pos), base_(base) {}
    int read_;
    int pos_;
    char base_;
}; // ~struct WordLocation


class StringGraph: public GraphBase {
public:
    StringGraph() : wordLength_(0), overlapLength_(0) {}

    StringGraph(unsigned w, unsigned o) : wordLength_(w), overlapLength_(o) {}

    virtual ~StringGraph() {}
    ///
    /// Build assembly graph from reads
    /// \return true if successful
    bool build(const AssemblyReadInput& reads, const AssemblyReadOutput& readInfo);

    void dumpToFile(const std::string& outFilePrefix) const;

    void getNeighbours(const std::string& v, AdjList& adj) const;

    unsigned size() const { return nodes_.size(); }


private:

    bool containsN_( const std::string& s );
    bool compareLabelsWithNs_( const std::string& label_small, const std::string& label_large);

    void collectWords_();

    bool
    getOverlap_(const int readLeft,
                const int readRight,
                const int offset,
                int& overlap,
                int& diffs,
                AllMatchToN& allMatchesToN);

    void
    computeOverlaps_();

    int
    countHits_(const unsigned readNum);

    void
    correctNs_(const AllMatchToN& allMatchesToN);

    void
    deduceMissingBases_(std::list<Edge>& edges );

    void
    reduceGraph_();

    unsigned wordLength_;
    unsigned overlapLength_;

    ReadStoreT allReads_;
    WordStoreT words_;
    std::vector<std::string> readSeqs_;
    std::vector<bool> isContained_;
    HitStoreT hits_;
    NodeStoreT nodes_;
    //std::vector<Walk> walks_;
    //std::vector<PathT > paths_;
}; // StringGraph

}

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

//#define DEBUG_DFS_WALKER
//#define DEBUG_GREEDY_WALKER

#pragma once

#include "DeBruijnGraph.hh"
#include "NodeChoiceFunction.hh"

#include "common/Contig.hh"

#include <iostream>

namespace rumovsky{

template <class Graph, class DecFun>
class GraphWalkerBase {

public:
    GraphWalkerBase(const Graph& g, const DecFun& f)
    : g_(g), f_(f) {}

    virtual ~GraphWalkerBase() {};

    const Graph& getGraph() const {return g_;}

    virtual void walk(const std::string& sourceKmer,
                      std::vector<Contig>& contigs) = 0;

protected:
    Graph g_;
    DecFun f_;
};

template <class Graph>
class GraphWalkerDFS : public GraphWalkerBase <Graph,NodeChoiceFunctionSimple<Graph> > {

public:
    GraphWalkerDFS(const Graph& g) : GraphWalkerBase<Graph,NodeChoiceFunctionSimple<Graph> >(g, NodeChoiceFunctionSimple<Graph>(g)) { }

    virtual ~GraphWalkerDFS() {};

    void walk(const std::string& sourceKmer,
              std::vector<Contig>& contigs)
    {
        std::string startVertex = "NA";
        str_uint_map_t::const_iterator ct;
        for (ct = this->g_.getKmerFreqTable().begin();
             ct != this->g_.getKmerFreqTable().end();
             ++ct)
        {
            if (ct->first.find(sourceKmer) != std::string::npos)
            {
                startVertex = ct->first;
                break;
            }
        }

        assert(ct!= this->g_.getKmerFreqTable().end());
        //std::cout << "Setting start vertex to " << startVertex << "\n";

        AssemblyContext ctxt;
        ctxt.ctg_.seq           = startVertex;
        ctxt.ctg_.numKmers      = 1;
        ctxt.ctg_.avgCoverage   = ct->second;
        ctxt.ctg_.seedReadCount = ct->second;

        seenVertices_[startVertex] = 1;

        walk_(ctxt,startVertex,contigs);
    }

private:
    str_uint_map_t seenVertices_;

    void walk_(AssemblyContext& context,
               std::string& vertex,
               std::vector<Contig>& contigs)
    {
#ifdef DEBUG_DFS_WALKER
        std::cerr << __PRETTY_FUNCTION__ << context.ctg_.seq << "\n";
        std::cerr << "vertex : " << vertex << "\n";
#endif

        AdjList neighbours;
        this->f_.getNeighbours(vertex,context,neighbours);
        const unsigned numNeighbours(neighbours.size());
        if (numNeighbours == 0)
        {
#ifdef DEBUG_DFS_WALKER
            std::cerr << "Branch end. Have contig " << context.ctg_.seq << "\n";
#endif
            // at the end of a branch, save contig
            contigs.push_back(context.ctg_);
            return;
        }

        bool movedForward(false);
        for (auto s: neighbours)
        {
            // don't walk back
            if (s.second != AFFIX_TYPE::SUFFIX) continue;
            if (!seenVertices_[s.first])
            {
                //std::cerr << "Testing neighbour " << s.first << " " << s.second << " of " << context.ctg_.seq << "\n";
                ++seenVertices_[s.first];
                // Create new context and copy contig sequence
                AssemblyContext newCtxt(context);
                vertex=s.first;
                // cut-off overlap
                std::string newSeq = s.first.substr(this->getGraph().getWordLength()-1);
                newCtxt.ctg_.seq += newSeq;
                str_uint_map_t::const_iterator ct = this->getGraph().getKmerFreqTable().find(vertex);
                assert(ct != this->getGraph().getKmerFreqTable().end());
                // update coverage
                newCtxt.ctg_.avgCoverage += ct->second;
                ++newCtxt.ctg_.numKmers;
#ifdef DEBUG_DFS_WALKER
                std::cerr << "Extending contig " << context.ctg_.seq << " into " << newCtxt.ctg_.seq << "\n";
#endif
                movedForward = true;
                walk_(newCtxt,vertex,contigs);
            }
        }

        if (!movedForward)
        {
#ifdef DEBUG_DFS_WALKER
            std::cerr << "Exhausted all neighbours. Have contig " << context.ctg_.seq << "\n";
#endif
            contigs.push_back(context.ctg_);
        }
    } // end of step()

};

template <class Graph>
class GraphWalkerGreedy: public GraphWalkerBase <Graph,NodeChoiceFunctionGreedy<Graph> > {

public:
    GraphWalkerGreedy(const Graph& g) : GraphWalkerBase<Graph,NodeChoiceFunctionGreedy<Graph> >(g, NodeChoiceFunctionGreedy<Graph>(g)) { }

    virtual ~GraphWalkerGreedy() {};

    void walk(const std::string& sourceKmer,
              std::vector<Contig>& contigs)
    {
        seenVertices_.clear();

        //std::cout << __PRETTY_FUNCTION__ << "\n";
        str_uint_map_t::const_iterator it = this->g_.getKmerFreqTable().find(sourceKmer);
        assert(it != this->g_.getKmerFreqTable().end());

        AssemblyContext ctxt;
        ctxt.ctg_.seq           = sourceKmer;
        ctxt.ctg_.numKmers      = 1;
        ctxt.ctg_.avgCoverage   = it->second;
        ctxt.ctg_.seedReadCount = it->second;

        seenVertices_[sourceKmer] = 1;

        // TODO : this is bad
        std::string startVertex(sourceKmer);

        bool forward(true);
        // contig forward extension
        walk_(ctxt,startVertex,forward);

        forward=false;
        // contig backward extension
        walk_(ctxt,startVertex,forward);

        // store contig
        contigs.push_back(ctxt.ctg_);
    }
private:
    str_uint_map_t seenVertices_;

    void walk_(AssemblyContext& context,
               std::string& vertex,
               const bool forward)
    {
#ifdef DEBUG_GREEDY_WALKER
        std::cerr << __PRETTY_FUNCTION__ << context.ctg_.seq << "\n";
        std::cerr << "vertex : " << vertex << "\n";
        std::cerr << "isForward " << forward << "\n";
#endif

        AdjList neighbours;
        this->f_.getNeighboursDirectional(vertex,context,forward,neighbours);
        if (neighbours.empty())
        {
#ifdef DEBUG_GREEDY_WALKER
            std::cerr << "Branch end. Have contig " << context.ctg_.seq << "\n";
#endif
            // at the end of a branch
            return;
        }

        for (auto s: neighbours)
        {
#ifdef DEBUG_GREEDY_WALKER
            std::cerr << "Testing neighbour " << s.first << " " << s.second << " of " << context.ctg_.seq << "\n";
            std::cerr << "numNeighbours=" << neighbours.size() << "\n";
            std::cerr << seenVertices_[s.first] << "\n";
#endif
            if (!seenVertices_[s.first])
            {
#ifdef DEBUG_GREEDY_WALKER
                std::cerr << "Walking along neighbour " << s.first << " " << s.second << " of " << context.ctg_.seq << "\n";
#endif
                ++seenVertices_[s.first];
                vertex=s.first;

                if (forward)
                {
                    std::string newSeq = s.first.substr(this->getGraph().getWordLength()-1);
                    context.ctg_.seq += newSeq;
                }
                else
                {
                    std::string newSeq = s.first.substr(0,1);
                    context.ctg_.seq = newSeq + context.ctg_.seq;
                }
                str_uint_map_t::const_iterator ct = this->getGraph().getKmerFreqTable().find(vertex);
                assert(ct != this->getGraph().getKmerFreqTable().end());
                // update coverage
                context.ctg_.avgCoverage += ct->second;
                ++context.ctg_.numKmers;

#ifdef DEBUG_GREEDY_WALKER
                std::cerr << "Extending contig into " << context.ctg_.seq << "\n";
#endif

                walk_(context,vertex,forward);
            } else {
#ifdef DEBUG_GREEDY_WALKER
                std::cout << "Seen this kmer (" <<  s.first  << ") already too often " << "\n";
                std::cout << seenVertices_[s.first] << " vs " << s.second << "\n";
#endif
            }
        }
    } // end of step()


};


}

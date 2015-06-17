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

#define DEBUG_CHOICE_FUNC

#pragma once

#include <string>

#include "DeBruijnGraph.hh"

namespace rumovsky {


///
/// Derive from this struct to pass info on current state of
/// assembly to choice function
struct AssemblyContext {
    ///
    /// Contig assembled so far
    Contig ctg_;
};

template <class Graph>
class NodeChoiceFunctionBase {
public:
    NodeChoiceFunctionBase(const Graph& g) : g_(g) {}
    virtual ~NodeChoiceFunctionBase() {};

    ///
    /// Give neighbours of \param vertex, given the current assembly state \param context
    virtual void getNeighbours(const std::string& vertex, const AssemblyContext& context, AdjList& neighbours) const = 0;

protected:
    Graph g_;
};

///
/// implements a choice functions for DFS
template <class Graph>
class NodeChoiceFunctionSimple : public NodeChoiceFunctionBase <Graph> {

    typedef NodeChoiceFunctionBase<Graph> Base;

public:
    NodeChoiceFunctionSimple(const Graph& g) : Base(g) { }

    virtual ~NodeChoiceFunctionSimple() {};

    ///
    /// Retrieve all neighbours
    void getNeighbours(const std::string& vertex, const AssemblyContext& /*context*/, AdjList& ns) const
    {
       Base::g_.getNeighbours(vertex,ns);
    }
};

///
/// Traverse all neighbours above minimum weight
template <class Graph>
class NodeChoiceFunctionWeighted : public NodeChoiceFunctionBase <Graph> {

    typedef NodeChoiceFunctionBase<Graph> Base;

public:
    NodeChoiceFunctionWeighted(const Graph& g, const unsigned minWeight) : Base(g), minWeight_(minWeight) { }

    virtual ~NodeChoiceFunctionWeighted() {};

    ///
    /// Retrieve all neighbours above minimum weight
    void getNeighbours(const std::string& vertex, const AssemblyContext& /*context*/, AdjList& ns) const
    {
       AdjList& tmp = Base::g_.getNeighbours(vertex,ns);
       for (auto s: tmp) {
           str_uint_map_t::const_iterator ct = Base::g_.getKmerFreqTable().find(vertex);
           assert(ct != Base::g_.getKmerFreqTable().end());
           if (ct->second > minWeight_)
           {
               ns.push_back(s);
           }
       }
    }

private:
    unsigned minWeight_;
};

///
/// Traverse all neighbours above minimum weight
template <class Graph>
class NodeChoiceFunctionGreedy : public NodeChoiceFunctionBase <Graph> {

    typedef NodeChoiceFunctionBase<Graph> Base;

public:
    NodeChoiceFunctionGreedy(const Graph& g) : Base(g) { }

    virtual ~NodeChoiceFunctionGreedy() {};

    ///
    /// Retrieve all neighbours above minimum weight
    void getNeighbours(const std::string& vertex, const AssemblyContext& /*context*/, AdjList& ns) const
    {
       unsigned maxWeight(0);
       Neighbour maxNeighbour;
       AdjList tmp;
       Base::g_.getNeighbours(vertex,tmp);
       for (auto s: tmp) {
           str_uint_map_t::const_iterator ct = Base::g_.getKmerFreqTable().find(vertex);
           assert(ct != Base::g_.getKmerFreqTable().end());
           if (ct->second > maxWeight)
           {
               maxNeighbour = s;
               maxWeight    = ct->second;
           }
       }
       ns.push_back(maxNeighbour);
    }

    ///
    /// Retrieve all neighbours above minimum weight
    void getNeighboursDirectional(const std::string& vertex,
                                  const AssemblyContext& /*context*/,
                                  const bool forward,
                                  AdjList& ns) const
    {
        unsigned maxWeight(0);
        Neighbour maxNeighbour;
        bool maxFound(false);
        AdjList tmp;
        Base::g_.getNeighbours(vertex,tmp);
        for (auto s: tmp) {

            if (forward && s.second != AFFIX_TYPE::SUFFIX)
                continue;
            if (!forward && s.second != AFFIX_TYPE::PREFIX)
                continue;

            str_uint_map_t::const_iterator ct = Base::g_.getKmerFreqTable().find(vertex);
            assert(ct != Base::g_.getKmerFreqTable().end());

            if (ct->second > maxWeight)
            {
                maxNeighbour = s;
                maxWeight    = ct->second;
                maxFound     = true;
            }
        }
        if (maxFound) ns.push_back(maxNeighbour);
    }



};

} // end of namespace

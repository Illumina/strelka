// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Starka
// Copyright (c) 2009-2014 Illumina, Inc.
//
// This software is provided under the terms and conditions of the
// Illumina Open Source Software License 1.
//
// You should have received a copy of the Illumina Open Source
// Software License 1 along with this program. If not, see
// <https://github.com/sequencing/licenses/>
//
/*
 *      Author: mkallberg
 */

#include "RandomForestModel.hh"

#include "blt_util/log.hh"
#include "blt_util/parse_util.hh"

#include <cassert>

using boost::property_tree::ptree;
using namespace illumina::blt_util;



namespace DTREE_NODE_TYPE
{
    enum index_t
    {
        TREE,
        VOTE,
        DECISION,
        SIZE
    };

    const char*
    get_label(const index_t i)
    {
        switch(i)
        {
        case TREE: return "tree";
        case VOTE: return "node_votes";
        case DECISION: return "decisions";
        default:
            assert(false && "Unknown node type");
            return nullptr;
        }
    }
}



template <typename L, typename R>
void
RandomForestModel::
parseTreeNode(
    const ptree::value_type& v,
    TreeNode<L,R>& val)
{
    assert(! val.isInit);
    assert(v.second.size() == 2);

    val.isInit = true;
    ptree::const_iterator viter(v.second.begin());
    val.left = viter->second.get_value<L>();
    ++viter;
    val.right = viter->second.get_value<R>();
}



void
RandomForestModel::
load(const ptree& pt)
{
    clear();

    // trees:
    for (const ptree::value_type& tree_pt : pt)
    {
        _forest.emplace_back();
        DecisionTree& dtree(_forest.back());

        using namespace DTREE_NODE_TYPE;

        // node types:
        for (int i(0); i<SIZE;++i)
        {
            const index_t nodeTypeIndex(static_cast<index_t>(i));
            // nodes:
            for (const ptree::value_type& v : tree_pt.second.get_child(get_label(nodeTypeIndex)))
            {
                const unsigned nodeIndex(parse_unsigned_str(v.first));
                if (dtree.data.size() <= nodeIndex)
                {
                    dtree.data.resize(nodeIndex+1);
                }
                DecisionTreeNode& node(dtree.data[nodeIndex]);

                switch(nodeTypeIndex)
                {
                case TREE:
                    parseTreeNode(v,node.tree);
                    break;
                case VOTE:
                    parseTreeNode(v,node.vote);
                    break;
                case DECISION:
                    parseTreeNode(v,node.decision);
                    if (node.decision.left >= static_cast<int>(_nFeatures))
                    {
                        _nFeatures=node.decision.left+1;
                    }
                    break;
                default:
                    assert(false && "Unknown node type");
                }
            }
        }
    }
}



double
RandomForestModel::
getDecisionTreeProb(
    const feature_type& features,
    const DecisionTree& dtree) const
{
    unsigned nodeIndex(0);

    //traverse a single tree
    while(true)
    {
        const DecisionTreeNode& node(dtree.getNode(nodeIndex));
        assert(node.tree.isInit);
        // test condition signifies we've reached a leaf node
        if (node.tree.left == -1)
        {
            assert(node.vote.isInit);
            const double total = node.vote.left + node.vote.right;
            return (node.vote.right / total);
        }

        assert(node.decision.isInit);
        if (features.at(node.decision.left) <= node.decision.right)
        {
            nodeIndex=node.tree.left;
        }
        else
        {
            nodeIndex=node.tree.right;
        }
    }
}



double
RandomForestModel::
getProb(
    const feature_type& features) const
{
    double retval(0);
    try
    {
        // get the probability for every tree and average them out.
        double prob(0);
        for (const DecisionTree& dtree : _forest)
        {
            prob += getDecisionTreeProb(features, dtree);
        }
        retval = (1-(prob/_forest.size()));
    }
    catch (...)
    {
        log_os << "Exception caught in random forest while scoring features:\n";
        for (const auto val : features)
        {
            log_os << "K:V " << val.first << " : " << val.second << "\n";
        }
        throw;
    }
    return retval;
}


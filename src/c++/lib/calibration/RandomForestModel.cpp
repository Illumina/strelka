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
 *      Author: Morten Kallberg
 */

#include "RandomForestModel.hh"

#include "blt_util/log.hh"
#include "blt_util/parse_util.hh"

#include <cstdlib>

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

static
const std::string
get_label(const index_t i)
{
    switch (i)
    {
    case TREE:
        return "tree";
    case VOTE:
        return "node_votes";
    case DECISION:
        return "decisions";
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
    const Json::Value& v,
    TreeNode<L,R>& val)
{
    assert(! val.isInit);
    assert(v.isArray());
    assert(v.size() == 2);

    val.isInit = true;
    val.left = static_cast<L>(v[0].asDouble());
    val.right = static_cast<R>(v[1].asDouble());
}

void RandomForestModel
::Deserialize( const Json::Value& root)
{
    clear();

    using namespace DTREE_NODE_TYPE;
    serialized_calibration_model::Deserialize(root);

    //TODO read in other RF specific values

    // Loop through all the trees
    Json::Value jmodels = root["Model"];
    for ( unsigned tree_count = 0; tree_count < jmodels.size(); ++tree_count)
    {
        _forest.emplace_back();
        DecisionTree& dtree(_forest.back());

        // loop through the three paramter categories (TREE,VOTE, DECISION) for each tree
        for (int i(0); i<SIZE; ++i)
        {
            const index_t nodeTypeIndex(static_cast<index_t>(i));
            Json::Value tree_cat = jmodels[tree_count][get_label(nodeTypeIndex)];

            // loop over all the nodes in the tree
            for (Json::Value::iterator it = tree_cat.begin(); it != tree_cat.end(); ++it)
            {
                unsigned nodeIndex(std::atoi(it.key().asCString()));
                if (dtree.data.size() <= nodeIndex)
                    dtree.data.resize(nodeIndex+1);

                Json::Value value = (*it);
                DecisionTreeNode& node(dtree.data[nodeIndex]);
                switch (nodeTypeIndex)
                {
                case TREE:
                    parseTreeNode(value,node.tree);
                    break;
                case VOTE:
                    parseTreeNode(value,node.vote);
                    break;
                case DECISION:
                    parseTreeNode(value,node.decision);
                    if (node.decision.left >= static_cast<int>(_nFeatures))
                    {
                        _nFeatures=node.decision.left+1;
                    }
                    break;
                default:
                    assert(false && "Unknown node type when reading in Random Forrest model");
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
    while (true)
    {
        const DecisionTreeNode& node(dtree.getNode(nodeIndex));
        assert(node.tree.isInit);
        // test condition signifies we've reached a leaf node
        if (node.tree.left == -1)
        {
            assert(node.vote.isInit);
            const double total = node.vote.left + node.vote.right;
            return (node.vote.left / total);
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
        retval = prob/_forest.size();
    }
    catch (...)
    {
//        log_os << "Exception caught in random forest while scoring features:\n";
//        for (const auto val : features)
//        {
//            log_os << "K:V " << val.first << " : " << val.second << "\n";
//        }
        throw;
    }
    return retval;
}

//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
/*
 *      Author: Morten Kallberg
 */

#include "RandomForestModel.hh"

#include "blt_util/log.hh"
#include "blt_util/parse_util.hh"
#include "common/Exceptions.hh"

#include <sstream>



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



void
RandomForestModel::
Deserialize(
    const unsigned expectedFeatureCount,
    const Json::Value& root)
{
    clear();

    using namespace DTREE_NODE_TYPE;

    //TODO read in other RF specific values

    // Loop through all the trees
    const Json::Value jmodels = root["Model"];
    const unsigned treeCount(jmodels.size());
    for (unsigned treeIndex = 0; treeIndex < treeCount; ++treeIndex)
    {
        _forest.emplace_back();
        DecisionTree& dtree(_forest.back());

        // loop through the three parameter categories (TREE,VOTE, DECISION) for each tree
        for (int i(0); i<SIZE; ++i)
        {
            const index_t nodeTypeIndex(static_cast<index_t>(i));
            const Json::Value tree_cat = jmodels[treeIndex][get_label(nodeTypeIndex)];

            // loop over all the nodes in the tree
            for (Json::Value::const_iterator it = tree_cat.begin(); it != tree_cat.end(); ++it)
            {
                using namespace illumina::blt_util;
                const unsigned nodeIndex(parse_unsigned_rvalue(it.key().asCString()));
                if (dtree.data.size() <= nodeIndex) dtree.data.resize(nodeIndex+1);

                const Json::Value value = (*it);
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
                    if (node.decision.left >= static_cast<int>(expectedFeatureCount))
                    {
                        using namespace illumina::common;

                        std::ostringstream oss;
                        oss << "ERROR: scoring model max feature index: " << node.decision.left
                            << " is inconsistent with expected feature count " << expectedFeatureCount;
                        BOOST_THROW_EXCEPTION(LogicException(oss.str()));
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
    const featureInput_t& features,
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
        if (features[node.decision.left] <= node.decision.right)
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
    const featureInput_t& features) const
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

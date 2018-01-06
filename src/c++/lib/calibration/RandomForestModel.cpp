//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2018 Illumina, Inc.
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
const char*
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



static
void
missingNodeError(
    const char* key)
{
    std::ostringstream oss;
    oss << "Can't find expected node '" << key << "' in  json scoring model file.";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
}



static
void
wrongValueTypeError(
    const char* key,
    const char* keyType)
{
    std::ostringstream oss;
    oss << "Node '" << key << "' in json scoring model file does not have expected type '" << keyType << "'.";
    BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
}



template <typename L, typename R>
void
RandomForestModel::
parseTreeNode(
    const rapidjson::Value& v,
    TreeNode<L, R>& val)
{
    assert(! val.isInit);
    assert(v.IsArray());
    assert(v.Size() == 2);

    val.isInit = true;
    val.left = static_cast<L>(v[0].GetDouble());
    val.right = static_cast<R>(v[1].GetDouble());
}



void
RandomForestModel::
Deserialize(
    const unsigned expectedFeatureCount,
    const rapidjson::Value& root)
{
    clear();

    using namespace DTREE_NODE_TYPE;

    auto getNodeMember = [](const rapidjson::Value& node, const char* label) -> const rapidjson::Value&
    {
        const rapidjson::Value::ConstMemberIterator iter(node.FindMember(label));
        if (iter == node.MemberEnd()) missingNodeError(label);
        return iter->value;
    };

    // Loop through all the trees
    static const char* modelLabel = "Model";
    const rapidjson::Value& rfTreeArray(getNodeMember(root, modelLabel));
    if (! rfTreeArray.IsArray()) wrongValueTypeError(modelLabel, "array");

    for (const auto& treeValue : rfTreeArray.GetArray())
    {
        _forest.emplace_back();
        DecisionTree& decisionTree(_forest.back());

        // loop through the three parameter categories (TREE, VOTE, DECISION) for each tree
        for (int i(0); i<SIZE; ++i)
        {
            const index_t nodeTypeIndex(static_cast<index_t>(i));
            const rapidjson::Value& treeCategoryValue(getNodeMember(treeValue, get_label(nodeTypeIndex)));

            // loop over all the nodes in the tree
            for (const auto& treeNode : treeCategoryValue.GetObject())
            {
                using namespace illumina::blt_util;
                const unsigned decisionTreeNodeIndex(parse_unsigned_rvalue(treeNode.name.GetString()));
                if (decisionTree.data.size() <= decisionTreeNodeIndex)
                {
                    decisionTree.data.resize(decisionTreeNodeIndex+1);
                }

                DecisionTreeNode& decisionTreeNode(decisionTree.data[decisionTreeNodeIndex]);
                switch (nodeTypeIndex)
                {
                case TREE:
                    parseTreeNode(treeNode.value, decisionTreeNode.tree);
                    break;
                case VOTE:
                    parseTreeNode(treeNode.value, decisionTreeNode.vote);
                    break;
                case DECISION:
                    parseTreeNode(treeNode.value, decisionTreeNode.decision);
                    if (decisionTreeNode.decision.left >= static_cast<int>(expectedFeatureCount))
                    {
                        std::ostringstream oss;
                        oss << "ERROR: scoring model max feature index: " << decisionTreeNode.decision.left
                            << " is inconsistent with expected feature count " << expectedFeatureCount;
                        BOOST_THROW_EXCEPTION(illumina::common::GeneralException(oss.str()));
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

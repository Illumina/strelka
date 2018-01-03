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

#pragma once

#include "VariantScoringModelBase.hh"

#include "rapidjson/document.h"

#include <cassert>

#include <map>
#include <vector>

struct RandomForestModel : public VariantScoringModelBase
{
    RandomForestModel() = default;

    double getProb(const featureInput_t& features) const override;

    void Deserialize(
        const unsigned expectedFeatureCount,
        const rapidjson::Value& root);

private:
    template <typename L, typename R>
    struct TreeNode
    {
        TreeNode() : isInit(false) {}

        bool isInit;
        L left;
        R right;
    };

    struct DecisionTreeNode
    {
        TreeNode<int,int> tree;
        TreeNode<double,double> vote;
        TreeNode<int,double> decision; // (feature index, feature value)
    };

    struct DecisionTree
    {
        const DecisionTreeNode&
        getNode(const unsigned i) const
        {
            assert(i<data.size());
            return data[i];
        }

        std::vector<DecisionTreeNode> data;
    };


    template <typename L, typename R>
    void
    parseTreeNode(
        const rapidjson::Value& v,
        TreeNode<L, R>& val);

    double
    getDecisionTreeProb(
        const featureInput_t& features,
        const DecisionTree& dtree) const;

    void
    clear()
    {
        _forest.clear();
    }

////////data:
    std::vector<DecisionTree> _forest;
};

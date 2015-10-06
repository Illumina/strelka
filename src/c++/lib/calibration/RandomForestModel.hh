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

#pragma once

#include <calibration/SerializedModel.hh>

#include <vector>
#include <cassert>


struct RandomForestModel : public serialized_calibration_model
{
    RandomForestModel() {}


    bool isInit() const
    {
        return (! _forest.empty());
    }

    unsigned expectedFeatureCount() const
    {
        return _nFeatures;
    }

    double getProb(const feature_type& features) const;

    void Deserialize( const Json::Value& root);

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
        const Json::Value& v,
        TreeNode<L,R>& val);

    double
    getDecisionTreeProb(
        const feature_type& features,
        const DecisionTree& dtree) const;

    void
    clear()
    {
        _nFeatures=0;
        _forest.clear();
    }

////////data:
    unsigned _nFeatures = 0;
    std::vector<DecisionTree> _forest;
};

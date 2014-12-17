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

#pragma once

#include "calibration/CalibrationModel.hh"

#include "boost/property_tree/ptree.hpp"

#include <vector>


struct RandomForestModel
{
    RandomForestModel(){};

    void load(const boost::property_tree::ptree& pt);

    bool isInit() const
    {
        return (! _forest.empty());
    }

    unsigned expectedFeatureCount() const
    {
        return _nFeatures;
    }

    double getProb(const feature_type& features) const;

private:
    template <typename L, typename R>
    struct TreeNode
    {
        bool isInit = false;
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
        const boost::property_tree::ptree::value_type& v,
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

#
# Starka
# Copyright (c) 2009-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

import os
import sys
import cPickle
import json


def dict_for_tree(tree):
    """ Make dictionary for tree """
    tree_dict = {}
    for i in range(tree.node_count):
        tree_dict[i] = [int(tree.children_left[i]), int(tree.children_right[i])]
    node_features = tree.feature
    node_thresholds = tree.threshold
    node_count = tree.node_count

    decision_dict = {}
    for i in range(node_count):
        if tree_dict[i][0] == -1 :
            node_features[i] = -1
            node_thresholds[i] = -1

        decision_dict[i] = [int(node_features[i]), float(node_thresholds[i])]
    node_values = tree.value

    node_values = node_values.tolist()
    node_values = [i[0] for i in node_values]
    node_values_dict = {}
    for i in range(len(node_values)):
        node_values_dict[i] = node_values[i]
    all_data = {}

    all_data['tree'] = tree_dict
    all_data['decisions'] = decision_dict
    all_data['node_votes'] = node_values_dict

    return all_data


def read_pickled_classifier(pickle_fn):
    """ Read classifier from pickle file """
    with open(pickle_fn, 'rb') as fid:
        clf = cPickle.load(fid)
    return clf


def write_classifier_pickle(clf, pickle_fn):
    """ store classifier in pickle file """
    with open(pickle_fn, 'wb') as fid:
        cPickle.dump(clf, fid)


def classifier_to_dict(clf):
    """ Convert random forest classifier to  """
    all_trees = []
    for i in range(len(clf.estimators_)):
        all_trees.append(dict_for_tree(clf.estimators_[i].tree_))
    return all_trees


def write_classifier_json(clf, fn):
    """ Write classifier as JSON file """
    all_trees = classifier_to_dict(clf)
    if not fn.endswith(".json"):
        fn += ".json"
    with open(fn, 'w') as fid:
        json.dump(all_trees, fid)

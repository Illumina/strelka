# coding=utf-8

import numpy as np
import pandas
from vqsr import VQSRModel

import vqsr.tools.io as io

from sklearn.ensemble import RandomForestClassifier


class StrelkaRF(VQSRModel):

    def __init__(self):
        self.clf = None


    def train(self, tp, fp, columns, *args, **kwargs):
        """ Train model from sets of TPs and FPs

        :param tp: data frame with rows of TP instances.
        :type tp: pandas.DataFrame
        :param fp: data frame with rows of FP instances.
        :type fp: pandas.DataFrame
        :param columns: the feature columns to use
        """

        tdf = pandas.DataFrame(tp)
        fdf = pandas.DataFrame(fp)

        # TODO try numeric labels
        tdf["tag"] = "TP"
        fdf["tag"] = "FP"

        allrows = pandas.concat([tdf, fdf])

        # TODO: parameters to try
        # {'max_depth' : range(1,20,1), 'n_estimators' : range(5, 30+1,5)},

        if not kwargs:
            kwargs = {"n_jobs": 8,
                      # "max_depth": 6,
                      # "min_samples_leaf": 100,
                      "n_estimators": 100,
                      # "max_features": None
                      }

        self.clf = RandomForestClassifier(**kwargs)
        self.clf.fit(allrows[columns].values, allrows["tag"].values)


    def classify(self, instances, columns, *args, **kwargs):
        """ Classify a set of instances after training

        :param instances: data frame with instances.
        :type instances: pandas.DataFrame
        :param columns: the feature columns to use

        :return: data frame with added column "tag" which gives the classification
        :rtype: pandas.DataFrame
        """
        instances = pandas.DataFrame(instances)
        preds = self.clf.predict(instances[columns].values)
        instances["ptag"] = preds

        cprobs = self.clf.predict_proba(instances[columns].values)
        tpcol = np.where(self.clf.classes_ == "TP")[0][0]
        instances["qual"] = cprobs[:, tpcol]

        return instances


    def save(self, filename):
        """ Save to file """
        if filename.endswith(".json"):
            io.write_classifier_json(self.clf, filename)
        else:
            io.write_classifier_pickle(self.clf, filename)

    def load(self, filename):
        """ Load from file """
        self.clf = io.read_pickled_classifier(filename)


    def plots(self, prefix, featurenames):
        """ Make diagnostic plots """
        importances = self.clf.feature_importances_
        std = np.std([tree.feature_importances_ for tree in self.clf.estimators_],
                      axis=0)
        indices = np.argsort(importances)[::-1]

        # Print the feature ranking
        print "Feature ranking:"

        for f in xrange(0, len(indices)):
            print "%d. feature %d:%s (%f +- %f)" % (f + 1, indices[f],
                                                    featurenames[indices[f]],
                                                    importances[indices[f]],
                                                    std[indices[f]])

VQSRModel.register("strelka.rf", StrelkaRF)

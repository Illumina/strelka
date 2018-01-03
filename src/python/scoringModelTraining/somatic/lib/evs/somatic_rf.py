#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

# coding=utf-8

import numpy as np
import pandas
from evs import EVSModel

import evs.tools.io as io

from sklearn.ensemble import RandomForestClassifier


class SomaticRF(EVSModel):

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
        tdf["tag"] = "TP"
        fdf["tag"] = "FP"
        allrows = pandas.concat([tdf, fdf])

        if not kwargs:
            kwargs = {"n_jobs": 8,
                      "max_depth": 6,
                      # "min_samples_leaf": 100,
                      "n_estimators": 100,
                      # "max_features": None
                      }

        self.clf = RandomForestClassifier(**kwargs)
        self.clf.fit(allrows[columns].values, allrows["tag"].values, allrows["weight"].values)

        # add audit trail into the model output:
        self.clf.columns = columns


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

    def save_json_strelka_format(self, filename, varianttype, threshold):
        """ Save to json including all strelka scoring model meta-data """
        import datetime
        import json

        date = datetime.datetime.utcnow().isoformat()
        meta = {
                "Date" : "%sZ" % (date),
                "Features" : self.clf.columns,
                "ModelType" : "RandomForest",
                "FilterCutoff" : threshold,
                "Calibration" : {"Power" : 1, "Scale" : 1}
                 }
        all_trees = io.classifier_to_dict(self.clf)
        full_model = meta
        full_model["Model"] = all_trees
        modelFile = {"CalibrationModels" : {"Somatic" : { varianttype : full_model }}}
        json.dump(modelFile, open(filename,"wb"))

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

EVSModel.register("somatic.rf", SomaticRF)

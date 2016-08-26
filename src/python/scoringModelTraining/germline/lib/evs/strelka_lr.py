#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2016 Illumina, Inc.
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

from sys import exit
import numpy as np
import pandas
from evs import EVSModel

import evs.tools.io as io

from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline

class StrelkaLR(EVSModel):

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

        allrows = pandas.concat([tdf, fdf]).dropna(axis=0,how='any',subset=columns)
 #       print(allrows.dtypes)
 #       scaler = preprocessing.StandardScaler().fit(allrows[columns].values)
 #       print("Data mean: ", allrows[columns].mean(axis=0))
 #       print("Scaler mean: ", scaler.mean_)
 #       print("Data std: ", allrows[columns].std(axis=0))
 #       print("Scaler scale: ", scaler.scale_)
#
#        dat = scaler.transform(allrows[columns].values)
#        print("Transformed mean: ", dat.mean(axis=0))
#        print("Transformed std: ", dat.std(axis=0))
#        raise SystemExit

        if not kwargs:
#            kwargs = {"loss": "log",
#                      "penalty": "l1",
#                      "l1_ratio": 1,
#                      "verbose": 1,
##                      "alpha": 0.0005,
#                      "alpha": 0.0001, # 10^-4: hardcoded value for SNP HET and SNP HOM on this dataset
#                      "n_iter": 100
#                      }

            kwargs = {"n_jobs": 8,
                      "penalty": "l1",
                      "dual": False,
                      "solver": "liblinear",
                      "verbose": 1,
#                      "C": 0.0005, 
                      "C": 0.0001, # 10^-4: hardcoded value for SNP HET and SNP HOM on this dataset                                                          
                      "tol": 0.01,
                      "max_iter": 100
                      }


#        self.clf = make_pipeline(StandardScaler(), LogisticRegression(**kwargs))
#        self.clf = make_pipeline(StandardScaler(), PolynomialFeatures(degree=2, interaction_only=True), SGDClassifier(**kwargs))
        self.clf = make_pipeline(StandardScaler(), PolynomialFeatures(degree=2, interaction_only=True, include_bias=False), LogisticRegression(**kwargs))
        self.clf.fit(allrows[columns].values, allrows["tag"].values)

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
        instances = pandas.DataFrame(instances).dropna(axis=0,how='any',subset=columns)
        preds = self.clf.predict(instances[columns].values)
        instances["ptag"] = preds

        cprobs = self.clf.predict_proba(instances[columns].values)
        tpcol = np.where(self.clf.classes_ == "TP")[0][0]
        instances["qual"] = cprobs[:, tpcol]

        return instances

    def save_json_strelka_format(self, filename):
        """ Save to json including all strelka scoring model meta-data """
        import datetime
        import json

        date = datetime.datetime.utcnow().isoformat()
        meta = {
                "Date" : "%sZ" % (date),
                "Features" : self.clf.columns,
#                "ModelType" : "SGDClassifier",
                "ModelType" : "LogisticRegression",
                "FilterCutoff" : 0.5
                }
        all_trees = io.classifier_to_dict(self.clf)
        #full_model = {"Meta" : meta, "Model" : all_trees }
        full_model = meta
        full_model["Model"] = all_trees
        modelFile = {"CalibrationModels" : {"germline_lr" : { "SNP" : full_model }}}
        json.dump(modelFile, open(filename,"wb"))

    def plots(self, prefix, featurenames):
        """ Make diagnostic plots """
#        coefficients = self.clf.named_steps['sgdclassifier'].coef_
        coefficients = self.clf.named_steps['logisticregression'].coef_
        indices = np.argsort(abs(coefficients[0]))[::-1]
        powers = self.clf.named_steps['polynomialfeatures'].powers_

        # Print the feature ranking
        print "Feature ranking:"
        for f in xrange(0, len(indices)):
            combinedfeatures = []
            for ff in xrange(0, len(featurenames)):
                if powers[indices[f]][ff] > 0 :
                    combinedfeatures.append(featurenames[ff])
            print "%d. feature %d:%s (%f)" % (f + 1, indices[f],
                                                    '.'.join(combinedfeatures),
                                                    coefficients[0][indices[f]])

#            print "%d. feature %d:%s (%f)" % (f + 1, indices[f],                   # 
#                                                    featurenames[f],     
#                                                    coefficients[0][indices[f]])


EVSModel.register("strelka.lr", StrelkaLR)

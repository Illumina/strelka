#!/illumina/thirdparty/python/python-2.7.5/bin/python
#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2017 Illumina, Inc.
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

########################################
# Tool for aggregating multiple indel error sampling result and 
#
# Author: Morten Kallberg
# Date: March 3rd, 2014
#
########################################

import os
import sys

import cPickle as pickle
import argparse
from dynamicCaseCounter import job,caseCounter
from numpy import mean, sqrt, square, arange
import numpy as np
import logging
import json
import glob
import pickle
from collections import defaultdict
import datetime
import time

DEBUG=0

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inputfolder', dest="tempFolder", help='Temp folder containing individual count outputs')
parser.add_argument('-o', '--outputModel', dest="outputModel", help='Destination for output model')
parser.add_argument('-m', '--inputModel', dest="inputModel", help='Input model JSON file')
parser.add_argument('-a', '--addNewModel', dest="matchExistingModel", action='store_false',
                    help='Add model to JSON instead of comparing to existing models')
parser.add_argument('-M', '--matchExistingModel', dest="matchExistingModel", action='store_true',
                    default=True,
                    help='Do not add model, just compare to existing models in input model JSON')
parser.add_argument('-n', "--modelName", dest="modelName",
                    help="When adding a new model, specify what it should be called")
parser.add_argument('-f', "--forceOverwriteModel", dest="forceOverwriteModel", action="store_true",
                    default=False,
                    help="If a model with modelName already exists, overwrite it")

args = parser.parse_args()
# setup logger
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logging.info("This is the indel aggregator script, we are in business")

if args.matchExistingModel and args.modelName is not None:
    logging.warning("Specified a model name even though no new model is being added")

if args.matchExistingModel and not os.path.isfile(args.inputModel):
    logging.error(args.inputModel + " is empty, cannot match to existing model")
    sys.exit("dynamicModelAggregator exited with an error")

if not args.matchExistingModel and args.modelName is None:
    logging.error("You must specify a model name if you want to add a new model")
    sys.exit("dynamicModelAggregator exited with an error")


# read in model file 
# read in a json file
if DEBUG:
    args.inputModel = "/home/mkallberg/workspace/starka/src/config/indel_models.json"
    args.outputModel = "/home/mkallberg/workspace/starka_debug/bams/error_model/test/indel_models_fit.json"
    args.tempFolder = "/home/mkallberg/workspace/starka_debug/bams/error_model/test/"


if args.inputModel is not None and os.path.isfile(args.inputModel):
    inputModelFile = open(args.inputModel)
    inputModelJson = json.load(inputModelFile) 
    inputModelFile.close()

# aggregate files results
# read in pickled count jobs from dynamicCaseCounter
# TODO eliminate 10% most extreme cases to avoid skewed estimate from unfortunate regions choices

masterCounter = None
for p in glob.glob(os.path.join(args.tempFolder,'*.pickle')):
    f = open(p)
    caseCounter = pickle.load(f)
    if not masterCounter: masterCounter = caseCounter
    else:   
        masterCounter.accumulateCounts(caseCounter)

# pickle the mastercounter to reference back to 
masterFile = open(os.path.join(args.tempFolder,'master.pickle'),'w')
pickle.dump(masterCounter, masterFile)
masterFile.close()

# model selection
# look at cumulative evidence from multiple region count and select model from input json
if args.inputModel is not None:
    if args.matchExistingModel and os.path.isfile(args.inputModel):
        # get master estimate for deletions, use first eight entries
        masterCounter.calcPars()
        # print masterCounter.parsSTR
        masterModel = masterCounter.parsSTR[1]['del']['simple']['props'][:8]

        # logging.info("Matching to existing models")
        # modelChoice = inputModelJson['IndelModels'].keys()[:1]     # pick first model as default
        # minRms = sys.maxint
        # for model in inputModelJson['IndelModels'].keys():
        #     inputModel = np.array([t[0] for t in inputModelJson['IndelModels'][model][:8]])
        #     diff = list(np.array(masterModel) - inputModel)
        #     rms = sqrt(mean(square(diff)))
        #     if rms<minRms:
        #         minRms = rms
        #         modelChoice = model
            
        # output model file
        # return chosen model by setting parameter in master pyflow....
        outputModelFile = open(opt.outputModel,'w')
        outputModelJson = inputModelJson
        outputModelJson['IndelModels']['bestFit'] = modelChoice
        data_string = json.dumps(outputModelJson)
        outputModelFile.write(data_string)
        outputModelFile.close()
    else:
        # check to see if model is already in inputModelJson and whether we've given
        # the go-ahead to overwrite it
        if os.path.isfile(args.inputModel):
            versionTag = "v1"
            if args.modelName in inputModelJson['IndelModels'].keys():
                ## need to increment version tag if model name is already in json
                versionTag = "v2"

        logging.info("Adding new model %s-%s" % (args.modelName, versionTag))

        masterCounter.calcPars()
        masterCounter.calcPars("ins")
        data = {}
        data["Name"] = args.modelName
        data["Version"] = versionTag

        ts = time.time()
        data["Date"] = datetime.datetime.fromtimestamp(ts).strftime("%m/%d/%Y")

        data["MaxMotifLength"] = len(masterCounter.parsSTR)
        data["MaxTractLength"] = 34
        data["Model"] = []
        tempErrorRates  = []
        for ml in masterCounter.parsSTR.keys():
            for indelType in sorted(masterCounter.parsSTR[ml].keys()):
                errorRate = masterCounter.parsSTR[ml][indelType]['simple']['props']
                currList = []
                for length, rate in enumerate(errorRate):
                    if rate != 0.1:
                        currList.append(rate/(length + 1))
                    else:
                        currList.append(0.)
                tempErrorRates.append(currList)
            data["Model"].append(map(list, zip(*tempErrorRates)))
            tempErrorRates = []

        print data

        if not os.path.isfile(args.inputModel):
            inputModelJson = defaultdict(dict)
        inputModelJson['IndelModels'][args.modelName] = data
        data_string = json.dumps(inputModelJson)
        # inputModelFile = open(args.inputModel, "w")
        outputModelFile = open("/home/mbekritsky/workspace/starka/error_model/test.json", "w")
        outputModelFile.write(data_string)
        outputModelFile.close()

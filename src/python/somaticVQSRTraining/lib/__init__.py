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

# coding=utf-8

import abc


class VQSRModel(object):
    """ Base class for VQSR models """

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def train(self, tp, fp, columns, *args, **kwargs):
        """ Train model from sets of TPs and FPs

        :param tp: data frame with rows of TP instances.
        :type tp: pandas.DataFrame
        :param fp: data frame with rows of FP instances.
        :type fp: pandas.DataFrame
        :param columns: the feature columns to use

        """
        pass

    # noinspection PyUnusedLocal
    @abc.abstractmethod
    def classify(self, instances, columns, *args, **kwargs):
        """ Classify a set of instances after training

        :param instances: data frame with instances.
        :type instances: pandas.DataFrame
        :param columns: the feature columns to use

        :return: data frame with added column "tag" which gives the classification
        :rtype: pandas.DataFrame
        """
        instances["tag"] = "FP"
        return instances

    @abc.abstractmethod
    def save(self, filename):
        """ Save to file """
        pass

    @abc.abstractmethod
    def load(self, filename):
        """ Load from file """
        pass

    # model factory
    _models = {}

    @classmethod
    def register(cls, mname, mclass):
        cls._models[mname] = mclass

    @classmethod
    def create(cls, mname):
        return cls._models[mname]()

    @classmethod
    def names(cls):
        return cls._models.keys()


import strelka_rf   # noqa
import strelka_rf_indel   # noqa

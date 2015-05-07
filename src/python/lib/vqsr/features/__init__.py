#!/illumina/development/haplocompare/hc-virtualenv/bin/python
# coding=utf-8
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/sequencing/licenses/blob/master/Simplified-BSD-License.txt

import abc


class FeatureSet(object):
    """ VCF paired Feature set for somatic comparison """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self.chr_depth = {}

    @abc.abstractmethod
    def collect(self, vcfname):
        """ Return a data frame with features collected from
            the given VCF, tagged by given type """
        pass

    @abc.abstractmethod
    def trainingfeatures(self):
        """ Return a list of columns that are features to use for VQSR training """
        pass

    sets = {}

    @staticmethod
    def register(name, xclass=None):
        if xclass is None:
            def fn(xclass2):
                FeatureSet.sets[name] = xclass2
            return fn
        else:
            FeatureSet.sets[name] = xclass

    @staticmethod
    def make(name):
        # noinspection PyCallingNonCallable
        return FeatureSet.sets[name]()

    def setChrDepths(self, cd):
        """ set depth normalisation factors (can come from VCF or BAM) """
        self.chr_depth = cd

import StrelkaSNV   # noqa
import StrelkaIndel  # noqa
import PosAndAlleles  # noqa

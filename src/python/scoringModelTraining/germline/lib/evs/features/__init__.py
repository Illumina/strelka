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

import abc


class FeatureSet(object):
    """ VCF paired Feature set for somatic comparison """

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        pass

    @abc.abstractmethod
    def collect(self, vcfname):
        """ Return a data frame with features collected from
            the given VCF, tagged by given type """
        pass

    @abc.abstractmethod
    def trainingfeatures(self):
        """ Return a list of columns that are features to use for EVS model training """
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


import GermlineSNV   # noqa
import GermlineIndel   # noqa
import RNASNV   # noqa
import RNAIndel   # noqa

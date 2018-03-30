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

import os,sys
import unittest

scriptDir=os.path.abspath(os.path.dirname(__file__))
targetDir=os.path.abspath(os.path.join(scriptDir,"@THIS_RELATIVE_PYTHON_LIBDIR@"))
sys.path.append(targetDir)

import configureUtil


class test_configureUtil(unittest.TestCase):

    def test_safeSetBool(self):
        class Foo :
            pass

        configureUtil.safeSetBool(Foo,"bar")
        self.assertFalse(Foo.bar)

    def test_checkForBamExtension(self):
        try :
            configureUtil.checkForBamExtension("foo.foo")
        except configureUtil.OptParseException:
            pass
        else :
            self.fail("Method did not raise OptParseException")

        try :
            configureUtil.checkForBamExtension("foo.bam")
        except :
            self.fail("Method raised unexpected exception")

        try :
            configureUtil.checkForBamExtension("foo.cram")
        except :
            self.fail("Method raised unexpected exception")

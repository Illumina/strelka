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

from checkChromSet import ordinalStr

class test_checkChromSet(unittest.TestCase):

    def test_ordinalStr(self):
        self.assertEqual(ordinalStr(1), '1st')
        self.assertEqual(ordinalStr(2), '2nd')
        self.assertEqual(ordinalStr(3), '3rd')
        for n in range(4, 10):
            self.assertEqual(ordinalStr(n), str(n)+'th')

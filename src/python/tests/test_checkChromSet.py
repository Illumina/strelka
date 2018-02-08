import os,sys
import unittest

scriptDir=os.path.abspath(os.path.dirname(__file__))
scriptName=os.path.basename(__file__)

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

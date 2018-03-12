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

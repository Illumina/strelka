########################################
# Pyflow for running errormodel estimator in parallele 
#
# Author: Morten Kallberg
# Date: September 16th, 2014
#
########################################


pyflowDir=os.path.join(scriptDir,"pyflow")
sys.path.append(os.path.abspath(pyflowDir))

from pyflow import WorkflowRunner

class indelErrorWorkflow(WorkflowRunner) :

    def __init__(self,params,paths) :
        self.params = params
        self.paths = paths

    def workflow(self) :
        sf = 1 
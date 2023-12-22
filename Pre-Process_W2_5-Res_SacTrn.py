
from hec.heclib.dss import HecDss
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.heclib.util import HecTime
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
import hec.hecmath.TimeSeriesMath as tsmath
from rma.util.RMAConst import MISSING_DOUBLE
import math
import sys
import datetime as dt
import os, sys

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project

sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import Forecast_preprocess as fpp
reload(fpp)

import Acc_Dep_ResSim_SacTrn
reload(Acc_Dep_ResSim_SacTrn)

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    data_preprocess = fpp.forecast_data_preprocess_W2_5Res(currentAlternative, computeOptions)

	# compute acc TODO
    #acc_dep = Acc_Dep_ResSim_SacTrn.computeAlternative(currentAlternative, computeOptions)

    if data_preprocess:# and acc_dep:
    	return True


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
#import hec.hecmath.TimeSeriesMath as tsmath
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))


from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

import Acc_Dep_ResSim_SacTrn
reload(Acc_Dep_ResSim_SacTrn)

import DMS_preprocess
reload(DMS_preprocess)


def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    data_preprocess = DMS_preprocess.preprocess_ResSim_5Res(currentAlternative, computeOptions)

    acc_dep = Acc_Dep_ResSim_SacTrn.computeAlternative(currentAlternative, computeOptions)

    if data_preprocess and acc_dep:
        return True

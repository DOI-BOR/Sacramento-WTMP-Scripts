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
import os, sys, time

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project

# print current path
print("Current paths: ", sys.path)

# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento", "American", "Stanislaus"]

# initialize and search for unwanted paths
matching_paths = []
for p in sys.path:
    if any(phrase in p for phrase in search_list):
        matching_paths.append(p)

# print paths containing unwanted phrases
print("Paths to be removed:")
for path in matching_paths:
    print(path)

# remove matching paths from sys.path
for path in matching_paths:
    if path in sys.path:
        sys.path.remove(path)

# append path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))


def computeAlternative(currentAlternative, computeOptions):
    """
    This just logs a message and pauses for 10 seconds to allow files to unlock

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `getOutputDataLocations()` for logging.
    computeOptions : object
        The compute options/settings object.

    Returns
    -------
    bool
        True once the waiting is over.

    """

    currentAlternative.addComputeMessage("Pausing to allow for DSS files to unlock")
    time.sleep(10)

    return True
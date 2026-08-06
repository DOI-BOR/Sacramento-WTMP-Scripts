"""
Pre-Process_ResSim_SacTrnPO
=============================

Compute a WAT Scripting Alternative that runs 5-reservoir ResSim
preprocessing followed by accretion/depletion calculations for the
Sacramento-Trinity (SacTrn) model.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `hec.heclib.dss`, `hec.hecmath`, `hec.heclib.util`,
  `hec.io`, `Acc_Dep_ResSim_SacTrn`, `DMS_preprocess`.
"""

from hec.heclib.dss import HecDss                     # HEC-DSS: open/read/write DSS files (imported for completeness; not directly used)
from hec.hecmath import HecMathException              # HEC math exception type (imported for completeness; not directly used)
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE   # HEC sentinel for undefined values (imported for completeness)
from hec.heclib.util import HecTime                   # HEC time object (imported for completeness; not directly used)
from hec.io import DSSIdentifier                      # DSS path identifier helper (not directly used)
from hec.io import TimeSeriesContainer                # Container class (imported for completeness; not directly used)
import hec.hecmath.TimeSeriesMath as tsmath           # Time series math wrapper (imported for completeness; not directly used)
from rma.util.RMAConst import MISSING_DOUBLE          # RMA sentinel for missing values (imported for completeness)
import math                                           # Standard math library (imported for completeness; not directly used)
import sys                                            # System utilities: path manipulation
import datetime as dt                                 # Python datetime (imported for completeness; not directly used)
import os, sys                                        # Filesystem and system path utilities (combined import per original)

from com.rma.io import DssFileManagerImpl             # RMA DSS manager (imported for completeness; not directly used)
from com.rma.model import Project                     # WAT project API: resolve workspace path for scripts folder
#import hec.hecmath.TimeSeriesMath as tsmath          # Original commented import (kept intact)

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
    # if the path is still present in sys.path, remove it
    if path in sys.path:
        sys.path.remove(path)

# append path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))


from com.rma.io import DssFileManagerImpl             # Duplicate import retained per original file content
from java.util import TimeZone                        # Java TimeZone (retained; not directly used)

import Acc_Dep_ResSim_SacTrn                          # Local module: reservoir balance-flow ("accretion/depletion") computations
reload(Acc_Dep_ResSim_SacTrn)                         # Jython: ensure latest version is loaded

import DMS_preprocess                                 # Local module: DMS unit/type standardization and derived-record generation
reload(DMS_preprocess)                                # Reload to ensure latest version


def computeAlternative(currentAlternative, computeOptions):
    """
    Compute a scripting alternative by running preprocessing and
    accretion/depletion calculations for the ResSim Sacramento
    Trinity (SacTrn) model.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging status messages during the
        compute process.
    computeOptions : object
        The compute options/settings object passed through to the
        preprocessing and accretion/depletion steps.

    Returns
    -------
    bool or None
        True if both the preprocessing step and the
        accretion/depletion step complete successfully. Returns
        None (implicitly) if either step fails, since no explicit
        False/None branch is defined.

    Notes
    -----
    This function logs a compute message to the alternative, then
    runs two dependent steps in sequence:

    1. Preprocessing of the 5-reservoir ResSim data via
       `DMS_preprocess.preprocess_ResSim_5Res`.
    2. Accretion/depletion computation via
       `Acc_Dep_ResSim_SacTrn.computeAlternative`.
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    data_preprocess = DMS_preprocess.preprocess_ResSim_5Res(currentAlternative, computeOptions)

    acc_dep = Acc_Dep_ResSim_SacTrn.computeAlternative(currentAlternative, computeOptions)

    if data_preprocess and acc_dep:
        return True
"""
Pre-Process_W2_5-Res_SacTrn
=============================

Compute a WAT Scripting Alternative that runs forecast data preprocessing for
the CE-QUAL-W2 5-reservoir Sacramento-Trinity model. A block of logic for
copying annual W2 model configuration files (mirroring the prescribed-run
scripts) and the accretion/depletion computation step are both present but
currently commented out.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `hec.heclib.dss`, `hec.hecmath`, `hec.heclib.util`,
  `hec.io`, `distutils.dir_util.copy_tree`, `Forecast_preprocess`,
  `Acc_Dep_ResSim_SacTrn`.
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
from distutils.dir_util import copy_tree              # Recursive directory copy utility (used only in commented-out block below)

from com.rma.io import DssFileManagerImpl             # RMA DSS manager (imported for completeness; not directly used)
from com.rma.model import Project                     # WAT project API: resolve workspace path for scripts folder

# print current path
print("Current paths: ", sys.path)

# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento", "American", "Stanislaus"]

# initialize and search for unwanted paths
matching_paths = []
# loop through every entry currently on sys.path
for p in sys.path:
    # if this path contains any of the unwanted phrases, mark it for removal
    if any(phrase in p for phrase in search_list):
        matching_paths.append(p)

# print paths containing unwanted phrases
print("Paths to be removed:")
for path in matching_paths:
    print(path)

# remove matching paths from sys.path
# loop through the matched paths again and remove them from sys.path
for path in matching_paths:
    if path in sys.path:
        sys.path.remove(path)

# append path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import Forecast_preprocess as fpp                     # Local module: forecast boundary-condition preprocessing steps
reload(fpp)                                            # Jython: ensure latest version is loaded

import Acc_Dep_ResSim_SacTrn                          # Local module: reservoir balance-flow computations (used only in commented-out block below)
reload(Acc_Dep_ResSim_SacTrn)                         # Reload to ensure latest version


# apparently, 2016 works with forecast, at least for Trinity, so we here make sure we are using the 2016
# "version" of the prescribed W2 models for forecast.
W2_models_for_input_copy = ['W2 Trinity Prescribed','W2 Lewiston Prescribed','W2 Keswick Prescribed']

def study_dir_from_run_dir(run_dir):
    """
    Derive the top-level study directory from a W2 run directory.

    Walks up three directory levels from the given run directory
    (run directory -> W2 simulation folder -> runs folder -> study
    folder) to find the root study directory.

    Parameters
    ----------
    run_dir : str
        Path to the W2 run directory for the current alternative.

    Returns
    -------
    str
        The path to the top-level study directory.
    """
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def annual_config_dirs_from_run_dir(run_dir,model_name,startyear_str):
    """
    Build the model, annual config, and base directory paths for a
    given W2 model and start year.

    Uses the study directory (derived from `run_dir`) along with the
    model name and start year to construct the paths to: the W2
    model's working directory, the annual configuration directory
    for the given start year, and the base directory used to store
    the original, unmodified model input files.

    Parameters
    ----------
    run_dir : str
        Path to the W2 run directory for the current alternative.
    model_name : str
        Name of the W2 model, e.g. `'W2 Trinity Prescribed'`.
    startyear_str : str
        The forecast start year as a string, e.g. `'2016'`.

    Returns
    -------
    tuple of (str, str, str)
        `(model_dir, annual_config_dir, base_dir)` paths.

    Notes
    -----
    This helper is currently unused, since the block of code that
    would call it (annual W2 configuration copying) is commented out
    in `computeAlternative` below.
    """
    study_dir = study_dir_from_run_dir(run_dir)
    model_dir = os.path.join(study_dir,'cequal-w2',model_name,model_name)  # don't know why this is two model_names deep!!
    annual_config_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,startyear_str)
    base_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,'base')
    return model_dir,annual_config_dir,base_dir

def computeAlternative(currentAlternative, computeOptions):
    """
    Compute the W2 5-Reservoir Sacramento Trinity pre-processing
    alternative.

    Parameters
    ----------
    currentAlternative : object
        The `ScriptingAlternative` being computed. Type:
        `hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt`.
    computeOptions : object
        The compute options for this run. Type:
        `hec.wat.model.ComputeOptions`.

    Returns
    -------
    bool or None
        True if the forecast data pre-processing step succeeded. If
        it fails, no explicit value is returned (implicit None),
        which HEC-WAT will typically treat as a failed alternative.

    Notes
    -----
    This function runs the forecast data pre-processing step for the
    W2 5-reservoir model (via
    `Forecast_preprocess.forecast_data_preprocess_W2_5Res`) and
    reports success only if that pre-processing step completes
    successfully.

    A block of logic for copying annual W2 model configuration files
    (based on `W2_models_for_input_copy` and `startyear_str`) is
    currently commented out below, as is the accumulated deposition
    computation step
    (`Acc_Dep_ResSim_SacTrn.computeAlternative`).
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    # startyear_str = '2016'
    
    # run_dir = computeOptions.getRunDirectory()
    # for W2_model in W2_models_for_input_copy:
        # model_dir,annual_config_dir,base_dir = annual_config_dirs_from_run_dir(run_dir,W2_model,startyear_str)
        # currentAlternative.addComputeMessage('model_dir: '+model_dir)
        # currentAlternative.addComputeMessage('annual_config_dir: '+annual_config_dir)
        # currentAlternative.addComputeMessage('base_dir: '+base_dir)
        # if not os.path.exists(annual_config_dir):
            # currentAlternative.addComputeMessage(W2_model+'- annual config not found; W2 may be configured incorrectly for this time window.')
        # else:        
            # # copy original W2 model alternative files to 'base' directory for safekeeping/later returning

            # #if not os.path.exists(base_dir):
            # #    os.makedir(base_dir)
            # #base_files = os.listdir(base_dir)
            # #if len(base_files) == 0:
            # #    currentAlternative.addComputeMessage(W2_model+'- base files not found; copying from model folder')
            # #    copy_tree(model_dir,base_dir)

            # # remove all W2 model input files EXCEPT the .w2alt file
            # for mfile in os.listdir(model_dir):
                # if not mfile.endswith('.w2Alt') and not mfile.endswith('.w2Alt.bak'):
                    # os.remove(os.path.join(model_dir,mfile))

            # # copy over annual config input files 
            # copy_tree(annual_config_dir,model_dir)
            # currentAlternative.addComputeMessage('Copied W2 inputs file for '+startyear_str+' to '+W2_model+' model alternative folder')
            
            # # now, W2 model alternative directory is configured for startyear, ready for the W2 plugin to
            # # work it's magic on the W2_con file and execute simulation   

    data_preprocess = fpp.forecast_data_preprocess_W2_5Res(currentAlternative, computeOptions)

    # compute acc TODO
    #acc_dep = Acc_Dep_ResSim_SacTrn.computeAlternative(currentAlternative, computeOptions)

    if data_preprocess:# and acc_dep:
        return True

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

import Forecast_preprocess as fpp
reload(fpp)

import DSS_Tools
reload(DSS_Tools)

model_output_and_target = [ 'model_output','target_temp']


def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    data_preprocess = fpp.forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions)

    # check if  W2 is a part of this simulation, and if so copy the target temperature to the linked record f-part
    # ------------------------------------------------------------------------------------------------------------------
    #scripting_run_dir = computeOptions.getRunDirectory()
    #run_dir,_ = os.path.split(scripting_run_dir)
    #_,run_name = os.path.split(run_dir)
    #print('run_name: ',run_name)
    #print('scripting_run_dir: ',scripting_run_dir)
    #if 'w2' in run_name.lower():
    #    currentAlternative.addComputeMessage("  'W2' found in run name, attempting to copy forecast target temperature to F-part of linked model record...")
    #    dss_file = computeOptions.getDssFilename()   
    #    locations_obj = currentAlternative.getInputDataLocations()
    #    locations_obj_2 = DSS_Tools.organizeLocations(currentAlternative, locations_obj, model_output_and_target, return_dss_paths=False)
    #    # locations_obj_2[0] is a linked output record from a model in the simulation
    #    # locations_obj_2[1] is a linked dss record, from a DSS file, NOT a model record from the simulation
    #    dss_path_tt,dss_file_tt = DSS_Tools.getDataLocationDSSInfo(locations_obj_2[1], currentAlternative, computeOptions)
    #    print(dss_file_tt,dss_path_tt)
    #    dss_model_path, _ = DSS_Tools.getDataLocationDSSInfo(locations_obj_2[0], currentAlternative, computeOptions)
    #    model_fpart = dss_model_path.split('/')[6]
    #    print('model_fpart: ',model_fpart)
    #    DSS_Tools.copy_dss_ts(dss_path_tt,
    #                          new_fpart=model_fpart,
    #                          dss_file_path=dss_file_tt,
    #                          dss_file_alt_outpath=dss_file)

    if data_preprocess: #and acc_dep:
        return True


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
from distutils.dir_util import copy_tree

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project

sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import Forecast_preprocess as fpp
reload(fpp)

import Acc_Dep_ResSim_SacTrn
reload(Acc_Dep_ResSim_SacTrn)


# apparently, 2016 works with forecast, at least for Trinity, so we here make sure we are using the 2016
# "version" of the prescribed W2 models for forecast.
W2_models_for_input_copy = ['W2 Trinity Prescribed','W2 Lewiston Prescribed','W2 Keswick Prescribed']

def study_dir_from_run_dir(run_dir):
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def annual_config_dirs_from_run_dir(run_dir,model_name,startyear_str):
    study_dir = study_dir_from_run_dir(run_dir)
    model_dir = os.path.join(study_dir,'cequal-w2',model_name,model_name)  # don't know why this is two model_names deep!!
    annual_config_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,startyear_str)
    base_dir = os.path.join(study_dir,'shared','W2_annual_configs',model_name,'base')
    return model_dir,annual_config_dir,base_dir

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    startyear_str = '2016'
    
    run_dir = computeOptions.getRunDirectory()
    for W2_model in W2_models_for_input_copy:
        model_dir,annual_config_dir,base_dir = annual_config_dirs_from_run_dir(run_dir,W2_model,startyear_str)
        currentAlternative.addComputeMessage('model_dir: '+model_dir)
        currentAlternative.addComputeMessage('annual_config_dir: '+annual_config_dir)
        currentAlternative.addComputeMessage('base_dir: '+base_dir)
        if not os.path.exists(annual_config_dir):
            currentAlternative.addComputeMessage(W2_model+'- annual config not found; W2 may be configured incorrectly for this time window.')
        else:        
            # copy original W2 model alternative files to 'base' directory for safekeeping/later returning

            #if not os.path.exists(base_dir):
            #    os.makedir(base_dir)
            #base_files = os.listdir(base_dir)
            #if len(base_files) == 0:
            #    currentAlternative.addComputeMessage(W2_model+'- base files not found; copying from model folder')
            #    copy_tree(model_dir,base_dir)

            # remove all W2 model input files EXCEPT the .w2alt file
            for mfile in os.listdir(model_dir):
                if not mfile.endswith('.w2Alt') and not mfile.endswith('.w2Alt.bak'):
                    os.remove(os.path.join(model_dir,mfile))

            # copy over annual config input files 
            copy_tree(annual_config_dir,model_dir)
            currentAlternative.addComputeMessage('Copied W2 inputs file for '+startyear_str+' to '+W2_model+' model alternative folder')
            
            # now, W2 model alternative directory is configured for startyear, ready for the W2 plugin to
            # work it's magic on the W2_con file and execute simulation   

    data_preprocess = fpp.forecast_data_preprocess_W2_5Res(currentAlternative, computeOptions)

    # compute acc TODO
    #acc_dep = Acc_Dep_ResSim_SacTrn.computeAlternative(currentAlternative, computeOptions)

    if data_preprocess:# and acc_dep:
        return True

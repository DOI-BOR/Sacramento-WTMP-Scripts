
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

import tz_offset
reload(tz_offset)

model_output_and_target = [ 'model_output','target_temp']


def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    data_preprocess = fpp.forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions)

    # check if  W2 is a part of this simulation, and if so copy the target temperature to the linked record f-part
    # ------------------------------------------------------------------------------------------------------------------	
    scripting_run_dir = computeOptions.getRunDirectory()
    run_dir,_ = os.path.split(scripting_run_dir)
    _,run_name = os.path.split(run_dir)
    print('run_name: ',run_name)
    print('scripting_run_dir: ',scripting_run_dir)
    if 'w2' in run_name.lower():
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

        # fix potential W2 output timing inaccuracies in linked W2 output temperature
        # --------------------------------------------------------------------------------------------------
        # Note - if the offset time series to be corrected is not consistent in hourly offset minutes, this
        # code below could fail, and create a discontinuous record. Remove when interpolation is fixed
        # in the W2 Plugin.
        dss_file = computeOptions.getDssFilename()
        locations_obj = currentAlternative.getInputDataLocations()
        locations_obj_2 = DSS_Tools.organizeLocations(currentAlternative, locations_obj, model_output_and_target, return_dss_paths=False)
        dss_model_path, _ = DSS_Tools.getDataLocationDSSInfo(locations_obj_2[0], currentAlternative, computeOptions)
        dssFm = HecDss.open(dss_file)
        tsc = dssFm.get(dss_model_path,True) # use get with True here to capture entire record, 'read' seems to leave off data randomly
        data_check = True
        if tsc is None:
            data_check = False
        elif tsc.times is None:
            data_check = False
        elif len(tsc.times) == 0:
            data_check = False
        if not data_check:
            currentAlternative.addComputeErrorMessage('W2 outflow temperature data not found: '+dss_model_path+' \n')
            return False
        dtt = DSS_Tools.datetimes_from_tsc(tsc)
        
        # new method --------------------------------------------------------------------------------
        # snap first/last dates to XX:00 time, add one hour before and after time period, 
        # and interpolate hourly between those datetimes
        td_1hour = dt.timedelta(hours=1)
        dStart = dtt[0].replace(minute=0) - td_1hour
        dEnd = dtt[-1].replace(minute=0) + td_1hour*2
        td_interp = dEnd-dStart
        hours_interp = int(DSS_Tools._timedelta_to_seconds(dEnd-dStart) / 3600.0)
        dtt_interp = [dStart + td_1hour*i for i in range(hours_interp)]
        values_interp = DSS_Tools.linear_interp_datetime(dtt, tsc.values, dtt_interp)
        hectimes_interp = [HecTime(d.strftime('%d%b%Y %H%M')).value() for d in dtt_interp]

        tsc_shift = TimeSeriesContainer() # have to create a new tsc as the start time changes (don't know another way)
        tsc_shift.startTime = hectimes_interp[0]
        tsc_shift.times = hectimes_interp
        tsc_shift.fullName = dss_model_path #[:-1] + '--/'
        tsc_shift.units = tsc.units
        tsc_shift.type = tsc.type
        tsc_shift.values = values_interp
        tsc_shift.numberValues = len(hectimes_interp)        
        
        # old method --------------------------------------------------------------------------------
#        fix = False
#        td_tz = dt.timedelta(hours=tz_offset.hours)
#        hectimes = []
#        values = []  # collect values into a list for later use
#        for i,d in enumerate(dtt):
#            if d.minute != 0:
#                print('Found non-hourly W2 output minute: ',d.year,d.month,d.day,d.hour,d.minute,' fixing -> 0')
#                if d.minute > 29:
#                    d = d.replace(minute=0) + dt.timedelta(hours=1)
#                else:
#                    d = d.replace(minute=0)
#            #d += td_tz
#            hectimes.append(HecTime(d.strftime('%d%b%Y %H%M')).value()) #; print('3',d.strftime('%d%b%Y %H%M'))
#            print('time snap:',tsc.times[i],hectimes[i])
#            values.append(tsc.values[i]) 
#    
#        # also, copy start and end values one hour before and after current endpoints of record, to capture
#        # issues where W2 output does not cover the whole run time window
#        hec_timestep = hectimes[1] - hectimes[0]
#        startTime = hectimes[0] - hec_timestep
#        endTime = hectimes[-1] + hec_timestep
#        times = [startTime] + hectimes + [endTime]    ; print(len(times),len(tsc.times),len(tsc.values))
#        values = [values[0]] + values + [values[-1]] ; print(len(values))    
#        print(tsc.times[0],hectimes[0],startTime,times[0])
#        
#        tsc_shift = TimeSeriesContainer() # have to create a new tsc as the start time changes (don't know another way)
#        tsc_shift.startTime = startTime
#        tsc_shift.times = times
#        tsc_shift.fullName = W2OutflowTempRec
#        tsc_shift.units = tsc.units
#        tsc_shift.type = tsc.type
#        tsc_shift.values = values
#        tsc_shift.numberValues = len(times)
    
        dssFm.put(tsc_shift)    
        dssFm.close()

    if data_preprocess: #and acc_dep:
        return True

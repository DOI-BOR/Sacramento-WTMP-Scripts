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
import os, sys, csv

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project
#import hec.hecmath.TimeSeriesMath as tsmath
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import DSS_Tools
reload(DSS_Tools)

import DMS_preprocess
reload(DMS_preprocess)

import equilibrium_temp
reload(equilibrium_temp)

import create_balance_flow_jython as cbfj
reload(cbfj)

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

def interp(x, xp, fp, left=None, right=None):
    """
    One-dimensional linear interpolation.

    Returns the one-dimensional piecewise linear interpolant to a function
    with given values at discrete data-points.

    Parameters
    ----------
    x : array_like
        The x-coordinates at which to evaluate the interpolated values.
    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.
    fp : 1-D sequence of floats
        The y-coordinates of the data points, same length as `xp`.
    left : float, optional
        Value to return for `x < xp[0]`, default is `fp[0]`.
    right : float, optional
        Value to return for `x > xp[-1]`, default is `fp[-1]`.

    Returns
    -------
    y : float or ndarray
        The interpolated values, same shape as `x`.
    """

    if isinstance(x, list):
        return [interp(point, xp, fp, left, right) for point in x]
    else:
        if left is None:
            left = fp[0]
        if right is None:
            right = fp[-1]

        if x < xp[0]:
            return left
        elif x > xp[-1]:
            return right
        else:
            for i in range(len(xp) - 1):
                if x >= xp[i] and x <= xp[i+1]:
                    # Perform the linear interpolation
                    t = (x - xp[i]) / (xp[i+1] - xp[i])
                    return fp[i] + t * (fp[i+1] - fp[i])

def eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out):
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()

	# get at data and times in the formats needed
    dssFm = HecDss.open(at[0])        
    tsc = dssFm.read(at[1], starttime_str, endtime_str, False).getData()
    tsc_int_times = tsc.times
    dtt = DSS_Tools.hectime_to_datetime(tsc)
    at_data = tsc.values
    dssFm.close()
	# get the rest of the data over the same period
    cl_data = DSS_Tools.data_from_dss(cl[0],cl[1],starttime_str,endtime_str)
    ws_data = DSS_Tools.data_from_dss(ws[0],ws[1],starttime_str,endtime_str)
    sr_data = DSS_Tools.data_from_dss(sr[0],sr[1],starttime_str,endtime_str)
    td_data = DSS_Tools.data_from_dss(td[0],td[1],starttime_str,endtime_str)
    
	# calc_equilibrium_temp(dtt, at, cl, sr, td, ws)
    Te = equilibrium_temp.calc_equilibrium_temp(dtt,at_data,cl_data,sr_data,td_data,ws_data)
    
    print('writing: ',eq_temp_out[1])
    tsc = TimeSeriesContainer()
    tsc.times = tsc_int_times
    tsc.fullName = eq_temp_out[1]
    tsc.values = Te
    tsc.units = 'C'
    tsc.type = 'INST-VAL'
    tsc.numberValues = len(tsc.values)

    tsm = tsmath(tsc)
    tsm_day = DSS_Tools.standardize_interval(tsm,'1day')
    tsm_wk = DSS_Tools.standardize_interval(tsm,'1week')
        
    dssFmOut = HecDss.open(eq_temp_out[0])
    dssFmOut.write(tsc)
    dssFmOut.write(tsm_day)
    dssFmOut.write(tsm_wk)
    dssFmOut.close()


def storage_to_elev(res_name,elev_stor_area,forecast_dss,storage_rec,conic=False):
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True) # read ALL data in record

    elev = []
    if conic:
        print('Conic interpolation of elevations from storage not supported yet.')
        sys.exit(-1)
    else:
        for j in range(tsc.numberValues):
            elev.append(cbfj.linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], tsc.values[j]))
            print('stor2elev: ',j,tsc.times[j],tsc.values[j])

    #if tsc.type.lower != 'inst-val':
    #    # going from PER-CUM -> INST-VAL, so need to move the times up one to account for the end-of-period reporting of PER-CUM
    #    time_delta = tsc.times[1] - tsc.times[0]
    #    times = [tsc.times[0]-time_delta]
    #    for j in range(tsc.numberValues-1):
    #        times.append(tsc.times[j])
    #    tsc.times = times

    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = elev
    dssFmRec.write(tsc)
    dssFmRec.close()

def invent_elevation(res_name,forecast_dss,storage_rec,elev_constant_ft):
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True)
    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = [elev_constant_ft for j in range(tsc.numberValues)]

    #if tsc.type.lower != 'inst-val':
    #    # going from PER-CUM -> INST-VAL, so need to move the times up one to account for the end-of-period reporting of PER-CUM
    #    time_delta = tsc.times[1] - tsc.times[0]
    #    times = [tsc.times[0]-time_delta]
    #    for j in range(tsc.numberValues-1):
    #        times.append(tsc.times[j])
    #    tsc.times = times
    
    dssFmRec.write(tsc)
    dssFmRec.close()

def write_forecast_elevations():

	# write an hourly forecast elevation based on starting elevation and flows
    DSS_Tools.resample_dss_ts(forecast_dss,'//FOLSOM/PUMPING (FP)//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'//Folsom-NF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'//Folsom-SF-in/FLOW-IN//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    DSS_Tools.resample_dss_ts(forecast_dss,'/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Day/AMER_BC_SCRIPT/',rtw,forecast_dss,'1HOUR')
    inflow_records = ['//Folsom-NF-in/FLOW-IN//1Hour/AMER_BC_SCRIPT/',
                      '//Folsom-SF-in/FLOW-IN//1Hour/AMER_BC_SCRIPT/',
                      '/AMERICAN RIVER/FOLSOM LAKE/FLOW-ACC-DEP//1Hour/AMER_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['//FOLSOM/FLOW-RELEASE//1Hour/AMER_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'//Folsom/ELEV//1Month/AMER_BC_SCRIPT/')
    print('starting_elevation ',starting_elevation)
    cbfj.predict_elevation(currentAlternative, rtw, 'Folsom', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Folsom/ELEV-FORECAST//1HOUR/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None)



def study_dir_from_run_dir(run_dir):
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def model_dir_from_run_dir(run_dir,model_place,model_name):
    study_dir = study_dir_from_run_dir(run_dir)
    model_dir = os.path.join(study_dir,'cequal-w2',model_place,model_name)
    return model_dir

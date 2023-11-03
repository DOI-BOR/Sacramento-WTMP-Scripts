"""
Created on 4/28/2022

@author: Stephen Andrews, Scott Burdick-Yahya
@organization: Resource Management Associates
@contact: steve@rmanet.com
@note:
modified for jython to be used in WAT by SBY on 8/3/2023
"""

# import datetime as dt
# import numpy as np
# from pydsstools.heclib.dss import HecDss
# from pydsstools.core import TimeSeriesContainer
# from scipy import     
# import pandas as pd

import math
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
import os
# import DSS_Tools
# reload(DSS_Tools)


from hec.io import DSSIdentifier
from hec.heclib.util import HecTime
from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

def linear_interpolation(x_values, y_values, x):
    if len(x_values) != len(y_values) or len(x_values) < 2:
        raise ValueError("Input lists must have the same length and contain at least 2 data points.")

    for i in range(1, len(x_values)):
        if x <= x_values[i]:
            x0, y0 = x_values[i - 1], y_values[i - 1]
            x1, y1 = x_values[i], y_values[i]

            # Perform linear interpolation
            y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)

            return y

    # If x is beyond the range of x_values, raise an error
    raise ValueError("Interpolation point is outside the range of provided data.")

def read_elev_storage_area_file(file_name, res_name):
    # These are in [elev, stor, area] with units [ft, acre-ft, acre]
    elevstorarea = {} #avoid lists doing weird things like mixing up order..
    elev = []
    stor = []
    area = []
    import os
    print('cwd: ' + os.getcwd())
    if res_name.lower() == 'natoma':
        with open(file_name, 'r') as fn:
            for line in fn:
                sline = line.strip().split(',')
                elev.append(float(sline[0]))
                area.append(float(sline[1]))
    else:
        with open(file_name, 'r') as fn:
            for line in fn:
                sline = line.strip().split(',')
                elev.append(float(sline[0]))
                stor.append(float(sline[1]))
                area.append(float(sline[2]))
    elevstorarea['elev'] = elev
    elevstorarea['stor'] = stor
    elevstorarea['area'] = area
    return elevstorarea

def build_conic_storage_array(elev, area, firstStorageValue=0.0):
    '''Find storage of slabs between measurement points on the elevation area curve,
    using a conic estimation.  Adapted from storage.java from HEC ResSim, 2022-06-17'''
    # calculate storage at each elevation using conic formula
    n_measures = len(elev)
    storage = []
    storage.append(firstStorageValue)
    for i in range(1, n_measures):
        h = elev[i] - elev[i-1]
        storage.append(h/3. * (area[i-1] + area[i] + math.sqrt(area[i-1] * area[i])) + storage[i-1])
    return storage


def conic_storage_interp(interpElev, elev, area, conicStorage, idx):
    '''Find storage between measurement points on the elevation area curve,
    using interpolation between conic layers.  Adapted from storage.java from
    HEC ResSim, 2022-06-17'''
    h = (interpElev - elev[idx]) / (elev[idx+1] - elev[idx])
    geomMean = math.sqrt(area[idx] * area[idx+1])
    areaInterp = area[idx] + 2.*(geomMean - area[idx])*h + (area[idx] + area[idx+1] - 2.*geomMean)*h*h
    storageInterp = (interpElev - elev[idx])/3. * (area[idx] + areaInterp + math.sqrt(area[idx] * areaInterp)) + conicStorage[idx]
    return storageInterp


def get_elev_layer_idx(elev, obs_elev, elev_stor_area):
    # find lower bounding index of where elevation lands in elev-stor-area table
    # idx = np.argmin(np.abs(elev-obs_elev))
    # if elev_stor_area[idx, 0] > obs_elev:
    #     idx -= 1
    # return idx

    idx = UNDEFINED_DOUBLE
    min_val = None
    for i in range(len(elev)):
        valchk = abs(elev[i]-obs_elev) #TODO: is this multidimensional?
        if math.isnan(valchk):
            min_val = valchk
            idx = i
        elif min_val == None:
            min_val = valchk
            idx = i
        elif valchk < min_val:
            min_val = valchk
            idx = i
    if idx != UNDEFINED_DOUBLE:
        if elev_stor_area['elev'][idx] > obs_elev: #TODO: is this multidimensional?
            idx -= 1
    else:
        idx = -1
    return idx

def get_balance_period(balance_period):
    if 'hour' in balance_period.lower():
        return float(balance_period.lower().replace('hour', ''))
    elif 'day' in balance_period.lower():
        return float(balance_period.lower().replace('day', '')) * 24
    elif 'min' in balance_period.lower():
        return float(balance_period.lower().replace('min', '')) / 60

def check_dss_intervals(records, balance_period, currentAlt):
    for r in records:
        if balance_period.lower() not in r.lower():
            currentAlt.addComputeMessage('DSS record {0} not matching time interval {1}'.format(r, balance_period))
            sys.exit(-1)


def read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str):
    '''pathname may contain the dss filepath additionally before the dss ts path, separated by '::'
       If so, use that dss file.'''
    if '::' in pathname:
        print('Splitting and reading:',pathname)
        alt_dss_file,pathname_clean = pathname.split('::')
        dssFmRec = HecDss.open(alt_dss_file)
        tsc = dssFmRec.read(pathname_clean, starttime_str, endtime_str, False).getData()
        dssFmRec.close()
    else:
        tsc = dssFm.read(pathname, starttime_str, endtime_str, False).getData()
    return tsc


def create_balance_flows(currentAlt, timewindow, res_name, inflow_records, outflow_records, stage_record, evap_record,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         storage_dss_record_name='', evap_dss_record_name='',
                         balance_period_str="1HOUR", use_conic=False, write_evap=False, write_storage=False,
                         alt_period=None,alt_period_string=None, lookback_padding=1440):


    check_dss_intervals(inflow_records, balance_period_str, currentAlt)
    check_dss_intervals(outflow_records, balance_period_str, currentAlt)
    check_dss_intervals([stage_record, evap_record], balance_period_str, currentAlt)
    
    balance_period = get_balance_period(balance_period_str) # convert to (float) hours
    print('balance_period ' + str(balance_period))
    
    cfs_2_acreft = balance_period * 3600. / 43559.9
    acreft_2_cfs = 1. / cfs_2_acreft

    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #01Jan2014 0000

    # add lookback padding to enable ResSim to have balance flows on 1st timestep
    # starttime_hectime_obj = HecTime(starttime_str).add(lookback_padding)
    # starttime_str = starttime_hectime_obj.date()


    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dss_file)

    inflows = []
    outflows = []
    times = []

    # Read inflows
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        print('\nreading: ' + str(pathname))
        try:
       
            print(starttime_str, endtime_str)
            print(dss_file)
            tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
            values = tsc.values
            hectimes = tsc.times
            units = tsc.units
            # print('num values {0}'.format(len(values)))
            # print('start {0}'.format(ts_data.getStartTime()))
            # print('end {0}'.format(ts_data.getEndTime()))
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        if len(inflows) == 0:
            inflows = values
            times = hectimes #TODO: check how this handles missing values
        else:
            for vi, v in enumerate(values):
                inflows[vi] += v

    
    
    # Read outflows
    print('Reading outflow records')
    for j, outflow_record in enumerate(outflow_records):  # for each of the dss paths in inflow_records    
        pathname = outflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
        try:
            values = tsc.values
            hectimes = tsc.times
            units = tsc.units
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
          
        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        if len(outflows) == 0:
            outflows = values
        else:
            for vi, v in enumerate(values):
                outflows[vi] += v


    # Inflow minus outflow record
    inflow_outflow = []
    for i in range(len(inflows[1:])):
        inflow_outflow.append(inflows[i+1] - outflows[i+1])
   # this is in cfs (period avg vals)

    # Read stage
    print('Reading stage')
    tsc = read_ts_rec_w_optional_fname(dssFm, stage_record, starttime_str, endtime_str)
    try:
        stage = tsc.values
        hectimes = tsc.times
        if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
            print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
            st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
            stage = stage[st_offset:]
            hectimes = hectimes[st_offset:]
        if hectimes[-1] > endtime_hectime:
            print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
            st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
            stage = stage[:(len(hectimes) - st_offset)]
            hectimes = hectimes[:(len(hectimes) - st_offset)]
        print('Number Stage Values: {0}'.format(len(stage)))
        
    except HecMathException:
        currentAlt.addComputeMessage('ERROR reading' + str(stage_record))
        sys.exit(-1)

    # Read evap
    print('Reading evap')
    tsc = read_ts_rec_w_optional_fname(dssFm, evap_record, starttime_str, endtime_str)
    try:
        evap = tsc.values
        hectimes = tsc.times
        if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
            print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
            st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
            evap = evap[st_offset:]
            hectimes = hectimes[st_offset:]
        if hectimes[-1] > endtime_hectime:
            print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
            st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
            evap = evap[:(len(hectimes) - st_offset)]
            hectimes = hectimes[:(len(hectimes) - st_offset)]
        print('Number Evap Values: {0}'.format(len(evap)))
    except HecMathException:
        currentAlt.addComputeMessage('ERROR reading' + str(evap_record))
        sys.exit(-1)

    # Build conic storage array for interpolation later
    conic_storage = build_conic_storage_array(elev_stor_area['elev'], elev_stor_area['area'])

    # Calculations
    n = len(stage) - 1
    flow_resid = []
    flow_evap = []
    # area_fnct =     .interp1d(elev_stor_area[:, 0], elev_stor_area[:, 2])
    # area_fnct = linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'])

    storage_record = []

    # if not use_conic:
        # stor_fnct =     .interp1d(elev_stor_area[:, 0], elev_stor_area[:, 1])
        # stor_fnct = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'])

    for k in range(n):
        stage_start = stage[k]
        stage_end = stage[k+1]

        if use_conic:
            idx1 = get_elev_layer_idx(elev_stor_area['elev'], stage_start, elev_stor_area)
            storage_start = conic_storage_interp(stage_start, elev_stor_area['elev'], elev_stor_area['area'], conic_storage, idx1)
            idx2 = get_elev_layer_idx(elev_stor_area['elev'], stage_end, elev_stor_area)
            storage_end = conic_storage_interp(stage_end, elev_stor_area['elev'], elev_stor_area['area'], conic_storage, idx2)
        else:
            storage_start = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], stage_start)
            storage_end = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], stage_end)

        delta_stor_from_stage = storage_end - storage_start  # in acre-ft
        delta_stor_flow = delta_stor_from_stage * acreft_2_cfs # in cfs
        inflow_minus_outflow = inflow_outflow[k]  # in cfs
        # area_avg = 0.5 * (area_fnct(stage_start) + area_fnct(stage_end))
        area_avg = 0.5 * (linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_start) +
                          linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_end))
        evap_flow_loss = (evap[k] * area_avg) * acreft_2_cfs  # in cfs

        resid = delta_stor_flow - (inflow_minus_outflow - evap_flow_loss)
        flow_resid.append(resid)
        flow_evap.append(evap_flow_loss)
        storage_record.append(storage_start)


    if True:
        print('Writing to CSV')
        # dump to CSV if DSS is mis-behaving or if flows are -999. etc.
        with open(os.path.join(shared_dir, "{0}_balance_flow.csv".format(res_name)), 'w') as opf:
            opf.write('date, balance_flow [cfs]\n')
            for i in range(len(flow_resid)):
                new_line = ','.join([str(times[i]), str(flow_resid[i]), '\n'])
                opf.write(new_line)
        # pd.DataFrame({'date':pd.to_datetime([tstart + model_time_step*i for i in range(len(flow_resid))]),'balance_flow [cfs]':flow_resid}).to_csv("%s balance flow.csv"%res_name)

    dssFm_out = HecDss.open(output_dss_file)
    
    # Output record
    # copy back 1st balance flow record 2 steps, instead of writing from 1st valid balance calc.
    # otherwise, time-averaging the balanece flows later leaves off the 1st time step needed for a ResSim run
    # best we can do I guess to make ResSim computes work

    steptime = times[1]-times[0]
    tsc = TimeSeriesContainer()
    #tsc.times = times[1:]
    tsc.startTime = times[0] - steptime
    tsc.interval = int(balance_period)*60
    tsc.fullName = output_dss_record_name
    #tsc.values = flow_resid
    tsc.values = [flow_resid[0],flow_resid[0]] + flow_resid
    #tsc.startTime = times[1]
    tsc.units = 'CFS'
    tsc.type = 'PER-AVER'
    #tsc.endTime = times[-1]
    # tsc.startHecTime = timewindow.getStartTime()
    # tsc.endHecTime = timewindow.getEndTime()
    tsc.numberValues = len(tsc.values)
    dssFm_out.write(tsc)

    if alt_period is not None:
        if alt_period_string.lower() != balance_period_str.lower():
            tsm = dssFm_out.read(output_dss_record_name)
            tsm_new_interval = tsm.transformTimeSeries(alt_period_string, "", "AVE")
            dssFm_out.write(tsm_new_interval)

    if write_evap:
        tsc = TimeSeriesContainer()
        tsc.times = times
        tsc.fullName = evap_dss_record_name
        tsc.values = flow_evap
        tsc.startTime = times[1]
        tsc.units = 'CFS'
        tsc.type = 'PER-AVER'
        tsc.endTime = times[-1]
        tsc.numberValues = len(flow_resid)
        tsc.startHecTime = timewindow.getStartTime()
        tsc.endHecTime = timewindow.getEndTime()
        dssFm_out.write(tsc)

    if write_storage:
        tsc = TimeSeriesContainer()
        tsc.times = times
        tsc.fullName = storage_dss_record_name
        tsc.values = storage_record
        tsc.startTime = times[1]
        tsc.units = "AC-FT"
        tsc.type = 'INST-VAL'  # is this right?
        tsc.endTime = times[-1]
        tsc.numberValues = len(flow_resid)
        tsc.startHecTime = timewindow.getStartTime()
        tsc.endHecTime = timewindow.getEndTime()
        dssFm_out.write(tsc)

    dssFm.close()
    dssFm_out.close()
    return True

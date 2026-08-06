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
    """
    Perform simple linear interpolation to estimate a y-value at a
    given x, based on a sorted list of known (x, y) points.

    Args:
        x_values (list of float): Known x-coordinates, assumed to
            be sorted in ascending order.
        y_values (list of float): Known y-coordinates, parallel to
            x_values.
        x (float): The x-value to interpolate a y-value for. Must
            fall within the range of x_values.

    Returns:
        float: The linearly interpolated y-value at x.

    Raises:
        ValueError: If x_values and y_values differ in length or
            contain fewer than 2 points, or if x falls outside the
            range covered by x_values.
    """
    if len(x_values) != len(y_values) or len(x_values) < 2:
        raise ValueError("Input lists must have the same length and contain at least 2 data points.")

    # search forward for the first known point at or beyond x, then interpolate
    # within the segment ending at that point
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
    """
    Read a CSV file containing an elevation-storage-area curve for
    a reservoir, and return it as a dictionary of parallel lists.

    NOTE: The 'natoma' reservoir uses a special-cased two-column
    file format (elevation, area only - no storage column), while
    all other reservoirs are expected to have a three-column format
    (elevation, storage, area). Calling this with res_name='natoma'
    against a three-column file (or vice versa) will misread the
    data.

    Args:
        file_name (str): Path to the CSV file to read. Expected to
            contain one row per elevation point, comma-separated.
        res_name (str): Reservoir name. If (case-insensitive) equal
            to 'natoma', the file is parsed as two columns
            (elev, area). Otherwise, it is parsed as three columns
            (elev, stor, area).

    Returns:
        dict: A dictionary with keys 'elev', 'stor', and 'area',
            each mapping to a list of floats parsed from the file,
            in units of feet, acre-feet, and acres respectively.
            For 'natoma', the 'stor' list will be empty, since that
            file format does not include a storage column.
    """
    # These are in [elev, stor, area] with units [ft, acre-ft, acre]
    elevstorarea = {} #avoid lists doing weird things like mixing up order..
    elev = []
    stor = []
    area = []
    import os
    print('cwd: ' + os.getcwd())
    # Natoma's file uses a different (2-column) format than other reservoirs
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
    # for each elevation layer, add the conic-estimated storage volume of that slab
    # onto the cumulative total from the layer below
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
    """
    Find the index of the elevation layer that lower-bounds a given
    observed elevation, within an elevation-storage-area table.

    Args:
        elev (list of float): Known elevation values to search
            within.
        obs_elev (float): The observed elevation to locate within
            elev.
        elev_stor_area (dict): The elevation-storage-area table
            (as produced by read_elev_storage_area_file), used to
            re-check whether the closest elevation is above or
            below obs_elev.

    Returns:
        int: The index of the elevation value that lower-bounds
            obs_elev (i.e. elev[idx] <= obs_elev). Returns -1 if no
            valid index could be determined (e.g. an empty elev
            list).
    """
    # find lower bounding index of where elevation lands in elev-stor-area table
    # idx = np.argmin(np.abs(elev-obs_elev))
    # if elev_stor_area[idx, 0] > obs_elev:
    #     idx -= 1
    # return idx

    idx = UNDEFINED_DOUBLE
    min_val = None
    # scan every elevation value to find the one closest to obs_elev
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
    # if a closest match was found, check whether it overshoots obs_elev and step back if so
    if idx != UNDEFINED_DOUBLE:
        if elev_stor_area['elev'][idx] > obs_elev: #TODO: is this multidimensional?
            idx -= 1
    else:
        idx = -1
    return idx

def get_balance_period(balance_period):
    """
    Convert a HEC-style balance period string into a number of
    hours.
    """
    # convert the period string to hours based on which unit substring it contains
    if 'hour' in balance_period.lower():
        return float(balance_period.lower().replace('hour', ''))
    elif 'day' in balance_period.lower():
        return float(balance_period.lower().replace('day', '')) * 24
    elif 'min' in balance_period.lower():
        return float(balance_period.lower().replace('min', '')) / 60

def check_dss_intervals(records, balance_period, currentAlt):
    """
    Validate that every DSS record path in a list mentions the
    expected balance period, logging an error and exiting if any do
    not.

    Args:
        records (list of str): DSS record paths to check.
        balance_period (str): The expected period string (e.g.
            '1Day') that should appear somewhere in each record
            path (typically the E-part).
        currentAlt: The alternative object being computed. Used to
            log an error message if a mismatch is found.

    Returns:
        None. Exits the process (sys.exit(-1)) if any record does
        not contain balance_period.
    """
    # check every record for the expected time interval substring
    for r in records:
        if balance_period.lower() not in r.lower():
            currentAlt.addComputeMessage('DSS record {0} not matching time interval {1}'.format(r, balance_period))
            sys.exit(-1)


def read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str):
    '''pathname may contain the dss filepath additionally before the dss ts path, separated by '::'
       If so, use that dss file.'''
    # if an alternate DSS file is embedded in the path string, open and read from that file instead
    if '::' in pathname:
        print('Splitting and reading:',pathname)
        alt_dss_file,pathname_clean = pathname.split('::')
        dssFmRec = HecDss.open(alt_dss_file)
        tsc = dssFmRec.read(pathname_clean, starttime_str, endtime_str, False).getData()
        dssFmRec.close()
    else:
        tsc = dssFm.read(pathname, starttime_str, endtime_str, False).getData()
    return tsc


def read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, starttime_str, endtime_str,
                          starttime_hectime, endtime_hectime):

     """
    Read and sum multiple inflow and outflow DSS records over a
    given time window, converting units to cfs as needed, and
    return the net inflow-minus-outflow time series.
    """
     dssFm = HecDss.open(dss_file)

    inflows = []
    outflows = []
    times = []

    # --- Read and sum all inflow records ---
    # Read inflows
    print('Reading inflows')
    # for each inflow record, read it, trim it to the time window, convert units if
    # needed, and add it to the running inflow total
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
            # if the record's data extends past the run window, trim the trailing values
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # if this record is in cms, convert it to cfs before summing
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting inflow to cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        if len(inflows) == 0:
            inflows = values
            times = hectimes #TODO: check how this handles missing values
        else:
            # add this record's values onto the running inflow total, timestep by timestep
            for vi, v in enumerate(values):
                inflows[vi] += v

    
    
    # Read outflows
    print('Reading outflow records')
    # for each outflow record, read it, trim it to the time window, convert units if
    # needed, and add it to the running outflow total
    for j, outflow_record in enumerate(outflow_records):  # for each of the dss paths in inflow_records    
        pathname = outflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
        try:
            values = tsc.values
            hectimes = tsc.times
            units = tsc.units
            # if the record's data starts before the run window, trim the leading values
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            # if the record's data extends past the run window, trim the trailing values
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
          
        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # if this record is in cms, convert it to cfs before summing
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting outflow cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # the first record seeds the running total; subsequent records are added on
        if len(outflows) == 0:
            outflows = values
        else:
            for vi, v in enumerate(values):
                outflows[vi] += v

    dssFm.close()

    # Inflow minus outflow record
    inflow_outflow = []
    for i in range(len(inflows[1:])):
        inflow_outflow.append(inflows[i+1] - outflows[i+1])
   # this is in cfs (period avg vals)

    currentAlt.addComputeMessage("Len inflow_outflow:"+str(len(inflow_outflow)))
    currentAlt.addComputeMessage("Len times:"+str(len(times)))

    return times,inflow_outflow


def predict_elevation(currentAlt, starttime_str, endtime_str, res_name, inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Hour'):
    '''From inflows/outflows, predict hourly elevation, useful for lookback/starting elevation for forecasts starting
    on arbitrary dates during forecast period
    '''
    balance_period = get_balance_period(balance_period_str) # convert to (float) hours
    
    check_dss_intervals(inflow_records, balance_period_str, currentAlt)
    check_dss_intervals(outflow_records, balance_period_str, currentAlt)
       
    cfs_2_acreft = balance_period * 3600. / 43559.9
    acreft_2_cfs = 1. / cfs_2_acreft

    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()

    times,inflow_outflow = read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, 
                                                 starttime_str, endtime_str, starttime_hectime, endtime_hectime)

    currentAlt.addComputeMessage("Len inflow_outflow:"+str(len(inflow_outflow)))
    currentAlt.addComputeMessage("Len times:"+str(len(times)))
    
    # TODO: support conic interpolation
    # TODO: support evap, but really that's just a positive outflow....
    storage = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], starting_elevation)
    storage = [storage,]
    elev_predicted = []
    for i in range(len(inflow_outflow)):
        storage.append( storage[-1] + inflow_outflow[i]*cfs_2_acreft )
        elev_predicted.append( linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], storage[-1]) )

    # Output record
    dssFm_out = HecDss.open(output_dss_file)
    steptime = times[1]-times[0]
    tsc = TimeSeriesContainer()
    #tsc.times = times[1:]
    tsc.startTime = times[0] #- steptime
    tsc.interval = int(balance_period)*60
    tsc.fullName = output_dss_record_name
    tsc.values = [starting_elevation] + elev_predicted
    #tsc.startTime = times[1]
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.numberValues = len(tsc.values)
    dssFm_out.write(tsc)

    recparts = output_dss_record_name.split('/')
    recparts[3] = 'STORAGE-PREDICTED'
    tsc.startTime = times[0]
    tsc.fullName = '/'.join(recparts)
    tsc.values = storage
    tsc.units = 'ac-ft'
    tsc.type = 'PER-CUM'
    dssFm_out.write(tsc)

    if alt_period is not None:
        if alt_period_string.lower() != balance_period_str.lower():
            tsm = dssFm_out.read(output_dss_record_name)
            tsm_new_interval = tsm.transformTimeSeries(alt_period_string, "", "AVE")
            dssFm_out.write(tsm_new_interval)

    dssFm_out.close()


def create_balance_flows(currentAlt, timewindow, res_name, inflow_records, outflow_records, stage_record, evap_record,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         storage_dss_record_name='', evap_dss_record_name='',
                         balance_period_str="1HOUR", use_conic=False, write_evap=False, write_storage=False,
                         alt_period=None,alt_period_string=None, lookback_padding=1440):


    """
    Compute a reservoir mass-balance ("balance flow") time series -
    the residual flow needed to reconcile observed storage change
    with recorded inflows, outflows, and evaporation - and write it
    to DSS, with optional companion evaporation and storage records
    and an optional resampled copy at an alternate time interval.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        timewindow: The run time window object, used to get the
            start and end time strings for reading input data.
        res_name (str): Reservoir name, used to label the debug CSV
            output file.
        inflow_records (list of str): DSS paths of records to sum
            as total inflow.
        outflow_records (list of str): DSS paths of records to sum
            as total outflow.
        stage_record (str): DSS path of the reservoir stage
            (elevation) record.
        evap_record (str): DSS path of the evaporation record.
        elev_stor_area (dict): Elevation-storage-area lookup table
            (as produced by read_elev_storage_area_file), with keys
            'elev', 'stor', and 'area'.
        dss_file (str): Path to the DSS file containing the input
            records.
        output_dss_record_name (str): DSS path to write the
            computed balance flow record to.
        output_dss_file (str): Path to the DSS file to write output
            records to.
        shared_dir (str): Directory to write the debug CSV file to.
        storage_dss_record_name (str, optional): DSS path to write
            the computed storage record to, if write_storage is
            True. Defaults to ''.
        evap_dss_record_name (str, optional): DSS path to write the
            computed evaporation-flow-loss record to, if write_evap
            is True. Defaults to ''.
        balance_period_str (str, optional): The expected time
            interval of all input records (e.g. '1HOUR'). Defaults
            to "1HOUR".
        use_conic (bool, optional): If True, use conic frustum
            interpolation (via get_elev_layer_idx and
            conic_storage_interp) to estimate storage from stage.
            If False, use simple linear interpolation. Defaults to
            False.
        write_evap (bool, optional): If True, also write the
            computed evaporation-flow-loss series to DSS. Defaults
            to False.
        write_storage (bool, optional): If True, also write the
            computed storage series to DSS. Defaults to False.
        alt_period (optional): If not None, triggers writing an
            additional resampled copy of the balance flow at
            alt_period_string's interval.
        alt_period_string (str, optional): The alternate time
            interval string to resample to, used only if alt_period
            is not None.
        lookback_padding (int, optional): Currently unused within
            the function body (referenced only in a commented-out
            block). Defaults to 1440.

    Returns:
        bool: True once the balance flow (and any requested
            companion records) have been written to DSS.
    """
    # --- Validate that all input records share the expected time interval ---
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

    # --- Read and sum inflow/outflow, then read stage and evaporation ---
    times,inflow_outflow = read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, 
                                                 starttime_str, endtime_str, starttime_hectime, endtime_hectime)
    print('len times:',len(times))
    print('len inflow_outflow:',len(inflow_outflow))

    dssFm = HecDss.open(dss_file)

    # Read stage
    print('Reading stage')
    tsc = read_ts_rec_w_optional_fname(dssFm, stage_record, starttime_str, endtime_str)
    try:
        stage = tsc.values
        hectimes = tsc.times
        # if the record's data starts before the run window, trim the leading values
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

        # if stage is recorded in meters, convert it to feet
        if tsc.units.lower() == 'm':
            currentAlt.addComputeMessage('Converting stage m to ft')
            print('Converting stage cms to cfs')
            convvals = []
            for elev in stage:
                convvals.append(elev * 3.280839895)
            stage = convvals
        
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

    # --- Compute the per-timestep mass balance residual (the balance flow) ---
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

    # for each timestep, compute storage change, net inflow/outflow, evaporation loss,
    # and the residual flow needed to balance them all
    for k in range(n):
        stage_start = stage[k]
        stage_end = stage[k+1]

        # use conic frustum interpolation if requested, otherwise simple linear interpolation
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


    # --- Write the balance flow series to a debug CSV file ---
    if True:
        print('Writing to CSV')
        # dump to CSV if DSS is mis-behaving or if flows are -999. etc.
        with open(os.path.join(shared_dir, "{0}_balance_flow.csv".format(res_name)), 'w') as opf:
            opf.write('date, balance_flow [cfs]\n')
            # write one CSV row per computed balance flow value
            for i in range(len(flow_resid)):
                new_line = ','.join([str(times[i]), str(flow_resid[i]), '\n'])
                opf.write(new_line)
        # pd.DataFrame({'date':pd.to_datetime([tstart + model_time_step*i for i in range(len(flow_resid))]),'balance_flow [cfs]':flow_resid}).to_csv("%s balance flow.csv"%res_name)

    dssFm_out = HecDss.open(output_dss_file)
    
    # Output record

    # sometimes ResSim does not include the start record in period average simulations, so if one flow or elevation data
    # record is missing, the calc can sometimes go way off.  Constrain to realistic values, set invalid to zero.
    # also, recs offset by timezone and/or daily records not at midnight can introduced bad values on the first or last days.
    # So, filter at least the first 24 hours, last 24 hours
    check_steps = 1
    bad_flow_bound = 1.e7
    if balance_period_str.lower() == '1hour':
        check_steps = 24
    for i in range(check_steps):
        for idx in [i,-1-i]:
            if math.isnan(flow_resid[idx]) or flow_resid[idx] > bad_flow_bound or flow_resid[idx] < -bad_flow_bound:
                flow_resid[idx] = 0.0

    steptime = times[1]-times[0]
    tsc = TimeSeriesContainer()
    #tsc.times = times[1:]
    tsc.startTime = times[0] - steptime
    tsc.interval = int(balance_period)*60
    tsc.fullName = output_dss_record_name
    #tsc.values = flow_resid

    # copy back 1st balance flow record 2 steps, instead of writing from 1st valid balance calc.
    # otherwise, time-averaging the balanece flows later leaves off the 1st time step needed for a ResSim run
    # best we can do I guess to make ResSim computes work
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

    # --- Optionally write companion evaporation and storage records ---
    # if requested, write the computed evaporation-flow-loss series to DSS
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

    # if requested, write the computed storage series to DSS
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

#version 2.0
#modified 03-28-2023 by Scott Burdick-Yahya

from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath
from com.rma.model import Project
import os,shutil,copy
from java.util import Vector, Date

from datetime import datetime
from hec.heclib.util import HecTime  # Replace with actual import


def standardize_interval(tsm, interval, makePerAver=True):
    tsc = tsm.getData()
    if interval.lower()=='1hour':
        intint = 60
    elif interval.lower()=='1day':
        intint=1440
    elif interval.lower()=='1week':
        intint=10800
    else:
        print('interval not supported:',interval)
        sys.exit(-1)

    if tsc.interval != intint:
        if makePerAver:
            #tsc.type = 'PER-AVER'  # make sure it's per-aver ... we are
            tsm.setType('PER-AVER')
        return tsm.transformTimeSeries(interval, "", "AVE")
    else:
        return tsm


def data_from_dss(dss_file,dss_rec,starttime_str, endtime_str,):
    dssFm = HecDss.open(dss_file)        
    tsc = dssFm.read(dss_rec, starttime_str, endtime_str, False).getData()
    dssFm.close()
    return tsc.values
 

def hectime_to_datetime(tsc):

    dtt = []
    for j in range(tsc.numberValues):
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        java_date = tsc.getHecTime(j).getJavaDate(0)  
	    
        # Convert Java Date to Python datetime
        timestamp = (java_date.getTime() / 1000)
        dtt.append(datetime.fromtimestamp(timestamp))

    return dtt

def fixInputLocationFpart(currentAlternative, tspath):
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])
    tspath = tspath.split('/')
    fpart = tspath[6]
    fpart_split = fpart.split(':')
    new_fpart = new_fpart_start + ':' + fpart_split[-1]
    tspath[6] = new_fpart
    tspath = '/'.join(tspath)
    return tspath

def appendAPart(current_path, ApartAppend):
    tspath = tspath.split('/')
    Apart = tspath[1]
    if len(Apart) == 0:
        new_Apart = ApartAppend
    else:
        new_Apart = Apart + '_' + ApartAppend
    tspath[1] = new_Apart
    tspath = '/'.join(tspath)
    return tspath

def getDataLocationDSSInfo(location, currentAlternative, computeOptions):
    if location.isLinkedToPreviousModel():
        tspath = str(currentAlternative.loadTimeSeries(location))
        tspath = fixInputLocationFpart(currentAlternative, tspath)
        dsspath = computeOptions.getDssFilename()
    else:
        tspath = location.getLinkedToLocation().getDssPath()
        rundir = Project.getCurrentProject().getProjectDirectory()
        dsspath = location.getLinkedToLocation().get_dssFile()
        dsspath = os.path.join(rundir, dsspath)
    return tspath, dsspath

def strip_templateID_and_rename_records(dssFilePath,currentAlt):

    # make copy of dss file
    shutil.copyfile(dssFilePath,dssFilePath+'.bak')

    # rename all records, stripping first 4 chars from f-part
    dss = HecDss.open(dssFilePath)
    rec_names = dss.getPathnameList()
    new_rec_names = Vector()
    #currentAlt.addComputeMessage(type(rec_names).__name__)
    for i,r in enumerate(rec_names):		
        #currentAlt.addComputeMessage(type(r).__name__)        
        parts = r.split('/')
        if not '-' in parts[-2]:
            return
        parts[-2] = parts[-2][4:]
        new_rec_names.add('/'.join(parts))
        currentAlt.addComputeMessage('Fixing path: '+r+' --> '+new_rec_names[-1])
        #currentAlt.addComputeMessage(type(new_rec_names).__name__)
        #currentAlt.addComputeMessage(type(new_rec_names[i]).__name__)
        #if i==2:
        #    break
    dss.renameRecords(rec_names, new_rec_names)
    #currentAlt.addComputeMessage(type(new_rec_names).__name__)
    dss.close()

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    output_data = []
    for dsspath in input_data:
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        if len(output_data) == 0:
            output_data = values
        else:
            for vi, val in enumerate(values):
                output_data[vi] += val
                
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_path
    tsc.values = output_data
    tsc.startTime = times[0]
    tsc.units = units
    tsc.endTime = times[-1]
    tsc.numberValues = len(output_data)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod):
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR'''
    dssFm = HecDss.open(inputDSSFile)
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")
    dssFm.close()

    dssFmout = HecDss.open(outputDSSFile)
    dssFmout.write(tsm_new)
    dssFmout.close()


def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    lapse = lapse_in_C
    if 'f' in tsm.getUnits().lower():
        lapse = lapse*9.0/5.0+32.0
    tsm = tsm.add(lapse)
    tsc = tsm.getData()
    dss.close()

    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()


def add_flows(currentAlt, timewindow, inflow_records, dss_file, output_dss_record_name, output_dss_file):
     
    #cfs_2_acreft = balance_period * 3600. / 43559.9
    #acreft_2_cfs = 1. / cfs_2_acreft

    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dss_file)

    inflows = []
    times = []

    # Read inflows
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        try:
       
            print(starttime_str, endtime_str)
            print(dss_file)
            ts = dssFm.read(pathname, starttime_str, endtime_str, False)
            ts_data = ts.getData()
            values = ts_data.values
            hectimes = ts_data.times
            units = ts_data.units
            tstype = ts_data.type
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

    # Output record
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    #tsc.startTime = times[1]
    tsc.units = 'CFS'
    tsc.type = tstype
    #tsc.endTime = times[-1]
    tsc.numberValues = len(inflows)
    #tsc.startHecTime = timewindow.getStartTime()
    #tsc.endHecTime = timewindow.getEndTime()
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    dssFm.close()
    dssFm_out.close()


def add_or_subtract_flows(currentAlt, timewindow, inflow_records, dss_file, operation,
                       output_dss_record_name, output_dss_file):
	# operation: list where True = add, False = subtract, e.g. [True,False,True] to substract the 2nd
	# record from the sun of the first and third records
     
    #cfs_2_acreft = balance_period * 3600. / 43559.9
    #acreft_2_cfs = 1. / cfs_2_acreft

    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dss_file)

    inflows = []
    times = []

    # Read inflows
    print('Reading inflows')
	
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        try:
       
            print(starttime_str, endtime_str)
            print(dss_file)
            ts = dssFm.read(pathname, starttime_str, endtime_str, False)
            ts_data = ts.getData()
            values = ts_data.values
            hectimes = ts_data.times
            units = ts_data.units
            tstype = ts_data.type
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
        	if operation[j]:
	            for vi, v in enumerate(values):
	                inflows[vi] += v
	        else:
	            for vi, v in enumerate(values):
	                inflows[vi] -= v
	                
    # Output record
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    #tsc.startTime = times[1]
    tsc.units = 'CFS'
    tsc.type = tstype
    #tsc.endTime = times[-1]
    tsc.numberValues = len(inflows)
    #tsc.startHecTime = timewindow.getStartTime()
    #tsc.endHecTime = timewindow.getEndTime()
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    dssFm.close()
    dssFm_out.close()


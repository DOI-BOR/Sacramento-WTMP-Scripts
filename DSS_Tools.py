#version 2.0
#modified 03-28-2023 by Scott Burdick-Yahya
#modifed Dec 2023 by Ben Saenz

from hec.heclib.dss import HecDss
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
import hec.hecmath.TimeSeriesMath as tsmath
from com.rma.model import Project
import os,shutil,copy,sys
from java.util import Vector, Date

import datetime
from hec.heclib.util import HecTime

def first_value(dss_file,dss_rec):
    dssFm = HecDss.open(dss_file)        
    tsc = dssFm.get(dss_rec,True)
    dssFm.close()
    return tsc.values[0]

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


def data_from_dss(dss_file,dss_rec,starttime_str, endtime_str):
    dssFm = HecDss.open(dss_file)
    if starttime_str is None and endtime_str is None:
        tsc = dssFm.get(dss_rec,True)
    else:
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
        dtt.append(datetime.datetime.fromtimestamp(timestamp))

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
        dsstype = ts.type
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
    tsc.type = dsstype
    tsc.endTime = times[-1]
    tsc.numberValues = len(output_data)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod):
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR, or downsample.  However, hecmath likes to
    clip of days that don't have the complete 24 hour cycle.  So, we pad here, but there is a chance we ask for data not
    available. The read gives garbage data and doesn't complain.  
    TODO: figure out how to check for bounds for non-midnight start and end times.
    '''
    dssFm = HecDss.open(inputDSSFile)
    if timewindow is not None:
        starttime_str = timewindow.getStartTimeString()
        endtime_str = timewindow.getEndTimeString()
        #if newPeriod.lower() == '1day':  # some computes don't end on 2400, causes problems when last day doesn't get produced in this func
        starttime_str = starttime_str[:-4] + '0000'
        endtime_str = endtime_str[:-4] + '2400' # clipped days don't work in computes ... hope the downloaded DMS data is long enough to do this.
        print('Resampling',newPeriod, inputRec,starttime_str,endtime_str)
        tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    else:
        print('Resampling',newPeriod, inputRec)
        tsm = dssFm.read(inputRec)  # caution - 'read' sometimes doesn't get whole record?  Need to use get?

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

def min_ts(dss_file,dss_rec,min_value,dss_outfile,f_part):
    dss = HecDss.open(dss_file)
    tsc = dss.get(dss_rec,True)
    dss.close()

    for vi, v in enumerate(tsc.values):
        tsc.values[vi] = max(v, min_value)

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
            if '::' in inflow_record:
                dss_file_alt,inflow_rec_alt = inflow_record.split('::')
                dssFm_alt = HecDss.open(dss_file_alt)
                ts = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False)
                dssFm_alt.close()
                print(dss_file_alt)
            else:
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
            if '::' in inflow_record:
                dss_file_alt,inflow_rec_alt = inflow_record.split('::')
                dssFm_alt = HecDss.open(dss_file_alt)
                ts = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False)
                dssFm_alt.close()
                print(dss_file_alt)
            else:            	
                ts = dssFm.read(pathname, starttime_str, endtime_str, False)                
                print(dss_file)
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


def hec_str_time_to_dt(hec_str_time):
    '''Convert HEC date time format to python datetime object'''

    dt_format = '%d%b%Y %H%M'
    add_day = False
    if hec_str_time.endswith('2400'):
        my_hec_str_time = hec_str_time[:-4] + '0000'
        add_day = True
    else:
        my_hec_str_time = hec_str_time

    dt = datetime.datetime.strptime(my_hec_str_time,dt_format)
    if add_day:
        dt = dt + datetime.timedelta(days=1)
    return dt


def create_constant_dss_rec(currentAlt, timewindow, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS', fpart='ZEROS'):
    '''Create and write a dss record with a constant in it for the given time windows.
       what={'flow','temp-water'}
       period={'1HOUR','1DAY'}
    '''

    if what.lower()=='flow':
        units = 'cfs'
        parameter = 'flow'
    elif what.lower()=='temp-water':
        units = 'C'
        parameter = 'temp-water'
    elif what.lower()=='gate':
        units = 'n/a'
        parameter = 'gate'
    elif what.lower()=='evap':
        units = 'ft'
        parameter = 'evap'
    elif what.lower()=='elev':
        units = 'ft'
        parameter = 'elev'
    else:
        currentAlt.addComputeMessage('create_zero_dss_rec: what not known: %s'%what)
        return False

    if period.lower()=='1hour':
        pass
    elif period.lower()=='1day':
        pass
    else:
        currentAlt.addComputeMessage('create_zero_dss_rec: period not known: %s'%period)
        return False

    dt_format = '%d%b%Y %H%M'
    
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #starttime_hectime = HecTime(starttime_str).value()
    #endtime_hectime = HecTime(endtime_str).value()

    # pad 1 day on records, in case these are used for lookbacks, or balance flow calcs, etc.
    starttime_dt = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=1)    
    endtime_dt = hec_str_time_to_dt(endtime_str) + datetime.timedelta(days=1)
    starttime_str_pad = starttime_dt.strftime(dt_format)
    endtime_str_pad = endtime_dt.strftime(dt_format)    
 
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    ########################
    # Zero-Flow Time Series
    ########################

    tsmath_zero_flow_day = tsmath.generateRegularIntervalTimeSeries(
        starttime_str_pad,
        endtime_str_pad,
        period, "0M", constant)
    tsmath_zero_flow_day.setUnits(units)
    tsmath_zero_flow_day.setType(dss_type)
    tsmath_zero_flow_day.setTimeInterval(period)
    tsmath_zero_flow_day.setLocation(cpart)
    tsmath_zero_flow_day.setParameterPart(parameter)
    tsmath_zero_flow_day.setVersion(fpart)

    dssFm = HecDss.open(output_dss_file)
    dssFm.write(tsmath_zero_flow_day)
    dssFm.close()

    return True


def calculate_relative_humidity(air_temp, dewpoint_temp):
    """
    Calculate Relative Humidity given the air temperature and dewpoint temperature - August-Roche-Magnus approximation

    :param air_temp: Air Temperature in degrees Celsius
    :param dewpoint_temp: Dew Point Temperature in degrees Celsius
    :return: Relative Humidity in percentage
    """
    numerator = (112.0 - 0.1 * dewpoint_temp + air_temp)
    denominator = (112.0 + 0.9 * air_temp)
    exponent = ((17.62 * dewpoint_temp) / (243.12 + dewpoint_temp)) - ((17.62 * air_temp) / (243.12 + air_temp))
    relative_humidity = 100.0 * (numerator / denominator) * math.exp(exponent)
    return max(0.01, min(100.0, relative_humidity))


def relhum_from_at_dp(met_dss_file, at_path, dp_path):
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    dp_data = DSS_Tools.data_from_dss(met_dss_file, dp_path, None, None)
    for i in range(tsc.numberValues):
        tsc.values[i] = calculate_relative_humidity(tsc.values[i], dp_data[i])
    parts = tsc.fullName.split('/')
    parts[2] = parts[2][:5]
    parts[3] = 'RELHUM-FROM-AT-DP'
    parts[6] = parts[6] + '-DERIVED'
    new_pathname = '/'.join(parts)
    tsc.fullName = new_pathname
    tsc.units = '%'
    print('writing: ', new_pathname)
    dss.write(tsc)
    dss.close()

def check_start_and_end(values, times, startime, endtime):
    if times[0] < startime:  # if startdate is before the timewindow..
        print('start date ({0}) from DSS before timewindow ({1})..'.format(times[0], startime))
        st_offset = (startime - times[0]) / (times[1] - times[0])
        values = values[st_offset:]
        times = times[st_offset:]
    if times[-1] > endtime:
        print('end date ({0}) from DSS after timewindow ({1})..'.format(times[-1], endtime))
        st_offset = (times[-1] - endtime) / (times[1] - times[0])
        values = values[:(len(times) - st_offset)]
        times = times[:(len(times) - st_offset)]
    return values, times


def replace_data(currentAlt, timewindow, pairs, dss_file, dss_outfile, months, standard_interval=None):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    # 01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    for pair in pairs:
        dssFm = HecDss.open(dss_file)
        currentAlt.addComputeMessage('Replacing data for {0} with {1} during {2}'.format(pair[0], pair[1], months))
        base = dssFm.read(pair[0], starttime_str, endtime_str, False)
        if standard_interval is not None:
            base = standardize_interval(base,standard_interval)
        base_data = base.getData()
        base_values = base_data.values
        base_hectimes = base_data.times
        base_units = base_data.units
        base_interval = base_data.interval
        base_type = base_data.type
        base_values, base_hectimes = check_start_and_end(base_values, base_hectimes, starttime_hectime, endtime_hectime)

        alt = dssFm.read(pair[1], starttime_str, endtime_str, False)
        if standard_interval is not None:
            alt = standardize_interval(alt,standard_interval)
        alt_data = alt.getData()
        alt_values = alt_data.values
        alt_hectimes = alt_data.times
        alt_units = alt_data.units
        alt_interval = alt_data.interval
        alt_values, alt_hectimes = check_start_and_end(alt_values, alt_hectimes, starttime_hectime, endtime_hectime)

        dssFm.close()

        if base_units != alt_units:
            currentAlt.addComputeMessage('Units do not match for {0} and {1}, skipping'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)
        if base_interval != alt_interval:
            currentAlt.addComputeMessage('Intervals do not match for {0} and {1}, changing interval...'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)

        for i in range(len(base_values)):
            if base_data.getHecTime(i).month() in months:
                base_values[i] = alt_values[i]
                #print('replaced {0}'.format(base_data.getHecTime(i)))

        new_pathname = base_data.fullName.split('/')
        alt_pathname = alt_data.fullName.split('/')
        new_pathname[-2] = 'MergedFrom_{0}'.format(alt_pathname[1])
        new_pathname = '/'.join(new_pathname)
        print('writing: ',new_pathname)
        tsc = TimeSeriesContainer()
        tsc.times = base_hectimes
        tsc.fullName = new_pathname
        tsc.values = base_values
        tsc.units = base_units
        tsc.type = base_type
        tsc.numberValues = len(base_values)

        dssFmOut = HecDss.open(dss_outfile)
        dssFmOut.write(tsc)

        dssFm.close()
        dssFmOut.close()

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
import os,shutil,copy,sys,math
from java.util import Vector, Date

import datetime
from hec.heclib.util import HecTime


# --- Utility Function for Python 2 timedelta.total_seconds() ---
# In Python 2, timedelta objects don't always have total_seconds(), so
# we implement it manually to ensure cross-version compatibility for float division.
def _timedelta_to_seconds(td):
    """
    Converts a datetime.timedelta object to a floating-point number of seconds.
    This ensures proper float division for time ratios in Python 2.

    Args:
        td (datetime.timedelta): The time delta to convert.

    Returns:
        float: The total duration represented by td, in seconds.
    """
    # 86400 seconds in a day (24 * 60 * 60)
    # 1,000,000 microseconds in a second
    return td.days * 86400.0 + td.seconds + td.microseconds / 1000000.0

# --- Main Interpolation Method ---
def linear_interp_datetime(time_data, value_data, query_times):
    """
    Performs linear interpolation on time-series data.

    Arguments:
    - time_data: List of ordered (ascending) datetime.datetime objects (X-axis).
    - value_data: List of float values corresponding to time_data (Y-axis).
    - query_times: List of datetime.datetime objects for which to interpolate values.

    Returns:
    - A list of interpolated float values corresponding to query_times.
    """
    # 1. Basic Validation
    if len(time_data) != len(value_data) or len(time_data) < 2:
        print "Error: time_data and value_data must have the same length (and at least 2 points)."
        return []

    # 2. Normalize Times to Numerical Values (Seconds from the first timestamp)
    # Using the first time point (t0) as the numerical zero reference (x_data[0] = 0.0)
    t0 = time_data[0]
    x_data = []
    # convert each known timestamp into seconds elapsed since t0
    for t in time_data:
        x_data.append(_timedelta_to_seconds(t - t0))

    x_query = []
    # convert each query timestamp into seconds elapsed since t0, using the same reference
    for t in query_times:
        x_query.append(_timedelta_to_seconds(t - t0))

    # 3. Perform Interpolation / Extrapolation
    results = []
    n = len(x_data)

    # for each query time, either extrapolate (if outside the known range) or interpolate
    for x in x_query:
        # A. Extrapolation: Before the first data point
        if x <= x_data[0]:
            # Use the closest float value (the first value)
            results.append(value_data[0])
            continue

        # B. Extrapolation: After the last data point
        if x >= x_data[n - 1]:
            # Use the closest float value (the last value)
            results.append(value_data[n - 1])
            continue

        # C. Interpolation: Find the interval [x_data[i], x_data[i+1]]
        # We search for the index 'i' such that x_data[i] <= x < x_data[i+1]
        i = 0
        while i < n - 1 and x >= x_data[i + 1]:
            i += 1

        # Define the segment (x0, y0) and (x1, y1)
        x0, y0 = x_data[i], value_data[i]
        x1, y1 = x_data[i + 1], value_data[i + 1]

        # Check for zero duration (in case of duplicate timestamps in time_data)
        dx = x1 - x0
        if dx == 0.0:
            y = y0
        else:
            # Linear interpolation formula: y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
            y = y0 + (y1 - y0) * (x - x0) / dx

        results.append(y)

    return results


def copy_dss_ts(dss_rec,new_fpart=None,new_dss_rec=None,
                dss_file_path=None,dss_file_handle=None,dss_file_alt_outpath=None, checkMakeCelsius=False):
    """
    Copy a DSS time series record to a new record path, either by
    substituting a new F-part or by using an entirely new record
    path, optionally normalizing degree units and converting
    Fahrenheit to Celsius, and optionally writing to a different
    DSS file.

    Args:
        dss_rec (str): DSS path of the source record to copy.
        new_fpart (str, optional): If provided (and new_dss_rec is
            not), the F-part to substitute into dss_rec's path to
            form the output path.
        new_dss_rec (str, optional): If provided, used directly as
            the output DSS path instead of substituting an F-part.
        dss_file_path (str, optional): Path to the DSS file to open
            for reading (and writing, unless
            dss_file_alt_outpath is given). Required if
            dss_file_handle is not provided.
        dss_file_handle (optional): An already-open DSS file handle
            to read from, instead of opening dss_file_path.
        dss_file_alt_outpath (str, optional): If provided, the
            copied record is written to this DSS file instead of
            the source file.
        checkMakeCelsius (bool, optional): If True and the source
            record's units are Fahrenheit, converts all values to
            Celsius before writing. Defaults to False.

    Returns:
        None. Writes the copied (and possibly converted) record to
        the appropriate output DSS file.

    Raises:
        ValueError: If neither dss_file_path nor dss_file_handle is
            provided, or if neither new_fpart nor new_dss_rec is
            provided.
    """
    # error check inputs - there are flexible way to copy record
    if dss_file_path is None and dss_file_handle is None:
        raise ValueError('copy_dss_rec_to_new_fpart: you must supply either a valid dss_file_path OR dss_file_handle')
    if new_fpart is None and new_dss_rec is None:
        raise ValueError('copy_dss_rec_to_new_fpart: you must supply either a new_fpart OR new_dss_rec')

    # (open dss) get record tsc
    if dss_file_handle is not None:
        dss_fm = dss_file_handle
        print('copy_dss_ts - reading input from open file:',dss_rec)
    else:
        dss_fm = HecDss.open(dss_file_path)
        print('copy_dss_ts - reading input:',dss_file_path,dss_rec)
    tsc = dss_fm.get(dss_rec,True)
    print('tsc.fullName',tsc.fullName)

    # new dss rec name
    # build the output path, either using the explicit new record name or by substituting a new F-part into the source path
    if new_dss_rec is not None:
        dss_rec_out = new_dss_rec
    else:
        rec_parts = tsc.fullName.split('/')
        print(rec_parts)
        rec_parts[6] = new_fpart
        dss_rec_out = '/'.join(rec_parts)    

    # fix some terrible units along the way (normalize non-standard degree unit labels to the standard C/F abbreviations)
    if tsc.units.lower() == 'degc':
        tsc.units = 'C'
    elif tsc.units.lower() == 'degf':
        tsc.units = 'F'

    if tsc.units.lower() == 'f' and checkMakeCelsius:
        T_values = tsc.values
        for i, TT in enumerate(T_values):
            T_values[i] = (TT-32.0)*5.0/9.0
        tsc.units = 'C'
        tsc.values = T_values
        
    # write
    tsc.fullName = dss_rec_out
    # write to the alternate output file if one was given, otherwise write back to the source file
    if dss_file_alt_outpath is not None:
        dss_fm_alt = HecDss.open(dss_file_alt_outpath)
        dss_fm_alt.put(tsc)
        dss_fm_alt.close()
    else:
        dss_fm.put(tsc)

    # only close the file handle if this function opened it itself
    if dss_file_handle is None:
        dss_fm.close()

def jday_from_tsc(tsc):
    """
    Compute the fractional (decimal) day-of-year for every timestamp
    in a time series container.

    Args:
        tsc: A DSS time series data container (TimeSeriesContainer)
            with time and value data.

    Returns:
        list of float: The decimal day-of-year for each timestamp
            in tsc, in the same order.
    """
    dtt = hectime_to_datetime(tsc)
    return [decimal_doy(dt) for dt in dtt]        


def decimal_doy(dt):
    """
    Compute the fractional (decimal) day-of-year for a single
    datetime, where the integer part is the calendar day-of-year
    and the fractional part represents the time of day.

    Args:
        dt (datetime.datetime): The date/time to convert.

    Returns:
        float: The decimal day-of-year (e.g. Jan 1 at noon would be
            approximately 1.5).
    """
    doy = dt.timetuple().tm_yday
    fractional_day = (dt.hour / 24.0) + (dt.minute / 1440.0) + (dt.second / 86400.0) + (dt.microsecond / 86400000000.0)
    return doy + fractional_day


def organizeLocations(curAlt, location_objs, loc_names, return_dss_paths=False):
    """
    Look up a list of named locations within a set of location
    objects, in the order given by loc_names, optionally resolving
    each to its DSS path instead of returning the location object
    itself.

    Args:
        curAlt: The alternative object being computed. Passed
            through to findLocationOrder and used for resolving
            linked time series paths.
        location_objs (list): The full list of location objects to
            search within (e.g. from
            currentAlternative.getInputDataLocations()).
        loc_names (list of str): Names to look up, in the desired
            output order.
        return_dss_paths (bool, optional): If True, resolves and
            returns each location's DSS path string instead of the
            location object itself. Defaults to False.

    Returns:
        list: If return_dss_paths is False, a list of location
            objects corresponding to loc_names, in that order. If
            True, a list of resolved DSS path strings instead.
            Exits the process (sys.exit(1)) via findLocationOrder if
            any name in loc_names cannot be found.
    """
    locations_list = []
    print('num_locs:',len(location_objs))
    # for each requested location name, find its matching location object and resolve it
    for name in loc_names:
        i_loc = findLocationOrder(curAlt,location_objs,name)
        print('name:',name,'i_loc:',i_loc)
        if return_dss_paths:
            lo1 = location_objs[i_loc]
            print(lo1)
            # if this location is linked to an upstream model, load its time series and fix
            # its F-part; otherwise, it already has a fixed DSS path to use directly
            if lo1.isLinkedToPreviousModel():
                tspath = str(curAlt.loadTimeSeries(lo1))
                tspath = fixInputLocationFpart(curAlt, tspath)
            else:
                tspath = str(lo1.getDssPath())
            print(tspath)
            locations_list.append(tspath)
        else:
            locations_list.append(location_objs[i_loc])
    return locations_list


def organizeLocationsPaired(curAlt, location_objs, loc_names_paired, return_dss_paths=False):   
    """
    Apply organizeLocations to multiple groups of location names at
    once, returning one resolved list per group.

    Args:
        curAlt: The alternative object being computed, passed
            through to organizeLocations.
        location_objs (list): The full list of location objects to
            search within.
        loc_names_paired (list of list of str): A list of name
            groups; each inner list is passed to organizeLocations
            independently.
        return_dss_paths (bool, optional): Passed through to
            organizeLocations for each group. Defaults to False.

    Returns:
        list of list: One resolved list (of location objects or DSS
            paths, depending on return_dss_paths) per input group in
            loc_names_paired, in the same order.
    """
    return [organizeLocations(curAlt, location_objs, pn, return_dss_paths) for pn in loc_names_paired]


def findLocationOrder(curAlt,location_objs,name):
    """
    Find the index of a location object within a list, matching by
    name.

    Args:
        curAlt: The alternative object being computed. Used to log
            an error message if the name cannot be found.
        location_objs (list): The list of location objects to
            search within. Each must support getName().
        name (str): The location name to search for.

    Returns:
        int: The index of the first location object whose
            getName() matches name. Exits the process
            (sys.exit(1)) if no match is found.
    """
    # search through every location for one whose name matches
    for i,loc in enumerate(location_objs):
        print(i,'Checking loc: ',loc.getName())
        print('loc:',loc)
        if name == loc.getName():
            return i
    # if we make it here, our input/output location name was not found
    curAlt.addComputeMessage("Scripting - Location name not found: "+name)
    sys.exit(1)

def first_value(dss_file,dss_rec,start_str=None,end_str=None):
    """
    Read a DSS record and return only its first value.

    Args:
        dss_file (str): Path to the DSS file to read from.
        dss_rec (str): DSS path of the record to read.
        start_str (str, optional): Start time string for a bounded
            read. If None (along with end_str), the entire record
            is read instead.
        end_str (str, optional): End time string for a bounded
            read. If None (along with start_str), the entire record
            is read instead.

    Returns:
        float: The first value in the record.
    """
    dssFm = HecDss.open(dss_file)        
    # if no time bounds were given, read the entire record; otherwise read only the given window
    if start_str is None and end_str is None:
        tsc = dssFm.get(dss_rec,True)
    else:
        tsc = dssFm.read(dss_rec,start_str,end_str,False).getData()
    dssFm.close()
    return tsc.values[0]


def standardize_interval(tsm, interval, makePerAver=True):
    """
    Resample a time series math object to a standard interval
    (hourly, daily, or weekly), if it is not already at that
    interval, averaging values within each new interval.

    Args:
        tsm: A HEC time series math object (tsmath) to standardize.
        interval (str): The target interval. Must be one of
            '1hour', '1day', or '1week' (case-insensitive).
        makePerAver (bool, optional): If True, sets the series type
            to 'PER-AVER' before transforming, ensuring the
            resample treats values as period averages. Defaults to
            True.

    Returns:
        A tsmath object at the requested interval. If tsm is
            already at that interval, tsm is returned unchanged
            (no averaging is performed).

    Note:
        Exits the process (sys.exit(-1)) if interval is not one of
        the recognized values.
    """
    tsc = tsm.getData()
    # map the requested interval name to its corresponding HEC interval-minutes value
    if interval.lower()=='1hour':
        intint = 60
    elif interval.lower()=='1day':
        intint=1440
    elif interval.lower()=='1week':
        intint=10800
    else:
        print('interval not supported:',interval)
        sys.exit(-1)

    # only resample if the current interval differs from the requested one
    if tsc.interval != intint:
        if makePerAver:
            #tsc.type = 'PER-AVER'  # make sure it's per-aver ... we are
            tsm.setType('PER-AVER')
        return tsm.transformTimeSeries(interval, "", "AVE")
    else:
        return tsm


def get_sanitized_record_list(dss_file_path):
    """
    The DSS library seems to return lists of paths with dates in them (getPathnameList()), and some of those
    dates don't even exist in the file or cannot be read and throw an error.  As of Jan 2024,
    this is an orregular problem, and the manual soluation is throwing away the DSS file but
    in many cases that is problematic. So, here we filter dates and check for duplicates.

    Args:
        dss_file_path (str): Path to the DSS file to list records
            from.

    Returns:
        list of str: A deduplicated list of DSS record paths with
            the date (D-part) field blanked out, so that records
            differing only by date are treated as the same record.
    """
    dss = HecDss.open(dss_file_path)
    recs = dss.getPathnameList()
    dss.close()
    sanitized_recs = []
    
    # for each raw record path, blank out its date field and skip if already seen
    for r in recs:
        rec_tokens = r.split('/')
        rec_tokens[4] = ''  # erase date string, if it exists
        r_sanitized = '/'.join(rec_tokens)
        if not r_sanitized in sanitized_recs:
            sanitized_recs.append(r_sanitized)
    return sanitized_recs  

def hectimes_from_tsm(tsm):
    """
    Extract a list of HecTime objects from a time series math
    object's underlying data container.

    Args:
        tsm: A HEC time series math object (tsmath).

    Returns:
        list of HecTime: One HecTime object per value in tsm,
            constructed from its raw integer time values.
    """
    times = tsm.getContainer().times
    htimes = []
    for i in range(tsm.getContainer().numberValues):
        htimes.append(HecTime())
        htimes[-1].set(times[i])
    return htimes

def hectimes_from_tsc(tsc):
    """
    Extract a list of HecTime objects from a time series data
    container.

    Args:
        tsc: A DSS time series data container (TimeSeriesContainer).

    Returns:
        list of HecTime: One HecTime object per value in tsc,
            constructed from its raw integer time values.
    """
    htimes = []
    # for each raw time value, build a corresponding HecTime object
    for i in range(tsc.numberValues):
        htimes.append(HecTime())
        htimes[-1].set(tsc.times[i])
    return htimes


def shift_pit_river_time(input_dss_file,dss_rec,output_dss_file,out_rec,start_date=None,end_date=None):
    """
    Read a time series, shift it 12 hours earlier in time, and
    write the result to a new record.

    Note: it is important that dss_rec and out_rec differ. If they
    are the same, running this function repeatedly will keep
    shifting the same record further back in time on every run.

    Args:
        input_dss_file (str): Path to the DSS file containing the
            source record.
        dss_rec (str): DSS path of the source record to shift.
        output_dss_file (str): Path to the DSS file to write the
            shifted record to.
        out_rec (str): DSS path for the shifted output record. Must
            differ from dss_rec.
        start_date (str, optional): Start time string for a bounded
            read. If None, the entire record is read.
        end_date (str, optional): End time string for a bounded
            read. If None, the entire record is read.

    Returns:
        None. Writes the shifted record to output_dss_file.
    """
    # important to not have dss_rec==out_rec, otherwise you will continually shift the record in time
    # whenever the calling script runs...    

    tsm_pit = dss_read_ts_safe(input_dss_file,dss_rec,start_date=start_date,end_date=end_date,returnTSM=True)

    tsm_pit = tsm_pit.shiftInTime("-12Hour")
    tsc_pit = tsm_pit.getData()
    tsc_pit.fullName = out_rec  

    dss_out = HecDss.open(output_dss_file)
    dss_out.put(tsc_pit)
    dss_out.close()


def shift_ts_time(input_dss_file,dss_rec,output_dss_file,out_rec,shift_str,start_date=None,end_date=None):
    """
    Read a time series, shift it in time by an arbitrary HEC-style
    offset, and write the result to a new record.

    Note: it is important that dss_rec and out_rec differ. If they
    are the same, running this function repeatedly will keep
    shifting the same record further in time on every run.

    Args:
        input_dss_file (str): Path to the DSS file containing the
            source record.
        dss_rec (str): DSS path of the source record to shift.
        output_dss_file (str): Path to the DSS file to write the
            shifted record to.
        out_rec (str): DSS path for the shifted output record. Must
            differ from dss_rec.
        shift_str (str): A HEC-recognized time shift string, with an
            optional leading minus sign (e.g. '-12Hour', '5Day').
        start_date (str, optional): Start time string for a bounded
            read. If None, the entire record is read.
        end_date (str, optional): End time string for a bounded
            read. If None, the entire record is read.

    Returns:
        None. Writes the shifted record to output_dss_file.
    """
    
    tsm = dss_read_ts_safe(input_dss_file,dss_rec,start_date=start_date,end_date=end_date,returnTSM=True)

    tsm = tsm.shiftInTime(shift_str)
    tsc = tsm.getData()
    tsc.fullName = out_rec  

    dss_out = HecDss.open(output_dss_file)
    dss_out.put(tsc)
    dss_out.close()


def dss_read_ts_safe(dssFilePath,dssRec,start_date=None,end_date=None,returnTSM=False,
                     returnPydatetimes=False,debug=False):
    """
    A read function that is date-flexible, and ensures the whole time range is returned if dates
    are None.

    Args:
        dssFilePath (str): Path to the DSS file to read from.
        dssRec (str): DSS path of the record to read.
        start_date (str, optional): Start time string for a bounded
            read. If both start_date and end_date are None, the
            entire record is read via get() instead.
        end_date (str, optional): End time string for a bounded
            read. If both start_date and end_date are None, the
            entire record is read via get() instead.
        returnTSM (bool, optional): If True, returns the tsmath
            wrapper object instead of the raw data container.
            Defaults to False.
        returnPydatetimes (bool, optional): Currently unused within
            this function body.
        debug (bool, optional): If True, prints extra diagnostic
            information about the read. Defaults to False.

    Returns:
        The time series data, either as a raw data container or as
        a tsmath object (depending on returnTSM), or None if neither
        the "both dates None" nor "both dates given" branch is
        satisfied (e.g. if only one of start_date/end_date is
        provided).
    """
    # this method is going to throw an error if the file doesn't exist
    dss = HecDss.open(dssFilePath,True)

    # recordExists() seems to check for the individual database "pages" or something, with a date, or
    # date range, required?  TODO: figure out how to check if a date-agnostic record exists.  Or, don't,
    # since HECLIB produces a pretty nice error if you try to open a non-existant record.
    #if not dss.recordExists(dssRec):
    #    print('DSS rec does not exist for reading:')
    #    print('File: '+dssFilePath)
    #    print('Rec: '+dssRec)
    #    sys.exit(-1)
    #else:
    if start_date is None and end_date is None:
        tsc = dss.get(dssRec,True) # use get with True here to capture entire record, 'read' seems to leave off data randomly
        dss.close()
        if debug:
            print('Reading DSS in script...')
            print('    file: '+dssFilePath)
            print('    record: '+dssRec)
        if returnTSM:
            return tsmath(tsc)
        else:
            return tsc
    # if both bounds were given, read only within that time window
    elif start_date is not None and end_date is not None:
        tsm = dss.read(dssRec,start_date,end_date,False) # 'read' allows time windows in call, but returns tsm
        dss.close()
        if debug:
            print('Reading DSS in script between '+start_date+' and ',+end_date)
            print('    file: '+dssFilePath)
            print('    record: '+dssRec)
        if returnTSM:
            return tsm
        else:
            return tsm.getData()


def data_from_dss(dss_file,dss_rec,starttime_str, endtime_str):
    """
    Read a DSS record's values, either over its entire span or
    within a specific time window.

    Args:
        dss_file (str): Path to the DSS file to read from.
        dss_rec (str): DSS path of the record to read.
        starttime_str (str): Start time string for a bounded read.
            If None (along with endtime_str), the entire record is
            read instead.
        endtime_str (str): End time string for a bounded read. If
            None (along with starttime_str), the entire record is
            read instead.

    Returns:
        list of float: The values from the requested record/time
            window.
    """
    dssFm = HecDss.open(dss_file)
    # if no time bounds were given, read the entire record; otherwise read only the given window
    if starttime_str is None and endtime_str is None:
        tsc = dssFm.get(dss_rec,True)
    else:
        tsc = dssFm.read(dss_rec, starttime_str, endtime_str, False).getData()
    dssFm.close()
    return tsc.values


def hectime_to_datetime(tsc):
    """
    Deprecated! This routine, and the javatime processing it depends on, is capable of rounding dates 
    to the wrong hour.  It also adds a timezone offset and returns GMT, which is not usually wanted.
    
    TODO: replace all instances of this with datetimes_from_tsc(), and rework/remove any extra processing
    previously needed for timezones, bad first/last values, etc.

    Args:
        tsc: A DSS time series data container (TimeSeriesContainer).

    Returns:
        list of datetime.datetime: One Python datetime per value in
            tsc, converted via Java Date/timestamp (in GMT, with
            potential rounding issues - see deprecation note above).
    """

    hectimes = hectimes_from_tsc(tsc)
    dtt = []
    # for each timestamp, convert via Java Date/timestamp into a Python datetime
    for j in range(tsc.numberValues):
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        java_date = tsc.getHecTime(j).getJavaDate(0)  
        
        # Convert Java Date to Python datetime
        timestamp = (java_date.getTime() / 1000)
        d = datetime.datetime.fromtimestamp(timestamp)
        #print('hectime_to_datetime:',tsc.getHecTime(j).toString(),timestamp,d.strftime('%d%b%Y %H%M'))
        dtt.append(datetime.datetime.fromtimestamp(timestamp))

    return dtt

def hecTime_to_datetime(hectimes):
    """
    Convert a list of HecTime objects to Python datetime objects,
    via their long-format string representation.

    Args:
        hectimes (list of HecTime): The HecTime objects to convert.

    Returns:
        list of datetime.datetime: One Python datetime per input
            HecTime object.
    """
    return [hec_str_time_to_dt(h.toString(),longStr=True) for h in hectimes]

def datetimes_from_tsc(tsc):
    """
    Convert all timestamps in a time series data container directly
    to Python datetime objects, without going through the
    Java-Date-based (and deprecated) hectime_to_datetime path.

    Args:
        tsc: A DSS time series data container (TimeSeriesContainer).

    Returns:
        list of datetime.datetime: One Python datetime per value in
            tsc.
    """
    return hecTime_to_datetime(hectimes_from_tsc(tsc))

def hec_str_time_to_dt(hec_str_time,longStr=False):
    """
    Convert HEC date time format to python datetime object

    Args:
        hec_str_time (str): The HEC-formatted date/time string to
            convert.
        longStr (bool, optional): If True, expects the long HEC
            string format ("%d %B %Y, %H:%M", with '24:00'/'00:00'
            markers). If False, expects the short format
            ('%d%b%Y %H%M', with '2400'/'0000' markers). Defaults
            to False.

    Returns:
        datetime.datetime: The parsed date/time. HEC's "24:00"
            end-of-day convention is converted to 00:00 of the
            following day.
    """
    add_day = False
    # select the expected string format and end-of-day markers based on longStr
    if longStr:
        dt_format = "%d %B %Y, %H:%M"
        check24 = '24:00'
        check00 = '00:00'
    else:
        dt_format = '%d%b%Y %H%M'
        check24 = '2400'
        check00 = '0000'

    checkLen = len(check24)*(-1)
    # if the time string uses HEC's "24:00" end-of-day convention, convert it to the
    # following day's 00:00 instead, since Python's datetime cannot represent hour 24
    if hec_str_time.endswith(check24):
        my_hec_str_time = hec_str_time[:checkLen] + check00
        add_day = True
    else:
        my_hec_str_time = hec_str_time

    dt = datetime.datetime.strptime(my_hec_str_time,dt_format)
    if add_day:
        dt = dt + datetime.timedelta(days=1)
    return dt

def fixInputLocationFpart(currentAlternative, tspath):
    """
    Rewrite a DSS path's F-part so that its version prefix matches
    the current alternative's configured input F-part, while
    preserving the original F-part's trailing scenario/model
    identifier.

    Args:
        currentAlternative: The alternative object being computed.
            Must support getInputFPart() to provide the expected
            F-part prefix.
        tspath (str): The DSS path whose F-part should be corrected.

    Returns:
        str: The DSS path with its F-part replaced by the
            alternative's configured F-part prefix, combined with
            the original path's trailing F-part segment.
    """
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])
    tspath = tspath.split('/')
    fpart = tspath[6]
    fpart_split = fpart.split(':')
    new_fpart = new_fpart_start + ':' + fpart_split[-1]
    tspath[6] = new_fpart
    tspath = '/'.join(tspath)
    return tspath

def appendAPart(current_path, ApartAppend):
    """
    Append a suffix onto a DSS path's A-part (river/basin name),
    separated by an underscore, or set the A-part directly if it
    was previously empty.

    NOTE: This function references the name `tspath` before it is
    ever assigned - the parameter is named `current_path`, but the
    first line splits `tspath`, which does not exist yet. As
    written, calling this function will raise a NameError. This has
    been left unchanged and is flagged here for visibility.

    Args:
        current_path (str): The DSS path whose A-part should be
            modified. (See NOTE above regarding a bug in the
            current implementation.)
        ApartAppend (str): The text to append to (or set as) the
            A-part.

    Returns:
        str: The DSS path with its A-part updated.
    """
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
    """
    Sum multiple DSS time series records together, value by value,
    and write the combined series to a new output record.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        dssFile (str): Path to the DSS file to read from and write
            to.
        timewindow: The run time window object, used to get the
            start and end time strings for reading input data.
        input_data (list of str): DSS paths of the records to sum
            together. All records are assumed to share the same
            timestamps, units, and type.
        output_path (str): DSS path to write the combined series to.

    Returns:
        int: Always returns 0 on completion.
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    output_data = []
    # for each input record, read it and add its values into the running total
    for dsspath in input_data:
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        dsstype = ts.type
        # the first record seeds the running total; subsequent records are added on
        if len(output_data) == 0:
            output_data = values
        else:
            # add this record's values onto the running total, timestep by timestep
            for vi, val in enumerate(values):
                output_data[vi] += val
                
    # build and write the combined time series container
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

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod, lookback_1mon=False, pad_start_days=0):
    '''Can upsample an even period DSS timeseries, e.g. go from 1DAY -> 1HOUR, or downsample.  However, hecmath likes to
    clip of days that don't have the complete 24 hour cycle.  So, we pad here, but there is a chance we ask for data not
    available. The read gives garbage data and doesn't complain.  
    TODO: figure out how to check for bounds for non-midnight start and end times.
    '''
    dssFm = HecDss.open(inputDSSFile)
    # if a time window was given, pad it out to full-day boundaries (and optionally
    # extra lookback days) before reading, since resampling clips incomplete days
    if timewindow is not None:
        starttime_str = timewindow.getStartTimeString()
        endtime_str = timewindow.getEndTimeString()
        #if newPeriod.lower() == '1day':  # some computes don't end on 2400, causes problems when last day doesn't get produced in this func
        starttime_str = starttime_str[:-4] + '0000'
        endtime_str = endtime_str[:-4] + '2400' # clipped days don't work in computes ... hope the downloaded DMS data is long enough to do this.
        # optionally pad extra days onto the start of the window
        if pad_start_days > 0:
            dt_start = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=pad_start_days)
            starttime_str = dt_start.strftime('%d%b%Y %H%M')
        if lookback_1mon:
            dt_start = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=-31)
            starttime_str = dt_start.strftime('%d%b%Y %H%M')      
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
    """
    Apply a fixed temperature lapse (offset) to an air temperature
    record, to approximate temperature change with elevation, and
    write the adjusted record under a new F-part.

    NOTE: This function is defined twice in this module (once here,
    and again later, identically). The later definition will
    silently override this one at import time. Consider removing
    the duplicate.

    Args:
        dss_file (str): Path to the DSS file containing the source
            air temperature record.
        dss_rec (str): DSS path of the source air temperature
            record.
        lapse_in_C (float): The temperature lapse/offset to apply,
            in degrees Celsius. If the source record's units are
            Fahrenheit, this is converted to an equivalent
            Fahrenheit offset before being applied.
        dss_outfile (str): Path to the DSS file to write the
            adjusted record to.
        f_part (str): The F-part to assign to the output record.

    Returns:
        None. Writes the lapse-adjusted record to dss_outfile.
    """
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    lapse = lapse_in_C
    # if the source series is in Fahrenheit, convert the lapse amount to an equivalent Fahrenheit offset
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
    """
    Clamp every value in a time series to a minimum value, and
    write the adjusted record under a new F-part.

    Args:
        dss_file (str): Path to the DSS file containing the source
            record.
        dss_rec (str): DSS path of the source record to clamp.
        min_value (float): The minimum allowed value; any value
            below this is raised to this value.
        dss_outfile (str): Path to the DSS file to write the
            adjusted record to.
        f_part (str): The F-part to assign to the output record.

    Returns:
        None. Writes the clamped record to dss_outfile.
    """
    dss = HecDss.open(dss_file)
    tsc = dss.get(dss_rec,True)
    dss.close()

    # for each value, raise it to the minimum if it falls below that threshold
    for vi, v in enumerate(tsc.values):
        tsc.values[vi] = max(v, min_value)

    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()

def add_flows(currentAlt, timewindow, inflow_records, dss_file, output_dss_record_name, output_dss_file):
     """
    Sum multiple flow time series together (with unit conversion
    from cms to cfs where needed, and trimming to the run time
    window), and write the combined flow record to a new output
    file.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        timewindow: The run time window object, used to determine
            the valid time range for trimming input records.
        inflow_records (list of str): DSS paths of the records to
            sum. A path may optionally use the form
            'other_dss_file::record_path' to read from a different
            DSS file than dss_file.
        dss_file (str): Path to the primary DSS file to read
            records from (used unless a record path specifies an
            alternate file via '::').
        output_dss_record_name (str): DSS path to write the summed
            result to.
        output_dss_file (str): Path to the DSS file to write the
            summed result to.

    Returns:
        None. Writes the combined flow series (in CFS) to
        output_dss_file.

    Raises:
        SystemExit: Calls sys.exit(-1) if a HecMathException occurs
            while reading any of the input records.
    """
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
    # for each inflow record, read its data (from an alternate file if specified),
    # trim it to the time window, convert units if needed, and add it to the running total
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        try:
       
            print(starttime_str, endtime_str)
            # if the record path specifies an alternate DSS file via '::', read from that file instead
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
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # the first record seeds the running total; subsequent records are added on
        if len(inflows) == 0:
            inflows = values
            times = hectimes #TODO: check how this handles missing values
        else:
            # add this record's values onto the running total, timestep by timestep
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
                       output_dss_record_name, output_dss_file, what="flow", prepend_n=0):
    """
    Combine multiple time series by adding or subtracting each one
    according to a matching list of True/False operations, with
    optional unit conversion (for flow) and optional prepending of
    lookback values.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        timewindow: The run time window object, used to determine
            the valid time range for trimming input records.
        inflow_records (list of str): DSS paths of the records to
            combine. A path may optionally use the form
            'other_dss_file::record_path' to read from a different
            DSS file than dss_file.
        dss_file (str): Path to the primary DSS file to read
            records from (used unless a record path specifies an
            alternate file via '::').
        operation (list of bool): One entry per record in
            inflow_records. True means add that record to the
            running total; False means subtract it. The first
            record's value is always used as-is to seed the total.
        output_dss_record_name (str): DSS path to write the combined
            result to.
        output_dss_file (str): Path to the DSS file to write the
            combined result to.
        what (str, optional): If "flow" (default), values are
            treated as flow and converted from cms to cfs as
            needed, and the output units are set to 'CFS'. If any
            other string, no unit conversion is performed and the
            output units are set to that string directly.
        prepend_n (int, optional): If greater than 0, prepends this
            many copies of the first value (and corresponding
            earlier timestamps) onto the front of the combined
            series, e.g. to provide lookback values for downstream
            models. Defaults to 0.

    Returns:
        None. Writes the combined series to output_dss_file.

    Raises:
        SystemExit: Calls sys.exit(-1) if a HecMathException occurs
            while reading any of the input records.
    """
    # operation: list where True = add, False = subtract, e.g. [True,False,True] to substract the 2nd
    # record from the sum of the first and third records
    # what == "flow": assume flows and rectify units
    # what == anything else: use what as units, and don't convert anything

    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('add_or_subtract_flows - Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dss_file)

    inflows = []
    times = []

    # Read inflows
    print('Reading inflows')
    
    # for each record, read its data (from an alternate file if specified), trim it to
    # the time window, convert units if needed, and add or subtract it from the running total
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        try:
       
            print(starttime_str, endtime_str)
            # if the record path specifies an alternate DSS file via '::', read from that file instead
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
            # if the record's data extends past the run window, trim the trailing values
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
                

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        # if treating values as flow and this record is in cms, convert it to cfs
        if what=="flow":
            if units.lower() == 'cms':
                currentAlt.addComputeMessage('Converting cms to cfs')
                convvals = []
                for flow in values:
                    convvals.append(flow * 35.314666213)
                values = convvals

        # the first record always seeds the running total directly; later records are
        # added or subtracted according to the matching operation flag
        if len(inflows) == 0:
            inflows = values
            times = hectimes #TODO: check how this handles missing values
        else:
            if operation[j]:
                # add this record's values onto the running total, timestep by timestep
                for vi, v in enumerate(values):
                    inflows[vi] += v
            else:
                # subtract this record's values from the running total, timestep by timestep
                for vi, v in enumerate(values):
                    inflows[vi] -= v

    # if requested, prepend lookback values (copies of the first value) onto the front of the series
    if prepend_n > 0:
        # add first value onto front of record
        # Sometimes ResSim needs some lookback values, or whatever
        p_times = [times[i] for i in range(len(times))] # convert to list - annoying
        time_delta = times[1] -times[0]
        times = [p_times[0] - time_delta*i for i in range(prepend_n,0,-1)] + p_times
        values = [inflows[i] for i in range(len(inflows))] # convert to list - annoying
        inflows = [values[0]]*prepend_n + values
                    
    # Output record
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    #tsc.startTime = times[1]
    # set the output units based on whether values represent flow or something else
    if what=="flow":
        tsc.units = 'CFS'
    else:
        tsc.units = what
    tsc.type = tstype
    #tsc.endTime = times[-1]
    tsc.numberValues = len(inflows)
    #tsc.startHecTime = timewindow.getStartTime()
    #tsc.endHecTime = timewindow.getEndTime()
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    dssFm.close()
    dssFm_out.close()


def create_constant_dss_rec(currentAlt, timewindow, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS', fpart='ZEROS'):
    """Create and write a dss record with a constant in it for the given time windows.
       what={'flow','temp-water'}
       period={'1HOUR','1DAY'}

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        timewindow: The run time window object, used to determine
            the time span of the constant record (padded by 1 day
            on each end).
        output_dss_file (str): Path to the DSS file to write the
            constant record to.
        constant (float, optional): The constant value to fill the
            record with. Defaults to 0.0.
        what (str, optional): The type of parameter being created.
            One of 'flow', 'temp-water', 'gate', 'evap', or 'elev'
            (case-insensitive), which determines the units and
            parameter label used. Defaults to 'flow'.
        dss_type (str, optional): The DSS record type (e.g.
            'PER-AVER', 'INST-VAL'). Defaults to 'PER-AVER'.
        period (str, optional): The record's time interval. Must be
            '1HOUR' or '1DAY' (case-insensitive). Defaults to
            '1HOUR'.
        cpart (str, optional): The C-part (location) label to assign
            to the record. Defaults to 'ZEROS'.
        fpart (str, optional): The F-part (version) label to assign
            to the record. Defaults to 'ZEROS'.

    Returns:
        bool: True once the constant record has been written
            successfully. False if `what` or `period` is not a
            recognized value (in which case nothing is written).
    """

    # select the units and parameter label based on the requested data type
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

    # validate that the requested period is supported
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

    # generate a regular-interval series filled entirely with the constant value
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

def calculate_dewpoint(air_temp, relative_humidity):
    """
    Calculate Dew Point Temperature given the air temperature and relative humidity,
    using the algebraic inversion of the simplified August-Roche-Magnus approximation.
    
    Parameters:
        air_temp (float): Air temperature in C.
        relative_humidity (float): Relative Humidity in percent (0-100).
    
    Returns:
        float: Dew Point Temperature in C.
    """
    gamma = math.log(relative_humidity / 100.0) + (17.62 * air_temp) / (243.12 + air_temp)
    dewpoint = 243.12 * gamma / (17.62 - gamma)
    return dewpoint


def relhum_from_at_dp(met_dss_file, at_path, dp_path):
    """
    Compute relative humidity from air temperature and dew point
    records, and write the result as a new derived DSS record.

    Args:
        met_dss_file (str): Path to the DSS file containing the
            source air temperature and dew point records, and where
            the derived humidity record will be written.
        at_path (str): DSS path of the air temperature record.
        dp_path (str): DSS path of the dew point record. Must align
            timestep-by-timestep with at_path.

    Returns:
        None. Writes the derived relative humidity record (with
        '-DERIVED' appended to the F-part and parameter set to
        'RELHUM-FROM-AT-DP') to met_dss_file.
    """
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    dp_data = data_from_dss(met_dss_file, dp_path, None, None)
    # for each timestep, compute relative humidity from the paired air temp and dew point
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


def dp_from_at_relhum(met_dss_file, at_path, rh_path):
    """
    Compute dew point temperature from air temperature and relative
    humidity records, and write the result as a new derived DSS
    record.

    Args:
        met_dss_file (str): Path to the DSS file containing the
            source air temperature and relative humidity records,
            and where the derived dew point record will be written.
        at_path (str): DSS path of the air temperature record.
        rh_path (str): DSS path of the relative humidity record.
            Must align timestep-by-timestep with at_path.

    Returns:
        None. Writes the derived dew point record (with '-DERIVED'
        appended to the F-part and parameter set to
        'temp-dewpoint') to met_dss_file.
    """
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    rh_data = data_from_dss(met_dss_file, rh_path, None, None)
    # for each timestep, compute dew point from the paired air temp and relative humidity
    for i in range(tsc.numberValues):
        #print('AT:',tsc.values[i], 'RH:',rh_data[i])
        tsc.values[i] = calculate_dewpoint(tsc.values[i], rh_data[i])
    parts = tsc.fullName.split('/')
    parts[2] = parts[2][:5]
    parts[3] = 'temp-dewpoint'
    parts[6] = parts[6] + '-DERIVED'
    new_pathname = '/'.join(parts)
    tsc.fullName = new_pathname
    #tsc.units = '%' - units should be C...
    print('writing: ', new_pathname)
    dss.write(tsc)
    dss.close()

def check_start_and_end(values, times, startime, endtime):
    """
    Trim a values/times pair so that the data falls entirely within
    a given start/end time window, dropping any leading or trailing
    values outside that range.

    Args:
        values (list): The data values to trim.
        times (list): The corresponding raw HEC time values, in the
            same order as values.
        startime: The window start time, in the same units as
            times.
        endtime: The window end time, in the same units as times.

    Returns:
        tuple: (values, times), each trimmed to fall within
            [startime, endtime].
    """
    if times[0] < startime:  # if startdate is before the timewindow..
        print('start date ({0}) from DSS before timewindow ({1})..'.format(times[0], startime))
        st_offset = (startime - times[0]) / (times[1] - times[0])
        values = values[st_offset:]
        times = times[st_offset:]
    # if the data extends past the window, trim the trailing values
    if times[-1] > endtime:
        print('end date ({0}) from DSS after timewindow ({1})..'.format(times[-1], endtime))
        st_offset = (times[-1] - endtime) / (times[1] - times[0])
        values = values[:(len(times) - st_offset)]
        times = times[:(len(times) - st_offset)]
    return values, times


def replace_data(currentAlt, timewindow, pairs, dss_file, dss_outfile, months, standard_interval=None):
    """
    For each (base, alternate) record pair, replace the base
    record's values with the alternate record's values during a
    specified set of months, and write the merged result to a new
    output record.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        timewindow: The run time window object, used to determine
            the read window for both records in each pair.
        pairs (list of list): A list of [base_path, alt_path] DSS
            path pairs. For each pair, values from alt_path replace
            values from base_path during the given months.
        dss_file (str): Path to the DSS file containing both the
            base and alternate records.
        dss_outfile (str): Path to the DSS file to write the merged
            result to.
        months (list of int): The calendar months (1-12) during
            which the alternate record's values should replace the
            base record's values.
        standard_interval (str, optional): If provided, both the
            base and alternate series are standardized to this
            interval (via standardize_interval) before merging.
            Defaults to None (no standardization).

    Returns:
        None. Writes one merged record per pair in pairs to
        dss_outfile, with the F-part/A-part renamed to indicate the
        merge source.

    Raises:
        SystemExit: Exits the process (sys.exit(1)) if the base and
            alternate records in a pair have mismatched units or
            mismatched intervals.
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    # 01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    # for each base/alternate record pair, merge them and write the result
    for pair in pairs:
        dssFm = HecDss.open(dss_file)
        currentAlt.addComputeMessage('Replacing data for {0} with {1} during {2}'.format(pair[0], pair[1], months))
        base = dssFm.read(pair[0], starttime_str, endtime_str, False)
        # optionally standardize the base series to a common interval before merging
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
        # optionally standardize the alternate series to the same common interval
        if standard_interval is not None:
            alt = standardize_interval(alt,standard_interval)
        alt_data = alt.getData()
        alt_values = alt_data.values
        alt_hectimes = alt_data.times
        alt_units = alt_data.units
        alt_interval = alt_data.interval
        alt_values, alt_hectimes = check_start_and_end(alt_values, alt_hectimes, starttime_hectime, endtime_hectime)

        dssFm.close()

        # the two records must share units to be merged meaningfully
        if base_units != alt_units:
            currentAlt.addComputeMessage('Units do not match for {0} and {1}, skipping'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)
        if base_interval != alt_interval:
            currentAlt.addComputeMessage('Intervals do not match for {0} and {1}, changing interval...'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)

        # for each timestep in the base record, if its month is in the replacement list,
        # swap in the alternate record's value instead
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

def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    """
    Apply a fixed temperature lapse (offset) to an air temperature
    record, to approximate temperature change with elevation, and
    write the adjusted record under a new F-part.

    NOTE: This is a duplicate definition of airtemp_lapse, identical
    to the one defined earlier in this module. Because Python/Jython
    executes definitions top-to-bottom, this later definition
    silently replaces the earlier one - any call to airtemp_lapse
    uses this version. Consider removing one of the two copies.

    Args:
        dss_file (str): Path to the DSS file containing the source
            air temperature record.
        dss_rec (str): DSS path of the source air temperature
            record.
        lapse_in_C (float): The temperature lapse/offset to apply,
            in degrees Celsius. If the source record's units are
            Fahrenheit, this is converted to an equivalent
            Fahrenheit offset before being applied.
        dss_outfile (str): Path to the DSS file to write the
            adjusted record to.
        f_part (str): The F-part to assign to the output record.

    Returns:
        None. Writes the lapse-adjusted record to dss_outfile.
    """
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    lapse = lapse_in_C
    # if the source series is in Fahrenheit, convert the lapse amount to an equivalent Fahrenheit offset
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

def preprend_first_value_on_ts(dss_file,dss_rec,prepend_n):
    '''Sometimes ResSim needs some lookback values, or whatever
    
    Be careful that the first record is where you want to start - sometimes these things change
    '''
    dss = HecDss.open(dss_file)
    tsc = dss.get(dss_rec,True)

    time_delta = tsc.times[1] - tsc.times[0]
    times = [tsc.times[i] for i in range(len(tsc.times))] # convert to list - annoying
    tsc.times = [times[0] - time_delta*i for i in range(prepend_n,0,-1)] + times
    tsc.startTime = tsc.times[0]
    values = [tsc.values[i] for i in range(len(tsc.values))] # convert to list - annoying
    tsc.values = [values[0]]*prepend_n + values
    tsc.numberValues = len(tsc.values)

    dss.put(tsc)
    dss.close()

def postprend_last_value_on_ts(dss_file,dss_rec,postpend_n):
    """Sometimes ResSim needs some lookback values, or whatever
    has to be a 'regular' record

    Args:
        dss_file (str): Path to the DSS file containing the source
            record, which is also modified in place.
        dss_rec (str): DSS path of the record to modify. Must be a
            regular-interval record.
        postpend_n (int): Number of copies of the last value (and
            corresponding later timestamps) to append onto the end
            of the record.

    Returns:
        None. Overwrites dss_rec in dss_file with the extended
        series.
    """
    dss = HecDss.open(dss_file)
    tsc = dss.get(dss_rec,True)

    time_delta = tsc.times[1] - tsc.times[0]
    times = [tsc.times[i] for i in range(len(tsc.times))] # convert to list - annoying
    tsc.times = times + [times[-1] + time_delta*i for i in range(1,postpend_n+1)]
    #tsc.startTime = tsc.times[0]
    values = [tsc.values[i] for i in range(len(tsc.values))] # convert to list - annoying
    tsc.values = values + [values[-1]]*postpend_n
    tsc.numberValues = len(tsc.values)

    dss.put(tsc)
    dss.close()   

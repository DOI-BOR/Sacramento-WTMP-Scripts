"""
Simple DSS + Interpolation Utilities

This module contains helper routines used within HEC-WAT/ResSim/W2 scripting
for the Sacramento WTMP workflow. It includes:

- DSS record utilities (copying, reading, resampling, unit fixes, flow math).
- Time conversion helpers between HEC `HecTime`/Java `Date` and Python `datetime`.
- Lightweight linear interpolation for Python `datetime` series.
- Convenience functions for path manipulation and alternative-based resolution.


Notes
-----

- Several helpers reference HEC Java APIs accessed via Jython. These calls
  assume the HEC-WAT environment is available at runtime.

"""

# --- Imports -----------------------------------------------------------------

from hec.heclib.dss import HecDss  # HEC-DSS file manager: open/read/write/rename records
from hec.io import DSSIdentifier  # DSS pathname identifier helper (A-F parts parsing/building)
from hec.io import TimeSeriesContainer  # Java container for DSS time series (times/values/units/type)
from rma.util.RMAConst import MISSING_DOUBLE  # RMA/HEC sentinel for missing double values
from hec.hecmath import HecMathException  # Exception raised by HEC Math operations (e.g., transform errors)
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE  # Sentinel for undefined doubles in Heclib utilities
import hec.hecmath.TimeSeriesMath as tsmath  # Time series math operations (generate, transform, add, shift)
from com.rma.model import Project  # Access to current project context (directories, settings) from Jython
import os, shutil, copy, sys, math  # Python stdlib: filesystem ops, copying, process control, math utilities
from java.util import Vector, Date  # Java utilities: Vector (mutable list), Date (epoch-based date-time)

import datetime  # Python datetime for conversions and arithmetic
from hec.heclib.util import HecTime  # HEC time representation used by DSS series



def _timedelta_to_seconds(td):
    """
    Convert a `datetime.timedelta` to seconds (float).

    Parameters
    ----------
    td : datetime.timedelta
        Time delta to convert.

    Returns
    -------
    float
        Total seconds represented by `td`.

    Notes
    -----
    Provides Python 2 compatibility where `timedelta.total_seconds()` may
    not be present or consistent. Uses explicit decomposition into days,
    seconds, and microseconds.
    """
    
    # 86400 seconds in a day (24 * 60 * 60)
    # 1,000,000 microseconds in a second
    return td.days * 86400.0 + td.seconds + td.microseconds / 1000000.0


def linear_interp_datetime(time_data, value_data, query_times):
    """
    Perform linear interpolation on datetime-indexed series.

    Parameters
    ----------
    time_data : list of datetime.datetime
        Ascending (ordered) timestamps for known data.
    value_data : list of float
        Values corresponding to `time_data`.
    query_times : list of datetime.datetime
        Timestamps at which to interpolate values.

    Returns
    -------
    list of float
        Interpolated (or edge extrapolated) values for `query_times`.

    Notes
    -----
    - Extrapolates using the nearest endpoint value when a query falls outside
      the known time range.
    - Handles duplicate adjacent timestamps by returning the lower endpoint value
      for that segment.
    """
    
    # 1. Basic Validation
    if len(time_data) != len(value_data) or len(time_data) < 2:
        print "Error: time_data and value_data must have the same length (and at least 2 points)."
        return []

    # 2. Normalize Times to Numerical Values (Seconds from the first timestamp)
    # Using the first time point (t0) as the numerical zero reference (x_data[0] = 0.0)
    t0 = time_data[0]
    x_data = []
    
    for t in time_data:
        x_data.append(_timedelta_to_seconds(t - t0))  # convert to seconds relative to t0

    x_query = []
    for t in query_times:
        x_query.append(_timedelta_to_seconds(t - t0))  # same normalization for queries

    # 3. Perform Interpolation / Extrapolation
    results = []
    n = len(x_data)

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
            y = y0  # degenerate segment: hold the left value
        
        else:
            # Linear interpolation formula: y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)
            y = y0 + (y1 - y0) * (x - x0) / dx

        # Append interpolated value
        results.append(y) 

    # Return final list of interpolated values
    return results  


def copy_dss_ts(dss_rec, new_fpart=None, new_dss_rec=None,
                dss_file_path=None, dss_file_handle=None, dss_file_alt_outpath=None, checkMakeCelsius=False):
    """
    Copy a DSS time series to a new record, optionally adjusting F-part and units.

    Parameters
    ----------
    dss_rec : str
        Source DSS pathname (e.g., "/A/B/C/D/E/F").
    new_fpart : str, optional
        New F-part to use (ignored if `new_dss_rec` is provided).
    new_dss_rec : str, optional
        Explicit destination DSS pathname (overrides `new_fpart`).
    dss_file_path : str, optional
        Path to DSS file to open (if a handle is not supplied).
    dss_file_handle : hec.heclib.dss.HecDss, optional
        Already-open DSS file handle to use for the read/write.
    dss_file_alt_outpath : str, optional
        If provided, writes the copied series to this alternate DSS file.
    checkMakeCelsius : bool, optional
        If True and the input units are Fahrenheit, convert to Celsius.

    Returns
    -------
    None

    Notes
    -----
    - Normalizes certain unit strings: 'DEGC' -> 'C', 'DEGF' -> 'F'.
    - If `checkMakeCelsius` is True, converts values from °F to °C and updates units.
    - Destination record name is constructed via F-part replacement or taken directly
      from `new_dss_rec`.
    """
    
    # Error check inputs - there are flexible way to copy record
    if dss_file_path is None and dss_file_handle is None:
        raise ValueError('copy_dss_rec_to_new_fpart: you must supply either a valid dss_file_path OR dss_file_handle')
    
    if new_fpart is None and new_dss_rec is None:
        raise ValueError('copy_dss_rec_to_new_fpart: you must supply either a new_fpart OR new_dss_rec')

    # Open dss get record tsc
    if dss_file_handle is not None:
        # reuse provided file handle
        dss_fm = dss_file_handle  

    else:
        # open by path if handle not supplied
        dss_fm = HecDss.open(dss_file_path)  

    # read full record
    tsc = dss_fm.get(dss_rec, True)  

    # new dss rec name
    if new_dss_rec is not None:
        # Explicit destination takes priority 
        dss_rec_out = new_dss_rec  
    
    else:
        # No pathname provided. Use default.
        rec_parts = tsc.fullName.split('/')  # tokenize pathname
        rec_parts[6] = new_fpart  # update F-part
        dss_rec_out = '/'.join(rec_parts)  # rebuild pathname

    # fix some terrible units along the way
    if tsc.units.lower() == 'degc':
        tsc.units = 'C'  # normalize to 'C'
    
    elif tsc.units.lower() == 'degf':
        tsc.units = 'F'  # normalize to 'F'

    # Convert from F to C if needed
    if tsc.units.lower() == 'f' and checkMakeCelsius:
        T_values = tsc.values  # get values
        
        for i, TT in enumerate(T_values):
            T_values[i] = (TT - 32.0) * 5.0 / 9.0  # convert °F -> °C in-place
        
        tsc.units = 'C'  # update units
        tsc.values = T_values  # reassign (explicit, per original)

    # Write back to the file
    tsc.fullName = dss_rec_out  # set output fullName
    if dss_file_alt_outpath is not None:
        dss_fm_alt = HecDss.open(dss_file_alt_outpath)  # open alternate output DSS
        dss_fm_alt.put(tsc)  # write to alt file
        dss_fm_alt.close()  # close alt
    
    else:
        dss_fm.put(tsc)  # write to same file

    # Close the DSS file if it was opened
    if dss_file_handle is None:
        dss_fm.close()  


def jday_from_tsc(tsc):
    """
    Compute decimal day-of-year for each time in a `TimeSeriesContainer`.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        DSS container with HEC times accessible via `getHecTime(i)`.

    Returns
    -------
    list of float
        Decimal DOY for each timestamp in `tsc`.
    """
    
    # Convert to Python datetimes
    dtt = hectime_to_datetime(tsc)  
    
    # Compute and return fractional DOY
    return [decimal_doy(dt) for dt in dtt]  


def decimal_doy(dt):
    """
    Convert a Python `datetime` to decimal day-of-year (DOY).

    Parameters
    ----------
    dt : datetime.datetime
        Timestamp to convert.

    Returns
    -------
    float
        Day-of-year with fractional component from time-of-day.
    """
    
    # integer day-of-year from struct_time
    doy = dt.timetuple().tm_yday  
    
    # Compute the current fraction of a day
    fractional_day = (dt.hour / 24.0) + (dt.minute / 1440.0) + (dt.second / 86400.0) + (dt.microsecond / 86400000000.0)
    
    # Sum and return integer and fractional parts
    return doy + fractional_day  


def organizeLocations(curAlt, location_objs, loc_names, return_dss_paths=False):
    """
    Order location objects or paths according to the provided name list.

    Parameters
    ----------
    curAlt : object
        Alternative context with `loadTimeSeries()` and logging.
    location_objs : list
        List of location objects supporting `getName()` / DSS linkage APIs.
    loc_names : list of str
        Desired ordering of location names.
    return_dss_paths : bool, optional
        If True, return string DSS paths; else return location objects.

    Returns
    -------
    list
        Ordered list of location objects or DSS paths.
    """
    
    # Define a list to hold the output
    locations_list = []  # accumulator for outputs
    
    # Loop and process over each location 
    for name in loc_names:
        # Get the index of the location
        i_loc = findLocationOrder(curAlt, location_objs, name)
        
        # Extract the DSS paths if requested
        if return_dss_paths:
            # Get the location information
            lo1 = location_objs[i_loc]  
            
            # Handle the request differently if the data is linked to the previous model
            if lo1.isLinkedToPreviousModel():
                # load from previous model linkage
                tspath = str(curAlt.loadTimeSeries(lo1))  
                
                # adjust F-part
                tspath = fixInputLocationFpart(curAlt, tspath)  
            
            else:
                # Data is not linked. Use a direct path.
                tspath = str(lo1.getDssPath())  
            
            # Append into the list
            locations_list.append(tspath)
        
        else:
            # Use the location information directly
            locations_list.append(location_objs[i_loc])
    
    # Return the list to the calling function
    return locations_list


def organizeLocationsPaired(curAlt, location_objs, loc_names_paired, return_dss_paths=False):
    """
    Apply `organizeLocations` to groups (pairs) of location names.

    Parameters
    ----------
    curAlt : object
        Alternative context for DSS path resolution.
    location_objs : list
        Location objects available for lookup.
    loc_names_paired : list of list of str
        Grouped names (e.g., [[A, B], [C, D]]) to process together.
    return_dss_paths : bool, optional
        If True, return DSS paths; else return location objects.

    Returns
    -------
    list of list
        Grouped results preserving input structure.
    """
    
    return [organizeLocations(curAlt, location_objs, pn, return_dss_paths) for pn in loc_names_paired]


def findLocationOrder(curAlt, location_objs, name):
    """
    Find the index of a location object by its name.

    Parameters
    ----------
    curAlt : object
        Alternative context for log messages on failure.
    location_objs : list
        Iterable of objects with `getName()` method.
    name : str
        Target location name to match.

    Returns
    -------
    int
        Index of the matched location.

    Notes
    -----
    Exits the script (`sys.exit(1)`) if the name is not found.
    """
    
    # Loop over each of the provided locations
    for i, loc in enumerate(location_objs):
        # Compare the name to the location name
        if name == loc.getName():
            # return match position
            return i  
    
    # if we make it here, our input/output location name was not found
    curAlt.addComputeMessage("Scripting - Location name not found: " + name)  # log
    
    # abort to indicate error
    sys.exit(1)  


def first_value(dss_file, dss_rec, start_str=None, end_str=None):
    """
    Return the first value from a DSS record (optionally within a time window).

    Parameters
    ----------
    dss_file : str
        Path to the DSS file.
    dss_rec : str
        Full DSS pathname to read.
    start_str : str or None, optional
        Window start in HEC format ("DDMonYYYY HHMM"). If None, read full record.
    end_str : str or None, optional
        Window end in HEC format ("DDMonYYYY HHMM"). If None, read full record.

    Returns
    -------
    float
        The first value in the record or in the specified window.
    """
    
    # open input DSS
    dssFm = HecDss.open(dss_file)  
    
    # Handle the read based on the sepcified times
    if start_str is None and end_str is None:
        # read full record
        tsc = dssFm.get(dss_rec, True)  
    
    else:
        # Start and end times are provided. Read the partial series.
        tsc = dssFm.read(dss_rec, start_str, end_str, False).getData()  
    
    # close DSS
    dssFm.close()  
    
    # Return first value
    return tsc.values[0]  


def standardize_interval(tsm, interval, makePerAver=True):
    """
    Ensure a `TimeSeriesMath` series is at the requested interval (transform if needed).

    Parameters
    ----------
    tsm : hec.hecmath.TimeSeriesMath
        Math wrapper for the time series to standardize.
    interval : {'1hour', '1day', '1week'}
        Target interval name used by HEC transforms.
    makePerAver : bool, optional
        If True, set type to 'PER-AVER' prior to transformation.

    Returns
    -------
    hec.hecmath.TimeSeriesMath
        Original `tsm` if already at interval; else transformed series.

    Raises
    ------
    SystemExit
        If an unsupported `interval` string is provided.
    """
    
    # Get the underlying data from the container
    tsc = tsm.getData()  
    
    # Determine the time interval
    if interval.lower() == '1hour':
        intint = 60  # minutes per hour
    
    elif interval.lower() == '1day':
        intint = 1440  # minutes per day
    
    elif interval.lower() == '1week':
        intint = 10800  # minutes per week
    
    else:
        # Invalid condition. Error and exit
        sys.exit(-1)

    # Check that the time series interval is the same as the series
    if tsc.interval != intint:
        if makePerAver:
            tsm.setType('PER-AVER')  # set series type on math object
        return tsm.transformTimeSeries(interval, "", "AVE")  # perform transform
    
    else:
        return tsm  # already standardized


def get_sanitized_record_list(dss_file_path):
    """
    Produce a de-duplicated list of DSS pathnames with date (E-part) cleared.

    The DSS library seems to return lists of paths with dates in them (getPathnameList()), and some of those
    dates don't even exist in the file or cannot be read and throw an error.  As of Jan 2024,
    this is an orregular problem, and the manual soluation is throwing away the DSS file but
    in many cases that is problematic. So, here we filter dates and check for duplicates.

    Parameters
    ----------
    dss_file_path : str
        Path to the DSS file.

    Returns
    -------
    list of str
        Sanitized pathnames with E-part blank and duplicates removed.

    Notes
    -----
    Workaround for occasional stale/unreadable date entries returned by
    `getPathnameList()`; stripping E-part avoids read errors downstream.
    """
    
    # Open DSS
    dss = HecDss.open(dss_file_path)  
    
    # Get all pathnames
    recs = dss.getPathnameList()  
    
    # Close the DSS file
    dss.close()  
    
    # Define a list to hold the values
    sanitized_recs = []
    
    # Loop over each of the records
    for r in recs:
        rec_tokens = r.split('/')  # split into parts
        rec_tokens[4] = ''  # erase date string, if it exists (E-part)
        r_sanitized = '/'.join(rec_tokens)  # rebuild pathname
    
        if not r_sanitized in sanitized_recs:
            sanitized_recs.append(r_sanitized)  # dedupe
    
    # Return cleaned list to the calling function
    return sanitized_recs  


def hectimes_from_tsm(tsm):
    """
    Build a list of `HecTime` objects from a `TimeSeriesMath` container.

    Parameters
    ----------
    tsm : hec.hecmath.TimeSeriesMath
        Math wrapper whose underlying container has `.times`.

    Returns
    -------
    list of hec.heclib.util.HecTime
        HEC time objects corresponding to each entry in the series.
    """
    
    # Get numeric HEC time values
    times = tsm.getContainer().times  
    htimes = []  # accumulator
    
    # Loop over timeseries and process each value
    for i in range(tsm.getContainer().numberValues):
        htimes.append(HecTime())  # construct new HecTime per value
        htimes[-1].set(times[i])  # assign numeric time
    
    # Return list of HecTime
    return htimes  


def hectimes_from_tsc(tsc):
    """
    Build `HecTime` objects from a `TimeSeriesContainer`.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        Container with `.times` (HEC numeric) and `.numberValues`.

    Returns
    -------
    list of hec.heclib.util.HecTime
        HEC time objects, one per series element.
    """
    
    # Create a list to hold data
    htimes = [] 
    
    # Loop over each timeseries value and process
    for i in range(tsc.numberValues):
        htimes.append(HecTime())  # new HecTime
        htimes[-1].set(tsc.times[i])  # assign from container times
    
    # Return the output list of values
    return htimes


def shift_pit_river_time(input_dss_file, dss_rec, output_dss_file, out_rec, start_date=None, end_date=None):
    """
    Shift a Pit River time series by -12 hours and write to a new record.

    Parameters
    ----------
    input_dss_file : str
        Source DSS file path.
    dss_rec : str
        Source record pathname.
    output_dss_file : str
        Destination DSS file path.
    out_rec : str
        Destination record pathname (must differ from source).
    start_date : str or None, optional
        HEC-format start time for windowed read.
    end_date : str or None, optional
        HEC-format end time for windowed read.

    Returns
    -------
    None

    Notes
    -----
    Important: do not set `dss_rec == out_rec` to avoid repeatedly re-shifting on subsequent runs.
    """
   
    # Read the series of values
    tsm_pit = dss_read_ts_safe(input_dss_file, dss_rec, start_date=start_date, end_date=end_date, returnTSM=True)

    # Shift the values by -12 hours
    tsm_pit = tsm_pit.shiftInTime("-12Hour")
    
    # Request the data and path
    tsc_pit = tsm_pit.getData()  # unwrap container
    tsc_pit.fullName = out_rec  # set destination pathname

    # Output the shifted series
    dss_out = HecDss.open(output_dss_file)  # open destination
    dss_out.put(tsc_pit)  # write shifted series
    dss_out.close()  # close output DSS


def shift_ts_time(input_dss_file, dss_rec, output_dss_file, out_rec, shift_str, start_date=None, end_date=None):
    """
    Shift a time series by a HEC-recognized offset string and write to a new record.

    Parameters
    ----------
    input_dss_file : str
        Source DSS file path.
    dss_rec : str
        Source record pathname.
    output_dss_file : str
        Destination DSS file path.
    out_rec : str
        Destination record pathname (must differ from source).
    shift_str : str
        HEC time-shift string (e.g., "-12Hour", "5Day").
    start_date : str or None, optional
        Windowed read start time in HEC format.
    end_date : str or None, optional
        Windowed read end time in HEC format.

    Returns
    -------
    None
    
    Notes
    -----
    Important to not have dss_rec==out_rec, otherwise you will continually shift the record in time
    whenever the calling script runs...

    shift_str is HEC-known string, with a possible negative on the front. E.g., '-12Hour','5Day'
    """

    # Read the timeseries
    tsm = dss_read_ts_safe(input_dss_file, dss_rec, start_date=start_date, end_date=end_date, returnTSM=True)

    # Shift the values
    tsm = tsm.shiftInTime(shift_str)  # perform time shift
    
    # Request the data and the path
    tsc = tsm.getData()  # unwrap to container
    tsc.fullName = out_rec  # set destination pathname

    # Output the shifted values
    dss_out = HecDss.open(output_dss_file)  # open destination DSS
    dss_out.put(tsc)  # write shifted series
    dss_out.close()  # close output


def dss_read_ts_safe(dssFilePath, dssRec, start_date=None, end_date=None, returnTSM=False,
                     returnPydatetimes=False, debug=False):
    """
    Read a DSS time series record safely, optionally with a time window.

    Parameters
    ----------
    dssFilePath : str
        Path to the DSS file.
    dssRec : str
        Full DSS pathname of the record.
    start_date : str or None, optional
        Window start in HEC format. If None, read entire record.
    end_date : str or None, optional
        Window end in HEC format. If None, read entire record.
    returnTSM : bool, optional
        If True, return `TimeSeriesMath`; else return `TimeSeriesContainer`.
    returnPydatetimes : bool, optional
        (Unused in original) Placeholder for returning Python datetimes.
    debug : bool, optional
        If True, print diagnostics.

    Returns
    -------
    hec.io.TimeSeriesContainer or hec.hecmath.TimeSeriesMath
        Data container or math wrapper depending on `returnTSM`.
    """
    '''A read function that is date-flexible, and ensures the whole time range is returned if dates
    are None.'''

    # This method is going to throw an error if the file doesn't exist
    dss = HecDss.open(dssFilePath, True)  # open file

    # recordExists() seems to check for the individual database "pages" or something, with a date, or
    # date range, required?  TODO: figure out how to check if a date-agnostic record exists.  Or, don't,
    # since HECLIB produces a pretty nice error if you try to open a non-existant record.
    if start_date is None and end_date is None:
        # full record read
        tsc = dss.get(dssRec, True)  
        
        # close DSS
        dss.close()  
        
        # Print debug information if enabled
        if debug:
            print('Reading DSS in script...')
            print('    file: ' + dssFilePath)
            print('    record: ' + dssRec)
        
        # Wrap in math container if requested
        if returnTSM:
            return tsmath(tsc)  
        
        else:
            return tsc  # return data container
    
    elif start_date is not None and end_date is not None:
        # Windowed read 
        tsm = dss.read(dssRec, start_date, end_date, False)  
        
        # close DSS
        dss.close()  
        
        # Print debug information if enabled
        if debug:
            print('Reading DSS in script between ' + start_date + ' and ', + end_date)
            print('    file: ' + dssFilePath)
            print('    record: ' + dssRec)
        
        # Wrap in matho container if requested
        if returnTSM:
            return tsm  
        
        else:
            return tsm.getData()  # unwrap to data container


def data_from_dss(dss_file, dss_rec, starttime_str, endtime_str):
    """
    Read values from a DSS record (optionally within a time window).

    Parameters
    ----------
    dss_file : str
        DSS file path.
    dss_rec : str
        Full DSS pathname of the record.
    starttime_str : str or None
        Window start (HEC format) or None for full record.
    endtime_str : str or None
        Window end (HEC format) or None for full record.

    Returns
    -------
    list
        Values array for the record or the windowed subset.
    """
    
    # Open input DSS
    dssFm = HecDss.open(dss_file)  
    
    # Determine what portion of the timeseries to provide
    if starttime_str is None and endtime_str is None:
        # Full record
        tsc = dssFm.get(dss_rec, True)  
    
    else:
        # Read a portion of the timeseries
        tsc = dssFm.read(dss_rec, starttime_str, endtime_str, False).getData()  
    
    # close DSS
    dssFm.close()  
    
    # Return values to the calling function
    return tsc.values  


def hectime_to_datetime(tsc):
    """
    **Deprecated**: Convert HEC times in a `TimeSeriesContainer` to Python `datetime`.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        Container with HEC times accessible via `getHecTime(i)`.

    Returns
    -------
    list of datetime.datetime
        Python datetimes corresponding to each container time.

    Notes
    -----
    This routine can round hours incorrectly and applies a timezone offset
    returning GMT. Prefer `datetimes_from_tsc()` with `hecTime_to_datetime()`
    for safer conversion in new code.
    """
    """
    Deprecated! This routine, and the javatime processing it depends on, is capable of rounding dates 
    to the wrong hour.  It also adds a timezone offset and returns GMT, which is not usually wanted.
    
    TODO: replace all instances of this with datetimes_from_tsc(), and rework/remove any extra processing
    previously needed for timezones, bad first/last values, etc.
    """
    
    # Get a list of hectimes
    hectimes = hectimes_from_tsc(tsc)
    
    # Define a list to hold values
    dtt = []
    
    # Loop and process each value
    for j in range(tsc.numberValues):
        # Assuming hectime can be converted to Java Date or has method to get the equivalent
        java_date = tsc.getHecTime(j).getJavaDate(0)  # Java date with timezone index

        # Convert Java Date to Python datetime
        timestamp = (java_date.getTime() / 1000)  # milliseconds to seconds
        d = datetime.datetime.fromtimestamp(timestamp)  # local conversion (unused local var)
        
        # Append converted time
        dtt.append(datetime.datetime.fromtimestamp(timestamp))  

    # Return list of datetimes
    return dtt  


def hecTime_to_datetime(hectimes):
    """
    Convert a list of `HecTime` objects to Python `datetime`.

    Parameters
    ----------
    hectimes : list of hec.heclib.util.HecTime
        HEC times to convert.

    Returns
    -------
    list of datetime.datetime
        Converted Python datetimes via string parsing helper.
    """
    
    return [hec_str_time_to_dt(h.toString(), longStr=True) for h in hectimes]


def datetimes_from_tsc(tsc):
    """
    Convert a `TimeSeriesContainer` to Python datetimes via `HecTime` conversion.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        Source container to convert.

    Returns
    -------
    list of datetime.datetime
        Converted Python datetimes.
    """
    
    return hecTime_to_datetime(hectimes_from_tsc(tsc)) 


def hec_str_time_to_dt(hec_str_time, longStr=False):
    """
    Convert a HEC time string to Python `datetime`, handling "24:00"/"2400" rollover.

    Parameters
    ----------
    hec_str_time : str
        HEC-formatted time string. If `longStr` is True, uses verbose format
        like "DD Month YYYY, HH:MM"; otherwise "DDMonYYYY HHMM".
    longStr : bool, optional
        Flag to select the verbose format with `24:00` rollover.

    Returns
    -------
    datetime.datetime
        Converted Python `datetime`.

    Notes
    -----
    If the string ends with "2400" or "24:00" (depending on format), it is
    converted to "0000"/"00:00" and one day is added.
    """

    # Define a flag to indicate rollover
    add_day = False
    
    # Definte the format of the time string
    if longStr:
        dt_format = "%d %B %Y, %H:%M"  # verbose format
        check24 = '24:00'
        check00 = '00:00'
    
    else:
        dt_format = '%d%b%Y %H%M'  # compact HEC format
        check24 = '2400'
        check00 = '0000'

    # Determine the length of the hour string
    checkLen = len(check24) * (-1) 
    
    # Determine if the hour time needs to be adjusted
    if hec_str_time.endswith(check24):
        my_hec_str_time = hec_str_time[:checkLen] + check00  # replace ending time
        add_day = True  # add one day after parse
    
    else:
        my_hec_str_time = hec_str_time  # unchanged

    # Parse the time with the correct format
    dt = datetime.datetime.strptime(my_hec_str_time, dt_format)
    
    # Add the day for the rollover
    if add_day:
        dt = dt + datetime.timedelta(days=1)
    
    # Return the jython timestamp
    return dt 


def fixInputLocationFpart(currentAlternative, tspath):
    """
    Adjust the F-part of an input location DSS path based on the Alternative's input F-part.

    Parameters
    ----------
    currentAlternative : object
        Alternative context providing `getInputFPart()`.
    tspath : str
        Full DSS pathname to adjust.

    Returns
    -------
    str
        Adjusted pathname with updated F-part preserving suffix.
    """
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])  # prefix sans last token
    tspath = tspath.split('/')  # tokenize pathname
    fpart = tspath[6]  # F-part index
    fpart_split = fpart.split(':')  # split by ':'
    new_fpart = new_fpart_start + ':' + fpart_split[-1]  # reconstructed F-part
    tspath[6] = new_fpart  # assign
    tspath = '/'.join(tspath)  # rebuild full path
    return tspath  # return adjusted pathname


def appendAPart(current_path, ApartAppend):
    """
    Append a suffix to the A-part of a DSS pathname.

    Parameters
    ----------
    current_path : str
        Existing DSS pathname string.
    ApartAppend : str
        Text to append (with underscore) to the A-part.

    Returns
    -------
    str
        Updated DSS pathname with modified A-part.

    """

    # Split the path
    tspath = tspath.split('/')

    # Extract the A part of the path
    Apart = tspath[1]  
    
    # Determine how to handle the A part
    if len(Apart) == 0:
        new_Apart = ApartAppend  # replace empty A-part
    
    else:
        new_Apart = Apart + '_' + ApartAppend  # append with underscore
    
    # Assign new A-part
    tspath[1] = new_Apart

    # Rebuild pathname
    tspath = '/'.join(tspath)  
    
    # Return updated path
    return tspath  


def getDataLocationDSSInfo(location, currentAlternative, computeOptions):
    """
    Resolve DSS path and file location for a model location, considering linkage.

    Parameters
    ----------
    location : object
        Location object supporting linkage queries.
    currentAlternative : object
        Alternative context for `loadTimeSeries()` and F-part adjustments.
    computeOptions : object
        Provides compute-time DSS filename.

    Returns
    -------
    (str, str)
        Tuple of (tspath, dsspath) for inputs/outputs.
    """
    
    # Determine if the data is linked to the previous model
    if location.isLinkedToPreviousModel():
        # Data is linked. Pull information from the previous model.
        tspath = str(currentAlternative.loadTimeSeries(location))  # path from previous model
        tspath = fixInputLocationFpart(currentAlternative, tspath)  # adjust F-part
        dsspath = computeOptions.getDssFilename()  # compute DSS file
    
    else:
        # Data is not linked. Create the direct path.
        tspath = location.getLinkedToLocation().getDssPath()  # direct path
        rundir = Project.getCurrentProject().getProjectDirectory()  # project directory
        dsspath = location.getLinkedToLocation().get_dssFile()  # relative DSS path
        dsspath = os.path.join(rundir, dsspath)  # absolute DSS path
    
    # Return the path information to the calling function
    return tspath, dsspath 


def strip_templateID_and_rename_records(dssFilePath, currentAlt):
    """
    Strip the first 4 characters from the F-part (when it contains '-') for all records, then rename.

    Parameters
    ----------
    dssFilePath : str
        Path to the DSS file to modify.
    currentAlt : object
        Alternative context for logging messages.

    Returns
    -------
    None

    Notes
    -----
    Creates a `.bak` copy, builds new pathnames in a `Vector`, and performs
    a bulk rename via `HecDss.renameRecords()`. If the F-part does not contain
    a '-', returns early (mirroring original code behavior).
    """
    
    # make copy of dss file
    shutil.copyfile(dssFilePath, dssFilePath + '.bak')  # backup original file

    # Rename all records, stripping first 4 chars from f-part
    dss = HecDss.open(dssFilePath)  # open DSS
    rec_names = dss.getPathnameList()  # enumerate pathnames
    new_rec_names = Vector()  # Java Vector for new names
    
    # Loop and process each record
    for i, r in enumerate(rec_names):
        # split into components
        parts = r.split('/')  
        
        # early exit if F-part does not contain '-'
        if not '-' in parts[-2]:
            return  
        
        # strip first 4 chars from F-part
        parts[-2] = parts[-2][4:]  
        
        # collect new pathname
        new_rec_names.add('/'.join(parts))  
        
        # log rename
        currentAlt.addComputeMessage('Fixing path: ' + r + ' --&gt; ' + new_rec_names[-1])  
    
    # Rename all records concurrently
    dss.renameRecords(rec_names, new_rec_names)  
    
    # Close the DSS file
    dss.close()  


def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    """
    Sum multiple DSS series over a window and write the result to `output_path`.

    Parameters
    ----------
    currentAlt : object
        Alternative context for logging.
    dssFile : str
        DSS file containing inputs and destination.
    timewindow : object
        Provides HEC-format start/end strings + HEC times.
    input_data : list of str
        DSS pathnames to read and sum.
    output_path : str
        Destination DSS pathname for the summed series.

    Returns
    -------
    int
        Always returns 0 after writing output.
    """
    
    # Get the temporal properties for the read
    starttime_str = timewindow.getStartTimeString()  # window start string
    endtime_str = timewindow.getEndTimeString()  # window end string
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))  # log window
    
    # Open the DSS file
    dssFm = HecDss.open(dssFile)
    
    # Create the data holder
    output_data = []
    
    # Loop and proccess each of the records
    for dsspath in input_data:
        print('reading', str(dsspath))  # diagnostic
        # Read the time series
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)  
        
        # Get the container and the data
        ts = ts.getData()  
        values = ts.values
        
        # Get the times and metadata
        times = ts.times  # HEC times
        units = ts.units  # units
        dsstype = ts.type  # series type
        
        # Append the data into the list
        if len(output_data) == 0:
            output_data = values  # seed with first series
        else:
            for vi, val in enumerate(values):
                output_data[vi] += val  # element-wise sum

    # Output the timeseries
    tsc = TimeSeriesContainer()  # prepare output container
    tsc.times = times  # assign times
    tsc.fullName = output_path  # destination path
    tsc.values = output_data  # summed values
    tsc.startTime = times[0]  # start time index
    tsc.units = units  # inherit units
    tsc.type = dsstype  # inherit type
    tsc.endTime = times[-1]  # end time index
    tsc.numberValues = len(output_data)  # series length
    tsc.startHecTime = timewindow.getStartTime()  # HEC start
    tsc.endHecTime = timewindow.getEndTime()  # HEC end
    dssFm.write(tsc)  # write output
    dssFm.close()  # close DSS
    
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))  # log count
    
    # Return that the write was sucessful
    return 0


def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod, lookback_1mon=False, pad_start_days=0):
    """
    Resample a DSS series to a new period (upsample or downsample), with optional padding/lookback.

    Parameters
    ----------
    inputDSSFile : str
        Input DSS file path.
    inputRec : str
        Input record pathname.
    timewindow : object or None
        If provided, supplies HEC-format start/end strings to constrain reads.
    outputDSSFile : str
        Output DSS file path.
    newPeriod : str
        Target period string for transform (e.g., '1HOUR', '1DAY').
    lookback_1mon : bool, optional
        If True, attempt a one-month lookback for the start time.
    pad_start_days : int, optional
        Additional days to subtract from the start for padding.

    Returns
    -------
    None

    Notes
    -----
    Coerces start to "0000" and end to "2400" for complete-day transforms. Can upsample an even period DSS timeseries, 
    e.g. go from 1DAY -> 1HOUR, or downsample.  However, hecmath likes to clip of days that don't have the complete 24 hour cycle.  
    So, we pad here, but there is a chance we ask for data not available. The read gives garbage data and doesn't complain.
    when a `timewindow` is supplied.
    """
    # TODO: figure out how to check for bounds for non-midnight start and end times.
    
    # Open input DSS
    dssFm = HecDss.open(inputDSSFile)  
    
    # Determine how to do the read
    if timewindow is not None:
        # Perform resampling on the specificed portion of the timeseries
        # Set the full record duration
        starttime_str = timewindow.getStartTimeString()  # initial start
        endtime_str = timewindow.getEndTimeString()  # initial end
        
        # Reformat the time
        starttime_str = starttime_str[:-4] + '0000'  # coerce to midnight start
        endtime_str = endtime_str[:-4] + '2400'  # coerce to midnight end
        
        # Adjust the padding
        if pad_start_days > 0:
            dt_start = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=pad_start_days)  # pad backwards
            starttime_str = dt_start.strftime('%d%b%Y %H%M')  # reformat to HEC string
        
        if lookback_1mon:
            dt_start = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=-31)  # original code (minus negative)
            starttime_str = dt_start.strftime('%d%b%Y %H%M')  # format back
        
        # Read the sieres
        print('Resampling', newPeriod, inputRec, starttime_str, endtime_str)  # diagnostic
        tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)  # windowed read
    
    else:
        # Resample on the whole series
        print('Resampling', newPeriod, inputRec)  # diagnostic
        tsm = dssFm.read(inputRec)  # full-record read (caution per original)

    # Perform the resampling action
    tsm_new = tsm.transformTimeSeries(newPeriod, "", "AVE")
    
    # Close the DSS file
    dssFm.close()  # close input DSS

    # Open, write to, and close the file
    dssFmout = HecDss.open(outputDSSFile)  # open output DSS
    dssFmout.write(tsm_new)  # write resampled series
    dssFmout.close()  # close output


def airtemp_lapse(dss_file, dss_rec, lapse_in_C, dss_outfile, f_part):
    """
    Apply a constant lapse to an air-temperature series and write to new record.

    Parameters
    ----------
    dss_file : str
        Input DSS file path.
    dss_rec : str
        Input air-temperature record pathname.
    lapse_in_C : float
        Lapse amount to add (°C). Converted to °F if units are Fahrenheit.
    dss_outfile : str
        Output DSS file path.
    f_part : str
        F-part for the destination record.

    Returns
    -------
    None
    """
    
    # Open the input DSS
    dss = HecDss.open(dss_file)  
    
    # Read the timeseries
    tsm = dss.read(dss_rec)
    
    # Create a copy of the series
    lapse = lapse_in_C  # local copy
    
    # Convert the units from F to C
    if 'f' in tsm.getUnits().lower():
        lapse = lapse * 9.0 / 5.0 + 32.0  # convert to °F if series units are Fahrenheit
    
    # Add the values to the series
    tsm = tsm.add(lapse) 
    
    # Get the values to write to the new series
    tsc = tsm.getData()  
    
    # Clsoe the input file
    dss.close()  

    # Write the output to a new file
    # Separate out the current path
    pathparts = dss_rec.split('/')  # split pathname
    pathparts[-2] = f_part  # change F-part
    
    # Rebuild fullName
    tsc.fullName = '/'.join(pathparts)  
    
    # Write
    dss_out = HecDss.open(dss_outfile)  # open output DSS
    dss_out.write(tsc)  # write adjusted series
    dss_out.close()  # close output DSS


def min_ts(dss_file, dss_rec, min_value, dss_outfile, f_part):
    """
    Clamp a time series to a minimum value and write to new record.

    Parameters
    ----------
    dss_file : str
        Input DSS file path.
    dss_rec : str
        Input record pathname.
    min_value : float
        Minimum value floor to apply.
    dss_outfile : str
        Output DSS file path.
    f_part : str
        F-part for the destination record.

    Returns
    -------
    None
    """
    
    # Open, read, and close the input file
    dss = HecDss.open(dss_file)  # open input DSS
    tsc = dss.get(dss_rec, True)  # full-record read
    dss.close()  # close input

    # Loop and compare each value in the timesereies
    for vi, v in enumerate(tsc.values):
        tsc.values[vi] = max(v, min_value)  # enforce minimum per element

    # Write the output to the output file
    # Separate out the current path
    pathparts = dss_rec.split('/')  # split pathname
    pathparts[-2] = f_part  # update F-part
    
    # Rebuild the path name
    tsc.fullName = '/'.join(pathparts)
    
    # Write
    dss_out = HecDss.open(dss_outfile)  # open output DSS
    dss_out.write(tsc)  # write clamped series
    dss_out.close()  # close output DSS


def add_flows(currentAlt, timewindow, inflow_records, dss_file, output_dss_record_name, output_dss_file):
    """
    Sum multiple inflow series (unit-converted if needed) within a window and write output.

    Parameters
    ----------
    currentAlt : object
        Alternative context for logging.
    timewindow : object
        Provides HEC-format start/end strings.
    inflow_records : list of str
        DSS paths for inflows; may include "file::path" for alternate files.
    dss_file : str
        Default DSS file to read when not using "file::path".
    output_dss_record_name : str
        Destination record pathname for summed flows.
    output_dss_file : str
        DSS file to write the output.

    Returns
    -------
    None
    """
    
    # Create the time for the series
    starttime_str = timewindow.getStartTimeString()  # window start
    endtime_str = timewindow.getEndTimeString()  # window end
    
    # Convert it to the hec format
    starttime_hectime = HecTime(starttime_str).value()  # numeric HEC start
    endtime_hectime = HecTime(endtime_str).value()  # numeric HEC end
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))  # log
    
    # Open the input file
    dssFm = HecDss.open(dss_file)  # open primary DSS

    # Define lists to hold the data
    inflows = []  # accumulator for sums
    times = []  # reference times

    # Read inflows
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records):  # iterate paths
        # Assign the path name
        pathname = inflow_record  # alias
        currentAlt.addComputeMessage('reading' + str(pathname))  # log
        print('\nreading' + str(pathname))  # diagnostic
        
        # Attempt to read the data from the file
        try:
            print(starttime_str, endtime_str)  # diagnostic times
            if '::' in inflow_record:
                dss_file_alt, inflow_rec_alt = inflow_record.split('::')  # alternate DSS file syntax
                dssFm_alt = HecDss.open(dss_file_alt)  # open alt DSS
                ts = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False)  # windowed read
                dssFm_alt.close()  # close alt DSS
                print(dss_file_alt)  # echo file
            
            else:
                print(dss_file)  # echo default file
                ts = dssFm.read(pathname, starttime_str, endtime_str, False)  # windowed read
            
            ts_data = ts.getData()  # unwrap to container
            values = ts_data.values  # series values
            hectimes = ts_data.times  # HEC times
            units = ts_data.units  # units string
            tstype = ts_data.type  # series type
            
            if hectimes[0] < starttime_hectime:  # if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])  # compute offset steps
                values = values[st_offset:]  # trim start values
                hectimes = hectimes[st_offset:]  # trim start times
            
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])  # compute tail steps
                values = values[:(len(hectimes) - st_offset)]  # trim end values
                hectimes = hectimes[:(len(hectimes) - st_offset)]  # trim end times

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))  # log failure
            sys.exit(-1)  # abort

        # Convert the units to imperial
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')  # log conversion
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)  # cms -> cfs
            values = convvals  # replace with converted

        if len(inflows) == 0:
            inflows = values  # seed sum
            times = hectimes  # store times (TODO: missing values handling)
        
        else:
            for vi, v in enumerate(values):
                inflows[vi] += v  # element-wise sum

    # Output record
    tsc = TimeSeriesContainer()  # build output container
    tsc.times = times  # assign times
    tsc.fullName = output_dss_record_name  # destination path
    tsc.values = inflows  # summed flows
    tsc.units = 'CFS'  # output in CFS
    tsc.type = tstype  # carry type
    tsc.numberValues = len(inflows)  # series length
    dssFm_out = HecDss.open(output_dss_file)  # open output DSS
    dssFm_out.write(tsc)  # write output

    dssFm.close()  # close input DSS
    dssFm_out.close()  # close output DSS


def add_or_subtract_flows(currentAlt, timewindow, inflow_records, dss_file, operation,
                          output_dss_record_name, output_dss_file, what="flow", prepend_n=0):
    """
    Combine flows (add or subtract) across records, with optional unit handling and preprend.

    Parameters
    ----------
    currentAlt : object
        Alternative context for logging.
    timewindow : object
        Provides HEC-format start/end strings.
    inflow_records : list of str
        DSS pathnames; may include "file::path" syntax for alternate files.
    dss_file : str
        Default DSS file path.
    operation : list of bool
        For each record: True to add, False to subtract.
    output_dss_record_name : str
        Destination DSS record pathname.
    output_dss_file : str
        DSS file to write the output.
    what : str, optional
        If "flow", convert cms->cfs and label as CFS; otherwise use `what` as units.
    prepend_n : int, optional
        Number of time steps to prepend the first value (lookback padding).

    Returns
    -------
    None
    """

    # Create the time for the series
    starttime_str = timewindow.getStartTimeString()  # window start
    endtime_str = timewindow.getEndTimeString()  # window end

    # Convert it to the hec format
    starttime_hectime = HecTime(starttime_str).value()  # numeric start
    endtime_hectime = HecTime(endtime_str).value()  # numeric end
    currentAlt.addComputeMessage('add_or_subtract_flows - Looking from {0} to {1}'.format(starttime_str, endtime_str))  # log
    
    # Open the input file
    dssFm = HecDss.open(dss_file)

    # Create lists to hold the data
    inflows = []  
    times = []  

    # Read inflows
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records):  # iterate records
        # Assign the pathanme
        pathname = inflow_record  # alias
        currentAlt.addComputeMessage('reading' + str(pathname))  # log
        print('\nreading' + str(pathname))  # diagnostic
        try:
            if '::' in inflow_record:
                dss_file_alt, inflow_rec_alt = inflow_record.split('::')  # alternate DSS file syntax
                dssFm_alt = HecDss.open(dss_file_alt)  # open alternate file
                ts = dssFm_alt.read(inflow_rec_alt, starttime_str, endtime_str, False)  # windowed read
                dssFm_alt.close()  # close alt file
                print(dss_file_alt)  # echo file
            
            else:
                ts = dssFm.read(pathname, starttime_str, endtime_str, False)  # primary read
                print(dss_file)  # echo file
            
            ts_data = ts.getData()  # unwrap to container
            values = ts_data.values  # numeric values
            hectimes = ts_data.times  # HEC times
            units = ts_data.units  # units
            tstype = ts_data.type  # type

            if hectimes[0] < starttime_hectime:  # trim start if before window
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])  # steps
                values = values[st_offset:]  # trim leading values
                hectimes = hectimes[st_offset:]  # trim leading times
            
            if hectimes[-1] > endtime_hectime:  # trim end if after window
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])  # steps
                values = values[:(len(hectimes) - st_offset)]  # trim trailing values
                hectimes = hectimes[:(len(hectimes) - st_offset)]  # trim trailing times

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))  # log error
            sys.exit(-1)  # abort

        if what == "flow":
            if units.lower() == 'cms':
                currentAlt.addComputeMessage('Converting cms to cfs')  # log conversion
                convvals = []

                for flow in values:
                    convvals.append(flow * 35.314666213)  # cms -> cfs

                values = convvals  # use converted values

        if len(inflows) == 0:
            inflows = values  # seed series
            times = hectimes  # reference times

        else:
            if operation[j]:
                for vi, v in enumerate(values):
                    inflows[vi] += v  # add series element-wise

            else:
                for vi, v in enumerate(values):
                    inflows[vi] -= v  # subtract element-wise

    if prepend_n > 0:
        # add first value onto front of record
        # Sometimes ResSim needs some lookback values, or whatever
        p_times = [times[i] for i in range(len(times))]  # copy to list
        time_delta = times[1] - times[0]  # one-step delta
        times = [p_times[0] - time_delta * i for i in range(prepend_n, 0, -1)] + p_times  # prepend times
        values = [inflows[i] for i in range(len(inflows))]  # copy values list
        inflows = [values[0]] * prepend_n + values  # prepend first value

    # Output record
    tsc = TimeSeriesContainer()  # build output container
    tsc.times = times  # assign times
    tsc.fullName = output_dss_record_name  # destination path
    tsc.values = inflows  # combined values
    
    if what == "flow":
        tsc.units = 'CFS'  # set units
    else:
        tsc.units = what  # custom units
    
    tsc.type = tstype  # carry type
    tsc.numberValues = len(inflows)  # series length
    dssFm_out = HecDss.open(output_dss_file)  # open output DSS
    dssFm_out.write(tsc)  # write series

    dssFm.close()  # close input DSS
    dssFm_out.close()  # close output DSS


def create_constant_dss_rec(currentAlt, timewindow, output_dss_file, constant=0.0, what='flow',
                            dss_type='PER-AVER', period='1HOUR', cpart='ZEROS', fpart='ZEROS'):
    """
    Create a constant-valued DSS series over a window and write it.

    Parameters
    ----------
    currentAlt : object
        Alternative for logging.
    timewindow : object
        Provides HEC-format start/end strings and HEC times.
    output_dss_file : str
        Destination DSS file path.
    constant : float, optional
        Constant value to fill the series (default 0.0).
    what : {'flow', 'temp-water', 'gate', 'evap', 'elev'}, optional
        Determines units and parameter part.
    dss_type : str, optional
        DSS type string (default 'PER-AVER').
    period : {'1HOUR', '1DAY'}, optional
        Time interval to use.
    cpart : str, optional
        Location (C-part) for output path.
    fpart : str, optional
        Version (F-part) for output path.

    Returns
    -------
    bool
        True on success; False for unsupported `what` or `period`.
    """

    # Determine what variable type is being processed
    if what.lower() == 'flow':
        units = 'cfs'  # cubic feet per second
        parameter = 'flow'  # parameter part
    
    elif what.lower() == 'temp-water':
        units = 'C'  # degrees Celsius
        parameter = 'temp-water'  # parameter part
    
    elif what.lower() == 'gate':
        units = 'n/a'  # unitless or not applicable
        parameter = 'gate'  # parameter part
    
    elif what.lower() == 'evap':
        units = 'ft'  # feet
        parameter = 'evap'  # parameter part
    
    elif what.lower() == 'elev':
        units = 'ft'  # feet
        parameter = 'elev'  # parameter part
    
    else:
        currentAlt.addComputeMessage('create_zero_dss_rec: what not known: %s' % what)  # log
        return False  # unsupported category

    # Define the time increment of the series
    if period.lower() == '1hour':
        pass
    
    elif period.lower() == '1day':
        pass
    
    else:
        currentAlt.addComputeMessage('create_zero_dss_rec: period not known: %s' % period)  # log
        return False  # unsupported period

    # Set the format string
    dt_format = '%d%b%Y %H%M'  # HEC format

    # Get the extents of the time series
    starttime_str = timewindow.getStartTimeString()  # start string
    endtime_str = timewindow.getEndTimeString()  # end string

    # pad 1 day on records, in case these are used for lookbacks, or balance flow calcs, etc.
    starttime_dt = hec_str_time_to_dt(starttime_str) - datetime.timedelta(days=1)  # pad -1 day
    endtime_dt = hec_str_time_to_dt(endtime_str) + datetime.timedelta(days=1)  # pad +1 day
    starttime_str_pad = starttime_dt.strftime(dt_format)  # format back
    endtime_str_pad = endtime_dt.strftime(dt_format)  # format back

    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))  # log

    ########################
    # Zero-Flow Time Series
    ########################

    tsmath_zero_flow_day = tsmath.generateRegularIntervalTimeSeries(
        starttime_str_pad,
        endtime_str_pad,
        period, "0M", constant)  # build constant series
    
    tsmath_zero_flow_day.setUnits(units)  # set units
    tsmath_zero_flow_day.setType(dss_type)  # set type
    tsmath_zero_flow_day.setTimeInterval(period)  # set interval
    tsmath_zero_flow_day.setLocation(cpart)  # set C-part
    tsmath_zero_flow_day.setParameterPart(parameter)  # set parameter part
    tsmath_zero_flow_day.setVersion(fpart)  # set F-part

    # Write the series to file
    dssFm = HecDss.open(output_dss_file)  # open output DSS
    dssFm.write(tsmath_zero_flow_day)  # write constant series
    dssFm.close()  # close output DSS

    # Return that the operation was successful
    return True


def calculate_relative_humidity(air_temp, dewpoint_temp):
    """
    Compute relative humidity (%) from air temperature (°C) and dew point (°C)
    using the August-Roche-Magnus approximation.

    Parameters
    ----------
    air_temp : float
        Air temperature in °C.
    dewpoint_temp : float
        Dew point temperature in °C.

    Returns
    -------
    float
        Relative humidity in percent, bounded to [0.01, 100].

    Notes
    -----
    Uses a simplified, commonly applied form suitable for typical ranges.
    """
  
 
    numerator = (112.0 - 0.1 * dewpoint_temp + air_temp)  # leading term
    denominator = (112.0 + 0.9 * air_temp)  # normalization
    
    exponent = ((17.62 * dewpoint_temp) / (243.12 + dewpoint_temp)) - ((17.62 * air_temp) / (243.12 + air_temp))  # exponent
    relative_humidity = 100.0 * (numerator / denominator) * math.exp(exponent)  # compute RH
    
    return max(0.01, min(100.0, relative_humidity))  # bound result


def calculate_dewpoint(air_temp, relative_humidity):
    """
    Calculate Dew Point Temperature given the air temperature and relative humidity, using the algebraic inversion of the 
    simplified August-Roche-Magnus approximation.

    Parameters
    ----------
    air_temp : float
        Air temperature in °C.
    relative_humidity : float
        Relative humidity (0-100%).

    Returns
    -------
    float
        Dew point temperature in °C.
    """

    gamma = math.log(relative_humidity / 100.0) + (17.62 * air_temp) / (243.12 + air_temp)  # intermediate term
    dewpoint = 243.12 * gamma / (17.62 - gamma)  # invert ARM
    
    return dewpoint  # dew point in °C


def relhum_from_at_dp(met_dss_file, at_path, dp_path):
    """
    Derive relative humidity (%) from air temperature and dew point records; write to derived DSS path.

    Parameters
    ----------
    met_dss_file : str
        DSS file containing inputs and destination.
    at_path : str
        Air temperature record pathname.
    dp_path : str
        Dew point temperature record pathname.

    Returns
    -------
    None
    """

    # Open the DSS file
    dss = HecDss.open(met_dss_file)
    
    # read AT container
    tsc = dss.read(at_path).getData()  
    
    # read DP values
    dp_data = data_from_dss(met_dss_file, dp_path, None, None)  
    
    # Loop over the series and compute the relative humidity for each step
    for i in range(tsc.numberValues):
        tsc.values[i] = calculate_relative_humidity(tsc.values[i], dp_data[i])
    
    # Format the path name for the new series
    parts = tsc.fullName.split('/')  # split pathname
    parts[2] = parts[2][:5]  # shorten C-part
    parts[3] = 'RELHUM-FROM-AT-DP'  # set parameter
    parts[6] = parts[6] + '-DERIVED'  # mark as derived
    new_pathname = '/'.join(parts)  # rebuild path
    
    # Create the output series
    tsc.fullName = new_pathname  # assign
    tsc.units = '%'  # set units
    dss.write(tsc)  # write derived series
    dss.close()  # close DSS


def dp_from_at_relhum(met_dss_file, at_path, rh_path):
    """
    Derive dew point temperature (°C) from air temperature and relative humidity records; write to derived DSS path.

    Parameters
    ----------
    met_dss_file : str
        DSS file containing inputs and destination.
    at_path : str
        Air temperature record pathname.
    rh_path : str
        Relative humidity record pathname.

    Returns
    -------
    None
    """
    
    # Open the DSS file
    dss = HecDss.open(met_dss_file)  # open DSS
    
    # Read the timeseries
    tsc = dss.read(at_path).getData()
    rh_data = data_from_dss(met_dss_file, rh_path, None, None)  # read RH values
    
    # Loop and calculate the dew point for each value
    for i in range(tsc.numberValues):
        tsc.values[i] = calculate_dewpoint(tsc.values[i], rh_data[i])  # compute dew point
    
    # Format the path name for the new series
    parts = tsc.fullName.split('/')  # split pathname
    parts[2] = parts[2][:5]  # shorten C-part
    parts[3] = 'temp-dewpoint'  # set parameter
    parts[6] = parts[6] + '-DERIVED'  # mark derived
    new_pathname = '/'.join(parts)  # rebuild
    
    # Create the output series
    tsc.fullName = new_pathname  # assign
    dss.write(tsc)  # write series
    dss.close()  # close DSS


def check_start_and_end(values, times, startime, endtime):
    """
    Trim arrays to fit within the [startime, endtime] boundaries.

    Parameters
    ----------
    values : list
        Numeric series values.
    times : list
        Numeric HEC time values (uniform step).
    startime : int
        HEC numeric start boundary.
    endtime : int
        HEC numeric end boundary.

    Returns
    -------
    (list, list)
        Tuple of (trimmed_values, trimmed_times).
    """
    
    # if startdate is before the timewindow
    if times[0] < startime:  
        print('start date ({0}) from DSS before timewindow ({1})..'.format(times[0], startime))
        st_offset = (startime - times[0]) / (times[1] - times[0])  # number of steps to trim
        values = values[st_offset:]  # drop leading values
        times = times[st_offset:]  # drop leading times
    
    # if end date is after the window..
    if times[-1] > endtime:  
        print('end date ({0}) from DSS after timewindow ({1})..'.format(times[-1], endtime))
        st_offset = (times[-1] - endtime) / (times[1] - times[0])  # number of steps to trim
        values = values[:(len(times) - st_offset)]  # drop trailing values
        times = times[:(len(times) - st_offset)]  # drop trailing times
    
    # Return trimmed outputs to the calling function
    return values, times


def replace_data(currentAlt, timewindow, pairs, dss_file, dss_outfile, months, standard_interval=None):
    """
    Replace base series values with an alternate series for selected months and write merged output.

    Parameters
    ----------
    currentAlt : object
        Alternative context for logging.
    timewindow : object
        Provides HEC-format start/end strings and HEC time values.
    pairs : list of (str, str)
        Tuples of (base_path, alt_path) for replacement.
    dss_file : str
        DSS file containing both input records.
    dss_outfile : str
        DSS file to write the merged output.
    months : list of int
        Month numbers (1-12) for which to replace base values.
    standard_interval : {'1hour', '1day', '1week'} or None, optional
        If provided, standardize both series before merging.

    Returns
    -------
    None
    """
    
    # Get the normal start and end times 
    starttime_str = timewindow.getStartTimeString()  # window start
    endtime_str = timewindow.getEndTimeString()  # window end
    
    # Convert the normal times to the HEC format
    starttime_hectime = HecTime(starttime_str).value()  # numeric start
    endtime_hectime = HecTime(endtime_str).value()  # numeric end
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))  # log window
    
    # Replace the data in each of the series
    for pair in pairs:
        # Get the values from the series to be filled
        # Open the DSS file
        dssFm = HecDss.open(dss_file)
        currentAlt.addComputeMessage('Replacing data for {0} with {1} during {2}'.format(pair[0], pair[1], months))  # log pair
        
        # Read the data from the file
        base = dssFm.read(pair[0], starttime_str, endtime_str, False)  # read base
        
        # Regularize the interval across the whole series
        if standard_interval is not None:
            base = standardize_interval(base, standard_interval)  # standardize base
        
        # Get the information from the series
        base_data = base.getData()  # container
        base_values = base_data.values  # values
        base_hectimes = base_data.times  # times
        base_units = base_data.units  # units
        base_interval = base_data.interval  # interval minutes
        base_type = base_data.type  # type
        base_values, base_hectimes = check_start_and_end(base_values, base_hectimes, starttime_hectime, endtime_hectime)  # trim base

        # Read the donor time series values
        alt = dssFm.read(pair[1], starttime_str, endtime_str, False)  # read alternate
        
        # Regularize the interval across the whole series
        if standard_interval is not None:
            alt = standardize_interval(alt, standard_interval)  # standardize alt
        
        # Get the information from the series
        alt_data = alt.getData()  # alt container
        alt_values = alt_data.values  # alt values
        alt_hectimes = alt_data.times  # alt times
        alt_units = alt_data.units  # alt units
        alt_interval = alt_data.interval  # alt interval
        alt_values, alt_hectimes = check_start_and_end(alt_values, alt_hectimes, starttime_hectime, endtime_hectime)  # trim alt

        # Close the DSS files
        dssFm.close()

        # Compare the units between series and enforce that they are the same
        if base_units != alt_units:
            currentAlt.addComputeMessage('Units do not match for {0} and {1}, skipping'.format(pair[0], pair[1]))  # log mismatch
            dssFm.close()
            sys.exit(1)  # abort
        
        # Compare the intervals between series and enforce that they are the same
        if base_interval != alt_interval:
            currentAlt.addComputeMessage('Intervals do not match for {0} and {1}, changing interval...'.format(pair[0], pair[1]))  # log mismatch
            dssFm.close()
            sys.exit(1)  # abort

        # Replace the values in the receipient series
        for i in range(len(base_values)):
            if base_data.getHecTime(i).month() in months:  # month selection
                base_values[i] = alt_values[i]  # replace value
                
        # Construct the new pathname
        new_pathname = base_data.fullName.split('/')  # base path tokens
        alt_pathname = alt_data.fullName.split('/')  # alt path tokens
        new_pathname[-2] = 'MergedFrom_{0}'.format(alt_pathname[1])  # set F-part to indicate source
        new_pathname = '/'.join(new_pathname)  # rebuild path
        
        # Write the series to the output file
        # Create the container
        tsc = TimeSeriesContainer()  # output container
        tsc.times = base_hectimes  # times from base
        tsc.fullName = new_pathname  # destination path
        tsc.values = base_values  # merged values
        tsc.units = base_units  # units preserved
        tsc.type = base_type  # type preserved
        tsc.numberValues = len(base_values)  # series length

        # Open the DSS file
        dssFmOut = HecDss.open(dss_outfile)
        
        # Write the container
        dssFmOut.write(tsc)

        # Close all DSS files
        dssFm.close()  
        dssFmOut.close()


def airtemp_lapse(dss_file, dss_rec, lapse_in_C, dss_outfile, f_part):
    """
    (Duplicate definition) Apply a constant lapse to an air-temperature series.

    Parameters
    ----------
    dss_file : str
        Input DSS file path.
    dss_rec : str
        Input air-temperature record pathname.
    lapse_in_C : float
        Lapse amount to add (°C). Converted to °F if units are Fahrenheit.
    dss_outfile : str
        Output DSS file path.
    f_part : str
        F-part for the destination record.

    Returns
    -------
    None
    """
    
    # Open the DSS file
    dss = HecDss.open(dss_file)
    
    # Read the timeseries
    tsm = dss.read(dss_rec)
    
    # Create a copy of the data
    lapse = lapse_in_C
    
    # Convert to °F if necessary
    if 'f' in tsm.getUnits().lower():
        lapse = lapse * 9.0 / 5.0 + 32.0  
    
    # add lapse
    tsm = tsm.add(lapse)  
    
    # unwrap container
    tsc = tsm.getData()  
    
    # Close the dss file
    dss.close()  

    # Write the new series to file
    pathparts = dss_rec.split('/')  # split pathname
    pathparts[-2] = f_part  # update F-part
    tsc.fullName = '/'.join(pathparts)  # rebuild pathname
    dss_out = HecDss.open(dss_outfile)  # open output DSS
    dss_out.write(tsc)  # write series
    dss_out.close()  # close output DSS


def preprend_first_value_on_ts(dss_file, dss_rec, prepend_n):
    """
    Prepend the first value of a series `prepend_n` times (for lookbacks) and write back.

    Parameters
    ----------
    dss_file : str
        DSS file path containing the record.
    dss_rec : str
        DSS pathname to modify.
    prepend_n : int
        Number of steps to prepend.

    Returns
    -------
    None
    """
    
    # Open the DSS file
    dss = HecDss.open(dss_file)
    
    # Read the record
    tsc = dss.get(dss_rec, True)

    # Calculate the interval of the timeseries
    time_delta = tsc.times[1] - tsc.times[0]  
    
    # Convert the time vector to a Python list
    times = [tsc.times[i] for i in range(len(tsc.times))]
    
    # Shift everything by a single timestep
    tsc.times = [times[0] - time_delta * i for i in range(prepend_n, 0, -1)] + times
    
    # Reset the start time
    tsc.startTime = tsc.times[0]
    
    # Copy the values into the new list
    values = [tsc.values[i] for i in range(len(tsc.values))]  # copy to Python list
    
    # Reset the values in the container
    tsc.values = [values[0]] * prepend_n + values  # prepend values
    tsc.numberValues = len(tsc.values)  # update length

    # Put values into the DSS file and close
    dss.put(tsc)  # write back to same record
    dss.close()  # close DSS


def postprend_last_value_on_ts(dss_file, dss_rec, postpend_n):
    """
    Append the last value of a regular series `postpend_n` times and write back.

    Parameters
    ----------
    dss_file : str
        DSS file path containing the record.
    dss_rec : str
        DSS pathname to modify.
    postpend_n : int
        Number of steps to append.

    Returns
    -------
    None

    Notes
    -----
    Assumes a regular (uniform-interval) record to compute time deltas.
    """
    
    # Open the DSS file
    dss = HecDss.open(dss_file)
    
    # Read the record from the file
    tsc = dss.get(dss_rec, True)

    # Calculate the interval of the timeseries
    time_delta = tsc.times[1] - tsc.times[0]
    
    # Convert the time vector to a Python list
    times = [tsc.times[i] for i in range(len(tsc.times))]
    
    # Shift everything by a single timestep
    tsc.times = times + [times[-1] + time_delta * i for i in range(1, postpend_n + 1)]
    
    # Copy the values to a new list
    values = [tsc.values[i] for i in range(len(tsc.values))]
    
    # Reset the values in the container
    tsc.values = values + [values[-1]] * postpend_n
    tsc.numberValues = len(tsc.values)

    # Put the values into the DSS file
    dss.put(tsc) 

"""
Simple DSS + Interpolation Utilities

This module contains helper routines used within HEC-WAT/ResSim/W2 scripting
for the Sacramento WTMP workflow. It includes:

- DSS record utilities (copying, reading, resampling, unit fixes, flow math).
- Time conversion helpers between HEC `HecTime`/Java `Date` and Python `datetime`.
- Lightweight linear interpolation for Python `datetime` series.
- Convenience functions for path manipulation and alternative-based resolution.
"""


# --- Utility Function for Python 2 timedelta.total_seconds() ---
# In Python 2, timedelta objects don't always have total_seconds(), so
# we implement it manually to ensure cross-version compatibility for float division.
def _timedelta_to_seconds(td):
    """
    Convert a `datetime.timedelta` to a floating-point number of seconds.

    This ensures proper float division for time ratios in Python 2,
    where `timedelta.total_seconds()` is not always available.

    Parameters
    ----------
    td : datetime.timedelta
        The time delta to convert.

    Returns
    -------
    float
        The total duration represented by `td`, in seconds.

    Notes
    -----
    Computed as ``days * 86400 + seconds + microseconds / 1e6``, i.e.
    86,400 seconds per day and 1,000,000 microseconds per second.
    """
    # 86400 seconds in a day (24 * 60 * 60)
    # 1,000,000 microseconds in a second
    return td.days * 86400.0 + td.seconds + td.microseconds / 1000000.0

# --- Main Interpolation Method ---
def linear_interp_datetime(time_data, value_data, query_times):
    """
    Perform linear interpolation on time-series data keyed by
    `datetime` objects.

    Parameters
    ----------
    time_data : list of datetime.datetime
        Ordered (ascending) timestamps, used as the interpolation
        x-axis.
    value_data : list of float
        Values corresponding to `time_data` (y-axis). Must be the
        same length as `time_data`.
    query_times : list of datetime.datetime
        Timestamps for which to interpolate values.

    Returns
    -------
    list of float
        Interpolated values corresponding to `query_times`. Returns
        an empty list if validation fails (see Notes).

    Notes
    -----
    - Requires at least 2 points in `time_data`/`value_data`; if this
      is not satisfied, an error is printed (Python 2 `print`
      statement) and an empty list is returned rather than raising.
    - Query times before the first known timestamp are clamped to the
      first value; query times after the last known timestamp are
      clamped to the last value (i.e. no true extrapolation is
      performed, despite the naming in the code comments).
    - Internally converts all timestamps to elapsed seconds since
      `time_data[0]` via `_timedelta_to_seconds` before interpolating.
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
    Copy a DSS time series record to a new record path.

    Copies either by substituting a new F-part or by using an
    entirely new record path, optionally normalizing degree units and
    converting Fahrenheit to Celsius, and optionally writing to a
    different DSS file.

    Parameters
    ----------
    dss_rec : str
        DSS path of the source record to copy.
    new_fpart : str, optional
        If provided (and `new_dss_rec` is not), the F-part to
        substitute into `dss_rec`'s path to form the output path.
    new_dss_rec : str, optional
        If provided, used directly as the output DSS path instead of
        substituting an F-part.
    dss_file_path : str, optional
        Path to the DSS file to open for reading (and writing, unless
        `dss_file_alt_outpath` is given). Required if
        `dss_file_handle` is not provided.
    dss_file_handle : object, optional
        An already-open DSS file handle to read from, instead of
        opening `dss_file_path`.
    dss_file_alt_outpath : str, optional
        If provided, the copied record is written to this DSS file
        instead of the source file.
    checkMakeCelsius : bool, optional
        If True and the source record's units are Fahrenheit,
        converts all values to Celsius before writing. Defaults to
        False.

    Returns
    -------
    None
        Writes the copied (and possibly converted) record to the
        appropriate output DSS file.

    Raises
    ------
    ValueError
        If neither `dss_file_path` nor `dss_file_handle` is provided,
        or if neither `new_fpart` nor `new_dss_rec` is provided.
    """
    # error check inputs - there are flexible way to copy record
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
    # build the output path, either using the explicit new record name or by substituting a new F-part into the source path
    if new_dss_rec is not None:
        # Explicit destination takes priority 
        dss_rec_out = new_dss_rec  
    
    else:
        # No pathname provided. Use default.
        rec_parts = tsc.fullName.split('/')  # tokenize pathname
        rec_parts[6] = new_fpart  # update F-part
        dss_rec_out = '/'.join(rec_parts)  # rebuild pathname

    # fix some terrible units along the way (normalize non-standard degree unit labels to the standard C/F abbreviations)
    if tsc.units.lower() == 'degc':
        tsc.units = 'C'  # normalize to 'C'
    
    elif tsc.units.lower() == 'degf':
        tsc.units = 'F'  # normalize to 'F'

    # Convert from F to C if needed
    if tsc.units.lower() == 'f' and checkMakeCelsius:
        T_values = tsc.values  # get values
        
        for i, TT in enumerate(T_values):
            T_values[i] = (TT - 32.0) * 5.0 / 9.0  # convert F -> C in-place
            
        tsc.units = 'C'
        tsc.values = T_values
        
    # write
    tsc.fullName = dss_rec_out
    # write to the alternate output file if one was given, otherwise write back to the source file
    if dss_file_alt_outpath is not None:
        dss_fm_alt = HecDss.open(dss_file_alt_outpath)  # open alternate output DSS
        dss_fm_alt.put(tsc)  # write to alt file
        dss_fm_alt.close()  # close alt
    
    else:
        dss_fm.put(tsc)  # write to same file

    # only close the file handle if this function opened it itself
    if dss_file_handle is None:
        dss_fm.close()  


def jday_from_tsc(tsc):
    """
    Compute the fractional (decimal) day-of-year for every timestamp
    in a time series container.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        A DSS time series data container with time and value data.

    Returns
    -------
    list of float
        The decimal day-of-year for each timestamp in `tsc`, in the
        same order.
    """
    dtt = hectime_to_datetime(tsc)
    return [decimal_doy(dt) for dt in dtt]        


def decimal_doy(dt):
    """
    Compute the fractional (decimal) day-of-year for a single datetime.

    The integer part is the calendar day-of-year and the fractional
    part represents the time of day.

    Parameters
    ----------
    dt : datetime.datetime
        The date/time to convert.

    Returns
    -------
    float
        The decimal day-of-year (e.g. Jan 1 at noon would be
        approximately 1.5).
    """
    doy = dt.timetuple().tm_yday
    fractional_day = (dt.hour / 24.0) + (dt.minute / 1440.0) + (dt.second / 86400.0) + (dt.microsecond / 86400000000.0)
    
    # Sum and return integer and fractional parts
    return doy + fractional_day  


def organizeLocations(curAlt, location_objs, loc_names, return_dss_paths=False):
    """
    Look up a list of named locations within a set of location
    objects, in the order given by `loc_names`.

    Optionally resolves each to its DSS path instead of returning the
    location object itself.

    Parameters
    ----------
    curAlt : object
        The alternative object being computed. Passed through to
        `findLocationOrder` and used for resolving linked time series
        paths.
    location_objs : list
        The full list of location objects to search within (e.g. from
        `currentAlternative.getInputDataLocations()`).
    loc_names : list of str
        Names to look up, in the desired output order.
    return_dss_paths : bool, optional
        If True, resolves and returns each location's DSS path string
        instead of the location object itself. Defaults to False.

    Returns
    -------
    list
        If `return_dss_paths` is False, a list of location objects
        corresponding to `loc_names`, in that order. If True, a list
        of resolved DSS path strings instead.

    Raises
    ------
    SystemExit
        Via `findLocationOrder`, if any name in `loc_names` cannot be
        found (calls `sys.exit(1)`).
    """
    locations_list = []
    print('num_locs:',len(location_objs))
    # for each requested location name, find its matching location object and resolve it
    for name in loc_names:
        # Get the index of the location
        i_loc = findLocationOrder(curAlt, location_objs, name)
        
        # Extract the DSS paths if requested
        if return_dss_paths:
            lo1 = location_objs[i_loc]
            print(lo1)
            # if this location is linked to an upstream model, load its time series and fix
            # its F-part; otherwise, it already has a fixed DSS path to use directly
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
    Apply `organizeLocations` to multiple groups of location names at
    once.

    Parameters
    ----------
    curAlt : object
        The alternative object being computed, passed through to
        `organizeLocations`.
    location_objs : list
        The full list of location objects to search within.
    loc_names_paired : list of list of str
        A list of name groups; each inner list is passed to
        `organizeLocations` independently.
    return_dss_paths : bool, optional
        Passed through to `organizeLocations` for each group.
        Defaults to False.

    Returns
    -------
    list of list
        One resolved list (of location objects or DSS paths,
        depending on `return_dss_paths`) per input group in
        `loc_names_paired`, in the same order.
    """
    return [organizeLocations(curAlt, location_objs, pn, return_dss_paths) for pn in loc_names_paired]


def findLocationOrder(curAlt,location_objs,name):
    """
    Find the index of a location object within a list, matching by
    name.

    Parameters
    ----------
    curAlt : object
        The alternative object being computed. Used to log an error
        message if the name cannot be found.
    location_objs : list
        The list of location objects to search within. Each must
        support `getName()`.
    name : str
        The location name to search for.

    Returns
    -------
    int
        The index of the first location object whose `getName()`
        matches `name`.

    Raises
    ------
    SystemExit
        Calls `sys.exit(1)` if no match is found, after logging an
        error message via `curAlt.addComputeMessage()`.
    """
    # search through every location for one whose name matches
    for i,loc in enumerate(location_objs):
        print(i,'Checking loc: ',loc.getName())
        print('loc:',loc)
        if name == loc.getName():
            # return match position
            return i  
    
    # if we make it here, our input/output location name was not found
    curAlt.addComputeMessage("Scripting - Location name not found: " + name)  # log
    
    # abort to indicate error
    sys.exit(1)  

def first_value(dss_file,dss_rec,start_str=None,end_str=None):
    """
    Read a DSS record and return only its first value.

    Parameters
    ----------
    dss_file : str
        Path to the DSS file to read from.
    dss_rec : str
        DSS path of the record to read.
    start_str : str, optional
        Start time string for a bounded read. If None (along with
        `end_str`), the entire record is read instead.
    end_str : str, optional
        End time string for a bounded read. If None (along with
        `start_str`), the entire record is read instead.

    Returns
    -------
    float
        The first value in the record.
    """
    dssFm = HecDss.open(dss_file)        
    # if no time bounds were given, read the entire record; otherwise read only the given window
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
    Resample a time series math object to a standard interval,
    averaging values within each new interval.

    Standard intervals are hourly, daily, or weekly.

    Parameters
    ----------
    tsm : object
        A HEC time series math object (`tsmath`) to standardize.
    interval : str
        The target interval. Must be one of `'1hour'`, `'1day'`, or
        `'1week'` (case-insensitive).
    makePerAver : bool, optional
        If True, sets the series type to `'PER-AVER'` before
        transforming, ensuring the resample treats values as period
        averages. Defaults to True.

    Returns
    -------
    object
        A `tsmath` object at the requested interval. If `tsm` is
        already at that interval, `tsm` is returned unchanged (no
        averaging is performed).

    Raises
    ------
    SystemExit
        Calls `sys.exit(-1)` if `interval` is not one of the
        recognized values.
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
        # Invalid condition. Error and exit
        sys.exit(-1)

    # only resample if the current interval differs from the requested one
    if tsc.interval != intint:
        if makePerAver:
            tsm.setType('PER-AVER')  # set series type on math object
        return tsm.transformTimeSeries(interval, "", "AVE")  # perform transform
    
    else:
        return tsm  # already standardized


def get_sanitized_record_list(dss_file_path):
    """
    Return a deduplicated list of DSS record paths with the date
    field blanked out.

    The DSS library seems to return lists of paths with dates in them
    (`getPathnameList()`), and some of those dates don't even exist
    in the file or cannot be read and throw an error. As of Jan 2024,
    this is a recurring problem, and the manual solution is throwing
    away the DSS file, but in many cases that is problematic. So,
    here we filter dates and check for duplicates.

    Parameters
    ----------
    dss_file_path : str
        Path to the DSS file to list records from.

    Returns
    -------
    list of str
        A deduplicated list of DSS record paths with the date
        (D-part) field blanked out, so that records differing only by
        date are treated as the same record.
    """
    dss = HecDss.open(dss_file_path)
    recs = dss.getPathnameList()
    dss.close()
    sanitized_recs = []
    
    # for each raw record path, blank out its date field and skip if already seen
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
    Extract a list of `HecTime` objects from a time series math
    object's underlying data container.

    Parameters
    ----------
    tsm : object
        A HEC time series math object (`tsmath`).

    Returns
    -------
    list of HecTime
        One `HecTime` object per value in `tsm`, constructed from its
        raw integer time values.
    """
    times = tsm.getContainer().times
    htimes = []
    for i in range(tsm.getContainer().numberValues):
        htimes.append(HecTime())  # construct new HecTime per value
        htimes[-1].set(times[i])  # assign numeric time
    
    # Return list of HecTime
    return htimes  


def hectimes_from_tsc(tsc):
    """
    Extract a list of `HecTime` objects from a time series data
    container.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        A DSS time series data container.

    Returns
    -------
    list of HecTime
        One `HecTime` object per value in `tsc`, constructed from its
        raw integer time values.
    """
    htimes = []
    # for each raw time value, build a corresponding HecTime object
    for i in range(tsc.numberValues):
        htimes.append(HecTime())  # new HecTime
        htimes[-1].set(tsc.times[i])  # assign from container times
    
    # Return the output list of values
    return htimes


def shift_pit_river_time(input_dss_file,dss_rec,output_dss_file,out_rec,start_date=None,end_date=None):
    """
    Read a time series, shift it 12 hours earlier in time, and write
    the result to a new record.

    Parameters
    ----------
    input_dss_file : str
        Path to the DSS file containing the source record.
    dss_rec : str
        DSS path of the source record to shift.
    output_dss_file : str
        Path to the DSS file to write the shifted record to.
    out_rec : str
        DSS path for the shifted output record. Must differ from
        `dss_rec`.
    start_date : str, optional
        Start time string for a bounded read. If None, the entire
        record is read.
    end_date : str, optional
        End time string for a bounded read. If None, the entire
        record is read.

    Returns
    -------
    None
        Writes the shifted record to `output_dss_file`.

    Notes
    -----
    It is important that `dss_rec` and `out_rec` differ. If they are
    the same, running this function repeatedly will keep shifting the
    same record further back in time on every run.
    """
    # important to not have dss_rec==out_rec, otherwise you will continually shift the record in time
    # whenever the calling script runs...    

    tsm_pit = dss_read_ts_safe(input_dss_file,dss_rec,start_date=start_date,end_date=end_date,returnTSM=True)

    # Shift the values by -12 hours
    tsm_pit = tsm_pit.shiftInTime("-12Hour")
    
    # Request the data and path
    tsc_pit = tsm_pit.getData()  # unwrap container
    tsc_pit.fullName = out_rec  # set destination pathname

    # Output the shifted series
    dss_out = HecDss.open(output_dss_file)  # open destination
    dss_out.put(tsc_pit)  # write shifted series
    dss_out.close()  # close output DSS


def shift_ts_time(input_dss_file,dss_rec,output_dss_file,out_rec,shift_str,start_date=None,end_date=None):
    """
    Read a time series, shift it in time by an arbitrary HEC-style
    offset, and write the result to a new record.

    Parameters
    ----------
    input_dss_file : str
        Path to the DSS file containing the source record.
    dss_rec : str
        DSS path of the source record to shift.
    output_dss_file : str
        Path to the DSS file to write the shifted record to.
    out_rec : str
        DSS path for the shifted output record. Must differ from
        `dss_rec`.
    shift_str : str
        A HEC-recognized time shift string, with an optional leading
        minus sign (e.g. `'-12Hour'`, `'5Day'`).
    start_date : str, optional
        Start time string for a bounded read. If None, the entire
        record is read.
    end_date : str, optional
        End time string for a bounded read. If None, the entire
        record is read.

    Returns
    -------
    None
        Writes the shifted record to `output_dss_file`.

    Notes
    -----
    It is important that `dss_rec` and `out_rec` differ. If they are
    the same, running this function repeatedly will keep shifting the
    same record further in time on every run.
    """
    
    tsm = dss_read_ts_safe(input_dss_file,dss_rec,start_date=start_date,end_date=end_date,returnTSM=True)

    # Shift the values
    tsm = tsm.shiftInTime(shift_str)  # perform time shift
    
    # Request the data and the path
    tsc = tsm.getData()  # unwrap to container
    tsc.fullName = out_rec  # set destination pathname

    # Output the shifted values
    dss_out = HecDss.open(output_dss_file)  # open destination DSS
    dss_out.put(tsc)  # write shifted series
    dss_out.close()  # close output


def dss_read_ts_safe(dssFilePath,dssRec,start_date=None,end_date=None,returnTSM=False,
                     returnPydatetimes=False,debug=False):
    """
    Read a DSS record in a date-flexible way, ensuring the whole time
    range is returned if dates are None.

    Parameters
    ----------
    dssFilePath : str
        Path to the DSS file to read from.
    dssRec : str
        DSS path of the record to read.
    start_date : str, optional
        Start time string for a bounded read. If both `start_date`
        and `end_date` are None, the entire record is read via
        `get()` instead.
    end_date : str, optional
        End time string for a bounded read. If both `start_date` and
        `end_date` are None, the entire record is read via `get()`
        instead.
    returnTSM : bool, optional
        If True, returns the `tsmath` wrapper object instead of the
        raw data container. Defaults to False.
    returnPydatetimes : bool, optional
        Currently unused within this function body.
    debug : bool, optional
        If True, prints extra diagnostic information about the read.
        Defaults to False.

    Returns
    -------
    object or None
        The time series data, either as a raw data container or as a
        `tsmath` object (depending on `returnTSM`), or `None` if
        neither the "both dates None" nor "both dates given" branch
        is satisfied (e.g. if only one of `start_date`/`end_date` is
        provided).
    """
    # this method is going to throw an error if the file doesn't exist
    dss = HecDss.open(dssFilePath,True)

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
            return tsc
    # if both bounds were given, read only within that time window
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


def data_from_dss(dss_file,dss_rec,starttime_str, endtime_str):
    """
    Read a DSS record's values, either over its entire span or within
    a specific time window.

    Parameters
    ----------
    dss_file : str
        Path to the DSS file to read from.
    dss_rec : str
        DSS path of the record to read.
    starttime_str : str
        Start time string for a bounded read. If None (along with
        `endtime_str`), the entire record is read instead.
    endtime_str : str
        End time string for a bounded read. If None (along with
        `starttime_str`), the entire record is read instead.

    Returns
    -------
    list of float
        The values from the requested record/time window.
    """
    dssFm = HecDss.open(dss_file)
    # if no time bounds were given, read the entire record; otherwise read only the given window
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
    **Deprecated**: Convert HEC times in a `TimeSeriesContainer` to
    Python `datetime`.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        Container with HEC times accessible via `getHecTime(i)`.

    Returns
    -------
    list of datetime.datetime
        Python datetimes corresponding to each container time,
        converted via Java Date/timestamp (in GMT, with potential
        rounding issues - see Notes below).

    Notes
    -----
    Deprecated! This routine, and the underlying Java-time processing
    it depends on, is capable of rounding dates to the wrong hour. It
    also adds a timezone offset and returns GMT, which is not usually
    wanted.

    TODO (from original source): replace all instances of this with
    `datetimes_from_tsc()`, and rework/remove any extra processing
    previously needed for timezones, bad first/last values, etc.
    """
    
    # Get a list of hectimes
    hectimes = hectimes_from_tsc(tsc)
    
    # Define a list to hold values
    dtt = []
    # for each timestamp, convert via Java Date/timestamp into a Python datetime
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
    Convert a list of `HecTime` objects to Python `datetime` objects.

    Conversion goes via each `HecTime`'s long-format string
    representation.

    Parameters
    ----------
    hectimes : list of hec.heclib.util.HecTime
        HEC times to convert.

    Returns
    -------
    list of datetime.datetime
        One Python datetime per input `HecTime` object.

    Notes
    -----
    This function is defined twice in this module (identically).
    Because Python/Jython executes definitions top-to-bottom, the
    later definition (further down in this file) silently replaces
    this one at import time; any call to `hecTime_to_datetime` uses
    the later version. Both are documented here for completeness.
    """
    return [hec_str_time_to_dt(h.toString(), longStr=True) for h in hectimes]


def datetimes_from_tsc(tsc):
    """
    Convert a `TimeSeriesContainer`'s timestamps to Python
    `datetime` objects via `HecTime` conversion.

    This is the preferred, non-deprecated alternative to
    `hectime_to_datetime`.

    Parameters
    ----------
    tsc : hec.io.TimeSeriesContainer
        Source container to convert.

    Returns
    -------
    list of datetime.datetime
        Converted Python datetimes.

    Notes
    -----
    This function is defined twice in this module (identically). The
    later definition (further down in this file) silently replaces
    this one at import time.
    """
    return hecTime_to_datetime(hectimes_from_tsc(tsc)) 


def hec_str_time_to_dt(hec_str_time, longStr=False):
    """
    Convert a HEC time string to a Python `datetime`, handling the
    `"24:00"`/`"2400"` end-of-day rollover.

    Parameters
    ----------
    hec_str_time : str
        HEC-formatted time string. If `longStr` is True, uses verbose
        format like `"DD Month YYYY, HH:MM"`; otherwise
        `"DDMonYYYY HHMM"`.
    longStr : bool, optional
        Flag to select the verbose format with `24:00` rollover.
        Defaults to False (compact format with `2400` rollover).

    Returns
    -------
    datetime.datetime
        The parsed date/time. HEC's "24:00" end-of-day convention is
        converted to 00:00 of the following day, since Python's
        `datetime` cannot represent hour 24 directly.

    Notes
    -----
    This function is defined twice in this module (identically). The
    later definition (further down in this file) silently replaces
    this one at import time.
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


def hecTime_to_datetime(hectimes):
    """
    Convert a list of `HecTime` objects to Python `datetime` objects,
    via their long-format string representation.

    This is the active (final) definition of `hecTime_to_datetime`
    used at runtime, superseding the earlier identical definition
    above.

    Parameters
    ----------
    hectimes : list of HecTime
        The `HecTime` objects to convert.

    Returns
    -------
    list of datetime.datetime
        One Python datetime per input `HecTime` object.
    """
    return [hec_str_time_to_dt(h.toString(),longStr=True) for h in hectimes]

def fixInputLocationFpart(currentAlternative, tspath):
    """
    Rewrite a DSS path's F-part to match the current alternative's
    configured input F-part.

    Preserves the original F-part's trailing scenario/model
    identifier while replacing the version prefix.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `getInputFPart()` to provide the expected F-part prefix.
    tspath : str
        The DSS path whose F-part should be corrected.

    Returns
    -------
    str
        The DSS path with its F-part replaced by the alternative's
        configured F-part prefix, combined with the original path's
        trailing F-part segment.
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
    Append a suffix onto a DSS path's A-part (river/basin name).

    Separated by an underscore, or sets the A-part directly if it
    was previously empty.

    Parameters
    ----------
    current_path : str
        The DSS path whose A-part should be modified. (See Notes
        below regarding a bug in the current implementation.)
    ApartAppend : str
        The text to append to (or set as) the A-part.

    Returns
    -------
    str
        The DSS path with its A-part updated.

    Raises
    ------
    NameError
        This function references the name `tspath` before it is ever
        assigned - the parameter is named `current_path`, but the
        first line splits `tspath`, which does not exist yet. As
        written, calling this function will always raise a
        `NameError`. This has been left unchanged and is flagged here
        for visibility.
    """
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
    Resolve DSS path and file location for a model location,
    considering linkage.

    Parameters
    ----------
    location : object
        Location object supporting linkage queries.
    currentAlternative : object
        Alternative context for `loadTimeSeries()` and F-part
        adjustments.
    computeOptions : object
        Provides compute-time DSS filename.

    Returns
    -------
    tuple of (str, str)
        `(tspath, dsspath)` for inputs/outputs.
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
    Strip the first 4 characters from the F-part for all records,
    then rename.

    Only applies when the F-part contains a `'-'` character.

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
    Creates a `.bak` copy, builds new pathnames in a `Vector`, and
    performs a bulk rename via `HecDss.renameRecords()`. If the
    F-part does not contain a `'-'`, returns early (mirroring original
    code behavior) - note this means the function will exit on the
    *first* record examined if that record's F-part lacks a `'-'`,
    rather than skipping just that record.
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
        currentAlt.addComputeMessage('Fixing path: ' + r + ' --> ' + new_rec_names[-1])  
    
    # Rename all records concurrently
    dss.renameRecords(rec_names, new_rec_names)  
    
    # Close the DSS file
    dss.close()  


def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    """
    Sum multiple DSS time series records together, value by value,
    and write the combined series to a new output record.

    Parameters
    ----------
    currentAlt : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging.
    dssFile : str
        Path to the DSS file to read from and write to.
    timewindow : object
        The run time window object, used to get the start and end
        time strings for reading input data.
    input_data : list of str
        DSS paths of the records to sum together. All records are
        assumed to share the same timestamps, units, and type.
    output_path : str
        DSS path to write the combined series to.

    Returns
    -------
    int
        Always returns 0 on completion.
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    
    # Create the data holder
    output_data = []
    # for each input record, read it and add its values into the running total
    for dsspath in input_data:
        print('reading', str(dsspath))  # diagnostic
        # Read the time series
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)  
        
        # Get the container and the data
        ts = ts.getData()  
        values = ts.values
        times = ts.times
        units = ts.units
        dsstype = ts.type
        # the first record seeds the running total; subsequent records are added on
        if len(output_data) == 0:
            output_data = values  # seed with first series
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
    """
    Resample a DSS time series to a new (regular) period, upsampling
    or downsampling as needed.

    Can upsample an even-period DSS time series, e.g. go from `1DAY`
    to `1HOUR`, or downsample. However, `hecmath` tends to clip off
    days that don't have a complete 24-hour cycle, so this function
    pads the read window to full-day boundaries first.

    Parameters
    ----------
    inputDSSFile : str
        Path to the DSS file containing the source record.
    inputRec : str
        DSS path of the record to resample.
    timewindow : object or None
        The run time window object, used to determine the read window
        (padded to full-day boundaries). If None, the entire record
        is read and resampled without windowing.
    outputDSSFile : str
        Path to the DSS file to write the resampled record to.
    newPeriod : str
        The target DSS interval string (e.g. `'1HOUR'`, `'1DAY'`).
    lookback_1mon : bool, optional
        If True, pads the start of the read window backward using a
        `timedelta(days=-31)` (note: this subtracts a negative,
        effectively adding 31 days - see Notes). Defaults to False.
    pad_start_days : int, optional
        Number of extra days to pad onto the start of the read
        window, applied before `lookback_1mon`. Defaults to 0.

    Returns
    -------
    None
        Writes the resampled series to `outputDSSFile`.

    Notes
    -----
    - Since incomplete days are clipped during transform, the read
      window's start/end times are forced to `0000`/`2400`
      respectively when `timewindow` is given.
    - There is a chance this requests data outside what is actually
      available in the source record; in that case the read may
      return garbage/undefined values without raising an error.
    - TODO (from original source): figure out how to check bounds for
      non-midnight start and end times.
    - The `lookback_1mon` branch computes
      ``hec_str_time_to_dt(starttime_str) - timedelta(days=-31)``,
      which is equivalent to *adding* 31 days rather than looking
      back a month; this appears to be a pre-existing quirk in the
      original code and has been left unchanged.
    """
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


def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    """
    Apply a fixed temperature lapse (offset) to an air temperature
    record.

    Approximates temperature change with elevation, writing the
    adjusted record under a new F-part.

    Parameters
    ----------
    dss_file : str
        Path to the DSS file containing the source air temperature
        record.
    dss_rec : str
        DSS path of the source air temperature record.
    lapse_in_C : float
        The temperature lapse/offset to apply, in degrees Celsius. If
        the source record's units are Fahrenheit, this is converted
        to an equivalent Fahrenheit offset before being applied.
    dss_outfile : str
        Path to the DSS file to write the adjusted record to.
    f_part : str
        The F-part to assign to the output record.

    Returns
    -------
    None
        Writes the lapse-adjusted record to `dss_outfile`.

    Notes
    -----
    This function is defined twice in this module (once here, and
    again later, identically). Because Python/Jython executes
    definitions top-to-bottom, the later definition silently
    overrides this one at import time; any call to `airtemp_lapse`
    uses that version. Consider removing one of the two copies.
    """
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    lapse = lapse_in_C
    # if the source series is in Fahrenheit, convert the lapse amount to an equivalent Fahrenheit offset
    if 'f' in tsm.getUnits().lower():
        lapse = lapse * 9.0 / 5.0 + 32.0  # convert to F if series units are Fahrenheit
    
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

def min_ts(dss_file,dss_rec,min_value,dss_outfile,f_part):
    """
    Clamp every value in a time series to a minimum value.

    Writes the adjusted record under a new F-part.

    Parameters
    ----------
    dss_file : str
        Path to the DSS file containing the source record.
    dss_rec : str
        DSS path of the source record to clamp.
    min_value : float
        The minimum allowed value; any value below this is raised to
        this value.
    dss_outfile : str
        Path to the DSS file to write the adjusted record to.
    f_part : str
        The F-part to assign to the output record.

    Returns
    -------
    None
        Writes the clamped record to `dss_outfile`.
    """
    dss = HecDss.open(dss_file)
    tsc = dss.get(dss_rec,True)
    dss.close()

    # for each value, raise it to the minimum if it falls below that threshold
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
    Sum multiple flow time series together and write the combined
    flow record to a new output file.

    Performs unit conversion from cms to cfs where needed, and trims
    each input record to the run time window.

    Parameters
    ----------
    currentAlt : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging.
    timewindow : object
        The run time window object, used to determine the valid time
        range for trimming input records.
    inflow_records : list of str
        DSS paths of the records to sum. A path may optionally use
        the form `'other_dss_file::record_path'` to read from a
        different DSS file than `dss_file`.
    dss_file : str
        Path to the primary DSS file to read records from (used
        unless a record path specifies an alternate file via
        `'::'`).
    output_dss_record_name : str
        DSS path to write the summed result to.
    output_dss_file : str
        Path to the DSS file to write the summed result to.

    Returns
    -------
    None
        Writes the combined flow series (in CFS) to
        `output_dss_file`.

    Raises
    ------
    SystemExit
        Calls `sys.exit(-1)` if a `HecMathException` occurs while
        reading any of the input records.

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

    # Define lists to hold the data
    inflows = []  # accumulator for sums
    times = []  # reference times

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
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            # if the record's data extends past the run window, trim the trailing values
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])  # compute tail steps
                values = values[:(len(hectimes) - st_offset)]  # trim end values
                hectimes = hectimes[:(len(hectimes) - st_offset)]  # trim end times

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))  # log failure
            sys.exit(-1)  # abort

        # if this record is in cms, convert it to cfs before summing
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')  # log conversion
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)  # cms -> cfs
            values = convvals  # replace with converted

        # the first record seeds the running total; subsequent records are added on
        if len(inflows) == 0:
            inflows = values  # seed sum
            times = hectimes  # store times (TODO: missing values handling)
        
        else:
            # add this record's values onto the running total, timestep by timestep
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
    Combine multiple time series by adding or subtracting each one
    according to a matching list of True/False operations.

    Supports optional unit conversion (for flow) and optional
    prepending of lookback values.

    Parameters
    ----------
    currentAlt : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging.
    timewindow : object
        The run time window object, used to determine the valid time
        range for trimming input records.
    inflow_records : list of str
        DSS paths of the records to combine. A path may optionally
        use the form `'other_dss_file::record_path'` to read from a
        different DSS file than `dss_file`.
    dss_file : str
        Path to the primary DSS file to read records from (used
        unless a record path specifies an alternate file via
        `'::'`).
    operation : list of bool
        One entry per record in `inflow_records`. True means add that
        record to the running total; False means subtract it. The
        first record's value is always used as-is to seed the total.
    output_dss_record_name : str
        DSS path to write the combined result to.
    output_dss_file : str
        Path to the DSS file to write the combined result to.
    what : str, optional
        If `"flow"` (default), values are treated as flow and
        converted from cms to cfs as needed, and the output units are
        set to `'CFS'`. If any other string, no unit conversion is
        performed and the output units are set to that string
        directly.
    prepend_n : int, optional
        If greater than 0, prepends this many copies of the first
        value (and corresponding earlier timestamps) onto the front
        of the combined series, e.g. to provide lookback values for
        downstream models. Defaults to 0.

    Returns
    -------
    None
        Writes the combined series to `output_dss_file`.

    Raises
    ------
    SystemExit
        Calls `sys.exit(-1)` if a `HecMathException` occurs while
        reading any of the input records.
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

    # Create lists to hold the data
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
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            # if the record's data extends past the run window, trim the trailing values
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])  # steps
                values = values[:(len(hectimes) - st_offset)]  # trim trailing values
                hectimes = hectimes[:(len(hectimes) - st_offset)]  # trim trailing times

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))  # log error
            sys.exit(-1)  # abort

        # if treating values as flow and this record is in cms, convert it to cfs
        if what=="flow":
            if units.lower() == 'cms':
                currentAlt.addComputeMessage('Converting cms to cfs')  # log conversion
                convvals = []

                for flow in values:
                    convvals.append(flow * 35.314666213)  # cms -> cfs

                values = convvals  # use converted values

        # the first record always seeds the running total directly; later records are
        # added or subtracted according to the matching operation flag
        if len(inflows) == 0:
            inflows = values  # seed series
            times = hectimes  # reference times

        else:
            if operation[j]:
                # add this record's values onto the running total, timestep by timestep
                for vi, v in enumerate(values):
                    inflows[vi] += v  # add series element-wise

            else:
                # subtract this record's values from the running total, timestep by timestep
                for vi, v in enumerate(values):
                    inflows[vi] -= v  # subtract element-wise

    # if requested, prepend lookback values (copies of the first value) onto the front of the series
    if prepend_n > 0:
        # add first value onto front of record
        # Sometimes ResSim needs some lookback values, or whatever
        p_times = [times[i] for i in range(len(times))]  # copy to list
        time_delta = times[1] - times[0]  # one-step delta
        times = [p_times[0] - time_delta * i for i in range(prepend_n, 0, -1)] + p_times  # prepend times
        values = [inflows[i] for i in range(len(inflows))]  # copy values list
        inflows = [values[0]] * prepend_n + values  # prepend first value

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
    """
    Create and write a DSS record with a constant value in it, for
    the given time window.

    Parameters
    ----------
    currentAlt : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging.
    timewindow : object
        The run time window object, used to determine the time span
        of the constant record (padded by 1 day on each end).
    output_dss_file : str
        Path to the DSS file to write the constant record to.
    constant : float, optional
        The constant value to fill the record with. Defaults to 0.0.
    what : str, optional
        The type of parameter being created. One of `'flow'`,
        `'temp-water'`, `'gate'`, `'evap'`, or `'elev'`
        (case-insensitive), which determines the units and parameter
        label used. Defaults to `'flow'`.
    dss_type : str, optional
        The DSS record type (e.g. `'PER-AVER'`, `'INST-VAL'`).
        Defaults to `'PER-AVER'`.
    period : str, optional
        The record's time interval. Must be `'1HOUR'` or `'1DAY'`
        (case-insensitive). Defaults to `'1HOUR'`.
    cpart : str, optional
        The C-part (location) label to assign to the record. Defaults
        to `'ZEROS'`.
    fpart : str, optional
        The F-part (version) label to assign to the record. Defaults
        to `'ZEROS'`.

    Returns
    -------
    bool
        True once the constant record has been written successfully.
        False if `what` or `period` is not a recognized value (in
        which case nothing is written).
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
        currentAlt.addComputeMessage('create_zero_dss_rec: what not known: %s' % what)  # log
        return False  # unsupported category

    # validate that the requested period is supported
    if period.lower()=='1hour':
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

    # generate a regular-interval series filled entirely with the constant value
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
    Compute relative humidity (%) from air temperature and dew point,
    using the August-Roche-Magnus approximation.

    Parameters
    ----------
    air_temp : float
        Air temperature in C.
    dewpoint_temp : float
        Dew point temperature in C.

    Returns
    -------
    float
        Relative humidity in percent, bounded to `[0.01, 100]`.

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
    Calculate dew point temperature given air temperature and
    relative humidity.

    Uses the algebraic inversion of the simplified
    August-Roche-Magnus approximation.

    Parameters
    ----------
    air_temp : float
        Air temperature in C.
    relative_humidity : float
        Relative humidity (0-100%).

    Returns
    -------
    float
        Dew point temperature in C.
    """

    gamma = math.log(relative_humidity / 100.0) + (17.62 * air_temp) / (243.12 + air_temp)  # intermediate term
    dewpoint = 243.12 * gamma / (17.62 - gamma)  # invert ARM
    
    return dewpoint  # dew point in C


def relhum_from_at_dp(met_dss_file, at_path, dp_path):
    """
    Compute relative humidity from air temperature and dew point
    records, writing the result as a new derived DSS record.

    Parameters
    ----------
    met_dss_file : str
        Path to the DSS file containing the source air temperature
        and dew point records, and where the derived humidity record
        will be written.
    at_path : str
        DSS path of the air temperature record.
    dp_path : str
        DSS path of the dew point record. Must align
        timestep-by-timestep with `at_path`.

    Returns
    -------
    None
        Writes the derived relative humidity record (with
        `'-DERIVED'` appended to the F-part and parameter set to
        `'RELHUM-FROM-AT-DP'`) to `met_dss_file`.
    """
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    dp_data = data_from_dss(met_dss_file, dp_path, None, None)
    # for each timestep, compute relative humidity from the paired air temp and dew point
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
    Compute dew point temperature from air temperature and relative
    humidity records, writing the result as a new derived DSS record.

    Parameters
    ----------
    met_dss_file : str
        Path to the DSS file containing the source air temperature
        and relative humidity records, and where the derived dew
        point record will be written.
    at_path : str
        DSS path of the air temperature record.
    rh_path : str
        DSS path of the relative humidity record. Must align
        timestep-by-timestep with `at_path`.

    Returns
    -------
    None
        Writes the derived dew point record (with `'-DERIVED'`
        appended to the F-part and parameter set to
        `'temp-dewpoint'`) to `met_dss_file`.
    """
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    rh_data = data_from_dss(met_dss_file, rh_path, None, None)
    # for each timestep, compute dew point from the paired air temp and relative humidity
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
    Trim a values/times pair so that the data falls entirely within
    a given start/end time window.

    Drops any leading or trailing values outside that range.

    Parameters
    ----------
    values : list
        The data values to trim.
    times : list
        The corresponding raw HEC time values, in the same order as
        `values`.
    startime : int or float
        The window start time, in the same units as `times`.
    endtime : int or float
        The window end time, in the same units as `times`.

    Returns
    -------
    tuple of (list, list)
        `(values, times)`, each trimmed to fall within
        `[startime, endtime]`.
    """
    if times[0] < startime:  # if startdate is before the timewindow..
        print('start date ({0}) from DSS before timewindow ({1})..'.format(times[0], startime))
        st_offset = (startime - times[0]) / (times[1] - times[0])
        values = values[st_offset:]
        times = times[st_offset:]
    # if the data extends past the window, trim the trailing values
    if times[-1] > endtime:
        print('end date ({0}) from DSS after timewindow ({1})..'.format(times[-1], endtime))
        st_offset = (times[-1] - endtime) / (times[1] - times[0])  # number of steps to trim
        values = values[:(len(times) - st_offset)]  # drop trailing values
        times = times[:(len(times) - st_offset)]  # drop trailing times
    
    # Return trimmed outputs to the calling function
    return values, times


def replace_data(currentAlt, timewindow, pairs, dss_file, dss_outfile, months, standard_interval=None):
    """
    Replace base record values with alternate record values during
    specified months, for each of several record pairs.

    For each (base, alternate) record pair, replace the base record's
    values with the alternate record's values during a specified set
    of months, and write the merged result to a new output record.

    Parameters
    ----------
    currentAlt : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging.
    timewindow : object
        The run time window object, used to determine the read window
        for both records in each pair.
    pairs : list of list
        A list of `[base_path, alt_path]` DSS path pairs. For each
        pair, values from `alt_path` replace values from `base_path`
        during the given months.
    dss_file : str
        Path to the DSS file containing both the base and alternate
        records.
    dss_outfile : str
        Path to the DSS file to write the merged result to.
    months : list of int
        The calendar months (1-12) during which the alternate
        record's values should replace the base record's values.
    standard_interval : str, optional
        If provided, both the base and alternate series are
        standardized to this interval (via `standardize_interval`)
        before merging. Defaults to None (no standardization).

    Returns
    -------
    None
        Writes one merged record per pair in `pairs` to `dss_outfile`,
        with the F-part/A-part renamed to indicate the merge source.

    Raises
    ------
    SystemExit
        Exits the process (`sys.exit(1)`) if the base and alternate
        records in a pair have mismatched units or mismatched
        intervals.
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    # 01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    # for each base/alternate record pair, merge them and write the result
    for pair in pairs:
        # Get the values from the series to be filled
        # Open the DSS file
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

        # the two records must share units to be merged meaningfully
        if base_units != alt_units:
            currentAlt.addComputeMessage('Units do not match for {0} and {1}, skipping'.format(pair[0], pair[1]))  # log mismatch
            dssFm.close()
            sys.exit(1)  # abort
        
        # Compare the intervals between series and enforce that they are the same
        if base_interval != alt_interval:
            currentAlt.addComputeMessage('Intervals do not match for {0} and {1}, changing interval...'.format(pair[0], pair[1]))  # log mismatch
            dssFm.close()
            sys.exit(1)  # abort

        # for each timestep in the base record, if its month is in the replacement list,
        # swap in the alternate record's value instead
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

def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    """
    Apply a fixed temperature lapse (offset) to an air temperature
    record.

    Approximates temperature change with elevation, writing the
    adjusted record under a new F-part. This is the active (final)
    definition of `airtemp_lapse` used at runtime, superseding the
    earlier identical definition above.

    Parameters
    ----------
    dss_file : str
        Path to the DSS file containing the source air temperature
        record.
    dss_rec : str
        DSS path of the source air temperature record.
    lapse_in_C : float
        The temperature lapse/offset to apply, in degrees Celsius. If
        the source record's units are Fahrenheit, this is converted
        to an equivalent Fahrenheit offset before being applied.
    dss_outfile : str
        Path to the DSS file to write the adjusted record to.
    f_part : str
        The F-part to assign to the output record.

    Returns
    -------
    None
        Writes the lapse-adjusted record to `dss_outfile`.

    Notes
    -----
    This is a duplicate definition of `airtemp_lapse`, identical to
    the one defined earlier in this module. Because Python/Jython
    executes definitions top-to-bottom, this later definition
    silently replaces the earlier one - any call to `airtemp_lapse`
    uses this version. Consider removing one of the two copies.
    """
    dss = HecDss.open(dss_file)
    
    # Read the timeseries
    tsm = dss.read(dss_rec)
    
    # Create a copy of the data
    lapse = lapse_in_C
    # if the source series is in Fahrenheit, convert the lapse amount to an equivalent Fahrenheit offset
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
    Prepend the first value of a series `prepend_n` times, for
    lookbacks, and write back to the same record.

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
        Writes the extended series back to `dss_rec` in `dss_file`.
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

def postprend_last_value_on_ts(dss_file,dss_rec,postpend_n):
    """
    Append copies of the last value onto the end of a time series.

    Sometimes ResSim needs some lookback values, or whatever has to
    be a "regular" record.

    Parameters
    ----------
    dss_file : str
        Path to the DSS file containing the source record, which is
        also modified in place.
    dss_rec : str
        DSS path of the record to modify. Must be a regular-interval
        record.
    postpend_n : int
        Number of copies of the last value (and corresponding later
        timestamps) to append onto the end of the record.

    Returns
    -------
    None
        Overwrites `dss_rec` in `dss_file` with the extended series.
    """
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
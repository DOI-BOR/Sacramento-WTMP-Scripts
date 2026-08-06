"""
Balance Flow & Elevation Utilities (WTMP / HEC-WAT)
===================================================

Purpose
-------
Provide helper functions to:
- Read and aggregate inflow/outflow records from DSS.
- Convert units and trim series to a specified runtime window.
- Compute predicted elevation and storage from inflow/outflow and stage.
- Create reservoir balance-flow time series suitable for ResSim/W2 workflows.

The module supports conic storage interpolation (consistent with HEC ResSim)
and linear interpolation as needed. It writes outputs back to DSS with
appropriate metadata (interval, units, type) and can optionally transform
to a different time step for alternate model contexts.

Notes
-----
- **Environment:** Jython in HEC-WAT using Python 2.7 syntax; relies on HEC-DSS APIs, 
  HecMath tools, and WAT runtime constructs (e.g., time windows, project directories).
- **Units/Conversions:** Flow conversion uses CMS→CFS factor `35.314666213`.
  Stage may be converted from meters to feet; storage conversions use
  `cfs_2_acreft`/`acreft_2_cfs` derived from the balance period.
- **Paths:** DSS paths may be provided as `file::/A/B/C/D/E/F/` to read from a
  specific DSS file; otherwise the primary DSS file is used.

See Also
--------
HecDss
    Open/read/write DSS files.
TimeSeriesContainer
    Container for time series writes.
HecTime
    Utility for HEC time values and conversions.
"""

import math                                   # Standard library math utilities (sqrt, etc.) used in storage/elevation calcs
from hec.heclib.dss import HecDss             # HEC-DSS API: open/read/write .dss files
from hec.hecmath import HecMathException      # HEC math exceptions: catch read/transform failures
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE  # HEC sentinel for undefined/missing numeric values
from hec.heclib.util import HecTime           # HEC time utility: string↔HEC time value conversions
from hec.io import DSSIdentifier              # HEC I/O identifier (imported for completeness; not directly used here)
from hec.io import TimeSeriesContainer        # HEC I/O container used to assemble series for writing
import hec.hecmath.TimeSeriesMath as tsmath   # HEC math utilities for transforming time series (e.g., intervals)
from rma.util.RMAConst import MISSING_DOUBLE  # RMA constant for missing values (imported for consistency)
import math                                   # Duplicate import retained intentionally (original code unchanged)
import sys                                    # Standard library: sys.exit on fatal conditions (e.g., interval mismatch)
import datetime as dt                         # Standard library datetime (not directly used but retained)
import os                                     # Standard library: filesystem operations for CSV output and paths
from hec.io import DSSIdentifier              # Duplicate import retained to preserve original file content
from hec.heclib.util import HecTime           # Duplicate import retained to preserve original file content
from com.rma.io import DssFileManagerImpl     # RMA DSS manager (not used directly here but kept for context)
from java.util import TimeZone                # Java TimeZone (not used directly; retained for completeness)


def linear_interpolation(x_values, y_values, x):
    """
    Linearly interpolate `y` at position `x` given ordered (x, y) vectors.

    Parameters
    ----------
    x_values : list of float
        Monotonic sequence of x-values (independent variable).
    y_values : list of float
        Corresponding y-values (dependent variable).
    x : float
        Interpolation point within the domain of `x_values`.

    Returns
    -------
    float
        Interpolated y-value at `x`.

    Raises
    ------
    ValueError
        If input lists differ in length or have fewer than two points, or if
        `x` lies outside the range of `x_values`.

    Notes
    -----
    - Performs simple linear interpolation over the first segment where
      `x <= x_values[i]` is satisfied.
    """

    # Peform initial data validation to confirm series length
    if len(x_values) != len(y_values) or len(x_values) < 2:
        raise ValueError("Input lists must have the same length and contain at least 2 data points.")

    # Scan segments to find the interval containing x; assumes x_values are ordered
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
    Read elevation-storage-area file into lists for interpolation and conic methods.

    Parameters
    ----------
    file_name : str
        Path to CSV-like data file. Expected columns differ by reservoir:
        - For ``Natoma``: `[elev, area]`
        - Otherwise: `[elev, stor, area]`
    res_name : str
        Reservoir name; used to select parsing behavior (`Natoma` special case).

    Returns
    -------
    dict
        Dictionary with keys: ``'elev'``, ``'stor'``, ``'area'``, each a list of floats.

    Notes
    -----
    - Units expected: elevation (ft), storage (acre-ft), area (acre).
    - Order of rows is preserved; caller should ensure monotonic elevation.
    """
    
    # These are in [elev, stor, area] with units [ft, acre-ft, acre]
    elevstorarea = {} #avoid lists doing weird things like mixing up order..
    elev = []
    stor = []
    area = []
    
    print('cwd: ' + os.getcwd())             # Useful diagnostic for file relative paths

    # Natoma has only elev, area; others include storage
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

    # Package into dict (lists aligned by row order)
    elevstorarea['elev'] = elev
    elevstorarea['stor'] = stor
    elevstorarea['area'] = area
    
    # Return the storage elevation curve
    return elevstorarea


def build_conic_storage_array(elev, area, firstStorageValue=0.0):
    """
    Compute cumulative storage array using conic estimation between elevation-area points.

    Adapted from HEC ResSim (storage.java, 2022-06-17). Each step integrates a
    conic segment between successive elevation-area points and accumulates.

    Parameters
    ----------
    elev : list of float
        Elevation measurement points (ft).
    area : list of float
        Surface area at each elevation (acre).
    firstStorageValue : float, optional
        Initial cumulative storage (acre-ft), by default ``0.0``.

    Returns
    -------
    list of float
        Cumulative storage at each elevation point (acre-ft).
    """

    # Calculate storage at each elevation using conic formula
    n_measures = len(elev)
    storage = []
    storage.append(firstStorageValue)
    
    # Loop on each of the points in the
    for i in range(1, n_measures):
        h = elev[i] - elev[i-1]                     # Elevation difference between successive points
        storage.append(h/3. * (area[i-1] + area[i] + math.sqrt(area[i-1] * area[i])) + storage[i-1])
    
    # Return the storage curve
    return storage


def conic_storage_interp(interpElev, elev, area, conicStorage, idx):
    """
    Interpolate storage between measurement points using conic layer approach.

    Parameters
    ----------
    interpElev : float
        Elevation at which to interpolate storage (ft).
    elev : list of float
        Elevation measurement points (ft).
    area : list of float
        Surface areas corresponding to `elev` (acre).
    conicStorage : list of float
        Cumulative conic storage array from ``build_conic_storage_array``.
    idx : int
        Lower-bounding index for `interpElev` within `elev`.

    Returns
    -------
    float
        Interpolated storage at `interpElev` (acre-ft).

    Notes
    -----
    - Uses geometric mean of adjacent areas per conic formulation and accumulates
      within the interval `[idx, idx+1]`.
    """

    # Normalize elevation position within layer
    h = (interpElev - elev[idx]) / (elev[idx+1] - elev[idx])  
    
    # Geometric mean area for conic
    geomMean = math.sqrt(area[idx] * area[idx+1])             
    
    # Interpolate for the area and the storage
    areaInterp = area[idx] + 2.*(geomMean - area[idx])*h + (area[idx] + area[idx+1] - 2.*geomMean)*h*h
    storageInterp = (interpElev - elev[idx])/3. * (area[idx] + areaInterp + math.sqrt(area[idx] * areaInterp)) + conicStorage[idx]
    
    # Return the interpolated storage
    return storageInterp


def get_elev_layer_idx(elev, obs_elev, elev_stor_area):
    """
    Find lower-bounding index in the elevation array for a given observed elevation.

    Parameters
    ----------
    elev : list of float
        Elevation reference points (ft).
    obs_elev : float
        Observed elevation to index (ft).
    elev_stor_area : dict
        Elevation-storage-area dictionary (see ``read_elev_storage_area_file``).

    Returns
    -------
    int
        Lower-bounding index for `obs_elev` within `elev`; may be decremented if
        the closest point is above `obs_elev`. Returns ``-1`` if not found.

    Notes
    -----
    - Uses absolute difference to find the nearest elevation then adjusts downward
      if that elevation exceeds the observed elevation.
    """

    # Initialize with HEC sentinel
    idx = UNDEFINED_DOUBLE  

    # Track minimum absolute difference 
    min_val = None                           
    
    # Loop over each point in the elevation
    for i in range(len(elev)):
        # Absolute difference to find nearest layer
        valchk = abs(elev[i]-obs_elev)       
        
        # Take action mased on the magnitude of the value
        if math.isnan(valchk):
            min_val = valchk
            idx = i
        
        elif min_val == None:
            min_val = valchk
            idx = i
        
        elif valchk < min_val:
            min_val = valchk
            idx = i

    # Adjust index if the nearest elevation is above observed elevation
    if idx != UNDEFINED_DOUBLE:
        if elev_stor_area['elev'][idx] > obs_elev: #TODO: is this multidimensional?
            idx -= 1
    
    else:
        # Fallback if index could not be determined
        idx = -1                                
    
    # Return to the calling function
    return idx


def get_balance_period(balance_period):
    """
    Convert a ResSim/WTMP balance-period string to hours (float).

    Parameters
    ----------
    balance_period : str
        Time step string (e.g., ``'1Hour'``, ``'1Day'``, ``'30Min'``).

    Returns
    -------
    float
        Duration in **hours**.

    Notes
    -----
    - Case-insensitive parsing; valid tokens include ``hour``, ``day``, ``min``.
    """
    
    if 'hour' in balance_period.lower():
        return float(balance_period.lower().replace('hour', ''))
    
    elif 'day' in balance_period.lower():
        return float(balance_period.lower().replace('day', '')) * 24
    
    elif 'min' in balance_period.lower():
        return float(balance_period.lower().replace('min', '')) / 60


def check_dss_intervals(records, balance_period, currentAlt):
    """
    Validate that each DSS record pathname contains the expected balance-period token.

    Parameters
    ----------
    records : list of str
        DSS pathnames to inspect (e.g., ``'/...//1Hour/...'``).
    balance_period : str
        Expected time-step token (e.g., ``'1Hour'``).
    currentAlt : object
        WAT scripting alternative used for logging and termination.

    Raises
    ------
    SystemExit
        If any record does not contain the expected token (case-insensitive).
    """
    
    for r in records:
        if balance_period.lower() not in r.lower():
            currentAlt.addComputeMessage('DSS record {0} not matching time interval {1}'.format(r, balance_period))
            sys.exit(-1)                      # Terminate compute if intervals mismatch


def read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str):
    """
    Read a time series, optionally from an alternate DSS file prefix.

    Parameters
    ----------
    dssFm : HecDss
        Open HecDss instance for default reads.
    pathname : str
        DSS pathname; may include prefix like ``'file.dss::/A/B/C/D/E/F/'``.
    starttime_str : str
        Start time string for the read (HEC date format).
    endtime_str : str
        End time string for the read (HEC date format).

    Returns
    -------
    hec.io.TimeSeriesContainer
        Data container with ``times``, ``values``, ``units``, etc.

    Notes
    -----
    - When an alternate file prefix is present, a new HecDss handle is opened,
      the series is read, and the handle is closed before returning the data.
    """

    if '::' in pathname:
        print('Splitting and reading:',pathname)        # Diagnostic for alternate-file reads
        alt_dss_file,pathname_clean = pathname.split('::')
        dssFmRec = HecDss.open(alt_dss_file)            # Open alternate DSS
        tsc = dssFmRec.read(pathname_clean, starttime_str, endtime_str, False).getData()
        dssFmRec.close()                                 # Close alternate DSS handle
    
    else:
        tsc = dssFm.read(pathname, starttime_str, endtime_str, False).getData()
    
    return tsc


def read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, starttime_str, endtime_str,
                          starttime_hectime, endtime_hectime):
    """
    Read and aggregate inflows and outflows, trimming series to the runtime window.

    Parameters
    ----------
    currentAlt : object
        WAT scripting alternative for logging.
    dss_file : str
        Primary DSS filename used when no alternate prefix is present.
    inflow_records : list of str
        DSS pathnames (or prefixed paths) for inflow components.
    outflow_records : list of str
        DSS pathnames (or prefixed paths) for outflow components.
    starttime_str : str
        Start time string for reads.
    endtime_str : str
        End time string for reads.
    starttime_hectime : int
        HEC time (int) corresponding to ``starttime_str``.
    endtime_hectime : int
        HEC time (int) corresponding to ``endtime_str``.

    Returns
    -------
    (list, list)
        ``times`` (HEC time vector) and ``inflow_outflow`` (cfs), where
        ``inflow_outflow[i] = inflows[i] - outflows[i]`` after unit conversion.

    Notes
    -----
    - CMS flows are converted to CFS.
    - Series are trimmed if their range exceeds the runtime window.
    """
    
    # Open primary DSS for reads
    dssFm = HecDss.open(dss_file)            

    inflows = []
    outflows = []
    times = []                                # Shared time vector from the first inflow series

    # Read inflows
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        print('\nreading: ' + str(pathname))
        try:
       
            print(starttime_str, endtime_str)           # Diagnostic prints for window
            print(dss_file)                             # Diagnostic: file being read
            tsc = read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str)
            values = tsc.values                         # Flow values (may be CMS)
            hectimes = tsc.times                        # HEC time stamps
            units = tsc.units                           # Units for possible conversion

            # Trim leading segment if it starts before the time window
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]

            # Trim trailing segment if it ends after the time window
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)                                 # Fatal read error

        # Convert CMS→CFS if needed
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting inflow to cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # Aggregate inflows across records; capture time vector from the first record
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
            values = tsc.values                           # Outflow values (may be CMS)
            hectimes = tsc.times                          # HEC time stamps
            units = tsc.units                             # Units for possible conversion

            # Trim leading segment if necessary
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], endtime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]

            # Trim trailing segment if necessary
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
          
        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)                                  # Fatal read error

        # Convert CMS→CFS if needed
        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            print('Converting outflow cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        # Aggregate outflows across records
        if len(outflows) == 0:
            outflows = values
        else:
            for vi, v in enumerate(values):
                outflows[vi] += v

    # Close DSS handle after reads
    dssFm.close()                                        

    # Inflow minus outflow record
    inflow_outflow = []
    for i in range(len(inflows[1:])):
        inflow_outflow.append(inflows[i+1] - outflows[i+1])  # Element-wise inflow - outflow (cfs)

    # Log diagnostic lengths to WAT
    currentAlt.addComputeMessage("Len inflow_outflow:"+str(len(inflow_outflow)))
    currentAlt.addComputeMessage("Len times:"+str(len(times)))

    # Return to the calling function
    return times,inflow_outflow


def predict_elevation(currentAlt, starttime_str, endtime_str, res_name, inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Hour'):
    """
    Predict elevation and storage over a period from inflow/outflow balance.

    Useful for generating lookback/starting elevations for forecast runs that
    begin mid-period; integrates inflow-outflow into storage and maps storage
    back to elevation via interpolation.

    Parameters
    ----------
    currentAlt : object
        WAT scripting alternative for logging.
    starttime_str : str
        Start time string for DSS reads.
    endtime_str : str
        End time string for DSS reads.
    res_name : str
        Reservoir name (used in CSV output naming and context).
    inflow_records : list of str
        DSS pathnames for inflow components (must match `balance_period_str`).
    outflow_records : list of str
        DSS pathnames for outflow components (must match `balance_period_str`).
    starting_elevation : float
        Initial elevation at `starttime_str` (ft).
    elev_stor_area : dict
        Elevation-storage-area dataset (lists) for interpolation/conic methods.
    dss_file : str
        DSS file for reading inflow/outflow.
    output_dss_record_name : str
        DSS pathname for predicted elevation output.
    output_dss_file : str
        DSS file where predicted elevation/storage will be written.
    shared_dir : str
        Directory used for auxiliary CSV outputs.
    use_conic : bool, optional
        If True, use conic storage interpolation; otherwise linear.
    alt_period : int, optional
        Alternate period length to transform output time step (e.g., minutes).
    alt_period_string : str, optional
        Alternate period string (e.g., ``'1Day'``) for transform.
    balance_period_str : str, optional
        Balance period token (default ``'1Hour'``).

    Returns
    -------
    None
        Writes predicted elevation and storage records to DSS.

    Notes
    -----
    - Storage integration uses `cfs_2_acreft` derived from balance period.
    - Elevation mapping uses linear interpolation on storage→elevation curve.
    """
    '''From inflows/outflows, predict hourly elevation, useful for lookback/starting elevation for forecasts starting
    on arbitrary dates during forecast period
    '''

    # Get the balance period of the calculation
    balance_period = get_balance_period(balance_period_str) # convert to (float) hours
    
    # Confirm the DSS values support the adjustment
    check_dss_intervals(inflow_records, balance_period_str, currentAlt)  # Validate inflow intervals
    check_dss_intervals(outflow_records, balance_period_str, currentAlt) # Validate outflow intervals
    
    # Perform unit conversions
    cfs_2_acreft = balance_period * 3600. / 43559.9         # Convert cfs-period to acre-ft
    acreft_2_cfs = 1. / cfs_2_acreft                        # Inverse conversion

    # Convert the datetimes into HEC format
    starttime_hectime = HecTime(starttime_str).value()       # HEC integer time for start
    endtime_hectime = HecTime(endtime_str).value()           # HEC integer time for end

    # Read inflow-outflow net series over the same window
    times,inflow_outflow = read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, 
                                                 starttime_str, endtime_str, starttime_hectime, endtime_hectime)

    # Log diagnostics for sanity-check
    currentAlt.addComputeMessage("Len inflow_outflow:"+str(len(inflow_outflow)))
    currentAlt.addComputeMessage("Len times:"+str(len(times)))
    
    # TODO: support conic interpolation
    # TODO: support evap, but really that's just a positive outflow....
    storage = linear_interpolation(elev_stor_area['elev'], elev_stor_area['stor'], starting_elevation)  # Initial storage from starting elev
    storage = [storage,]                                           # Begin cumulative storage sequence
    elev_predicted = []                                            # Predicted elevations to be written

    # Integrate inflow-outflow over the period and map storage→elevation
    for i in range(len(inflow_outflow)):
        storage.append( storage[-1] + inflow_outflow[i]*cfs_2_acreft )
        elev_predicted.append( linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], storage[-1]) )

    # Output record
    dssFm_out = HecDss.open(output_dss_file)                       # DSS for writing predicted series
    steptime = times[1]-times[0]                                   # Step size in HEC time units

    # Format the timeseries container
    tsc = TimeSeriesContainer()
    tsc.startTime = times[0] #- steptime                           # Start at the first time
    tsc.interval = int(balance_period)*60                          # Interval in minutes
    tsc.fullName = output_dss_record_name                          # Predicted elevation path
    tsc.values = [starting_elevation] + elev_predicted             # Include initial elevation at t0
    #tsc.startTime = times[1]
    tsc.units = 'ft'                                               # Elevation units
    tsc.type = 'INST-VAL'                                          # Instantaneous values
    tsc.numberValues = len(tsc.values)
    dssFm_out.write(tsc)                                           # Write elevation prediction

    # Also write predicted storage under a modified record name
    recparts = output_dss_record_name.split('/')
    recparts[3] = 'STORAGE-PREDICTED'                              # Replace E-part with storage label
    tsc.startTime = times[0]                                       # Storage start time
    tsc.fullName = '/'.join(recparts)
    tsc.values = storage                                           # Cumulative storage
    tsc.units = 'ac-ft'                                            # Storage units
    tsc.type = 'PER-CUM'                                           # Cumulative per-period
    dssFm_out.write(tsc)

    # Optional: transform to alternate period if requested
    if alt_period is not None:
        if alt_period_string.lower() != balance_period_str.lower():
            tsm = dssFm_out.read(output_dss_record_name)
            tsm_new_interval = tsm.transformTimeSeries(alt_period_string, "", "AVE")
            dssFm_out.write(tsm_new_interval)

    dssFm_out.close()                                             # Close DSS handle after writes


def create_balance_flows(currentAlt, timewindow, res_name, inflow_records, outflow_records, stage_record, evap_record,
                         elev_stor_area, dss_file, output_dss_record_name, output_dss_file, shared_dir,
                         storage_dss_record_name='', evap_dss_record_name='',
                         balance_period_str="1HOUR", use_conic=False, write_evap=False, write_storage=False,
                         alt_period=None,alt_period_string=None, lookback_padding=1440):
    """
    Compute reservoir balance flow series from inflow/outflow, stage, and evaporation.

    Parameters
    ----------
    currentAlt : object
        WAT scripting alternative for logging.
    timewindow : object
        Runtime window object; provides start/end strings and HEC times.
    res_name : str
        Reservoir name used in CSV output and context.
    inflow_records : list of str
        DSS inflow pathnames (may include file prefixes).
    outflow_records : list of str
        DSS outflow pathnames (may include file prefixes).
    stage_record : str
        DSS pathname for stage/elevation series (INST-VAL).
    evap_record : str
        DSS pathname for evaporation depth (per-period accumulation).
    elev_stor_area : dict
        Elevation-storage-area dataset from file reader.
    dss_file : str
        Primary DSS filename for reading series.
    output_dss_record_name : str
        DSS pathname for balance-flow output.
    output_dss_file : str
        DSS filename where outputs will be written.
    shared_dir : str
        Directory used for CSV diagnostics.
    storage_dss_record_name : str, optional
        DSS pathname for optional storage output.
    evap_dss_record_name : str, optional
        DSS pathname for optional evaporation-flow output.
    balance_period_str : str, optional
        Balance period token (default ``"1HOUR"``).
    use_conic : bool, optional
        If True, use conic method; else linear interpolation for storage.
    write_evap : bool, optional
        If True, write evaporation flow series.
    write_storage : bool, optional
        If True, write storage series.
    alt_period : int, optional
        Alternate period length for time transform.
    alt_period_string : str, optional
        Alternate period string (e.g., ``"1Day"``).
    lookback_padding : int, optional
        Padding (minutes) for lookback (currently commented out).

    Returns
    -------
    bool
        True on success.

    Notes
    -----
    - Trims stage/evap series to runtime window; converts stage from meters to feet if needed.
    - Applies CMS→CFS conversion to inflow/outflow in reader helper.
    - Writes CSV diagnostics and DSS balance-flow record; filters spurious values at edges.
    """

    # Validate intervals before processing (defensive check)
    check_dss_intervals(inflow_records, balance_period_str, currentAlt)
    check_dss_intervals(outflow_records, balance_period_str, currentAlt)
    check_dss_intervals([stage_record, evap_record], balance_period_str, currentAlt)
    
    balance_period = get_balance_period(balance_period_str) # convert to (float) hours
    print('balance_period ' + str(balance_period))          # Diagnostic in console
    
    cfs_2_acreft = balance_period * 3600. / 43559.9         # Convert period-avg cfs to acre-ft
    acreft_2_cfs = 1. / cfs_2_acreft                        # Inverse conversion for later use

    starttime_str = timewindow.getStartTimeString()         # Window start string for reads
    endtime_str = timewindow.getEndTimeString()             # Window end string for reads

    starttime_hectime = HecTime(starttime_str).value()      # HEC int start time
    endtime_hectime = HecTime(endtime_str).value()          # HEC int end time
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    # Read net inflow-outflow series aligned to time window
    times,inflow_outflow = read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, 
                                                 starttime_str, endtime_str, starttime_hectime, endtime_hectime)
    print('len times:',len(times))                          # Diagnostic
    print('len inflow_outflow:',len(inflow_outflow))        # Diagnostic

    dssFm = HecDss.open(dss_file)                           # DSS for stage/evap reads

    # Read stage
    print('Reading stage')
    tsc = read_ts_rec_w_optional_fname(dssFm, stage_record, starttime_str, endtime_str)
    
    try:
        stage = tsc.values                                  # Elevation series
        hectimes = tsc.times                                # Elevation timestamps
        
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

        # Convert meters→feet if required
        if tsc.units.lower() == 'm':
            currentAlt.addComputeMessage('Converting stage m to ft')
            print('Converting stage cms to cfs')            # (Message retained as in original)
            convvals = []
            for elev in stage:
                convvals.append(elev * 3.280839895)
            stage = convvals
        
    except HecMathException:
        currentAlt.addComputeMessage('ERROR reading' + str(stage_record))
        sys.exit(-1)                                        # Fatal read error

    # Read evap
    print('Reading evap')
    tsc = read_ts_rec_w_optional_fname(dssFm, evap_record, starttime_str, endtime_str)
    
    try:
        evap = tsc.values                                   # Evaporation depth per period
        hectimes = tsc.times                                # Timestamps for evaporation
        
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
        sys.exit(-1)                                        # Fatal read error

    # Build conic storage array for interpolation later
    conic_storage = build_conic_storage_array(elev_stor_area['elev'], elev_stor_area['area'])

    # Calculations
    n = len(stage) - 1                                      # Number of intervals
    flow_resid = []                                         # Residual flow (cfs)
    flow_evap = []                                          # Evaporation flow (cfs)
    storage_record = []                                     # Optional storage record (start of each interval)

    for k in range(n):
        stage_start = stage[k]
        stage_end = stage[k+1]

        # Determine storage change using conic or linear interpolation
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

        # Average area across endpoints for evap flow calculation
        area_avg = 0.5 * (linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_start) +
                          linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_end))
        evap_flow_loss = (evap[k] * area_avg) * acreft_2_cfs  # in cfs

        # Residual required to balance storage change with net inflow minus evap
        resid = delta_stor_flow - (inflow_minus_outflow - evap_flow_loss)
        flow_resid.append(resid)
        flow_evap.append(evap_flow_loss)
        storage_record.append(storage_start)                  # Record storage at interval start


    if True:
        print('Writing to CSV')
        # dump to CSV if DSS is mis-behaving or if flows are -999. etc.
        with open(os.path.join(shared_dir, "{0}_balance_flow.csv".format(res_name)), 'w') as opf:
            opf.write('date, balance_flow [cfs]\n')
            
            for i in range(len(flow_resid)):
                new_line = ','.join([str(times[i]), str(flow_resid[i]), '\n'])
                opf.write(new_line)

    dssFm_out = HecDss.open(output_dss_file)                 # DSS for balance-flow writes
    
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
                flow_resid[idx] = 0.0                        # Zero-out unrealistic edge values

    steptime = times[1]-times[0]                             # HEC step delta
    tsc = TimeSeriesContainer()
    tsc.startTime = times[0] - steptime                      # Start one step earlier (ResSim compatibility)
    tsc.interval = int(balance_period)*60                    # Interval in minutes
    tsc.fullName = output_dss_record_name                    # Balance-flow output path

    # copy back 1st balance flow record 2 steps, instead of writing from 1st valid balance calc.
    # otherwise, time-averaging the balanece flows later leaves off the 1st time step needed for a ResSim run
    # best we can do I guess to make ResSim computes work
    tsc.values = [flow_resid[0],flow_resid[0]] + flow_resid  # Prepend duplicates for ResSim averaging behavior
    tsc.units = 'CFS'                                        # Flow units
    tsc.type = 'PER-AVER'                                    # Period-averaged values
    tsc.numberValues = len(tsc.values)
    dssFm_out.write(tsc)                                     # Write balance flow series

    # Optional alternate time transform if requested
    if alt_period is not None:
        if alt_period_string.lower() != balance_period_str.lower():
            tsm = dssFm_out.read(output_dss_record_name)
            tsm_new_interval = tsm.transformTimeSeries(alt_period_string, "", "AVE")
            dssFm_out.write(tsm_new_interval)

    if write_evap:
        tsc = TimeSeriesContainer()
        tsc.times = times                                    # Use original times for evap flow output
        tsc.fullName = evap_dss_record_name
        tsc.values = flow_evap
        tsc.startTime = times[1]                             # Align with first valid interval
        tsc.units = 'CFS'
        tsc.type = 'PER-AVER'
        tsc.endTime = times[-1]
        tsc.numberValues = len(flow_resid)
        tsc.startHecTime = timewindow.getStartTime()
        tsc.endHecTime = timewindow.getEndTime()
        dssFm_out.write(tsc)

    if write_storage:
        tsc = TimeSeriesContainer()
        tsc.times = times                                    # Use original times for storage output
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

    dssFm.close()                                            # Close input DSS
    dssFm_out.close()                                        # Close output DSS
    return True                                              # Indicate success
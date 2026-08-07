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
- **Units/Conversions:** Flow conversion uses CMS to CFS factor `35.314666213`.
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
from hec.heclib.util import HecTime           # HEC time utility: string to HEC time value conversions
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
    Estimate a y-value at a given x by linear interpolation between
    known (x, y) points.

    Parameters
    ----------
    x_values : list of float
        Known x-coordinates, assumed to be sorted in ascending order.
    y_values : list of float
        Known y-coordinates, parallel to `x_values`.
    x : float
        The x-value to interpolate a y-value for. Must fall within the
        range of `x_values`.

    Returns
    -------
    float
        The linearly interpolated y-value at `x`.

    Raises
    ------
    ValueError
        If `x_values` and `y_values` differ in length or contain fewer
        than 2 points, or if `x` falls outside the range covered by
        `x_values`.

    Notes
    -----
    Performs a linear forward scan (not a binary search) to find the
    bracketing segment, which is fine for the small lookup tables used
    in this module but would be inefficient for very large tables.
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
    Read a CSV elevation-storage-area curve for a reservoir into
    parallel lists.

    Parameters
    ----------
    file_name : str
        Path to the CSV file to read. Expected to contain one row per
        elevation point, comma-separated.
    res_name : str
        Reservoir name. If (case-insensitive) equal to ``'natoma'``,
        the file is parsed as two columns (elev, area). Otherwise, it
        is parsed as three columns (elev, stor, area).

    Returns
    -------
    dict
        Dictionary with keys ``'elev'``, ``'stor'``, and ``'area'``,
        each mapping to a list of floats parsed from the file, in
        units of feet, acre-feet, and acres respectively. For
        ``'natoma'``, the ``'stor'`` list will be empty, since that
        file format does not include a storage column.

    Notes
    -----
    The ``'natoma'`` reservoir uses a special-cased two-column file
    format (elevation, area only - no storage column), while all
    other reservoirs are expected to have a three-column format
    (elevation, storage, area). Calling this with ``res_name='natoma'``
    against a three-column file (or vice versa) will misread the data.
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

    # Package into dict (lists aligned by row order)
    elevstorarea['elev'] = elev
    elevstorarea['stor'] = stor
    elevstorarea['area'] = area
    
    # Return the storage elevation curve
    return elevstorarea


def build_conic_storage_array(elev, area, firstStorageValue=0.0):
    """
    Compute cumulative storage array using conic estimation between
    elevation-area points.

    Adapted from HEC ResSim (storage.java, 2022-06-17). Each step
    integrates a conic segment between successive elevation-area
    points and accumulates.

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

    Notes
    -----
    The conic (frustum) formula estimates the volume of the slab
    between two elevation-area points as
    ``h/3 * (A0 + A1 + sqrt(A0*A1))``, which is exact for a linearly
    tapering (conic) shape and is the same approach used internally
    by HEC-ResSim.
    """

    # Calculate storage at each elevation using conic formula
    n_measures = len(elev)
    storage = []
    storage.append(firstStorageValue)
    # for each elevation layer, add the conic-estimated storage volume of that slab
    # onto the cumulative total from the layer below
    for i in range(1, n_measures):
        h = elev[i] - elev[i-1]                     # Elevation difference between successive points
        storage.append(h/3. * (area[i-1] + area[i] + math.sqrt(area[i-1] * area[i])) + storage[i-1])
    
    # Return the storage curve
    return storage


def conic_storage_interp(interpElev, elev, area, conicStorage, idx):
    """
    Interpolate storage between measurement points using conic layer
    approach.

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
    Find the index of the elevation layer that lower-bounds a given
    observed elevation, within an elevation-storage-area table.

    Parameters
    ----------
    elev : list of float
        Known elevation values to search within.
    obs_elev : float
        The observed elevation to locate within `elev`.
    elev_stor_area : dict
        The elevation-storage-area table (as produced by
        `read_elev_storage_area_file`), used to re-check whether the
        closest elevation is above or below `obs_elev`.

    Returns
    -------
    int
        The index of the elevation value that lower-bounds `obs_elev`
        (i.e. ``elev[idx] <= obs_elev``). Returns -1 if no valid index
        could be determined (e.g. an empty `elev` list).

    Notes
    -----
    Equivalent (commented out above) to a vectorized numpy
    ``argmin(abs(elev - obs_elev))`` with a step-back correction, but
    implemented as an explicit scan since numpy is not available in
    this Jython environment.
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
    # if a closest match was found, check whether it overshoots obs_elev and step back if so
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
    Convert a HEC-style balance period string into a number of hours.

    Parameters
    ----------
    balance_period : str
        HEC-style interval string containing one of the substrings
        ``'hour'``, ``'day'``, or ``'min'`` (case-insensitive), e.g.
        ``'1Hour'``, ``'1Day'``, ``'15Min'``.

    Returns
    -------
    float
        The equivalent period length in hours. Returns ``None``
        implicitly if none of the recognized substrings are present.

    Notes
    -----
    The numeric portion is extracted by stripping the recognized unit
    substring and converting the remainder to a float, so the input
    must be of the form ``'<number><unit>'``.
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
    Validate that every DSS record path mentions the expected balance
    period, aborting the compute if any do not.

    Parameters
    ----------
    records : list of str
        DSS record paths to check.
    balance_period : str
        The expected period string (e.g. ``'1Day'``) that should
        appear somewhere in each record path (typically the E-part).
    currentAlt : object
        The alternative object being computed. Used to log an error
        message if a mismatch is found.

    Returns
    -------
    None
        Exits the process (``sys.exit(-1)``) if any record does not
        contain `balance_period`.
    """
    # check every record for the expected time interval substring
    for r in records:
        if balance_period.lower() not in r.lower():
            currentAlt.addComputeMessage('DSS record {0} not matching time interval {1}'.format(r, balance_period))
            sys.exit(-1)                      # Terminate compute if intervals mismatch


def read_ts_rec_w_optional_fname(dssFm, pathname, starttime_str, endtime_str):
    """
    Read a DSS time series record, optionally from an alternate DSS
    file specified inline in the pathname string.

    Parameters
    ----------
    dssFm : object
        An already-open HEC-DSS file handle to use if `pathname` does
        not specify an alternate file.
    pathname : str
        DSS record path. May optionally contain the DSS filepath
        before the DSS ts path, separated by ``'::'`` (e.g.
        ``'C:/data/alt.dss::/A/B/C/D/E/F/'``). If present, that file
        is opened, read, and closed instead of using `dssFm`.
    starttime_str : str
        Start time string for the windowed read.
    endtime_str : str
        End time string for the windowed read.

    Returns
    -------
    object
        The `TimeSeriesContainer` data for the requested record and
        time window.
    """
    # if an alternate DSS file is embedded in the path string, open and read from that file instead
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
    Read and sum multiple inflow and outflow DSS records over a given
    time window, converting units to cfs as needed, and return the
    net inflow-minus-outflow time series.

    Parameters
    ----------
    currentAlt : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging.
    dss_file : str
        Path to the primary DSS file to read records from (unless a
        record path specifies an alternate file via ``'::'``).
    inflow_records : list of str
        DSS pathnames of records to sum together as total inflow.
    outflow_records : list of str
        DSS pathnames of records to sum together as total outflow.
    starttime_str : str
        Start time string for the windowed reads.
    endtime_str : str
        End time string for the windowed reads.
    starttime_hectime : int
        HEC integer time corresponding to `starttime_str`, used for
        trimming records that start earlier than the window.
    endtime_hectime : int
        HEC integer time corresponding to `endtime_str`, used for
        trimming records that extend later than the window.

    Returns
    -------
    tuple of (list, list)
        ``(times, inflow_outflow)`` where `times` are the HEC times
        of the first inflow record read, and `inflow_outflow` is the
        elementwise net flow (inflow minus outflow) in cfs, aligned
        to those times (offset by one, see Notes).

    Raises
    ------
    SystemExit
        Calls ``sys.exit(-1)`` if a `HecMathException` occurs while
        reading any of the input records.

    Notes
    -----
    - The returned `inflow_outflow` list is built from
      ``inflows[1:]`` and ``outflows[1:]``, i.e. it skips the first
      timestep of each accumulated series.
    - If a record's units are ``'cms'``, values are converted to cfs
      using the factor ``35.314666213``.
    """
    dssFm = HecDss.open(dss_file)

    inflows = []
    outflows = []
    times = []                                # Shared time vector from the first inflow series

    # --- Read and sum all inflow records ---
    # Read inflows
    print('Reading inflows')
    # for each inflow record, read it, trim it to the time window, convert units if
    # needed, and add it to the running inflow total
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
            # if the record's data extends past the run window, trim the trailing values
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)                                 # Fatal read error

        # if this record is in cms, convert it to cfs before summing
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
            sys.exit(-1)                                  # Fatal read error

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
    Predict elevation and storage over a period from inflow/outflow
    balance.

    Useful for generating lookback/starting elevations for forecast
    runs that begin mid-period; integrates inflow-outflow into
    storage and maps storage back to elevation via interpolation.

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
    - Elevation mapping uses linear interpolation on storage to elevation curve.
    - Conic interpolation and evaporation are not yet supported here
      (see TODO comments in the function body).
    """

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

    # Integrate inflow-outflow over the period and map storage to elevation
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
    Compute a reservoir mass-balance ("balance flow") time series and
    write it to DSS.

    The balance flow is the residual flow needed to reconcile
    observed storage change with recorded inflows, outflows, and
    evaporation, with optional companion evaporation and storage
    records and an optional resampled copy at an alternate time
    interval.

    Parameters
    ----------
    currentAlt : object
        The alternative object being computed. Must support
        `addComputeMessage()` for logging.
    timewindow : object
        The run time window object, used to get the start and end
        time strings for reading input data.
    res_name : str
        Reservoir name, used to label the debug CSV output file.
    inflow_records : list of str
        DSS paths of records to sum as total inflow.
    outflow_records : list of str
        DSS paths of records to sum as total outflow.
    stage_record : str
        DSS path of the reservoir stage (elevation) record.
    evap_record : str
        DSS path of the evaporation record.
    elev_stor_area : dict
        Elevation-storage-area lookup table (as produced by
        `read_elev_storage_area_file`), with keys `'elev'`, `'stor'`,
        and `'area'`.
    dss_file : str
        Path to the DSS file containing the input records.
    output_dss_record_name : str
        DSS path to write the computed balance flow record to.
    output_dss_file : str
        Path to the DSS file to write output records to.
    shared_dir : str
        Directory to write the debug CSV file to.
    storage_dss_record_name : str, optional
        DSS path to write the computed storage record to, if
        `write_storage` is True. Defaults to ``''``.
    evap_dss_record_name : str, optional
        DSS path to write the computed evaporation-flow-loss record
        to, if `write_evap` is True. Defaults to ``''``.
    balance_period_str : str, optional
        The expected time interval of all input records (e.g.
        ``'1HOUR'``). Defaults to ``"1HOUR"``.
    use_conic : bool, optional
        If True, use conic frustum interpolation (via
        `get_elev_layer_idx` and `conic_storage_interp`) to estimate
        storage from stage. If False, use simple linear
        interpolation. Defaults to False.
    write_evap : bool, optional
        If True, also write the computed evaporation-flow-loss series
        to DSS. Defaults to False.
    write_storage : bool, optional
        If True, also write the computed storage series to DSS.
        Defaults to False.
    alt_period : optional
        If not None, triggers writing an additional resampled copy of
        the balance flow at `alt_period_string`'s interval.
    alt_period_string : str, optional
        The alternate time interval string to resample to, used only
        if `alt_period` is not None.
    lookback_padding : int, optional
        Currently unused within the function body (referenced only in
        a commented-out block). Defaults to 1440.

    Returns
    -------
    bool
        True once the balance flow (and any requested companion
        records) have been written to DSS.

    Notes
    -----
    - Edge-value safeguard: the first ``check_steps`` and last
      ``check_steps`` computed residual-flow values are zeroed out if
      NaN or unrealistically large in magnitude (beyond
      ``bad_flow_bound``), since timezone/period-boundary effects on
      the first/last days can otherwise produce large spurious
      values. `check_steps` is 24 for an hourly balance period, else 1.
    - The written balance-flow series duplicates its first computed
      value twice at the front (``[flow_resid[0], flow_resid[0]] +
      flow_resid``) as a workaround so that later time-averaging
      doesn't drop the first needed timestep for a ResSim run.
    """
    # --- Validate that all input records share the expected time interval ---
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

    # --- Read and sum inflow/outflow, then read stage and evaporation ---
    times,inflow_outflow = read_inflows_outflows(currentAlt, dss_file, inflow_records, outflow_records, 
                                                 starttime_str, endtime_str, starttime_hectime, endtime_hectime)
    print('len times:',len(times))                          # Diagnostic
    print('len inflow_outflow:',len(inflow_outflow))        # Diagnostic

    dssFm = HecDss.open(dss_file)                           # DSS for stage/evap reads

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

    # --- Compute the per-timestep mass balance residual (the balance flow) ---
    # Calculations
    n = len(stage) - 1                                      # Number of intervals
    flow_resid = []                                         # Residual flow (cfs)
    flow_evap = []                                          # Evaporation flow (cfs)
    storage_record = []                                     # Optional storage record (start of each interval)

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

        # Average area across endpoints for evap flow calculation
        area_avg = 0.5 * (linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_start) +
                          linear_interpolation(elev_stor_area['elev'], elev_stor_area['area'], stage_end))
        evap_flow_loss = (evap[k] * area_avg) * acreft_2_cfs  # in cfs

        # Residual required to balance storage change with net inflow minus evap
        resid = delta_stor_flow - (inflow_minus_outflow - evap_flow_loss)
        flow_resid.append(resid)
        flow_evap.append(evap_flow_loss)
        storage_record.append(storage_start)                  # Record storage at interval start


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

    # --- Optionally write companion evaporation and storage records ---
    # if requested, write the computed evaporation-flow-loss series to DSS
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

    # if requested, write the computed storage series to DSS
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
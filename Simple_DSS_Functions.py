"""
Simple_DSS_Functions
=====================

Lightweight DSS helper functions for summing multiple time series together
and resampling a time series to a new interval, used within HEC-WAT / WTMP
scripting workflows.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `hec.heclib.dss`, `hec.io`, `rma.util.RMAConst`,
  `hec.hecmath`, `hec.heclib.util.Heclib`.
"""

from hec.heclib.dss import HecDss                     # HEC-DSS: open/read/write DSS files
from hec.io import DSSIdentifier                      # DSS path identifier helper (imported for completeness; not directly used)
from hec.io import TimeSeriesContainer                # Container class for writing summed time series to DSS
from rma.util.RMAConst import MISSING_DOUBLE          # RMA sentinel for missing values (imported for completeness; not directly used)
from hec.hecmath import HecMathException              # HEC math exception type (imported for completeness; not directly used)
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE   # HEC sentinel for undefined values (imported for completeness; not directly used)
import hec.hecmath.TimeSeriesMath as tsmath           # Time series math wrapper (imported for completeness; not directly used)

def add_DSS_Data(currentAlt, dssFile, timewindow, input_data, output_path):
    """
    Read multiple DSS time series, sum their values together, and
    write the combined result out as a single new time series.

    For each DSS path in `input_data`, this function reads the time
    series over the given time window and adds its values
    element-wise into a running total. The final summed series is
    written to `output_path` in the same DSS file.

    Parameters
    ----------
    currentAlt : object
        The current alternative object, used for logging compute
        messages via `addComputeMessage`.
    dssFile : str
        Path to the DSS file to read from and write to.
    timewindow : object
        Object providing `getStartTimeString()`,
        `getEndTimeString()`, `getStartTime()`, and `getEndTime()`
        for the compute time window.
    input_data : list of str
        A list of DSS paths (records) to read and sum together.
    output_path : str
        The DSS path to write the summed output time series to.

    Returns
    -------
    int
        Always returns 0.

    Notes
    -----
    All input records are assumed to share the same timestamps and
    units; only the values from the first-read record are used to
    seed the output container's `units`, and only the last-read
    record's `times` is used for the output container's time vector
    (since `times` and `units` are simply overwritten on each loop
    iteration rather than validated for consistency).
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    output_data = []
    # loop through each input DSS path, read its time series, and sum
    # its values into output_data
    for dsspath in input_data:
        print('reading', str(dsspath))
        ts = dssFm.read(dsspath, starttime_str, endtime_str, False)
        ts = ts.getData()
        values = ts.values
        times = ts.times
        units = ts.units
        # if this is the first time series read, use its values as the
        # starting point for the running total
        if len(output_data) == 0:
            output_data = values
        else:
            # otherwise, add this time series' values into the running
            # total value by value
            for vi, val in enumerate(values):
                output_data[vi] += val
                
    # build a new TimeSeriesContainer to hold the summed output 
    tsc = TimeSeriesContainer()
    tsc.times = times
    # set the DSS output path
    tsc.fullName = output_path
    tsc.values = output_data
    # set up series 
    tsc.startTime = times[0]
    tsc.units = units
    tsc.endTime = times[-1]
    tsc.numberValues = len(output_data)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    # write the completed time series container to the DSS file, then close it
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(output_data)))
    return 0

def resample_dss_ts(inputDSSFile, inputRec, timewindow, outputDSSFile, newPeriod):
    """
    Resample a DSS time series to a new time interval.

    Can upsample an even period DSS timeseries, e.g. go from 1DAY ->
    1HOUR. Reads the input time series from `inputDSSFile` over the
    given time window, transforms it to the new period using an
    averaging method, and writes the resampled series to
    `outputDSSFile`.

    Parameters
    ----------
    inputDSSFile : str
        Path to the DSS file to read the input time series from.
    inputRec : str
        The DSS path (record) of the time series to resample.
    timewindow : object
        Object providing `getStartTimeString()` and
        `getEndTimeString()` for the read time window.
    outputDSSFile : str
        Path to the DSS file to write the resampled time series to.
    newPeriod : str
        The new time interval to resample to, e.g. `"1HOUR"`.

    Returns
    -------
    None
        Writes the resampled series to `outputDSSFile`.
    """
    dssFm = HecDss.open(inputDSSFile)
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    tsm = dssFm.read(inputRec, starttime_str, endtime_str, False)
    tsm_new = tsm.transformTimeSeries(newPeriod,"","AVE")
    dssFm.close()
    # write to dss file
    dssFmout = HecDss.open(outputDSSFile)
    dssFmout.write(tsm_new)
    dssFmout.close()
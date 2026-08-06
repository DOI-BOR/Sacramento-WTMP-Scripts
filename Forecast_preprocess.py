"""
Forecast Data Preprocess

This module orchestrates several preprocessing steps needed before running
WTMP/ResSim/W2 forecast workflows. It includes:

- **Meteorology and equilibrium temperature**: Reading met data, computing
  equilibrium water temperature (`equilibrium_temp`), and writing timeseries at
  multiple intervals.
- **Reservoir elevation derivations**: Converting storage to elevation using
  tabulated elevation-storage-area data, inventing elevation records where
  timing surrogates are needed, and predicting future elevations from inflows/
  outflows.
- **Target temperature routing**: Building upstream temperature targets from
  downstream targets, including river travel-time and atmospheric exchange
  effects.
- **Utility helpers**: Lightweight 1-D linear interpolation, fractional-month
  weights, and path/integration helpers.

This file is designed for Jython within HEC-WAT, relying on Java-backed HEC
libraries and project context objects.
"""

from hec.heclib.dss import HecDss  # HEC-DSS file manager: open/read/write time series records
from hec.hecmath import HecMathException  # Exception class for HEC math operations (e.g., transforms)
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE  # Heclib sentinel for undefined doubles
from hec.heclib.util import HecTime  # HEC time object representing DSS-friendly dates/times
from hec.io import DSSIdentifier  # DSS pathname identifier helper (A-F part parsing/building)
from hec.io import TimeSeriesContainer  # Container for DSS time series (times/values/units/type)
from hec.io import DataContainer  # Base container type used across HEC library data structures
import hec.hecmath.TimeSeriesMath as tsmath  # Time series math wrapper for transforms/regularization
from rma.util.RMAConst import MISSING_DOUBLE  # RMA/HEC sentinel representing missing double values
import math  # Standard math utilities (logarithms, trig, etc.)
import sys  # System-level utilities (path, exit, etc.)
import datetime as dt  # Python datetime for timestamp arithmetic and formatting
import os, sys, csv, calendar  # Filesystem ops, duplicate sys import (preserved), CSV I/O, month-day utilities

from com.rma.io import DssFileManagerImpl  # Java DSS file manager (not directly used but preserved import)
from com.rma.model import Project  # Project context: workspace paths, directories for current HEC-WAT project

# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento", "American", "Stanislaus"]  # Terms indicating study-specific paths

# initialize and search for unwanted paths
matching_paths = []
# for each path in sys.path, check if it contains any unwanted phrase
for p in sys.path:
    # if the path matches one of the unwanted phrases, mark it for removal
    if any(phrase in p for phrase in search_list):
        matching_paths.append(p)

# print paths containing unwanted phrases
print("Paths to be removed:")

for path in matching_paths:
    print(path)  # Log each matching path prior to removal

# remove matching paths from sys.path
for path in matching_paths:
    if path in sys.path:
        sys.path.remove(path)  # Remove to avoid import collisions or stale modules

# append path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))  # Ensure 'scripts' is importable

# --- Local module imports (reload ensured) ------------------------------------

import DSS_Tools  # Local DSS utility functions (read/write helpers, unit handling, path fixes)
reload(DSS_Tools)  # Reload to ensure latest definitions in Jython runtime

import DMS_preprocess  # Preprocessing utilities for DMS output types/units
reload(DMS_preprocess)  # Reload to pick up any changes during development

import equilibrium_temp  # Equilibrium temperature computations (radiation + flux models)
reload(equilibrium_temp)  # Reload to ensure updated constants/solvers

import create_balance_flow_jython as cbfj  # Balance-flow helpers and linear interpolation utilities
reload(cbfj)  # Reload to ensure latest behavior

from com.rma.io import DssFileManagerImpl  # (Duplicate import in original; preserved)
from java.util import TimeZone  # Java TimeZone class (not directly used here but preserved)




def interp(x, xp, fp, left=None, right=None):
    """
    One-dimensional linear interpolation (piecewise).

    Parameters
    ----------
    x : float or list of float
        The x-coordinate(s) at which to evaluate the interpolant.
    xp : 1-D sequence of float
        Monotonic increasing x-coordinates for known data points.
    fp : 1-D sequence of float
        y-coordinates for known data points (same length as `xp`).
    left : float, optional
        Return value for `x < xp[0]` (default: `fp[0]`).
    right : float, optional
        Return value for `x > xp[-1]` (default: `fp[-1]`).

    Returns
    -------
    float or list of float
        Interpolated value(s) at `x`, with edge extrapolation using `left`/`right`.

    Notes
    -----
    - This is a lightweight alternative to `numpy.interp`, implemented to
      keep compatibility with Jython environments lacking NumPy/SciPy.
    - When `x` is a list, interpolation is applied element-wise.
    """

    # if x is a list, interpolate each point recursively and return a list
    if isinstance(x, list):
        # Input type is a list. Recursively interp each of the list entries
        return [interp(point, xp, fp, left, right) for point in x] 
    
    else:
        # Element is not a list. Interpolate the specific values
        if left is None:
            # Default left extrapolation value
            left = fp[0]  
        
        if right is None:
            # Default right extrapolation value
            right = fp[-1]  

        # if x is below the known range, clamp to the left value
        if x < xp[0]:
            return left
        # if x is above the known range, clamp to the right value
        elif x > xp[-1]:
            # Right extrapolation
            return right  
        
        else:
            # find the bracketing interval in xp and linearly interpolate within it
            for i in range(len(xp) - 1):
                
                # Find bracketing interval
                if x >= xp[i] and x <= xp[i + 1]:  
                    # Perform the linear interpolation
                    t = (x - xp[i]) / (xp[i + 1] - xp[i])  # Fraction within interval
                    
                    # Return linear interpolation result
                    return fp[i] + t * (fp[i + 1] - fp[i])  


def eq_temp(rtw, at, cl, ws, sr, td, eq_temp_out):
    """
    Compute equilibrium water temperature series and write multiple
    intervals to DSS.

    Computes equilibrium water temperature over the run time window
    from meteorological inputs, and writes the hourly result plus
    daily and weekly standardized versions to DSS.

    Parameters
    ----------
    rtw : object
        The run time window object, used to get the start and end
        time strings for reading input data.
    at : list
        `[dss_file_path, dss_record_path]` for air temperature.
    cl : list
        `[dss_file_path, dss_record_path]` for cloud cover.
    ws : list
        `[dss_file_path, dss_record_path]` for wind speed.
    sr : list
        `[dss_file_path, dss_record_path]` for solar radiation.
    td : list
        `[dss_file_path, dss_record_path]` for dew point temperature.
    eq_temp_out : list
        `[dss_file_path, dss_record_path]` specifying where to write
        the resulting equilibrium temperature record.

    Returns
    -------
    None
        Writes the hourly equilibrium temperature record, and its
        1-day and 1-week standardized-interval versions, to the DSS
        file specified in `eq_temp_out`.
    """
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()

    # get at data and times in the formats needed
    dssFm = HecDss.open(at[0])        
    tsc = dssFm.read(at[1], starttime_str, endtime_str, False).getData()
    tsc_int_times = tsc.times
    dtt = DSS_Tools.hectime_to_datetime(tsc)
    at_data = tsc.values
    dssFm.close()
    
    # get the rest of the data over the same period
    cl_data = DSS_Tools.data_from_dss(cl[0],cl[1],starttime_str,endtime_str)
    ws_data = DSS_Tools.data_from_dss(ws[0],ws[1],starttime_str,endtime_str)
    sr_data = DSS_Tools.data_from_dss(sr[0],sr[1],starttime_str,endtime_str)
    td_data = DSS_Tools.data_from_dss(td[0],td[1],starttime_str,endtime_str)

    # calc_equilibrium_temp(dtt, at, cl, sr, td, ws)
    Te = equilibrium_temp.calc_equilibrium_temp(dtt,at_data,cl_data,sr_data,td_data,ws_data)
    
    print('writing: ',eq_temp_out[1])
    # build the time series container for the computed equilibrium temperature
    tsc = TimeSeriesContainer()
    tsc.times = tsc_int_times
    tsc.fullName = eq_temp_out[1]
    tsc.values = Te
    tsc.units = 'C'
    tsc.type = 'INST-VAL'
    tsc.numberValues = len(tsc.values)

    # also produce daily and weekly standardized-interval versions of the same series
    tsm = tsmath(tsc)
    tsm_day = DSS_Tools.standardize_interval(tsm,'1day')
    tsm_wk = DSS_Tools.standardize_interval(tsm,'1week')
        
    # write all three versions (hourly, daily, weekly) to the output DSS file
    dssFmOut = HecDss.open(eq_temp_out[0])
    dssFmOut.write(tsc)
    dssFmOut.write(tsm_day)
    dssFmOut.write(tsm_wk)
    dssFmOut.close()


def storage_to_elev(res_name,elev_stor_area,forecast_dss,storage_rec,conic=False):
    """
    Convert a reservoir storage time series into an elevation time
    series.

    Uses a storage-elevation-area lookup table, and writes the result
    back to DSS under the same F-part.

    Parameters
    ----------
    res_name : str
        Reservoir name to place in the B-part of the output DSS path.
    elev_stor_area : dict
        Lookup table with `'stor'` and `'elev'` keys (parallel
        lists), used for interpolation.
    forecast_dss : str
        Path to the DSS file containing the storage record and where
        the elevation record will be written.
    storage_rec : str
        DSS path of the input storage record to convert.
    conic : bool, optional
        If True, use conic interpolation (not yet implemented - the
        function will exit). Defaults to False (linear
        interpolation).

    Returns
    -------
    None
        Writes the converted elevation record back to `forecast_dss`,
        reusing the same F-part as `storage_rec` but with B-part set
        to `res_name` and C-part set to `'ELEV'`.

    Raises
    ------
    SystemExit
        Calls `sys.exit(-1)` if `conic` is True, since conic
        interpolation from storage is not yet supported.
    """
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True) # read ALL data in record

    elev = []
    # if conic interpolation was requested, it is not supported yet, so exit
    if conic:
        # Use conic interpolation
        print('Conic interpolation of elevations from storage not supported yet.')
        sys.exit(-1)  # Preserve original error path
    else:
        # for each stored value, linearly interpolate its corresponding elevation
        for j in range(tsc.numberValues):
            elev.append(cbfj.linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], tsc.values[j]))
            print('stor2elev: ',j,tsc.times[j],tsc.values[j])

    # rewrite the DSS path's B-part (location) and C-part (parameter) for the elevation record
    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = elev
    dssFmRec.write(tsc)
    dssFmRec.close()

def invent_elevation(res_name,forecast_dss,storage_rec,elev_constant_ft):
    """
    Create a constant-value elevation time series and write it to
    DSS.

    Reuses the same timestamps as an existing storage record. This is
    used to fabricate a placeholder elevation record purely for
    timing purposes, when an actual computed elevation is not needed
    or not yet available.

    Parameters
    ----------
    res_name : str
        Reservoir name to place in the B-part of the output DSS path.
    forecast_dss : str
        Path to the DSS file containing the reference storage record
        and where the elevation record will be written.
    storage_rec : str
        DSS path of an existing record whose timestamps will be
        reused for the constant elevation series.
    elev_constant_ft : float
        The constant elevation value, in feet, to assign to every
        timestamp.

    Returns
    -------
    None
        Writes the constant elevation record to `forecast_dss`.
    """
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True)
    
    # rewrite the DSS path's B-part (location) and C-part (parameter) for the elevation record
    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = [elev_constant_ft for j in range(tsc.numberValues)]
    
    # Write and close the series
    dssFmRec.write(tsc)  # Write to DSS
    dssFmRec.close()  # Close file


def write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir):
    """
    Compute and write forecasted elevation time series for the
    5-reservoir Sacramento/Trinity system.

    Computes elevations for Shasta Lake, Trinity Lake, and
    Whiskeytown Lake (plus placeholder elevations for Keswick and
    Lewiston reservoirs), based on storage-elevation-area
    relationships and mass-balance flow routing.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed, passed through to
        `cbfj.predict_elevation` for logging.
    rtw : object
        The run time window object, used to determine the forecast
        start and end dates.
    forecast_dss : str
        Path to the shared forecast DSS file containing storage and
        flow records, and where elevation records will be written.
    shared_dir : str
        Directory containing the storage-elevation-area CSV lookup
        files.

    Returns
    -------
    None
        Writes multiple elevation and forecast records to
        `forecast_dss` for Shasta, Trinity, and Whiskeytown.

    Notes
    -----
    Workflow (repeated per reservoir, with reservoir-specific
    inflow/outflow records):

    1. Determines the starting date for the elevation forecast (the
       end of the prior month, or Dec 31 of the previous year if the
       run starts in January).
    2. Loads the storage-elevation-area lookup table for the
       reservoir from a CSV file.
    3. Converts both the monthly and daily (CVP-forecast) storage
       records into elevation records via `storage_to_elev`.
    4. For Keswick and Lewiston, invents placeholder elevation
       records (constant values) from the Shasta/Trinity storage
       records, since actual elevation isn't computed for these but a
       record is needed for timing.
    5. Resamples relevant flow-release records to daily resolution.
    6. Looks up the starting elevation from the monthly elevation
       record.
    7. Calls `cbfj.predict_elevation` to forecast elevation forward
       using inflow and outflow records and a daily mass balance.
    """
    
    # get date for starting elevation - look for end-of-month before start time
    ht = HecTime(rtw.getStartTimeString())
    # if the run starts in January, the prior month-end is Dec 31 of the previous year
    if ht.month() == 1:
        start_str = dt.datetime(ht.year() - 1, 12, 31).strftime('%d%b%Y') + ' 2400'  # Previous year-end
    else:
        start_dt = dt.datetime(ht.year(),ht.month(),1)
        start_dt = start_dt - dt.timedelta(days=1)
        start_str = start_dt.strftime('%d%b%Y')+ ' 2400'
    end_str = rtw.getEndTimeString()   
    
    # Shasta
    # covert storage to monthly elevation
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_shasta.csv'), 'Shasta')  # Load curve
    storage_to_elev('Shasta Lake', elev_stor_area, forecast_dss, '/SACRAMENTO RIVER/SHASTA LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/', conic=False)
    
    # convert bc script predicted storage
    storage_to_elev('Shasta Lake', elev_stor_area, forecast_dss, '/SACRAMENTO RIVER/SHASTA LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/', conic=False)

    # invent flow-reg reservoir elevation record from shasta storage rec (used for timing only)
    invent_elevation('Keswick Reservoir', forecast_dss, '/SACRAMENTO RIVER/SHASTA LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/', 582.0)
    invent_elevation('Lewiston Reservoir', forecast_dss, '/TRINITY RIVER/TRINITY LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/', 1901.0)

    # also make a one day step, to see if that solves some issues (used for timing only)
    invent_elevation('Keswick Reservoir', forecast_dss, '/SACRAMENTO RIVER/SHASTA LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/', 582.0)
    invent_elevation('Lewiston Reservoir', forecast_dss, '/TRINITY RIVER/TRINITY LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/', 1901.0)

    # write an hourly forecast elevation based on starting elevation and flows
    DSS_Tools.resample_dss_ts(forecast_dss, '/SACRAMENTO RIVER/SHASTA LAKE/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/', None, forecast_dss, '1DAY')

    inflow_records = ['//SHASTA-PIT-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SHASTA-SAC-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SHASTA-SULANHARAS-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SHASTA-MCCLOUD-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '/SACRAMENTO RIVER/SHASTA LAKE/FLOW-ACC-DEP//1DAY/SACTRN_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['/SACRAMENTO RIVER/SHASTA LAKE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)
    # forecast Shasta elevation forward using the daily inflow/outflow mass balance
    cbfj.predict_elevation(currentAlternative, start_str,end_str, 'Shasta Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Shasta Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')  # Predict daily elevation

    # Trinity
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_trinity.csv'), 'trinity')  # Load curve
    storage_to_elev('Trinity Lake', elev_stor_area, forecast_dss, '/TRINITY RIVER/TRINITY LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/', conic=False)
    
    # convert bc script predicted storage
    storage_to_elev('Trinity Lake', elev_stor_area, forecast_dss, '/TRINITY RIVER/TRINITY LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/', conic=False)

    DSS_Tools.resample_dss_ts(forecast_dss, '/TRINITY RIVER/TRINITY LAKE/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/', None, forecast_dss, '1DAY')
 
    inflow_records = ['//EF TRINITY/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//STUART FORK/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SWIFT CR/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//TRINITY RIVER/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '/TRINITY RIVER/TRINITY LAKE/FLOW-ACC-DEP//1DAY/SACTRN_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['/TRINITY RIVER/TRINITY LAKE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'/TRINITY RIVER/TRINITY LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)
    # forecast Trinity elevation forward using the daily inflow/outflow mass balance
    cbfj.predict_elevation(currentAlternative, start_str,end_str, 'Trinity Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Trinity Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')  # Predict daily elevation

    # Whiskeytown
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_whiskeytown.csv'), 'Whiskeytown')  # Load curve
    storage_to_elev('Whiskeytown Lake', elev_stor_area, forecast_dss, '/CLEAR CREEK/WHISKEYTOWN LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/', conic=False)

    # convert bc script predicted storage
    storage_to_elev('Whiskeytown Lake', elev_stor_area, forecast_dss, '/CLEAR CREEK/WHISKEYTOWN LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/', conic=False)

    DSS_Tools.resample_dss_ts(forecast_dss, '/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1HOUR/SACTRN_BC_SCRIPT/', None, forecast_dss, '1DAY')
    DSS_Tools.resample_dss_ts(forecast_dss, '/CLEAR CREEK/WHISKEYTOWN DAM/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/', None, forecast_dss, '1DAY')
    DSS_Tools.resample_dss_ts(forecast_dss, '/CLEAR CREEK/CARR POWERHOUSE/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/', None, forecast_dss, '1DAY')

    inflow_records = ['/USBR-LINEARINTERP/CLEAR CR ABOVE JCR INFLOW/FLOW//1DAY/SACTRN_BC_SCRIPT/',
                      '/CLEAR CREEK/CARR POWERHOUSE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/',
                      '/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-ACC-DEP//1DAY/SACTRN_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow

    outflow_records = ['/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1DAY/SACTRN_BC_SCRIPT/',
                       '/CLEAR CREEK/WHISKEYTOWN DAM/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'/CLEAR CREEK/WHISKEYTOWN LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)
    # forecast Whiskeytown elevation forward using the daily inflow/outflow mass balance
    cbfj.predict_elevation(currentAlternative, start_str,end_str, 'Whiskeytown Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Whiskeytown Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')  # Predict daily elevation


def splice_met(currentAlternative, rtw, forecast_dss, output_dss):
    """
    Splice Redding (KRDD) meteorological data into the Lewiston
    Reservoir meteorological records for January-March.

    Lewiston is still dependent on Redding data during those months.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed, passed through to
        `DSS_Tools.replace_data` for logging.
    rtw : object
        The run time window object, defining the period over which
        data may be replaced.
    forecast_dss : str
        Path to the DSS file containing the source (Redding) and
        target (Lewiston) records.
    output_dss : str
        Path to the DSS file where the spliced records should be
        written.

    Returns
    -------
    None
        Delegates the actual replacement to `DSS_Tools.replace_data`
        for each Lewiston/Redding record pair, limited to the
        specified months.
    """
    # Lewiston is still dependent on using Met data from Redding during Jan-Feb-Mar.  Create those spliced Met data records...
    pairs = [
            ["/MR Sac.-Lewiston Res./TCAC1/Temp-Air//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Lewiston Res./TCAC1/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR SAC.-LEWISTON RES./TCAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/",
             "/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            #["/MR Sac.-Lewiston Res./TCAC1/Dir-Wind-radians//1Hour/SACTRN_BC_SCRIPT/", # removing this, b/c this is not used and the unit diff is causing forecast problems
            #"/MR Sac.-Clear Cr. to Sac R./KRDD/Dir-Wind//1Hour/SACTRN_BC_SCRIPT/"],  # already in radians
            ["/MR Sac.-Lewiston Res./TCAC1/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Trinity River/TCAC1/%-Cloud Cover//1Day/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Trinity River/TCAC1/%-Cloud Cover-FRAC//1Day/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"]
             ]
    
    # Define the months over which data is being replaced with pair #2
    months = [1, 2, 3]
    
    # Execute the replacement
    DSS_Tools.replace_data(currentAlternative, rtw, pairs, forecast_dss, output_dss, months, standard_interval='1HOUR')


def study_dir_from_run_dir(run_dir):
    """
    Derive the top-level study directory from a compute run
    directory.

    Walks up three directory levels from `run_dir`.

    Parameters
    ----------
    run_dir : str
        The full path to the current run's directory.

    Returns
    -------
    str
        The path to the study directory (three levels above
        `run_dir`).
    """
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def model_dir_from_run_dir(run_dir,model_place,model_name):
    """
    Build the path to a specific CE-QUAL-W2 model directory.

    Parameters
    ----------
    run_dir : str
        The full path to the current run's directory, used to locate
        the parent study directory.
    model_place : str
        Subdirectory under `'cequal-w2'` where the model lives.
    model_name : str
        Name of the specific model directory.

    Returns
    -------
    str
        The full path to the requested W2 model directory.
    """
    study_dir = study_dir_from_run_dir(run_dir)
    model_dir = os.path.join(study_dir,'cequal-w2',model_place,model_name)
    return model_dir

def forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions):
    """
    Preprocess forecast boundary condition data for the ResSim
    5-Reservoir Sacramento/Trinity system.

    Fixes data types/units, applies temperature lapse and met-data
    splicing, forecasts reservoir elevations, computes equilibrium
    temperature, creates a set of constant reference DSS records, and
    computes an upstream (Shasta) target temperature from the
    downstream target.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()` and `getTimeStep()` for logging and
        retrieving the balance period.
    computeOptions : object
        The compute options/settings object. Must support
        `getDssFilename()`, `getRunTimeWindow()`, and
        `getRunDirectory()`.

    Returns
    -------
    bool
        True once all preprocessing steps have completed
        successfully.

    Notes
    -----
    Workflow:

    1. Retrieves the DSS filename, run time window, run directory,
       and shared data directory for this compute.
    2. Fixes data types/units in the shared forecast DSS file via
       `DMS_preprocess.fix_DMS_types_units`.
    3. Applies an air temperature lapse rate correction for the
       elevation at Shasta Lake.
    4. Splices in Redding meteorological data for Lewiston during
       Jan-Mar via `splice_met`.
    5. Forecasts elevations for Shasta, Trinity, and Whiskeytown (and
       placeholders for Keswick/Lewiston) via
       `write_forecast_elevations`.
    6. Computes equilibrium water temperature via `eq_temp`.
    7. Creates a set of constant-value reference DSS records (tiny
       flow, zero flow/gate, constant temperature targets, etc.) used
       elsewhere in the model as fixed boundary conditions.
    8. Computes relative humidity from air temperature and dew point.
    9. Resamples the Keswick flow-release and target temperature
       records to daily resolution.
    10. Reads the configured downstream control location.
    11. If the downstream location is at Shasta Dam itself
        (`location == 0`), copies the target temperature directly.
        Otherwise, back-calculates the required upstream (Shasta)
        target temperature from the downstream target via
        `upstream_target`.
    """
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    forecast_dss = os.path.join(shared_dir,'WTMP_SacTrn_Forecast.dss')
    DMS_preprocess.fix_DMS_types_units(forecast_dss)
    
    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+forecast_dss)
    currentAlternative.addComputeMessage('lapse outfile: '+forecast_dss)
    DSS_Tools.airtemp_lapse(forecast_dss, "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-AIR//1HOUR/SACTRN_BC_SCRIPT/",
                  0.7, forecast_dss, "Shasta_Lapse")
    
    # splice in Redding met data for Lewiston during Jan-Mar
    splice_met(currentAlternative, rtw, forecast_dss, forecast_dss)

    # forecast reservoir elevations for Shasta, Trinity, and Whiskeytown
    write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")
    # eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out)
    eq_temp(rtw,
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Hour/sactrn_bc_script/"]
           )

    #  Create a set of fixed constant-value reference records 
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')   
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=1, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ONES',fpart='ONES')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=13.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='WHI-target-13',fpart='WHI-target-13')

    # Compute/Write EQ temperature
    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")  # Notice
    
    eq_temp(rtw,
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Hour/sactrn_bc_script/"]
           )  

    # Create the records
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow',
                        dss_type='PER-AVER', period='1DAY', cpart='TinyFlow', fpart='TinyFlow')  # Tiny daily flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow',
                        dss_type='PER-AVER', period='1HOUR', cpart='TinyFlow', fpart='TinyFlow')  # Tiny hourly flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water',
                        dss_type='PER-AVER', period='1DAY', cpart='TENS', fpart='TENS')  # Daily temp 10C
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow',
                        dss_type='PER-AVER', period='1DAY', cpart='ZEROS', fpart='ZEROS')  # Daily zero flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow',
                        dss_type='PER-AVER', period='1HOUR', cpart='ZEROS', fpart='ZEROS')  # Hourly zero flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0, what='gate',
                        dss_type='INST-VAL', period='1HOUR', cpart='ZEROS', fpart='ZEROS')  # Hourly gate zero
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=1, what='gate',
                        dss_type='INST-VAL', period='1HOUR', cpart='ONES', fpart='ONES')  # Hourly gate ones
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=13.0, what='temp-water',
                        dss_type='PER-AVER', period='1DAY', cpart='WHI-target-13', fpart='WHI-target-13')  # WHI target

    # Estimate the relative humidity
    DSS_Tools.relhum_from_at_dp(forecast_dss,
                      "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-AIR//1HOUR/sactrn_bc_script/",
                      "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-DEWPOINT//1HOUR/sactrn_bc_script/")  # RH derived from AT/DP

    # TODO: Perhaps generate tributary flows/temps based on exceedence and/or temp regressions?

    # compute W2 regression downstream target temps, in case we want to try/use them in ResSim
    # Keswick need daily record.
    DSS_Tools.resample_dss_ts(forecast_dss, '//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY', pad_start_days=1)  # Hour to Day
    DSS_Tools.resample_dss_ts(forecast_dss, '/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY')  # Hour to Day targets

    # read location from DSS
    location = get_downstream_loc(forecast_dss)  # 0=Dam, others downstream points

    # Set the target paths
    TT_rec = "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/"  # Downstream TT
    TT_W2_rec = "/USBR/SHASTA/TEMP-WATER-TARGET-W2-UPSTREAM//1Day/SACTRN_BC_SCRIPT/"  # Upstream TT out
    
    # Apply the regression to account for warming between the release and the 
    if location == 0:
        # At Shasta Dam, use exact TT without adjustment
        DSS_Tools.copy_dss_ts(TT_rec, new_dss_rec=TT_W2_rec, dss_file_path=forecast_dss, checkMakeCelsius=True)  # Copy/convert to C
    
    else:
        # otherwise, back-calculate the upstream (Shasta) target temperature needed
        # to meet the downstream target, accounting for travel time and heating
        upstream_target(forecast_dss,rtw,
                        "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/",
                        "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Day/sactrn_bc_script/",
                        "//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Day/SACTRN_BC_SCRIPT/",
                        "/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1Day/SACTRN_BC_SCRIPT/",
                        location, TT_W2_rec, ResSimRiver=False)  # Build upstream TT via regression/back-routing

    # Return that the setup is successful
    return True  

def route_downstream(tKeswick,keswickFlowDaily,eqTempDaily,step,hour_of_day,loc):
    """
    Route a starting temperature downstream hour-by-hour.

    Allows it to approach the equilibrium temperature at a fixed
    hourly exchange rate, for the travel time appropriate to the
    given location and flow.

    Parameters
    ----------
    tKeswick : float
        Starting temperature (typically at Keswick) to route
        downstream.
    keswickFlowDaily : list of float
        Daily flow values, used to determine travel time via
        `travel_time_hrs`.
    eqTempDaily : list of float
        Daily equilibrium temperature values, which the routed
        temperature drifts toward each hour.
    step : int
        Index into `keswickFlowDaily`/`eqTempDaily` representing the
        starting day.
    hour_of_day : int
        Starting hour of day (0-23) for the routing simulation.
    loc : int
        Downstream location index, passed to `travel_time_hrs` to
        determine travel distance/time.

    Returns
    -------
    float
        The temperature after being routed downstream for the
        computed travel time.
    """
    exchCoef = 0.015 # hourly exchange rate between atmosphere and river temp

    hrs = travel_time_hrs(loc,keswickFlowDaily[step]) # loc 2 = CCR
    t = tKeswick
    i = step # + keswickResidenceTime ?
    imax = len(eqTempDaily)
    h = hour_of_day
    # step forward hour by hour for the computed travel time, letting temperature
    # drift toward the equilibrium temperature at the exchange rate
    for k in range(hrs):
        deltaTemp = (eqTempDaily[i] - t) * exchCoef
        t += deltaTemp
        # advance to the next day's equilibrium temperature index once past hour 23
        if h == 23:
            h = 0  # Wrap hour
            i = min(imax - 1, i + 1)  # Advance day with bound

        else:
            h += 1  # Next hour

    return t  # Routed downstream temperature


def travel_time_hrs(loc,keswickFlow):
    """
    Estimate integer travel time (hours) from Keswick to a
    downstream location.

    Uses a power-law approximation of river velocity based on flow.

    Parameters
    ----------
    loc : int
        Downstream location index. Must be one of: 1 (Highway 44),
        2 (CCR), or 3 (Ball's Ferry).
    keswickFlow : float
        Flow rate at Keswick, in cfs, used to estimate velocity.

    Returns
    -------
    int
        Estimated travel time in whole hours. Returns 0 if the
        computed velocity is not positive (with a warning printed).

    Raises
    ------
    NotImplementedError
        If `loc` is not one of the recognized location indices (1, 2,
        or 3).
    """
    # determine the downstream distance in feet based on the requested location
    if loc == 1:  # Highway 44
        downstreamDistance = 30000.  # in feet
    
    elif loc == 2:  # CCR
        downstreamDistance = 53000.
    
    elif loc == 3:  # Ball's Ferry
        downstreamDistance = 137000.
    
    else:
        raise NotImplementedError('Downstream location index ' + str(loc) + ' not recognized.')

    # Power law approximation for velocity in the Sacramento River
    Kcoef = 2.  # Coefficient for velocity calculation
    alpha = 0.33  # Exponent in power-law relation
    velocity = Kcoef * (keswickFlow / 1000) ** alpha  # power law approximation
    
    # Calculate travel time in seconds
    if velocity <= 0.0:
        print('WARNING: Keswick flow <= 0 in Forecast TCD script. Specified flows may be incorrect.')
        travTime = 0.0  # Zero travel time in degenerate case
    
    else:
        travTime = downstreamDistance / velocity  # Seconds to reach location

    # Return hrs
    return int(travTime / (60.0 * 60.0))  # Convert seconds to whole hours


def fractional_month(date_obj):
    """
    Compute month-position blending weights for a given date.

    Parameters
    ----------
    date_obj : datetime.datetime
        The date to evaluate.

    Returns
    -------
    tuple of float
        A 4-tuple:

        - `fractional_month` : Position of the day within the month,
          from 0.0 (first day) to just under 1.0 (last day).
        - `fractional_previous` : Blending weight toward the previous
          month. Non-zero only when `date_obj` falls in the first
          half of the month.
        - `fractional_current` : Blending weight toward the current
          month. Always the largest of the three weights, peaking
          mid-month.
        - `fractional_next` : Blending weight toward the next month.
          Non-zero only when `date_obj` falls in the second half of
          the month.
    """

    year = date_obj.year  # Year of date
    month = date_obj.month  # Month number
    day = date_obj.day  # Day of month

    # Get the number of days in the current month using the calendar module
    _, days_in_month = calendar.monthrange(year, month)  # (weekday_of_first, days)

    # Calculate the fractional month. Use 1.0 to force float division.
    fractional_month = (day - 1.0) / days_in_month
    fractional_previous = 0.0
    fractional_next = 0.0
    
    # determine whether the date falls in the first or second half of the month,
    # which controls whether it blends toward the previous or next month
    if fractional_month > 0.5:
        fractional_current = (1.0 - fractional_month) + 0.5  # Late month: more current+next
        fractional_next = 1.0 - fractional_current  # Remaining weight goes to next

    else:
        fractional_current = fractional_month + 0.5  # Early month: more previous+current
        fractional_previous = 1.0 - fractional_current  # Remaining weight goes to previous

    # Return to the calling function
    return fractional_month, fractional_previous, fractional_current, fractional_next  # Weights tuple


def get_step_future_and_RiverHrs(wqTargetDaily,keswickFlowDaily,step,loc):
    """
    Determine the future daily time step corresponding to water
    released today.

    Accounts for downstream river travel time and Keswick reservoir
    residence time.

    Parameters
    ----------
    wqTargetDaily : list
        Daily water quality target series, used only for its length
        to clamp the future step index.
    keswickFlowDaily : list of float
        Daily flow values at Keswick, used to estimate travel time.
    step : int
        Index of the current daily time step.
    loc : int
        Downstream location index, passed to `travel_time_hrs`.

    Returns
    -------
    tuple of (int, int)
        `(step_future, hrs)` where `step_future` is the clamped
        future daily index accounting for travel and residence time,
        and `hrs` is the estimated river travel time in hours.
    """
    # Power law approximation for velocity in the Sacramento River
    hrs = travel_time_hrs(loc, keswickFlowDaily[step])  # Compute travel hours

    # Get Keswick pool information - from forecast TCD script
    flowVol = keswickFlowDaily[step] * 86400.  # Daily volume [ft3]
    kesConPoolVol = 20100. * 43560.  # cubic feet, assumed this is top of conservation
    kesFraction = flowVol / kesConPoolVol  # Residence proxy
    multiplier = 0.14  # Calibration factor
    kesFraction = min(kesFraction * multiplier, 1.)  # Bound multiplier effect

    travel_days = int(round(hrs / 24 + kesFraction))  # Combined travel days + residence effect
    step_future = min(step + travel_days, len(wqTargetDaily))  # Bound index within available targets

    # Return future index and travel hours
    return step_future, hrs  


#######################################################################################################
# Backcalculate the temperature required at Shasta Dam from the downstream temperature target
def backRouteWQTarget2(eqTempDaily, targetTempFuture, sha2kes_diff, hrs, step, loc, hour_of_day=10):
    """
    Back-calculate the outlet (Shasta Dam) temperature required to
    meet a downstream target.

    Uses an iterative bisection-style search so that, after Keswick
    heating/cooling and downstream travel time, the resulting
    temperature meets a given future downstream target.

    Parameters
    ----------
    eqTempDaily : list of float
        Daily equilibrium temperature values used during downstream
        routing.
    targetTempFuture : float
        The desired downstream temperature target at the future time
        step.
    sha2kes_diff : float
        Temperature difference between Shasta outlet and Keswick
        (negative if heating occurs in Keswick reservoir).
    hrs : int
        Downstream travel time in hours, used for routing.
    step : int
        Current daily time step index, used as the starting point for
        downstream routing.
    loc : int
        Downstream location index (used only for
        informational/logging purposes here; routing itself uses
        `hrs`).
    hour_of_day : int, optional
        Starting hour of day for the routing simulation. Defaults to
        10.

    Returns
    -------
    float
        The estimated required outlet temperature at Shasta Dam.

    Raises
    ------
    ValueError
        If the search range does not bracket the target temperature
        and the target also cannot be met at the lowest tested outlet
        temperature.
    """
    #exchCoef = 0.0098  # exchange rate between atmosphere and river temp
    exchCoef = 0.015 # exchange rate between atmosphere and river temp
    tSearchMin = targetTempFuture + sha2kes_diff - 6.
    tSearchMax = tSearchMin + 18.
    numIters = 51
    bracketed = False
    cantBeMet = False
    #network.printMessage('Keswick vars ' + str(keswickResAvgTemp) + ', ' + str(kesFraction))
    #network.printMessage('Travel time steps ' + str(travTimeSteps))
    
    # sweep candidate outlet temperatures across the search range, routing each
    # downstream, to find where the routed result brackets the target temperature
    for j in range(numIters):
        # Candidate outlet temp
        outletTemp = tSearchMin + float(j) / float(numIters + 1) * (tSearchMax - tSearchMin)  
        
        # Impact of Keswick
        t = outletTemp - sha2kes_diff  # sha2kes_diff is negative if heating in Keswick
        t1 = 1  # Placeholder variable retained from original
        
        # Route downstream
        i = step 
        imax = len(eqTempDaily)
        h = hour_of_day
        # step forward hour by hour for the travel time, letting temperature drift
        # toward the equilibrium temperature at the exchange rate
        for k in range(hrs):
            deltaTemp = (eqTempDaily[i] - t) * exchCoef  # Relax toward equilibrium
            t += deltaTemp  # Update routed temp

            if h == 23:
                h = 0  # Wrap hour
                i = min(imax - 1, i + 1)  # Advance day

            else:
                h += 1  # Next hour

        print('--', j, t1, t)  # Diagnostic print

        # Check for hte convergence
        if j == 0:
            prevT = t
            prevOutletT = outletTemp
        # check whether the routed temperature crossed the target between this
        # iteration and the previous one, in either direction, to bracket the answer
        if t > targetTempFuture and prevT < targetTempFuture:
            upperOutletT = outletTemp
            upperT = t
            lowerOutletT = prevOutletT
            lowerT = prevT
            bracketed = True
            print('Break loop 1 ' + str(prevT) + ', ' + str(t) + ', ' + str(targetTempFuture))
            break

        elif prevT > targetTempFuture and t < targetTempFuture:
            lowerOutletT = outletTemp
            lowerT = t
            upperOutletT = prevOutletT
            upperT = prevT
            bracketed = True
            print('Break loop 2 ' + str(prevT) + ', ' + str(t) + ', ' + str(targetTempFuture))
            break

        elif j == 0 and t > targetTempFuture:
            cantBeMet = True  # Immediately infeasible with first candidate
            break
        
        prevT = t  # Update previous
        prevOutletT = outletTemp  # Update previous outlet temp

    # resolve the final outlet temperature based on how the search resolved
    if bracketed:
        # Linear interpolation
        targetTemp = (targetTempFuture - lowerT) / (upperT - lowerT) * (upperOutletT - lowerOutletT) + lowerOutletT  # Backrouted TT

    elif cantBeMet:
        targetTemp = outletTemp  # Use last candidate as fallback

    else:
        if t < targetTempFuture:
            targetTemp = outletTemp  # Use last candidate below target

        else:
            # network.printMessage('Target Temperature Downstream' + str(targetTempFuture))
            raise ValueError('Outlet temperature not bracketed')  # No safe interpolation

def upstream_target(forecastDSS,rtw,downstreamTT_rec,eqTemp_rec,kesFlow_rec,sppFlow_rec,loc,TT_W2_rec,ResSimRiver=True):
    """
    Back-calculate the required Shasta Dam (upstream) target
    temperature from a downstream temperature target.

    Uses location-specific regression coefficients and (for
    ResSim-based routing) hour-by-hour downstream routing with
    equilibrium temperature.

    Parameters
    ----------
    forecastDSS : str
        Path to the DSS file containing the downstream target,
        equilibrium temperature, and flow records.
    rtw : object
        The run time window object, used to get the start and end
        time strings for reading input data.
    downstreamTT_rec : str
        DSS path of the downstream target temperature record.
    eqTemp_rec : str
        DSS path of the equilibrium temperature record.
    kesFlow_rec : str
        DSS path of the Keswick flow record.
    sppFlow_rec : str
        DSS path of the Spring Creek Powerplant (or similar) flow
        record.
    loc : int
        Downstream location index. Currently only `loc == 2` is
        supported.
    TT_W2_rec : str
        DSS path where the resulting upstream (Shasta) target
        temperature record should be written.
    ResSimRiver : bool, optional
        If True, use the ResSim-river regression coefficients and
        full hour-by-hour back-routing via `backRouteWQTarget2`. If
        False, use the CCR regression coefficients and a simpler
        direct offset. Defaults to True.

    Returns
    -------
    None
        Writes the computed upstream target temperature record to
        `forecastDSS` at `TT_W2_rec`.

    Raises
    ------
    ValueError
        If `loc` is not 2 (the only currently supported location).

    Notes
    -----
    **Known issue:** `TTDiff` is referenced inside the per-timestep
    loop (in both the `ResSimRiver` and non-`ResSimRiver` branches)
    but is never assigned anywhere in this function - it appears the
    monthly-blended coefficients (`a0, b0, c0, d0`, `a1, b1, c1, d1`,
    `a2, b2, c2, d2`, and the fractional weights `pFrac`, `cFrac`,
    `nFrac`) were intended to be combined to compute `TTDiff` before
    use, but that computation is missing from the source. As written,
    this function will raise a `NameError` when it reaches the first
    call to `backRouteWQTarget2` (or the non-`ResSimRiver` offset
    calculation).
    """
    kes2sha_coeffs = [
        [ 0.5026658277531562 , -0.023436069646011252 , -0.06557997671807375 , -4.685462463211685e-06 ],
        [ -1.0945789953889893 , -0.03547324181218445 , 0.3634668604715895 , -9.493485475182071e-05 ],
        [ -2.2979699638623208 , -0.024911009946374824 , 0.6188147891878977 , -0.00011857144326738345 ],
        [ -3.166418078853781 , -0.031958344183594084 , 0.8406859050306963 , -0.00012072892037607769 ],
        [ -4.013361936990245 , -0.038181594634872126 , 1.0869790557430794 , -0.00015667409032667326 ],
        [ -5.926346258972109 , -0.031964680347665884 , 1.5240375730724862 , -0.00016933328522669477 ],
        [ -7.90033309005282 , -0.020825728153169958 , 1.894831160970044 , -0.00010526577327102076 ],
        [ -7.152475999113026 , -0.015734400267369567 , 1.7102290556150588 , -0.0002686871465276408 ],
        [ 2.13934407431362 , -0.0320109908976759 , -0.5277579029594971 , -0.00025007383389563745 ],
        [ 1.8607414748455695 , -0.05343220981005764 , -0.33435348319941816 , -4.181194104351927e-05 ],
        [ 3.0620743202772447 , -0.031718714577417005 , -0.720332573077426 , 4.455859225551441e-05 ],
        [ 0.8715248291341835 , -0.014522140866670415 , -0.14027205100034504 , -8.962822645247351e-06 ],
    ]

    # Set the regression coefficient for clear creek
    ccr2sha_coeffs = [
        [ 1.792959080538937 , -0.059354492077861414 , -0.2960320376671016 , 5.046191161081086e-07 ],
        [ -0.7512647356887353 , -0.08260980765521911 , 0.3989728525297206 , -0.00010475304117792841 ],
        [ -2.2942016504506095 , -0.07724045046488413 , 0.7261448433127441 , -0.000105287304891731 ],
        [ -4.4991291842472005 , -0.07450124827763506 , 1.2838524899115868 , -0.00015759183944165977 ],
        [ -5.891695527938446 , -0.06622045217750946 , 1.5982456393750977 , -0.00015336390252308166 ],
        [ -9.012911351050086 , -0.05819889006557895 , 2.3237120979677863 , -0.00017729324265147367 ],
        [ -11.807214560294568 , -0.04586684354549874 , 2.88558316872268 , -9.275162445499182e-05 ],
        [ -10.475158908615029 , -0.04277016383124388 , 2.5906412962866536 , -0.00027885911475739506 ],
        [ 0.7463268269374388 , -0.06582341829199213 , -0.08213960515083435 , -0.0002556556527696652 ],
        [ 3.3789888674223043 , -0.08818106994859856 , -0.6214925066947776 , -2.910878547596829e-05 ],
        [ 8.334250042494999 , -0.06392276058932216 , -2.003215573079848 , 6.311125793283819e-05 ],
        [ 4.138461583069988 , -0.036555149612575485 , -0.8940381389461702 , 3.245473186093828e-05 ],
    ]

    # select the correct regression coefficient set based on location and river type
    if loc == 2:
        if ResSimRiver:
            coeffs = kes2sha_coeffs  # Use Keswick to Shasta coefficients when routing via river model
        
        else:
            coeffs = ccr2sha_coeffs  # Use CCR to Shasta coefficients for regression-only path
    
    else:
        raise ValueError("upstream_target: loc unknown, currently must be one of  {2,}")

    # Load datasets over rtw, must be same lengths
    # TODO: should check units, expecting CFS and C
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    dssFm = HecDss.open(forecastDSS)        
    tsc = dssFm.read(downstreamTT_rec, starttime_str, endtime_str, False).getData()
    tsc_int_times = tsc.times
    dtt = DSS_Tools.hectime_to_datetime(tsc)
    downstreamTT = tsc.values
    TT_units = str(tsc.units)
    
    # if the downstream target is in Fahrenheit, convert it to Celsius
    if TT_units.lower() == 'f' or TT_units.lower() == 'degF':
        for i, TT in enumerate(downstreamTT):
            downstreamTT[i] = (TT - 32.0) * 5.0 / 9.0  # Convert to C

        # Mark as C
        tsc.units = 'C'  
    
    # Close file
    dssFm.close()  

    # get the rest of the data over the same period
    eqTemp = DSS_Tools.data_from_dss(forecastDSS, eqTemp_rec, starttime_str, endtime_str)  # Daily EQ temps
    kesFlow = DSS_Tools.data_from_dss(forecastDSS, kesFlow_rec, starttime_str, endtime_str)  # Keswick daily flows
    sppFlow = DSS_Tools.data_from_dss(forecastDSS, sppFlow_rec, starttime_str, endtime_str)  # Spring Creek diversion daily flows

    # first daily eqTemp seem to often be bad, maybe timezone/DSS reading issue
    eqTemp[0] = eqTemp[1]  # Replace first value with second (original behavior preserved)

    shastaTT = []
    # for each downstream target temperature value, back-calculate the required
    # upstream (Shasta) temperature using month-blended regression coefficients
    for i, TT in enumerate(downstreamTT):
        mo_i = dtt[i].month - 1  # Current month index [0..11]
        mo_i_prev = mo_i - 1 if mo_i - 1 >= 0 else 11  # Previous month
        mo_i_next = mo_i + 1 if mo_i + 1 < 12 else 0  # Next month
        _, pFrac, cFrac, nFrac = fractional_month(dtt[i])  # Month weights

        a0, b0, c0, d0 = coeffs[mo_i_prev]  # Previous month coefficients
        a1, b1, c1, d1 = coeffs[mo_i]  # Current month coefficients
        a2, b2, c2, d2 = coeffs[mo_i]  # (Original code: next uses current coefficients; preserved)

        # if using the full ResSim river model, back-route hour-by-hour to solve for
        # the outlet temp; otherwise apply the simpler direct regression offset
        if ResSimRiver:
            i_future, hrs = get_step_future_and_RiverHrs(downstreamTT, kesFlow, i, loc)  # Step ahead and hours
            upstreamTT = backRouteWQTarget2(eqTemp, downstreamTT[i_future], TTDiff, hrs, i, loc, hour_of_day=10)  # Back-route to outlet temp
        else:
            upstreamTT = TT + TTDiff  # Regression-only path without routing
        print('result : ', TT, TTDiff, TT + TTDiff, upstreamTT)  # Diagnostic

        # Append the value into the time series
        shastaTT.append(upstreamTT)  # Collect result

    # copy over first record, which always fails maybe due to time zone read issues
    print('len(shastaTT)', len(shastaTT))  # Diagnostic length print
    shastaTT[0] = shastaTT[1]  # Replace first element with second (original behavior)

    # Write the regressed target to the DSS file
    dssFmRec = HecDss.open(forecastDSS)  # Open DSS for writing
    tsc.fullName = TT_W2_rec  # Set output path
    tsc.values = shastaTT  # Assign computed upstream TT values
    dssFmRec.write(tsc)  # Write to DSS
    dssFmRec.close()  # Close file


def get_downstream_loc(forecastDSS):
    """
    Read the configured downstream control location index.

    Reads from a dedicated integer DSS record.

    Parameters
    ----------
    forecastDSS : str
        Path to the DSS file containing the downstream control
        location record.

    Returns
    -------
    int
        The downstream control location index (e.g. 0 for Shasta
        Dam, or another value for a downstream point).
    """
    dssFm = HecDss.open(forecastDSS)        
    tsc = dssFm.get('//DOWNSTREAM_CONTROL_LOC///INTEGER/SACTRN_TARGET_TEMP/', True) # this should be passed in a linked record at some point
    loc = int(str(tsc.getText()).strip())
    print('Downstream Loc: ',str(tsc))
    dssFm.close()
    
    # Return the location to the calling function
    return loc


def forecast_data_preprocess_W2_5Res(currentAlternative, computeOptions):
    """
    Preprocess forecast boundary condition data for the CE-QUAL-W2
    5-Reservoir Sacramento/Trinity system.

    Fixes data types/units, applies met-data splicing, forecasts
    reservoir elevations, computes equilibrium temperature, creates a
    set of constant reference DSS records, and computes an upstream
    (Shasta) target temperature from the downstream target.

    This is the W2-specific counterpart to
    `forecast_data_preprocess_ResSim_5Res`: it performs a similar set
    of steps but omits the ResSim-specific air temperature lapse
    correction and gate-related constant records, and always uses the
    simpler (non-ResSim-river) upstream target calculation.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()` and `getTimeStep()` for logging and
        retrieving the balance period.
    computeOptions : object
        The compute options/settings object. Must support
        `getDssFilename()`, `getRunTimeWindow()`, and
        `getRunDirectory()`.

    Returns
    -------
    bool
        True once all preprocessing steps have completed
        successfully.

    Notes
    -----
    Workflow:

    1. Retrieves the DSS filename, run time window, run directory,
       and shared data directory for this compute.
    2. Fixes data types/units in the shared forecast DSS file via
       `DMS_preprocess.fix_DMS_types_units`.
    3. Splices in Redding meteorological data for Lewiston during
       Jan-Mar via `splice_met`.
    4. Computes equilibrium water temperature via `eq_temp`.
    5. Forecasts elevations for Shasta, Trinity, and Whiskeytown via
       `write_forecast_elevations`.
    6. Creates a set of constant-value reference DSS records used
       elsewhere in the model as fixed boundary conditions.
    7. Resamples the Keswick flow-release and target temperature
       records to daily resolution.
    8. Reads the configured downstream control location.
    9. If the downstream location is at Shasta Dam itself
       (`location == 0`), copies the target temperature directly.
       Otherwise, back-calculates the required upstream (Shasta)
       target temperature from the downstream target via
       `upstream_target`, using the simpler (non-ResSim-river)
       regression path.
    """
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    # Setup the anlysis items
    dss_file = computeOptions.getDssFilename()  # Compute DSS file
    rtw = computeOptions.getRunTimeWindow()  # Window for data

    run_dir = computeOptions.getRunDirectory()  # Run directory
    project_dir = Project.getCurrentProject().getProjectDirectory()  # Project directory
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)  # Log
    currentAlternative.addComputeMessage('run dir: ' + run_dir)  # Log
    balance_period = currentAlternative.getTimeStep()  # Time step
    shared_dir = os.path.join(project_dir, 'shared')  # Shared dir

    # output from scripting all goes to the <study>_Forecast.dss file.
    forecast_dss = os.path.join(shared_dir,'WTMP_SacTrn_Forecast.dss')
    DMS_preprocess.fix_DMS_types_units(forecast_dss)
    
    # splice in Redding met data for Lewiston during Jan-Mar
    splice_met(currentAlternative, rtw, forecast_dss, forecast_dss)

    # Compute/write EQ temperature
    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")
    eq_temp(rtw,
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Hour/sactrn_bc_script/"]
           )

    # forecast reservoir elevations for Shasta, Trinity, and Whiskeytown
    write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    # Create a set of fixed constant-value reference records 
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=13.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='WHI-target-13',fpart='WHI-target-13')

    # Keswick need daily record.
    DSS_Tools.resample_dss_ts(forecast_dss, '//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY', pad_start_days=1)  # Hour to Day
    DSS_Tools.resample_dss_ts(forecast_dss, '/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY')  # Hour to Day

    # read location from DSS
    location = get_downstream_loc(forecast_dss)  # Read downstream control location

    TT_rec = "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/"
    TT_W2_rec = "/USBR/SHASTA/TEMP-WATER-TARGET-W2-UPSTREAM//1Day/SACTRN_BC_SCRIPT/"
    # if the control point is at Shasta Dam itself, no back-routing is needed - copy directly
    if location == 0: 
        # @ Shasta Dam, use exact TT
        DSS_Tools.copy_dss_ts(TT_rec,new_dss_rec=TT_W2_rec,dss_file_path=forecast_dss,checkMakeCelsius=True)
    # otherwise, back-calculate the upstream (Shasta) target temperature needed to meet the downstream target
    else:
        # At any location other than Shasta. Apply the regression to move the location
        upstream_target(forecast_dss, rtw,
                        "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/",
                        "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Day/sactrn_bc_script/",
                        "//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Day/SACTRN_BC_SCRIPT/",
                        "/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1Day/SACTRN_BC_SCRIPT/",
                        location, TT_W2_rec, ResSimRiver=False)  # Regression/back-routing path

    # Return that the setup was successful
    return True
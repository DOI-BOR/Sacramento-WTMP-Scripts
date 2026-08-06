"""
Forecast Data Preprocess

This module orchestrates several preprocessing steps needed before running
WTMP/ResSim/W2 forecast workflows. It includes:

- **Meteorology and equilibrium temperature**: Reading met data, computing
  equilibrium water temperature (`equilibrium_temp`), and writing timeseries at
  multiple intervals.
- **Reservoir elevation derivations**: Converting storage to elevation using
  tabulated elevation–storage–area data, inventing elevation records where
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
from hec.io import DSSIdentifier  # DSS pathname identifier helper (A–F part parsing/building)
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
for p in sys.path:
    if any(phrase in p for phrase in search_list):  # Flag entries containing any study phrase
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
    
    # Handle behavior based on the input type
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

        if x < xp[0]:
            # Left extrapolation
            return left  
        
        elif x > xp[-1]:
            # Right extrapolation
            return right  
        
        else:
            # Point is in the middle. Look for the point that brackets it.
            for i in range(len(xp) - 1):
                
                # Find bracketing interval
                if x >= xp[i] and x <= xp[i + 1]:  
                    # Perform the linear interpolation
                    t = (x - xp[i]) / (xp[i + 1] - xp[i])  # Fraction within interval
                    
                    # Return linear interpolation result
                    return fp[i] + t * (fp[i + 1] - fp[i])  


def eq_temp(rtw, at, cl, ws, sr, td, eq_temp_out):
    """
    Compute equilibrium water temperature series and write multiple intervals.

    Parameters
    ----------
    rtw : object
        Run time window with `getStartTimeString()` / `getEndTimeString()`.
    at : [str, str]
        Pair `[dss_file_path, at_path]` for air temperature (°C).
    cl : [str, str]
        Pair `[dss_file_path, cloud_path]` for cloud fraction (0–1).
    ws : [str, str]
        Pair `[dss_file_path, wind_speed_path]` for wind speed (m/s).
    sr : [str, str]
        Pair `[dss_file_path, solar_radiation_path]` for shortwave irradiance (W/m²).
    td : [str, str]
        Pair `[dss_file_path, dewpoint_path]` for dew point (°C).
    eq_temp_out : [str, str]
        Pair `[output_dss_file_path, output_path]` for the equilibrium temperature series.

    Returns
    -------
    None

    Notes
    -----
    - Computes equilibrium temperature with `equilibrium_temp.calc_equilibrium_temp`.
    - Writes instantaneous values and standardized 1-day / 1-week transforms to DSS.
    """
    
    # Get the time window 
    starttime_str = rtw.getStartTimeString()  # Window start string
    endtime_str = rtw.getEndTimeString()      # Window end string

    # Get at data and times in the formats needed
    dssFm = HecDss.open(at[0])  # Open DSS file for air temperature
    tsc = dssFm.read(at[1], starttime_str, endtime_str, False).getData()  # Read AT data container
    tsc_int_times = tsc.times  # HEC numeric times for direct writeback
    dtt = DSS_Tools.hectime_to_datetime(tsc)  # Python datetimes for solver inputs
    at_data = tsc.values  # Air temperature values
    dssFm.close()  # Close input file

    # get the rest of the data over the same period
    cl_data = DSS_Tools.data_from_dss(cl[0], cl[1], starttime_str, endtime_str)  # Cloud fraction
    ws_data = DSS_Tools.data_from_dss(ws[0], ws[1], starttime_str, endtime_str)  # Wind speed
    sr_data = DSS_Tools.data_from_dss(sr[0], sr[1], starttime_str, endtime_str)  # Solar radiation
    td_data = DSS_Tools.data_from_dss(td[0], td[1], starttime_str, endtime_str)  # Dew point

    # calc_equilibrium_temp(dtt, at, cl, sr, td, ws)
    Te = equilibrium_temp.calc_equilibrium_temp(dtt, at_data, cl_data, sr_data, td_data, ws_data)  # Solver returns °C

    # Output the container
    tsc = TimeSeriesContainer()  # Prepare container for instantaneous series
    tsc.times = tsc_int_times  # Assign times (HEC numeric)
    tsc.fullName = eq_temp_out[1]  # Destination path
    tsc.values = Te  # Equilibrium temperature values
    tsc.units = 'C'  # Degrees Celsius
    tsc.type = 'INST-VAL'  # Instantaneous values
    tsc.numberValues = len(tsc.values)  # Count

    tsm = tsmath(tsc)  # Wrap in TimeSeriesMath for transforms
    tsm_day = DSS_Tools.standardize_interval(tsm, '1day')  # Standardize to 1-day interval
    tsm_wk = DSS_Tools.standardize_interval(tsm, '1week')  # Standardize to 1-week interval

    dssFmOut = HecDss.open(eq_temp_out[0])  # Open output DSS file
    dssFmOut.write(tsc)  # Write instantaneous series
    dssFmOut.write(tsm_day)  # Write daily series
    dssFmOut.write(tsm_wk)  # Write weekly series
    dssFmOut.close()  # Close output


def storage_to_elev(res_name, elev_stor_area, forecast_dss, storage_rec, conic=False):
    """
    Convert a storage time series to elevation using tabulated curves.

    Parameters
    ----------
    res_name : str
        Reservoir name for the output C-part.
    elev_stor_area : dict
        Curve data with keys `'stor'` and `'elev'` (lists of breakpoints).
    forecast_dss : str
        DSS file path containing the storage record.
    storage_rec : str
        DSS pathname for the input storage record.
    conic : bool, optional
        Placeholder for conic interpolation (currently unsupported).

    Returns
    -------
    None

    Raises
    ------
    SystemExit
        If `conic=True` (not implemented).
    """
    
    # Open and read the content
    dssFmRec = HecDss.open(forecast_dss)  # Open DSS file
    tsc = dssFmRec.get(storage_rec, True)  # read ALL data in record

    # Create a holder for the data
    elev = []  
    
    # Handle the interpolation
    if conic:
        # Use conic interpolation
        print('Conic interpolation of elevations from storage not supported yet.')
        sys.exit(-1)  # Preserve original error path
    else:
        # Use linear interpolation
        for j in range(tsc.numberValues):
            elev.append(cbfj.linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], tsc.values[j]))  # storage→elev
            print('stor2elev: ', j, tsc.times[j], tsc.values[j])  # Diagnostic

    # Write the output to the dss file
    recparts = tsc.fullName.split('/')  # Tokenize pathname for edits
    recparts[2] = res_name  # Set C-part to reservoir name
    recparts[3] = 'ELEV'  # Set parameter to ELEV
    tsc.fullName = '/'.join(recparts)  # Reassemble pathname
    tsc.units = 'ft'  # Feet
    tsc.type = 'INST-VAL'  # Instantaneous values
    tsc.values = elev  # Assign computed elevations
    dssFmRec.write(tsc)  # Write to DSS
    dssFmRec.close()  # Close file


def invent_elevation(res_name, forecast_dss, storage_rec, elev_constant_ft):
    """
    Invent an elevation time series from a storage record by constant fill.

    Parameters
    ----------
    res_name : str
        Reservoir name for C-part.
    forecast_dss : str
        DSS file path containing input record.
    storage_rec : str
        DSS pathname of a storage record whose timebase is reused.
    elev_constant_ft : float
        Constant elevation (ft) used to fill all timesteps.

    Returns
    -------
    None
    """
    
    # Open and read the data
    dssFmRec = HecDss.open(forecast_dss)  # Open DSS
    tsc = dssFmRec.get(storage_rec, True)  # Read entire record
    recparts = tsc.fullName.split('/')  # Tokenize pathname
    recparts[2] = res_name  # Set reservoir name
    recparts[3] = 'ELEV'  # Parameter as ELEV
    tsc.fullName = '/'.join(recparts)  # Reassemble pathname
    tsc.units = 'ft'  # Feet
    tsc.type = 'INST-VAL'  # Instantaneous
    
    # Fill with a constant value
    tsc.values = [elev_constant_ft for j in range(tsc.numberValues)]
    
    # Write and close the series
    dssFmRec.write(tsc)  # Write to DSS
    dssFmRec.close()  # Close file


def write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir):
    """
    Generate and write forecast reservoir elevations for Shasta, Trinity, Whiskeytown.

    Parameters
    ----------
    currentAlternative : object
        Alternative used for logging messages.
    rtw : object
        Run time window object with start/end strings.
    forecast_dss : str
        Forecast DSS file path for reads/writes.
    shared_dir : str
        Path to shared directory containing elevation–storage–area CSVs.

    Returns
    -------
    None

    Notes
    -----
    - Converts monthly and daily storage to elevation using tabulated curves.
    - Invents hourly/daily elevation placeholders for timing.
    - Predicts daily elevation using a balance-flow approach from inflows/outflows.
    """
    
    # Get date for starting elevation - look for end-of-month before start time
    ht = HecTime(rtw.getStartTimeString())  # Parse HEC start time
    
    if ht.month() == 1:
        start_str = dt.datetime(ht.year() - 1, 12, 31).strftime('%d%b%Y') + ' 2400'  # Previous year-end
    else:
        start_dt = dt.datetime(ht.year(), ht.month(), 1)  # First of current month
        start_dt = start_dt - dt.timedelta(days=1)  # Back up one day (EOM)
        start_str = start_dt.strftime('%d%b%Y') + ' 2400'  # End-of-month time string
    end_str = rtw.getEndTimeString()  # End of window

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
    
    outflow_records = ['/SACRAMENTO RIVER/SHASTA LAKE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']  # Release as outflow
    
    starting_elevation = DSS_Tools.first_value(forecast_dss, '/SACRAMENTO RIVER/SHASTA LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/', start_str, end_str)  # First elev
    print('starting_elevation ', starting_elevation)  # Diagnostic
    
    cbfj.predict_elevation(currentAlternative, start_str, end_str, 'Shasta Lake', inflow_records, outflow_records, starting_elevation,
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

    outflow_records = ['/TRINITY RIVER/TRINITY LAKE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']  # Outflow

    starting_elevation = DSS_Tools.first_value(forecast_dss, '/TRINITY RIVER/TRINITY LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/', start_str, end_str)
    print('starting_elevation ', starting_elevation)  # Diagnostic

    cbfj.predict_elevation(currentAlternative, start_str, end_str, 'Trinity Lake', inflow_records, outflow_records, starting_elevation,
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

    starting_elevation = DSS_Tools.first_value(forecast_dss, '/CLEAR CREEK/WHISKEYTOWN LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/', start_str, end_str)
    print('starting_elevation ', starting_elevation)  # Diagnostic

    cbfj.predict_elevation(currentAlternative, start_str, end_str, 'Whiskeytown Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Whiskeytown Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')  # Predict daily elevation


def splice_met(currentAlternative, rtw, forecast_dss, output_dss):
    """
    Splice meteorological records for Lewiston using Redding data in Jan–Mar.

    Parameters
    ----------
    currentAlternative : object
        Alternative used for logging messages.
    rtw : object
        Run time window.
    forecast_dss : str
        Input DSS file containing met records.
    output_dss : str
        Output DSS file (may be same as input) for spliced results.

    Returns
    -------
    None

    Notes
    -----
    - Replaces Lewiston met records for the months `[1, 2, 3]` with selected
      Redding records using `DSS_Tools.replace_data`.
    - Standardizes interval to `'1HOUR'` for consistent splicing coverage.
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
    Derive study directory path from a run directory.

    Parameters
    ----------
    run_dir : str
        Path to a specific run directory.

    Returns
    -------
    str
        Path to the top-level study directory.
    """
    
    # Remove last segment
    w2sim, _ = os.path.split(run_dir)  
    
    # Remove 'w2sim' segment
    runs_dir, _ = os.path.split(w2sim)  
    
    # Remove 'runs' segment
    study_dir, _ = os.path.split(runs_dir)  
    
    # Return the top-level study dir
    return study_dir  


def model_dir_from_run_dir(run_dir, model_place, model_name):
    """
    Construct a CE-QUAL-W2 model directory from the run directory.

    Parameters
    ----------
    run_dir : str
        Path to a specific run directory.
    model_place : str
        Subfolder under 'cequal-w2' indicating the model place (e.g., 'SacTrn').
    model_name : str
        Model name folder under `model_place`.

    Returns
    -------
    str
        Full path to the W2 model directory.
    """
    
    # Ascend to study root
    study_dir = study_dir_from_run_dir(run_dir)  
    
    # Build path
    model_dir = os.path.join(study_dir, 'cequal-w2', model_place, model_name)  
    
    # Return model directory path
    return model_dir  


def forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions):
    """
    Preprocess data for the 5-reservoir ResSim forecast workflow.

    Parameters
    ----------
    currentAlternative : object
        Alternative context for logging and options retrieval.
    computeOptions : object
        Contains DSS filename, run directory, time window, etc.

    Returns
    -------
    bool
        True on successful preprocessing.

    Notes
    -----
    Steps include DMS type/unit fixup, met lapse, met splicing, elevation writes,
    equilibrium temperature computation, constant record creation, relative humidity
    derivation, resampling for Keswick releases and targets, and upstream target
    assembly based on downstream location.
    """
    
    # Get the setup information for the analysis
    dss_file = computeOptions.getDssFilename()  # Current compute DSS file
    rtw = computeOptions.getRunTimeWindow()     # Runtime window

    run_dir = computeOptions.getRunDirectory()  # Path to this compute’s run directory
    project_dir = Project.getCurrentProject().getProjectDirectory()  # Project dir for study
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)  # Log path
    currentAlternative.addComputeMessage('run dir: ' + run_dir)  # Log path
    balance_period = currentAlternative.getTimeStep()  # Current time step (not directly used)
    shared_dir = os.path.join(project_dir, 'shared')  # Shared artifacts folder

    forecast_dss = os.path.join(shared_dir, 'WTMP_SacTrn_Forecast.dss')  # Forecast DSS path
    DMS_preprocess.fix_DMS_types_units(forecast_dss)  # Ensure types/units are expected formats

    # Calculate meteorological airtemp lapse for the elevation at Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: ' + forecast_dss)
    currentAlternative.addComputeMessage('lapse outfile: ' + forecast_dss)
    
    # Apply lapse rate to AT
    DSS_Tools.airtemp_lapse(forecast_dss, "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-AIR//1HOUR/SACTRN_BC_SCRIPT/",
                  0.7, forecast_dss, "Shasta_Lapse")  

    # Splice met records for Lewiston
    splice_met(currentAlternative, rtw, forecast_dss, forecast_dss)  

    # Generate forecast elevations
    write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)  

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
                        dss_type='PER-AVER', period='1DAY', cpart='TENS', fpart='TENS')  # Daily temp 10°C
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
    DSS_Tools.resample_dss_ts(forecast_dss, '//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY', pad_start_days=1)  # Hour→Day
    DSS_Tools.resample_dss_ts(forecast_dss, '/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY')  # Hour→Day targets

    # read location from DSS
    location = get_downstream_loc(forecast_dss)  # 0=Dam, others downstream points

    # Set the target paths
    TT_rec = "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/"  # Downstream TT
    TT_W2_rec = "/USBR/SHASTA/TEMP-WATER-TARGET-W2-UPSTREAM//1Day/SACTRN_BC_SCRIPT/"  # Upstream TT out
    
    # Apply the regression to account for warming between the release and the 
    if location == 0:
        # At Shasta Dam, use exact TT without adjustment
        DSS_Tools.copy_dss_ts(TT_rec, new_dss_rec=TT_W2_rec, dss_file_path=forecast_dss, checkMakeCelsius=True)  # Copy/convert to °C
    
    else:
        # The compliacne location is downstream. Apply the regression.
        upstream_target(forecast_dss, rtw,
                        "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/",
                        "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Day/sactrn_bc_script/",
                        "//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Day/SACTRN_BC_SCRIPT/",
                        "/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1Day/SACTRN_BC_SCRIPT/",
                        location, TT_W2_rec, ResSimRiver=False)  # Build upstream TT via regression/back-routing

    # Return that the setup is successful
    return True  


def route_downstream(tKeswick, keswickFlowDaily, eqTempDaily, step, hour_of_day, loc):
    """
    Forward route temperature downstream using a simple exchange-with-atmosphere model.

    Parameters
    ----------
    tKeswick : float
        Temperature at Keswick [°C].
    keswickFlowDaily : list of float
        Daily Keswick flows [CFS].
    eqTempDaily : list of float
        Daily equilibrium temperatures [°C].
    step : int
        Current day index.
    hour_of_day : int
        Hour within day at which to start routing (0–23).
    loc : int
        Downstream location index (1=Hwy44, 2=CCR, 3=Ball's Ferry).

    Returns
    -------
    float
        Routed temperature at the downstream location [°C].

    Notes
    -----
    Applies a fixed exchange coefficient (0.015 per hour) and steps through
    `hrs` travel time computed by `travel_time_hrs`.
    """
    
    exchCoef = 0.015  # hourly exchange rate between atmosphere and river temp

    hrs = travel_time_hrs(loc, keswickFlowDaily[step])  # loc 2 = CCR
    t = tKeswick  # Initialize routed temp
    i = step  # Start at current day
    imax = len(eqTempDaily)  # Max index bound
    h = hour_of_day  # Current hour
    
    for k in range(hrs):
        deltaTemp = (eqTempDaily[i] - t) * exchCoef  # Relax toward eq temp
        t += deltaTemp  # Update

        if h == 23:
            h = 0  # Wrap hour
            i = min(imax - 1, i + 1)  # Advance day with bound

        else:
            h += 1  # Next hour

    return t  # Routed downstream temperature


def travel_time_hrs(loc, keswickFlow):
    """
    Estimate integer travel time (hours) from Keswick to a downstream location.

    Parameters
    ----------
    loc : int
        Location index: 1=Highway 44, 2=CCR, 3=Ball's Ferry.
    keswickFlow : float
        Keswick release flow [CFS].

    Returns
    -------
    int
        Travel time in whole hours.

    Raises
    ------
    NotImplementedError
        For unsupported location indices.
    """
    
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
    Compute fractional month and weights for previous/current/next months.

    Parameters
    ----------
    date_obj : datetime.datetime
        Timestamp to analyze.

    Returns
    -------
    tuple
        `(fractional_month, fractional_previous, fractional_current, fractional_next)`:
        - fractional_month : float
            Fraction into month from 0 to <1.
        - fractional_previous : float
            Weight applied to previous month.
        - fractional_current : float
            Weight applied to current month.
        - fractional_next : float
            Weight applied to next month.

    Notes
    -----
    Weights sum to 1.0 with a pivot at 0.5 (mid-month).
    """

    year = date_obj.year  # Year of date
    month = date_obj.month  # Month number
    day = date_obj.day  # Day of month

    # Get the number of days in the current month using the calendar module
    _, days_in_month = calendar.monthrange(year, month)  # (weekday_of_first, days)

    # Calculate the fractional month. Use 1.0 to force float division.
    fractional_month = (day - 1.0) / days_in_month  # Fraction [0, 1)
    fractional_previous = 0.0  # Default previous weight
    fractional_next = 0.0  # Default next weight

    if fractional_month > 0.5:
        fractional_current = (1.0 - fractional_month) + 0.5  # Late month: more current+next
        fractional_next = 1.0 - fractional_current  # Remaining weight goes to next

    else:
        fractional_current = fractional_month + 0.5  # Early month: more previous+current
        fractional_previous = 1.0 - fractional_current  # Remaining weight goes to previous

    # Return to the calling function
    return fractional_month, fractional_previous, fractional_current, fractional_next  # Weights tuple


def get_step_future_and_RiverHrs(wqTargetDaily, keswickFlowDaily, step, loc):
    """
    Determine future step index and river travel hours given flow and location.

    Parameters
    ----------
    wqTargetDaily : list of float
        Daily downstream WQ target temperatures [°C].
    keswickFlowDaily : list of float
        Daily Keswick flows [CFS].
    step : int
        Current step (day index).
    loc : int
        Downstream location index (e.g., 2=CCR).

    Returns
    -------
    (int, int)
        `(step_future, hrs)` where `step_future` is the future day index after
        travel/residence and `hrs` is the integer river travel time in hours.

    Notes
    -----
    Combines `travel_time_hrs` with a Keswick pool fraction representing residence
    time effects (bounded by calibration multiplier).
    """

    # Power law approximation for velocity in the Sacramento River
    hrs = travel_time_hrs(loc, keswickFlowDaily[step])  # Compute travel hours

    # Get Keswick pool information - from forecast TCD script
    flowVol = keswickFlowDaily[step] * 86400.  # Daily volume [ft³]
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
    Back-route downstream target temperature to required outlet temperature at Shasta.

    Parameters
    ----------
    eqTempDaily : list of float
        Daily equilibrium temperatures [°C].
    targetTempFuture : float
        Downstream target temperature in the future [°C].
    sha2kes_diff : float
        Temperature difference term between Shasta and Keswick (diff may be negative if heating).
    hrs : int
        River travel time in hours.
    step : int
        Current day index.
    loc : int
        Downstream location index.
    hour_of_day : int, optional
        Hour within day to start routing (default: 10).

    Returns
    -------
    float
        Required outlet temperature at Shasta [°C] to meet downstream target.

    Raises
    ------
    ValueError
        If bracketing fails and the target cannot be interpolated safely (except certain cases).

    Notes
    -----
    - Searches over a bracket [tSearchMin, tSearchMax] using a fixed number of iterations.
    - Switches to simple assignment when `cantBeMet=True`; otherwise linear interpolation of bracket.
    """
    
    # Set the initial search information
    exchCoef = 0.015  # exchange rate between atmosphere and river temp
    tSearchMin = targetTempFuture + sha2kes_diff - 6.  # Lower search bound
    tSearchMax = tSearchMin + 18.  # Upper search bound
    numIters = 51  # Number of samples
    bracketed = False  # Flag for successful bracketing
    cantBeMet = False  # Flag for infeasible target
    
    # Loop and search until converged or the number of iterations are met
    for j in range(numIters):
        # Candidate outlet temp
        outletTemp = tSearchMin + float(j) / float(numIters + 1) * (tSearchMax - tSearchMin)  
        
        # Impact of Keswick
        t = outletTemp - sha2kes_diff  # sha2kes_diff is negative if heating in Keswick
        t1 = 1  # Placeholder variable retained from original
        
        # Route downstream
        i = step  # Day index
        imax = len(eqTempDaily)  # Bounds
        h = hour_of_day  # Hour to start

        # Loop over the hours 
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
            prevT = t  # Establish previous
            prevOutletT = outletTemp  # Track previous candidate

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

    # Validate the estimate from the search
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

    # Return the required outlet temp at Shasta
    return targetTemp  


def upstream_target(forecastDSS, rtw, downstreamTT_rec, eqTemp_rec, kesFlow_rec, sppFlow_rec, loc, TT_W2_rec, ResSimRiver=True):
    """
    Construct upstream target temperature at Shasta from downstream targets.

    Parameters
    ----------
    forecastDSS : str
        DSS file path for forecast datasets.
    rtw : object
        Run time window.
    downstreamTT_rec : str
        DSS pathname for downstream target temperature (daily).
    eqTemp_rec : str
        DSS pathname for daily equilibrium temperatures used in routing.
    kesFlow_rec : str
        DSS pathname for Keswick daily flows [CFS].
    sppFlow_rec : str
        DSS pathname for Spring Creek diversion flows (used in regressions).
    loc : int
        Downstream location index (currently must be 2).
    TT_W2_rec : str
        Output DSS pathname for upstream target temperature.
    ResSimRiver : bool, optional
        If True, apply ResSim river routing; else use regression direct difference.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If `loc` is not among supported indices (presently `{2}` only).

    Notes
    -----
    - Applies month-dependent regression coefficients to compute temperature
      difference (TTDiff) between downstream and upstream.
    - When `ResSimRiver` is True, back-routes downstream targets accounting
      for travel time/atmospheric exchange; otherwise adds TTDiff directly.
    """
    
    # Set the regression coeffients for keswick
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

    # todo: generalize this so that it handles all model coupling
    if loc == 2:
        if ResSimRiver:
            coeffs = kes2sha_coeffs  # Use Keswick→Shasta coefficients when routing via river model
        
        else:
            coeffs = ccr2sha_coeffs  # Use CCR→Shasta coefficients for regression-only path
    
    else:
        raise ValueError("upstream_target: loc unknown, currently must be one of  {2,}")

    # Load datasets over rtw, must be same lengths
    # TODO: should check units, expecting CFS and C
    # Get the times of the series
    starttime_str = rtw.getStartTimeString()  # Window start
    endtime_str = rtw.getEndTimeString()  # Window end
    
    # Open the file
    dssFm = HecDss.open(forecastDSS)  # Open forecast DSS
    
    # Read the data
    tsc = dssFm.read(downstreamTT_rec, starttime_str, endtime_str, False).getData()  # Read downstream TT
    tsc_int_times = tsc.times  # Save times for writing
    
    # Convert the datetimes
    dtt = DSS_Tools.hectime_to_datetime(tsc)  # Convert to Python datetimes
    
    # Read the values and units
    downstreamTT = tsc.values  # Values
    TT_units = str(tsc.units)  # Units string
    
    # Convert from F to C, if necessary
    if TT_units.lower() == 'f' or TT_units.lower() == 'degF':
        for i, TT in enumerate(downstreamTT):
            downstreamTT[i] = (TT - 32.0) * 5.0 / 9.0  # Convert to °C

        # Mark as °C
        tsc.units = 'C'  
    
    # Close file
    dssFm.close()  

    # get the rest of the data over the same period
    eqTemp = DSS_Tools.data_from_dss(forecastDSS, eqTemp_rec, starttime_str, endtime_str)  # Daily EQ temps
    kesFlow = DSS_Tools.data_from_dss(forecastDSS, kesFlow_rec, starttime_str, endtime_str)  # Keswick daily flows
    sppFlow = DSS_Tools.data_from_dss(forecastDSS, sppFlow_rec, starttime_str, endtime_str)  # Spring Creek diversion daily flows

    # first daily eqTemp seem to often be bad, maybe timezone/DSS reading issue
    eqTemp[0] = eqTemp[1]  # Replace first value with second (original behavior preserved)

    # Upstream target temperature series
    shastaTT = []  
    for i, TT in enumerate(downstreamTT):
        mo_i = dtt[i].month - 1  # Current month index [0..11]
        mo_i_prev = mo_i - 1 if mo_i - 1 >= 0 else 11  # Previous month
        mo_i_next = mo_i + 1 if mo_i + 1 < 12 else 0  # Next month
        _, pFrac, cFrac, nFrac = fractional_month(dtt[i])  # Month weights

        a0, b0, c0, d0 = coeffs[mo_i_prev]  # Previous month coefficients
        a1, b1, c1, d1 = coeffs[mo_i]  # Current month coefficients
        a2, b2, c2, d2 = coeffs[mo_i]  # (Original code: next uses current coefficients; preserved)

        print(i, ' : ', dtt[i].year, dtt[i].month, dtt[i].day, dtt[i].hour, ' : ', TT, eqTemp[i], kesFlow[i], sppFlow[i])  # Diagnostic
        upstreamTT = -1  # Placeholder initial value

        TTDiff = pFrac * (a0 + b0 * eqTemp[i] + c0 * math.log10(kesFlow[i]) + d0 * sppFlow[i]) + \
                 cFrac * (a1 + b1 * eqTemp[i] + c1 * math.log10(kesFlow[i]) + d1 * sppFlow[i]) + \
                 nFrac * (a2 + b2 * eqTemp[i] + c2 * math.log10(kesFlow[i]) + d2 * sppFlow[i])  # Month-weighted difference

        # Adjust based on the target location
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
    Read downstream control location index from DSS.

    Parameters
    ----------
    forecastDSS : str
        DSS file containing the downstream control location record.

    Returns
    -------
    int
        Location index parsed from the DSS record text.

    Notes
    -----
    Expects the linked record `//DOWNSTREAM_CONTROL_LOC///INTEGER/SACTRN_TARGET_TEMP/`
    to contain an integer in its text field.
    """
    
    # Open the DSS file
    dssFm = HecDss.open(forecastDSS)
    
    # Get teh compliance location
    tsc = dssFm.get('//DOWNSTREAM_CONTROL_LOC///INTEGER/SACTRN_TARGET_TEMP/', True)  # this should be passed in a linked record at some point
    
    # Parse out the values
    loc = int(str(tsc.getText()).strip())  # Parse integer from text
    print('Downstream Loc: ', str(tsc))  # Diagnostic print of container
    
    # Close the DSS file
    dssFm.close()
    
    # Return the location to the calling function
    return loc


def forecast_data_preprocess_W2_5Res(currentAlternative, computeOptions):
    """
    Preprocess data for the 5-reservoir W2 forecast workflow.

    Parameters
    ----------
    currentAlternative : object
        Alternative context for logging and options retrieval.
    computeOptions : object
        Contains DSS filename, run directory, time window, etc.

    Returns
    -------
    bool
        True on successful preprocessing.

    Notes
    -----
    Similar to the ResSim preprocessing, but tailored for W2 with specific
    resampling and upstream TT assembly flows.
    """
    
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
    forecast_dss = os.path.join(shared_dir, 'WTMP_SacTrn_Forecast.dss')  # Forecast DSS
    DMS_preprocess.fix_DMS_types_units(forecast_dss)  # Fix types/units

    splice_met(currentAlternative, rtw, forecast_dss, forecast_dss)  # Splice met

    # Compute/write EQ temperature
    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")
    eq_temp(rtw,
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss, "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Hour/sactrn_bc_script/"]
           )  

    # Write the elevation
    write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    # Create the constant DSS records as placeholders
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow',
                        dss_type='PER-AVER', period='1DAY', cpart='TinyFlow', fpart='TinyFlow')  # Tiny daily flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow',
                        dss_type='PER-AVER', period='1HOUR', cpart='TinyFlow', fpart='TinyFlow')  # Tiny hourly flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water',
                        dss_type='PER-AVER', period='1DAY', cpart='TENS', fpart='TENS')  # Daily 10°C temp
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow',
                        dss_type='PER-AVER', period='1DAY', cpart='ZEROS', fpart='ZEROS')  # Daily zero flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow',
                        dss_type='PER-AVER', period='1HOUR', cpart='ZEROS', fpart='ZEROS')  # Hourly zero flow
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=13.0, what='temp-water',
                        dss_type='PER-AVER', period='1DAY', cpart='WHI-target-13', fpart='WHI-target-13')  # WHI 13°C target

    # Keswick need daily record.
    DSS_Tools.resample_dss_ts(forecast_dss, '//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY', pad_start_days=1)  # Hour→Day
    DSS_Tools.resample_dss_ts(forecast_dss, '/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/', rtw, forecast_dss, '1DAY')  # Hour→Day

    # read location from DSS
    location = get_downstream_loc(forecast_dss)  # Read downstream control location

    # Define the locations for the records
    TT_rec = "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/"  # Downstream TT record
    TT_W2_rec = "/USBR/SHASTA/TEMP-WATER-TARGET-W2-UPSTREAM//1Day/SACTRN_BC_SCRIPT/"  # Upstream TT out path
    
    # Create the target timeseries
    if location == 0:
        # At Shasta Dam, use exact TT without modification
        DSS_Tools.copy_dss_ts(TT_rec, new_dss_rec=TT_W2_rec, dss_file_path=forecast_dss, checkMakeCelsius=True)
    
    else:
        # At any location other than Shasta. Apply the regression to move the location
        upstream_target(forecast_dss, rtw,
                        "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/",
                        "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Day/sactrn_bc_script/",
                        "//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Day/SACTRN_BC_SCRIPT/",
                        "/CLEAR CREEK/WHOW-DIVERSION-SPRING-CR//1Day/SACTRN_BC_SCRIPT/",
                        location, TT_W2_rec, ResSimRiver=False)  # Regression/back-routing path

    # Return that the setup was successful
    return True
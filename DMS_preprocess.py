"""
Pre-Process: ResSim/W2 (Sacramento–Trinity WTMP)
===============================================

Purpose
-------
Utilities and workflow steps to:
- Sanitize and standardize DMS-provided DSS time series (types/units).
- Create derived outflow, plotting, and river-balance records for ResSim/W2.
- Prepare CE-QUAL-W2-specific withdrawal flows for Shasta TCD.
- Splice Lewiston meteorological data with Redding records for winter months.

Notes
-----
- **Environment:** Jython in HEC‑WAT (Python 2.7 semantics), using HEC‑DSS APIs.
- **I/O:** DSS paths may be specified as `file::/A/B/C/D/E/F/` to target
  a specific file; otherwise the primary file handle is used.
- **Units:** Flow CMS→CFS factor `35.314666213`; wind conversions (kph↔m/s);
  temperature standardized to Celsius where required for model linking.

See Also
--------
DSS_Tools
    Local utilities used throughout (safe reads, resampling, transformations).
Acc_Dep_ResSim_SacTrn
    Accretion/depletion helpers for Sacramento/Trinity context (imported).
"""

# --- Imports ---------------------------------------------------------------------------
from hec.heclib.dss import HecDss                     # HEC-DSS: open/read/write time series and metadata
from hec.hecmath import HecMathException              # HEC math exception type for transform/read failures
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE   # Sentinel for undefined numeric values within HEC containers
from hec.heclib.util import HecTime                   # HEC time conversion utilities (string ↔ HEC int time)
from hec.io import DSSIdentifier                      # DSS pathname identifier (retained per original code)
from hec.io import TimeSeriesContainer                # Container class for writing time series back to DSS
import hec.hecmath.TimeSeriesMath as tsmath           # Time-series transform helpers (e.g., interval changes)
from rma.util.RMAConst import MISSING_DOUBLE          # RMA constant for missing double values (import retained)
import math                                           # Standard math library; used across formulas and checks
import sys                                            # System utilities; used for path edits and sys.exit on errors
import datetime as dt                                 # Datetime utilities (retained; not used in all functions)
import os, sys                                        # Filesystem and system path utilities (combined import per original)

from com.rma.io import DssFileManagerImpl             # RMA DSS manager (retained for environment context)
from com.rma.model import Project                     # WAT project API to resolve workspace and project directories
#import hec.hecmath.TimeSeriesMath as tsmath          # Original commented import (kept intact)

from com.rma.io import DssFileManagerImpl             # Duplicate import retained per original file content
from java.util import TimeZone                        # Java TimeZone (retained; not directly used)

# --- sys.path hygiene for WAT workspace scripts ----------------------------------------
# create list of unwanted folders in sys.path
search_list = ["SacTrn", "Sacramento", "American", "Stanislaus"]  # Keywords to exclude from module search

# initialize and search for unwanted paths
matching_paths = []
for p in sys.path:
    if any(phrase in p for phrase in search_list):    # Detect undesired paths matching keywords
        matching_paths.append(p)

# remove matching paths from sys.path
for path in matching_paths:
    if path in sys.path:
        sys.path.remove(path)                         # Prune workspace search paths to avoid collisions

# append path
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))  # Add WAT scripts folder

# --- Local module imports (workspace scripts) -------------------------------------------

import Acc_Dep_ResSim_SacTrn                           # Accretion/depletion helpers for ResSim (Sac/Trn)
reload(Acc_Dep_ResSim_SacTrn)                          # Jython: ensure the latest version is loaded

import DSS_Tools                                       # Local DSS utilities used heavily in preprocessing
reload(DSS_Tools)                                      # Reload to pick up recent edits during interactive runs

# Units that need conversion/standardization for model linking
units_need_fixing = ['tenths','deg','kph'] #'radians',]  # 'radians' handled in legacy function below

def fix_DMS_types_units(dss_file):
    """
    Standardize DMS record types and units in-place.

    This method was implemented to change data types to `PER-AVER`
    for records that are not coming from the DMS that way, and to
    normalize a handful of oddball unit labels into standard/derived
    forms.

    Parameters
    ----------
    dss_file : str
        Path to the source DSS file (DMS hydro/met time series).

    Returns
    -------
    None
        Writes corrected records back to `dss_file`.

    Notes
    -----
    - Sets `PER-AVER` for `/flow` or `/1Day/` (non-storage) records.
    - Converts unit codes:
      * `tenths` → `FRAC` (values ÷ 10)
      * `radians` → `deg` (values × 360 / 2π)
      * `deg` → `radians` (values × 2π / 360)
      * `kph` → `m/s` (values ÷ 3.6)
    - Skips paired, scalar/text, and special mapping/control paths.
    - Unlike `fix_DMS_types_units_old`, this version reads/writes via
      the fast `dss.get()`/`dss.put()` `TimeSeriesContainer` path
      rather than the `tsmath` wrapper.
    """

    # Retrieve cleaned list of DSS pathnames
    recs = DSS_Tools.get_sanitized_record_list(dss_file)     

    # Open DSS file for read/write operations
    dss = HecDss.open(dss_file)                              
    
    # Loop on each of the records
    for r in recs:
        # Move the record name to lower case.
        rlow = r.lower()
        
        # things not to read: paired data, integer/scalar/text vars and some other things that are causing trouble.
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow and not 'temp-water-target' in rlow:
            # Read as TimeSeriesContainer (fast path)
            tsc = dss.get(r,True)                             

            # Handle the flow records separately
            if "/flow" in rlow or "/1day/" in rlow:
                if not "/storage" and not "/stor" in rlow:    # Guard against storage series
                    tsc.type = 'PER-AVER'                     # Standardize type for period-averaged flows
                    #tsc.setStoreAsDoubles(True)
                    dss.write(tsc)                            # Write back updated type

            # Convert the units to lower case for comparison
            units = str(tsc.units).lower()  
            
            # Handle the other types of values
            if units in units_need_fixing:
                if units == 'tenths':
                    # save off a copy of cloud record in 0-1 for ResSim
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'                   # Annotate E-part with FRAC
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'                        # Standard fractional unit
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 10.0  # Convert tenths → fraction
                    #tsc.setStoreAsDoubles(True)         
                    dss.put(tsc)                              # Write converted series
                
                if units == 'radians':
                    # save off a copy in deg
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-DEG'                    # Tag E-part with DEG
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'deg'                         # Convert to degrees
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0  # rad → deg
                    #tsc.setStoreAsDoubles(True)         
                    dss.put(tsc)
                
                if units == 'deg':
                    # save off a copy in redians
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-RADIANS'                # Tag E-part with RADIANS
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'radians'                     # Convert to radians
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)  # deg → rad
                    #tsc.setStoreAsDoubles(True)            
                    dss.put(tsc)
                
                if units == 'kph':
                    # convert to m/s 
                    tsc.units = 'm/s'                         # Target unit for wind speed
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6   # kph → m/s
                    #tsc.setStoreAsDoubles(True)
                    dss.put(tsc)

    # Close DSS file handle
    dss.close()                                               


def fix_DMS_types_units_old(dss_file):
    """
    Standardize DMS records (legacy implementation).

    This method was implemented to change data types to `PER-AVER`
    for records that are not coming from the DMS that way. It is the
    older counterpart to `fix_DMS_types_units`, using the `tsmath`
    wrapper for reads/transforms rather than reading raw
    `TimeSeriesContainer` objects directly.

    Parameters
    ----------
    dss_file : str
        Path to source DSS file.

    Returns
    -------
    None
        Writes corrected records back to `dss_file`.

    Notes
    -----
    - Uses `dss.read` (math container) to set type and units, then writes
      the underlying `TimeSeriesContainer`.
    - Includes additional `m/s → W2link` copy for wind speed workaround,
      which the newer `fix_DMS_types_units` does not perform.
    """
    
    # Gather sanitized DSS record list
    recs = DSS_Tools.get_sanitized_record_list(dss_file)      

    # Open DSS file for transformations
    dss = HecDss.open(dss_file)                                
    
    # Loop over each DSS record
    for r in recs:
        # Convert the record name to lower case
        rlow = r.lower()
        
        # things not to read: paired data, integer/scalar/text vars and some other things that are causing trouble.
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow:
        
            # Read math container for path
            tsm = dss.read(r)                                  

            # Handle flow records separately
            if "/flow" in rlow or "/1day/" in rlow:
                tsm.setType('PER-AVER')                        # Standardize to period-averaged
                tsc = tsm.getData()
                dss.write(tsc)                                 # Write back updated type
            
            # Handle other types of values
            if tsm.getUnits().lower() in units_need_fixing:
                if tsm.getUnits() == 'tenths':
                    # save off a copy of cloud record in 0-1 for ResSim
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 10.0         
                    dss.write(tsc)
                
                if tsm.getUnits() == 'radians':
                    # save off a copy in deg
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-DEG'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'deg'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0       
                    dss.write(tsc)
                
                if tsm.getUnits() == 'deg':
                    # save off a copy in redians
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-RADIANS'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'radians'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)
                    dss.write(tsc)
                
                if tsm.getUnits() == 'kph':
                    # convert to m/s 
                    tsc = tsm.getData()
                    tsc.units = 'm/s'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    dss.write(tsc)

                if tsm.getUnits() == 'm/s':
                    # make a copy divied by kph conversion as a hack to get W2 linking the wind speed correctly 
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    
                    if not "w2link" in rec_parts[3].lower():
                        rec_parts[3] += '-W2link'              # Annotate E-part for W2 link workaround
                        tsc.fullName = '/'.join(rec_parts)
                        
                        for i in range(len(tsc.values)) :
                            tsc.values[i] = tsc.values[i] / 3.6 # Adjust values for W2 linkage
                        dss.write(tsc)
    
    # Close DSS handle
    dss.close()                                               


def standardize_bc_temp_water_to_C(dss_file,output_dss_file):
    """
    Create Celsius copies of `temp-water` records for ResSim linking.

    Parameters
    ----------
    dss_file : str
        Source DSS file containing `temp-water` series.
    output_dss_file : str
        Destination DSS file to write standardized Celsius copies.

    Returns
    -------
    None
        Writes derivative records with `-C` E-part and `C` units.

    Notes
    -----
    - If source and destination files are the same, writes in-place using
      the same handle.
    - Fahrenheit inputs (`F`/`degF`) are converted to Celsius.
    """
    
    # Collect relevant DSS pathnames
    recs = DSS_Tools.get_sanitized_record_list(dss_file)    

    # Open source DSS
    dss = HecDss.open(dss_file)                               

    # Check where the data is going to be written
    if dss_file == output_dss_file:
        # Reuse handle if writing in-place
        dss_out = dss                                         
    
    else:
        # Separate output DSS
        dss_out = HecDss.open(output_dss_file)                
    
    # Loop and handle each record in the file
    for r in recs:
        # Convert the dss path to lower case for comparison
        rlow = r.lower()
        
        # Filter only temp-water series
        if '/temp-water' in rlow:                             
            # Read series as container
            tsc = dss.get(r,True)                             

            # Track original units for conversion
            incoming_units = tsc.units.lower()                
        
            # Get fresh container for output name edit
            tsc = dss.get(r,True)                             
            rec_parts = tsc.fullName.split('/')
            rec_parts[3] += '-C'                              # E-part annotation for Celsius
            tsc.fullName = '/'.join(rec_parts)
            tsc.units = 'C'                                   # Set output units to Celsius
            
            # Convert units if necessary
            if incoming_units == 'f' or incoming_units == 'degf':                
                for i in range(len(tsc.values)) :
                    tsc.values[i] = (tsc.values[i] - 32.0)*5.0/9.0  # F → C conversion             

            # Write standardized series
            dss_out.put(tsc)                                  

    # Close the DSS file
    dss.close()
    
    # Close output handle if separate
    if dss_file != output_dss_file:
        dss_out.close()                                       


def DMS_fix_units_types(hydro_dss,met_dss_file):
    """
    Convenience wrapper to standardize hydro and met DSS files.

    Parameters
    ----------
    hydro_dss : str
        Path to DMS hydro DSS file.
    met_dss_file : str
        Path to DMS meteorological DSS file.

    Returns
    -------
    None
        Calls `fix_DMS_types_units` on both inputs.
    """
    
    # Sanitize hydro time series
    fix_DMS_types_units(hydro_dss)

    # Sanitize meteorological time series   
    fix_DMS_types_units(met_dss_file)                         


def splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file, months=[1,2,3]):
    """
    Replace Lewiston met data with Redding met data for specified months.

    Parameters
    ----------
    currentAlternative : object
        WAT scripting alternative used for logging.
    rtw : object
        Runtime window (start/end strings, HEC times).
    met_dss_file : str
        Input meteorological DSS file.
    output_dss_file : str
        Destination DSS file for spliced outputs.
    months : list of int, optional
        Month numbers to splice (default `[1, 2, 3]` → Jan–Mar).

    Returns
    -------
    None
        Writes updated met records to `output_dss_file`.

    Notes
    -----
    - Uses `DSS_Tools.replace_data` with a fixed set of Lewiston↔Redding pairs.
    """
    
    # Lewiston is still dependent on using Met data from Redding during Jan-Feb-Mar.  Create those spliced Met data records.
    pairs = [["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Air temperature/Temp-Air//1Hour/232.6.53.1.1/",
              "/MR Sac.-Clear Cr. to Sac R./KRDD-Air temperature/Temp-Air//1Hour/235.40.53.1.1/"],
             ["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Dew Point/Temp-DewPoint//1Hour/232.6.51.1.1/",
              "/MR Sac.-Clear Cr. to Sac R./KRDD-Dew Point/Temp-DewPoint//1Hour/235.40.51.1.1/"],
             ["/MR SAC.-LEWISTON RES./TCAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/232.5.135.1.1/",
              "/MR SAC.-CLEAR CR. TO SAC R./RRAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/235.41.135.1.1/"],
             ["/MR Sac.-Lewiston Res./TCAC1-Wind Direction/Dir-Wind/0/1Hour/232.5.133.1.2/",
              "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Direction/Dir-Wind//1Hour/235.40.133.1.2/"],
             ["/MR Sac.-Lewiston Res./TCAC1-Wind Speed/Speed-Wind//1Hour/232.5.133.1.1/",
              "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Speed/Speed-Wind//1Hour/235.40.133.1.1/"],
             ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover//1Day/236.9.129.1.1/",
              "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover//1Hour/235.41.129.1.1/"],
             ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover-FRAC//1Day/236.9.129.1.1/",
              "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover-FRAC//1Hour/235.41.129.1.1/"]]

    # Execute splice operation
    DSS_Tools.replace_data(currentAlternative, rtw, pairs, met_dss_file, output_dss_file, months, standard_interval='1HOUR')  
   

def compute_river_balance_flows(currentAlternative, rtw, hydro_dss, obs_dss_file, output_dss_file):
    """
    Compute river balance flows (IGO and Bend Bridge) and write results.

    Parameters
    ----------
    currentAlternative : object
        WAT scripting alternative for logging.
    rtw : object
        Runtime window.
    hydro_dss : str
        Hydro DSS file (DMS inputs).
    obs_dss_file : str
        Observational DSS file (USGS records).
    output_dss_file : str
        Destination DSS file for derived river balance flows.

    Returns
    -------
    None
        Writes daily and hourly balance series to `output_dss_file`.

    Notes
    -----
    - Aligns and shifts Bend Bridge records (−1 day) to better balance with releases.
    """
    
    # balance at IGO
    flow_records = [obs_dss_file + "::/USGS CLEAR CR/11372000 NR IGO/FLOW//1HOUR/USGS-MERGED-CROP/",
                    output_dss_file+'::/MR Sac.-Whiskeytown Lake/WHI-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/']
    out_rec = "/CLEAR CR/IGO BALANCE FLOW/FLOW//1HOUR/ResSim_PreProcess/"
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_records, hydro_dss,
                              [True, False],
                              out_rec, output_dss_file, what="flow",prepend_n=25)
    DSS_Tools.resample_dss_ts(output_dss_file,out_rec,rtw,output_dss_file,'1DAY')  # Write daily aggregate

    # balance at Bend Bridge
    DSS_Tools.resample_dss_ts(hydro_dss,'/MR Sac.-Keswick Res./KES-Dam Total Release/Flow//1Hour/234.1.125.1.1/',rtw,output_dss_file,'1DAY',pad_start_days=1)
    DSS_Tools.resample_dss_ts(hydro_dss,'/MR Sac.-Sac River/11377100 Sacramento R at Bend Bridge-15min Flow/Flow//15Minute/237.64.125.1.1/',rtw,output_dss_file,'1DAY',pad_start_days=1)

    # bend bridge is pretty far downstream from the main release at Keswick - shifting backward by one day
    # makes for a better balance
    in_rec = '/MR Sac.-Sac River/11377100 Sacramento R at Bend Bridge-15min Flow/Flow//1Day/237.64.125.1.1/'
    
    # not using out_rec, time shifting is not working
    out_rec = '/MR Sac.-Sac River/11377100 Sacramento R at Bend Bridge-15min Flow/Flow//1Day/237.64.125.1.1-1Day/'
    DSS_Tools.shift_ts_time(output_dss_file,in_rec,output_dss_file,out_rec,'-1Day',start_date=None,end_date=None)    
    DSS_Tools.postprend_last_value_on_ts(output_dss_file,out_rec,1)

    # Construct the flow record paths
    flow_records = [output_dss_file + "::" + in_rec,
                    obs_dss_file + "::/USGS CLEAR CR/11372000 NR IGO/FLOW//1DAY/USGS/",
                    "/MR Sac.-Sac River/11370700 ACID-Dly Flow/Flow//1Day/237.60.125.2.1/",
                    "/MR Sac.-Sac River/11374000 Cow Creek-Dly Flow/Flow//1Day/237.61.125.2.1/",
                    "/MR Sac.-Sac River/11376000 Cottonwood-Dly Flow/Flow//1Day/237.62.125.2.1/",
                    "/MR Sac.-Sac River/11376550 Battle Creek-Dly Flow/Flow//1Day/237.63.125.2.1/",
                    output_dss_file+"::/MR Sac.-Keswick Res./KES-Dam Total Release/Flow//1Day/234.1.125.1.1/",]
    
    # Define the proceprocessing path
    out_rec = "/SACRAMENTO R/BEND BR BALANCE FLOW/FLOW//1DAY/ResSim_PreProcess/"
    
    # Apply the flow adjustment
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_records, hydro_dss,
                              [True, False, False, False, False, False, False],
                              out_rec, output_dss_file, what="flow",prepend_n=2)  # Final Bend Bridge balance


def compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file):
    """
    Create combined outflow records used by ResSim/TCD scripts across the five northern reservoirs.

    Parameters
    ----------
    currentAlternative : object
        WAT scripting alternative for logging.
    rtw : object
        Runtime window.
    hydro_dss : str
        DMS hydro DSS file.
    output_dss_file : str
        Destination DSS file for derived outflows.

    Returns
    -------
    None
        Writes combined outflow series (hourly) for Trinity and Lewiston, then
        delegates Shasta gate combination to `combine_shasta_gates_flows`.
    """
    

    # Trinity: add Generation, G1 &amp; G2, which are powerplant and jet-valve (bypass) flows from the powerplant intake (G3 is the low-level bypass)
    inflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/231.5.125.8.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Trinity Lake/TRN-GenerationG1G2_Sum/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Lewiston: add Generation and outlet flows which come from the same level
    inflow_records = ['/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/232.12.125.2.1/',
                      '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/232.12.125.3.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Lewiston Res./LEW-Gen_plus_Outlet Release/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Shasta gate sums & outlet combinations
    combine_shasta_gates_flows(currentAlternative,rtw,hydro_dss,output_dss_file)  


def compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file):
    """
    Generate hourly total-dam-outflow records for plotting if proper sums do not exist.

    Parameters
    ----------
    currentAlternative : object
        WAT scripting alternative.
    rtw : object
        Runtime window.
    hydro_dss : str
        DMS hydro DSS file.
    output_dss_file : str
        Destination DSS file.

    Returns
    -------
    None
        Writes total outflow records for Whiskeytown, Shasta, and Lewiston.
    """
    
    # Trinity: calc'd total outflow is included in DMS template
    # Keswick: outlet gase only, no summing needed
    # Whiskeytown: outlet + spill = total dam outflow
    inflow_records = ['/MR Sac.-Whiskeytown Lake/WHI-Outlet Release/Flow//1Hour/233.14.125.2.1/',
                      '/MR Sac.-Whiskeytown Lake/WHI-Spill Release/Flow//1Hour/233.14.125.5.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Whiskeytown Lake/WHI-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Shasta: gen + spill + river outlets = total dam release
    inflow_records = ['/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/',
                      '/MR Sac.-Shasta Lake/SHA-Outlet Release/Flow//1Hour/230.11.125.5.1/',
                      '/MR Sac.-Shasta Lake/SHA-Spill Release/Flow//1Hour/230.11.125.4.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Shasta Lake/SHA-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Lewiston: outlet + gen + hatchery + spill = total dam outflow
    inflow_records = ['/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/232.12.125.2.1/',
                      '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/232.12.125.3.1/',
                      '/MR Sac.-Lewiston Res./LEW-Fish Hatchery Release/Flow//1Hour/232.12.125.1.1/',
                      '/MR Sac.-Lewiston Res./LEW-Spill Release Hrly/Flow//1Hour/232.12.125.5.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Lewiston Res./LEW-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)


def combine_shasta_gates_flows(currentAlt,timewindow,hydro_dss,output_dss_file):
    """
    Sum Shasta TCD gate counts by level and combine outlet flows by elevation.

    Parameters
    ----------
    currentAlt : object
        WAT scripting alternative for logging.
    timewindow : object
        Runtime window.
    hydro_dss : str
        DMS hydro DSS file for reads.
    output_dss_file : str
        Destination DSS file for derived sums.

    Returns
    -------
    None
        Writes gate-count sums and outlet-flow sums (750/850/950) to DSS.
    """

    # Construct the gate record DSS path components
    g_rec_1 = "/MR Sac.-Shasta Lake/SHA-TCD-Gate Position "
    g_rec_2 = "/Gate Position "
    g_rec_3 = "//1Hour/230.13.224."
    
    # Loop on each of the gate elevations
    for level in ['Top','Middle','Bottom','Side']:
        # Set the number of gates that are open
        ngate = 2 if level=='Side' else 5
        
        # Create the gate record list
        gate_recs = []
        for ng in range(1,ngate+1):
             sg = str(ng)
             lg = level[0]
             gate_recs.append(g_rec_1 + sg+"_"+lg + g_rec_2 + sg + g_rec_3 + lg+"."+sg)
        
        # Define the output DSS path
        out_rec = "/MR Sac.-Shasta Lake/SHA-TCD-"+level+" Gate Sum/count-gate//1Hour/Derived/"
        
        # Create the flows for each gate elevation
        DSS_Tools.add_or_subtract_flows(currentAlt, timewindow, gate_recs, hydro_dss, 
                       [True for i in range(ngate)],
                       out_rec, output_dss_file, what="n/a")  # Sum counts (treated as “flows” for add function)

    # Combine river outlet flows
    # Example record:
    # "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G1/Flow//1Hour/230.11.125.14.1/"

    # just writing these out - seems a bit confusing to code them like above
    ro_flows = {
    '750': [
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G1/Flow//1Hour/230.11.125.14.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G2/Flow//1Hour/230.11.125.15.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G3/Flow//1Hour/230.11.125.16.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G4/Flow//1Hour/230.11.125.17.1/",
    ],
    '850': [
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G1/Flow//1Hour/230.11.125.18.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G2/Flow//1Hour/230.11.125.19.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G3/Flow//1Hour/230.11.125.20.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G4/Flow//1Hour/230.11.125.21.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G5/Flow//1Hour/230.11.125.22.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G6/Flow//1Hour/230.11.125.23.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G7/Flow//1Hour/230.11.125.24.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G8/Flow//1Hour/230.11.125.25.1/",
    ],
    '950' : [
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G1/Flow//1Hour/230.11.125.26.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G2/Flow//1Hour/230.11.125.27.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G3/Flow//1Hour/230.11.125.28.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G4/Flow//1Hour/230.11.125.29.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G5/Flow//1Hour/230.11.125.30.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G6/Flow//1Hour/230.11.125.31.1/",
    ]}
    
    # Create the output series for each river outlet
    for level,gate_recs in ro_flows.items():
        out_rec = "/MR Sac.-Shasta Lake/SHA-Outlet Flow "+level+" Sum/Flow//1Hour/Derived/"
        DSS_Tools.add_or_subtract_flows(currentAlt, timewindow, gate_recs, hydro_dss,
                              [True for i in range(len(gate_recs))],
                              out_rec, output_dss_file, what="flow")  # Sum outlet flows per elevation band
    
def leakage(gen_flow,nMidGates,nLowGates,nSideGates,WSE,pre2010=False):
    """
    Estimate leakage flow through a dam's gate levels.

    Total leakage is first estimated as a fraction of generation
    flow, then split across six leakage points (LKG1-LKG6) using
    empirically-derived percentage splits that depend on how many
    gates are open at the mid, low, and side gate levels. A
    different set of empirical percentages is used depending on
    whether the dam configuration is pre-2010 or post-2010. Finally,
    leakage at the uppermost points (LKG1, LKG2, LKG3) is reduced or
    zeroed out if the water surface elevation has dropped below
    their withdrawal elevation, since a leakage point cannot leak
    water it cannot reach.

    Parameters
    ----------
    gen_flow : float
        Total generation (powerhouse) flow, used as the basis for
        estimating total leakage flow.
    nMidGates : int
        Number of open gates at the mid gate level. Only used in the
        pre-2010 calculation.
    nLowGates : int
        Number of open gates at the low gate level.
    nSideGates : int
        Number of open gates at the side gate level.
    WSE : float
        Current water surface elevation, in feet, used to determine
        which leakage points are still submerged.
    pre2010 : bool, optional
        If True, use the pre-2010 empirical leakage percentage
        splits. If False (default), use the post-2010 splits.

    Returns
    -------
    list of float
        `[LKG1, LKG2, LKG3, LKG4, LKG5, LKG6]`, the estimated leakage
        flow at each of the six leakage points, in the same units as
        `gen_flow`.

    Notes
    -----
    The elevation bands used to taper/zero the upper leakage points
    (945-1000 ft, 900-945 ft, 831-900 ft) are hardcoded to the
    specific dam being modeled.
    """
    # select the empirical leakage percentage splits based on dam configuration era
    if pre2010:
        # Apply the pre2010 flow ratios
        leakage_flow = gen_flow * 0.2
        LKG1 = (13.09 / 100) * leakage_flow
        LKG2 = ((8.05 + (5 - nMidGates) / 5 * 11.65) / 100) * leakage_flow
        LKG3 = (9.34 + (2 - nSideGates) / 2 * 3.31) / 100 * leakage_flow
        LKG4 = (1.03 + (5 - nLowGates) / 5 * 10.01) / 100 * leakage_flow
        LKG5 = (3.84 + 31.12) / 100 * leakage_flow
        LKG6 = (1.79 + 6.77) / 100 * leakage_flow
    
    else:
        # Apply the post 2010 adjustment
        leakage_flow = gen_flow * 0.1606
        LKG1 = 16.3 / 100 * leakage_flow
        LKG2 = 0.0
        LKG3 = (11.63 + (2 - nSideGates) / 2 * 4.12) / 100 * leakage_flow
        LKG4 = (1.28 + (5 - nLowGates) / 5 * 12.47) / 100 * leakage_flow
        LKG5 = (4.78 + 38.76) / 100 * leakage_flow
        LKG6 = (2.23 + 8.44) / 100 * leakage_flow

    # don't leak if leakage withdraw elevation is out of the water
    # taper or zero out the upper leakage points as WSE drops through their elevation bands
    if 945 <= WSE <1000:
        LKG1 *= (WSE-945)/(1000-945)
    elif 900 <= WSE < 945:
        LKG1 = 0.0
        LKG2 *= (WSE-900)/(945-900)                 # Scale LKG2 by depth fraction
    
    elif 831 <= WSE < 900:
        LKG1 = 0.0
        LKG2 = 0.0
        LKG3 *= (WSE-831)/(900-831)                 # Scale LKG3 by depth fraction

    # Return leakage components
    return [LKG1,LKG2,LKG3,LKG4,LKG5,LKG6]         


def get_w2_withdraw_points(nUpperGates,nMidGates,nLowGates,nSideGates,WSE,withdraw_elevs):
    """
    Determine W2 withdrawal points based on open gates and water surface.

    Parameters
    ----------
    nUpperGates : int
        Count of open upper gates.
    nMidGates : int
        Count of open middle gates.
    nLowGates : int
        Count of open low gates.
    nSideGates : int
        Count of open side gates.
    WSE : float
        Water surface elevation (ft).
    withdraw_elevs : dict
        Map of gate groups (`TCDU`, `TCDM`, `TCDL`, `TCDS`) to elevations.

    Returns
    -------
    list of str
        List of withdrawal point identifiers (e.g., `['TCDU1','TCDU2']`).

    Raises
    ------
    SystemExit
        Calls `sys.exit(-1)` if no qualifying withdrawal points can be
        determined for any gate level, even after the fallback logic.

    Notes
    -----
    - Withdraw points must lie at least 3 ft below WSE.
    - Falls back to highest available level if no open-gate points qualify.
    """

    # Define the gate levels and the number of gates that are open
    gates = ['TCDU','TCDM','TCDL','TCDS']
    n_gates_open = [nUpperGates, nMidGates, nLowGates, nSideGates]

    # Loop and determine which the points are active based on the gates
    w_points = []
    for gate, n in zip(gates,n_gates_open):
        if n > 0 and not w_points:
            for wl,elev in withdraw_elevs[gate].items():
                if elev < WSE - 3: # withdraws must be at least 3 ft below WSE
                    w_points.append(gate + str(wl))

    # Check that the point assignment has produced a valid configuration
    if not w_points:
        # Above failed, usually because gates are erroneously only open above WSE.
        # In this case, pick highest level that allows flow (without considering if gates are open or not)
        for gate in gates:
            if not w_points:
                for wl, elev in withdraw_elevs[gate].items():
                    if elev < WSE - 3: # withdraws must be at least 3 ft below WSE
                        w_points.append(gate + str(wl))

    # Recheck that the above correction fixed the issue. If not, error for safe handling.
    if not w_points:
        # there is some bad data
        print("get_w2_withdraw_points fail!")
        print("Gates:",nUpperGates,nMidGates,nLowGates,nSideGates)
        print("WSE:",WSE)
        sys.exit(-1)                                        # Halt if no valid withdrawal points found

    # Return candidate withdraw points
    return w_points                                         

def repeat_annual_daily_over_timewindow(timewindow,annual_dss_file,annual_dss_rec,output_dss_file,out_rec):
    """
    Repeat a single 365-value annual daily pattern across every year
    spanned by a run time window.

    For each year in the window's range, the same 365 daily values
    are copied in as-is, with timestamps generated for that specific
    year. For leap years, an extra value is appended at the end of
    that year's block (a repeat of the year's last value) so that
    366 daily timestamps are produced instead of 365.

    Parameters
    ----------
    timewindow : object
        The run time window object, used to determine the starting
        and ending calendar years to cover.
    annual_dss_file : str
        Path to the DSS file containing the source annual daily
        pattern record.
    annual_dss_rec : str
        DSS path of the source record. Must contain exactly 365
        daily values representing one calendar year's pattern.
    output_dss_file : str
        Path to the DSS file to write the repeated multi-year series
        to.
    out_rec : str
        DSS path to write the resulting record to.

    Returns
    -------
    None
        Writes a daily time series spanning from the start year
        through the end year of `timewindow` (inclusive) to `out_rec`
        in `output_dss_file`.
    """
    start_year = HecTime(timewindow.getStartTimeString()).year()
    end_year = HecTime(timewindow.getEndTimeString()).year()
    annual_tsc = DSS_Tools.dss_read_ts_safe(annual_dss_file,annual_dss_rec) # must be 365 values
    values = [annual_tsc.values[i] for i in range(annual_tsc.numberValues)] # convert to list, dumb

    # Create the minutes to day conversion factor
    day_interval = 1440                                           # Minutes per day for DSS interval

    # Loop and repeat the series across years
    values_out = []
    times_out = []
    # for each year spanned by the time window, repeat the annual pattern with that year's timestamps
    for i,y in enumerate(range(start_year,end_year+1)):
        values_out += values                                      # Append 365 daily values
        ht = HecTime("01Jan%i"%y,"2400").value()                  # HEC time for end of Jan 1
        ytimes = [ht+i*day_interval for i in range(len(values))]  # Build daily HEC times
        times_out += ytimes        
        if HecTime.isLeap(y):
            values_out.append(values_out[-1])                     # Repeat last value on leap day
            times_out.append(times_out[-1]+day_interval)          # Append leap day time

    # build and write the combined multi-year daily time series
    tsc_result = TimeSeriesContainer()
    tsc_result.fullName = out_rec
    tsc_result.units = annual_tsc.units
    tsc_result.type = annual_tsc.type
    tsc_result.interval = day_interval
    tsc_result.times = times_out
    tsc_result.values = values_out
    tsc_result.startTime = times_out[0]
    tsc_result.numberValues = len(values_out)

    # Put the record into the container file
    dss_outfile = HecDss.open(output_dss_file)
    dss_outfile.put(tsc_result)                                   # Write repeated daily series
    dss_outfile.close()


def W2_shasta_TCD_flow(timewindow,hydro_dss,output_dss_file):
    """
    Compute CE-QUAL-W2 TCD withdrawal flows for Shasta based on
    gates and WSE.

    Parameters
    ----------
    timewindow : object
        Runtime window for all reads/writes.
    hydro_dss : str
        DMS hydro DSS file (elevation & generation inputs).
    output_dss_file : str
        Destination DSS file for W2 TCD point flow outputs.

    Returns
    -------
    None
        Writes hourly flows for each W2 sink (`TCDU/M/L/S1–3`, leakage `LKG1–6`,
        and `TCD_down`) under `/SHA-W2-TCD-*/Flow//1Hour/Derived/`.

    Raises
    ------
    SystemExit
        Calls `sys.exit(-1)` if the source elevation record's units
        are neither feet nor meters.

    Notes
    -----
    Algorithm, per timestep:
      1. Assign all generation flow initially to the highest open
         gate level.
      2. Determine active withdrawal points via
         `get_w2_withdraw_points` and leakage components via
         `leakage`.
      3. Route 35% of the (non-leaked) flow to `TCD_down` if any side
         gates are open, then split the remainder evenly across the
         active withdrawal points.
      4. Perform a mass-balance sanity check, printing diagnostics if
         the sum of computed sink flows does not match `gen_flow`
         (rounded to the nearest integer).
    """
    
    # W2 TCD point sink elevations - hard code those here
    withdraw_elevs = {'TCDU':{1:1042.0,
                          2:1021.0,
                          3:1000.0},
                  'TCDM':{1:942.0,
                          2:921.0,
                          3:900.0},
                  'TCDL':{1:830.0,
                          2:816.0,
                          3:802.0},
                  'TCDS':{1:800.0,
                          2:760.0,
                          3:720.0}
                   }

    # load TCD gates summary from DSS, make regular hourly if needed and fill
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    tsc_u_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Top Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)
    tsc_m_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Middle Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)
    tsc_l_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Bottom Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)
    tsc_s_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Side Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)

    # get year of all records to determine what leakage to use 
    hectimes = DSS_Tools.hectimes_from_tsc(tsc_u_gates)
    years = [hectimes[i].year() for i in range(len(hectimes))]    # Build year vector from HEC times

    # get WSE data, convert to ft if neccessary
    tsc_wse = DSS_Tools.dss_read_ts_safe(hydro_dss,'/MR Sac.-Shasta Lake/SHA-Elevation/Elev//1Hour/230.11.145.2.1/',starttime_str,endtime_str)
    if tsc_wse.units.lower() != 'ft':
        if tsc_wse.units.lower() == 'm':
            for i in range(tsc_wse.numberValues):
               tsc_wse.values[i] = tsc_wse.values[i] * 3.28084 # so lame we have to do this
            tsc_wse.units = 'ft'
        else:
            print('W2_shasta_TCD_flow: elevation units not understood:',tsc_wse.units)
            sys.exit(-1)                                       # Abort if elevation units unknown

    # get generation flow data, also need hourly, only make calcs where we
    # have both gate and flow data and WSE data
    tsc_gen_flow = DSS_Tools.dss_read_ts_safe(hydro_dss,'/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/',starttime_str,endtime_str)

    # move data to lists
    WSE = tsc_wse.values
    nUpperGates = tsc_u_gates.values
    nMidGates = tsc_m_gates.values
    nLowerGates = tsc_l_gates.values
    nSideGates = tsc_s_gates.values
    gen_flow = tsc_gen_flow.values

    # setup results dict
    withdraws_in_order = ['TCDU1', 'TCDU2', 'TCDU3', 'TCDM1', 'TCDM2', 'TCDM3', 'TCDL1', 'TCDL2', 'TCDL3', 'TCDS1', 'TCDS2', 'TCDS3',
                          'LKG1', 'LKG2', 'LKG3', 'LKG4', 'LKG5', 'LKG6', 'TCD_down']
    
    # Define the holder for the TCD flows
    TCD_flows = {}
    for w in withdraws_in_order:
        TCD_flows[w] = []                                      # Initialize per-sink series list

    # init flow at each of 4 levels to 0
    # 1) assign all flow to highest gate level open
    # 2) divide flow at each level by 3 and assign as the point source
    #    sinks for each level. At this point, all flow is still on the highest
    #    open level
    # 3) apply withdraw depth logic (≥3 ft below WSE) and leakage components
    for i in range(len(gen_flow)):
        # Get the withdrawal points associated with each elevation
        wp = get_w2_withdraw_points(nUpperGates[i],nMidGates[i],nLowerGates[i],nSideGates[i],WSE[i],withdraw_elevs)
        
        # Estimate the leakage that is occuring from each elevation
        leak = leakage(gen_flow[i],nMidGates[i],nLowerGates[i],nSideGates[i],WSE[i],pre2010=years[i] < 2010)

        # Partition flows among the available withdrawal points
        step_flows = {w:0.0 for w in withdraws_in_order}       # Initialize step flow dict
        sum_leak = sum(leak)                                   # Total leakage at this step
        step_flows['TCD_down'] = 0.35 * (gen_flow[i]-sum_leak) if nSideGates[i] > 0 else 0.0
        wp_flow = (gen_flow[i] - sum_leak - step_flows['TCD_down']) / len(wp)  # Distribute among withdraws
        
        # Assign equal share per withdrawal point
        for p in wp:
            step_flows[p] = wp_flow                            

        # Populate leakage sinks
        for j in range(6): # six leakage withddraw points
            step_flows['LKG'+str(j+1)] = leak[j]              

        # Calculate the total flows
        total_flow = 0.0
        for p in withdraws_in_order:
            total_flow += step_flows[p]
            TCD_flows[p].append(step_flows[p])                # Append to time series

        # Sanity check on mass balance
        if int(round(total_flow)) != int(round(gen_flow[i])): 
            print('DEBUG W2 TCD FLOWS-----------------------------------------------------')
            print('WSE: ', WSE[i])
            print(nUpperGates[i],nMidGates[i],nLowerGates[i],nSideGates[i])
            print(wp)
            print(i,total_flow,gen_flow[i])
            print(step_flows)

    # Write the flows to the DSS files
    dss_outfile = HecDss.open(output_dss_file)
    for p,flow in TCD_flows.items():
        # re-use gen_flow tsc to write TCD flows (including original units)
        tsc_gen_flow.fullName = "/MR Sac.-Shasta Lake/SHA-W2-TCD-"+p+"/Flow//1Hour/Derived/"
        tsc_gen_flow.values = flow
        dss_outfile.put(tsc_gen_flow)                          # Write each sink flow series
    
    # Close the DSS file
    dss_outfile.close()


def preprocess_W2_5Res(currentAlternative, computeOptions):
    """
    Preprocess DMS hydrology and meteorology data for the CE-QUAL-W2
    5-reservoir Sacramento/Trinity model.

    Fixes data types/units, creates constant reference records,
    splices in Redding met data for Lewiston, computes reservoir
    outflows and TCD flow, enforces minimum flows to avoid
    zero-flow issues in W2, and corrects a known Pit River timestamp
    offset.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()` and `getTimeStep()` for logging and
        retrieving the balance period.
    computeOptions : object
        The compute options/settings object. Must support
        `getRunDirectory()` to locate the current run.

    Returns
    -------
    bool
        True once preprocessing has completed. Writes multiple
        preprocessed and derived records to the shared
        `DMS_SacTrn_ResSim_Pre-Process.dss` file (and related shared
        DSS files).

    Notes
    -----
    Order-dependent steps: `splice_lewiston_met_data` and
    `compute_5Res_outflows` must run before `W2_shasta_TCD_flow`,
    since the latter depends on outputs of the former.
    """
    rtw = computeOptions.getRunTimeWindow()
    
    # Runtime window for all operations
    rtw = computeOptions.getRunTimeWindow()                    
    
    # Setup the paths for the operations
    run_dir = computeOptions.getRunDirectory()                       # Run directory path (diagnostic)
    project_dir = Project.getCurrentProject().getProjectDirectory()  # Project directory root
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()          # Alternative timestep string (diagnostic)
    shared_dir = os.path.join(project_dir, 'shared')           # Shared assets folder

    output_dss_file = os.path.join(shared_dir,'DMS_SacTrn_ResSim_Pre-Process.dss')  # Preprocess outputs

    # --- Fix source data types/units before any further processing ---
    currentAlternative.addComputeMessage('Rectifying units - this may take a while if the length of DMS data is large...')
    hydro_dss = os.path.join(shared_dir, 'DMS_SacTrnHydroTS.dss')
    fix_DMS_types_units(hydro_dss)                             # Sanitize hydro DSS
    met_dss_file = os.path.join(shared_dir,'DMS_SacTrnMet.dss')
    fix_DMS_types_units(met_dss_file)                          # Sanitize meteorological DSS

    # --- Create constant-value reference records used as fixed boundary conditions ---
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')

    # --- Splice, compute outflows, and compute TCD flow (order-dependent) ---
    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file, months=[1,2,3,12])
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)
    W2_shasta_TCD_flow(rtw,hydro_dss,output_dss_file) # depends on compute_5Res_outflows having completed

    # link these flows for W2, to avoid zero-dam-flow situations - needs to be above 2.0 cfs for the flowweightaverage script to recognize it :/
    # --- Enforce a minimum flow floor to avoid zero-flow issues downstream ---
    DSS_Tools.min_ts(output_dss_file, '/MR Sac.-Trinity Lake/TRN-GenerationG1G2_Sum/Flow//1Hour/ResSim_PreProcess/', 1.1, output_dss_file, 'min_flow')
    DSS_Tools.min_ts(output_dss_file, '/MR Sac.-Lewiston Res./LEW-Gen_plus_Outlet Release/Flow//1Hour/ResSim_PreProcess/', 1.1, output_dss_file, 'min_flow')

    # pit river is recorded at 12:00 daily, W2 plugin can't handle that, shift backward 12 hours
    DSS_Tools.shift_pit_river_time(hydro_dss,"/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Day/230.6.125.1.1/",
                         output_dss_file,"/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Day/Shifted-12H/",
                         start_date=None,end_date=None)

    # Return success flag for WAT scripting alternative
    return True


def preprocess_ResSim_5Res(currentAlternative, computeOptions):
    """
    Preprocess DMS hydrology and meteorology data for the ResSim
    5-reservoir Sacramento/Trinity model.

    Fixes data types/units, standardizes water temperature boundary
    conditions to Celsius, creates constant reference records,
    generates a repeating annual temperature pattern, splices in
    Redding met data for Lewiston, computes reservoir
    outflows/plotting records/river balance flows, and derives an
    air temperature lapse and several relative humidity / dew point
    records.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()` and `getTimeStep()` for logging and
        retrieving the balance period.
    computeOptions : object
        The compute options/settings object. Must support
        `getDssFilename()` and `getRunDirectory()`.

    Returns
    -------
    bool
        True once preprocessing has completed. Writes multiple
        preprocessed and derived records to the shared
        `DMS_SacTrn_ResSim_Pre-Process.dss` file (and related shared
        DSS files).

    Notes
    -----
    Order-dependent steps: met splicing, outflow computation, and
    plotting records must run before `compute_river_balance_flows`,
    since the latter depends on the Whiskeytown total dam flow record
    created in `compute_plotting_records`.
    """
    rtw = computeOptions.getRunTimeWindow()
    
    # Setup the paths for the script
    run_dir = computeOptions.getRunDirectory()                 # Run directory path
    project_dir = Project.getCurrentProject().getProjectDirectory()  # Project directory
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()          # Alternative timestep string
    shared_dir = os.path.join(project_dir, 'shared')           # Shared folder path

    output_dss_file = os.path.join(shared_dir,'DMS_SacTrn_ResSim_Pre-Process.dss')  # Preprocess outputs

    # --- Fix source data types/units, and standardize water temperature units ---
    currentAlternative.addComputeMessage('Rectifying units - this may take a while if the length of DMS data is large...')
    hydro_dss = os.path.join(shared_dir, 'DMS_SacTrnHydroTS.dss')
    fix_DMS_types_units(hydro_dss)                             # Standardize hydro
    met_dss_file = os.path.join(shared_dir,'DMS_SacTrnMet.dss')
    fix_DMS_types_units(met_dss_file)                          # Standardize met
    
    # ressim can't handle different units under model linking
    standardize_bc_temp_water_to_C(hydro_dss,output_dss_file)  # Ensure temp-water in Celsius

    # --- Create constant-value reference records used as fixed boundary conditions ---
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=1, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ONES',fpart='ONES')

    # --- Generate a repeating annual temperature pattern for Spring Creek Debris Dam ---
    repeat_annual_daily_over_timewindow(rtw,os.path.join(shared_dir,'PatternHydro.dss'),
                                        "//SPRING_CREEK_DEBRIS_DAM/TEMP-WATER//1Day/ANNUAL_TEMPLATE/",
                                        output_dss_file,"//SPRING_CREEK_DEBRIS_DAM/TEMP-WATER//1Day/ANNUAL/")
                        
    # if template IDs exist still, remove them
    #DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    #DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)

    # --- Splice met data, then compute outflows, plotting records, and balance flows (order-dependent) ---
    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file,months=[1,2,3])
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)    
    compute_river_balance_flows(currentAlternative, rtw, hydro_dss, 
        os.path.join(shared_dir,"WTMP_SacTrn_Historical.dss"), output_dss_file)  # depends on WHI dam flow, created in compute_plotting_records

    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    # --- Apply air temperature lapse correction for Shasta Lake's elevation ---
    currentAlternative.addComputeMessage('lapse infile: '+met_dss_file)
    currentAlternative.addComputeMessage('lapse outfile: '+output_dss_file)
    
    # Estimate the meteorological properties
    DSS_Tools.airtemp_lapse(met_dss_file, "/MR SAC.-CLEAR CR. TO SAC R./KRDD-AIR TEMPERATURE/TEMP-AIR//1HOUR/235.40.53.1.1/",
                  -0.7, output_dss_file, "Shasta_Lapse")                                         # Lapse at Shasta elevation

    DSS_Tools.relhum_from_at_dp(met_dss_file,
                      "/MR Sac.-Clear Cr. to Sac R./KRDD-Air temperature/Temp-Air//1Hour/235.40.53.1.1/",
                      "/MR Sac.-Clear Cr. to Sac R./KRDD-Dew Point/Temp-DewPoint//1Hour/235.40.51.1.1/")
    DSS_Tools.relhum_from_at_dp(met_dss_file,
                      "/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Air temperature/Temp-Air//1Hour/232.6.53.1.1/",
                      "/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Dew Point/Temp-DewPoint//1Hour/232.6.51.1.1/")
    DSS_Tools.dp_from_at_relhum(met_dss_file,
                      "/MR SAC.-CLEAR CR. TO SAC R./KRDD-AIR TEMPERATURE/TEMP-AIR//1Hour/Shasta_Lapse/",
                      "/MR Sac.-Clear Cr. to Sac R./KRDD-/RELHUM-FROM-AT-DP//1Hour/235.40.53.1.1-DERIVED/")  # Dewpoint from lapse + RH
   
    # Signal success to WAT
    return True
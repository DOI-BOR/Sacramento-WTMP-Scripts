
"""
OutputLink_W2-Whiskeytown-Downstream_SacTrn
===========================================

Compose downstream output link datasets for Whiskeytown within the Sacramento
WTMP workflow, using Jython and HEC-DSS APIs. This script performs two primary
operations during a WAT Scripting Alternative compute:

1. **Flow Combination (Whiskeytown downstream link)**  
   Reads a base time series and adds one additional time series (first two input
   locations only), converting units (CMS → CFS) as needed, and writes the
   combined time series to the model's DSS file. A log message reports the
   number of written values.

2. **Fixed Heating for Spring Creek Tunnel**  
   Applies a constant temperature increase to the Spring Creek tunnel time
   series and writes the heated series to DSS. It also writes a duplicate using
   the input model F-part to support plotting with the input model.

The script expects to be called by HEC-WAT during the compute of a
ScriptingAlternative, and uses runtime window and DSS path information provided
via WAT/HEC libraries.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `hec.heclib.dss`, `hec.io`, `DSS_Tools`, `rma.util.RMAConst`,
  `com.rma.model`. These are part of HEC/USBR modeling toolchains.
- **Unit Handling:** If an input time series is in `cms`, it is converted to
  `cfs` using the constant `35.314666212661`; otherwise it is treated as-is.
- **Error Handling:** The script logs messages to `currentAlternative` and uses
  `sys.exit(1)` in cases where input data locations are missing or insufficient.
- **F-part Fixing:** The helper `fixFpartToInput` replaces the F-part of an
  output path with the F-part from the input path to align series for plotting.

See Also
--------
DSS_Tools
    Utility module (reloaded at import) that provides helpers such as
    `getDataLocationDSSInfo` and `fixInputLocationFpart` used for path
    resolution and F-part alignment.
"""

import sys  # Standard library: process/system utilities; used here for sys.exit on invalid inputs

from hec.heclib.dss import HecDss  # HEC-DSS API: open/read/write DSS files and time series
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE  # HEC utility constant (unused in this script but available)
from hec.io import DSSIdentifier  # HEC I/O: DSS path identifier helper (not directly used in this script)
from hec.io import TimeSeriesContainer  # HEC I/O: container class for time series data written/read to/from DSS

from rma.util.RMAConst import MISSING_DOUBLE  # RMA/HEC constant representing missing numeric values in series
from com.rma.model import Project  # RMA model package (not directly used here; available in WAT environment)

import math  # Standard library math utilities (not used directly in this script)
import datetime as dt  # Standard library datetime; alias dt (not used directly in this script)

import DSS_Tools  # Local/WTMP utility module: DSS path resolution, F-part alignment, etc.
reload(DSS_Tools)  # Ensure the latest version of DSS_Tools is loaded into the Jython runtime


def fixFpartToInput(inputpath, outpath):
    """
    Replace the F-part of a DSS output path with the F-part taken
    from a given input DSS path.

    This keeps an output record tagged with the same F-part
    (typically a model/run identifier) as the input it was derived
    from, rather than whatever F-part was assigned when the output
    path was created.

    Args:
        inputpath (str): DSS path string for the input location.
            Its F-part (index 6 after splitting on '/') is used.
        outpath (str): The DSS output path whose F-part should be
            replaced.

    Returns:
        str: The output DSS path with its F-part replaced by the
            F-part from inputpath.
    """
    # get F-part from input locations
    location_fpart = inputpath.split('/')[6]  # Extract F-part from the input path
    out_parts = outpath.split('/')            # Split the output path into parts

    out_parts[6] = location_fpart             # Replace output F-part with input F-part
    return '/'.join(out_parts)                # Recompose the full DSS path


def computeAlternative(currentAlternative, computeOptions):
    """
    Compute a scripting alternative that performs two independent
    operations for the Sacramento River system:

    1. Keswick flow balancing: sums the first two configured input
       locations (converting units to cfs as needed) and writes the
       combined flow to the first output location.
    2. Clear Creek tunnel heating: adds a fixed temperature offset
       to the Spring Creek (last input location) time series and
       writes the heated result to the second output location, both
       under its own F-part and under the input model's F-part (for
       easier plotting alongside the source model).

    Workflow:
      1. Logs a compute status message to the alternative.
      2. Validates that at least 2 input data locations are
         configured; exits if fewer than 2 are found.
      3. Reads the first input location as the base time series.
      4. Determines the base series' units and computes a
         cms-to-cfs conversion factor if needed.
      5. For the second input location, reads its time series,
         converts units if needed, and adds it value-by-value to
         the base series, substituting MISSING_DOUBLE for any
         index where a value is unavailable.
      6. Writes the combined flow time series to the first output
         location in cfs.
      7. Re-reads the input data locations and selects the last one
         (Spring Creek) as the source for the tunnel heating step.
      8. Adds a fixed heating offset to every value in that series.
      9. Writes the heated series to the second output location,
         then writes a second copy under the F-part of the original
         input model, so both the model output and a
         plotting-friendly, F-part-matched copy exist in DSS.

    Args:
        currentAlternative: The alternative object being computed.
            Must support addComputeMessage(), getInputDataLocations(),
            getOutputDataLocations(), createOutputTimeSeries(), and
            loadTimeSeries() for logging and resolving linked data
            locations.
        computeOptions: The compute options/settings object. Must
            support getDssFilename() and getRunTimeWindow() to
            provide the target DSS file and the time window to
            compute over.

    Returns:
        bool: True once both the flow balancing and tunnel heating
            steps have completed and their results have been
            written to DSS. Calls sys.exit(1) instead of returning
            if fewer than 2 input data locations are found.
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    
    # Retrieve the list of input data locations provided by the ScriptingAlternative
    locations = currentAlternative.getInputDataLocations()
    # validate that enough input data locations are configured
    if len(locations) == 1:
        currentAlternative.addComputeMessage("Found only 1 datapath locations. Need at least 2.")
        sys.exit(1)
    
    elif len(locations) == 0:
        currentAlternative.addComputeMessage("Found no datapath locations. Exiting.")
        sys.exit(1)
    
    # Base time series: use the first location as the "base" to which others are added
    base_ts = locations[0]  # this is just the first, we'll add everything to this one 
    
    # Resolve DSS path information (tspath and file path) for the base input
    base_tspath, base_DSSPath = DSS_Tools.getDataLocationDSSInfo(base_ts, currentAlternative, computeOptions)   
    currentAlternative.addComputeMessage('Found DSS path: {0}'.format(base_tspath))
    
    currentAlternative.addComputeMessage('\n')  # Visual break in logs
    rtw = computeOptions.getRunTimeWindow()     # Obtain runtime window object for start/end times

    # Adding flow to keep keswick balanced, coming from independent DSS files, so can't
    # use DSS_Tools (currently).  Only add first two input locations!!
    # ------------------------------------------------------------------------------------ 

    # Prepare output: use the first output location to write the combined time series
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])    
    
    # if extra output locations were configured, note that only the first is used
    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))

    # Align F-part in output path to match the input model conventions via DSS_Tools
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    # ---- Read the base time series over the runtime window --------------------------------
    currentAlternative.addComputeMessage("\n##### Adding Timeseries #####")
    starttime_str = rtw.getStartTimeString()  # Runtime window start as string
    endtime_str = rtw.getEndTimeString()      # Runtime window end as string
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    # read the base series' data over the run time window, then close the source file
    dss_base = HecDss.open(base_DSSPath)
    base_TS = dss_base.read(base_tspath, starttime_str, endtime_str, False)
    dss_base.close()                          # Close DSS handle to avoid leaks

    base_TS = base_TS.getData()               # Extract the data container (TimeSeriesContainer-like)
    hecstarttimes = base_TS.times             # HEC times for the time series
    all_values = base_TS.values               # Values to be combined (will be updated)
    units = base_TS.units                     # Units of the base series
    dsstype = base_TS.type                    # Data type (e.g., INST-VAL)

    # we rely here on the input from W2 Whisketytown and the debris dam flow from the DMS, which might be in
    # different units

    cms2cfs_base = 1.0
    # if the base series is in cms, convert it to cfs
    if str(units).lower() == 'cms':
        cms2cfs_base = 35.314666212661        # Convert base series from CMS to CFS if needed

    # for the second input location, read its series and add it to the running total
    for other_loc in locations[1:2]:
        new_values = []  # Accumulator for combined values

        # Resolve DSS path info for the additional location
        other_tspath, other_DSSPath = DSS_Tools.getDataLocationDSSInfo(other_loc, currentAlternative, computeOptions)
        currentAlternative.addComputeMessage("Adding {0}".format(other_tspath))
        currentAlternative.addComputeMessage('DSS path: {0}'.format(other_DSSPath))

        # Read the additional series over the same runtime window
        dss_other = HecDss.open(other_DSSPath)
        other_TS = dss_other.read(other_tspath, starttime_str, endtime_str, False)
        dss_other.close()
        other_TS = other_TS.getData()  # Extract data container for values/units

        # Determine conversion factor for the additional series
        cms2cfs_other = 1.0
        # if the other series is in cms, convert it to cfs
        if str(other_TS.units).lower() == 'cms':
            cms2cfs_other = 35.314666212661
        
        # Combine values element-wise (with unit conversion), handling missing entries
        other_values = other_TS.values
        
        # for each value in the base series, add the corresponding converted value from the other series
        for vi, val in enumerate(all_values):
            try:
                new_values.append(val*cms2cfs_base + other_values[vi]*cms2cfs_other)
            except:
                # If the other series is missing at index vi, log and mark as missing
                currentAlternative.addComputeMessage("No value for location {0} at idx {1} {2}".format(other_loc, vi, hecstarttimes[vi]))
                new_values.append(MISSING_DOUBLE)

        # Update the running "all_values" with the newly combined series
        all_values = new_values
        
    # build a new time series container holding the combined (balanced) flow values, in cfs
    dssfn = computeOptions.getDssFilename()
    dssFm = HecDss.open(dssfn)
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputpath
    tsc.values = all_values
    tsc.startTime = hecstarttimes[0]
    tsc.units = 'cfs'
    tsc.type = dsstype
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(all_values)
    tsc.startHecTime = rtw.getStartTime()
    tsc.endHecTime = rtw.getEndTime()
    
    # write the combined flow time series to the first output location and close the file
    dssFm.write(tsc)
    dssFm.close()
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(new_values)))

    
    # Add tunnel heating to Clear Creek tunnel temperatures
    # ------------------------------------------------------------------------------------ 

    # Refresh input locations (Spring Creek tunnel expected as the last location)
    locations = currentAlternative.getInputDataLocations()   
    SpringCreek = locations[-1]               # Take the final input location as Spring Creek
    
    # resolve the Spring Creek input location's DSS path, correcting its F-part
    SpringCreek = locations[-1]
   
    tspath =str(currentAlternative.loadTimeSeries(SpringCreek))
    tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)
            
    currentAlternative.addComputeMessage('Found DSS path: {0}'.format(tspath))
    
    currentAlternative.addComputeMessage('\n')  # Break in the log
    dssFile = computeOptions.getDssFilename()   # DSS file used by the model
    rtw = computeOptions.getRunTimeWindow()     # Runtime window (reused)

    # second outpath is CC tuneel w/ heating
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[1])
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    # ---- Read Spring Creek series and apply fixed heating --------------------------------
    currentAlternative.addComputeMessage("\n##### PERFORMING HEATING #####")
    starttime_str = rtw.getStartTimeString()  # Get start time string
    endtime_str = rtw.getEndTimeString()      # Get end time string
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    # read the Spring Creek series over the run time window, to be heated below
    dssFm = HecDss.open(dssFile)
    TS = dssFm.read(tspath, starttime_str, endtime_str, False)
    TS = TS.getData()                         # Extract data container
    hecstarttimes = TS.times                  # Time vector for Spring Creek

    heat_amount_c = 0.32
    new_values = []
    
    # for each value in the Spring Creek series, add the fixed heating offset
    for val in TS.values:
        new_values.append(val + heat_amount_c)  # Apply constant temperature increase
        
    # build a new time series container holding the heated Spring Creek values
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputpath
    tsc.values = new_values
    tsc.startTime = hecstarttimes[0]
    tsc.units = TS.units
    tsc.type = TS.type
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(new_values)
    tsc.startHecTime = rtw.getStartTime()
    tsc.endHecTime = rtw.getEndTime()
    # write the heated series to the second output location
    dssFm.write(tsc)

    dssFm.write(tsc)  # Write the heated Spring Creek series

    # Write a copy out using the input model f-part, so that plotting will be able to use it with the input model
    tsc.fullName = fixFpartToInput(tspath, outputpath)  # Mirror input F-part for plot compatibility
    dssFm.write(tsc)    
    
    dssFm.close()  # Close DSS file after writes
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(new_values)))

    # Exit
    return True

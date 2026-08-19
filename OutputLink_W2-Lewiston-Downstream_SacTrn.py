"""
OutputLink_W2-Lewiston-Downstream_SacTrn
=========================================

Compute a WAT Scripting Alternative that performs two operations for the
Lewiston/Clear Creek system:

1. Flow-weighted averaging of Lewiston outflow temperatures for linking to
   downstream models, written both under a cleaned-up F-part and under the
   F-part of the first input location (for plotting alongside that model).
2. Clear Creek tunnel heating: adds a month-specific temperature offset to
   the Clear Creek tunnel time series and writes the heated result, both
   under its own F-part and under the input model's F-part.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `hec.heclib.dss`, `hec.io`, `rma.util.RMAConst`,
  `BoundaryFixes`, `flowweightaverage`, `DSS_Tools`.
"""

import sys                                            # Standard library: used here only to print sys.path for diagnostics
print(sys.path)                                       # Diagnostic: show module search path at script load time
import BoundaryFixes                                  # Local module: threshold/NaN replacement helpers (imported but not directly called in this file)
reload(BoundaryFixes)                                 # Jython: ensure latest version is loaded
from hec.heclib.dss import HecDss                     # HEC-DSS: open/read/write DSS files
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE   # HEC sentinel for undefined values (imported for completeness)
from hec.io import DSSIdentifier                      # DSS path identifier helper (not directly used)
from hec.io import TimeSeriesContainer                # Container for writing time series to DSS
from rma.util.RMAConst import MISSING_DOUBLE          # RMA sentinel for missing values (imported for completeness)
import math                                           # Standard math library (not directly used in this file)
import datetime as dt                                 # Python datetime, used for month-based heating lookups
import flowweightaverage                              # Local module: flow-weighted average temperature computation (FWA2)
reload(flowweightaverage)                             # Reload to ensure latest version
import DSS_Tools                                      # Local module: DSS path resolution / F-part fixing helpers
reload(DSS_Tools)                                     # Reload to ensure latest version

##
#
# computeAlternative function is called when the ScriptingAlternative is computed.
# Arguments:
#   currentAlternative - the ScriptingAlternative. hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt
#   computeOptions     - the compute options.  hec.wat.model.ComputeOptions
#
# return True if the script was successful, False if not.
# no explicit return will be treated as a successful return
#
##

def fixFpartToInput(inputpath, outpath):
    """
    Replace the F-part of a DSS output path with the F-part taken
    from a given input DSS path.

    This keeps an output record tagged with the same F-part
    (typically a model/run identifier) as the input it was derived
    from, rather than whatever F-part was assigned when the output
    path was created.

    Parameters
    ----------
    inputpath : str
        DSS path string for the input location. Its F-part (index 6
        after splitting on `'/'`) is used.
    outpath : str
        The DSS output path whose F-part should be replaced.

    Returns
    -------
    str
        The output DSS path with its F-part replaced by the F-part
        from `inputpath`.
    """
    # get F-part from input locations
    location_fpart = inputpath.split('/')[6]
    out_parts = outpath.split('/')
    out_parts[6] = location_fpart
    return '/'.join(out_parts)

def computeAlternative(currentAlternative, computeOptions):
    """
    Compute a scripting alternative for the Lewiston/Clear Creek
    system.

    Performs two operations:

    1. Flow-weighted averaging of Lewiston outflow temperatures
       (excluding the last input location) for linking, writing the
       result both under a cleaned-up F-part and under the F-part
       of the first input location (for plotting alongside that
       model).
    2. Clear Creek tunnel heating: adds a month-specific temperature
       offset to the last input location's time series and writes
       the heated result to the second output location, both under
       its own F-part and under the input model's F-part.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()`, `getInputDataLocations()`,
        `getOutputDataLocations()`, `createOutputTimeSeries()`, and
        `loadTimeSeries()` for logging and resolving linked data
        locations.
    computeOptions : object
        The compute options/settings object. Must support
        `getDssFilename()` and `getRunTimeWindow()` to provide the
        target DSS file and the time window to compute over.

    Returns
    -------
    bool
        True once the flow-weighted average and tunnel heating
        steps have completed and their results have been written to
        DSS.

    Notes
    -----
    Workflow:

    1. Logs a compute status message to the alternative.
    2. Retrieves all input data locations except the last one, and
       organizes them into record pairs via
       `flowweightaverage.organizeLocations`.
    3. Logs the resolved DSS paths for each organized location pair.
    4. Retrieves the DSS filename and run time window.
    5. Resolves the first output location's path and strips any
       prefix before a `'|'` character out of its F-part, to get a
       clean F-part for the flow-weighted output.
    6. Runs `flowweightaverage.FWA2` to compute the flow-weighted
       average and write it to that cleaned output path.
    7. Runs `flowweightaverage.FWA2` a second time, writing an
       additional copy under the F-part of the first input
       location, so plotting tools can associate it with that
       model.
    8. Re-reads all input data locations and selects the last one
       (Clear Creek) as the source for the tunnel heating step.
    9. Reads that series over the run time window, and for each
       timestamp, converts the HEC time to a Python datetime
       (correcting the 24:00 end-of-day convention to the following
       day's 00:00/first hour) to determine which calendar month it
       falls in.
    10. Looks up that month's heating offset from `monthly_heating`
        and adds it to the value.
    11. Writes the heated series to the second output location, then
        writes a second copy under the F-part of the original input
        model, so both the model output and a plotting-friendly,
        F-part-matched copy exist in DSS.
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
 
    locations = currentAlternative.getInputDataLocations()[:-1]
    locations = flowweightaverage.organizeLocations(currentAlternative, locations)
    currentAlternative.addComputeMessage('Found DSS paths:')

    # for each organized location pair, log every resolved DSS path
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    # Flow-weight average Lewiston outflow temps for linking
    # ------------------------------------------------------------------------------------

    # Set up the output path for the flow-weighted dam temperature
    # first outpath is flow-weighted dam temperature
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    #if len(outputlocations) > 1:
    #    currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    currentAlternative.addComputeMessage("Original output location 0 {0}".format(outputlocations[0]))
    currentAlternative.addComputeMessage("Original tsc fullname 0 {0}".format(outputpath.fullName))
    tspath =str(outputpath)
    
    currentAlternative.addComputeMessage("Original outpath 0 {0}".format(outputpath))
    tspath = tspath.split('/')
    fpart = tspath[6]
#    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
#    new_fpart += ':cequalw2-' #replace scripting with W2
#    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that

    # if the F-part contains a '|' separator, strip everything up to and including it
    if '|' in fpart:
        # remove everything before | including |
        tspath[6] = fpart[fpart.find('|')+1:]
    outputpath = '/'.join(tspath)
    #outputpath = str(outputpath)
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    # Compute the flow-weighted average and write it out twice
    # First write: the cleaned-up output path (main linking record).
    # Second write: a duplicate copy tagged with the first input location's F-part, purely so plotting tools can find it alongside that model's data.
    cfs_limit = 1.0 #float
    # compute and write the flow-weighted average temperature to the cleaned output path
    flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit, 10.0)

    # write a copy out using the input model f-part, so that plotting will be able to use it with the input model
    #currentAlternative.addComputeMessage("Len of locations: {0}".format(len(locations)))
    # ineffcient, could add extra copy in sommewhere else....
    outputpath = fixFpartToInput(str(locations[0]), outputpath)
    flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit, 10.0)

    

    # Add tunnel heating to Clear Creek tunnel temperatures
    # ------------------------------------------------------------------------------------ 

    # Re-fetch input locations and grab the Clear Creek record
    # we should have a flow and a temperature clear creek path, separate them
    if 'Flow' in str(clear_creek_locations[1]):
        ClearCreek_flow = clear_creek_locations[1]
        ClearCreek_temperature = clear_creek_locations[0]
    else:
        ClearCreek_flow = clear_creek_locations[0]
        ClearCreek_temperature = clear_creek_locations[1]

    # get the paths from the locations
    tspath_temperature =str(currentAlternative.loadTimeSeries(ClearCreek_temperature))
    tspath_temperature = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath_temperature)
    tspath_flow = str(currentAlternative.loadTimeSeries(ClearCreek_flow))
    tspath_flow = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath_flow)

    currentAlternative.addComputeMessage('Found DSS path: {0}'.format(tspath_temperature))
    currentAlternative.addComputeMessage('Found DSS path: {0}'.format(tspath_flow))

    currentAlternative.addComputeMessage('\n')
    dssFile = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    # get the temp and flow clear creek output paths
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath_temperature, outputpath_flow= DSS_Tools.organizeLocations(currentAlternative, outputlocations, ['ClearCreek_heated', 'ClearCreek_flow'], return_dss_paths=True)

    outputpath_temperature = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath_temperature))
    outputpath_flow = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath_flow))

    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath_temperature))
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath_flow))

    # Define the fixed monthly heating offsets (in Celsius)
    #heat in C
    monthly_heating = {1: 0.81, 2: 0.71, 3: 0.74, 4: 0.73, 5: 0.93, 6: 0.71, 7: 0.75, 8: 0.74, 9: 0.76, 10: 0.68, 11: 0.77, 12: 0.82}

    currentAlternative.addComputeMessage("\n##### PERFORMING HEATING #####")
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    # Read the source data and convert its timestamps to real dates
    dssFm = HecDss.open(dssFile)

    # read temperature data
    ts_temperature = dssFm.read(tspath_temperature, starttime_str, endtime_str, False)
    ts_temperature = ts_temperature.getData()
    hecstarttimes = ts_temperature.times
    readabledates = []

    # for each timestamp, convert the HEC time to a Python datetime, correcting for the 24:00 convention
    for i in range(len(ts_temperature.times)):
        hecdate = ts_temperature.getHecTime(i).dateAndTime(4) #delete later
        if hecdate.split(':')[0][-2:] == '24':
            hecdate = hecdate.split(',')[0] + ', 23:00'
            dthecdate = dt.datetime.strptime(hecdate, '%d%b%Y, %H:%M') #31Dec2018, 24:00
            dthecdate += dt.timedelta(hours=1)
        else:
            dthecdate = dt.datetime.strptime(hecdate, '%d%b%Y, %H:%M') #31Dec2018, 24:00
        readabledates.append(dthecdate)

    # read flow data
    ts_flow = dssFm.read(tspath_flow, starttime_str, endtime_str, False)
    ts_flow = ts_flow.getData()

    new_values = []
    flow_values = []

    for vi, val in enumerate(ts_temperature.values):

        # if the temperature value is missing it will be a very negative number
        # detect this and set temperature and flow to zero
        if val <= -100:
            new_values.append(0)
            flow_values.append(0)
        else:
            # Apply the monthly heating offset to every value
            # for each value, look up its month's heating offset and add it to the value
            valmonth = readabledates[vi].month
            heat_amount_c = monthly_heating[valmonth]
            new_values.append(val + heat_amount_c)
            flow_values.append(ts_flow.values[vi])

    # Package the heated results and write them to DSS
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputpath_temperature
    tsc.values = new_values
    tsc.startTime = hecstarttimes[0]
    tsc.units = ts_temperature.units
    tsc.type = ts_temperature.type
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(new_values)
    tsc.startHecTime = rtw.getStartTime()
    tsc.endHecTime = rtw.getEndTime()
    dssFm.write(tsc)

    # write a copy out using the input model f-part, so that plotting will be able to use it with the input model
    #currentAlternative.addComputeMessage("Len of locations: {0}".format(len(locations)))
    tsc.fullName = fixFpartToInput(tspath_temperature, str(outputpath_temperature))
    dssFm.write(tsc)
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(new_values)))

    # same for flow data
    tsc_flow = TimeSeriesContainer()
    tsc_flow.times = hecstarttimes
    tsc_flow.fullName = outputpath_flow
    tsc_flow.values = flow_values
    tsc_flow.startTime = hecstarttimes[0]
    tsc_flow.units = ts_flow.units
    tsc_flow.type = ts_flow.type
    tsc_flow.endTime = hecstarttimes[-1]
    tsc_flow.numberValues = len(flow_values)
    tsc_flow.startHecTime = rtw.getStartTime()
    tsc_flow.endHecTime = rtw.getEndTime()
    dssFm.write(tsc_flow)

    # write a copy out using the input model f-part, so that plotting will be able to use it with the input model
    # currentAlternative.addComputeMessage("Len of locations: {0}".format(len(locations)))
    tsc_flow.fullName = fixFpartToInput(tspath_flow, str(outputpath_flow))
    dssFm.write(tsc_flow)

    dssFm.close()
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(flow_values)))

    # exit
    return True

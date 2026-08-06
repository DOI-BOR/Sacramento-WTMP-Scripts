import sys
from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
from com.rma.model import Project
import math
import datetime as dt
import DSS_Tools
reload(DSS_Tools)

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
    location_fpart = inputpath.split('/')[6]
    out_parts = outpath.split('/')
    out_parts[6] = location_fpart
    return '/'.join(out_parts)

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
    
    locations = currentAlternative.getInputDataLocations()
    # validate that enough input data locations are configured
    if len(locations) == 1:
        currentAlternative.addComputeMessage("Found only 1 datapath locations. Need at least 2.")
        sys.exit(1)
    elif len(locations) == 0:
        currentAlternative.addComputeMessage("Found no datapath locations. Exiting.")
        sys.exit(1)
    
    base_ts = locations[0] #this is just the first, we'll add everything to this one 
    base_tspath, base_DSSPath = DSS_Tools.getDataLocationDSSInfo(base_ts, currentAlternative, computeOptions)   
    currentAlternative.addComputeMessage('Found DSS path: {0}'.format(base_tspath))
    
    currentAlternative.addComputeMessage('\n')
    rtw = computeOptions.getRunTimeWindow()

    # Adding flow to keep keswick balanced, coming from independent DSS files, so can't
    # use DSS_Tools (currently).  Only add first two input locations!!
    # ------------------------------------------------------------------------------------ 

    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])    
    
    # if extra output locations were configured, note that only the first is used
    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    currentAlternative.addComputeMessage("\n##### Adding Timeseries #####")
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    # read the base series' data over the run time window, then close the source file
    dss_base = HecDss.open(base_DSSPath)
    base_TS = dss_base.read(base_tspath, starttime_str, endtime_str, False)
    dss_base.close()
    base_TS = base_TS.getData()
    hecstarttimes = base_TS.times
    all_values = base_TS.values
    units = base_TS.units
    dsstype = base_TS.type

    # we rely here on the input from W2 Whisketytown and the debris dam flow from the DMS, which might be in
    # different units

    cms2cfs_base = 1.0
    # if the base series is in cms, convert it to cfs
    if str(units).lower() == 'cms':
        cms2cfs_base = 35.314666212661

    # for the second input location, read its series and add it to the running total
    for other_loc in locations[1:2]:
        new_values = []
        other_tspath, other_DSSPath = DSS_Tools.getDataLocationDSSInfo(other_loc, currentAlternative, computeOptions)
        currentAlternative.addComputeMessage("Adding {0}".format(other_tspath))
        currentAlternative.addComputeMessage('DSS path: {0}'.format(other_DSSPath))
        dss_other = HecDss.open(other_DSSPath)
        other_TS = dss_other.read(other_tspath, starttime_str, endtime_str, False)
        dss_other.close()
        other_TS = other_TS.getData()

        cms2cfs_other = 1.0
        # if the other series is in cms, convert it to cfs
        if str(other_TS.units).lower() == 'cms':
            cms2cfs_other = 35.314666212661
        
        other_values = other_TS.values
        
        # for each value in the base series, add the corresponding converted value from the other series
        for vi, val in enumerate(all_values):
            try:
                new_values.append(val*cms2cfs_base + other_values[vi]*cms2cfs_other)
            except:
                currentAlternative.addComputeMessage("No value for location {0} at idx {1} {2}".format(other_loc, vi, hecstarttimes[vi]))
                new_values.append(MISSING_DOUBLE)
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

    locations = currentAlternative.getInputDataLocations()
    #if len(locations) > 1:
    #    currentAlternative.addComputeMessage("Found more than 1 datapath locations. Using the first, {0}".format(outputlocations[0]))
    #elif len(locations) == 0:
    #    currentAlternative.addComputeMessage("Found no datapath locations. Exiting.")
    #    sys.exit(1)
    
    # resolve the Spring Creek input location's DSS path, correcting its F-part
    SpringCreek = locations[-1]
   
    tspath =str(currentAlternative.loadTimeSeries(SpringCreek))
    tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)
            
    currentAlternative.addComputeMessage('Found DSS path: {0}'.format(tspath))
    
    currentAlternative.addComputeMessage('\n')
    dssFile = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    # second outpath is CC tuneel w/ heating
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[1])
    
    #if len(outputlocations) > 1:
    #    currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))

    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    currentAlternative.addComputeMessage("\n##### PERFORMING HEATING #####")
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    # read the Spring Creek series over the run time window, to be heated below
    dssFm = HecDss.open(dssFile)
    TS = dssFm.read(tspath, starttime_str, endtime_str, False)
    TS = TS.getData()    
    hecstarttimes = TS.times

    heat_amount_c = 0.32
    new_values = []
    
    # for each value in the Spring Creek series, add the fixed heating offset
    for val in TS.values:
        new_values.append(val + heat_amount_c)
        
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

    # write a copy out using the input model f-part, so that plotting will be able to use it with the input model
    tsc.fullName = fixFpartToInput(tspath, outputpath)
    dssFm.write(tsc)    
    
    dssFm.close()
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(new_values)))

    # exit
    return True

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
    # get F-part from input locations
    location_fpart = inputpath.split('/')[6]
    out_parts = outpath.split('/')
    out_parts[6] = location_fpart
    return '/'.join(out_parts)

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    
    locations = currentAlternative.getInputDataLocations()
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
    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    currentAlternative.addComputeMessage("\n##### Adding Timeseries #####")
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dss_base = HecDss.open(base_DSSPath)
    base_TS = dss_base.read(base_tspath, starttime_str, endtime_str, False)
    dss_base.close()
    base_TS = base_TS.getData()
    hecstarttimes = base_TS.times
    all_values = base_TS.values
    units = base_TS.units
    dsstype = base_TS.type

    for other_loc in locations[1:2]:
        new_values = []
        other_tspath, other_DSSPath = DSS_Tools.getDataLocationDSSInfo(other_loc, currentAlternative, computeOptions)
        currentAlternative.addComputeMessage("Adding {0}".format(other_tspath))
        currentAlternative.addComputeMessage('DSS path: {0}'.format(other_DSSPath))
        dss_other = HecDss.open(other_DSSPath)
        other_TS = dss_other.read(other_tspath, starttime_str, endtime_str, False)
        dss_other.close()
        other_TS = other_TS.getData()
        other_values = other_TS.values
        for vi, val in enumerate(all_values):
            try:
                new_values.append(val + other_values[vi])
            except:
                currentAlternative.addComputeMessage("No value for location {0} at idx {1} {2}".format(other_loc, vi, hecstarttimes[vi]))
                new_values.append(MISSING_DOUBLE)
        all_values = new_values
        
    dssfn = computeOptions.getDssFilename()
    dssFm = HecDss.open(dssfn)
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputpath
    tsc.values = all_values
    tsc.startTime = hecstarttimes[0]
    tsc.units = units
    tsc.type = dsstype
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(all_values)
    tsc.startHecTime = rtw.getStartTime()
    tsc.endHecTime = rtw.getEndTime()
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

    #heat in C
    monthly_heating = {1: 0.81, 2: 0.71, 3: 0.74, 4: 0.73, 5: 0.93, 6: 0.71, 7: 0.75, 8: 0.74, 9: 0.76, 10: 0.68, 11: 0.77, 12: 0.82}

    currentAlternative.addComputeMessage("\n##### PERFORMING HEATING #####")
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    TS = dssFm.read(tspath, starttime_str, endtime_str, False)
    TS = TS.getData()    
    hecstarttimes = TS.times

    heat_amount_c = 0.32
    new_values = []
    for val in TS.values:
        new_values.append(val + heat_amount_c)
        
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
    dssFm.write(tsc)

    # write a copy out using the input model f-part, so that plotting will be able to use it with the input model
    tsc.fullName = fixFpartToInput(tspath, outputpath)
    dssFm.write(tsc)    
    
    dssFm.close()
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(new_values)))

    # exit
    return True

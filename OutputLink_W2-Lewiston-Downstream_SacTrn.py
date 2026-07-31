import sys
print(sys.path)
import BoundaryFixes
reload(BoundaryFixes)
from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
import math
import datetime as dt
import flowweightaverage
reload(flowweightaverage)
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
 
    all_locations = currentAlternative.getInputDataLocations()
    locations = []
    clear_creek_locations = []

    # sort into clear creek (seg 111) and the rest
    for location in all_locations:
        if '111' in str(location):
            clear_creek_locations.append(location)
        else:
            locations.append(location)

    locations = flowweightaverage.organizeLocations(currentAlternative, locations)
    currentAlternative.addComputeMessage('Found DSS paths:')
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    # Flow-weight average Lewiston outflow temps for linking
    # ------------------------------------------------------------------------------------

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
    if '|' in fpart:
        # remove everything before | including |
        tspath[6] = fpart[fpart.find('|')+1:]
    outputpath = '/'.join(tspath)
    #outputpath = str(outputpath)
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    cfs_limit = 1.0 #float
    flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit, 10.0)

    # write a copy out using the input model f-part, so that plotting will be able to use it with the input model
    #currentAlternative.addComputeMessage("Len of locations: {0}".format(len(locations)))
    # ineffcient, could add extra copy in sommewhere else....
    outputpath = fixFpartToInput(str(locations[0]), outputpath)
    flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit, 10.0)

    

    # Add tunnel heating to Clear Creek tunnel temperatures
    # ------------------------------------------------------------------------------------ 

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

    #heat in C
    monthly_heating = {1: 0.81, 2: 0.71, 3: 0.74, 4: 0.73, 5: 0.93, 6: 0.71, 7: 0.75, 8: 0.74, 9: 0.76, 10: 0.68, 11: 0.77, 12: 0.82}

    currentAlternative.addComputeMessage("\n##### PERFORMING HEATING #####")
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    currentAlternative.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)

    # read temperature data
    ts_temperature = dssFm.read(tspath_temperature, starttime_str, endtime_str, False)
    ts_temperature = ts_temperature.getData()
    hecstarttimes = ts_temperature.times
    readabledates = []

    # get timeseries times
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
            valmonth = readabledates[vi].month
            heat_amount_c = monthly_heating[valmonth]
            new_values.append(val + heat_amount_c)
            flow_values.append(ts_flow.values[vi])

    # write temperature data
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

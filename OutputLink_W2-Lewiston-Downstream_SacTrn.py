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

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
 
    locations = currentAlternative.getInputDataLocations()[:-1]
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
    flowweightaverage.FWA(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit)
    

    # Add tunnel heating to Clear Creek tunnel temperatures
    # ------------------------------------------------------------------------------------ 

    locations = currentAlternative.getInputDataLocations()
    #if len(locations) > 1:
    #    currentAlternative.addComputeMessage("Found more than 1 datapath locations. Using the first, {0}".format(outputlocations[0]))
    #elif len(locations) == 0:
    #    currentAlternative.addComputeMessage("Found no datapath locations. Exiting.")
    #    sys.exit(1)
    
    ClearCreek = locations[-1]
   
    tspath =str(currentAlternative.loadTimeSeries(ClearCreek))
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
    readabledates = []
    for i in range(len(TS.times)):
        hecdate = TS.getHecTime(i).dateAndTime(4) #delete later
        if hecdate.split(':')[0][-2:] == '24':
            hecdate = hecdate.split(',')[0] + ', 23:00'
            dthecdate = dt.datetime.strptime(hecdate, '%d%b%Y, %H:%M') #31Dec2018, 24:00
            dthecdate += dt.timedelta(hours=1)
        else:
            dthecdate = dt.datetime.strptime(hecdate, '%d%b%Y, %H:%M') #31Dec2018, 24:00
        readabledates.append(dthecdate)
        
    new_values = []
    for vi, val in enumerate(TS.values):
        valmonth = readabledates[vi].month
        heat_amount_c = monthly_heating[valmonth]
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
    dssFm.close()
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(new_values)))

    # exit
    return True

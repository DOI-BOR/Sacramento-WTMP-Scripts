import sys
import BoundaryFixes
reload(BoundaryFixes)
from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
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

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    
    locations = currentAlternative.getInputDataLocations()
    if len(locations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 datapath locations. Using the first, {0}".format(outputlocations[0]))
    elif len(locations) == 0:
        currentAlternative.addComputeMessage("Found no datapath locations. Exiting.")
        sys.exit(1)
    
    SpringCreek = locations[0]
   
    tspath =str(currentAlternative.loadTimeSeries(SpringCreek))
    tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)
            
    currentAlternative.addComputeMessage('Found DSS path: {0}'.format(tspath))
    
    currentAlternative.addComputeMessage('\n')
    dssFile = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    
    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))

    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    #heat in C

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
    tsc.units = 'c'
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(new_values)
    tsc.startHecTime = rtw.getStartTime()
    tsc.endHecTime = rtw.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlternative.addComputeMessage("Number of Written values: {0}".format(len(new_values)))
    return True


            


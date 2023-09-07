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
    locations = flowweightaverage.organizeLocations(currentAlternative, locations)
    currentAlternative.addComputeMessage('Found DSS paths:')
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    outputFpart = 'PostProcessed'
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    tspath =str(outputpath)
    tspath = tspath.split('/')
    fpart = tspath[6]
#    new_fpart = fpart.lower().split(':scripting-')[0] #fpart will have scripting in it, but we want w2 version
#    new_fpart += ':cequalw2-' #replace scripting with W2
#    new_fpart += '-'.join(computeOptions.getSimulationName().split('-')[:-1]) #sim name will have the sim group, so snip that
    tspath[6] = outputFpart
    outputpath = '/'.join(tspath)
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    cfs_limit = 1.0 #float
    flowweightaverage.FWA(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit)
    
    return True


            


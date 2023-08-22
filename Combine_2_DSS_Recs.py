import sys
print(sys.path)
import Simple_DSS_Functions
reload(Simple_DSS_Functions)

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

def fixInputLocationFpart(currentAlternative, tspath):
    new_fpart_start = ':'.join(currentAlternative.getInputFPart().split(':')[:-1])
    tspath = tspath.split('/')
    fpart = tspath[6]
    fpart_split = fpart.split(':')
    new_fpart = new_fpart_start + ':' + fpart_split[-1]
    tspath[6] = new_fpart
    tspath = '/'.join(tspath)
    return tspath

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    locations = currentAlternative.getInputDataLocations()
    
    currentAlternative.addComputeMessage('Found DSS paths:')
    location_paths = []
    for location in locations:
        currentAlternative.addComputeMessage(str(location))
        tspath =str(currentAlternative.loadTimeSeries(location))
        loc_path = fixInputLocationFpart(currentAlternative, tspath)
        location_paths.append(loc_path)
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    outputFpart = 'PostProcessed'
    outputlocations = currentAlternative.getOutputDataLocations()

    # We combine 2 DSS records, so (input) locations should be in sequential pairs, and outputpaths should = locations/2
    for i in range(len(outputlocations)):
        outputpath = currentAlternative.createOutputTimeSeries(outputlocations[i])
        if len(outputlocations) > 1:
            currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
        tspath =str(outputpath)
        tspath = tspath.split('/')
        fpart = tspath[6]
        tspath[6] = outputFpart
        outputpath = '/'.join(tspath)
        currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))
        lpi = 2*i
        Simple_DSS_Functions.add_DSS_Data(currentAlternative, dss_file, rtw, location_paths[lpi:lpi+1], outputpath)

    return True


            


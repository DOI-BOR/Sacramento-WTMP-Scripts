import sys
print(sys.path)
import DSS_Tools
reload(DSS_Tools)
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

def fixFpartToInput(locations_paths, outpath):
    # get F-part from input locations
    location_fpart = locations_paths[0][0].split('/')[6]
    out_parts = outpath.split('/')
    out_parts[6] = location_fpart
    return '/'.join(out_parts)

def getOutputPaths(locations_paths,currentAlternative):
    outputlocations = currentAlternative.getOutputDataLocations()
    outputPaths = []
    for opl in outputlocations:
        op_ts = currentAlternative.createOutputTimeSeries(opl)
        outputPaths.append(fixFpartToInput(locations_paths,str(op_ts)))
    return outputPaths
    
def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    locations = currentAlternative.getInputDataLocations()
    locations_paths = flowweightaverage.organizeLocations(currentAlternative, locations)    
    output_paths = getOutputPaths(locations_paths,currentAlternative)
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    # Clear Creek Tunnel Heating
    DSS_Tools.add_DSS_Data(currentAlternative, dss_file, rtw, locations_paths[0], output_paths[0])

    # Spring creek tunnel heating
    DSS_Tools.add_DSS_Data(currentAlternative, dss_file, rtw, locations_paths[1], output_paths[1])
 
    return True


            


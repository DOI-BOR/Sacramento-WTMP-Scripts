import sys
print(sys.path)
from hec.heclib.dss import HecDss

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

input_names = [
'WQ_BOUND_TEMP_SRC_CCTUNNEL',
'LEWISTON RESERVOIR-TUNNEL',
'WQ_BOUND_TEMP_SRC_SCTUNNEL',
'WHISKEYTOWN LAKE-TUNNEL'
]

output_names = [
'LewistonTunnel_With_Heating',
'WhiskeyTownTunnel_With_Heating'
]

#copy_to_resim_names = [
#'TargetTemp',
#]

ResSim_linked_rec = 'LEWISTON RESERVOIR-TUNNEL' # only used for get f part, could be any W2 output

def fixFpartToInput(locations_paths, outpath):
    # get F-part from input locations
    location_fpart = locations_paths[0].split('/')[6]
    out_parts = outpath.split('/')
    out_parts[6] = location_fpart
    return '/'.join(out_parts)

def getOutputPaths(locations_paths,currentAlternative):
    outputlocations_obs = currentAlternative.getOutputDataLocations()    
    outputPaths = []
    for opl in outputlocations_obs:        
        path = str(currentAlternative.createOutputTimeSeries(opl))
        outputPaths.append(fixFpartToInput(locations_paths,path))
    return outputPaths
    
def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    locations_obj = currentAlternative.getInputDataLocations()
    locations_paths = DSS_Tools.organizeLocations(currentAlternative, locations_obj, input_names, return_dss_paths=True)

    ResSim_FPart = DSS_Tools.organizeLocations(currentAlternative, locations_obj, [ResSim_linked_rec], return_dss_paths=True)[0].split('/')[6]
    print('ResSim FPart:',ResSim_FPart)

    # this is choking on linked DSS data - can't get linking to pass in the right thing.  Hardcoding ressim_copy_paths
    #ressim_copy_paths = DSS_Tools.organizeLocations(currentAlternative, locations_obj, copy_to_resim_names, return_dss_paths=True)
    ressim_copy_paths = []
    
    output_paths = getOutputPaths(locations_paths,currentAlternative)
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()

    # Clear Creek Tunnel Heating
    print('Clear Creek',locations_paths[:2],output_paths[0])
    DSS_Tools.add_DSS_Data(currentAlternative, dss_file, rtw, locations_paths[:2], output_paths[0])

    # Spring creek tunnel heating
    print('Spring Creek',locations_paths[2:],output_paths[1])
    DSS_Tools.add_DSS_Data(currentAlternative, dss_file, rtw, locations_paths[2:], output_paths[1])

  
    # copy other recs to ResSim fpart - there is potential for overwriting here, but seems
    # unlikely that two models would have the same rec name that you want to copy?
    for loc_path in ressim_copy_paths:
        DSS_Tools.copy_dss_ts(loc_path,new_fpart=ResSim_FPart,dss_file_path=dss_file)

 
    return True


            


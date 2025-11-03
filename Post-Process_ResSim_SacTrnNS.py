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

model_recs = [ 'alt_model_output1','alt_model_output2'] # linked to records from upstream simulation models
copy_rec = ['target_temp','target_temp_upstream'] # Linked to a DSS record in DSS file, not a model result

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
    ressim_copy_paths = ["/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/",]
    
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

    
    # Copy the opy_rec (target temp) to thoer f-parts for easier reporting
    copy_obj = DSS_Tools.organizeLocations(currentAlternative, locations_obj, copy_rec, return_dss_paths=False)
    copy_dss_rec,copy_dss_filepath = DSS_Tools.getDataLocationDSSInfo(copy_obj[0], currentAlternative, computeOptions)
    print(copy_dss_filepath,copy_dss_rec)    

    # get the f-part from the linked model record, copy the copy_rec to model-fpart in the simulation DSS
#    for model in model_recs:
#        model_obj = DSS_Tools.organizeLocations(currentAlternative, locations_obj, [model], return_dss_paths=False)
#        dss_model_path, _ = DSS_Tools.getDataLocationDSSInfo(model_obj[0], currentAlternative, computeOptions)
#        model_fpart = dss_model_path.split('/')[6]
#        print('model_fpart: ',model_fpart)
#        DSS_Tools.copy_dss_ts(copy_dss_rec,
#                              new_fpart=model_fpart,
#                              dss_file_path=copy_dss_filepath,
#                              dss_file_alt_outpath=dss_file)
                              
    # get the f-part from the linked model record, copy the copy_rec to model-fpart in the simulation DSS
    for co in copy_obj:
        copy_dss_rec,copy_dss_filepath = DSS_Tools.getDataLocationDSSInfo(co, currentAlternative, computeOptions)
        print(copy_dss_filepath,copy_dss_rec)    
    
        # get the f-part from the linked model record, copy the copy_rec to model-fpart in the simulation DSS
        for model in model_recs:
            model_obj = DSS_Tools.organizeLocations(currentAlternative, locations_obj, [model], return_dss_paths=False)
            dss_model_path, _ = DSS_Tools.getDataLocationDSSInfo(model_obj[0], currentAlternative, computeOptions)
            model_fpart = dss_model_path.split('/')[6]
            print('model_fpart: ',model_fpart)
            DSS_Tools.copy_dss_ts(copy_dss_rec,
                                  new_fpart=model_fpart,
                                  dss_file_path=copy_dss_filepath,
                                  dss_file_alt_outpath=dss_file)
 
    return True


            


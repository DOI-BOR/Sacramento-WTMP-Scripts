import sys,os
print(sys.path)
from hec.heclib.dss import HecDss

from com.rma.model import Project
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import DSS_Tools
reload(DSS_Tools)

model_recs = [ 'alt_model_output1','alt_model_output2'] # linked to records from upstream simulation models
copy_rec = ['target_temp','target_temp_upstream'] # Linked to a DSS record in DSS file, not a model result
    
def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
   
    dss_file = computeOptions.getDssFilename()  
    
    locations_obj = currentAlternative.getInputDataLocations()
    copy_obj = DSS_Tools.organizeLocations(currentAlternative, locations_obj, copy_rec, return_dss_paths=False)
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


            


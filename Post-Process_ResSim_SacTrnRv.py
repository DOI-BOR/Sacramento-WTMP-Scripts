"""
Post-Process_ResSim_SacTrnRv
=============================

Compute a WAT Scripting Alternative that copies linked target temperature
DSS records into the simulation DSS file, tagged with the F-part of each
linked upstream model output record.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `hec.heclib.dss`, `com.rma.model.Project`, `DSS_Tools`.
"""

import sys,os                                          # Standard library: sys.path diagnostics and path joins
print(sys.path)                                        # Diagnostic: show module search path at script load time
from hec.heclib.dss import HecDss                      # HEC-DSS: open/read/write DSS files (imported for completeness; not directly used)

from com.rma.model import Project                      # WAT project API: resolve workspace path for scripts folder
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))  # Ensure 'scripts' folder is importable

import DSS_Tools                                       # Local module: DSS path resolution, F-part fixing, copy_dss_ts
reload(DSS_Tools)                                      # Jython: ensure latest version is loaded

# linked to records from upstream simulation models
model_recs = [ 'alt_model_output1','alt_model_output2']
# Linked to a DSS record in DSS file, not a model result 
copy_rec = ['target_temp','target_temp_upstream'] 
    
def computeAlternative(currentAlternative, computeOptions):
    """
    Compute a scripting alternative that copies linked target
    temperature records to each linked model's F-part.

    Copies linked target temperature DSS records into the
    simulation DSS file, tagged with the F-part of each linked
    upstream model output record.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()` and `getInputDataLocations()` for
        logging and retrieving linked data locations.
    computeOptions : object
        The compute options/settings object, used to retrieve the
        DSS filename for this compute.

    Returns
    -------
    bool
        True once all target temperature records have been copied
        for all linked model records.

    Notes
    -----
    Workflow:

    1. Logs a compute status message to the alternative.
    2. Retrieves the DSS filename for this compute and the
       alternative's input data locations.
    3. Organizes the locations tied to `copy_rec` (the target
       temperature records to be copied).
    4. For each of those target temperature records:

       - Resolves its source DSS path and file.
       - For each linked upstream model record in `model_recs`:

         - Resolves the model record's DSS path to extract its
           F-part.
         - Copies the target temperature record into the simulation
           DSS file, renamed with that model's F-part.
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
   
    dss_file = computeOptions.getDssFilename()  
    
    locations_obj = currentAlternative.getInputDataLocations()
    copy_obj = DSS_Tools.organizeLocations(currentAlternative, locations_obj, copy_rec, return_dss_paths=False)
    
    # for each target temperature record to be copied
    for co in copy_obj:
        copy_dss_rec,copy_dss_filepath = DSS_Tools.getDataLocationDSSInfo(co, currentAlternative, computeOptions)
        print(copy_dss_filepath,copy_dss_rec)    
    
        # get the f-part from the linked model record, copy the copy_rec to model-fpart in the simulation DSS
        # for each linked upstream model output record
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
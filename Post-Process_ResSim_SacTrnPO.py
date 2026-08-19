"""
Post-Process_ResSim_SacTrnPO
=============================

Compute a WAT Scripting Alternative that applies tunnel heating adjustments
to Clear Creek and Spring Creek tunnel temperature data, and optionally
copies additional records into ResSim using a shared F-part.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `hec.heclib.dss`, `DSS_Tools`, `flowweightaverage`.
"""

import sys                       # Standard library: used here only to print sys.path for diagnostics
print(sys.path)                  # Diagnostic: show module search path at script load time
from hec.heclib.dss import HecDss  # HEC-DSS: open/read/write DSS files (imported for completeness; not directly used)

import DSS_Tools                 # Local module: DSS path resolution, F-part fixing, add_DSS_Data, copy_dss_ts
reload(DSS_Tools)                # Jython: ensure latest version is loaded
import flowweightaverage         # Local module: flow-weighted average utilities (imported but not directly used here)
reload(flowweightaverage)        # Reload to ensure latest version

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

# set input names
input_names = [
'WQ_BOUND_TEMP_SRC_CCTUNNEL',
'LEWISTON RESERVOIR-TUNNEL',
'WQ_BOUND_TEMP_SRC_SCTUNNEL',
'WHISKEYTOWN LAKE-TUNNEL'
]

# set output names
output_names = [
'LewistonTunnel_With_Heating',
'WhiskeyTownTunnel_With_Heating'
]

#copy_to_resim_names = [
#'TargetTemp',
#]

ResSim_linked_rec = 'LEWISTON RESERVOIR-TUNNEL' # only used for get f part, could be any W2 output

def fixFpartToInput(locations_paths, outpath):
    """
    Replace the F-part of a DSS output path with the F-part taken
    from the first input location path.

    This keeps output records tagged with the same F-part (typically
    a model/run identifier) as the inputs they were derived from,
    rather than whatever default F-part was assigned when the output
    path was created.

    Parameters
    ----------
    locations_paths : list of str
        DSS path strings for the input locations. Only the first
        entry's F-part (index 6 after splitting on `'/'`) is used.
    outpath : str
        The DSS output path whose F-part should be replaced.

    Returns
    -------
    str
        The output DSS path with its F-part replaced by the F-part
        from `locations_paths[0]`.
    """
    # get F-part from input locations
    location_fpart = locations_paths[0].split('/')[6]
    out_parts = outpath.split('/')
    out_parts[6] = location_fpart
    return '/'.join(out_parts)

def getOutputPaths(locations_paths,currentAlternative):
    """
    Build the list of output DSS paths for this alternative, each
    tagged with the F-part from the given input location paths.

    For every output data location configured on the alternative,
    this creates the corresponding output time series path and then
    rewrites its F-part to match the input locations, via
    `fixFpartToInput`.

    Parameters
    ----------
    locations_paths : list of str
        DSS path strings for the input locations, used to source the
        F-part to apply.
    currentAlternative : object
        The alternative object being computed. Must support
        `getOutputDataLocations()` and `createOutputTimeSeries()`.

    Returns
    -------
    list of str
        DSS output paths, one per output data location, each with
        its F-part fixed to match the input locations.
    """
    outputlocations_obs = currentAlternative.getOutputDataLocations()    
    outputPaths = []
    # for each configured output data location, build and fix its output path
    for opl in outputlocations_obs:        
        path = str(currentAlternative.createOutputTimeSeries(opl))
        outputPaths.append(fixFpartToInput(locations_paths,path))
    return outputPaths
    
def computeAlternative(currentAlternative, computeOptions):
    """
    Compute a scripting alternative for tunnel heating with optional
    ResSim record copying.

    Applies tunnel heating adjustments to Clear Creek and Spring
    Creek tunnel temperature data, and optionally copies additional
    records into ResSim using a shared F-part.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()`, `getInputDataLocations()`,
        `getOutputDataLocations()`, and `createOutputTimeSeries()`
        for logging and resolving linked data locations.
    computeOptions : object
        The compute options/settings object. Must support
        `getDssFilename()` and `getRunTimeWindow()` to provide the
        target DSS file and the time window to compute over.

    Returns
    -------
    bool
        True once the tunnel heating calculations have been applied
        and any ResSim copy records have been written.

    Notes
    -----
    Workflow:

    1. Logs a compute status message to the alternative.
    2. Resolves DSS paths for all configured input locations
       (`input_names`).
    3. Resolves the F-part associated with the ResSim-linked record
       (`ResSim_linked_rec`), for use when copying any additional
       records into ResSim.
    4. Builds the list of output DSS paths, each tagged with the
       F-part from the input locations (via `getOutputPaths`).
    5. Retrieves the DSS filename and run time window for this
       compute.
    6. Applies the Clear Creek tunnel heating calculation using the
       first two input paths, writing to the first output path.
    7. Applies the Spring Creek tunnel heating calculation using the
       remaining input paths, writing to the second output path.
    8. Copies any additional records in `ressim_copy_paths` into the
       simulation DSS file, tagged with the ResSim F-part.
       (Currently disabled/hardcoded to an empty list due to an
       issue linking DSS data for this step - see inline note.)
    """
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
"""
OutputLink_W2-Trinity-W2-Lewiston_SacTrn
==========================================

Compute a WAT Scripting Alternative that applies a flow-weighted average
(FWA) to configured input data locations for the Trinity-to-Lewiston W2
link, using the last input record pair as an override source for values
above a small cfs threshold.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `flowweightaverage`, `BoundaryFixes`, `DSS_Tools`.
"""

import sys                       # Standard library: retained (not directly used beyond potential diagnostics)
import flowweightaverage         # Local module: flow-weighted average temperature computation (FWA2)
reload(flowweightaverage)        # Jython: ensure latest version is loaded
import BoundaryFixes             # Local module: threshold/NaN replacement helpers (used only in unreachable code below)
reload(BoundaryFixes)            # Reload to ensure latest version
import DSS_Tools                 # Local module: DSS path resolution / F-part fixing helpers
reload(DSS_Tools)                # Reload to ensure latest version

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
    """
    Compute a scripting alternative that applies a flow weighted
    average (FWA) to the configured input data locations.

    Uses the last record pair as an override source for values
    above a small cfs threshold, and writes the result to the first
    output location.

    Parameters
    ----------
    currentAlternative : object
        The alternative object being computed. Must support
        `addComputeMessage()`, `getInputDataLocations()`,
        `getOutputDataLocations()`, `createOutputTimeSeries()`, and
        `loadTimeSeries()` for logging and resolving linked data
        locations.
    computeOptions : object
        The compute options/settings object. Must support
        `getDssFilename()` and `getRunTimeWindow()` to provide the
        target DSS file and the time window to compute over.

    Returns
    -------
    bool
        True once the flow weighted average has been computed and
        written to the first output location.

    Notes
    -----
    Workflow:

    1. Logs a compute status message to the alternative.
    2. Retrieves and organizes the input data locations into record
       pairs via `flowweightaverage.organizeLocations`.
    3. Logs the resolved DSS paths for each organized location pair.
    4. Retrieves the DSS filename and run time window.
    5. Resolves the first output data location's DSS path,
       correcting its F-part, and warns if more than one output
       location is configured (only the first is used).
    6. Runs `flowweightaverage.FWA2` with a very low `cfs_limit`
       (1.0), using the last record pair in `locations` as an
       override source (`last_override=True`) rather than treating
       it as another value to average in.
    7. Returns True immediately after this step.

    **Known issue (left unchanged per instructions):** This function
    contains a large block of code after the initial `return True`
    statement (`structure1Flow_location`, `totalTemperature_location`
    filtering, a second `FWA2` call, and the "REPLACING VALUES" /
    "REPLACING REMAINING VALUES" sections). Because of the early
    return, none of that code can ever execute - it is dead code as
    written. It has been left in place and documented below in case
    it is intended to be reachable in a future version of this
    script.
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )

    locations = currentAlternative.getInputDataLocations()
    locations = flowweightaverage.organizeLocations(currentAlternative, locations)
    currentAlternative.addComputeMessage('Found DSS paths:')
    
    # for each organized location pair, log every resolved DSS path
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    
    # if extra output locations were configured, note that only the first is used
    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))
    cfs_limit = 1.0 #float
#    currentAlternative.addComputeMessage("\n##### PERFORMING FLOW WEIGHTED AVERAGE FOR UNDER 1CFS #####")
    # last record-pair in locations, if over cfs_limit, replaces FWA
    flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit, 10.0, last_override=True)

    return True


    # NOTE: everything below this point is unreachable, since the function already returned True above
    # ------------------------------------------------------------------------------------------------
    # The following block (structure1Flow_location / totalTemperature_location filtering, a second
    # FWA2 call at a much lower cfs_limit, and the commented-out "replace values over threshold" /
    # "replace remaining NaN values" steps) is dead code: it can never execute because of the
    # unconditional `return True` immediately above. It is preserved and documented here, unmodified,
    # in case it is meant to become reachable again in a future revision of this script.
    #
    # Intended purpose of the dead code (based on its structure and variable names):
    #   1. Filter out the 'STR1_Flow' and 'Total_Temp' input locations from the full input list,
    #      leaving only genuine flow/temperature pairs.
    #   2. Organize the filtered locations into pairs as before.
    #   3. Compute a flow-weighted average with a much smaller cfs_limit (0.001) than the reachable
    #      code above (1.0), writing to the same first output location.
    #   4. (Commented out) Replace any values in the output above a 1 cfs threshold using
    #      BoundaryFixes.replaceValuesOverThresh, sourced from the 'STR1_Flow' and 'Total_Temp'
    #      locations.
    #   5. (Commented out) Replace any remaining NaN/undefined values in the output using
    #      BoundaryFixes.replaceNaNValues, sourced from a 'WD2_Temp' location.

    structure1Flow_location = 'STR1_Flow' #used later to filter
    totalTemperature_location = 'Total_Temp' #used later to filter
    locations = currentAlternative.getInputDataLocations()

    filtered_locations = []
    # for each input location, keep only those not matching the structure flow or temperature locations
    for location in locations:
        if str(location) not in [structure1Flow_location, totalTemperature_location]: #filter out the later locations
            filtered_locations.append(location)
            
    locations = flowweightaverage.organizeLocations(currentAlternative, filtered_locations)

    currentAlternative.addComputeMessage('Found DSS paths:')
    # for each organized location pair, log every resolved DSS path
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    
    # if extra output locations were configured, note that only the first is used
    if len(outputlocations) > 1:
        currentAlternative.addComputeMessage("Found more than 1 output datapath locations. Using the first, {0}".format(outputlocations[0]))
    outputpath = DSS_Tools.fixInputLocationFpart(currentAlternative, str(outputpath))
#    tspath =str(outputpath)
#    tspath = tspath.split('/')
#    fpart = tspath[6]
#    new_fpart = 'PostProcessed'
#    tspath[6] = new_fpart
#    outputpath = '/'.join(tspath)
    currentAlternative.addComputeMessage("Outputting to {0}".format(outputpath))

    cfs_limit = 0.001 #float
#    currentAlternative.addComputeMessage("\n##### PERFORMING FLOW WEIGHTED AVERAGE FOR UNDER 1CFS #####")
    flowweightaverage.FWA2(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit, 10.0)

#    currentAlternative.addComputeMessage("\n##### REPLACING VALUES OVER 1CFS #####")
    threshold = 1.0
    # find and resolve the DSS paths for the structure flow and total temperature locations
    for location in currentAlternative.getInputDataLocations():
        if str(location) == structure1Flow_location:
            structure1Flow_dsspath = str(currentAlternative.loadTimeSeries(location))
            structure1Flow_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, structure1Flow_dsspath)
        elif str(location) == totalTemperature_location:
            totalTemperature_dsspath = str(currentAlternative.loadTimeSeries(location))
            totalTemperature_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, totalTemperature_dsspath)
    
#    BoundaryFixes.replaceValuesOverThresh(currentAlternative, dss_file, rtw, structure1Flow_dsspath, totalTemperature_dsspath, outputpath, threshold)
    
    currentAlternative.addComputeMessage("\n##### REPLACING REMAINING VALUES #####")
    Wd2_temps_location = 'WD2_Temp'
    # find and resolve the DSS path for the WD2 temperature location
    for location in currentAlternative.getInputDataLocations():
        if str(location) == Wd2_temps_location:
            Wd2_temps_dsspath = str(currentAlternative.loadTimeSeries(location))
            Wd2_temps_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, Wd2_temps_dsspath)
#    BoundaryFixes.replaceNaNValues(currentAlternative, dss_file, rtw, outputpath, Wd2_temps_dsspath)
    return True
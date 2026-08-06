"""
Trinity2Lewiston
=================

Compute a WAT Scripting Alternative that performs a three-stage
post-processing workflow on boundary condition time series for the
Trinity-to-Lewiston link: flow-weighted averaging under a low-flow
threshold, replacement of values over that threshold with a dedicated
source record, and replacement of any remaining NaN values with a
fallback record.

Notes
-----
- **Environment:** Jython (Python 2.7 semantics) within HEC-WAT.
- **Dependencies:** `flowweightaverage`, `BoundaryFixes`, `DSS_Tools`.
"""

import sys                       # Standard library: retained (not directly used beyond potential diagnostics)
import flowweightaverage         # Local module: flow-weighted average temperature computation (FWA)
reload(flowweightaverage)        # Jython: ensure latest version is loaded
import BoundaryFixes             # Local module: threshold/NaN replacement helpers
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
    Compute the Trinity to Lewiston scripting alternative.

    This function performs a three stage post-processing workflow on
    boundary condition time series data:

    1. Filters out the structure 1 flow and total temperature input
       locations, then organizes the remaining input locations and
       computes a flow weighted average for values under 1 CFS using
       `flowweightaverage.FWA`.
    2. Replaces output values over the 1 CFS threshold using
       `BoundaryFixes.replaceValuesOverThresh`, based on the
       structure 1 flow and total temperature input data.
    3. Replaces any remaining NaN values in the output time series
       using `BoundaryFixes.replaceNaNValues`, based on the WD2
       temperature input data.

    Parameters
    ----------
    currentAlternative : object
        The `ScriptingAlternative` being computed. Type:
        `hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt`.
    computeOptions : object
        The compute options for this run. Type:
        `hec.wat.model.ComputeOptions`.

    Returns
    -------
    bool
        True if the script completed successfully. An implicit
        return (no explicit return statement) is treated as a
        successful return.

    Notes
    -----
    This function is structurally very similar to
    `OutputLink_W2-Trinity-W2-Lewiston_SacTrn.py`'s
    `computeAlternative`, but here the full three-stage workflow
    (FWA, threshold replacement, and NaN filling) actually executes,
    since there is no early `return` statement preceding it.
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    
    structure1Flow_location = 'STR1_Flow' #used later to filter
    totalTemperature_location = 'Total_Temp' #used later to filter
    locations = currentAlternative.getInputDataLocations()

    filtered_locations = []
    # loop through all input data locations so the structure 1 flow and
    # total temperature locations can be excluded from this list
    for location in locations:
        # keep only the locations that are not the structure 1 flow or
        # total temperature locations, since those are handled separately below
        if str(location) not in [structure1Flow_location, totalTemperature_location]: #filter out the later locations
            filtered_locations.append(location)
            
    locations = flowweightaverage.organizeLocations(currentAlternative, filtered_locations)

    currentAlternative.addComputeMessage('Found DSS paths:')
    # loop through each organized location group to log its DSS paths
    for location in locations:
        # loop through each path within the current location group
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    # if more than one output location was found, warn that only the
    # first one will be used
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

    cfs_limit = 1.0 #float
    currentAlternative.addComputeMessage("\n##### PERFORMING FLOW WEIGHTED AVERAGE FOR UNDER 1CFS #####")
    flowweightaverage.FWA(currentAlternative, dss_file, rtw, locations, outputpath, cfs_limit)

    currentAlternative.addComputeMessage("\n##### REPLACING VALUES OVER 1CFS #####")
    threshold = 1.0
    # loop through all input data locations to find and load the
    # structure 1 flow and total temperature time series
    for location in currentAlternative.getInputDataLocations():
        # if this location is the structure 1 flow location, load its
        # time series and fix its DSS F-part
        if str(location) == structure1Flow_location:
            structure1Flow_dsspath = str(currentAlternative.loadTimeSeries(location))
            structure1Flow_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, structure1Flow_dsspath)
        # otherwise, if this location is the total temperature location,
        # load its time series and fix its DSS F-part
        elif str(location) == totalTemperature_location:
            totalTemperature_dsspath = str(currentAlternative.loadTimeSeries(location))
            totalTemperature_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, totalTemperature_dsspath)
    
    BoundaryFixes.replaceValuesOverThresh(currentAlternative, dss_file, rtw, structure1Flow_dsspath, totalTemperature_dsspath, outputpath, threshold)
    
    currentAlternative.addComputeMessage("\n##### REPLACING REMAINING VALUES #####")
    Wd2_temps_location = 'WD2_Temp'
    # loop through all input data locations to find and load the
    # WD2 temperature time series
    for location in currentAlternative.getInputDataLocations():
        # if this location is the WD2 temperature location, load its
        # time series and fix its DSS F-part
        if str(location) == Wd2_temps_location:
            Wd2_temps_dsspath = str(currentAlternative.loadTimeSeries(location))
            Wd2_temps_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, Wd2_temps_dsspath)
    BoundaryFixes.replaceNaNValues(currentAlternative, dss_file, rtw, outputpath, Wd2_temps_dsspath)
    return True
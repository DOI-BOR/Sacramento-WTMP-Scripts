import sys
import flowweightaverage
reload(flowweightaverage)
import BoundaryFixes
reload(BoundaryFixes)
import DSS_Tools
reload(DSS_Tools)

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
    
    structure1Flow_location = 'STR1_Flow' #used later to filter
    totalTemperature_location = 'Total_Temp' #used later to filter
    locations = currentAlternative.getInputDataLocations()

    filtered_locations = []
    for location in locations:
        if str(location) not in [structure1Flow_location, totalTemperature_location]: #filter out the later locations
            filtered_locations.append(location)
            
    locations = flowweightaverage.organizeLocations(currentAlternative, filtered_locations)

    currentAlternative.addComputeMessage('Found DSS paths:')
    for location in locations:
        for path in location:
            currentAlternative.addComputeMessage(str(path))
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    outputlocations = currentAlternative.getOutputDataLocations()
    outputpath = currentAlternative.createOutputTimeSeries(outputlocations[0])
    
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
    for location in currentAlternative.getInputDataLocations():
        if str(location) == structure1Flow_location:
            structure1Flow_dsspath = str(currentAlternative.loadTimeSeries(location))
            structure1Flow_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, structure1Flow_dsspath)
        elif str(location) == totalTemperature_location:
            totalTemperature_dsspath = str(currentAlternative.loadTimeSeries(location))
            totalTemperature_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, totalTemperature_dsspath)
    
    BoundaryFixes.replaceValuesOverThresh(currentAlternative, dss_file, rtw, structure1Flow_dsspath, totalTemperature_dsspath, outputpath, threshold)
    
    currentAlternative.addComputeMessage("\n##### REPLACING REMAINING VALUES #####")
    Wd2_temps_location = 'WD2_Temp'
    for location in currentAlternative.getInputDataLocations():
        if str(location) == Wd2_temps_location:
            Wd2_temps_dsspath = str(currentAlternative.loadTimeSeries(location))
            Wd2_temps_dsspath = DSS_Tools.fixInputLocationFpart(currentAlternative, Wd2_temps_dsspath)
    BoundaryFixes.replaceNaNValues(currentAlternative, dss_file, rtw, outputpath, Wd2_temps_dsspath)
    return True


            


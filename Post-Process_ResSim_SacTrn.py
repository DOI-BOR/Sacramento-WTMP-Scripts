import random
from com.rma.io import DssFileManagerImpl

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
	write_example_dl(currentAlternative)
	return True


##
#
# Output data locations can be retrieved from the alternative.
# The method createOutputTimeSeries can be used to create the associated
# TimeSeriesContainer for each output data location.
# Changes to the returned TimeSeriesContainer must be written to disk.
#
##
def write_example_dl(currentAlternative):
	dfm = DssFileManagerImpl.getDssFileManager()
	odls = currentAlternative.getOutputDataLocations()
	val = 0
	for odl in odls:
		val = val + 1
		tsc = currentAlternative.createOutputTimeSeries(odl)
		# currentAlternative.addComputeMessage("tsc.fileName:" + tsc.fileName + " fullName:" + tsc.fullName)
		end = len(tsc.values)
		for i in range(0,end):
			tsc.values[i] = i+val
		dfm.write(tsc)


##
#
# computeOutputVariable function is called when Scripting Output Variables use the 'Other' statistic or
# when variables are not linked to output data locations.
# The computed value needs to be passed to the variable via currentVariable.setValue().
# The computeOutputVariable function does not need to be defined if output variables are not used
# or if all output variables are linked to output data locations and do not use "Other".
#
##
def computeOutputVariable(currentAlternative, currentVariable):
	demo_val = random.randint(10, 100)
	currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName()
		+ " Variable:'" + currentVariable.getName() + "' setValue(" + str(demo_val) + ")")
	currentVariable.setValue(demo_val)
	return True



from hec.heclib.dss import HecDss
from hec.hecmath import HecMathException
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.heclib.util import HecTime
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
import hec.hecmath.TimeSeriesMath as tsmath
from rma.util.RMAConst import MISSING_DOUBLE
import math
import sys
import datetime as dt
import os, sys

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project

sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import Forecast_preprocess as fpp
reload(fpp)

def computeAlternative(currentAlternative, computeOptions):
    """
    Compute a forecast scripting alternative that distributes total
    powerhouse flow across the Trinity/Shasta TCD (Temperature
    Control Device) gate levels and point sinks, based on gate
    status, water surface elevation (WSE), and generation flow.

    NOTE: This function currently only sets up static reference
    data (outlet order and sink elevations) and contains a detailed
    plan, written as comments, for the flow-distribution algorithm.
    The actual computation steps (loading DSS data, computing
    gate/flow/WSE logic, and returning a result) are not yet
    implemented below. This docstring documents both what is
    implemented and the intended behavior described in those
    planning comments, so the gap is clear.

    Intended workflow (per in-line comments, not yet coded):
      1. Load the TCD gate change summary from DSS. Resample to
         regular hourly values if needed, and fill gaps.
      2. Create summed totals of the number of open gates at the
         upper, mid, lower, and side gate levels.
      3. Load WSE (water surface elevation) data, converting units
         to feet if necessary.
      4. Load generation flow data, resampled to hourly. Only
         perform calculations for timestamps where gate, flow, and
         WSE data are all available.
      5. Initialize flow at each of the 4 gate levels to 0.
      6. Distribute flow across levels and point sinks:
         - Assign all flow to the highest open gate level.
         - Split that level's flow evenly across its 3 point sinks
           (while all flow is still nominally on the highest open
           level).
         - Compare WSE minus 3 ft against each sink's elevation:
           if the sink is submerged, record its flow, then check
           if remaining flow is greater than zero.
         - Remove any sinks above the WSE from further
           consideration.
         - If any submerged sinks remain in the current level,
           assign all flow to those remaining 1 or 2 sinks.
         - Otherwise, move to the next gate level down and
           determine how many of ITS point sinks are submerged,
           using a lookahead check against the following level's
           first sink (index into first_next_gate).

    Args:
        currentAlternative: The alternative object being computed.
            Must support addComputeMessage() for logging status
            messages during the compute process.
        computeOptions: The compute options/settings object, used
            to control preprocessing behavior for this run.

    Returns:
        None: No return value is currently implemented. Based on
            the pattern used elsewhere in this codebase, this
            function is expected to eventually return True on
            success (and possibly False on data validation
            failure), once the flow-distribution logic is written.
    """
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')


    outlet_name_order = ['TCDU1','TCDU2','TCDU3','TCDM1','TCDM2','TCDM3','TCDL1','TCDL2','TCDL3','TCDS1','TCDS2','TCDS3']

    # W2 TCD point sink elevations - hard code those here
    sink_elevs = {'TCDU':{1:1042.0,
                          2:1021.0,
                          3:1000.0}
                  'TCDM':{1:942.0,
                          2:921.0,
                          3:900.0}
                  'TCDL':{1:830.0,
                          2:816.0,
                          3:802.0}
                  'TCDS':{1:800.0,
                          2:760.0,
                          3:720.0}
                   }

    # load TCD change summary from DSS, make regular hourly if needed and fill


    # create summed total number gates at upper, mid, lower, and side gates


    # get WSE data, convert to ft if neccessary



    # get generation flow data, also need hourly, only make calcs where we
    # have both gate and flow data and WSE data




    # init flow at each of 4 levels to 0

    # 1) assign all flow to highest gate level open

    # 2) divide flow at each level by 3 and assign as the point source
    # sinks for each level. At this point, all flow is still on the highest
    # open level

    # 3) - If WSE-3 (ft) is less than the withdraw sink, record the flow.
    #     - Then test if that flow is greater than zero
    #	 - Find all the withdraw points with flow again above WSE, and drop them from
    #	 the list of withdraw points with flow.
    #	 - If there are still withdraw points (i.e. the whole gate is not out of the
    #	   water, assing all flow to remaining 1 or 2 withdraw points
    #	 - Else
    #	   go to next level, and figure out how many withdraw points in the NEXT
    #	   active level are below water, using a funny test of the 1st member of the 
    #	   NEXT NEXT level (4th index of first_next_gate, which should be first_next_sink)
    #      *** should the index to first_next_gate actually be 3?



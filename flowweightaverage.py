"""
Flow-Weighted Average Temperature Utilities 

This module provides helper functions to organize paired flow/temperature
DSS pathnames and to compute flow-weighted average (FWA) temperatures
over a specified time window. It is used within HEC-WAT / WTMP workflows
and relies on HEC-DSS APIs available in Jython.

**What’s included**
- Pair-organization for locations linked to prior models.
- Unit-normalization helpers for flow (CFS/CMS) and temperature (C/F).
- Two FWA implementations:
  * `FWA2`: Robust, pair-wise accumulation with data filtering and optional
    backfill/override behavior.
  * `FWA`: Original approach using per-record dictionaries and aggregation.


Notes
-----
- This file expects HEC-WAT/Jython runtime with HEC-DSS libraries on classpath.
"""

from hec.heclib.dss import HecDss  # HEC-DSS file manager (open/read/write/rename records)
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE  # Heclib sentinel value representing undefined doubles
from hec.io import DSSIdentifier  # Optional helper for DSS pathname parsing/building (A–F parts)
from hec.io import TimeSeriesContainer  # Container class holding times/values/units/type for a time series
from rma.util.RMAConst import MISSING_DOUBLE  # RMA/HEC sentinel value representing missing doubles
import math, sys  # Standard Python: math functions and system utilities (e.g., exit)

import DSS_Tools  # Project-local module with DSS helper functions (e.g., fixInputLocationFpart)
reload(DSS_Tools)  # Ensure we have the latest version of DSS_Tools in Jython during runtime


def organizeLocations(currentAlternative, locations):
    """
    Resolve a flat list of input data locations into paired
    (flow, temperature) DSS path pairs.

    Locations are assumed to alternate flow, temperature, flow,
    temperature, and so on, with each flow location immediately
    followed by its matching temperature location. Each location is
    resolved to its actual DSS path and has its F-part corrected.

    Args:
        currentAlternative: The alternative object being computed.
            Must support addComputeMessage() and loadTimeSeries()
            for logging and resolving linked data locations.
        locations (list): Flat list of input data locations,
            alternating flow and temperature entries. Must contain
            an even number of entries.

    Returns:
        list of list: A list of [flow_path, temp_path] pairs, one
            per two input locations, in the order they were given.
            Exits the process (sys.exit(1)) if an odd number of
            locations is provided.
    """
    locations_list = []
    # a flow/temperature pairing requires an even number of locations
    if len(locations) % 2 != 0:
        currentAlternative.addComputeMessage("Uneven amount of Flow/Temp pairings. Check inputs.")
        sys.exit(1)
    # for each location, resolve its DSS path and group it into a flow/temp pair
    for li, location in enumerate(locations):
        tspath =str(currentAlternative.loadTimeSeries(location))
        tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)
        # even index starts a new pair (flow); odd index completes it (temperature)
        if li % 2 == 0: 
            current_pair = [tspath]
        else:
            current_pair.append(tspath)  # complete the pair with temp path
            locations_list.append(current_pair)  # collect finalized pair

    # Return the pair to the calling function
    return locations_list 


def flow_in_cfs(units, flows):
    """
    Normalize flow values to CFS.

    Parameters
    ----------
    units : str
        Units string from DSS (e.g., 'CFS', 'CMS').
    flows : list of float
        Flow values in the given units.

    Returns
    -------
    list of float
        Flow values converted to CFS (or unchanged if already in CFS).

def flow_in_cfs(units,flows):
    """
    Convert a list of flow values to cubic feet per second (cfs),
    if not already in that unit.

    Args:
        units (str): The unit label of the input flow values.
            Recognized values (case-insensitive) are 'cfs' and
            'cms'.
        flows (list of float): Flow values to convert.

    Returns:
        list of float: Flow values in cfs. If units is already
            'cfs', the original list is returned unchanged. If
            units is unrecognized, the process exits
            (sys.exit(-1)) after printing an error.
    """
    if units.lower()=='cfs':
        return flows
    elif units.lower()=='cms':
        values_converted = []

        # Loop and convert values 
        for f in flows:
            values_converted.append(f * 35.314666213)  # conversion: cms → cfs
    
        # Return the converted values
        return values_converted
    
    else:
        # The unit values aren't understood. Fail because that's safe.
        print('FWA2: flow units not known:', units)  # diagnostic
        sys.exit(-1)  # abort on unknown units

def temperature_in_C(units,temps):
    """
    Convert a list of temperature values to degrees Celsius, if not
    already in that unit.

    Args:
        units (str): The unit label of the input temperature
            values. Recognized values (case-insensitive) are 'c',
            'deg c', 'f', and 'deg f'.
        temps (list of float): Temperature values to convert.

    Returns:
        list of float: Temperature values in degrees Celsius. If
            units is already Celsius, the original list is returned
            unchanged. If units is unrecognized, the process exits
            (sys.exit(-1)) after printing an error.
    """
    if units.lower()=='c' or units.lower()=='deg c':
        return temps
    elif units.lower()=='f' or units.lower()=='deg f':
        values_converted = []
        
        # Loop and convert the values
        for t in temps:
            values_converted.append((t - 32.0) * 5.0 / 9.0)  # conversion: °F → °C
        
        # Return the converted values
        return values_converted
    
    else:
        # The unit values aren't understood. Fail because that's safe.
        print('FWA2: temperature units not known:', units)  # diagnostic
        sys.exit(-1)  # abort on unknown units


def FWA2(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None, bad_data_fill_tempC=None, last_override=False):
    """
    Compute a flow-weighted average temperature across multiple
    flow/temperature record pairs, with optional filtering of low
    flows and optional override behavior using the last record
    pair's values directly when its flow exceeds the cfs_limit.

    This is a rewritten version of FWA, intended to fix odd
    behavior observed in the original implementation.

    Workflow:
      1. Opens the DSS file and reads each flow/temperature pair
         from DSSPaths_list.
      2. Converts each pair's units to cfs and Celsius.
      3. For each timestep, if both flow and temperature pass data
         quality checks (not NaN, flow above cfs_limit and below a
         sanity ceiling, temperature within 0-80 C), accumulates
         flow and flow*temperature totals to compute the weighted
         average.
      4. After processing all pairs, computes the flow-weighted
         average temperature at each timestep, using fill_value
         where no valid pairs contributed.
      5. If last_override is True and the last record pair was read
         successfully, overrides the weighted average at any
         timestep where the last pair's own flow/temperature values
         pass the same data quality checks, using its temperature
         value directly instead of the weighted average.
      6. Writes the resulting time series to outputname in dssFile,
         reusing the last temperature record's container.

    NOTE: If reading the very first record pair (dspi == 0) fails
    inside the try/except block, nrecs and temp_type will not have
    been set, and the function will raise a NameError further down.
    This edge case is not currently handled.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        dssFile (str): Path to the DSS file to read from and write
            to.
        timewindow: The run time window object, used to get the
            start and end time strings for reading input data.
        DSSPaths_list (list of list): A list of [flow_path,
            temp_path] pairs, as produced by organizeLocations.
        outputname (str): DSS path to write the resulting
            flow-weighted average temperature record to.
        cfs_limit (float, optional): Minimum flow (in cfs) required
            for a flow/temperature pair to be included in the
            average. Defaults to 0.0 if not provided.
        bad_data_fill_tempC (float, optional): Value to use at
            timesteps where no valid flow/temperature pairs are
            available. Defaults to UNDEFINED_DOUBLE if not provided.
        last_override (bool, optional): If True, use the last
            record pair's own temperature value directly (instead
            of the weighted average) at timesteps where its flow
            exceeds cfs_limit and its temperature is valid. Defaults
            to False.

    Returns:
        int: Always returns 0 on completion.
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)

    Parameters
    ----------
    currentAlt : object
        Alternative providing start/end times and logging.
    dssFile : str
        Path to DSS file for reading and writing.
    timewindow : object
        Provides `getStartTimeString()` / `getEndTimeString()` (HEC format).
    DSSPaths_list : list of [str, str]
        List of DSS path pairs: `[flow_path, temp_path]`.
    outputname : str
        Destination DSS pathname to write the FWA temperature series.
    cfs_limit : float or None, optional
        Minimum flow threshold (CFS); flows below are excluded from weighting.
    bad_data_fill_tempC : float or None, optional
        Fill value (°C) used when no valid pairs exist at a timestep; if None, uses `UNDEFINED_DOUBLE`.
    last_override : bool, optional
        If True and the last record is valid, override FWA at each index with raw temp
        when the concurrent flow/temperature pass all filters.

    Returns
    -------
    int
        Returns 0 after writing the output series.

    Notes
    -----
    - Enforces multiple data checks: non-NaN values, flow within [cfs_limit, 9e6],
      and temperature within [0, 80] °C.
    - Uses temperature type from the first temp record for output `tsc.type`.
    """

    # Get the times for the analysis
    starttime_str = timewindow.getStartTimeString()  # start time (HEC string)
    endtime_str = timewindow.getEndTimeString()      # end time (HEC string)
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))  # log window
    
    # Read each flow/temperature pair and accumulate weighted totals 
    # for each flow/temperature record pair, read and validate the data, then
    # accumulate flow and flow*temperature totals at each timestep
    for dspi, dsspaths in enumerate(DSSPaths_list):
        # Unpack the provided paths
        flow_dss_path = dsspaths[0]  # flow record path
        temp_dss_path = dsspaths[1]  # temperature record path
        currentAlt.addComputeMessage(str(flow_dss_path))  # log flow path
        print('FWA2 Reading:', flow_dss_path)  # diagnostic

        # Read and convert the flow information
        tsc_flow = dssFm.read(flow_dss_path, starttime_str, endtime_str, False).getData()  # read flow data
        flows = flow_in_cfs(tsc_flow.units, tsc_flow.values)  # normalize to CFS

        print('FWA2 Reading:', temp_dss_path)  # diagnostic
        last_rec_valid = False  # test to see if we can use override values

        # Attempt to process the temperature information
        try:
	        tsc_temp = dssFm.read(temp_dss_path, starttime_str, endtime_str, False).getData()
	        temps = temperature_in_C(tsc_temp.units,tsc_temp.values)
	        print('tscf',tsc_flow.values[0])
	        print(flows[0])
	        print('tsct',tsc_temp.values[0])
	        print(temps[0])
	
	        # use type of 1st temp record
	        if dspi==0:
	            nrecs = len(flows)
	            temp_type = tsc_temp.type
	
	        # all pairs must have the same number of records to align by timestep
	        if len(flows) != nrecs or len(temps) != nrecs:
	            currentAlt.addComputeMessage("FWA2: record lengths do not match!")
	            print("FWA2: record lengths do not match!",nrecs,len(flows),len(temps))
	            sys.exit(-1)
	
	        # for each timestep, validate the flow/temp values and accumulate weighted totals
	        for i in range(nrecs):
	            if dspi==0:
	                n_pairs.append(0) # init counter for number of flow/temp pairs in weighted average
	                flow_total.append(0.0)
	                flowtemp_total.append(0.0)
	            # perform a lot of checks on data
	            #print(i,flows[i],temps[i])
	            if not math.isnan(flows[i]) and not math.isnan(temps[i]):
	                if flows[i] > flow_limit and flows[i] < 9.0e6: # could lower upper limit to something relevant to watershed
	                    if temps[i] >= 0.0 and temps[i] <= 80.0:
	                        # passed the data checks
	                        
	                        n_pairs[i] += 1
	                        flow_total[i] += flows[i]
	                        flowtemp_total[i] += flows[i]*temps[i]
	
	                        #print(dspi,i,n_pairs[i],flows[i],temps[i],flow_total[i],flowtemp_total[i])
	        last_rec_valid = True
        except:
	        currentAlt.addComputeMessage('FWA2: data not addeded for record: '+temp_dss_path)
	        last_rec_valid = False
		
    # Compute the final flow-weighted average at each timestep 
    fwat = []
    print('nrecs:',nrecs)
    # for each timestep, compute the weighted average, or apply the last-pair override if requested
    for i in range(nrecs):
        if n_pairs[i] > 0:
            fwat.append(flowtemp_total[i] / flow_total[i])  # compute FWA when pairs exist
        else:
            fwat.append(fill_value)
        # if requested and the last record pair was read successfully, override with
        # its own value directly whenever it passes the same data quality checks
        if last_override and last_rec_valid:
            if flows[i] > flow_limit and flows[i] < 9.0e6:  # same flow filter
                if temps[i] >= 0.0 and temps[i] <= 80.0:    # same temp filter
                    fwat[i] = temps[i]  # override with raw temp (per request)

        # print(i,fwat[i])

    # use last temp container to write output series back to DSS
    tsc_temp.type = temp_type          # enforce type from first temp record
    tsc_temp.fullName = outputname     # destination pathname
    tsc_temp.values = fwat             # computed FWA series

    dssFm.write(tsc_temp)  # write to DSS
    dssFm.close()          # close file
    return 0               # success code


def FWA(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None):
    """
    Compute a flow-weighted average temperature across multiple
    flow/temperature record pairs.

    This is the original flow-weighted average implementation.
    Superseded by FWA2, which was written to address inconsistent
    behavior observed in this version, but this function is kept
    for reference/backward compatibility.

    Workflow:
      1. Opens the DSS file and reads each flow/temperature pair
         from DSSPaths_list into a per-pair dictionary (dss_data).
      2. Converts flow to cfs if needed, and zeroes out any flow
         values below cfs_limit.
      3. Computes flow*temperature for each pair at each timestep.
      4. Sums flow across all pairs at each timestep to get total
         flow.
      5. Sums flow*temperature across all pairs at each timestep
         (skipping NaN contributions) to get total weighted
         temperature.
      6. Divides total weighted temperature by total flow at each
         timestep to get the flow-weighted average, using
         UNDEFINED_DOUBLE where total flow is zero.
      7. Writes the resulting time series to outputname in dssFile.

    NOTE: This function relies on dss_data.keys()[0] to pick an
    arbitrary reference pair, which assumes a stable key/insertion
    order for the dss_data dictionary. Depending on the Python/
    Jython version in use, this may not be guaranteed.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        dssFile (str): Path to the DSS file to read from and write
            to.
        timewindow: The run time window object, used to get the
            start and end time strings for reading input data.
        DSSPaths_list (list of list): A list of [flow_path,
            temp_path] pairs, as produced by organizeLocations.
        outputname (str): DSS path to write the resulting
            flow-weighted average temperature record to.
        cfs_limit (float, optional): Minimum flow (in cfs) required
            for a flow value to be included; flows below this limit
            are zeroed out rather than excluded. If None, no
            filtering is applied.

    Returns:
        int: Always returns 0 on completion.
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    dss_data = {}
    # Read each flow/temperature pair into a per-pair data dictionary 
    # for each flow/temperature record pair, read the data and store it by pair index
    for dspi, dsspaths in enumerate(DSSPaths_list):
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))
        flowTS = dssFm.read(flow_dss_path, starttime_str, endtime_str, False)

    Parameters
    ----------
    currentAlt : object
        Alternative providing window times and logging.
    dssFile : str
        Path to DSS file for reading and writing.
    timewindow : object
        Provides `getStartTimeString()` / `getEndTimeString()` (HEC format).
    DSSPaths_list : list of [str, str]
        List of DSS path pairs: `[flow_path, temp_path]`.
    outputname : str
        Destination DSS pathname to write the FWA temperature series.
    cfs_limit : float or None, optional
        Minimum flow threshold (CFS); flows below are set to zero.

    Returns
    -------
    int
        Returns 0 after writing the output series.

    Notes
    -----
    - Converts CMS to CFS when needed and sets sub-threshold flows to zero.
    - Aggregates per-record flow and flow*temp, then computes FWA.
    """
    
    # Get the time window of the analysis
    starttime_str = timewindow.getStartTimeString()  # window start string
    endtime_str = timewindow.getEndTimeString()      # window end string
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))  # log
    
    # open DSS file
    dssFm = HecDss.open(dssFile)  
    
    # Define a holder for the data
    dss_data = {}  # dict keyed by record index holding flow/temp arrays

    # Loop and process each set of paths
    for dspi, dsspaths in enumerate(DSSPaths_list):
        # Read the flow series
        flow_dss_path = dsspaths[0]  # path for flow record
        currentAlt.addComputeMessage(str(flow_dss_path))  # log path
        flowTS = dssFm.read(flow_dss_path, starttime_str, endtime_str, False)  # windowed read (math obj)

        # Remove the data from the container
        flowTS = flowTS.getData()
        
        flow_units = flowTS.units
        # if the flow series is in cms, convert it to cfs
        if flow_units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')  # log conversion
            flowvals = []
            # convert each cms value to cfs
            for flow in flowTS.values:
                flowvals.append(flow * 35.314666213)  # cms → cfs
            dss_data[dspi] = {'flow': flowvals}  # start dict with converted flows
        else:
            dss_data[dspi] = {'flow': flowTS.values}  # start dict with raw flows

        
        # if a minimum flow limit was provided, zero out any flow below it
        if cfs_limit != None:
            for fi, flow in enumerate(dss_data[dspi]['flow']):
                if flow < cfs_limit:
                    dss_data[dspi]['flow'][fi] = 0  # set to zero

        temp_dss_path = dsspaths[1]  # path for temp record
        currentAlt.addComputeMessage(str(temp_dss_path))  # log path
        TempTS = dssFm.read(temp_dss_path, starttime_str, endtime_str, False)  # windowed read (math obj)
        TempTS = TempTS.getData()  # unwrap to data container
        tempunits = TempTS.units  # store units for output
        temptype = TempTS.type    # store type for output
        dss_data[dspi]['temp'] = TempTS.values  # attach temperature series

#    print(readabledates)
        
    # Compute flow*temperature for each pair
    for dspi in dss_data.keys():
        flowtemps = []
        offset = 0
        
        # for each timestep, multiply this pair's flow by its temperature
        for i, flow in enumerate(dss_data[dspi]['flow']):
            
            temp = dss_data[dspi]['temp'][i]
            flowtemps.append(flow*temp)
            
        dss_data[dspi]['flowtemp'] = flowtemps
        
    # Sum total flow across all pairs at each timestep
    total_flows = []
    dspi = dss_data.keys()[0]
    # for each timestep, sum this pair's flow with every other pair's flow
    for i, flow in enumerate(dss_data[dspi]['flow']):
        temptotalflow = flow
        # add in the flow from every other pair at this same timestep
        for key in dss_data.keys():
            if key != dspi:
                temptotalflow += dss_data[key]['flow'][i]
        total_flows.append(temptotalflow)
        
    # Sum total flow*temperature across all pairs at each timestep 
    total_flowtemp = []
    dspi = dss_data.keys()[0]
    # for each timestep, sum this pair's flow*temperature with every other pair's, skipping NaNs
    for i, flowtemp in enumerate(dss_data[dspi]['flowtemp']):
        temptotalflowtemp = flowtemp
        # add in the flow*temperature from every other pair at this same timestep
        for key in dss_data.keys():
            if key != dspi:
                if not math.isnan(dss_data[key]['flowtemp'][i]):
                    if math.isnan(temptotalflowtemp):
                        temptotalflowtemp = dss_data[key]['flowtemp'][i]  # replace NaN with other record
                    else:
                        temptotalflowtemp += dss_data[key]['flowtemp'][i]
        total_flowtemp.append(temptotalflowtemp)
    
    # Divide to get the final flow-weighted average at each timestep 
    FW_Avg_vals = []
    # for each timestep, divide total flow*temperature by total flow to get the weighted average
    for i, flow in enumerate(total_flows):
        flowtemp = total_flowtemp[i]  # paired sum of flow*temp
        if flow == 0:
            FW_Avg_vals.append(UNDEFINED_DOUBLE)  # sentinel when no flow
        else:
            FW_Avg_vals.append(flowtemp / flow)
    
    # Package the result and write it to DSS 
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputname
    tsc.values = FW_Avg_vals
    tsc.startTime = hecstarttimes[0]
    tsc.units = tempunits
    tsc.type = temptype
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(FW_Avg_vals)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))
    return 0
   
           

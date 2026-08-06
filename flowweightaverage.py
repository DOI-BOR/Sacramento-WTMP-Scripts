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
    Organize paired Flow/Temp locations into DSS path pairs.

    Parameters
    ----------
    currentAlternative : object
        HEC-WAT Alternative providing `loadTimeSeries(location)` and logging.
    locations : list
        List of location objects; expected in alternating Flow/Temp order.

    Returns
    -------
    list of [str, str]
        List of DSS pathname pairs: `[flow_path, temp_path]`.

    Raises
    ------
    SystemExit
        If the number of locations is odd (uneven Flow/Temp pairings).

    Notes
    -----
    - Uses `DSS_Tools.fixInputLocationFpart` to align input F-part
      to the Alternative’s expected F-part structure.
    """
    
    # Define a list to hold output
    locations_list = []
    
    # Confirm that there are an event number of series
    if len(locations) % 2 != 0:
        currentAlternative.addComputeMessage("Uneven amount of Flow/Temp pairings. Check inputs.")
        sys.exit(1)  # abort if not paired

    # Loop and process each location
    for li, location in enumerate(locations):
        tspath = str(currentAlternative.loadTimeSeries(location))  # derive linked DSS path
        tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)  # adjust F-part
        
        if li % 2 == 0:
            current_pair = [tspath]  # start a new flow/temp pair
    
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

    Raises
    ------
    SystemExit
        If units are not recognized ('cfs' or 'cms').
    """
    
    # Determine if values need to be converted
    if units.lower() == 'cfs':
        # If the units are in cfs, return them unchanged
        return flows  
    
    elif units.lower() == 'cms':
        # Units are in cms. Make the conversion.
        # Define a holder for the values
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


def temperature_in_C(units, temps):
    """
    Normalize temperature values to degrees Celsius.

    Parameters
    ----------
    units : str
        Units string from DSS (e.g., 'C', 'F', 'DEG C', 'DEG F').
    temps : list of float
        Temperature values in the given units.

    Returns
    -------
    list of float
        Temperature values converted to Celsius (or unchanged if already °C).

    Raises
    ------
    SystemExit
        If units are not recognized variants of Celsius or Fahrenheit.
    """
    
    # Determine if the units need to be converted
    if units.lower() == 'c' or units.lower() == 'deg c':
        # Do not convert since they're already in C
        return temps  # already in °C
    
    elif units.lower() == 'f' or units.lower() == 'deg f':
        # Convert from F to C
        # Define a holder for the values
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


def FWA2(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None,
         bad_data_fill_tempC=None, last_override=False):
    """
    Compute flow-weighted average temperature using pair-wise accumulation and validation.

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
    
    # Open the DSS file
    dssFm = HecDss.open(dssFile)  

    # Create holders for the values
    flow_total = []          # per-time-step sum of valid flows
    flowtemp_total = []      # per-time-step sum of flow*temp across valid pairs
    n_pairs = []             # per-time-step count of valid flow/temp pairs

    # Set placeholder values
    flow_limit = 0.0 if cfs_limit is None else cfs_limit  # default threshold
    fill_value = UNDEFINED_DOUBLE if bad_data_fill_tempC is None else bad_data_fill_tempC  # default fill

    # Loop and process each of the DSS paths
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
            # Read the temperature series
            tsc_temp = dssFm.read(temp_dss_path, starttime_str, endtime_str, False).getData()  # read temp data
            temps = temperature_in_C(tsc_temp.units, tsc_temp.values)  # normalize to °C
            print('tscf', tsc_flow.values[0])  # debug print of first flow
            print(flows[0])                    # normalized first flow
            print('tsct', tsc_temp.values[0])  # debug print of first temp
            print(temps[0])                    # normalized first temp
    
            # use type of 1st temp record
            if dspi == 0:
                nrecs = len(flows)   # store expected length from first flow series
                temp_type = tsc_temp.type  # capture type for output
    
            if len(flows) != nrecs or len(temps) != nrecs:
                currentAlt.addComputeMessage("FWA2: record lengths do not match!")
                print("FWA2: record lengths do not match!", nrecs, len(flows), len(temps))
                sys.exit(-1)  # enforce equal length per timestep
    
            for i in range(nrecs):
                if dspi == 0:
                    n_pairs.append(0)          # init counter for number of flow/temp pairs in weighted average
                    flow_total.append(0.0)     # init total flow per index
                    flowtemp_total.append(0.0) # init total flow*temp per index
                
                # perform a lot of checks on data
                if not math.isnan(flows[i]) and not math.isnan(temps[i]):
                    if flows[i] > flow_limit and flows[i] < 9.0e6:  # could lower upper limit to something relevant to watershed
                        if temps[i] >= 0.0 and temps[i] <= 80.0:    # plausible temperature bounds
                            # passed the data checks
                            n_pairs[i] += 1
                            flow_total[i] += flows[i]
                            flowtemp_total[i] += flows[i] * temps[i]
    
                            # print(dspi,i,n_pairs[i],flows[i],temps[i],flow_total[i],flowtemp_total[i])
            last_rec_valid = True  # mark last read successful
        
        except:
            currentAlt.addComputeMessage('FWA2: data not addeded for record: ' + temp_dss_path)
            last_rec_valid = False  # mark read failure
    
    # output FWA temperature series
    fwat = []  

    for i in range(nrecs):
        if n_pairs[i] > 0:
            fwat.append(flowtemp_total[i] / flow_total[i])  # compute FWA when pairs exist
        else:
            fwat.append(fill_value)  # use fill/sentinel when no valid pairs

        # Optional override using last record’s raw temperature (if valid and enabled)
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
    Original flow-weighted average temperature computation using dictionaries.

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
        
        # Get the time series
        hecstarttimes = flowTS.times  # store time array for output series alignment

        # Convert the flow units, if necessary
        flow_units = flowTS.units  # input units
        if flow_units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')  # log conversion
            flowvals = []
            for flow in flowTS.values:
                flowvals.append(flow * 35.314666213)  # cms → cfs
            dss_data[dspi] = {'flow': flowvals}  # start dict with converted flows
        else:
            dss_data[dspi] = {'flow': flowTS.values}  # start dict with raw flows

        # Optionally zero-out flows below threshold to exclude them from weighting
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

    # Calculate the flow weighted temperature
    for dspi in dss_data.keys():
        flowtemps = []  # element-wise flow*temp products
        offset = 0      # retained from original (unused but preserved)
        for i, flow in enumerate(dss_data[dspi]['flow']):
            temp = dss_data[dspi]['temp'][i]  # matching temp value
            flowtemps.append(flow * temp)     # product for weighting
        dss_data[dspi]['flowtemp'] = flowtemps  # store products

    # Sum the flow components
    total_flows = []  # element-wise sum of flows across records
    dspi = dss_data.keys()[0]  # base key for iteration
    for i, flow in enumerate(dss_data[dspi]['flow']):
        temptotalflow = flow  # start with base record's flow
        for key in dss_data.keys():
            if key != dspi:
                temptotalflow += dss_data[key]['flow'][i]  # add other records' flows
        total_flows.append(temptotalflow)  # collect sum per index

    # Perform QA/QC on the flow-temp series
    total_flowtemp = []  # element-wise sum of flow*temp across records
    dspi = dss_data.keys()[0]  # base key again
    for i, flowtemp in enumerate(dss_data[dspi]['flowtemp']):
        temptotalflowtemp = flowtemp  # start with base record's flow*temp
        for key in dss_data.keys():
            if key != dspi:
                if not math.isnan(dss_data[key]['flowtemp'][i]):
                    if math.isnan(temptotalflowtemp):
                        temptotalflowtemp = dss_data[key]['flowtemp'][i]  # replace NaN with other record
                    else:
                        temptotalflowtemp += dss_data[key]['flowtemp'][i]  # sum valid products
        total_flowtemp.append(temptotalflowtemp)  # collect sum per index

    # Calculate the flow weighted temperature
    FW_Avg_vals = []  # final flow-weighted average temperatures
    for i, flow in enumerate(total_flows):
        flowtemp = total_flowtemp[i]  # paired sum of flow*temp
        if flow == 0:
            FW_Avg_vals.append(UNDEFINED_DOUBLE)  # sentinel when no flow
        else:
            FW_Avg_vals.append(flowtemp / flow)   # compute average

    # Define a timeseries for output
    tsc = TimeSeriesContainer()  # build container for output series
    tsc.times = hecstarttimes    # assign original times
    tsc.fullName = outputname    # destination pathname
    tsc.values = FW_Avg_vals     # computed FWA series
    tsc.startTime = hecstarttimes[0]  # start time index
    tsc.units = tempunits        # keep temperature units from input series
    tsc.type = temptype          # keep series type from input series
    tsc.endTime = hecstarttimes[-1]   # end time index
    tsc.numberValues = len(FW_Avg_vals)  # count of values
    tsc.startHecTime = timewindow.getStartTime()  # window (HEC time)
    tsc.endHecTime = timewindow.getEndTime()      # window (HEC time)

    # Write to the file and close
    dssFm.write(tsc)  # write output series to DSS
    dssFm.close()     # close file
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))  # log size

    # Return that writing has been successful
    return 0
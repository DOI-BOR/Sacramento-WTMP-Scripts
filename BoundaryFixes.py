from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE

def replaceValuesOverThresh(currentAlt, dssFile, timewindow, primary_data_dsspath, secondary_data_dsspath, tertiary_data_dsspath, threshold):
    """
    When Primary file under threshold, use established data values from tertiary_data_dsspath record
    when over threshold, use values from secondary data dsspath

    For each timestep where the primary record's value exceeds the
    given threshold, the value in the tertiary ("existing") record
    at that same timestamp is overwritten with the corresponding
    value from the secondary record. Timesteps where the primary
    value stays at or below the threshold are left untouched in the
    tertiary record. The (possibly modified) tertiary record is
    then written back to the same DSS path.
    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        dssFile (str): Path to the DSS file to read all three
            records from and write the result to.
        timewindow: The run time window object, used to get the
            start and end time strings for reading input data.
        primary_data_dsspath (str): DSS path of the record whose
            values are compared against threshold. Converted from
            cms to cfs if necessary.
        secondary_data_dsspath (str): DSS path of the record to
            pull replacement values from, at timesteps where the
            primary value exceeds threshold.
        tertiary_data_dsspath (str): DSS path of the "existing"
            record that gets selectively overwritten and then
            rewritten back to this same path.
        threshold (float): The threshold value (in the primary
            record's units, after cfs conversion if applicable)
            above which replacement occurs.

    Returns:
        int: Always returns 0 on completion. Writes the modified
            tertiary record back to tertiary_data_dsspath in
            dssFile. If a matching timestamp cannot be found in
            either the tertiary or secondary record for a given
            primary timestep over threshold, a message is printed
            (but no error is raised) and that timestep is left
            unmodified.
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()

    # Log the inspection window for traceability in WAT messages
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    # Open the target DSS file containing all series needed for this operation
    dssFm = HecDss.open(dssFile)
    
    # --- Read the primary record and convert units if necessary ---
    PrimaryTS = dssFm.read(primary_data_dsspath, starttime_str, endtime_str, False)
    PrimaryTS = PrimaryTS.getData()  # Extract the internal data container

    # Capture metadata and values from the primary read
    Primary_units = PrimaryTS.units
    PrimaryTS_values = PrimaryTS.values
    PrimaryTS_times = PrimaryTS.times
    # if the primary series is in cms, convert it to cfs before comparing against threshold
    if Primary_units.lower() == 'cms':
        currentAlt.addComputeMessage('Converting cms to cfs')
        PrimaryTS_values = []  # Rebuild values with unit conversion applied
        for val in PrimaryTS.values:
            PrimaryTS_values.append(val * 35.314666213)  # CMS → CFS

    # --- Read the secondary (replacement source) and tertiary (existing/target) records ---
    SecondaryTS = dssFm.read(secondary_data_dsspath, starttime_str, endtime_str, False)
    SecondaryTS = SecondaryTS.getData()
    SecondaryTS_values = SecondaryTS.values
    SecondaryTS_times = SecondaryTS.times

    # ---- Read the existing/tertiary series that will be modified and re-written ----------
    ExistingTS = dssFm.read(tertiary_data_dsspath, starttime_str, endtime_str, False)
    ExistingTS = ExistingTS.getData()
    ExistingTS_values = ExistingTS.values
    ExistingTS_times = ExistingTS.times
    ExistingTS_units = ExistingTS.units
    ExistingTS_type = ExistingTS.type

    # --- Replace values in the existing record wherever the primary exceeds the threshold ---
    # for each primary value, if it exceeds the threshold, substitute the secondary
    # record's value into the existing record at the matching timestamp
    for i, primary_val in enumerate(PrimaryTS_values):
        if primary_val > threshold:
            primarytime = PrimaryTS_times[i]
            # only substitute if the same timestamp exists in both the existing and secondary records
            if primarytime in ExistingTS_times and primarytime in SecondaryTS_times:
                existing_time_index = ExistingTS_times.index(primarytime)   # Find index in existing series
                secondary_time_index = SecondaryTS_times.index(primarytime) # Find index in secondary series
                ExistingTS_values[existing_time_index] = SecondaryTS_values[secondary_time_index]  # Perform replacement
            
            else:
                # report which record(s) are missing this timestamp, and leave the value unmodified
                if primarytime not in ExistingTS_times:
                    print('Unable to find time {0} in {1}'.format(primarytime, tertiary_data_dsspath))
                if primarytime not in SecondaryTS_times:
                    print('Unable to find time {0} in {1}'.format(primarytime, secondary_data_dsspath))

    # --- Write the modified existing record back to DSS ---
    tsc = TimeSeriesContainer()
    tsc.times = ExistingTS_times            # Use the original time vector
    tsc.fullName = tertiary_data_dsspath    # Write back to the tertiary (existing) DSS path
    tsc.values = ExistingTS_values          # Modified values after replacements
    tsc.startTime = ExistingTS_times[0]     # Start time marker for container
    tsc.units = ExistingTS_units            # Preserve original units
    tsc.type = ExistingTS_type              # Preserve original type (e.g., INST-VAL)
    tsc.endTime = ExistingTS_times[-1]      # End time marker for container
    tsc.numberValues = len(ExistingTS_values)  # Count of values

    # Provide HEC-time boundaries based on the input time window for completeness
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()

    # Write the modified series and close the DSS file to release resources
    dssFm.write(tsc)
    dssFm.close()

    # Log the number of values written for transparency in the WAT compute log
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(ExistingTS_values)))
    
    # Return successful
    return 0


def replaceNaNValues(currentAlt, dssFile, timewindow, existing_dsspath, fill_dsspath):
    """
    Fill in undefined (missing) values in an existing time series
    record using matching timestamps from a separate fill record,
    and write the result back to the same DSS path.

    NOTE: The error-message branch below references the variables
    `primarytime` and `tertiary_data_dsspath`, neither of which
    exist in this function (they belong to the similarly-structured
    replaceValuesOverThresh function above). As written, this branch
    will raise a NameError if it is ever reached (i.e. if a missing
    value's timestamp cannot be found in the fill record). This
    appears to be leftover from copying that function's logic and
    has been left unchanged here for visibility.

    Args:
        currentAlt: The alternative object being computed. Must
            support addComputeMessage() for logging.
        dssFile (str): Path to the DSS file to read both records
            from and write the result to.
        timewindow: The run time window object, used to get the
            start and end time strings for reading input data.
        existing_dsspath (str): DSS path of the record to fill
            missing values in, and to write the result back to.
        fill_dsspath (str): DSS path of the record to pull
            replacement values from, at timestamps where the
            existing record has an undefined value.

    Returns:
        int: Always returns 0 on completion. Writes the filled
            record back to existing_dsspath in dssFile.

    Raises:
        NameError: If a missing value's timestamp cannot be found
            in the fill record (see NOTE above).
    """
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    # Open the DSS file holding both target and fill series
    dssFm = HecDss.open(dssFile)

    # --- Read the existing record to be filled ---
    ExistingTS = dssFm.read(existing_dsspath, starttime_str, endtime_str, False)
    ExistingTS = ExistingTS.getData()      # Extract container
    ExistingTS_values = ExistingTS.values  # Values to scan for UNDEFINED_DOUBLE
    ExistingTS_times = ExistingTS.times    # Time vector
    ExistingTS_units = ExistingTS.units    # Preserve units on write
    ExistingTS_type = ExistingTS.type      # Preserve type on write

    # --- Read the fill record to pull replacement values from ---
    FillTS = dssFm.read(fill_dsspath, starttime_str, endtime_str, False)
    FillTS = FillTS.getData()             # Extract container
    FillTS_values = FillTS.values         # Donor values
    FillTS_times = FillTS.times           # Donor time vector
    FillTS_units = FillTS.units           # Units (typically same domain; not used directly)

    # --- Fill any undefined values using the matching timestamp in the fill record ---
    # for each value in the existing record, if it is undefined, look up the
    # matching timestamp in the fill record and substitute that value in
    for i, value in enumerate(ExistingTS_values):
        if value == UNDEFINED_DOUBLE:  # Check for sentinel indicating undefined/missing
            existingtime = ExistingTS_times[i]  # Timestamp corresponding to the missing value

            # Replace only when the timestamp is present in the fill series
            if existingtime in FillTS_times:
                FillTS_time_index = FillTS_times.index(existingtime)  # Locate donor index
                ExistingTS_values[i] = FillTS_values[FillTS_time_index]  # Perform replacement
                # print('Filled at {0}-{1}'.format(existingtime, FillTS_times[FillTS_time_index]))
            else:
                # Diagnostic print if the timestamp is not found in the fill series
                # NOTE: The message below references variables not defined in this function
                # (e.g., primarytime, tertiary_data_dsspath). Left unchanged intentionally.
                print('Unable to find time {0} in {1}'.format(primarytime, tertiary_data_dsspath))
                
    # --- Write the filled record back to DSS ---
    tsc = TimeSeriesContainer()
    tsc.times = ExistingTS_times             # Original time vector
    tsc.fullName = existing_dsspath          # Write back to the existing series path
    tsc.values = ExistingTS_values           # Updated values after replacements
    tsc.startTime = ExistingTS_times[0]      # Start marker
    tsc.units = ExistingTS_units             # Preserve units
    tsc.type = ExistingTS_type               # Preserve type    
    tsc.endTime = ExistingTS_times[-1]       # End marker
    tsc.numberValues = len(ExistingTS_values)  # Number of values

    # Provide HEC-time boundaries from the input time window
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()

    # Write the cleaned series, close the DSS file, and log the operation
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(ExistingTS_values)))
    
    # Return successful
    return 0 
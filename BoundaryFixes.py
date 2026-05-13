"""
DSS Value Replacement Utilities (WTMP / HEC-WAT)
================================================

Purpose
-------
Provide helper functions to modify/write time series in a DSS file for the
Sacramento–Trinity WTMP workflow:

1. ``replaceValuesOverThresh``  
   For a given primary series, when a value is **over** a threshold, replace
   the corresponding value in an existing (tertiary) series with the value
   from a secondary series. When a value is **under** the threshold, retain
   the existing (tertiary) value.

2. ``replaceNaNValues``  
   Replace undefined/missing values (``UNDEFINED_DOUBLE`` sentinel) in an
   existing series using values from a fill series at matching timestamps.

Notes
-----
- **Environment:** Jython within HEC‑WAT; uses HEC‑DSS APIs and containers.
- **Unit handling:** ``replaceValuesOverThresh`` converts primary series
  from CMS → CFS if units are reported as ``'cms'``.
- **Time matching:** Replacement operations require that timestamps in the
  target series exist in the donor series; otherwise a message is printed
  and the value is not replaced.
- **Write behavior:** Results are written back to the DSS file using the
  original path of the modified (existing/tertiary or existing) series.

See Also
--------
hec.heclib.dss.HecDss
    API to open, read, and write DSS time series.
hec.io.TimeSeriesContainer
    Container used to write updated time series back to DSS.
"""

from hec.heclib.dss import HecDss            # HEC-DSS API: open/read/write operations on .dss files
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE  # HEC utility sentinel for undefined/missing values
from hec.io import TimeSeriesContainer        # HEC I/O container used to assemble time series for writing
from rma.util.RMAConst import MISSING_DOUBLE  # RMA constant for missing values (imported for consistency; not used here)


def replaceValuesOverThresh(currentAlt, dssFile, timewindow,
                            primary_data_dsspath, secondary_data_dsspath,
                            tertiary_data_dsspath, threshold):
    """
    Replace values in an existing (tertiary) series when the primary series exceeds a threshold.

    When the primary value at a timestamp is **over** the threshold, the
    corresponding value in the tertiary (existing) series is replaced with
    the value from the secondary series. When the primary value is **under**
    the threshold, the tertiary value remains unchanged.

    Parameters
    ----------
    currentAlt : object
        WAT scripting alternative used for logging (e.g., ``addComputeMessage``).
    dssFile : str
        Path to the DSS file containing all referenced time series.
    timewindow : object
        Runtime window object providing start/end times and HEC time conversions.
    primary_data_dsspath : str
        Full DSS pathname for the primary series (driver for threshold logic).
    secondary_data_dsspath : str
        Full DSS pathname for the secondary (donor) series used when threshold is exceeded.
    tertiary_data_dsspath : str
        Full DSS pathname for the existing (target) series to be modified and re-written.
    threshold : float
        Numeric threshold; replacement occurs when ``primary_val > threshold``.

    Returns
    -------
    int
        ``0`` on completion 

    Notes
    -----
    - If the primary units are ``cms``, values are converted to ``cfs`` using
      factor ``35.314666213`` before threshold comparison.
    - Replacement is performed only when the timestamp exists in **both**
      the tertiary and secondary series; otherwise an informational message
      is printed and the value at that timestamp is left unchanged.
    """
    
    # Extract start/end strings from the runtime window for DSS read operations
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()

    # Log the inspection window for traceability in WAT messages
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    # Open the target DSS file containing all series needed for this operation
    dssFm = HecDss.open(dssFile)
    
    # ---- Read the primary series that drives the threshold logic -------------------------
    PrimaryTS = dssFm.read(primary_data_dsspath, starttime_str, endtime_str, False)
    PrimaryTS = PrimaryTS.getData()  # Extract the internal data container

    # Capture metadata and values from the primary read
    Primary_units = PrimaryTS.units
    PrimaryTS_values = PrimaryTS.values
    PrimaryTS_times = PrimaryTS.times

    # If the primary series is reported in 'cms', convert to 'cfs' for consistent comparison
    if Primary_units.lower() == 'cms':
        currentAlt.addComputeMessage('Converting cms to cfs')
        PrimaryTS_values = []  # Rebuild values with unit conversion applied
        for val in PrimaryTS.values:
            PrimaryTS_values.append(val * 35.314666213)  # CMS → CFS

    # ---- Read the secondary (donor) series used when threshold is exceeded ---------------
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

    # ---- Iterate over primary values, performing replacements when threshold is exceeded -
    for i, primary_val in enumerate(PrimaryTS_values):
        if primary_val &gt; threshold:  # Use HTML entity as provided; replacement triggers when primary exceeds threshold
            primarytime = PrimaryTS_times[i]  # Timestamp at the current index

            # Only replace if the timestamp exists in both the tertiary and secondary series
            if primarytime in ExistingTS_times and primarytime in SecondaryTS_times:
                existing_time_index = ExistingTS_times.index(primarytime)   # Find index in existing series
                secondary_time_index = SecondaryTS_times.index(primarytime) # Find index in secondary series
                ExistingTS_values[existing_time_index] = SecondaryTS_values[secondary_time_index]  # Perform replacement
            
            else:
                # Print diagnostics if either timestamp is missing; values remain unchanged at this time
                if primarytime not in ExistingTS_times:
                    print('Unable to find time {0} in {1}'.format(primarytime, tertiary_data_dsspath))
                if primarytime not in SecondaryTS_times:
                    print('Unable to find time {0} in {1}'.format(primarytime, secondary_data_dsspath))

    # ---- Assemble a TimeSeriesContainer and write the modified tertiary series back ------
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
    Replace undefined values in an existing series using a fill series at matching timestamps.

    This function scans the existing series for values equal to the
    ``UNDEFINED_DOUBLE`` sentinel and replaces those values with the
    corresponding values from a fill series (when the timestamp is present).

    Parameters
    ----------
    currentAlt : object
        WAT scripting alternative used for logging.
    dssFile : str
        Path to the DSS file containing both existing and fill series.
    timewindow : object
        Runtime window object providing start/end time strings and HEC times.
    existing_dsspath : str
        Full DSS pathname for the existing series to be checked/updated.
    fill_dsspath : str
        Full DSS pathname for the fill series providing replacement values.

    Returns
    -------
    int
        ``0`` on completion.

    Notes
    -----
    - Replacement occurs only when the timestamp exists in the fill series.
    - A message is printed if a matching timestamp is not found in the fill series.
    """
    
    # Retrieve human-readable start/end strings for DSS queries and log the window
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))

    # Open the DSS file holding both target and fill series
    dssFm = HecDss.open(dssFile)

    # ---- Read the existing series to be cleaned -----------------------------------------
    ExistingTS = dssFm.read(existing_dsspath, starttime_str, endtime_str, False)
    ExistingTS = ExistingTS.getData()      # Extract container
    ExistingTS_values = ExistingTS.values  # Values to scan for UNDEFINED_DOUBLE
    ExistingTS_times = ExistingTS.times    # Time vector
    ExistingTS_units = ExistingTS.units    # Preserve units on write
    ExistingTS_type = ExistingTS.type      # Preserve type on write

    # ---- Read the fill series providing replacement values ------------------------------
    FillTS = dssFm.read(fill_dsspath, starttime_str, endtime_str, False)
    FillTS = FillTS.getData()             # Extract container
    FillTS_values = FillTS.values         # Donor values
    FillTS_times = FillTS.times           # Donor time vector
    FillTS_units = FillTS.units           # Units (typically same domain; not used directly)

    # ---- Iterate and replace undefined values when a matching timestamp exists ----------
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
                
    # ---- Assemble container and write the updated existing series back to DSS -----------
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
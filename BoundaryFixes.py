from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE

def replaceValuesOverThresh(currentAlt, dssFile, timewindow, primary_data_dsspath, secondary_data_dsspath, tertiary_data_dsspath, threshold):
    '''
    When Primary file under threshold, use established data values from tertiary_data_dsspath record
    when over threshold, use values from secondary data dsspath
    '''
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    
    PrimaryTS = dssFm.read(primary_data_dsspath, starttime_str, endtime_str, False)
    PrimaryTS = PrimaryTS.getData()

    Primary_units = PrimaryTS.units
    PrimaryTS_values = PrimaryTS.values
    PrimaryTS_times = PrimaryTS.times
    if Primary_units.lower() == 'cms':
        currentAlt.addComputeMessage('Converting cms to cfs')
        PrimaryTS_values = []
        for val in PrimaryTS.values:
            PrimaryTS_values.append(val * 35.314666213)

    SecondaryTS = dssFm.read(secondary_data_dsspath, starttime_str, endtime_str, False)
    SecondaryTS = SecondaryTS.getData()
    SecondaryTS_values = SecondaryTS.values
    SecondaryTS_times = SecondaryTS.times

    ExistingTS = dssFm.read(tertiary_data_dsspath, starttime_str, endtime_str, False)
    ExistingTS = ExistingTS.getData()
    ExistingTS_values = ExistingTS.values
    ExistingTS_times = ExistingTS.times
    ExistingTS_units = ExistingTS.units
    ExistingTS_type = ExistingTS.type

    for i, primary_val in enumerate(PrimaryTS_values):
        if primary_val > threshold:
            primarytime = PrimaryTS_times[i]
            if primarytime in ExistingTS_times and primarytime in SecondaryTS_times:
                existing_time_index = ExistingTS_times.index(primarytime)
                secondary_time_index = SecondaryTS_times.index(primarytime)
                ExistingTS_values[existing_time_index] = SecondaryTS_values[secondary_time_index]
            else:
                if primarytime not in ExistingTS_times:
                    print('Unable to find time {0} in {1}'.format(primarytime, tertiary_data_dsspath))
                if primarytime not in SecondaryTS_times:
                    print('Unable to find time {0} in {1}'.format(primarytime, secondary_data_dsspath))

    tsc = TimeSeriesContainer()
    tsc.times = ExistingTS_times
    tsc.fullName = tertiary_data_dsspath
    tsc.values = ExistingTS_values
    tsc.startTime = ExistingTS_times[0]
    tsc.units = ExistingTS_units
    tsc.type = ExistingTS_type
    tsc.endTime = ExistingTS_times[-1]
    tsc.numberValues = len(ExistingTS_values)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(ExistingTS_values)))
    return 0    

def replaceNaNValues(currentAlt, dssFile, timewindow, existing_dsspath, fill_dsspath):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)

    ExistingTS = dssFm.read(existing_dsspath, starttime_str, endtime_str, False)
    ExistingTS = ExistingTS.getData()
    ExistingTS_values = ExistingTS.values
    ExistingTS_times = ExistingTS.times
    ExistingTS_units = ExistingTS.units
    ExistingTS_type = ExistingTS.type

    FillTS = dssFm.read(fill_dsspath, starttime_str, endtime_str, False)
    FillTS = FillTS.getData()
    FillTS_values = FillTS.values
    FillTS_times = FillTS.times
    FillTS_units = FillTS.units

    for i, value in enumerate(ExistingTS_values):
        if value == UNDEFINED_DOUBLE:
            existingtime = ExistingTS_times[i]
            if existingtime in FillTS_times:
                FillTS_time_index = FillTS_times.index(existingtime)
                ExistingTS_values[i] = FillTS_values[FillTS_time_index]
#                print('Filled at {0}-{1}'.format(existingtime, FillTS_times[FillTS_time_index]))
            else:
                print('Unable to find time {0} in {1}'.format(primarytime, tertiary_data_dsspath))
                
    tsc = TimeSeriesContainer()
    tsc.times = ExistingTS_times
    tsc.fullName = existing_dsspath
    tsc.values = ExistingTS_values
    tsc.startTime = ExistingTS_times[0]
    tsc.units = ExistingTS_units
    tsc.type = ExistingTS_type    
    tsc.endTime = ExistingTS_times[-1]
    tsc.numberValues = len(ExistingTS_values)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(ExistingTS_values)))
    return 0    

         
         

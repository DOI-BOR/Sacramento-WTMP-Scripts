
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
#import hec.hecmath.TimeSeriesMath as tsmath
from Simple_DSS_Functions import resample_dss_ts

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

units_need_fixing = ['radians','tenths',r'langley/min']

def fix_DMS_types_units(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    dss = HecDss.open(dss_file)
    recs = dss.getPathnameList()
    for r in recs:
        tsm = dss.read(r)
        rlow = r.lower()
        if "/flow" in rlow or "/1day/" in rlow:
            tsm.setType('PER-AVER')
            dss.write(tsm)
        if tsm.getUnits() in units_need_fixing:
            tsm = standardize_units_tsm(tsm)
            dss.write(tsm)
    dss.close()

def standardize_units_tsm(tsm):
    tsc = tsm.getData()
    tsc = standardize_units_tsc(tsc)
    tsm.setData(tsc)
    return tsm

def standardize_units_tsc(tsc):
    if tsc.units == 'radians':
        for i in range(len(tsc.values)) :
            tsc.values[i] = tsc.values[i] / 2*3.141592653589793 * 360.0
        tsc.units = 'deg'
    if tsc.units == 'tenths':
        # this is not really 'tenths', but in the data from DMS, it descirbes data ranging from 0-10.
        # ResSim expects 0-1 for cloud cover
        for i in range(len(tsc.values)) :
            tsc.values[i] = tsc.values[i] / 10.0
        tsc.units = 'FRAC'
    if tsc.units == r'langley/min':
        conv = 41840.0 / 60.0  # lang/min * 41840 j/m2 * 1 min/60 s  j/m2/s = W/m2
        for i in range(len(tsc.values)) :
            tsc.values[i] = tsc.values[i] * conv
        tsc.units = r'W/m2'
    return tsc


def standardize_interval(tsm, interval, makePerAver=True):
    tsc = tsm.getData()
    if interval.lower()=='1hour':
        intint = 60
    elif interval.lower()=='1day':
        intint=1440
    else:
        print('internal not supported:',interval)
        sys.exit(-1)

    if tsc.interval != intint:
        if makePerAver:
            #tsc.type = 'PER-AVER'  # make sure it's per-aver ... we are
            tsm.setType('PER-AVER')
        return tsm.transformTimeSeries(interval, "", "AVE")
    else:
        return tsm

def replace_data(currentAlt, timewindow, pairs, dss_file, dss_outfile, months, standard_interval=None):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    # 01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    for pair in pairs:
        dssFm = HecDss.open(dss_file)
        currentAlt.addComputeMessage('Replacing data for {0} with {1} during {2}'.format(pair[0], pair[1], months))
        base = dssFm.read(pair[0], starttime_str, endtime_str, False)
        if standard_interval is not None:
            base = standardize_interval(base,standard_interval)
        base_data = base.getData()
        base_values = base_data.values
        base_hectimes = base_data.times
        base_units = base_data.units
        base_interval = base_data.interval
        base_values, base_hectimes = check_start_and_end(base_values, base_hectimes, starttime_hectime, endtime_hectime)

        alt = dssFm.read(pair[1], starttime_str, endtime_str, False)
        if standard_interval is not None:
            alt = standardize_interval(alt,standard_interval)
        alt_data = alt.getData()
        alt_values = alt_data.values
        alt_hectimes = alt_data.times
        alt_units = alt_data.units
        alt_interval = alt_data.interval
        alt_values, alt_hectimes = check_start_and_end(alt_values, alt_hectimes, starttime_hectime, endtime_hectime)

        dssFm.close()

        if base_units != alt_units:
            currentAlt.addComputeMessage('Units do not match for {0} and {1}, skipping'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)
        if base_interval != alt_interval:
            currentAlt.addComputeMessage('Intervals do not match for {0} and {1}, changing interval...'.format(pair[0], pair[1]))
            dssFm.close()
            sys.exit(1)

        for i in range(len(base_values)):
            if base_data.getHecTime(i).month() in months:
                base_values[i] = alt_values[i]
                print('replaced {0}'.format(base_data.getHecTime(i)))

        new_pathname = base_data.fullName.split('/')
        alt_pathname = alt_data.fullName.split('/')
        new_pathname[-2] = 'MergedFrom_{0}'.format(alt_pathname[1])
        new_pathname = '/'.join(new_pathname)
        print('writing: ',new_pathname)
        tsc = TimeSeriesContainer()
        tsc.times = base_hectimes
        tsc.fullName = new_pathname
        tsc.values = base_values
        tsc.units = base_units
        tsc.numberValues = len(base_values)

        dssFmOut = HecDss.open(dss_outfile)
        dssFmOut.write(tsc)

        dssFm.close()
        dssFmOut.close()


def check_start_and_end(values, times, startime, endtime):
    if times[0] < startime:  # if startdate is before the timewindow..
        print('start date ({0}) from DSS before timewindow ({1})..'.format(times[0], startime))
        st_offset = (startime - times[0]) / (times[1] - times[0])
        values = values[st_offset:]
        times = times[st_offset:]
    if times[-1] > endtime:
        print('end date ({0}) from DSS after timewindow ({1})..'.format(times[-1], endtime))
        st_offset = (times[-1] - endtime) / (times[1] - times[0])
        values = values[:(len(times) - st_offset)]
        times = times[:(len(times) - st_offset)]
    return values, times

def airtemp_lapse(dss_file,dss_rec,lapse_in_C,dss_outfile,f_part):
    dss = HecDss.open(dss_file)
    tsm = dss.read(dss_rec)
    lapse = lapse_in_C
    if 'f' in tsm.getUnits().lower():
        lapse = lapse*9.0/5.0+32.0
    tsm = tsm.add(lapse)
    tsc = tsm.getData()
    dss.close()

    pathparts = dss_rec.split('/')
    pathparts[-2] = f_part
    tsc.fullName = '/'.join(pathparts)
    dss_out = HecDss.open(dss_outfile)
    dss_out.write(tsc)
    dss_out.close()


def add_flows(currentAlt, timewindow, inflow_records, dss_file, output_dss_record_name, output_dss_file):
     
    #cfs_2_acreft = balance_period * 3600. / 43559.9
    #acreft_2_cfs = 1. / cfs_2_acreft

    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    #01Jan2014 0000
    starttime_hectime = HecTime(starttime_str).value()
    endtime_hectime = HecTime(endtime_str).value()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dss_file)

    inflows = []
    times = []

    # Read inflows
    print('Reading inflows')
    for j, inflow_record in enumerate(inflow_records): #for each of the dss paths in inflow_records
        pathname = inflow_record
        currentAlt.addComputeMessage('reading' + str(pathname))
        print('\nreading' + str(pathname))
        try:
       
            print(starttime_str, endtime_str)
            print(dss_file)
            ts = dssFm.read(pathname, starttime_str, endtime_str, False)
            ts_data = ts.getData()
            values = ts_data.values
            hectimes = ts_data.times
            units = ts_data.units
            # print('num values {0}'.format(len(values)))
            # print('start {0}'.format(ts_data.getStartTime()))
            # print('end {0}'.format(ts_data.getEndTime()))
            if hectimes[0] < starttime_hectime: #if startdate is before the timewindow..
                print('start date ({0}) from DSS before timewindow ({1})..'.format(hectimes[0], starttime_hectime))
                st_offset = (starttime_hectime - hectimes[0]) / (hectimes[1] - hectimes[0])
                values = values[st_offset:]
                hectimes = hectimes[st_offset:]
            if hectimes[-1] > endtime_hectime:
                print('end date ({0}) from DSS after timewindow ({1})..'.format(hectimes[-1], endtime_hectime))
                st_offset = (hectimes[-1] - endtime_hectime) / (hectimes[1] - hectimes[0])
                values = values[:(len(hectimes) - st_offset)]
                hectimes = hectimes[:(len(hectimes) - st_offset)]
                

        except HecMathException:
            currentAlt.addComputeMessage('ERROR reading' + str(pathname))
            sys.exit(-1)

        if units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            convvals = []
            for flow in values:
                convvals.append(flow * 35.314666213)
            values = convvals

        if len(inflows) == 0:
            inflows = values
            times = hectimes #TODO: check how this handles missing values
        else:
            for vi, v in enumerate(values):
                inflows[vi] += v

    # Output record
    tsc = TimeSeriesContainer()
    tsc.times = times
    tsc.fullName = output_dss_record_name
    tsc.values = inflows
    #tsc.startTime = times[1]
    tsc.units = 'CFS'
    #tsc.endTime = times[-1]
    tsc.numberValues = len(inflows)
    #tsc.startHecTime = timewindow.getStartTime()
    #tsc.endHecTime = timewindow.getEndTime()
    dssFm_out = HecDss.open(output_dss_file)
    dssFm_out.write(tsc)

    dssFm.close()
    dssFm_out.close()


def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'DMS_SacTrn_ResSim_Pre-Process.dss')

    hydro_dss = os.path.join(shared_dir, 'DMS_SacTrnHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_SacTrnMet.dss')
    fix_DMS_types_units(met_dss_file)

    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+met_dss_file)
    currentAlternative.addComputeMessage('lapse outfile: '+output_dss_file)
    airtemp_lapse(met_dss_file, "/MR SAC.-CLEAR CR. TO SAC R./KRDD-AIR TEMPERATURE/TEMP-AIR//1HOUR/238-235.40.53.1.1/",
                  0.7, output_dss_file, "Shasta_Lapse")

    # add PG flows 1-5 to create PG_SUM record used by ResSim/TCD scripts
    # Hmm someone has done that for us, hooray!
    
    # TODO: maybe we need to generate combined river outlet flows.  

    # TODO: Before, the Debris Dam flow needed to be generated by subtracting Spring Creek tunnel flows from the WHI Gen
    # time series, or something like that (also, we capped the tunnel flow at the max allowed by tunnel, and made the
    # subtracted debris dam flow non-negative).  We may need to do soething like that, when those records are available in
    # in the template.

    # Trinity: add Generation, G1 & G2, which are powerplant and jet-valve (bypass) flows from the powerplant intake (G3 is the low-level bypass)
    inflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/225-231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/225-231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/225-231.5.125.8.1/']
    add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Trinity Lake/TRN-GenerationG1G2_Sum/Flow//1Hour/ResSim_PreProcess/', output_dss_file)


    # Lewiston: add Generation and outlet flows which come from the same level
    inflow_records = ['/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/226-232.12.125.2.1/',
                      '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/226-232.12.125.3.1/']
    add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Lewiston Res./LEW-Gen_plus_Outlet Release/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Lewiston is still dependent on using Met data from Redding during Jan-Feb-Mar.  Create those spliced Met data records...
    pairs = [
            ["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Air temperature/Temp-Air//1Hour/242-232.6.53.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Air temperature/Temp-Air//1Hour/238-235.40.53.1.1/"],
            ["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Dew Point/Temp-DewPoint//1Hour/242-232.6.51.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Dew Point/Temp-DewPoint//1Hour/238-235.40.51.1.1/"],
            ["/MR SAC.-LEWISTON RES./TCAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/242-232.5.135.1.1/",
             "/MR SAC.-CLEAR CR. TO SAC R./RRAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/238-235.41.135.1.1/"],
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Direction/Dir-Wind/0/1Hour/242-232.5.133.1.2/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Direction/Dir-Wind//1Hour/238-235.40.133.1.2/"],
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Speed/Speed-Wind//1Hour/242-232.5.133.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Speed/Speed-Wind//1Hour/238-235.40.133.1.1/"],
            ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover//1Day/242-236.9.129.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover//1Hour/238-235.41.129.1.1/"]
             ]

    months = [1,2,3] #these are the months we are replacing with pair #2
    replace_data(currentAlternative, rtw, pairs, met_dss_file, output_dss_file, months, standard_interval='1HOUR')


    # Generate Flow records needed for plotting - mostly, sums of dam outlet flows, if proper sums do not exist
    # Trinity: calc'd total outflow is included in DMS template
    # Keswick: outlet gase only, no summing needed
    # Whiskeytown: outlet + spill = total dam outflow
    inflow_records = ['/MR Sac.-Whiskeytown Lake/WHI-Outlet Release/Flow//1Hour/227-233.14.125.2.1/',
                      '/MR Sac.-Whiskeytown Lake/WHI-Generation Release/Flow//1Hour/229-233.14.125.1.1/']
    add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Whiskeytown Lake/WHI-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Shasta: gen + spill + river outlets = total dam release
    inflow_records = ['/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/229-230.11.125.3.1/',
                      '/MR Sac.-Shasta Lake/SHA-Outlet Release/Flow//1Hour/229-230.11.125.5.1/',
                      '/MR Sac.-Shasta Lake/SHA-Spill Release/Flow//1Hour/229-230.11.125.4.1/']
    add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Shasta Lake/SHA-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Lewiston: outlet + gen + hatchery + spill = total dam outflow
    inflow_records = ['/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/226-232.12.125.2.1/',
                      '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/226-232.12.125.3.1/',
                      '/MR Sac.-Lewiston Res./LEW-Fish Hatchery Release/Flow//1Hour/226-232.12.125.1.1/',
                      '/MR Sac.-Lewiston Res./LEW-Spill Release Hrly/Flow//1Hour/226-232.12.125.5.1/']
    add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Lewiston Res./LEW-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

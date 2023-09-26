
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

import DSS_Tools
reload(DSS_Tools)

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

from Acc_Dep_ResSim_SacTrn import computeAlternative as acc_dep_5Res_ResSim

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
        base_type = base_data.type
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
        tsc.type = base_type
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


def forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'forecast_SacTrn_ResSim_Pre-Process.dss')
    forecast_dss = os.path.join(shared_dir,'WTMP_SacTrn_Forecast.dss')

    hydro_dss = os.path.join(shared_dir, 'DMS_SacTrnHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_SacTrnMet.dss')
    fix_DMS_types_units(met_dss_file)

    
    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+met_dss_file)
    currentAlternative.addComputeMessage('lapse outfile: '+output_dss_file)
    DSS_Tools.airtemp_lapse(forecast_dss, "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-AIR//1HOUR/SACTRN_BC_SCRIPT/",
                  0.7, output_dss_file, "Shasta_Lapse")
    
    # Lewiston is still dependent on using Met data from Redding during Jan-Feb-Mar.  Create those spliced Met data records...
    pairs = [
            ["/MR Sac.-Lewiston Res./TCAC1/Temp-Air//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Lewiston Res./TCAC1/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR SAC.-LEWISTON RES./TCAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/",
             "/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Lewiston Res./TCAC1/Dir-Wind/0/1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Dir-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Lewiston Res./TCAC1/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Trinity River/TCAC1/%-Cloud Cover//1Day/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover//1Hour/SACTRN_BC_SCRIPT/"]
             ]
    months = [1,2,3] #these are the months we are replacing with pair #2
    replace_data(currentAlternative, rtw, pairs, forecast_dss, output_dss_file, months, standard_interval='1HOUR')

	# Do we need to convert the evap (used in place of balance flows?) perhaps we at least need to make them negative?

    # TODO: Generate Flow records needed for plotting
    return True

def fix_DMS_parts(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    hydro_dss = os.path.join(shared_dir, 'DMS_SacTrnHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_SacTrnMet.dss')
    fix_DMS_types_units(met_dss_file)

    DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)


def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    # --------- TEMPORARY TO FIX DMS F-PARTS
    fix_DMS_parts(currentAlternative, computeOptions)
    # --------- TEMPORARY TO FIX DMS F-PARTS

    data_preprocess = forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions)

    #acc_dep = acc_dep_5Res_ResSim(currentAlternative, computeOptions)

    if data_preprocess: #and acc_dep:
    	return True

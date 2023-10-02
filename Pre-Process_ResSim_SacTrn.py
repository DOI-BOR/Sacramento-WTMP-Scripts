
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
import equilibrium_temp
reload(equilibrium_temp)
import DMS_preprocess
reload(DMS_preprocess)

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

from Acc_Dep_ResSim_SacTrn import computeAlternative as acc_dep_5Res_ResSim


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
            base = DSS_Tools.standardize_interval(base,standard_interval)
        base_data = base.getData()
        base_values = base_data.values
        base_hectimes = base_data.times
        base_units = base_data.units
        base_interval = base_data.interval
        base_type = base_data.type
        base_values, base_hectimes = check_start_and_end(base_values, base_hectimes, starttime_hectime, endtime_hectime)

        alt = dssFm.read(pair[1], starttime_str, endtime_str, False)
        if standard_interval is not None:
            alt = DSS_Tools.standardize_interval(alt,standard_interval)
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
       
def eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out):
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()

	# get at data and times in the formats needed
    dssFm = HecDss.open(at[0])        
    tsc = dssFm.read(at[1], starttime_str, endtime_str, False).getData()
    tsc_int_times = tsc.times
    dtt = DSS_Tools.hectime_to_datetime(tsc)
    at_data = tsc.values
    dssFm.close()
	# get the rest of the data over the same period
    cl_data = DSS_Tools.data_from_dss(cl[0],cl[1],starttime_str,endtime_str)
    ws_data = DSS_Tools.data_from_dss(ws[0],ws[1],starttime_str,endtime_str)
    sr_data = DSS_Tools.data_from_dss(sr[0],sr[1],starttime_str,endtime_str)
    td_data = DSS_Tools.data_from_dss(td[0],td[1],starttime_str,endtime_str)
    
	# calc_equilibrium_temp(dtt, at, cl, sr, td, ws)
    Te = equilibrium_temp.calc_equilibrium_temp(dtt,at_data,cl_data,sr_data,td_data,ws_data)
    
    print('writing: ',eq_temp_out[1])
    tsc = TimeSeriesContainer()
    tsc.times = tsc_int_times
    tsc.fullName = eq_temp_out[1]
    tsc.values = Te
    tsc.units = 'C'
    tsc.type = 'INST-VAL'
    tsc.numberValues = len(tsc.values)

    tsm = tsmath(tsc)
    tsm_day = DSS_Tools.standardize_interval(tsm,'1day')
    tsm_wk = DSS_Tools.standardize_interval(tsm,'1week')
        
    dssFmOut = HecDss.open(eq_temp_out[0])
    dssFmOut.write(tsc)
    dssFmOut.write(tsm_day)
    dssFmOut.write(tsm_wk)
    dssFmOut.close()


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
    DMS_preprocess.fix_DMS_types_units(forecast_dss)
    
    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+forecast_dss)
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
            ["/MR Sac.-Lewiston Res./TCAC1/Dir-Wind-radians/0/1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Dir-Wind-radians//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Lewiston Res./TCAC1/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Trinity River/TCAC1/%-Cloud Cover//1Day/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Trinity River/TCAC1/%-Cloud Cover-FRAC//1Day/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"]
             ]
    months = [1,2,3] #these are the months we are replacing with pair #2
    replace_data(currentAlternative, rtw, pairs, forecast_dss, output_dss_file, months, standard_interval='1HOUR')

    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")
    # eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out)
    eq_temp(rtw,
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
			[forecast_dss,"/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"],
			[forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
			[forecast_dss,"/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
			[forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            [output_dss_file,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Hour/sactrn_bc_script/"]
		   )

	# TODO: Perhaps generate tributary flows/temps based on exceedence and/or temp regressions?

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
    DMS_preprocess.fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_SacTrnMet.dss')
    DMS_preprocess.fix_DMS_types_units(met_dss_file)

    DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)


def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')

    # --------- TEMPORARY TO FIX DMS F-PARTS
    #fix_DMS_parts(currentAlternative, computeOptions)
    # --------- TEMPORARY TO FIX DMS F-PARTS

    data_preprocess = forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions)

    #acc_dep = acc_dep_5Res_ResSim(currentAlternative, computeOptions)

    if data_preprocess: #and acc_dep:
    	return True

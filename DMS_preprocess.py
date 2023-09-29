
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

sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

import Acc_Dep_ResSim_SacTrn
reload(Acc_Dep_ResSim_SacTrn)

import DSS_Tools
reload(DSS_Tools)

units_need_fixing = ['radians','tenths']

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
            if tsm.getUnits() == 'tenths':
                # save off a copy of cloud record in 0-1 for ResSim
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-FRAC'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'FRAC'
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 10.0                
                dss.write(tsc)				
            if tsm.getUnits() == 'radians':
                # save off a copy in deg
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-DEG'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'deg'
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0                
                dss.write(tsc)	
    dss.close()



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


def DMS_fix_units_types(hydro_dss,met_dss_file):
    fix_DMS_types_units(hydro_dss)
    fix_DMS_types_units(met_dss_file)


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
                #print('replaced {0}'.format(base_data.getHecTime(i)))

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


def splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file):

    # Lewiston is still dependent on using Met data from Redding during Jan-Feb-Mar.  Create those spliced Met data records...
    pairs = [
            ["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Air temperature/Temp-Air//1Hour/232.6.53.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Air temperature/Temp-Air//1Hour/235.40.53.1.1/"],
            ["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Dew Point/Temp-DewPoint//1Hour/232.6.51.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Dew Point/Temp-DewPoint//1Hour/235.40.51.1.1/"],
            ["/MR SAC.-LEWISTON RES./TCAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/232.5.135.1.1/",
             "/MR SAC.-CLEAR CR. TO SAC R./RRAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/235.41.135.1.1/"],
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Direction/Dir-Wind/0/1Hour/232.5.133.1.2/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Direction/Dir-Wind-Deg//1Hour/235.40.133.1.2/"],
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Speed/Speed-Wind//1Hour/232.5.133.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Speed/Speed-Wind//1Hour/235.40.133.1.1/"],
            ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover//1Day/236.9.129.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover//1Hour/235.41.129.1.1/"],
            ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover-FRAC//1Day/236.9.129.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover-FRAC//1Hour/235.41.129.1.1/"]
             ]

    months = [1,2,3] #these are the months we are replacing with pair #2
    replace_data(currentAlternative, rtw, pairs, met_dss_file, output_dss_file, months, standard_interval='1HOUR')
   

def compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file):
    # add PG flows 1-5 to create PG_SUM record used by ResSim/TCD scripts
    # Hmm someone has done that for us, hooray!
    
    # TODO: maybe we need to generate combined river outlet flows.  

    # TODO: Before, the Debris Dam flow needed to be generated by subtracting Spring Creek tunnel flows from the WHI Gen
    # time series, or something like that (also, we capped the tunnel flow at the max allowed by tunnel, and made the
    # subtracted debris dam flow non-negative).  We may need to do soething like that, when those records are available in
    # in the template.

    # Trinity: add Generation, G1 & G2, which are powerplant and jet-valve (bypass) flows from the powerplant intake (G3 is the low-level bypass)
    inflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/231.5.125.8.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Trinity Lake/TRN-GenerationG1G2_Sum/Flow//1Hour/ResSim_PreProcess/', output_dss_file)


    # Lewiston: add Generation and outlet flows which come from the same level
    inflow_records = ['/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/232.12.125.2.1/',
                      '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/232.12.125.3.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Lewiston Res./LEW-Gen_plus_Outlet Release/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

def compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file):

    # Generate Flow records needed for plotting - mostly, sums of dam outlet flows, if proper sums do not exist
    # Trinity: calc'd total outflow is included in DMS template
    # Keswick: outlet gase only, no summing needed
    # Whiskeytown: outlet + spill = total dam outflow
    inflow_records = ['/MR Sac.-Whiskeytown Lake/WHI-Outlet Release/Flow//1Hour/233.14.125.2.1/',
                      '/MR Sac.-Whiskeytown Lake/WHI-Generation Release/Flow//1Hour/233.14.125.1.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Whiskeytown Lake/WHI-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Shasta: gen + spill + river outlets = total dam release
    inflow_records = ['/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/',
                      '/MR Sac.-Shasta Lake/SHA-Outlet Release/Flow//1Hour/230.11.125.5.1/',
                      '/MR Sac.-Shasta Lake/SHA-Spill Release/Flow//1Hour/230.11.125.4.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Shasta Lake/SHA-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)

    # Lewiston: outlet + gen + hatchery + spill = total dam outflow
    inflow_records = ['/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/232.12.125.2.1/',
                      '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/232.12.125.3.1/',
                      '/MR Sac.-Lewiston Res./LEW-Fish Hatchery Release/Flow//1Hour/232.12.125.1.1/',
                      '/MR Sac.-Lewiston Res./LEW-Spill Release Hrly/Flow//1Hour/232.12.125.5.1/']
    DSS_Tools.add_flows(currentAlternative, rtw, inflow_records, hydro_dss,
              '/MR Sac.-Lewiston Res./LEW-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/', output_dss_file)


def preprocess_W2_5Res(currentAlternative, computeOptions):
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

    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file)
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)


def preprocess_ResSim_5Res(currentAlternative, computeOptions):
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
    
    # if template IDs exist still, remove them
    #DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    #DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)

    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file)
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+met_dss_file)
    currentAlternative.addComputeMessage('lapse outfile: '+output_dss_file)
    airtemp_lapse(met_dss_file, "/MR SAC.-CLEAR CR. TO SAC R./KRDD-AIR TEMPERATURE/TEMP-AIR//1HOUR/235.40.53.1.1/",
                  0.7, output_dss_file, "Shasta_Lapse")

    return True



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

units_need_fixing = ['tenths','m/s','deg','kph'] #'radians',]

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
            if tsm.getUnits() == 'deg':
                # save off a copy in redians
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                rec_parts[3] += '-RADIANS'
                tsc.fullName = '/'.join(rec_parts)
                tsc.units = 'radians'
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)                
                dss.write(tsc)
            if tsm.getUnits() == 'kph':
                # convert to m/s 
                tsc = tsm.getData()
                tsc.units = 'm/s'
                for i in range(len(tsc.values)) :
                    tsc.values[i] = tsc.values[i] / 3.6
                dss.write(tsc)

                # also, add w2link
                rec_parts = tsc.fullName.split('/')
                if not "w2link" in rec_parts[3].lower():
                    rec_parts[3] += '-W2link'
                    tsc.fullName = '/'.join(rec_parts)
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    dss.write(tsc)
 
            if tsm.getUnits() == 'm/s':
                # make a copy divied by kph conversion as a hack to get W2 linking the wind speed correctly 
                tsc = tsm.getData()
                rec_parts = tsc.fullName.split('/')
                if not "w2link" in rec_parts[3].lower():
                    rec_parts[3] += '-W2link'
                    tsc.fullName = '/'.join(rec_parts)
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
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

def calculate_relative_humidity(air_temp, dewpoint_temp):
    """
    Calculate Relative Humidity given the air temperature and dewpoint temperature - August-Roche-Magnus approximation
    
    :param air_temp: Air Temperature in degrees Celsius
    :param dewpoint_temp: Dew Point Temperature in degrees Celsius
    :return: Relative Humidity in percentage
    """
    numerator = (112.0 - 0.1 * dewpoint_temp + air_temp)
    denominator = (112.0 + 0.9 * air_temp)
    exponent = ((17.62 * dewpoint_temp) / (243.12 + dewpoint_temp)) - ((17.62 * air_temp) / (243.12 + air_temp))
    relative_humidity = 100.0 * (numerator / denominator) * math.exp(exponent)	
    return max(0.01,min(100.0,relative_humidity))

def relhum_from_at_dp(met_dss_file,at_path,dp_path):
    dss = HecDss.open(met_dss_file)
    tsc = dss.read(at_path).getData()
    dp_data = DSS_Tools.data_from_dss(met_dss_file,dp_path,None,None)
    for i in range(tsc.numberValues):
        tsc.values[i] = calculate_relative_humidity(tsc.values[i],dp_data[i])
    parts = tsc.fullName.split('/')
    parts[2] = parts[2][:5]
    parts[3] = 'RELHUM-FROM-AT-DP'
    parts[6] = parts[6]+'-DERIVED'
    new_pathname = '/'.join(parts)
    tsc.fullName = new_pathname
    tsc.units = '%'
    print('writing: ',new_pathname)
    dss.write(tsc)
    dss.close()


def DMS_fix_units_types(hydro_dss,met_dss_file):
    fix_DMS_types_units(hydro_dss)
    fix_DMS_types_units(met_dss_file)


def splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file, months=[1,2,3]):

    # Lewiston is still dependent on using Met data from Redding during Jan-Feb-Mar.  Create those spliced Met data records...
    pairs = [
            ["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Air temperature/Temp-Air//1Hour/232.6.53.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Air temperature/Temp-Air//1Hour/235.40.53.1.1/"],
            ["/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Dew Point/Temp-DewPoint//1Hour/232.6.51.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Dew Point/Temp-DewPoint//1Hour/235.40.51.1.1/"],
            ["/MR SAC.-LEWISTON RES./TCAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/232.5.135.1.1/",
             "/MR SAC.-CLEAR CR. TO SAC R./RRAC1-SOLAR RADIATION/IRRAD-SOLAR//1HOUR/235.41.135.1.1/"],
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Direction/Dir-Wind-RADIANS/0/1Hour/232.5.133.1.2/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Direction/Dir-Wind//1Hour/235.40.133.1.2/"],
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Speed/Speed-Wind-W2link//1Hour/232.5.133.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Speed/Speed-Wind-W2link//1Hour/235.40.133.1.1/"],
            ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover//1Day/236.9.129.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover//1Hour/235.41.129.1.1/"],
            ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover-FRAC//1Day/236.9.129.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover-FRAC//1Hour/235.41.129.1.1/"]
             ]

    DSS_Tools.replace_data(currentAlternative, rtw, pairs, met_dss_file, output_dss_file, months, standard_interval='1HOUR')
   

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

    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')

    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file, months=[1,2,3,12])
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)

	# link these flows for W2, to avoid zero-dam-flow situations - needs to be above 2.0 cfs for the flowweightaverage script to recognize it :/
    DSS_Tools.min_ts(output_dss_file, '/MR Sac.-Trinity Lake/TRN-GenerationG1G2_Sum/Flow//1Hour/ResSim_PreProcess/', 2.0, output_dss_file, 'min_flow')
    DSS_Tools.min_ts(output_dss_file, '/MR Sac.-Lewiston Res./LEW-Gen_plus_Outlet Release/Flow//1Hour/ResSim_PreProcess/', 2.0, output_dss_file, 'min_flow')


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

    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0.0, what='temp-water', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=0, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, output_dss_file, constant=1, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ONES',fpart='ONES')


                        
    # if template IDs exist still, remove them
    #DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    #DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)

    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file,months=[1,2,3])
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)

    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+met_dss_file)
    currentAlternative.addComputeMessage('lapse outfile: '+output_dss_file)
    airtemp_lapse(met_dss_file, "/MR SAC.-CLEAR CR. TO SAC R./KRDD-AIR TEMPERATURE/TEMP-AIR//1HOUR/235.40.53.1.1/",
                  -0.7, output_dss_file, "Shasta_Lapse")

    relhum_from_at_dp(met_dss_file,
	                  "/MR Sac.-Clear Cr. to Sac R./KRDD-Air temperature/Temp-Air//1Hour/235.40.53.1.1/",
	                  "/MR Sac.-Clear Cr. to Sac R./KRDD-Dew Point/Temp-DewPoint//1Hour/235.40.51.1.1/")
    relhum_from_at_dp(met_dss_file,
	                  "/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Air temperature/Temp-Air//1Hour/232.6.53.1.1/",
	                  "/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Dew Point/Temp-DewPoint//1Hour/232.6.51.1.1/")

    return True


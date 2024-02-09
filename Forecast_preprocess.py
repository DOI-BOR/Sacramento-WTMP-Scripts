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
import os, sys, csv

from com.rma.io import DssFileManagerImpl
from com.rma.model import Project
#import hec.hecmath.TimeSeriesMath as tsmath
sys.path.append(os.path.join(Project.getCurrentProject().getWorkspacePath(), "scripts"))

import DSS_Tools
reload(DSS_Tools)

import DMS_preprocess
reload(DMS_preprocess)

import equilibrium_temp
reload(equilibrium_temp)

import create_balance_flow_jython as cbfj
reload(cbfj)

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

def interp(x, xp, fp, left=None, right=None):
    """
    One-dimensional linear interpolation.

    Returns the one-dimensional piecewise linear interpolant to a function
    with given values at discrete data-points.

    Parameters
    ----------
    x : array_like
        The x-coordinates at which to evaluate the interpolated values.
    xp : 1-D sequence of floats
        The x-coordinates of the data points, must be increasing.
    fp : 1-D sequence of floats
        The y-coordinates of the data points, same length as `xp`.
    left : float, optional
        Value to return for `x < xp[0]`, default is `fp[0]`.
    right : float, optional
        Value to return for `x > xp[-1]`, default is `fp[-1]`.

    Returns
    -------
    y : float or ndarray
        The interpolated values, same shape as `x`.
    """

    if isinstance(x, list):
        return [interp(point, xp, fp, left, right) for point in x]
    else:
        if left is None:
            left = fp[0]
        if right is None:
            right = fp[-1]

        if x < xp[0]:
            return left
        elif x > xp[-1]:
            return right
        else:
            for i in range(len(xp) - 1):
                if x >= xp[i] and x <= xp[i+1]:
                    # Perform the linear interpolation
                    t = (x - xp[i]) / (xp[i+1] - xp[i])
                    return fp[i] + t * (fp[i+1] - fp[i])



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


def storage_to_elev(res_name,elev_stor_area,forecast_dss,storage_rec,conic=False):
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True) # read ALL data in record

    elev = []
    if conic:
        print('Conic interpolation of elevations from storage not supported yet.')
        sys.exit(-1)
    else:
        for j in range(tsc.numberValues):
            elev.append(cbfj.linear_interpolation(elev_stor_area['stor'], elev_stor_area['elev'], tsc.values[j]))
            print('stor2elev: ',j,tsc.times[j],tsc.values[j])

    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = elev
    dssFmRec.write(tsc)
    dssFmRec.close()

def invent_elevation(res_name,forecast_dss,storage_rec,elev_constant_ft):
    dssFmRec = HecDss.open(forecast_dss)
    tsc = dssFmRec.get(storage_rec,True)
    recparts = tsc.fullName.split('/')
    recparts[2] = res_name
    recparts[3] = 'ELEV'
    tsc.fullName = '/'.join(recparts)
    tsc.units = 'ft'
    tsc.type = 'INST-VAL'
    tsc.values = [elev_constant_ft for j in range(tsc.numberValues)]
    dssFmRec.write(tsc)
    dssFmRec.close()

def write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir):

    # get date for starting elevation - look for end-of-month before start time
    ht = HecTime(rtw.getStartTimeString())
    if ht.month() == 1:
        start_str = dt.datetime(ht.year()-1,12,31).strftime('%d%b%Y')+ ' 2400'
    else:
        start_dt = dt.datetime(ht.year(),ht.month(),1)
        start_dt = start_dt - dt.timedelta(days=1)
        start_str = start_dt.strftime('%d%b%Y')+ ' 2400'
    end_str = rtw.getEndTimeString()   
    

    # covert storage to monthly elevation
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_shasta.csv'), 'Shasta')
    storage_to_elev('Shasta Lake',elev_stor_area,forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/',conic=False)

    # invent flow-reg reservoir elevation record from shasta storage rec (used for timing only)
    invent_elevation('Keswick Reservoir',forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/',582.0)
    invent_elevation('Lewiston Reservoir',forecast_dss,'/TRINITY RIVER/TRINITY LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/',1901.0)
    
    # write an hourly forecast elevation based on starting elevation and flows
    DSS_Tools.resample_dss_ts(forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/',None,forecast_dss,'1DAY')
    inflow_records = ['//SHASTA-PIT-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SHASTA-SAC-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SHASTA-SULANHARAS-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SHASTA-MCCLOUD-IN/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '/SACRAMENTO RIVER/SHASTA LAKE/FLOW-ACC-DEP//1DAY/SACTRN_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['/SACRAMENTO RIVER/SHASTA LAKE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)
    cbfj.predict_elevation(currentAlternative, start_str,end_str, 'Shasta Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Shasta Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')

    # Trinity
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_trinity.csv'), 'trinity')
    storage_to_elev('Trinity Lake',elev_stor_area,forecast_dss,'/TRINITY RIVER/TRINITY LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/',conic=False)

    DSS_Tools.resample_dss_ts(forecast_dss,'/TRINITY RIVER/TRINITY LAKE/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/',None,forecast_dss,'1DAY')
    inflow_records = ['//EF TRINITY/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//STUART FORK/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//SWIFT CR/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '//TRINITY RIVER/FLOW-IN//1DAY/SACTRN_BC_SCRIPT/',
                      '/TRINITY RIVER/TRINITY LAKE/FLOW-ACC-DEP//1DAY/SACTRN_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['/TRINITY RIVER/TRINITY LAKE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'/TRINITY RIVER/TRINITY LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)
    cbfj.predict_elevation(currentAlternative, start_str,end_str, 'Trinity Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Trinity Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')

    # Whiskeytown
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_whiskeytown.csv'), 'Whiskeytown')
    storage_to_elev('Whiskeytown Lake',elev_stor_area,forecast_dss,'/CLEAR CREEK/WHISKEYTOWN LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/',conic=False)

    DSS_Tools.resample_dss_ts(forecast_dss,'/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1HOUR/SACTRN_BC_SCRIPT/',None,forecast_dss,'1DAY')
    DSS_Tools.resample_dss_ts(forecast_dss,'/CLEAR CREEK/WHISKEYTOWN DAM/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/',None,forecast_dss,'1DAY')
    DSS_Tools.resample_dss_ts(forecast_dss,'/CLEAR CREEK/CARR POWERHOUSE/FLOW-RELEASE//1HOUR/SACTRN_BC_SCRIPT/',None,forecast_dss,'1DAY')
    inflow_records = ['/USBR-LINEARINTERP/CLEAR CR ABOVE JCR INFLOW/FLOW//1DAY/SACTRN_BC_SCRIPT/',
                      '/CLEAR CREEK/CARR POWERHOUSE/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/',
                      '/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-ACC-DEP//1DAY/SACTRN_BC_SCRIPT/']  # this actually evap, but negative already, so it goes as inflow
    outflow_records = ['/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1DAY/SACTRN_BC_SCRIPT/',
                       '/CLEAR CREEK/WHISKEYTOWN DAM/FLOW-RELEASE//1DAY/SACTRN_BC_SCRIPT/']
    starting_elevation = DSS_Tools.first_value(forecast_dss,'/CLEAR CREEK/WHISKEYTOWN LAKE/ELEV//1MON/SACTRN_BC_SCRIPT/',start_str,end_str)
    print('starting_elevation ',starting_elevation)
    cbfj.predict_elevation(currentAlternative, start_str,end_str, 'Whiskeytown Lake', inflow_records, outflow_records, starting_elevation,
                         elev_stor_area, forecast_dss, '//Whiskeytown Lake/ELEV-FORECAST//1DAY/AMER_BC_SCRIPT/', forecast_dss, shared_dir,
                         use_conic=False, alt_period=None, alt_period_string=None, balance_period_str='1Day')

def splice_met(currentAlternative, rtw, forecast_dss, output_dss):
    # Lewiston is still dependent on using Met data from Redding during Jan-Feb-Mar.  Create those spliced Met data records...
    pairs = [
            ["/MR Sac.-Lewiston Res./TCAC1/Temp-Air//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Lewiston Res./TCAC1/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR SAC.-LEWISTON RES./TCAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/",
             "/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            #["/MR Sac.-Lewiston Res./TCAC1/Dir-Wind-radians//1Hour/SACTRN_BC_SCRIPT/", # removing this, b/c this is not used and the unit diff is causing forecast problems
            #"/MR Sac.-Clear Cr. to Sac R./KRDD/Dir-Wind//1Hour/SACTRN_BC_SCRIPT/"],  # already in radians
            ["/MR Sac.-Lewiston Res./TCAC1/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Trinity River/TCAC1/%-Cloud Cover//1Day/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover//1Hour/SACTRN_BC_SCRIPT/"],
            ["/MR Sac.-Trinity River/TCAC1/%-Cloud Cover-FRAC//1Day/SACTRN_BC_SCRIPT/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"]
             ]
    months = [1,2,3] #these are the months we are replacing with pair #2
    DSS_Tools.replace_data(currentAlternative, rtw, pairs, forecast_dss, output_dss, months, standard_interval='1HOUR')


def study_dir_from_run_dir(run_dir):
    w2sim,_ = os.path.split(run_dir)
    runs_dir,_ = os.path.split(w2sim)
    study_dir,_ = os.path.split(runs_dir)
    return study_dir

def model_dir_from_run_dir(run_dir,model_place,model_name):
    study_dir = study_dir_from_run_dir(run_dir)
    model_dir = os.path.join(study_dir,'cequal-w2',model_place,model_name)
    return model_dir

def forecast_data_preprocess_ResSim_5Res(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    forecast_dss = os.path.join(shared_dir,'WTMP_SacTrn_Forecast.dss')
    DMS_preprocess.fix_DMS_types_units(forecast_dss)
    
    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+forecast_dss)
    currentAlternative.addComputeMessage('lapse outfile: '+forecast_dss)
    DSS_Tools.airtemp_lapse(forecast_dss, "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-AIR//1HOUR/SACTRN_BC_SCRIPT/",
                  0.7, forecast_dss, "Shasta_Lapse")
    
    splice_met(currentAlternative, rtw, forecast_dss, forecast_dss)

    write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    currentAlternative.addComputeMessage("Computing equilibrium temperature, this may take a while...")
    # eq_temp(rtw,at,cl,ws,sr,td,eq_temp_out)
    eq_temp(rtw,
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Air//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./RRAC1/%-Cloud Cover-FRAC//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Speed-Wind//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR SAC.-CLEAR CR. TO SAC R./RRAC1/IRRAD-SOLAR//1HOUR/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-DewPoint//1Hour/SACTRN_BC_SCRIPT/"],
            [forecast_dss,"/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Hour/sactrn_bc_script/"]
           )

    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')   
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=1, what='gate', 
                        dss_type='INST-VAL', period='1HOUR',cpart='ONES',fpart='ONES')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=13.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='WHI-target-13',fpart='WHI-target-13')

    DSS_Tools.relhum_from_at_dp(forecast_dss,
                      "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-AIR//1HOUR/sactrn_bc_script/",
                      "/MR SAC.-CLEAR CR. TO SAC R./KRDD/TEMP-DEWPOINT//1HOUR/sactrn_bc_script/")

    # TODO: Perhaps generate tributary flows/temps based on exceedence and/or temp regressions?

    # TODO: Generate Flow records needed for plotting
    return True

def forecast_data_preprocess_W2_5Res(currentAlternative, computeOptions):
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    # output from scripting all goes to the <study>_Forecast.dss file.
    forecast_dss = os.path.join(shared_dir,'WTMP_SacTrn_Forecast.dss')
    DMS_preprocess.fix_DMS_types_units(forecast_dss)
    
    splice_met(currentAlternative, rtw, forecast_dss, forecast_dss)

    write_forecast_elevations(currentAlternative, rtw, forecast_dss, shared_dir)

    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.001, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='TinyFlow',fpart='TinyFlow')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=10.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='TENS',fpart='TENS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1DAY',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=0.0, what='flow', 
                        dss_type='PER-AVER', period='1HOUR',cpart='ZEROS',fpart='ZEROS')
    DSS_Tools.create_constant_dss_rec(currentAlternative, rtw, forecast_dss, constant=13.0, what='temp-water', 
                        dss_type='PER-AVER', period='1DAY',cpart='WHI-target-13',fpart='WHI-target-13')

    # Keswick need daily record.
    DSS_Tools.resample_dss_ts(forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/FLOW-KESWICK-CFS//1MON/SACTRN_BC_SCRIPT/',rtw,forecast_dss,'1DAY')

    # TODO: Generate Flow records needed for plotting
    return True

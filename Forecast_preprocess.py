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
import os, sys, csv, calendar

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
    # convert bc script predicted storage
    storage_to_elev('Shasta Lake',elev_stor_area,forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/',conic=False)

    # invent flow-reg reservoir elevation record from shasta storage rec (used for timing only)
    invent_elevation('Keswick Reservoir',forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/',582.0)
    invent_elevation('Lewiston Reservoir',forecast_dss,'/TRINITY RIVER/TRINITY LAKE/STORAGE//1MON/SACTRN_BC_SCRIPT/',1901.0)

    # also make a one day step, to see if that solves some issues (used for timing only)
    invent_elevation('Keswick Reservoir',forecast_dss,'/SACRAMENTO RIVER/SHASTA LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/',582.0)
    invent_elevation('Lewiston Reservoir',forecast_dss,'/TRINITY RIVER/TRINITY LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/',1901.0)

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
    # convert bc script predicted storage
    storage_to_elev('Trinity Lake',elev_stor_area,forecast_dss,'/TRINITY RIVER/TRINITY LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/',conic=False)

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
    # convert bc script predicted storage
    storage_to_elev('Whiskeytown Lake',elev_stor_area,forecast_dss,'/CLEAR CREEK/WHISKEYTOWN LAKE/STORAGE-CVP//1Day/SACTRN_BC_SCRIPT/',conic=False)

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

	# compute W2 regression downstream target temps, in case we want to try/use them in ResSim
    # Keswick need daily record.
    DSS_Tools.resample_dss_ts(forecast_dss,'//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Hour/SACTRN_BC_SCRIPT/',rtw,forecast_dss,'1DAY',pad_start_days=1)
    DSS_Tools.resample_dss_ts(forecast_dss,'/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/',rtw,forecast_dss,'1DAY')

    location = 2
    TT_rec = "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/"
    TT_W2_rec = "/USBR/SHASTA/TEMP-WATER-TARGET-W2-UPSTREAM//1Day/SACTRN_BC_SCRIPT/"
    if location == 0: 
        # @ Shasta Dam, use exact TT
        DSS_tools.copy_dss_ts(TT_rec,new_dss_rec=TT_W2_rec,dss_file_path=forecast_dss)
    else:
        upstream_target(forecast_dss,rtw,
                        "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/",
                        "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Day/sactrn_bc_script/",
                        "//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Day/SACTRN_BC_SCRIPT/",
                        "/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1Day/SACTRN_BC_SCRIPT/",
                        location,TT_W2_rec,ResSimRiver=True)
    
    return True

def route_downstream(tKeswick,keswickFlowDaily,eqTempDaily,step,hour_of_day,loc):
    exchCoef = 0.015 # hourly exchange rate between atmosphere and river temp

    hrs = travel_time_hrs(loc,keswickFlowDaily[step]) # loc 2 = CCR
    t = tKeswick
    i = step # + keswickResidenceTime ?
    imax = len(eqTempDaily)
    h = hour_of_day
    for k in range(hrs):
        deltaTemp = (eqTempDaily[i] - t) * exchCoef
        t += deltaTemp
        if h == 23:
            h = 0
            i = min(imax-1, i+1)
        else:
            h += 1
    return t


def travel_time_hrs(loc,keswickFlow):

    if loc == 1:  # Highway 44
        downstreamDistance = 30000.  # in feet
    elif loc == 2:  # CCR
        downstreamDistance = 53000.
    elif loc == 3:  # Ball's Ferry
        downstreamDistance = 137000.
    else:
        raise NotImplementedError('Downstream location index ' + str(loc) + ' not recognized.')

    # Power law approximation for velocity in the Sacramento River
    Kcoef = 2.
    alpha = 0.33
    velocity = Kcoef * (keswickFlow / 1000) ** alpha  # power law approximation
    # Calculate travel time in seconds
    if velocity <= 0.0:
        print('WARNING: Keswick flow <= 0 in Forecast TCD script. Specified flows may be incorrect.')
        travTime = 0.0
    else:
        travTime = downstreamDistance / velocity

    # return hrs
    return int(travTime/(60.0*60.0))

def fractional_month(date_obj):
    """
    Args:
        dates: A  datetime.datetime object.

    Returns:
        fractional month & multipliers.
    """
    year = date_obj.year
    month = date_obj.month
    day = date_obj.day

    # Get the number of days in the current month using the calendar module
    _, days_in_month = calendar.monthrange(year, month)

    # Calculate the fractional month. Use 1.0 to force float division.
    fractional_month = (day - 1.0) / days_in_month
    fractional_previous = 0.0
    fractional_next = 0.0
    if fractional_month > 0.5:
        fractional_current = (1.0 - fractional_month)+0.5
        fractional_next = 1.0 - fractional_current
    else:
        fractional_current = fractional_month + 0.5
        fractional_previous = 1.0 - fractional_current

    return fractional_month,fractional_previous,fractional_current,fractional_next



def get_step_future_and_RiverHrs(wqTargetDaily,keswickFlowDaily,step,loc):

    # Power law approximation for velocity in the Sacramento River
    hrs = travel_time_hrs(loc,keswickFlowDaily[step])
    
    # Get Keswick pool information - from forecast TCD script
    flowVol = keswickFlowDaily[step] * 86400.
    kesConPoolVol = 20100. * 43560.  # cubic feet, assumed this is top of conservation
    kesFraction = flowVol / kesConPoolVol
    multiplier = 0.14  # Calibration factor
    kesFraction = min(kesFraction * multiplier, 1.)

    travel_days = int(round(hrs/24 + kesFraction)) # TODO: check is this is a good representaion of keswick travel/residence
    step_future = min(step+travel_days,len(wqTargetDaily))
    
    return step_future,hrs

#######################################################################################################
# Backcalculate the temperature required at Shasta Dam from the downstream temperature target
def backRouteWQTarget2(eqTempDaily, targetTempFuture, sha2kes_diff, hrs, step, loc, hour_of_day=10):
    
    #exchCoef = 0.0098  # exchange rate between atmosphere and river temp
    exchCoef = 0.015 # exchange rate between atmosphere and river temp
    tSearchMin = targetTempFuture + sha2kes_diff - 6.
    tSearchMax = tSearchMin + 18.
    numIters = 51
    bracketed = False
    cantBeMet = False
    #network.printMessage('Keswick vars ' + str(keswickResAvgTemp) + ', ' + str(kesFraction))
    #network.printMessage('Travel time steps ' + str(travTimeSteps))
    for j in range(numIters):
        outletTemp = tSearchMin + float(j) / float(numIters+1) * (tSearchMax - tSearchMin)
        # Impact of Keswick
        t = outletTemp - sha2kes_diff # sha2kes_diff is negative if heating in Keswick
        t1 = 1
        # Route downstream
        i = step 
        imax = len(eqTempDaily)
        h = hour_of_day
        for k in range(hrs):
            deltaTemp = (eqTempDaily[i] - t) * exchCoef
            t += deltaTemp
            if h == 23:
                h = 0
                i = min(imax-1, i+1)
            else:
                h += 1

        print('--',j,t1,t)

        if j == 0:
            prevT = t
            prevOutletT = outletTemp
        if t > targetTempFuture and prevT < targetTempFuture:
            upperOutletT = outletTemp
            upperT = t
            lowerOutletT = prevOutletT
            lowerT = prevT
            bracketed = True
            print('Break loop 1 ' + str(prevT) + ', ' + str(t) + ', ' + str(targetTempFuture))
            break
        elif prevT > targetTempFuture and t < targetTempFuture:
            lowerOutletT = outletTemp
            lowerT = t
            upperOutletT = prevOutletT
            upperT = prevT
            bracketed = True
            print('Break loop 2 ' + str(prevT) + ', ' + str(t) + ', ' + str(targetTempFuture))
            break
        elif j == 0 and t > targetTempFuture:
            cantBeMet = True
            break
        prevT = t
        prevOutletT = outletTemp

    if bracketed:
        # Linear interpolation
        targetTemp = (targetTempFuture - lowerT) / (upperT - lowerT) * (upperOutletT - lowerOutletT) + lowerOutletT
        #network.printMessage('Upper, lower T ' + str(upperT) + ', ' + str(lowerT))
        #network.printMessage('Upper, lower outlet T ' + str(upperOutletT) + ', ' + str(lowerOutletT))
        #network.printMessage('Backrouted target temp ' + str(targetTemp))
    elif cantBeMet:
        targetTemp = outletTemp
    else:
        if t < targetTempFuture:
            targetTemp = outletTemp
        else:
            #network.printMessage('Target Temperature Downstream' + str(targetTempFuture))
            raise ValueError('Outlet temperature not bracketed')
    
    return targetTemp

def upstream_target(forecastDSS,rtw,downstreamTT_rec,eqTemp_rec,kesFlow_rec,sppFlow_rec,loc,TT_W2_rec,ResSimRiver=True):

    kes2sha_coeffs = [
        [ 0.5026658277531562 , -0.023436069646011252 , -0.06557997671807375 , -4.685462463211685e-06 ],
        [ -1.0945789953889893 , -0.03547324181218445 , 0.3634668604715895 , -9.493485475182071e-05 ],
        [ -2.2979699638623208 , -0.024911009946374824 , 0.6188147891878977 , -0.00011857144326738345 ],
        [ -3.166418078853781 , -0.031958344183594084 , 0.8406859050306963 , -0.00012072892037607769 ],
        [ -4.013361936990245 , -0.038181594634872126 , 1.0869790557430794 , -0.00015667409032667326 ],
        [ -5.926346258972109 , -0.031964680347665884 , 1.5240375730724862 , -0.00016933328522669477 ],
        [ -7.90033309005282 , -0.020825728153169958 , 1.894831160970044 , -0.00010526577327102076 ],
        [ -7.152475999113026 , -0.015734400267369567 , 1.7102290556150588 , -0.0002686871465276408 ],
        [ 2.13934407431362 , -0.0320109908976759 , -0.5277579029594971 , -0.00025007383389563745 ],
        [ 1.8607414748455695 , -0.05343220981005764 , -0.33435348319941816 , -4.181194104351927e-05 ],
        [ 3.0620743202772447 , -0.031718714577417005 , -0.720332573077426 , 4.455859225551441e-05 ],
        [ 0.8715248291341835 , -0.014522140866670415 , -0.14027205100034504 , -8.962822645247351e-06 ],
    ]

    ccr2sha_coeffs = [
        [ 1.792959080538937 , -0.059354492077861414 , -0.2960320376671016 , 5.046191161081086e-07 ],
        [ -0.7512647356887353 , -0.08260980765521911 , 0.3989728525297206 , -0.00010475304117792841 ],
        [ -2.2942016504506095 , -0.07724045046488413 , 0.7261448433127441 , -0.000105287304891731 ],
        [ -4.4991291842472005 , -0.07450124827763506 , 1.2838524899115868 , -0.00015759183944165977 ],
        [ -5.891695527938446 , -0.06622045217750946 , 1.5982456393750977 , -0.00015336390252308166 ],
        [ -9.012911351050086 , -0.05819889006557895 , 2.3237120979677863 , -0.00017729324265147367 ],
        [ -11.807214560294568 , -0.04586684354549874 , 2.88558316872268 , -9.275162445499182e-05 ],
        [ -10.475158908615029 , -0.04277016383124388 , 2.5906412962866536 , -0.00027885911475739506 ],
        [ 0.7463268269374388 , -0.06582341829199213 , -0.08213960515083435 , -0.0002556556527696652 ],
        [ 3.3789888674223043 , -0.08818106994859856 , -0.6214925066947776 , -2.910878547596829e-05 ],
        [ 8.334250042494999 , -0.06392276058932216 , -2.003215573079848 , 6.311125793283819e-05 ],
        [ 4.138461583069988 , -0.036555149612575485 , -0.8940381389461702 , 3.245473186093828e-05 ],
    ]

    if loc == 2:
        if ResSimRiver:
            coeffs = kes2sha_coeffs
        else:
            coeffs = ccr2sha_coeffs
    else:
        raise ValueError("upstream_target: loc unknown, currently must be one of  {2,}")

    # load datasets over rtw, must be same lengths
    # TODO: should check units, expecting CFS and C
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()
    dssFm = HecDss.open(forecastDSS)        
    tsc = dssFm.read(downstreamTT_rec, starttime_str, endtime_str, False).getData()
    tsc_int_times = tsc.times
    dtt = DSS_Tools.hectime_to_datetime(tsc)
    downstreamTT = tsc.values
    dssFm.close()
    # get the rest of the data over the same period
    eqTemp = DSS_Tools.data_from_dss(forecastDSS,eqTemp_rec,starttime_str,endtime_str)
    kesFlow = DSS_Tools.data_from_dss(forecastDSS,kesFlow_rec,starttime_str,endtime_str)
    sppFlow = DSS_Tools.data_from_dss(forecastDSS,sppFlow_rec,starttime_str,endtime_str)

    # first daily eqTemp seem to often be bad, maybe timezone/DSS reading issue
    eqTemp[0] = eqTemp[1]

    shastaTT = []
    for i, TT in enumerate(downstreamTT):
        mo_i = dtt[i].month - 1
        mo_i_prev = mo_i-1 if mo_i-1 >= 0 else 11
        mo_i_next = mo_i+1 if mo_i+1 < 12 else 0
        _,pFrac,cFrac,nFrac = fractional_month(dtt[i])
        
        a0,b0,c0,d0 = coeffs[mo_i_prev]
        a1,b1,c1,d1 = coeffs[mo_i]
        a2,b2,c2,d2 = coeffs[mo_i]
        
        print(i,' : ',dtt[i].year,dtt[i].month,dtt[i].day,dtt[i].hour,' : ',TT,eqTemp[i],kesFlow[i],sppFlow[i])
        upstreamTT = -1
        #try:
        TTDiff = pFrac * (a0 + b0*eqTemp[i] + c0*math.log10(kesFlow[i]) + d0*sppFlow[i]) + \
                 cFrac * (a1 + b1*eqTemp[i] + c1*math.log10(kesFlow[i]) + d1*sppFlow[i]) + \
                 nFrac * (a2 + b2*eqTemp[i] + c2*math.log10(kesFlow[i]) + d2*sppFlow[i])

        if ResSimRiver:
            i_future,hrs = get_step_future_and_RiverHrs(downstreamTT,kesFlow,i,loc)
            upstreamTT = backRouteWQTarget2(eqTemp, downstreamTT[i_future], TTDiff, hrs, i, loc, hour_of_day=10)
        else:
            upstreamTT = TT+TTDiff
        print('result : ',TT,TTDiff,TT+TTDiff,upstreamTT)
        #except:
        #    print('failed!')
        shastaTT.append(upstreamTT)

    # copy over first record, which always fails maybe due to time zone read issues
    print('len(shastaTT)',len(shastaTT))
    shastaTT[0] = shastaTT[1]

    dssFmRec = HecDss.open(forecastDSS)
    tsc.fullName = TT_W2_rec
    tsc.values = shastaTT
    dssFmRec.write(tsc)
    dssFmRec.close()


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
    DSS_Tools.resample_dss_ts(forecast_dss,'//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Hour/SACTRN_BC_SCRIPT/',rtw,forecast_dss,'1DAY',pad_start_days=1)
    DSS_Tools.resample_dss_ts(forecast_dss,'/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/',rtw,forecast_dss,'1DAY')

    # read location form DSS somehow?
    location = 2 

    TT_rec = "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/"
    TT_W2_rec = "/USBR/SHASTA/TEMP-WATER-TARGET-W2-UPSTREAM//1Day/SACTRN_BC_SCRIPT/"
    if location == 0: 
        # @ Shasta Dam, use exact TT
        DSS_tools.copy_dss_ts(TT_rec,new_dss_rec=TT_W2_rec,dss_file_path=forecast_dss)
    else:
        upstream_target(forecast_dss,rtw,
                        "/USBR/SHASTA/TEMP-WATER-TARGET//1Day/SACTRN_BC_SCRIPT/",
                        "/MR Sac.-Clear Cr. to Sac R./KRDD/Temp-Equil//1Day/sactrn_bc_script/",
                        "//SHASTA/FLOW-RELEASE-KESWICK-CFS//1Day/SACTRN_BC_SCRIPT/",
                        "/CLEAR CREEK/WHISKEYTOWN LAKE/FLOW-DIVERSION-SPRING-CR//1Day/SACTRN_BC_SCRIPT/",
                        location,TT_W2_rec,ResSimRiver=False)
    
    

    # TODO: Generate Flow records needed for plotting
    return True

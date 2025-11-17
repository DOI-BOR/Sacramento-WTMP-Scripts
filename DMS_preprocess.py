
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

from com.rma.io import DssFileManagerImpl
from java.util import TimeZone

import Acc_Dep_ResSim_SacTrn
reload(Acc_Dep_ResSim_SacTrn)

import DSS_Tools
reload(DSS_Tools)

units_need_fixing = ['tenths','deg','kph'] #'radians',]

def fix_DMS_types_units(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    recs = DSS_Tools.get_sanitized_record_list(dss_file)

    dss = HecDss.open(dss_file)
    
    for r in recs:
        rlow = r.lower()
        # things not to read: paired data, integer/scalar/text vars and some
        # other things that are causing trouble.
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow and not 'temp-water-target' in rlow:
        
            tsc = dss.get(r,True)

            if "/flow" in rlow or "/1day/" in rlow:
                if not "/storage" and not "/stor" in rlow:
                    tsc.type = 'PER-AVER'
                    #tsc.setStoreAsDoubles(True)
                    dss.write(tsc)

            units = str(tsc.units).lower()  # just to make sure
            
            if units in units_need_fixing:
                if units == 'tenths':
                    # save off a copy of cloud record in 0-1 for ResSim
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 10.0
                    #tsc.setStoreAsDoubles(True)         
                    dss.put(tsc)
                if units == 'radians':
                    # save off a copy in deg
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-DEG'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'deg'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / (2*3.141592653589793) * 360.0
                    #tsc.setStoreAsDoubles(True)         
                    dss.put(tsc)
                if units == 'deg':
                    # save off a copy in redians
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-RADIANS'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'radians'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 360.0 * (2*3.141592653589793)
                    #tsc.setStoreAsDoubles(True)            
                    dss.put(tsc)
                if units == 'kph':
                    # convert to m/s 
                    tsc.units = 'm/s'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    #tsc.setStoreAsDoubles(True)
                    dss.put(tsc)

    dss.close()


def fix_DMS_types_units_old(dss_file):
    '''This method was implemented to change data types to PER-AVER that are not coming from the DMS that way'''
    recs = DSS_Tools.get_sanitized_record_list(dss_file)

    dss = HecDss.open(dss_file)
    
    for r in recs:
        rlow = r.lower()
        # things not to read: paired data, integer/scalar/text vars and some
        # other things that are causing trouble.
        if not '/location info' in rlow and not '/temp-equil' in rlow and \
          not '/depth-temp' in rlow and not 'icpathsmap' in rlow and \
          not '/downstream_control_loc' in rlow:
        
            tsm = dss.read(r)

            if "/flow" in rlow or "/1day/" in rlow:
                tsm.setType('PER-AVER')
                tsc = tsm.getData()
                #tsc.setStoreAsDoubles(True)
                dss.write(tsc)
            
            if tsm.getUnits().lower() in units_need_fixing:
                if tsm.getUnits() == 'tenths':
                    # save off a copy of cloud record in 0-1 for ResSim
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    rec_parts[3] += '-FRAC'
                    tsc.fullName = '/'.join(rec_parts)
                    tsc.units = 'FRAC'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 10.0
                    #tsc.setStoreAsDoubles(True)         
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
                    #tsc.setStoreAsDoubles(True)         
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
                    #tsc.setStoreAsDoubles(True)            
                    dss.write(tsc)
                if tsm.getUnits() == 'kph':
                    # convert to m/s 
                    tsc = tsm.getData()
                    tsc.units = 'm/s'
                    for i in range(len(tsc.values)) :
                        tsc.values[i] = tsc.values[i] / 3.6
                    #tsc.setStoreAsDoubles(True)
                    dss.write(tsc)

                    # also, add w2link
                    #rec_parts = tsc.fullName.split('/')
                    #if not "w2link" in rec_parts[3].lower():
                    #    rec_parts[3] += '-W2link'
                    #    tsc.fullName = '/'.join(rec_parts)
                    #    for i in range(len(tsc.values)) :
                    #        tsc.values[i] = tsc.values[i] / 3.6
                    #    dss.write(tsc)

                if tsm.getUnits() == 'm/s':
                    # make a copy divied by kph conversion as a hack to get W2 linking the wind speed correctly 
                    tsc = tsm.getData()
                    rec_parts = tsc.fullName.split('/')
                    if not "w2link" in rec_parts[3].lower():
                        rec_parts[3] += '-W2link'
                        tsc.fullName = '/'.join(rec_parts)
                        for i in range(len(tsc.values)) :
                            tsc.values[i] = tsc.values[i] / 3.6
                        #tsc.setStoreAsDoubles(True)
                        dss.write(tsc)
    dss.close()


def standardize_bc_temp_water_to_C(dss_file,output_dss_file):
    '''Make copies of temp-water records in C (standardizing on C) for ResSim linking'''
    recs = DSS_Tools.get_sanitized_record_list(dss_file)
    dss = HecDss.open(dss_file)

    if dss_file == output_dss_file:
        dss_out = dss
    else:
        dss_out = HecDss.open(output_dss_file)
    
    for r in recs:
        rlow = r.lower()
        if '/temp-water' in rlow:

            tsc = dss.get(r,True)

            incoming_units = tsc.units.lower()
        
            tsc = dss.get(r,True)
            rec_parts = tsc.fullName.split('/')
            rec_parts[3] += '-C'
            tsc.fullName = '/'.join(rec_parts)
            tsc.units = 'C'
                        
            if incoming_units == 'f' or incoming_units == 'degf':                
                for i in range(len(tsc.values)) :
                    tsc.values[i] = (tsc.values[i] - 32.0)*5.0/9.0             

            dss_out.put(tsc)

    dss.close()
    if dss_file != output_dss_file:
        dss_out.close()


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
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Direction/Dir-Wind/0/1Hour/232.5.133.1.2/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Direction/Dir-Wind//1Hour/235.40.133.1.2/"],
            ["/MR Sac.-Lewiston Res./TCAC1-Wind Speed/Speed-Wind//1Hour/232.5.133.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./KRDD-Wind Speed/Speed-Wind//1Hour/235.40.133.1.1/"],
            ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover//1Day/236.9.129.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover//1Hour/235.41.129.1.1/"],
            ["/MR Sac.-Trinity River/TCAC1 - Calc Data-Cloud Cover/%-Cloud Cover-FRAC//1Day/236.9.129.1.1/",
             "/MR Sac.-Clear Cr. to Sac R./RRAC1-Cloud Cover/%-Cloud Cover-FRAC//1Hour/235.41.129.1.1/"]
             ]

    DSS_Tools.replace_data(currentAlternative, rtw, pairs, met_dss_file, output_dss_file, months, standard_interval='1HOUR')
   

def compute_river_balance_flows(currentAlternative, rtw, hydro_dss, obs_dss_file, output_dss_file):

    # balance at IGO
    flow_records = [obs_dss_file + "::/USGS CLEAR CR/11372000 NR IGO/FLOW//1HOUR/USGS-MERGED-CROP/",
                    output_dss_file+'::/MR Sac.-Whiskeytown Lake/WHI-Total Dam Outflow/Flow//1Hour/ResSim_PreProcess/']
    out_rec = "/CLEAR CR/IGO BALANCE FLOW/FLOW//1HOUR/ResSim_PreProcess/"
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_records, hydro_dss,
                              [True, False],
                              out_rec, output_dss_file, what="flow",prepend_n=25)
    DSS_Tools.resample_dss_ts(output_dss_file,out_rec,rtw,output_dss_file,'1DAY')

    # balance at Bend Bridge
    DSS_Tools.resample_dss_ts(hydro_dss,'/MR Sac.-Keswick Res./KES-Dam Total Release/Flow//1Hour/234.1.125.1.1/',rtw,output_dss_file,'1DAY',pad_start_days=1)
    DSS_Tools.resample_dss_ts(hydro_dss,'/MR Sac.-Sac River/11377100 Sacramento R at Bend Bridge-15min Flow/Flow//15Minute/237.64.125.1.1/',rtw,output_dss_file,'1DAY',pad_start_days=1)

    # bend bridge is pretty far downstream from the main release at Keswick - shifting backward by one day
    # makes for a better balance
    in_rec = '/MR Sac.-Sac River/11377100 Sacramento R at Bend Bridge-15min Flow/Flow//1Day/237.64.125.1.1/'
    
    # not using out_rec, time shifting is not working
    out_rec = '/MR Sac.-Sac River/11377100 Sacramento R at Bend Bridge-15min Flow/Flow//1Day/237.64.125.1.1-1Day/'
    DSS_Tools.shift_ts_time(output_dss_file,in_rec,output_dss_file,out_rec,'-1Day',start_date=None,end_date=None)    
    DSS_Tools.postprend_last_value_on_ts(output_dss_file,out_rec,1)

    #flow_records = [obs_dss_file + "::/USGS SACRAMENTO R/11377100 AB BEND BRDIGE/FLOW//1DAY/USGS/",
    flow_records = [output_dss_file + "::" + in_rec,
                    obs_dss_file + "::/USGS CLEAR CR/11372000 NR IGO/FLOW//1DAY/USGS/",
                    "/MR Sac.-Sac River/11370700 ACID-Dly Flow/Flow//1Day/237.60.125.2.1/",
                    "/MR Sac.-Sac River/11374000 Cow Creek-Dly Flow/Flow//1Day/237.61.125.2.1/",
                    "/MR Sac.-Sac River/11376000 Cottonwood-Dly Flow/Flow//1Day/237.62.125.2.1/",
                    "/MR Sac.-Sac River/11376550 Battle Creek-Dly Flow/Flow//1Day/237.63.125.2.1/",
                    output_dss_file+"::/MR Sac.-Keswick Res./KES-Dam Total Release/Flow//1Day/234.1.125.1.1/",]
    out_rec = "/SACRAMENTO R/BEND BR BALANCE FLOW/FLOW//1DAY/ResSim_PreProcess/"
    DSS_Tools.add_or_subtract_flows(currentAlternative, rtw, flow_records, hydro_dss,
                              [True, False, False, False, False, False, False],
                              out_rec, output_dss_file, what="flow",prepend_n=2)

def compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file):
    # add PG flows 1-5 to create PG_SUM record used by ResSim/TCD scripts
    # Hmm someone has done that for us, hooray!

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

    combine_shasta_gates_flows(currentAlternative,rtw,hydro_dss,output_dss_file)


def compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file):

    # Generate Flow records needed for plotting - mostly, sums of dam outlet flows, if proper sums do not exist
    # Trinity: calc'd total outflow is included in DMS template
    # Keswick: outlet gase only, no summing needed
    # Whiskeytown: outlet + spill = total dam outflow
    inflow_records = ['/MR Sac.-Whiskeytown Lake/WHI-Outlet Release/Flow//1Hour/233.14.125.2.1/',
                      '/MR Sac.-Whiskeytown Lake/WHI-Spill Release/Flow//1Hour/233.14.125.5.1/']
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


def combine_shasta_gates_flows(currentAlt,timewindow,hydro_dss,output_dss_file):

    # construct and add gate openings - old record names
    # Example record:
    # /MR Sac.-Shasta Lake/SHA-TCD-Bottom Gate 1/Gate Position//1Hour/230.13.223.B1.1/ --- old rec name
    #g_rec_1 = "/MR Sac.-Shasta Lake/SHA-TCD-"
    #g_rec_2 = "/Gate Position//1Hour/230.13.223."
    #for level in ['Top','Middle','Bottom','Side']:
    #    ngate = 2 if level=='Side' else 5
    #    gate_recs = []
    #    for ng in range(1,ngate+1):
    #         sg = str(ng)
    #         gate_recs.append(g_rec_1+level+" Gate "+sg+g_rec_2+level[0]+sg+".1/")
    #    out_rec = "/MR Sac.-Shasta Lake/SHA-TCD-"+level+" Gate Sum/count-gate//1Hour/Derived/"
    #    DSS_Tools.add_or_subtract_flows(currentAlt, timewindow, gate_recs, hydro_dss, 
    #                   [True for i in range(ngate)],
    #                   out_rec, output_dss_file, what="n/a")

    # construct and add gate openings
    # Example record:
    # /MR Sac.-Shasta Lake/SHA-TCD-Gate Position 1_B/Gate Position 1//1Hour/230.13.224.B.1/
    g_rec_1 = "/MR Sac.-Shasta Lake/SHA-TCD-Gate Position "
    g_rec_2 = "/Gate Position "
    g_rec_3 = "//1Hour/230.13.224."
    for level in ['Top','Middle','Bottom','Side']:
        ngate = 2 if level=='Side' else 5
        gate_recs = []
        for ng in range(1,ngate+1):
             sg = str(ng)
             lg = level[0]
             gate_recs.append(g_rec_1 + sg+"_"+lg + g_rec_2 + sg + g_rec_3 + lg+"."+sg)
        out_rec = "/MR Sac.-Shasta Lake/SHA-TCD-"+level+" Gate Sum/count-gate//1Hour/Derived/"
        DSS_Tools.add_or_subtract_flows(currentAlt, timewindow, gate_recs, hydro_dss, 
                       [True for i in range(ngate)],
                       out_rec, output_dss_file, what="n/a")

    # Combine river outlet flows
    # Example record:
    # "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G1/Flow//1Hour/230.11.125.14.1/"

    # just writing these out - seems a bit confusing to code them like above
    ro_flows = {
    '750': [
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G1/Flow//1Hour/230.11.125.14.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G2/Flow//1Hour/230.11.125.15.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G3/Flow//1Hour/230.11.125.16.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 G4/Flow//1Hour/230.11.125.17.1/",
    ],
    '850': [
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G1/Flow//1Hour/230.11.125.18.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G2/Flow//1Hour/230.11.125.19.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G3/Flow//1Hour/230.11.125.20.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G4/Flow//1Hour/230.11.125.21.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G5/Flow//1Hour/230.11.125.22.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G6/Flow//1Hour/230.11.125.23.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G7/Flow//1Hour/230.11.125.24.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 G8/Flow//1Hour/230.11.125.25.1/",
    ],
    '950' : [
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G1/Flow//1Hour/230.11.125.26.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G2/Flow//1Hour/230.11.125.27.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G3/Flow//1Hour/230.11.125.28.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G4/Flow//1Hour/230.11.125.29.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G5/Flow//1Hour/230.11.125.30.1/",
        "/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 G6/Flow//1Hour/230.11.125.31.1/",
    ]}
    for level,gate_recs in ro_flows.items():
        out_rec = "/MR Sac.-Shasta Lake/SHA-Outlet Flow "+level+" Sum/Flow//1Hour/Derived/"
        DSS_Tools.add_or_subtract_flows(currentAlt, timewindow, gate_recs, hydro_dss,
                              [True for i in range(len(gate_recs))],
                              out_rec, output_dss_file, what="flow")
    
def leakage(gen_flow,nMidGates,nLowGates,nSideGates,WSE,pre2010=False):
    if pre2010:
        leakage_flow = gen_flow * 0.2
        LKG1 = (13.09 / 100) * leakage_flow
        LKG2 = ((8.05 + (5 - nMidGates) / 5 * 11.65) / 100) * leakage_flow
        LKG3 = (9.34 + (2 - nSideGates) / 2 * 3.31) / 100 * leakage_flow
        LKG4 = (1.03 + (5 - nLowGates) / 5 * 10.01) / 100 * leakage_flow
        LKG5 = (3.84 + 31.12) / 100 * leakage_flow
        LKG6 = (1.79 + 6.77) / 100 * leakage_flow
    else:
        leakage_flow = gen_flow * 0.1606
        LKG1 = 16.3 / 100 * leakage_flow
        LKG2 = 0.0
        LKG3 = (11.63 + (2 - nSideGates) / 2 * 4.12) / 100 * leakage_flow
        LKG4 = (1.28 + (5 - nLowGates) / 5 * 12.47) / 100 * leakage_flow
        LKG5 = (4.78 + 38.76) / 100 * leakage_flow
        LKG6 = (2.23 + 8.44) / 100 * leakage_flow

    # don't leak if leakage withdraw elevation is out of the water
    if 945 <= WSE <1000:
        LKG1 *= (WSE-945)/(1000-945)
    elif 900 <= WSE < 945:
        LKG1 = 0.0
        LKG2 *= (WSE-900)/(945-900)
    elif 831 <= WSE < 900:
        LKG1 = 0.0
        LKG2 = 0.0
        LKG3 *= (WSE-831)/(900-831)

    return [LKG1,LKG2,LKG3,LKG4,LKG5,LKG6]


def get_w2_withdraw_points(nUpperGates,nMidGates,nLowGates,nSideGates,WSE,withdraw_elevs):
    '''Returns a list of withdraw points in the format ['TCDU1','TCDU2','TCDU3'] base on
    open gates and water surface elevation [in feet] where W2 by default withdraws from
    the Shasta TCD'''

    gates = ['TCDU','TCDM','TCDL','TCDS']
    n_gates_open = [nUpperGates,nMidGates,nLowGates,nSideGates]

    w_points = []
    for gate,n in zip(gates,n_gates_open):
        if n > 0 and not w_points:
            for wl,elev in withdraw_elevs[gate].items():
                if elev < WSE - 3: # withdraws must be at least 3 ft below WSE
                    w_points.append(gate+str(wl))

    if not w_points:
        # Above failed, usually because gates are erroneously only open above WSE.
        # In this case, pick highest level that allows flow (without considering if gates are open or not)
        for gate in gates:
            if not w_points:
                for wl, elev in withdraw_elevs[gate].items():
                    if elev < WSE - 3: # withdraws must be at least 3 ft below WSE
                        w_points.append(gate + str(wl))

    if not w_points:
        # there is some bad data
        print("get_w2_withdraw_points fail!")
        print("Gates:",nUpperGates,nMidGates,nLowGates,nSideGates)
        print("WSE:",WSE)
        sys.exit(-1)        

    return w_points

def repeat_annual_daily_over_timewindow(timewindow,annual_dss_file,annual_dss_rec,output_dss_file,out_rec):
    start_year = HecTime(timewindow.getStartTimeString()).year()
    end_year = HecTime(timewindow.getEndTimeString()).year()
    annual_tsc = DSS_Tools.dss_read_ts_safe(annual_dss_file,annual_dss_rec) # must be 365 values
    values = [annual_tsc.values[i] for i in range(annual_tsc.numberValues)] # convert to list, dumb

    day_interval = 1440

    values_out = []
    times_out = []
    for i,y in enumerate(range(start_year,end_year+1)):
        values_out += values
        ht = HecTime("01Jan%i"%y,"2400").value()
        ytimes = [ht+i*day_interval for i in range(len(values))]
        times_out += ytimes        
        if HecTime.isLeap(y):
            values_out.append(values_out[-1])
            times_out.append(times_out[-1]+day_interval)

    tsc_result = TimeSeriesContainer()
    tsc_result.fullName = out_rec
    tsc_result.units = annual_tsc.units
    tsc_result.type = annual_tsc.type
    tsc_result.interval = day_interval
    tsc_result.times = times_out
    tsc_result.values = values_out
    tsc_result.startTime = times_out[0]
    tsc_result.numberValues = len(values_out)

    dss_outfile = HecDss.open(output_dss_file)
    dss_outfile.put(tsc_result)
    dss_outfile.close()


def W2_shasta_TCD_flow(timewindow,hydro_dss,output_dss_file):

    # W2 TCD point sink elevations - hard code those here
    withdraw_elevs = {'TCDU':{1:1042.0,
                          2:1021.0,
                          3:1000.0},
                  'TCDM':{1:942.0,
                          2:921.0,
                          3:900.0},
                  'TCDL':{1:830.0,
                          2:816.0,
                          3:802.0},
                  'TCDS':{1:800.0,
                          2:760.0,
                          3:720.0}
                   }

    # load TCD gates summary from DSS, make regular hourly if needed and fill
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    tsc_u_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Top Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)
    tsc_m_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Middle Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)
    tsc_l_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Bottom Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)
    tsc_s_gates = DSS_Tools.dss_read_ts_safe(output_dss_file,'/MR Sac.-Shasta Lake/SHA-TCD-Side Gate Sum/count-gate//1Hour/Derived/',starttime_str,endtime_str)

    # get year of all records to determine what leakage to use 
    hectimes = DSS_Tools.hectimes_from_tsc(tsc_u_gates)
    years = [hectimes[i].year() for i in range(len(hectimes))]

    # get WSE data, convert to ft if neccessary
    tsc_wse = DSS_Tools.dss_read_ts_safe(hydro_dss,'/MR Sac.-Shasta Lake/SHA-Elevation/Elev//1Hour/230.11.145.2.1/',starttime_str,endtime_str)
    if tsc_wse.units.lower() != 'ft':
        if tsc_wse.units.lower() == 'm':
            for i in range(tsc_wse.numberValues):
               tsc_wse.values[i] = tsc_wse.values[i] * 3.28084 # so lame we have to do this
            tsc_wse.units = 'ft'
        else:
            print('W2_shasta_TCD_flow: elevation units not understood:',tsc_wse.units)
            sys.exit(-1)

    # get generation flow data, also need hourly, only make calcs where we
    # have both gate and flow data and WSE data
    tsc_gen_flow = DSS_Tools.dss_read_ts_safe(hydro_dss,'/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/',starttime_str,endtime_str)

    # move data to lists
    WSE = tsc_wse.values
    nUpperGates = tsc_u_gates.values
    nMidGates = tsc_m_gates.values
    nLowerGates = tsc_l_gates.values
    nSideGates = tsc_s_gates.values
    gen_flow = tsc_gen_flow.values

    # setup results dict
    withdraws_in_order = ['TCDU1', 'TCDU2', 'TCDU3', 'TCDM1', 'TCDM2', 'TCDM3', 'TCDL1', 'TCDL2', 'TCDL3', 'TCDS1', 'TCDS2', 'TCDS3',
                          'LKG1', 'LKG2', 'LKG3', 'LKG4', 'LKG5', 'LKG6', 'TCD_down']
    TCD_flows = {}
    for w in withdraws_in_order:
        TCD_flows[w] = []

    # init flow at each of 4 levels to 0

    # 1) assign all flow to highest gate level open

    # 2) divide flow at each level by 3 and assign as the point source
    # sinks for each level. At this point, all flow is still on the highest
    # open level

    # 3) - If WSE-3 (ft) is less than the withdraw sink, record the flow.
    #     - Then test if that flow is greater than zero
    #     - Find all the withdraw points with flow again above WSE, and drop them from
    #     the list of withdraw points with flow.
    #     - If there are still withdraw points (i.e. the whole gate is not out of the
    #       water, assing all flow to remaining 1 or 2 withdraw points
    #     - Else
    #       go to next level, and figure out how many withdraw points in the NEXT
    #       active level are below water, using a funny test of the 1st member of the 
    #       NEXT NEXT level (4th index of first_next_gate, which should be first_next_sink)
    #      *** should the index to first_next_gate actually be 3?
    for i in range(len(gen_flow)):
        #print(dates[i],' U:',nUpperGates[i],' M:',nMidGates[i],' L:',nLowerGates[i],' S:',nSideGates[i],WSE[i])
        wp = get_w2_withdraw_points(nUpperGates[i],nMidGates[i],nLowerGates[i],nSideGates[i],WSE[i],withdraw_elevs)
        leak = leakage(gen_flow[i],nMidGates[i],nLowerGates[i],nSideGates[i],WSE[i],pre2010=years[i]<2010)

        step_flows = {w:0.0 for w in withdraws_in_order}
        sum_leak = sum(leak)
        step_flows['TCD_down'] = 0.35 * (gen_flow[i]-sum_leak) if nSideGates[i] > 0 else 0.0
        wp_flow = (gen_flow[i] - sum_leak - step_flows['TCD_down']) / len(wp)
        for p in wp:
            step_flows[p] = wp_flow

        for j in range(6): # six leakage withddraw points
            step_flows['LKG'+str(j+1)] = leak[j]

        total_flow = 0.0
        for p in withdraws_in_order:
            total_flow += step_flows[p]
            TCD_flows[p].append(step_flows[p])

        if int(round(total_flow)) != int(round(gen_flow[i])):
            print('DEBUG W2 TCD FLOWS-----------------------------------------------------')
            print('WSE: ', WSE[i])
            print(nUpperGates[i],nMidGates[i],nLowerGates[i],nSideGates[i])
            print(wp)
            print(i,total_flow,gen_flow[i])
            print(step_flows)

    dss_outfile = HecDss.open(output_dss_file)
    for p,flow in TCD_flows.items():
        # re-use gen_flow tsc to write TCD flows (including original units)
        tsc_gen_flow.fullName = "/MR Sac.-Shasta Lake/SHA-W2-TCD-"+p+"/Flow//1Hour/Derived/"
        tsc_gen_flow.values = flow
        dss_outfile.put(tsc_gen_flow)
    dss_outfile.close()


def preprocess_W2_5Res(currentAlternative, computeOptions):
    rtw = computeOptions.getRunTimeWindow()
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    output_dss_file = os.path.join(shared_dir,'DMS_SacTrn_ResSim_Pre-Process.dss') 

    currentAlternative.addComputeMessage('Rectifying units - this may take a while if the length of DMS data is large...')
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
    W2_shasta_TCD_flow(rtw,hydro_dss,output_dss_file) # depends on compute_5Res_outflows having completed

    # link these flows for W2, to avoid zero-dam-flow situations - needs to be above 2.0 cfs for the flowweightaverage script to recognize it :/
    DSS_Tools.min_ts(output_dss_file, '/MR Sac.-Trinity Lake/TRN-GenerationG1G2_Sum/Flow//1Hour/ResSim_PreProcess/', 1.1, output_dss_file, 'min_flow')
    DSS_Tools.min_ts(output_dss_file, '/MR Sac.-Lewiston Res./LEW-Gen_plus_Outlet Release/Flow//1Hour/ResSim_PreProcess/', 1.1, output_dss_file, 'min_flow')

    # pit river is recorded at 12:00 daily, W2 plugin can't handle that, shift backward 12 hours
    DSS_Tools.shift_pit_river_time(hydro_dss,"/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Day/230.6.125.1.1/",
                         output_dss_file,"/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Day/Shifted-12H/",
                         start_date=None,end_date=None)

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

    currentAlternative.addComputeMessage('Rectifying units - this may take a while if the length of DMS data is large...')
    hydro_dss = os.path.join(shared_dir, 'DMS_SacTrnHydroTS.dss')
    fix_DMS_types_units(hydro_dss)
    met_dss_file = os.path.join(shared_dir,'DMS_SacTrnMet.dss')
    fix_DMS_types_units(met_dss_file)
    # ressim can't handle different units under model linking
    standardize_bc_temp_water_to_C(hydro_dss,output_dss_file)

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

    repeat_annual_daily_over_timewindow(rtw,os.path.join(shared_dir,'PatternHydro.dss'),
                                        "//SPRING_CREEK_DEBRIS_DAM/TEMP-WATER//1Day/ANNUAL_TEMPLATE/",
                                        output_dss_file,"//SPRING_CREEK_DEBRIS_DAM/TEMP-WATER//1Day/ANNUAL/")
                        
    # if template IDs exist still, remove them
    #DSS_Tools.strip_templateID_and_rename_records(hydro_dss,currentAlternative)
    #DSS_Tools.strip_templateID_and_rename_records(met_dss_file,currentAlternative)

    splice_lewiston_met_data(currentAlternative, rtw, met_dss_file, output_dss_file,months=[1,2,3])
    compute_5Res_outflows(currentAlternative, rtw, hydro_dss, output_dss_file)
    compute_plotting_records(currentAlternative, rtw, hydro_dss, output_dss_file)    
    compute_river_balance_flows(currentAlternative, rtw, hydro_dss, 
        os.path.join(shared_dir,"WTMP_SacTrn_Historical.dss"), output_dss_file)  # depends on WHI dam flow, created in compute_plotting_records

    # calculate meteorological airtemp lapse for the elevation @ Shasta Lake
    currentAlternative.addComputeMessage('lapse infile: '+met_dss_file)
    currentAlternative.addComputeMessage('lapse outfile: '+output_dss_file)
    DSS_Tools.airtemp_lapse(met_dss_file, "/MR SAC.-CLEAR CR. TO SAC R./KRDD-AIR TEMPERATURE/TEMP-AIR//1HOUR/235.40.53.1.1/",
                  -0.7, output_dss_file, "Shasta_Lapse")

    DSS_Tools.relhum_from_at_dp(met_dss_file,
                      "/MR Sac.-Clear Cr. to Sac R./KRDD-Air temperature/Temp-Air//1Hour/235.40.53.1.1/",
                      "/MR Sac.-Clear Cr. to Sac R./KRDD-Dew Point/Temp-DewPoint//1Hour/235.40.51.1.1/")
    DSS_Tools.relhum_from_at_dp(met_dss_file,
                      "/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Air temperature/Temp-Air//1Hour/232.6.53.1.1/",
                      "/MR Sac.-Lewiston Res./TCAC1 - Calc Data-Dew Point/Temp-DewPoint//1Hour/232.6.51.1.1/")
    DSS_Tools.dp_from_at_relhum(met_dss_file,
                      "/MR SAC.-CLEAR CR. TO SAC R./KRDD-AIR TEMPERATURE/TEMP-AIR//1Hour/Shasta_Lapse/",
                      "/MR Sac.-Clear Cr. to Sac R./KRDD-/RELHUM-FROM-AT-DP//1Hour/235.40.53.1.1-DERIVED/")
   
    

    return True


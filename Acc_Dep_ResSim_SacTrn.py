'''
Created on 8/7/2023
@author: Scott Burdick-Yahya
@organization: Resource Management Associates
@contact: scott@rmanet.com
@note:
'''

import create_balance_flow_jython as cbfj
reload(cbfj)
from com.rma.model import Project
import os
import Simple_DSS_Functions as sdf
reload(sdf)

def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()

    # Shasta Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values
    
    run_dir = computeOptions.getRunDirectory()
    project_dir = Project.getCurrentProject().getProjectDirectory()
    currentAlternative.addComputeMessage('project_dir: ' + project_dir)
    currentAlternative.addComputeMessage('run dir: ' + run_dir)
    balance_period_str = currentAlternative.getTimeStep()
    shared_dir = os.path.join(project_dir, 'shared')

    DMS_hydro_dss_file = os.path.join(shared_dir, "DMS_SacTrnHydroTS.dss")
    output_dss_file = os.path.join(shared_dir,'DMS_SacTrn_ResSim_Pre-Process.dss')
    fallback_dss_file = os.path.join(shared_dir,'WTMP_SacTrn_Historical.dss')

    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/Sacramento R. a Delta-Flow/Flow//1Day/230.9.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/McCloud River-Flow/Flow//1Day/230.8.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Day/230.6.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR SAC.-SHASTA LAKE/SULANHARAS CREEK-FLOW/Flow//1Day/230.7.125.1.1/',rtw,output_dss_file,'1HOUR')

    inflow_records = ['::'.join([output_dss_file,'/MR Sac.-Shasta Lake/McCloud River-Flow/Flow//1Hour/230.8.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/McCloud River-Flow/Flow//1Hour/230.8.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Hour/230.6.125.1.1/']),
                      '::'.join([output_dss_file,'/MR SAC.-SHASTA LAKE/SULANHARAS CREEK-FLOW/Flow//1Hour/230.7.125.1.1/']),]

    outflow_records = ['/MR Sac.-Shasta Lake/SHA-Spill Release/Flow//1Hour/230.11.125.4.1/',
                       '::'.join([fallback_dss_file,'/USBR/SHASTA_RRU/FLOW//1HOUR/SUPP/']),
                       '::'.join([fallback_dss_file,'/USBR/SHASTA_RRM/FLOW//1HOUR/SUPP/']),
                       '::'.join([fallback_dss_file,'/USBR/SHASTA_RRL/FLOW//1HOUR/SUPP/']),
                       '/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/'
                       ]

    stage_record = '::'.join([fallback_dss_file,'/USBR/SHASTA/ELEVATION//1HOUR/USBR_BLESSED/'])
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_shasta.csv'), 'Shasta') #TODO: check this    

    use_conic = False
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/SHASTA RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/SHASTA RESERVOIR/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/SHASTA RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/SHASTA RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/SHASTA RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    cbfj.create_balance_flows(currentAlternative, rtw, 'Shasta', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')


    # Keswick Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values

    inflow_records = ['::'.join([DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/']),
                      '::'.join([DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/SHA-Spill Release/Flow//1Hour/230.11.125.4.1/']),
                      '/USBR/SHASTA_RRU/FLOW//1HOUR/SUPP/',
                      '/USBR/SHASTA_RRM/FLOW//1HOUR/SUPP/',
                      '/USBR/SHASTA_RRL/FLOW//1HOUR/SUPP/',                      
                      '::'.join([DMS_hydro_dss_file,'/MR SAC.-WHISKEYTOWN LAKE/WHI-GENERATION RELEASE/FLOW//1HOUR/233.14.125.1.1/']),
                      '/SPC DEBRIS DAM FLOW CALCULATION/SPRING_CREEK_SPC2017SHIFT-WHI_GEN_RELEASE-NONEG/FLOW//1HOUR/USBR_DERIVED/'
                      ]

    outflow_records = ['/USBR/KESWICK_QOUT_SUM/FLOW//1HOUR/USBR_DERIVED/']

    stage_record = '/USBR/KESWICK/ELEVATION//1HOUR/USBR_BLESSED/'
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_keswick.csv'), 'Keswick') #TODO: check this

    use_conic = False
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/KESWICK RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/KESWICK RESERVOIR/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/KESWICK RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/KESWICK RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/KESWICK RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    cbfj.create_balance_flows(currentAlternative, rtw, 'Keswick', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, fallback_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=60*3, alt_period_string='3Hour')

    # Lewiston Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values

    inflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/231.5.125.8.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G3/Flow//1Hour/231.5.125.9.1/',
                      '/MR Sac.-Trinity Lake/TRN-Spill Release/Flow//1Hour/231.5.125.3.1/']

    outflow_records = ['/MR Sac.-Whiskeytown Lake/JCR-Generation Release/Flow//1Hour/233.13.125.1.1/', #lewiston  CC diversion tunnel: Clear Creek Transfer operation
                       '/MR Sac.-Lewiston Res./LEW-Fish Hatchery Release/Flow//1Hour/232.12.125.1.1/',
                       '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/232.12.125.3.1/',
                       '/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/232.12.125.2.1/',
                       '/MR Sac.-Lewiston Res./LEW-Spill Release Hrly/Flow//1Hour/232.12.125.5.1/', #Lewiston Res. Dam at Trinity River - Powerhouse
                       ]

    stage_record = '::'.join([fallback_dss_file,'/USBR_CLEANED_FULLLINEARINTERP/LEW_ELEVATION_PROCESSED/ELEV//1HOUR/GRAB2/'])
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_lewiston.csv'), 'lewiston') #TODO: check this

    use_conic = False
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/LEWISTON RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/LEWISTON RESERVOIR/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    cbfj.create_balance_flows(currentAlternative, rtw, 'LEWISTON', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=60*6, alt_period_string='6Hour')

    # Trinity Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values

    # resample the Trinity inflow data to hourly
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Trinity Lake/EF Trinity River-Inflow/Flow//1Day/231.7.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Trinity Lake/Stuart Fork-Inflow/Flow//1Day/231.8.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Trinity Lake/Swift Creek-Inflow/Flow//1Day/231.9.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Trinity Lake/Trinity River-Inflow/Flow//1Day/231.6.125.1.1/',rtw,output_dss_file,'1HOUR')

    inflow_records = ['::'.join([output_dss_file,'/MR Sac.-Trinity Lake/EF Trinity River-Inflow/Flow//1Hour/231.7.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Trinity Lake/Stuart Fork-Inflow/Flow//1Hour/231.8.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Trinity Lake/Swift Creek-Inflow/Flow//1Hour/231.9.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Trinity Lake/Trinity River-Inflow/Flow//1Hour/231.6.125.1.1/']),]

    outflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/231.5.125.8.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G3/Flow//1Hour/231.5.125.9.1/',
                      '/MR Sac.-Trinity Lake/TRN-Spill Release/Flow//1Hour/231.5.125.3.1/']

    stage_record = '::'.join([fallback_dss_file,'/USBR_RECLEANED_FULLLINEARINTERP/TRN_ELEVATION/ELEV//1HOUR/GRAB2/'])
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_trinity.csv'), 'trinity') #TODO: check this

    use_conic = False
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/TRINITY RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/TRINITY RESERVOIR/CONIC_FROM_ELEV/STORAGE//1HOUR/DERIVED/"
    output_dss_record_name = "/TRINITY RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/TRINITY RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/TRINITY RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"


    cbfj.create_balance_flows(currentAlternative, rtw, 'Trinity', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')

    # whiskeytown Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values
    
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR SAC.-WHISKEYTOWN LAKE/CLEAR CREEK AB WHISKEYTOWN LAKE-FLOW/FLOW//1DAY/233.131.125.1.1/',rtw,output_dss_file,'1HOUR')
    inflow_records = ['/MR Sac.-Whiskeytown Lake/JCR-Generation Release/Flow//1Hour/233.13.125.1.1/', 
                      '::'.join([output_dss_file,'/MR SAC.-WHISKEYTOWN LAKE/CLEAR CREEK AB WHISKEYTOWN LAKE-FLOW/FLOW//1Hour/233.131.125.1.1/'])
                      ]

    outflow_records = ['/MR SAC.-WHISKEYTOWN LAKE/WHI-OUTLET RELEASE/FLOW//1HOUR/233.14.125.2.1/', #Whiskeytown Lake - Upper Dam Outlet
                       '/MR SAC.-WHISKEYTOWN LAKE/WHI-GENERATION RELEASE/FLOW//1HOUR/233.14.125.1.1/',#Whiskeytown Lake - Spring CreekTunnel diversion
                       '/MR SAC.-WHISKEYTOWN LAKE/WHI-SPILL RELEASE/FLOW//1HOUR/233.14.125.5.1/', #spill release
                       ]

    stage_record = '::'.join([fallback_dss_file,'/USBR+23HrShift/WHI_ELEVATION/ELEV//1HOUR/GRAB2-MANUAL-FILTER-INTERPFILL/'])
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_whiskeytown.csv'), 'Whiskeytown') #TODO: check this

    use_conic = False
    write_evap = False
    write_storage = False

    evap_dss_record_name = "/WHISKEYTOWN RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/WHISKEYTOWN RESERVOIR/CONIC_FROM_ELEV/STORAGE//1HOUR/DERIVED/"
    output_dss_record_name = "/WHISKEYTOWN RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    if use_conic:
        output_dss_record_name = "/WHISKEYTOWN RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/WHISKEYTOWN RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    cbfj.create_balance_flows(currentAlternative, rtw, 'whiskeytown', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')


    #######################################################################################
    # TODO: Calculate River balances
    #######################################################################################
    # Trinity River: Limekiln Gulch

    # Trinity River: Douglas City

    # Trinity River: Junction City

    # Clear Creek at South Fork junction (IGO)

    # Sacramento River at Bend Bridge

    return True

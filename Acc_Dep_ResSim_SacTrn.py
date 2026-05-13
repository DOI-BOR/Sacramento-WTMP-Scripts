
"""
create_balance_flows for WTMP (Sacramento-Trinity System)
=========================================================

Purpose
-------
Prepare and write reservoir balance-flow time series for multiple reservoirs
(Shasta, Keswick, Lewiston, Trinity, Whiskeytown) used by the WTMP workflow.
This script relies on DSS records provided by DMS/WTMP and resamples/derives
inputs as needed, then calls helper routines to compute mass-balance flows.

Two entry points are provided:

- ``computeAlternative``  
  Builds balance flows for Shasta, Keswick, Lewiston, Trinity, and Whiskeytown
  using the runtime window defined by the WAT compute options and writes outputs
  to shared DSS files.

- ``compute_W2_forecast_balance``  
  Builds Lewiston balance flows for a W2 forecast context and writes outputs to
  the forecast DSS file.

Notes
-----
- **Environment:** Jython within HEC-WAT using Python 2.7 syntax. The script 
  expects WAT to supply the current project context and compute options 
  (DSS filename, run directory, runtime window, etc.).
- **Time-step assumptions:**  
  Flows are period-averaged; evaporation is period-accumulated (e.g., feet);
  stage/elevation is instantaneous. Balance period is derived from the WAT
  alternative time step string.
- **DSS paths:**  
  Many inputs/outputs are specified as DSS pathnames; some are referenced with
  a prefix (``file::/A/B/C/D/E/F/``) to read/write from a specific DSS file.

See Also
--------
create_balance_flow_jython
    Provides helpers: ``read_elev_storage_area_file`` (read elevation-storage-area
    lookups) and ``create_balance_flows`` (compute/write reservoir balance flows).
Simple_DSS_Functions
    Provides ``resample_dss_ts`` for resampling DSS time series between time steps.
"""

import os                                  # Standard library: filesystem path joins and file operations

import create_balance_flow_jython as cbfj  # Local helper: balance-flow computation & elev-storage-area readers
reload(cbfj)                               # Jython: ensure latest module version is loaded at runtime

from com.rma.model import Project          # WAT/RMA project API: access current project directory/context

import Simple_DSS_Functions as sdf         # Local DSS utilities: resampling time series across time steps
reload(sdf)                                # Jython: reload to pick up any recent changes


def computeAlternative(currentAlternative, computeOptions):
    """
    Main WAT scripting entry point for balance flows (WTMP simulation context).

    Computes and writes reservoir balance-flow time series for:
    Shasta, Keswick, Lewiston, Trinity, and Whiskeytown.

    Parameters
    ----------
    currentAlternative : object
        WAT scripting alternative instance providing logging and data-location
        helpers (e.g., ``addComputeMessage``, ``getTimeStep``).
    computeOptions : object
        WAT compute options providing runtime window, DSS filenames, and run
        directory (e.g., ``getRunTimeWindow``, ``getRunDirectory``).

    Returns
    -------
    bool
        True on completion (conventional success indicator for WAT scripting).

    Notes
    -----
    - Uses period assumptions consistent with ResSim hydro timestep.
    - Resamples select DMS inflows to hourly before computing balances.
    - Writes derived balance flows and (optionally) derived evaporation/storage.
    """
    
    # Log the start of the script
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')  # Log spacer for readability

    # Get the model DSS output
    dss_file = computeOptions.getDssFilename()  

    # Get the current runtime window
    rtw = computeOptions.getRunTimeWindow()  
    
    # Setup run information for the rest of the script
    run_dir = computeOptions.getRunDirectory()                  # Directory of this run execution
    project_dir = Project.getCurrentProject().getProjectDirectory()  # WAT project directory root

    currentAlternative.addComputeMessage('project_dir: ' + project_dir)  # Diagnostic: project dir
    currentAlternative.addComputeMessage('run dir: ' + run_dir)          # Diagnostic: run dir

    balance_period_str = currentAlternative.getTimeStep()       # Use WAT alt timestep string (e.g., '1Hour')

    shared_dir = os.path.join(project_dir, 'shared')            # Shared folder holding common DSS assets

    DMS_hydro_dss_file = os.path.join(shared_dir, "DMS_SacTrnHydroTS.dss")             # DMS hydro inputs
    output_dss_file = os.path.join(shared_dir,'DMS_SacTrn_ResSim_Pre-Process.dss')     # Pre-process outputs
    fallback_dss_file = os.path.join(shared_dir,'WTMP_SacTrn_Historical.dss')          # Historical fallback DSS

    # Shasta Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values

    # Resample Shasta inflow components from daily to hourly for balance computation
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/Sacramento R. a Delta-Flow/Flow//1Day/230.9.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/McCloud River-Flow/Flow//1Day/230.8.125.1.1/',rtw,output_dss_file,'1HOUR')

    # Pit river is recorded at GMT midnight, and so when requesting the run time windows from DSS, the last day is left off the record.
    # This generates some bad value in the balance flow, which hopefully we catch and make zero, and which hopefully don't effect the run much
    # because it's only the last day.
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Day/230.6.125.1.1/',rtw,output_dss_file,'1HOUR')
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR SAC.-SHASTA LAKE/SULANHARAS CREEK-FLOW/Flow//1Day/230.7.125.1.1/',rtw,output_dss_file,'1HOUR')
        
    # Construct inflow list (file::path syntax ensures a specific DSS file is used)
    inflow_records = ['::'.join([output_dss_file,'/MR Sac.-Shasta Lake/McCloud River-Flow/Flow//1Hour/230.8.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/Sacramento R. a Delta-Flow/Flow//1Hour/230.9.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/Pit R. Branch-Flow/Flow//1Hour/230.6.125.1.1/']),
                      '::'.join([output_dss_file,'/MR SAC.-SHASTA LAKE/SULANHARAS CREEK-FLOW/Flow//1Hour/230.7.125.1.1/']),]

    # Outflows include spill, outlet sums, and generation release; some derived in pre-process DSS
    outflow_records = ['/MR Sac.-Shasta Lake/SHA-Spill Release/Flow//1Hour/230.11.125.4.1/',
                       '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 Sum/Flow//1Hour/Derived/']),
                       '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 Sum/Flow//1Hour/Derived/']),
                       '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 Sum/Flow//1Hour/Derived/']),
                       '/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/'
                       ]

    # Set the stage series
    stage_record = '/MR Sac.-Shasta Lake/SHA-Elevation/Elev//1Hour/230.11.145.2.1/'              # Hourly elevation

    # Set the evaporation series
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])                       # Placeholder evap (zeros)

    # Set the elevation/storage curve
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_shasta.csv'), 'Shasta') #TODO: check this    

    # Set the optional write flags
    use_conic = False     # Use conic interpolation (disabled)
    write_evap = False    # Write derived evaporation (disabled)
    write_storage = False # Write derived storage (disabled)

    # Output DSS record names for optional derived series
    evap_dss_record_name = "/SHASTA RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/SHASTA RESERVOIR/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/SHASTA RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    if use_conic:
        output_dss_record_name = "/SHASTA RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/SHASTA RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # ---- Compute and write Shasta balance flows -----------------------------------------
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

    # Construct inflows for Keswick: upstream Shasta releases, Whiskeytown generation, Spring Creek Debris Dam
    inflow_records = ['::'.join([DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/SHA-Generation Release/Flow//1Hour/230.11.125.3.1/']),
                      '::'.join([DMS_hydro_dss_file,'/MR Sac.-Shasta Lake/SHA-Spill Release/Flow//1Hour/230.11.125.4.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/SHA-Outlet Flow 950 Sum/Flow//1Hour/Derived/']),
                      '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/SHA-Outlet Flow 850 Sum/Flow//1Hour/Derived/']),
                      '::'.join([output_dss_file,'/MR Sac.-Shasta Lake/SHA-Outlet Flow 750 Sum/Flow//1Hour/Derived/']),
                      '::'.join([DMS_hydro_dss_file,'/MR SAC.-WHISKEYTOWN LAKE/WHI-GENERATION RELEASE/FLOW//1HOUR/233.14.125.1.1/']),
                      '::'.join([DMS_hydro_dss_file,'/MR Sac.-Keswick Res./Spring Creek Debris Dam (SPR)-Dam Total Release/Flow//1Hour/234.2.125.1.1/']),
                      ]

    # Set the outflow series
    outflow_records = ['::'.join([DMS_hydro_dss_file,'/MR Sac.-Keswick Res./KES-Dam Total Release/Flow//1Hour/234.1.125.1.1/'])]

    # Set the stage series
    stage_record = '::'.join([DMS_hydro_dss_file,'/MR Sac.-Keswick Res./KES-Reservoir Elevation/Elev//1Hour/234.1.145.1.1/'])  # Hourly elevation
    
    # Set the evaporation series
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])                                                    # Placeholder evap (zeros)

    # Set the elevation/storage curve
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_keswick.csv'), 'Keswick') #TODO: check this

    # Set the optional write flags
    use_conic = False     # Interpolation option disabled
    write_evap = False    # Do not write derived evaporation
    write_storage = False # Do not write derived storage

    # Output DSS record names specific to Keswick
    evap_dss_record_name = "/KESWICK RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/KESWICK RESERVOIR/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/KESWICK RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    if use_conic:
        output_dss_record_name = "/KESWICK RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/KESWICK RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # ---- Compute and write Keswick balance flows ----------------------------------------
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

    # Trinity upstream releases to Lewiston (inflows)
    inflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/231.5.125.8.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G3/Flow//1Hour/231.5.125.9.1/',
                      '/MR Sac.-Trinity Lake/TRN-Spill Release/Flow//1Hour/231.5.125.3.1/']

    # Lewiston outflows: tunnel generation/diversions, hatchery, generation, outlet, spill
    outflow_records = ['/MR Sac.-Whiskeytown Lake/JCR-Generation Release/Flow//1Hour/233.13.125.1.1/', # lewiston CC diversion tunnel: Clear Creek Transfer operation
                       '/MR Sac.-Lewiston Res./LEW-Fish Hatchery Release/Flow//1Hour/232.12.125.1.1/',
                       '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/232.12.125.3.1/',
                       '/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/232.12.125.2.1/',
                       '/MR Sac.-Lewiston Res./LEW-Spill Release Hrly/Flow//1Hour/232.12.125.5.1/', # Lewiston Res. Dam at Trinity River - Powerhouse
                       ]

    # Set the stage timeseries
    stage_record = '/MR Sac.-Lewiston Res./LEW-Reservoir Elevation Hrly/Elev//1Hour/232.12.145.1.1/'   # Hourly elevation

    # Set the evaporation timeseries
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])                            # Placeholder evap (zeros)

    # Set the elevation/storage curve
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_lewiston.csv'), 'lewiston') #TODO: check this

    # Set the optional write flags
    use_conic = False     # Interpolation option disabled
    write_evap = False    # Do not write derived evaporation
    write_storage = False # Do not write derived storage

    # Output DSS record names specific to Lewiston
    evap_dss_record_name = "/LEWISTON RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/LEWISTON RESERVOIR/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    if use_conic:
        output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # ---- Compute and write Lewiston balance flows ---------------------------------------
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
        
    # Construct resampled inflow records for Trinity
    inflow_records = ['::'.join([output_dss_file,'/MR Sac.-Trinity Lake/EF Trinity River-Inflow/Flow//1Hour/231.7.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Trinity Lake/Stuart Fork-Inflow/Flow//1Hour/231.8.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Trinity Lake/Swift Creek-Inflow/Flow//1Hour/231.9.125.1.1/']),
                      '::'.join([output_dss_file,'/MR Sac.-Trinity Lake/Trinity River-Inflow/Flow//1Hour/231.6.125.1.1/']),]

    # Trinity outflows: generation, outlets (G1–G3), and spill
    outflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/231.5.125.8.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G3/Flow//1Hour/231.5.125.9.1/',
                      '/MR Sac.-Trinity Lake/TRN-Spill Release/Flow//1Hour/231.5.125.3.1/']

    # Set the stage timeseries
    stage_record = '/MR Sac.-Trinity Lake/TRN-Reservoir Elevation/Elev//1Hour/231.5.145.1.1/'  # Hourly elevation

    # Set the evaporation timeseries
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])                    # Placeholder evap (zeros)

    # Set the elevation/storage curve
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_trinity.csv'), 'trinity') #TODO: check this

    # Set the optional write flags
    use_conic = False     # Interpolation option disabled
    write_evap = False    # Do not write derived evaporation
    write_storage = False # Do not write derived storage

    # Output DSS record names specific to Trinity
    evap_dss_record_name = "/TRINITY RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/TRINITY RESERVOIR/CONIC_FROM_ELEV/STORAGE//1HOUR/DERIVED/"
    output_dss_record_name = "/TRINITY RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    if use_conic:
        output_dss_record_name = "/TRINITY RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/TRINITY RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # ---- Compute and write Trinity balance flows ----------------------------------------
    cbfj.create_balance_flows(currentAlternative, rtw, 'Trinity', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=1440, alt_period_string='1Day')

    # Whiskeytown Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values
    
    # Resample Clear Creek above Whiskeytown inflow from daily to hourly
    sdf.resample_dss_ts(DMS_hydro_dss_file,'/MR SAC.-WHISKEYTOWN LAKE/CLEAR CREEK AB WHISKEYTOWN LAKE-FLOW/FLOW//1DAY/233.131.125.1.1/',rtw,output_dss_file,'1HOUR')

    inflow_records = ['/MR Sac.-Whiskeytown Lake/JCR-Generation Release/Flow//1Hour/233.13.125.1.1/', 
                      '::'.join([output_dss_file,'/MR SAC.-WHISKEYTOWN LAKE/CLEAR CREEK AB WHISKEYTOWN LAKE-FLOW/FLOW//1Hour/233.131.125.1.1/'])
                      ]

    # Whiskeytown outflows: outlet, generation (tunnel diversion), spill
    outflow_records = ['/MR SAC.-WHISKEYTOWN LAKE/WHI-OUTLET RELEASE/FLOW//1HOUR/233.14.125.2.1/', # Whiskeytown Lake - Upper Dam Outlet
                       '/MR SAC.-WHISKEYTOWN LAKE/WHI-GENERATION RELEASE/FLOW//1HOUR/233.14.125.1.1/',# Whiskeytown Lake - Spring CreekTunnel diversion
                       '/MR SAC.-WHISKEYTOWN LAKE/WHI-SPILL RELEASE/FLOW//1HOUR/233.14.125.5.1/',     # spill release
                       ]

    # Set the stage timeseries
    stage_record = '/MR Sac.-Whiskeytown Lake/WHI-Reservoir Elevation/Elev//1Hour/233.14.145.1.1/'    # Hourly elevation
    
    # Set the evaporation timeseries
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])                           # Placeholder evap (zeros)

    # Set the elevation/storage curve
    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_whiskeytown.csv'), 'Whiskeytown') #TODO: check this

    # Set the optional write flags
    use_conic = False     # Interpolation option disabled
    write_evap = False    # Do not write derived evaporation
    write_storage = False # Do not write derived storage

    # Output DSS record names specific to Whiskeytown
    evap_dss_record_name = "/WHISKEYTOWN RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/WHISKEYTOWN RESERVOIR/CONIC_FROM_ELEV/STORAGE//1HOUR/DERIVED/"
    output_dss_record_name = "/WHISKEYTOWN RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    if use_conic:
        output_dss_record_name = "/WHISKEYTOWN RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
        
        if 'ZEROS' in evap_record:
            output_dss_record_name = "/WHISKEYTOWN RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # ---- Compute and write Whiskeytown balance flows -------------------------------------
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


def compute_W2_forecast_balance(currentAlternative, computeOptions):
    """
    W2 forecast balance-flow computation (Lewiston only).

    Computes and writes **Lewiston** reservoir balance flows to the forecast DSS
    file (``WTMP_SacTrn_Forecast.dss``). This entry point mirrors the Lewiston
    portion of ``computeAlternative`` but targets the forecast context.

    Parameters
    ----------
    currentAlternative : object
        WAT scripting alternative instance used for logging and configuration.
    computeOptions : object
        WAT compute options providing runtime window and directories.

    Returns
    -------
    bool
        True on completion.

    Notes
    -----
    - Uses zero-evap placeholder path for simplicity (can be replaced if needed).
    - Time-step assumptions follow ResSim hydro timestep conventions.
    """
    
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName())
    currentAlternative.addComputeMessage('\n')  # Spacer for logs

    dss_file = computeOptions.getDssFilename()    # Current model DSS (not directly used for forecast writes)
    rtw = computeOptions.getRunTimeWindow()       # Forecast runtime window

    # Shasta Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values
    
    run_dir = computeOptions.getRunDirectory()                  # Forecast run directory
    project_dir = Project.getCurrentProject().getProjectDirectory()  # Project directory

    currentAlternative.addComputeMessage('project_dir: ' + project_dir)  # Diagnostic: project dir
    currentAlternative.addComputeMessage('run dir: ' + run_dir)          # Diagnostic: run dir

    balance_period_str = currentAlternative.getTimeStep()       # Alt timestep string for balance period
    shared_dir = os.path.join(project_dir, 'shared')            # Shared folder reference

    output_dss_file = os.path.join(shared_dir,'WTMP_SacTrn_Forecast.dss')  # Forecast DSS for outputs

    # Lewiston Inputs **********************************************************************
    # ******* Use same time resolution as ResSim hydro model time step ************
    # Flows are assumed to be period averaged
    # Evap assumed to be period accumulated length (e.g., ft)
    # Stage assumed to be instantaneous values

    # Trinity upstream releases to Lewiston (inflows) for forecast
    inflow_records = ['/MR Sac.-Trinity Lake/TRN-Generation Release/Flow//1Hour/231.5.125.2.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G1/Flow//1Hour/231.5.125.7.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G2/Flow//1Hour/231.5.125.8.1/',
                      '/MR Sac.-Trinity Lake/TRN-Outlet Release G3/Flow//1Hour/231.5.125.9.1/',
                      '/MR Sac.-Trinity Lake/TRN-Spill Release/Flow//1Hour/231.5.125.3.1/']

    # Lewiston outflows in the forecast context
    outflow_records = ['/MR Sac.-Whiskeytown Lake/JCR-Generation Release/Flow//1Hour/233.13.125.1.1/', # lewiston CC diversion tunnel: Clear Creek Transfer operation
                       '/MR Sac.-Lewiston Res./LEW-Fish Hatchery Release/Flow//1Hour/232.12.125.1.1/',
                       '/MR Sac.-Lewiston Res./LEW-Generation Release Hrly/Flow//1Hour/232.12.125.3.1/',
                       '/MR Sac.-Lewiston Res./LEW-Outlet Release Hrly/Flow//1Hour/232.12.125.2.1/',
                       '/MR Sac.-Lewiston Res./LEW-Spill Release Hrly/Flow//1Hour/232.12.125.5.1/', # Lewiston Res. Dam at Trinity River - Powerhouse
                       ]

    stage_record = '/MR Sac.-Lewiston Res./LEW-Reservoir Elevation Hrly/Elev//1Hour/232.12.145.1.1/'  # Hourly elevation (forecast)
    evap_record = '::'.join([output_dss_file,'//ZEROS/FLOW//1HOUR/ZEROS/'])                           # Zero-evap placeholder

    elev_stor_area = cbfj.read_elev_storage_area_file(os.path.join(shared_dir, 'AMR_scratch_lewiston.csv'), 'lewiston') #TODO: check this

    use_conic = False     # Interpolation option disabled
    write_evap = False    # Evaporation write disabled
    write_storage = False # Storage write disabled

    # Forecast output record names (Lewiston)
    evap_dss_record_name = "/LEWISTON RESERVOIR/EVAP FLOW/FLOW//1HOUR/DERIVED/"
    storage_dss_record_name = "/LEWISTON RESERVOIR/STORAGE/FLOW//1HOUR/DERIVED/"
    output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED/"
    
    if use_conic:
        output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP/"
    
    if 'ZEROS' in evap_record:
            output_dss_record_name = "/LEWISTON RESERVOIR/BALANCE FLOW/FLOW//1HOUR/DERIVED-CONIC INTERP NO EVAP/"

    # ---- Compute and write Lewiston forecast balance flows -------------------------------
    cbfj.create_balance_flows(currentAlternative, rtw, 'LEWISTON', inflow_records, outflow_records, stage_record, evap_record,
                                elev_stor_area, DMS_hydro_dss_file, output_dss_record_name, output_dss_file, shared_dir,
                                evap_dss_record_name=evap_dss_record_name, storage_dss_record_name=storage_dss_record_name,
                                balance_period_str=balance_period_str, use_conic=use_conic, write_evap=write_evap, write_storage=write_storage,
                                alt_period=60*6, alt_period_string='6Hour')

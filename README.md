# Sacramento WTMP Scripts
The HEC-WAT scripts folder for the Sacramento River within the Bureau of Reclamation (USBR) Water Temperature Modeling Platform (WTMP).
This repository is a dependency of the Sacramento WTMP Study repository.
These scripts are called at different points in the WTMP workflow.
These scripts are all in Jython, the implementation of Python in Java.

## W2 files
### Forecast
### Hindcast
### Planning
### Shared or not specified
* Lewiston-W2.py
* OutputLink_W2-Lewiston-Downstream_SacTrn.py
* OutputLink_W2-Trinity-W2-Lewiston_SacTrn.py
* OutputLink_W2-Whiskeytown-Downstream_SacTrn.py
* Pre-Process_W2_5-Res_Prescribed_SacTrn.py
  * The pre-processing script for W2
* Pre-Process_W2_5-Res_Prescribed_V45_SacTrn.py
  * The pre-processing script for W2
* Pre-Process_W2_5-Res_SacTrn.py
  * The pre-processing script for W2
* w2_shasta_TCD_flow.py

## ResSim files
### Forecast
### Hindcast
### Planning
### Shared or not specified
* Lewiston_Tunnel_AddHeating_5ResFcst.py
* Lewiston_Tunnel_AddHeating_5ResTemp.py
* Post-Process_ResSim_SacTrn.py
  * The post-processing script for ResSim
* Post-Process_ResSim_SacTrnNS.py
  * The post-processing script for ResSim
* Post-Process_ResSim_SacTrnPO.py
  * The post-processing script for ResSim
* Post-Process_ResSim_SacTrnRv.py
  * The post-processing script for ResSim
* Pre-Process_ResSim_SacTrn.py
  * The pre-processing script for ResSim
* Pre-Process_ResSim_SacTrnPO.py
  * The pre-processing script for ResSim

## Shared files between W2 and ResSim
These files contain various function that are shared by different models.
* Acc_Dep_ResSim_SacTrn.py
* AddTwoTimeSeries.py
* BoundaryFixes.py
* ClearCreekHeating.py
* Combine_2_DSS_Recs.py
* create_balance_flow_jython.py
* DMS_preprocess.py
* DSS_Tools.py
* equilibrium_temp.py
* flowweightaverage.py
* Simple_DSS_Functions.py
* SpringCreekTunnelHeating.py
* Trinity2Lewiston.py
* tz_offset.py
* UpperSac_RSS_CombineTemp.py

### Forecast
Forecast_preprocess.py

## Dependencies
There are dependencies that come from the WAT and from external locations that must be brought in. 

- hec
  - hec.heclib
  - hec.io
  - hec.hecmath
- rma
  - com.rma
  - rma.util.RMAConst

## Usage
### Usage withing the build process
To be added later.

### Post build implementation
After making any desired changes to the code, files must be replaced in the scripts folder in a WTMP Sacramento Study folder. 
Scripts will be triggered as various points of the modeling workflow.

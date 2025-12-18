import sys
print(sys.path)
from hec.heclib.dss import HecDss

import DSS_Tools
reload(DSS_Tools)
import flowweightaverage
reload(flowweightaverage)

# TESTING
from java.io import File
from hec.heclib.dss import HecDss
from hec.io import TimeSeriesContainer
from hec.heclib.util import HecTime
from java.util import Vector
from java.text import SimpleDateFormat
from java.util import Date
from jarray import array
from java.lang import Double
from java.io import BufferedReader, FileReader

##
#
# computeAlternative function is called when the ScriptingAlternative is computed.
# Arguments:
#   currentAlternative - the ScriptingAlternative. hec2.wat.plugin.java.impl.scripting.model.ScriptPluginAlt
#   computeOptions     - the compute options.  hec.wat.model.ComputeOptions
#
# return True if the script was successful, False if not.
# no explicit return will be treated as a successful return
#
##

input_names = [
'WQ_BOUND_TEMP_SRC_CCTUNNEL',
'LEWISTON RESERVOIR-TUNNEL',
'WQ_BOUND_TEMP_SRC_SCTUNNEL',
'WHISKEYTOWN LAKE-TUNNEL'
]

output_names = [
'LewistonTunnel_With_Heating',
'WhiskeyTownTunnel_With_Heating'
]

copy_to_resim_names = [
'TargetTemp',
]

ResSim_linked_rec = 'LEWISTON RESERVOIR-TUNNEL' # only used for get f part, could be any W2 output

def fixFpartToInput(locations_paths, outpath):
    # get F-part from input locations
    location_fpart = locations_paths[0].split('/')[6]
    out_parts = outpath.split('/')
    out_parts[6] = location_fpart
    return '/'.join(out_parts)

def getOutputPaths(locations_paths,currentAlternative):
    outputlocations_obs = currentAlternative.getOutputDataLocations()    
    outputPaths = []
    for opl in outputlocations_obs:        
        path = str(currentAlternative.createOutputTimeSeries(opl))
        outputPaths.append(fixFpartToInput(locations_paths,path))
    return outputPaths
    
def computeAlternative(currentAlternative, computeOptions):
    currentAlternative.addComputeMessage("Computing ScriptingAlternative:" + currentAlternative.getName() )
    locations_obj = currentAlternative.getInputDataLocations()
    locations_paths = DSS_Tools.organizeLocations(currentAlternative, locations_obj, input_names, return_dss_paths=True)

    ResSim_FPart = DSS_Tools.organizeLocations(currentAlternative, locations_obj, [ResSim_linked_rec], return_dss_paths=True)[0].split('/')[6]
    print('ResSim FPart:',ResSim_FPart)

    # this is choking on linked DSS data - can't get linking to pass in the right thing.  Hardcoding ressim_copy_paths
    #ressim_copy_paths = DSS_Tools.organizeLocations(currentAlternative, locations_obj, copy_to_resim_names, return_dss_paths=True)
    ressim_copy_paths = ["/USBR/SHASTA/TEMP-WATER-TARGET//1Hour/SACTRN_BC_SCRIPT/",]
    
    output_paths = getOutputPaths(locations_paths,currentAlternative)
    
    currentAlternative.addComputeMessage('\n')
    dss_file = computeOptions.getDssFilename()
    rtw = computeOptions.getRunTimeWindow()
    starttime_str = rtw.getStartTimeString()
    endtime_str = rtw.getEndTimeString()

    # Clear Creek Tunnel Heating
    print('Clear Creek',locations_paths[:2],output_paths[0])
    DSS_Tools.add_DSS_Data(currentAlternative, dss_file, rtw, locations_paths[:2], output_paths[0])

    # Spring creek tunnel heating
    print('Spring Creek',locations_paths[2:],output_paths[1])
    DSS_Tools.add_DSS_Data(currentAlternative, dss_file, rtw, locations_paths[2:], output_paths[1])
  
    # copy other recs to ResSim fpart - there is potential for overwriting here, but seems
    # unlikely that two models would have the same rec name that you want to copy?
    for loc_path in ressim_copy_paths:
        DSS_Tools.copy_dss_ts(loc_path,new_fpart=ResSim_FPart,dss_file_path=dss_file)


    # TESTING STARTS HERE #

# Configure
csv_file_path = "D:/rollout_package_update/WTMP_SacramentoTrinity_10012025/WTMP_SacramentoTrinity_10012025/runs/W2_-_ResSim-October2024Fcst/CeQual-W2/W2/W2 Shasta 5/spr1.csv"
dss_file_path = "D:/rollout_package_update/WTMP_SacramentoTrinity_10012025/WTMP_SacramentoTrinity_10012025/runs/ResSim-October2024Fcst/ResSim-October2024Fcst.dss"
start_year = 2024

# Julian day to date
def julian_to_date(julian_day, year):
    sdf = SimpleDateFormat("yyyy-MM-dd HH:mm")
    base_date = sdf.parse("%d-01-01 00:00" % year)
    millis_per_day = 86400000
    return Date(base_date.getTime() + (int(float(julian_day)) - 1) * millis_per_day)

# Prep to read csv
reader = BufferedReader(FileReader(csv_file_path))
data_by_constituent = {}

line = reader.readLine()  # Skip header
print("First line:", line)

# Read csv and organize by constituent
while True:
    line = reader.readLine()
    if line is None:
        break
    parts = [p.strip().replace('"', '') for p in line.split(',')]
    if len(parts) < 7:
        continue

    constituent = parts[0]
    julian_day = parts[1]
    value = parts[4]  # Seg value

    if constituent not in data_by_constituent:
        data_by_constituent[constituent] = {"dates": [], "values": []}

    try:
        date = julian_to_date(julian_day, start_year)
        val = float(value)
        data_by_constituent[constituent]["dates"].append(date)
        data_by_constituent[constituent]["values"].append(val)
    except Exception as e:
        print(e)

reader.close()

# Write to dss
dss = HecDss.open(dss_file_path)

for constituent, data in data_by_constituent.items():
    if len(data["values"]) == 0 or len(data["dates"]) == 0:
        print("Skipping", constituent, "due to empty data")
        continue

    # Clean C-part
    constituent_clean = constituent.upper().replace("(", "").replace(")", "").replace("/", "-").replace(" ", "_").replace("__________________", "")
    constituent_clean = constituent_clean[:64]

    # Format D-part as date range
    date_format = SimpleDateFormat("ddMMMyyyy")
    start_date_str = date_format.format(data["dates"][0])
    end_date_str = date_format.format(data["dates"][-1])
    d_part = "%s-%s" % (start_date_str, end_date_str)

    # Convert date to DSS time integers
    sdf = SimpleDateFormat("ddMMMyyyy HH:mm")
    time_ints = []

    # E-part and interval
    e_part = "IR-DAY"
    interval = -1

    # Build pathname
    pathname = "/USBR/SHASTA/%s/%s/%s/SPR1/" % (constituent_clean, d_part, e_part)

    # Convert dates to dss time integers
    time_ints = []
    for date in data["dates"]:
        ht = HecTime()
        ht.set(sdf.format(date))  # pass a string
        time_ints.append(ht.value())

    # Create TimeSeriesContainer
    tsc = TimeSeriesContainer()
    tsc.fullName = pathname
    tsc.values = array(data["values"], 'd')
    tsc.times = array(time_ints, 'i')
    tsc.numberValues = len(data["values"])
    tsc.interval = interval
    tsc.units = "CMS" if "flow" in constituent.lower() else "M/S" if "velocity" in constituent.lower() else "DEG C"
    tsc.type = "INST-VAL"

    dss.put(tsc)
    print("Wrote:", pathname)

dss.done()







            


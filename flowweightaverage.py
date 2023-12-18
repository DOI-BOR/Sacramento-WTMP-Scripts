from hec.heclib.dss import HecDss
from hec.heclib.util.Heclib import UNDEFINED_DOUBLE
from hec.io import DSSIdentifier
from hec.io import TimeSeriesContainer
from rma.util.RMAConst import MISSING_DOUBLE
import math,sys
import DSS_Tools
reload(DSS_Tools)

#version 2.1
#modified 11-30-2022 by Scott Burdick-Yahya

def organizeLocations(currentAlternative, locations):
    locations_list = []
    if len(locations) % 2 != 0:
        currentAlternative.addComputeMessage("Uneven amount of Flow/Temp pairings. Check inputs.")
        sys.exit(1)
    for li, location in enumerate(locations):
        tspath =str(currentAlternative.loadTimeSeries(location))
        tspath = DSS_Tools.fixInputLocationFpart(currentAlternative, tspath)
        if li % 2 == 0: 
            current_pair = [tspath]
        else:
            current_pair.append(tspath)
            locations_list.append(current_pair)
    return locations_list


def flow_in_cfs(units,flows):
    if units.lower()=='cfs':
        return flows
    elif units.lower()=='cms':
        values_converted = []
        for f in flows:
            values_converted.append(f * 35.314666213)
        return values_converted
    else:
        print('FWA2: flow units not known:',units)
        sys.exit(-1)

def temperature_in_C(units,temps):
    if units.lower()=='c' or units.lower()=='deg c':
        return temps
    elif units.lower()=='f' or units.lower()=='deg f':
        values_converted = []
        for t in temps:
            values_converted.append((t - 32.0)*5.0/9.0)
        return values_converted
    else:
        print('FWA2: temperature units not known:',units)
        sys.exit(-1)

def FWA2(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None, bad_data_fill_tempC=None, last_override=False):
    '''Made a new flow-weighted average temperature function; other one was producing weirdness '''
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)

    flow_total = []
    flowtemp_total = []
    n_pairs = []

    flow_limit = 0.0 if cfs_limit is None else cfs_limit
    fill_value = UNDEFINED_DOUBLE if bad_data_fill_tempC is None else bad_data_fill_tempC
    
    for dspi, dsspaths in enumerate(DSSPaths_list):
        flow_dss_path = dsspaths[0]
        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(flow_dss_path))
        print('FWA2 Reading:',flow_dss_path)
        tsc_flow = dssFm.read(flow_dss_path, starttime_str, endtime_str, False).getData()
        flows = flow_in_cfs(tsc_flow.units,tsc_flow.values)
        print('FWA2 Reading:',temp_dss_path)
        last_rec_valid = False  # test to see if we can use override values
        try:
	        tsc_temp = dssFm.read(temp_dss_path, starttime_str, endtime_str, False).getData()
	        temps = temperature_in_C(tsc_temp.units,tsc_temp.values)
	        print('tscf',tsc_flow.values[0])
	        print(flows[0])
	        print('tsct',tsc_temp.values[0])
	        print(temps[0])
	
	        # use type of 1st temp record
	        if dspi==0:
	            nrecs = len(flows)
	            temp_type = tsc_temp.type
	
	        if len(flows) != nrecs or len(temps) != nrecs:
	            currentAlt.addComputeMessage("FWA2: record lengths do not match!")
	            print("FWA2: record lengths do not match!",nrecs,len(flows),len(temps))
	            sys.exit(-1)
	
	        for i in range(nrecs):
	            if dspi==0:
	                n_pairs.append(0) # init counter for number of flow/temp pairs in weighted average
	                flow_total.append(0.0)
	                flowtemp_total.append(0.0)
	            # perform a lot of checks on data
	            #print(i,flows[i],temps[i])
	            if not math.isnan(flows[i]) and not math.isnan(temps[i]):
	                if flows[i] > flow_limit and flows[i] < 9.0e6: # could lower upper limit to something relevant to watershed
	                    if temps[i] >= 0.0 and temps[i] <= 80.0:
	                        # passed the data checks
	                        
	                        n_pairs[i] += 1
	                        flow_total[i] += flows[i]
	                        flowtemp_total[i] += flows[i]*temps[i]
	
	                        #print(dspi,i,n_pairs[i],flows[i],temps[i],flow_total[i],flowtemp_total[i])
	        last_rec_valid = True
        except:
	        currentAlt.addComputeMessage('FWA2: data not addeded for record: '+temp_dss_path)
	        last_rec_valid = False
		
    fwat = []
    print('nrecs:',nrecs)
    for i in range(nrecs):
        if n_pairs[i] > 0:
            fwat.append(flowtemp_total[i]/flow_total[i])		
        else:
            fwat.append(fill_value)
        if last_override and last_rec_valid:
            if flows[i] > flow_limit and flows[i] < 9.0e6: # could lower upper limit to something relevant to watershed
                if temps[i] >= 0.0 and temps[i] <= 80.0:
                    fwat[i] = temps[i]
        #print(i,fwat[i])

    # use last temp container to write
    tsc_temp.type = temp_type
    tsc_temp.fullName = outputname
    tsc_temp.values = fwat
    dssFm.write(tsc_temp)
    dssFm.close()
    return 0


def FWA(currentAlt, dssFile, timewindow, DSSPaths_list, outputname, cfs_limit=None):
    starttime_str = timewindow.getStartTimeString()
    endtime_str = timewindow.getEndTimeString()
    currentAlt.addComputeMessage('Looking from {0} to {1}'.format(starttime_str, endtime_str))
    dssFm = HecDss.open(dssFile)
    dss_data = {}
    for dspi, dsspaths in enumerate(DSSPaths_list):
        flow_dss_path = dsspaths[0]
        currentAlt.addComputeMessage(str(flow_dss_path))
        flowTS = dssFm.read(flow_dss_path, starttime_str, endtime_str, False)

#        readabledates = []
#        for i in range(len(flowTS.getData().times)):
#            hecdate = flowTS.getData().getHecTime(i).dateAndTime(4) #delete later
#            readabledates.append(hecdate)

        
        flowTS = flowTS.getData()
        hecstarttimes = flowTS.times
        
        flow_units = flowTS.units
        if flow_units.lower() == 'cms':
            currentAlt.addComputeMessage('Converting cms to cfs')
            flowvals = []
            for flow in flowTS.values:
                flowvals.append(flow * 35.314666213)
            dss_data[dspi] = {'flow': flowvals} #start the dict
        else:
            dss_data[dspi] = {'flow': flowTS.values} #start the dict

        
        if cfs_limit != None:
            for fi, flow in enumerate(dss_data[dspi]['flow']):
                if flow < cfs_limit:
#                    print('Flow of {0} removed for being under limit: {1}'.format(flow, cfs_limit)) #noisy..
                    dss_data[dspi]['flow'][fi] = 0
                

        temp_dss_path = dsspaths[1]
        currentAlt.addComputeMessage(str(temp_dss_path))
        TempTS = dssFm.read(temp_dss_path, starttime_str, endtime_str, False)
        TempTS = TempTS.getData()
        tempunits = TempTS.units
        temptype = TempTS.type
        dss_data[dspi]['temp'] = TempTS.values

#    print(readabledates)
        
    for dspi in dss_data.keys():
        flowtemps = []
        offset = 0
        
        for i, flow in enumerate(dss_data[dspi]['flow']):
            
            temp = dss_data[dspi]['temp'][i]
            flowtemps.append(flow*temp)
            
        dss_data[dspi]['flowtemp'] = flowtemps
        
    total_flows = []
    dspi = dss_data.keys()[0]
    for i, flow in enumerate(dss_data[dspi]['flow']):
        temptotalflow = flow
        for key in dss_data.keys():
            if key != dspi:
                temptotalflow += dss_data[key]['flow'][i]
        total_flows.append(temptotalflow)
        
    total_flowtemp = []
    dspi = dss_data.keys()[0]
    for i, flowtemp in enumerate(dss_data[dspi]['flowtemp']):
        temptotalflowtemp = flowtemp
        for key in dss_data.keys():
            if key != dspi:             
                if not math.isnan(dss_data[key]['flowtemp'][i]):
                    if math.isnan(temptotalflowtemp):
                        temptotalflowtemp = dss_data[key]['flowtemp'][i]
                    else:
                        temptotalflowtemp += dss_data[key]['flowtemp'][i]
        total_flowtemp.append(temptotalflowtemp)
    
    FW_Avg_vals = []
    for i, flow in enumerate(total_flows):
        flowtemp = total_flowtemp[i]
        if flow == 0:
            FW_Avg_vals.append(UNDEFINED_DOUBLE)
        else:
            FW_Avg_vals.append(flowtemp / flow)
    
    tsc = TimeSeriesContainer()
    tsc.times = hecstarttimes
    tsc.fullName = outputname
    tsc.values = FW_Avg_vals
    tsc.startTime = hecstarttimes[0]
    tsc.units = tempunits
    tsc.type = temptype
    tsc.endTime = hecstarttimes[-1]
    tsc.numberValues = len(FW_Avg_vals)
    tsc.startHecTime = timewindow.getStartTime()
    tsc.endHecTime = timewindow.getEndTime()
    dssFm.write(tsc)
    dssFm.close()
    currentAlt.addComputeMessage("Number of Written values: {0}".format(len(FW_Avg_vals)))
    return 0
   
           

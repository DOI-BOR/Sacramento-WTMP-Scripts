
outlet_name_order = ['TCDU1','TCDU2','TCDU3','TCDM1','TCDM2','TCDM3','TCDL1','TCDL2','TCDL3','TCDS1','TCDS2','TCDS3']

# W2 TCD point sink elevations - hard code those here
sink_elevs = {'TCDU':{1:1042.0,
                      2:1021.0,
                      3:1000.0}
              'TCDM':{1:942.0,
                      2:921.0,
                      3:900.0}
              'TCDL':{1:830.0,
                      2:816.0,
                      3:802.0}
              'TCDS':{1:800.0,
                      2:760.0,
                      3:720.0}
               }

# load TCD change summary from DSS, make regular hourly if needed and fill


# create summed total number gates at upper, mid, lower, and side gates


# get WSE data, convert to ft if neccessary



# get generation flow data, also need hourly, only make calcs where we
# have both gate and flow data and WSE data




# init flow at each of 4 levels to 0

# 1) assign all flow to highest gate level open

# 2) divide flow at each level by 3 and assign as the point source
# sinks for each level. At this point, all flow is still on the highest
# open level

# 3) - If WSE-3 (ft) is less than the withdraw sink, record the flow.
#     - Then test if that flow is greater than zero
#	 - Find all the withdraw points with flow again above WSE, and drop them from
#	 the list of withdraw points with flow.
#	 - If there are still withdraw points (i.e. the whole gate is not out of the
#	   water, assing all flow to remaining 1 or 2 withdraw points
#	 - Else
#	   go to next level, and figure out how many withdraw points in the NEXT
#	   active level are below water, using a funny test of the 1st member of the 
#	   NEXT NEXT level (4th index of first_next_gate, which should be first_next_sink)
#      *** should the index to first_next_gate actually be 3?



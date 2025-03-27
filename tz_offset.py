# version 1.0
# Mar 2025 by Ben Saenz
import os, sys
from com.rma.model import Project
from hec.heclib.util import HecTime

# It seems there is global timezone set in WAT that is pacific time -8 hours
# That means that in times/DSS records read by WAT in scripting, the times
# come out as shifted -8 hours, even though there is not tz set in the DSS data, 
# and the final results DSS files seem to be written out without tz with the correct
# time.  This was all very confusing, and something is wrong. Scripts that rely on the
# hour of day, or interface with data not read by the model alternative, may need to
# use this offset. 
#
# for WTMP_American is was needed multiple places, so it is imported from here

hours = 8.0
days = hours/24.0

"""
tz_offset
=========

Shared timezone offset constants used by scripts that need to correct for a
timezone discrepancy observed in HEC-WAT scripting.

Background
----------
There appears to be a global timezone setting in WAT that is Pacific Time
(-8 hours). This means that in times/DSS records read by WAT in scripting,
the times come out as shifted -8 hours, even though there is no timezone set
in the DSS data itself, and the final results DSS files seem to be written
out without a timezone, with the correct time. This behavior was confusing
to track down, and something about it seems incorrect at a deeper level.

Scripts that rely on the hour of day, or that interface with data not read
through the model alternative (e.g., data read directly via `HecDss` rather
than through WAT's own time-series loading), may need to apply this offset
to correct for the discrepancy.

For WTMP_American, this offset was needed in multiple places, so it is
centralized here and imported from other scripts (e.g., `equilibrium_temp`,
`Pre-Process_ResSim_SacTrn`) rather than being hardcoded repeatedly.

Attributes
----------
hours : float
    The timezone offset, in hours, to apply as a correction (8.0).
days : float
    The same offset, expressed in fractional days (`hours / 24.0`), for use
    in day-of-year style calculations (e.g., `equilibrium_temp.get_decimal_day_of_year`).
"""

import os, sys                        # Standard library: retained (not directly used in this module)
from com.rma.model import Project     # WAT project API (retained; not directly used in this module)
from hec.heclib.util import HecTime   # HEC time utility (retained; not directly used in this module)

# It seems there is global timezone set in WAT that is pacific time -8 hours
# That means that in times/DSS records read by WAT in scripting, the times
# come out as shifted -8 hours, even though there is not tz set in the DSS data, 
# and the final results DSS files seem to be written out without tz with the correct
# time.  This was all very confusing, and something is wrong. Scripts that rely on the
# hour of day, or interface with data not read by the model alternative, may need to
# use this offset. 
#
# for WTMP_American is was needed multiple places, so it is imported from here

hours = 8.0        # Timezone correction offset, in hours
days = hours/24.0  # Same offset expressed in fractional days, for day-of-year calculations
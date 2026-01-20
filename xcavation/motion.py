
#-----------------------------------------------------------------------#
# xcavation.motion v0.3.1
# By Hunter Brooks, at UToledo, Toledo: Jan. 20, 2026
#
# Purpose: Propogate celestial coordinates using proper motion
#-----------------------------------------------------------------------#



# Import Data Management
# ------------------------------------------------------ #
import numpy as np
from astropy.time import Time
from datetime import datetime
# ------------------------------------------------------ #




# Calculate Decimal Year
#-----------------------------------------------------------------------#
def decimal_year(t):
    year_start = datetime(t.year, 1, 1)
    next_year_start = datetime(t.year + 1, 1, 1)

    year_length = (next_year_start - year_start).total_seconds()
    seconds_into_year = (t - year_start).total_seconds()

    return t.year + seconds_into_year / year_length
#-----------------------------------------------------------------------#



# Calculate Time Passed since CW Observation
#-----------------------------------------------------------------------#
def time_mjd(mjd):
    t = Time(mjd, format='mjd')
    decimal_year_cw = t.decimalyear

    now = datetime.now()
    time_passed = decimal_year(now) - decimal_year_cw
    return time_passed
#-----------------------------------------------------------------------#



# Adjust RA and Dec for Proper Motion
#-----------------------------------------------------------------------#
def proper_motion(ra, dec, pmra, pmdec, time_passed):
    dRA  = (pmra  * time_passed) / (3600 * np.cos(np.deg2rad(dec)))
    dDec = (pmdec * time_passed) / 3600

    ra_deg  = ra  + dRA
    dec_deg = dec + dDec
    return ra_deg, dec_deg
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
# xcavation.motion v1.0.0
# By Hunter Brooks, at UToledo, Toledo: Marc. 14, 2026
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
    """
    Convert a calendar date to a decimal year.

    Parameters
    ----------
      t: Calendar date and time to be converted.

    Returns
    -------
      decimal_year:  Date expressed as a decimal year (float).
    """



    # ----- Get Years ----- #
    # The Years Closest to Inputted Years
    year_start = datetime(t.year, 1, 1) # The Start of the Inputted year
    next_year_start = datetime(t.year + 1, 1, 1) # The Start of the Next Year
    # --------------------- #



    # ----- Get Difference ----- #
    # The Length of Years From Above
    year_length = (next_year_start - year_start).total_seconds()
    seconds_into_year = (t - year_start).total_seconds()
    # -------------------------- #



    return t.year + seconds_into_year / year_length
#-----------------------------------------------------------------------#



# Calculate Time Passed Since Inputted MJD
#-----------------------------------------------------------------------#
def time_mjd(mjd, mjd_obs):
    """
    Compute the elapsed time in years between two epochs given in
    Modified Julian Date (MJD).

    Parameters
    ----------
      mjd: Reference MJD corresponding to the inputted epoch (float)
      mjd_obs: Observed MJD corresponding to the observed epoch (float)

    Returns
    -------
      time_passed: Time elapsed from `mjd` to `mjd_obs` in decimal years (float)
    """



    # ----- MJD Input ----- #
    # Convert MJD (input) to decimal year (i.e. 2026.2315)
    t = Time(mjd, format='mjd')
    decimal_year_inputted = t.decimalyear
    # --------------------- #



    # ----- MJD Observed ----- #
    # Convert MJD (observed) to decimal year (i.e. 2026.2315)
    t = Time(mjd_obs, format='mjd')
    decimal_year_obs = t.decimalyear
    # ------------------------ #



    # Get the Time of Final Time Propagated
    time_passed = decimal_year_obs - decimal_year_inputted # Calc. Time Passed
    return time_passed
#-----------------------------------------------------------------------#



# Adjust RA and Dec for Proper Motion
#-----------------------------------------------------------------------#
def proper_motion(ra, dec, pmra, pmdec, time_passed):
    """
    Propagate right ascension and declination using proper motion.

    Parameters
    ----------
      ra: Initial R.A. in degrees (float)
      dec : Initial Decl. in degrees (float)
      pmra: Proper motion in R.A. in arcseconds per year (float)
      pmdec: Proper motion in Decl. in arcseconds per year (float)
      time_passed: Time interval over which to propagate the position (float)

    Returns
    -------
      ra_deg: Propagated R.A. in degrees (float)
      dec_deg: Propagated Decl. in degrees (float)
    """



    # ----- Delta R.A. and Decl. ----- #
    dRA  = (pmra  * time_passed) / (3600 * np.cos(np.deg2rad(dec))) # Delta R.A.
    dDec = (pmdec * time_passed) / 3600 # Delta Decl.
    # -------------------------------- #



    # ----- Final R.A. and Decl. ----- #
    ra_deg  = ra  + dRA # R.A.
    dec_deg = dec + dDec # Decl.
    # -------------------------------- #



    return ra_deg, dec_deg
#-----------------------------------------------------------------------#


#-----------------------------------------------------------------------#
# xcavation.genspec v1.0.1
# By Hunter Brooks, at UToledo, Toledo: Apr. 08, 2026
#
# Purpose: Main API Function for SphereX Data Retrieval and Photometry
#-----------------------------------------------------------------------#



# Import Internal Modules
# ------------------------------------------------------ #
from .aperture import spherex_aperature_phot
from .motion import decimal_year, time_mjd, proper_motion
from .quality import finder_chart, spectra_plot
# ------------------------------------------------------ #


# Import Data Management
# ------------------------------------------------------ #
import time
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
from datetime import datetime
from astropy.table import Table
from astroquery.ipac.irsa import Irsa
from dataclasses import dataclass, field
# ------------------------------------------------------ #



# Import PYVO Query Features
# ------------------------------------------------------ #
import pyvo
import http.client
import urllib.error
from astropy.utils.data import conf
# ------------------------------------------------------ #



# Import WCS, Photometry, and Plotting
# ------------------------------------------------------ #
from astropy.coordinates import SkyCoord
# ------------------------------------------------------ #



# Import Multithreading and Watching Packages
# ------------------------------------------------------ #
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
# ------------------------------------------------------ #



# User Configuration Data
# ------------------------------------------------------ #
def variable_verify(ra, dec, # Core Inputs
                    r_fwhm, r_annulus_in, r_annulus_out, # Photometry Radii
                    pmra, pmdec, mjd_epoch, mjd_query, # Proper Motion Propagation
                    save_data, output_path, # Saved Data
                    threads, # Number of Multi-Threads
                    enable_print, ram_download, retry_count, # User Control
                    clean_type, bad_bits, # Cleaning Flux Data
                    background_type, # Background Type
                    cutout_size, # Cutout Size
                    zodi_subtract, # ZODI Light Subtraction
                    sigclip_sigma, sigclip_maxiters): # Astropy.SigmaClip
  """
  Verify user input parameters for SPHEREx aperture photometry.

  Parameters
  ----------
    ra: R.A. of the target in degrees (0 <= ra <= 360) (float)
    dec: Decl. of the target in degrees (-90 <= dec <= 90) (float)
    r_fwhm: Aperture radius in units of PSF FWHM (float)
    r_annulus_in: Inner radius of background annulus (float)
    r_annulus_out: Outer radius of background annulus (float)
    pmra: Proper motion in R.A. (float)
    pmdec: Proper motion in Decl. (float)
    mjd: MJD corresponding to the input coordinates (float)
    save_data: Whether to save the photometry output (bool)
    output_path: Path for saving the output file (str)
    threads: Number of CPU threads to allocate (int)
    enable_print: Enable or disable console print statements (bool)
    ram_download: If True, download FITS entirely into memory (bool)
    retry_count: Number of times to retry downloading a file (int)
    clean_type: Method for cleaning flux data (str)
    bad_bits: List of bit indices in FLAG extension (list)
    cutout_size: Cutout Size in arcsec
    zodi_subtract: Whether to subtract ZODI light
    sigclip_sigma: Astropy.SigmaClip(sigma=5.0)
    sigclip_maxiters: Astropy.SigmaClip(maxiters=5)

  Returns
  -------
    bool
        True if all input variables are valid, False otherwise
  """

  # Checks R.A. Decl. Inputs
  if (isinstance(ra, (int, float)) or isinstance(dec, (int, float))
                or (0 <= ra <= 360) or (-90 <= dec <= 90)) == False:
    print('Please Input Valid: R.A. (deg, float) or Decl. (deg, float)')
    return False

  # Checks PMRA, PMDEC, and MJD Inputs
  if (isinstance(pmra, (int, float))
      or isinstance(pmdec, (int, float))
               or isinstance(mjd_epoch, (int, float))
                  or isinstance(mjd_query, (int, float))) == False:
    print('Please Input Valid: PMRA (float), PMDEC (float), or MJD (float)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Aperture Radius Input
  if (isinstance(r_fwhm, (int, float)) or (r_fwhm) > 0) == False:
    print('Please Input Valid: Aperture Radius')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Inner Annulus Radius Input
  if (isinstance(r_annulus_in, (int, float)) or ((r_annulus_in) > 0)) == False:
    print('Please Input Valid: Inner Annulus Radius')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Outer Annulus Radius Input
  if (isinstance(r_annulus_out, (int, float)) or ((r_annulus_out) > 0)) == False:
    print('Please Input Valid: Outer Annulus Radius')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Saved Data Inputted Variable
  if type(save_data) is not bool:
    print('Please Input Valid: Save_data (True or False)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Output Path Inputted
  if (type(output_path) is not str):
    print('Please Input Valid: Spectrum Table (string)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Inputted Threads
  if (type(threads) is not int) or (1 >= threads >= 50):
    print('Please Input a Valid: Amount of Allocated CPU Threads')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Inputted Print Variables
  if (type(enable_print) is not bool):
    print('Please Input a Valid: Enable Print (True or False)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks RAM Download Variable
  if (type(ram_download) is not bool):
    print('Please Input a Valid: Ram Download (True or False)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Download Retry Count Variables
  if (type(retry_count) is not int) or (retry_count < 1):
    print('Please Input a Valid: Retry Count (integer)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Clean Type Variable
  if ((clean_type != 'none') and (clean_type != 'mask')
      and (clean_type != 'median_mask') and (clean_type != 'interp_mask')):
    print('Please Input a Valid: Clean Tyle (none, laplacian, mask)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Background Type
  if ((background_type != 'mean') and (background_type != 'median')
      and (background_type != 'mode')):
    print('Please Input a Valid: Background Type (mean, median, mode)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Flagged Variable List
  if (type(bad_bits) is not list) or (len(bad_bits) < 1):
    print('Please Input a Valid: Bad Bits (list)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Cutout Size
  if (type(cutout_size) is not int) or (cutout_size < 1):
    print('Please Input a Valid: Cutout Size (int)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks ZODI Variable
  if (type(zodi_subtract) is not bool):
    print('Please Input a Valid: ZODI Subtraction (bool)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks sigclip_sigma Variable
  if ((type(sigclip_sigma) is not float) and (type(sigclip_sigma) is not int)) or (sigclip_sigma < 0):
    print('Please Input a Valid: SigmaClip (float)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks sigclip_maxiters Variable
  if ((type(sigclip_maxiters) is not float) and (type(sigclip_maxiters) is not int)) or (sigclip_maxiters < 0):
    print('Please Input a Valid: Max Iterations (float)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Returns True if Every Variable is Good
  return True
# ------------------------------------------------------ #



# User Configuration Data
# ------------------------------------------------------ #
mjd_now = Time(datetime.now()).mjd # Current MJD
@dataclass
class genspec_profile:
    r_fwhm: float = 2 # Aperture Radius FWHM
    r_annulus_in: float = 3.5 # Inner Annulus Radius FWHM
    r_annulus_out: float = 10 # Outer Annulus Radius FWHM
    pmra: float = 0 # Proper Motion R.A. (arcsec/yr)
    pmdec: float = 0 # Proper Motion Decl. (arcsec/yr)
    mjd_epoch: float = 61000 # Modified Julian Date of Epoch
    mjd_query: float = mjd_now # Modified Julian Date of Query
    save_data: bool = False # Save Q.A., Table, and Spectral Data
    output_path: str = 'xcavation_save' # Saved Location
    threads: int = 1 # Number of Threads for Multi-Threading
    enable_print: bool = True # Whether to Print
    ram_download: bool = False # Download Images to RAM
    retry_count: int = 10 # How Many HTML Retries
    clean_type: str = 'none' # What type of Image Cleaning
    bad_bits: list = field(default_factory=lambda: [0,1,10,11]) # Flags Removed
    background_type: str = 'mean' # What type of background subtraction
    cutout_size: int = 150 # Cutout Size in arcsec
    zodi_subtract: bool = True # Whether ZODI Light is subtracted
    sigclip_sigma: float = 5 # Astropy.SigmaClip(sigma=5.0)
    sigclip_maxiters: float = 5 # Astropy.SigmaClip(maxiters=5)
# ------------------------------------------------------ #



# Creates Multi-Try Query
# ------------------------------------------------------ #
def retry(func, *args, retries=10, delay=1, **kwargs):
    """
    Execute a function with automatic retries on failure.

    Parameters
    ----------
      func: Called function
      *args: All relavent data for called function
      retries : Max number of retries (int)
      delay: Delay between retries (float)

    Returns
    -------
      result: object or None
          Return value of `func` if successful; otherwise None if all
          retries fail.
    """



    # ----- Repeat Calling Loop ----- #
    for attempt in range(1, retries + 1):
        try:  # Call main function (aperture photometry)
            return func(*args, **kwargs)

        except Exception as e:
            err_msg = str(e)

            # Special handling for 502 Bad Gateway
            if "502 Server Error" in err_msg and "Bad Gateway" in err_msg:
                print(f"[Attempt {attempt}/{retries}] Query failed with 502 Bad Gateway, Retrying..... ")
            else:
                print(f"[Attempt {attempt}/{retries}] Error occurred: {err_msg}")

            if attempt < retries:  # Keep going if retries remain
                time.sleep(delay)
            else:  # Fail if max retries reached
                print("Maximum retries reached. Returning None.")
                return None
    # ------------------------------- #
# ------------------------------------------------------ #



# Multi-Thread the SphereX Query Function
# ------------------------------------------------------ #
def genspec(ra, dec, config: genspec_profile):
    """
    Generate a SPHEREx spectrophotometric spectrum at a given sky position.
    This is the main function for the xcavation SPHEREx
    spectrophotometry tool. It performs multi-threaded SPHEREx aperture
    photometry calculations while optionally produces finder charts, spectral
    plots, and spectral table.

    Parameters
    ----------
      ra: R.A. of the target in degrees (float).
      dec: Decl. of the target in degrees (float).
      config: User generated configuration for queried data

    Returns
    -------
      output : pandas.DataFrame
          - wavelength: float
              Wavelength at queried location
          - delta_lambda: float
              Spectral resolution using R value solution
          - flux: float
              Aperture flux in μJy
          - flux_err: float
              Propagated flux uncertainty in μJy
          - flag_count: str
              Dictionary of flagged pixel counts per bit
          - SNR: float
              Signal-to-noise ratio of given fluxes
          - flux_cutout: array
              Flux cutout
          - flag_cutout: array
              FLAG cutout
          - aperture: array
              Boolean aperture cutout
          - annulus: array
              Boolean annulus cutout
          - x_loc, y_loc: float
              Pixel coordinates of the target in the cutout
          - ap_radius: float
              Aperture radius
          - inner_annulus, outer_annulus : float
              Annulus radii
          - url: str
              Source FITS URL.
          - mjd: float
              Observation MJD
          - cutout_size
              Cutout Size in Arcsec
          - zodi_subtract
              Whether to subtract ZODI light
          - sigclip_sigma
              Astropy.SigmaClip(sigma=5.0)
          - sigclip_maxiters
              Astropy.SigmaClip(maxiters=5)

    Raises
    ------
      ValueError
          If input variables fail validation or network failure
    """



    # ----- Obtain User Config. Data ----- #
    r_fwhm = config.r_fwhm # Radius in FWHM
    r_annulus_in = config.r_annulus_in # Inner Annulus Radius in FWHM
    r_annulus_out = config.r_annulus_out # Outer Annulus Radius in FWHM
    pmra = config.pmra # Proper Motion of R.A. in arcsec/yr
    pmdec = config.pmdec # Proper Motion of Decl. in arcsec/yr
    mjd = config.mjd_epoch # R.A. and Decl. Start MJD
    mjd_query = config.mjd_query # User Defined Query MJD
    save_data = config.save_data # Boolean on Whether to Save Data
    output_path = config.output_path # Path to Save Data
    threads = config.threads # Number of Multi-Threads
    enable_print = config.enable_print # Boolean on Whether to Print Results
    ram_download = config.ram_download # Download FITS Images to RAM
    bad_bits = config.bad_bits # Which Flags to Mask
    retry_count = config.retry_count # Number of Retries
    clean_type = config.clean_type # Which Method of Bad Pixel Fixes to USE
    background_type = config.background_type # Which type of background subtraction
    cutout_size = config.cutout_size # Cutout Size in Arcsec
    zodi_subtract = config.zodi_subtract # Whether ZODI Light is subtracted
    sigclip_sigma = config.sigclip_sigma # Astropy.SigmaClip(sigma=5.0)
    sigclip_maxiters = config.sigclip_maxiters # Astropy.SigmaClip(maxiters=5)
    # ------------------------------------ #



    # ----- Verify User Inputted Variables ----- #
    verify = variable_verify(ra, dec, # Positional Data
                             r_fwhm, r_annulus_in, r_annulus_out, # Radii Data
                             pmra, pmdec, mjd, mjd_query, # Proper Motion Data
                             save_data, output_path, # Saved Data
                             threads, enable_print, ram_download, retry_count, # Misc.
                             clean_type, bad_bits, # Masking Data
                             background_type, # Background Type
                             cutout_size, # Cutout Size
                             zodi_subtract, # ZODI Light Subtraction
                             sigclip_sigma, sigclip_maxiters) # Astropy.SigmaClip

    # Raise Error of User Inputted Variables is Bad
    if verify is False:
      raise ValueError('Poorly Assigned Variable. Please Read Documentation.')
    # ------------------------------------------ #



    # ----- Proper Motion Propagation ----- #
    # Calculate Propagated Proper Motion
    mjd_now = mjd_query
    time_passed = time_mjd(mjd, mjd_now) # Years since observation
    ra_deg, dec_deg = proper_motion(ra, dec,
                          pmra, pmdec, time_passed) # Adjust for Proper Motion

    # Print Current R.A. and Decl. Data if Enable Print is True
    if enable_print is True:
      print('\n+====================================================================+')
      print('       Starting xcavation.genspec SPHEREx Spectrophotometry Tool')
      print('+====================================================================+')
      print(f'''Queried Coordinates: R.A. = {round(ra_deg, 6)} (deg), Decl. = {round(dec_deg, 6)} (deg)''') # Print Adjusted Coordinates
      print(f"Query Time: {datetime.now()}")
      print('+====================================================================+\n')
    # ------------------------------------- #



    # ----- SphereX Query ----- #
    # Query For All SPHEREx URLs in a Nearby Search
    ra_deg, dec_deg = ra_deg * u.degree, dec_deg * u.degree # Convert to astropy Units
    size = cutout_size/3600 * u.degree
    service = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")

    # Query SphereX QR2 (Will Need Updates)
    query = f"""
              SELECT
                  'https://irsa.ipac.caltech.edu/' || a.uri || '?center={ra_deg.value},{dec_deg.value}d&size={size.value}' AS uri,
                  p.time_bounds_lower
              FROM spherex.artifact a
              JOIN spherex.plane p ON a.planeid = p.planeid
              WHERE 1 = CONTAINS(POINT('ICRS', {ra_deg.value}, {dec_deg.value}), p.poly)
            """

    results = service.search(query)
    urls = results['uri'].tolist() # Gets the Image URLs

    # Read in the Number of Each Detectors
    D1, D2, D3, D4, D5, D6 = 0, 0, 0, 0, 0, 0 # Creates Empty Numbers
    for url in urls:
      # Splits URL to Obtain Detector Number
      num_detector = int(url.split('D')[len(url.split('D')) - 1][0])
      if num_detector == 1: # Detector 1
        D1 += 1
      if num_detector == 2: # Detector 2
        D2 += 1
      if num_detector == 3: # Detector 3
        D3 += 1
      if num_detector == 4: # Detector 4
        D4 += 1
      if num_detector == 5: # Detector 5
        D5 += 1
      if num_detector == 6: # Detector 6
        D6 += 1
    # ------------------------- #



    # ----- Multi-Thread Wavelength and Aperature Calc. ----- #
    output = [] # Creates Empty List
    with ThreadPoolExecutor(max_workers=threads) as executor:

        # Multi-Thread w/ tqdm
        futures = [executor.submit(retry, spherex_aperature_phot, url,
                          [ra, dec], [pmra, pmdec], mjd,
                            r_fwhm, r_annulus_in, r_annulus_out,
                              ram_download, clean_type, bad_bits,
                                   background_type,
                                   cutout_size,
                                   zodi_subtract,
                                   sigclip_sigma, sigclip_maxiters,
                                   retries = retry_count) for url in urls]

        # Prints Start of Query
        if enable_print is True:
          print('+====================================================================+')
          print(f'{len(urls)} Images Were Found in SPHEREx')

        # Sets Up Progress Bar
        for future in tqdm(futures, desc="Calculating", ncols=70, colour = '#C5D89D',
                           bar_format='''{l_bar}{bar}|{n_fmt}/{total_fmt} [{elapsed}<{remaining}]''',
                           ascii="░▒█", disable=not enable_print):

            # Starts Multi-Thread
            result = future.result()

            # Raises Error if Network is Down
            if result is None:
              raise ValueError('Network Connectivity Issue. Please Try Again')

            # Adds to Empty Array For Total Dictionary
            output.append({
                "wavelength": result["wavelength"], # Wavelength
                "delta_lambda": result["delta_lambda"], # Resolving Power
                "flux": result["flux"], # Aperture Flux
                "flux_err": result["flux_err"], # Aperture Flux Error
                "flag_count": str(result["flag_count"]), # Dict. of Flags
                "flag": result["flag"], # Bad Bit Count Flag
                "SNR": result["SNR"], # Signal-to-Noise Ratio
                "flux_cutout": result["flux_cutout"], # Flux Cutout
                "flag_cutout": result["flag_cutout"], # Bool. Flag Cutout
                "aperture": result["aperture"], # Aperture Mask
                "annulus": result["annulus"], # Annulus Mask
                "x_loc": result["x_loc"], # X Pixel Location in Cutout
                "y_loc": result["y_loc"], # Y Pixel Location in Cutout
                "ap_radius": result['ap_radius'], # Aperture Radius
                "inner_annulus": result['inner_annulus'], # Inner Annulus Radius
                "outer_annulus": result['outer_annulus'], # Outer Annulus Radius
                "url": result["url"], # Image URL
                "mjd": result["mjd"] # Observation MJD
            })

        # Prints Finished Query
        if enable_print is True:
          print('Spectral Calulcation Complete')
          print('+====================================================================+\n')

    # Converts Final Dictionary to pandas DataFrame
    output = pd.DataFrame(output)
    # ------------------------------------------------------- #



    # ----- Finder Chart and Spectral Plots ----- #
    # Only Queries if Selected
    if save_data is True:

      # Start Q.A. Finder Chart Creation
      finder_chart(output, output_path + '.pdf')

      # Prints Finished Finder Chart
      if enable_print is True:
        print('+====================================================================+')
        print('Finder Chart Was Successfully Saved')
        print(f'Spectral Plot of: R.A. = {round(ra, 3)} (deg) and Decl. = {round(dec, 3)} (deg)')

      # Gets Relevant Data for Table
      lambda_list = np.array(output['wavelength']) # Wavelength
      flux_list = np.array(output['flux']) # Flux
      flux_err_list = np.array(output['flux_err']) # Flux Error
      delta_lambda = np.array(output['delta_lambda']) # Resolving Power
      flag = np.array(output['flag']) # Flag Count
      mjd_list = np.array(output['mjd']) # MJD List
      flag_count = np.array(output['flag_count']) # Flag Dictionary

      # Sort Relevant Data Based on Wavelength
      sort_idx = np.argsort(lambda_list) # Sort Based on Wavelength
      lambda_list = lambda_list[sort_idx] # Sort Wavelength
      flux_list = flux_list[sort_idx] # Sort Flux
      flux_err_list = flux_err_list[sort_idx] # Sort Flux Error
      delta_lambda = delta_lambda[sort_idx] # Sort Resolving Power
      flag = flag[sort_idx] # Flag Count
      mjd_list = mjd_list[sort_idx] # Sorted MJD List
      flag_count = flag_count[sort_idx] # Sort Flag Dictionary

      # Creates Table with Header
      t = Table([lambda_list, flux_list, flux_err_list,
              delta_lambda, flag, mjd_list, flag_count],
          names=('lambda', 'flux', 'flux_err', 'delta_lambda', 'flag', 'mjd', 'flag_count'))

      # Opens Created Files To Edit
      with open(output_path + '.txt', 'w') as f:

          # Generates Relevant Table Header
          f.write(f'# Description: Generated Using xcavation.genspec SPHEREx Spectrophotometry Tool\n')
          f.write(f'# Publisher: IPAC/IRSA\n')
          f.write(f'# Contact: Mr. Hunter Brooks\n')
          f.write(f'# Email: hbrooks8@rockets.utoledo.edu\n')
          f.write(f'# Data Source: SPHEREx QR2 Images via IRSA\n')
          f.write(f'# Title: Spectrum (ASCII)\n')
          f.write(f'# RA: {ra}\n')
          f.write(f'# Dec: {dec}\n')
          f.write(f'# Wavelength Range: {round(np.nanmin(output['mjd']), 3)} to {round(np.nanmax(output['mjd']), 3)}\n')
          f.write(f'# Spectral resolution (Δλ/λ): 35 – 130\n')
          f.write(f'# Number of images measured: {len(urls)}\n')
          f.write(f'# Photometry Type: Aperture\n')
          f.write(f'# Average SNR: {round(np.nanmean(output['flux']/output['flux_err']), 2)}\n')
          f.write(f'# Average Flux Error: {round(np.nanmean(output['flux_err']))}\n')
          f.write(f'# Images per detector: D1={D1}, D2={D2}, D3={D3}, D4={D4}, D5={D5}, D6={D6}\n')
          f.write(f'# References: https://spherex.caltech.edu/page/publications; https://github.com/huntbrooks85/Xcavation\n')
          f.write(f'# Column 1: Wavelength (micron)\n')
          f.write(f'# Column 2: Flux (uJy)\n')
          f.write(f'# Column 3: Flux Error (uJy)\n')
          f.write(f'# Column 4: Delta Lambda (microns)\n')
          f.write(f'# Column 5: Number of Flags (int)\n')
          f.write(f'# Column 6: MJD of Observation (float)\n')
          f.write(f'# Column 7: Aperture Flag Bit Count (dictionary)\n')
          t.write(f, format='ascii.commented_header', overwrite=True)

      # Print Finished Saving Data
      if enable_print is True:
        print('Spectrum Table Was Successfully Saved')
        print('+====================================================================+\n')
    # ------------------------------------------- #



    # ----- Print xcavation Summary ----- #
    if enable_print is True:
      print('+====================================================================+')
      print('              xcavation Spectrophotometry Tool Summary')
      print('+====================================================================+')

      # Plots Spectral Plot
      print(f"Spectra Plot:")
      spectra_plot(output)

      # Print Coordinate Data
      print(f"Queried Coordinates: R.A. = {round(ra_deg, 6)} (deg), Decl. = {round(dec_deg, 6)} (deg)") # Print Adjusted Coordinates
      print(f"SIMBAD Query: https://simbad.u-strasbg.fr/simbad/sim-coo?Coord={ra}+{dec}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=0.5&Radius.unit=arcmin&submit=submit+query&CoordList=")
      print(f"Query Finished: {datetime.now()}\n")

      # Print Number of Image Data
      print(f"Number of Images Found: {len(urls)}")
      print(f"Images Per Detector: D1={D1}, D2={D2}, D3={D3}, D4={D4}, D5={D5}, D6={D6}\n")

      # Print Flux Data
      print(f'Average SNR: {round(np.nanmean(output['flux']/output['flux_err']), 2)}')
      print(f'Average Flux Error: {round(np.nanmean(output['flux_err']))}\n')

      # Print Spectral Data
      print(f'Observation Range: {round(np.nanmin(output['mjd']), 1)} to {round(np.nanmax(output['mjd']), 1)}')
      print(f'Wavelength Range: {round(np.nanmin(output['wavelength']), 3)} to {round(np.nanmax(output['wavelength']), 3)}')
      print(f"Flux Range: {round(np.nanmin(output['flux']), 3)} to {round(np.nanmax(output['flux']), 3)}")

      # Print Where Data Was Saved
      if save_data == True:
        print(f'\nData Saved To: {output_path}')
      print('+====================================================================+')
    # ----------------------------------- #



    # Returns The Output Dictionary
    return output
# ------------------------------------------------------ #

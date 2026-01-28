
#-----------------------------------------------------------------------#
# xcavation.genspec v0.4.0
# By Hunter Brooks, at UToledo, Toledo: Jan. 28, 2026
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
                    pmra, pmdec, mjd, # Proper Motion Propagation
                    save_data, output_path, # Saved Data
                    threads, # Number of Multi-Threads
                    enable_print, ram_download, retry_count, # User Control
                    clean_type, bad_bits, # Cleaning Flux Data
                    background_type):
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

  Returns
  -------
    bool
        True if all input variables are valid, False otherwise
  """

  # Checks R.A. Decl. Inputs
  if ((type(ra) is not float) or (type(dec) is not float)
                or (0 >= ra >= 360) or (-90 >= dec >= 90)):
    print('Please Input Valid: R.A. (deg, float) or Decl. (deg, float)')
    return False

  # Checks PMRA, PMDEC, and MJD Inputs
  if (((type(pmra) is not float) and (type(pmra) is not int))
      or ((type(pmdec) is not float) and (type(pmdec) is not int))
               or (type(mjd) is not float and type(mjd) is not int)):
    print('Please Input Valid: PMRA (float), PMDEC (float), or MJD (float)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Aperture Radius Input
  if (type(r_fwhm) is not float and type(r_fwhm) is not int) or ((r_fwhm) < 0):
    print('Please Input Valid: Aperture Radius')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Inner Annulus Radius Input
  if ((type(r_annulus_in) is not float and type(r_annulus_in) is not int)
                                                  or ((r_annulus_in) < 0)):
    print('Please Input Valid: Inner Annulus Radius')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  # Checks Outer Annulus Radius Input
  if ((type(r_annulus_out) is not float and type(r_annulus_out) is not int)
                                                    or ((r_annulus_out) < 0)):
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
  if ((clean_type != 'none') and (clean_type != 'lacosmic')
      and (clean_type != 'mask') and (clean_type != 'median_mask')
        and (clean_type != 'interp_mask')):
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

  # Returns True if Every Variable is Good
  return True
# ------------------------------------------------------ #



# User Configuration Data
# ------------------------------------------------------ #
@dataclass
class genspec_profile:
    r_fwhm: float = 2 # Aperture Radius FWHM
    r_annulus_in: float = 7.5 # Inner Annulus Radius FWHM
    r_annulus_out: float = 15 # Outer Annulus Radius FWHM
    pmra: float = 0 # Proper Motion R.A. (arcsec/yr)
    pmdec: float = 0 # Proper Motion Decl. (arcsec/yr)
    mjd: float = 61000 # Modified Julian Date
    save_data: bool = False # Save Q.A., Table, and Spectral Data
    output_path: str = 'xcavation_save' # Saved Location
    threads: int = 1 # Number of Threads for Multi-Threading
    enable_print: bool = True # Whether to Print
    ram_download: bool = False # Download Images to RAM
    retry_count: int = 10 # How Many HTML Retries
    clean_type: str = 'none' # What type of Image Cleaning
    bad_bits: list = field(default_factory=lambda: [0,1,10,11]) # Flags Removed
    background_type: str = 'mean' # What type of background subtraction
# ------------------------------------------------------ #



# Creates Multi-Try Query
# ------------------------------------------------------ #
def retry(func, *args, retries=10, delay=0.5, **kwargs):
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
        try: # Call main function (aperture photometry)
            return func(*args, **kwargs)
        except Exception as e: # If failed retry (usual fail is HTML)
            if attempt < retries: # Keep Going if more retries
                time.sleep(delay)
            else: # Fail if max retries reached
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
    mjd = config.mjd # R.A. and Decl. Start MJD
    save_data = config.save_data # Boolean on Whether to Save Data
    output_path = config.output_path # Path to Save Data
    threads = config.threads # Number of Multi-Threads
    enable_print = config.enable_print # Boolean on Whether to Print Results
    ram_download = config.ram_download # Download FITS Images to RAM
    bad_bits = config.bad_bits # Which Flags to Mask
    retry_count = config.retry_count # Number of Retries
    clean_type = config.clean_type # Which Method of Bad Pixel Fixes to USE
    background_type = config.background_type # Which type of background subtraction
    # ------------------------------------ #



    # ----- Verify User Inputted Variables ----- #
    verify = variable_verify(ra, dec, # Positional Data
                             r_fwhm, r_annulus_in, r_annulus_out, # Radii Data
                             pmra, pmdec, mjd, # Proper Motion Data
                             save_data, output_path, # Saved Data
                             threads, enable_print, ram_download, retry_count, # Misc.
                             clean_type, bad_bits, # Masking Data
                             background_type) # Background Type

    # Raise Error of User Inputted Variables is Bad
    if verify is False:
      raise ValueError('Poorly Assigned Variable. Please Read Documentation.')
    # ------------------------------------------ #



    # ----- Proper Motion Propagation ----- #
    # Calculate Propagated Proper Motion
    mjd_now = Time(datetime.now()).mjd # Current MJD
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
    coord = SkyCoord(ra_deg, dec_deg, unit='deg') # Convert to astropy Units
    results = Irsa.query_sia(pos=(coord, 10 * u.arcsec),
              collection='spherex_qr2') # Query SphereX QR2 (Will Need Updates)
    urls = results["access_url"] # Gets the Image URLs

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

      # Plots Spectral Plot
      spectra_plot(output)

      # Gets Relevant Data for Table
      lambda_list = np.array(output['wavelength']) # Wavelength
      flux_list = np.array(output['flux']) # Flux
      flux_err_list = np.array(output['flux_err']) # Flux Error
      delta_lambda = np.array(output['delta_lambda']) # Resolving Power
      flag_count = np.array(output['flag_count']) # Flag Dictionary

      # Sort Relevant Data Based on Wavelength
      sort_idx = np.argsort(lambda_list) # Sort Based on Wavelength
      lambda_list = lambda_list[sort_idx] # Sort Wavelength
      flux_list = flux_list[sort_idx] # Sort Flux
      flux_err_list = flux_err_list[sort_idx] # Sort Flux Error
      delta_lambda = delta_lambda[sort_idx] # Sort Resolving Power
      flag_count = flag_count[sort_idx] # Sort Flag Dictionary

      # Creates Table with Header
      t = Table([lambda_list, flux_list, flux_err_list,
              delta_lambda, flag_count],
          names=('lambda', 'flux', 'flux_err', 'delta_lambda', 'flag_count'))

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
          f.write(f'# Column 5: Aperture Flag Bit Count (dictionary)\n')
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


#-----------------------------------------------------------------------#
# xcavation.genspec v0.3.1
# By Hunter Brooks, at UToledo, Toledo: Jan. 20, 2026
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



# Verify Input Variables
# ------------------------------------------------------ #
def variable_verify(ra, dec, style, r_fwhm, pmra, pmdec, mjd, save_data, output_path, threads, enable_print, ram_download, bad_bits):
  if (type(ra) is not float) or (type(dec) is not float) or (0 >= ra >= 360) or (-90 >= dec >= 90):
    print('Please Input Valid: R.A. (deg, float) or Decl. (deg, float)')
    return False

  if ((type(pmra) is not float) and (type(pmra) is not int)) or ((type(pmdec) is not float) and (type(pmdec) is not int)) or (type(mjd) is not float and type(mjd) is not int):
    print('Please Input Valid: PMRA (arcsec/yr, float), PMDEC (arcsec/yr, float), or MJD (float)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  if (type(style) is not str) or (style != 'psf' and style != 'aperture'):
    print('Please Input Valid: Photometry Style ("aperture" or "psf")')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  if style == 'aperture':
    if (type(r_fwhm) is not float and type(r_fwhm) is not int) or ((r_fwhm) < 0):
      print('Please Input Valid: Aperture Radius')
      print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
      return False

  if type(save_data) is not bool:
    print('Please Input Valid: Save_data (True or False)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  if (type(output_path) is not str):
    print('Please Input Valid: Spectrum Table (string)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  if (type(threads) is not int) or (1 >= threads >= 50):
    print('Please Input a Valid: Amount of Allocated CPU Threads (between 1 and 50)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  if (type(enable_print) is not bool):
    print('Please Input a Valid: Enable Print (True or False)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  if (type(ram_download) is not bool):
    print('Please Input a Valid: Ram Download (True or False)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  if (type(bad_bits) is not list) or (len(bad_bits) < 1):
    print('Please Input a Valid: Bad Bits (list)')
    print('Read Documentation: https://github.com/huntbrooks85/Xcavation')
    return False

  return True
# ------------------------------------------------------ #



# DataClass Management
# ------------------------------------------------------ #
@dataclass
class genspec_profile:
    r_fwhm: float = 2
    pmra: float = 0
    pmdec: float = 0
    mjd: float = 61000
    save_data: bool = False
    output_path: str = 'xcavation_save'
    threads: int = 1
    enable_print: bool = True
    ram_download: bool = False
    bad_bits: list = field(default_factory=lambda: [0,1,2,4,6,7,9,10,11,14,15,17,19])
    
# HIERARCH MP_TRANSIENT = 0 / Transient detected during SUR                       
# HIERARCH MP_OVERFLOW = 1 / Overflow detected during SUR                         
# HIERARCH MP_SUR_ERROR = 2 / Error in onboard processing                         
# HIERARCH MP_PHANTOM = 4 / Phantom pixel                                         
# HIERARCH MP_REFERENCE = 5 / Reference pixel                                     
# HIERARCH MP_NONFUNC = 6 / Permanently unusable                                  
# HIERARCH MP_DICHROIC = 7 / Low efficiency due to dichroic                       
# HIERARCH MP_MISSING_DATA = 9 / Onboard data lost                                
# MP_HOT  =                   10 / Hot pixel                                      
# MP_COLD =                   11 / Anomalously low signal                         
# HIERARCH MP_FULLSAMPLE = 12 / Pixel full sample history is available            
# HIERARCH MP_PHANMISS = 14 / Phantom correction was not applied                  
# HIERARCH MP_NONLINEAR = 15 / Nonlinearity correction cannot be applied reliably
# HIERARCH MP_PERSIST = 17 / Persistent charge above threshold                    
# HIERARCH MP_OUTLIER = 19 / Pixel flagged by Detect Outliers                     
# HIERARCH MP_SOURCE = 21 / Pixel mapped to a known source       
# ------------------------------------------------------ #
    
    

# SphereX 402 Error Retry Function
# ------------------------------------------------------ #
def retry(func, *args, retries=10, delay=0.5, **kwargs):
    for attempt in range(1, retries + 1):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if attempt < retries:
                time.sleep(delay)
            else:
                return None
# ------------------------------------------------------ #



# Multi-Thread the SphereX Query Function
# ------------------------------------------------------ #
def genspec(ra, dec, style, config: genspec_profile):

    r_fwhm = config.r_fwhm
    pmra = config.pmra
    pmdec = config.pmdec
    mjd = config.mjd
    save_data = config.save_data
    output_path = config.output_path
    threads = config.threads
    enable_print = config.enable_print
    ram_download = config.ram_download
    bad_bits = config.bad_bits



    verify = variable_verify(ra, dec, style, r_fwhm, pmra, pmdec, mjd, save_data, output_path, threads, enable_print, ram_download, bad_bits)
    if verify is False:
      raise ValueError('Variables Assigned Incorrectly. Please Read Documentation.')




    # ----- Proper Motion Propagation ----- #
    time_passed = time_mjd(mjd) # years since observation
    ra_deg, dec_deg = proper_motion(ra, dec, pmra, pmdec, time_passed) # Adjust for Proper Motion
    if enable_print is True:
      print('\n+====================================================================+')
      print('       Starting xcavation.genspec SPHEREx Spectrophotometry Tool')
      print('+====================================================================+')
      print(f"Queried Coordinates: R.A. = {round(ra_deg, 6)} (deg), Decl. = {round(dec_deg, 6)} (deg)") # Print Adjusted Coordinates
      print(f"Query Time: {datetime.now()}")
      print('+====================================================================+\n')
    # ------------------------------------- #

    # ----- SphereX Query ----- #
    coord = SkyCoord(ra_deg, dec_deg, unit='deg') # Convert to astropy Units
    results = Irsa.query_sia(pos=(coord, 10 * u.arcsec), collection='spherex_qr2') # Query SphereX QR2 (Will Need Updates)
    urls = results["access_url"] # Gets the Image URLs
    D1, D2, D3, D4, D5, D6 = 0, 0, 0, 0, 0, 0
    for url in urls:
      num_detector = int(url.split('D')[len(url.split('D')) - 1][0])
      if num_detector == 1:
        D1 += 1
      if num_detector == 2:
        D2 += 1
      if num_detector == 3:
        D3 += 1
      if num_detector == 4:
        D4 += 1
      if num_detector == 5:
        D5 += 1
      if num_detector == 6:
        D6 += 1
    # ------------------------- #



    # ----- Multi-Thread Wavelength and Aperature Calc. ----- #
    output = []
    with ThreadPoolExecutor(max_workers=threads) as executor:

        # Multi-Thread w/ tqdm
        if style == 'aperture':
            futures = [executor.submit(retry, spherex_aperature_phot, url, coord, r_fwhm, ram_download, bad_bits) for url in urls]
        elif style == 'psf':
            print('Please Wait Till Next Major Release')

        if enable_print is True:
          print('+====================================================================+')
          print(f'{len(urls)} Images Were Found in SPHEREx')
        for future in tqdm(futures, desc="Calculating", ncols=70, colour = '#C5D89D', bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]", ascii="░▒█", disable=not enable_print):
            result = future.result()

            if result is None: 
              raise ValueError('Network Connectivity Issue. Please Try Again')

            output.append({
                "wavelength": result["wavelength"],
                "delta_lambda": result["delta_lambda"],
                "flux": result["flux"],
                "flux_err": result["flux_err"],
                "SNR": result["SNR"],
                "flux_cutout": result["flux_cutout"],
                "aperture": result["aperture"],
                "annulus": result["annulus"],
                "x_loc": result["x_loc"],
                "y_loc": result["y_loc"],
                "ap_radius": result['ap_radius'],
                "inner_annulus": result['inner_annulus'],
                "outer_annulus": result['outer_annulus'],
                "url": result["url"],
                "mjd": result["mjd"]
            })
        if enable_print is True:
          print('Spectral Calulcation Complete')
          print('+====================================================================+\n')
    # ------------------------------------------------------- #
    output = pd.DataFrame(output)



    # ----- Finder Chart and Spectral Plots ----- #
    if save_data is True:

      finder_chart(output, output_path + '.pdf')

      if enable_print is True:
        print('+====================================================================+')
        print('Finder Chart Was Successfully Saved')


        print(f'Spectral Plot of: R.A. = {round(ra, 3)} (deg) and Decl. = {round(dec, 3)} (deg)')
      spectra_plot(output)

      lambda_list = np.array(output['wavelength'])
      flux_list = np.array(output['flux'])
      flux_err_list = np.array(output['flux_err'])
      delta_lambda = np.array(output['delta_lambda'])

      sort_idx = np.argsort(lambda_list)
      lambda_list = lambda_list[sort_idx]
      flux_list = flux_list[sort_idx]
      flux_err_list = flux_err_list[sort_idx]
      delta_lambda = delta_lambda[sort_idx]

      t = Table([lambda_list, flux_list, flux_err_list, delta_lambda], names=("lambda", "flux", "flux_err", "delta_lambda"))
      with open(output_path + '.txt', 'w') as f:
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
          f.write(f'# Photometry Type: {style}\n')
          f.write(f'# Average SNR: {round(np.nanmean(output['flux']/output['flux_err']), 2)}\n')
          f.write(f'# Average Flux Error: {round(np.nanmean(output['flux_err']))}\n')
          f.write(f'# Images per detector: D1={D1}, D2={D2}, D3={D3}, D4={D4}, D5={D5}, D6={D6}\n')
          f.write(f'# References: https://spherex.caltech.edu/page/publications; https://github.com/huntbrooks85/Xcavation\n')
          f.write(f'# Column 1: Wavelength (micron)\n')
          f.write(f'# Column 2: Flux (uJy)\n')
          f.write(f'# Column 3: Flux Error (uJy)\n')
          f.write(f'# Column 4: Delta Lambda (microns)\n')

          # Now write the table
          t.write(f, format='ascii', overwrite=True)
      if enable_print is True:
        print('Spectrum Table Was Successfully Saved')
        print('+====================================================================+\n')
    # # ------------------------------------------- #


    if enable_print is True:
      print('+====================================================================+')
      print('              xcavation Spectrophotometry Tool Summary')
      print('+====================================================================+')

      print(f"Queried Coordinates: R.A. = {round(ra_deg, 6)} (deg), Decl. = {round(dec_deg, 6)} (deg)") # Print Adjusted Coordinates
      print(f"SIMBAD Query: https://simbad.u-strasbg.fr/simbad/sim-coo?Coord={ra}+{dec}&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&CooDefinedFrames=none&Radius=0.5&Radius.unit=arcmin&submit=submit+query&CoordList=")
      print(f"Query Finished: {datetime.now()}\n")

      print(f"Number of Images Found: {len(urls)}")
      print(f"Images Per Detector: D1={D1}, D2={D2}, D3={D3}, D4={D4}, D5={D5}, D6={D6}\n")

      print(f'Average SNR: {round(np.nanmean(output['flux']/output['flux_err']), 2)}')
      print(f'Average Flux Error: {round(np.nanmean(output['flux_err']))}\n')

      print(f'Wavelength Range: {round(np.nanmin(output['mjd']), 1)} to {round(np.nanmax(output['mjd']), 1)}')
      print(f"Flux Range: {round(np.nanmin(output['flux']), 3)} to {round(np.nanmax(output['flux']), 3)}")
      print(f'Observation Range: {round(np.nanmin(output['wavelength']), 3)} to {round(np.nanmax(output['wavelength']), 3)}')

      if save_data == True:
        print(f'\nData Saved To: {output_path}')
      print('+====================================================================+')


    # Returns The Output Dictionary
    return output
# ------------------------------------------------------ #

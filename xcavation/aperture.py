
#-----------------------------------------------------------------------#
# xcavation.aperture v1.0.1
# By Hunter Brooks, at UToledo, Toledo: Apr. 08, 2026
#
# Purpose: Perform Aperture Photometry on SphereX Data
#-----------------------------------------------------------------------#



# Import Data Management
# ------------------------------------------------------ #
import numpy as np
import astropy.units as u
from astropy.io import fits
from scipy.interpolate import griddata
# ------------------------------------------------------ #



# Import Internal Modules
# ------------------------------------------------------ #
from .motion import time_mjd, proper_motion
# ------------------------------------------------------ #



# Import WCS, Photometry, and Plotting
# ------------------------------------------------------ #
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from scipy.ndimage import median_filter
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry, ApertureStats
# ------------------------------------------------------ #



# Import API Mangament Tools
# ------------------------------------------------------ #
import requests
from io import BytesIO
# ------------------------------------------------------ #




# Table of Resolving Powers Based on Wavelength
# ------------------------------------------------------ #
def resolving_table(wave):
  """
  Return the spectral resolution element (Δλ) for a given wavelength
  based on the SPHEREx instrument resolving powers. The resolving power,
  R = λ / Δλ, is taken from Akeson et al. (2025):
  https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp_QR_v1.0.pdf.

  Parameters
  ----------
    wave: Wavelength inputted to find corresponding R value (float)

  Returns
  -------
    delta_lambda: Spectral wavelength resolution at a given wavelength (float)
  """


  # ----- Table of R Values ----- #
  # Use Inputted Wavelength to Find Relevant R value
  if wave < 1.11:
    R = 39
  if 1.11 <= wave < 2.42:
    R = 41
  elif 2.42 <= wave < 3.82:
    R = 35
  elif 3.82 <= wave < 4.42:
    R = 110
  elif 4.42 <= wave:
    R = 130
  # ----------------------------- #



  # Calculate Delta Lambda Using Wavelength and R Value
  return wave/R
# ------------------------------------------------------ #




# Perform SphereX Wavelength and Aperature Calculation
# ------------------------------------------------------ #
def spherex_aperature_phot(url_og, coord, pm, mjd,  # Coords
                           r_fwhm, r_annulus_in, r_annulus_out, # Radii
                           ram_download, # Download Types
                           clean_type, bad_bits, # Bad Pixel Fixes
                           background_type, # Background Subraction
                           cutout_size, # Cutout Size
                           zodi_subtract, # ZODI Subtraction
                           sigclip_sigma, sigclip_maxiters): # SigmaClip Valus
    """
    Perform SPHEREx aperture photometry for a given sky position and epoch.
    Downloads SPHEREx FITS file, propagates the input coordinates for
    proper motion, creates cutouts of flux, variance, and flags. Performs
    aperture photometry with a local background annulus. Spectral wavelength
    and resolution are also calulcated alongside counting bad pixels within
    the aperture.

    Parameters
    ----------
      url: URL of SPHEREx FITS file (str)
      coord : Input sky coordinates (RA, Dec) in degrees (list)
      pm: Proper motion in arcseconds per year (list)
      mjd: MJD corresponding to the input coordinates (float)
      r_fwhm: Aperture radius in units of PSF FWHM (float)
      r_annulus_in: Inner radius of background annulus in units of FWHM (float)
      r_annulus_out:  Outer radius of background annulus in units of FWHM (float)
      ram_download: Download the FITS entirely into memory (bool)
      bad_bits: FLAG extension to mask as bad pixels (list)
      background_type: Mean, Median, or Mode of Background Annulus (str)
      cutout_size: Cutout Size in pixels (int)
      zodi_subtract: Whether ZODI light is subtracted (bool)
      sigclip_sigma: Background Astropy.SigmaClip sigma level (float)
      sigclip_maxiters: Background Astropy.SigmaClip max iterations (float)

    Returns
    -------
      result : dict
          Dictionary containing the following keys:
          - 'wavelength': Wavelength at target position (float)
          - 'delta_lambda': Spectral resolution element Δλ (float)
          - 'flux': Aperture flux (float)
          - 'flux_err': Propagated flux uncertainty (float)
          - 'flag_count': Counts of flagged pixels per bad bit (dict)
          - 'SNR': Signal-to-noise ratio of the aperture flux (float)
          - 'flux_cutout': Flux cutout (array)
          - 'flag_cutout': FLAG cutout (array)
          - 'aperture': Binary aperture mask (bool array)
          - 'annulus': Binary annulus mask (bool array)
          - 'x_loc', 'y_loc': Pixel coordinates of target in cutout (float)
          - 'ap_radius': Aperture radius (float)
          - 'inner_annulus', 'outer_annulus': Annulus radii (float)
          - 'url': Input URL (str)
          - 'mjd': MJD of observation (float)
    """



    # ----- Download FITS entirely in RAM ----- #
    if ram_download is True:
      resp = requests.get(url_og)
      resp.raise_for_status()
      url = BytesIO(resp.content)
    else:
      url = url_og
    # ----------------------------------------- #



    # ---------------- Open File ---------------- #
    with fits.open(url, lazy_load_hdus=True) as hdul:

        # ---------------- WCS ---------------- #
        # Gets Image Header and States Pixel Scale
        spectral_header = hdul["IMAGE"].header.copy()
        pixel_scale_deg=0.001708333 # Pixel Scale From Documentation

        # Blocks SIP warning message doesnt show
        for key in list(spectral_header.keys()):
            if key.startswith(("A_", "B_", "AP_", "BP_")):
                del spectral_header[key]

        # Gets Flag Data and Finds All Selected Bad Bits
        flags = hdul['FLAGS'].data # Flag Data
        bad_bitmask = sum(1 << bit for bit in bad_bits) # Bad Bits Data

        # Apply bad pixel mask and sets to 0
        bad_pixel_mask = (flags & bad_bitmask) != 0

        # Obtain Header Data
        celestial_wcs = WCS(hdul[1].header.copy()) # Normal WCS
        spectral_wcs = WCS(spectral_header, fobj=hdul, key="W") # Wavelength WCS

        # Obtain Flux Image
        flux_image_temp =  hdul[1].data.astype(float) # Flux Image
        if zodi_subtract == True:
          zodi_image_temp = hdul[4].data.astype(float) # ZODI Image
          flux_image = flux_image_temp - zodi_image_temp # Image - ZODI
        else:
          flux_image = flux_image_temp # Image
        flux_image = flux_image.astype(np.float32, copy=True)

        # Obtain Variance Image
        variance_image = hdul['VARIANCE'].data.astype(float) # VARIANCE Image
        variance_image = variance_image.astype(np.float32, copy=True)
        # ------------------------------------- #



        # ----- Proper Motion Propagation ----- #
        ra_in, dec_in = coord[0], coord[1] # Inputted R.A.
        pmra_in, pmdec_in = pm[0], pm[1] # Inputted Decl.
        mjd_obs = spectral_header['MJD-OBS'] # Observed MJD
        time_passed = time_mjd(mjd, mjd_obs) # Years since observation
        ra_deg, dec_deg = proper_motion(ra_in, dec_in, # Proper Motion Prop.
                                      pmra_in, pmdec_in, time_passed)
        coord = SkyCoord(ra_deg, dec_deg, unit='deg') # Convert to astropy Units
        # ------------------------------------- #



        # ---------------- Interpolated Masking ---------------- #
        if clean_type == 'interp_mask':
          ny, nx = flux_image.shape
          y, x = np.mgrid[0:ny, 0:nx]

          # Good pixels only
          good = (~bad_pixel_mask) & np.isfinite(flux_image)

          # Coordinates and values of good pixels
          points = np.column_stack((x[good], y[good]))
          values_flux = flux_image[good]
          values_var = variance_image[good]

          # Interpolate over full grid
          interp_flux = griddata(points, values_flux, (x, y), method='nearest')
          interp_var = griddata(points, values_var, (x, y), method='nearest')
          flux_image[bad_pixel_mask] = interp_flux[bad_pixel_mask]
          variance_image[bad_pixel_mask] = interp_var[bad_pixel_mask]
        # ------------------------------------------------------ #



        # ---------------- Median Masking ---------------- #
        if clean_type == 'median_mask':
          bad = flags != 0
          med = median_filter(flux_image, size=5, mode='mirror')
          flux_image[bad] = med[bad]
        # ------------------------------------------------ #
        del bad_pixel_mask # Delete Mask to Save on RAM




        # ---------------- Convert Units ---------------- #
        # Converts Units of Flux and Variance Cutout
        pixel_area_sr = (pixel_scale_deg * np.pi / 180)**2 # Pixel^2 to Sr
        flux_ujy = flux_image * 1e6 * pixel_area_sr * 1e6 # Flux Convert
        var_ujy = variance_image * (1e6 * pixel_area_sr * 1e6)**2 # Variance Convert
        # ----------------------------------------------- #




        # ----- Set-Up Aperture Photometry ----- #
        # Gets x, y Pixel Locations and FWHM
        x_query, y_query = celestial_wcs.world_to_pixel(coord) # Cutouts Coords.
        fwhm = hdul[1].header['PSF_FWHM'] # Image FWHM

        # Aperature Set-Up
        ap_radius_arcsec = (r_fwhm * fwhm * u.arcsec).to(u.deg).value
        ap_radius_pix = ap_radius_arcsec / pixel_scale_deg # Aperture Radius
        aperture = CircularAperture((x_query, y_query), # Query Positon
                                    r=ap_radius_pix) # Radius

        # Annulus Set-up
        an_inner_arcsec = (r_annulus_in * fwhm * u.arcsec).to(u.deg).value
        an_outer_arcsec = (r_annulus_out * fwhm * u.arcsec).to(u.deg).value
        an_inner = an_inner_arcsec / pixel_scale_deg # Annulus Inner Radius
        an_outer = an_outer_arcsec / pixel_scale_deg # Annulus Outer Radius
        annulus_ap = CircularAnnulus((x_query, y_query), # Query Position
                                     r_in=an_inner, r_out=an_outer) # Radii
        # -------------------------------------- #



        # ----- Background Frame ----- #
        # Calculates Background Frame
        sigclip = SigmaClip(sigma=sigclip_sigma, maxiters=sigclip_maxiters) # Sigma Clip For Background
        if background_type == 'mean':
          bkg_per_pix = (ApertureStats(flux_ujy, annulus_ap,
                          sigma_clip=sigclip)).mean # Average Annulus Background
        elif background_type == 'median':
          bkg_per_pix = (ApertureStats(flux_ujy, annulus_ap,
                          sigma_clip=sigclip)).median # Median Annulus Background
        elif background_type == 'mode':
          bkg_per_pix = (ApertureStats(flux_ujy, annulus_ap,
                          sigma_clip=sigclip)).mode # Mode Annulus Background

        ap_mask = aperture.to_mask(method='exact').to_image(flux_ujy.shape)

        # Calculate Aperture Area (w/ or w/o mask)
        if clean_type == 'mask':
          good_ap_mask = ap_mask * (~bad_pixel_mask) # w/ mask
        else:
          good_ap_mask = ap_mask # w/o mask

        # Calcualtes Final Background Frame
        area = np.sum(good_ap_mask) # Total Aperture Area
        if str(area) == 'None':
          return {
                  "wavelength": np.nan,
                  "delta_lambda": np.nan,
                  "flux": np.nan,
                  "flux_err": np.nan,
                  "flag_count": np.nan,
                  "flag": np.nan,
                  "SNR": np.nan,
                  "flux_cutout": np.nan,
                  "flag_cutout": np.nan,
                  "aperture": np.nan,
                  "annulus": np.nan,
                  "x_loc": np.nan,
                  "y_loc": np.nan,
                  "ap_radius": np.nan,
                  "inner_annulus": np.nan,
                  "outer_annulus": np.nan,
                  "url": url_og,
                  "mjd": np.nan
                  }
        total_bkg = bkg_per_pix * area # Background Frame in Aperture Area
        # ---------------------------- #



        # ----- Final Aperture Photometry ----- #
        # Calculates the Aperature Photometry (w/ or w/o mask)
        if clean_type == 'mask':
          phot_table = aperture_photometry(flux_ujy, aperture,
                                           mask = bad_pixel_mask) # w/ mask
        else:
          phot_table = aperture_photometry(flux_ujy, aperture) # w/o mask
        flux_ap = phot_table['aperture_sum'][0] - total_bkg # Final Aperture
        # ------------------------------------ #



        # ----- Aperture Error Photometry ----- #
        # Aperture Error (w/ or w/o mask)
        if clean_type == 'mask':
          phot_var_table = aperture_photometry(var_ujy, aperture,
                                                mask = bad_pixel_mask) # w/ mask
        else:
          phot_var_table = aperture_photometry(var_ujy, aperture) # w/o mask

        # Calculate Progated Errors
        var_flux_ap = phot_var_table['aperture_sum'][0] # Aperture Error
        top_var = ((ApertureStats(flux_ujy, annulus_ap,
                                  sigma_clip=sigclip)).std)**2
        bot_var = (annulus_ap.area_overlap(var_ujy))
        var_bkg_per_pix = top_var/bot_var # Background Error
        flux_err = np.sqrt(var_flux_ap+(area**2)*var_bkg_per_pix) # Total Error
        # ------------------------------------- #



        # ----- Flag Aperture Count ----- #
        # Gets the Pixels in the Aperture
        ap_pixels = ap_mask > 0
        ap_flags = flags[ap_pixels]

        # Queried Array
        total_bit_counts = {} # Empty Array
        bit_counts = 0
        total_bad_bits = [0, 1, 2, 4, 5, 6, 7,
                          9, 10, 11, 12, 14, 15,
                          17, 19, 21]

        # Counts Each Flag in the Aperture
        for bit in total_bad_bits:
            bitmask = 1 << bit # Convert to bit
            total_bit_counts[bit] = np.count_nonzero(ap_flags & bitmask) # Count

        # Counts Each Flag in the Aperture
        for bit in bad_bits:
            bitmask = 1 << bit # Convert to bit
            bit_counts += np.count_nonzero(ap_flags & bitmask) # Count
        # ------------------------------- #



        # ----- Wavelength ----- #
        # General Set-Up
        y, x = np.indices(hdul["IMAGE"].data.shape) # Obtain Detector Size
        wave, *_ = spectral_wcs.pixel_to_world(x, y) # Wavelength WCS
        wave = np.asarray(wave) # Make Wavelength Arrary

        # Obtains Wavelength Average
        x_query_total, y_query_total = celestial_wcs.world_to_pixel(coord)
        y = int(np.clip(y_query_total, 0, wave.shape[0] - 1))
        x = int(np.clip(x_query_total, 0, wave.shape[1] - 1))
        wave_point = wave[y, x] # Wavelength
        # ---------------------- #
    # ------------------------------------------- #



    # ----- Output Aperture/Annulus Array ----- #
    ap_mask_obj = aperture.to_mask(method='exact') # Aperture Mask
    an_mask_obj = annulus_ap.to_mask(method='exact') # Annulus Mask

    # Convert to image-sized arrays
    ap_mask = ap_mask_obj.to_image(flux_image.shape) # Aperture
    an_mask = an_mask_obj.to_image(flux_image.shape) # Annulus
    # ----------------------------------------- #



    # ----- Output Dictionary ----- #
    return {
        "wavelength": wave_point,
        "delta_lambda": resolving_table(wave_point),
        "flux": flux_ap,
        "flux_err": flux_err,
        "flag_count": total_bit_counts,
        "flag": bit_counts,
        "SNR": (flux_ap/flux_err),
        "flux_cutout": flux_ujy,
        "flag_cutout": flags.data,
        "aperture": ap_mask,
        "annulus": an_mask,
        "x_loc": x_query,
        "y_loc": y_query,
        "ap_radius": r_fwhm * fwhm,
        "inner_annulus": r_annulus_in * fwhm,
        "outer_annulus": r_annulus_out * fwhm,
        "url": url_og,
        "mjd": spectral_header['MJD-OBS']
    }
    # ----------------------------- #
# ------------------------------------------------------ #

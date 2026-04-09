
#-----------------------------------------------------------------------#
# xcavation.quality v1.0.2
# By Hunter Brooks, at UToledo, Toledo: Apr. 08, 2026
#
# Purpose: Plotting Tools for SphereX Data
#-----------------------------------------------------------------------#



# Import Data Management
# ------------------------------------------------------ #
import math
import numpy as np
# ------------------------------------------------------ #



# Import WCS, Photometry, and Plotting
# ------------------------------------------------------ #
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import ListedColormap
# ------------------------------------------------------ #



# Plot Style
# ------------------------------------------------------ #
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
  'font.family': 'STIXGeneral',
  'font.size': 14,
  'axes.labelsize': 20,
  'axes.titlesize': 25,
  'xtick.labelsize': 15,
  'ytick.labelsize': 15,
  'axes.linewidth': 1.2,
  'legend.fontsize': 12
                  })
# ------------------------------------------------------ #



# Function to Generate Q.A. Finder Charts
# ------------------------------------------------------ #
def finder_chart(output, pdf_path):
  """
  Generate a multi-panel finder chart PDF from SPHEREx aperture photometry results.

  This function creates a grid of cutout images for multiple wavelengths,
  overlays the aperture and annulus regions, and marks the target position.
  The panels are sorted by wavelength. The final figure is saved as a PDF.

  Parameters
  ----------
    output: dict
        Dictionary containing SPHEREx photometry:
        - 'flux_cutout': list of 2D flux arrays (μJy)
        - 'aperture': list of 2D aperture masks
        - 'annulus': list of 2D annulus masks
        - 'x_loc', 'y_loc': lists of pixel coordinates of the target in cutouts
        - 'wavelength': list or array of wavelengths (μm)
        - 'flux': list of aperture flux values (μJy)
        - 'ap_radius', 'inner_annulus', 'outer_annulus': lists of radii in arcsec
    pdf_path: File path where the PDF finder chart will be saved (str)

  Returns
  -------
    None: Figure save to the specified path
  """



  # ------- Set-Up Subplots ------- #
  n_images = len(output['flux_cutout']) # Total number of flux arrays
  ncols = 5 # Number of Columns (pre-set)
  nrows = math.ceil(n_images / ncols) # Number of Rows (based on total)
  # ------------------------------- #



  # ------- Get Required Data From Input Array ------- #
  # Sort By Wavelength
  wavelengths = np.array(output['wavelength']) # Original Wavelength List
  sorted_idx = np.argsort(wavelengths) # Sorted Wavelength List

  # Sorts All Needed Variables
  flux_cutouts_sorted = [output['flux_cutout'][i] for i in sorted_idx] # Cutouts
  apers_sorted = [output['aperture'][i] for i in sorted_idx] # Aperture Mask
  annuls_sorted = [output['annulus'][i] for i in sorted_idx] # Annulus Mask
  x_locs_sorted = [output['x_loc'][i] for i in sorted_idx] # X Pixel Location
  y_locs_sorted = [output['y_loc'][i] for i in sorted_idx] # Y Pixel Location
  wavelengths_sorted = wavelengths[sorted_idx] # Sorted Wavelength
  fluxes_sorted = [output['flux'][i] for i in sorted_idx] # Flux
  aperture_sorted = [output['ap_radius'][i] for i in sorted_idx] # Aperture Radius
  inner_annulus_sorted = [output['inner_annulus'][i] for i in sorted_idx] # Inner Annulus Radius
  outer_annulus_sorted = [output['outer_annulus'][i] for i in sorted_idx] # Outer Annulus Radius
  # -------------------------------------------------- #



  # ------- Create Subplot ------- #
  fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows))
  # ------------------------------ #



  # ------------ Plot Each Cutout ------------ #
  axes = axes.ravel()
  valid_count = 0
  for i in range(n_images):
      # Get Relavent Data
      cutouts = flux_cutouts_sorted[i] # Cutout
      wavelength = wavelengths_sorted[i] # Wavelength
      flux = fluxes_sorted[i] # Flux

      # Verify that each Cutout, Wavelength, and Flux is Valid
      if (cutouts is not None) and (np.isfinite(wavelength)) and (np.isfinite(flux)):
        ax = axes[valid_count] # Counts Valid Axes

        # Get the Valid Data
        apers = apers_sorted[i] # Aperture Mask
        annuls = annuls_sorted[i] # Annulus Mask
        x_loc = x_locs_sorted[i] # Queried X Location
        y_loc = y_locs_sorted[i] # Queried Y Location
        aperture = aperture_sorted[i] # Aperture Radius
        inner_annulus = inner_annulus_sorted[i] # Inner Annulus Radius
        outer_annulus = outer_annulus_sorted[i] # Outer Annulus Radius

        # Plots the Image Data
        ax.scatter(x_loc, y_loc, c='yellow',marker='*', s = 100) # Queried Point
        ax.imshow(cutouts, cmap='Greys',
          vmin=np.nanpercentile(cutouts, 1), vmax=np.nanpercentile(cutouts, 99))
        ax.imshow(np.ma.masked_where(apers == 0, apers),
                  cmap=ListedColormap(['red']), alpha=0.2) # Aperture Mask
        ax.imshow(np.ma.masked_where(annuls == 0, annuls),
                  cmap=ListedColormap(['blue']), alpha=0.2) # Annulus Mask

        # Creates Matplotlib Circles for Aperture and Annulus Radii
        ap = Circle((x_loc, y_loc), aperture/6.15,
                    color='red', lw = 2, fill = False) # Aperture Radius
        in_an = Circle((x_loc, y_loc), inner_annulus/6.15,
                    color='blue', lw = 2, fill = False) # Inner Annulus Radius
        out_an = Circle((x_loc, y_loc), outer_annulus/6.15,
                    color='blue', lw = 2, fill = False) # Outer Annulus Radius

        # Plots Each Circle
        ax.add_patch(ap) # Aperture
        ax.add_patch(in_an) # Inner Annulus
        ax.add_patch(out_an) # Outer Annulus

        # General Subplot Looks
        ax.set_title(f"$λ$: {wavelength:.4f} μm and $F_ν$: {flux:.0f} μJy", fontsize=12)
        ax.axis('off') # Turn off axis
        valid_count += 1 # Counts Valid Plots
  # ------------------------------------------ #



  # ------------ Fix Plots ------------ #
  # Remove Unused Plots
  for ax in axes[valid_count:]:
      fig.delaxes(ax)

  # Adjust Each Subplot
  plt.subplots_adjust(left=0.02, right=0.98,
                      bottom=0.02, top=0.95,
                      wspace=0.05, hspace=0.15)
  # -------------------------------------- #




  # ------------ Save Figure ------------ #
  plt.savefig(pdf_path, bbox_inches='tight') # Save Figure
  plt.close() # Close Figure to Save Memory
  # ------------------------------------- #
# ------------------------------------------------------ #



# Plots Spectra Without Saving
# ------------------------------------------------------ #
def spectra_plot(output):
  """
  Plot a SPHEREx spectrum with flux uncertainties.

  Parameters
  ----------
    output : dict
        Dictionary containing SPHEREx photometry results. Expected keys:
        - 'wavelength': list or array of wavelengths in microns (μm)
        - 'flux': list or array of flux values in μJy
        - 'flux_err': list or array of flux uncertainties in μJy

  Returns
  -------
    None
        Displays the spectrum plot using matplotlib. Does not return data.
  """



  # ------- Create Plot ------- #
  plt.figure(figsize=(6,4))
  # --------------------------- #



  # ------- Get Relevant Data ------- #
  wavelengths = np.array(output['wavelength']) # Wavelength
  fluxes = np.array(output['flux']) # Flux
  errors = np.array(output['flux_err']) # Flux Error
  flags = np.array(output['flag']) # Flag Count
  # --------------------------------- #



  # ------- Sort Data Based On Wavelength ------- #
  sort_idx = np.argsort(wavelengths) # Sort Based on Wavelength
  wavelengths = wavelengths[sort_idx] # Wavelength Sort
  fluxes = fluxes[sort_idx] # Flux Sort
  errors = errors[sort_idx] # Flux Error Sort
  flags = flags[sort_idx] # Flag Count Sort
  finite = np.isfinite(fluxes) # Ensures that all Flus Data is Real
  # --------------------------------------------- #



  # ------- Plot Flux and Errors ------- #
  plt.errorbar(wavelengths[finite], fluxes[finite], # Plot Flux
               yerr=errors[finite], fmt='o', color='slategray',
               ecolor='slategray', elinewidth=1.5, capsize=3,
               markersize=1, alpha = 0.5)
  plt.step(wavelengths[finite], fluxes[finite], # Plot Errors
           color='k', alpha=0.6, where = 'mid')
  sc = plt.scatter(wavelengths[finite], fluxes[finite],
                 c=flags[finite], cmap='jet',
                 alpha=0.8, s=15) # Plot Flags
  # ------------------------------------ #



  # ------- Make Plot Pretty ------- #
  plt.colorbar(sc, label='Flag Count') # Colorbar

  plt.xlabel("Wavelength [μm]", fontsize=14) # X Label
  plt.ylabel("Flux [μJy]", fontsize=14) # Y Label

  plt.grid(alpha=0.2) # Plot Grid
  plt.minorticks_on() # Plot Minor Tickmarks
  plt.tight_layout() # Forces Plot to Be Tight on all Paddings

  # Make Y Limit Using Percentiles
  bot = np.nanpercentile(fluxes[finite], 1) - 100
  top =  np.nanpercentile(fluxes[finite], 95) + 100
  plt.ylim(bot, top)

  # Show Plot
  plt.show()
  # -------------------------------- #
# ------------------------------------------------------ #

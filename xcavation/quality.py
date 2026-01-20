
#-----------------------------------------------------------------------#
# xcavation.quality v0.3.1
# By Hunter Brooks, at UToledo, Toledo: Jan. 20, 2026
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



# Generate Finder Chart
#----------------------------------------------------------------------- #
def finder_chart(output, pdf_path):
  n_images = len(output['flux_cutout'])
  ncols = 5
  nrows = math.ceil(n_images / ncols)

  wavelengths = np.array(output['wavelength'])
  sorted_idx = np.argsort(wavelengths)
  flux_cutouts_sorted = [output['flux_cutout'][i] for i in sorted_idx]
  apers_sorted = [output['aperture'][i] for i in sorted_idx]
  annuls_sorted = [output['annulus'][i] for i in sorted_idx]
  x_locs_sorted = [output['x_loc'][i] for i in sorted_idx]
  y_locs_sorted = [output['y_loc'][i] for i in sorted_idx]
  wavelengths_sorted = wavelengths[sorted_idx]
  fluxes_sorted = [output['flux'][i] for i in sorted_idx]

  aperture_sorted = [output['ap_radius'][i] for i in sorted_idx]
  inner_annulus_sorted = [output['inner_annulus'][i] for i in sorted_idx]
  outer_annulus_sorted = [output['outer_annulus'][i] for i in sorted_idx]

  fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows))

  axes = axes.ravel()
  valid_count = 0
  for i in range(n_images):
    cutouts = flux_cutouts_sorted[i]
    wavelength = wavelengths_sorted[i]
    flux = fluxes_sorted[i]

    if (cutouts is not None) and (np.isfinite(wavelength)) and (np.isfinite(flux)):
      ax = axes[valid_count]

      apers = apers_sorted[i]
      annuls = annuls_sorted[i]
      x_loc = x_locs_sorted[i]
      y_loc = y_locs_sorted[i]
      aperture = aperture_sorted[i]
      inner_annulus = inner_annulus_sorted[i]
      outer_annulus = outer_annulus_sorted[i]

      ax.scatter(x_loc, y_loc, c='yellow', marker='*', s = 100)
      ax.imshow(cutouts, cmap='Greys', vmin=np.nanpercentile(cutouts, 1), vmax=np.nanpercentile(cutouts, 99))
      ax.imshow(np.ma.masked_where(apers == 0, apers), cmap=ListedColormap(['red']), alpha=0.2)
      ax.imshow(np.ma.masked_where(annuls == 0, annuls), cmap=ListedColormap(['blue']), alpha=0.2)

      ap = Circle((x_loc, y_loc), aperture/6.15, color='red', lw = 2, fill = False)
      in_an = Circle((x_loc, y_loc), inner_annulus/6.15, color='blue', lw = 2, fill = False)
      out_an = Circle((x_loc, y_loc), outer_annulus/6.15, color='blue', lw = 2, fill = False)

      ax.add_patch(ap)
      ax.add_patch(in_an)
      ax.add_patch(out_an)

      ax.set_title(f"$λ$: {wavelength:.4f} μm and $F_λ$: {flux:.0f}", fontsize=12)
      ax.axis('off')
      valid_count += 1

  for ax in axes[valid_count:]:
    fig.delaxes(ax)

  plt.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.95, wspace=0.05, hspace=0.15)

  plt.savefig(pdf_path, bbox_inches='tight')
  plt.close()
# ----------------------------------------------------------------------- #



# Plot Spectra
#----------------------------------------------------------------------- #
def spectra_plot(output):
  plt.figure(figsize=(6,4))

  wavelengths = np.array(output['wavelength'])
  fluxes = np.array(output['flux'])
  errors = np.array(output['flux_err'])

  sort_idx = np.argsort(wavelengths)
  wavelengths = wavelengths[sort_idx]
  fluxes = fluxes[sort_idx]
  errors = errors[sort_idx]
  finite = np.isfinite(fluxes)

  plt.errorbar(wavelengths[finite], fluxes[finite], yerr=errors[finite], fmt='o', color='slategray', ecolor='slategray', elinewidth=1.5, capsize=3, markersize=1, alpha = 0.5)

  plt.step(wavelengths[finite], fluxes[finite], color='#656D3F', alpha=0.6, where = 'mid')

  plt.xlabel("Wavelength [μm]", fontsize=14)
  plt.ylabel("Flux [μJy]", fontsize=14)

  plt.grid(alpha=0.2)
  plt.minorticks_on()
  plt.tight_layout()

  plt.ylim(-100, np.nanpercentile(fluxes[finite], 99))
  plt.show()
# ----------------------------------------------------------------------- #

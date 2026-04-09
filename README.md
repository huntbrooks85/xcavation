
<p align="center">
    <a href="https://ibb.co/pvx3rgVs"><img src="/misc/xcavation-logo.png" width="180%"></a> <br>
</p>

<div align="center">
<p align="center">

[![SPHEREx QR2](https://img.shields.io/badge/SPHEREx-QR2-green)](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/overview_qr.html)
[![PyPI version](https://img.shields.io/pypi/v/xcavation)](https://pypi.org/project/xcavation/)
[![PyPI downloads](https://img.shields.io/pypi/dm/xcavation)](https://pypistats.org/packages/xcavation)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.19477227-blue)](https://doi.org/10.5281/zenodo.19477227)
</p>
</div>

<div align="center">
  <p id="description" > <b>Xcavation</b> is a Python package for rapid forced-photometry calculations of <a href="https://spherex.caltech.edu/">SPHEREx</a> QR2 spectral data. It uses ASQL PyVo to access <a href="https://spherex.caltech.edu/">SPHEREx</a> images through the IRSA API. Aperture photometry are supported and include multi-threading for improved deeper analysis.  </p>
</div>

<div align="center">
  <h2 style="font-size: 2em;"> 🔍 Example 🔎 </h2>
</div>

<div align="center">
  <p style="font-size: 1.2em;"> <b>
      The example notebook is provided through <a href="https://drive.google.com/file/d/1S0iLap2IoNdk4ErpdEyMyipBvJeaJ3ay/view?usp=sharing">Google Colab</a>, as it has excellent internet speed for the image API. Additionally, there is an example notebook provided on the GitHub repository. 
    </b> </p>
</div>

<div align="center">
  <h2 style="font-size: 2em;">🛠️ Installation 🛠️</h2>
</div>

<div align="center">
<pp><b> pip Installation </b><pp>
</div>
<div align="center">
</div>

1. **Download Python:** Visit [here](https://www.python.org/downloads/) to install Python 
2. **Download pip:** Visit [here](https://pip.pypa.io/en/stable/installation/) to install pip
3. **Run Install Command:** Run the command in terminal:
   ```bash
   pip install xcavation

<div align="center">
  <h2 style="font-size: 2em;">🏗️ How to Use xcavation 🏗️</h2>
</div>

<code>genspec</code> is the primary function for generating a forced-photometry spectrum, providing a full spectrum in a single line of code. User controlled options are provided to aid in problems such as: proper motion, background subtraction, masking, plotting, etc. Multi-threading for efficient processing of large datasets is provided to speed the process further. We recommend starting with <code>genspec</code> rather than calling internal functions directly. <code>genspec</code> returns a dictionary containing all relevant information from the SPHEREx calibrated images; further information on outputted data is provided below.

### How to Use <code> genspec </code>

1. After Xcavation is installed, verify the installation by running the following command: ```from xcavation.genspec import *```. If you encounter any issues during installation, please reach out to Hunter Brooks for assistance.
2. Assign the relavent variables as described below.
3. Execute the command: ```genspec(ra, dec)```. These are the minimum required parameters for Xcavation to run. You can include optional variables if needed.


### Relavent Variables For <code> genspec </code>

- **Required Variables:**
  - **ra:** Right Accension in Degrees: *float*:
     - *example:* ```131.123```

  - **dec:** Declination in Degrees: *float*:
     - *example:* ```-12.31254```

- **User Configuration:**

### Relavent Variables For <code> genspec </code>

- **Required Variables:**
  - **ra:** Right Accension in Degrees: *float*:
     - *example:* ```131.123```

  - **dec:** Declination in Degrees: *float*:
     - *example:* ```-12.31254```

- **User Configuration:**
  - **r_fwhm:** The Number Multiplied by the Image FWHM to Get The Aperture Radius: *float*
    - *example:* ```2.5```, default=```2```
  - **r_annulus_in:** The Number Multiplied by the Image FWHM to Get the Inner Annulus Radius: *float*
    - *example:* ```5```, default= ```3.5```
  - **r_annulus_out:** The Number Multiplied by the Image FWHM to Get the Outer Annulus Radius: *float*
    - *example:* ```10```, default= ```10```
  - **pmra:** Proper Motion in Right Accension (in arcsec/year): *float*
    - *example:* ```-0.981```, default=```0```
  - **pmdec:** Proper Motion in Declination (in arcsec/year): *float*
    - *example:* ```0.123```, default=```0```
  - **mjd_epoch:** Modified Julian Date of Epoch Position: *float*
    - *example:* ```57500```, default=```61000```
  - **mjd_query:** Modified Julian Date of Query Position: *float*
    -  *example:* ```61500```, default=```Time(datetime.now()).mjd```
  - **mjd:** Modified Julian Date of inputed R.A. and Decl. from Above: *float*
    - *example:* ```57170```, default=```61000```
  - **save_data:** Save Q.A. Finder Chart, Spectrum, and ASCII Table if True: *boolean*
    - *example:* ```True```, default=```False```
  - **output_path:** Saved Path (not including extention) for Q.A. Finder Chart, Spectrum, and ASCII Table if ```save_data``` is True: *string*
    - *example:* ```/content/drive/MyDrive/Research/example_object```, default=```xcavation_save```
  - **threads:** Number of Threads for Multi-Threading: *int*
    - *example:* ```2```, default=```1```
  - **enable_print:** Enable ```xcavation.genspec``` Progress Printing: *boolean*
    - *example:* ```False```, default=```True```
  - **ram_download:** Enable ```xcavation.aperture``` Downloading Image Data on RAM Instead of Local Harddrive (best for computers with small amounts of storage): *boolean*
    - *example:* ```True```, default=```False```
  - **retry_count:** The Number of HTML 504 Retries Before xcavation.genspec Crashes: *int*
    - *example:* ```1```, default=```10```
  - **clean_type:** The Type of Flagged Pixel Cleaning Done to the Calibrated image: *string*
    - *options:* ```none```, ```mask```, ```median_mask```, ```interp_mask```, ```lacosmic```
    - *example:* ```mask```, default=```none```
  - **bad_bits:** Masks Out Flags Seen in [SPHEREx Documentation](https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/#id-9-explore-the-second-extension-flags) and Performs Laplace Transforms Using the ```astroscrappy.detect_cosmics``` to Replace These Pixels: *list*
    - *options:* Go to Section 6 of [HERE](https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/#id-7-visualize-a-spherex-spectral-image-mef-using-the-firefly-python-client) for each flag type.
    - *example:* ```[0, 1, 5, 6, 10, 11]```, default=```[0, 1, 10, 11]```
  - **background_type:** How the Background Field is Calculated: *string*
    - *options:* ```mean```, ```median```, ```mode```
    - *example:* ```mode```, default=```mean```
  - **cutout_size:** The Size of the Cutout in Arcsec: *int*
    - *example:* ```500```, default=```150```
  - **zodi_subtract:** Whether or not the Zodiacal light is subtracted: *bool*
    - *example:* ```False```, default=```True```
  - **sigclip_sigma:** Background Annulus [Astropy.SigmaClip](https://docs.astropy.org/en/stable/api/astropy.stats.SigmaClip.html) Sigma Variable: *float*
    - *example:* ```3```, default=```5```
  - **sigclip_maxiters:** Background Annulus [Astropy.SigmaClip](https://docs.astropy.org/en/stable/api/astropy.stats.SigmaClip.html) Max Iterations Variable: *float*
    - *example:* ```10```, default=```5```

### Output Variables:
1.  <code>wavelength</code>: Wavelength of a queried position: $\mu m$: <code> float </code>
2.  <code>delta_lambda</code>: The resolving power solution ($\Delta \lambda = \frac{\lambda}{R}$) from Table 2 of [Akeson et al. (2025)](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHEREx_Expsupp_QR_v1.0.pdf): $\mu m$: <code> float </code>
3.  <code>flux</code>: Flux from the given aperture photometry: $\mu Jy$: <code> float </code>
4.  <code>flux_err</code>: Flux error from the given aperture photometry: $\mu Jy$: <code> float </code>
5. <code>flag_count</code>: Dictionary of every bit measured in the aperture
6.  <code>SNR</code>: Signal-to-Noise Ratio ($\frac{F}{\sigma_{F}}$): <code> float </code>
7.  <code>flux_cutout</code>: The cutout flux image in $\mu Jy$: <code> numpy array</code>
8.  <code>flag_cutout</code>: The cutout flag image in bits: <code> numpy array</code>
9.  <code>aperture</code>: The aperture array in cutout size in booleans: <code> numpy array</code>
10.  <code>annulus</code>: The annulus array in cutout size in booleans: <code> numpy array</code>
11.  <code>x_loc</code>: The x_loc of the centeroid in the flux cutout: <code> float </code>
12.  <code>y_loc</code>: The y_loc of the centeroid in the flux cutout: <code> float </code>
13. <code> ap_radius</code>: The radius of the aperture: 2.5 $\times$ FWHM: $arcsec$: <code> float </code>
14. <code> inner_annulus</code>: The inner radius of the annulus: $arcsec$: <code> float </code>
15. <code> outer_annulus</code>: The outer radius of the annulus: $arcsec$: <code> float </code>
16. <code> url</code>: The IPAC IRSA SPHEREx API url: <code> string </code>
17. <code> mjd</code>: The Modified Julian Date of the given SPHEREx image: <code> float </code>

### Example
```
from xcavation.genspec import *

# Required Variables
ra = 12.7955492 
dec = -15.7382839 

# User Configuration
user_config = genspec_profile(pmra = 0, pmdec = 0,
                              mjd_epoch = 61000, mjd_query = 61500,
                              r_fwhm = 2, r_annulus_in = 5, r_annulus_out = 10,
                              save_data = True, output_path = 'test', 
                              threads = 10,
                              enable_print = True, ram_download = False, 
                              retry_count = 10,
                              clean_type = 'mask', bad_bits = [0, 1, 10, 11],
                              background_type = 'mean',
                              cutout_size = 150,
                              zodi_subtract = True,
                              sigclip_sigma = 5, sigclip_maxiters = 5)

# Call function with test input parameters
df = genspec(ra, dec, config = user_config)
```

<div align="center">
  <h2 style="font-size: 2em;">📞 Support Team 📞</h2>
</div>

- **Mr. Hunter Brooks**
  - hbrooks8 (at) rockets.utoledo.edu
- **Dr. Michael Cushing**

<div align="center">
  <h2 style="font-size: 2em;">📖 Acknowledgments 📖</h2>
</div>

1. If you intend to publish any calculations done by xcavation, please reference Brooks et al. (in prep.).

2. Please reference the relevant [SPHEREx citations](https://spherex.caltech.edu/page/publications).


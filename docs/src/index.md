# UVITTools.jl
*A package for AstroSat/UVIT grating spectroscopy and aperture photometry.*

## Package Features
- Extract flux calibrated UVIT FUV/NUV grating spectra.
- Generate source and background PHA spectral data compatible with XSPEC/Sherpa.
- Perform aperture photometry of FUV/NUV sources in any of the broadband filters.
- Saturation correction for point sources.

## UVIT Grating Spectroscopy
The Indian multi-wavelength astronomy satellite carries four co-aligned instruments. One of them, the Ultra-Violet Imaging Telescope ([UVIT](http://astrosat-ssc.iucaa.in/uvitData); Tandon et al. ([2017](https://iopscience.iop.org/article/10.3847/1538-3881/aa8451),[2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3))) carries two gratings in the FUV channel and a single grating in the NUV channel. These gratings are useful for low resolution, slitless spectroscopy in the far and near UV bands of a variety of cosmic sources such as hot stars, interacting binaries, active galactic nuclei, etc. The calibration of the gratings are described in [Dewangan (2021)](https://link.springer.com/article/10.1007/s12036-021-09691-w). This package includes tools for the UVIT grating spectroscopy based on updated calibration following the method described in [Dewangan (2021)](https://link.springer.com/article/10.1007/s12036-021-09691-w). The calibration files are included in this package in the `caldata` directory and are directly used by the tools.

The `UVITTools.jl` package has been designed to use the processed UVIT data generated with the [CCDLAB package](https://github.com/user29A/CCDLAB) ([Postma & Leahy 2017](https://ui.adsabs.harvard.edu/abs/2017PASP..129k5002P/abstract)), hence the users will first need to process the level1 UVIT data with the CCDLAB. The `UVITTools.jl` package provides functions for each of the steps in generating the calibrated spectra (i. e., $f_\lambda$ Vs $\lambda$) as well as a single tool to generate the final calibrated spectra. The package also provides functions to generate source and background PHA spectral files and the associated instrument response file (response matrix + effective area) for each of the gratings. These PHA spectral data can directly be used in X-ray spectral fitting packages such as [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) and [SHERPA](https://cxc.cfa.harvard.edu/sherpa/).  The grating PHA spectral data are treated in the same as the X-ray spectral data, and are useful for simultaneous UV/X-ray spectral modelling.

### Extraction of 1d count spectra

The [CCDLAB pipeline](https://ui.adsabs.harvard.edu/abs/2017PASP..129k5002P/abstract) has the option of generating images from the grating observations either in the detector coordinates left unrotated or rotated to align the dispersion axis with the X-axis.  The 1d spectral extraction functions for the FUV gratings (`fuv_grating1_count_spec` and `fuv_grating2_count_spec`) include two parameters `disp_aligned_to_xaxis` and `angle_xaxis_disp_deg` ($\theta$) that can be used to specify if the grating images are rotated or not (for more details, please see the function documentation). The default values of the angles should work fine, otherwise users can measure the angle of the dispersion direction with respect to the X-axis and specify in the command line. 

The 1d spectral extraction tools construct a linear relation $y = mx +c$   ($m = \tan{\theta}$) that defines the  spectral trace. The coordinates along the dispersion direction are defined relative to the zero order position provided by the use in the form of a ds9 region file. The 1d count spectrum is then extracted by integrating the counts along the spatial direction within a cross dispersion width provided by the user.

### Wavelength Calibration

The relative pixel coordinates along the dispersion direction are converted using linear dispersion relations 
$$
\lambda = c_0 + c_1 \times pixel,
$$
where the coefficients $c_0$ and $c_1$ are derived from wavelength calibration based on observations of the planetary nebula NGC40. The UVITToools.jl functions  ` fuv_grating1_pixel2lamA`, `fuv_grating1_wavelength_calib` for FUV-Grating1 and similar functions for other gratings use updated wavelength calibrations compared to those provided in [Dewangan (2021)](https://link.springer.com/article/10.1007/s12036-021-09691-w).

### Flux calibration

The count rates at each wavelength are converted to flux using the effective area derived from the grating observations of spectrophotometric standards. The effective area files for three grating and calibrated orders are provided in the `caldata` directory, these files are used by the flux calibration functions `fuv_grating1_flux_calib`, `fuv_grating2_flux_calib` and `nuv_grating1_m1_flux_calib` (see the function documentation for more details). 

### Generating wavelengh and flux calibrated spectra in one step

The `UVITTools.jl` package provides functions for direct generation of background-corrected, wavelength and flux calibrated spectra in one shot without going through the individual steps. Users can use the functions `fuv_grating1_fluxed_spec`, `fuv_grating2_fluxed_spec`, and `nuv_grating_m1_fluxed_spec` for this purpose (see function documentation for more details).

### XSPEC/Sherpa-compatible PHA spectral data

Working with data from different instruments onboard AstroSat require tools and techniques to facilitate joint analysis of multi-wavelength data. In particular, the broadband spectral coverage of AstroSat, from near UV to hard X-rays, requires tools for simultaneous fitting of spectral models to the multi-wavelength data.  For this purpose, `UVITTools.jl` package provides functions to generate grating PHA spectral data for the source and background; it also includes the associated response files in the caldata directory. The functions available for the extraction of  PHA spectral data are  `fuv_grating1_phafile`, `fuv_grating2_phafile` and `nuv_grating_phafile`. These functions will automatically use the appropriate response files from the caldata directory. See the function documentation for more details.

### Redistribution matrix and ancillary response files

The `UVITTools.jl` package includes updated response files generated using the method described in [Dewangan (2021)](https://link.springer.com/article/10.1007/s12036-021-09691-w). The response files are available in the caldata directory. Users are not required to download these files as the PHA spectral functions will automatically download the relevant response file. Users can utilise the response files to simulate grating spectal data in the usual way with XSPEC or SHERPA.

## UVIT Aperture Photometry

The UVIT is primarily a high resolution imaging instrument in a number of broadband filters (see <https://uvit.iiap.res.in/Instrument/Filters>). The photometric calibration of all the broadband filters are performed in Tandon et al. ([2017](https://iopscience.iop.org/article/10.3847/1538-3881/aa8451), [2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3)). The `UVITTools.jl` provides  functions for aperture photometry (see `uvit_aphot`) using simple apertures such as circular, elliptical, annular regions. The aperture photometry tasks output count rate, flux and magnitude based on the photometric calibration.

### Saturation correction

The FUV/NUV channels operate in the photon counting mode. If multiple photon events fall within $3\times 3$ pixels in a frame, then these photon events can be detected individually but recorded as arising due to a single photon. Such saturation of photon events cannot be avoided unless the photon rate per frame is much less than 1. Tandon et al. ([2017](https://iopscience.iop.org/article/10.3847/1538-3881/aa8451),[2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3)) have devised an algorithm to correct for the saturation.  This algorithm is implemented in the `UVITTools.jl` function `uvit_saturation_corr`, which is called by the aperture photometry task `uvit_aphot` if the option for saturation correction is provided. The algorithm for the saturation correction is applicable to point sources only, and users need to use a circular extraction region of radius about 12 arcsec (see [Tandon et al. 2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3) for details).

### XSPEC/Sherpa compatible single channel spectral data and responses

The `UVITTools.jl` package also includes function to extract single channel PHA spectral data from FUV/NUV images in any of the broadband filters. This function also uses appropriate response files available in the caldata directory, the users do not need to select or separately download the response files.  These response files can also be used to simulate count rate in any of the filters using XSPEC or SHERPA. All the filter response files were created from the updated calibration provided in [Tandon et al. (2020)](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3). The single channel spectra can be used along with X-ray spectral data for broadband spectral fitting in XSPEC or SHERPA.



## Function Documentation
#### FUV-Grating1
```@docs
fuv_grating1_count_spec
```

```@docs
fuv_grating1_net_countrate_spec
```

```@docs
fuv_grating1_pixel2lamA
```

```@docs
fuv_grating1_wavelength_calib
```

```@docs
fuv_grating1_flux_calib
```

```@docs
fuv_grating1_fluxed_spec
```

```@docs
fuv_grating1_phafile
```

```@docs
fuv_grating1_ea
```

#### FUV-Grating2

```@docs
fuv_grating2_count_spec
```

```@docs
fuv_grating2_net_countrate_spec
```

```@docs
fuv_grating2_pixel2lamA
```
```@docs
fuv_grating2_wavelength_calib
```

```@docs
fuv_grating2_flux_calib
```

```@docs
fuv_grating2_fluxed_spec
```

```@docs
fuv_grating2_phafile
```

```@docs
fuv_grating2_ea
```

#### NUV-Grating

```@docs
nuv_grating_m1_count_spec
```

```@docs
nuv_grating_m1_net_countrate_spec
```

```@docs
nuv_grating_m1_wavelength_calib
```
```@docs
nuv_grating_m1_flux_calib
```

```@docs
nuv_grating_m1_fluxed_spec
```

```@docs
nuv_grating_phafile
```

### Aperture Photometry

```@docs
read_ds9reg
```

```@docs
uvit_aphot
```

```@docs
uvit_countrate2flux
```

```@docs
uvit_filter2pha
```

```@docs
uvit_saturation_corr
```

```@docs
uvit_zp_uc
```

### Miscellaneous

```@docs
xycen_from_ds9reg
```

```@docs
write_uvit_grating_phafile
```
```@docs
lambdaA2keV
```

```@docs
lambdaA2ergs
```

```@docs
flux_density_photons2cgs
```

```@docs
flux_density_cgs2photons
```
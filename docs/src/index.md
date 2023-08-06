# UVITTools.jl
*A package for AstroSat/UVIT grating spectroscopy and aperture photometry.*

## Package Features
- Extract flux calibrated UVIT FUV/NUV grating spectra.
- Generate  source and background FUV/NUV grating spectral PHA data compatible with XSPEC/Sherpa.
- Perform aperture photometry of FUV/NUV sources in any of the broadband filters.
- Generate single channel PHA source and background spectral data for any of the  FUV/NUV broadband filters.
- Saturation correction for point sources.

## UVIT Grating Spectroscopy
The Indian multi-wavelength astronomy satellite carries four co-aligned instruments. One of them, the Ultra-Violet Imaging Telescope ([UVIT](http://astrosat-ssc.iucaa.in/uvitData); Tandon et al. ([2017](https://iopscience.iop.org/article/10.3847/1538-3881/aa8451),[2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3))) is equipped with two gratings in the FUV channel and a single grating in the NUV channel in addition to a number of broadband filters. The slitless gratings are useful for low resolution far and near UV spectroscopy  of a variety of cosmic sources such as hot stars, interacting binaries, active galactic nuclei, etc. The `UVITTools.jl` package includes tools for the UVIT grating spectroscopy based on updated calibration following the method described in [Dewangan (2021)](https://link.springer.com/article/10.1007/s12036-021-09691-w). The calibration files are included in this package in the `caldata` directory and are directly used by the tools within the package.

The `UVITTools.jl` package has been designed to use the processed UVIT data generated with the [CCDLAB pipeline](https://github.com/user29A/CCDLAB) ([Postma & Leahy 2017](https://ui.adsabs.harvard.edu/abs/2017PASP..129k5002P/abstract)), hence the users first need to process the level-1 UVIT data with the CCDLAB. The `UVITTools.jl` package provides functions for each of the steps in generating the calibrated spectra (i. e., $f_\lambda$ Vs $\lambda$) as well as single tools to generate the final calibrated spectra for different gratings and orders. The package also provides functions to generate source and background PHA spectral files and the associated instrument response file (response matrix + effective area) for each of the gratings.  The grating PHA spectral data are treated in the same as the X-ray spectral data, and are useful for simultaneous UV/X-ray spectral modelling using  [XSPEC](https://heasarc.gsfc.nasa.gov/xanadu/xspec/) and [SHERPA](https://cxc.cfa.harvard.edu/sherpa/).

### Extraction of 1d count spectra from 2d images

The [CCDLAB pipeline](https://ui.adsabs.harvard.edu/abs/2017PASP..129k5002P/abstract) has the option of generating images from the grating observations either in the detector coordinates left unrotated or rotated to align the dispersion axis with the X-axis.  The 1d spectral extraction functions for the FUV gratings (`fuv_grating1_count_spec` and `fuv_grating2_count_spec`) include two parameters `disp_aligned_to_xaxis` and `angle_xaxis_disp_deg` ($\theta$) that can be used to specify if the grating images are rotated and the angle between the X-axis and the dispersion direction towards the negative orders (for more details, please see the function documentation). The default values of the angles should work fine, otherwise users can measure the angle of the dispersion direction towards negative orders with respect to the X-axis and specify in the command line. 

The 1d spectral extraction tools construct a linear relation $y = mx +c$  with  ($m = \tan{\theta}$) that defines the  spectral trace. The coordinates along the dispersion direction are defined relative to the zero order position provided by the user in the form of a [SAOImage/DS9](https://sites.google.com/cfa.harvard.edu/saoimageds9/download)  region file. The 1d count spectrum is then extracted by integrating the counts along the spatial direction within a cross dispersion width provided by the user.

The following steps can be followed to generate 1d count spectrum of a source from a grating image generated with the CCDLAB pipeline.

-  Create region file for a source of interest using DS9 in the detector coordinates. Use a circular region centred at the zero order position, the size of the circle is immaterial.

- To extract 1d count spectrum of NGC40 in the -2 order of FUV-Grating1 using a 40 pixel width along the spatial direction, issue the following command.  

  ```julia
  julia> (rel_pixels, count_rate) = fuv_grating1_count_spec("NGC40_FUV_Grating1_IMAGE.fits","src.reg", order=-2, cross_disp_width_pixels = 40, angle_xaxis_disp_deg=0.0, rate=true)
  ([-629, -628, -627, -626, -625, -624, -623, -622, -621, -620  …  -422, -421, -420, -419, -418, -417, -416, -415, -414, -413], Measurements.Measurement{Float64}[0.0068 ± 0.0022, 0.0091 ± 0.0026, 0.0068 ± 0.0022, 0.0105 ± 0.0028, 0.009 ± 0.0026, 0.0083 ± 0.0025, 0.0105 ± 0.0028, 0.0105 ± 0.0028, 0.0067 ± 0.0022, 0.0075 ± 0.0024  …  0.006 ± 0.0021, 0.0067 ± 0.0022, 0.006 ± 0.0021, 0.0045 ± 0.0018, 0.006 ± 0.0021, 0.0067 ± 0.0022, 0.0037 ± 0.0017, 0.0075 ± 0.0023, 0.0074 ± 0.0023, 0.0075 ± 0.0023])
  
  ```

  In addition to the screen output, a DS9 region file used for the grating order and the cpunt spectrum are written to files `src_fuv_grating1_m2.reg` and `src_fuv_grating1_m2_count_spec.dat`. The region files can be loaded in DS9 to cross-check the used grating order. 

  Similarly, one can use the function `fuv_grating1_net_countrate_spec` to extract the background-corrected net count rate spectrum of a source by additionally using a DS9 region background file.

  ```julia
  julia> (rel_pixels, net_count_rate) = fuv_grating1_net_countrate_spec("NGC40_FUV_Grating1_IMAGE.fits","src.reg", "bgd.reg",order=-2, cross_disp_width_pixels = 40, angle_xaxis_disp_deg=0.0)
  ([-629, -628, -627, -626, -625, -624, -623, -622, -621, -620  …  -422, -421, -420, -419, -418, -417, -416, -415, -414, -413], Measurements.Measurement{Float64}[0.0053 ± 0.0025, 0.0068 ± 0.0029, 0.0045 ± 0.0026, 0.006 ± 0.0033, 0.0075 ± 0.0028, 0.006 ± 0.0028, 0.0083 ± 0.0031, 0.0098 ± 0.0029, 0.003 ± 0.0028, 0.0045 ± 0.0028  …  0.003 ± 0.0026, 0.003 ± 0.0028, 0.0045 ± 0.0023, 0.003 ± 0.0021, 0.0045 ± 0.0023, 0.0022 ± 0.0029, -0.00073 ± 0.0025, 0.006 ± 0.0026, 0.0045 ± 0.0028, 0.0052 ± 0.0027])
  ```

  

### Wavelength Calibration

The relative pixel coordinates along the dispersion direction are converted using linear dispersion relations 
$$
\lambda = c_0 + c_1 \times pixel,
$$
where the coefficients $c_0$ and $c_1$ are derived from wavelength calibration based on observations of the planetary nebula NGC40. The UVITToools.jl functions  ` fuv_grating1_pixel2lamA`, `fuv_grating1_wavelength_calib` for FUV-Grating1 and similar functions for other gratings use updated wavelength calibrations compared to those provided in [Dewangan (2021)](https://link.springer.com/article/10.1007/s12036-021-09691-w).

One can convert the relative pixel coordinates to wavelength using the command

```julia
julia> fuv_grating1_pixel2lamA(-500, order=-2)
1440.3065041673246

julia> fuv_grating1_pixel2lamA(-250, order=-1)
1440.2665
```

and perform wavelength calibration of a count spectrum as follows.

```julia
julia> (rel_pixels, net_count_rate) = fuv_grating1_net_countrate_spec("NGC40_FUV_Grating1_IMAGE.fits","src.reg", "bgd.reg",order=-2, cross_disp_width_pixels = 40, angle_xaxis_disp_deg=0.0)
([-629, -628, -627, -626, -625, -624, -623, -622, -621, -620  …  -422, -421, -420, -419, -418, -417, -416, -415, -414, -413], Measurements.Measurement{Float64}[0.0053 ± 0.0025, 0.0068 ± 0.0029, 0.0045 ± 0.0026, 0.006 ± 0.0033, 0.0075 ± 0.0028, 0.006 ± 0.0028, 0.0083 ± 0.0031, 0.0098 ± 0.0029, 0.003 ± 0.0028, 0.0045 ± 0.0028  …  0.003 ± 0.0026, 0.003 ± 0.0028, 0.0045 ± 0.0023, 0.003 ± 0.0021, 0.0045 ± 0.0023, 0.0022 ± 0.0029, -0.00073 ± 0.0025, 0.006 ± 0.0026, 0.0045 ± 0.0028, 0.0052 ± 0.0027])
```

### Flux calibration

The count rates at each wavelength are converted to flux using the effective area derived from the grating observations of spectrophotometric standards. The effective area files for three grating and calibrated orders are provided in the `caldata` directory, these files are used by the flux calibration functions `fuv_grating1_flux_calib`, `fuv_grating2_flux_calib` and `nuv_grating1_m1_flux_calib` (see the function documentation for more details). These functions assume that the extraction width along the spatial direction include all the source flux as the aperture correction is not built-in.

Flux calibration of a wavelength-calibrated count rate spectrum for a -2 order FUV-Grating1 can be performed as follows.

```julia
julia> (λ, f_λ)=fuv_grating1_flux_calib(λ, net_count_rate, order=-2)
([1800.1595912456048, 1797.3700324310446, 1794.5804736164841, 1791.790914801924, 1789.0013559873637, 1786.2117971728032, 1783.422238358243, 1780.6326795436828, 1777.8431207291223, 1775.0535619145621  …  1222.7209166316204, 1219.9313578170602, 1217.1417990024997, 1214.3522401879395, 1211.5626813733793, 1208.7731225588188, 1205.9835637442586, 1203.1940049296982, 1200.404446115138, 1197.6148873005777], Measurements.Measurement{Float64}[6.5e-14 ± 3.0e-14, 8.3e-14 ± 3.5e-14, 5.5e-14 ± 3.1e-14, 7.1e-14 ± 4.0e-14, 8.8e-14 ± 3.3e-14, 6.9e-14 ± 3.2e-14, 9.2e-14 ± 3.4e-14, 1.06e-13 ± 3.1e-14, 3.1e-14 ± 2.9e-14, 4.6e-14 ± 2.8e-14  …  5.4e-14 ± 4.6e-14, 6.3e-14 ± 5.9e-14, 1.11e-13 ± 5.8e-14, 8.7e-14 ± 6.1e-14, 1.49e-13 ± 7.8e-14, 8.3e-14 ± 1.1e-13, -2.9e-14 ± 9.6e-14, 2.3e-13 ± 9.9e-14, 1.59e-13 ± 9.9e-14, 1.64e-13 ± 8.4e-14])
```

### Generating wavelength and flux calibrated spectra in one step

The `UVITTools.jl` package provides functions for direct generation of background-corrected, wavelength and flux calibrated spectra in one shot without going through the individual steps. Users can use the functions `fuv_grating1_fluxed_spec`, `fuv_grating2_fluxed_spec`, and `nuv_grating_m1_fluxed_spec` for this purpose (see function documentation for more details). For example, to extract background-correced, wavelength and flux-calibrated spectrum in the -2 order of FUV-Grating1, one can use the following command.

```julia
julia> (lam, flam, err_flam) = fuv_grating1_fluxed_spec("ngc40", "NGC40_FUV_Grating1_IMAGE.fits","src.reg","bgd.reg", order=-2, cross_disp_width_pixels = 40, angle_xaxis_disp_deg=0.0)
-----------------------------------
target=ngc40
UVIT channel=FUV
Grating=Grating1
OBS_ID=C07_015T01_9000005256
order=-2
Exposure time=1353.258 seconds
---------------------------------------
Wrote spectral ascii file: ngc40_C07_015T01_9000005256_FUV_Grating1m2_crossdisp40pix_xax_disp_0.0deg_spec.dat
```

Similar functions for FUV-Grating2 and NUV-Grating are also available.

### XSPEC/Sherpa-compatible PHA spectral data

Working with data from different instruments onboard AstroSat require tools and techniques to facilitate joint analysis of multi-wavelength data. In particular, the broadband spectral coverage of AstroSat, from near UV to hard X-rays, requires tools for simultaneous fitting of spectral models to the multi-wavelength data.  For this purpose, `UVITTools.jl` package provides functions to generate grating PHA spectral data for the source and background; it also includes the associated response files in the `caldata` directory. The functions available for the extraction of  PHA spectral data are  `fuv_grating1_phafile`, `fuv_grating2_phafile` and `nuv_grating_phafile`. These functions will automatically use the appropriate response files from the caldata directory. See the function documentation for more details. To generate PHA spectral dataset (source and background spectral data, and response file) for a source from -2 order of FUV-Grating1 image, run the following command. The background and response files are automatically written in the header of the spurce spectral file.

```julia
julia> fuv_grating1_phafile("ngc40", "NGC40_FUV_Grating1_IMAGE.fits","src.reg","bgd.reg", order=-2, cross_disp_width_pixels = 40, angle_xaxis_disp_deg=0.0)
-----------------------------------
target=ngc40
UVIT channel=FUV
Grating=Grating1
OBS_ID=C07_015T01_9000005256
order=-2
Exposure time=1353.258 seconds
---------------------------------------
Writing source and background PHA files..
Using respfile/home/gulabd/work/julia_dev/UVITTools/caldata/fuv_grating1_m2_12nov22.rmf
Using /home/gulabd/work/julia_dev/UVITTools/caldata/fuv_grating1_m2_12nov22.rmf
("ngc40_C07_015T01_9000005256_FUV_Grating1_m2_crossdisp40pix_xax_disp_0.0deg_src.pha", "ngc40_C07_015T01_9000005256_FUV_Grating1_m2_crossdisp40pix_xax_disp_0.0deg_bgd.pha")
```



### Redistribution matrix and ancillary response files

The `UVITTools.jl` package includes updated response files generated using the method described in [Dewangan (2021)](https://link.springer.com/article/10.1007/s12036-021-09691-w). The response files are available in the caldata directory. Users are not required to download these files as the PHA spectral functions will automatically download the relevant response file. Users can utilise the response files to simulate grating spectal data in the usual way with XSPEC or SHERPA.

## UVIT Aperture Photometry

The UVIT is primarily a high resolution imaging instrument in a number of broadband filters (see <https://uvit.iiap.res.in/Instrument/Filters>). The photometric calibration of all the broadband filters are performed in Tandon et al. ([2017](https://iopscience.iop.org/article/10.3847/1538-3881/aa8451), [2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3)). The `UVITTools.jl` provides  functions for aperture photometry (see `uvit_aphot`) using simple apertures such as circular, elliptical, annular regions. The aperture photometry task output count rate, flux and magnitude based on the photometric calibration. Currently, `uvit_aphot` does not account for aperture correction, hence users will need to use large enough extraction region (a circular region with radius at least 25 pixels)  to include almost all counts from point sources.

To perform aperture photometry including saturation correction (see below) of a source in FUV/BaF2 image, one can run the following command. 

```julia
julia> uvit_aphot("NGC4593_FUV_BaF2___MASTER.fits","src.reg","bgd.reg",satu_corr=true)
---------count rates, flux and magnitude------------
Detector: FUV
Filter: BaF2
Exposure time: 17129.3 seconds
source+background rate = 4.2 ± 0.016 counts/s
background rate = 0.051 ± 0.0 counts/s
net source count rate = 4.149 ± 0.016 counts/s
Saturation corrected net source count rate= 4.452 ± 0.016 counts/s
Saturation corrected f_λ [BaF2]= 1.591e-14 ± 1.6e-16 ergs/cm2/s/A
Saturation corrected  magnitude[BaF2] (AB system) = 16.15 ± 0.011
------------------------------------------------------
Mean BJD: 2.457588076583e6
------------------------------------------------------
(2.457588076583e6, 4.451869565968372, 0.016218003910273074)
```

### Saturation correction

The FUV/NUV channels operate in the photon counting mode. If multiple photon events fall within $3\times 3$ pixels in a frame, then these photon events can be detected individually but recorded as arising due to a single photon. Such saturation of photon events cannot be avoided unless the photon rate per frame is much less than 1. Tandon et al. ([2017](https://iopscience.iop.org/article/10.3847/1538-3881/aa8451),[2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3)) have devised an algorithm to correct for the saturation.  This algorithm is implemented in the `UVITTools.jl` function `uvit_saturation_corr`, which is called by the aperture photometry task `uvit_aphot` if the option for saturation correction is provided. The algorithm for the saturation correction is applicable to point sources only, and users need to use a circular extraction region of radius about 12 arcsec (see [Tandon et al. 2020](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3) for details).

### XSPEC/Sherpa compatible single channel spectral data and responses

The `UVITTools.jl` package also includes function to extract single channel PHA spectral data from FUV/NUV images in any of the broadband filters. This function also uses appropriate response files available in the caldata directory, the users do not need to select or separately download the response files.  These response files can also be used to simulate count rate in any of the filters using XSPEC or SHERPA. All the filter response files were created from the updated calibration provided in [Tandon et al. (2020)](https://iopscience.iop.org/article/10.3847/1538-3881/ab72a3). The single channel spectra can be used along with X-ray spectral data for broadband spectral fitting in XSPEC or SHERPA.

Single channel source and background count spectra (and the corresponding response) can be generated as follows.

```julia
julia> uvit_filter2pha("ngc4593","NGC4593_FUV_BaF2___MASTER.fits","src.reg","bgd.reg")
---------count rates------------
Detector: FUV
Filter: F2
Exposure time: 17129.3 seconds
source+background rate = 4.2+/-0.0157 counts/s
background rate = 0.051+/-0.0004 counts/s
net source count rate = 4.149+/-0.0157 counts/s
Saturation corrected net source count rate= 4.452 ± 0.016 counts/s
--------------------------------
---Writing PHA file-----
Dict{String, Vector}("GROUPING" => Int16[1], "COUNTS" => [77132.30772409608], "CHANNEL" => Int32[1], "QUALITY" => Int16[0])
Dict{String, Vector}("GROUPING" => Int16[1], "COUNTS" => [874.8137822323066], "CHANNEL" => Int32[1], "QUALITY" => Int16[0])
("ngc4593_G05_219T01_9000000_FUV_F2_spec_src.pha", "ngc4593_G05_219T01_9000000_FUV_F2_spec_bgd.pha")
```

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
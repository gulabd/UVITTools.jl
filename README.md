# UVITTools

A tool to analyze data from the [AstroSat](http://astrosat-ssc.iucaa.in/)/[UVIT](https://uvit.iiap.res.in/) payload.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gulabd.github.io/UVITTools.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gulabd.github.io/UVITTools.jl/dev)
[![Build Status](https://github.com/gulabd/UVITTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gulabd/UVITTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/gulabd/UVITTools.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gulabd/UVITTools.jl)

## Installation

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/gulabd/UVITTools.jl")
```

## Usage
UVITTools can perform photometery, spectroscopy of astronomical sources from UVIT FUV/NUV images acquired with the AstroSat/UVIT payload.  The package is meant for images generated with the [CCDLAB package](https://iopscience.iop.org/article/10.1088/1538-3873/aa8800). The tools are not compatible with the data generated from the  [UVIT level2 pipeline software](http://astrosat-ssc.iucaa.in/uvitData).
```julia
julia> using UVITTools
```
To find the various tools available, type

```julia
julia> ?UVITTools
```
List of tool to analyze data from the AstroSat/UVIT payload.


## Grating Spectroscopy
- `xycen_from_ds9reg`
- `write_uvit_grating_phafile`
- `fuv_grating1_count_spec `
- `fuv_grating1_net_countrate_spec`
- `fuv_grating1_pixel2lamA ` 
- `fuv_grating1_wavelength_calib `
- `fuv_grating1_flux_calib `
- `fuv_grating1_fluxed_spec `
- `fuv_grating1_phafile `
- `fuv_grating1_ea`
- `fuv_grating2_count_spec `
- `fuv_grating2_net_countrate_spec `
- `fuv_grating2_pixel2lamA `
- `fuv_grating2_wavelength_calib `
- `fuv_grating2_flux_calib`
- `fuv_grating2_fluxed_spec`
- `fuv_grating2_phafile `
- `fuv_grating2_ea`
- `nuv_grating_m1_count_spec `
- `nuv_grating_m1_net_countrate_spec `
- `nuv_grating_m1_wavelength_calib `
- `nuv_grating_m1_flux_calib `
- `nuv_grating_m1_fluxed_spec`
- `nuv_grating_phafile`


##  UVIT Photometry 
- `read_ds9reg`
- `uvit_aphot `
- `uvit_countrate2flux `
- `uvit_filter2pha`
- `uvit_saturation_corr`
- `uvit_zp_uc`

## Unit conversions
- `lambdaA2keV`
- `lambdaA2ergs`
- `flux_density_photons2cgs`
- `flux_density_cgs2photons`



To see the help for any tool (e.g., `nuv_grating_phafile`), type

```julia
julia>?
help?> nuv_grating_phafile
```

### Photometry
UVITTools can perform photometry of individual sources present in FUV or NUV image in any filter. Users will need to create source and background region file using DS9 and input to the photometry tool. 

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

### XSPEC-style PHA spectral file for broadband filters

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
### Grating Spectroscopy - fluxed spectra

  To extract fluxed spectrum in the -2 order of FUV-Grating1, type 
```julia
 julia> (lam, flam, err_flam) = fuv_grating1_fluxed_spec("ngc40", "NGC40_FUV_Grating1_IMAGE.fits","src.reg","bgd.reg", order=-2, cross_disp_width_pixels = 40, angle_xaxis_disp_deg=0.0)
-----------------------------------
target=ngc40
UVIT channel=FUV
Grating=Grating1
OBS_ID=C02_010T01_9000000
order=-2
Exposure time=1194.007 seconds
---------------------------------------
Wrote spectral ascii file: ngc40_C02_010T01_9000000_FUV_Grating1m2_crossdisp_50pixels_spec.dat

([1196.1544204935722, 1198.9454527675634, 1201.7364850415547, 1204.527517315546, 1207.3185495895373, 1210.1095818635285, 1212.9006141375198, 1215.691646411511, 1218.4826786855024, 1221.2737109594932  …  1773.898101209758, 1776.6891334837492, 1779.4801657577405, 1782.2711980317317, 1785.062230305723, 1787.8532625797143, 1790.6442948537056, 1793.4353271276968, 1796.2263594016881, 1799.0173916756794], [5.045466629574188e-13, 4.3150185620641407e-13, 3.20002438057737e-13, 2.1473337940718712e-13, 3.1256136935946915e-13, 3.309620601584206e-13, 2.330335115352435e-13, 2.8446859045530054e-13, 2.969928077009422e-13, 1.8559317187843672e-13  …  4.78007919733574e-13, 4.1855639102711703e-13, 6.104105792024475e-13, 5.793213715591102e-13, 5.984233168019304e-13, 6.895818454682478e-13, 6.118620245824862e-13, 7.615785936600397e-13, 8.907679878495615e-13, 7.61419751498956e-13], [1.0267947190347392e-13, 8.637264104461592e-14, 6.972295297643998e-14, 5.840516168789386e-14, 5.6833639705142054e-14, 5.42522362473877e-14, 4.587473745925019e-14, 4.4072190863771767e-14, 4.2681889008873683e-14, 3.599499296873876e-14  …  5.96764312139293e-14, 6.233276080702569e-14, 7.279889935967496e-14, 7.101598452391445e-14, 7.784002780312587e-14, 8.738266189464572e-14, 8.352764933553134e-14, 1.0110334674900172e-13, 1.1547586393436977e-13, 1.0992802034895543e-13])
```

To plot the spectrum,

julia> plot(lam,flam,yerr=err_flam)



To extract XSPEC-style PHA spectral files, type

```julia
julia> fuv_grating1_phafile("ngc40", "NGC40_FUV_Grating1_IMAGE.fits","src.reg","bgd.reg", order=-2, cross_disp_width_pixel = 40, angle_xaxis_disp_deg=0.0)
-----------------------------------
target=ngc40
UVIT channel=FUV
Grating=Grating1
OBS_ID=C02_010T01_9000000
order=-2
Exposure time=1194.007 seconds
---------------------------------------
Writing source and background PHA files..
Using /soft/astrosat/responses/uvit/fuv_grating1_m2_4oct19.rmf
("ngc40_C02_010T01_9000000_FUV_Grating1_m2_spec_src.pha", "ngc40_C02_010T01_9000000_FUV_Grating1_m2_spec_bgd.pha")
```


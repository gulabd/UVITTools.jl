# UVITTools.jl
*Tools for AstroSat/UVIT grating spectroscopy and photometry.*
## Package Features
- Extract flux calbrated UVIT FUV/NUV grating spectra.
- generate source and background PHA spectral data compatible with XSPEC/Sherpa.
- Perform aperture photometry of FUV/NUV sources in any of the broadband filters.
- Saturation correction for point sources.

## UVIT Grating Spectoscopy
### Extraction of 1d count spectra
### Wavelength Calibration
### Flux calibration
### XSPEC/Sherpa-compatible PHA spectral data
### Redistribution matrix and ancillary response files

## UVIT Aperture Photometry
### Count rate, flux and magnitudes
### Saturation correction
### XSPEC/Sherpa compatible single channel spectral data and responses

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
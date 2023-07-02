"""
List of tools to analyze data from the AstroSat/UVIT payload processed with CCDLAB.

##  UVIT Photometry 
- `read_ds9reg` : 	Read DS9 region file.
- `uvit_aphot ` : Perform aper photometry on UVIT images in any of the broadband filters.
- `uvit_countrate2flux`: Convert count rate to flux density for any filter.
- `uvit_filter2pha`: Generate source and PHA spectral files compatible with XSPEC from UVIT images in any of the filters.
- `uvit_saturation_corr` : Perform saturation correction.
- `uvit_zp_uc`: Print magnitude zero point and unit conversion factor using Tandon et al. 2020.

## Grating Spectroscopy
- `xycen_from_ds9reg`: Read x,y centre for DS9 circular region file.
- `write_uvit_grating_phafile`: Write PHA spectral file for UVIT gratings.
- `fuv_grating1_count_spec`: Extract 1d count spectrum from FUV-Grating1 image file.
- `fuv_grating1_net_countrate_spec`: Derive background-corrected net count rate spectrum from FUV-Grating1 image.
- `fuv_grating1_pixel2lamA` : Convert FUV-Grating1 pixel numbers relative to zero order position into wavelength in Å.
- `fuv_grating1_wavelength_calib`: Perform wavelength calibration of FUV-Grating1 count sectrum.
- `fuv_grating1_flux_calib`: Perform flux calibration of FUV-Grating1 wavelength-calibrated count spectrum.
- `fuv_grating1_fluxed_spec`: Generate wavelength and flux calibrated spectrum from FUV-Grating1 images.
- `fuv_grating1_phafile`: Generate PHA spectral files from FUV-Grating1 images.
- `fuv_grating1_ea` : Print FUV-Grating1 effective area.
- `fuv_grating2_count_spec` : Extract 1d count spectrum from FUV-Grating2 image file.
- `fuv_grating2_net_countrate_spec`: Derive background-corrected net count rate spectrum from FUV-Grating2 image.
- `fuv_grating2_pixel2lamA`:  Convert FUV-Grating2 pixel numbers relative to zero order position into wavelength in Å.
- `fuv_grating2_wavelength_calib`: Perform wavelength calibration of FUV-Grating2 count sectrum.
- `fuv_grating2_flux_calib`: Perform flux calibration of FUV-Grating2 wavelength-calibrated count spectrum.
- `fuv_grating2_fluxed_spec`: Generate wavelength and flux calibrated spectrum from FUV-Grating2 images.
- `fuv_grating2_phafile`: Generate PHA spectral files from FUV-Grating2 images.
- `fuv_grating2_ea`: Print FUV-Grating1 effective area.
- `nuv_grating_m1_count_spec`: Extract 1d count spectrum from NUV-Grating image file.
- `nuv_grating_m1_net_countrate_spec`: Derive background-corrected net count rate spectrum from NUV-Grating image.
- `nuv_grating_m1_wavelength_calib`: Perform wavelength calibration of NUV-Grating count sectrum.
- `nuv_grating_m1_flux_calib`: Perform flux calibration of NUV-Grating wavelength-calibrated count spectrum.
- `nuv_grating_m1_fluxed_spec`: Generate wavelength and flux calibrated spectrum from NUV-Grating images.
- `nuv_grating_phafile`: Generate PHA spectral files from NUV-Grating images.

## Unit conversions
- `lambdaA2keV` : Convert wavelengh in Å to energy in keV.
- `lambdaA2ergs` : Convert λ in Å to energy in ergs.
- `flux_density_photons2cgs` : Convert photon flux density n_λ (photons/cm2/s/Å) to energy flux density f_λ (ergs/cm2/s/Å).
- `flux_density_cgs2photons` : Convert flux density f_λ (ergs/cm2/s/Å) to photon number density n_λ (photons/cm2/s/Å).


To see detailed help for any tool (e.g., `nuv_grating_phafile`), type

```julia
julia>?
help?> nuv_grating_phafile
```
"""
module UVITTools

using DelimitedFiles, Photometry
export DelimitedFiles, readdlm, writedlm
using Measurements
using DataFrames

# General purpose tools

"""
   Convert wavelength in Angstrom to Energy in keV.
"""
lambdaA2keV(lam_A::Float64) = 12.39937/lam_A

"""
   Convert wavelength in Angstrom to Energy in keV.
"""
lambdaA2keV(lam_A::Int64) = 12.39937/lam_A

"""
   Convert Energy in keV to wavelength in Angstrom.
"""
keV2lambdaA(kev::Real) = 12.39937/kev

"""
   Convert wavelength in Å to eargs.
"""
lambdaA2ergs(lam_A::Float64) = (6.6260e-27 * 3e10 / (lam_A * 1e-8))
lambdaA2ergs(lam_A::Int64) = (6.6260e-27 * 3e10 / (lam_A * 1e-8))

"""
   Convert  photon number density n_λ to flux density f_λ.
"""
flux_density_photons2cgs(f_photon_cm2_s_A::Float64, lam_A) = f_photon_cm2_s_A * lambdaA2ergs(lam_A)

"""
   Convert flux density f_λ (ergs/cm2/s/Å) to photon number density n_λ (photons/cm2/s/Å).
"""
flux_density_cgs2photons(f_ergs_cm2_s_A::Float64, lam_A) = f_ergs_cm2_s_A / lambdaA2ergs(lam_A)


#Unit conversions
export lambdaA2keV, lambdaA2ergs, flux_density_photons2cgs, flux_density_cgs2photons
# read/write useful files


# UVIT Photometry

#	include("AperturePhotometry.jl") 
#	export sum_circle, sum_circann, sum_ellipse

	include("read_ds9reg.jl")
	include("uvit_aphot.jl")
	include("write_uvit_phafile.jl")
	include("uvit_filter2pha.jl")

	export read_ds9reg, uvit_aphot, uvit_countrate2flux, uvit_filter2pha, uvit_saturation_corr, uvit_zp_uc


	# UVIT Grating Spectroscopy


	include("xycen_from_ds9reg.jl")
	include("write_uvit_grating_phafile.jl")
	include("fuv_grating1_spectroscopy.jl")
	include("fuv_grating2_spectroscopy.jl")
	include("nuv_grating_spectroscopy.jl")


	export xycen_from_ds9reg
	export write_uvit_grating_phafile
	export fuv_grating1_count_spec, fuv_grating1_net_countrate_spec, fuv_grating1_pixel2lamA, fuv_grating1_wavelength_calib, fuv_grating1_flux_calib, fuv_grating1_fluxed_spec, fuv_grating1_phafile, fuv_grating1_ea
	export fuv_grating2_count_spec, fuv_grating2_net_countrate_spec, fuv_grating2_pixel2lamA, fuv_grating2_wavelength_calib, fuv_grating2_flux_calib, fuv_grating2_fluxed_spec, fuv_grating2_phafile, fuv_grating2_ea
	export nuv_grating_m1_count_spec, nuv_grating_m1_net_countrate_spec, nuv_grating_m1_wavelength_calib, nuv_grating_m1_flux_calib, nuv_grating_m1_fluxed_spec, nuv_grating_phafile

end

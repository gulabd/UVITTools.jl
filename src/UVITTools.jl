"""
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
"""
module UVITTools

using DataFrames, FITSIO, Measurements,  DelimitedFiles, SmoothingSplines, Dates, CSV
# export DelimitedFiles, readdlm, writedlm

# using Plots
#using DataFrames

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

	# include("AperturePhotometry.jl") 
	# export sum_circle, sum_circann, sum_ellipse

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
	include("uvit_lc_from_orbitwise_images.jl")

	export xycen_from_ds9reg
	export write_uvit_grating_phafile
	export fuv_grating1_count_spec, fuv_grating1_net_countrate_spec, fuv_grating1_pixel2lamA, fuv_grating1_wavelength_calib, fuv_grating1_flux_calib, fuv_grating1_fluxed_spec, fuv_grating1_phafile, fuv_grating1_ea
	export fuv_grating2_count_spec, fuv_grating2_net_countrate_spec, fuv_grating2_pixel2lamA, fuv_grating2_wavelength_calib, fuv_grating2_flux_calib, fuv_grating2_fluxed_spec, fuv_grating2_phafile, fuv_grating2_ea
	export nuv_grating_m1_count_spec, nuv_grating_m1_net_countrate_spec, nuv_grating_m1_wavelength_calib, nuv_grating_m1_flux_calib, nuv_grating_m1_fluxed_spec, nuv_grating_phafile, uvit_lc_from_orbitwise_images

end

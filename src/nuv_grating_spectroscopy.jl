# Author: Gulab Dewangan


# For Julia 1.0 
using FITSIO, FITSIO.Libcfitsio, Dierckx, Measurements, DelimitedFiles,  DataFrames


"""
   `nuv_grating_m1_count_spec(nuv_grating_image_file, ds9regfile[, cross_disp_width_pixels = 50, rate = true, outfile="nuv_grating_m1_count_spec.dat"])` 

Extract count rate spectrum from AstroSat/UVIT NUV-Grating dispersed image generated from CCDLAB processing pipeline.

...
# Arguments
## Required parameters
- `nuv_grating_image_file::String`: Name of the NUV-Grating image file in FITS format generated using CCDLAB.
- `ds9regfile::String`: Name of the ds9 region file with center as the zero order position.
## Optional parameters
- `order::Int`: -2 (default), Grating order to be used to extract the spectrum. 
                 Allowed orders=-1 and -2.
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
- `rate::Bool`: true (default) for count rate spectrum, otherwise false for count spectrum.
- `outfile::String`: Name of ascii output file name. Default file name: "`nuv_grating_m1_count_spec.dat`".

...
"""
function nuv_grating_m1_count_spec(nuv_grating_image_file::String, ds9regfile::String; cross_disp_width_pixels::Int = 50, rate::Bool = false, outfile::String = "nuv_grating_m1_count_spec.dat")

	# Read grating image file
   	fb = FITS(nuv_grating_image_file)
   	gimg = read(fb[1])
   	exposure_time_sec = float(read_key(fb[1], "RDCDTIME")[1])
   	uvit_detector = read_key(fb[1], "DETECTOR")[1]
   	uvit_grating = read_key(fb[1], "FILTERID")[1]
   	naxis1 = read_key(fb[1], "NAXIS1")[1]

# Extract 1d spectrum

#Extract source x/y center
   	(xzero, yzero) = xycen_from_ds9reg(ds9regfile)
   	println(xzero)


# Determine the trace and cross-dispersion range for spectral extraction

# Trace is linear line : y =mx+c with m=tand(358.983) for NUV grating
#The slope and intercept are
    m = tand(358.904)
    c = yzero - m * xzero

# Need to find x-coordinate along the trace i.e., for each y.
# Number of rows i.e., number of y values is naxis2.
#xvals = range(1,step=1,length=naxis1)
    xvals = collect(1:naxis1)
    yvals = round.(Int, (m * xvals .+ c))
#println(xvals)



   	ylo = round.(Int, yvals .- cross_disp_width_pixels / 2)
	#println(xlo)
   	yhi = round.(Int, yvals .+ cross_disp_width_pixels / 2)
	#println(xhi)

	# Sum the counts from xlo to xhi for each xvalue

   	spec1d = sum.([gimg[i, ylo[i]:yhi[i]] for i in 1:naxis1])

	#display(plot(src_spec))

	#pixel_numbers = linearindices(src_spec)
  #pixel_numbers = range(1,1,length(spec1d))
    pixel_numbers = collect(1:length(spec1d))

   	pixel_num_wrt_zero_order =  pixel_numbers .- round(Int, xzero)


# Select the range of pixel numbers wrt zero order appropriate for -1 order

   	pixels_m1_order = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -546) .& (pixel_num_wrt_zero_order .< -335)]
   	spec_m1_order = spec1d[(pixel_num_wrt_zero_order .> -546) .& (pixel_num_wrt_zero_order .< -335)]
   	spec_counts_m1_order = measurement.(spec_m1_order, sqrt.(spec_m1_order))
   	spec_m1_order_counts_per_s =  spec_counts_m1_order / exposure_time_sec

   	if rate == true
      		return pixels_m1_order, spec_m1_order_counts_per_s
   	else
     		 return pixels_m1_order, spec_counts_m1_order
   	end
end


"""

   `nuv_grating_m1_net_countrate_spec(nuv_grating_image_file, ds9srcregfile, ds9bgdregfile[, cross_disp_width_pixels = 50,  outfile="nuv_grating_m1_net_countrate_spec.dat"])` 

Extract background corrected, net count rate spectrum from AstroSat/UVIT NUV-Grating dispersed image generated from CCDLAB processing pipeline.

...
# Arguments
## Required parameters
- `nuv_grating_image_file::String`: Name of the NUV-Grating image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file with source center as the zero order position.
- `ds9bgdregfile::String`: Name of the ds9 region file with  center in a source-free region of the image.
## Optional parameters
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
- `outfile::String`: Name of ascii output file name. Default file name: "`nuv_grating_m1_net_countrate_spec.dat`".
...
"""
function nuv_grating_m1_net_countrate_spec(nuv_grating_image_file::String, ds9srcregfile::String, ds9bgdregfile::String; cross_disp_width_pixels::Int = 60, rate::Bool = true,outfile::String = "nuv_grating_m1_net_count_spec.dat")
   	(pixels_m1_order, src_spec_m1_order_counts_per_s) = nuv_grating_m1_count_spec(nuv_grating_image_file, ds9srcregfile,  cross_disp_width_pixels = cross_disp_width_pixels, rate = true, outfile = "nuv_grating_m1_src_count_spec.dat")
   	(pixels_m1_order, bgd_spec_m1_order_counts_per_s) = nuv_grating_m1_count_spec(nuv_grating_image_file,  ds9bgdregfile, cross_disp_width_pixels = cross_disp_width_pixels, rate = true, outfile = "nuv_grating_m1_bgd_count_spec.dat")
   	netsrc_spec_m1_order_counts_per_s = (src_spec_m1_order_counts_per_s .- bgd_spec_m1_order_counts_per_s)
 #  	display(plot(pixels_m1_order, Measurements.value.(netsrc_spec_m1_order_counts_per_s), yerr = Measurements.uncertainty.(netsrc_spec_m1_order_counts_per_s), xlabel = "Pixel numbers wrt zero order", ylabel = "counts/s"))
   	return pixels_m1_order, netsrc_spec_m1_order_counts_per_s
end



function nuv_grating_m1_extr_count_spec(nuv_grating_image_file::String, ds9srcregfile::String, ds9bgdregfile::String, cross_disp_width_pixels::Int = 50)

  #  Use rbinned and rotated nuv grating image

   	fb = FITS(nuv_grating_image_file)
   	gimg = read(fb[1])
   	exposure_time_sec = float(read_key(fb[1], "RDCDTIME")[1])
   	uvit_detector = read_key(fb[1], "DETECTOR")[1]
   	uvit_grating = read_key(fb[1], "FILTERID")[1]
   	naxis1 = read_key(fb[1], "NAXIS1")[1]

# Extract 1d spectrum



#Extract source x/y center
   	(src_cenx, src_ceny) = xycen_from_ds9reg(ds9srcregfile)
   	println(src_cenx)
   	println(src_ceny)

# Determine the trace and cross-dispersion range for spectral extraction

# Trace is linear line : y =mx+c with m=tand(358.983) for NUV grating
#The slope and intercept are
    m = tand(358.904)
    c = src_ceny - m * src_cenx

# Need to find x-coordinate along the trace i.e., for each y.
# Number of rows i.e., number of y values is naxis2.
    xvals = range(1, 1, naxis1)
    yvals = round.(Int, (m * xvals + c))
#println(xvals)



   	ylo = round.(Int, yvals - cross_disp_width_pixels / 2)
	#println(xlo)
   	yhi = round.(Int, yvals + cross_disp_width_pixels / 2)
	#println(xhi)

	# Sum the counts from xlo to xhi for each xvalue

   	src_spec = sum.([gimg[i, ylo[i]:yhi[i]] for i in 1:naxis1])



	#display(plot(src_spec))

	#pixel_numbers = linearindices(src_spec)
    pixel_numbers = range(1, step=1, stop=length(src_spec))
   	pixel_num_wrt_zero_order =  pixel_numbers - round(Int, src_cenx)



# Select the range of pixel numbers wrt zero order appropriate for -1 order

   	pixels_m1_order = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -546) .& (pixel_num_wrt_zero_order .< -335)]
   	src_spec_m1_order = src_spec[(pixel_num_wrt_zero_order .> -546) .& (pixel_num_wrt_zero_order .< -335)]
   	src_spec_counts_m1_order = measurement.(src_spec_m1_order, sqrt.(src_spec_m1_order))



# Extract background spectrum from a source-free region
    (bgd_x, bgd_y) = xycen_from_ds9reg(ds9bgdregfile)
#    new_bgd_x = round(Int,bgd_x/2)
#    new_bgd_y = round(Int,bgd_y/2)

# Determine background range of pixels in cross-dispersion direction
   	bgd_ylo = round(Int, bgd_y - cross_disp_width_pixels / 2)
   	bgd_yhi = round(Int, bgd_y + cross_disp_width_pixels / 2)
# Extract bgd image from background extraction area with the same size as the source extraction size in x and y

    bgd_gimg = gimg[:, bgd_ylo:bgd_yhi]
    bgd_spec = sum(bgd_gimg, 2)[:,1]
    bgd_spec_m1_order = float.(bgd_spec[(pixel_num_wrt_zero_order .> -546) .& (pixel_num_wrt_zero_order .< -335)])
    bgd_spec_counts_m1_order = measurement.(bgd_spec_m1_order, sqrt.(bgd_spec_m1_order))
    netsrc_spec_m1_order_counts_per_s = (src_spec_counts_m1_order .- bgd_spec_counts_m1_order) / exposure_time_sec

# 	 display(plot(pixel_num_wrt_zero_order, src_spec, xlabel = "Pixel number w.r.t. zero order", ylabel = "Source Counts"))
#    display(plot!(pixel_num_wrt_zero_order, bgd_spec, xlabel = "Pixel number w.r.t. zero order", ylabel = "Source Counts"))
#    display(plot!(pixels_m1_order, src_spec_m1_order, yerr = sqrt.(src_spec_m1_order), ms = 1.0))

#    display(plot(pixels_m1_order, Measurements.value.(netsrc_spec_m1_order_counts_per_s), yerr = Measurements.uncertainty.(netsrc_spec_m1_order_counts_per_s), ms = 1.0))

    writedlm("ngc40_nuv_grating_1dspec_1by8pixel.dat", zip(pixels_m1_order, Measurements.value.(netsrc_spec_m1_order_counts_per_s), Measurements.uncertainty.(netsrc_spec_m1_order_counts_per_s)))
   	return pixels_m1_order, netsrc_spec_m1_order_counts_per_s
end



"""
   `nuv_grating_pixel2lamA(pixel_num_wrt_zero_order[, order = -1])`

Convert NUV-Grating (m=-1) pixel number relative to zero order to wavelength in Angstrom.

This function is used for wavelength calibration of FUV-Grating2 count spectrum.
	
## Required parameters
- `pixel_num_wrt_zero_order::Int`: Pixel numbers relative to zero order.
## Optional parameters
- `order::Int`: -1 (default), Grating order. Allowed orders=-1. Order=-2 is not yet calibrated.
...
"""
function nuv_grating_pixel2lamA(pixel_num_wrt_zero_order;order=-1)
	pixels_m1_order = pixel_num_wrt_zero_order
	if order == -1
		# (c0,c1,period,offset,ampl)=(20.806073952157906,-5.5812026743932979,9.0298864493480622,11.385668016104884,6.1811810466779304)
		# nuv_lambdaA =  c1 * pixels_m1_order + c0 + ampl * sin(2 * π * (pixels_m1_order - offset) / period)
		nuv_lambdaA =  -5.5229149729172899 * pixels_m1_order .+ 45.083829084835834
		return nuv_lambdaA
	else
		println("NUV Grating order=$order not calibrated.")
		return 0
	end
end


#=
Depricated function.

function lamA2lohi(lamA::Array{Float64,1})
    sort!(lamA)
    nbins = length(lamA)
    lamA_lo = Array{Float64}(undef,nbins)
    lamA_hi = Array{Float64}(undef,nbins)

    for i=1:nbins
        if i==1
            lamA_lo[i] = lamA[i] - (lamA[i+1] - lamA[i])/2.0
            lamA_hi[i] = lamA[i] + (lamA[i+1] - lamA[i])/2.0
        elseif i < nbins
            lamA_hi[i] = lamA[i] + (lamA[i+1] - lamA[i])/2.0
            lamA_lo[i] = lamA_hi[i-1]
        else
            lamA_lo[i] = lamA_hi[i-1]
            lamA_hi[i] = lamA[i] + (lamA[i] - lamA[i-1])/2.0
    	end
    end
    return lamA_lo, lamA_hi
end
=#


"""
   `nuv_grating_m1_wavelength_calib(pixels, netsrc_spec_counts_per_s)`

Convert pixel numbers relative to zero to wavelengths.

This function is used for wavelength calibration of NUV-Grating order=-1 count spectrum.

...
# Arguments
## Required parameters
- `pixels::Array`: An array of pixel numbers relative to zero order.
- `netsrc_spec_counts_per_s::Array`: An array of net count rates corresponding to the relative pixel numbers.
## Optional parameters
- None
...
"""
function nuv_grating_m1_wavelength_calib(pixels_m1_order, netsrc_spec_m1_order_counts_per_s)
#Wavelength calibration analysis of NGC40 data in 1/8pixel resoultion
#    nuv_lambdaA = -5.585562438209716 * pixels_m1_order .+ 18.055833282132344
#	nuv_lambdaA =  -5.5229149729172899 * pixels_m1_order .+ 45.083829084835834

	nuv_lambdaA=nuv_grating_pixel2lamA.(pixels_m1_order)
#	counts_per_s = Measurement.value.(netsrc_spec_m1_order_counts_per_s)
#	counts_per_s_err = Measurement.uncertainty.(netsrc_spec_m1_order_counts_per_s)

#	spec_df=DataFrame(nuv_lambdaA,counts_per_s,counts_per_s_err)
#	sort!(spec_df,:nuv_lambdaA)
#	lamA_sorted=spec_df[!,1]
#	counts_per_s_sorted=spec_df[!,2]
#	counts_per_s_err_sorted=spec_df[!,2]
	# (lamAlo, lamAhi) = lamA2lohi(nuv_lambdaA)
	# delta_lamA = lamAhi .- lamAlo
	# print(delta_lamA)
    netsrc_spec_m1_order_counts_per_s_A = netsrc_spec_m1_order_counts_per_s ./ 5.5229149729172899
    return nuv_lambdaA, netsrc_spec_m1_order_counts_per_s_A
end


"""
    nuv_grating_m1_ea(lamA)

Calculate NUV-Grating order=-1 effective area in cm^2 at a desired wavelength (Angstrom).

The calculation of the effective area is based on the grating calibration peformed in [Dewangan (2021)](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract). 
This function is used for flux calibration of count spectrum.

...
# Arguments
## Required
- `lam::Number`: Wavelength in Angstrom.
## Optional 
- None
...

# Example
```jldoctest
julia> nuv_grating_m1_ea(2100)
13.136720959979547
```
"""
function nuv_grating_m1_ea(lamA)
   	#ea = 900363.87080438435 - 2548.8671167487901 * l + 3.0751080479561437 * l^2  - 0.0020498923173653013 * l^3 + 8.1553032340433166e-07 * l^4 - 1.9365502778395436e-10 * l^5 + 2.5415821035881287e-14 * l^6  - 1.4222943466581621e-18 * l^7
	#ea1 = 644892.65780534293 -1834.0276469235359 * l + 2.2215485385805573 * l^2 -0.0014860733979006186*l^3  + 5.9300905470747156e-07 *l^4 -1.4117972257081421e-10 * l^5 + 1.8568903309695907e-14 * l^6 -1.0409481307291975e-18 * l^7
	(c0,c1,c2,c3,c4,c5,c6,c7)=
	(2.367139749998207,
	-0.036201927466790476,
	0.0011382352427383288,
	-5.1847986134512304e-06,
	1.0855006153897135e-08,
	-1.1914174989628051e-11,
	6.5939696397322293e-15,
	-1.4485941469673373e-18)

  # account for fixed offset parameter in polynomial fit
	l = lamA - 1900.0

	ea =  c0 + c1 * l  + c2 * l^2 + c3*l^3  + c4* l^4 + c5 * l^5 + c6 *l^6 + c7 * l^7
#
	return ea
end



"""
   `nuv_grating_m1_flux_calib(lamA, netsrc_spec_counts_per_s_A)`

Flux calibrate the wavelength-calibrated count spectrum from NUV-Grating1 order=-1.

This function calculates the effective area at each wavelength of the count spectrum, 
then uses the effective areas to convert net count rates to f_λ in CGS units.

...
# Arguments
## Required
-`lamA::Array{Float64}`: Array of wavelengths in Å.
-`netsrc_spec_counts_per_s_A::Array{Float64}`: Array of background corrected counts/s/Å corresponding to wavelength array.
## Optional
- None.
...
"""
function nuv_grating_m1_flux_calib(nuv_lambdaA, netsrc_spec_m1_order_counts_per_s_A)

  	 nuv_g_m1_ea_cm2_at_nuv_lambdaA = nuv_grating_m1_ea.(nuv_lambdaA)

# Calculate flux density
   	nuv_g_m1_n_λ = netsrc_spec_m1_order_counts_per_s_A ./ nuv_g_m1_ea_cm2_at_nuv_lambdaA
   	nuv_g_m1_f_λ = nuv_g_m1_n_λ .* lambdaA2ergs.(nuv_lambdaA)

    return nuv_lambdaA, nuv_g_m1_f_λ
end


"""
    nuv_grating_m1_fluxed_spec(target,nuv_grating_image_file, ds9srcregfile, ds9bgdregfile[,order=-1, cross_disp_width_pixels= 50, outfile="default"])

Extract flux calibrated spectrum from AstroSat/UVIT NUV-Grating order=-1 dispersed image generated from CCDLAB processing pipeline.

This is a main function that uses other functions for extraction of source and background spectra, wavelenth 
	and flux calibrations, and outputs fluxed spectrum. For details on grating orders, wavelength and flux calibrations, 
	see [Dewangan (2021)](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract).
...
# Arguments
## Required parameters
- `target::String`: Name of the target available in the observed UVIT field.
- `nuv_grating_image_file::String`: Name of the NUV-Grating image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file with source center as the zero order position.
- `ds9bgdregfile::String`: Name of the ds9 region file with  center in a source-free region of the image.
## Optional parameters
- `order::Int`: -1 (default), Grating order to be used to extract the spectrum. Allowed orders=-1. Order=-2 is not yet calibrated.
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
-`outfile::String`:Name of the output file. A default file name based on the target name/grating will be generated if outfile="default".
## Output
- `(λ, f_λ, err_f_λ)`
- Fluxed spectrum saved in an ascii file.
...
"""
function nuv_grating_m1_fluxed_spec(target::String,nuv_grating_image_file::String, ds9srcregfile::String, ds9bgdregfile::String; order::Number=-1, cross_disp_width_pixels::Int = 50, outfile::String="default")
    (pixels_m1_order, netsrc_spec_m1_order_counts_per_s) = nuv_grating_m1_net_countrate_spec(nuv_grating_image_file, ds9srcregfile, ds9bgdregfile, cross_disp_width_pixels = cross_disp_width_pixels)
   	(nuv_lambdaA, netsrc_spec_m1_order_counts_per_s_A) = nuv_grating_m1_wavelength_calib(pixels_m1_order, netsrc_spec_m1_order_counts_per_s)
   	(nuv_lambdaA, f_lambda_with_error) = nuv_grating_m1_flux_calib(nuv_lambdaA, netsrc_spec_m1_order_counts_per_s_A)
   	f_λ = Measurements.value.(f_lambda_with_error)
   	err_f_λ =  Measurements.uncertainty.(f_lambda_with_error)

	println("-----------------------------------")
	println("target=$target")
	fff = FITS(nuv_grating_image_file)
	uvit_detector = read_key(fff[1], "DETECTOR")[1]
	println("UVIT channel=$uvit_detector")
	uvit_grating = read_key(fff[1], "FILTERID")[1]
	println("Grating=$uvit_grating")
	obsid = read_key(fff[1], "OBS_ID")[1]
	println("OBS_ID=$obsid")
	println("order=$order")
	exposure_time_sec = float(read_key(fff[1], "RDCDTIME")[1])
	println("Exposure time=$exposure_time_sec seconds")
#	tstart=read_key(fff[1], "TSTART")[1]
	close(fff)
	println("---------")

	if order==-1
		gorder="m1"
	elseif order==-2
		gorder="m2"
	else
		println("Grating order not calibrated")
	end
	if outfile == "default"
		outfilename = target * "_" * obsid * "_" * uvit_detector * "_" * uvit_grating * gorder *  "_cross_disp_" * string(cross_disp_width_pixels) * "pixels_spec.dat"
	else
		outfilename=outfile
	end

	reverse!(nuv_lambdaA)
	reverse!(f_λ)
	reverse(err_f_λ)


   # display(plot(nuv_lambdaA, f_λ, yerr = err_f_λ, xlabel = L"\rm{Wavelength (\AA)}", ylabel = L"\rm{f_{\lambda} (ergs~cm^{-2}~s^{-1}~\AA^{-1})}", ms = 1, label = "NUV Grating"))
   	writedlm(outfilename, zip(nuv_lambdaA, f_λ, err_f_λ))
   	return nuv_lambdaA, f_λ, err_f_λ
end




"""
    nuv_grating_phafile(target,nuv_grating1_image_file, ds9srcregfile, ds9bgdregfile[, cross_disp_width_pixels= 50])

Extract XSPEC/Sherpa compatible source and background PHA spectral files from AstroSat/UVIT NUV-Grating  dispersed image generated from CCDLAB processing pipeline.

This function extracts source and background count spectra using the zero order positions provided 
in the DS9 region files and converts them into PHA spectral files. For details on grating orders 
and spectral responses, see [Dewangan (2021)](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract).

...
# Arguments
## Required parameters
- `target::String`: Name of the target available in the observed UVIT field.
- `nuv_grating_image_file::String`: Name of the NUV-Grating image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file with source center as the zero order position.
- `ds9bgdregfile::String`: Name of the ds9 region file with  center in a source-free region of the image.
## Optional parameters
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
- `respdir::String`: Name of the directory containing the response matrices in the local machine. The response files can be downloaded from the GitHub page.
#read necessary keywords`
## Output
- Source and background PHA files.
- Some relevant information are also printed on the screen.
...
"""
function nuv_grating_phafile(target::String,nuv_grating_image_file::String, ds9srcregfile::String, ds9bgdregfile::String; cross_disp_width_pixels = 50, respdir::String = "/soft/astrosat/responses/uvit/")
	#read necessary keywords
   	fb = FITS(nuv_grating_image_file)
   	exposure_time_sec = float(read_key(fb[1], "RDCDTIME")[1])
   	uvit_detector = read_key(fb[1], "DETECTOR")[1]
   	uvit_grating = read_key(fb[1], "FILTERID")[1]
	obsid=read_key(fb[1], "OBS_ID")[1]
	close(fb)


	# Define filenames

		srcphafile = target  * "_" * uvit_detector * "_" * uvit_grating * "m1_cross_disp_" * string(cross_disp_width_pixels) * "pixels_src"* ".pha"
		bgdphafile = target *  "_" * uvit_detector * "_" * uvit_grating * "m1_cross_disp_" * string(cross_disp_width_pixels) * "pixels_bgd" * ".pha"
	
	#Extract 1d count spectrum
   	(pixels_m1_order, src_count_spec) = nuv_grating_m1_count_spec(nuv_grating_image_file, ds9srcregfile,  cross_disp_width_pixels = cross_disp_width_pixels, rate = false, outfile = "nuv_grating_src_count_spec.dat")
   	(pixels_m1_order, bgd_count_spec) = nuv_grating_m1_count_spec(nuv_grating_image_file,  ds9bgdregfile, cross_disp_width_pixels = cross_disp_width_pixels, rate = false, outfile = "nuv_grating_bgd_count_spec.dat")

   	src_count_spec_vals = Measurements.value.(src_count_spec)
   	bgd_count_spec_vals = Measurements.value.(bgd_count_spec)
  # Define spectral channels
  #	pha_channels = convert.(Int64, linearindices(src_count_spec))
    pha_channels = collect(1:1:length(src_count_spec))
	# Write phafiles


   	src_pha_file_written = write_uvit_grating_phafile(uvit_detector, uvit_grating, pha_channels, reverse(src_count_spec_vals), exposure_time_sec, phafile = srcphafile)
	  println("Source PHA file written: $src_pha_file_written")
	  bgd_pha_file_written = write_uvit_grating_phafile(uvit_detector, uvit_grating, pha_channels, reverse(bgd_count_spec_vals), exposure_time_sec, phafile = bgdphafile)
	  println("Background PHA file written: $bgd_pha_file_written")

	# Update Respfile in source spectrum
	# rmffile=respdir * "nuv_grating_m1_9oct19.rmf"
  
	rmffile= respdir * "nuv_grating_m1_9oct19.rmf"
	println("Using respfile $rmffile")

	f=fits_open_file(srcphafile, +1)
	fits_movabs_hdu(f,2)
	fits_update_key(f,"BACKFILE",bgdphafile,"Background pha file")
	println("Background phafile updated in the header of srcpha file.")
	fits_update_key(f,"RESPFILE",rmffile,"Response matrix with effective area")
	println("rmffile updated in the header of srcpha file.")
	fits_close_file(f)
	files=(src_pha_file_written,bgd_pha_file_written)
	return files
end





using FITSIO, FITSIO.Libcfitsio, Dierckx, Measurements, DelimitedFiles,  SmoothingSplines, Dates

"""
   `fuv_grating2_count_spec(fuv_grating2_image_file, ds9regfile[, order = -2, disp_aligned_to_xaxis = false, angle_xaxis_disp_deg = 267.479, cross_disp_width_pixels = 40, rate = true])` 

Extract count rate spectrum from AstroSat/UVIT FUV-Grating2 dispersed image generated from CCDLAB processing pipeline.

...
# Arguments
## Required parameters
- `fuv_grating2_image_file::String`: Name of the FUV-Grating2 image file in FITS format generated using CCDLAB.
- `ds9regfile::String`: Name of the ds9 region file with center as the zero order position.
## Optional parameters
- `order::Int`: -2 (default), Grating order to be used to extract the spectrum. 
                 Allowed orders=-1 and -2.
- `disp_aligned_to_xaxis`: false (default), true if the dispersion axis is rotated to align to x-axis in the detector coordinates, false if dispersion axis left unrotated.
- `angle_xaxis_disp_deg`: 267.479 (default), the angle in degrees between dispersion direction (from zero to minus orders) and the x-axis. 
		If the dispersion axis is left unrotated in the detector coordiantes, the default value should be appropriate. 
		If the dispersion axis is rotated to align to x-axis, the angle should be set to 180 or zero degrees.
- `cross_disp_width_pixels::Int`: 40 (default), width in pixels in the cross-dispersion direction.
- `rate::Bool`: true (default) for count rate spectrum, otherwise false for count spectrum.

## Output
- count spectrum file
- region file compatible with ds9 that can be used to verify the extraction region.
- Count spectrum as (pixel numbers relative to zero order, counts/s or counts, errors)
...
"""
function fuv_grating2_count_spec(fuv_grating2_image_file::String,ds9regfile::String; order=-2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=50, rate::Bool=true)
	fb = FITS(fuv_grating2_image_file)
	gimg = read(fb[1])
	exposure_time_sec = float(read_key(fb[1],"RDCDTIME")[1])
	uvit_detector = read_key(fb[1],"DETECTOR")[1]
	uvit_grating = read_key(fb[1],"FILTERID")[1]
	naxis2=read_key(fb[1],"NAXIS2")[1]
	naxis1=read_key(fb[1],"NAXIS1")[1]
	close(fb)

# Extract 1d spectrum

#Extract source x/y center
	(cenx,ceny) = xycen_from_ds9reg(ds9regfile)
	

#=
 Determine the trace and cross-dispersion range for spectral extraction
Angle of dispersion direction i.e. zero order to grating order being consisdered to x axis = 267.479 degrees (can be slightly different due to registration)
Trace is linear line : y =mx+c with m=tand(267.479) for FUV grating2
 The slope and intercept are
=#
	m=tand(angle_xaxis_disp_deg)
	c = ceny - m*cenx

if disp_aligned_to_xaxis == false

# Need to find x-coordinate along the trace i.e., for each y.
# Number of rows i.e., number of y values is naxis2.

	yvals = range(1,step=1,length=naxis2)
	xvals = round.(Int, (yvals .- c)/m)

	xlo=round.(Int, xvals .- cross_disp_width_pixels/2)
	xhi=round.(Int, xvals .+ cross_disp_width_pixels/2)

	# Sum the counts from xlo to xhi for each xvalue
	spec =sum.([gimg[xlo[i]:xhi[i], i] for i in 1:naxis2])
	pixel_numbers = range(1,step=1,length=length(spec))
	pixel_num_wrt_zero_order =  pixel_numbers .- round(Int,ceny)

elseif disp_aligned_to_xaxis == true
#	println("dispersion axis aligned to x axis.")
	xvals = collect(1:naxis1)
    yvals = round.(Int, (m * xvals .+ c))
	# Define y-range along spatial direction
   	ylo = round.(Int, yvals .- cross_disp_width_pixels / 2)
   	yhi = round.(Int, yvals .+ cross_disp_width_pixels / 2)
	# Sum the counts along spatial direction
   	spec = sum.([gimg[i, ylo[i]:yhi[i]] for i in 1:naxis1])
	pixel_numbers = collect(1:length(spec))
   	pixel_num_wrt_zero_order =  pixel_numbers .- round(Int, cenx)
else
	println("""Dispersion axis is either aligned to x-axis (rotated) or kept unchanged in detector coordinates.""")
end

# Select the range of pixel numbers wrt zero order appropriate for -1 order
if order==-2
	pixels = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -641) .& (pixel_num_wrt_zero_order .< -442)]
	grating_spec = spec[(pixel_num_wrt_zero_order .> -641) .& (pixel_num_wrt_zero_order .< -442)]
	grating_spec_counts = measurement.(grating_spec, sqrt.(grating_spec))
	grating_spec_counts_per_s =  grating_spec_counts / exposure_time_sec

	# Write region file
	if disp_aligned_to_xaxis == true
		xcen_m2 = (-640 - 443)/2 + cenx
		ycen_m2 = m * xcen_m2 +  c
		xsize_m2 = -443 + 640 + 1
		ysize_m2 = cross_disp_width_pixels
	elseif disp_aligned_to_xaxis == false
		ycen_m2 = (-640 - 443)/2 + ceny
		xcen_m2 = (ycen_m2 -  c)/m
		xsize_m2 = -443 + 640 + 1
		ysize_m2 = cross_disp_width_pixels
	end
	println([xcen_m2, ycen_m2, xsize_m2, ysize_m2])
	reg = """# Region file format: DS9 version 4.1
	global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
	physical
	box($xcen_m2,$ycen_m2,$xsize_m2,$ysize_m2,$angle_xaxis_disp_deg)
	"""
# 	Construct and write region filename from input zero order region file
	outregfile = split(ds9regfile, ".")[1] * "_fuv_grating2_m2" * ".reg"
	open(outregfile, "w") do file
		write(file, reg)
	end
	outspecfile = split(ds9regfile, ".")[1] * "_fuv_grating2_m2_count_spec" * ".dat"
	if rate == true
		#display(plot(pixels_m2_order, Measurements.value.(grating_spec_counts_per_s), yerr=Measurements.uncertainty.(grating_spec_counts_per_s),xlabel="Pixel numbers wrt zero order", ylabel="counts/s"))
		writedlm(outspecfile,zip(pixels, Measurements.value.(grating_spec_counts_per_s), Measurements.uncertainty.(grating_spec_counts_per_s)))
		return pixels, Measurements.value.(grating_spec_counts_per_s), Measurements.uncertainty.(grating_spec_counts_per_s)
	else
		#display(plot(pixels_m2_order, Measurements.value.(grating_spec_counts), yerr=Measurements.uncertainty.(grating_spec),xlabel="Pixel numbers wrt zero order", ylabel="counts"))
		writedlm(outfile,zip(pixels, Measurements.value.(grating_spec_counts), Measurements.uncertainty.(grating_spec_counts)))
		return pixels, Measurements.value.(grating_spec_counts), Measurements.uncertainty.(grating_spec_counts)
	end


elseif order==-1
	# original code 
	pixels = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -334) .& (pixel_num_wrt_zero_order .< -237)]
	grating_spec = spec[(pixel_num_wrt_zero_order .> -334) .& (pixel_num_wrt_zero_order .< -237)]
	

#= Modified to check red leak
	pixels = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -414) .& (pixel_num_wrt_zero_order .< -227)]
	grating_spec = spec[(pixel_num_wrt_zero_order .> -414) .& (pixel_num_wrt_zero_order .< -227)]
=#

	grating_spec_counts = measurement.(grating_spec, sqrt.(grating_spec))
	grating_spec_counts_per_s =  grating_spec_counts / exposure_time_sec
	# Write region file
	if disp_aligned_to_xaxis == true
		xcen_m1 = (-333 - 238)/2 + cenx
		ycen_m1 = m * xcen_m1 +  c
		xsize_m1 = -238 + 333 + 1
		ysize_m1 = cross_disp_width_pixels
	elseif disp_aligned_to_xaxis == false
		ycen_m1 = (-333 - 238)/2 + ceny
		xcen_m1 = (ycen_m1 -  c)/m
		xsize_m1 = -238 + 333 + 1
		ysize_m1 = cross_disp_width_pixels
	end
	reg_m1 = """# Region file format: DS9 version 4.1
	global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
	physical
	box($xcen_m1,$ycen_m1,$xsize_m1,$ysize_m1,$angle_xaxis_disp_deg)
	"""
	# Construct output DS9 region file 
	outregfile = split(ds9regfile, ".")[1] * "_fuv_grating2_m1" * ".reg"
	open(outregfile, "w") do file
		write(file, reg_m1)
	end	
	
	outspecfile = split(ds9regfile, ".")[1] * "_fuv_grating2_m1_count_spec" * ".dat"
	if rate == true
		#display(plot(pixels_m2_order, Measurements.value.(grating_spec_counts_per_s), yerr=Measurements.uncertainty.(grating_spec_counts_per_s),xlabel="Pixel numbers wrt zero order", ylabel="counts/s"))
		writedlm(outspecfile,zip(pixels, Measurements.value.(grating_spec_counts_per_s), Measurements.uncertainty.(grating_spec_counts_per_s)))
		return pixels, grating_spec_counts_per_s
	else
		#display(plot(pixels_m2_order, Measurements.value.(grating_spec_counts), yerr=Measurements.uncertainty.(grating_spec),xlabel="Pixel numbers wrt zero order", ylabel="counts"))
		writedlm(outfile,zip(pixels, Measurements.value.(grating_spec_counts), Measurements.uncertainty.(grating_spec_counts)))
		return pixels, grating_spec_counts
	end
else
	println("Grating order $order not calibrated")

# Determine if rate or counts to be returned.


end
end


"""

   `fuv_grating2_net_countrate_spec(fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile[, order = -2, disp_aligned_to_xaxis=false, angle_xaxis_disp_deg=267.479, cross_disp_width_pixels = 40])` 

Extract background corrected, net count rate spectrum from AstroSat/UVIT FUV-Grating2 dispersed image generated from CCDLAB processing pipeline.

...
# Arguments
## Required parameters
- `fuv_grating2_image_file::String`: Name of the FUV-Grating2 image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file with source center as the zero order position.
- `ds9bgdregfile::String`: Name of the ds9 region file with  center in a source-free region of the image.
## Optional parameters
- `order::Int`: -2 (default), Grating order to be used to extract the spectrum. 
                 Allowed orders=-1 and -2.
- `disp_aligned_to_xaxis`: false (default), true if the dispersion axis is rotated to align to x-axis in the detector coordinates, false if dispersion axis left unrotated.
- `angle_xaxis_disp_deg`: 267.479 (default), the angle in degrees between dispersion direction (from zero to minus orders) and the x-axis. 
		If the dispersion axis is left unrotated in the detector coordiantes, the default value should be appropriate. 
		If the dispersion axis is rotated to align to x-axis, the angle should be set to 180 or zero degrees.
- `cross_disp_width_pixels::Int`: 40 (default), width in pixels in the cross-dispersion direction.
...
"""
function fuv_grating2_net_countrate_spec(fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order::Int = -2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=40)
	if order == -2 || order == -1
		(pixels, src_spec_counts_per_s) = fuv_grating2_count_spec(fuv_grating2_image_file, ds9srcregfile,  order=order, disp_aligned_to_xaxis=disp_aligned_to_xaxis, angle_xaxis_disp_deg=angle_xaxis_disp_deg, cross_disp_width_pixels=cross_disp_width_pixels,rate=true)
		(pixels, bgd_spec_counts_per_s) = fuv_grating2_count_spec(fuv_grating2_image_file, ds9bgdregfile, order=order, disp_aligned_to_xaxis=disp_aligned_to_xaxis, angle_xaxis_disp_deg=angle_xaxis_disp_deg, cross_disp_width_pixels=cross_disp_width_pixels,rate=true)
		netsrc_spec_counts_per_s = (src_spec_counts_per_s .- bgd_spec_counts_per_s)
	#	outfile="fuv_grating2_order" * string(order) * "_net_countrate_spec.dat"
		if order == -2
			outfile = "net_fuv_grating2_" * "m2" * "_countrate_spec" * "_cross_disp_with_" * string(cross_disp_width_pixels) * "pixels_" * string(angle_xaxis_disp_deg) * "deg.dat"
		else
			outfile = "net_fuv_grating2_" * "m1" * "_countrate_spec" * "_cross_disp_with_" * string(cross_disp_width_pixels) * "pixels_" * string(angle_xaxis_disp_deg) * "deg.dat"
		end
		writedlm(outfile, zip(pixels, Measurements.value.(netsrc_spec_counts_per_s), Measurements.uncertainty.(netsrc_spec_counts_per_s)))
	#display(plot(pixels_m2_order, Measurements.value.(netsrc_spec_m2_order_counts_per_s), yerr=Measurements.uncertainty.(netsrc_spec_m2_order_counts_per_s), xlabel="Pixel numbers wrt zero order", ylabel="counts/s"))
		println("Extracted spectrum (Net source) for order=$order , and written to file $outfile")
		return pixels, netsrc_spec_counts_per_s
	else
		println("Grating order $order not calibrated")
		return 0.0
 	 end
end


"""
   `fuv_grating2_net_countrate(fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile[, order = -2, cross_disp_width_pixels = 50])` 

Calculate background corrected, net count rate from AstroSat/UVIT FUV-Grating2 dispersed image generated from CCDLAB processing pipeline.
if order == -2
	outfile = "net_fuv_grating1_" * "m2" * "_countrate_spec" * "_cross_disp_with_" * string(cross_disp_width_pixels) * "pixels_" * string(angle_xaxis_disp_deg) * "deg.dat"
else
	outfile = "net_fuv_grating1_" * "m1" * "_countrate_spec" * "_cross_disp_with_" * string(cross_disp_width_pixels) * "pixels_" * string(angle_xaxis_disp_deg) * "deg.dat"
end
 writedlm(outfile, zip(pixels, Measurements.value.(netsrc_spec_counts_per_s), Measurements.uncertainty.(netsrc_spec_counts_per_s)))
else
 println("Grating order $order not calibrated")
end
...
# Arguments
## Required parameters
- `fuv_grating2_image_file::String`: Name of the FUV-Grating2 image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file with source center as the zero order position.
- `ds9bgdregfile::String`: Name of the ds9 region file with  center in a source-free region of the image.
## Optional parameters
- `order::Int`: -2 (default), Grating order to be used to extract the spectrum. 
                 Allowed orders=-1 and -2.
- `disp_aligned_to_xaxis`: false (default), true if the dispersion axis is rotated to align to x-axis in the detector coordinates, false if dispersion axis left unrotated.
- `angle_xaxis_disp_deg`: 267.479 (default), the angle in degrees between dispersion direction (from zero to minus orders) and the x-axis. 
	If the dispersion axis is left unrotated in the detector coordiantes, the default value should be appropriate. 
	If the dispersion axis is rotated to align to x-axis, the angle should be set to 180 or zero degrees.
- `cross_disp_width_pixels::String`: 40 (default), width in pixels in the cross-dispersion direction.
- `mst_or_bjd::String`: Print mission time (mst) or barycentric julain data (bjd). Default: "mst".
...
"""
function fuv_grating2_net_countrate(fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order=-2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=40,mst_or_bjd="mst")
	fb = FITS(fuv_grating2_image_file)
#	gimg = read(fb[1])
	exposure_time_sec = float(read_key(fb[1],"RDCDTIME")[1])
	uvit_detector = read_key(fb[1],"DETECTOR")[1]
	uvit_grating = read_key(fb[1],"FILTERID")[1]
	if mst_or_bjd=="mst"
		tstart = read_key(fb[1],"TSTART")[1]
		tstop = read_key(fb[1],"TSTOP")[1]
		meanmst=(tstart + tstop) / 2.0
	elseif mst_or_bjd=="bjd"
		meanbjd=read_key(fb[1],"MEANBJD")[1]
	else
		println("The value of keyword mst_or_vjd is not one of mst or bjd.")
	end

	close(fb)
	if order == -2 || order == -1
		(pixels, src_spec_counts) = fuv_grating2_count_spec(fuv_grating2_image_file, ds9srcregfile,  order=order, disp_aligned_to_xaxis=disp_aligned_to_xaxis, angle_xaxis_disp_deg=angle_xaxis_disp_deg, cross_disp_width_pixels=cross_disp_width_pixels,rate=false)
		(pixels, bgd_spec_counts) = fuv_grating2_count_spec(fuv_grating2_image_file,  ds9bgdregfile, order=order, disp_aligned_to_xaxis=disp_aligned_to_xaxis, angle_xaxis_disp_deg=angle_xaxis_disp_deg, cross_disp_width_pixels=cross_disp_width_pixels,rate=false)
		total_src_plus_bgd_counts = sum(src_spec_counts)
		total_bgd_counts=sum(bgd_spec_counts)
		netsrc_count_rate = (total_src_plus_bgd_counts - total_bgd_counts)/exposure_time_sec
	#	outfile="fuv_grating2_order" * string(order) * "_net_countrate_spec.dat"
	#	writedlm(outfile, zip(pixels, Measurements.value.(netsrc_spec_counts_per_s), Measurements.uncertainty.(netsrc_spec_counts_per_s)))
	#display(plot(pixels_m2_order, Measurements.value.(netsrc_spec_m2_order_counts_per_s), yerr=Measurements.uncertainty.(netsrc_spec_m2_order_counts_per_s), xlabel="Pixel numbers wrt zero order", ylabel="counts/s"))
	#	println("Extracted spectrum (Net source) for order=$order , and written to file $outfile")
		println("")
		println("=====Calculating net count rate=====")
		println("UVIT Channel: $uvit_detector")
		println("Grating: $uvit_grating")
		println("Grating order= $order")
		println("Net source count rate= $netsrc_count_rate")
		if mst_or_bjd=="mst"
			println("Mission time: $meanmst")
			println("===================================")
			return meanmst, netsrc_count_rate
		else mst_or_bjd=="bjd"
			println("BJD= $meanbjd")
			println("===================================")
			return meanbjd, netsrc_count_rate
		end
		
	else
		println("Grating order $order not calibrated")
		return 0.0
 	 end
end


"""
   `fuv_grating2_pixel2lamA(pixel_num_wrt_zero_order[, order = -2])`

Convert FUV-Grating2 pixel number relative to zero order to wavelength in Angstrom.

This function is used for wavelength calibration of FUV-Grating2 count spectrum.
	
## Required parameters
- `pixel_num_wrt_zero_order::Int`: Pixel numbers relative to zero order.
## Optional parameters
- `order::Int`: -2 (default), Grating order. Allowed orders=-1 and -2.
...
"""
function fuv_grating2_pixel2lamA(pixel_num_wrt_zero_order;order=-2)
	pixels = pixel_num_wrt_zero_order
	if order == -2
	# Updated calibration based on linear fit (seems sufficient)

		(c0,c1) = (52.471254140178075, -2.789654664753731)	
		fuv_lambdaA = c0 + c1 * pixels

	# Updated calibration based on observations performed in August 2022 (quadratic fit)
	#fuv_lambdaA = 33.0415509809 - 2.861676144524709 * pixel_num_wrt_zero_order - 6.6585613550555e-05 * pixel_num_wrt_zero_order^2
	###	
		return fuv_lambdaA
	elseif order == -1
		 fuv_lambdaA = -5.6254676*pixels .+ 45
		return fuv_lambdaA
	else
		println("Grating order=$order is not calibrated, wavelength calibration not available.")
		return -1
	end
end

"""
   `fuv_grating2_wavelength_calib(pixels, netsrc_spec_counts_per_s[, order = -2])`

Convert pixel numbers relative to zero to wavelengths.

This function is used for wavelength calibration of FUV-Grating2 count spectrum.

...
# Arguments
## Required parameters
- `pixels::Array`: An array of pixel numbers relative to zero order.
- `netsrc_spec_counts_per_s::Array`: An array of net count rates corresponding to the relative pixel numbers.
## Optional parameters
- `order::Int`: -2 (default), Grating order. Allowed orders=-1 and -2.
...
"""
function fuv_grating2_wavelength_calib(pixels, netsrc_spec_counts_per_s; order=-2)
#Wavelength calibration analysis of NGC40 data in 1/8pixel resoultion

 	fuv_lambdaA = fuv_grating2_pixel2lamA.(pixels,order=order)
	if order == -2
		netsrc_spec_counts_per_s_A = netsrc_spec_counts_per_s ./ 2.789654664753731
		return fuv_lambdaA, netsrc_spec_counts_per_s_A
	elseif order == -1
		netsrc_spec_counts_per_s_A = netsrc_spec_counts_per_s ./ 5.6046
		return fuv_lambdaA, netsrc_spec_counts_per_s_A
	else
		println("Grating order=$order is not calibrated, wavelength calibration not available.")
		return -1
	end
end


"""
    fuv_grating2_ea(lamA[,order = -2])

Calculate FUV-Grating2 effective area in cm^2 at a desired wavelength (Angstrom) and grating order -1 or -2.

The calculation of the effective area is based on the grating calibration peformed in [Dewangan (2021)](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract). 
This function is used for flux calibration of count spectrum.

...
# Arguments
## Required
- `lam::Number`: Wavelength in Angstrom.
## Optional 
-`order::Int`: Grating order -1 or -2.
...

# Example
```jldoctest
julia> fuv_grating2_ea(1450.4,order=-2)
4.076758827954109
```
"""
function fuv_grating2_ea(l; order=-2)
	l = convert.(Float64, l)
	if order == -2    

# Read effective area based on WD0308_crreject data

#f = FITS("/home/gulabd/work/julia_dev/UVITTools.jl/caldata/fuv_grating2m2_effarea_12nov22.fits")
f = FITS(joinpath(dirname(dirname(pathof(UVITTools))),"caldata", "fuv_grating2m2_effarea_12nov22.fits"))
ea_lamA = read(f[2], "X")
ea_cm2 = read(f[2],"MODEL")

spl_ea = fit(SmoothingSpline, ea_lamA, ea_cm2, 1.0)
ea = SmoothingSplines.predict(spl_ea, l)



#= Older effective area : Oct 2019
		(c0,c1,c2,c3,c4,c5)=(15052.403902313306,
			-49.434425415335248,
			0.064493346882384228,
			-4.1805856614204243e-05,
			1.3476394915698584e-08,
			-1.7298353285086891e-12)

		
		#	(c0,c1,c2,c3,c4,c5)=(14130.283637884739, -46.4858825372079, 0.06074028601196389, -3.9428721534637826e-05, 1.2727091791971908e-08, -1.6357902088292345e-12)
		ea =  c0 + c1 * l  + c2 * l^2 + c3*l^3  + c4* l^4 + c5 * l^5
		
		return ea
=#
	elseif order == -1
		# Coeeficients of best-fitting polynomial to the m1 effective area

		c_0=-1.2975923842910477
 		c_1=0.041839500032446297
 		c_2=-0.00036139043396815656
 		c_3=1.3619592871773961e-06
 		c_4=-2.2111705274426864e-09
 		c_5=1.2713230746820525e-12
		offset=1250.0
		l=l-offset

		ea =  c_0 + c_1 * l  + c_2 * l^2 + c_3*l^3  + c_4* l^4 + c_5 * l^5
		return ea
	else
		println("Grating order=$order not calibrated, effective area not available.")
		return  -1
	end
end


"""
   `fuv_grating2_flux_calib(lamA, netsrc_spec_counts_per_s_A[,order=-2])`

Flux calibrate the wavelength-calibrated count spectrum from FUV-Grating2.

This function calculates the effective area at each wavelength of the count spectrum, 
then uses the effective areas to convert net count rates to f_λ in CGS units.

...
# Arguments
## Required
-`lamA::Array{Float64}`: Array of wavelengths in Å.
-`netsrc_spec_counts_per_s_A::Array{Float64}`: Array of background corrected counts/s/Å corresponding to wavelength array.
## Optional
-`Order::Number`: Grating order -2 (default) or -1.
...
"""
function fuv_grating2_flux_calib(fuv_lambdaA, netsrc_spec_counts_per_s_A;order=-2)
	#= older code
  # read effective area file and evaluate to match the wavelength grid of spectrum
 	#fuv_g2_ea_lambdaA = convert.(Float64, readdlm("/soft/julia_load/UVITTools/responses/EA/fuv_grating2_smoothed_effarea.dat")[:,1][2:end])
 	# fuv_g2_ea_cm2 =  convert.(Float64, readdlm("/soft/julia_load/UVITTools/responses/EA/fuv_grating2_smoothed_effarea.dat")[:,2][2:end])
	# fuv_g2_ea_cm2_intrp = Spline1D(fuv_g2_ea_lambdaA, fuv_g2_ea_cm2)
  fuv_g2_ea_cm2_at_fuv_lambdaA = fuv_g2_ea_cm2_intrp(fuv_lambdaA)
=#
fuv_g2_ea_cm2_at_fuv_lambdaA = fuv_grating2_ea.(fuv_lambdaA, order=order)

# Calculate flux density
	fuv_g2_n_λ = netsrc_spec_counts_per_s_A ./fuv_g2_ea_cm2_at_fuv_lambdaA
	fuv_g2_f_λ = fuv_g2_n_λ .* lambdaA2ergs.(fuv_lambdaA)
  return fuv_lambdaA, fuv_g2_f_λ
end



"""
    fuv_grating2_fluxed_spec(target,fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile[,order = -2, disp_aligned_to_xaxis=false, angle_xaxis_disp_deg=267.479, cross_disp_width_pixels= 40])

Extract flux calibrated spectrum from AstroSat/UVIT FUV-Grating2 dispersed image generated from CCDLAB processing pipeline.

This is a main function that uses other functions for extraction of source and background spectra, wavelenth 
	and flux calibrations, and outputs fluxed spectrum. For details on grating orders, wavelength and flux calibrations, 
	see [Dewangan (2021)](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract).
...
# Arguments
## Required parameters
- `target::String`: Name of the target available in the observed UVIT field.
- `fuv_grating2_image_file::String`: Name of the FUV-Grating2 image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file with source center as the zero order position.
- `ds9bgdregfile::String`: Name of the ds9 region file with  center in a source-free region of the image.
## Optional parameters
- `order::Int`: -2 (default), Grating order to be used to extract the spectrum. Allowed orders=-1 and -2.
- `disp_aligned_to_xaxis`: false (default), true if the dispersion axis is rotated to align to x-axis in the detector coordinates, false if dispersion axis left unrotated.
- `angle_xaxis_disp_deg`: 267.479 (default), the angle in degrees between dispersion direction (from zero to minus orders) and the x-axis. 
		If the dispersion axis is left unrotated in the detector coordiantes, the default value should be appropriate. 
		If the dispersion axis is rotated to align to x-axis, the angle should be set to 180 or zero degrees.
- `cross_disp_width_pixels::String`: 40 (default), width in pixels in the cross-dispersion direction.
## Output
- `(λ, f_λ, err_f_λ)`
- Fluxed spectrum saved in an ascii file.
...
"""
function fuv_grating2_fluxed_spec(target::String, fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order=-2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=40)
	(pixels, netsrc_spec_counts_per_s)=fuv_grating2_net_countrate_spec(fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile, order=order, disp_aligned_to_xaxis=disp_aligned_to_xaxis, angle_xaxis_disp_deg=angle_xaxis_disp_deg, cross_disp_width_pixels=cross_disp_width_pixels)
	(fuv_lambdaA, netsrc_spec_counts_per_s_A) = fuv_grating2_wavelength_calib(pixels, netsrc_spec_counts_per_s,order=order)
	(fuv_lambdaA, f_lambda_with_error) = fuv_grating2_flux_calib(fuv_lambdaA,netsrc_spec_counts_per_s_A,order=order)
	f_λ = Measurements.value.(f_lambda_with_error)
	err_f_λ =  Measurements.uncertainty.(f_lambda_with_error)
  	#display(plot(fuv_lambdaA,f_λ,yerr=err_f_λ,xlabel="Wavelength (Å)", ylabel="f_λ (ergs cm^-2 s^-1 Å^-1)",ms=1))

	println("-----------------------------------")
	println("target=$target")
	fff = FITS(fuv_grating2_image_file)
	uvit_detector = read_key(fff[1], "DETECTOR")[1]
	println("UVIT channel=$uvit_detector")
	uvit_grating = read_key(fff[1], "FILTERID")[1]
	println("Grating=$uvit_grating")
	obsid = read_key(fff[1], "OBS_ID")[1]
	println("OBS_ID=$obsid")
	println("order=$order")
	exposure_time_sec = float(read_key(fff[1], "RDCDTIME")[1])
	println("Exposure time=$exposure_time_sec seconds")
	close(fff)
	println("---------------------------------------")
	if order==-1
		gorder="m1"
	elseif order==-2
		gorder="m2"
	else
		println("Grating order not calibrated")
	end

	#Reverse the spectrum to order in wavelength
	reverse!(fuv_lambdaA)
	reverse!(f_λ)
	reverse!(err_f_λ)
	outfile = target * "_" * obsid * "_" * uvit_detector * "_" * uvit_grating * gorder * "_crossdisp" * string(cross_disp_width_pixels) * "pix_" * "xax_disp_" * string(angle_xaxis_disp_deg) * "deg_spec.dat"
#	outfile = target * "_" * obsid * "_" * uvit_detector * "_" * uvit_grating * gorder *  "_cross_disp_" * string(cross_disp_width_pixels) * "pixels_spec.dat"

	writedlm(outfile,zip(fuv_lambdaA,f_λ, err_f_λ))
	return fuv_lambdaA,f_λ,err_f_λ
end



"""
    fuv_grating2_phafile(target,fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile[, order = -2, disp_aligned_to_xaxis=false, angle_xaxis_disp_deg=267.479, cross_disp_width_pixels= 40])

Extract XSPEC/Sherpa compatible source and background PHA spectral files from AstroSat/UVIT FUV-Grating2 dispersed image generated from CCDLAB processing pipeline.

This function extracts source and background count spectra using the zero order positions provided 
in the DS9 region files and converts them into PHA spectral files. For details on grating orders 
and spectral responses, see [Dewangan (2021)](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract).

...
# Arguments
## Required parameters
- `target::String`: Name of the target available in the observed UVIT field.
- `fuv_grating2_image_file::String`: Name of the FUV-Grating1 image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file with source center as the zero order position.
- `ds9bgdregfile::String`: Name of the ds9 region file with  center in a source-free region of the image.
## Optional parameters
- `order::Int`: -2 (default), Grating order to be used to extract the spectrum. Allowed orders=-1 and -2.
- `disp_aligned_to_xaxis`: false (default), true if the dispersion axis is rotated to align to x-axis in the detector coordinates, false if dispersion axis left unrotated.
- `angle_xaxis_disp_deg`: 267.479 (default), the angle in degrees between dispersion direction (from zero to minus orders) and the x-axis. 
		If the dispersion axis is left unrotated in the detector coordiantes, the default value should be appropriate. 
		If the dispersion axis is rotated to align to x-axis, the angle should be set to 180 or zero degrees.
- `cross_disp_width_pixels::Int`: 40 (default), width in pixels in the cross-dispersion direction.
## Output
- Source and background PHA files.
- Some relevant information are also printed on the screen.
...
"""
function fuv_grating2_phafile(target::String,fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order=-2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=40)
	#read necessary keywords
	fb = FITS(fuv_grating2_image_file)
	exposure_time_sec = float(read_key(fb[1],"RDCDTIME")[1])
	uvit_detector = read_key(fb[1],"DETECTOR")[1]
	println("Detectors is $uvit_detector")
	obsid = read_key(fb[1], "OBS_ID")[1]
	uvit_grating = read_key(fb[1],"FILTERID")[1]
	println("Grating is $uvit_grating")
	println("Grating order is $order")
	#Extract 1d count spectrum
	(pixels, src_count_spec) = fuv_grating2_count_spec(fuv_grating2_image_file, ds9srcregfile,  order=order, disp_aligned_to_xaxis=disp_aligned_to_xaxis, angle_xaxis_disp_deg=angle_xaxis_disp_deg, cross_disp_width_pixels=cross_disp_width_pixels,rate=false)
	(pixels, bgd_count_spec) = fuv_grating2_count_spec(fuv_grating2_image_file, ds9bgdregfile, order=order, disp_aligned_to_xaxis=disp_aligned_to_xaxis, angle_xaxis_disp_deg=angle_xaxis_disp_deg, cross_disp_width_pixels=cross_disp_width_pixels,rate=false)

	src_count_spec_vals = Measurements.value.(src_count_spec)
	bgd_count_spec_vals = Measurements.value.(bgd_count_spec)
  # Define spectral channels
	pha_channels = collect(1:length(src_count_spec))
	# Create filenames & Write phafiles
	if order==-2
		gorder="m2"
	elseif order==-1
		gorder="m1"
	else
		println("Grating order = $order not calibrated")
	end
	srcphafile = target * "_" * obsid * "_" * uvit_detector * "_" * uvit_grating * "_" * gorder * "_"  * "crossdisp" * string(cross_disp_width_pixels) *  "pix_" * "xax_disp_" * string(angle_xaxis_disp_deg) * "deg_src.pha"
	bgdphafile = target * "_" * obsid * "_" * uvit_detector * "_" * uvit_grating * "_" * gorder * "_" * "crossdisp" * string(cross_disp_width_pixels) *  "pix_" * "xax_disp_" * string(angle_xaxis_disp_deg) * "deg_bgd.pha"
#	srcphafile = target * "_" * obsid * "_" * uvit_detector * "_" * uvit_grating  * gorder * "_"  * string(cross_disp_width_pixels) * "_src_" * string(now())[1:10] * ".pha"
#	bgdphafile = target * "_" * obsid * "_" * uvit_detector * "_" * uvit_grating * gorder * "_" * string(cross_disp_width_pixels) * "_bgd_" * string(now())[1:10] * ".pha"

	#srcphafile=target * "fuv_grating2_" * gorder * "_crossdispwidth_" * string(cross_disp_width_pixels) * "pixels_" * string(now())[1:10] * "_src.pha"
	#bgdphafile=target * "fuv_grating2_" * gorder * "_crossdispwidth_" * string(cross_disp_width_pixels) * "pixels_" * string(now())[1:10] * "_bgd.pha"
	
	src_pha_file_written= write_uvit_grating_phafile(uvit_detector,uvit_grating,pha_channels,reverse(src_count_spec_vals),exposure_time_sec; phafile=srcphafile)
	bgd_pha_file_written=write_uvit_grating_phafile(uvit_detector,uvit_grating,pha_channels,reverse(bgd_count_spec_vals),exposure_time_sec; phafile=bgdphafile)
	
# Find correct response file
	respdir = joinpath(dirname(dirname(pathof(UVITTools))), "caldata")
	if  order==-2
        	rmffile=respdir * "fuv_grating2_m2_12nov22.rmf"
	elseif  order==-1
        	rmffile=respdir *"fuv_grating2_m1_3oct19.rmf"
	else
    	print("Detector/Grating not recognised, see http://uvit.iiap.res.in/Instrument")
    	print("rmf/arf filenames not updated in the PHA header.")
    	rmffile="NONE"
    #	arffile="NONE"
    end
	println("Using $rmffile")
	# Update necessary keywords
	f=fits_open_file(srcphafile, +1)
	fits_movabs_hdu(f,2)
	fits_update_key(f,"BACKFILE",bgdphafile,"Background pha file")
	fits_update_key(f,"RESPFILE",rmffile,"Response matrix with effective area")
	fits_close_file(f)

	return src_pha_file_written, bgd_pha_file_written
end

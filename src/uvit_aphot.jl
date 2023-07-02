# UVIT Photometry Tool

#  15 Feb 2018
# Gulab Dewangan (gulabd@iucaa.in)

# 1st Update
#  Instead of Filter name (Silica, etc.), filter number e.g., F3, F4 is now used.
#  The FILTER keyword in the image file store the filter number. Documentation added.


using FITSIO, Measurements, Photometry

"""
    uvit_aphot(imagefile,ds9srcregion, ds9bgdregion[, satu_corr=true, mst_or_bjd="mst"])

Perform aperture photometry of a source present in an UVIT FUV/NUV FITS image
 obtained from CCDLAB pipeline.

This function uses DS9 source and background region files and extracts count rates, corrects 
	for any saturation, calculates the background net count rate, and then converts to 
		flux density f_λ and AB magnitude using calibration  information available in 
		[Tandon et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..158T/abstract).

If the source region is a circle and satu_corr is true, then the source count rate is also corrected for saturation 
using the function `uvit_saturation_corr.jl` (see below). For details, 
	[see Tandon et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017AJ....154..128T/abstract).

...
# Arguments
## Required parameters
- `imagefile::String`: Name of the FUV or NUV image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file for the source.
- `ds9bgdregfile::String`: Name of the ds9 region file for background.
## Optional parameters
- `satu_corr::Bool`: true (default), Apply saturation correction only if source extraction region is circular. 
   Valid for point sources only and if the redius of the circular extraction region is 
    about 29 pixels or 12 arcsec (see Tandon et al. 2020).
- `mst_or_bjd::String`: Mission time or barycentric julian day, default: "mst".
...
"""
function uvit_aphot(imagefile::String,ds9srcregion::String, ds9bgdregion::String; satu_corr::Bool=true, mst_or_bjd::String="bjd")
	im_id = FITS(imagefile)
#	im_data = transpose(read(im_id[1]))
	im_data = read(im_id[1])
	im_header = read_header(im_id[1])
	uvit_filter = im_header["FILTERID"]
#	println(uvit_filter)
	detector = im_header["DETECTOR"]
	channel = detector
	frampers = im_header["FRAMPERS"]
	exposure_time_sec = float(read_key(im_id[1],"RDCDTIME")[1])

# read source and background region files
	srcreg = read_ds9reg(ds9srcregion)
	bgdreg = read_ds9reg(ds9bgdregion)

	if srcreg[1]=="circle"
		xcen = srcreg[2]
		ycen = srcreg[3]
		radius = srcreg[4]
		aper = CircularAperture(xcen, ycen, radius)
#		src_counts = photometry(aper, im_data)
#		err_src_counts = sqrt(src_counts)
	#		(src_counts, err_src_counts) = 	sum_circle(im_data,xcen,ycen, radius)
	# Apply saturation correction
	#		src_counts_satu_corr = uvit_saturation_corr(src_counts/exposure_time_sec, frames_per_sec=frampers) * exposure_time_sec		
	#    println(src_counts)
	    src_area = pi * radius^2
	elseif srcreg[1]=="annulus"
		xcen = srcreg[2]
		ycen = srcreg[3]
		r_in = srcreg[4]
		r_out = srcreg[5]
		aper = CircularAnnulus(xcen, ycen, r_in, r_out)
	#	(src_counts, err_src_counts) = 	sum_circann(im_data,xcen,ycen, inner_radius, outer_radius)
		src_area = pi * (outer_radius^2 - inner_radius^2)
	elseif srcreg[1]=="ellipse"
		xcen = srcreg[2]
		ycen = srcreg[3]
		major_radius = srcreg[4]
		minor_radius = srcreg[5]
		theta = srcreg[6]
		aper =  EllipticalAperture(xcen, ycen, major_radius, minor_radius, theta)
		#=  OLDDER code
		AperPhotometry requires theta measured w.r.t. to +X-axis in radians 
		and in the range of -pi/2 <= theta <= pi/2

		The angle in the ds9 ellipse region is measured from the positive axis.
		=#

		#=  Older code that utlised AperPhotometry
		if theta <= 90.0
			theta_rad = theta * π/180.0
		elseif theta > 90.0 && theta <= 180.0
			theta_rad = -(180.0-theta) * π/180.0
		elseif theta > 180.0 && theta <= 270.0
			theta_rad = (theta - 180.0) * π/180.0
		else theta > 270.0
			theta_rad = -(360.0 - theta) * π/180.0
		end
		println([xcen,ycen,major_radius,minor_radius,theta_rad])
		=#
	#	(src_counts, err_src_counts) = sum_ellipse(im_data,xcen,ycen,major_radius,minor_radius,theta_rad)
		src_area = pi * major_radius * minor_radius
	else
		println("Source region not recognized")
		return 0
	end

	# Perform Aperture Photometry
	src_counts = photometry(aper, im_data).aperture_sum
	err_src_counts = sqrt(src_counts)

	#if satu_corr == true
		# Saturation correction should be applied only for point sources, and use circular region of 12arcsec radius or about 29 pixels.
	# Apply saturation correction
		src_counts_satu_corr = uvit_saturation_corr(src_counts/exposure_time_sec, frames_per_sec=frampers) * exposure_time_sec
	#	err_src_counts_satu_corr = 	err_src_counts * (src_counts_satu_corr / src_counts)
	#else
	#	println("Saturation correction not applied")
	#end

	if bgdreg[1]=="circle"
		xcen_b = bgdreg[2]
		ycen_b = bgdreg[3]
		radius_b = bgdreg[4]
		aper_b = CircularAperture(xcen_b, ycen_b, radius_b)
	#    (bgd_counts, err_bgd_counts) = 	sum_circle(im_data,xcen_b,ycen_b, radius_b)
	    bgd_area = pi * radius_b^2
	elseif bgdreg[1]=="annulus"
		xcen_b = bgdreg[2]
		ycen_b = bgdreg[3]
		inner_radius_b = bgdreg[4]
		outer_radius_b = bgdreg[5]
		aper_b = CircularAnnulus(xcen_b, ycen_b, inner_radius_b, outer_radius_b)
	#	(bgd_counts, err_bgd_counts) = 	sum_circann(im_data,xcen_b,ycen_b, inner_radius_b,outer_radius_b)
		bgd_area = pi * (outer_radius_b^2  - inner_radius_b^2)
	elseif bgdreg[1]=="ellipse"
		xcen_b = bgdreg[2]
		ycen_b = bgdreg[3]
		major_radius_b = bgdreg[4]
		minor_radius_b = bgdreg[5]
		theta_b = bgdreg[6]
		aper_b =  EllipticalAperture(xcen_b, ycen_b, major_radius_b, minor_radius_b, theta_b)
	#ds9 ellipse region theta in degrees to Julia region angle in radians.
	#=	if theta_b <= 90.0
			theta_b_rad = theta_b * π/180.0
		elseif theta_b > 90.0 && theta_b <= 180.0
			theta_b_rad = -(180.0-theta_b) * π/180.0
		elseif theta_b > 180.0 && theta_b <= 270.0
			theta_b_rad = (theta_b - 180.0) * π/180.0
		else theta_b > 270.0
			theta_b_rad = -(360.0 - theta_b) * π/180.0
		end
	=#
	#	(bgd_counts, err_bgd_counts) = sum_ellipse(im_data,xcen_b,ycen_b,major_radius_b,minor_radius_b,theta_b_rad)
		bgd_area = pi * major_radius_b * minor_radius_b
	else
		print("Background region not recognized")
		return 0
	end

	bgd_counts = photometry(aper_b, im_data).aperture_sum
	# err_bgd_counts = sqrt(bgd_counts)

    scl_bgd_counts = bgd_counts * (src_area / bgd_area)


	



	src_count_rate = measurement(src_counts/exposure_time_sec, sqrt(src_counts)/exposure_time_sec)
	bgd_count_rate = measurement(scl_bgd_counts/exposure_time_sec, sqrt(bgd_counts) * (src_area / bgd_area) / exposure_time_sec)
	#err_bgd_count_rate = sqrt(bgd_counts) * (src_area / bgd_area) / exposure_time_sec
	net_count_rate = measurement((src_counts - scl_bgd_counts)/exposure_time_sec, sqrt((Measurements.uncertainty(src_count_rate))^2 + (Measurements.uncertainty(bgd_count_rate))^2))
	if srcreg[1]=="circle" && satu_corr == true
		src_counts_satu_corr = uvit_saturation_corr(src_counts/exposure_time_sec, frames_per_sec=frampers) * exposure_time_sec
		src_count_rate_satu_corr = measurement(src_counts_satu_corr/exposure_time_sec, sqrt(src_counts_satu_corr)/exposure_time_sec)
		net_count_rate_satu_corr = measurement((src_counts_satu_corr - scl_bgd_counts)/exposure_time_sec, sqrt((Measurements.uncertainty(src_count_rate_satu_corr))^2 + (Measurements.uncertainty(bgd_count_rate))^2))
	end

	if mst_or_bjd=="mst"
		tstart_mst = float(read_key(im_id[1],"TSTART")[1])
		tstop_mst = float(read_key(im_id[1],"TSTOP")[1])
		time_meanmst = (tstart_mst + tstop_mst)/2
	elseif mst_or_bjd=="bjd"
		time_meanbjd = float(read_key(im_id[1],"MEANBJD")[1])
	else 
		println("Timing info NOT found in the header")
	end

# Calculate source flux and magnitude

#= Now UC and ZP are calculated from the new routine uvit_zp_uc below, and 
the following line of code based on 2017 calibration are not required.

    if uvit_filter=="F148W" || uvit_filter=="CaF2-1" 
    	UC=measurement(3.09e-15, 2.9e-17)
    	ZP=measurement(18.016, 0.01)

   	 	elseif uvit_filter=="F154W" || uvit_filter=="BaF2"

    		UC=measurement(3.55e-15, 4.0e-17)
    		ZP=measurement(17.778, 0.01)

    	elseif uvit_filter=="F169M" || uvit_filter=="Sapphire"
    		UC=measurement(4.392e-15, 3.7e-17)
    		ZP=measurement(17.455, 0.01)

    	elseif uvit_filter=="F172M" || uvit_filter=="Silica"
    		UC=measurement(1.074e-14, 1.6e-16)

    		ZP=measurement(16.342,0.02)

    	elseif uvit_filter=="F148Wa" || uvit_filter=="CaF2-2" || uvit_filter=="CaF2"
    		UC=measurement(3.28e-15, 2.5e-17)

    		ZP=measurement(17.994, 0.01)


    	elseif uvit_filter=="N242W" || uvit_filter=="Silica-1" || uvit_filter=="Silica15"
    		UC=measurement(2.220e-16, 6.5e-19)

    		ZP=measurement(19.81, 0.002)

    	elseif uvit_filter=="N219M" || uvit_filter=="NUVB15"
    		UC=measurement(5.25e-15, 8.2e-17)
    		ZP=measurement(16.59, 0.02)

    	elseif uvit_filter=="N245M" || uvit_filter=="NUVB13"
    		UC=measurement(7.25e-16, 3.6e-18)

    		ZP=measurement(18.50, 0.07)


   		elseif uvit_filter=="N263M" || uvit_filter=="NUVB4"
			UC=measurement(8.44e-16, 9.6e-18)

    		ZP=measurement(18.18, 0.01)


    	elseif uvit_filter=="N279N" || uvit_filter=="NUVN2"
			UC=measurement(3.50e-15, 3.5e-17)
    		ZP=measurement(16.50, 0.01)

    	elseif uvit_filter=="N242Wa" || uvit_filter=="Silica-2"
			print("Photometric calibration not available for", uvit_filter)

			elseif uvit_filter=="Grating1" || uvit_filter=="Grating2"
				print("Photometric calibration not available for", uvit_filter)
    	else
    		print("Filter name not recognised, see http://uvit.iiap.res.in/Instrument/Filters")
    	return 0
		end
		
=#

	UC = uvit_zp_uc(detector, uvit_filter)[4]
	ZP= uvit_zp_uc(detector, uvit_filter)[5]

	if satu_corr == true  && srcreg[1] == "circle"
     	f_lambda = net_count_rate_satu_corr * UC
     	magnitude = -2.5log10(net_count_rate_satu_corr) + ZP
	else
		f_lambda = net_count_rate * UC
     	magnitude = -2.5log10(net_count_rate) + ZP
	end


	
		#Correct count rate for saturation using uvit_saturation_corr.

	println("---------count rates, flux and magnitude------------")
	println("Detector: ",detector)
	println("Filter: ", uvit_filter)
	println("Exposure time: ", round(exposure_time_sec; digits=1), " seconds")
#	println("Mean MST: ", time_meanmst)
	println("source+background rate = ", round(src_count_rate,digits=3)," counts/s")
	println("background rate = ", round(bgd_count_rate,digits=3),  " counts/s")
	
	if srcreg[1]=="circle" && satu_corr == true
		println("net source count rate = ", round(net_count_rate,digits=3), " counts/s")
		println("Saturation corrected net source count rate= ",round(net_count_rate_satu_corr,digits=3), " counts/s")
		println("Saturation corrected f_λ [$uvit_filter]= ", f_lambda," ergs/cm2/s/A")
		println("Saturation corrected  magnitude[$uvit_filter] (AB system) = ", round(magnitude,digits=3))
	else
		println("----No saturation correction-----")
		println("net source count rate = ", round(net_count_rate,digits=3), " counts/s")
	#	println("Net source count rate= ",round(net_count_rate_satu_corr,digits=3), " counts/s")
		println("Flux density, f_λ [$uvit_filter]= ", f_lambda," ergs/cm2/s/A")
		println("Magnitude[$uvit_filter] (AB system) = ", round(magnitude,digits=3))
	end
	println("------------------------------------------------------")

 if mst_or_bjd=="mst"
	println("Mean MST: ", time_meanmst)
	if srcreg[1]=="circle"
		return  time_meanmst, Measurements.value(net_count_rate_satu_corr), Measurements.uncertainty(net_count_rate_satu_corr)
	else
		return  time_meanmst, Measurements.value(net_count_rate), Measurements.uncertainty(net_count_rate)
	end
 elseif mst_or_bjd=="bjd"
	println("Mean BJD: ", time_meanbjd)
	println("------------------------------------------------------")
	if srcreg[1]=="circle" && satu_corr==true
		return  time_meanbjd, Measurements.value(net_count_rate_satu_corr), Measurements.uncertainty(net_count_rate_satu_corr)
	else
		return  time_meanbjd, Measurements.value(net_count_rate), Measurements.uncertainty(net_count_rate)
	end
 else
	println("Requested time col not available, returnin time in bjd")
	if srcreg[1]=="circle" && satu_corr == true
		return  time_meanbjd, Measurements.value(net_count_rate_satu_corr), Measurements.uncertainty(net_count_rate_satu_corr)
	else
		return  time_meanbjd, Measurements.value(net_count_rate), Measurements.uncertainty(net_count_rate)
	end
end


end

"""
    uvit_saturation_corr(count_rate[, frames_per_sec=28.717])

Correct UVIT source count rate for saturation effect.
	
This function is based on the algorithm provided in Tandon et al. (2017), 
and used by `uvit_aphot.jl` for circular extraction region only.

"""
function uvit_saturation_corr(count_rate::Float64;frames_per_sec::Float64=28.717)
  cpf5 = count_rate / frames_per_sec
  # cpf5 = 1-exp(-icpf5)
  icpf5 =  -log(1 - cpf5)
  icorr = icpf5 - cpf5
  rcorr = icorr * (0.89 - 0.30 * icorr^2)
  corr_cpf = (rcorr + cpf5 ) 
  corr_cps = corr_cpf * frames_per_sec
  return  corr_cps
end

"""
    uvit_zp_uc(channel, uvit_filter)

Obtain the magnitude zero point and unit conversation factor for FUV or NUV channel in any filter.

This function is based on the calibration information from Tandon et al. (2017) and (2020). 
	The function is used by `uvit_aphot.jl` to count rates to flux densities and magnitudes.

...
# Arguments
- `channel::String`: UVIT channel either FUV or NUV.
- `uvit_filter::String`: UVIT filter name. Different names can be used for the same filter e.g., "F1", "F148W", "CaF2-1".

...
"""
function uvit_zp_uc(channel::String, uvit_filter::String)
 if channel=="FUV"
    if uvit_filter=="F1" || uvit_filter=="F148W" || uvit_filter=="CaF2"
			zp=measurement(18.097,0.010)
			λ_mean=1481

		elseif uvit_filter== "F2" || uvit_filter=="F154W" || uvit_filter=="BaF2"
			zp=measurement(17.771,0.010)
			λ_mean=1541
		elseif uvit_filter=="F3" || uvit_filter=="F169M" || uvit_filter=="Sapphire"
      zp=measurement(17.410, 0.010)
			λ_mean=1608

		elseif uvit_filter=="F5" || uvit_filter=="F172M" || uvit_filter=="Silica"
			zp=measurement(16.274,0.020)
			λ_mean=1717
    end
  elseif channel=="NUV"
    if uvit_filter=="F1" || uvit_filter=="N242W" || uvit_filter=="Silica-1"
		zp=measurement(19.763,0.002)
			λ_mean=2418
		elseif uvit_filter=="F2" || uvit_filter=="N219M" || uvit_filter=="NUVB15"
			zp=measurement(16.654, 0.020)
			λ_mean=2196
   
		elseif uvit_filter=="F3" || uvit_filter=="N245M" || uvit_filter=="NUVB13"
			zp=measurement(18.452, 0.005)
			λ_mean=2447
    
		elseif uvit_filter=="F5" || uvit_filter=="N263M" || uvit_filter=="NUVB4"
      zp=measurement(18.146,0.010)
			λ_mean=2632
      
		elseif uvit_filter=="F6" || uvit_filter=="N279N" || uvit_filter=="NUVN2"
			zp=measurement(16.416,0.010)
			λ_mean=2792
    end
  end

	uc=10^(-0.4 *  (zp + 2.407))/λ_mean^2
	return channel, uvit_filter, λ_mean, uc, zp
end


"""
   uvit_countrate2flux(channel, filter, cps)

Convert net source count rate to flux density or magnitude 
for any of FUV/NUV broadband filters. This function is used by `uvit_aphot.jl`.

...
# Arguments
- `channel::String`: UVIT channel FUV or NUV
- `filter::String` : FUV or NUV filter e.g., "BaF2".
- `cps::Float64` : Counts per second.
...

"""
function uvit_countrate2flux(channel::String, filter::String, cps::Number)
	UC = uvit_zp_uc(channel, filter)[4]
	ZP= uvit_zp_uc(channel, filter)[5]

	f_lambda = cps * UC
	magnitude = -2.5log10(cps) + ZP
		 
	return f_lambda, magnitude
end

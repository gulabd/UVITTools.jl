# Tool to create source and background pha spectral files from UVIT broadband images

#  15 Feb 2018
# Gulab Dewangan (gulabd@iucaa.in)

# 1st Update
#  Instead of Filter name (Silica, etc.), filter number e.g., F3, F4 is now used.
#  The FILTER keyword in the image file store the filter number. Documentation added.

# ------------------------------------------------------------------

using FITSIO, FITSIO.Libcfitsio

"""
    uvit_filter2pha(target,imagefile, ds9srcregfile, ds9bgdregfile[,satu_corr=true, respdir=""/soft/astrosat/resp/uvit/"])

Extract XSPEC/Sherpa compatible source and background PHA spectral files from AstroSat/UVIT broadband image generated from CCDLAB processing pipeline.

This function extracts source and background counts using the DS9 region files and converts them into PHA spectral files. For details on grating orders 
and spectral responses, see [Dewangan (2021)](https://ui.adsabs.harvard.edu/abs/2021JApA...42...49D/abstract).

...
# Arguments
## Required parameters
- `target::String`: Name of the target available in the observed UVIT field.
- `imagefile::String`: Name of the FUV or NUV image file in FITS format generated using CCDLAB.
- `ds9srcregfile::String`: Name of the ds9 region file for source.
- `ds9bgdregfile::String`: Name of the ds9 region file for background.
## Optional parameters
- `satu_corr::Bool`: Apply saturation correction for point sources only if true (default). 
- `respdir::String`: Name of the directory containing the response matrices in the local machine. The response files can be downloaded from the GitHub page.
## Output
- Source and background PHA files.
- Some relevant information are also printed on the screen.
...
"""
function uvit_filter2pha(target::String,imagefile::String,ds9srcregion::String, ds9bgdregion::String; satu_corr::Bool=true, respdir::String="/soft/astrosat/resp/uvit/")
	im_id = FITS(imagefile)
#	im_data = transpose(read(im_id[1]))
	im_data = read(im_id[1])
	im_header = read_header(im_id[1])
	filterid = im_header["FILTERID"]
	filter = im_header["FILTER"]
	detector = im_header["DETECTOR"]
	obsid = im_header["OBS_ID"]
	frampers = im_header["FRAMPERS"]
	exposure_time_sec = float(read_key(im_id[1],"RDCDTIME")[1])

# Read source and background region files
	srcreg = read_ds9reg(ds9srcregion)
	bgdreg = read_ds9reg(ds9bgdregion)

	#= OLDER CODE that utilised AperturePhotometry.jl
	if srcreg[1]=="circle"
		xcen = srcreg[2]
		ycen = srcreg[3]
		radius = srcreg[4]
	    (src_counts, err_src_counts) = 	sum_circle(im_data,xcen,ycen, radius)
	#    println(src_counts)
	    src_area = pi * radius^2
	elseif srcreg[1]=="annulus"
		xcen = srcreg[2]
		ycen = srcreg[3]
		inner_radius = srcreg[4]
		outer_radius = srcreg[5]
		(src_counts, err_src_counts) = 	sum_circann(im_data,xcen,ycen, inner_radius, outer_radius)
		src_area = pi * (outer_radius^2 - inner_radius^2)
	elseif srcreg[1]=="ellipse"
		xcen = srcreg[2]
		ycen = srcreg[3]
		major_radius = srcreg[4]
		minor_radius = srcreg[5]
		theta = srcreg[6]
		(src_counts, err_src_counts) = sum_ellipse(im_data,xcen,ycen,major_radius,minor_radius,theta)
		src_area = pi * major_radius * minor_radius
	else
		println("Source region not recognized")
		return 0
	end

	if bgdreg[1]=="circle"
		xcen_b = bgdreg[2]
		ycen_b = bgdreg[3]
		radius_b = bgdreg[4]
	    (bgd_counts, err_bgd_counts) = 	sum_circle(im_data,xcen_b,ycen_b, radius_b)
	    bgd_area = pi * radius_b^2
	elseif bgdreg[1]=="annulus"
		xcen_b = bgdreg[2]
		ycen_b = bgdreg[3]
		inner_radius_b = bgdreg[4]
		outer_radius_b = bgdreg[5]
		(bgd_counts, err_bgd_counts) = 	sum_circann(im_data,xcen_b,ycen_b, inner_radius_b,outer_radius_b)
		bgd_area = pi * (outer_radius_b^2  - inner_radius_b^2)
	elseif bgdreg[1]=="ellipse"
		xcen_b = bgdreg[2]
		ycen_b = bgdreg[3]
		major_radius_b = bgdreg[4]
		minor_radius_b = bgdreg[5]
		theta_b = bgdreg[6]
		(bgd_counts, err_bgd_counts) = sum_ellipse(im_data,xcen_b,ycen_b,major_radius_b,minor_radius_b,theta_b)
		bgd_area = pi * major_radius_b * minor_radius_b
	else
		print("Background region not recognized")
		return 0
	end
    scl_bgd_counts = bgd_counts * (src_area / bgd_area)
	=#

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
	#	src_counts_satu_corr = uvit_saturation_corr(src_counts/exposure_time_sec, frames_per_sec=frampers) * exposure_time_sec
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
    scl_bgd_counts = bgd_counts * (src_area / bgd_area)

# Get count rates
	src_count_rate = src_counts/exposure_time_sec
	err_src_count_rate = sqrt(src_counts)/exposure_time_sec
	bgd_count_rate = scl_bgd_counts/exposure_time_sec
	err_bgd_count_rate = sqrt(bgd_counts) * (src_area / bgd_area) / exposure_time_sec
	net_count_rate = (src_counts - scl_bgd_counts)/exposure_time_sec
	err_net_count_rate = sqrt(err_src_count_rate^2 + err_bgd_count_rate^2)
	meanbjd = float(read_key(im_id[1],"MEANBJD")[1])
	println("---------count rates------------")
	println("Detector: ",detector)
	println("Filter: ", filter)
	println("Exposure time: ", round(exposure_time_sec; digits=1), " seconds")
	println("source+background rate = ", round(src_count_rate;digits=3),"+/-",round(err_src_count_rate;digits=4)," counts/s")
	println("background rate = ",round(bgd_count_rate;digits=3), "+/-",round(err_bgd_count_rate;digits=4)," counts/s")
	println("net source count rate = ",round(net_count_rate;digits=3), "+/-", round(err_net_count_rate;digits=4), " counts/s")
	if srcreg[1]=="circle" && satu_corr == true
		src_counts_satu_corr = uvit_saturation_corr(src_counts/exposure_time_sec, frames_per_sec=frampers) * exposure_time_sec
		src_count_rate_satu_corr = measurement(src_counts_satu_corr/exposure_time_sec, sqrt(src_counts_satu_corr)/exposure_time_sec)
		net_count_rate_satu_corr = measurement((src_counts_satu_corr - scl_bgd_counts)/exposure_time_sec, sqrt((Measurements.uncertainty(src_count_rate_satu_corr))^2 + (Measurements.uncertainty(bgd_count_rate))^2))
		println("Saturation corrected net source count rate= ",round(net_count_rate_satu_corr,digits=3), " counts/s")
		# Treat saturatation corrected source counts as the true source counts
		src_counts = src_counts_satu_corr
	end
	println("--------------------------------")
	println("---Writing PHA file-----")
	srcphafile = target * "_" * obsid * "_" * detector * "_" * filter *  "_spec_src.pha"
    bgdphafile = target * "_" * obsid * "_" * detector * "_" * filter * "_spec_bgd.pha"
	srcphafile_written = write_uvit_phafile(detector,filter,src_counts,exposure_time_sec;phafile=srcphafile)
	# Use scaled background counts for the same source and background extraction area
 	bgdphafile_written = write_uvit_phafile(detector,filter,scl_bgd_counts,exposure_time_sec;phafile=bgdphafile)
 # Update background phafile in the header of source phafile
 

		srcpha = fits_open_file(srcphafile,+1)
		fits_movabs_hdu(srcpha,2)
		#fits_write_key(srcpha, "BACKFILE", bgdphafile, "Associated background file")
		fits_update_key(srcpha,"BACKFILE",bgdphafile,"Background pha file")
    fits_close_file(srcpha)
	return srcphafile_written, bgdphafile_written
end

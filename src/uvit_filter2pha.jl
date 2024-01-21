# Tool to create source and background pha spectral files from UVIT broadband images

#  15 Feb 2018
# Gulab Dewangan (gulabd@iucaa.in)

# 1st Update
#  Instead of Filter name (Silica, etc.), filter number e.g., F3, F4 is now used.
#  The FILTER keyword in the image file store the filter number. Documentation added.

# ------------------------------------------------------------------

# using FITSIO, FITSIO.Libcfitsio

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
function uvit_filter2pha(target::String,imagefile::String,ds9srcregion::String, ds9bgdregion::String; satu_corr::Bool=true)
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
		src_area = pi * major_radius * minor_radius
	else
		println("Source region not recognized")
		return 0
	end

	# Perform Aperture Photometry
	src_counts = photometry(aper, im_data).aperture_sum
	err_src_counts = sqrt(src_counts)


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
	

# Indeitify the correct response file for a particular filter
 	if filter=="F1" || filter=="F148W" || filter=="CaF2-1"
		respfile="F148W_effarea_Tandon_etal2020.rsp"
	elseif filter== "F2" || filter=="F154W" || filter=="BaF2"
		respfile="F154W_effarea_Tandon_etal2020.rsp"
	elseif filter=="F3" || filter=="F169M" || filter=="Sapphire"
		respfile="F169M_effarea_Tandon_etal2020.rsp"
	elseif filter=="F5" || filter=="F172M" || filter=="Silica"
		respfile="F172M_effarea_Tandon_etal2020.rsp"
	elseif filter=="F7" || filter=="F148Wa" || filter=="CaF2-2"
		respfile="NONE"
	elseif filter=="F1" || filter=="N242W" || filter=="Silica-1"
		respfile="N242W_effarea_Tandon_etal2020.rsp"
	elseif filter=="F2" || filter=="N219M" || filter=="NUVB15"
		respfile="N219M_effarea_Tandon_etal2020.rsp"
	elseif filter=="F3" || filter=="N245M" || filter=="NUVB13"
		respfile="N245M_effarea_Tandon_etal2020.rsp"
	elseif filter=="F5" || filter=="N263M" || filter=="NUVB4"
		respfile="N263M_effarea_Tandon_etal2020.rsp"
	elseif filter=="F6" || filter=="N279N" || filter=="NUVN2"
		respfile="N279N_effarea_Tandon_etal2020.rsp"
	elseif filter=="F7" || filter=="N242Wa" || filter=="Silica-2"
		respfile="NONE"
	else 
		print("Filter name not recognised, see http://uvit.iiap.res.in/Instrument/Filters")
		print("rmf/arf filenames not updated in the PHA header.")
		respfile="NONE"
	end

respdir = joinpath(dirname(dirname(pathof(UVITTools))), "caldata")
if respfile=="NONE" 
	respfilefullpath=" "
else
	respfilefullpath = joinpath(respdir, respfile)
end

println("--------------------------------")
println("---Writing PHA file-----")
srcphafile = target * "_" * obsid * "_" * detector * "_" * filter *  "_spec_src.pha"
bgdphafile = target * "_" * obsid * "_" * detector * "_" * filter * "_spec_bgd.pha"
srcphafile_written = write_uvit_phafile(detector,filter,src_counts,exposure_time_sec;phafile=srcphafile)
# Use scaled background counts for the same source and background extraction area
bgdphafile_written = write_uvit_phafile(detector,filter,scl_bgd_counts,exposure_time_sec;phafile=bgdphafile)

# Update background phafile and response file in the header of source phafile
srcpha = fits_open_file(srcphafile,+1)
fits_movabs_hdu(srcpha,2)
fits_update_key(srcpha,"BACKFILE",bgdphafile,"Background pha file")
cp(respfilefullpath, joinpath(pwd(),  respfile), force=true)
fits_update_key(srcpha, "RESPFILE", respfile, "Response file, includes effective area" )
fits_close_file(srcpha)
return srcphafile_written, bgdphafile_written
end

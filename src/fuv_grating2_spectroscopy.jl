using FITSIO, FITSIO.Libcfitsio, Dierckx, Measurements, DelimitedFiles,  SmoothingSplines, Dates

"""
    `fuv_grating2_count_spec(fuv_grating2_image_file, ds9regfile[, order=-2, disp_aligned_to_xaxis=false, angle_xaxis_disp_deg=267.479, cross_disp_width_pixels=40, rate=true])`

Extract count rate spectrum from AstroSat/UVIT FUV-Grating2 dispersed image generated from CCDLAB processing pipeline.

...
# Arguments
## Required parameters
- `fuv_grating2_image_file::String`: Name of the FUV-Grating2 image file in FITS format generated using CCDLAB.
- `ds9regfile::String`: Name of the ds9 region file with center as the zero order position.
## Optional parameters
- `order::Int`: -2 (default), Grating order to be used to extract the spectrum. 
                 Allowed orders=-1 and -2.
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
- `rate::Bool`: true (default) for count rate spectrum, otherwise false for count spectrum.
- `outfile::String`: Name of ascii output file name. Default file name: "`fuv_grating2_count_spec.dat`".

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
	xvals = collect(1:naxis1)
	yvals = round.(Int, (m * xvals .+ c))

	ylo=round.(Int, yvals .- cross_disp_width_pixels/2)
	yhi=round.(Int, yvals .+ cross_disp_width_pixels/2)

	# Sum the counts from xlo to xhi for each xvalue
	spec =sum.([gimg[ylo[i]:yhi[i], i] for i in 1:naxis1])
	pixel_numbers = range(1,step=1,length=length(spec))
	pixel_num_wrt_zero_order =  pixel_numbers .- round(Int,cenx)
else
	println("""Dispersion axis is either aligned to x-axis (rotated) or kept unchanged in detector coordinates.""")
end

  	

# Select the range of pixel numbers wrt zero order appropriate for -1 order
if order==-2
	pixels = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -625) .& (pixel_num_wrt_zero_order .< -425)]
	grating_spec = spec[(pixel_num_wrt_zero_order .> -625) .& (pixel_num_wrt_zero_order .< -425)]
	grating_spec_counts = measurement.(grating_spec, sqrt.(grating_spec))
	grating_spec_counts_per_s =  grating_spec_counts / exposure_time_sec
	xcen_m2 = (-624 - 426)/2 + cenx
	ycen_m2 = m * xcen_m2 +  c
	xsize_m2 = -426 + 624 + 1
	ysize_m2 = cross_disp_width_pixels
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
		return pixels, grating_spec_counts_per_s
	else
		#display(plot(pixels_m2_order, Measurements.value.(grating_spec_counts), yerr=Measurements.uncertainty.(grating_spec),xlabel="Pixel numbers wrt zero order", ylabel="counts"))
		writedlm(outfile,zip(pixels, Measurements.value.(grating_spec_counts), Measurements.uncertainty.(grating_spec_counts)))
		return pixels, grating_spec_counts
	end


elseif order==-1
	# original code 
	pixels = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -314) .& (pixel_num_wrt_zero_order .< -227)]
	grating_spec = spec[(pixel_num_wrt_zero_order .> -314) .& (pixel_num_wrt_zero_order .< -227)]
	

#= Modified to check red leak
	pixels = pixel_num_wrt_zero_order[(pixel_num_wrt_zero_order .> -414) .& (pixel_num_wrt_zero_order .< -227)]
	grating_spec = spec[(pixel_num_wrt_zero_order .> -414) .& (pixel_num_wrt_zero_order .< -227)]
=#

	grating_spec_counts = measurement.(grating_spec, sqrt.(grating_spec))
	grating_spec_counts_per_s =  grating_spec_counts / exposure_time_sec
	xcen_m1 = (-313 - 228)/2 + cenx
	ycen_m1 = m * xcen_m1 +  c
	xsize_m1 = -228 + 313 + 1
	ysize_m1 = cross_disp_width_pixels
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
   `fuv_grating2_net_countrate_spec(fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile[, order = -2, cross_disp_width_pixels = 50,  outfile="fuv_grating2_net_countrate_spec.dat"])` 

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
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
- `outfile::String`: Name of ascii output file name. Default file name: "`fuv_grating2_net_countrate_spec.dat`".
...
"""
function fuv_grating2_net_countrate_spec(fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order::Int = -2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=50)
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
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
- `mst_or_bjd::String`: Print mission time (mst) or barycentric julain data (bjd). Default: "mst".
...
"""
function fuv_grating2_net_countrate(fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order=-2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=60,mst_or_bjd="mst")
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
"""
function fuv_grating2_ea(l; order=-2)
	l = convert.(Float64, l)
	if order == -2

	#= Updated effective area in Feb 2022
		ea_lamA = [1200.0, 1200.5, 1201.0, 1201.5, 1202.0, 1202.5, 1203.0, 1203.5, 1204.0, 1204.5, 1205.0, 1205.5, 1206.0, 1206.5, 1207.0, 1207.5, 1208.0, 1208.5, 1209.0, 1209.5, 1210.0, 1210.5, 1211.0, 1211.5, 1212.0, 1212.5, 1213.0, 1213.5, 1214.0, 1214.5, 1215.0, 1215.5, 1216.0, 1216.5, 1217.0, 1217.5, 1218.0, 1218.5, 1219.0, 1219.5, 1220.0, 1220.5, 1221.0, 1221.5, 1222.0, 1222.5, 1223.0, 1223.5, 1224.0, 1224.5, 1225.0, 1225.5, 1226.0, 1226.5, 1227.0, 1227.5, 1228.0, 1228.5, 1229.0, 1229.5, 1230.0, 1230.5, 1231.0, 1231.5, 1232.0, 1232.5, 1233.0, 1233.5, 1234.0, 1234.5, 1235.0, 1235.5, 1236.0, 1236.5, 1237.0, 1237.5, 1238.0, 1238.5, 1239.0, 1239.5, 1240.0, 1240.5, 1241.0, 1241.5, 1242.0, 1242.5, 1243.0, 1243.5, 1244.0, 1244.5, 1245.0, 1245.5, 1246.0, 1246.5, 1247.0, 1247.5, 1248.0, 1248.5, 1249.0, 1249.5, 1250.0, 1250.5, 1251.0, 1251.5, 1252.0, 1252.5, 1253.0, 1253.5, 1254.0, 1254.5, 1255.0, 1255.5, 1256.0, 1256.5, 1257.0, 1257.5, 1258.0, 1258.5, 1259.0, 1259.5, 1260.0, 1260.5, 1261.0, 1261.5, 1262.0, 1262.5, 1263.0, 1263.5, 1264.0, 1264.5, 1265.0, 1265.5, 1266.0, 1266.5, 1267.0, 1267.5, 1268.0, 1268.5, 1269.0, 1269.5, 1270.0, 1270.5, 1271.0, 1271.5, 1272.0, 1272.5, 1273.0, 1273.5, 1274.0, 1274.5, 1275.0, 1275.5, 1276.0, 1276.5, 1277.0, 1277.5, 1278.0, 1278.5, 1279.0, 1279.5, 1280.0, 1280.5, 1281.0, 1281.5, 1282.0, 1282.5, 1283.0, 1283.5, 1284.0, 1284.5, 1285.0, 1285.5, 1286.0, 1286.5, 1287.0, 1287.5, 1288.0, 1288.5, 1289.0, 1289.5, 1290.0, 1290.5, 1291.0, 1291.5, 1292.0, 1292.5, 1293.0, 1293.5, 1294.0, 1294.5, 1295.0, 1295.5, 1296.0, 1296.5, 1297.0, 1297.5, 1298.0, 1298.5, 1299.0, 1299.5, 1300.0, 1300.5, 1301.0, 1301.5, 1302.0, 1302.5, 1303.0, 1303.5, 1304.0, 1304.5, 1305.0, 1305.5, 1306.0, 1306.5, 1307.0, 1307.5, 1308.0, 1308.5, 1309.0, 1309.5, 1310.0, 1310.5, 1311.0, 1311.5, 1312.0, 1312.5, 1313.0, 1313.5, 1314.0, 1314.5, 1315.0, 1315.5, 1316.0, 1316.5, 1317.0, 1317.5, 1318.0, 1318.5, 1319.0, 1319.5, 1320.0, 1320.5, 1321.0, 1321.5, 1322.0, 1322.5, 1323.0, 1323.5, 1324.0, 1324.5, 1325.0, 1325.5, 1326.0, 1326.5, 1327.0, 1327.5, 1328.0, 1328.5, 1329.0, 1329.5, 1330.0, 1330.5, 1331.0, 1331.5, 1332.0, 1332.5, 1333.0, 1333.5, 1334.0, 1334.5, 1335.0, 1335.5, 1336.0, 1336.5, 1337.0, 1337.5, 1338.0, 1338.5, 1339.0, 1339.5, 1340.0, 1340.5, 1341.0, 1341.5, 1342.0, 1342.5, 1343.0, 1343.5, 1344.0, 1344.5, 1345.0, 1345.5, 1346.0, 1346.5, 1347.0, 1347.5, 1348.0, 1348.5, 1349.0, 1349.5, 1350.0, 1350.5, 1351.0, 1351.5, 1352.0, 1352.5, 1353.0, 1353.5, 1354.0, 1354.5, 1355.0, 1355.5, 1356.0, 1356.5, 1357.0, 1357.5, 1358.0, 1358.5, 1359.0, 1359.5, 1360.0, 1360.5, 1361.0, 1361.5, 1362.0, 1362.5, 1363.0, 1363.5, 1364.0, 1364.5, 1365.0, 1365.5, 1366.0, 1366.5, 1367.0, 1367.5, 1368.0, 1368.5, 1369.0, 1369.5, 1370.0, 1370.5, 1371.0, 1371.5, 1372.0, 1372.5, 1373.0, 1373.5, 1374.0, 1374.5, 1375.0, 1375.5, 1376.0, 1376.5, 1377.0, 1377.5, 1378.0, 1378.5, 1379.0, 1379.5, 1380.0, 1380.5, 1381.0, 1381.5, 1382.0, 1382.5, 1383.0, 1383.5, 1384.0, 1384.5, 1385.0, 1385.5, 1386.0, 1386.5, 1387.0, 1387.5, 1388.0, 1388.5, 1389.0, 1389.5, 1390.0, 1390.5, 1391.0, 1391.5, 1392.0, 1392.5, 1393.0, 1393.5, 1394.0, 1394.5, 1395.0, 1395.5, 1396.0, 1396.5, 1397.0, 1397.5, 1398.0, 1398.5, 1399.0, 1399.5, 1400.0, 1400.5, 1401.0, 1401.5, 1402.0, 1402.5, 1403.0, 1403.5, 1404.0, 1404.5, 1405.0, 1405.5, 1406.0, 1406.5, 1407.0, 1407.5, 1408.0, 1408.5, 1409.0, 1409.5, 1410.0, 1410.5, 1411.0, 1411.5, 1412.0, 1412.5, 1413.0, 1413.5, 1414.0, 1414.5, 1415.0, 1415.5, 1416.0, 1416.5, 1417.0, 1417.5, 1418.0, 1418.5, 1419.0, 1419.5, 1420.0, 1420.5, 1421.0, 1421.5, 1422.0, 1422.5, 1423.0, 1423.5, 1424.0, 1424.5, 1425.0, 1425.5, 1426.0, 1426.5, 1427.0, 1427.5, 1428.0, 1428.5, 1429.0, 1429.5, 1430.0, 1430.5, 1431.0, 1431.5, 1432.0, 1432.5, 1433.0, 1433.5, 1434.0, 1434.5, 1435.0, 1435.5, 1436.0, 1436.5, 1437.0, 1437.5, 1438.0, 1438.5, 1439.0, 1439.5, 1440.0, 1440.5, 1441.0, 1441.5, 1442.0, 1442.5, 1443.0, 1443.5, 1444.0, 1444.5, 1445.0, 1445.5, 1446.0, 1446.5, 1447.0, 1447.5, 1448.0, 1448.5, 1449.0, 1449.5, 1450.0, 1450.5, 1451.0, 1451.5, 1452.0, 1452.5, 1453.0, 1453.5, 1454.0, 1454.5, 1455.0, 1455.5, 1456.0, 1456.5, 1457.0, 1457.5, 1458.0, 1458.5, 1459.0, 1459.5, 1460.0, 1460.5, 1461.0, 1461.5, 1462.0, 1462.5, 1463.0, 1463.5, 1464.0, 1464.5, 1465.0, 1465.5, 1466.0, 1466.5, 1467.0, 1467.5, 1468.0, 1468.5, 1469.0, 1469.5, 1470.0, 1470.5, 1471.0, 1471.5, 1472.0, 1472.5, 1473.0, 1473.5, 1474.0, 1474.5, 1475.0, 1475.5, 1476.0, 1476.5, 1477.0, 1477.5, 1478.0, 1478.5, 1479.0, 1479.5, 1480.0, 1480.5, 1481.0, 1481.5, 1482.0, 1482.5, 1483.0, 1483.5, 1484.0, 1484.5, 1485.0, 1485.5, 1486.0, 1486.5, 1487.0, 1487.5, 1488.0, 1488.5, 1489.0, 1489.5, 1490.0, 1490.5, 1491.0, 1491.5, 1492.0, 1492.5, 1493.0, 1493.5, 1494.0, 1494.5, 1495.0, 1495.5, 1496.0, 1496.5, 1497.0, 1497.5, 1498.0, 1498.5, 1499.0, 1499.5, 1500.0, 1500.5, 1501.0, 1501.5, 1502.0, 1502.5, 1503.0, 1503.5, 1504.0, 1504.5, 1505.0, 1505.5, 1506.0, 1506.5, 1507.0, 1507.5, 1508.0, 1508.5, 1509.0, 1509.5, 1510.0, 1510.5, 1511.0, 1511.5, 1512.0, 1512.5, 1513.0, 1513.5, 1514.0, 1514.5, 1515.0, 1515.5, 1516.0, 1516.5, 1517.0, 1517.5, 1518.0, 1518.5, 1519.0, 1519.5, 1520.0, 1520.5, 1521.0, 1521.5, 1522.0, 1522.5, 1523.0, 1523.5, 1524.0, 1524.5, 1525.0, 1525.5, 1526.0, 1526.5, 1527.0, 1527.5, 1528.0, 1528.5, 1529.0, 1529.5, 1530.0, 1530.5, 1531.0, 1531.5, 1532.0, 1532.5, 1533.0, 1533.5, 1534.0, 1534.5, 1535.0, 1535.5, 1536.0, 1536.5, 1537.0, 1537.5, 1538.0, 1538.5, 1539.0, 1539.5, 1540.0, 1540.5, 1541.0, 1541.5, 1542.0, 1542.5, 1543.0, 1543.5, 1544.0, 1544.5, 1545.0, 1545.5, 1546.0, 1546.5, 1547.0, 1547.5, 1548.0, 1548.5, 1549.0, 1549.5, 1550.0, 1550.5, 1551.0, 1551.5, 1552.0, 1552.5, 1553.0, 1553.5, 1554.0, 1554.5, 1555.0, 1555.5, 1556.0, 1556.5, 1557.0, 1557.5, 1558.0, 1558.5, 1559.0, 1559.5, 1560.0, 1560.5, 1561.0, 1561.5, 1562.0, 1562.5, 1563.0, 1563.5, 1564.0, 1564.5, 1565.0, 1565.5, 1566.0, 1566.5, 1567.0, 1567.5, 1568.0, 1568.5, 1569.0, 1569.5, 1570.0, 1570.5, 1571.0, 1571.5, 1572.0, 1572.5, 1573.0, 1573.5, 1574.0, 1574.5, 1575.0, 1575.5, 1576.0, 1576.5, 1577.0, 1577.5, 1578.0, 1578.5, 1579.0, 1579.5, 1580.0, 1580.5, 1581.0, 1581.5, 1582.0, 1582.5, 1583.0, 1583.5, 1584.0, 1584.5, 1585.0, 1585.5, 1586.0, 1586.5, 1587.0, 1587.5, 1588.0, 1588.5, 1589.0, 1589.5, 1590.0, 1590.5, 1591.0, 1591.5, 1592.0, 1592.5, 1593.0, 1593.5, 1594.0, 1594.5, 1595.0, 1595.5, 1596.0, 1596.5, 1597.0, 1597.5, 1598.0, 1598.5, 1599.0, 1599.5, 1600.0, 1600.5, 1601.0, 1601.5, 1602.0, 1602.5, 1603.0, 1603.5, 1604.0, 1604.5, 1605.0, 1605.5, 1606.0, 1606.5, 1607.0, 1607.5, 1608.0, 1608.5, 1609.0, 1609.5, 1610.0, 1610.5, 1611.0, 1611.5, 1612.0, 1612.5, 1613.0, 1613.5, 1614.0, 1614.5, 1615.0, 1615.5, 1616.0, 1616.5, 1617.0, 1617.5, 1618.0, 1618.5, 1619.0, 1619.5, 1620.0, 1620.5, 1621.0, 1621.5, 1622.0, 1622.5, 1623.0, 1623.5, 1624.0, 1624.5, 1625.0, 1625.5, 1626.0, 1626.5, 1627.0, 1627.5, 1628.0, 1628.5, 1629.0, 1629.5, 1630.0, 1630.5, 1631.0, 1631.5, 1632.0, 1632.5, 1633.0, 1633.5, 1634.0, 1634.5, 1635.0, 1635.5, 1636.0, 1636.5, 1637.0, 1637.5, 1638.0, 1638.5, 1639.0, 1639.5, 1640.0, 1640.5, 1641.0, 1641.5, 1642.0, 1642.5, 1643.0, 1643.5, 1644.0, 1644.5, 1645.0, 1645.5, 1646.0, 1646.5, 1647.0, 1647.5, 1648.0, 1648.5, 1649.0, 1649.5, 1650.0, 1650.5, 1651.0, 1651.5, 1652.0, 1652.5, 1653.0, 1653.5, 1654.0, 1654.5, 1655.0, 1655.5, 1656.0, 1656.5, 1657.0, 1657.5, 1658.0, 1658.5, 1659.0, 1659.5, 1660.0, 1660.5, 1661.0, 1661.5, 1662.0, 1662.5, 1663.0, 1663.5, 1664.0, 1664.5, 1665.0, 1665.5, 1666.0, 1666.5, 1667.0, 1667.5, 1668.0, 1668.5, 1669.0, 1669.5, 1670.0, 1670.5, 1671.0, 1671.5, 1672.0, 1672.5, 1673.0, 1673.5, 1674.0, 1674.5, 1675.0, 1675.5, 1676.0, 1676.5, 1677.0, 1677.5, 1678.0, 1678.5, 1679.0, 1679.5, 1680.0, 1680.5, 1681.0, 1681.5, 1682.0, 1682.5, 1683.0, 1683.5, 1684.0, 1684.5, 1685.0, 1685.5, 1686.0, 1686.5, 1687.0, 1687.5, 1688.0, 1688.5, 1689.0, 1689.5, 1690.0, 1690.5, 1691.0, 1691.5, 1692.0, 1692.5, 1693.0, 1693.5, 1694.0, 1694.5, 1695.0, 1695.5, 1696.0, 1696.5, 1697.0, 1697.5, 1698.0, 1698.5, 1699.0, 1699.5, 1700.0, 1700.5, 1701.0, 1701.5, 1702.0, 1702.5, 1703.0, 1703.5, 1704.0, 1704.5, 1705.0, 1705.5, 1706.0, 1706.5, 1707.0, 1707.5, 1708.0, 1708.5, 1709.0, 1709.5, 1710.0, 1710.5, 1711.0, 1711.5, 1712.0, 1712.5, 1713.0, 1713.5, 1714.0, 1714.5, 1715.0, 1715.5, 1716.0, 1716.5, 1717.0, 1717.5, 1718.0, 1718.5, 1719.0, 1719.5, 1720.0, 1720.5, 1721.0, 1721.5, 1722.0, 1722.5, 1723.0, 1723.5, 1724.0, 1724.5, 1725.0, 1725.5, 1726.0, 1726.5, 1727.0, 1727.5, 1728.0, 1728.5, 1729.0, 1729.5, 1730.0, 1730.5, 1731.0, 1731.5, 1732.0, 1732.5, 1733.0, 1733.5, 1734.0, 1734.5, 1735.0, 1735.5, 1736.0, 1736.5, 1737.0, 1737.5, 1738.0, 1738.5, 1739.0, 1739.5, 1740.0, 1740.5, 1741.0, 1741.5, 1742.0, 1742.5, 1743.0, 1743.5, 1744.0, 1744.5, 1745.0, 1745.5, 1746.0, 1746.5, 1747.0, 1747.5, 1748.0, 1748.5, 1749.0, 1749.5, 1750.0, 1750.5, 1751.0, 1751.5, 1752.0, 1752.5, 1753.0, 1753.5, 1754.0, 1754.5, 1755.0, 1755.5, 1756.0, 1756.5, 1757.0, 1757.5, 1758.0, 1758.5, 1759.0, 1759.5, 1760.0, 1760.5, 1761.0, 1761.5, 1762.0, 1762.5, 1763.0, 1763.5, 1764.0, 1764.5, 1765.0, 1765.5, 1766.0, 1766.5, 1767.0, 1767.5, 1768.0, 1768.5, 1769.0, 1769.5, 1770.0, 1770.5, 1771.0, 1771.5, 1772.0, 1772.5, 1773.0, 1773.5, 1774.0, 1774.5, 1775.0, 1775.5, 1776.0, 1776.5, 1777.0, 1777.5, 1778.0, 1778.5, 1779.0, 1779.5, 1780.0, 1780.5, 1781.0, 1781.5, 1782.0, 1782.5, 1783.0, 1783.5, 1784.0, 1784.5, 1785.0, 1785.5, 1786.0, 1786.5, 1787.0, 1787.5, 1788.0, 1788.5, 1789.0, 1789.5, 1790.0, 1790.5, 1791.0, 1791.5, 1792.0, 1792.5, 1793.0, 1793.5, 1794.0, 1794.5, 1795.0, 1795.5, 1796.0, 1796.5, 1797.0, 1797.5, 1798.0, 1798.5, 1799.0, 1799.5, 1800.0]
		ea_cm2 = [1.1893044191444615, 1.1672265475024037, 1.1456906600414547, 1.1246806883155915, 1.1041812064530345, 1.0841774034980876, 1.0646550482125403, 1.0456004647146493, 1.0270004996252113, 1.0088425042663836, 0.9911143017226299, 0.973804170878594, 0.9569008241787235, 0.9403933835896687, 0.9242713630641562, 0.9085246546101113, 0.8931435068739301, 0.8781185090676398, 0.8634405771398055, 0.849100940709263, 0.8350911254597723, 0.8214029446030066, 0.8080284824977085, 0.7949600877661792, 0.7821903562920174, 0.7697121254784189, 0.7575184640440705, 0.7456026585520473, 0.73395820904956, 0.7225788168199861, 0.7114583798576565, 0.7005909823788554, 0.6899708881151211, 0.6795925338677394, 0.6694505223509345, 0.6595396143233488, 0.6498547255510754, 0.6403909201877523, 0.6311434033144615, 0.6221075152545886, 0.6132787320815486, 0.6046526521106153, 0.5962249993448808, 0.5879916137764774, 0.5799484510637738, 0.5720915738116522, 0.56441715141249, 0.5569214565385747, 0.5496008586089962, 0.5424518224283917, 0.5354709055104077, 0.5286547531648547, 0.5220000962026426, 0.5155037491712711, 0.5091626054756534, 0.5029736376660914, 0.4969338912203836, 0.49104048597517197, 0.485290609284705, 0.47968151914230595, 0.4742105370156538, 0.468875050044639, 0.46367250486326905, 0.45860040946032643, 0.45365632833008723, 0.44883819874303144, 0.4441476071618776, 0.4395882389920643, 0.4351634813775167, 0.4308764122886827, 0.4267298341325281, 0.42272673289862206, 0.41887137075893016, 0.4151680354324712, 0.4116207338063314, 0.4082332342502963, 0.40500909977413024, 0.4019515926236682, 0.3990636851802471, 0.39634827493703995, 0.3938082368271634, 0.39144645820008694, 0.38926569538783484, 0.3872674402349888, 0.38545273693869697, 0.38382284582607745, 0.38237926383443577, 0.38112374817356764, 0.3800574864689228, 0.3791799060364599, 0.37849048695597154, 0.3779890876007985, 0.37767594359056983, 0.37755158787993587, 0.377615321491849, 0.37786538107969, 0.37830032493427296, 0.3789190739573794, 0.3797208951574975, 0.3807051045941589, 0.3818698331300166, 0.3832130074918529, 0.3847327974538444, 0.3864275922633765, 0.3882959790206419, 0.39033649684528055, 0.3925474397633505, 0.39492717561362467, 0.39747415433344996, 0.40018688784965323, 0.4030639527726415, 0.40610428251403663, 0.4093070426815261, 0.41267139356923094, 0.4161964755614095, 0.4198813775718463, 0.4237252913607715, 0.42772793620417554, 0.4318891017121496, 0.43620849806204126, 0.4406857393638875, 0.445320319265712, 0.45011158367597265, 0.4550587025494151, 0.46016065664170513, 0.4654162136782854, 0.4708239024796998, 0.4763819102141887, 0.48208756062326674, 0.48793753917320676, 0.49392810358979744, 0.5000550450784633, 0.5063136555328883, 0.5126986503970894, 0.5192040735076999, 0.525823336138174, 0.5325491918392695, 0.5393737039865654, 0.5462882675217872, 0.5532845754442165, 0.5603545733648644, 0.567489603422169, 0.5746803576698286, 0.5819168644055117, 0.5891889576487456, 0.5964882717778165, 0.6038067706314859, 0.6111360201518182, 0.6184671938945223, 0.6257911033130239, 0.6330995032250876, 0.6403861653222764, 0.6476449872907153, 0.6548697946498083, 0.6620543576263043, 0.6691925692850237, 0.6762802215406015, 0.6833144453510055, 0.6902925744884774, 0.6972121403398778, 0.7040708790763884, 0.7108671875401665, 0.717601291762278, 0.7242740772495164, 0.7308867578072655, 0.7374408685943016, 0.7439382790666312, 0.7503817178746645, 0.7567749216333771, 0.7631220139538238, 0.7694274563509704, 0.7756960404372611, 0.781932890950087, 0.7881435416004035, 0.7943338779473493, 0.8005100839750533, 0.8066786360820872, 0.812846295322929, 0.8190198397607763, 0.8252056937887989, 0.8314104294806369, 0.8376408603508637, 0.8439040427818222, 0.8502071980025333, 0.8565565175070783, 0.8629572713661839, 0.869414862865019, 0.875934874191797, 0.8825230675471293, 0.8891848657478959, 0.8959232307755077, 0.9027403903595435, 0.9096386534989257, 0.9166204051757411, 0.9236880966655096, 0.9308428736337514, 0.9380835530631549, 0.9454086908308355, 0.952816798575565, 0.9603063432305867, 0.9678756355969097, 0.9755217256439674, 0.9832407937453714, 0.9910288920115483, 0.9988819381428343, 1.006795716433429, 1.0147659683032393, 1.0227886062902989, 1.0308594786606293, 1.038974312064992, 1.0471286986753816, 1.055318127720081, 1.063538784722249, 1.071787760440462, 1.0800621336250085, 1.0883589146592205, 1.0966750436736075, 1.1050075737355038, 1.1133548319441693, 1.1217157255351078, 1.1300891588132311, 1.138474039533887, 1.1468692715987492, 1.1552741045658521, 1.1636886115875789, 1.1721130225615102, 1.1805476059683204, 1.1889926636096273, 1.1974485328063607, 1.205915453708995, 1.2143935841392401, 1.2228831096052242, 1.2313842460368334, 1.2398972420853995, 1.248422254891361, 1.2569588819826871, 1.265506542484868, 1.2740646680368803, 1.2826326962469745, 1.2912100723484534, 1.2997962792715638, 1.3083908596201985, 1.3169933667002953, 1.3256033699205625, 1.3342204461559655, 1.3428442701579817, 1.3514755106430605, 1.36011548274511, 1.3687655647925012, 1.3774271851185376, 1.386101825086248, 1.394791485429971, 1.4034998176108664, 1.4122309131964825, 1.4209889713682002, 1.4297783039878238, 1.4386033553549749, 1.44746944184533, 1.4563829213123431, 1.4653503616140768, 1.474378494086221, 1.4834742181539111, 1.492644560564328, 1.5018964194936604, 1.5112367366335768, 1.5206726261951962, 1.5302113911714261, 1.5398605236002887, 1.5496270129133138, 1.55951636805692, 1.569534018818806, 1.579685548589727, 1.589976708872719, 1.6004133167295567, 1.6109994677637727, 1.6217377437524236, 1.6326307632163604, 1.6436812279799986, 1.6548919333788115, 1.6662652014096995, 1.6778006815702216, 1.6894971612950902, 1.701353386218803, 1.7133680576130017, 1.725539802615036, 1.737865849304283, 1.7503410684062088, 1.7629599241888465, 1.7757166664614974, 1.788605319408486, 1.8016195258733747, 1.8147510656656256, 1.8279904015768182, 1.8413276092436626, 1.8547523854291066, 1.8682540317289968, 1.8818213528237648, 1.8954424189965897, 1.9091048009971532, 1.9227956243872681, 1.9365015630015479, 1.950208880550111, 1.963904655748303, 1.9775770710082814, 1.991214041551063, 2.0048031382977274, 2.0183315867805893, 2.0317868548554996, 2.0451602102201316, 2.0584445588250024, 2.0716327012764957, 2.08471734310066, 2.0976910998807647, 2.110548468164532, 2.123288510331395, 2.1359110493922775, 2.1484160346730476, 2.1608035412319224, 2.173073913815602, 2.185230215528346, 2.197277921701692, 2.2092228175756286, 2.221070918243559, 2.232828472128663, 2.244502241670039, 2.256100629525429, 2.2676327190151278, 2.279107849516709, 2.2905356076821426, 2.3019258280365307, 2.313288614523499, 2.3246343734188626, 2.3359737493297072, 2.34731762887382, 2.358677137810625, 2.370063517746889, 2.3814869377605166, 2.392956965628996, 2.4044833691045575, 2.4160761224928966, 2.427745411446308, 2.439500638319622, 2.451348055370474, 2.4632933662749927, 2.4753423964610892, 2.487501107197141, 2.499775498764881, 2.5121685940553724, 2.5246798146709715, 2.53730836315652, 2.5500534236134804, 2.5629141576395122, 2.575888985829522, 2.5889711530784334, 2.602151515316494, 2.6154206846743335, 2.628769033156646, 2.6421866749144587, 2.655661625676304, 2.669177331306422, 2.6827162211486115, 2.6962602776566613, 2.70979103777297, 2.7232895377250648, 2.7367355069202217, 2.750107479714277, 2.763383488641766, 2.7765411035648717, 2.789557442935259, 2.8024098726882, 2.8150786289037963, 2.827544625296315, 2.8397885145885056, 2.8517907231479045, 2.8635315151282104, 2.8749946464992804, 2.8861697355035423, 2.8970469936882823, 2.907616702737505, 2.917869230255953, 2.9277957359954323, 2.9373948010267, 2.946669651035417, 2.9556238429653194, 2.964261203964079, 2.9725858156973084, 2.980604147182682, 2.9883299222043394, 2.9957785814662654, 3.002965829318342, 3.0099076095921324, 3.0166201632975596, 3.0231227503078206, 3.0294380818496522, 3.035589218885526, 3.0415993587537313, 3.0474918177463, 3.053290229140832, 3.0590197943076705, 3.0647064326762985, 3.070376143429842, 3.0760550032378897, 3.0817691663641824, 3.0875445737184584, 3.093406589361731, 3.0993805921030315, 3.1054920915835478, 3.111766749474195, 3.1182302507302637, 3.1249058231685964, 3.131814694428293, 3.1389782749149244, 3.146418254875502, 3.1541566489016897, 3.1622144497520384, 3.1706064887866305, 3.1793460260536808, 3.1884466931206727, 3.197922516462609, 3.2077878771700754, 3.2180524597302886, 3.228718075081094, 3.239786065672773, 3.251258054469895, 3.2631359395504314, 3.2754210554097547, 3.2881062189759755, 3.301179045337068, 3.314626997521302, 3.3284374206060825, 3.342597497635754, 3.357092376070579, 3.3719009182393935, 3.3870002981083993, 3.4023671401993276, 3.41797747013815, 3.433806686838294, 3.4498297637583426, 3.4660212741764442, 3.482355066248713, 3.4988042179907843, 3.5153410078613145, 3.531937588155796, 3.5485700780699148, 3.565215926398028, 3.581851940626883, 3.5984542775294224, 3.6149984480234543, 3.6314620505983046, 3.647828309313166, 3.6640808956491515, 3.680203107386708, 3.696177886429626, 3.711988142525477, 3.7276219794450043, 3.7430719098509644, 3.7583305186873055, 3.7733903379316915, 3.788243851316109, 3.8028846954876387, 3.817312100016129, 3.831527054676938, 3.8455307084619545, 3.8593243567873907, 3.8729094755510998, 3.886288941763903, 3.8994678221576633, 3.9124515648981912, 3.925245830659788, 3.9378564791652293, 3.9502893990533496, 3.9625489855862206, 3.9746387801010776, 3.986562450401983, 3.998323785297107, 4.0099267058561, 4.021374045854534, 4.032664849544065, 4.0437974283353695, 4.054770122992369, 4.065581292372619, 4.076229242632671, 4.086709730589422, 4.097015561108973, 4.107139337639362, 4.11707362077688, 4.126810929731954, 4.136343454364696, 4.14566128039544, 4.154753591106116, 4.163609512671855, 4.172218133889753, 4.180568516146076, 4.188649930214281, 4.196452158142277, 4.203965059587182, 4.211178521621482, 4.218082466299532, 4.224666996638198, 4.2309244914351245, 4.236849222541093, 4.242435626795214, 4.247678252708546, 4.252571771150013, 4.2571117750556695, 4.261297646968964, 4.265129973542168, 4.268609494306927, 4.27173710607927, 4.274513881176863, 4.276942957712089, 4.279030624493392, 4.28078357016767, 4.282208602255779, 4.283312655889372, 4.284103027008252, 4.284589771620789, 4.284784581043241, 4.284699218385518, 4.28434549338915, 4.283735242518843, 4.282881196505264, 4.28179885041858, 4.280504211207644, 4.279013241321174, 4.27734184142186, 4.2755058739652, 4.273522047812197, 4.271408022069777, 4.269181392687884, 4.266859635404354, 4.264460099310447, 4.261999922034904, 4.259495563068294, 4.256963124712673, 4.254418580299847, 4.251877777767614, 4.249356433892087, 4.246869408770448, 4.2444298252502675, 4.242050504160675, 4.239744174570696, 4.23752348589533, 4.23540086441478, 4.233386341581095, 4.231488018750844, 4.229713915161976, 4.228072038191698, 4.22657037269768, 4.225216112148949, 4.224012845954393, 4.222963163179267, 4.222069689974626, 4.221335087324935, 4.22076202636901, 4.220351720484131, 4.220103097347491, 4.2200149414453145, 4.220086084190061, 4.220315413308819, 4.220701775683864, 4.221243180368596, 4.221937167897028, 4.22278132569375, 4.22377328255592, 4.224910712273222, 4.226191530455021, 4.227614321468625, 4.229177830460259, 4.230880838590652, 4.232722166780825, 4.234700699166868, 4.2368162158890925, 4.239069489338507, 4.2414613879354555, 4.24399281350242, 4.246664706681919, 4.249478230085728, 4.252435769613542, 4.255540236605479, 4.2587945965414065, 4.262201876516941, 4.265765155994701, 4.269488082881454, 4.273375482515747, 4.277432402676247, 4.281663971837931, 4.286075396819741, 4.290672014050634, 4.295459792907266, 4.30044523310151, 4.305634966651895, 4.311035741309275, 4.316654436734763, 4.322497889830009, 4.328572307236653, 4.334883808827115, 4.341438662961033, 4.3482432911393145, 4.355304242788721, 4.362626905391594, 4.370214758641636, 4.378071244255478, 4.386199930285894, 4.394604531429188, 4.403288614704303, 4.4122530666950075, 4.421497247972388, 4.431020570913924, 4.440822517520494, 4.450902630406415, 4.4612591506958355, 4.471886096705284, 4.482776612175203, 4.49392376510524, 4.505320527944802, 4.516959656334487, 4.528830710315242, 4.540919629289184, 4.553211888653011, 4.565692659127541, 4.578346788188206, 4.591158287432504, 4.6041076251775275, 4.617173512077738, 4.630334147691608, 4.6435672090859015, 4.65684982896757, 4.670158162724894, 4.683466876551761, 4.696749912633001, 4.70998060178722, 4.723131658752282, 4.736175243833715, 4.749083803737985, 4.761829926446765, 4.77438565907723, 4.786722505450079, 4.798811448799843, 4.810623654280443, 4.822132923931623, 4.833313545233968, 4.844139474388925, 4.8545843713955765, 4.864621656522513, 4.874226677879301, 4.883377906455132, 4.892054029886585, 4.900233711285234, 4.907895621216566, 4.915018691877441, 4.921584084178689, 4.9275743215536645, 4.9329721525047585, 4.937760569938032, 4.941922845151281, 4.945442616827566, 4.948304045185453, 4.950491706308134, 4.951990619223776, 4.952786278899506, 4.952864695577834, 4.9522127756513505, 4.950818374942613, 4.948669996991617, 4.945756801343818, 4.942068645765273, 4.93759694653558, 4.932339204990731, 4.9262957972452535, 4.919467837442484, 4.911857173150329, 4.903466409703501, 4.894303191767214, 4.884384961567724, 4.8737308398052654, 4.862360394783141, 4.850293589778906, 4.837551281776449, 4.824162851974283, 4.810164178762128, 4.795591191443758, 4.780479597235405, 4.764864837345608, 4.748783935504041, 4.732281973514573, 4.715405623155855, 4.6982006900930315, 4.680712086509223, 4.662983810280774, 4.645062111537348, 4.626996976690208, 4.608837516610285, 4.590631542134771, 4.572425561857399, 4.554265022586867, 4.536196297675371, 4.518265684321596, 4.500518169499145, 4.482997456036988, 4.465746028114724, 4.448804701217601, 4.432211693603229, 4.416003970326379, 4.40021760114407, 4.384887854178672, 4.370049103468438, 4.355730972909366, 4.341958185982745, 4.328754976307896, 4.316145377179613, 4.304153291163404, 4.2928010697989905, 4.282101473473299, 4.272063751662304, 4.262697536106197, 4.25401288018504, 4.246020277957069, 4.238725150517906, 4.232121817438004, 4.2262038237681105, 4.220965432677376, 4.2164015990907515, 4.212507289725468, 4.209267940369125, 4.206661751411975, 4.2046673275596556, 4.203263845793377, 4.202430992997161, 4.202146774497972, 4.202380370624297, 4.203098723562877, 4.204268952264872, 4.205858294277914, 4.207833976099675, 4.210160828077264, 4.212800034092558, 4.215712278457005, 4.218857976457537, 4.222197214539519, 4.2256898532387765, 4.229296640892457, 4.232978637995626, 4.236696467292306, 4.240410304582089, 4.24407986849314, 4.247666385275297, 4.251136598230013, 4.2544579860742395, 4.257597701179321, 4.260522605787039, 4.263199473737165, 4.265600216776555, 4.26770235064507, 4.269483599281713, 4.270921630320802, 4.27199410986703, 4.272680057864286, 4.272967049008975, 4.272846100505074, 4.272308421866817, 4.271345414615352, 4.269948712814302, 4.268114404331926, 4.26584759171355, 4.263154734634072, 4.260042532983404, 4.256517917898078, 4.252588511850343, 4.2482691136856126, 4.243580049059514, 4.238541913665049, 4.233175371670989, 4.227501138604207, 4.221540748977241, 4.215319042044216, 4.208861709925239, 4.202194351856068, 4.195342435841526, 4.188331272941007, 4.181184424524819, 4.173922952641248, 4.166567582844004, 4.159138873577683, 4.151657215905516, 4.14414217135199, 4.136606994258203, 4.12906147032647, 4.121515259563893, 4.1139779356104205, 4.106458979836944, 4.098965512352464, 4.091497909058593, 4.0840552999004505, 4.0766368323742475, 4.06924165731399, 4.0618688419785425, 4.054515115672026, 4.0471747576359425, 4.039841977716841, 4.032511036227844, 4.025176244545308, 4.0178320587597876, 4.010473635367523, 4.003096417628423, 3.9956958887684104, 3.9882675762826483, 3.980807051250716, 3.97331099742087, 3.96577832652389, 3.9582082292038283, 3.95059990859232, 3.9429525785851607, 3.935265591475744, 3.9275401768226983, 3.919778968745611, 3.9119846233342916, 3.9041597713353235, 3.896307014603368, 3.888429579069324, 3.8805333665368558, 3.872624945613637, 3.864710812897813, 3.8567974082901073, 3.848891142440428, 3.8410005873704334, 3.833137478758211, 3.8253136962335543, 3.817540983921617, 3.8098309901538014, 3.802195573209203, 3.7946496373341816, 3.7872096505664956, 3.779891958073288, 3.7727127662022713, 3.765688177204128, 3.7588342982542176, 3.7521675316204286, 3.74570425673045, 3.7394607929959154, 3.7334534180753076, 3.727698253381694, 3.722207859343246, 3.716991117427279, 3.712056795690564, 3.7074137658143616, 3.7030710063026553, 3.6990360881861752, 3.6953073945828865, 3.6918798927849923, 3.688748707691248, 3.685909113413805, 3.683356486569185, 3.6810821343095084, 3.679068856771425, 3.677298496462221, 3.6757529416325516, 3.674414104024911, 3.673263673658423, 3.6722800102655686, 3.67143891277466, 3.6707160227715634, 3.670086898943093, 3.6695269465719123, 3.6690116286413206, 3.668516911995629, 3.66801881665599, 3.6674932025437514, 3.6669157886545123, 3.666262188422107, 3.665509940899483, 3.6646394606383694, 3.6636313023452227, 3.662465961344933, 3.661123849403062, 3.6595859367524475, 3.657838548464777, 3.6558709097635895, 3.6536723215034574, 3.6512321204291447, 3.648539722493118, 3.64558692171879, 3.6423722702132024, 3.63889562946292, 3.6351569502885677, 3.6311563105212676, 3.626894049850805, 3.6223748423875675, 3.61760785418234, 3.612602527398497, 3.6073683614385232, 3.6019149097885, 3.5962525281054027, 3.5903962959161513, 3.584363056353361, 3.5781696036085466, 3.571832689778545, 3.565368995000886, 3.5587963247657677, 3.552134788283002, 3.5454046642592165, 3.538626089219678, 3.5318190560554275, 3.525003359376313, 3.5181979043912075, 3.511420884734179, 3.504690354536419, 3.4980242389596494, 3.4914403229848388, 3.4849553014327723, 3.4785816330889485, 3.47233064909504, 3.466213630444083, 3.46024182525684, 3.4544263753409963, 3.4487758109447983, 3.4432948788081976, 3.4379880954882087, 3.4328599813153637, 3.4279151124256786, 3.423157924159992, 3.4185913025723544, 3.414217342531489, 3.4100381953214636, 3.406056073217442, 3.4022732559760827, 3.398692127577926, 3.395315190409169, 3.392145050552847, 3.389184383787298, 3.3864359594064206, 3.3839026044556655, 3.3815863384374207, 3.379488367126825, 3.3776099412715372, 3.3759524162540675, 3.374517237660656, 3.373305248550887, 3.3723130990717434, 3.3715359025243328, 3.3709688265713007, 3.370607082976864, 3.3704459049758886, 3.370477602446871, 3.3706885969486664, 3.3710645252659894, 3.3715909392117505, 3.372253279229991, 3.373036525069204, 3.373920894576579, 3.3748829387690695, 3.375898875680086, 3.376944640755931, 3.3779958964951873, 3.379027206728619, 3.3800096366276917, 3.380913060601067, 3.381707018074029, 3.382360702565492, 3.3828430150791577, 3.383123495594029, 3.3831727841003, 3.3829613923093502, 3.3824596727530833, 3.3816378575741926, 3.3804668140690453, 3.378923902580492, 3.3769899883101284, 3.374646111347684, 3.371873491367995, 3.3686535824626977, 3.3649717760593867, 3.3608241079380132, 3.3562087945850183, 3.3511244145675882, 3.345569923319486, 3.339544923414551, 3.3330558182169496, 3.3261158824378136, 3.3187389638913136, 3.3109391655006943, 3.3027308017166934, 3.2941291292575334, 3.2851539884624343, 3.2758269180850346, 3.2661694698387054, 3.2562031428345213, 3.2459493648930096, 3.235430077095297, 3.2246682949960155, 3.213686970595221, 3.202508802664294, 3.1911562369161226, 3.179651343942523, 3.1680148239558674, 3.1562663045597987, 3.1444250834717207, 3.132510140051296, 3.1205401208426164, 3.108532549417599, 3.0965014156661748, 3.0844596259339654, 3.0724198516054697, 3.0603945095411795, 3.048395718952562, 3.036432822926277, 3.0245114321521918, 3.012636780769486, 3.0008139908905953, 2.9890480696742276, 2.9773434685839417, 2.9657003732314755, 2.954116824701327, 2.942590867338216, 2.931120572292994, 2.9197040363942746, 2.9083383698381895, 2.897017835784297, 2.885736327418666, 2.8744878215295953, 2.863266409222878, 2.852066234668112, 2.8408803825981823, 2.829700906158505, 2.818519953146419, 2.8073298010247436, 2.7961228766533197, 2.784891470776775, 2.7736263201602345, 2.7623177506311354, 2.7509562669833096, 2.739532574862049, 2.7280375882173513, 2.7164614692840723, 2.704792818646009, 2.69302028709162, 2.681132819006826, 2.669119643917352, 2.6569702241887327, 2.644673493794231, 2.6322181119324317, 2.6195931017158034, 2.606787860083055, 2.5937921882853257, 2.580596737914703, 2.5671943824888626, 2.553578813147443, 2.539744108260125, 2.5256847365787083, 2.511395638602193, 2.4968755486312046, 2.4821281670117896, 2.467157698205591, 2.4519685084670515, 2.4365650894268733, 2.4209528763130232, 2.4051441639231186, 2.3891545568648005, 2.372999394792389, 2.3566936860209298, 2.340252113288651, 2.323691106047454, 2.3070323372875583, 2.290297751905915, 2.2735085642848896, 2.25668522077732, 2.239847500809759, 2.223016153769243, 2.2062128051645944, 2.1894582554001008, 2.17277243876032, 2.1561743998753076, 2.1396821283461787, 2.123311431103474, 2.107076912191324, 2.0909924892282303, 2.0750714304040514, 2.0593263552508576, 2.0437671126895585, 2.0283989088206438, 2.013226199668328, 1.998253198365033, 1.9834838840529199, 1.9689216679814137, 1.9545649136570635, 1.9404085922294954, 1.9264477195307828, 1.9126774791802705, 1.8990931881339197, 1.885689139203471, 1.8724552827620093, 1.8593807850471626, 1.8464551568046692, 1.8336682416231491, 1.821010183677138, 1.808470555364633, 1.7960380629683634, 1.7837016753731627, 1.7714506995012653, 1.7592747749942055, 1.7471640347355262, 1.7351104190696303, 1.7231068800131073, 1.7111465808253263, 1.6992229047549121, 1.6873294399134386, 1.6754611421017482, 1.6636162927023843, 1.6517937468304658, 1.6399924004841633, 1.6282111774180783, 1.6164491748823249, 1.6047086425397874, 1.5929948238855949, 1.581312897465798, 1.5696678555197037, 1.558064509953889, 1.54650823653928, 1.5350085040146293, 1.523575948601564, 1.5122208258081347, 1.5009530034191993, 1.489782005157829, 1.4787183440986458, 1.4677746817693156, 1.4569635339811675, 1.4462969927761264, 1.435786744435658, 1.4254441065218229, 1.415279873556771, 1.405304424724414, 1.3955278709672447, 1.3859601059570954, 1.3766108365003298, 1.3674887965345148, 1.3585993799001073, 1.3499472155592656, 1.3415369968055102, 1.3333735152467552, 1.3254615992717926, 1.3178030141267305, 1.310395360766446, 1.3032362379413533, 1.296323530660785, 1.2896554490221155, 1.2832299146194872, 1.2770401962992242, 1.2710774135129344, 1.2653329814506336, 1.2597986304302178, 1.2544663579921145, 1.2493272708192753, 1.2443696330390839, 1.2395813636833939, 1.2349505520299644, 1.2304654081059365, 1.2261142061956851, 1.221884954065583, 1.217765328015086, 1.213743008100218, 1.2098056180262418, 1.2059407330813694, 1.2021360221446638, 1.198380088920987, 1.194661798267481, 1.1909698820747778, 1.1872929364414166, 1.1836193983955778, 1.1799386815133344, 1.1762422099462868, 1.172521508113262, 1.1687679542884222, 1.1649727781152626, 1.1611272448019139, 1.1572251021866482, 1.1532618364726843, 1.1492328747139742, 1.1451335642621796, 1.14095915144213, 1.136705383472132, 1.132370221048148, 1.127952150264261, 1.123449587386786, 1.118860937779801, 1.1141845480366217, 1.1094187052297881, 1.1045616726358412, 1.099611654263532, 1.0945667961264856, 1.0894251883821717, 1.0841848775939018, 1.078843840871176, 1.073400001526169, 1.0678512085960017, 1.0621952532334442, 1.0564298567866832, 1.0505526566885863, 1.0445612271552487, 1.0384530540113783, 1.0322255417274875, 1.025876007063647]

		spl_ea = fit(SmoothingSpline, ea_lamA, ea_cm2, 1.0)
		ea = SmoothingSplines.predict(spl_ea, l)
 
=#
    #Update effective area on 21 August 2022
	#=
	(c0,c1,c2,c3,c4,c5, c6, c7)= (-1203075.293157474, 5642.681983417952, -11.301034138859688, 0.012529016935390545, -8.304935026402432e-06, 3.2916444437952605e-09, -7.223770113658085e-13, 6.772151844554497e-17)
	ea =  c0 + c1 * l  + c2 * l^2 + c3*l^3  + c4* l^4 + c5 * l^5 + c6 *l^6 + c7*l^7
  =#

  #= Updated effective area based on WD0308 21 Aug 2022.

   ea_lamA = [1240.86414133, 1243.65379599, 1246.44345065, 1249.23310532,
   1252.02275998, 1254.81241465, 1257.60206931, 1260.39172398,
   1263.18137864, 1265.97103331, 1268.76068797, 1271.55034264,
   1274.3399973 , 1277.12965197, 1279.91930663, 1282.7089613 ,
   1285.49861596, 1288.28827063, 1291.07792529, 1293.86757996,
   1296.65723462, 1299.44688929, 1302.23654395, 1305.02619861,
   1307.81585328, 1310.60550794, 1313.39516261, 1316.18481727,
   1318.97447194, 1321.7641266 , 1324.55378127, 1327.34343593,
   1330.1330906 , 1332.92274526, 1335.71239993, 1338.50205459,
   1341.29170926, 1344.08136392, 1346.87101859, 1349.66067325,
   1352.45032792, 1355.23998258, 1358.02963724, 1360.81929191,
   1363.60894657, 1366.39860124, 1369.1882559 , 1371.97791057,
   1374.76756523, 1377.5572199 , 1380.34687456, 1383.13652923,
   1385.92618389, 1388.71583856, 1391.50549322, 1394.29514789,
   1397.08480255, 1399.87445722, 1402.66411188, 1405.45376655,
   1408.24342121, 1411.03307588, 1413.82273054, 1416.6123852 ,
   1419.40203987, 1422.19169453, 1424.9813492 , 1427.77100386,
   1430.56065853, 1433.35031319, 1436.13996786, 1438.92962252,
   1441.71927719, 1444.50893185, 1447.29858652, 1450.08824118,
   1452.87789585, 1455.66755051, 1458.45720518, 1461.24685984,
   1464.03651451, 1466.82616917, 1469.61582384, 1472.4054785 ,
   1475.19513316, 1477.98478783, 1480.77444249, 1483.56409716,
   1486.35375182, 1489.14340649, 1491.93306115, 1494.72271582,
   1497.51237048, 1500.30202515, 1503.09167981, 1505.88133448,
   1508.67098914, 1511.46064381, 1514.25029847, 1517.03995314,
   1519.8296078 , 1522.61926247, 1525.40891713, 1528.19857179,
   1530.98822646, 1533.77788112, 1536.56753579, 1539.35719045,
   1542.14684512, 1544.93649978, 1547.72615445, 1550.51580911,
   1553.30546378, 1556.09511844, 1558.88477311, 1561.67442777,
   1564.46408244, 1567.2537371 , 1570.04339177, 1572.83304643,
   1575.6227011 , 1578.41235576, 1581.20201043, 1583.99166509,
   1586.78131975, 1589.57097442, 1592.36062908, 1595.15028375,
   1597.93993841, 1600.72959308, 1603.51924774, 1606.30890241,
   1609.09855707, 1611.88821174, 1614.6778664 , 1617.46752107,
   1620.25717573, 1623.0468304 , 1625.83648506, 1628.62613973,
   1631.41579439, 1634.20544906, 1636.99510372, 1639.78475839,
   1642.57441305, 1645.36406771, 1648.15372238, 1650.94337704,
   1653.73303171, 1656.52268637, 1659.31234104, 1662.1019957 ,
   1664.89165037, 1667.68130503, 1670.4709597 , 1673.26061436,
   1676.05026903, 1678.83992369, 1681.62957836, 1684.41923302,
   1687.20888769, 1689.99854235, 1692.78819702, 1695.57785168,
   1698.36750634, 1701.15716101, 1703.94681567, 1706.73647034,
   1709.526125  , 1712.31577967, 1715.10543433, 1717.895089  ,
   1720.68474366, 1723.47439833, 1726.26405299, 1729.05370766,
   1731.84336232, 1734.63301699, 1737.42267165, 1740.21232632,
   1743.00198098, 1745.79163565, 1748.58129031, 1751.37094498,
   1754.16059964, 1756.9502543 , 1759.73990897, 1762.52956363,
   1765.3192183 , 1768.10887296, 1770.89852763, 1773.68818229,
   1776.47783696, 1779.26749162, 1782.05714629, 1784.84680095,
   1787.63645562, 1790.42611028, 1793.21576495]

   ea_cm2 = [0.21043198, 0.21412498, 0.21971267, 0.22715562, 0.23642387,
   0.24749582, 0.26035693, 0.27499881, 0.29141805, 0.30961542,
   0.32959484, 0.35136261, 0.37492668, 0.40029586, 0.42747923,
   0.4564855 , 0.4873224 , 0.51999626, 0.55451157, 0.59087042,
   0.62907228, 0.66911361, 0.71098752, 0.75468364, 0.80018778,
   0.84748183, 0.89654359, 0.94734667, 0.99986039, 1.05404974,
   1.10987533, 1.16729344, 1.22625607, 1.28671081, 1.34860109,
   1.41186632, 1.47644171, 1.54225878, 1.60924521, 1.67732513,
   1.74641929, 1.81644524, 1.88731746, 1.95894772, 2.0312452 ,
   2.10411671, 2.17746703, 2.25119895, 2.32521373, 2.39941133,
   2.47369054, 2.54794936, 2.62208517, 2.69599508, 2.76957621,
   2.84272579, 2.91534157, 2.98732209, 3.05856691, 3.12897669,
   3.1984538 , 3.26690212, 3.33422774, 3.40033881, 3.46514595,
   3.5285625 , 3.59050454, 3.6508913 , 3.70964527, 3.7666923 ,
   3.82196193, 3.87538739, 3.92690589, 3.97645866, 4.02399114,
   4.06945308, 4.11279858, 4.15398634, 4.19297958, 4.22974616,
   4.26425869, 4.29649448, 4.32643563, 4.35406913, 4.37938665,
   4.40238468, 4.4230645 , 4.44143211, 4.45749819, 4.47127805,
   4.48279157, 4.49206299, 4.49912104, 4.50399866, 4.50673282,
   4.50736449, 4.50593873, 4.50250389, 4.49711213, 4.48981879,
   4.48068233, 4.4697642 , 4.45712867, 4.44284225, 4.42697396,
   4.4095947 , 4.39077716, 4.37059577, 4.3491259 , 4.32644421,
   4.30262799, 4.27775497, 4.25190306, 4.22515016, 4.19757363,
   4.16925036, 4.14025609, 4.1106654 , 4.08055129, 4.04998505,
   4.01903574, 3.98776998, 3.95625196, 3.92454271, 3.89270006,
   3.86077867, 3.82882908, 3.79689821, 3.76502851, 3.73325832,
   3.70162111, 3.67014579, 3.63885619, 3.60777086, 3.57690334,
   3.54626155, 3.51584803, 3.48565966, 3.45568762, 3.42591766,
   3.39632941, 3.36689714, 3.33758938, 3.30836916, 3.27919395,
   3.25001597, 3.22078223, 3.19143496, 3.16191148, 3.13214497,
   3.10206438, 3.07159508, 3.04065891, 3.00917524, 2.97706094,
   2.94423117, 2.91060006, 2.8760811 , 2.84058809, 2.80403607,
   2.76634171, 2.72742457, 2.68720803, 2.64562007, 2.60259468,
   2.5580729 , 2.51200394, 2.46434667, 2.41507092, 2.36415893,
   2.31160702, 2.25742709, 2.20164832, 2.14431887, 2.08550839,
   2.02530898, 1.96383822, 1.90124077, 1.83769074, 1.77339424,
   1.70859168, 1.64356054, 1.57861766, 1.51412267, 1.45048031,
   1.38814375, 1.32761755, 1.26946149, 1.21429315, 1.16279199,
   1.1157029 , 1.07383977, 1.03808972, 1.00941702, 0.98886705,
   0.97757141, 0.97675128, 0.98772302, 1.01190256]
=#

#= 
Update effective area based on WD0308 data after rejecting cosmic-ray events

ea_lamA = [1240.86414133, 1243.65379599, 1246.44345065, 1249.23310532,
1252.02275998, 1254.81241465, 1257.60206931, 1260.39172398,
1263.18137864, 1265.97103331, 1268.76068797, 1271.55034264,
1274.3399973 , 1277.12965197, 1279.91930663, 1282.7089613 ,
1285.49861596, 1288.28827063, 1291.07792529, 1293.86757996,
1296.65723462, 1299.44688929, 1302.23654395, 1305.02619861,
1307.81585328, 1310.60550794, 1313.39516261, 1316.18481727,
1318.97447194, 1321.7641266 , 1324.55378127, 1327.34343593,
1330.1330906 , 1332.92274526, 1335.71239993, 1338.50205459,
1341.29170926, 1344.08136392, 1346.87101859, 1349.66067325,
1352.45032792, 1355.23998258, 1358.02963724, 1360.81929191,
1363.60894657, 1366.39860124, 1369.1882559 , 1371.97791057,
1374.76756523, 1377.5572199 , 1380.34687456, 1383.13652923,
1385.92618389, 1388.71583856, 1391.50549322, 1394.29514789,
1397.08480255, 1399.87445722, 1402.66411188, 1405.45376655,
1408.24342121, 1411.03307588, 1413.82273054, 1416.6123852 ,
1419.40203987, 1422.19169453, 1424.9813492 , 1427.77100386,
1430.56065853, 1433.35031319, 1436.13996786, 1438.92962252,
1441.71927719, 1444.50893185, 1447.29858652, 1450.08824118,
1452.87789585, 1455.66755051, 1458.45720518, 1461.24685984,
1464.03651451, 1466.82616917, 1469.61582384, 1472.4054785 ,
1475.19513316, 1477.98478783, 1480.77444249, 1483.56409716,
1486.35375182, 1489.14340649, 1491.93306115, 1494.72271582,
1497.51237048, 1500.30202515, 1503.09167981, 1505.88133448,
1508.67098914, 1511.46064381, 1514.25029847, 1517.03995314,
1519.8296078 , 1522.61926247, 1525.40891713, 1528.19857179,
1530.98822646, 1533.77788112, 1536.56753579, 1539.35719045,
1542.14684512, 1544.93649978, 1547.72615445, 1550.51580911,
1553.30546378, 1556.09511844, 1558.88477311, 1561.67442777,
1564.46408244, 1567.2537371 , 1570.04339177, 1572.83304643,
1575.6227011 , 1578.41235576, 1581.20201043, 1583.99166509,
1586.78131975, 1589.57097442, 1592.36062908, 1595.15028375,
1597.93993841, 1600.72959308, 1603.51924774, 1606.30890241,
1609.09855707, 1611.88821174, 1614.6778664 , 1617.46752107,
1620.25717573, 1623.0468304 , 1625.83648506, 1628.62613973,
1631.41579439, 1634.20544906, 1636.99510372, 1639.78475839,
1642.57441305, 1645.36406771, 1648.15372238, 1650.94337704,
1653.73303171, 1656.52268637, 1659.31234104, 1662.1019957 ,
1664.89165037, 1667.68130503, 1670.4709597 , 1673.26061436,
1676.05026903, 1678.83992369, 1681.62957836, 1684.41923302,
1687.20888769, 1689.99854235, 1692.78819702, 1695.57785168,
1698.36750634, 1701.15716101, 1703.94681567, 1706.73647034,
1709.526125  , 1712.31577967, 1715.10543433, 1717.895089  ,
1720.68474366, 1723.47439833, 1726.26405299, 1729.05370766,
1731.84336232, 1734.63301699, 1737.42267165, 1740.21232632,
1743.00198098, 1745.79163565, 1748.58129031, 1751.37094498,
1754.16059964, 1756.9502543 , 1759.73990897, 1762.52956363,
1765.3192183 , 1768.10887296, 1770.89852763, 1773.68818229,
1776.47783696, 1779.26749162, 1782.05714629, 1784.84680095,
1787.63645562, 1790.42611028, 1793.21576495]

ea_cm2 = [0.075906  , 0.10787247, 0.13749315, 0.16509385, 0.19098686,
0.2154708 , 0.23883057, 0.26133732, 0.28324834, 0.30480724,
0.32624382, 0.34777425, 0.36960107, 0.39191335, 0.41488685,
0.4386841 , 0.46345463, 0.48933515, 0.51644974, 0.54491013,
0.57481589, 0.60625472, 0.63930271, 0.67402464, 0.71047429,
0.7486947 , 0.78871855, 0.83056844, 0.87425734, 0.91978871,
0.96715716, 1.01634852, 1.06734042, 1.12010253, 1.17459701,
1.23077886, 1.28859631, 1.34799117, 1.40889926, 1.47125078,
1.5349707 , 1.59997912, 1.66619163, 1.73351982, 1.80187147,
1.87115103, 1.94126009, 2.0120975 , 2.08355991, 2.15554215,
2.2279374 , 2.30063777, 2.37353441, 2.44651804, 2.51947912,
2.59230827, 2.66489648, 2.73713553, 2.80891817, 2.88013846,
2.95069201, 3.02047624, 3.08939069, 3.15733719, 3.2242201 ,
3.28994655, 3.35442666, 3.41757366, 3.4793042 , 3.53953842,
3.59820011, 3.65521694, 3.71052047, 3.7640465 , 3.81573488,
3.86552983, 3.91337997, 3.95923833, 4.00306256, 4.04481481,
4.08446189, 4.12197522, 4.15733093, 4.19050977, 4.22149709,
4.250283  , 4.27686208, 4.30123348, 4.32340087, 4.34337234,
4.36116025, 4.3767813 , 4.39025627, 4.40160999, 4.4108712 ,
4.41807243, 4.42324982, 4.42644303, 4.42769502, 4.42705188,
4.42456271, 4.42027944, 4.41425648, 4.4065508 , 4.39722143,
4.38632942, 4.37393765, 4.36011039, 4.34491339, 4.32841335,
4.31067789, 4.29177515, 4.27177374, 4.25074227, 4.22874929,
4.20586298, 4.18215085, 4.15767955, 4.13251465, 4.10672034,
4.08035922, 4.05349201, 4.02617742, 3.99847186, 3.97042919,
3.94210052, 3.91353411, 3.88477493, 3.85586472, 3.82684169,
3.79774024, 3.76859107, 3.73942081, 3.71025187, 3.68110254,
3.65198661, 3.62291352, 3.59388803, 3.56491055, 3.53597657,
3.50707717, 3.47819868, 3.44932289, 3.42042695, 3.39148361,
3.36246117, 3.33332374, 3.30403126, 3.27453986, 3.24480185,
3.21476613, 3.18437851, 3.15358184, 3.12231645, 3.09052067,
3.05813097, 3.02508261, 2.99131023, 2.95674817, 2.92133115,
2.88499493, 2.847677  , 2.80931723, 2.76985871, 2.72924851,
2.68743869, 2.64438705, 2.60005822, 2.55442469, 2.50746782,
2.45917914, 2.4095615 , 2.35863036, 2.30641499, 2.25296016,
2.19832741, 2.14259658, 2.0858676 , 2.02826198, 1.96992471,
1.91102603, 1.85176342, 1.7923634 , 1.7330839 , 1.67421614,
1.61608702, 1.55906149, 1.50354484, 1.44998515, 1.39887616,
1.35075977, 1.30622885, 1.26592999, 1.23056689, 1.2009029 ,
1.17776469, 1.16204522, 1.15470746, 1.15678744]

=#

# Read effective area based on WD0308_crreject data

f = FITS("./caldata/fuv_grating2m2_effarea_12nov22.fits")
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
    fuv_grating2_fluxed_spec(target,fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile[,order = -2, cross_disp_width_pixels= 50])

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
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
## Output
- `(λ, f_λ, err_f_λ)`
- Fluxed spectrum saved in an ascii file.
...
"""
function fuv_grating2_fluxed_spec(target::String, fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order=-2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels::Int=50)
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
    fuv_grating2_phafile(target,fuv_grating2_image_file, ds9srcregfile, ds9bgdregfile[,order = -2, cross_disp_width_pixels= 50])

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
- `cross_disp_width_pixels::String`: 50 (default), width in pixels in the cross-dispersion direction.
- `respdir::String`: Name of the directory containing the response matrices in the local machine. The response files can be downloaded from the GitHub page.
#read necessary keywords`
## Output
- Source and background PHA files.
- Some relevant information are also printed on the screen.
...
"""
function fuv_grating2_phafile(target::String,fuv_grating2_image_file::String,ds9srcregfile::String, ds9bgdregfile::String; order=-2, disp_aligned_to_xaxis::Bool=false, angle_xaxis_disp_deg::Float64=267.479, cross_disp_width_pixels=50, respdir::String="/soft/astrosat/responses/uvit/")
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

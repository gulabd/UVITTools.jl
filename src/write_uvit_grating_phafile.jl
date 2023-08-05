
#using FITSIO
"""
   `write_uvit_grating_phafile(uvit_detector, uvit_grating, channels, counts, exptime_sec[, phafile="phafile.pha"])`

Write OGIP compatible PHA spectral file from grating count spectra - channel vs counts.

This function is used by grating spectroscopy tools.
...
# Arguments
- `uvit_detector::String`: FUV or NUV.
- `uvit_grating::String`: One of the UVIT gratings - FUV-Grating1, FUV-Grating2, NUV-Grating.
- `channels::Array{Int64}`: An array of spectral channels.
- `counts::Array{Float64}`: An array of counts corresponding to the spectral channel array.
- `exposure_time_sec::Float64`: Exposure time (s).
# Optional
- `phafile::String`: Name of the PHA spectral file to be generated.
...

"""
function write_uvit_grating_phafile(uvit_detector::String, uvit_grating::String,channels::Array{Int64,1}, counts::Array{Float64,1}, exptime_sec::Float64; phafile="phafile.pha")


# Define FITS pha keyword

	kys=String["XTENSION", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS", "TTYPE1", "TFORM1", "TTYPE2", "TFORM2", "TUNIT2", "TTYPE3", "TFORM3", "TTYPE4", "TFORM4", "EXTNAME", "HDUCLASS", "HDUCLAS1", "HDUVERS1", "HDUVERS", "HDUCLAS2", "HDUCLAS3", "TLMIN1", "TLMAX1", "TELESCOP", "INSTRUME", "FILTER", "EXPOSURE", "AREASCAL", "BACKFILE", "BACKSCAL", "CORRFILE", "CORRSCAL", "RESPFILE", "ANCRFILE", "PHAVERSN", "DETCHANS", "CHANTYPE", "POISSERR", "STAT_ERR", "SYS_ERR"]

	vals=Any["BINTABLE", Int8(8), Int8(2), Int8(10), Int8(1), Int8(0), Int8(1), Int8(4), "CHANNEL", "I", "COUNTS", "J", "count", "QUALITY", "I", "GROUPING", "I", "SPECTRUM", "OGIP", "SPECTRUM", "1.1.0", "1.1.0", " ", "COUNT", Int8(0), Int8(0), "ASTROSAT/UVIT", "UVIT", " ", 0.0, 1.0, "NONE", 1.0, "NONE", -1.0, "", "", "1992a", Int8(1), "PI", true, Int8(0), Int8(0)]

	comments=String["binary table extension", "8-bit bytes", "2-dimensional binary table", "width of table in bytes", "number of rows in table", "size of special data area", "one data group (required keyword)", "number of fields in each row", "Pulse Invariant (PI) Channel", "data format of field: 2-byte INTEGER", "Counts per channel", "data format of field: 4-byte INTEGER", "physical unit of field", "Quality flag of this channel (0=good)", "data format of field: 2-byte INTEGER", "Grouping flag for channel (0=undefined)", "data format of field: 2-byte INTEGER", "name of this binary table ex", "format conforms to OGIP standard", "PHA dataset (OGIP memo OGIP-92-007)", "Obsolete - included for backwards compatibility", "Version of format (OGIP memo OGIP-92-007)", "WARNING This is NOT an OGIP-approved value", "PHA data stored as Counts (not count/s)", "Lowest legal channel number", "Highest legal channel number", "Telescope (mission) name", "Instrument name", "Instrument filter in use", " Exposure time", " Area scaling factor", "Background FITS file", " Background scale factor", "Correlation FITS file", " Correlation scale factor", "redistribution matrix", "ancillary response", "obsolete", "total number possible channels", "channel type (PHA, PI etc)", "Poisson errors to be assumed", "no statistical error specified", "no systematic error specified"]


	pha_header = FITSHeader(kys,vals,comments)
	pha_header["EXPOSURE"]=exptime_sec
	pha_header["DETECTOR"]=uvit_detector
	pha_header["FILTER"]=uvit_grating
# Set number of valid detector channels
	pha_header["DETCHANS"]=length(channels)


# Open new file to write
	pha = FITS(phafile,"w")

# Convert to same data type as used by FTOOLS
	channels = convert.(Int32,round.(channels))	
	counts = convert.(Int32,round.(counts))

	data = Dict("CHANNEL"=>channels,"COUNTS"=>counts,"QUALITY"=>convert.(Int16, zeros(length(channels))),"GROUPING"=>convert.(Int16, ones(length(channels))))
#	println(data)
	write(pha,data,hdutype=TableHDU, name="SPECTRUM", header=pha_header)
	close(pha)

	return phafile
end

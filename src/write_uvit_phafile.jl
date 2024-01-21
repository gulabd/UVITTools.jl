#= This function writes OGIP pha file for a UVIT filter. The inputs are filter name either HST style name or original name, number of counts, exposure time in seconds. The parameters phafile and respdir are optional.
=#

#Date : 13 Feb 2018, G. C. Dewangan
#
# First Update:  18 October 2020 
#				Included new response files based on Tandon et al. 2020 calibration
#       Earlier rmf/arf files replaced with .rsp file that includes both rmf/arf

#using FITSIO

function write_uvit_phafile(uvit_detector::String, uvit_filter::String,counts, exptime_sec::Float64; phafile::String="phafile.pha")


# Define FITS pha keyword

	kys=String["XTENSION", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS", "TTYPE1", "TFORM1", "TTYPE2", "TFORM2", "TUNIT2", "TTYPE3", "TFORM3", "TTYPE4", "TFORM4", "EXTNAME", "HDUCLASS", "HDUCLAS1", "HDUVERS1", "HDUVERS", "HDUCLAS2", "HDUCLAS3", "TLMIN1", "TLMAX1", "TELESCOP", "INSTRUME", "FILTER", "EXPOSURE", "AREASCAL", "BACKFILE", "BACKSCAL", "CORRFILE", "CORRSCAL", "RESPFILE", "ANCRFILE", "PHAVERSN", "DETCHANS", "CHANTYPE", "POISSERR", "STAT_ERR", "SYS_ERR"]

	vals=Any["BINTABLE", Int8(8), Int8(2), Int8(10), Int8(1), Int8(0), Int8(1), Int8(4), "CHANNEL", "I", "COUNTS", "J", "count", "QUALITY", "I", "GROUPING", "I", "SPECTRUM", "OGIP", "SPECTRUM", "1.1.0", "1.1.0", " ", "COUNT", Int8(0), Int8(0), "ASTROSAT", "UVIT", " ", 0.0, 1.0, "NONE", 1.0, "NONE", -1.0, "", "", "1992a", Int8(1), "PI", true, Int8(0), Int8(0)]

	comments=String["binary table extension", "8-bit bytes", "2-dimensional binary table", "width of table in bytes", "number of rows in table", "size of special data area", "one data group (required keyword)", "number of fields in each row", "Pulse Invariant (PI) Channel", "data format of field: 2-byte INTEGER", "Counts per channel", "data format of field: 4-byte INTEGER", "physical unit of field", "Quality flag of this channel (0=good)", "data format of field: 2-byte INTEGER", "Grouping flag for channel (0=undefined)", "data format of field: 2-byte INTEGER", "name of this binary table ex", "format conforms to OGIP standard", "PHA dataset (OGIP memo OGIP-92-007)", "Obsolete - included for backwards compatibility", "Version of format (OGIP memo OGIP-92-007)", "WARNING This is NOT an OGIP-approved value", "PHA data stored as Counts (not count/s)", "Lowest legal channel number", "Highest legal channel number", "Telescope (mission) name", "Instrument name", "Instrument filter in use", " Exposure time", " Area scaling factor", "Background FITS file", " Background scale factor", "Correlation FITS file", " Correlation scale factor", "redistribution matrix", "ancillary response", "obsolete", "total number possible channels", "channel type (PHA, PI etc)", "Poisson errors to be assumed", "no statistical error specified", "no systematic error specified"]


	pha_header = FITSHeader(kys,vals,comments)
	pha_header["EXPOSURE"]=exptime_sec
	pha_header["DETECTOR"]=uvit_detector
	pha_header["FILTER"]=uvit_filter
#	pha_header["DATE"]=Dates.now(Dates.UTC)
#	set_comment!(pha_header, "DATE", "file creation date (YYYY-MM-DDThh:mm:ss UTC")

# Open new file to write
	pha = FITS(phafile,"w")
	data = Dict("CHANNEL"=>convert.(Int32,[1]),"COUNTS"=>[counts],"QUALITY"=>convert.(Int16,[0]),"GROUPING"=>convert.(Int16,[1]))
	println(data)
	write(pha,data;header=pha_header,hdutype=TableHDU, name="SPECTRUM")
	close(pha)
	return phafile
end

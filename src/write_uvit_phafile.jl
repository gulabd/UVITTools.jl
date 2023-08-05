#= This function writes OGIP pha file for a UVIT filter. The inputs are filter name either HST style name or original name, number of counts, exposure time in seconds. The parameters phafile and respdir are optional.
=#

#Date : 13 Feb 2018, G. C. Dewangan
#
# First Update:  18 October 2020 
#				Included new response files based on Tandon et al. 2020 calibration
#       Earlier rmf/arf files replaced with .rsp file that includes both rmf/arf

#using FITSIO

function write_uvit_phafile(uvit_detector::String, uvit_filter::String,counts, exptime_sec::Float64; phafile::String="phafile.pha",respdir::String="/soft/astrosat/responses/uvit/")


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
# Update rmf/arf files in the pha header
   
		if uvit_filter=="F1" || uvit_filter=="F148W" || uvit_filter=="CaF2-1"
			respfile="F148W_effarea_Tandon_etal2020.rsp"
   # 	rmffile="uvit_fuv_caf21.rmf"
   # 	arffile="uvit_fuv_caf21_updated.arf"

		elseif uvit_filter== "F2" || uvit_filter=="F154W" || uvit_filter=="BaF2"
			respfile="F154W_effarea_Tandon_etal2020.rsp"
    #	rmffile="uvit_fuv_baf2.rmf"
    #	arffile="uvit_fuv_baf2_updated.arf"
		elseif uvit_filter=="F3" || uvit_filter=="F169M" || uvit_filter=="Sapphire"
			respfile="F169M_effarea_Tandon_etal2020.rsp"
   # 	rmffile="uvit_fuv_saph.rmf"
   # 	arffile="uvit_fuv_saph_updated.arf"
		elseif uvit_filter=="F5" || uvit_filter=="F172M" || uvit_filter=="Silica"
			respfile="F172M_effarea_Tandon_etal2020.rsp"
   # 	rmffile="uvit_fuv_sil.rmf"
   # 	arffile="uvit_fuv_sil_updated.arf"
		elseif uvit_filter=="F7" || uvit_filter=="F148Wa" || uvit_filter=="CaF2-2"
			respfile="NONE"
    #	rmffile="uvit_fuv_caf22.rmf"
    #	arffile="uvit_fuv_caf22_updated.arf"
    

		elseif uvit_filter=="F1" || uvit_filter=="N242W" || uvit_filter=="Silica-1"
			respfile="N242W_effarea_Tandon_etal2020.rsp"
    #	rmffile="NONE"
    #	arffile="NONE"
		elseif uvit_filter=="F2" || uvit_filter=="N219M" || uvit_filter=="NUVB15"
			respfile="N219M_effarea_Tandon_etal2020.rsp"
    #	rmffile="uvit_nuv_b15.rmf"
    #	arffile="uvit_nuv_b15_updated.arf"
   
		elseif uvit_filter=="F3" || uvit_filter=="N245M" || uvit_filter=="NUVB13"
			respfile="N245M_effarea_Tandon_etal2020.rsp"
    #	rmffile="uvit_nuv_b13.rmf"
    #	arffile="uvit_nuv_b13_updated.arf"
    

		elseif uvit_filter=="F5" || uvit_filter=="N263M" || uvit_filter=="NUVB4"
			respfile="N263M_effarea_Tandon_etal2020.rsp"
	#	rmffile="uvit_nuv_b4.rmf"
    #	arffile="uvit_nuv_b4_updated.arf"
    
		elseif uvit_filter=="F6" || uvit_filter=="N279N" || uvit_filter=="NUVN2"
			respfile="N279N_effarea_Tandon_etal2020.rsp"
	#	rmffile="uvit_nuv_n2.rmf"
   # 	arffile="uvit_nuv_n2_updated.arf"
    elseif uvit_filter=="F7" || uvit_filter=="N242Wa" || uvit_filter=="Silica-2"
		rmffile="NONE"
    	arffile="NONE"
    else 
    	print("Filter name not recognised, see http://uvit.iiap.res.in/Instrument/Filters")
    	print("rmf/arf filenames not updated in the PHA header.")
    	rmffile="NONE"
    	arffile="NONE"
    end
	respdir = joinpath(dirname(dirname(pathof(UVITTools))), "caldata")
    if rmffile=="NONE" 
    	pha_header["RESPFILE"]=" "
    #	pha_header["ANCRFILE"]=" "
    else
    	pha_header["RESPFILE"]=respdir * respfile
   # 	pha_header["ANCRFILE"]=respdir * arffile
    end
    pha_header

# Open new file to write
	pha = FITS(phafile,"w")
	data = Dict("CHANNEL"=>convert.(Int32,[1]),"COUNTS"=>[counts],"QUALITY"=>convert.(Int16,[0]),"GROUPING"=>convert.(Int16,[1]))
	println(data)
	write(pha,data;header=pha_header,hdutype=TableHDU, name="SPECTRUM")
	close(pha)
	return phafile
end

"""
    xycen_from_ds9reg(ds9regfile)

Extract center coordinates (x,y) from a DS9 circular region file.
"""
function xycen_from_ds9reg(ds9regfile::String)
	freg=open(ds9regfile)
	freg_content=readlines(freg)[4]
	src_region=split(split(freg_content,"(")[2], ")")[1]

#Extract source x/y center 
	src_cenx_string = split(src_region, ",")[1]
	src_ceny_string = split(src_region, ",")[2]

	src_cenx = parse(Float64,src_cenx_string)
	src_ceny = parse(Float64,src_ceny_string)
	return src_cenx,src_ceny
end
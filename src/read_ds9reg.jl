# Feb 13, 2018 : Written by G. C. Dewangan
#

#= This code reads ds9 region files (circular, annular, and elliptical) files and outputs
  region parameters circle : (xcen, ycen, radius), annulus: (xcen, ycen, inner_radius, outer_radius),
  ellipse : (xcen, ycen, major_radius, minor_radius, theta).

  The code is appropriate for single region, if there are multiple regions in the region file, the first region will be read.
=#


# export read_ds9circular_reg, read_ds9cirannu_reg, read_ds9elliptical_reg
"""
	read_ds9reg(ds9regfile)

Read a DS9 region file.

This function reads a region file created with the DS9 tool, and outputs the region parameters.
...
# Arguments
## Required parameters
- `ds9regfile::String` : DS9 region file for a circular, annular or elliptical region.
...
"""
function read_ds9reg(ds9regfile::String)
	f=open(ds9regfile)
	region=readlines(f)[4]
    reg_type = split(region,"(")[1]
	if reg_type =="circle"
		xcen = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[1])
		ycen  = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[2])
		radius = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[3])
		return reg_type, xcen, ycen, radius
	elseif reg_type =="annulus"
		xcen = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[1])
		ycen  = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[2])
		inner_radius = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[3])
		outer_radius = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[4])
		return reg_type, xcen, ycen, inner_radius, outer_radius
	elseif reg_type =="ellipse"
		xcen = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[1])
		ycen  = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[2])
		major_radius = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[3])
		minor_radius = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[4])
		theta = parse(Float64,split(split(split(region,"(")[2], ")")[1],",")[5])
		return reg_type, xcen, ycen, major_radius, minor_radius, theta
	else
        print("Region is not one Circle, Annulus or Ellipse")
        return 1
	end
end

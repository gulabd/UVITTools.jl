# test/runtests.jl

using UVITTools
using Test

@testset "Wavelength calculations" begin
	@test fuv_grating2_pixel2lamA(-500,order=-2) == 1447.2985865170435
	@test fuv_grating1_pixel2lamA(-500,order=-2) == 1440.3065041673246
end
# test/runtests.jl

using UVITTools
using Test


@testset "Effective Area calculation" begin
    @test fuv_grating2_ea(1400.0,order=-2) == 2.7117826719319313
    @test fuv_grating1_ea(1500.0,order=-2) == 4.23520016495318
end

@testset "Wavelength calculations" begin
	@test fuv_grating2_pixel2lamA(-500,order=-2) == 1447.2985865170435
	@test fuv_grating1_pixel2lamA(-500,order=-2) == 1440.3065041673246
end
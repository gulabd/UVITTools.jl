
"""
    lamA2lohi(lamA)
Generate wavelength bins (A_lo, A_hi) and energy bins (E_lo, E_hi) from an array of wavelengths in Å.

...
# Arguments
- `lamA::Array{Float64}`: An array of wavelengths in Å.

...
"""
function lamA2lohi(lamA::Array{Float64,1})
    sort!(lamA)
    nbins = length(lamA)
    lamA_lo = Array{Float64}(undef,nbins)
    lamA_hi = Array{Float64}(undef,nbins)


    for i=1:nbins
        if i==1
            lamA_lo[i] = lamA[i] - (lamA[i+1] - lamA[i])/2.0
            lamA_hi[i] = lamA[i] + (lamA[i+1] - lamA[i])/2.0
        elseif i < nbins
            lamA_hi[i] = lamA[i] + (lamA[i+1] - lamA[i])/2.0
            lamA_lo[i] = lamA_hi[i-1]
        else
            lamA_lo[i] = lamA_hi[i-1]
            lamA_hi[i] = lamA[i] + (lamA[i] - lamA[i-1])/2.0
    end

        println(lamA_lo[i])
       println(lamA_hi[i])
    end
    enkev_lo = reverse(lambdaA2keV.(lamA_hi))
    enkev_hi = reverse(lambdaA2keV.(lamA_lo))
    return lamA_lo, lamA_hi, enkev_lo, enkev_hi
end



module SpinHilbertSpaceMod

using ..HilbertSpaceMod
using ..BasisBuildersMod: build_particle_basis

mutable struct SpinHilbertSpace<:AbstractHilbertSpace
    nalpha::Int64
    nbeta::Int64
    norbital::Int64
    hilbert::Array{Int64,1}
end

function BuildHilbert(nparticle,hilbertspace::SpinHilbertSpace)
    norbital = hilbertspace.norbital
    hilbertspace.hilbert = build_particle_basis(Int64, nparticle, norbital; bitstep=2)
    return hilbertspace.hilbert
    
end

function BuildSpinHilbert(hilbertspace::SpinHilbertSpace)
    nalpha = hilbertspace.nalpha
    nbeta = hilbertspace.nbeta
    halpha = BuildHilbert(nalpha,hilbertspace)
    hbeta = BuildHilbert(nbeta,hilbertspace)

    for i in eachindex(hbeta)
        hbeta[i] = hbeta[i] << 1
    end
    
    h_tot = zeros(Int64,length(halpha)*length(hbeta))

    for i in eachindex(halpha)
        for j in eachindex(hbeta)
            h_tot[(i-1)*length(hbeta)+j] = halpha[i] | hbeta[j]
        end
    end

    return h_tot
end

end

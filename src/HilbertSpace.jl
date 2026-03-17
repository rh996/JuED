module HilbertSpaceMod

export AbstractHilbertSpace, ToDict

using ..BasisBuildersMod: build_particle_basis
using ..IndexTypesMod: state_eltype

abstract type AbstractHilbertSpace end

mutable struct GeneralHilbertSpace{Ti}<:AbstractHilbertSpace
    nparticle::Int64
    norbital::Int64
    hilbert::Array{Ti,1}
end

function BuildHilbert(hilbertspace::GeneralHilbertSpace{Ti}) where {Ti}
    nparticle = hilbertspace.nparticle
    norbital = hilbertspace.norbital
    hilbertspace.hilbert = build_particle_basis(Ti, nparticle, norbital)
    return hilbertspace.hilbert
    
end

function ToDict(hilbertspace::AbstractHilbertSpace, ::Type{Tv}=Int) where {Tv<:Integer}
    dict = Dict{state_eltype(hilbertspace),Tv}()
    for i in eachindex(hilbertspace.hilbert)
        dict[hilbertspace.hilbert[i]] = Tv(i)
    end
    return dict
    
end


end

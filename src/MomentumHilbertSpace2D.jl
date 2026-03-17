


module MomentumHilbertSpace2DMod
export MomentumHilbertSpace2D, BuildHilbert, momentum_add_2d, momentum_sub_2d

using ..HilbertSpaceMod
using ..BasisBuildersMod: build_momentum_basis
using ..MomentumUtilsMod: momentum_add_2d, momentum_sub_2d

mutable struct MomentumHilbertSpace2D{Ti}<:AbstractHilbertSpace
    nparticle::Int64
    Nkx::Int64
    Nky::Int64
    momentum::Int64
    hilbert::Array{Ti,1}
end

function BuildHilbert(hilbertspace::MomentumHilbertSpace2D{Ti}; use_cache::Bool=true) where {Ti}
    nparticle = hilbertspace.nparticle
    norbital = hilbertspace.Nkx*hilbertspace.Nky
    k = hilbertspace.momentum
    add_momentum = (k1, k2, systemsize) -> momentum_add_2d(k1, k2, hilbertspace.Nkx, hilbertspace.Nky)
    sub_momentum = (k1, k2, systemsize) -> momentum_sub_2d(k1, k2, hilbertspace.Nkx, hilbertspace.Nky)
    hilbertspace.hilbert = build_momentum_basis(Ti, nparticle, norbital, k, norbital, add_momentum, sub_momentum; use_cache)
    return hilbertspace.hilbert
    
end
end

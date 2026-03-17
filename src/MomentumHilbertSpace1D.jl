


module MomentumHilbertSpace1DMod

export MomentumHilbertSpace1D, BuildHilbert, momentum_add_1d, momentum_sub_1d

using ..HilbertSpaceMod
using ..BasisBuildersMod: build_momentum_basis
using ..MomentumUtilsMod: momentum_add_1d, momentum_sub_1d

mutable struct MomentumHilbertSpace1D{Ti}<:AbstractHilbertSpace
    nparticle::Int64
    norbital::Int64
    momentum::Int64
    hilbert::Array{Ti,1}
end

function BuildHilbert(hilbertspace::MomentumHilbertSpace1D{Ti}; use_cache::Bool=true) where {Ti}
    nparticle = hilbertspace.nparticle
    norbital = hilbertspace.norbital
    k = hilbertspace.momentum
    hilbertspace.hilbert = build_momentum_basis(Ti, nparticle, norbital, k, norbital, momentum_add_1d, momentum_sub_1d; use_cache)
    return hilbertspace.hilbert
    
end


end

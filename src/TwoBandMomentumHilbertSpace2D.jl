module TwoBandMomentumHilbertSpace2DMod

export TwoBandMomentumHilbertSpace2D, BuildHilbert, BuildTwoBandHilbert

using ..HilbertSpaceMod
using ..BasisBuildersMod: build_momentum_basis
using ..IndexTypesMod: state_eltype
using ..MomentumUtilsMod: momentum_add_2d, momentum_sub_2d

mutable struct TwoBandMomentumHilbertSpace2D{Ti}<:AbstractHilbertSpace
    nparticle::Int64
    Nkx::Int64
    Nky::Int64
    momentum::Int64
    hilbert::Array{Ti,1}
end


function BuildHilbert(nparticle,momentum,hilbertspace::TwoBandMomentumHilbertSpace2D; use_cache::Bool=true)
    norbital = hilbertspace.Nkx*hilbertspace.Nky
    k = momentum
    Ti = state_eltype(hilbertspace)
    add_momentum = (k1, k2, systemsize) -> momentum_add_2d(k1, k2, hilbertspace.Nkx, hilbertspace.Nky)
    sub_momentum = (k1, k2, systemsize) -> momentum_sub_2d(k1, k2, hilbertspace.Nkx, hilbertspace.Nky)
    hilbert = build_momentum_basis(Ti, nparticle, norbital, k, norbital, add_momentum, sub_momentum; bitstep=2, use_cache)
    return hilbert
    
end

function BuildTwoBandHilbert(hilbertspace::TwoBandMomentumHilbertSpace2D; use_cache::Bool=true)
    nparticle = hilbertspace.nparticle
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    norbital = hilbertspace.Nkx*hilbertspace.Nky
    momentum = hilbertspace.momentum
    Ti = state_eltype(hilbertspace)
    result = Vector{Ti}()

    for nalpha in 0:nparticle
        nbeta = nparticle - nalpha
        for k1 in 0:norbital-1
            for k2 in 0:norbital-1
                k_tot = momentum_add_2d(k1,k2,Nkx,Nky)
                if k_tot == momentum
                    hilbert_alpha = BuildHilbert(nalpha,k1,hilbertspace; use_cache)
                    hilbert_beta = BuildHilbert(nbeta,k2,hilbertspace; use_cache)
                    for i in eachindex(hilbert_beta)
                        hilbert_beta[i] = hilbert_beta[i] << 1
                    end
                    hilbert_tot = zeros(Ti,length(hilbert_alpha)*length(hilbert_beta))
                    for i in eachindex(hilbert_alpha)
                        for j in eachindex(hilbert_beta)
                            hilbert_tot[(i-1)*length(hilbert_beta)+j] = hilbert_alpha[i] | hilbert_beta[j]
                        end
                    end
                    result = vcat(result,hilbert_tot)
                end
            end
        end
    end




    hilbertspace.hilbert = result
    return result
end





end

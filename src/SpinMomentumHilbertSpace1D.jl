
module SpinMomentumHilbertSpace1DMod

using ..HilbertSpaceMod
using ..BasisBuildersMod: build_momentum_basis
using ..MomentumUtilsMod: momentum_add_1d, momentum_sub_1d

mutable struct SpinMomentumHilbertSpace1D{Ti}<:AbstractHilbertSpace
    nalpha::Int64
    nbeta::Int64
    norbital::Int64
    momentum::Int64
    hilbert::Array{Ti,1}
end

function BuildHilbert(nparticle,momentum,hilbertspace::SpinMomentumHilbertSpace1D)
    
    norbital = hilbertspace.norbital
    k = momentum
    Ti = typeof(hilbertspace.hilbert).parameters[1]
    hilbert = build_momentum_basis(Ti, nparticle, norbital, k, norbital, momentum_add_1d, momentum_sub_1d; bitstep=2)
    return hilbert
    
end

function BuildSpinHilbert(hilbertspace::SpinMomentumHilbertSpace1D)
    nalpha = hilbertspace.nalpha
    nbeta = hilbertspace.nbeta
    norbital = hilbertspace.norbital
    momentum = hilbertspace.momentum

    Ti = typeof(hilbertspace.hilbert).parameters[1]
    result = Vector{Ti}()
    for k1 in 0:norbital-1
        for k2 in 0:norbital-1
            k_tot = momentum_add_1d(k1,k2,norbital)
            if k_tot == momentum
                hilbert_alpha = BuildHilbert(nalpha,k1,hilbertspace)
                hilbert_beta = BuildHilbert(nbeta,k2,hilbertspace)
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
    hilbertspace.hilbert = result
    return result


    
    
end

export SpinMomentumHilbertSpace1D, momentum_add_1d, momentum_sub_1d

end

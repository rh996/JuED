
module SpinMomentumHilbertSpace1DMod

using ..HilbertSpaceMod

mutable struct SpinMomentumHilbertSpace1D{Ti}<:AbstractHilbertSpace
    nalpha::Int64
    nbeta::Int64
    norbital::Int64
    momentum::Int64
    hilbert::Array{Ti,1}
end


function momentum_add_1d(k1,k2,systemsize)
    return mod((k1+k2),systemsize)
    
end

function momentum_sub_1d(k1,k2,systemsize)
    return mod((k1-k2),systemsize)
    
end

function _dfs(ne::Int64,no::Int64,k_curr::Int64,hilbertspace::SpinMomentumHilbertSpace1D)
    systemsize = hilbertspace.norbital
    if ne==no
        k_new = 0
        for i in 0:no-1
            k1 = (systemsize - 1-i)
            k_new = momentum_add_1d(k_new,k1,systemsize)
        end
        if k_new == k_curr
            a::Int64 = 0
            for i in 1:no
                a |= 1<<((i-1)*2)
            end
            return [a]
        else
            return []
        end
    elseif ne==0
        if k_curr == 0
            return [0]
        else
            return []
        end
    elseif ne>no
        return []
    else
        kind = systemsize-no
        k1 = k_curr
        k2 = kind
        k_new = momentum_sub_1d(k1,k2,systemsize)
        
        left = _dfs(ne,no-1,k_curr,hilbertspace)
        right = _dfs(ne-1,no-1,k_new,hilbertspace)
        shifted_right = right .+ (1<<((no-1)*2))
        curr = vcat(left,shifted_right)
        return curr
    end
end

function BuildHilbert(nparticle,momentum,hilbertspace::SpinMomentumHilbertSpace1D)
    
    norbital = hilbertspace.norbital
    k = momentum
    hilbert = _dfs(nparticle,norbital,k,hilbertspace)
    return hilbert
    
end

function BuildSpinHilbert(hilbertspace::SpinMomentumHilbertSpace1D)
    nalpha = hilbertspace.nalpha
    nbeta = hilbertspace.nbeta
    norbital = hilbertspace.norbital
    momentum = hilbertspace.momentum

    result :: typeof(hilbertspace.hilbert)= []
    for k1 in 0:norbital-1
        for k2 in 0:norbital-1
            k_tot = momentum_add_1d(k1,k2,norbital)
            if k_tot == momentum
                hilbert_alpha = BuildHilbert(nalpha,k1,hilbertspace)
                hilbert_beta = BuildHilbert(nbeta,k2,hilbertspace)
                for i in eachindex(hilbert_beta)
                    hilbert_beta[i] = hilbert_beta[i] << 1
                end
                hilbert_tot = zeros(typeof(hilbertspace.hilbert).parameters[1],length(hilbert_alpha)*length(hilbert_beta))
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

export SpinMomentumHilbertSpace1D

end
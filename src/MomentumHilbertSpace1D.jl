


module MomentumHilbertSpace1DMod

using ..HilbertSpaceMod

mutable struct MomentumHilbertSpace1D{Ti}<:AbstractHilbertSpace
    nparticle::Int64
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

function _dfs(ne::Int64,no::Int64,k_curr::Int64,hilbertspace::MomentumHilbertSpace1D)
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
                a |= 1<<(i-1)
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
        shifted_right = right .+ (1<<(no-1))
        curr = vcat(left,shifted_right)
        return curr
    end
end

function BuildHilbert(hilbertspace::MomentumHilbertSpace1D)
    nparticle = hilbertspace.nparticle
    norbital = hilbertspace.norbital
    k = hilbertspace.momentum
    hilbertspace.hilbert = _dfs(nparticle,norbital,k,hilbertspace)
    return hilbertspace.hilbert
    
end


end
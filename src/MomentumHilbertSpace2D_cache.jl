

module MomentumHilbertSpace2DMod
export momentum_add_2d, momentum_sub_2d

using ..HilbertSpaceMod
mutable struct MomentumHilbertSpace2D{Ti}<:AbstractHilbertSpace
    nparticle::Int64
    Nkx::Int64
    Nky::Int64
    momentum::Int64
    hilbert::Array{Ti,1}
end


@inline function momentum_add_2d(k1,k2,Nkx,Nky)
    k1x = mod(k1,Nkx)
    k1y = fld(k1,Nkx)
    k2x = mod(k2,Nkx)
    k2y = fld(k2,Nkx)
    return mod((k1x+k2x),Nkx) + mod((k1y+k2y),Nky)*Nkx
    
end


@inline function momentum_sub_2d(k1,k2,Nkx,Nky)
    k1x = mod(k1,Nkx)
    k1y = fld(k1,Nkx)
    k2x = mod(k2,Nkx)
    k2y = fld(k2,Nkx)
    return mod((k1x-k2x),Nkx) + mod((k1y-k2y),Nky)*Nkx
    
end

function _dfs(ne::Int64,no::Int64,k_curr::Int64,hilbertspace::MomentumHilbertSpace2D,cache::Dict)
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    systemsize = Nkx*Nky

    if ne>no
        return []
    end
    if ne==no
        k_new = 0
        for i in 0:no-1
            k1 = (systemsize - 1-i)
            k_new = momentum_add_2d(k_new,k1,Nkx,Nky)
        end
        if k_new == k_curr
            a::Int64 = 0
            for i in 1:ne
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
    end

    key = (ne,no,k_curr)
    if haskey(cache,key)
        return cache[key]
        
    else
        kind = systemsize-no
        k1 = k_curr
        k2 = kind
        k_new = momentum_sub_2d(k1,k2,Nkx,Nky)
        
        left = _dfs(ne,no-1,k_curr,hilbertspace,cache)
        right = _dfs(ne-1,no-1,k_new,hilbertspace,cache)
        shifted_right = right .+ (1<<(no-1))
        curr = vcat(left,shifted_right)
        cache[key] = curr
        return curr
    end
end

function BuildHilbert(hilbertspace::MomentumHilbertSpace2D)
    cache = Dict()
    nparticle = hilbertspace.nparticle
    norbital = hilbertspace.Nkx*hilbertspace.Nky
    k = hilbertspace.momentum
    hilbertspace.hilbert = _dfs(nparticle,norbital,k,hilbertspace,cache)
    return hilbertspace.hilbert
    
end
end
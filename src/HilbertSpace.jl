module HilbertSpaceMod

export AbstractHilbertSpace, ToDict

abstract type AbstractHilbertSpace end

mutable struct GeneralHilbertSpace{Ti}<:AbstractHilbertSpace
    nparticle::Int64
    norbital::Int64
    hilbert::Array{Ti,1}
end


function _dfs(ne::Int64,no::Int64,cache::Dict,hilbertspace::GeneralHilbertSpace)
    if ne == 0
        return [0]
    end
    if ne==no
        a::Int64 = 0
        for i in 1:no
            a |= 1<<(i-1)
        end
        return [a]
    end

    key = (ne,no)
    if haskey(cache,key)
        return cache[key]
    else

        left = _dfs(ne,no-1,cache,hilbertspace)
        right = _dfs(ne-1,no-1,cache,hilbertspace)
        shifted_right = right .+ (1<<(no-1))
        curr = vcat(left,shifted_right)

        cache[key] = curr
        return curr
    end
end

function BuildHilbert(hilbertspace::GeneralHilbertSpace)
    cache = Dict()
    nparticle = hilbertspace.nparticle
    norbital = hilbertspace.norbital
    hilbertspace.hilbert = _dfs(nparticle,norbital,cache,hilbertspace)
    return hilbertspace.hilbert
    
end

function ToDict(hilbertspace::AbstractHilbertSpace)
    dict = Dict()
    for i in eachindex(hilbertspace.hilbert)
        dict[hilbertspace.hilbert[i]] = i
    end
    return dict
    
end


end


module SpinHilbertSpaceMod

using ..HilbertSpaceMod

mutable struct SpinHilbertSpace<:AbstractHilbertSpace
    nalpha::Int64
    nbeta::Int64
    norbital::Int64
    hilbert::Array{Int64,1}
end



function _dfs(ne::Int64,no::Int64,cache,hilbertspace::SpinHilbertSpace)
    if ne == 0
        return [0]
    end
    if ne==no
        a::Int64 = 0
        for i in 0:no-1
            a |= 1<< (i*2)
        end
        return [a]
    end

    key = (ne,no)
    if haskey(cache,key)
        return cache[key]
    else

        left = _dfs(ne,no-1,cache,hilbertspace)
        right = _dfs(ne-1,no-1,cache,hilbertspace)
        shifted_right = right .+ (1<<((no-1)*2))
        curr = vcat(left,shifted_right)

        cache[key] = curr
        return curr
    end
end

function BuildHilbert(nparticle,hilbertspace::SpinHilbertSpace)
    cache = Dict()
    
    norbital = hilbertspace.norbital
    hilbertspace.hilbert = _dfs(nparticle,norbital,cache,hilbertspace)
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
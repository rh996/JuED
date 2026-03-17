module BasisBuildersMod

export build_particle_basis, build_momentum_basis

function _filled_state(::Type{Ti}, no::Int, bitstep::Int) where {Ti<:Integer}
    state = zero(Ti)
    for i in 0:(no - 1)
        state |= one(Ti) << (i * bitstep)
    end
    return state
end

function _particle_dfs(::Type{Ti}, ne::Int, no::Int, bitstep::Int, cache::Union{Nothing,Dict{Tuple{Int,Int},Vector{Ti}}}) where {Ti<:Integer}
    if ne == 0
        return Ti[0]
    end
    if ne == no
        return Ti[_filled_state(Ti, no, bitstep)]
    end

    key = (ne, no)
    if cache !== nothing && haskey(cache, key)
        return cache[key]
    end

    left = _particle_dfs(Ti, ne, no - 1, bitstep, cache)
    right = _particle_dfs(Ti, ne - 1, no - 1, bitstep, cache)
    shifted_right = right .+ (one(Ti) << ((no - 1) * bitstep))
    curr = vcat(left, shifted_right)
    if cache !== nothing
        cache[key] = curr
    end
    return curr
end

function build_particle_basis(::Type{Ti}, ne::Int, no::Int; bitstep::Int=1, use_cache::Bool=true) where {Ti<:Integer}
    cache = use_cache ? Dict{Tuple{Int,Int},Vector{Ti}}() : nothing
    return _particle_dfs(Ti, ne, no, bitstep, cache)
end

function _momentum_dfs(
    ::Type{Ti},
    ne::Int,
    no::Int,
    k_curr::Int,
    systemsize::Int,
    add_momentum,
    sub_momentum,
    bitstep::Int,
    cache::Union{Nothing,Dict{NTuple{3,Int},Vector{Ti}}},
) where {Ti<:Integer}
    if ne > no
        return Vector{Ti}()
    end
    if ne == no
        k_new = 0
        for i in 0:(no - 1)
            k_new = add_momentum(k_new, systemsize - 1 - i, systemsize)
        end
        if k_new == k_curr
            return Ti[_filled_state(Ti, ne, bitstep)]
        end
        return Vector{Ti}()
    end
    if ne == 0
        if k_curr == 0
            return Ti[0]
        end
        return Vector{Ti}()
    end

    key = (ne, no, k_curr)
    if cache !== nothing && haskey(cache, key)
        return cache[key]
    end

    k_new = sub_momentum(k_curr, systemsize - no, systemsize)
    left = _momentum_dfs(Ti, ne, no - 1, k_curr, systemsize, add_momentum, sub_momentum, bitstep, cache)
    right = _momentum_dfs(Ti, ne - 1, no - 1, k_new, systemsize, add_momentum, sub_momentum, bitstep, cache)
    shifted_right = right .+ (one(Ti) << ((no - 1) * bitstep))
    curr = vcat(left, shifted_right)
    if cache !== nothing
        cache[key] = curr
    end
    return curr
end

function build_momentum_basis(::Type{Ti}, ne::Int, no::Int, k::Int, systemsize::Int, add_momentum, sub_momentum; bitstep::Int=1, use_cache::Bool=true) where {Ti<:Integer}
    cache = use_cache ? Dict{NTuple{3,Int},Vector{Ti}}() : nothing
    return _momentum_dfs(Ti, ne, no, k, systemsize, add_momentum, sub_momentum, bitstep, cache)
end

end

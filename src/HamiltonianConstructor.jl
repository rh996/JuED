module HamiltonianConstructorMod

export Hamiltonian_list_constructor, Hamiltonian_momentum_constructor, HamiltonianAction

using ..MomentumHilbertSpace2DMod
using ..FermionOperatorMod
using ..IndexTypesMod: choose_pointer_type
using SparseArrays

AbstractHilbertSpace = MomentumHilbertSpace2DMod.AbstractHilbertSpace
ToDict = MomentumHilbertSpace2DMod.ToDict

struct OneBodyChannel{T}
    annihilate::Int
    create::Int
    value::T
end

struct TwoBodyChannel{T}
    annihilate1::Int
    annihilate2::Int
    create1::Int
    create2::Int
    value::T
end

function GeneralHamiltonian(OneBody::Array{Float64,2}, TwoBody::Array{Float64,4}, hilbertspace::AbstractHilbertSpace)
    hilbert = hilbertspace.hilbert
    ind_dict = ToDict(hilbertspace)
    ndim = size(OneBody)[1]
    data::Array{Float64,1} = []
    column::Array{Int64,1} = []
    row::Array{Int64,1} = []
    for i in 1:ndim
        for j in 1:ndim
            for state_ind in eachindex(hilbert)
                fermion = FermionOperator(hilbert[state_ind], 1)
                AnnihilationOperator!(fermion, j)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, i)
                if fermion.fermion_sign == 0
                    continue
                end
                push!(data, OneBody[i, j] * fermion.fermion_sign)
                push!(column, state_ind)
                push!(row, ind_dict[fermion.state])
            end
        end
    end

    for i in 1:ndim
        for j in 1:ndim
            for k in 1:ndim
                for l in 1:ndim
                    for state_ind in eachindex(hilbert)
                        fermion1 = FermionOperator(hilbert[state_ind], 1)
                        AnnihilationOperator!(fermion1, k)
                        if fermion1.fermion_sign == 0
                            continue
                        end
                        AnnihilationOperator!(fermion1, l)
                        if fermion1.fermion_sign == 0
                            continue
                        end
                        CreationOperator!(fermion1, j)
                        if fermion1.fermion_sign == 0
                            continue
                        end
                        CreationOperator!(fermion1, i)
                        if fermion1.fermion_sign == 0
                            continue
                        end
                        push!(data, TwoBody[i, j, k, l] * fermion1.fermion_sign)
                        push!(column, state_ind)
                        push!(row, ind_dict[fermion1.state])
                    end
                end
            end
        end
    end

    return sparse(row, column, data, length(hilbert), length(hilbert))
end

@inline function _chunk_bounds(t::Int, chunk_size::Int, total_size::Int)
    start_idx = (t - 1) * chunk_size + 1
    end_idx = min(t * chunk_size, total_size)
    return start_idx, end_idx
end

function _onebody_channels(OneBody)
    norbital = size(OneBody, 1)
    channels = Vector{OneBodyChannel{eltype(OneBody)}}()
    for ind in findall(x -> abs(x) > 1e-10, OneBody)
        i = ind[1]
        j = ind[2]
        push!(channels, OneBodyChannel((norbital - j) + 1, (norbital - i) + 1, OneBody[i, j]))
    end
    return channels
end

function _list_twobody_channels(Eri)
    norbital = size(Eri, 1)
    channels = Vector{TwoBodyChannel{eltype(Eri)}}()
    for ind in findall(x -> abs(x) > 1e-10, Eri)
        i = ind[1]
        j = ind[2]
        k = ind[3]
        l = ind[4]
        push!(channels, TwoBodyChannel((norbital - k) + 1, (norbital - l) + 1, (norbital - j) + 1, (norbital - i) + 1, Eri[i, j, k, l]))
    end
    return channels
end

function _momentum_twobody_channels(Eri, hilbertspace)
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    nk = Nkx * Nky
    channels = Vector{TwoBodyChannel{eltype(Eri)}}()

    for ik in 0:(nk - 1)
        ikx = mod(ik, Nkx)
        iky = fld(ik, Nkx)
        for ikp in 0:(nk - 1)
            ikpx = mod(ikp, Nkx)
            ikpy = fld(ikp, Nkx)
            for iq in 0:(nk - 1)
                iqx = mod(iq, Nkx)
                iqy = fld(iq, Nkx)
                value = Eri[ikx+1, iky+1, ikpx+1, ikpy+1, iqx+1, iqy+1]
                if abs(value) <= 1e-10
                    continue
                end

                create2 = nk - momentum_add_2d(ik, iq, Nkx, Nky)
                annihilate1 = nk - ik
                create1 = nk - momentum_sub_2d(ikp, iq, Nkx, Nky)
                annihilate2 = nk - ikp
                push!(channels, TwoBodyChannel(annihilate1, annihilate2, create1, create2, value))
            end
        end
    end

    return channels
end

function _hamiltonian_channels(OneBody, Eri::Array{T,4}, hilbertspace) where {T}
    return _onebody_channels(OneBody), _list_twobody_channels(Eri)
end

function _hamiltonian_channels(OneBody, Eri::Array{T,6}, hilbertspace) where {T}
    return _onebody_channels(OneBody), _momentum_twobody_channels(Eri, hilbertspace)
end

@inline function _apply_onebody_channel(state, channel::OneBodyChannel)
    fermion = FermionOperator(state, 1)
    AnnihilationOperator!(fermion, channel.annihilate)
    if fermion.fermion_sign == 0
        return 0, state
    end
    CreationOperator!(fermion, channel.create)
    if fermion.fermion_sign == 0
        return 0, state
    end
    return fermion.fermion_sign, fermion.state
end

@inline function _apply_twobody_channel(state, channel::TwoBodyChannel)
    fermion = FermionOperator(state, 1)
    AnnihilationOperator!(fermion, channel.annihilate1)
    if fermion.fermion_sign == 0
        return 0, state
    end
    AnnihilationOperator!(fermion, channel.annihilate2)
    if fermion.fermion_sign == 0
        return 0, state
    end
    CreationOperator!(fermion, channel.create1)
    if fermion.fermion_sign == 0
        return 0, state
    end
    CreationOperator!(fermion, channel.create2)
    if fermion.fermion_sign == 0
        return 0, state
    end
    return fermion.fermion_sign, fermion.state
end

@inline function _count_connected_terms(state, onebody_channels, twobody_channels)
    local_count = 0

    for channel in onebody_channels
        sign, _ = _apply_onebody_channel(state, channel)
        if sign != 0
            local_count += 1
        end
    end

    for channel in twobody_channels
        sign, _ = _apply_twobody_channel(state, channel)
        if sign != 0
            local_count += 1
        end
    end

    return local_count
end

@inline function _fill_connected_terms!(state, writepos, data, row, ind_dict, onebody_channels, twobody_channels)
    for channel in onebody_channels
        sign, next_state = _apply_onebody_channel(state, channel)
        if sign == 0
            continue
        end
        data[writepos] = channel.value * sign
        row[writepos] = ind_dict[next_state]
        writepos += 1
    end

    for channel in twobody_channels
        sign, next_state = _apply_twobody_channel(state, channel)
        if sign == 0
            continue
        end
        data[writepos] = channel.value * sign
        row[writepos] = ind_dict[next_state]
        writepos += 1
    end

    return nothing
end

function _assemble_hamiltonian_csc(::Type{T}, hilbertspace, onebody_channels, twobody_channels) where {T}
    hilbert = hilbertspace.hilbert
    dim_hilbert = length(hilbert)
    nthreads_total = Threads.nthreads()
    chunk_size = ceil(Int, dim_hilbert / nthreads_total)
    nnz_col = zeros(Int, dim_hilbert)

    Threads.@threads for t in 1:nthreads_total
        start_idx, end_idx = _chunk_bounds(t, chunk_size, dim_hilbert)
        @inbounds for state_ind in start_idx:end_idx
            nnz_col[state_ind] = _count_connected_terms(hilbert[state_ind], onebody_channels, twobody_channels)
        end
    end

    pointertype = choose_pointer_type(sum(nnz_col) + 1)
    indptr = Array{pointertype}(undef, dim_hilbert + 1)
    indptr[1] = 1
    for col in 1:dim_hilbert
        indptr[col+1] = indptr[col] + nnz_col[col]
    end

    total_nnz = indptr[end] - 1
    data = Array{T}(undef, total_nnz)
    row = Array{pointertype}(undef, total_nnz)
    ind_dict = ToDict(hilbertspace, pointertype)

    Threads.@threads for t in 1:nthreads_total
        start_idx, end_idx = _chunk_bounds(t, chunk_size, dim_hilbert)
        @inbounds for state_ind in start_idx:end_idx
            _fill_connected_terms!(hilbert[state_ind], indptr[state_ind], data, row, ind_dict, onebody_channels, twobody_channels)
        end
    end

    return data, row, indptr, dim_hilbert
end

function _hamiltonian_action(coefftype::Type, hilbertspace, onebody_channels, twobody_channels)
    hilbert = hilbertspace.hilbert
    dim = length(hilbert)
    ind_dict = ToDict(hilbertspace)

    function action(v::AbstractVector)
        result_type = promote_type(eltype(v), coefftype)
        local_results = [zeros(result_type, dim) for _ in 1:Threads.nthreads()]
        nthreads_total = Threads.nthreads()
        chunk_size = ceil(Int, dim / nthreads_total)

        Threads.@threads for t in 1:nthreads_total
            start_idx, end_idx = _chunk_bounds(t, chunk_size, dim)
            local_result = local_results[Threads.threadid()]
            @inbounds for state_ind in start_idx:end_idx
                amplitude = v[state_ind]
                if iszero(amplitude)
                    continue
                end

                state = hilbert[state_ind]
                for channel in onebody_channels
                    sign, next_state = _apply_onebody_channel(state, channel)
                    if sign == 0
                        continue
                    end
                    local_result[ind_dict[next_state]] += channel.value * sign * amplitude
                end

                for channel in twobody_channels
                    sign, next_state = _apply_twobody_channel(state, channel)
                    if sign == 0
                        continue
                    end
                    local_result[ind_dict[next_state]] += channel.value * sign * amplitude
                end
            end
        end

        result = zeros(result_type, dim)
        for local_result in local_results
            result .+= local_result
        end
        return result
    end

    return action, dim
end

function Hamiltonian_momentum_constructor(OneBody, Eri::Array{T,6}, hilbertspace) where {T}
    onebody_channels, twobody_channels = _hamiltonian_channels(OneBody, Eri, hilbertspace)
    return _assemble_hamiltonian_csc(T, hilbertspace, onebody_channels, twobody_channels)
end

function Hamiltonian_list_constructor(OneBody, Eri::Array{T,4}, hilbertspace) where {T}
    onebody_channels, twobody_channels = _hamiltonian_channels(OneBody, Eri, hilbertspace)
    return _assemble_hamiltonian_csc(T, hilbertspace, onebody_channels, twobody_channels)
end

function HamiltonianAction(OneBody, Eri::Array{T,4}, hilbertspace) where {T}
    onebody_channels, twobody_channels = _hamiltonian_channels(OneBody, Eri, hilbertspace)
    coefftype = promote_type(eltype(OneBody), T)
    return _hamiltonian_action(coefftype, hilbertspace, onebody_channels, twobody_channels)
end

function HamiltonianAction(OneBody, Eri::Array{T,6}, hilbertspace) where {T}
    onebody_channels, twobody_channels = _hamiltonian_channels(OneBody, Eri, hilbertspace)
    coefftype = promote_type(eltype(OneBody), T)
    return _hamiltonian_action(coefftype, hilbertspace, onebody_channels, twobody_channels)
end

end

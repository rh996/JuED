module DensityMatricesMod

using ..FermionOperatorMod
using ..HilbertSpaceMod: ToDict
using ..IndexTypesMod: choose_pointer_type, choose_state_type
using ..ModelTypesMod
using ..MomentumHilbertSpace2DMod
using ..MomentumUtilsMod: momentum_add_2d
using ..TwoBandMomentumHilbertSpace2DMod
using JLD2

export RDMWorkspace, CompactRDM2, CompactRDM3, RDM1, RDM2, RDM3, RDM2Compact, RDM3Compact, RDM3_single_thread, RDM2_cache, RDM3_naive, RDM2_naive, todense

const PairTuple = NTuple{2,Int}
const TripleTuple = NTuple{3,Int}
const QuadTuple = NTuple{4,Int}
const SextTuple = NTuple{6,Int}
const RDM3_LEFT_PERMS = ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
const RDM3_RIGHT_PERMS = ((4, 5, 6), (4, 6, 5), (5, 4, 6), (5, 6, 4), (6, 4, 5), (6, 5, 4))
const RDM2_CACHE_SCHEMA_VERSION = 1

struct RDMSectorWorkspace{TH,TS<:Integer,TI<:Integer}
    kind::Symbol
    nparticle::Int
    momentum::Int
    Nkx::Int
    Nky::Int
    norb::Int
    hilbertspace::TH
    hilbert::Vector{TS}
    ind_dict::Dict{TS,TI}
    rdm2_elements::Vector{QuadTuple}
    rdm3_elements::Vector{SextTuple}
end

struct CompactRDM2{T}
    norb::Int
    elements::Vector{QuadTuple}
    values::Vector{T}
end

struct CompactRDM3{T}
    norb::Int
    elements::Vector{SextTuple}
    values::Vector{T}
end

function orbital_pairs(norb::Int)
    pairs = Vector{PairTuple}()
    sizehint!(pairs, binomial(norb, 2))
    for i in 1:norb
        for j in (i + 1):norb
            push!(pairs, (i, j))
        end
    end
    return pairs
end

function orbital_triples(norb::Int)
    triples = Vector{TripleTuple}()
    sizehint!(triples, binomial(norb, 3))
    for i in 1:norb
        for j in (i + 1):norb
            for k in (j + 1):norb
                push!(triples, (i, j, k))
            end
        end
    end
    return triples
end

@inline spinless_orbital_momentum(i::Int) = i - 1
@inline twoband_orbital_momentum(i::Int) = fld(i - 1, 2)

function momentum_sum_2d(indices::Tuple, orbital_momentum, Nkx::Int, Nky::Int)
    total = 0
    for index in indices
        total = momentum_add_2d(total, orbital_momentum(index), Nkx, Nky)
    end
    return total
end

function momentum_conserving_rdm2_elements(
    pairs::Vector{PairTuple},
    Nkx::Int,
    Nky::Int;
    orbital_momentum,
    reorder_second_pair::Function=identity,
)
    independent_elements = Vector{QuadTuple}()
    number_of_pairs = length(pairs)
    for ind1 in eachindex(pairs)
        i, j = pairs[ind1]
        left_momentum = momentum_sum_2d((i, j), orbital_momentum, Nkx, Nky)
        for ind2 in ind1:number_of_pairs
            k0, l0 = pairs[ind2]
            k, l = reorder_second_pair((k0, l0))
            right_momentum = momentum_sum_2d((k, l), orbital_momentum, Nkx, Nky)
            if left_momentum == right_momentum
                push!(independent_elements, (i, j, k, l))
            end
        end
    end
    return independent_elements
end

function momentum_conserving_rdm3_elements(
    triples::Vector{TripleTuple},
    Nkx::Int,
    Nky::Int;
    orbital_momentum,
)
    independent_elements = Vector{SextTuple}()
    number_of_triples = length(triples)
    for ind1 in eachindex(triples)
        left_triple = triples[ind1]
        left_momentum = momentum_sum_2d(left_triple, orbital_momentum, Nkx, Nky)
        for ind2 in ind1:number_of_triples
            right_triple = triples[ind2]
            if left_momentum == momentum_sum_2d(right_triple, orbital_momentum, Nkx, Nky)
                i, j, k = left_triple
                l, m, n = right_triple
                push!(independent_elements, (i, j, k, l, m, n))
            end
        end
    end
    return independent_elements
end

function RDMWorkspace(model::ModelParams2DSpinlessList, momentum::Int; use_cache::Bool=true)
    Nkx = model.Nkx
    Nky = model.Nky
    norb = Nkx * Nky
    state_type = choose_state_type(norb)
    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{state_type}(model.nparticle, Nkx, Nky, momentum, state_type[])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace; use_cache)
    pointertype = choose_pointer_type(length(hilbert))
    pairs = orbital_pairs(norb)
    triples = orbital_triples(norb)
    return RDMSectorWorkspace(
        :spinless_2d,
        model.nparticle,
        momentum,
        Nkx,
        Nky,
        norb,
        hilbertspace,
        hilbert,
        ToDict(hilbertspace, pointertype),
        momentum_conserving_rdm2_elements(
            pairs,
            Nkx,
            Nky;
            orbital_momentum=spinless_orbital_momentum,
            reorder_second_pair=pair -> (pair[2], pair[1]),
        ),
        momentum_conserving_rdm3_elements(triples, Nkx, Nky; orbital_momentum=spinless_orbital_momentum),
    )
end

function RDMWorkspace(model::ModelParams2DTwoBand, momentum::Int; use_cache::Bool=true)
    Nkx = model.Nkx
    Nky = model.Nky
    norb = Nkx * Nky * 2
    state_type = choose_state_type(norb)
    hilbertspace = TwoBandMomentumHilbertSpace2DMod.TwoBandMomentumHilbertSpace2D{state_type}(model.nparticle, Nkx, Nky, momentum, state_type[])
    hilbert = TwoBandMomentumHilbertSpace2DMod.BuildTwoBandHilbert(hilbertspace; use_cache)
    pointertype = choose_pointer_type(length(hilbert))
    pairs = orbital_pairs(norb)
    return RDMSectorWorkspace(
        :twoband_2d,
        model.nparticle,
        momentum,
        Nkx,
        Nky,
        norb,
        hilbertspace,
        hilbert,
        ToDict(hilbertspace, pointertype),
        momentum_conserving_rdm2_elements(pairs, Nkx, Nky; orbital_momentum=twoband_orbital_momentum),
        SextTuple[],
    )
end

function _check_coefficients(workspace::RDMSectorWorkspace, coeffs::AbstractVector)
    if length(coeffs) != length(workspace.hilbert)
        throw(ArgumentError("Coefficient vector length $(length(coeffs)) does not match Hilbert-space dimension $(length(workspace.hilbert))."))
    end
end

function _lookup_transition(workspace::RDMSectorWorkspace{<:Any,TS}, state::TS, creation_sites::Tuple, annihilation_sites::Tuple) where {TS<:Integer}
    internal_creation_sites = map(site -> basis_site_index(workspace.norb, site), creation_sites)
    internal_annihilation_sites = map(site -> basis_site_index(workspace.norb, site), annihilation_sites)
    result = apply_operator_string(state, internal_creation_sites, internal_annihilation_sites)
    result === nothing && return nothing
    left_state, sign = result
    if !haskey(workspace.ind_dict, left_state)
        return nothing
    end
    return workspace.ind_dict[left_state], sign
end

@inline _lookup_rdm1_transition(workspace::RDMSectorWorkspace, state, i::Int, j::Int) = _lookup_transition(workspace, state, (i,), (j,))
@inline _lookup_rdm2_transition(workspace::RDMSectorWorkspace, state, i::Int, j::Int, k::Int, l::Int) = _lookup_transition(workspace, state, (j, i), (k, l))
@inline _lookup_rdm3_transition(workspace::RDMSectorWorkspace, state, i::Int, j::Int, k::Int, l::Int, m::Int, n::Int) = _lookup_transition(workspace, state, (k, j, i), (l, m, n))

@inline function _transition_value(coeffs::AbstractVector{T}, right, left, sign) where {T}
    return sign * coeffs[right] * conj(coeffs[left])
end

@inline function _chunk_bounds(dim::Int, nthreads_total::Int, t::Int)
    chunk_size = cld(dim, nthreads_total)
    start_idx = (t - 1) * chunk_size + 1
    end_idx = min(t * chunk_size, dim)
    return start_idx, end_idx
end

function _antisymmetrize_rdm2!(rdm2)
    rdm2 .-= permutedims(rdm2, (2, 1, 3, 4))
    rdm2 .-= permutedims(rdm2, (1, 2, 4, 3))
    return rdm2
end

function permutation_sign(perm)
    sign = 1
    n = length(perm)
    for i in 1:n
        for j in (i + 1):n
            if perm[i] > perm[j]
                sign *= -1
            end
        end
    end
    return sign
end

function _antisymmetrize_rdm3(rdm3)
    result = zeros(eltype(rdm3), size(rdm3)...)
    for eta in RDM3_LEFT_PERMS
        for sigma in RDM3_RIGHT_PERMS
            sign = permutation_sign(eta) * permutation_sign(sigma)
            result .+= sign .* permutedims(rdm3, (eta..., sigma...))
        end
    end
    return result
end

function todense(rdm::CompactRDM2{T}) where {T}
    dense = zeros(T, rdm.norb, rdm.norb, rdm.norb, rdm.norb)
    for idx in eachindex(rdm.elements)
        i, j, k, l = rdm.elements[idx]
        value = rdm.values[idx]
        dense[i, j, k, l] += value
        if (i, j) != (k, l)
            dense[k, l, i, j] += conj(value)
        end
    end
    return _antisymmetrize_rdm2!(dense)
end

function todense(rdm::CompactRDM3{T}) where {T}
    dense = zeros(T, rdm.norb, rdm.norb, rdm.norb, rdm.norb, rdm.norb, rdm.norb)
    for idx in eachindex(rdm.elements)
        i, j, k, l, m, n = rdm.elements[idx]
        value = rdm.values[idx]
        dense[i, j, k, l, m, n] += value
        if (i, j, k) != (l, m, n)
            dense[l, m, n, i, j, k] += conj(value)
        end
    end
    return _antisymmetrize_rdm3(dense)
end

function _rdm_representation(representation::Symbol, compact_value, dense_builder::Function)
    if representation == :compact
        return compact_value
    elseif representation == :dense
        return dense_builder()
    else
        throw(ArgumentError("Unsupported RDM representation $(representation). Use :dense or :compact."))
    end
end

function RDM1(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}) where {T}
    _check_coefficients(workspace, coeffs)
    norb = workspace.norb
    rdm1 = zeros(T, norb, norb)
    for i in 1:norb
        for j in 1:i
            for right in eachindex(workspace.hilbert)
                transition = _lookup_rdm1_transition(workspace, workspace.hilbert[right], i, j)
                transition === nothing && continue
                left, sign = transition
                value = _transition_value(coeffs, right, left, sign)
                rdm1[i, j] += value
                if i != j
                    rdm1[j, i] += conj(value)
                end
            end
        end
    end
    return rdm1
end

function RDM2Compact(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}) where {T}
    _check_coefficients(workspace, coeffs)
    values = zeros(T, length(workspace.rdm2_elements))
    nthreads_total = Threads.nthreads()
    thread_values = [zeros(T, length(workspace.rdm2_elements)) for _ in 1:nthreads_total]
    Threads.@threads for t in 1:nthreads_total
        tid = Threads.threadid()
        start_idx, end_idx = _chunk_bounds(length(workspace.hilbert), nthreads_total, t)
        local_values = thread_values[tid]
        @inbounds for right in start_idx:end_idx
            state = workspace.hilbert[right]
            for element_idx in eachindex(workspace.rdm2_elements)
                i, j, k, l = workspace.rdm2_elements[element_idx]
                transition = _lookup_rdm2_transition(workspace, state, i, j, k, l)
                transition === nothing && continue
                left, sign = transition
                local_values[element_idx] += _transition_value(coeffs, right, left, sign)
            end
        end
    end
    for local_values in thread_values
        values .+= local_values
    end
    return CompactRDM2(workspace.norb, workspace.rdm2_elements, values)
end

function RDM2(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}; representation::Symbol=:dense) where {T}
    compact = RDM2Compact(workspace, coeffs)
    return _rdm_representation(representation, compact, () -> todense(compact))
end

function _load_rdm2_cache_payload(file::String)
    return load(file)
end

function _validate_rdm2_cache_payload!(workspace::RDMSectorWorkspace, payload)
    if get(payload, "schema_version", Int32(-1)) != Int32(RDM2_CACHE_SCHEMA_VERSION)
        throw(ArgumentError("RDM2 cache schema version mismatch for workspace $(workspace.kind)."))
    end
    if get(payload, "cache_kind", "") != "rdm2"
        throw(ArgumentError("Unsupported cache payload kind for RDM2 cache file."))
    end
    if get(payload, "model_kind", "") != String(workspace.kind)
        throw(ArgumentError("RDM2 cache model kind does not match the requested workspace."))
    end
    if get(payload, "Nkx", Int32(-1)) != Int32(workspace.Nkx) || get(payload, "Nky", Int32(-1)) != Int32(workspace.Nky)
        throw(ArgumentError("RDM2 cache lattice dimensions do not match the requested workspace."))
    end
    if get(payload, "momentum", Int32(-1)) != Int32(workspace.momentum)
        throw(ArgumentError("RDM2 cache momentum does not match the requested workspace."))
    end
    if get(payload, "nparticle", Int32(-1)) != Int32(workspace.nparticle)
        throw(ArgumentError("RDM2 cache particle number does not match the requested workspace."))
    end
    return nothing
end

function _load_rdm2_cache_entries(workspace::RDMSectorWorkspace, file::String)
    payload = _load_rdm2_cache_payload(file)
    if haskey(payload, "entries")
        _validate_rdm2_cache_payload!(workspace, payload)
        return payload["entries"], false
    end
    return payload["k_$(workspace.momentum)"], true
end

function RDM2Compact(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}, file::String) where {T}
    _check_coefficients(workspace, coeffs)
    values = zeros(T, length(workspace.rdm2_elements))
    element_to_index = Dict{QuadTuple,Int}()
    for idx in eachindex(workspace.rdm2_elements)
        element_to_index[workspace.rdm2_elements[idx]] = idx
    end
    entries, is_legacy_layout = _load_rdm2_cache_entries(workspace, file)
    for (key, value) in entries
        i, j, k, l, right = key
        left, sign = value
        element = is_legacy_layout ? (Int(i), Int(j), Int(l), Int(k)) : (Int(i), Int(j), Int(k), Int(l))
        if haskey(element_to_index, element)
            values[element_to_index[element]] += _transition_value(coeffs, right, left, sign)
        end
    end
    return CompactRDM2(workspace.norb, workspace.rdm2_elements, values)
end

function RDM2(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}, file::String; representation::Symbol=:dense) where {T}
    compact = RDM2Compact(workspace, coeffs, file)
    return _rdm_representation(representation, compact, () -> todense(compact))
end

function rdm2_cache_filename(workspace::RDMSectorWorkspace; dir::String="./data")
    model_tag = String(workspace.kind)
    return joinpath(dir, "rdm2_cache_v$(RDM2_CACHE_SCHEMA_VERSION)_$(model_tag)_n$(workspace.nparticle)_K$(workspace.Nkx)x$(workspace.Nky)_k$(workspace.momentum).jld2")
end

function RDM2_cache(workspace::RDMSectorWorkspace; file::String=rdm2_cache_filename(workspace))
    indtype = valtype(typeof(workspace.ind_dict))
    cache_key_type = NTuple{5,indtype}
    cache_value_type = Tuple{indtype,Int8}
    thread_cache = [Dict{cache_key_type,cache_value_type}() for _ in 1:Threads.nthreads()]
    Threads.@threads for right in eachindex(workspace.hilbert)
        tid = Threads.threadid()
        local_cache = thread_cache[tid]
        state = workspace.hilbert[right]
        for element in workspace.rdm2_elements
            i, j, k, l = element
            transition = _lookup_rdm2_transition(workspace, state, i, j, k, l)
            transition === nothing && continue
            left, sign = transition
            local_cache[(indtype(i), indtype(j), indtype(k), indtype(l), indtype(right))] = (left, sign)
        end
    end
    entries = Dict{cache_key_type,cache_value_type}()
    for local_cache in thread_cache
        merge!(entries, local_cache)
    end
    mkpath(dirname(file))
    save(
        file,
        "schema_version", Int32(RDM2_CACHE_SCHEMA_VERSION),
        "cache_kind", "rdm2",
        "model_kind", String(workspace.kind),
        "Nkx", Int32(workspace.Nkx),
        "Nky", Int32(workspace.Nky),
        "momentum", Int32(workspace.momentum),
        "nparticle", Int32(workspace.nparticle),
        "entries", entries,
        "k_$(workspace.momentum)", entries,
    )
    return file
end

function RDM2(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int, file::String; representation::Symbol=:dense) where {T}
    return RDM2(RDMWorkspace(model, momentum), coeffs, file; representation)
end

function RDM2(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int; representation::Symbol=:dense) where {T}
    return RDM2(RDMWorkspace(model, momentum), coeffs; representation)
end

function RDM2Compact(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM2Compact(RDMWorkspace(model, momentum), coeffs)
end

function RDM2Compact(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int, file::String) where {T}
    return RDM2Compact(RDMWorkspace(model, momentum), coeffs, file)
end

function RDM2_cache(model::ModelParams2DSpinlessList, momentum::Int)
    return RDM2_cache(RDMWorkspace(model, momentum))
end

function RDM1(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM1(RDMWorkspace(model, momentum), coeffs)
end

function RDM1(model::ModelParams2DTwoBand, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM1(RDMWorkspace(model, momentum), coeffs)
end

function RDM2(model::ModelParams2DTwoBand, coeffs::AbstractVector{T}, momentum::Int; representation::Symbol=:dense) where {T}
    return RDM2(RDMWorkspace(model, momentum), coeffs; representation)
end

function RDM2Compact(model::ModelParams2DTwoBand, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM2Compact(RDMWorkspace(model, momentum), coeffs)
end

function RDM3Compact(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}) where {T}
    _check_coefficients(workspace, coeffs)
    if workspace.kind != :spinless_2d
        throw(ArgumentError("RDM3 is currently only implemented for the spinless 2D workspace."))
    end
    values = zeros(T, length(workspace.rdm3_elements))
    nthreads_total = Threads.nthreads()
    thread_values = [zeros(T, length(workspace.rdm3_elements)) for _ in 1:nthreads_total]
    Threads.@threads for t in 1:nthreads_total
        tid = Threads.threadid()
        start_idx, end_idx = _chunk_bounds(length(workspace.hilbert), nthreads_total, t)
        local_values = thread_values[tid]
        @inbounds for right in start_idx:end_idx
            state = workspace.hilbert[right]
            for element_idx in eachindex(workspace.rdm3_elements)
                i, j, k, l, m, n = workspace.rdm3_elements[element_idx]
                transition = _lookup_rdm3_transition(workspace, state, i, j, k, l, m, n)
                transition === nothing && continue
                left, sign = transition
                local_values[element_idx] += _transition_value(coeffs, right, left, sign)
            end
        end
    end
    for local_values in thread_values
        values .+= local_values
    end
    return CompactRDM3(workspace.norb, workspace.rdm3_elements, values)
end

function RDM3(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}; representation::Symbol=:dense) where {T}
    compact = RDM3Compact(workspace, coeffs)
    return _rdm_representation(representation, compact, () -> todense(compact))
end

function RDM3(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int; representation::Symbol=:dense) where {T}
    return RDM3(RDMWorkspace(model, momentum), coeffs; representation)
end

function RDM3Compact(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM3Compact(RDMWorkspace(model, momentum), coeffs)
end

function RDM3_single_thread(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}) where {T}
    _check_coefficients(workspace, coeffs)
    if workspace.kind != :spinless_2d
        throw(ArgumentError("RDM3 is currently only implemented for the spinless 2D workspace."))
    end
    norb = workspace.norb
    rdm3 = zeros(T, norb, norb, norb, norb, norb, norb)
    @inbounds for right in eachindex(workspace.hilbert)
        state = workspace.hilbert[right]
        for element in workspace.rdm3_elements
            i, j, k, l, m, n = element
            transition = _lookup_rdm3_transition(workspace, state, i, j, k, l, m, n)
            transition === nothing && continue
            left, sign = transition
            value = _transition_value(coeffs, right, left, sign)
            rdm3[i, j, k, l, m, n] += value
            if (i, j, k) != (l, m, n)
                rdm3[l, m, n, i, j, k] += conj(value)
            end
        end
    end
    return _antisymmetrize_rdm3(rdm3)
end

function RDM3_single_thread(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM3_single_thread(RDMWorkspace(model, momentum), coeffs)
end

function RDM3_naive(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}) where {T}
    _check_coefficients(workspace, coeffs)
    if workspace.kind != :spinless_2d
        throw(ArgumentError("RDM3 is currently only implemented for the spinless 2D workspace."))
    end
    norb = workspace.norb
    rdm3 = zeros(T, norb, norb, norb, norb, norb, norb)
    for right in eachindex(workspace.hilbert)
        state = workspace.hilbert[right]
        for i in 1:norb, j in 1:norb, k in 1:norb, l in 1:norb, m in 1:norb, n in 1:norb
            transition = _lookup_rdm3_transition(workspace, state, i, j, k, l, m, n)
            transition === nothing && continue
            left, sign = transition
            rdm3[i, j, k, l, m, n] += _transition_value(coeffs, right, left, sign)
        end
    end
    return rdm3
end

function RDM3_naive(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM3_naive(RDMWorkspace(model, momentum), coeffs)
end

function RDM2_naive(workspace::RDMSectorWorkspace, coeffs::AbstractVector{T}) where {T}
    _check_coefficients(workspace, coeffs)
    norb = workspace.norb
    rdm2 = zeros(T, norb, norb, norb, norb)
    for right in eachindex(workspace.hilbert)
        state = workspace.hilbert[right]
        for i in 1:norb, j in 1:norb, k in 1:norb, l in 1:norb
            transition = _lookup_rdm2_transition(workspace, state, i, j, l, k)
            transition === nothing && continue
            left, sign = transition
            rdm2[i, j, l, k] += _transition_value(coeffs, right, left, sign)
        end
    end
    return rdm2
end

function RDM2_naive(model::ModelParams2DSpinlessList, coeffs::AbstractVector{T}, momentum::Int) where {T}
    return RDM2_naive(RDMWorkspace(model, momentum), coeffs)
end

end

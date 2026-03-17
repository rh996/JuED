module FermionOperatorMod

export FermionOperator, CreationOperator!, AnnihilationOperator!
export creation_kernel, annihilation_kernel, apply_operator_string, basis_site_index

"""
    basis_site_index(norb, orbital)

Map a 1-based orbital index from the code's matrix/tensor ordering to the
underlying Fock-basis bit position. JuED stores site `1` in the least
significant bit, so orbital `1` in an `norb`-orbital tensor lives at bit
position `norb`.
"""
@inline basis_site_index(norb::Integer, orbital::Integer) = Int(norb - orbital + 1)

@inline function _max_supported_site(::Type{T}) where {T<:Integer}
    return T <: Signed ? (8 * sizeof(T)) - 1 : 8 * sizeof(T)
end

@inline function _validate_site(::Type{T}, isite::Integer) where {T<:Integer}
    max_site = _max_supported_site(T)
    if isite < 1 || isite > max_site
        throw(ArgumentError("isite=$(isite) is out of range for $(T); supported range is 1:$(max_site)."))
    end
    return nothing
end

@inline function _bitmask(::Type{T}, isite::Integer) where {T<:Integer}
    return one(T) << (isite - 1)
end

@inline function _shifted_state(state::T, isite::Integer) where {T<:Integer}
    return state >> (isite - 1)
end

@inline function creation_kernel(state::T, isite::Integer) where {T<:Integer}
    _validate_site(T, isite)
    shifted = _shifted_state(state, isite)
    if (shifted & one(T)) != zero(T)
        return nothing
    end
    sign = isodd(count_ones(unsigned(shifted))) ? Int8(-1) : Int8(1)
    return state | _bitmask(T, isite), sign
end

@inline function annihilation_kernel(state::T, isite::Integer) where {T<:Integer}
    _validate_site(T, isite)
    shifted = _shifted_state(state, isite)
    if (shifted & one(T)) == zero(T)
        return nothing
    end
    sign = iseven(count_ones(unsigned(shifted))) ? Int8(-1) : Int8(1)
    return state & ~_bitmask(T, isite), sign
end

"""
    apply_operator_string(state, creation_sites, annihilation_sites)

Apply all annihilations first, then all creations, accumulating the fermionic
sign. The `creation_sites` and `annihilation_sites` tuples are already expected
to use the internal bit-position convention described by [`basis_site_index`](@ref).
Returns `nothing` when the operator string annihilates the state.
"""
@inline function apply_operator_string(state::T, creation_sites::Tuple, annihilation_sites::Tuple) where {T<:Integer}
    next_state = state
    sign = Int8(1)

    for site in annihilation_sites
        result = annihilation_kernel(next_state, site)
        result === nothing && return nothing
        next_state, local_sign = result
        sign = Int8(sign * local_sign)
    end

    for site in creation_sites
        result = creation_kernel(next_state, site)
        result === nothing && return nothing
        next_state, local_sign = result
        sign = Int8(sign * local_sign)
    end

    return next_state, sign
end

mutable struct FermionOperator{T<:Integer}
    state::T
    fermion_sign::Int8
end

FermionOperator(state::T, fermion_sign::Integer=1) where {T<:Integer} = FermionOperator{T}(state, Int8(fermion_sign))

function CreationOperator!(fermion_operator::FermionOperator{T}, isite::Integer) where {T<:Integer}
    fermion_operator.fermion_sign == 0 && return fermion_operator
    result = creation_kernel(fermion_operator.state, isite)
    if result === nothing
        fermion_operator.fermion_sign = 0
    else
        fermion_operator.state, local_sign = result
        fermion_operator.fermion_sign = Int8(fermion_operator.fermion_sign * local_sign)
    end
    return fermion_operator
end

function AnnihilationOperator!(fermion_operator::FermionOperator{T}, isite::Integer) where {T<:Integer}
    fermion_operator.fermion_sign == 0 && return fermion_operator
    result = annihilation_kernel(fermion_operator.state, isite)
    if result === nothing
        fermion_operator.fermion_sign = 0
    else
        fermion_operator.state, local_sign = result
        fermion_operator.fermion_sign = Int8(fermion_operator.fermion_sign * local_sign)
    end
    return fermion_operator
end

end

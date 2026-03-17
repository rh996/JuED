module IndexTypesMod

export choose_state_type, choose_pointer_type, state_eltype

function choose_state_type(nbits::Integer)
    if nbits > 31
        return Int64
    end
    return Int32
end

function choose_pointer_type(nentries::Integer)
    if nentries <= typemax(Int32)
        return Int32
    end
    return Int64
end

state_eltype(hilbert) = eltype(hilbert.hilbert)

end

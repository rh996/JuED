
module HamiltonianConstructorMod


export Hamiltonian_list_constructor, Hamiltonian_momentum_constructor


using ..MomentumHilbertSpace2DMod
using ..FermionOperatorMod
using ..IndexTypesMod: choose_pointer_type
using SparseArrays

AbstractHilbertSpace = MomentumHilbertSpace2DMod.AbstractHilbertSpace
ToDict = MomentumHilbertSpace2DMod.ToDict
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
    # return data, row, column
end

@inline function _chunk_bounds(t::Int, chunk_size::Int, total_size::Int)
    start_idx = (t - 1) * chunk_size + 1
    end_idx = min(t * chunk_size, total_size)
    return start_idx, end_idx
end

function _assemble_hamiltonian_csc(::Type{T}, hilbertspace, count_column!, fill_column!) where {T}
    hilbert = hilbertspace.hilbert
    dim_hilbert = length(hilbert)
    nthreads_total = Threads.nthreads()
    chunk_size = ceil(Int, dim_hilbert / nthreads_total)
    nnz_col = zeros(Int, dim_hilbert)

    Threads.@threads for t in 1:nthreads_total
        start_idx, end_idx = _chunk_bounds(t, chunk_size, dim_hilbert)
        @inbounds for state_ind in start_idx:end_idx
            nnz_col[state_ind] = count_column!(state_ind)
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
            fill_column!(state_ind, indptr[state_ind], data, row, ind_dict)
        end
    end

    return data, row, indptr, dim_hilbert
end

function Hamiltonian_momentum_constructor(OneBody, Eri::Array{T,6}, hilbertspace) where {T}
    onebodyind = findall(x -> abs(x) > 1e-10, OneBody)
    hilbert = hilbertspace.hilbert
    norbital = size(OneBody)[1]
    Nkx, Nky = hilbertspace.Nkx, hilbertspace.Nky
    Nk = Nkx * Nky

    function count_column!(state_ind)
        local_count = 0

        for ind in onebodyind
            i = ind[1]
            j = ind[2]
            inversei = (norbital - i) + 1
            inversej = (norbital - j) + 1

            fermion = FermionOperator(hilbert[state_ind], 1)
            AnnihilationOperator!(fermion, inversej)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversei)
            if fermion.fermion_sign == 0
                continue
            end

            local_count += 1
        end

        for ik in 0:(Nk-1)
            for ikp in 0:(Nk-1)
                for iq in 0:(Nk-1)
                    ind1 = norbital - momentum_add_2d(ik, iq, Nkx, Nky)
                    ind2 = norbital - ik
                    ind3 = norbital - momentum_sub_2d(ikp, iq, Nkx, Nky)
                    ind4 = norbital - ikp

                    fermion = FermionOperator(hilbert[state_ind], 1)
                    AnnihilationOperator!(fermion, ind2)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    AnnihilationOperator!(fermion, ind4)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind3)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind1)
                    if fermion.fermion_sign == 0
                        continue
                    end

                    local_count += 1
                end
            end
        end

        return local_count
    end

    function fill_column!(state_ind, writepos, data, row, ind_dict)
        for ind in onebodyind
            i = ind[1]
            j = ind[2]
            inversei = (norbital - i) + 1
            inversej = (norbital - j) + 1
            fermion = FermionOperator(hilbert[state_ind], 1)

            AnnihilationOperator!(fermion, inversej)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversei)
            if fermion.fermion_sign == 0
                continue
            end

            data[writepos] = OneBody[i, j] * fermion.fermion_sign
            row[writepos] = ind_dict[fermion.state]
            writepos += 1
        end

        for ik in 0:(norbital-1)
            for ikp in 0:(norbital-1)
                for iq in 0:(norbital-1)
                    ind1 = norbital - momentum_add_2d(ik, iq, Nkx, Nky)
                    ind2 = norbital - ik
                    ind3 = norbital - momentum_sub_2d(ikp, iq, Nkx, Nky)
                    ind4 = norbital - ikp

                    fermion = FermionOperator(hilbert[state_ind], 1)
                    AnnihilationOperator!(fermion, ind2)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    AnnihilationOperator!(fermion, ind4)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind3)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind1)
                    if fermion.fermion_sign == 0
                        continue
                    end

                    ikx = mod(ik, Nkx)
                    iky = fld(ik, Nkx)
                    ikpx = mod(ikp, Nkx)
                    ikpy = fld(ikp, Nkx)
                    iqx = mod(iq, Nkx)
                    iqy = fld(iq, Nkx)

                    data[writepos] = Eri[ikx+1, iky+1, ikpx+1, ikpy+1, iqx+1, iqy+1] * fermion.fermion_sign
                    row[writepos] = ind_dict[fermion.state]
                    writepos += 1
                end
            end
        end

        return nothing
    end

    return _assemble_hamiltonian_csc(T, hilbertspace, count_column!, fill_column!)
end











function Hamiltonian_list_constructor(OneBody, Eri::Array{T,4}, hilbertspace) where {T}
    onebodyind = findall(x -> abs(x) > 1e-10, OneBody)
    eriind = findall(x -> abs(x) > 1e-10, Eri)
    hilbert = hilbertspace.hilbert
    norbital = size(OneBody)[1]
    function count_column!(state_ind)
        local_count = 0

        for ind1 in onebodyind
            i = ind1[1]
            j = ind1[2]
            fermion = FermionOperator(hilbert[state_ind], 1)
            inversei = (norbital - i) + 1
            inversej = (norbital - j) + 1

            AnnihilationOperator!(fermion, inversej)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversei)
            if fermion.fermion_sign == 0
                continue
            end

            local_count += 1
        end

        for ind2 in eriind
            i = ind2[1]
            j = ind2[2]
            k = ind2[3]
            l = ind2[4]

            inversei = (norbital - i) + 1
            inversej = (norbital - j) + 1
            inversek = (norbital - k) + 1
            inversel = (norbital - l) + 1

            fermion = FermionOperator(hilbert[state_ind], 1)
            AnnihilationOperator!(fermion, inversek)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, inversel)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversej)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversei)
            if fermion.fermion_sign == 0
                continue
            end

            local_count += 1
        end

        return local_count
    end

    function fill_column!(state_ind, writepos, data, row, ind_dict)
        for ind1 in onebodyind
            i = ind1[1]
            j = ind1[2]
            inversei = (norbital - i) + 1
            inversej = (norbital - j) + 1
            fermion = FermionOperator(hilbert[state_ind], 1)

            AnnihilationOperator!(fermion, inversej)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversei)
            if fermion.fermion_sign == 0
                continue
            end

            data[writepos] = OneBody[i, j] * fermion.fermion_sign
            row[writepos] = ind_dict[fermion.state]
            writepos += 1
        end

        for ind2 in eriind
            i = ind2[1]
            j = ind2[2]
            k = ind2[3]
            l = ind2[4]

            inversei = (norbital - i) + 1
            inversej = (norbital - j) + 1
            inversek = (norbital - k) + 1
            inversel = (norbital - l) + 1

            fermion = FermionOperator(hilbert[state_ind], 1)
            AnnihilationOperator!(fermion, inversek)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, inversel)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversej)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, inversei)
            if fermion.fermion_sign == 0
                continue
            end

            data[writepos] = Eri[i, j, k, l] * fermion.fermion_sign
            row[writepos] = ind_dict[fermion.state]
            writepos += 1
        end

        return nothing
    end

    return _assemble_hamiltonian_csc(T, hilbertspace, count_column!, fill_column!)

end



end

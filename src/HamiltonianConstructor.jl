
module HamiltonianConstructorMod


export Hamiltonian_list_constructor, Hamiltonian_momentum_constructor


using ProgressMeter
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



function Hamiltonian_momentum_constructor(OneBody, Eri::Array{T,6}, hilbertspace) where {T}


    onebodyind = findall(x -> abs(x) > 1e-10, OneBody)

    # Build your spin Hilbert space
    hilbert = hilbertspace.hilbert
    DimHilbert = length(hilbert)
    norbital = size(OneBody)[1]
    Nkx, Nky = hilbertspace.Nkx, hilbertspace.Nky
    Nk = Nkx * Nky
    # Number of threads and chunk size for parallel loop.
    nthreads_total = Threads.nthreads()
    chunk_size = ceil(Int, DimHilbert / nthreads_total)

    # We'll store the nonzeros-per-column here.
    nnz_col = zeros(Int, DimHilbert)

    # Optional progress bar, showing total columns processed
    p = Progress(DimHilbert, desc="Counting nnz_col")

    # 1) ------------------- First pass: compute nnz_col[state_ind] -------------------
    Threads.@threads for t in 1:nthreads_total
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)

        @inbounds for state_ind in start_idx:end_idx
            next!(p)  # update progress bar
            local_count = 0

            # ==================== Kinetic part ====================
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



            # ==================== Interaction part =================
            for ik in 0:(Nk-1)
                for ikp in 0:(Nk-1)
                    for iq in 0:(Nk-1)
                        ind1 = momentum_add_2d(ik, iq, Nkx, Nky)
                        ind1 = norbital - ind1
                        ind2 = ik
                        ind2 = norbital - ind2
                        ind3 = momentum_sub_2d(ikp, iq, Nkx, Nky)
                        ind3 = norbital - ind3
                        ind4 = ikp
                        ind4 = norbital - ind4

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

            nnz_col[state_ind] = local_count
        end
    end

    # 2) ------------------- Build indptr array via prefix sum of nnz_col --------------
    # By CSC convention, indptr has length = DimHilbert + 1,
    # with indptr[1] = 1, indptr[col+1] = indptr[col] + nnz_col[col].
    pointertype = choose_pointer_type(sum(nnz_col) + 1)
    indptr = Array{pointertype}(undef, DimHilbert + 1)
    indptr[1] = 1
    for col in 1:DimHilbert
        indptr[col+1] = indptr[col] + nnz_col[col]
    end

    total_nnz = indptr[end] - 1
    # If you want the same floating type as U, do `Array{typeof(U)}`:
    data = Array{T}(undef, total_nnz)
    row = Array{pointertype}(undef, total_nnz)
    ind_dict = ToDict(hilbertspace, pointertype)
    # 3) ------------------- Second pass: fill in data[] and row[] in parallel ----------
    # We'll do the same chunk approach, but this time we know exactly where each
    # column's data belongs: starting at indptr[col].
    p = Progress(DimHilbert, desc="Filling CSC data")
    Threads.@threads for t in 1:nthreads_total
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)

        @inbounds for state_ind in start_idx:end_idx
            next!(p)  # update progress
            writepos = indptr[state_ind]  # start of this column in data/row

            # ==================== Kinetic part ====================
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

                # Kinetic energy factor

                eps = OneBody[i, j]

                data[writepos] = eps * fermion.fermion_sign
                row[writepos] = ind_dict[fermion.state]
                writepos += 1
            end


            # ==================== Interaction part =================
            for ik in 0:(norbital-1)
                for ikp in 0:(norbital-1)
                    for iq in 0:(norbital-1)
                        ind1 = momentum_add_2d(ik, iq, Nkx, Nky)
                        ind1 = norbital - ind1
                        ind2 = ik
                        ind2 = norbital - ind2
                        ind3 = momentum_sub_2d(ikp, iq, Nkx, Nky)
                        ind3 = norbital - ind3
                        ind4 = ikp
                        ind4 = norbital - ind4

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
        end
    end

    return data, row, indptr, DimHilbert
end












function Hamiltonian_list_constructor(OneBody, Eri::Array{T,4}, hilbertspace) where {T}


    onebodyind = findall(x -> abs(x) > 1e-10, OneBody)
    eriind = findall(x -> abs(x) > 1e-10, Eri)


    # Build Hilbert space
    hilbert = hilbertspace.hilbert
    DimHilbert = length(hilbert)
    norbital = size(OneBody)[1]
    # Nkx, Nky     = hilbertspace.Nkx, hilbertspace.Nky

    # Number of threads and chunk size for parallel loop.
    nthreads_total = Threads.nthreads()
    chunk_size = ceil(Int, DimHilbert / nthreads_total)

    # We'll store the nonzeros-per-column here.
    nnz_col = zeros(Int, DimHilbert)

    # Optional progress bar, showing total columns processed
    p = Progress(DimHilbert, desc="Counting nnz_col")

    # 1) ------------------- First pass: compute nnz_col[state_ind] -------------------
    Threads.@threads for t in 1:nthreads_total
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)

        @inbounds for state_ind in start_idx:end_idx
            next!(p)  # update progress bar
            local_count = 0

            # ==================== Kinetic part ====================
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


            # ==================== Interaction part =================
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

            nnz_col[state_ind] = local_count
        end
    end

    # 2) ------------------- Build indptr array via prefix sum of nnz_col --------------
    # By CSC convention, indptr has length = DimHilbert + 1,
    # with indptr[1] = 1, indptr[col+1] = indptr[col] + nnz_col[col].
    totaldata = sum(nnz_col)
    pointertype = choose_pointer_type(totaldata + 1)


    indptr = Array{pointertype}(undef, DimHilbert + 1)
    indptr[1] = 1
    for col in 1:DimHilbert
        indptr[col+1] = indptr[col] + nnz_col[col]
    end

    total_nnz = indptr[end] - 1
    # If you want the same floating type as U, do `Array{typeof(Eri)}`:
    data = Array{T}(undef, total_nnz)


    row = Array{pointertype}(undef, total_nnz)
    ind_dict = ToDict(hilbertspace, pointertype)
    # 3) ------------------- Second pass: fill in data[] and row[] in parallel ----------
    # We'll do the same chunk approach, but this time we know exactly where each
    # column's data belongs: starting at indptr[col].
    p = Progress(DimHilbert, desc="Filling CSC data")
    Threads.@threads for t in 1:nthreads_total
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)

        @inbounds for state_ind in start_idx:end_idx
            next!(p)  # update progress
            writepos = indptr[state_ind]  # start of this column in data/row

            # ==================== Kinetic part ====================
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

                # Kinetic energy factor

                eps = OneBody[i, j]

                data[writepos] = eps * fermion.fermion_sign
                row[writepos] = ind_dict[fermion.state]
                writepos += 1
            end

            # ==================== Interaction part =================
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
        end
    end

    return data, row, indptr, DimHilbert

end



end

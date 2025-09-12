include("../src/FermionOperator.jl")
include("../src/SpinMomentumHilbertSpace2D.jl")

using .SpinMomentumHilbertSpace2DMod
using SparseArrays
using KrylovKit
using Arpack
using ProgressMeter
using .FermionOperatorMod

nalpha = 6
nbeta = 6
Nkx = 2
Nky = 6

AbstractHilbertSpace = SpinMomentumHilbertSpace2DMod.AbstractHilbertSpace


function count_data_size(hilbertspace::AbstractHilbertSpace)
    hilbert = SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)
    m = length(hilbert)
    ptrlenth = m + 1
    datalenth = 0
    norbital = hilbertspace.Nkx * hilbertspace.Nky
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    p = Progress(m)
    @inbounds for state_ind in eachindex(hilbert)
        next!(p)
        for i in 0:norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = norbital - i
            position = 2 * position
            AnnihilationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end

            datalenth += 1

        end

        for i in 0:norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = norbital - i
            position = 2 * position - 1
            AnnihilationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end

            datalenth += 1
        end

        for ik in 0:norbital-1
            for ikp in 0:norbital-1
                for iq in 0:norbital-1
                    ind1 = 2 * momentum_add_2d(ik, iq, Nkx, Nky)
                    ind1 = 2 * norbital - ind1
                    ind2 = 2 * ik
                    ind2 = 2 * norbital - ind2
                    ind3 = 2 * momentum_sub_2d(ikp, iq, Nkx, Nky) + 1
                    ind3 = 2 * norbital - ind3
                    ind4 = 2 * ikp + 1
                    ind4 = 2 * norbital - ind4
                    fermion = FermionOperator(hilbert[state_ind], 1)
                    AnnihilationOperator!(fermion, ind4)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind3)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    AnnihilationOperator!(fermion, ind2)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind1)
                    if fermion.fermion_sign == 0
                        continue
                    end

                    datalenth += 1
                end
            end
        end

    end

    return datalenth, ptrlenth
end


function count_data_size_paralell(hilbertspace::AbstractHilbertSpace)
    nthreads_total = Threads.nthreads()
    partial_sums = Vector{Int64}(undef, nthreads_total)


    hilbert = SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)
    m = length(hilbert)
    ptrlenth = m + 1
    chunk_size = ceil(Int64, m / nthreads_total)


    norbital = hilbertspace.Nkx * hilbertspace.Nky
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    p = Progress(m)


    Threads.@threads for t in 1:nthreads_total
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, m)

        local_sum = zero(Int64)

        @inbounds for state_ind in start_idx:end_idx
            next!(p)
            for i in 0:norbital-1

                fermion = FermionOperator(hilbert[state_ind], 1)
                position = norbital - i
                position = 2 * position
                AnnihilationOperator!(fermion, position)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, position)
                if fermion.fermion_sign == 0
                    continue
                end

                local_sum += 1

            end

            for i in 0:norbital-1

                fermion = FermionOperator(hilbert[state_ind], 1)
                position = norbital - i
                position = 2 * position - 1
                AnnihilationOperator!(fermion, position)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, position)
                if fermion.fermion_sign == 0
                    continue
                end

                local_sum += 1
            end

            for ik in 0:norbital-1
                for ikp in 0:norbital-1
                    for iq in 0:norbital-1
                        ind1 = 2 * momentum_add_2d(ik, iq, Nkx, Nky)
                        ind1 = 2 * norbital - ind1
                        ind2 = 2 * ik
                        ind2 = 2 * norbital - ind2
                        ind3 = 2 * momentum_sub_2d(ikp, iq, Nkx, Nky) + 1
                        ind3 = 2 * norbital - ind3
                        ind4 = 2 * ikp + 1
                        ind4 = 2 * norbital - ind4
                        fermion = FermionOperator(hilbert[state_ind], 1)
                        AnnihilationOperator!(fermion, ind4)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        CreationOperator!(fermion, ind3)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        AnnihilationOperator!(fermion, ind2)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        CreationOperator!(fermion, ind1)
                        if fermion.fermion_sign == 0
                            continue
                        end

                        local_sum += 1
                    end
                end
            end
        end
        partial_sums[t] = local_sum
    end

    # Finally, sum up the partial sums in a single pass (main thread).
    total_sum = zero(Int64)
    for t in 1:nthreads_total
        total_sum += partial_sums[t]
    end

    return total_sum, ptrlenth
end


function Hubbard2DMatrix(U::Real, datasize::Int64, hilbertspace::AbstractHilbertSpace)
    hilbert = hilbertspace.hilbert
    ind_dict = SpinMomentumHilbertSpace2DMod.ToDict(hilbertspace)
    DimHilbert = length(hilbert)
    # println("Sector Dim: ", m)
    norbital = hilbertspace.Nkx * hilbertspace.Nky
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky


    datatp = typeof(U)
    data = Array{Float64}(undef, datasize)
    row = Array{Int32}(undef, datasize)
    indptr = Array{Int32}(undef, DimHilbert + 1)
    indptr[1] = 1
    COUNT = 0


    @inbounds for state_ind in eachindex(hilbert)
        ncolcount = 0
        for i in 0:norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = norbital - i
            position = 2 * position
            AnnihilationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end
            COUNT += 1
            ikx = mod(i, Nkx)
            iky = fld(i, Nkx)
            data[COUNT] = (-2 * cos(2 * pi * ikx / Nkx) - 2 * cos(2 * pi * iky / Nky)) * fermion.fermion_sign

            row[COUNT] = ind_dict[fermion.state]
            ncolcount += 1

        end

        for i in 0:norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = norbital - i
            position = 2 * position - 1
            AnnihilationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, position)
            if fermion.fermion_sign == 0
                continue
            end
            COUNT += 1
            ikx = mod(i, Nkx)
            iky = fld(i, Nkx)
            data[COUNT] = (-2 * cos(2 * pi * ikx / Nkx) - 2 * cos(2 * pi * iky / Nky)) * fermion.fermion_sign
            row[COUNT] = ind_dict[fermion.state]
            ncolcount += 1
        end

        for ik in 0:norbital-1
            for ikp in 0:norbital-1
                for iq in 0:norbital-1
                    ind1 = 2 * momentum_add_2d(ik, iq, Nkx, Nky)
                    ind1 = 2 * norbital - ind1
                    ind2 = 2 * ik
                    ind2 = 2 * norbital - ind2
                    ind3 = 2 * momentum_sub_2d(ikp, iq, Nkx, Nky) + 1
                    ind3 = 2 * norbital - ind3
                    ind4 = 2 * ikp + 1
                    ind4 = 2 * norbital - ind4
                    fermion = FermionOperator(hilbert[state_ind], 1)
                    AnnihilationOperator!(fermion, ind4)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind3)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    AnnihilationOperator!(fermion, ind2)
                    if fermion.fermion_sign == 0
                        continue
                    end
                    CreationOperator!(fermion, ind1)
                    if fermion.fermion_sign == 0
                        continue
                    end

                    COUNT += 1
                    data[COUNT] = (U / norbital) * fermion.fermion_sign
                    row[COUNT] = ind_dict[fermion.state]
                    ncolcount += 1
                end
            end
        end
        # push!(indptr, ncolcount+indptr[end])
        indptr[state_ind+1] = ncolcount + indptr[state_ind]
    end

    return data, row, indptr, DimHilbert

end





function buildHubbard2DMatrixCSC_parallel(U::Real, hilbertspace)
    # Build your spin Hilbert space (as you did in count_data_size_paralell).
    hilbert = SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)
    DimHilbert = length(hilbert)
    norbital = hilbertspace.Nkx * hilbertspace.Nky
    Nkx, Nky = hilbertspace.Nkx, hilbertspace.Nky

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

            # ==================== Kinetic (spin-up) part ====================
            for i in 0:norbital-1
                fermion = FermionOperator(hilbert[state_ind], 1)
                pos = 2 * (norbital - i)
                AnnihilationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end

                local_count += 1
            end

            # ==================== Kinetic (spin-down) part ==================
            for i in 0:norbital-1
                fermion = FermionOperator(hilbert[state_ind], 1)
                pos = 2 * (norbital - i) - 1
                AnnihilationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end

                local_count += 1
            end

            # ==================== Interaction (U-term) part =================
            for ik in 0:(norbital-1)
                for ikp in 0:(norbital-1)
                    for iq in 0:(norbital-1)
                        ind1 = 2 * momentum_add_2d(ik, iq, Nkx, Nky)
                        ind1 = 2 * norbital - ind1
                        ind2 = 2 * ik
                        ind2 = 2 * norbital - ind2
                        ind3 = 2 * momentum_sub_2d(ikp, iq, Nkx, Nky) + 1
                        ind3 = 2 * norbital - ind3
                        ind4 = 2 * ikp + 1
                        ind4 = 2 * norbital - ind4

                        fermion = FermionOperator(hilbert[state_ind], 1)
                        AnnihilationOperator!(fermion, ind4)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        CreationOperator!(fermion, ind3)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        AnnihilationOperator!(fermion, ind2)
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
    indptr = Array{Int32}(undef, DimHilbert + 1)
    indptr[1] = 1
    for col in 1:DimHilbert
        indptr[col+1] = indptr[col] + nnz_col[col]
    end

    total_nnz = indptr[end] - 1
    # If you want the same floating type as U, do `Array{typeof(U)}`:
    data = Array{Float64}(undef, total_nnz)
    row = Array{Int32}(undef, total_nnz)
    ind_dict = SpinMomentumHilbertSpace2DMod.ToDict(hilbertspace)
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

            # ==================== Kinetic (spin-up) part ====================
            for i in 0:norbital-1
                fermion = FermionOperator(hilbert[state_ind], 1)
                pos = 2 * (norbital - i)
                AnnihilationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end

                # Kinetic energy factor
                ikx = mod(i, Nkx)
                iky = fld(i, Nkx)
                eps = -2 * cos(2π * ikx / Nkx) - 2 * cos(2π * iky / Nky)

                data[writepos] = eps * fermion.fermion_sign
                row[writepos] = Int32(ind_dict[fermion.state])
                writepos += 1
            end

            # ==================== Kinetic (spin-down) part ==================
            for i in 0:norbital-1
                fermion = FermionOperator(hilbert[state_ind], 1)
                pos = 2 * (norbital - i) - 1
                AnnihilationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, pos)
                if fermion.fermion_sign == 0
                    continue
                end

                ikx = mod(i, Nkx)
                iky = fld(i, Nkx)
                eps = -2 * cos(2π * ikx / Nkx) - 2 * cos(2π * iky / Nky)

                data[writepos] = eps * fermion.fermion_sign
                row[writepos] = Int32(ind_dict[fermion.state])
                writepos += 1
            end

            # ==================== Interaction (U-term) part =================
            for ik in 0:(norbital-1)
                for ikp in 0:(norbital-1)
                    for iq in 0:(norbital-1)
                        ind1 = 2 * momentum_add_2d(ik, iq, Nkx, Nky)
                        ind1 = 2 * norbital - ind1
                        ind2 = 2 * ik
                        ind2 = 2 * norbital - ind2
                        ind3 = 2 * momentum_sub_2d(ikp, iq, Nkx, Nky) + 1
                        ind3 = 2 * norbital - ind3
                        ind4 = 2 * ikp + 1
                        ind4 = 2 * norbital - ind4

                        fermion = FermionOperator(hilbert[state_ind], 1)
                        AnnihilationOperator!(fermion, ind4)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        CreationOperator!(fermion, ind3)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        AnnihilationOperator!(fermion, ind2)
                        if fermion.fermion_sign == 0
                            continue
                        end
                        CreationOperator!(fermion, ind1)
                        if fermion.fermion_sign == 0
                            continue
                        end

                        data[writepos] = (U / norbital) * fermion.fermion_sign
                        row[writepos] = Int32(ind_dict[fermion.state])
                        writepos += 1
                    end
                end
            end
        end
    end

    return data, row, indptr, DimHilbert
end






function test(nalpha, nbeta, Nkx, Nky)
    elist = []
    for k in 0:Nkx*Nky-1
        hilbertspace = SpinMomentumHilbertSpace2DMod.SpinMomentumHilbertSpace2D{Int64}(nalpha, nbeta, Nkx, Nky, k, [])
        # @show datasize,ptrsize = count_data_size(hilbertspace)

        # @show datasize,ptrsize = count_data_size_paralell(hilbertspace)
        # data,row,indptr, dim = Hubbard2DMatrix(8.0,datasize,hilbertspace)

        data, row, indptr, dim = buildHubbard2DMatrixCSC_parallel(4.0, hilbertspace)

        @show dim
        # vals,vecs = eigs(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), nev=5, which=:SR)
        vals, vecs, info = eigsolve(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), ones(Float64, dim), 10, :SR, Lanczos())
        push!(elist, vals[1:5])

    end
    elist
end



@time test(nalpha, nbeta, Nkx, Nky)

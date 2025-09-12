module DensityMatricesMod
using ..FermionOperatorMod
using ..ModelTypesMod
using ..MomentumHilbertSpace2DMod
using ..TwoBandMomentumHilbertSpace2DMod
using JLD2
export RDM1, RDM2, RDM3, RDM3_single_thread, RDM2_cache, RDM3_naive, RDM2_naive

function RDM1(ModelParams2DSpinless::ModelParams2DSpinlessList, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky
    if Nkx * Nky > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams2DSpinless.nparticle

    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)
    ind_dict = MomentumHilbertSpace2DMod.ToDict(hilbertspace)

    rdm1 = zeros(T, norb, norb)

    for i in 0:norb-1
        for j in 0:i
            ireverse = norb - i
            jreverse = norb - j
            for k in eachindex(hilbert)

                state = hilbert[k]
                fermion = FermionOperator(state, 1)

                AnnihilationOperator!(fermion, jreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, ireverse)
                if fermion.fermion_sign == 0
                    continue
                end
                # println(fermion.state,ind_dict[fermion.state])
                if haskey(ind_dict, fermion.state)
                    left = ind_dict[fermion.state]
                else
                    continue
                end

                value = fermion.fermion_sign * coeffs[k] * conj(coeffs[left])
                rdm1[i+1, j+1] += value

                if i != j
                    rdm1[j+1, i+1] += conj(value)
                end

            end

        end
    end

    return rdm1

end

function RDM1(ModelParams::ModelParams2DTwoBand, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams.Nkx
    Nky = ModelParams.Nky
    norb = Nkx * Nky * 2
    if norb > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams.nparticle

    hilbertspace = TwoBandMomentumHilbertSpace2DMod.TwoBandMomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = TwoBandMomentumHilbertSpace2DMod.BuildTwoBandHilbert(hilbertspace)
    ind_dict = TwoBandMomentumHilbertSpace2DMod.ToDict(hilbertspace)

    rdm1 = zeros(T, norb, norb)

    for i in 0:norb-1
        for j in 0:i
            ireverse = norb - i
            jreverse = norb - j
            for k in eachindex(hilbert)

                state = hilbert[k]
                fermion = FermionOperator(state, 1)

                AnnihilationOperator!(fermion, jreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, ireverse)
                if fermion.fermion_sign == 0
                    continue
                end
                # println(fermion.state,ind_dict[fermion.state])
                if haskey(ind_dict, fermion.state)
                    left = ind_dict[fermion.state]
                else
                    continue
                end

                value = fermion.fermion_sign * coeffs[k] * conj(coeffs[left])
                rdm1[i+1, j+1] += value

                if i != j
                    rdm1[j+1, i+1] += conj(value)
                end

            end

        end
    end

    return rdm1

end

function RDM2(ModelParams2DSpinless::ModelParams2DSpinlessList, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky
    if Nkx * Nky > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams2DSpinless.nparticle

    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)
    ind_dict = MomentumHilbertSpace2DMod.ToDict(hilbertspace)
    pairs = []

    for i in 1:norb
        for j in i+1:norb
            push!(pairs, (i, j))
        end
    end

    independent_elements = []
    number_of_pairs = length(pairs)

    for ind in eachindex(pairs)
        i, j = pairs[ind]
        for ind2 in ind:number_of_pairs
            k, l = pairs[ind2]
            if momentum_add_2d(i - 1, j - 1, Nkx, Nky) == momentum_add_2d(l - 1, k - 1, Nkx, Nky)
                push!(independent_elements, (i, j, k, l))
            end
        end
    end

    # @show independent_elements
    # @show length(independent_elements)

    rdm2_final = zeros(T, norb, norb, norb, norb)


    DimHilbert = length(hilbert)

    nthreads_total = Threads.nthreads()
    chunk_size = ceil(Int, DimHilbert / nthreads_total)

    rdm2_thread = [zeros(T, norb, norb, norb, norb) for _ in 1:nthreads_total]

    Threads.@threads for t in 1:nthreads_total
        tid = Threads.threadid()
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)
        local_rdm2 = rdm2_thread[tid]
        for m in start_idx:end_idx
            state = hilbert[m]
            for ind in eachindex(independent_elements)
                i, j, l, k = independent_elements[ind]
                ireverse = norb - i + 1
                jreverse = norb - j + 1
                kreverse = norb - k + 1
                lreverse = norb - l + 1


                fermion = FermionOperator(state, 1)

                AnnihilationOperator!(fermion, lreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                AnnihilationOperator!(fermion, kreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, jreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, ireverse)
                if fermion.fermion_sign == 0
                    continue
                end

                left = ind_dict[fermion.state]

                value = fermion.fermion_sign * coeffs[m] * conj(coeffs[left])
                local_rdm2[i, j, l, k] += value

                if i != l || j != k
                    local_rdm2[l, k, i, j] += conj(value)
                end

            end


        end

    end

    # Combine the thread-local arrays into the final result.
    for t in 1:nthreads_total
        rdm2_final .+= rdm2_thread[t]
    end

    rdm2_final -= permutedims(rdm2_final, (2, 1, 3, 4))
    rdm2_final -= permutedims(rdm2_final, (1, 2, 4, 3))

    return rdm2_final

end


function RDM2(ModelParams2DSpinless::ModelParams2DSpinlessList, coeffs::Array{T,1}, momentum::Int, file::String) where {T}
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky


    ijkl_cache = load(file)["k_$(momentum)"]

    rdm2 = zeros(T, norb, norb, norb, norb)
    for (key, value) in ijkl_cache
        i, j, k, l, right = key
        left, sign = value

        value = sign * coeffs[right] * conj(coeffs[left])
        rdm2[i, j, l, k] += value
        if (i, j) != (l, k)
            rdm2[l, k, i, j] += conj(value)
        end
    end


    rdm2 -= permutedims(rdm2, (2, 1, 3, 4))
    rdm2 -= permutedims(rdm2, (1, 2, 4, 3))
    return rdm2

end

function RDM2_cache(ModelParams2DSpinless::ModelParams2DSpinlessList, momentum::Int)

    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky
    if Nkx * Nky > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams2DSpinless.nparticle

    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)
    ind_dict = MomentumHilbertSpace2DMod.ToDict(hilbertspace)
    pairs = []

    for i in 1:norb
        for j in i+1:norb
            push!(pairs, (i, j))
        end
    end

    independent_elements = []
    number_of_pairs = length(pairs)

    for ind in eachindex(pairs)
        i, j = pairs[ind]
        for ind2 in ind:number_of_pairs
            k, l = pairs[ind2]
            if momentum_add_2d(i - 1, j - 1, Nkx, Nky) == momentum_add_2d(l - 1, k - 1, Nkx, Nky)
                push!(independent_elements, (i, j, k, l))
            end

        end
    end




    ijkl_cache = Dict{NTuple{5,Int32}}{NTuple{2,Int32}}()


    # Create a thread-local dictionary for each thread.
    thread_cache = [Dict{NTuple{5,Int32},NTuple{2,Int32}}() for _ in 1:Threads.nthreads()]

    Threads.@threads for m in eachindex(hilbert)
        tid = Threads.threadid()
        local_cache = thread_cache[tid]
        state = hilbert[m]
        for ind in eachindex(independent_elements)
            i, j, l, k = independent_elements[ind]




            ireverse = norb - i + 1
            jreverse = norb - j + 1
            kreverse = norb - k + 1
            lreverse = norb - l + 1

            fermion = FermionOperator(state, 1)

            AnnihilationOperator!(fermion, lreverse)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, kreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, jreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, ireverse)
            if fermion.fermion_sign == 0
                continue
            end


            left = ind_dict[fermion.state]


            local_cache[(i, j, k, l, m)] = (left, fermion.fermion_sign)
        end
    end

    # Merge thread-local dictionaries into the global cache.
    for local_cache in thread_cache
        merge!(ijkl_cache, local_cache)
    end

    save("./data/rdm2_chahe_K$(Nkx)x$(Nky).jld2", "k_$(momentum)", ijkl_cache)
end

function RDM2(ModelParams::ModelParams2DTwoBand, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams.Nkx
    Nky = ModelParams.Nky

    norb = Nkx * Nky * 2
    if norb > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams.nparticle

    hilbertspace = TwoBandMomentumHilbertSpace2DMod.TwoBandMomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = TwoBandMomentumHilbertSpace2DMod.BuildTwoBandHilbert(hilbertspace)
    ind_dict = TwoBandMomentumHilbertSpace2DMod.ToDict(hilbertspace)
    pairs = []

    for i in 1:norb
        for j in i+1:norb
            push!(pairs, (i, j))
        end
    end

    independent_elements = []
    number_of_pairs = length(pairs)

    for ind in eachindex(pairs)
        i, j = pairs[ind]
        for ind2 in ind:number_of_pairs
            k, l = pairs[ind2]
            if momentum_add_2d(fld(i - 1, 2), fld(j - 1, 2), Nkx, Nky) == momentum_add_2d(fld(l - 1, 2), fld(k - 1, 2), Nkx, Nky)
                push!(independent_elements, (i, j, k, l))
            end
            # push!(independent_elements, (i,j,k,l))
        end
    end

    # @show independent_elements
    @show length(independent_elements)


    rdm2_final = zeros(T, norb, norb, norb, norb)


    DimHilbert = length(hilbert)

    nthreads_total = Threads.nthreads()
    chunk_size = ceil(Int, DimHilbert / nthreads_total)

    rdm2_thread = [zeros(T, norb, norb, norb, norb) for _ in 1:nthreads_total]

    Threads.@threads for t in 1:nthreads_total
        tid = Threads.threadid()
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)
        local_rdm2 = rdm2_thread[tid]
        for m in start_idx:end_idx
            state = hilbert[m]
            for ind in eachindex(independent_elements)
                i, j, k, l = independent_elements[ind]

                ireverse = norb - (i - 1)
                jreverse = norb - (j - 1)
                kreverse = norb - (k - 1)
                lreverse = norb - (l - 1)


                fermion = FermionOperator(state, 1)

                AnnihilationOperator!(fermion, kreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                AnnihilationOperator!(fermion, lreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, jreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, ireverse)
                if fermion.fermion_sign == 0
                    continue
                end


                left = ind_dict[fermion.state]


                value = fermion.fermion_sign * coeffs[m] * conj(coeffs[left])
                local_rdm2[i, j, k, l] += value

                if (i, j) != (k, l)
                    local_rdm2[k, l, i, j] += conj(value)
                end

            end


        end

    end

    # Combine the thread-local arrays into the final result.
    for t in 1:nthreads_total
        rdm2_final .+= rdm2_thread[t]
    end

    rdm2_final -= permutedims(rdm2_final, (2, 1, 3, 4))
    rdm2_final -= permutedims(rdm2_final, (1, 2, 4, 3))

    return rdm2_final

end




#rmd3

function permutation_sign(perm)
    sign = 1
    n = length(perm)
    for i in 1:n
        for j in i+1:n
            if perm[i] > perm[j]
                sign *= -1
            end
        end
    end
    return sign
end


function RDM3(ModelParams2DSpinless::ModelParams2DSpinlessList, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky
    if Nkx * Nky > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams2DSpinless.nparticle

    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)
    ind_dict = MomentumHilbertSpace2DMod.ToDict(hilbertspace)


    # Generate all unique triples (i, j, k) with i < j < k
    triples = []
    for i in 1:norb
        for j in i+1:norb
            for k in j+1:norb
                push!(triples, (i, j, k))
            end
        end
    end

    independent_elements = []
    number_of_triples = length(triples)
    for ind1 in 1:number_of_triples
        i, j, k = triples[ind1]
        for ind2 in ind1:number_of_triples
            l, m, n = triples[ind2]
            k1 = momentum_add_2d(i - 1, j - 1, Nkx, Nky)
            k1 = momentum_add_2d(k1, k - 1, Nkx, Nky)
            k2 = momentum_add_2d(l - 1, m - 1, Nkx, Nky)
            k2 = momentum_add_2d(n - 1, k2, Nkx, Nky)
            if k1 == k2
                push!(independent_elements, (i, j, k, l, m, n))
            end
        end
    end


    @show length(independent_elements)
    rdm3_final = zeros(T, norb, norb, norb, norb, norb, norb)

    DimHilbert = length(hilbert)
    nthreads_total = Threads.nthreads()
    chunk_size = ceil(Int, DimHilbert / nthreads_total)

    rdm3_thread = [zeros(T, norb, norb, norb, norb, norb, norb) for _ in 1:nthreads_total]

    Threads.@threads for t in 1:nthreads_total
        tid = Threads.threadid()
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)
        local_rdm3 = rdm3_thread[tid]
        @inbounds for stateid in start_idx:end_idx
            state = hilbert[stateid]
            for ind in eachindex(independent_elements)
                i, j, k, l, m, n = independent_elements[ind]
                ireverse = norb - i + 1
                jreverse = norb - j + 1
                kreverse = norb - k + 1
                lreverse = norb - l + 1
                mreverse = norb - m + 1
                nreverse = norb - n + 1

                fermion = FermionOperator(state, 1)

                AnnihilationOperator!(fermion, lreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                AnnihilationOperator!(fermion, mreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                AnnihilationOperator!(fermion, nreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, kreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, jreverse)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, ireverse)
                if fermion.fermion_sign == 0
                    continue
                end


                left = ind_dict[fermion.state]



                value = fermion.fermion_sign * coeffs[stateid] * conj(coeffs[left])
                local_rdm3[i, j, k, l, m, n] += value

                if (i, j, k) != (l, m, n)
                    local_rdm3[l, m, n, i, j, k] += conj(value)
                end

            end
        end
    end

    for t in 1:nthreads_total
        rdm3_final .+= rdm3_thread[t]
    end




    perms1 = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
    perms2 = [[4, 5, 6], [4, 6, 5], [5, 4, 6], [5, 6, 4], [6, 4, 5], [6, 5, 4]]

    rdm3_result = zeros(T, norb, norb, norb, norb, norb, norb)
    for eta in perms1
        for sigma in perms2
            sign = permutation_sign(eta) * permutation_sign(sigma)
            #concatenate eta and sigma
            neworder = Tuple([eta..., sigma...])
            rdm3_result += sign * permutedims(rdm3_final, neworder)
        end
    end

    return rdm3_result
end


function RDM3_single_thread(ModelParams2DSpinless::ModelParams2DSpinlessList, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky
    if Nkx * Nky > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams2DSpinless.nparticle

    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)
    ind_dict = MomentumHilbertSpace2DMod.ToDict(hilbertspace)

    # Generate all unique triples (i, j, k) with i < j < k
    triples = []
    for i in 1:norb
        for j in i+1:norb
            for k in j+1:norb
                push!(triples, (i, j, k))
            end
        end
    end

    independent_elements = []
    number_of_triples = length(triples)
    for ind1 in 1:number_of_triples
        i, j, k = triples[ind1]
        for ind2 in ind1:number_of_triples
            l, m, n = triples[ind2]
            push!(independent_elements, (i, j, k, l, m, n))
        end
    end

    rdm3 = zeros(T, norb, norb, norb, norb, norb, norb)

    @inbounds for stateid in eachindex(hilbert)
        state = hilbert[stateid]
        for ind in eachindex(independent_elements)
            i, j, k, l, m, n = independent_elements[ind]
            ireverse = norb - i + 1
            jreverse = norb - j + 1
            kreverse = norb - k + 1
            lreverse = norb - l + 1
            mreverse = norb - m + 1
            nreverse = norb - n + 1

            fermion = FermionOperator(state, 1)

            AnnihilationOperator!(fermion, lreverse)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, mreverse)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, nreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, kreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, jreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, ireverse)
            if fermion.fermion_sign == 0
                continue
            end

            if haskey(ind_dict, fermion.state)
                left = ind_dict[fermion.state]
            else
                continue
            end


            value = fermion.fermion_sign * coeffs[stateid] * conj(coeffs[left])
            rdm3[i, j, k, l, m, n] += value

            if (i, j, k) != (l, m, n)
                rdm3[l, m, n, i, j, k] += conj(value)
            end
        end

    end


    perms1 = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
    perms2 = [[4, 5, 6], [4, 6, 5], [5, 4, 6], [5, 6, 4], [6, 4, 5], [6, 5, 4]]

    rdm3_result = zeros(T, norb, norb, norb, norb, norb, norb)
    for eta in perms1
        for sigma in perms2
            sign = permutation_sign(eta) * permutation_sign(sigma)
            #concatenate eta and sigma
            neworder = Tuple([eta..., sigma...])
            rdm3_result += sign * permutedims(rdm3, neworder)
        end
    end

    return rdm3_result


end


function RDM3_naive(ModelParams2DSpinless::ModelParams2DSpinlessList, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky
    if Nkx * Nky > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams2DSpinless.nparticle

    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)
    ind_dict = MomentumHilbertSpace2DMod.ToDict(hilbertspace)



    rdm3 = zeros(T, norb, norb, norb, norb, norb, norb)

    for stateid in eachindex(hilbert)
        state = hilbert[stateid]

        for i in 1:norb, j in 1:norb, k in 1:norb, l in 1:norb, m in 1:norb, n in 1:norb


            ireverse = norb - i + 1
            jreverse = norb - j + 1
            kreverse = norb - k + 1
            lreverse = norb - l + 1
            mreverse = norb - m + 1
            nreverse = norb - n + 1

            fermion = FermionOperator(state, 1)

            AnnihilationOperator!(fermion, lreverse)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, mreverse)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, nreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, kreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, jreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, ireverse)
            if fermion.fermion_sign == 0
                continue
            end

            if haskey(ind_dict, fermion.state)
                left = ind_dict[fermion.state]
            else
                continue
            end


            value = fermion.fermion_sign * coeffs[stateid] * conj(coeffs[left])
            rdm3[i, j, k, l, m, n] += value


        end

    end

    return rdm3
end


function RDM2_naive(ModelParams2DSpinless::ModelParams2DSpinlessList, coeffs::Array{T,1}, momentum::Int) where {T}
    Nkx = ModelParams2DSpinless.Nkx
    Nky = ModelParams2DSpinless.Nky
    norb = Nkx * Nky
    if Nkx * Nky > 31
        indtype = Int64
    else
        indtype = Int32
    end

    nparticle = ModelParams2DSpinless.nparticle

    hilbertspace = MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{indtype}(nparticle, Nkx, Nky, momentum, [])
    hilbert = MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)
    ind_dict = MomentumHilbertSpace2DMod.ToDict(hilbertspace)

    rdm2 = zeros(T, norb, norb, norb, norb)
    for stateid in eachindex(hilbert)
        state = hilbert[stateid]

        for i in 1:norb, j in 1:norb, k in 1:norb, l in 1:norb


            ireverse = norb - i + 1
            jreverse = norb - j + 1
            kreverse = norb - k + 1
            lreverse = norb - l + 1

            fermion = FermionOperator(state, 1)

            AnnihilationOperator!(fermion, lreverse)
            if fermion.fermion_sign == 0
                continue
            end
            AnnihilationOperator!(fermion, kreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, jreverse)
            if fermion.fermion_sign == 0
                continue
            end
            CreationOperator!(fermion, ireverse)
            if fermion.fermion_sign == 0
                continue
            end

            if haskey(ind_dict, fermion.state)
                left = ind_dict[fermion.state]
            else
                continue
            end
            value = fermion.fermion_sign * coeffs[stateid] * conj(coeffs[left])
            rdm2[i, j, l, k] += value

        end
    end

    return rdm2

end



end

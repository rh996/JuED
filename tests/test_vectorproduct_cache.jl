# include("../src/FermionOperator.jl")
include("../src/EDMain.jl")

using .EDMod
using SparseArrays
using KrylovKit
using Arpack
using ProgressMeter

using .EDMod.FermionOperatorMod
momentum_add_2d = EDMod.MomentumHilbertSpace2DMod.momentum_add_2d
momentum_sub_2d = EDMod.MomentumHilbertSpace2DMod.momentum_sub_2d



function make_links(ind_dict, hilbertspace)
    norbital = hilbertspace.Nkx * hilbertspace.Nky
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    hilbert = hilbertspace.hilbert
    links1 = Vector{Any}(undef, length(hilbert))
    links2 = Vector{Any}(undef, length(hilbert))
    Threads.@threads for state_ind in eachindex(hilbert)
        innerlinks1 = []
        innerlinks2 = []



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
            push!(innerlinks1, i)

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
            push!(innerlinks1, i)
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
                    push!(innerlinks2, (ik, ikp, iq, fermion.fermion_sign * ind_dict[fermion.state]))
                end
            end
        end
        links2[state_ind] = innerlinks2
        links1[state_ind] = innerlinks1
    end
    return links1, links2
end



function vectormap(coeff, links1, links2, hilbertspace; U=4.0)
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    norbital = hilbertspace.Nkx * hilbertspace.Nky

    result = zeros(Float64, length(coeff))
    n_threads = Threads.nthreads()
    local_results = [zeros(Float64, length(coeff)) for _ in 1:n_threads]

    Threads.@threads for state0 in eachindex(links2)
        thread_id = Threads.threadid()
        innerlinks2 = links2[state0]
        for (ik, ikp, iq, state1) in innerlinks2
            local_results[thread_id][abs(state1)] += U / norbital * coeff[state0] * sign(state1)
        end

        innerlinks1 = links1[state0]
        for i in innerlinks1
            ikx = mod(i, Nkx)
            iky = fld(i, Nkx)
            eps = -2 * cos(2π * ikx / Nkx) - 2 * cos(2π * iky / Nky)
            local_results[thread_id][state0] += eps * coeff[state0]
        end


    end

    # Combine the results from all threads
    for t_res in local_results
        result .+= t_res
    end

    return result
end





let
    nalpha = 6
    nbeta = 6

    Nkx = 6
    Nky = 2
    k = 0
    hilbertspace = EDMod.SpinMomentumHilbertSpace2DMod.SpinMomentumHilbertSpace2D{Int64}(nalpha, nbeta, Nkx, Nky, k, [])
    EDMod.SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)
    ind_dict = EDMod.SpinMomentumHilbertSpace2DMod.ToDict(hilbertspace)
    @show dim = length(hilbertspace.hilbert)

    links1, links2 = make_links(ind_dict, hilbertspace)

    @show length(links2)
    @show length(links2[1])
    @show length(links2[2])
    @show length(links2[3])


    f = (v::Vector{Float64}) -> vectormap(v, links1, links2, hilbertspace; U=4.0)

    # f(rand(dim))
    @time energy, ψ = eigsolve(f, dim, 5, :SR; maxiter=1000, tol=1e-8, issymmetric=true, verbosity=3)

    #the memory usage of links2
    @show mem = sizeof(links2) / 1024^2
    energy

end

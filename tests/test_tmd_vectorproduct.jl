include("../src/EDMain.jl")
include("../src/FermionOperator.jl")





using .EDMod
using SparseArrays
using KrylovKit
using ProgressMeter
using HDF5
using .FermionOperatorMod
using PythonCall

using .EDMod.MomentumHilbertSpace2DMod

plt = pyimport("matplotlib.pyplot")

@inline function momentum_add_2d(k1::Int, k2::Int, Nkx::Int, Nky::Int)
    k1x = mod(k1, Nkx)
    k1y = fld(k1, Nkx)
    k2x = mod(k2, Nkx)
    k2y = fld(k2, Nkx)
    return mod(k1x + k2x, Nkx) + mod(k1y + k2y, Nky) * Nkx
end
@inline function momentum_sub_2d(k1::Int, k2::Int, Nkx::Int, Nky::Int)
    k1x = mod(k1, Nkx)
    k1y = fld(k1, Nkx)
    k2x = mod(k2, Nkx)
    k2y = fld(k2, Nkx)
    return mod(k1x - k2x, Nkx) + mod(k1y - k2y, Nky) * Nkx
end


function make_links(hilbertspace)
    norbital = hilbertspace.Nkx * hilbertspace.Nky
    nparticle = hilbertspace.nparticle
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky
    hilbert = hilbertspace.hilbert
    DimHilbert = length(hilbert)

    links1 = [zeros(Int64, nparticle) for _ in 1:DimHilbert]



    nthreads_total = Threads.nthreads()
    # @show nthreads_total
    chunk_size = ceil(Int, DimHilbert / nthreads_total)
    p = Progress(DimHilbert, desc="Counting links2")
    links2Count = zeros(UInt32, DimHilbert)
    Threads.@threads for t in 1:nthreads_total
        start_idx = (t - 1) * chunk_size + 1
        end_idx = min(t * chunk_size, DimHilbert)

        @inbounds for state_ind in start_idx:end_idx
            next!(p)

            count1 = 0
            for i in 0:norbital-1

                fermion = FermionOperator(hilbert[state_ind], 1)
                position = norbital - i

                AnnihilationOperator!(fermion, position)
                if fermion.fermion_sign == 0
                    continue
                end
                CreationOperator!(fermion, position)
                if fermion.fermion_sign == 0
                    continue
                end
                count1 += 1
                links1[state_ind][count1] = i

            end







            local_count = 0
            for ik in 0:norbital-1
                for ikp in 0:norbital-1
                    for iq in 0:norbital-1
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
            links2Count[state_ind] = UInt32(local_count)
        end

    end
    @show links2Count[1]
    # links2 = [zeros(Int32,4,links2Count[i]) for i in 1:DimHilbert]

    # p = Progress(DimHilbert, desc = "Constructing links2")
    # @inbounds Threads.@threads for state_ind in eachindex(hilbert)
    #     next!(p)
    #     localCount = 0
    #     for ik in 0:norbital-1
    #         for ikp in 0:norbital-1
    #             for iq in 0:norbital-1
    #                 ind1 = momentum_add_2d(ik,iq,Nkx,Nky)
    #                 ind1 = norbital - ind1
    #                 ind2 = ik
    #                 ind2 = norbital - ind2
    #                 ind3 = momentum_sub_2d(ikp,iq,Nkx,Nky)
    #                 ind3 = norbital - ind3
    #                 ind4 = ikp
    #                 ind4 = norbital - ind4
    #                 fermion = FermionOperator(hilbert[state_ind], 1)
    #                 AnnihilationOperator!(fermion, ind2)
    #                 if fermion.fermion_sign ==0
    #                     continue
    #                 end
    #                 AnnihilationOperator!(fermion, ind4)
    #                 if fermion.fermion_sign ==0
    #                     continue
    #                 end
    #                 CreationOperator!(fermion, ind3)
    #                 if fermion.fermion_sign ==0
    #                     continue
    #                 end
    #                 CreationOperator!(fermion, ind1)
    #                 if fermion.fermion_sign ==0
    #                     continue
    #                 end
    #                 localCount += 1
    #                 links2[state_ind][:,localCount] = [ik,ikp,iq,fermion.fermion_sign*ind_dict[fermion.state]]

    #             end
    #         end
    #     end

    # end



    # return links1,links2
end



function vectormap(OneBody, TwoBody, coeff, links1, links2)


    result = zeros(ComplexF64, length(coeff))
    n_threads = Threads.nthreads()
    local_results = [zeros(ComplexF64, length(coeff)) for _ in 1:n_threads]

    Threads.@threads for state0 in eachindex(links2)
        thread_id = Threads.threadid()
        innerlinks2 = links2[state0]
        # @show innerlinks2
        for excitations in eachcol(innerlinks2)
            ik = excitations[1]
            ikp = excitations[2]
            iq = excitations[3]
            state1 = excitations[4]

            local_results[thread_id][abs(state1)] += TwoBody[ik+1, ikp+1, iq+1] * coeff[state0] * sign(state1)
        end

        innerlinks1 = links1[state0]
        for i in innerlinks1
            eps = OneBody[i+1]
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
    nparticle = 8

    Nkx = 4
    Nky = 8

    filename = "./data/twoband_TMD_degree_3.5_gate_10.0_nkx_$(Nkx)_nky_$(Nky).h5"
    interactions = h5read(filename, "interactions")
    band = h5read(filename, "band")

    OneBody = zeros(ComplexF64, Nkx * Nky)
    TwoBody = zeros(ComplexF64, Nkx * Nky, Nkx * Nky, Nkx * Nky)
    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            ik = ikx + iky * Nkx
            OneBody[ik+1] = band[ikx+1, iky+1, 1]
        end
    end

    for ikx in 0:Nkx-1
        for iky in 0:Nky-1
            for ikpx in 0:Nkx-1
                for ikpy in 0:Nky-1
                    for iqx in 0:Nkx-1
                        for iqy in 0:Nky-1
                            ik = ikx + iky * Nkx
                            ikp = ikpx + ikpy * Nkx
                            iq = iqx + iqy * Nkx
                            TwoBody[ik+1, ikp+1, iq+1] = interactions[1, 1, 1, 1, ikx+1, iky+1, ikpx+1, ikpy+1, iqx+1, iqy+1]
                        end
                    end
                end
            end
        end
    end

    elist = []
    for k in 0:Nkx*Nky-1
        # k = 0
        hilbertspace = EDMod.MomentumHilbertSpace2DMod.MomentumHilbertSpace2D{Int64}(nparticle, Nkx, Nky, k, [])
        EDMod.MomentumHilbertSpace2DMod.BuildHilbert(hilbertspace)

        ind_dict = EDMod.MomentumHilbertSpace2DMod.ToDict(hilbertspace)

        # @show sizeof(hilbertspace.hilbert)/1024^2

        @show dim = length(hilbertspace.hilbert)

        make_links(hilbertspace)
        links1, links2 = make_links(hilbertspace)
        # println("Total memory consumption: ", Base.summarysize(links2)/1024^2, " MB")
        # @show size(links)

        # @show links2[1][:,1]

        # @show length(links2)
        # @show length(links2[1])
        # @show length(links2[2])
        # @show length(links2[3])

        # @show mem = sizeof(links2)/1024^2


        f = (v::Vector{ComplexF64}) -> vectormap(OneBody, TwoBody, v, links1, links2)



        @time energy, ψ = eigsolve(f, dim, 5, :SR, ComplexF64; maxiter=1000, tol=1e-6, ishermitian=true, verbosity=3)
        push!(elist, energy[1:5])

    end
    elist = hcat(elist...)
    X = [[i for i in 1:Nkx*Nky] for j in 1:5]
    X = hcat(X...)

    plt.scatter(X, elist')
    plt.show()
    plt.close()


end

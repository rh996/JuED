include("../src/FermionOperator.jl")
include("../src/SpinMomentumHilbertSpace1D.jl")

using .SpinMomentumHilbertSpace1DMod
using SparseArrays
using KrylovKit
using Arpack
using ProgressMeter
using .FermionOperatorMod


nalpha = 6
nbeta = 6
orbital = 12

AbstractHilbertSpace = SpinMomentumHilbertSpace1DMod.AbstractHilbertSpace


function count_data_size(hilbertspace::AbstractHilbertSpace)
    hilbert = SpinMomentumHilbertSpace1DMod.BuildSpinHilbert(hilbertspace)
    m = length(hilbert)
    ptrlenth = m + 1
    datalenth = 0

    p = Progress(m)
    @inbounds for state_ind in eachindex(hilbert)
        next!(p)
        for i in 0:hilbertspace.norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = hilbertspace.norbital - i
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

        for i in 0:hilbertspace.norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = hilbertspace.norbital - i
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

        for ik in 0:hilbertspace.norbital-1
            for ikp in 0:hilbertspace.norbital-1
                for iq in 0:hilbertspace.norbital-1
                    ind1 = 2 * mod((ik + iq), hilbertspace.norbital)
                    ind1 = 2 * hilbertspace.norbital - ind1
                    ind2 = 2 * ik
                    ind2 = 2 * hilbertspace.norbital - ind2
                    ind3 = 2 * mod((ikp - iq), hilbertspace.norbital) + 1
                    ind3 = 2 * hilbertspace.norbital - ind3
                    ind4 = 2 * ikp + 1
                    ind4 = 2 * hilbertspace.norbital - ind4
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




function Hubbard1D(U::Real, datasize::Int64, hilbertspace::AbstractHilbertSpace)
    hilbert = hilbertspace.hilbert
    ind_dict = SpinMomentumHilbertSpace1DMod.ToDict(hilbertspace)
    m = length(hilbert)
    # println("Sector Dim: ", m)

    datatp = typeof(U)
    data::Array{datatp,1} = zeros(datatp, datasize)
    row::Array{Int32,1} = zeros(Int32, datasize)
    indptr::Array{Int32,1} = zeros(Int32, m + 1)
    indptr[1] = 1
    COUNT = 0
    @inbounds for state_ind in eachindex(hilbert)
        ncolcount = 0
        for i in 0:hilbertspace.norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = hilbertspace.norbital - i
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
            # push!(data, -2*cos(2*pi*i/hilbertspace.norbital)*fermion.fermion_sign)
            # push!(row, ind_dict[fermion.state])
            data[COUNT] = -2 * cos(2 * pi * i / hilbertspace.norbital) * fermion.fermion_sign
            row[COUNT] = ind_dict[fermion.state]
            ncolcount += 1

        end

        for i in 0:hilbertspace.norbital-1

            fermion = FermionOperator(hilbert[state_ind], 1)
            position = hilbertspace.norbital - i
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
            data[COUNT] = -2 * cos(2 * pi * i / hilbertspace.norbital) * fermion.fermion_sign
            row[COUNT] = ind_dict[fermion.state]
            # push!(data, -2*cos(2*pi*i/hilbertspace.norbital)*fermion.fermion_sign)
            # push!(row, ind_dict[fermion.state])
            ncolcount += 1
        end

        for ik in 0:hilbertspace.norbital-1
            for ikp in 0:hilbertspace.norbital-1
                for iq in 0:hilbertspace.norbital-1
                    ind1 = 2 * mod((ik + iq), hilbertspace.norbital)
                    ind1 = 2 * hilbertspace.norbital - ind1
                    ind2 = 2 * ik
                    ind2 = 2 * hilbertspace.norbital - ind2
                    ind3 = 2 * mod((ikp - iq), hilbertspace.norbital) + 1
                    ind3 = 2 * hilbertspace.norbital - ind3
                    ind4 = 2 * ikp + 1
                    ind4 = 2 * hilbertspace.norbital - ind4
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
                    # push!(data, (U/hilbertspace.norbital)*fermion.fermion_sign)
                    # push!(row, ind_dict[fermion.state])
                    COUNT += 1
                    data[COUNT] = (U / hilbertspace.norbital) * fermion.fermion_sign
                    row[COUNT] = ind_dict[fermion.state]
                    ncolcount += 1
                end
            end
        end
        # push!(indptr, ncolcount+indptr[end])
        indptr[state_ind+1] = ncolcount + indptr[state_ind]
    end

    return data, row, indptr, m

end


function test(nalpha, nbeta)
    elist = []
    for k in 0:orbital-1
        hilbertspace = SpinMomentumHilbertSpace1DMod.SpinMomentumHilbertSpace1D{Int64}(nalpha, nbeta, orbital, k, [])
        datasize, _ = count_data_size(hilbertspace)
        data, row, indptr, dim = Hubbard1D(4.0, datasize, hilbertspace)
        # vals,vecs = eigs(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), nev=5, which=:SR)
        vals, _, _ = eigsolve(SparseMatrixCSC{Float64,Int32}(dim, dim, indptr, row, data), ones(Float64, dim), 5, :SR, Lanczos())
        push!(elist, vals[1:5])

    end
    println(elist)
end



@time test(nalpha, nbeta)

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


function vectormap(coeff, ind_dict, hilbert, hilbertspace; U=4.0)
    result = zeros(Float64, length(coeff))


    norbital = hilbertspace.Nkx * hilbertspace.Nky
    Nkx = hilbertspace.Nkx
    Nky = hilbertspace.Nky

    @inbounds for state_ind in eachindex(hilbert)

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

            ikx = mod(i, Nkx)
            iky = fld(i, Nkx)
            eps = -2 * cos(2π * ikx / Nkx) - 2 * cos(2π * iky / Nky)
            result[ind_dict[fermion.state]] += eps * fermion.fermion_sign * coeff[state_ind]
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
            ikx = mod(i, Nkx)
            iky = fld(i, Nkx)
            eps = -2 * cos(2π * ikx / Nkx) - 2 * cos(2π * iky / Nky)

            result[ind_dict[fermion.state]] += eps * fermion.fermion_sign * coeff[state_ind]
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
                    result[ind_dict[fermion.state]] += (U / norbital) * fermion.fermion_sign * coeff[state_ind]
                end
            end
        end

    end

    return result
end


let
    nalpha = 4
    nbeta = 4

    Nkx = 4
    Nky = 2
    k = 0
    hilbertspace = EDMod.SpinMomentumHilbertSpace2DMod.SpinMomentumHilbertSpace2D{Int64}(nalpha, nbeta, Nkx, Nky, k, [])
    hilbert = EDMod.SpinMomentumHilbertSpace2DMod.BuildSpinHilbert(hilbertspace)
    ind_dict = EDMod.SpinMomentumHilbertSpace2DMod.ToDict(hilbertspace)
    @show dim = length(hilbert)

    f = (v::Vector{Float64}) -> vectormap(v, ind_dict, hilbert, hilbertspace; U=4.0)

    # f(rand(dim))
    @time energy, ψ = eigsolve(f, dim, 5, :SR; maxiter=1000, tol=1e-8, issymmetric=true, verbosity=3)

    energy
end

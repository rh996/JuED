using LinearAlgebra
using Random
using Test

if isdefined(Main, :JuED)
    const EDMod = Main.JuED.EDMod
else
    include("../src/EDMain.jl")
    const EDMod = Main.EDMod
end

function find_nonempty_spinless_workspace(model)
    for momentum in 0:(model.Nkx * model.Nky - 1)
        workspace = EDMod.RDMWorkspace(model, momentum)
        if !isempty(workspace.hilbert)
            return workspace
        end
    end
    error("No nonempty momentum sector found.")
end

function random_state_density_invariants(dim::Int)
    coeffs = randn(ComplexF64, dim)
    coeffs ./= norm(coeffs)
    return coeffs
end

function contract_rdm2(rdm2, nparticle::Int)
    norbital = size(rdm2, 1)
    contracted = zeros(eltype(rdm2), norbital, norbital)
    for i in 1:norbital
        for j in 1:norbital
            for k in 1:norbital
                contracted[i, j] += rdm2[i, k, j, k] / nparticle
            end
        end
    end
    return contracted
end

@testset "Density-matrix traces and contractions preserve particle number" begin
    Random.seed!(20260316)
    nparticle = 2
    model = EDMod.SpinlessListModel(nparticle, 2, 2, zeros(ComplexF64, 4, 4), zeros(ComplexF64, 4, 4, 4, 4))
    workspace = find_nonempty_spinless_workspace(model)
    coeffs = random_state_density_invariants(length(workspace.hilbert))

    rdm1 = EDMod.RDM1(workspace, coeffs)
    rdm2 = EDMod.RDM2(workspace, coeffs)

    @test isapprox(tr(rdm1), nparticle; atol=1e-10)
    @test isapprox(contract_rdm2(rdm2, nparticle), rdm1; atol=1e-10)
    @test isapprox(rdm1, rdm1'; atol=1e-10)
    @test isapprox(rdm2, permutedims(conj(rdm2), (3, 4, 1, 2)); atol=1e-10)
end

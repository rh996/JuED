@testset "Density-matrix traces and contractions preserve particle number" begin
    Random.seed!(20260316)
    nparticle = 2
    model = SpinlessListModel(nparticle, 2, 2, zeros(ComplexF64, 4, 4), zeros(ComplexF64, 4, 4, 4, 4))
    workspace = find_nonempty_workspace(model)
    coeffs = random_state(length(workspace.hilbert))

    rdm1 = RDM1(workspace, coeffs)
    rdm2 = RDM2(workspace, coeffs)

    @test isapprox(tr(rdm1), nparticle; atol=1e-10)
    @test isapprox(contract_rdm2(rdm2, nparticle), rdm1; atol=1e-10)
    @test isapprox(rdm1, rdm1'; atol=1e-10)
    @test isapprox(rdm2, permutedims(conj(rdm2), (3, 4, 1, 2)); atol=1e-10)
end
